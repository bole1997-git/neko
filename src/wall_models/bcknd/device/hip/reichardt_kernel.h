#ifndef __COMMON_REICHARDT_KERNEL_H__
#define __COMMON_REICHARDT_KERNEL_H__
/*
 Copyright (c) 2025, The Neko Authors
 All rights reserved.

 Redistribution and use in source and binary forms, with or without
 modification, are permitted provided that the following conditions
 are met:

   * Redistributions of source code must retain the above copyright
     notice, this list of conditions and the following disclaimer.

   * Redistributions in binary form must reproduce the above
     copyright notice, this list of conditions and the following
     disclaimer in the documentation and/or other materials provided
     with the distribution.

   * Neither the name of the authors nor the names of its
     contributors may be used to endorse or promote products derived
     from this software without specific prior written permission.

 THIS SOFTWARE IS PROVIDED BY THE COPYRIGHT HOLDERS AND CONTRIBUTORS
 "AS IS" AND ANY EXPRESS OR IMPLIED WARRANTIES, INCLUDING, BUT NOT
 LIMITED TO, THE IMPLIED WARRANTIES OF MERCHANTABILITY AND FITNESS
 FOR A PARTICULAR PURPOSE ARE DISCLAIMED. IN NO EVENT SHALL THE
 COPYRIGHT OWNER OR CONTRIBUTORS BE LIABLE FOR ANY DIRECT, INDIRECT,
 INCIDENTAL, SPECIAL, EXEMPLARY, OR CONSEQUENTIAL DAMAGES (INCLUDING,
 BUT NOT LIMITED TO, PROCUREMENT OF SUBSTITUTE GOODS OR SERVICES;
 LOSS OF USE, DATA, OR PROFITS; OR BUSINESS INTERRUPTION) HOWEVER
 CAUSED AND ON ANY THEORY OF LIABILITY, WHETHER IN CONTRACT, STRICT
 LIABILITY, OR TORT (INCLUDING NEGLIGENCE OR OTHERWISE) ARISING IN
 ANY WAY OUT OF THE USE OF THIS SOFTWARE, EVEN IF ADVISED OF THE
 POSSIBILITY OF SUCH DAMAGE.
*/

/**
 * Device kernel for reichardt_compute
 *
 * Implements the original two-term Reichardt (1951) universal law of the wall:
 *
 *   u⁺ = (1/κ)*ln(1 + κ*y⁺) + 7.8*[1 - exp(-y⁺/11) - (y⁺/11)*exp(-y⁺/3)]
 *
 * Reference:
 *   Reichardt, H. (1951). "Vollständige Darstellung der turbulenten
 *   Geschwindigkeitsverteilung in Rohren." Zeitschrift für angewandte
 *   Mathematik und Mechanik, 31(7-8), 208-219.
 */
#include <cmath>
#include <algorithm>

// Reichardt (1951) constants
#define A_DAMP   11.0   // Damping length scale
#define B_EXP     3.0   // Exponential decay scale
#define EXP_COEFF 7.8   // Exponential amplitude

/**
 * Dimensionless velocity u⁺ from dimensionless distance y⁺
 * using the original two-term Reichardt (1951) formula.
 *
 *   u⁺ = (1/κ)*ln(1 + κ*y⁺) + 7.8*[1 - exp(-y⁺/11) - (y⁺/11)*exp(-y⁺/3)]
 *
 * Limiting behaviour:
 *   - y⁺ → 0 : u⁺ → y⁺  (viscous sublayer recovered analytically)
 *   - y⁺ → ∞ : u⁺ → (1/κ)*ln(y⁺) + B  (log law)
 */
template<typename T>
__device__ T reichardt_u_plus(const T y_plus, const T kappa) {
    if (y_plus < 1.0e-12) {
        return y_plus;
    }

    // Logarithmic term: (1/κ)*ln(1 + κ*y⁺)
    T log_arg = 1.0 + kappa * y_plus;
    T log_term = (1.0 / kappa) * log(log_arg);

    // Exponential correction: 7.8*[1 - exp(-y⁺/11) - (y⁺/11)*exp(-y⁺/3)]
    T exp_term1 = exp(-y_plus / A_DAMP);
    T exp_term2 = (y_plus / A_DAMP) * exp(-y_plus / B_EXP);

    return log_term + EXP_COEFF * (1.0 - exp_term1 - exp_term2);
}

/**
 * Analytical derivative du⁺/dy⁺ of the original Reichardt (1951) formula.
 *
 * Derived from:
 *   u⁺ = (1/κ)*ln(1 + κ*y⁺) + 7.8*[1 - exp(-y⁺/11) - (y⁺/11)*exp(-y⁺/3)]
 *
 * Giving:
 *   du⁺/dy⁺ = 1/(1 + κ*y⁺)
 *            + 7.8*[ (1/11)*exp(-y⁺/11)  -  (1/11)*exp(-y⁺/3)  +  (y⁺/33)*exp(-y⁺/3) ]
 *
 * Note: the sign of the last two exponential sub-terms is MINUS (they are
 * grouped as exp_term2 = (1/11)*exp(-y/3) - (y/33)*exp(-y/3), so
 * d_exp = 7.8*(exp_term1 - exp_term2)).
 *
 * Used in Newton-Raphson via:
 *   dF/du_τ = u⁺ + u_τ * (du⁺/dy⁺) * (y/ν)
 */
template<typename T>
__device__ T reichardt_du_plus_dy(const T y_plus, const T kappa) {
    if (y_plus < 1.0e-12) {
        return 1.0;  // viscous sublayer: du⁺/dy⁺ = 1
    }

    // Derivative of log term: d/dy[(1/κ)*ln(1 + κ*y)] = 1/(1 + κ*y)
    T log_arg = 1.0 + kappa * y_plus;
    T d_log = 1.0 / log_arg;

    // Derivative of exponential correction:
    //   d/dy[7.8*(1 - exp(-y/11) - (y/11)*exp(-y/3))]
    // = 7.8*[ (1/11)*exp(-y/11)  -  (1/11)*exp(-y/3)  +  (y/33)*exp(-y/3) ]
    //
    // Grouping the last two terms into exp_term2:
    //   exp_term2 = (1/11)*exp(-y/3) - (y/33)*exp(-y/3)
    // so  d_exp = 7.8*(exp_term1 - exp_term2)   <-- MINUS sign
    T exp_term1 = exp(-y_plus / A_DAMP) / A_DAMP;
    T exp_term2 = (1.0 / A_DAMP) * exp(-y_plus / B_EXP)
                - (y_plus / A_DAMP) * (1.0 / B_EXP) * exp(-y_plus / B_EXP);
    T d_exp = EXP_COEFF * (exp_term1 - exp_term2);

    return d_log + d_exp;
}

/**
 * Newton-Raphson solver for friction velocity u_τ using the Reichardt law.
 *
 * Solves the implicit equation:
 *   F(u_τ) = u_τ * u⁺(y⁺) - U = 0
 * where y⁺ = y * u_τ / ν and U is the tangential velocity magnitude.
 *
 * Jacobian (chain rule through y⁺):
 *   dF/du_τ = u⁺ + u_τ * (du⁺/dy⁺) * (y/ν)
 */
template<typename T>
__device__ T solve_reichardt(const T u, const T y, const T guess, const T nu,
                             const T kappa) {
    T utau = guess;
    T y_plus, u_plus, du_dy;
    T error, f, df, utau_old;
    const int maxiter = 100;

    for (int k = 0; k < maxiter; ++k) {
        utau_old = utau;

        y_plus = y * utau / nu;
        u_plus = reichardt_u_plus(y_plus, kappa);
        du_dy  = reichardt_du_plus_dy(y_plus, kappa);

        // Residual and Jacobian
        f  = utau * u_plus - u;
        df = u_plus + utau * du_dy * (y / nu);

        // Guard against near-zero Jacobian
        if (fabs(df) < 1.0e-14) {
            utau = utau * 0.99;
            continue;
        }

        // Newton step
        utau = utau - f / df;

        // Keep u_τ strictly positive
        if (utau <= 0.0) {
            utau = utau_old * 0.5;
        }

        // Convergence check on relative change
        error = fabs((utau - utau_old) / (utau + 1.0e-16));
        if (error < 1.0e-8) {
            break;
        }
    }

    return utau;
}

/**
 * CUDA/HIP kernel for Reichardt's wall model.
 */
template<typename T>
__global__ void reichardt_compute(const T * __restrict__ u_d,
                                  const T * __restrict__ v_d,
                                  const T * __restrict__ w_d,
                                  const int * __restrict__ ind_r_d,
                                  const int * __restrict__ ind_s_d,
                                  const int * __restrict__ ind_t_d,
                                  const int * __restrict__ ind_e_d,
                                  const T * __restrict__ n_x_d,
                                  const T * __restrict__ n_y_d,
                                  const T * __restrict__ n_z_d,
                                  const T * __restrict__ nu_d,
                                  const T * __restrict__ h_d,
                                  T * __restrict__ tau_x_d,
                                  T * __restrict__ tau_y_d,
                                  T * __restrict__ tau_z_d,
                                  const int n_nodes,
                                  const int lx,
                                  const T kappa,
                                  const T B,
                                  const int tstep) {

    const int idx = blockIdx.x * blockDim.x + threadIdx.x;
    const int str = blockDim.x * gridDim.x;

    for (int i = idx; i < n_nodes; i += str) {
        // Sample velocity at the off-wall point (Fortran-1 index offset)
        const int index = (ind_e_d[i] - 1) * lx * lx * lx +
                          (ind_t_d[i] - 1) * lx * lx +
                          (ind_s_d[i] - 1) * lx +
                          (ind_r_d[i] - 1);

        T ui = u_d[index];
        T vi = v_d[index];
        T wi = w_d[index];

        const T nx = n_x_d[i];
        const T ny = n_y_d[i];
        const T nz = n_z_d[i];
        const T h  = h_d[i];
        const T nu = nu_d[i];

        // Remove wall-normal component to get tangential velocity
        T normu = ui * nx + vi * ny + wi * nz;
        ui -= normu * nx;
        vi -= normu * ny;
        wi -= normu * nz;

        // Magnitude of tangential velocity
        T magu = sqrt(ui * ui + vi * vi + wi * wi);

        // Initial guess for Newton solver
        T guess;
        if (tstep == 1) {
            // First time-step: laminar sublayer estimate
            guess = sqrt(magu * nu / h);
        } else {
            // Warm start from previous shear stress magnitude
            T tau_mag_sq = tau_x_d[i] * tau_x_d[i] +
                           tau_y_d[i] * tau_y_d[i] +
                           tau_z_d[i] * tau_z_d[i];
            guess = sqrt(sqrt(tau_mag_sq));
        }

        // Solve for friction velocity u_tau
        T utau = solve_reichardt(magu, h, guess, nu, kappa);

        // Distribute shear stress in the tangential velocity direction
        // tau_wall = -u_tau² * (u_tangential / |u_tangential|)
        if (magu > 1.0e-14) {
            tau_x_d[i] = -utau * utau * ui / magu;
            tau_y_d[i] = -utau * utau * vi / magu;
            tau_z_d[i] = -utau * utau * wi / magu;
        } else {
            tau_x_d[i] = 0.0;
            tau_y_d[i] = 0.0;
            tau_z_d[i] = 0.0;
        }
    }
}

#endif // __COMMON_REICHARDT_KERNEL_H__