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
 * Implements the Reichardt (1951) universal law of the wall:
 *   u⁺ = y⁺/(1 + y⁺/11) + (1/κ)*ln(1 + 0.41*y⁺) + 7.8*[1 - exp(-y⁺/11) - (y⁺/11)*exp(-y⁺/3)]
 * 
 * Reference:
 *   Reichardt, H. (1951). "Vollständige Darstellung der turbulenten 
 *   Geschwindigkeitsverteilung in Rohren." Zeitschrift für angewandte 
 *   Mathematik und Mechanik, 31(7-8), 208-219.
 */
#include <cmath>
#include <algorithm>

// Reichardt constants
#define C_REICH 0.41
#define A_DAMP 11.0
#define B_EXP 3.0
#define EXP_COEFF 7.8

/**
 * Dimensionless velocity u⁺ from dimensionless distance y⁺ 
 * using the original Reichardt (1951) formula.
 */
template<typename T>
__device__ T reichardt_u_plus(const T y_plus, const T kappa) {
    if (y_plus < 1.0e-12) {
        return y_plus;
    }

    // Linear term with damping: y⁺/(1 + y⁺/11)
    T linear_term = y_plus / (1.0 + y_plus / A_DAMP);

    // Logarithmic term: (1/κ)*ln(1 + 0.41*y⁺)
    T log_arg = 1.0 + C_REICH * y_plus;
    if (log_arg <= 0.0) {
        return linear_term;
    }
    T log_term = (1.0 / kappa) * log(log_arg);

    // Exponential correction term: 7.8*[1 - exp(-y⁺/11) - (y⁺/11)*exp(-y⁺/3)]
    T exp_term1 = exp(-y_plus / A_DAMP);
    T exp_term2 = (y_plus / A_DAMP) * exp(-y_plus / B_EXP);

    // Total dimensionless velocity
    T u_plus = linear_term + log_term + EXP_COEFF * (1.0 - exp_term1 - exp_term2);

    return u_plus;
}

/**
 * Derivative du⁺/dy⁺ of the original Reichardt (1951) formula.
 * Used for Newton-Raphson iteration.
 */
template<typename T>
__device__ T reichardt_du_plus_dy(const T y_plus, const T kappa) {
    if (y_plus < 1.0e-12) {
        return 1.0;
    }

    // Derivative of linear term: d/dy[y/(1 + y/11)] = 1/(1 + y/11)²
    T d_linear = 1.0 / ((1.0 + y_plus / A_DAMP) * (1.0 + y_plus / A_DAMP));

    // Derivative of log term: d/dy[(1/κ)*ln(1 + 0.41*y)]
    T log_arg = 1.0 + C_REICH * y_plus;
    if (log_arg <= 0.0) {
        return d_linear;
    }
    T d_log = (1.0 / kappa) * C_REICH / log_arg;

    // Derivative of exponential term:
    // d/dy[7.8*(1 - exp(-y/11) - (y/11)*exp(-y/3))]
    T exp_term1 = exp(-y_plus / A_DAMP) / A_DAMP;
    T exp_term2 = (1.0 / A_DAMP) * exp(-y_plus / B_EXP) - 
                  (y_plus / A_DAMP) * (1.0 / B_EXP) * exp(-y_plus / B_EXP);
    T d_exp = EXP_COEFF * (exp_term1 + exp_term2);

    T du_dy = d_linear + d_log + d_exp;

    return du_dy;
}

/**
 * Newton-Raphson solver for friction velocity using original Reichardt law.
 * Solves the implicit equation: F(u_τ) = u_τ * u⁺(y⁺) - U = 0
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

        // Current dimensionless coordinates
        y_plus = y * utau / nu;
        u_plus = reichardt_u_plus(y_plus, kappa);
        du_dy = reichardt_du_plus_dy(y_plus, kappa);

        // Residual: F(u_τ) = u_τ * u⁺ - U
        f = utau * u_plus - u;

        // Sensitivity: dF/du_τ = u⁺ + u_τ * du⁺/dy⁺ * (y/ν)
        df = u_plus + utau * du_dy * (y / nu);

        // Safeguard against singular Jacobian
        if (fabs(df) < 1.0e-14) {
            utau = utau * 0.99;
            continue;
        }

        // Newton step
        utau = utau - f / df;

        // Safeguard: prevent u_τ from becoming negative or zero
        if (utau <= 0.0) {
            utau = utau_old * 0.5;
        }

        // Relative error check for convergence
        error = fabs((utau - utau_old) / (utau + 1.0e-16));

        // Convergence tolerance
        if (error < 1.0e-8) {
            break;
        }
    }

    return utau;
}

/**
 * HIP/CUDA kernel for Reichardt's wall model.
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
        // Sample the velocity at the off-wall point
        const int index = (ind_e_d[i] - 1) * lx * lx * lx +
                          (ind_t_d[i] - 1) * lx * lx +
                          (ind_s_d[i] - 1) * lx +
                          (ind_r_d[i] - 1);

        T ui = u_d[index];
        T vi = v_d[index];
        T wi = w_d[index];

        // Load normal vectors and wall distance
        T nx = n_x_d[i];
        T ny = n_y_d[i];
        T nz = n_z_d[i];
        T h = h_d[i];
        T nu = nu_d[i];

        // Project velocity onto tangential plane (remove normal component)
        T normu = ui * nx + vi * ny + wi * nz;

        ui -= normu * nx;
        vi -= normu * ny;
        wi -= normu * nz;

        // Magnitude of tangential velocity
        T magu = sqrt(ui * ui + vi * vi + wi * wi);

        // Get initial guess for Newton solver
        T guess;
        if (tstep == 1) {
            // First time-step: use simple estimate
            guess = sqrt(magu * nu / h);
        } else {
            // Use previous solution as starting point
            T tau_mag_sq = tau_x_d[i] * tau_x_d[i] +
                           tau_y_d[i] * tau_y_d[i] +
                           tau_z_d[i] * tau_z_d[i];
            guess = sqrt(sqrt(tau_mag_sq));
        }

        // Solve for friction velocity u_tau
        T utau = solve_reichardt(magu, h, guess, nu, kappa);

        // Distribute shear stress according to velocity direction
        // tau_wall = -ρ * u_tau² * (u_tangential / |u_tangential|)
        if (magu > 1.0e-14) {
            tau_x_d[i] = -utau * utau * ui / magu;
            tau_y_d[i] = -utau * utau * vi / magu;
            tau_z_d[i] = -utau * utau * wi / magu;
        } else {
            // Avoid division by zero when velocity is zero
            tau_x_d[i] = 0.0;
            tau_y_d[i] = 0.0;
            tau_z_d[i] = 0.0;
        }
    }
}

#endif // __COMMON_REICHARDT_KERNEL_H__