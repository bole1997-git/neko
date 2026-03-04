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

#include <hip/hip_runtime.h>
#include <device/device_config.h>
#include <device/hip/check.h>
#include "reichardt_kernel.h"

/**
 * HIP C wrapper for the Reichardt wall model computation.
 * 
 * This function launches the HIP kernel that computes wall shear stress
 * using the Reichardt (1951) universal law of the wall on AMD GPUs.
 * 
 * @param u_d GPU pointer to velocity x-component
 * @param v_d GPU pointer to velocity y-component
 * @param w_d GPU pointer to velocity z-component
 * @param ind_r_d GPU pointer to r-direction indices
 * @param ind_s_d GPU pointer to s-direction indices
 * @param ind_t_d GPU pointer to t-direction indices
 * @param ind_e_d GPU pointer to element indices
 * @param n_x_d GPU pointer to x-component of wall normal
 * @param n_y_d GPU pointer to y-component of wall normal
 * @param n_z_d GPU pointer to z-component of wall normal
 * @param nu_d GPU pointer to kinematic viscosity
 * @param h_d GPU pointer to wall-normal distance
 * @param tau_x_d GPU pointer to output x-component of wall shear stress
 * @param tau_y_d GPU pointer to output y-component of wall shear stress
 * @param tau_z_d GPU pointer to output z-component of wall shear stress
 * @param n_nodes Number of boundary nodes
 * @param lx Polynomial order in element
 * @param kappa Von Kármán constant
 * @param B Log-law intercept (not used in Reichardt formula)
 * @param tstep Current time-step
 */
extern "C" {
  void hip_reichardt_compute(void *u_d, void *v_d, void *w_d,
          void *ind_r_d, void *ind_s_d, void *ind_t_d, void *ind_e_d,
          void *n_x_d, void *n_y_d, void *n_z_d, void *nu_d, void *h_d,
          void *tau_x_d, void *tau_y_d, void *tau_z_d,
          int *n_nodes, int *lx, real *kappa, real *B, int *tstep) {

    // Define HIP grid and block dimensions
    // Use 1024 threads per block for good GPU occupancy
    const dim3 nthrds(1024, 1, 1);
    const dim3 nblcks(((*n_nodes) + 1024 - 1) / 1024, 1, 1);
    
    // Launch kernel only if there are nodes to process
    if (*n_nodes > 0) {
      hipLaunchKernelGGL(HIP_KERNEL_NAME(reichardt_compute<real>),
                         nblcks, nthrds, 0, (hipStream_t) glb_cmd_queue,
                         (real *) u_d, (real *) v_d, (real *) w_d,
                         (int *) ind_r_d, (int *) ind_s_d,
                         (int *) ind_t_d, (int *) ind_e_d,
                         (real *) n_x_d, (real *) n_y_d,
                         (real *) n_z_d, (real *) nu_d, (real *) h_d,
                         (real *) tau_x_d, (real *) tau_y_d,
                         (real *) tau_z_d,
                         *n_nodes, *lx, *kappa, *B, *tstep);
      HIP_CHECK(hipGetLastError());
    }
  }
}
