! Copyright (c) 2025, The Neko Authors
! All rights reserved.
!
! Redistribution and use in source and binary forms, with or without
! modification, are permitted provided that the following conditions
! are met:
!
!   * Redistributions of source code must retain the above copyright
!     notice, this list of conditions and the following disclaimer.
!
!   * Redistributions in binary form must reproduce the above
!     copyright notice, this list of conditions and the following
!     disclaimer in the documentation and/or other materials provided
!     with the distribution.
!
!   * Neither the name of the authors nor the names of its
!     contributors may be used to endorse or promote products derived
!     from this software without specific prior written permission.
!
! THIS SOFTWARE IS PROVIDED BY THE COPYRIGHT HOLDERS AND CONTRIBUTORS
! "AS IS" AND ANY EXPRESS OR IMPLIED WARRANTIES, INCLUDING, BUT NOT
! LIMITED TO, THE IMPLIED WARRANTIES OF MERCHANTABILITY AND FITNESS
! FOR A PARTICULAR PURPOSE ARE DISCLAIMED. IN NO EVENT SHALL THE
! COPYRIGHT OWNER OR CONTRIBUTORS BE LIABLE FOR ANY DIRECT, INDIRECT,
! INCIDENTAL, SPECIAL, EXEMPLARY, OR CONSEQUENTIAL DAMAGES (INCLUDING,
! BUT NOT LIMITED TO, PROCUREMENT OF SUBSTITUTE GOODS OR SERVICES;
! LOSS OF USE, DATA, OR PROFITS; OR BUSINESS INTERRUPTION) HOWEVER
! CAUSED AND ON ANY THEORY OF LIABILITY, WHETHER IN CONTRACT, STRICT
! LIABILITY, OR TORT (INCLUDING NEGLIGENCE OR OTHERWISE) ARISING IN
! ANY WAY OUT OF THE USE OF THIS SOFTWARE, EVEN IF ADVISED OF THE
! POSSIBILITY OF SUCH DAMAGE.
!
!> Fortran-C interface for the device (GPU) kernel of `reichardt_t`.
!!
!! This module provides the Fortran binding to the CUDA/HIP C kernel that
!! computes wall shear stress using the original two-term Reichardt (1951)
!! law of the wall:
!!
!!   u⁺ = (1/κ)*ln(1 + κ*y⁺) + 7.8*[1 - exp(-y⁺/11) - (y⁺/11)*exp(-y⁺/3)]
!!
!! All formula logic lives in reichardt_kernel.h. This module contains only
!! the Fortran interface declaration and the dispatch wrapper.
!!
module reichardt_device
  use num_types, only : rp, c_rp
  use, intrinsic :: iso_c_binding, only : c_ptr
  use utils, only : neko_error
  implicit none
  private

#ifdef HAVE_HIP
  interface
     subroutine hip_reichardt_compute(u_d, v_d, w_d, &
          ind_r_d, ind_s_d, ind_t_d, ind_e_d, &
          n_x_d, n_y_d, n_z_d, nu_d, h_d, &
          tau_x_d, tau_y_d, tau_z_d, n_nodes, lx, &
          kappa, B, tstep) &
          bind(c, name = 'hip_reichardt_compute')
       use, intrinsic :: iso_c_binding, only : c_ptr, c_int
       use num_types, only : c_rp
       implicit none
       type(c_ptr), value :: u_d, v_d, w_d
       type(c_ptr), value :: ind_r_d, ind_s_d, ind_t_d, ind_e_d
       type(c_ptr), value :: n_x_d, n_y_d, n_z_d, h_d, nu_d
       real(c_rp) :: kappa, B
       type(c_ptr), value :: tau_x_d, tau_y_d, tau_z_d
       integer(c_int) :: n_nodes, lx, tstep
     end subroutine hip_reichardt_compute
  end interface
#elif HAVE_CUDA
  interface
     subroutine cuda_reichardt_compute(u_d, v_d, w_d, &
          ind_r_d, ind_s_d, ind_t_d, ind_e_d, &
          n_x_d, n_y_d, n_z_d, nu_d, h_d, &
          tau_x_d, tau_y_d, tau_z_d, n_nodes, lx, &
          kappa, B, tstep) &
          bind(c, name = 'cuda_reichardt_compute')
       use, intrinsic :: iso_c_binding, only : c_ptr, c_int
       use num_types, only : c_rp
       implicit none
       type(c_ptr), value :: u_d, v_d, w_d
       type(c_ptr), value :: ind_r_d, ind_s_d, ind_t_d, ind_e_d
       type(c_ptr), value :: n_x_d, n_y_d, n_z_d, h_d, nu_d
       real(c_rp) :: kappa, B
       type(c_ptr), value :: tau_x_d, tau_y_d, tau_z_d
       integer(c_int) :: n_nodes, lx, tstep
     end subroutine cuda_reichardt_compute
  end interface
#elif HAVE_OPENCL
#endif
  public :: reichardt_compute_device

contains

  !> Dispatch the Reichardt wall-model computation to the active GPU backend.
  !!
  !! Calls the HIP or CUDA C kernel (defined in reichardt_kernel.h) which
  !! implements the two-term Reichardt (1951) law. OpenCL is not yet supported.
  !!
  !! @param u_d       GPU pointer to velocity x-component.
  !! @param v_d       GPU pointer to velocity y-component.
  !! @param w_d       GPU pointer to velocity z-component.
  !! @param ind_r_d   GPU pointer to r-direction indices.
  !! @param ind_s_d   GPU pointer to s-direction indices.
  !! @param ind_t_d   GPU pointer to t-direction indices.
  !! @param ind_e_d   GPU pointer to element indices.
  !! @param n_x_d     GPU pointer to x-component of wall-normal vector.
  !! @param n_y_d     GPU pointer to y-component of wall-normal vector.
  !! @param n_z_d     GPU pointer to z-component of wall-normal vector.
  !! @param nu_d      GPU pointer to kinematic viscosity.
  !! @param h_d       GPU pointer to wall-normal distance.
  !! @param tau_x_d   GPU pointer to output x-component of wall shear stress.
  !! @param tau_y_d   GPU pointer to output y-component of wall shear stress.
  !! @param tau_z_d   GPU pointer to output z-component of wall shear stress.
  !! @param n_nodes   Number of boundary nodes.
  !! @param lx        Polynomial order in element.
  !! @param kappa     Von Kármán constant.
  !! @param B         Not used in the Reichardt formula; kept for API consistency.
  !! @param tstep     Current time-step (used for initial guess selection).
  subroutine reichardt_compute_device(u_d, v_d, w_d, &
       ind_r_d, ind_s_d, ind_t_d, ind_e_d, &
       n_x_d, n_y_d, n_z_d, nu_d, h_d, tau_x_d, tau_y_d, tau_z_d, &
       n_nodes, lx, kappa, B, tstep)
    integer, intent(in) :: n_nodes, lx, tstep
    type(c_ptr), intent(in) :: u_d, v_d, w_d
    type(c_ptr), intent(in) :: ind_r_d, ind_s_d, ind_t_d, ind_e_d
    type(c_ptr), intent(in) :: n_x_d, n_y_d, n_z_d, h_d, nu_d
    type(c_ptr), intent(inout) :: tau_x_d, tau_y_d, tau_z_d
    real(kind=rp), intent(in) :: kappa, B

#if HAVE_HIP
    call hip_reichardt_compute(u_d, v_d, w_d, &
         ind_r_d, ind_s_d, ind_t_d, ind_e_d, &
         n_x_d, n_y_d, n_z_d, nu_d, h_d, &
         tau_x_d, tau_y_d, tau_z_d, n_nodes, lx, kappa, B, tstep)
#elif HAVE_CUDA
    call cuda_reichardt_compute(u_d, v_d, w_d, &
         ind_r_d, ind_s_d, ind_t_d, ind_e_d, &
         n_x_d, n_y_d, n_z_d, nu_d, h_d, &
         tau_x_d, tau_y_d, tau_z_d, n_nodes, lx, kappa, B, tstep)
#elif HAVE_OPENCL
    call neko_error("OPENCL is not implemented for Reichardt's model")
#else
    call neko_error('No device backend configured')
#endif

  end subroutine reichardt_compute_device

end module reichardt_device