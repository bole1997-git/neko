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
!> Implements the GPU device kernel for the `reichardt_t` type.
module reichardt_device
  use num_types, only : rp
  implicit none
  private

  public :: reichardt_compute_device

contains

  !> Placeholder for GPU device kernel (not yet implemented).
  !!
  !! This subroutine is a placeholder for future GPU implementation
  !! using CUDA or HIP. For now, it simply returns an error message.
  !!
  subroutine reichardt_compute_device(u, v, w, ind_r, ind_s, ind_t, ind_e, &
       n_x, n_y, n_z, nu, h, tau_x, tau_y, tau_z, n_nodes, lx, &
       kappa, B, tstep)
    real(kind=rp), intent(in) :: u(:,:,:,:), v(:,:,:,:), w(:,:,:,:)
    integer, intent(in) :: ind_r(:), ind_s(:), ind_t(:), ind_e(:)
    real(kind=rp), intent(in) :: n_x(:), n_y(:), n_z(:)
    real(kind=rp), intent(in) :: nu(:), h(:)
    real(kind=rp), intent(inout) :: tau_x(:), tau_y(:), tau_z(:)
    integer, intent(in) :: n_nodes, lx
    real(kind=rp), intent(in) :: kappa, B
    integer, intent(in) :: tstep

    ! Placeholder: Device kernel not yet implemented
    ! This would contain CUDA/HIP code for GPU acceleration
    error stop "Reichardt GPU kernel not yet implemented"

  end subroutine reichardt_compute_device

end module reichardt_device