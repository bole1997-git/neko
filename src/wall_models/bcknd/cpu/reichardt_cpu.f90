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
!> Implements the CPU kernel for the `reichardt_t` type.
module reichardt_cpu
  use num_types, only : rp
  use logger, only : neko_log, NEKO_LOG_DEBUG, LOG_SIZE
  implicit none
  private

  public :: reichardt_compute_cpu

  !> Damping length scale in the Reichardt (1951) formula.
  real(kind=rp), parameter :: A_DAMP    = 11.0_rp
  !> Exponential decay scale in the Reichardt (1951) formula.
  real(kind=rp), parameter :: B_EXP     = 3.0_rp
  !> Exponential-term amplitude in the Reichardt (1951) formula.
  real(kind=rp), parameter :: EXP_COEFF = 7.8_rp

contains

  !> Compute the wall shear stress on cpu using Reichardt's model.
  !! @param rho_w The fluid density at the boundary.
  !! @param tstep The current time-step.
  subroutine reichardt_compute_cpu(u, v, w, ind_r, ind_s, ind_t, ind_e, &
       n_x, n_y, n_z, nu, rho_w, h, tau_x, tau_y, tau_z, n_nodes, lx, nelv, &
       kappa, B, tstep)
    integer, intent(in) :: n_nodes, lx, nelv, tstep
    real(kind=rp), dimension(lx, lx, lx, nelv), intent(in) :: u, v, w
    integer, intent(in), dimension(n_nodes) :: ind_r, ind_s, ind_t, ind_e
    real(kind=rp), dimension(n_nodes), intent(in) :: n_x, n_y, n_z, h, nu
    real(kind=rp), dimension(n_nodes), intent(in) :: rho_w
    real(kind=rp), dimension(n_nodes), intent(inout) :: tau_x, tau_y, tau_z
    real(kind=rp), intent(in) :: kappa, B
    integer :: i
    real(kind=rp) :: ui, vi, wi, magu, utau, normu, guess, tau_mag

    do i = 1, n_nodes
       ! Sample the velocity
       ui = u(ind_r(i), ind_s(i), ind_t(i), ind_e(i))
       vi = v(ind_r(i), ind_s(i), ind_t(i), ind_e(i))
       wi = w(ind_r(i), ind_s(i), ind_t(i), ind_e(i))

       ! Project on tangential direction
       normu = ui * n_x(i) + vi * n_y(i) + wi * n_z(i)
       ui = ui - normu * n_x(i)
       vi = vi - normu * n_y(i)
       wi = wi - normu * n_z(i)

       magu = sqrt(ui**2 + vi**2 + wi**2)

       ! Get initial guess for Newton solver
       if (tstep .eq. 1) then
          guess = sqrt(magu * nu(i) / h(i))
       else
          tau_mag = sqrt(tau_x(i)**2 + tau_y(i)**2 + tau_z(i)**2)
          guess   = sqrt(tau_mag / max(rho_w(i), 1.0e-14_rp))
          guess   = max(guess, 1.0e-10_rp)
       end if

       utau = solve_reichardt_cpu(magu, h(i), guess, nu(i), kappa)

       ! Distribute according to the velocity vector
       if (magu > 1.0e-14_rp) then
          tau_x(i) = -rho_w(i) * utau**2 * ui / magu
          tau_y(i) = -rho_w(i) * utau**2 * vi / magu
          tau_z(i) = -rho_w(i) * utau**2 * wi / magu
       else
          tau_x(i) = 0.0_rp
          tau_y(i) = 0.0_rp
          tau_z(i) = 0.0_rp
       end if
    end do

  end subroutine reichardt_compute_cpu

  !> Dimensionless velocity u+ from dimensionless distance y+.
  pure function reichardt_u_plus(y_plus, kappa) result(u_plus)
    real(kind=rp), intent(in) :: y_plus, kappa
    real(kind=rp) :: u_plus
    real(kind=rp) :: log_term, exp_term1, exp_term2

    if (y_plus < 1.0e-12_rp) then
       u_plus = y_plus
       return
    end if

    log_term = (1.0_rp / kappa) * log(1.0_rp + kappa * y_plus)

    exp_term1 = exp(-y_plus / A_DAMP)
    exp_term2 = (y_plus / A_DAMP) * exp(-y_plus / B_EXP)

    u_plus = log_term + EXP_COEFF * (1.0_rp - exp_term1 - exp_term2)

  end function reichardt_u_plus

  !> Derivative du+/dy+, used in the Newton Jacobian.
  pure function reichardt_du_plus_dy(y_plus, kappa) result(du_dy)
    real(kind=rp), intent(in) :: y_plus, kappa
    real(kind=rp) :: du_dy
    real(kind=rp) :: d_log, d_exp, exp_term1, exp_term2

    if (y_plus < 1.0e-12_rp) then
       du_dy = 1.0_rp
       return
    end if

    d_log = 1.0_rp / (1.0_rp + kappa * y_plus)

    exp_term1 = exp(-y_plus / A_DAMP) / A_DAMP
    exp_term2 = (1.0_rp / A_DAMP) * exp(-y_plus / B_EXP) - &
                (y_plus / A_DAMP) * (1.0_rp / B_EXP) * exp(-y_plus / B_EXP)
    d_exp = EXP_COEFF * (exp_term1 - exp_term2)

    du_dy = d_log + d_exp

  end function reichardt_du_plus_dy

  !> Newton solver for the friction velocity using Reichardt's model.
  function solve_reichardt_cpu(u, y, guess, nu, kappa) result(utau)
    real(kind=rp), intent(in) :: u, y, guess, nu, kappa
    real(kind=rp) :: utau
    real(kind=rp) :: y_plus, u_plus, du_dy
    real(kind=rp) :: error, f, df, utau_old
    integer :: k, maxiter
    character(len=LOG_SIZE) :: log_msg

    utau    = guess
    maxiter = 100
    error   = 0.0_rp

    do k = 1, maxiter
       utau_old = utau

       y_plus = y * utau / nu
       u_plus = reichardt_u_plus(y_plus, kappa)
       du_dy  = reichardt_du_plus_dy(y_plus, kappa)

       f  = utau * u_plus - u
       df = u_plus + utau * du_dy * (y / nu)

       if (abs(df) < 1.0e-14_rp) then
          utau = utau * 0.99_rp
          cycle
       end if

       utau = utau - f / df

       if (utau <= 0.0_rp) utau = utau_old * 0.5_rp

       error = abs((utau - utau_old) / (utau + 1.0e-16_rp))
       if (error < 1.0e-8_rp) exit

    end do

    if (k .eq. maxiter + 1) then
       write(log_msg, '(A,E10.3,A,E10.3,A,E10.3)') &
            "Reichardt NC: err=", error, " utau=", utau, " u=", u
       call neko_log%message(log_msg, NEKO_LOG_DEBUG)
    end if

  end function solve_reichardt_cpu

end module reichardt_cpu