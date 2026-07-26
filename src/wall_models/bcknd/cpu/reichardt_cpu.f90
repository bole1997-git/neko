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
!!
!! Uses the original Reichardt (1951) two-term formula:
!!
!!   u+ = (1/kappa)*ln(1 + kappa*y+)
!!      + 7.8*[1 - exp(-y+/11) - (y+/11)*exp(-y+/3)]
!!
!! Reference:
!!   Reichardt, H. (1951). "Vollstandige Darstellung der turbulenten
!!   Geschwindigkeitsverteilung in Rohren." Zeitschrift fur angewandte
!!   Mathematik und Mechanik, 31(7-8), 208-219.
!!
module reichardt_cpu
  use num_types, only : rp
  use logger, only : neko_log, NEKO_LOG_DEBUG, LOG_SIZE
  implicit none
  private

  public :: reichardt_compute_cpu

  ! Reichardt (1951) constants
  real(kind=rp), parameter :: A_DAMP    = 11.0_rp   ! Damping length scale
  real(kind=rp), parameter :: B_EXP     = 3.0_rp    ! Exponential decay scale
  real(kind=rp), parameter :: EXP_COEFF = 7.8_rp    ! Exponential amplitude

contains

  !> Compute the wall shear stress on CPU using the original Reichardt (1951) law.
  !!
  !! @param u,v,w      Velocity components on the full (lx,lx,lx,nelv) mesh.
  !! @param ind_r/s/t/e Off-wall sampling-point indices into the 4D array.
  !! @param n_x/y/z    Wall-normal unit vector at each boundary node.
  !! @param nu         Kinematic viscosity at boundary nodes.
  !! @param h          Wall-normal distance to the sampling point.
  !! @param tau_x/y/z  Wall shear stress components (inout).
  !! @param n_nodes    Number of boundary nodes.
  !! @param lx         GLL polynomial order.
  !! @param nelv       Number of elements.
  !! @param kappa      Von Karman constant.
  !! @param B          Log-law intercept (not used; kept for API compatibility).
  !! @param tstep      Current time step.
  subroutine reichardt_compute_cpu(u, v, w, ind_r, ind_s, ind_t, ind_e, &
       n_x, n_y, n_z, nu, h, tau_x, tau_y, tau_z, n_nodes, lx, nelv, &
       kappa, B, tstep)
    integer, intent(in) :: n_nodes, lx, nelv, tstep
    real(kind=rp), dimension(lx, lx, lx, nelv), intent(in) :: u, v, w
    integer, intent(in), dimension(n_nodes) :: ind_r, ind_s, ind_t, ind_e
    real(kind=rp), dimension(n_nodes), intent(in) :: n_x, n_y, n_z, h, nu
    real(kind=rp), dimension(n_nodes), intent(inout) :: tau_x, tau_y, tau_z
    real(kind=rp), intent(in) :: kappa, B
    integer :: i
    real(kind=rp) :: ui, vi, wi, magu, utau, normu, guess

    do i = 1, n_nodes
       ! Sample velocity at the off-wall point.
       ui = u(ind_r(i), ind_s(i), ind_t(i), ind_e(i))
       vi = v(ind_r(i), ind_s(i), ind_t(i), ind_e(i))
       wi = w(ind_r(i), ind_s(i), ind_t(i), ind_e(i))

       ! Remove wall-normal component to get tangential velocity.
       normu = ui * n_x(i) + vi * n_y(i) + wi * n_z(i)
       ui = ui - normu * n_x(i)
       vi = vi - normu * n_y(i)
       wi = wi - normu * n_z(i)

       magu = sqrt(ui**2 + vi**2 + wi**2)

       ! Initial guess for the Newton solver.
       if (tstep .eq. 1) then
          ! First timestep: laminar sublayer estimate u_tau ~ sqrt(U*nu/y).
          guess = sqrt(magu * nu(i) / h(i))
       else
          ! Warm start from previous shear stress magnitude.
          guess = sqrt(sqrt(tau_x(i)**2 + tau_y(i)**2 + tau_z(i)**2))
       end if

       ! Solve for friction velocity u_tau via Newton-Raphson.
       utau = solve_reichardt_cpu(magu, h(i), guess, nu(i), kappa)

       ! Distribute shear stress in the tangential velocity direction:
       !   tau_wall = -u_tau^2 * (u_tang / |u_tang|)
       if (magu > 1.0e-14_rp) then
          tau_x(i) = -utau**2 * ui / magu
          tau_y(i) = -utau**2 * vi / magu
          tau_z(i) = -utau**2 * wi / magu
       else
          tau_x(i) = 0.0_rp
          tau_y(i) = 0.0_rp
          tau_z(i) = 0.0_rp
       end if
    end do

  end subroutine reichardt_compute_cpu

  !> Dimensionless velocity u+ from dimensionless distance y+
  !! using the original two-term Reichardt (1951) formula.
  !!
  !!   u+ = (1/kappa)*ln(1 + kappa*y+)
  !!      + 7.8*[1 - exp(-y+/11) - (y+/11)*exp(-y+/3)]
  !!
  !! Limiting behaviour:
  !!   y+ -> 0:  u+ -> y+          (viscous sublayer)
  !!   y+ -> inf: u+ -> (1/kappa)*ln(y+) + B   (log law)
  pure function reichardt_u_plus(y_plus, kappa) result(u_plus)
    real(kind=rp), intent(in) :: y_plus, kappa
    real(kind=rp) :: u_plus
    real(kind=rp) :: log_term, exp_term1, exp_term2

    if (y_plus < 1.0e-12_rp) then
       u_plus = y_plus
       return
    end if

    ! Logarithmic term: (1/kappa)*ln(1 + kappa*y+)
    log_term = (1.0_rp / kappa) * log(1.0_rp + kappa * y_plus)

    ! Exponential correction: 7.8*[1 - exp(-y+/11) - (y+/11)*exp(-y+/3)]
    exp_term1 = exp(-y_plus / A_DAMP)
    exp_term2 = (y_plus / A_DAMP) * exp(-y_plus / B_EXP)

    u_plus = log_term + EXP_COEFF * (1.0_rp - exp_term1 - exp_term2)

  end function reichardt_u_plus

  !> Analytical derivative du+/dy+ of the Reichardt (1951) formula.
  !! Used in the Newton-Raphson iteration via the chain rule:
  !!   dF/du_tau = u+(y+) + u_tau * (du+/dy+) * (y/nu)
  pure function reichardt_du_plus_dy(y_plus, kappa) result(du_dy)
    real(kind=rp), intent(in) :: y_plus, kappa
    real(kind=rp) :: du_dy
    real(kind=rp) :: d_log, d_exp, exp_term1, exp_term2

    if (y_plus < 1.0e-12_rp) then
       du_dy = 1.0_rp
       return
    end if

    ! Derivative of log term: 1/(1 + kappa*y+)
    d_log = 1.0_rp / (1.0_rp + kappa * y_plus)

    ! Derivative of exponential correction:
    !   d/dy[7.8*(1 - exp(-y/11) - (y/11)*exp(-y/3))]
    ! = 7.8*[ (1/11)*exp(-y/11)  -  (1/11)*exp(-y/3)  +  (y/33)*exp(-y/3) ]
    exp_term1 = exp(-y_plus / A_DAMP) / A_DAMP
    exp_term2 = (1.0_rp / A_DAMP) * exp(-y_plus / B_EXP) - &
                (y_plus / A_DAMP) * (1.0_rp / B_EXP) * exp(-y_plus / B_EXP)
    d_exp = EXP_COEFF * (exp_term1 - exp_term2)

    du_dy = d_log + d_exp

  end function reichardt_du_plus_dy

  !> Newton-Raphson solver for friction velocity u_tau using the Reichardt law.
  !!
  !! Solves: F(u_tau) = u_tau * u+(y+) - U = 0
  !! Jacobian: dF/du_tau = u+ + u_tau * (du+/dy+) * (y/nu)
  !!
  !! Non-convergence is logged at DEBUG level, matching the Spalding standard.
  !!
  !! @param u     Tangential velocity magnitude U.
  !! @param y     Wall-normal distance to the sampling point.
  !! @param guess Initial guess for u_tau.
  !! @param nu    Kinematic viscosity.
  !! @param kappa Von Karman constant.
  !! @return      Friction velocity u_tau.
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