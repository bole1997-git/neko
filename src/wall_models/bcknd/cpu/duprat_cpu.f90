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
!> Implements the CPU kernel for the `duprat_t` type.
!!
!! Uses the Duprat et al. (2011) extended law of the wall with
!! pressure gradient correction via Simpson velocity scale.
!!
!! Reference:
!!   Duprat, C., et al. (2011). "A wall-layer model for large-eddy 
!!   simulations of turbulent flows with/without pressure gradient."
!!   Physics of Fluids, 23(015101).
!!
module duprat_cpu
  use num_types, only: rp
  implicit none
  private

  public :: duprat_compute_cpu

  ! Duprat constants
  real(kind=rp), parameter :: YPLUS_THRESHOLD = 5.0_rp  ! Log region threshold

contains

  !> Compute the wall shear stress on CPU using the Duprat (2011) law.
  !!
  !! Handles both zero pressure gradient (ZPG) and adverse pressure gradient (APG)
  !! through extended scaling with Simpson velocity scale.
  !!
  !! @param u Velocity component in x-direction (4D array).
  !! @param v Velocity component in y-direction (4D array).
  !! @param w Velocity component in z-direction (4D array).
  !! @param ind_r Radial indices of sampling points.
  !! @param ind_s S-direction indices of sampling points.
  !! @param ind_t T-direction indices of sampling points.
  !! @param ind_e Element indices of sampling points.
  !! @param n_x X-component of wall-normal vector.
  !! @param n_y Y-component of wall-normal vector.
  !! @param n_z Z-component of wall-normal vector.
  !! @param nu Kinematic viscosity at boundary nodes.
  !! @param dpdx Pressure gradient (streamwise) at boundary nodes.
  !! @param h Distance from wall to sampling point.
  !! @param tau_x Output: X-component of wall shear stress.
  !! @param tau_y Output: Y-component of wall shear stress.
  !! @param tau_z Output: Z-component of wall shear stress.
  !! @param n_nodes Number of boundary nodes.
  !! @param lx Polynomial order in element.
  !! @param nelv Number of elements.
  !! @param kappa Von Kármán constant.
  !! @param beta Duprat damping coefficient.
  !! @param A Van Driest constant.
  !! @param tstep Current time-step.
  subroutine duprat_compute_cpu(u, v, w, ind_r, ind_s, ind_t, ind_e, &
       n_x, n_y, n_z, nu, dpdx, h, tau_x, tau_y, tau_z, n_nodes, lx, nelv, &
       kappa, beta, A, tstep)
    integer, intent(in) :: n_nodes, lx, nelv, tstep
    real(kind=rp), dimension(lx, lx, lx, nelv), intent(in) :: u, v, w
    integer, intent(in), dimension(n_nodes) :: ind_r, ind_s, ind_t, ind_e
    real(kind=rp), dimension(n_nodes), intent(in) :: n_x, n_y, n_z, nu, dpdx, h
    real(kind=rp), dimension(n_nodes), intent(inout) :: tau_x, tau_y, tau_z
    real(kind=rp), intent(in) :: kappa, beta, A
    integer :: i
    real(kind=rp) :: ui, vi, wi, magu, utau, normu, guess
    real(kind=rp) :: rho

    ! Normalized density (unit density for this problem)
    rho = 1.0_rp

    do i = 1, n_nodes
       ! Sample the velocity at the off-wall point
       ui = u(ind_r(i), ind_s(i), ind_t(i), ind_e(i))
       vi = v(ind_r(i), ind_s(i), ind_t(i), ind_e(i))
       wi = w(ind_r(i), ind_s(i), ind_t(i), ind_e(i))

       ! Remove wall-normal component to get tangential velocity
       normu = ui * n_x(i) + vi * n_y(i) + wi * n_z(i)
       ui = ui - normu * n_x(i)
       vi = vi - normu * n_y(i)
       wi = wi - normu * n_z(i)

       ! Magnitude of tangential velocity
       magu = sqrt(ui**2 + vi**2 + wi**2)

       ! Initial guess for the Newton solver
       if (tstep .eq. 1) then
          ! First time-step: laminar sublayer estimate u_tau ~ sqrt(U*nu/y)
          guess = sqrt(magu * nu(i) / h(i))
       else
          ! Warm start from previous shear stress magnitude
          guess = sqrt(sqrt(tau_x(i)**2 + tau_y(i)**2 + tau_z(i)**2))
       end if

       ! Solve for friction velocity u_tau via Newton-Raphson
       utau = solve_duprat_cpu(magu, h(i), guess, nu(i), dpdx(i), &
                               kappa, beta, A, rho)

       ! Distribute shear stress in the tangential velocity direction
       ! tau_wall = -u_tau² * (u_tangential / |u_tangential|)
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

  end subroutine duprat_compute_cpu

  !> Newton-Raphson solver for friction velocity u_τ using the Duprat law.
  !!
  !! Solves the implicit equation:
  !!   F(u_τ) = u_τ * u⁺(y⁺) - U = 0
  !! where y⁺ = y * u_τ / ν and u⁺ follows the Duprat extended scaling law.
  !!
  !! @param u Tangential velocity magnitude U.
  !! @param y Wall-normal distance to the sampling point.
  !! @param guess Initial guess for u_τ.
  !! @param nu Kinematic viscosity.
  !! @param dpdx Pressure gradient (streamwise component).
  !! @param kappa Von Kármán constant.
  !! @param beta Duprat damping coefficient.
  !! @param A Van Driest constant.
  !! @param rho Density (normalized).
  !! @return Friction velocity u_τ.
  function solve_duprat_cpu(u, y, guess, nu, dpdx, kappa, beta, A, rho) &
       result(utau)
    real(kind=rp), intent(in) :: u, y, guess, nu, dpdx
    real(kind=rp), intent(in) :: kappa, beta, A, rho
    real(kind=rp) :: utau
    
    real(kind=rp) :: y_plus, u_plus, du_dy
    real(kind=rp) :: error, f, df, utau_old
    real(kind=rp) :: u_p, alpha, sign_dpdx_val
    integer :: k

    utau = guess

    do k = 1, 100
       utau_old = utau

       ! Compute Simpson velocity scale
       if (abs(dpdx) < 1.0e-10_rp) then
          ! ZPG case
          u_p = 1.0_rp
          alpha = 1.0_rp
          sign_dpdx_val = 0.0_rp
       else
          ! APG case
          u_p = (rho / 2.0_rp * abs(dpdx))**(1.0_rp / 3.0_rp)
          alpha = utau**2 / (u_p**2 + 1.0e-15_rp)
          alpha = max(0.0_rp, min(1.0_rp, alpha))  ! Clamp to [0, 1]
          sign_dpdx_val = sign(1.0_rp, dpdx)
       end if

       ! Dimensionless coordinates at current iterate
       y_plus = y * utau / nu
       u_plus = duprat_u_plus(y_plus, kappa, alpha, A, sign_dpdx_val)
       du_dy = duprat_du_plus_dy(y_plus, kappa, alpha, A, sign_dpdx_val)

       ! Residual F and Jacobian dF/du_τ
       f = utau * u_plus - u
       df = u_plus + utau * du_dy * (y / nu)

       ! Guard against near-zero Jacobian
       if (abs(df) < 1.0e-14_rp) then
          utau = utau * 0.99_rp
          cycle
       end if

       ! Newton step
       utau = utau - f / df

       ! Keep u_τ strictly positive
       if (utau <= 0.0_rp) then
          utau = utau_old * 0.5_rp
       end if

       ! Convergence check on relative change
       error = abs((utau - utau_old) / (utau + 1.0e-16_rp))
       if (error < 1.0e-8_rp) exit

    end do

  end function solve_duprat_cpu

  !> Dimensionless velocity u⁺ from dimensionless distance y⁺ using Duprat law.
  !!
  !! Extended scaling law with pressure gradient correction:
  !!   U* = U / u_p,  y* = y * u_p / ν
  !!   α = u_τ² / u_p² (0 = separation, 1 = ZPG)
  !!
  !! Velocity profile:
  !!   χ = κ*y* + y*(1-α)^1.5 * [1 - exp(-y*/(1+A*α³))]²
  !!   u* = χ + (1/κ)*ln(1 + κ*χ) * [1 - exp(-y*/26)]  (approximate)
  !!
  !! @param y_star Dimensionless wall distance in extended scaling.
  !! @param kappa Von Kármán constant.
  !! @param alpha Pressure gradient parameter.
  !! @param A Van Driest constant.
  !! @param sign_dpdx Sign of pressure gradient.
  !! @return Dimensionless velocity u⁺.
  pure function duprat_u_plus(y_star, kappa, alpha, A, sign_dpdx) result(u_plus)
    real(kind=rp), intent(in) :: y_star, kappa, alpha, A, sign_dpdx
    real(kind=rp) :: u_plus
    real(kind=rp) :: chi, viscous_part, log_part, damp_factor, log_arg
    real(kind=rp) :: alpha_clamped

    ! Clamp alpha to [0, 1]
    alpha_clamped = max(0.0_rp, min(1.0_rp, alpha))

    ! Avoid numerical issues at very small y_star
    if (y_star < 1.0e-12_rp) then
       u_plus = y_star
       return
    end if

    ! Van Driest damping with PG correction
    chi = kappa * y_star + y_star * (1.0_rp - alpha_clamped)**1.5_rp * &
          (1.0_rp - exp(-y_star / (1.0_rp + A * alpha_clamped**3)))**2

    viscous_part = chi

    ! Logarithmic contribution
    if (y_star > YPLUS_THRESHOLD) then
       log_arg = 1.0_rp + kappa * chi
       if (log_arg > 0.0_rp) then
          log_part = (1.0_rp / kappa) * log(log_arg)
       else
          log_part = viscous_part
       end if
    else
       log_part = viscous_part
    end if

    ! Damping function
    damp_factor = 1.0_rp - exp(-y_star / 26.0_rp)

    u_plus = viscous_part + (log_part - viscous_part) * damp_factor

  end function duprat_u_plus

  !> Analytical derivative du⁺/dy⁺ of the Duprat velocity profile.
  !!
  !! Computed via the chain rule through the profile equation.
  !!
  !! @param y_star Dimensionless wall distance.
  !! @param kappa Von Kármán constant.
  !! @param alpha Pressure gradient parameter.
  !! @param A Van Driest constant.
  !! @param sign_dpdx Sign of pressure gradient.
  !! @return Derivative du⁺/dy⁺.
  pure function duprat_du_plus_dy(y_star, kappa, alpha, A, sign_dpdx) &
       result(du_dy)
    real(kind=rp), intent(in) :: y_star, kappa, alpha, A, sign_dpdx
    real(kind=rp) :: du_dy
    real(kind=rp) :: chi, dchi_dy, d_damp, log_arg, d_log
    real(kind=rp) :: viscous_part, log_part, damp_factor
    real(kind=rp) :: d_viscous, d_log_part, d_damp_dy
    real(kind=rp) :: alpha_clamped, exp_term, A_factor

    ! Clamp alpha
    alpha_clamped = max(0.0_rp, min(1.0_rp, alpha))

    ! Limiting case
    if (y_star < 1.0e-12_rp) then
       du_dy = 1.0_rp
       return
    end if

    ! Compute chi and its derivative
    A_factor = 1.0_rp + A * alpha_clamped**3
    exp_term = exp(-y_star / A_factor)
    
    chi = kappa * y_star + y_star * (1.0_rp - alpha_clamped)**1.5_rp * &
          (1.0_rp - exp_term)**2

    ! dχ/dy*
    dchi_dy = kappa + (1.0_rp - alpha_clamped)**1.5_rp * &
              ((1.0_rp - exp_term)**2 + &
               y_star * 2.0_rp * (1.0_rp - exp_term) * &
               (1.0_rp / A_factor) * exp_term)

    viscous_part = chi
    d_viscous = dchi_dy

    ! Log part and its derivative
    if (y_star > YPLUS_THRESHOLD) then
       log_arg = 1.0_rp + kappa * chi
       if (log_arg > 0.0_rp) then
          log_part = (1.0_rp / kappa) * log(log_arg)
          d_log = dchi_dy / log_arg
       else
          log_part = viscous_part
          d_log = d_viscous
       end if
       d_log_part = d_log
    else
       log_part = viscous_part
       d_log_part = d_viscous
    end if

    ! Damping function and its derivative
    damp_factor = 1.0_rp - exp(-y_star / 26.0_rp)
    d_damp_dy = (1.0_rp / 26.0_rp) * exp(-y_star / 26.0_rp)

    ! Final derivative: d(u⁺)/dy⁺
    du_dy = d_viscous + &
            (d_log_part - d_viscous) * damp_factor + &
            (log_part - viscous_part) * d_damp_dy

  end function duprat_du_plus_dy

end module duprat_cpu