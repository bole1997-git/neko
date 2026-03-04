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
!! Uses the original Reichardt (1951) formula as presented in:
!!   Brill, S. (2022). "An Enriched-Basis High-Order Method for 
!!   Wall-Modeled Large Eddy Simulation". Stanford University PhD Thesis.
!!   Equation 2.42
!!
module reichardt_cpu
  use num_types, only : rp
  implicit none
  private

  public :: reichardt_compute_cpu

  ! Reichardt law constants (Original Reichardt 1951)
  real(kind=rp), parameter :: C_REICH = 0.41_rp         ! Log coefficient
  real(kind=rp), parameter :: A_DAMP = 11.0_rp          ! Damping constant
  real(kind=rp), parameter :: B_EXP = 3.0_rp            ! Second exponential constant
  real(kind=rp), parameter :: EXP_COEFF = 7.8_rp        ! Exponential term coefficient

contains

  !> Compute the wall shear stress on CPU using Reichardt's universal law.
  !!
  !! The original Reichardt (1951) law provides an explicit formula covering 
  !! the entire inner layer (viscous sublayer through logarithmic region):
  !!
  !!   u⁺ = y⁺/(1 + y⁺/11) + (1/κ)*ln(1 + 0.41*y⁺) + 7.8*[1 - exp(-y⁺/11) - (y⁺/11)*exp(-y⁺/3)]
  !!
  !! Reference:
  !!   Reichardt, H. (1951). "Vollständige Darstellung der turbulenten 
  !!   Geschwindigkeitsverteilung in Rohren." Zeitschrift für angewandte 
  !!   Mathematik und Mechanik, 31(7-8), 208-219.
  !!
  !! Also presented in:
  !!   Brill, S. (2022). "An Enriched-Basis High-Order Method for 
  !!   Wall-Modeled Large Eddy Simulation." Stanford University PhD Thesis,
  !!   Equation 2.42.
  !!
  !! @param u Velocity component in x-direction.
  !! @param v Velocity component in y-direction.
  !! @param w Velocity component in z-direction.
  !! @param ind_r Radial indices of sampling points.
  !! @param ind_s S-direction indices of sampling points.
  !! @param ind_t T-direction indices of sampling points.
  !! @param ind_e Element indices of sampling points.
  !! @param n_x X-component of wall-normal vector.
  !! @param n_y Y-component of wall-normal vector.
  !! @param n_z Z-component of wall-normal vector.
  !! @param nu Kinematic viscosity at boundary nodes.
  !! @param h Distance from wall to sampling point.
  !! @param tau_x Output: X-component of wall shear stress.
  !! @param tau_y Output: Y-component of wall shear stress.
  !! @param tau_z Output: Z-component of wall shear stress.
  !! @param n_nodes Number of boundary nodes.
  !! @param lx Polynomial order in element.
  !! @param nelv Number of elements.
  !! @param kappa Von Kármán constant.
  !! @param B Log-law intercept (not used in Reichardt formula).
  !! @param tstep Current time-step.
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

    do i=1, n_nodes
       ! Sample the velocity at the off-wall point
       ui = u(ind_r(i), ind_s(i), ind_t(i), ind_e(i))
       vi = v(ind_r(i), ind_s(i), ind_t(i), ind_e(i))
       wi = w(ind_r(i), ind_s(i), ind_t(i), ind_e(i))

       ! Project velocity onto tangential plane (remove normal component)
       ! Wall-normal component
       normu = ui * n_x(i) + vi * n_y(i) + wi * n_z(i)

       ! Tangential velocity components
       ui = ui - normu * n_x(i)
       vi = vi - normu * n_y(i)
       wi = wi - normu * n_z(i)

       ! Magnitude of tangential velocity
       magu = sqrt(ui**2 + vi**2 + wi**2)

       ! Get initial guess for Newton solver
       if (tstep .eq. 1) then
          ! First time-step: use simple estimate
          guess = sqrt(magu * nu(i) / h(i))
       else
          ! Use previous solution as starting point
          guess = sqrt(sqrt(tau_x(i)**2 + tau_y(i)**2 + tau_z(i)**2))
       end if

       ! Solve for friction velocity u_tau
       utau = solve_reichardt_cpu(magu, h(i), guess, nu(i), kappa)

       ! Distribute shear stress according to velocity direction
       ! tau_wall = -ρ * u_tau² * (u_tangential / |u_tangential|)
       if (magu > 1.0e-14_rp) then
          tau_x(i) = -utau**2 * ui / magu
          tau_y(i) = -utau**2 * vi / magu
          tau_z(i) = -utau**2 * wi / magu
       else
          ! Avoid division by zero when velocity is zero
          tau_x(i) = 0.0_rp
          tau_y(i) = 0.0_rp
          tau_z(i) = 0.0_rp
       end if
    end do

  end subroutine reichardt_compute_cpu

  !> Dimensionless velocity u⁺ from dimensionless distance y⁺ 
  !! using the original Reichardt (1951) formula.
  !!
  !! Formula (from Reichardt 1951, Equation 2.42 in Brill thesis):
  !!   u⁺ = y⁺/(1 + y⁺/11) + (1/κ)*ln(1 + 0.41*y⁺) + 7.8*[1 - exp(-y⁺/11) - (y⁺/11)*exp(-y⁺/3)]
  !!
  !! This formula smoothly covers the entire inner layer without piecewise definitions.
  !!
  !! @param y_plus Dimensionless distance y⁺ = y*u_τ/ν.
  !! @param kappa Von Kármán constant (typically 0.41).
  !! @return Dimensionless velocity u⁺.
  pure function reichardt_u_plus(y_plus, kappa) result(u_plus)
    real(kind=rp), intent(in) :: y_plus, kappa
    real(kind=rp) :: u_plus
    real(kind=rp) :: linear_term, log_term, exp_term1, exp_term2, log_arg

    ! Avoid computation for very small y⁺ (linear region: u⁺ ≈ y⁺)
    if (y_plus < 1.0e-12_rp) then
       u_plus = y_plus
       return
    end if

    ! Linear term with damping: y⁺/(1 + y⁺/11)
    linear_term = y_plus / (1.0_rp + y_plus / A_DAMP)

    ! Logarithmic term: (1/κ)*ln(1 + 0.41*y⁺)
    log_arg = 1.0_rp + C_REICH * y_plus
    if (log_arg <= 0.0_rp) then
       ! Safeguard: use linear term only
       u_plus = linear_term
       return
    end if
    log_term = (1.0_rp / kappa) * log(log_arg)

    ! Exponential correction term: 7.8*[1 - exp(-y⁺/11) - (y⁺/11)*exp(-y⁺/3)]
    exp_term1 = exp(-y_plus / A_DAMP)
    exp_term2 = (y_plus / A_DAMP) * exp(-y_plus / B_EXP)

    ! Total dimensionless velocity
    u_plus = linear_term + log_term + EXP_COEFF * (1.0_rp - exp_term1 - exp_term2)

  end function reichardt_u_plus

  !> Derivative du⁺/dy⁺ of the original Reichardt (1951) formula.
  !!
  !! Used for Newton-Raphson iteration. Computed analytically from:
  !!   u⁺ = y⁺/(1 + y⁺/11) + (1/κ)*ln(1 + 0.41*y⁺) + 7.8*[1 - exp(-y⁺/11) - (y⁺/11)*exp(-y⁺/3)]
  !!
  !! @param y_plus Dimensionless distance y⁺.
  !! @param kappa Von Kármán constant.
  !! @return Derivative du⁺/dy⁺.
  pure function reichardt_du_plus_dy(y_plus, kappa) result(du_dy)
    real(kind=rp), intent(in) :: y_plus, kappa
    real(kind=rp) :: du_dy
    real(kind=rp) :: d_linear, d_log, d_exp, exp_term1, exp_term2, log_arg

    if (y_plus < 1.0e-12_rp) then
       du_dy = 1.0_rp  ! Linear region: du⁺/dy⁺ = 1
       return
    end if

    ! Derivative of linear term: d/dy[y/(1 + y/11)] = 1/(1 + y/11)²
    d_linear = 1.0_rp / (1.0_rp + y_plus / A_DAMP)**2

    ! Derivative of log term: d/dy[(1/κ)*ln(1 + 0.41*y)]
    log_arg = 1.0_rp + C_REICH * y_plus
    if (log_arg <= 0.0_rp) then
       du_dy = d_linear
       return
    end if
    d_log = (1.0_rp / kappa) * C_REICH / log_arg

    ! Derivative of exponential term:
    ! d/dy[7.8*(1 - exp(-y/11) - (y/11)*exp(-y/3))]
    ! = 7.8*[exp(-y/11)/11 - (1/11)*exp(-y/3) + (y/11)*(1/3)*exp(-y/3)]
    exp_term1 = exp(-y_plus / A_DAMP) / A_DAMP
    exp_term2 = (1.0_rp / A_DAMP) * exp(-y_plus / B_EXP) - &
                (y_plus / A_DAMP) * (1.0_rp / B_EXP) * exp(-y_plus / B_EXP)
    d_exp = EXP_COEFF * (exp_term1 + exp_term2)

    du_dy = d_linear + d_log + d_exp

  end function reichardt_du_plus_dy

  !> Newton-Raphson solver for friction velocity using original Reichardt law.
  !!
  !! Solves the implicit equation:
  !!   F(u_τ) = u_τ * u⁺(y⁺) - U = 0
  !! where y⁺ = y * u_τ / ν and U is the tangential velocity magnitude.
  !!
  !! The Newton update is:
  !!   u_τ^(n+1) = u_τ^(n) - F(u_τ^(n)) / F'(u_τ^(n))
  !!
  !! where F'(u_τ) = u⁺ + u_τ * du⁺/dy⁺ * (y/ν)
  !!
  !! @param u Tangential velocity magnitude.
  !! @param y Wall-normal distance.
  !! @param guess Initial guess for u_τ.
  !! @param nu Kinematic viscosity.
  !! @param kappa Von Kármán constant.
  !! @return Friction velocity u_τ.
  function solve_reichardt_cpu(u, y, guess, nu, kappa) result(utau)
    real(kind=rp), intent(in) :: u
    real(kind=rp), intent(in) :: y
    real(kind=rp), intent(in) :: guess
    real(kind=rp), intent(in) :: nu, kappa
    real(kind=rp) :: y_plus, u_plus, du_dy
    real(kind=rp) :: error, f, df, utau_old, utau
    integer :: k, maxiter, niter

    utau = guess

    ! Newton-Raphson iteration parameters
    maxiter = 100

    do k = 1, maxiter
       utau_old = utau

       ! Current dimensionless coordinates
       y_plus = y * utau / nu
       u_plus = reichardt_u_plus(y_plus, kappa)
       du_dy = reichardt_du_plus_dy(y_plus, kappa)

       ! Residual: F(u_τ) = u_τ * u⁺ - U
       f = utau * u_plus - u

       ! Sensitivity: dF/du_τ = u⁺ + u_τ * du⁺/dy⁺ * dy⁺/du_τ
       ! where dy⁺/du_τ = y / ν
       df = u_plus + utau * du_dy * (y / nu)

       ! Safeguard against singular Jacobian
       if (abs(df) < 1.0e-14_rp) then
          ! Jacobian is nearly zero, use simpler update
          utau = utau * 0.99_rp
          cycle
       end if

       ! Newton step
       utau = utau - f / df

       ! Safeguard: prevent u_τ from becoming negative or zero
       if (utau <= 0.0_rp) then
          utau = utau_old * 0.5_rp
       end if

       ! Relative error check for convergence
       error = abs((utau - utau_old) / (utau + 1.0e-16_rp))
       niter = k

       ! Convergence tolerance
       if (error < 1.0e-8_rp) then
          exit
       end if

    end do

  end function solve_reichardt_cpu

end module reichardt_cpu