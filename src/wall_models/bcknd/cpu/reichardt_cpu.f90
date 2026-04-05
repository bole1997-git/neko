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
!!   u⁺ = (1/κ)*ln(1 + κ*y⁺) + 7.8*[1 - exp(-y⁺/11) - (y⁺/11)*exp(-y⁺/3)]
!!
!! Reference:
!!   Reichardt, H. (1951). "Vollständige Darstellung der turbulenten
!!   Geschwindigkeitsverteilung in Rohren." Zeitschrift für angewandte
!!   Mathematik und Mechanik, 31(7-8), 208-219.
!!
!! Note: Some secondary sources (e.g. Brill 2022 eq. 2.42) add a spurious
!! linear term y⁺/(1 + y⁺/11) that is NOT part of the original Reichardt
!! formula and causes significant over-prediction of u⁺ in the log region.
!!
module reichardt_cpu
  use num_types, only : rp
  implicit none
  private

  public :: reichardt_compute_cpu

  ! Reichardt (1951) constants
  real(kind=rp), parameter :: A_DAMP   = 11.0_rp   ! Damping length scale
  real(kind=rp), parameter :: B_EXP    = 3.0_rp    ! Exponential decay scale
  real(kind=rp), parameter :: EXP_COEFF = 7.8_rp   ! Exponential amplitude

contains

  !> Compute the wall shear stress on CPU using the original Reichardt (1951) law.
  !!
  !! The two-term formula covers the entire inner layer continuously:
  !!
  !!   u⁺ = (1/κ)*ln(1 + κ*y⁺) + 7.8*[1 - exp(-y⁺/11) - (y⁺/11)*exp(-y⁺/3)]
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
  !! @param B Log-law intercept (not used in Reichardt formula, kept for API compatibility).
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
       utau = solve_reichardt_cpu(magu, h(i), guess, nu(i), kappa)

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

  end subroutine reichardt_compute_cpu

  !> Dimensionless velocity u⁺ from dimensionless distance y⁺
  !! using the original two-term Reichardt (1951) formula.
  !!
  !! Formula:
  !!   u⁺ = (1/κ)*ln(1 + κ*y⁺) + 7.8*[1 - exp(-y⁺/11) - (y⁺/11)*exp(-y⁺/3)]
  !!
  !! Limiting behaviour:
  !!   - As y⁺ → 0:  log_term → y⁺,  exp correction → 0,  so u⁺ → y⁺  (viscous sublayer)
  !!   - As y⁺ → ∞:  recovers u⁺ = (1/κ)*ln(y⁺) + B  (log law)
  !!
  !! @param y_plus Dimensionless wall distance y⁺ = y*u_τ/ν.
  !! @param kappa Von Kármán constant (typically 0.41).
  !! @return Dimensionless velocity u⁺.
  pure function reichardt_u_plus(y_plus, kappa) result(u_plus)
    real(kind=rp), intent(in) :: y_plus, kappa
    real(kind=rp) :: u_plus
    real(kind=rp) :: log_term, exp_term1, exp_term2, log_arg

    ! For very small y⁺ the formula reduces analytically to u⁺ ≈ y⁺
    if (y_plus < 1.0e-12_rp) then
       u_plus = y_plus
       return
    end if

    ! Logarithmic term: (1/κ)*ln(1 + κ*y⁺)
    log_arg = 1.0_rp + kappa * y_plus
    log_term = (1.0_rp / kappa) * log(log_arg)

    ! Exponential correction: 7.8*[1 - exp(-y⁺/11) - (y⁺/11)*exp(-y⁺/3)]
    exp_term1 = exp(-y_plus / A_DAMP)
    exp_term2 = (y_plus / A_DAMP) * exp(-y_plus / B_EXP)

    u_plus = log_term + EXP_COEFF * (1.0_rp - exp_term1 - exp_term2)

  end function reichardt_u_plus

  !> Analytical derivative du⁺/dy⁺ of the original Reichardt (1951) formula.
  !!
  !! Derived from:
  !!   u⁺ = (1/κ)*ln(1 + κ*y⁺) + 7.8*[1 - exp(-y⁺/11) - (y⁺/11)*exp(-y⁺/3)]
  !!
  !! Giving:
  !!   du⁺/dy⁺ = κ / (κ*(1 + κ*y⁺))
  !!            + 7.8*[ exp(-y⁺/11)/11  -  (1/11)*exp(-y⁺/3)  +  (y⁺/33)*exp(-y⁺/3) ]
  !!
  !! Used in the Newton-Raphson iteration via the chain rule:
  !!   dF/du_τ = u⁺(y⁺) + u_τ * (du⁺/dy⁺) * (y/ν)
  !!
  !! @param y_plus Dimensionless wall distance y⁺.
  !! @param kappa Von Kármán constant.
  !! @return Derivative du⁺/dy⁺.
  pure function reichardt_du_plus_dy(y_plus, kappa) result(du_dy)
    real(kind=rp), intent(in) :: y_plus, kappa
    real(kind=rp) :: du_dy
    real(kind=rp) :: d_log, d_exp, exp_term1, exp_term2, log_arg

    ! Limiting case: du⁺/dy⁺ → 1 as y⁺ → 0 (viscous sublayer slope)
    if (y_plus < 1.0e-12_rp) then
       du_dy = 1.0_rp
       return
    end if

    ! Derivative of log term: d/dy[(1/κ)*ln(1 + κ*y⁺)] = 1/(1 + κ*y⁺)
    log_arg = 1.0_rp + kappa * y_plus
    d_log = 1.0_rp / log_arg

    ! Derivative of exponential correction:
    !   d/dy[7.8*(1 - exp(-y/11) - (y/11)*exp(-y/3))]
    ! = 7.8*[ (1/11)*exp(-y/11)  -  (1/11)*exp(-y/3)  +  (y/33)*exp(-y/3) ]
    !
    ! Grouping the last two terms:
    !   exp_term2 = (1/11)*exp(-y/3) - (y/33)*exp(-y/3)
    ! so d_exp = 7.8*(exp_term1 - exp_term2)  [note the MINUS sign]
    exp_term1 = exp(-y_plus / A_DAMP) / A_DAMP
    exp_term2 = (1.0_rp / A_DAMP) * exp(-y_plus / B_EXP) - &
                (y_plus / A_DAMP) * (1.0_rp / B_EXP) * exp(-y_plus / B_EXP)
    d_exp = EXP_COEFF * (exp_term1 - exp_term2)

    du_dy = d_log + d_exp

  end function reichardt_du_plus_dy

  !> Newton-Raphson solver for friction velocity u_τ using the Reichardt law.
  !!
  !! Solves the implicit equation:
  !!   F(u_τ) = u_τ * u⁺(y⁺) - U = 0
  !! where y⁺ = y * u_τ / ν and U is the tangential velocity magnitude.
  !!
  !! The Jacobian (by chain rule through y⁺):
  !!   dF/du_τ = u⁺ + u_τ * (du⁺/dy⁺) * (y/ν)
  !!
  !! @param u Tangential velocity magnitude U.
  !! @param y Wall-normal distance to the sampling point.
  !! @param guess Initial guess for u_τ.
  !! @param nu Kinematic viscosity.
  !! @param kappa Von Kármán constant.
  !! @return Friction velocity u_τ.
  function solve_reichardt_cpu(u, y, guess, nu, kappa) result(utau)
    real(kind=rp), intent(in) :: u, y, guess, nu, kappa
    real(kind=rp) :: utau
    real(kind=rp) :: y_plus, u_plus, du_dy
    real(kind=rp) :: error, f, df, utau_old
    integer :: k

    utau = guess

    do k = 1, 100
       utau_old = utau

       ! Dimensionless coordinates at current iterate
       y_plus = y * utau / nu
       u_plus = reichardt_u_plus(y_plus, kappa)
       du_dy  = reichardt_du_plus_dy(y_plus, kappa)

       ! Residual F and Jacobian dF/du_τ
       f  = utau * u_plus - u
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

  end function solve_reichardt_cpu

end module reichardt_cpu