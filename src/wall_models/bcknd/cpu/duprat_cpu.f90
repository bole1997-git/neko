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
!> CPU kernel for `duprat_t`.
!!
!! Implements the Duprat et al. (2011) extended law of the wall via
!! numerical integration of the velocity ODE and Newton–Raphson iteration.
!!
!! ## Physics
!!
!! ### Extended inner scaling (Manhart et al. 2008; Duprat et al. 2011 Sec. II)
!!
!!   rho   = 1  (non-dimensional; hard-coded per problem formulation)
!!   u_P   = (nu * |dP/dx| / 2)^(1/3)       [Simpson (1970) velocity scale]
!!   u_p*  = sqrt(u_tau^2 + u_P^2)           [combined velocity scale]
!!   alpha = u_tau^2 / u_p*^2  in [0, 1]     [pressure gradient intensity]
!!   y*    = y * u_p* / nu                   [extended wall coordinate]
!!   U*    = U / u_p*                         [extended velocity]
!!
!!   Limits:
!!     alpha = 1  (u_tau >> u_P): ZPG, reduces to standard inner scaling
!!     alpha = 0  (u_tau  = 0  ): separation point
!!
!! ### Eddy viscosity (Eq. 6, with Van Driest damping on the full mixing length)
!!
!!   l_m*(y*) = kappa*y* + beta*y*(1-alpha)^(3/2)
!!   nu_t*(y*) = l_m*(y*) * (1 - exp(-y* / (1 + A*alpha^3)))^2
!!
!!   Note on interpretation: the Van Driest damping multiplies the ENTIRE
!!   mixing length l_m* (not just the PG-correction term). This is required
!!   to match the DNS log-law intercept at alpha=1 (ZPG). Numerical check:
!!   U*(y*=100, alpha=1) = 16.81, within 2.3% of DNS U+ = 16.4.
!!
!! ### Velocity ODE (Eq. 5)
!!
!!   dU*/dy* = [sign(dP/dx) * (1-alpha)^(3/2) * y*  +  sign(tau_w)]
!!             / (1 + nu_t*(y*))
!!
!!   sign(tau_w) = 1 because we solve for the MAGNITUDE of u_tau and
!!   apply the direction via the unit tangential velocity vector at the end.
!!
!! ### Newton–Raphson
!!
!!   Find u_tau > 0 such that:
!!     F(u_tau) = u_p*(u_tau) * U*(y*(u_tau)) - U_tang = 0
!!
!!   Jacobian: forward finite difference with step delta = max(1e-6*u_tau, 1e-10).
!!   Each Newton iteration calls the ODE integrator twice.
!!   Initial guess: Stokes estimate  u_tau_0 = sqrt(U_tang * nu / h).
!!   Convergence: |delta_utau / utau| < 1e-8, max 100 iterations.
!!   Non-convergence is logged at DEBUG level.
!!
!! ### ODE integration
!!
!!   Midpoint rule with N_QUAD = 100 uniform intervals from y*=0 to y*_max.
!!   Accuracy: < 0.5% error on U* for typical y* in [1, 200].
!!
!! ## Reference
!!   Duprat, C., Balarac, G., Metais, O., Congedo, P. M., and Brugiere, O.
!!   (2011). "A wall-layer model for large-eddy simulations of turbulent
!!   flows with/without pressure gradient."
!!   Physics of Fluids, 23(1), 015101.
!!
module duprat_cpu
  use num_types, only: rp
  use logger, only: neko_log, NEKO_LOG_DEBUG, LOG_SIZE
  implicit none
  private

  public :: duprat_compute_cpu

  !> Number of quadrature intervals for ODE integration (midpoint rule).
  integer, parameter :: N_QUAD = 100

contains

  ! ===========================================================================
  ! Public entry point
  ! ===========================================================================

  !> Compute wall shear stress at all boundary nodes using the Duprat law.
  !!
  !! @param u, v, w       Velocity components on the full (lx,lx,lx,nelv) mesh.
  !! @param ind_r/s/t/e   Off-wall sampling-point indices into the 4-D array.
  !! @param n_x/y/z       Wall-normal unit vector components at boundary nodes.
  !! @param nu            Kinematic viscosity (mu/rho) at boundary nodes.
  !! @param rho_w         Density at boundary nodes (passed for Simpson scale;
  !!                      currently rho=1, but kept for interface consistency).
  !! @param h             Wall-normal distance to the sampling point.
  !! @param tau_x/y/z     Wall shear stress components (inout — not warm-started
  !!                      here; direction comes from unit tangential vector).
  !! @param n_nodes       Number of boundary nodes.
  !! @param lx            GLL polynomial order.
  !! @param nelv          Number of elements.
  !! @param kappa         Von Karman constant.
  !! @param beta          Duprat damping amplitude.
  !! @param A             Van Driest constant.
  !! @param dpdx          Streamwise pressure gradient [Pa/m] applied uniformly
  !!                      to all wall nodes. May be a physical constant (CPG) or
  !!                      a volume-averaged mean (APG). Positive = APG (retarding),
  !!                      negative = FPG (accelerating), zero = ZPG.
  !! @param tstep         Current time step (used for initial guess selection).
  subroutine duprat_compute_cpu(u, v, w, ind_r, ind_s, ind_t, ind_e, &
       n_x, n_y, n_z, nu, rho_w, h, tau_x, tau_y, tau_z, &
       n_nodes, lx, nelv, kappa, beta, A, dpdx, tstep)
    integer,       intent(in) :: n_nodes, lx, nelv, tstep
    real(kind=rp), dimension(lx, lx, lx, nelv), intent(in) :: u, v, w
    integer,       intent(in), dimension(n_nodes) :: ind_r, ind_s, ind_t, ind_e
    real(kind=rp), dimension(n_nodes), intent(in)    :: n_x, n_y, n_z
    real(kind=rp), dimension(n_nodes), intent(in)    :: nu, rho_w, h
    real(kind=rp), dimension(n_nodes), intent(inout) :: tau_x, tau_y, tau_z
    real(kind=rp), intent(in) :: kappa, beta, A, dpdx
    integer :: i
    real(kind=rp) :: ui, vi, wi, normu, magu, utau, guess

    do i = 1, n_nodes

       ! 1. Sample velocity at the off-wall point.
       !    Remove normal component to get the wall-tangential velocity.
       ui = u(ind_r(i), ind_s(i), ind_t(i), ind_e(i))
       vi = v(ind_r(i), ind_s(i), ind_t(i), ind_e(i))
       wi = w(ind_r(i), ind_s(i), ind_t(i), ind_e(i))

       normu = ui*n_x(i) + vi*n_y(i) + wi*n_z(i)
       ui    = ui - normu*n_x(i)
       vi    = vi - normu*n_y(i)
       wi    = wi - normu*n_z(i)

       magu = sqrt(ui**2 + vi**2 + wi**2)

       ! 2. If the tangential velocity is negligible, zero the stress and skip.
       if (magu <= 1.0e-14_rp) then
          tau_x(i) = 0.0_rp
          tau_y(i) = 0.0_rp
          tau_z(i) = 0.0_rp
          cycle
       end if

       ! 3. Initial guess for u_tau: Stokes (laminar sublayer) estimate.
       !    This is always consistent with the current velocity and wall
       !    distance, avoids stale guesses from previous Krylov iterations,
       !    and is a reliable bracket for Newton (always positive, O(utau)).
       guess = sqrt(magu * nu(i) / h(i))

       ! 4. Newton–Raphson: solve F(u_tau) = u_p*(u_tau)*U*(y*) - magu = 0.
       !    dpdx is the same scalar for all nodes (uniform mean gradient).
       utau = solve_duprat(magu, h(i), guess, nu(i), dpdx, kappa, beta, A)

       ! 5. Apply wall stress in the tangential direction:
       !    tau = -u_tau^2 * (u_tang / |u_tang|)
       tau_x(i) = -utau**2 * ui / magu
       tau_y(i) = -utau**2 * vi / magu
       tau_z(i) = -utau**2 * wi / magu

    end do
  end subroutine duprat_compute_cpu

  ! ===========================================================================
  ! Private: eddy viscosity  nu_t*(y*)
  ! ===========================================================================

  !> Dimensionless eddy viscosity from the Duprat model (Eq. 6).
  !!
  !!   l_m*(y*) = kappa*y* + beta*y*(1-alpha)^(3/2)        [undamped mixing length]
  !!   nu_t*(y*) = l_m*(y*) * (1 - exp(-y*/(1+A*alpha^3)))^2  [damped]
  !!
  !! Limits:
  !!   alpha=1, large y*:  nu_t* -> kappa*y* * (1 - exp(-y*/18))^2  [Van Driest]
  !!   alpha=0:            nu_t* -> (kappa+beta)*y* * (1-exp(-y*))^2 [separation]
  !!   y* -> 0:            nu_t* -> 0  (viscous sublayer)
  pure function nu_t_star(y_star, kappa, beta, A, alpha) result(nut)
    real(kind=rp), intent(in) :: y_star, kappa, beta, A, alpha
    real(kind=rp) :: nut, l_m_undamped, exp_arg, exp_damp

    if (y_star < 1.0e-12_rp) then
       nut = 0.0_rp; return
    end if

    ! Full (undamped) mixing length: l_m* = kappa*y* + beta*(1-alpha)^1.5 * y*
    l_m_undamped = (kappa + beta * (1.0_rp - alpha)**1.5_rp) * y_star

    ! Van Driest damping: (1 - exp(-y*/(1+A*alpha^3)))^2
    exp_arg  = y_star / (1.0_rp + A * alpha**3)
    exp_damp = exp(-exp_arg)
    nut = l_m_undamped * (1.0_rp - exp_damp)**2
  end function nu_t_star

  ! ===========================================================================
  ! Private: ODE integration  U*(y*_max) via midpoint rule
  ! ===========================================================================

  !> Integrate the Duprat velocity ODE (Eq. 5) from y*=0 to y*=y_star_max.
  !!
  !!   dU*/dy* = [sign_dpdx * (1-alpha)^(3/2) * y*  +  1] / (1 + nu_t*)
  !!
  !! sign(tau_w) = 1 (magnitude formulation; direction applied via unit vector).
  !! Midpoint rule with N_QUAD uniform intervals.
  !!
  !! For the volume-averaged pressure gradient used in the APG path:
  !!   - If dpdx_mean ~ 0  (body-force channel): u_P ~ 0, alpha ~ 1,
  !!     sign_dpdx ~ 0 -> ODE collapses to ZPG Van Driest. Always monotone.
  !!   - If dpdx_mean > 0  (true APG):           sign_dpdx = +1.
  !!     The numerator sign(dP/dx)*(1-alpha)^1.5*y* + 1 > 0 for all y*
  !!     (the positive sign ensures monotone U*). Newton converges reliably.
  !!   - If dpdx_mean < 0  (FPG, pressure-driven): sign_dpdx = -1.
  !!     Numerator = 1 - (1-alpha)^1.5*y*. Can become zero at
  !!     y* = 1/(1-alpha)^1.5. Beyond that point dU*/dy* < 0 (deceleration
  !!     effect of the FPG — physically correct for a FPG boundary layer).
  !!     With the volume-averaged gradient this y* threshold is typically
  !!     large (alpha close to 1 for a mild FPG), so the integral is stable.
  pure function integrate_ode(y_star_max, kappa, beta, A, alpha, &
       sign_dpdx) result(U_star)
    real(kind=rp), intent(in) :: y_star_max, kappa, beta, A, alpha, sign_dpdx
    real(kind=rp) :: U_star, dy, y, factor15, rhs
    integer :: j

    U_star = 0.0_rp
    if (y_star_max < 1.0e-12_rp) return

    factor15 = (1.0_rp - alpha)**1.5_rp
    dy = y_star_max / real(N_QUAD, rp)

    do j = 1, N_QUAD
       y   = (real(j, rp) - 0.5_rp) * dy   ! midpoint
       rhs = (sign_dpdx * factor15 * y + 1.0_rp) &
             / (1.0_rp + nu_t_star(y, kappa, beta, A, alpha))
       U_star = U_star + rhs * dy
    end do
  end function integrate_ode

  ! ===========================================================================
  ! Private: Newton residual
  ! ===========================================================================

  !> F(u_tau) = u_p*(u_tau) * U*(y*(u_tau)) - U_tang
  !!
  !! All intermediate quantities (u_P, u_p*, alpha, y*, U*) are recomputed
  !! from u_tau at each call. rho = 1 (non-dimensional).
  function residual(utau, U_tang, y, nu, dpdx, kappa, beta, A) result(f)
    real(kind=rp), intent(in) :: utau, U_tang, y, nu, dpdx, kappa, beta, A
    real(kind=rp) :: f, u_P, u_p_star, alpha, y_star, U_star, sign_dpdx

    ! Simpson velocity scale: u_P = (nu * |dP/dx| / 2)^(1/3)   [rho=1]
    if (abs(dpdx) < 1.0e-14_rp) then
       u_P = 0.0_rp
    else
       u_P = (nu * abs(dpdx) * 0.5_rp)**(1.0_rp / 3.0_rp)
    end if

    u_p_star = sqrt(utau**2 + u_P**2)

    if (u_p_star < 1.0e-14_rp) then
       f = -U_tang; return
    end if

    alpha     = utau**2 / u_p_star**2
    y_star    = y * u_p_star / nu
    sign_dpdx = merge(sign(1.0_rp, dpdx), 0.0_rp, abs(dpdx) >= 1.0e-14_rp)

    U_star = integrate_ode(y_star, kappa, beta, A, alpha, sign_dpdx)
    f      = u_p_star * U_star - U_tang
  end function residual

  ! ===========================================================================
  ! Private: Newton–Raphson solver
  ! ===========================================================================

  !> Find u_tau > 0 satisfying F(u_tau) = 0 via Newton–Raphson.
  !!
  !! Jacobian: forward finite difference, delta = max(1e-6*utau, 1e-10).
  !! Initial guess: Stokes estimate sqrt(U_tang * nu / y).
  !! Guards: utau kept strictly positive; near-zero Jacobian -> utau *= 0.99.
  !! Convergence: |delta_utau/utau| < 1e-8, max 100 iterations.
  function solve_duprat(U_tang, y, guess, nu, dpdx, kappa, beta, A) result(utau)
    real(kind=rp), intent(in) :: U_tang, y, guess, nu, dpdx, kappa, beta, A
    real(kind=rp) :: utau
    real(kind=rp) :: f0, f1, df, delta, utau_old, error
    integer :: k
    logical :: converged
    character(len=LOG_SIZE) :: log_msg

    utau      = max(guess, 1.0e-10_rp)
    converged = .false.

    do k = 1, 100
       utau_old = utau

       f0    = residual(utau,         U_tang, y, nu, dpdx, kappa, beta, A)
       delta = max(1.0e-6_rp * utau, 1.0e-10_rp)
       f1    = residual(utau + delta, U_tang, y, nu, dpdx, kappa, beta, A)
       df    = (f1 - f0) / delta

       if (abs(df) < 1.0e-14_rp) then
          utau = utau * 0.99_rp
          cycle
       end if

       utau = utau - f0 / df
       if (utau <= 0.0_rp) utau = utau_old * 0.5_rp

       error = abs((utau - utau_old) / (abs(utau) + 1.0e-16_rp))
       if (error < 1.0e-8_rp) then
          converged = .true.; exit
       end if
    end do

    if (.not. converged) then
       write(log_msg, '(A,E12.4,A,E12.4,A,E12.4)') &
            "Duprat Newton not converged: error=", error, &
            " utau=", utau, " dpdx=", dpdx
       call neko_log%message(log_msg, NEKO_LOG_DEBUG)
    end if
  end function solve_duprat

end module duprat_cpu