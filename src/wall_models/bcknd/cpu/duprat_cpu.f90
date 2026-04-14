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
!> CPU kernels for `duprat_t`.
!!
!! Two public entry points:
!!   duprat_compute_cpu     — CPG mode (uniform scalar dP/dx).
!!   duprat_compute_apg_cpu — APG mode (per-node IIR-filtered dP/dx array).
!!
!! ## Physics (Duprat et al. 2011, Phys. Fluids 23, 015101)
!!
!! Extended inner scaling (Manhart et al. 2008):
!!   rho   = 1  (non-dimensional; hard-coded per problem formulation)
!!   u_P   = (nu * |dP/dx| / 2)^(1/3)       [Simpson (1970) velocity scale]
!!   u_p*  = sqrt(u_tau^2 + u_P^2)           [combined velocity scale]
!!   alpha = u_tau^2 / u_p*^2  in [0, 1]     [pressure gradient intensity]
!!   y*    = y * u_p*/nu                     [extended wall coordinate]
!!   U*    = U / u_p*                         [extended velocity]
!!
!! Eddy viscosity (Eq. 6, Van Driest damping on FULL mixing length):
!!   l_m*(y*) = kappa*y* + beta*(1-alpha)^(3/2) * y*
!!   nu_t*(y*) = l_m*(y*) * (1 - exp(-y*/(1+A*alpha^3)))^2
!!
!!   At alpha=1: nu_t* = kappa*y*(1-exp(-y*/18))^2  [standard Van Driest, A=17]
!!   At alpha=0: nu_t* = (kappa+beta)*y*(1-exp(-y*))^2  [separation limit]
!!
!! Velocity ODE (Eq. 5):
!!   dU*/dy* = [sign(dP/dx)*(1-alpha)^(3/2)*y* + 1] / (1 + nu_t*)
!!   sign(tau_w) = 1 (magnitude solve; direction applied via unit tangent vector)
!!
!! Newton-Raphson:
!!   F(u_tau) = u_p*(u_tau) * U*(y*(u_tau)) - U_tang = 0
!!   Jacobian: forward finite difference, delta = max(1e-6*utau, 1e-10).
!!   Initial guess: Stokes estimate sqrt(magu*nu/h).
!!   Convergence: |delta_utau/utau| < 1e-8, max 100 iterations.
!!
!! ODE integration: midpoint rule, N_QUAD=100 uniform intervals.
!!
!! ## Self-consistent Stokes clamp (APG kernel only)
!!
!!   max_dpdx = 2 * (magu*nu/h)^(3/2) / nu
!!
!!   Derivation: u_tau_stokes = sqrt(magu*nu/h); set u_P_max = u_tau_stokes;
!!   invert u_P = (nu*|dpdx|/2)^(1/3) -> max_dpdx above.
!!   When clamped: u_P <= u_tau_stokes -> alpha >= 0.5 (GUARANTEED).
!!   ODE numerator = sign*0.354*y* + 1 > 0 for all y* -> monotone -> Newton OK.
!!   This is a SECONDARY defence; the primary defence is t_filter_start in duprat.f90.
!!
module duprat_cpu
  use num_types, only: rp
  use logger, only: neko_log, NEKO_LOG_DEBUG, LOG_SIZE
  implicit none
  private

  public :: duprat_compute_cpu, duprat_compute_apg_cpu

  !> Quadrature intervals for the midpoint-rule ODE integration.
  integer, parameter :: N_QUAD = 100

contains

  ! ===========================================================================
  ! Public: CPG / ZPG entry point (uniform scalar dpdx_const)
  ! ===========================================================================

  subroutine duprat_compute_cpu(u, v, w, ind_r, ind_s, ind_t, ind_e, &
       n_x, n_y, n_z, nu, rho_w, h, tau_x, tau_y, tau_z, &
       n_nodes, lx, nelv, kappa, beta, A, dpdx_const)
    integer,       intent(in) :: n_nodes, lx, nelv
    real(kind=rp), dimension(lx,lx,lx,nelv), intent(in) :: u, v, w
    integer,       intent(in), dimension(n_nodes) :: ind_r, ind_s, ind_t, ind_e
    real(kind=rp), dimension(n_nodes), intent(in)    :: n_x, n_y, n_z
    real(kind=rp), dimension(n_nodes), intent(in)    :: nu, rho_w, h
    real(kind=rp), dimension(n_nodes), intent(inout) :: tau_x, tau_y, tau_z
    real(kind=rp), intent(in) :: kappa, beta, A, dpdx_const
    integer :: i
    real(kind=rp) :: ui, vi, wi, normu, magu, utau, guess

    do i = 1, n_nodes
       ui    = u(ind_r(i), ind_s(i), ind_t(i), ind_e(i))
       vi    = v(ind_r(i), ind_s(i), ind_t(i), ind_e(i))
       wi    = w(ind_r(i), ind_s(i), ind_t(i), ind_e(i))
       normu = ui*n_x(i) + vi*n_y(i) + wi*n_z(i)
       ui    = ui - normu*n_x(i)
       vi    = vi - normu*n_y(i)
       wi    = wi - normu*n_z(i)
       magu  = sqrt(ui**2 + vi**2 + wi**2)

       if (magu <= 1.0e-14_rp) then
          tau_x(i) = 0.0_rp; tau_y(i) = 0.0_rp; tau_z(i) = 0.0_rp
          cycle
       end if

       guess = sqrt(magu * nu(i) / h(i))
       utau  = solve_duprat(magu, h(i), guess, nu(i), rho_w(i), dpdx_const, &
                            kappa, beta, A)

       tau_x(i) = -utau**2 * ui / magu
       tau_y(i) = -utau**2 * vi / magu
       tau_z(i) = -utau**2 * wi / magu
    end do
  end subroutine duprat_compute_cpu

  ! ===========================================================================
  ! Public: APG entry point (per-node filtered dpdx_filt array)
  ! ===========================================================================

  !> @param dpdx_filt  Per-node IIR-filtered streamwise dP/dx [Pa/m], n_nodes.
  !!   Zero for t < t_filter_start (ZPG). Converges to local mean for t >= t_filter_start.
  !!   Computed in duprat_compute using the PREVIOUS step's gradient (double-buffer).
  subroutine duprat_compute_apg_cpu(u, v, w, ind_r, ind_s, ind_t, ind_e, &
       n_x, n_y, n_z, nu, rho_w, h, tau_x, tau_y, tau_z, &
       n_nodes, lx, nelv, kappa, beta, A, dpdx_filt)
    integer,       intent(in) :: n_nodes, lx, nelv
    real(kind=rp), dimension(lx,lx,lx,nelv), intent(in) :: u, v, w
    integer,       intent(in), dimension(n_nodes) :: ind_r, ind_s, ind_t, ind_e
    real(kind=rp), dimension(n_nodes), intent(in)    :: n_x, n_y, n_z
    real(kind=rp), dimension(n_nodes), intent(in)    :: nu, rho_w, h
    real(kind=rp), dimension(n_nodes), intent(inout) :: tau_x, tau_y, tau_z
    real(kind=rp), dimension(n_nodes), intent(in)    :: dpdx_filt
    real(kind=rp), intent(in) :: kappa, beta, A
    integer :: i
    real(kind=rp) :: ui, vi, wi, normu, magu, utau, guess, dpdx_i, max_dpdx

    do i = 1, n_nodes
       ui    = u(ind_r(i), ind_s(i), ind_t(i), ind_e(i))
       vi    = v(ind_r(i), ind_s(i), ind_t(i), ind_e(i))
       wi    = w(ind_r(i), ind_s(i), ind_t(i), ind_e(i))
       normu = ui*n_x(i) + vi*n_y(i) + wi*n_z(i)
       ui    = ui - normu*n_x(i)
       vi    = vi - normu*n_y(i)
       wi    = wi - normu*n_z(i)
       magu  = sqrt(ui**2 + vi**2 + wi**2)

       if (magu <= 1.0e-14_rp) then
          tau_x(i) = 0.0_rp; tau_y(i) = 0.0_rp; tau_z(i) = 0.0_rp
          cycle
       end if

       ! Use the pre-filtered gradient (already clamped in duprat_compute).
       ! Apply the Stokes clamp here too as a secondary safety net.
       dpdx_i  = dpdx_filt(i)
       max_dpdx = 2.0_rp * (magu * nu(i) / h(i))**1.5_rp / nu(i)
       dpdx_i   = sign(min(abs(dpdx_i), max_dpdx), dpdx_i)

       guess = sqrt(magu * nu(i) / h(i))
       utau  = solve_duprat(magu, h(i), guess, nu(i), rho_w(i), dpdx_i, &
                            kappa, beta, A)

       tau_x(i) = -utau**2 * ui / magu
       tau_y(i) = -utau**2 * vi / magu
       tau_z(i) = -utau**2 * wi / magu
    end do
  end subroutine duprat_compute_apg_cpu

  ! ===========================================================================
  ! Private: eddy viscosity  nu_t*(y*)
  ! ===========================================================================

  pure function nu_t_star(y_star, kappa, beta, A, alpha) result(nut)
    real(kind=rp), intent(in) :: y_star, kappa, beta, A, alpha
    real(kind=rp) :: nut, l_m, exp_damp

    if (y_star < 1.0e-12_rp) then
       nut = 0.0_rp; return
    end if
    l_m      = (kappa + beta*(1.0_rp - alpha)**1.5_rp) * y_star
    exp_damp = exp(-y_star / (1.0_rp + A*alpha**3))
    nut      = l_m * (1.0_rp - exp_damp)**2
  end function nu_t_star

  ! ===========================================================================
  ! Private: ODE integration
  ! ===========================================================================

  pure function integrate_ode(y_star_max, kappa, beta, A, alpha, &
       sign_dpdx) result(U_star)
    real(kind=rp), intent(in) :: y_star_max, kappa, beta, A, alpha, sign_dpdx
    real(kind=rp) :: U_star, dy, y, factor15
    integer :: j

    U_star = 0.0_rp
    if (y_star_max < 1.0e-12_rp) return
    factor15 = (1.0_rp - alpha)**1.5_rp
    dy = y_star_max / real(N_QUAD, rp)
    do j = 1, N_QUAD
       y      = (real(j, rp) - 0.5_rp) * dy
       U_star = U_star + &
            (sign_dpdx * factor15 * y + 1.0_rp) / &
            (1.0_rp + nu_t_star(y, kappa, beta, A, alpha)) * dy
    end do
  end function integrate_ode

  ! ===========================================================================
  ! Private: Newton residual
  ! ===========================================================================

  function residual(utau, U_tang, y, nu, rho_w, dpdx, kappa, beta, A) result(f)
    real(kind=rp), intent(in) :: utau, U_tang, y, nu, rho_w, dpdx
    real(kind=rp), intent(in) :: kappa, beta, A
    real(kind=rp) :: f, u_P, u_p_star, alpha, y_star, U_star, sign_dpdx

    if (abs(dpdx) < 1.0e-14_rp) then
       u_P = 0.0_rp
    else
       u_P = (nu * abs(dpdx) / (2.0_rp * rho_w))**(1.0_rp/3.0_rp)
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
  ! Private: Newton-Raphson solver
  ! ===========================================================================

  function solve_duprat(U_tang, y, guess, nu, rho_w, dpdx, kappa, beta, A) &
       result(utau)
    real(kind=rp), intent(in) :: U_tang, y, guess, nu, rho_w, dpdx
    real(kind=rp), intent(in) :: kappa, beta, A
    real(kind=rp) :: utau, f0, f1, df, delta, utau_old, error
    integer :: k
    logical :: converged
    character(len=LOG_SIZE) :: log_msg

    utau      = max(guess, 1.0e-10_rp)
    converged = .false.

    do k = 1, 100
       utau_old = utau
       f0    = residual(utau,         U_tang, y, nu, rho_w, dpdx, kappa, beta, A)
       delta = max(1.0e-6_rp * utau, 1.0e-10_rp)
       f1    = residual(utau + delta, U_tang, y, nu, rho_w, dpdx, kappa, beta, A)
       df    = (f1 - f0) / delta

       if (abs(df) < 1.0e-14_rp) then
          utau = utau * 0.99_rp; cycle
       end if

       utau = utau - f0 / df
       if (utau <= 0.0_rp) utau = utau_old * 0.5_rp

       error = abs((utau - utau_old) / (abs(utau) + 1.0e-16_rp))
       if (error < 1.0e-8_rp) then
          converged = .true.; exit
       end if
    end do

    if (.not. converged) then
       write(log_msg, '(A,E10.3,A,E10.3,A,E10.3)') &
            "Duprat NC: err=", error, " utau=", utau, " dp=", dpdx
       call neko_log%message(log_msg, NEKO_LOG_DEBUG)
    end if
  end function solve_duprat

end module duprat_cpu