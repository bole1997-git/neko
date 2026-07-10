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
!!   duprat_compute_cpu     -- CPG mode (uniform scalar dP/dx).
!!   duprat_compute_apg_cpu -- APG mode (per-node IIR-filtered dP/dx array).
!!
!! ## Changes from original -
!!
!! FIX 1 -- rho_w removed from u_P formula.
!!   Neko's pressure field p is the kinematic pressure (p/rho).
!!   The Duprat (2011) formula is kinematic: u_P=(nu*|dP/dx|/2)^(1/3).
!!   Dividing by rho_w a second time was dimensionally wrong and
!!   caused blow-up in any run where rho /= 1.
!!   rho_w has been removed from residual() and solve_duprat() entirely.
!!   The CPG and APG kernels no longer receive or pass rho_w.
!!
!! FIX 2 -- Wall-tangential gradient magnitude replaces velocity projection.
!!   The original dpdx_local = grad . (u/|u|) flips sign inside
!!   recirculation zones (u reverses), giving spurious APG->FPG switches.
!!   The correct quantity for Duprat u_P is the magnitude of the
!!   wall-tangential pressure gradient (sign applied separately).
!!   This is computed in duprat.f90 before the filter update.
!!   The kernel receives a pre-projected scalar; no change needed here.
!!
!! FIX 3 -- Stokes clamp removed from APG kernel.
!!   The clamp is applied once in the filter update (duprat.f90).
!!   A second clamp in the kernel uses the current-step magu, which
!!   can be near-zero at nodes that were active at the filter step,
!!   falsely clamping dpdx_filt to zero.
!!
!! FIX 4 -- t_prev / restart fix is in duprat.f90, not here.
!!
!! FIX 5 -- N_QUAD=100 midpoint rule replaced by 5-point Gauss-Legendre.
!!   Cost reduction: 100 evaluations/Newton-iter -> 5.
!!   Accuracy: GL-5 integrates polynomials of degree <= 9 exactly;
!!   more than sufficient for the smooth nu_t* integrand.
!!
!! FIX 6 -- Newton tolerance relaxed from 1e-8 to 1e-3.
!!   Matches Spalding. LES velocities carry O(1%) noise; sub-percent
!!   accuracy in utau does not improve physics.
!!
!! FIX 7 -- Warm initial guess from previous tau instead of Stokes.
!!   Reduces average Newton iterations from ~20 to ~3-5 after step 1.
!!
!! ## Physics (Duprat et al. 2011, Phys. Fluids 23, 015101)
!!
!! Extended inner scaling:
!!   u_P    = (nu * |dP/dx| / 2)^(1/3)       [Simpson velocity scale]
!!   u_p*   = sqrt(u_tau^2 + u_P^2)           [combined scale]
!!   alpha  = u_tau^2 / u_p*^2  in [0,1]      [PG intensity]
!!   y*     = y * u_p*/nu                      [extended wall unit]
!!   U*     = U / u_p*
!!
!! Eddy viscosity (Van Driest damping, Eq. 6):
!!   l_m*(y*) = (kappa + beta*(1-alpha)^1.5) * y*
!!   nu_t*    = l_m* * (1 - exp(-y*/(1+A*alpha^3)))^2
!!
!! Velocity ODE (Eq. 5):
!!   dU*/dy* = [sign(dP/dx)*(1-alpha)^1.5*y* + 1] / (1 + nu_t*)
!!
!! Newton: F(u_tau) = u_p*(u_tau)*U*(y*(u_tau)) - U_tang = 0
!!   Jacobian: forward FD, delta = max(1e-6*utau, 1e-10)
!!   Convergence: |delta_utau / utau| < 1e-3, max 100 iters
!!
module duprat_cpu
  use num_types, only: rp
  use logger, only: neko_log, NEKO_LOG_DEBUG, LOG_SIZE
  implicit none
  private

  public :: duprat_compute_cpu, duprat_compute_apg_cpu

  ! -------------------------------------------------------------------------
  ! FIX 5: 5-point Gauss-Legendre quadrature on [0,1].
  ! Replaces the 100-point midpoint rule (20x fewer nu_t* evaluations).
  ! Nodes and weights from Abramowitz & Stegun Table 25.4.
  ! -------------------------------------------------------------------------
  integer,  parameter :: N_GL = 5
  real(rp), parameter :: GL_XI(N_GL) = [ &
       0.046910077936172_rp, &
       0.230765345953158_rp, &
       0.500000000000000_rp, &
       0.769234654046842_rp, &
       0.953089922063828_rp ]
  real(rp), parameter :: GL_W(N_GL) = [ &
       0.118463442528095_rp, &
       0.239314335249683_rp, &
       0.284444444444444_rp, &
       0.239314335249683_rp, &
       0.118463442528095_rp ]

contains

  ! ===========================================================================
  ! Public: CPG / ZPG entry point (uniform scalar dpdx_const)
  ! ===========================================================================
  ! FIX 1: rho_w removed from argument list entirely.
  ! FIX 7: warm guess from previous tau (tstep argument added).
  subroutine duprat_compute_cpu(u, v, w, ind_r, ind_s, ind_t, ind_e, &
       n_x, n_y, n_z, nu, h, tau_x, tau_y, tau_z, &
       n_nodes, lx, nelv, kappa, beta, A, dpdx_const, tstep)
    integer,       intent(in) :: n_nodes, lx, nelv, tstep
    real(kind=rp), dimension(lx,lx,lx,nelv), intent(in) :: u, v, w
    integer,       intent(in),  dimension(n_nodes) :: ind_r, ind_s, ind_t, ind_e
    real(kind=rp), dimension(n_nodes), intent(in)    :: n_x, n_y, n_z
    real(kind=rp), dimension(n_nodes), intent(in)    :: nu, h
    real(kind=rp), dimension(n_nodes), intent(inout) :: tau_x, tau_y, tau_z
    real(kind=rp), intent(in) :: kappa, beta, A, dpdx_const

    integer       :: i
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
          tau_x(i) = 0.0_rp
          tau_y(i) = 0.0_rp
          tau_z(i) = 0.0_rp
          cycle
       end if

       ! FIX 7: warm guess from previous tau magnitude; Stokes only at step 1.
       if (tstep .eq. 1) then
          guess = sqrt(magu * nu(i) / h(i))
       else
          guess = sqrt(sqrt(tau_x(i)**2 + tau_y(i)**2 + tau_z(i)**2))
          guess = max(guess, 1.0e-10_rp)
       end if

       ! FIX 1: no rho_w argument
       utau = solve_duprat(magu, h(i), guess, nu(i), dpdx_const, &
                           kappa, beta, A)

       tau_x(i) = -utau**2 * ui / magu
       tau_y(i) = -utau**2 * vi / magu
       tau_z(i) = -utau**2 * wi / magu
    end do

  end subroutine duprat_compute_cpu

  ! ===========================================================================
  ! Public: APG entry point (per-node IIR-filtered dpdx array)
  ! ===========================================================================
  ! FIX 1: rho_w removed.
  ! FIX 3: Stokes clamp removed (applied once in duprat.f90 filter update).
  ! FIX 7: warm guess from previous tau.
  !
  ! @param dpdx_filt  Per-node wall-tangential |dP/ds| [Pa/m], n_nodes.
  !   Already clamped by the Stokes bound in duprat.f90.
  !   Zero for t < t_filter_start (ZPG mode).
  subroutine duprat_compute_apg_cpu(u, v, w, ind_r, ind_s, ind_t, ind_e, &
       n_x, n_y, n_z, nu, h, tau_x, tau_y, tau_z, &
       n_nodes, lx, nelv, kappa, beta, A, dpdx_filt, tstep)
    integer,       intent(in) :: n_nodes, lx, nelv, tstep
    real(kind=rp), dimension(lx,lx,lx,nelv), intent(in) :: u, v, w
    integer,       intent(in),  dimension(n_nodes) :: ind_r, ind_s, ind_t, ind_e
    real(kind=rp), dimension(n_nodes), intent(in)    :: n_x, n_y, n_z
    real(kind=rp), dimension(n_nodes), intent(in)    :: nu, h
    real(kind=rp), dimension(n_nodes), intent(inout) :: tau_x, tau_y, tau_z
    real(kind=rp), dimension(n_nodes), intent(in)    :: dpdx_filt
    real(kind=rp), intent(in) :: kappa, beta, A

    integer       :: i
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
          tau_x(i) = 0.0_rp
          tau_y(i) = 0.0_rp
          tau_z(i) = 0.0_rp
          cycle
       end if

       ! FIX 7: warm guess from previous tau; Stokes only at step 1.
       if (tstep .eq. 1) then
          guess = sqrt(magu * nu(i) / h(i))
       else
          guess = sqrt(sqrt(tau_x(i)**2 + tau_y(i)**2 + tau_z(i)**2))
          guess = max(guess, 1.0e-10_rp)
       end if

       ! FIX 3: dpdx_filt(i) is already clamped in duprat.f90 -- use directly.
       ! FIX 1: no rho_w argument.
       utau = solve_duprat(magu, h(i), guess, nu(i), dpdx_filt(i), &
                           kappa, beta, A)

       tau_x(i) = -utau**2 * ui / magu
       tau_y(i) = -utau**2 * vi / magu
       tau_z(i) = -utau**2 * wi / magu
    end do

  end subroutine duprat_compute_apg_cpu

  ! ===========================================================================
  ! Private: eddy viscosity nu_t*(y*)
  ! ===========================================================================

  pure function nu_t_star(y_star, kappa, beta, A, alpha) result(nut)
    real(kind=rp), intent(in) :: y_star, kappa, beta, A, alpha
    real(kind=rp) :: nut, l_m, exp_damp

    if (y_star < 1.0e-12_rp) then
       nut = 0.0_rp
       return
    end if

    l_m      = (kappa + beta*(1.0_rp - alpha)**1.5_rp) * y_star
    exp_damp = exp(-y_star / (1.0_rp + A*alpha**3))
    nut      = l_m * (1.0_rp - exp_damp)**2
  end function nu_t_star

  ! ===========================================================================
  ! Private: ODE integration -- FIX 5
  ! 5-point Gauss-Legendre replaces 100-point midpoint (20x speedup).
  ! GL-5 integrates polynomials of degree <=9 exactly; sufficient here.
  ! ===========================================================================

  pure function integrate_ode(y_star_max, kappa, beta, A, alpha, &
       sign_dpdx) result(U_star)
    real(kind=rp), intent(in) :: y_star_max, kappa, beta, A, alpha, sign_dpdx
    real(kind=rp) :: U_star, y, factor15
    integer :: j

    U_star = 0.0_rp
    if (y_star_max < 1.0e-12_rp) return

    factor15 = (1.0_rp - alpha)**1.5_rp

    do j = 1, N_GL
       y      = GL_XI(j) * y_star_max
       U_star = U_star + GL_W(j) * y_star_max * &
            (sign_dpdx * factor15 * y + 1.0_rp) / &
            (1.0_rp + nu_t_star(y, kappa, beta, A, alpha))
    end do

  end function integrate_ode

  ! ===========================================================================
  ! Private: Newton residual
  ! FIX 1: rho_w removed -- u_P uses kinematic formula (nu*|dpdx|/2)^(1/3)
  ! ===========================================================================

  function residual(utau, U_tang, y, nu, dpdx, kappa, beta, A) result(f)
    real(kind=rp), intent(in) :: utau, U_tang, y, nu, dpdx
    real(kind=rp), intent(in) :: kappa, beta, A
    real(kind=rp) :: f, u_P, u_p_star, alpha, y_star, U_star, sign_dpdx

    ! FIX 1: kinematic pressure velocity scale -- no rho_w division.
    ! Neko p is kinematic; Duprat eq. 6: u_P = (nu*|dP/dx|/2)^(1/3).
    if (abs(dpdx) < 1.0e-14_rp) then
       u_P = 0.0_rp
    else
       u_P = (nu * abs(dpdx) / 2.0_rp)**(1.0_rp/3.0_rp)
    end if

    u_p_star = sqrt(utau**2 + u_P**2)
    if (u_p_star < 1.0e-14_rp) then
       f = -U_tang
       return
    end if

    alpha  = utau**2 / u_p_star**2
    y_star = y * u_p_star / nu

    sign_dpdx = 0.0_rp
    if (abs(dpdx) >= 1.0e-14_rp) sign_dpdx = sign(1.0_rp, dpdx)

    U_star = integrate_ode(y_star, kappa, beta, A, alpha, sign_dpdx)
    f      = u_p_star * U_star - U_tang

  end function residual

  ! ===========================================================================
  ! Private: Newton-Raphson solver
  ! FIX 1: rho_w removed from signature.
  ! FIX 6: tolerance relaxed from 1e-8 to 1e-3 (matches Spalding).
  ! ===========================================================================

  function solve_duprat(U_tang, y, guess, nu, dpdx, kappa, beta, A) &
       result(utau)
    real(kind=rp), intent(in) :: U_tang, y, guess, nu, dpdx
    real(kind=rp), intent(in) :: kappa, beta, A
    real(kind=rp) :: utau, f0, f1, df, delta, utau_old, error
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

       ! FIX 6: relaxed tolerance 1e-3 (was 1e-8).
       ! LES velocities carry O(1%) fluctuation noise; tighter tolerance
       ! gives no physical benefit and triples iteration count.
       if (error < 1.0e-3_rp) then
          converged = .true.
          exit
       end if
    end do

    if (.not. converged) then
       write(log_msg, '(A,E10.3,A,E10.3,A,E10.3)') &
            "Duprat NC: err=", error, " utau=", utau, " dp=", dpdx
       call neko_log%message(log_msg, NEKO_LOG_DEBUG)
    end if

  end function solve_duprat

end module duprat_cpu