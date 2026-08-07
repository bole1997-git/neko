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
module duprat_cpu
  use num_types, only: rp
  use logger, only: neko_log, NEKO_LOG_DEBUG, LOG_SIZE
  implicit none
  private

  public :: duprat_compute_cpu, duprat_compute_apg_cpu

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

  integer, parameter :: N_PANEL = 10

contains

  subroutine duprat_compute_cpu(u, v, w, ind_r, ind_s, ind_t, ind_e, &
       n_x, n_y, n_z, nu, rho_w, h, tau_x, tau_y, tau_z, &
       n_nodes, lx, nelv, kappa, beta, A, dpdx_const, tstep)
    integer,       intent(in) :: n_nodes, lx, nelv, tstep
    real(kind=rp), dimension(lx,lx,lx,nelv), intent(in) :: u, v, w
    integer,       intent(in),  dimension(n_nodes) :: ind_r, ind_s, ind_t, ind_e
    real(kind=rp), dimension(n_nodes), intent(in)    :: n_x, n_y, n_z
    real(kind=rp), dimension(n_nodes), intent(in)    :: nu, rho_w, h
    real(kind=rp), dimension(n_nodes), intent(inout) :: tau_x, tau_y, tau_z
    real(kind=rp), intent(in) :: kappa, beta, A, dpdx_const

    integer       :: i
    real(kind=rp) :: ui, vi, wi, normu, magu, utau, guess, tau_mag

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

       if (tstep .eq. 1) then
          guess = sqrt(magu * nu(i) / h(i))
       else
          tau_mag = sqrt(tau_x(i)**2 + tau_y(i)**2 + tau_z(i)**2)
          guess   = sqrt(tau_mag / max(rho_w(i), 1.0e-14_rp))
          guess   = max(guess, 1.0e-10_rp)
       end if

       utau = solve_duprat(magu, h(i), guess, nu(i), rho_w(i), dpdx_const, &
                           kappa, beta, A)

       tau_x(i) = -rho_w(i) * utau**2 * ui / magu
       tau_y(i) = -rho_w(i) * utau**2 * vi / magu
       tau_z(i) = -rho_w(i) * utau**2 * wi / magu
    end do

  end subroutine duprat_compute_cpu

  subroutine duprat_compute_apg_cpu(u, v, w, ind_r, ind_s, ind_t, ind_e, &
       n_x, n_y, n_z, nu, rho_w, h, tau_x, tau_y, tau_z, &
       n_nodes, lx, nelv, kappa, beta, A, dpdx_filt, tstep)
    integer,       intent(in) :: n_nodes, lx, nelv, tstep
    real(kind=rp), dimension(lx,lx,lx,nelv), intent(in) :: u, v, w
    integer,       intent(in),  dimension(n_nodes) :: ind_r, ind_s, ind_t, ind_e
    real(kind=rp), dimension(n_nodes), intent(in)    :: n_x, n_y, n_z
    real(kind=rp), dimension(n_nodes), intent(in)    :: nu, rho_w, h
    real(kind=rp), dimension(n_nodes), intent(inout) :: tau_x, tau_y, tau_z
    real(kind=rp), dimension(n_nodes), intent(in)    :: dpdx_filt
    real(kind=rp), intent(in) :: kappa, beta, A

    integer       :: i
    real(kind=rp) :: ui, vi, wi, normu, magu, utau, guess, tau_mag

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

       if (tstep .eq. 1) then
          guess = sqrt(magu * nu(i) / h(i))
       else
          tau_mag = sqrt(tau_x(i)**2 + tau_y(i)**2 + tau_z(i)**2)
          guess   = sqrt(tau_mag / max(rho_w(i), 1.0e-14_rp))
          guess   = max(guess, 1.0e-10_rp)
       end if

       utau = solve_duprat(magu, h(i), guess, nu(i), rho_w(i), dpdx_filt(i), &
                           kappa, beta, A)

       tau_x(i) = -rho_w(i) * utau**2 * ui / magu
       tau_y(i) = -rho_w(i) * utau**2 * vi / magu
       tau_z(i) = -rho_w(i) * utau**2 * wi / magu
    end do

  end subroutine duprat_compute_apg_cpu

  pure function nu_t_star(y_star, kappa, beta, A, alpha) result(nut)
    real(kind=rp), intent(in) :: y_star, kappa, beta, A, alpha
    real(kind=rp) :: nut, bracket, exp_damp

    if (y_star < 1.0e-12_rp) then
       nut = 0.0_rp
       return
    end if

    bracket  = alpha + y_star * (1.0_rp - alpha)**1.5_rp
    exp_damp = exp(-y_star / (1.0_rp + A*alpha**3))
    nut      = kappa * y_star * bracket**beta * (1.0_rp - exp_damp)**2
  end function nu_t_star

  pure function integrate_ode(y_star_max, kappa, beta, A, alpha, &
       sign_dpdx) result(U_star)
    real(kind=rp), intent(in) :: y_star_max, kappa, beta, A, alpha, sign_dpdx
    real(kind=rp) :: U_star, factor15, s_max, ds, s0, s, y, weight
    integer :: panel, j

    U_star = 0.0_rp
    if (y_star_max < 1.0e-12_rp) return

    factor15 = (1.0_rp - alpha)**1.5_rp
    s_max    = log(1.0_rp + y_star_max)
    ds       = s_max / real(N_PANEL, rp)

    do panel = 1, N_PANEL
       s0 = real(panel - 1, rp) * ds
       do j = 1, N_GL
          s      = s0 + GL_XI(j) * ds
          y      = exp(s) - 1.0_rp
          weight = GL_W(j) * ds * exp(s)
          U_star = U_star + weight * &
               (alpha + sign_dpdx * factor15 * y) / &
               (1.0_rp + nu_t_star(y, kappa, beta, A, alpha))
       end do
    end do

  end function integrate_ode

  function residual(utau, U_tang, y, nu, rho_w, dpdx, kappa, beta, A) &
       result(f)
    real(kind=rp), intent(in) :: utau, U_tang, y, nu, rho_w, dpdx
    real(kind=rp), intent(in) :: kappa, beta, A
    real(kind=rp) :: f, u_P, u_p_star, alpha, y_star, U_star, sign_dpdx

    if (abs(dpdx) < 1.0e-14_rp) then
       u_P = 0.0_rp
    else
       u_P = (nu * abs(dpdx) / max(rho_w, 1.0e-14_rp))**(1.0_rp/3.0_rp)
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

       f0    = residual(utau, U_tang, y, nu, rho_w, dpdx, kappa, beta, A)
       delta = max(1.0e-6_rp * utau, 1.0e-10_rp)
       f1    = residual(utau + delta, U_tang, y, nu, rho_w, dpdx, kappa, &
                        beta, A)
       df    = (f1 - f0) / delta

       if (abs(df) < 1.0e-14_rp) then
          utau = utau * 0.99_rp
          cycle
       end if

       utau = utau - f0 / df
       if (utau <= 0.0_rp) utau = utau_old * 0.5_rp

       error = abs((utau - utau_old) / (abs(utau) + 1.0e-16_rp))

       if (error < 1.0e-8_rp) then
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