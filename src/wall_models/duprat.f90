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
!> Implements `duprat_t`.
!!
!! Wall model based on Duprat et al. (2011) extended law of the wall.
!!
!! ## Physics
!!
!!   rho   = 1  (non-dimensional; hard-coded per problem setup)
!!   u_P   = (nu * |dP/dx| / 2)^(1/3)           [Simpson velocity scale]
!!   u_p*  = sqrt(u_tau^2 + u_P^2)
!!   alpha = u_tau^2 / u_p*^2  in [0, 1]
!!   y*    = y * u_p*/nu,   U* = U/u_p*
!!   nu_t* = (kappa*y* + beta*y*(1-alpha)^1.5)*(1-exp(-y*/(1+A*alpha^3)))^2
!!   dU*/dy* = [sign(dP/dx)*(1-alpha)^1.5*y* + 1] / (1 + nu_t*)
!!
!! ## Pressure gradient modes
!!
!! **CPG** (`use_constant_dpdx: true`, default):
!!   Uniform scalar `dpdx_constant` [Pa/m] at every wall node. Set to 0 for ZPG.
!!
!! **APG** (`use_constant_dpdx: false`):
!!   Per-node temporally filtered local streamwise pressure gradient.
!!
!! ## APG framework design
!!
!! Three structural problems caused all previous blow-ups, each addressed here:
!!
!! ### Problem A — Gradient timing
!!
!!   `compute()` is called from the Momentum Krylov solver. At step N's
!!   first call, the pnpn ordering (Pressure solve -> Velocity correct ->
!!   Momentum solve) means `p` in the registry is the CONVERGED step-N
!!   pressure. Using `dudxyz(p_N)` during step N gives the current-step's
!!   gradient, which can be large during turbulent transition even after
!!   the previous step was stable.
!!
!!   FIX: double-buffer the gradient.
!!     grad_px/py/pz     = gradient from p_{N-1} (PREVIOUS step, safe to use)
!!     grad_px_buf/...   = gradient from p_N (CURRENT step, stored for next step)
!!
!!   At step N's first call:
!!     1. Use grad_px (= p_{N-1} gradient) to update the IIR filter.
!!     2. Compute dudxyz(p_N) -> grad_px_buf.
!!     3. Copy buf -> active: grad_px <- grad_px_buf.
!!     4. Record grad_tstep = tstep.
!!
!!   This guarantees the filter always operates on a one-step-lagged,
!!   previously-converged pressure gradient.
!!
!! ### Problem B — Filter warm-up during transition
!!
!!   For t < t_filter_start, the flow is in turbulent transition. The
!!   pressure field has physically meaningless transient values. Accumulating
!!   the IIR filter during this period seeds it with garbage, which compounds
!!   over many steps even with the Stokes clamp active.
!!
!!   FIX: the `t_filter_start` parameter (JSON, default 5.0 outer time units).
!!     t <  t_filter_start: dpdx_filt(:) = 0 (ZPG), filter NOT updated.
!!     t >= t_filter_start: IIR filter accumulates normally.
!!
!!   The double-buffer IS updated every step regardless, so at t_filter_start
!!   the filter immediately has access to a correct, recent p_{N-1} gradient.
!!
!!   Default t_filter_start = 5.0 ensures transition is complete before the
!!   APG correction activates. The model runs as ZPG Duprat (Van Driest)
!!   during transition — physically correct and numerically stable.
!!
!! ### Problem C — Inexact dt
!!
!!   Previous code estimated dt = t/(tstep-1), giving the MEAN dt, not the
!!   current dt. For variable-timestep runs, this gives wrong beta_f.
!!
!!   FIX: store `t_prev` as a type member. Exact dt = t - t_prev.
!!
!! ### Clamp (secondary defence)
!!
!!   Stokes-based self-consistent clamp applied to dpdx_local before the
!!   filter update:
!!
!!     max_dpdx = 2 * (magu * nu / h)^(3/2) / nu
!!
!!   Derived by setting u_P_max = u_tau_stokes = sqrt(magu*nu/h) and
!!   inverting the Simpson scale. When clamped, alpha >= 0.5 is guaranteed,
!!   keeping the ODE numerator positive and Newton convergent.
!!
!! ## JSON parameters
!!
!!   kappa            : real    (default 0.41)
!!   beta             : real    (default 0.78)
!!   A                : real    (default 17.0)
!!   use_constant_dpdx: logical (default .true.)
!!   dpdx_constant    : real    (default 0.0  [Pa/m])
!!   t_filter_start   : real    (default 5.0  [outer time units])
!!   t_filter         : real    (default 1.0  [outer time units])
!!
!! ## Reference
!!   Duprat, C., Balarac, G., Metais, O., Congedo, P. M., and Brugiere, O.
!!   (2011). Physics of Fluids, 23(1), 015101.
!!
module duprat
  use field, only: field_t
  use num_types, only: rp
  use json_module, only: json_file
  use coefs, only: coef_t
  use neko_config, only: NEKO_BCKND_DEVICE
  use wall_model, only: wall_model_t
  use registry, only: neko_registry
  use json_utils, only: json_get_or_lookup, json_get_or_default
  use field_math, only: field_invcol3
  use vector, only: vector_t
  use math, only: masked_gather_copy_0
  use device_math, only: device_masked_gather_copy_0
  use scratch_registry, only: neko_scratch_registry
  use operators, only: dudxyz
  use duprat_cpu, only: duprat_compute_cpu, duprat_compute_apg_cpu
  use utils, only: neko_error

  implicit none
  private

  type, public, extends(wall_model_t) :: duprat_t
     ! Model constants
     real(kind=rp) :: kappa = 0.41_rp
     real(kind=rp) :: beta  = 0.78_rp
     real(kind=rp) :: A     = 17.0_rp
     ! Mode
     logical        :: use_constant_dpdx = .true.
     real(kind=rp) :: dpdx_const = 0.0_rp
     ! APG filter parameters
     real(kind=rp) :: t_filter_start = 5.0_rp  ! start filter at this time
     real(kind=rp) :: t_filter       = 1.0_rp  ! IIR time scale
     ! Nu and rho at wall nodes
     type(vector_t) :: nu
     type(vector_t) :: rho_w
     ! Per-node IIR filtered streamwise pressure gradient [Pa/m]
     ! Zero until t >= t_filter_start. Converges to local mean afterwards.
     real(kind=rp), allocatable :: dpdx_filt(:)
     ! ACTIVE gradient arrays: hold grad(p_{N-1}) — the PREVIOUS step's gradient.
     ! Safe to use for the filter update at step N.
     real(kind=rp), allocatable :: grad_px(:,:,:,:)
     real(kind=rp), allocatable :: grad_py(:,:,:,:)
     real(kind=rp), allocatable :: grad_pz(:,:,:,:)
     ! BUFFER gradient arrays: hold grad(p_N) — the CURRENT step's gradient.
     ! Computed at step N's first call. Promoted to active at step N+1.
     real(kind=rp), allocatable :: grad_px_buf(:,:,:,:)
     real(kind=rp), allocatable :: grad_py_buf(:,:,:,:)
     real(kind=rp), allocatable :: grad_pz_buf(:,:,:,:)
     ! Krylov guard: fires update on FIRST call per timestep only
     integer        :: grad_tstep = 0
     ! Previous time for exact dt computation
     real(kind=rp) :: t_prev = -1.0_rp   ! -1 signals "first call ever"
   contains
     procedure, pass(this) :: init              => duprat_init
     procedure, pass(this) :: partial_init      => duprat_partial_init
     procedure, pass(this) :: finalize          => duprat_finalize
     procedure, pass(this) :: init_from_components => duprat_init_from_components
     procedure, pass(this) :: free              => duprat_free
     procedure, pass(this) :: compute_nu        => duprat_compute_nu
     procedure, pass(this) :: compute           => duprat_compute
  end type duprat_t

contains

  ! ===========================================================================
  ! Constructors / destructor
  ! ===========================================================================

  subroutine duprat_init(this, scheme_name, coef, msk, facet, h_index, json)
    class(duprat_t), intent(inout) :: this
    character(len=*), intent(in)   :: scheme_name
    type(coef_t),     intent(in)   :: coef
    integer,          intent(in)   :: msk(:), facet(:), h_index
    type(json_file),  intent(inout):: json
    real(kind=rp) :: kappa, beta, A, dpdx_const, t_filter_start, t_filter
    logical       :: use_constant_dpdx

    call json_get_or_lookup(json, "kappa", kappa)
    call json_get_or_lookup(json, "beta",  beta)
    call json_get_or_lookup(json, "A",     A)
    call json_get_or_default(json, "use_constant_dpdx", use_constant_dpdx, .true.)
    call json_get_or_default(json, "dpdx_constant",     dpdx_const,      0.0_rp)
    call json_get_or_default(json, "t_filter_start",    t_filter_start,  5.0_rp)
    call json_get_or_default(json, "t_filter",          t_filter,        1.0_rp)

    call this%init_from_components(scheme_name, coef, msk, facet, h_index, &
         kappa, beta, A, use_constant_dpdx, dpdx_const, t_filter_start, t_filter)
  end subroutine duprat_init

  subroutine duprat_partial_init(this, coef, json)
    class(duprat_t), intent(inout) :: this
    type(coef_t),    intent(in)    :: coef
    type(json_file), intent(inout) :: json

    call this%partial_init_base(coef, json)
    call json_get_or_lookup(json, "kappa", this%kappa)
    call json_get_or_lookup(json, "beta",  this%beta)
    call json_get_or_lookup(json, "A",     this%A)
    call json_get_or_default(json, "use_constant_dpdx", this%use_constant_dpdx, .true.)
    call json_get_or_default(json, "dpdx_constant",     this%dpdx_const,      0.0_rp)
    call json_get_or_default(json, "t_filter_start",    this%t_filter_start,  5.0_rp)
    call json_get_or_default(json, "t_filter",          this%t_filter,        1.0_rp)
  end subroutine duprat_partial_init

  subroutine duprat_finalize(this, msk, facet)
    class(duprat_t), intent(inout) :: this
    integer,         intent(in)    :: msk(:), facet(:)
    integer :: lx, nelv

    call this%finalize_base(msk, facet)
    call this%nu%init(this%n_nodes)
    call this%rho_w%init(this%n_nodes)
    this%grad_tstep = 0
    this%t_prev     = -1.0_rp

    if (.not. this%use_constant_dpdx) then
       lx   = this%coef%Xh%lx
       nelv = this%coef%msh%nelv
       allocate(this%dpdx_filt   (this%n_nodes))
       allocate(this%grad_px     (lx, lx, lx, nelv))
       allocate(this%grad_py     (lx, lx, lx, nelv))
       allocate(this%grad_pz     (lx, lx, lx, nelv))
       allocate(this%grad_px_buf (lx, lx, lx, nelv))
       allocate(this%grad_py_buf (lx, lx, lx, nelv))
       allocate(this%grad_pz_buf (lx, lx, lx, nelv))
       this%dpdx_filt    = 0.0_rp
       this%grad_px      = 0.0_rp
       this%grad_py      = 0.0_rp
       this%grad_pz      = 0.0_rp
       this%grad_px_buf  = 0.0_rp
       this%grad_py_buf  = 0.0_rp
       this%grad_pz_buf  = 0.0_rp
    end if
  end subroutine duprat_finalize

  ! NOTE: init_from_components must allocate APG arrays too because
  ! init_base calls finalize_base only, NOT the child duprat_finalize.
  subroutine duprat_init_from_components(this, scheme_name, coef, msk, &
       facet, h_index, kappa, beta, A, use_constant_dpdx, dpdx_const, &
       t_filter_start, t_filter)
    class(duprat_t), intent(inout) :: this
    character(len=*), intent(in)   :: scheme_name
    type(coef_t),     intent(in)   :: coef
    integer,          intent(in)   :: msk(:), facet(:), h_index
    real(kind=rp),    intent(in)   :: kappa, beta, A, dpdx_const
    real(kind=rp),    intent(in)   :: t_filter_start, t_filter
    logical,          intent(in)   :: use_constant_dpdx
    integer :: lx, nelv

    call this%free()
    call this%init_base(scheme_name, coef, msk, facet, h_index)

    this%kappa             = kappa
    this%beta              = beta
    this%A                 = A
    this%use_constant_dpdx = use_constant_dpdx
    this%dpdx_const        = dpdx_const
    this%t_filter_start    = t_filter_start
    this%t_filter          = t_filter
    this%grad_tstep        = 0
    this%t_prev            = -1.0_rp

    call this%nu%init(this%n_nodes)
    call this%rho_w%init(this%n_nodes)

    if (.not. this%use_constant_dpdx) then
       lx   = this%coef%Xh%lx
       nelv = this%coef%msh%nelv
       allocate(this%dpdx_filt   (this%n_nodes))
       allocate(this%grad_px     (lx, lx, lx, nelv))
       allocate(this%grad_py     (lx, lx, lx, nelv))
       allocate(this%grad_pz     (lx, lx, lx, nelv))
       allocate(this%grad_px_buf (lx, lx, lx, nelv))
       allocate(this%grad_py_buf (lx, lx, lx, nelv))
       allocate(this%grad_pz_buf (lx, lx, lx, nelv))
       this%dpdx_filt    = 0.0_rp
       this%grad_px      = 0.0_rp
       this%grad_py      = 0.0_rp
       this%grad_pz      = 0.0_rp
       this%grad_px_buf  = 0.0_rp
       this%grad_py_buf  = 0.0_rp
       this%grad_pz_buf  = 0.0_rp
    end if
  end subroutine duprat_init_from_components

  subroutine duprat_free(this)
    class(duprat_t), intent(inout) :: this
    call this%free_base()
    call this%nu%free()
    call this%rho_w%free()
    if (allocated(this%dpdx_filt   )) deallocate(this%dpdx_filt   )
    if (allocated(this%grad_px     )) deallocate(this%grad_px     )
    if (allocated(this%grad_py     )) deallocate(this%grad_py     )
    if (allocated(this%grad_pz     )) deallocate(this%grad_pz     )
    if (allocated(this%grad_px_buf )) deallocate(this%grad_px_buf )
    if (allocated(this%grad_py_buf )) deallocate(this%grad_py_buf )
    if (allocated(this%grad_pz_buf )) deallocate(this%grad_pz_buf )
  end subroutine duprat_free

  ! ===========================================================================
  ! Nu / rho gather
  ! ===========================================================================

  subroutine duprat_compute_nu(this)
    class(duprat_t), intent(inout) :: this
    type(field_t), pointer :: temp
    integer :: idx

    call neko_scratch_registry%request_field(temp, idx, .false.)
    call field_invcol3(temp, this%mu, this%rho)

    if (NEKO_BCKND_DEVICE .eq. 1) then
       call device_masked_gather_copy_0(this%nu%x_d, temp%x_d, this%msk_d, &
            temp%size(), this%nu%size())
       call device_masked_gather_copy_0(this%rho_w%x_d, this%rho%x_d, &
            this%msk_d, this%rho%size(), this%rho_w%size())
    else
       call masked_gather_copy_0(this%nu%x, temp%x, this%msk, &
            temp%size(), this%nu%size())
       call masked_gather_copy_0(this%rho_w%x, this%rho%x, this%msk, &
            this%rho%size(), this%rho_w%size())
    end if

    call neko_scratch_registry%relinquish_field(idx)
  end subroutine duprat_compute_nu

  ! ===========================================================================
  ! Main compute
  ! ===========================================================================

  !> Compute wall shear stress using the Duprat (2011) model.
  !!
  !! **CPG path** (`use_constant_dpdx = .true.`):
  !!   Passes `dpdx_const` directly to the kernel. No field access.
  !!
  !! **APG path** (`use_constant_dpdx = .false.`):
  !!
  !!   Gated by `grad_tstep` to fire EXACTLY ONCE per timestep.
  !!   On the first call of step N:
  !!
  !!   Step 1 — compute exact dt:
  !!     if t_prev < 0: dt = t  (first ever call)
  !!     else:          dt = t - t_prev  (exact, not mean)
  !!
  !!   Step 2 — IIR filter update (only if t >= t_filter_start):
  !!     For each wall node i, using the ACTIVE gradient (= grad(p_{N-1})):
  !!       dpdx_local = grad_px[ind_r,ind_s,ind_t,ind_e] * (ui/magu)
  !!                  + grad_py[...] * (vi/magu)
  !!                  + grad_pz[...] * (wi/magu)
  !!       Clamp: |dpdx_local| <= 2*(magu*nu/h)^1.5/nu   [Stokes bound]
  !!       beta_f = min(dt/t_filter, 0.5)
  !!       dpdx_filt[i] = (1-beta_f)*dpdx_filt[i] + beta_f*dpdx_local
  !!
  !!     If t < t_filter_start: dpdx_filt stays at 0 (ZPG, no update).
  !!
  !!   Step 3 — update double-buffer (ALWAYS, regardless of t_filter_start):
  !!     Compute dudxyz(p_N) -> grad_px_buf/py_buf/pz_buf   [current step's grad]
  !!     grad_px <- grad_px_buf                              [promote: now p_{N-1} for N+1]
  !!     grad_py <- grad_py_buf
  !!     grad_pz <- grad_pz_buf
  !!     t_prev = t,  grad_tstep = tstep
  !!
  !!   Note: `compute_mag_field()` is NOT called here; wall_model_bc does it.
  subroutine duprat_compute(this, t, tstep)
    class(duprat_t), intent(inout) :: this
    real(kind=rp),   intent(in)    :: t
    integer,         intent(in)    :: tstep
    type(field_t), pointer :: u, v, w, p
    integer       :: i
    real(kind=rp) :: ui, vi, wi, normu, magu, dpdx_local, beta_f, dt
    real(kind=rp) :: max_dpdx_local

    call this%compute_nu()

    u => neko_registry%get_field("u")
    v => neko_registry%get_field("v")
    w => neko_registry%get_field("w")

    if (NEKO_BCKND_DEVICE .eq. 1) then
       call neko_error("Duprat GPU kernel not yet implemented")
    end if

    if (this%use_constant_dpdx) then

       ! ---- CPG / ZPG path --------------------------------------------------
       call duprat_compute_cpu( &
            u%x, v%x, w%x, &
            this%ind_r, this%ind_s, this%ind_t, this%ind_e, &
            this%n_x%x, this%n_y%x, this%n_z%x, &
            this%nu%x, this%rho_w%x, this%h%x, &
            this%tau_x%x, this%tau_y%x, this%tau_z%x, &
            this%n_nodes, u%Xh%lx, u%msh%nelv, &
            this%kappa, this%beta, this%A, this%dpdx_const)

    else

       ! ---- APG path --------------------------------------------------------
       !
       ! Gate: fire ONCE per timestep (first Krylov call only).
       ! On subsequent calls: reuse dpdx_filt unchanged.
       if (tstep .ne. this%grad_tstep) then

          ! --- Step 1: exact dt ---
          if (this%t_prev < 0.0_rp) then
             dt = t             ! very first call ever; t ≈ first dt
          else
             dt = t - this%t_prev
          end if
          dt = max(dt, 1.0e-14_rp)   ! guard against zero or negative dt

          ! --- Step 2: IIR filter update ---
          ! Uses the ACTIVE gradient arrays = grad(p_{N-1}).
          ! At the very first call (tstep=1), the active arrays are zero
          ! (initialised in finalize/init_from_components), so dpdx_local=0
          ! and dpdx_filt stays 0. This is correct: ZPG on the first step.
          if (t >= this%t_filter_start) then

             beta_f = min(dt / max(this%t_filter, 1.0e-14_rp), 0.5_rp)

             do i = 1, this%n_nodes
                ! Sample velocity at off-wall point; project out normal component.
                ui = u%x(this%ind_r(i), this%ind_s(i), this%ind_t(i), this%ind_e(i))
                vi = v%x(this%ind_r(i), this%ind_s(i), this%ind_t(i), this%ind_e(i))
                wi = w%x(this%ind_r(i), this%ind_s(i), this%ind_t(i), this%ind_e(i))
                normu = ui*this%n_x%x(i) + vi*this%n_y%x(i) + wi*this%n_z%x(i)
                ui = ui - normu*this%n_x%x(i)
                vi = vi - normu*this%n_y%x(i)
                wi = wi - normu*this%n_z%x(i)
                magu = sqrt(ui**2 + vi**2 + wi**2)

                if (magu <= 1.0e-14_rp) cycle  ! leave dpdx_filt unchanged

                ! Project PREVIOUS step's gradient onto tangential direction.
                dpdx_local = &
                     this%grad_px( &
                          this%ind_r(i),this%ind_s(i),this%ind_t(i),this%ind_e(i)) &
                     * (ui/magu) &
                     + this%grad_py( &
                          this%ind_r(i),this%ind_s(i),this%ind_t(i),this%ind_e(i)) &
                     * (vi/magu) &
                     + this%grad_pz( &
                          this%ind_r(i),this%ind_s(i),this%ind_t(i),this%ind_e(i)) &
                     * (wi/magu)

                ! Stokes self-consistent clamp (last-resort guard):
                ! Ensures u_P <= u_tau_stokes -> alpha >= 0.5 when clamped.
                max_dpdx_local = 2.0_rp &
                     * (magu * this%nu%x(i) / this%h%x(i))**1.5_rp &
                     / this%nu%x(i)
                dpdx_local = sign(min(abs(dpdx_local), max_dpdx_local), dpdx_local)

                ! IIR exponential moving average.
                this%dpdx_filt(i) = (1.0_rp - beta_f) * this%dpdx_filt(i) &
                                  + beta_f * dpdx_local
             end do

          else
             ! t < t_filter_start: zero the filter (ZPG, no update).
             this%dpdx_filt(:) = 0.0_rp
          end if

          ! --- Step 3: update double-buffer (ALWAYS) ---
          !
          ! Compute grad(p_N) -> buffer. Then promote buffer -> active.
          ! After this, grad_px/py/pz = grad(p_N) which becomes
          ! grad(p_{N-1}) for the next timestep's filter update.
          !
          ! This uses the CURRENT step's converged pressure (p_N), stored
          ! as the lagged gradient for step N+1.
          p => neko_registry%get_field("p")

          call dudxyz(this%grad_px_buf, p%x, &
               this%coef%drdx, this%coef%dsdx, this%coef%dtdx, this%coef)
          call dudxyz(this%grad_py_buf, p%x, &
               this%coef%drdy, this%coef%dsdy, this%coef%dtdy, this%coef)
          call dudxyz(this%grad_pz_buf, p%x, &
               this%coef%drdz, this%coef%dsdz, this%coef%dtdz, this%coef)

          ! Promote buffer -> active (simple array copy; avoids pointer aliasing).
          this%grad_px = this%grad_px_buf
          this%grad_py = this%grad_py_buf
          this%grad_pz = this%grad_pz_buf

          this%t_prev     = t
          this%grad_tstep = tstep

       end if  ! end of first-call-only block

       ! Call APG kernel with per-node filtered gradient.
       call duprat_compute_apg_cpu( &
            u%x, v%x, w%x, &
            this%ind_r, this%ind_s, this%ind_t, this%ind_e, &
            this%n_x%x, this%n_y%x, this%n_z%x, &
            this%nu%x, this%rho_w%x, this%h%x, &
            this%tau_x%x, this%tau_y%x, this%tau_z%x, &
            this%n_nodes, u%Xh%lx, u%msh%nelv, &
            this%kappa, this%beta, this%A, &
            this%dpdx_filt)

    end if

  end subroutine duprat_compute

end module duprat