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
!! Reference: Duprat, C., Balarac, G., Metais, O., Congedo, P. M., and
!!   Brugiere, O. (2011). Physics of Fluids, 23(1), 015101.
!!
!! ## Pressure gradient modes
!!
!! CPG (`use_constant_dpdx: true`, default):
!!   Uniform scalar `dpdx_constant` [Pa/m] at every wall node.
!!   Set to 0.0 for ZPG warmup.
!!
!! APG (`use_constant_dpdx: false`):
!!   Per-node temporally filtered wall-tangential pressure gradient
!!   magnitude with sign, computed from dudxyz(p) corrected by jac.
!!
!! ## APG double-buffer design
!!
!!   grad_px/py/pz     = physical gradient dp/dx from p_{N-1}
!!   grad_px_buf/...   = physical gradient dp/dx from p_N (buffer)
!!
!!   Both arrays are zero for t < t_filter_start - t_filter.
!!
!!   At step N's first Krylov call (t >= t_filter_start - t_filter):
!!     1. Use grad_px (= p_{N-1}) for IIR filter update.
!!     2. Compute dudxyz(p_N)*jac -> grad_px_buf.
!!     3. Promote: grad_px <- grad_px_buf.
!!     4. Set grad_tstep = tstep, t_prev = t.
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
     real(kind=rp) :: t_filter_start = 5.0_rp
     real(kind=rp) :: t_filter       = 1.0_rp
     ! Nu and rho at wall nodes
     type(vector_t) :: nu
     type(vector_t) :: rho_w   ! kept for completeness; not passed to kernel
     ! Per-node IIR filtered wall-tangential |dP/ds| with sign [Pa/m].
     ! Zero until t >= t_filter_start.
     real(kind=rp), allocatable :: dpdx_filt(:)
     ! ACTIVE gradient: physical dp/dx from p_{N-1}, used at step N's filter.
     ! Zero for t < t_filter_start - t_filter
     real(kind=rp), allocatable :: grad_px(:,:,:,:)
     real(kind=rp), allocatable :: grad_py(:,:,:,:)
     real(kind=rp), allocatable :: grad_pz(:,:,:,:)
     ! BUFFER gradient: physical dp/dx from p_N, promoted at step N+1.
     ! Zero for t < t_filter_start - t_filter
     real(kind=rp), allocatable :: grad_px_buf(:,:,:,:)
     real(kind=rp), allocatable :: grad_py_buf(:,:,:,:)
     real(kind=rp), allocatable :: grad_pz_buf(:,:,:,:)
     ! Krylov guard: fires once per timestep only.
     integer        :: grad_tstep = 0
     ! t_prev for exact dt. -1 = "first call ever".
     real(kind=rp) :: t_prev = -1.0_rp
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
    call json_get_or_default(json, "use_constant_dpdx", use_constant_dpdx, &
         .true.)
    call json_get_or_default(json, "dpdx_constant",  dpdx_const,     0.0_rp)
    call json_get_or_default(json, "t_filter_start", t_filter_start, 5.0_rp)
    call json_get_or_default(json, "t_filter",       t_filter,       1.0_rp)

    call this%init_from_components(scheme_name, coef, msk, facet, h_index, &
         kappa, beta, A, use_constant_dpdx, dpdx_const, t_filter_start, &
         t_filter)
  end subroutine duprat_init

  ! ---------------------------------------------------------------------------

  subroutine duprat_partial_init(this, coef, json)
    class(duprat_t), intent(inout) :: this
    type(coef_t),    intent(in)    :: coef
    type(json_file), intent(inout) :: json

    call this%partial_init_base(coef, json)
    call json_get_or_lookup(json, "kappa", this%kappa)
    call json_get_or_lookup(json, "beta",  this%beta)
    call json_get_or_lookup(json, "A",     this%A)
    call json_get_or_default(json, "use_constant_dpdx", &
         this%use_constant_dpdx, .true.)
    call json_get_or_default(json, "dpdx_constant",  this%dpdx_const,     &
         0.0_rp)
    call json_get_or_default(json, "t_filter_start", this%t_filter_start, &
         5.0_rp)
    call json_get_or_default(json, "t_filter",       this%t_filter,       &
         1.0_rp)
  end subroutine duprat_partial_init

  ! ---------------------------------------------------------------------------

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
       this%dpdx_filt   = 0.0_rp
       this%grad_px     = 0.0_rp
       this%grad_py     = 0.0_rp
       this%grad_pz     = 0.0_rp
       this%grad_px_buf = 0.0_rp
       this%grad_py_buf = 0.0_rp
       this%grad_pz_buf = 0.0_rp
    end if
  end subroutine duprat_finalize

  ! ---------------------------------------------------------------------------

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
       this%dpdx_filt   = 0.0_rp
       this%grad_px     = 0.0_rp
       this%grad_py     = 0.0_rp
       this%grad_pz     = 0.0_rp
       this%grad_px_buf = 0.0_rp
       this%grad_py_buf = 0.0_rp
       this%grad_pz_buf = 0.0_rp
    end if
  end subroutine duprat_init_from_components

  ! ---------------------------------------------------------------------------

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
       call device_masked_gather_copy_0(this%nu%x_d, temp%x_d, &
            this%msk_d, temp%size(), this%nu%size())
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

  subroutine duprat_compute(this, t, tstep)
    class(duprat_t), intent(inout) :: this
    real(kind=rp),   intent(in)    :: t
    integer,         intent(in)    :: tstep

    type(field_t), pointer :: u, v, w, p
    integer       :: i
    real(kind=rp) :: ui, vi, wi, normu, magu
    real(kind=rp) :: dpx, dpy, dpz, dp_n
    real(kind=rp) :: dp_tx, dp_ty, dp_tz
    real(kind=rp) :: dpdx_local, max_dpdx, beta_f, dt
    ! linear ramp factor for smooth APG activation.
    real(kind=rp) :: ramp

    call this%compute_nu()

    u => neko_registry%get_field("u")
    v => neko_registry%get_field("v")
    w => neko_registry%get_field("w")

    if (NEKO_BCKND_DEVICE .eq. 1) then
       call neko_error("Duprat GPU kernel not yet implemented")
    end if

    ! =========================================================================
    ! CPG / ZPG path
    ! =========================================================================
    if (this%use_constant_dpdx) then

       call duprat_compute_cpu( &
            u%x, v%x, w%x, &
            this%ind_r, this%ind_s, this%ind_t, this%ind_e, &
            this%n_x%x, this%n_y%x, this%n_z%x, &
            this%nu%x, this%h%x, &
            this%tau_x%x, this%tau_y%x, this%tau_z%x, &
            this%n_nodes, u%Xh%lx, u%msh%nelv, &
            this%kappa, this%beta, this%A, this%dpdx_const, tstep)

    ! =========================================================================
    ! APG path
    ! =========================================================================
    else

       ! Krylov guard: first call per timestep only.
       if (tstep .ne. this%grad_tstep) then

          ! -------------------------------------------------------------------
          ! restart-safe exact dt.
          ! t_prev = -1 on very first call -> dt = 0 -> skip filter update.
          ! dpdx_filt stays 0 (ZPG), buffer populated for next step.
          ! -------------------------------------------------------------------
          if (this%t_prev < 0.0_rp) then
             dt = 0.0_rp
          else
             dt = t - this%t_prev
          end if
          dt = max(dt, 0.0_rp)

          ! -------------------------------------------------------------------
          ! IIR filter update.
          ! Fires only when t >= t_filter_start AND dt > 0.
          ! Uses ACTIVE gradient = physical dp/dx from p_{N-1}.
          ! -------------------------------------------------------------------
          if (t >= this%t_filter_start .and. dt > 1.0e-14_rp) then

             ! linear ramp of beta_f over first t_filter window.
             ! ramp = 0 at t = t_filter_start  (no update yet)
             ! ramp = 1 at t = t_filter_start + t_filter  (full weight)
             ramp   = min((t - this%t_filter_start) &
                          / max(this%t_filter, 1.0e-14_rp), 1.0_rp)
             ramp   = max(ramp, 0.0_rp)
             beta_f = min(dt / max(this%t_filter, 1.0e-14_rp), 0.5_rp) &
                    * ramp

             do i = 1, this%n_nodes

                ! Sample velocity at sampling point; remove wall-normal part.
                ui    = u%x(this%ind_r(i), this%ind_s(i), &
                            this%ind_t(i), this%ind_e(i))
                vi    = v%x(this%ind_r(i), this%ind_s(i), &
                            this%ind_t(i), this%ind_e(i))
                wi    = w%x(this%ind_r(i), this%ind_s(i), &
                            this%ind_t(i), this%ind_e(i))
                normu = ui*this%n_x%x(i) + vi*this%n_y%x(i) + &
                        wi*this%n_z%x(i)
                ui    = ui - normu*this%n_x%x(i)
                vi    = vi - normu*this%n_y%x(i)
                wi    = wi - normu*this%n_z%x(i)
                magu  = sqrt(ui**2 + vi**2 + wi**2)

                if (magu <= 1.0e-14_rp) cycle

                ! Sample PREVIOUS-step physical gradient at sampling point.
                ! grad_px contains dp/dx (after Jacobian correction).
                dpx = this%grad_px( &
                     this%ind_r(i), this%ind_s(i), &
                     this%ind_t(i), this%ind_e(i))
                dpy = this%grad_py( &
                     this%ind_r(i), this%ind_s(i), &
                     this%ind_t(i), this%ind_e(i))
                dpz = this%grad_pz( &
                     this%ind_r(i), this%ind_s(i), &
                     this%ind_t(i), this%ind_e(i))

                ! wall-tangential gradient magnitude + physical sign.
                !
                ! Step 2a: remove wall-normal component from grad_p.
                dp_n  = dpx*this%n_x%x(i) + dpy*this%n_y%x(i) + &
                        dpz*this%n_z%x(i)
                dp_tx = dpx - dp_n*this%n_x%x(i)
                dp_ty = dpy - dp_n*this%n_y%x(i)
                dp_tz = dpz - dp_n*this%n_z%x(i)

                ! Step 2b: magnitude of wall-tangential gradient vector.
                dpdx_local = sqrt(dp_tx**2 + dp_ty**2 + dp_tz**2)

                ! Step 2c: sign from dot of tangential gradient with
                ! wall-parallel velocity direction.
                ! Positive = APG (decelerating), negative = FPG.
                ! Physically correct even in recirculation zones.
                dpdx_local = dpdx_local * &
                     sign(1.0_rp, dp_tx*(ui/magu) + dp_ty*(vi/magu) + &
                                  dp_tz*(wi/magu))

                ! Stokes clamp applied ONCE here (not in kernel).
                ! Ensures u_P <= u_tau_stokes -> alpha >= 0.5 when clamped.
                max_dpdx = 2.0_rp * &
                     (magu * this%nu%x(i) / this%h%x(i))**1.5_rp &
                     / this%nu%x(i)
                dpdx_local = sign( &
                     min(abs(dpdx_local), max_dpdx), dpdx_local)

                ! IIR exponential moving average.
                this%dpdx_filt(i) = (1.0_rp - beta_f)*this%dpdx_filt(i) &
                                  + beta_f * dpdx_local

             end do

          else if (t < this%t_filter_start) then
             ! Before filter start: keep ZPG (zero), no update.
             this%dpdx_filt(:) = 0.0_rp
          end if
          ! If t >= t_filter_start but dt = 0 (first call after restart):
          ! skip update, reuse last valid dpdx_filt unchanged.

          ! -------------------------------------------------------------------
          ! Double-buffer update with Jacobian correction.
          !
          ! only run dudxyz within one t_filter window before
          ! t_filter_start and thereafter. Before that window, keep zero.
          !
          ! dudxyz returns (1/J)*dp/dx due to internal jacinv
          ! multiplication in cpu_dudxyz. Multiply by coef%jac to recover
          ! the true physical gradient dp/dx. Without this correction, on
          ! curved walls (periodic hill) the gradient is wrong by J ~ 2-5x,
          ! causing wrong u_P and CFL blow-up when APG activates.
          ! -------------------------------------------------------------------
          if (t >= this%t_filter_start - this%t_filter) then

             p => neko_registry%get_field("p")

             ! Compute weak gradient (1/J)*dp/dx over entire domain
             call dudxyz(this%grad_px_buf, p%x, &
                  this%coef%drdx, this%coef%dsdx, this%coef%dtdx, this%coef)
             call dudxyz(this%grad_py_buf, p%x, &
                  this%coef%drdy, this%coef%dsdy, this%coef%dtdy, this%coef)
             call dudxyz(this%grad_pz_buf, p%x, &
                  this%coef%drdz, this%coef%dsdz, this%coef%dtdz, this%coef)

             ! multiply by J to get true physical gradient dp/dx.
             ! coef%jac has shape (lx, ly, lz, nelv) — same as grad buffers.
             this%grad_px_buf = this%grad_px_buf * this%coef%jac
             this%grad_py_buf = this%grad_py_buf * this%coef%jac
             this%grad_pz_buf = this%grad_pz_buf * this%coef%jac

             ! Promote buffer -> active.
             this%grad_px = this%grad_px_buf
             this%grad_py = this%grad_py_buf
             this%grad_pz = this%grad_pz_buf

          else
             ! Before pre-activation window: zero all gradient arrays.
             ! First filter update sees grad_p = 0 and dpdx_filt ramps
             ! up cleanly from zero.
             this%grad_px     = 0.0_rp
             this%grad_py     = 0.0_rp
             this%grad_pz     = 0.0_rp
             this%grad_px_buf = 0.0_rp
             this%grad_py_buf = 0.0_rp
             this%grad_pz_buf = 0.0_rp
          end if

          ! record current time for next step's exact dt.
          this%t_prev     = t
          this%grad_tstep = tstep

       end if  ! end of first-call-only block

       ! -----------------------------------------------------------------------
       ! APG kernel: per-node pre-clamped filtered gradient.
       ! dpdx_filt is zero before t_filter_start and ramps up smoothly after.
       ! -----------------------------------------------------------------------
       call duprat_compute_apg_cpu( &
            u%x, v%x, w%x, &
            this%ind_r, this%ind_s, this%ind_t, this%ind_e, &
            this%n_x%x, this%n_y%x, this%n_z%x, &
            this%nu%x, this%h%x, &
            this%tau_x%x, this%tau_y%x, this%tau_z%x, &
            this%n_nodes, u%Xh%lx, u%msh%nelv, &
            this%kappa, this%beta, this%A, &
            this%dpdx_filt, tstep)

    end if

  end subroutine duprat_compute

end module duprat