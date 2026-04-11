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
!! Wall model based on Duprat et al. (2011) extended law of the wall
!! for turbulent flows with/without streamwise pressure gradient.
!!
!! ## Physics
!!
!!   u_P   = (nu * |dP/dx| / (2*rho))^(1/3)     [Simpson velocity scale; rho=1]
!!   u_p*  = sqrt(u_tau^2 + u_P^2)
!!   alpha = u_tau^2 / u_p*^2  in [0, 1]
!!   y*    = y * u_p*/nu,   U* = U/u_p*
!!   nu_t* = (kappa*y* + beta*y*(1-alpha)^1.5)*(1-exp(-y*/(1+A*alpha^3)))^2
!!   dU*/dy* = [sign(dP/dx)*(1-alpha)^1.5*y* + 1] / (1 + nu_t*)
!!
!! ## Pressure gradient modes
!!
!! **CPG** (`use_constant_dpdx: true`, default):
!!   Uses `dpdx_constant` [Pa/m] uniformly. Set to 0 for ZPG.
!!
!! **APG** (`use_constant_dpdx: false`):
!!   Extracts the VOLUME-AVERAGED streamwise pressure gradient from the
!!   resolved pressure field `p` at every timestep.
!!
!!   ### Why volume-average and not pointwise?
!!
!!   The pressure field p in LES contains both the mean gradient AND
!!   turbulent pressure fluctuations of magnitude p_rms ~ 3 u_tau^2.
!!   The fluctuation gradient nabla(p_fluct) has random signs per GLL node.
!!   For nodes where sign(dp_fluct/dx) = -1 (apparent FPG), the Duprat ODE
!!   numerator sign(dP/dx)*(1-alpha)^1.5*y* + 1 becomes negative for
!!   y* > 1/(1-alpha)^1.5 ~ 8, making U*(y*) non-monotone. The Newton
!!   solver then either diverges or converges to a spurious root, producing
!!   a wrong utau -> wrong tau -> velocity blow-up.
!!
!!   The MEAN pressure gradient has a definite, physically correct sign
!!   and does not cause this problem.
!!
!!   ### Formula (MPI-parallel, exact)
!!
!!   opgrad returns the weak-form (B-weighted) gradient:
!!     gx(i) = B_i * (dp/dx)_i,   B_i = J_i * w_i  (= coef%B)
!!
!!   The volume-weighted mean streamwise gradient is then:
!!     dP_mean/dx = sum_i(gx_i) / sum_i(B_i)
!!               = glsum(gx, n) / glsum(coef%B, n)
!!
!!   glsum performs MPI_Allreduce, giving the GLOBAL volume average.
!!   No pointwise division by B is needed — the B factors cancel.
!!
!!   ### Behaviour in a body-force-driven channel
!!
!!   When flow_rate_force drives the channel, the mean dP/dx is absorbed
!!   into the body force, so opgrad(p) gives fluctuations only ->
!!   dP_mean/dx ~ 0 -> u_P ~ 0 -> alpha ~ 1 -> ZPG Duprat (Van Driest).
!!   This is physically correct for a statistically stationary channel.
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
  use math, only: masked_gather_copy_0, glsum
  use device_math, only: device_masked_gather_copy_0
  use scratch_registry, only: neko_scratch_registry
  use operators, only: opgrad
  use duprat_cpu, only: duprat_compute_cpu
  use utils, only: neko_error

  implicit none
  private

  !> Duprat (2011) wall model.
  !!
  !! JSON keys:
  !!   kappa              (real, default 0.41)
  !!   beta               (real, default 0.78)
  !!   A                  (real, default 17.0)
  !!   use_constant_dpdx  (logical, default .true.)
  !!   dpdx_constant      (real, default 0.0 [Pa/m])
  !!
  type, public, extends(wall_model_t) :: duprat_t
     real(kind=rp) :: kappa = 0.41_rp
     real(kind=rp) :: beta  = 0.78_rp
     real(kind=rp) :: A     = 17.0_rp
     logical        :: use_constant_dpdx = .true.
     real(kind=rp) :: dpdx_const = 0.0_rp
     type(vector_t) :: nu
     type(vector_t) :: rho_w
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
    real(kind=rp) :: kappa, beta, A, dpdx_const
    logical       :: use_constant_dpdx

    call json_get_or_lookup(json, "kappa", kappa)
    call json_get_or_lookup(json, "beta",  beta)
    call json_get_or_lookup(json, "A",     A)
    call json_get_or_default(json, "use_constant_dpdx", use_constant_dpdx, .true.)
    call json_get_or_default(json, "dpdx_constant",     dpdx_const, 0.0_rp)

    call this%init_from_components(scheme_name, coef, msk, facet, h_index, &
         kappa, beta, A, use_constant_dpdx, dpdx_const)
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
    call json_get_or_default(json, "dpdx_constant",     this%dpdx_const, 0.0_rp)
  end subroutine duprat_partial_init

  subroutine duprat_finalize(this, msk, facet)
    class(duprat_t), intent(inout) :: this
    integer,         intent(in)    :: msk(:), facet(:)

    call this%finalize_base(msk, facet)
    call this%nu%init(this%n_nodes)
    call this%rho_w%init(this%n_nodes)
  end subroutine duprat_finalize

  subroutine duprat_init_from_components(this, scheme_name, coef, msk, &
       facet, h_index, kappa, beta, A, use_constant_dpdx, dpdx_const)
    class(duprat_t), intent(inout) :: this
    character(len=*), intent(in)   :: scheme_name
    type(coef_t),     intent(in)   :: coef
    integer,          intent(in)   :: msk(:), facet(:), h_index
    real(kind=rp),    intent(in)   :: kappa, beta, A, dpdx_const
    logical,          intent(in)   :: use_constant_dpdx

    call this%free()
    call this%init_base(scheme_name, coef, msk, facet, h_index)

    this%kappa             = kappa
    this%beta              = beta
    this%A                 = A
    this%use_constant_dpdx = use_constant_dpdx
    this%dpdx_const        = dpdx_const

    call this%nu%init(this%n_nodes)
    call this%rho_w%init(this%n_nodes)
  end subroutine duprat_init_from_components

  subroutine duprat_free(this)
    class(duprat_t), intent(inout) :: this
    call this%free_base()
    call this%nu%free()
    call this%rho_w%free()
  end subroutine duprat_free

  ! ===========================================================================
  ! Nu / rho gather
  ! ===========================================================================

  subroutine duprat_compute_nu(this)
    class(duprat_t), intent(inout) :: this
    type(field_t), pointer :: temp
    integer :: idx

    call neko_scratch_registry%request_field(temp, idx, .false.)
    call field_invcol3(temp, this%mu, this%rho)   ! temp = mu/rho = nu

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
  !! **CPG path**: passes `dpdx_const` directly to the kernel. Fast, no
  !! pressure field access.
  !!
  !! **APG path**: computes the global volume-averaged streamwise pressure
  !! gradient from the resolved pressure field `p`:
  !!
  !!   1. opgrad(gx, gy_dummy, gy_dummy, p, coef) -> gx = B*(dp/dx)
  !!      (we reuse gy_dummy for y and z outputs; those values are not used)
  !!   2. dP_mean/dx = glsum(gx, n) / glsum(coef%B, n)
  !!      Both glsum calls perform MPI_Allreduce => global result.
  !!   3. This single scalar is passed uniformly to every wall node.
  !!
  !! The volume-average eliminates turbulent fluctuation gradients that have
  !! random signs per node and would make the Duprat ODE non-monotone.
  !!
  !! Note: compute_mag_field() is NOT called here; wall_model_bc does it.
  subroutine duprat_compute(this, t, tstep)
    class(duprat_t), intent(inout) :: this
    real(kind=rp),   intent(in)    :: t
    integer,         intent(in)    :: tstep
    type(field_t), pointer :: u, v, w, p
    type(field_t), pointer :: gx, gy_dummy
    integer       :: idx_gx, idx_gy, n_total
    real(kind=rp) :: dpdx_apg

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
            this%kappa, this%beta, this%A, this%dpdx_const, tstep)

    else

       ! ---- APG path: volume-averaged mean streamwise pressure gradient ------
       !
       ! opgrad(gx, gy, gz, p, coef) computes the WEAK-form gradient:
       !   gx_i = B_i * (dp/dx)_i
       !   gy_i = B_i * (dp/dy)_i
       !   gz_i = B_i * (dp/dz)_i
       !
       ! We only need gx for the volume average. To avoid allocating a third
       ! scratch field for gz, we pass gy_dummy for both the y and z outputs.
       ! The Neko CPU opgrad backend writes each component in a separate
       ! loop, so the final content of gy_dummy is B*(dp/dz) (last write),
       ! which we discard. gx is unaffected.
       !
       ! Volume-averaged mean:
       !   dP_mean/dx = sum_i[B_i*(dp/dx)_i] / sum_i[B_i]
       !              = glsum(gx, n) / glsum(B, n)
       !
       ! glsum performs MPI_Allreduce: result is globally consistent.
       ! The B factors cancel exactly — no pointwise division needed.

       p => neko_registry%get_field("p")
       n_total = this%coef%Xh%lxyz * this%coef%msh%nelv

       call neko_scratch_registry%request_field(gx,       idx_gx,  .false.)
       call neko_scratch_registry%request_field(gy_dummy, idx_gy,   .false.)

       ! Weak-form gradient: gx <- B*(dp/dx),  gy_dummy <- B*(dp/dz) [discarded]
       call opgrad(gx%x, gy_dummy%x, gy_dummy%x, p%x, this%coef)

       ! Global volume-weighted mean streamwise pressure gradient
       dpdx_apg = glsum(gx%x, n_total) / glsum(this%coef%B, n_total)

       call neko_scratch_registry%relinquish_field(idx_gx)
       call neko_scratch_registry%relinquish_field(idx_gy)

       ! Pass the single scalar to the kernel (same as CPG path).
       call duprat_compute_cpu( &
            u%x, v%x, w%x, &
            this%ind_r, this%ind_s, this%ind_t, this%ind_e, &
            this%n_x%x, this%n_y%x, this%n_z%x, &
            this%nu%x, this%rho_w%x, this%h%x, &
            this%tau_x%x, this%tau_y%x, this%tau_z%x, &
            this%n_nodes, u%Xh%lx, u%msh%nelv, &
            this%kappa, this%beta, this%A, dpdx_apg, tstep)

    end if

  end subroutine duprat_compute

end module duprat