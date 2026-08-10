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

  !> Wall model based on Duprat et al. (2011) extended law of the wall.
  !! Reference: https://doi.org/10.1063/1.3529358
  type, public, extends(wall_model_t) :: duprat_t
     !> The von Karman coefficient.
     real(kind=rp) :: kappa = 0.41_rp
     !> Pressure-gradient mixing-length exponent.
     real(kind=rp) :: beta  = 0.78_rp
     !> Van Driest damping constant.
     real(kind=rp) :: A     = 17.0_rp
     !> Use a constant streamwise dP/dx instead of the filtered field value.
     logical        :: use_constant_dpdx = .true.
     !> Constant streamwise pressure gradient, used if use_constant_dpdx.
     real(kind=rp) :: dpdx_const = 0.0_rp
     !> Simulation time at which the pressure-gradient filter starts.
     real(kind=rp) :: t_filter_start = 5.0_rp
     !> Time constant of the pressure-gradient filter.
     real(kind=rp) :: t_filter       = 1.0_rp
     !> The kinematic viscosity.
     type(vector_t) :: nu
     !> The fluid density at the boundary.
     type(vector_t) :: rho_w
     !> Per-node filtered wall-tangential pressure gradient.
     real(kind=rp), allocatable :: dpdx_filt(:)
     real(kind=rp), allocatable :: grad_px(:,:,:,:)
     real(kind=rp), allocatable :: grad_py(:,:,:,:)
     real(kind=rp), allocatable :: grad_pz(:,:,:,:)
     real(kind=rp), allocatable :: grad_px_buf(:,:,:,:)
     real(kind=rp), allocatable :: grad_py_buf(:,:,:,:)
     real(kind=rp), allocatable :: grad_pz_buf(:,:,:,:)
     integer        :: grad_tstep = 0
     real(kind=rp) :: t_prev = -1.0_rp
   contains
     !> Constructor from JSON.
     procedure, pass(this) :: init              => duprat_init
     !> Partial constructor from JSON.
     procedure, pass(this) :: partial_init      => duprat_partial_init
     !> Finalize the construction using the mask and facet arrays of the bc.
     procedure, pass(this) :: finalize          => duprat_finalize
     !> Constructor from components.
     procedure, pass(this) :: init_from_components => duprat_init_from_components
     !> Destructor.
     procedure, pass(this) :: free              => duprat_free
     !> Compute the kinematic viscosity and density at the wall.
     procedure, pass(this) :: compute_nu        => duprat_compute_nu
     !> Compute the wall shear stress.
     procedure, pass(this) :: compute           => duprat_compute
  end type duprat_t

contains

  !> Constructor from JSON.
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

  !> Partial constructor from JSON.
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

  !> Finalize the construction using the mask and facet arrays of the bc.
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

  !> Constructor from components.
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

  !> Destructor.
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

  !> Compute the kinematic viscosity and density at the wall.
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

  !> Compute the wall shear stress. GPU backend not yet implemented.
  !! @param t The time value.
  !! @param tstep The current time-step.
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
    real(kind=rp) :: ramp

    call this%compute_nu()

    u => neko_registry%get_field("u")
    v => neko_registry%get_field("v")
    w => neko_registry%get_field("w")

    if (NEKO_BCKND_DEVICE .eq. 1) then
       call neko_error("Duprat GPU kernel not yet implemented")
    end if

    if (this%use_constant_dpdx) then

       call duprat_compute_cpu( &
            u%x, v%x, w%x, &
            this%ind_r, this%ind_s, this%ind_t, this%ind_e, &
            this%n_x%x, this%n_y%x, this%n_z%x, &
            this%nu%x, this%rho_w%x, this%h%x, &
            this%tau_x%x, this%tau_y%x, this%tau_z%x, &
            this%n_nodes, u%Xh%lx, u%msh%nelv, &
            this%kappa, this%beta, this%A, this%dpdx_const, tstep)

    else

       ! Update the filtered pressure gradient once per time-step.
       if (tstep .ne. this%grad_tstep) then

          if (this%t_prev < 0.0_rp) then
             dt = 0.0_rp
          else
             dt = t - this%t_prev
          end if
          dt = max(dt, 0.0_rp)

          if (t >= this%t_filter_start .and. dt > 1.0e-14_rp) then

             ramp   = min((t - this%t_filter_start) &
                          / max(this%t_filter, 1.0e-14_rp), 1.0_rp)
             ramp   = max(ramp, 0.0_rp)
             beta_f = min(dt / max(this%t_filter, 1.0e-14_rp), 0.5_rp) &
                    * ramp

             do i = 1, this%n_nodes

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

                dpx = this%grad_px( &
                     this%ind_r(i), this%ind_s(i), &
                     this%ind_t(i), this%ind_e(i))
                dpy = this%grad_py( &
                     this%ind_r(i), this%ind_s(i), &
                     this%ind_t(i), this%ind_e(i))
                dpz = this%grad_pz( &
                     this%ind_r(i), this%ind_s(i), &
                     this%ind_t(i), this%ind_e(i))

                ! Remove the wall-normal component
                dp_n  = dpx*this%n_x%x(i) + dpy*this%n_y%x(i) + &
                        dpz*this%n_z%x(i)
                dp_tx = dpx - dp_n*this%n_x%x(i)
                dp_ty = dpy - dp_n*this%n_y%x(i)
                dp_tz = dpz - dp_n*this%n_z%x(i)

                dpdx_local = sqrt(dp_tx**2 + dp_ty**2 + dp_tz**2)

                ! Sign from the velocity-aligned component
                dpdx_local = dpdx_local * &
                     sign(1.0_rp, dp_tx*(ui/magu) + dp_ty*(vi/magu) + &
                                  dp_tz*(wi/magu))

                ! Stokes clamp on the pressure gradient
                max_dpdx = this%rho_w%x(i) * &
                     (magu * this%nu%x(i) / this%h%x(i))**1.5_rp &
                     / this%nu%x(i)
                dpdx_local = sign( &
                     min(abs(dpdx_local), max_dpdx), dpdx_local)

                this%dpdx_filt(i) = (1.0_rp - beta_f)*this%dpdx_filt(i) &
                                  + beta_f * dpdx_local

             end do

          else if (t < this%t_filter_start) then
             this%dpdx_filt(:) = 0.0_rp
          end if

          if (t >= this%t_filter_start - this%t_filter) then

             p => neko_registry%get_field("p")

             call dudxyz(this%grad_px_buf, p%x, &
                  this%coef%drdx, this%coef%dsdx, this%coef%dtdx, this%coef)
             call dudxyz(this%grad_py_buf, p%x, &
                  this%coef%drdy, this%coef%dsdy, this%coef%dtdy, this%coef)
             call dudxyz(this%grad_pz_buf, p%x, &
                  this%coef%drdz, this%coef%dsdz, this%coef%dtdz, this%coef)

             this%grad_px = this%grad_px_buf
             this%grad_py = this%grad_py_buf
             this%grad_pz = this%grad_pz_buf

          else
             this%grad_px     = 0.0_rp
             this%grad_py     = 0.0_rp
             this%grad_pz     = 0.0_rp
             this%grad_px_buf = 0.0_rp
             this%grad_py_buf = 0.0_rp
             this%grad_pz_buf = 0.0_rp
          end if

          this%t_prev     = t
          this%grad_tstep = tstep

       end if

       call duprat_compute_apg_cpu( &
            u%x, v%x, w%x, &
            this%ind_r, this%ind_s, this%ind_t, this%ind_e, &
            this%n_x%x, this%n_y%x, this%n_z%x, &
            this%nu%x, this%rho_w%x, this%h%x, &
            this%tau_x%x, this%tau_y%x, this%tau_z%x, &
            this%n_nodes, u%Xh%lx, u%msh%nelv, &
            this%kappa, this%beta, this%A, &
            this%dpdx_filt, tstep)

    end if

  end subroutine duprat_compute

end module duprat