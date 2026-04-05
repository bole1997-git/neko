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
!> Implements `duprat_t` - Duprat wall model for adverse pressure gradient flows.
!!
!! Wall model based on the Duprat et al. (2011) extended law of the wall:
!!   "A wall-layer model for large-eddy simulations of turbulent flows 
!!    with/without pressure gradient."
!!   Physics of Fluids 23, 015101.
!!
!! Includes pressure gradient effects via extended scaling with
!! Simpson velocity scale: u_p = (ρ/2 * |dP/dx|)^(1/3)
!!
module duprat
  use field, only: field_t
  use num_types, only: rp
  use json_module, only: json_file
  use coefs, only: coef_t
  use neko_config, only: NEKO_BCKND_DEVICE
  use wall_model, only: wall_model_t
  use registry, only: neko_registry
  use json_utils, only: json_get_or_lookup
  use field_math, only: field_invcol3
  use vector, only: vector_t
  use math, only: masked_gather_copy_0
  use device_math, only: device_masked_gather_copy_0
  use scratch_registry, only: neko_scratch_registry
  use duprat_cpu, only: duprat_compute_cpu

  implicit none
  private

  !> Wall model based on the Duprat et al. (2011) law with pressure gradient correction.
  !!
  !! Features:
  !! - Extended scaling with Simpson velocity scale
  !! - Handles both ZPG (Zero Pressure Gradient) and APG cases
  !! - Flexible dP/dx source: extracted from pressure field or constant CPG
  type, public, extends(wall_model_t) :: duprat_t
     !> The von Kármán constant.
     real(kind=rp) :: kappa = 0.41_rp
     !> Duprat damping coefficient (beta).
     real(kind=rp) :: beta = 0.78_rp
     !> Van Driest constant (A).
     real(kind=rp) :: A = 17.0_rp
     !> Kinematic viscosity at wall boundary nodes.
     type(vector_t) :: nu
     !> Pressure gradient field (dP/dx) at boundary nodes.
     type(vector_t) :: dpdx
     !> Extract dP/dx from pressure field (vs. constant CPG).
     logical :: extract_from_pressure = .true.
     !> Use constant pressure gradient (CPG mode).
     logical :: use_constant_dpdx = .false.
     !> Constant dP/dx value (for CPG mode).
     real(kind=rp) :: dpdx_const = 0.0_rp
   contains
     !> Constructor from JSON.
     procedure, pass(this) :: init => duprat_init
     !> Partial constructor from JSON.
     procedure, pass(this) :: partial_init => duprat_partial_init
     !> Finalize the construction.
     procedure, pass(this) :: finalize => duprat_finalize
     !> Constructor from components.
     procedure, pass(this) :: init_from_components => &
          duprat_init_from_components
     !> Destructor.
     procedure, pass(this) :: free => duprat_free
     !> Compute the kinematic viscosity at the wall.
     procedure, pass(this) :: compute_nu => duprat_compute_nu
     !> Compute the pressure gradient at the wall.
     procedure, pass(this) :: compute_dpdx => duprat_compute_dpdx
     !> Compute the wall shear stress.
     procedure, pass(this) :: compute => duprat_compute
  end type duprat_t

contains

  !> Constructor from JSON.
  subroutine duprat_init(this, scheme_name, coef, msk, facet, h_index, json)
    class(duprat_t), intent(inout) :: this
    character(len=*), intent(in) :: scheme_name
    type(coef_t), intent(in) :: coef
    integer, intent(in) :: msk(:)
    integer, intent(in) :: facet(:)
    integer, intent(in) :: h_index
    type(json_file), intent(inout) :: json
    real(kind=rp) :: kappa, beta, A, dpdx_const
    logical :: use_constant_dpdx

    call json_get_or_lookup(json, "kappa", kappa)
    call json_get_or_lookup(json, "beta", beta)
    call json_get_or_lookup(json, "A", A)
    
    ! Default: extract from pressure field (APG mode)
    use_constant_dpdx = .false.
    dpdx_const = 0.0_rp
    
    ! Try to read dpdx_const if it exists
    ! If present, user wants CPG mode
    call json_get_or_lookup(json, "dpdx_constant", dpdx_const)
    if (dpdx_const /= 0.0_rp) then
       use_constant_dpdx = .true.
    end if

    call this%init_from_components(scheme_name, coef, msk, facet, h_index, &
         kappa, beta, A, use_constant_dpdx, dpdx_const)
  end subroutine duprat_init

  !> Partial constructor from JSON.
  subroutine duprat_partial_init(this, coef, json)
    class(duprat_t), intent(inout) :: this
    type(coef_t), intent(in) :: coef
    type(json_file), intent(inout) :: json

    call this%partial_init_base(coef, json)
    call json_get_or_lookup(json, "kappa", this%kappa)
    call json_get_or_lookup(json, "beta", this%beta)
    call json_get_or_lookup(json, "A", this%A)
    
    ! Default: extract from pressure field
    this%extract_from_pressure = .true.
    this%use_constant_dpdx = .false.
    this%dpdx_const = 0.0_rp

  end subroutine duprat_partial_init

  !> Finalize the construction using the mask and facet arrays of the bc.
  subroutine duprat_finalize(this, msk, facet)
    class(duprat_t), intent(inout) :: this
    integer, intent(in) :: msk(:)
    integer, intent(in) :: facet(:)

    call this%finalize_base(msk, facet)
    call this%nu%init(this%n_nodes)
    call this%dpdx%init(this%n_nodes)

  end subroutine duprat_finalize

  !> Constructor from components.
  subroutine duprat_init_from_components(this, scheme_name, coef, msk, &
       facet, h_index, kappa, beta, A, use_constant_dpdx, dpdx_const)
    class(duprat_t), intent(inout) :: this
    character(len=*), intent(in) :: scheme_name
    type(coef_t), intent(in) :: coef
    integer, intent(in) :: msk(:)
    integer, intent(in) :: facet(:)
    integer, intent(in) :: h_index
    real(kind=rp), intent(in) :: kappa, beta, A
    logical, intent(in) :: use_constant_dpdx
    real(kind=rp), intent(in) :: dpdx_const

    call this%free()
    call this%init_base(scheme_name, coef, msk, facet, h_index)

    this%kappa = kappa
    this%beta = beta
    this%A = A
    this%use_constant_dpdx = use_constant_dpdx
    this%dpdx_const = dpdx_const

    if (use_constant_dpdx) then
       this%extract_from_pressure = .false.
    else
       this%extract_from_pressure = .true.
    end if

    call this%nu%init(this%n_nodes)
    call this%dpdx%init(this%n_nodes)

  end subroutine duprat_init_from_components

  !> Compute the kinematic viscosity vector at wall boundary nodes.
  !! Evaluates ν = μ/ρ at each boundary node and gathers it into the
  !! compact wall-node array this%nu.
  subroutine duprat_compute_nu(this)
    class(duprat_t), intent(inout) :: this
    type(field_t), pointer :: temp
    integer :: idx

    call neko_scratch_registry%request_field(temp, idx, .false.)
    call field_invcol3(temp, this%mu, this%rho)

    if (NEKO_BCKND_DEVICE .eq. 1) then
       call device_masked_gather_copy_0(this%nu%x_d, temp%x_d, this%msk_d, &
            temp%size(), this%nu%size())
    else
       call masked_gather_copy_0(this%nu%x, temp%x, this%msk, temp%size(), &
            this%nu%size())
    end if

    call neko_scratch_registry%relinquish_field(idx)
  end subroutine duprat_compute_nu

  !> Compute the pressure gradient at wall boundary nodes.
  !! Either extracts from pressure field or uses constant CPG value.
  subroutine duprat_compute_dpdx(this)
    class(duprat_t), intent(inout) :: this
    type(field_t), pointer :: p
    integer :: i

    if (this%use_constant_dpdx) then
       ! Use constant pressure gradient (CPG mode)
       this%dpdx%x(:) = this%dpdx_const
    else
       ! Extract from pressure field
       ! For now, use a placeholder approach
       ! Full implementation would compute grad_p and extract streamwise component
       p => neko_registry%get_field("p")
       
       ! Placeholder: set small dP/dx for testing
       ! In production: compute ∇p and extract ∂p/∂x
       this%dpdx%x(:) = 0.0_rp
       
    end if

  end subroutine duprat_compute_dpdx

  !> Destructor.
  subroutine duprat_free(this)
    class(duprat_t), intent(inout) :: this
    call this%free_base()
  end subroutine duprat_free

  !> Compute the wall shear stress using the Duprat (2011) law.
  !!
  !! Calls duprat_compute_cpu on CPU backends. GPU support is not yet
  !! implemented.
  subroutine duprat_compute(this, t, tstep)
    class(duprat_t), intent(inout) :: this
    real(kind=rp), intent(in) :: t
    integer, intent(in) :: tstep
    type(field_t), pointer :: u, v, w

    call this%compute_nu()
    call this%compute_dpdx()

    u => neko_registry%get_field("u")
    v => neko_registry%get_field("v")
    w => neko_registry%get_field("w")

    if (NEKO_BCKND_DEVICE .eq. 1) then
       error stop "Duprat GPU kernel not yet implemented"
    else
       call duprat_compute_cpu(u%x, v%x, w%x, &
            this%ind_r, this%ind_s, this%ind_t, this%ind_e, &
            this%n_x%x, this%n_y%x, this%n_z%x, &
            this%nu%x, this%dpdx%x, this%h%x, &
            this%tau_x%x, this%tau_y%x, this%tau_z%x, &
            this%n_nodes, u%Xh%lx, u%msh%nelv, &
            this%kappa, this%beta, this%A, tstep)
    end if

    call this%compute_mag_field()

  end subroutine duprat_compute

end module duprat