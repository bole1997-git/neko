! Copyright (c) 2024, The Neko Authors
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
!> Implements `reichardt_t`.
!! Wall model based on Reichardt's universal law of the wall (1951).
module reichardt
  use field, only: field_t
  use num_types, only : rp
  use json_module, only : json_file
  use coefs, only : coef_t
  use neko_config, only : NEKO_BCKND_DEVICE
  use wall_model, only : wall_model_t
  use registry, only : neko_registry
  use json_utils, only : json_get_or_lookup
  use field_math, only: field_invcol3
  use vector, only : vector_t
  use math, only: masked_gather_copy_0
  use device_math, only: device_masked_gather_copy_0
  use scratch_registry, only : neko_scratch_registry
  use reichardt_cpu, only : reichardt_compute_cpu

  implicit none
  private

  !> Wall model based on Reichardt's universal law of the wall.
  type, public, extends(wall_model_t) :: reichardt_t
     !> The von Karman coefficient.
     real(kind=rp) :: kappa = 0.41_rp
     !> The log-law intercept.
     real(kind=rp) :: B = 5.2_rp
     !> The kinematic viscosity.
     type(vector_t) :: nu
   contains
     !> Constructor from JSON.
     procedure, pass(this) :: init => reichardt_init
     !> Partial constructor from JSON.
     procedure, pass(this) :: partial_init => reichardt_partial_init
     !> Finalize the construction.
     procedure, pass(this) :: finalize => reichardt_finalize
     !> Constructor from components.
     procedure, pass(this) :: init_from_components => &
          reichardt_init_from_components
     !> Destructor.
     procedure, pass(this) :: free => reichardt_free
     !> Compute the kinematic viscosity at the wall.
     procedure, pass(this) :: compute_nu => reichardt_compute_nu
     !> Compute the wall shear stress.
     procedure, pass(this) :: compute => reichardt_compute
  end type reichardt_t

contains

  !> Constructor from JSON.
  subroutine reichardt_init(this, scheme_name, coef, msk, facet, h_index, json)
    class(reichardt_t), intent(inout) :: this
    character(len=*), intent(in) :: scheme_name
    type(coef_t), intent(in) :: coef
    integer, intent(in) :: msk(:)
    integer, intent(in) :: facet(:)
    integer, intent(in) :: h_index
    type(json_file), intent(inout) :: json
    real(kind=rp) :: kappa, B

    call json_get_or_lookup(json, "kappa", kappa)
    call json_get_or_lookup(json, "B", B)

    call this%init_from_components(scheme_name, coef, msk, facet, h_index, &
         kappa, B)
  end subroutine reichardt_init

  !> Partial constructor from JSON.
  subroutine reichardt_partial_init(this, coef, json)
    class(reichardt_t), intent(inout) :: this
    type(coef_t), intent(in) :: coef
    type(json_file), intent(inout) :: json

    call this%partial_init_base(coef, json)
    call json_get_or_lookup(json, "kappa", this%kappa)
    call json_get_or_lookup(json, "B", this%B)

  end subroutine reichardt_partial_init

  !> Finalize the construction using the mask and facet arrays of the bc.
  subroutine reichardt_finalize(this, msk, facet)
    class(reichardt_t), intent(inout) :: this
    integer, intent(in) :: msk(:)
    integer, intent(in) :: facet(:)

    call this%finalize_base(msk, facet)
    call this%nu%init(this%n_nodes)
  end subroutine reichardt_finalize

  !> Constructor from components.
  subroutine reichardt_init_from_components(this, scheme_name, coef, msk, &
       facet, h_index, kappa, B)
    class(reichardt_t), intent(inout) :: this
    character(len=*), intent(in) :: scheme_name
    type(coef_t), intent(in) :: coef
    integer, intent(in) :: msk(:)
    integer, intent(in) :: facet(:)
    integer, intent(in) :: h_index
    real(kind=rp), intent(in) :: kappa
    real(kind=rp), intent(in) :: B

    call this%free()
    call this%init_base(scheme_name, coef, msk, facet, h_index)

    this%kappa = kappa
    this%B = B

    call this%nu%init(this%n_nodes)
  end subroutine reichardt_init_from_components

  !> Compute the kinematic viscosity vector.
  subroutine reichardt_compute_nu(this)
    class(reichardt_t), intent(inout) :: this
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
  end subroutine reichardt_compute_nu

  !> Destructor for the reichardt_t class.
  subroutine reichardt_free(this)
    class(reichardt_t), intent(inout) :: this
    call this%free_base()
  end subroutine reichardt_free

  !> Compute the wall shear stress using Reichardt law.
  subroutine reichardt_compute(this, t, tstep)
    class(reichardt_t), intent(inout) :: this
    real(kind=rp), intent(in) :: t
    integer, intent(in) :: tstep
    type(field_t), pointer :: u
    type(field_t), pointer :: v
    type(field_t), pointer :: w

    call this%compute_nu()

    u => neko_registry%get_field("u")
    v => neko_registry%get_field("v")
    w => neko_registry%get_field("w")

    if (NEKO_BCKND_DEVICE .eq. 1) then
       ! GPU implementation - placeholder
       error stop "Reichardt GPU kernel not yet implemented"
    else
       ! CPU implementation - call reichardt_compute_cpu
       call reichardt_compute_cpu(u%x, v%x, w%x, &
            this%ind_r, this%ind_s, this%ind_t, this%ind_e, &
            this%n_x%x, this%n_y%x, this%n_z%x, &
            this%nu%x, this%h%x, &
            this%tau_x%x, this%tau_y%x, this%tau_z%x, &
            this%n_nodes, u%Xh%lx, u%msh%nelv, &
            this%kappa, this%B, tstep)
    end if

    call this%compute_mag_field()

  end subroutine reichardt_compute

end module reichardt