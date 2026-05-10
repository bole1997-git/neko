! Taylor-Green vortex (TGV)
!
! Time-integration of the Taylor-Green vortex up to specified time, evaluating
! total kinetic energy and enstrophy. The resulting flow fields may be
! visualised using Paraview or VisIt by opening the field0.nek5000 file.
!
module user
  use neko
  implicit none

  ! Global user variables
  type(field_t) :: w1
  real(kind=rp), parameter :: llx = 9.0_rp, lly = 3.035_rp, llz=4.5_rp
  real(kind=rp), parameter :: h_hill = 1.0_rp, w_hill = 1.929_rp

contains

  ! Register user-defined functions (see user_intf.f90)
  subroutine user_setup(user)
    type(user_t), intent(inout) :: user
    user%initial_conditions => initial_conditions
    user%mesh_setup => user_mesh_scale
    user%compute => user_calc_quantities
  end subroutine user_setup

  function math_ran_dst(ix, iy, iz, ieg, xl, fcoeff) result(rng)
      integer :: ix,iy,iz,ieg
      real(kind=rp) :: xl(3)
      real(kind=rp) :: fcoeff(3), rng

      rng = fcoeff(1)*(ieg + xl(1)*sin(xl(2))) + &
           fcoeff(2)*ix*iy + fcoeff(3)*ix
      rng = &
           fcoeff(1)*(ieg + xl(2)*sin(rng)) + &
           fcoeff(2)*iz*ix + fcoeff(3)*iz
      rng = 1.e3*sin(rng)
      rng = 1.e3*sin(rng)
      rng = cos(rng)
  end function math_ran_dst

  function hill_step(x, w, h)
    real(kind=rp), intent(in) :: x, w, h
    real(kind=rp) :: y, hill_step, xs

    xs = x/w

    if (xs.le.0) then
       y = h
    elseif (xs.gt.0.and.xs.le.9./54.) then
       y = h*min(1.0_rp,1.0_rp+7.05575248e-1*xs**2-1.1947737203e1*xs**3)
    elseif (xs.gt.9./54.and.xs.le.14./54.) then
       y = h*(0.895484248+1.881283544*xs-10.582126017*xs**2 &
            +10.627665327*xs**3)
    elseif (xs.gt.14./54.and.xs.le.20./54.) then
       y = h*(0.92128609+1.582719366*xs-9.430521329*xs**2 &
            +9.147030728*xs**3)
    elseif (xs.gt.20./54..and.xs.le.30./54.) then
       y = h*(1.445155365-2.660621763*xs+2.026499719*xs**2 &
            -1.164288215*xs**3)
    elseif (xs.gt.30./54..and.xs.le.40./54.) then
       y = h*(0.640164762+1.6863274926*xs-5.798008941*xs**2 &
            +3.530416981*xs**3)
    elseif (xs.gt.40./54..and.xs.le.1.) then
       y = h*(2.013932568-3.877432121*xs+1.713066537*xs**2 &
            +0.150433015*xs**3)
    else
       y = 0.
    endif

    hill_step = y
  end

  function hill_height(x, llx, w, h)
      real(kind=rp), intent(in) :: x, llx, w, h
      real(kind=rp) :: xx, hill_height

      if (x .lt. 0) then
         xx = llx + mod(x, llx)
      elseif (x .gt. llx) then
         xx = mod(x, llx)
      else
         xx = x
      endif
      hill_height = hill_step(xx, w, h) + hill_step(llx - xx, w, h)
  end

  function shift(x, y, llx, lly, w_hill) result(dx)
      real(kind=rp), intent(in) :: x, y, llx, lly, w_hill
      real(kind=rp) :: dx, xfac, yfac

      yfac = (1.0_rp - y/lly) ** 3
      xfac = 0.0_rp

      if (x .le. w_hill/2) then
         xfac = -2.0_rp / w_hill * x
      elseif (x .gt. w_hill/2 .and. x .le. llx - w_hill/2) then
         xfac = 2.0_rp / (llx - w_hill) * x -1.0_rp - w_hill / (llx - w_hill)
      elseif (x .gt. llx - w_hill/2) then
         xfac = -2.0_rp / w_hill * x + 2.0_rp * llx / w_hill
      endif

      dx = xfac*yfac
  end function

  ! Rescale mesh
  subroutine user_mesh_scale(mesh, time)
    !import mesh_t, time_state_t
    type(mesh_t), intent(inout) :: mesh
    type(time_state_t), intent(in) :: time
    integer :: i, p, nvert
    real(kind=rp) :: beta_x, beta_y, amp
    real(kind=rp) :: x, y, z, y_hill

    beta_x = 2.0
    beta_y = 1.5

    amp = 0.45

    nvert = size(mesh%points)


    do i = 1, nvert
       x = mesh%points(i)%x(1)
       y = mesh%points(i)%x(2)
       z = mesh%points(i)%x(3)

       ! Compress x resolution towards the center
       x = &
           0.5_rp * ( sinh(beta_x *(x - 0.5_rp)) / sinh(beta_x*0.5) + 1.0_rp)
       ! Compress y resolution towards the walls
       y = &
           0.5 * (tanh(beta_y * (2*y - 1.0_rp)) / tanh(beta_y) + 1.0_rp)

       ! Rescale
       x = x * llx
       y = y * lly
       z = z * llz

       ! Shift in x
       x = x + amp*shift(x, y, llx, lly, w_hill)

       mesh%apply_deform => phill



       !y_hill = hill_height(x, llx, w_hill, h_hill)
       !y = y_hill + y * (1 - y_hill / ly)

       mesh%points(i)%x(1) = x
       mesh%points(i)%x(2) = y
       mesh%points(i)%x(3) = z
    end do

  end subroutine user_mesh_scale

  subroutine phill(msh, x, y, z, lx, ly, lz)
    class(mesh_t) :: msh
    integer, intent(in) :: lx, ly, lz
    real(kind=rp), intent(inout) :: x(lx, lx, lx, msh%nelv)
    real(kind=rp), intent(inout) :: y(lx, lx, lx, msh%nelv)
    real(kind=rp), intent(inout) :: z(lx, lx, lx, msh%nelv)
    real(kind=rp) :: y_hill

    integer :: e, i, j ,k, l

    do e = 1,msh%nelv
       do k = 1, lz
          do j = 1, ly
              do i = 1, lx
                 y_hill = hill_height(x(i,j,k,e), llx, w_hill, h_hill)
                 y(i,j,k,e) = y_hill + y(i,j,k,e) * (1 - y_hill / lly)
              end do
          end do
       end do
    end do

  end subroutine phill

  ! User-defined initial condition
  ! subroutine initial_conditions(u, v, w, p, params)
  !  type(field_t), intent(inout) :: u
  !  type(field_t), intent(inout) :: v
  !  type(field_t), intent(inout) :: w
  !  type(field_t), intent(inout) :: p
  !  type(json_file), intent(inout) :: params
  !  integer :: i, ntot, idx(4), ix,iy,iz,ieg
  !  real(kind=rp) :: y_hill, x, y, z, fcoeff(3), xl(3), amp, ran
  !  real(kind=rp), parameter :: eps = 1e-3

  subroutine initial_conditions(scheme_name, fields)
    character(len=*), intent(in) :: scheme_name
    type(field_list_t), intent(inout) :: fields
    integer :: i, ntot, idx(4), ix,iy,iz,ieg
    real(kind=rp) :: uvw(3)
    real(kind=rp) :: y_hill, x, y, z, fcoeff(3), xl(3), amp, ran
    real(kind=rp), parameter :: eps = 1e-3
    type(dofmap_t), pointer :: dof
    type (field_t), pointer :: u, v, w, p

    dof => fields%dof(1)
    u => fields%get_by_name("u")
    v => fields%get_by_name("v")
    w => fields%get_by_name("w")
    p => fields%get_by_name("p")

    ntot = dof%size()
    
    !ntot = u%dof%size()
    do i = 1, ntot
       idx = nonlinear_index(i, u%dof%Xh%lx, u%dof%Xh%ly, u%dof%Xh%lz)
       x = u%dof%x(i,1,1,1)
       y = u%dof%y(i,1,1,1)
       z = u%dof%z(i,1,1,1)
       ix = idx(1)
       iy = idx(2)
       iz = idx(3)
       ieg = idx(4)

       xl = [x, y, z]

       !y_hill = hill_height(x, lx, w_hill, h_hill)
       !y = y_hill + y * (1 - y_hill / ly)
       !u%dof%y(i,1,1,1) = y

       amp = 0.2

       ran = 3.e4*(ieg+x*sin(y)+z*cos(y)) &
           + 4.7e2*ix*iy*iz - 1.5e3*ix*iy + .5e5*ix
       ran = 6.e3*sin(ran)
       ran = 3.e3*sin(ran)
       ran = cos(ran)
       u%x(i,1,1,1)= 1. + ran*amp

       ran = (2+ran)*1.e4*(ieg+y*sin(z)+x*cos(z)) &
           + 1.5e3*ix*iy*iz - 2.5e3*ix*iy + 8.9e4*ix
       ran = 2.e3*sin(ran)
       ran = 7.e3*sin(ran)
       ran = cos(ran)
       v%x(i,1,1,1) = ran*amp

       ran = (4+ran)*5.1e4*(ieg+z*sin(x)+y*cos(x)) &
           + 4.6e3*ix*iy*iz - 2.9e4*ix*iy + 3.7e3*ix
       ran = 9.e3*sin(ran)
       ran = 4.e3*sin(ran)
       ran = cos(ran)
       w%x(i,1,1,1) = ran*amp

       !fcoeff(1)=  3.0e4
       !fcoeff(2)= -1.5e3
       !fcoeff(3)=  0.5e5

       !u%x(i,1,1,1) = 1.0_rp + &
       !    eps * math_ran_dst(idx(1), idx(2), idx(3), idx(4), xl, fcoeff)

       !fcoeff(1)=  2.3e4
       !fcoeff(2)=  2.3e3
       !fcoeff(3)= -2.0e5

       !v%x(i,1,1,1) = &
       !    eps * math_ran_dst(idx(1), idx(2), idx(3), idx(4), xl, fcoeff)


       !fcoeff(1)= 2.e4
       !fcoeff(2)= 1.e3
       !fcoeff(3)= 1.e5

       !w%x(i,1,1,1) = &
       !    eps * math_ran_dst(idx(1), idx(2), idx(3), idx(4), xl, fcoeff)

    end do

    !p = 0._rp
    call field_rzero(p)
  end subroutine initial_conditions

  ! User-defined routine called at the end of every time step
  !subroutine user_calc_quantities(t, tstep, u, v, w, p, coef, params)
   ! real(kind=rp), intent(in) :: t
   ! integer, intent(in) :: tstep
   ! type(coef_t), intent(inout) :: coef
   ! type(json_file), intent(inout) :: params
   ! type(field_t), intent(inout) :: u
   ! type(field_t), intent(inout) :: v
   ! type(field_t), intent(inout) :: w
   ! type(field_t), intent(inout) :: p
   ! type(field_t), pointer :: omega_x, omega_y, omega_z
   ! integer :: ntot, i
   ! real(kind=rp) :: vv, sum_e1(1), e1, e2, sum_e2(1), oo, e3

  !end subroutine user_calc_quantities

    ! User-defined routine called at the end of every time step
  subroutine user_calc_quantities(time)
    type(time_state_t), intent(in) :: time

  end subroutine user_calc_quantities
     
end module user
