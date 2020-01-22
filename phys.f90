module phys

  use constant
  
  implicit none
  
contains

  subroutine f_transport(u,f)
    real(dp), dimension(:), intent(in) :: u
    real(dp), dimension(:,:), intent(inout) :: f
    real(dp), dimension(:,:), allocatable :: a

    allocate(a(size(u),2))
    
    ! Transport première variable
    a(1,1)=1.0_dp
    a(1,2)=1.0_dp
    ! Transport deuxième variable
    !a(2,1)=1.0_dp
    !a(2,2)=1.0_dp
    
    f(1,:)=a(1,:)*u(1)
    !f(2,:)=a(2,:)*u(2)

    deallocate(a)
    
    return
  end subroutine f_transport

  subroutine f_burgers(u,f)
    real(dp), dimension(:), intent(in) :: u
    real(dp), dimension(:,:), intent(inout) :: f
    
    f(1,:)=u(1)*u(1)/2.0_dp
    
    return
  end subroutine f_burgers

  subroutine f_euler(u,f)
    real(dp), dimension(:), intent(in) :: u     !u=(rho,rhou,rhov,E)
    real(dp), dimension(:,:), intent(inout) :: f
    real(dp) :: p,ux,uy

    call unconserv(u,"euler               ",2,ux)
    call unconserv(u,"euler               ",3,uy)
    call unconserv(u,"euler               ",4,p)

    ! x direction
    
    f(1,1)=u(2)
    f(2,1)=u(2)*ux+p
    f(3,1)=u(2)*uy
    f(4,1)=ux*(u(4)+p)

    ! y direction

    f(1,2)=u(3)
    f(2,2)=u(3)*ux
    f(3,2)=u(3)*uy+p
    f(4,2)=uy*(u(4)+p)

    return
  end subroutine f_euler

  subroutine f_euler_is(u,f)
    real(dp), dimension(:), intent(in) :: u     !u=(rho,rhou,rhov,E)
    real(dp), dimension(:,:), intent(inout) :: f
    real(dp) :: p,ux,uy

    call unconserv(u,"euler_is            ",2,ux)
    call unconserv(u,"euler_is            ",3,uy)
    call unconserv(u,"euler_is            ",4,p)

    ! x direction
    
    f(1,1)=u(2)
    f(2,1)=u(2)*ux+p
    f(3,1)=u(2)*uy

    ! y direction

    f(1,2)=u(3)
    f(2,2)=u(3)*ux
    f(3,2)=u(3)*uy+p

    return
  end subroutine f_euler_is

  subroutine f_M1(u,f)
    real(dp), dimension(:), intent(in) :: u     !u=(E,Fx,Fy)
    real(dp), dimension(:,:), intent(inout) :: f
    real(dp) :: x,f2

    f2=(u(2)/(c*u(1)))**2
    x=(3.0_dp+4.0_dp*f2)/(5.0_dp+2.0_dp*sqrt(4.0_dp-3.0_dp*f2))

    ! x direction
    
    f(1,1)=u(2)
    f(2,1)=(c**2)*u(1)*x
    f(3,1)=(c**2)*u(1)*(3.0_dp*x-1.0_dp)/2.0_dp

    ! y direction

    f(1,1)=u(3)
    f(2,1)=(c**2)*u(1)*(3.0_dp*x-1.0_dp)/2.0_dp
    f(3,1)=(c**2)*u(1)*x

    return
  end subroutine f_M1

  subroutine conserv(U,str_equa,isol,uc)
    real(dp), dimension(:), intent(in) :: U
    character(len=20), intent(in) :: str_equa
    integer, intent(in) :: isol
    real(dp), intent(out) :: uc

    select case (trim(str_equa))
    case("advection")
       uc=U(isol)
    case("euler")
       select case (isol)
       case(1)
          uc=U(1)
       case(2)
          uc=U(1)*U(2)
       case(3)
          uc=U(1)*U(3)
       case(4)
          uc=U(1)*(0.5_dp*(U(2)**2+U(3)**2)+U(4)/((gamma-1)*U(1)))
       case default
          print*, "Change conserv function"
          call exit()
       end select
    case("euler_is")
       select case (isol)
       case(1)
          uc=U(1)
       case(2)
          uc=U(1)*U(2)
       case(3)
          uc=U(1)*U(3)
       case(4)
          uc=U(1)*(0.5_dp*(U(2)**2+U(3)**2)+(U(1)**gamma)/((gamma-1)*U(1)))
       case default
          print*, "Change conserv function"
          call exit()
       end select
    case("M1")
       uc=U(isol)
    case default
       print*,"Conserv function not implemented"
       call exit()
    end select

    return
  end subroutine conserv

  subroutine unconserv(U,str_equa,isol,uu)
    real(dp), dimension(:), intent(in) :: U
    character(len=20), intent(in) :: str_equa
    integer, intent(in) :: isol
    real(dp), intent(out) :: uu

    select case (trim(str_equa))
    case("transport")
       uu=U(isol)
    case("euler")
       select case (isol)
       case(1)
          uu=U(1)
       case(2)
          uu=U(2)/U(1)
       case(3)
          uu=U(3)/U(1)
       case(4)
          uu=(gamma-1)*(U(4)-0.5_dp*(U(2)**2+U(3)**2)/U(1))
       case default
          print*, "Change unconserv function"
          call exit()
       end select
    case("euler_is")
       select case (isol)
       case(1)
          uu=U(1)
       case(2)
          uu=U(2)/U(1)
       case(3)
          uu=U(3)/U(1)
       case(4)
          uu=U(1)**gamma
       case default
          print*, "Change unconserv function"
          call exit()
       end select
    case("M1")
       uu=U(isol)
    case default
       print*,"Unconserv function not implemented"
       call exit()
    end select

    return
  end subroutine unconserv
    
end module phys
