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

  subroutine conserv(U,str_equa,isol,uc)
    real(dp), dimension(:), intent(in) :: U
    character(len=20), intent(in) :: str_equa
    integer, intent(in) :: isol
    real(dp), intent(out) :: uc

    select case (trim(str_equa))
    case("transport")
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
    case default
       print*,"Unconserv function not implemented"
       call exit()
    end select

    return
  end subroutine unconserv
    
end module phys
