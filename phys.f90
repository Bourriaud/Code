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
    a(2,1)=1.0_dp
    a(2,2)=0.0_dp
    
    f(1,:)=a(1,:)*u(1)
    f(2,:)=a(2,:)*u(2)

    deallocate(a)
    
    return
  end subroutine f_transport

  subroutine f_euler(u,f)
    real(dp), dimension(:), intent(in) :: u     !u=(rho,rhou,rhov,E)
    real(dp), dimension(:,:), intent(inout) :: f
    real(dp) :: p,gamma,ux,uy

    gamma=1.4_dp
    ux=u(2)/u(1)
    uy=u(3)/u(1)
    p=(u(4)-u(1)*(0.5*(ux**2+uy**2)))*(gamma-1.0_dp)

    ! x direction
    
    f(1,1)=u(2)
    f(2,1)=ux*u(2)+p
    f(3,1)=ux*u(3)
    f(4,1)=ux*(u(4)+p)

    ! y direction

    f(1,2)=u(3)
    f(2,2)=uy*u(2)
    f(3,2)=uy*u(3)+p
    f(4,2)=uy*(u(4)+p)

    return
  end subroutine f_euler
    

end module phys
