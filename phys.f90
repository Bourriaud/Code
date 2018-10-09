module phys

  implicit none

contains

  subroutine f_transport(u,f)
    real, dimension(:), intent(in) :: u
    real, dimension(:,:), intent(inout) :: f
    real, dimension(:,:), allocatable :: a

    allocate(a(size(u),2))
    
    ! Transport première variable
    a(1,1)=1.
    a(1,2)=0.
    ! Transport deuxième variable
    a(2,1)=1.
    a(2,2)=1.
    
    f(1,:)=a(1,:)*u(1)
    f(2,:)=a(2,:)*u(2)

    deallocate(a)
    
    return
  end subroutine f_transport

end module phys
