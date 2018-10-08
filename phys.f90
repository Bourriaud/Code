module phys

  implicit none

contains

  subroutine f_transport(u,f)
    real, dimension(:), intent(in) :: u
    real, dimension(:,:), intent(inout) :: f
    real, dimension(2) :: a

    a(1)=1.
    a(2)=1.
    f(:,1)=a(1)*u
    f(:,2)=a(2)*u

    return
  end subroutine f_transport

end module phys
