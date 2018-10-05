module phys

  implicit none

contains

  subroutine f_transport(u,f)
    real, intent(in) :: u
    real, dimension(2), intent(out) :: f
    real, dimension(2) :: a
    
    a=1.
    f=u*a
    
    return
  end subroutine f_transport

end module phys
