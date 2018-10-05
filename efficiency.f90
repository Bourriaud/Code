module efficiency

  implicit none

contains

  subroutine exactSol(Ix,Iy,isol)
    real, dimension(0:), intent(in) :: Ix,Iy
    real, dimension(0:,0:), intent(out) :: isol
    integer :: i,j

    do i=0,size(Ix)-1
       do j=0,size(Iy)-1
          isol(i,j)=i+j
       enddo
    enddo
    return
  end subroutine exactSol

end module efficiency
