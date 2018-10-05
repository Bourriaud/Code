module efficiency

  use types
  
  implicit none

contains

  subroutine exactSol(mesh,isol)
    type (meshStruct) :: mesh
    real, dimension(0:,0:), intent(out) :: isol
    integer :: i,j

    do i=0,mesh%nx
       do j=0,mesh%ny
          isol(i,j)=i+j
       enddo
    enddo
    
    return
  end subroutine exactSol

end module efficiency
