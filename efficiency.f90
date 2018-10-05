module efficiency

  use types
  
  implicit none

contains

  subroutine exactSol(mesh,tabSol)
    type (meshStruct) :: mesh
    real, dimension(0:,0:,:), intent(out) :: tabSol
    integer :: i,j

    do i=0,mesh%nx
       do j=0,mesh%ny
          tabSol(i,j,1)=real(i+j)
       enddo
    enddo
    
    return
  end subroutine exactSol

  subroutine userSol(mesh,sol)
    type(meshStruct), intent(in) :: mesh
    type(solStruct), intent(inout) :: sol
    integer :: i,j
    
    sol%nsolUser=1
    
    allocate(sol%user(0:mesh%nx,0:mesh%ny,sol%nsolUser),sol%nameUser(sol%nsolUser))

    sol%nameUser(1)='SolAnal'
    call exactSol(mesh,sol%user(:,:,1:sol%nsolUser))
    
    return
  end subroutine userSol

end module efficiency
