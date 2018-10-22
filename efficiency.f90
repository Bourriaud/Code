module efficiency

  use types
  
  implicit none

contains

  subroutine exactSol(t,mesh,tabSol)
    real, intent(in) :: t
    type (meshStruct) :: mesh
    real, dimension(:,:), intent(inout) :: tabSol
    integer :: k
    real :: a1,a2,x,y

    a1=1.
    a2=1.
    do k=1,mesh%nc
       x=mesh%cell(k)%xc
       y=mesh%cell(k)%yc
       tabsol(k,1)=exp(-(x-a1*t-5.)**2-(y-a2*t-5.)**2)
    enddo
    
    return
  end subroutine exactSol

  subroutine userSol(t,mesh,sol)
    real, intent(in) :: t
    type(meshStruct), intent(in) :: mesh
    type(solStruct), intent(inout) :: sol

    call exactSol(t,mesh,sol%user(:,1:1))

    return
  end subroutine userSol

  subroutine errorL1(mesh,sol,exactsol,eL1)
    type(meshStruct), intent(in) :: mesh
    real, dimension(:), intent(in) :: sol,exactsol
    real, intent(out) :: eL1
    integer :: k
    real :: dx,dy

    eL1=0.
    do k=1,size(sol(:))
       dx=mesh%cell(k)%dx
       dy=mesh%cell(k)%dy
       eL1=eL1+abs(sol(k)-exactsol(k))*dx*dy
    enddo

    return
  end subroutine errorL1

  subroutine errorL2(mesh,sol,exactsol,eL2)
    type(meshStruct), intent(in) :: mesh
    real, dimension(:), intent(in) :: sol,exactsol
    real, intent(out) :: eL2
    integer :: k
    real :: dx,dy

    eL2=0.
    do k=1,size(sol(:))
       dx=mesh%cell(k)%dx
       dy=mesh%cell(k)%dy
       eL2=eL2+dx*dy*(sol(k)-exactsol(k))**2
    enddo
    eL2=eL2**0.5

    return
  end subroutine errorL2

end module efficiency
