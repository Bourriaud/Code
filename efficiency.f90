module efficiency

  use types
  
  implicit none

contains

  subroutine exactSol(t,mesh,tabSol)
    real, intent(in) :: t
    type (meshStruct) :: mesh
    real, dimension(:,:,:), intent(inout) :: tabSol
    integer :: i,j
    real :: a1,a2

    a1=1.
    a2=0.
    do i=1,mesh%nx
       do j=1,mesh%ny
          if (((mesh%CX(i,j)-a1*t)**2+(mesh%CY(i,j)-a2*t)**2)**0.5<5) then
             if((mesh%CX(i,j)-a1*t>0.).and.(mesh%CY(i,j)-a2*t>0.)) then
                tabsol(i,j,1)=1.
             else
                tabsol(i,j,1)=0.
             endif
          else
             tabsol(i,j,1)=0.
          endif
       enddo
    enddo
    
    return
  end subroutine exactSol

  subroutine userSol(t,mesh,sol)
    real, intent(in) :: t
    type(meshStruct), intent(in) :: mesh
    type(solStruct), intent(inout) :: sol
    integer :: i,j

    call exactSol(t,mesh,sol%user(:,:,1:1))

    return
  end subroutine userSol

  subroutine errorL1(sol,exactsol,eL1)
    real, dimension(:,:), intent(in) :: sol,exactsol
    real, intent(out) :: eL1
    integer :: i,j

    eL1=0.
    do i=1,size(sol(:,1))
       do j=1,size(sol(1,:))
          eL1=eL1+abs(sol(i,j)-exactsol(i,j))
       enddo
    enddo

    return
  end subroutine errorL1

  subroutine errorL2(sol,exactsol,eL2)
    real, dimension(:,:), intent(in) :: sol,exactsol
    real, intent(out) :: eL2
    integer :: i,j

    eL2=0.
    do i=1,size(sol(:,1))
       do j=1,size(sol(1,:))
          eL2=eL2+(sol(i,j)-exactsol(i,j))**2
       enddo
    enddo
    eL2=eL2**0.5

    return
  end subroutine errorL2

end module efficiency
