module efficiency

  use types
  
  implicit none

contains

  subroutine exactSol(t,mesh,tabSol)
    real, intent(in) :: t
    type (meshStruct) :: mesh
    real, dimension(:,:), intent(inout) :: tabSol
    integer :: k
    real :: a1,a2

    a1=1.
    a2=0.
    do k=1,mesh%nc
       if (((mesh%cell(k)%xc-a1*t)**2+(mesh%cell(k)%yc-a2*t)**2)**0.5<5) then
          if((mesh%cell(k)%xc-a1*t>0.).and.(mesh%cell(k)%yc-a2*t>0.)) then
             tabsol(k,1)=1.
          else
             tabsol(k,1)=0.
          endif
       else
          tabsol(k,1)=0.
       endif
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

  subroutine errorL1(sol,exactsol,eL1)
    real, dimension(:), intent(in) :: sol,exactsol
    real, intent(out) :: eL1
    integer :: k

    eL1=0.
    do k=1,size(sol(:))
          eL1=eL1+abs(sol(k)-exactsol(k))
    enddo

    return
  end subroutine errorL1

  subroutine errorL2(sol,exactsol,eL2)
    real, dimension(:), intent(in) :: sol,exactsol
    real, intent(out) :: eL2
    integer :: k

    eL2=0.
    do k=1,size(sol(:))
       eL2=eL2+(sol(k)-exactsol(k))**2
    enddo
    eL2=eL2**0.5

    return
  end subroutine errorL2

end module efficiency
