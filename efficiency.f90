module efficiency

  use constant
  use types
  
  implicit none
  
contains

  subroutine exactSol(x,y,t,s)
    real(dp), intent(in) :: x,y,t
    real(dp), intent(out) :: s
    real(dp) :: a1,a2

    a1=1.0_dp
    a2=1.0_dp

    s=cos((x-a1*t-5.0_dp)*pi/5.0_dp)+cos((y-a2*t-5.0_dp)*pi/5.0_dp)
    if(s>1.5_dp)then
       s=1.0_dp
    else
       s=0.0_dp
    endif

    return
  end subroutine exactSol
  
  subroutine exactTab(t,mesh,tab,gauss_weight)
    real(dp), intent(in) :: t
    type (meshStruct), intent(in) :: mesh
    real(dp), dimension(:,:), intent(inout) :: tab
    real(dp), dimension(:), intent(in) :: gauss_weight
    integer :: k,p1,p2
    real(dp) :: s

    do k=1,mesh%nc
       tab(k,1)=0.0_dp
       do p1=1,size(mesh%cell(k)%X_gauss)
          do p2=1,size(mesh%cell(k)%Y_gauss)
             call exactSol(mesh%cell(k)%X_gauss(p1),mesh%cell(k)%Y_gauss(p2),t,s)
             tab(k,1)=tab(k,1)+s*gauss_weight(p1)*gauss_weight(p2)/4.0_dp
          enddo
       enddo
    enddo
    
    return
  end subroutine exactTab

  subroutine userSol(t,mesh,sol,gauss_weight)
    real(dp), intent(in) :: t
    type(meshStruct), intent(in) :: mesh
    type(solStruct), intent(inout) :: sol
    real(dp), dimension(:), intent(in) :: gauss_weight

    call exactTab(t,mesh,sol%user(:,1:1),gauss_weight)
    sol%user(:,2:2)=abs(sol%user(:,1:1)-sol%val(:,1:1))

    return
  end subroutine userSol

  subroutine errorL1(mesh,sol,exactsol,eL1)
    type(meshStruct), intent(in) :: mesh
    real(dp), dimension(:), intent(in) :: sol,exactsol
    real(dp), intent(out) :: eL1
    integer :: k
    real(dp) :: dx,dy

    eL1=0.0_dp
    do k=1,size(sol(:))
       dx=mesh%cell(k)%dx
       dy=mesh%cell(k)%dy
       eL1=eL1+abs(sol(k)-exactsol(k))*dx*dy
    enddo

    print*, "errorL1 = ",eL1

    return
  end subroutine errorL1

  subroutine errorL2(mesh,sol,exactsol,eL2)
    type(meshStruct), intent(in) :: mesh
    real(dp), dimension(:), intent(in) :: sol,exactsol
    real(dp), intent(out) :: eL2
    integer :: k
    real(dp) :: dx,dy

    eL2=0.0_dp
    do k=1,size(sol(:))
       dx=mesh%cell(k)%dx
       dy=mesh%cell(k)%dy
       eL2=eL2+dx*dy*(sol(k)-exactsol(k))**2
    enddo
    eL2=sqrt(eL2)

    print*, "errorL2 = ",eL2

    return
  end subroutine errorL2

  subroutine check_conservativity(mesh,sol)
    type(meshStruct), intent(in) :: mesh
    type(solStruct), intent(in) :: sol
    real(dp) :: dx,dy
    real(dp), dimension(:), allocatable :: total
    integer :: isol,k

    allocate(total(sol%nvar))
    
    total=0.0_dp
    do k=1,mesh%nc
       dx=mesh%cell(k)%dx
       dy=mesh%cell(k)%dy
       do isol=1,sol%nvar
          total(isol)=total(isol)+sol%val(k,isol)*dx*dy
       enddo
    enddo

    print*,"Total quantities : ",total

    deallocate(total)

    return
  end subroutine check_conservativity

end module efficiency
