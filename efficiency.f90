module efficiency

  use constant
  use types
  use phys
  
  implicit none
  
contains
  
  subroutine exactTab(t,mesh,tab,exactSol,gauss_weight,order)
    real(dp), intent(in) :: t
    type (meshStruct), intent(in) :: mesh
    real(dp), dimension(:,:), intent(inout) :: tab
    procedure (sub_exactsol), pointer, intent(in) :: exactSol
    real(dp), dimension(:), intent(in) :: gauss_weight
    integer, intent(in) :: order
    integer :: k,p1,p2
    real(dp) :: s

    do k=1,mesh%nc
       tab(k,1)=0.0_dp
       do p1=1,order
          do p2=1,order
             call exactSol(mesh%cell(k)%X_gauss(p1),mesh%cell(k)%Y_gauss(p2),t,s)
             tab(k,1)=tab(k,1)+s*gauss_weight(p1)*gauss_weight(p2)/4.0_dp
          enddo
       enddo
    enddo
    
    return
  end subroutine exactTab

  subroutine userSol(t,mesh,sol,str_equa,exactSol,gauss_weight,order)
    real(dp), intent(in) :: t
    type(meshStruct), intent(in) :: mesh
    type(solStruct), intent(inout) :: sol
    character(len=20), intent(in) :: str_equa
    procedure (sub_exactsol), pointer, intent(in) :: exactSol
    real(dp), dimension(:), intent(in) :: gauss_weight
    integer, intent(in) :: order
    integer :: i,k
    character(len=20) :: str

    do i=1,sol%nsolUser
       select case (sol%var_user(i))
       case(1)
          sol%name_user(i)="SolAnal"
          call exactTab(t,mesh,sol%user(:,i:i),exactSol,gauss_weight,order)
       case(2)
          sol%name_user(i)="Error"
          call exactTab(t,mesh,sol%user(:,i:i),exactSol,gauss_weight,order)
          sol%user(:,i:i)=abs(sol%user(:,i:i)-sol%val(:,1:1))
       case(3)
          if (trim(str_equa)=="euler") then
             sol%name_user(i)="IntE"
             do k=1,mesh%nc
                call unconserv(sol%val(k,:),str_equa,4,sol%user(k,i))
                sol%user(k,i)=sol%user(k,i)/((gamma-1.0_dp)*sol%val(k,1))
             enddo
          else
             print*,"No internal energy computed for your equation"
             call exit()
          endif
       case(101:109)
          write (str,"(I1)")sol%var_user(i)-100
          sol%name_user(i)="Uncons"//trim(str)
          do k=1,mesh%nc
             call unconserv(sol%val(k,:),str_equa,sol%var_user(i)-100,sol%user(k,i))
          enddo
       case default
          print*,"Wrong number of user_sol"
          call exit()
       end select
    enddo

    return
  end subroutine userSol

  subroutine errorL1(mesh,sol,exacttab,eL1)
    type(meshStruct), intent(in) :: mesh
    real(dp), dimension(:), intent(in) :: sol,exacttab
    real(dp), intent(out) :: eL1
    integer :: k
    real(dp) :: dx,dy

    eL1=0.0_dp
    do k=1,size(sol(:))
       dx=mesh%cell(k)%dx
       dy=mesh%cell(k)%dy
       eL1=eL1+abs(sol(k)-exacttab(k))*dx*dy
    enddo

    print*, "errorL1 = ",eL1

    return
  end subroutine errorL1

  subroutine errorL2(mesh,sol,exacttab,eL2)
    type(meshStruct), intent(in) :: mesh
    real(dp), dimension(:), intent(in) :: sol,exacttab
    real(dp), intent(out) :: eL2
    integer :: k
    real(dp) :: dx,dy

    eL2=0.0_dp
    do k=1,size(sol(:))
       dx=mesh%cell(k)%dx
       dy=mesh%cell(k)%dy
       eL2=eL2+dx*dy*(sol(k)-exacttab(k))**2
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

!!!!!!!!!!!!!!!!!!!!!!!!!!!!!! Exact solutions !!!!!!!!!!!!!!!!!!!!!!!!!!!!!!

  subroutine exactSol_sinus(x,y,t,s)
    real(dp), intent(in) :: x,y,t
    real(dp), intent(out) :: s
    real(dp) :: a1,a2

    a1=1.0_dp
    a2=1.0_dp

    s=cos((x-a1*t-5.0_dp)*pi/5.0_dp)+cos((y-a2*t-5.0_dp)*pi/5.0_dp)

    return
  end subroutine exactSol_sinus

  subroutine exactSol_sinus_dis(x,y,t,s)
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
  end subroutine exactSol_sinus_dis

  subroutine exactSol_vortex(x,y,t,s)
    real(dp), intent(in) :: x,y,t
    real(dp), intent(out) :: s
    real(dp) :: x0,y0,xmin,ymin,beta,r,Temp

    x0=-5.0_dp+mod(int(t+5.0_dp),10)+t-int(t)
    y0=-5.0_dp+mod(int(t+5.0_dp),10)+t-int(t)
    beta=5.0_dp
    xmin=min(abs(x-x0),abs(x-x0-10.0_dp),abs(x-x0+10.0_dp))
    ymin=min(abs(y-y0),abs(y-y0-10.0_dp),abs(y-y0+10.0_dp))
    r=sqrt(xmin**2+ymin**2)
    Temp=1.0_dp-(gamma-1.0_dp)*beta**2*exp(1-r**2)/(8.0_dp*gamma*pi**2)
    s=Temp**(1.0_dp/(gamma-1.0_dp))

    return
  end subroutine exactSol_vortex
    

  subroutine exactSol_none(x,y,t,s)
    real(dp), intent(in) :: x,y,t
    real(dp), intent(out) :: s
    if(.false.)print*,x,y,t,s
    
    print*,"There is no exact solution for this configuration"
    call exit()

    return
  end subroutine exactSol_none

end module efficiency
