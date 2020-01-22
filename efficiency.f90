module efficiency

  use constant
  use types
  use phys
  
  implicit none
  
contains
  
  subroutine exactTab(t,mesh,tab,exactSol)
    real(dp), intent(in) :: t
    type (meshStruct), intent(in) :: mesh
    real(dp), dimension(:,:), intent(inout) :: tab
    procedure (sub_exactsol), pointer, intent(in) :: exactSol
    integer :: k,p1,p2
    real(dp) :: s,Xg,Yg

    do k=1,mesh%nc
       tab(k,1)=0.0_dp
       do p1=1,6
          do p2=1,6
             Xg=mesh%cell(k)%xc+gauss_point6(p1)*mesh%cell(k)%dx/2.0_dp
             Yg=mesh%cell(k)%yc+gauss_point6(p2)*mesh%cell(k)%dy/2.0_dp
             call exactSol(Xg,Yg,t,s)
             tab(k,1)=tab(k,1)+s*gauss_weight6(p1)*gauss_weight6(p2)/4.0_dp
          enddo
       enddo
    enddo
    
    return
  end subroutine exactTab

  subroutine exactTab2(mesh,exact_file,dim,tab)
    type (meshStruct), intent(in) :: mesh
    character (len=20), intent(in) :: exact_file
    integer, intent(in) :: dim
    real(dp), dimension(:,:), intent(inout) :: tab
    character (len=20) :: namefile
    integer :: n,i,k,p1,p2
    real(dp), dimension(:), allocatable :: x
    real(dp), dimension(:), allocatable :: solAnal
    real(dp) :: Xg,Yg,dx,c,r

    namefile="sol/"//trim(exact_file)
    open(14,file=namefile,form="formatted")

    read(14,*)n
    allocate(x(n),solAnal(n))
    
    do k=1,n
       read(14,*)x(k),solAnal(k)
    enddo
    dx=(x(n)-x(1))/real(n-1)
    if (dim==2) then
       c=(x(1)+x(n))/2
    endif

    do k=1,mesh%nc
       tab(k,1)=0.0_dp
       do p1=1,6
          do p2=1,6
             Xg=mesh%cell(k)%xc+gauss_point6(p1)*mesh%cell(k)%dx/2.0_dp
             if (dim==1) then
                r=Xg
             else
                Yg=mesh%cell(k)%yc+gauss_point6(p2)*mesh%cell(k)%dy/2.0_dp
                r=sqrt((Xg-c)**2+(Yg-c)**2)
             endif
             r=r-x(1)
             i=ceiling(r/dx)
             if (i>n) i=n
             tab(k,1)=tab(k,1)+solAnal(i)*gauss_weight6(p1)*gauss_weight6(p2)/4.0_dp
          enddo
       enddo
    enddo

    deallocate(x,solAnal)

    close(14)
    
    
  end subroutine exactTab2

  subroutine userSol(t,mesh,sol,str_equa,exactSol,exact_file,dim)
    real(dp), intent(in) :: t
    type(meshStruct), intent(in) :: mesh
    type(solStruct), intent(inout) :: sol
    character(len=20), intent(in) :: str_equa,exact_file
    procedure (sub_exactsol), pointer, intent(in) :: exactSol
    integer, intent(in) :: dim
    integer :: i,k
    character(len=20) :: str

    do i=1,sol%nsolUser
       select case (sol%var_user(i))
       case(1)
          sol%name_user(i)="SolAnal"
          if (trim(exact_file)=="none") then
             call exactTab(t,mesh,sol%user(:,i:i),exactSol)
          else
             call exactTab2(mesh,exact_file,dim,sol%user(:,i:i))
          endif
       case(2)
          sol%name_user(i)="Error"
          if (trim(exact_file)=="none") then
             call exactTab(t,mesh,sol%user(:,i:i),exactSol)
          else
             call exactTab2(mesh,exact_file,dim,sol%user(:,i:i))
          endif
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
       case(4)
          if (trim(str_equa)=="M1") then
             sol%name_user(i)="F_anis"
             do k=1,mesh%nc
                sol%user(k,i)=sol%val(k,2)/(sol%val(k,1)*c)
             enddo
          else
             print*,"No internal energy computed for your equation"
             call exit()
          endif
       case(5)
          if (trim(str_equa)=="M1") then
             sol%name_user(i)="TR"
             do k=1,mesh%nc
                sol%user(k,i)=(sol%val(k,1)/7.56573085e-16_dp)**(0.25_dp)
             enddo
          else
             print*,"No internal energy computed for your equation"
             call exit()
          endif
       case(99)
          sol%name_user(i)="User"
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

  subroutine check_conservativity(mesh,sol,total)
    type(meshStruct), intent(in) :: mesh
    type(solStruct), intent(in) :: sol
    real(dp), dimension(:), intent(inout) :: total
    real(dp) :: dx,dy
    integer :: isol,k
    
    total=0.0_dp
    do k=1,mesh%nc
       dx=mesh%cell(k)%dx
       dy=mesh%cell(k)%dy
       do isol=1,sol%nvar
          total(isol)=total(isol)+sol%val(k,isol)*dx*dy
       enddo
    enddo

    print*,"Total quantities : ",total

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

  subroutine exactSol_test(x,y,t,s)
    real(dp), intent(in) :: x,y,t
    real(dp), intent(out) :: s
    real(dp) :: a1,a2,x_period,scale
    integer :: n
    if(.false.)print*,y

    a1=1.0_dp
    a2=1.0_dp
    scale=10.0_dp

    x_period=x-a1*t-int((x-a1*t)/scale)*scale
    if (x-a1*t<0.0_dp) x_period=x_period+scale
    n=floor(2.0_dp*x_period)
    select case (n)
    case(3)
       s=2.0_dp*(x_period)-3.0_dp
    case(4)
       s=-2.0_dp*(x_period)+5.0_dp
    case(7:8)
       s=1.0_dp
    case(11:12)
       s=2.0_dp*sqrt(0.25_dp-(x_period-6.0_dp)**2)
    case (15)
       s=-2.0_dp*sqrt(0.25_dp-(x_period-7.5_dp)**2)+1.0_dp
    case (16)
       s=-2.0_dp*sqrt(0.25_dp-(x_period-8.5_dp)**2)+1.0_dp
    case default
       s=0.0_dp
    end select

    return
  end subroutine exactSol_test

  subroutine exactSol_test2(x,y,t,s)
    real(dp), intent(in) :: x,y,t
    real(dp), intent(out) :: s
    real(dp) :: a1,a2
    if(.false.)print*,y

    a1=1.0_dp
    a2=1.0_dp

    if (x-a1*t<1.5_dp) then
       s=1.0_dp
    else
       s=0.0_dp
    endif

    return
  end subroutine exactSol_test2

  subroutine exactSol_vortex(x,y,t,s)
    real(dp), intent(in) :: x,y,t
    real(dp), intent(out) :: s
    real(dp) :: x0,y0,xmin,ymin,beta,r,Temp
    if (.false.) print*,t

    x0=0.0_dp!-5.0_dp+mod(int(t+5.0_dp),10)+t-int(t)
    y0=0.0_dp!-5.0_dp+mod(int(t+5.0_dp),10)+t-int(t)
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
    if(.false.)print*,x,y,t

    s=0.0_dp
    print*,"There is no exact solution for this configuration"
    call exit()

    return
  end subroutine exactSol_none

end module efficiency
