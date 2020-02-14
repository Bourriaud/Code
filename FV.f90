module FV

  use constant
  use types
  use phys
  use inout
  use efficiency
  use reconstruction
  
  implicit none
  
contains

  subroutine boundary(flux,f_equa,cell,neigh,u,mesh,sol,edge,p,boundtype,bound,normal,order,F)
    procedure (sub_flux), pointer, intent(in) :: flux
    procedure (sub_f), pointer, intent(in) :: f_equa
    integer, intent(in) :: cell,neigh,edge,p
    real(dp), dimension(:), intent(in) :: u
    type(meshStruct), intent(inout) :: mesh
    type(solStruct), intent(in) :: sol
    character(len=20), intent(in) :: boundtype
    real(dp), dimension(:), intent(in) :: bound
    integer, intent(in) :: normal      !1=left 2=bottom 3=right 4=top
    integer, intent(in) :: order
    real(dp), dimension(:), intent(inout) :: F
    real(dp), dimension(:), allocatable :: v
    integer :: k,period,isol
    
    select case (trim(boundtype))
    case ('DIRICHLET')
       select case (normal)
       case (1)
          call flux(bound,u,f_equa,1,F)
       case (2)
          call flux(bound,u,f_equa,2,F)
       case (3)
          call flux(u,bound,f_equa,1,F)
       case (4)
          call flux(u,bound,f_equa,2,F)
       end select
       
    case ('NEUMANN')
       allocate(v(sol%nvar))
       select case (normal)
       case (1)
          call evaluate(mesh,sol,mesh%cell(cell)%polCoef,cell,order,mesh%cell(cell)%xc,mesh%edge(edge)%Y_gauss(p),v)
          do isol=1,sol%nvar
             v(isol)=sol%val(cell,isol)-bound(isol)*mesh%cell(cell)%dx/2.0_dp
          enddo
          call flux(v,v,f_equa,1,F)
       case (2)
          call evaluate(mesh,sol,mesh%cell(cell)%polCoef,cell,order,mesh%edge(edge)%X_gauss(p),mesh%cell(cell)%yc,v)
          do isol=1,sol%nvar
             v(isol)=sol%val(cell,isol)-bound(isol)*mesh%cell(cell)%dy/2.0_dp
          enddo
          call flux(v,v,f_equa,2,F)
       case (3)
          call evaluate(mesh,sol,mesh%cell(cell)%polCoef,cell,order,mesh%cell(cell)%xc,mesh%edge(edge)%Y_gauss(p),v)
          do isol=1,sol%nvar
             v(isol)=sol%val(cell,isol)+bound(isol)*mesh%cell(cell)%dx/2.0_dp
          enddo
          call flux(v,v,f_equa,1,F)
       case (4)
          call evaluate(mesh,sol,mesh%cell(cell)%polCoef,cell,order,mesh%edge(edge)%X_gauss(p),mesh%cell(cell)%yc,v)
          do isol=1,sol%nvar
             v(isol)=sol%val(cell,isol)+bound(isol)*mesh%cell(cell)%dy/2.0_dp
          enddo
          call flux(v,v,f_equa,2,F)
       end select
       
       deallocate(v)          
       
    case ('TRANSMISSIVE')
       select case (normal)
       case (1)
          call flux(u,u,f_equa,1,F)
       case (2)
          call flux(u,u,f_equa,2,F)
       case (3)
          call flux(u,u,f_equa,1,F)
       case (4)
          call flux(u,u,f_equa,2,F)
       end select
       
    case ('WALL')
       allocate(v(sol%nvar))
       v=u
       select case (normal)
       case (1)
          v(2)=-u(2)
          call flux(v,u,f_equa,1,F)
       case (2)
          v(3)=-u(3)
          call flux(v,u,f_equa,2,F)
       case (3)
          v(2)=-u(2)
          call flux(u,v,f_equa,1,F)
       case (4)
          v(3)=-u(3)
          call flux(u,v,f_equa,2,F)
       end select
       deallocate(v)

    case ('PERIODIC')
       allocate(v(sol%nvar))
       k=abs(neigh)
       period=mesh%edge(edge)%period
       call evaluate(mesh,sol,mesh%cell(k)%polCoef,k,order,mesh%edge(period)%X_gauss(p),mesh%edge(period)%Y_gauss(p),v)
       select case (normal)
       case (1)
          call flux(v,u,f_equa,1,F)
       case (2)
          call flux(v,u,f_equa,2,F)
       case (3)
          call flux(u,v,f_equa,1,F)
       case (4)
          call flux(u,v,f_equa,2,F)
       end select
       deallocate(v)
       
    case default
       print*,"Boundary condition ",trim(boundtype)," not implemented"
       call exit()
    end select
    
  end subroutine boundary

  subroutine RS_advection(u1,u2,f_equa,ustar,dir)
    real(dp), dimension(:), intent(in) :: u1,u2
    procedure (sub_f), pointer, intent(in) :: f_equa
    real(dp), dimension(:), intent(inout) :: ustar
    integer, intent(in) :: dir
    real(dp) :: s
    real(dp), dimension(:,:), allocatable :: F1vect,F2vect
    integer :: i

    allocate(F1vect(size(u1),2),F2vect(size(u1),2))
    
    call f_equa(u1,F1vect(:,:))
    call f_equa(u2,F2vect(:,:))
    
    do i=1,size(u1)
       if (u1(i)/=u2(i)) then
          s=(F2vect(i,dir)-F1vect(i,dir))/(u2(i)-u1(i))
          if (s<0.0_dp) then
             ustar(i)=u2(i)
          else
             ustar(i)=u1(i)
          endif
       else
          ustar(i)=u1(i)
       endif
    enddo

    deallocate(F1vect,F2vect)

    return
  end subroutine RS_advection

  subroutine RS_burgers(u1,u2,f_equa,ustar,dir)
    real(dp), dimension(:), intent(in) :: u1,u2
    procedure (sub_f), pointer, intent(in) :: f_equa
    real(dp), dimension(:), intent(inout) :: ustar
    integer, intent(in) :: dir
    real(dp) :: s
    real(dp), dimension(:,:), allocatable :: F1vect,F2vect
    integer :: i

    allocate(F1vect(size(u1),2),F2vect(size(u1),2))
    
    call f_equa(u1,F1vect(:,:))
    call f_equa(u2,F2vect(:,:))
    
    do i=1,size(u1)
       if (u1(i)<u2(i)) then
          if (u1(i)>0.0_dp) then
             ustar(i)=u1(i)
          elseif (u2(i)<0.0_dp) then
             ustar(i)=u2(i)
          else
             ustar(i)=0.0_dp
          endif
       elseif (u1(i)>u2(i)) then
          s=(F2vect(i,dir)-F1vect(i,dir))/(u2(i)-u1(i))   !(u1(i)+u2(i))/2.0_dp
          if (s<0.0_dp) then
             ustar(i)=u2(i)
          else
             ustar(i)=u1(i)
          endif
       else
          ustar(i)=u1(i)
       endif
    enddo

    deallocate(F1vect,F2vect)

    return
  end subroutine RS_burgers
  
  subroutine flux_godunov_adv(u1,u2,f_equa,dir,F)
    real(dp), dimension(:), intent(in) :: u1,u2
    procedure (sub_f), pointer, intent(in) :: f_equa
    integer, intent(in) :: dir     !1=vertical 2=horizontal
    real(dp), dimension(:), intent(inout) :: F
    real(dp), dimension(:), allocatable :: ustar
    real(dp), dimension(:,:), allocatable :: Fvect

    allocate (ustar(size(u1)),Fvect(size(u1),2))
    
    call RS_advection(u1,u2,f_equa,ustar,dir)
    call f_equa(ustar,Fvect)

    F(:)=Fvect(:,dir)

    deallocate(ustar,Fvect)
    
    return
  end subroutine flux_godunov_adv

  subroutine flux_godunov_bur(u1,u2,f_equa,dir,F)
    real(dp), dimension(:), intent(in) :: u1,u2
    procedure (sub_f), pointer, intent(in) :: f_equa
    integer, intent(in) :: dir     !1=vertical 2=horizontal
    real(dp), dimension(:), intent(inout) :: F
    real(dp), dimension(:), allocatable :: ustar
    real(dp), dimension(:,:), allocatable :: Fvect

    allocate (ustar(size(u1)),Fvect(size(u1),2))
    
    call RS_burgers(u1,u2,f_equa,ustar,dir)
    call f_equa(ustar,Fvect)

    F(:)=Fvect(:,dir)

    deallocate(ustar,Fvect)
    
    return
  end subroutine flux_godunov_bur

  subroutine speed_godunov(u1,u2,f_equa,Smax)
    real(dp), dimension(:), intent(in) :: u1,u2
    procedure (sub_f), pointer, intent(in) :: f_equa
    real(dp), dimension(2), intent(out) :: Smax
    real(dp) :: s
    real(dp), dimension(:,:), allocatable :: f1vect,f2vect
    integer :: i

    allocate(F1vect(size(u1),2),F2vect(size(u1),2))

    Smax=0.0_dp
    
    call f_equa(u1,F1vect(:,:))
    call f_equa(u2,F2vect(:,:))
    
    do i=1,size(u1)
       if (u1(i)/=u2(i)) then
          s=(F2vect(i,1)-F1vect(i,1))/(u2(i)-u1(i))
          Smax(1)=max(abs(s),Smax(1))
          s=(F2vect(i,2)-F1vect(i,2))/(u2(i)-u1(i))
          Smax(2)=max(abs(s),Smax(2))
       endif
    enddo

    deallocate(F1vect,F2vect)

    return
  end subroutine speed_godunov

  subroutine flux_HLL(u1,u2,f_equa,dir,F)
    real(dp), dimension(:), intent(in) :: u1,u2
    procedure (sub_f), pointer, intent(in) :: f_equa
    integer, intent(in) :: dir     !1=vertical 2=horizontal
    real(dp), dimension(:), intent(inout) :: F
    real(dp), dimension(:,:), allocatable :: F1vect,F2vect
    real(dp) :: SL,SR,p1,p2,a1,a2

    allocate (F1vect(size(u1),2),F2vect(size(u2),2))

    if (u1(1)<=0.0_dp) then
       F=0.0_dp
       return
    endif
    if (u2(1)<=0.0_dp) then
       F=0.0_dp
       return
    endif

    call unconserv(u1,"euler               ",4,p1)
    call unconserv(u2,"euler               ",4,p2)

    if (p1<=0.0_dp) then
       F=0.0_dp
       return
    endif
    if (p2<=0.0_dp) then
       F=0.0_dp
       return
    endif
    
    a1=sqrt(gamma*abs(p1/u1(1)))
    a2=sqrt(gamma*abs(p2/u2(1)))
    SL=min(-abs(u1(1+dir)/u1(1))-a1,-abs(u2(1+dir)/u2(1))-a2)
    SR=max(abs(u1(1+dir)/u1(1))+a1,abs(u2(1+dir)/u2(1))+a2)

    if (SL>=0.0_dp) then
       call f_equa(u1,F1vect)
       F(:)=F1vect(:,dir)
    else if (SR<=0.0_dp) then
       call f_equa(u2,F2vect)
       F(:)=F2vect(:,dir)
    else
       call f_equa(u1,F1vect)
       call f_equa(u2,F2vect)
       F(:)=(SR*F1vect(:,dir)-SL*F2vect(:,dir)+SL*SR*(u2(:)-u1(:)))/(SR-SL)
    endif

    deallocate(F1vect,F2vect)

    return
  end subroutine flux_HLL

  subroutine speed_HLL(u1,u2,f_equa,Smax)
    use constant
    real(dp), dimension(:), intent(in) :: u1,u2
    procedure (sub_f), pointer, intent(in) :: f_equa
    real(dp), dimension(2), intent(out) :: Smax
    real(dp) :: SL,SR,p1,p2,a1,a2
    real(dp), dimension(:,:), allocatable :: F
    if(.false.)call f_equa(u1,F)

    call unconserv(u1,"euler               ",4,p1)
    call unconserv(u2,"euler               ",4,p2)
    a1=sqrt(gamma*abs(p1/u1(1)))
    a2=sqrt(gamma*abs(p2/u2(1)))
    
    SL=min(-abs(u1(2)/u1(1))-a1,-abs(u2(2)/u2(1))-a2)
    SR=max(abs(u1(2)/u1(1))+a1,abs(u2(2)/u2(1))+a2)    
    Smax(1)=max(abs(SL),abs(SR))

    SL=min(-abs(u1(3)/u1(1))-a1,-abs(u2(3)/u2(1))-a2)
    SR=max(abs(u1(3)/u1(1))+a1,abs(u2(3)/u2(1))+a2)
    Smax(2)=max(abs(SL),abs(SR))

    return
  end subroutine speed_HLL

  subroutine flux_HLL_is(u1,u2,f_equa,dir,F)
    real(dp), dimension(:), intent(in) :: u1,u2
    procedure (sub_f), pointer, intent(in) :: f_equa
    integer, intent(in) :: dir     !1=vertical 2=horizontal
    real(dp), dimension(:), intent(inout) :: F
    real(dp), dimension(:,:), allocatable :: F1vect,F2vect
    real(dp) :: SL,SR,p1,p2,a1,a2

    allocate (F1vect(size(u1),2),F2vect(size(u2),2))

    if (u1(1)<=0.0_dp) then
       F=0.0_dp
       return
    endif
    if (u2(1)<=0.0_dp) then
       F=0.0_dp
       return
    endif

    call unconserv(u1,"euler_is            ",4,p1)
    call unconserv(u2,"euler_is            ",4,p2)

    if (p1<=0.0_dp) then
       F=0.0_dp
       return
    endif
    if (p2<=0.0_dp) then
       F=0.0_dp
       return
    endif
    
    a1=sqrt(gamma*abs(p1/u1(1)))
    a2=sqrt(gamma*abs(p2/u2(1)))
    SL=min(-abs(u1(1+dir)/u1(1))-a1,-abs(u2(1+dir)/u2(1))-a2)
    SR=max(abs(u1(1+dir)/u1(1))+a1,abs(u2(1+dir)/u2(1))+a2)

    if (SL>=0.0_dp) then
       call f_equa(u1,F1vect)
       F(:)=F1vect(:,dir)
    else if (SR<=0.0_dp) then
       call f_equa(u2,F2vect)
       F(:)=F2vect(:,dir)
    else
       call f_equa(u1,F1vect)
       call f_equa(u2,F2vect)
       F(:)=(SR*F1vect(:,dir)-SL*F2vect(:,dir)+SL*SR*(u2(:)-u1(:)))/(SR-SL)
    endif

    deallocate(F1vect,F2vect)

    return
  end subroutine flux_HLL_is

  subroutine speed_HLL_is(u1,u2,f_equa,Smax)
    use constant
    real(dp), dimension(:), intent(in) :: u1,u2
    procedure (sub_f), pointer, intent(in) :: f_equa
    real(dp), dimension(2), intent(out) :: Smax
    real(dp) :: SL,SR,p1,p2,a1,a2
    real(dp), dimension(:,:), allocatable :: F
    if(.false.)call f_equa(u1,F)

    call unconserv(u1,"euler_is            ",4,p1)
    call unconserv(u2,"euler_is            ",4,p2)
    a1=sqrt(gamma*abs(p1/u1(1)))
    a2=sqrt(gamma*abs(p2/u2(1)))
    
    SL=min(-abs(u1(2)/u1(1))-a1,-abs(u2(2)/u2(1))-a2)
    SR=max(abs(u1(2)/u1(1))+a1,abs(u2(2)/u2(1))+a2)    
    Smax(1)=max(abs(SL),abs(SR))

    SL=min(-abs(u1(3)/u1(1))-a1,-abs(u2(3)/u2(1))-a2)
    SR=max(abs(u1(3)/u1(1))+a1,abs(u2(3)/u2(1))+a2)
    Smax(2)=max(abs(SL),abs(SR))

    return
  end subroutine speed_HLL_is

  subroutine flux_HLLC(u1,u2,f_equa,dir,F)
    real(dp), dimension(:), intent(in) :: u1,u2
    procedure (sub_f), pointer, intent(in) :: f_equa
    integer, intent(in) :: dir     !1=vertical 2=horizontal
    real(dp), dimension(:), intent(inout) :: F
    real(dp) :: aL,aR,aM,rhoM,pL,pR,ppvrs,rhoL,rhoR,uL,uR,vL,vR,ps
    real(dp) :: qL,qR,SL,SR,Ss
    real(dp), dimension(:,:), allocatable :: FL,FR
    real(dp), dimension(:), allocatable :: ULs,URs

    allocate (ULs(size(u1)),URs(size(u2)),FL(size(u1),2),FR(size(u2),2))

    if (u1(1)<=0.0_dp) then
       F=0.0_dp
       return
    endif
    if (u2(1)<=0.0_dp) then
       F=0.0_dp
       return
    endif
    
    call unconserv(u1,"euler               ",1,rhoL)
    call unconserv(u2,"euler               ",1,rhoR)
    call unconserv(u1,"euler               ",1+dir,uL)
    call unconserv(u2,"euler               ",1+dir,uR)
    call unconserv(u1,"euler               ",4-dir,vL)
    call unconserv(u2,"euler               ",4-dir,vR)
    call unconserv(u1,"euler               ",4,pL)
    call unconserv(u2,"euler               ",4,pR)

    if (pL<=0.0_dp) then
       F=0.0_dp
       return
    endif
    if (pR<=0.0_dp) then
       F=0.0_dp
       return
    endif
    
    rhoM=(u1(1)+u2(1))/2.0_dp
    aL=sqrt(gamma*abs(pL/rhoL))
    aR=sqrt(gamma*abs(pR/rhoR))
    aM=(aL+aR)/2.0_dp
    ppvrs=(pL+pR)/2.0_dp-((uR-uL)*rhoM*aM)/2.0_dp
    ps=max(0.0_dp,ppvrs)

    if (ps<=pL) then
       qL=1.0_dp
    else
       qL=sqrt(1.0_dp+((gamma+1.0_dp)/2.0_dp)*(ps/pL-1.0_dp))
    endif

    if (ps<=pR) then
       qR=1.0_dp
    else
       qR=sqrt(1.0_dp+((gamma+1.0_dp)/2.0_dp)*(ps/pR-1.0_dp))
    endif
    SL=uL-aL*qL
    SR=uR+aR*qR
    Ss=(pR-pL+rhoL*uL*(SL-uL)-rhoR*uR*(SR-uR))/(rhoL*(SL-uL)-rhoR*(SR-uR))

    if (SL>=0.0_dp) then
       call f_equa(u1,FL)
       F=FL(:,dir)
    elseif (Ss>=0.0_dp) then
       call f_equa(u1,FL)
       ULs(1)=1.0_dp
       ULs(1+dir)=Ss
       ULs(4-dir)=vL
       ULs(4)=u1(3)/rhoL+(Ss-uL)*(Ss+pL/(rhoL*(SL-uL)))
       ULs=ULs*rhoL*(SL-uL)/(SL-Ss)
       F=FL(:,dir)+SL*(ULs-u1)
    elseif (SR>0.0_dp) then
       call f_equa(u2,FR)
       URs(1)=1.0_dp
       URs(1+dir)=Ss
       URs(4-dir)=vR
       URs(4)=u2(3)/rhoR+(Ss-uR)*(Ss+pR/(rhoR*(SR-uR)))
       URs=URs*rhoR*(SR-uR)/(SR-Ss)
       F=FR(:,dir)+SR*(URs-u2)
    else
       call f_equa(u2,FR)
       F=FR(:,dir)
    endif

    deallocate(ULs,URs,FL,FR)

    return
  end subroutine flux_HLLC

  subroutine speed_HLLC(u1,u2,f_equa,Smax)
    real(dp), dimension(:), intent(in) :: u1,u2
    procedure (sub_f), pointer, intent(in) :: f_equa
    real(dp), dimension(2), intent(out) :: Smax
    real(dp) :: aL,aR,aM,rhoM,pL,pR,ppvrs,rhoL,rhoR,uL,uR,vL,vR,ps
    real(dp) :: qL,qR,SL,SR
    real(dp), dimension(:,:), allocatable :: F
    if(.false.)call f_equa(u1,F)

    call unconserv(u1,"euler               ",1,rhoL)
    call unconserv(u2,"euler               ",1,rhoR)
    call unconserv(u1,"euler               ",2,uL)
    call unconserv(u2,"euler               ",2,uR)
    call unconserv(u1,"euler               ",3,vL)
    call unconserv(u2,"euler               ",3,vR)
    call unconserv(u1,"euler               ",4,pL)
    call unconserv(u2,"euler               ",4,pR)
    
    rhoM=(u1(1)+u2(1))/2.0_dp
    aL=sqrt(gamma*abs(pL/rhoL))
    aR=sqrt(gamma*abs(pR/rhoR))
    aM=(aL+aR)/2.0_dp
    
    ppvrs=(pL+pR)/2.0_dp-((uR-uL)*rhoM*aM)/2.0_dp
    ps=max(0.0_dp,ppvrs)
    if (ps<=pL) then
       qL=1.0_dp
    else
       qL=sqrt(1.0_dp+((gamma+1.0_dp)/2.0_dp)*(ps/pL-1.0_dp))
    endif
    if (ps<=pR) then
       qR=1.0_dp
    else
       qR=sqrt(1.0_dp+((gamma+1.0_dp)/2.0_dp)*(ps/pR-1.0_dp))
    endif
    SL=uL-aL*qL
    SR=uR+aR*qR
    Smax(1)=max(abs(SL),abs(SR))

    ppvrs=(pL+pR)/2.0_dp-((vR-vL)*rhoM*aM)/2.0_dp
    ps=max(0.0_dp,ppvrs)
    if (ps<=pL) then
       qL=1.0_dp
    else
       qL=sqrt(1.0_dp+((gamma+1.0_dp)/2.0_dp)*(ps/pL-1.0_dp))
    endif
    if (ps<=pR) then
       qR=1.0_dp
    else
       qR=sqrt(1.0_dp+((gamma+1.0_dp)/2.0_dp)*(ps/pR-1.0_dp))
    endif
    SL=vL-aL*qL
    SR=vR+aR*qR
    Smax(2)=max(abs(SL),abs(SR))
       
    return
  end subroutine speed_HLLC

  subroutine flux_HLLC_is(u1,u2,f_equa,dir,F)
    real(dp), dimension(:), intent(in) :: u1,u2
    procedure (sub_f), pointer, intent(in) :: f_equa
    integer, intent(in) :: dir     !1=vertical 2=horizontal
    real(dp), dimension(:), intent(inout) :: F
    real(dp) :: aL,aR,aM,rhoM,pL,pR,ppvrs,rhoL,rhoR,uL,uR,vL,vR,ps
    real(dp) :: qL,qR,SL,SR,Ss
    real(dp), dimension(:,:), allocatable :: FL,FR
    real(dp), dimension(:), allocatable :: ULs,URs

    allocate (ULs(size(u1)),URs(size(u2)),FL(size(u1),2),FR(size(u2),2))

    if (u1(1)<=0.0_dp) then
       F=0.0_dp
       return
    endif
    if (u2(1)<=0.0_dp) then
       F=0.0_dp
       return
    endif
    
    call unconserv(u1,"euler_is            ",1,rhoL)
    call unconserv(u2,"euler_is            ",1,rhoR)
    call unconserv(u1,"euler_is            ",1+dir,uL)
    call unconserv(u2,"euler_is            ",1+dir,uR)
    call unconserv(u1,"euler_is            ",4-dir,vL)
    call unconserv(u2,"euler_is            ",4-dir,vR)
    call unconserv(u1,"euler_is            ",4,pL)
    call unconserv(u2,"euler_is            ",4,pR)

    if (pL<=0.0_dp) then
       F=0.0_dp
       return
    endif
    if (pR<=0.0_dp) then
       F=0.0_dp
       return
    endif
    
    rhoM=(u1(1)+u2(1))/2.0_dp
    aL=sqrt(gamma*abs(pL/rhoL))
    aR=sqrt(gamma*abs(pR/rhoR))
    aM=(aL+aR)/2.0_dp
    ppvrs=(pL+pR)/2.0_dp-((uR-uL)*rhoM*aM)/2.0_dp
    ps=max(0.0_dp,ppvrs)

    if (ps<=pL) then
       qL=1.0_dp
    else
       qL=sqrt(1.0_dp+((gamma+1.0_dp)/2.0_dp)*(ps/pL-1.0_dp))
    endif

    if (ps<=pR) then
       qR=1.0_dp
    else
       qR=sqrt(1.0_dp+((gamma+1.0_dp)/2.0_dp)*(ps/pR-1.0_dp))
    endif
    SL=uL-aL*qL
    SR=uR+aR*qR
    Ss=(pR-pL+rhoL*uL*(SL-uL)-rhoR*uR*(SR-uR))/(rhoL*(SL-uL)-rhoR*(SR-uR))

    if (SL>=0.0_dp) then
       call f_equa(u1,FL)
       F=FL(:,dir)
    elseif (Ss>=0.0_dp) then
       call f_equa(u1,FL)
       ULs(1)=1.0_dp
       ULs(1+dir)=Ss
       ULs(4-dir)=vL
       ULs(4)=u1(3)/rhoL+(Ss-uL)*(Ss+pL/(rhoL*(SL-uL)))
       ULs=ULs*rhoL*(SL-uL)/(SL-Ss)
       F=FL(:,dir)+SL*(ULs-u1)
    elseif (SR>0.0_dp) then
       call f_equa(u2,FR)
       URs(1)=1.0_dp
       URs(1+dir)=Ss
       URs(4-dir)=vR
       URs(4)=u2(3)/rhoR+(Ss-uR)*(Ss+pR/(rhoR*(SR-uR)))
       URs=URs*rhoR*(SR-uR)/(SR-Ss)
       F=FR(:,dir)+SR*(URs-u2)
    else
       call f_equa(u2,FR)
       F=FR(:,dir)
    endif

    deallocate(ULs,URs,FL,FR)

    return
  end subroutine flux_HLLC_is

  subroutine speed_HLLC_is(u1,u2,f_equa,Smax)
    real(dp), dimension(:), intent(in) :: u1,u2
    procedure (sub_f), pointer, intent(in) :: f_equa
    real(dp), dimension(2), intent(out) :: Smax
    real(dp) :: aL,aR,aM,rhoM,pL,pR,ppvrs,rhoL,rhoR,uL,uR,vL,vR,ps
    real(dp) :: qL,qR,SL,SR
    real(dp), dimension(:,:), allocatable :: F
    if(.false.)call f_equa(u1,F)

    call unconserv(u1,"euler_is            ",1,rhoL)
    call unconserv(u2,"euler_is            ",1,rhoR)
    call unconserv(u1,"euler_is            ",2,uL)
    call unconserv(u2,"euler_is            ",2,uR)
    call unconserv(u1,"euler_is            ",3,vL)
    call unconserv(u2,"euler_is            ",3,vR)
    call unconserv(u1,"euler_is            ",4,pL)
    call unconserv(u2,"euler_is            ",4,pR)
    
    rhoM=(u1(1)+u2(1))/2.0_dp
    aL=sqrt(gamma*abs(pL/rhoL))
    aR=sqrt(gamma*abs(pR/rhoR))
    aM=(aL+aR)/2.0_dp
    
    ppvrs=(pL+pR)/2.0_dp-((uR-uL)*rhoM*aM)/2.0_dp
    ps=max(0.0_dp,ppvrs)
    if (ps<=pL) then
       qL=1.0_dp
    else
       qL=sqrt(1.0_dp+((gamma+1.0_dp)/2.0_dp)*(ps/pL-1.0_dp))
    endif
    if (ps<=pR) then
       qR=1.0_dp
    else
       qR=sqrt(1.0_dp+((gamma+1.0_dp)/2.0_dp)*(ps/pR-1.0_dp))
    endif
    SL=uL-aL*qL
    SR=uR+aR*qR
    Smax(1)=max(abs(SL),abs(SR))

    ppvrs=(pL+pR)/2.0_dp-((vR-vL)*rhoM*aM)/2.0_dp
    ps=max(0.0_dp,ppvrs)
    if (ps<=pL) then
       qL=1.0_dp
    else
       qL=sqrt(1.0_dp+((gamma+1.0_dp)/2.0_dp)*(ps/pL-1.0_dp))
    endif
    if (ps<=pR) then
       qR=1.0_dp
    else
       qR=sqrt(1.0_dp+((gamma+1.0_dp)/2.0_dp)*(ps/pR-1.0_dp))
    endif
    SL=vL-aL*qL
    SR=vR+aR*qR
    Smax(2)=max(abs(SL),abs(SR))
       
    return
  end subroutine speed_HLLC_is

  subroutine flux_rusanov(u1,u2,f_equa,dir,F)
    real(dp), dimension(:), intent(in) :: u1,u2
    procedure (sub_f), pointer, intent(in) :: f_equa
    integer, intent(in) :: dir     !1=vertical 2=horizontal
    real(dp), dimension(:), intent(inout) :: F
    real(dp) :: aL,aR,aM,rhoM,pL,pR,rhoL,rhoR,uL,uR
    real(dp) :: Sp
    real(dp), dimension(:,:), allocatable :: FL,FR
    real(dp), dimension(:), allocatable :: ULs,URs

    allocate (ULs(size(u1)),URs(size(u2)),FL(size(u1),2),FR(size(u2),2))

    if (u1(1)<=0.0_dp) then
       F=0.0_dp
       return
    endif
    if (u2(1)<=0.0_dp) then
       F=0.0_dp
       return
    endif
    
    call unconserv(u1,"euler               ",1,rhoL)
    call unconserv(u2,"euler               ",1,rhoR)
    call unconserv(u1,"euler               ",1+dir,uL)
    call unconserv(u2,"euler               ",1+dir,uR)
    call unconserv(u1,"euler               ",4,pL)
    call unconserv(u2,"euler               ",4,pR)

    if (pL<=0.0_dp) then
       F=0.0_dp
       return
    endif
    if (pR<=0.0_dp) then
       F=0.0_dp
       return
    endif
    
    rhoM=(u1(1)+u2(1))/2.0_dp
    aL=sqrt(gamma*abs(pL/rhoL))
    aR=sqrt(gamma*abs(pR/rhoR))
    aM=(aL+aR)/2.0_dp
    
    !Sp=max(abs(uL-aL),abs(uR-aR),abs(uL+aL),abs(uR+aR))
    Sp=max(abs(uL)+aL,abs(uR)+aR)

    call f_equa(u1,FL)
    call f_equa(u2,FR)
    
    F=0.5_dp*(FL(:,dir)+FR(:,dir))-0.5_dp*Sp*(u2-u1)

    deallocate(ULs,URs,FL,FR)

    return
  end subroutine flux_rusanov

  subroutine speed_rusanov(u1,u2,f_equa,Smax)
    real(dp), dimension(:), intent(in) :: u1,u2
    procedure (sub_f), pointer, intent(in) :: f_equa
    real(dp), dimension(2), intent(out) :: Smax
    real(dp) :: aL,aR,aM,rhoM,pL,pR,rhoL,rhoR,uL,uR,vL,vR
    real(dp) :: Sp
    real(dp), dimension(4,2) :: test
    if(.false.)call f_equa(u1,test)

    call unconserv(u1,"euler               ",1,rhoL)
    call unconserv(u2,"euler               ",1,rhoR)
    call unconserv(u1,"euler               ",2,uL)
    call unconserv(u2,"euler               ",2,uR)
    call unconserv(u1,"euler               ",3,vL)
    call unconserv(u2,"euler               ",3,vR)
    call unconserv(u1,"euler               ",4,pL)
    call unconserv(u2,"euler               ",4,pR)
    
    rhoM=(u1(1)+u2(1))/2.0_dp
    aL=sqrt(gamma*abs(pL/rhoL))
    aR=sqrt(gamma*abs(pR/rhoR))
    aM=(aL+aR)/2.0_dp
    
    !Sp=max(abs(uL-aL),abs(uR-aR),abs(uL+aL),abs(uR+aR))
    Sp=max(abs(uL)+aL,abs(uR)+aR)
    Smax(1)=Sp

    !Sp=max(abs(uL-aL),abs(uR-aR),abs(uL+aL),abs(uR+aR))
    Sp=max(abs(vL)+aL,abs(vR)+aR)
    Smax(2)=Sp

    return
  end subroutine speed_rusanov

  subroutine flux_rusanov_is(u1,u2,f_equa,dir,F)
    real(dp), dimension(:), intent(in) :: u1,u2
    procedure (sub_f), pointer, intent(in) :: f_equa
    integer, intent(in) :: dir     !1=vertical 2=horizontal
    real(dp), dimension(:), intent(inout) :: F
    real(dp) :: aL,aR,aM,rhoM,pL,pR,rhoL,rhoR,uL,uR
    real(dp) :: Sp
    real(dp), dimension(:,:), allocatable :: FL,FR
    real(dp), dimension(:), allocatable :: ULs,URs

    allocate (ULs(size(u1)),URs(size(u2)),FL(size(u1),2),FR(size(u2),2))

    if (u1(1)<=0.0_dp) then
       F=0.0_dp
       return
    endif
    if (u2(1)<=0.0_dp) then
       F=0.0_dp
       return
    endif
    
    call unconserv(u1,"euler_is            ",1,rhoL)
    call unconserv(u2,"euler_is            ",1,rhoR)
    call unconserv(u1,"euler_is            ",1+dir,uL)
    call unconserv(u2,"euler_is            ",1+dir,uR)
    call unconserv(u1,"euler_is            ",4,pL)
    call unconserv(u2,"euler_is            ",4,pR)

    if (pL<=0.0_dp) then
       F=0.0_dp
       return
    endif
    if (pR<=0.0_dp) then
       F=0.0_dp
       return
    endif
    
    rhoM=(u1(1)+u2(1))/2.0_dp
    aL=sqrt(gamma*abs(pL/rhoL))
    aR=sqrt(gamma*abs(pR/rhoR))
    aM=(aL+aR)/2.0_dp
    
    !Sp=max(abs(uL-aL),abs(uR-aR),abs(uL+aL),abs(uR+aR))
    Sp=max(abs(uL)+aL,abs(uR)+aR)

    call f_equa(u1,FL)
    call f_equa(u2,FR)
    
    F=0.5_dp*(FL(:,dir)+FR(:,dir))-0.5_dp*Sp*(u2-u1)

    deallocate(ULs,URs,FL,FR)

    return
  end subroutine flux_rusanov_is

  subroutine speed_rusanov_is(u1,u2,f_equa,Smax)
    real(dp), dimension(:), intent(in) :: u1,u2
    procedure (sub_f), pointer, intent(in) :: f_equa
    real(dp), dimension(2), intent(out) :: Smax
    real(dp) :: aL,aR,aM,rhoM,pL,pR,rhoL,rhoR,uL,uR,vL,vR
    real(dp) :: Sp
    real(dp), dimension(3,2) :: test
    if(.false.)call f_equa(u1,test)

    call unconserv(u1,"euler_is            ",1,rhoL)
    call unconserv(u2,"euler_is            ",1,rhoR)
    call unconserv(u1,"euler_is            ",2,uL)
    call unconserv(u2,"euler_is            ",2,uR)
    call unconserv(u1,"euler_is            ",3,vL)
    call unconserv(u2,"euler_is            ",3,vR)
    call unconserv(u1,"euler_is            ",4,pL)
    call unconserv(u2,"euler_is            ",4,pR)
    
    rhoM=(u1(1)+u2(1))/2.0_dp
    aL=sqrt(gamma*abs(pL/rhoL))
    aR=sqrt(gamma*abs(pR/rhoR))
    aM=(aL+aR)/2.0_dp
    
    !Sp=max(abs(uL-aL),abs(uR-aR),abs(uL+aL),abs(uR+aR))
    Sp=max(abs(uL)+aL,abs(uR)+aR)
    Smax(1)=Sp

    !Sp=max(abs(uL-aL),abs(uR-aR),abs(uL+aL),abs(uR+aR))
    Sp=max(abs(vL)+aL,abs(vR)+aR)
    Smax(2)=Sp

    return
  end subroutine speed_rusanov_is

  subroutine flux_rusanov_M1(u1,u2,f_equa,dir,F)
    real(dp), dimension(:), intent(in) :: u1,u2
    procedure (sub_f), pointer, intent(in) :: f_equa
    integer, intent(in) :: dir
    real(dp), dimension(:), intent(inout) :: F
    real(dp) :: Sp
    real(dp), dimension(:,:), allocatable :: FL,FR

    allocate (FL(size(u1),2),FR(size(u2),2))   
    Sp=c

    call f_equa(u1,FL)
    call f_equa(u2,FR)
    
    F=0.5_dp*(FL(:,dir)+FR(:,dir))-0.5_dp*Sp*(u2-u1)

    deallocate(FL,FR)

    return
  end subroutine flux_rusanov_M1

  subroutine speed_rusanov_M1(u1,u2,f_equa,Smax)
    real(dp), dimension(:), intent(in) :: u1,u2
    procedure (sub_f), pointer, intent(in) :: f_equa
    real(dp), dimension(2), intent(out) :: Smax
    real(dp), dimension(4,2) :: test
    if(.false.)call f_equa(u1,test)
    if(.false.)call f_equa(u2,test)

    Smax(1)=c

    return
  end subroutine speed_rusanov_M1

end module FV
