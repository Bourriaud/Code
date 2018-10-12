module FV

  use types
  use phys
  use inout
  use efficiency
  
  implicit none

  abstract interface
     subroutine sub_flux (u1,u2,f_equa,dir,F)
       real, dimension(:), intent(in) :: u1,u2
       procedure (sub_f), pointer, intent(in) :: f_equa
       integer, intent(in) :: dir
       real, dimension(:), intent(inout) :: F
     end subroutine sub_flux
  end interface
  
contains

  subroutine calculation(mesh,sol,dt,tf,fs,namefile)
    type(meshStruct), intent(inout) :: mesh
    type(solStruct), intent(inout) :: sol
    real, intent(in) :: dt,tf
    integer, intent(in) :: fs
    type(solStruct) :: sol2
    character(len=20),intent(in) :: namefile
    integer :: k,i,n,neigh,k1,k2,normal
    real, dimension(:,:), allocatable :: F
    real, dimension(:), allocatable :: Ftemp
    real :: t,dx,dy
    procedure (sub_flux), pointer :: flux => null()
    procedure (sub_f), pointer :: f_equa => null()

    allocate(F(sol%nvar,4),Ftemp(sol%nvar))
    allocate(sol2%val(mesh%nc,sol%nvar))
    flux => flux_HLL
    f_equa => f_euler
    
    t=0.
    n=1
    do while (t<tf)
       do k=1,mesh%nc
          F=0.
          dx=mesh%cell(k)%dx
          dy=mesh%cell(k)%dy
          do i=1,size(mesh%cell(k)%edge)
             neigh=mesh%cell(k)%edge(i)%neigh
             normal=mesh%cell(k)%edge(i)%normal
             if (neigh==-1) then
                call boundary(flux,f_equa,sol%val(k,:),mesh%cell(k)%edge(i)%boundType,mesh%cell(k)%edge(i)%bound,normal,Ftemp(:))
             else
                k1=min(k,neigh)
                k2=max(k,neigh)
                call flux(sol%val(k1,:),sol%val(k2,:),f_equa,normal,Ftemp(:))
             endif
             F(:,normal)=F(:,normal)+Ftemp(:)
          enddo
          sol2%val(k,:)=sol%val(k,:)-dt/dx*(F(:,3)-F(:,1))-dt/dy*(F(:,4)-F(:,2))
       enddo

       sol%val=sol2%val
       
       if (mod(n,fs)==0) then
          call userSol(t,mesh,sol)
          call writeSol(mesh,sol,namefile,n/fs)
          call print(t,n)
       endif
       t=t+dt
       n=n+1
    enddo

    deallocate(F,Ftemp)
    deallocate(sol2%val)

  end subroutine calculation

  subroutine boundary(flux_ptr,f_equa,u,boundtype,bound,normal,F)
    procedure (sub_flux), pointer, intent(in) :: flux_ptr
    procedure (sub_f), pointer, intent(in) :: f_equa
    real, dimension(:), intent(in) :: u
    character(len=20), intent(in) :: boundtype
    real, dimension(:), intent(in) :: bound
    integer, intent(in) :: normal      !1=left 2=bottom 3=right 4=top
    real, dimension(:), intent(inout) :: F
    real, dimension(:), allocatable :: v

    if (trim(boundtype)=='DIRICHLET') then
       if (normal==1) then
          call flux_ptr(bound,u,f_equa,1,F)
       elseif (normal==2) then
          call flux_ptr(bound,u,f_equa,2,F)
       elseif (normal==3) then
          call flux_ptr(u,bound,f_equa,1,F)
       else
          call flux_ptr(u,bound,f_equa,2,F)
       endif
       
    elseif (trim(boundtype)=='TRANSMISSIVE') then
       call flux_ptr(u,u,f_equa,normal,F)
       
    elseif (trim(boundtype)=='WALL') then
       allocate(v(size(u)))
       v=u
       if (normal==1) then
          v(2)=-u(2)
          call flux_ptr(v,u,f_equa,1,F)
       elseif (normal==2) then
          v(3)=-u(3)
          call flux_ptr(v,u,f_equa,2,F)
       elseif (normal==3) then
          v(2)=-u(2)
          call flux_ptr(u,v,f_equa,1,F)
       else
          v(3)=-u(3)
          call flux_ptr(u,v,f_equa,2,F)
       endif
       deallocate(v)
       
    else
       print*,"Boundary condition ",trim(boundtype)," not implemented"
       call exit()
    endif

  end subroutine boundary
  
  subroutine flux_godunov(u1,u2,f_equa,normal,F)
    real, dimension(:), intent(in) :: u1,u2
    procedure (sub_f), pointer, intent(in) :: f_equa
    integer, intent(in) :: normal     !1=left 2=bottom 3=right 4=top
    real, dimension(:), intent(inout) :: F
    real, dimension(:), allocatable :: ustar
    real, dimension(:,:), allocatable :: Fvect
    integer :: dir

    allocate (ustar(size(u1)),Fvect(size(u1),2))

    if ((normal==1).or.(normal)==3) then
       dir=1
    else
       dir=2
    endif
    
    call RS_advection(u1,u2,f_equa,ustar,dir)
    call f_equa(ustar,Fvect)
    
    F(:)=Fvect(:,dir)

    deallocate(ustar,Fvect)
    
    return
  end subroutine flux_godunov

  subroutine RS_advection(u1,u2,f_equa,ustar,dir)
    real, dimension(:), intent(in) :: u1,u2
    procedure (sub_f), pointer, intent(in) :: f_equa
    real, dimension(:), intent(inout) :: ustar
    integer, intent(in) :: dir
    real :: s
    real, dimension(:,:), allocatable :: f1vect,f2vect
    real, dimension(:), allocatable :: f1,f2
    integer :: i

    allocate(f1vect(size(u1),2),f2vect(size(u1),2),f1(size(u1)),f2(size(u1)))
    
    call f_equa(u1,f1vect(:,:))
    call f_equa(u2,f2vect(:,:))
    f1=f1vect(:,dir)
    f2=f2vect(:,dir)

    do i=1,size(u1)
       s=(f2(i)-f1(i))/(u2(i)-u1(i))
       if (s<0.) then
          ustar(i)=u2(i)
       else
          ustar(i)=u1(i)
       endif
    enddo

    return
  end subroutine RS_advection

  subroutine flux_HLL(u1,u2,f_equa,normal,F)
    real, dimension(:), intent(in) :: u1,u2
    procedure (sub_f), pointer, intent(in) :: f_equa
    integer, intent(in) :: normal     !1=left 2=bottom 3=right 4=top
    real, dimension(:), intent(inout) :: F
    real, dimension(:,:), allocatable :: F1vect,F2vect
    integer :: dir
    real :: SL,SR,p1,p2,a1,a2,gamma

    allocate (F1vect(size(u1),2),F2vect(size(u2),2))

    if ((normal==1).or.(normal)==3) then
       dir=1
    else
       dir=2
    endif

    gamma=1.4
    p1=(u1(4)-u1(1)*(0.5*((u1(2)/u1(1))**2+(u1(3)/u1(1))**2)))*(gamma-1)
    p2=(u2(4)-u2(1)*(0.5*((u2(2)/u2(1))**2+(u2(3)/u2(1))**2)))*(gamma-1)
    a1=(gamma*p1/u1(1))**0.5
    a2=(gamma*p2/u2(1))**0.5
    SL=min(u1(2+dir-1)-a1,u2(2+dir-1)-a2)
    SR=max(u1(2+dir-1)+a1,u2(2+dir-1)+a2)
    
    if (SL>0.) then
       call f_equa(u1,F1vect)
       F(:)=F1vect(:,dir)
    else if (SR<0.) then
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
       

end module FV
