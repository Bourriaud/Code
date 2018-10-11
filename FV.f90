module FV

  use types
  use phys
  use inout
  use efficiency
  
  implicit none

  abstract interface
     subroutine sub (u1,u2,dir,F)
       real, dimension(:), intent(in) :: u1,u2
       integer, intent(in) :: dir
       real, dimension(:), intent(inout) :: F
     end subroutine sub
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
    procedure (sub), pointer :: flux => null()

    allocate(F(sol%nvar,4),Ftemp(sol%nvar))
    allocate(sol2%val(mesh%nc,sol%nvar))
    flux => flux_godunov
    
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
                call boundary(flux,sol%val(k,:),mesh%cell(k)%edge(i)%boundType,mesh%cell(k)%edge(i)%bound,normal,Ftemp(:))
             else
                k1=min(k,neigh)
                k2=max(k,neigh)
                call flux(sol%val(k1,:),sol%val(k2,:),normal,Ftemp(:))
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

  subroutine boundary(flux_ptr,u,boundtype,bound,normal,F)
    procedure (sub), pointer, intent(in) :: flux_ptr
    real, dimension(:), intent(in) :: u
    character(len=20), intent(in) :: boundtype
    real, dimension(:), intent(in) :: bound
    integer, intent(in) :: normal      !1=left 2=bottom 3=right 4=top
    real, dimension(:), intent(inout) :: F

    if (trim(boundtype)=='DIRICHLET') then
       if (normal==1) then
          call flux_ptr(bound,u,1,F)
       elseif (normal==2) then
          call flux_ptr(bound,u,2,F)
       elseif (normal==3) then
          call flux_ptr(u,bound,1,F)
       else
          call flux_ptr(u,bound,2,F)
       endif
    else
       print*,"Boundary condition ",trim(boundtype)," not implemented"
       call exit()
    endif

  end subroutine boundary
  
  subroutine flux_godunov(u1,u2,normal,F)
    real, dimension(:), intent(in) :: u1,u2
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
    call RP(u1,u2,ustar,dir)
    call f_transport(ustar,Fvect)
    F(:)=Fvect(:,dir)

    deallocate(ustar,Fvect)
    
    return
  end subroutine flux_godunov

  subroutine RP(u1,u2,ustar,dir)
    real, dimension(:), intent(in) :: u1,u2
    real, dimension(:), intent(inout) :: ustar
    integer, intent(in) :: dir
    real :: s
    real, dimension(:,:), allocatable :: f1vect,f2vect
    real, dimension(:), allocatable :: f1,f2
    integer :: i

    ! Seulement advection linéaire

    allocate(f1vect(size(u1),2),f2vect(size(u1),2),f1(size(u1)),f2(size(u1)))
    call f_transport(u1,f1vect(:,:))
    call f_transport(u2,f2vect(:,:))
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
  end subroutine RP
       

end module FV
