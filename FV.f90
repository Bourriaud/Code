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

  subroutine calculation(mesh,sol,dx,dy,dt,tf,fs,namefile)
    type(meshStruct), intent(inout) :: mesh
    type(solStruct), intent(inout) :: sol
    real, intent(in) :: dx,dy,dt,tf
    integer, intent(in) :: fs
    type(solStruct) :: sol2
    character(len=20),intent(in) :: namefile
    integer :: i,j,n
    real, dimension(:), allocatable :: Fp,Fm,Gp,Gm
    real :: t
    procedure (sub), pointer :: flux => null()

    allocate(Fp(sol%nvar),Fm(sol%nvar),Gp(sol%nvar),Gm(sol%nvar))
    allocate(sol2%val(mesh%nx,mesh%ny,sol%nvar))
    flux => flux_godunov
    
    t=0.
    n=1
    do while (t<tf)
       do j=2,mesh%ny-1
          call boundary(flux,sol%val(1,j,:),mesh%boundtype(1,j),mesh%bound(1,j,:),1,Fm)
          call flux(sol%val(1,j,:),sol%val(2,j,:),1,Fp)
          call flux(sol%val(1,j-1,:),sol%val(1,j,:),2,Gm)
          call flux(sol%val(1,j,:),sol%val(1,j+1,:),2,Gp)
          sol2%val(1,j,:)=sol%val(1,j,:)-dt/dx*(Fp-Fm)-dt/dy*(Gp-Gm)
          call boundary(flux,sol%val(mesh%nx,j,:),mesh%boundtype(mesh%nx,j),mesh%bound(mesh%nx,j,:),3,Fp)
          call flux(sol%val(mesh%nx-1,j,:),sol%val(mesh%nx,j,:),1,Fm)
          call flux(sol%val(mesh%nx,j-1,:),sol%val(mesh%nx,j,:),2,Gm)
          call flux(sol%val(mesh%nx,j,:),sol%val(mesh%nx,j+1,:),2,Gp)
          sol2%val(mesh%nx,j,:)=sol%val(mesh%nx,j,:)-dt/dx*(Fp-Fm)-dt/dy*(Gp-Gm)
       enddo
       do i=2,mesh%nx-1
          call boundary(flux,sol%val(i,1,:),mesh%boundtype(i,1),mesh%bound(i,1,:),2,Gm)
          call flux(sol%val(i-1,1,:),sol%val(i,1,:),1,Fm)
          call flux(sol%val(i,1,:),sol%val(i+1,1,:),1,Fp)
          call flux(sol%val(i,1,:),sol%val(i,2,:),2,Gp)
          sol2%val(i,1,:)=sol%val(i,1,:)-dt/dx*(Fp-Fm)-dt/dy*(Gp-Gm)
          call boundary(flux,sol%val(i,mesh%ny,:),mesh%boundtype(i,mesh%ny),mesh%bound(i,mesh%ny,:),4,Gp)
          call flux(sol%val(i-1,mesh%ny,:),sol%val(i,mesh%ny,:),1,Fm)
          call flux(sol%val(i,mesh%ny,:),sol%val(i+1,mesh%ny,:),1,Fp)
          call flux(sol%val(i,mesh%ny-1,:),sol%val(i,mesh%ny,:),2,Gm)
          sol2%val(i,mesh%ny,:)=sol%val(i,mesh%ny,:)-dt/dx*(Fp-Fm)-dt/dy*(Gp-Gm)
       enddo
       call boundary(flux,sol%val(1,1,:),mesh%boundtype(1,1),mesh%bound(1,1,:),1,Fm)
       call boundary(flux,sol%val(1,1,:),mesh%boundtype(1,1),mesh%bound(1,1,:),2,Gm)
       call flux(sol%val(1,1,:),sol%val(2,1,:),1,Fp)
       call flux(sol%val(1,1,:),sol%val(1,2,:),2,Gp)
       sol2%val(1,1,:)=sol%val(1,1,:)-dt/dx*(Fp-Fm)-dt/dy*(Gp-Gm)
       call boundary(flux,sol%val(1,mesh%ny,:),mesh%boundtype(1,mesh%ny),mesh%bound(1,mesh%ny,:),1,Fm)
       call boundary(flux,sol%val(1,mesh%ny,:),mesh%boundtype(1,mesh%ny),mesh%bound(1,mesh%ny,:),4,Gp)
       call flux(sol%val(1,mesh%ny,:),sol%val(2,mesh%ny,:),1,Fp)
       call flux(sol%val(1,mesh%ny-1,:),sol%val(1,mesh%ny,:),2,Gm)
       sol2%val(1,mesh%ny,:)=sol%val(1,mesh%ny,:)-dt/dx*(Fp-Fm)-dt/dy*(Gp-Gm)
       call boundary(flux,sol%val(mesh%nx,1,:),mesh%boundtype(mesh%nx,1),mesh%bound(mesh%nx,1,:),3,Fp)
       call boundary(flux,sol%val(mesh%nx,1,:),mesh%boundtype(mesh%nx,1),mesh%bound(mesh%nx,1,:),2,Gm)
       call flux(sol%val(mesh%nx-1,1,:),sol%val(mesh%nx,1,:),1,Fm)
       call flux(sol%val(mesh%nx,1,:),sol%val(mesh%nx,2,:),2,Gp)
       sol2%val(mesh%nx,1,:)=sol%val(mesh%nx,1,:)-dt/dx*(Fp-Fm)-dt/dy*(Gp-Gm)
       call boundary(flux,sol%val(mesh%nx,mesh%ny,:),mesh%boundtype(mesh%nx,mesh%ny),mesh%bound(mesh%nx,mesh%ny,:),3,Fp)
       call boundary(flux,sol%val(mesh%nx,mesh%ny,:),mesh%boundtype(mesh%nx,mesh%ny),mesh%bound(mesh%nx,mesh%ny,:),4,Gp)
       call flux(sol%val(mesh%nx-1,mesh%ny,:),sol%val(mesh%nx,mesh%ny,:),1,Fm)
       call flux(sol%val(mesh%nx,mesh%ny-1,:),sol%val(mesh%nx,mesh%ny,:),2,Gm)
       sol2%val(mesh%nx,mesh%ny,:)=sol%val(mesh%nx,mesh%ny,:)-dt/dx*(Fp-Fm)-dt/dy*(Gp-Gm)
       do i=2,mesh%nx-1
          do j=2,mesh%ny-1
             call flux(sol%val(i-1,j,:),sol%val(i,j,:),1,Fm)
             call flux(sol%val(i,j,:),sol%val(i+1,j,:),1,Fp)
             call flux(sol%val(i,j-1,:),sol%val(i,j,:),2,Gm)
             call flux(sol%val(i,j,:),sol%val(i,j+1,:),2,Gp)
             sol2%val(i,j,:)=sol%val(i,j,:)-dt/dx*(Fp-Fm)-dt/dy*(Gp-Gm)
          enddo
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

    deallocate(Fp,Fm,Gp,Gm)
    deallocate(sol2%val)
    
  end subroutine calculation

  subroutine boundary(flux_ptr,u,boundtype,bound,side,F)
    procedure (sub), pointer, intent(in) :: flux_ptr
    real, dimension(:), intent(in) :: u
    character(len=20), intent(in) :: boundtype
    real, dimension(:), intent(in) :: bound
    integer, intent(in) :: side      !1=left 2=bottom 3=right 4=top
    real, dimension(:), intent(inout) :: F

    if (trim(boundtype)=='DIRICHLET') then
       if (side==1) then
          call flux_ptr(bound,u,1,F)
       elseif (side==2) then
          call flux_ptr(bound,u,2,F)
       elseif (side==3) then
          call flux_ptr(u,bound,1,F)
       else
          call flux_ptr(u,bound,2,F)
       endif
    else
       print*,"Boundary condition ",trim(boundtype)," not implemented"
       call exit()
    endif

  end subroutine boundary
  
  subroutine flux_godunov(u1,u2,dir,F)
    real, dimension(:), intent(in) :: u1,u2
    integer, intent(in) :: dir
    real, dimension(:), intent(inout) :: F
    real, dimension(:), allocatable :: ustar
    real, dimension(:,:), allocatable :: Fvect

    allocate (ustar(size(u1)),Fvect(size(u1),2))

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
