module FV

  use types
  use phys
  use inout
  
  implicit none

contains

  subroutine calculation(mesh,sol,dx,dy,dt,tf,namefile)
    type(meshStruct), intent(inout) :: mesh
    type(solStruct), intent(inout) :: sol
    real, intent(in) :: dx,dy,dt,tf
    type(solStruct) :: sol2
    character(len=20),intent(in) :: namefile
    integer :: i,j,n
    real, dimension(:), allocatable :: Fp,Fm,Gp,Gm
    real :: t

    allocate(Fp(sol%nvar),Fm(sol%nvar),Gp(sol%nvar),Gm(sol%nvar))
    allocate(sol2%val(mesh%nx,mesh%ny,sol%nvar))

    t=0.
    n=1
    do while (t<tf)
       do j=2,mesh%ny-1
          call boundary(mesh%boundtype(1,j),mesh%bound(1,j),Fm)
          call flux_godunov(sol%val(1,j,:),sol%val(2,j,:),1,Fp)
          call flux_godunov(sol%val(1,j-1,:),sol%val(1,j,:),2,Gm)
          call flux_godunov(sol%val(1,j,:),sol%val(1,j+1,:),2,Gp)
          sol2%val(1,j,:)=sol%val(1,j,:)-dt/dx*(Fp-Fm)-dt/dy*(Gp-Gm)
          call boundary(mesh%boundtype(mesh%nx,j),mesh%bound(mesh%nx,j),Fp)
          call flux_godunov(sol%val(mesh%nx-1,j,:),sol%val(mesh%nx,j,:),1,Fm)
          call flux_godunov(sol%val(mesh%nx,j-1,:),sol%val(mesh%nx,j,:),2,Gm)
          call flux_godunov(sol%val(mesh%nx,j,:),sol%val(mesh%nx,j+1,:),2,Gp)
          sol2%val(mesh%nx,j,:)=sol%val(mesh%nx,j,:)-dt/dx*(Fp-Fm)-dt/dy*(Gp-Gm)
       enddo
       do i=2,mesh%nx-1
          call boundary(mesh%boundtype(i,1),mesh%bound(i,1),Gm)
          call flux_godunov(sol%val(i-1,1,:),sol%val(i,1,:),1,Fm)
          call flux_godunov(sol%val(i,1,:),sol%val(i+1,1,:),1,Fp)
          call flux_godunov(sol%val(i,1,:),sol%val(i,2,:),2,Gp)
          sol2%val(i,1,:)=sol%val(i,1,:)-dt/dx*(Fp-Fm)-dt/dy*(Gp-Gm)
          call boundary(mesh%boundtype(i,mesh%ny),mesh%bound(i,mesh%ny),Gp)
          call flux_godunov(sol%val(i-1,mesh%ny,:),sol%val(i,mesh%ny,:),1,Fm)
          call flux_godunov(sol%val(i,mesh%ny,:),sol%val(i+1,mesh%ny,:),1,Fp)
          call flux_godunov(sol%val(i,mesh%ny-1,:),sol%val(i,mesh%ny,:),2,Gm)
          sol2%val(i,mesh%ny,:)=sol%val(i,mesh%ny,:)-dt/dx*(Fp-Fm)-dt/dy*(Gp-Gm)
       enddo
       call boundary(mesh%boundtype(1,1),mesh%bound(1,1),Fm)
       call boundary(mesh%boundtype(1,1),mesh%bound(1,1),Gm)
       call flux_godunov(sol%val(1,1,:),sol%val(2,1,:),1,Fp)
       call flux_godunov(sol%val(1,1,:),sol%val(1,2,:),2,Gp)
       sol2%val(1,1,:)=sol%val(1,1,:)-dt/dx*(Fp-Fm)-dt/dy*(Gp-Gm)
       call boundary(mesh%boundtype(1,mesh%ny),mesh%bound(1,mesh%ny),Fm)
       call boundary(mesh%boundtype(1,mesh%ny),mesh%bound(1,mesh%ny),Gp)
       call flux_godunov(sol%val(1,mesh%ny,:),sol%val(2,mesh%ny,:),1,Fp)
       call flux_godunov(sol%val(1,mesh%ny-1,:),sol%val(1,mesh%ny,:),2,Gm)
       sol2%val(1,mesh%ny,:)=sol%val(1,mesh%ny,:)-dt/dx*(Fp-Fm)-dt/dy*(Gp-Gm)
       call boundary(mesh%boundtype(mesh%nx,1),mesh%bound(mesh%nx,1),Fp)
       call boundary(mesh%boundtype(mesh%nx,1),mesh%bound(mesh%nx,1),Gm)
       call flux_godunov(sol%val(mesh%nx-1,1,:),sol%val(mesh%nx,1,:),1,Fm)
       call flux_godunov(sol%val(mesh%nx,1,:),sol%val(mesh%nx,2,:),2,Gp)
       sol2%val(mesh%nx,1,:)=sol%val(mesh%nx,1,:)-dt/dx*(Fp-Fm)-dt/dy*(Gp-Gm)
       call boundary(mesh%boundtype(mesh%nx,mesh%ny),mesh%bound(mesh%nx,mesh%ny),Fp)
       call boundary(mesh%boundtype(mesh%nx,mesh%ny),mesh%bound(mesh%nx,mesh%ny),Gp)
       call flux_godunov(sol%val(mesh%nx-1,mesh%ny,:),sol%val(mesh%nx,mesh%ny,:),1,Fm)
       call flux_godunov(sol%val(mesh%nx,mesh%ny-1,:),sol%val(mesh%nx,mesh%ny,:),2,Gm)
       sol2%val(mesh%nx,mesh%ny,:)=sol%val(mesh%nx,mesh%ny,:)-dt/dx*(Fp-Fm)-dt/dy*(Gp-Gm)
       do i=2,mesh%nx-1
          do j=2,mesh%ny-1
             call flux_godunov(sol%val(i-1,j,:),sol%val(i,j,:),1,Fm)
             call flux_godunov(sol%val(i,j,:),sol%val(i+1,j,:),1,Fp)
             call flux_godunov(sol%val(i,j-1,:),sol%val(i,j,:),2,Gm)
             call flux_godunov(sol%val(i,j,:),sol%val(i,j+1,:),2,Gp)
             sol2%val(i,j,:)=sol%val(i,j,:)-dt/dx*(Fp-Fm)-dt/dy*(Gp-Gm)
          enddo
       enddo
       sol%val=sol2%val
       call writeSol(mesh,sol,namefile,n)
       t=t+dt
       n=n+1
    enddo

    deallocate(Fp,Fm,Gp,Gm)
    deallocate(sol2%val)
    
  end subroutine calculation

  subroutine boundary(boundtype,bound,F)
    character(len=20), intent(in) :: boundtype
    real, intent(in) :: bound
    real, dimension(:), intent(inout) :: F

    F=0.

  end subroutine boundary
  
  subroutine flux_godunov(u1,u2,dir,F)
    real, dimension(:), intent(in) :: u1,u2
    real, dimension(:), intent(inout) :: F
    integer, intent(in) :: dir
    real, dimension(:), allocatable :: ustar
    real, dimension(:,:), allocatable :: Fvect

    allocate (ustar(size(u1)),Fvect(size(u1),2))

    call RP1D(u1,u2,ustar,dir)
    call f_transport(ustar,Fvect)
    F(:)=Fvect(:,dir)

    deallocate(ustar,Fvect)

  end subroutine flux_godunov

  subroutine RP1D(u1,u2,ustar,dir)
    real, dimension(:), intent(in) :: u1,u2
    real, dimension(:), intent(inout) :: ustar
    integer, intent(in) :: dir
    real :: s
    real, dimension(:,:), allocatable :: f1vect,f2vect
    real, dimension(:), allocatable :: f1,f2

    ! Seulement advection scalaire linéaire

    allocate(f1vect(size(u1),2),f2vect(size(u1),2),f1(size(u1)),f2(size(u1)))
    call f_transport(u1,f1vect(:,:))
    call f_transport(u2,f2vect(:,:))
    f1=f1vect(:,1)
    f2=f2vect(:,2)
    
    s=(f2(1)-f1(1))/(u2(1)-u1(1))

    if (s<0.) then
       ustar=u2
    else
       ustar=u1
    endif
  end subroutine RP1D
       

end module FV
