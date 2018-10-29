program test

  use constant
  use reconstruction
  use types
  use inout
  use efficiency
  use phys
  use FV
  use time

  implicit none

  type(meshStruct) :: mesh
  type(solStruct) :: sol
  real(dp) :: xL,xR,yL,yR,cfl,tf,error,error1,error2
  integer :: nx,ny,nx2,ny2,nvar,fs
  character(len=20) :: namefile,str_equa,str_flux,str_time_scheme
  real(dp), dimension(:), allocatable :: X,U,Xstencil,Ustencil
  real(dp), dimension(3) :: X2
  real(dp), dimension(2) :: U2
  integer :: i,j,k,kpos,order,dir,neigh
  real(dp), dimension(2) :: ul,ur
  real(dp) :: val

  
  call init(xL,xR,yL,yR,nx,ny,nvar,cfl,tf,fs,namefile,mesh,sol,str_equa,str_flux,str_time_scheme,order)
  call buildmesh(xL,xR,yL,yR,nx,ny,mesh)
  call IC(mesh,sol)
  call BC(nx,ny,nvar,mesh)
  call userSol(0.0_dp,mesh,sol)
  call writeSol(mesh,sol,namefile,0)
  call calculation(mesh,sol,str_equa,str_flux,str_time_scheme,order,cfl,tf,fs,namefile)
  call errorL1(mesh,sol%val(:,2),sol%user(:,1),error1)

  nx2=nx*2
  ny2=ny*2
  mesh%nc=nx2*ny2
  mesh%np=(nx2+1)*(ny2+1)

  deallocate(mesh%node,mesh%cell,sol%val,sol%user)
  allocate(mesh%node(mesh%np),mesh%cell(mesh%nc))
  allocate(sol%val(mesh%nc,sol%nvar),sol%user(mesh%nc,sol%nsolUser))

  call buildmesh(xL,xR,yL,yR,nx2,ny2,mesh)
  call IC(mesh,sol)
  call BC(nx2,ny2,nvar,mesh)
  call userSol(0.0_dp,mesh,sol)
  call writeSol(mesh,sol,namefile,0)
  call calculation(mesh,sol,str_equa,str_flux,str_time_scheme,order,cfl,tf,fs,namefile)
  call errorL1(mesh,sol%val(:,2),sol%user(:,1),error2)

  print*,"ordre : ",(log(error1)-log(error2))/(log((xR-xL)/nx)-log((xR-xL)/nx2))

  
!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!

  !do i=1,100
     !do j=1,100
        !sol%val(i+(j-1)*100,1)=((50.5_dp-i)*0.1_dp)**1
        !sol%val(i+(j-1)*100,1)=((i-0.5_dp)*0.1_dp)**10
        !sol%val(i+(j-1)*100,2)=((i-0.5_dp)*0.1_dp)**1
        !if(i<100)then
           !sol%val(i+(j-1)*100,1)=1.0_dp
        !else
           !sol%val(i+(j-1)*100,1)=0.0_dp
        !endif
     !enddo
  !enddo

  !k=1
  !order=2
  !dir=1
  !neigh=51
  !call extractDirection(mesh,sol,k,order,dir,1,X,U,kpos)
  !call buildStencil(X,U,kpos,order,Xstencil,Ustencil)
  !print*,Xstencil
  !print*,Ustencil
  !print*,U(kpos),kpos,U(98),98
  !do i=2,99
     !k=i
     !neigh=i+1
     !dir=3
     !call reconstruct(mesh,sol,k,neigh,order,dir,ul,ur)
     !print*,"k=",k
     !print*,ul(2),ur(2)
     !print*,sol%val(k,2),sol%val(k+1,2)
     !print*,"--------------------------------------------"
  !enddo
  !call reconstruct_boundary(mesh,sol,k,order,dir,ul)

!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!

  !X2(1)=0.0_dp
  !X2(2)=1.0_dp
  !X2(3)=2.0_dp
  !U2(1)=0.5_dp
  !U2(2)=1.5_dp

  !call evaluate(3.0_dp,1,X2,U2,val)
  !print*,val


  contains

  subroutine init_FV(str_equa,str_flux,str_time_scheme,f_equa,flux,time_scheme)
    character(len=20), intent(in) :: str_equa,str_flux,str_time_scheme
    procedure (sub_f), pointer, intent(out) :: f_equa
    procedure (sub_flux), pointer, intent(out) :: flux
    procedure (sub_time), pointer, intent(out) :: time_scheme

    if (trim(str_equa)=='transport') then
       f_equa => f_transport
    elseif (trim(str_equa)=='euler') then
       f_equa => f_euler
    else
       print*,trim(str_equa)," equation not implemented"
       call exit()
    endif
    
    if (trim(str_flux)=='godunov') then
       flux => flux_godunov
    elseif (trim(str_flux)=='HLL') then
       flux => flux_HLL
    else
       print*,trim(str_flux)," flux not implemented"
       call exit()
    endif

    if (trim(str_time_scheme)=='euler_exp') then
       time_scheme => euler_exp
    elseif (trim(str_time_scheme)=='SSPRK2') then
       time_scheme => SSPRK2
    elseif (trim(str_time_scheme)=='SSPRK3') then
       time_scheme => SSPRK3
    else
       print*,trim(str_time_scheme)," time scheme not implemented"
       call exit()
    endif

    return
  end subroutine init_FV

  subroutine calculation(mesh,sol,str_equa,str_flux,str_time_scheme,order,cfl,tf,fs,namefile)
    type(meshStruct), intent(inout) :: mesh
    type(solStruct), intent(inout) :: sol
    character(len=20), intent(in) :: str_equa,str_flux,str_time_scheme
    integer, intent(in) :: order
    real(dp), intent(in) :: cfl,tf
    integer, intent(in) :: fs
    character(len=20),intent(in) :: namefile
    integer :: n
    real(dp) :: t
    procedure (sub_f), pointer :: f_equa
    procedure (sub_flux), pointer :: flux
    procedure (sub_time), pointer :: time_scheme

    call init_FV(str_equa,str_flux,str_time_scheme,f_equa,flux,time_scheme)
    
    t=0.0_dp
    n=1
    do while (t<tf)
       call time_scheme(mesh,sol,f_equa,flux,order,cfl,t)
       if (mod(n,fs)==0) then
          call userSol(t,mesh,sol)
          call writeSol(mesh,sol,namefile,n/fs)
          call print(t,n)
       endif
       n=n+1
    enddo

  end subroutine calculation

end program test
