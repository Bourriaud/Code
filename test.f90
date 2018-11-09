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
  real(dp) :: xL,xR,yL,yR,cfl,tf,error,error1,error2,dx,dy
  integer :: nx,ny,nx2,ny2,nvar,fs
  character(len=20) :: namefile,str_equa,str_flux,str_time_scheme
  real(dp), dimension(:), allocatable :: X,U,Xstencil,Ustencil
  real(dp), dimension(3) :: X2
  real(dp), dimension(2) :: U2
  integer :: i,j,k,kpos,order,dir,neigh
  real(dp), dimension(2) :: ul,ur
  real(dp) :: val
  real(dp), dimension(3,2) :: A
  real(dp), dimension(3) :: b
  real(dp), dimension(2) :: vect
  real(dp), dimension(:,:), allocatable :: R
  procedure (sub_reconstruction), pointer :: func
  procedure (sub_f), pointer :: f_equa
  procedure (sub_flux), pointer :: flux
  procedure (sub_speed), pointer :: speed
  procedure (sub_time), pointer :: time_scheme
  procedure (quadrature_t), pointer :: quad_t
  procedure (quadrature_c_alpha), pointer :: quad_c_alpha
  procedure (quadrature_reconstruction), pointer :: quad_reconstruct

  
  call init(xL,xR,yL,yR,nx,ny,nvar,cfl,tf,fs,namefile,mesh,sol,str_equa,str_flux,str_time_scheme,order)
  call init_FV(str_equa,str_flux,str_time_scheme,order,f_equa,flux,speed,time_scheme,quad_t,quad_c_alpha,quad_reconstruct)
  call buildmesh(xL,xR,yL,yR,nx,ny,mesh)
  call IC(mesh,sol,quad_t)
  call BC(nx,ny,nvar,mesh)
  call userSol(0.0_dp,mesh,sol,quad_t)
  call writeSol(mesh,sol,namefile,0)
  call calculation(mesh,sol,order,cfl,tf,fs,namefile,f_equa,flux,speed,time_scheme,quad_t,quad_c_alpha,quad_reconstruct)
  
  call errorL1(mesh,sol%val(:,2),sol%user(:,1),error1)
  print*, "errorL1 = ",error1

  nx2=nx*2
  ny2=ny*2
  mesh%nc=nx2*ny2
  mesh%ne=nx2*(nx2+1)+ny2*(ny2+1)
  mesh%np=(nx2+1)*(ny2+1)

  deallocate(mesh%node,mesh%edge,mesh%cell,sol%val,sol%user)
  allocate(mesh%node(mesh%np),mesh%edge(mesh%ne),mesh%cell(mesh%nc))
  allocate(sol%val(mesh%nc,sol%nvar),sol%user(mesh%nc,sol%nsolUser))

  call init_FV(str_equa,str_flux,str_time_scheme,order,f_equa,flux,speed,time_scheme,quad_t,quad_c_alpha,quad_reconstruct)
  call buildmesh(xL,xR,yL,yR,nx2,ny2,mesh)
  call IC(mesh,sol,quad_t)
  call BC(nx2,ny2,nvar,mesh)
  call userSol(0.0_dp,mesh,sol,quad_t)
  call writeSol(mesh,sol,namefile,0)
  call calculation(mesh,sol,order,cfl,tf,fs,namefile,f_equa,flux,speed,time_scheme,quad_t,quad_c_alpha,quad_reconstruct)
  
  call errorL1(mesh,sol%val(:,2),sol%user(:,1),error2)
  print*, "errorL1 = ",error2

  print*,"ordre : ",(log(error1)-log(error2))/(log((xR-xL)/nx)-log((xR-xL)/nx2))

  deallocate(mesh%node,mesh%edge,mesh%cell,sol%val,sol%user,sol%name,sol%nameUser)
  
!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!

  !do i=1,10
     !do j=1,10
        !sol%val(i+(j-1)*100,1)=((50.5_dp-i)*0.1_dp)**1
        !sol%val(i+(j-1)*100,1)=((i-0.5_dp)*0.1_dp)**2
        !sol%val(i+(j-1)*100,2)=((j-0.5_dp)*0.1_dp)**2
        !sol%val(i+(j-1)*10,1)=(i**3/3.0_dp-(i-1)**3/3.0_dp)
        !sol%val(i+(j-1)*10,2)=(j**3/3.0_dp-(j-1)**3/3.0_dp)
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

!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!

  !call init(xL,xR,yL,yR,nx,ny,nvar,cfl,tf,fs,namefile,mesh,sol,str_equa,str_flux,str_time_scheme,order)
  !call buildmesh(xL,xR,yL,yR,nx,ny,mesh)
  !call IC(mesh,sol)

  !order=3
  !k=11
  !dx=mesh%cell(k)%dx
  !dy=mesh%cell(k)%dy
  !func => evaluate2
  
  !call reconstruct2(mesh,sol,k,order)
  !print*,mesh%cell(k)%polCoef(:,1)
  !call quadrature3(func,mesh,sol,order,1,k,ul)
  !call evaluate2(mesh,sol,k,order,mesh%cell(k)%xc-dx/2.0_dp,mesh%cell(k)%yc,ul)

  !print*,"En x = ",mesh%cell(k)%xc-dx/2.0_dp," y = ",mesh%cell(k)%yc
  !print*,ul

  !deallocate(mesh%cell(k)%polCoef)
  
  !A(1,1)=0.0_dp
  !A(1,2)=-10.0_dp
  !A(2,1)=10.0_dp
  !A(2,2)=-10.0_dp
  !A(3,1)=10.0_dp
  !A(3,2)=0.0_dp
  !b(1)=10.0_dp
  !b(2)=10.0_dp
  !b(3)=0.0_dp

  !call QR(A,b,vect)
  !print*,vect

  contains

    subroutine init_FV(str_equa,str_flux,str_time_scheme,order,f_equa,flux,speed,time_scheme,quad_t,quad_c_alpha,quad_reconstruct)
    character(len=20), intent(in) :: str_equa,str_flux,str_time_scheme
    integer, intent(in) :: order
    procedure (sub_f), pointer, intent(out) :: f_equa
    procedure (sub_flux), pointer, intent(out) :: flux
    procedure (sub_speed), pointer, intent(out) :: speed
    procedure (sub_time), pointer, intent(out) :: time_scheme
    procedure (quadrature_t), pointer, intent(out) :: quad_t
    procedure (quadrature_c_alpha), pointer, intent(out) :: quad_c_alpha
    procedure (quadrature_reconstruction), pointer, intent(out) :: quad_reconstruct

    select case (trim(str_equa))
    case ('transport')
       f_equa => f_transport
    case ('euler')
       f_equa => f_euler
    case default
       print*,trim(str_equa)," equation not implemented"
       call exit()
    end select

    select case (trim(str_flux))
    case ('godunov')
       flux => flux_godunov
       speed => speed_godunov
    case ('HLL')
       flux => flux_HLL
       speed => speed_HLL
    case default
       print*,trim(str_flux)," flux not implemented"
       call exit()
    end select

    select case (trim(str_time_scheme))
    case ('euler_exp')
       time_scheme => euler_exp
    case ('SSPRK2')
       time_scheme => SSPRK2
    case ('SSPRK3')
       time_scheme => SSPRK3
    case default
       print*,trim(str_time_scheme)," time scheme not implemented"
       call exit()
    end select

    quad_t => quadrature3_t
    select case (order)
    case (1:3)
       quad_c_alpha => quadrature3_c_alpha
       quad_reconstruct => quadrature3_reconstruction
    case default
       print*,"Space order too high, no good enough quadrature implemented"
       call exit()
    end select

    return
  end subroutine init_FV

  subroutine calculation(mesh,sol,order,cfl,tf,fs,namefile,f_equa,flux,speed,time_scheme,quad_t,quad_c_alpha,quad_reconstruct)
    type(meshStruct), intent(inout) :: mesh
    type(solStruct), intent(inout) :: sol
    integer, intent(in) :: order
    real(dp), intent(in) :: cfl,tf
    integer, intent(in) :: fs
    character(len=20),intent(in) :: namefile
    procedure (sub_f), pointer, intent(in) :: f_equa
    procedure (sub_flux), pointer, intent(in) :: flux
    procedure (sub_speed), pointer, intent(in) :: speed
    procedure (sub_time), pointer, intent(in) :: time_scheme
    procedure (quadrature_t), pointer, intent(in) :: quad_t
    procedure (quadrature_c_alpha), pointer, intent(in) :: quad_c_alpha
    procedure (quadrature_reconstruction), pointer, intent(in) :: quad_reconstruct
    integer :: n
    real(dp) :: t
    
    t=0.0_dp
    n=1
    do while (t<tf)
       call time_scheme(mesh,sol,f_equa,flux,speed,order,cfl,t,tf,quad_c_alpha,quad_reconstruct)
       if (mod(n,fs)==0) then
          call userSol(t,mesh,sol,quad_t)
          call writeSol(mesh,sol,namefile,n/fs)
          call print(t,n)
       endif
       n=n+1
    enddo

  end subroutine calculation
    
  end program test
