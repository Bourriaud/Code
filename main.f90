program main

  use constant
  use inout
  use efficiency
  use types
  use phys
  use FV
  use time
  use reconstruction
  
  implicit none

  type(meshStruct) :: mesh
  type(solStruct) :: sol
  real(dp) :: xL,xR,yL,yR,cfl,tf,error
  integer :: nx,ny,nvar,fs,order
  integer, dimension(:), allocatable :: L_var_criteria
  character(len=20) :: namefile,str_equa,str_flux,str_time_scheme
  character(len=20), dimension(:), allocatable :: L_str_criteria
  procedure (sub_f), pointer :: f_equa
  procedure (sub_flux), pointer :: flux
  procedure (sub_speed), pointer :: speed
  procedure (sub_time), pointer :: time_scheme
  procedure (quadrature_t), pointer :: quad_t
  procedure (quadrature_c_alpha), pointer :: quad_c_alpha
  procedure (quadrature_reconstruction), pointer :: quad_reconstruct

  call init(xL,xR,yL,yR,nx,ny,nvar,cfl,tf,fs,namefile,mesh,sol,str_equa,str_flux,str_time_scheme, &
       order,L_str_criteria,L_var_criteria)
  call init_FV(str_equa,str_flux,str_time_scheme,order,f_equa,flux,speed,time_scheme, &
       quad_t,quad_c_alpha,quad_reconstruct)
  call buildmesh(xL,xR,yL,yR,nx,ny,mesh)
  call IC(mesh,sol,quad_t)
  call BC(nx,ny,nvar,mesh)
  call userSol(0.0_dp,mesh,sol,quad_t)
  call writeSol(mesh,sol,namefile,0)
  call calculation(mesh,sol,order,cfl,tf,fs,namefile,f_equa,flux,speed,time_scheme, &
       L_str_criteria,L_var_criteria,quad_t,quad_c_alpha,quad_reconstruct)
  
  call errorL1(mesh,sol%val(:,2),sol%user(:,1),error)
  print*, "errorL1 = ",error
  call errorL2(mesh,sol%val(:,2),sol%user(:,1),error)
  print*, "errorL2 = ",error
  
  deallocate(mesh%node,mesh%edge,mesh%cell,sol%val,sol%user,sol%name,sol%nameUser,L_str_criteria,L_var_criteria)

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

  subroutine calculation(mesh,sol,order,cfl,tf,fs,namefile,f_equa,flux,speed,time_scheme, &
       L_str_criteria,L_var_criteria,quad_t,quad_c_alpha,quad_reconstruct)
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
    character(len=20), dimension(:), intent(in) :: L_str_criteria
    integer, dimension(:), intent(in) :: L_var_criteria
    procedure (quadrature_t), pointer, intent(in) :: quad_t
    procedure (quadrature_c_alpha), pointer, intent(in) :: quad_c_alpha
    procedure (quadrature_reconstruction), pointer, intent(in) :: quad_reconstruct
    integer :: n
    real(dp) :: t
    
    t=0.0_dp
    n=1
    do while (t<tf)
       call time_scheme(mesh,sol,f_equa,flux,speed,order,cfl,t,tf,quad_c_alpha,quad_reconstruct,L_str_criteria,L_var_criteria)
       if (mod(n,fs)==0) then
          call userSol(t,mesh,sol,quad_t)
          call writeSol(mesh,sol,namefile,n/fs)
          call print(t,n)
       endif
       n=n+1
    enddo

  end subroutine calculation
  
end program main
