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
  real(dp), dimension(:), allocatable :: gauss_point,gauss_weight

  call init(xL,xR,yL,yR,nx,ny,nvar,cfl,tf,fs,namefile,mesh,sol,str_equa,str_flux,str_time_scheme, &
       order,L_str_criteria,L_var_criteria,gauss_point,gauss_weight)
  call init_FV(str_equa,str_flux,str_time_scheme,f_equa,flux,speed,time_scheme)
  call buildmesh(xL,xR,yL,yR,nx,ny,gauss_point,mesh)
  call IC(mesh,sol,gauss_weight)
  call BC(nx,ny,nvar,mesh)
  call userSol(0.0_dp,mesh,sol,gauss_weight)
  call writeSol(mesh,sol,namefile,0)
  call calculation(mesh,sol,order,cfl,tf,fs,namefile,f_equa,flux,speed,time_scheme, &
       L_str_criteria,L_var_criteria,gauss_weight)

  call userSol(tf,mesh,sol,gauss_weight)
  call errorL1(mesh,sol%val(:,2),sol%user(:,1),error)
  call errorL2(mesh,sol%val(:,2),sol%user(:,1),error)
  
  deallocate(mesh%node,mesh%edge,mesh%cell)
  deallocate(sol%val,sol%user,sol%name,sol%nameUser)
  deallocate(L_str_criteria,L_var_criteria,gauss_point,gauss_weight)

contains

  subroutine init_FV(str_equa,str_flux,str_time_scheme,f_equa,flux,speed,time_scheme)
    character(len=20), intent(in) :: str_equa,str_flux,str_time_scheme
    procedure (sub_f), pointer, intent(out) :: f_equa
    procedure (sub_flux), pointer, intent(out) :: flux
    procedure (sub_speed), pointer, intent(out) :: speed
    procedure (sub_time), pointer, intent(out) :: time_scheme

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

    return
  end subroutine init_FV

  subroutine calculation(mesh,sol,order,cfl,tf,fs,namefile,f_equa,flux,speed,time_scheme, &
       L_str_criteria,L_var_criteria,gauss_weight)
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
    real(dp), dimension(:), intent(in) :: gauss_weight
    integer :: n
    real(dp) :: t
    
    t=0.0_dp
    n=1
    call check_conservativity(mesh,sol)
    do while (t<tf)
       call time_scheme(mesh,sol,f_equa,flux,speed,order,cfl,t,tf,L_str_criteria,L_var_criteria,gauss_weight)
       if (mod(n,fs)==0) then
          call userSol(t,mesh,sol,gauss_weight)
          call writeSol(mesh,sol,namefile,n/fs)
          call print(mesh,sol,t,n)
       endif
       n=n+1
    enddo

  end subroutine calculation
  
end program main
