program main

  use constant
  use inout
  use efficiency
  use types
  use phys
  use FV
  use time
  use reconstruction
  use ICBC
  use AMR
  
  implicit none

  type(meshStruct) :: mesh
  type(solStruct) :: sol
  real(dp) :: xL,xR,yL,yR,cfl,tf,error
  integer :: level,nvar,fs,order,i,f_adapt
  integer, dimension(:), allocatable :: L_var_criteria
  real(dp), dimension(:), allocatable :: L_eps
  character(len=20) :: config_file,test_case,namefile,str_equa,str_flux,str_time_scheme
  character(len=20) :: str_fn_adapt
  character(len=20), dimension(:), allocatable :: L_str_criteria
  procedure (sub_IC), pointer :: IC_func
  procedure (sub_BC), pointer :: BC
  procedure (sub_exactsol), pointer :: exactSol
  procedure (sub_f), pointer :: f_equa
  procedure (sub_flux), pointer :: flux
  procedure (sub_speed), pointer :: speed
  procedure (sub_time), pointer :: time_scheme
  procedure (sub_adapt), pointer :: fn_adapt
  real(dp), dimension(:), allocatable :: gauss_point,gauss_weight
  type(c_ptr) :: p4est,quadrants
  logical :: bool_AMR

  call get_config(config_file)
  call init(config_file,test_case,xL,xR,yL,yR,level,nvar,cfl,tf,fs,namefile,sol, &
       str_equa,str_flux,str_time_scheme,order,L_str_criteria,L_var_criteria,L_eps, &
       gauss_point,gauss_weight,bool_AMR,str_fn_adapt,f_adapt)
  call buildP4EST(level,p4est)
  call buildMesh_P4EST(p4est,xL,xR,yL,yR,gauss_point,order,mesh,sol,quadrants)
  !call buildmesh(xL,xR,yL,yR,64,64,gauss_point,order,mesh,sol)
  call init_FV(test_case,str_equa,str_flux,str_time_scheme,str_fn_adapt,IC_func, &
       BC,exactSol,f_equa,flux,speed,time_scheme,sol,fn_adapt)
  call IC(IC_func,mesh,sol,order)
  call BC(nvar,mesh)

  if (bool_AMR) then
     do i=1,1
        call adapt(fn_adapt,p4est,quadrants,mesh,sol)
        call buildMesh_P4EST(p4est,xL,xR,yL,yR,gauss_point,order,mesh,sol,quadrants)
        call IC(IC_func,mesh,sol,order)
        call BC(nvar,mesh)
     enddo
  endif
  
  call userSol(0.0_dp,mesh,sol,str_equa,exactSol)
  call writeSol(mesh,sol,namefile,0)
  call calculation(mesh,sol,order,cfl,tf,fs,namefile,str_equa, &
       f_equa,flux,speed,time_scheme,exactSol, &
       L_str_criteria,L_var_criteria,L_eps,gauss_weight, &
       bool_AMR,fn_adapt,f_adapt)

  call userSol(tf,mesh,sol,str_equa,exactSol)
  select case (trim(test_case))
  case ('Sinus')
     call errorL1(mesh,sol%val(:,1),sol%user(:,1),error)
     call errorL2(mesh,sol%val(:,1),sol%user(:,1),error)
  case ('Sinus_dis')
     call errorL1(mesh,sol%val(:,1),sol%user(:,1),error)
     call errorL2(mesh,sol%val(:,1),sol%user(:,1),error)
  case ('Vortex')
     call errorL1(mesh,sol%val(:,1),sol%user(:,1),error)
     call errorL2(mesh,sol%val(:,1),sol%user(:,1),error)
  case default
     print*,"No analytical solution for this configuration"
  end select
  
  deallocate(mesh%node,mesh%edge,mesh%cell)
  deallocate(sol%val,sol%user,sol%name,sol%var_user,sol%name_user,sol%conserv_var)
  deallocate(L_str_criteria,L_var_criteria,L_eps)
  deallocate(gauss_point,gauss_weight)

contains

  subroutine init_FV(test_case,str_equa,str_flux,str_time_scheme,str_fn_adapt,IC_func, &
       BC,exactSol,f_equa,flux,speed,time_scheme,sol,fn_adapt)
    character(len=20), intent(in) :: test_case,str_equa,str_flux,str_time_scheme
    character(len=20), intent(in) :: str_fn_adapt
    procedure (sub_IC), pointer, intent(out) :: IC_func
    procedure (sub_BC), pointer, intent(out) :: BC
    procedure (sub_exactsol), pointer, intent(out) :: exactSol
    procedure (sub_f), pointer, intent(out) :: f_equa
    procedure (sub_flux), pointer, intent(out) :: flux
    procedure (sub_speed), pointer, intent(out) :: speed
    procedure (sub_time), pointer, intent(out) :: time_scheme
    type(solStruct), intent(inout) :: sol
    procedure (sub_adapt), pointer, intent(out) :: fn_adapt
    integer :: i

    select case (trim(test_case))
    case ('Sinus')
       IC_func => IC_func_sinus
       BC => BC_sinus
       exactSol => exactSol_sinus
    case ('Sinus_dis')
       IC_func => IC_func_sinus_dis
       BC => BC_sinus_dis
       exactSol => exactSol_sinus_dis
    case ('Sod')
       IC_func => IC_func_sod
       BC => BC_sod
       exactSol => exactSol_none
    case ('Sod_2D')
       IC_func => IC_func_sod_2D
       BC => BC_sod_2D
       exactSol => exactSol_none
    case ('Sod_mod')
       IC_func => IC_func_sod_mod
       BC => BC_sod_mod
       exactSol => exactSol_none
    case ('Shu')
       IC_func => IC_func_shu
       BC => BC_shu
       exactSol => exactSol_none
    case ('123')
       IC_func => IC_func_123
       BC => BC_123
       exactSol => exactSol_none
    case ('Vortex')
       IC_func => IC_func_vortex
       BC => BC_vortex
       exactSol => exactSol_vortex
    case ('RP2D_3')
       IC_func => IC_func_RP2D_3
       BC => BC_RP2D_3
       exactSol => exactSol_none
    case default
       print*,"Test case ",trim(test_case)," not implemented"
       call exit()
    end select
    
    select case (trim(str_equa))
    case ('transport')
       f_equa => f_transport
       do i=1,sol%nvar
          sol%conserv_var(i,1:2)=i
       enddo
    case ('euler')
       f_equa => f_euler
       sol%conserv_var(1,1:2)=1
       sol%conserv_var(2,1)=2
       sol%conserv_var(2,2)=3
       sol%conserv_var(3,1)=2
       sol%conserv_var(3,2)=3
       sol%conserv_var(4,1:2)=4
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

    select case (trim(str_fn_adapt))
    case ('sinus')
       fn_adapt => adapt_sinus
    case ('vortex')
       fn_adapt => adapt_vortex
    case default
       print*,trim(str_fn_adapt)," adaptation function not implemented"
       call exit()
    end select

    return
  end subroutine init_FV

  subroutine calculation(mesh,sol,order,cfl,tf,fs,namefile,str_equa, &
       f_equa,flux,speed,time_scheme,exactSol, &
       L_str_criteria,L_var_criteria,L_eps,gauss_weight, &
       bool_AMR,fn_adapt,f_adapt)
    type(meshStruct), intent(inout) :: mesh
    type(solStruct), intent(inout) :: sol
    integer, intent(in) :: order
    real(dp), intent(in) :: cfl,tf
    integer, intent(in) :: fs,f_adapt
    character(len=20),intent(in) :: namefile,str_equa
    procedure (sub_f), pointer, intent(in) :: f_equa
    procedure (sub_flux), pointer, intent(in) :: flux
    procedure (sub_speed), pointer, intent(in) :: speed
    procedure (sub_time), pointer, intent(in) :: time_scheme
    procedure (sub_exactsol), pointer, intent(in) :: exactSol
    procedure (sub_adapt), pointer, intent(in) :: fn_adapt
    character(len=20), dimension(:), intent(in) :: L_str_criteria
    integer, dimension(:), intent(in) :: L_var_criteria
    real(dp), dimension(:), intent(in) :: L_eps
    real(dp), dimension(:), intent(in) :: gauss_weight
    logical, intent(in) :: bool_AMR
    integer :: n
    real(dp) :: t
    
    t=0.0_dp
    n=1
    call check_conservativity(mesh,sol)
    do while (t<tf)
       call time_scheme(mesh,sol,str_equa,f_equa,flux,speed,order,cfl,t,n,tf, &
            L_str_criteria,L_var_criteria,L_eps,gauss_weight)
       if (bool_AMR.and.mod(n,f_adapt)==0) then
          call adapt(fn_adapt,p4est,quadrants,mesh,sol)
          call buildMesh_P4EST(p4est,xL,xR,yL,yR,gauss_point,order,mesh,sol,quadrants)
          call new_sol(mesh,quadrants,sol)
          call BC(nvar,mesh)
       endif
       if (mod(n,fs)==0.or.t>=tf) then
          call userSol(t,mesh,sol,str_equa,exactSol)
          call writeSol(mesh,sol,namefile,n/fs)
          call print(mesh,sol,t,n)
       endif
       n=n+1
    enddo

  end subroutine calculation
  
end program main
