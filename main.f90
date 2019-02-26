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
  real(dp) :: xL,xR,yL,yR,cfl,tf,eL1,eL2,tstart,tfinish
  integer :: nvar,fs,verbosity,order,i,f_adapt,recursivity,total_cell,average_cell
  integer :: level,minlevel,maxlevel,dim
  integer, dimension(:), allocatable :: L_var_criteria
  real(dp), dimension(:), allocatable :: L_eps
  character(len=20) :: config_file,test_case,namefile,str_equa,str_flux,str_time_scheme
  character(len=20) :: str_fn_adapt,str_exactSol,exact_file
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
  type(c_ptr) :: p4est,connectivity,quadrants
  logical :: bool_AMR

  call CPU_TIME(tstart)
  
  call get_config(config_file)
  call init(config_file,test_case,xL,xR,yL,yR,level,nvar,cfl,tf,fs,namefile,verbosity,sol, &
       str_equa,str_flux,str_time_scheme,order,L_str_criteria,L_var_criteria,L_eps, &
       gauss_point,gauss_weight,str_exactSol,exact_file,bool_AMR,str_fn_adapt,f_adapt,recursivity)
  call buildP4EST(level,connectivity,p4est)
  call buildMesh_P4EST(p4est,xL,xR,yL,yR,gauss_point,order,mesh,sol,quadrants)
  !call buildmesh(xL,xR,yL,yR,64,64,gauss_point,order,mesh,sol)
  call init_FV(test_case,str_equa,str_flux,str_time_scheme,str_fn_adapt,IC_func, &
       BC,exactSol,f_equa,flux,speed,time_scheme,sol,fn_adapt,cfl,str_exactSol,exact_file,dim)
  call IC(IC_func,mesh,sol,order)
  call BC(nvar,mesh)
  
  if (bool_AMR) then
     do i=0,recursivity-1
        call adapt(fn_adapt,p4est,quadrants,mesh,sol,level+i,order,gauss_weight,gauss_point,minlevel,maxlevel)
        call buildMesh_P4EST(p4est,xL,xR,yL,yR,gauss_point,order,mesh,sol,quadrants)
        call new_sol(mesh,quadrants,sol)
        call IC(IC_func,mesh,sol,order)
        call BC(nvar,mesh)
     enddo
  endif

  if(verbosity>0) then
     call write_output_header(test_case,xL,xR,yL,yR,level,cfl,tf,namefile, &
       str_equa,str_flux,str_time_scheme,order,L_str_criteria, &
       bool_AMR,str_fn_adapt,f_adapt,recursivity,minlevel,maxlevel)
  endif
  call userSol(0.0_dp,mesh,sol,str_equa,exactSol,exact_file,dim)
  call writeSol(mesh,sol,namefile,0)
  call calculation(mesh,sol,level,order,cfl,tf,fs,namefile,verbosity,str_equa, &
       f_equa,flux,speed,time_scheme,exactSol, &
       L_str_criteria,L_var_criteria,L_eps,gauss_weight,gauss_point, &
       bool_AMR,fn_adapt,f_adapt,recursivity,total_cell,average_cell,exact_file,dim)

  select case (trim(str_exactSol))
  case ('sinus')
     call errorL1(mesh,sol%val(:,1),sol%user(:,1),eL1)
     call errorL2(mesh,sol%val(:,1),sol%user(:,1),eL2)
  case ('sinus_dis')
     call errorL1(mesh,sol%val(:,1),sol%user(:,1),eL1)
     call errorL2(mesh,sol%val(:,1),sol%user(:,1),eL2)
  case ('vortex')
     call errorL1(mesh,sol%val(:,1),sol%user(:,1),eL1)
     call errorL2(mesh,sol%val(:,1),sol%user(:,1),eL2)
  case ('file')
     call errorL1(mesh,sol%val(:,1),sol%user(:,1),eL1)
     call errorL2(mesh,sol%val(:,1),sol%user(:,1),eL2)
  case default
     print*,"No analytical solution for this configuration"
     eL1=-1.0_dp
     eL2=-1.0_dp
  end select
  
  deallocate(mesh%node,mesh%edge,mesh%cell)
  deallocate(sol%val,sol%user,sol%name,sol%var_user,sol%name_user,sol%conserv_var)
  deallocate(L_str_criteria,L_var_criteria,L_eps)
  deallocate(gauss_point,gauss_weight)
  call p4_destroy(connectivity,p4est)
  
  call CPU_TIME(tfinish)
  print*,"Calculation completed, time of execution : ",tfinish-tstart,"seconds"
  if (verbosity>0) then
     call write_output_summary(tfinish-tstart,eL1,eL2,total_cell,average_cell)
  endif

contains

  subroutine init_FV(test_case,str_equa,str_flux,str_time_scheme,str_fn_adapt,IC_func, &
       BC,exactSol,f_equa,flux,speed,time_scheme,sol,fn_adapt,cfl,str_exactSol,exact_file,dim)
    character(len=20), intent(in) :: test_case,str_equa,str_flux,str_time_scheme
    character(len=20), intent(in) :: str_fn_adapt,str_exactSol
    procedure (sub_IC), pointer, intent(out) :: IC_func
    procedure (sub_BC), pointer, intent(out) :: BC
    procedure (sub_exactsol), pointer, intent(out) :: exactSol
    procedure (sub_f), pointer, intent(out) :: f_equa
    procedure (sub_flux), pointer, intent(out) :: flux
    procedure (sub_speed), pointer, intent(out) :: speed
    procedure (sub_time), pointer, intent(out) :: time_scheme
    type(solStruct), intent(inout) :: sol
    procedure (sub_adapt), pointer, intent(out) :: fn_adapt
    real(dp), intent(inout) :: cfl
    character(len=20), intent(inout) :: exact_file
    integer, intent(out) :: dim
    integer :: i

    select case (trim(test_case))
    case ('Sinus')
       IC_func => IC_func_sinus
       BC => BC_sinus
       dim=2
    case ('Sinus_dis')
       IC_func => IC_func_sinus_dis
       BC => BC_sinus_dis
       dim=2
    case ('Sod')
       IC_func => IC_func_sod
       BC => BC_sod
       dim=1
    case ('Sod_2D')
       IC_func => IC_func_sod_2D
       BC => BC_sod_2D
       dim=2
    case ('Sod_mod')
       IC_func => IC_func_sod_mod
       BC => BC_sod_mod
       dim=1
    case ('Shu')
       IC_func => IC_func_shu
       BC => BC_shu
       dim=2
    case ('123')
       IC_func => IC_func_123
       BC => BC_123
       dim=2
    case ('Vortex')
       IC_func => IC_func_vortex
       BC => BC_vortex
       dim=2
    case ('RP2D_3')
       IC_func => IC_func_RP2D_3
       BC => BC_RP2D_3
       dim=2
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
       cfl=cfl*0.5_dp
    case ('SSPRK2')
       time_scheme => SSPRK2
       cfl=cfl*0.5_dp
    case ('SSPRK3')
       time_scheme => SSPRK3
       cfl=cfl*0.5_dp
    case ('SSPRK4')
       time_scheme => SSPRK4
       cfl=cfl*1.50818004918983_dp*0.5_dp
    case ('SSPRK5')
       time_scheme => SSPRK5
       cfl=cfl*3.39533683277420_dp*0.5_dp
    case default
       print*,trim(str_time_scheme)," time scheme not implemented"
       call exit()
    end select

    select case (trim(str_fn_adapt))
    case ('zone')
       fn_adapt => adapt_zone
    case ('sinus')
       fn_adapt => adapt_sinus
    case ('vortex')
       fn_adapt => adapt_vortex
    case ('sod')
       fn_adapt => adapt_sod
    case default
       print*,trim(str_fn_adapt)," adaptation function not implemented"
       call exit()
    end select

    select case (trim(str_exactSol))
    case ('sinus')
       exactSol => exactSol_sinus
       exact_file="none"
    case ('sinus_dis')
       exactSol => exactSol_sinus_dis
       exact_file="none"
    case ('vortex')
       exactSol => exactSol_vortex
       exact_file="none"
    case ('file')
       exactSol => exactSol_none
    case ('none')
       exactSol => exactSol_none
       exact_file="none"
    case default
       print*,"Exact solution ",trim(str_exactSol)," not implemented"
       call exit()
    end select

    return
  end subroutine init_FV

  subroutine calculation(mesh,sol,level,order,cfl,tf,fs,namefile,verbosity,str_equa, &
       f_equa,flux,speed,time_scheme,exactSol, &
       L_str_criteria,L_var_criteria,L_eps,gauss_weight,gauss_point, &
       bool_AMR,fn_adapt,f_adapt,recursivity,total_cell,average_cell,exact_file,dim)
    type(meshStruct), intent(inout) :: mesh
    type(solStruct), intent(inout) :: sol
    integer, intent(in) :: level,order,fs,verbosity,f_adapt,recursivity,dim
    real(dp), intent(in) :: cfl,tf
    character(len=20),intent(in) :: namefile,str_equa,exact_file
    procedure (sub_f), pointer, intent(in) :: f_equa
    procedure (sub_flux), pointer, intent(in) :: flux
    procedure (sub_speed), pointer, intent(in) :: speed
    procedure (sub_time), pointer, intent(in) :: time_scheme
    procedure (sub_exactsol), pointer, intent(in) :: exactSol
    procedure (sub_adapt), pointer, intent(in) :: fn_adapt
    character(len=20), dimension(:), intent(in) :: L_str_criteria
    integer, dimension(:), intent(in) :: L_var_criteria
    real(dp), dimension(:), intent(in) :: L_eps
    real(dp), dimension(:), intent(in) :: gauss_weight,gauss_point
    logical, intent(in) :: bool_AMR
    integer, intent(out) :: total_cell,average_cell
    integer :: i,n,nout,minlevel,maxlevel
    real(dp) :: t
    real(dp), dimension(:), allocatable :: total_quantities
    character(len=30) :: AMRfile

    allocate(total_quantities(sol%nvar))
    
    t=0.0_dp
    n=1
    nout=1
    total_cell=mesh%nc
    call print(mesh,sol,0.0_dp,0,total_quantities)
    if (verbosity>0) then
       call write_output_calculation(0.0_dp,0,mesh%nc,total_quantities)
       write(AMRfile,'(A,A,A)')'./results/',trim(namefile),'_AMR.txt'
       open(13,file=AMRfile,form="formatted")
       write(13,'(i8,i8)')0,mesh%nc
    endif
    
    do while (t<tf)
       call time_scheme(mesh,sol,str_equa,f_equa,flux,speed,order,cfl,t,n,tf, &
            L_str_criteria,L_var_criteria,L_eps,gauss_weight,verbosity)
       if (bool_AMR.and.mod(n,f_adapt)==0) then
          do i=0,recursivity-1
             call adapt(fn_adapt,p4est,quadrants,mesh,sol,level+i,order,gauss_weight,gauss_point,minlevel,maxlevel)
             call buildMesh_P4EST(p4est,xL,xR,yL,yR,gauss_point,order,mesh,sol,quadrants)
             call new_sol(mesh,quadrants,sol)
             call BC(nvar,mesh)
          enddo
       endif
       if (mod(n,fs)==0.or.t>=tf) then
          call userSol(t,mesh,sol,str_equa,exactSol,exact_file,dim)
          call writeSol(mesh,sol,namefile,nout)
          call print(mesh,sol,t,n,total_quantities)
          if (verbosity>0) call write_output_calculation(t,n,mesh%nc,total_quantities)
          nout=nout+1
       endif
       if (verbosity>0) write(13,'(i8,i8)')n,mesh%nc
       total_cell=total_cell+mesh%nc
       n=n+1
    enddo
    average_cell=total_cell/n

    if (verbosity>0) close(13)
    deallocate(total_quantities)
    
    return
  end subroutine calculation
  
end program main
