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
  integer :: nvar,fs,fp,verbosity,order,i,k,f_adapt,recursivity,total_cell,average_cell
  integer :: level,minlevel,maxlevel,dim
  logical :: period
  integer, dimension(:), allocatable :: L_var_criteria,order_pc
  real(dp), dimension(:), allocatable :: L_eps
  character(len=20) :: config_file,test_case,namefile,str_equa,str_flux,str_time_scheme
  character(len=20) :: str_fn_adapt,str_exactSol,exact_file,restart_file
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
  call init(config_file,test_case,restart_file,xL,xR,yL,yR,level,nvar,cfl,tf,fs,fp,namefile,verbosity,sol, &
       str_equa,str_flux,str_time_scheme,order,period,L_str_criteria,L_var_criteria,L_eps, &
       gauss_point,gauss_weight,str_exactSol,exact_file,bool_AMR,str_fn_adapt,f_adapt,recursivity,order_pc)
  call buildP4EST(level,connectivity,p4est)
  call buildMesh_P4EST(p4est,xL,xR,yL,yR,gauss_point,order,mesh,sol,quadrants,period)
  !call buildmesh(xL,xR,yL,yR,64,64,gauss_point,order,mesh,sol)
  call init_FV(test_case,str_equa,str_flux,str_time_scheme,str_fn_adapt,IC_func, &
       BC,exactSol,f_equa,flux,speed,time_scheme,sol,fn_adapt,cfl,str_exactSol,exact_file,dim)
  if (restart_file=="none") then
     call IC(IC_func,mesh,sol,order)
  else
     call IC_restart(restart_file,mesh,sol)
  endif
  call BC(nvar,mesh)

  if (bool_AMR) then
     do i=0,recursivity-1
        do k=1,mesh%nc
           call reconstruct1(mesh,sol,k,order,gauss_weight,period,mesh%cell(k)%polCoef)
        enddo
        call adapt(fn_adapt,p4est,quadrants,mesh,sol,level+i,order,gauss_weight,gauss_point, &
             period,minlevel,maxlevel)
        call buildMesh_P4EST(p4est,xL,xR,yL,yR,gauss_point,order,mesh,sol,quadrants,period)
        call new_sol(mesh,order,quadrants,sol)
        if (restart_file=="none") then
           call IC(IC_func,mesh,sol,order)
        else
           !call IC_restart(restart_file,mesh,sol)
           print*,"Restart not possible with AMR"
           call exit()
        endif
        call BC(nvar,mesh)
     enddo
  else
     do k=1,mesh%nc
        call reconstruct1(mesh,sol,k,order,gauss_weight,period,mesh%cell(k)%polCoef)
     enddo
  endif

  if(verbosity>0) then
     call write_output_header(test_case,xL,xR,yL,yR,level,cfl,tf,namefile, &
       str_equa,str_flux,str_time_scheme,order,L_str_criteria, &
       bool_AMR,str_fn_adapt,f_adapt,recursivity,minlevel,maxlevel)
     call write_orders_header()
  endif
  call userSol(0.0_dp,mesh,sol,str_equa,exactSol,exact_file,dim)
  call writeSol(mesh,sol,namefile,0)
  call writeSol2(mesh,sol,namefile,0)
  call calculation(mesh,sol,level,order,cfl,tf,fs,fp,namefile,verbosity,str_equa, &
       f_equa,flux,speed,time_scheme,exactSol,order_pc, &
       L_str_criteria,L_var_criteria,L_eps,gauss_weight,gauss_point,period, &
       bool_AMR,fn_adapt,f_adapt,recursivity,total_cell,average_cell,exact_file,dim)

  select case (trim(str_exactSol))
  case ('sinus','sinus_dis','vortex','test','file')
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
     close(12)
     close(17)
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
    case ('Test')
       IC_func => IC_func_test
       BC => BC_test
       dim=1
    case ('Test2')
       IC_func => IC_func_test
       BC => BC_test
       dim=1
    case ('Step_burgers')
       IC_func => IC_func_step_burgers
       BC => BC_step_burgers
       dim=1
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
    case ('Sod_is')
       IC_func => IC_func_sod_is
       BC => BC_sod_is
       dim=1
    case ('Shu')
       IC_func => IC_func_shu
       BC => BC_shu
       dim=2
    case ('123')
       IC_func => IC_func_123
       BC => BC_123
       dim=2
    case ('Blastwave')
       IC_func => IC_func_blastwave
       BC => BC_blastwave
       dim=1
    case ('Vortex')
       IC_func => IC_func_vortex
       BC => BC_vortex
       dim=2
    case ('RP2D_3')
       IC_func => IC_func_RP2D_3
       BC => BC_RP2D_3
       dim=2
    case ('Riemann_M1')
       IC_func => IC_func_riemann_M1
       BC => BC_riemann_M1
       dim=1
    case default
       print*,"Test case ",trim(test_case)," not implemented"
       call exit()
    end select
    
    select case (trim(str_equa))
    case ('advection')
       f_equa => f_transport
       do i=1,sol%nvar
          sol%conserv_var(i,1:2)=i
       enddo
    case ('burgers')
       f_equa => f_burgers
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
    case ('euler_is')
       f_equa => f_euler_is
       sol%conserv_var(1,1:2)=1
       sol%conserv_var(2,1)=2
       sol%conserv_var(2,2)=3
       sol%conserv_var(3,1)=2
       sol%conserv_var(3,2)=3
    case ('M1')
       f_equa => f_M1
       sol%conserv_var(1,1:2)=1
       sol%conserv_var(2,1)=2
       sol%conserv_var(2,2)=3
       sol%conserv_var(3,1)=2
       sol%conserv_var(3,2)=3
    case default
       print*,trim(str_equa)," equation not implemented"
       call exit()
    end select

    select case (trim(str_flux))
    case ('godunov')
       select case (trim(str_equa))
       case('advection')
          flux => flux_godunov_adv
       case('burgers')
          flux => flux_godunov_bur
       case default
          print*,trim(str_flux)," is not valid for equation ",trim(str_equa)
          call exit()
       end select
       speed => speed_godunov
    case ('HLL')
       select case (trim(str_equa))
       case('euler')
          flux => flux_HLL
          speed => speed_HLL
       case('euler_is')
          flux => flux_HLL_is
          speed => speed_HLL_is
       case default
          print*,trim(str_flux)," is not valid for equation ",trim(str_equa)
          call exit()
       end select
    case ('HLLc')
       select case (trim(str_equa))
       case('euler')
          flux => flux_HLLc
          speed => speed_HLLc
       case('euler_is')
          flux => flux_HLLc_is
          speed => speed_HLLc_is
       case default
          print*,trim(str_flux)," is not valid for equation ",trim(str_equa)
          call exit()
       end select
    case ('rusanov')
       select case (trim(str_equa))
       case('euler')
          flux => flux_rusanov
          speed => speed_rusanov
       case('euler_is')
          flux => flux_rusanov_is
          speed => speed_rusanov_is
       case('M1')
          flux => flux_rusanov_M1
          speed => speed_rusanov_M1
       case default
          print*,trim(str_flux)," is not valid for equation ",trim(str_equa)
          call exit()
       end select
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
       cfl=cfl*0.5_dp
       !cfl=cfl*1.50818004918983_dp*0.5_dp
    case ('SSPRK5')
       time_scheme => SSPRK5
       cfl=cfl*0.5_dp
       !cfl=cfl*3.39533683277420_dp*0.5_dp
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
    case ('sod2D')
       fn_adapt => adapt_sod2D
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
    case ('test')
       exactSol => exactSol_test
       exact_file="none"
    case ('test2')
       exactSol => exactSol_test2
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

  subroutine calculation(mesh,sol,level,order,cfl,tf,fs,fp,namefile,verbosity,str_equa, &
       f_equa,flux,speed,time_scheme,exactSol,order_pc, &
       L_str_criteria,L_var_criteria,L_eps,gauss_weight,gauss_point,period, &
       bool_AMR,fn_adapt,f_adapt,recursivity,total_cell,average_cell,exact_file,dim)
    type(meshStruct), intent(inout) :: mesh
    type(solStruct), intent(inout) :: sol
    integer, intent(in) :: level,order,fs,fp,verbosity,f_adapt,recursivity,dim
    real(dp), intent(in) :: cfl,tf
    character(len=20),intent(in) :: namefile,str_equa,exact_file
    procedure (sub_f), pointer, intent(in) :: f_equa
    procedure (sub_flux), pointer, intent(in) :: flux
    procedure (sub_speed), pointer, intent(in) :: speed
    procedure (sub_time), pointer, intent(in) :: time_scheme
    procedure (sub_exactsol), pointer, intent(in) :: exactSol
    integer, dimension(:), intent(inout) :: order_pc
    procedure (sub_adapt), pointer, intent(in) :: fn_adapt
    character(len=20), dimension(:), intent(in) :: L_str_criteria
    integer, dimension(:), intent(in) :: L_var_criteria
    real(dp), dimension(:), intent(in) :: L_eps
    real(dp), dimension(:), intent(in) :: gauss_weight,gauss_point
    logical, intent(in) :: period,bool_AMR
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
            L_str_criteria,L_var_criteria,L_eps,gauss_weight,period,verbosity,order_pc)
       if (bool_AMR.and.mod(n,f_adapt)==0) then
          do i=0,recursivity-1
             call adapt(fn_adapt,p4est,quadrants,mesh,sol,level+i,order,gauss_weight,gauss_point, &
                  period,minlevel,maxlevel)
             call buildMesh_P4EST(p4est,xL,xR,yL,yR,gauss_point,order,mesh,sol,quadrants,period)
             call new_sol(mesh,order,quadrants,sol)
             call BC(nvar,mesh)
          enddo
       endif
       if (mod(n,fs)==0.or.t>=tf) then
          call userSol(t,mesh,sol,str_equa,exactSol,exact_file,dim)
          call writeSol(mesh,sol,namefile,nout)
          call writeSol2(mesh,sol,namefile,nout)
          if (verbosity>0) then
             call write_output_calculation(t,n,mesh%nc,total_quantities)
             call write_orders(mesh,order,n)
          endif
          print*,"Solution saved, iteration",n
          nout=nout+1
       endif
       if (mod(n,fp)==0.or.t>=tf) then
          call print(mesh,sol,t,n,total_quantities)
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
