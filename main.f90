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
  character(len=20) :: namefile,str_equa,str_flux,str_time_scheme

  call init(xL,xR,yL,yR,nx,ny,nvar,cfl,tf,fs,namefile,mesh,sol,str_equa,str_flux,str_time_scheme,order)
  call buildmesh(xL,xR,yL,yR,nx,ny,mesh)
  call IC(mesh,sol)
  call BC(nx,ny,nvar,mesh)
  call userSol(0.0_dp,mesh,sol)
  call writeSol(mesh,sol,namefile,0)
  call calculation(mesh,sol,str_equa,str_flux,str_time_scheme,order,cfl,tf,fs,namefile)

  call errorL1(mesh,sol%val(:,2),sol%user(:,1),error)
  print*, "errorL1 = ",error
  call errorL2(mesh,sol%val(:,2),sol%user(:,1),error)
  print*, "errorL2 = ",error
  
  deallocate(mesh%node,mesh%cell,sol%val,sol%user,sol%name,sol%nameUser)

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
  
end program main
