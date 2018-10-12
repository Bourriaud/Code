program main
  
  use inout
  use efficiency
  use types
  use phys
  use FV
  
  implicit none

  type(meshStruct) :: mesh
  type(solStruct) :: sol
  real :: xL,xR,yL,yR,dt,tf,error
  integer :: nx,ny,nvar,fs
  character(len=20) :: namefile

  call init(xL,xR,yL,yR,nx,ny,nvar,dt,tf,fs,namefile,mesh,sol)
  call buildmesh(xL,xR,yL,yR,nx,ny,mesh)
  call IC(mesh,sol)
  call BC(nx,ny,nvar,mesh)
  call userSol(0.,mesh,sol)
  call writeSol(mesh,sol,namefile,0)
  
  call calculation(mesh,sol,dt,tf,fs,namefile)
  
  call errorL1(sol%val(:,1),sol%user(:,1),error)
  print*, "errorL1 = ",error
  call errorL2(sol%val(:,1),sol%user(:,1),error)
  print*, "errorL2 = ",error

  do fs=1,mesh%nc
     if (mesh%cell(fs)%xc==0.505) then
        print*,sol%val(fs,1),mesh%cell(fs)%yc
     endif
  enddo
  
  deallocate(mesh%node,mesh%cell,sol%val,sol%user,sol%name,sol%nameUser)
  
end program main
