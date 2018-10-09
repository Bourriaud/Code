program main
  
  use inout
  use efficiency
  use types
  use phys
  use FV
  
  implicit none

  type(meshStruct) :: mesh
  type(solStruct) :: sol
  real :: xL,xR,yL,yR,dx,dy,dt,tf,error
  integer :: fs
  character(len=20) :: namefile

  call init(xL,xR,yL,yR,dx,dy,dt,tf,fs,namefile,mesh,sol)
  call buildmesh(xL,xR,yL,yR,mesh)
  call IC(mesh,sol)
  call BC(xL,xR,yL,yR,mesh)
  call userSol(0.,mesh,sol)
  call writeSol(mesh,sol,namefile,0)
  
  call calculation(mesh,sol,dx,dy,dt,tf,fs,namefile)
  
  call errorL1(sol%val(:,:,1),sol%user(:,:,1),error)
  print*, "errorL1 = ",error
  call errorL2(sol%val(:,:,1),sol%user(:,:,1),error)
  print*, "errorL2 = ",error

  deallocate(mesh%X,mesh%Y,mesh%CX,mesh%CY,mesh%boundType,mesh%bound,sol%val,sol%user,sol%name,sol%nameUser)
  
end program main
