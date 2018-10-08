program main
  
  use inout
  use efficiency
  use types
  use phys
  use FV
  
  implicit none

  type(meshStruct) :: mesh
  type(solStruct) :: sol
  real :: xL,xR,yL,yR,dx,dy,dt,tf
  character(len=20) :: namefile

  call init(xL,xR,yL,yR,dx,dy,dt,tf,namefile,mesh,sol)
  call buildmesh(xL,xR,yL,yR,mesh)
  call IC(mesh,sol)
  call BC(xL,xR,yL,yR,mesh)
  call userSol(mesh,sol)
  call writeSol(mesh,sol,namefile,0)
  
  call calculation(mesh,sol,dx,dy,dt,tf,namefile)

  deallocate(mesh%X,mesh%Y,mesh%boundType,mesh%bound,sol%val,sol%user,sol%name,sol%nameUser)
  
end program main
