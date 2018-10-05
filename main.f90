program main
  
  use inout
  use efficiency
  use types
  
  implicit none

  type(meshStruct) :: mesh
  type(solStruct) :: sol
  real :: xL,xR,yL,yR
  real :: dx,dy
  character(len=20) :: namefile
  integer :: i,j
  
  call init(xL,xR,yL,yR,dx,dy,namefile,mesh,sol)
  allocate(mesh%X(0:mesh%nx,0:mesh%ny),mesh%Y(0:mesh%nx,0:mesh%ny),sol%val(0:mesh%nx,0:mesh%ny,sol%nsol))

  call buildmesh(xL,xR,yL,yR,mesh)

  do i=0,mesh%nx
     do j=0,mesh%ny
        sol%val(i,j,1)=1.
     enddo
  enddo
  call exactsol(mesh,sol%val(:,:,2))

  call writeSol(mesh,sol,namefile)

  deallocate(mesh%X,mesh%Y,sol%val,sol%name)
  
end program main
