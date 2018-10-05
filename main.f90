program main
  
  use inout
  use efficiency
  
  implicit none

  real, dimension(:), allocatable :: Ix,Iy
  real, dimension(:,:,:), allocatable :: sol
  real :: xL,xR,yL,yR
  integer :: nx,ny,nsol
  real :: dx,dy
  character(len=20) :: namefile
  character(len=20), dimension(:), allocatable :: solname
  integer :: i,j
  
  call init(xL,xR,yL,yR,nx,ny,dx,dy,nsol,namefile,solname)

  allocate(Ix(0:nx),Iy(0:ny),sol(0:nx,0:ny,nsol))

  call buildmesh(xL,xR,yL,yR,nx,ny,Ix,Iy)

  do i=0,nx
     do j=0,ny
        sol(i,j,1)=1.
     enddo
  enddo
  call exactsol(Ix,Iy,sol(:,:,2))

  call writeSol(Ix,Iy,nsol,sol,namefile,solname)

  deallocate(Ix,Iy,sol,solname)
  
end program main
