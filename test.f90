program test

  use constant
  use reconstruction
  use types
  use inout

  implicit none

  type(meshStruct) :: mesh
  type(solStruct) :: sol
  real(dp) :: xL,xR,yL,yR,cfl,tf,error
  integer :: nx,ny,nvar,fs
  character(len=20) :: namefile,str_equa,str_flux,str_time_scheme
  real(dp), dimension(:), allocatable :: X,U,Xstencil,Ustencil
  real(dp), dimension(3) :: X2
  real(dp), dimension(2) :: U2
  integer :: i,j,k,kpos,order,dir,neigh
  real(dp), dimension(2) :: ul,ur
  real(dp) :: val
  

  call init(xL,xR,yL,yR,nx,ny,nvar,cfl,tf,fs,namefile,mesh,sol,str_equa,str_flux,str_time_scheme,order)
  call buildmesh(xL,xR,yL,yR,nx,ny,mesh)

  do i=1,100
     do j=1,100
        !sol%val(i+(j-1)*100,1)=((50.5_dp-i)*0.1_dp)**1
        sol%val(i+(j-1)*100,1)=((i-0.5_dp)*0.1_dp)**10
        sol%val(i+(j-1)*100,2)=((i-0.5_dp)*0.1_dp)**1
        !if(i<100)then
           !sol%val(i+(j-1)*100,1)=1.0_dp
        !else
           !sol%val(i+(j-1)*100,1)=0.0_dp
        !endif
     enddo
  enddo

  k=1
  order=1
  dir=1
  neigh=51
  !call extractDirection(mesh,sol,k,order,dir,1,X,U,kpos)
  !call buildStencil(X,U,kpos,order,Xstencil,Ustencil)
  !print*,Xstencil
  !print*,Ustencil
  !print*,U(kpos),kpos,U(98),98
  !call reconstruct(mesh,sol,k,neigh,order,dir,ul,ur)
  call reconstruct_boundary(mesh,sol,k,order,dir,ul)

!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!

  !X2(1)=0.0_dp
  !X2(2)=1.0_dp
  !X2(3)=2.0_dp
  !U2(1)=0.5_dp
  !U2(2)=1.5_dp

  !call evaluate(4.0_dp,2,X2,U2,val)
  !print*,val


end program test
