program test

  use constant
  use inout
  use efficiency
  use types
  use phys
  use FV
  use time
  use reconstruction
  use ICBC
  use limit
  
  implicit none

  type(meshStruct) :: mesh
  type(solStruct) :: sol
  real(dp) :: xL,xR,yL,yR,cfl,tf,error
  integer :: nx,ny,nvar,fs,order
  integer, dimension(:), allocatable :: L_var_criteria
  real(dp), dimension(:), allocatable :: L_eps
  character(len=20) :: config_file,test_case,namefile,str_equa,str_flux,str_time_scheme
  character(len=20), dimension(:), allocatable :: L_str_criteria
  procedure (sub_IC), pointer :: IC_func
  procedure (sub_BC), pointer :: BC
  procedure (sub_exactsol), pointer :: exactSol
  procedure (sub_f), pointer :: f_equa
  procedure (sub_flux), pointer :: flux
  procedure (sub_speed), pointer :: speed
  procedure (sub_time), pointer :: time_scheme
  real(dp), dimension(:), allocatable :: gauss_point,gauss_weight

  call get_config(config_file)
  call init(config_file,test_case,xL,xR,yL,yR,nx,ny,nvar,cfl,tf,fs,namefile,mesh,sol, &
       str_equa,str_flux,str_time_scheme,order,L_str_criteria,L_var_criteria,L_eps, &
       gauss_point,gauss_weight)
  call init_FV(test_case,str_equa,str_flux,str_time_scheme,IC_func,BC,exactSol, &
       f_equa,flux,speed,time_scheme,sol)
  call buildmesh2(xL,xR,yL,yR,nx,ny,gauss_point,mesh)
  call IC(IC_func,mesh,sol,order,gauss_point6,gauss_weight6)
  call BC(nx,ny,nvar,mesh)
  call userSol(0.0_dp,mesh,sol,str_equa,exactSol,gauss_weight)
  call writeSol(mesh,sol,namefile,0)
  call calculation(mesh,sol,order,cfl,tf,fs,namefile,str_equa, &
       f_equa,flux,speed,time_scheme,exactSol, &
       L_str_criteria,L_var_criteria,L_eps,gauss_weight)

  call userSol(tf,mesh,sol,str_equa,exactSol,gauss_weight)
  
  deallocate(mesh%node,mesh%edge,mesh%cell)
  deallocate(sol%val,sol%user,sol%name,sol%var_user,sol%name_user,sol%conserv_var)
  deallocate(L_str_criteria,L_var_criteria,L_eps)
  deallocate(gauss_point,gauss_weight)

contains

  subroutine buildMesh2(xL,xR,yL,yR,nx,ny,gauss_point,mesh)
    real(dp), intent(in) :: xL,xR,yL,yR
    integer, intent(in) :: nx,ny
    real(dp), dimension(:), intent(in) :: gauss_point
    type(meshStruct), intent(inout) :: mesh
    integer :: i,j,k,p,dir
    real(dp) :: dx,dy,a,b,c,center,diff,xc,yc
    type(edgeStruct) :: edge

    dx=(xR-xL)/nx
    dy=(yR-yL)/ny

    !Initialisation des points
    
    do j=0,ny
       do i=0,nx
          k=j*(nx+1)+i+1
          mesh%node(k)%x=xL+i*dx
          mesh%node(k)%y=yL+j*dy

          allocate(mesh%node(k)%connect(4))

          mesh%node(k)%connect(1)=j*(nx+1)+i
          mesh%node(k)%connect(2)=(j-1)*(nx+1)+i+1
          mesh%node(k)%connect(3)=j*(nx+1)+i+2
          mesh%node(k)%connect(4)=(j+1)*(nx+1)+i+1
          if (i==0) then
             mesh%node(k)%connect(1)=-1
          endif
          if (j==0) then
             mesh%node(k)%connect(2)=-1
          endif
          if (i==nx) then
             mesh%node(k)%connect(3)=-1
          endif
          if (j==ny) then
             mesh%node(k)%connect(4)=-1
          endif
       enddo
    enddo

!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
    allocate(mesh%node(17)%connect(3),mesh%node(18)%connect(3),mesh%node(19)%connect(4))
    allocate(mesh%node(20)%connect(3),mesh%node(21)%connect(3))
    mesh%node(17)%x=1.5_dp
    mesh%node(17)%y=1.0_dp
    mesh%node(18)%x=1.0_dp
    mesh%node(18)%y=1.5_dp
    mesh%node(19)%x=1.5_dp
    mesh%node(19)%y=1.5_dp
    mesh%node(20)%x=2.0_dp
    mesh%node(20)%y=1.5_dp
    mesh%node(21)%x=1.5_dp
    mesh%node(21)%y=2.0_dp
    
    mesh%node(17)%connect(1)=6
    mesh%node(17)%connect(2)=7
    mesh%node(17)%connect(3)=19
    mesh%node(18)%connect(1)=6
    mesh%node(18)%connect(2)=19
    mesh%node(18)%connect(3)=10
    mesh%node(19)%connect(1)=18
    mesh%node(19)%connect(2)=17
    mesh%node(19)%connect(3)=20
    mesh%node(19)%connect(4)=21
    mesh%node(20)%connect(1)=19
    mesh%node(20)%connect(2)=7
    mesh%node(20)%connect(3)=11
    mesh%node(21)%connect(1)=10
    mesh%node(21)%connect(2)=19
    mesh%node(21)%connect(3)=11
!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!    

    !Initialisation des cells
    
    do j=1,ny
       do i=1,nx
          k=(j-1)*nx+i
          mesh%cell(k)%dx=dx
          mesh%cell(k)%dy=dy
          mesh%cell(k)%xc=xL+i*dx-dx/2.0_dp
          mesh%cell(k)%yc=yL+j*dy-dy/2.0_dp

          allocate(mesh%cell(k)%edge(4),mesh%cell(k)%neigh(8))
          allocate(mesh%cell(k)%X_gauss(size(gauss_point)),mesh%cell(k)%Y_gauss(size(gauss_point)))

          mesh%cell(k)%corner(1)=(j-1)*(nx+1)+i
          mesh%cell(k)%corner(2)=(j-1)*(nx+1)+i+1
          mesh%cell(k)%corner(3)=j*(nx+1)+i+1
          mesh%cell(k)%corner(4)=j*(nx+1)+i

          mesh%cell(k)%edge(1)=k+j-1
          mesh%cell(k)%edge(2)=(nx+1)*ny+j+(i-1)*(ny+1)
          mesh%cell(k)%edge(3)=k+j
          mesh%cell(k)%edge(4)=(nx+1)*ny+j+(i-1)*(ny+1)+1
          
          mesh%cell(k)%neigh(1)=k-nx-1
          mesh%cell(k)%neigh(2)=k-1
          mesh%cell(k)%neigh(3)=k+nx-1
          mesh%cell(k)%neigh(4)=k+nx
          mesh%cell(k)%neigh(5)=k+nx+1
          mesh%cell(k)%neigh(6)=k+1
          mesh%cell(k)%neigh(7)=k-nx+1
          mesh%cell(k)%neigh(8)=k-nx

          xc=mesh%cell(k)%xc
          yc=mesh%cell(k)%yc
          do p=1,size(gauss_point)
             mesh%cell(k)%X_gauss(p)=xc+dx*gauss_point(p)/2.0_dp
             mesh%cell(k)%Y_gauss(p)=yc+dy*gauss_point(p)/2.0_dp
          enddo
          
          if (i==1) then
             mesh%cell(k)%neigh(1)=-((j-1)*nx)
             mesh%cell(k)%neigh(2)=-(j*nx)
             mesh%cell(k)%neigh(3)=-((j+1)*nx)
          endif

          if (j==1) then
             mesh%cell(k)%neigh(7)=-(i+(ny-1)*nx+1)
             mesh%cell(k)%neigh(8)=-(i+(ny-1)*nx)
             mesh%cell(k)%neigh(1)=-(i+(ny-1)*nx-1)
          endif
          
          if (i==nx) then
             mesh%cell(k)%neigh(5)=-(j*nx+1)
             mesh%cell(k)%neigh(6)=-((j-1)*nx+1)
             mesh%cell(k)%neigh(7)=-((j-2)*nx+1)
          endif

          if (j==ny) then
             mesh%cell(k)%neigh(3)=-(i-1)
             mesh%cell(k)%neigh(4)=-i
             mesh%cell(k)%neigh(5)=-(i+1)
          endif
       enddo
    enddo
    mesh%cell(1)%neigh(1)=-nx*ny
    mesh%cell(nx*ny)%neigh(5)=-1
    mesh%cell(nx)%neigh(7)=-((ny-1)*nx+1)
    mesh%cell((ny-1)*nx+1)%neigh(3)=-nx

!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
    mesh%cell(5)%dx=dx/2.0_dp
    mesh%cell(5)%dy=dy/2.0_dp
    mesh%cell(5)%xc=1.25_dp
    mesh%cell(5)%yc=1.25_dp
    mesh%cell(10)%dx=dx/2.0_dp
    mesh%cell(10)%dy=dy/2.0_dp
    mesh%cell(10)%xc=1.75_dp
    mesh%cell(10)%yc=1.25_dp
    mesh%cell(11)%dx=dx/2.0_dp
    mesh%cell(11)%dy=dy/2.0_dp
    mesh%cell(11)%xc=1.25_dp
    mesh%cell(11)%yc=1.75_dp
    mesh%cell(12)%dx=dx/2.0_dp
    mesh%cell(12)%dy=dy/2.0_dp
    mesh%cell(12)%xc=1.75_dp
    mesh%cell(12)%yc=1.75_dp

    allocate(mesh%cell(10)%edge(4),mesh%cell(11)%edge(4),mesh%cell(12)%edge(4))
    deallocate(mesh%cell(2)%edge,mesh%cell(4)%edge,mesh%cell(6)%edge,mesh%cell(8)%edge)
    allocate(mesh%cell(2)%edge(5),mesh%cell(4)%edge(5),mesh%cell(6)%edge(5),mesh%cell(8)%edge(5))
    deallocate(mesh%cell(2)%neigh,mesh%cell(4)%neigh,mesh%cell(5)%neigh,mesh%cell(6)%neigh,mesh%cell(8)%neigh)
    allocate(mesh%cell(2)%neigh(9),mesh%cell(4)%neigh(9),mesh%cell(5)%neigh(6),mesh%cell(6)%neigh(9),mesh%cell(8)%neigh(9))
    allocate(mesh%cell(10)%neigh(6),mesh%cell(11)%neigh(6),mesh%cell(12)%neigh(6))
    allocate(mesh%cell(10)%X_gauss(size(gauss_point)),mesh%cell(10)%Y_gauss(size(gauss_point)))
    allocate(mesh%cell(11)%X_gauss(size(gauss_point)),mesh%cell(11)%Y_gauss(size(gauss_point)))
    allocate(mesh%cell(12)%X_gauss(size(gauss_point)),mesh%cell(12)%Y_gauss(size(gauss_point)))

    mesh%cell(5)%corner(1)=6
    mesh%cell(5)%corner(2)=17
    mesh%cell(5)%corner(3)=19
    mesh%cell(5)%corner(4)=18
    mesh%cell(10)%corner(1)=17
    mesh%cell(10)%corner(2)=7
    mesh%cell(10)%corner(3)=20
    mesh%cell(10)%corner(4)=19
    mesh%cell(11)%corner(1)=18
    mesh%cell(11)%corner(2)=19
    mesh%cell(11)%corner(3)=21
    mesh%cell(11)%corner(4)=10
    mesh%cell(12)%corner(1)=19
    mesh%cell(12)%corner(2)=20
    mesh%cell(12)%corner(3)=11
    mesh%cell(12)%corner(4)=21

    mesh%cell(5)%edge(1)=6
    mesh%cell(5)%edge(2)=27
    mesh%cell(5)%edge(3)=7
    mesh%cell(5)%edge(4)=28
    mesh%cell(10)%edge(1)=7
    mesh%cell(10)%edge(2)=30
    mesh%cell(10)%edge(3)=18
    mesh%cell(10)%edge(4)=31
    mesh%cell(11)%edge(1)=19
    mesh%cell(11)%edge(2)=28
    mesh%cell(11)%edge(3)=25
    mesh%cell(11)%edge(4)=29
    mesh%cell(12)%edge(1)=25
    mesh%cell(12)%edge(2)=31
    mesh%cell(12)%edge(3)=26
    mesh%cell(12)%edge(4)=32

    mesh%cell(2)%edge(1)=2
    mesh%cell(2)%edge(2)=17
    mesh%cell(2)%edge(3)=3
    mesh%cell(2)%edge(4)=30
    mesh%cell(2)%edge(5)=27
    mesh%cell(4)%edge(1)=5
    mesh%cell(4)%edge(2)=14
    mesh%cell(4)%edge(3)=6
    mesh%cell(4)%edge(4)=19
    mesh%cell(4)%edge(5)=15
    mesh%cell(6)%edge(1)=26
    mesh%cell(6)%edge(2)=18
    mesh%cell(6)%edge(3)=22
    mesh%cell(6)%edge(4)=8
    mesh%cell(6)%edge(5)=23
    mesh%cell(8)%edge(1)=10
    mesh%cell(8)%edge(2)=29
    mesh%cell(8)%edge(3)=32
    mesh%cell(8)%edge(4)=11
    mesh%cell(8)%edge(5)=20

    mesh%cell(5)%neigh(1)=1
    mesh%cell(5)%neigh(2)=4
    mesh%cell(5)%neigh(3)=11
    mesh%cell(5)%neigh(4)=12
    mesh%cell(5)%neigh(5)=10
    mesh%cell(5)%neigh(6)=2
    mesh%cell(10)%neigh(1)=2
    mesh%cell(10)%neigh(2)=5
    mesh%cell(10)%neigh(3)=11
    mesh%cell(10)%neigh(4)=12
    mesh%cell(10)%neigh(5)=6
    mesh%cell(10)%neigh(6)=3
    mesh%cell(11)%neigh(1)=4
    mesh%cell(11)%neigh(2)=7
    mesh%cell(11)%neigh(3)=8
    mesh%cell(11)%neigh(4)=12
    mesh%cell(11)%neigh(5)=10
    mesh%cell(11)%neigh(6)=5
    mesh%cell(12)%neigh(1)=5
    mesh%cell(12)%neigh(2)=11
    mesh%cell(12)%neigh(3)=8
    mesh%cell(12)%neigh(4)=9
    mesh%cell(12)%neigh(5)=6
    mesh%cell(12)%neigh(6)=10

    mesh%cell(2)%neigh(1)=-7
    mesh%cell(2)%neigh(2)=1
    mesh%cell(2)%neigh(3)=4
    mesh%cell(2)%neigh(4)=5
    mesh%cell(2)%neigh(5)=10
    mesh%cell(2)%neigh(6)=6
    mesh%cell(2)%neigh(7)=3
    mesh%cell(2)%neigh(8)=-9
    mesh%cell(2)%neigh(9)=-8
    mesh%cell(4)%neigh(1)=-3
    mesh%cell(4)%neigh(2)=-6
    mesh%cell(4)%neigh(3)=-9
    mesh%cell(4)%neigh(4)=7
    mesh%cell(4)%neigh(5)=8
    mesh%cell(4)%neigh(6)=11
    mesh%cell(4)%neigh(7)=5
    mesh%cell(4)%neigh(8)=2
    mesh%cell(4)%neigh(9)=1
    mesh%cell(6)%neigh(1)=2
    mesh%cell(6)%neigh(2)=10
    mesh%cell(6)%neigh(3)=12
    mesh%cell(6)%neigh(4)=8
    mesh%cell(6)%neigh(5)=9
    mesh%cell(6)%neigh(6)=-7
    mesh%cell(6)%neigh(7)=-4
    mesh%cell(6)%neigh(8)=-1
    mesh%cell(6)%neigh(9)=3
    mesh%cell(8)%neigh(1)=4
    mesh%cell(8)%neigh(2)=7
    mesh%cell(8)%neigh(3)=-1
    mesh%cell(8)%neigh(4)=-2
    mesh%cell(8)%neigh(5)=-3
    mesh%cell(8)%neigh(6)=9
    mesh%cell(8)%neigh(7)=6
    mesh%cell(8)%neigh(8)=12
    mesh%cell(8)%neigh(9)=11
    
    do p=1,size(gauss_point)
       mesh%cell(5)%X_gauss(p)=mesh%cell(5)%xc+dx*gauss_point(p)/2.0_dp
       mesh%cell(5)%Y_gauss(p)=mesh%cell(5)%yc+dy*gauss_point(p)/2.0_dp
       mesh%cell(10)%X_gauss(p)=mesh%cell(10)%xc+dx*gauss_point(p)/2.0_dp
       mesh%cell(10)%Y_gauss(p)=mesh%cell(10)%yc+dy*gauss_point(p)/2.0_dp
       mesh%cell(11)%X_gauss(p)=mesh%cell(11)%xc+dx*gauss_point(p)/2.0_dp
       mesh%cell(11)%Y_gauss(p)=mesh%cell(11)%yc+dy*gauss_point(p)/2.0_dp
       mesh%cell(12)%X_gauss(p)=mesh%cell(12)%xc+dx*gauss_point(p)/2.0_dp
       mesh%cell(12)%Y_gauss(p)=mesh%cell(12)%yc+dy*gauss_point(p)/2.0_dp
    enddo
!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!

    !Initialisation des edges
    
    do j=1,ny
       do i=1,nx+1
          k=(j-1)*(nx+1)+i
          mesh%edge(k)%node1=k
          mesh%edge(k)%node2=k+nx+1
          mesh%edge(k)%cell1=(j-1)*nx+i-1
          mesh%edge(k)%cell2=(j-1)*nx+i
          mesh%edge(k)%dir=1
          mesh%edge(k)%length=dy
          mesh%edge(k)%lengthN=dx
          mesh%edge(k)%period=k
          mesh%edge(k)%sub=0
       enddo
       mesh%edge((j-1)*(nx+1)+1)%cell1=-j*nx
       mesh%edge((j-1)*(nx+1)+1)%period=(j-1)*(nx+1)+nx+1
       mesh%edge((j-1)*(nx+1)+nx+1)%cell2=-((j-1)*nx+1)
       mesh%edge((j-1)*(nx+1)+nx+1)%period=(j-1)*(nx+1)+1
    enddo
    do i=1,nx
       do j=1,ny+1
          k=(i-1)*(ny+1)+j+(nx+1)*ny
          mesh%edge(k)%node1=i+(j-1)*(nx+1)
          mesh%edge(k)%node2=i+(j-1)*(nx+1)+1
          mesh%edge(k)%cell1=i+(j-1)*nx-nx
          mesh%edge(k)%cell2=i+(j-1)*nx
          mesh%edge(k)%dir=2
          mesh%edge(k)%length=dx
          mesh%edge(k)%lengthN=dy
          mesh%edge(k)%period=k
          mesh%edge(k)%sub=0
       enddo
       mesh%edge((i-1)*(ny+1)+1+(nx+1)*ny)%cell1=-(nx*(ny-1)+i)
       mesh%edge((i-1)*(ny+1)+1+(nx+1)*ny)%period=(i-1)*(ny+1)+ny+1+(nx+1)*ny
       mesh%edge((i-1)*(ny+1)+ny+1+(nx+1)*ny)%cell2=-i
       mesh%edge((i-1)*(ny+1)+ny+1+(nx+1)*ny)%period=(i-1)*(ny+1)+1+(nx+1)*ny
    enddo

    do i=1,mesh%ne
       edge=mesh%edge(i)
       dir=edge%dir
       allocate(mesh%edge(i)%X_gauss(size(gauss_point)),mesh%edge(i)%Y_gauss(size(gauss_point)))
       allocate(mesh%edge(i)%flux_acc(size(gauss_point)))
       mesh%edge(i)%flux_acc=.false.
       select case (dir)
       case(1)
          a=mesh%node(edge%node1)%y
          b=mesh%node(edge%node2)%y
          c=mesh%node(edge%node1)%x
          center=(a+b)/2.0_dp
          diff=(b-a)/2.0_dp
          do p=1,size(gauss_point)
             mesh%edge(i)%X_gauss(p)=c
             mesh%edge(i)%Y_gauss(p)=center+diff*gauss_point(p)
          enddo
       case(2)
          a=mesh%node(edge%node1)%x
          b=mesh%node(edge%node2)%x
          c=mesh%node(edge%node1)%y
          center=(a+b)/2.0_dp
          diff=(b-a)/2.0_dp
          do p=1,size(gauss_point)
             mesh%edge(i)%X_gauss(p)=center+diff*gauss_point(p)
             mesh%edge(i)%Y_gauss(p)=c
          enddo
       end select
    enddo

!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
    mesh%edge(6)%node1=6
    mesh%edge(6)%node2=18
    mesh%edge(6)%cell1=4
    mesh%edge(6)%cell2=5
    mesh%edge(6)%dir=1
    mesh%edge(6)%length=dy
    mesh%edge(6)%lengthN=dx
    mesh%edge(6)%period=6
    mesh%edge(6)%sub(1)=1
    mesh%edge(6)%sub(2)=0
    
    mesh%edge(7)%node1=17
    mesh%edge(7)%node2=19
    mesh%edge(7)%cell1=5
    mesh%edge(7)%cell2=10
    mesh%edge(7)%dir=1
    mesh%edge(7)%length=dy
    mesh%edge(7)%lengthN=dx
    mesh%edge(7)%period=7
    mesh%edge(7)%sub(1)=0
    mesh%edge(7)%sub(2)=0

    mesh%edge(18)%node1=7
    mesh%edge(18)%node2=20
    mesh%edge(18)%cell1=10
    mesh%edge(18)%cell2=6
    mesh%edge(18)%dir=1
    mesh%edge(18)%length=dy
    mesh%edge(18)%lengthN=dx
    mesh%edge(18)%period=18
    mesh%edge(18)%sub(1)=0
    mesh%edge(18)%sub(2)=1

    mesh%edge(19)%node1=18
    mesh%edge(19)%node2=10
    mesh%edge(19)%cell1=4
    mesh%edge(19)%cell2=11
    mesh%edge(19)%dir=1
    mesh%edge(19)%length=dy
    mesh%edge(19)%lengthN=dx
    mesh%edge(19)%period=19
    mesh%edge(19)%sub(1)=1
    mesh%edge(19)%sub(2)=0

    mesh%edge(25)%node1=19
    mesh%edge(25)%node2=21
    mesh%edge(25)%cell1=11
    mesh%edge(25)%cell2=12
    mesh%edge(25)%dir=1
    mesh%edge(25)%length=dy
    mesh%edge(25)%lengthN=dx
    mesh%edge(25)%period=25
    mesh%edge(25)%sub(1)=0
    mesh%edge(25)%sub(2)=0

    mesh%edge(26)%node1=20
    mesh%edge(26)%node2=11
    mesh%edge(26)%cell1=12
    mesh%edge(26)%cell2=6
    mesh%edge(26)%dir=1
    mesh%edge(26)%length=dy
    mesh%edge(26)%lengthN=dx
    mesh%edge(26)%period=26
    mesh%edge(26)%sub(1)=0
    mesh%edge(26)%sub(2)=1

    mesh%edge(27)%node1=6
    mesh%edge(27)%node2=17
    mesh%edge(27)%cell1=2
    mesh%edge(27)%cell2=5
    mesh%edge(27)%dir=2
    mesh%edge(27)%length=dx
    mesh%edge(27)%lengthN=dy
    mesh%edge(27)%period=27
    mesh%edge(27)%sub(1)=1
    mesh%edge(27)%sub(2)=0

    mesh%edge(28)%node1=18
    mesh%edge(28)%node2=19
    mesh%edge(28)%cell1=5
    mesh%edge(28)%cell2=11
    mesh%edge(28)%dir=2
    mesh%edge(28)%length=dx
    mesh%edge(28)%lengthN=dy
    mesh%edge(28)%period=28
    mesh%edge(28)%sub(1)=0
    mesh%edge(28)%sub(2)=0

    mesh%edge(29)%node1=10
    mesh%edge(29)%node2=21
    mesh%edge(29)%cell1=11
    mesh%edge(29)%cell2=8
    mesh%edge(29)%dir=2
    mesh%edge(29)%length=dx
    mesh%edge(29)%lengthN=dy
    mesh%edge(29)%period=29
    mesh%edge(29)%sub(1)=0
    mesh%edge(29)%sub(2)=1

    mesh%edge(30)%node1=17
    mesh%edge(30)%node2=7
    mesh%edge(30)%cell1=2
    mesh%edge(30)%cell2=10
    mesh%edge(30)%dir=2
    mesh%edge(30)%length=dx
    mesh%edge(30)%lengthN=dy
    mesh%edge(30)%period=30
    mesh%edge(30)%sub(1)=1
    mesh%edge(30)%sub(2)=0

    mesh%edge(31)%node1=19
    mesh%edge(31)%node2=20
    mesh%edge(31)%cell1=10
    mesh%edge(31)%cell2=12
    mesh%edge(31)%dir=2
    mesh%edge(31)%length=dx
    mesh%edge(31)%lengthN=dy
    mesh%edge(31)%period=31
    mesh%edge(31)%sub(1)=0
    mesh%edge(31)%sub(2)=0

    mesh%edge(32)%node1=21
    mesh%edge(32)%node2=11
    mesh%edge(32)%cell1=12
    mesh%edge(32)%cell2=8
    mesh%edge(32)%dir=2
    mesh%edge(32)%length=dx
    mesh%edge(32)%lengthN=dy
    mesh%edge(32)%period=32
    mesh%edge(32)%sub(1)=0
    mesh%edge(32)%sub(2)=1

    do p=1,size(gauss_point)
       mesh%edge(6)%X_gauss(p)=1.0_dp
       mesh%edge(6)%Y_gauss(p)=1.25_dp+dy/4.0_dp*gauss_point(p)
       mesh%edge(7)%X_gauss(p)=1.5_dp
       mesh%edge(7)%Y_gauss(p)=1.25_dp+dy/4.0_dp*gauss_point(p)
       mesh%edge(18)%X_gauss(p)=2.0_dp
       mesh%edge(18)%Y_gauss(p)=1.25_dp+dy/4.0_dp*gauss_point(p)
       mesh%edge(19)%X_gauss(p)=1.0_dp
       mesh%edge(19)%Y_gauss(p)=1.75_dp+dy/4.0_dp*gauss_point(p)
       mesh%edge(25)%X_gauss(p)=1.5_dp
       mesh%edge(25)%Y_gauss(p)=1.75_dp+dy/4.0_dp*gauss_point(p)
       mesh%edge(26)%X_gauss(p)=2.0_dp
       mesh%edge(26)%Y_gauss(p)=1.75_dp+dy/4.0_dp*gauss_point(p)
       mesh%edge(27)%X_gauss(p)=1.25_dp+dx/4.0_dp*gauss_point(p)
       mesh%edge(27)%Y_gauss(p)=1.0_dp
       mesh%edge(28)%X_gauss(p)=1.25_dp+dx/4.0_dp*gauss_point(p)
       mesh%edge(28)%Y_gauss(p)=1.5_dp
       mesh%edge(29)%X_gauss(p)=1.25_dp+dx/4.0_dp*gauss_point(p)
       mesh%edge(29)%Y_gauss(p)=2.0_dp
       mesh%edge(30)%X_gauss(p)=1.75_dp+dx/4.0_dp*gauss_point(p)
       mesh%edge(30)%Y_gauss(p)=1.0_dp
       mesh%edge(31)%X_gauss(p)=1.75_dp+dx/4.0_dp*gauss_point(p)
       mesh%edge(31)%Y_gauss(p)=1.5_dp
       mesh%edge(32)%X_gauss(p)=1.75_dp+dx/4.0_dp*gauss_point(p)
       mesh%edge(32)%Y_gauss(p)=2.0_dp
    enddo
    

!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
    
    return
  end subroutine buildMesh2

  subroutine init_FV(test_case,str_equa,str_flux,str_time_scheme,IC_func,BC,exactSol, &
       f_equa,flux,speed,time_scheme,sol)
    character(len=20), intent(in) :: test_case,str_equa,str_flux,str_time_scheme
    procedure (sub_IC), pointer, intent(out) :: IC_func
    procedure (sub_BC), pointer, intent(out) :: BC
    procedure (sub_exactsol), pointer, intent(out) :: exactSol
    procedure (sub_f), pointer, intent(out) :: f_equa
    procedure (sub_flux), pointer, intent(out) :: flux
    procedure (sub_speed), pointer, intent(out) :: speed
    procedure (sub_time), pointer, intent(out) :: time_scheme
    type(solStruct), intent(inout) :: sol
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
    case ('Test')
       IC_func => IC_func_test
       BC => BC_test
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

    return
  end subroutine init_FV

  subroutine calculation(mesh,sol,order,cfl,tf,fs,namefile,str_equa, &
       f_equa,flux,speed,time_scheme,exactSol, &
       L_str_criteria,L_var_criteria,L_eps,gauss_weight)
    type(meshStruct), intent(inout) :: mesh
    type(solStruct), intent(inout) :: sol
    integer, intent(in) :: order
    real(dp), intent(in) :: cfl,tf
    integer, intent(in) :: fs
    character(len=20),intent(in) :: namefile,str_equa
    procedure (sub_f), pointer, intent(in) :: f_equa
    procedure (sub_flux), pointer, intent(in) :: flux
    procedure (sub_speed), pointer, intent(in) :: speed
    procedure (sub_time), pointer, intent(in) :: time_scheme
    procedure (sub_exactsol), pointer, intent(in) :: exactSol
    character(len=20), dimension(:), intent(in) :: L_str_criteria
    integer, dimension(:), intent(in) :: L_var_criteria
    real(dp), dimension(:), intent(in) :: L_eps
    real(dp), dimension(:), intent(in) :: gauss_weight
    integer :: n
    real(dp) :: t
    
    t=0.0_dp
    n=1
    call check_conservativity(mesh,sol)
    do while (t<tf)
       call time_scheme(mesh,sol,str_equa,f_equa,flux,speed,order,cfl,t,n,tf, &
            L_str_criteria,L_var_criteria,L_eps,gauss_weight)
       if (mod(n,fs)==0.or.t>=tf) then
          call userSol(t,mesh,sol,str_equa,exactSol,gauss_weight)
          call writeSol(mesh,sol,namefile,n/fs)
          call print(mesh,sol,t,n)
       endif
       n=n+1
    enddo

  end subroutine calculation
    
end program test
