module inout

  use constant
  use types
  use efficiency
  use phys
  
  implicit none

contains

  subroutine get_config(config_file)
    character(len=20), intent(out) :: config_file
    
    print*,"Configuration file : "
    read*,config_file

    return
  end subroutine get_config
  
  subroutine init(config_file,test_case,xL,xR,yL,yR,nx,ny,nvar,cfl,tf,fs,namefile,mesh,sol, &
       str_equa,str_flux,str_time_scheme,order,L_str_criteria,L_var_criteria,L_eps,gauss_point,gauss_weight)
    character(len=20), intent(out) :: config_file,test_case,namefile,str_equa,str_flux,str_time_scheme
    real(dp), intent(out) :: xL,xR,yL,yR,cfl,tf
    integer, intent(out) :: nx,ny,nvar,fs,order
    type(meshStruct), intent(out) :: mesh
    type(solStruct), intent(out) :: sol
    character(len=20), dimension(:), allocatable, intent(out) :: L_str_criteria
    integer, dimension(:), allocatable, intent(out) :: L_var_criteria
    real(dp), dimension(:), allocatable, intent(out) :: L_eps
    real(dp), dimension(:), allocatable, intent(out) :: gauss_point,gauss_weight
    integer :: i,ncriteria

    config_file="config/"//trim(config_file)
    open(11,file=config_file,form="formatted")
    read(11,*)test_case
    read(11,*)xL
    read(11,*)xR
    read(11,*)yL
    read(11,*)yR
    read(11,*)nx
    read(11,*)ny
    read(11,*)cfl
    read(11,*)tf
    read(11,*)str_equa
    read(11,*)str_flux
    read(11,*)str_time_scheme
    read(11,*)order
    read(11,*)ncriteria
    allocate(L_str_criteria(ncriteria),L_var_criteria(ncriteria),L_eps(ncriteria))
    do i=1,ncriteria
       read(11,*)L_str_criteria(i),L_var_criteria(i),L_eps(i)
    enddo
    read(11,*)fs
    read(11,*)namefile
    read(11,*)nvar
    allocate(sol%name(nvar))
    do i=1,nvar
       read(11,*)sol%name(i)
    enddo
    read(11,*)sol%nsolUser
    allocate(sol%var_user(sol%nsolUser),sol%name_user(sol%nsolUser))
    do i=1,sol%nsolUser
       read(11,*)sol%var_user(i)
    enddo
    close(11)
    
    mesh%nc=nx*ny
    mesh%ne=ny*(nx+1)+nx*(ny+1)
    mesh%np=(nx+1)*(ny+1)
    sol%nvar=nvar
    
    allocate(mesh%node(mesh%np),mesh%edge(mesh%ne),mesh%cell(mesh%nc))
    allocate(sol%val(mesh%nc,sol%nvar),sol%user(mesh%nc,sol%nsolUser),sol%conserv_var(sol%nvar,2))
    do i=1,mesh%ne
       allocate(mesh%edge(i)%flux(order,nvar))
    enddo

    allocate(gauss_point(order),gauss_weight(order))
    select case (order)
    case (1)
       gauss_point=gauss_point1
       gauss_weight=gauss_weight1
    case(2)
       gauss_point=gauss_point2
       gauss_weight=gauss_weight2
    case(3)
       gauss_point=gauss_point3
       gauss_weight=gauss_weight3
    case(4)
       gauss_point=gauss_point4
       gauss_weight=gauss_weight4
    case(5)
       gauss_point=gauss_point5
       gauss_weight=gauss_weight5
    case(6)
       gauss_point=gauss_point5
       gauss_weight=gauss_weight5
    case default
       print*,"Space order too high, no good enough quadrature implemented"
       call exit()
    end select

    return
  end subroutine init
  
  subroutine IC(IC_func,mesh,sol,order,gauss_point,gauss_weight)
    procedure (sub_IC), pointer, intent(in) :: IC_func
    type(meshStruct), intent(in) :: mesh
    type(solStruct), intent(inout) :: sol
    integer, intent(in) :: order
    real(dp), dimension(:), intent(in) :: gauss_point,gauss_weight
    real(dp), dimension(:), allocatable :: U0,S
    integer :: k,p1,p2
    real(dp) :: Xg,Yg

    allocate(U0(sol%nvar),S(sol%nvar))
    
    do k=1,mesh%nc
       U0=0.0_dp
       do p1=1,6
          do p2=1,6
             Xg=mesh%cell(k)%xc+gauss_point(p1)*mesh%cell(k)%dx/2.0_dp
             Yg=mesh%cell(k)%yc+gauss_point(p2)*mesh%cell(k)%dy/2.0_dp
             call IC_func(Xg,Yg,S)
             U0=U0+S*gauss_weight(p1)*gauss_weight(p2)/4.0_dp
          enddo
       enddo
       sol%val(k,:)=U0
    enddo

    if (order>6) then
       print*,"The order of initial condition's quadrature is to low"
       call exit()
    endif

    deallocate(U0,S)
    
    return
  end subroutine IC
  
  subroutine buildMesh(xL,xR,yL,yR,nx,ny,gauss_point,mesh)
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
    
    return
  end subroutine buildMesh
  
  subroutine writeSol(mesh,sol,namefile,nfile)
    type(meshStruct), intent(in) :: mesh
    type(solStruct), intent(in) :: sol
    character(len=20), intent(in) :: namefile
    integer, intent(in) :: nfile
    integer :: k,n
    integer :: i1,i2,i3,i4
    character(len=34) :: completenamefile

    write(completenamefile,'(A,A,I3.3,A)')'./results/',trim(namefile),nfile,'.vtk'
    open(11,file=completenamefile,form="formatted")
    
    write(11,'(a)')"# vtk DataFile Version 2.0"
    write(11,'(a)')"Results of the calculation"
    write(11,'(a)')"ASCII"
    write(11,'(a)')"DATASET UNSTRUCTURED_GRID"
    
    write(11,'(a,i8,a)')"POINTS ",mesh%np," float"
    do k=1,mesh%np
       write(11,'(e15.8,a,e15.8,a,e15.8)')mesh%node(k)%x," ",mesh%node(k)%y," ",0.0_dp
    enddo

    write(11,'(a,i8,i9)')"CELLS ",mesh%nc,5*mesh%nc
    do k=1,mesh%nc
       i1=mesh%cell(k)%corner(1)-1
       i2=mesh%cell(k)%corner(2)-1
       i3=mesh%cell(k)%corner(3)-1
       i4=mesh%cell(k)%corner(4)-1
       write(11,'(i1,i8,i8,i8,i8)')4,i1,i2,i3,i4
    enddo

    write(11,'(a,i8)')"CELL_TYPES ",mesh%nc
    do k=1,mesh%nc
       write(11,'(i1)')9
    enddo
    
    write(11,'(a,i8)')"CELL_DATA ",mesh%nc
    do n=1,sol%nvar
       write(11,'(a,a,a)')"SCALARS ",sol%name(n)," float 1"
       write(11,'(a)')"LOOKUP_TABLE default"
       do k=1,mesh%nc
          if (abs(sol%val(k,n))>10.0**(-dp)) then
             write(11,'(e15.8)')sol%val(k,n)
          else
             write(11,'(e15.8)')0.0_dp
          endif
       enddo
    enddo
    
    do n=1,sol%nsolUser
       write(11,'(a,a,a)')"SCALARS ",sol%name_user(n)," float 1"
       write(11,'(a)')"LOOKUP_TABLE default"
       do k=1,mesh%nc
             write(11,'(e15.8)')sol%user(k,n)
       enddo
    enddo
    
    close(11)
    
    return
  end subroutine writeSol

  subroutine print(mesh,sol,t,n)
    type(meshStruct), intent(in) :: mesh
    type(solStruct), intent(in) :: sol
    real(dp), intent(in) :: t
    integer, intent(in) :: n

    print*,"t=",t,"itération ",n
    call check_conservativity(mesh,sol)
    print*,"-----------------------------------------"
    
    return
  end subroutine print

  subroutine write_accept(mesh,NOT_ACCEPTED_CELL,n,count)
    type(meshStruct), intent(in) :: mesh
    integer, dimension(:), intent(in) :: NOT_ACCEPTED_CELL
    integer, intent(in) :: n,count
    integer :: k
    integer :: i1,i2,i3,i4
    character(len=34) :: completenamefile

    write(completenamefile,'(A,A,I3.3,I3.3,A)')'./results/','accept',n,count,'.vtk'
    open(11,file=completenamefile,form="formatted")
    
    write(11,'(a)')"# vtk DataFile Version 2.0"
    write(11,'(a)')"Results of the calculation"
    write(11,'(a)')"ASCII"
    write(11,'(a)')"DATASET UNSTRUCTURED_GRID"
    
    write(11,'(a,i8,a)')"POINTS ",mesh%np," float"
    do k=1,mesh%np
       write(11,'(e15.8,a,e15.8,a,e15.8)')mesh%node(k)%x," ",mesh%node(k)%y," ",0.0_dp
    enddo

    write(11,'(a,i8,i9)')"CELLS ",mesh%nc,5*mesh%nc
    do k=1,mesh%nc
       i1=mesh%edge(mesh%cell(k)%edge(1))%node1-1
       i2=mesh%edge(mesh%cell(k)%edge(3))%node1-1
       i3=mesh%edge(mesh%cell(k)%edge(3))%node2-1
       i4=mesh%edge(mesh%cell(k)%edge(1))%node2-1
       write(11,'(i1,i8,i8,i8,i8)')4,i1,i2,i3,i4
    enddo

    write(11,'(a,i8)')"CELL_TYPES ",mesh%nc
    do k=1,mesh%nc
       write(11,'(i1)')9
    enddo
    
    write(11,'(a,i8)')"CELL_DATA ",mesh%nc
    write(11,'(a,a,a)')"SCALARS ","accept"," float 1"
    write(11,'(a)')"LOOKUP_TABLE default"
    do k=1,mesh%nc
       if (all(k/=NOT_ACCEPTED_CELL)) then
          write(11,'(e15.8)')0.0_dp
       else
          write(11,'(e15.8)')1.0_dp
       endif
    enddo
    
    close(11)
    
    return
  end subroutine write_accept

end module inout
