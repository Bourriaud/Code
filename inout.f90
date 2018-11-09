module inout

  use constant
  use types
  use efficiency
  
  implicit none

contains

  subroutine init(xL,xR,yL,yR,nx,ny,nvar,cfl,tf,fs,namefile,mesh,sol,str_equa,str_flux,str_time_scheme,order)
    real(dp), intent(out) :: xL,xR,yL,yR,cfl,tf
    integer, intent(out) :: nx,ny,nvar,fs,order
    character(len=20), intent(out) :: namefile,str_equa,str_flux,str_time_scheme
    type(meshStruct), intent(out) :: mesh
    type(solStruct), intent(out) :: sol
    integer :: i

    open(11,file="configuration",form="formatted")
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
    read(11,*)fs
    read(11,*)namefile
    read(11,*)nvar
    allocate(sol%name(nvar))
    do i=1,nvar
       read(11,*)sol%name(i)
    enddo
    read(11,*)sol%nsolUser
    allocate(sol%nameUser(sol%nsolUser))
    do i=1,sol%nsolUser
       read(11,*)sol%nameUser(i)
    enddo
    close(11)
    
    mesh%nc=nx*ny
    mesh%ne=nx*(nx+1)+ny*(ny+1)
    mesh%np=(nx+1)*(ny+1)
    sol%nvar=nvar
    
    allocate(mesh%node(mesh%np),mesh%edge(mesh%ne),mesh%cell(mesh%nc))
    allocate(sol%val(mesh%nc,sol%nvar),sol%user(mesh%nc,sol%nsolUser))

    return
  end subroutine init

  subroutine IC_func(x,y,t,s)
    real(dp), intent(in) :: x,y,t
    real(dp), intent(out) :: s

    s=cos((x-5.0_dp)*pi/5.0_dp)+cos((y-5.0_dp)*pi/5.0_dp)
    
    return
  end subroutine IC_func
  
  subroutine IC(mesh,sol,quad_t)
    type(meshStruct), intent(in) :: mesh
    type(solStruct), intent(inout) :: sol
    procedure (quadrature_t), pointer, intent(in) :: quad_t
    integer :: k
    real(dp) :: x,y,dx,dy,rho
    !real(dp) :: u,v,p,gamma

    !gamma=1.4_dp
    do k=1,mesh%nc
       x=mesh%cell(k)%xc
       y=mesh%cell(k)%yc
       dx=mesh%cell(k)%dx
       dy=mesh%cell(k)%dy
       call quad_t(IC_func,mesh,0.0_dp,k,rho)
       rho=rho/(dx*dy)
       sol%val(k,1)=rho
       sol%val(k,2)=rho
       !sol%val(k,3)=rho*v
       !sol%val(k,4)=rho*0.5*(u**2+v**2)+p/(gamma-1)
    enddo
    
    return
  end subroutine IC

  subroutine BC(nx,ny,nvar,mesh)
    integer, intent(in) :: nx,ny,nvar
    type(meshStruct), intent(inout) :: mesh
    integer :: i,j

    do i=1,mesh%ne
       mesh%edge(i)%boundType='NOT A BOUNDARY'
       allocate(mesh%edge(i)%bound(nvar))
    enddo
    
    do j=1,ny
       mesh%edge((j-1)*(nx+1)+1)%boundType='PERIODIC'
       mesh%edge((j-1)*(nx+1)+1)%bound(:)=0.0_dp
       
       mesh%edge(j*(nx+1))%boundType='PERIODIC'
       mesh%edge(j*(nx+1))%bound(:)=0.0_dp
    enddo
    
    do i=1,nx
       mesh%edge((i-1)*(ny+1)+1+(nx+1)*ny)%boundType='PERIODIC'
       mesh%edge((i-1)*(ny+1)+1+(nx+1)*ny)%bound(:)=0.0_dp
       
       mesh%edge(i*(ny+1)+(nx+1)*ny)%boundType='PERIODIC'
       mesh%edge(i*(ny+1)+(nx+1)*ny)%bound(:)=0.0_dp
    enddo
    return   
  end subroutine BC
  
  subroutine buildMesh(xL,xR,yL,yR,nx,ny,mesh)
    real(dp), intent(in) :: xL,xR,yL,yR
    integer, intent(in) :: nx,ny
    type(meshStruct), intent(inout) :: mesh
    integer :: i,j,k
    real(dp) :: dx,dy

    dx=(xR-xL)/nx
    dy=(yR-yL)/ny

    !Initialisation des points
    
    do j=0,ny
       do i=0,nx
          k=j*(nx+1)+i+1
          mesh%node(k)%x=xL+i*dx
          mesh%node(k)%y=yL+j*dy

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

          mesh%cell(k)%neigh(1)=k-nx-1
          mesh%cell(k)%neigh(2)=k-1
          mesh%cell(k)%neigh(3)=k+nx-1
          mesh%cell(k)%neigh(4)=k+nx
          mesh%cell(k)%neigh(5)=k+nx+1
          mesh%cell(k)%neigh(6)=k+1
          mesh%cell(k)%neigh(7)=k-nx+1
          mesh%cell(k)%neigh(8)=k-nx

          mesh%cell(k)%edge(1)%node1=(j-1)*(nx+1)+i
          mesh%cell(k)%edge(1)%node2=j*(nx+1)+i
          if (i==1) then
             mesh%cell(k)%neigh(1)=-1
             mesh%cell(k)%neigh(2)=-1
             mesh%cell(k)%neigh(3)=-1
          endif

          mesh%cell(k)%edge(2)%node1=(j-1)*(nx+1)+i
          mesh%cell(k)%edge(2)%node2=(j-1)*(nx+1)+i+1
          if (j==1) then
             mesh%cell(k)%neigh(7)=-1
             mesh%cell(k)%neigh(8)=-1
             mesh%cell(k)%neigh(1)=-1
          endif
          
          mesh%cell(k)%edge(3)%node1=(j-1)*(nx+1)+i+1
          mesh%cell(k)%edge(3)%node2=j*(nx+1)+i+1
          if (i==nx) then
             mesh%cell(k)%neigh(5)=-1
             mesh%cell(k)%neigh(6)=-1
             mesh%cell(k)%neigh(7)=-1
          endif

          mesh%cell(k)%edge(4)%node1=j*(nx+1)+i
          mesh%cell(k)%edge(4)%node2=j*(nx+1)+i+1
          if (j==ny) then
             mesh%cell(k)%neigh(3)=-1
             mesh%cell(k)%neigh(4)=-1
             mesh%cell(k)%neigh(5)=-1
          endif
       enddo
    enddo

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
       enddo
       mesh%edge((j-1)*(nx+1)+1)%cell1=-j*nx
       mesh%edge((j-1)*(nx+1)+nx+1)%cell2=-((j-1)*nx+1)
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
       enddo
       mesh%edge((i-1)*(ny+1)+1+(nx+1)*ny)%cell1=-(nx*(ny-1)+i)
       mesh%edge((i-1)*(ny+1)+ny+1+(nx+1)*ny)%cell2=-i
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
       i1=mesh%cell(k)%edge(1)%node1-1
       i2=mesh%cell(k)%edge(3)%node1-1
       i3=mesh%cell(k)%edge(3)%node2-1
       i4=mesh%cell(k)%edge(1)%node2-1
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
          write(11,'(e15.8)')sol%val(k,n)
       enddo
    enddo
    
    do n=1,sol%nsolUser
       write(11,'(a,a,a)')"SCALARS ",sol%nameUser(n)," float 1"
       write(11,'(a)')"LOOKUP_TABLE default"
       do k=1,mesh%nc
             write(11,'(e15.8)')sol%user(k,n)
       enddo
    enddo
    
    close(11)
    
    return
  end subroutine writeSol

  subroutine print(t,n)
    real(dp), intent(in) :: t
    integer, intent(in) :: n

    print*,"t=",t,"itération ",n
    
    return
  end subroutine print

end module inout
