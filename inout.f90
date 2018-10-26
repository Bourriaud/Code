module inout

  use constant
  use types
  
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
    mesh%np=(nx+1)*(ny+1)
    sol%nvar=nvar
    
    allocate(mesh%node(mesh%np),mesh%cell(mesh%nc))
    allocate(sol%val(mesh%nc,sol%nvar),sol%user(mesh%nc,sol%nsolUser))

    return
  end subroutine init

  subroutine IC(mesh,sol)
    type(meshStruct), intent(in) :: mesh
    type(solStruct), intent(inout) :: sol
    integer :: k
    real(dp) :: x,y,rho
    !real(dp) :: u,v,p,gamma

    !gamma=1.4_dp
    do k=1,mesh%nc
       x=mesh%cell(k)%xc
       y=mesh%cell(k)%yc
       rho=cos((x-5.0_dp)*pi/5.0_dp)+cos((y-5.0_dp)*pi/5.0_dp)
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
    integer :: i,j,k

    do k=1,mesh%nc
       do i=1,size(mesh%cell(k)%edge)
          mesh%cell(k)%edge(i)%boundType='NOT A BOUNDARY'
          allocate(mesh%cell(k)%edge(i)%bound(nvar))
       enddo
    enddo
    
    do j=1,ny
       mesh%cell((j-1)*nx+1)%edge(1)%boundType='PERIODIC'
       mesh%cell((j-1)*nx+1)%edge(1)%bound(:)=0.0_dp
       
       mesh%cell((j-1)*nx+nx)%edge(3)%boundType='PERIODIC'
       mesh%cell((j-1)*nx+nx)%edge(3)%bound(:)=0.0_dp
    enddo
    
    do i=1,nx
       mesh%cell(i)%edge(2)%boundType='PERIODIC'
       mesh%cell(i)%edge(2)%bound(:)=0.0_dp
       
       mesh%cell((ny-1)*nx+i)%edge(4)%boundType='PERIODIC'
       mesh%cell((ny-1)*nx+i)%edge(4)%bound(:)=0.0_dp
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

    do j=1,ny
       do i=1,nx
          k=(j-1)*nx+i
          mesh%cell(k)%dx=dx
          mesh%cell(k)%dy=dy
          mesh%cell(k)%xc=xL+i*dx-dx/2.0_dp
          mesh%cell(k)%yc=yL+j*dy-dy/2.0_dp

          allocate(mesh%cell(k)%edge(4))

          mesh%cell(k)%edge(1)%node1=(j-1)*(nx+1)+i
          mesh%cell(k)%edge(1)%node2=j*(nx+1)+i
          mesh%cell(k)%edge(1)%normal=1
          mesh%cell(k)%edge(1)%neigh=k-1
          if (i==1) then
             mesh%cell(k)%edge(1)%neigh=-(j*nx)
          endif

          mesh%cell(k)%edge(2)%node1=(j-1)*(nx+1)+i
          mesh%cell(k)%edge(2)%node2=(j-1)*(nx+1)+i+1
          mesh%cell(k)%edge(2)%normal=2
          mesh%cell(k)%edge(2)%neigh=k-nx
          if (j==1) then
             mesh%cell(k)%edge(2)%neigh=-(nx*(ny-1)+i)
          endif
          
          mesh%cell(k)%edge(3)%node1=(j-1)*(nx+1)+i+1
          mesh%cell(k)%edge(3)%node2=j*(nx+1)+i+1
          mesh%cell(k)%edge(3)%normal=3
          mesh%cell(k)%edge(3)%neigh=k+1
          if (i==nx) then
             mesh%cell(k)%edge(3)%neigh=-((j-1)*nx+1)
          endif

          mesh%cell(k)%edge(4)%node1=j*(nx+1)+i
          mesh%cell(k)%edge(4)%node2=j*(nx+1)+i+1
          mesh%cell(k)%edge(4)%normal=4
          mesh%cell(k)%edge(4)%neigh=k+nx
          if (j==ny) then
             mesh%cell(k)%edge(4)%neigh=-i
          endif
       enddo
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
