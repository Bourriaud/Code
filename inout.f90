module inout

  use types
  
  implicit none

contains

  subroutine init(xL,xR,yL,yR,dx,dy,dt,tf,fs,namefile,mesh,sol)
    real, intent(out) :: xL,xR,yL,yR,dx,dy,dt,tf
    character(len=20), intent(out) :: namefile
    type(meshStruct), intent(out) :: mesh
    type(solStruct), intent(out) :: sol
    integer, intent(out) :: fs
    integer :: i

    open(11,file="configuration",form="formatted")
    read(11,*)xL
    read(11,*)xR
    read(11,*)yL
    read(11,*)yR
    read(11,*)mesh%nx
    read(11,*)mesh%ny
    read(11,*)dt
    read(11,*)tf
    read(11,*)fs
    read(11,*)namefile
    read(11,*)sol%nvar
    allocate(sol%name(sol%nvar))
    do i=1,sol%nvar
       read(11,*)sol%name(i)
    enddo
    read(11,*)sol%nsolUser
    allocate(sol%nameUser(sol%nsolUser))
    do i=1,sol%nsolUser
       read(11,*)sol%nameUser(i)
    enddo
    close(11)

    dx=(xR-xL)/mesh%nx
    dy=(yR-yL)/mesh%ny
    allocate(mesh%X(0:mesh%nx,0:mesh%ny),mesh%Y(0:mesh%nx,0:mesh%ny))
    allocate(mesh%CX(mesh%nx,mesh%ny),mesh%CY(mesh%nx,mesh%ny))
    allocate(mesh%boundType(mesh%nx,mesh%ny),mesh%bound(mesh%nx,mesh%ny,sol%nvar))
    allocate(sol%val(mesh%nx,mesh%ny,sol%nvar),sol%user(mesh%nx,mesh%ny,sol%nsolUser))

    return
  end subroutine init

  subroutine IC(mesh,sol)
    type(meshStruct), intent(in) :: mesh
    type(solStruct), intent(inout) :: sol
    integer :: i,j,n

    do i=1,mesh%nx
       do j=1,mesh%ny
          if ((mesh%CX(i,j)**2+mesh%CY(i,j)**2)**0.5<5) then
             sol%val(i,j,1)=1.
             sol%val(i,j,2)=1.
          else
             sol%val(i,j,1)=0.
             sol%val(i,j,2)=0.
          endif
       enddo
    enddo
    return
  end subroutine IC

  subroutine BC(xL,xR,yL,yR,mesh)
    real, intent(in) :: xL,xR,yL,yR
    type(meshStruct), intent(inout) :: mesh
    integer :: i,j

    mesh%boundType(:,:)='NOT A BOUNDARY'
    mesh%boundType(1,:)='DIRICHLET'
    mesh%boundType(mesh%nx,:)='DIRICHLET'
    mesh%boundType(:,1)='DIRICHLET'
    mesh%boundType(:,mesh%ny)='DIRICHLET'
    mesh%bound(:,:,:)=0.
    
    return   
  end subroutine BC
  
  subroutine buildMesh(xL,xR,yL,yR,mesh)
    real, intent(in) :: xL,xR,yL,yR
    type(meshStruct), intent(inout) :: mesh
    integer :: i
    real :: dx,dy

    dx=(xR-xL)/mesh%nx
    dy=(yR-yL)/mesh%ny

    do i=1,mesh%nx
       mesh%X(i,:)=xL+i*dx
       mesh%CX(i,:)=xL+i*dx-dx/2.
    enddo
    do i=1,mesh%ny
       mesh%Y(:,i)=yL+i*dy
       mesh%CY(:,i)=yL+i*dy-dy/2
    enddo
    mesh%X(0,:)=0.
    mesh%Y(:,0)=0.
    
    return
  end subroutine buildMesh
  
  subroutine writeSol(mesh,sol,namefile,nfile)
    type(meshStruct), intent(in) :: mesh
    type(solStruct), intent(in) :: sol
    character(len=20), intent(in) :: namefile
    integer, intent(in) :: nfile
    integer :: i,j,n
    integer :: nx,ny
    character(len=34) :: completenamefile

    write(completenamefile,'(A,A,I3.3,A)')'./results/',trim(namefile),nfile,'.vtk'
    open(11,file=completenamefile,form="formatted")
    write(11,'(a)')"# vtk DataFile Version 2.0"
    write(11,'(a)')"Results of the calculation"
    write(11,'(a)')"ASCII"
    write(11,'(a)')"DATASET STRUCTURED_GRID"
    write(11,'(a,i4,a,i4,a)')"DIMENSIONS ",mesh%nx+1," ",mesh%ny+1," 1"
    write(11,'(a,i8,a)')"POINTS ",(mesh%nx+1)*(mesh%ny+1)," float"
    do j=0,mesh%ny
       do i=0,mesh%nx
          write(11,'(e15.8,a,e15.8,a,e15.8)')mesh%X(i,j)," ",mesh%Y(i,j)," ",0.
       enddo
    enddo

    write(11,'(a,i8)')"CELL_DATA ",(mesh%nx)*(mesh%ny)
    do n=1,sol%nvar
       write(11,'(a,a,a)')"SCALARS ",sol%name(n)," float 1"
       write(11,'(a)')"LOOKUP_TABLE default"
       do j=1,mesh%ny
          do i=1,mesh%nx
             write(11,'(e15.8)')sol%val(i,j,n)
          enddo
       enddo
    enddo
    
    do n=1,sol%nsolUser
       write(11,'(a,a,a)')"SCALARS ",sol%nameUser(n)," float 1"
       write(11,'(a)')"LOOKUP_TABLE default"
       do j=1,mesh%ny
          do i=1,mesh%nx
             write(11,'(e15.8)')sol%user(i,j,n)
          enddo
       enddo
    enddo
    close(11)
    
    return
  end subroutine writeSol

  subroutine print(t,n)
    real, intent(in) :: t
    integer, intent(in) :: n

    print*,"t=",t,"itération ",n
    
    return
  end subroutine print

end module inout
