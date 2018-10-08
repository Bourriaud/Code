module inout

  use types
  
  implicit none

contains

  subroutine init(xL,xR,yL,yR,dx,dy,dt,tf,namefile,mesh,sol)
    real, intent(out) :: xL,xR,yL,yR,dx,dy,dt,tf
    character(len=20), intent(out) :: namefile
    type(meshStruct), intent(out) :: mesh
    type(solStruct), intent(out) :: sol
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
    read(11,*)namefile
    read(11,*)sol%nvar
    allocate(sol%name(sol%nvar))
    do i=1,sol%nvar
       read(11,*)sol%name(i)
    enddo
    close(11)

    dx=(xR-xL)/mesh%nx
    dy=(yR-yL)/mesh%ny
    allocate(mesh%X(mesh%nx,mesh%ny),mesh%Y(mesh%nx,mesh%ny))
    allocate(mesh%boundType(mesh%nx,mesh%ny),mesh%bound(mesh%nx,mesh%ny))
    allocate(sol%val(mesh%nx,mesh%ny,sol%nvar))

    return
  end subroutine init

  subroutine IC(mesh,sol)
    type(meshStruct), intent(in) :: mesh
    type(solStruct), intent(inout) :: sol
    integer :: i,j,n

    do i=1,mesh%nx
       do j=1,mesh%ny
          if ((mesh%X(i,j)**2+mesh%Y(i,j)**2)**0.5<5) then
             sol%val(i,j,1)=1.
          else
             sol%val(i,j,1)=0.
          endif
       enddo
    enddo
    return
  end subroutine IC

  subroutine BC(xL,xR,yL,yR,mesh)
    real, intent(in) :: xL,xR,yL,yR
    type(meshStruct), intent(inout) :: mesh
    integer :: i,j

    do i=1,mesh%nx
       do j=1,mesh%ny
          if ((mesh%X(i,j)==xL).or.(mesh%Y(i,j)==yL)) then
             mesh%boundType(i,j)='DIRICHLET'
             mesh%bound(i,j)=0.
          elseif ((mesh%X(i,j)==xR).or.(mesh%Y(i,j)==yR)) then
             mesh%boundType(i,j)='DIRICHLET'
             mesh%bound(i,j)=0.
          else
             mesh%boundType(i,j)='NOT A BOUNDARY'
             mesh%bound(i,j)=0.
          endif
       enddo
    enddo
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
       mesh%X(i,:)=xL+i*dx-dx/2.
    enddo
    do i=1,mesh%ny
       mesh%Y(:,i)=yL+i*dy-dy/2
    enddo
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
    write(11,'(a,i4,a,i4,a)')"DIMENSIONS ",mesh%nx," ",mesh%ny," 1"
    write(11,'(a,i8,a)')"POINTS ",(mesh%nx)*(mesh%ny)," float"
    do j=1,mesh%ny
       do i=1,mesh%nx
          write(11,'(e15.8,a,e15.8,a,e15.8)')mesh%X(i,j)," ",mesh%Y(i,j)," ",0.
       enddo
    enddo

    write(11,'(a,i8)')"POINT_DATA ",(mesh%nx)*(mesh%ny)
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

end module inout
