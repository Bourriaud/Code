module inout
  
  implicit none
  
contains

  subroutine init(xL,xR,yL,yR,nx,ny,dx,dy,nsol,namefile,solname)
    real, intent(out) :: xL,xR,yL,yR,dx,dy
    integer, intent(out) :: nx,ny,nsol
    character(len=20), intent(out) :: namefile
    character(len=20), dimension(:), allocatable, intent(out) :: solname
    integer :: i

    open(11,file="configuration",form="formatted")
    read(11,*)xL
    read(11,*)xR
    read(11,*)yL
    read(11,*)yR
    read(11,*)nx
    read(11,*)ny
    read(11,*)namefile
    read(11,*)nsol
    allocate(solname(nsol))
    do i=1,nsol
       read(11,*)solname(i)
    enddo
    close(11)
    dx=(xR-xL)/nx
    dy=(yR-yL)/ny
    return
  end subroutine init
  
  subroutine buildMesh(xL,xR,yL,yR,nx,ny,Ix,Iy)
    real, intent(in) :: xL,xR,yL,yR
    integer, intent(in) :: nx,ny
    real, dimension(0:nx), intent(out) :: Ix
    real, dimension(0:ny), intent(out) :: Iy
    integer :: i
    real :: dx,dy

    dx=(xR-xL)/nx
    dy=(yR-yL)/ny
    do i=0,nx
       Ix(i)=xL+i*dx
    enddo
    do i=0,ny
       Iy(i)=yL+i*dy
    enddo
    return
  end subroutine buildMesh
  
  subroutine writeSol(Ix,Iy,nsol,sol,namefile,solname)
    integer, intent(in) :: nsol
    real, dimension(0:), intent(in) :: Ix,Iy
    real, dimension(0:,0:,:), intent(in) :: sol
    character(len=20), intent(in) :: namefile
    character(len=20), dimension(:), intent(in) :: solname
    integer :: i,j,n
    integer :: nx,ny

    nx=size(Ix)-1
    ny=size(Iy)-1

    open(11,file=namefile,form="formatted")
    write(11,*)"# vtk DataFile Version 2.0"
    write(11,*)"Results of the calculation"
    write(11,*)"ASCII"
    write(11,*)"DATASET STRUCTURED_GRID"
    write(11,*)"DIMENSIONS ",nx+1," ",ny+1," 1"
    write(11,*)"POINTS ",(nx+1)*(ny+1)," float"
    do i=0,nx
       do j=0,ny
          write(11,*)Ix(i)," ",Iy(j)," ",0.
       enddo
    enddo
    write(11,*)"POINT_DATA ",(nx+1)*(ny+1)
    do n=1,nsol
       write(11,*)"SCALARS ",solname(n)," float 1"
       write(11,*)"LOOKUP_TABLE default"
       do j=0,ny
          do i=0,nx
             write(11,*)sol(i,j,n)
          enddo
       enddo
    enddo
    close(11)
  return
end subroutine writeSol

end module inout
