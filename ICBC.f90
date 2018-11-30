module ICBC

  use constant
  use types
  use phys

  implicit none

contains

!!!!!!!!!!!!!!!!!!!!!!!!!!!!! Sinus !!!!!!!!!!!!!!!!!!!!!!!!!!!!!

  subroutine IC_func_sinus(x,y,S)
    real(dp), intent(in) :: x,y
    real(dp), dimension(:), intent(inout) :: S

    S(1)=cos((x-5.0_dp)*pi/5.0_dp)+cos((y-5.0_dp)*pi/5.0_dp)
    
    return
  end subroutine IC_func_sinus

  subroutine BC_sinus(nx,ny,nvar,mesh)
    integer, intent(in) :: nx,ny,nvar
    type(meshStruct), intent(inout) :: mesh
    integer :: i,j
    
    do i=1,mesh%ne
       allocate(mesh%edge(i)%bound(nvar))
       mesh%edge(i)%boundType='NOT A BOUNDARY'
       mesh%edge(i)%bound=0.0_dp
    enddo
    
    do j=1,ny
       mesh%edge((j-1)*(nx+1)+1)%boundType='PERIODIC'
       mesh%edge((j-1)*(nx+1)+1)%bound=0.0_dp
       
       mesh%edge(j*(nx+1))%boundType='PERIODIC'
       mesh%edge(j*(nx+1))%bound=0.0_dp
       
    enddo
    
    do i=1,nx
       mesh%edge((i-1)*(ny+1)+1+(nx+1)*ny)%boundType='PERIODIC'
       mesh%edge((i-1)*(ny+1)+1+(nx+1)*ny)%bound(:)=0.0_dp
       
       mesh%edge(i*(ny+1)+(nx+1)*ny)%boundType='PERIODIC'
       mesh%edge(i*(ny+1)+(nx+1)*ny)%bound(:)=0.0_dp
    enddo
    
    return   
  end subroutine BC_sinus

!!!!!!!!!!!!!!!!!!!!!!!!!!!!! Sinus !!!!!!!!!!!!!!!!!!!!!!!!!!!!!

  subroutine IC_func_sinus_dis(x,y,S)
    real(dp), intent(in) :: x,y
    real(dp), dimension(:), intent(inout) :: S

    S(1)=cos((x-5.0_dp)*pi/5.0_dp)+cos((y-5.0_dp)*pi/5.0_dp)

    if (S(1)>1.5_dp) then
       S(1)=1.0_dp
    else
       S(1)=0.0_dp
    endif
    
    return
  end subroutine IC_func_sinus_dis

  subroutine BC_sinus_dis(nx,ny,nvar,mesh)
    integer, intent(in) :: nx,ny,nvar
    type(meshStruct), intent(inout) :: mesh
    integer :: i,j
    
    do i=1,mesh%ne
       allocate(mesh%edge(i)%bound(nvar))
       mesh%edge(i)%boundType='NOT A BOUNDARY'
       mesh%edge(i)%bound=0.0_dp
    enddo
    
    do j=1,ny
       mesh%edge((j-1)*(nx+1)+1)%boundType='PERIODIC'
       mesh%edge((j-1)*(nx+1)+1)%bound=0.0_dp
       
       mesh%edge(j*(nx+1))%boundType='PERIODIC'
       mesh%edge(j*(nx+1))%bound=0.0_dp
       
    enddo
    
    do i=1,nx
       mesh%edge((i-1)*(ny+1)+1+(nx+1)*ny)%boundType='PERIODIC'
       mesh%edge((i-1)*(ny+1)+1+(nx+1)*ny)%bound(:)=0.0_dp
       
       mesh%edge(i*(ny+1)+(nx+1)*ny)%boundType='PERIODIC'
       mesh%edge(i*(ny+1)+(nx+1)*ny)%bound(:)=0.0_dp
    enddo
    
    return   
  end subroutine BC_sinus_dis
  
!!!!!!!!!!!!!!!!!!!!!!!!!!!!! Sod !!!!!!!!!!!!!!!!!!!!!!!!!!!!!

  subroutine IC_func_sod(x,y,S)
    real(dp), intent(in) :: x,y
    real(dp), dimension(:), intent(inout) :: S
    real(dp), dimension(:), allocatable :: U
    integer :: i
    real(dp) :: x0
    if(.false.)print*,y

    allocate(U(size(S)))

    x0=0.5_dp
    
    if (x<x0) then
       U(1)=1.0_dp
       U(2)=0.0_dp
       U(3)=0.0_dp
       U(4)=1.0_dp
    else
       U(1)=0.125_dp
       U(2)=0.0_dp
       U(3)=0.0_dp
       U(4)=0.1_dp
    endif
    
    do i=1,4
       call conserv(U,"euler               ",i,S(i))
    enddo

    deallocate(U)
    
    return
  end subroutine IC_func_sod

  subroutine BC_sod(nx,ny,nvar,mesh)
    integer, intent(in) :: nx,ny,nvar
    type(meshStruct), intent(inout) :: mesh
    integer :: i,j,isol
    real(dp), dimension(:), allocatable :: bound

    allocate(bound(nvar))
    
    do i=1,mesh%ne
       allocate(mesh%edge(i)%bound(nvar))
       mesh%edge(i)%boundType='NOT A BOUNDARY'
       mesh%edge(i)%bound(:)=0.0_dp
    enddo
    
    do j=1,ny
       mesh%edge((j-1)*(nx+1)+1)%boundType='DIRICHLET'
       bound(1)=1.0_dp
       bound(2)=0.0_dp
       bound(3)=0.0_dp
       bound(4)=1.0_dp
       do isol=1,nvar
          call conserv(bound(:),"euler               ",isol,mesh%edge((j-1)*(nx+1)+1)%bound(isol))
       enddo
       
       mesh%edge(j*(nx+1))%boundType='DIRICHLET'
       bound(1)=0.125_dp
       bound(2)=0.0_dp
       bound(3)=0.0_dp
       bound(4)=0.1_dp
       do isol=1,nvar
          call conserv(bound(:),"euler               ",isol,mesh%edge(j*(nx+1))%bound(isol))
       enddo
       
    enddo
    
    do i=1,nx
       mesh%edge((i-1)*(ny+1)+1+(nx+1)*ny)%boundType='PERIODIC'
       mesh%edge((i-1)*(ny+1)+1+(nx+1)*ny)%bound(:)=0.0_dp
       
       mesh%edge(i*(ny+1)+(nx+1)*ny)%boundType='PERIODIC'
       mesh%edge(i*(ny+1)+(nx+1)*ny)%bound(:)=0.0_dp
    enddo

    deallocate(bound)
    
    return   
  end subroutine BC_sod

!!!!!!!!!!!!!!!!!!!!!!!!!!!!! Sod_mod !!!!!!!!!!!!!!!!!!!!!!!!!!!!!

  subroutine IC_func_sod_mod(x,y,S)
    real(dp), intent(in) :: x,y
    real(dp), dimension(:), intent(inout) :: S
    real(dp), dimension(:), allocatable :: U
    integer :: i
    real(dp) :: x0
    if(.false.)print*,y

    allocate(U(size(S)))

    x0=0.3_dp
    
    if (x<x0) then
       U(1)=1.0_dp
       U(2)=0.75_dp
       U(3)=0.0_dp
       U(4)=1.0_dp
    else
       U(1)=0.125_dp
       U(2)=0.0_dp
       U(3)=0.0_dp
       U(4)=0.1_dp
    endif
    
    do i=1,4
       call conserv(U,"euler               ",i,S(i))
    enddo

    deallocate(U)
    
    return
  end subroutine IC_func_sod_mod

  subroutine BC_sod_mod(nx,ny,nvar,mesh)
    integer, intent(in) :: nx,ny,nvar
    type(meshStruct), intent(inout) :: mesh
    integer :: i,j,isol
    real(dp), dimension(:), allocatable :: bound

    allocate(bound(nvar))
    
    do i=1,mesh%ne
       allocate(mesh%edge(i)%bound(nvar))
       mesh%edge(i)%boundType='NOT A BOUNDARY'
       mesh%edge(i)%bound(:)=0.0_dp
    enddo
    
    do j=1,ny
       mesh%edge((j-1)*(nx+1)+1)%boundType='DIRICHLET'
       bound(1)=1.0_dp
       bound(2)=0.75_dp
       bound(3)=0.0_dp
       bound(4)=1.0_dp
       do isol=1,nvar
          call conserv(bound(:),"euler               ",isol,mesh%edge((j-1)*(nx+1)+1)%bound(isol))
       enddo
       
       mesh%edge(j*(nx+1))%boundType='DIRICHLET'
       bound(1)=0.125_dp
       bound(2)=0.0_dp
       bound(3)=0.0_dp
       bound(4)=0.1_dp
       do isol=1,nvar
          call conserv(bound(:),"euler               ",isol,mesh%edge(j*(nx+1))%bound(isol))
       enddo
       
    enddo
    
    do i=1,nx
       mesh%edge((i-1)*(ny+1)+1+(nx+1)*ny)%boundType='PERIODIC'
       mesh%edge((i-1)*(ny+1)+1+(nx+1)*ny)%bound(:)=0.0_dp
       
       mesh%edge(i*(ny+1)+(nx+1)*ny)%boundType='PERIODIC'
       mesh%edge(i*(ny+1)+(nx+1)*ny)%bound(:)=0.0_dp
    enddo

    deallocate(bound)
    
    return   
  end subroutine BC_sod_mod
  
!!!!!!!!!!!!!!!!!!!!!!!!!!!!! Shu !!!!!!!!!!!!!!!!!!!!!!!!!!!!!
  
  subroutine IC_func_shu(x,y,S)
    real(dp), intent(in) :: x,y
    real(dp), dimension(:), intent(inout) :: S
    real(dp), dimension(:), allocatable :: U
    integer :: i
    real(dp) :: x0
    if(.false.)print*,y

    allocate(U(size(S)))

    x0=-4.0_dp
    
    if (x<x0) then
       U(1)=3.857143_dp
       U(2)=2.629369_dp
       U(3)=0.0_dp
       U(4)=10.33333_dp
    else
       U(1)=1.0_dp+0.2_dp*sin(5.0_dp*x)
       U(2)=0.0_dp
       U(3)=0.0_dp
       U(4)=1.0_dp
    endif
    
    do i=1,4
       call conserv(U,"euler               ",i,S(i))
    enddo

    deallocate(U)
    
    return
  end subroutine IC_func_shu

  subroutine BC_shu(nx,ny,nvar,mesh)
    integer, intent(in) :: nx,ny,nvar
    type(meshStruct), intent(inout) :: mesh
    integer :: i,j,isol
    real(dp), dimension(:), allocatable :: bound

    allocate(bound(nvar))
    
    do i=1,mesh%ne
       allocate(mesh%edge(i)%bound(nvar))
       mesh%edge(i)%boundType='NOT A BOUNDARY'
       mesh%edge(i)%bound(:)=0.0_dp
    enddo
    
    do j=1,ny
       mesh%edge((j-1)*(nx+1)+1)%boundType='DIRICHLET'
       bound(1)=3.857143_dp
       bound(2)=2.629369_dp
       bound(3)=0.0_dp
       bound(4)=10.33333_dp
       do isol=1,nvar
          call conserv(bound(:),"euler               ",isol,mesh%edge((j-1)*(nx+1)+1)%bound(isol))
       enddo
       
       mesh%edge(j*(nx+1))%boundType='DIRICHLET'
       bound(1)=1.0_dp+0.2_dp*sin(25.0_dp)
       bound(2)=0.0_dp
       bound(3)=0.0_dp
       bound(4)=1.0_dp
       do isol=1,nvar
          call conserv(bound(:),"euler               ",isol,mesh%edge(j*(nx+1))%bound(isol))
       enddo
       
    enddo
    
    do i=1,nx
       mesh%edge((i-1)*(ny+1)+1+(nx+1)*ny)%boundType='PERIODIC'
       mesh%edge((i-1)*(ny+1)+1+(nx+1)*ny)%bound(:)=0.0_dp
       
       mesh%edge(i*(ny+1)+(nx+1)*ny)%boundType='PERIODIC'
       mesh%edge(i*(ny+1)+(nx+1)*ny)%bound(:)=0.0_dp
    enddo

    deallocate(bound)
    
    return   
  end subroutine BC_shu

!!!!!!!!!!!!!!!!!!!!!!!!!!!!! 123 !!!!!!!!!!!!!!!!!!!!!!!!!!!!!

  subroutine IC_func_123(x,y,S)
    real(dp), intent(in) :: x,y
    real(dp), dimension(:), intent(inout) :: S
    real(dp), dimension(:), allocatable :: U
    integer :: i
    real(dp) :: x0
    if(.false.)print*,y

    allocate(U(size(S)))

    x0=0.5_dp
    
    if (x<x0) then
       U(1)=1.0_dp
       U(2)=-2.0_dp
       U(3)=0.0_dp
       U(4)=0.4_dp
    else
       U(1)=1.0_dp
       U(2)=2.0_dp
       U(3)=0.0_dp
       U(4)=0.4_dp
    endif
    
    do i=1,4
       call conserv(U,"euler               ",i,S(i))
    enddo

    deallocate(U)
    
    return
  end subroutine IC_func_123

  subroutine BC_123(nx,ny,nvar,mesh)
    integer, intent(in) :: nx,ny,nvar
    type(meshStruct), intent(inout) :: mesh
    integer :: i,j,isol
    real(dp), dimension(:), allocatable :: bound

    allocate(bound(nvar))
    
    do i=1,mesh%ne
       allocate(mesh%edge(i)%bound(nvar))
       mesh%edge(i)%boundType='NOT A BOUNDARY'
       mesh%edge(i)%bound(:)=0.0_dp
    enddo
    
    do j=1,ny
       mesh%edge((j-1)*(nx+1)+1)%boundType='DIRICHLET'
       bound(1)=1.0_dp
       bound(2)=-2.0_dp
       bound(3)=0.0_dp
       bound(4)=0.4_dp
       do isol=1,nvar
          call conserv(bound(:),"euler               ",isol,mesh%edge((j-1)*(nx+1)+1)%bound(isol))
       enddo
       
       mesh%edge(j*(nx+1))%boundType='DIRICHLET'
       bound(1)=1_dp
       bound(2)=2.0_dp
       bound(3)=0.0_dp
       bound(4)=0.4_dp
       do isol=1,nvar
          call conserv(bound(:),"euler               ",isol,mesh%edge(j*(nx+1))%bound(isol))
       enddo
       
    enddo
    
    do i=1,nx
       mesh%edge((i-1)*(ny+1)+1+(nx+1)*ny)%boundType='PERIODIC'
       mesh%edge((i-1)*(ny+1)+1+(nx+1)*ny)%bound(:)=0.0_dp
       
       mesh%edge(i*(ny+1)+(nx+1)*ny)%boundType='PERIODIC'
       mesh%edge(i*(ny+1)+(nx+1)*ny)%bound(:)=0.0_dp
    enddo

    deallocate(bound)
    
    return   
  end subroutine BC_123

end module ICBC
