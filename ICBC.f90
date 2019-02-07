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

  subroutine BC_sinus(nvar,mesh)
    integer, intent(in) :: nvar
    type(meshStruct), intent(inout) :: mesh
    integer :: i
    
    do i=1,mesh%ne
       allocate(mesh%edge(i)%bound(nvar))
       if (mesh%edge(i)%cell1<0.or.mesh%edge(i)%cell2<0) then
          mesh%edge(i)%boundType='PERIODIC'
          mesh%edge(i)%bound=0.0_dp
       else
          mesh%edge(i)%boundType='NOT A BOUNDARY'
          mesh%edge(i)%bound=0.0_dp
       endif
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

  subroutine BC_sinus_dis(nvar,mesh)
    integer, intent(in) :: nvar
    type(meshStruct), intent(inout) :: mesh
    integer :: i
    
    do i=1,mesh%ne
       allocate(mesh%edge(i)%bound(nvar))
       if (mesh%edge(i)%cell1<0.or.mesh%edge(i)%cell2<0) then
          mesh%edge(i)%boundType='PERIODIC'
          mesh%edge(i)%bound=0.0_dp
       else
          mesh%edge(i)%boundType='NOT A BOUNDARY'
          mesh%edge(i)%bound=0.0_dp
       endif
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

  subroutine BC_sod(nvar,mesh)
    integer, intent(in) :: nvar
    type(meshStruct), intent(inout) :: mesh
    integer :: i,isol
    real(dp), dimension(:), allocatable :: bound

    allocate(bound(nvar))
    
    do i=1,mesh%ne
       allocate(mesh%edge(i)%bound(nvar))
       if (mesh%edge(i)%cell1<0.and.mesh%edge(i)%dir==1) then
          mesh%edge(i)%boundType='DIRICHLET'
          bound(1)=1.0_dp
          bound(2)=0.0_dp
          bound(3)=0.0_dp
          bound(4)=1.0_dp
          do isol=1,nvar
             call conserv(bound(:),"euler               ",isol,mesh%edge(i)%bound(isol))
          enddo
       elseif (mesh%edge(i)%cell2<0.and.mesh%edge(i)%dir==1) then
          mesh%edge(i)%boundType='DIRICHLET'
          bound(1)=0.125_dp
          bound(2)=0.0_dp
          bound(3)=0.0_dp
          bound(4)=0.1_dp
          do isol=1,nvar
             call conserv(bound(:),"euler               ",isol,mesh%edge(i)%bound(isol))
          enddo
       elseif (mesh%edge(i)%cell1<0.and.mesh%edge(i)%dir==2) then
          mesh%edge(i)%boundType='PERIODIC'
          mesh%edge(i)%bound(:)=0.0_dp
       elseif (mesh%edge(i)%cell2<0.and.mesh%edge(i)%dir==2) then
          mesh%edge(i)%boundType='PERIODIC'
          mesh%edge(i)%bound(:)=0.0_dp
       else
          mesh%edge(i)%boundType='NOT A BOUNDARY'
          mesh%edge(i)%bound=0.0_dp
       endif
    enddo

    deallocate(bound)
    
    return   
  end subroutine BC_sod

!!!!!!!!!!!!!!!!!!!!!!!!!!!!! Sod_2D !!!!!!!!!!!!!!!!!!!!!!!!!!!!!

  subroutine IC_func_sod_2D(x,y,S)
    real(dp), intent(in) :: x,y
    real(dp), dimension(:), intent(inout) :: S
    real(dp), dimension(:), allocatable :: U
    integer :: i
    real(dp) :: r,r0

    allocate(U(size(S)))

    r0=0.25_dp
    r=sqrt(x**2+y**2)
    
    if (r>r0) then
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
  end subroutine IC_func_sod_2D

  subroutine BC_sod_2D(nvar,mesh)
    integer, intent(in) :: nvar
    type(meshStruct), intent(inout) :: mesh
    integer :: i,isol
    real(dp), dimension(:), allocatable :: bound

    allocate(bound(nvar))
    
    do i=1,mesh%ne
       allocate(mesh%edge(i)%bound(nvar))
       if (mesh%edge(i)%cell1<0.or.mesh%edge(i)%cell2<0) then
          mesh%edge(i)%boundType='DIRICHLET'
          bound(1)=1.0_dp
          bound(2)=0.0_dp
          bound(3)=0.0_dp
          bound(4)=1.0_dp
          do isol=1,nvar
             call conserv(bound(:),"euler               ",isol,mesh%edge(i)%bound(isol))
          enddo
       else
          mesh%edge(i)%boundType='NOT A BOUNDARY'
          mesh%edge(i)%bound=0.0_dp
       endif
    enddo

    deallocate(bound)
    
    return   
  end subroutine BC_sod_2D

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

  subroutine BC_sod_mod(nvar,mesh)
    integer, intent(in) :: nvar
    type(meshStruct), intent(inout) :: mesh
    integer :: i,isol
    real(dp), dimension(:), allocatable :: bound

    allocate(bound(nvar))
    
    do i=1,mesh%ne
       allocate(mesh%edge(i)%bound(nvar))
       if (mesh%edge(i)%cell1<0.and.mesh%edge(i)%dir==1) then
          mesh%edge(i)%boundType='DIRICHLET'
          bound(1)=1.0_dp
          bound(2)=0.75_dp
          bound(3)=0.0_dp
          bound(4)=1.0_dp
          do isol=1,nvar
             call conserv(bound(:),"euler               ",isol,mesh%edge(i)%bound(isol))
          enddo
       elseif (mesh%edge(i)%cell2<0.and.mesh%edge(i)%dir==1) then
          mesh%edge(i)%boundType='DIRICHLET'
          bound(1)=0.125_dp
          bound(2)=0.0_dp
          bound(3)=0.0_dp
          bound(4)=0.1_dp
          do isol=1,nvar
             call conserv(bound(:),"euler               ",isol,mesh%edge(i)%bound(isol))
          enddo
       elseif (mesh%edge(i)%cell1<0.and.mesh%edge(i)%dir==2) then
          mesh%edge(i)%boundType='PERIODIC'
          mesh%edge(i)%bound(:)=0.0_dp
       elseif (mesh%edge(i)%cell2<0.and.mesh%edge(i)%dir==2) then
          mesh%edge(i)%boundType='PERIODIC'
          mesh%edge(i)%bound(:)=0.0_dp
       else
          mesh%edge(i)%boundType='NOT A BOUNDARY'
          mesh%edge(i)%bound=0.0_dp
       endif
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

  subroutine BC_shu(nvar,mesh)
    integer, intent(in) :: nvar
    type(meshStruct), intent(inout) :: mesh
    integer :: i,isol
    real(dp), dimension(:), allocatable :: bound

    allocate(bound(nvar))
    
    do i=1,mesh%ne
       allocate(mesh%edge(i)%bound(nvar))
       if (mesh%edge(i)%cell1<0.and.mesh%edge(i)%dir==1) then
          mesh%edge(i)%boundType='DIRICHLET'
          bound(1)=3.857143_dp
          bound(2)=2.629369_dp
          bound(3)=0.0_dp
          bound(4)=10.33333_dp
          do isol=1,nvar
             call conserv(bound(:),"euler               ",isol,mesh%edge(i)%bound(isol))
          enddo
       elseif (mesh%edge(i)%cell2<0.and.mesh%edge(i)%dir==1) then
          mesh%edge(i)%boundType='DIRICHLET'
          bound(1)=1.0_dp+0.2_dp*sin(25.0_dp)
          bound(2)=0.0_dp
          bound(3)=0.0_dp
          bound(4)=1.0_dp
          do isol=1,nvar
             call conserv(bound(:),"euler               ",isol,mesh%edge(i)%bound(isol))
          enddo
       elseif (mesh%edge(i)%cell1<0.and.mesh%edge(i)%dir==2) then
          mesh%edge(i)%boundType='PERIODIC'
          mesh%edge(i)%bound(:)=0.0_dp
       elseif (mesh%edge(i)%cell2<0.and.mesh%edge(i)%dir==2) then
          mesh%edge(i)%boundType='PERIODIC'
          mesh%edge(i)%bound(:)=0.0_dp
       else
          mesh%edge(i)%boundType='NOT A BOUNDARY'
          mesh%edge(i)%bound=0.0_dp
       endif
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

  subroutine BC_123(nvar,mesh)
    integer, intent(in) :: nvar
    type(meshStruct), intent(inout) :: mesh
    integer :: i,isol
    real(dp), dimension(:), allocatable :: bound

    allocate(bound(nvar))
    
    do i=1,mesh%ne
       allocate(mesh%edge(i)%bound(nvar))
       if (mesh%edge(i)%cell1<0.and.mesh%edge(i)%dir==1) then
          mesh%edge(i)%boundType='DIRICHLET'
          bound(1)=1.0_dp
          bound(2)=-2.0_dp
          bound(3)=0.0_dp
          bound(4)=0.4_dp
          do isol=1,nvar
             call conserv(bound(:),"euler               ",isol,mesh%edge(i)%bound(isol))
          enddo
       elseif (mesh%edge(i)%cell2<0.and.mesh%edge(i)%dir==1) then
          mesh%edge(i)%boundType='DIRICHLET'
          bound(1)=1_dp
          bound(2)=2.0_dp
          bound(3)=0.0_dp
          bound(4)=0.4_dp
          do isol=1,nvar
             call conserv(bound(:),"euler               ",isol,mesh%edge(i)%bound(isol))
          enddo
       elseif (mesh%edge(i)%cell1<0.and.mesh%edge(i)%dir==2) then
          mesh%edge(i)%boundType='PERIODIC'
          mesh%edge(i)%bound(:)=0.0_dp
       elseif (mesh%edge(i)%cell2<0.and.mesh%edge(i)%dir==2) then
          mesh%edge(i)%boundType='PERIODIC'
          mesh%edge(i)%bound(:)=0.0_dp
       else
          mesh%edge(i)%boundType='NOT A BOUNDARY'
          mesh%edge(i)%bound=0.0_dp
       endif
    enddo

    deallocate(bound)
    
    return   
  end subroutine BC_123

!!!!!!!!!!!!!!!!!!!!!!!!!!!!! Vortex !!!!!!!!!!!!!!!!!!!!!!!!!!!!!
  
  subroutine IC_func_vortex(x,y,S)
    real(dp), intent(in) :: x,y
    real(dp), dimension(:), intent(inout) :: S
    real(dp), dimension(:), allocatable :: U
    integer :: i
    real(dp) :: beta,r,T

    allocate(U(size(S)))

    beta=5.0_dp
    r=sqrt(x**2+y**2)
    T=1.0_dp-(gamma-1.0_dp)*beta**2*exp(1-r**2)/(8.0_dp*gamma*(pi**2))
    U(1)=T**(1.0_dp/(gamma-1.0_dp))
    U(2)=0.0_dp-y*beta*exp(0.5_dp*(1-r**2))/(2.0_dp*pi)
    U(3)=0.0_dp+x*beta*exp(0.5_dp*(1-r**2))/(2.0_dp*pi)
    U(4)=U(1)**gamma

    do i=1,4
       call conserv(U,"euler               ",i,S(i))
    enddo

    deallocate(U)
    
    return
  end subroutine IC_func_vortex

  subroutine BC_vortex(nvar,mesh)
    integer, intent(in) :: nvar
    type(meshStruct), intent(inout) :: mesh
    integer :: i
    
    do i=1,mesh%ne
       allocate(mesh%edge(i)%bound(nvar))
       if (mesh%edge(i)%cell1<0.or.mesh%edge(i)%cell2<0) then
          mesh%edge(i)%boundType='PERIODIC'
          mesh%edge(i)%bound=0.0_dp
       else
          mesh%edge(i)%boundType='NOT A BOUNDARY'
          mesh%edge(i)%bound=0.0_dp
       endif
    enddo
    
    return   
  end subroutine BC_vortex

!!!!!!!!!!!!!!!!!!!!!!!!!!!!! RP2D_3 !!!!!!!!!!!!!!!!!!!!!!!!!!!!!

  subroutine IC_func_RP2D_3(x,y,S)
    real(dp), intent(in) :: x,y
    real(dp), dimension(:), intent(inout) :: S
    real(dp), dimension(:), allocatable :: U
    integer :: i
    real(dp) :: x0,y0

    allocate(U(size(S)))

    x0=0.5_dp
    y0=0.5_dp
    
    if (x<=x0.and.y<=y0) then
       U(1)=0.138_dp
       U(2)=1.206_dp
       U(3)=1.206_dp
       U(4)=0.029_dp
    elseif (x<=x0.and.y>y0) then
       U(1)=0.5323_dp
       U(2)=1.206_dp
       U(3)=0.0_dp
       U(4)=0.3_dp
    elseif (x>x0.and.y<=y0) then
       U(1)=0.5323_dp
       U(2)=0.0_dp
       U(3)=1.206_dp
       U(4)=0.3_dp
    else
       U(1)=1.5_dp
       U(2)=0.0_dp
       U(3)=0.0_dp
       U(4)=1.5_dp
    endif
    
    do i=1,4
       call conserv(U,"euler               ",i,S(i))
    enddo

    deallocate(U)
    
    return
  end subroutine IC_func_RP2D_3

  subroutine BC_RP2D_3(nvar,mesh)
    integer, intent(in) :: nvar
    type(meshStruct), intent(inout) :: mesh
    integer :: i
    real(dp), dimension(:), allocatable :: bound

    allocate(bound(nvar))
    
    do i=1,mesh%ne
       allocate(mesh%edge(i)%bound(nvar))
       if (mesh%edge(i)%cell1<0.or.mesh%edge(i)%cell2<0) then
          mesh%edge(i)%boundType='NEUMANN'
          mesh%edge(i)%bound=0.0_dp
       else
          mesh%edge(i)%boundType='NOT A BOUNDARY'
          mesh%edge(i)%bound=0.0_dp
       endif
    enddo

    deallocate(bound)
    
    return   
  end subroutine BC_RP2D_3

!!!!!!!!!!!!!!!!!!!!!!!!!!!!! Sinus !!!!!!!!!!!!!!!!!!!!!!!!!!!!!

  subroutine IC_func_test(x,y,S)
    real(dp), intent(in) :: x,y
    real(dp), dimension(:), intent(inout) :: S
    if(.false.)print*,y

    if (x<1.5_dp) then
       S(1)=1.0_dp
    else
       S(1)=0.0_dp
    endif
    
    return
  end subroutine IC_func_test

  subroutine BC_test(nvar,mesh)
    integer, intent(in) :: nvar
    type(meshStruct), intent(inout) :: mesh
    integer :: i
    
    do i=1,mesh%ne
       allocate(mesh%edge(i)%bound(nvar))
       if (mesh%edge(i)%cell1<0.or.mesh%edge(i)%cell2<0) then
          mesh%edge(i)%boundType='NEUMANN'
          mesh%edge(i)%bound=0.0_dp
       else
          mesh%edge(i)%boundType='NOT A BOUNDARY'
          mesh%edge(i)%bound=0.0_dp
       endif
    enddo
    
    return   
  end subroutine BC_test


end module ICBC
