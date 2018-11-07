module efficiency

  use constant
  use types
  
  implicit none

  interface quadrature1
     module procedure quadrature1_t,quadrature1_c_alpha,quadrature1_reconstruction
  end interface quadrature1
  
  interface quadrature3
     module procedure quadrature3_t,quadrature3_c_alpha,quadrature3_reconstruction
  end interface quadrature3
  
contains

  subroutine exactSol(x,y,t,s)
    real(dp), intent(in) :: x,y,t
    real(dp), intent(out) :: s
    real(dp) :: a1,a2

    a1=1.0_dp
    a2=1.0_dp

    s=cos((x-a1*t-5.0_dp)*pi/5.0_dp)+cos((y-a2*t-5.0_dp)*pi/5.0_dp)

    return
  end subroutine exactSol
  
  subroutine exactTab(t,mesh,tab)
    real(dp), intent(in) :: t
    type (meshStruct), intent(in) :: mesh
    real(dp), dimension(:,:), intent(inout) :: tab
    integer :: k
    real(dp) :: dx,dy

    do k=1,mesh%nc
       dx=mesh%cell(k)%dx
       dy=mesh%cell(k)%dy
       call quadrature3(exactSol,mesh,t,k,tab(k,1))
       tab(k,1)=tab(k,1)/(dx*dy)
    enddo
    
    return
  end subroutine exactTab

  subroutine userSol(t,mesh,sol)
    real(dp), intent(in) :: t
    type(meshStruct), intent(in) :: mesh
    type(solStruct), intent(inout) :: sol

    call exactTab(t,mesh,sol%user(:,1:1))
    sol%user(:,2:2)=abs(sol%user(:,1:1)-sol%val(:,2:2))

    return
  end subroutine userSol

  subroutine errorL1(mesh,sol,exactsol,eL1)
    type(meshStruct), intent(in) :: mesh
    real(dp), dimension(:), intent(in) :: sol,exactsol
    real(dp), intent(out) :: eL1
    integer :: k
    real(dp) :: dx,dy

    eL1=0.0_dp
    do k=1,size(sol(:))
       dx=mesh%cell(k)%dx
       dy=mesh%cell(k)%dy
       eL1=eL1+abs(sol(k)-exactsol(k))*dx*dy
    enddo

    return
  end subroutine errorL1

  subroutine errorL2(mesh,sol,exactsol,eL2)
    type(meshStruct), intent(in) :: mesh
    real(dp), dimension(:), intent(in) :: sol,exactsol
    real(dp), intent(out) :: eL2
    integer :: k
    real(dp) :: dx,dy

    eL2=0.0_dp
    do k=1,size(sol(:))
       dx=mesh%cell(k)%dx
       dy=mesh%cell(k)%dy
       eL2=eL2+dx*dy*(sol(k)-exactsol(k))**2
    enddo
    eL2=sqrt(eL2)

    return
  end subroutine errorL2

  subroutine quadrature1_t(func,mesh,t,k,int)
    procedure(sub_quadra_t) :: func
    type(meshStruct), intent(in) :: mesh
    real(dp), intent(in) :: t
    integer, intent(in) :: k
    real(dp), intent(out) :: int
    real(dp) :: x,y

    x=mesh%cell(k)%xc
    y=mesh%cell(k)%yc
    call func(x,y,t,int)
    int=int*mesh%cell(k)%dx*mesh%cell(k)%dy

    return
  end subroutine quadrature1_t

  subroutine quadrature1_c_alpha(func,mesh,c,alpha,k,int)
    procedure(sub_quadra_c_alpha) :: func
    type(meshStruct), intent(in) :: mesh
    real(dp), dimension(2), intent(in) :: c
    integer, dimension(2), intent(in) :: alpha
    integer, intent(in) :: k
    real(dp), intent(out) :: int
    real(dp) :: x,y

    x=mesh%cell(k)%xc
    y=mesh%cell(k)%yc
    call func(x,y,c,alpha,int)
    int=int*mesh%cell(k)%dx*mesh%cell(k)%dy

    return
  end subroutine quadrature1_c_alpha

  subroutine quadrature1_reconstruction(func,mesh,sol,order,normal,k,int)
    procedure(sub_reconstruction) :: func
    type(meshStruct), intent(in) :: mesh
    type(solStruct), intent(in) :: sol
    integer, intent(in) :: order,normal,k
    real(dp), dimension(:), intent(inout) :: int
    real(dp) :: x,y,dx,dy

    x=mesh%cell(k)%xc
    y=mesh%cell(k)%yc
    dx=mesh%cell(k)%dx
    dy=mesh%cell(k)%dy

    select case (normal)
    case(1)
       call func(mesh,sol,k,order,x-dx/2.0_dp,y,int)
       print*,int,dy
       int=int*dy
    case(2)
       call func(mesh,sol,k,order,x,y-dy/2.0_dp,int)
       int=int*dx
    case(3)
       call func(mesh,sol,k,order,x+dx/2.0_dp,y,int)
       int=int*dy
    case(4)
       call func(mesh,sol,k,order,x,y+dy/2.0_dp,int)
       int=int*dx
    end select
    
    return
  end subroutine quadrature1_reconstruction

  subroutine quadrature3_t(func,mesh,t,k,int)
    procedure(sub_quadra_t) :: func
    type(meshStruct), intent(in) :: mesh
    real(dp), intent(in) :: t
    integer, intent(in) :: k
    real(dp), intent(out) :: int
    real(dp) :: x,y,dx,dy,ax,ay,sq
    real(dp) :: f1,f2,f3,f4,f5,f6,f7,f8,f9

    x=mesh%cell(k)%xc
    y=mesh%cell(k)%yc
    dx=mesh%cell(k)%dx
    dy=mesh%cell(k)%dy
    ax=x-dx/2.0_dp
    ay=y-dy/2.0_dp
    sq=sqrt(0.6_dp)
    
    call func(ax+dx/2.0_dp*(1.0_dp-sq),ay+dy/2.0_dp*(1.0_dp-sq),t,f1)
    int=25.0_dp*f1
    call func(ax+dx/2.0_dp*(1.0_dp),ay+dy/2.0_dp*(1.0_dp-sq),t,f2)
    int=int+40.0_dp*f2
    call func(ax+dx/2.0_dp*(1.0_dp+sq),ay+dy/2.0_dp*(1.0_dp-sq),t,f3)
    int=int+25.0_dp*f3
    call func(ax+dx/2.0_dp*(1.0_dp-sq),ay+dy/2.0_dp*(1.0_dp),t,f4)
    int=int+40.0_dp*f4
    call func(ax+dx/2.0_dp*(1.0_dp),ay+dy/2.0_dp*(1.0_dp),t,f5)
    int=int+64.0_dp*f5
    call func(ax+dx/2.0_dp*(1.0_dp+sq),ay+dy/2.0_dp*(1.0_dp),t,f6)
    int=int+40.0_dp*f6
    call func(ax+dx/2.0_dp*(1.0_dp-sq),ay+dy/2.0_dp*(1.0_dp+sq),t,f7)
    int=int+25.0_dp*f7
    call func(ax+dx/2.0_dp*(1.0_dp),ay+dy/2.0_dp*(1.0_dp+sq),t,f8)
    int=int+40.0_dp*f8
    call func(ax+dx/2.0_dp*(1.0_dp+sq),ay+dy/2.0_dp*(1.0_dp+sq),t,f9)
    int=int+25.0_dp*f9

    int=int/(4.0_dp*81.0_dp)
    int=int*mesh%cell(k)%dx*mesh%cell(k)%dy

    return
  end subroutine quadrature3_t

  subroutine quadrature3_c_alpha(func,mesh,c,alpha,k,int)
    procedure(sub_quadra_c_alpha) :: func
    type(meshStruct), intent(in) :: mesh
    real(dp), dimension(2), intent(in) :: c
    integer, dimension(2), intent(in) :: alpha
    integer, intent(in) :: k
    real(dp), intent(out) :: int
    real(dp) :: x,y,dx,dy,ax,ay,sq
    real(dp) :: f1,f2,f3,f4,f5,f6,f7,f8,f9

    x=mesh%cell(k)%xc
    y=mesh%cell(k)%yc
    dx=mesh%cell(k)%dx
    dy=mesh%cell(k)%dy
    ax=x-dx/2.0_dp
    ay=y-dy/2.0_dp
    sq=sqrt(0.6_dp)
    
    call func(ax+dx/2.0_dp*(1.0_dp-sq),ay+dy/2.0_dp*(1.0_dp-sq),c,alpha,f1)
    int=25.0_dp*f1
    call func(ax+dx/2.0_dp*(1.0_dp),ay+dy/2.0_dp*(1.0_dp-sq),c,alpha,f2)
    int=int+40.0_dp*f2
    call func(ax+dx/2.0_dp*(1.0_dp+sq),ay+dy/2.0_dp*(1.0_dp-sq),c,alpha,f3)
    int=int+25.0_dp*f3
    call func(ax+dx/2.0_dp*(1.0_dp-sq),ay+dy/2.0_dp*(1.0_dp),c,alpha,f4)
    int=int+40.0_dp*f4
    call func(ax+dx/2.0_dp*(1.0_dp),ay+dy/2.0_dp*(1.0_dp),c,alpha,f5)
    int=int+64.0_dp*f5
    call func(ax+dx/2.0_dp*(1.0_dp+sq),ay+dy/2.0_dp*(1.0_dp),c,alpha,f6)
    int=int+40.0_dp*f6
    call func(ax+dx/2.0_dp*(1.0_dp-sq),ay+dy/2.0_dp*(1.0_dp+sq),c,alpha,f7)
    int=int+25.0_dp*f7
    call func(ax+dx/2.0_dp*(1.0_dp),ay+dy/2.0_dp*(1.0_dp+sq),c,alpha,f8)
    int=int+40.0_dp*f8
    call func(ax+dx/2.0_dp*(1.0_dp+sq),ay+dy/2.0_dp*(1.0_dp+sq),c,alpha,f9)
    int=int+25.0_dp*f9

    int=int/(4.0_dp*81.0_dp)
    int=int*mesh%cell(k)%dx*mesh%cell(k)%dy

    return
  end subroutine quadrature3_c_alpha

  ! /!\ quadrature3_reconstruction est une moyenne
  subroutine quadrature3_reconstruction(func,mesh,sol,order,normal,k,int)
    procedure(sub_reconstruction) :: func
    type(meshStruct), intent(in) :: mesh
    type(solStruct), intent(in) :: sol
    integer, intent(in) :: order,normal,k
    real(dp), dimension(:), intent(inout) :: int
    real(dp) :: x,y,dx,dy,ax,ay,bx,by,sq
    real(dp), dimension(:), allocatable :: f

    allocate(f(size(int)))

    x=mesh%cell(k)%xc
    y=mesh%cell(k)%yc
    dx=mesh%cell(k)%dx
    dy=mesh%cell(k)%dy
    ax=x-dx/2.0_dp
    ay=y-dy/2.0_dp
    bx=x+dx/2.0_dp
    by=y+dy/2.0_dp
    sq=sqrt(0.6_dp)
    int=0.0_dp

    select case (normal)
    case(1)
       call func(mesh,sol,k,order,ax,ay+dy/2.0_dp*(1.0_dp-sq),f)
       int=int+5.0_dp*f
       call func(mesh,sol,k,order,ax,ay+dy/2.0_dp*(1.0_dp),f)
       int=int+8.0_dp*f
       call func(mesh,sol,k,order,ax,ay+dy/2.0_dp*(1.0_dp+sq),f)
       int=int+5.0_dp*f
    case(2)
       call func(mesh,sol,k,order,ax+dx/2.0_dp*(1.0_dp-sq),ay,f)
       int=5.0_dp*f
       call func(mesh,sol,k,order,ax+dx/2.0_dp*(1.0_dp),ay,f)
       int=int+8.0_dp*f
       call func(mesh,sol,k,order,ax+dx/2.0_dp*(1.0_dp+sq),ay,f)
       int=int+5.0_dp*f
    case(3)
       call func(mesh,sol,k,order,bx,ay+dy/2.0_dp*(1.0_dp-sq),f)
       int=int+5.0_dp*f
       call func(mesh,sol,k,order,bx,ay+dy/2.0_dp*(1.0_dp),f)
       int=int+8.0_dp*f
       call func(mesh,sol,k,order,bx,ay+dy/2.0_dp*(1.0_dp+sq),f)
       int=int+5.0_dp*f
    case(4)
       call func(mesh,sol,k,order,ax+dx/2.0_dp*(1.0_dp-sq),by,f)
       int=int+5.0_dp*f
       call func(mesh,sol,k,order,ax+dx/2.0_dp*(1.0_dp),by,f)
       int=int+8.0_dp*f
       call func(mesh,sol,k,order,ax+dx/2.0_dp*(1.0_dp+sq),by,f)
       int=int+5.0_dp*f       
    end select
    int=int/18.0_dp

    deallocate(f)

    return
  end subroutine quadrature3_reconstruction

end module efficiency
