module efficiency

  use constant
  use types
  
  implicit none

contains

  subroutine exactSol(x,y,t,s)
    real(dp), intent(in) :: x,y,t
    real(dp), intent(out) :: s
    real(dp) :: a1,a2

    a1=-1.0_dp
    a2=-1.0_dp

    s=cos((x-a1*t-5.0_dp)*pi/5.0_dp)+cos((y-a2*t-5.0_dp)*pi/5.0_dp)

    return
  end subroutine exactSol
  
  subroutine exactTab(t,mesh,tab)
    real(dp), intent(in) :: t
    type (meshStruct) :: mesh
    real(dp), dimension(:,:), intent(inout) :: tab
    integer :: k
    real(dp) :: x,y

    do k=1,mesh%nc
       x=mesh%cell(k)%xc
       y=mesh%cell(k)%yc
       call quadrature3(exactSol,mesh,t,k,tab(k,1))
    enddo
    
    return
  end subroutine exactTab

  subroutine userSol(t,mesh,sol)
    real(dp), intent(in) :: t
    type(meshStruct), intent(in) :: mesh
    type(solStruct), intent(inout) :: sol

    call exactTab(t,mesh,sol%user(:,1:1))

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

  subroutine quadrature1(func,mesh,t,k,int)
    procedure(sub_fquadra) :: func
    type(meshStruct), intent(in) :: mesh
    real(dp), intent(in) :: t
    integer, intent(in) :: k
    real(dp), intent(out) :: int
    real(dp) :: x,y

    x=mesh%cell(k)%xc
    y=mesh%cell(k)%yc
    call func(x,y,t,int)

    return
  end subroutine quadrature1

  subroutine quadrature3(func,mesh,t,k,int)
    procedure(sub_fquadra) :: func
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

    return
  end subroutine quadrature3

end module efficiency
