module FV

  use constant
  use types
  use phys
  use inout
  use efficiency
  use reconstruction
  
  implicit none
  
contains

  subroutine boundary(flux_ptr,f_equa,neigh,u,mesh,sol,boundtype,bound,normal,order,F,S)
    procedure (sub_flux), pointer, intent(in) :: flux_ptr
    procedure (sub_f), pointer, intent(in) :: f_equa
    integer, intent(in) :: neigh
    real(dp), dimension(:), intent(in) :: u
    type(meshStruct), intent(in) :: mesh
    type(solStruct), intent(in) :: sol
    character(len=20), intent(in) :: boundtype
    real(dp), dimension(:), intent(in) :: bound
    integer, intent(in) :: normal      !1=left 2=bottom 3=right 4=top
    integer, intent(in) :: order
    real(dp), dimension(:), intent(inout) :: F
    real(dp), intent(out) :: S
    real(dp), dimension(:), allocatable :: v
    integer :: k
    procedure (sub_reconstruction), pointer :: func

    select case (trim(boundtype))
    case ('DIRICHLET')
       select case (normal)
       case (1)
          call flux_ptr(bound,u,f_equa,1,F,S)
       case (2)
          call flux_ptr(bound,u,f_equa,2,F,S)
       case (3)
          call flux_ptr(u,bound,f_equa,1,F,S)
       case (4)
          call flux_ptr(u,bound,f_equa,2,F,S)
       end select
       
    case ('TRANSMISSIVE')
       call flux_ptr(u,u,f_equa,normal,F,S)
       
    case ('WALL')
       allocate(v(size(u)))
       v=u
       select case (normal)
       case (1)
          v(2)=-u(2)
          call flux_ptr(v,u,f_equa,1,F,S)
       case (2)
          v(3)=-u(3)
          call flux_ptr(v,u,f_equa,2,F,S)
       case (3)
          v(2)=-u(2)
          call flux_ptr(u,v,f_equa,1,F,S)
       case (4)
          v(3)=-u(3)
          call flux_ptr(u,v,f_equa,2,F,S)
       end select
       deallocate(v)

    case ('PERIODIC')
       allocate(v(size(u)))
       func => evaluate2
       k=abs(neigh)
       select case (normal)
       case (1)
          call quadrature3(func,mesh,sol,order,3,k,v)
          call flux_ptr(v,u,f_equa,1,F,S)
       case (2)
          call quadrature3(func,mesh,sol,order,4,k,v)
          call flux_ptr(v,u,f_equa,2,F,S)
       case (3)
          call quadrature3(func,mesh,sol,order,1,k,v)
          call flux_ptr(u,v,f_equa,1,F,S)
       case (4)
          call quadrature3(func,mesh,sol,order,2,k,v)
          call flux_ptr(u,v,f_equa,2,F,S)
       end select
       deallocate(v)
       
    case default
       print*,"Boundary condition ",trim(boundtype)," not implemented"
       call exit()
    end select

  end subroutine boundary
  
  subroutine flux_godunov(u1,u2,f_equa,normal,F,Smax)
    real(dp), dimension(:), intent(in) :: u1,u2
    procedure (sub_f), pointer, intent(in) :: f_equa
    integer, intent(in) :: normal     !1=left 2=bottom 3=right 4=top
    real(dp), dimension(:), intent(inout) :: F
    real(dp), intent(out) :: Smax
    real(dp), dimension(:), allocatable :: ustar
    real(dp), dimension(:,:), allocatable :: Fvect
    integer :: dir

    allocate (ustar(size(u1)),Fvect(size(u1),2))

    select case (normal)
    case (1)
       dir=1
    case (2)
       dir=2
    case (3)
       dir=1
    case (4)
       dir=2
    end select
    
    call RS_advection(u1,u2,f_equa,ustar,dir,Smax)
    call f_equa(ustar,Fvect)

    F(:)=Fvect(:,dir)

    deallocate(ustar,Fvect)
    
    return
  end subroutine flux_godunov

  subroutine RS_advection(u1,u2,f_equa,ustar,dir,Smax)
    real(dp), dimension(:), intent(in) :: u1,u2
    procedure (sub_f), pointer, intent(in) :: f_equa
    real(dp), dimension(:), intent(inout) :: ustar
    integer, intent(in) :: dir
    real(dp), intent(out) :: Smax
    real(dp) :: s
    real(dp), dimension(:,:), allocatable :: f1vect,f2vect
    real(dp), dimension(:), allocatable :: f1,f2
    integer :: i

    allocate(f1vect(size(u1),2),f2vect(size(u1),2),f1(size(u1)),f2(size(u1)))

    Smax=0.
    
    call f_equa(u1,f1vect(:,:))
    call f_equa(u2,f2vect(:,:))
    f1=f1vect(:,dir)
    f2=f2vect(:,dir)
    
    do i=1,size(u1)
       if (u1(i)/=u2(i)) then
          s=(f2(i)-f1(i))/(u2(i)-u1(i))
          Smax=max(abs(s),Smax)
          if (s<0.0_dp) then
             ustar(i)=u2(i)
          else
             ustar(i)=u1(i)
          endif
       else
          ustar(i)=u1(i)
       endif
    enddo

    return
  end subroutine RS_advection

  subroutine flux_HLL(u1,u2,f_equa,normal,F,Smax)
    real(dp), dimension(:), intent(in) :: u1,u2
    procedure (sub_f), pointer, intent(in) :: f_equa
    integer, intent(in) :: normal     !1=left 2=bottom 3=right 4=top
    real(dp), dimension(:), intent(inout) :: F
    real(dp), intent(out) :: Smax
    real(dp), dimension(:,:), allocatable :: F1vect,F2vect
    integer :: dir
    real(dp) :: SL,SR,p1,p2,a1,a2,gamma

    allocate (F1vect(size(u1),2),F2vect(size(u2),2))

    select case (normal)
    case (1)
       dir=1
    case (2)
       dir=2
    case (3)
       dir=1
    case (4)
       dir=2
    end select

    gamma=1.4
    p1=(u1(4)-u1(1)*(0.5_dp*((u1(2)/u1(1))**2+(u1(3)/u1(1))**2)))*(gamma-1)
    p2=(u2(4)-u2(1)*(0.5_dp*((u2(2)/u2(1))**2+(u2(3)/u2(1))**2)))*(gamma-1)
    a1=sqrt(gamma*p1/u1(1))
    a2=sqrt(gamma*p2/u2(1))
    SL=min(u1(1+dir)/u1(1)-a1,u2(1+dir)/u2(1)-a2)
    SR=max(u1(1+dir)/u1(1)+a1,u2(1+dir)/u2(1)+a2)
    
    if (SL>0.0_dp) then
       call f_equa(u1,F1vect)
       F(:)=F1vect(:,dir)
    else if (SR<0.0_dp) then
       call f_equa(u2,F2vect)
       F(:)=F2vect(:,dir)
    else
       call f_equa(u1,F1vect)
       call f_equa(u2,F2vect)
       F(:)=(SR*F1vect(:,dir)-SL*F2vect(:,dir)+SL*SR*(u2(:)-u1(:)))/(SR-SL)
    endif
    Smax=max(abs(SL),abs(SR))
    
    deallocate(F1vect,F2vect)
    
    return
  end subroutine flux_HLL
       

end module FV
