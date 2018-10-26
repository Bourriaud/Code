module FV

  use constant
  use types
  use phys
  use inout
  use efficiency
  
  implicit none
  
contains

  subroutine boundary(flux_ptr,f_equa,neigh,u,sol,boundtype,bound,normal,F,S)
    procedure (sub_flux), pointer, intent(in) :: flux_ptr
    procedure (sub_f), pointer, intent(in) :: f_equa
    integer, intent(in) :: neigh
    real(dp), dimension(:), intent(in) :: u
    type(solStruct), intent(in) :: sol
    character(len=20), intent(in) :: boundtype
    real(dp), dimension(:), intent(in) :: bound
    integer, intent(in) :: normal      !1=left 2=bottom 3=right 4=top
    real(dp), dimension(:), intent(inout) :: F
    real(dp), intent(out) :: S
    real(dp), dimension(:), allocatable :: v

    if (trim(boundtype)=='DIRICHLET') then
       if (normal==1) then
          call flux_ptr(bound,u,f_equa,1,F,S)
       elseif (normal==2) then
          call flux_ptr(bound,u,f_equa,2,F,S)
       elseif (normal==3) then
          call flux_ptr(u,bound,f_equa,1,F,S)
       else
          call flux_ptr(u,bound,f_equa,2,F,S)
       endif
       
    elseif (trim(boundtype)=='TRANSMISSIVE') then
       call flux_ptr(u,u,f_equa,normal,F,S)
       
    elseif (trim(boundtype)=='WALL') then
       allocate(v(size(u)))
       v=u
       if (normal==1) then
          v(2)=-u(2)
          call flux_ptr(v,u,f_equa,1,F,S)
       elseif (normal==2) then
          v(3)=-u(3)
          call flux_ptr(v,u,f_equa,2,F,S)
       elseif (normal==3) then
          v(2)=-u(2)
          call flux_ptr(u,v,f_equa,1,F,S)
       else
          v(3)=-u(3)
          call flux_ptr(u,v,f_equa,2,F,S)
       endif
       deallocate(v)

    elseif (trim(boundtype)=='PERIODIC') then
       allocate(v(size(u)))
       v=sol%val(abs(neigh),:)
       if (normal==1) then
          call flux_ptr(v,u,f_equa,1,F,S)
       elseif (normal==2) then
          call flux_ptr(v,u,f_equa,2,F,S)
       elseif (normal==3) then
          call flux_ptr(u,v,f_equa,1,F,S)
       else
          call flux_ptr(u,v,f_equa,2,F,S)
       endif
       deallocate(v)
       
    else
       print*,"Boundary condition ",trim(boundtype)," not implemented"
       call exit()
    endif

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

    if ((normal==1).or.(normal)==3) then
       dir=1
    else
       dir=2
    endif
    
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

    if ((normal==1).or.(normal)==3) then
       dir=1
    else
       dir=2
    endif

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
