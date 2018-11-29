module FV

  use constant
  use types
  use phys
  use inout
  use efficiency
  use reconstruction
  
  implicit none
  
contains

  subroutine boundary(flux,f_equa,gauss_weight,cell,neigh,u,mesh,sol,edge,p,boundtype,bound,normal,order,F)
    procedure (sub_flux), pointer, intent(in) :: flux
    procedure (sub_f), pointer, intent(in) :: f_equa
    real(dp), dimension(:), intent(in) :: gauss_weight
    integer, intent(in) :: cell,neigh,edge,p
    real(dp), dimension(:), intent(in) :: u
    type(meshStruct), intent(inout) :: mesh
    type(solStruct), intent(in) :: sol
    character(len=20), intent(in) :: boundtype
    real(dp), dimension(:), intent(in) :: bound
    integer, intent(in) :: normal      !1=left 2=bottom 3=right 4=top
    integer, intent(in) :: order
    real(dp), dimension(:), intent(inout) :: F
    real(dp), dimension(:), allocatable :: v
    integer :: k,period,isol
    
    select case (trim(boundtype))
    case ('DIRICHLET')
       select case (normal)
       case (1)
          call flux(bound,u,f_equa,1,F)
       case (2)
          call flux(bound,u,f_equa,2,F)
       case (3)
          call flux(u,bound,f_equa,1,F)
       case (4)
          call flux(u,bound,f_equa,2,F)
       end select
       
    case ('NEUMANN')
       allocate(v(sol%nvar))
       select case (normal)
       case (1)
          call evaluate(mesh,sol,cell,order,gauss_weight,mesh%cell(cell)%xc,mesh%edge(edge)%Y_gauss(p),v)
          do isol=1,sol%nvar
             v(isol)=sol%val(cell,isol)-bound(isol)*mesh%cell(cell)%dx/2.0_dp
          enddo
          call flux(v,v,f_equa,1,F)
       case (2)
          call evaluate(mesh,sol,cell,order,gauss_weight,mesh%edge(edge)%X_gauss(p),mesh%cell(cell)%yc,v)
          do isol=1,sol%nvar
             v(isol)=sol%val(cell,isol)-bound(isol)*mesh%cell(cell)%dy/2.0_dp
          enddo
          call flux(v,v,f_equa,2,F)
       case (3)
          call evaluate(mesh,sol,cell,order,gauss_weight,mesh%cell(cell)%xc,mesh%edge(edge)%Y_gauss(p),v)
          do isol=1,sol%nvar
             v(isol)=sol%val(cell,isol)+bound(isol)*mesh%cell(cell)%dx/2.0_dp
          enddo
          call flux(v,v,f_equa,1,F)
       case (4)
          call evaluate(mesh,sol,cell,order,gauss_weight,mesh%edge(edge)%X_gauss(p),mesh%cell(cell)%yc,v)
          do isol=1,sol%nvar
             v(isol)=sol%val(cell,isol)+bound(isol)*mesh%cell(cell)%dy/2.0_dp
          enddo
          call flux(v,v,f_equa,2,F)
       end select
       
       deallocate(v)          
       
    case ('TRANSMISSIVE')
       select case (normal)
       case (1)
          call flux(u,u,f_equa,1,F)
       case (2)
          call flux(u,u,f_equa,2,F)
       case (3)
          call flux(u,u,f_equa,1,F)
       case (4)
          call flux(u,u,f_equa,2,F)
       end select
       
    case ('WALL')
       allocate(v(sol%nvar))
       v=u
       select case (normal)
       case (1)
          v(2)=-u(2)
          call flux(v,u,f_equa,1,F)
       case (2)
          v(3)=-u(3)
          call flux(v,u,f_equa,2,F)
       case (3)
          v(2)=-u(2)
          call flux(u,v,f_equa,1,F)
       case (4)
          v(3)=-u(3)
          call flux(u,v,f_equa,2,F)
       end select
       deallocate(v)

    case ('PERIODIC')
       allocate(v(sol%nvar))
       k=abs(neigh)
       period=mesh%edge(edge)%period
       call evaluate(mesh,sol,k,order,gauss_weight,mesh%edge(period)%X_gauss(p),mesh%edge(period)%Y_gauss(p),v)
       select case (normal)
       case (1)
          call flux(v,u,f_equa,1,F)
       case (2)
          call flux(v,u,f_equa,2,F)
       case (3)
          call flux(u,v,f_equa,1,F)
       case (4)
          call flux(u,v,f_equa,2,F)
       end select
       deallocate(v)
       
    case default
       print*,"Boundary condition ",trim(boundtype)," not implemented"
       call exit()
    end select
    
  end subroutine boundary

  subroutine RS_advection(u1,u2,f_equa,ustar,dir)
    real(dp), dimension(:), intent(in) :: u1,u2
    procedure (sub_f), pointer, intent(in) :: f_equa
    real(dp), dimension(:), intent(inout) :: ustar
    integer, intent(in) :: dir
    real(dp) :: s
    real(dp), dimension(:,:), allocatable :: F1vect,F2vect
    integer :: i

    allocate(F1vect(size(u1),2),F2vect(size(u1),2))
    
    call f_equa(u1,F1vect(:,:))
    call f_equa(u2,F2vect(:,:))
    
    do i=1,size(u1)
       if (u1(i)/=u2(i)) then
          s=(F2vect(i,dir)-F1vect(i,dir))/(u2(i)-u1(i))
          if (s<0.0_dp) then
             ustar(i)=u2(i)
          else
             ustar(i)=u1(i)
          endif
       else
          ustar(i)=u1(i)
       endif
    enddo

    deallocate(F1vect,F2vect)

    return
  end subroutine RS_advection
  
  subroutine flux_godunov(u1,u2,f_equa,dir,F)
    real(dp), dimension(:), intent(in) :: u1,u2
    procedure (sub_f), pointer, intent(in) :: f_equa
    integer, intent(in) :: dir     !1=vertical 2=horizontal
    real(dp), dimension(:), intent(inout) :: F
    real(dp), dimension(:), allocatable :: ustar
    real(dp), dimension(:,:), allocatable :: Fvect

    allocate (ustar(size(u1)),Fvect(size(u1),2))
    
    call RS_advection(u1,u2,f_equa,ustar,dir)
    call f_equa(ustar,Fvect)

    F(:)=Fvect(:,dir)

    deallocate(ustar,Fvect)
    
    return
  end subroutine flux_godunov

  subroutine speed_godunov(u1,u2,f_equa,Smax)
    real(dp), dimension(:), intent(in) :: u1,u2
    procedure (sub_f), pointer, intent(in) :: f_equa
    real(dp), dimension(2), intent(out) :: Smax
    real(dp) :: s
    real(dp), dimension(:,:), allocatable :: f1vect,f2vect
    integer :: i

    allocate(F1vect(size(u1),2),F2vect(size(u1),2))

    Smax=0.0_dp
    
    call f_equa(u1,F1vect(:,:))
    call f_equa(u2,F2vect(:,:))
    
    do i=1,size(u1)
       if (u1(i)/=u2(i)) then
          s=(F2vect(i,1)-F1vect(i,1))/(u2(i)-u1(i))
          Smax(1)=max(abs(s),Smax(1))
          s=(F2vect(i,2)-F1vect(i,2))/(u2(i)-u1(i))
          Smax(2)=max(abs(s),Smax(2))
       endif
    enddo

    deallocate(F1vect,F2vect)

    return
  end subroutine speed_godunov

  subroutine flux_HLL(u1,u2,f_equa,dir,F)
    real(dp), dimension(:), intent(in) :: u1,u2
    procedure (sub_f), pointer, intent(in) :: f_equa
    integer, intent(in) :: dir     !1=vertical 2=horizontal
    real(dp), dimension(:), intent(inout) :: F
    real(dp), dimension(:,:), allocatable :: F1vect,F2vect
    real(dp) :: SL,SR,p1,p2,a1,a2

    allocate (F1vect(size(u1),2),F2vect(size(u2),2))

    call unconserv(u1,4,p1)
    call unconserv(u2,4,p2)
    a1=sqrt(gamma*abs(p1/u1(1)))
    a2=sqrt(gamma*abs(p2/u2(1)))
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
    
    deallocate(F1vect,F2vect)

    return
  end subroutine flux_HLL

  subroutine speed_HLL(u1,u2,f_equa,Smax)
    use constant
    real(dp), dimension(:), intent(in) :: u1,u2
    procedure (sub_f), pointer, intent(in) :: f_equa
    real(dp), dimension(2), intent(out) :: Smax
    real(dp) :: SL,SR,p1,p2,a1,a2

    call unconserv(u1,4,p1)
    call unconserv(u2,4,p2)
    a1=sqrt(gamma*abs(p1/u1(1)))
    a2=sqrt(gamma*abs(p2/u2(1)))
    SL=min(u1(2)/u1(1)-a1,u2(2)/u2(1)-a2)
    SR=max(u1(2)/u1(1)+a1,u2(2)/u2(1)+a2)
    
    Smax(1)=max(abs(SL),abs(SR))

    SL=min(u1(3)/u1(1)-a1,u2(3)/u2(1)-a2)
    SR=max(u1(3)/u1(1)+a1,u2(3)/u2(1)+a2)

    Smax(2)=max(abs(SL),abs(SR))

    return
  end subroutine speed_HLL
       

end module FV
