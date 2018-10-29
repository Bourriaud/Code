module reconstruction

  use constant
  use types

  implicit none

contains

  recursive subroutine DivDiff(X,U,i1,i2,dd)
    real(dp), dimension(:), intent(in) :: X,U
    integer, intent(in) :: i1,i2
    real(dp), intent(out) :: dd
    real(dp) :: dd1,dd2

    if (i2==i1+1) then
       dd=U(i1)
    else
       call DivDiff(X,U,i1,i2-1,dd1)
       call DivDiff(X,U,i1+1,i2,dd2)
       dd=(dd2-dd1)/(X(i2)-X(i1))
    endif

    return
  end subroutine DivDiff

  subroutine evaluate(x,deg,stencil,U,val)
    real(dp), intent(in) :: x
    integer, intent(in) :: deg
    real(dp), dimension(:), intent(in) :: stencil,U
    real(dp), intent(out) :: val
    integer :: j,m,l
    real(dp) :: prod,sum,dd

    val=0.0_dp
    do j=1,deg
       sum=0.0_dp
       do m=0,j-1
          prod=1.0_dp
          do l=0,m-1
             prod=prod*(x-stencil(l+1))
          enddo
          do l=m+1,j-1
             prod=prod*(x-stencil(l+1))
          enddo
          sum=sum+prod
       enddo
       call DivDiff(stencil,U,1,j+1,dd)
       val=val+dd*sum
    enddo

    return
  end subroutine evaluate

  subroutine buildStencil(X,U,k,order,stencil,Ux)
    real(dp), dimension(:), intent(in) :: X,U
    integer, intent(in) :: k,order
    real(dp), dimension(:), allocatable, intent(out) :: stencil,Ux
    integer :: i,nl,nr
    real(dp) :: dd1,dd2

    if (k+1<=ubound(X,1)) then
       nl=k
       nr=k+1
    else
       nl=k-1
       nr=k
    endif

    do i=2,order
       if (nl-1<=0) then
          nr=nr+1
       elseif (nr+1>size(X)) then
          nl=nl-1
       else
          call DivDiff(X,U,nl-1,nr,dd1)
          call DivDiff(X,U,nl,nr+1,dd2)
          if (abs(dd1) < abs(dd2)) then
             nl=nl-1
          else
             nr=nr+1   
          endif
       endif
    enddo

    allocate(stencil(nr-nl+1))
    stencil(1:nr-nl+1)=X(nl:nr)

    allocate(Ux(nr-nl))
    Ux(1:nr-nl)=U(nl:nr-1)   

    return
  end subroutine buildStencil

  subroutine extractDirection(mesh,sol,k,order,normal,isol,X,U,kpos)
    type(meshStruct), intent(in) :: mesh
    type(solStruct), intent(in) :: sol
    integer, intent(in) :: k,order,normal,isol
    real(dp), dimension(:), allocatable, intent(out) :: X,U
    integer, intent(out) :: kpos
    real(dp), dimension(:), allocatable :: X2,U2
    integer :: i,j,neigh1,neigh2
    
    allocate(X2(k-2*order:k+2*order),U2(k-2*order:k+2*order-1))

    if ((normal==1).or.(normal==3)) then
       X2(k)=mesh%cell(k)%xc-mesh%cell(k)%dx/2.0_dp
       neigh1=mesh%cell(k)%edge(1)%neigh
       neigh2=k
       i=1
       do while ((i<=2*order).and.(neigh1>-1))
          X2(k-i)=X2(k-i+1)-mesh%cell(neigh1)%dx
          U2(k-i)=sol%val(neigh1,isol)
          neigh1=mesh%cell(neigh1)%edge(1)%neigh
          i=i+1
       enddo
       j=1
       do while ((j<=2*order).and.(neigh2>-1))
          X2(k+j)=X2(k+j-1)+mesh%cell(neigh2)%dx
          U2(k+j-1)=sol%val(neigh2,isol)
          neigh2=mesh%cell(neigh2)%edge(3)%neigh
          j=j+1
       enddo
    else
       X2(k)=mesh%cell(k)%yc-mesh%cell(k)%dy/2.
       neigh1=mesh%cell(k)%edge(2)%neigh
       neigh2=k
       i=1
       do while ((i<=2*order).and.(neigh1>-1))
          X2(k-i)=X2(k-i+1)-mesh%cell(neigh1)%dy
          U2(k-i)=sol%val(neigh1,isol)
          neigh1=mesh%cell(neigh1)%edge(2)%neigh
          i=i+1
       enddo
       j=1
       do while ((j<=2*order).and.(neigh2>-1))
          X2(k+j)=X2(k+j-1)+mesh%cell(neigh2)%dy
          U2(k+j-1)=sol%val(neigh2,isol)
          neigh2=mesh%cell(neigh2)%edge(4)%neigh
          j=j+1
       enddo
    endif
    
    kpos=i

    allocate(X(1:i+j-1),U(1:i+j-2))
    X(:)=X2(k-i+1:k+j-1)
    U(:)=U2(k-i+1:k+j-2)

    deallocate(X2,U2)
    
    return
  end subroutine extractDirection
          
  
  subroutine reconstruct(mesh,sol,k,neigh,order,normal,ul,ur)
    type(meshStruct), intent(in) :: mesh
    type(solStruct), intent(in) :: sol
    integer, intent(in) :: k,neigh,order,normal
    real(dp), dimension(:), intent(inout) :: ul,ur
    real(dp), dimension(:), allocatable :: X,U,Xstencil,Ustencil
    integer :: kpos,isol
    real(dp) :: xbound

    select case (normal)
    case (1)
       xbound=mesh%cell(k)%xc-mesh%cell(k)%dx/2.0_dp
    case (2)
       xbound=mesh%cell(k)%yc-mesh%cell(k)%dy/2.0_dp
    case (3)
       xbound=mesh%cell(k)%xc+mesh%cell(k)%dx/2.0_dp
    case (4)
       xbound=mesh%cell(k)%yc+mesh%cell(k)%dy/2.0_dp
    end select

    do isol=1,sol%nvar
       if ((normal==3).or.(normal==4)) then
          call extractDirection(mesh,sol,k,order,normal,isol,X,U,kpos)
          call buildStencil(X,U,kpos,order,Xstencil,Ustencil)
          call evaluate(xbound,order,Xstencil,Ustencil,ul(isol))
          call extractDirection(mesh,sol,neigh,order,normal,isol,X,U,kpos)
          call buildStencil(X,U,kpos,order,Xstencil,Ustencil)
          call evaluate(xbound,order,Xstencil,Ustencil,ur(isol))
       else
          call extractDirection(mesh,sol,k,order,normal,isol,X,U,kpos)
          call buildStencil(X,U,kpos,order,Xstencil,Ustencil)
          call evaluate(xbound,order,Xstencil,Ustencil,ur(isol))
          call extractDirection(mesh,sol,neigh,order,normal,isol,X,U,kpos)
          call buildStencil(X,U,kpos,order,Xstencil,Ustencil)
          call evaluate(xbound,order,Xstencil,Ustencil,ul(isol))
       endif
    enddo
    
    deallocate(X,U,Xstencil,Ustencil)

    return
  end subroutine reconstruct

  subroutine reconstruct_boundary(mesh,sol,k,order,normal,ubound)
    type(meshStruct), intent(in) :: mesh
    type(solStruct), intent(in) :: sol
    integer, intent(in) :: k,order,normal
    real(dp), dimension(:), intent(inout) :: ubound
    real(dp), dimension(:), allocatable :: X,U,Xstencil,Ustencil
    integer :: kpos,isol
    real(dp) :: xbound

    select case (normal)
    case (1)
       xbound=mesh%cell(k)%xc-mesh%cell(k)%dx/2.0_dp
    case (2)
       xbound=mesh%cell(k)%yc-mesh%cell(k)%dy/2.0_dp
    case (3)
       xbound=mesh%cell(k)%xc+mesh%cell(k)%dx/2.0_dp
    case (4)
       xbound=mesh%cell(k)%yc+mesh%cell(k)%dy/2.0_dp
    end select

    do isol=1,sol%nvar
       call extractDirection(mesh,sol,k,order,normal,isol,X,U,kpos)
       call buildStencil(X,U,kpos,order,Xstencil,Ustencil)
       call evaluate(xbound,order,Xstencil,Ustencil,ubound(isol))
    enddo
    
    deallocate(X,U,Xstencil,Ustencil)

    return
  end subroutine reconstruct_boundary
    
  
end module reconstruction
