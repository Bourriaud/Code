module reconstruction

  use types

  implicit none

contains

  recursive subroutine DivDiff(X,U,i1,i2,dd)
    real, dimension(:), intent(in) :: X,U
    integer, intent(in) :: i1,i2
    real, intent(out) :: dd
    real :: dd1,dd2

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
    real, intent(in) :: x
    integer, intent(in) :: deg
    real, dimension(:), intent(in) :: stencil,U
    real, intent(out) :: val
    integer :: j,m,l
    real :: prod,sum,dd

    val=0.
    do j=1,deg
       sum=0.
       do m=0,j-1
          prod=1.
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
    real, dimension(:), intent(in) :: X,U
    integer, intent(in) :: k,order
    real, dimension(:), allocatable, intent(out) :: stencil,Ux
    integer :: i,nl,nr
    real :: dd1,dd2

    if (k+1<size(X)) then
       nl=k
       nr=k+1
    else
       nl=k-1
       nr=k
    endif

    do i=2,order
       call DivDiff(X,U,nl-1,nr,dd1)
       call DivDiff(X,U,nl,nr+1,dd2)
       if (((abs(dd1) < abs(dd2)).and.(nl-1>0)).or.(nr+1>size(X))) then
          nl=nl-1
       else
          nr=nr+1   
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
    real, dimension(:), allocatable, intent(out) :: X,U
    integer, intent(out) :: kpos
    real, dimension(:), allocatable :: X2,U2
    integer :: i,j,neigh1,neigh2

    allocate(X2(k-2*order:k+2*order),U2(k-2*order:k+2*order-1))

    if ((normal==1).or.(normal==3)) then
       X2(k)=mesh%cell(k)%xc-mesh%cell(k)%dx/2.
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
    real, dimension(:), intent(inout) :: ul,ur
    real, dimension(:), allocatable :: X,U,Xstencil,Ustencil
    integer :: kpos,isol
    real :: xbound

    if (normal==1) then
       xbound=mesh%cell(k)%xc-mesh%cell(k)%dx/2.
    elseif (normal==2) then
       xbound=mesh%cell(k)%yc-mesh%cell(k)%dy/2.
    elseif (normal==3) then
       xbound=mesh%cell(k)%xc+mesh%cell(k)%dx/2.
    else
       xbound=mesh%cell(k)%yc+mesh%cell(k)%dy/2.
    endif

    do isol=1,sol%nvar
       if (k<neigh) then
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

    !print*,"xbound : ",mesh%cell(k)%xc,xbound
    !print*,"Ugauche : ",U(kpos-1),"Uk : ",U(kpos),"Udroite : ",U(kpos+1)
    !print*,"    ul : ",ul,"    ur : ",ur
    !print*,"Bon u  : ",xbound**10,xbound**1
    
    deallocate(X,U,Xstencil,Ustencil)

    return
  end subroutine reconstruct

  subroutine reconstruct_boundary(mesh,sol,k,order,normal,ubound)
    type(meshStruct), intent(in) :: mesh
    type(solStruct), intent(in) :: sol
    integer, intent(in) :: k,order,normal
    real, dimension(:), intent(inout) :: ubound
    real, dimension(:), allocatable :: X,U,Xstencil,Ustencil
    integer :: kpos,isol
    real :: xbound

    if (normal==1) then
       xbound=mesh%cell(k)%xc-mesh%cell(k)%dx/2.
    elseif (normal==2) then
       xbound=mesh%cell(k)%yc-mesh%cell(k)%dy/2.
    elseif (normal==3) then
       xbound=mesh%cell(k)%xc+mesh%cell(k)%dx/2.
    else
       xbound=mesh%cell(k)%yc+mesh%cell(k)%dy/2.
    endif

    do isol=1,sol%nvar
       call extractDirection(mesh,sol,k,order,normal,isol,X,U,kpos)
       call buildStencil(X,U,kpos,order,Xstencil,Ustencil)
       call evaluate(xbound,order,Xstencil,Ustencil,ubound(isol))
    enddo

    !print*,"xbound : ",mesh%cell(k)%xc,xbound
    !print*,"Ugauche : ",U(kpos-1),"Uk : ",U(kpos),"Udroite : ",U(kpos+1)
    !print*,"    ubound : ",ubound
    !print*,"Bon ubound  : ",xbound**10,xbound**1
    
    deallocate(X,U,Xstencil,Ustencil)

    return
  end subroutine reconstruct_boundary
    
  
end module reconstruction
