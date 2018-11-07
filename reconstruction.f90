module reconstruction

  use constant
  use types
  use efficiency

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

!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!

  subroutine evaluate2(mesh,sol,k,order,x,y,u)
    type(meshStruct), intent(in) :: mesh
    type(solStruct), intent(in) :: sol
    integer, intent(in) :: k,order
    real(dp), intent(in) :: x,y
    real(dp), dimension(:), intent(inout) :: u
    integer :: i1,i2,i,d,isol
    real(dp) :: xc,yc,Kk,intk
    real(dp), dimension(2) :: c
    integer, dimension(2) :: alpha    
    procedure(sub_quadra_c_alpha), pointer :: func
    
    u=sol%val(k,:)
    xc=mesh%cell(k)%xc
    yc=mesh%cell(k)%yc
    c(1)=xc
    c(2)=yc
    func => polynomialProduct
    Kk=mesh%cell(k)%dx*mesh%cell(k)%dy
    d=order-1

    do isol=1,sol%nvar
       i=1
       do i2=1,d
          do i1=0,d-i2
             alpha(1)=i1
             alpha(2)=i2
             call quadrature3(func,mesh,c,alpha,k,intk)
             intk=intk/Kk
             u(isol)=u(isol)+mesh%cell(k)%polCoef(i,isol)*((x-xc)**alpha(1)*(y-yc)**alpha(2)-intk)
             i=i+1
          enddo
       enddo
       
       do i1=1,d
          alpha(1)=i1
          alpha(2)=0
          call quadrature3(func,mesh,c,alpha,k,intk)
          intk=intk/Kk
          u(isol)=u(isol)+mesh%cell(k)%polCoef(i,isol)*((x-xc)**alpha(1)*(y-yc)**alpha(2)-intk)
          i=i+1
       enddo
    enddo

    return
  end subroutine evaluate2
  
  subroutine reconstruct2(mesh,sol,k,order)
    type(meshStruct), intent(inout) :: mesh
    type(solStruct), intent(in) :: sol
    integer, intent(in) :: k,order
    integer, dimension(:), allocatable :: stencil
    real(dp), dimension(:,:), allocatable :: X,U
    integer :: d,N,Ni,Nj,i,i1,i2,j,isol
    real(dp) :: Kk,Kj,intk,intj
    real(dp), dimension(2) :: c
    integer, dimension(2) :: alpha
    procedure(sub_quadra_c_alpha), pointer :: func

    d=order-1
    call Nequa(d,N)
    call buildStencil2(mesh,k,N,stencil)
    Ni=size(stencil)-1
    Nj=d*(d+1)/2+d

    allocate(X(Ni,Nj),mesh%cell(k)%polCoef(Nj,sol%nvar),U(Ni,sol%nvar))

    func => polynomialProduct
    Kk=mesh%cell(k)%dx*mesh%cell(k)%dy
    c(1)=mesh%cell(k)%xc
    c(2)=mesh%cell(k)%yc
    i=1

    do i2=1,d
       do i1=0,d-i2
          alpha(1)=i1
          alpha(2)=i2
          call quadrature3(func,mesh,c,alpha,k,intk)
          intk=intk/Kk
          do j=2,Ni+1
             Kj=mesh%cell(stencil(j))%dx*mesh%cell(stencil(j))%dy
             call quadrature3(func,mesh,c,alpha,stencil(j),intj)
             intj=intj/Kj
             X(j-1,i)=intj-intk
          enddo
          i=i+1
       enddo
    enddo

    do i1=1,d
       alpha(1)=i1
       alpha(2)=0
       call quadrature3(func,mesh,c,alpha,k,intk)
       intk=intk/Kk
       do j=2,Ni+1
          Kj=mesh%cell(stencil(j))%dx*mesh%cell(stencil(j))%dy
          call quadrature3(func,mesh,c,alpha,stencil(j),intj)
          intj=intj/Kj
          X(j-1,i)=intj-intk
       enddo
       i=i+1
    enddo

    do i=1,Ni
       do isol=1,sol%nvar
          U(i,isol)=sol%val(stencil(i+1),isol)-sol%val(k,isol)
       enddo
    enddo
    
    do isol=1,sol%nvar
       call solve(X,U(:,isol),mesh%cell(k)%polCoef(:,isol))
       !call QR(X,U(:,isol),mesh%cell(k)%polCoef(:,isol))
    enddo
    
    deallocate(stencil,X,U)

    return
  end subroutine reconstruct2
  
  subroutine Nequa(d,N)
    integer, intent(in) :: d
    integer, intent(out) :: N

    N=(d+1)*(d+2)/2-1
    
    return
  end subroutine Nequa

  subroutine couronne(mesh,stencil,n)
    type(meshStruct), intent(in) :: mesh
    integer, dimension(:), allocatable, intent(inout) :: stencil
    integer, intent(in) :: n
    integer, dimension(:), allocatable :: stencil2
    integer :: i,j,k,neigh

    allocate(stencil2(2*(2*n+1)**2))
    stencil2=-1
    stencil2(1:size(stencil))=stencil
    k=size(stencil)

    do i=1,size(stencil)
       do j=1,size(mesh%cell(stencil(i))%neigh)
          neigh=mesh%cell(stencil(i))%neigh(j)
          if ((all(stencil2/=neigh)).and.(neigh>0)) then
             k=k+1
             stencil2(k)=neigh
          endif
       enddo
    enddo
    
    deallocate(stencil)
    allocate(stencil(k))

    stencil=stencil2(1:k)

    deallocate(stencil2)
    
    return
  end subroutine couronne
  
  subroutine buildStencil2(mesh,k,N,stencil)
    type(meshStruct), intent(in) :: mesh
    integer, intent(in) :: k,N
    integer, dimension(:), allocatable, intent(out) :: stencil
    integer :: i

    allocate(stencil(1))
    stencil(1)=k
    i=1
    
    do while (size(stencil)<ceiling(1.5*N))
       call couronne(mesh,stencil,i)
       i=i+1
    enddo

    return
  end subroutine buildStencil2

  subroutine polynomialProduct(x,y,c,alpha,s)
    real(dp), intent(in) :: x,y
    real(dp), dimension(2), intent(in) :: c
    integer, dimension(2), intent(in) :: alpha
    real(dp), intent(out) :: s

    s=(x-c(1))**alpha(1)*(y-c(2))**alpha(2)

    return
  end subroutine polynomialProduct

  subroutine solve(X,U,R)
    real(dp), dimension(:,:), intent(in) :: X
    real(dp), dimension(:), intent(in) :: U
    real(dp), dimension(:), intent(inout) :: R
    real(dp), dimension(:,:), allocatable :: A
    real(dp), dimension(:), allocatable :: b
    integer :: i,j,k

    allocate(A(size(R),size(R)),b(size(R)))
    A=0.0_dp
    b=0.0_dp

    do i=1,size(R)
       do j=1,size(R)
          do k=1,size(U)
             A(i,j)=A(i,j)+X(k,i)*X(k,j)
          enddo
       enddo
       do k=1,size(U)
          b(i)=b(i)+X(k,i)*U(k)
       enddo
    enddo

    call cholesky(A,b,R)

    deallocate(A,b)

    return
  end subroutine solve

  subroutine cholesky(A,b,x)
    real(dp), dimension(:,:), intent(in) :: A
    real(dp), dimension(:), intent(in) :: b
    real(dp), dimension(:), intent(inout) :: x
    real(dp), dimension(:,:), allocatable :: L
    integer :: i,j,k,n
    real(dp), dimension(:), allocatable :: x1

    n=size(b)
    allocate(L(n,n),x1(n))
    
    L=0.0_dp
    do j=1,n
       do i=1,j-1
          L(i,j)=A(i,j)
          do k=1,i-1
             L(i,j)=L(i,j)-L(k,i)*L(k,j)
          enddo
          L(i,j)=L(i,j)/L(i,i)
       enddo
       L(j,j)=A(j,j)
       do k=1,j-1
          L(j,j)=L(j,j)-L(k,j)*L(k,j)
       enddo
       L(j,j)=sqrt(L(j,j))
    enddo

    do j=1,n
       x1(j)=b(j)
       do k=1,j-1
          x1(j)=x1(j)-x1(k)*L(k,j)
       enddo
       x1(j)=x1(j)/L(j,j)
    enddo

    do j=1,n
       x(n-j+1)=x1(n-j+1)
       do k=1,j-1
          x(n-j+1)=x(n-j+1)-x(n-k+1)*L(n-j+1,n-k+1)
       enddo
       x(n-j+1)=x(n-j+1)/L(n-j+1,n-j+1)
    enddo

    deallocate(L,x1)
    
    return
  end subroutine cholesky

  subroutine QR(A,b,x)
    real(dp), dimension(:,:), intent(in) :: A
    real(dp), dimension(:), intent(in) :: b
    real(dp), dimension(:), intent(inout) :: x
    real(dp), dimension(:,:), allocatable :: Q,R
    real(dp), dimension(:), allocatable :: x1,x2
    integer :: n,m,i,j,k
    real(dp) :: alpha,beta,c
    real(dp), dimension(:), allocatable :: v

    n=size(A(1,:))
    m=size(A(:,1))
    allocate(Q(m,m),R(m,n),x1(m),x2(n),v(m))

    Q=0.0_dp
    R=A

    do i=1,m
       do j=1,m
          if (i==j) then
             Q(i,j)=1.0_dp
          endif
       enddo
    enddo

    do k=1,n-1
       alpha=0.0_dp
       do i=k,m
          alpha=alpha+R(i,k)**2
       enddo
       alpha=sqrt(alpha)
       alpha=sign(alpha,-R(k,k))
       if (alpha/=0.0_dp) then
          beta=alpha**2-alpha*R(k,k)
          v(k)=R(k,k)-alpha
          do i=k+1,m
             v(i)=R(i,k)
          enddo
          do j=k,n
             c=0.0_dp
             do i=k,m
                c=c+v(i)*R(i,j)
             enddo
             c=c/beta
             do i=k,m
                R(i,j)=R(i,j)-c*v(i)
             enddo
          enddo
          do j=1,m
             c=0.0_dp
             do i=k,m
                c=c+v(i)*Q(i,j)
             enddo
             c=c/beta
             do i=k,m
                Q(i,j)=Q(i,j)-c*v(i)
             enddo
          enddo
       endif
    enddo
    
    !do i=1,m
       !x1(i)=0.0_dp
       !do k=1,m
          !x1(i)=x1(i)+Q(i,k)*b(k)
       !enddo
    !enddo

    !do j=1,n
       !x(n-j+1)=x1(n-j+1)
       !do k=1,j-1
          !x(n-j+1)=x(n-j+1)-x(n-k+1)*R(n-j+1,n-k+1)
       !enddo
       !x(n-j+1)=x(n-j+1)/R(n-j+1,n-j+1)
    !enddo

    do i=1,n
       x1(i)=0.0_dp
       do k=1,m
          x1(i)=x1(i)+A(k,i)*b(k)
       enddo
    enddo

    do j=1,n
       x2(j)=x1(j)
       do k=1,j-1
          x2(j)=x2(j)-x2(k)*R(k,j)
       enddo
       x2(j)=x2(j)/R(j,j)
    enddo

    do j=1,n
       x(n-j+1)=x2(n-j+1)
       do k=1,j-1
          x(n-j+1)=x(n-j+1)-x(n-k+1)*R(n-j+1,n-k+1)
       enddo
       x(n-j+1)=x(n-j+1)/R(n-j+1,n-j+1)
    enddo

    deallocate(Q,R,x1,x2,v)
    
    return
  end subroutine QR
  
end module reconstruction
