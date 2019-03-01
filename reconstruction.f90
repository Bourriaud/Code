module reconstruction

  use constant
  use types
  use efficiency

  implicit none

contains
  
  subroutine evaluate(mesh,sol,k,order,gauss_weight,x,y,u)
    type(meshStruct), intent(in) :: mesh
    type(solStruct), intent(in) :: sol
    integer, intent(in) :: k,order
    real(dp), dimension(:), intent(in) :: gauss_weight
    real(dp), intent(in) :: x,y
    real(dp), dimension(:), intent(inout) :: u
    integer :: i1,i2,i,d,isol
    real(dp) :: xc,yc,Kk,intk
    real(dp), dimension(2) :: c
    integer, dimension(2) :: alpha    
    real(dp), dimension(:,:), allocatable :: pol

    allocate(pol(size(mesh%cell(k)%polCoef(:,1)),size(mesh%cell(k)%polCoef(1,:))))
    pol=mesh%cell(k)%polCoef

    xc=mesh%cell(k)%xc
    yc=mesh%cell(k)%yc
    c(1)=xc
    c(2)=yc
    Kk=mesh%cell(k)%dx*mesh%cell(k)%dy
    d=order-1
    u(:)=sol%val(k,:)
    i=1

    do i2=1,d
       do i1=0,d-i2
          alpha(1)=i1
          alpha(2)=i2
          call quad_monome(mesh,c,alpha,k,gauss_weight,intk)
          do isol=1,sol%nvar
             u(isol)=u(isol)+pol(i,isol)*((x-xc)**alpha(1)*(y-yc)**alpha(2)-intk)
          enddo
          i=i+1
       enddo
    enddo

    do i1=1,d
       alpha(1)=i1
       alpha(2)=0
       call quad_monome(mesh,c,alpha,k,gauss_weight,intk)
       do isol=1,sol%nvar
          u(isol)=u(isol)+pol(i,isol)*((x-xc)**alpha(1)*(y-yc)**alpha(2)-intk)
       enddo
       i=i+1
    enddo

    deallocate(pol)

    return
  end subroutine evaluate
  
  subroutine reconstruct(mesh,sol,k,order,gauss_weight)
    type(meshStruct), intent(inout) :: mesh
    type(solStruct), intent(in) :: sol
    integer, intent(in) :: k,order
    real(dp), dimension(:), intent(in) :: gauss_weight
    integer, dimension(:), allocatable :: stencil
    real(dp), dimension(:,:), allocatable :: X,U
    integer :: d,N,Ni,Nj,i,i1,i2,j,isol
    real(dp) :: Kk,Kj,intk,intj,pond
    real(dp), dimension(2) :: c
    integer, dimension(2) :: alpha
    
    d=order-1
    call Nequa(d,N)
    call buildStencil(mesh,k,N,order,stencil)
    Ni=size(stencil)
    Nj=d*(d+1)/2+d
    
    allocate(X(Ni,Nj),U(Ni,sol%nvar))

    Kk=mesh%cell(k)%dx*mesh%cell(k)%dy
    c(1)=mesh%cell(k)%xc
    c(2)=mesh%cell(k)%yc
    pond=1.5_dp
    i=1

    do i2=1,d
       do i1=0,d-i2
          alpha(1)=i1
          alpha(2)=i2
          call quad_monome(mesh,c,alpha,k,gauss_weight,intk)
          do j=1,Ni
             Kj=mesh%cell(stencil(j))%dx*mesh%cell(stencil(j))%dy
             call quad_monome(mesh,c,alpha,stencil(j),gauss_weight,intj)
             X(j,i)=(intj-intk)/(((mesh%cell(stencil(j))%xc-c(1))**2+(mesh%cell(stencil(j))%yc-c(2))**2)**(pond/4.0_dp))
          enddo
          i=i+1
       enddo
    enddo

    do i1=1,d
       alpha(1)=i1
       alpha(2)=0
       call quad_monome(mesh,c,alpha,k,gauss_weight,intk)
       do j=1,Ni
          Kj=mesh%cell(stencil(j))%dx*mesh%cell(stencil(j))%dy
          call quad_monome(mesh,c,alpha,stencil(j),gauss_weight,intj)
          X(j,i)=(intj-intk)/(((mesh%cell(stencil(j))%xc-c(1))**2+(mesh%cell(stencil(j))%yc-c(2))**2)**(pond/4.0_dp))
       enddo
       i=i+1
    enddo
    
    do i=1,Ni
       do isol=1,sol%nvar
          U(i,isol)=(sol%val(stencil(i),isol)-sol%val(k,isol))/ &
               (((mesh%cell(stencil(i))%xc-c(1))**2+(mesh%cell(stencil(i))%yc-c(2))**2)**(pond/4.0_dp))
       enddo
    enddo
    
    if (allocated(mesh%cell(k)%polCoef)) deallocate(mesh%cell(k)%polCoef)
    allocate(mesh%cell(k)%polCoef(Nj,sol%nvar))

    do isol=1,sol%nvar
       call solve(X,U(:,isol),mesh%cell(k)%polCoef(:,isol))
    enddo

    deallocate(stencil,X,U)

    return
  end subroutine reconstruct
  
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
  
  subroutine buildStencil(mesh,k,N,order,stencil)
    type(meshStruct), intent(in) :: mesh
    integer, intent(in) :: k,N,order
    integer, dimension(:), allocatable, intent(out) :: stencil
    integer, dimension(:), allocatable :: stencil2
    integer :: i,s

    allocate(stencil2(1))
    stencil2(1)=k
    i=1
    s=ceiling(size_stencil(order)*N)
    
    do while (size(stencil2)<s)
       call couronne(mesh,stencil2,i)
       i=i+1
    enddo

    s=size(stencil2)-1
    allocate(stencil(s))
    stencil(1:s)=stencil2(2:s+1)
    
    deallocate(stencil2)

    return
  end subroutine buildStencil

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
    !call QR(X,b,R)

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
    real(dp), dimension(:,:), allocatable :: R
    real(dp), dimension(:), allocatable :: w,w2
    integer :: k,k2,j,m,n
    real(dp) :: norm,s,u1,tau

    m=size(A(:,1))
    n=size(A(1,:))
    allocate(R(m,n),w(m),w2(n))
    R=A

    do j=1,n
       norm=0
       do k=j,m
          norm=norm+R(k,j)**2
       enddo
       norm=sqrt(norm)
       s=sign(1.0_dp,-R(j,j))
       u1=R(j,j)-s*norm
       w(j:m)=R(j:m,j)/u1
       w(j)=1
       tau=-s*u1/norm
       w2=0.0_dp
       do k2=1,n
          do k=j,m
             w2(k2)=w2(k2)+w(k)*R(k,k2)
          enddo
       enddo
       do k=j,m
          do k2=1,n
             R(k,k2)=R(k,k2)-tau*w(k)*w2(k2)
          enddo
       enddo
    enddo
    
    do j=1,n
       w2(j)=b(j)
       do k=1,j-1
          w2(j)=w2(j)-w2(k)*R(k,j)
       enddo
       w2(j)=w2(j)/R(j,j)
    enddo
    
    do j=1,n
       x(n-j+1)=w2(n-j+1)
       do k=1,j-1
          x(n-j+1)=x(n-j+1)-x(n-k+1)*R(n-j+1,n-k+1)
       enddo
       x(n-j+1)=x(n-j+1)/R(n-j+1,n-j+1)
    enddo

    deallocate(R,w,w2)

    return
  end subroutine QR

  subroutine quad_monome(mesh,c,alpha,k,gauss_weight,int)
    type(meshStruct), intent(in) :: mesh
    real(dp), dimension(2), intent(in) :: c
    integer, dimension(2), intent(in) :: alpha
    integer, intent(in) :: k
    real(dp), dimension(:), intent(in) :: gauss_weight
    real(dp), intent(out) :: int
    integer :: p1,p2
    real(dp) :: s
    
    int=0.0_dp
    do p1=1,size(gauss_weight)
       do p2=1,size(gauss_weight)
          if (size(gauss_weight)==2) then
             call polynomialProduct(mesh%cell(k)%X_gauss2(p1),mesh%cell(k)%Y_gauss2(p2),c,alpha,s)
          else
             call polynomialProduct(mesh%cell(k)%X_gauss(p1),mesh%cell(k)%Y_gauss(p2),c,alpha,s)
          endif
          int=int+s*gauss_weight(p1)*gauss_weight(p2)/4.0_dp
       enddo
    enddo

    return
  end subroutine quad_monome
  
end module reconstruction
