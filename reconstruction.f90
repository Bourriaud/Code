module reconstruction

  use constant
  use types
  use efficiency

  implicit none

contains
  
  subroutine evaluate(mesh,sol,k,order,gauss_weight,X,Y,U)
    type(meshStruct), intent(in) :: mesh
    type(solStruct), intent(in) :: sol
    integer, intent(in) :: k,order
    real(dp), dimension(:), intent(in) :: gauss_weight,X,Y
    real(dp), dimension(:,:), intent(inout) :: U   !U(p,var)
    integer :: i1,i2,i,d,isol,p
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
    
    do p=1,order
       u(p,:)=sol%val(k,:)
    enddo
    i=1

    do i2=1,d
       do i1=0,d-i2
          alpha(1)=i1
          alpha(2)=i2
          call quad_monome(mesh,c,alpha,k,gauss_weight,intk)
          do isol=1,sol%nvar
             do p=1,order
                u(p,isol)=u(p,isol)+pol(i,isol)*((X(p)-xc)**alpha(1)*(Y(p)-yc)**alpha(2)-intk)
             enddo
          enddo
          i=i+1
       enddo
    enddo

    do i1=1,d
       alpha(1)=i1
       alpha(2)=0
       call quad_monome(mesh,c,alpha,k,gauss_weight,intk)
       do isol=1,sol%nvar
          do p=1,order
             u(p,isol)=u(p,isol)+pol(i,isol)*((X(p)-xc)**alpha(1)*(Y(p)-yc)**alpha(2)-intk)
          enddo
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
    real(dp) :: Kk,Kj,intk,intj
    real(dp), dimension(2) :: c
    integer, dimension(2) :: alpha
    
    d=order-1
    call Nequa(d,N)
    call buildStencil(mesh,k,N,stencil)
    Ni=size(stencil)-1
    Nj=d*(d+1)/2+d

    allocate(X(Ni,Nj),U(Ni,sol%nvar))

    Kk=mesh%cell(k)%dx*mesh%cell(k)%dy
    c(1)=mesh%cell(k)%xc
    c(2)=mesh%cell(k)%yc
    i=1

    do i2=1,d
       do i1=0,d-i2
          alpha(1)=i1
          alpha(2)=i2
          call quad_monome(mesh,c,alpha,k,gauss_weight,intk)
          do j=2,Ni+1
             Kj=mesh%cell(stencil(j))%dx*mesh%cell(stencil(j))%dy
             call quad_monome(mesh,c,alpha,stencil(j),gauss_weight,intj)
             X(j-1,i)=intj-intk
          enddo
          i=i+1
       enddo
    enddo

    do i1=1,d
       alpha(1)=i1
       alpha(2)=0
       call quad_monome(mesh,c,alpha,k,gauss_weight,intk)
       do j=2,Ni+1
          Kj=mesh%cell(stencil(j))%dx*mesh%cell(stencil(j))%dy
          call quad_monome(mesh,c,alpha,stencil(j),gauss_weight,intj)
          X(j-1,i)=intj-intk
       enddo
       i=i+1
    enddo

    do i=1,Ni
       do isol=1,sol%nvar
          U(i,isol)=sol%val(stencil(i+1),isol)-sol%val(k,isol)
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
  
  subroutine buildStencil(mesh,k,N,stencil)
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
    do p1=1,size(mesh%cell(k)%X_gauss)
       do p2=1,size(mesh%cell(k)%Y_gauss)
          call polynomialProduct(mesh%cell(k)%X_gauss(p1),mesh%cell(k)%Y_gauss(p2),c,alpha,s)
          int=int+s*gauss_weight(p1)*gauss_weight(p2)/4.0_dp
       enddo
    enddo

    return
  end subroutine quad_monome
  
end module reconstruction
