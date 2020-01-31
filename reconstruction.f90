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
  
  subroutine reconstruct(mesh,sol,k,order,gauss_weight,period)
    type(meshStruct), intent(inout) :: mesh
    type(solStruct), intent(in) :: sol
    integer, intent(in) :: k,order
    real(dp), dimension(:), intent(in) :: gauss_weight
    logical, intent(in) :: period
    integer, dimension(:), allocatable :: stencil,stencil_type
    real(dp), dimension(:,:), allocatable :: X,U
    integer :: d,N,Ni,Nj,i,isol
    real(dp) :: Kk,pond,dist
    real(dp), dimension(2) :: c
    character(len=99) :: id_char
    integer(16) :: id
    
    d=order-1
    call Nequa(d,N)
    if (period) then
       call buildStencil_period(mesh,k,N,order,stencil,stencil_type)
    else
       call buildStencil(mesh,k,N,order,stencil,stencil_type)
    endif
    Ni=size(stencil)
    Nj=d*(d+1)/2+d
    
    allocate(X(Ni,Nj),U(Ni,sol%nvar))

    Kk=mesh%cell(k)%dx*mesh%cell(k)%dy
    c(1)=mesh%cell(k)%xc
    c(2)=mesh%cell(k)%yc
    pond=1.5_dp

    write(id_char,'(99I1)')order,0,stencil_type
    read(id_char, '(I99)' )id
    
    select case (id)
    case (205555_dp)
       X=X_2_5555
       call adjust_X(mesh%cell(k)%dx,mesh%cell(k)%dy,d,X)
    case (30555555555555_dp)
       X=X_3_555555555555
       call adjust_X(mesh%cell(k)%dx,mesh%cell(k)%dy,d,X)
    case (40555555555555555555555555_16)
       X=X_4_555555555555555555555555
       call adjust_X(mesh%cell(k)%dx,mesh%cell(k)%dy,d,X)
    case default
       call calculate_X(mesh,stencil,gauss_weight,k,c,d,Ni,X)
    end select
    
    !if (k==100.and.order==4) then
       !print*,order
       !print*,stencil_type
       !do i=1,size(X(1,:))
          !print*,X(:,i) !/(mesh%cell(k)%dx)**3
       !enddo
       !print*,"---------------"
       !do i=1,size(X_3_000000000000(:,1))
          !print*,X_3_000000000000(i,:)
       !enddo
       !call exit()
    !endif
    
    do i=1,Ni
       dist=((mesh%cell(stencil(i))%xc-c(1))**2+(mesh%cell(stencil(i))%yc-c(2))**2)**(pond/4.0_dp)
       do isol=1,sol%nvar
          U(i,isol)=(sol%val(stencil(i),isol)-sol%val(k,isol))/dist
       enddo
       X(i,:)=X(i,:)/dist
    enddo
    
    if (allocated(mesh%cell(k)%polCoef)) deallocate(mesh%cell(k)%polCoef)
    allocate(mesh%cell(k)%polCoef(Nj,sol%nvar))
    
    do isol=1,sol%nvar
       call solve(X,U(:,isol),mesh%cell(k)%polCoef(:,isol))
    enddo

    deallocate(stencil,X,U)

    return
  end subroutine reconstruct

  subroutine adjust_X(dx,dy,d,X)
    real(dp), intent(in) :: dx,dy
    integer, intent(in) :: d
    real(dp), dimension(:,:), intent(inout) :: X
    integer :: i1,i2,i

    i=1
    do i2=1,d
       do i1=0,d-i2
          X(:,i)=X(:,i)*(dx**i1)*(dy**i2)
          i=i+1
       enddo
    enddo
    do i1=1,d
       X(:,i)=X(:,i)*(dx**i1)
       i=i+1
    enddo

  end subroutine adjust_X

  subroutine calculate_X(mesh,stencil,gauss_weight,k,c,d,Ni,X)
    type(meshStruct), intent(in) :: mesh
    integer, dimension(:), intent(in) :: stencil
    real(dp), dimension(:), intent(in) :: gauss_weight
    integer, intent(in) :: k,d,Ni
    real(dp), dimension(2), intent(in) :: c
    real(dp), dimension(:,:), intent(inout) :: X
    integer :: i1,i2,i,j
    integer, dimension(2) :: alpha
    real(dp) :: intj,intk,Kj

    i=1
    do i2=1,d
       do i1=0,d-i2
          alpha(1)=i1
          alpha(2)=i2
          call quad_monome(mesh,c,alpha,k,gauss_weight,intk)
          do j=1,Ni
             Kj=mesh%cell(stencil(j))%dx*mesh%cell(stencil(j))%dy
             call quad_monome(mesh,c,alpha,stencil(j),gauss_weight,intj)
             X(j,i)=intj-intk
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
          X(j,i)=intj-intk
       enddo
       i=i+1
    enddo

  end subroutine calculate_X
  
  subroutine Nequa(d,N)
    integer, intent(in) :: d
    integer, intent(out) :: N

    N=(d+1)*(d+2)/2-1

    return
  end subroutine Nequa

  subroutine couronne(mesh,stencil,stencil_type,n)
    type(meshStruct), intent(in) :: mesh
    integer, dimension(:), allocatable, intent(inout) :: stencil,stencil_type
    integer, intent(in) :: n
    integer, dimension(:), allocatable :: stencil2,stencil2_type
    integer :: i,j,k,neigh,s

    s=3*(2*n+1)**2
    allocate(stencil2(s),stencil2_type(s))
    stencil2=-1
    stencil2(1:size(stencil))=stencil
    stencil2_type=0
    stencil2_type(1:size(stencil_type))=stencil_type
    k=size(stencil)

    do i=1,size(stencil)
       do j=1,size(mesh%cell(stencil(i))%edge)
          neigh=mesh%cell(stencil(i))%neigh(j)
          if ((all(stencil2/=neigh)).and.(neigh>0)) then
             k=k+1
             stencil2(k)=neigh
             stencil2_type(k)=mesh%cell(stencil(1))%level-mesh%cell(abs(neigh))%level+5
          endif
       enddo
    enddo
    
    deallocate(stencil,stencil_type)
    allocate(stencil(k),stencil_type(k))

    stencil=stencil2(1:k)
    stencil_type=stencil2_type(1:k)

    deallocate(stencil2,stencil2_type)
    
    return
  end subroutine couronne

  subroutine couronne_period(mesh,stencil,stencil_type,n)
    type(meshStruct), intent(in) :: mesh
    integer, dimension(:), allocatable, intent(inout) :: stencil,stencil_type
    integer, intent(in) :: n
    integer, dimension(:), allocatable :: stencil2,stencil2_type
    integer :: i,j,k,neigh,s

    s=3*(2*n+1)**2
    allocate(stencil2(s),stencil2_type(s))
    stencil2=-1
    stencil2(1:size(stencil))=stencil
    stencil2_type=0
    stencil2_type(1:size(stencil_type))=stencil_type
    k=size(stencil)

    do i=1,size(stencil)
       do j=1,size(mesh%cell(stencil(i))%edge)
          neigh=mesh%cell(stencil(i))%neigh(j)
          if (all(stencil2/=abs(neigh))) then
             k=k+1
             stencil2(k)=abs(neigh)
             stencil2_type(k)=mesh%cell(stencil(1))%level-mesh%cell(abs(neigh))%level+5
          endif
       enddo
    enddo

    deallocate(stencil,stencil_type)
    allocate(stencil(k),stencil_type(k))

    stencil=stencil2(1:k)
    stencil_type=stencil2_type(1:k)

    deallocate(stencil2,stencil2_type)
    
    return
  end subroutine couronne_period
  
  subroutine buildStencil(mesh,k,N,order,stencil,stencil_type)
    type(meshStruct), intent(in) :: mesh
    integer, intent(in) :: k,N,order
    integer, dimension(:), allocatable, intent(out) :: stencil,stencil_type
    integer, dimension(:), allocatable :: stencil2,stencil2_type
    integer :: i,s

    allocate(stencil2(1),stencil2_type(1))
    stencil2(1)=k
    stencil2_type(1)=0
    i=1
    s=ceiling(size_stencil(order)*N)

    do while (size(stencil2)<s)
       call couronne(mesh,stencil2,stencil2_type,i)
       i=i+1
    enddo

    s=size(stencil2)-1
    allocate(stencil(s),stencil_type(s))
    stencil(1:s)=stencil2(2:s+1)
    stencil_type(1:s)=stencil2_type(2:s+1)
    
    deallocate(stencil2,stencil2_type)

    return
  end subroutine buildStencil

  subroutine buildStencil_period(mesh,k,N,order,stencil,stencil_type)
    type(meshStruct), intent(in) :: mesh
    integer, intent(in) :: k,N,order
    integer, dimension(:), allocatable, intent(out) :: stencil,stencil_type
    integer, dimension(:), allocatable :: stencil2,stencil2_type
    integer :: i,s

    allocate(stencil2(1),stencil2_type(1))
    stencil2(1)=k
    stencil2_type(1)=0
    i=1
    s=ceiling(size_stencil(order)*N)
    
    do while (size(stencil2)<s)
       call couronne_period(mesh,stencil2,stencil2_type,i)
       i=i+1
    enddo

    s=size(stencil2)-1
    allocate(stencil(s),stencil_type(s))
    stencil(1:s)=stencil2(2:s+1)
    stencil_type(1:s)=stencil2_type(2:s+1)
    
    deallocate(stencil2,stencil2_type)

    return
  end subroutine buildStencil_period

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
    !integer :: i,j,k

    allocate(A(size(R),size(R)),b(size(R)))
    !A=0.0_dp
    !b=0.0_dp

    !do i=1,size(R)
       !do j=1,size(R)
          !do k=1,size(U)
             !A(i,j)=A(i,j)+X(k,i)*X(k,j)
          !enddo
       !enddo
       !do k=1,size(U)
          !b(i)=b(i)+X(k,i)*U(k)
       !enddo
    !enddo
    A=matmul(transpose(X),X)
    b=matmul(transpose(X),U)
    
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
    real(dp) :: s,s2

    if (size(gauss_weight)==2) then
       int=0.0_dp
       do p1=1,size(gauss_weight)
          s2=0.0_dp
          do p2=1,size(gauss_weight)
             call polynomialProduct(mesh%cell(k)%X_gauss2(p1),mesh%cell(k)%Y_gauss2(p2),c,alpha,s)
             s2=s2+s*gauss_weight(p2)
          enddo
          int=int+gauss_weight(p1)*s2
       enddo
       int=int*0.25_dp
    else
       int=0.0_dp
       do p1=1,size(gauss_weight)
          s2=0.0_dp
          do p2=1,size(gauss_weight)
             call polynomialProduct(mesh%cell(k)%X_gauss(p1),mesh%cell(k)%Y_gauss(p2),c,alpha,s)
             s2=s2+s*gauss_weight(p2)
          enddo
          int=int+gauss_weight(p1)*s2
       enddo
       int=int*0.25_dp
    endif

    return
  end subroutine quad_monome
  
end module reconstruction
