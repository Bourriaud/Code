module reconstruction

  use constant
  use types
  use efficiency

  implicit none

contains

  subroutine evaluate_cell(mesh,sol,k,xc,yc,dx,dy,gauss_weight,gauss_point,u)
    type(meshStruct), intent(in) :: mesh
    type(solStruct), intent(in) :: sol
    integer, intent(in) :: k
    real(dp), intent(in) :: xc,yc,dx,dy
    real(dp), dimension(:), intent(in) :: gauss_weight,gauss_point
    real(dp), dimension(:), intent(inout) :: u
    integer :: p1,p2,order
    real(dp), dimension(:), allocatable :: utemp
    real(dp), dimension(:,:), allocatable :: pol2

    allocate(utemp(size(u)))

    call reconstruct1(mesh,sol,k,1,gauss_weight,(/.true.,.true./),pol2)
    order=size(gauss_point)
    
    u=0.0_dp
    do p1=1,order
       do p2=1,order
          call evaluate(mesh,sol,pol2,k,xc+gauss_point(p1)*dx/2.0_dp,yc+gauss_point(p2)*dy/2.0_dp,utemp)
          u=u+utemp*gauss_weight(p1)*gauss_weight(p2)/4.0_dp
       enddo
    enddo

    deallocate(utemp)

    return
  end subroutine evaluate_cell
  
  subroutine evaluate(mesh,sol,pol,k,x,y,u)
    type(meshStruct), intent(in) :: mesh
    type(solStruct), intent(in) :: sol
    real(dp), dimension(:,:), intent(in) :: pol
    integer, intent(in) :: k
    real(dp), intent(in) :: x,y
    real(dp), dimension(:), intent(inout) :: u
    integer :: i1,i2,i,d,isol,order
    real(dp) :: xc,yc,Kk,intk,temp
    integer, dimension(2) :: alpha
    
    xc=mesh%cell(k)%xc
    yc=mesh%cell(k)%yc
    Kk=mesh%cell(k)%dx*mesh%cell(k)%dy
    u(:)=sol%val(k,:)

    select case (size(pol(:,1)))
    case (0)
       order=1
    case (2)
       order=2
    case (5)
       order=3
    case (9)
       order=4
    case (14)
       order=5
    case default
       print*,"Order not available in evaluate"
       call exit()
    end select
    d=order-1
    i=1

    do i2=1,d
       do i1=0,d-i2
          alpha(1)=i1
          alpha(2)=i2
          call quad_monome2(mesh,alpha,k,intk)
          temp=(x-xc)**alpha(1)*(y-yc)**alpha(2)-intk
          do isol=1,sol%nvar
             u(isol)=u(isol)+pol(i,isol)*temp
          enddo
          i=i+1
       enddo
    enddo

    do i1=1,d
       alpha(1)=i1
       alpha(2)=0
       call quad_monome2(mesh,alpha,k,intk)
       temp=(x-xc)**alpha(1)*(y-yc)**alpha(2)-intk
       do isol=1,sol%nvar
          u(isol)=u(isol)+pol(i,isol)*temp
       enddo
       i=i+1
    enddo
    
    return
  end subroutine evaluate

  subroutine reconstruct(mesh,sol,k,order,gauss_weight,period,pol_OUT)
    type(meshStruct), intent(inout) :: mesh
    type(solStruct), intent(in) :: sol
    integer, intent(in) :: k,order
    real(dp), dimension(:), intent(in) :: gauss_weight
    real(dp), dimension(:), allocatable :: gauss_point
    logical, dimension(2), intent(in) :: period
    real(dp), dimension(:,:), allocatable, intent(inout) :: pol_OUT

    allocate (gauss_point(size(gauss_weight)))
    select case (order)
    case (1)
       gauss_point=gauss_point1
    case (2)
       gauss_point=gauss_point2
    case (3)
       gauss_point=gauss_point3
    case (4)
       gauss_point=gauss_point4
    case (5)
       gauss_point=gauss_point5
    end select
    
    if (mesh%cell(k)%bound<=2) then
       if (mesh%cell(k)%bound==0.or.period(max(mesh%cell(k)%bound,1))) then
          select case (order)
          case (4)
             call reconstruct1(mesh,sol,k,order,gauss_weight,period,pol_OUT)
             !call reconstruct2(mesh,sol,k,order,pol_OUT)
          case default
             !call reconstruct1(mesh,sol,k,order,gauss_weight,period,pol_OUT)
             call reconstruct2(mesh,sol,k,order,pol_OUT)
          end select
       else
          call reconstruct1(mesh,sol,k,order,gauss_weight,period,pol_OUT)
       endif
    else
       if (period(1).and.period(2)) then
          select case (order)
          case (4)
             call reconstruct1(mesh,sol,k,order,gauss_weight,period,pol_OUT)
             !call reconstruct2(mesh,sol,k,order,pol_OUT)
          case default
             !call reconstruct1(mesh,sol,k,order,gauss_weight,period,pol_OUT)
             call reconstruct2(mesh,sol,k,order,pol_OUT)
          end select
       else
          call reconstruct1(mesh,sol,k,order,gauss_weight,period,pol_OUT)
       endif
    endif

    deallocate (gauss_point)

    return
  end subroutine reconstruct
  
  subroutine reconstruct1(mesh,sol,k,order,gauss_weight,period,pol)
    type(meshStruct), intent(in) :: mesh
    type(solStruct), intent(in) :: sol
    integer, intent(in) :: k,order
    real(dp), dimension(:), intent(in) :: gauss_weight
    logical, dimension(2), intent(in) :: period
    real(dp), dimension(:,:), allocatable, intent(inout) :: pol
    integer, dimension(:), allocatable :: stencil,stencil_type
    real(dp), dimension(:,:), allocatable :: stencil_bound
    real(dp), dimension(:,:), allocatable :: X,U
    integer :: d,N,Ni,Nj,i,isol
    real(dp) :: Kk,pond,dist
    real(dp), dimension(2) :: c
    character(len=99) :: id_char
    integer(16) :: id
    
    d=order-1
    call Nequa(d,N)
    call buildStencil(mesh,k,N,order,period,stencil,stencil_type,stencil_bound)
    
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
    case (805555_dp)
       X=X_2_5555
       call adjust_X(mesh%cell(k)%dx,mesh%cell(k)%dy,d,X)
    case (80555555555555_dp)
       X=X_3_555555555555
       call adjust_X(mesh%cell(k)%dx,mesh%cell(k)%dy,d,X)
    case (80555555555555555555555555_16)
       X=X_4_555555555555555555555555
       call adjust_X(mesh%cell(k)%dx,mesh%cell(k)%dy,d,X)
    case default
       call calculate_X(mesh,stencil,stencil_bound,gauss_weight,k,c,d,Ni,X)
    end select

    do i=1,Ni
       dist=((mesh%cell(stencil(i))%xc-c(1)+stencil_bound(i,1))**2+ &
            (mesh%cell(stencil(i))%yc-c(2)+stencil_bound(i,2))**2)**(pond/4.0_dp)
       do isol=1,sol%nvar
          U(i,isol)=(sol%val(stencil(i),isol)-sol%val(k,isol))/dist
       enddo
       X(i,:)=X(i,:)/dist
    enddo

    if (allocated(pol)) deallocate(pol)
    allocate(pol(Nj,sol%nvar))
    
    do isol=1,sol%nvar
       call solve(X,U(:,isol),pol(:,isol))
    enddo

    deallocate(stencil,X,U)

    return
  end subroutine reconstruct1

  subroutine reconstruct2(mesh,sol,k,order,pol_OUT)
    type(meshStruct), intent(inout) :: mesh
    type(solStruct), intent(in) :: sol
    integer, intent(in) :: k,order
    real(dp), dimension(:,:), allocatable, intent(inout) :: pol_OUT
    integer, dimension(:), allocatable :: stencil
    real(dp), dimension(:,:), allocatable :: X,U,pol
    integer :: d,Ni,Nj,isol
    real(dp) :: pond,dist,dx,dy,dx2,dy2
    real(dp), dimension(2) :: c

    Ni=16
    d=order-1
    Nj=d*(d+1)/2+d
    
    allocate(X(Ni,Nj),U(Ni,sol%nvar),stencil(12),pol(0,0))

    stencil=abs(mesh%cell(k)%stencil2)
    c(1)=mesh%cell(k)%xc
    c(2)=mesh%cell(k)%yc
    dx=mesh%cell(k)%dx
    dy=mesh%cell(k)%dy
    pond=2.0_dp

    select case (order)
    !case (4)
       !X=X_4
       !call adjust_X(dx,dy,d,X)
    case (3)
       X=X_3
       call adjust_X(dx,dy,d,X)
    case (2)
       X=X_2
       call adjust_X(dx,dx,d,X)
    case default
       !call calculate_X(mesh,stencil,mesh%cell(k)%stencil2_bound,gauss_weight,k,c,d,Ni,X)
       if(allocated(pol_OUT))deallocate(pol_OUT)
       allocate(pol_OUT(Nj,sol%nvar))
       pol_OUT(:,:)=0.0_dp
       return
    end select
    
    dx2=0.5*dx
    dy2=0.5*dy
    !dist=((dx)**2+(dy)**2)**(pond/2.0_dp)
    dist=0.1_dp
    call evaluate(mesh,sol,pol,stencil(1), &
         c(1)-1.0_dp*dx-mesh%cell(k)%stencil2_bound(1,1),c(2)+1.0_dp*dy-mesh%cell(k)%stencil2_bound(1,2),U(1,:))
    U(1,:)=(U(1,:)-sol%val(k,:))/dist
    X(1,:)=X(1,:)/dist
    call evaluate(mesh,sol,pol,stencil(4), &
         c(1)-1.0_dp*dx-mesh%cell(k)%stencil2_bound(4,1),c(2)-1.0_dp*dy-mesh%cell(k)%stencil2_bound(4,2),U(5,:))
    U(5,:)=(U(5,:)-sol%val(k,:))/dist
    X(5,:)=X(5,:)/dist
    call evaluate(mesh,sol,pol,stencil(7), &
         c(1)+1.0_dp*dx-mesh%cell(k)%stencil2_bound(7,1),c(2)-1.0_dp*dy-mesh%cell(k)%stencil2_bound(7,2),U(9,:))
    U(9,:)=(U(9,:)-sol%val(k,:))/dist
    X(9,:)=X(9,:)/dist
    call evaluate(mesh,sol,pol,stencil(10), &
         c(1)+1.0_dp*dx-mesh%cell(k)%stencil2_bound(10,1),c(2)+1.0_dp*dy-mesh%cell(k)%stencil2_bound(10,2),U(13,:))
    U(13,:)=(U(13,:)-sol%val(k,:))/dist
    X(13,:)=X(13,:)/dist

    !dist=((1.0_dp*dx)**2+(0.5*dy)**2)**(pond/2.0_dp)
    dist=10.0_dp
    call evaluate(mesh,sol,pol,stencil(2), &
         c(1)-1.0_dp*dx-mesh%cell(k)%stencil2_bound(2,1),c(2)+0.5_dp*dy-mesh%cell(k)%stencil2_bound(2,2),U(2,:))
    U(2,:)=(U(2,:)-sol%val(k,:))/dist
    X(2,:)=X(2,:)/dist
    call evaluate(mesh,sol,pol,stencil(3), &
         c(1)-1.0_dp*dx-mesh%cell(k)%stencil2_bound(3,1),c(2)-0.5_dp*dy-mesh%cell(k)%stencil2_bound(3,2),U(4,:))
    U(4,:)=(U(4,:)-sol%val(k,:))/dist
    X(4,:)=X(4,:)/dist
    call evaluate(mesh,sol,pol,stencil(8), &
         c(1)+1.0_dp*dx-mesh%cell(k)%stencil2_bound(8,1),c(2)-0.5_dp*dy-mesh%cell(k)%stencil2_bound(8,2),U(10,:))
    U(10,:)=(U(10,:)-sol%val(k,:))/dist
    X(10,:)=X(10,:)/dist
    call evaluate(mesh,sol,pol,stencil(9), &
         c(1)+1.0_dp*dx-mesh%cell(k)%stencil2_bound(9,1),c(2)+0.5_dp*dy-mesh%cell(k)%stencil2_bound(9,2),U(12,:))
    U(12,:)=(U(12,:)-sol%val(k,:))/dist
    X(12,:)=X(12,:)/dist

    !dist=((0.5_dp*dx)**2+(1.0*dy)**2)**(pond/2.0_dp)
    dist=10.0_dp
    call evaluate(mesh,sol,pol,stencil(5), &
         c(1)-0.5_dp*dx-mesh%cell(k)%stencil2_bound(5,1),c(2)-1.0_dp*dy-mesh%cell(k)%stencil2_bound(5,2),U(6,:))
    U(6,:)=(U(6,:)-sol%val(k,:))/dist
    X(6,:)=X(6,:)/dist
    call evaluate(mesh,sol,pol,stencil(6), &
         c(1)+0.5_dp*dx-mesh%cell(k)%stencil2_bound(6,1),c(2)-1.0_dp*dy-mesh%cell(k)%stencil2_bound(6,2),U(8,:))
    U(8,:)=(U(8,:)-sol%val(k,:))/dist
    X(8,:)=X(8,:)/dist
    call evaluate(mesh,sol,pol,stencil(11), &
         c(1)+0.5_dp*dx-mesh%cell(k)%stencil2_bound(11,1),c(2)+1.0_dp*dy-mesh%cell(k)%stencil2_bound(11,2),U(14,:))
    U(14,:)=(U(14,:)-sol%val(k,:))/dist
    X(14,:)=X(14,:)/dist
    call evaluate(mesh,sol,pol,stencil(12), &
         c(1)-0.5_dp*dx-mesh%cell(k)%stencil2_bound(12,1),c(2)+1.0_dp*dy-mesh%cell(k)%stencil2_bound(12,2),U(16,:))
    U(16,:)=(U(16,:)-sol%val(k,:))/dist
    X(16,:)=X(16,:)/dist

    dist=1.0_dp
    call evaluate(mesh,sol,pol,stencil(2), &
         c(1)-1.0_dp*dx-mesh%cell(k)%stencil2_bound(2,1),c(2)-0.0_dp*dy-mesh%cell(k)%stencil2_bound(2,2),U(3,:))
    U(3,:)=(U(3,:)-sol%val(k,:))/dist
    X(3,:)=X(3,:)/dist
    call evaluate(mesh,sol,pol,stencil(5), &
         c(1)+0.0_dp*dx-mesh%cell(k)%stencil2_bound(5,1),c(2)-1.0_dp*dy-mesh%cell(k)%stencil2_bound(5,2),U(7,:))
    U(7,:)=(U(7,:)-sol%val(k,:))/dist
    X(7,:)=X(7,:)/dist
    call evaluate(mesh,sol,pol,stencil(8), &
         c(1)+1.0_dp*dx-mesh%cell(k)%stencil2_bound(8,1),c(2)+0.0_dp*dy-mesh%cell(k)%stencil2_bound(8,2),U(11,:))
    U(11,:)=(U(11,:)-sol%val(k,:))/dist
    X(11,:)=X(11,:)/dist
    call evaluate(mesh,sol,pol,stencil(11), &
         c(1)-0.0_dp*dx-mesh%cell(k)%stencil2_bound(11,1),c(2)+1.0_dp*dy-mesh%cell(k)%stencil2_bound(11,2),U(15,:))
    U(15,:)=(U(15,:)-sol%val(k,:))/dist
    X(15,:)=X(15,:)/dist

    if (allocated(pol_OUT)) deallocate(pol_OUT)
    allocate(pol_OUT(Nj,sol%nvar))
    
    do isol=1,sol%nvar
       call solve(X,U(:,isol),pol_OUT(:,isol))
    enddo

    deallocate(stencil,X,U)
    
    return
  end subroutine reconstruct2

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

  subroutine calculate_X(mesh,stencil,stencil_bound,gauss_weight,k,c,d,Ni,X)
    type(meshStruct), intent(in) :: mesh
    integer, dimension(:), intent(in) :: stencil
    real(dp), dimension(:,:), intent(in) :: stencil_bound
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
          call quad_monome2(mesh,alpha,k,intk)
          do j=1,Ni
             Kj=mesh%cell(stencil(j))%dx*mesh%cell(stencil(j))%dy
             call quad_monome(mesh,c-stencil_bound(j,:),alpha,stencil(j),gauss_weight,intj)
             X(j,i)=intj-intk
          enddo
          i=i+1
       enddo
    enddo
    do i1=1,d
       alpha(1)=i1
       alpha(2)=0
       call quad_monome2(mesh,alpha,k,intk)
       do j=1,Ni
          Kj=mesh%cell(stencil(j))%dx*mesh%cell(stencil(j))%dy
          call quad_monome(mesh,c-stencil_bound(j,:),alpha,stencil(j),gauss_weight,intj)
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

  subroutine couronne(mesh,period,stencil,stencil_type,n)
    type(meshStruct), intent(in) :: mesh
    logical, dimension(2), intent(in) :: period
    integer, dimension(:), allocatable, intent(inout) :: stencil,stencil_type
    integer, intent(in) :: n
    type(cellStruct) :: cell
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
       cell=mesh%cell(stencil(i))
       do j=1,size(cell%edge)
          neigh=cell%neigh(j)
          if (neigh>0.or.period(mesh%edge(cell%edge(j))%dir)) then
             if (all(stencil2/=abs(neigh))) then
                k=k+1
                stencil2(k)=abs(neigh)
                stencil2_type(k)=mesh%cell(stencil(1))%level-mesh%cell(abs(neigh))%level+5
             endif
          endif
       enddo

       neigh=cell%corner_cell(1)
       if (neigh>0) then
          if (all(stencil2/=abs(neigh))) then
             k=k+1
             stencil2(k)=abs(neigh)
             stencil2_type(k)=mesh%cell(stencil(1))%level-mesh%cell(abs(neigh))%level+5
          endif
       else if (period(mesh%edge(cell%edge(1))%dir).and.mesh%cell(abs(neigh))%yc<mesh%cell(stencil(i))%yc) then
          if (all(stencil2/=abs(neigh))) then
             k=k+1
             stencil2(k)=abs(neigh)
             stencil2_type(k)=mesh%cell(stencil(1))%level-mesh%cell(abs(neigh))%level+5
          endif
       else if (period(mesh%edge(cell%edge(3))%dir).and.mesh%cell(abs(neigh))%xc<mesh%cell(stencil(i))%xc) then
          if (all(stencil2/=abs(neigh))) then
             k=k+1
             stencil2(k)=abs(neigh)
             stencil2_type(k)=mesh%cell(stencil(1))%level-mesh%cell(abs(neigh))%level+5
          endif 
       endif

       neigh=cell%corner_cell(2)
       if (neigh>0) then
          if (all(stencil2/=abs(neigh))) then
             k=k+1
             stencil2(k)=abs(neigh)
             stencil2_type(k)=mesh%cell(stencil(1))%level-mesh%cell(abs(neigh))%level+5
          endif
       else if (period(mesh%edge(cell%edge(2))%dir).and.mesh%cell(abs(neigh))%yc<mesh%cell(stencil(i))%yc) then
          if (all(stencil2/=abs(neigh))) then
             k=k+1
             stencil2(k)=abs(neigh)
             stencil2_type(k)=mesh%cell(stencil(1))%level-mesh%cell(abs(neigh))%level+5
          endif
       else if (period(mesh%edge(cell%edge(3))%dir).and.mesh%cell(abs(neigh))%xc>mesh%cell(stencil(i))%xc) then
          if (all(stencil2/=abs(neigh))) then
             k=k+1
             stencil2(k)=abs(neigh)
             stencil2_type(k)=mesh%cell(stencil(1))%level-mesh%cell(abs(neigh))%level+5
          endif 
       endif

       neigh=cell%corner_cell(3)
       if (neigh>0) then
          if (all(stencil2/=abs(neigh))) then
             k=k+1
             stencil2(k)=abs(neigh)
             stencil2_type(k)=mesh%cell(stencil(1))%level-mesh%cell(abs(neigh))%level+5
          endif
       else if (period(mesh%edge(cell%edge(1))%dir).and.mesh%cell(abs(neigh))%yc>mesh%cell(stencil(i))%yc) then
          if (all(stencil2/=abs(neigh))) then
             k=k+1
             stencil2(k)=abs(neigh)
             stencil2_type(k)=mesh%cell(stencil(1))%level-mesh%cell(abs(neigh))%level+5
          endif
       else if (period(mesh%edge(cell%edge(4))%dir).and.mesh%cell(abs(neigh))%xc<mesh%cell(stencil(i))%xc) then
          if (all(stencil2/=abs(neigh))) then
             k=k+1
             stencil2(k)=abs(neigh)
             stencil2_type(k)=mesh%cell(stencil(1))%level-mesh%cell(abs(neigh))%level+5
          endif 
       endif

       neigh=cell%corner_cell(4)
       if (neigh>0) then
          if (all(stencil2/=abs(neigh))) then
             k=k+1
             stencil2(k)=abs(neigh)
             stencil2_type(k)=mesh%cell(stencil(1))%level-mesh%cell(abs(neigh))%level+5
          endif
       else if (period(mesh%edge(cell%edge(2))%dir).and.mesh%cell(abs(neigh))%yc>mesh%cell(stencil(i))%yc) then
          if (all(stencil2/=abs(neigh))) then
             k=k+1
             stencil2(k)=abs(neigh)
             stencil2_type(k)=mesh%cell(stencil(1))%level-mesh%cell(abs(neigh))%level+5
          endif
       else if (period(mesh%edge(cell%edge(4))%dir).and.mesh%cell(abs(neigh))%xc>mesh%cell(stencil(i))%xc) then
          if (all(stencil2/=abs(neigh))) then
             k=k+1
             stencil2(k)=abs(neigh)
             stencil2_type(k)=mesh%cell(stencil(1))%level-mesh%cell(abs(neigh))%level+5
          endif 
       endif

    enddo

    deallocate(stencil,stencil_type)
    allocate(stencil(k),stencil_type(k))

    stencil=stencil2(1:k)
    stencil_type=stencil2_type(1:k)

    deallocate(stencil2,stencil2_type)
    
    return
  end subroutine couronne

  subroutine buildStencil(mesh,k,N,order,period,stencil,stencil_type,stencil_bound)
    type(meshStruct), intent(in) :: mesh
    integer, intent(in) :: k,N,order
    logical, dimension(2), intent(in) :: period
    integer, dimension(:), allocatable, intent(out) :: stencil,stencil_type
    real(dp), dimension(:,:), allocatable, intent(out) :: stencil_bound
    integer, dimension(:), allocatable :: stencil2,stencil2_type
    integer :: i,s

    allocate(stencil2(1),stencil2_type(1))
    stencil2(1)=k
    stencil2_type(1)=0
    i=1
    s=ceiling(size_stencil(order)*N)
    
    do while (size(stencil2)<s)
       call couronne(mesh,period,stencil2,stencil2_type,i)
       i=i+1
    enddo

    s=size(stencil2)-1
    allocate(stencil(s),stencil_type(s),stencil_bound(s,2))
    stencil(1:s)=stencil2(2:s+1)
    stencil_type(1:s)=stencil2_type(2:s+1)
    stencil_bound=0.0_dp
    do i=1,s
       if (mesh%cell(k)%xc-mesh%cell(abs(stencil(i)))%xc>mesh%Lx/2.0_dp) stencil_bound(i,1)=stencil_bound(i,1)+mesh%Lx
       if (mesh%cell(k)%xc-mesh%cell(abs(stencil(i)))%xc<-mesh%Lx/2.0_dp) stencil_bound(i,1)=stencil_bound(i,1)-mesh%Lx
       if (mesh%cell(k)%yc-mesh%cell(abs(stencil(i)))%yc>mesh%Ly/2.0_dp) stencil_bound(i,2)=stencil_bound(i,2)+mesh%Ly
       if (mesh%cell(k)%yc-mesh%cell(abs(stencil(i)))%yc<-mesh%Ly/2.0_dp) stencil_bound(i,2)=stencil_bound(i,2)-mesh%Ly
    enddo
    
    deallocate(stencil2,stencil2_type)

    return
  end subroutine buildStencil

  subroutine buildStencil2(mesh,k)
    type(meshStruct), intent(inout) :: mesh
    integer, intent(in) :: k
    integer :: i

    i=1

    if (mesh%cell(k)%level<mesh%cell(abs(mesh%cell(k)%neigh(i)))%level) then
       mesh%cell(k)%stencil2(2)=mesh%cell(k)%neigh(i+1)
       mesh%cell(k)%stencil2(3)=mesh%cell(k)%neigh(i)
       i=i+2
    else
       mesh%cell(k)%stencil2(2)=mesh%cell(k)%neigh(i)
       mesh%cell(k)%stencil2(3)=mesh%cell(k)%neigh(i)
       i=i+1
    endif

    if (mesh%cell(k)%level<mesh%cell(abs(mesh%cell(k)%neigh(i)))%level) then
       mesh%cell(k)%stencil2(8)=mesh%cell(k)%neigh(i)
       mesh%cell(k)%stencil2(9)=mesh%cell(k)%neigh(i+1)
       i=i+2
    else
       mesh%cell(k)%stencil2(8)=mesh%cell(k)%neigh(i)
       mesh%cell(k)%stencil2(9)=mesh%cell(k)%neigh(i)
       i=i+1
    endif

    if (mesh%cell(k)%level<mesh%cell(abs(mesh%cell(k)%neigh(i)))%level) then
       mesh%cell(k)%stencil2(5)=mesh%cell(k)%neigh(i)
       mesh%cell(k)%stencil2(6)=mesh%cell(k)%neigh(i+1)
       i=i+2
    else
       mesh%cell(k)%stencil2(5)=mesh%cell(k)%neigh(i)
       mesh%cell(k)%stencil2(6)=mesh%cell(k)%neigh(i)
       i=i+1
    endif

    if (mesh%cell(k)%level<mesh%cell(abs(mesh%cell(k)%neigh(i)))%level) then
       mesh%cell(k)%stencil2(11)=mesh%cell(k)%neigh(i+1)
       mesh%cell(k)%stencil2(12)=mesh%cell(k)%neigh(i)
       i=i+2
    else
       mesh%cell(k)%stencil2(11)=mesh%cell(k)%neigh(i)
       mesh%cell(k)%stencil2(12)=mesh%cell(k)%neigh(i)
       i=i+1
    endif

    mesh%cell(k)%stencil2(1)=mesh%cell(k)%corner_cell(3)
    mesh%cell(k)%stencil2(4)=mesh%cell(k)%corner_cell(1)
    mesh%cell(k)%stencil2(7)=mesh%cell(k)%corner_cell(2)
    mesh%cell(k)%stencil2(10)=mesh%cell(k)%corner_cell(4)

    mesh%cell(k)%stencil2_bound=0
    do i=1,12
       if (mesh%cell(k)%xc-mesh%cell(abs(mesh%cell(k)%stencil2(i)))%xc>mesh%Lx/2.0_dp) then
          mesh%cell(k)%stencil2_bound(i,1)=mesh%cell(k)%stencil2_bound(i,1)+mesh%Lx
       endif
       if (mesh%cell(k)%xc-mesh%cell(abs(mesh%cell(k)%stencil2(i)))%xc<-mesh%Lx/2.0_dp)  then
          mesh%cell(k)%stencil2_bound(i,1)=mesh%cell(k)%stencil2_bound(i,1)-mesh%Lx
       endif
       if (mesh%cell(k)%yc-mesh%cell(abs(mesh%cell(k)%stencil2(i)))%yc>mesh%Ly/2.0_dp) then
          mesh%cell(k)%stencil2_bound(i,2)=mesh%cell(k)%stencil2_bound(i,2)+mesh%Ly
       endif
       if (mesh%cell(k)%yc-mesh%cell(abs(mesh%cell(k)%stencil2(i)))%yc<-mesh%Ly/2.0_dp) then
          mesh%cell(k)%stencil2_bound(i,2)=mesh%cell(k)%stencil2_bound(i,2)-mesh%Ly
       endif
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

   subroutine quad_monome2(mesh,alpha,k,int)
    type(meshStruct), intent(in) :: mesh
    integer, dimension(2), intent(in) :: alpha
    integer, intent(in) :: k
    real(dp), intent(out) :: int

    if (mod(alpha(1),2)==0.and.mod(alpha(2),2)==0) then
       int=mesh%cell(k)%dx**(alpha(1))*mesh%cell(k)%dy**(alpha(2))/((alpha(1)+1)*(alpha(2)+1)*2**(alpha(1)+alpha(2)))
    else
       int=0.0_dp
    endif

    return
  end subroutine quad_monome2
  
end module reconstruction
