module time

  use constant
  use types
  use FV
  use reconstruction
  
  implicit none

contains

  subroutine advance(mesh,sol,sol2,f_ptr,flux_ptr,order,cfl,t,tf,quad_t,quad_c_alpha,quad_reconstruct)
    type(meshStruct), intent(inout) :: mesh
    type(solStruct), intent(in) :: sol
    type(solStruct), intent(inout) :: sol2
    procedure (sub_f), pointer, intent(in) :: f_ptr
    procedure (sub_flux), pointer, intent(in) :: flux_ptr
    integer, intent(in) :: order
    real(dp), intent(in) :: cfl,tf
    real(dp), intent(inout) :: t
    procedure (quadrature_t), pointer, intent(in) :: quad_t
    procedure (quadrature_c_alpha), pointer, intent(in) :: quad_c_alpha
    procedure (quadrature_reconstruction), pointer, intent(in) :: quad_reconstruct
    integer :: k,i,neigh,normal
    real(dp), dimension(:,:,:), allocatable :: F   !F(k,var,dir)
    real(dp), dimension(:), allocatable :: Ftemp,u1,u2
    type(cellStruct) :: cell
    type(edgeStruct) :: edge
    real(dp) :: xi,yi,dx,dy,dt,S
    procedure (sub_reconstruction), pointer :: func

    allocate(F(mesh%nc,sol%nvar,4),Ftemp(sol%nvar),u1(sol%nvar),u2(sol%nvar))

    dt=10.0_dp**20
    F=0.0_dp
    func => evaluate

    do k=1,mesh%nc
       call reconstruct(mesh,sol,k,order,quad_c_alpha)
    enddo
    
    do k=1,mesh%nc
       cell=mesh%cell(k)
       dx=cell%dx
       dy=cell%dy
       do i=1,size(cell%edge)
          edge=cell%edge(i)
          neigh=edge%neigh
          normal=edge%normal
          xi=(mesh%node(cell%edge(i)%node1)%x+mesh%node(cell%edge(i)%node2)%x)/2.0_dp
          yi=(mesh%node(cell%edge(i)%node1)%y+mesh%node(cell%edge(i)%node2)%y)/2.0_dp
          call quad_reconstruct(func,mesh,sol,order,quad_c_alpha,normal,k,u1)
          if (neigh<0) then  
             call boundary(flux_ptr,f_ptr,quad_c_alpha,quad_reconstruct,neigh, &
                  u1,mesh,sol,edge%boundType,edge%bound,normal,order,Ftemp(:),S)
          else
             if ((normal==1).or.(normal==2)) then
                call quad_reconstruct(func,mesh,sol,order,quad_c_alpha,normal+2,neigh,u2)
                call flux_ptr(u2,u1,f_ptr,normal,Ftemp(:),S)
             else
                call quad_reconstruct(func,mesh,sol,order,quad_c_alpha,normal-2,neigh,u2)
                call flux_ptr(u1,u2,f_ptr,normal,Ftemp(:),S)
             endif
          endif
          F(k,:,normal)=F(k,:,normal)+Ftemp(:)
          dt=min(dt,cfl*min(dx,dy)/S)
       enddo
    enddo

    if (t+dt>tf) then
       dt=tf-t
    endif
    
    sol2%val(1:mesh%nc,:)=sol%val(1:mesh%nc,:)-dt/dx*(F(1:mesh%nc,:,3)-F(1:mesh%nc,:,1))-dt/dy*(F(1:mesh%nc,:,4)-F(1:mesh%nc,:,2))
    t=t+dt

    do k=1,mesh%nc
       deallocate(mesh%cell(k)%polCoef)
    enddo
    
    deallocate(F,Ftemp,u1,u2)
    
    return
  end subroutine advance

  subroutine euler_exp(mesh,sol,f_ptr,flux_ptr,order,cfl,t,tf,quad_t,quad_c_alpha,quad_reconstruct)
    type(meshStruct), intent(inout) :: mesh
    type(solStruct), intent(inout) :: sol
    procedure (sub_f), pointer, intent(in) :: f_ptr
    procedure (sub_flux), pointer, intent(in) :: flux_ptr
    integer, intent(in) :: order
    real(dp), intent(in) :: cfl,tf
    real(dp), intent(inout) :: t
    procedure (quadrature_t), pointer, intent(in) :: quad_t
    procedure (quadrature_c_alpha), pointer, intent(in) :: quad_c_alpha
    procedure (quadrature_reconstruction), pointer, intent(in) :: quad_reconstruct
    type(solStruct) :: sol1

    allocate(sol1%val(mesh%nc,sol%nvar))
    sol1=sol

    call advance(mesh,sol,sol1,f_ptr,flux_ptr,order,cfl,t,tf,quad_t,quad_c_alpha,quad_reconstruct)
    sol%val=sol1%val

    deallocate(sol1%val)

    return
  end subroutine euler_exp

  subroutine SSPRK2(mesh,sol,f_ptr,flux_ptr,order,cfl,t,tf,quad_t,quad_c_alpha,quad_reconstruct)
    type(meshStruct), intent(inout) :: mesh
    type(solStruct), intent(inout) :: sol
    procedure (sub_f), pointer, intent(in) :: f_ptr
    procedure (sub_flux), pointer, intent(in) :: flux_ptr
    integer, intent(in) :: order
    real(dp), intent(in) :: cfl,tf
    real(dp), intent(inout) :: t
    procedure (quadrature_t), pointer, intent(in) :: quad_t
    procedure (quadrature_c_alpha), pointer, intent(in) :: quad_c_alpha
    procedure (quadrature_reconstruction), pointer, intent(in) :: quad_reconstruct
    real(dp) :: t1
    type(solStruct) :: sol1,sol2

    allocate(sol1%val(mesh%nc,sol%nvar),sol2%val(mesh%nc,sol%nvar))
    sol1=sol
    sol2=sol

    t1=t
    call advance(mesh,sol,sol1,f_ptr,flux_ptr,order,cfl,t1,tf,quad_t,quad_c_alpha,quad_reconstruct)
    call advance(mesh,sol1,sol2,f_ptr,flux_ptr,order,cfl,t,tf,quad_t,quad_c_alpha,quad_reconstruct)
    sol%val=0.5_dp*(sol%val+sol2%val)

    deallocate(sol1%val,sol2%val)

    return
  end subroutine SSPRK2

  subroutine SSPRK3(mesh,sol,f_ptr,flux_ptr,order,cfl,t,tf,quad_t,quad_c_alpha,quad_reconstruct)
    type(meshStruct), intent(inout) :: mesh
    type(solStruct), intent(inout) :: sol
    procedure (sub_f), pointer, intent(in) :: f_ptr
    procedure (sub_flux), pointer, intent(in) :: flux_ptr
    integer, intent(in) :: order
    real(dp), intent(in) :: cfl,tf
    real(dp), intent(inout) :: t
    procedure (quadrature_t), pointer, intent(in) :: quad_t
    procedure (quadrature_c_alpha), pointer, intent(in) :: quad_c_alpha
    procedure (quadrature_reconstruction), pointer, intent(in) :: quad_reconstruct
    real(dp) :: t1,t2
    type(solStruct) :: sol1,sol2,sol3

    allocate(sol1%val(mesh%nc,sol%nvar),sol2%val(mesh%nc,sol%nvar),sol3%val(mesh%nc,sol%nvar))
    sol1=sol
    sol2=sol
    sol3=sol

    t1=t
    call advance(mesh,sol,sol1,f_ptr,flux_ptr,order,cfl,t1,tf,quad_t,quad_c_alpha,quad_reconstruct)
    t2=t
    call advance(mesh,sol1,sol2,f_ptr,flux_ptr,order,cfl,t2,tf,quad_t,quad_c_alpha,quad_reconstruct)
    sol2%val=0.75_dp*sol%val+0.25_dp*sol2%val    
    call advance(mesh,sol2,sol3,f_ptr,flux_ptr,order,cfl,t,tf,quad_t,quad_c_alpha,quad_reconstruct)
    sol%val=1.0_dp/3.0_dp*sol%val+2.0_dp/3.0_dp*sol3%val

    deallocate(sol1%val,sol2%val,sol3%val)

    return
  end subroutine SSPRK3

end module time
