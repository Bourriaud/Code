module time

  use constant
  use types
  use FV
  use reconstruction
  
  implicit none

contains

  subroutine advance(mesh,sol,sol2,f_ptr,flux_ptr,order,cfl,t)
    type(meshStruct), intent(in) :: mesh
    type(solStruct), intent(in) :: sol
    type(solStruct), intent(inout) :: sol2
    procedure (sub_f), pointer, intent(in) :: f_ptr
    procedure (sub_flux), pointer, intent(in) :: flux_ptr
    integer, intent(in) :: order
    real(dp), intent(in) :: cfl
    real(dp), intent(inout) :: t
    integer :: k,i,neigh,normal
    real(dp), dimension(:,:,:), allocatable :: F   !F(k,var,dir)
    real(dp), dimension(:), allocatable :: Ftemp,ul,ur
    type(cellStruct) :: cell
    type(edgeStruct) :: edge
    real(dp) :: dx,dy,dt,S

    allocate(F(mesh%nc,sol%nvar,4),Ftemp(sol%nvar),ul(sol%nvar),ur(sol%nvar))

    dt=10.0_dp**20
    F=0.0_dp
    do k=1,mesh%nc
       cell=mesh%cell(k)
       dx=cell%dx
       dy=cell%dy
       do i=1,size(cell%edge)
          edge=cell%edge(i)
          neigh=edge%neigh
          normal=edge%normal
          if (neigh<0) then
             call reconstruct_boundary(mesh,sol,k,order,normal,ul)
             call boundary(flux_ptr,f_ptr,neigh,ul,mesh,sol,edge%boundType,edge%bound,normal,order,Ftemp(:),S)
          else
             call reconstruct(mesh,sol,k,neigh,order,normal,ul,ur)
             call flux_ptr(ul,ur,f_ptr,normal,Ftemp(:),S)
          endif
          F(k,:,normal)=F(k,:,normal)+Ftemp(:)
          dt=min(dt,cfl*min(dx,dy)/S)
       enddo
    enddo
    
    sol2%val(1:mesh%nc,:)=sol%val(1:mesh%nc,:)-dt/dx*(F(1:mesh%nc,:,3)-F(1:mesh%nc,:,1))-dt/dy*(F(1:mesh%nc,:,4)-F(1:mesh%nc,:,2))
    t=t+dt

    deallocate(F,Ftemp,ul,ur)
    
    return
  end subroutine advance

  subroutine euler_exp(mesh,sol,f_ptr,flux_ptr,order,cfl,t)
    type(meshStruct), intent(in) :: mesh
    type(solStruct), intent(inout) :: sol
    procedure (sub_f), pointer, intent(in) :: f_ptr
    procedure (sub_flux), pointer, intent(in) :: flux_ptr
    integer, intent(in) :: order
    real(dp), intent(in) :: cfl
    real(dp), intent(inout) :: t
    type(solStruct) :: sol1

    allocate(sol1%val(mesh%nc,sol%nvar))
    sol1=sol
    
    call advance(mesh,sol,sol1,f_ptr,flux_ptr,order,cfl,t)
    sol%val=sol1%val

    deallocate(sol1%val)

    return
  end subroutine euler_exp

  subroutine SSPRK2(mesh,sol,f_ptr,flux_ptr,order,cfl,t)
    type(meshStruct), intent(in) :: mesh
    type(solStruct), intent(inout) :: sol
    procedure (sub_f), pointer, intent(in) :: f_ptr
    procedure (sub_flux), pointer, intent(in) :: flux_ptr
    integer, intent(in) :: order
    real(dp), intent(in) :: cfl
    real(dp), intent(inout) :: t
    real(dp) :: t1
    type(solStruct) :: sol1,sol2

    allocate(sol1%val(mesh%nc,sol%nvar),sol2%val(mesh%nc,sol%nvar))
    sol1=sol
    sol2=sol

    call advance(mesh,sol,sol1,f_ptr,flux_ptr,order,cfl,t)
    t1=t
    call advance(mesh,sol1,sol2,f_ptr,flux_ptr,order,cfl,t1)
    sol%val=0.5_dp*(sol%val+sol2%val)

    deallocate(sol1%val,sol2%val)

    return
  end subroutine SSPRK2

  subroutine SSPRK3(mesh,sol,f_ptr,flux_ptr,order,cfl,t)
    type(meshStruct), intent(in) :: mesh
    type(solStruct), intent(inout) :: sol
    procedure (sub_f), pointer, intent(in) :: f_ptr
    procedure (sub_flux), pointer, intent(in) :: flux_ptr
    integer, intent(in) :: order
    real(dp), intent(in) :: cfl
    real(dp), intent(inout) :: t
    real(dp) :: t1,t2
    type(solStruct) :: sol1,sol2,sol3

    allocate(sol1%val(mesh%nc,sol%nvar),sol2%val(mesh%nc,sol%nvar),sol3%val(mesh%nc,sol%nvar))
    sol1=sol
    sol2=sol
    sol3=sol
    
    call advance(mesh,sol,sol1,f_ptr,flux_ptr,order,cfl,t)
    t1=t
    call advance(mesh,sol1,sol2,f_ptr,flux_ptr,order,cfl,t1)
    t2=t1
    sol2%val=0.75_dp*sol%val+0.25_dp*sol2%val    
    call advance(mesh,sol2,sol3,f_ptr,flux_ptr,order,cfl,t2)   
    sol%val=1.0_dp/3.0_dp*sol%val+2.0_dp/3.0_dp*sol3%val

    deallocate(sol1%val,sol2%val,sol3%val)

    return
  end subroutine SSPRK3

end module time
