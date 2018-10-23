module time

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
    real, intent(in) :: cfl
    real, intent(inout) :: t
    integer :: k,i,neigh,normal
    real, dimension(:,:), allocatable :: F
    real, dimension(:), allocatable :: Ftemp,ul,ur
    type(cellStruct) :: cell
    type(edgeStruct) :: edge
    real :: dx,dy,dt,S

    allocate(F(sol%nvar,4),Ftemp(sol%nvar),ul(sol%nvar),ur(sol%nvar))

    dt=10.**20
    do k=1,mesh%nc
       F=0.
       cell=mesh%cell(k)
       dx=cell%dx
       dy=cell%dy
       do i=1,size(cell%edge)
          edge=cell%edge(i)
          neigh=edge%neigh
          normal=edge%normal
          if (neigh==-1) then
             call reconstruct_boundary(mesh,sol,k,order,normal,ul)
             call boundary(flux_ptr,f_ptr,ul,edge%boundType,edge%bound,normal,Ftemp(:),S)
          else
             call reconstruct(mesh,sol,k,neigh,order,normal,ul,ur)
             call flux_ptr(ul,ur,f_ptr,normal,Ftemp(:),S)
          endif
          F(:,normal)=F(:,normal)+Ftemp(:)
          dt=min(dt,cfl*min(dx,dy)/S)
       enddo
       sol2%val(k,:)=sol%val(k,:)-dt/dx*(F(:,3)-F(:,1))-dt/dy*(F(:,4)-F(:,2))
    enddo
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
    real, intent(in) :: cfl
    real, intent(inout) :: t
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
    real, intent(in) :: cfl
    real, intent(inout) :: t
    real :: t1
    type(solStruct) :: sol1,sol2

    allocate(sol1%val(mesh%nc,sol%nvar),sol2%val(mesh%nc,sol%nvar))
    sol1=sol
    sol2=sol

    call advance(mesh,sol,sol1,f_ptr,flux_ptr,order,cfl,t)
    t1=t
    call advance(mesh,sol1,sol2,f_ptr,flux_ptr,order,cfl,t1)
    sol%val=(sol%val+sol2%val)/2.

    deallocate(sol1%val,sol2%val)

    return
  end subroutine SSPRK2

  subroutine SSPRK3(mesh,sol,f_ptr,flux_ptr,order,cfl,t)
    type(meshStruct), intent(in) :: mesh
    type(solStruct), intent(inout) :: sol
    procedure (sub_f), pointer, intent(in) :: f_ptr
    procedure (sub_flux), pointer, intent(in) :: flux_ptr
    integer, intent(in) :: order
    real, intent(in) :: cfl
    real, intent(inout) :: t
    real :: t1,t2
    type(solStruct) :: sol1,sol2,sol3

    allocate(sol1%val(mesh%nc,sol%nvar),sol2%val(mesh%nc,sol%nvar),sol3%val(mesh%nc,sol%nvar))
    sol1=sol
    sol2=sol
    sol3=sol
    
    call advance(mesh,sol,sol1,f_ptr,flux_ptr,order,cfl,t)
    t1=t
    call advance(mesh,sol1,sol2,f_ptr,flux_ptr,order,cfl,t1)
    t2=t1
    sol2%val=3./4.*sol%val+1./4.*sol2%val    
    call advance(mesh,sol2,sol3,f_ptr,flux_ptr,order,cfl,t2)   
    sol%val=1./3.*sol%val+2./3.*sol3%val

    deallocate(sol1%val,sol2%val,sol3%val)

    return
  end subroutine SSPRK3

end module time
