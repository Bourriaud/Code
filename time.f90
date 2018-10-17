module time

  use types
  use FV
  
  implicit none

contains

  subroutine advance(mesh,sol,sol2,f_ptr,flux_ptr,cfl,t)
    type(meshStruct), intent(in) :: mesh
    type(solStruct), intent(in) :: sol
    type(solStruct), intent(inout) :: sol2
    procedure (sub_f), pointer, intent(in) :: f_ptr
    procedure (sub_flux), pointer, intent(in) :: flux_ptr
    real, intent(in) :: cfl
    real, intent(inout) :: t
    integer :: k,i,neigh,normal,k1,k2
    real, dimension(:,:), allocatable :: F
    real, dimension(:), allocatable :: Ftemp
    type(cellStruct) :: cell
    type(edgeStruct) :: edge
    real :: dx,dy,dt,S
    
    allocate(F(sol%nvar,4),Ftemp(sol%nvar))

    dt=1.
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
             call boundary(flux_ptr,f_ptr,sol%val(k,:),edge%boundType,edge%bound,normal,Ftemp(:),S)
          else
             k1=min(k,neigh)
             k2=max(k,neigh)
             call flux_ptr(sol%val(k1,:),sol%val(k2,:),f_ptr,normal,Ftemp(:),S)
          endif
          F(:,normal)=F(:,normal)+Ftemp(:)
          dt=min(dt,cfl*min(dx,dy)/S)
       enddo
       sol2%val(k,:)=sol%val(k,:)-dt/dx*(F(:,3)-F(:,1))-dt/dy*(F(:,4)-F(:,2))
    enddo
    t=t+dt

    deallocate(F,Ftemp)
    
    return
  end subroutine advance

  subroutine euler_exp(mesh,sol,f_ptr,flux_ptr,cfl,t)
    type(meshStruct), intent(in) :: mesh
    type(solStruct), intent(inout) :: sol
    procedure (sub_f), pointer, intent(in) :: f_ptr
    procedure (sub_flux), pointer, intent(in) :: flux_ptr
    real, intent(in) :: cfl
    real, intent(inout) :: t
    type(solStruct) :: sol1

    allocate(sol1%val(mesh%nc,sol%nvar))
    
    call advance(mesh,sol,sol1,f_ptr,flux_ptr,cfl,t)
    sol%val=sol1%val

    deallocate(sol1%val)

    return
  end subroutine euler_exp

  subroutine SSPRK2(mesh,sol,f_ptr,flux_ptr,cfl,t)
    type(meshStruct), intent(in) :: mesh
    type(solStruct), intent(inout) :: sol
    procedure (sub_f), pointer, intent(in) :: f_ptr
    procedure (sub_flux), pointer, intent(in) :: flux_ptr
    real, intent(in) :: cfl
    real, intent(inout) :: t
    real :: t1
    type(solStruct) :: sol1,sol2

    allocate(sol1%val(mesh%nc,sol%nvar),sol2%val(mesh%nc,sol%nvar))

    call advance(mesh,sol,sol1,f_ptr,flux_ptr,cfl,t)
    t1=t
    call advance(mesh,sol1,sol2,f_ptr,flux_ptr,cfl,t1)
    sol%val=(sol%val+sol2%val)/2.

    deallocate(sol1%val,sol2%val)

    return
  end subroutine SSPRK2

  subroutine SSPRK3(mesh,sol,f_ptr,flux_ptr,cfl,t)
    type(meshStruct), intent(in) :: mesh
    type(solStruct), intent(inout) :: sol
    procedure (sub_f), pointer, intent(in) :: f_ptr
    procedure (sub_flux), pointer, intent(in) :: flux_ptr
    real, intent(in) :: cfl
    real, intent(inout) :: t
    real :: t1,t2
    type(solStruct) :: sol1,sol2

    allocate(sol1%val(mesh%nc,sol%nvar),sol2%val(mesh%nc,sol%nvar))
    
    call advance(mesh,sol,sol1,f_ptr,flux_ptr,cfl,t)
    t1=t
    call advance(mesh,sol1,sol2,f_ptr,flux_ptr,cfl,t1)
    t2=t1
    sol1%val=3./4.*sol%val+1./4.*sol1%val    
    call advance(mesh,sol1,sol2,f_ptr,flux_ptr,cfl,t2)   
    sol%val=1./3.*sol%val+2./3.*sol2%val

    deallocate(sol1%val,sol2%val)

    return
  end subroutine SSPRK3

end module time
