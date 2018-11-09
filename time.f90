module time

  use constant
  use types
  use FV
  use reconstruction
  
  implicit none

contains

  subroutine compute_timestep(mesh,sol,f_equa,speed,cfl,tf,t,dt)
    type(meshStruct), intent(in) :: mesh
    type(solStruct), intent(in) :: sol
    procedure (sub_f), pointer, intent(in) :: f_equa
    procedure (sub_speed), pointer, intent(in) :: speed
    real(dp), intent(in) :: cfl,tf,t
    real(dp), intent(out) :: dt
    real(dp) :: dx,dy
    real(dp), dimension(2) :: S
    integer :: i,cell1,cell2

    dt=tf-t

    do i=1,mesh%ne
       cell1=mesh%edge(i)%cell1
       cell2=mesh%edge(i)%cell2
       if ((cell1>0).and.(cell2>0)) then
          call speed(sol%val(cell1,:),sol%val(cell2,:),f_equa,S)
          dx=min(mesh%cell(cell1)%dx,mesh%cell(cell2)%dx)
          dy=min(mesh%cell(cell1)%dy,mesh%cell(cell2)%dy)
          dt=min(dt,cfl/(S(1)/dx+S(2)/dy))
       end if
    end do       
    
    return
  end subroutine compute_timestep

  subroutine advance(mesh,sol,sol2,f_equa,flux,order,dt,quad_c_alpha,quad_reconstruct)
    type(meshStruct), intent(inout) :: mesh
    type(solStruct), intent(in) :: sol
    type(solStruct), intent(inout) :: sol2
    procedure (sub_f), pointer, intent(in) :: f_equa
    procedure (sub_flux), pointer, intent(in) :: flux
    integer, intent(in) :: order
    real(dp), intent(in) :: dt
    procedure (quadrature_c_alpha), pointer, intent(in) :: quad_c_alpha
    procedure (quadrature_reconstruction), pointer, intent(in) :: quad_reconstruct
    integer :: k,i,cell1,cell2,dir
    real(dp), dimension(:), allocatable :: Ftemp,u1,u2
    type(edgeStruct) :: edge
    procedure (sub_reconstruction), pointer :: func

    allocate(Ftemp(sol%nvar),u1(sol%nvar),u2(sol%nvar))
    
    func => evaluate
    sol2%val=sol%val
    
    do k=1,mesh%nc
       call reconstruct(mesh,sol,k,order,quad_c_alpha)
    enddo

    do i=1,mesh%ne
       edge=mesh%edge(i)
       cell1=edge%cell1
       cell2=edge%cell2
       dir=edge%dir
       if (cell1<0) then
          call quad_reconstruct(func,mesh,sol,order,quad_c_alpha,dir,cell2,u2)
          call boundary(flux,f_equa,quad_c_alpha,quad_reconstruct,cell1, &
               u2,mesh,sol,edge%boundType,edge%bound,dir,order,Ftemp(:))
          sol2%val(cell2,:)=sol2%val(cell2,:)+Ftemp(:)*dt/edge%length
       elseif (cell2<0) then
          call quad_reconstruct(func,mesh,sol,order,quad_c_alpha,dir+2,cell1,u1)
          call boundary(flux,f_equa,quad_c_alpha,quad_reconstruct,cell2, &
               u1,mesh,sol,edge%boundType,edge%bound,dir+2,order,Ftemp(:))
          sol2%val(cell1,:)=sol2%val(cell1,:)-Ftemp(:)*dt/edge%length
       else
          call quad_reconstruct(func,mesh,sol,order,quad_c_alpha,dir+2,cell1,u1)
          call quad_reconstruct(func,mesh,sol,order,quad_c_alpha,dir,cell2,u2)
          call flux(u1,u2,f_equa,dir,Ftemp(:))
          sol2%val(cell1,:)=sol2%val(cell1,:)-Ftemp(:)*dt/edge%length
          sol2%val(cell2,:)=sol2%val(cell2,:)+Ftemp(:)*dt/edge%length
       endif
    enddo

    do k=1,mesh%nc
       deallocate(mesh%cell(k)%polCoef)
    enddo
    
    deallocate(Ftemp,u1,u2)
    
    return
  end subroutine advance

  subroutine euler_exp(mesh,sol,f_equa,flux,speed,order,cfl,t,tf,quad_c_alpha,quad_reconstruct)
    type(meshStruct), intent(inout) :: mesh
    type(solStruct), intent(inout) :: sol
    procedure (sub_f), pointer, intent(in) :: f_equa
    procedure (sub_flux), pointer, intent(in) :: flux
    procedure (sub_speed), pointer, intent(in) :: speed
    integer, intent(in) :: order
    real(dp), intent(in) :: cfl,tf
    real(dp), intent(inout) :: t
    procedure (quadrature_c_alpha), pointer, intent(in) :: quad_c_alpha
    procedure (quadrature_reconstruction), pointer, intent(in) :: quad_reconstruct
    type(solStruct) :: sol1
    real(dp) :: dt

    allocate(sol1%val(mesh%nc,sol%nvar))
    sol1=sol

    call compute_timestep(mesh,sol,f_equa,speed,cfl,tf,t,dt)

    call advance(mesh,sol,sol1,f_equa,flux,order,dt,quad_c_alpha,quad_reconstruct)
    sol%val=sol1%val

    t=t+dt

    deallocate(sol1%val)

    return
  end subroutine euler_exp

  subroutine SSPRK2(mesh,sol,f_equa,flux,speed,order,cfl,t,tf,quad_c_alpha,quad_reconstruct)
    type(meshStruct), intent(inout) :: mesh
    type(solStruct), intent(inout) :: sol
    procedure (sub_f), pointer, intent(in) :: f_equa
    procedure (sub_flux), pointer, intent(in) :: flux
    procedure (sub_speed), pointer, intent(in) :: speed
    integer, intent(in) :: order
    real(dp), intent(in) :: cfl,tf
    real(dp), intent(inout) :: t
    procedure (quadrature_c_alpha), pointer, intent(in) :: quad_c_alpha
    procedure (quadrature_reconstruction), pointer, intent(in) :: quad_reconstruct
    type(solStruct) :: sol1,sol2
    real(dp) :: dt

    allocate(sol1%val(mesh%nc,sol%nvar),sol2%val(mesh%nc,sol%nvar))
    sol1=sol
    sol2=sol

    call compute_timestep(mesh,sol,f_equa,speed,cfl,tf,t,dt)

    call advance(mesh,sol,sol1,f_equa,flux,order,dt,quad_c_alpha,quad_reconstruct)
    call advance(mesh,sol1,sol2,f_equa,flux,order,dt,quad_c_alpha,quad_reconstruct)
    sol%val=0.5_dp*(sol%val+sol2%val)

    t=t+dt

    deallocate(sol1%val,sol2%val)

    return
  end subroutine SSPRK2

  subroutine SSPRK3(mesh,sol,f_equa,flux,speed,order,cfl,t,tf,quad_c_alpha,quad_reconstruct)
    type(meshStruct), intent(inout) :: mesh
    type(solStruct), intent(inout) :: sol
    procedure (sub_f), pointer, intent(in) :: f_equa
    procedure (sub_flux), pointer, intent(in) :: flux
    procedure (sub_speed), pointer, intent(in) :: speed
    integer, intent(in) :: order
    real(dp), intent(in) :: cfl,tf
    real(dp), intent(inout) :: t
    procedure (quadrature_c_alpha), pointer, intent(in) :: quad_c_alpha
    procedure (quadrature_reconstruction), pointer, intent(in) :: quad_reconstruct
    type(solStruct) :: sol1,sol2,sol3
    real(dp) :: dt

    allocate(sol1%val(mesh%nc,sol%nvar),sol2%val(mesh%nc,sol%nvar),sol3%val(mesh%nc,sol%nvar))
    sol1=sol
    sol2=sol
    sol3=sol

    call compute_timestep(mesh,sol,f_equa,speed,cfl,tf,t,dt)

    call advance(mesh,sol,sol1,f_equa,flux,order,dt,quad_c_alpha,quad_reconstruct)
    call advance(mesh,sol1,sol2,f_equa,flux,order,dt,quad_c_alpha,quad_reconstruct)
    sol2%val=0.75_dp*sol%val+0.25_dp*sol2%val    
    call advance(mesh,sol2,sol3,f_equa,flux,order,dt,quad_c_alpha,quad_reconstruct)
    sol%val=1.0_dp/3.0_dp*sol%val+2.0_dp/3.0_dp*sol3%val

    t=t+dt

    deallocate(sol1%val,sol2%val,sol3%val)

    return
  end subroutine SSPRK3

end module time
