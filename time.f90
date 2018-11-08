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
    real(dp) :: S
    integer :: k,i,normal,neigh

    dt=tf-t
    do k=1,mesh%nc
       do i=1,size(mesh%cell(k)%edge)
          normal=mesh%cell(k)%edge(i)%normal
          neigh=mesh%cell(k)%edge(i)%neigh
          if (neigh>0) then
             select case (normal)
             case (1)
                call speed(sol%val(neigh,:),sol%val(k,:),f_equa,1,S)
             case (2)
                call speed(sol%val(neigh,:),sol%val(k,:),f_equa,2,S)
             case (3)
                call speed(sol%val(k,:),sol%val(neigh,:),f_equa,1,S)
             case (4)
                call speed(sol%val(k,:),sol%val(neigh,:),f_equa,2,S)
             end select
             dt=min(dt,cfl*min(mesh%cell(k)%dx,mesh%cell(k)%dy)/S)
          endif
       enddo
    enddo

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
    integer :: k,i,neigh,normal
    real(dp), dimension(:,:,:), allocatable :: F   !F(k,var,dir)
    real(dp), dimension(:), allocatable :: Ftemp,u1,u2
    type(cellStruct) :: cell
    type(edgeStruct) :: edge
    real(dp) :: xi,yi,dx,dy
    procedure (sub_reconstruction), pointer :: func

    allocate(F(mesh%nc,sol%nvar,4),Ftemp(sol%nvar),u1(sol%nvar),u2(sol%nvar))

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
             call boundary(flux,f_equa,quad_c_alpha,quad_reconstruct,neigh, &
                  u1,mesh,sol,edge%boundType,edge%bound,normal,order,Ftemp(:))
          else
             if ((normal==1).or.(normal==2)) then
                call quad_reconstruct(func,mesh,sol,order,quad_c_alpha,normal+2,neigh,u2)
                call flux(u2,u1,f_equa,normal,Ftemp(:))
             else
                call quad_reconstruct(func,mesh,sol,order,quad_c_alpha,normal-2,neigh,u2)
                call flux(u1,u2,f_equa,normal,Ftemp(:))
             endif
          endif
          F(k,:,normal)=F(k,:,normal)+Ftemp(:)
       enddo
    enddo
    
    sol2%val(1:mesh%nc,:)=sol%val(1:mesh%nc,:)-dt/dx*(F(1:mesh%nc,:,3)-F(1:mesh%nc,:,1))-dt/dy*(F(1:mesh%nc,:,4)-F(1:mesh%nc,:,2))

    do k=1,mesh%nc
       deallocate(mesh%cell(k)%polCoef)
    enddo
    
    deallocate(F,Ftemp,u1,u2)
    
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
