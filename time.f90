module time

  use constant
  use types
  use FV
  use reconstruction
  use limit
  
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

  subroutine advance(mesh,sol,sol2,f_equa,flux,order,dt,quad_c_alpha,quad_reconstruct,L_str_criteria,L_var_criteria)
    type(meshStruct), intent(inout) :: mesh
    type(solStruct), intent(in) :: sol
    type(solStruct), intent(inout) :: sol2
    procedure (sub_f), pointer, intent(in) :: f_equa
    procedure (sub_flux), pointer, intent(in) :: flux
    integer, intent(in) :: order
    real(dp), intent(in) :: dt
    procedure (quadrature_c_alpha), pointer, intent(in) :: quad_c_alpha
    procedure (quadrature_reconstruction), pointer, intent(in) :: quad_reconstruct
    character(len=20), dimension(:), intent(in) :: L_str_criteria
    integer, dimension(:), intent(in) :: L_var_criteria
    integer :: k,i,j,cell1,cell2,dir,count,deg
    real(dp), dimension(:), allocatable :: u1,u2
    type(edgeStruct) :: edge
    procedure (sub_reconstruction), pointer :: func
    integer, dimension(:), allocatable :: NOT_ACCEPTED_CELL,NOT_ACCEPTED_EDGE
    type(solStruct) :: soltemp
    integer, dimension(2) :: L_deg

    allocate(u1(sol%nvar),u2(sol%nvar))
    allocate(NOT_ACCEPTED_CELL(mesh%nc),NOT_ACCEPTED_EDGE(mesh%ne),soltemp%val(mesh%nc,sol%nvar))
    
    func => evaluate
    sol2=sol
    soltemp=sol
    
    do k=1,mesh%nc
       mesh%cell(k)%deg=order-1
       NOT_ACCEPTED_CELL(k)=k
    enddo
    do k=1,mesh%ne
       NOT_ACCEPTED_EDGE(k)=k
    enddo

    L_deg(1)=1
    L_deg(2)=0
    
    count=0
    do while (size(NOT_ACCEPTED_CELL)>0)
       count=count+1
       do i=1,size(NOT_ACCEPTED_EDGE)
          j=NOT_ACCEPTED_EDGE(i)
          dir=mesh%edge(j)%dir
          cell1=mesh%edge(j)%cell1
          cell2=mesh%edge(j)%cell2
          mesh%edge(j)%deg=min(mesh%cell(abs(cell1))%deg,mesh%cell(abs(cell2))%deg)
          edge=mesh%edge(j)
          call reconstruct(mesh,sol,abs(cell1),dir+2,mesh%edge(j)%deg+1,quad_c_alpha)
          call reconstruct(mesh,sol,abs(cell2),dir,mesh%edge(j)%deg+1,quad_c_alpha)
          if (cell1<0) then
             call quad_reconstruct(func,mesh,sol,edge%deg+1,quad_c_alpha,dir,cell2,u2)
             call boundary(flux,f_equa,quad_c_alpha,quad_reconstruct,cell1, &
                  u2,mesh,sol,edge%boundType,edge%bound,dir,edge%deg+1,mesh%edge(j)%flux(:))
             soltemp%val(cell2,:)=soltemp%val(cell2,:)+mesh%edge(j)%flux(:)*dt/edge%length
          elseif (cell2<0) then
             call quad_reconstruct(func,mesh,sol,edge%deg+1,quad_c_alpha,dir+2,cell1,u1)
             call boundary(flux,f_equa,quad_c_alpha,quad_reconstruct,cell2, &
                  u1,mesh,sol,edge%boundType,edge%bound,dir+2,edge%deg+1,mesh%edge(j)%flux(:))
             soltemp%val(cell1,:)=soltemp%val(cell1,:)-mesh%edge(j)%flux(:)*dt/edge%length
          else
             call quad_reconstruct(func,mesh,sol,edge%deg+1,quad_c_alpha,dir+2,cell1,u1)
             call quad_reconstruct(func,mesh,sol,edge%deg+1,quad_c_alpha,dir,cell2,u2)
             call flux(u1,u2,f_equa,dir,mesh%edge(j)%flux(:))
             soltemp%val(cell1,:)=soltemp%val(cell1,:)-mesh%edge(j)%flux(:)*dt/edge%length
             soltemp%val(cell2,:)=soltemp%val(cell2,:)+mesh%edge(j)%flux(:)*dt/edge%length
          endif
       enddo

       do k=1,size(NOT_ACCEPTED_EDGE)
          cell1=mesh%edge(NOT_ACCEPTED_EDGE(k))%cell1
          cell2=mesh%edge(NOT_ACCEPTED_EDGE(k))%cell2
          if (cell2>0) then
             if(allocated(mesh%cell(cell2)%polCoefL))deallocate(mesh%cell(cell2)%polCoefL)
             if(allocated(mesh%cell(cell2)%polCoefB))deallocate(mesh%cell(cell2)%polCoefB)
          endif
          if (cell1>0) then
             if(allocated(mesh%cell(cell1)%polCoefR))deallocate(mesh%cell(cell1)%polCoefR)
             if(allocated(mesh%cell(cell1)%polCoefT))deallocate(mesh%cell(cell1)%polCoefT)
          endif
       enddo

       deg=L_deg(min(count,size(L_deg)))
       call decrement(mesh,sol,soltemp,deg,dt,L_str_criteria,L_var_criteria,NOT_ACCEPTED_CELL,NOT_ACCEPTED_EDGE)

       !if(count>19)call exit()      
       !print*,NOT_ACCEPTED_CELL
       !print*,"-----------------------------------"

    enddo

    sol2%val=soltemp%val
    
    deallocate(u1,u2,NOT_ACCEPTED_CELL,NOT_ACCEPTED_EDGE,soltemp%val)

    return
  end subroutine advance

  subroutine euler_exp(mesh,sol,f_equa,flux,speed,order,cfl,t,tf,quad_c_alpha,quad_reconstruct,L_str_criteria,L_var_criteria)
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
    character(len=20), dimension(:), intent(in) :: L_str_criteria
    integer, dimension(:), intent(in) :: L_var_criteria
    type(solStruct) :: sol1
    real(dp) :: dt

    allocate(sol1%val(mesh%nc,sol%nvar))
    sol1=sol

    call compute_timestep(mesh,sol,f_equa,speed,cfl,tf,t,dt)

    call advance(mesh,sol,sol1,f_equa,flux,order,dt,quad_c_alpha,quad_reconstruct,L_str_criteria,L_var_criteria)
    sol%val=sol1%val

    t=t+dt

    deallocate(sol1%val)

    return
  end subroutine euler_exp

  subroutine SSPRK2(mesh,sol,f_equa,flux,speed,order,cfl,t,tf,quad_c_alpha,quad_reconstruct,L_str_criteria,L_var_criteria)
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
    character(len=20), dimension(:), intent(in) :: L_str_criteria
    integer, dimension(:), intent(in) :: L_var_criteria
    type(solStruct) :: sol1,sol2
    real(dp) :: dt

    allocate(sol1%val(mesh%nc,sol%nvar),sol2%val(mesh%nc,sol%nvar))
    sol1=sol
    sol2=sol

    call compute_timestep(mesh,sol,f_equa,speed,cfl,tf,t,dt)

    call advance(mesh,sol,sol1,f_equa,flux,order,dt,quad_c_alpha,quad_reconstruct,L_str_criteria,L_var_criteria)
    call advance(mesh,sol1,sol2,f_equa,flux,order,dt,quad_c_alpha,quad_reconstruct,L_str_criteria,L_var_criteria)
    sol%val=0.5_dp*(sol%val+sol2%val)

    t=t+dt

    deallocate(sol1%val,sol2%val)

    return
  end subroutine SSPRK2

  subroutine SSPRK3(mesh,sol,f_equa,flux,speed,order,cfl,t,tf,quad_c_alpha,quad_reconstruct,L_str_criteria,L_var_criteria)
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
    character(len=20), dimension(:), intent(in) :: L_str_criteria
    integer, dimension(:), intent(in) :: L_var_criteria
    type(solStruct) :: sol1,sol2,sol3
    real(dp) :: dt

    allocate(sol1%val(mesh%nc,sol%nvar),sol2%val(mesh%nc,sol%nvar),sol3%val(mesh%nc,sol%nvar))
    sol1=sol
    sol2=sol
    sol3=sol

    call compute_timestep(mesh,sol,f_equa,speed,cfl,tf,t,dt)

    call advance(mesh,sol,sol1,f_equa,flux,order,dt,quad_c_alpha,quad_reconstruct,L_str_criteria,L_var_criteria)
    call advance(mesh,sol1,sol2,f_equa,flux,order,dt,quad_c_alpha,quad_reconstruct,L_str_criteria,L_var_criteria)
    sol2%val=0.75_dp*sol%val+0.25_dp*sol2%val    
    call advance(mesh,sol2,sol3,f_equa,flux,order,dt,quad_c_alpha,quad_reconstruct,L_str_criteria,L_var_criteria)
    sol%val=1.0_dp/3.0_dp*sol%val+2.0_dp/3.0_dp*sol3%val

    t=t+dt

    deallocate(sol1%val,sol2%val,sol3%val)

    return
  end subroutine SSPRK3

end module time
