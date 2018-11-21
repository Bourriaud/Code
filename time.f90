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

  subroutine advance(mesh,sol,sol2,f_equa,flux,order,dt,n,L_str_criteria,L_var_criteria,gauss_weight)
    type(meshStruct), intent(inout) :: mesh
    type(solStruct), intent(in) :: sol
    type(solStruct), intent(inout) :: sol2
    procedure (sub_f), pointer, intent(in) :: f_equa
    procedure (sub_flux), pointer, intent(in) :: flux
    integer, intent(in) :: order,n
    real(dp), intent(in) :: dt
    character(len=20), dimension(:), intent(in) :: L_str_criteria
    integer, dimension(:), intent(in) :: L_var_criteria
    real(dp), dimension(:), intent(in) :: gauss_weight
    integer :: k,i,j,p,cell1,cell2,dir,count,deg
    real(dp), dimension(:), allocatable :: u1,u2
    type(edgeStruct) :: edge
    integer, dimension(:), allocatable :: NOT_ACCEPTED_CELL,NOT_ACCEPTED_EDGE
    type(solStruct) :: soltemp
    integer, dimension(2) :: L_deg

    allocate(u1(sol%nvar),u2(sol%nvar))
    allocate(NOT_ACCEPTED_CELL(mesh%nc),NOT_ACCEPTED_EDGE(mesh%ne),soltemp%val(mesh%nc,sol%nvar))
    
    sol2=sol
    soltemp=sol
    
    do k=1,mesh%nc
       mesh%cell(k)%deg=order-1
       NOT_ACCEPTED_CELL(k)=k
    enddo
    do k=1,mesh%ne
       NOT_ACCEPTED_EDGE(k)=k
    enddo

    L_deg(1)=0
    L_deg(2)=0
    
    count=0
    do while (size(NOT_ACCEPTED_CELL)>0)
       count=count+1
       do i=1,size(NOT_ACCEPTED_EDGE)
          j=NOT_ACCEPTED_EDGE(i)
          edge=mesh%edge(j)
          cell1=edge%cell1
          cell2=edge%cell2
          dir=edge%dir
          deg=min(mesh%cell(abs(cell1))%deg,mesh%cell(abs(cell2))%deg)
          mesh%edge(j)%deg=deg
          call reconstruct(mesh,sol,abs(cell1),deg+1,gauss_weight)
          call reconstruct(mesh,sol,abs(cell2),deg+1,gauss_weight)
          do p=1,order
             if (.not.edge%flux_acc(p)) then
                if (cell1<0) then
                   call evaluate(mesh,sol,cell2,deg+1,gauss_weight,edge%X_gauss(p),edge%Y_gauss(p),u2)
                   call boundary(flux,f_equa,gauss_weight,cell1, &
                        u2,mesh,sol,j,p,edge%boundType,edge%bound,dir,deg+1,mesh%edge(j)%flux(p,:))
                   soltemp%val(cell2,:)=soltemp%val(cell2,:)+gauss_weight(p)*mesh%edge(j)%flux(p,:)*dt/(mesh%edge(j)%length*2.0_dp)
                elseif (cell2<0) then
                   call evaluate(mesh,sol,cell1,deg+1,gauss_weight,edge%X_gauss(p),edge%Y_gauss(p),u1)
                   call boundary(flux,f_equa,gauss_weight,cell2, &
                        u1,mesh,sol,j,p,edge%boundType,edge%bound,dir+2,deg+1,mesh%edge(j)%flux(p,:))
                   soltemp%val(cell1,:)=soltemp%val(cell1,:)-gauss_weight(p)*mesh%edge(j)%flux(p,:)*dt/(mesh%edge(j)%length*2.0_dp)
                else
                   call evaluate(mesh,sol,cell1,deg+1,gauss_weight,edge%X_gauss(p),edge%Y_gauss(p),u1)
                   call evaluate(mesh,sol,cell2,deg+1,gauss_weight,edge%X_gauss(p),edge%Y_gauss(p),u2)
                   call flux(u1,u2,f_equa,dir,mesh%edge(j)%flux(p,:))
                   soltemp%val(cell1,:)=soltemp%val(cell1,:)-gauss_weight(p)*mesh%edge(j)%flux(p,:)*dt/(mesh%edge(j)%length*2.0_dp)
                   soltemp%val(cell2,:)=soltemp%val(cell2,:)+gauss_weight(p)*mesh%edge(j)%flux(p,:)*dt/(mesh%edge(j)%length*2.0_dp)
                endif
             endif
          enddo
       enddo

       do k=1,size(NOT_ACCEPTED_EDGE)
          cell1=mesh%edge(NOT_ACCEPTED_EDGE(k))%cell1
          cell2=mesh%edge(NOT_ACCEPTED_EDGE(k))%cell2
          if (cell2>0) then
             if(allocated(mesh%cell(cell2)%polCoef))deallocate(mesh%cell(cell2)%polCoef)
          endif
          if (cell1>0) then
             if(allocated(mesh%cell(cell1)%polCoef))deallocate(mesh%cell(cell1)%polCoef)
          endif
       enddo
       
       deg=L_deg(min(count,size(L_deg)))
       call decrement(mesh,sol,soltemp,deg,dt,L_str_criteria,L_var_criteria,gauss_weight,NOT_ACCEPTED_CELL,NOT_ACCEPTED_EDGE)

       !if(count>10)call exit()
       !print*,size(NOT_ACCEPTED_CELL)
       !print*,"-----------------------------------"
       !call write_accept(mesh,sol,NOT_ACCEPTED_CELL,n,count)

    enddo

    sol2%val=soltemp%val
    deallocate(u1,u2,NOT_ACCEPTED_CELL,NOT_ACCEPTED_EDGE,soltemp%val)

    return
  end subroutine advance

  subroutine euler_exp(mesh,sol,f_equa,flux,speed,order,cfl,t,n,tf,L_str_criteria,L_var_criteria,gauss_weight)
    type(meshStruct), intent(inout) :: mesh
    type(solStruct), intent(inout) :: sol
    procedure (sub_f), pointer, intent(in) :: f_equa
    procedure (sub_flux), pointer, intent(in) :: flux
    procedure (sub_speed), pointer, intent(in) :: speed
    integer, intent(in) :: order,n
    real(dp), intent(in) :: cfl,tf
    real(dp), intent(inout) :: t
    character(len=20), dimension(:), intent(in) :: L_str_criteria
    integer, dimension(:), intent(in) :: L_var_criteria
    real(dp), dimension(:), intent(in) :: gauss_weight
    type(solStruct) :: sol1
    real(dp) :: dt

    allocate(sol1%val(mesh%nc,sol%nvar))
    sol1=sol

    call compute_timestep(mesh,sol,f_equa,speed,cfl,tf,t,dt)

    call advance(mesh,sol,sol1,f_equa,flux,order,dt,n,L_str_criteria,L_var_criteria,gauss_weight)
    sol%val=sol1%val

    t=t+dt

    deallocate(sol1%val)

    return
  end subroutine euler_exp

  subroutine SSPRK2(mesh,sol,f_equa,flux,speed,order,cfl,t,n,tf,L_str_criteria,L_var_criteria,gauss_weight)
    type(meshStruct), intent(inout) :: mesh
    type(solStruct), intent(inout) :: sol
    procedure (sub_f), pointer, intent(in) :: f_equa
    procedure (sub_flux), pointer, intent(in) :: flux
    procedure (sub_speed), pointer, intent(in) :: speed
    integer, intent(in) :: order,n
    real(dp), intent(in) :: cfl,tf
    real(dp), intent(inout) :: t
    character(len=20), dimension(:), intent(in) :: L_str_criteria
    integer, dimension(:), intent(in) :: L_var_criteria
    real(dp), dimension(:), intent(in) :: gauss_weight
    type(solStruct) :: sol1,sol2
    real(dp) :: dt

    allocate(sol1%val(mesh%nc,sol%nvar),sol2%val(mesh%nc,sol%nvar))
    sol1=sol
    sol2=sol

    call compute_timestep(mesh,sol,f_equa,speed,cfl,tf,t,dt)

    call advance(mesh,sol,sol1,f_equa,flux,order,dt,n,L_str_criteria,L_var_criteria,gauss_weight)
    call advance(mesh,sol1,sol2,f_equa,flux,order,dt,n,L_str_criteria,L_var_criteria,gauss_weight)
    sol%val=0.5_dp*(sol%val+sol2%val)

    t=t+dt

    deallocate(sol1%val,sol2%val)

    return
  end subroutine SSPRK2

  subroutine SSPRK3(mesh,sol,f_equa,flux,speed,order,cfl,t,n,tf,L_str_criteria,L_var_criteria,gauss_weight)
    type(meshStruct), intent(inout) :: mesh
    type(solStruct), intent(inout) :: sol
    procedure (sub_f), pointer, intent(in) :: f_equa
    procedure (sub_flux), pointer, intent(in) :: flux
    procedure (sub_speed), pointer, intent(in) :: speed
    integer, intent(in) :: order,n
    real(dp), intent(in) :: cfl,tf
    real(dp), intent(inout) :: t
    character(len=20), dimension(:), intent(in) :: L_str_criteria
    integer, dimension(:), intent(in) :: L_var_criteria
    real(dp), dimension(:), intent(in) :: gauss_weight
    type(solStruct) :: sol1,sol2,sol3
    real(dp) :: dt

    allocate(sol1%val(mesh%nc,sol%nvar),sol2%val(mesh%nc,sol%nvar),sol3%val(mesh%nc,sol%nvar))
    sol1=sol
    sol2=sol
    sol3=sol

    call compute_timestep(mesh,sol,f_equa,speed,cfl,tf,t,dt)

    call advance(mesh,sol,sol1,f_equa,flux,order,dt,n,L_str_criteria,L_var_criteria,gauss_weight)
    call advance(mesh,sol1,sol2,f_equa,flux,order,dt,n,L_str_criteria,L_var_criteria,gauss_weight)
    sol2%val=0.75_dp*sol%val+0.25_dp*sol2%val    
    call advance(mesh,sol2,sol3,f_equa,flux,order,dt,n,L_str_criteria,L_var_criteria,gauss_weight)
    sol%val=1.0_dp/3.0_dp*sol%val+2.0_dp/3.0_dp*sol3%val

    t=t+dt

    deallocate(sol1%val,sol2%val,sol3%val)

    return
  end subroutine SSPRK3

end module time
