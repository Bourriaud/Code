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

  subroutine timestep(mesh,sol,str_equa,f_equa,flux,speed,time_scheme,order,nrk,cfl,t,n,tf, &
            cascade,L_str_criteria,L_var_criteria,L_eps,gauss_weight,period,verbosity,order_pc)
    type(meshStruct), intent(inout) :: mesh
    type(solStruct), intent(inout) :: sol
    character(len=20), intent(in) :: str_equa
    procedure (sub_f), pointer, intent(in) :: f_equa
    procedure (sub_flux), pointer, intent(in) :: flux
    procedure (sub_speed), pointer, intent(in) :: speed
    procedure (sub_time), pointer, intent(in) :: time_scheme
    integer, intent(in) :: order,nrk,n,verbosity
    real(dp), intent(in) :: cfl,tf
    real(dp), intent(inout) :: t
    character(len=20), dimension(:), intent(in) :: L_str_criteria
    integer, dimension(:), intent(in) :: cascade,L_var_criteria
    real(dp), dimension(:), intent(in) :: L_eps
    real(dp), dimension(:), intent(in) :: gauss_weight
    logical, dimension(2), intent(in) :: period
    integer, dimension(:), intent(inout) :: order_pc
    integer, dimension(:), allocatable :: NOT_ACCEPTED_CELL,NOT_ACCEPTED_EDGE,NOT_ACCEPTED_CELL_OLD,NOT_ACCEPTED_EDGE_OLD
    integer, dimension(:), allocatable :: NAC_cycle,NAC_reason
    type(solStruct) :: soltemp,soltemp2
    integer :: count,deg,k,j,cell1,cell2
    real(dp) :: dt
    if(.false.)print*,cell1,cell2

    call compute_timestep(mesh,sol,f_equa,speed,cfl,tf,t,dt)
    
    allocate(NOT_ACCEPTED_CELL(mesh%nc),NOT_ACCEPTED_CELL_OLD(mesh%nc),NOT_ACCEPTED_EDGE(mesh%ne),NOT_ACCEPTED_EDGE_OLD(mesh%ne))
    allocate(soltemp%val(mesh%nc,sol%nvar),NAC_cycle(mesh%nc),NAC_reason(mesh%nc))
    
    soltemp=sol
    soltemp2=sol
    NAC_cycle=0

    do k=1,mesh%nc
       mesh%cell(k)%accept=.false.
       mesh%cell(k)%updated=.false.
       mesh%cell(k)%deg=order-1
       NOT_ACCEPTED_CELL(k)=k
    enddo
    do k=1,mesh%ne
       NOT_ACCEPTED_EDGE(k)=k
       mesh%edge(k)%done=.false.
    enddo
    NOT_ACCEPTED_CELL_OLD=NOT_ACCEPTED_CELL
    NOT_ACCEPTED_EDGE_OLD=NOT_ACCEPTED_EDGE

    count=0
    
    do while (size(NOT_ACCEPTED_EDGE)>0)
       count=count+1

       do k=1,size(NOT_ACCEPTED_CELL)
          j=NOT_ACCEPTED_CELL(k)
       enddo

       do k=1,mesh%nc
          mesh%cell(k)%accept=.true.
          mesh%cell(k)%accept_temp=.true.
       enddo

       do k=1,size(NOT_ACCEPTED_EDGE)
          j=NOT_ACCEPTED_EDGE(k)
          cell1=mesh%edge(j)%cell1
          cell2=mesh%edge(j)%cell2
          if (mesh%edge(j)%done) then
             if (cell1>0) soltemp2%val(cell1,:)=soltemp2%val(cell1,:)-mesh%edge(j)%flux(nrk+1,1,:)
             if (cell2>0) soltemp2%val(cell2,:)=soltemp2%val(cell2,:)-mesh%edge(j)%flux(nrk+1,2,:)
          endif
       enddo
       
       call time_scheme(mesh,sol,soltemp2,soltemp,f_equa,flux,order,t,dt,n,gauss_weight,period, &
            order_pc,NOT_ACCEPTED_EDGE,str_equa,L_str_criteria,L_var_criteria,L_eps,NOT_ACCEPTED_CELL)

       if (size(L_str_criteria)>0.and.order>1) then
          deg=min(cascade(min(count,order-1)),order)-1
          call decrement(mesh,sol,soltemp,str_equa,deg,nrk,dt,L_str_criteria,L_var_criteria,L_eps, &
               gauss_weight,period,NOT_ACCEPTED_CELL,NOT_ACCEPTED_EDGE,NAC_reason,verbosity)
       else
          deg=order-1
          deallocate(NOT_ACCEPTED_CELL,NOT_ACCEPTED_EDGE)
          allocate(NOT_ACCEPTED_CELL(0),NOT_ACCEPTED_EDGE(0))
       endif

       do k=1,size(NOT_ACCEPTED_EDGE_OLD)
          j=NOT_ACCEPTED_EDGE_OLD(k)
          cell1=mesh%edge(j)%cell1
          cell2=mesh%edge(j)%cell2
          if (cell1>0) soltemp2%val(cell1,:)=soltemp2%val(cell1,:)+mesh%edge(j)%flux(nrk+1,1,:)
          if (cell2>0) soltemp2%val(cell2,:)=soltemp2%val(cell2,:)+mesh%edge(j)%flux(nrk+1,2,:)
          mesh%edge(j)%done=.true.
       enddo

       deallocate(NOT_ACCEPTED_CELL_OLD,NOT_ACCEPTED_EDGE_OLD)
       allocate(NOT_ACCEPTED_CELL_OLD(size(NOT_ACCEPTED_CELL)),NOT_ACCEPTED_EDGE_OLD(size(NOT_ACCEPTED_EDGE)))
       NOT_ACCEPTED_CELL_OLD=NOT_ACCEPTED_CELL
       NOT_ACCEPTED_EDGE_OLD=NOT_ACCEPTED_EDGE

       if (verbosity>1) then
          do k=1,size(NOT_ACCEPTED_CELL)
             NAC_cycle(NOT_ACCEPTED_CELL(k))=NAC_cycle(NOT_ACCEPTED_CELL(k))+1
          enddo
          if (size(NOT_ACCEPTED_EDGE)==0) then
             call write_accept(mesh,NAC_cycle,NAC_reason,n)
          endif
       endif
       
    enddo
    t=t+dt

    do k=1,mesh%nc
       if (allocated(mesh%cell(k)%polCoef)) deallocate(mesh%cell(k)%polCoef)
       allocate(mesh%cell(k)%polCoef(size(mesh%cell(k)%polTest(:,1)),size(mesh%cell(k)%polTest(1,:))))
       mesh%cell(k)%polCoef=mesh%cell(k)%polTest
    enddo

    sol=soltemp2

    deallocate(NOT_ACCEPTED_CELL,NOT_ACCEPTED_CELL_OLD,NOT_ACCEPTED_EDGE,NOT_ACCEPTED_EDGE_OLD)
    deallocate(soltemp%val,soltemp2%val,NAC_cycle,NAC_reason)

    return
  end subroutine timestep   

  subroutine advance(mesh,sol,sol2,f_equa,flux,order,dt,t,n,irk,gauss_weight,period, &
            order_pc,NOT_ACCEPTED_EDGE,str_equa,L_str_criteria,L_var_criteria,L_eps,NOT_ACCEPTED_CELL)
    type(meshStruct), intent(inout) :: mesh
    type(solStruct), intent(in) :: sol
    type(solStruct), intent(inout) :: sol2
    procedure (sub_f), pointer, intent(in) :: f_equa
    procedure (sub_flux), pointer, intent(in) :: flux
    integer, intent(in) :: order,n,irk
    real(dp), intent(in) :: dt,t
    real(dp), dimension(:), intent(in) :: gauss_weight
    logical, dimension(2), intent(in) :: period
    integer, dimension(:), intent(inout) :: order_pc
    integer, dimension(:), intent(in) :: NOT_ACCEPTED_EDGE
    character(len=20), intent(in) :: str_equa
    integer, dimension(:), intent(in) :: L_var_criteria
    character(len=20), dimension(:), intent(in) :: L_str_criteria
    real(dp), dimension(:), intent(in) :: L_eps
    integer, dimension(:), intent(inout) :: NOT_ACCEPTED_CELL
    integer :: i,j,p,cell1,cell2,dir,deg,sub1,sub2,ac1,ac2
    real(dp), dimension(:), allocatable :: u1,u2
    type(edgeStruct) :: edge
    type(solStruct) :: soltemp
    real(dp) :: lengthN1,lengthN2
    real(dp), dimension(:,:), allocatable :: flux_temp,flux_temp1,flux_temp2
    if(.false.)print*,irk,order_pc,t,n

    allocate(u1(sol%nvar),u2(sol%nvar),flux_temp(order,sol%nvar),flux_temp1(order,sol%nvar),flux_temp2(order,sol%nvar))
    
    sol2=sol
    soltemp=sol

    do i=1,size(NOT_ACCEPTED_EDGE)
       j=NOT_ACCEPTED_EDGE(i)
       edge=mesh%edge(j)
       cell1=edge%cell1
       cell2=edge%cell2
       ac1=abs(cell1)
       ac2=abs(cell2)
       dir=edge%dir
       lengthN1=mesh%cell(ac1)%dx
       lengthN2=mesh%cell(ac2)%dx
       sub1=edge%sub(1)
       sub2=edge%sub(2)

       deg=min(mesh%cell(ac1)%deg,mesh%cell(ac2)%deg)

       mesh%edge(j)%flux(irk,:,:)=0.0_dp
       mesh%edge(j)%deg=deg
       call reconstruct(mesh,sol,ac1,deg+1,gauss_weight,period,mesh%cell(ac1)%polTest)
       call reconstruct(mesh,sol,ac2,deg+1,gauss_weight,period,mesh%cell(ac2)%polTest)
          
       do p=1,order
          if (.not.edge%flux_acc(p)) then
             if (cell1<0) then
                call evaluate(mesh,sol,mesh%cell(cell2)%polTest,cell2,edge%X_gauss(p),edge%Y_gauss(p),u2)
                call boundary(flux,f_equa,cell2,cell1, &
                     u2,mesh,sol,j,p,edge%boundType,edge%bound,dir,flux_temp(p,:))
                flux_temp2(p,:)=gauss_weight(p)*flux_temp(p,:)*dt/(lengthN2*2.0_dp**(sub2+1))
                mesh%edge(j)%flux(irk,2,:)=mesh%edge(j)%flux(irk,2,:)+flux_temp2(p,:)
                soltemp%val(cell2,:)=soltemp%val(cell2,:)+flux_temp2(p,:)                       
             elseif (cell2<0) then
                call evaluate(mesh,sol,mesh%cell(cell1)%polTest,cell1,edge%X_gauss(p),edge%Y_gauss(p),u1)
                call boundary(flux,f_equa,cell1,cell2, &
                     u1,mesh,sol,j,p,edge%boundType,edge%bound,dir+2,flux_temp(p,:))
                flux_temp1(p,:)=-gauss_weight(p)*flux_temp(p,:)*dt/(lengthN1*2.0_dp**(sub1+1))
                mesh%edge(j)%flux(irk,1,:)=mesh%edge(j)%flux(irk,1,:)+flux_temp1(p,:)
                soltemp%val(cell1,:)=soltemp%val(cell1,:)+flux_temp1(p,:)
             else
                call evaluate(mesh,sol,mesh%cell(cell1)%polTest,cell1,edge%X_gauss(p),edge%Y_gauss(p),u1)
                call evaluate(mesh,sol,mesh%cell(cell2)%polTest,cell2,edge%X_gauss(p),edge%Y_gauss(p),u2)
                call flux(u1,u2,f_equa,dir,flux_temp(p,:))
                flux_temp1(p,:)=-gauss_weight(p)*flux_temp(p,:)*dt/(lengthN1*2.0_dp**(sub1+1))
                flux_temp2(p,:)=gauss_weight(p)*flux_temp(p,:)*dt/(lengthN2*2.0_dp**(sub2+1))
                mesh%edge(j)%flux(irk,1,:)=mesh%edge(j)%flux(irk,1,:)+flux_temp1(p,:)
                mesh%edge(j)%flux(irk,2,:)=mesh%edge(j)%flux(irk,2,:)+flux_temp2(p,:)
                soltemp%val(cell1,:)=soltemp%val(cell1,:)+flux_temp1(p,:)
                soltemp%val(cell2,:)=soltemp%val(cell2,:)+flux_temp2(p,:)
             endif
          endif
       enddo

    enddo

    call detect(mesh,sol,sol2,str_equa,L_str_criteria,L_var_criteria,L_eps,gauss_weight,period,NOT_ACCEPTED_CELL)

    sol2%val=soltemp%val-sol%val

    deallocate(u1,u2,flux_temp,flux_temp1,flux_temp2)
    
    return
  end subroutine advance

  subroutine euler_exp(mesh,sol0,sol,soltemp,f_equa,flux,order,t,dt,n, &
       gauss_weight,period,order_pc,NOT_ACCEPTED_EDGE,str_equa, &
       L_str_criteria,L_var_criteria,L_eps,NOT_ACCEPTED_CELL)
    type(meshStruct), intent(inout) :: mesh
    type(solStruct), intent(in) :: sol0,sol
    type(solStruct), intent(inout) :: soltemp
    procedure (sub_f), pointer, intent(in) :: f_equa
    procedure (sub_flux), pointer, intent(in) :: flux
    integer, intent(in) :: order,n
    real(dp), intent(in) :: dt
    real(dp), intent(inout) :: t
    real(dp), dimension(:), intent(in) :: gauss_weight
    logical, dimension(2), intent(in) :: period
    integer, dimension(:), intent(inout) :: order_pc
    integer, dimension(:), intent(in) :: NOT_ACCEPTED_EDGE
    character(len=20), intent(in) :: str_equa
    integer, dimension(:), intent(in) :: L_var_criteria
    character(len=20), dimension(:), intent(in) :: L_str_criteria
    real(dp), dimension(:), intent(in) :: L_eps
    integer, dimension(:), intent(inout) :: NOT_ACCEPTED_CELL
    type(solStruct) :: Fsol
    integer :: i

    allocate(Fsol%val(mesh%nc,sol%nvar))   

    call advance(mesh,sol0,Fsol,f_equa,flux,order,dt,t,n,1, &
         gauss_weight,period,order_pc,NOT_ACCEPTED_EDGE,str_equa, &
         L_str_criteria,L_var_criteria,L_eps,NOT_ACCEPTED_CELL)
    soltemp%user=Fsol%user
    soltemp%val=sol%val+Fsol%val

    do i=1,size(mesh%edge)
       mesh%edge(i)%flux(2,:,:)=mesh%edge(i)%flux(1,:,:)
    enddo

    deallocate(Fsol%val)

    return
  end subroutine euler_exp

  subroutine SSPRK2(mesh,sol0,sol,soltemp,f_equa,flux,order,t,dt,n, &
       gauss_weight,period,order_pc,NOT_ACCEPTED_EDGE,str_equa, &
       L_str_criteria,L_var_criteria,L_eps,NOT_ACCEPTED_CELL)
    type(meshStruct), intent(inout) :: mesh
    type(solStruct), intent(in) :: sol0,sol
    type(solStruct), intent(inout) :: soltemp
    procedure (sub_f), pointer, intent(in) :: f_equa
    procedure (sub_flux), pointer, intent(in) :: flux
    integer, intent(in) :: order,n
    real(dp), intent(in) :: dt
    real(dp), intent(inout) :: t
    real(dp), dimension(:), intent(in) :: gauss_weight
    logical, dimension(2), intent(in) :: period
    integer, dimension(:), intent(inout) :: order_pc
    integer, dimension(:), intent(in) :: NOT_ACCEPTED_EDGE
    character(len=20), intent(in) :: str_equa
    integer, dimension(:), intent(in) :: L_var_criteria
    character(len=20), dimension(:), intent(in) :: L_str_criteria
    real(dp), dimension(:), intent(in) :: L_eps
    integer, dimension(:), intent(inout) :: NOT_ACCEPTED_CELL
    type(solStruct) :: sol1,Fsol,Fsol1
    real(dp), dimension(:,:), allocatable :: f1
    integer :: i

    allocate(sol1%val(mesh%nc,sol%nvar),Fsol%val(mesh%nc,sol%nvar),Fsol1%val(mesh%nc,sol%nvar))
    allocate(f1(2,sol%nvar))
    
    sol1=sol

    call advance(mesh,sol0,Fsol,f_equa,flux,order,dt,t,n,1, &
         gauss_weight,period,order_pc,NOT_ACCEPTED_EDGE,str_equa, &
         L_str_criteria,L_var_criteria,L_eps,NOT_ACCEPTED_CELL)
    sol1%user=Fsol%user
    sol1%val=sol%val+Fsol%val
    call advance(mesh,sol1,Fsol1,f_equa,flux,order,dt,t,n,2, &
         gauss_weight,period,order_pc,NOT_ACCEPTED_EDGE,str_equa, &
         L_str_criteria,L_var_criteria,L_eps,NOT_ACCEPTED_CELL)
    soltemp%user=Fsol1%user
    soltemp%val=0.5_dp*sol%val+0.5_dp*sol1%val+0.5_dp*Fsol1%val

    do i=1,size(mesh%edge)
       f1=mesh%edge(i)%flux(1,:,:)
       mesh%edge(i)%flux(3,:,:)=0.5_dp*f1+0.5_dp*mesh%edge(i)%flux(2,:,:)
    enddo

    deallocate(sol1%val,Fsol%val,Fsol1%val,f1)

    return
  end subroutine SSPRK2

  subroutine SSPRK3(mesh,sol0,sol,soltemp,f_equa,flux,order,t,dt,n, &
       gauss_weight,period,order_pc,NOT_ACCEPTED_EDGE,str_equa, &
       L_str_criteria,L_var_criteria,L_eps,NOT_ACCEPTED_CELL)
    type(meshStruct), intent(inout) :: mesh
    type(solStruct), intent(in) :: sol0,sol
    type(solStruct), intent(inout) :: soltemp
    procedure (sub_f), pointer, intent(in) :: f_equa
    procedure (sub_flux), pointer, intent(in) :: flux
    integer, intent(in) :: order,n
    real(dp), intent(in) :: dt
    real(dp), intent(inout) :: t
    real(dp), dimension(:), intent(in) :: gauss_weight
    logical, dimension(2), intent(in) :: period
    integer, dimension(:), intent(inout) :: order_pc
    integer, dimension(:), intent(in) :: NOT_ACCEPTED_EDGE
    character(len=20), intent(in) :: str_equa
    integer, dimension(:), intent(in) :: L_var_criteria
    character(len=20), dimension(:), intent(in) :: L_str_criteria
    real(dp), dimension(:), intent(in) :: L_eps
    integer, dimension(:), intent(inout) :: NOT_ACCEPTED_CELL
    type(solStruct) :: sol1,sol2,Fsol,Fsol1,Fsol2
    real(dp), dimension(:,:), allocatable :: f1,f2
    integer :: i

    allocate(sol1%val(mesh%nc,sol%nvar),sol2%val(mesh%nc,sol%nvar))
    allocate(Fsol%val(mesh%nc,sol%nvar),Fsol1%val(mesh%nc,sol%nvar),Fsol2%val(mesh%nc,sol%nvar))
    allocate(f1(2,sol%nvar),f2(2,sol%nvar))
    
    sol1=sol
    sol2=sol

    call advance(mesh,sol0,Fsol,f_equa,flux,order,dt,t,n,1, &
         gauss_weight,period,order_pc,NOT_ACCEPTED_EDGE,str_equa, &
         L_str_criteria,L_var_criteria,L_eps,NOT_ACCEPTED_CELL)
    sol1%user=Fsol%user
    sol1%val=sol%val+Fsol%val
    call advance(mesh,sol1,Fsol1,f_equa,flux,order,dt,t,n,2, &
         gauss_weight,period,order_pc,NOT_ACCEPTED_EDGE,str_equa, &
         L_str_criteria,L_var_criteria,L_eps,NOT_ACCEPTED_CELL)
    sol2%user=Fsol1%user
    sol2%val=0.75_dp*sol%val+0.25_dp*sol1%val+0.25*Fsol1%val
    call advance(mesh,sol2,Fsol2,f_equa,flux,order,dt,t,n,3, &
         gauss_weight,period,order_pc,NOT_ACCEPTED_EDGE,str_equa, &
         L_str_criteria,L_var_criteria,L_eps,NOT_ACCEPTED_CELL)
    soltemp%user=Fsol2%user
    soltemp%val=1.0_dp/3.0_dp*sol%val+2.0_dp/3.0_dp*sol2%val+2.0_dp/3.0_dp*Fsol2%val

    do i=1,size(mesh%edge)
       f1=mesh%edge(i)%flux(1,:,:)
       f2=0.25_dp*f1+0.25_dp*mesh%edge(i)%flux(2,:,:)
       mesh%edge(i)%flux(4,:,:)=2.0_dp/3.0_dp*f2+2.0_dp/3.0_dp*mesh%edge(i)%flux(3,:,:)
    enddo

    deallocate(sol1%val,sol2%val,Fsol%val,Fsol1%val,Fsol2%val,f1,f2)

    return
  end subroutine SSPRK3

  subroutine SSPRK4(mesh,sol0,sol,soltemp,f_equa,flux,order,t,dt,n, &
       gauss_weight,period,order_pc,NOT_ACCEPTED_EDGE,str_equa, &
       L_str_criteria,L_var_criteria,L_eps,NOT_ACCEPTED_CELL)
    type(meshStruct), intent(inout) :: mesh
    type(solStruct), intent(in) :: sol0,sol
    type(solStruct), intent(inout) :: soltemp
    procedure (sub_f), pointer, intent(in) :: f_equa
    procedure (sub_flux), pointer, intent(in) :: flux
    integer, intent(in) :: order,n
    real(dp), intent(in) :: dt
    real(dp), intent(inout) :: t
    real(dp), dimension(:), intent(in) :: gauss_weight
    logical, dimension(2), intent(in) :: period
    integer, dimension(:), intent(inout) :: order_pc
    integer, dimension(:), intent(in) :: NOT_ACCEPTED_EDGE
    character(len=20), intent(in) :: str_equa
    integer, dimension(:), intent(in) :: L_var_criteria
    character(len=20), dimension(:), intent(in) :: L_str_criteria
    real(dp), dimension(:), intent(in) :: L_eps
    integer, dimension(:), intent(inout) :: NOT_ACCEPTED_CELL
    type(solStruct) :: sol1,sol2,sol3,sol4,Fsol,Fsol1,Fsol2,Fsol3,Fsol4
    real(dp), dimension(:,:), allocatable :: f1,f2,f3,f4
    integer :: i
    
    allocate(sol1%val(mesh%nc,sol%nvar),sol2%val(mesh%nc,sol%nvar),sol3%val(mesh%nc,sol%nvar),sol4%val(mesh%nc,sol%nvar))
    allocate(Fsol%val(mesh%nc,sol%nvar),Fsol1%val(mesh%nc,sol%nvar),Fsol2%val(mesh%nc,sol%nvar))
    allocate(Fsol3%val(mesh%nc,sol%nvar),Fsol4%val(mesh%nc,sol%nvar))
    allocate(f1(2,sol%nvar),f2(2,sol%nvar),f3(2,sol%nvar),f4(2,sol%nvar))

    sol1=sol
    sol2=sol
    sol3=sol
    sol4=sol

    call advance(mesh,sol0,Fsol,f_equa,flux,order,dt,t,n,1, &
         gauss_weight,period,order_pc,NOT_ACCEPTED_EDGE,str_equa, &
         L_str_criteria,L_var_criteria,L_eps,NOT_ACCEPTED_CELL)
    sol1%user=Fsol%user
    sol1%val=sol%val+0.391752226571890_dp*Fsol%val
    call advance(mesh,sol1,Fsol1,f_equa,flux,order,dt,t,n,2, &
         gauss_weight,period,order_pc,NOT_ACCEPTED_EDGE,str_equa, &
         L_str_criteria,L_var_criteria,L_eps,NOT_ACCEPTED_CELL)
    sol2%user=Fsol1%user
    sol2%val=0.444370493651235_dp*sol%val+0.555629506348765_dp*sol1%val+0.368410593050371_dp*Fsol1%val   
    call advance(mesh,sol2,Fsol2,f_equa,flux,order,dt,t,n,3, &
         gauss_weight,period,order_pc,NOT_ACCEPTED_EDGE,str_equa, &
         L_str_criteria,L_var_criteria,L_eps,NOT_ACCEPTED_CELL)
    sol3%user=Fsol2%user
    sol3%val=0.620101851488403_dp*sol%val+0.379898148511597_dp*sol2%val+0.251891774271694_dp*Fsol2%val
    call advance(mesh,sol3,Fsol3,f_equa,flux,order,dt,t,n,4, &
         gauss_weight,period,order_pc,NOT_ACCEPTED_EDGE,str_equa, &
         L_str_criteria,L_var_criteria,L_eps,NOT_ACCEPTED_CELL)
    sol4%user=Fsol3%user
    sol4%val=0.178079954393132_dp*sol%val+0.821920045606868_dp*sol3%val+0.544974750228521_dp*Fsol3%val
    call advance(mesh,sol4,Fsol4,f_equa,flux,order,dt,t,n,5, &
         gauss_weight,period,order_pc,NOT_ACCEPTED_EDGE,str_equa, &
         L_str_criteria,L_var_criteria,L_eps,NOT_ACCEPTED_CELL)
    soltemp%user=Fsol4%user
    soltemp%val=0.517231671970585_dp*sol2%val+0.096059710526147_dp*sol3%val+0.386708617503269_dp*sol4%val+ &
         0.063692468666290_dp*Fsol3%val+0.226007483236906_dp*Fsol4%val

    do i=1,size(mesh%edge)
       f1=0.391752226571890_dp*mesh%edge(i)%flux(1,:,:)
       f2=0.555629506348765_dp*f1+0.368410593050371_dp*mesh%edge(i)%flux(2,:,:)
       f3=0.379898148511597_dp*f2+0.251891774271694_dp*mesh%edge(i)%flux(3,:,:)
       f4=0.821920045606868_dp*f3+0.544974750228521_dp*mesh%edge(i)%flux(4,:,:)
       mesh%edge(i)%flux(6,:,:)=0.517231671970585_dp*f2+0.096059710526147_dp*f3+0.386708617503269_dp*f4+ &
         0.063692468666290_dp*mesh%edge(i)%flux(4,:,:)+0.226007483236906_dp*mesh%edge(i)%flux(5,:,:)
    enddo

    deallocate(sol1%val,sol2%val,sol3%val,sol4%val)
    deallocate(Fsol%val,Fsol1%val,Fsol2%val,Fsol3%val,Fsol4%val)
    deallocate(f1,f2,f3,f4)

    return
  end subroutine SSPRK4

  subroutine SSPRK5(mesh,sol0,sol,soltemp,f_equa,flux,order,t,dt,n, &
       gauss_weight,period,order_pc,NOT_ACCEPTED_EDGE,str_equa, &
       L_str_criteria,L_var_criteria,L_eps,NOT_ACCEPTED_CELL)
    type(meshStruct), intent(inout) :: mesh
    type(solStruct), intent(in) :: sol0,sol
    type(solStruct), intent(inout) :: soltemp
    procedure (sub_f), pointer, intent(in) :: f_equa
    procedure (sub_flux), pointer, intent(in) :: flux
    integer, intent(in) :: order,n
    real(dp), intent(in) :: dt
    real(dp), intent(inout) :: t
    real(dp), dimension(:), intent(in) :: gauss_weight
    logical, dimension(2), intent(in) :: period
    integer, dimension(:), intent(inout) :: order_pc
    integer, dimension(:), intent(in) :: NOT_ACCEPTED_EDGE
    character(len=20), intent(in) :: str_equa
    integer, dimension(:), intent(in) :: L_var_criteria
    character(len=20), dimension(:), intent(in) :: L_str_criteria
    real(dp), dimension(:), intent(in) :: L_eps
    integer, dimension(:), intent(inout) :: NOT_ACCEPTED_CELL
    type(solStruct) :: sol1,sol2,sol3,sol4,sol5,sol6,sol7,sol8,sol9
    type(solStruct) :: Fsol,Fsol1,Fsol2,Fsol3,Fsol4,Fsol5,Fsol6,Fsol7,Fsol8,Fsol9
    real(dp), dimension(:,:), allocatable :: f1,f2,f3,f4,f5,f6,f7,f8,f9
    integer :: i
    
    allocate(sol1%val(mesh%nc,sol%nvar),sol2%val(mesh%nc,sol%nvar),sol3%val(mesh%nc,sol%nvar))
    allocate(sol4%val(mesh%nc,sol%nvar),sol5%val(mesh%nc,sol%nvar),sol6%val(mesh%nc,sol%nvar))
    allocate(sol7%val(mesh%nc,sol%nvar),sol8%val(mesh%nc,sol%nvar),sol9%val(mesh%nc,sol%nvar))
    allocate(Fsol%val(mesh%nc,sol%nvar),Fsol1%val(mesh%nc,sol%nvar),Fsol2%val(mesh%nc,sol%nvar))
    allocate(Fsol3%val(mesh%nc,sol%nvar),Fsol4%val(mesh%nc,sol%nvar),Fsol5%val(mesh%nc,sol%nvar))
    allocate(Fsol6%val(mesh%nc,sol%nvar),Fsol7%val(mesh%nc,sol%nvar),Fsol8%val(mesh%nc,sol%nvar),Fsol9%val(mesh%nc,sol%nvar))
    allocate(f1(2,sol%nvar),f2(2,sol%nvar),f3(2,sol%nvar),f4(2,sol%nvar),f5(2,sol%nvar))
    allocate(f6(2,sol%nvar),f7(2,sol%nvar),f8(2,sol%nvar),f9(2,sol%nvar))

    sol1=sol
    sol2=sol
    sol3=sol
    sol4=sol
    sol5=sol
    sol6=sol
    sol7=sol
    sol8=sol
    sol9=sol

    call advance(mesh,sol0,Fsol,f_equa,flux,order,dt,t,n,1, &
         gauss_weight,period,order_pc,NOT_ACCEPTED_EDGE,str_equa, &
         L_str_criteria,L_var_criteria,L_eps,NOT_ACCEPTED_CELL)
    sol1%user=Fsol%user
    sol1%val=sol%val+0.173586107937995_dp*Fsol%val
    call advance(mesh,sol1,Fsol1,f_equa,flux,order,dt,t,n,2, &
         gauss_weight,period,order_pc,NOT_ACCEPTED_EDGE,str_equa, &
         L_str_criteria,L_var_criteria,L_eps,NOT_ACCEPTED_CELL)
    sol2%user=Fsol1%user
    sol2%val=0.258168167463650_dp*sol%val+0.741831832536350_dp*sol1%val+0.218485490268790_dp*Fsol1%val  
    call advance(mesh,sol2,Fsol2,f_equa,flux,order,dt,t,n,3, &
         gauss_weight,period,order_pc,NOT_ACCEPTED_EDGE,str_equa, &
         L_str_criteria,L_var_criteria,L_eps,NOT_ACCEPTED_CELL)
    sol3%user=Fsol2%user
    sol3%val=0.037493531856076_dp*sol1%val+0.962506468143924_dp*sol2%val+ &
         0.011042654588541_dp*Fsol1%val+0.283478934653295_dp*Fsol2%val
    call advance(mesh,sol3,Fsol3,f_equa,flux,order,dt,t,n,4, &
         gauss_weight,period,order_pc,NOT_ACCEPTED_EDGE,str_equa, &
         L_str_criteria,L_var_criteria,L_eps,NOT_ACCEPTED_CELL)
    sol4%user=Fsol3%user
    sol4%val=0.595955269449077_dp*sol%val+0.404044730550923_dp*sol2%val+0.118999896166647_dp*Fsol2%val
    call advance(mesh,sol4,Fsol4,f_equa,flux,order,dt,t,n,5, &
         gauss_weight,period,order_pc,NOT_ACCEPTED_EDGE,str_equa, &
         L_str_criteria,L_var_criteria,L_eps,NOT_ACCEPTED_CELL)
    sol5%user=Fsol4%user
    sol5%val=0.331848124368345_dp*sol%val+0.008466192609453_dp*sol3%val+0.659685683022202_dp*sol4%val+ &
         0.025030881091201_dp*Fsol%val-0.002493476502164_dp*Fsol3%val+0.194291675763785_dp*Fsol4%val
    call advance(mesh,sol5,Fsol5,f_equa,flux,order,dt,t,n,6, &
         gauss_weight,period,order_pc,NOT_ACCEPTED_EDGE,str_equa, &
         L_str_criteria,L_var_criteria,L_eps,NOT_ACCEPTED_CELL)
    sol6%user=Fsol5%user
    sol6%val=0.086976414344414_dp*sol%val+0.913023585655586_dp*sol5%val+0.268905157462563_dp*Fsol5%val
    call advance(mesh,sol6,Fsol6,f_equa,flux,order,dt,t,n,7, &
         gauss_weight,period,order_pc,NOT_ACCEPTED_EDGE,str_equa, &
         L_str_criteria,L_var_criteria,L_eps,NOT_ACCEPTED_CELL)
    sol7%user=Fsol6%user
    sol7%val=0.075863700003186_dp*sol%val+0.267513039663395_dp*sol2%val+0.656623260333419_dp*sol6%val+ &
         0.066115378914543_dp*Fsol2%val+0.193389726166555_dp*Fsol6%val
    call advance(mesh,sol7,Fsol7,f_equa,flux,order,dt,t,n,8, &
         gauss_weight,period,order_pc,NOT_ACCEPTED_EDGE,str_equa, &
         L_str_criteria,L_var_criteria,L_eps,NOT_ACCEPTED_CELL)
    sol8%user=Fsol7%user
    sol8%val=0.005212058095597_dp*sol%val+0.407430107306541_dp*sol3%val+0.587357834597862_dp*sol7%val- &
         0.119996962708895_dp*Fsol3%val+0.172989562899406_dp*Fsol7%val
    call advance(mesh,sol8,Fsol8,f_equa,flux,order,dt,t,n,9, &
         gauss_weight,period,order_pc,NOT_ACCEPTED_EDGE,str_equa, &
         L_str_criteria,L_var_criteria,L_eps,NOT_ACCEPTED_CELL)
    sol9%user=Fsol8%user
    sol9%val=0.122832051947995_dp*sol%val+0.877167948052005_dp*sol8%val+ &
         0.000000000000035_dp*Fsol%val+0.258344898092277_dp*Fsol8%val
    call advance(mesh,sol9,Fsol9,f_equa,flux,order,dt,t,n,10, &
         gauss_weight,period,order_pc,NOT_ACCEPTED_EDGE,str_equa, &
         L_str_criteria,L_var_criteria,L_eps,NOT_ACCEPTED_CELL)
    soltemp%user=Fsol9%user
    soltemp%val=0.075346276482673_dp*sol%val+0.000425904246091_dp*sol1%val+0.064038648145995_dp*sol5%val+ &
         0.354077936287492_dp*sol6%val+0.506111234837749_dp*sol9%val+0.016982542367506_dp*Fsol%val+ &
         0.018860764424857_dp*Fsol5%val+0.098896719553054_dp*Fsol6%val+0.149060685217562_dp*Fsol9%val

    do i=1,size(mesh%edge)
       f1=0.173586107937995_dp*mesh%edge(i)%flux(1,:,:)
       f2=0.741831832536350_dp*f1+0.218485490268790_dp*mesh%edge(i)%flux(2,:,:)
       f3=0.037493531856076_dp*f1+0.962506468143924_dp*f2+0.011042654588541_dp*mesh%edge(i)%flux(2,:,:)+ &
         0.283478934653295_dp*mesh%edge(i)%flux(3,:,:)
       f4=0.404044730550923_dp*f2+0.118999896166647_dp*mesh%edge(i)%flux(3,:,:)
       f5=0.008466192609453_dp*f3+0.659685683022202_dp*f4+0.025030881091201_dp*mesh%edge(i)%flux(1,:,:)- &
         0.002493476502164_dp*mesh%edge(i)%flux(4,:,:)+0.194291675763785_dp*mesh%edge(i)%flux(5,:,:)
       f6=0.913023585655586_dp*f5+0.268905157462563_dp*mesh%edge(i)%flux(6,:,:)
       f7=0.267513039663395_dp*f2+0.656623260333419_dp*f6+0.066115378914543_dp*mesh%edge(i)%flux(3,:,:)+ &
         0.193389726166555_dp*mesh%edge(i)%flux(7,:,:)
       f8=0.407430107306541_dp*f3+0.587357834597862_dp*f7-0.119996962708895_dp*mesh%edge(i)%flux(4,:,:)+ &
         0.172989562899406_dp*mesh%edge(i)%flux(8,:,:)
       f9=0.877167948052005_dp*f8+0.000000000000035_dp*mesh%edge(i)%flux(1,:,:)+ &
         0.258344898092277_dp*mesh%edge(i)%flux(9,:,:)
       mesh%edge(i)%flux(11,:,:)=0.000425904246091_dp*f1+0.064038648145995_dp*f5+0.354077936287492_dp*f6+ &
         0.506111234837749_dp*f9+0.016982542367506_dp*mesh%edge(i)%flux(1,:,:)+ &
         0.018860764424857_dp*mesh%edge(i)%flux(6,:,:)+0.098896719553054_dp*mesh%edge(i)%flux(7,:,:)+ &
         0.149060685217562_dp*mesh%edge(i)%flux(10,:,:)
    enddo

    deallocate(sol1%val,sol2%val,sol3%val,sol4%val,sol5%val,sol6%val,sol7%val,sol8%val,sol9%val)
    deallocate(Fsol%val,Fsol1%val,Fsol2%val,Fsol3%val,Fsol4%val,Fsol5%val,Fsol6%val,Fsol7%val,Fsol8%val,Fsol9%val)
    deallocate(f1,f2,f3,f4,f5,f6,f7,f8,f9)
    
    return
  end subroutine SSPRK5

end module time
