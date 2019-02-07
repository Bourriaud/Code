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

  subroutine advance(mesh,sol,sol2,str_equa,f_equa,flux,order,dt,n, &
       L_str_criteria,L_var_criteria,L_eps,gauss_weight,verbosity)
    type(meshStruct), intent(inout) :: mesh
    type(solStruct), intent(in) :: sol
    type(solStruct), intent(inout) :: sol2
    character(len=20), intent(in) :: str_equa
    procedure (sub_f), pointer, intent(in) :: f_equa
    procedure (sub_flux), pointer, intent(in) :: flux
    integer, intent(in) :: order,n,verbosity
    real(dp), intent(in) :: dt
    character(len=20), dimension(:), intent(in) :: L_str_criteria
    integer, dimension(:), intent(in) :: L_var_criteria
    real(dp), dimension(:), intent(in) :: L_eps
    real(dp), dimension(:), intent(in) :: gauss_weight
    integer :: k,i,j,p,cell1,cell2,dir,count,deg,sub1,sub2,ac1,ac2
    real(dp), dimension(:), allocatable :: u1,u2
    type(edgeStruct) :: edge
    integer, dimension(:), allocatable :: NOT_ACCEPTED_CELL,NOT_ACCEPTED_EDGE
    integer, dimension(:), allocatable :: NAC_cycle,NAC_reason
    type(solStruct) :: soltemp
    integer, dimension(2) :: L_deg
    real(dp) :: lengthN1,lengthN2
    if(.false.)print*,n

    allocate(u1(sol%nvar),u2(sol%nvar))
    allocate(NOT_ACCEPTED_CELL(mesh%nc),NOT_ACCEPTED_EDGE(mesh%ne),soltemp%val(mesh%nc,sol%nvar))
    allocate(NAC_cycle(mesh%nc),NAC_reason(mesh%nc))
    
    sol2=sol
    soltemp=sol
    NAC_cycle=0
    
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
          ac1=abs(cell1)
          ac2=abs(cell2)
          dir=edge%dir
          lengthN1=mesh%cell(ac1)%dx
          lengthN2=mesh%cell(ac2)%dx
          sub1=edge%sub(1)
          sub2=edge%sub(2)
          deg=min(mesh%cell(ac1)%deg,mesh%cell(ac2)%deg)
          mesh%edge(j)%deg=deg
          call reconstruct(mesh,sol,ac1,deg+1,gauss_weight)
          call reconstruct(mesh,sol,ac2,deg+1,gauss_weight)
          if (deg+1==order) then
             mesh%cell(ac1)%polMax=mesh%cell(ac1)%polCoef
             mesh%cell(ac2)%polMax=mesh%cell(ac2)%polCoef
          endif
          do p=1,order
             if (.not.edge%flux_acc(p)) then
                if (cell1<0) then
                   call evaluate(mesh,sol,cell2,deg+1,gauss_weight,edge%X_gauss(p),edge%Y_gauss(p),u2)
                   call boundary(flux,f_equa,gauss_weight,cell2,cell1, &
                        u2,mesh,sol,j,p,edge%boundType,edge%bound,dir,deg+1,mesh%edge(j)%flux(p,:))
                   soltemp%val(cell2,:)=soltemp%val(cell2,:)+gauss_weight(p)*mesh%edge(j)%flux(p,:)*dt/(lengthN2*2.0_dp**(sub2+1))
                elseif (cell2<0) then
                   call evaluate(mesh,sol,cell1,deg+1,gauss_weight,edge%X_gauss(p),edge%Y_gauss(p),u1)
                   call boundary(flux,f_equa,gauss_weight,cell1,cell2, &
                        u1,mesh,sol,j,p,edge%boundType,edge%bound,dir+2,deg+1,mesh%edge(j)%flux(p,:))
                   soltemp%val(cell1,:)=soltemp%val(cell1,:)-gauss_weight(p)*mesh%edge(j)%flux(p,:)*dt/(lengthN1*2.0_dp**(sub1+1))
                else
                   call evaluate(mesh,sol,cell1,deg+1,gauss_weight,edge%X_gauss(p),edge%Y_gauss(p),u1)
                   call evaluate(mesh,sol,cell2,deg+1,gauss_weight,edge%X_gauss(p),edge%Y_gauss(p),u2)
                   call flux(u1,u2,f_equa,dir,mesh%edge(j)%flux(p,:))
                   soltemp%val(cell1,:)=soltemp%val(cell1,:)-gauss_weight(p)*mesh%edge(j)%flux(p,:)*dt/(lengthN1*2.0_dp**(sub1+1))
                   soltemp%val(cell2,:)=soltemp%val(cell2,:)+gauss_weight(p)*mesh%edge(j)%flux(p,:)*dt/(lengthN2*2.0_dp**(sub2+1))
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
       call decrement(mesh,sol,soltemp,str_equa,deg,dt,L_str_criteria,L_var_criteria,L_eps, &
            gauss_weight,NOT_ACCEPTED_CELL,NOT_ACCEPTED_EDGE,NAC_reason,verbosity)

       if (verbosity>1) then
          do k=1,size(NOT_ACCEPTED_CELL)
             NAC_cycle(NOT_ACCEPTED_CELL(k))=NAC_cycle(NOT_ACCEPTED_CELL(k))+1
          enddo
          if (size(NOT_ACCEPTED_CELL)==0) then
             call write_accept(mesh,NAC_cycle,NAC_reason,n)
          endif
       endif

    enddo

    sol2%val=soltemp%val-sol%val

    deallocate(u1,u2,NOT_ACCEPTED_CELL,NOT_ACCEPTED_EDGE,soltemp%val,NAC_cycle,NAC_reason)
    
    return
  end subroutine advance

  subroutine euler_exp(mesh,sol,str_equa,f_equa,flux,speed,order,cfl,t,n,tf, &
       L_str_criteria,L_var_criteria,L_eps,gauss_weight,verbosity)
    type(meshStruct), intent(inout) :: mesh
    type(solStruct), intent(inout) :: sol
    character(len=20), intent(in) :: str_equa
    procedure (sub_f), pointer, intent(in) :: f_equa
    procedure (sub_flux), pointer, intent(in) :: flux
    procedure (sub_speed), pointer, intent(in) :: speed
    integer, intent(in) :: order,n,verbosity
    real(dp), intent(in) :: cfl,tf
    real(dp), intent(inout) :: t
    character(len=20), dimension(:), intent(in) :: L_str_criteria
    integer, dimension(:), intent(in) :: L_var_criteria
    real(dp), dimension(:), intent(in) :: L_eps
    real(dp), dimension(:), intent(in) :: gauss_weight
    type(solStruct) :: Fsol
    real(dp) :: dt

    allocate(Fsol%val(mesh%nc,sol%nvar))

    call compute_timestep(mesh,sol,f_equa,speed,cfl,tf,t,dt)

    call advance(mesh,sol,Fsol,str_equa,f_equa,flux,order,dt,n, &
         L_str_criteria,L_var_criteria,L_eps,gauss_weight,verbosity)
    sol%val=sol%val+Fsol%val

    t=t+dt

    deallocate(Fsol%val)

    return
  end subroutine euler_exp

  subroutine SSPRK2(mesh,sol,str_equa,f_equa,flux,speed,order,cfl,t,n,tf, &
       L_str_criteria,L_var_criteria,L_eps,gauss_weight,verbosity)
    type(meshStruct), intent(inout) :: mesh
    type(solStruct), intent(inout) :: sol
    character(len=20), intent(in) :: str_equa
    procedure (sub_f), pointer, intent(in) :: f_equa
    procedure (sub_flux), pointer, intent(in) :: flux
    procedure (sub_speed), pointer, intent(in) :: speed
    integer, intent(in) :: order,n,verbosity
    real(dp), intent(in) :: cfl,tf
    real(dp), intent(inout) :: t
    character(len=20), dimension(:), intent(in) :: L_str_criteria
    integer, dimension(:), intent(in) :: L_var_criteria
    real(dp), dimension(:), intent(in) :: L_eps
    real(dp), dimension(:), intent(in) :: gauss_weight
    type(solStruct) :: sol1,Fsol,Fsol1
    real(dp) :: dt

    allocate(sol1%val(mesh%nc,sol%nvar),Fsol%val(mesh%nc,sol%nvar),Fsol1%val(mesh%nc,sol%nvar))
    sol1=sol

    call compute_timestep(mesh,sol,f_equa,speed,cfl,tf,t,dt)

    call advance(mesh,sol,Fsol,str_equa,f_equa,flux,order,dt,n, &
         L_str_criteria,L_var_criteria,L_eps,gauss_weight,verbosity)
    sol1%val=sol%val+Fsol%val
    call advance(mesh,sol1,Fsol1,str_equa,f_equa,flux,order,dt,n, &
         L_str_criteria,L_var_criteria,L_eps,gauss_weight,verbosity)
    sol%val=0.5_dp*sol%val+0.5_dp*sol1%val+0.5_dp*Fsol1%val

    t=t+dt

    deallocate(sol1%val,Fsol%val,Fsol1%val)

    return
  end subroutine SSPRK2

  subroutine SSPRK3(mesh,sol,str_equa,f_equa,flux,speed,order,cfl,t,n,tf,L_str_criteria,L_var_criteria,L_eps,gauss_weight,verbosity)
    type(meshStruct), intent(inout) :: mesh
    type(solStruct), intent(inout) :: sol
    character(len=20), intent(in) :: str_equa
    procedure (sub_f), pointer, intent(in) :: f_equa
    procedure (sub_flux), pointer, intent(in) :: flux
    procedure (sub_speed), pointer, intent(in) :: speed
    integer, intent(in) :: order,n,verbosity
    real(dp), intent(in) :: cfl,tf
    real(dp), intent(inout) :: t
    character(len=20), dimension(:), intent(in) :: L_str_criteria
    integer, dimension(:), intent(in) :: L_var_criteria
    real(dp), dimension(:), intent(in) :: L_eps
    real(dp), dimension(:), intent(in) :: gauss_weight
    type(solStruct) :: sol1,sol2,Fsol,Fsol1,Fsol2
    real(dp) :: dt

    allocate(sol1%val(mesh%nc,sol%nvar),sol2%val(mesh%nc,sol%nvar))
    allocate(Fsol%val(mesh%nc,sol%nvar),Fsol1%val(mesh%nc,sol%nvar),Fsol2%val(mesh%nc,sol%nvar))
    sol1=sol
    sol2=sol

    call compute_timestep(mesh,sol,f_equa,speed,cfl,tf,t,dt)

    call advance(mesh,sol,Fsol,str_equa,f_equa,flux,order,dt,n, &
         L_str_criteria,L_var_criteria,L_eps,gauss_weight,verbosity)
    sol1%val=sol%val+Fsol%val
    call advance(mesh,sol1,Fsol1,str_equa,f_equa,flux,order,dt,n, &
         L_str_criteria,L_var_criteria,L_eps,gauss_weight,verbosity)
    sol2%val=0.75_dp*sol%val+0.25_dp*sol1%val+0.25*Fsol1%val
    call advance(mesh,sol2,Fsol2,str_equa,f_equa,flux,order,dt,n, &
         L_str_criteria,L_var_criteria,L_eps,gauss_weight,verbosity)
    sol%val=1.0_dp/3.0_dp*sol%val+2.0_dp/3.0_dp*sol2%val+2.0_dp/3.0_dp*Fsol2%val

    t=t+dt

    deallocate(sol1%val,sol2%val,Fsol%val,Fsol1%val,Fsol2%val)

    return
  end subroutine SSPRK3

  subroutine SSPRK4(mesh,sol,str_equa,f_equa,flux,speed,order,cfl,t,n,tf,L_str_criteria,L_var_criteria,L_eps,gauss_weight,verbosity)
    type(meshStruct), intent(inout) :: mesh
    type(solStruct), intent(inout) :: sol
    character(len=20), intent(in) :: str_equa
    procedure (sub_f), pointer, intent(in) :: f_equa
    procedure (sub_flux), pointer, intent(in) :: flux
    procedure (sub_speed), pointer, intent(in) :: speed
    integer, intent(in) :: order,n,verbosity
    real(dp), intent(in) :: cfl,tf
    real(dp), intent(inout) :: t
    character(len=20), dimension(:), intent(in) :: L_str_criteria
    integer, dimension(:), intent(in) :: L_var_criteria
    real(dp), dimension(:), intent(in) :: L_eps
    real(dp), dimension(:), intent(in) :: gauss_weight
    type(solStruct) :: sol1,sol2,sol3,sol4,Fsol,Fsol1,Fsol2,Fsol3,Fsol4
    real(dp) :: dt

    allocate(sol1%val(mesh%nc,sol%nvar),sol2%val(mesh%nc,sol%nvar),sol3%val(mesh%nc,sol%nvar),sol4%val(mesh%nc,sol%nvar))
    allocate(Fsol%val(mesh%nc,sol%nvar),Fsol1%val(mesh%nc,sol%nvar),Fsol2%val(mesh%nc,sol%nvar))
    allocate(Fsol3%val(mesh%nc,sol%nvar),Fsol4%val(mesh%nc,sol%nvar))

    call compute_timestep(mesh,sol,f_equa,speed,cfl,tf,t,dt)
    sol1=sol
    sol2=sol
    sol3=sol
    sol4=sol

    call advance(mesh,sol,Fsol,str_equa,f_equa,flux,order,dt,n, &
         L_str_criteria,L_var_criteria,L_eps,gauss_weight,verbosity)
    sol1%val=sol%val+0.391752226571890_dp*Fsol%val
    call advance(mesh,sol1,Fsol1,str_equa,f_equa,flux,order,dt,n, &
         L_str_criteria,L_var_criteria,L_eps,gauss_weight,verbosity)
    sol2%val=0.444370493651235_dp*sol%val+0.555629506348765_dp*sol1%val+0.368410593050371_dp*Fsol1%val   
    call advance(mesh,sol2,Fsol2,str_equa,f_equa,flux,order,dt,n, &
         L_str_criteria,L_var_criteria,L_eps,gauss_weight,verbosity)
    sol3%val=0.620101851488403_dp*sol%val+0.379898148511597_dp*sol2%val+0.251891774271694_dp*Fsol2%val
    call advance(mesh,sol3,Fsol3,str_equa,f_equa,flux,order,dt,n, &
         L_str_criteria,L_var_criteria,L_eps,gauss_weight,verbosity)
    sol4%val=0.178079954393132_dp*sol%val+0.821920045606868_dp*sol3%val+0.544974750228521_dp*Fsol3%val
    call advance(mesh,sol4,Fsol4,str_equa,f_equa,flux,order,dt,n, &
         L_str_criteria,L_var_criteria,L_eps,gauss_weight,verbosity)
    sol%val=0.517231671970585_dp*sol2%val+0.096059710526147_dp*sol3%val+0.386708617503269_dp*sol4%val+ &
         0.063692468666290_dp*Fsol3%val+0.226007483236906_dp*Fsol4%val

    t=t+dt

    deallocate(sol1%val,sol2%val,sol3%val,sol4%val)
    deallocate(Fsol%val,Fsol1%val,Fsol2%val,Fsol3%val,Fsol4%val)

    return
  end subroutine SSPRK4

end module time
