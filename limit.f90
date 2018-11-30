module limit

  use constant
  use types
  use reconstruction

  implicit none

contains

  subroutine decrement(mesh,sol,soltemp,str_equa,deg,dt,L_str_criteria,L_var_criteria,L_eps, &
       gauss_weight,NOT_ACCEPTED_CELL,NOT_ACCEPTED_EDGE)
    type(meshStruct), intent(inout) :: mesh
    type(solStruct), intent(in) :: sol
    type(solStruct), intent(inout) :: soltemp
    character(len=20), intent(in) :: str_equa
    integer, intent(in) :: deg
    real(dp), intent(in) :: dt
    integer, dimension(:), intent(in) :: L_var_criteria
    character(len=20), dimension(:), intent(in) :: L_str_criteria
    real(dp), dimension(:), intent(in) :: L_eps
    real(dp), dimension(:), intent(in) :: gauss_weight
    type(solStruct) :: sol2
    integer, dimension(:), allocatable, intent(inout) :: NOT_ACCEPTED_CELL,NOT_ACCEPTED_EDGE
    integer, dimension(:), allocatable :: NAC,NAE
    integer :: i,j,k,n,p,isol,nc,ne,cell1,cell2,edge
    procedure (sub_criteria), pointer :: criteria
    logical :: accept

    allocate(NAC(mesh%nc),NAE(mesh%ne))

    nc=0
    ne=0
    NAC=0
    NAE=0
    sol2=soltemp
    
    do n=1,size(L_str_criteria)
       
       isol=L_var_criteria(n)
       
       select case (trim(L_str_criteria(n)))
       case ('DMP')
          criteria => DMP
       case ('DMPu2')
          criteria => DMPu2
       case('PAD')
          criteria => PAD
       case default
          print*,trim(L_str_criteria(n))," criteria not implemented"
          call exit()
       end select

       do i=1,size(NOT_ACCEPTED_CELL)
          k=NOT_ACCEPTED_CELL(i)
          call criteria(mesh,sol,sol2,k,isol,L_eps(n),gauss_weight,str_equa,accept)
          if (.not.accept) then
             do j=1,size(mesh%cell(k)%edge)
                edge=mesh%cell(k)%edge(j)
                cell1=mesh%edge(edge)%cell1
                cell2=mesh%edge(edge)%cell2
                if (all(NAE/=edge).and.mesh%edge(edge)%deg>0) then
                   
                   select case (trim(mesh%edge(edge)%boundtype))
                      
                   case('PERIODIC')
                      ne=ne+1
                      NAE(ne)=edge
                      if (cell1*cell2<0) then
                         ne=ne+1
                         NAE(ne)=mesh%edge(edge)%period
                      endif
                      if (all(NAC/=abs(cell1))) then
                         mesh%cell(abs(cell1))%deg=deg
                         nc=nc+1
                         NAC(nc)=abs(cell1)
                      endif
                      if (all(NAC/=abs(cell2))) then
                         mesh%cell(abs(cell2))%deg=deg
                         nc=nc+1
                         NAC(nc)=abs(cell2)
                      endif
                      do p=1,size(gauss_weight)
                         call criteria_flux(mesh%edge(edge)%flux(p,:),mesh%edge(edge)%flux_acc(p))
                         if (.not.mesh%edge(edge)%flux_acc(p)) then
                            soltemp%val(abs(cell1),:)=soltemp%val(abs(cell1),:)+ &
                                 gauss_weight(p)*mesh%edge(edge)%flux(p,:)*dt/(mesh%edge(edge)%lengthN*2.0_dp)
                            soltemp%val(abs(cell2),:)=soltemp%val(abs(cell2),:)- &
                                 gauss_weight(p)*mesh%edge(edge)%flux(p,:)*dt/(mesh%edge(edge)%lengthN*2.0_dp)
                         endif
                      enddo
                      
                   case default
                      ne=ne+1
                      NAE(ne)=edge
                      if (all(NAC/=cell1).and.cell1>0) then
                         mesh%cell(cell1)%deg=deg
                         nc=nc+1
                         NAC(nc)=cell1
                      endif
                      if (all(NAC/=cell2).and.cell2>0) then
                         mesh%cell(cell2)%deg=deg
                         nc=nc+1
                         NAC(nc)=cell2
                      endif
                      do p=1,size(gauss_weight)
                         call criteria_flux(mesh%edge(edge)%flux(p,:),mesh%edge(edge)%flux_acc(p))
                         if (.not.mesh%edge(edge)%flux_acc(p)) then
                            if (cell1>0) then
                               soltemp%val(cell1,:)=soltemp%val(cell1,:)+ &
                                    gauss_weight(p)*mesh%edge(edge)%flux(p,:)*dt/(mesh%edge(edge)%lengthN*2.0_dp)
                            endif
                            if (cell2>0) then
                               soltemp%val(cell2,:)=soltemp%val(cell2,:)- &
                                    gauss_weight(p)*mesh%edge(edge)%flux(p,:)*dt/(mesh%edge(edge)%lengthN*2.0_dp)
                            endif
                         endif
                      enddo
                   end select
                   
                endif
             enddo
          endif
       enddo

    enddo

    deallocate(NOT_ACCEPTED_CELL,NOT_ACCEPTED_EDGE)
    allocate(NOT_ACCEPTED_CELL(nc),NOT_ACCEPTED_EDGE(ne))

    NOT_ACCEPTED_CELL=NAC(1:nc)
    NOT_ACCEPTED_EDGE=NAE(1:ne)

    deallocate(NAC,NAE)
    
    return
  end subroutine decrement

  subroutine DMP(mesh,sol,sol2,k,isol,eps,gauss_weight,str_equa,accept)
    type(meshStruct), intent(inout) :: mesh
    type(solStruct), intent(in) :: sol,sol2
    integer, intent(in) :: k,isol
    real(dp), intent(in) :: eps
    real(dp), dimension(:), intent(in) :: gauss_weight
    character(len=20), intent(in) :: str_equa
    logical, intent(out) :: accept
    integer :: j,neigh
    real(dp) :: mini,maxi,test
    if(.false.)print*,gauss_weight,str_equa

    call norme2(sol,k,isol,mini)
    call norme2(sol,k,isol,maxi)
    do j=1,size(mesh%cell(k)%neigh)
       neigh=abs(mesh%cell(k)%neigh(j))
       call norme2(sol,neigh,isol,test)
       if (test<mini) then
          mini=test
       else if (test>maxi) then
          maxi=test
       endif
    enddo

    call norme2(sol2,k,isol,test)
    if ((test-mini>=eps).and.(test-maxi<=eps)) then
       accept=.true.
    else
       accept=.false.
    endif

    return
  end subroutine DMP

  subroutine DMPu2(mesh,sol,sol2,k,isol,eps,gauss_weight,str_equa,accept)
    type(meshStruct), intent(inout) :: mesh
    type(solStruct), intent(in) :: sol,sol2
    integer, intent(in) :: k,isol
    real(dp), intent(in) :: eps
    real(dp), dimension(:), intent(in) :: gauss_weight
    character(len=20), intent(in) :: str_equa
    logical, intent(out) :: accept
    integer :: i,j,neigh
    real(dp) :: mini,maxi,test
    logical :: extrema
    if(.false.)print*,str_equa

    call norme2(sol,k,isol,mini)
    call norme2(sol,k,isol,maxi)
    do j=1,size(mesh%cell(k)%neigh)
       neigh=abs(mesh%cell(k)%neigh(j))
       call norme2(sol,neigh,isol,test)
       if (test<mini) then
          mini=test
       else if (test>maxi) then
          maxi=test
       endif
    enddo

    call norme2(sol2,k,isol,test)
    if ((test-mini>=-eps).and.(test-maxi<=eps)) then
       accept=.true.
    else
       accept=.false.
       do i=sol%conserv_var(isol,1),sol%conserv_var(isol,2)
          call u2(mesh,sol2,k,i,gauss_weight,extrema)
          if (extrema) then
             accept=.true.
          endif
       enddo
    endif

    return
  end subroutine DMPu2

  subroutine u2(mesh,sol,k,isol,gauss_weight,extrema)
    type(meshStruct), intent(inout) :: mesh
    type(solStruct), intent(in) :: sol
    integer, intent(in) :: k,isol
    real(dp), dimension(:), intent(in) :: gauss_weight
    logical, intent(out) :: extrema
    integer :: j,neigh
    real(dp) :: Xmin,Xmax,Ymin,Ymax

    call reconstruct(mesh,sol,k,3,gauss_weight)
    do j=1,size(mesh%cell(k)%neigh)
       neigh=mesh%cell(k)%neigh(j)
       if (neigh>0) then
          call reconstruct(mesh,sol,neigh,3,gauss_weight)
       endif
    enddo

    Xmin=mesh%cell(k)%polCoef(5,isol)
    Xmax=mesh%cell(k)%polCoef(5,isol)
    Ymin=mesh%cell(k)%polCoef(3,isol)
    Ymax=mesh%cell(k)%polCoef(3,isol)
    do j=1,size(mesh%cell(k)%neigh)
       neigh=mesh%cell(k)%neigh(j)
       if (neigh>0) then
          if (mesh%cell(neigh)%polCoef(5,isol)<Xmin) then
             Xmin=mesh%cell(neigh)%polCoef(5,isol)
          else if (mesh%cell(neigh)%polCoef(5,isol)>Xmax) then
             Xmax=mesh%cell(neigh)%polCoef(5,isol)
          endif
          if (mesh%cell(neigh)%polCoef(3,isol)<Ymin) then
             Ymin=mesh%cell(neigh)%polCoef(3,isol)
          else if (mesh%cell(neigh)%polCoef(3,isol)>Ymax) then
             Ymax=mesh%cell(neigh)%polCoef(3,isol)
          endif
       endif
    enddo

    if (Xmax==0.0_dp) then
       extrema=.true.
    else
       if ((Xmin/Xmax>0.5_dp).and.(Ymin/Ymax>0.5_dp)) then
          extrema=.true.
       else
          extrema=.false.
       endif
    endif

    deallocate(mesh%cell(k)%polCoef)
    do j=1,size(mesh%cell(k)%neigh)
       neigh=mesh%cell(k)%neigh(j)
       if (neigh>0) then
          deallocate(mesh%cell(neigh)%polCoef)
       endif
    enddo

    return
  end subroutine u2

  subroutine PAD(mesh,sol,sol2,k,isol,eps,gauss_weight,str_equa,accept)
    type(meshStruct), intent(inout) :: mesh
    type(solStruct), intent(in) :: sol,sol2
    integer, intent(in) :: k,isol
    real(dp), intent(in) :: eps
    real(dp), dimension(:), intent(in) :: gauss_weight
    character(len=20), intent(in) :: str_equa
    logical, intent(out) :: accept
    real(dp) :: rho,p
    if(.false.)print*,eps,gauss_weight,isol,mesh%nc,sol%nsolUser

    call unconserv(sol2%val(k,:),str_equa,1,rho)
    call unconserv(sol2%val(k,:),str_equa,4,p)
    if (rho>=0.0_dp.and.p>=0.0_dp) then
       accept=.true.
    else
       accept=.false.
    endif

    return
  end subroutine PAD

  subroutine criteria_flux(flux,accept)
    real(dp), dimension(:), intent(in) :: flux
    logical, intent(out) :: accept
    if(.false.)print*,flux
    
    accept=.false.
    
    return
  end subroutine criteria_flux

  subroutine norme2(sol,k,isol,norm)
    type(solStruct), intent(in) :: sol
    integer, intent(in) :: k,isol
    real(dp), intent(out) :: norm
    integer :: i

    select case (sol%conserv_var(isol,2)-sol%conserv_var(isol,1))
    case(0)
       norm=sol%val(k,isol)
    case default
       norm=0.0_dp
       do i=sol%conserv_var(isol,1),sol%conserv_var(isol,2)
          norm=norm+sol%val(k,i)**2
       enddo
       norm=sqrt(norm)
    end select

    return
  end subroutine norme2
    
end module limit
