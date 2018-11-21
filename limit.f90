module limit

  use constant
  use types
  use reconstruction

  implicit none

contains

  subroutine decrement(mesh,sol,soltemp,deg,dt,L_str_criteria,L_var_criteria,gauss_weight,NOT_ACCEPTED_CELL,NOT_ACCEPTED_EDGE)
    type(meshStruct), intent(inout) :: mesh
    type(solStruct), intent(in) :: sol
    type(solStruct), intent(inout) :: soltemp
    integer, intent(in) :: deg
    real(dp), intent(in) :: dt
    integer, dimension(:), intent(in) :: L_var_criteria
    character(len=20), dimension(:), intent(in) :: L_str_criteria
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
       case default
          print*,trim(L_str_criteria(n))," criteria not implemented"
          call exit()
       end select

       do i=1,size(NOT_ACCEPTED_CELL)
          k=NOT_ACCEPTED_CELL(i)
          call criteria(mesh,sol,sol2,k,isol,gauss_weight,accept)
          if (.not.accept) then
             do j=1,size(mesh%cell(k)%edge)
                edge=mesh%cell(k)%edge(j)
                cell1=mesh%edge(edge)%cell1
                cell2=mesh%edge(edge)%cell2
                if (all(NAE/=edge)) then
                   ne=ne+1
                   NAE(ne)=edge
                   if (cell1*cell2<0) then
                      ne=ne+1
                      NAE(ne)=mesh%edge(edge)%period
                   endif
                   do p=1,size(gauss_weight)
                      call criteria_flux(mesh%edge(edge)%flux(p,:),mesh%edge(edge)%flux_acc(p))
                      if (.not.mesh%edge(edge)%flux_acc(p)) then
                         soltemp%val(abs(cell1),:)=soltemp%val(abs(cell1),:)+ &
                              gauss_weight(p)*mesh%edge(edge)%flux(p,:)*dt/(mesh%edge(edge)%length*2.0_dp)
                         soltemp%val(abs(cell2),:)=soltemp%val(abs(cell2),:)- &
                              gauss_weight(p)*mesh%edge(edge)%flux(p,:)*dt/(mesh%edge(edge)%length*2.0_dp)
                      endif
                   enddo
                endif
                cell1=abs(cell1)
                cell2=abs(cell2)
                if (all(NAC/=cell1)) then
                   mesh%cell(cell1)%deg=deg
                   nc=nc+1
                   NAC(nc)=cell1
                endif
                if (all(NAC/=cell2)) then
                   mesh%cell(cell2)%deg=deg
                   nc=nc+1
                   NAC(nc)=cell2
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

  subroutine DMP(mesh,sol,sol2,k,isol,gauss_weight,accept)
    type(meshStruct), intent(inout) :: mesh
    type(solStruct), intent(in) :: sol,sol2
    integer, intent(in) :: k,isol
    real(dp), dimension(:), intent(in) :: gauss_weight
    logical, intent(out) :: accept
    integer :: j,neigh
    real(dp) :: mini,maxi,eps2

    mini=sol%val(k,isol)
    maxi=sol%val(k,isol)
    do j=1,size(mesh%cell(k)%neigh)
       neigh=abs(mesh%cell(k)%neigh(j))
       if (sol%val(neigh,isol)<mini) then
          mini=sol%val(neigh,isol)
       else if (sol%val(neigh,isol)>maxi) then
          maxi=sol%val(neigh,isol)
       endif
    enddo

    eps2=max((maxi-mini)*1.0e-6_dp,eps)
    if ((sol2%val(k,isol)-mini>=-eps2).and.(sol2%val(k,isol)-maxi<=eps2)) then
       accept=.true.
    else
       accept=.false.
    endif

    return
  end subroutine DMP

  subroutine DMPu2(mesh,sol,sol2,k,isol,gauss_weight,accept)
    type(meshStruct), intent(inout) :: mesh
    type(solStruct), intent(in) :: sol,sol2
    integer, intent(in) :: k,isol
    real(dp), dimension(:), intent(in) :: gauss_weight
    logical, intent(out) :: accept
    integer :: j,neigh
    real(dp) :: mini,maxi,eps2
    logical :: extrema

    mini=sol%val(k,isol)
    maxi=sol%val(k,isol)
    do j=1,size(mesh%cell(k)%neigh)
       neigh=abs(mesh%cell(k)%neigh(j))
       if (sol%val(neigh,isol)<mini) then
          mini=sol%val(neigh,isol)
       else if (sol%val(neigh,isol)>maxi) then
          maxi=sol%val(neigh,isol)
       endif
    enddo

    eps2=max((maxi-mini)*1.0e-6_dp,eps)
    if ((sol2%val(k,isol)-mini>=-eps2).and.(sol2%val(k,isol)-maxi<=eps2)) then
       accept=.true.
    else
       call u2(mesh,sol2,k,isol,gauss_weight,extrema)
       if (extrema) then
          accept=.true.
       else
          accept=.false.
       endif
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

  subroutine criteria_flux(flux,accept)
    real(dp), dimension(:), intent(in) :: flux
    logical, intent(out) :: accept

    accept=.false.
    
    return
  end subroutine criteria_flux
    
end module limit
