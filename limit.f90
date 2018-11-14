module limit

  use constant
  use types
  use reconstruction

  implicit none

contains

  subroutine decrement(mesh,sol,soltemp,dt,L_str_criteria,L_var_criteria,NOT_ACCEPTED_CELL,NOT_ACCEPTED_EDGE)
    type(meshStruct), intent(inout) :: mesh
    type(solStruct), intent(in) :: sol
    type(solStruct), intent(inout) :: soltemp
    real(dp), intent(in) :: dt
    integer, dimension(:), intent(in) :: L_var_criteria
    character(len=20), dimension(:), intent(in) :: L_str_criteria
    integer, dimension(:), allocatable, intent(inout) :: NOT_ACCEPTED_CELL,NOT_ACCEPTED_EDGE
    integer, dimension(:), allocatable :: NAC,NAE
    integer :: i,j,k,n,isol,nc,ne,cell1,cell2,edge
    procedure (sub_criteria), pointer :: criteria

    allocate(NAC(mesh%nc),NAE(mesh%ne))

    nc=0
    ne=0
    NAC=0
    NAE=0
    
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
          call criteria(mesh,sol,soltemp,k,isol)
          if (.not.mesh%cell(k)%accept) then
             do j=1,size(mesh%cell(k)%edge)
                edge=mesh%cell(k)%edge(j)
                cell1=mesh%edge(edge)%cell1
                cell2=mesh%edge(edge)%cell2
                if (all(NAE/=edge)) then
                   ne=ne+1
                   NAE(ne)=edge
                   if (cell1>0) then
                      soltemp%val(cell1,:)=soltemp%val(cell1,:)+mesh%edge(edge)%flux(:)*dt/mesh%edge(edge)%length
                   endif
                   if (cell2>0) then
                      soltemp%val(cell2,:)=soltemp%val(cell2,:)-mesh%edge(edge)%flux(:)*dt/mesh%edge(edge)%length
                   endif
                endif
                cell1=abs(cell1)
                cell2=abs(cell2)
                if (all(NAC/=cell1)) then
                   mesh%cell(cell1)%deg=0
                   nc=nc+1
                   NAC(nc)=cell1
                endif
                if (all(NAC/=cell2)) then
                   mesh%cell(cell2)%deg=0
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

  subroutine DMP(mesh,sol,sol2,k,isol)
    type(meshStruct), intent(inout) :: mesh
    type(solStruct), intent(in) :: sol,sol2
    integer, intent(in) :: k,isol
    integer :: j,neigh
    real(dp) :: min,max

    min=sol%val(k,isol)
    max=sol%val(k,isol)
    do j=1,size(mesh%cell(k)%neigh)
       neigh=abs(mesh%cell(k)%neigh(j))
       if (sol%val(neigh,isol)<min) then
          min=sol%val(neigh,isol)
       else if (sol%val(neigh,isol)>max) then
          max=sol%val(neigh,isol)
       endif
    enddo

    if ((sol2%val(k,isol)-min>=-eps).and.(sol2%val(k,isol)-max<=eps)) then
       mesh%cell(k)%accept=.true.
    endif

    return
  end subroutine DMP

  subroutine DMPu2(mesh,sol,sol2,k,isol)
    type(meshStruct), intent(inout) :: mesh
    type(solStruct), intent(in) :: sol,sol2
    integer, intent(in) :: k,isol
    integer :: j,neigh
    real(dp) :: min,max
    logical :: extrema

    min=sol%val(k,isol)
    max=sol%val(k,isol)
    do j=1,size(mesh%cell(k)%neigh)
       neigh=abs(mesh%cell(k)%neigh(j))
       if (sol%val(neigh,isol)<min) then
          min=sol%val(neigh,isol)
       else if (sol%val(neigh,isol)>max) then
          max=sol%val(neigh,isol)
       endif
    enddo

    if ((sol2%val(k,isol)-min>=-eps).and.(sol2%val(k,isol)-max<=eps)) then
       mesh%cell(k)%accept=.true.
    else
       call u2(mesh,sol2,k,isol,extrema)
       if (extrema) then
          mesh%cell(k)%accept=.true.
       endif
    endif

    return
  end subroutine DMPu2

  subroutine u2(mesh,sol,k,isol,extrema)
    type(meshStruct), intent(inout) :: mesh
    type(solStruct), intent(in) :: sol
    integer, intent(in) :: k,isol
    integer :: j,neigh
    real(dp) :: Xmin,Xmax,Ymin,Ymax
    logical :: extrema
    procedure (quadrature_c_alpha), pointer :: quad_c_alpha

    quad_c_alpha => quadrature3_c_alpha

    call reconstruct(mesh,sol,k,1,3,quad_c_alpha)
    do j=1,size(mesh%cell(k)%neigh)
       neigh=mesh%cell(k)%neigh(j)
       if (neigh>0) then
          call reconstruct(mesh,sol,neigh,1,3,quad_c_alpha)
       endif
    enddo

    Xmin=mesh%cell(k)%polCoefL(5,isol)
    Xmax=mesh%cell(k)%polCoefL(5,isol)
    Ymin=mesh%cell(k)%polCoefL(3,isol)
    Ymax=mesh%cell(k)%polCoefL(3,isol)
    do j=1,size(mesh%cell(k)%neigh)
       neigh=mesh%cell(k)%neigh(j)
       if (neigh>0) then
          if (mesh%cell(neigh)%polCoefL(5,isol)<Xmin) then
             Xmin=mesh%cell(neigh)%polCoefL(5,isol)
          else if (mesh%cell(neigh)%polCoefL(5,isol)>Xmax) then
             Xmax=mesh%cell(neigh)%polCoefL(5,isol)
          endif
          if (mesh%cell(neigh)%polCoefL(3,isol)<Ymin) then
             Ymin=mesh%cell(neigh)%polCoefL(3,isol)
          else if (mesh%cell(neigh)%polCoefL(3,isol)>Ymax) then
             Ymax=mesh%cell(neigh)%polCoefL(3,isol)
          endif
       endif
    enddo

    if ((Xmin/Xmax>0.5_dp).and.(Ymin/Ymax>0.5_dp)) then
       extrema=.true.
    else
       extrema=.false.
    endif

    deallocate(mesh%cell(k)%polCoefL)
    do j=1,size(mesh%cell(k)%neigh)
       neigh=mesh%cell(k)%neigh(j)
       if (neigh>0) then
          deallocate(mesh%cell(neigh)%polCoefL)
       endif
    enddo

    return
  end subroutine u2
    
end module limit
