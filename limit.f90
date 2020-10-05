module limit

  use constant
  use types
  use reconstruction

  implicit none

contains

  subroutine detect(mesh,sol,sol2,str_equa,L_str_criteria,L_var_criteria,L_eps,gauss_weight,period,NOT_ACCEPTED_CELL)
     type(meshStruct), intent(inout) :: mesh
     type(solStruct), intent(in) :: sol
     type(solStruct), intent(in) :: sol2
     character(len=20), intent(in) :: str_equa
     integer, dimension(:), intent(in) :: L_var_criteria
     character(len=20), dimension(:), intent(in) :: L_str_criteria
     real(dp), dimension(:), intent(in) :: L_eps
     real(dp), dimension(:), intent(in) :: gauss_weight
     logical, intent(in) :: period
     integer, dimension(:), intent(inout) :: NOT_ACCEPTED_CELL
     integer :: n,isol,i,k
     procedure (sub_criteria), pointer :: criteria
     logical :: accept

     do n=1,size(L_str_criteria)
        isol=L_var_criteria(n)
        select case (trim(L_str_criteria(n)))
        case ('DMP')
           criteria => DMP
        case ('DMPu2')
           criteria => DMPu2
        case('PAD')
           select case (trim(str_equa))
           case('euler')
              criteria => PAD_euler
           case('euler_is')
              criteria => PAD_euler
           case('M1')
              criteria => PAD_M1
           case default
              print*,trim(L_str_criteria(n))," is not valid for equation ",trim(str_equa)
              call exit()
           end select
        case default
           print*,trim(L_str_criteria(n))," criteria not implemented"
           call exit()
        end select

        do i=1,size(NOT_ACCEPTED_CELL)
          k=NOT_ACCEPTED_CELL(i)

          if (mesh%cell(k)%deg==0) then
             accept=.true.
          else
             call criteria(mesh,sol,sol2,k,isol,L_eps(n),gauss_weight,period,str_equa,accept)
          endif

          if (.not.accept) mesh%cell(k)%accept=.false.
       enddo

     enddo

  end subroutine detect

  subroutine decrement(mesh,sol,soltemp,str_equa,deg,nrk,dt,L_str_criteria,L_var_criteria,L_eps, &
       gauss_weight,period,NOT_ACCEPTED_CELL,NOT_ACCEPTED_EDGE,NAC_reason,verbosity)
    type(meshStruct), intent(inout) :: mesh
    type(solStruct), intent(in) :: sol
    type(solStruct), intent(inout) :: soltemp
    character(len=20), intent(in) :: str_equa
    integer, intent(in) :: deg,nrk,verbosity
    real(dp), intent(in) :: dt
    integer, dimension(:), intent(in) :: L_var_criteria
    character(len=20), dimension(:), intent(in) :: L_str_criteria
    real(dp), dimension(:), intent(in) :: L_eps
    real(dp), dimension(:), intent(in) :: gauss_weight
    logical, intent(in) :: period
    type(solStruct) :: sol2
    integer, dimension(:), allocatable, intent(inout) :: NOT_ACCEPTED_CELL,NOT_ACCEPTED_EDGE,NAC_reason
    integer, dimension(:), allocatable :: NAC,NAE
    integer :: i,j,k,n,p,isol,nc,ne,cell1,cell2,edge,sub1,sub2,neigh
    procedure (sub_criteria), pointer :: criteria
    logical :: accept
    real(dp) :: lengthN1,lengthN2
    if(.false.)print*,dt,neigh,p

    allocate(NAC(mesh%nc),NAE(mesh%ne))

    nc=0
    ne=0
    NAC=0
    NAE=0
    sol2=soltemp
    if (verbosity>1.and.size(NOT_ACCEPTED_CELL)==mesh%nc) NAC_reason=0
    
    do n=1,size(L_str_criteria)
       
       isol=L_var_criteria(n)
       
       select case (trim(L_str_criteria(n)))
       case ('DMP')
          criteria => DMP
       case ('DMPu2')
          criteria => DMPu2
       case('PAD')
          select case (trim(str_equa))
          case('euler')
             criteria => PAD_euler
          case('euler_is')
             criteria => PAD_euler
          case('M1')
             criteria => PAD_M1
          case default
             print*,trim(L_str_criteria(n))," is not valid for equation ",trim(str_equa)
             call exit()
          end select
       case default
          print*,trim(L_str_criteria(n))," criteria not implemented"
          call exit()
       end select

       do i=1,size(NOT_ACCEPTED_CELL)
          k=NOT_ACCEPTED_CELL(i)

          if (mesh%cell(k)%deg==0) then
             accept=.true.
          else if (mesh%cell(k)%accept) then
             call criteria(mesh,sol,sol2,k,isol,L_eps(n),gauss_weight,period,str_equa,accept)
          else
             accept=.false.
          endif

          if (.not.accept) then
             if (verbosity>1) then
                if (size(NOT_ACCEPTED_CELL)==mesh%nc.and.NAC_reason(k)==0) NAC_reason(k)=n
             endif
             mesh%cell(k)%deg=deg
             mesh%cell(k)%accept=.false.
             if (all(NAC/=k)) then
                nc=nc+1
                NAC(nc)=k
             endif
             call crown(mesh,k,nc,1,nrk,deg,NAC)
          endif
       enddo

       do i=1,nc
          k=NAC(i)
          do j=1,size(mesh%cell(k)%edge)
             edge=mesh%cell(k)%edge(j)
             cell1=mesh%edge(edge)%cell1
             cell2=mesh%edge(edge)%cell2
             lengthN1=mesh%cell(abs(cell1))%dx
             lengthN2=mesh%cell(abs(cell2))%dx
             sub1=mesh%edge(edge)%sub(1)
             sub2=mesh%edge(edge)%sub(2)
             
             if (all(NAE/=edge)) then

                select case (trim(mesh%edge(edge)%boundtype))

                case('PERIODIC')
                   ne=ne+1
                   NAE(ne)=edge
                   if (cell1*cell2<0) then
                      ne=ne+1
                      NAE(ne)=mesh%edge(edge)%period
                   endif
                   !call criteria_flux(mesh%edge(edge)%flux(1,:),mesh%edge(edge)%flux_acc(1))

                case default
                   ne=ne+1
                   NAE(ne)=edge
                   !call criteria_flux(mesh%edge(edge)%flux(1,:),mesh%edge(edge)%flux_acc(1))

                end select

                !do p=1,size(gauss_weight)
                   !call criteria_flux(mesh%edge(edge)%flux(p,:),mesh%edge(edge)%flux_acc(p))
                !enddo
                
             endif
          enddo
       enddo
    enddo

    deallocate(NOT_ACCEPTED_CELL,NOT_ACCEPTED_EDGE)
    allocate(NOT_ACCEPTED_CELL(nc),NOT_ACCEPTED_EDGE(ne))

    NOT_ACCEPTED_CELL=NAC(1:nc)
    NOT_ACCEPTED_EDGE=NAE(1:ne)

    deallocate(NAC,NAE)

    return
  end subroutine decrement

  recursive subroutine crown(mesh,k,nc,i,nrk,deg,NAC)
    type(meshStruct), intent(inout) :: mesh
    integer, intent(in) :: k,i,nrk,deg
    integer, intent(inout) :: nc
    integer, dimension(:), intent(inout) :: NAC
    integer :: j,neigh

    do j=1,size(mesh%cell(k)%stencil)
       neigh=mesh%cell(k)%stencil(j)
       if (i<=nrk) then
          mesh%cell(neigh)%accept=.false.
          !mesh%cell(neigh)%deg=deg
       endif
       if (all(NAC/=neigh)) then
          nc=nc+1
          NAC(nc)=neigh
       endif
    enddo
    if (i<nrk) then
       do j=1,size(mesh%cell(k)%stencil)
          neigh=mesh%cell(k)%stencil(j)
          call crown(mesh,neigh,nc,i+1,nrk,deg,NAC)
       enddo
    endif
  
    return
  end subroutine crown

  subroutine DMP(mesh,sol,sol2,k,isol,eps,gauss_weight,period,str_equa,accept)
    type(meshStruct), intent(inout) :: mesh
    type(solStruct), intent(in) :: sol,sol2
    integer, intent(in) :: k,isol
    real(dp), intent(in) :: eps
    real(dp), dimension(:), intent(in) :: gauss_weight
    logical, intent(in) :: period
    character(len=20), intent(in) :: str_equa
    logical, intent(inout) :: accept
    integer :: j,neigh
    real(dp) :: mini,maxi,test,eps2
    if(.false.)print*,gauss_weight,str_equa,period

    call norme2(sol,k,isol,mini)
    call norme2(sol,k,isol,maxi)
    do j=1,size(mesh%cell(k)%stencil)     !neigh ?
       neigh=abs(mesh%cell(k)%stencil(j))
       call norme2(sol,neigh,isol,test)
       if (test<mini) then
          mini=test
       else if (test>maxi) then
          maxi=test
       endif
    enddo

    eps2=max(eps,eps*(maxi-mini))
    call norme2(sol2,k,isol,test)
    if ((test-mini>=-eps2).and.(test-maxi<=eps2)) then
       accept=.true.
    else
       accept=.false.
    endif

    return
  end subroutine DMP

  subroutine DMPu2(mesh,sol,sol2,k,isol,eps,gauss_weight,period,str_equa,accept)
    type(meshStruct), intent(inout) :: mesh
    type(solStruct), intent(in) :: sol,sol2
    integer, intent(in) :: k,isol
    real(dp), intent(in) :: eps
    real(dp), dimension(:), intent(in) :: gauss_weight
    logical, intent(in) :: period
    character(len=20), intent(in) :: str_equa
    logical, intent(inout) :: accept
    integer :: i,j,neigh
    real(dp) :: mini,maxi,test,eps2
    logical :: extrema
    if(.false.)print*,str_equa

    call norme2(sol,k,isol,mini)
    call norme2(sol,k,isol,maxi)
    do j=1,size(mesh%cell(k)%stencil)     !neigh ?
       neigh=abs(mesh%cell(k)%stencil(j))
       call norme2(sol,neigh,isol,test)
       if (test<mini) then
          mini=test
       else if (test>maxi) then
          maxi=test
       endif
    enddo

    eps2=max(eps,eps*(maxi-mini))
    call norme2(sol2,k,isol,test)
    if ((test-mini>=-eps2).and.(test-maxi<=eps2)) then
       accept=.true.
    else
       accept=.false.
       do i=sol%conserv_var(isol,1),sol%conserv_var(isol,2)
          call u2(mesh,sol2,k,i,gauss_weight,period,extrema)
          if (extrema) then
             accept=.true.
          endif
       enddo
    endif
       
    return
  end subroutine DMPu2

  subroutine u2(mesh,sol,k,isol,gauss_weight,period,extrema)
    type(meshStruct), intent(inout) :: mesh
    type(solStruct), intent(in) :: sol
    integer, intent(in) :: k,isol
    real(dp), dimension(:), intent(in) :: gauss_weight
    logical, intent(in) :: period
    logical, intent(out) :: extrema
    integer :: j,neigh,test_min,test_max
    real(dp) :: Xmin,Xmax,Ymin,Ymax

    call reconstruct1(mesh,sol,k,3,gauss_weight,period,mesh%cell(k)%polCoef3)
    do j=1,size(mesh%cell(k)%edge)
       neigh=mesh%cell(k)%neigh(j)
       if (neigh>0) then
          call reconstruct1(mesh,sol,neigh,3,gauss_weight,period,mesh%cell(neigh)%polCoef3)
       endif
    enddo

    Xmin=mesh%cell(k)%polCoef3(5,isol)
    Xmax=mesh%cell(k)%polCoef3(5,isol)
    Ymin=mesh%cell(k)%polCoef3(3,isol)
    Ymax=mesh%cell(k)%polCoef3(3,isol)
    test_min=k
    test_max=k
    do j=1,size(mesh%cell(k)%edge)
       neigh=mesh%cell(k)%neigh(j)
       if (neigh>0) then
          if (abs(mesh%cell(neigh)%polCoef3(5,isol))<abs(Xmin)) then
             Xmin=mesh%cell(neigh)%polCoef3(5,isol)
             test_min=neigh
          else if (abs(mesh%cell(neigh)%polCoef3(5,isol))>abs(Xmax)) then
             Xmax=mesh%cell(neigh)%polCoef3(5,isol)
             test_max=neigh
          endif
          if (abs(mesh%cell(neigh)%polCoef3(3,isol))<abs(Ymin)) then
             Ymin=mesh%cell(neigh)%polCoef3(3,isol)
          else if (abs(mesh%cell(neigh)%polCoef3(3,isol))>abs(Ymax)) then
             Ymax=mesh%cell(neigh)%polCoef3(3,isol)
          endif
       endif
    enddo

    if (Xmax==0.0_dp.or.Ymax==0.0_dp) then
       extrema=.true.
    else
       if ((Xmin/Xmax>0.5_dp).and.(Ymin/Ymax>0.5_dp)) then
          extrema=.true.
       else
          extrema=.false.
       endif
    endif

    deallocate(mesh%cell(k)%polCoef3)
    do j=1,size(mesh%cell(k)%edge)
       neigh=mesh%cell(k)%neigh(j)
       if (neigh/=k.and.neigh>0) then
          deallocate(mesh%cell(neigh)%polCoef3)
       endif
    enddo

    return
  end subroutine u2

  subroutine PAD_euler(mesh,sol,sol2,k,isol,eps,gauss_weight,period,str_equa,accept)
    type(meshStruct), intent(inout) :: mesh
    type(solStruct), intent(in) :: sol,sol2
    integer, intent(in) :: k,isol
    real(dp), intent(in) :: eps
    real(dp), dimension(:), intent(in) :: gauss_weight
    logical, intent(in) :: period
    character(len=20), intent(in) :: str_equa
    logical, intent(inout) :: accept
    real(dp) :: rho,p
    if(.false.)print*,eps,gauss_weight,isol,mesh%nc,sol%nsolUser,period

    call unconserv(sol2%val(k,:),str_equa,1,rho)
    accept=.false.
    if (rho>eps) then
       call unconserv(sol2%val(k,:),str_equa,4,p)
       if (p>eps) then
          accept=.true.
       endif
    endif

    return
  end subroutine PAD_euler

  subroutine PAD_M1(mesh,sol,sol2,k,isol,eps,gauss_weight,period,str_equa,accept)
    type(meshStruct), intent(inout) :: mesh
    type(solStruct), intent(in) :: sol,sol2
    integer, intent(in) :: k,isol
    real(dp), intent(in) :: eps
    real(dp), dimension(:), intent(in) :: gauss_weight
    logical, intent(in) :: period
    character(len=20), intent(in) :: str_equa
    logical, intent(inout) :: accept
    real(dp) :: test
    if(.false.)print*,gauss_weight,isol,mesh%nc,sol%nsolUser,period
    if(.false.)call unconserv(sol2%val(k,:),str_equa,1,test)

    accept=.false.
    if (sol2%val(k,1)<eps) then
       accept=.true.
    endif
    if ((sol2%val(k,2)<eps.or.sol2%val(k,3)<eps).or.(sol2%val(k,2)>1.0_dp-eps.or.sol2%val(k,3)>1.0_dp-eps)) then
       accept=.true.
    endif

    return
  end subroutine PAD_M1

  subroutine criteria_flux(flux,accept)
    real(dp), dimension(:), intent(in) :: flux
    logical, intent(inout) :: accept
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
