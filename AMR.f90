module AMR
  
  use constant
  use types
  use reconstruction
  use ISO_C_BINDING

contains

  subroutine adapt(fn_adapt,p4est,quadrants,mesh,sol,level,order,gauss_weight,gauss_point)
    procedure (sub_adapt), pointer, intent(in) :: fn_adapt
    type(c_ptr), intent(in) :: p4est,quadrants
    type(meshStruct), intent(inout) :: mesh
    type(solStruct), intent(in) :: sol
    integer, intent(in) :: level,order
    real(dp), dimension(:), intent(in) :: gauss_weight,gauss_point
    type(c_ptr) :: C_sol,C_sol_coarsen,C_sol_refine,C_sol_interp
    real(dp), dimension(:,:), allocatable, target :: F_sol,F_sol_interp
    integer, dimension(:), allocatable, target :: sol_coarsen,sol_refine
    integer :: maxlevel,coarsen_recursive,refine_recursive
    
    allocate(F_sol(mesh%nc,sol%nvar),F_sol_interp(4*mesh%nc,sol%nvar))
    allocate(sol_coarsen(mesh%nc),sol_refine(mesh%nc))
    
    call fn_adapt(mesh,sol,level,maxlevel,coarsen_recursive,refine_recursive, &
         sol_coarsen,sol_refine)
    call interp(mesh,sol,order,gauss_weight,gauss_point,sol_refine,F_sol_interp)
    C_sol_interp=c_loc(F_sol_interp)
    
    F_sol=sol%val
    C_sol=c_loc(F_sol)
    C_sol_coarsen=c_loc(sol_coarsen)
    C_sol_refine=c_loc(sol_refine)
    
    call p4_adapt(p4est,quadrants,C_sol,C_sol_interp,sol%nvar,C_sol_coarsen,C_sol_refine,maxlevel, &
         coarsen_recursive,refine_recursive)
    
    deallocate(F_sol,F_sol_interp,sol_coarsen,sol_refine)

    return
  end subroutine adapt

  subroutine interp(mesh,sol,order,gauss_weight,gauss_point,sol_refine,sol_interp)
    type(meshStruct), intent(inout) :: mesh
    type(solStruct), intent(in) :: sol
    integer, intent(in) :: order
    real(dp), dimension(:), intent(in) :: gauss_weight,gauss_point
    integer, dimension(:), intent(in) :: sol_refine
    real(dp), dimension(:,:), intent(inout):: sol_interp
    integer :: k,p1,p2
    real(dp), dimension(:), allocatable :: u1,u2,u3,u4
    real(dp) :: xc,yc,dx,dy,Xg,Yg

    allocate(u1(sol%nvar),u2(sol%nvar),u3(sol%nvar),u4(sol%nvar))

    sol_interp=0.0_dp
    do k=1,mesh%nc
       select case (sol_refine(k))
       case(0)
          sol_interp(k,:)=sol%val(k,:)
          sol_interp(k+mesh%nc,:)=sol%val(k,:)
          sol_interp(k+2*mesh%nc,:)=sol%val(k,:)
          sol_interp(k+3*mesh%nc,:)=sol%val(k,:)
       case(1)
          if (allocated(mesh%cell(k)%polCoef)) deallocate(mesh%cell(k)%polCoef)
          allocate(mesh%cell(k)%polcoef(order*(order-1)/2+order-1,sol%nvar))
          call reconstruct(mesh,sol,k,order,gauss_weight)
          xc=mesh%cell(k)%xc
          yc=mesh%cell(k)%yc
          dx=mesh%cell(k)%dx/4.0_dp
          dy=mesh%cell(k)%dy/4.0_dp
          u1=sol%val(k,:)
          u2=sol%val(k,:)
          u3=sol%val(k,:)
          u4=sol%val(k,:)
          do p1=1,order
             do p2=1,order
                Xg=xc-dx+gauss_point(p1)*dx
                Yg=yc-dy+gauss_point(p2)*dy
                call evaluate(mesh,sol,k,order,gauss_weight,Xg,Yg,u1)
                sol_interp(k,:)=sol_interp(k,:)+u1*gauss_weight(p1)*gauss_weight(p2)/4.0_dp
                Xg=xc+dx+gauss_point(p1)*dx
                Yg=yc-dy+gauss_point(p2)*dy
                call evaluate(mesh,sol,k,order,gauss_weight,Xg,Yg,u2)
                sol_interp(k+mesh%nc,:)=sol_interp(k+mesh%nc,:)+u2*gauss_weight(p1)*gauss_weight(p2)/4.0_dp
                Xg=xc-dx+gauss_point(p1)*dx
                Yg=yc+dy+gauss_point(p2)*dy
                call evaluate(mesh,sol,k,order,gauss_weight,Xg,Yg,u3)
                sol_interp(k+2*mesh%nc,:)=sol_interp(k+2*mesh%nc,:)+u3*gauss_weight(p1)*gauss_weight(p2)/4.0_dp
                Xg=xc+dx+gauss_point(p1)*dx
                Yg=yc+dy+gauss_point(p2)*dy
                call evaluate(mesh,sol,k,order,gauss_weight,Xg,Yg,u4)
                sol_interp(k+3*mesh%nc,:)=sol_interp(k+3*mesh%nc,:)+u4*gauss_weight(p1)*gauss_weight(p2)/4.0_dp
             enddo
          enddo
       end select
    enddo

    deallocate(u1,u2,u3,u4)
          
    return
  end subroutine interp

  subroutine new_sol(mesh,quadrants,sol)
    type(meshStruct), intent(in) :: mesh
    type(c_ptr), intent(in) :: quadrants
    type(solStruct), intent(inout) :: sol
    type(c_ptr) :: C_sol
    real(dp), dimension(:), pointer :: F_sol
    integer :: k,isol

    call p4_new_sol(quadrants,C_sol)
    call c_f_pointer(C_sol,F_sol,(/sol%nvar*mesh%nc/))
    do k=1,mesh%nc
       do isol=1,sol%nvar
          sol%val(k,isol)=F_sol((isol-1)*mesh%nc+k)
       enddo
    enddo
    call p4_free(C_sol)

    return
  end subroutine new_sol

  subroutine adapt_zone(mesh,sol,level,maxlevel,coarsen_recursive,refine_recursive,sol_coarsen,sol_refine)
    type(meshStruct), intent(inout) :: mesh
    type(solStruct), intent(in) :: sol
    integer, intent(in) :: level
    integer, intent(out) :: maxlevel,coarsen_recursive,refine_recursive
    integer, dimension(:), intent(inout) :: sol_coarsen,sol_refine
    integer :: k,isol
    real(dp) :: xmin,xmax,ymin,ymax
    if(.false.)print*,sol%nvar

    maxlevel=level+1
    coarsen_recursive=0
    refine_recursive=0
    isol=1
    xmin=4.0_dp
    xmax=6.0_dp
    ymin=4.0_dp
    ymax=6.0_dp
    do k=1,mesh%nc
       sol_coarsen(k)=0
       sol_refine(k)=0
       if (mesh%cell(k)%xc>xmin.and.mesh%cell(k)%xc<xmax) then
          if (mesh%cell(k)%yc>ymin.and.mesh%cell(k)%yc<ymax) then
             sol_refine(k)=1
          endif
       endif
    enddo
    
    return
  end subroutine adapt_zone

  subroutine adapt_sinus(mesh,sol,level,maxlevel,coarsen_recursive,refine_recursive,sol_coarsen,sol_refine)
    type(meshStruct), intent(inout) :: mesh
    type(solStruct), intent(in) :: sol
    integer, intent(in) :: level
    integer, intent(out) :: maxlevel,coarsen_recursive,refine_recursive
    integer, dimension(:), intent(inout) :: sol_coarsen,sol_refine
    integer :: k,isol

    maxlevel=level+1
    coarsen_recursive=0
    refine_recursive=0
    isol=1
    do k=1,mesh%nc
       sol_coarsen(k)=0
       sol_refine(k)=0
       if (sol%val(k,isol)<0.5_dp.and.mesh%cell(k)%level>level) sol_coarsen(k)=1
       if (sol%val(k,isol)>0.5_dp) sol_refine(k)=1
    enddo
    
    return
  end subroutine adapt_sinus

  subroutine adapt_vortex(mesh,sol,level,maxlevel,coarsen_recursive,refine_recursive,sol_coarsen,sol_refine)
    type(meshStruct), intent(inout) :: mesh
    type(solStruct), intent(in) :: sol
    integer, intent(in) :: level
    integer, intent(out) :: maxlevel,coarsen_recursive,refine_recursive
    integer, dimension(:), intent(inout) :: sol_coarsen,sol_refine
    integer :: k,isol
    real(dp) :: grad

    maxlevel=level+1
    coarsen_recursive=0
    refine_recursive=0
    isol=1
    do k=1,mesh%nc
       sol_coarsen(k)=0
       sol_refine(k)=0
       call reconstruct(mesh,sol,k,2,gauss_weight2)
       grad=sqrt(mesh%cell(k)%polCoef(1,isol)**2+mesh%cell(k)%polCoef(2,isol)**2)
       if (grad<0.001_dp.and.mesh%cell(k)%level>level-2) sol_coarsen(k)=1
       if (grad>0.01_dp) sol_refine(k)=1
    enddo
    
    return
  end subroutine adapt_vortex
    

end module AMR
