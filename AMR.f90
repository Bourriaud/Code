module AMR
  
  use constant
  use types
  use reconstruction
  use ISO_C_BINDING

contains

  subroutine adapt(fn_adapt,p4est,quadrants,mesh,sol,level,order,gauss_weight,gauss_point,period,minlevel,maxlevel)
    procedure (sub_adapt), pointer, intent(in) :: fn_adapt
    type(c_ptr), intent(in) :: p4est,quadrants
    type(meshStruct), intent(inout) :: mesh
    type(solStruct), intent(in) :: sol
    integer, intent(in) :: level,order
    real(dp), dimension(:), intent(in) :: gauss_weight,gauss_point
    logical, dimension(2), intent(in) :: period
    integer, intent(out) :: minlevel,maxlevel
    type(c_ptr) :: C_sol,C_sol_coarsen,C_sol_refine,C_sol_interp,C_pol_interp
    real(dp), dimension(:,:), allocatable, target :: F_sol,F_sol_interp
    real(dp), dimension(:,:,:), allocatable, target :: F_pol_interp
    integer, dimension(:), allocatable, target :: sol_coarsen,sol_refine
    integer :: coarsen_recursive,refine_recursive,size_pol

    size_pol=order*(order-1)/2+order-1
    allocate(F_sol(mesh%nc,sol%nvar),F_sol_interp(4*mesh%nc,sol%nvar))
    allocate(F_pol_interp(mesh%nc,size_pol,sol%nvar))
    allocate(sol_coarsen(mesh%nc),sol_refine(mesh%nc))

    call fn_adapt(mesh,sol,level,minlevel,maxlevel,coarsen_recursive,refine_recursive, &
         period,sol_coarsen,sol_refine)
    call interp(mesh,sol,order,gauss_weight,gauss_point,period,sol_refine,F_sol_interp,F_pol_interp)
    C_sol_interp=c_loc(F_sol_interp)
    C_pol_interp=c_loc(F_pol_interp)
    
    F_sol=sol%val
    C_sol=c_loc(F_sol)
    C_sol_coarsen=c_loc(sol_coarsen)
    C_sol_refine=c_loc(sol_refine)

    call p4_adapt(p4est,quadrants,C_sol,C_sol_interp,C_pol_interp,size_pol,sol%nvar,C_sol_coarsen,C_sol_refine, &
         maxlevel,coarsen_recursive,refine_recursive)

    deallocate(F_sol,F_sol_interp,F_pol_interp,sol_coarsen,sol_refine)
    
    return
  end subroutine adapt

  subroutine interp(mesh,sol,order,gauss_weight,gauss_point,period,sol_refine,sol_interp,pol_interp)
    type(meshStruct), intent(inout) :: mesh
    type(solStruct), intent(in) :: sol
    integer, intent(in) :: order
    real(dp), dimension(:), intent(in) :: gauss_weight,gauss_point
    logical, dimension(2), intent(in) :: period
    integer, dimension(:), intent(in) :: sol_refine
    real(dp), dimension(:,:), intent(inout):: sol_interp
    real(dp), dimension(:,:,:), intent(inout):: pol_interp
    integer :: k,p1,p2,scheme
    real(dp), dimension(:), allocatable :: u1,u2,u3,u4
    real(dp) :: xc,yc,dx,dy,Xg,Yg

    allocate(u1(sol%nvar),u2(sol%nvar),u3(sol%nvar),u4(sol%nvar))

    scheme=2
    sol_interp=0.0_dp
    do k=1,mesh%nc
       pol_interp(k,1:size(mesh%cell(k)%polCoef(:,1)),1:sol%nvar)=mesh%cell(k)%polCoef
       select case (sol_refine(k))
       case(0)
          sol_interp(k,:)=sol%val(k,:)
          sol_interp(k+mesh%nc,:)=sol%val(k,:)
          sol_interp(k+2*mesh%nc,:)=sol%val(k,:)
          sol_interp(k+3*mesh%nc,:)=sol%val(k,:)
       case(1)
          if (allocated(mesh%cell(k)%polCoef)) deallocate(mesh%cell(k)%polCoef)
          allocate(mesh%cell(k)%polcoef(order*(order-1)/2+order-1,sol%nvar))
          call reconstruct(mesh,sol,k,order,scheme,gauss_weight,period,mesh%cell(k)%polTest)
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
                call evaluate(mesh,sol,mesh%cell(k)%polTest,k,Xg,Yg,u1)
                sol_interp(k,:)=sol_interp(k,:)+u1*gauss_weight(p1)*gauss_weight(p2)/4.0_dp
                Xg=xc+dx+gauss_point(p1)*dx
                Yg=yc-dy+gauss_point(p2)*dy
                call evaluate(mesh,sol,mesh%cell(k)%polTest,k,Xg,Yg,u2)
                sol_interp(k+mesh%nc,:)=sol_interp(k+mesh%nc,:)+u2*gauss_weight(p1)*gauss_weight(p2)/4.0_dp
                Xg=xc-dx+gauss_point(p1)*dx
                Yg=yc+dy+gauss_point(p2)*dy
                call evaluate(mesh,sol,mesh%cell(k)%polTest,k,Xg,Yg,u3)
                sol_interp(k+2*mesh%nc,:)=sol_interp(k+2*mesh%nc,:)+u3*gauss_weight(p1)*gauss_weight(p2)/4.0_dp
                Xg=xc+dx+gauss_point(p1)*dx
                Yg=yc+dy+gauss_point(p2)*dy
                call evaluate(mesh,sol,mesh%cell(k)%polTest,k,Xg,Yg,u4)
                sol_interp(k+3*mesh%nc,:)=sol_interp(k+3*mesh%nc,:)+u4*gauss_weight(p1)*gauss_weight(p2)/4.0_dp
             enddo
          enddo
       end select
    enddo

    deallocate(u1,u2,u3,u4)
          
    return
  end subroutine interp

  subroutine new_sol(mesh,order,quadrants,sol)
    type(meshStruct), intent(inout) :: mesh
    integer, intent(in) :: order
    type(c_ptr), intent(in) :: quadrants
    type(solStruct), intent(inout) :: sol
    type(c_ptr) :: C_sol,C_pol
    real(dp), dimension(:), pointer :: F_sol,F_pol
    integer :: k,isol,i,size_pol

    size_pol=order*(order-1)/2+order-1
    call p4_new_sol(quadrants,C_sol,C_pol)
    call c_f_pointer(C_sol,F_sol,(/sol%nvar*mesh%nc/))
    call c_f_pointer(C_pol,F_pol,(/sol%nvar*mesh%nc*size_pol/))
    do k=1,mesh%nc
       allocate(mesh%cell(k)%polCoef(size_pol,sol%nvar))
       do isol=1,sol%nvar
          sol%val(k,isol)=F_sol((isol-1)*mesh%nc+k)
          do i=1,size_pol
             mesh%cell(k)%polCoef(i,isol)=F_pol((isol-1)*size_pol*mesh%nc+(i-1)*mesh%nc+k)
          enddo
       enddo
    enddo
    call p4_free(C_sol)
    call p4_free(C_pol)

    return
  end subroutine new_sol

  subroutine adapt_zone(mesh,sol,level,minlevel,maxlevel,coarsen_recursive,refine_recursive,period,sol_coarsen,sol_refine)
    type(meshStruct), intent(inout) :: mesh
    type(solStruct), intent(in) :: sol
    integer, intent(in) :: level
    integer, intent(out) :: minlevel,maxlevel,coarsen_recursive,refine_recursive
    logical, dimension(2), intent(in) :: period
    integer, dimension(:), intent(inout) :: sol_coarsen,sol_refine
    integer :: k,isol
    real(dp) :: xmin,xmax,ymin,ymax
    if(.false.)print*,sol%nvar,period

    minlevel=level
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

  subroutine adapt_sinus(mesh,sol,level,minlevel,maxlevel,coarsen_recursive,refine_recursive,period,sol_coarsen,sol_refine)
    type(meshStruct), intent(inout) :: mesh
    type(solStruct), intent(in) :: sol
    integer, intent(in) :: level
    integer, intent(out) :: minlevel,maxlevel,coarsen_recursive,refine_recursive
    logical, dimension(2), intent(in) :: period
    integer, dimension(:), intent(inout) :: sol_coarsen,sol_refine
    integer :: k,isol
    if(.false.)print*,period

    minlevel=level
    maxlevel=level+1
    coarsen_recursive=0
    refine_recursive=0
    isol=1
    
    do k=1,mesh%nc
       sol_coarsen(k)=0
       sol_refine(k)=0
       if (sol%val(k,isol)<0.5_dp.and.mesh%cell(k)%level>minlevel) sol_coarsen(k)=1
       if (sol%val(k,isol)>0.5_dp) sol_refine(k)=1
    enddo
    
    return
  end subroutine adapt_sinus

  subroutine adapt_vortex(mesh,sol,level,minlevel,maxlevel,coarsen_recursive,refine_recursive,period,sol_coarsen,sol_refine)
    type(meshStruct), intent(inout) :: mesh
    type(solStruct), intent(in) :: sol
    integer, intent(in) :: level
    integer, intent(out) :: minlevel,maxlevel,coarsen_recursive,refine_recursive
    logical, dimension(2), intent(in) :: period
    integer, dimension(:), intent(inout) :: sol_coarsen,sol_refine
    integer :: k,isol,scheme
    real(dp) :: grad

    minlevel=level-1
    maxlevel=level+1
    coarsen_recursive=0
    refine_recursive=0
    isol=1
    scheme=2

    do k=1,mesh%nc
       sol_coarsen(k)=0
       sol_refine(k)=0
       call reconstruct(mesh,sol,k,2,scheme,gauss_weight2,period,mesh%cell(k)%polTest)
       grad=sqrt(mesh%cell(k)%polTest(1,isol)**2+mesh%cell(k)%polTest(2,isol)**2)
       if (grad<0.03_dp.and.mesh%cell(k)%level>minlevel) sol_coarsen(k)=1   !0.002
       if (grad>0.05_dp) sol_refine(k)=1   !0.004
       deallocate(mesh%cell(k)%polTest)
    enddo
    
    return
  end subroutine adapt_vortex

  subroutine adapt_sod(mesh,sol,level,minlevel,maxlevel,coarsen_recursive,refine_recursive,period,sol_coarsen,sol_refine)
    type(meshStruct), intent(inout) :: mesh
    type(solStruct), intent(in) :: sol
    integer, intent(in) :: level
    integer, intent(out) :: minlevel,maxlevel,coarsen_recursive,refine_recursive
    logical, dimension(2), intent(in) :: period
    integer, dimension(:), intent(inout) :: sol_coarsen,sol_refine
    integer :: k,isol,scheme
    real(dp) :: grad

    minlevel=level-1
    maxlevel=level+1
    coarsen_recursive=0
    refine_recursive=0
    isol=1
    scheme=2
    
    do k=1,mesh%nc
       sol_coarsen(k)=0
       sol_refine(k)=0
       call reconstruct(mesh,sol,k,2,scheme,gauss_weight2,period,mesh%cell(k)%polTest)
       grad=sqrt(mesh%cell(k)%polTest(1,isol)**2+mesh%cell(k)%polTest(2,isol)**2)
       if (grad<0.8_dp.and.mesh%cell(k)%level>minlevel) sol_coarsen(k)=1
       if (grad>1.0_dp) sol_refine(k)=1
       deallocate(mesh%cell(k)%polTest)
    enddo
    
    return
  end subroutine adapt_sod

  subroutine adapt_sod2D(mesh,sol,level,minlevel,maxlevel,coarsen_recursive,refine_recursive,period,sol_coarsen,sol_refine)
    type(meshStruct), intent(inout) :: mesh
    type(solStruct), intent(in) :: sol
    integer, intent(in) :: level
    integer, intent(out) :: minlevel,maxlevel,coarsen_recursive,refine_recursive
    logical, dimension(2), intent(in) :: period
    integer, dimension(:), intent(inout) :: sol_coarsen,sol_refine
    integer :: k,isol,scheme
    real(dp) :: lap

    minlevel=level-1
    maxlevel=level+1
    coarsen_recursive=0
    refine_recursive=0
    isol=1
    scheme=2
    
    do k=1,mesh%nc
       sol_coarsen(k)=0
       sol_refine(k)=0
       call reconstruct(mesh,sol,k,3,scheme,gauss_weight3,period,mesh%cell(k)%polTest)
       lap=abs(mesh%cell(k)%polTest(3,isol)+mesh%cell(k)%polTest(5,isol))
       if (lap<3.0_dp.and.mesh%cell(k)%level>minlevel) sol_coarsen(k)=1
       if (lap>6.0_dp) sol_refine(k)=1
       deallocate(mesh%cell(k)%polTest)
    enddo
    
    return
  end subroutine adapt_sod2D
    

end module AMR
