module AMR
  
  use constant
  use types
  use reconstruction
  use ISO_C_BINDING

contains

  subroutine adapt(fn_adapt,p4est,quadrants,mesh,sol)
    procedure (sub_adapt), pointer, intent(in) :: fn_adapt
    type(c_ptr), intent(in) :: p4est,quadrants
    type(meshStruct), intent(inout) :: mesh
    type(solStruct), intent(in) :: sol
    type(c_ptr) :: C_coarsen_fn,C_refine_fn,C_init_fn
    character(len=20), pointer :: f_coarsen,f_refine,f_init
    type(c_ptr) :: C_sol,C_sol_coarsen,C_sol_refine
    real(dp), dimension(:,:), pointer :: F_sol
    integer, dimension(:), pointer :: sol_coarsen,sol_refine
    integer :: maxlevel,coarsen_recursive,refine_recursive
    
    allocate(f_coarsen,f_refine,f_init)
    allocate(F_sol(mesh%nc,sol%nvar))
    allocate(sol_coarsen(mesh%nc),sol_refine(mesh%nc))
    
    call fn_adapt(mesh,sol,maxlevel,coarsen_recursive,refine_recursive,sol_coarsen,sol_refine)
    
    f_init='null'//c_null_char
    C_coarsen_fn=c_loc(f_coarsen)
    C_refine_fn=c_loc(f_refine)
    C_init_fn=c_loc(f_init)
    F_sol=sol%val
    C_sol=c_loc(F_sol)
    C_sol_coarsen=c_loc(sol_coarsen)
    C_sol_refine=c_loc(sol_refine)
    
    call p4_adapt(p4est,quadrants,C_sol,sol%nvar,C_sol_coarsen,C_sol_refine,maxlevel, &
         coarsen_recursive,refine_recursive)
    
    deallocate(f_coarsen,f_refine,f_init)

    return
  end subroutine adapt

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

  subroutine adapt_sinus(mesh,sol,maxlevel,coarsen_recursive,refine_recursive,sol_coarsen,sol_refine)
    type(meshStruct), intent(inout) :: mesh
    type(solStruct), intent(in) :: sol
    integer, intent(out) :: maxlevel,coarsen_recursive,refine_recursive
    integer, dimension(:), intent(inout) :: sol_coarsen,sol_refine
    integer :: k,isol

    maxlevel=7
    coarsen_recursive=0
    refine_recursive=0
    isol=1
    do k=1,mesh%nc
       sol_coarsen(k)=0
       sol_refine(k)=0
       if (sol%val(k,isol)<0.5_dp.and.mesh%cell(k)%level>6) sol_coarsen(k)=1
       if (sol%val(k,isol)>0.5_dp) sol_refine(k)=1
    enddo
    
    return
  end subroutine adapt_sinus

  subroutine adapt_vortex(mesh,sol,maxlevel,coarsen_recursive,refine_recursive,sol_coarsen,sol_refine)
    type(meshStruct), intent(inout) :: mesh
    type(solStruct), intent(in) :: sol
    integer, intent(out) :: maxlevel,coarsen_recursive,refine_recursive
    integer, dimension(:), intent(inout) :: sol_coarsen,sol_refine
    integer :: k,isol
    real(dp) :: grad

    maxlevel=7
    coarsen_recursive=0
    refine_recursive=0
    isol=1
    do k=1,mesh%nc
       sol_coarsen(k)=0
       sol_refine(k)=0
       call reconstruct(mesh,sol,k,2,gauss_weight2)
       grad=sqrt(mesh%cell(k)%polCoef(1,isol)**2+mesh%cell(k)%polCoef(2,isol)**2)
       if (grad<0.01_dp.and.mesh%cell(k)%level>5) sol_coarsen(k)=1
       if (grad>0.1_dp) sol_refine(k)=1
    enddo
    
    return
  end subroutine adapt_vortex
    

end module AMR
