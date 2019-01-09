module AMR
  
  use constant
  use types
  use ISO_C_BINDING

contains

  subroutine adapt(p4est,quadrants,mesh,sol,maxlevel,coarsen_recursive,refine_recursive, &
       coarsen_fn,refine_fn)
    type(c_ptr), intent(in) :: p4est,quadrants
    type(meshStruct), intent(in) :: mesh
    type(solStruct) :: sol
    integer, intent(in) :: maxlevel,coarsen_recursive,refine_recursive
    character(len=20), intent(in) :: coarsen_fn,refine_fn
    type(c_ptr) :: C_coarsen_fn,C_refine_fn,C_init_fn
    character(len=20), pointer :: f_coarsen,f_refine,f_init
    type(c_ptr) :: C_sol
    real(dp), dimension(:,:), pointer :: F_sol
    
    allocate(f_coarsen,f_refine,f_init)
    allocate(F_sol(mesh%nc,sol%nvar))

    f_coarsen=trim(coarsen_fn)//c_null_char
    f_refine=trim(refine_fn)//c_null_char
    f_init='null'//c_null_char
    C_coarsen_fn=c_loc(f_coarsen)
    C_refine_fn=c_loc(f_refine)
    C_init_fn=c_loc(f_init)
    F_sol=sol%val
    C_sol=c_loc(F_sol)
    call p4_adapt(p4est,quadrants,C_sol,sol%nvar,maxlevel,coarsen_recursive,C_coarsen_fn, &
         refine_recursive,C_refine_fn,C_init_fn)
    
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

    return
  end subroutine new_sol
    

end module AMR
