module AMR
  
  use constant
  use types
  use ISO_C_BINDING

contains

  subroutine adapt(p4est,quadrants,mesh,sol)
    type(c_ptr), intent(in) :: p4est,quadrants
    type(meshStruct), intent(in) :: mesh
    type(solStruct) :: sol
    type(c_ptr) :: refine_fn,init_fn
    character(len=20), pointer :: f_refine,f_init
    integer :: refine_recursive
    type(c_ptr) :: C_sol
    real(dp), dimension(:,:), pointer :: F_sol
    
    allocate(f_refine,f_init)
    allocate(F_sol(mesh%nc,sol%nvar))
    
    f_refine='refine_test'//c_null_char
    f_init='null'//c_null_char
    refine_recursive=0
    refine_fn=c_loc(f_refine)
    init_fn=c_loc(f_init)
    F_sol=sol%val
    C_sol=c_loc(F_sol)
    call p4_refine(p4est,quadrants,C_sol,sol%nvar,refine_recursive,refine_fn,init_fn)
    
    deallocate(f_refine,f_init)

    return
  end subroutine adapt

end module AMR
