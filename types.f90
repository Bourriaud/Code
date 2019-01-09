module types

  use constant
  
  implicit none
  
  type :: nodeStruct
     real(dp) :: x,y
  end type nodeStruct

  type :: edgeStruct
     integer :: cell1,cell2
     integer :: dir
     integer :: period
     real(dp) :: length,lengthN
     character(len=20) :: boundType
     real(dp), dimension(:), allocatable :: bound
     real(dp), dimension(:), allocatable :: X_gauss,Y_gauss
     real(dp), dimension(:,:), allocatable :: flux   !flux(gauss_point,var)
     logical, dimension(:), allocatable :: flux_acc
     integer :: deg
     integer, dimension(2) :: sub
  end type edgeStruct

  type :: cellStruct
     real(dp) :: dx,dy
     real(dp) :: xc,yc
     real(dp), dimension(:,:), allocatable :: polCoef
     integer, dimension(4) :: corner
     integer, dimension(:), allocatable :: node
     integer, dimension(:), allocatable :: neigh
     integer, dimension(:), allocatable :: edge
     real(dp), dimension(:), allocatable :: X_gauss,Y_gauss
     logical :: accept
     integer :: deg
  end type cellStruct
  
  type :: meshStruct
     integer :: np,ne,nc
     type(nodeStruct), dimension(:), allocatable :: node
     type(edgeStruct), dimension(:), allocatable :: edge
     type(cellStruct), dimension(:), allocatable :: cell
  end type meshStruct

  type :: solStruct
     integer :: nvar,nsolUser
     character(len=20), dimension(:), allocatable :: name,name_user
     integer, dimension(:), allocatable :: var_user
     real(dp), dimension(:,:), allocatable :: val,user
     integer, dimension(:,:), allocatable :: conserv_var   !conserv_var(isol,1/2)=index of the first/last component of the conservative vector associated to isol
  end type solStruct

  abstract interface

     subroutine sub_IC (x,y,S)
       use constant
       real(dp), intent(in) :: x,y
       real(dp), dimension(:), intent(inout) :: S
     end subroutine sub_IC

     subroutine sub_BC (nvar,mesh)
       use constant
       import meshStruct
       integer, intent(in) :: nvar
       type(meshStruct), intent(inout) :: mesh
     end subroutine sub_BC
     
     subroutine sub_f (u,f)
       use constant
       real(dp), dimension(:), intent(in) :: u
       real(dp), dimension(:,:), intent(inout) :: f
     end subroutine sub_f
     
     subroutine sub_flux (u1,u2,f_equa,dir,F)
       use constant
       real(dp), dimension(:), intent(in) :: u1,u2
       procedure (sub_f), pointer, intent(in) :: f_equa
       integer, intent(in) :: dir
       real(dp), dimension(:), intent(inout) :: F
     end subroutine sub_flux

     subroutine sub_speed (u1,u2,f_equa,Smax)
       use constant
       real(dp), dimension(:), intent(in) :: u1,u2
       procedure (sub_f), pointer, intent(in) :: f_equa
       real(dp), dimension(2), intent(out) :: Smax
     end subroutine sub_speed

     subroutine sub_time (mesh,sol,str_equa,f_equa,flux,speed,order,cfl,t,n,tf, &
          L_str_criteria,L_var_criteria,L_eps,gauss_weight)
       use constant
       import meshStruct
       import solStruct
       type(meshStruct), intent(inout) :: mesh
       type(solStruct), intent(inout) :: sol
       character(len=20), intent(in) :: str_equa
       procedure (sub_f), pointer, intent(in) :: f_equa
       procedure (sub_flux), pointer, intent(in) :: flux
       procedure (sub_speed), pointer, intent(in) :: speed
       integer, intent(in) :: order,n
       real(dp), intent(in) :: cfl,tf
       real(dp), intent(inout) :: t
       character(len=20), dimension(:), intent(in) :: L_str_criteria
       integer, dimension(:), intent(in) :: L_var_criteria
       real(dp), dimension(:), intent(in) :: L_eps
       real(dp), dimension(:), intent(in) :: gauss_weight
     end subroutine sub_time

     subroutine sub_criteria(mesh,sol,sol2,k,isol,eps,gauss_weight,str_equa,accept)
       use constant
       import meshStruct
       import solStruct
       type(meshStruct), intent(inout) :: mesh
       type(solStruct), intent(in) :: sol,sol2
       integer, intent(in) :: k,isol
       real(dp), intent(in) :: eps
       real(dp), dimension(:), intent(in) :: gauss_weight
       character(len=20), intent(in) :: str_equa
       logical, intent(out) :: accept
     end subroutine sub_criteria

     subroutine sub_exactsol(x,y,t,s)
       use constant
       real(dp), intent(in) :: x,y,t
       real(dp), intent(out) :: s
     end subroutine sub_exactsol
     
  end interface

  interface

     subroutine p4_new(level,p4est) bind(C)
       use, intrinsic :: ISO_C_BINDING
       integer(c_int), value, intent(in) :: level
       type(c_ptr), intent(out) :: p4est
     end subroutine p4_new
     
     subroutine p4_build_mesh(p4est,tt,mesh,quadrants,nodes,edges,np,nc,ne) bind(C)
       use, intrinsic :: ISO_C_BINDING
       type(c_ptr), value, intent(in) :: p4est
       type(c_ptr), intent(out) :: mesh,quadrants,nodes,edges
       integer(c_int), intent(out) :: tt
       integer(c_int), intent(out) :: np,nc,ne
     end subroutine p4_build_mesh

     subroutine p4_get_node(p4est,tt,nodes,k,X,Y) bind(C)
       use, intrinsic :: ISO_C_BINDING
       type(c_ptr), value, intent(in) :: p4est,nodes
       integer(c_int), value, intent(in) :: tt
       integer(c_int), value, intent(in) :: k
       real(c_double), intent(out) :: X,Y
     end subroutine p4_get_node

     subroutine p4_get_cell(p4est,mesh,tt,quadrants,nodes,edges,k,xc,yc,dx,dy, &
          corners,Nneigh,neighbors,Nnodes,cell_nodes,N_edge,level) bind(C)
       use, intrinsic :: ISO_C_BINDING
       type(c_ptr), value, intent(in) :: p4est,mesh,quadrants,nodes,edges
       integer(c_int), value, intent(in) :: tt
       integer(c_int), value, intent(in) :: k
       real(c_double), intent(out) :: xc,yc,dx,dy
       type(c_ptr), intent(out) :: corners,neighbors,cell_nodes
       integer(c_int), intent(out) :: Nneigh,Nnodes,N_edge,level
     end subroutine p4_get_cell

     subroutine p4_get_edge(p4est,mesh,quadrants,edges,k,i,iedge_out,nedge,cell1,cell2,sub,period) bind(C)
       use, intrinsic :: ISO_C_BINDING
       type(c_ptr), value, intent(in) :: p4est,mesh,quadrants,edges
       integer(c_int), value, intent(in) :: k,i
       integer(c_int), intent(out) :: nedge
       type(c_ptr), intent(out) :: iedge_out,cell1,cell2,sub,period
     end subroutine p4_get_edge

     subroutine p4_adapt(p4est,quadrants,sol,nsol,maxlevel,coarsen_recursive,coarsen_fn, &
          refine_recursive,refine_fn,init_fn) bind(C)
       use, intrinsic :: ISO_C_BINDING
       type(c_ptr), value, intent(in) :: p4est,quadrants,sol,coarsen_fn,refine_fn,init_fn
       integer(c_int), value, intent(in) :: nsol,maxlevel,coarsen_recursive,refine_recursive
     end subroutine p4_adapt

     subroutine p4_new_sol(quadrants,sol) bind(C)
       use, intrinsic :: ISO_C_BINDING
       type(c_ptr), value, intent(in) :: quadrants
       type(c_ptr), intent(out) :: sol
     end subroutine p4_new_sol

     subroutine p4_free(ptr) bind(C)
       use, intrinsic :: ISO_C_BINDING
       type(c_ptr), value, intent(in) :: ptr
     end subroutine p4_free
     
  end interface

end module types
