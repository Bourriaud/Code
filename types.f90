module types

  use constant
  
  implicit none
  
  type :: nodeStruct
     real(dp) :: x,y
     integer, dimension(:), allocatable :: connect
  end type nodeStruct

  type :: edgeStruct
     integer :: node1,node2
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

     subroutine sub_BC (nx,ny,nvar,mesh)
       use constant
       import meshStruct
       integer, intent(in) :: nx,ny,nvar
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

end module types
