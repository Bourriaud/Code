module types

  use constant
  
  implicit none
  
  type :: nodeStruct
     real(dp) :: x,y
     integer, dimension(4) :: connect
  end type nodeStruct

  type :: edgeStruct
     integer :: node1,node2
     integer :: cell1,cell2
     integer :: dir
     real(dp) :: length
     character(len=20) :: boundType
     real(dp), dimension(:), allocatable :: bound
  end type edgeStruct

  type :: cellStruct
     real(dp) :: dx,dy
     real(dp) :: xc,yc
     real(dp), dimension(:,:), allocatable :: polCoef
     integer, dimension(:), allocatable :: neigh
     type(edgeStruct), dimension(:), allocatable :: edge
  end type cellStruct
  
  type :: meshStruct
     integer :: np,ne,nc
     type(nodeStruct), dimension(:), allocatable :: node
     type(edgeStruct), dimension(:), allocatable :: edge
     type(cellStruct), dimension(:), allocatable :: cell
  end type meshStruct

  type :: solStruct
     integer :: nvar,nsolUser
     character(len=20), dimension(:), allocatable :: name,nameUser
     real(dp), dimension(:,:), allocatable :: val,user
  end type solStruct

  abstract interface
     
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

     subroutine sub_time (mesh,sol,f_equa,flux,speed,order,cfl,t,tf,quad_c_alpha,quad_reconstruct)
       use constant
       import meshStruct
       import solStruct
       type(meshStruct), intent(inout) :: mesh
       type(solStruct), intent(inout) :: sol
       procedure (sub_f), pointer, intent(in) :: f_equa
       procedure (sub_flux), pointer, intent(in) :: flux
       procedure (sub_speed), pointer, intent(in) :: speed
       integer, intent(in) :: order
       real(dp), intent(in) :: cfl,tf
       real(dp), intent(inout) :: t
       procedure (quadrature_c_alpha), pointer, intent(in) :: quad_c_alpha
       procedure (quadrature_reconstruction), pointer, intent(in) :: quad_reconstruct
     end subroutine sub_time

     subroutine sub_quadra_t (x,y,t,s)
       use constant
       real(dp), intent(in) :: x,y,t
       real(dp), intent(out) :: s
     end subroutine sub_quadra_t

     subroutine sub_quadra_c_alpha (x,y,c,alpha,s)
       use constant
       real(dp), intent(in) :: x,y
       real(dp), dimension(2), intent(in) :: c
       integer, dimension(2), intent(in) :: alpha
       real(dp), intent(out) :: s
     end subroutine sub_quadra_c_alpha

     subroutine sub_reconstruction (mesh,sol,k,order,quad_c_alpha,x,y,u)
       use constant
       import meshStruct
       import solStruct
       type(meshStruct), intent(in) :: mesh
       type(solStruct), intent(in) :: sol
       integer, intent(in) :: k,order
       procedure (quadrature_c_alpha), pointer, intent(in) :: quad_c_alpha
       real(dp), intent(in) :: x,y
       real(dp), dimension(:), intent(inout) :: u
     end subroutine sub_reconstruction

     subroutine quadrature_t(func,mesh,t,k,int)
       use constant
       import meshStruct
       procedure(sub_quadra_t) :: func
       type(meshStruct), intent(in) :: mesh
       real(dp), intent(in) :: t
       integer, intent(in) :: k
       real(dp), intent(out) :: int
       real(dp) :: x,y
     end subroutine quadrature_t

     subroutine quadrature_c_alpha(func,mesh,c,alpha,k,int)
       use constant
       import meshStruct
       procedure(sub_quadra_c_alpha) :: func
       type(meshStruct), intent(in) :: mesh
       real(dp), dimension(2), intent(in) :: c
       integer, dimension(2), intent(in) :: alpha
       integer, intent(in) :: k
       real(dp), intent(out) :: int
       real(dp) :: x,y
     end subroutine quadrature_c_alpha

     subroutine quadrature_reconstruction(func,mesh,sol,order,quad_c_alpha,normal,k,int)
       use constant
       import meshStruct
       import solStruct
       procedure(sub_reconstruction) :: func
       type(meshStruct), intent(in) :: mesh
       type(solStruct), intent(in) :: sol
       integer, intent(in) :: order,normal,k
       procedure (quadrature_c_alpha), pointer, intent(in) :: quad_c_alpha
       real(dp), dimension(:), intent(inout) :: int
       real(dp) :: x,y,dx,dy
     end subroutine quadrature_reconstruction
     
  end interface

end module types
