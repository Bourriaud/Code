module types

  implicit none
  
  type :: nodeStruct
     real :: x,y
     integer, dimension(4) :: connect
  end type nodeStruct

  type :: edgeStruct
     integer :: node1,node2
     integer :: neigh
     integer :: normal    !1=left 2=bottom 3=right 4=top
     character(len=20) :: boundType
     real, dimension(:), allocatable :: bound
  end type edgeStruct

  type :: cellStruct
     real :: dx,dy
     real :: xc,yc
     type(edgeStruct), dimension(:), allocatable :: edge
  end type cellStruct
  
  type :: meshStruct
     integer :: np,nc
     type(nodeStruct), dimension(:), allocatable :: node
     type(cellStruct), dimension(:), allocatable :: cell
  end type meshStruct

  type :: solStruct
     integer :: nvar,nsolUser
     character(len=20), dimension(:), allocatable :: name,nameUser
     real, dimension(:,:), allocatable :: val,user
  end type solStruct

  abstract interface
     
     subroutine sub_f (u,f)
       real, dimension(:), intent(in) :: u
       real, dimension(:,:), intent(inout) :: f
     end subroutine sub_f
     
     subroutine sub_flux (u1,u2,f_equa,dir,F,Smax)
       real, dimension(:), intent(in) :: u1,u2
       procedure (sub_f), pointer, intent(in) :: f_equa
       integer, intent(in) :: dir
       real, dimension(:), intent(inout) :: F
       real, intent(out) :: Smax
     end subroutine sub_flux

     subroutine sub_time (mesh,sol,f_ptr,flux_ptr,cfl,t)
       import meshStruct
       import solStruct
       type(meshStruct), intent(in) :: mesh
       type(solStruct), intent(inout) :: sol
       procedure (sub_f), pointer, intent(in) :: f_ptr
       procedure (sub_flux), pointer, intent(in) :: flux_ptr
       real, intent(in) :: cfl
       real, intent(inout) :: t
     end subroutine sub_time
     
  end interface

end module types
