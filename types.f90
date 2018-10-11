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

end module types
