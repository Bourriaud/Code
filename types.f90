module types

  implicit none

  type :: meshStruct
     integer :: nx,ny
     real, dimension(:,:), allocatable :: X,Y
  end type meshStruct

  type :: solStruct
     integer :: nsol
     character(len=20), dimension(:), allocatable :: name
     real, dimension(:,:,:), allocatable :: val
  end type solStruct

end module types
