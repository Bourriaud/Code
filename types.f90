module types

  implicit none

  type :: meshStruct
     integer :: nx,ny
     real, dimension(:,:), allocatable :: X,Y
     real, dimension(:,:), allocatable :: CX,CY
     character(len=20), dimension(:,:), allocatable :: boundType
     real, dimension(:,:,:), allocatable :: bound
  end type meshStruct

  type :: solStruct
     integer :: nvar,nsolUser
     character(len=20), dimension(:), allocatable :: name, nameUser
     real, dimension(:,:,:), allocatable :: val, user
  end type solStruct

end module types
