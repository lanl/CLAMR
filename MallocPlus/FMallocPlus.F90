      module FMallocPlus_mod

#ifndef INT8_KIND_DIGITS
#define INT8_KIND_DIGITS 16
#endif

#ifndef REAL4_KIND_DIGITS
#define REAL4_KIND_DIGITS 6
#endif

#ifndef REAL8_KIND_DIGITS
#define REAL8_KIND_DIGITS 12
#endif

      use, intrinsic :: ISO_C_Binding, only : C_PTR, C_NULL_ptr

      implicit none

      integer, parameter :: INT8   = SELECTED_INT_KIND(INT8_KIND_DIGITS)
      integer, parameter :: REAL4  = SELECTED_REAL_KIND(REAL4_KIND_DIGITS)
      integer, parameter :: REAL8  = SELECTED_REAL_KIND(REAL8_KIND_DIGITS)

      integer(INT8), public, parameter :: HOST_REGULAR_MEMORY   = 000
      integer(INT8), public, parameter :: HOST_MANAGED_MEMORY   = 001
      integer(INT8), public, parameter :: DEVICE_REGULAR_MEMORY = 002
      integer(INT8), public, parameter :: INDEX_ARRAY_MEMORY    = 004
      integer(INT8), public, parameter :: LOAD_BALANCE_MEMORY   = 008
      integer(INT8), public, parameter :: RESTART_DATA          = 016
      integer(INT8), public, parameter :: FORTRAN_DATA          = 032

      public FMallocPlus_type, FMallocPlus_new
      public FMallocPlus_malloc, FMallocPlus_add
      public FMallocPlus_report
      public FMallocPlus_get_mem_size

      interface FMallocPlus_add
        module procedure FMallocPlus_add_1D_int
        module procedure FMallocPlus_add_2D_int
        module procedure FMallocPlus_add_3D_int
        module procedure FMallocPlus_add_4D_int
        module procedure FMallocPlus_add_1D_dble
        module procedure FMallocPlus_add_2D_dble
        module procedure FMallocPlus_add_3D_dble
        module procedure FMallocPlus_add_4D_dble
        module procedure FMallocPlus_add_5D_dble
        module procedure FMallocPlus_add_6D_dble
      end interface FMallocPlus_add

      interface FMallocPlus_malloc
        module procedure FMallocPlus_malloc_1D_int
        module procedure FMallocPlus_malloc_2D_int
        module procedure FMallocPlus_malloc_3D_int
        module procedure FMallocPlus_malloc_4D_int
        module procedure FMallocPlus_malloc_1D_dble
        module procedure FMallocPlus_malloc_2D_dble
        module procedure FMallocPlus_malloc_3D_dble
        module procedure FMallocPlus_malloc_4D_dble
        module procedure FMallocPlus_malloc_5D_dble
        module procedure FMallocPlus_malloc_6D_dble
      end interface

      type FMallocPlus_type
        type(C_PTR) :: object = C_NULL_ptr
      end type FMallocPlus_type

      integer(INT8) :: mem_size = 0

      private

      interface
        function C_MallocPlus_new() result(mem_object) &
            bind(C,name="MallocPlus_new")
          use, intrinsic :: ISO_C_Binding, only : C_PTR
          type(C_PTR) :: mem_object
        end function C_MallocPlus_new

        subroutine C_MallocPlus_report(mem_object) &
            bind(C,name="MallocPlus_memory_report")
          use, intrinsic :: ISO_C_Binding, only : C_PTR
          type(C_PTR), value :: mem_object
        end subroutine C_MallocPlus_report

        subroutine C_MallocPlus_add(mem_object, intptr, nelem, &
            elsize, name, flags) bind(C,NAME="MallocPlus_memory_add")
          use, intrinsic :: ISO_C_Binding, only : C_PTR, C_INT, &
            C_CHAR, C_SIZE_T
          type(C_PTR),       value   :: mem_object
          type(C_PTR),       value   :: intptr
          integer(C_SIZE_T), value   :: nelem, elsize
          character(kind=C_CHAR, len=1) :: name(*)
          integer(C_SIZE_T), value   :: flags
        end subroutine C_MallocPlus_add

        subroutine C_MallocPlus_add_multi(mem_object, intptr, ndim, nelem, &
            elsize, name, flags) bind(C,NAME="MallocPlus_memory_add_nD")
          use, intrinsic :: ISO_C_Binding, only : C_PTR, C_INT, &
            C_CHAR, C_SIZE_T
          type(C_PTR),       value   :: mem_object
          type(C_PTR),       value   :: intptr
          integer(C_INT),    value   :: ndim
          type(C_PTR),       value   :: nelem
          integer(C_SIZE_T), value   :: elsize
          character(kind=C_CHAR, len=1) :: name(*)
          integer(C_SIZE_T), value   :: flags
        end subroutine C_MallocPlus_add_multi
      end interface

      contains

!=======================================================================
      subroutine FMallocPlus_new(mem_object)
        type(FMallocPlus_type), intent(inout) :: mem_object
        mem_object%object = C_MallocPlus_new()
!       print *,"DEBUG FMallocPlus_new ",mem_object%object
      end subroutine FMallocPlus_new
!=======================================================================
      subroutine FMallocPlus_report(mem_object)
        type(FMallocPlus_type), intent(in) :: mem_object
        call C_MallocPlus_report(mem_object%object)
      end subroutine FMallocPlus_report
!=======================================================================
      function FMallocPlus_get_mem_size()
        integer(INT8) :: FMallocPlus_get_mem_size
        FMallocPlus_get_mem_size = mem_size
      end function FMallocPlus_get_mem_size
!=======================================================================
      subroutine FMallocPlus_add_1D_int(mem_object, intptr, nelem, elsize, &
          name, flags)
        use iso_c_binding
        implicit none

        ! Argument definitions
        type(FMallocPlus_type), intent(inout) :: mem_object
        integer, dimension(:), pointer, intent(in) :: intptr
        integer(INT8), intent(in) :: nelem, elsize
        integer(INT8), intent(in) :: flags
        character(*), intent(in) :: name

        ! Local variables
        integer :: low1

        ! Executable code
        low1 = lbound(intptr,1)
        call C_MallocPlus_add(mem_object%object, c_loc(intptr(low1)), nelem, &
          elsize, trim(name)//char(0), flags)
      end subroutine FMallocPlus_add_1D_int
!=======================================================================
      subroutine FMallocPlus_add_2D_int(mem_object, intptr, ndim, nelem, &
          elsize, name, flags)

        use iso_c_binding, only : c_loc

        implicit none

        ! Argument definitions
        type(FMallocPlus_type), intent(inout) :: mem_object
        integer, dimension(:,:), pointer, intent(in) :: intptr
        integer, intent(in) :: ndim
        integer(INT8), pointer, intent(in) :: nelem(:)
        integer(INT8), intent(in) :: elsize
        integer(INT8), intent(in) :: flags
        character(*), intent(in) :: name

        ! Local variables
        integer :: low1, low2

        ! Executable code
        low1 = lbound(intptr,1)
        low2 = lbound(intptr,2)
        call C_MallocPlus_add_multi(mem_object%object, &
          c_loc(intptr(low1,low2)), ndim, &
          c_loc(nelem(1)), elsize, trim(name)//char(0), flags)
      end subroutine FMallocPlus_add_2D_int
!=======================================================================
      subroutine FMallocPlus_add_3D_int(mem_object, intptr, ndim, nelem, &
          elsize, name, flags)

        use iso_c_binding, only : c_loc

        implicit none

        ! Argument definitions
        type(FMallocPlus_type), intent(inout) :: mem_object
        integer, dimension(:,:,:), pointer, intent(in) :: intptr
        integer, intent(in) :: ndim
        integer(INT8), pointer, intent(in) :: nelem(:)
        integer(INT8), intent(in) :: elsize
        integer(INT8), intent(in) :: flags
        character(*), intent(in) :: name

        ! Local variables
        integer :: low1, low2, low3

        ! Executable code
        low1 = lbound(intptr,1)
        low2 = lbound(intptr,2)
        low3 = lbound(intptr,3)
        call C_MallocPlus_add_multi(mem_object%object, &
          c_loc(intptr(low1,low2,low3)), ndim, &
          c_loc(nelem(1)), elsize, trim(name)//char(0), flags)
      end subroutine FMallocPlus_add_3D_int
!=======================================================================
      subroutine FMallocPlus_add_4D_int(mem_object, intptr, ndim, nelem, &
          elsize, name, flags)

        use iso_c_binding, only : c_loc

        implicit none

        ! Argument definitions
        type(FMallocPlus_type), intent(inout) :: mem_object
        integer, dimension(:,:,:,:), pointer, intent(in) :: intptr
        integer, intent(in) :: ndim
        integer(INT8), pointer, intent(in) :: nelem(:)
        integer(INT8), intent(in) :: elsize
        integer(INT8), intent(in) :: flags
        character(*), intent(in) :: name

        ! Local variables
        integer :: low1,low2,low3,low4

        ! Executable code
        low1 = lbound(intptr,1)
        low2 = lbound(intptr,2)
        low3 = lbound(intptr,3)
        low4 = lbound(intptr,4)
        call C_MallocPlus_add_multi(mem_object%object, &
          c_loc(intptr(low1,low2,low3,low4)), ndim, &
          c_loc(nelem(1)), elsize, trim(name)//char(0), flags)
      end subroutine FMallocPlus_add_4D_int
!=======================================================================
      subroutine FMallocPlus_add_1D_dble(mem_object, dbleptr, nelem, elsize, &
          name, flags)

        use iso_c_binding, only : c_loc

        implicit none

        ! Argument definitions
        type(FMallocPlus_type), intent(inout) :: mem_object
        real(REAL8), dimension(:), pointer, intent(in) :: dbleptr
        integer(INT8), intent(in) :: nelem, elsize
        integer(INT8), intent(in) :: flags
        character(*), intent(in) :: name

        ! Local variables
        integer :: low1

        ! Executable code
        low1 = lbound(dbleptr,1)
        ! print *,"DEBUG mem ptr entering add is ",loc(dbleptr),mem_object%object,nelem,elsize
        call C_MallocPlus_add(mem_object%object, c_loc(dbleptr(low1)), nelem, &
          elsize, trim(name)//char(0), flags)
      end subroutine FMallocPlus_add_1D_dble
!=======================================================================
      subroutine FMallocPlus_add_2D_dble(mem_object, dbleptr, ndim, nelem, &
          elsize, name, flags)

        use iso_c_binding, only : c_loc

        implicit none

        ! Argument definitions
        type(FMallocPlus_type), intent(inout) :: mem_object
        real(REAL8), dimension(:,:), pointer, intent(in) :: dbleptr
        integer, intent(in) :: ndim
        integer(INT8), intent(in) :: elsize
        integer(INT8), pointer, intent(in) :: nelem(:)
        integer(INT8), intent(in) :: flags
        character(*), intent(in) :: name

        ! Local variables
        integer :: low1, low2

        ! Executable code
        low1 = lbound(dbleptr,1)
        low2 = lbound(dbleptr,2)
        call C_MallocPlus_add_multi(mem_object%object, &
          c_loc(dbleptr(low1,low2)), ndim, &
          c_loc(nelem(1)), elsize, trim(name)//char(0), flags)
      end subroutine FMallocPlus_add_2D_dble
!=======================================================================
      subroutine FMallocPlus_add_3D_dble(mem_object, dbleptr, ndim, nelem, &
          elsize, name, flags)

        use iso_c_binding, only : c_loc

        implicit none

        ! Argument definitions
        type(FMallocPlus_type), intent(inout) :: mem_object
        real(REAL8), dimension(:,:,:), pointer, intent(in) :: dbleptr
        integer, intent(in) :: ndim
        integer(INT8), pointer, intent(in) :: nelem(:)
        integer(INT8), intent(in) :: elsize
        integer(INT8), intent(in) :: flags
        character(*), intent(in) :: name

        ! Local variables
        integer :: low1, low2, low3

        ! Executable code
        low1 = lbound(dbleptr,1)
        low2 = lbound(dbleptr,2)
        low3 = lbound(dbleptr,3)
        call C_MallocPlus_add_multi(mem_object%object, &
          c_loc(dbleptr(low1,low2,low3)), ndim, &
          c_loc(nelem(1)), elsize, trim(name)//char(0), flags)
      end subroutine FMallocPlus_add_3D_dble
!=======================================================================
      subroutine FMallocPlus_add_4D_dble(mem_object, dbleptr, ndim, nelem, &
          elsize, name, flags)

        use iso_c_binding, only : c_loc

        implicit none

        ! Argument definitions
        type(FMallocPlus_type), intent(inout) :: mem_object
        real(REAL8), dimension(:,:,:,:), pointer, intent(in) :: dbleptr
        integer, intent(in) :: ndim
        integer(INT8), pointer, intent(in) :: nelem(:)
        integer(INT8), intent(in) :: elsize
        integer(INT8), intent(in) :: flags
        character(*), intent(in) :: name

        ! Local variables
        integer :: low1, low2, low3, low4

        ! Executable code
        low1 = lbound(dbleptr,1)
        low2 = lbound(dbleptr,2)
        low3 = lbound(dbleptr,3)
        low4 = lbound(dbleptr,4)
        call C_MallocPlus_add_multi(mem_object%object, &
          c_loc(dbleptr(low1,low2,low3,low4)), ndim, &
          c_loc(nelem(1)), elsize, trim(name)//char(0), flags)
      end subroutine FMallocPlus_add_4D_dble
!=======================================================================
      subroutine FMallocPlus_add_5D_dble(mem_object, dbleptr, ndim, nelem, &
          elsize, name, flags)

        use iso_c_binding, only : c_loc

        implicit none

        ! Argument definitions
        type(FMallocPlus_type), intent(inout) :: mem_object
        real(REAL8), dimension(:,:,:,:,:), pointer, intent(in) :: dbleptr
        integer, intent(in) :: ndim
        integer(INT8), pointer, intent(in) :: nelem(:)
        integer(INT8), intent(in) :: elsize
        integer(INT8), intent(in) :: flags
        character(*), intent(in) :: name

        ! Local variables
        integer :: low1, low2, low3, low4, low5

        ! Executable code
        low1 = lbound(dbleptr,1)
        low2 = lbound(dbleptr,2)
        low3 = lbound(dbleptr,3)
        low4 = lbound(dbleptr,4)
        low5 = lbound(dbleptr,5)
        call C_MallocPlus_add_multi(mem_object%object, &
          c_loc(dbleptr(low1,low2,low3,low4,low5)), ndim, &
          c_loc(nelem(1)), elsize, trim(name)//char(0), flags)
      end subroutine FMallocPlus_add_5D_dble
!=======================================================================
      subroutine FMallocPlus_add_6D_dble(mem_object, dbleptr, ndim, nelem, &
          elsize, name, flags)

        use iso_c_binding, only : c_loc

        implicit none

        ! Argument definitions
        type(FMallocPlus_type), intent(inout) :: mem_object
        real(REAL8), dimension(:,:,:,:,:,:), pointer, intent(in) :: dbleptr
        integer, intent(in) :: ndim
        integer(INT8), pointer, intent(in) :: nelem(:)
        integer(INT8), intent(in) :: elsize
        integer(INT8), intent(in) :: flags
        character(*), intent(in) :: name

        ! Local variables
        integer :: low1, low2, low3, low4, low5, low6

        ! Executable code
        low1 = lbound(dbleptr,1)
        low2 = lbound(dbleptr,2)
        low3 = lbound(dbleptr,3)
        low4 = lbound(dbleptr,4)
        low5 = lbound(dbleptr,5)
        low6 = lbound(dbleptr,6)
        call C_MallocPlus_add_multi(mem_object%object, &
          c_loc(dbleptr(low1,low2,low3,low4,low5,low6)), ndim, &
          c_loc(nelem(1)), elsize, trim(name)//char(0), flags)
      end subroutine FMallocPlus_add_6D_dble
!=======================================================================
      subroutine FMallocPlus_malloc_1D_int(mem_object, intptr, &
          nlbound, nubound, &
          name, flags)

        use iso_c_binding, only : c_loc

        implicit none

        ! Argument definitions
        type(FMallocPlus_type), intent(inout) :: mem_object
        integer, dimension(:), pointer, intent(inout) :: intptr
        integer, intent(in) :: nlbound, nubound
        integer(INT8), intent(in) :: flags
        character(*), intent(in) :: name

        ! Local variables
        integer :: ierr
        integer(INT8) :: msize
        integer(INT8), parameter :: elsize = 4

        ! Executable code

        msize = (nubound-nlbound+1)
        mem_size = mem_size + msize * elsize
        allocate(intptr(nlbound:nubound), stat=ierr)
        call allocate_error(ierr, name, mem_size)

        call C_MallocPlus_add(mem_object%object, c_loc(intptr(nlbound)), msize, &
          elsize, trim(name)//char(0), flags)

      end subroutine FMallocPlus_malloc_1D_int
      
!=======================================================================
      subroutine FMallocPlus_malloc_2D_int(mem_object, intptr, &
          nlbound, nubound, mlbound, mubound, &
          name, flags)

        use iso_c_binding, only : c_loc

        implicit none

        ! Argument definitions
        type(FMallocPlus_type), intent(inout) :: mem_object
        integer, dimension(:,:), pointer, intent(inout) :: intptr
        integer, intent(in) :: nlbound, nubound, mlbound, mubound
        integer(INT8), intent(in) :: flags
        character(*), intent(in) :: name

        ! Local variables
        integer :: ierr
        integer(INT8), pointer :: msize(:)
        integer, parameter :: ndim = 2
        integer(INT8), parameter :: elsize = 4

        ! Executable code

        allocate(msize(ndim), stat=ierr)
        call allocate_error(ierr, name, mem_size)

        msize(1) = nubound-nlbound+1
        msize(2) = mubound-mlbound+1
        mem_size = mem_size + msize(1)*msize(2)*elsize

        allocate(intptr(nlbound:nubound, mlbound:mubound), stat=ierr)
        call allocate_error(ierr, name, mem_size)

        call C_MallocPlus_add_multi(mem_object%object, &
          c_loc(intptr(nlbound,mlbound)), ndim, &
          c_loc(msize(1)), elsize, trim(name)//char(0), flags)

        deallocate(msize)

      end subroutine FMallocPlus_malloc_2D_int
!=======================================================================
      subroutine FMallocPlus_malloc_3D_int(mem_object, intptr, &
          nlbound, nubound, mlbound, mubound, &
          llbound, lubound, &
          name, flags)

        use iso_c_binding, only : c_loc

        implicit none

        ! Argument definitions
        type(FMallocPlus_type), intent(inout) :: mem_object
        integer, dimension(:,:,:), pointer, intent(inout) :: intptr
        integer, intent(in) :: nlbound, nubound, mlbound, mubound, &
                               llbound, lubound
        integer(INT8), intent(in) :: flags
        character(*), intent(in) :: name

        ! Local variables
        integer :: ierr
        integer(INT8), pointer :: msize(:)
        integer, parameter :: ndim = 3
        integer(INT8), parameter :: elsize = 4

        ! Executable code

        allocate(msize(ndim), stat=ierr)
        call allocate_error(ierr, name, mem_size)

        msize(1) = nubound-nlbound+1
        msize(2) = mubound-mlbound+1
        msize(3) = lubound-llbound+1
        mem_size = mem_size + msize(1)*msize(2)*msize(3)*elsize

        allocate(intptr(nlbound:nubound, mlbound:mubound, &
                        llbound:lubound), stat=ierr)
        call allocate_error(ierr, name, mem_size)

        call C_MallocPlus_add_multi(mem_object%object, &
          c_loc(intptr(mlbound,mlbound,llbound)), ndim, &
          c_loc(msize(1)), elsize, trim(name)//char(0), flags)

        deallocate(msize)

      end subroutine FMallocPlus_malloc_3D_int
!=======================================================================
      subroutine FMallocPlus_malloc_4D_int(mem_object, intptr, &
          nlbound, nubound, mlbound, mubound, &
          llbound, lubound, klbound, kubound, &
          name, flags)

        use iso_c_binding, only : c_loc

        implicit none

        ! Argument definitions
        type(FMallocPlus_type), intent(inout) :: mem_object
        integer, dimension(:,:,:,:), pointer, intent(inout) :: intptr
        integer, intent(in) :: nlbound, nubound, mlbound, mubound, &
                               llbound, lubound, klbound, kubound
        integer(INT8), intent(in) :: flags
        character(*), intent(in) :: name

        ! Local variables
        integer :: ierr
        integer(INT8), pointer :: msize(:)
        integer, parameter :: ndim = 4
        integer(INT8), parameter :: elsize = 4

        ! Executable code

        allocate(msize(ndim), stat=ierr)
        call allocate_error(ierr, name, mem_size)

        msize(1) = nubound-nlbound+1
        msize(2) = mubound-mlbound+1
        msize(3) = lubound-llbound+1
        msize(4) = kubound-klbound+1
        mem_size = mem_size + msize(1)*msize(2)*msize(3)*msize(4)*elsize

        allocate(intptr(nlbound:nubound, mlbound:mubound, &
                        llbound:lubound, klbound:kubound), stat=ierr)
        call allocate_error(ierr, name, mem_size)

        call C_MallocPlus_add_multi(mem_object%object, &
          c_loc(intptr(nlbound,mlbound,llbound,klbound)), ndim, &
          c_loc(msize(1)), elsize, trim(name)//char(0), flags)

        deallocate(msize)

      end subroutine FMallocPlus_malloc_4D_int
!=======================================================================
      subroutine FMallocPlus_malloc_1D_dble(mem_object, dbleptr, &
          nlbound, nubound, &
          name, flags)

        use iso_c_binding, only : c_loc

        implicit none

        ! Argument definitions
        type(FMallocPlus_type), intent(inout) :: mem_object
        real(REAL8), dimension(:), pointer, intent(inout) :: dbleptr
        integer, intent(in) :: nlbound, nubound
        integer(INT8), intent(in) :: flags
        character(*), intent(in) :: name

        ! Local variables
        integer :: ierr
        integer(INT8) :: msize
        integer(INT8), parameter :: elsize = 8

        ! Executable code

        msize = (nubound-nlbound+1)
        mem_size = mem_size + msize * elsize
        allocate(dbleptr(nlbound:nubound), stat=ierr)
        call allocate_error(ierr, name, mem_size)

        !call FMallocPlus_add(mem_object, dbleptr, msize, elsize, name, flags)
        call C_MallocPlus_add(mem_object%object, &
          c_loc(dbleptr(nlbound)), msize, &
          elsize, trim(name)//char(0), flags)

      end subroutine FMallocPlus_malloc_1D_dble
      
!=======================================================================
      subroutine FMallocPlus_malloc_2D_dble(mem_object, dbleptr, &
          nlbound, nubound, mlbound, mubound, &
          name, flags)

        use iso_c_binding, only : c_loc

        implicit none

        ! Argument definitions
        type(FMallocPlus_type), intent(inout) :: mem_object
        real(REAL8), dimension(:,:), pointer, intent(inout) :: dbleptr
        integer, intent(in) :: nlbound, nubound, mlbound, mubound
        integer(INT8), intent(in) :: flags
        character(*), intent(in) :: name

        ! Local variables
        integer :: ierr
        integer(INT8), pointer :: msize(:)
        integer, parameter :: ndim = 2
        integer(INT8), parameter :: elsize = 8

        ! Executable code

        allocate(msize(ndim), stat=ierr)
        call allocate_error(ierr, name, mem_size)

        msize(1) = nubound-nlbound+1
        msize(2) = mubound-mlbound+1
        mem_size = mem_size + msize(1)*msize(2)*elsize

        allocate(dbleptr(nlbound:nubound, mlbound:mubound), stat=ierr)
        call allocate_error(ierr, name, mem_size)

        call C_MallocPlus_add_multi(mem_object%object, &
          c_loc(dbleptr(nlbound,mlbound)), ndim, &
          c_loc(msize(1)), elsize, trim(name)//char(0), flags)

        deallocate(msize)

      end subroutine FMallocPlus_malloc_2D_dble
!=======================================================================
      subroutine FMallocPlus_malloc_3D_dble(mem_object, dbleptr, &
          nlbound, nubound, mlbound, mubound, &
          llbound, lubound, &
          name, flags)

        use iso_c_binding, only : c_loc

        implicit none

        ! Argument definitions
        type(FMallocPlus_type), intent(inout) :: mem_object
        real(REAL8), dimension(:,:,:), pointer, intent(inout) :: dbleptr
        integer, intent(in) :: nlbound, nubound
        integer, intent(in) :: mlbound, mubound
        integer, intent(in) :: llbound, lubound
        integer(INT8), intent(in) :: flags
        character(*), intent(in) :: name

        ! Local variables
        integer :: ierr
        integer(INT8), pointer :: msize(:)
        integer, parameter :: ndim = 3
        integer(INT8), parameter :: elsize = 8

        ! Executable code

        allocate(msize(ndim), stat=ierr)
        call allocate_error(ierr, name, mem_size)

        msize(1) = nubound-nlbound+1
        msize(2) = mubound-mlbound+1
        msize(3) = lubound-llbound+1
        mem_size = mem_size + msize(1)*msize(2)*msize(3)*elsize

        allocate(dbleptr(nlbound:nubound, mlbound:mubound, &
                         llbound:lubound), stat=ierr)
        call allocate_error(ierr, name, mem_size)

        call C_MallocPlus_add_multi(mem_object%object, &
          c_loc(dbleptr(nlbound,mlbound,llbound)), ndim, &
          c_loc(msize(1)), elsize, trim(name)//char(0), flags)

        deallocate(msize)

      end subroutine FMallocPlus_malloc_3D_dble
!=======================================================================
      subroutine FMallocPlus_malloc_4D_dble(mem_object, dbleptr, &
          nlbound, nubound, mlbound, mubound, &
          llbound, lubound, klbound, kubound, &
          name, flags)

        use iso_c_binding, only : c_loc

        implicit none

        ! Argument definitions
        type(FMallocPlus_type), intent(inout) :: mem_object
        real(REAL8), dimension(:,:,:,:), pointer, intent(inout) :: dbleptr
        integer, intent(in) :: nlbound, nubound
        integer, intent(in) :: mlbound, mubound
        integer, intent(in) :: llbound, lubound
        integer, intent(in) :: klbound, kubound
        integer(INT8), intent(in) :: flags
        character(*), intent(in) :: name

        ! Local variables
        integer :: ierr
        integer(INT8), pointer :: msize(:)
        integer, parameter :: ndim = 4
        integer(INT8), parameter :: elsize = 8

        ! Executable code

        allocate(msize(ndim), stat=ierr)
        call allocate_error(ierr, name, mem_size)

        msize(1) = nubound-nlbound+1
        msize(2) = mubound-mlbound+1
        msize(3) = lubound-llbound+1
        msize(4) = kubound-klbound+1
        mem_size = mem_size + msize(1)*msize(2)*msize(3)*msize(4)*elsize

        allocate(dbleptr(nlbound:nubound, mlbound:mubound, &
                         llbound:lubound, klbound:kubound), stat=ierr)
        call allocate_error(ierr, name, mem_size)

        call C_MallocPlus_add_multi(mem_object%object, &
          c_loc(dbleptr(nlbound,mlbound,llbound,klbound)), ndim, &
          c_loc(msize(1)), elsize, trim(name)//char(0), flags)

        deallocate(msize)

      end subroutine FMallocPlus_malloc_4D_dble
!=======================================================================
      subroutine FMallocPlus_malloc_5D_dble(mem_object, dbleptr, &
          nlbound, nubound, mlbound, mubound, &
          llbound, lubound, klbound, kubound, &
          jlbound, jubound, &
          name, flags)

        use iso_c_binding, only : c_loc

        implicit none

        ! Argument definitions
        type(FMallocPlus_type), intent(inout) :: mem_object
        real(REAL8), dimension(:,:,:,:,:), pointer, intent(inout) :: dbleptr
        integer, intent(in) :: nlbound, nubound
        integer, intent(in) :: mlbound, mubound
        integer, intent(in) :: llbound, lubound
        integer, intent(in) :: klbound, kubound
        integer, intent(in) :: jlbound, jubound
        integer(INT8), intent(in) :: flags
        character(*), intent(in) :: name

        ! Local variables
        integer :: ierr
        integer(INT8), pointer :: msize(:)
        integer, parameter :: ndim = 5
        integer(INT8), parameter :: elsize = 8

        ! Executable code

        allocate(msize(ndim), stat=ierr)
        call allocate_error(ierr, name, mem_size)

        msize(1) = nubound-nlbound+1
        msize(2) = mubound-mlbound+1
        msize(3) = lubound-llbound+1
        msize(4) = kubound-klbound+1
        msize(5) = jubound-jlbound+1
        mem_size = mem_size + msize(1)*msize(2)*msize(3)*msize(4)*msize(5)*elsize

        allocate(dbleptr(nlbound:nubound, mlbound:mubound, &
                         llbound:lubound, klbound:kubound, &
                         jlbound:jubound), stat=ierr)
        call allocate_error(ierr, name, mem_size)

        call C_MallocPlus_add_multi(mem_object%object, &
          c_loc(dbleptr(nlbound,mlbound,llbound,klbound,jlbound)), ndim, &
          c_loc(msize(1)), elsize, trim(name)//char(0), flags)

        deallocate(msize)

      end subroutine FMallocPlus_malloc_5D_dble
!=======================================================================
      subroutine FMallocPlus_malloc_6D_dble(mem_object, dbleptr, &
          nlbound, nubound, mlbound, mubound, &
          llbound, lubound, klbound, kubound, &
          jlbound, jubound, ilbound, iubound, &
          name, flags)

        use iso_c_binding, only : c_loc

        implicit none

        ! Argument definitions
        type(FMallocPlus_type), intent(inout) :: mem_object
        real(REAL8), dimension(:,:,:,:,:,:), pointer, intent(inout) :: dbleptr
        integer, intent(in) :: nlbound, nubound
        integer, intent(in) :: mlbound, mubound
        integer, intent(in) :: llbound, lubound
        integer, intent(in) :: klbound, kubound
        integer, intent(in) :: jlbound, jubound
        integer, intent(in) :: ilbound, iubound
        integer(INT8), intent(in) :: flags
        character(*), intent(in) :: name

        ! Local variables
        integer :: ierr
        integer(INT8), pointer :: msize(:)
        integer, parameter :: ndim = 6
        integer(INT8), parameter :: elsize = 8

        ! Executable code

        allocate(msize(ndim), stat=ierr)
        call allocate_error(ierr, name, mem_size)

        msize(1) = nubound-nlbound+1
        msize(2) = mubound-mlbound+1
        msize(3) = lubound-llbound+1
        msize(4) = kubound-klbound+1
        msize(5) = jubound-jlbound+1
        msize(6) = iubound-ilbound+1
        mem_size = mem_size + msize(1)*msize(2)*msize(3)*msize(4)*msize(5)*msize(6)*elsize

        allocate(dbleptr(nlbound:nubound, mlbound:mubound, &
                         llbound:lubound, klbound:kubound, &
                         jlbound:jubound, ilbound:iubound), stat=ierr)
        call allocate_error(ierr, name, mem_size)

        call C_MallocPlus_add_multi(mem_object%object, &
          c_loc(dbleptr(nlbound,mlbound,llbound,klbound,jlbound,ilbound)), ndim, &
          c_loc(msize(1)), elsize, trim(name)//char(0), flags)

        deallocate(msize)

      end subroutine FMallocPlus_malloc_6D_dble
!=======================================================================
      subroutine allocate_error(ierr, message, mem_size)

         implicit none

#ifdef HAVE_MPI
#include "mpif.h"
#endif

         ! Argument definitions
         integer, intent(in) :: ierr
         character*(*), intent(in) :: message
         integer(INT8), intent(in) :: mem_size

         ! Local declarations

#ifdef HAVE_MPI
         integer :: ierror_reduce
#endif
         integer :: ierr_global, mpi_rank = 0

         ! Executable statements

#ifdef HAVE_MPI
         call MPI_COMM_RANK(MPI_COMM_WORLD, mpi_rank, ierr)

         call MPI_ALLREDUCE(ierr, ierr_global, 1, MPI_INTEGER, MPI_MAX, &
                            MPI_COMM_WORLD, ierror_reduce)
#else
         ierr_global = ierr
#endif

         if (ierr_global .eq. 0) return

         if (mem_size .lt. 10000000) then
            write (*,'(i3,":ERROR allocating ",a12,". Memsize is ",i6," KB") ' ) &
                      mpi_rank,message, mem_size/1000
         else
            write (*,'(i3,":ERROR allocating ",a12,". Memsize is ",i6," MB") ' ) &
                      mpi_rank,message, mem_size/1000000
         endif

#ifdef HAVE_MPI
         call mpi_finalize(ierr)
#endif
         stop

      end subroutine allocate_error

!=======================================================================

      end module FMallocPlus_mod

