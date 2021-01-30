!* Copyright 2015-2019.  Triad National Security, LLC. This material was produced
!* under U.S. Government contract 89233218CNA000001 for Los Alamos National 
!* Laboratory (LANL), which is operated by Triad National Security, LLC
!* for the U.S. Department of Energy. The U.S. Government has rights to use,
!* reproduce, and distribute this software.  NEITHER THE GOVERNMENT NOR
!* TRIAD NATIONAL SECURITY, LLC MAKES ANY WARRANTY, EXPRESS OR IMPLIED, OR
!* ASSUMES ANY LIABILITY FOR THE USE OF THIS SOFTWARE.  If software is modified
!* to produce derivative works, such modified software should be clearly marked,
!* so as not to confuse it with the version available from LANL.
!*
!* Licensed under the Apache License, Version 2.0 (the "License"); you may not
!* use this file except in compliance with the License. You may obtain a copy
!* of the License at
!*
!* http://www.apache.org/licenses/LICENSE-2.0
!*
!* Unless required by applicable law or agreed to in writing, software distributed
!* under the License is distributed on an "AS IS" BASIS, WITHOUT WARRANTIES OR
!* CONDITIONS OF ANY KIND, either express or implied. See the License for the
!* specific language governing permissions and limitations under the License.
!*
!* Under this license, it is required to include a reference to this work. We
!* request that each derivative work contain a reference to LANL Copyright
!* Disclosure C15076/LA-CC-15-054 so that this work's impact can be roughly
!* measured.
!*
!* This is LANL Copyright Disclosure C15076/LA-CC-15-054
!*/

!*
!*  PowerParser is a general purpose input file parser for software applications.
!*
!*  Authors: Chuck Wingate   XCP-2   caw@lanl.gov
!*           Robert Robey    XCP-2   brobey@lanl.gov
!*/

! ------------------------------------------------------------------------------
! This module contains a collection of subroutines and related data that
! is used to read input data from an input file.
! ------------------------------------------------------------------------------
  module FParser_module
 
    implicit none
    save
    private

#ifndef INT4_KIND_DIGITS
#define INT4_KIND_DIGITS 6
#endif

#ifndef INT8_KIND_DIGITS
#define INT8_KIND_DIGITS 16
#endif

#ifndef REAL4_KIND_DIGITS
#define REAL4_KIND_DIGITS 6
#endif

#ifndef REAL8_KIND_DIGITS
#define REAL8_KIND_DIGITS 12
#endif

    integer, parameter :: INT4   = SELECTED_INT_KIND(INT4_KIND_DIGITS)
    integer, parameter :: INT32  = SELECTED_INT_KIND(INT4_KIND_DIGITS)

    integer, parameter :: INT8   = SELECTED_INT_KIND(INT8_KIND_DIGITS)
    integer, parameter :: INT64  = SELECTED_INT_KIND(INT8_KIND_DIGITS)

    integer, parameter :: REAL4  = SELECTED_REAL_KIND(REAL4_KIND_DIGITS)
    integer, parameter :: REAL32 = SELECTED_REAL_KIND(REAL4_KIND_DIGITS)

    integer, parameter :: REAL8  = SELECTED_REAL_KIND(REAL8_KIND_DIGITS)
    integer, parameter :: REAL64 = SELECTED_REAL_KIND(REAL8_KIND_DIGITS)

    integer :: mype
    integer :: numpe
    integer :: iope

    integer :: FParser_allostat

    public QueryFParser               ! Interface to get command values.
    public FParser_size               ! Interface to get sizes before getting values.
    public FParser_sizeb              ! Interface to get all sizes before getting values.
 
    public FParser_initialize         ! Initialize the parser and read in the input deck.
    public FParser_dictionary_add     ! Add a constant to the Parser dictionary
    public FParser_dictionary_env_add ! Add an environment variable constant to the Parser dictionary
    public FParser_compile_buffer     ! Compile the input buffer
    public FParser_initialized_check  ! Interface to check if FParser has been initialized
    public FParser_check_processed    ! Check processed flags.
    public FParser_disable_check      ! Disable check processed.
    public FParser_terminate          ! Last FParser call, check processed flags, clean.
    public FParser_cmd_in_input       ! Check for command in the user input.
    public FParser_cmd_set_processed  ! Set the command as being processed.
    public FParser_echo_user_input    ! Echo user input to a file.
    public FParser_count_var          ! Count numer of internal input variables
    public FParser_echo_fvf           ! Echo functions, variables, final buffer
    public FParser_get_num_include_files ! The number of include files that have been processed
    public FParser_get_include_file    ! Retrieve the names of any input files included in the main input file
 
    private FParser_list_funcs        ! List functions to file.
    private FParser_list_vars         ! List variables to file.
    private FParser_echo_fbuffer      ! Echo final buffer to a file.
 
    ! Define an interface block that allows several subroutines to be called
    ! using a common name, i.e. QueryFParser. This is for getting values.
    interface QueryFParser
       module procedure QueryFParser_c_00   ! Single character
       module procedure QueryFParser_c_0    ! Character string
       module procedure QueryFParser_c_1    ! 1d array of character strings
       module procedure QueryFParser_c_2    ! 2d array of character strings
       module procedure QueryFParser_c_3    ! 3d array of character strings
       module procedure QueryFParser_c_4    ! 4d array of character strings
 
       module procedure QueryFParser_i_0    ! Single integer
       module procedure QueryFParser_i_1    ! 1d integer array
       module procedure QueryFParser_i_2    ! 2d integer array
       module procedure QueryFParser_i_3    ! 3d integer array
       module procedure QueryFParser_i_4    ! 4d integer array
 
       module procedure QueryFParser_i8_0   ! Single long int (int64_t)
       module procedure QueryFParser_i8_1   ! Single long int (int64_t)
 
       module procedure QueryFParser_l_0    ! Single logical
       module procedure QueryFParser_l_1    ! 1d logical array
       module procedure QueryFParser_l_2    ! 2d logical array
       module procedure QueryFParser_l_3    ! 3d logical array
       module procedure QueryFParser_l_4    ! 4d logical array
 
       module procedure QueryFParser_r8_0   ! Single real
       module procedure QueryFParser_r8_1   ! 1d real array
       module procedure QueryFParser_r8_2   ! 2d real array
       module procedure QueryFParser_r8_3   ! 3d real array
       module procedure QueryFParser_r8_4   ! 4d real array
    end interface
 
    ! Define an interface block that allows several subroutines to be called
    ! using a common name, i.e. FParser_size. This is for getting sizes.
    ! For variables that do not exist, a size of 0 must be returned.
    ! This size is used as a check and passed into other routines.
    ! So returning a negative number for the size of non-existant vars
    ! does cause problems.
    interface FParser_size
       module procedure FParser_size_1    ! Single dimensional array
       module procedure FParser_size_2    ! 2d array
       module procedure FParser_size_3    ! 3d array
       module procedure FParser_size_4    ! 4d array
       module procedure FParser_size_5    ! 5d array
       module procedure FParser_size_6    ! 6d array
    end interface
    ! Define an interface block that allows several subroutines to be called
    ! using a common name, i.e. FParser_sizeb. This is for getting all sizes.
    interface FParser_sizeb
       module procedure FParser_sizeb_2   ! 2d array
    end interface
 
    ! restart_block subroutines.
    public FParser_rb_check           ! Do the restart blocks check
    public FParser_echo_rb_info       ! Echo restart block info for one block
 
 
    ! When...then subroutines.
    public FParser_whenthen_num       ! Get number of when...thens
    public FParser_whenthen_check     ! Do the whenthen check
    public FParser_whenthen_setcfp
    public FParser_whenthen_reset     ! Reset cmds final buffer pointer.
    public FParser_whenthen_casize
    public FParser_whenthen_ca
    public FParser_whenthen_satsize
    public FParser_whenthen_getsat
    public FParser_whenthen_setsat
    public FParser_whenthen_getprocessed
    public FParser_whenthen_setprocessed
    public FParser_whenthen_getseq
    public FParser_whenthen_setseq
 
    ! Private module data.
 
 
    logical :: FParser_initialized = .false.
    logical :: FParser_do_check    = .true.
 
! ==============================================================================
  contains
! ------------------------------------------------------------------------------
! ******************************************************************************
! ******************************************************************************
! Get logical values
! ******************************************************************************
! ******************************************************************************
 
    ! ==========================================================================
    ! Get a single logical value.
    ! ==========================================================================
    subroutine QueryFParser_l_0(cmdname, cmdvalue, noskip)
 
      character(*), intent(in)    :: cmdname
      logical,      intent(inout) :: cmdvalue
      logical,      intent(in)    :: noskip
 
      character :: cmd_array(len(cmdname))
      integer :: ival, iskip, c
 
      if(.not.FParser_initialized)return
 
      ! Pass an array of single chars rather than a character string.
      do c = 1,len(cmdname)
         cmd_array(c) = cmdname(c:c)
      enddo
 
      ! Default ival. Note that cmdvalue only changes if the command is found
      ! in the input file.
      ival = 0
      if(cmdvalue) ival = 1
 
      ! Cannot pass logicals to C++ as booleans, so use integer.
      iskip = 1
      if (noskip) iskip = 0
 
      ! This is a C++ funtion. Note that we have to use integers because fortran
      ! logical and C++ bool are not compatible.
      call get_logical0(cmd_array, ival, len(cmdname), iskip)
 
      if (noskip) then
         cmdvalue = .false.
         if (ival .eq. 1) cmdvalue = .true.
      endif
 
    end subroutine QueryFParser_l_0
 
 
    ! ==========================================================================
    ! Get a one dimensional logical array.
    ! ==========================================================================
    subroutine QueryFParser_l_1(cmdname, array, noskip, array_size)
 
      character(*), intent(in)    :: cmdname
      integer,      intent(in)    :: array_size
      logical,      intent(inout) :: array(array_size)
      logical,      intent(in)    :: noskip
 
      character :: cmd_array(len(cmdname))
      integer :: ival_array(array_size), i, iskip, c
 
      if(.not.FParser_initialized)return
 
      ! Pass an array of single chars rather than a character string.
      do c = 1,len(cmdname)
         cmd_array(c) = cmdname(c:c)
      enddo
 
      ! Default ival_array. We have to do this because the array values only
      ! change if the command is found in the input file and the array as
      ! input already has default values.
      ival_array = 0;
      if (noskip) then
         do i=1,array_size
            ival_array(i) = 0
            if (array(i)) ival_array(i) = 1
         enddo
      endif
 
      ! Cannot pass logicals to C++ as booleans, so use integer.
      iskip = 1
      if (noskip) iskip = 0
 
      ! This is a C++ funtion. Note that we have to use integers because fortran
      ! logical and C++ bool are not compatible.
      call get_logical1(cmd_array, ival_array, array_size, len(cmdname), iskip)
 
      if (noskip) then
         do i=1,array_size
            array(i) = .false.
            if (ival_array(i) .eq. 1) array(i) = .true.
         enddo
      endif
 
    end subroutine QueryFParser_l_1
 
 
    ! ==========================================================================
    ! Get two dimensional logical array values.
    ! ==========================================================================
    subroutine QueryFParser_l_2(cmdname, array, noskip, size1, size2)
 
      character(*), intent(in)    :: cmdname
      integer,      intent(in)    :: size1, size2
      logical,      intent(inout) :: array(size1, size2)
      logical,      intent(in)    :: noskip
 
      character :: cmd_array(len(cmdname))
      integer :: val_array(size1*size2)
      integer :: i, j, i1d, iskip, c
 
      if(.not.FParser_initialized)return
 
      ! Pass an array of single chars rather than a character string.
      do c = 1,len(cmdname)
         cmd_array(c) = cmdname(c:c)
      enddo
 
      ! Default val_array. We have to do this because the array values only
      ! change if the command is found in the input file and the array as
      ! input already has default values.
      ! We pass a 1d array to C++ because it is easier.
      val_array = 0;
      if (noskip) then
         i1d = 1
         do j = 1, size2
            do i = 1, size1
               !i1d = i + (j-1)*size1
               val_array(i1d) = 0
               if(array(i,j)) val_array(i1d) = 1
               i1d = i1d + 1
            enddo
         enddo
      endif
 
      ! Cannot pass logicals to C++ as booleans, so use integer.
      iskip = 1
      if (noskip) iskip = 0
 
      ! This is a C++ funtion. Note that we have to use integers because fortran
      ! logical and C++ bool are not compatible.
      call get_logical2(cmd_array, val_array, size1, size2, len(cmdname), iskip)
 
      ! Copy values back in to the 2d array.
      if (noskip) then
         i1d = 1
         do j = 1, size2
            do i = 1, size1
               !i1d = i + (j-1)*size1
               array(i,j) = .false.
               if (val_array(i1d) .eq. 1) array(i,j) = .true.
               i1d = i1d + 1
            enddo
         enddo
      endif
 
    end subroutine QueryFParser_l_2
 
 
    ! ==========================================================================
    ! Get three dimensional logical array values.
    ! ==========================================================================
    subroutine QueryFParser_l_3(cmdname, array, noskip, size1, size2, size3)
 
      character(*), intent(in)    :: cmdname
      integer,      intent(in)    :: size1, size2, size3
      logical,      intent(inout) :: array(size1, size2, size3)
      logical,      intent(in)    :: noskip
 
      character :: cmd_array(len(cmdname))
      integer :: val_array(size1*size2*size3)
      integer :: i, j, k, i1d, iskip, c
 
      if(.not.FParser_initialized)return
 
      ! Pass an array of single chars rather than a character string.
      do c = 1,len(cmdname)
         cmd_array(c) = cmdname(c:c)
      enddo
 
      ! Default val_array. We have to do this because the array values only
      ! change if the command is found in the input file and the array as
      ! input already has default values.
      ! We pass a 1d array to C++ because it is easier.
      val_array = 0;
      if (noskip) then
         i1d = 1
         do k = 1, size3
            do j = 1, size2
               do i = 1, size1
                  !i1d = i + (j-1)*size1 + (k-1)*size1*size2
                  val_array(i1d) = 0
                  if (array(i,j,k)) val_array(i1d) = 1
                  i1d = i1d + 1
               enddo
            enddo
         enddo
      endif
 
      ! Cannot pass logicals to C++ as booleans, so use integer.
      iskip = 1
      if (noskip) iskip = 0
 
      ! This is a C++ funtion. Note that we have to use integers because fortran
      ! logical and C++ bool are not compatible.
      call get_logical3(cmd_array, val_array, size1, size2, size3, len(cmdname), &
                        iskip)
 
      ! Copy values back in to the 3d array.
      if (noskip) then
         i1d = 1
         do k = 1, size3
            do j = 1, size2
               do i = 1, size1
                  !i1d = i + (j-1)*size1 + (k-1)*size1*size2
                  array(i,j,k) = .false.
                  if (val_array(i1d) .eq. 1) array(i,j,k) = .true.
                  i1d = i1d + 1
               enddo
            enddo
         enddo
      endif
 
    end subroutine QueryFParser_l_3
 
 
    ! ==========================================================================
    ! Get four dimensional logical array values.
    ! ==========================================================================
    subroutine QueryFParser_l_4(cmdname, array, noskip, size1, size2, size3, size4)
 
      character(*), intent(in)    :: cmdname
      integer,      intent(in)    :: size1, size2, size3, size4
      logical,      intent(inout) :: array(size1, size2, size3, size4)
      logical,      intent(in)    :: noskip
 
      character :: cmd_array(len(cmdname))
      integer :: val_array(size1*size2*size3*size4)
      integer :: i, j, k, l, i1d, iskip, c
 
      if(.not.FParser_initialized)return
 
      ! Pass an array of single chars rather than a character string.
      do c = 1,len(cmdname)
         cmd_array(c) = cmdname(c:c)
      enddo
 
      ! Default val_array. We have to do this because the array values only
      ! change if the command is found in the input file and the array as
      ! input already has default values.
      ! We pass a 1d array to C++ because it is easier.
      val_array = 0;
      if (noskip) then
         i1d = 1
         do l = 1, size4
            do k = 1, size3
               do j = 1, size2
                  do i = 1, size1
                     !i1d = i + (j-1)*size1 + (k-1)*size1*size2 + &
                     !     (l-1)*size1*size2*size3
                     val_array(i1d) = 0
                     if (array(i,j,k,l)) val_array(i1d) = 1
                     i1d = i1d + 1
                  enddo
               enddo
            enddo
         enddo
      endif
 
      ! Cannot pass logicals to C++ as booleans, so use integer.
      iskip = 1
      if (noskip) iskip = 0
 
      ! This is a C++ funtion. Note that we have to use integers because fortran
      ! logical and C++ bool are not compatible.
      call get_logical4(cmd_array, val_array, size1, size2, size3, size4, &
                        len(cmdname), iskip)
 
      ! Copy values back in to the 4d array.
      if (noskip) then
         i1d = 1
         do l = 1, size4
            do k = 1, size3
               do j = 1, size2
                  do i = 1, size1
                     !i1d = i + (j-1)*size1 + (k-1)*size1*size2 + &
                     !     (l-1)*size1*size2*size3
                     array(i,j,k,l) = .false.
                     if (val_array(i1d) .eq. 1) array(i,j,k,l) = .true.
                     i1d = i1d + 1
                  enddo
               enddo
            enddo
         enddo
      endif
 
    end subroutine QueryFParser_l_4
 
 
 
 
! ******************************************************************************
! ******************************************************************************
! Get integer values
! ******************************************************************************
! ******************************************************************************
 
    ! ==========================================================================
    ! Get a single integer value.
    ! ==========================================================================
    subroutine QueryFParser_i_0(cmdname, cmdvalue, noskip)
 
      character(*), intent(in)    :: cmdname
      integer,      intent(inout) :: cmdvalue
      logical,      intent(in)    :: noskip
 
      character :: cmd_array(len(cmdname))
      integer :: iskip, c
 
      if(.not.FParser_initialized)return
 
      ! Pass an array of single chars rather than a character string.
      do c = 1,len(cmdname)
         cmd_array(c) = cmdname(c:c)
      enddo
 
      ! Cannot pass logicals to C++ as booleans, so use integer.
      iskip = 1
      if (noskip) iskip = 0
 
      ! This is a C++ funtion. Note that ival (and thus cmdvalue) only changes
      ! if the command is found in the input file.
      call get_integer0(cmd_array, cmdvalue, len(cmdname), iskip)
 
    end subroutine QueryFParser_i_0
 
 
    ! ==========================================================================
    ! Get a one dimensional integer array.
    ! ==========================================================================
    subroutine QueryFParser_i_1(cmdname, array, noskip, array_size)
 
      character(*), intent(in)    :: cmdname
      integer,      intent(in)    :: array_size
      integer,      intent(inout) :: array(array_size)
      logical,      intent(in)    :: noskip
 
      character :: cmd_array(len(cmdname))
      integer :: iskip, c
 
      if(.not.FParser_initialized)return
 
      ! Pass an array of single chars rather than a character string.
      do c = 1,len(cmdname)
         cmd_array(c) = cmdname(c:c)
      enddo
 
      ! Cannot pass logicals to C++ as booleans, so use integer.
      iskip = 1
      if (noskip) iskip = 0
 
      ! This is a C++ funtion.
      call get_integer1(cmd_array, array, array_size, len(cmdname), iskip)
 
    end subroutine QueryFParser_i_1
 
 
    ! ==========================================================================
    ! Get two dimensional integer array values.
    ! ==========================================================================
    subroutine QueryFParser_i_2(cmdname, array, noskip, size1, size2)
 
      character(*), intent(in)    :: cmdname
      integer,      intent(in)    :: size1, size2
      integer,      intent(inout) :: array(size1, size2)
      logical,      intent(in)    :: noskip
 
      character :: cmd_array(len(cmdname))
      integer :: iskip, c
 
      if(.not.FParser_initialized)return
 
      ! Pass an array of single chars rather than a character string.
      do c = 1,len(cmdname)
         cmd_array(c) = cmdname(c:c)
      enddo
 
      ! Cannot pass logicals to C++ as booleans, so use integer.
      iskip = 1
      if (noskip) iskip = 0
 
      ! This is a C++ funtion.
      call get_integer2(cmd_array, array, size1, size2, len(cmdname), iskip)
 
    end subroutine QueryFParser_i_2
 
 
    ! ==========================================================================
    ! Get three dimensional integer array values.
    ! ==========================================================================
    subroutine QueryFParser_i_3(cmdname, array, noskip, size1, size2, size3)
 
      character(*), intent(in)    :: cmdname
      integer,      intent(in)    :: size1, size2, size3
      integer,      intent(inout) :: array(size1, size2, size3)
      logical,      intent(in)    :: noskip
 
      character :: cmd_array(len(cmdname))
      integer :: iskip, c
 
      if(.not.FParser_initialized)return
 
      ! Pass an array of single chars rather than a character string.
      do c = 1,len(cmdname)
         cmd_array(c) = cmdname(c:c)
      enddo
 
      ! Cannot pass logicals to C++ as booleans, so use integer.
      iskip = 1
      if (noskip) iskip = 0
 
      ! This is a C++ funtion.
      call get_integer3(cmd_array, array, size1, size2, size3, len(cmdname), &
                        iskip)
 
    end subroutine QueryFParser_i_3
 
 
    ! ==========================================================================
    ! Get four dimensional integer array values.
    ! ==========================================================================
    subroutine QueryFParser_i_4(cmdname, array, noskip, size1, size2, size3, size4)
 
      character(*), intent(in)    :: cmdname
      integer,      intent(in)    :: size1, size2, size3, size4
      integer,      intent(inout) :: array(size1, size2, size3, size4)
      logical,      intent(in)    :: noskip
 
      character :: cmd_array(len(cmdname))
      integer :: iskip, c
 
      if(.not.FParser_initialized)return
 
      ! Pass an array of single chars rather than a character string.
      do c = 1,len(cmdname)
         cmd_array(c) = cmdname(c:c)
      enddo
 
      ! Cannot pass logicals to C++ as booleans, so use integer.
      iskip = 1
      if (noskip) iskip = 0
 
      ! This is a C++ funtion.
      call get_integer4(cmd_array, array, size1, size2, size3, size4, &
                        len(cmdname), iskip)
 
    end subroutine QueryFParser_i_4
 
    ! ==========================================================================
    ! Get a single long integer (int64_t) value.
    ! ==========================================================================
    subroutine QueryFParser_i8_0(cmdname, cmdvalue, noskip)
 
      character(*),  intent(in)    :: cmdname
      integer(INT8), intent(inout) :: cmdvalue
      logical,       intent(in)    :: noskip
 
      character :: cmd_array(len(cmdname))
      integer :: iskip, c
 
      if(.not.FParser_initialized)return
 
      ! Pass an array of single chars rather than a character string.
      do c = 1,len(cmdname)
         cmd_array(c) = cmdname(c:c)
      enddo
 
      ! Cannot pass logicals to C++ as booleans, so use integer.
      iskip = 1
      if (noskip) iskip = 0
 
      ! This is a C++ funtion. Note that ival (and thus cmdvalue) only changes
      ! if the command is found in the input file.
      call get_int8_0(cmd_array, cmdvalue, len(cmdname), iskip)
 
    end subroutine QueryFParser_i8_0

    ! ==========================================================================
    ! Get a one dimensional integer (int64_t) array.
    ! ==========================================================================
    subroutine QueryFParser_i8_1(cmdname, array, noskip, array_size)
 
      character(*),  intent(in)    :: cmdname
      integer,       intent(in)    :: array_size
      integer(INT8), intent(inout) :: array(array_size)
      logical,       intent(in)    :: noskip
 
      character :: cmd_array(len(cmdname))
      integer :: iskip, c
 
      if(.not.FParser_initialized)return
 
      ! Pass an array of single chars rather than a character string.
      do c = 1,len(cmdname)
         cmd_array(c) = cmdname(c:c)
      enddo
 
      ! Cannot pass logicals to C++ as booleans, so use integer.
      iskip = 1
      if (noskip) iskip = 0
 
      ! This is a C++ funtion.
      call get_int8_1(cmd_array, array, array_size, len(cmdname), iskip)
 
    end subroutine QueryFParser_i8_1
 
 
 
 
! ******************************************************************************
! ******************************************************************************
! Get real values
! ******************************************************************************
! ******************************************************************************
 
    ! ==========================================================================
    ! Get a single real value.
    ! ==========================================================================
    subroutine QueryFParser_r8_0(cmdname, cmdvalue, noskip)
 
      character(*), intent(in)    :: cmdname
      real(REAL8),  intent(inout) :: cmdvalue
      logical,      intent(in)    :: noskip
 
      character :: cmd_array(len(cmdname))
      integer :: iskip, c
 
      if(.not.FParser_initialized)return
 
      ! Pass an array of single chars rather than a character string.
      do c = 1,len(cmdname)
         cmd_array(c) = cmdname(c:c)
      enddo
 
      ! Cannot pass logicals to C++ as booleans, so use integer.
      iskip = 1
      if (noskip) iskip = 0
 
      ! This is a C++ funtion. Note that rval (and thus cmdvalue) only changes
      ! if the command is found in the input file.
      call get_real0(cmd_array, cmdvalue, len(cmdname), iskip)
 
    end subroutine QueryFParser_r8_0
 
 
    ! ==========================================================================
    ! Get one dimensional real array values.
    ! ==========================================================================
    subroutine QueryFParser_r8_1(cmdname, array, noskip, array_size)
 
      character(*), intent(in)    :: cmdname
      integer,      intent(in)    :: array_size
      real(REAL8),  intent(inout) :: array(array_size)
      logical,      intent(in)    :: noskip
 
      character :: cmd_array(len(cmdname))
      integer :: iskip, c
 
      if(.not.FParser_initialized)return
 
      ! Pass an array of single chars rather than a character string.
      do c = 1,len(cmdname)
         cmd_array(c) = cmdname(c:c)
      enddo
 
      ! Cannot pass logicals to C++ as booleans, so use integer.
      iskip = 1
      if (noskip) iskip = 0
 
      ! This is a C++ funtion.
      call get_real1(cmd_array, array, array_size, len(cmdname), iskip)
 
    end subroutine QueryFParser_r8_1
 
 
    ! ==========================================================================
    ! Get two dimensional real array values.
    ! ==========================================================================
    subroutine QueryFParser_r8_2(cmdname, array, noskip, size1, size2)
 
      character(*), intent(in)    :: cmdname
      integer,      intent(in)    :: size1, size2
      real(REAL8),  intent(inout) :: array(size1, size2)
      logical,      intent(in)    :: noskip
 
      character :: cmd_array(len(cmdname))
      integer :: iskip, c
 
      if(.not.FParser_initialized)return
 
      ! Pass an array of single chars rather than a character string.
      do c = 1,len(cmdname)
         cmd_array(c) = cmdname(c:c)
      enddo
 
      ! Cannot pass logicals to C++ as booleans, so use integer.
      iskip = 1
      if (noskip) iskip = 0
 
      ! This is a C++ funtion.
      call get_real2(cmd_array, array, size1, size2, len(cmdname), &
                     iskip)
 
    end subroutine QueryFParser_r8_2
 
 
    ! ==========================================================================
    ! Get three dimensional real array values.
    ! ==========================================================================
    subroutine QueryFParser_r8_3(cmdname, array, noskip, size1, size2, size3)
 
      character(*), intent(in)    :: cmdname
      integer,      intent(in)    :: size1, size2, size3
      real(REAL8),  intent(inout) :: array(size1, size2, size3)
      logical,      intent(in)    :: noskip
 
      character :: cmd_array(len(cmdname))
      integer :: iskip, c
 
      if(.not.FParser_initialized)return
 
      ! Pass an array of single chars rather than a character string.
      do c = 1,len(cmdname)
         cmd_array(c) = cmdname(c:c)
      enddo
 
      ! Cannot pass logicals to C++ as booleans, so use integer.
      iskip = 1
      if (noskip) iskip = 0
 
      ! This is a C++ funtion.
      call get_real3(cmd_array, array, size1, size2, size3, len(cmdname), &
                     iskip)
 
    end subroutine QueryFParser_r8_3
 
 
    ! ==========================================================================
    ! Get four dimensional real array values.
    ! ==========================================================================
    subroutine QueryFParser_r8_4(cmdname, array, noskip, size1, size2, size3, size4)
 
      character(*), intent(in)    :: cmdname
      integer,      intent(in)    :: size1, size2, size3, size4
      real(REAL8),  intent(inout) :: array(size1, size2, size3, size4)
      logical,      intent(in)    :: noskip
 
      character :: cmd_array(len(cmdname))
      integer :: iskip, c
 
      if(.not.FParser_initialized)return
 
      ! Pass an array of single chars rather than a character string.
      do c = 1,len(cmdname)
         cmd_array(c) = cmdname(c:c)
      enddo
 
      ! Cannot pass logicals to C++ as booleans, so use integer.
      iskip = 1
      if (noskip) iskip = 0
 
      ! This is a C++ funtion.
      call get_real4(cmd_array, array, size1, size2, size3, size4, &
                     len(cmdname), iskip)
 
    end subroutine QueryFParser_r8_4
 
 
! ******************************************************************************
! ******************************************************************************
! Get character values
! ******************************************************************************
! ******************************************************************************
 
    ! ==========================================================================
    ! Get a single character value.
    ! ==========================================================================
    subroutine QueryFParser_c_00(cmdname, cmdvalue, noskip)
 
      character(*), intent(in)    :: cmdname
      character(*), intent(inout) :: cmdvalue
      logical,      intent(in)    :: noskip
 
      character :: cmd_array(len(cmdname))
      integer :: iskip, c
 
      if(.not.FParser_initialized)return
 
      ! Pass an array of single chars rather than a character string.
      do c = 1,len(cmdname)
         cmd_array(c) = cmdname(c:c)
      enddo
 
      ! Cannot pass logicals to C++ as booleans, so use integer.
      iskip = 1
      if (noskip) iskip = 0
 
      ! This is a C++ funtion. Note that cval (and thus cmdvalue) only changes
      ! if the command is found in the input file.
      call get_char00(cmd_array, cmdvalue(1:1), len(cmdname), iskip)
 
    end subroutine QueryFParser_c_00
 
 
    ! ==========================================================================
    ! Get a character string.
    ! ==========================================================================
    subroutine QueryFParser_c_0(cmdname, cmdvalue, noskip, nchar, nchar_actual)
 
      character(*),     intent(in)            :: cmdname
      integer,          intent(in)            :: nchar
      character(nchar), intent(inout)         :: cmdvalue
      logical,          intent(in)            :: noskip
      integer,          intent(out), optional :: nchar_actual ! the actual max len
 
      character :: cmd_array(len(cmdname)), val_array(nchar)
      integer :: iskip, c
      integer :: nchar_actual_val
 
      if(.not.FParser_initialized)return
 
      ! Pass an array of single chars rather than a character string.
      do c = 1,len(cmdname)
         cmd_array(c) = cmdname(c:c)
      enddo
      do c = 1,nchar
         val_array(c) = cmdvalue(c:c)
      enddo
 
      ! Cannot pass logicals to C++ as booleans, so use integer.
      iskip = 1
      if (noskip) iskip = 0
 
      ! This is a C++ funtion. Note that cval (and thus cmdvalue) only changes
      ! if the command is found in the input file.
      call get_char0(cmd_array, val_array, len(cmdname), nchar, iskip, nchar_actual_val)
      if( present(nchar_actual) ) then
         nchar_actual = nchar_actual_val
      endif
 
      ! Unpack the results into the output array of strings.
      if (noskip) then
         do c = 1,nchar
            cmdvalue(c:c) = val_array(c)
         enddo
      endif
 
    end subroutine QueryFParser_c_0
 
 
    ! ==========================================================================
    ! Get one dimensional array of character strings.
    ! ==========================================================================
    subroutine QueryFParser_c_1(cmdname, array, noskip, nchar, size1)
 
      character(*),     intent(in)    :: cmdname
      integer,          intent(in)    :: size1
      integer,          intent(in)    :: nchar
      character(nchar), intent(inout) :: array(size1)
      logical,          intent(in)    :: noskip
 
      character :: cmd_array(len(cmdname)), val_array(nchar*size1)
      integer :: i, c, i1d, iskip
 
      if(.not.FParser_initialized)return
 
      ! Pass an array of single chars rather than a character string.
      do c = 1,len(cmdname)
         cmd_array(c) = cmdname(c:c)
      enddo
 
      ! We pack all the fortran strings into one long character array, this
      ! makes it easier to send to C++.
      val_array = ' '
      if (noskip) then
         i1d = 1
         do i = 1,size1
            do c = 1,nchar
               val_array(i1d) = array(i)(c:c)
               i1d = i1d+1
            enddo
         enddo
      endif
 
      ! Cannot pass logicals to C++ as booleans, so use integer.
      iskip = 1
      if (noskip) iskip = 0
 
      ! This is a C++ funtion.
      call get_char1(cmd_array, val_array, nchar, size1, len(cmdname), iskip)
 
      ! Unpack the results into the output array of strings.
      if (noskip) then
         i1d = 1
         do i = 1,size1
            do c = 1,nchar
               array(i)(c:c) = val_array(i1d)
               i1d = i1d+1
            enddo
         enddo
      endif
 
    end subroutine QueryFParser_c_1
 
 
    ! ==========================================================================
    ! Get two dimensional array of character strings.
    ! ==========================================================================
    subroutine QueryFParser_c_2(cmdname, array, noskip, nchar, size1, size2)
 
      character(*),     intent(in)    :: cmdname
      integer,          intent(in)    :: size1, size2
      integer,          intent(in)    :: nchar
      character(nchar), intent(inout) :: array(size1, size2)
      logical,          intent(in)    :: noskip
 
      character :: cmd_array(len(cmdname)), val_array(nchar*size1*size2)
      integer :: i, j, c, i1d, iskip
 
      if(.not.FParser_initialized)return
 
      ! Pass an array of single chars rather than a character string.
      do c = 1,len(cmdname)
         cmd_array(c) = cmdname(c:c)
      enddo
 
      ! We pack all the fortran strings into one long character array, this
      ! makes it easier to send to C++.
      val_array = ' '
      if (noskip) then
         i1d = 1
         do j = 1,size2
            do i = 1,size1
               do c = 1,nchar
                  val_array(i1d) = array(i,j)(c:c)
                  i1d = i1d + 1
               enddo
            enddo
         enddo
      endif
 
      ! Cannot pass logicals to C++ as booleans, so use integer.
      iskip = 1
      if (noskip) iskip = 0
 
      ! This is a C++ funtion.
      call get_char2(cmd_array, val_array, nchar, size1, size2, len(cmdname), &
                     iskip)
 
      ! Unpack the results into the output array of strings.
      if (noskip) then
         i1d = 1
         do j = 1,size2
            do i = 1,size1
               do c = 1,nchar
                  array(i,j)(c:c) = val_array(i1d)
                  i1d = i1d+1
               enddo
            enddo
         enddo
      endif
 
    end subroutine QueryFParser_c_2
 
 
    ! ==========================================================================
    ! Get three dimensional array of character strings.
    ! ==========================================================================
    subroutine QueryFParser_c_3(cmdname, array, noskip, nchar, size1, size2, size3)
 
      character(*),     intent(in)    :: cmdname
      integer,          intent(in)    :: size1, size2, size3
      integer,          intent(in)    :: nchar
      character(nchar), intent(inout) :: array(size1, size2, size3)
      logical,          intent(in)    :: noskip
 
      character :: cmd_array(len(cmdname)), val_array(nchar*size1*size2*size3)
      integer :: i, j, k, c, i1d, iskip
 
      if(.not.FParser_initialized)return
 
      ! Pass an array of single chars rather than a character string.
      do c = 1,len(cmdname)
         cmd_array(c) = cmdname(c:c)
      enddo
 
      ! We pack all the fortran strings into one long character array, this
      ! makes it easier to send to C++.
      val_array = ' '
      if (noskip) then
         i1d = 1
         do k = 1,size3
            do j = 1,size2
               do i = 1,size1
                  do c = 1,nchar
                     val_array(i1d) = array(i,j,k)(c:c)
                     i1d = i1d + 1
                  enddo
               enddo
            enddo
         enddo
      endif
 
      ! Cannot pass logicals to C++ as booleans, so use integer.
      iskip = 1
      if (noskip) iskip = 0
 
      ! This is a C++ funtion.
      call get_char3(cmd_array, val_array, nchar, size1, size2, size3, &
                     len(cmdname), iskip)
 
      ! Unpack the results into the output array of strings.
      if (noskip) then
         i1d = 1
         do k = 1,size3
            do j = 1,size2
               do i = 1,size1
                  do c = 1,nchar
                     array(i,j,k)(c:c) = val_array(i1d)
                     i1d = i1d+1
                  enddo
               enddo
            enddo
         enddo
      endif
 
    end subroutine QueryFParser_c_3
 
 
    ! ==========================================================================
    ! Get four dimensional array of character strings.
    ! ==========================================================================
    subroutine QueryFParser_c_4(cmdname, array, noskip, nchar, size1, size2, size3, &
                          size4)
 
      character(*),     intent(in)    :: cmdname
      integer,          intent(in)    :: size1, size2, size3, size4
      integer,          intent(in)    :: nchar
      character(nchar), intent(inout) :: array(size1, size2, size3, size4)
      logical,          intent(in)    :: noskip
 
      character :: cmd_array(len(cmdname))
      character :: val_array(nchar*size1*size2*size3*size4)
      integer :: i,j,k,l,c,i1d, iskip
 
      if(.not.FParser_initialized)return
 
      ! Pass an array of single chars rather than a character string.
      do c = 1,len(cmdname)
         cmd_array(c) = cmdname(c:c)
      enddo
 
      ! We pack all the fortran strings into one long character array, this
      ! makes it easier to send to C++.
      val_array = ' '
      if (noskip) then
         i1d = 1
         do l = 1,size4
            do k = 1,size3
               do j = 1,size2
                  do i = 1,size1
                     do c = 1,nchar
                        val_array(i1d) = array(i,j,k,l)(c:c)
                        i1d = i1d + 1
                     enddo
                  enddo
               enddo
            enddo
         enddo
      endif
 
      ! Cannot pass logicals to C++ as booleans, so use integer.
      iskip = 1
      if (noskip) iskip = 0
 
      ! This is a C++ funtion.
      call get_char4(cmd_array, val_array, nchar, size1, size2, size3, size4, &
                     len(cmdname), iskip)
 
      ! Unpack the results into the output array of strings.
      if (noskip) then
         i1d = 1
         do l = 1,size4
            do k = 1,size3
               do j = 1,size2
                  do i = 1,size1
                     do c = 1,nchar
                        array(i,j,k,l)(c:c) = val_array(i1d)
                        i1d = i1d+1
                     enddo
                  enddo
               enddo
            enddo
         enddo
      endif
 
    end subroutine QueryFParser_c_4
 
 
 
 
 
! ******************************************************************************
! ******************************************************************************
! Get sizes
! ******************************************************************************
! ******************************************************************************
 
    ! ==========================================================================
    ! Get input sizes before actually reading values - 1d array.
    ! ==========================================================================
    subroutine FParser_size_1(cmdname, array_size)
 
      character(*), intent(in)    :: cmdname
      integer,      intent(inout) :: array_size
 
      character :: cmd_array(len(cmdname))
      integer :: c
 
      ! Important to initalize since array_size is accumulated.
      array_size = 0
 
      if(.not.FParser_initialized)return
 
      ! Pass an array of single chars rather than a character string.
      do c = 1,len(cmdname)
         cmd_array(c) = cmdname(c:c)
      enddo
 
      ! One dimensional arrays.
      call get_size1(cmd_array, array_size, len(cmdname))
 
    end subroutine FParser_size_1
 
 
    ! ==========================================================================
    ! Get input sizes before actually reading values - 2d array.
    ! ==========================================================================
    subroutine FParser_size_2(cmdname, size1, size2)
 
      character(*), intent(in)    :: cmdname
      integer,      intent(in)    :: size1
      integer,      intent(inout) :: size2
 
      character :: cmd_array(len(cmdname))
      integer :: c
 
      ! Important to initalize since array_size is accumulated.
      size2 = 0
 
      if(.not.FParser_initialized)return
 
      ! Pass an array of single chars rather than a character string.
      do c = 1,len(cmdname)
         cmd_array(c) = cmdname(c:c)
      enddo
 
      ! 2d arrays.
      call get_size2(cmd_array, size1, size2, len(cmdname))
 
    end subroutine FParser_size_2
 
 
    ! ==========================================================================
    ! Get all input sizes before actually reading values - 2d array.
    ! ==========================================================================
    subroutine FParser_sizeb_2(cmdname, size1, size2)
 
      character(*), intent(in)    :: cmdname
      integer,      intent(inout) :: size1, size2
 
      character :: cmd_array(len(cmdname))
      integer :: c
 
      ! Important to initalize since array_size is accumulated.
      size1 = 0
      size2 = 0
 
      if(.not.FParser_initialized)return
 
      ! Pass an array of single chars rather than a character string.
      do c = 1,len(cmdname)
         cmd_array(c) = cmdname(c:c)
      enddo
 
      ! 2d arrays.
      call get_sizeb2(cmd_array, size1, size2, len(cmdname))
 
    end subroutine FParser_sizeb_2
 
 
    ! ==========================================================================
    ! Get input sizes before actually reading values - 3d array.
    ! ==========================================================================
    subroutine FParser_size_3(cmdname, size1, size2, size3)
 
      character(*), intent(in)    :: cmdname
      integer,      intent(in)    :: size1, size2
      integer,      intent(inout) :: size3
 
      character :: cmd_array(len(cmdname))
      integer :: c
 
      ! Important to initalize since sizes are accumulated.
      size3 = 0
 
      if(.not.FParser_initialized)return
 
      ! Pass an array of single chars rather than a character string.
      do c = 1,len(cmdname)
         cmd_array(c) = cmdname(c:c)
      enddo
 
 
      ! 3d arrays.
      call get_size3(cmd_array, size1, size2, size3, len(cmdname))
 
    end subroutine FParser_size_3
 
 
    ! ==========================================================================
    ! Get input sizes before actually reading values - 4d array.
    ! ==========================================================================
    subroutine FParser_size_4(cmdname, size1, size2, size3, size4)
 
      character(*), intent(in)    :: cmdname
      integer,      intent(in)    :: size1, size2, size3
      integer,      intent(inout) :: size4
 
      character :: cmd_array(len(cmdname))
      integer :: c
 
      ! Important to initalize since sizes are accumulated.
      size4 = 0
 
      if(.not.FParser_initialized)return
 
      ! Pass an array of single chars rather than a character string.
      do c = 1,len(cmdname)
         cmd_array(c) = cmdname(c:c)
      enddo
 
!     ! 4d arrays.
!     call get_size4(cmd_array, size1, size2, size3, size4, len(cmdname))
 
    end subroutine FParser_size_4
 
 
    ! ==========================================================================
    ! Get input sizes before actually reading values - 5d array.
    ! ==========================================================================
    subroutine FParser_size_5(cmdname, size1, size2, size3, size4, size5)
 
      character(*), intent(in)    :: cmdname
      integer,      intent(in)    :: size1, size2, size3, size4
      integer,      intent(inout) :: size5
 
      character :: cmd_array(len(cmdname))
      integer :: c
 
      ! Important to initalize since sizes are accumulated.
      size5 = 0
 
      if(.not.FParser_initialized)return
 
      ! Pass an array of single chars rather than a character string.
      do c = 1,len(cmdname)
         cmd_array(c) = cmdname(c:c)
      enddo
 
 
      ! 5d arrays.
      call get_size4(cmd_array, size1, size2, size3, size4, size5, len(cmdname))
 
    end subroutine FParser_size_5
 
 
    ! ==========================================================================
    ! Get input sizes before actually reading values - 6d array.
    ! ==========================================================================
    subroutine FParser_size_6(cmdname, size1, size2, size3, size4, size5, size6)
 
      character(*), intent(in)    :: cmdname
      integer,      intent(in)    :: size1, size2, size3, size4, size5
      integer,      intent(inout) :: size6
 
      character :: cmd_array(len(cmdname))
      integer :: c
 
      ! Important to initalize since sizes are accumulated.
      size6 = 0
 
      if(.not.FParser_initialized)return
 
      ! Pass an array of single chars rather than a character string.
      do c = 1,len(cmdname)
         cmd_array(c) = cmdname(c:c)
      enddo
 
      ! 6d arrays.
      call get_size6(cmd_array, size1, size2, size3, size4, size5, size6, &
                     len(cmdname))
 
    end subroutine FParser_size_6
 
 
! ******************************************************************************
! ******************************************************************************
! Parser setup and meta-operations
! ******************************************************************************
! ******************************************************************************
 
    ! ==========================================================================
    ! Initialize the parser
    ! ==========================================================================
    subroutine FParser_initialize(first,deckname,oargs)
 
      character(len=*), intent(in) :: deckname,oargs
      logical         , intent(in) :: first
 
      call parser_create(oargs, len(oargs), trim(deckname), len_trim(deckname))

      call parser_comm_info(mype,numpe,iope)
      FParser_initialized = .true.
 
    end subroutine FParser_initialize

    ! ==========================================================================
    ! Check for other packages to make sure FParser has been initialized
    ! ==========================================================================
    logical function FParser_initialized_check()
       FParser_initialized_check = .false.
       if (FParser_initialized) FParser_initialized_check = .true.
    end function FParser_initialized_check

! ******************************************************************************
! ******************************************************************************
! Dictionary operations
! ******************************************************************************
! ******************************************************************************

   subroutine FParser_dictionary_add(name, value, pred, vdesc)
      use ISO_C_BINDING, only : C_NULL_CHAR

      character(len=*), intent(in) :: name
      real(REAL64),     intent(in) :: value
      logical,          intent(in) :: pred 
      character(len=*), intent(in) :: vdesc

      integer :: c_pred = 0

      if (pred) c_pred = 1
      call parser_dictionary_add(trim(name)//C_NULL_CHAR, value, c_pred,   &
                                 trim(vdesc)//C_NULL_CHAR)

   end subroutine FParser_dictionary_add

   subroutine FParser_dictionary_env_add(name, pred)
      use ISO_C_BINDING, only : C_NULL_CHAR

      character(len=*), intent(in) :: name
      logical,          intent(in) :: pred 

      integer :: c_pred = 0

      if (pred) c_pred = 1
      call parser_dictionary_env_add(trim(name)//C_NULL_CHAR, c_pred)

   end subroutine FParser_dictionary_env_add

! ******************************************************************************
! ******************************************************************************
! Compile buffer
! ******************************************************************************
! ******************************************************************************
   subroutine FParser_compile_buffer(value_returned)
      integer, intent(out), optional   :: value_returned

      integer  :: valreturn

      if (present(value_returned)) value_returned = 0
      valreturn      = 0
 
      call parser_compile_buffer(valreturn)
      if (present(value_returned)) value_returned = valreturn

   end subroutine FParser_compile_buffer
 
! ******************************************************************************
! ******************************************************************************
! Handle processed flags.
! ******************************************************************************
! ******************************************************************************
 
    ! ==========================================================================
    ! Check that every command/word has been processed.
    ! ==========================================================================
    subroutine FParser_check_processed(good)
      logical, intent(inout) :: good
 
      integer :: good_int = 1

      if(.not.FParser_initialized)then
         good=.false.
         return
      endif

      if (.not.  FParser_do_check) then
         good=.true.
         return
      endif
 
      call check_processed(good_int)
 
      good = .false.
      if (good_int .eq. 1) good = .true.
 
    end subroutine FParser_check_processed
 
 
    ! ==========================================================================
    ! Check that every command/word has been processed.
    ! We are done with the parser and could clean up here.
    ! ==========================================================================
    subroutine FParser_terminate(good, value_returned)

      ! The argument value_returned will be necessary to print the errors in the
      ! data base after parser has failed:

      logical, intent(inout)           :: good
      integer, intent(out), optional   :: value_returned

      integer  :: valreturn

      value_returned = 0
      valreturn      = 0

      if(.not.FParser_initialized)then
         good = .false.
         return
      endif

      call process_error_final(valreturn)
      if (present(value_returned)) value_returned = valreturn
      call FParser_check_processed(good)
      call parser_destroy

    end subroutine FParser_terminate

    
    ! ==========================================================================
    ! Disable the FParser check operation
    ! ==========================================================================
    subroutine FParser_disable_check()

       FParser_do_check = .false.
       
    end subroutine FParser_disable_check
 
 
    ! ==========================================================================
    ! Set the processed flag for all words for all commands that match cmdname.
    ! The value to set the processed flag to is bval.
    ! This sets the processed flag for commands in the final buffer and in the
    ! when...then final buffers.
    ! ==========================================================================
    subroutine FParser_cmd_set_processed(cmdname, bval)
      character(*), intent(in) :: cmdname
      logical     , intent(in) :: bval
 
      character :: cmd_array(len(cmdname))
      integer :: bval_i, c
 
      if(.not.FParser_initialized)return
 
      ! Pass an array of single chars rather than a character string.
      do c = 1,len(cmdname)
         cmd_array(c) = cmdname(c:c)
      enddo
 
      bval_i = 0
      if (bval) bval_i = 1
 
      call cmd_set_processed(cmd_array, len(cmdname), bval_i)
 
    end subroutine FParser_cmd_set_processed
 
 
 
! ******************************************************************************
! ******************************************************************************
! Check for command in final, parsed user input.
! ******************************************************************************
! ******************************************************************************
 
    ! ==========================================================================
    ! Check if the input command, cname, appears in the final, parsed user input.
    !
    ! The two outputs are the logicals in_input and in_whenthen
    !    in_input     command is in (or not) the main part of the input, i.e.
    !                 everything except the when...then statements.
    !    in_whenthen  command is in (or not) at least one when...then statement.
    ! ==========================================================================
    subroutine FParser_cmd_in_input(cmdname, in_input, in_whenthen)
      character(*), intent(in)    :: cmdname
      logical     , intent(inout) :: in_input, in_whenthen
 
      character :: cmd_array(len(cmdname))
      integer   :: in_input_i = 0, in_whenthen_i = 0, c
 
 
      in_input    = .false.
      in_whenthen = .false.
 
      if(.not.FParser_initialized)return
 
      ! Pass an array of single chars rather than a character string.
      do c = 1,len(cmdname)
         cmd_array(c) = cmdname(c:c)
      enddo
 
      call cmd_in_input(cmd_array, len(cmdname), in_input_i, in_whenthen_i)
 
      if (in_input_i    .eq. 1) in_input    = .true.
      if (in_whenthen_i .eq. 1) in_whenthen = .true.
 
    end subroutine FParser_cmd_in_input
 
 
 
! ******************************************************************************
! ******************************************************************************
! List variables, functions, etc, to the output file.
! ******************************************************************************
! ******************************************************************************
 
    ! ==========================================================================
    ! Echo user input to a file.
    ! ==========================================================================
    subroutine FParser_echo_user_input(kw)
 
      integer, intent(in) :: kw
      integer, PARAMETER :: MAX_LENGTH=65536 ! overkill I hope
                                             ! should get line length and alloc
                                             ! but no dynamic string alloc
                                             ! until Fortran2003
      character(MAX_LENGTH) :: fline
      integer :: nchar, bint
 
      if (mype .ne. iope) return
 
      if(.not.FParser_initialized)then
 
         write(kw, '(/,a)') "THE INPUT DECK HAS NOT BEEN OPENED, SO"
         write(kw, '(/,a)') "IT CANNOT BE PRINTED OUT"
         return
      endif
 
      write(kw, '(/,a)') "BEGIN INPUT FILE -------------------------"
 
      call echo_ui_start
      nchar = MAX_LENGTH
      do
         call get_output_line(fline, nchar, bint)
         if (bint .eq. 0) exit;
         write(kw,'(a)') trim(fline)
      enddo
 
      write(kw, '(a,/)') "END INPUT FILE -------------------------"
 
    end subroutine FParser_echo_user_input
 
    ! ==========================================================================
    ! Find the number of input variable names
    ! ==========================================================================
    subroutine FParser_count_var(kw,vnszmx)
 
      integer, intent(in)  :: kw
      integer, intent(out) :: vnszmx

      character(len=256)  :: fline
      character(len=128)  :: lineshortened
      integer             :: nchar, bint
      integer             :: i
      logical             :: ido
 
      if (mype .ne. iope) return

      ido   = .false.
 
      if(.not.FParser_initialized)then
         write(kw, '(/,a)') "THE INPUT DECK HAS NOT BEEN OPENED, SO"
         write(kw, '(/,a)') "WE CANNOT COUNT THE NUMBER OF VARIABLES"
         return
      endif
 
      write(kw, '(/,a)') "*****************************************************************"
      write(kw, '(  a)') "********** This is what the code uses to count internal variables"
 
      call FParser_list_vars (kw)

      call echo_final_buffer
      nchar = 256
      i     = 0
      do
         call get_output_line(fline, nchar, bint)
         if (bint .eq. 0) exit;
         !write(kw,'(a)') trim(fline)
         
         call  shorten_string (fline, lineshortened, ido)

         if (ido) then
           i        = i+1
         endif
      enddo

      if (i.gt.0)  then
        vnszmx = i
      else 
        vnszmx = 1
        write(kw, '(a)') "******* Could not find any variables.        ******"
        write(kw, '(a)') "******* End of counting internal variables.  ******"
      endif
 
    end subroutine FParser_count_var

    ! ==========================================================================
    ! Echo the parser functions, variables, and final buffer to file unit
    ! number kw.
    ! ==========================================================================
    subroutine FParser_echo_fvf(kw, vname, vnsize)
 
      integer, intent(in)  :: kw
      integer, intent(out) :: vnsize

      character(len=*), intent(out) :: vname(*)
 
      if (mype .ne. iope) return
 
      if(.not.FParser_initialized)then
         write(kw, '(/,a)') "THE INPUT DECK HAS NOT BEEN OPENED, SO ECHOING"
         write(kw, '(/,a)') "IT TO THE SCREEN CANNOT BE PERFORMED"
         return
      endif
 
      call FParser_list_funcs  (kw)
      call FParser_list_vars   (kw)
      call FParser_echo_fbuffer(kw, vname, vnsize)
 
    end subroutine FParser_echo_fvf

    ! ==========================================================================
    ! List available functions to a file.
    ! ==========================================================================
    subroutine FParser_list_funcs(kw)
 
      integer, intent(in) :: kw
      character(256) :: fline
      integer :: nchar, bint
 
      if (mype .ne. iope) return
 
      if(.not.FParser_initialized)then
         write(kw, '(/,a)') "THE INPUT DECK HAS NOT BEEN OPENED, SO"
         write(kw, '(/,a)') "A LIST OF ITS FUNCTIONS CANNOT BE PRINTED OUT"
         return
      endif
 
      write(kw, '(/,/,a)') "********************************************************************************"
      write(kw, '(    a)') "********** This is a list of functions known in the user input file."
 
      call list_functions_start
      nchar = 256
      do
         call get_output_line(fline, nchar, bint)
         if (bint .eq. 0) exit;
         write(kw,'(a)') trim(fline)
      enddo
 
      write(kw, '(a)') "********** End of list of functions."
      write(kw, '(a,/,/)') "********************************************************************************"
 
    end subroutine FParser_list_funcs
 
 
    ! ==========================================================================
    ! List available variables to a file.
    ! ==========================================================================
    subroutine FParser_list_vars(kw)
 
      integer, intent(in) :: kw
      character(256) :: fline
      integer :: nchar, bint
 
      if (mype .ne. iope) return
 
      if(.not.FParser_initialized)then
         write(kw, '(/,a)') "THE INPUT DECK HAS NOT BEEN OPENED, SO"
         write(kw, '(/,a)') "A LIST OF ITS VARIABLES CANNOT BE PRINTED OUT"
         return
      endif
 
      write(kw, '(/,a)') "********************************************************************************"
      write(kw, '(  a)') "********** List of input file variables."
      write(kw, '(  a)') "********** The first list is only the parser pre-defined variables."
      write(kw, '(  a)') "********** The second list is all variables at the end of processing the input"
      write(kw, '(  a)') "********** file including variables defined on the execute line."
      write(kw, '(  a)') " "
 
      call list_variables_start
      nchar = 256
      do
         call get_output_line(fline, nchar, bint)
         if (bint .eq. 0) exit;
         write(kw,'(a)') trim(fline)
      enddo
 
      write(kw, '(a)') "********** End of variable list."
      write(kw, '(a,/,/)') "********************************************************************************"
 
    end subroutine FParser_list_vars
 
 
    ! ==========================================================================
    ! Echo the final buffer to a file.
    ! ==========================================================================
    subroutine FParser_echo_fbuffer(kw, vname, vnsize)
 
      integer, intent(in)  :: kw
      integer, intent(out) :: vnsize

      character(len=*), intent(out) :: vname(*)

      character(len=256)  :: fline
      character(len=128)  :: lineshortened
      integer             :: nchar, bint
      integer             :: i
      logical             :: ido
 
      if (mype .ne. iope) return

      ido   = .false.
 
      if(.not.FParser_initialized)then
         write(kw, '(/,a)') "THE INPUT DECK HAS NOT BEEN OPENED, SO"
         write(kw, '(/,a)') "A PARSED LIST OF ITS CONTENTS CANNOT BE PRINTED OUT"
         return
      endif
 
      write(kw, '(/,a)') "********************************************************************************"
      write(kw, '(  a)') "********** Echo final parser buffer, this is what the code uses to set internal "
      write(kw, '(  a)') "********** code variables."
 
      call echo_final_buffer
      nchar = 256
      i     = 0
      do
         call get_output_line(fline, nchar, bint)
         if (bint .eq. 0) exit;
         write(kw,'(a)') trim(fline)
         
         call  shorten_string (fline, lineshortened, ido)

         if (ido) then
           i        = i+1
           vname(i) = lineshortened
         endif
      enddo
      if (i.gt.0)  vnsize = i

      ! eliminate double pname
      if (vname(1) .eq. vname(2)) then
        do i = 1, vnsize-1
           vname(i) = vname(i+1)
        enddo
        vnsize = vnsize - 1
      endif
 
      write(kw, '(a)') "********** End of echo final parser buffer."
      write(kw, '(a,/,/)') "********************************************************************************"
 
 
      write(kw, '(/,a)') "********************************************************************************"
      write(kw, '(a)')   "********** Echo final when...then parser buffers, this is what the code uses "
      write(kw, '(a)')   "********** to set internal code variables when processing when...then commands."
 
      call echo_wt_final_buffer
      nchar = 256
      do
         call get_output_line(fline, nchar, bint)
         if (bint .eq. 0) exit;
         write(kw,'(a)') trim(fline)
      enddo
 
      write(kw, '(a)')     "********** End of echo final when...then parser buffers."
      write(kw, '(a,/,/)') "********************************************************************************"
 
 
      write(kw, '(/,a)') "********************************************************************************"
      write(kw, '(a)')   "********** Echo restart block information."
 
      call echo_rb_info
      nchar = 256
      do
         call get_output_line(fline, nchar, bint)
         if (bint .eq. 0) exit;
         write(kw,'(a)') trim(fline)
      enddo
 
      write(kw, '(a)')     "********** End of echo restart block information."
      write(kw, '(a,/,/)') "********************************************************************************"
 
    end subroutine FParser_echo_fbuffer
 
! ******************************************************************************
! ******************************************************************************
! Shorten input lines for listing just input variable names
! ******************************************************************************
! ******************************************************************************

      subroutine shorten_string (fline, lineshortened, ido)

      character(*),       intent(in)    :: fline
      character(len=128), intent(out)   :: lineshortened
      logical,            intent(out)   :: ido

      character(len=128) :: lineintermed
      character(len=128) :: lineempty = " "
      integer            :: inde1, inde2
      character(len= 1)  :: deleq = "="
      character(len= 1)  :: delco = "!"
      character(len= 1)  :: deldi = "("
      character(len= 1)  :: delde = "$"

      ido   = .false.
      inde1 = 0
      inde2 = 0
      lineintermed  = trim(fline)
      lineshortened = " "

      ! check first for ( to see if it is not a dimensioned variable
      ! since I do not want to print the (dim ... on the list
      inde1 = scan(lineintermed,deldi)
      if (inde1.le.1) then
         ! Then look for the = sign
         inde2          = scan(lineintermed,deleq)
         if(inde2 .gt. 2) then
          lineshortened = lineintermed(1:inde2-1)
          ido           = .true.
         endif
      else if (inde1.gt.1) then 
        ! In this case, the equal sign comes after the (:
        ido = .true.
        lineshortened = lineintermed(1:inde1-1)
        lineshortened = trim(lineshortened)
      endif

      ! Now, check to see if the line is a comment line (comment before the =):
      inde1 = 0
      inde1 = scan(lineshortened,delco)
      if (inde1.gt.0) then
        !it is a comment
        write(*,*) "found a comment before the = or the ("
        ido = .false.
        lineshortened = lineempty
        return
      endif

      ! Now, check to see if the line is a definition line ("$" before the =):
      inde1 = 0
      inde1 = scan(lineshortened,delde)
      if (inde1.gt.0) then
        !it is a comment
        write(*,*) "found a comment before the = or the ("
        ido = .false.
        lineshortened = lineempty
        return
      endif
       
      end subroutine shorten_string
 
! ******************************************************************************
! ******************************************************************************
! Handle restart_block commands.
! ******************************************************************************
! ******************************************************************************
 
    ! ==========================================================================
    ! Do the restart_block checks
    ! ==========================================================================
    subroutine FParser_rb_check(code_varnames, code_values, code_vv_active, &
                               nchar, nvv, rb_check, rb_ntriggered, &
                               rb_num, rb_triggered_indices)
 
      integer,                             intent(in)    :: nchar, nvv, rb_num
      character(nchar), dimension(nvv),    intent(in)    :: code_varnames, code_values
      integer,          dimension(nvv),    intent(in)    :: code_vv_active
      logical,                             intent(  out) :: rb_check
      integer,                             intent(inout) :: rb_ntriggered
      integer,          dimension(rb_num), intent(inout) :: rb_triggered_indices
 
      character :: val_array(nchar*nvv), varname_array(nchar*nvv)
      integer :: i,c,i1d, rbci
 
      rb_check = .false.
      rb_ntriggered = 0
 
      if(.not.FParser_initialized)return
 
      ! We pack all the fortran strings into one long character array, this
      ! makes it easier to send to C++.
      i1d = 1
      do i = 1,nvv
         do c = 1,nchar
            val_array(i1d) = code_values(i)(c:c)
            varname_array(i1d) = code_varnames(i)(c:c)
            i1d = i1d+1
         enddo
      enddo
 
      call parser_rb_check(varname_array, val_array, code_vv_active, &
                           nchar, nvv, rbci, &
                           rb_ntriggered, rb_triggered_indices)
 
      if (rbci .eq. 1) rb_check = .true.
 
    end subroutine FParser_rb_check
 
 
    ! ==========================================================================
    ! Echo restart block data to a file.
    !     kw = unit number for the file to write to
    !     rb = index of the restart block to write. Note that this is an index
    !          which starts from 0
    ! ==========================================================================
    subroutine FParser_echo_rb_info(kw, rb)
 
      integer, intent(in) :: kw, rb
      character(256) :: fline
      integer :: nchar, bint
 
      if (mype .ne. iope) return
 
      if(.not.FParser_initialized)then
         write(kw, '(/,a)') "THE INPUT DECK HAS NOT BEEN OPENED, SO"
         write(kw, '(/,a)') "A PARSED LIST OF ITS CONTENTS CANNOT BE PRINTED OUT"
         return
      endif
 
      ! C++ function which echos the restart block info to a stringstream
      ! internal to the Parser package.
      call echo_rb1_info(rb)
 
      ! Go through every line in the internal stringstream, copy the line
      ! to a fortran string, and write the fortran string to the file.
      nchar = 256
      do
         call get_output_line(fline, nchar, bint)
         if (bint .eq. 0) exit;
         if (kw .eq. 6) then
            write(*,'(a)') trim(fline)
         else
            write(kw,'(a)') trim(fline)
         endif
      enddo
 
    end subroutine FParser_echo_rb_info
 
 
! ******************************************************************************
! ******************************************************************************
! Handle when...then commands
! ******************************************************************************
! ******************************************************************************
 
    ! ==========================================================================
    ! Get the number of when...then commands.
    ! ==========================================================================
    subroutine FParser_whenthen_num(wtnum)
      integer, intent(out) :: wtnum
      wtnum = 0
      if(.not.FParser_initialized)return
      call parser_get_wtnum(wtnum)
    end subroutine FParser_whenthen_num
 
 
    ! ==========================================================================
    ! Do the when...then checks
    ! ==========================================================================
    subroutine FParser_whenthen_check(wti, code_varnames, code_values,  &
                                     code_vv_active, nchar, nvv, wt_check)
 
      integer,                          intent(in)    :: wti, nchar, nvv
      character(nchar), dimension(nvv), intent(in)    :: code_varnames, code_values
      integer,          dimension(nvv), intent(in)    :: code_vv_active
      logical,                          intent(  out) :: wt_check
 
      character :: val_array(nchar*nvv), varname_array(nchar*nvv)
      integer :: i,c,i1d, wtci
 
      wt_check = .false.
 
      if(.not.FParser_initialized)return
 
      ! We pack all the fortran strings into one long character array, this
      ! makes it easier to send to C++.
      i1d = 1
      do i = 1,nvv
         do c = 1,nchar
            val_array(i1d) = code_values(i)(c:c)
            varname_array(i1d) = code_varnames(i)(c:c)
            i1d = i1d+1
         enddo
      enddo
 
      call parser_wt_check2(wti, varname_array, val_array, &
                            code_vv_active, nchar, nvv, wtci)
 
      if (wtci .eq. 1) wt_check = .true.
 
    end subroutine FParser_whenthen_check
 
 
    ! ==========================================================================
    ! Set the final commands pointer to the whenthen command buffer.
    ! This is also done in the whenthen check routine if the check is true.
    ! ==========================================================================
    subroutine FParser_whenthen_setcfp(wti)
      integer, intent(in)    :: wti
      if(.not.FParser_initialized)return
      call parser_wt_set_cmdsfp(wti)
    end subroutine FParser_whenthen_setcfp
 
 
    ! ==========================================================================
    ! Reset the commands final buffer pointer.
    ! ==========================================================================
    subroutine FParser_whenthen_reset
      if(.not.FParser_initialized)return
      call parser_wt_reset
    end subroutine FParser_whenthen_reset
 
 
    ! ==========================================================================
    ! Get the character array and its size describing a whenthen.
    ! ==========================================================================
    subroutine FParser_whenthen_ca(wti, wt_ca, wt_casize)
      integer, intent(in)    :: wti, wt_casize
      character(1), dimension(wt_casize), intent(inout) :: wt_ca
      wt_ca = ' '
      if(.not.FParser_initialized)return
      call parser_wt_ca(wti, wt_ca, wt_casize)
    end subroutine FParser_whenthen_ca
 
    subroutine FParser_whenthen_casize(wti, wt_casize)
      integer, intent(in)    :: wti
      integer, intent(inout) :: wt_casize
      wt_casize=0
      if(.not.FParser_initialized)return
      call parser_wt_casize(wti, wt_casize)
    end subroutine FParser_whenthen_casize
 
 
    ! ==========================================================================
    ! Get the number of satisfied flags for a whenthen. This is equal to the
    ! number of sub-conditions in the whenthen condition.
    ! ==========================================================================
    subroutine FParser_whenthen_satsize(wti, wt_satsize)
      integer, intent(in)    :: wti
      integer, intent(inout) :: wt_satsize
      wt_satsize=0
      if(.not.FParser_initialized)return
      call parser_wt_satsize(wti, wt_satsize)
    end subroutine FParser_whenthen_satsize
 
 
    ! ==========================================================================
    ! Get and Set the whenthen satisfied flags for a whenthen.
    ! ==========================================================================
    subroutine FParser_whenthen_getsat(wti, wt_sat, wt_satsize)
      integer, intent(in)    :: wti, wt_satsize
      integer, dimension(wt_satsize), intent(inout) :: wt_sat
      wt_sat=0
      if(.not.FParser_initialized)return
      call parser_wt_getsat(wti, wt_sat, wt_satsize)
    end subroutine FParser_whenthen_getsat
 
    subroutine FParser_whenthen_setsat(wti, wt_sat, wt_satsize)
      integer, intent(in)    :: wti, wt_satsize
      integer, dimension(wt_satsize), intent(inout) :: wt_sat
      wt_sat=0
      if(.not.FParser_initialized)return
      call parser_wt_setsat(wti, wt_sat, wt_satsize)
    end subroutine FParser_whenthen_setsat
 
 
    ! ==========================================================================
    ! Get and Set the processed flag for a whenthen.
    ! ==========================================================================
    subroutine FParser_whenthen_getprocessed(wti, wtp)
      integer, intent(in)    :: wti
      integer, intent(inout) :: wtp
      wtp=0
      if(.not.FParser_initialized)return
      call parser_wt_getprocessed(wti, wtp)
    end subroutine FParser_whenthen_getprocessed
 
    subroutine FParser_whenthen_setprocessed(wti, wtp)
      integer, intent(in) :: wti
      integer, intent(in) :: wtp
      if(.not.FParser_initialized)return
      call parser_wt_setprocessed(wti, wtp)
    end subroutine FParser_whenthen_setprocessed
 
 
    ! ==========================================================================
    ! Get and Set the sequence number for a whenthen.
    ! ==========================================================================
    subroutine FParser_whenthen_getseq(wti, wtseq)
      integer, intent(in)    :: wti
      integer, intent(inout) :: wtseq
      wtseq=0
      if(.not.FParser_initialized)return
      call parser_wt_getseq(wti, wtseq)
    end subroutine FParser_whenthen_getseq
 
    subroutine FParser_whenthen_setseq(wti, wtseq)
      integer, intent(in) :: wti
      integer, intent(in) :: wtseq
      if(.not.FParser_initialized)return
      call parser_wt_setseq(wti, wtseq)
    end subroutine FParser_whenthen_setseq
 
    integer function FParser_get_num_include_files()
      interface
        subroutine parser_get_num_include_files(num_include) bind(C,NAME="parser_get_num_include_files")
          use iso_c_binding, only : C_INT
          integer(C_INT), intent(out) :: num_include
        end subroutine parser_get_num_include_files
      end interface
      call parser_get_num_include_files(FParser_get_num_include_files)
    end function FParser_get_num_include_files

    subroutine FParser_list_include_files() !Mainly for debugging
      use iso_c_binding
      interface
        subroutine parser_list_include_files() bind(C,NAME="parser_list_include_files")
        end subroutine parser_list_include_files
      end interface
      call parser_list_include_files()
    end subroutine FParser_list_include_files

    subroutine FParser_get_include_file(kw,fname,i)
      use iso_c_binding
      interface
        integer(C_INT) function parser_get_include_file_name_length(n) bind(C,NAME="parser_get_include_file_name_length")
          use iso_c_binding, only : C_INT
          integer(C_INT), intent(in) :: n
        end function parser_get_include_file_name_length
      end interface
      interface
        subroutine parser_get_include_file_name(cfile_name,cfile_name_len,n) bind(C,NAME="parser_get_include_file_name")
          use iso_c_binding, only : C_INT, C_CHAR
          integer(C_INT), intent(in) :: cfile_name_len, n
          character(kind=C_CHAR,len=1), dimension(cfile_name_len), intent(in) :: cfile_name
        end subroutine parser_get_include_file_name
      end interface
      character(len=120), intent(out) :: fname
      integer, intent(in) :: i, kw

      integer(C_INT) :: n, file_name_length, cfile_name_length, ii
      character(kind=c_char,len=1), dimension(:), allocatable  :: cfile_name

      n = i-1
      file_name_length = parser_get_include_file_name_length(n)
      cfile_name_length = file_name_length + 1
      allocate(cfile_name(cfile_name_length),stat=FParser_allostat)
      call FParser_test_allostat(kw, "cfile_name allocate")
      call parser_get_include_file_name(cfile_name,cfile_name_length,n)
!     allocate(character(len=file_name_length)::fname)
      do ii = 1,file_name_length
        fname(ii:ii) = cfile_name(ii)
      enddo
      deallocate(cfile_name,stat=FParser_allostat)
      call FParser_test_allostat(kw, "cfile_name deallocate")
    end subroutine FParser_get_include_file

    subroutine FParser_test_allostat(kw, comment)
    !*******************************************************************************
    !                                                                              *
    ! Test allostat and abort in not zero                                          *
    !                                                                              *
    !*******************************************************************************
    character(*), intent(in) :: comment
    integer, intent(in) :: kw

      if(FParser_allostat.ne.0)then
         write(kw,"(a,i6,2a)")'TEST_ALLOSTAT: FParser_allostat = ',FParser_allostat,' - ',trim(comment)

!        call MPI_Abort(MPI_COMM_WORLD,1,mpierror)

         stop
      endif

    end subroutine FParser_test_allostat
 
! ------------------------------------------------------------------------------
  end module FParser_module
! ==============================================================================
