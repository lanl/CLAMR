!*
!*  Copyright (c) 2010-2014, Los Alamos National Security, LLC.
!*  All rights Reserved.
!*
!*  Copyright 2010-2014. Los Alamos National Security, LLC. This software was produced 
!*  under U.S. Government contract DE-AC52-06NA25396 for Los Alamos National 
!*  Laboratory (LANL), which is operated by Los Alamos National Security, LLC 
!*  for the U.S. Department of Energy. The U.S. Government has rights to use, 
!*  reproduce, and distribute this software.  NEITHER THE GOVERNMENT NOR LOS 
!*  ALAMOS NATIONAL SECURITY, LLC MAKES ANY WARRANTY, EXPRESS OR IMPLIED, OR 
!*  ASSUMES ANY LIABILITY FOR THE USE OF THIS SOFTWARE.  If software is modified
!*  to produce derivative works, such modified software should be clearly marked,
!*  so as not to confuse it with the version available from LANL.
!*
!*  Additionally, redistribution and use in source and binary forms, with or without
!*  modification, are permitted provided that the following conditions are met:
!*     * Redistributions of source code must retain the above copyright
!*       notice, this list of conditions and the following disclaimer.
!*     * Redistributions in binary form must reproduce the above copyright
!*       notice, this list of conditions and the following disclaimer in the
!*       documentation and/or other materials provided with the distribution.
!*     * Neither the name of the Los Alamos National Security, LLC, Los Alamos 
!*       National Laboratory, LANL, the U.S. Government, nor the names of its 
!*       contributors may be used to endorse or promote products derived from 
!*       this software without specific prior written permission.
!*  
!*  THIS SOFTWARE IS PROVIDED BY THE LOS ALAMOS NATIONAL SECURITY, LLC AND 
!*  CONTRIBUTORS "AS IS" AND ANY EXPRESS OR IMPLIED WARRANTIES, INCLUDING, BUT 
!*  NOT LIMITED TO, THE IMPLIED WARRANTIES OF MERCHANTABILITY AND FITNESS FOR
!*  A PARTICULAR PURPOSE ARE DISCLAIMED. IN NO EVENT SHALL LOS ALAMOS NATIONAL
!*  SECURITY, LLC OR CONTRIBUTORS BE LIABLE FOR ANY DIRECT, INDIRECT, INCIDENTAL,
!*  SPECIAL, EXEMPLARY, OR CONSEQUENTIAL DAMAGES (INCLUDING, BUT NOT LIMITED TO,
!*  PROCUREMENT OF SUBSTITUTE GOODS OR SERVICES; LOSS OF USE, DATA, OR PROFITS;
!*  OR BUSINESS INTERRUPTION) HOWEVER CAUSED AND ON ANY THEORY OF LIABILITY,
!*  WHETHER IN CONTRACT, STRICT LIABILITY, OR TORT (INCLUDING NEGLIGENCE OR OTHERWISE)
!*  ARISING IN ANY WAY OUT OF THE USE OF THIS SOFTWARE, EVEN IF ADVISED OF THE
!*  POSSIBILITY OF SUCH DAMAGE.
!*  
!*  CLAMR -- LA-CC-11-094
!*  This research code is being developed as part of the 
!*  2011 X Division Summer Workshop for the express purpose
!*  of a collaborative code for development of ideas in
!*  the implementation of AMR codes for Exascale platforms
!*  
!*  Authors: Chuck Wingate   XCP-2   caw@lanl.gov
!*           Robert Robey    XCP-2   brobey@lanl.gov
!*

! ------------------------------------------------------------------------------
! This module contains a collection of subroutines and related data that
! is used to read input data from an input file.
! ------------------------------------------------------------------------------
  module PPinput_module
 
    implicit none
    save
    private

    integer, parameter :: INT4   = SELECTED_INT_KIND(6)
    integer, parameter :: INT32  = SELECTED_INT_KIND(6)
    integer, parameter :: INT8   = SELECTED_INT_KIND(16)
    integer, parameter :: INT64  = SELECTED_INT_KIND(16)

    integer, parameter :: REAL4  = SELECTED_REAL_KIND(6)
    integer, parameter :: REAL8  = SELECTED_REAL_KIND(12)

    integer, parameter :: REAL32 = SELECTED_REAL_KIND(6)
    integer, parameter :: REAL64 = SELECTED_REAL_KIND(12)

    integer :: mype
    integer :: numpe
    integer :: iope

    public PPinput                    ! Interface to get command values.
    public PPinput_size               ! Interface to get sizes before getting values.
    public PPinput_sizeb              ! Interface to get all sizes before getting values.
 
    public PPinput_initialize         ! Initialize the parser and read in the input deck.
    public PPinput_dictionary_add     ! Add a constant to the Parser dictionary
    public PPinput_dictionary_env_add ! Add an environment variable constant to the Parser dictionary
    public PPinput_compile_buffer     ! Compile the input buffer
    public PPinput_initialized_check  ! Interface to check if PPinput has been initialized
    public PPinput_check_processed    ! Check processed flags.
    public PPinput_disable_check      ! Disable check processed.
    public PPinput_terminate          ! Last PPinput call, check processed flags, clean.
    public PPinput_cmd_in_input       ! Check for command in the user input.
    public PPinput_cmd_set_processed  ! Set the command as being processed.
    public PPinput_echo_user_input    ! Echo user input to a file.
    public PPinput_echo_fvf           ! Echo functions, variables, final buffer
 
    private PPinput_list_funcs        ! List functions to file.
    private PPinput_list_vars         ! List variables to file.
    private PPinput_echo_fbuffer      ! Echo final buffer to a file.
 
    ! Define an interface block that allows several subroutines to be called
    ! using a common name, i.e. PPinput. This is for getting values.
    interface PPinput
       module procedure PPinput_c_00   ! Single character
       module procedure PPinput_c_0    ! Character string
       module procedure PPinput_c_1    ! 1d array of character strings
       module procedure PPinput_c_2    ! 2d array of character strings
       module procedure PPinput_c_3    ! 3d array of character strings
       module procedure PPinput_c_4    ! 4d array of character strings
 
       module procedure PPinput_i_0    ! Single integer
       module procedure PPinput_i_1    ! 1d integer array
       module procedure PPinput_i_2    ! 2d integer array
       module procedure PPinput_i_3    ! 3d integer array
       module procedure PPinput_i_4    ! 4d integer array
 
       module procedure PPinput_i8_0   ! Single long int (int64_t)
       module procedure PPinput_i8_1   ! Single long int (int64_t)
 
       module procedure PPinput_l_0    ! Single logical
       module procedure PPinput_l_1    ! 1d logical array
       module procedure PPinput_l_2    ! 2d logical array
       module procedure PPinput_l_3    ! 3d logical array
       module procedure PPinput_l_4    ! 4d logical array
 
       module procedure PPinput_r8_0   ! Single real
       module procedure PPinput_r8_1   ! 1d real array
       module procedure PPinput_r8_2   ! 2d real array
       module procedure PPinput_r8_3   ! 3d real array
       module procedure PPinput_r8_4   ! 4d real array
    end interface
 
    ! Define an interface block that allows several subroutines to be called
    ! using a common name, i.e. PPinput_size. This is for getting sizes.
    ! For variables that do not exist, a size of 0 must be returned.
    ! This size is used as a check and passed into other routines.
    ! So returning a negative number for the size of non-existant vars
    ! does cause problems.
    interface PPinput_size
       module procedure PPinput_size_1    ! Single dimensional array
       module procedure PPinput_size_2    ! 2d array
       module procedure PPinput_size_3    ! 3d array
       module procedure PPinput_size_4    ! 4d array
       module procedure PPinput_size_5    ! 5d array
       module procedure PPinput_size_6    ! 6d array
    end interface
    ! Define an interface block that allows several subroutines to be called
    ! using a common name, i.e. PPinput_sizeb. This is for getting all sizes.
    interface PPinput_sizeb
       module procedure PPinput_sizeb_2   ! 2d array
    end interface
 
    ! restart_block subroutines.
    public PPinput_rb_check           ! Do the restart blocks check
    public PPinput_echo_rb_info       ! Echo restart block info for one block
 
 
    ! When...then subroutines.
    public PPinput_whenthen_num       ! Get number of when...thens
    public PPinput_whenthen_check     ! Do the whenthen check
    public PPinput_whenthen_setcfp
    public PPinput_whenthen_reset     ! Reset cmds final buffer pointer.
    public PPinput_whenthen_casize
    public PPinput_whenthen_ca
    public PPinput_whenthen_satsize
    public PPinput_whenthen_getsat
    public PPinput_whenthen_setsat
    public PPinput_whenthen_getprocessed
    public PPinput_whenthen_setprocessed
    public PPinput_whenthen_getseq
    public PPinput_whenthen_setseq
 
    ! Private module data.
 
 
    logical :: PPinput_initialized = .false.
    logical :: PPinput_do_check    = .true.
 
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
    subroutine PPinput_l_0(cmdname, cmdvalue, noskip)
 
      character(*), intent(in)    :: cmdname
      logical,      intent(inout) :: cmdvalue
      logical,      intent(in)    :: noskip
 
      character :: cmd_array(len(cmdname))
      integer :: ival, iskip, c
 
      if(.not.PPinput_initialized)return
 
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
 
    end subroutine PPinput_l_0
 
 
    ! ==========================================================================
    ! Get a one dimensional logical array.
    ! ==========================================================================
    subroutine PPinput_l_1(cmdname, array, noskip, array_size)
 
      character(*), intent(in)    :: cmdname
      integer,      intent(in)    :: array_size
      logical,      intent(inout) :: array(array_size)
      logical,      intent(in)    :: noskip
 
      character :: cmd_array(len(cmdname))
      integer :: ival_array(array_size), i, iskip, c
 
      if(.not.PPinput_initialized)return
 
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
 
    end subroutine PPinput_l_1
 
 
    ! ==========================================================================
    ! Get two dimensional logical array values.
    ! ==========================================================================
    subroutine PPinput_l_2(cmdname, array, noskip, size1, size2)
 
      character(*), intent(in)    :: cmdname
      integer,      intent(in)    :: size1, size2
      logical,      intent(inout) :: array(size1, size2)
      logical,      intent(in)    :: noskip
 
      character :: cmd_array(len(cmdname))
      integer :: val_array(size1*size2)
      integer :: i, j, i1d, iskip, c
 
      if(.not.PPinput_initialized)return
 
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
 
    end subroutine PPinput_l_2
 
 
    ! ==========================================================================
    ! Get three dimensional logical array values.
    ! ==========================================================================
    subroutine PPinput_l_3(cmdname, array, noskip, size1, size2, size3)
 
      character(*), intent(in)    :: cmdname
      integer,      intent(in)    :: size1, size2, size3
      logical,      intent(inout) :: array(size1, size2, size3)
      logical,      intent(in)    :: noskip
 
      character :: cmd_array(len(cmdname))
      integer :: val_array(size1*size2*size3)
      integer :: i, j, k, i1d, iskip, c
 
      if(.not.PPinput_initialized)return
 
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
 
    end subroutine PPinput_l_3
 
 
    ! ==========================================================================
    ! Get four dimensional logical array values.
    ! ==========================================================================
    subroutine PPinput_l_4(cmdname, array, noskip, size1, size2, size3, size4)
 
      character(*), intent(in)    :: cmdname
      integer,      intent(in)    :: size1, size2, size3, size4
      logical,      intent(inout) :: array(size1, size2, size3, size4)
      logical,      intent(in)    :: noskip
 
      character :: cmd_array(len(cmdname))
      integer :: val_array(size1*size2*size3*size4)
      integer :: i, j, k, l, i1d, iskip, c
 
      if(.not.PPinput_initialized)return
 
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
 
    end subroutine PPinput_l_4
 
 
 
 
! ******************************************************************************
! ******************************************************************************
! Get integer values
! ******************************************************************************
! ******************************************************************************
 
    ! ==========================================================================
    ! Get a single integer value.
    ! ==========================================================================
    subroutine PPinput_i_0(cmdname, cmdvalue, noskip)
 
      character(*), intent(in)    :: cmdname
      integer,      intent(inout) :: cmdvalue
      logical,      intent(in)    :: noskip
 
      character :: cmd_array(len(cmdname))
      integer :: iskip, c
 
      if(.not.PPinput_initialized)return
 
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
 
    end subroutine PPinput_i_0
 
 
    ! ==========================================================================
    ! Get a one dimensional integer array.
    ! ==========================================================================
    subroutine PPinput_i_1(cmdname, array, noskip, array_size)
 
      character(*), intent(in)    :: cmdname
      integer,      intent(in)    :: array_size
      integer,      intent(inout) :: array(array_size)
      logical,      intent(in)    :: noskip
 
      character :: cmd_array(len(cmdname))
      integer :: iskip, c
 
      if(.not.PPinput_initialized)return
 
      ! Pass an array of single chars rather than a character string.
      do c = 1,len(cmdname)
         cmd_array(c) = cmdname(c:c)
      enddo
 
      ! Cannot pass logicals to C++ as booleans, so use integer.
      iskip = 1
      if (noskip) iskip = 0
 
      ! This is a C++ funtion.
      call get_integer1(cmd_array, array, array_size, len(cmdname), iskip)
 
    end subroutine PPinput_i_1
 
 
    ! ==========================================================================
    ! Get two dimensional integer array values.
    ! ==========================================================================
    subroutine PPinput_i_2(cmdname, array, noskip, size1, size2)
 
      character(*), intent(in)    :: cmdname
      integer,      intent(in)    :: size1, size2
      integer,      intent(inout) :: array(size1, size2)
      logical,      intent(in)    :: noskip
 
      character :: cmd_array(len(cmdname))
      integer :: iskip, c
 
      if(.not.PPinput_initialized)return
 
      ! Pass an array of single chars rather than a character string.
      do c = 1,len(cmdname)
         cmd_array(c) = cmdname(c:c)
      enddo
 
      ! Cannot pass logicals to C++ as booleans, so use integer.
      iskip = 1
      if (noskip) iskip = 0
 
      ! This is a C++ funtion.
      call get_integer2(cmd_array, array, size1, size2, len(cmdname), iskip)
 
    end subroutine PPinput_i_2
 
 
    ! ==========================================================================
    ! Get three dimensional integer array values.
    ! ==========================================================================
    subroutine PPinput_i_3(cmdname, array, noskip, size1, size2, size3)
 
      character(*), intent(in)    :: cmdname
      integer,      intent(in)    :: size1, size2, size3
      integer,      intent(inout) :: array(size1, size2, size3)
      logical,      intent(in)    :: noskip
 
      character :: cmd_array(len(cmdname))
      integer :: iskip, c
 
      if(.not.PPinput_initialized)return
 
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
 
    end subroutine PPinput_i_3
 
 
    ! ==========================================================================
    ! Get four dimensional integer array values.
    ! ==========================================================================
    subroutine PPinput_i_4(cmdname, array, noskip, size1, size2, size3, size4)
 
      character(*), intent(in)    :: cmdname
      integer,      intent(in)    :: size1, size2, size3, size4
      integer,      intent(inout) :: array(size1, size2, size3, size4)
      logical,      intent(in)    :: noskip
 
      character :: cmd_array(len(cmdname))
      integer :: iskip, c
 
      if(.not.PPinput_initialized)return
 
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
 
    end subroutine PPinput_i_4
 
    ! ==========================================================================
    ! Get a single long integer (int64_t) value.
    ! ==========================================================================
    subroutine PPinput_i8_0(cmdname, cmdvalue, noskip)
 
      character(*),  intent(in)    :: cmdname
      integer(INT8), intent(inout) :: cmdvalue
      logical,       intent(in)    :: noskip
 
      character :: cmd_array(len(cmdname))
      integer :: iskip, c
 
      if(.not.PPinput_initialized)return
 
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
 
    end subroutine PPinput_i8_0

    ! ==========================================================================
    ! Get a one dimensional integer (int64_t) array.
    ! ==========================================================================
    subroutine PPinput_i8_1(cmdname, array, noskip, array_size)
 
      character(*),  intent(in)    :: cmdname
      integer,       intent(in)    :: array_size
      integer(INT8), intent(inout) :: array(array_size)
      logical,       intent(in)    :: noskip
 
      character :: cmd_array(len(cmdname))
      integer :: iskip, c
 
      if(.not.PPinput_initialized)return
 
      ! Pass an array of single chars rather than a character string.
      do c = 1,len(cmdname)
         cmd_array(c) = cmdname(c:c)
      enddo
 
      ! Cannot pass logicals to C++ as booleans, so use integer.
      iskip = 1
      if (noskip) iskip = 0
 
      ! This is a C++ funtion.
      call get_int8_1(cmd_array, array, array_size, len(cmdname), iskip)
 
    end subroutine PPinput_i8_1
 
 
 
 
! ******************************************************************************
! ******************************************************************************
! Get real values
! ******************************************************************************
! ******************************************************************************
 
    ! ==========================================================================
    ! Get a single real value.
    ! ==========================================================================
    subroutine PPinput_r8_0(cmdname, cmdvalue, noskip)
 
      character(*), intent(in)    :: cmdname
      real(REAL8),  intent(inout) :: cmdvalue
      logical,      intent(in)    :: noskip
 
      character :: cmd_array(len(cmdname))
      integer :: iskip, c
 
      if(.not.PPinput_initialized)return
 
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
 
    end subroutine PPinput_r8_0
 
 
    ! ==========================================================================
    ! Get one dimensional real array values.
    ! ==========================================================================
    subroutine PPinput_r8_1(cmdname, array, noskip, array_size)
 
      character(*), intent(in)    :: cmdname
      integer,      intent(in)    :: array_size
      real(REAL8),  intent(inout) :: array(array_size)
      logical,      intent(in)    :: noskip
 
      character :: cmd_array(len(cmdname))
      integer :: iskip, c
 
      if(.not.PPinput_initialized)return
 
      ! Pass an array of single chars rather than a character string.
      do c = 1,len(cmdname)
         cmd_array(c) = cmdname(c:c)
      enddo
 
      ! Cannot pass logicals to C++ as booleans, so use integer.
      iskip = 1
      if (noskip) iskip = 0
 
      ! This is a C++ funtion.
      call get_real1(cmd_array, array, array_size, len(cmdname), iskip)
 
    end subroutine PPinput_r8_1
 
 
    ! ==========================================================================
    ! Get two dimensional real array values.
    ! ==========================================================================
    subroutine PPinput_r8_2(cmdname, array, noskip, size1, size2)
 
      character(*), intent(in)    :: cmdname
      integer,      intent(in)    :: size1, size2
      real(REAL8),  intent(inout) :: array(size1, size2)
      logical,      intent(in)    :: noskip
 
      character :: cmd_array(len(cmdname))
      integer :: iskip, c
 
      if(.not.PPinput_initialized)return
 
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
 
    end subroutine PPinput_r8_2
 
 
    ! ==========================================================================
    ! Get three dimensional real array values.
    ! ==========================================================================
    subroutine PPinput_r8_3(cmdname, array, noskip, size1, size2, size3)
 
      character(*), intent(in)    :: cmdname
      integer,      intent(in)    :: size1, size2, size3
      real(REAL8),  intent(inout) :: array(size1, size2, size3)
      logical,      intent(in)    :: noskip
 
      character :: cmd_array(len(cmdname))
      integer :: iskip, c
 
      if(.not.PPinput_initialized)return
 
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
 
    end subroutine PPinput_r8_3
 
 
    ! ==========================================================================
    ! Get four dimensional real array values.
    ! ==========================================================================
    subroutine PPinput_r8_4(cmdname, array, noskip, size1, size2, size3, size4)
 
      character(*), intent(in)    :: cmdname
      integer,      intent(in)    :: size1, size2, size3, size4
      real(REAL8),  intent(inout) :: array(size1, size2, size3, size4)
      logical,      intent(in)    :: noskip
 
      character :: cmd_array(len(cmdname))
      integer :: iskip, c
 
      if(.not.PPinput_initialized)return
 
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
 
    end subroutine PPinput_r8_4
 
 
! ******************************************************************************
! ******************************************************************************
! Get character values
! ******************************************************************************
! ******************************************************************************
 
    ! ==========================================================================
    ! Get a single character value.
    ! ==========================================================================
    subroutine PPinput_c_00(cmdname, cmdvalue, noskip)
 
      character(*), intent(in)    :: cmdname
      character(*), intent(inout) :: cmdvalue
      logical,      intent(in)    :: noskip
 
      character :: cmd_array(len(cmdname))
      integer :: iskip, c
 
      if(.not.PPinput_initialized)return
 
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
 
    end subroutine PPinput_c_00
 
 
    ! ==========================================================================
    ! Get a character string.
    ! ==========================================================================
    subroutine PPinput_c_0(cmdname, cmdvalue, noskip, nchar, nchar_actual)
 
      character(*),     intent(in)            :: cmdname
      integer,          intent(in)            :: nchar
      character(nchar), intent(inout)         :: cmdvalue
      logical,          intent(in)            :: noskip
      integer,          intent(out), optional :: nchar_actual ! the actual max len
 
      character :: cmd_array(len(cmdname)), val_array(nchar)
      integer :: iskip, c
      integer :: nchar_actual_val
 
      if(.not.PPinput_initialized)return
 
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
 
    end subroutine PPinput_c_0
 
 
    ! ==========================================================================
    ! Get one dimensional array of character strings.
    ! ==========================================================================
    subroutine PPinput_c_1(cmdname, array, noskip, nchar, size1)
 
      character(*),     intent(in)    :: cmdname
      integer,          intent(in)    :: size1
      integer,          intent(in)    :: nchar
      character(nchar), intent(inout) :: array(size1)
      logical,          intent(in)    :: noskip
 
      character :: cmd_array(len(cmdname)), val_array(nchar*size1)
      integer :: i, c, i1d, iskip
 
      if(.not.PPinput_initialized)return
 
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
 
    end subroutine PPinput_c_1
 
 
    ! ==========================================================================
    ! Get two dimensional array of character strings.
    ! ==========================================================================
    subroutine PPinput_c_2(cmdname, array, noskip, nchar, size1, size2)
 
      character(*),     intent(in)    :: cmdname
      integer,          intent(in)    :: size1, size2
      integer,          intent(in)    :: nchar
      character(nchar), intent(inout) :: array(size1, size2)
      logical,          intent(in)    :: noskip
 
      character :: cmd_array(len(cmdname)), val_array(nchar*size1*size2)
      integer :: i, j, c, i1d, iskip
 
      if(.not.PPinput_initialized)return
 
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
 
    end subroutine PPinput_c_2
 
 
    ! ==========================================================================
    ! Get three dimensional array of character strings.
    ! ==========================================================================
    subroutine PPinput_c_3(cmdname, array, noskip, nchar, size1, size2, size3)
 
      character(*),     intent(in)    :: cmdname
      integer,          intent(in)    :: size1, size2, size3
      integer,          intent(in)    :: nchar
      character(nchar), intent(inout) :: array(size1, size2, size3)
      logical,          intent(in)    :: noskip
 
      character :: cmd_array(len(cmdname)), val_array(nchar*size1*size2*size3)
      integer :: i, j, k, c, i1d, iskip
 
      if(.not.PPinput_initialized)return
 
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
 
    end subroutine PPinput_c_3
 
 
    ! ==========================================================================
    ! Get four dimensional array of character strings.
    ! ==========================================================================
    subroutine PPinput_c_4(cmdname, array, noskip, nchar, size1, size2, size3, &
                          size4)
 
      character(*),     intent(in)    :: cmdname
      integer,          intent(in)    :: size1, size2, size3, size4
      integer,          intent(in)    :: nchar
      character(nchar), intent(inout) :: array(size1, size2, size3, size4)
      logical,          intent(in)    :: noskip
 
      character :: cmd_array(len(cmdname))
      character :: val_array(nchar*size1*size2*size3*size4)
      integer :: i,j,k,l,c,i1d, iskip
 
      if(.not.PPinput_initialized)return
 
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
 
    end subroutine PPinput_c_4
 
 
 
 
 
! ******************************************************************************
! ******************************************************************************
! Get sizes
! ******************************************************************************
! ******************************************************************************
 
    ! ==========================================================================
    ! Get input sizes before actually reading values - 1d array.
    ! ==========================================================================
    subroutine PPinput_size_1(cmdname, array_size)
 
      character(*), intent(in)    :: cmdname
      integer,      intent(inout) :: array_size
 
      character :: cmd_array(len(cmdname))
      integer :: c
 
      ! Important to initalize since array_size is accumulated.
      array_size = 0
 
      if(.not.PPinput_initialized)return
 
      ! Pass an array of single chars rather than a character string.
      do c = 1,len(cmdname)
         cmd_array(c) = cmdname(c:c)
      enddo
 
      ! One dimensional arrays.
      call get_size1(cmd_array, array_size, len(cmdname))
 
    end subroutine PPinput_size_1
 
 
    ! ==========================================================================
    ! Get input sizes before actually reading values - 2d array.
    ! ==========================================================================
    subroutine PPinput_size_2(cmdname, size1, size2)
 
      character(*), intent(in)    :: cmdname
      integer,      intent(in)    :: size1
      integer,      intent(inout) :: size2
 
      character :: cmd_array(len(cmdname))
      integer :: c
 
      ! Important to initalize since array_size is accumulated.
      size2 = 0
 
      if(.not.PPinput_initialized)return
 
      ! Pass an array of single chars rather than a character string.
      do c = 1,len(cmdname)
         cmd_array(c) = cmdname(c:c)
      enddo
 
      ! 2d arrays.
      call get_size2(cmd_array, size1, size2, len(cmdname))
 
    end subroutine PPinput_size_2
 
 
    ! ==========================================================================
    ! Get all input sizes before actually reading values - 2d array.
    ! ==========================================================================
    subroutine PPinput_sizeb_2(cmdname, size1, size2)
 
      character(*), intent(in)    :: cmdname
      integer,      intent(inout) :: size1, size2
 
      character :: cmd_array(len(cmdname))
      integer :: c
 
      ! Important to initalize since array_size is accumulated.
      size1 = 0
      size2 = 0
 
      if(.not.PPinput_initialized)return
 
      ! Pass an array of single chars rather than a character string.
      do c = 1,len(cmdname)
         cmd_array(c) = cmdname(c:c)
      enddo
 
      ! 2d arrays.
      call get_sizeb2(cmd_array, size1, size2, len(cmdname))
 
    end subroutine PPinput_sizeb_2
 
 
    ! ==========================================================================
    ! Get input sizes before actually reading values - 3d array.
    ! ==========================================================================
    subroutine PPinput_size_3(cmdname, size1, size2, size3)
 
      character(*), intent(in)    :: cmdname
      integer,      intent(in)    :: size1, size2
      integer,      intent(inout) :: size3
 
      character :: cmd_array(len(cmdname))
      integer :: c
 
      ! Important to initalize since sizes are accumulated.
      size3 = 0
 
      if(.not.PPinput_initialized)return
 
      ! Pass an array of single chars rather than a character string.
      do c = 1,len(cmdname)
         cmd_array(c) = cmdname(c:c)
      enddo
 
 
      ! 3d arrays.
      call get_size3(cmd_array, size1, size2, size3, len(cmdname))
 
    end subroutine PPinput_size_3
 
 
    ! ==========================================================================
    ! Get input sizes before actually reading values - 4d array.
    ! ==========================================================================
    subroutine PPinput_size_4(cmdname, size1, size2, size3, size4)
 
      character(*), intent(in)    :: cmdname
      integer,      intent(in)    :: size1, size2, size3
      integer,      intent(inout) :: size4
 
      character :: cmd_array(len(cmdname))
      integer :: c
 
      ! Important to initalize since sizes are accumulated.
      size4 = 0
 
      if(.not.PPinput_initialized)return
 
      ! Pass an array of single chars rather than a character string.
      do c = 1,len(cmdname)
         cmd_array(c) = cmdname(c:c)
      enddo
 
      ! 4d arrays.
      call get_size4(cmd_array, size1, size2, size3, size4, len(cmdname))
 
    end subroutine PPinput_size_4
 
 
    ! ==========================================================================
    ! Get input sizes before actually reading values - 5d array.
    ! ==========================================================================
    subroutine PPinput_size_5(cmdname, size1, size2, size3, size4, size5)
 
      character(*), intent(in)    :: cmdname
      integer,      intent(in)    :: size1, size2, size3, size4
      integer,      intent(inout) :: size5
 
      character :: cmd_array(len(cmdname))
      integer :: c
 
      ! Important to initalize since sizes are accumulated.
      size5 = 0
 
      if(.not.PPinput_initialized)return
 
      ! Pass an array of single chars rather than a character string.
      do c = 1,len(cmdname)
         cmd_array(c) = cmdname(c:c)
      enddo
 
 
      ! 5d arrays.
      call get_size4(cmd_array, size1, size2, size3, size4, size5, len(cmdname))
 
    end subroutine PPinput_size_5
 
 
    ! ==========================================================================
    ! Get input sizes before actually reading values - 6d array.
    ! ==========================================================================
    subroutine PPinput_size_6(cmdname, size1, size2, size3, size4, size5, size6)
 
      character(*), intent(in)    :: cmdname
      integer,      intent(in)    :: size1, size2, size3, size4, size5
      integer,      intent(inout) :: size6
 
      character :: cmd_array(len(cmdname))
      integer :: c
 
      ! Important to initalize since sizes are accumulated.
      size6 = 0
 
      if(.not.PPinput_initialized)return
 
      ! Pass an array of single chars rather than a character string.
      do c = 1,len(cmdname)
         cmd_array(c) = cmdname(c:c)
      enddo
 
      ! 6d arrays.
      call get_size6(cmd_array, size1, size2, size3, size4, size5, size6, &
                     len(cmdname))
 
    end subroutine PPinput_size_6
 
 
! ******************************************************************************
! ******************************************************************************
! Parser setup and meta-operations
! ******************************************************************************
! ******************************************************************************
 
    ! ==========================================================================
    ! Initialize the parser
    ! ==========================================================================
    subroutine PPinput_initialize(first,deckname,oargs)
 
      character(len=*), intent(in) :: deckname,oargs
      logical         , intent(in) :: first
 
      call parser_create(oargs, len(oargs), trim(deckname), len_trim(deckname))

      call parser_comm_info(mype,numpe,iope)
      PPinput_initialized = .true.
 
    end subroutine PPinput_initialize

    ! ==========================================================================
    ! Check for other packages to make sure PPinput has been initialized
    ! ==========================================================================
    logical function PPinput_initialized_check()
       PPinput_initialized_check = .false.
       if (PPinput_initialized) PPinput_initialized_check = .true.
    end function PPinput_initialized_check

! ******************************************************************************
! ******************************************************************************
! Dictionary operations
! ******************************************************************************
! ******************************************************************************

   subroutine PPinput_dictionary_add(name, value, pred, vdesc)
      use ISO_C_BINDING, only : C_NULL_CHAR

      character(len=*), intent(in) :: name
      real(REAL64),     intent(in) :: value
      logical,          intent(in) :: pred 
      character(len=*), intent(in) :: vdesc

      integer :: c_pred = 0

      if (pred) c_pred = 1
      call parser_dictionary_add(trim(name)//C_NULL_CHAR, value, c_pred,   &
                                 trim(vdesc)//C_NULL_CHAR)

   end subroutine PPinput_dictionary_add

   subroutine PPinput_dictionary_env_add(name, pred)
      use ISO_C_BINDING, only : C_NULL_CHAR

      character(len=*), intent(in) :: name
      logical,          intent(in) :: pred 

      integer :: c_pred = 0

      if (pred) c_pred = 1
      call parser_dictionary_env_add(trim(name)//C_NULL_CHAR, c_pred)

   end subroutine PPinput_dictionary_env_add

! ******************************************************************************
! ******************************************************************************
! Compile buffer
! ******************************************************************************
! ******************************************************************************
   subroutine PPinput_compile_buffer()
      call parser_compile_buffer()
   end subroutine PPinput_compile_buffer
 
! ******************************************************************************
! ******************************************************************************
! Handle processed flags.
! ******************************************************************************
! ******************************************************************************
 
    ! ==========================================================================
    ! Check that every command/word has been processed.
    ! ==========================================================================
    subroutine PPinput_check_processed(good)
      logical, intent(inout) :: good
 
      integer :: good_int = 1
 
      if(.not.PPinput_initialized)then
         good=.false.
         return
      endif

      if (.not.  PPinput_do_check) then
         good=.true.
         return
      endif
 
      call check_processed(good_int)
 
      good = .false.
      if (good_int .eq. 1) good = .true.
 
    end subroutine PPinput_check_processed
 
 
    ! ==========================================================================
    ! Check that every command/word has been processed.
    ! We are done with the parser and could clean up here.
    ! ==========================================================================
    subroutine PPinput_terminate(good)
      logical, intent(inout) :: good
 
 
      if(.not.PPinput_initialized)then
         good = .false.
         return
      endif

      call process_error_final
      call PPinput_check_processed(good)
      call parser_destroy
    end subroutine PPinput_terminate

    
    ! ==========================================================================
    ! Disable the PPinput check operation
    ! ==========================================================================
    subroutine PPinput_disable_check()

       PPinput_do_check = .false.
       
    end subroutine PPinput_disable_check
 
 
    ! ==========================================================================
    ! Set the processed flag for all words for all commands that match cmdname.
    ! The value to set the processed flag to is bval.
    ! This sets the processed flag for commands in the final buffer and in the
    ! when...then final buffers.
    ! ==========================================================================
    subroutine PPinput_cmd_set_processed(cmdname, bval)
      character(*), intent(in) :: cmdname
      logical     , intent(in) :: bval
 
      character :: cmd_array(len(cmdname))
      integer :: bval_i, c
 
      if(.not.PPinput_initialized)return
 
      ! Pass an array of single chars rather than a character string.
      do c = 1,len(cmdname)
         cmd_array(c) = cmdname(c:c)
      enddo
 
      bval_i = 0
      if (bval) bval_i = 1
 
      call cmd_set_processed(cmd_array, len(cmdname), bval_i)
 
    end subroutine PPinput_cmd_set_processed
 
 
 
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
    subroutine PPinput_cmd_in_input(cmdname, in_input, in_whenthen)
      character(*), intent(in)    :: cmdname
      logical     , intent(inout) :: in_input, in_whenthen
 
      character :: cmd_array(len(cmdname))
      integer   :: in_input_i = 0, in_whenthen_i = 0, c
 
 
      in_input    = .false.
      in_whenthen = .false.
 
      if(.not.PPinput_initialized)return
 
      ! Pass an array of single chars rather than a character string.
      do c = 1,len(cmdname)
         cmd_array(c) = cmdname(c:c)
      enddo
 
      call cmd_in_input(cmd_array, len(cmdname), in_input_i, in_whenthen_i)
 
      if (in_input_i    .eq. 1) in_input    = .true.
      if (in_whenthen_i .eq. 1) in_whenthen = .true.
 
    end subroutine PPinput_cmd_in_input
 
 
 
! ******************************************************************************
! ******************************************************************************
! List variables, functions, etc, to the output file.
! ******************************************************************************
! ******************************************************************************
 
    ! ==========================================================================
    ! Echo user input to a file.
    ! ==========================================================================
    subroutine PPinput_echo_user_input(kw)
 
      integer, intent(in) :: kw
      integer, PARAMETER :: MAX_LENGTH=65536 ! overkill I hope
                                             ! should get line length and alloc
                                             ! but no dynamic string alloc
                                             ! until Fortran2003
      character(MAX_LENGTH) :: fline
      integer :: nchar, bint
 
      if (mype .ne. iope) return
 
      if(.not.PPinput_initialized)then
 
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
 
    end subroutine PPinput_echo_user_input
 
 
    ! ==========================================================================
    ! Echo the parser functions, variables, and final buffer to file unit
    ! number kw.
    ! ==========================================================================
    subroutine PPinput_echo_fvf(kw)
 
      integer, intent(in) :: kw
 
      if (mype .ne. iope) return
 
      if(.not.PPinput_initialized)then
         write(kw, '(/,a)') "THE INPUT DECK HAS NOT BEEN OPENED, SO ECHOING"
         write(kw, '(/,a)') "IT TO THE SCREEN CANNOT BE PERFORMED"
         return
      endif
 
      call PPinput_list_funcs  (kw)
      call PPinput_list_vars   (kw)
      call PPinput_echo_fbuffer(kw)
 
    end subroutine PPinput_echo_fvf

    ! ==========================================================================
    ! List available functions to a file.
    ! ==========================================================================
    subroutine PPinput_list_funcs(kw)
 
      integer, intent(in) :: kw
      character(256) :: fline
      integer :: nchar, bint
 
      if (mype .ne. iope) return
 
      if(.not.PPinput_initialized)then
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
 
    end subroutine PPinput_list_funcs
 
 
    ! ==========================================================================
    ! List available variables to a file.
    ! ==========================================================================
    subroutine PPinput_list_vars(kw)
 
      integer, intent(in) :: kw
      character(256) :: fline
      integer :: nchar, bint
 
      if (mype .ne. iope) return
 
      if(.not.PPinput_initialized)then
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
 
    end subroutine PPinput_list_vars
 
 
    ! ==========================================================================
    ! Echo the final buffer to a file.
    ! ==========================================================================
    subroutine PPinput_echo_fbuffer(kw)
 
      integer, intent(in) :: kw
      character(256) :: fline
      integer :: nchar, bint
 
      if (mype .ne. iope) return
 
      if(.not.PPinput_initialized)then
         write(kw, '(/,a)') "THE INPUT DECK HAS NOT BEEN OPENED, SO"
         write(kw, '(/,a)') "A PARSED LIST OF ITS CONTENTS CANNOT BE PRINTED OUT"
         return
      endif
 
      write(kw, '(/,a)') "********************************************************************************"
      write(kw, '(  a)') "********** Echo final parser buffer, this is what the code uses to set internal "
      write(kw, '(  a)') "********** code variables."
 
      call echo_final_buffer
      nchar = 256
      do
         call get_output_line(fline, nchar, bint)
         if (bint .eq. 0) exit;
         write(kw,'(a)') trim(fline)
      enddo
 
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
 
    end subroutine PPinput_echo_fbuffer
 
 
! ******************************************************************************
! ******************************************************************************
! Handle restart_block commands.
! ******************************************************************************
! ******************************************************************************
 
    ! ==========================================================================
    ! Do the restart_block checks
    ! ==========================================================================
    subroutine PPinput_rb_check(code_varnames, code_values, code_vv_active, &
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
 
      if(.not.PPinput_initialized)return
 
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
 
    end subroutine PPinput_rb_check
 
 
    ! ==========================================================================
    ! Echo restart block data to a file.
    !     kw = unit number for the file to write to
    !     rb = index of the restart block to write. Note that this is an index
    !          which starts from 0
    ! ==========================================================================
    subroutine PPinput_echo_rb_info(kw, rb)
 
      integer, intent(in) :: kw, rb
      character(256) :: fline
      integer :: nchar, bint
 
      if (mype .ne. iope) return
 
      if(.not.PPinput_initialized)then
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
 
    end subroutine PPinput_echo_rb_info
 
 
! ******************************************************************************
! ******************************************************************************
! Handle when...then commands
! ******************************************************************************
! ******************************************************************************
 
    ! ==========================================================================
    ! Get the number of when...then commands.
    ! ==========================================================================
    subroutine PPinput_whenthen_num(wtnum)
      integer, intent(out) :: wtnum
      wtnum = 0
      if(.not.PPinput_initialized)return
      call parser_get_wtnum(wtnum)
    end subroutine PPinput_whenthen_num
 
 
    ! ==========================================================================
    ! Do the when...then checks
    ! ==========================================================================
    subroutine PPinput_whenthen_check(wti, code_varnames, code_values,  &
                                     code_vv_active, nchar, nvv, wt_check)
 
      integer,                          intent(in)    :: wti, nchar, nvv
      character(nchar), dimension(nvv), intent(in)    :: code_varnames, code_values
      integer,          dimension(nvv), intent(in)    :: code_vv_active
      logical,                          intent(  out) :: wt_check
 
      character :: val_array(nchar*nvv), varname_array(nchar*nvv)
      integer :: i,c,i1d, wtci
 
      wt_check = .false.
 
      if(.not.PPinput_initialized)return
 
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
 
    end subroutine PPinput_whenthen_check
 
 
    ! ==========================================================================
    ! Set the final commands pointer to the whenthen command buffer.
    ! This is also done in the whenthen check routine if the check is true.
    ! ==========================================================================
    subroutine PPinput_whenthen_setcfp(wti)
      integer, intent(in)    :: wti
      if(.not.PPinput_initialized)return
      call parser_wt_set_cmdsfp(wti)
    end subroutine PPinput_whenthen_setcfp
 
 
    ! ==========================================================================
    ! Reset the commands final buffer pointer.
    ! ==========================================================================
    subroutine PPinput_whenthen_reset
      if(.not.PPinput_initialized)return
      call parser_wt_reset
    end subroutine PPinput_whenthen_reset
 
 
    ! ==========================================================================
    ! Get the character array and its size describing a whenthen.
    ! ==========================================================================
    subroutine PPinput_whenthen_ca(wti, wt_ca, wt_casize)
      integer, intent(in)    :: wti, wt_casize
      character(1), dimension(wt_casize), intent(inout) :: wt_ca
      wt_ca = ' '
      if(.not.PPinput_initialized)return
      call parser_wt_ca(wti, wt_ca, wt_casize)
    end subroutine PPinput_whenthen_ca
 
    subroutine PPinput_whenthen_casize(wti, wt_casize)
      integer, intent(in)    :: wti
      integer, intent(inout) :: wt_casize
      wt_casize=0
      if(.not.PPinput_initialized)return
      call parser_wt_casize(wti, wt_casize)
    end subroutine PPinput_whenthen_casize
 
 
    ! ==========================================================================
    ! Get the number of satisfied flags for a whenthen. This is equal to the
    ! number of sub-conditions in the whenthen condition.
    ! ==========================================================================
    subroutine PPinput_whenthen_satsize(wti, wt_satsize)
      integer, intent(in)    :: wti
      integer, intent(inout) :: wt_satsize
      wt_satsize=0
      if(.not.PPinput_initialized)return
      call parser_wt_satsize(wti, wt_satsize)
    end subroutine PPinput_whenthen_satsize
 
 
    ! ==========================================================================
    ! Get and Set the whenthen satisfied flags for a whenthen.
    ! ==========================================================================
    subroutine PPinput_whenthen_getsat(wti, wt_sat, wt_satsize)
      integer, intent(in)    :: wti, wt_satsize
      integer, dimension(wt_satsize), intent(inout) :: wt_sat
      wt_sat=0
      if(.not.PPinput_initialized)return
      call parser_wt_getsat(wti, wt_sat, wt_satsize)
    end subroutine PPinput_whenthen_getsat
 
    subroutine PPinput_whenthen_setsat(wti, wt_sat, wt_satsize)
      integer, intent(in)    :: wti, wt_satsize
      integer, dimension(wt_satsize), intent(inout) :: wt_sat
      wt_sat=0
      if(.not.PPinput_initialized)return
      call parser_wt_setsat(wti, wt_sat, wt_satsize)
    end subroutine PPinput_whenthen_setsat
 
 
    ! ==========================================================================
    ! Get and Set the processed flag for a whenthen.
    ! ==========================================================================
    subroutine PPinput_whenthen_getprocessed(wti, wtp)
      integer, intent(in)    :: wti
      integer, intent(inout) :: wtp
      wtp=0
      if(.not.PPinput_initialized)return
      call parser_wt_getprocessed(wti, wtp)
    end subroutine PPinput_whenthen_getprocessed
 
    subroutine PPinput_whenthen_setprocessed(wti, wtp)
      integer, intent(in) :: wti
      integer, intent(in) :: wtp
      if(.not.PPinput_initialized)return
      call parser_wt_setprocessed(wti, wtp)
    end subroutine PPinput_whenthen_setprocessed
 
 
    ! ==========================================================================
    ! Get and Set the sequence number for a whenthen.
    ! ==========================================================================
    subroutine PPinput_whenthen_getseq(wti, wtseq)
      integer, intent(in)    :: wti
      integer, intent(inout) :: wtseq
      wtseq=0
      if(.not.PPinput_initialized)return
      call parser_wt_getseq(wti, wtseq)
    end subroutine PPinput_whenthen_getseq
 
    subroutine PPinput_whenthen_setseq(wti, wtseq)
      integer, intent(in) :: wti
      integer, intent(in) :: wtseq
      if(.not.PPinput_initialized)return
      call parser_wt_setseq(wti, wtseq)
    end subroutine PPinput_whenthen_setseq
 
 
 
! ------------------------------------------------------------------------------
  end module PPinput_module
! ==============================================================================
