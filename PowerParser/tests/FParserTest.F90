program FParserTest

   use FParser_module

   implicit none


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
   integer :: maxpe

   real(REAL64), parameter :: PI        = 3.14159265358979323846_REAL64
   ! Speed of light in a vacuum (exact)
   real(REAL64), parameter :: SOL           = 2.997924580E+10_REAL64 ! [cm/s]

   character(24) :: Package = " Package:: Parser:"
   character(24) :: Module  = " Module:: QueryFParser:"
   character(120) :: msg

   character(20) :: deckname_arg = "parsetest.in"
   character(20) :: oargs = " "
   integer :: check_input = 0

   character(32) :: special_version
   integer       :: modver, modext

   ! local variables
   integer :: i, n, exearg_cmd01

   logical :: some_logical_cmd
   real(REAL8) :: exp_val, denom, rdiff, volume_cmd

   logical :: sp_logical_array(10), mult_logical_array(6)
   integer :: math_result1, math_result3, math_result4, math_result5
   integer :: math_result5B
   integer :: math_result6, math_result7, math_result8, math_result9
   integer :: math_result10, upm01, upm02
   real(REAL8) :: math_result2, math_result11, math_result12
   integer :: int_array(4), pint_nested(10)
   real(REAL8) :: negnum, xcenter, cmdml(5), f01, f02
   real(REAL8) :: acmd, acmd2, acmd3(5), acmd4, acmd5(8)
   real(REAL8) :: depcmd01(5)
   real(REAL8), dimension(:), allocatable :: acontline
   integer, dimension(:), allocatable :: skip1d
   integer :: acontline_size, skip1d_size
   integer :: oddstrs_size
   logical :: good, idefvar_cmd01, idefvar_cmd02
   real(REAL8) :: idefvar_cmd03, idefvar_cmd04

   real(REAL8) :: skip_check_cmd(2)

   real(REAL8) :: var1d_res, var2d_res, vnc_cmd
   integer :: var8d_cmd, var8d_cmd2
   logical :: log1d_cmd
   character(24) :: vchar3d_cmd

   logical :: in_input, in_whenthen

   logical :: math_result13, math_result14, math_result15
   logical :: math_result16, math_result17, math_result18

   logical, dimension(:),       allocatable :: logical_array
   logical, dimension(:,:),     allocatable :: log2d
   logical, dimension(:,:,:),   allocatable :: log3d
   logical, dimension(:,:,:,:), allocatable :: log4d
   integer :: logical_array_size, log2d_size, log3d_size, log4d_size

   real(REAL8)                                  :: a2d(3,2)
   real(REAL8), dimension(:,:,:), allocatable   :: a3d
   real(REAL8), dimension(:,:,:,:), allocatable :: a4d
   integer :: a2d_size, a3d_size, a4d_size

   integer, dimension(:,:),     allocatable :: i2d
   integer, dimension(:,:,:),   allocatable :: i3d
   integer, dimension(:,:,:,:), allocatable :: i4d
   integer :: i2d_size, i3d_size, i4d_size

   character :: single_char, single_charq
   character(24) :: title, char1d(6), c1d_mult(5)
   character(24), dimension(:), allocatable :: oddstrs
   character(24), dimension(:,:), allocatable     :: char2d
   character(24), dimension(:,:,:), allocatable   :: char3d
   character(24), dimension(:,:,:,:), allocatable :: char4d
   integer :: char2d_size, char3d_size, char4d_size

   real(REAL8) :: delta_y_cmd01, delta_y_cmd02, delta_y_cmd03
   integer :: delta_y_cmd04
   real(REAL8) :: delta_y_cmd05, delta_y_cmd06, delta_y_cmd07
   real(REAL8) :: delta_x_cmd01, delta_x_cmd02, delta_x_cmd03

   integer :: do_sum_cmd01, do_sum_cmd02, do_sum_cmd03, do_sum_cmd04
   integer :: do_sum_cmd05, do_sum_cmd06, do_sum_cmd07

   integer :: sub_cmd01, sub_cmd02, sub_cmd03, sub_cmd04, sub_cmd05
   integer :: sub_cmd06(4), sub_cmd07, sub_cmd08, sub_cmd09
   integer :: sub_cmd10

   real(REAL8) :: inc_cmd01, inc_cmd02, inc_cmd03, inc_cmd04

   real(REAL8), dimension(:), allocatable     :: asm1d
   real(REAL8), dimension(:,:), allocatable   :: asm2d, mults
   integer :: asm1d_size, asm2d_size
   integer :: mults_size1, mults_size2

   integer :: iarith_cmd01, iarith_cmd02, iarith_cmd03

   integer :: strlen_cmd01, strlen_cmd02
   character(24) :: strcat_cmd01, strcat_cmd02, strerase_cmd01
   character(24) :: strerase_cmd02, strinsert_cmd01
   character(24) :: strsubstr_cmd01, strtrim_cmd01

   real(REAL8) :: matdef(200, 99)
   integer :: matreg(99), matreg2(99), nummat, numreg, numreg2

   integer :: ppmm_cmd01, ppmm_cmd02, ppmm_cmd03, ppmm_cmd04
   integer :: ppmm_cmd05, ppmm_cmd06, ppmm_cmd07, ppmm_cmd08
   integer :: ppmm_cmd09, ppmm_cmd10, ppmm_cmd11
   real(REAL8) :: ppmm_cmd12, ppmm_cmd13
   integer :: ppmm_cmd14
   real(REAL8) :: ppmm_cmd15, ppmm_cmd16

   logical :: rb_check
   integer :: rb_ntriggered, rb_num
   integer, dimension(:), allocatable :: rb_triggered_indices

   integer :: shortmodcyc, wtnum, modcyc, ncycle, wt_cmd04, wt_cmd05
   integer :: num_wt_cyc, wt_cmd06
   logical :: wt_check, wttf_c01, wt_cmd01, wt_cmd03
   real(REAL8) :: sim_time, wt_cmd02, sim_pressure(5)
   character(24) :: wttf_c02
   character(24), dimension(5) :: code_varnames, code_values
   integer, dimension(5) :: code_vv_active
   integer :: wt_casize = 0, max_casize
   character(1), dimension(:), allocatable :: wt_ca
   integer :: wt_satsize = 0, max_satsize
   integer, dimension(:), allocatable :: wt_sat

   integer :: quad_root1, quad_root2

   ! Passing C strings to Fortran90 is fraught with danger, so pass
   ! them in as character arrays and copy the contents into Fortran
   ! strings

   ! Initialize the QueryFParser user input system - for getting user input.
   !call QueryFParsera(fname, .false.)
   call FParser_initialize(.true.,deckname_arg,oargs)

   ! Add dictionary entries that are "pre-defined" -- third argument is
   ! true and called before the compile_buffer operations
   call FParser_dictionary_add("$sol",SOL,.true.,"Speed of light (cm/s)")
   ! PI is defined in define_kind
   call FParser_dictionary_add("$pi", PI, .true.,"Circle circumference/diameter")

   call FParser_compile_buffer()

   ! Get processor info.
   !     mype  = id for this processor (0 to numpe-1)
   !     iope  = id for the i/o processor (0 to numpe-1 -- normally 0)
   !     numpe = number of processors
   mype = 0
   iope = 0
   numpe = 1
   call parser_comm_info(mype, numpe, iope)
   maxpe = numpe-1

!  call set_version_info()

   ! Open the output file.
   open(12, FILE="test_output", STATUS='REPLACE')

   ! Check for duplicate scalar values.
   ! Echo the input file, vars, functions, final buffer.
   call parser_chk_scalar_dup
   call FParser_echo_user_input(12)
  !call FParser_echo_fvf(12)

!  modver = v_modver
!  modext = v_modext
!  special_version = v_special_version

   ! Initialize the unit test package.

   ! ---------------------------------------------------------------------------
   ! ---------------------------------------------------------------------------
   ! Logical data.

   ! Test simple logical command.
   some_logical_cmd = .false.
   call QueryFParser('some_logical_cmd', some_logical_cmd, .true.)
   !!!!!!@assertTrue(some_logical_cmd, trim(Package)//trim(Module)//" simple logical command")

   ! ******************************************************************************
   ! ******************************************************************************
   ! Parser - Check Processed
   ! At this point, all the parsing has been done, check that everything has been
   ! processed. We could also terminate the parser, i.e. free memory, clean up.
   ! ******************************************************************************
   ! ******************************************************************************
   !call FParser_terminate(good)

   call parser_destroy

   ! ******************************************************************************
   ! ******************************************************************************
   ! Parser - Final Summary
   ! ******************************************************************************
   ! ******************************************************************************

   ! Close the output file.
   close(12)

end program FParserTest
