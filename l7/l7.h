/*
 *  Copyright (c) 2011-2012, Los Alamos National Security, LLC.
 *  All rights Reserved.
 *
 *  Copyright 2011-2012. Los Alamos National Security, LLC. This software was produced 
 *  under U.S. Government contract DE-AC52-06NA25396 for Los Alamos National 
 *  Laboratory (LANL), which is operated by Los Alamos National Security, LLC 
 *  for the U.S. Department of Energy. The U.S. Government has rights to use, 
 *  reproduce, and distribute this software.  NEITHER THE GOVERNMENT NOR LOS 
 *  ALAMOS NATIONAL SECURITY, LLC MAKES ANY WARRANTY, EXPRESS OR IMPLIED, OR 
 *  ASSUMES ANY LIABILITY FOR THE USE OF THIS SOFTWARE.  If software is modified
 *  to produce derivative works, such modified software should be clearly marked,
 *  so as not to confuse it with the version available from LANL.
 *
 *  Additionally, redistribution and use in source and binary forms, with or without
 *  modification, are permitted provided that the following conditions are met:
 *     * Redistributions of source code must retain the above copyright
 *       notice, this list of conditions and the following disclaimer.
 *     * Redistributions in binary form must reproduce the above copyright
 *       notice, this list of conditions and the following disclaimer in the
 *       documentation and/or other materials provided with the distribution.
 *     * Neither the name of the Los Alamos National Security, LLC, Los Alamos 
 *       National Laboratory, LANL, the U.S. Government, nor the names of its 
 *       contributors may be used to endorse or promote products derived from 
 *       this software without specific prior written permission.
 *  
 *  THIS SOFTWARE IS PROVIDED BY THE LOS ALAMOS NATIONAL SECURITY, LLC AND 
 *  CONTRIBUTORS "AS IS" AND ANY EXPRESS OR IMPLIED WARRANTIES, INCLUDING, BUT 
 *  NOT LIMITED TO, THE IMPLIED WARRANTIES OF MERCHANTABILITY AND FITNESS FOR
 *  A PARTICULAR PURPOSE ARE DISCLAIMED. IN NO EVENT SHALL LOS ALAMOS NATIONAL
 *  SECURITY, LLC OR CONTRIBUTORS BE LIABLE FOR ANY DIRECT, INDIRECT, INCIDENTAL,
 *  SPECIAL, EXEMPLARY, OR CONSEQUENTIAL DAMAGES (INCLUDING, BUT NOT LIMITED TO,
 *  PROCUREMENT OF SUBSTITUTE GOODS OR SERVICES; LOSS OF USE, DATA, OR PROFITS;
 *  OR BUSINESS INTERRUPTION) HOWEVER CAUSED AND ON ANY THEORY OF LIABILITY,
 *  WHETHER IN CONTRACT, STRICT LIABILITY, OR TORT (INCLUDING NEGLIGENCE OR OTHERWISE)
 *  ARISING IN ANY WAY OUT OF THE USE OF THIS SOFTWARE, EVEN IF ADVISED OF THE
 *  POSSIBILITY OF SUCH DAMAGE.
 */  
#ifndef L7_H_
#define L7_H_
#ifdef HAVE_OPENCL
#include "ezcl/ezcl.h"
#endif

//#define _L7_DEBUG

#ifdef __cplusplus
extern "C"
{
#endif
	
/*
 * Some L7 parameters.
 */

#define L7_OK   0 /* Successful return. */

enum  L7_Datatype
{
   L7_GENERIC8  = 0,
   L7_BYTE,
   L7_PACKED,
   
   L7_CHAR,
   L7_INT,
   L7_LONG,
   L7_LONG_LONG_INT,
   L7_FLOAT,
   L7_DOUBLE,
   
   L7_CHARACTER,
   L7_LOGICAL,
   L7_INTEGER4,
   L7_INTEGER8,
   L7_REAL4,
   L7_REAL8,

   L7_DATATYPE_MIN = L7_GENERIC8,
   L7_DATATYPE_MAX = L7_REAL8
};

/* Processor/disk patterns */
enum L7_DiskPatternType
{
   L7_ONE_PROC           = 101,
   L7_ALL_PROCS,
   L7_BUFFERED_ONE_PROC,
   L7_BUFFERED_ALL_PROCS,
   
   L7_DISK_PATTERN_MIN = L7_ONE_PROC,
   L7_DISK_PATTERN_MAX = L7_BUFFERED_ALL_PROCS
};

/* Data patterns
 * Starting these at 1001 so they don't get crossed
 * up with the disk patterns
 */
enum L7_DataPatternType
{
   L7_REPLICATED  = 1001,
   L7_DISTRIBUTED,
   
   L7_DATA_PATTERN_MIN = L7_REPLICATED,
   L7_DATA_PATTERN_MAX = L7_DISTRIBUTED   
};

/* IO Profiling Verbosity Levels */
enum L7_IO_ProfilingLevel
{
   L7_IO_PROF_OFF     = 0,
   L7_IO_PROF_SIMPLE,
   L7_IO_PROF_VERBOSE,

   L7_IO_PROF_LEVEL_MIN = L7_IO_PROF_OFF,
   L7_IO_PROF_LEVEL_MAX = L7_IO_PROF_VERBOSE   
};

/*
 * C Prototypes.
 * 
 * Note that scalar variables that are input are usually made
 * const. Large input arrays are not made const to avoid copying
 * by the compiler.  
 */

int L7_Init (
	int *mype,
	int *numpes,
	int *argc,
	char **argv,
        int do_quo_setup,
        int lttrace_on
	);

void L7_Dev_Init(void);

int L7_Terminate(void);

int L7_Broadcast (
		void                   *data_buffer,
		const int              count,
		const enum L7_Datatype l7_datatype,
		const int              root_pe
		);

void L7_Abort (void); 

void L7_Sum(
		void                    *input,
		const int               count,
		const enum L7_Datatype  l7_datatype,
		void                    *output
		);

int L7_Max(
		void                    *input,
		const int               count,
		const enum L7_Datatype  l7_datatype,
		void                    *output
		);

int L7_Min(
		void                    *input,
		const int               count,
		const enum L7_Datatype  l7_datatype,
		void                    *output
		);

int  L7_Any(
		void                    *input,
		const int               count,
		const enum L7_Datatype  l7_datatype,
		void                    *output
		);

int L7_All(
		void                    *input,
		const int               count,
		const enum L7_Datatype  l7_datatype,
		void                    *output
		);

int L7_Array_Sum(
		void                    *input,
		const int               count,
		const enum L7_Datatype  l7_datatype,
		void                    *output
		);

int L7_Array_Max(
		void                    *input,
		const int               count,
		const enum L7_Datatype  l7_datatype,
		void                    *output
		);

int L7_Array_Min(
		void                    *input,
		const int               count,
		const enum L7_Datatype  l7_datatype,
		void                    *output
		);

int L7_MaxLoc(
		void                    *input,
		const int               count,
		const enum L7_Datatype  l7_datatype,
		int                     *output
		);

int L7_MinLoc(
      void                    *input,
      const int               count,
      const enum L7_Datatype  l7_datatype,
      int                     *output
      );

int L7_MaxValLoc(
		void                    *input,
		const int               count,
		const enum L7_Datatype  l7_datatype,
		void                    *val,
		int                     *loc
		);

int L7_MinValLoc(
      void                    *input,
      const int               count,
      const enum L7_Datatype  l7_datatype,
      void                    *val,
      int                     *loc
      );

/*
int L7_MaxValLocLoc(
		void       *input,
		void       *loc2input,
		const int  count,
		const int  type,
		void       *val,
		void       *loc,
		void       *loc2
		);
*/

int L7_MaxValLocLoc_Int4(
		void                    *input,
		void                    *loc2input,
		const int               count,
		const enum L7_Datatype  l7_datatype,
		void                    *val,
		void                    *loc,
		void                    *loc2
		);

int L7_MaxValLocLoc_Int8(
		void                    *input,
		void                    *loc2input,
		const int               count,
		const enum L7_Datatype  l7_datatype,
		void                    *val,
		void                    *loc,
		void                    *loc2
		);

void L7_Global_Offset(
                const int count,
                int *output
                );

int L7_GetGlobal(
		void                    *array,
		const int               count,
		const enum L7_Datatype  l7_datatype,
		const int               index,
		void                    *val
		);

int L7_Scatter(
		void                    *global_buffer,
		void                    *local_buffer,
		const int               count,
		const enum L7_Datatype  l7_datatype,
		const int               root_pe
		);

int L7_Allscatter(
		void                    *global_buffer,
		void                    *local_buffer,
		const int               count,
		const enum L7_Datatype  l7_datatype
		);

int L7_AllScatter_Sizes_Known(
		void                    *global_buffer,
		void                    *local_buffer,
		const int               *nsizearray,
		const enum L7_Datatype  l7_datatype
		);

int L7_Gather(
		void                    *local_buffer,
		const int               count,
		void                    *global_buffer,
		const enum L7_Datatype  l7_datatype,
		const int               root_pe
		);

int L7_Allgather(
		void                    *local_buffer,
		const int               count,
		void                    *global_buffer,
		const enum L7_Datatype  l7_datatype
		);

int L7_Get_Rank(void);

int L7_Get_Numpes(void);

double L7_Wtime(void);

int L7_Get_Timeout_Signal(void);

int L7_Setup(
        const int       num_base,
        const int       my_start_index,
        const int       num_indices_owned,
        int             *indices_needed,
        const int       num_indices_needed,
        /* Remove for now until can be made portable and is needed */
        //const MPI_Comm  *mpi_comm_user, 
        int             *l7_id
        );

int L7_Dev_Setup(
        const int       num_base,
        const int       my_start_index,
        const int       num_indices_owned,
        int             *indices_needed,
        const int       num_indices_needed,
        /* Remove for now until can be made portable and is needed */
        //const MPI_Comm  *mpi_comm_user, 
        int             *l7_id
        );

int L7_Free(
        const int   *l7_id
        );

int L7_Dev_Free(
        void
        );

int L7_Update(
      void                    *data_buffer,
      const enum L7_Datatype  l7_datatype,
      const int               l7_id
      );

#ifdef HAVE_OPENCL
int L7_Dev_Update(
      cl_mem                  dev_data_buffer,
      const enum L7_Datatype  l7_datatype,
      const int               l7_id
      );
#endif

int L7_Update_Check(
      void                    *data_buffer,
      const enum L7_Datatype  l7_datatype,
      const int               l7_id
      );

int L7_Get_Num_Indices(
      const int               l7_id
      );

int L7_Get_Local_Indices(
      const int               l7_id,
      int                     *local_indices
      );

int L7_Push_Setup(
      const int               num_comm_partners,
      const int               *comm_partner,
      const int               *send_buffer_count,
      int                     **send_database, 
      int                     *receive_count_total,
      int                     *l7_push_id
      );

int L7_Push_Update(
      const int               *array,
      int                     *return_array,
      const int               l7_push_id
      );

int L7_Push_Free(
      const int               *l7_push_id
      );

/*
 * L7 File Prototypes.
 */

int L7_File_Open(
      const char *fdesc,
      const char *type,
      const enum L7_DiskPatternType  l7_disk_pattern
      );

int L7_File_Close(
      const int fid
      );

void L7_File_Read(
      const int                      fid,
      long long                      *disk_loc,
      void                           *buf,
      const long long                nwords,
      const enum L7_Datatype         l7_datatype,
      const enum L7_DataPatternType  l7_data_pattern
      );

void L7_File_Write(
      const int                      fid,
      long long                      *disk_loc,
      void                           *buf,
      const long long                nwords,
      const enum L7_Datatype         l7_datatype,
      const enum L7_DataPatternType  l7_data_pattern
      );

int L7_File_Inquire(
      const char                     *fildes,
      const enum L7_DiskPatternType  l7_disk_pattern
      );

int L7_File_Unlink(
      const char *name,
      const enum L7_DiskPatternType  l7_disk_pattern
      );

long long L7_File_Size(
      const char *fildes,
      const enum L7_DiskPatternType  l7_disk_pattern
      );

int L7_File_Link(
      const char *oldname,
      const char *newname,
      const enum L7_DiskPatternType  l7_disk_pattern
      );

int L7_File_Unlink(
      const char *name,
      const enum L7_DiskPatternType  l7_disk_pattern
      );

int L7_File_Symlink(
      const char *oldname,
      const char *newname,
      const enum L7_DiskPatternType  l7_disk_pattern
      );

void L7_Sort(
      void                   *array_in,
      const int              nsize,
      const enum L7_Datatype l7_datatype
      );

void L7_Sort_2Arrays(
      void *                 array_in1,
      void *                 array_in2,
      const int              nsize,
      const enum L7_Datatype l7_datatype
      );

void L7_Index_Sort(
      void *                 array_in,
      const int              nsize,
      const enum L7_Datatype l7_datatype,
      int *                  index
      );

void L7_Indexi8_Sort(
      void *                 array_in,
      const int              nsize,
      const enum L7_Datatype l7_datatype,
      long *                 index
      );


long long L7_Address(void *var);
void L7_Index_sort_real8(const int n,double array_in[],int index[]);
void L7_Index_sort_int8(const int n,long long iarray_in[], int index[]);
void L7_Index_sort_int4(const int n, int iarray_in[], int index[]);
void L7_Index_sort_real8_int8(const int n,double array_in[],long long index[]);
void L7_INDEX_SORT_REAL8_INT8_F(const int *n, double *array_in, long long *index);
void L7_Index_sort_int8_int8(const int n,long long iarray_in[], long long index[]);
void L7_Index_sort_int4_int8(const int n, int iarray_in[], long long index[]);
void L7_Sort_real8(const int n,double array_in[]);
void L7_Sort_int8(const int n,long long array_in[]);
void L7_Sort_int4(const int n,int array_in[]);
void L7_Sort_real8_real8(const int n,double array_in[],double array_in2[]);
void L7_Sort_int8_int8(const int n,long long array_in[],long long array_in2[]);
void L7_Sort_int4_int4(const int n,int array_in[],int array_in2[]);














/*
 * End prototypes.
 */

/*
 * remove typesafe linkage if compiling under c++
 */

#ifdef __cplusplus
}
#endif

#endif /* L7_H */
