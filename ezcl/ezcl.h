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
 *  
 *  CLAMR -- LA-CC-11-094
 *  This research code is being developed as part of the 
 *  2011 X Division Summer Workshop for the express purpose
 *  of a collaborative code for development of ideas in
 *  the implementation of AMR codes for Exascale platforms
 *  
 *  AMR implementation of the Wave code previously developed
 *  as a demonstration code for regular grids on Exascale platforms
 *  as part of the Supercomputing Challenge and Los Alamos 
 *  National Laboratory
 *  
 *  Authors: Bob Robey       XCP-2   brobey@lanl.gov
 *           Neal Davis              davis68@lanl.gov, davis68@illinois.edu
 *           David Nicholaeff        dnic@lanl.gov, mtrxknight@aol.com
 *           Dennis Trujillo         dptrujillo@lanl.gov, dptru10@gmail.com
 * 
 */
#define HAVE_OPENCL 1

#ifdef HAVE_CONFIG_H
#include "config.h"
#endif

#ifdef HAVE_OPENCL

#ifdef __APPLE_CC__
#include <OpenCL/OpenCL.h>
#else
#include <CL/cl.h>
#endif

#ifdef HAVE_CL_DOUBLE
typedef cl_double cl_real;
typedef cl_double4 cl_real4;
#else
typedef cl_float cl_real;
typedef cl_float4 cl_real4;
#endif

#define PROGRAM_OBJECT       1
#define KERNEL_OBJECT        2
#define COMMAND_QUEUE_OBJECT 3
#define CONTEXT_OBJECT       4
#define EVENT_OBJECT         5

#define EZCL_SUCCESS         0
#define EZCL_NODEVICE        1

#define EZCL_PINNED_MEMORY   1

#ifdef __cplusplus
extern "C"
{
#endif

// Profiling interfaces -- file and line number are inserted at these interfaces to give users
// better error messages
/* init and end routines */
#define ezcl_init(  ezcl_gpu_context, ezcl_cpu_context, ezcl_accelerator_context) \
      ( ezcl_init_p(ezcl_gpu_context, ezcl_cpu_context, ezcl_accelerator_context, __FILE__, __LINE__) )
#define ezcl_devtype_init(  device_type, return_context, return_command_queue, mype) \
      ( ezcl_devtype_init_p(device_type, return_context, return_command_queue, mype, __FILE__, __LINE__) )
#define ezcl_finalize() \
      ( ezcl_finalize_p(__FILE__,__LINE__) )

/* device based routines */
#define ezcl_get_device(  context) \
      ( ezcl_get_device_p(context, __FILE__, __LINE__) ) 
#define ezcl_device_info(  device) \
      ( ezcl_device_info_p(device, __FILE__, __LINE__) )

/* context based routines*/
#define ezcl_create_command_queue(  context, mype) \
      ( ezcl_create_command_queue_p(context, mype, __FILE__, __LINE__) )
#define ezcl_get_command_queue() \
      ( ezcl_get_command_queue_p(__FILE__, __LINE__) )
#define ezcl_get_context() \
      ( ezcl_get_context_p(__FILE__, __LINE__) )

/* memory routines */
#define ezcl_malloc(  host_mem_ptr, dims, elsize, flags, ezcl_flags) \
      ( ezcl_malloc_p(host_mem_ptr, dims, elsize, flags, ezcl_flags, __FILE__, __LINE__) )
#define ezcl_malloc_memory_add(  malloc_mem_ptr, size) \
      ( ezcl_malloc_memory_add_p(malloc_mem_ptr, size, __FILE__, __LINE__) )
#define ezcl_device_memory_add(  dev_mem_ptr, size) \
      ( ezcl_device_memory_add_p(dev_mem_ptr, size, __FILE__, __LINE__) )
#define ezcl_mapped_memory_add(  map_mem_ptr, size) \
      ( ezcl_mapped_memory_add_p(map_mem_ptr, size, __FILE__, __LINE__) )
#define ezcl_device_memory_remove(  dev_mem_ptr) \
      ( ezcl_device_memory_remove_p(dev_mem_ptr, __FILE__, __LINE__) )
#define ezcl_mapped_memory_remove(  map_mem_ptr) \
      ( ezcl_mapped_memory_remove_p(map_mem_ptr, __FILE__, __LINE__) )
#define ezcl_malloc_memory_remove(  malloc_mem_ptr) \
      ( ezcl_malloc_memory_remove_p(malloc_mem_ptr, __FILE__, __LINE__) )
#define ezcl_mem_free_all() \
      ( ezcl_mem_free_all_p(__FILE__, __LINE__) ) 
#define genvector(  inum, elsize) \
      ( genvector_p(inum, elsize, __FILE__, __LINE__) )
#define genvectorfree(  var) \
      ( genvectorfree_p(var, __FILE__, __LINE__) )
#define genmatrix(  jnum, inum, elsize) \
      ( genmatrix_p(jnum, inum, elsize, __FILE__, __LINE__) )
#define gentrimatrix(  knum, jnum, inum, elsize) \
      ( gentrimatrix_p(knum, jnum, inum, elsize, __FILE__, __LINE__) )
#define genmatrixfree(  var) \
      ( genmatrixfree_p(var, __FILE__, __LINE__) )
#define gentrimatrixfree(  var) \
      ( gentrimatrixfree_p(var, __FILE__, __LINE__) )

/* kernel and program routines */
#define ezcl_create_kernel(  context, filename, kernel_name, flags) \
      ( ezcl_create_kernel_p(context, filename, kernel_name, flags, __FILE__, __LINE__) )

#define ezcl_program_release(  program) \
      ( ezcl_program_release_p(program, __FILE__, __LINE__) )
#define ezcl_kernel_release(  kernel) \
      ( ezcl_kernel_release_p(kernel, __FILE__, __LINE__) )
#define ezcl_command_queue_release(  command_queue) \
      ( ezcl_command_queue_release_p(command_queue, __FILE__, __LINE__) )
#define ezcl_context_release(  context) \
      ( ezcl_context_release_p(context, __FILE__, __LINE__) )

#define ezcl_enqueue_write_buffer(  command_queue, mem_buffer, blocking_flag, offset, size, ptr, event) \
      ( ezcl_enqueue_write_buffer_p(command_queue, mem_buffer, blocking_flag, offset, size, ptr, event, __FILE__, __LINE__) )
#define ezcl_enqueue_read_buffer(  command_queue, mem_buffer, blocking_flag, offset, size, ptr, event) \
      ( ezcl_enqueue_read_buffer_p(command_queue, mem_buffer, blocking_flag, offset, size, ptr, event, __FILE__, __LINE__) )
#define ezcl_enqueue_ndrange_kernel(  command_queue, kernel, work_dim, global_work_offset, global_work_size, local_work_size, event) \
      ( ezcl_enqueue_ndrange_kernel_p(command_queue, kernel, work_dim, global_work_offset, global_work_size, local_work_size, event, __FILE__, __LINE__) )
#define ezcl_set_kernel_arg(  kernel, arg_index, arg_size, arg_value) \
      ( ezcl_set_kernel_arg_p(kernel, arg_index, arg_size, arg_value, __FILE__, __LINE__) )
#define ezcl_flush(  command_queue) \
      ( ezcl_flush_p(command_queue, __FILE__, __LINE__) )
#define ezcl_finish(  command_queue) \
      ( ezcl_finish_p(command_queue, __FILE__, __LINE__) )
#define ezcl_get_event_profiling_info(  event, param_name, param_value_size, param_value) \
      ( ezcl_get_event_profiling_info_p(event, param_name, param_value_size, param_value, __FILE__, __LINE__) )

#define ezcl_timer_calc(  start_read_event, end_read_event) \
      ( ezcl_timer_calc_p(start_read_event, end_read_event, __FILE__, __LINE__) )
   

/* init and finish routines -- use either ezcl_init or ezcl_devtype_init */
/*   the devtype_init is the simpler and recommended path */
cl_int ezcl_init_p(cl_context *ezcl_gpu_context, cl_context *ezcl_cpu_context, cl_context *ezcl_accelerator_context, const char *file, const int line);
cl_int ezcl_devtype_init_p(cl_device_type device_type, cl_context *return_context, cl_command_queue *return_command_queue, const int mype, const char *file, const int line);
cl_int ezcl_finalize_p(const char *file, const int line);

/* Error reporting */
void ezcl_print_error(const int ierr, const char *routine, const char *cl_routine, const char *file, const int line);

/* device based routines */
cl_device_id ezcl_get_device_p(cl_context context, const char *file, const int line);
void ezcl_device_info_p(cl_device_id device, const char *file, const int line);

/* context based routines*/
cl_command_queue ezcl_create_command_queue_p(cl_context context, const int mype, const char *file, const int line);
cl_command_queue ezcl_get_command_queue_p(const char *file, const int line);
cl_context ezcl_get_context_p(const char *file, const int line);

/* memory routines */
cl_mem ezcl_malloc_p(void *host_mem_ptr, size_t dims[], size_t elsize, size_t flags, int ezcl_flags, const char *file, const int line);
void *ezcl_malloc_memory_add_p(void *malloc_mem_ptr, size_t size, const char *file, const int line);
void ezcl_device_memory_add_p(cl_mem dev_mem_ptr, size_t size, const char *file, const int line);
void ezcl_mapped_memory_add_p(cl_mem map_mem_ptr, size_t size, const char *file, const int line);
void ezcl_device_memory_remove_p(void *dev_mem_ptr, const char *file, const int line);   
void ezcl_mapped_memory_remove_p(void *map_mem_ptr, const char *file, const int line);   
void ezcl_malloc_memory_remove_p(void *malloc_mem_ptr, const char *file, const int line);   
void ezcl_mem_free_all_p(const char *file, const int line);
void *genvector_p(int inum, size_t elsize, const char *file, const int line);
void genvectorfree_p(void *var, const char *file, const int line);
void **genmatrix_p(int jnum, int inum, size_t elsize, const char *file, const int line);
void ***gentrimatrix_p(int knum, int jnum, int inum, size_t elsize, const char *file, const int line);
void genmatrixfree_p(void **var, const char *file, const int line);
void gentrimatrixfree_p(void ***var, const char *file, const int line);

/* kernel and program routines */
cl_kernel ezcl_create_kernel_p(cl_context context, const char *filename, const char *kernel_name, int flags, const char *file, const int line);
//cl_command_queue ezcl_create_command_queue(void);

void ezcl_program_release_p(cl_program program, const char *file, const int line);
void ezcl_kernel_release_p(cl_kernel kernel, const char *file, const int line);
void ezcl_command_queue_release_p(cl_command_queue command_queue, const char *file, const int line);
void ezcl_context_release_p(cl_context context, const char *file, const int line);

void ezcl_enqueue_write_buffer_p(cl_command_queue command_queue, cl_mem mem_buffer, cl_bool blocking_flag,
                               size_t offset, size_t size, const void *ptr, cl_event *event, const char *file, const int line);
void ezcl_enqueue_read_buffer_p(cl_command_queue command_queue, cl_mem mem_buffer, cl_bool blocking_flag,
                              size_t offset, size_t size, void *ptr, cl_event *event, const char *file, const int line);
void ezcl_enqueue_ndrange_kernel_p(cl_command_queue command_queue, cl_kernel kernel, cl_uint work_dim,
                                 const size_t *global_work_offset, const size_t *global_work_size, const size_t *local_work_size, cl_event *event, const char *file, const int line);
void ezcl_set_kernel_arg_p(cl_kernel kernel, cl_uint arg_index, size_t arg_size, const void *arg_value, const char *file, const int line);
void ezcl_flush_p(cl_command_queue command_queue, const char *file, const int line);
void ezcl_finish_p(cl_command_queue command_queue, const char *file, const int line);
void ezcl_get_event_profiling_info_p(cl_event event, cl_profiling_info param_name, size_t param_value_size, void *param_value, const char *file, const int line);

void ezcl_set_timing_flag(void);
long ezcl_timer_calc_p(cl_event *start_read_event, cl_event *end_read_event, const char *file, const int line);
   
#ifdef __cplusplus
}
#endif

#endif
      
