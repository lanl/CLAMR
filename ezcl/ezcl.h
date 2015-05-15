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
#ifndef _EZCL_H
#define _EZCL_H

#ifdef HAVE_CONFIG_H
#include "config.h"
#endif

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

#define EZCL_REGULAR_MEMORY  0x0001
#define EZCL_MANAGED_MEMORY  0x0002

#define COMPUTE_DEVICE_NVIDIA   1
#define COMPUTE_DEVICE_ATI      2
#define COMPUTE_DEVICE_INTEL    3

#ifdef __cplusplus
extern "C"
{
#endif

// Profiling interfaces -- file and line number are inserted at these interfaces to give users
// better error messages
/* init and end routines */
#define ezcl_init(  ezcl_gpu_context, ezcl_cpu_context, ezcl_accelerator_context) \
      ( ezcl_init_p(ezcl_gpu_context, ezcl_cpu_context, ezcl_accelerator_context, __FILE__, __LINE__) )

/****************************************************************//**
 * \brief
 * Detects and sets up an OpenCL device of the specified type
 *
 * **Parameters**
 * * cl_device_type device type -- device type specified can be
 *     CL_DEVICE_TYPE_GPU
 *     CL_DEVICE_TYPE_ACCELERATOR
 *     CL_DEVICE_TYPE_CPU
 * Returns an int with the possible values:
 *   EZCL_NODEVICE
 *   EZCL_SUCCESS
 *
 * Typical Usage
 *
 *     int ierr = ezcl_devtype_init(CL_DEVICE_TYPE_GPU);
 *******************************************************************/
#define ezcl_devtype_init(  device_type) \
      ( ezcl_devtype_init_p(device_type, __FILE__, __LINE__) )

#define ezcl_finalize() \
      ( ezcl_finalize_p(__FILE__,__LINE__) )
#define ezcl_terminate() \
      ( ezcl_terminate_p(__FILE__,__LINE__) )

/* device based routines */
#define ezcl_get_device(  context) \
      ( ezcl_get_device_p(context, __FILE__, __LINE__) ) 
#define ezcl_get_compute_device() \
      ( ezcl_get_compute_device_p(__FILE__, __LINE__) ) 
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
#define ezcl_malloc(  host_mem_ptr, name, dims, elsize, flags, ezcl_flags) \
      ( ezcl_malloc_p(host_mem_ptr, name, dims, elsize, flags, ezcl_flags, __FILE__, __LINE__) )
#define ezcl_malloc_memory_add(  malloc_mem_ptr, name, size) \
      ( ezcl_malloc_memory_add_p(malloc_mem_ptr, name, size, __FILE__, __LINE__) )
#define ezcl_device_memory_add(  dev_mem_ptr, name, num_elements, elsize, flags, ezcl_flags) \
      ( ezcl_device_memory_add_p(dev_mem_ptr, name, num_elements, elsize, flags, ezcl_flags, __FILE__, __LINE__) )
#define ezcl_mapped_memory_add(  map_mem_ptr, name, num_elements, elsize) \
      ( ezcl_mapped_memory_add_p(map_mem_ptr, name, num_elements, elsize, __FILE__, __LINE__) )
#define ezcl_device_memory_malloc(  context, host_mem_ptr, name, num_elements, elsize, flags, ezcl_flags) \
      ( ezcl_device_memory_malloc_p(context, host_mem_ptr, name, num_elements, elsize, flags, ezcl_flags, __FILE__, __LINE__) )
#define ezcl_device_memory_realloc(  dev_mem_ptr, num_elements) \
      ( ezcl_device_memory_realloc_p(dev_mem_ptr, num_elements, __FILE__, __LINE__) )
#define ezcl_device_memory_request(  dev_mem_ptr, capacity) \
      ( ezcl_device_memory_request_p(dev_mem_ptr, capacity, __FILE__, __LINE__) )
#define ezcl_device_memory_swap(  dev_mem_ptr_old, dev_mem_ptr_new) \
      ( ezcl_device_memory_swap_p(dev_mem_ptr_old, dev_mem_ptr_new, __FILE__, __LINE__) )
#define ezcl_device_memory_replace(  dev_mem_ptr_old, dev_mem_ptr_new) \
      ( ezcl_device_memory_replace_p(dev_mem_ptr_old, dev_mem_ptr_new, __FILE__, __LINE__) )
#define ezcl_device_memory_remove(  dev_mem_ptr) \
      ( ezcl_device_memory_remove_p(dev_mem_ptr, __FILE__, __LINE__) )
#define ezcl_mapped_memory_remove(  map_mem_ptr) \
      ( ezcl_mapped_memory_remove_p(map_mem_ptr, __FILE__, __LINE__) )
#define ezcl_malloc_memory_remove(  malloc_mem_ptr) \
      ( ezcl_malloc_memory_remove_p(malloc_mem_ptr, __FILE__, __LINE__) )
#define ezcl_device_memory_delete(  dev_mem_ptr) \
      ( ezcl_device_memory_delete_p(dev_mem_ptr, __FILE__, __LINE__) )
#define ezcl_mapped_memory_delete(  map_mem_ptr) \
      ( ezcl_mapped_memory_delete_p(map_mem_ptr, __FILE__, __LINE__) )
#define ezcl_malloc_memory_delete(  malloc_mem_ptr) \
      ( ezcl_malloc_memory_delete_p(malloc_mem_ptr, __FILE__, __LINE__) )
#define ezcl_mem_free_all() \
      ( ezcl_mem_free_all_p(__FILE__, __LINE__) ) 
#define ezcl_mem_walk_all() \
      ( ezcl_mem_walk_all_p(__FILE__, __LINE__) ) 
#define ezcl_mem_walk_one(mem_buffer) \
      ( ezcl_mem_walk_one_p(mem_buffer, __FILE__, __LINE__) ) 
#define ezcl_get_device_mem_nelements(  mem_buffer) \
      ( ezcl_get_device_mem_nelements_p(mem_buffer, __FILE__, __LINE__) ) 
#define ezcl_get_device_mem_elsize(  mem_buffer) \
      ( ezcl_get_device_mem_elsize_p(mem_buffer, __FILE__, __LINE__) ) 
#define ezcl_get_device_mem_capacity(  mem_buffer) \
      ( ezcl_get_device_mem_capacity_p(mem_buffer, __FILE__, __LINE__) ) 

/* kernel and program routines */
#define ezcl_create_kernel(  context, filename, kernel_name) \
      ( ezcl_create_kernel_p(context, filename, kernel_name, __FILE__, __LINE__) )
#define ezcl_create_kernel_wsource(  context, source, kernel_name) \
      ( ezcl_create_kernel_wsource_p(context, source, kernel_name, __FILE__, __LINE__) )
#define ezcl_create_program_wsource(  context, defines, source) \
      ( ezcl_create_program_wsource_p(context, defines, source, __FILE__, __LINE__) )
#define ezcl_create_kernel_wprogram(  program, kernel_name) \
      ( ezcl_create_kernel_wprogram_p(program, kernel_name, __FILE__, __LINE__) )

#define ezcl_program_release(  program) \
      ( ezcl_program_release_p(program, __FILE__, __LINE__) )
#define ezcl_kernel_release(  kernel) \
      ( ezcl_kernel_release_p(kernel, __FILE__, __LINE__) )
#define ezcl_command_queue_release(  command_queue) \
      ( ezcl_command_queue_release_p(command_queue, __FILE__, __LINE__) )
#define ezcl_context_release(  context) \
      ( ezcl_context_release_p(context, __FILE__, __LINE__) )
#define ezcl_event_release(  event) \
      ( ezcl_event_release_p(event, __FILE__, __LINE__) )

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
#define ezcl_wait_for_events(  num_events, events) \
      ( ezcl_wait_for_events_p(num_events, events, __FILE__, __LINE__) )

#define ezcl_timer_calc(  start_read_event, end_read_event) \
      ( ezcl_timer_calc_p(start_read_event, end_read_event, __FILE__, __LINE__) )
   

/* init and finish routines -- use either ezcl_init or ezcl_devtype_init */
/*   the devtype_init is the simpler and recommended path */
cl_int ezcl_init_p(cl_context *ezcl_gpu_context, cl_context *ezcl_cpu_context, cl_context *ezcl_accelerator_context, const char *file, const int line);
cl_int ezcl_devtype_init_p(cl_device_type device_type, const char *file, const int line);
cl_int ezcl_finalize_p(const char *file, const int line);
void ezcl_terminate_p(const char *file, const int line);

/* Error reporting */
void ezcl_print_error(const int ierr, const char *routine, const char *cl_routine, const char *file, const int line);

/* device based routines */
cl_device_id ezcl_get_device_p(cl_context context, const char *file, const int line);
int ezcl_get_compute_device_p(const char *file, const int line);
void ezcl_device_info_p(cl_device_id device, const char *file, const int line);

/* context based routines*/
cl_command_queue ezcl_create_command_queue_p(cl_context context, const int mype, const char *file, const int line);
cl_command_queue ezcl_get_command_queue_p(const char *file, const int line);
cl_context ezcl_get_context_p(const char *file, const int line);

/* memory routines */
cl_mem ezcl_malloc_p(void *host_mem_ptr, const char *name, size_t dims[], size_t elsize, size_t flags, int ezcl_flags, const char *file, const int line);
void *ezcl_malloc_memory_add_p(void *malloc_mem_ptr, const char *name, size_t size, const char *file, const int line);
void ezcl_device_memory_add_p(cl_mem dev_mem_ptr, const char *name, size_t num_elements, size_t elsize, size_t flags, int ezcl_flags, const char *file, const int line);
void ezcl_mapped_memory_add_p(cl_mem map_mem_ptr, const char *name, size_t num_elements, size_t elsize, const char *file, const int line);
cl_mem ezcl_device_memory_malloc_p(cl_context context, void *host_mem_ptr, const char *name, size_t num_elements, size_t elsize, size_t flags, int ezcl_flags, const char *file, const int line);
cl_mem ezcl_device_memory_realloc_p(cl_mem dev_mem_ptr, size_t num_elements, const char *file, const int line);   
cl_mem ezcl_device_memory_request_p(cl_mem dev_mem_ptr, size_t capacity, const char *file, const int line);   
void ezcl_device_memory_swap_p(cl_mem *dev_mem_ptr_old, cl_mem *dev_mem_ptr_new, const char *file, const int line);   
void ezcl_device_memory_replace_p(void **dev_mem_ptr_old, void **dev_mem_ptr_new, const char *file, const int line);   
void ezcl_device_memory_remove_p(void *dev_mem_ptr, const char *file, const int line);   
void ezcl_mapped_memory_remove_p(void *map_mem_ptr, const char *file, const int line);   
void ezcl_malloc_memory_remove_p(void *malloc_mem_ptr, const char *file, const int line);   
void ezcl_device_memory_delete_p(void *dev_mem_ptr, const char *file, const int line);   
void ezcl_mapped_memory_delete_p(void *map_mem_ptr, const char *file, const int line);   
void ezcl_malloc_memory_delete_p(void *malloc_mem_ptr, const char *file, const int line);   
void ezcl_mem_free_all_p(const char *file, const int line);
void ezcl_mem_walk_all_p(const char *file, const int line);
void ezcl_mem_walk_one_p(cl_mem mem_buffer, const char *file, const int line);
int ezcl_get_device_mem_nelements_p(cl_mem mem_buffer, const char *file, const int line);
int ezcl_get_device_mem_elsize_p(cl_mem mem_buffer, const char *file, const int line);
size_t ezcl_get_device_mem_capacity_p(cl_mem mem_buffer, const char *file, const int line);

/* kernel and program routines */
cl_kernel ezcl_create_kernel_p(cl_context context, const char *filename, const char *kernel_name, const char *file, const int line);
cl_kernel ezcl_create_kernel_wsource_p(cl_context context, const char *source, const char *kernel_name, const char *file, const int line);
cl_program ezcl_create_program_wsource_p(cl_context context, const char *defines, const char *source, const char *file, const int line);
cl_kernel ezcl_create_kernel_wprogram_p(cl_program program, const char *kernel_name, const char *file, const int line);
//cl_command_queue ezcl_create_command_queue(void);

void ezcl_program_release_p(cl_program program, const char *file, const int line);
void ezcl_kernel_release_p(cl_kernel kernel, const char *file, const int line);
void ezcl_command_queue_release_p(cl_command_queue command_queue, const char *file, const int line);
void ezcl_context_release_p(cl_context context, const char *file, const int line);
void ezcl_event_release_p(cl_event event, const char *file, const int line);

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
void ezcl_wait_for_events_p(unsigned int num_events, cl_event *events, const char *file, const int line);

void ezcl_set_timing_flag(void);
long ezcl_timer_calc_p(cl_event *start_read_event, cl_event *end_read_event, const char *file, const int line);
   
#ifdef __cplusplus
}
#endif

#endif  /* _EZCL_H */
