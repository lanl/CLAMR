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
#ifndef DEBUG
#define DEBUG 0
#endif

#ifndef DEVICE_DETECT_DEBUG
#define DEVICE_DETECT_DEBUG 0
#endif

#define EZCL_MEM_FACTOR 1.6
#define SUPPRESS_WARNING 1

#if !defined(FULL_PRECISION) && !defined(MIXED_PRECISION) && !defined(MINIMUM_PRECISION)
#define FULL_PRECISION
#endif
#ifdef NO_CL_DOUBLE
#undef  FULL_PRECISION
#undef  MIXED_PRECISION
#define MINIMUM_PRECISION
#endif

#include <stdlib.h>
#include <stdio.h>
#include <string.h>
#include <sys/queue.h>
#include <sys/stat.h>
#include <execinfo.h>

#include "ezcl.h"
#ifdef HAVE_MPI
#include "mpi.h"
#endif

typedef unsigned int uint;

// Control flags
static struct {
   unsigned int timing     : 1;
} ezcl_flags={1};

static cl_uint          nPlatforms = 0;
static cl_uint          nDevices = 0;
static cl_platform_id  *platforms;
static cl_device_id    *devices;
static cl_command_queue command_queue;
static cl_context       context;
static cl_program       program;
static cl_platform_id   platform = NULL;
static int              compute_device = 0;

SLIST_HEAD(slist_device_memory_head, device_memory_entry) device_memory_head = SLIST_HEAD_INITIALIZER(device_memory_head);
struct slist_device_memory_head *device_memory_headp;
struct device_memory_entry {
   cl_mem cl_mem_ptr;
   size_t num_elements;
   size_t el_size;
   size_t capacity;
   size_t flags;
   int    ezcl_flags;
   char   name[31];
   int    line;
   char   file[41];
   SLIST_ENTRY(device_memory_entry) device_memory_entries;
} *device_memory_item, *device_memory_item_old, *device_memory_item_new;

SLIST_HEAD(slist_mapped_memory_head, mapped_memory_entry) mapped_memory_head = SLIST_HEAD_INITIALIZER(mapped_memory_head);
struct slist_mapped_memory_head *mapped_memory_headp;
struct mapped_memory_entry {
   cl_mem cl_mem_ptr;
   size_t num_elements;
   size_t el_size;
   char   name[31];
   int    line;
   char   file[41];
   SLIST_ENTRY(mapped_memory_entry) mapped_memory_entries;
} *mapped_memory_item;

SLIST_HEAD(slist_malloc_memory_head, malloc_memory_entry) malloc_memory_head = SLIST_HEAD_INITIALIZER(malloc_memory_head);
struct slist_malloc_memory_head *malloc_memory_headp;
struct malloc_memory_entry {
   void *mem_ptr;
   size_t mem_size;
   char   name[31];
   int    line;
   char   file[41];
   SLIST_ENTRY(malloc_memory_entry) malloc_memory_entries;
} *malloc_memory_item;

SLIST_HEAD(slist_object_head, object_entry) object_head = SLIST_HEAD_INITIALIZER(object_head);
struct slist_object_head *object_headp;
struct object_entry {
   int object_type;
   cl_program program;
   cl_kernel kernel;
   cl_command_queue command_queue;
   cl_context context;
   cl_event event;
   char *filename;
   char name[31];
   int  line;
   char file[41];
   SLIST_ENTRY(object_entry) object_entries;
} *object_item;


cl_int ezcl_init_p(cl_context *ezcl_gpu_context, cl_context *ezcl_cpu_context, cl_context *ezcl_accelerator_context, const char *file, int line)
{
   int ierr;
   //cl_uint nPlatforms;
   //cl_uint nDevices = 0;

   ierr = clGetPlatformIDs(0, 0, &nPlatforms);
   if (ierr != CL_SUCCESS){
      printf("EZCL_INIT: Error with clGetPlatformIDs call in file %s at line %d\n", file, line);
      if (ierr == CL_INVALID_VALUE){
         printf("Invalid value in clGetPlatformID call\n");
      }
      exit(-1);
   }
   if (nPlatforms == 0) {
      printf("EZCL_INIT: Error -- No opencl platforms detected in file %s at line %d\n", file, line);
      exit(-1);
   }
   if (DEVICE_DETECT_DEBUG){
      printf("\n\nEZCL_INIT: %d opencl platform(s) detected\n",nPlatforms);
   }
   
   cl_uint num_platforms_requested = nPlatforms;
   platforms = (cl_platform_id *)malloc(nPlatforms*sizeof(cl_platform_id));
   ezcl_malloc_memory_add(platforms, "PLATFORMS", nPlatforms*sizeof(cl_platform_id));
   
   ierr = clGetPlatformIDs(num_platforms_requested, platforms, &nPlatforms);
   if (ierr != CL_SUCCESS || nPlatforms != num_platforms_requested){
      printf("EZCL_INIT: Error with clGetPlatformIDs call in file %s at line %d\n", file, line);
        if (ierr == CL_INVALID_VALUE){
          printf("Invalid value in clGetPlatformID call\n");
        }
   }

   if (DEVICE_DETECT_DEBUG){
      char info[1024];
      for (uint iplatform=0; iplatform<nPlatforms; iplatform++){
         printf("  Platform %d:\n",iplatform+1);
      
         //clGetPlatformInfo(platforms[iplatform],CL_PLATFORM_PROFILE,   1024L,info,0);
         //printf("    CL_PLATFORM_PROFILE    : %s\n",info);
      
         clGetPlatformInfo(platforms[iplatform],CL_PLATFORM_VERSION,   1024L,info,0);
         printf("    CL_PLATFORM_VERSION    : %s\n",info);
      
         clGetPlatformInfo(platforms[iplatform],CL_PLATFORM_NAME,      1024L,info,0);
         printf("    CL_PLATFORM_NAME       : %s\n",info);
      
         clGetPlatformInfo(platforms[iplatform],CL_PLATFORM_VENDOR,    1024L,info,0);
         printf("    CL_PLATFORM_VENDOR     : %s\n",info);
      
         //clGetPlatformInfo(platforms[iplatform],CL_PLATFORM_EXTENSIONS,1024L,info,0);
         //printf("    CL_PLATFORM_EXTENSIONS : %s\n",info);
      }
      printf("\n");
   }
   
   clGetDeviceIDs(platforms[0],(size_t)CL_DEVICE_TYPE_ALL,0,0,&nDevices);
   if (nDevices == 0) {
      printf("EZCL_INIT: Error -- No opencl devices detected in file %s at line %d\n", file, line);
      //exit(-1);
   }
   if (ierr != CL_SUCCESS) {
     /* Possible Errors
      *  CL_INVALID_PLATFORM:
      *  CL_INVALID_DEVICE_TYPE:
      *  CL_INVALID_VALUE:
      *  CL_DEVICE_NOT_FOUND:
      */
     ezcl_print_error(ierr, "EZCL_INIT", "clGetDeviceIDs", file, line);
   }
   if (DEVICE_DETECT_DEBUG){
      printf("EZCL_INIT: %d opencl devices(s) detected at line\n",nDevices);
   }
   
   cl_uint num_devices_requested = nDevices;
   devices = (cl_device_id *)malloc(nDevices*sizeof(cl_device_id));
   ezcl_malloc_memory_add(devices, "DEVICES", nDevices*sizeof(cl_device_id));

   ierr = clGetDeviceIDs(platforms[0],(size_t)CL_DEVICE_TYPE_ALL,num_devices_requested,devices,&nDevices);
   if (ierr != CL_SUCCESS) {
     /* Possible Errors
      *  CL_INVALID_PLATFORM:
      *  CL_INVALID_DEVICE_TYPE:
      *  CL_INVALID_VALUE:
      *  CL_DEVICE_NOT_FOUND:
      */
     ezcl_print_error(ierr, "EZCL_INIT", "clGetDeviceIDs", file, line);
   }
   
   if (DEVICE_DETECT_DEBUG){
      for (uint idevice=0; idevice<nDevices; idevice++){
         printf(  "  Device %d:\n", idevice+1);
         ezcl_device_info(devices[idevice]);
      }
   }
   
   *ezcl_gpu_context = clCreateContextFromType(NULL, (size_t)CL_DEVICE_TYPE_GPU, NULL, NULL, &ierr);
   if (ierr == CL_INVALID_VALUE){
      printf("Invalid value in clCreateContext call\n");
      if (devices == NULL) printf("Devices is NULL\n");
   }
   if (ierr != CL_SUCCESS) {
     /* Possible Errors
      *  CL_INVALID_PLATFORM:
      *  CL_INVALID_VALUE:
      *  CL_INVALID_DEVICE_TYPE:
      *  CL_DEVICE_NOT_AVAILABLE:
      *  CL_DEVICE_NOT_FOUND:
      *  CL_OUT_OF_HOST_MEMORY:
      */
     ezcl_print_error(ierr, "EZCL_INIT", "clCreateContextFromType", file, line);
   }
   
   if (*ezcl_gpu_context != NULL){
      printf("EZCL_INIT: GPU device context created\n");
   } else {
      if (DEVICE_DETECT_DEBUG) printf("EZCL_INIT: No gpu device found in file %s at line %d\n", file, line);
   }
   
   *ezcl_cpu_context = clCreateContextFromType(0, (size_t)CL_DEVICE_TYPE_CPU, NULL, NULL, &ierr);
   if (ierr == CL_INVALID_VALUE){
      printf("Invalid value in clCreateContext call\n");
      if (devices == NULL) printf("Devices is NULL\n");
   }
   if (ierr != CL_SUCCESS) {
     /* Possible Errors
      *  CL_INVALID_PLATFORM:
      *  CL_INVALID_VALUE:
      *  CL_INVALID_DEVICE_TYPE:
      *  CL_DEVICE_NOT_AVAILABLE:
      *  CL_DEVICE_NOT_FOUND:
      *  CL_OUT_OF_HOST_MEMORY:
      */
     ezcl_print_error(ierr, "EZCL_INIT", "clCreateContextFromType", file, line);
   }
   if (*ezcl_cpu_context != NULL){
      printf("EZCL_INIT: CPU device context created\n");
   } else {
      if (DEBUG) printf("EZCL_INIT: No cpu device found\n");
   }

   *ezcl_accelerator_context = clCreateContextFromType(0, (size_t)CL_DEVICE_TYPE_ACCELERATOR, NULL, NULL, &ierr);
   if (ierr == CL_INVALID_VALUE){
      printf("Invalid value in clCreateContext call\n");
      if (devices == NULL) printf("Devices is NULL\n");
   }
   if (ierr != CL_SUCCESS) {
     /* Possible Errors
      *  CL_INVALID_PLATFORM:
      *  CL_INVALID_VALUE:
      *  CL_INVALID_DEVICE_TYPE:
      *  CL_DEVICE_NOT_AVAILABLE:
      *  CL_DEVICE_NOT_FOUND:
      *  CL_OUT_OF_HOST_MEMORY:
      */
     ezcl_print_error(ierr, "EZCL_INIT", "clCreateContextFromType", file, line);
   }
   if (*ezcl_accelerator_context != NULL){
      printf("EZCL_INIT: Accelerator device context created\n");
   } else {
      if (DEVICE_DETECT_DEBUG) printf("EZCL_INIT: No accelerator device found in file %s at line %d\n", file, line);
   }
    printf("\n");
   
   return(0);
}

cl_int ezcl_devtype_init_p(cl_device_type device_type, const char *file, int line){
   cl_device_id device;
   int ierr;
   cl_uint nDevices_selected=0;

   //cl_uint nPlatforms;
   //cl_uint nDevices;
   cl_int platform_selected = -1;

   // Get the number of platforms first, then allocate and get the platform
   ierr = clGetPlatformIDs(0, NULL, &nPlatforms);
   if (ierr != CL_SUCCESS){
      printf("EZCL_DEVTYPE_INIT: Error with clGetDeviceIDs call in file %s at line %d\n", file, line);
      if (ierr == CL_INVALID_VALUE){
         printf("Invalid value in clGetPlatformID call\n");
      }
      exit(-1);
   }
   if (nPlatforms == 0) {
      printf("EZCL_DEVTYPE_INIT: Error -- No opencl platforms detected in file %s at line %d\n", file, line);
      exit(-1);
   }
   if (DEVICE_DETECT_DEBUG){
      printf("\n\nEZCL_DEVTYPE_INIT: %d opencl platform(s) detected\n",nPlatforms);
   }

   platforms = (cl_platform_id *)malloc(nPlatforms*sizeof(cl_platform_id));
   ezcl_malloc_memory_add(platforms, "PLATFORMS", nPlatforms*sizeof(cl_platform_id));

   ierr = clGetPlatformIDs(nPlatforms, platforms, NULL);
   if (ierr != CL_SUCCESS){
      printf("EZCL_DEVTYPE_INIT: Error with clGetPlatformIDs call in file %s at line %d\n", file, line);
      if (ierr == CL_INVALID_VALUE){
         printf("Invalid value in clGetPlatformID call\n");
      }
   }

   if (DEVICE_DETECT_DEBUG){
      char info[1024];
      for (uint iplatform=0; iplatform<nPlatforms; iplatform++){
         printf("  Platform %d:\n",iplatform+1);

         //clGetPlatformInfo(platforms[iplatform],CL_PLATFORM_PROFILE,   1024L,info,0);
         //printf("    CL_PLATFORM_PROFILE    : %s\n",info);

         clGetPlatformInfo(platforms[iplatform],CL_PLATFORM_VERSION,   1024L,info,0);
         printf("    CL_PLATFORM_VERSION    : %s\n",info);

         clGetPlatformInfo(platforms[iplatform],CL_PLATFORM_NAME,      1024L,info,0);
         printf("    CL_PLATFORM_NAME       : %s\n",info);

         clGetPlatformInfo(platforms[iplatform],CL_PLATFORM_VENDOR,    1024L,info,0);
         printf("    CL_PLATFORM_VENDOR     : %s\n",info);

         //clGetPlatformInfo(platforms[iplatform],CL_PLATFORM_EXTENSIONS,1024L,info,0);
         //printf("    CL_PLATFORM_EXTENSIONS : %s\n",info);
      }
      printf("\n");
   }

   // Get the number of devices, allocate, and get the devices
   for (uint iplatform=0; iplatform<nPlatforms; iplatform++){
      ierr = clGetDeviceIDs(platforms[iplatform],device_type,0,NULL,&nDevices);
      if (ierr == CL_DEVICE_NOT_FOUND) {
         if (DEVICE_DETECT_DEBUG) {
           printf("Warning: Device of requested type not found for platform %d in clGetDeviceID call\n",iplatform);
         }
         continue;
      }
      if (ierr != CL_SUCCESS) {
        /* Possible Errors
         *  CL_INVALID_PLATFORM:
         *  CL_INVALID_DEVICE_TYPE:
         *  CL_INVALID_VALUE:
         *  CL_DEVICE_NOT_FOUND:
         */
        ezcl_print_error(ierr, "EZCL_DEVTYPE_INIT", "clGetDeviceIDs", file, line);
      }
      if (DEVICE_DETECT_DEBUG){
         printf("EZCL_DEVTYPE_INIT: %d opencl devices(s) detected\n",nDevices);
      }
      platform_selected = iplatform;
      platform = platforms[iplatform];
      nDevices_selected = nDevices;
   }

   if (platform_selected == -1){
      printf("Warning: Device of requested type not found in clGetDeviceID call\n");
      ezcl_malloc_memory_delete(platforms);
      return(EZCL_NODEVICE);
   }

   nDevices = nDevices_selected;

   devices = (cl_device_id *)malloc(nDevices*sizeof(cl_device_id));
   ezcl_malloc_memory_add(devices, "DEVICES", nDevices*sizeof(cl_device_id));

   ierr = clGetDeviceIDs(platforms[platform_selected],device_type,nDevices,devices,NULL);
   if (ierr != CL_SUCCESS) {
     /* Possible Errors
      *  CL_INVALID_PLATFORM:
      *  CL_INVALID_DEVICE_TYPE:
      *  CL_INVALID_VALUE:
      *  CL_DEVICE_NOT_FOUND:
      */
     ezcl_print_error(ierr, "EZCL_DEVTYPE_INIT", "clGetDeviceIDs", file, line);
   }

   if (DEVICE_DETECT_DEBUG){
      for (uint idevice=0; idevice<nDevices; idevice++){
         printf(  "  Device %d:\n", idevice+1);
         ezcl_device_info(devices[idevice]);
      }
   }

   cl_context_properties context_properties[3]=
   {
     CL_CONTEXT_PLATFORM,
     (cl_context_properties)platform,
     0
   };

   context = clCreateContext(context_properties, nDevices, devices, NULL, NULL, &ierr);
   if (ierr == CL_INVALID_VALUE){
      printf("Invalid value in clCreateContext call\n");
      if (devices == NULL) printf("Devices is NULL\n");
   }
   if (ierr != CL_SUCCESS) {
     /* Possible Errors
      *  CL_INVALID_PLATFORM:
      *  CL_INVALID_VALUE:
      *  CL_INVALID_DEVICE_TYPE:
      *  CL_DEVICE_NOT_AVAILABLE:
      *  CL_DEVICE_NOT_FOUND:
      *  CL_OUT_OF_HOST_MEMORY:
      */
     ezcl_print_error(ierr, "EZCL_DEVTYPE_INIT", "clCreateContext", file, line);
   }
   if (DEVICE_DETECT_DEBUG){
      if (context != NULL){
         if(device_type & CL_DEVICE_TYPE_CPU )
            printf("EZCL_DEVTYPE_INIT: CPU device context created\n");
         else if (device_type & CL_DEVICE_TYPE_GPU){
            printf("EZCL_DEVTYPE_INIT: GPU device context created\n");
         }
         else if (device_type & CL_DEVICE_TYPE_ACCELERATOR){
            printf("EZCL_DEVTYPE_INIT: ACCELERATOR device context created\n");
         }
         else if (device_type & CL_DEVICE_TYPE_DEFAULT){
            printf("EZCL_DEVTYPE_INIT: Default device context created\n");
         }
      } else {
         if (DEBUG == 2) printf("EZCL_DEVTYPE_INIT: No device of type specified found\n");
      }
   }

   object_item = malloc(sizeof(struct object_entry));
   object_item->object_type = CONTEXT_OBJECT;
   object_item->context = context;
   object_item->line = line;
   snprintf(object_item->file,(size_t)30,"%-30s",file);
   if (DEBUG) printf("EZCL_DEVTYPE_INIT: DEBUG -- context is %p\n",context);
   SLIST_INSERT_HEAD(&object_head, object_item, object_entries);

   clGetContextInfo(context, CL_CONTEXT_DEVICES, sizeof(device), &device, NULL);
   if (context == NULL){
      printf("EZCL_DEVTYPE_INIT: Failed to find device and setup context in file %s at line %d\n", file, line);
      exit(-1); /* No device is available, something is wrong */
   }
   if (DEVICE_DETECT_DEBUG == 2){
      ezcl_device_info(device);
   }

   char info[1024];

   clGetDeviceInfo(devices[0], CL_DEVICE_VENDOR, sizeof(info), &info, NULL);
   if (! strncmp(info,"NVIDIA",(size_t)6) ) compute_device = COMPUTE_DEVICE_NVIDIA;
   if (! strncmp(info,"Advanced Micro Devices",(size_t)6) ) compute_device = COMPUTE_DEVICE_ATI;
   if (! strncmp(info,"Intel(R) Corporation",(size_t)6) ) compute_device = COMPUTE_DEVICE_INTEL;
   //printf("DEBUG -- device vendor is |%s|, compute_device %d\n",info,compute_device);

   int mype = 0;
#ifdef HAVE_MPI
   int numpe = 1;
   int numpe_node = 1;
   int mpi_initialized = 0;
   MPI_Initialized(&mpi_initialized);
   if (mpi_initialized) {
      MPI_Comm_rank(MPI_COMM_WORLD, &mype);
      MPI_Comm_size(MPI_COMM_WORLD, &numpe);
      numpe_node = numpe;

#if OMPI_MAJOR_VERSION >= 1 && OMPI_MINOR_VERSION >= 7
      MPI_Comm   mpi_comm_node;
      MPI_Comm_split_type(MPI_COMM_WORLD, MPI_COMM_TYPE_SHARED, 0, MPI_INFO_NULL, &mpi_comm_node);
      MPI_Comm_size(mpi_comm_node, &numpe_node);
      //printf("EZCL_DEVTYPE_INIT: Info -- major_version %d minor_version %d\n",
      //       OMPI_MAJOR_VERSION,OMPI_MINOR_VERSION);
#endif
   }
#endif

#ifdef XXX
//#ifndef __APPLE_CC__
   if (numpe_node > (int)nDevices) {
      printf("%d:EZCL_DEVTYPE_INIT: Error -- not enough GPUs for mpi ranks. nDevices %d numpe_node %d\n",
             mype,nDevices,numpe_node);
#ifdef HAVE_MPI
      MPI_Abort(MPI_COMM_WORLD, -1);
#endif
      exit(-1); /* Not enough devices for mpi ranks */
   }
//#endif
#endif

   command_queue = ezcl_create_command_queue(context, mype);

   return(EZCL_SUCCESS);
}

cl_device_id ezcl_get_device_p(cl_context context, const char *file, const int line){
   cl_device_id device;
   int ierr;
   ierr = clGetContextInfo(context, CL_CONTEXT_DEVICES, sizeof(device), &device, NULL);
   if (ierr != CL_SUCCESS) {
     /* Possible Errors
      *  CL_INVALID_CONTEXT:
      *  CL_INVALID_VALUE:
      */
     ezcl_print_error(ierr, "EZCL_GET_DEVICE", "clGetContextInfo", file, line);
   }
   return(device);
}

int ezcl_get_compute_device_p(const char *file, const int line)
{
   if (! SUPPRESS_WARNING) {
      printf(" Error with ezcl_get_compute_device from file %s line %d\n",file,line);
   }
   return(compute_device);
}

void ezcl_device_info_p(cl_device_id device, const char *file, const int line){
   if (device == NULL) {
      printf(" Error with device in ezcl_device_info called from file %s line %d\n",file,line);
   }
   char info[1024];
   cl_bool iflag;
   cl_uint inum;
   size_t isize;
   cl_ulong ilong;
   cl_device_type device_type;
   cl_command_queue_properties iprop;

   clGetDeviceInfo(device,CL_DEVICE_TYPE,sizeof(device_type),&device_type,0);
   if( device_type & CL_DEVICE_TYPE_CPU )
      printf("    CL_DEVICE_TYPE                       : %s\n", "CL_DEVICE_TYPE_CPU");
   if( device_type & CL_DEVICE_TYPE_GPU )
      printf("    CL_DEVICE_TYPE                       : %s\n", "CL_DEVICE_TYPE_GPU");
   if( device_type & CL_DEVICE_TYPE_ACCELERATOR )
      printf("    CL_DEVICE_TYPE                       : %s\n", "CL_DEVICE_TYPE_ACCELERATOR");
   if( device_type & CL_DEVICE_TYPE_DEFAULT )
      printf("    CL_DEVICE_TYPE                       : %s\n", "CL_DEVICE_TYPE_DEFAULT");
               
   clGetDeviceInfo(device,CL_DEVICE_AVAILABLE,sizeof(iflag),&iflag,0);
   if (iflag == CL_TRUE) {
      printf(  "    CL_DEVICE_AVAILABLE                  : TRUE\n");
   } else {
      printf(  "    CL_DEVICE_AVAILABLE                  : FALSE\n");
   }
   
   clGetDeviceInfo(device, CL_DEVICE_VENDOR, sizeof(info), &info, NULL);
   printf(  "    CL_DEVICE_VENDOR                     : %s\n", info);
   
   clGetDeviceInfo(device, CL_DEVICE_NAME, sizeof(info), &info, NULL);
   printf(  "    CL_DEVICE_NAME                       : %s\n", info);
   
   clGetDeviceInfo(device, CL_DRIVER_VERSION, sizeof(info), &info, NULL);
   printf(  "    CL_DRIVER_VERSION                    : %s\n", info);
   
   clGetDeviceInfo(device, CL_DEVICE_VERSION, sizeof(info), &info, NULL);
   printf(  "    CL_DEVICE_VERSION                    : %s\n", info);
   
   clGetDeviceInfo(device,CL_DEVICE_MAX_COMPUTE_UNITS,sizeof(inum),&inum,0);
   printf(  "    CL_DEVICE_MAX_COMPUTE_UNITS          : %d\n", inum);
   
   clGetDeviceInfo(device,CL_DEVICE_MAX_WORK_ITEM_DIMENSIONS,sizeof(inum),&inum,0);
   printf(  "    CL_DEVICE_MAX_WORK_ITEM_DIMENSIONS   : %d\n", inum);
   
   size_t *item_sizes = (size_t *)malloc(inum*sizeof(size_t));
   clGetDeviceInfo(device,CL_DEVICE_MAX_WORK_ITEM_SIZES,sizeof(item_sizes),item_sizes,0);
   printf(  "    CL_DEVICE_MAX_WORK_ITEM_SIZES        : %ld %ld %ld\n",
         item_sizes[0], item_sizes[1], item_sizes[2]);
   free(item_sizes);
   
   clGetDeviceInfo(device,CL_DEVICE_MAX_WORK_GROUP_SIZE,sizeof(isize),&isize,0);
   printf(  "    CL_DEVICE_MAX_WORK_GROUP_SIZE        : %ld\n", isize);
   
   clGetDeviceInfo(device,CL_DEVICE_MAX_CLOCK_FREQUENCY,sizeof(inum),&inum,0);
   printf(  "    CL_DEVICE_MAX_CLOCK_FREQUENCY        : %d\n", inum);
   
   clGetDeviceInfo(device,CL_DEVICE_MAX_MEM_ALLOC_SIZE,sizeof(inum),&inum,0);
   printf(  "    CL_DEVICE_MAX_MEM_ALLOC_SIZE         : %d\n", inum);
   
#ifdef __APPLE_CC__
   clGetDeviceInfo(device,CL_DEVICE_GLOBAL_MEM_SIZE,sizeof(ilong),&ilong,0);
   printf(  "    CL_DEVICE_GLOBAL_MEM_SIZE            : %llu\n", ilong);
   
   clGetDeviceInfo(device,CL_DEVICE_GLOBAL_MEM_CACHE_SIZE,sizeof(ilong),&ilong,0);
   printf(  "    CL_DEVICE_GLOBAL_MEM_CACHE_SIZE      : %llu\n", ilong);
#else
   clGetDeviceInfo(device,CL_DEVICE_GLOBAL_MEM_SIZE,sizeof(ilong),&ilong,0);
   printf(  "    CL_DEVICE_GLOBAL_MEM_SIZE            : %lu\n", ilong);

   clGetDeviceInfo(device,CL_DEVICE_GLOBAL_MEM_CACHE_SIZE,sizeof(ilong),&ilong,0);
   printf(  "    CL_DEVICE_GLOBAL_MEM_CACHE_SIZE      : %lu\n", ilong);
#endif
   clGetDeviceInfo(device,CL_DEVICE_GLOBAL_MEM_CACHELINE_SIZE,sizeof(inum),&inum,0);
   printf(  "    CL_DEVICE_GLOBAL_MEM_CACHELINE_SIZE  : %d\n", inum);
   
   clGetDeviceInfo(device,CL_DEVICE_MAX_CONSTANT_ARGS,sizeof(inum),&inum,0);
   printf(  "    CL_DEVICE_GLOBAL_MAX_CONSTANT_ARGS   : %d\n", inum);
   
   clGetDeviceInfo(device,CL_DEVICE_ERROR_CORRECTION_SUPPORT,sizeof(iflag),&iflag,0);
   if (iflag == CL_TRUE) {
      printf(  "    CL_DEVICE_ERROR_CORRECTION_SUPPORT   : TRUE\n");
   } else {
      printf(  "    CL_DEVICE_ERROR_CORRECTION_SUPPORT   : FALSE\n");
   }
   
   clGetDeviceInfo(device,CL_DEVICE_PROFILING_TIMER_RESOLUTION,sizeof(isize),&isize,0);
   printf(  "    CL_DEVICE_PROFILING_TIMER_RESOLUTION : %ld nanosecs\n", isize);
   
   clGetDeviceInfo(device,CL_DEVICE_QUEUE_PROPERTIES,sizeof(iprop),&iprop,0);
   if (iprop & CL_QUEUE_PROFILING_ENABLE) {
      printf(  "    CL_DEVICE_QUEUE PROFILING            : AVAILABLE\n");
   } else {
      printf(  "    CL_DEVICE_QUEUE PROFILING            : NOT AVAILABLE\n");
   }
   
   clGetDeviceInfo(device, CL_DEVICE_EXTENSIONS, sizeof(info), &info, NULL);
   printf(  "    CL_DEVICE_EXTENSIONS                 : %s\n\n", info);
   
}

cl_command_queue ezcl_create_command_queue_p(cl_context context, const int mype, const char *file, const int line){
   cl_device_id *device;
   int ierr;
   cl_command_queue_properties queueProps;

   device = (cl_device_id *)malloc(nDevices*sizeof(cl_device_id));

   ierr = clGetContextInfo(context, CL_CONTEXT_DEVICES, sizeof(cl_device_id)*nDevices, device, NULL);
   if (ierr != CL_SUCCESS) {
     /* Possible Errors
      *  CL_INVALID_CONTEXT:
      *  CL_INVALID_VALUE:
      */
     ezcl_print_error(ierr, "EZCL_CREATE_COMMAND_QUEUE", "clGetContextInfo", file, line);
   }

   if (DEVICE_DETECT_DEBUG) printf("%d: DEBUG -- running on device %d\n",mype,mype%nDevices);

   queueProps = ezcl_flags.timing ? CL_QUEUE_PROFILING_ENABLE : 0;
   command_queue = clCreateCommandQueue(context, device[mype%nDevices], queueProps, &ierr);
   if (ierr != CL_SUCCESS) {
     /* Possible Errors
      *  CL_INVALID_CONTEXT:
      *  CL_INVALID_DEVICE:
      *  CL_INVALID_VALUE:
      *  CL_OUT_OF_HOST_MEMORY:
      */
     ezcl_print_error(ierr, "EZCL_CREATE_COMMAND_QUEUE", "clCreateCommandQueue", file, line);
   }

   free(device);

   object_item = malloc(sizeof(struct object_entry));
   object_item->object_type = COMMAND_QUEUE_OBJECT;
   object_item->command_queue = command_queue;
   object_item->line = line;
   snprintf(object_item->file,(size_t)30,"%-30s",file);
   if (DEBUG) printf("EZCL_CREATE_COMMAND_QUEUE: DEBUG -- command_queue is %p\n",command_queue);
   SLIST_INSERT_HEAD(&object_head, object_item, object_entries);

   return(command_queue);
}

cl_command_queue ezcl_get_command_queue_p(const char *file, const int line){
   if (command_queue == NULL) {
      printf(" Error with command_queue in ezcl_get_command_queue called from file %s line %d\n",file,line);
   }
   return(command_queue);
}

cl_context ezcl_get_context_p(const char *file, const int line){
   if (context == NULL) {
      printf(" Error with context in ezcl_get_context called from file %s line %d\n",file,line);
   }
   return(context);
}

cl_mem ezcl_malloc_p(void *host_mem_ptr, const char *name, size_t dims[], size_t elsize, size_t flags, int ezcl_flags, const char *file, const int line){
   void *buf_mem_ptr = NULL;
   cl_mem mem_ptr = NULL;
   size_t size=0;

   size=dims[0]*elsize;
   if (DEBUG) printf("EZCL_MALLOC: variable %20s dims %ld elsize is %ld called from file %s line %d\n",name,dims[0],elsize,file,line);

   if (flags & CL_MEM_ALLOC_HOST_PTR){
      if (ezcl_flags & EZCL_PINNED_MEMORY){
         buf_mem_ptr = ezcl_device_memory_malloc_p(context, host_mem_ptr, name, dims[0], elsize, flags, ezcl_flags, file, line);

         if (flags & CL_MEM_READ_ONLY){
           mem_ptr = clEnqueueMapBuffer(command_queue, buf_mem_ptr, CL_TRUE, (size_t)CL_MAP_WRITE, 0L, size, 0, NULL, NULL, NULL);
         }else if (flags & CL_MEM_WRITE_ONLY){
           mem_ptr = clEnqueueMapBuffer(command_queue, buf_mem_ptr, CL_TRUE, (size_t)CL_MAP_READ, 0L, size, 0, NULL, NULL, NULL);
         } else {
           mem_ptr = clEnqueueMapBuffer(command_queue, buf_mem_ptr, CL_TRUE, 0L, 0L, size, 0, NULL, NULL, NULL);
         }
         ezcl_mapped_memory_add_p(mem_ptr, name, dims[0], elsize, file, line);

      } else {
         mem_ptr = malloc(size);
         ezcl_malloc_memory_add(mem_ptr, name, size);
      }
   } else {
      //printf("DEBUG line 691 ezcl_flags %d\n",ezcl_flags);
      mem_ptr = ezcl_device_memory_malloc_p(context, host_mem_ptr, name, dims[0], elsize, flags, ezcl_flags, file, line);
   }

   return(mem_ptr);
}

void ezcl_device_memory_add_p(cl_mem dev_mem_ptr, const char *name, size_t num_elements, size_t elsize, size_t flags, int ezcl_flags, const char *file, const int line){
   if (SLIST_EMPTY(&device_memory_head)) SLIST_INIT(&device_memory_head);

   device_memory_item = malloc(sizeof(struct device_memory_entry));
   snprintf(device_memory_item->name,(size_t)30,"%s",name);
   device_memory_item->cl_mem_ptr = dev_mem_ptr;
   device_memory_item->num_elements=num_elements;
   device_memory_item->flags=flags;
   device_memory_item->ezcl_flags=ezcl_flags;
   device_memory_item->capacity=num_elements;
   device_memory_item->el_size=elsize;
   device_memory_item->line=line;
   snprintf(device_memory_item->file,(size_t)40,"%s",file);
   if (DEBUG) printf("EZCL_DEVICE_MEMORY_ADD: DEBUG -- cl memory pointer is %p\n",dev_mem_ptr);
   SLIST_INSERT_HEAD(&device_memory_head, device_memory_item, device_memory_entries);
}

cl_mem ezcl_device_memory_malloc_p(cl_context context, void *host_mem_ptr, const char *name, size_t num_elements, size_t elsize, size_t flags, int ezcl_flags, const char *file, const int line){
   if (SLIST_EMPTY(&device_memory_head)) SLIST_INIT(&device_memory_head);
   device_memory_item = malloc(sizeof(struct device_memory_entry));
   snprintf(device_memory_item->name,(size_t)30,"%s",name);

   int ierr=0;
   size_t size = 0;
   //printf("DEBUG ezcl device memory malloc %d test %d\n",ezcl_flags,ezcl_flags & EZCL_MANAGED_MEMORY);
   if (ezcl_flags & EZCL_MANAGED_MEMORY) {
      size = EZCL_MEM_FACTOR*num_elements*elsize;
      device_memory_item->capacity=EZCL_MEM_FACTOR*num_elements;
   } else {
      size = num_elements*elsize;
      device_memory_item->capacity=num_elements;
   }
   cl_mem dev_mem_ptr = clCreateBuffer(context, flags, size, host_mem_ptr, &ierr);
   if (ierr != CL_SUCCESS) {
       /* Possible Errors
        *  CL_INVALID_CONTEXT:
        *  CL_INVALID_VALUE:
        *  CL_INVALID_BUFFER_SIZE:
        *  CL_INVALID_HOST_PTR:
        *  CL_MEM_OBJECT_ALLOCATION_FAILURE:
        *  CL_OUT_OF_HOST_MEMORY:
        */
      if (compute_device == COMPUTE_DEVICE_NVIDIA) {
         system("nvidia-smi -q -d MEMORY");
      }
      ezcl_print_error(ierr, "EZCL_DEVICE_MEMORY_MALLOC", "clCreateBuffer", file, line);
   }

   device_memory_item->cl_mem_ptr = dev_mem_ptr;
   device_memory_item->num_elements=num_elements;
   device_memory_item->flags=flags;
   device_memory_item->ezcl_flags=ezcl_flags;
   device_memory_item->el_size=elsize;
   device_memory_item->line=line;
   snprintf(device_memory_item->file,(size_t)40,"%s",file);
   if (DEBUG) printf("EZCL_DEVICE_MEMORY_MALLOC: DEBUG -- cl memory pointer is %p\n",dev_mem_ptr);
   SLIST_INSERT_HEAD(&device_memory_head, device_memory_item, device_memory_entries);

   return(dev_mem_ptr);
}

cl_mem ezcl_device_memory_realloc_p(cl_mem dev_mem_ptr, size_t num_elements, const char *file, const int line){

   SLIST_FOREACH(device_memory_item, &device_memory_head, device_memory_entries){
      if (device_memory_item->cl_mem_ptr == dev_mem_ptr) {

         if (device_memory_item->ezcl_flags & EZCL_MANAGED_MEMORY){
            if (num_elements > device_memory_item->capacity){

               int ierr=0;
               size_t size = EZCL_MEM_FACTOR*num_elements*device_memory_item->el_size;
               device_memory_item->capacity=EZCL_MEM_FACTOR*num_elements;

               cl_mem dev_mem_ptr = clCreateBuffer(context, device_memory_item->flags, size, NULL, &ierr);
               if (ierr != CL_SUCCESS) {
                   /* Possible Errors
                    *  CL_INVALID_CONTEXT:
                    *  CL_INVALID_VALUE:
                    *  CL_INVALID_BUFFER_SIZE:
                    *  CL_INVALID_HOST_PTR:
                    *  CL_MEM_OBJECT_ALLOCATION_FAILURE:
                    *  CL_OUT_OF_HOST_MEMORY:
                    */
                  ezcl_print_error(ierr, "EZCL_DEVICE_MEMORY_REALLOC", "clCreateBuffer", file, line);
               }

               device_memory_item->cl_mem_ptr = dev_mem_ptr;
               device_memory_item->num_elements=num_elements;
               device_memory_item->line=line;
               snprintf(device_memory_item->file,(size_t)40,"%s",file);
            } else {
               device_memory_item->num_elements = num_elements;
            }
         } else {
            int ierr=0;
            size_t size = num_elements*device_memory_item->el_size;

            cl_mem dev_mem_ptr = clCreateBuffer(context, device_memory_item->flags, size, NULL, &ierr);
            if (ierr != CL_SUCCESS) {
                /* Possible Errors
                 *  CL_INVALID_CONTEXT:
                 *  CL_INVALID_VALUE:
                 *  CL_INVALID_BUFFER_SIZE:
                 *  CL_INVALID_HOST_PTR:
                 *  CL_MEM_OBJECT_ALLOCATION_FAILURE:
                 *  CL_OUT_OF_HOST_MEMORY:
                 */
               ezcl_print_error(ierr, "EZCL_DEVICE_MEMORY_REALLOC", "clCreateBuffer", file, line);
            }
            device_memory_item->cl_mem_ptr = dev_mem_ptr;
            device_memory_item->num_elements=num_elements;
            device_memory_item->capacity=num_elements;
            device_memory_item->line=line;
            snprintf(device_memory_item->file,(size_t)40,"%s",file);
         }
         if (DEBUG) printf("EZCL_DEVICE_MEMORY_REALLOC: DEBUG -- cl memory pointer is %p\n",dev_mem_ptr);
      }
   }
   return(dev_mem_ptr);
}

cl_mem ezcl_device_memory_request_p(cl_mem dev_mem_ptr, size_t capacity, const char *file, const int line){

   SLIST_FOREACH(device_memory_item, &device_memory_head, device_memory_entries){
      if (device_memory_item->cl_mem_ptr == dev_mem_ptr) {

         int ierr=0;
         size_t size = capacity*device_memory_item->el_size;
         device_memory_item->capacity=capacity;

         cl_mem dev_mem_ptr = clCreateBuffer(context, device_memory_item->flags, size, NULL, &ierr);
         if (ierr != CL_SUCCESS) {
             /* Possible Errors
              *  CL_INVALID_CONTEXT:
              *  CL_INVALID_VALUE:
              *  CL_INVALID_BUFFER_SIZE:
              *  CL_INVALID_HOST_PTR:
              *  CL_MEM_OBJECT_ALLOCATION_FAILURE:
              *  CL_OUT_OF_HOST_MEMORY:
              */
            ezcl_print_error(ierr, "EZCL_DEVICE_MEMORY_REQUEST", "clCreateBuffer", file, line);
         }

         device_memory_item->cl_mem_ptr = dev_mem_ptr;
         device_memory_item->line=line;
         snprintf(device_memory_item->file,(size_t)40,"%s",file);

         if (DEBUG) printf("EZCL_DEVICE_MEMORY_REQUEST: DEBUG -- cl memory pointer is %p\n",dev_mem_ptr);
      }
   }
   return(dev_mem_ptr);
}

void ezcl_mapped_memory_add_p(cl_mem map_mem_ptr, const char *name, size_t num_elements, size_t elsize, const char *file, const int line){
   if (SLIST_EMPTY(&mapped_memory_head)) SLIST_INIT(&mapped_memory_head);

   mapped_memory_item = malloc(sizeof(struct mapped_memory_entry));
   snprintf(mapped_memory_item->name,(size_t)30,"%s",name);
   mapped_memory_item->cl_mem_ptr = map_mem_ptr;
   mapped_memory_item->num_elements=num_elements;
   mapped_memory_item->el_size=elsize;
   mapped_memory_item->line=line;
   snprintf(mapped_memory_item->file,(size_t)40,"%s",file);
   if (DEBUG) printf("EZCL_MAPPED_MEMORY_ADD: DEBUG -- cl memory pointer is %p\n",map_mem_ptr);
   SLIST_INSERT_HEAD(&mapped_memory_head, mapped_memory_item, mapped_memory_entries);
}

void *ezcl_malloc_memory_add_p(void *malloc_mem_ptr, const char *name, size_t size, const char *file, const int line){
   if (SLIST_EMPTY(&malloc_memory_head)) SLIST_INIT(&malloc_memory_head);
   
   malloc_memory_item = malloc(sizeof(struct malloc_memory_entry));
   snprintf(malloc_memory_item->name,(size_t)30,"%s",name);
   malloc_memory_item->mem_ptr = malloc_mem_ptr;
   malloc_memory_item->mem_size = size;
   malloc_memory_item->line     = line;
   snprintf(malloc_memory_item->file,(size_t)40,"%s",file);
   if (DEBUG) printf("EZCL_MALLOC_MEMORY_ADD: DEBUG -- malloc memory pointer is %p\n",malloc_mem_ptr);

   SLIST_INSERT_HEAD(&malloc_memory_head, malloc_memory_item, malloc_memory_entries);

   return(malloc_mem_ptr);
}

void ezcl_device_memory_swap_p(cl_mem *dev_mem_ptr_old, cl_mem *dev_mem_ptr_new, const char *file, const int line){
   if (*dev_mem_ptr_old == NULL || *dev_mem_ptr_new == NULL) {
      printf(" Error with dev_mem_ptr in ezcl_device_memory_swap from file %s line %d\n",file,line);
   }
   SLIST_FOREACH(device_memory_item_old, &device_memory_head, device_memory_entries){
      if (device_memory_item_old->cl_mem_ptr == *dev_mem_ptr_old) {
         break;
      }
   }
   SLIST_FOREACH(device_memory_item_new, &device_memory_head, device_memory_entries){
      if (device_memory_item_new->cl_mem_ptr == *dev_mem_ptr_new) {
         break;
      }
   }
   char   name_tmp[31];
   strncpy(name_tmp,device_memory_item_old->name,31);
   strncpy(device_memory_item_old->name,device_memory_item_new->name,31);
   strncpy(device_memory_item_new->name,name_tmp,31);

   cl_mem dev_mem_ptr_tmp = *dev_mem_ptr_old;
          *dev_mem_ptr_old = *dev_mem_ptr_new;
          *dev_mem_ptr_new = dev_mem_ptr_tmp;
}

void ezcl_device_memory_replace_p(void **dev_mem_ptr_old, void **dev_mem_ptr_new, const char *file, const int line){
   if (*dev_mem_ptr_old == NULL || *dev_mem_ptr_new == NULL) {
      printf(" Error with dev_mem_ptr in ezcl_device_memory_swap from file %s line %d\n",file,line);
   }
   SLIST_FOREACH(device_memory_item_old, &device_memory_head, device_memory_entries){
      if (device_memory_item_old->cl_mem_ptr == *dev_mem_ptr_old) {
         break;
      }
   }
   SLIST_FOREACH(device_memory_item_new, &device_memory_head, device_memory_entries){
      if (device_memory_item_new->cl_mem_ptr == *dev_mem_ptr_new) {
         break;
      }
   }

   //printf("DEBUG -- oldname %s newname %s\n",device_memory_item_old->name,device_memory_item_new->name);
   // Move the old name to the new entry
   strncpy(device_memory_item_new->name,device_memory_item_old->name,31);

   // Delete the old entry
   ezcl_device_memory_delete(*dev_mem_ptr_old);
   *dev_mem_ptr_old = NULL;

   // Return the new pointer in the old location. Now only the new entry exists with the old name and as
   // the old pointer
   *dev_mem_ptr_old = *dev_mem_ptr_new;
}

void ezcl_device_memory_delete_p(void *dev_mem_ptr, const char *file, const int line){
   if (dev_mem_ptr == NULL) {
      printf(" Error with dev_mem_ptr in ezcl_device_memory_delete from file %s line %d\n",file,line);
   }
   SLIST_FOREACH(device_memory_item, &device_memory_head, device_memory_entries){
      if (device_memory_item->cl_mem_ptr == dev_mem_ptr) {
         if (DEBUG) printf("EZCL_DEVICE_MEMORY_DELETE: DEBUG -- freeing device memory pointer %p\n",dev_mem_ptr);
         clReleaseMemObject((void *)device_memory_item->cl_mem_ptr);
         SLIST_REMOVE(&device_memory_head, device_memory_item, device_memory_entry, device_memory_entries);
         free(device_memory_item);
         break;         
      }
   }
}   

void ezcl_mapped_memory_delete_p(void *map_mem_ptr, const char *file, const int line){
   if (map_mem_ptr == NULL) {
      printf(" Error with map_mem_ptr in ezcl_mapped_memory_delete from file %s line %d\n",file,line);
   }
   SLIST_FOREACH(mapped_memory_item, &mapped_memory_head, mapped_memory_entries){
      if (mapped_memory_item->cl_mem_ptr == map_mem_ptr) {
         if (DEBUG) printf("EZCL_MAPPED_MEMORY_DELETE: DEBUG -- freeing mapped memory pointer %p\n",map_mem_ptr);
         clReleaseMemObject((void *)mapped_memory_item->cl_mem_ptr);
         SLIST_REMOVE(&mapped_memory_head, mapped_memory_item, mapped_memory_entry, mapped_memory_entries);
         free(mapped_memory_item);
         break;
      }
   }
}

void ezcl_malloc_memory_delete_p(void *malloc_mem_ptr, const char *file, const int line){
   if (malloc_mem_ptr == NULL) {
      printf(" Error with malloc_mem_ptr in ezcl_malloc_memory_delete from file %s line %d\n",file,line);
   }
   SLIST_FOREACH(malloc_memory_item, &malloc_memory_head, malloc_memory_entries){
      if (malloc_memory_item->mem_ptr == malloc_mem_ptr) {
         if (DEBUG) printf("EZCL_MALLOC_MEMORY_DELETE: DEBUG -- freeing malloc memory pointer %p\n",malloc_mem_ptr);
         free(malloc_mem_ptr);
         SLIST_REMOVE(&malloc_memory_head, malloc_memory_item, malloc_memory_entry, malloc_memory_entries);
         free(malloc_memory_item);
         break;         
      }
   }
}   

void ezcl_device_memory_remove_p(void *dev_mem_ptr, const char *file, const int line){
   if (dev_mem_ptr == NULL) {
      printf(" Error with dev_mem_ptr in ezcl_device_memory_remove from file %s line %d\n",file,line);
   }
   SLIST_FOREACH(device_memory_item, &device_memory_head, device_memory_entries){
      if (device_memory_item->cl_mem_ptr == dev_mem_ptr) {
         if (DEBUG) printf("EZCL_DEVICE_MEMORY_REMOVE: DEBUG -- freeing device memory pointer %p\n",dev_mem_ptr);
         clReleaseMemObject((void *)device_memory_item->cl_mem_ptr);
         SLIST_REMOVE(&device_memory_head, device_memory_item, device_memory_entry, device_memory_entries);
         break;         
      }
   }
}   

void ezcl_mapped_memory_remove_p(void *map_mem_ptr, const char *file, const int line){
   if (map_mem_ptr == NULL) {
      printf(" Error with map_mem_ptr in ezcl_mapped_memory_remove from file %s line %d\n",file,line);
   }
   SLIST_FOREACH(mapped_memory_item, &mapped_memory_head, mapped_memory_entries){
      if (mapped_memory_item->cl_mem_ptr == map_mem_ptr) {
         if (DEBUG) printf("EZCL_MAPPED_MEMORY_REMOVE: DEBUG -- freeing mapped memory pointer %p\n",map_mem_ptr);
         clReleaseMemObject((void *)mapped_memory_item->cl_mem_ptr);
         SLIST_REMOVE(&mapped_memory_head, mapped_memory_item, mapped_memory_entry, mapped_memory_entries);
         break;
      }
   }
}

void ezcl_malloc_memory_remove_p(void *malloc_mem_ptr, const char *file, const int line){
   if (malloc_mem_ptr == NULL) {
      printf(" Error with malloc_mem_ptr in ezcl_malloc_memory_remove from file %s line %d\n",file,line);
   }
   SLIST_FOREACH(malloc_memory_item, &malloc_memory_head, malloc_memory_entries){
      if (malloc_memory_item->mem_ptr == malloc_mem_ptr) {
         if (DEBUG) printf("EZCL_MALLOC_MEMORY_REMOVE: DEBUG -- freeing malloc memory pointer %p\n",malloc_mem_ptr);
         free(malloc_mem_ptr);
         SLIST_REMOVE(&malloc_memory_head, malloc_memory_item, malloc_memory_entry, malloc_memory_entries);
         break;         
      }
   }
}   

void ezcl_program_release_p(cl_program program, const char *file, const int line){
   if (program == NULL) {
      printf(" Error with program in ezcl_program_release from file %s line %d\n",file,line);
   }
    SLIST_FOREACH(object_item, &object_head, object_entries){
        if (object_item->object_type == PROGRAM_OBJECT && object_item->program == program) {
            if (DEBUG) printf("EZCL_PROGRAM_RELEASE: DEBUG -- releasing program %p\n",program);
            clReleaseProgram(program);
            SLIST_REMOVE(&object_head, object_item, object_entry, object_entries);
            free(object_item);
            break;
        }
    }
}

void ezcl_kernel_release_p(cl_kernel kernel, const char *file, const int line){
   if (kernel == NULL) {
      printf(" Error with kernel in ezcl_kernel_release from file %s line %d\n",file,line);
   }
    SLIST_FOREACH(object_item, &object_head, object_entries){
        if (object_item->object_type == KERNEL_OBJECT && object_item->kernel == kernel) {
            if (DEBUG) printf("EZCL_KERNEL_RELEASE: DEBUG -- releasing kernel %p\n",kernel);
            clReleaseKernel(kernel);
            SLIST_REMOVE(&object_head, object_item, object_entry, object_entries);
            free(object_item);
            break;
        }
    }
}

void ezcl_command_queue_release_p(cl_command_queue command_queue, const char *file, const int line){
    if (command_queue == NULL) {
       printf(" Error with command_queue in ezcl_command_queue_release from file %s line %d\n",file,line);
    }
    SLIST_FOREACH(object_item, &object_head, object_entries){
        if (object_item->object_type == COMMAND_QUEUE_OBJECT && object_item->command_queue == command_queue) {
            if (DEBUG) printf("EZCL_COMMAND_QUEUE_RELEASE: DEBUG -- releasing command queue %p\n",command_queue);
            clReleaseCommandQueue(command_queue);
            SLIST_REMOVE(&object_head, object_item, object_entry, object_entries);
            free(object_item);
            break;
        }
    }
}


void ezcl_context_release_p(cl_context context, const char *file, const int line){
    if (context == NULL) {
       printf(" Error with context in ezcl_context_release from file %s line %d\n",file,line);
    }
    SLIST_FOREACH(object_item, &object_head, object_entries){
        if (object_item->object_type == CONTEXT_OBJECT && object_item->context == context) {
            if (DEBUG) printf("EZCL_CONTEXT_RELEASE: DEBUG -- releasing context %p\n",context);
            clReleaseContext(context);
            SLIST_REMOVE(&object_head, object_item, object_entry, object_entries);
            free(object_item);
            break;
        }
    }
}

void ezcl_event_release_p(cl_event event, const char *file, const int line){
    if (event == NULL) {
       printf(" Error with event in ezcl_event_release from file %s line %d\n",file,line);
    }
    SLIST_FOREACH(object_item, &object_head, object_entries){
        if (object_item->object_type == EVENT_OBJECT && object_item->event == event) {
            if (DEBUG) printf("EZCL_EVENT_RELEASE: DEBUG -- releasing event %p\n",event);
            clReleaseEvent(event);
            SLIST_REMOVE(&object_head, object_item, object_entry, object_entries);
            free(object_item);
            break;
        }
    }
}

void ezcl_mem_free_all_p(const char *file, const int line){
    if (! SUPPRESS_WARNING) {
       printf(" Error with ezcl_mem_free_all from file %s line %d\n",file,line);
    }
   while (!SLIST_EMPTY(&device_memory_head)) {
      device_memory_item = SLIST_FIRST(&device_memory_head);
      if (DEBUG) printf("EZCL_MEM_FREE_ALL: DEBUG -- freeing cl memory %p\n",device_memory_item->cl_mem_ptr);
      clReleaseMemObject((cl_mem)device_memory_item->cl_mem_ptr);
      SLIST_REMOVE_HEAD(&device_memory_head, device_memory_entries);
      free(device_memory_item);
   }

   while (!SLIST_EMPTY(&malloc_memory_head)) {
      malloc_memory_item = SLIST_FIRST(&malloc_memory_head);
      if (DEBUG) printf("EZCL_MEM_FREE_ALL: DEBUG -- freeing malloc memory %p\n",malloc_memory_item->mem_ptr);
      free(malloc_memory_item->mem_ptr);
      SLIST_REMOVE_HEAD(&malloc_memory_head, malloc_memory_entries);
      free(malloc_memory_item);
   }
   
   while (!SLIST_EMPTY(&object_head)) {
      object_item = SLIST_FIRST(&object_head);
      switch(object_item->object_type){
      case PROGRAM_OBJECT:
        if (DEBUG) printf("EZCL_MEM_FREE_ALL: DEBUG -- releasing program %p %s\n",object_item->program, object_item->filename);
        free(object_item->filename);
        clReleaseProgram(object_item->program);
        break;
      case KERNEL_OBJECT:
        if (DEBUG) printf("EZCL_MEM_FREE_ALL: DEBUG -- releasing kernel %p\n",object_item->kernel);
        clReleaseKernel(object_item->kernel);
        break;
      case COMMAND_QUEUE_OBJECT:
        if (DEBUG) printf("EZCL_MEM_FREE_ALL: DEBUG -- releasing command_queue %p\n",object_item->command_queue);
        clReleaseCommandQueue(object_item->command_queue);
        break;
      case CONTEXT_OBJECT:
        if (DEBUG) printf("EZCL_MEM_FREE_ALL: DEBUG -- releasing context %p\n",object_item->context);
        clReleaseContext(object_item->context);
        break;
      case EVENT_OBJECT:
        if (DEBUG) printf("EZCL_MEM_FREE_ALL: DEBUG -- releasing event %p\n",object_item->event);
        clReleaseEvent(object_item->event);
        break;
      }
      SLIST_REMOVE_HEAD(&object_head, object_entries);
      free(object_item);
   }
   
}

int ezcl_get_device_mem_nelements_p(cl_mem mem_buffer, const char *file, const int line)
{
   if (mem_buffer == NULL) {
      printf(" Error with mem_buffer in ezcl_get_device_mem_nelements from file %s line %d\n",file,line);
   }
   int nelements = -1;
   if (! SLIST_EMPTY(&device_memory_head)){
      SLIST_FOREACH(device_memory_item, &device_memory_head, device_memory_entries){
         if (mem_buffer == device_memory_item->cl_mem_ptr){
            nelements = device_memory_item->num_elements;
            break;
         }
      }
   }
   return(nelements);
}

int ezcl_get_device_mem_elsize_p(cl_mem mem_buffer, const char *file, const int line)
{
   if (mem_buffer == NULL) {
      printf(" Error with mem_buffer in ezcl_get_device_mem_elsize from file %s line %d\n",file,line);
   }
   int elsize = -1;
   if (! SLIST_EMPTY(&device_memory_head)){
      SLIST_FOREACH(device_memory_item, &device_memory_head, device_memory_entries){
         if (mem_buffer == device_memory_item->cl_mem_ptr){
            elsize = device_memory_item->el_size;
            break;
         }
      }
   }
   return(elsize);
}

size_t ezcl_get_device_mem_capacity_p(cl_mem mem_buffer, const char *file, const int line)
{
   if (mem_buffer == NULL) {
      printf(" Error with mem_buffer in ezcl_get_device_mem_capacity from file %s line %d\n",file,line);
   }
   size_t capacity = 0;
   if (! SLIST_EMPTY(&device_memory_head)){
      SLIST_FOREACH(device_memory_item, &device_memory_head, device_memory_entries){
         if (mem_buffer == device_memory_item->cl_mem_ptr){
            capacity = device_memory_item->capacity;
            break;
         }
      }
   }
   return(capacity);
}

void ezcl_mem_walk_one_p(cl_mem mem_buffer, const char *file, const int line){
   if (mem_buffer == NULL) {
      printf(" Error with mem_buffer in ezcl_mem_walk_one from file %s line %d\n",file,line);
   }
   if (! SLIST_EMPTY(&device_memory_head)){
      printf("\n ------ DEVICE MEMORY ALLOCATIONS -----\n");
      SLIST_FOREACH(device_memory_item, &device_memory_head, device_memory_entries){
         if (mem_buffer == device_memory_item->cl_mem_ptr){
            printf("EZCL_MEM_WALK_ONE: cl device memory %p name %-30s num %ld elsize %ld at line %d in file %s\n",
               device_memory_item->cl_mem_ptr,
               device_memory_item->name,
               device_memory_item->num_elements,
               device_memory_item->el_size,
               device_memory_item->line,
               device_memory_item->file);
            break;
         }
      }
   }
}


void ezcl_mem_walk_all_p(const char *file, const int line){
    if (! SUPPRESS_WARNING) {
       printf(" Error with ezcl_mem_walk_all from file %s line %d\n",file,line);
    }
   if (! SLIST_EMPTY(&device_memory_head)){
      printf("\n ------ DEVICE MEMORY ALLOCATIONS -----\n");
      SLIST_FOREACH(device_memory_item, &device_memory_head, device_memory_entries){
         printf("EZCL_MEM_WALK_ALL: cl device memory %p name %-30s num %ld elsize %ld at line %d in file %s\n",
            device_memory_item->cl_mem_ptr,
            device_memory_item->name,
            device_memory_item->num_elements,
            device_memory_item->el_size,
            device_memory_item->line,
            device_memory_item->file);
      }
   }

   if (! SLIST_EMPTY(&malloc_memory_head)){
      printf("\n ------ MALLOC MEMORY ALLOCATIONS -----\n");
      SLIST_FOREACH(malloc_memory_item, &malloc_memory_head, malloc_memory_entries){
         printf("EZCL_MEM_WALK_ALL: cl malloc memory %p name %-30s size %ld at line %d in file %s\n",
            malloc_memory_item->mem_ptr,
            malloc_memory_item->name,
            malloc_memory_item->mem_size,
            malloc_memory_item->line,
            malloc_memory_item->file);
      }
   }

   if (! SLIST_EMPTY(&mapped_memory_head)){
      printf("\n ------ MAPPED MEMORY ALLOCATIONS -----\n");
      SLIST_FOREACH(mapped_memory_item, &mapped_memory_head, mapped_memory_entries){
         printf("EZCL_MEM_WALK_ALL: cl mapped memory %p name %-30s num %ld elsize %ld at line %d in file %s\n",
            mapped_memory_item->cl_mem_ptr,
            mapped_memory_item->name,
            mapped_memory_item->num_elements,
            mapped_memory_item->el_size,
            mapped_memory_item->line,
            mapped_memory_item->file);
      }
   }

   if (! SLIST_EMPTY(&object_head)){
      int nprogram_objects = 0;
      int nkernel_objects = 0;
      int ncommand_queue_objects = 0;
      int ncontext_objects = 0;
      int nevent_objects = 0;
      SLIST_FOREACH(object_item, &object_head, object_entries){
         switch (object_item->object_type){
         case PROGRAM_OBJECT:
            nprogram_objects++;
            break;
         case KERNEL_OBJECT:
            nkernel_objects++;
            break;
         case COMMAND_QUEUE_OBJECT:
            ncommand_queue_objects++;
            break;
         case CONTEXT_OBJECT:
            ncontext_objects++;
            break;
         case EVENT_OBJECT:
            nevent_objects++;
            break;
         }
      }

      printf("\n ====== OBJECT ALLOCATIONS ======\n");

      if (nprogram_objects) {
         printf("\n ------ PROGRAM OBJECTS -----\n");
         SLIST_FOREACH(object_item, &object_head, object_entries){
            if (object_item->object_type == PROGRAM_OBJECT){
               printf("EZCL_MEM_WALK_ALL: program object %p name %-30s at line %d in file %s\n",
                  object_item->program,
                  object_item->name,
                  object_item->line,
                  object_item->file);
            }
         }
      }

      if (nkernel_objects) {
         printf("\n ------ KERNEL OBJECTS -----\n");
         SLIST_FOREACH(object_item, &object_head, object_entries){
            if (object_item->object_type == KERNEL_OBJECT){
               printf("EZCL_MEM_WALK_ALL: kernel object %p name %-30s at line %d in file %s\n",
                  object_item->kernel,
                  object_item->name,
                  object_item->line,
                  object_item->file);
            }
         }
      }

      if (ncommand_queue_objects) {
         printf("\n ------ COMMAND_QUEUE OBJECTS -----\n");
         SLIST_FOREACH(object_item, &object_head, object_entries){
            if (object_item->object_type == COMMAND_QUEUE_OBJECT){
               printf("EZCL_MEM_WALK_ALL: command_queue object %p at line %d in file %s\n",
                  object_item->command_queue,
                  object_item->line,
                  object_item->file);
            }
         }
      }

      if (ncontext_objects) {
         printf("\n ------ CONTEXT OBJECTS -----\n");
         SLIST_FOREACH(object_item, &object_head, object_entries){
            if (object_item->object_type == CONTEXT_OBJECT){
               printf("EZCL_MEM_WALK_ALL: context object %p at line %d in file %s\n",
                  object_item->context,
                  object_item->line,
                  object_item->file);
            }
         }
      }
   
      if (nevent_objects) {
         printf("\n ------ EVENT OBJECTS -----\n");
         SLIST_FOREACH(object_item, &object_head, object_entries){
            if (object_item->object_type == EVENT_OBJECT){
               printf("EZCL_MEM_WALK_ALL: event object %p name %-30s at line %d in file %s\n",
                  object_item->event,
                  object_item->name,
                  object_item->line,
                  object_item->file);
            }
         }
      }
   }
}

cl_kernel ezcl_create_kernel_p(cl_context context, const char *filename, const char *kernel_name, const char *file, const int line){
   cl_kernel kernel;
   struct stat statbuf;
   FILE *fh;
   char *source;
   cl_int ierr;
   size_t nReportSize;
   char *BuildReport;
   char *filename_copy;

   fh=fopen(filename, "r");
   if (fh == 0){
      printf("EZCL_CREATE_KERNEL: Error -- cannot open %s in file %s at line %d\n",filename, file, line);
      exit(-1);
   }
   
   stat(filename, &statbuf);
   source = (char *)malloc((size_t)(statbuf.st_size + 1));
   ezcl_malloc_memory_add(source, kernel_name, (size_t)(statbuf.st_size + 1));
   
   if (fread(source, (size_t)statbuf.st_size, 1L, fh) != 1) {
      printf("ERROR: problem reading program source file %s\n",filename);
      exit(-1);
   }
   source[statbuf.st_size] = '\0';

   program = clCreateProgramWithSource(context, 1, (const char **)&source, NULL, &ierr);
   if (ierr != CL_SUCCESS){
      printf("EZCL_CREATE_KERNEL: clCreateProgramWithSource returned an error %d in file %s at line %d\n",ierr, file, line);
      switch (ierr){
      case CL_INVALID_CONTEXT:
       printf("Invalid context in clCreateProgramWithSource\n");
       break;
      case CL_INVALID_VALUE:
        printf("Invalid value in clCreateProgramWithSource\n");
        break;
      case CL_OUT_OF_HOST_MEMORY:
        printf("Out of host memory in clCreateProgramWithSource\n");
        break;
     }
   }

   ezcl_malloc_memory_delete(source);

   char CompileString[80];
   CompileString[0] ='\0';

#ifdef HAVE_CL_DOUBLE
   strcat(CompileString,"-DHAVE_CL_DOUBLE");
#else
   strcat(CompileString,"-DNO_CL_DOUBLE -cl-single-precision-constant");
#endif
#ifdef FULL_PRECISION
   strcat(CompileString," -DFULL_PRECISION");
#endif
#ifdef MIXED_PRECISION
   strcat(CompileString," -DMIXED_PRECISION");
#endif
#ifdef MINIMUM_PRECISION
   strcat(CompileString," -DMINIMUM_PRECISION");
#endif

   if (compute_device == COMPUTE_DEVICE_NVIDIA) {
      strcat(CompileString," -DIS_NVIDIA");
   }

   ierr = clBuildProgram(program, 0, NULL, CompileString, NULL, NULL);

   if (ierr != CL_SUCCESS){
      printf("EZCL_CREATE_KERNEL: clBuildProgram returned an error %d in file %s at line %d\n",ierr, file, line);
        switch (ierr){
        case CL_INVALID_PROGRAM:
          printf("Invalid program in clBuildProgram\n");
          break;
        case CL_INVALID_VALUE:
          printf("Invalid value in clBuildProgram\n");
          break;
        case CL_INVALID_DEVICE:
          printf("Invalid device in clBuildProgram\n");
          break;
        case CL_INVALID_BUILD_OPTIONS:
          printf("Invalid build options in clBuildProgram\n");
          break;
        case CL_INVALID_OPERATION:
          printf("Invalid operation in clBuildProgram\n");
          break;
        case CL_COMPILER_NOT_AVAILABLE:
          printf("CL compiler not available in clBuildProgram\n");
          break;
        case CL_BUILD_PROGRAM_FAILURE:
          printf("Build program failure in clBuildProgram\n");
          ierr = clGetProgramBuildInfo(program, devices[0], CL_PROGRAM_BUILD_LOG, 0L, NULL, &nReportSize);
          if (ierr != CL_SUCCESS) {
            switch (ierr){
               case CL_INVALID_DEVICE:
                  printf("Invalid device in clProgramBuildInfo\n");
                  break;
               case CL_INVALID_VALUE:
                  printf("Invalid value in clProgramBuildInfo\n");
                  break;
               case CL_INVALID_PROGRAM:
                  printf("Invalid program in clProgramBuildInfo\n");
                  break;
               }
            }
                 
            BuildReport = (char *)malloc(nReportSize);
                 
            ierr = clGetProgramBuildInfo(program, devices[0], CL_PROGRAM_BUILD_LOG, nReportSize, BuildReport, NULL);
            if (ierr != CL_SUCCESS) {
               switch (ierr){
                  case CL_INVALID_DEVICE:
                     printf("Invalid device in clProgramBuildInfo\n");
                     break;
                  case CL_INVALID_VALUE:
                     printf("Invalid value in clProgramBuildInfo\n");
                     break;
                  case CL_INVALID_PROGRAM:
                     printf("Invalid program in clProgramBuildInfo\n");
                     break;
               }
            }
              
            printf("EZCL_CREATE_KERNEL: Build Log: %s\n",BuildReport);
            free(BuildReport);
            exit(-1);
            break;
        case CL_OUT_OF_HOST_MEMORY:
          printf("Out of host memory in clBuildProgram\n");
          break;
        }
   } else {
      if (DEBUG)
         printf("EZCL_CREATE_KERNEL: Build is SUCCESSFUL with no errors\n");
   }

   object_item = malloc(sizeof(struct object_entry));
   object_item->object_type = PROGRAM_OBJECT;
   object_item->program = program;
   snprintf(object_item->name,(size_t)30,"%-30s",kernel_name);
   object_item->line = line;
   snprintf(object_item->file,(size_t)30,"%-30s",file);
   filename_copy = (char *)malloc(strlen(filename) + 1);
   strcpy(filename_copy, filename);
   object_item->filename = filename_copy;
   if (DEBUG) printf("EZCL_CREATE_KERNEL: DEBUG -- program is %p\n",program);
   SLIST_INSERT_HEAD(&object_head, object_item, object_entries);

   kernel = clCreateKernel(program, kernel_name, &ierr);
   if (ierr != CL_SUCCESS){
      /* Possible Errors
        case CL_INVALID_PROGRAM:
        case CL_INVALID_PROGRAM_EXECUTABLE:
        case CL_INVALID_KERNEL_NAME:
        case CL_INVALID_KERNEL_DEFINITION:
        case CL_INVALID_VALUE:
        case CL_OUT_OF_HOST_MEMORY:
      */
      ezcl_print_error(ierr, "EZCL_CREATE_KERNEL", "clCreateKernel", file, line);
      exit(-1);
   }

   object_item = malloc(sizeof(struct object_entry));
   object_item->object_type = KERNEL_OBJECT;
   object_item->kernel = kernel;
   snprintf(object_item->name,(size_t)30,"%-30s",kernel_name);
   object_item->line = line;
   snprintf(object_item->file,(size_t)30,"%-30s",file);
   if (DEBUG) printf("EZCL_CREATE_KERNEL: DEBUG -- kernel is %p\n",kernel);
   SLIST_INSERT_HEAD(&object_head, object_item, object_entries);

   ezcl_program_release(program);

   return(kernel);
}

cl_kernel ezcl_create_kernel_wsource_p(cl_context context, const char *source, const char *kernel_name, const char *file, const int line){
   cl_kernel kernel;
   cl_int ierr;
   size_t nReportSize;
   char *BuildReport;
   char *filename_copy;
   
   program = clCreateProgramWithSource(context, 1, (const char **)&source, NULL, &ierr);
   if (ierr != CL_SUCCESS){
      printf("EZCL_CREATE_KERNEL_WSOURCE: clBuildProgramWithSource returned an error %d in file %s at line %d\n",ierr, file, line);
      switch (ierr){
      case CL_INVALID_CONTEXT:
       printf("Invalid context in clBuildProgramWithSource\n");
       break;
      case CL_INVALID_VALUE:
        printf("Invalid value in clBuildProgramWithSource\n");
        break;
      case CL_OUT_OF_HOST_MEMORY:
        printf("Out of host memory in clBuildProgramWithSource\n");
        break;
     }
   }

   char CompileString[80];
   CompileString[0] ='\0';

#ifdef HAVE_CL_DOUBLE
   strcat(CompileString,"-DHAVE_CL_DOUBLE");
#else
   strcat(CompileString,"-DNO_CL_DOUBLE -cl-single-precision-constant");
#endif
#ifdef FULL_PRECISION
   strcat(CompileString," -DFULL_PRECISION");
#endif
#ifdef MIXED_PRECISION
   strcat(CompileString," -DMIXED_PRECISION");
#endif
#ifdef MINIMUM_PRECISION
   strcat(CompileString," -DMINIMUM_PRECISION");
#endif

   if (compute_device == COMPUTE_DEVICE_NVIDIA) {
      strcat(CompileString," -DIS_NVIDIA");
   }

   ierr = clBuildProgram(program, 0, NULL, CompileString, NULL, NULL);

   if (ierr != CL_SUCCESS){
      printf("EZCL_CREATE_KERNEL_WSOURCE: clBuildProgram returned an error %d in file %s at line %d\n",ierr, file, line);
        switch (ierr){
        case CL_INVALID_PROGRAM:
          printf("Invalid program in clBuildProgram\n");
          break;
        case CL_INVALID_VALUE:
          printf("Invalid value in clBuildProgram\n");
          break;
        case CL_INVALID_DEVICE:
          printf("Invalid device in clBuildProgram\n");
          break;
        case CL_INVALID_BUILD_OPTIONS:
          printf("Invalid build options in clBuildProgram\n");
          break;
        case CL_INVALID_OPERATION:
          printf("Invalid operation in clBuildProgram\n");
          break;
        case CL_COMPILER_NOT_AVAILABLE:
          printf("CL compiler not available in clBuildProgram\n");
          break;
        case CL_BUILD_PROGRAM_FAILURE:
          printf("Build program failure in clBuildProgram\n");
          ierr = clGetProgramBuildInfo(program, devices[0], CL_PROGRAM_BUILD_LOG, 0L, NULL, &nReportSize);
          if (ierr != CL_SUCCESS) {
            switch (ierr){
               case CL_INVALID_DEVICE:
                  printf("Invalid device in clProgramBuildInfo\n");
                  break;
               case CL_INVALID_VALUE:
                  printf("Invalid value in clProgramBuildInfo\n");
                  break;
               case CL_INVALID_PROGRAM:
                  printf("Invalid program in clProgramBuildInfo\n");
                  break;
               }
            }
                 
            BuildReport = (char *)malloc(nReportSize);
                 
            ierr = clGetProgramBuildInfo(program, devices[0], CL_PROGRAM_BUILD_LOG, nReportSize, BuildReport, NULL);
            if (ierr != CL_SUCCESS) {
               switch (ierr){
                  case CL_INVALID_DEVICE:
                     printf("Invalid device in clProgramBuildInfo\n");
                     break;
                  case CL_INVALID_VALUE:
                     printf("Invalid value in clProgramBuildInfo\n");
                     break;
                  case CL_INVALID_PROGRAM:
                     printf("Invalid program in clProgramBuildInfo\n");
                     break;
               }
            }
              
            printf("EZCL_CREATE_KERNEL_WSOURCE: Build Log: %s\n",BuildReport);
            free(BuildReport);
            exit(-1);
            break;
        case CL_OUT_OF_HOST_MEMORY:
          printf("Out of host memory in clBuildProgram\n");
          break;
        }
   } else {
      if (DEBUG)
         printf("EZCL_CREATE_KERNEL_WSOURCE: Build is SUCCESSFUL with no errors\n");
   }

   object_item = malloc(sizeof(struct object_entry));
   object_item->object_type = PROGRAM_OBJECT;
   object_item->program = program;
   snprintf(object_item->name,(size_t)30,"%-30s",kernel_name);
   object_item->line = line;
   snprintf(object_item->file,(size_t)30,"%-30s",file);
   char *filename = "dummy";
   filename_copy = (char *)malloc(strlen(filename) + 1);
   strcpy(filename_copy, filename);
   object_item->filename = filename_copy;
   if (DEBUG) printf("EZCL_CREATE_KERNEL_WSOURCE: DEBUG -- program is %p\n",program);
   SLIST_INSERT_HEAD(&object_head, object_item, object_entries);

   kernel = clCreateKernel(program, kernel_name, &ierr);
   if (ierr != CL_SUCCESS){
      /* Possible Errors
        case CL_INVALID_PROGRAM:
        case CL_INVALID_PROGRAM_EXECUTABLE:
        case CL_INVALID_KERNEL_NAME:
        case CL_INVALID_KERNEL_DEFINITION:
        case CL_INVALID_VALUE:
        case CL_OUT_OF_HOST_MEMORY:
      */
      ezcl_print_error(ierr, "EZCL_CREATE_KERNEL_WSOURCE", "clCreateKernel", file, line);
      exit(-1);
   }

   object_item = malloc(sizeof(struct object_entry));
   object_item->object_type = KERNEL_OBJECT;
   object_item->kernel = kernel;
   snprintf(object_item->name,(size_t)30,"%-30s",kernel_name);
   object_item->line = line;
   snprintf(object_item->file,(size_t)30,"%-30s",file);
   if (DEBUG) printf("EZCL_CREATE_KERNEL_WSOURCE: DEBUG -- kernel is %p\n",kernel);
   SLIST_INSERT_HEAD(&object_head, object_item, object_entries);

   ezcl_program_release(program);

   return(kernel);
}

cl_program ezcl_create_program_wsource_p(cl_context context, const char *defines, const char *source, const char *file, const int line){
   cl_int ierr;
   size_t nReportSize;
   char *BuildReport;
   char *filename_copy;
   
   program = clCreateProgramWithSource(context, 1, (const char **)&source, NULL, &ierr);
   if (ierr != CL_SUCCESS){
      printf("EZCL_CREATE_PROGRAM_WSOURCE: clCreateProgramWithSource returned an error %d in file %s at line %d\n",ierr, file, line);
      switch (ierr){
      case CL_INVALID_CONTEXT:
       printf("Invalid context in clCreatProgramWithSource\n");
       break;
      case CL_INVALID_VALUE:
        printf("Invalid value in clCreateProgramWithSource\n");
        break;
      case CL_OUT_OF_HOST_MEMORY:
        printf("Out of host memory in clCreateProgramWithSource\n");
        break;
     }
   }

   int define_length = 63;
   if (defines != NULL) define_length = strlen(defines)+63;
   char *CompileString = (char *)malloc( define_length*sizeof(char) );
   CompileString[0]='\0';
   if (defines != NULL) strcat(CompileString, defines);

#ifdef HAVE_CL_DOUBLE
   strcat(CompileString, " -DHAVE_CL_DOUBLE");
#else
   strcat(CompileString, "-DNO_CL_DOUBLE -cl-single-precision-constant");
#endif
#ifdef FULL_PRECISION
   strcat(CompileString, " -DFULL_PRECISION");
#endif
#ifdef MIXED_PRECISION
   strcat(CompileString, " -DMIXED_PRECISION");
#endif
#ifdef MINIMUM_PRECISION
   strcat(CompileString, " -DMINIMUM_PRECISION");
#endif

   if (compute_device == COMPUTE_DEVICE_NVIDIA) {
      strcat(CompileString, " -DIS_NVIDIA");
   }

   //printf("DEBUG file %s line %d CompileString %s\n",__FILE__,__LINE__,CompileString);
   ierr = clBuildProgram(program, 0, NULL, CompileString, NULL, NULL);

   free(CompileString);

   if (ierr != CL_SUCCESS){
      printf("EZCL_CREATE_PROGRAM_WSOURCE: clBuildProgram returned an error %d in file %s at line %d\n",ierr, file, line);
        switch (ierr){
        case CL_INVALID_PROGRAM:
          printf("Invalid program in clBuildProgram\n");
          break;
        case CL_INVALID_VALUE:
          printf("Invalid value in clBuildProgram\n");
          break;
        case CL_INVALID_DEVICE:
          printf("Invalid device in clBuildProgram\n");
          break;
        case CL_INVALID_BUILD_OPTIONS:
          printf("Invalid build options in clBuildProgram\n");
          break;
        case CL_INVALID_OPERATION:
          printf("Invalid operation in clBuildProgram\n");
          break;
        case CL_COMPILER_NOT_AVAILABLE:
          printf("CL compiler not available in clBuildProgram\n");
          break;
        case CL_BUILD_PROGRAM_FAILURE:
          printf("Build program failure in clBuildProgram\n");
          ierr = clGetProgramBuildInfo(program, devices[0], CL_PROGRAM_BUILD_LOG, 0L, NULL, &nReportSize);
          if (ierr != CL_SUCCESS) {
            switch (ierr){
               case CL_INVALID_DEVICE:
                  printf("Invalid device in clProgramBuildInfo\n");
                  break;
               case CL_INVALID_VALUE:
                  printf("Invalid value in clProgramBuildInfo\n");
                  break;
               case CL_INVALID_PROGRAM:
                  printf("Invalid program in clProgramBuildInfo\n");
                  break;
               }
            }
                 
            BuildReport = (char *)malloc(nReportSize);
                 
            ierr = clGetProgramBuildInfo(program, devices[0], CL_PROGRAM_BUILD_LOG, nReportSize, BuildReport, NULL);
            if (ierr != CL_SUCCESS) {
               switch (ierr){
                  case CL_INVALID_DEVICE:
                     printf("Invalid device in clProgramBuildInfo\n");
                     break;
                  case CL_INVALID_VALUE:
                     printf("Invalid value in clProgramBuildInfo\n");
                     break;
                  case CL_INVALID_PROGRAM:
                     printf("Invalid program in clProgramBuildInfo\n");
                     break;
               }
            }
              
            printf("EZCL_CREATE_KERNEL: Build Log: %s\n",BuildReport);
            free(BuildReport);
            exit(-1);
            break;
        case CL_OUT_OF_HOST_MEMORY:
          printf("Out of host memory in clBuildProgram\n");
          break;
      }
   } else {
      if (DEBUG)
         printf("EZCL_CREATE_PROGRAM_WSOURCE: Build is SUCCESSFUL with no errors\n");

      if (compute_device == COMPUTE_DEVICE_INTEL) {
         ierr = clGetProgramBuildInfo(program, devices[0], CL_PROGRAM_BUILD_LOG, 0L, NULL, &nReportSize);
         BuildReport = (char *)malloc(nReportSize);
         ierr = clGetProgramBuildInfo(program, devices[0], CL_PROGRAM_BUILD_LOG, nReportSize, BuildReport, NULL);
         printf("EZCL_CREATE_PROGRAM_WSOURCE: Build Log: %s\n",BuildReport);
         free(BuildReport);
      }
   }

   object_item = malloc(sizeof(struct object_entry));
   object_item->object_type = PROGRAM_OBJECT;
   object_item->program = program;
   strcpy(object_item->name, "\0");
   object_item->line = line;
   snprintf(object_item->file,(size_t)30,"%-30s",file);
   char *filename = "dummy";
   filename_copy = (char *)malloc(strlen(filename) + 1);
   strcpy(filename_copy, filename);
   object_item->filename = filename_copy;
   if (DEBUG) printf("EZCL_CREATE_PROGRAM_WSOURCE: DEBUG -- program is %p\n",program);
   SLIST_INSERT_HEAD(&object_head, object_item, object_entries);

   return(program);
}

cl_kernel ezcl_create_kernel_wprogram_p(cl_program program, const char *kernel_name, const char *file, const int line){
   cl_kernel kernel;
   cl_int ierr;
   
   kernel = clCreateKernel(program, kernel_name, &ierr);
   if (ierr != CL_SUCCESS){
      /* Possible Errors
        case CL_INVALID_PROGRAM:
        case CL_INVALID_PROGRAM_EXECUTABLE:
        case CL_INVALID_KERNEL_NAME:
        case CL_INVALID_KERNEL_DEFINITION:
        case CL_INVALID_VALUE:
        case CL_OUT_OF_HOST_MEMORY:
      */
      ezcl_print_error(ierr, "EZCL_CREATE_PROGRAM", "clCreateKernel", file, line);
      exit(-1);
   }

   object_item = malloc(sizeof(struct object_entry));
   object_item->object_type = KERNEL_OBJECT;
   object_item->kernel = kernel;
   snprintf(object_item->name,(size_t)30,"%-30s",kernel_name);
   object_item->line = line;
   snprintf(object_item->file,(size_t)30,"%-30s",file);
   if (DEBUG) printf("EZCL_CREATE_KERNEL_WPROGRAM: DEBUG -- kernel is %p\n",kernel);
   SLIST_INSERT_HEAD(&object_head, object_item, object_entries);

   return(kernel);
}


void ezcl_enqueue_write_buffer_p(cl_command_queue command_queue, cl_mem mem_buffer, cl_bool blocking_flag, 
                            size_t offset, size_t size, const void *ptr, cl_event *event, const char *file, int line)
{
   int ierr;
   ierr=clEnqueueWriteBuffer(command_queue, mem_buffer, blocking_flag, offset, size, ptr, 0, NULL, event);
   if (ierr != CL_SUCCESS) {
     /* Possible Errors
      *  CL_INVALID_COMMAND_QUEUE:
      *  CL_INVALID_CONTEXT:
      *  CL_INVALID_MEM_OBJECT:
      *  CL_INVALID_VALUE:
      *  CL_INVALID_EVENT_WAIT_LIST:
      *  CL_MEM_OBJECT_ALLOCATION_FAILURE:
      *  CL_OUT_OF_HOST_MEMORY:
      */
     int num_elements = ezcl_get_device_mem_nelements_p(mem_buffer, file, line);
     int elsize = ezcl_get_device_mem_elsize_p(mem_buffer, file, line);
     //ezcl_mem_walk_one_p(mem_buffer, file, line);
     if (size + offset > (size_t)(num_elements*elsize)){
        printf("\nERROR: EZCL_ENQUEUE_WRITE_BUFFER -- buffer %p has %d elements with elsize %d for total of %d bytes but trying to write %ld bytes\n",mem_buffer, num_elements, elsize, num_elements*elsize, size);
     }
     ezcl_print_error(ierr, "EZCL_ENQUEUE_WRITE_BUFFER", "clEnqueueWriteBuffer", file, line);
   }
   if (event != NULL) {
      object_item = malloc(sizeof(struct object_entry));
      object_item->object_type = EVENT_OBJECT;
      object_item->event = *event;
      sprintf(object_item->name,"Write Buffer");
      object_item->line = line;
      snprintf(object_item->file,(size_t)30,"%-30s",file);
      if (DEBUG) printf("EZCL_ENQUEUE_WRITE_BUFFER: DEBUG -- event is %p\n",*event);
      SLIST_INSERT_HEAD(&object_head, object_item, object_entries);
   }
}

void ezcl_enqueue_read_buffer_p(cl_command_queue command_queue, cl_mem mem_buffer, cl_bool blocking_flag, 
                               size_t offset, size_t size, void *ptr, cl_event *event, const char *file, int line)
{
   int ierr;
   ierr=clEnqueueReadBuffer(command_queue, mem_buffer, blocking_flag, offset, size, ptr, 0, NULL, event);
   if (ierr != CL_SUCCESS) {
     /* Possible Errors
      *  CL_INVALID_COMMAND_QUEUE:
      *  CL_INVALID_CONTEXT:
      *  CL_INVALID_MEM_OBJECT:
      *  CL_INVALID_VALUE:
      *  CL_INVALID_EVENT_WAIT_LIST:
      *  CL_MEM_OBJECT_ALLOCATION_FAILURE:
      *  CL_OUT_OF_HOST_MEMORY:
      */
     int num_elements = ezcl_get_device_mem_nelements_p(mem_buffer, file, line);
     int elsize = ezcl_get_device_mem_elsize_p(mem_buffer, file, line);
     //ezcl_mem_walk_one_p(mem_buffer, file, line);
     if (size+offset > (size_t)(num_elements*elsize)){
        printf("\nERROR: EZCL_ENQUEUE_READ_BUFFER -- buffer %p has %d elements with element size %d for total of %d bytes but trying to read %ld bytes with an offset of %ld\n",mem_buffer, num_elements, elsize, num_elements*elsize, size, offset);
     }
     ezcl_print_error(ierr, "EZCL_ENQUEUE_READ_BUFFER", "clEnqueueReadBuffer", file, line);
   }
   if (event != NULL) {
      object_item = malloc(sizeof(struct object_entry));
      object_item->object_type = EVENT_OBJECT;
      object_item->event = *event;
      sprintf(object_item->name,"Read Buffer");
      object_item->line = line;
      snprintf(object_item->file,(size_t)30,"%-30s",file);
      if (DEBUG) printf("EZCL_ENQUEUE_READ_BUFFER: DEBUG -- event is %p\n",*event);
      SLIST_INSERT_HEAD(&object_head, object_item, object_entries);
   }
}

void ezcl_enqueue_ndrange_kernel_p(cl_command_queue command_queue, cl_kernel kernel, cl_uint work_dim, 
                                 const size_t *global_work_offset, const size_t *global_work_size, const size_t *local_work_size, cl_event *event,
                                 const char *file, const int line)
{
   int ierr; 
   ierr=clEnqueueNDRangeKernel(command_queue, kernel, work_dim, global_work_offset, global_work_size, local_work_size, 0, NULL, event);
   if (ierr != CL_SUCCESS) {
     /* Possible Errors
      *  CL_INVALID_PROGRAM_EXECUTABLE:
      *  CL_INVALID_COMMAND_QUEUE:
      *  CL_INVALID_KERNEL:
      *  CL_INVALID_CONTEXT:
      *  CL_INVALID_KERNEL_ARGS:
      *  CL_INVALID_WORK_DIMENSION:
      *  CL_INVALID_GLOBAL_WORK_SIZE:
      *  CL_INVALID_WORK_GROUP_SIZE:
      *  CL_INVALID_WORK_ITEM_SIZE:
      *  CL_INVALID_GLOBAL_OFFSET:
      *  CL_OUT_OF_RESOURCES:
      *  CL_MEM_OBJECT_ALLOCATION_FAILURE:
      *  CL_INVALID_EVENT_WAIT_LIST:
      *  CL_OUT_OF_HOST_MEMORY:
      */
     if (compute_device == COMPUTE_DEVICE_NVIDIA) {
         system("nvidia-smi -q -d MEMORY");
     }
     ezcl_print_error(ierr, "EZCL_ENQUEUE_NDRANGE_KERNEL", "clEnqueueNDRangeKernel", file, line);
   }
   if (event != NULL) {
      object_item = malloc(sizeof(struct object_entry));
      object_item->object_type = EVENT_OBJECT;
      object_item->event = *event;
      sprintf(object_item->name,"NDRange Kernel");
      object_item->line = line;
      snprintf(object_item->file,(size_t)30,"%-30s",file);
      if (DEBUG) printf("EZCL_ENQUEUE_NDRANGE_KERNEL: DEBUG -- event is %p\n",*event);
      SLIST_INSERT_HEAD(&object_head, object_item, object_entries);
   }
}

void ezcl_set_kernel_arg_p(cl_kernel kernel, cl_uint arg_index, size_t arg_size, const void *arg_value, const char *file, const int line)
{
   int ierr;
   ierr=clSetKernelArg(kernel, arg_index, arg_size, arg_value);
   if (ierr != CL_SUCCESS) {
     /* Possible Errors
      *  CL_INVALID_KERNEL:
      *  CL_INVALID_ARG_INDEX:
      *  CL_INVALID_ARG_VALUE:
      *  CL_INVALID_MEM_OBJECT:
      *  CL_INVALID_SAMPLER:
      *  CL_INVALID_ARG_SIZE:
      */
     ezcl_print_error(ierr, "EZCL_SET_KERNEL_ARG", "clSetKernelArg", file, line);
   }
}

void ezcl_flush_p(cl_command_queue command_queue, const char *file, const int line)
{
   int ierr;
   ierr=clFlush(command_queue);
   if (ierr != CL_SUCCESS) {
     /* Possible Errors
      *  CL_INVALID_COMMAND_QUEUE:
      *  CL_OUT_OF_HOST_MEMORY:
      */
     ezcl_print_error(ierr, "EZCL_FLUSH", "clFlush", file, line);
   }
}

void ezcl_finish_p(cl_command_queue command_queue, const char *file, const int line)
{
   int ierr;
   ierr=clFinish(command_queue);
   if (ierr != CL_SUCCESS) {
     /* Possible Errors
      *  CL_INVALID_COMMAND_QUEUE:
      *  CL_OUT_OF_HOST_MEMORY:
      */
     ezcl_print_error(ierr, "EZCL_FINISH", "clFinish", file, line);
   }
}

void ezcl_get_event_profiling_info_p(cl_event event, cl_profiling_info param_name, size_t param_value_size, void *param_value, const char *file, const int line)
{
   int ierr;
   ierr=clGetEventProfilingInfo(event, param_name, param_value_size, param_value, NULL);
   if (ierr != CL_SUCCESS) {
     /* Possible Errors
      *  CL_PROFILING_INFO_NOT_AVAILABLE:
      *  CL_INVALID_VALUE:
      *  CL_INVALID_EVENT:
      */
     ezcl_print_error(ierr, "EZCL_GET_EVENT_PROFILING_INFO", "clGetEventProfilingInfo", file, line);
   }
}

void ezcl_wait_for_events_p(unsigned int num_events, cl_event *events, const char *file, const int line)
{
   int ierr;
   ierr=clWaitForEvents(num_events, events);
   if (ierr != CL_SUCCESS) {
     /* Possible Errors
      *  CL_INVALID_VALUE:
      *  CL_INVALID_CONTEXT;
      *  CL_INVALID_EVENT:
      */
     ezcl_print_error(ierr, "EZCL_WAIT_FOR_EVENTS", "clWaitForEvents", file, line);
   }
}

void ezcl_set_timing_flag(void){
   ezcl_flags.timing = 1;
}

long ezcl_timer_calc_p(cl_event *start_read_event, cl_event *end_read_event, const char *file, const int line){
   long result;
   long gpu_time_start, gpu_time_end;

   int ierr = clWaitForEvents(1, end_read_event);
   if (ierr != CL_SUCCESS) {
     /* Possible Errors
      *  CL_INVALID_VALUE:
      *  CL_INVALID_CONTEXT:
      *  CL_INVALID_EVENT:
      */
      ezcl_print_error(ierr, "EZCL_TIMER_CALC", "clWaitForEvents", file, line);
   }

   ezcl_get_event_profiling_info_p(*start_read_event, CL_PROFILING_COMMAND_START, sizeof(gpu_time_start), &gpu_time_start, file, line);
   ezcl_get_event_profiling_info_p(*end_read_event,   CL_PROFILING_COMMAND_END,   sizeof(gpu_time_end),   &gpu_time_end,   file, line);
   result = gpu_time_end - gpu_time_start;
   if (*start_read_event != *end_read_event) {
      ezcl_event_release_p(*start_read_event, file, line);
   }
   ezcl_event_release_p(*end_read_event, file, line);

   return(result);
}


cl_int ezcl_finalize_p(const char *file, const int line){
   if (! SUPPRESS_WARNING) {
      printf(" Error with ezcl_finalize from file %s line %d\n",file,line);
   }

   ezcl_mem_free_all();
   
   return(0);
}

void ezcl_terminate_p(const char *file, const int line){
   if (! SUPPRESS_WARNING) {
      printf(" Error with ezcl_terminate from file %s line %d\n",file,line);
   }

   ezcl_malloc_memory_delete(devices);
   ezcl_malloc_memory_delete(platforms);

   ezcl_command_queue_release(command_queue);
   ezcl_context_release(context);
}

void ezcl_print_error(const int ierr, const char *routine, const char *cl_routine, const char *file, const int line)
{
   switch (ierr){
      case CL_DEVICE_NOT_FOUND:                //#define CL_DEVICE_NOT_FOUND                 -1
         printf("\nERROR: %s -- Device not found in %s at line %d in file %s\n", routine, cl_routine, line, file);
         break;
      case CL_DEVICE_NOT_AVAILABLE:            //#define CL_DEVICE_NOT_AVAILABLE             -2
         printf("\nERROR: %s -- Device not available in %s at line %d in file %s\n", routine, cl_routine, line, file);
         break;
      case CL_COMPILER_NOT_AVAILABLE:          //#define CL_COMPILER_NOT_AVAILABLE           -3
         printf("\nERROR: %s -- CL compiler not available failure in %s at line %d in file %s\n", routine, cl_routine, line, file);
         break;
      case CL_MEM_OBJECT_ALLOCATION_FAILURE:   //#define CL_MEM_OBJECT_ALLOCATION_FAILURE    -4
         printf("\nERROR: %s -- Mem object allocation failure in %s at line %d in file %s\n", routine, cl_routine, line, file);
         if (compute_device == COMPUTE_DEVICE_NVIDIA) {
            system("nvidia-smi -q -d MEMORY");
         }
         break;
      case CL_OUT_OF_RESOURCES:                //#define CL_OUT_OF_RESOURCES                 -5
         printf("\nERROR: %s -- Out of resources in %s at line %d in file %s\n", routine, cl_routine, line, file);
         break;
      case CL_OUT_OF_HOST_MEMORY:              //#define CL_OUT_OF_HOST_MEMORY               -6
         printf("\nERROR: %s -- Out of host memory in %s at line %d in file %s\n", routine, cl_routine, line, file);
         break;
      case CL_PROFILING_INFO_NOT_AVAILABLE:    //#define CL_PROFILING_INFO_NOT_AVAILABLE     -7
         printf("\nERROR: %s -- Profiling info not available in %s at line %d in file %s\n", routine, cl_routine, line, file);
         break;
      case CL_MEM_COPY_OVERLAP:                //#define CL_MEM_COPY_OVERLAP                 -8
         printf("\nERROR: %s -- Mem copy overlap in %s at line %d in file %s\n", routine, cl_routine, line, file);
         break;
      case CL_IMAGE_FORMAT_MISMATCH:           //#define CL_IMAGE_FORMAT_MISMATCH            -9
         printf("\nERROR: %s -- Image format mismatch in %s at line %d in file %s\n", routine, cl_routine, line, file);
         break;
      case CL_IMAGE_FORMAT_NOT_SUPPORTED:      //#define CL_IMAGE_FORMAT_NOT_SUPPORTED      -10
         printf("\nERROR: %s -- Image format not supported in %s at line %d in file %s\n", routine, cl_routine, line, file);
         break;
      case CL_BUILD_PROGRAM_FAILURE:           //#define CL_BUILD_PROGRAM_FAILURE           -11
         printf("\nERROR: %s -- Build program failure in %s at line %d in file %s\n", routine, cl_routine, line, file);
         break;
      case CL_MAP_FAILURE:                     //#define CL_MAP_FAILURE                     -12
         printf("\nERROR: %s -- Map failure in %s at line %d in file %s\n", routine, cl_routine, line, file);
         break;
#ifdef CL_MISALIGNED_SUB_BUFFER_OFFSET
      case CL_MISALIGNED_SUB_BUFFER_OFFSET:    //#define CL_MISALIGNED_SUB_BUFFER_OFFSET    -13
         printf("\nERROR: %s -- Misaligned sub buffer offset in %s at line %d in file %s\n", routine, cl_routine, line, file);
         break;
#endif
#ifdef CL_EXEC_STATUS_ERROR_FOR_EVENTS_IN_WAIT_LIST
      case CL_EXEC_STATUS_ERROR_FOR_EVENTS_IN_WAIT_LIST:  //#define CL_EXEC_STATUS_ERROR_FOR_EVENTS_IN_WAIT_LIST    -14
         printf("\nERROR: %s -- Error for events in wait list in %s at line %d in file %s\n", routine, cl_routine, line, file);
         break;
#endif
#ifdef CL_COMPILE_PROGRAM_FAILURE
      case CL_COMPILE_PROGRAM_FAILURE:         //#define CL_COMPILE_PROGRAM_FAILURE         -15
         printf("\nERROR: %s -- Compile program failure in %s at line %d in file %s\n", routine, cl_routine, line, file);
         break;
#endif
#ifdef CL_LINKER_NOT_AVAILABLE
      case CL_LINKER_NOT_AVAILABLE:            //#define CL_LINKER_NOT_AVAILABLE            -16
         printf("\nERROR: %s -- Linker not available in %s at line %d in file %s\n", routine, cl_routine, line, file);
         break;
#endif
#ifdef CL_LINK_PROGRAM_FAILURE
      case CL_LINK_PROGRAM_FAILURE:            //#define CL_LINK_PROGRAM_FAILURE            -17
         printf("\nERROR: %s -- Link program failure in %s at line %d in file %s\n", routine, cl_routine, line, file);
         break;
#endif
#ifdef CL_DEVICE_PARTITION_FAILED
      case CL_DEVICE_PARTITION_FAILED:         //#define CL_DEVICE_PARTITION_FAILED         -18
         printf("\nERROR: %s -- Device partition failed in %s at line %d in file %s\n", routine, cl_routine, line, file);
         break;
#endif
#ifdef CL_KERNEL_ARG_INFO_NOT_AVAILABLE
      case CL_KERNEL_ARG_INFO_NOT_AVAILABLE:   //#define CL_KERNEL_ARG_INFO_NOT_AVAILABLE   -19
         printf("\nERROR: %s -- Kernel arg info not available in %s at line %d in file %s\n", routine, cl_routine, line, file);
         break;
#endif

      case CL_INVALID_VALUE:                   //#define CL_INVALID_VALUE                   -30
         printf("\nERROR: %s -- Invalid value in %s at line %d in file %s\n", routine, cl_routine, line, file);
         break;
      case CL_INVALID_DEVICE_TYPE:             //#define CL_INVALID_DEVICE_TYPE             -31
         printf("\nERROR: %s -- Invalid device type in %s at line %d in file %s\n", routine, cl_routine, line, file);
         break;
      case CL_INVALID_PLATFORM:                //#define CL_INVALID_PLATFORM                -32
         printf("\nERROR: %s -- Invalid platform in %s at line %d in file %s\n", routine, cl_routine, line, file);
         break;
      case CL_INVALID_DEVICE:                  //#define CL_INVALID_DEVICE                  -33
         printf("\nERROR: %s -- Invalid device in %s at line %d in file %s\n", routine, cl_routine, line, file);
         break;
      case CL_INVALID_CONTEXT:                 //#define CL_INVALID_CONTEXT                 -34
         printf("\nERROR: %s -- Invalid context in %s at line %d in file %s\n", routine, cl_routine, line, file);
         break;
      case CL_INVALID_QUEUE_PROPERTIES:        //#define CL_INVALID_QUEUE_PROPERTIES        -35
         printf("\nERROR: %s -- Invalid queue properties in %s at line %d in file %s\n", routine, cl_routine, line, file);
         break;
      case CL_INVALID_COMMAND_QUEUE:           //#define CL_INVALID_COMMAND_QUEUE           -36
         printf("\nERROR: %s -- Invalid command queue in %s at line %d in file %s\n", routine, cl_routine, line, file);
         break;
      case CL_INVALID_HOST_PTR:                //#define CL_INVALID_HOST_PTR                -37
         printf("\nERROR: %s -- Invalid host pointer in %s at line %d in file %s\n", routine, cl_routine, line, file);
         break;
      case CL_INVALID_MEM_OBJECT:              //#define CL_INVALID_MEM_OBJECT              -38
         printf("\nERROR: %s -- Invalid memory object in %s at line %d in file %s\n", routine, cl_routine, line, file);
         break;
      case CL_INVALID_IMAGE_FORMAT_DESCRIPTOR: //#define CL_INVALID_IMAGE_FORMAT_DESCRIPTOR -39
         printf("\nERROR: %s -- Invalid image format descriptor in %s at line %d in file %s\n", routine, cl_routine, line, file);
         break;
      case CL_INVALID_IMAGE_SIZE:              //#define CL_INVALID_IMAGE_SIZE              -40
         printf("\nERROR: %s -- Invalid image size in %s at line %d in file %s\n", routine, cl_routine, line, file);
         break;
      case CL_INVALID_SAMPLER:                 //#define CL_INVALID_SAMPLER                 -41
         printf("\nERROR: %s -- Invalid sampler in %s at line %d in file %s\n", routine, cl_routine, line, file);
         break;
      case CL_INVALID_BINARY:                  //#define CL_INVALID_BINARY                  -42
         printf("\nERROR: %s -- Invalid binary in %s at line %d in file %s\n", routine, cl_routine, line, file);
         break;
      case CL_INVALID_BUILD_OPTIONS:           //#define CL_INVALID_BUILD_OPTIONS           -43
         printf("\nERROR: %s -- Invalid build options in %s at line %d in file %s\n", routine, cl_routine, line, file);
         break;
      case CL_INVALID_PROGRAM:                 //#define CL_INVALID_PROGRAM                 -44
         printf("\nERROR: %s -- Invalid program in %s at line %d in file %s\n", routine, cl_routine, line, file);
         break;
      case CL_INVALID_PROGRAM_EXECUTABLE:      //#define CL_INVALID_PROGRAM_EXECUTABLE      -45
         printf("\nERROR: %s -- Invalid program executable in %s at line %d in file %s\n", routine, cl_routine, line, file);
         break;
      case CL_INVALID_KERNEL_NAME:             //#define CL_INVALID_KERNEL_NAME             -46
         printf("\nERROR: %s -- Invalid kernel name in %s at line %d in file %s\n", routine, cl_routine, line, file);
         break;
      case CL_INVALID_KERNEL_DEFINITION:       //#define CL_INVALID_KERNEL_DEFINITION       -47
         printf("\nERROR: %s -- Invalid kernel definition in %s at line %d in file %s\n", routine, cl_routine, line, file);
         break;
      case CL_INVALID_KERNEL:                  //#define CL_INVALID_KERNEL                  -48
         printf("\nERROR: %s -- Invalid kernel in %s at line %d in file %s\n", routine, cl_routine, line, file);
         break;
      case CL_INVALID_ARG_INDEX:               //#define CL_INVALID_ARG_INDEX               -49
         printf("\nERROR: %s -- Invalid arg index in %s at line %d in file %s\n", routine, cl_routine, line, file);
         break;
      case CL_INVALID_ARG_VALUE:               //#define CL_INVALID_ARG_VALUE               -50
         printf("\nERROR: %s -- Invalid arg value in %s at line %d in file %s\n", routine, cl_routine, line, file);
         break;
      case CL_INVALID_ARG_SIZE:                //#define CL_INVALID_ARG_SIZE                -51
         printf("\nERROR: %s -- Invalid arg size in %s at line %d in file %s\n", routine, cl_routine, line, file);
         break;
      case CL_INVALID_KERNEL_ARGS:             //#define CL_INVALID_KERNEL_ARGS             -52
         printf("\nERROR: %s -- Invalid kernel args in %s at line %d in file %s\n", routine, cl_routine, line, file);
         break;
      case CL_INVALID_WORK_DIMENSION:          //#define CL_INVALID_WORK_DIMENSION          -53
         printf("\nERROR: %s -- Invalid work dimension in %s at line %d in file %s\n", routine, cl_routine, line, file);
         break;
      case CL_INVALID_WORK_GROUP_SIZE:         //#define CL_INVALID_WORK_GROUP_SIZE         -54
         printf("\nERROR: %s -- Invalid work group size in %s at line %d in file %s\n", routine, cl_routine, line, file);
         break;
      case CL_INVALID_WORK_ITEM_SIZE:          //#define CL_INVALID_WORK_ITEM_SIZE          -55
         printf("\nERROR: %s -- Invalid work item size in %s at line %d in file %s\n", routine, cl_routine, line, file);
         break;
      case CL_INVALID_GLOBAL_OFFSET:           //#define CL_INVALID_GLOBAL_OFFSET           -56
         printf("\nERROR: %s -- Invalid global offset in %s at line %d in file %s\n", routine, cl_routine, line, file);
         break;
      case CL_INVALID_EVENT_WAIT_LIST:         //#define CL_INVALID_EVENT_WAIT_LIST         -57
         printf("\nERROR: %s -- Invalid event wait list in %s at line %d in file %s\n", routine, cl_routine, line, file);
         break;
      case CL_INVALID_EVENT:                   //#define CL_INVALID_EVENT                   -58
         printf("\nERROR: %s -- Invalid event in %s at line %d in file %s\n", routine, cl_routine, line, file);
         break;
      case CL_INVALID_OPERATION:               //#define CL_INVALID_OPERATION               -59
         printf("\nERROR: %s -- Invalid operation in %s at line %d in file %s\n", routine, cl_routine, line, file);
         break;
      case CL_INVALID_GL_OBJECT:               //#define CL_INVALID_GL_OBJECT               -60
         printf("\nERROR: %s -- Invalid GL object in %s at line %d in file %s\n", routine, cl_routine, line, file);
         break;
      case CL_INVALID_BUFFER_SIZE:             //#define CL_INVALID_BUFFER_SIZE             -61
         printf("\nERROR: %s -- Invalid buffer size in %s at line %d in file %s\n", routine, cl_routine, line, file);
         break;
      case CL_INVALID_MIP_LEVEL:               //#define CL_INVALID_MIP_LEVEL               -62
         printf("\nERROR: %s -- Invalid mip level in %s at line %d in file %s\n", routine, cl_routine, line, file);
         break;
      case CL_INVALID_GLOBAL_WORK_SIZE:        //#define CL_INVALID_GLOBAL_WORK_SIZE        -63
         printf("\nERROR: %s -- Invalid global work size in %s at line %d in file %s\n", routine, cl_routine, line, file);
         break;

      default:
        printf("\nERROR: %s -- %d in %s at line %d in file %s\n", routine, ierr, cl_routine, line, file);
        break;
   }
   void* callstack[128];
   int frames = backtrace(callstack, 128);
   if (frames > 2) {
#ifdef HAVE_ADDR2LINE
      char hex_address[21];
      const char command_string[80];
#endif
      char** strs = backtrace_symbols(callstack, frames);
      fprintf(stderr,"\n  =============== Backtrace ===============\n");
      for (int i = 1; i < frames-1; ++i) {
          fprintf(stderr,"   %s    \t", strs[i]);
#ifdef HAVE_ADDR2LINE
          sscanf(strs[i],"%*s [%[^]]s]",hex_address);
          //printf("DEBUG addr2line -e clamr -f -s -i -p %s\n",hex_address);
          sprintf(command_string,"addr2line -e clamr -f -s -i -p %s",hex_address);
          system(command_string);
          // on mac, need to install binutils using macports "port install binutils"
#endif
          fprintf(stderr,"\n");
      }
      fprintf(stderr,"  =============== Backtrace ===============\n\n");
      free(strs);
   }

#ifdef HAVE_MPI
   MPI_Finalize();
#endif
   exit(-1);
}

