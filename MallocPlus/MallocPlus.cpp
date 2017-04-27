/*
 *  Copyright (c) 2011-2014, Los Alamos National Security, LLC.
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

// SKG TODO op realloc (similar to managed)

#include "MallocPlus.h"
#include <stdlib.h>
#include <unistd.h>
#include <stdio.h>
#include <algorithm>
#include <queue>
#include <string.h>
#ifdef HAVE_OPENCL
#include "ezcl/ezcl.h"
#endif

#ifndef DEBUG
#define DEBUG 0
#endif
#define WARNING_SUPPRESSION 0

#ifdef HAVE_CL_DOUBLE
#ifdef HAVE_OPENCL
typedef cl_double2  cl_real2;
#endif
#else
#ifdef HAVE_OPENCL
typedef cl_float2   cl_real2;
#endif
#endif

#ifndef MIN
#define MIN(a,b) ((a) < (b) ? (a) : (b))
#endif

#ifndef SWAP_PTR
#define SWAP_PTR(xnew,xold,xtmp) (xtmp=xnew, xnew=xold, xold=xtmp)
#endif

typedef unsigned int uint;
list<malloc_plus_memory_entry>::iterator it_save;

#if defined(HAVE_MPI)
void
MallocPlus::pinit(MPI_Comm smComm, std::size_t memPoolSize)
{
#if defined(HAVE_J7)
    try {
        j7 = new J7(smComm, memPoolSize);
    }
    catch(...) {
        std::cerr << "*** pinit failure ***" << std::endl;
        throw;
    }
#else
    // Just to suppress compiler warnings
    if (WARNING_SUPPRESSION) printf("DEBUG memPoolSize = %lu smComm = %p\n",memPoolSize,smComm);
#endif
}

void
MallocPlus::pfini(void)
{
#if defined(HAVE_J7)
    try {
        delete j7;
        j7 = NULL;
    }
    catch(...) {
        std::cerr << "*** pfini failure ***" << std::endl;
        throw;
    }
#endif
}
#endif // if defined(HAVE_MPI)

void *MallocPlus::memory_malloc(size_t nelem, size_t elsize, const char *name, int flags){
   malloc_plus_memory_entry memory_item;

   memory_item.mem_nelem    = (size_t *)malloc(1*sizeof(size_t));
   memory_item.mem_nelem[0] = nelem;
   memory_item.mem_ndims    = 1;
   memory_item.mem_elsize   = elsize;
   memory_item.mem_flags    = flags;
   if ((flags & DEVICE_REGULAR_MEMORY) != 0){
#ifdef HAVE_OPENCL
      cl_context context = ezcl_get_context();
      memory_item.mem_capacity = nelem;
      memory_item.mem_ptr      = ezcl_device_memory_malloc(context, NULL, name, nelem, elsize, CL_MEM_READ_WRITE, 0);
#endif
   }
   else if ((flags & HOST_MANAGED_MEMORY) != 0){
      memory_item.mem_capacity = 2 * nelem;
      memory_item.mem_ptr      = malloc(2* nelem*elsize);
   }
#ifdef HAVE_J7
   else if (flags & LOAD_BALANCE_MEMORY) {
      memory_item.mem_capacity = nelem;
      memory_item.mem_ptr      = j7->memAlloc(nelem * elsize);
   }
#endif
   else {
      memory_item.mem_capacity = nelem;
      memory_item.mem_ptr      = malloc(nelem*elsize);
   }
   char *mem_name = (char *)malloc(MIN(strlen(name)+1,20));
   strncpy(mem_name,name,MIN(strlen(name),19));
   mem_name[strlen(name)]='\0';
   memory_item.mem_name = mem_name;
   memory_list.push_back(memory_item);
   if (DEBUG) printf("MALLOC_PLUS_MEMORY_MALLOC: DEBUG -- malloc plus memory pointer for %s is %p nelements %ld elsize is %ld\n",mem_name,memory_item.mem_ptr,memory_item.mem_nelem[0],memory_item.mem_elsize);

   return(memory_item.mem_ptr);
}

void *MallocPlus::memory_realloc(size_t nelem, void *malloc_mem_ptr){
   list<malloc_plus_memory_entry>::iterator it;
   void *mem_ptr=NULL;

   for ( it=memory_list.begin(); it != memory_list.end(); it++){
      if (DEBUG) printf("Testing it ptr %p ptr in %p name %s\n",it->mem_ptr,malloc_mem_ptr,it->mem_name);
      if (malloc_mem_ptr == it->mem_ptr) break;
   }
   if (it != memory_list.end()){
      if (DEBUG) printf("MALLOC_PLUS_MEMORY_REALLOC: DEBUG -- reallocated memory pointer %p\n",it->mem_ptr);
      if (it->mem_flags & HOST_MANAGED_MEMORY){
         // Check to see if memory needs to be expanded
         if (nelem > it->mem_capacity) {
            // Need to realloc memory. Allocate extra for growth of array.
            mem_ptr=realloc(it->mem_ptr, 2*nelem*it->mem_elsize);
            it->mem_capacity = 2*nelem;
            it->mem_nelem[0] = nelem;
            it->mem_ptr      = mem_ptr;
         } else {
            // Just move size to use more of memory buffer
            it->mem_nelem[0] = nelem;
         }
      }
#ifdef HAVE_J7
      else if (it->mem_flags & LOAD_BALANCE_MEMORY) {
         mem_ptr = j7->memRealloc(it->mem_ptr, nelem * it->mem_elsize);
         it->mem_capacity = nelem;
         it->mem_nelem[0] = nelem;
         it->mem_ptr      = mem_ptr;
      }
#endif
      else {
         mem_ptr=realloc(it->mem_ptr, nelem*it->mem_elsize);
         it->mem_capacity = nelem;
         it->mem_nelem[0] = nelem;
         it->mem_ptr      = mem_ptr;
      }
   } else {
      if (DEBUG) printf("Warning -- memory pointer %p not found\n",malloc_mem_ptr);
   }
   return(mem_ptr);
}

void *MallocPlus::memory_realloc(size_t nelem, const char *name){
   list<malloc_plus_memory_entry>::iterator it;
   void *mem_ptr=NULL;

   for ( it=memory_list.begin(); it != memory_list.end(); it++){
      if (DEBUG) printf("Testing it name %s input name %s\n",it->mem_name,name);
      if (! strcmp(name,it->mem_name)) break;
   }
   if (it != memory_list.end()){
      if (DEBUG) printf("MALLOC_PLUS_MEMORY_REALLOC: DEBUG -- "
                        "reallocated memory pointer %p\n",it->mem_ptr);
      if (it->mem_flags & HOST_MANAGED_MEMORY) {
         // Check to see if memory needs to be expanded
         if (nelem > it->mem_capacity) {
            // Need to realloc memory. Allocate extra for growth of array.
            mem_ptr=realloc(it->mem_ptr, 2*nelem*it->mem_elsize);
            it->mem_capacity = 2*nelem;
            it->mem_nelem[0] = nelem;
            it->mem_ptr      = mem_ptr;
         } else {
            // Just move size to use more of memory buffer
            it->mem_nelem[0] = nelem;
         }
      }
#ifdef HAVE_J7
      else if (it->mem_flags & LOAD_BALANCE_MEMORY) {
         mem_ptr = j7->memRealloc(it->mem_ptr, nelem * it->mem_elsize);
         it->mem_capacity = nelem;
         it->mem_nelem[0] = nelem;
         it->mem_ptr      = mem_ptr;
      }
#endif
      else {
         mem_ptr=realloc(it->mem_ptr, nelem*it->mem_elsize);
         it->mem_capacity = nelem;
         it->mem_nelem[0] = nelem;
         it->mem_ptr      = mem_ptr;
      }
   } else {
      if (DEBUG) printf("Warning -- memory named %s not found\n",name);
   }
   return(mem_ptr);
}

void *MallocPlus::memory_request(size_t new_capacity, void *malloc_mem_ptr){
   list<malloc_plus_memory_entry>::iterator it;
   void *mem_ptr=NULL;

   for ( it=memory_list.begin(); it != memory_list.end(); it++){
      if (DEBUG) printf("Testing it ptr %p ptr in %p name %s\n",it->mem_ptr,malloc_mem_ptr,it->mem_name);
      if (malloc_mem_ptr == it->mem_ptr) break;
   }
   if (it != memory_list.end()){
      if (DEBUG) printf("MALLOC_PLUS_MEMORY_REQUEST: DEBUG -- reallocated memory pointer %p\n",it->mem_ptr);
      mem_ptr=realloc(it->mem_ptr, new_capacity*it->mem_elsize);
      it->mem_capacity = new_capacity;
      it->mem_ptr      = mem_ptr;
   } else {
      if (DEBUG) printf("Warning -- memory pointer %p not found\n",malloc_mem_ptr);
   }
   return(mem_ptr);
}

void *MallocPlus::memory_request(size_t new_capacity, const char *name){
   list<malloc_plus_memory_entry>::iterator it;
   void *mem_ptr=NULL;

   for ( it=memory_list.begin(); it != memory_list.end(); it++){
      if (DEBUG) printf("Testing it name %s input name %s\n",it->mem_name,name);
      if (! strcmp(name,it->mem_name)) break;
   }
   if (it != memory_list.end()){
      if (DEBUG) printf("MALLOC_PLUS_MEMORY_REQUEST: DEBUG -- reallocated memory pointer %p\n",it->mem_ptr);
      mem_ptr=realloc(it->mem_ptr, new_capacity*it->mem_elsize);
      it->mem_capacity = new_capacity;
      it->mem_ptr      = mem_ptr;
   } else {
      if (DEBUG) printf("Warning -- memory named %s not found\n",name);
   }
   return(mem_ptr);
}

void MallocPlus::memory_realloc_all(size_t nelem){
   list<malloc_plus_memory_entry>::iterator it;
   void *mem_ptr=NULL;

   for ( it=memory_list.begin(); it != memory_list.end(); it++){
      if (it->mem_flags & HOST_MANAGED_MEMORY) {
         if (nelem > it->mem_capacity) {
            mem_ptr=realloc(it->mem_ptr, nelem*it->mem_elsize);
            if (DEBUG) printf("MALLOC_PLUS_MEMORY_REALLOC_ALL: DEBUG -- reallocated memory pointer %p new pointer %p\n",it->mem_ptr,mem_ptr);
            it->mem_capacity = nelem;
            it->mem_nelem[0] = nelem;
            it->mem_ptr      = mem_ptr;
         } else {
            it->mem_nelem[0] = nelem;
         }
      }
#ifdef HAVE_J7
      else if (it->mem_flags & LOAD_BALANCE_MEMORY) {
         mem_ptr = j7->memRealloc(it->mem_ptr, nelem * it->mem_elsize);
         it->mem_capacity = nelem;
         it->mem_nelem[0] = nelem;
         it->mem_ptr      = mem_ptr;
      }
#endif
      else {
         mem_ptr=realloc(it->mem_ptr, nelem*it->mem_elsize);
         if (DEBUG) printf("MALLOC_PLUS_MEMORY_REALLOC_ALL: DEBUG -- reallocated memory pointer %p new pointer %p\n",it->mem_ptr,mem_ptr);
         it->mem_capacity = nelem;
         it->mem_nelem[0] = nelem;
         it->mem_ptr      = mem_ptr;
      }
   }
}

void MallocPlus::memory_request_all(size_t new_capacity){
   list<malloc_plus_memory_entry>::iterator it;
   void *mem_ptr=NULL;

   for ( it=memory_list.begin(); it != memory_list.end(); it++){
      mem_ptr=realloc(it->mem_ptr, new_capacity*it->mem_elsize);
      if (DEBUG) printf("MALLOC_PLUS_MEMORY_REQUEST_ALL: DEBUG -- reallocated memory pointer %p new pointer %p\n",it->mem_ptr,mem_ptr);
      it->mem_capacity = new_capacity;
      it->mem_ptr      = mem_ptr;
   }
}

// This routine is for memory allocated by the host program and added to the database
void *MallocPlus::memory_add(void *malloc_mem_ptr, size_t nelem, size_t elsize, const char *name, int flags){
   malloc_plus_memory_entry memory_item;

   memory_item.mem_nelem    = (size_t *)malloc(1*sizeof(size_t));
   memory_item.mem_nelem[0] = nelem;
   memory_item.mem_ndims    = 1;
   memory_item.mem_capacity = nelem;
   memory_item.mem_elsize   = elsize;
   memory_item.mem_flags    = flags;
   memory_item.mem_ptr      = malloc_mem_ptr;
   char *mem_name = (char *)malloc(MIN(strlen(name)+1,20));
   strncpy(mem_name,name,MIN(strlen(name),19));
   mem_name[MIN(strlen(name),20)]='\0';
   memory_item.mem_name = mem_name;
   memory_list.push_front(memory_item);
   if (DEBUG) printf("MALLOC_PLUS_MEMORY_ADD: DEBUG -- added memory pointer for %s is %p\n",mem_name,malloc_mem_ptr);

   return(malloc_mem_ptr);
}

// This routine is for memory allocated by the host program and added to the database
void *MallocPlus::memory_add(void *malloc_mem_ptr, int ndim, size_t *nelem, size_t elsize, const char *name, int flags){
   malloc_plus_memory_entry memory_item;

   memory_item.mem_nelem    = (size_t *)malloc(ndim*sizeof(size_t));
   for (int i=0; i<ndim; i++){
     memory_item.mem_nelem[i] = nelem[i];
   }
   memory_item.mem_ndims    = ndim;
   memory_item.mem_capacity = 0;
   memory_item.mem_elsize   = elsize;
   memory_item.mem_flags    = flags;
   memory_item.mem_ptr      = malloc_mem_ptr;
   memory_item.mem_name = strdup(name); // mallocs memory
   char *mem_name = (char *)malloc(MIN(strlen(name)+1,20));
   strncpy(mem_name,name,MIN(strlen(name),19));
   mem_name[MIN(strlen(name),20)]='\0';
   memory_item.mem_name = mem_name;
   if (DEBUG) printf("MALLOC_PLUS_MEMORY_ADD: DEBUG -- added memory pointer for %s is %p\n",name,malloc_mem_ptr);

   return(malloc_mem_ptr);
}

double *MallocPlus::memory_reorder(double *malloc_mem_ptr, int *iorder){
   list<malloc_plus_memory_entry>::iterator it;
   double *ptr;

   for ( it=memory_list.begin(); it != memory_list.end(); it++){
      if (DEBUG) printf("Testing it ptr %p ptr in %p name %s\n",it->mem_ptr,malloc_mem_ptr,it->mem_name);
      if (malloc_mem_ptr == it->mem_ptr) break;
   }
   if (it != memory_list.end() ){
      if (DEBUG) printf("Found it ptr %p name %s\n",it->mem_ptr,it->mem_name);
      double *tmp = (double *)malloc(it->mem_nelem[0]*it->mem_elsize);
#ifdef _OPENMP
#pragma omp parallel for
#endif
      for (uint ic = 0; ic < it->mem_nelem[0]; ic++){
         tmp[ic] = malloc_mem_ptr[iorder[ic]];
      }
      SWAP_PTR(malloc_mem_ptr, tmp, ptr);
      free(tmp);
      it->mem_ptr = malloc_mem_ptr;
   } else {
      if (DEBUG) printf("Warning -- memory pointer %p not found\n",malloc_mem_ptr);
   }

   return(malloc_mem_ptr);
}

float *MallocPlus::memory_reorder(float *malloc_mem_ptr, int *iorder){
   list<malloc_plus_memory_entry>::iterator it;
   float *ptr;

   for ( it=memory_list.begin(); it != memory_list.end(); it++){
      if (DEBUG) printf("Testing it ptr %p ptr in %p name %s\n",it->mem_ptr,malloc_mem_ptr,it->mem_name);
      if (malloc_mem_ptr == it->mem_ptr) break;
   }
   if (it != memory_list.end() ){
      if (DEBUG) printf("Found it ptr %p name %s\n",it->mem_ptr,it->mem_name);
      float *tmp = (float *)malloc(it->mem_nelem[0]*it->mem_elsize);
#ifdef _OPENMP
#pragma omp parallel for
#endif
      for (uint ic = 0; ic < it->mem_nelem[0]; ic++){
         tmp[ic] = malloc_mem_ptr[iorder[ic]];
      }
      SWAP_PTR(malloc_mem_ptr, tmp, ptr);
      free(tmp);
      it->mem_ptr = malloc_mem_ptr;
   } else {
      if (DEBUG) printf("Warning -- memory pointer %p not found\n",malloc_mem_ptr);
   }

   return(malloc_mem_ptr);
}

int *MallocPlus::memory_reorder(int *malloc_mem_ptr, int *iorder){
   list<malloc_plus_memory_entry>::iterator it;
   int *ptr;

   for ( it=memory_list.begin(); it != memory_list.end(); it++){
      if (DEBUG) printf("Testing it ptr %p ptr in %p name %s\n",it->mem_ptr,malloc_mem_ptr,it->mem_name);
      if (malloc_mem_ptr == it->mem_ptr) break;
   }
   if (it != memory_list.end() ){
      if (DEBUG) printf("Found it ptr %p name %s\n",it->mem_ptr,it->mem_name);
      int *tmp = (int *)malloc(it->mem_nelem[0]*it->mem_elsize);
#ifdef _OPENMP
#pragma omp parallel for
#endif
      for (uint ic = 0; ic < it->mem_nelem[0]; ic++){
         tmp[ic] = malloc_mem_ptr[iorder[ic]];
      }
      SWAP_PTR(malloc_mem_ptr, tmp, ptr);
      free(tmp);
      it->mem_ptr = malloc_mem_ptr;
   } else {
      if (DEBUG) printf("Warning -- memory pointer %p not found\n",malloc_mem_ptr);
   }

   return(malloc_mem_ptr);
}

int *MallocPlus::memory_reorder_indexarray(int *malloc_mem_ptr, int *iorder, int *inv_iorder){
   list<malloc_plus_memory_entry>::iterator it;
   int *ptr;

   for ( it=memory_list.begin(); it != memory_list.end(); it++){
      if (DEBUG) printf("Testing it ptr %p ptr in %p name %s\n",it->mem_ptr,malloc_mem_ptr,it->mem_name);
      if (malloc_mem_ptr == it->mem_ptr) break;
   }
   if (it != memory_list.end() ){
      if (DEBUG) printf("Found it ptr %p name %s\n",it->mem_ptr,it->mem_name);
      int *tmp = (int *)malloc(it->mem_nelem[0]*it->mem_elsize);
      for (uint ic = 0; ic < it->mem_nelem[0]; ic++){
         tmp[ic] = inv_iorder[malloc_mem_ptr[iorder[ic]]];
      }
      SWAP_PTR(malloc_mem_ptr, tmp, ptr);
      free(tmp);
      it->mem_ptr = malloc_mem_ptr;
   } else {
      if (DEBUG) printf("Warning -- memory pointer %p not found\n",malloc_mem_ptr);
   }

   return(malloc_mem_ptr);
}

void MallocPlus::memory_reorder_all(int *iorder){
   list<malloc_plus_memory_entry>::iterator it;
   vector<int> inv_iorder;

   for ( it=memory_list.begin(); it != memory_list.end(); it++){
      if (it->mem_flags & 0x100) {
         if (inv_iorder.size() < it->mem_nelem[0]) {
            inv_iorder.resize(it->mem_nelem[0]);
            for (int ic = 0; ic < (int)it->mem_nelem[0]; ic++){
               inv_iorder[iorder[ic]] = ic;
            }
         }
         int *ptr;
         int *malloc_mem_ptr = (int *)it->mem_ptr;
         int *tmp = (int *)malloc(it->mem_nelem[0]*it->mem_elsize);
         for (uint ic = 0; ic < it->mem_nelem[0]; ic++){
            tmp[ic] = inv_iorder[malloc_mem_ptr[iorder[ic]]];
         }
         memory_replace(malloc_mem_ptr, tmp);
         SWAP_PTR(malloc_mem_ptr, tmp, ptr);
         free(tmp);
         it->mem_ptr = malloc_mem_ptr;
      } else if (it->mem_elsize == 8){
         double *ptr;
         double *malloc_mem_ptr = (double *)it->mem_ptr;
         double *tmp = (double *)malloc(it->mem_nelem[0]*it->mem_elsize);
         for (uint ic = 0; ic < it->mem_nelem[0]; ic++){
            tmp[ic] = malloc_mem_ptr[iorder[ic]];
         }
         SWAP_PTR(malloc_mem_ptr, tmp, ptr);
         free(tmp);
         it->mem_ptr = malloc_mem_ptr;
      } else {
         float *ptr;
         float *malloc_mem_ptr = (float *)it->mem_ptr;
         float *tmp = (float *)malloc(it->mem_nelem[0]*it->mem_elsize);
         for (uint ic = 0; ic < it->mem_nelem[0]; ic++){
            tmp[ic] = malloc_mem_ptr[iorder[ic]];
         }
         memory_replace(malloc_mem_ptr, tmp);
         SWAP_PTR(malloc_mem_ptr, tmp, ptr);
         free(tmp);
         it->mem_ptr = malloc_mem_ptr;
      }
   }
   inv_iorder.clear();
}

void MallocPlus::memory_report(void){
   list<malloc_plus_memory_entry>::iterator it;

   for ( it=memory_list.begin(); it != memory_list.end(); it++){
      printf("MallocPlus %10s ptr %p dims %lu nelem (",
            it->mem_name,it->mem_ptr,it->mem_ndims);

      char nelemstring[80];
      char *str_ptr = nelemstring;
      str_ptr += sprintf(str_ptr,"%lu", it->mem_nelem[0]);
      for (uint i = 1; i < it->mem_ndims; i++){
         str_ptr += sprintf(str_ptr,", %lu", it->mem_nelem[i]);
      }
      printf("%12s",nelemstring);

      printf(") elsize %lu flags %d capacity %lu\n",
            it->mem_elsize,it->mem_flags,it->mem_capacity);
   }
}

void *MallocPlus::memory_delete(void *malloc_mem_ptr){
   list<malloc_plus_memory_entry>::iterator it;

   for ( it=memory_list.begin(); it != memory_list.end(); it++){
      if (DEBUG) printf("Testing it ptr %p ptr in %p name %s\n",it->mem_ptr,malloc_mem_ptr,it->mem_name);
      if (malloc_mem_ptr == it->mem_ptr) break;
   }
   if (it != memory_list.end()){
      if (DEBUG) printf("MALLOC_PLUS_MEMORY_REMOVE: DEBUG -- removed memory pointer %p\n",it->mem_ptr);
      free(it->mem_name);
      if ((it->mem_flags & DEVICE_REGULAR_MEMORY) != 0){
#ifdef HAVE_OPENCL
         //printf("MALLOC_PLUS_MEMORY_REMOVE: DEBUG -- removed memory pointer %p\n",it->mem_ptr);
         ezcl_device_memory_delete(it->mem_ptr);
#endif
      }
#ifdef HAVE_J7
      else if (it->mem_flags & LOAD_BALANCE_MEMORY) {
         j7->memFree(it->mem_ptr);
      }
#endif
      else {
         free(it->mem_ptr);
      }
      memory_list.erase(it);
   } else {
      if (DEBUG) printf("Warning -- memory pointer %p not found\n",malloc_mem_ptr);
   }

   return(NULL);
}

void *MallocPlus::memory_delete(const char *name){
   list<malloc_plus_memory_entry>::iterator it;

   for ( it=memory_list.begin(); it != memory_list.end(); it++){
      if (DEBUG) printf("Testing it name %s input name %s\n",it->mem_name,name);
      if (! strcmp(name,it->mem_name)) break;
   }
   if (it != memory_list.end()){
      if (DEBUG) printf("MALLOC_PLUS_MEMORY_REMOVE: DEBUG -- removed memory pointer %p\n",it->mem_ptr);
      free(it->mem_name);
      if ((it->mem_flags & DEVICE_REGULAR_MEMORY) != 0){
#ifdef HAVE_OPENCL
         ezcl_device_memory_delete(it->mem_ptr);
#endif
      }
#ifdef HAVE_J7
      else if (it->mem_flags & LOAD_BALANCE_MEMORY) {
         j7->memFree(it->mem_ptr);
      }
#endif
      else {
         free(it->mem_ptr);
      }
      memory_list.erase(it);
   } else {
      if (DEBUG) printf("Warning -- memory named %s not found\n",name);
   }

   return(NULL);
}

void MallocPlus::memory_delete_all(void){
   list<malloc_plus_memory_entry>::iterator it;
   list<malloc_plus_memory_entry> memory_list_old = memory_list;

   for ( it=memory_list_old.begin(); it != memory_list_old.end(); it++){
      if (DEBUG) printf("MALLOC_PLUS_MEMORY_REMOVE: DEBUG -- removed memory pointer %p name %s\n",it->mem_ptr,it->mem_name);
      free(it->mem_name);
      if ((it->mem_flags & DEVICE_REGULAR_MEMORY) != 0){
#ifdef HAVE_OPENCL
         ezcl_device_memory_delete(it->mem_ptr);
#endif
      } else {
         free(it->mem_ptr);
      }
   }
   memory_list.clear();
}

void MallocPlus::memory_remove(void *malloc_mem_ptr){
   list<malloc_plus_memory_entry>::iterator it;

   for ( it=memory_list.begin(); it != memory_list.end(); it++){
      if (DEBUG) printf("Testing it ptr %p ptr in %p name %s\n",it->mem_ptr,malloc_mem_ptr,it->mem_name);
      if (malloc_mem_ptr == it->mem_ptr) break;
   }
   if (it != memory_list.end()){
      if (DEBUG) printf("MALLOC_PLUS_MEMORY_REMOVE: DEBUG -- removed memory pointer %p\n",it->mem_ptr);
      free(it->mem_name);
      memory_list.erase(it);
   } else {
      if (DEBUG) printf("Warning -- memory pointer %p not found\n",malloc_mem_ptr);
   }
}

void MallocPlus::memory_remove(const char *name){
   list<malloc_plus_memory_entry>::iterator it;

   for ( it=memory_list.begin(); it != memory_list.end(); it++){
      if (DEBUG) printf("Testing it name %s input name %s\n",it->mem_name,name);
      if (! strcmp(name,it->mem_name)) break;
   }
   if (it != memory_list.end()){
      if (DEBUG) printf("MALLOC_PLUS_MEMORY_REMOVE: DEBUG -- removed memory pointer %p\n",it->mem_ptr);
      free(it->mem_name);
      memory_list.erase(it);
   } else {
      if (DEBUG) printf("Warning -- memory named %s not found\n",name);
   }
}

void *MallocPlus::memory_begin(void){
   it_save = memory_list.begin();
   return(it_save->mem_ptr);
}

void *MallocPlus::memory_next(void){
   list<malloc_plus_memory_entry>::iterator it;

   it = it_save;
   it++;

   if (it != memory_list.end()){
      if (DEBUG) printf("Found it ptr %p name %s\n",it->mem_ptr,it->mem_name);
      it_save = it;
      return(it->mem_ptr);
   } else {
      if (DEBUG) printf("Warning -- memory not found\n");
      return(NULL);
   }
}

list<malloc_plus_memory_entry>::iterator MallocPlus::memory_entry_begin(void){
   it_save = memory_list.begin();
   return(it_save);
}

list<malloc_plus_memory_entry>::iterator MallocPlus::memory_entry_next(void){
   list<malloc_plus_memory_entry>::iterator it;

   it = it_save;
   it++;

   if (it != memory_list.end()){
      if (DEBUG) printf("Found it ptr %p name %s\n",it->mem_ptr,it->mem_name);
      it_save = it;
      return(it);
   } else {
      if (DEBUG) printf("Warning -- memory not found\n");
      return((list<malloc_plus_memory_entry>::iterator) NULL);
   }
}

size_t MallocPlus::get_memory_size(void *malloc_mem_ptr){
   list<malloc_plus_memory_entry>::iterator it;

   for ( it=memory_list.begin(); it != memory_list.end(); it++){
      if (DEBUG) printf("Testing it ptr %p ptr in %p name %s\n",it->mem_ptr,malloc_mem_ptr,it->mem_name);
      if (malloc_mem_ptr == it->mem_ptr) break;
   }
   if (it != memory_list.end()){
      if (DEBUG) printf("Found it ptr %p name %s\n",it->mem_ptr,it->mem_name);
      return(it->mem_nelem[0]);
   } else {
      if (DEBUG) printf("Warning -- memory not found\n");
   }
   return(0);
}

int MallocPlus::get_memory_elemsize(void *malloc_mem_ptr){
   list<malloc_plus_memory_entry>::iterator it;

   for ( it=memory_list.begin(); it != memory_list.end(); it++){
      if (DEBUG) printf("Testing it ptr %p ptr in %p name %s\n",it->mem_ptr,malloc_mem_ptr,it->mem_name);
      if (malloc_mem_ptr == it->mem_ptr) break;
   }
   if (it != memory_list.end()){
      if (DEBUG) printf("Found it ptr %p name %s\n",it->mem_ptr,it->mem_name);
      return(it->mem_elsize);
   } else {
      if (DEBUG) printf("Warning -- memory not found\n");
   }
   return(0);
}

int MallocPlus::get_memory_flags(void *malloc_mem_ptr){
   list<malloc_plus_memory_entry>::iterator it;

   for ( it=memory_list.begin(); it != memory_list.end(); it++){
      if (DEBUG) printf("Testing it ptr %p ptr in %p name %s\n",it->mem_ptr,malloc_mem_ptr,it->mem_name);
      if (malloc_mem_ptr == it->mem_ptr) break;
   }
   if (it != memory_list.end()){
      if (DEBUG) printf("Found it ptr %p name %s attribute %d\n",it->mem_ptr,it->mem_name,it->mem_flags);
      return(it->mem_flags);
   } else {
      if (DEBUG) printf("Warning -- memory not found\n");
   }
   return(0);
}

size_t MallocPlus::get_memory_capacity(void *malloc_mem_ptr){
   list<malloc_plus_memory_entry>::iterator it;

   for ( it=memory_list.begin(); it != memory_list.end(); it++){
      if (DEBUG) printf("Testing it ptr %p ptr in %p name %s\n",it->mem_ptr,malloc_mem_ptr,it->mem_name);
      if (malloc_mem_ptr == it->mem_ptr) break;
   }
   if (it != memory_list.end()){
      if (DEBUG) printf("Found it ptr %p name %s\n",it->mem_ptr,it->mem_name);
      return(it->mem_capacity);
   } else {
      if (DEBUG) printf("Warning -- memory not found\n");
   }
   return(0);
}

const char * MallocPlus::get_memory_name(void *malloc_mem_ptr){
   list<malloc_plus_memory_entry>::iterator it;

   for ( it=memory_list.begin(); it != memory_list.end(); it++){
      if (DEBUG) printf("Testing it ptr %p ptr in %p name %s\n",it->mem_ptr,malloc_mem_ptr,it->mem_name);
      if (malloc_mem_ptr == it->mem_ptr) break;
   }
   if (it != memory_list.end()){
      if (DEBUG) printf("Found it ptr %p name %s\n",it->mem_ptr,it->mem_name);
      return(it->mem_name);
   } else {
      if (DEBUG) printf("Warning -- memory not found\n");
   }
   return(NULL);
}

void *MallocPlus::memory_replace(void *malloc_mem_ptr_old, void * const malloc_mem_ptr_new){
   list<malloc_plus_memory_entry>::iterator it, it_old=memory_list.end(), it_new=memory_list.end();

   for ( it=memory_list.begin(); it != memory_list.end(); it++){
      if (DEBUG) printf("Testing it ptr %p ptr old %p ptr new %p name %s\n",it->mem_ptr,malloc_mem_ptr_old,malloc_mem_ptr_new,it->mem_name);
      if (malloc_mem_ptr_old == it->mem_ptr) it_old = it;
      if (malloc_mem_ptr_new == it->mem_ptr) it_new = it;
      if (it_new != memory_list.end() && it_old != memory_list.end() ) break;
   }
   if (it != memory_list.end()){
      if (DEBUG) printf("Found it ptr_old %p name %s ptr_new %p name %s\n",it_old->mem_ptr,it_old->mem_name,it_new->mem_ptr,it_new->mem_name);
      if ((it_old->mem_flags & DEVICE_REGULAR_MEMORY) != 0){
#ifdef HAVE_OPENCL
         if (DEBUG) printf("Deleting device memory name %s pointer %p\n",it_old->mem_name,it_old->mem_ptr);
         ezcl_device_memory_delete(it_old->mem_ptr);
#endif
      }
#ifdef HAVE_J7
      else if (it->mem_flags & LOAD_BALANCE_MEMORY) {
         j7->memFree(it_old->mem_ptr);
      }
#endif
      else {
         free(it_old->mem_ptr);
      }
      it_old->mem_ptr      = it_new->mem_ptr;
      it_old->mem_nelem[0] = it_new->mem_nelem[0];
      it_old->mem_capacity = it_new->mem_capacity;
      it_old->mem_elsize   = it_new->mem_elsize;
      it_old->mem_flags    = it_new->mem_flags;
      malloc_mem_ptr_old = (void *)malloc_mem_ptr_new;
      free(it_new->mem_name);
      memory_list.erase(it_new);
      return(it_old->mem_ptr);
   } else {
      if (DEBUG) printf("Warning -- memory not found\n");
   }
   return(NULL);
}

void MallocPlus::memory_swap(int **malloc_mem_ptr_old, int **malloc_mem_ptr_new){
   list<malloc_plus_memory_entry>::iterator it, it_old=memory_list.end(), it_new=memory_list.end();
   malloc_plus_memory_entry it_tmp;

   for ( it=memory_list.begin(); it != memory_list.end(); it++){
      if (DEBUG) printf("Testing it ptr %p ptr old %p ptr new %p name %s\n",it->mem_ptr,malloc_mem_ptr_old,malloc_mem_ptr_new,it->mem_name);
      if (*malloc_mem_ptr_old == it->mem_ptr) it_old = it;
      if (*malloc_mem_ptr_new == it->mem_ptr) it_new = it;
      if (it_new != memory_list.end() && it_old != memory_list.end() ) break;
   }
   if (it != memory_list.end()){
      if (DEBUG) printf("Found it ptr_old %p name %s ptr_new %p name %s\n",it_old->mem_ptr,it_old->mem_name,it_new->mem_ptr,it_new->mem_name);

      it_tmp.mem_name  = it_old->mem_name;
      it_old->mem_name = it_new->mem_name;
      it_new->mem_name = it_tmp.mem_name;

      *malloc_mem_ptr_old = (int *)it_new->mem_ptr;
      *malloc_mem_ptr_new = (int *)it_old->mem_ptr;
   } else {
      if (DEBUG) printf("Warning -- memory not found\n");
   }
}

void MallocPlus::memory_swap(float **malloc_mem_ptr_old, float **malloc_mem_ptr_new){
   list<malloc_plus_memory_entry>::iterator it, it_old=memory_list.end(), it_new=memory_list.end();
   malloc_plus_memory_entry it_tmp;

   for ( it=memory_list.begin(); it != memory_list.end(); it++){
      if (DEBUG) printf("Testing it ptr %p ptr old %p ptr new %p name %s\n",it->mem_ptr,malloc_mem_ptr_old,malloc_mem_ptr_new,it->mem_name);
      if (*malloc_mem_ptr_old == it->mem_ptr) it_old = it;
      if (*malloc_mem_ptr_new == it->mem_ptr) it_new = it;
      if (it_new != memory_list.end() && it_old != memory_list.end() ) break;
   }
   if (it != memory_list.end()){
      if (DEBUG) printf("Found it ptr_old %p name %s ptr_new %p name %s\n",it_old->mem_ptr,it_old->mem_name,it_new->mem_ptr,it_new->mem_name);

      it_tmp.mem_name  = it_old->mem_name;
      it_old->mem_name = it_new->mem_name;
      it_new->mem_name = it_tmp.mem_name;

      *malloc_mem_ptr_old = (float *)it_new->mem_ptr;
      *malloc_mem_ptr_new = (float *)it_old->mem_ptr;
   } else {
      if (DEBUG) printf("Warning -- memory not found\n");
   }
}

void MallocPlus::memory_swap(double **malloc_mem_ptr_old, double **malloc_mem_ptr_new){
   list<malloc_plus_memory_entry>::iterator it, it_old=memory_list.end(), it_new=memory_list.end();
   malloc_plus_memory_entry it_tmp;

   for ( it=memory_list.begin(); it != memory_list.end(); it++){
      if (DEBUG) printf("Testing it ptr %p ptr old %p ptr new %p name %s\n",it->mem_ptr,malloc_mem_ptr_old,malloc_mem_ptr_new,it->mem_name);
      if (*malloc_mem_ptr_old == it->mem_ptr) it_old = it;
      if (*malloc_mem_ptr_new == it->mem_ptr) it_new = it;
      if (it_new != memory_list.end() && it_old != memory_list.end() ) break;
   }
   if (it != memory_list.end()){
      if (DEBUG) printf("Found it ptr_old %p name %s ptr_new %p name %s\n",it_old->mem_ptr,it_old->mem_name,it_new->mem_ptr,it_new->mem_name);

      it_tmp.mem_name  = it_old->mem_name;
      it_old->mem_name = it_new->mem_name;
      it_new->mem_name = it_tmp.mem_name;

      *malloc_mem_ptr_old = (double *)it_new->mem_ptr;
      *malloc_mem_ptr_new = (double *)it_old->mem_ptr;
   } else {
      if (DEBUG) printf("Warning -- memory not found\n");
   }
}

void *MallocPlus::memory_duplicate(void *malloc_mem_ptr, const char *addname){
   list<malloc_plus_memory_entry>::iterator it;
   void *mem_ptr_dup;

   for ( it=memory_list.begin(); it != memory_list.end(); it++){
      if (DEBUG) printf("Testing it ptr %p ptr %p name %s\n",it->mem_ptr,malloc_mem_ptr,it->mem_name);
      if (malloc_mem_ptr == it->mem_ptr) break;
   }
   if (it != memory_list.end()){
      int str_size = strlen(it->mem_name) + strlen(addname)+2;
      char *new_name = (char *)malloc(str_size*sizeof(char));
      new_name = strcat(new_name, addname);
      new_name = strcat(new_name, it->mem_name);
      mem_ptr_dup   = memory_malloc(it->mem_nelem[0], it->mem_elsize, new_name, it->mem_flags);
      free(new_name);
      mem_ptr_dup = memcpy(mem_ptr_dup, malloc_mem_ptr, it->mem_nelem[0]*it->mem_elsize);
      return(mem_ptr_dup);
   } else {
      if (DEBUG) printf("Warning -- memory not found\n");
   }
   return(NULL);
}

void *MallocPlus::get_memory_ptr(const char *name){
   list<malloc_plus_memory_entry>::iterator it;

   for ( it=memory_list.begin(); it != memory_list.end(); it++){
      if (DEBUG) printf("Testing it name %s input name %s\n",it->mem_name,name);
      if (! strcmp(name,it->mem_name) ) break;
   }
   if (it != memory_list.end()){
      if (DEBUG) printf("Found it ptr %p name %s\n",it->mem_ptr,it->mem_name);
      return(it->mem_ptr);
   } else {
      if (DEBUG) printf("Warning -- memory not found\n");
   }
   return(NULL);
}

bool MallocPlus::check_memory_attribute(void *malloc_mem_ptr, int attribute){
   list<malloc_plus_memory_entry>::iterator it;

   for ( it=memory_list.begin(); it != memory_list.end(); it++){
      if (DEBUG) printf("Testing it ptr %p ptr in %p name %s\n",it->mem_ptr,malloc_mem_ptr,it->mem_name);
      if (malloc_mem_ptr == it->mem_ptr) break;
   }
   if (it == memory_list.end()){
      printf("Error -- memory not found\n");
      exit(1);
   }

   if (DEBUG) printf("Found it ptr %p name %s attribute %d\n",it->mem_ptr,it->mem_name,it->mem_flags);
   bool bvalue = false;
   if (it->mem_flags & attribute) bvalue = true;

   return bvalue;
}

void MallocPlus::set_memory_attribute(void *malloc_mem_ptr, int attribute){
   list<malloc_plus_memory_entry>::iterator it;

   for ( it=memory_list.begin(); it != memory_list.end(); it++){
      if (DEBUG) printf("Testing it ptr %p ptr in %p name %s\n",it->mem_ptr,malloc_mem_ptr,it->mem_name);
      if (malloc_mem_ptr == it->mem_ptr) break;
   }
   if (it != memory_list.end()){
      if (DEBUG) printf("Found it ptr %p name %s attribute %d\n",it->mem_ptr,it->mem_name,it->mem_flags);
      it->mem_flags |= attribute;
   } else {
      if (DEBUG) printf("Warning -- memory not found\n");
   }
}

void MallocPlus::clear_memory_attribute(void *malloc_mem_ptr, int attribute){
   list<malloc_plus_memory_entry>::iterator it;

   for ( it=memory_list.begin(); it != memory_list.end(); it++){
      if (DEBUG) printf("Testing it ptr %p ptr in %p name %s\n",it->mem_ptr,malloc_mem_ptr,it->mem_name);
      if (malloc_mem_ptr == it->mem_ptr) break;
   }
   if (it != memory_list.end()){
      if (DEBUG) printf("Found it ptr %p name %s attribute %d\n",it->mem_ptr,it->mem_name,it->mem_flags);
      it->mem_flags &= attribute;
   } else {
      if (DEBUG) printf("Warning -- memory not found\n");
   }
}

extern "C" {
   MallocPlus *MallocPlus_new(){
     return new MallocPlus;
   }

   void MallocPlus_memory_report(MallocPlus *mem_object) {
      mem_object->memory_report();
   }

   void MallocPlus_memory_add(MallocPlus *mem_object, void *dbleptr, size_t nelem,
       size_t elsize, char *name, unsigned long long flags){
//   printf("DEBUG -- nelem %lu elsize %lu\n", nelem, elsize);
     mem_object->memory_add(dbleptr, nelem, elsize, name,
       (unsigned long long)flags);
   }
   void MallocPlus_memory_add_nD(MallocPlus *mem_object, void *dbleptr, int ndim, size_t *nelem,
       size_t elsize, char *name, unsigned long long flags){
//   printf("DEBUG -- ndim %d nelem[0] %lu elsize %lu\n",ndim, nelem[0], elsize);
     mem_object->memory_add(dbleptr, ndim, nelem, elsize, name,
       (unsigned long long)flags);
   }
}
