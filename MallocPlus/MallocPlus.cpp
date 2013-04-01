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
#include <unistd.h>
#include <stdio.h>
#include <algorithm>
#include <queue>
#include "string.h"
#ifdef HAVE_OPENCL
#include "ezcl/ezcl.h"
#endif
#include "MallocPlus.h"

//#undef DEBUG
//#define DEBUG 1
#ifndef DEBUG
#define DEBUG 0
#endif

#ifdef HAVE_CL_DOUBLE
typedef double      real;
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

void *MallocPlus::memory_malloc(size_t nelem, size_t elsize, const char *name){
   malloc_plus_memory_entry memory_item;

   memory_item.mem_nelem    = nelem;
   memory_item.mem_capacity = nelem;
   memory_item.mem_elsize   = elsize;
   memory_item.mem_flags    = HOST_REGULAR_MEMORY;
   memory_item.mem_ptr      = malloc(nelem*elsize);
   char *mem_name = (char *)malloc(MIN(strlen(name)+1,20));
   strncpy(mem_name,name,MIN(strlen(name),19));
   mem_name[strlen(name)]='\0';
   memory_item.mem_name = mem_name;
   memory_list.push_back(memory_item);
   if (DEBUG) printf("MALLOC_PLUS_MEMORY_MALLOC: DEBUG -- malloc plus memory pointer for %s is %p nelements %ld elsize is %ld\n",mem_name,memory_item.mem_ptr,memory_item.mem_nelem,memory_item.mem_elsize);

   return(memory_item.mem_ptr);
}

void *MallocPlus::memory_malloc(size_t nelem, size_t elsize, int flags, const char *name){
   malloc_plus_memory_entry memory_item;

   memory_item.mem_nelem  = nelem;
   memory_item.mem_elsize = elsize;
   memory_item.mem_flags  = flags;
   if ((flags & DEVICE_REGULAR_MEMORY) != 0){
#ifdef HAVE_OPENCL
      cl_context context = ezcl_get_context();
      memory_item.mem_capacity = nelem;
      memory_item.mem_ptr      = ezcl_device_memory_malloc(context, NULL, name, nelem, elsize, CL_MEM_READ_WRITE, 0);
#endif
   } else if ((flags & HOST_MANAGED_MEMORY) != 0){
      memory_item.mem_capacity = nelem;
      memory_item.mem_ptr      = malloc(nelem*elsize);
   } else {
      memory_item.mem_capacity = 2*nelem;
      memory_item.mem_ptr      = malloc(2*nelem*elsize);
   }
   char *mem_name = (char *)malloc(MIN(strlen(name)+1,20));
   strncpy(mem_name,name,MIN(strlen(name),19));
   mem_name[strlen(name)]='\0';
   memory_item.mem_name = mem_name;
   memory_list.push_back(memory_item);
   if (DEBUG) printf("MALLOC_PLUS_MEMORY_MALLOC: DEBUG -- malloc plus memory pointer for %s is %p nelements %ld elsize is %ld\n",mem_name,memory_item.mem_ptr,memory_item.mem_nelem,memory_item.mem_elsize);

   return(memory_item.mem_ptr);
}

void *MallocPlus::memory_realloc(size_t nelem, size_t elsize, void *malloc_mem_ptr){
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
            mem_ptr=realloc(it->mem_ptr, 2*nelem*elsize);
            it->mem_capacity = 2*nelem;
            it->mem_nelem    = nelem;
            it->mem_elsize   = elsize;
            it->mem_ptr      = mem_ptr;
         } else {
            // Just move size to use more of memory buffer
            it->mem_nelem  = nelem;
         }
      } else {
         mem_ptr=realloc(it->mem_ptr, nelem*elsize);
         it->mem_capacity = nelem;
         it->mem_nelem    = nelem;
         it->mem_elsize   = elsize;
         it->mem_ptr      = mem_ptr;
      }
   } else {
      if (DEBUG) printf("Warning -- memory pointer %p not found\n",malloc_mem_ptr);
   }
   return(mem_ptr);
}

void *MallocPlus::memory_realloc(size_t nelem, size_t elsize, const char *name){
   list<malloc_plus_memory_entry>::iterator it;
   void *mem_ptr=NULL;

   for ( it=memory_list.begin(); it != memory_list.end(); it++){
      if (DEBUG) printf("Testing it name %s input name %s\n",it->mem_name,name);
      if (! strcmp(name,it->mem_name)) break;
   }
   if (it != memory_list.end()){
      if (DEBUG) printf("MALLOC_PLUS_MEMORY_REALLOC: DEBUG -- reallocated memory pointer %p\n",it->mem_ptr);
      if (it->mem_flags & HOST_MANAGED_MEMORY){
         // Check to see if memory needs to be expanded
         if (nelem > it->mem_capacity) {
            // Need to realloc memory. Allocate extra for growth of array.
            mem_ptr=realloc(it->mem_ptr, 2*nelem*elsize);
            it->mem_capacity = 2*nelem;
            it->mem_nelem    = nelem;
            it->mem_elsize   = elsize;
            it->mem_ptr      = mem_ptr;
         } else {
            // Just move size to use more of memory buffer
            it->mem_nelem  = nelem;
         }
      } else {
         mem_ptr=realloc(it->mem_ptr, nelem*elsize);
         it->mem_capacity = nelem;
         it->mem_nelem    = nelem;
         it->mem_elsize   = elsize;
         it->mem_ptr      = mem_ptr;
      }
   } else {
      if (DEBUG) printf("Warning -- memory named %s not found\n",name);
   }
   return(mem_ptr);
}

void *MallocPlus::memory_request(size_t nelem, size_t elsize, void *malloc_mem_ptr){
   list<malloc_plus_memory_entry>::iterator it;
   void *mem_ptr=NULL;

   for ( it=memory_list.begin(); it != memory_list.end(); it++){
      if (DEBUG) printf("Testing it ptr %p ptr in %p name %s\n",it->mem_ptr,malloc_mem_ptr,it->mem_name);
      if (malloc_mem_ptr == it->mem_ptr) break;
   }
   if (it != memory_list.end()){
      if (DEBUG) printf("MALLOC_PLUS_MEMORY_REQUEST: DEBUG -- reallocated memory pointer %p\n",it->mem_ptr);
      mem_ptr=realloc(it->mem_ptr, nelem*elsize);
      it->mem_capacity = nelem;
      it->mem_elsize   = elsize;
      it->mem_ptr      = mem_ptr;
   } else {
      if (DEBUG) printf("Warning -- memory pointer %p not found\n",malloc_mem_ptr);
   }
   return(mem_ptr);
}

void *MallocPlus::memory_request(size_t nelem, size_t elsize, const char *name){
   list<malloc_plus_memory_entry>::iterator it;
   void *mem_ptr=NULL;

   for ( it=memory_list.begin(); it != memory_list.end(); it++){
      if (DEBUG) printf("Testing it name %s input name %s\n",it->mem_name,name);
      if (! strcmp(name,it->mem_name)) break;
   }
   if (it != memory_list.end()){
      if (DEBUG) printf("MALLOC_PLUS_MEMORY_REQUEST: DEBUG -- reallocated memory pointer %p\n",it->mem_ptr);
      mem_ptr=realloc(it->mem_ptr, nelem*elsize);
      it->mem_capacity = nelem;
      it->mem_elsize   = elsize;
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
            it->mem_nelem = nelem;
            it->mem_ptr   = mem_ptr;
         } else {
            it->mem_nelem = nelem;
         }
      } else {
         mem_ptr=realloc(it->mem_ptr, nelem*it->mem_elsize);
         if (DEBUG) printf("MALLOC_PLUS_MEMORY_REALLOC_ALL: DEBUG -- reallocated memory pointer %p new pointer %p\n",it->mem_ptr,mem_ptr);
         it->mem_nelem = nelem;
         it->mem_ptr   = mem_ptr;
      }
   }
}

void MallocPlus::memory_request_all(size_t nelem){
   list<malloc_plus_memory_entry>::iterator it;
   void *mem_ptr=NULL;

   for ( it=memory_list.begin(); it != memory_list.end(); it++){
      mem_ptr=realloc(it->mem_ptr, nelem*it->mem_elsize);
      if (DEBUG) printf("MALLOC_PLUS_MEMORY_REQUEST_ALL: DEBUG -- reallocated memory pointer %p new pointer %p\n",it->mem_ptr,mem_ptr);
      it->mem_capacity = nelem;
      it->mem_ptr      = mem_ptr;
   }
}

void *MallocPlus::memory_add(void *malloc_mem_ptr, size_t nelem, size_t elsize, const char *name){
   malloc_plus_memory_entry memory_item;

   memory_item.mem_nelem    = nelem;
   memory_item.mem_capacity = nelem;
   memory_item.mem_elsize   = elsize;
   memory_item.mem_flags    = HOST_REGULAR_MEMORY;
   memory_item.mem_ptr      = malloc_mem_ptr;
   char *mem_name = (char *)malloc(MIN(strlen(name)+1,20));
   strncpy(mem_name,name,MIN(strlen(name),19));
   mem_name[MIN(strlen(name),20)]='\0';
   memory_item.mem_name = mem_name;
   memory_list.push_front(memory_item);
   if (DEBUG) printf("MALLOC_PLUS_MEMORY_ADD: DEBUG -- added memory pointer for %s is %p\n",mem_name,malloc_mem_ptr);

   return(malloc_mem_ptr);
}

void *MallocPlus::memory_add(void *malloc_mem_ptr, size_t nelem, size_t elsize, int flags, const char *name){
   malloc_plus_memory_entry memory_item;

   memory_item.mem_nelem    = nelem;
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

real *MallocPlus::memory_reorder(real *malloc_mem_ptr, int *iorder){
   list<malloc_plus_memory_entry>::iterator it;
   real *ptr;

   for ( it=memory_list.begin(); it != memory_list.end(); it++){
      if (DEBUG) printf("Testing it ptr %p ptr in %p name %s\n",it->mem_ptr,malloc_mem_ptr,it->mem_name);
      if (malloc_mem_ptr == it->mem_ptr) break;
   }
   if (it != memory_list.end() ){
      if (DEBUG) printf("Found it ptr %p name %s\n",it->mem_ptr,it->mem_name);
      real *tmp = (real *)malloc(it->mem_nelem*it->mem_elsize);
      for (uint ic = 0; ic < it->mem_nelem; ic++){
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

void MallocPlus::memory_report(void){
   list<malloc_plus_memory_entry>::iterator it;
   for ( it=memory_list.begin(); it != memory_list.end(); it++){
      printf("MallocPlus %10s ptr %p nelements %lu elsize %lu flags %d capacity %lu\n",it->mem_name,it->mem_ptr,it->mem_nelem,it->mem_elsize,it->mem_flags,it->mem_capacity);
   }
}

void MallocPlus::memory_delete(void *malloc_mem_ptr){
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
      } else {
         free(it->mem_ptr);
      }
      memory_list.erase(it);
   } else {
      if (DEBUG) printf("Warning -- memory pointer %p not found\n",malloc_mem_ptr);
   }
}

void MallocPlus::memory_delete(const char *name){
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
      } else {
         free(it->mem_ptr);
      }
      memory_list.erase(it);
   } else {
      if (DEBUG) printf("Warning -- memory named %s not found\n",name);
   }
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
   return(memory_list.begin());
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
      return(it->mem_nelem);
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
      } else {
         free(it_old->mem_ptr);
      }
      it_old->mem_ptr      = it_new->mem_ptr;
      it_old->mem_nelem    = it_new->mem_nelem;
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

