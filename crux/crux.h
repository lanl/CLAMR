/*
 *  Copyright (c) 2014, Los Alamos National Security, LLC.
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
 *  
 *  Authors: Brian Atkinson          bwa@g.clemson.edu
             Bob Robey        XCP-2  brobey@lanl.gov
 */

#ifndef CRUX_H_
#define CRUX_H_

#include <stdio.h>

enum crux_types{
   CRUX_NONE,
   CRUX_DISK,
   CRUX_IN_MEMORY
};

class Crux
{
   int num_of_rollback_states;
   int crux_type;
   int checkpoint_counter;

public:

   Crux(int crux_type_in, int num_of_rollback_states_in, bool restart);
   ~Crux();

   void store_begin(size_t nsize, int ncycle);
   void store_ints(int *int_vals, size_t nelem);
   void store_longs(long long *long_vals, size_t nelem);
   void store_bools(bool *bool_vals, size_t nelem);
   void store_doubles(double *double_vals, size_t nelem);
   void store_int_array(int *int_array, size_t nelem);
   void store_long_array(long long *long_array, size_t nelem);
   void store_float_array(float *float_array, size_t nelem);
   void store_double_array(double *double_array, size_t nelem);
   void store_end(void);

   void       restore_begin(char *restart_file, int rollback_counter);
   void       restore_ints(int *int_vals, size_t nelem);
   void       restore_bools(bool *bool_vals, size_t nelem);
   void       restore_longs(long long *long_vals, size_t nelem);
   void       restore_doubles(double *double_vals, size_t nelem);
   int       *restore_int_array(int *int_array, size_t nsize);
   long long *restore_long_array(long long *long_array, size_t nsize);
   float     *restore_float_array(float *float_array, size_t nsize);
   double    *restore_double_array(double *double_array, size_t nsize);
   void       restore_end(void);

   int get_rollback_number();

};
#endif // CRUX_H_
