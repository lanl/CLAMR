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

__constant ulong prime = 4294967291;
__constant uint hash_jump_prime = 41;

enum choose_hash_method
{  METHOD_UNSET = 0,            //  use 0 for no method set
   PERFECT_HASH,                //  perfect hash 1
   LINEAR,                      //  linear hash 2
   QUADRATIC,                   //  quadratic hash 3
   PRIME_JUMP  };               //  prime_jump hash 4

void write_hash(
            const int   hash_method,
            const ulong hash_table_size,
            const ulong AA,
            const ulong BB,
            const uint  giX,
            const ulong hashkey,
   __global       int   *hash);

int read_hash(
            const int   hash_method,
            const ulong hash_table_size,
            const ulong AA,
            const ulong BB,
            const ulong hashkey,
   __global const int   *hash);

void write_hash(
            const int   hash_method,
            const ulong hash_table_size,
            const ulong AA,
            const ulong BB,
            const uint  giX,
            const ulong hashkey,
   __global       int   *hash)
{
   uint hashloc;
   int icount = 0;
   uint jump;
   int old_key;
   int MaxTries = 1000;

   switch (hash_method) {
   case PERFECT_HASH:
      hash[hashkey] = giX;
      break;
   case LINEAR:
      hashloc = (hashkey*AA+BB)%prime%hash_table_size;
      old_key = atomic_cmpxchg(&hash[2*hashloc],-1,hashkey);

      for (icount = 1; old_key != hashkey && old_key != -1 && icount < MaxTries; icount++){
         hashloc++;
         hashloc %= hash_table_size;
      
         old_key = atomic_cmpxchg(&hash[2*hashloc],-1,hashkey);
      }

      if (icount < MaxTries) hash[2*hashloc+1] = giX;
      break;
   case QUADRATIC:
      hashloc = (hashkey*AA+BB)%prime%hash_table_size;
      old_key = atomic_cmpxchg(&hash[2*hashloc],-1,hashkey);

      for (icount = 1; old_key != hashkey && old_key != -1 && icount < MaxTries; icount++){
         hashloc+=(icount*icount);
         hashloc %= hash_table_size;
      
         old_key = atomic_cmpxchg(&hash[2*hashloc],-1,hashkey);
      }

      if (icount < MaxTries) hash[2*hashloc+1] = giX;
      break;
   case PRIME_JUMP:
      jump = 1+hashkey%hash_jump_prime;
      hashloc = (hashkey*AA+BB)%prime%hash_table_size;
      old_key = atomic_cmpxchg(&hash[2*hashloc],-1,hashkey);

      for (icount = 1; old_key != hashkey && old_key != -1 && icount < MaxTries; icount++){
         hashloc += (icount*jump);
         hashloc %= hash_table_size;
      
         old_key = atomic_cmpxchg(&hash[2*hashloc],-1,hashkey);
      }

      if (icount < MaxTries) hash[2*hashloc+1] = giX;
      break;
   }
}

int read_hash(
            const int   hash_method,
            const ulong hash_table_size,
            const ulong AA,
            const ulong BB,
            const ulong hashkey,
   __global const int   *hash)
{
   int hashval = -1;
   uint hashloc;
   int icount = 0;
   uint jump;

   switch (hash_method) {
   case PERFECT_HASH:
      return(hash[hashkey]);
      break;
   case LINEAR:
      for (hashloc = (hashkey*AA+BB)%prime%hash_table_size; hash[2*hashloc] != hashkey && hash[2*hashloc] != -1; hashloc++,hashloc %= hash_table_size);
      if (hash[2*hashloc] != -1) hashval = hash[2*hashloc+1];
      return(hashval);
      break;
   case QUADRATIC:
      for (hashloc = (hashkey*AA+BB)%prime%hash_table_size; hash[2*hashloc] != hashkey && hash[2*hashloc] != -1; hashloc+=(icount*icount),hashloc %= hash_table_size){
         icount++;
      }
      if (hash[2*hashloc] != -1) hashval = hash[2*hashloc+1];
      return(hashval);
      break;
   case PRIME_JUMP:
      jump = 1+hashkey%hash_jump_prime;
      for (hashloc = (hashkey*AA+BB)%prime%hash_table_size; hash[2*hashloc] != hashkey && hash[2*hashloc] != -1; hashloc+=(icount*jump),hashloc %= hash_table_size){
         icount++;
      }
      if (hash[2*hashloc] != -1) hashval = hash[2*hashloc+1];
      return(hashval);
      break;
   default:
      return(-1);
      break;
   }
}

