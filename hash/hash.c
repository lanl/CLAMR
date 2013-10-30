#include <stdio.h>
#define __USE_XOPEN
#include <stdlib.h>
#include "hash.h"
#include "genmalloc/genmalloc.h"


static ulong AA;
static ulong BB;
ulong prime=4294967291;
uint hashtablesize;
uint hash_stride;
uint hash_ncells;
uint write_hash_collisions;
uint read_hash_collisions;
uint hash_report_level = 2;
uint hash_queries;
int hash_method = METHOD_UNSET;
uint hash_jump_prime = 41;
double hash_mult = 3.0;

float mem_opt_factor;

int   choose_hash_method = METHOD_UNSET;

int *compact_hash_init(int ncells, uint isize, uint jsize, uint report_level){
   hash_ncells = 0;
   write_hash_collisions = 0;
   read_hash_collisions = 0;
   hash_queries = 0;
   hash_report_level = report_level;
   hash_stride = isize;
   int *hash = NULL;

   if (choose_hash_method != METHOD_UNSET) hash_method = choose_hash_method;

   uint compact_hash_size = (uint)((double)ncells*hash_mult);
   uint perfect_hash_size = (uint)(isize*jsize);

   if (hash_method == METHOD_UNSET){
      float hash_mem_factor = 20.0;
      float hash_mem_ratio = (double)perfect_hash_size/(double)compact_hash_size;
      if (mem_opt_factor != 1.0) hash_mem_factor /= (mem_opt_factor*0.2); 
      hash_method = (hash_mem_ratio < hash_mem_factor) ? PERFECT_HASH : QUADRATIC;

      if (hash_report_level >= 2) printf("DEBUG hash_method %d hash_mem_ratio %f hash_mem_factor %f mem_opt_factor %f perfect_hash_size %u compact_hash_size %u\n",
         hash_method,hash_mem_ratio,hash_mem_factor,mem_opt_factor,perfect_hash_size,compact_hash_size);
   }

   int do_compact_hash = (hash_method == PERFECT_HASH) ? 0 : 1;

   if (hash_report_level >= 2) printf("DEBUG do_compact_hash %d hash_method %d perfect_hash_size %u compact_hash_size %u\n",
      do_compact_hash,hash_method,perfect_hash_size,compact_hash_size);

   if (do_compact_hash) {
      hashtablesize = compact_hash_size;
      AA = (ulong)(1.0+(double)(prime-1)*drand48());
      BB = (ulong)(0.0+(double)(prime-1)*drand48());
      if (AA > prime-1 || BB > prime-1) exit(0);
      if (hash_report_level > 1) printf("Factors AA %lu BB %lu\n",AA,BB);

      hash = (int *)genvector(2*hashtablesize,sizeof(int));
      for (uint ii = 0; ii<2*hashtablesize; ii+=2){
         hash[ii] = -1;
      }
   } else {
      hashtablesize = perfect_hash_size;

      hash = (int *)genvector(hashtablesize,sizeof(int));
      for (uint ii = 0; ii<hashtablesize; ii++){
         hash[ii] = -1;
      }
   }

   if (hash_report_level >= 2) {
      printf("Hash table size %u perfect hash table size %u memory savings %u by percentage %lf\n",
        hashtablesize,isize*jsize,isize*jsize-hashtablesize,
        (double)hashtablesize/(double)(isize*jsize));
   }

   return(hash);
}

void write_hash(uint ic, ulong hashkey, int *hash){
   int icount = 0;
   uint hashloc;
   if (hash_method == PERFECT_HASH) {
      hash[hashkey] = ic;
      return;
   }
   if (hash_method == LINEAR){
      if (hash_report_level == 0) {
         for (hashloc = (hashkey*AA+BB)%prime%hashtablesize; hash[2*hashloc] != -1 && hash[2*hashloc]!= (int)hashkey; hashloc++,hashloc = hashloc%hashtablesize);
      } else if (hash_report_level == 1) {
         hash_ncells++;
         for (hashloc = (hashkey*AA+BB)%prime%hashtablesize; hash[2*hashloc] != -1 && hash[2*hashloc]!= (int)hashkey; hashloc++,hashloc = hashloc%hashtablesize){
            write_hash_collisions++;
         }
      } else if (hash_report_level == 2) {
         hash_ncells++;
         for (hashloc = (hashkey*AA+BB)%prime%hashtablesize; hash[2*hashloc] != -1 && hash[2*hashloc]!= (int)hashkey; hashloc++,hashloc = hashloc%hashtablesize){
            write_hash_collisions++;
         }
      } else if (hash_report_level == 3) {
         hash_ncells++;
         hashloc = (hashkey*AA+BB)%prime%hashtablesize;
         printf("%d: cell %d hashloc is %d hash[2*hashloc] = %d hashkey %lu ii %lu jj %lu\n",icount,ic,hashloc,hash[2*hashloc],hashkey,hashkey%hash_stride,hashkey/hash_stride);
         for (hashloc = (hashkey*AA+BB)%prime%hashtablesize; hash[2*hashloc] != -1 && hash[2*hashloc]!= (int)hashkey; hashloc++,hashloc = hashloc%hashtablesize){
            int hashloctmp = hashloc+1;
            hashloctmp = hashloctmp%hashtablesize;
            printf("%d: cell %d hashloc is %d hash[2*hashloc] = %d hashkey %lu ii %lu jj %lu\n",icount,ic,hashloctmp,hash[2*hashloctmp],hashkey,hashkey%hash_stride,hashkey/hash_stride);
            icount++;
         }
         write_hash_collisions += icount;
      }
   } else if (hash_method == QUADRATIC){
      if (hash_report_level == 0) {
         for (hashloc = (hashkey*AA+BB)%prime%hashtablesize; hash[2*hashloc] != -1 && hash[2*hashloc]!= (int)hashkey; hashloc+=(icount*icount),hashloc = hashloc%hashtablesize) {
            icount++;
         }
      } else if (hash_report_level == 1) {
         hash_ncells++;
         for (hashloc = (hashkey*AA+BB)%prime%hashtablesize; hash[2*hashloc] != -1 && hash[2*hashloc]!= (int)hashkey; hashloc+=(icount*icount),hashloc = hashloc%hashtablesize){
            icount++;
         }
         write_hash_collisions += icount;
      } else if (hash_report_level == 2) {
         hash_ncells++;
         for (hashloc = (hashkey*AA+BB)%prime%hashtablesize; hash[2*hashloc] != -1 && hash[2*hashloc]!= (int)hashkey; hashloc+=(icount*icount),hashloc = hashloc%hashtablesize){
            icount++;
         }
         write_hash_collisions += icount;
      } else if (hash_report_level == 3) {
         hash_ncells++;
         hashloc = (hashkey*AA+BB)%prime%hashtablesize;
         printf("%d: cell %d hashloc is %d hash[2*hashloc] = %d hashkey %lu ii %lu jj %lu\n",icount,ic,hashloc,hash[2*hashloc],hashkey,hashkey%hash_stride,hashkey/hash_stride);
         for (hashloc = (hashkey*AA+BB)%prime%hashtablesize; hash[2*hashloc] != -1 && hash[2*hashloc]!= (int)hashkey; hashloc+=(icount*icount),hashloc = hashloc%hashtablesize){
            icount++;
            int hashloctmp = hashloc+icount*icount;
            hashloctmp = hashloctmp%hashtablesize;
            printf("%d: cell %d hashloc is %d hash[2*hashloc] = %d hashkey %lu ii %lu jj %lu\n",icount,ic,hashloctmp,hash[2*hashloctmp],hashkey,hashkey%hash_stride,hashkey/hash_stride);
         }
         write_hash_collisions += icount;
      }
   } else if (hash_method == PRIME_JUMP){
      uint jump = 1+hashkey%hash_jump_prime;
      if (hash_report_level == 0) {
         for (hashloc = (hashkey*AA+BB)%prime%hashtablesize; hash[2*hashloc] != -1 && hash[2*hashloc]!= (int)hashkey; hashloc+=(icount*jump),hashloc = hashloc%hashtablesize) {
            icount++;
         }
      } else if (hash_report_level == 1) {
         hash_ncells++;
         for (hashloc = (hashkey*AA+BB)%prime%hashtablesize; hash[2*hashloc] != -1 && hash[2*hashloc]!= (int)hashkey; hashloc+=(icount*jump),hashloc = hashloc%hashtablesize){
            icount++;
         }
         write_hash_collisions += icount;
      } else if (hash_report_level == 2) {
         hash_ncells++;
         for (hashloc = (hashkey*AA+BB)%prime%hashtablesize; hash[2*hashloc] != -1 && hash[2*hashloc]!= (int)hashkey; hashloc+=(icount*jump),hashloc = hashloc%hashtablesize){
            icount++;
         }
         write_hash_collisions += icount;
      } else if (hash_report_level == 3) {
         hash_ncells++;
         hashloc = (hashkey*AA+BB)%prime%hashtablesize;
         printf("%d: cell %d hashloc is %d hash[2*hashloc] = %d hashkey %lu ii %lu jj %lu\n",icount,ic,hashloc,hash[2*hashloc],hashkey,hashkey%hash_stride,hashkey/hash_stride);
         for (hashloc = (hashkey*AA+BB)%prime%hashtablesize; hash[2*hashloc] != -1 && hash[2*hashloc]!= (int)hashkey; hashloc+=(icount*jump),hashloc = hashloc%hashtablesize){
            icount++;
            int hashloctmp = hashloc+1;
            hashloctmp = hashloctmp%hashtablesize;
            printf("%d: cell %d hashloc is %d hash[2*hashloc] = %d hashkey %lu ii %lu jj %lu\n",icount,ic,hashloctmp,hash[2*hashloctmp],hashkey,hashkey%hash_stride,hashkey/hash_stride);
         }
         write_hash_collisions += icount;
      }
   }

   hash[2*hashloc] = hashkey;
   hash[2*hashloc+1] = ic;
}

int read_hash(ulong hashkey, int *hash){
   int max_collisions_allowed = 1000;
   int hashval = -1;
   uint hashloc;
   int icount=0;
   if (hash_method == PERFECT_HASH) {
      return(hash[hashkey]);
   }
   if (hash_method == LINEAR) {
      if (hash_report_level == 0) {
         for (hashloc = (hashkey*AA+BB)%prime%hashtablesize; hash[2*hashloc] != (int)hashkey && hash[2*hashloc] != -1; hashloc++,hashloc = hashloc%hashtablesize){
            icount++;
         }
      } else if (hash_report_level == 1) {
         hash_queries++;
         for (hashloc = (hashkey*AA+BB)%prime%hashtablesize; hash[2*hashloc] != (int)hashkey && hash[2*hashloc] != -1; hashloc++,hashloc = hashloc%hashtablesize){
            icount++;
         }
         read_hash_collisions += icount;
      } else if (hash_report_level == 2) {
         hash_queries++;
         for (hashloc = (hashkey*AA+BB)%prime%hashtablesize; hash[2*hashloc] != (int)hashkey && hash[2*hashloc] != -1; hashloc++,hashloc = hashloc%hashtablesize){
            icount++;
            if (icount > max_collisions_allowed) {
               printf("Error -- too many read hash collisions\n");
               exit(0);
            }
         }
         read_hash_collisions += icount;
      } else if (hash_report_level == 3) {
         hash_queries++;
         hashloc = (hashkey*AA+BB)%prime%hashtablesize;
         printf("%d: hashloc is %d hash[2*hashloc] = %d hashkey %lu ii %lu jj %lu\n",icount,hashloc,hash[2*hashloc],hashkey,hashkey%hash_stride,hashkey/hash_stride);
         for (hashloc = (hashkey*AA+BB)%prime%hashtablesize; hash[2*hashloc] != (int)hashkey && hash[2*hashloc] != -1; hashloc++,hashloc = hashloc%hashtablesize){
            icount++;
            uint hashloctmp = hashloc+1;
            hashloctmp = hashloctmp%hashtablesize;
            printf("%d: hashloc is %d hash[2*hashloc] = %d hashkey %lu ii %lu jj %lu\n",icount,hashloctmp,hash[2*hashloctmp],hashkey,hashkey%hash_stride,hashkey/hash_stride);
            if (icount > max_collisions_allowed) {
               printf("Error -- too many read hash collisions\n");
               exit(0);
            }
         }
         read_hash_collisions += icount;
      }
   } else if (hash_method == QUADRATIC) {
      if (hash_report_level == 0) {
         for (hashloc = (hashkey*AA+BB)%prime%hashtablesize; hash[2*hashloc] != (int)hashkey && hash[2*hashloc] != -1; hashloc+=(icount*icount),hashloc = hashloc%hashtablesize){
            icount++;
         }
      } else if (hash_report_level == 1) {
         hash_queries++;
         for (hashloc = (hashkey*AA+BB)%prime%hashtablesize; hash[2*hashloc] != (int)hashkey && hash[2*hashloc] != -1; hashloc+=(icount*icount),hashloc = hashloc%hashtablesize){
            icount++;
         }
         read_hash_collisions += icount;
      } else if (hash_report_level == 2) {
         hash_queries++;
         for (hashloc = (hashkey*AA+BB)%prime%hashtablesize; hash[2*hashloc] != (int)hashkey && hash[2*hashloc] != -1; hashloc+=(icount*icount),hashloc = hashloc%hashtablesize){
            icount++;
            if (icount > max_collisions_allowed) {
               printf("Error -- too many read hash collisions\n");
               exit(0);
            }
         }
         read_hash_collisions += icount;
      } else if (hash_report_level == 3) {
         hash_queries++;
         hashloc = (hashkey*AA+BB)%prime%hashtablesize;
         printf("%d: hashloc is %d hash[2*hashloc] = %d hashkey %lu ii %lu jj %lu\n",icount,hashloc,hash[2*hashloc],hashkey,hashkey%hash_stride,hashkey/hash_stride);
         for (hashloc = (hashkey*AA+BB)%prime%hashtablesize; hash[2*hashloc] != (int)hashkey && hash[2*hashloc] != -1; hashloc+=(icount*icount),hashloc = hashloc%hashtablesize){
            icount++;
            uint hashloctmp = hashloc+1;
            hashloctmp = hashloctmp%hashtablesize;
            printf("%d: hashloc is %d hash[2*hashloc] = %d hashkey %lu ii %lu jj %lu\n",icount,hashloctmp,hash[2*hashloctmp],hashkey,hashkey%hash_stride,hashkey/hash_stride);
            if (icount > max_collisions_allowed) {
               printf("Error -- too many read hash collisions\n");
               exit(0);
            }
         }
         read_hash_collisions += icount;
      }
   } else if (hash_method == PRIME_JUMP) {
      uint jump = 1+hashkey%hash_jump_prime;
      if (hash_report_level == 0) {
         for (hashloc = (hashkey*AA+BB)%prime%hashtablesize; hash[2*hashloc] != (int)hashkey && hash[2*hashloc] != -1; hashloc+=(icount*jump),hashloc = hashloc%hashtablesize){
            icount++;
         }
      } else if (hash_report_level == 1) {
         hash_queries++;
         for (hashloc = (hashkey*AA+BB)%prime%hashtablesize; hash[2*hashloc] != (int)hashkey && hash[2*hashloc] != -1; hashloc+=(icount*jump),hashloc = hashloc%hashtablesize){
            icount++;
         }
         read_hash_collisions += icount;
      } else if (hash_report_level == 2) {
         hash_queries++;
         for (hashloc = (hashkey*AA+BB)%prime%hashtablesize; hash[2*hashloc] != (int)hashkey && hash[2*hashloc] != -1; hashloc+=(icount*jump),hashloc = hashloc%hashtablesize){
            icount++;
            if (icount > max_collisions_allowed) {
               printf("Error -- too many read hash collisions\n");
               exit(0);
            }
         }
         read_hash_collisions += icount;
      } else if (hash_report_level == 3) {
         hash_queries++;
         hashloc = (hashkey*AA+BB)%prime%hashtablesize;
         printf("%d: hashloc is %d hash[2*hashloc] = %d hashkey %lu ii %lu jj %lu\n",icount,hashloc,hash[2*hashloc],hashkey,hashkey%hash_stride,hashkey/hash_stride);
         for (hashloc = (hashkey*AA+BB)%prime%hashtablesize; hash[2*hashloc] != (int)hashkey && hash[2*hashloc] != -1; hashloc+=(icount*jump),hashloc = hashloc%hashtablesize){
            icount++;
            uint hashloctmp = hashloc+1;
            hashloctmp = hashloctmp%hashtablesize;
            printf("%d: hashloc is %d hash[2*hashloc] = %d hashkey %lu ii %lu jj %lu\n",icount,hashloctmp,hash[2*hashloctmp],hashkey,hashkey%hash_stride,hashkey/hash_stride);
            if (icount > max_collisions_allowed) {
               printf("Error -- too many read hash collisions\n");
               exit(0);
            }
         }
         read_hash_collisions += icount;
      }
   }

   if (hash[2*hashloc] != -1) hashval = hash[2*hashloc+1];
   return(hashval);
}

void compact_hash_delete(int *hash){
      genvectorfree((void *)hash);
}

