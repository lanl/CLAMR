#include <stdio.h>
#define __USE_XOPEN
#include <stdlib.h>
#include "hash.h"
#include "genmalloc/genmalloc.h"
#ifdef HAVE_OPENCL
#include "hashlib_kern.inc"
#include "hashlib_source_kern.inc"
#endif

static ulong AA;
static ulong BB;
static ulong prime=4294967291;
static uint hashtablesize;
static uint hash_stride;
static uint hash_ncells;
static uint write_hash_collisions;
static uint read_hash_collisions;
static double write_hash_collisions_runsum = 0.0;
static double read_hash_collisions_runsum = 0.0;
static uint write_hash_collisions_count = 0;
static uint read_hash_collisions_count = 0;
static uint hash_report_level = 2;
static uint hash_queries;
static int hash_method = METHOD_UNSET;
static uint hash_jump_prime = 41;
static double hash_mult = 3.0;

size_t hash_header_size = 16;

#ifdef HAVE_OPENCL
cl_mem dev_hash_header = NULL;
#endif

float mem_opt_factor;

int   choose_hash_method = METHOD_UNSET;

#define MIN(a,b) ((a) < (b) ? (a) : (b))

int (*read_hash)(ulong, int *);
void (*write_hash)(uint, ulong, int *);

int get_hash_method(void) {
  return(hash_method);
}

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

      if (hash_method == LINEAR){
         if (hash_report_level == 0){
            read_hash  = read_hash_linear;
            write_hash = write_hash_linear;
         } else if (hash_report_level == 1){
            read_hash  = read_hash_linear_report_level_1;
            write_hash = write_hash_linear_report_level_1;
         } else if (hash_report_level == 2){
            read_hash  = read_hash_linear_report_level_2;
            write_hash = write_hash_linear_report_level_2;
         } else if (hash_report_level == 3){
            read_hash  = read_hash_linear_report_level_3;
            write_hash = write_hash_linear_report_level_3;
         }
      } else if (hash_method == QUADRATIC) {
         if (hash_report_level == 0){
            read_hash  = read_hash_quadratic;
            write_hash = write_hash_quadratic;
         } else if (hash_report_level == 1){
            read_hash  = read_hash_quadratic_report_level_1;
            write_hash = write_hash_quadratic_report_level_1;
         } else if (hash_report_level == 2){
            read_hash  = read_hash_quadratic_report_level_2;
            write_hash = write_hash_quadratic_report_level_2;
         } else if (hash_report_level == 3){
            read_hash  = read_hash_quadratic_report_level_3;
            write_hash = write_hash_quadratic_report_level_3;
         }
      } else if (hash_method == PRIME_JUMP) {
         if (hash_report_level == 0){
            read_hash  = read_hash_primejump;
            write_hash = write_hash_primejump;
         } else if (hash_report_level == 1){
            read_hash  = read_hash_primejump_report_level_1;
            write_hash = write_hash_primejump_report_level_1;
         } else if (hash_report_level == 2){
            read_hash  = read_hash_primejump_report_level_2;
            write_hash = write_hash_primejump_report_level_2;
         } else if (hash_report_level == 3){
            read_hash  = read_hash_primejump_report_level_3;
            write_hash = write_hash_primejump_report_level_3;
         }
      }
   } else {
      hashtablesize = perfect_hash_size;

      hash = (int *)genvector(hashtablesize,sizeof(int));
#ifdef _OPENMP
#pragma omp parallel for
#endif
      for (uint ii = 0; ii<hashtablesize; ii++){
         hash[ii] = -1;
      }

      read_hash  = read_hash_perfect;
      write_hash = write_hash_perfect;
   }

   if (hash_report_level >= 2) {
      printf("Hash table size %u perfect hash table size %u memory savings %u by percentage %lf\n",
        hashtablesize,isize*jsize,isize*jsize-hashtablesize,
        (double)hashtablesize/(double)(isize*jsize));
   }

   return(hash);
}

void write_hash_perfect(uint ic, ulong hashkey, int *hash){
   hash[hashkey] = ic;
}

void write_hash_linear(uint ic, ulong hashkey, int *hash){
   uint hashloc;

   for (hashloc = (hashkey*AA+BB)%prime%hashtablesize; hash[2*hashloc] != -1 && hash[2*hashloc]!= (int)hashkey; hashloc++,hashloc = hashloc%hashtablesize);

   hash[2*hashloc] = hashkey;
   hash[2*hashloc+1] = ic;
}

void write_hash_linear_report_level_1(uint ic, ulong hashkey, int *hash){
   uint hashloc;

   hash_ncells++;
   for (hashloc = (hashkey*AA+BB)%prime%hashtablesize; hash[2*hashloc] != -1 && hash[2*hashloc]!= (int)hashkey; hashloc++,hashloc = hashloc%hashtablesize){
      write_hash_collisions++;
   }

   hash[2*hashloc] = hashkey;
   hash[2*hashloc+1] = ic;
}

void write_hash_linear_report_level_2(uint ic, ulong hashkey, int *hash){
   uint hashloc;

   hash_ncells++;
   for (hashloc = (hashkey*AA+BB)%prime%hashtablesize; hash[2*hashloc] != -1 && hash[2*hashloc]!= (int)hashkey; hashloc++,hashloc = hashloc%hashtablesize){
      write_hash_collisions++;
   }

   hash[2*hashloc] = hashkey;
   hash[2*hashloc+1] = ic;
}

void write_hash_linear_report_level_3(uint ic, ulong hashkey, int *hash){
   int icount = 0;
   uint hashloc;

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

   hash[2*hashloc] = hashkey;
   hash[2*hashloc+1] = ic;
}

void write_hash_quadratic(uint ic, ulong hashkey, int *hash){
   int icount = 0;
   uint hashloc;

   for (hashloc = (hashkey*AA+BB)%prime%hashtablesize; hash[2*hashloc] != -1 && hash[2*hashloc]!= (int)hashkey; hashloc+=(icount*icount),hashloc = hashloc%hashtablesize) {
      icount++;
   }

   hash[2*hashloc] = hashkey;
   hash[2*hashloc+1] = ic;
}

void write_hash_quadratic_report_level_1(uint ic, ulong hashkey, int *hash){
   int icount = 0;
   uint hashloc;

   hash_ncells++;
   for (hashloc = (hashkey*AA+BB)%prime%hashtablesize; hash[2*hashloc] != -1 && hash[2*hashloc]!= (int)hashkey; hashloc+=(icount*icount),hashloc = hashloc%hashtablesize){
      icount++;
   }
   write_hash_collisions += icount;

   hash[2*hashloc] = hashkey;
   hash[2*hashloc+1] = ic;
}

void write_hash_quadratic_report_level_2(uint ic, ulong hashkey, int *hash){
   int icount = 0;
   uint hashloc;

   hash_ncells++;
   for (hashloc = (hashkey*AA+BB)%prime%hashtablesize; hash[2*hashloc] != -1 && hash[2*hashloc]!= (int)hashkey; hashloc+=(icount*icount),hashloc = hashloc%hashtablesize){
      icount++;
   }
   write_hash_collisions += icount;

   hash[2*hashloc] = hashkey;
   hash[2*hashloc+1] = ic;
}

void write_hash_quadratic_report_level_3(uint ic, ulong hashkey, int *hash){
   int icount = 0;
   uint hashloc;

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

   hash[2*hashloc] = hashkey;
   hash[2*hashloc+1] = ic;
}

void write_hash_primejump(uint ic, ulong hashkey, int *hash){
   int icount = 0;
   uint hashloc;

   uint jump = 1+hashkey%hash_jump_prime;
   for (hashloc = (hashkey*AA+BB)%prime%hashtablesize; hash[2*hashloc] != -1 && hash[2*hashloc]!= (int)hashkey; hashloc+=(icount*jump),hashloc = hashloc%hashtablesize) {
      icount++;
   }

   hash[2*hashloc] = hashkey;
   hash[2*hashloc+1] = ic;
}

void write_hash_primejump_report_level_1(uint ic, ulong hashkey, int *hash){
   int icount = 0;
   uint hashloc;

   uint jump = 1+hashkey%hash_jump_prime;
   hash_ncells++;
   for (hashloc = (hashkey*AA+BB)%prime%hashtablesize; hash[2*hashloc] != -1 && hash[2*hashloc]!= (int)hashkey; hashloc+=(icount*jump),hashloc = hashloc%hashtablesize){
      icount++;
   }
   write_hash_collisions += icount;

   hash[2*hashloc] = hashkey;
   hash[2*hashloc+1] = ic;
}

void write_hash_primejump_report_level_2(uint ic, ulong hashkey, int *hash){
   int icount = 0;
   uint hashloc;

   uint jump = 1+hashkey%hash_jump_prime;
   hash_ncells++;
   for (hashloc = (hashkey*AA+BB)%prime%hashtablesize; hash[2*hashloc] != -1 && hash[2*hashloc]!= (int)hashkey; hashloc+=(icount*jump),hashloc = hashloc%hashtablesize){
      icount++;
   }
   write_hash_collisions += icount;

   hash[2*hashloc] = hashkey;
   hash[2*hashloc+1] = ic;
}

void write_hash_primejump_report_level_3(uint ic, ulong hashkey, int *hash){
   int icount = 0;
   uint hashloc;

   uint jump = 1+hashkey%hash_jump_prime;
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

   hash[2*hashloc] = hashkey;
   hash[2*hashloc+1] = ic;
}

int read_hash_perfect(ulong hashkey, int *hash){
   return(hash[hashkey]);
}

int read_hash_linear(ulong hashkey, int *hash){
   int hashval = -1;
   uint hashloc;
   int icount=0;

   for (hashloc = (hashkey*AA+BB)%prime%hashtablesize; hash[2*hashloc] != (int)hashkey && hash[2*hashloc] != -1; hashloc++,hashloc = hashloc%hashtablesize){
      icount++;
   }

   if (hash[2*hashloc] != -1) hashval = hash[2*hashloc+1];
   return(hashval);
}

int read_hash_linear_report_level_1(ulong hashkey, int *hash){
   int hashval = -1;
   uint hashloc;
   int icount=0;

   hash_queries++;
   for (hashloc = (hashkey*AA+BB)%prime%hashtablesize; hash[2*hashloc] != (int)hashkey && hash[2*hashloc] != -1; hashloc++,hashloc = hashloc%hashtablesize){
      icount++;
   }
   read_hash_collisions += icount;

   if (hash[2*hashloc] != -1) hashval = hash[2*hashloc+1];
   return(hashval);
}

int read_hash_linear_report_level_2(ulong hashkey, int *hash){
   int max_collisions_allowed = 1000;
   int hashval = -1;
   uint hashloc;
   int icount=0;

   hash_queries++;
   for (hashloc = (hashkey*AA+BB)%prime%hashtablesize; hash[2*hashloc] != (int)hashkey && hash[2*hashloc] != -1; hashloc++,hashloc = hashloc%hashtablesize){
      icount++;
      if (icount > max_collisions_allowed) {
         printf("Error -- too many read hash collisions\n");
         exit(0);
      }
   }
   read_hash_collisions += icount;

   if (hash[2*hashloc] != -1) hashval = hash[2*hashloc+1];
   return(hashval);
}

int read_hash_linear_report_level_3(ulong hashkey, int *hash){
   int max_collisions_allowed = 1000;
   int hashval = -1;
   uint hashloc;
   int icount=0;

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

   if (hash[2*hashloc] != -1) hashval = hash[2*hashloc+1];
   return(hashval);
}

int read_hash_quadratic(ulong hashkey, int *hash){
   int hashval = -1;
   uint hashloc;
   int icount=0;

   for (hashloc = (hashkey*AA+BB)%prime%hashtablesize; hash[2*hashloc] != (int)hashkey && hash[2*hashloc] != -1; hashloc+=(icount*icount),hashloc = hashloc%hashtablesize){
      icount++;
   }

   if (hash[2*hashloc] != -1) hashval = hash[2*hashloc+1];
   return(hashval);
}

int read_hash_quadratic_report_level_1(ulong hashkey, int *hash){
   int hashval = -1;
   uint hashloc;
   int icount=0;

   hash_queries++;
   for (hashloc = (hashkey*AA+BB)%prime%hashtablesize; hash[2*hashloc] != (int)hashkey && hash[2*hashloc] != -1; hashloc+=(icount*icount),hashloc = hashloc%hashtablesize){
      icount++;
   }
   read_hash_collisions += icount;

   if (hash[2*hashloc] != -1) hashval = hash[2*hashloc+1];
   return(hashval);
}

int read_hash_quadratic_report_level_2(ulong hashkey, int *hash){
   int max_collisions_allowed = 1000;
   int hashval = -1;
   uint hashloc;
   int icount=0;

   hash_queries++;
   for (hashloc = (hashkey*AA+BB)%prime%hashtablesize; hash[2*hashloc] != (int)hashkey && hash[2*hashloc] != -1; hashloc+=(icount*icount),hashloc = hashloc%hashtablesize){
      icount++;
      if (icount > max_collisions_allowed) {
         printf("Error -- too many read hash collisions\n");
         exit(0);
      }
   }
   read_hash_collisions += icount;

   if (hash[2*hashloc] != -1) hashval = hash[2*hashloc+1];
   return(hashval);
}

int read_hash_quadratic_report_level_3(ulong hashkey, int *hash){
   int max_collisions_allowed = 1000;
   int hashval = -1;
   uint hashloc;
   int icount=0;

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

   if (hash[2*hashloc] != -1) hashval = hash[2*hashloc+1];
   return(hashval);
}

int read_hash_primejump(ulong hashkey, int *hash){
   int hashval = -1;
   uint hashloc;
   int icount=0;

   uint jump = 1+hashkey%hash_jump_prime;
   for (hashloc = (hashkey*AA+BB)%prime%hashtablesize; hash[2*hashloc] != (int)hashkey && hash[2*hashloc] != -1; hashloc+=(icount*jump),hashloc = hashloc%hashtablesize){
      icount++;
   }

   if (hash[2*hashloc] != -1) hashval = hash[2*hashloc+1];
   return(hashval);
}

int read_hash_primejump_report_level_1(ulong hashkey, int *hash){
   int hashval = -1;
   uint hashloc;
   int icount=0;

   uint jump = 1+hashkey%hash_jump_prime;
   hash_queries++;
   for (hashloc = (hashkey*AA+BB)%prime%hashtablesize; hash[2*hashloc] != (int)hashkey && hash[2*hashloc] != -1; hashloc+=(icount*jump),hashloc = hashloc%hashtablesize){
      icount++;
   }
   read_hash_collisions += icount;

   if (hash[2*hashloc] != -1) hashval = hash[2*hashloc+1];
   return(hashval);
}

int read_hash_primejump_report_level_2(ulong hashkey, int *hash){
   int max_collisions_allowed = 1000;
   int hashval = -1;
   uint hashloc;
   int icount=0;

   uint jump = 1+hashkey%hash_jump_prime;
   hash_queries++;
   for (hashloc = (hashkey*AA+BB)%prime%hashtablesize; hash[2*hashloc] != (int)hashkey && hash[2*hashloc] != -1; hashloc+=(icount*jump),hashloc = hashloc%hashtablesize){
      icount++;
      if (icount > max_collisions_allowed) {
         printf("Error -- too many read hash collisions\n");
         exit(0);
      }
   }
   read_hash_collisions += icount;

   if (hash[2*hashloc] != -1) hashval = hash[2*hashloc+1];
   return(hashval);
}

int read_hash_primejump_report_level_3(ulong hashkey, int *hash){
   int max_collisions_allowed = 1000;
   int hashval = -1;
   uint hashloc;
   int icount=0;

   uint jump = 1+hashkey%hash_jump_prime;
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

   if (hash[2*hashloc] != -1) hashval = hash[2*hashloc+1];
   return(hashval);
}

void compact_hash_delete(int *hash){
   read_hash = NULL;
   genvectorfree((void *)hash);
   hash_method = METHOD_UNSET;
}

void write_hash_collision_report(void){
   if (hash_method == PERFECT_HASH) return;
   if (hash_report_level == 1) {
      write_hash_collisions_runsum += (double)write_hash_collisions/(double)hash_ncells;
      write_hash_collisions_count++;
   } else if (hash_report_level >= 2) {
      printf("Write hash collision report -- collisions per cell %lf, collisions %d cells %d\n",(double)write_hash_collisions/(double)hash_ncells,write_hash_collisions,hash_ncells);
   }
}

void read_hash_collision_report(void){
   //printf("hash table size  bytes %ld\n",hashtablesize*sizeof(int));
   if (hash_method == PERFECT_HASH) return;
   if (hash_report_level == 1) {
      read_hash_collisions_runsum += (double)read_hash_collisions/(double)hash_queries;
      read_hash_collisions_count++;
   } else if (hash_report_level >= 2) {
      printf("Read hash collision report -- collisions per cell %lf, collisions %d cells %d\n",(double)read_hash_collisions/(double)hash_queries,read_hash_collisions,hash_queries);
      hash_queries = 0;
      read_hash_collisions = 0;
   }
}

void final_hash_collision_report(void){
   printf("hash table size  bytes %ld\n",hashtablesize*sizeof(int));
   if (hash_report_level >= 1 && read_hash_collisions_count > 0) { 
      printf("Final hash collision report -- write/read collisions per cell %lf/%lf\n",write_hash_collisions_runsum/(double)write_hash_collisions_count,read_hash_collisions_runsum/(double)read_hash_collisions_count);
   }
}

#ifdef HAVE_OPENCL
const char *get_hash_kernel_source_string(void)
{
   return(hashlib_source_kern_source);
}
#endif

#ifdef HAVE_OPENCL
static cl_kernel kernel_hash_init;
void hash_lib_init(void){

   cl_context context = ezcl_get_context();

   const char *defines = NULL;
   cl_program program = ezcl_create_program_wsource(context, defines, hashlib_kern_source);

   kernel_hash_init = ezcl_create_kernel_wprogram(program, "hash_init_cl");

   ezcl_program_release(program);
}

void hash_lib_terminate(void){
   ezcl_kernel_release(kernel_hash_init);
}

cl_mem gpu_compact_hash_init(ulong ncells, int imaxsize, int jmaxsize, int gpu_hash_method, uint hash_report_level_in,
   ulong *gpu_hash_table_size, ulong *hashsize, cl_mem *dev_hash_header_in)
{
   hash_report_level = hash_report_level_in;

   uint gpu_compact_hash_size = (uint)((double)ncells*hash_mult);
   uint gpu_perfect_hash_size = (uint)(imaxsize*jmaxsize);

   if (gpu_hash_method == METHOD_UNSET) {
      float gpu_hash_mem_factor = 20.0;
      float gpu_hash_mem_ratio = (double)gpu_perfect_hash_size/(double)gpu_compact_hash_size;
      if (mem_opt_factor != 1.0) gpu_hash_mem_factor /= (mem_opt_factor*0.2);
      gpu_hash_method = (gpu_hash_mem_ratio < gpu_hash_mem_factor) ? PERFECT_HASH : QUADRATIC;
   }

   int gpu_do_compact_hash = (gpu_hash_method == PERFECT_HASH) ? 0 : 1;

   ulong gpu_AA = 1;
   ulong gpu_BB = 0;
   if (gpu_do_compact_hash){
      (*gpu_hash_table_size) = gpu_compact_hash_size;
      gpu_AA = (ulong)(1.0+(double)(prime-1)*drand48());
      gpu_BB = (ulong)(0.0+(double)(prime-1)*drand48());
      //if ( gpu_AA > prime-1 || gpu_BB > prime-1) exit(0);
      (*hashsize) = 2*gpu_compact_hash_size;
   } else {
      (*gpu_hash_table_size) = gpu_perfect_hash_size;
      (*hashsize) = gpu_perfect_hash_size;
   }

   hashtablesize = (*hashsize);

   const uint TILE_SIZE = 128;

   cl_command_queue command_queue = ezcl_get_command_queue();

   cl_mem dev_hash = ezcl_malloc(NULL, "dev_hash", hashsize, sizeof(cl_int),  CL_MEM_READ_WRITE, 0);
   ulong *gpu_hash_header = (ulong *)genvector(hash_header_size, sizeof(ulong));
   gpu_hash_header[0] = (ulong)gpu_hash_method; 
   gpu_hash_header[1] =        (*gpu_hash_table_size);
   gpu_hash_header[2] =        gpu_AA;
   gpu_hash_header[3] =        gpu_BB;
   dev_hash_header = ezcl_malloc(NULL, "dev_hash_header", &hash_header_size, sizeof(cl_ulong),  CL_MEM_READ_WRITE, 0);
   ezcl_enqueue_write_buffer(command_queue, dev_hash_header, CL_TRUE, 0, hash_header_size*sizeof(cl_ulong), &gpu_hash_header[0], NULL);

   genvectorfree(gpu_hash_header);

   (*dev_hash_header_in) = dev_hash_header;

   size_t hash_local_work_size  = MIN((*hashsize), TILE_SIZE);
   size_t hash_global_work_size = (((*hashsize)+hash_local_work_size - 1) /hash_local_work_size) * hash_local_work_size;

   ezcl_set_kernel_arg(kernel_hash_init, 0, sizeof(cl_int),  (void *)hashsize);
   ezcl_set_kernel_arg(kernel_hash_init, 1, sizeof(cl_mem),  (void *)&dev_hash);
   ezcl_enqueue_ndrange_kernel(command_queue, kernel_hash_init,   1, NULL, &hash_global_work_size, &hash_local_work_size, NULL);

   return(dev_hash);
}

void gpu_compact_hash_delete(cl_mem dev_hash, cl_mem dev_hash_header){
   ezcl_device_memory_delete(dev_hash);
   ezcl_device_memory_delete(dev_hash_header);
   hash_method = METHOD_UNSET;
}

cl_mem gpu_get_hash_header(void){
   return(dev_hash_header);
}
#endif

int read_dev_hash(int hash_method, ulong hashtablesize, ulong AA, ulong BB, ulong hashkey, int *hash){
   //int hash_report_level = 3;
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
      } else {
         printf("Error -- Illegal value of hash_report_level %d\n",hash_report_level);
         exit(1);
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
      } else {
         printf("Error -- Illegal value of hash_report_level %d\n",hash_report_level);
         exit(1);
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
      } else {
         printf("Error -- Illegal value of hash_report_level %d\n",hash_report_level);
         exit(1);
      }
   } else {
      printf("Error -- Illegal value of hash_method %d\n",hash_method);
      exit(1);
   }

   if (hash[2*hashloc] != -1) hashval = hash[2*hashloc+1];
   return(hashval);
}

