#ifndef _HASH_H
#define _HASH_H

#include "ezcl/ezcl.h"

enum choose_hash_method
{  METHOD_UNSET = 0,            //  use 0 for no method set
   PERFECT_HASH,                //  perfect hash 1
   LINEAR,                      //  linear hash 2
   QUADRATIC,                   //  quadratic hash 3
   PRIME_JUMP  };               //  prime_jump hash 4

typedef unsigned int uint;
typedef unsigned long ulong;

#ifdef __cplusplus
extern "C"
{
#endif

int *compact_hash_init(int ncells, uint isize, uint jsize, uint report_level);
void write_hash(uint ic, ulong hashkey, int *hash);
int read_hash(ulong hashkey, int *hash);
void compact_hash_delete(int *hash);

void write_hash_collision_report(void);
void read_hash_collision_report(void);
void final_hash_collision_report(void);

const char *get_hash_kernel_source_string(void);
void hash_lib_init(void);
void hash_lib_terminate(void);

cl_mem gpu_compact_hash_init(ulong ncells, int imaxsize, int jmaxsize, int gpu_hash_method, uint hash_report_level_in,
   ulong *gpu_hash_table_size, ulong *hashsize, cl_mem *dev_hash_header_in);
cl_mem gpu_get_hash_header(void);
void gpu_compact_hash_delete(cl_mem dev_hash, cl_mem dev_hash_header);
int read_dev_hash(int hash_method, ulong hashtablesize, ulong AA, ulong BB, ulong hashkey, int *hash);

#ifdef __cplusplus
}
#endif


#endif // _HASH_H

