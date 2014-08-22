#ifndef CRUX_H_
#define CRUX_H_

#include <stdio.h>

enum rollback_types{
   ROLLBACK_NONE,
   ROLLBACK_DISK,
   ROLLBACK_IN_MEMORY
};

class Crux
{
   int num_of_rollback_states;
   int rollback_type;
   int checkpoint_counter;

public:

   Crux(int rollback_type_in, int num_of_rollback_states_in, bool restart);
   ~Crux();

   void store_begin(size_t nsize, int ncycle);
   void store_ints(int *int_vals, size_t nelem);
   void store_longs(long long *long_vals, size_t nelem);
   void store_bools(bool *bool_vals, size_t nelem);
   void store_doubles(double *double_vals, size_t nelem);
   void store_int_array(int *int_array, size_t nelem);
   void store_float_array(float *float_array, size_t nelem);
   void store_double_array(double *double_array, size_t nelem);
   void store_end(void);

   void    restore_begin(char *restart_file, int rollback_counter);
   void    restore_ints(int *int_vals, size_t nelem);
   void    restore_bools(bool *bool_vals, size_t nelem);
   void    restore_longs(long long *long_vals, size_t nelem);
   void    restore_doubles(double *double_vals, size_t nelem);
   int    *restore_int_array(int *int_array, size_t nsize);
   float  *restore_float_array(float *float_array, size_t nsize);
   double *restore_double_array(double *double_array, size_t nsize);
   void    restore_end(void);

};
#endif // CRUX_H_
