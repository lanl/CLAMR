#include <float.h>
#include <stdlib.h>
#include <limits.h>
#include "l7.h"
#include "l7p.h"

#define L7_LOCATION "L7_REDUCE"

void L7_Sum(
      void                    *input,
      const int               count,
      const enum L7_Datatype  l7_datatype,
      void                    *output
      )
{
   
	double    local_sum_double;
	int       local_sum_int;
	long      local_sum_long;
	long long local_sum_long_long;
	int       i;
	
#ifdef HAVE_MPI
	double out_double;
	
	switch (l7_datatype){
	case L7_INT:
	case L7_INTEGER4:
		local_sum_int = 0;
		for (i=0; i<count; i++){
			local_sum_int += ((int*)input)[i];
		}
		if (l7.initialized_mpi){
			MPI_Allreduce(&local_sum_int, output, 1, MPI_INT, MPI_SUM,
					MPI_COMM_WORLD);
		}
		else{
			*((int *)output) = local_sum_int;
		}
		break;
   case L7_LONG:
      local_sum_long = 0L;
      for (i=0; i<count; i++){
         local_sum_long += ((long *)input)[i];
      }
      if (l7.initialized_mpi){
         MPI_Allreduce(&local_sum_long, output, 1, MPI_LONG,
               MPI_SUM, MPI_COMM_WORLD);
      }
      else{
         *((long *)output) = local_sum_long;
      }
      break;
	case L7_LONG_LONG_INT:
	case L7_INTEGER8:
		local_sum_long_long = 0L;
		for (i=0; i<count; i++){
			local_sum_long_long += ((long long*)input)[i];
		}
		if (l7.initialized_mpi){
			MPI_Allreduce(&local_sum_long_long, output, 1, MPI_LONG_LONG_INT,
					MPI_SUM, MPI_COMM_WORLD);
		}
		else{
			*((long long *)output) = local_sum_long_long;
		}
		break;
	case L7_FLOAT:
	case L7_REAL4:
		local_sum_double = 0.0;
		for (i=0; i<count; i++){
			local_sum_double += ((float*)input)[i];
		}
		if (l7.initialized_mpi){
			MPI_Allreduce(&local_sum_double, &out_double, 1, MPI_DOUBLE_PRECISION,
					MPI_SUM, MPI_COMM_WORLD);
			*((float*)output) = (float)out_double;
		}
		else{
			*((float*)output) = (float)local_sum_double;
		}
		break;
	case L7_DOUBLE:
	case L7_REAL8:
		local_sum_double = 0.0;
		for (i=0; i<count; i++){
			local_sum_double += ((double*)input)[i];
		}
		if (l7.initialized_mpi){
			MPI_Allreduce(&local_sum_double, output, 1, MPI_DOUBLE_PRECISION,
					MPI_SUM, MPI_COMM_WORLD);
		}
		else{
			*((double*)output) = local_sum_double;
		}
		break;
	default:
	   printf("Error -- L7_DATATYPE not supported in L7_Sum\n");
	   exit(1);
	   break;
	}
	return;
#else
   switch (l7_datatype){
   case L7_INT:
   case L7_INTEGER4:
      local_sum_int = 0;
      for (i=0; i<count; i++){
         local_sum_int += ((int*)input)[i];
      }
      *((int *)output) = local_sum_int;
      break;
   case L7_LONG:
      local_sum_long = 0L;
      for (i=0; i<count; i++){
         local_sum_long += ((long *)input)[i];
      }
      *((long *)output) = local_sum_long;
      break;
   case L7_LONG_LONG_INT:
   case L7_INTEGER8:
      local_sum_long_long = 0L;
      for (i=0; i<count; i++){
         local_sum_long_long += ((long long*)input)[i];
      }
      *((long long *)output) = local_sum_long_long;
      break;
   case L7_FLOAT:
   case L7_REAL4:
      local_sum_double = 0.0;
      for (i=0; i<count; i++){
         local_sum_double += ((float*)input)[i];
      }
      *((float*)output) = (float)local_sum_double;
      break;
   case L7_DOUBLE:
   case L7_REAL8:
      local_sum_double = 0.0;
      for (i=0; i<count; i++){
         local_sum_double += ((double*)input)[i];
      }
      *((double*)output) = local_sum_double;
      break;
   default:
      printf("Error -- L7_DATATYPE not supported in L7_Sum\n");
      exit(1);
      break;
   }
   return;
#endif
} /* End L7_Sum */

int L7_Max(
      void                    *input,
      const int               count,
      const enum L7_Datatype  l7_datatype,
      void                    *output
      )
{
	double    local_max_double;
	int       local_max_int;
   long      local_max_long;
	long long local_max_long_long;
	int       i;

#ifdef HAVE_MPI
	double out_double;
	
	switch (l7_datatype){
	case L7_INT:
	case L7_INTEGER4:
		local_max_int = INT_MIN;
		for (i=0; i<count; i++){
			if (((int *)input)[i] > local_max_int){
			   local_max_int = ((int*)input)[i];
			}
		}
		if (l7.initialized_mpi){
			MPI_Allreduce(&local_max_int, output, 1, MPI_INT, MPI_MAX,
					MPI_COMM_WORLD);
		}
		else{
			*((int *)output) = local_max_int;
		}
		break;
   case L7_LONG:
      local_max_long = LONG_MIN;
      for (i=0; i<count; i++){
         if (((long *)input)[i] > local_max_long){
            local_max_long = ((long *)input)[i];
         }
      }
      if (l7.initialized_mpi){
         MPI_Allreduce(&local_max_long, output, 1, MPI_LONG,
               MPI_MAX, MPI_COMM_WORLD);
      }
      else{
         *((long *)output) = local_max_long;
      }
      break;
	case L7_LONG_LONG_INT:
	case L7_INTEGER8:
		local_max_long_long = LONG_MIN;
		for (i=0; i<count; i++){
			if (((long long *)input)[i] > local_max_long_long){
			   local_max_long_long = ((long long *)input)[i];
			}
		}
		if (l7.initialized_mpi){
			MPI_Allreduce(&local_max_long_long, output, 1, MPI_LONG_LONG_INT,
					MPI_MAX, MPI_COMM_WORLD);
		}
		else{
			*((long long *)output) = local_max_long_long;
		}
		break;
	case L7_FLOAT:
	case L7_REAL4:
		local_max_double = -FLT_MAX;
		for (i=0; i<count; i++){
			if ((double)((float *)input)[i] > local_max_double){
			   local_max_double = (double)((float*)input)[i];
			}
		}
		if (l7.initialized_mpi){
			MPI_Allreduce(&local_max_double, &out_double, 1, MPI_DOUBLE_PRECISION,
					MPI_MAX, MPI_COMM_WORLD);
         *((float*)output) = (float)out_double;
		}
		else{
			*((float*)output) = (float)local_max_double;
		}
		break;
	case L7_DOUBLE:
	case L7_REAL8:
		local_max_double = -DBL_MAX;
		for (i=0; i<count; i++){
			if (((double *)input)[i] > local_max_double){
			   local_max_double = ((double*)input)[i];
			}
		}
		if (l7.initialized_mpi){
			MPI_Allreduce(&local_max_double, output, 1, MPI_DOUBLE_PRECISION,
					MPI_MAX, MPI_COMM_WORLD);
		}
		else{
			*((double*)output) = local_max_double;
		}
		break;
   default:
      printf("Error -- L7_DATATYPE not supported in L7_Max\n");
      exit(1);
      break;
	}
	return(0);
#else
   switch (l7_datatype){
   case L7_INT:
   case L7_INTEGER4:
      local_max_int = INT_MIN;
      for (i=0; i<count; i++){
         if (((int *)input)[i] > local_max_int){
            local_max_int = ((int*)input)[i];
         }
      }
      *((int *)output) = local_max_int;
      break;
   case L7_LONG:
      local_max_long = LONG_MIN;
      for (i=0; i<count; i++){
         if (((long *)input)[i] > local_max_long){
            local_max_long = ((long *)input)[i];
         }
      }
      *((long *)output) = local_max_long;
      break;
   case L7_LONG_LONG_INT:
   case L7_INTEGER8:
      local_max_long_long = LONG_MIN;
      for (i=0; i<count; i++){
         if (((long long *)input)[i] > local_max_long_long){
            local_max_long_long = ((long long *)input)[i];
         }
      }
      *((long long *)output) = local_max_long_long;
      break;
   case L7_FLOAT:
   case L7_REAL4:
      local_max_double = -FLT_MAX;
      for (i=0; i<count; i++){
         if ((double)((float *)input)[i] > local_max_double){
            local_max_double = (double)((float*)input)[i];
         }
      }
      *((float*)output) = (float)local_max_double;
      break;
   case L7_DOUBLE:
   case L7_REAL8:
      local_max_double = -DBL_MAX;
      for (i=0; i<count; i++){
         if (((double *)input)[i] > local_max_double){
            local_max_double = ((double*)input)[i];
         }
      }
      *((double*)output) = local_max_double;
      break;
   default:
      printf("Error -- L7_DATATYPE not supported in L7_Max\n");
      exit(1);
      break;
   }
   return(0);
#endif
} /* End L7_Max */

void l7_max_(
      void                    *input,
      const int               *count,
      const enum L7_Datatype  *l7_datatype,
      void                    *output,
      int                     *ierr
      )
{
	*ierr = L7_Max(input, *count, *l7_datatype, output);
} /* End l7_max_ */

int L7_Min(
      void                    *input,
      const int               count,
      const enum L7_Datatype  l7_datatype,
      void                    *output
      )
{
   double    local_min_double;
   int       local_min_int;
   long      local_min_long;
   long long local_min_long_long;
   int       i;
   
#ifdef HAVE_MPI
   double out_double;
   
   switch (l7_datatype){
   case L7_INT:
   case L7_INTEGER4:
      local_min_int = INT_MAX;
      for (i=0; i<count; i++){
         if (((int *)input)[i] < local_min_int){
            local_min_int = ((int*)input)[i];
         }
      }
      if (l7.initialized_mpi){
         MPI_Allreduce(&local_min_int, output, 1, MPI_INT, MPI_MIN,
               MPI_COMM_WORLD);
      }
      else{
         *((int *)output) = local_min_int;
      }
      break;
   case L7_LONG:
      local_min_long = LONG_MAX;
      for (i=0; i<count; i++){
         if (((long *)input)[i] < local_min_long){
            local_min_long = ((long *)input)[i];
         }
      }
      if (l7.initialized_mpi){
         MPI_Allreduce(&local_min_long, output, 1, MPI_LONG,
               MPI_MIN, MPI_COMM_WORLD);
      }
      else{
         *((long *)output) = local_min_long;
      }
      break;
   case L7_LONG_LONG_INT:
   case L7_INTEGER8:
      local_min_long_long = LONG_MAX;
      for (i=0; i<count; i++){
         if (((long long *)input)[i] < local_min_long_long){
            local_min_long_long = ((long long *)input)[i];
         }
      }
      if (l7.initialized_mpi){
         MPI_Allreduce(&local_min_long_long, output, 1, MPI_LONG_LONG_INT,
               MPI_MIN, MPI_COMM_WORLD);
      }
      else{
         *((long long *)output) = local_min_long_long;
      }
      break;
   case L7_FLOAT:
   case L7_REAL4:
      local_min_double = FLT_MAX;
      for (i=0; i<count; i++){
         if ((double)((float *)input)[i] < local_min_double){
            local_min_double = (double)((float*)input)[i];
         }
      }
      if (l7.initialized_mpi){
         MPI_Allreduce(&local_min_double, &out_double, 1, MPI_DOUBLE_PRECISION,
               MPI_MIN, MPI_COMM_WORLD);
         *((float*)output) = (float)out_double;
      }
      else{
         *((float*)output) = (float)local_min_double;
      }
      break;
   case L7_DOUBLE:
   case L7_REAL8:
      local_min_double = DBL_MAX;
      for (i=0; i<count; i++){
         if (((double *)input)[i] < local_min_double){
            local_min_double = ((double*)input)[i];
         }
      }
      if (l7.initialized_mpi){
         MPI_Allreduce(&local_min_double, output, 1, MPI_DOUBLE_PRECISION,
               MPI_MIN, MPI_COMM_WORLD);
      }
      else{
         *((double*)output) = local_min_double;
      }
      break;
   default:
      printf("Error -- L7_DATATYPE not supported in L7_Min\n");
      exit(1);
      break;
   }
   return(0);
#else
   switch (l7_datatype){
   case L7_INT:
   case L7_INTEGER4:
      local_min_int = INT_MAX;
      for (i=0; i<count; i++){
         if (((int *)input)[i] < local_min_int){
            local_min_int = ((int*)input)[i];
         }
      }
      *((int *)output) = local_min_int;
      break;
   case L7_LONG:
      local_min_long = LONG_MAX;
      for (i=0; i<count; i++){
         if (((long *)input)[i] < local_min_long){
            local_min_long = ((long *)input)[i];
         }
      }
      *((long *)output) = local_min_long;
      break;
   case L7_LONG_LONG_INT:
   case L7_INTEGER8:
      local_min_long_long = LONG_MAX;
      for (i=0; i<count; i++){
         if (((long long *)input)[i] < local_min_long_long){
            local_min_long_long = ((long long *)input)[i];
         }
      }
      *((long long *)output) = local_min_long_long;
      break;
   case L7_FLOAT:
   case L7_REAL4:
      local_min_double = FLT_MAX;
      for (i=0; i<count; i++){
         if ((double)((float *)input)[i] < local_min_double){
            local_min_double = (double)((float*)input)[i];
         }
      }
      *((float*)output) = (float)local_min_double;
      break;
   case L7_DOUBLE:
   case L7_REAL8:
      local_min_double = DBL_MAX;
      for (i=0; i<count; i++){
         if (((double *)input)[i] < local_min_double){
            local_min_double = ((double*)input)[i];
         }
      }
      *((double*)output) = local_min_double;
      break;
   default:
      printf("Error -- L7_DATATYPE not supported in L7_Min\n");
      exit(1);
      break;
   }
   return(0);
#endif
} /* End L7_Min */

void l7_min_(
      void                    *input,
      const int               *count,
      const enum L7_Datatype  *l7_datatype,
      void                    *output,
      int                     *ierr
      )
{
   *ierr = L7_Min(input, *count, *l7_datatype, output);
} /* End l7_min_ */

int L7_Any(
      void                    *input,
      const int               count,
      const enum L7_Datatype  l7_datatype,
      void                    *output
      )
{
	int       local_any_int;
	long      local_any_long;
	long long local_any_long_long;
	int       i;
	
#ifdef HAVE_MPI
	switch (l7_datatype){
	case L7_INT:
	case L7_LOGICAL:
		local_any_int = 0;
		for (i=0; i<count; i++){
			if (((int *)input)[i]){
			   local_any_int = ((int*)input)[i];
			}
		}
		if (l7.initialized_mpi){
			MPI_Allreduce(&local_any_int, output, 1, MPI_INT, MPI_SUM,
					MPI_COMM_WORLD);
			if (*((int*)output)){
				if (*((int *)output) < 0){
					*((int *)output) = -1;
				}
				else {
				    *((int *)output) = 1;	
				}
			}
		}
		else{
			*((int *)output) = local_any_int;
		}
		break;
   case L7_LONG:
      local_any_long = 0L;
      for (i=0; i<count; i++){
         if (((long *)input)[i]){
            local_any_long = ((long*)input)[i];
         }
      }
      if (l7.initialized_mpi){
         MPI_Allreduce(&local_any_long, output, 1, MPI_LONG, MPI_SUM,
               MPI_COMM_WORLD);
         if (*((long*)output)){
            if (*((long *)output) < 0){
               *((long *)output) = -1;
            }
            else {
                *((long *)output) = 1;  
            }
         }
      }
      else{
         *((long *)output) = local_any_long;
      }
      break;
	case L7_LONG_LONG_INT:
		local_any_long_long = 0L;
		for (i=0; i<count; i++){
			if (((long long *)input)[i]){
			   local_any_long_long = ((long long*)input)[i];
			}
		}
		if (l7.initialized_mpi){
			MPI_Allreduce(&local_any_long_long, output, 1, MPI_LONG_LONG_INT, MPI_SUM,
					MPI_COMM_WORLD);
			if (*((long long*)output)){
				if (*((long long *)output) < 0){
					*((long long *)output) = -1;
				}
				else {
				    *((long long *)output) = 1;	
				}
			}
		}
		else{
			*((long long *)output) = local_any_long_long;
		}
		break;
	 default:
	      printf("Error -- L7_DATATYPE not supported in L7_Any\n");
	      exit(1);
	      break;
	}
	return(0);
#else
   switch (l7_datatype){
   case L7_INT:
   case L7_LOGICAL:
      local_any_int = 0;
      for (i=0; i<count; i++){
         if (((int *)input)[i]){
            local_any_int = ((int*)input)[i];
         }
      }
      *((int *)output) = local_any_int;
      break;
   case L7_LONG:
      local_any_long = 0L;
      for (i=0; i<count; i++){
         if (((long *)input)[i]){
            local_any_long = ((long*)input)[i];
         }
      }
      *((long *)output) = local_any_long;
      break;
   case L7_LONG_LONG_INT:
      local_any_long_long = 0L;
      for (i=0; i<count; i++){
         if (((long long *)input)[i]){
            local_any_long_long = ((long long*)input)[i];
         }
      }
      *((long long *)output) = local_any_long_long;
      break;
    default:
      printf("Error -- L7_DATATYPE not supported in L7_Any\n");
      exit(1);
      break;
   }
   return(0);
#endif
} /* End L7_Any */

void l7_any_(
      void                   *input,
      const int              *count,
      const enum L7_Datatype *l7_datatype,
      void                   *output,
      int                    *ierr
      )
{
	*ierr = L7_Any(input, *count, *l7_datatype, output);
} /* End l7_any_ */

int L7_All(
      void                    *input,
      const int               count,
      const enum L7_Datatype  l7_datatype,
      void                    *output
      )
{
	int       local_all_int;
   long      local_all_long, sign_long;
	long long local_all_long_long, sign_long_long;
	int       i, sign;
	
#ifdef HAVE_MPI
	switch (l7_datatype){
	case L7_INT:
	case L7_LOGICAL:
		local_all_int = 0;
		sign = ((int *)input)[0];
		for (i=0; i<count; i++){
			if (! ((int *)input)[i]){
			   local_all_int = 0;
			}
		}
		if (l7.initialized_mpi){
			MPI_Allreduce(&local_all_int, output, 1, MPI_INT, MPI_LAND,
					MPI_COMM_WORLD);
			if (*((int *)output)){
				*((int *)output) = sign;
			}
		}
		else{
			*((int *)output) = local_all_int;
		}
		break;
   case L7_LONG:
      local_all_long = 0L;
      sign_long=((long *)input)[0];
      for (i=0; i<count; i++){
         if (! ((long *)input)[i]){
            local_all_long = 0L;
         }
      }
      if (l7.initialized_mpi){
         MPI_Allreduce(&local_all_long, output, 1, MPI_LONG, MPI_LAND,
               MPI_COMM_WORLD);
         if (*((long *)output)){
            *((long *)output) = sign_long;
         }
      }
      else{
         *((long *)output) = local_all_long;
      }
      break;
	case L7_LONG_LONG_INT:
		local_all_long_long = 0L;
		sign_long_long=((long long *)input)[0];
		for (i=0; i<count; i++){
			if (! ((long long *)input)[i]){
			   local_all_long_long = 0L;
			}
		}
		if (l7.initialized_mpi){
			MPI_Allreduce(&local_all_long_long, output, 1, MPI_LONG_LONG_INT, MPI_LAND,
					MPI_COMM_WORLD);
			if (*((long long *)output)){
				*((long long *)output) = sign_long_long;
			}
		}
		else{
			*((long long *)output) = local_all_long_long;
		}
		break;
   default:
        printf("Error -- L7_DATATYPE not supported in L7_All\n");
        exit(1);
        break;
	}
	return(0);
#else
   switch (l7_datatype){
   case L7_INT:
   case L7_LOGICAL:
      local_all_int = 0;
      sign = ((int *)input)[0];
      for (i=0; i<count; i++){
         if (! ((int *)input)[i]){
            local_all_int = 0;
         }
      }
      *((int *)output) = local_all_int;
      break;
   case L7_LONG:
      local_all_long = 0L;
      sign_long=((long *)input)[0];
      for (i=0; i<count; i++){
         if (! ((long *)input)[i]){
            local_all_long = 0L;
         }
      }
      *((long *)output) = local_all_long;
      break;
   case L7_LONG_LONG_INT:
      local_all_long_long = 0L;
      sign_long_long=((long long *)input)[0];
      for (i=0; i<count; i++){
         if (! ((long long *)input)[i]){
            local_all_long_long = 0L;
         }
      }
      *((long long *)output) = local_all_long_long;
      break;
   default:
        printf("Error -- L7_DATATYPE not supported in L7_All\n");
        exit(1);
        break;
   }
   return(0);
#endif
} /* End L7_All */

void l7_all_(
      void                    *input,
      const int               *count,
      const enum L7_Datatype  *l7_datatype,
      void                    *output,
      int                     *ierr
      )
{
	*ierr = L7_All(input, *count, *l7_datatype, output);
} /* End l7_all_ */

int L7_Array_Sum(
      void                    *input,
      const int               count,
      const enum L7_Datatype  l7_datatype,
      void                    *output
      )
{
#ifdef HAVE_MPI
	MPI_Datatype
	  mpi_type;
	
	switch (l7_datatype){
	case L7_INT:
	case L7_INTEGER4:
	case L7_LONG_LONG_INT:
	case L7_INTEGER8:
	case L7_FLOAT:
	case L7_REAL4:
	case L7_DOUBLE:
	case L7_REAL8:
		if (l7.initialized_mpi){
			mpi_type = l7p_mpi_type (l7_datatype);
			MPI_Allreduce(input, output, count, mpi_type, MPI_SUM,
					MPI_COMM_WORLD);
		}
		break;
   default:
        printf("Error -- L7_DATATYPE not supported in L7_Array_Sum\n");
        exit(1);
        break;
	}
#else
	int i;
	
   switch (l7_datatype){
   case L7_INT:
   case L7_INTEGER4:
      for (i=0; i<count; i++){
         ((int*)output)[i]=((int*)input)[i];
      }
   case L7_LONG:
      for (i=0; i<count; i++){
         ((long*)output)[i]=((long*)input)[i];
      }
   case L7_LONG_LONG_INT:
   case L7_INTEGER8:
      for (i=0; i<count; i++){
         ((long long*)output)[i]=((long long*)input)[i];
      }
   case L7_FLOAT:
   case L7_REAL4:
      for (i=0; i<count; i++){
         ((float*)output)[i]=((float*)input)[i];
      }
   case L7_DOUBLE:
   case L7_REAL8:
      for (i=0; i<count; i++){
         ((double*)output)[i]=((double*)input)[i];
      }
      break;
   default:
        printf("Error -- L7_DATATYPE not supported in L7_Array_Sum\n");
        exit(1);
        break;
   }
#endif
   return(0);
} /* End L7_Array_Sum */

void l7_array_sum_(
      void                    *input,
      const int               *count,
      const enum L7_Datatype  *l7_datatype,
      void                    *output,
      int                     *ierr
      )
{
	*ierr = L7_Array_Sum(input, *count, *l7_datatype, output);
} /* End l7_array_sum_ */

int L7_Array_Max(
      void                    *input,
      const int               count,
      const enum L7_Datatype  l7_datatype,
      void                    *output
      )
{
#ifdef HAVE_MPI
	MPI_Datatype
	  mpi_type;
	
	switch (l7_datatype){
	case L7_INT:
	case L7_INTEGER4:
	case L7_LONG_LONG_INT:
	case L7_INTEGER8:
	case L7_FLOAT:
	case L7_REAL4:
	case L7_DOUBLE:
	case L7_REAL8:
		if (l7.initialized_mpi){
			mpi_type = l7p_mpi_type (l7_datatype);
			MPI_Allreduce(input, output, count, mpi_type, MPI_MAX,
					MPI_COMM_WORLD);
		}
		break;
   default:
        printf("Error -- L7_DATATYPE not supported in L7_Array_Max\n");
        exit(1);
        break;
	}
#else
   int i;
   
   switch (l7_datatype){
   case L7_INT:
   case L7_INTEGER4:
      for (i=0; i<count; i++){
         ((int*)output)[i]=((int*)input)[i];
      }
   case L7_LONG:
      for (i=0; i<count; i++){
         ((long*)output)[i]=((long*)input)[i];
      }
   case L7_LONG_LONG_INT:
   case L7_INTEGER8:
      for (i=0; i<count; i++){
         ((long long*)output)[i]=((long long*)input)[i];
      }
   case L7_FLOAT:
   case L7_REAL4:
      for (i=0; i<count; i++){
         ((float*)output)[i]=((float*)input)[i];
      }
   case L7_DOUBLE:
   case L7_REAL8:
      for (i=0; i<count; i++){
         ((double*)output)[i]=((double*)input)[i];
      }
      break;
   default:
        printf("Error -- L7_DATATYPE not supported in L7_Array_Sum\n");
        exit(1);
        break;
   }
#endif
	return(0);
} /* End L7_Array_Max */

void l7_array_max_(
      void                    *input,
      const int               *count,
      const enum L7_Datatype  *l7_datatype,
      void                    *output,
      int                     *ierr
      )
{
	*ierr = L7_Array_Max(input, *count, *l7_datatype, output);
} /* End l7_array_max_ */

int L7_Array_Min(
      void                    *input,
      const int               count,
      const enum L7_Datatype  l7_datatype,
      void                    *output
      )
{
#ifdef HAVE_MPI
	MPI_Datatype
	  mpi_type;
	
	switch (l7_datatype){
	case L7_INT:
	case L7_INTEGER4:
	case L7_LONG_LONG_INT:
	case L7_INTEGER8:
	case L7_FLOAT:
	case L7_REAL4:
	case L7_DOUBLE:
	case L7_REAL8:
		if (l7.initialized_mpi){
			mpi_type = l7p_mpi_type (l7_datatype);
			MPI_Allreduce(input, output, count, mpi_type, MPI_MIN,
					MPI_COMM_WORLD);
		}
		break;
   default:
        printf("Error -- L7_DATATYPE not supported in L7_Array_Min\n");
        exit(1);
        break;
	}
#else
   int i;
   
   switch (l7_datatype){
   case L7_INT:
   case L7_INTEGER4:
      for (i=0; i<count; i++){
         ((int*)output)[i]=((int*)input)[i];
      }
   case L7_LONG:
      for (i=0; i<count; i++){
         ((long*)output)[i]=((long*)input)[i];
      }
   case L7_LONG_LONG_INT:
   case L7_INTEGER8:
      for (i=0; i<count; i++){
         ((long long*)output)[i]=((long long*)input)[i];
      }
   case L7_FLOAT:
   case L7_REAL4:
      for (i=0; i<count; i++){
         ((float*)output)[i]=((float*)input)[i];
      }
   case L7_DOUBLE:
   case L7_REAL8:
      for (i=0; i<count; i++){
         ((double*)output)[i]=((double*)input)[i];
      }
      break;
   default:
        printf("Error -- L7_DATATYPE not supported in L7_Array_Sum\n");
        exit(1);
        break;
   }
#endif
   return(0);
} /* End L7_Array_Min */

void l7_array_min_(
      void                    *input,
      const int               *count,
      const enum L7_Datatype  *l7_datatype,
      void                    *output,
      int                     *ierr
      )
{
	*ierr = L7_Array_Min(input, *count, *l7_datatype, output);
} /* End l7_array_min_ */

int L7_MaxLoc(
      void                    *input,
      const int               count,
      const enum L7_Datatype  l7_datatype,
      int                     *output
      )
{
   int       int_cur_max;
   long      long_cur_max;
   //long long longlong_cur_max;
   float     float_cur_max;
   double real_cur_max;
   
   int    i, local_maxloc;

   int *iptrinp = NULL;
   long *lptrinp = NULL;
   //long long *llptrinp=NULL;
   float *fptrinp = NULL;
   double *rptrinp = NULL;

#ifdef HAVE_MPI
	int    iprocs, istart;
	int    nprocs, mype;
	int    local_count[1];
	int    *counts = NULL;
	
   struct {
      int    value;
      int    index;
   } in_int, out_int;
   struct {
      long    value;
      int     index;
   } in_long, out_long;
//   struct {
//      long long   value;
//      int         index;
//   } in_longlong, out_longlong;
   struct {
      float    value;
      int      index;
   } in_float, out_float;
	struct {
		double value;
		int    index;
	} in_double, out_double;
	
	nprocs = l7.numpes;
	mype = l7.penum;
	
   istart = 0;
   if (l7.initialized_mpi){
     counts = (int *)malloc(nprocs*sizeof(int));
     local_count[0] = count;
	  MPI_Allgather(local_count, 1, MPI_INT, counts, 1, MPI_INT,
			MPI_COMM_WORLD);
	  for (iprocs = 0; iprocs < mype; iprocs++){
	    istart+=counts[iprocs];
	  }
	  free(counts);
   }
	
	switch (l7_datatype){
	   case L7_INT:
	   case L7_INTEGER4:
         int_cur_max=INT_MIN;
         local_maxloc = 0;
         iptrinp=(int *)input;
         for (i=0;i<count;i++){
            if (iptrinp[i] > int_cur_max){
               local_maxloc=i;
               int_cur_max=iptrinp[i];
            }
         }
         local_maxloc += istart;
         if (l7.initialized_mpi){
           in_int.value = int_cur_max;
           in_int.index = local_maxloc;
           MPI_Allreduce(&in_int, &out_int, 1, MPI_2INT,
               MPI_MAXLOC, MPI_COMM_WORLD);
           *output = out_int.index;
         }
         else {
            *output = local_maxloc;
         }
         break;
	   case L7_LONG:
         long_cur_max=LONG_MIN;
         local_maxloc = 0;
         lptrinp=(long *)input;
         for (i=0;i<count;i++){
            if (lptrinp[i] > long_cur_max){
               local_maxloc=i;
               long_cur_max=lptrinp[i];
            }
         }
         local_maxloc += istart;
         if (l7.initialized_mpi){
           in_long.value = long_cur_max;
           in_long.index = local_maxloc;
           MPI_Allreduce(&in_long, &out_long, 1, MPI_LONG_INT,
               MPI_MAXLOC, MPI_COMM_WORLD);
           *output = out_long.index;
         }
         else {
            *output = local_maxloc;
         }
         break;
	   case L7_FLOAT:
	   case L7_REAL4:
         float_cur_max=-FLT_MAX;
         local_maxloc = 0;
         fptrinp=(float *)input;
         for (i=0;i<count;i++){
            if (fptrinp[i] > float_cur_max){
               local_maxloc=i;
               float_cur_max=fptrinp[i];
            }
         }
         local_maxloc += istart;
         if (l7.initialized_mpi){
           in_float.value = float_cur_max;
           in_float.index = local_maxloc;
           MPI_Allreduce(&in_float, &out_float, 1, MPI_FLOAT_INT,
               MPI_MAXLOC, MPI_COMM_WORLD);
           *output = out_float.index;
         }
         else {
            *output = local_maxloc;
         }
         break;
	   case L7_DOUBLE:
	   case L7_REAL8:
    		real_cur_max=-DBL_MAX;
	   	local_maxloc = 0;
		   rptrinp=(double *)input;
		   for (i=0;i<count;i++){
			   if (rptrinp[i] > real_cur_max){
				   local_maxloc=i;
				   real_cur_max=rptrinp[i];
   			}
	   	}
		   local_maxloc += istart;
         if (l7.initialized_mpi){
           in_double.value = real_cur_max;
	   	  in_double.index = local_maxloc;
		     MPI_Allreduce(&in_double, &out_double, 1, MPI_DOUBLE_INT,
			   	MPI_MAXLOC, MPI_COMM_WORLD);
	   	  *output = out_double.index;
         }
         else {
            *output = local_maxloc;
         }
		   break;
   	 default:
	      printf("Error -- L7_DATATYPE not supported in L7_MaxLoc\n");
	      exit(1);
	      break;
	}
#else
   switch (l7_datatype){
      case L7_INT:
      case L7_INTEGER4:
         int_cur_max=INT_MIN;
         local_maxloc = 0;
         iptrinp=(int *)input;
         for (i=0;i<count;i++){
            if (iptrinp[i] > int_cur_max){
               local_maxloc=i;
               int_cur_max=iptrinp[i];
            }
         }
         *output = local_maxloc;
         break;
      case L7_LONG:
         long_cur_max=LONG_MIN;
         local_maxloc = 0;
         lptrinp=(long *)input;
         for (i=0;i<count;i++){
            if (lptrinp[i] > long_cur_max){
               local_maxloc=i;
               long_cur_max=lptrinp[i];
            }
         }
         *output = local_maxloc;
         break;
      case L7_FLOAT:
      case L7_REAL4:
         float_cur_max=-FLT_MAX;
         local_maxloc = 0;
         fptrinp=(float *)input;
         for (i=0;i<count;i++){
            if (fptrinp[i] > float_cur_max){
               local_maxloc=i;
               float_cur_max=fptrinp[i];
            }
         }
         *output = local_maxloc;
         break;
      case L7_DOUBLE:
      case L7_REAL8:
         real_cur_max=-DBL_MAX;
         local_maxloc = 0;
         rptrinp=(double *)input;
         for (i=0;i<count;i++){
            if (rptrinp[i] > real_cur_max){
               local_maxloc=i;
               real_cur_max=rptrinp[i];
            }
         }
         *output = local_maxloc;
         break;
       default:
         printf("Error -- L7_DATATYPE not supported in L7_MaxLoc\n");
         exit(1);
         break;
   }
#endif
	return(0);
} /* End L7_MaxLoc */

void l7_maxloc_(
      void                    *input,
      const int               *count,
      const enum L7_Datatype  *l7_datatype,
      int                     *output,
      int                     *ierr
      )
{
	*ierr=L7_MaxLoc(input, *count, *l7_datatype, output);
} /* End l7_maxloc_ */

int L7_MinLoc(
      void                    *input,
      const int               count,
      const enum L7_Datatype  l7_datatype,
      int                     *output
      )
{
   int       int_cur_min;
   long      long_cur_min;
   //long long longlong_cur_min;
   float     float_cur_min;
   double real_cur_min;
   
   int    i, local_minloc;
   
   int *iptrinp = NULL;
   long *lptrinp = NULL;
   //long long *llptrinp=NULL;
   float *fptrinp = NULL;
   double *rptrinp = NULL;

#ifdef HAVE_MPI
   int    iprocs, istart;
   int    nprocs, mype;
   int    local_count[1];
   int    *counts = NULL;
   
   struct {
      int    value;
      int    index;
   } in_int, out_int;
   struct {
      long    value;
      int     index;
   } in_long, out_long;
//   struct {
//      long long   value;
//      int         index;
//   } in_longlong, out_longlong;
   struct {
      float    value;
      int      index;
   } in_float, out_float;
   struct {
      double value;
      int    index;
   } in_double, out_double;
   
   nprocs = l7.numpes;
   mype = l7.penum;
   
   istart = 0;
   if (l7.initialized_mpi){
     counts = (int *)malloc(nprocs*sizeof(int));
     local_count[0] = count;
     MPI_Allgather(local_count, 1, MPI_INT, counts, 1, MPI_INT,
         MPI_COMM_WORLD);
     for (iprocs = 0; iprocs < mype; iprocs++){
        istart+=counts[iprocs];
     }
     free(counts);
   }
   
   switch (l7_datatype){
      case L7_INT:
      case L7_INTEGER4:
         int_cur_min=INT_MAX;
         local_minloc = 0;
         iptrinp=(int *)input;
         for (i=0;i<count;i++){
            if (iptrinp[i] < int_cur_min){
               local_minloc=i;
               int_cur_min=iptrinp[i];
            }
         }
         local_minloc += istart;
         if (l7.initialized_mpi){
           in_int.value = int_cur_min;
           in_int.index = local_minloc;
           MPI_Allreduce(&in_int, &out_int, 1, MPI_2INT,
               MPI_MINLOC, MPI_COMM_WORLD);
           *output = out_int.index;
         }
         else {
            *output = local_minloc;
         }
         break;
      case L7_LONG:
         long_cur_min=LONG_MAX;
         local_minloc = 0;
         lptrinp=(long *)input;
         for (i=0;i<count;i++){
            if (lptrinp[i] < long_cur_min){
               local_minloc=i;
               long_cur_min=lptrinp[i];
            }
         }
         local_minloc += istart;
         if (l7.initialized_mpi){
           in_long.value = long_cur_min;
           in_long.index = local_minloc;
           MPI_Allreduce(&in_long, &out_long, 1, MPI_LONG_INT,
               MPI_MINLOC, MPI_COMM_WORLD);
           *output = out_long.index;
         }
         else {
            *output = local_minloc;
         }
         break;
      case L7_FLOAT:
      case L7_REAL4:
         float_cur_min=FLT_MAX;
         local_minloc = 0;
         fptrinp=(float *)input;
         for (i=0;i<count;i++){
            if (fptrinp[i] < float_cur_min){
               local_minloc=i;
               float_cur_min=fptrinp[i];
            }
         }
         local_minloc += istart;
         if (l7.initialized_mpi){
           in_float.value = float_cur_min;
           in_float.index = local_minloc;
           MPI_Allreduce(&in_float, &out_float, 1, MPI_FLOAT_INT,
               MPI_MINLOC, MPI_COMM_WORLD);
           *output = out_float.index;
         }
         else {
            *output = local_minloc;
         }
         break;
      case L7_DOUBLE:
      case L7_REAL8:
         real_cur_min=DBL_MAX;
         local_minloc = 0;
         rptrinp=(double *)input;
         for (i=0;i<count;i++){
            if (rptrinp[i] < real_cur_min){
               local_minloc=i;
               real_cur_min=rptrinp[i];
            }
         }
         local_minloc += istart;
         if (l7.initialized_mpi){
           in_double.value = real_cur_min;
           in_double.index = local_minloc;
           MPI_Allreduce(&in_double, &out_double, 1, MPI_DOUBLE_INT,
               MPI_MINLOC, MPI_COMM_WORLD);
           *output = out_double.index;
         }
         else {
            *output = local_minloc;
         }
         break;
       default:
         printf("Error -- L7_DATATYPE not supported in L7_MinLoc\n");
         exit(1);
         break;
   }
#else
   switch (l7_datatype){
      case L7_INT:
      case L7_INTEGER4:
         int_cur_min=INT_MAX;
         local_minloc = 0;
         iptrinp=(int *)input;
         for (i=0;i<count;i++){
            if (iptrinp[i] < int_cur_min){
               local_minloc=i;
               int_cur_min=iptrinp[i];
            }
         }
         *output = local_minloc;
         break;
      case L7_LONG:
         long_cur_min=LONG_MAX;
         local_minloc = 0;
         lptrinp=(long *)input;
         for (i=0;i<count;i++){
            if (lptrinp[i] < long_cur_min){
               local_minloc=i;
               long_cur_min=lptrinp[i];
            }
         }
         *output = local_minloc;
         break;
      case L7_FLOAT:
      case L7_REAL4:
         float_cur_min=FLT_MAX;
         local_minloc = 0;
         fptrinp=(float *)input;
         for (i=0;i<count;i++){
            if (fptrinp[i] < float_cur_min){
               local_minloc=i;
               float_cur_min=fptrinp[i];
            }
         }
         *output = local_minloc;
         break;
      case L7_DOUBLE:
      case L7_REAL8:
         real_cur_min=DBL_MAX;
         local_minloc = 0;
         rptrinp=(double *)input;
         for (i=0;i<count;i++){
            if (rptrinp[i] < real_cur_min){
               local_minloc=i;
               real_cur_min=rptrinp[i];
            }
         }
         *output = local_minloc;
         break;
       default:
         printf("Error -- L7_DATATYPE not supported in L7_MinLoc\n");
         exit(1);
         break;
   }
#endif
   return(0);
} /* End L7_MinLoc */

void l7_minloc_(
      void                    *input,
      const int               *count,
      const enum L7_Datatype  *l7_datatype,
      int                     *output,
      int                     *ierr
      )
{
   *ierr=L7_MinLoc(input, *count, *l7_datatype, output);
} /* End l7_minloc_ */

int L7_MaxValLoc(
      void                    *input,
      const int               count,
      const enum L7_Datatype  l7_datatype,
      void                    *val,
      int                     *loc
      )
{
   int       int_cur_max;
   long      long_cur_max;
   //long long longlong_cur_max;
   float     float_cur_max;
   double real_cur_max;
   
   int    i, local_maxloc;

   int *iptrinp = NULL;
   long *lptrinp = NULL;
   //long long *llptrinp=NULL;
   float *fptrinp = NULL;
   double *rptrinp = NULL;
   
#ifdef HAVE_MPI
   
   int    iprocs, istart;
   int    nprocs, mype;
   int    local_count[1];
   int    *counts = NULL;
   
   struct {
      int    value;
      int    index;
   } in_int, out_int;
   struct {
      long    value;
      int     index;
   } in_long, out_long;
//   struct {
//      long long   value;
//      int         index;
//   } in_longlong, out_longlong;
   struct {
      float    value;
      int      index;
   } in_float, out_float;
   struct {
      double value;
      int    index;
   } in_double, out_double;
   
   nprocs = l7.numpes;
   mype = l7.penum;
   
   istart = 0;
   if (l7.initialized_mpi){
     counts = (int *)malloc(nprocs*sizeof(int));
     local_count[0] = count;
     MPI_Allgather(local_count, 1, MPI_INT, counts, 1, MPI_INT,
         MPI_COMM_WORLD);
     for (iprocs = 0; iprocs < mype; iprocs++){
        istart+=counts[iprocs];
     }
     free(counts);
   }
   
   switch (l7_datatype){
      case L7_INT:
      case L7_INTEGER4:
         int_cur_max=INT_MIN;
         local_maxloc = 0;
         iptrinp=(int *)input;
         for (i=0;i<count;i++){
            if (iptrinp[i] > int_cur_max){
               local_maxloc=i;
               int_cur_max=iptrinp[i];
            }
         }
         if (l7.initialized_mpi){
           local_maxloc += istart;
           in_int.value = int_cur_max;
           in_int.index = local_maxloc;
           MPI_Allreduce(&in_int, &out_int, 1, MPI_2INT,
               MPI_MAXLOC, MPI_COMM_WORLD);
           *((int *)val) = out_int.value;
           *loc = out_int.index;
         }
         else {
            *((int *)val) = int_cur_max;
            *loc = local_maxloc;
         }
         break;
      case L7_LONG:
         long_cur_max=LONG_MIN;
         local_maxloc = 0;
         lptrinp=(long *)input;
         for (i=0;i<count;i++){
            if (lptrinp[i] > long_cur_max){
               local_maxloc=i;
               long_cur_max=lptrinp[i];
            }
         }
         local_maxloc += istart;
         if (l7.initialized_mpi){
           in_long.value = long_cur_max;
           in_long.index = local_maxloc;
           MPI_Allreduce(&in_long, &out_long, 1, MPI_LONG_INT,
               MPI_MAXLOC, MPI_COMM_WORLD);
           *((long *)val) = out_long.value;
           *loc = out_long.index;
         }
         else {
            *((long *)val) = long_cur_max;
            *loc = local_maxloc;
         }
         break;
      case L7_FLOAT:
      case L7_REAL4:
         float_cur_max=-FLT_MAX;
         local_maxloc = 0;
         fptrinp=(float *)input;
         for (i=0;i<count;i++){
            if (fptrinp[i] > float_cur_max){
               local_maxloc=i;
               float_cur_max=fptrinp[i];
            }
         }
         local_maxloc += istart;
         if (l7.initialized_mpi){
           in_float.value = float_cur_max;
           in_float.index = local_maxloc;
           MPI_Allreduce(&in_float, &out_float, 1, MPI_FLOAT_INT,
               MPI_MAXLOC, MPI_COMM_WORLD);
           *((float *)val) = out_float.value;
           *loc = out_float.index;
         }
         else {
            *((float *)val) = float_cur_max;
            *loc = local_maxloc;
            
         }
         break;
      case L7_DOUBLE:
      case L7_REAL8:
         real_cur_max=-DBL_MAX;
         local_maxloc = 0;
         rptrinp=(double *)input;
         for (i=0;i<count;i++){
            if (rptrinp[i] > real_cur_max){
               local_maxloc=i;
               real_cur_max=rptrinp[i];
            }
         }
         if (l7.initialized_mpi){
            local_maxloc += istart;
            in_double.value = real_cur_max;
            in_double.index = local_maxloc;
            MPI_Allreduce(&in_double, &out_double, 1, MPI_DOUBLE_INT,
               MPI_MAXLOC, MPI_COMM_WORLD);
         
            *((double *)val) = out_double.value;
            *loc = out_double.index;
         }
         else {
            *((double *)val) = real_cur_max;
            *loc = local_maxloc;
         }

         break;
       default:
         printf("Error -- L7_DATATYPE not supported in L7_MaxValLoc\n");
         exit(1);
         break;
   }
#else
   switch (l7_datatype){
      case L7_INT:
      case L7_INTEGER4:
         int_cur_max=INT_MIN;
         local_maxloc = 0;
         iptrinp=(int *)input;
         for (i=0;i<count;i++){
            if (iptrinp[i] > int_cur_max){
               local_maxloc=i;
               int_cur_max=iptrinp[i];
            }
         }
         *((int *)val) = int_cur_max;
         *loc = local_maxloc;
         break;
      case L7_LONG:
         long_cur_max=LONG_MIN;
         local_maxloc = 0;
         lptrinp=(long *)input;
         for (i=0;i<count;i++){
            if (lptrinp[i] > long_cur_max){
               local_maxloc=i;
               long_cur_max=lptrinp[i];
            }
         }
         *((long *)val) = long_cur_max;
         *loc = local_maxloc;
         break;
      case L7_FLOAT:
      case L7_REAL4:
         float_cur_max=-FLT_MAX;
         local_maxloc = 0;
         fptrinp=(float *)input;
         for (i=0;i<count;i++){
            if (fptrinp[i] > float_cur_max){
               local_maxloc=i;
               float_cur_max=fptrinp[i];
            }
         }
         *((float *)val) = float_cur_max;
         *loc = local_maxloc;
         break;
      case L7_DOUBLE:
      case L7_REAL8:
         real_cur_max=-DBL_MAX;
         local_maxloc = 0;
         rptrinp=(double *)input;
         for (i=0;i<count;i++){
            if (rptrinp[i] > real_cur_max){
               local_maxloc=i;
               real_cur_max=rptrinp[i];
            }
         }
         *((double *)val) = real_cur_max;
         *loc = local_maxloc;
         break;
       default:
         printf("Error -- L7_DATATYPE not supported in L7_MaxValLoc\n");
         exit(1);
         break;
   }
#endif
   return(0);
} /* End L7_MaxValLoc */


void l7_maxvalloc_(
      void                    *input,
      const int               *count,
      const enum L7_Datatype  *l7_datatype,
      void                    *val,
      int                     *loc,
      int                     *ierr
      )
{
	*ierr=L7_MaxValLoc(input, *count, *l7_datatype, val, loc);
} /* End l7_maxvalloc_ */

int L7_MinValLoc(
      void                    *input,
      const int               count,
      const enum L7_Datatype  l7_datatype,
      void                    *val,
      int                     *loc
      )
{
   int       int_cur_min;
   long      long_cur_min;
   //long long longlong_cur_min;
   float     float_cur_min;
   double real_cur_min;
   
   int    i, local_minloc;

   int *iptrinp = NULL;
   long *lptrinp = NULL;
   //long long *llptrinp=NULL;
   float *fptrinp = NULL;
   double *rptrinp = NULL;

#ifdef HAVE_MPI
   int    iprocs, istart;
   int    nprocs, mype;
   int    local_count[1];
   int    *counts = NULL;
   
   struct {
      int    value;
      int    index;
   } in_int, out_int;
   struct {
      long    value;
      int     index;
   } in_long, out_long;
//   struct {
//      long long   value;
//      int         index;
//   } in_longlong, out_longlong;
   struct {
      float    value;
      int      index;
   } in_float, out_float;
   struct {
      double value;
      int    index;
   } in_double, out_double;
   
   nprocs = l7.numpes;
   mype = l7.penum;
   
   if (l7.initialized_mpi){
     counts = (int *)malloc(nprocs*sizeof(int));
     local_count[0] = count;
     MPI_Allgather(local_count, 1, MPI_INT, counts, 1, MPI_INT,
         MPI_COMM_WORLD);
     istart = 0;
     for (iprocs = 0; iprocs < mype; iprocs++){
        istart+=counts[iprocs];
     }
     free(counts);
   }
   
   switch (l7_datatype){
      case L7_INT:
      case L7_INTEGER4:
         int_cur_min=INT_MAX;
         local_minloc = 0;
         iptrinp=(int *)input;
         for (i=0;i<count;i++){
            if (iptrinp[i] < int_cur_min){
               local_minloc=i;
               int_cur_min=iptrinp[i];
            }
         }
         if (l7.initialized_mpi){
           local_minloc += istart;
           in_int.value = int_cur_min;
           in_int.index = local_minloc;
           MPI_Allreduce(&in_int, &out_int, 1, MPI_2INT,
               MPI_MINLOC, MPI_COMM_WORLD);
           *((int *)val) = out_int.value;
           *loc = out_int.index;
         }
         else {
            *((int *)val) = int_cur_min;
            *loc = local_minloc;
         }
         break;
      case L7_LONG:
         long_cur_min=LONG_MAX;
         local_minloc = 0;
         lptrinp=(long *)input;
         for (i=0;i<count;i++){
            if (lptrinp[i] < long_cur_min){
               local_minloc=i;
               long_cur_min=lptrinp[i];
            }
         }
         if (l7.initialized_mpi){
           local_minloc += istart;
           in_long.value = long_cur_min;
           in_long.index = local_minloc;
           MPI_Allreduce(&in_long, &out_long, 1, MPI_LONG_INT,
               MPI_MINLOC, MPI_COMM_WORLD);
           *((long *)val) = out_long.value;
           *loc = out_long.index;
         }
         else {
            *((long *)val) = long_cur_min;
            *loc = local_minloc;
         }
         break;
      case L7_FLOAT:
      case L7_REAL4:
         float_cur_min=FLT_MAX;
         local_minloc = 0;
         fptrinp=(float *)input;
         for (i=0;i<count;i++){
            if (fptrinp[i] < float_cur_min){
               local_minloc=i;
               float_cur_min=fptrinp[i];
            }
         }
         if (l7.initialized_mpi){
           local_minloc += istart;
           in_float.value = float_cur_min;
           in_float.index = local_minloc;
           MPI_Allreduce(&in_float, &out_float, 1, MPI_FLOAT_INT,
               MPI_MINLOC, MPI_COMM_WORLD);
           *((float *)val) = out_float.value;
           *loc = out_float.index;
         }
         else {
            *((float *)val) = float_cur_min;
            *loc = local_minloc;
         }
         break;
      case L7_DOUBLE:
      case L7_REAL8:
         real_cur_min=DBL_MAX;
         local_minloc = 0;
         rptrinp=(double *)input;
         for (i=0;i<count;i++){
            if (rptrinp[i] < real_cur_min){
               local_minloc=i;
               real_cur_min=rptrinp[i];
            }
         }
         if (l7.initialized_mpi){
            local_minloc += istart;
            in_double.value = real_cur_min;
            in_double.index = local_minloc;
            MPI_Allreduce(&in_double, &out_double, 1, MPI_DOUBLE_INT,
               MPI_MINLOC, MPI_COMM_WORLD);
         
            *((double *)val) = out_double.value;
            *loc = out_double.index;
         }
         else {
            *((double *)val) = real_cur_min;
            *loc = local_minloc;
         }

         break;
       default:
         printf("Error -- L7_DATATYPE not supported in L7_MinValLoc\n");
         exit(1);
         break;
   }
#else
   switch (l7_datatype){
      case L7_INT:
      case L7_INTEGER4:
         int_cur_min=INT_MAX;
         local_minloc = 0;
         iptrinp=(int *)input;
         for (i=0;i<count;i++){
            if (iptrinp[i] < int_cur_min){
               local_minloc=i;
               int_cur_min=iptrinp[i];
            }
         }
         *((int *)val) = int_cur_min;
         *loc = local_minloc;
         break;
      case L7_LONG:
         long_cur_min=LONG_MAX;
         local_minloc = 0;
         lptrinp=(long *)input;
         for (i=0;i<count;i++){
            if (lptrinp[i] < long_cur_min){
               local_minloc=i;
               long_cur_min=lptrinp[i];
            }
         }
         *((long *)val) = long_cur_min;
         *loc = local_minloc;
         break;
      case L7_FLOAT:
      case L7_REAL4:
         float_cur_min=FLT_MAX;
         local_minloc = 0;
         fptrinp=(float *)input;
         for (i=0;i<count;i++){
            if (fptrinp[i] < float_cur_min){
               local_minloc=i;
               float_cur_min=fptrinp[i];
            }
         }
         *((float *)val) = float_cur_min;
         *loc = local_minloc;
         break;
      case L7_DOUBLE:
      case L7_REAL8:
         real_cur_min=DBL_MAX;
         local_minloc = 0;
         rptrinp=(double *)input;
         for (i=0;i<count;i++){
            if (rptrinp[i] < real_cur_min){
               local_minloc=i;
               real_cur_min=rptrinp[i];
            }
         }
         *((double *)val) = real_cur_min;
         *loc = local_minloc;
         break;
       default:
         printf("Error -- L7_DATATYPE not supported in L7_MinValLoc\n");
         exit(1);
         break;
   }
#endif
   return(0);
} /* End L7_MinValLoc */


void l7_minvalloc_(
      void                    *input,
      const int               *count,
      const enum L7_Datatype  *l7_datatype,
      void                    *val,
      int                     *loc,
      int                     *ierr
      )
{
   *ierr=L7_MinValLoc(input, *count, *l7_datatype, val, loc);
} /* End l7_minvalloc_ */

int L7_MaxValLocLoc_Int4(
      void                    *input,
      void                    *loc2input,
      const int               count,
      const enum L7_Datatype  l7_datatype,
      void                    *val,
      void                    *loc,
      void                    *loc2
      )
{
   double real_cur_max;
   int    i, local_maxloc, local_maxloc2;
   double *rptrinp = NULL;
   int    *iptrinp = NULL;

#ifdef HAVE_MPI
	int    iprocs, istart;
	int    nprocs, mype;
	int    local_count[1];
	int    *counts = NULL;
	
	struct {
		double value;
		int    index;
	} in_double, out_double;
	
	nprocs = l7.numpes;
	mype = l7.penum;
	istart = 0;
	
	if (l7.initialized_mpi){
	  counts = (int *)malloc(nprocs*sizeof(int));
	  local_count[0] = count;
	  MPI_Allgather(local_count, 1, MPI_INT, counts, 1, MPI_INT,
			MPI_COMM_WORLD);
	  istart = 0;
	  for (iprocs = 0; iprocs < mype; iprocs++){
		  istart+=counts[iprocs];
	  }
	  free(counts);
	}
	
	switch (l7_datatype){
	case L7_REAL8:
		real_cur_max=-DBL_MAX;
		local_maxloc = 0;
		rptrinp=(double *)input;
		iptrinp=(int *)loc2input;
		for (i=0;i<count;i++){
			if (rptrinp[i] > real_cur_max){
				local_maxloc=i;
				local_maxloc2=iptrinp[i];
				real_cur_max=rptrinp[i];
			}
		}
		if (l7.initialized_mpi){
		  local_maxloc += istart;
		  in_double.value = real_cur_max;
		  in_double.index = local_maxloc;
		  MPI_Allreduce(&in_double, &out_double, 1, MPI_DOUBLE_INT,
				MPI_MAXLOC, MPI_COMM_WORLD);
		  *((double *)val) = out_double.value;
		  *((int *)loc) = out_double.index;
		  
		  in_double.index = local_maxloc2;
		  MPI_Allreduce(&in_double, &out_double, 1, MPI_DOUBLE_INT,
				  MPI_MAXLOC, MPI_COMM_WORLD);
		  *((int*)loc2) = out_double.index;
		}
		else {
			  *((double *)val) = real_cur_max;
			  *((int *)loc) = local_maxloc;
			  *((int *)loc2) = local_maxloc2;
			
		}
		break;
   default:
        printf("Error -- L7_DATATYPE not supported in L7_MaxValLocLoc_Int4\n");
        exit(1);
        break;
	}
#else
   switch (l7_datatype){
   case L7_REAL8:
      real_cur_max=-DBL_MAX;
      local_maxloc = 0;
      rptrinp=(double *)input;
      iptrinp=(int *)loc2input;
      for (i=0;i<count;i++){
         if (rptrinp[i] > real_cur_max){
            local_maxloc=i;
            local_maxloc2=iptrinp[i];
            real_cur_max=rptrinp[i];
         }
      }
      *((double *)val) = real_cur_max;
      *((int *)loc) = local_maxloc;
      *((int *)loc2) = local_maxloc2;
      break;
   default:
        printf("Error -- L7_DATATYPE not supported in L7_MaxValLocLoc_Int4\n");
        exit(1);
        break;
   }
#endif
	return(0);
} /* End L7_MaxValLocLoc_Int4 */

void l7_maxvallocloc_int4_(
      void                    *input,
      void                    *loc2input,
      const int               *count,
      const enum L7_Datatype  *l7_datatype,
      void                    *val,
      void                    *loc,
      void                    *loc2,
      int                     *ierr
      )
{
	*ierr=L7_MaxValLocLoc_Int4(input, loc2input, *count, *l7_datatype, val, loc, loc2);
} /* End l7_maxvallocloc_int4_ */

int L7_MaxValLocLoc_Int8(
      void                    *input,
      void                    *loc2input,
      const int               count,
      const enum L7_Datatype  l7_datatype,
      void                    *val,
      void                    *loc,
      void                    *loc2
      )
{
	double    real_cur_max;
	int       i, local_maxloc, local_maxloc2;
   double    *rptrinp = NULL;
   long long *iptrinp = NULL;
	
#ifdef HAVE_MPI
	int       iprocs, istart;
	int       nprocs, mype;
	int       local_count[1];
	int       *counts = NULL;
	
	struct {
		double value;
		int    index;
	} in_double, out_double;
	
	nprocs = l7.numpes;
	mype = l7.penum;
	istart = 0;
	
	if (l7.initialized_mpi){
	  counts = (int *)malloc(nprocs*sizeof(int));
	  local_count[0] = count;
	  MPI_Allgather(local_count, 1, MPI_INT, counts, 1, MPI_INT,
			MPI_COMM_WORLD);
	  istart = 0;
	  for (iprocs = 0; iprocs < mype; iprocs++){
		  istart+=counts[iprocs];
	  }
	  free(counts);
	}
	
	switch (l7_datatype){
	case L7_REAL8:
		real_cur_max=-DBL_MAX;
		local_maxloc = 0;
		rptrinp=(double *)input;
		iptrinp=(long long *)loc2input;
		for (i=0;i<count;i++){
			if (rptrinp[i] > real_cur_max){
				local_maxloc=i;
				local_maxloc2=(int)iptrinp[i];
				real_cur_max=rptrinp[i];
			}
		}
		if (l7.initialized_mpi){
		  local_maxloc += istart;
		  in_double.value = real_cur_max;
		  in_double.index = local_maxloc;
		  MPI_Allreduce(&in_double, &out_double, 1, MPI_DOUBLE_INT,
				MPI_MAXLOC, MPI_COMM_WORLD);
		  *((double *)val) = out_double.value;
		  *((int *)loc) = out_double.index;
		  
		  in_double.index = local_maxloc2;
		  MPI_Allreduce(&in_double, &out_double, 1, MPI_DOUBLE_INT,
				  MPI_MAXLOC, MPI_COMM_WORLD);
		  *((int*)loc2) = out_double.index;
		}
		else {
			  *((double *)val) = real_cur_max;
			  *((int *)loc) = local_maxloc;
			  *((int *)loc2) = local_maxloc2;
			
		}
		break;
   default:
        printf("Error -- L7_DATATYPE not supported in L7_MaxValLocLoc_Int8\n");
        exit(1);
        break;
	}
#else
   switch (l7_datatype){
   case L7_REAL8:
      real_cur_max=-DBL_MAX;
      local_maxloc = 0;
      rptrinp=(double *)input;
      iptrinp=(long long *)loc2input;
      for (i=0;i<count;i++){
         if (rptrinp[i] > real_cur_max){
            local_maxloc=i;
            local_maxloc2=iptrinp[i];
            real_cur_max=rptrinp[i];
         }
      }
      *((double *)val) = real_cur_max;
      *((int *)loc) = local_maxloc;
      *((int *)loc2) = local_maxloc2;
      break;
   default:
        printf("Error -- L7_DATATYPE not supported in L7_MaxValLocLoc_Int8\n");
        exit(1);
        break;
   }
#endif
	return(0);
} /* End L7_MaxValLocLoc_Int8 */

void l7_maxvallocloc_int8_(
      void                    *input,
      void                    *loc2input,
      const int               *count,
      const enum L7_Datatype  *l7_datatype,
      void                    *val,
      void                    *loc,
      void                    *loc2,
      int                     *ierr
      )
{
	*ierr=L7_MaxValLocLoc_Int8(input, loc2input, *count, *l7_datatype, val, loc, loc2);
} /* End l7_maxvallocloc_int8_ */

int L7_GetGlobal(
      void                    *array,
      const int               count,
      const enum L7_Datatype  l7_datatype,
      const int               index,
      void                    *val
      )
{
#ifdef HAVE_MPI
	int   nprocs, mype;
	int   iprocs, istart, index_start, index_end, root_pe;
	int   found_value;
	
	int   local_count[1];
	int   *counts = NULL;
	int   *iptrinp = NULL;
	
	nprocs = l7.numpes;
	mype   = l7.penum;
	
	istart = 0;
	index_start = 0;
	index_end   = -1;
	root_pe     = 0;
	found_value = 0;
	
	iptrinp = (int *)array;
	
	if (l7.initialized_mpi){
		  counts = (int *)malloc(nprocs*sizeof(int));
		  local_count[0] = count;
		  MPI_Allgather(local_count, 1, MPI_INT, counts, 1, MPI_INT,
				MPI_COMM_WORLD);
		  istart = 0;
		  for (iprocs = 0; iprocs < mype; iprocs++){
			  istart+=counts[iprocs];
		  }
		  
		  for (iprocs = 0; iprocs < mype; iprocs++){
			  index_end += counts[iprocs];
			  if (index >= index_start && index <= index_end){
				  root_pe = iprocs;
			  }
			  index_start+=counts[iprocs];
		  }
		  
		  free(counts);
		  
		  if (mype == root_pe){
			  found_value = iptrinp[index-istart];
		  }
	}
	else {
		found_value = iptrinp[index];
	}

	L7_Broadcast(&found_value, 1, l7_datatype, root_pe);
	
	*((int *)val) = found_value;
	
#else
   *((int *)val) = ((int *)array)[index];
#endif
	return(0);
}

void l7_getglobal_(
      void                   *array,
      const int              *count,
      const enum L7_Datatype *l7_datatype,
		const int              *index,
		void                   *val,
		int                    *ierr
		)
{
	*ierr=L7_GetGlobal(array, *count, *l7_datatype, *index, val);
}

