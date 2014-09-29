#include "mesh.h"
#include "l7.h"
#include "mpi.h"

#if !defined(FULL_PRECISION) && !defined(MIXED_PRECISION) && !defined(MINIMUM_PRECISION)
#define FULL_PRECISION
#endif
#ifdef NO_CL_DOUBLE
#undef  FULL_PRECISION
#undef  MIXED_PRECISION
#define MINIMUM_PRECISION
#endif

#if defined(MINIMUM_PRECISION)
   typedef float real_t; // this is used for intermediate calculations
#elif defined(MIXED_PRECISION) // intermediate values calculated high precision and stored as floats
   typedef double real_t;
#elif defined(FULL_PRECISION)
   typedef double real_t;
#endif

enum partition_method initial_order = HILBERT_SORT;
enum partition_method cycle_reorder = ORIGINAL_ORDER;
bool localStencil;

static Mesh *mesh;

int main (int argc, char **argv)
{
   //  Process command-line arguments, if any.
   int mype=0;
   int numpe=0;
   int do_quo_setup = 0;
   int lttrace_on = 0;
   L7_Init(&mype, &numpe, &argc, argv, do_quo_setup, lttrace_on);

   int nx = 4;
   int ny = 4;
   int levmx = 2;
   int ndim = 2;
   int boundary = 1;
   int parallel_in = 1;
   int do_gpu_calc = 0;
   double circ_radius = 6.0;
   circ_radius *= (double)nx / 128.0; // scaling radius for problem size

   double deltax_in = 1.0;
   double deltay_in = 1.0;
   mesh = new Mesh(nx, ny, levmx, ndim, deltax_in, deltay_in, boundary, parallel_in, do_gpu_calc);

   mesh->init(nx, ny, circ_radius, initial_order, do_gpu_calc);

   MallocPlus state_memory;

   real_t *density = (real_t *)state_memory.memory_malloc(mesh->ncells, sizeof(real_t), "density");
   for (int ic=0; ic<mesh->ncells; ic++){
      density[ic]=1.0;
   }

   mesh->do_load_balance_local(mesh->ncells, NULL, state_memory);

   int ncells_test = mesh->ncells_global/mesh->numpe;
   if (mesh->ncells_global%mesh->numpe > mype) ncells_test++;

   int ierr = 0;
   if (ncells_test != mesh->ncells) ierr = 1;

   int ierr_global = 0;
   MPI_Allreduce(&ierr, &ierr_global, 1, MPI_INT, MPI_SUM, MPI_COMM_WORLD);

   //printf("%d: DEBUG -- ncells %ld ncells_global %ld ncells_test %d ierr %d ierr_global %d\n",mype,mesh->ncells,mesh->ncells_global,ncells_test,ierr,ierr_global);

   if (mype == 0){
      if (ierr_global){
          printf("  Error with load balance\n");
      }
      else{
          printf("  PASSED load balance\n");
      }
   }


   L7_Terminate();
   exit(0);
}
