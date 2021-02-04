#include "mesh/mesh.h"
#include <unistd.h>
#include <stdio.h>
#include <assert.h>
#include <algorithm>
#include <queue>
#include "state.h"
#include "timer/timer.h"
#include "genmalloc/genmalloc.h"

#ifdef HAVE_MPI
#include <mpi.h>
#endif
#include <omp.h>
#undef DEBUG
//#define DEBUG 0
#undef DEBUG_RESTORE_VALS
#define TIMING_LEVEL 2

#if defined(HALF_PRECISION)
#define ZERO 0.0f
#define ONE 1.0f
#define HALF 0.5f
#define EPSILON 1.0f-30
#define STATE_EPS        15.0
// calc refine is done in single precision
#define REFINE_GRADIENT  0.10f
#define COARSEN_GRADIENT 0.05f
#define REFINE_HALF 0.5f
#define REFINE_NEG_THOUSAND -1000.0f

#elif defined(MINIMUM_PRECISION)
#define ZERO 0.0f
#define ONE 1.0f
#define HALF 0.5f
#define EPSILON 1.0f-30
#define STATE_EPS        15.0
// calc refine is done in single precision
#define REFINE_GRADIENT  0.10f
#define COARSEN_GRADIENT 0.05f
#define REFINE_HALF 0.5f
#define REFINE_NEG_THOUSAND -1000.0f

#elif defined(MIXED_PRECISION) // intermediate values calculated high precision and stored as floats
#define ZERO 0.0
#define ONE 1.0
#define HALF 0.5
#define EPSILON 1.0e-30
#define STATE_EPS        .02
// calc refine is done in single precision
#define REFINE_GRADIENT  0.10f
#define COARSEN_GRADIENT 0.05f
#define REFINE_HALF 0.5f
#define REFINE_NEG_THOUSAND -1000.0f

#elif defined(FULL_PRECISION)
#define ZERO 0.0
#define ONE 1.0
#define HALF 0.5
#define EPSILON 1.0e-30
#define STATE_EPS        .001
// calc refine is done in single precision
#define REFINE_GRADIENT  0.10
#define COARSEN_GRADIENT 0.05
#define REFINE_HALF 0.5
#define REFINE_NEG_THOUSAND -1000.0

#endif

#ifdef _OPENMP
//static bool iversion_flag = false;
#endif

bool phantom_debug = false;

typedef unsigned int uint;

static const char *state_timer_descriptor[STATE_TIMER_SIZE] = {
   "state_timer_apply_BCs",
   "state_timer_set_timestep",
   "state_timer_finite_difference",
   "state_timer_finite_diff_part1",
   "state_timer_finite_diff_part2",
   "state_timer_finite_diff_part3",
   "state_timer_finite_diff_part4",
   "state_timer_finite_diff_part5",
   "state_timer_finite_diff_part6",
   "state_timer_refine_potential",
   "state_timer_calc_mpot",
   "state_timer_rezone_all",
   "state_timer_mass_sum",
   "state_timer_read",
   "state_timer_write"
};

const int CRUX_STATE_VERSION = 102;
const int num_int_vals       = 1;

int int_vals[num_int_vals] = {CRUX_STATE_VERSION};

#ifdef HAVE_OPENCL
#include "state_kernel.inc"
#endif

struct esum_type{
   double sum;
   double correction;
};
#ifdef HAVE_MPI
MPI_Datatype MPI_TWO_DOUBLES;
MPI_Op KNUTH_SUM;
int commutative = 1;
void knuth_sum(struct esum_type *in, struct esum_type *inout, int *len, MPI_Datatype *MPI_TWO_DOUBLES);
#endif

int save_ncells;

#define CONSERVED_EQNS

//#define PRECISION_CHECK 1.0e-7
//#define PRECISION_CHECK_STATS 1

#ifdef PRECISION_CHECK_STATS
static int fail_prec_count = 0;
static double fail_F_plus_sum    = 0.0;
static double fail_F_minus_sum   = 0.0;
static double fail_G_plus_sum    = 0.0;
static double fail_G_minus_sum   = 0.0;
static double fail_wminusx_H_sum = 0.0;
static double fail_wplusx_H_sum  = 0.0;
static double fail_wminusy_H_sum = 0.0;
static double fail_wplusy_H_sum  = 0.0;
static int prec_count = 0;
static double F_plus_sum    = 0.0;
static double F_minus_sum   = 0.0;
static double G_plus_sum    = 0.0;
static double G_minus_sum   = 0.0;
static double wminusx_H_sum = 0.0;
static double wplusx_H_sum  = 0.0;
static double wminusy_H_sum = 0.0;
static double wplusy_H_sum  = 0.0;

static int fail_prec_avg_count = 0;
static double fail_F_plus_avg    = 0.0;
static double fail_F_minus_avg   = 0.0;
static double fail_G_plus_avg    = 0.0;
static double fail_G_minus_avg   = 0.0;
static double fail_wminusx_H_avg = 0.0;
static double fail_wplusx_H_avg  = 0.0;
static double fail_wminusy_H_avg = 0.0;
static double fail_wplusy_H_avg  = 0.0;
static int prec_avg_count = 0;
static double F_plus_avg    = 0.0;
static double F_minus_avg   = 0.0;
static double G_plus_avg    = 0.0;
static double G_minus_avg   = 0.0;
static double wminusx_H_avg = 0.0;
static double wplusx_H_avg  = 0.0;
static double wminusy_H_avg = 0.0;
static double wplusy_H_avg  = 0.0;
#endif

#ifdef PRECISION_CHECK
FILE *fprecise;
#endif

#define SQR(x) ( x*x )
#define MIN3(x,y,z) ( min( min(x,y), z) )

void doubleToHex(FILE *fp, double val){
	fprintf(fp, "0x%02x", *(((unsigned char*)(&val))+sizeof(double)-1)  );
 		for(int i=sizeof(double)-2; i >=0; i--) {
			fprintf(fp, "%02x", *(((unsigned char*)(&val))+i)  );
		}
}

#ifdef HAVE_OPENCL
cl_kernel kernel_set_timestep;
cl_kernel kernel_reduction_min;
cl_kernel kernel_copy_state_data;
cl_kernel kernel_copy_state_ghost_data;
cl_kernel kernel_apply_boundary_conditions;
cl_kernel kernel_apply_boundary_conditions_local;
cl_kernel kernel_apply_boundary_conditions_ghost;
cl_kernel kernel_calc_finite_difference;
cl_kernel kernel_calc_finite_difference_via_faces_face;
cl_kernel kernel_calc_finite_difference_via_faces_cell;
cl_kernel kernel_calc_finite_difference_in_place_cell_comps;
cl_kernel kernel_calc_finite_difference_in_place_fixup;
cl_kernel kernel_calc_finite_difference_in_place_fill_new;
cl_kernel kernel_calc_finite_difference_via_face_in_place_face_comps;
cl_kernel kernel_calc_finite_difference_via_face_in_place_fixup;
cl_kernel kernel_calc_finite_difference_via_face_in_place_fill_new;
cl_kernel kernel_calc_finite_difference_regular_cells_comps;
cl_kernel kernel_calc_finite_difference_regular_cells_fill;
cl_kernel kernel_calc_finite_difference_regular_cells_face_comps;
cl_kernel kernel_calc_finite_difference_regular_cells_by_faces_fill;
cl_kernel kernel_refine_potential;
cl_kernel kernel_reduce_sum_mass_stage1of2;
cl_kernel kernel_reduce_sum_mass_stage2of2;
cl_kernel kernel_reduce_epsum_mass_stage1of2;
cl_kernel kernel_reduce_epsum_mass_stage2of2;
#endif


#ifdef _OPENMP_SIMD
#pragma omp declare simd
#endif
inline real_t U_halfstep(// XXX Fix the subindices to be more intuitive XXX
        real_t    deltaT,     // Timestep
        real_t    U_i,        // Initial cell's (downwind's) state variable
        real_t    U_n,        // Next cell's    (upwind's)   state variable
        real_t    F_i,        // Initial cell's (downwind's) state variable flux
        real_t    F_n,        // Next cell's    (upwind's)   state variable flux
        real_t    r_i,        // Initial cell's (downwind's) center to face distance
        real_t    r_n,        // Next cell's    (upwind's)   center to face distance
        real_t    A_i,        // Cell's            face surface area
        real_t    A_n,        // Cell's neighbor's face surface area
        real_t    V_i,        // Cell's            volume
        real_t    V_n) {      // Cell's neighbor's volume

   return (( r_i*U_n + r_n*U_i ) / ( r_i + r_n )) 
          - HALF*deltaT*(( F_n*A_n*min(ONE, A_i/A_n) - F_i*A_i*min(ONE, A_n/A_i) )
                    / ( V_n*min(HALF, V_i/V_n) + V_i*min(HALF, V_n/V_i) ));

}

inline real_t U_fullstep(
        real_t    deltaT,
        real_t    dr,
        real_t    U,
        real_t    F_plus,
        real_t    F_minus,
        real_t    G_plus,
        real_t    G_minus) {

#ifdef PRECISION_CHECK_BEST_PARENTHESIS
   //"best" parentheses version
   //return (U - (deltaT / dr)*((F_plus - F_minus) + (G_plus - G_minus)));
   return (U + (-(deltaT/dr)*((F_plus-F_minus)+(G_plus-G_minus))));
#else
   //original, no parentheses
   return (U - (deltaT / dr)*(F_plus - F_minus + G_plus - G_minus));
#endif

}

#ifdef PRECISION_CHECK

inline void U_fullstep_precision_check(
        int       ic,
        real_t    deltaT_in,
        real_t    dr_in,
        real_t    U_in,
        real_t    U_new_in,
        real_t    F_plus_in,
        real_t    F_minus_in,
        real_t    G_plus_in,
        real_t    G_minus_in,
        real_t    wminusx_H_in,
        real_t    wplusx_H_in,
        real_t    wminusy_H_in,
        real_t    wplusy_H_in,
        int       *fail) {

   float U       = (float)U_in;
   float deltaT  = (float)deltaT_in;
   float dr      = (float)dr_in;
   float F_plus  = (float)F_plus_in;
   float F_minus = (float)F_minus_in;
   float G_plus  = (float)G_plus_in;
   float G_minus = (float)G_minus_in;
   float wminusx_H = (float)wminusx_H_in;
   float wplusx_H = (float)wplusx_H_in;
   float wminusy_H = (float)wminusy_H_in;
   float wplusy_H = (float)wplusy_H_in;

#ifdef PRECISION_CHECK_WITH_PARENTHESIS
   //Some parentheses
   double U_new = U - (deltaT / dr)*(F_plus - F_minus + G_plus - G_minus)
                   +( -wminusx_H + wplusx_H - wminusy_H + wplusy_H);
#else

#ifdef PRECISION_CHECK_BEST_PARENTHESIS
   //"best" parentheses version
   double U_new = U - (deltaT / dr)*((F_plus - F_minus) + (G_plus - G_minus))
                  // + (((-wminusx_H - wminusy_H) + wplusy_H) + wplusx_H);
                   + (((wplusx_H + wplusy_H) - wminusy_H) - wminusx_H);
#else
   //original, no parentheses
   double U_new = U - (deltaT / dr)*(F_plus - F_minus + G_plus - G_minus)
                   + -wminusx_H + wplusx_H - wminusy_H + wplusy_H;
#endif
#endif

   *fail = 0;

   if (fabs(U_new - U_new_in)/U_new > PRECISION_CHECK) {
      fprintf(fprecise, "DEBUG -- found one at ic %d precision diff is %12.6lg relative %12.6lg\n",ic,fabs(U_new - U_new_in), fabs(U_new - U_new_in)/U_new);

      *fail = 1;

#ifdef PRECISION_CHECK_STATS
      fail_prec_count++;
      fail_F_plus_sum    += fabs(F_plus_in);
      fail_F_minus_sum   += fabs(F_minus_in);
      fail_G_plus_sum    += fabs(G_plus_in);
      fail_G_minus_sum   += fabs(G_minus_in);
      fail_wminusx_H_sum += fabs(wminusx_H_in);
      fail_wplusx_H_sum  += fabs(wplusx_H_in);
      fail_wminusy_H_sum += fabs(wminusy_H_in);
      fail_wplusy_H_sum  += fabs(wplusy_H_in);
#endif
   }

#ifdef PRECISION_CHECK_STATS
   prec_count++;
   F_plus_sum    += fabs(F_plus_in);
   F_minus_sum   += fabs(F_minus_in);
   G_plus_sum    += fabs(G_plus_in);
   G_minus_sum   += fabs(G_minus_in);
   wminusx_H_sum += fabs(wminusx_H_in);
   wplusx_H_sum  += fabs(wplusx_H_in);
   wminusy_H_sum += fabs(wminusy_H_in);
   wplusy_H_sum  += fabs(wplusy_H_in);
#endif
}

#endif

#ifdef _OPENMP_SIMD
#pragma omp declare simd
#endif
inline real_t w_corrector(
        real_t    deltaT,       // Timestep
        real_t    dr,           // Cell's center to face distance
        real_t    U_eigen,      // State variable's eigenvalue (speed)
        real_t    grad_half,    // Centered gradient
        real_t    grad_minus,   // Downwind gradient
        real_t    grad_plus) {  // Upwind gradient

   real_t nu     = HALF * U_eigen * deltaT / dr;
   nu          = nu * (ONE - nu);

   real_t rdenom = ONE / max(SQR(grad_half), EPSILON);
   real_t rplus  = (grad_plus  * grad_half) * rdenom;
   real_t rminus = (grad_minus * grad_half) * rdenom;

   return HALF*nu*(ONE- max(MIN3(ONE, rplus, rminus), ZERO));
}

inline real_t U_reggrid_halfstep(// XXX Fix the subindices to be more intuitive XXX
        real_t    deltaT,     // Timestep
        real_t    deltax,     // Cell size in direction of flux
        real_t    U_i,        // Initial cell's (downwind's) state variable
        real_t    U_n,        // Next cell's    (upwind's)   state variable
        real_t    F_i,        // Initial cell's (downwind's) state variable flux
        real_t    F_n) {      // Next cell's    (upwind's)   state variable flux

   return ( HALF*( ((U_n) + (U_i)) - (deltaT)/(deltax)*((F_n) - (F_i)) ) );
}


State::State(Mesh *mesh_in)
{
   state_memory.memory_add(int_vals,   (size_t)num_int_vals,     4, "state_int_vals",   RESTART_DATA | REPLICATED_DATA);
   state_memory.memory_add(cpu_timers, (size_t)STATE_TIMER_SIZE, 8, "state_cpu_timers", RESTART_DATA);
   state_memory.memory_add(gpu_timers, (size_t)STATE_TIMER_SIZE, 8, "state_gpu_timers", RESTART_DATA);
   for (int i = 0; i < STATE_TIMER_SIZE; i++){
      cpu_timers[i] = 0.0;
   }
   for (int i = 0; i < STATE_TIMER_SIZE; i++){
      gpu_timers[i] = 0L;
   }

   mesh = mesh_in;

#ifdef HAVE_MPI
   int mpi_init;
   MPI_Initialized(&mpi_init);
   if (mpi_init){
      MPI_Type_contiguous(2, MPI_DOUBLE, &MPI_TWO_DOUBLES);
      MPI_Type_commit(&MPI_TWO_DOUBLES);
      MPI_Op_create((MPI_User_function *)knuth_sum, commutative, &KNUTH_SUM);
      // FIXME add fini and set size
      if (mesh->parallel) state_memory.pinit(MPI_COMM_WORLD, 2L * 1024 * 1024 * 1024);
   }
#endif
}

void State::init(int do_gpu_calc)
{
#ifdef PRECISION_CHECK
   fprecise = fopen("precision.out","w");
#endif
   if (do_gpu_calc) {
#ifdef HAVE_OPENCL
      cl_context context = ezcl_get_context();

      if (mesh->mype == 0) printf("Starting compile of kernels in state\n");
      const char *defines = NULL;
      cl_program program                 = ezcl_create_program_wsource(context, defines, state_kern_source);

      kernel_set_timestep                    = ezcl_create_kernel_wprogram(program, "set_timestep_cl");
      kernel_reduction_min                   = ezcl_create_kernel_wprogram(program, "finish_reduction_min_cl");
      kernel_copy_state_data                 = ezcl_create_kernel_wprogram(program, "copy_state_data_cl");
      kernel_copy_state_ghost_data           = ezcl_create_kernel_wprogram(program, "copy_state_ghost_data_cl");
      kernel_apply_boundary_conditions       = ezcl_create_kernel_wprogram(program, "apply_boundary_conditions_cl");
      kernel_apply_boundary_conditions_local = ezcl_create_kernel_wprogram(program, "apply_boundary_conditions_local_cl");
      kernel_apply_boundary_conditions_ghost = ezcl_create_kernel_wprogram(program, "apply_boundary_conditions_ghost_cl");
      kernel_calc_finite_difference          = ezcl_create_kernel_wprogram(program, "calc_finite_difference_cl");
      kernel_calc_finite_difference_via_faces_face = ezcl_create_kernel_wprogram(program, "calc_finite_difference_via_faces_face_comps_cl");
      kernel_calc_finite_difference_via_faces_cell = ezcl_create_kernel_wprogram(program, "calc_finite_difference_via_faces_cell_comps_cl");
      kernel_calc_finite_difference_in_place_cell_comps = ezcl_create_kernel_wprogram(program, "calc_finite_difference_in_place_cell_comps_cl");
      kernel_calc_finite_difference_in_place_fixup = ezcl_create_kernel_wprogram(program, "calc_finite_difference_in_place_fixup_cl");
      kernel_calc_finite_difference_in_place_fill_new = ezcl_create_kernel_wprogram(program, "calc_finite_difference_in_place_fill_new_cl");
      kernel_calc_finite_difference_via_face_in_place_face_comps = ezcl_create_kernel_wprogram(program, "calc_finite_difference_via_face_in_place_face_comps_cl");
      kernel_calc_finite_difference_via_face_in_place_fixup = ezcl_create_kernel_wprogram(program, "calc_finite_difference_via_face_in_place_fixup_cl");
      kernel_calc_finite_difference_via_face_in_place_fill_new = ezcl_create_kernel_wprogram(program, "calc_finite_difference_via_face_in_place_fill_new_cl");
      kernel_calc_finite_difference_regular_cells_comps = ezcl_create_kernel_wprogram(program, "calc_finite_difference_regular_cells_comps_cl");
      kernel_calc_finite_difference_regular_cells_fill = ezcl_create_kernel_wprogram(program, "calc_finite_difference_regular_cells_fill_cl");
      kernel_calc_finite_difference_regular_cells_face_comps = ezcl_create_kernel_wprogram(program, "calc_finite_difference_regular_cells_face_comps_cl");
      kernel_calc_finite_difference_regular_cells_by_faces_fill = ezcl_create_kernel_wprogram(program, "calc_finite_difference_regular_cells_by_faces_fill_cl");
      kernel_refine_potential                = ezcl_create_kernel_wprogram(program, "refine_potential_cl");
      kernel_reduce_sum_mass_stage1of2       = ezcl_create_kernel_wprogram(program, "reduce_sum_mass_stage1of2_cl");
      kernel_reduce_sum_mass_stage2of2       = ezcl_create_kernel_wprogram(program, "reduce_sum_mass_stage2of2_cl");
      kernel_reduce_epsum_mass_stage1of2     = ezcl_create_kernel_wprogram(program, "reduce_epsum_mass_stage1of2_cl");
      kernel_reduce_epsum_mass_stage2of2     = ezcl_create_kernel_wprogram(program, "reduce_epsum_mass_stage2of2_cl");

      ezcl_program_release(program);
      if (mesh->mype == 0) printf("Finishing compile of kernels in state\n");
#endif

   }

   //printf("\nDEBUG -- Calling state memory memory malloc at line %d\n",__LINE__);
   allocate(mesh->ncells);
   //state_memory.memory_report();
   //printf("DEBUG -- Finished state memory memory malloc at line %d\n\n",__LINE__);

}

void State::allocate(size_t ncells)
{
   int flags = 0;
   flags = (RESTART_DATA | REZONE_DATA | LOAD_BALANCE_MEMORY);
   //if (mesh->parallel) flags = (flags | LOAD_BALANCE_MEMORY);

   H = (state_t *)state_memory.memory_malloc(ncells, sizeof(state_t), "H", flags);
   U = (state_t *)state_memory.memory_malloc(ncells, sizeof(state_t), "U", flags);
   V = (state_t *)state_memory.memory_malloc(ncells, sizeof(state_t), "V", flags);
#ifdef PRECISION_CHECK_GRAPHICS
   PCHECK = (state_t *)state_memory.memory_malloc(ncells, sizeof(state_t), "PCHECK", flags);
#endif
}

void State::resize(size_t new_ncells){
   size_t current_size = state_memory.get_memory_size(H);
   if (new_ncells > current_size) state_memory.memory_realloc_all(new_ncells);

   //printf("\nDEBUG -- Calling state memory resize at line %d\n",__LINE__);
   //state_memory.memory_report();
   //printf("DEBUG -- Finished state memory resize at line %d\n\n",__LINE__);
}

void State::memory_reset_ptrs(void){
   H = (state_t *)state_memory.get_memory_ptr("H");
   U = (state_t *)state_memory.get_memory_ptr("U");
   V = (state_t *)state_memory.get_memory_ptr("V");
#ifdef PRECISION_CHECK_GRAPHICS
   PCHECK = (state_t *)state_memory.get_memory_ptr("PCHECK");
#endif

   //printf("\nDEBUG -- Calling state memory reset_ptrs at line %d\n",__LINE__);
   //state_memory.memory_report();
   //printf("DEBUG -- Finished state memory reset_ptrs at line %d\n\n",__LINE__);
}

#ifdef HAVE_OPENCL
void State::gpu_memory_reset_ptrs(void)
{
   dev_H = (cl_mem)gpu_state_memory.get_memory_ptr("dev_H");
   dev_U = (cl_mem)gpu_state_memory.get_memory_ptr("dev_U");
   dev_V = (cl_mem)gpu_state_memory.get_memory_ptr("dev_V");
}
#endif

void State::terminate(void)
{
#ifdef PRECISION_CHECK
   fclose(fprecise);
#endif
#ifdef PRECISION_CHECK_STATS
   printf("Stats are Fplus %lf Fminus %lf Gplus %lf Gminus %lf wminusx %lf wplusx %lf wminusy %lf wplusy %lf\n",
      F_plus_avg/(double)prec_avg_count,
      F_minus_avg/(double)prec_avg_count,
      G_plus_avg/(double)prec_avg_count,
      G_minus_avg/(double)prec_avg_count,
      wminusx_H_avg/(double)prec_avg_count,
      wplusx_H_avg/(double)prec_avg_count,
      wminusy_H_avg/(double)prec_avg_count,
      wplusy_H_avg/(double)prec_avg_count);
   printf("Stats for fails are Fplus %lf Fminus %lf Gplus %lf Gminus %lf wminusx %lf wplusx %lf wminusy %lf wplusy %lf\n",
      fail_F_plus_avg/(double)fail_prec_avg_count,
      fail_F_minus_avg/(double)fail_prec_avg_count,
      fail_G_plus_avg/(double)fail_prec_avg_count,
      fail_G_minus_avg/(double)fail_prec_avg_count,
      fail_wminusx_H_avg/(double)fail_prec_avg_count,
      fail_wplusx_H_avg/(double)fail_prec_avg_count,
      fail_wminusy_H_avg/(double)fail_prec_avg_count,
      fail_wplusy_H_avg/(double)fail_prec_avg_count);
#endif

   state_memory.memory_delete(H);
   state_memory.memory_delete(U);
   state_memory.memory_delete(V);
   state_memory.memory_remove(int_vals);
   state_memory.memory_remove(cpu_timers);
   state_memory.memory_remove(gpu_timers);

#ifdef HAVE_OPENCL
   ezcl_device_memory_delete(dev_deltaT);

   gpu_state_memory.memory_delete(dev_H);
   gpu_state_memory.memory_delete(dev_U);
   gpu_state_memory.memory_delete(dev_V);

   ezcl_kernel_release(kernel_set_timestep);
   ezcl_kernel_release(kernel_reduction_min);
   ezcl_kernel_release(kernel_copy_state_data);
   ezcl_kernel_release(kernel_copy_state_ghost_data);
   ezcl_kernel_release(kernel_apply_boundary_conditions);
   ezcl_kernel_release(kernel_apply_boundary_conditions_local);
   ezcl_kernel_release(kernel_apply_boundary_conditions_ghost);
   ezcl_kernel_release(kernel_calc_finite_difference);
   ezcl_kernel_release(kernel_calc_finite_difference_via_faces_face);
   ezcl_kernel_release(kernel_calc_finite_difference_via_faces_cell);
   ezcl_kernel_release(kernel_calc_finite_difference_in_place_cell_comps);
   ezcl_kernel_release(kernel_calc_finite_difference_in_place_fixup);
   ezcl_kernel_release(kernel_calc_finite_difference_in_place_fill_new);
   ezcl_kernel_release(kernel_calc_finite_difference_via_face_in_place_face_comps);
   ezcl_kernel_release(kernel_calc_finite_difference_via_face_in_place_fixup);
   ezcl_kernel_release(kernel_calc_finite_difference_via_face_in_place_fill_new);
   ezcl_kernel_release(kernel_calc_finite_difference_regular_cells_comps);
   ezcl_kernel_release(kernel_calc_finite_difference_regular_cells_fill);
   ezcl_kernel_release(kernel_calc_finite_difference_regular_cells_face_comps);
   ezcl_kernel_release(kernel_calc_finite_difference_regular_cells_by_faces_fill);

   ezcl_kernel_release(kernel_refine_potential);
   ezcl_kernel_release(kernel_reduce_sum_mass_stage1of2);
   ezcl_kernel_release(kernel_reduce_sum_mass_stage2of2);
   ezcl_kernel_release(kernel_reduce_epsum_mass_stage1of2);
   ezcl_kernel_release(kernel_reduce_epsum_mass_stage2of2);
#endif
#ifdef HAVE_MPI
   if (mesh->parallel) state_memory.pfini();
#endif
}

#ifdef HAVE_MPI
void knuth_sum(struct esum_type *in, struct esum_type *inout, int *len, MPI_Datatype *MPI_TWO_DOUBLES)
{
   double u, v, upt, up, vpp;
   u = inout->sum;
   v = in->sum + (in->correction+inout->correction);
   upt = u + v;
   up = upt - v;
   vpp = upt - up;
   inout->sum = upt;
   inout->correction = (u - up) + (v - vpp);

   // Just to block compiler warnings
   //if (1==2) printf("DEBUG len %d datatype %lld\n",*len,(long long)(*MPI_TWO_DOUBLES) );
}
#endif

void State::add_boundary_cells(void)
{
   struct timespec tstart_cpu;
   cpu_timer_start(&tstart_cpu);

   // This is for a mesh with no boundary cells -- they are added and
   // the mesh sizes increased
   size_t &ncells        = mesh->ncells;
   vector<int>  &index    = mesh->index;
   vector<spatial_t> &x        = mesh->x;
   vector<spatial_t> &dx       = mesh->dx;
   vector<spatial_t> &y        = mesh->y;
   vector<spatial_t> &dy       = mesh->dy;

   int *i        = mesh->i;
   int *j        = mesh->j;
   uchar_t *level    = mesh->level;
   char_t *celltype = mesh->celltype;
   int *nlft     = mesh->nlft;
   int *nrht     = mesh->nrht;
   int *nbot     = mesh->nbot;
   int *ntop     = mesh->ntop;

   vector<int> &lev_ibegin = mesh->lev_ibegin;
   vector<int> &lev_iend   = mesh->lev_iend;
   vector<int> &lev_jbegin = mesh->lev_jbegin;
   vector<int> &lev_jend   = mesh->lev_jend;

   // Pre-count number of cells to add
   int icount = 0;
   for (uint ic=0; ic<ncells; ic++) {
      if (i[ic] == lev_ibegin[level[ic]]) icount++; // Left boundary
      if (i[ic] == lev_iend[level[ic]])   icount++; // Right boundary
      if (j[ic] == lev_jbegin[level[ic]]) icount++; // Bottom boundary
      if (j[ic] == lev_jend[level[ic]])   icount++; // Top boundary
   }
      
   int new_ncells = ncells + icount;
   // Increase the arrays for the new boundary cells
   H=(state_t *)state_memory.memory_realloc(new_ncells, H);
   U=(state_t *)state_memory.memory_realloc(new_ncells, U);
   V=(state_t *)state_memory.memory_realloc(new_ncells, V);
   //printf("\nDEBUG add_boundary cells\n"); 
   //state_memory.memory_report();
   //printf("DEBUG end add_boundary cells\n\n"); 

   mesh->i        =(int *)mesh->mesh_memory.memory_realloc(new_ncells, i);
   mesh->j        =(int *)mesh->mesh_memory.memory_realloc(new_ncells, j);
   mesh->level    =(uchar_t *)mesh->mesh_memory.memory_realloc(new_ncells, level);
   // needs to cast char_t to void so doesn't mistake it for string
   mesh->celltype =(char_t *)mesh->mesh_memory.memory_realloc(new_ncells, (void *)celltype);
   mesh->nlft     =(int *)mesh->mesh_memory.memory_realloc(new_ncells, nlft);
   mesh->nrht     =(int *)mesh->mesh_memory.memory_realloc(new_ncells, nrht);
   mesh->nbot     =(int *)mesh->mesh_memory.memory_realloc(new_ncells, nbot);
   mesh->ntop     =(int *)mesh->mesh_memory.memory_realloc(new_ncells, ntop);
   //memory_reset_ptrs();
   i        = mesh->i;
   j        = mesh->j;
   level    = mesh->level;
   celltype = mesh->celltype;
   nlft     = mesh->nlft;
   nrht     = mesh->nrht;
   nbot     = mesh->nbot;
   ntop     = mesh->ntop;

   index.resize(new_ncells);
   x.resize(new_ncells);
   dx.resize(new_ncells);
   y.resize(new_ncells);
   dy.resize(new_ncells);

#ifdef _OPENMP_SIMD
#pragma omp simd
#endif
   for (int nc=ncells; nc<new_ncells; nc++) {
      nlft[nc] = -1;
      nrht[nc] = -1;
      nbot[nc] = -1;
      ntop[nc] = -1;
   }
      
   // In the first pass, set two of the neighbor indices and all
   // the other data to be brought across. Set the inverse of the
   // the velocity to enforce the reflective boundary condition
   uint nc=ncells;
   for (uint ic=0; ic<ncells; ic++) {
      if (i[ic] == lev_ibegin[level[ic]]) {
         nlft[ic] = nc;
         nlft[nc] = nc;
         nrht[nc] = ic;
         i[nc] = lev_ibegin[level[ic]]-1;
         j[nc] = j[ic];
         level[nc] = level[ic];
         dx[nc] = dx[ic];
         dy[nc] = dy[ic];
         x[nc] = x[ic]-dx[ic];
         y[nc] = y[ic];
         H[nc] =  H[ic];
         U[nc] = -U[ic];
         V[nc] =  V[ic];
         nc++;
      }
      if (i[ic] == lev_iend[level[ic]]) {
         nrht[ic] = nc;
         nrht[nc] = nc;
         nlft[nc] = ic;
         i[nc] = lev_iend[level[ic]]+1;
         j[nc] = j[ic];
         level[nc] = level[ic];
         dx[nc] = dx[ic];
         dy[nc] = dy[ic];
         x[nc] = x[ic]+dx[ic];
         y[nc] = y[ic];
         H[nc] =  H[ic];
         U[nc] = -U[ic];
         V[nc] =  V[ic];
         nc++;
      }
      if (j[ic] == lev_jbegin[level[ic]]) {
         nbot[ic] = nc;
         nbot[nc] = nc;
         ntop[nc] = ic;
         i[nc] = i[ic];
         j[nc] = lev_jbegin[level[ic]]-1;
         level[nc] = level[ic];
         dx[nc] = dx[ic];
         dy[nc] = dy[ic];
         x[nc] = x[ic];
         y[nc] = y[ic]-dy[ic];
         H[nc] =  H[ic];
         U[nc] =  U[ic];
         V[nc] = -V[ic];
         nc++;
      }
      if (j[ic] == lev_jend[level[ic]]) {
         ntop[ic] = nc;
         ntop[nc] = nc;
         nbot[nc] = ic;
         i[nc] = i[ic];
         j[nc] = lev_jend[level[ic]]+1;
         level[nc] = level[ic];
         dx[nc] = dx[ic];
         dy[nc] = dy[ic];
         x[nc] = x[ic];
         y[nc] = y[ic]+dy[ic];
         H[nc] =  H[ic];
         U[nc] =  U[ic];
         V[nc] = -V[ic];
         nc++;
      }
   }

   // Now set the other two neighbor indices
   for (int nc=ncells; nc<new_ncells; nc++) {
      if (i[nc] == lev_ibegin[level[nc]]-1) {
         // Need to check if also a bottom boundary cell
         if (j[nc] == lev_jbegin[level[nc]]){
           nbot[nc] = nc;
         } else {
           nbot[nc] = nlft[nbot[nrht[nc]]];
         }
         if (j[nc] == lev_jend[level[nc]]){
           ntop[nc] = nc;
         } else {
           ntop[nc] = nlft[ntop[nrht[nc]]];
         }
      }
      if (i[nc] == lev_iend[level[nc]]+1)   {
         if (level[nc] <= level[nbot[nlft[nc]]]){
            if (j[nc] == lev_jbegin[level[nc]]){
               nbot[nc] = nc;
            } else {
               nbot[nc] = nrht[nbot[nlft[nc]]];
            }
            if (j[nc] == lev_jend[level[nc]]){
               ntop[nc] = nc;
            } else {
               ntop[nc] = nrht[ntop[nlft[nc]]];
            }
         // calculation is a little different if going through a
         // finer zoned region
         } else {
            nbot[nc] = nrht[nrht[nbot[nlft[nc]]]];
            ntop[nc] = nrht[nrht[ntop[nlft[nc]]]];
         }
      }
      if (j[nc] == lev_jbegin[level[nc]]-1) {
         if (i[nc] == lev_ibegin[level[nc]]){
            nlft[nc] = nc;
         } else {
            nlft[nc] = nbot[nlft[ntop[nc]]];
         }
         if (i[nc] == lev_iend[level[nc]]){
            nrht[nc] = nc;
         } else {
            nrht[nc] = nbot[nrht[ntop[nc]]];
         }
      }
      if (j[nc] == lev_jend[level[nc]]+1)   {
         if (level[nc] <= level[nlft[nbot[nc]]]){
            if (i[nc] == lev_ibegin[level[nc]]){
               nlft[nc] = nc;
            } else {
               nlft[nc] = ntop[nlft[nbot[nc]]];
            }
            if (i[nc] == lev_iend[level[nc]]){
               nrht[nc] = nc;
            } else {
               nrht[nc] = ntop[nrht[nbot[nc]]];
            }
         } else {
            nlft[nc] = ntop[ntop[nlft[nbot[nc]]]];
            nrht[nc] = ntop[ntop[nrht[nbot[nc]]]];
         }
      }
   }
   save_ncells = ncells;
   ncells = new_ncells;

   cpu_timers[STATE_TIMER_APPLY_BCS] += cpu_timer_stop(tstart_cpu);
}


void State::apply_boundary_conditions(void)
{
   int *nlft = mesh->nlft;
   int *nrht = mesh->nrht;
   int *nbot = mesh->nbot;
   int *ntop = mesh->ntop;

#ifdef _OPENMP
#pragma omp master
   {
#endif
   if (mesh->ncells_ghost < mesh->ncells) mesh->ncells_ghost = mesh->ncells;
#ifdef _OPENMP
      }    
#pragma omp barrier
#endif

   // This is for a mesh with boundary cells
   int lowerBound, upperBound;
   mesh->get_bounds(lowerBound, upperBound);

   for (int ic=lowerBound; ic<upperBound; ic++) {
      if (mesh->is_left_boundary(ic)) {
         int nr = nrht[ic];
         if (nr < (int)mesh->ncells) {
            H[ic] =  H[nr];
            U[ic] = -U[nr];
            V[ic] =  V[nr];
         }
      }
      if (mesh->is_right_boundary(ic))  {
         int nl = nlft[ic];
         if (nl < (int)mesh->ncells) {
            H[ic] =  H[nl];
            U[ic] = -U[nl];
            V[ic] =  V[nl];
         }
      }
      if (mesh->is_bottom_boundary(ic)) {
         int nt = ntop[ic];
         if (nt < (int)mesh->ncells) {
            H[ic] =  H[nt];
            U[ic] =  U[nt];
            V[ic] = -V[nt];
         }
      }
      if (mesh->is_top_boundary(ic)) {
         int nb = nbot[ic];
         if (nb < (int)mesh->ncells) {
            H[ic] =  H[nb];
            U[ic] =  U[nb];
            V[ic] = -V[nb];
         }
      }
   }

   if (mesh->numpe > 1) {

#ifdef HAVE_MPI

#ifdef _OPENMP
#pragma omp barrier
#pragma omp master
      {    
#endif
      H=(state_t *)state_memory.memory_realloc(mesh->ncells_ghost, H);
      U=(state_t *)state_memory.memory_realloc(mesh->ncells_ghost, U);
      V=(state_t *)state_memory.memory_realloc(mesh->ncells_ghost, V);

      L7_Update(&H[0], L7_STATE_T, mesh->cell_handle);
      L7_Update(&U[0], L7_STATE_T, mesh->cell_handle);
      L7_Update(&V[0], L7_STATE_T, mesh->cell_handle);
#ifdef _OPENMP
      }    
#pragma omp barrier
#endif

#endif

      // This is for a mesh with boundary cells
      for (int ic=lowerBound; ic<upperBound; ic++) {
         if (mesh->is_left_boundary(ic)) {
            int nr = nrht[ic];
            if (nr >= (int)mesh->ncells) {
               H[ic] =  H[nr];
               U[ic] = -U[nr];
               V[ic] =  V[nr];
            }
         }
         if (mesh->is_right_boundary(ic))  {
            int nl = nlft[ic];
            if (nl >= (int)mesh->ncells) {
               H[ic] =  H[nl];
               U[ic] = -U[nl];
               V[ic] =  V[nl];
            }
         }
         if (mesh->is_bottom_boundary(ic)) {
            int nt = ntop[ic];
            if (nt >= (int)mesh->ncells) {
               H[ic] =  H[nt];
               U[ic] =  U[nt];
               V[ic] = -V[nt];
            }
         }
         if (mesh->is_top_boundary(ic)) {
            int nb = nbot[ic];
            if (nb >= (int)mesh->ncells) {
               H[ic] =  H[nb];
               U[ic] =  U[nb];
               V[ic] = -V[nb];
            }
         }
      }
   }
}

void State::remove_boundary_cells(void)
{
   if(! mesh->have_boundary) {

#ifdef _OPENMP
#pragma omp barrier
#pragma omp master
      {
#endif
         size_t &ncells = mesh->ncells;

         // Resize to drop all the boundary cells
         ncells = save_ncells;
         H=(state_t *)state_memory.memory_realloc(save_ncells, H);
         U=(state_t *)state_memory.memory_realloc(save_ncells, U);
         V=(state_t *)state_memory.memory_realloc(save_ncells, V);
         //printf("\nDEBUG remove_boundary cells\n"); 
         //state_memory.memory_report();
         //printf("DEBUG end remove_boundary cells\n\n"); 

         mesh->i        = (int *)mesh->mesh_memory.memory_realloc(save_ncells, mesh->i);
         mesh->j        = (int *)mesh->mesh_memory.memory_realloc(save_ncells, mesh->j);
         mesh->level    = (uchar_t *)mesh->mesh_memory.memory_realloc(save_ncells, mesh->level);
         // needs to cast char_t to void so doesn't mistake it for string
         mesh->celltype = (char_t *)mesh->mesh_memory.memory_realloc(save_ncells, (void *)mesh->celltype);
         mesh->nlft     = (int *)mesh->mesh_memory.memory_realloc(save_ncells, mesh->nlft);
         mesh->nrht     = (int *)mesh->mesh_memory.memory_realloc(save_ncells, mesh->nrht);
         mesh->nbot     = (int *)mesh->mesh_memory.memory_realloc(save_ncells, mesh->nbot);
         mesh->ntop     = (int *)mesh->mesh_memory.memory_realloc(save_ncells, mesh->ntop);

         // Reset the neighbors due to the dropped boundary cells
         mesh->index.resize(save_ncells);
         mesh->x.resize(save_ncells);
         mesh->dx.resize(save_ncells);
         mesh->y.resize(save_ncells);
         mesh->dy.resize(save_ncells);
#ifdef _OPENMP
      }
#pragma omp barrier
#endif

      mesh->set_bounds(mesh->ncells);

      int lowerBound, upperBound;
      mesh->get_bounds(lowerBound, upperBound);
      for (int ic=lowerBound; ic<upperBound; ic++) {
         if (mesh->i[ic] == mesh->lev_ibegin[mesh->level[ic]]) mesh->nlft[ic] = ic;
         if (mesh->i[ic] == mesh->lev_iend[mesh->level[ic]])   mesh->nrht[ic] = ic;
         if (mesh->j[ic] == mesh->lev_jbegin[mesh->level[ic]]) mesh->nbot[ic] = ic;
         if (mesh->j[ic] == mesh->lev_jend[mesh->level[ic]])   mesh->ntop[ic] = ic;
      }

   } // if have_boundary
}

double State::set_timestep(double g, double sigma)
{
   double globalmindeltaT;

   struct timespec tstart_cpu;
   cpu_timer_start(&tstart_cpu);

   static double mindeltaT;

   int lowerBounds, upperBounds;
   mesh->set_bounds(mesh->ncells);
   mesh->get_bounds(lowerBounds, upperBounds);

#ifdef _OPENMP
#pragma omp barrier
#pragma omp master
   {
#endif
      mindeltaT = 1000;
#ifdef _OPENMP
   }
#pragma omp barrier
#endif

   double mymindeltaT = 1000.0; // private for each thread
#ifdef _OPENMP_SIMD
#pragma omp simd reduction(min:mymindeltaT)
#endif
   for (int ic=lowerBounds; ic<upperBounds; ic++) {
      if (mesh->celltype[ic] == REAL_CELL) {
         uchar_t lev = mesh->level[ic];
         double wavespeed = sqrt(g*H[ic]);
         double xspeed = (fabs(U[ic])+wavespeed)/mesh->lev_deltax[lev];
         double yspeed = (fabs(V[ic])+wavespeed)/mesh->lev_deltay[lev];
         double deltaT=sigma/(xspeed+yspeed);
         if (deltaT < mymindeltaT) mymindeltaT = deltaT;
      }
   }

#ifdef _OPENMP
#pragma omp critical
   {
#endif
      if (mymindeltaT < mindeltaT) mindeltaT = mymindeltaT;
#ifdef _OPENMP
   } // End critical region
#pragma omp barrier
#endif

#ifdef _OPENMP
#pragma omp master
   {
#endif


   globalmindeltaT = mindeltaT;
#ifdef HAVE_MPI
      if (mesh->parallel) MPI_Allreduce(&mindeltaT, &globalmindeltaT, 1, MPI_DOUBLE, MPI_MIN, MPI_COMM_WORLD);
#endif

      cpu_timers[STATE_TIMER_SET_TIMESTEP] += cpu_timer_stop(tstart_cpu);
#ifdef _OPENMP
   } // End master region
#pragma omp barrier
#endif

   return(globalmindeltaT);
}

#ifdef HAVE_OPENCL
double State::gpu_set_timestep(double sigma)
{
   double deltaT, globalmindeltaT;

   struct timespec tstart_cpu;
   cpu_timer_start(&tstart_cpu);

   cl_command_queue command_queue = ezcl_get_command_queue();

   size_t &ncells       = mesh->ncells;
#ifdef HAVE_MPI
   int &parallel        = mesh->parallel;
#endif
   cl_mem &dev_level    = mesh->dev_level;
   cl_mem &dev_celltype = mesh->dev_celltype;
   cl_mem &dev_levdx    = mesh->dev_levdx;
   cl_mem &dev_levdy    = mesh->dev_levdy;

   assert(dev_H);
   assert(dev_U);
   assert(dev_V);
   assert(dev_level);
   assert(dev_celltype);
   assert(dev_levdx);
   assert(dev_levdy);

   size_t local_work_size = 128;
   size_t global_work_size = ((ncells+local_work_size - 1) /local_work_size) * local_work_size;
   size_t block_size     = global_work_size/local_work_size;

   cl_mem dev_redscratch = ezcl_malloc(NULL, const_cast<char *>("dev_redscratch"), &block_size, sizeof(cl_real_t), CL_MEM_READ_WRITE, 0);

      /*
      __kernel void set_timestep_cl(
                       const int       ncells,     // 0  Total number of cells.
                       const real_t    sigma,      // 1
              __global const state_t  *H,          // 2
              __global const state_t  *U,          // 3
              __global const state_t  *V,          // 4
              __global const uchar_t  *level,      // 5  Array of level information.
              __global const char_t   *celltype,   // 6
              __global const real_t   *lev_dx,     // 7
              __global const real_t   *lev_dy,     // 8
              __global       real_t   *redscratch, // 9
              __global       real_t   *deltaT,     // 10
              __local        real_t   *tile)       // 11
      */

   real_t sigma_local = sigma;
   ezcl_set_kernel_arg(kernel_set_timestep,  0, sizeof(cl_int),  (void *)&ncells);
   ezcl_set_kernel_arg(kernel_set_timestep,  1, sizeof(cl_real_t), (void *)&sigma_local);
   ezcl_set_kernel_arg(kernel_set_timestep,  2, sizeof(cl_mem),  (void *)&dev_H);
   ezcl_set_kernel_arg(kernel_set_timestep,  3, sizeof(cl_mem),  (void *)&dev_U);
   ezcl_set_kernel_arg(kernel_set_timestep,  4, sizeof(cl_mem),  (void *)&dev_V);
   ezcl_set_kernel_arg(kernel_set_timestep,  5, sizeof(cl_mem),  (void *)&dev_level);
   ezcl_set_kernel_arg(kernel_set_timestep,  6, sizeof(cl_mem),  (void *)&dev_celltype);
   ezcl_set_kernel_arg(kernel_set_timestep,  7, sizeof(cl_mem),  (void *)&dev_levdx);
   ezcl_set_kernel_arg(kernel_set_timestep,  8, sizeof(cl_mem),  (void *)&dev_levdy);
   ezcl_set_kernel_arg(kernel_set_timestep,  9, sizeof(cl_mem),  (void *)&dev_redscratch);
   ezcl_set_kernel_arg(kernel_set_timestep, 10, sizeof(cl_mem),  (void *)&dev_deltaT);
   ezcl_set_kernel_arg(kernel_set_timestep, 11, local_work_size*sizeof(cl_real_t),  NULL);

   ezcl_enqueue_ndrange_kernel(command_queue, kernel_set_timestep, 1, NULL, &global_work_size, &local_work_size, NULL);

   if (block_size > 1){
         /*
         __kernel void finish_reduction_min_cl(
           const    int      isize,
           __global real_t  *redscratch,
           __global real_t  *deltaT,
           __local  real_t  *tile)
         */
      ezcl_set_kernel_arg(kernel_reduction_min, 0, sizeof(cl_int),  (void *)&block_size);
      ezcl_set_kernel_arg(kernel_reduction_min, 1, sizeof(cl_mem),  (void *)&dev_redscratch);
      ezcl_set_kernel_arg(kernel_reduction_min, 2, sizeof(cl_mem),  (void *)&dev_deltaT);
      ezcl_set_kernel_arg(kernel_reduction_min, 3, local_work_size*sizeof(cl_real_t), NULL);

     ezcl_enqueue_ndrange_kernel(command_queue, kernel_reduction_min, 1, NULL, &local_work_size, &local_work_size, NULL);
   }

   real_t deltaT_local;
   ezcl_enqueue_read_buffer(command_queue, dev_deltaT, CL_TRUE,  0, sizeof(cl_real_t), &deltaT_local, NULL);
   deltaT = deltaT_local;

   globalmindeltaT = deltaT;
#ifdef HAVE_MPI
   if (parallel) MPI_Allreduce(&deltaT, &globalmindeltaT, 1, MPI_DOUBLE, MPI_MIN, MPI_COMM_WORLD);
#endif

   ezcl_device_memory_delete(dev_redscratch);

   gpu_timers[STATE_TIMER_SET_TIMESTEP] += (long)(cpu_timer_stop(tstart_cpu)*1.0e9);

   return(globalmindeltaT);
}
#endif

void State::fill_circle(double  circ_radius,//  Radius of circle in grid units.
                        double  fill_value, //  Circle height for shallow water.
                        double  background) //  Background height for shallow water.
{  
   size_t &ncells = mesh->ncells;
   vector<spatial_t> &x  = mesh->x;
   vector<spatial_t> &dx = mesh->dx;
   vector<spatial_t> &y  = mesh->y;
   vector<spatial_t> &dy = mesh->dy;

   for (uint ic = 0; ic < ncells; ic++)
   {  H[ic] = background;
      U[ic] = V[ic] = 0.0; }
   
   //  Clear the old k-D tree and generate new data (slow but necessary here).
   //KDTree_Destroy(&mesh->tree);
   mesh->kdtree_setup();
   
   int nez;
   vector<int>    ind(ncells);
   vector<double> weight(ncells);
   
#ifdef FULL_PRECISION
   KDTree_QueryCircleInterior_Double(&mesh->tree, &nez, &(ind[0]), circ_radius, ncells,
                                     &x[0], &dx[0],
                                     &y[0], &dy[0]);
#else
   KDTree_QueryCircleInterior_Float(&mesh->tree, &nez, &(ind[0]), circ_radius, ncells,
                                    &x[0], &dx[0],
                                    &y[0], &dy[0]);
#endif
   for (int ic = 0; ic < nez; ++ic)
   {  H[ind[ic]] = fill_value; }
   
#ifdef FULL_PRECISION
   KDTree_QueryCircleIntersectWeighted_Double(&mesh->tree, &nez, &(ind[0]), &(weight[0]),
                              circ_radius, ncells,
                              &x[0], &dx[0],
                              &y[0], &dy[0]);
#else
   KDTree_QueryCircleIntersectWeighted_Float(&mesh->tree, &nez, &(ind[0]), &(weight[0]),
                              circ_radius, ncells,
                              &x[0], &dx[0],
                              &y[0], &dy[0]);
#endif

   for (int ic = 0; ic < nez; ++ic)
   {  H[ind[ic]] = background + (fill_value - background) * weight[ic]; }

   KDTree_Destroy(&mesh->tree);
}

void State::state_reorder(vector<int> iorder)
{
   H = state_memory.memory_reorder(H, &iorder[0]);
   U = state_memory.memory_reorder(U, &iorder[0]);
   V = state_memory.memory_reorder(V, &iorder[0]);
   //printf("\nDEBUG reorder cells\n"); 
   //state_memory.memory_report();
   //printf("DEBUG end reorder cells\n\n"); 
}

void State::rezone_all(int icount, int jcount, vector<char_t> mpot)
{
   struct timespec tstart_cpu;
   cpu_timer_start(&tstart_cpu);

   mesh->rezone_all(icount, jcount, mpot, 1, state_memory);

#ifdef _OPENMP
#pragma omp master
   {
#endif
   memory_reset_ptrs();

   cpu_timers[STATE_TIMER_REZONE_ALL] += cpu_timer_stop(tstart_cpu);
#ifdef _OPENMP
   } // end master region
#endif
}


#ifdef HAVE_OPENCL
void State::gpu_rezone_all(int icount, int jcount, bool localStencil)
{
   struct timespec tstart_cpu;
   cpu_timer_start(&tstart_cpu);

   // Just to get rid of compiler warnings
   //if (1 == 2) printf("DEBUG -- localStencil is %d\n",localStencil);

   mesh->gpu_rezone_all(icount, jcount, dev_mpot, gpu_state_memory);
   dev_H = (cl_mem)gpu_state_memory.get_memory_ptr("dev_H");
   dev_U = (cl_mem)gpu_state_memory.get_memory_ptr("dev_U");
   dev_V = (cl_mem)gpu_state_memory.get_memory_ptr("dev_V");

   gpu_timers[STATE_TIMER_REZONE_ALL] += (long)(cpu_timer_stop(tstart_cpu)*1.0e9);
}
#endif

//define macro for squaring a number
#define SQ(x) ((x)*(x))
//define macro to find minimum of 3 values
//#define MIN3(a,b,c) (min(min((a),(b)),(c)))

#define HXFLUX(ic)  ( U[ic] )
#define UXFLUX(ic)  ( SQ(U[ic])/H[ic] + ghalf*SQ(H[ic]) )
#define UVFLUX(ic)  ( U[ic]*V[ic]/H[ic] )

#define HXFLUXIC ( Uic )
#define HXFLUXNL ( Ul )
#define HXFLUXNR ( Ur )
#define HXFLUXNB ( Ub )
#define HXFLUXNT ( Ut )

#define UXFLUXIC ( SQ(Uic)/Hic + ghalf*SQ(Hic) )
#define UXFLUXNL ( SQ(Ul)/Hl + ghalf*SQ(Hl) )
#define UXFLUXNR ( SQ(Ur)/Hr + ghalf*SQ(Hr) )
#define UXFLUXNB ( SQ(Ub)/Hb + ghalf*SQ(Hb) )
#define UXFLUXNT ( SQ(Ut)/Ht + ghalf*SQ(Ht) )

#define UVFLUXIC ( Uic*Vic/Hic )
#define UVFLUXNL ( Ul*Vl/Hl )
#define UVFLUXNR ( Ur*Vr/Hr )
#define UVFLUXNB ( Ub*Vb/Hb )
#define UVFLUXNT ( Ut*Vt/Ht )

#define HYFLUX(ic)  ( V[ic] )
#define VUFLUX(ic)  ( V[ic]*U[ic]/H[ic] )
#define VYFLUX(ic)  ( SQ(V[ic])/H[ic] + ghalf*SQ(H[ic]) )

#define HYFLUXIC ( Vic )
#define HYFLUXNL ( Vl )
#define HYFLUXNR ( Vr )
#define HYFLUXNB ( Vb )
#define HYFLUXNT ( Vt )

#define VUFLUXIC  ( Vic*Uic/Hic )
#define VUFLUXNL  ( Vl*Ul/Hl )
#define VUFLUXNR  ( Vr*Ur/Hr )
#define VUFLUXNB  ( Vb*Ub/Hb )
#define VUFLUXNT  ( Vt*Ut/Ht )

#define VYFLUXIC  ( SQ(Vic)/Hic + ghalf*SQ(Hic) )
#define VYFLUXNL  ( SQ(Vl)/Hl + ghalf*SQ(Hl) )
#define VYFLUXNR  ( SQ(Vr)/Hr + ghalf*SQ(Hr) )
#define VYFLUXNB  ( SQ(Vb)/Hb + ghalf*SQ(Hb) )
#define VYFLUXNT  ( SQ(Vt)/Ht + ghalf*SQ(Ht) )


#define HNEWXFLUXMINUS  ( Uxminus )
#define HNEWXFLUXPLUS   ( Uxplus )
#define UNEWXFLUXMINUS  ( SQ(Uxminus)/Hxminus + ghalf*SQ(Hxminus) )
#define UNEWXFLUXPLUS   ( SQ(Uxplus) /Hxplus +  ghalf*SQ(Hxplus)  )
#define UVNEWFLUXMINUS  ( Uxminus*Vxminus/Hxminus )
#define UVNEWFLUXPLUS   ( Uxplus *Vxplus /Hxplus  )

#define HNEWYFLUXMINUS  ( Vyminus )
#define HNEWYFLUXPLUS   ( Vyplus  )
#define VNEWYFLUXMINUS  ( SQ(Vyminus)/Hyminus + ghalf*SQ(Hyminus) )
#define VNEWYFLUXPLUS   ( SQ(Vyplus) /Hyplus  + ghalf*SQ(Hyplus)  )
#define VUNEWFLUXMINUS  ( Vyminus*Uyminus/Hyminus )
#define VUNEWFLUXPLUS   ( Vyplus *Uyplus /Hyplus )

#define HXFLUXFACE (Ux)
#define UXFLUXFACE (SQ(Ux)/Hx + ghalf*SQ(Hx))
#define VXFLUXFACE (Ux*Vx/Hx)

#define HYFLUXFACE (Vy)
#define UYFLUXFACE (Vy*Uy/Hy)
#define VYFLUXFACE (SQ(Vy)/Hy + ghalf*SQ(Hy))

// XXX ADDED XXX
#define HXFLUXNLT ( Ult )
#define HXFLUXNRT ( Urt )
#define UXFLUXNLT ( SQR(Ult)/Hlt + ghalf*SQR(Hlt) )
#define UXFLUXNRT ( SQR(Urt)/Hrt + ghalf*SQR(Hrt) )
#define UVFLUXNLT ( Ult*Vlt/Hlt )
#define UVFLUXNRT ( Urt*Vrt/Hrt )
#define HYFLUXNBR ( Vbr )
#define HYFLUXNTR ( Vtr )
#define VUFLUXNBR  ( Vbr*Ubr/Hbr )
#define VUFLUXNTR  ( Vtr*Utr/Htr )
#define VYFLUXNBR  ( SQR(Vbr)/Hbr + ghalf*SQR(Hbr) )
#define VYFLUXNTR  ( SQR(Vtr)/Htr + ghalf*SQR(Htr) )
#define HNEWXFLUXMINUS2  ( Uxminus2 )
#define HNEWXFLUXPLUS2   ( Uxplus2 )
#define UNEWXFLUXMINUS2  ( SQR(Uxminus2)/Hxminus2 + ghalf*SQR(Hxminus2) )
#define UNEWXFLUXPLUS2   ( SQR(Uxplus2) /Hxplus2 +  ghalf*SQR(Hxplus2)  )
#define UVNEWFLUXMINUS2  ( Uxminus2*Vxminus2/Hxminus2 )
#define UVNEWFLUXPLUS2   ( Uxplus2 *Vxplus2 /Hxplus2  )
#define HNEWYFLUXMINUS2  ( Vyminus2 )
#define HNEWYFLUXPLUS2   ( Vyplus2  )
#define VNEWYFLUXMINUS2  ( SQR(Vyminus2)/Hyminus2 + ghalf*SQR(Hyminus2) )
#define VNEWYFLUXPLUS2   ( SQR(Vyplus2) /Hyplus2  + ghalf*SQR(Hyplus2)  )
#define VUNEWFLUXMINUS2  ( Vyminus2*Uyminus2/Hyminus2 )
#define VUNEWFLUXPLUS2   ( Vyplus2 *Uyplus2 /Hyplus2 )

void State::calc_finite_difference(double deltaT)
{
   real_t   g     = 9.80;   // gravitational constant
   real_t   ghalf = 0.5*g;

   //printf("\nDEBUG finite diff\n"); 

   // We need to populate the ghost regions since the calc neighbors has just been
   // established for the mesh shortly before
   apply_boundary_conditions();

   static state_t *H_new, *U_new, *V_new;
#ifdef PRECISION_CHECK_GRAPHICS
   static state_t *PCHECK_new;
#endif
   int *nlft, *nrht, *nbot, *ntop;
   uchar_t *level;

   nlft  = mesh->nlft;
   nrht  = mesh->nrht;
   nbot  = mesh->nbot;
   ntop  = mesh->ntop;
   level = mesh->level;

   vector<real_t> &lev_deltax = mesh->lev_deltax;
   vector<real_t> &lev_deltay = mesh->lev_deltay;

   int flags = 0;
   flags = (RESTART_DATA | REZONE_DATA | LOAD_BALANCE_MEMORY);

#ifdef _OPENMP
#pragma omp master
#endif
   {
      H_new = (state_t *)state_memory.memory_malloc(mesh->ncells_ghost, sizeof(state_t), "H_new", flags);
      U_new = (state_t *)state_memory.memory_malloc(mesh->ncells_ghost, sizeof(state_t), "U_new", flags);
      V_new = (state_t *)state_memory.memory_malloc(mesh->ncells_ghost, sizeof(state_t), "V_new", flags);
#ifdef PRECISION_CHECK_GRAPHICS
      PCHECK_new = (state_t *)state_memory.memory_malloc(mesh->ncells_ghost, sizeof(state_t), "PCHECK_new", flags);
#endif
   }
#ifdef _OPENMP
#pragma omp barrier
#endif

   struct timespec tstart_cpu;
   cpu_timer_start(&tstart_cpu);

   int lowerBound, upperBound;
   mesh->get_bounds(lowerBound, upperBound);

#ifdef PRECISION_CHECK_STATS
   fail_prec_count = 0;
   fail_F_plus_sum    = 0.0;
   fail_F_minus_sum   = 0.0;
   fail_G_plus_sum    = 0.0;
   fail_G_minus_sum   = 0.0;
   fail_wminusx_H_sum = 0.0;
   fail_wplusx_H_sum  = 0.0;
   fail_wminusy_H_sum = 0.0;
   fail_wplusy_H_sum  = 0.0;
   prec_count = 0;
   F_plus_sum    = 0.0;
   F_minus_sum   = 0.0;
   G_plus_sum    = 0.0;
   G_minus_sum   = 0.0;
   wminusx_H_sum = 0.0;
   wplusx_H_sum  = 0.0;
   wminusy_H_sum = 0.0;
   wplusy_H_sum  = 0.0;
#endif

#ifdef _OPENMP_SIMD
#pragma omp simd
#endif
   for(int gix = lowerBound; gix < upperBound; gix++) {

      uchar_t lvl     = level[gix];
      int nl      = nlft[gix];
      int nr      = nrht[gix];
      int nt      = ntop[gix];
      int nb      = nbot[gix];

      real_t Hic     = H[gix];
      real_t Uic     = U[gix];
      real_t Vic     = V[gix];

      int nll     = nlft[nl];
      real_t Hl      = H[nl];
      real_t Ul      = U[nl];
      real_t Vl      = V[nl];

      int nrr     = nrht[nr];
      real_t Hr      = H[nr];
      real_t Ur      = U[nr];
      real_t Vr      = V[nr];

      int ntt     = ntop[nt];
      real_t Ht      = H[nt];
      real_t Ut      = U[nt];
      real_t Vt      = V[nt];

      int nbb     = nbot[nb];
      real_t Hb      = H[nb];
      real_t Ub      = U[nb];
      real_t Vb      = V[nb];

      int nlt     = ntop[nl];
      int nrt     = ntop[nr];
      int ntr     = nrht[nt];
      int nbr     = nrht[nb];

      real_t Hll     = H[nll];
      real_t Ull     = U[nll];
      //real_t Vll     = V[nll];

      real_t Hrr     = H[nrr];
      real_t Urr     = U[nrr];
      //real_t Vrr     = V[nrr];

      real_t Htt     = H[ntt];
      //real_t Utt     = U[ntt];
      real_t Vtt     = V[ntt];

      real_t Hbb     = H[nbb];
      //real_t Ubb     = U[nbb];
      real_t Vbb     = V[nbb];

      real_t dxic    = lev_deltax[lvl];
      real_t dyic    = lev_deltay[lvl];

      real_t dxl     = lev_deltax[level[nl]];
      real_t dxr     = lev_deltax[level[nr]];

      real_t dyt     = lev_deltay[level[nt]];
      real_t dyb     = lev_deltay[level[nb]];

      real_t drl     = dxl;
      real_t drr     = dxr;
      real_t drt     = dyt;
      real_t drb     = dyb;

      real_t dric    = dxic;

      int nltl = 0;
      real_t Hlt = 0.0, Ult = 0.0, Vlt = 0.0;
      real_t Hll2 = 0.0;
      real_t Ull2 = 0.0;
      if(lvl < level[nl]) {
         Hlt  = H[ ntop[nl] ];
         Ult  = U[ ntop[nl] ];
         Vlt  = V[ ntop[nl] ];
         nltl = nlft[nlt];
         Hll2 = H[nltl];
         Ull2 = U[nltl];
      }

      int nrtr = 0;
      real_t Hrt = 0.0, Urt = 0.0, Vrt = 0.0;
      real_t Hrr2 = 0.0;
      real_t Urr2 = 0.0;
      if(lvl < level[nr]) {
         Hrt  = H[ ntop[nr] ];
         Urt  = U[ ntop[nr] ];
         Vrt  = V[ ntop[nr] ];
         nrtr = nrht[nrt];
         Hrr2 = H[nrtr];
         Urr2 = U[nrtr];
      }

      int nbrb = 0;
      real_t Hbr = 0.0, Ubr = 0.0, Vbr = 0.0;
      real_t Hbb2 = 0.0;
      real_t Vbb2 = 0.0;
      if(lvl < level[nb]) {
         Hbr  = H[ nrht[nb] ];
         Ubr  = U[ nrht[nb] ];
         Vbr  = V[ nrht[nb] ];
         nbrb = nbot[nbr];
         Hbb2 = H[nbrb];
         Vbb2 = V[nbrb];
      }

      int ntrt = 0;
      real_t Htr = 0.0, Utr = 0.0, Vtr = 0.0;
      real_t Htt2 = 0.0;
      real_t Vtt2 = 0.0;
      if(lvl < level[nt]) {
         Htr  = H[ nrht[nt] ];
         Utr  = U[ nrht[nt] ];
         Vtr  = V[ nrht[nt] ];
         ntrt = ntop[ntr];
         Htt2 = H[ntrt];
         Vtt2 = V[ntrt];
      }


      real_t Hxminus = U_halfstep(deltaT, Hl, Hic, HXFLUXNL, HXFLUXIC,
                           dxl, dxic, dxl, dxic, SQR(dxl), SQR(dxic));
      real_t Uxminus = U_halfstep(deltaT, Ul, Uic, UXFLUXNL, UXFLUXIC,
                           dxl, dxic, dxl, dxic, SQR(dxl), SQR(dxic));
      real_t Vxminus = U_halfstep(deltaT, Vl, Vic, UVFLUXNL, UVFLUXIC,
                           dxl, dxic, dxl, dxic, SQR(dxl), SQR(dxic));

      real_t Hxplus  = U_halfstep(deltaT, Hic, Hr, HXFLUXIC, HXFLUXNR,
                           dxic, dxr, dxic, dxr, SQR(dxic), SQR(dxr));
      real_t Uxplus  = U_halfstep(deltaT, Uic, Ur, UXFLUXIC, UXFLUXNR,
                           dxic, dxr, dxic, dxr, SQR(dxic), SQR(dxr));
      real_t Vxplus  = U_halfstep(deltaT, Vic, Vr, UVFLUXIC, UVFLUXNR,
                           dxic, dxr, dxic, dxr, SQR(dxic), SQR(dxr));

      real_t Hyminus = U_halfstep(deltaT, Hb, Hic, HYFLUXNB, HYFLUXIC,
                           dyb, dyic, dyb, dyic, SQR(dyb), SQR(dyic));
      real_t Uyminus = U_halfstep(deltaT, Ub, Uic, VUFLUXNB, VUFLUXIC,
                           dyb, dyic, dyb, dyic, SQR(dyb), SQR(dyic));
      real_t Vyminus = U_halfstep(deltaT, Vb, Vic, VYFLUXNB, VYFLUXIC,
                           dyb, dyic, dyb, dyic, SQR(dyb), SQR(dyic));

      real_t Hyplus  = U_halfstep(deltaT, Hic, Ht, HYFLUXIC, HYFLUXNT,
                           dyic, dyt, dyic, dyt, SQR(dyic), SQR(dyt));
      real_t Uyplus  = U_halfstep(deltaT, Uic, Ut, VUFLUXIC, VUFLUXNT,
                           dyic, dyt, dyic, dyt, SQR(dyic), SQR(dyt));
      real_t Vyplus  = U_halfstep(deltaT, Vic, Vt, VYFLUXIC, VYFLUXNT,
                           dyic, dyt, dyic, dyt, SQR(dyic), SQR(dyt));

      real_t Hxfluxminus = HNEWXFLUXMINUS;
      real_t Uxfluxminus = UNEWXFLUXMINUS;
      real_t Vxfluxminus = UVNEWFLUXMINUS;

      real_t Hxfluxplus  = HNEWXFLUXPLUS;
      real_t Uxfluxplus  = UNEWXFLUXPLUS;
      real_t Vxfluxplus  = UVNEWFLUXPLUS;

      real_t Hyfluxminus = HNEWYFLUXMINUS;
      real_t Uyfluxminus = VUNEWFLUXMINUS;
      real_t Vyfluxminus = VNEWYFLUXMINUS;

      real_t Hyfluxplus  = HNEWYFLUXPLUS;
      real_t Uyfluxplus  = VUNEWFLUXPLUS;
      real_t Vyfluxplus  = VNEWYFLUXPLUS;

      real_t Hxminus2 = 0.0;
      real_t Uxminus2 = 0.0;
      real_t Vxminus2 = 0.0;
      if(lvl < level[nl]) {

         Hxminus2 = U_halfstep(deltaT, Hlt, Hic, HXFLUXNLT, HXFLUXIC,
                               drl, dric, drl, dric, SQR(drl), SQR(dric));
         Uxminus2 = U_halfstep(deltaT, Ult, Uic, UXFLUXNLT, UXFLUXIC,
                               drl, dric, drl, dric, SQR(drl), SQR(dric));
         Vxminus2 = U_halfstep(deltaT, Vlt, Vic, UVFLUXNLT, UVFLUXIC,
                               drl, dric, drl, dric, SQR(drl), SQR(dric));

         Hxfluxminus = (Hxfluxminus + HNEWXFLUXMINUS2) * HALF;
         Uxfluxminus = (Uxfluxminus + UNEWXFLUXMINUS2) * HALF;
         Vxfluxminus = (Vxfluxminus + UVNEWFLUXMINUS2) * HALF;

      }

      real_t Hxplus2 = 0.0;
      real_t Uxplus2 = 0.0;
      real_t Vxplus2 = 0.0;
      if(lvl < level[nr]) {

         Hxplus2  = U_halfstep(deltaT, Hic, Hrt, HXFLUXIC, HXFLUXNRT,
                               dric, drr, dric, drr, SQR(dric), SQR(drr));
         Uxplus2  = U_halfstep(deltaT, Uic, Urt, UXFLUXIC, UXFLUXNRT,
                               dric, drr, dric, drr, SQR(dric), SQR(drr));
         Vxplus2  = U_halfstep(deltaT, Vic, Vrt, UVFLUXIC, UVFLUXNRT,
                               dric, drr, dric, drr, SQR(dric), SQR(drr));

         Hxfluxplus  = (Hxfluxplus + HNEWXFLUXPLUS2) * HALF;
         Uxfluxplus  = (Uxfluxplus + UNEWXFLUXPLUS2) * HALF;
         Vxfluxplus  = (Vxfluxplus + UVNEWFLUXPLUS2) * HALF;

      }

      real_t Hyminus2 = 0.0;
      real_t Uyminus2 = 0.0;
      real_t Vyminus2 = 0.0;
      if(lvl < level[nb]) {

         Hyminus2 = U_halfstep(deltaT, Hbr, Hic, HYFLUXNBR, HYFLUXIC,
                               drb, dric, drb, dric, SQR(drb), SQR(dric));
         Uyminus2 = U_halfstep(deltaT, Ubr, Uic, VUFLUXNBR, VUFLUXIC,
                               drb, dric, drb, dric, SQR(drb), SQR(dric));
         Vyminus2 = U_halfstep(deltaT, Vbr, Vic, VYFLUXNBR, VYFLUXIC,
                               drb, dric, drb, dric, SQR(drb), SQR(dric));

         Hyfluxminus = (Hyfluxminus + HNEWYFLUXMINUS2) * HALF;
         Uyfluxminus = (Uyfluxminus + VUNEWFLUXMINUS2) * HALF;
         Vyfluxminus = (Vyfluxminus + VNEWYFLUXMINUS2) * HALF;

      }

      real_t Hyplus2 = 0.0;
      real_t Uyplus2 = 0.0;
      real_t Vyplus2 = 0.0;
      if(lvl < level[nt]) {

         Hyplus2  = U_halfstep(deltaT, Hic, Htr, HYFLUXIC, HYFLUXNTR,
                               dric, drt, dric, drt, SQR(dric), SQR(drt));
         Uyplus2  = U_halfstep(deltaT, Uic, Utr, VUFLUXIC, VUFLUXNTR,
                               dric, drt, dric, drt, SQR(dric), SQR(drt));
         Vyplus2  = U_halfstep(deltaT, Vic, Vtr, VYFLUXIC, VYFLUXNTR,
                               dric, drt, dric, drt, SQR(dric), SQR(drt));

         Hyfluxplus  = (Hyfluxplus + HNEWYFLUXPLUS2) * HALF;
         Uyfluxplus  = (Uyfluxplus + VUNEWFLUXPLUS2) * HALF;
         Vyfluxplus  = (Vyfluxplus + VNEWYFLUXPLUS2) * HALF;

      }

      ////////////////////////////////////////
      /// Artificial Viscosity corrections ///
      ////////////////////////////////////////


      if(level[nl] < level[nll]) {
         Hll = (Hll + H[ ntop[nll] ]) * HALF;
         Ull = (Ull + U[ ntop[nll] ]) * HALF;
      }

      real_t Hr2 = Hr;
      real_t Ur2 = Ur;
      if(lvl < level[nr]) {
         Hr2 = (Hr2 + Hrt) * HALF;
         Ur2 = (Ur2 + Urt) * HALF;
      }

      real_t wminusx_H = w_corrector(deltaT, (dric+dxl)*HALF, fabs(Uxminus/Hxminus) + sqrt(g*Hxminus),
                              Hic-Hl, Hl-Hll, Hr2-Hic);

      wminusx_H *= Hic - Hl;

      if(lvl < level[nl]) {
         if(level[nlt] < level[nltl])
            Hll2 = (Hll2 + H[ ntop[nltl] ]) * HALF;
         wminusx_H = ((w_corrector(deltaT, (dric+dxl)*HALF, fabs(Uxminus2/Hxminus2) +
                                  sqrt(g*Hxminus2), Hic-Hlt, Hlt-Hll2, Hr2-Hic) *
                      (Hic - Hlt)) + wminusx_H)*HALF*HALF;
      }


      if(level[nr] < level[nrr]) {
         Hrr = (Hrr + H[ ntop[nrr] ]) * HALF;
         Urr = (Urr + U[ ntop[nrr] ]) * HALF;
      }

      real_t Hl2 = Hl;
      real_t Ul2 = Ul;
      if(lvl < level[nl]) {
         Hl2 = (Hl2 + Hlt) * HALF;
         Ul2 = (Ul2 + Ult) * HALF;
      }

      real_t wplusx_H = w_corrector(deltaT, (dric+dxr)*HALF, fabs(Uxplus/Hxplus) + sqrt(g*Hxplus),
                           Hr-Hic, Hic-Hl2, Hrr-Hr);

      wplusx_H *= Hr - Hic;

      if(lvl < level[nr]) {
         if(level[nrt] < level[nrtr])
            Hrr2 = (Hrr2 + H[ ntop[nrtr] ]) * HALF;
         wplusx_H = ((w_corrector(deltaT, (dric+dxr)*HALF, fabs(Uxplus2/Hxplus2) +
                                  sqrt(g*Hxplus2), Hrt-Hic, Hic-Hl2, Hrr2-Hrt) *
                      (Hrt - Hic))+wplusx_H)*HALF*HALF;
      }


      real_t wminusx_U = w_corrector(deltaT, (dric+dxl)*HALF, fabs(Uxminus/Hxminus) + sqrt(g*Hxminus),
                              Uic-Ul, Ul-Ull, Ur2-Uic);

      wminusx_U *= Uic - Ul;

      if(lvl < level[nl]) {
         if(level[nlt] < level[nltl])
            Ull2 = (Ull2 + U[ ntop[nltl] ]) * HALF;
         wminusx_U = ((w_corrector(deltaT, (dric+dxl)*HALF, fabs(Uxminus2/Hxminus2) +
                                  sqrt(g*Hxminus2), Uic-Ult, Ult-Ull2, Ur2-Uic) *
                      (Uic - Ult))+wminusx_U)*HALF*HALF;
      }


      real_t wplusx_U = w_corrector(deltaT, (dric+dxr)*HALF, fabs(Uxplus/Hxplus) + sqrt(g*Hxplus),
                              Ur-Uic, Uic-Ul2, Urr-Ur);

      wplusx_U *= Ur - Uic;

      if(lvl < level[nr]) {
         if(level[nrt] < level[nrtr])
            Urr2 = (Urr2 + U[ ntop[nrtr] ]) * HALF;
         wplusx_U = ((w_corrector(deltaT, (dric+dxr)*HALF, fabs(Uxplus2/Hxplus2) +
                                  sqrt(g*Hxplus2), Urt-Uic, Uic-Ul2, Urr2-Urt) *
                      (Urt - Uic))+wplusx_U)*HALF*HALF;
      }


      if(level[nb] < level[nbb]) {
         Hbb = (Hbb + H[ nrht[nbb] ]) * HALF;
         Vbb = (Vbb + V[ nrht[nbb] ]) * HALF;
      }

      real_t Ht2 = Ht;
      real_t Vt2 = Vt;
      if(lvl < level[nt]) {
         Ht2 = (Ht2 + Htr) * HALF;
         Vt2 = (Vt2 + Vtr) * HALF;
      }

      real_t wminusy_H = w_corrector(deltaT, (dric+dyb)*HALF, fabs(Vyminus/Hyminus) + sqrt(g*Hyminus),
                              Hic-Hb, Hb-Hbb, Ht2-Hic);

      wminusy_H *= Hic - Hb;

      if(lvl < level[nb]) {
         if(level[nbr] < level[nbrb])
            Hbb2 = (Hbb2 + H[ nrht[nbrb] ]) * HALF;
         wminusy_H = ((w_corrector(deltaT, (dric+dyb)*HALF, fabs(Vyminus2/Hyminus2) +
                                  sqrt(g*Hyminus2), Hic-Hbr, Hbr-Hbb2, Ht2-Hic) *
                      (Hic - Hbr))+wminusy_H)*HALF*HALF;
      }


      if(level[nt] < level[ntt]) {
         Htt = (Htt + H[ nrht[ntt] ]) * HALF;
         Vtt = (Vtt + V[ nrht[ntt] ]) * HALF;
      }

      real_t Hb2 = Hb;
      real_t Vb2 = Vb;
      if(lvl < level[nb]) {
         Hb2 = (Hb2 + Hbr) * HALF;
         Vb2 = (Vb2 + Vbr) * HALF;
      }

      real_t wplusy_H = w_corrector(deltaT, (dric+dyt)*HALF, fabs(Vyplus/Hyplus) + sqrt(g*Hyplus),
                             Ht-Hic, Hic-Hb2, Htt-Ht);

      wplusy_H *= Ht - Hic;

      if(lvl < level[nt]) {
         if(level[ntr] < level[ntrt])
            Htt2 = (Htt2 + H[ nrht[ntrt] ]) * HALF;
         wplusy_H = ((w_corrector(deltaT, (dric+dyt)*HALF, fabs(Vyplus2/Hyplus2) +
                                  sqrt(g*Hyplus2), Htr-Hic, Hic-Hb2, Htt2-Htr) *
                      (Htr - Hic))+wplusy_H)*HALF*HALF;
      }

      real_t wminusy_V = w_corrector(deltaT, (dric+dyb)*HALF, fabs(Vyminus/Hyminus) + sqrt(g*Hyminus),
                              Vic-Vb, Vb-Vbb, Vt2-Vic);

      wminusy_V *= Vic - Vb;

      if(lvl < level[nb]) {
         if(level[nbr] < level[nbrb])
            Vbb2 = (Vbb2 + V[ nrht[nbrb] ]) * HALF;
         wminusy_V = ((w_corrector(deltaT, (dric+dyb)*HALF, fabs(Vyminus2/Hyminus2) +
                                  sqrt(g*Hyminus2), Vic-Vbr, Vbr-Vbb2, Vt2-Vic) *
                      (Vic - Vbr))+wminusy_V)*HALF*HALF;
      }

      real_t wplusy_V = w_corrector(deltaT, (dric+dyt)*HALF, fabs(Vyplus/Hyplus) + sqrt(g*Hyplus),
                           Vt-Vic, Vic-Vb2, Vtt-Vt);

      wplusy_V *= Vt - Vic;

      if(lvl < level[nt]) {
         if(level[ntr] < level[ntrt])
            Vtt2 = (Vtt2 + V[ nrht[ntrt] ]) * HALF;
         wplusy_V = ((w_corrector(deltaT, (dric+dyt)*HALF, fabs(Vyplus2/Hyplus2) +
                                  sqrt(g*Hyplus2), Vtr-Vic, Vic-Vb2, Vtt2-Vtr) *
                      (Vtr - Vic))+wplusy_V)*HALF*HALF;
      }

#ifdef PRECISION_CHECK_WITH_PARENTHESIS
      //Basic parentheses
      H_new[gix] = U_fullstep(deltaT, dxic, Hic,
                       Hxfluxplus, Hxfluxminus, Hyfluxplus, Hyfluxminus)
                  +( - wminusx_H + wplusx_H - wminusy_H + wplusy_H );
#else
#ifdef PRECISION_CHECK_BEST_PARENTHESIS
      //Version with "best" parentheses
      H_new[gix] = U_fullstep(deltaT, dxic, Hic,
                       Hxfluxplus, Hxfluxminus, Hyfluxplus, Hyfluxminus)
                   //+ (((-wminusx_H - wminusy_H) + wplusy_H) + wplusx_H);
                   + (((wplusx_H + wplusy_H) - wminusy_H) - wminusx_H);
#else
      //Original version - no parentheses
      H_new[gix] = U_fullstep(deltaT, dxic, Hic,
                       Hxfluxplus, Hxfluxminus, Hyfluxplus, Hyfluxminus)
                  - wminusx_H + wplusx_H - wminusy_H + wplusy_H;
#endif
#endif

#ifdef PRECISION_CHECK
      int fail;
      U_fullstep_precision_check(gix, deltaT, dxic, Hic, H_new[gix],
                      Hxfluxplus, Hxfluxminus, Hyfluxplus, Hyfluxminus,
                      wminusx_H, wplusx_H, wminusy_H, wplusy_H, &fail); 


#ifdef PRECISION_CHECK_GRAPHICS
      PCHECK_new[gix] = 100.0*(double)fail;
#endif
#endif

      U_new[gix] = U_fullstep(deltaT, dxic, Uic,
                       Uxfluxplus, Uxfluxminus, Uyfluxplus, Uyfluxminus)
                  - wminusx_U + wplusx_U;
      V_new[gix] = U_fullstep(deltaT, dxic, Vic,
                       Vxfluxplus, Vxfluxminus, Vyfluxplus, Vyfluxminus)
                  - wminusy_V + wplusy_V;

   } // cell loop

#ifdef PRECISION_CHECK_STATS
   fail_F_plus_sum    /= (double)fail_prec_count;
   fail_F_minus_sum   /= (double)fail_prec_count;
   fail_G_plus_sum    /= (double)fail_prec_count;
   fail_G_minus_sum   /= (double)fail_prec_count;
   fail_wminusx_H_sum /= (double)fail_prec_count;
   fail_wplusx_H_sum  /= (double)fail_prec_count;
   fail_wminusy_H_sum /= (double)fail_prec_count;
   fail_wplusy_H_sum  /= (double)fail_prec_count;
   F_plus_sum    /= (double)prec_count;
   F_minus_sum   /= (double)prec_count;
   G_plus_sum    /= (double)prec_count;
   G_minus_sum   /= (double)prec_count;
   wminusx_H_sum /= (double)prec_count;
   wplusx_H_sum  /= (double)prec_count;
   wminusy_H_sum /= (double)prec_count;
   wplusy_H_sum  /= (double)prec_count;

   fail_prec_avg_count++;
   fail_F_plus_avg    += fail_F_plus_sum;
   fail_F_minus_avg   += fail_F_minus_sum;
   fail_G_plus_avg    += fail_G_plus_sum;
   fail_G_minus_avg   += fail_G_minus_sum;
   fail_wminusx_H_avg += fail_wminusx_H_sum;
   fail_wplusx_H_avg  += fail_wplusx_H_sum;
   fail_wminusy_H_avg += fail_wminusy_H_sum;
   fail_wplusy_H_avg  += fail_wplusy_H_sum;
   prec_avg_count++;
   F_plus_avg    += F_plus_sum;
   F_minus_avg   += F_minus_sum;
   G_plus_avg    += G_plus_sum;
   G_minus_avg   += G_minus_sum;
   wminusx_H_avg += wminusx_H_sum;
   wplusx_H_avg  += wplusx_H_sum;
   wminusy_H_avg += wminusy_H_sum;
   wplusy_H_avg  += wplusy_H_sum;
#endif


#ifdef _OPENMP
#pragma omp barrier
#pragma omp master
   {
#endif
      // Replace H with H_new and deallocate H. New memory will have the characteristics
      // of the new memory and the name of the old. Both return and arg1 will be reset to new memory
      H = (state_t *)state_memory.memory_replace(H, H_new);
      U = (state_t *)state_memory.memory_replace(U, U_new);
      V = (state_t *)state_memory.memory_replace(V, V_new);
#ifdef PRECISION_CHECK_GRAPHICS
      PCHECK = (state_t *)state_memory.memory_replace(PCHECK, PCHECK_new);
#endif

      //state_memory.memory_report();
      //printf("DEBUG end finite diff\n\n"); 

      cpu_timers[STATE_TIMER_FINITE_DIFFERENCE] += cpu_timer_stop(tstart_cpu);
#ifdef _OPENMP
   }
#pragma omp barrier
#endif
}

void State::calc_finite_difference_cell_in_place(double deltaT)
{
   real_t   g     = 9.80;   // gravitational constant
   real_t   ghalf = HALF*g;

   struct timespec tstart_cpu;
   cpu_timer_start(&tstart_cpu);

   struct timespec tstart_cpu_part;
   cpu_timer_start(&tstart_cpu_part);

   // We need to populate the ghost regions since the calc neighbors has just been
   // established for the mesh shortly before
   apply_boundary_conditions();

   static state_t *H_new, *U_new, *V_new;
#ifdef PRECISION_CHECK_GRAPHICS
   static state_t *PCHECK_new;
#endif
   static real_t *Hxfluxminus, *Uxfluxminus, *Vxfluxminus, *Hxfluxplus, *Uxfluxplus, *Vxfluxplus;
   static real_t *Hyfluxminus, *Uyfluxminus, *Vyfluxminus, *Hyfluxplus, *Uyfluxplus, *Vyfluxplus;
   static real_t *wminusx_H, *wminusx_U, *wplusx_H, *wplusx_U, *wminusy_H, *wminusy_V, *wplusy_H, *wplusy_V;

#ifdef _OPENMP
#pragma omp barrier
#pragma omp master
   {
#endif
      mesh->calc_face_list_wbidirmap_phantom(state_memory, deltaT);
      memory_reset_ptrs(); //reset the pointers H,U,V that were recently reallocated in wbidirmap call

      int flags = (RESTART_DATA | REZONE_DATA | LOAD_BALANCE_MEMORY);

      H_new = (state_t *)state_memory.memory_malloc(mesh->ncells_ghost, sizeof(state_t), "H_new", flags);
      U_new = (state_t *)state_memory.memory_malloc(mesh->ncells_ghost, sizeof(state_t), "U_new", flags);
      V_new = (state_t *)state_memory.memory_malloc(mesh->ncells_ghost, sizeof(state_t), "V_new", flags);
#ifdef PRECISION_CHECK_GRAPHICS
      PCHECK_new = (state_t *)state_memory.memory_malloc(mesh->ncells_ghost, sizeof(state_t), "PCHECK_new", flags);
#endif

      Hxfluxminus = (real_t *)malloc(mesh->ncells*sizeof(real_t));
      Uxfluxminus = (real_t *)malloc(mesh->ncells*sizeof(real_t));
      Vxfluxminus = (real_t *)malloc(mesh->ncells*sizeof(real_t));

      Hxfluxplus  = (real_t *)malloc(mesh->ncells*sizeof(real_t));
      Uxfluxplus  = (real_t *)malloc(mesh->ncells*sizeof(real_t));
      Vxfluxplus  = (real_t *)malloc(mesh->ncells*sizeof(real_t));

      Hyfluxminus = (real_t *)malloc(mesh->ncells*sizeof(real_t));
      Uyfluxminus = (real_t *)malloc(mesh->ncells*sizeof(real_t));
      Vyfluxminus = (real_t *)malloc(mesh->ncells*sizeof(real_t));

      Hyfluxplus  = (real_t *)malloc(mesh->ncells*sizeof(real_t));
      Uyfluxplus  = (real_t *)malloc(mesh->ncells*sizeof(real_t));
      Vyfluxplus  = (real_t *)malloc(mesh->ncells*sizeof(real_t));

      wminusx_H   = (real_t *)malloc(mesh->ncells*sizeof(real_t));
      wminusx_U   = (real_t *)malloc(mesh->ncells*sizeof(real_t));

      wplusx_H    = (real_t *)malloc(mesh->ncells*sizeof(real_t));
      wplusx_U    = (real_t *)malloc(mesh->ncells*sizeof(real_t));

      wminusy_H   = (real_t *)malloc(mesh->ncells*sizeof(real_t));
      wminusy_V   = (real_t *)malloc(mesh->ncells*sizeof(real_t));

      wplusy_H    = (real_t *)malloc(mesh->ncells*sizeof(real_t));
      wplusy_V    = (real_t *)malloc(mesh->ncells*sizeof(real_t));

#ifdef _OPENMP
   }
#pragma omp barrier
#endif

   cpu_timers[STATE_TIMER_FINITE_DIFFERENCE_PART1] += cpu_timer_stop(tstart_cpu_part);
   cpu_timer_start(&tstart_cpu_part);

   int lowerBound, upperBound;

   mesh->get_bounds(lowerBound, upperBound);

#if defined(__GNUC_MINOR__)
   static real_t * dxcell;
   static real_t * dycell;
#ifdef _OPENMP
#pragma omp barrier
#pragma omp master
   {
#endif
   dxcell = (real_t *)malloc(sizeof(real_t) * mesh->ncells);
   dycell = (real_t *)malloc(sizeof(real_t) * mesh->ncells);
#ifdef _OPENMP
   }
#pragma omp barrier
#endif
   for (int ic = lowerBound; ic < upperBound; ic++) {
      uchar_t lev = mesh->level[ic];
      dxcell[ic]  = mesh->lev_deltax[lev];
      dycell[ic]  = mesh->lev_deltay[lev];
   }
   
   state_t *H_loc = H;
   state_t *U_loc = U;
   state_t *V_loc = V;
#endif
   
#ifdef _OPENMP
#ifdef _OPENMP_SIMD
#pragma omp for simd
#else
#pragma omp for
#endif
#else
#ifdef _OPENMP_SIMD
#pragma omp simd
#endif
#endif
   for (int ic = lowerBound; ic < upperBound; ic++) {
   //for (int ic = 0; ic < mesh->ncells; ic++) {
      if (mesh->celltype[ic] != REAL_CELL) continue;
#ifdef _OPENMP
      //printf("%d) %d %d\n", omp_get_thread_num(), lowerBound, upperBound);
#endif

#if defined(__GNUC_MINOR__)
      real_t dxic = dxcell[ic];
      real_t dyic = dycell[ic];
#else
      uchar_t lev = mesh->level[ic];
      real_t dxic    = mesh->lev_deltax[lev];
      real_t dyic    = mesh->lev_deltay[lev];
#endif

      real_t Cxhalf = 0.5*deltaT/dxic;
      real_t Cyhalf = 0.5*deltaT/dyic;

      // set the four faces
      int fl = mesh->map_xcell2face_left1[ic];
      int fr = mesh->map_xcell2face_right1[ic];
      int fb = mesh->map_ycell2face_bot1[ic];
      int ft = mesh->map_ycell2face_top1[ic];

      // set the four neighboring cells
      int nl = mesh->map_xface2cell_lower[fl];
      int nr = mesh->map_xface2cell_upper[fr];
      int nb = mesh->map_yface2cell_lower[fb];
      int nt = mesh->map_yface2cell_upper[ft];
      //printf("%d) %d\n", ic, nl);

      real_t Hic     = H[ic];
      real_t Uic     = U[ic];
      real_t Vic     = V[ic];

      int nll     = mesh->map_xface2cell_lower[mesh->map_xcell2face_left1[nl]];
      real_t Hl      = H[nl];
      real_t Ul      = U[nl];
      real_t Vl      = V[nl];

      int nrr     = mesh->map_xface2cell_upper[mesh->map_xcell2face_right1[nr]];
      real_t Hr      = H[nr];
      real_t Ur      = U[nr];
      real_t Vr      = V[nr];

      int ntt     = mesh->map_yface2cell_upper[mesh->map_ycell2face_top1[nt]];
      real_t Ht      = H[nt];
      real_t Ut      = U[nt];
      real_t Vt      = V[nt];

      int nbb     = mesh->map_yface2cell_lower[mesh->map_ycell2face_bot1[nb]];
      real_t Hb      = H[nb];
      real_t Ub      = U[nb];
      real_t Vb      = V[nb];

      real_t Hll     = H[nll];
      real_t Ull     = U[nll];

      real_t Hrr     = H[nrr];
      real_t Urr     = U[nrr];

      real_t Htt     = H[ntt];
      real_t Vtt     = V[ntt];

      real_t Hbb     = H[nbb];
      real_t Vbb     = V[nbb];

      // Halfstep in space and time
      real_t Hxminus = HALF*(Hic+Hl)-Cxhalf*(HXFLUXIC-HXFLUXNL);
      real_t Uxminus = HALF*(Uic+Ul)-Cxhalf*(UXFLUXIC-UXFLUXNL);
      real_t Vxminus = HALF*(Vic+Vl)-Cxhalf*(UVFLUXIC-UVFLUXNL);

      real_t Hxplus = HALF*(Hr+Hic)-Cxhalf*(HXFLUXNR-HXFLUXIC);
      real_t Uxplus = HALF*(Ur+Uic)-Cxhalf*(UXFLUXNR-UXFLUXIC);
      real_t Vxplus = HALF*(Vr+Vic)-Cxhalf*(UVFLUXNR-UVFLUXIC);

      real_t Hyminus = HALF*(Hic+Hb)-Cyhalf*(HYFLUXIC-HYFLUXNB);
      real_t Uyminus = HALF*(Uic+Ub)-Cyhalf*(VUFLUXIC-VUFLUXNB);
      real_t Vyminus = HALF*(Vic+Vb)-Cyhalf*(VYFLUXIC-VYFLUXNB);

      real_t Hyplus = HALF*(Ht+Hic)-Cyhalf*(HYFLUXNT-HYFLUXIC);
      real_t Uyplus = HALF*(Ut+Uic)-Cyhalf*(VUFLUXNT-VUFLUXIC);
      real_t Vyplus = HALF*(Vt+Vic)-Cyhalf*(VYFLUXNT-VYFLUXIC);

      Hxfluxminus[ic] = HNEWXFLUXMINUS;
      Uxfluxminus[ic] = UNEWXFLUXMINUS;
      Vxfluxminus[ic] = UVNEWFLUXMINUS;

      Hxfluxplus[ic]  = HNEWXFLUXPLUS;
      Uxfluxplus[ic]  = UNEWXFLUXPLUS;
      Vxfluxplus[ic]  = UVNEWFLUXPLUS;

      Hyfluxminus[ic] = HNEWYFLUXMINUS;
      Uyfluxminus[ic] = VUNEWFLUXMINUS;
      Vyfluxminus[ic] = VNEWYFLUXMINUS;

      Hyfluxplus[ic]  = HNEWYFLUXPLUS;
      Uyfluxplus[ic]  = VUNEWFLUXPLUS;
      Vyfluxplus[ic]  = VNEWYFLUXPLUS;

      ////////////////////////////////////////
      /// Artificial Viscosity corrections ///
      ////////////////////////////////////////

      real_t U_eigen = fabs(Uxminus/Hxminus) + sqrt(g*Hxminus);
      wminusx_H[ic] = w_corrector(deltaT, dxic, U_eigen, Hic-Hl, Hl-Hll, Hr-Hic) * (Hic - Hl);
      wminusx_U[ic] = w_corrector(deltaT, dxic, U_eigen, Uic-Ul, Ul-Ull, Ur-Uic) * (Uic - Ul);

      U_eigen = fabs(Uxplus/Hxplus) + sqrt(g*Hxplus);
      wplusx_H[ic] = w_corrector(deltaT, dxic, U_eigen, Hr-Hic, Hic-Hl, Hrr-Hr) * (Hr - Hic);
      wplusx_U[ic] = w_corrector(deltaT, dxic, U_eigen, Ur-Uic, Uic-Ul, Urr-Ur) * (Ur - Uic);

      U_eigen = fabs(Vyminus/Hyminus) + sqrt(g*Hyminus);
      wminusy_H[ic] = w_corrector(deltaT, dyic, U_eigen, Hic-Hb, Hb-Hbb, Ht-Hic) * (Hic - Hb);
      wminusy_V[ic] = w_corrector(deltaT, dyic, U_eigen, Vic-Vb, Vb-Vbb, Vt-Vic) * (Vic - Vb);

      U_eigen = fabs(Vyplus/Hyplus) + sqrt(g*Hyplus);
      wplusy_H[ic] = w_corrector(deltaT, dyic, U_eigen, Ht-Hic, Hic-Hb, Htt-Ht) * (Ht - Hic);
      wplusy_V[ic] = w_corrector(deltaT, dyic, U_eigen, Vt-Vic, Vic-Vb, Vtt-Vt) * (Vt - Vic);
   }

   cpu_timers[STATE_TIMER_FINITE_DIFFERENCE_PART2] += cpu_timer_stop(tstart_cpu_part);
   cpu_timer_start(&tstart_cpu_part);

#ifdef _OPENMP
#ifdef _OPENMP_SIMD
#pragma omp for simd
#else
#pragma omp for
#endif
#else
#ifdef _OPENMP_SIMD
#pragma omp simd
#endif
#endif
   for (int ifix = 0; ifix < mesh->nxfixup; ifix++){
      int ic = mesh->xrecvCIdx[ifix];

      if (mesh->xplusCell2Idx[ic] > -1) {
         int ifixup = mesh->xplusCell2Idx[ic];

         // set the sending cells
         int ns1 = mesh->map_xface2cell_upper[mesh->xsendIdx1[ifixup]];
         int ns2 = mesh->map_xface2cell_upper[mesh->xsendIdx2[ifixup]];

         Hxfluxplus[ic] = (Hxfluxminus[ns1] + Hxfluxminus[ns2]) * HALF;
         Uxfluxplus[ic] = (Uxfluxminus[ns1] + Uxfluxminus[ns2]) * HALF;
         Vxfluxplus[ic] = (Vxfluxminus[ns1] + Vxfluxminus[ns2]) * HALF;
         wplusx_H[ic] = (wminusx_H[ns1] + wminusx_H[ns2]) * 0.25;
         wplusx_U[ic] = (wminusx_U[ns1] + wminusx_U[ns2]) * 0.25;
      }

      if (mesh->xminusCell2Idx[ic] > -1) {
         int ifixup = mesh->xminusCell2Idx[ic];

         // set the sending cells
         int ns1 = mesh->map_xface2cell_lower[mesh->xsendIdx1[ifixup]];
         int ns2 = mesh->map_xface2cell_lower[mesh->xsendIdx2[ifixup]];

         Hxfluxminus[ic] = (Hxfluxplus[ns1] + Hxfluxplus[ns2]) * HALF;
         Uxfluxminus[ic] = (Uxfluxplus[ns1] + Uxfluxplus[ns2]) * HALF;
         Vxfluxminus[ic] = (Vxfluxplus[ns1] + Vxfluxplus[ns2]) * HALF;
         //printf("(%d) %f\n", ic, Vxfluxminus[ic]);
         wminusx_H[ic] = (wplusx_H[ns1] + wplusx_H[ns2]) * 0.25;
         wminusx_U[ic] = (wplusx_U[ns1] + wplusx_U[ns2]) * 0.25;
      }
   }

#ifdef _OPENMP
#ifdef _OPENMP_SIMD
#pragma omp for simd
#else
#pragma omp for
#endif
#else
#ifdef _OPENMP_SIMD
#pragma omp simd
#endif
#endif
   for (int ifix = 0; ifix < mesh->nyfixup; ifix++){
      int ic = mesh->yrecvCIdx[ifix];

      if (mesh->yplusCell2Idx[ic] > -1) {
         int ifixup = mesh->yplusCell2Idx[ic];

         // set the sending cells
         int ns1 = mesh->map_yface2cell_upper[mesh->ysendIdx1[ifixup]];
         int ns2 = mesh->map_yface2cell_upper[mesh->ysendIdx2[ifixup]];

         Hyfluxplus[ic] = (Hyfluxminus[ns1] + Hyfluxminus[ns2]) * HALF;
         Uyfluxplus[ic] = (Uyfluxminus[ns1] + Uyfluxminus[ns2]) * HALF;
         Vyfluxplus[ic] = (Vyfluxminus[ns1] + Vyfluxminus[ns2]) * HALF;
         wplusy_H[ic] = (wminusy_H[ns1] + wminusy_H[ns2]) * 0.25;
         wplusy_V[ic] = (wminusy_V[ns1] + wminusy_V[ns2]) * 0.25;
      }

      if (mesh->yminusCell2Idx[ic] > -1) {
         int ifixup = mesh->yminusCell2Idx[ic];

         // set the sending cells
         int ns1 = mesh->map_yface2cell_lower[mesh->ysendIdx1[ifixup]];
         int ns2 = mesh->map_yface2cell_lower[mesh->ysendIdx2[ifixup]];

         Hyfluxminus[ic] = (Hyfluxplus[ns1] + Hyfluxplus[ns2]) * HALF;
         Uyfluxminus[ic] = (Uyfluxplus[ns1] + Uyfluxplus[ns2]) * HALF;
         Vyfluxminus[ic] = (Vyfluxplus[ns1] + Vyfluxplus[ns2]) * HALF;
         wminusy_H[ic] = (wplusy_H[ns1] + wplusy_H[ns2]) * 0.25;
         wminusy_V[ic] = (wplusy_V[ns1] + wplusy_V[ns2]) * 0.25;
      }
   }

   cpu_timers[STATE_TIMER_FINITE_DIFFERENCE_PART3] += cpu_timer_stop(tstart_cpu_part);
   cpu_timer_start(&tstart_cpu_part);

#if defined(__GNUC_MINOR__)
   H_loc = H;
   U_loc = U;
   V_loc = V;
#endif

#ifdef _OPENMP
#pragma omp barrier
#endif
#ifdef _OPENMP_SIMD
#pragma omp simd
#endif
   for (int ic = lowerBound; ic < upperBound; ic++) {
      if (mesh->celltype[ic] != REAL_CELL) continue;

#if defined(__GNUC_MINOR__)
      real_t dxic = dxcell[ic];
      //real_t dyic = dycell[ic];
#else
      uchar_t lev = mesh->level[ic];
      real_t dxic    = mesh->lev_deltax[lev];
      //real_t dyic    = mesh->lev_deltay[lev];
#endif

      H_new[ic] = U_fullstep(deltaT, dxic, H[ic],
                       Hxfluxplus[ic], Hxfluxminus[ic], Hyfluxplus[ic], Hyfluxminus[ic])
                  - wminusx_H[ic] + wplusx_H[ic] - wminusy_H[ic] + wplusy_H[ic];
      U_new[ic] = U_fullstep(deltaT, dxic, U[ic],
                       Uxfluxplus[ic], Uxfluxminus[ic], Uyfluxplus[ic], Uyfluxminus[ic])
                  - wminusx_U[ic] + wplusx_U[ic];
      V_new[ic] = U_fullstep(deltaT, dxic, V[ic],
                       Vxfluxplus[ic], Vxfluxminus[ic], Vyfluxplus[ic], Vyfluxminus[ic])
                  - wminusy_V[ic] + wplusy_V[ic];
   } // cell loop

   cpu_timers[STATE_TIMER_FINITE_DIFFERENCE_PART4] += cpu_timer_stop(tstart_cpu_part);

#ifdef _OPENMP
#pragma omp barrier
#pragma omp master
   {
#endif
      free(Hxfluxminus);
      free(Uxfluxminus);
      free(Vxfluxminus);

      free(Hxfluxplus);
      free(Uxfluxplus);
      free(Vxfluxplus);

      free(Hyfluxminus);
      free(Uyfluxminus);
      free(Vyfluxminus);

      free(Hyfluxplus);
      free(Uyfluxplus);
      free(Vyfluxplus);

      free(wminusx_H);
      free(wminusx_U);

      free(wplusx_H);
      free(wplusx_U);

      free(wminusy_H);
      free(wminusy_V);

      free(wplusy_H);
      free(wplusy_V);

      // Replace H with H_new and deallocate H. New memory will have the characteristics
      // of the new memory and the name of the old. Both return and arg1 will be reset to new memory
      H = (state_t *)state_memory.memory_replace(H, H_new);
      U = (state_t *)state_memory.memory_replace(U, U_new);
      V = (state_t *)state_memory.memory_replace(V, V_new);
#ifdef PRECISION_CHECK_GRAPHICS
      PCHECK = (state_t *)state_memory.memory_replace(PCHECK, PCHECK_new);
#endif

#if defined(__GNUC_MINOR__)
      free(dxcell);
      free(dycell);
#endif

      cpu_timers[STATE_TIMER_FINITE_DIFFERENCE] += cpu_timer_stop(tstart_cpu);
#ifdef _OPENMP
   }
#pragma omp barrier
#endif
}

void State::calc_finite_difference_face_in_place(double deltaT)
{
   real_t   g     = 9.80;   // gravitational constant
   real_t   ghalf = HALF*g;

   struct timespec tstart_cpu;
   cpu_timer_start(&tstart_cpu);

   struct timespec tstart_cpu_part;
   cpu_timer_start(&tstart_cpu_part);

   // We need to populate the ghost regions since the calc neighbors has just been
   // established for the mesh shortly before
   apply_boundary_conditions();

   static int xfaceSize; //new "update" nxface inc. phantoms
   static int yfaceSize;

   static state_t *HxFlux, *UxFlux, *VxFlux, *Wx_H, *Wx_U;
   static state_t *HyFlux, *UyFlux, *VyFlux, *Wy_H, *Wy_V;


#ifdef _OPENMP
#pragma omp barrier
#pragma omp master
   {
#endif
      mesh->calc_face_list_wbidirmap_phantom(state_memory, deltaT);
      xfaceSize = mesh->pxface; //new "update" nxface inc. phantoms
      yfaceSize = mesh->pyface; //new "update" nyface inc. phantoms
      memory_reset_ptrs(); //reset the pointers H,U,V that were recently reallocated in wbidirmap call

      /*
      Wx_H.clear();
      Wx_H.resize(xfaceSize, (state_t)0.0);
      Wx_U.clear();
      Wx_U.resize(xfaceSize, (state_t)0.0);
      HxFlux.clear();
      HxFlux.resize(xfaceSize, (state_t)0.0);
      UxFlux.clear();
      UxFlux.resize(xfaceSize, (state_t)0.0);
      VxFlux.clear();
      VxFlux.resize(xfaceSize, (state_t)0.0);

      Wy_H.clear();
      Wy_H.resize(yfaceSize, (state_t)0.0);
      Wy_V.clear();
      Wy_V.resize(yfaceSize, (state_t)0.0);
      HyFlux.clear();
      HyFlux.resize(yfaceSize, (state_t)0.0);
      UyFlux.clear();
      UyFlux.resize(yfaceSize, (state_t)0.0);
      VyFlux.clear();
      VyFlux.resize(yfaceSize, (state_t)0.0);
      */
    
      HxFlux = (state_t *)malloc(xfaceSize*sizeof(state_t));
      UxFlux = (state_t *)malloc(xfaceSize*sizeof(state_t));
      VxFlux = (state_t *)malloc(xfaceSize*sizeof(state_t));
      Wx_H = (state_t *)malloc(xfaceSize*sizeof(state_t));
      Wx_U = (state_t *)malloc(xfaceSize*sizeof(state_t));

      HyFlux = (state_t *)malloc(yfaceSize*sizeof(state_t));
      UyFlux = (state_t *)malloc(yfaceSize*sizeof(state_t));
      VyFlux = (state_t *)malloc(yfaceSize*sizeof(state_t));
      Wy_H = (state_t *)malloc(yfaceSize*sizeof(state_t));
      Wy_V = (state_t *)malloc(yfaceSize*sizeof(state_t));

#ifdef _OPENMP
   }
#pragma omp barrier
#endif
   cpu_timers[STATE_TIMER_FINITE_DIFFERENCE_PART1] += cpu_timer_stop(tstart_cpu_part);
   cpu_timer_start(&tstart_cpu_part);

#ifdef _OPENMP
#ifdef _OPENMP_SIMD
#pragma omp for simd
#else
#pragma omp for
#endif
#else
#ifdef _OPENMP_SIMD
#pragma omp simd
#endif
#endif
   //normally use xfaceSize
   for (int iface = 0; iface < mesh->nxface; iface++){
      uint cell_lower = mesh->map_xface2cell_lower[iface];
      uint cell_upper = mesh->map_xface2cell_upper[iface];
      if (cell_lower >= mesh->ncells && cell_upper >= mesh->ncells) continue;

      // set the two faces
      int fl = mesh->map_xcell2face_left1[cell_lower];
      int fr = mesh->map_xcell2face_right1[cell_upper];
#ifndef _OPENMP
      if (fl == -1 || fr == -1) continue;
#endif
      real_t Hx, Ux, Vx;

      // set the two cells away
      int nll = mesh->map_xface2cell_lower[fl];
      int nrr = mesh->map_xface2cell_upper[fr];

      uchar_t lev = mesh->level[cell_lower];
      real_t dxic    = mesh->lev_deltax[lev];
      real_t Cxhalf = 0.5*deltaT/dxic;

      real_t Hic = H[cell_lower];
      real_t Hr  = H[cell_upper];
      real_t Hl  = H[nll];
      real_t Hrr = H[nrr];
      real_t Uic = U[cell_lower];
      real_t Ur  = U[cell_upper];
      real_t Ul  = U[nll];
      real_t Urr = U[nrr];

      Hx=HALF*(H[cell_upper]+H[cell_lower]) - Cxhalf*( HXFLUX(cell_upper)-HXFLUX(cell_lower) );
      Ux=HALF*(U[cell_upper]+U[cell_lower]) - Cxhalf*( UXFLUX(cell_upper)-UXFLUX(cell_lower) );
      Vx=HALF*(V[cell_upper]+V[cell_lower]) - Cxhalf*( UVFLUX(cell_upper)-UVFLUX(cell_lower) );

      real_t U_eigen = fabs(Ux/Hx) + sqrt(g*Hx);

      Wx_H[iface] = w_corrector(deltaT, dxic, U_eigen, Hr-Hic, Hic-Hl, Hrr-Hr) * (Hr - Hic);
      Wx_U[iface] = w_corrector(deltaT, dxic, U_eigen, Ur-Uic, Uic-Ul, Urr-Ur) * (Ur - Uic);

      HxFlux[iface] = HXFLUXFACE;
      UxFlux[iface] = UXFLUXFACE;
      VxFlux[iface] = VXFLUXFACE;

   }

#ifdef _OPENMP
#ifdef _OPENMP_SIMD
#pragma omp for simd
#else
#pragma omp for
#endif
#else
#ifdef _OPENMP_SIMD
#pragma omp simd
#endif
#endif
   for (int ifixup = 0; ifixup < mesh->nxfixup; ifixup++){
      int ir  = mesh->xrecvIdx[ifixup];
      int is1 = mesh->xsendIdx1[ifixup];
      int is2 = mesh->xsendIdx2[ifixup];
      HxFlux[ir] = (HxFlux[is1] + HxFlux[is2]) * HALF;
      UxFlux[ir] = (UxFlux[is1] + UxFlux[is2]) * HALF;
      VxFlux[ir] = (VxFlux[is1] + VxFlux[is2]) * HALF;
      Wx_H[ir] = (Wx_H[is1] + Wx_H[is2]) * 0.25;
      Wx_U[ir] = (Wx_U[is1] + Wx_U[is2]) * 0.25;
   }
   cpu_timers[STATE_TIMER_FINITE_DIFFERENCE_PART2] += cpu_timer_stop(tstart_cpu_part);
   cpu_timer_start(&tstart_cpu_part);


#ifdef _OPENMP
#ifdef _OPENMP_SIMD
#pragma omp for simd
#else
#pragma omp for
#endif
#else
#ifdef _OPENMP_SIMD
#pragma omp simd
#endif
#endif
   //normally use yfaceSize
   for (int iface = 0; iface < mesh->nyface; iface++){
      int cell_lower = mesh->map_yface2cell_lower[iface];
      int cell_upper = mesh->map_yface2cell_upper[iface];
      //if (cell_lower >= mesh->ncells && cell_upper >= mesh->ncells) continue;

      // set the two faces
      int fb = mesh->map_ycell2face_bot1[cell_lower];
      int ft = mesh->map_ycell2face_top1[cell_upper];
#ifndef _OPENMP
//      if (fb == -1 || ft == -1) continue;
#endif
      real_t Hy, Uy, Vy;

      // set the two cells away
      int nbb = mesh->map_yface2cell_lower[fb];
      int ntt = mesh->map_yface2cell_upper[ft];

      uchar_t lev = mesh->level[cell_lower];
      real_t dyic    = mesh->lev_deltay[lev];
      real_t Cyhalf = 0.5*deltaT/dyic;

      real_t Hic = H[cell_lower];
      real_t Ht  = H[cell_upper];
      real_t Hb  = H[nbb];
      real_t Htt = H[ntt];
      real_t Vic = V[cell_lower];
      real_t Vt  = V[cell_upper];
      real_t Vb  = V[nbb];
      real_t Vtt = V[ntt];

      Hy=HALF*(H[cell_upper]+H[cell_lower]) - Cyhalf*( HYFLUX(cell_upper)-HYFLUX(cell_lower) );
      Uy=HALF*(U[cell_upper]+U[cell_lower]) - Cyhalf*( UVFLUX(cell_upper)-UVFLUX(cell_lower) );
      Vy=HALF*(V[cell_upper]+V[cell_lower]) - Cyhalf*( VYFLUX(cell_upper)-VYFLUX(cell_lower) );

      real_t U_eigen = fabs(Vy/Hy) + sqrt(g*Hy);

      Wy_H[iface] = w_corrector(deltaT, dyic, U_eigen, Ht-Hic, Hic-Hb, Htt-Ht) * (Ht - Hic);
      Wy_V[iface] = w_corrector(deltaT, dyic, U_eigen, Vt-Vic, Vic-Vb, Vtt-Vt) * (Vt - Vic);

      HyFlux[iface] = HYFLUXFACE;
      UyFlux[iface] = UYFLUXFACE;
      VyFlux[iface] = VYFLUXFACE;
   }

#ifdef _OPENMP
#ifdef _OPENMP_SIMD
#pragma omp for simd
#else
#pragma omp for
#endif
#else
#ifdef _OPENMP_SIMD
#pragma omp simd
#endif
#endif
   for (int ifixup = 0; ifixup < mesh->nyfixup; ifixup++){
      int ir  = mesh->yrecvIdx[ifixup];
      int is1 = mesh->ysendIdx1[ifixup];
      int is2 = mesh->ysendIdx2[ifixup];
      HyFlux[ir] = (HyFlux[is1] + HyFlux[is2]) * HALF;
      UyFlux[ir] = (UyFlux[is1] + UyFlux[is2]) * HALF;
      VyFlux[ir] = (VyFlux[is1] + VyFlux[is2]) * HALF;
      Wy_H[ir] = (Wy_H[is1] + Wy_H[is2]) * 0.25;
      Wy_V[ir] = (Wy_V[is1] + Wy_V[is2]) * 0.25;
    }
   cpu_timers[STATE_TIMER_FINITE_DIFFERENCE_PART3] += cpu_timer_stop(tstart_cpu_part);
   cpu_timer_start(&tstart_cpu_part);

   static state_t *H_new, *U_new, *V_new;
#ifdef PRECISION_CHECK_GRAPHICS
   static state_t *PCHECK_new;
#endif

#ifdef _OPENMP
#pragma omp barrier
#pragma omp master
   {
#endif
      int flags = (RESTART_DATA | REZONE_DATA | LOAD_BALANCE_MEMORY);

      H_new = (state_t *)state_memory.memory_malloc(mesh->ncells_ghost, sizeof(state_t), "H_new", flags);
      U_new = (state_t *)state_memory.memory_malloc(mesh->ncells_ghost, sizeof(state_t), "U_new", flags);
      V_new = (state_t *)state_memory.memory_malloc(mesh->ncells_ghost, sizeof(state_t), "V_new", flags);
#ifdef PRECISION_CHECK_GRAPHICS
      PCHECK_new = (state_t *)state_memory.memory_malloc(mesh->ncells_ghost, sizeof(state_t), "PCHECK_new", flags);
#endif
#ifdef _OPENMP
   }
#pragma omp barrier
#endif


   int lowerBound, upperBound;

   mesh->get_bounds(lowerBound, upperBound);

#if defined(__GNUC_MINOR__)
   static real_t * dxcell;
#ifdef _OPENMP
#pragma omp barrier
#pragma omp master
   {
#endif
   dxcell = (real_t *)malloc(sizeof(real_t) * mesh->ncells);
#ifdef _OPENMP
   }
#pragma omp barrier
#endif
   for (int ic = lowerBound; ic < upperBound; ic++) {
      uchar_t lev = mesh->level[ic];
      dxcell[ic]  = mesh->lev_deltax[lev];
   }

   state_t *H_loc = H;
   state_t *U_loc = U;
   state_t *V_loc = V;

   int *map_xcell2face_left1_loc = mesh->map_xcell2face_left1;
   int *map_xcell2face_right1_loc = mesh->map_xcell2face_right1;
   int *map_ycell2face_bot1_loc = mesh->map_ycell2face_bot1;
   int *map_ycell2face_top1_loc = mesh->map_ycell2face_top1;
#endif

#ifdef _OPENMP_SIMD
#pragma omp simd
#endif
   for (int ic = lowerBound; ic < upperBound; ic++){
      if (mesh->celltype[ic] != REAL_CELL) continue;

#if defined(__GNUC_MINOR__)
      real_t dxic = dxcell[ic];
#else
      real_t dxic = mesh->lev_deltax[mesh->level[ic]];
#endif
      // set the four faces
      int fl = mesh->map_xcell2face_left1[ic];
      int fr = mesh->map_xcell2face_right1[ic];
      int fb = mesh->map_ycell2face_bot1[ic];
      int ft = mesh->map_ycell2face_top1[ic];

      H_new[ic] = U_fullstep(deltaT,dxic,H[ic],
                  HxFlux[fr], HxFlux[fl], HyFlux[ft], HyFlux[fb])
                - Wx_H[fl] + Wx_H[fr] - Wy_H[fb] + Wy_H[ft];
      U_new[ic] = U_fullstep(deltaT,dxic,U[ic],
                  UxFlux[fr], UxFlux[fl], UyFlux[ft], UyFlux[fb])
                - Wx_U[fl] + Wx_U[fr];
      V_new[ic] = U_fullstep(deltaT,dxic,V[ic],
                  VxFlux[fr], VxFlux[fl], VyFlux[ft], VyFlux[fb])
                - Wy_V[fb] + Wy_V[ft];

   } // cell loop

   cpu_timers[STATE_TIMER_FINITE_DIFFERENCE_PART4] += cpu_timer_stop(tstart_cpu_part);
   cpu_timer_start(&tstart_cpu_part);

#ifdef _OPENMP
#pragma omp barrier
#pragma omp master
   {
#endif
      free(HxFlux);
      free(UxFlux);
      free(VxFlux);
      free(Wx_H);
      free(Wx_U);
      free(HyFlux);
      free(UyFlux);
      free(VyFlux);
      free(Wy_H);
      free(Wy_V);

      // Replace H with H_new and deallocate H. New memory will have the characteristics
      // of the new memory and the name of the old. Both return and arg1 will be reset to new memory
      H = (state_t *)state_memory.memory_replace(H, H_new);
      U = (state_t *)state_memory.memory_replace(U, U_new);
      V = (state_t *)state_memory.memory_replace(V, V_new);
#ifdef PRECISION_CHECK_GRAPHICS
      PCHECK = (state_t *)state_memory.memory_replace(PCHECK, PCHECK_new);
#endif

#if defined(__GNUC_MINOR__)
      free(dxcell);
#endif

      cpu_timers[STATE_TIMER_FINITE_DIFFERENCE] += cpu_timer_stop(tstart_cpu);
#ifdef _OPENMP
   }
#pragma omp barrier
#endif
}

void State::calc_finite_difference_via_faces(double deltaT)
{
   real_t   g     = 9.80;   // gravitational constant
   real_t   ghalf = HALF*g;

   struct timespec tstart_cpu;
   cpu_timer_start(&tstart_cpu);

   struct timespec tstart_cpu_part;
   cpu_timer_start(&tstart_cpu_part);

   // We need to populate the ghost regions since the calc neighbors has just been
   // established for the mesh shortly before
   apply_boundary_conditions();

   //No longer used, can remove
   //int *nlft, *nbot;
   int *nrht, *ntop;
   uchar_t *level;

   //nlft  = mesh->nlft;
   nrht  = mesh->nrht;
   //nbot  = mesh->nbot;
   ntop  = mesh->ntop;
   level = mesh->level;

   vector<real_t> &lev_deltax = mesh->lev_deltax;
   vector<real_t> &lev_deltay = mesh->lev_deltay;

   static int xfaceSize;
   static int yfaceSize;

   static vector<state_t> HxFlux, UxFlux, VxFlux, Wx_H, Wx_U;
   static vector<state_t> HyFlux, UyFlux, VyFlux, Wy_H, Wy_V;

#ifdef _OPENMP
#pragma omp barrier
#pragma omp master
   {
#endif
      mesh->calc_face_list_wbidirmap();
      xfaceSize = mesh->nxface;//new "update" nxface inc. phantoms
      yfaceSize = mesh->nyface; //new "update" nyface inc. phantoms
      memory_reset_ptrs(); //reset the pointers H,U,V that were recently reallocated in wbidirmap call

      HxFlux.resize(xfaceSize, (state_t)-999999);
      UxFlux.resize(xfaceSize, (state_t)-999999);
      VxFlux.resize(xfaceSize, (state_t)-999999);
      Wx_H.resize(xfaceSize, (state_t)-999999);
      Wx_U.resize(xfaceSize, (state_t)-999999);


      HyFlux.resize(yfaceSize, (state_t)-999999);
      UyFlux.resize(yfaceSize, (state_t)-999999);
      VyFlux.resize(yfaceSize, (state_t)-999999);
      Wy_H.resize(yfaceSize, (state_t)-999999);
      Wy_V.resize(yfaceSize, (state_t)-999999);
#ifdef _OPENMP
   }
#pragma omp barrier
#endif

   cpu_timers[STATE_TIMER_FINITE_DIFFERENCE_PART1] += cpu_timer_stop(tstart_cpu_part);
   cpu_timer_start(&tstart_cpu_part);


#ifdef _OPENMP
#ifdef _OPENMP_SIMD
#pragma omp for simd
#else
#pragma omp for
#endif
#else
#ifdef _OPENMP_SIMD
#pragma omp simd
#endif
#endif
   for (int iface = 0; iface < mesh->nxface; iface++){
      int cell_lower = mesh->map_xface2cell_lower[iface];
      int cell_upper = mesh->map_xface2cell_upper[iface];
      //printf("%d) %d %d\n", iface, cell_lower, cell_upper);
      real_t Hx, Ux, Vx;
      if (level[cell_lower] == level[cell_upper]) {

         // set the two faces
         //int fl = mesh->map_xcell2face_left1[cell_lower];
         //int fr = mesh->map_xcell2face_right1[cell_upper];
         // set the two cells away
         //int nll = mesh->map_xface2cell_lower[fl];
         int nll = mesh->nlft[cell_lower];
         int nrr = mesh->nrht[cell_upper];
 
         int lev = level[cell_lower];
         real_t dxic = lev_deltax[lev];
         real_t Cxhalf = 0.5*deltaT/dxic;
 
         real_t Hic = H[cell_lower];
         real_t Hr  = H[cell_upper];
         real_t Hl  = H[nll];
         real_t Hrr = H[nrr];
         real_t Uic = U[cell_lower];
         real_t Ur  = U[cell_upper];
         real_t Ul  = U[nll];
         real_t Urr = U[nrr];
 
         Hx=HALF*(H[cell_upper]+H[cell_lower]) - Cxhalf*( HXFLUX(cell_upper)-HXFLUX(cell_lower) );
         Ux=HALF*(U[cell_upper]+U[cell_lower]) - Cxhalf*( UXFLUX(cell_upper)-UXFLUX(cell_lower) );
         Vx=HALF*(V[cell_upper]+V[cell_lower]) - Cxhalf*( UVFLUX(cell_upper)-UVFLUX(cell_lower) );
 
         real_t U_eigen = fabs(Ux/Hx) + sqrt(g*Hx);
 
         Wx_H[iface] = w_corrector(deltaT, dxic, U_eigen, Hr-Hic, Hic-Hl, Hrr-Hr) * (Hr - Hic);
 
         Wx_U[iface] = w_corrector(deltaT, dxic, U_eigen, Ur-Uic, Uic-Ul, Urr-Ur) * (Ur - Uic);

         HxFlux[iface] = HXFLUXFACE;
         UxFlux[iface] = UXFLUXFACE;
         VxFlux[iface] = VXFLUXFACE;
      } else {

         real_t dx_lower = lev_deltax[level[cell_lower]];
         real_t dx_upper = lev_deltax[level[cell_upper]];

         real_t FA_lower = dx_lower;
         real_t FA_upper = dx_upper;
         real_t FA_lolim = FA_lower*min(ONE, FA_upper/FA_lower);
         real_t FA_uplim = FA_upper*min(ONE, FA_lower/FA_upper);

         real_t CV_lower = SQ(dx_lower);
         real_t CV_upper = SQ(dx_upper);
         real_t CV_lolim = CV_lower*min(HALF, CV_upper/CV_lower);
         real_t CV_uplim = CV_upper*min(HALF, CV_lower/CV_upper);

         // Weighted half-step calculation
         //
         // (dx_lower*H[cell_upper]+dx_upper*H[cell_lower])
         // -----------------------------------------------   -
         //             (dx_lower+dx_upper)
         //
         //                ( (FA_uplim*HXFLUX(cell_upper))-(FA_lolim*HXFLUX(cell_lower)) )
         // 0.5*deltaT  *  ----------------------------------------------------------------
         //                                    (CV_uplim+CV_lolim)
         //

         Hx=(dx_lower*H[cell_upper]+dx_upper*H[cell_lower])/(dx_lower+dx_upper) -
                   HALF*deltaT*( (FA_uplim*HXFLUX(cell_upper))-(FA_lolim*HXFLUX(cell_lower)) )/
                   (CV_uplim+CV_lolim);
         Ux=(dx_lower*U[cell_upper]+dx_upper*U[cell_lower])/(dx_lower+dx_upper) -
                   HALF*deltaT*( (FA_uplim*UXFLUX(cell_upper))-(FA_lolim*UXFLUX(cell_lower)) )/
                   (CV_uplim+CV_lolim);
         Vx=(dx_lower*V[cell_upper]+dx_upper*V[cell_lower])/(dx_lower+dx_upper) -
                   HALF*deltaT*( (FA_uplim*UVFLUX(cell_upper))-(FA_lolim*UVFLUX(cell_lower)) )/
                   (CV_uplim+CV_lolim); 

         // set the two faces
         //int fl = mesh->map_xcell2face_left1[cell_lower];
         //int fr = mesh->map_xcell2face_right1[cell_upper];
         // set the two cells away
         //int nll = mesh->map_xface2cell_lower[fl];
         //int nrr = mesh->map_xface2cell_upper[fr];
         int nll = mesh->nlft[cell_lower];
         int nrr = mesh->nrht[cell_upper];

         uchar_t lev = level[cell_lower];
         uchar_t levr = level[cell_upper];
         real_t dxic = lev_deltax[lev]; 
         real_t dxr = lev_deltax[levr];

         real_t Hic = H[cell_lower];
         real_t Hr  = H[cell_upper];
         real_t Hl  = H[nll];
         real_t Hrr = H[nrr];
         real_t Uic = U[cell_lower];
         real_t Ur  = U[cell_upper];
         real_t Ul  = U[nll];
         real_t Urr = U[nrr];

         real_t U_eigen = fabs(Ux/Hx) + sqrt(g*Hx);
         real_t dx_avg = (dxic+dxr)*HALF;

         if(level[cell_upper] < level[nrr]) {
            Hrr = (Hrr + H[ntop[nrr]]) * HALF;
            Urr = (Urr + U[ntop[nrr]]) * HALF;
         }

         real_t Hl2 = Hl;
         real_t Ul2 = Ul;
         if(lev < level[nll]) {
            Hl2 = (Hl2 + H[ntop[nll]]) * HALF;
            Ul2 = (Ul2 + U[ntop[nll]]) * HALF;
         }

         Wx_H[iface] = w_corrector(deltaT, dx_avg, U_eigen, Hr-Hic, Hic-Hl2, Hrr-Hr) * (Hr - Hic);
         Wx_U[iface] = w_corrector(deltaT, dx_avg, U_eigen, Ur-Uic, Uic-Ul2, Urr-Ur) * (Ur - Uic);
 
         HxFlux[iface] = HXFLUXFACE;
         UxFlux[iface] = UXFLUXFACE;
         VxFlux[iface] = VXFLUXFACE;
      }
      //printf("\t%d) %f | %f | %f\n", iface, HxFlux[iface], UxFlux[iface], VxFlux[iface]);
   }

   cpu_timers[STATE_TIMER_FINITE_DIFFERENCE_PART2] += cpu_timer_stop(tstart_cpu_part);
   cpu_timer_start(&tstart_cpu_part);


#ifdef _OPENMP
#ifdef _OPENMP_SIMD
#pragma omp for simd
#else
#pragma omp for
#endif
#else
#ifdef _OPENMP_SIMD
#pragma omp simd
#endif
#endif
   for (int iface = 0; iface < mesh->nyface; iface++){
      int cell_lower = mesh->map_yface2cell_lower[iface];
      int cell_upper = mesh->map_yface2cell_upper[iface];
      real_t Hy, Uy, Vy;
      if (level[cell_lower] == level[cell_upper]) {

         // set the two faces
         //int fb = mesh->map_ycell2face_bot1[cell_lower];
         //int ft = mesh->map_ycell2face_top1[cell_upper];
         // set the two cells away
         //int nbb = mesh->map_yface2cell_lower[fb];
         //int ntt = mesh->map_yface2cell_upper[ft];
         int nbb = mesh->nbot[cell_lower];
         int ntt = mesh->ntop[cell_upper];
	
         int lev = level[cell_lower];
         real_t dyic    = lev_deltay[lev];
         real_t Cyhalf = 0.5*deltaT/lev_deltay[lev];

         real_t Hic = H[cell_lower];
         real_t Ht  = H[cell_upper];
         real_t Hb  = H[nbb];
         real_t Htt = H[ntt];
         real_t Vic = V[cell_lower];
         real_t Vt  = V[cell_upper];
         real_t Vb  = V[nbb];
         real_t Vtt = V[ntt];

         Hy=HALF*(H[cell_upper]+H[cell_lower]) - Cyhalf*( HYFLUX(cell_upper)-HYFLUX(cell_lower) );
         Uy=HALF*(U[cell_upper]+U[cell_lower]) - Cyhalf*( UVFLUX(cell_upper)-UVFLUX(cell_lower) );
         Vy=HALF*(V[cell_upper]+V[cell_lower]) - Cyhalf*( VYFLUX(cell_upper)-VYFLUX(cell_lower) );

         real_t U_eigen = fabs(Vy/Hy) + sqrt(g*Hy);
 
         Wy_H[iface] = w_corrector(deltaT, dyic, U_eigen, Ht-Hic, Hic-Hb, Htt-Ht) * (Ht - Hic);
         Wy_V[iface] = w_corrector(deltaT, dyic, U_eigen, Vt-Vic, Vic-Vb, Vtt-Vt) * (Vt - Vic);

         HyFlux[iface] = HYFLUXFACE;
         UyFlux[iface] = UYFLUXFACE;
         VyFlux[iface] = VYFLUXFACE;
      } else {
         real_t dy_lower = lev_deltay[level[cell_lower]];
         real_t dy_upper = lev_deltay[level[cell_upper]];

         real_t FA_lower = dy_lower;
         real_t FA_upper = dy_upper;
         real_t FA_lolim = FA_lower*min(ONE, FA_upper/FA_lower);
         real_t FA_uplim = FA_upper*min(ONE, FA_lower/FA_upper);

         real_t CV_lower = SQ(dy_lower);
         real_t CV_upper = SQ(dy_upper);
         real_t CV_lolim = CV_lower*min(HALF, CV_upper/CV_lower);
         real_t CV_uplim = CV_upper*min(HALF, CV_lower/CV_upper);

         // Weighted half-step calculation
         //
         // (dy_lower*H[cell_upper]+dy_upper*H[cell_lower])
         // -----------------------------------------------   -
         //             (dy_lower+dy_upper)
         //
         //                ( (FA_uplim*HYFLUX(cell_upper))-(FA_lolim*HYFLUX(cell_lower)) )
         // 0.5*deltaT  *  ----------------------------------------------------------------
         //                                    (CV_uplim+CV_lolim)
         //

         Hy=(dy_lower*H[cell_upper]+dy_upper*H[cell_lower])/(dy_lower+dy_upper) -
                   HALF*deltaT*( (FA_uplim*HYFLUX(cell_upper))-(FA_lolim*HYFLUX(cell_lower)) )/
                   (CV_uplim+CV_lolim);
         Uy=(dy_lower*U[cell_upper]+dy_upper*U[cell_lower])/(dy_lower+dy_upper) -
                   HALF*deltaT*( (FA_uplim*UVFLUX(cell_upper))-(FA_lolim*UVFLUX(cell_lower)) )/
                   (CV_uplim+CV_lolim);
         Vy=(dy_lower*V[cell_upper]+dy_upper*V[cell_lower])/(dy_lower+dy_upper) -
                   HALF*deltaT*( (FA_uplim*VYFLUX(cell_upper))-(FA_lolim*VYFLUX(cell_lower)) )/
                   (CV_uplim+CV_lolim);

         // set the two faces
         //int fb = mesh->map_ycell2face_bot1[cell_lower];
         //int ft = mesh->map_ycell2face_top1[cell_upper];
         // set the two cells away
         //int nbb = mesh->map_yface2cell_lower[fb];
         //int ntt = mesh->map_yface2cell_upper[ft];
         int nbb = mesh->nbot[cell_lower];
         int ntt = mesh->ntop[cell_upper];
 
         uchar_t lev = level[cell_lower];
         uchar_t levt = level[cell_upper];
         real_t dyic = lev_deltay[lev];
         real_t dyt = lev_deltay[levt];

         real_t Hic = H[cell_lower];
         real_t Ht  = H[cell_upper];
         real_t Hb  = H[nbb];
         real_t Htt = H[ntt];
         real_t Vic = V[cell_lower];
         real_t Vt  = V[cell_upper];
         real_t Vb  = V[nbb];
         real_t Vtt = V[ntt];

         real_t V_eigen = fabs(Vy/Hy) + sqrt(g*Hy);
         real_t dy_avg = (dyic+dyt)*HALF; 

         if(level[cell_upper] < level[ntt]) {
            Htt = (Htt + H[nrht[ntt]]) * HALF;
            Vtt = (Vtt + V[nrht[ntt]]) * HALF;
         } 

         real_t Hb2 = Hb;
         real_t Vb2 = Vb;
         if(lev < level[nbb]) {
            Hb2 = (Hb2 + H[nrht[nbb]]) * HALF;
            Vb2 = (Vb2 + V[nrht[nbb]]) * HALF;
         }

         Wy_H[iface] = w_corrector(deltaT, dy_avg, V_eigen, Ht-Hic, Hic-Hb2, Htt-Ht) * (Ht - Hic);
         Wy_V[iface] = w_corrector(deltaT, dy_avg, V_eigen, Vt-Vic, Vic-Vb2, Vtt-Vt) * (Vt - Vic);
 
         HyFlux[iface] = HYFLUXFACE;
         UyFlux[iface] = UYFLUXFACE;
         VyFlux[iface] = VYFLUXFACE;
      }
   }

   cpu_timers[STATE_TIMER_FINITE_DIFFERENCE_PART3] += cpu_timer_stop(tstart_cpu_part);
   cpu_timer_start(&tstart_cpu_part);

   static state_t *H_new, *U_new, *V_new;
#ifdef PRECISION_CHECK_GRAPHICS
   static state_t *PCHECK_new;
#endif

#ifdef _OPENMP
#pragma omp barrier
#pragma omp master
   {
#endif
      int flags = (RESTART_DATA | REZONE_DATA | LOAD_BALANCE_MEMORY);

      H_new = (state_t *)state_memory.memory_malloc(mesh->ncells_ghost, sizeof(state_t), "H_new", flags);
      U_new = (state_t *)state_memory.memory_malloc(mesh->ncells_ghost, sizeof(state_t), "U_new", flags);
      V_new = (state_t *)state_memory.memory_malloc(mesh->ncells_ghost, sizeof(state_t), "V_new", flags);
#ifdef PRECISION_CHECK_GRAPHICS
      PCHECK_new = (state_t *)state_memory.memory_malloc(mesh->ncells_ghost, sizeof(state_t), "PCHECK_new", flags);
#endif
#ifdef _OPENMP
   }
#pragma omp barrier
#endif

   int lowerBound, upperBound;

   mesh->get_bounds(lowerBound, upperBound);
#ifdef PRECISION_CHECK_STATS
   fail_prec_count = 0;
   fail_F_plus_sum    = 0.0;
   fail_F_minus_sum   = 0.0;
   fail_G_plus_sum    = 0.0;
   fail_G_minus_sum   = 0.0;
   fail_wminusx_H_sum = 0.0;
   fail_wplusx_H_sum  = 0.0;
   fail_wminusy_H_sum = 0.0;
   fail_wplusy_H_sum  = 0.0;
   prec_count = 0;
   F_plus_sum    = 0.0;
   F_minus_sum   = 0.0;
   G_plus_sum    = 0.0;
   G_minus_sum   = 0.0;
   wminusx_H_sum = 0.0;
   wplusx_H_sum  = 0.0;
   wminusy_H_sum = 0.0;
   wplusy_H_sum  = 0.0;
#endif

#ifdef _OPENMP_SIMD
#pragma omp simd
#endif
   for (int ic = lowerBound; ic < upperBound; ic++){
      real_t dxic    = lev_deltax[level[ic]];
      // set the four faces
      int fl = mesh->map_xcell2face_left1[ic];
      int fr = mesh->map_xcell2face_right1[ic];
      int fb = mesh->map_ycell2face_bot1[ic];
      int ft = mesh->map_ycell2face_top1[ic];
      int fl2 = mesh->map_xcell2face_left2[ic];
      int fr2 = mesh->map_xcell2face_right2[ic];
      int fb2 = mesh->map_ycell2face_bot2[ic];
      int ft2 = mesh->map_ycell2face_top2[ic];

      //printf("%d) %d %d %d %d\n", ic, fl, fr, fb, ft);

      // set the four neighboring cells
      int nl = mesh->map_xface2cell_lower[fl];
      int nr = mesh->map_xface2cell_upper[fr];
      int nb = mesh->map_yface2cell_lower[fb];
      int nt = mesh->map_yface2cell_upper[ft];

      if (nb == ic  || nt == ic || nl == ic || nr == ic) continue;

      real_t Hic     = H[ic];
      real_t Uic     = U[ic];
      real_t Vic     = V[ic];

      real_t Hxfluxminus = HxFlux[fl];
      real_t Uxfluxminus = UxFlux[fl];
      real_t Vxfluxminus = VxFlux[fl];

      real_t Hxfluxplus  = HxFlux[fr];
      real_t Uxfluxplus  = UxFlux[fr];
      real_t Vxfluxplus  = VxFlux[fr];

      real_t Hyfluxminus = HyFlux[fb];
      real_t Uyfluxminus = UyFlux[fb];
      real_t Vyfluxminus = VyFlux[fb];

      real_t Hyfluxplus  = HyFlux[ft];
      real_t Uyfluxplus  = UyFlux[ft];
      real_t Vyfluxplus  = VyFlux[ft];

      real_t wminusx_H = Wx_H[fl];
      real_t wminusx_U = Wx_U[fl];

      real_t wplusx_H = Wx_H[fr];
      real_t wplusx_U = Wx_U[fr];

      real_t wminusy_H = Wy_H[fb];
      real_t wminusy_V = Wy_V[fb];

      real_t wplusy_H = Wy_H[ft];
      real_t wplusy_V = Wy_V[ft];

      if (level[ic] < level[nl]) {
         Hxfluxminus = (Hxfluxminus + HxFlux[fl2]) * HALF;
         Uxfluxminus = (Uxfluxminus + UxFlux[fl2]) * HALF;
         Vxfluxminus = (Vxfluxminus + VxFlux[fl2]) * HALF;
         wminusx_H = (wminusx_H + Wx_H[fl2]) * HALF * HALF;
         wminusx_U = (wminusx_U + Wx_U[fl2]) * HALF * HALF;
      }
   
      if (level[ic] < level[nr]) {
         Hxfluxplus = (Hxfluxplus + HxFlux[fr2]) * HALF;
         Uxfluxplus = (Uxfluxplus + UxFlux[fr2]) * HALF;
         Vxfluxplus = (Vxfluxplus + VxFlux[fr2]) * HALF;
         wplusx_H = (wplusx_H + Wx_H[fr2]) * HALF * HALF;
         wplusx_U = (wplusx_U + Wx_U[fr2]) * HALF * HALF;
      }
   
      if (level[ic] < level[nb]) {
         Hyfluxminus = (Hyfluxminus + HyFlux[fb2]) * HALF;
         Uyfluxminus = (Uyfluxminus + UyFlux[fb2]) * HALF;
         Vyfluxminus = (Vyfluxminus + VyFlux[fb2]) * HALF;
         wminusy_H = (wminusy_H + Wy_H[fb2]) * HALF * HALF;
         wminusy_V = (wminusy_V + Wy_V[fb2]) * HALF * HALF;
      }
   
      if (level[ic] < level[nt]) {
         Hyfluxplus = (Hyfluxplus + HyFlux[ft2]) * HALF;
         Uyfluxplus = (Uyfluxplus + UyFlux[ft2]) * HALF;
         Vyfluxplus = (Vyfluxplus + VyFlux[ft2]) * HALF;
         wplusy_H = (wplusy_H + Wy_H[ft2]) * HALF * HALF;
         wplusy_V = (wplusy_V + Wy_V[ft2]) * HALF * HALF;
      }


#ifdef PRECISION_CHECK_WITH_PARENTHESIS
   //Some parentheses
      H_new[ic] = U_fullstep(deltaT, dxic, Hic,
                      Hxfluxplus, Hxfluxminus, Hyfluxplus, Hyfluxminus)
                 + (- wminusx_H + wplusx_H - wminusy_H + wplusy_H);
#else

#ifdef PRECISION_CHECK_BEST_PARENTHESIS
   //"best" parentheses version
      H_new[ic] = U_fullstep(deltaT, dxic, Hic,
                      Hxfluxplus, Hxfluxminus, Hyfluxplus, Hyfluxminus)
                   //+ (((-wminusx_H - wminusy_H) + wplusy_H) + wplusx_H);
                   + (((wplusx_H + wplusy_H) - wminusy_H) - wminusx_H);
#else
   //original, no parentheses
      H_new[ic] = U_fullstep(deltaT, dxic, Hic,
                      Hxfluxplus, Hxfluxminus, Hyfluxplus, Hyfluxminus)
                 - wminusx_H + wplusx_H - wminusy_H + wplusy_H;
#endif
#endif

#ifdef PRECISION_CHECK
      int fail;
      U_fullstep_precision_check(ic, deltaT, dxic, Hic, H_new[ic],
                      Hxfluxplus, Hxfluxminus, Hyfluxplus, Hyfluxminus,
                      wminusx_H, wplusx_H, wminusy_H, wplusy_H, &fail);


#ifdef PRECISION_CHECK_GRAPHICS
      PCHECK_new[ic] = 100.0*(double)fail;
#endif
#endif

      U_new[ic] = U_fullstep(deltaT, dxic, Uic,
                      Uxfluxplus, Uxfluxminus, Uyfluxplus, Uyfluxminus)
                 - wminusx_U + wplusx_U;
      V_new[ic] = U_fullstep(deltaT, dxic, Vic,
                      Vxfluxplus, Vxfluxminus, Vyfluxplus, Vyfluxminus)
                 - wminusy_V + wplusy_V;
      //printf("%d) %f | %f | %f\n", ic, H_new[ic], U_new[ic], V_new[ic]);
#ifdef HAVE_MPI
         if (mesh->mype == 1) {
             int nb, nbb, nt, ntt;
             real_t Hb, Hbb, Ht, Htt;
             nb = mesh->nbot[ic];
             nbb = mesh->nbot[nb];
             nt = mesh->ntop[ic];
             ntt = mesh->ntop[nt];
             Hb = H[nb];
             Hb = H[nb];
             Hbb = H[nbb];
             Ht = H[nt];
             Htt = H[ntt];
             //if (ic == 48) printf("%d %d %d %d %lf %lf %lf %lf %lf\n", nb, nbb, nt, ntt, Hb, Hbb, Ht, Htt, wplusy_H);
#endif
            //printf("%d) %d %d %f %f %f %f %f %f %f %f %f\n", ic, mesh->i[ic], mesh->j[ic], Hic, Hxfluxminus, Hxfluxplus, Hyfluxminus, Hyfluxplus, wminusx_H, wplusx_H, wminusy_H, wplusy_H); 
            //printf("%d) %d %d %d %d %d %d\n", ic, mesh->i[ic], mesh->j[ic], mesh->nlft[ic], mesh->nrht[ic], mesh->nbot[ic], mesh->ntop[ic]);
#ifdef HAVE_MPI
            }
#endif
            


   } // cell loop

#ifdef PRECISION_CHECK_STATS
   fail_F_plus_sum    /= (double)fail_prec_count;
   fail_F_minus_sum   /= (double)fail_prec_count;
   fail_G_plus_sum    /= (double)fail_prec_count;
   fail_G_minus_sum   /= (double)fail_prec_count;
   fail_wminusx_H_sum /= (double)fail_prec_count;
   fail_wplusx_H_sum  /= (double)fail_prec_count;
   fail_wminusy_H_sum /= (double)fail_prec_count;
   fail_wplusy_H_sum  /= (double)fail_prec_count;
   F_plus_sum    /= (double)prec_count;
   F_minus_sum   /= (double)prec_count;
   G_plus_sum    /= (double)prec_count;
   G_minus_sum   /= (double)prec_count;
   wminusx_H_sum /= (double)prec_count;
   wplusx_H_sum  /= (double)prec_count;
   wminusy_H_sum /= (double)prec_count;
   wplusy_H_sum  /= (double)prec_count;

   fail_prec_avg_count++;
   fail_F_plus_avg    += fail_F_plus_sum;
   fail_F_minus_avg   += fail_F_minus_sum;
   fail_G_plus_avg    += fail_G_plus_sum;
   fail_G_minus_avg   += fail_G_minus_sum;
   fail_wminusx_H_avg += fail_wminusx_H_sum;
   fail_wplusx_H_avg  += fail_wplusx_H_sum;
   fail_wminusy_H_avg += fail_wminusy_H_sum;
   fail_wplusy_H_avg  += fail_wplusy_H_sum;
   prec_avg_count++;
   F_plus_avg    += F_plus_sum;
   F_minus_avg   += F_minus_sum;
   G_plus_avg    += G_plus_sum;
   G_minus_avg   += G_minus_sum;
   wminusx_H_avg += wminusx_H_sum;
   wplusx_H_avg  += wplusx_H_sum;
   wminusy_H_avg += wminusy_H_sum;
   wplusy_H_avg  += wplusy_H_sum;
#endif

   cpu_timers[STATE_TIMER_FINITE_DIFFERENCE_PART4] += cpu_timer_stop(tstart_cpu_part);
   cpu_timer_start(&tstart_cpu_part);


#ifdef _OPENMP
#pragma omp barrier
#pragma omp master
   {
#endif
      // Replace H with H_new and deallocate H. New memory will have the characteristics
      // of the new memory and the name of the old. Both return and arg1 will be reset to new memory
      H = (state_t *)state_memory.memory_replace(H, H_new);
      U = (state_t *)state_memory.memory_replace(U, U_new);
      V = (state_t *)state_memory.memory_replace(V, V_new);
#ifdef PRECISION_CHECK_GRAPHICS
      PCHECK = (state_t *)state_memory.memory_replace(PCHECK, PCHECK_new);
#endif

      cpu_timers[STATE_TIMER_FINITE_DIFFERENCE] += cpu_timer_stop(tstart_cpu);
#ifdef _OPENMP
   }
#pragma omp barrier
#endif
}

#define HXRGFLUXIC ( U_reg_lev[ll][jj][ii] )
#define HXRGFLUXNL ( U_reg_lev[ll][jj][ii-1] )
#define HXRGFLUXNR ( U_reg_lev[ll][jj][ii+1] )
#define HXRGFLUXNB ( U_reg_lev[ll][jj-1][ii] )
#define HXRGFLUXNT ( U_reg_lev[ll][jj+1][ii] )

#define UXRGFLUXIC ( SQ(U_reg_lev[ll][jj][ii])  /H_reg_lev[ll][jj][ii]   + ghalf*SQ(H_reg_lev[ll][jj][ii]) )
#define UXRGFLUXNL ( SQ(U_reg_lev[ll][jj][ii-1])/H_reg_lev[ll][jj][ii-1] + ghalf*SQ(H_reg_lev[ll][jj][ii-1]) )
#define UXRGFLUXNR ( SQ(U_reg_lev[ll][jj][ii+1])/H_reg_lev[ll][jj][ii+1] + ghalf*SQ(H_reg_lev[ll][jj][ii+1]) )
#define UXRGFLUXNB ( SQ(U_reg_lev[ll][jj-1][ii])/H_reg_lev[ll][jj-1][ii] + ghalf*SQ(H_reg_lev[ll][jj-1][ii]) )
#define UXRGFLUXNT ( SQ(U_reg_lev[ll][jj+1][ii])/H_reg_lev[ll][jj+1][ii] + ghalf*SQ(H_reg_lev[ll][jj+1][ii]) )

#define VXRGFLUXIC ( U_reg_lev[ll][jj][ii]  *V_reg_lev[ll][jj][ii]  /H_reg_lev[ll][jj][ii] )
#define VXRGFLUXNL ( U_reg_lev[ll][jj][ii-1]*V_reg_lev[ll][jj][ii-1]/H_reg_lev[ll][jj][ii-1] )
#define VXRGFLUXNR ( U_reg_lev[ll][jj][ii+1]*V_reg_lev[ll][jj][ii+1]/H_reg_lev[ll][jj][ii+1] )
#define VXRGFLUXNB ( U_reg_lev[ll][jj-1][ii]*V_reg_lev[ll][jj-1][ii]/H_reg_lev[ll][jj-1][ii] )
#define VXRGFLUXNT ( U_reg_lev[ll][jj+1][ii]*V_reg_lev[ll][jj+1][ii]/H_reg_lev[ll][jj+1][ii] )

#define HYRGFLUXIC ( V_reg_lev[ll][jj][ii] )
#define HYRGFLUXNL ( V_reg_lev[ll][jj][ii-1] )
#define HYRGFLUXNR ( V_reg_lev[ll][jj][ii+1] )
#define HYRGFLUXNB ( V_reg_lev[ll][jj-1][ii] )
#define HYRGFLUXNT ( V_reg_lev[ll][jj+1][ii] )

#define UYRGFLUXIC  ( V_reg_lev[ll][jj][ii]  *U_reg_lev[ll][jj][ii]  /H_reg_lev[ll][jj][ii] )
#define UYRGFLUXNL  ( V_reg_lev[ll][jj][ii-1]*U_reg_lev[ll][jj][ii-1]/H_reg_lev[ll][jj][ii-1] )
#define UYRGFLUXNR  ( V_reg_lev[ll][jj][ii+1]*U_reg_lev[ll][jj][ii+1]/H_reg_lev[ll][jj][ii+1] )
#define UYRGFLUXNB  ( V_reg_lev[ll][jj-1][ii]*U_reg_lev[ll][jj-1][ii]/H_reg_lev[ll][jj-1][ii] )
#define UYRGFLUXNT  ( V_reg_lev[ll][jj+1][ii]*U_reg_lev[ll][jj+1][ii]/H_reg_lev[ll][jj+1][ii] )

#define VYRGFLUXIC  ( SQ(V_reg_lev[ll][jj][ii])  /H_reg_lev[ll][jj][ii]   + ghalf*SQ(H_reg_lev[ll][jj][ii]) )
#define VYRGFLUXNL  ( SQ(V_reg_lev[ll][jj][ii-1])/H_reg_lev[ll][jj][ii-1] + ghalf*SQ(H_reg_lev[ll][jj][ii-1]) )
#define VYRGFLUXNR  ( SQ(V_reg_lev[ll][jj][ii+1])/H_reg_lev[ll][jj][ii+1] + ghalf*SQ(H_reg_lev[ll][jj][ii+1]) )
#define VYRGFLUXNB  ( SQ(V_reg_lev[ll][jj-1][ii])/H_reg_lev[ll][jj-1][ii] + ghalf*SQ(H_reg_lev[ll][jj-1][ii]) )
#define VYRGFLUXNT  ( SQ(V_reg_lev[ll][jj+1][ii])/H_reg_lev[ll][jj+1][ii] + ghalf*SQ(H_reg_lev[ll][jj+1][ii]) )

#define HNEWXRGFLUXFL  ( Ux )
#define UNEWXRGFLUXFL  ( SQ(Ux)/Hx + ghalf*SQ(Hx) )
#define VNEWXRGFLUXFL  ( Ux*Vx/Hx )

#define HNEWYRGFLUXFB  ( Vy )
#define UNEWYRGFLUXFB  ( Vy*Uy/Hy )
#define VNEWYRGFLUXFB  ( SQ(Vy)/Hy + ghalf*SQ(Hy) )

#define HNEWXRGFLUXMINUS  ( Uxminus )
#define HNEWXRGFLUXPLUS   ( Uxplus )
#define UNEWXRGFLUXMINUS  ( SQ(Uxminus)/Hxminus + ghalf*SQ(Hxminus) )
#define UNEWXRGFLUXPLUS   ( SQ(Uxplus) /Hxplus  + ghalf*SQ(Hxplus) )
#define VNEWXRGFLUXMINUS  ( Uxminus*Vxminus/Hxminus )
#define VNEWXRGFLUXPLUS   ( Uxplus *Vxplus /Hxplus )

#define HNEWYRGFLUXMINUS  ( Vyminus )
#define HNEWYRGFLUXPLUS   ( Vyplus )
#define UNEWYRGFLUXMINUS  ( Vyminus*Uyminus/Hyminus )
#define UNEWYRGFLUXPLUS   ( Vyplus *Uyplus /Hyplus )
#define VNEWYRGFLUXMINUS  ( SQ(Vyminus)/Hyminus + ghalf*SQ(Hyminus) )
#define VNEWYRGFLUXPLUS   ( SQ(Vyplus) /Hyplus  + ghalf*SQ(Hyplus) )

/*****************************************************************/

void State::calc_finite_difference_regular_cells(double deltaT)
{
   real_t   g     = 9.80;   // gravitational constant
   real_t   ghalf = 0.5*g;

   struct timespec tstart_cpu;
   cpu_timer_start(&tstart_cpu);

   struct timespec tstart_cpu_part;
   cpu_timer_start(&tstart_cpu_part);

   //printf("\nDEBUG finite diff\n"); 

   // We need to populate the ghost regions since the calc neighbors has just been
   // established for the mesh shortly before

   //static state_t *H_new, *U_new, *V_new;
   static state_t ***H_reg_lev, ***U_reg_lev, ***V_reg_lev;
   static int ***mask_reg_lev;

   vector<real_t> &lev_deltax = mesh->lev_deltax;
   vector<real_t> &lev_deltay = mesh->lev_deltay;

   //static state ***varH, ***varU, ***varV, ****states_new; 
   static state_t ****states_new;
   static state_t ***Hxfluxplus, ***Hxfluxminus, ***Uxfluxplus, ***Uxfluxminus, ***Vxfluxplus, ***Vxfluxminus; 
   static state_t ***Hyfluxplus, ***Hyfluxminus, ***Uyfluxplus, ***Uyfluxminus, ***Vyfluxplus, ***Vyfluxminus; 
   static state_t ***wplusx_H, ***wminusx_H, ***wplusy_H, ***wminusy_H, ***wplusx_U, ***wminusx_U, ***wplusy_V, ***wminusy_V; 
   //static int    ***passFlag;
   apply_boundary_conditions();

#ifdef _OPENMP
#pragma omp barrier
#pragma omp master
   {
#endif

      mesh->calc_face_list_wbidirmap_phantom(state_memory, deltaT);

      memory_reset_ptrs(); //reset the pointers H,U,V that were recently reallocated in wbidirmap call
      mesh->generate_regular_cell_meshes(state_memory);

      H_reg_lev = (state_t ***)malloc((mesh->levmx+1)*sizeof(state_t **));
      U_reg_lev = (state_t ***)malloc((mesh->levmx+1)*sizeof(state_t **));
      V_reg_lev = (state_t ***)malloc((mesh->levmx+1)*sizeof(state_t **));
      mask_reg_lev = (int ***)malloc((mesh->levmx+1)*sizeof(int **));

      Hxfluxplus = (state_t ***)malloc((mesh->levmx+1)*sizeof(state_t **));
      Hxfluxminus = (state_t ***)malloc((mesh->levmx+1)*sizeof(state_t **));
      Hyfluxplus = (state_t ***)malloc((mesh->levmx+1)*sizeof(state_t **));
      Hyfluxminus = (state_t ***)malloc((mesh->levmx+1)*sizeof(state_t **));
      wplusx_H = (state_t ***)malloc((mesh->levmx+1)*sizeof(state_t **));
      wminusx_H = (state_t ***)malloc((mesh->levmx+1)*sizeof(state_t **));
      wplusy_H = (state_t ***)malloc((mesh->levmx+1)*sizeof(state_t **));
      wminusy_H = (state_t ***)malloc((mesh->levmx+1)*sizeof(state_t **));
      Uxfluxplus = (state_t ***)malloc((mesh->levmx+1)*sizeof(state_t **));
      Uxfluxminus = (state_t ***)malloc((mesh->levmx+1)*sizeof(state_t **));
      Uyfluxplus = (state_t ***)malloc((mesh->levmx+1)*sizeof(state_t **));
      Uyfluxminus = (state_t ***)malloc((mesh->levmx+1)*sizeof(state_t **));
      wplusx_U = (state_t ***)malloc((mesh->levmx+1)*sizeof(state_t **));
      wminusx_U = (state_t ***)malloc((mesh->levmx+1)*sizeof(state_t **));
      Vxfluxplus = (state_t ***)malloc((mesh->levmx+1)*sizeof(state_t **));
      Vxfluxminus = (state_t ***)malloc((mesh->levmx+1)*sizeof(state_t **));
      Vyfluxplus = (state_t ***)malloc((mesh->levmx+1)*sizeof(state_t **));
      Vyfluxminus = (state_t ***)malloc((mesh->levmx+1)*sizeof(state_t **));
      wplusy_V = (state_t ***)malloc((mesh->levmx+1)*sizeof(state_t **));
      wminusy_V = (state_t ***)malloc((mesh->levmx+1)*sizeof(state_t **));
      //varH = (double ***) malloc((mesh->levmx+1) * sizeof(double **));
      //varU = (double ***) malloc((mesh->levmx+1) * sizeof(double **));
      //varV = (double ***) malloc((mesh->levmx+1) * sizeof(double **));
      //passFlag = (int ***) malloc((mesh->levmx+1) * sizeof(int **));
      states_new = (state_t ****)malloc((mesh->levmx+1)*sizeof(double ***));

      for (int lev = 0; lev < mesh->levmx + 1; lev++){
          state_t ***pstate = mesh->meshes[lev].pstate;
          H_reg_lev[lev] = pstate[0];
          U_reg_lev[lev] = pstate[1];
          V_reg_lev[lev] = pstate[2];
          mask_reg_lev[lev] = mesh->meshes[lev].mask;

    
          Hxfluxplus[lev] = (state_t **)genmatrix(mesh->lev_jregsize[lev],mesh->lev_iregsize[lev],sizeof(double));
          Hxfluxminus[lev] = (state_t **)genmatrix(mesh->lev_jregsize[lev],mesh->lev_iregsize[lev],sizeof(double));
          Hyfluxplus[lev] = (state_t **)genmatrix(mesh->lev_jregsize[lev],mesh->lev_iregsize[lev],sizeof(double));
          Hyfluxminus[lev] = (state_t **)genmatrix(mesh->lev_jregsize[lev],mesh->lev_iregsize[lev],sizeof(double));
          wplusx_H[lev] = (state_t **)genmatrix(mesh->lev_jregsize[lev],mesh->lev_iregsize[lev],sizeof(double));
          wminusx_H[lev] = (state_t **)genmatrix(mesh->lev_jregsize[lev],mesh->lev_iregsize[lev],sizeof(double));
          wplusy_H[lev] = (state_t **)genmatrix(mesh->lev_jregsize[lev],mesh->lev_iregsize[lev],sizeof(double));
          wminusy_H[lev] = (state_t **)genmatrix(mesh->lev_jregsize[lev],mesh->lev_iregsize[lev],sizeof(double));
          Uxfluxplus[lev] = (state_t **)genmatrix(mesh->lev_jregsize[lev],mesh->lev_iregsize[lev],sizeof(double));
          Uxfluxminus[lev] = (state_t **)genmatrix(mesh->lev_jregsize[lev],mesh->lev_iregsize[lev],sizeof(double));
          Uyfluxplus[lev] = (state_t **)genmatrix(mesh->lev_jregsize[lev],mesh->lev_iregsize[lev],sizeof(double));
          Uyfluxminus[lev] = (state_t **)genmatrix(mesh->lev_jregsize[lev],mesh->lev_iregsize[lev],sizeof(double));
          wplusx_U[lev] = (state_t **)genmatrix(mesh->lev_jregsize[lev],mesh->lev_iregsize[lev],sizeof(double));
          wminusx_U[lev] = (state_t **)genmatrix(mesh->lev_jregsize[lev],mesh->lev_iregsize[lev],sizeof(double));
          Vxfluxplus[lev] = (state_t **)genmatrix(mesh->lev_jregsize[lev],mesh->lev_iregsize[lev],sizeof(double));
          Vxfluxminus[lev] = (state_t **)genmatrix(mesh->lev_jregsize[lev],mesh->lev_iregsize[lev],sizeof(double));
          Vyfluxplus[lev] = (state_t **)genmatrix(mesh->lev_jregsize[lev],mesh->lev_iregsize[lev],sizeof(double));
          Vyfluxminus[lev] = (state_t **)genmatrix(mesh->lev_jregsize[lev],mesh->lev_iregsize[lev],sizeof(double));
          wplusy_V[lev] = (state_t **)genmatrix(mesh->lev_jregsize[lev],mesh->lev_iregsize[lev],sizeof(double));
          wminusy_V[lev] = (state_t **)genmatrix(mesh->lev_jregsize[lev],mesh->lev_iregsize[lev],sizeof(double));
          //varH[lev] = (double **)genmatrix(mesh->lev_jregsize[lev],mesh->lev_iregsize[lev]*2,sizeof(double));
          //varU[lev] = (double **)genmatrix(mesh->lev_jregsize[lev],mesh->lev_iregsize[lev]*2,sizeof(double));
          //varV[lev] = (double **)genmatrix(mesh->lev_jregsize[lev],mesh->lev_iregsize[lev]*2,sizeof(double));
          //passFlag[lev] = (int **)genmatrix(mesh->lev_jregsize[lev],mesh->lev_iregsize[lev]*4,sizeof(int));
/*
         for(int jj=0; jj<mesh->lev_jregsize[lev]; jj++){
             for(int ii=0; ii<mesh->lev_iregsize[lev]; ii++){
                 Hxfluxplus[lev][jj][ii] = 0.0;
                 Hxfluxminus[lev][jj][ii] = 0.0;
                 Hyfluxplus[lev][jj][ii] = 0.0;
                 Hyfluxminus[lev][jj][ii] = 0.0;
                 Uxfluxplus[lev][jj][ii] = 0.0;
                 Uxfluxminus[lev][jj][ii] = 0.0;
                 Uyfluxplus[lev][jj][ii] = 0.0;
                 Uyfluxminus[lev][jj][ii] = 0.0;
                 Vxfluxplus[lev][jj][ii] = 0.0;
                 Vxfluxminus[lev][jj][ii] = 0.0;
                 Vyfluxplus[lev][jj][ii] = 0.0;
                 Vyfluxminus[lev][jj][ii] = 0.0;
                 wplusx_H[lev][jj][ii] = 0.0;
                 wminusx_H[lev][jj][ii] = 0.0;
                 wplusy_H[lev][jj][ii] = 0.0;
                 wminusy_H[lev][jj][ii] = 0.0;
                 wplusx_U[lev][jj][ii] = 0.0;
                 wminusx_U[lev][jj][ii] = 0.0;
                 wplusy_V[lev][jj][ii] = 0.0;
                 wminusy_V[lev][jj][ii] = 0.0;
             }
         }
*/
      }
   cpu_timers[STATE_TIMER_FINITE_DIFFERENCE_PART1] += cpu_timer_stop(tstart_cpu_part);
   cpu_timer_start(&tstart_cpu_part);

   for (int ll = 0; ll < mesh->levmx+1; ll++){
      int iimax = mesh->lev_iregsize[ll];
      int jjmax = mesh->lev_jregsize[ll];
      states_new[ll] = (state_t ***)gentrimatrix(3, jjmax, iimax, sizeof(state_t));
   }
#ifdef _OPENMP
   }
#pragma omp barrier
#endif

   //static state_t **H_reg, **U_reg, **V_reg;
   //int **mask_reg;

   //for (int ll=mesh->levmx; ll>-1; ll--){
   for (int ll=0; ll<mesh->levmx+1; ll++){
      //state_t ***pstate = mesh->meshes[ll].pstate;
      int iimax = mesh->lev_iregsize[ll];
      int jjmax = mesh->lev_jregsize[ll];
      int jj, ii;

      //H_reg = H_reg_lev[ll];
      //U_reg = U_reg_lev[ll];
      //V_reg = V_reg_lev[ll];
      //mask_reg = mask_reg_lev[ll];
      real_t dx = lev_deltax[ll];
      real_t dy = lev_deltay[ll];
      //real_t Cx = deltaT/dx;
      //real_t Cy = deltaT/dy;

      //state_t **H_reg_new = (state_t **)genmatrix(jjmax, iimax, sizeof(state_t));
      //state_t **U_reg_new = (state_t **)genmatrix(jjmax, iimax, sizeof(state_t));
      //state_t **V_reg_new = (state_t **)genmatrix(jjmax, iimax, sizeof(state_t));
      //

#ifdef _OPENMP
#pragma omp for
#endif
      for(jj=2; jj<jjmax-2; jj++){
#if defined(__GNUC_MINOR__)
         state_t * H_reg_loc = states_new[ll][0][jj];
         state_t * H_reg_lev_loc = H_reg_lev[ll][jj];
         state_t * U_reg_lev_loc = U_reg_lev[ll][jj];
         state_t * V_reg_lev_loc = V_reg_lev[ll][jj];
         state_t * Hxfluxminus_loc = Hxfluxminus[ll][jj];
         state_t * Uxfluxminus_loc = Uxfluxminus[ll][jj];
         state_t * Vxfluxminus_loc = Vxfluxminus[ll][jj];
         state_t * Hxfluxplus_loc = Hxfluxplus[ll][jj];
         state_t * Uxfluxplus_loc = Uxfluxplus[ll][jj];
         state_t * Vxfluxplus_loc = Vxfluxplus[ll][jj];
         state_t * Hyfluxminus_loc = Hyfluxminus[ll][jj];
         state_t * Uyfluxminus_loc = Uyfluxminus[ll][jj];
         state_t * Vyfluxminus_loc = Vyfluxminus[ll][jj];
         state_t * Hyfluxplus_loc = Hyfluxplus[ll][jj];
         state_t * Uyfluxplus_loc = Uyfluxplus[ll][jj];
         state_t * Vyfluxplus_loc = Vyfluxplus[ll][jj];
         state_t * H_reg_lev_minus_loc = H_reg_lev[ll][jj-1];
         state_t * H_reg_lev_minus2_loc = H_reg_lev[ll][jj-2];
         state_t * U_reg_lev_minus_loc = U_reg_lev[ll][jj-1];
         state_t * U_reg_lev_minus2_loc = U_reg_lev[ll][jj-2];
         state_t * V_reg_lev_minus_loc = V_reg_lev[ll][jj-1];
         state_t * V_reg_lev_minus2_loc = V_reg_lev[ll][jj-2];
         state_t * H_reg_lev_plus_loc = H_reg_lev[ll][jj+1];
         state_t * H_reg_lev_plus2_loc = H_reg_lev[ll][jj+2];
         state_t * U_reg_lev_plus_loc = U_reg_lev[ll][jj+1];
         state_t * U_reg_lev_plus2_loc = U_reg_lev[ll][jj+2];
         state_t * V_reg_lev_plus_loc = V_reg_lev[ll][jj+1];
         state_t * V_reg_lev_plus2_loc = V_reg_lev[ll][jj+2];
         state_t * wminusx_H_loc = wminusx_H[ll][jj];
         state_t * wplusx_H_loc = wplusx_H[ll][jj];
         state_t * wminusx_U_loc = wminusx_U[ll][jj];
         state_t * wplusx_U_loc = wplusx_U[ll][jj];
         state_t * wminusy_H_loc = wminusy_H[ll][jj];
         state_t * wplusy_H_loc = wplusy_H[ll][jj];
         state_t * wminusy_V_loc = wminusy_V[ll][jj];
         state_t * wplusy_V_loc = wplusy_V[ll][jj];
#endif
#ifdef _OPENMP_SIMD
#pragma omp simd
#endif
         for(ii=2; ii<iimax-2; ii++){
            if (mask_reg_lev[ll][jj][ii] == 1) {

            real_t Hxminus = HALF*( ((H_reg_lev[ll][jj][ii-1]) + (H_reg_lev[ll][jj][ii])) - (deltaT)/(dx)*((HXRGFLUXIC) - (HXRGFLUXNL)) );
            real_t Uxminus = HALF*( ((U_reg_lev[ll][jj][ii-1]) + (U_reg_lev[ll][jj][ii])) - (deltaT)/(dx)*((UXRGFLUXIC) - (UXRGFLUXNL)) );
            real_t Vxminus = HALF*( ((V_reg_lev[ll][jj][ii-1]) + (V_reg_lev[ll][jj][ii])) - (deltaT)/(dx)*((VXRGFLUXIC) - (VXRGFLUXNL)) );

            real_t Hxplus  = HALF*( ((H_reg_lev[ll][jj][ii]) + (H_reg_lev[ll][jj][ii+1])) - (deltaT)/(dx)*((HXRGFLUXNR) - (HXRGFLUXIC)) );
            real_t Uxplus  = HALF*( ((U_reg_lev[ll][jj][ii]) + (U_reg_lev[ll][jj][ii+1])) - (deltaT)/(dx)*((UXRGFLUXNR) - (UXRGFLUXIC)) );
            real_t Vxplus  = HALF*( ((V_reg_lev[ll][jj][ii]) + (V_reg_lev[ll][jj][ii+1])) - (deltaT)/(dx)*((VXRGFLUXNR) - (VXRGFLUXIC)) );

            real_t Hyminus = HALF*( ((H_reg_lev[ll][jj-1][ii]) + (H_reg_lev[ll][jj][ii])) - (deltaT)/(dy)*((HYRGFLUXIC) - (HYRGFLUXNB)) );
            real_t Uyminus = HALF*( ((U_reg_lev[ll][jj-1][ii]) + (U_reg_lev[ll][jj][ii])) - (deltaT)/(dy)*((UYRGFLUXIC) - (UYRGFLUXNB)) );
            real_t Vyminus = HALF*( ((V_reg_lev[ll][jj-1][ii]) + (V_reg_lev[ll][jj][ii])) - (deltaT)/(dy)*((VYRGFLUXIC) - (VYRGFLUXNB)) );

            real_t Hyplus  = HALF*( ((H_reg_lev[ll][jj][ii]) + (H_reg_lev[ll][jj+1][ii])) - (deltaT)/(dy)*((HYRGFLUXNT) - (HYRGFLUXIC)) );
            real_t Uyplus  = HALF*( ((U_reg_lev[ll][jj][ii]) + (U_reg_lev[ll][jj+1][ii])) - (deltaT)/(dy)*((UYRGFLUXNT) - (UYRGFLUXIC)) );
            real_t Vyplus  = HALF*( ((V_reg_lev[ll][jj][ii]) + (V_reg_lev[ll][jj+1][ii])) - (deltaT)/(dy)*((VYRGFLUXNT) - (VYRGFLUXIC)) );

            Hxfluxminus[ll][jj][ii] = HNEWXRGFLUXMINUS;
            Uxfluxminus[ll][jj][ii] = UNEWXRGFLUXMINUS;
            Vxfluxminus[ll][jj][ii] = VNEWXRGFLUXMINUS;

            Hxfluxplus[ll][jj][ii] = HNEWXRGFLUXPLUS;
            Uxfluxplus[ll][jj][ii] = UNEWXRGFLUXPLUS;
            Vxfluxplus[ll][jj][ii] = VNEWXRGFLUXPLUS;

            Hyfluxminus[ll][jj][ii] = HNEWYRGFLUXMINUS;
            Uyfluxminus[ll][jj][ii] = UNEWYRGFLUXMINUS;
            Vyfluxminus[ll][jj][ii] = VNEWYRGFLUXMINUS;

            Hyfluxplus[ll][jj][ii] = HNEWYRGFLUXPLUS;
            Uyfluxplus[ll][jj][ii] = UNEWYRGFLUXPLUS;
            Vyfluxplus[ll][jj][ii] = VNEWYRGFLUXPLUS;

            wminusx_H[ll][jj][ii] = w_corrector(deltaT, dx, fabs(Uxminus/Hxminus) + sqrt(g*Hxminus),
                              H_reg_lev[ll][jj][ii]-H_reg_lev[ll][jj][ii-1], H_reg_lev[ll][jj][ii-1]-H_reg_lev[ll][jj][ii-2], H_reg_lev[ll][jj][ii+1]-H_reg_lev[ll][jj][ii]);

            wminusx_H[ll][jj][ii] *= H_reg_lev[ll][jj][ii] - H_reg_lev[ll][jj][ii-1];

            wplusx_H[ll][jj][ii] = w_corrector(deltaT, dx, fabs(Uxplus/Hxplus) + sqrt(g*Hxplus),
                              H_reg_lev[ll][jj][ii+1]-H_reg_lev[ll][jj][ii], H_reg_lev[ll][jj][ii]-H_reg_lev[ll][jj][ii-1], H_reg_lev[ll][jj][ii+2]-H_reg_lev[ll][jj][ii+1]);

            wplusx_H[ll][jj][ii] *= H_reg_lev[ll][jj][ii+1] - H_reg_lev[ll][jj][ii];

            wminusx_U[ll][jj][ii] = w_corrector(deltaT, dx, fabs(Uxminus/Hxminus) + sqrt(g*Hxminus),
                              U_reg_lev[ll][jj][ii]-U_reg_lev[ll][jj][ii-1], U_reg_lev[ll][jj][ii-1]-U_reg_lev[ll][jj][ii-2], U_reg_lev[ll][jj][ii+1]-U_reg_lev[ll][jj][ii]);

            wminusx_U[ll][jj][ii] *= U_reg_lev[ll][jj][ii] - U_reg_lev[ll][jj][ii-1];

            wplusx_U[ll][jj][ii] = w_corrector(deltaT, dx, fabs(Uxplus/Hxplus) + sqrt(g*Hxplus),
                              U_reg_lev[ll][jj][ii+1]-U_reg_lev[ll][jj][ii], U_reg_lev[ll][jj][ii]-U_reg_lev[ll][jj][ii-1], U_reg_lev[ll][jj][ii+2]-U_reg_lev[ll][jj][ii+1]);

            wplusx_U[ll][jj][ii] *= U_reg_lev[ll][jj][ii+1] - U_reg_lev[ll][jj][ii];

            wminusy_H[ll][jj][ii] = w_corrector(deltaT, dy, fabs(Vyminus/Hyminus) + sqrt(g*Hyminus),
                              H_reg_lev[ll][jj][ii]-H_reg_lev[ll][jj-1][ii], H_reg_lev[ll][jj-1][ii]-H_reg_lev[ll][jj-2][ii], H_reg_lev[ll][jj+1][ii]-H_reg_lev[ll][jj][ii]);

            wminusy_H[ll][jj][ii] *= H_reg_lev[ll][jj][ii] - H_reg_lev[ll][jj-1][ii];

            wplusy_H[ll][jj][ii] = w_corrector(deltaT, dy, fabs(Vyplus/Hyplus) + sqrt(g*Hyplus),
                              H_reg_lev[ll][jj+1][ii]-H_reg_lev[ll][jj][ii], H_reg_lev[ll][jj][ii]-H_reg_lev[ll][jj-1][ii], H_reg_lev[ll][jj+2][ii]-H_reg_lev[ll][jj+1][ii]);

            wplusy_H[ll][jj][ii] *= H_reg_lev[ll][jj+1][ii] - H_reg_lev[ll][jj][ii];

            wminusy_V[ll][jj][ii] = w_corrector(deltaT, dy, fabs(Vyminus/Hyminus) + sqrt(g*Hyminus),
                              V_reg_lev[ll][jj][ii]-V_reg_lev[ll][jj-1][ii], V_reg_lev[ll][jj-1][ii]-V_reg_lev[ll][jj-2][ii], V_reg_lev[ll][jj+1][ii]-V_reg_lev[ll][jj][ii]);

            wminusy_V[ll][jj][ii] *= V_reg_lev[ll][jj][ii] - V_reg_lev[ll][jj-1][ii];

            wplusy_V[ll][jj][ii] = w_corrector(deltaT, dy, fabs(Vyplus/Hyplus) + sqrt(g*Hyplus),
                              V_reg_lev[ll][jj+1][ii]-V_reg_lev[ll][jj][ii], V_reg_lev[ll][jj][ii]-V_reg_lev[ll][jj-1][ii], V_reg_lev[ll][jj+2][ii]-V_reg_lev[ll][jj+1][ii]);

            wplusy_V[ll][jj][ii] *= V_reg_lev[ll][jj+1][ii] - V_reg_lev[ll][jj][ii];

/*
            if (passFlag[ll][jj][ii*4+0] == 1) {
               Hxfluxminus = 0.0;
               Uxfluxminus = 0.0;
               Vxfluxminus = 0.0;
                wminusx_H = 0.0;
                wminusx_U = 0.0;
            }
            if (passFlag[ll][jj][ii*4+1] == 1) {
               Hxfluxplus = 0.0;
               Uxfluxplus = 0.0;
               Vxfluxplus = 0.0;
                wplusx_H = 0.0;
                wplusx_U = 0.0;
            }
            if (passFlag[ll][jj][ii*4+2] == 1) {
               Hyfluxminus = 0.0;
               Uyfluxminus = 0.0;
               Vyfluxminus = 0.0;
                wminusy_H = 0.0;
                wminusy_V = 0.0;
            }
            if (passFlag[ll][jj][ii*4+3] == 1) {
               Hyfluxplus = 0.0;
               Uyfluxplus = 0.0;
               Vyfluxplus = 0.0;
                wplusy_H = 0.0;
                wplusy_V = 0.0;
            }
#ifdef _OPENMP
#pragma omp critical
{
#endif
            if ((mesh->phantomXFluxRG[ll][jj][ii] >= 0) && (mesh->phantomXFluxRG[ll][jj][ii] < mesh->ncells)) {
               int recvic = mesh->phantomXFluxRG[ll][jj][ii];
               int recvJ = mesh->j[recvic] - mesh->lev_jregmin[ll-1];
               int recvI = mesh->i[recvic] - mesh->lev_iregmin[ll-1];
               varH[ll-1][recvJ][recvI*2+0] += Hxfluxminus * HALF;
               varU[ll-1][recvJ][recvI*2+0] += Uxfluxminus * HALF;
               varV[ll-1][recvJ][recvI*2+0] += Vxfluxminus * HALF;
               varH[ll-1][recvJ][recvI*2+1] += wminusx_H / 4;
               varU[ll-1][recvJ][recvI*2+1] += wminusx_U / 4;
               passFlag[ll-1][recvJ][recvI*4+1] = 1;
           }
           else if ((mesh->phantomXFluxRG[ll][jj][ii] < 0) && (abs(mesh->phantomXFluxRG[ll][jj][ii]) < mesh->ncells)) {
               int recvic = abs(mesh->phantomXFluxRG[ll][jj][ii]);
               int recvJ = mesh->j[recvic] - mesh->lev_jregmin[ll-1];
               int recvI = mesh->i[recvic] - mesh->lev_iregmin[ll-1];
               varH[ll-1][recvJ][recvI*2+0] -= Hxfluxplus * HALF;
               varU[ll-1][recvJ][recvI*2+0] -= Uxfluxplus * HALF;
               varV[ll-1][recvJ][recvI*2+0] -= Vxfluxplus * HALF;
               varH[ll-1][recvJ][recvI*2+1] -= wplusx_H / 4;
               varU[ll-1][recvJ][recvI*2+1] -= wplusx_U / 4;
               passFlag[ll-1][recvJ][recvI*4+0] = 1;
           }
           if ((mesh->phantomYFluxRG[ll][jj][ii] >= 0) && (mesh->phantomYFluxRG[ll][jj][ii] < mesh->ncells)) {
               int recvic = mesh->phantomYFluxRG[ll][jj][ii];
               int recvJ = mesh->j[recvic] - mesh->lev_jregmin[ll-1];
               int recvI = mesh->i[recvic] - mesh->lev_iregmin[ll-1];
               varH[ll-1][recvJ][recvI*2+0] += Hyfluxminus * HALF;
               varU[ll-1][recvJ][recvI*2+0] += Uyfluxminus * HALF;
               varV[ll-1][recvJ][recvI*2+0] += Vyfluxminus * HALF;
               varH[ll-1][recvJ][recvI*2+1] += wminusy_H / 4;
               varV[ll-1][recvJ][recvI*2+1] += wminusy_V / 4;
               passFlag[ll-1][recvJ][recvI*4+3] = 1;
           }
           else if ((mesh->phantomYFluxRG[ll][jj][ii] < 0) && (abs(mesh->phantomYFluxRG[ll][jj][ii]) < mesh->ncells)) {
               int recvic = abs(mesh->phantomYFluxRG[ll][jj][ii]);
               int recvJ = mesh->j[recvic] - mesh->lev_jregmin[ll-1];
               int recvI = mesh->i[recvic] - mesh->lev_iregmin[ll-1];
               varH[ll-1][recvJ][recvI*2+0] -= Hyfluxplus * HALF;
               varU[ll-1][recvJ][recvI*2+0] -= Uyfluxplus * HALF;
               varV[ll-1][recvJ][recvI*2+0] -= Vyfluxplus * HALF;
               varH[ll-1][recvJ][recvI*2+1] -= wplusy_H / 4;
               varV[ll-1][recvJ][recvI*2+1] -= wplusy_V / 4;
               passFlag[ll-1][recvJ][recvI*4+2] = 1;
           }

           states_new[ll][0][jj][ii] = U_fullstep(deltaT, dx, H_reg_lev[ll][jj][ii],
           //H_reg_new[jj][ii] = U_fullstep(deltaT, dx, H_reg[jj][ii],
                       Hxfluxplus, Hxfluxminus, Hyfluxplus + varH[ll][jj][ii*2], Hyfluxminus)
                  - wminusx_H + wplusx_H - wminusy_H + wplusy_H + varH[ll][jj][ii*2+1];

           states_new[ll][1][jj][ii] = U_fullstep(deltaT, dx, U_reg_lev[ll][jj][ii],
           //U_reg_new[jj][ii] = U_fullstep(deltaT, dx, U_reg[jj][ii],
                       Uxfluxplus, Uxfluxminus, Uyfluxplus + varU[ll][jj][ii*2], Uyfluxminus)
                  - wminusx_U + wplusx_U + varU[ll][jj][ii*2+1];

           states_new[ll][2][jj][ii] = U_fullstep(deltaT, dx, V_reg_lev[ll][jj][ii],
           //V_reg_new[jj][ii] = U_fullstep(deltaT, dx, V_reg[jj][ii],
                       Vxfluxplus, Vxfluxminus, Vyfluxplus + varV[ll][jj][ii*2], Vyfluxminus)
                  - wminusy_V + wplusy_V + varV[ll][jj][ii*2+1];
#ifdef _OPENMP
}
#endif
*/
         }
         } // ii
      } // jj 


      /*for(int jj=2; jj<jjmax-2; jj++){
         for(int ii=2; ii<iimax-2; ii++){
             H_reg[jj][ii] = H_reg_new[jj][ii];
             U_reg[jj][ii] = U_reg_new[jj][ii];
             V_reg[jj][ii] = V_reg_new[jj][ii];
         }
      }*/
   } // ll


#ifdef _OPENMP
#pragma omp barrier
#pragma omp master
   {
#endif
      cpu_timers[STATE_TIMER_FINITE_DIFFERENCE_PART2] += cpu_timer_stop(tstart_cpu_part);
      cpu_timer_start(&tstart_cpu_part);
#ifdef _OPENMP
   }
#pragma omp barrier
#endif


#ifdef _OPENMP
#ifdef _OPENMP_SIMD
#pragma omp for simd
#else
#pragma omp for
#endif
#else
#ifdef _OPENMP_SIMD
#pragma omp simd
#endif
#endif
   for (int ifix = 0; ifix < mesh->nxfixup; ifix++){
      int ic = mesh->xrecvCIdx[ifix];
      int llr = mesh->level[ic];
      int icjj = mesh->j[ic] - mesh->lev_jregmin[llr];
      int icii = mesh->i[ic] - mesh->lev_iregmin[llr];
      if (mesh->xplusCell2Idx[ic] > -1) {
         int ifixup = mesh->xplusCell2Idx[ic];
         int ns1 = mesh->map_xface2cell_upper[mesh->xsendIdx1[ifixup]];
         int ns2 = mesh->map_xface2cell_upper[mesh->xsendIdx2[ifixup]];

         int ll = mesh->level[ns1];
         int ns1jj = mesh->j[ns1] - mesh->lev_jregmin[ll];
         int ns1ii = mesh->i[ns1] - mesh->lev_iregmin[ll];
         int ns2jj = mesh->j[ns2] - mesh->lev_jregmin[ll];
         int ns2ii = mesh->i[ns2] - mesh->lev_iregmin[ll];


         Hxfluxplus[llr][icjj][icii] = (Hxfluxminus[ll][ns1jj][ns1ii] + Hxfluxminus[ll][ns2jj][ns2ii]) * HALF;
         Uxfluxplus[llr][icjj][icii] = (Uxfluxminus[ll][ns1jj][ns1ii] + Uxfluxminus[ll][ns2jj][ns2ii]) * HALF;
         Vxfluxplus[llr][icjj][icii] = (Vxfluxminus[ll][ns1jj][ns1ii] + Vxfluxminus[ll][ns2jj][ns2ii]) * HALF;
         wplusx_H[llr][icjj][icii] = (wminusx_H[ll][ns1jj][ns1ii] + wminusx_H[ll][ns2jj][ns2ii]) * 0.25;
         wplusx_U[llr][icjj][icii] = (wminusx_U[ll][ns1jj][ns1ii] + wminusx_U[ll][ns2jj][ns2ii]) * 0.25;
         //Hxfluxplus[ic] = (Hxfluxminus[ns1] + Hxfluxminus[ns2]) * HALF;
         //Uxfluxplus[ic] = (Uxfluxminus[ns1] + Uxfluxminus[ns2]) * HALF;
         //Vxfluxplus[ic] = (Vxfluxminus[ns1] + Vxfluxminus[ns2]) * HALF;
         //wplusx_H[ic] = (wminusx_H[ns1] + wminusx_H[ns2]) * 0.25;
         //wplusx_U[ic] = (wminusx_U[ns1] + wminusx_U[ns2]) * 0.25;
      }

      if (mesh->xminusCell2Idx[ic] > -1) {
         int ifixup = mesh->xminusCell2Idx[ic];
         int ns1 = mesh->map_xface2cell_lower[mesh->xsendIdx1[ifixup]];
         int ns2 = mesh->map_xface2cell_lower[mesh->xsendIdx2[ifixup]];
         
         int ll = mesh->level[ns1];
         int ns1jj = mesh->j[ns1] - mesh->lev_jregmin[ll];
         int ns1ii = mesh->i[ns1] - mesh->lev_iregmin[ll];
         int ns2jj = mesh->j[ns2] - mesh->lev_jregmin[ll];
         int ns2ii = mesh->i[ns2] - mesh->lev_iregmin[ll];

         Hxfluxminus[llr][icjj][icii] = (Hxfluxplus[ll][ns1jj][ns1ii] + Hxfluxplus[ll][ns2jj][ns2ii]) * HALF;
         Uxfluxminus[llr][icjj][icii] = (Uxfluxplus[ll][ns1jj][ns1ii] + Uxfluxplus[ll][ns2jj][ns2ii]) * HALF;
         Vxfluxminus[llr][icjj][icii] = (Vxfluxplus[ll][ns1jj][ns1ii] + Vxfluxplus[ll][ns2jj][ns2ii]) * HALF;
         wminusx_H[llr][icjj][icii] = (wplusx_H[ll][ns1jj][ns1ii] + wplusx_H[ll][ns2jj][ns2ii]) * 0.25;
         wminusx_U[llr][icjj][icii] = (wplusx_U[ll][ns1jj][ns1ii] + wplusx_U[ll][ns2jj][ns2ii]) * 0.25;
         //Hxfluxminus[ic] = (Hxfluxplus[ns1] + Hxfluxplus[ns2]) * HALF;
         //Uxfluxminus[ic] = (Uxfluxplus[ns1] + Uxfluxplus[ns2]) * HALF;
         //Vxfluxminus[ic] = (Vxfluxplus[ns1] + Vxfluxplus[ns2]) * HALF;
         //wminusx_H[ic] = (wplusx_H[ns1] + wplusx_H[ns2]) * 0.25;
         //wminusx_U[ic] = (wplusx_U[ns1] + wplusx_U[ns2]) * 0.25;
      }
   }

#ifdef _OPENMP
#ifdef _OPENMP_SIMD
#pragma omp for simd
#else
#pragma omp for
#endif
#else
#ifdef _OPENMP_SIMD
#pragma omp simd
#endif
#endif
   for (int ifix = 0; ifix < mesh->nyfixup; ifix++){
      int ic = mesh->yrecvCIdx[ifix];
      int llr = mesh->level[ic];
      int icjj = mesh->j[ic] - mesh->lev_jregmin[llr];
      int icii = mesh->i[ic] - mesh->lev_iregmin[llr];
      if (mesh->yplusCell2Idx[ic] > -1) {
         int ifixup = mesh->yplusCell2Idx[ic];
         int ns1 = mesh->map_yface2cell_upper[mesh->ysendIdx1[ifixup]];
         int ns2 = mesh->map_yface2cell_upper[mesh->ysendIdx2[ifixup]];

         int ll = mesh->level[ns1];
         int ns1jj = mesh->j[ns1] - mesh->lev_jregmin[ll];
         int ns1ii = mesh->i[ns1] - mesh->lev_iregmin[ll];
         int ns2jj = mesh->j[ns2] - mesh->lev_jregmin[ll];
         int ns2ii = mesh->i[ns2] - mesh->lev_iregmin[ll];

         Hyfluxplus[llr][icjj][icii] = (Hyfluxminus[ll][ns1jj][ns1ii] + Hyfluxminus[ll][ns2jj][ns2ii]) * HALF;
         Uyfluxplus[llr][icjj][icii] = (Uyfluxminus[ll][ns1jj][ns1ii] + Uyfluxminus[ll][ns2jj][ns2ii]) * HALF;
         Vyfluxplus[llr][icjj][icii] = (Vyfluxminus[ll][ns1jj][ns1ii] + Vyfluxminus[ll][ns2jj][ns2ii]) * HALF;
         wplusy_H[llr][icjj][icii] = (wminusy_H[ll][ns1jj][ns1ii] + wminusy_H[ll][ns2jj][ns2ii]) * 0.25;
         wplusy_V[llr][icjj][icii] = (wminusy_V[ll][ns1jj][ns1ii] + wminusy_V[ll][ns2jj][ns2ii]) * 0.25;
         //Hyfluxplus[ic] = (Hyfluxminus[ns1] + Hyfluxminus[ns2]) * HALF;
         //Uyfluxplus[ic] = (Uyfluxminus[ns1] + Uyfluxminus[ns2]) * HALF;
         //Vyfluxplus[ic] = (Vyfluxminus[ns1] + Vyfluxminus[ns2]) * HALF;
         //wplusy_H[ic] = (wminusy_H[ns1] + wminusy_H[ns2]) * 0.25;
         //wplusy_V[ic] = (wminusy_V[ns1] + wminusy_V[ns2]) * 0.25;
      }

      if (mesh->yminusCell2Idx[ic] > -1) {
         int ifixup = mesh->yminusCell2Idx[ic];
         int ns1 = mesh->map_yface2cell_lower[mesh->ysendIdx1[ifixup]];
         int ns2 = mesh->map_yface2cell_lower[mesh->ysendIdx2[ifixup]];

         int ll = mesh->level[ns1];
         int ns1jj = mesh->j[ns1] - mesh->lev_jregmin[ll];
         int ns1ii = mesh->i[ns1] - mesh->lev_iregmin[ll];
         int ns2jj = mesh->j[ns2] - mesh->lev_jregmin[ll];
         int ns2ii = mesh->i[ns2] - mesh->lev_iregmin[ll];

         Hyfluxminus[llr][icjj][icii] = (Hyfluxplus[ll][ns1jj][ns1ii] + Hyfluxplus[ll][ns2jj][ns2ii]) * HALF;
         Uyfluxminus[llr][icjj][icii] = (Uyfluxplus[ll][ns1jj][ns1ii] + Uyfluxplus[ll][ns2jj][ns2ii]) * HALF;
         Vyfluxminus[llr][icjj][icii] = (Vyfluxplus[ll][ns1jj][ns1ii] + Vyfluxplus[ll][ns2jj][ns2ii]) * HALF;
         wminusy_H[llr][icjj][icii] = (wplusy_H[ll][ns1jj][ns1ii] + wplusy_H[ll][ns2jj][ns2ii]) * 0.25;
         wminusy_V[llr][icjj][icii] = (wplusy_V[ll][ns1jj][ns1ii] + wplusy_V[ll][ns2jj][ns2ii]) * 0.25;
         //Hyfluxminus[ic] = (Hyfluxplus[ns1] + Hyfluxplus[ns2]) * HALF;
         //Uyfluxminus[ic] = (Uyfluxplus[ns1] + Uyfluxplus[ns2]) * HALF;
         //Vyfluxminus[ic] = (Vyfluxplus[ns1] + Vyfluxplus[ns2]) * HALF;
         //wminusy_H[ic] = (wplusy_H[ns1] + wplusy_H[ns2]) * 0.25;
         //wminusy_V[ic] = (wplusy_V[ns1] + wplusy_V[ns2]) * 0.25;
      }
   }

#ifdef _OPENMP
#pragma omp barrier
#pragma omp master
   {
#endif
      cpu_timers[STATE_TIMER_FINITE_DIFFERENCE_PART3] += cpu_timer_stop(tstart_cpu_part);
      cpu_timer_start(&tstart_cpu_part);

#ifdef _OPENMP
   }
#pragma omp barrier
#endif



   for (int ll=0; ll<mesh->levmx+1; ll++){
      int iimax = mesh->lev_iregsize[ll];
      int jjmax = mesh->lev_jregsize[ll];
      int jj, ii;
      real_t dx = lev_deltax[ll];
      //real_t dy = lev_deltay[ll];

#ifdef _OPENMP
#pragma omp for
#endif
      for(jj=2; jj<jjmax-2; jj++){
#if defined(__GNUC_MINOR__)
         state_t *H_reg_new_loc = states_new[ll][0][jj];
         state_t *U_reg_new_loc = states_new[ll][1][jj];
         state_t *V_reg_new_loc = states_new[ll][2][jj];
         state_t *H_reg_loc = H_reg_lev[ll][jj];
         state_t *U_reg_loc = U_reg_lev[ll][jj];
         state_t *V_reg_loc = V_reg_lev[ll][jj];
         state_t *Hxfluxplus_loc = Hxfluxplus[ll][jj];
         state_t *Hxfluxminus_loc = Hxfluxminus[ll][jj];
         state_t *Hyfluxplus_loc = Hyfluxplus[ll][jj];
         state_t *Hyfluxminus_loc = Hyfluxminus[ll][jj];
         state_t *Uxfluxplus_loc = Uxfluxplus[ll][jj];
         state_t *Uxfluxminus_loc = Uxfluxminus[ll][jj];
         state_t *Uyfluxplus_loc = Uyfluxplus[ll][jj];
         state_t *Uyfluxminus_loc = Uyfluxminus[ll][jj];
         state_t *Vxfluxplus_loc = Vxfluxplus[ll][jj];
         state_t *Vxfluxminus_loc = Vxfluxminus[ll][jj];
         state_t *Vyfluxplus_loc = Vyfluxplus[ll][jj];
         state_t *Vyfluxminus_loc = Vyfluxminus[ll][jj];
         state_t *wplusx_H_loc = wplusx_H[ll][jj];
         state_t *wminusx_H_loc = wminusx_H[ll][jj];
         state_t *wplusy_H_loc = wplusy_H[ll][jj];
         state_t *wminusy_H_loc = wminusy_H[ll][jj];
         state_t *wplusx_U_loc = wplusx_U[ll][jj];
         state_t *wminusx_U_loc = wminusx_U[ll][jj];
         state_t *wplusy_V_loc = wplusy_V[ll][jj];
         state_t *wminusy_V_loc = wminusy_V[ll][jj];
#endif
#ifdef _OPENMP_SIMD
#pragma omp simd
#endif
         for(ii=2; ii<iimax-2; ii++){
            if (mask_reg_lev[ll][jj][ii] == 1) {
#ifdef PRECISION_CHECK_WITH_PARENTHESIS
           //Basic parentheses
           states_new[ll][0][jj][ii] = U_fullstep(deltaT, dx, H_reg_lev[ll][jj][ii],
                       Hxfluxplus[ll][jj][ii], Hxfluxminus[ll][jj][ii], Hyfluxplus[ll][jj][ii], Hyfluxminus[ll][jj][ii])
                  + (- wminusx_H[ll][jj][ii] + wplusx_H[ll][jj][ii] - wminusy_H[ll][jj][ii] + wplusy_H[ll][jj][ii]);
#else
#ifdef PRECISION_CHECK_BEST_PARENTHESIS
           //Best parentheses version
           states_new[ll][0][jj][ii] = U_fullstep(deltaT, dx, H_reg_lev[ll][jj][ii],
                       Hxfluxplus[ll][jj][ii], Hxfluxminus[ll][jj][ii], Hyfluxplus[ll][jj][ii], Hyfluxminus[ll][jj][ii])
                  //+ ((( -wminusx_H[ll][jj][ii] - wminusy_H[ll][jj][ii] ) + wplusy_H[ll][jj][ii] ) + wplusx_H[ll][jj][ii] );
                  + (((wplusx_H[ll][jj][ii] + wplusy_H[ll][jj][ii]) - wminusy_H[ll][jj][ii]) - wminusx_H[ll][jj][ii]);
#else
           //Original version with no parentheses
           states_new[ll][0][jj][ii] = U_fullstep(deltaT, dx, H_reg_lev[ll][jj][ii],
                       Hxfluxplus[ll][jj][ii], Hxfluxminus[ll][jj][ii], Hyfluxplus[ll][jj][ii], Hyfluxminus[ll][jj][ii])
                  - wminusx_H[ll][jj][ii] + wplusx_H[ll][jj][ii] - wminusy_H[ll][jj][ii] + wplusy_H[ll][jj][ii];
#endif
#endif

           states_new[ll][1][jj][ii] = U_fullstep(deltaT, dx, U_reg_lev[ll][jj][ii],
           //U_reg_new[jj][ii] = U_fullstep(deltaT, dx, U_reg[jj][ii],
                       Uxfluxplus[ll][jj][ii], Uxfluxminus[ll][jj][ii], Uyfluxplus[ll][jj][ii], Uyfluxminus[ll][jj][ii])
                  - wminusx_U[ll][jj][ii] + wplusx_U[ll][jj][ii];

           states_new[ll][2][jj][ii] = U_fullstep(deltaT, dx, V_reg_lev[ll][jj][ii],
           //V_reg_new[jj][ii] = U_fullstep(deltaT, dx, V_reg[jj][ii],
                       Vxfluxplus[ll][jj][ii], Vxfluxminus[ll][jj][ii], Vyfluxplus[ll][jj][ii], Vyfluxminus[ll][jj][ii])
                  - wminusy_V[ll][jj][ii] + wplusy_V[ll][jj][ii];
            }
         } // ii
      } // jj
   } // ll

#ifdef _OPENMP
#pragma omp barrier
#pragma omp master
   {
#endif
      cpu_timers[STATE_TIMER_FINITE_DIFFERENCE_PART4] += cpu_timer_stop(tstart_cpu_part);
      cpu_timer_start(&tstart_cpu_part);

   for (int ll=0; ll<mesh->levmx+1; ll++){
      //genmatrixfree((void **)varH[ll]);
      //genmatrixfree((void **)varU[ll]);
      //genmatrixfree((void **)varV[ll]);
      //genmatrixfree((void **)passFlag[ll]);
      genmatrixfree((void **)Hxfluxplus[ll]);
      genmatrixfree((void **)Hxfluxminus[ll]);
      genmatrixfree((void **)Hyfluxplus[ll]);
      genmatrixfree((void **)Hyfluxminus[ll]);
      genmatrixfree((void **)Uxfluxplus[ll]);
      genmatrixfree((void **)Uxfluxminus[ll]);
      genmatrixfree((void **)Uyfluxplus[ll]);
      genmatrixfree((void **)Uyfluxminus[ll]);
      genmatrixfree((void **)Vxfluxplus[ll]);
      genmatrixfree((void **)Vxfluxminus[ll]);
      genmatrixfree((void **)Vyfluxplus[ll]);
      genmatrixfree((void **)Vyfluxminus[ll]);
      genmatrixfree((void **)wplusx_H[ll]);
      genmatrixfree((void **)wminusx_H[ll]);
      genmatrixfree((void **)wplusy_H[ll]);
      genmatrixfree((void **)wminusy_H[ll]);
      genmatrixfree((void **)wplusx_U[ll]);
      genmatrixfree((void **)wminusx_U[ll]);
      genmatrixfree((void **)wplusy_V[ll]);
      genmatrixfree((void **)wminusy_V[ll]);
      state_t ***tmp_reg;
      //printf("tmp %x pstate %x new %x\n", tmp_reg, mesh->meshes[ll].pstate, states_new);
      SWAP_PTR(states_new[ll], mesh->meshes[ll].pstate, tmp_reg);
      //printf("tmp %x pstate %x new %x\n", tmp_reg, mesh->meshes[ll].pstate, states_new);
      gentrimatrixfree((void ***) states_new[ll]);
      //genmatrixfree((void **)H_reg_new);
      //genmatrixfree((void **)U_reg_new);
      //genmatrixfree((void **)V_reg_new);
   } // ll

   //free(varH);
   //free(varU);
   //free(varV);
   //free(passFlag);
   free(Hxfluxplus);
   free(Hxfluxminus);
   free(Hyfluxplus);
   free(Hyfluxminus);
   free(Uxfluxplus);
   free(Uxfluxminus);
   free(Uyfluxplus);
   free(Uyfluxminus);
   free(Vxfluxplus);
   free(Vxfluxminus);
   free(Vyfluxplus);
   free(Vyfluxminus);
   free(wplusx_H);
   free(wminusx_H);
   free(wplusy_H);
   free(wminusy_H);
   free(wplusx_U);
   free(wminusx_U);
   free(wplusy_V);
   free(wminusy_V);
   free(H_reg_lev);
   free(U_reg_lev);
   free(V_reg_lev);
   free(mask_reg_lev);
   free(states_new);

   mesh->destroy_regular_cell_meshes(state_memory);

   cpu_timers[STATE_TIMER_FINITE_DIFFERENCE_PART5] += cpu_timer_stop(tstart_cpu_part);
   cpu_timer_start(&tstart_cpu_part);

   cpu_timers[STATE_TIMER_FINITE_DIFFERENCE] += cpu_timer_stop(tstart_cpu);
#ifdef _OPENMP
      }
#endif

#ifdef _OPENMP
#pragma omp barrier
#endif
}

void State::calc_finite_difference_regular_cells_by_faces(double deltaT)
{
   real_t   g     = 9.80;   // gravitational constant
   real_t   ghalf = 0.5*g;

   struct timespec tstart_cpu;
   cpu_timer_start(&tstart_cpu);

   struct timespec tstart_cpu_part;
   cpu_timer_start(&tstart_cpu_part);

   //printf("\nDEBUG finite diff\n"); 

   // We need to populate the ghost regions since the calc neighbors has just been
   // established for the mesh shortly before
   apply_boundary_conditions();

   //static state_t *H_new, *U_new, *V_new;
   static state_t ****states_new;
   static state_t ***H_reg_lev, ***U_reg_lev, ***V_reg_lev;
   static int ***mask_reg_lev;

   static real_t ***HxFlux, ***UxFlux, ***VxFlux, ***Wx_H, ***Wx_U;
   static real_t ***HyFlux, ***UyFlux, ***VyFlux, ***Wy_H, ***Wy_V;
   //int ***passFlagX, ***passFlagY;

#ifdef _OPENMP
#pragma omp barrier
#pragma omp master
   {
#endif
      mesh->calc_face_list_wbidirmap_phantom(state_memory, deltaT);
      memory_reset_ptrs(); //reset the pointers H,U,V that were recently reallocated in wbidirmap call
      mesh->generate_regular_cell_meshes(state_memory);
      H_reg_lev = (state_t ***)malloc((mesh->levmx+1)*sizeof(state_t **));
      U_reg_lev = (state_t ***)malloc((mesh->levmx+1)*sizeof(state_t **));
      V_reg_lev = (state_t ***)malloc((mesh->levmx+1)*sizeof(state_t **));
      mask_reg_lev = (int ***)malloc((mesh->levmx+1)*sizeof(int **));

      HxFlux = (real_t ***)malloc((mesh->levmx+1)*sizeof(real_t**));
      UxFlux = (real_t ***)malloc((mesh->levmx+1)*sizeof(real_t**));
      VxFlux = (real_t ***)malloc((mesh->levmx+1)*sizeof(real_t**));
      Wx_H   = (real_t ***)malloc((mesh->levmx+1)*sizeof(real_t**));
      Wx_U   = (real_t ***)malloc((mesh->levmx+1)*sizeof(real_t**));
      HyFlux = (real_t ***)malloc((mesh->levmx+1)*sizeof(real_t**));
      UyFlux = (real_t ***)malloc((mesh->levmx+1)*sizeof(real_t**));
      VyFlux = (real_t ***)malloc((mesh->levmx+1)*sizeof(real_t**));
      Wy_H   = (real_t ***)malloc((mesh->levmx+1)*sizeof(real_t**));
      Wy_V   = (real_t ***)malloc((mesh->levmx+1)*sizeof(real_t**));
      states_new = (state_t ****)malloc((mesh->levmx+1)*sizeof(state_t ***));
      //passFlagX   = (int ***)malloc((mesh->levmx+1)*sizeof(int**));
      //passFlagY   = (int ***)malloc((mesh->levmx+1)*sizeof(int**));

      for (int lev = 0; lev < mesh->levmx+1; lev++){
         state_t ***pstate = mesh->meshes[lev].pstate;
         H_reg_lev[lev] = pstate[0];
         U_reg_lev[lev] = pstate[1];
         V_reg_lev[lev] = pstate[2];
         mask_reg_lev[lev] = mesh->meshes[lev].mask;

         HxFlux[lev] = (real_t **)genmatrix(mesh->lev_jregsize[lev],mesh->lev_iregsize[lev],sizeof(real_t));
         UxFlux[lev] = (real_t **)genmatrix(mesh->lev_jregsize[lev],mesh->lev_iregsize[lev],sizeof(real_t));
         VxFlux[lev] = (real_t **)genmatrix(mesh->lev_jregsize[lev],mesh->lev_iregsize[lev],sizeof(real_t));
         Wx_H[lev]   = (real_t **)genmatrix(mesh->lev_jregsize[lev],mesh->lev_iregsize[lev],sizeof(real_t));
         Wx_U[lev]   = (real_t **)genmatrix(mesh->lev_jregsize[lev],mesh->lev_iregsize[lev],sizeof(real_t));
         HyFlux[lev] = (real_t **)genmatrix(mesh->lev_jregsize[lev],mesh->lev_iregsize[lev],sizeof(real_t));
         UyFlux[lev] = (real_t **)genmatrix(mesh->lev_jregsize[lev],mesh->lev_iregsize[lev],sizeof(real_t));
         VyFlux[lev] = (real_t **)genmatrix(mesh->lev_jregsize[lev],mesh->lev_iregsize[lev],sizeof(real_t));
         Wy_H[lev]   = (real_t **)genmatrix(mesh->lev_jregsize[lev],mesh->lev_iregsize[lev],sizeof(real_t));
         Wy_V[lev]   = (real_t **)genmatrix(mesh->lev_jregsize[lev],mesh->lev_iregsize[lev],sizeof(real_t));
         //passFlagX[lev]   = (int **)genmatrix(mesh->lev_jregsize[lev],mesh->lev_iregsize[lev],sizeof(int));
         //passFlagY[lev]   = (int **)genmatrix(mesh->lev_jregsize[lev],mesh->lev_iregsize[lev],sizeof(int));

      }


      for (int ll = 0; ll < mesh->levmx+1; ll++){
          int iimax = mesh->lev_iregsize[ll];
          int jjmax = mesh->lev_jregsize[ll];
          states_new[ll] = (state_t ***)gentrimatrix(3, jjmax, iimax, sizeof(state_t));
      }

      cpu_timers[STATE_TIMER_FINITE_DIFFERENCE_PART1] += cpu_timer_stop(tstart_cpu_part);
      cpu_timer_start(&tstart_cpu_part);
#ifdef _OPENMP
   }
#pragma omp barrier
#endif

/*
   for (int ll = 0; ll < mesh->levmx; ll++){
      for (int jj = 0; jj < mesh->lev_jregsize[ll]; jj++) {
         for (int ii = 0; ii < mesh->lev_iregsize[ll]; ii++) {
            passFlagX[ll][jj][ii] = 0;
            passFlagY[ll][jj][ii] = 0;
         }
      }
   }
*/
   for (int ll = 0; ll < mesh->levmx+1; ll++){
      //state_t ***pstate = mesh->meshes[ll].pstate;
      int iimax = mesh->lev_iregsize[ll];
      int jjmax = mesh->lev_jregsize[ll];

      //state_t **H_reg = H_reg_lev[ll];
      //state_t **U_reg = U_reg_lev[ll];
      //state_t **V_reg = V_reg_lev[ll];
      //int **mask_reg = mask_reg_lev[ll];
      real_t dx = mesh->lev_deltax[ll];
      real_t dy = mesh->lev_deltay[ll];
      real_t Cxhalf = deltaT/dx;
      real_t Cyhalf = deltaT/dy;

#ifdef _OPENMP
#pragma omp for
#endif
      for(int jj=2; jj<jjmax-1; jj++){
#if defined(__GNUC_MINOR__)
         real_t * HxFlux_loc = HxFlux[ll][jj];
         real_t * UxFlux_loc = UxFlux[ll][jj];
         real_t * VxFlux_loc = VxFlux[ll][jj];
         real_t * Wx_H_loc = Wx_H[ll][jj];
         real_t * Wx_U_loc = Wx_U[ll][jj];
         state_t * H_reg_lev_loc = H_reg_lev[ll][jj];
         state_t * U_reg_lev_loc = U_reg_lev[ll][jj];
         state_t * V_reg_lev_loc = V_reg_lev[ll][jj];
         state_t * H_reg_lev_plus_loc = H_reg_lev[ll][jj+1];
         state_t * U_reg_lev_plus_loc = U_reg_lev[ll][jj+1];
         state_t * V_reg_lev_plus_loc = V_reg_lev[ll][jj+1];
         state_t * H_reg_lev_minus_loc = H_reg_lev[ll][jj-1];
         state_t * U_reg_lev_minus_loc = U_reg_lev[ll][jj-1];
         state_t * V_reg_lev_minus_loc = V_reg_lev[ll][jj-1];
         state_t * H_reg_lev_minus2_loc = H_reg_lev[ll][jj-2];
         state_t * U_reg_lev_minus2_loc = U_reg_lev[ll][jj-2];
         state_t * V_reg_lev_minus2_loc = V_reg_lev[ll][jj-2];
         real_t * HyFlux_loc = HyFlux[ll][jj];
         real_t * UyFlux_loc = UyFlux[ll][jj];
         real_t * VyFlux_loc = VyFlux[ll][jj];
         real_t * Wy_H_loc = Wy_H[ll][jj];
         real_t * Wy_V_loc = Wy_V[ll][jj];
#endif
#ifdef _OPENMP_SIMD
#pragma omp simd
#endif
         for(int ii=2; ii<iimax-1; ii++){
            if ((mask_reg_lev[ll][jj][ii-1] == 1 || mask_reg_lev[ll][jj][ii] == 1)){
               real_t Hx = HALF*( (H_reg_lev[ll][jj][ii-1] + H_reg_lev[ll][jj][ii]) - Cxhalf*(HXRGFLUXIC - HXRGFLUXNL) );
               real_t Ux = HALF*( (U_reg_lev[ll][jj][ii-1] + U_reg_lev[ll][jj][ii]) - Cxhalf*(UXRGFLUXIC - UXRGFLUXNL) );
               real_t Vx = HALF*( (V_reg_lev[ll][jj][ii-1] + V_reg_lev[ll][jj][ii]) - Cxhalf*(VXRGFLUXIC - VXRGFLUXNL) );
 
               real_t U_eigen = fabs(Ux/Hx) + sqrt(g*Hx);

               // Cell numbering is the same for the cell to right of the face -- ll l r rr is -2 -1 0 1
               real_t Hll = H_reg_lev[ll][jj][ii-2];
               real_t Hl  = H_reg_lev[ll][jj][ii-1];
               real_t Hr  = H_reg_lev[ll][jj][ii  ];
               real_t Hrr = H_reg_lev[ll][jj][ii+1];
               real_t Ull = U_reg_lev[ll][jj][ii-2];
               real_t Ul  = U_reg_lev[ll][jj][ii-1];
               real_t Ur  = U_reg_lev[ll][jj][ii  ];
               real_t Urr = U_reg_lev[ll][jj][ii+1];

               Wx_H[ll][jj][ii] = w_corrector(deltaT, dx, U_eigen, Hr-Hl, Hl-Hll, Hrr-Hr) * (Hr-Hl);
               Wx_U[ll][jj][ii] = w_corrector(deltaT, dx, U_eigen, Ur-Ul, Ul-Ull, Urr-Ur) * (Ur-Ul);

               HxFlux[ll][jj][ii] = HNEWXRGFLUXFL;
               UxFlux[ll][jj][ii] = UNEWXRGFLUXFL;
               VxFlux[ll][jj][ii] = VNEWXRGFLUXFL;
            }

            if ((mask_reg_lev[ll][jj-1][ii] == 1 || mask_reg_lev[ll][jj][ii] == 1)){
               real_t Hy = HALF*( (H_reg_lev[ll][jj-1][ii] + H_reg_lev[ll][jj][ii]) - Cyhalf*(HYRGFLUXIC - HYRGFLUXNB) );
               real_t Uy = HALF*( (U_reg_lev[ll][jj-1][ii] + U_reg_lev[ll][jj][ii]) - Cyhalf*(UYRGFLUXIC - UYRGFLUXNB) );
               real_t Vy = HALF*( (V_reg_lev[ll][jj-1][ii] + V_reg_lev[ll][jj][ii]) - Cyhalf*(VYRGFLUXIC - VYRGFLUXNB) );

               real_t U_eigen = fabs(Vy/Hy) + sqrt(g*Hy);

               // Cell numbering is the same for the cell to right of the face -- ll l r rr is -2 -1 0 1
               real_t Hbb = H_reg_lev[ll][jj-2][ii];
               real_t Hb  = H_reg_lev[ll][jj-1][ii];
               real_t Ht  = H_reg_lev[ll][jj  ][ii];
               real_t Htt = H_reg_lev[ll][jj+1][ii];
               real_t Vbb = V_reg_lev[ll][jj-2][ii];
               real_t Vb  = V_reg_lev[ll][jj-1][ii];
               real_t Vt  = V_reg_lev[ll][jj  ][ii];
               real_t Vtt = V_reg_lev[ll][jj+1][ii];

               Wy_H[ll][jj][ii] = w_corrector(deltaT, dy, U_eigen, Ht-Hb, Hb-Hbb, Htt-Ht) * (Ht-Hb);
               Wy_V[ll][jj][ii] = w_corrector(deltaT, dy, U_eigen, Vt-Vb, Vb-Vbb, Vtt-Vt) * (Vt-Vb);

               HyFlux[ll][jj][ii] = HNEWYRGFLUXFB;
               UyFlux[ll][jj][ii] = UNEWYRGFLUXFB;
               VyFlux[ll][jj][ii] = VNEWYRGFLUXFB;
            } 
         } // ii
      } // jj
   } // ll

#ifdef _OPENMP
#pragma omp barrier
#pragma omp master
   {
#endif
      cpu_timers[STATE_TIMER_FINITE_DIFFERENCE_PART2] += cpu_timer_stop(tstart_cpu_part);
      cpu_timer_start(&tstart_cpu_part);
#ifdef _OPENMP
   }
#pragma omp barrier
#endif


#ifdef _OPENMP
#ifdef _OPENMP_SIMD
#pragma omp for simd
#else
#pragma omp for
#endif
#else
#ifdef _OPENMP_SIMD
#pragma omp simd
#endif
#endif
   for (int ifixup = 0; ifixup < mesh->nxfixup; ifixup++){
      int ir  = mesh->xrecvIdx[ifixup];
      int ns1 = mesh->xsendIdx1[ifixup];
      int ns2 = mesh->xsendIdx2[ifixup];

      int llr = mesh->xface_level[ir];
      int irjj = mesh->xface_j[ir] - mesh->lev_jregmin[llr];
      int irii = mesh->xface_i[ir] - mesh->lev_iregmin[llr];
      int ll = mesh->xface_level[ns1];
      int ns1jj = mesh->xface_j[ns1] - mesh->lev_jregmin[ll];
      int ns1ii = mesh->xface_i[ns1] - mesh->lev_iregmin[ll];
      int ns2jj = mesh->xface_j[ns2] - mesh->lev_jregmin[ll];
      int ns2ii = mesh->xface_i[ns2] - mesh->lev_iregmin[ll];

      HxFlux[llr][irjj][irii] = (HxFlux[ll][ns1jj][ns1ii] + HxFlux[ll][ns2jj][ns2ii]) * HALF;
      UxFlux[llr][irjj][irii] = (UxFlux[ll][ns1jj][ns1ii] + UxFlux[ll][ns2jj][ns2ii]) * HALF;
      VxFlux[llr][irjj][irii] = (VxFlux[ll][ns1jj][ns1ii] + VxFlux[ll][ns2jj][ns2ii]) * HALF;
      Wx_H[llr][irjj][irii] = (Wx_H[ll][ns1jj][ns1ii] + Wx_H[ll][ns2jj][ns2ii]) * 0.25;
      Wx_U[llr][irjj][irii] = (Wx_U[ll][ns1jj][ns1ii] + Wx_U[ll][ns2jj][ns2ii]) * 0.25;
   }

#ifdef _OPENMP
#ifdef _OPENMP_SIMD
#pragma omp for simd
#else
#pragma omp for
#endif
#else
#ifdef _OPENMP_SIMD
#pragma omp simd
#endif
#endif
   for (int ifixup = 0; ifixup < mesh->nyfixup; ifixup++){
      int ir  = mesh->yrecvIdx[ifixup];
      int ns1 = mesh->ysendIdx1[ifixup];
      int ns2 = mesh->ysendIdx2[ifixup];

      int llr = mesh->yface_level[ir];
      int irjj = mesh->yface_j[ir] - mesh->lev_jregmin[llr];
      int irii = mesh->yface_i[ir] - mesh->lev_iregmin[llr];
      int ll = mesh->yface_level[ns1];
      int ns1jj = mesh->yface_j[ns1] - mesh->lev_jregmin[ll];
      int ns1ii = mesh->yface_i[ns1] - mesh->lev_iregmin[ll];
      int ns2jj = mesh->yface_j[ns2] - mesh->lev_jregmin[ll];
      int ns2ii = mesh->yface_i[ns2] - mesh->lev_iregmin[ll];

      HyFlux[llr][irjj][irii] = (HyFlux[ll][ns1jj][ns1ii] + HyFlux[ll][ns2jj][ns2ii]) * HALF;
      UyFlux[llr][irjj][irii] = (UyFlux[ll][ns1jj][ns1ii] + UyFlux[ll][ns2jj][ns2ii]) * HALF;
      VyFlux[llr][irjj][irii] = (VyFlux[ll][ns1jj][ns1ii] + VyFlux[ll][ns2jj][ns2ii]) * HALF;
      Wy_H[llr][irjj][irii] = (Wy_H[ll][ns1jj][ns1ii] + Wy_H[ll][ns2jj][ns2ii]) * 0.25;
      Wy_V[llr][irjj][irii] = (Wy_V[ll][ns1jj][ns1ii] + Wy_V[ll][ns2jj][ns2ii]) * 0.25;
    }

#ifdef _OPENMP
#pragma omp barrier
#pragma omp master
   {
#endif
      cpu_timers[STATE_TIMER_FINITE_DIFFERENCE_PART3] += cpu_timer_stop(tstart_cpu_part);
      cpu_timer_start(&tstart_cpu_part);
#ifdef _OPENMP
   }
#pragma omp barrier
#endif

   /*
      for (int jj = 2; jj < jjmax-1; jj++) {
         for (int ii = 2; ii < iimax-1; ii++) {
            if (mask_reg[jj][ii-1] == 1 || mask_reg[jj][ii] == 1) {
                int recvIdx = mesh->phantomXFluxRGFace[ll][jj][ii];
                if (recvIdx > -1) {
                    int recvJ = mesh->xface_j[recvIdx] - mesh->lev_jregmin[ll-1];
                    int recvI = mesh->xface_i[recvIdx] - mesh->lev_iregmin[ll-1];

                    if (jj % 2 == 0) { // we are the first to get here
                        HxFlux[ll-1][recvJ][recvI] = 0;
                        UxFlux[ll-1][recvJ][recvI] = 0;
                        VxFlux[ll-1][recvJ][recvI] = 0;
                        Wx_H[ll-1][recvJ][recvI] = 0;
                        Wx_U[ll-1][recvJ][recvI] = 0;
                        passFlagX[ll-1][recvJ][recvI] = 1;
                    }
                    HxFlux[ll-1][recvJ][recvI] += HxFlux[ll][jj][ii] * HALF;
                    UxFlux[ll-1][recvJ][recvI] += UxFlux[ll][jj][ii] * HALF;
                    VxFlux[ll-1][recvJ][recvI] += VxFlux[ll][jj][ii] * HALF;
                    Wx_H[ll-1][recvJ][recvI] += Wx_H[ll][jj][ii] / 4;
                    Wx_U[ll-1][recvJ][recvI] += Wx_U[ll][jj][ii] / 4;
                }
            }

            if (mask_reg[jj-1][ii] == 1 || mask_reg[jj][ii] == 1) {
                int recvIdx = mesh->phantomYFluxRGFace[ll][jj][ii];
                if (recvIdx > -1) {
                    int recvJ = mesh->yface_j[recvIdx] - mesh->lev_jregmin[ll-1];
                    int recvI = mesh->yface_i[recvIdx] - mesh->lev_iregmin[ll-1];
                    if (ii % 2 == 0) { // we are the first to get here
                        HyFlux[ll-1][recvJ][recvI] = 0;
                        UyFlux[ll-1][recvJ][recvI] = 0;
                        VyFlux[ll-1][recvJ][recvI] = 0;
                        Wy_H[ll-1][recvJ][recvI] = 0;
                        Wy_V[ll-1][recvJ][recvI] = 0;
                        passFlagY[ll-1][recvJ][recvI] = 1;
                    }
                    HyFlux[ll-1][recvJ][recvI] += HyFlux[ll][jj][ii] * HALF;
                    UyFlux[ll-1][recvJ][recvI] += UyFlux[ll][jj][ii] * HALF;
                    VyFlux[ll-1][recvJ][recvI] += VyFlux[ll][jj][ii] * HALF;
                    Wy_H[ll-1][recvJ][recvI] += Wy_H[ll][jj][ii] / 4;
                    Wy_V[ll-1][recvJ][recvI] += Wy_V[ll][jj][ii] / 4;
                }
            }
         }
      }*/

#ifdef _OPENMP
#pragma omp barrier
#endif

   for (int ll = 0; ll < mesh->levmx+1; ll++){
      int iimax = mesh->lev_iregsize[ll];
      int jjmax = mesh->lev_jregsize[ll];
      real_t dx = mesh->lev_deltax[ll];
      //real_t dy = mesh->lev_deltay[ll];

#ifdef _OPENMP
#pragma omp for
#endif
      for(int jj=2; jj<jjmax-1; jj++){
#if defined(__GNUC_MINOR__)
         state_t *H_reg_new_loc = states_new[ll][0][jj];
         state_t *U_reg_new_loc = states_new[ll][1][jj];
         state_t *V_reg_new_loc = states_new[ll][2][jj];
         state_t *H_reg_loc = H_reg_lev[ll][jj];
         state_t *U_reg_loc = U_reg_lev[ll][jj];
         state_t *V_reg_loc = V_reg_lev[ll][jj];
         real_t *Hxflux_loc = HxFlux[ll][jj];
         real_t *Hyfluxplus_loc = HyFlux[ll][jj+1];
         real_t *Hyfluxminus_loc = HyFlux[ll][jj];
         real_t *Uxflux_loc = UxFlux[ll][jj];
         real_t *Uyflux_loc = UyFlux[ll][jj];
         real_t *Uyfluxplus_loc = UyFlux[ll][jj+1];
         real_t *Vxflux_loc = VxFlux[ll][jj];
         real_t *Vyflux_loc = VyFlux[ll][jj];
         real_t *Vyfluxplus_loc = VyFlux[ll][jj+1];
         real_t *Wx_H_loc = Wx_H[ll][jj];
         real_t *Wy_H_loc = Wy_H[ll][jj];
         real_t *Wy_Hplus_loc = Wy_H[ll][jj+1];
         real_t *Wx_U_loc = Wx_U[ll][jj];
         real_t *Wy_V_loc = Wy_V[ll][jj];
         real_t *Wy_Vplus_loc = Wy_V[ll][jj+1];
#endif
#ifdef _OPENMP_SIMD
#pragma omp simd
#endif
         for(int ii=2; ii<iimax-1; ii++){
            if (mask_reg_lev[ll][jj][ii] != 1) continue;

            real_t Hxfluxminus = HxFlux[ll][jj][ii];
            real_t Uxfluxminus = UxFlux[ll][jj][ii];
            real_t Vxfluxminus = VxFlux[ll][jj][ii];

            real_t Hxfluxplus = HxFlux[ll][jj][ii+1];
            real_t Uxfluxplus = UxFlux[ll][jj][ii+1];
            real_t Vxfluxplus = VxFlux[ll][jj][ii+1];

            real_t Hyfluxminus = HyFlux[ll][jj][ii];
            real_t Uyfluxminus = UyFlux[ll][jj][ii];
            real_t Vyfluxminus = VyFlux[ll][jj][ii];

            real_t Hyfluxplus = HyFlux[ll][jj+1][ii];
            real_t Uyfluxplus = UyFlux[ll][jj+1][ii];
            real_t Vyfluxplus = VyFlux[ll][jj+1][ii];

            real_t wminusx_H = Wx_H[ll][jj][ii];
            real_t wplusx_H  = Wx_H[ll][jj][ii+1];
            real_t wminusx_U = Wx_U[ll][jj][ii];
            real_t wplusx_U  = Wx_U[ll][jj][ii+1];

            real_t wminusy_H = Wy_H[ll][jj][ii];
            real_t wplusy_H  = Wy_H[ll][jj+1][ii];
            real_t wminusy_V = Wy_V[ll][jj][ii];
            real_t wplusy_V  = Wy_V[ll][jj+1][ii];

            states_new[ll][0][jj][ii] = U_fullstep(deltaT, dx, H_reg_lev[ll][jj][ii],
                       Hxfluxplus, Hxfluxminus, Hyfluxplus, Hyfluxminus)
                  - wminusx_H + wplusx_H - wminusy_H + wplusy_H;

            states_new[ll][1][jj][ii] = U_fullstep(deltaT, dx, U_reg_lev[ll][jj][ii],
                       Uxfluxplus, Uxfluxminus, Uyfluxplus, Uyfluxminus)
                  - wminusx_U + wplusx_U;

            states_new[ll][2][jj][ii] = U_fullstep(deltaT, dx, V_reg_lev[ll][jj][ii],
                       Vxfluxplus, Vxfluxminus, Vyfluxplus, Vyfluxminus)
                  - wminusy_V + wplusy_V;

         } // ii
      } // jj 
   } // ll

#ifdef _OPENMP
#pragma omp barrier
#pragma omp master
   {
#endif
      cpu_timers[STATE_TIMER_FINITE_DIFFERENCE_PART4] += cpu_timer_stop(tstart_cpu_part);
      cpu_timer_start(&tstart_cpu_part);
#ifdef _OPENMP
   }
#pragma omp barrier
#endif

#ifdef _OPENMP
#pragma omp barrier
#pragma omp master
   {
#endif

   for (int ll = 0; ll < mesh->levmx+1; ll++){
      genmatrixfree((void **)HxFlux[ll]);
      genmatrixfree((void **)UxFlux[ll]);
      genmatrixfree((void **)VxFlux[ll]);
      genmatrixfree((void **)Wx_H[ll]);
      genmatrixfree((void **)Wx_U[ll]);
      genmatrixfree((void **)HyFlux[ll]);
      genmatrixfree((void **)UyFlux[ll]);
      genmatrixfree((void **)VyFlux[ll]);
      genmatrixfree((void **)Wy_H[ll]);
      genmatrixfree((void **)Wy_V[ll]);

      // Replace old state with new state
      state_t ***tmp_reg;
      SWAP_PTR(states_new[ll], mesh->meshes[ll].pstate, tmp_reg);
      gentrimatrixfree((void ***) states_new[ll]);
   } // ll


   //free memory
   free(HxFlux);
   free(UxFlux);
   free(VxFlux);
   free(Wx_H);
   free(Wx_U);
   free(HyFlux);
   free(UyFlux);
   free(VyFlux);
   free(Wy_H);
   free(Wy_V);
   free(H_reg_lev);
   free(U_reg_lev);
   free(V_reg_lev);
   free(mask_reg_lev);
   free(states_new);

   mesh->destroy_regular_cell_meshes(state_memory);

      cpu_timers[STATE_TIMER_FINITE_DIFFERENCE_PART5] += cpu_timer_stop(tstart_cpu_part);
      cpu_timer_start(&tstart_cpu_part);

   cpu_timers[STATE_TIMER_FINITE_DIFFERENCE] += cpu_timer_stop(tstart_cpu);
#ifdef _OPENMP
   }
#pragma omp barrier
#endif
}

/************************************************************************/


#ifdef HAVE_OPENCL
void State::gpu_calc_finite_difference(double deltaT)
{
    
   struct timespec tstart_cpu;
   cpu_timer_start(&tstart_cpu);

   cl_command_queue command_queue = ezcl_get_command_queue();

   //cl_mem dev_ptr = NULL;

   size_t &ncells    = mesh->ncells;
   size_t &ncells_ghost = mesh->ncells_ghost;
   if (ncells_ghost < ncells) ncells_ghost = ncells;

   int &levmx           = mesh->levmx;
   cl_mem &dev_celltype = mesh->dev_celltype;
   cl_mem &dev_nlft     = mesh->dev_nlft;
   cl_mem &dev_nrht     = mesh->dev_nrht;
   cl_mem &dev_nbot     = mesh->dev_nbot;
   cl_mem &dev_ntop     = mesh->dev_ntop;
   cl_mem &dev_level    = mesh->dev_level;
   cl_mem &dev_levdx    = mesh->dev_levdx;
   cl_mem &dev_levdy    = mesh->dev_levdy;


   assert(dev_H);
   assert(dev_U);
   assert(dev_V);
   assert(dev_nlft);
   assert(dev_nrht);
   assert(dev_nbot);
   assert(dev_ntop);
   assert(dev_level);
   assert(dev_levdx);
   assert(dev_levdy);

   cl_mem dev_H_new = (cl_mem)gpu_state_memory.memory_malloc(ncells_ghost, sizeof(cl_state_t), const_cast<char *>("dev_H_new"), DEVICE_REGULAR_MEMORY);
   cl_mem dev_U_new = (cl_mem)gpu_state_memory.memory_malloc(ncells_ghost, sizeof(cl_state_t), const_cast<char *>("dev_U_new"), DEVICE_REGULAR_MEMORY);
   cl_mem dev_V_new = (cl_mem)gpu_state_memory.memory_malloc(ncells_ghost, sizeof(cl_state_t), const_cast<char *>("dev_V_new"), DEVICE_REGULAR_MEMORY);
 
   size_t local_work_size = 128;
   size_t global_work_size = ((ncells+local_work_size - 1) /local_work_size) * local_work_size;

#ifdef HAVE_MPI
   if (mesh->numpe > 1) {
      ezcl_set_kernel_arg(kernel_apply_boundary_conditions_local, 0, sizeof(cl_int), &ncells);
      ezcl_set_kernel_arg(kernel_apply_boundary_conditions_local, 1, sizeof(cl_mem), &dev_celltype);
      ezcl_set_kernel_arg(kernel_apply_boundary_conditions_local, 2, sizeof(cl_mem), &dev_nlft);
      ezcl_set_kernel_arg(kernel_apply_boundary_conditions_local, 3, sizeof(cl_mem), &dev_nrht);
      ezcl_set_kernel_arg(kernel_apply_boundary_conditions_local, 4, sizeof(cl_mem), &dev_ntop);
      ezcl_set_kernel_arg(kernel_apply_boundary_conditions_local, 5, sizeof(cl_mem), &dev_nbot);
      ezcl_set_kernel_arg(kernel_apply_boundary_conditions_local, 6, sizeof(cl_mem), &dev_H);
      ezcl_set_kernel_arg(kernel_apply_boundary_conditions_local, 7, sizeof(cl_mem), &dev_U);
      ezcl_set_kernel_arg(kernel_apply_boundary_conditions_local, 8, sizeof(cl_mem), &dev_V);
      ezcl_enqueue_ndrange_kernel(command_queue, kernel_apply_boundary_conditions_local,   1, NULL, &global_work_size, &local_work_size, NULL);
    
        /*
        __kernel void copy_state_data_cl(
                         const int    isize,         // 0
                __global      state_t *H,            // 1
                __global      state_t *U,            // 2
                __global      state_t *V,            // 3
                __global      state_t *H_new,        // 4
                __global      state_t *U_new,        // 5
                __global      state_t *V_new)        // 6
        */

      ezcl_set_kernel_arg(kernel_copy_state_data, 0, sizeof(cl_int), (void *)&ncells);
      ezcl_set_kernel_arg(kernel_copy_state_data, 1, sizeof(cl_mem), (void *)&dev_H);
      ezcl_set_kernel_arg(kernel_copy_state_data, 2, sizeof(cl_mem), (void *)&dev_U);
      ezcl_set_kernel_arg(kernel_copy_state_data, 3, sizeof(cl_mem), (void *)&dev_V);
      ezcl_set_kernel_arg(kernel_copy_state_data, 4, sizeof(cl_mem), (void *)&dev_H_new);
      ezcl_set_kernel_arg(kernel_copy_state_data, 5, sizeof(cl_mem), (void *)&dev_U_new);
      ezcl_set_kernel_arg(kernel_copy_state_data, 6, sizeof(cl_mem), (void *)&dev_V_new);

      //ezcl_enqueue_ndrange_kernel(command_queue, kernel_copy_state_data,   1, NULL, &global_work_size, &local_work_size, &copy_state_data_event);
      ezcl_enqueue_ndrange_kernel(command_queue, kernel_copy_state_data,   1, NULL, &global_work_size, &local_work_size, NULL);

      dev_H = (cl_mem)gpu_state_memory.memory_replace(dev_H, dev_H_new);
      dev_U = (cl_mem)gpu_state_memory.memory_replace(dev_U, dev_U_new);
      dev_V = (cl_mem)gpu_state_memory.memory_replace(dev_V, dev_V_new);

      L7_Dev_Update(dev_H, L7_STATE_T, mesh->cell_handle);
      L7_Dev_Update(dev_U, L7_STATE_T, mesh->cell_handle);
      L7_Dev_Update(dev_V, L7_STATE_T, mesh->cell_handle);

      dev_H_new = (cl_mem)gpu_state_memory.memory_malloc(ncells_ghost, sizeof(cl_state_t), const_cast<char *>("dev_H_new"), DEVICE_REGULAR_MEMORY);
      dev_U_new = (cl_mem)gpu_state_memory.memory_malloc(ncells_ghost, sizeof(cl_state_t), const_cast<char *>("dev_U_new"), DEVICE_REGULAR_MEMORY);
      dev_V_new = (cl_mem)gpu_state_memory.memory_malloc(ncells_ghost, sizeof(cl_state_t), const_cast<char *>("dev_V_new"), DEVICE_REGULAR_MEMORY);

      ezcl_set_kernel_arg(kernel_apply_boundary_conditions_ghost, 0, sizeof(cl_int), &ncells);
      ezcl_set_kernel_arg(kernel_apply_boundary_conditions_ghost, 1, sizeof(cl_mem), &dev_celltype);
      ezcl_set_kernel_arg(kernel_apply_boundary_conditions_ghost, 2, sizeof(cl_mem), &dev_nlft);
      ezcl_set_kernel_arg(kernel_apply_boundary_conditions_ghost, 3, sizeof(cl_mem), &dev_nrht);
      ezcl_set_kernel_arg(kernel_apply_boundary_conditions_ghost, 4, sizeof(cl_mem), &dev_ntop);
      ezcl_set_kernel_arg(kernel_apply_boundary_conditions_ghost, 5, sizeof(cl_mem), &dev_nbot);
      ezcl_set_kernel_arg(kernel_apply_boundary_conditions_ghost, 6, sizeof(cl_mem), &dev_H);
      ezcl_set_kernel_arg(kernel_apply_boundary_conditions_ghost, 7, sizeof(cl_mem), &dev_U);
      ezcl_set_kernel_arg(kernel_apply_boundary_conditions_ghost, 8, sizeof(cl_mem), &dev_V);
      ezcl_enqueue_ndrange_kernel(command_queue, kernel_apply_boundary_conditions_ghost,   1, NULL, &global_work_size, &local_work_size, NULL);
   } else {
      ezcl_set_kernel_arg(kernel_apply_boundary_conditions, 0, sizeof(cl_int), &ncells);
      ezcl_set_kernel_arg(kernel_apply_boundary_conditions, 1, sizeof(cl_mem), &dev_celltype);
      ezcl_set_kernel_arg(kernel_apply_boundary_conditions, 2, sizeof(cl_mem), &dev_nlft);
      ezcl_set_kernel_arg(kernel_apply_boundary_conditions, 3, sizeof(cl_mem), &dev_nrht);
      ezcl_set_kernel_arg(kernel_apply_boundary_conditions, 4, sizeof(cl_mem), &dev_ntop);
      ezcl_set_kernel_arg(kernel_apply_boundary_conditions, 5, sizeof(cl_mem), &dev_nbot);
      ezcl_set_kernel_arg(kernel_apply_boundary_conditions, 6, sizeof(cl_mem), &dev_H);
      ezcl_set_kernel_arg(kernel_apply_boundary_conditions, 7, sizeof(cl_mem), &dev_U);
      ezcl_set_kernel_arg(kernel_apply_boundary_conditions, 8, sizeof(cl_mem), &dev_V);
      ezcl_enqueue_ndrange_kernel(command_queue, kernel_apply_boundary_conditions,   1, NULL, &global_work_size, &local_work_size, NULL);
   }
#else
   ezcl_set_kernel_arg(kernel_apply_boundary_conditions, 0, sizeof(cl_int), &ncells);
   ezcl_set_kernel_arg(kernel_apply_boundary_conditions, 1, sizeof(cl_mem), &dev_celltype);
   ezcl_set_kernel_arg(kernel_apply_boundary_conditions, 2, sizeof(cl_mem), &dev_nlft);
   ezcl_set_kernel_arg(kernel_apply_boundary_conditions, 3, sizeof(cl_mem), &dev_nrht);
   ezcl_set_kernel_arg(kernel_apply_boundary_conditions, 4, sizeof(cl_mem), &dev_ntop);
   ezcl_set_kernel_arg(kernel_apply_boundary_conditions, 5, sizeof(cl_mem), &dev_nbot);
   ezcl_set_kernel_arg(kernel_apply_boundary_conditions, 6, sizeof(cl_mem), &dev_H);
   ezcl_set_kernel_arg(kernel_apply_boundary_conditions, 7, sizeof(cl_mem), &dev_U);
   ezcl_set_kernel_arg(kernel_apply_boundary_conditions, 8, sizeof(cl_mem), &dev_V);
   ezcl_enqueue_ndrange_kernel(command_queue, kernel_apply_boundary_conditions,   1, NULL, &global_work_size, &local_work_size, NULL);
#endif


     /*
     __kernel void calc_finite_difference_cl(
                      const int     ncells,    // 0  Total number of cells.
                      const int     lvmax,     // 1  Maximum level
             __global       state_t *H,        // 2
             __global       state_t *U,        // 3
             __global       state_t *V,        // 4
             __global       state_t *H_new,    // 5
             __global       state_t *U_new,    // 6
             __global       state_t *V_new,    // 7
             __global const int     *nlft,     // 8  Array of left neighbors.
             __global const int     *nrht,     // 9  Array of right neighbors.
             __global const int     *ntop,     // 10  Array of bottom neighbors.
             __global const int     *nbot,     // 11  Array of top neighbors.
             __global const uchar_t *level,    // 12  Array of level information.
                      const real_t   deltaT,   // 13  Size of time step.
             __global const real_t  *lev_dx,   // 14
             __global const real_t  *lev_dy,   // 15
             __local        state4_t *tile,    // 16  Tile size in state4.
             __local        int8  *itile)      // 17  Tile size in int8.
     */
   cl_event calc_finite_difference_event;

   real_t deltaT_local = deltaT;
   ezcl_set_kernel_arg(kernel_calc_finite_difference, 0, sizeof(cl_int),  (void *)&ncells);
   ezcl_set_kernel_arg(kernel_calc_finite_difference, 1, sizeof(cl_int),  (void *)&levmx);
   ezcl_set_kernel_arg(kernel_calc_finite_difference, 2, sizeof(cl_mem),  (void *)&dev_H);
   ezcl_set_kernel_arg(kernel_calc_finite_difference, 3, sizeof(cl_mem),  (void *)&dev_U);
   ezcl_set_kernel_arg(kernel_calc_finite_difference, 4, sizeof(cl_mem),  (void *)&dev_V);
   ezcl_set_kernel_arg(kernel_calc_finite_difference, 5, sizeof(cl_mem),  (void *)&dev_H_new);
   ezcl_set_kernel_arg(kernel_calc_finite_difference, 6, sizeof(cl_mem),  (void *)&dev_U_new);
   ezcl_set_kernel_arg(kernel_calc_finite_difference, 7, sizeof(cl_mem),  (void *)&dev_V_new);
   ezcl_set_kernel_arg(kernel_calc_finite_difference, 8, sizeof(cl_mem),  (void *)&dev_nlft);
   ezcl_set_kernel_arg(kernel_calc_finite_difference, 9, sizeof(cl_mem),  (void *)&dev_nrht);
   ezcl_set_kernel_arg(kernel_calc_finite_difference,10, sizeof(cl_mem),  (void *)&dev_ntop);
   ezcl_set_kernel_arg(kernel_calc_finite_difference,11, sizeof(cl_mem),  (void *)&dev_nbot);
   ezcl_set_kernel_arg(kernel_calc_finite_difference,12, sizeof(cl_mem),  (void *)&dev_level);
   ezcl_set_kernel_arg(kernel_calc_finite_difference,13, sizeof(cl_real_t), (void *)&deltaT_local);
   ezcl_set_kernel_arg(kernel_calc_finite_difference,14, sizeof(cl_mem),  (void *)&dev_levdx);
   ezcl_set_kernel_arg(kernel_calc_finite_difference,15, sizeof(cl_mem),  (void *)&dev_levdy);
   ezcl_set_kernel_arg(kernel_calc_finite_difference,16, local_work_size*sizeof(cl_state4_t),    NULL);
   ezcl_set_kernel_arg(kernel_calc_finite_difference,17, local_work_size*sizeof(cl_int8),    NULL);

   ezcl_enqueue_ndrange_kernel(command_queue, kernel_calc_finite_difference,   1, NULL, &global_work_size, &local_work_size, &calc_finite_difference_event);

   ezcl_wait_for_events(1, &calc_finite_difference_event);
   ezcl_event_release(calc_finite_difference_event);

   dev_H = (cl_mem)gpu_state_memory.memory_replace(dev_H, dev_H_new);
   dev_U = (cl_mem)gpu_state_memory.memory_replace(dev_U, dev_U_new);
   dev_V = (cl_mem)gpu_state_memory.memory_replace(dev_V, dev_V_new);
/*
   vector<state_t>H_loc(ncells);
   ezcl_enqueue_read_buffer(command_queue, dev_H,     CL_TRUE, 0, ncells*sizeof(cl_real_t), &H_loc[0], NULL);
   vector<state_t>U_loc(ncells);
   ezcl_enqueue_read_buffer(command_queue, dev_U,     CL_TRUE, 0, ncells*sizeof(cl_real_t), &U_loc[0], NULL);
   vector<state_t>V_loc(ncells);
   ezcl_enqueue_read_buffer(command_queue, dev_V,     CL_TRUE, 0, ncells*sizeof(cl_real_t), &V_loc[0], NULL);
   printf("\n");
   for (int jello = 0; jello < ncells; jello++) { printf("%d) %f | %f | %f\n", jello, H_loc[jello], U_loc[jello], V_loc[jello]); }
*/
   gpu_timers[STATE_TIMER_FINITE_DIFFERENCE] += (long)(cpu_timer_stop(tstart_cpu)*1.0e9);
}

void State::gpu_faces_setup(size_t mem_requestx, size_t mem_requesty)
{
    dev_HxFlux = ezcl_malloc(NULL, const_cast<char *>("dev_HxFlux"), &mem_requestx, sizeof(cl_real_t), CL_MEM_READ_WRITE, 0);
    dev_UxFlux = ezcl_malloc(NULL, const_cast<char *>("dev_UxFlux"), &mem_requestx, sizeof(cl_real_t), CL_MEM_READ_WRITE, 0);
    dev_VxFlux = ezcl_malloc(NULL, const_cast<char *>("dev_VxFlux"), &mem_requestx, sizeof(cl_real_t), CL_MEM_READ_WRITE, 0);
    dev_HyFlux = ezcl_malloc(NULL, const_cast<char *>("dev_HyFlux"), &mem_requesty, sizeof(cl_real_t), CL_MEM_READ_WRITE, 0);
    dev_UyFlux = ezcl_malloc(NULL, const_cast<char *>("dev_UyFlux"), &mem_requesty, sizeof(cl_real_t), CL_MEM_READ_WRITE, 0);
    dev_VyFlux = ezcl_malloc(NULL, const_cast<char *>("dev_VyFlux"), &mem_requesty, sizeof(cl_real_t), CL_MEM_READ_WRITE, 0);
    dev_Wx_H = ezcl_malloc(NULL, const_cast<char *>("dev_Wx_H"), &mem_requestx, sizeof(cl_real_t), CL_MEM_READ_WRITE, 0);
    dev_Wx_U = ezcl_malloc(NULL, const_cast<char *>("dev_Wx_U"), &mem_requestx, sizeof(cl_real_t), CL_MEM_READ_WRITE, 0);
    dev_Wy_H = ezcl_malloc(NULL, const_cast<char *>("dev_Wy_H"), &mem_requesty, sizeof(cl_real_t), CL_MEM_READ_WRITE, 0);
    dev_Wy_V = ezcl_malloc(NULL, const_cast<char *>("dev_Wy_V"), &mem_requesty, sizeof(cl_real_t), CL_MEM_READ_WRITE, 0);
}

void State::gpu_faces_setup_phantom(size_t mem_request)
{
    dev_Hxfluxplus = ezcl_malloc(NULL, const_cast<char *>("dev_Hxfluxplus"), &mem_request, sizeof(cl_real_t), CL_MEM_READ_WRITE, 0);
    dev_Hxfluxminus = ezcl_malloc(NULL, const_cast<char *>("dev_Hxfluxminus"), &mem_request, sizeof(cl_real_t), CL_MEM_READ_WRITE, 0);
    dev_Uxfluxplus = ezcl_malloc(NULL, const_cast<char *>("dev_Uxfluxplus"), &mem_request, sizeof(cl_real_t), CL_MEM_READ_WRITE, 0);
    dev_Uxfluxminus = ezcl_malloc(NULL, const_cast<char *>("dev_Uxfluxminus"), &mem_request, sizeof(cl_real_t), CL_MEM_READ_WRITE, 0);
    dev_Vxfluxplus = ezcl_malloc(NULL, const_cast<char *>("dev_Vxfluxplus"), &mem_request, sizeof(cl_real_t), CL_MEM_READ_WRITE, 0);
    dev_Vxfluxminus = ezcl_malloc(NULL, const_cast<char *>("dev_Vxfluxminus"), &mem_request, sizeof(cl_real_t), CL_MEM_READ_WRITE, 0);
    dev_Hyfluxplus = ezcl_malloc(NULL, const_cast<char *>("dev_Hyfluxplus"), &mem_request, sizeof(cl_real_t), CL_MEM_READ_WRITE, 0);
    dev_Hyfluxminus = ezcl_malloc(NULL, const_cast<char *>("dev_Hyfluxminus"), &mem_request, sizeof(cl_real_t), CL_MEM_READ_WRITE, 0);
    dev_Uyfluxplus = ezcl_malloc(NULL, const_cast<char *>("dev_Uyfluxplus"), &mem_request, sizeof(cl_real_t), CL_MEM_READ_WRITE, 0);
    dev_Uyfluxminus = ezcl_malloc(NULL, const_cast<char *>("dev_Uyfluxminus"), &mem_request, sizeof(cl_real_t), CL_MEM_READ_WRITE, 0);
    dev_Vyfluxplus = ezcl_malloc(NULL, const_cast<char *>("dev_Vyfluxplus"), &mem_request, sizeof(cl_real_t), CL_MEM_READ_WRITE, 0);
    dev_Vyfluxminus = ezcl_malloc(NULL, const_cast<char *>("dev_Vyfluxminus"), &mem_request, sizeof(cl_real_t), CL_MEM_READ_WRITE, 0);
    dev_Wplusx_H = ezcl_malloc(NULL, const_cast<char *>("dev_Wplusx_H"), &mem_request, sizeof(cl_real_t), CL_MEM_READ_WRITE, 0);
    dev_Wminusx_H = ezcl_malloc(NULL, const_cast<char *>("dev_Wminusx_H"), &mem_request, sizeof(cl_real_t), CL_MEM_READ_WRITE, 0);
    dev_Wplusx_U = ezcl_malloc(NULL, const_cast<char *>("dev_Wplusx_U"), &mem_request, sizeof(cl_real_t), CL_MEM_READ_WRITE, 0);
    dev_Wminusx_U = ezcl_malloc(NULL, const_cast<char *>("dev_Wminusx_U"), &mem_request, sizeof(cl_real_t), CL_MEM_READ_WRITE, 0);
    dev_Wplusy_H = ezcl_malloc(NULL, const_cast<char *>("dev_Wplusy_H"), &mem_request, sizeof(cl_real_t), CL_MEM_READ_WRITE, 0);
    dev_Wminusy_H = ezcl_malloc(NULL, const_cast<char *>("dev_Wminusy_H"), &mem_request, sizeof(cl_real_t), CL_MEM_READ_WRITE, 0);
    dev_Wplusy_V = ezcl_malloc(NULL, const_cast<char *>("dev_Wplusy_V"), &mem_request, sizeof(cl_real_t), CL_MEM_READ_WRITE, 0);
    dev_Wminusy_V = ezcl_malloc(NULL, const_cast<char *>("dev_Wminusy_V"), &mem_request, sizeof(cl_real_t), CL_MEM_READ_WRITE, 0);
}

void State::gpu_faces_delete()
{
   ezcl_device_memory_delete(dev_HxFlux);
   ezcl_device_memory_delete(dev_UxFlux);
   ezcl_device_memory_delete(dev_VxFlux);
   ezcl_device_memory_delete(dev_HyFlux);
   ezcl_device_memory_delete(dev_UyFlux);
   ezcl_device_memory_delete(dev_VyFlux);
   ezcl_device_memory_delete(dev_Wx_H);
   ezcl_device_memory_delete(dev_Wx_U);
   ezcl_device_memory_delete(dev_Wy_H);
   ezcl_device_memory_delete(dev_Wy_V);
}

void State::gpu_faces_delete_phantom()
{
   ezcl_device_memory_delete(dev_Hxfluxplus);
   ezcl_device_memory_delete(dev_Hxfluxminus);
   ezcl_device_memory_delete(dev_Uxfluxplus);
   ezcl_device_memory_delete(dev_Uxfluxminus);
   ezcl_device_memory_delete(dev_Vxfluxplus);
   ezcl_device_memory_delete(dev_Vxfluxminus);
   ezcl_device_memory_delete(dev_Hyfluxplus);
   ezcl_device_memory_delete(dev_Hyfluxminus);
   ezcl_device_memory_delete(dev_Uyfluxplus);
   ezcl_device_memory_delete(dev_Uyfluxminus);
   ezcl_device_memory_delete(dev_Vyfluxplus);
   ezcl_device_memory_delete(dev_Vyfluxminus);
   ezcl_device_memory_delete(dev_Wplusx_H);
   ezcl_device_memory_delete(dev_Wminusx_H);
   ezcl_device_memory_delete(dev_Wplusx_U);
   ezcl_device_memory_delete(dev_Wminusx_U);
   ezcl_device_memory_delete(dev_Wplusy_H);
   ezcl_device_memory_delete(dev_Wminusy_H);
   ezcl_device_memory_delete(dev_Wplusy_V);
   ezcl_device_memory_delete(dev_Wminusy_V);
}

void State::gpu_calc_finite_difference_via_faces(double deltaT)
{

    
   //struct timespec tstart_cpu_part;
   //cpu_timer_start(&tstart_cpu_part);


   cl_command_queue command_queue = ezcl_get_command_queue();

   //cl_mem dev_ptr = NULL;
   mesh->gpu_wbidirmap_setup();

   size_t &ncells    = mesh->ncells;
   size_t &ncells_ghost = mesh->ncells_ghost;
   if (ncells_ghost < ncells) ncells_ghost = ncells;
   int &levmx           = mesh->levmx;
   cl_mem &dev_nface    = mesh->dev_nface;
   cl_mem &dev_celltype = mesh->dev_celltype;
   cl_mem &dev_nlft     = mesh->dev_nlft;
   cl_mem &dev_nrht     = mesh->dev_nrht;
   cl_mem &dev_nbot     = mesh->dev_nbot;
   cl_mem &dev_ntop     = mesh->dev_ntop;
   cl_mem &dev_level    = mesh->dev_level;
   cl_mem &dev_levdx    = mesh->dev_levdx;
   cl_mem &dev_levdy    = mesh->dev_levdy;
   cl_mem &dev_map_xface2cell_lower = mesh->dev_map_xface2cell_lower;
   cl_mem &dev_map_xface2cell_upper = mesh->dev_map_xface2cell_upper;
   cl_mem &dev_map_xcell2face_left1 = mesh->dev_map_xcell2face_left1;
   cl_mem &dev_map_xcell2face_left2 = mesh->dev_map_xcell2face_left2;
   cl_mem &dev_map_xcell2face_right1 = mesh->dev_map_xcell2face_right1;
   cl_mem &dev_map_xcell2face_right2 = mesh->dev_map_xcell2face_right2;
   cl_mem &dev_map_yface2cell_lower = mesh->dev_map_yface2cell_lower;
   cl_mem &dev_map_yface2cell_upper = mesh->dev_map_yface2cell_upper;
   cl_mem &dev_map_ycell2face_bot1 = mesh->dev_map_ycell2face_bot1;
   cl_mem &dev_map_ycell2face_bot2 = mesh->dev_map_ycell2face_bot2;
   cl_mem &dev_map_ycell2face_top1 = mesh->dev_map_ycell2face_top1;
   cl_mem &dev_map_ycell2face_top2 = mesh->dev_map_ycell2face_top2;
    cl_mem dev_xface_level = mesh->dev_xface_level;
    cl_mem dev_xface_i = mesh->dev_xface_i;
    cl_mem dev_xface_j = mesh->dev_xface_j;
    cl_mem dev_ixmin_level = mesh->dev_ixmin_level;
    cl_mem dev_ixmax_level = mesh->dev_ixmax_level;
    cl_mem dev_jxmin_level = mesh->dev_jxmin_level;
    cl_mem dev_jxmax_level = mesh->dev_jxmax_level;
    cl_mem dev_yface_level = mesh->dev_yface_level;
    cl_mem dev_yface_i = mesh->dev_yface_i;
    cl_mem dev_yface_j = mesh->dev_yface_j;
    cl_mem dev_iymin_level = mesh->dev_iymin_level;
    cl_mem dev_iymax_level = mesh->dev_iymax_level;
    cl_mem dev_jymin_level = mesh->dev_jymin_level;
    cl_mem dev_jymax_level = mesh->dev_jymax_level;
   //cl_mem &dev_Hx = mesh->dev_Hx;
   //cl_mem &dev_Ux = mesh->dev_Ux;
   //cl_mem &dev_Vx = mesh->dev_Vx;
   //cl_mem &dev_Hy = mesh->dev_Hy;
   //cl_mem &dev_Uy = mesh->dev_Uy;
   //cl_mem &dev_Vy = mesh->dev_Vy;

   assert(dev_nface);
   assert(dev_H);
   assert(dev_U);
   assert(dev_V);
   assert(dev_nlft);
   assert(dev_nrht);
   assert(dev_nbot);
   assert(dev_ntop);
   assert(dev_level);
   assert(dev_levdx);
   assert(dev_levdy);
   assert(dev_map_xface2cell_lower);
   assert(dev_map_xface2cell_upper);
   assert(dev_map_xcell2face_left1);
   assert(dev_map_xcell2face_left2);
   assert(dev_map_xcell2face_right1);
   assert(dev_map_xcell2face_right2);
   assert(dev_map_yface2cell_lower);
   assert(dev_map_yface2cell_upper);
   assert(dev_map_ycell2face_bot1);
   assert(dev_map_ycell2face_bot2);
   assert(dev_map_ycell2face_top1);
   assert(dev_map_ycell2face_top2);
   assert(dev_xface_level);
   assert(dev_xface_i);
   assert(dev_xface_j);
   assert(dev_ixmin_level);
   assert(dev_ixmax_level);
   assert(dev_jxmin_level);
   assert(dev_jxmax_level);
   assert(dev_yface_level);
   assert(dev_yface_i);
   assert(dev_yface_j);
   assert(dev_iymin_level);
   assert(dev_iymax_level);
   assert(dev_jymin_level);
   assert(dev_jymax_level);
   //assert(dev_Hx);
   //assert(dev_Ux);
   //assert(dev_Vx);
   //assert(dev_Hy);
   //assert(dev_Uy);
   //assert(dev_Vy);

   cl_mem dev_H_new = (cl_mem)gpu_state_memory.memory_malloc(ncells_ghost, sizeof(cl_state_t), const_cast<char *>("dev_H_new"), DEVICE_REGULAR_MEMORY);
   cl_mem dev_U_new = (cl_mem)gpu_state_memory.memory_malloc(ncells_ghost, sizeof(cl_state_t), const_cast<char *>("dev_U_new"), DEVICE_REGULAR_MEMORY);
   cl_mem dev_V_new = (cl_mem)gpu_state_memory.memory_malloc(ncells_ghost, sizeof(cl_state_t), const_cast<char *>("dev_V_new"), DEVICE_REGULAR_MEMORY);
 
   size_t local_work_size = 128;
   size_t global_work_size = ((ncells+local_work_size - 1) /local_work_size) * local_work_size;

#ifdef HAVE_MPI
   if (mesh->numpe > 1) {
      ezcl_set_kernel_arg(kernel_apply_boundary_conditions_local, 0, sizeof(cl_int), &ncells);
      ezcl_set_kernel_arg(kernel_apply_boundary_conditions_local, 1, sizeof(cl_mem), &dev_celltype);
      ezcl_set_kernel_arg(kernel_apply_boundary_conditions_local, 2, sizeof(cl_mem), &dev_nlft);
      ezcl_set_kernel_arg(kernel_apply_boundary_conditions_local, 3, sizeof(cl_mem), &dev_nrht);
      ezcl_set_kernel_arg(kernel_apply_boundary_conditions_local, 4, sizeof(cl_mem), &dev_ntop);
      ezcl_set_kernel_arg(kernel_apply_boundary_conditions_local, 5, sizeof(cl_mem), &dev_nbot);
      ezcl_set_kernel_arg(kernel_apply_boundary_conditions_local, 6, sizeof(cl_mem), &dev_H);
      ezcl_set_kernel_arg(kernel_apply_boundary_conditions_local, 7, sizeof(cl_mem), &dev_U);
      ezcl_set_kernel_arg(kernel_apply_boundary_conditions_local, 8, sizeof(cl_mem), &dev_V);
      ezcl_enqueue_ndrange_kernel(command_queue, kernel_apply_boundary_conditions_local,   1, NULL, &global_work_size, &local_work_size, NULL);
    
        /*
        __kernel void copy_state_data_cl(
                         const int    isize,         // 0
                __global      state_t *H,            // 1
                __global      state_t *U,            // 2
                __global      state_t *V,            // 3
                __global      state_t *H_new,        // 4
                __global      state_t *U_new,        // 5
                __global      state_t *V_new)        // 6
        */

      ezcl_set_kernel_arg(kernel_copy_state_data, 0, sizeof(cl_int), (void *)&ncells);
      ezcl_set_kernel_arg(kernel_copy_state_data, 1, sizeof(cl_mem), (void *)&dev_H);
      ezcl_set_kernel_arg(kernel_copy_state_data, 2, sizeof(cl_mem), (void *)&dev_U);
      ezcl_set_kernel_arg(kernel_copy_state_data, 3, sizeof(cl_mem), (void *)&dev_V);
      ezcl_set_kernel_arg(kernel_copy_state_data, 4, sizeof(cl_mem), (void *)&dev_H_new);
      ezcl_set_kernel_arg(kernel_copy_state_data, 5, sizeof(cl_mem), (void *)&dev_U_new);
      ezcl_set_kernel_arg(kernel_copy_state_data, 6, sizeof(cl_mem), (void *)&dev_V_new);

      //ezcl_enqueue_ndrange_kernel(command_queue, kernel_copy_state_data,   1, NULL, &global_work_size, &local_work_size, &copy_state_data_event);
      ezcl_enqueue_ndrange_kernel(command_queue, kernel_copy_state_data,   1, NULL, &global_work_size, &local_work_size, NULL);

      dev_H = (cl_mem)gpu_state_memory.memory_replace(dev_H, dev_H_new);
      dev_U = (cl_mem)gpu_state_memory.memory_replace(dev_U, dev_U_new);
      dev_V = (cl_mem)gpu_state_memory.memory_replace(dev_V, dev_V_new);

      L7_Dev_Update(dev_H, L7_STATE_T, mesh->cell_handle);
      L7_Dev_Update(dev_U, L7_STATE_T, mesh->cell_handle);
      L7_Dev_Update(dev_V, L7_STATE_T, mesh->cell_handle);

      dev_H_new = (cl_mem)gpu_state_memory.memory_malloc(ncells_ghost, sizeof(cl_state_t), const_cast<char *>("dev_H_new"), DEVICE_REGULAR_MEMORY);
      dev_U_new = (cl_mem)gpu_state_memory.memory_malloc(ncells_ghost, sizeof(cl_state_t), const_cast<char *>("dev_U_new"), DEVICE_REGULAR_MEMORY);
      dev_V_new = (cl_mem)gpu_state_memory.memory_malloc(ncells_ghost, sizeof(cl_state_t), const_cast<char *>("dev_V_new"), DEVICE_REGULAR_MEMORY);

      ezcl_set_kernel_arg(kernel_apply_boundary_conditions_ghost, 0, sizeof(cl_int), &ncells);
      ezcl_set_kernel_arg(kernel_apply_boundary_conditions_ghost, 1, sizeof(cl_mem), &dev_celltype);
      ezcl_set_kernel_arg(kernel_apply_boundary_conditions_ghost, 2, sizeof(cl_mem), &dev_nlft);
      ezcl_set_kernel_arg(kernel_apply_boundary_conditions_ghost, 3, sizeof(cl_mem), &dev_nrht);
      ezcl_set_kernel_arg(kernel_apply_boundary_conditions_ghost, 4, sizeof(cl_mem), &dev_ntop);
      ezcl_set_kernel_arg(kernel_apply_boundary_conditions_ghost, 5, sizeof(cl_mem), &dev_nbot);
      ezcl_set_kernel_arg(kernel_apply_boundary_conditions_ghost, 6, sizeof(cl_mem), &dev_H);
      ezcl_set_kernel_arg(kernel_apply_boundary_conditions_ghost, 7, sizeof(cl_mem), &dev_U);
      ezcl_set_kernel_arg(kernel_apply_boundary_conditions_ghost, 8, sizeof(cl_mem), &dev_V);
      ezcl_enqueue_ndrange_kernel(command_queue, kernel_apply_boundary_conditions_ghost,   1, NULL, &global_work_size, &local_work_size, NULL);
   } else {
      ezcl_set_kernel_arg(kernel_apply_boundary_conditions, 0, sizeof(cl_int), &ncells);
      ezcl_set_kernel_arg(kernel_apply_boundary_conditions, 1, sizeof(cl_mem), &dev_celltype);
      ezcl_set_kernel_arg(kernel_apply_boundary_conditions, 2, sizeof(cl_mem), &dev_nlft);
      ezcl_set_kernel_arg(kernel_apply_boundary_conditions, 3, sizeof(cl_mem), &dev_nrht);
      ezcl_set_kernel_arg(kernel_apply_boundary_conditions, 4, sizeof(cl_mem), &dev_ntop);
      ezcl_set_kernel_arg(kernel_apply_boundary_conditions, 5, sizeof(cl_mem), &dev_nbot);
      ezcl_set_kernel_arg(kernel_apply_boundary_conditions, 6, sizeof(cl_mem), &dev_H);
      ezcl_set_kernel_arg(kernel_apply_boundary_conditions, 7, sizeof(cl_mem), &dev_U);
      ezcl_set_kernel_arg(kernel_apply_boundary_conditions, 8, sizeof(cl_mem), &dev_V);
      ezcl_enqueue_ndrange_kernel(command_queue, kernel_apply_boundary_conditions,   1, NULL, &global_work_size, &local_work_size, NULL);
   }
#else
   ezcl_set_kernel_arg(kernel_apply_boundary_conditions, 0, sizeof(cl_int), &ncells);
   ezcl_set_kernel_arg(kernel_apply_boundary_conditions, 1, sizeof(cl_mem), &dev_celltype);
   ezcl_set_kernel_arg(kernel_apply_boundary_conditions, 2, sizeof(cl_mem), &dev_nlft);
   ezcl_set_kernel_arg(kernel_apply_boundary_conditions, 3, sizeof(cl_mem), &dev_nrht);
   ezcl_set_kernel_arg(kernel_apply_boundary_conditions, 4, sizeof(cl_mem), &dev_ntop);
   ezcl_set_kernel_arg(kernel_apply_boundary_conditions, 5, sizeof(cl_mem), &dev_nbot);
   ezcl_set_kernel_arg(kernel_apply_boundary_conditions, 6, sizeof(cl_mem), &dev_H);
   ezcl_set_kernel_arg(kernel_apply_boundary_conditions, 7, sizeof(cl_mem), &dev_U);
   ezcl_set_kernel_arg(kernel_apply_boundary_conditions, 8, sizeof(cl_mem), &dev_V);
   ezcl_enqueue_ndrange_kernel(command_queue, kernel_apply_boundary_conditions,   1, NULL, &global_work_size, &local_work_size, NULL);
#endif


   //cl_mem dev_Hx = (cl_mem)gpu_state_memory.memory_malloc(nxface, sizeof(cl_state_t), const_cast<char *>("dev_Hx"), DEVICE_REGULAR_MEMORY);
   //cl_mem dev_Ux = (cl_mem)gpu_state_memory.memory_malloc(nxface, sizeof(cl_state_t), const_cast<char *>("dev_Ux"), DEVICE_REGULAR_MEMORY);
   //cl_mem dev_Vx = (cl_mem)gpu_state_memory.memory_malloc(nxface, sizeof(cl_state_t), const_cast<char *>("dev_Vx"), DEVICE_REGULAR_MEMORY);
   //cl_mem dev_Hy = (cl_mem)gpu_state_memory.memory_malloc(nyface, sizeof(cl_state_t), const_cast<char *>("dev_Hy"), DEVICE_REGULAR_MEMORY);
   //cl_mem dev_Uy = (cl_mem)gpu_state_memory.memory_malloc(nyface, sizeof(cl_state_t), const_cast<char *>("dev_Uy"), DEVICE_REGULAR_MEMORY);
   //cl_mem dev_Vy = (cl_mem)gpu_state_memory.memory_malloc(nyface, sizeof(cl_state_t), const_cast<char *>("dev_Vy"), DEVICE_REGULAR_MEMORY);

    //gpu_timers[STATE_TIMER_FINITE_DIFFERENCE_PART1] += (long)(cpu_timer_stop(tstart_cpu_part)*1.0e9);
    //cpu_timer_start(&tstart_cpu_part);

    mesh->gpu_calc_face_list_wbidirmap();

    //gpu_timers[STATE_TIMER_FINITE_DIFFERENCE_PART2] += (long)(cpu_timer_stop(tstart_cpu_part)*1.0e9);
    //cpu_timer_start(&tstart_cpu_part);

    size_t mem_requestx, mem_requesty;
    ezcl_enqueue_read_buffer(command_queue, dev_nface,     CL_TRUE, 0, sizeof(cl_int), &mem_requestx, NULL);
    ezcl_enqueue_read_buffer(command_queue, dev_nface,     CL_TRUE, 1*sizeof(cl_int), sizeof(cl_int), &mem_requesty, NULL);
    //printf("\nMem requests %d and %d\n", mem_requestx, mem_requesty);
    //printf("\n%d\n", ncells);
    
    gpu_faces_setup(mem_requestx, mem_requesty);

   struct timespec tstart_cpu;
   cpu_timer_start(&tstart_cpu);

    /*__kernel void calc_finite_difference_via_faces_face_comps_cl(
            __global          int       *nface,                     // 0 Number array of faces
                        const int       levmx,                      // 1 Maximum level
            __global    const state_t   *H,                         // 2
            __global    const state_t   *U,                         // 3
            __global    const state_t   *V,                         // 4
            __global    const uchar_t   *level,                     // 5 Array of level information
                        const real_t    deltaT,                     // 6 Size of time step
            __global    const real_t    *lev_deltax,                // 7
            __global    const real_t    *lev_deltay,                // 8
            __local           state4_t  *tile,                      // 9Tile size in state4_t
            __local           int8      *itile,                     // 10 Tile size in int8
            __local           int8      *xface,                     // 11 xFace size in int8
            __local           int8      *yface,                     // 12 yFace size in int8 
            __global    const int       *map_xface2cell_lower,      // 13 A face's left cell 
            __global    const int       *map_xface2cell_upper,      // 14 A face's left cell 
            __global    const int       *map_yface2cell_lower,      // 15 A face's below cell 
            __global    const int       *map_yface2cell_upper,      // 16 A face's above cell 
            __global    const int       *map_xcell2face_left1,      // 17 
            __global    const int       *map_xcell2face_left2,      // 18
            __global    const int       *map_xcell2face_right1,     // 19 
            __global    const int       *map_xcell2face_right2,     // 20 
            __global    const int       *map_ycell2face_bot1,       // 21 
            __global    const int       *map_ycell2face_bot2,       // 22 
            __global    const int       *map_ycell2face_top1,       // 23 
            __global    const int       *map_ycell2face_top2,       // 24 
            __global          state_t   *HxFlux,                    // 25
            __global          state_t   *UxFlux,                    // 26
            __global          state_t   *VxFlux,                    // 27
            __global          state_t   *HyFlux,                    // 28
            __global          state_t   *UyFlux,                    // 29
            __global          state_t   *VyFlux,                    // 30
            __global          state_t   *Wx_H,                      // 31
            __global          state_t   *Wx_U,                      // 32
            __global          state_t   *Wy_H,                      // 33
            __global          state_t   *Wy_V,                      // 34
            __global          int       *nlft,                      // 35
            __global          int       *nrht,                      // 36
            __global          int       *nbot,                      // 37
            __global          int       *ntop) {                    // 38
            */

   cl_event calc_finite_difference_via_faces_face_event, calc_finite_difference_via_faces_cell_event;

   //size_t local_face_work = 128;
   //size_t global_face_work = ((MAX(mem_requestx, mem_requesty)+local_face_work - 1) /local_face_work) * local_face_work;
   //printf("\nglobal face work %d\n", global_face_work);

   real_t deltaT_local = deltaT;
   ezcl_set_kernel_arg(kernel_calc_finite_difference_via_faces_face, 0, sizeof(cl_mem), (void *)&dev_nface); 
   ezcl_set_kernel_arg(kernel_calc_finite_difference_via_faces_face, 1, sizeof(cl_int), (void *)&levmx); 
   ezcl_set_kernel_arg(kernel_calc_finite_difference_via_faces_face, 2, sizeof(cl_mem), (void *)&dev_H); 
   ezcl_set_kernel_arg(kernel_calc_finite_difference_via_faces_face, 3, sizeof(cl_mem), (void *)&dev_U); 
   ezcl_set_kernel_arg(kernel_calc_finite_difference_via_faces_face, 4, sizeof(cl_mem), (void *)&dev_V); 
   ezcl_set_kernel_arg(kernel_calc_finite_difference_via_faces_face, 5, sizeof(cl_mem), (void *)&dev_level); 
   ezcl_set_kernel_arg(kernel_calc_finite_difference_via_faces_face, 6, sizeof(cl_real_t), (void *)&deltaT_local); 
   ezcl_set_kernel_arg(kernel_calc_finite_difference_via_faces_face, 7, sizeof(cl_mem), (void *)&dev_levdx); 
   ezcl_set_kernel_arg(kernel_calc_finite_difference_via_faces_face, 8, sizeof(cl_mem), (void *)&dev_levdy); 
   ezcl_set_kernel_arg(kernel_calc_finite_difference_via_faces_face, 9, sizeof(cl_state4_t), NULL); 
   ezcl_set_kernel_arg(kernel_calc_finite_difference_via_faces_face, 10, sizeof(cl_int8), NULL); 
   ezcl_set_kernel_arg(kernel_calc_finite_difference_via_faces_face, 11, sizeof(cl_int8), NULL); 
   ezcl_set_kernel_arg(kernel_calc_finite_difference_via_faces_face, 12, sizeof(cl_int8), NULL); 
   ezcl_set_kernel_arg(kernel_calc_finite_difference_via_faces_face, 13, sizeof(cl_mem), (void *)&dev_map_xface2cell_lower); 
   ezcl_set_kernel_arg(kernel_calc_finite_difference_via_faces_face, 14, sizeof(cl_mem), (void *)&dev_map_xface2cell_upper); 
   ezcl_set_kernel_arg(kernel_calc_finite_difference_via_faces_face, 15, sizeof(cl_mem), (void *)&dev_map_yface2cell_lower); 
   ezcl_set_kernel_arg(kernel_calc_finite_difference_via_faces_face, 16, sizeof(cl_mem), (void *)&dev_map_yface2cell_upper); 
   ezcl_set_kernel_arg(kernel_calc_finite_difference_via_faces_face, 17, sizeof(cl_mem), (void *)&dev_map_xcell2face_left1); 
   ezcl_set_kernel_arg(kernel_calc_finite_difference_via_faces_face, 18, sizeof(cl_mem), (void *)&dev_map_xcell2face_left2); 
   ezcl_set_kernel_arg(kernel_calc_finite_difference_via_faces_face, 19, sizeof(cl_mem), (void *)&dev_map_xcell2face_right1); 
   ezcl_set_kernel_arg(kernel_calc_finite_difference_via_faces_face, 20, sizeof(cl_mem), (void *)&dev_map_xcell2face_right2); 
   ezcl_set_kernel_arg(kernel_calc_finite_difference_via_faces_face, 21, sizeof(cl_mem), (void *)&dev_map_ycell2face_bot1); 
   ezcl_set_kernel_arg(kernel_calc_finite_difference_via_faces_face, 22, sizeof(cl_mem), (void *)&dev_map_ycell2face_bot2); 
   ezcl_set_kernel_arg(kernel_calc_finite_difference_via_faces_face, 23, sizeof(cl_mem), (void *)&dev_map_ycell2face_top1); 
   ezcl_set_kernel_arg(kernel_calc_finite_difference_via_faces_face, 24, sizeof(cl_mem), (void *)&dev_map_ycell2face_top2); 
   ezcl_set_kernel_arg(kernel_calc_finite_difference_via_faces_face, 25, sizeof(cl_mem), (void *)&dev_HxFlux); 
   ezcl_set_kernel_arg(kernel_calc_finite_difference_via_faces_face, 26, sizeof(cl_mem), (void *)&dev_UxFlux); 
   ezcl_set_kernel_arg(kernel_calc_finite_difference_via_faces_face, 27, sizeof(cl_mem), (void *)&dev_VxFlux); 
   ezcl_set_kernel_arg(kernel_calc_finite_difference_via_faces_face, 28, sizeof(cl_mem), (void *)&dev_HyFlux); 
   ezcl_set_kernel_arg(kernel_calc_finite_difference_via_faces_face, 29, sizeof(cl_mem), (void *)&dev_UyFlux); 
   ezcl_set_kernel_arg(kernel_calc_finite_difference_via_faces_face, 30, sizeof(cl_mem), (void *)&dev_VyFlux); 
   ezcl_set_kernel_arg(kernel_calc_finite_difference_via_faces_face, 31, sizeof(cl_mem), (void *)&dev_Wx_H); 
   ezcl_set_kernel_arg(kernel_calc_finite_difference_via_faces_face, 32, sizeof(cl_mem), (void *)&dev_Wx_U); 
   ezcl_set_kernel_arg(kernel_calc_finite_difference_via_faces_face, 33, sizeof(cl_mem), (void *)&dev_Wy_H); 
   ezcl_set_kernel_arg(kernel_calc_finite_difference_via_faces_face, 34, sizeof(cl_mem), (void *)&dev_Wy_V); 
   ezcl_set_kernel_arg(kernel_calc_finite_difference_via_faces_face, 35, sizeof(cl_mem), (void *)&dev_nlft);
   ezcl_set_kernel_arg(kernel_calc_finite_difference_via_faces_face, 36, sizeof(cl_mem), (void *)&dev_nrht);
   ezcl_set_kernel_arg(kernel_calc_finite_difference_via_faces_face, 37, sizeof(cl_mem), (void *)&dev_nbot);
   ezcl_set_kernel_arg(kernel_calc_finite_difference_via_faces_face, 38, sizeof(cl_mem), (void *)&dev_ntop);

   ezcl_enqueue_ndrange_kernel(command_queue, kernel_calc_finite_difference_via_faces_face, 1, NULL, &global_work_size, &local_work_size, &calc_finite_difference_via_faces_face_event);

   ezcl_wait_for_events(1, &calc_finite_difference_via_faces_face_event);
   ezcl_event_release(calc_finite_difference_via_faces_face_event);

   int nxface;
   ezcl_enqueue_read_buffer(command_queue, dev_nface,     CL_TRUE, 0, sizeof(cl_int), &nxface, NULL);
   /* 
   vector<int>lefty(nxface);
   ezcl_enqueue_read_buffer(command_queue, dev_map_xface2cell_lower,     CL_TRUE, 0, nxface*sizeof(cl_int), &lefty[0], NULL);
   vector<int>righty(nxface);
   ezcl_enqueue_read_buffer(command_queue, dev_map_xface2cell_upper,     CL_TRUE, 0, nxface*sizeof(cl_int), &righty[0], NULL);
   vector<state_t>Hx_loc(nxface);
   ezcl_enqueue_read_buffer(command_queue, dev_HxFlux,     CL_TRUE, 0, nxface*sizeof(cl_real_t), &Hx_loc[0], NULL);
   vector<state_t>Ux_loc(nxface);
   ezcl_enqueue_read_buffer(command_queue, dev_UxFlux,     CL_TRUE, 0, nxface*sizeof(cl_real_t), &Ux_loc[0], NULL);
   vector<state_t>Vx_loc(nxface);
   ezcl_enqueue_read_buffer(command_queue, dev_VxFlux,     CL_TRUE, 0, nxface*sizeof(cl_real_t), &Vx_loc[0], NULL);
   printf("\n");
   for (int jello = 0; jello < nxface; jello++) { printf("%d) %d %d\n\t%d) %f | %f | %f\n", jello, lefty[jello], righty[jello], jello, Hx_loc[jello], Ux_loc[jello], Vx_loc[jello]); }
   */

    /*
    __kernel void calc_finite_difference_via_faces_cell_comps_cl (
            __global    const state_t   *H,                         // 0
            __global    const state_t   *U,                         // 1
            __global    const state_t   *V,                         // 2
            __global    const uchar_t   *level,                     // 3 Array of level information
            __global    const int       *map_xface2cell_lower,      // 4 A face's left cell 
            __global    const int       *map_xface2cell_upper,      // 5 A face's left cell 
            __global    const int       *map_yface2cell_lower,      // 6 A face's below cell 
            __global    const int       *map_yface2cell_upper,      // 7 A face's above cell 
            __global    const int       *map_xcell2face_left1,      // 8
            __global    const int       *map_xcell2face_left2,      // 9 
            __global    const int       *map_xcell2face_right1,      // 10 
            __global    const int       *map_xcell2face_right2,      // 11 
            __global    const int       *map_ycell2face_bot1,      // 12 
            __global    const int       *map_ycell2face_bot2,      // 13 
            __global    const int       *map_ycell2face_top1,      // 14 
            __global    const int       *map_ycell2face_top2,      // 15 
            __global          state_t   *HxFlux,                        // 16
            __global          state_t   *UxFlux,                        // 17
            __global          state_t   *VxFlux,                        // 18
            __global          state_t   *HyFlux,                        // 19
            __global          state_t   *UyFlux,                        // 20
            __global          state_t   *VyFlux,                        // 21
            __global          state_t   *Wx_H,                      // 22
            __global          state_t   *Wx_U,                      // 23
            __global          state_t   *Wy_H,                      // 24
            __global          state_t   *Wy_V,                      // 25
            __global          state_t   *H_new,                         // 26
            __global          state_t   *U_new,                         // 27
            __global          state_t   *V_new,                         // 28
                        const int       ncells,                     // 29  Total number of cells
                        const real_t    deltaT,                     // 30 Size of time step
            __global    const real_t    *lev_deltax,                    // 31
            __global    const real_t    *lev_deltay) {                    // 32
     * */ 

   ezcl_set_kernel_arg(kernel_calc_finite_difference_via_faces_cell, 0, sizeof(cl_mem), (void *)&dev_H); 
   ezcl_set_kernel_arg(kernel_calc_finite_difference_via_faces_cell, 1, sizeof(cl_mem), (void *)&dev_U); 
   ezcl_set_kernel_arg(kernel_calc_finite_difference_via_faces_cell, 2, sizeof(cl_mem), (void *)&dev_V); 
   ezcl_set_kernel_arg(kernel_calc_finite_difference_via_faces_cell, 3, sizeof(cl_mem), (void *)&dev_level); 
   ezcl_set_kernel_arg(kernel_calc_finite_difference_via_faces_cell, 4, sizeof(cl_mem), (void *)&dev_map_xface2cell_lower); 
   ezcl_set_kernel_arg(kernel_calc_finite_difference_via_faces_cell, 5, sizeof(cl_mem), (void *)&dev_map_xface2cell_upper); 
   ezcl_set_kernel_arg(kernel_calc_finite_difference_via_faces_cell, 6, sizeof(cl_mem), (void *)&dev_map_yface2cell_lower); 
   ezcl_set_kernel_arg(kernel_calc_finite_difference_via_faces_cell, 7, sizeof(cl_mem), (void *)&dev_map_yface2cell_upper); 
   ezcl_set_kernel_arg(kernel_calc_finite_difference_via_faces_cell, 8, sizeof(cl_mem), (void *)&dev_map_xcell2face_left1); 
   ezcl_set_kernel_arg(kernel_calc_finite_difference_via_faces_cell, 9, sizeof(cl_mem), (void *)&dev_map_xcell2face_left2); 
   ezcl_set_kernel_arg(kernel_calc_finite_difference_via_faces_cell, 10, sizeof(cl_mem), (void *)&dev_map_xcell2face_right1); 
   ezcl_set_kernel_arg(kernel_calc_finite_difference_via_faces_cell, 11, sizeof(cl_mem), (void *)&dev_map_xcell2face_right2); 
   ezcl_set_kernel_arg(kernel_calc_finite_difference_via_faces_cell, 12, sizeof(cl_mem), (void *)&dev_map_ycell2face_bot1); 
   ezcl_set_kernel_arg(kernel_calc_finite_difference_via_faces_cell, 13, sizeof(cl_mem), (void *)&dev_map_ycell2face_bot2); 
   ezcl_set_kernel_arg(kernel_calc_finite_difference_via_faces_cell, 14, sizeof(cl_mem), (void *)&dev_map_ycell2face_top1); 
   ezcl_set_kernel_arg(kernel_calc_finite_difference_via_faces_cell, 15, sizeof(cl_mem), (void *)&dev_map_ycell2face_top2); 
   ezcl_set_kernel_arg(kernel_calc_finite_difference_via_faces_cell, 16, sizeof(cl_mem), (void *)&dev_HxFlux); 
   ezcl_set_kernel_arg(kernel_calc_finite_difference_via_faces_cell, 17, sizeof(cl_mem), (void *)&dev_UxFlux); 
   ezcl_set_kernel_arg(kernel_calc_finite_difference_via_faces_cell, 18, sizeof(cl_mem), (void *)&dev_VxFlux); 
   ezcl_set_kernel_arg(kernel_calc_finite_difference_via_faces_cell, 19, sizeof(cl_mem), (void *)&dev_HyFlux); 
   ezcl_set_kernel_arg(kernel_calc_finite_difference_via_faces_cell, 20, sizeof(cl_mem), (void *)&dev_UyFlux); 
   ezcl_set_kernel_arg(kernel_calc_finite_difference_via_faces_cell, 21, sizeof(cl_mem), (void *)&dev_VyFlux); 
   ezcl_set_kernel_arg(kernel_calc_finite_difference_via_faces_cell, 22, sizeof(cl_mem), (void *)&dev_Wx_H); 
   ezcl_set_kernel_arg(kernel_calc_finite_difference_via_faces_cell, 23, sizeof(cl_mem), (void *)&dev_Wx_U); 
   ezcl_set_kernel_arg(kernel_calc_finite_difference_via_faces_cell, 24, sizeof(cl_mem), (void *)&dev_Wy_H); 
   ezcl_set_kernel_arg(kernel_calc_finite_difference_via_faces_cell, 25, sizeof(cl_mem), (void *)&dev_Wy_V); 
   ezcl_set_kernel_arg(kernel_calc_finite_difference_via_faces_cell, 26, sizeof(cl_mem), (void *)&dev_H_new); 
   ezcl_set_kernel_arg(kernel_calc_finite_difference_via_faces_cell, 27, sizeof(cl_mem), (void *)&dev_U_new); 
   ezcl_set_kernel_arg(kernel_calc_finite_difference_via_faces_cell, 28, sizeof(cl_mem), (void *)&dev_V_new); 
   ezcl_set_kernel_arg(kernel_calc_finite_difference_via_faces_cell, 29, sizeof(cl_int), (void *)&ncells); 
   ezcl_set_kernel_arg(kernel_calc_finite_difference_via_faces_cell, 30, sizeof(cl_real_t), (void *)&deltaT_local); 
   ezcl_set_kernel_arg(kernel_calc_finite_difference_via_faces_cell, 31, sizeof(cl_mem), (void *)&dev_levdx); 
   ezcl_set_kernel_arg(kernel_calc_finite_difference_via_faces_cell, 32, sizeof(cl_mem), (void *)&dev_levdy); 

   ezcl_enqueue_ndrange_kernel(command_queue, kernel_calc_finite_difference_via_faces_cell, 1, NULL, &global_work_size, &local_work_size, &calc_finite_difference_via_faces_cell_event);

   ezcl_wait_for_events(1, &calc_finite_difference_via_faces_cell_event);
   ezcl_event_release(calc_finite_difference_via_faces_cell_event);

   dev_H = (cl_mem)gpu_state_memory.memory_replace(dev_H, dev_H_new);
   dev_U = (cl_mem)gpu_state_memory.memory_replace(dev_U, dev_U_new);
   dev_V = (cl_mem)gpu_state_memory.memory_replace(dev_V, dev_V_new);
   /*  
   vector<state_t>H_loc(ncells);
   ezcl_enqueue_read_buffer(command_queue, dev_H,     CL_TRUE, 0, ncells*sizeof(cl_real_t), &H_loc[0], NULL);
   vector<state_t>U_loc(ncells);
   ezcl_enqueue_read_buffer(command_queue, dev_U,     CL_TRUE, 0, ncells*sizeof(cl_real_t), &U_loc[0], NULL);
   vector<state_t>V_loc(ncells);
   ezcl_enqueue_read_buffer(command_queue, dev_V,     CL_TRUE, 0, ncells*sizeof(cl_real_t), &V_loc[0], NULL);
   
   printf("\n");
   for (int jello = 0; jello < ncells; jello++) { printf("%d) %f | %f | %f\n", jello, H_loc[jello], U_loc[jello], V_loc[jello]); }
   */

   gpu_timers[STATE_TIMER_FINITE_DIFFERENCE] += (long)(cpu_timer_stop(tstart_cpu)*1.0e9);
    
   gpu_faces_delete();
   mesh->gpu_wbidirmap_delete();
}

void State::gpu_calc_finite_difference_in_place(double deltaT)
{

   //struct timespec tstart_cpu_part;
   //cpu_timer_start(&tstart_cpu_part);

    
   cl_command_queue command_queue = ezcl_get_command_queue();


   //cl_mem dev_ptr = NULL;

   size_t &ncells    = mesh->ncells;
   size_t &ncells_ghost = mesh->ncells_ghost;
   if (ncells_ghost < ncells) ncells_ghost = ncells;
   int &levmx           = mesh->levmx;
   real_t deltaT_local = deltaT;
   cl_mem &dev_nface    = mesh->dev_nface;
   cl_mem &dev_celltype = mesh->dev_celltype;
   cl_mem &dev_nlft     = mesh->dev_nlft;
   cl_mem &dev_nrht     = mesh->dev_nrht;
   cl_mem &dev_nbot     = mesh->dev_nbot;
   cl_mem &dev_ntop     = mesh->dev_ntop;
   cl_mem &dev_level    = mesh->dev_level;
   cl_mem &dev_levdx    = mesh->dev_levdx;
   cl_mem &dev_levdy    = mesh->dev_levdy;
   cl_mem &dev_xrecvCIdx = mesh->dev_xrecvCIdx;
   cl_mem &dev_xminusCell2Idx = mesh->dev_xminusCell2Idx;
   cl_mem &dev_xplusCell2Idx = mesh->dev_xplusCell2Idx;
   cl_mem &dev_xsendIdx1 = mesh->dev_xsendIdx1;
   cl_mem &dev_xsendIdx2 = mesh->dev_xsendIdx2;
   cl_mem &dev_yrecvCIdx = mesh->dev_yrecvCIdx;
   cl_mem &dev_yminusCell2Idx = mesh->dev_yminusCell2Idx;
   cl_mem &dev_yplusCell2Idx = mesh->dev_yplusCell2Idx;
   cl_mem &dev_ysendIdx1 = mesh->dev_ysendIdx1;
   cl_mem &dev_ysendIdx2 = mesh->dev_ysendIdx2;
   int &nxfixup = mesh->nxfixup;
   int &nyfixup = mesh->nyfixup;

   assert(dev_nface);
   assert(dev_H);
   assert(dev_U);
   assert(dev_V);
   assert(dev_nlft);
   assert(dev_nrht);
   assert(dev_nbot);
   assert(dev_ntop);
   assert(dev_level);
   assert(dev_levdx);
   assert(dev_levdy);

    mesh->gpu_wbidirmap_setup();
    mesh->gpu_calc_face_list_wbidirmap_phantom(gpu_state_memory, deltaT_local);
    gpu_memory_reset_ptrs();

    //gpu_timers[STATE_TIMER_FINITE_DIFFERENCE_PART1] += (long)(cpu_timer_stop(tstart_cpu_part)*1.0e9);
    //cpu_timer_start(&tstart_cpu_part);

   cl_mem dev_H_new = (cl_mem)gpu_state_memory.memory_malloc(ncells, sizeof(cl_state_t), const_cast<char *>("dev_H_new"), DEVICE_REGULAR_MEMORY);
   cl_mem dev_U_new = (cl_mem)gpu_state_memory.memory_malloc(ncells, sizeof(cl_state_t), const_cast<char *>("dev_U_new"), DEVICE_REGULAR_MEMORY);
   cl_mem dev_V_new = (cl_mem)gpu_state_memory.memory_malloc(ncells, sizeof(cl_state_t), const_cast<char *>("dev_V_new"), DEVICE_REGULAR_MEMORY);
 
    //gpu_timers[STATE_TIMER_WRITE] += (long)(cpu_timer_stop(tstart_cpu_part)*1.0e9);
    //cpu_timer_start(&tstart_cpu_part);

   size_t local_work_size = 128;
   size_t global_work_size = ((ncells+local_work_size - 1) /local_work_size) * local_work_size;
   //size_t global_work_size = MAX(MAX(mesh->nxface, mesh->nyface), ncells);

#ifdef HAVE_MPI
   if (mesh->numpe > 1) {
      ezcl_set_kernel_arg(kernel_apply_boundary_conditions_local, 0, sizeof(cl_int), &ncells);
      ezcl_set_kernel_arg(kernel_apply_boundary_conditions_local, 1, sizeof(cl_mem), &dev_celltype);
      ezcl_set_kernel_arg(kernel_apply_boundary_conditions_local, 2, sizeof(cl_mem), &dev_nlft);
      ezcl_set_kernel_arg(kernel_apply_boundary_conditions_local, 3, sizeof(cl_mem), &dev_nrht);
      ezcl_set_kernel_arg(kernel_apply_boundary_conditions_local, 4, sizeof(cl_mem), &dev_ntop);
      ezcl_set_kernel_arg(kernel_apply_boundary_conditions_local, 5, sizeof(cl_mem), &dev_nbot);
      ezcl_set_kernel_arg(kernel_apply_boundary_conditions_local, 6, sizeof(cl_mem), &dev_H);
      ezcl_set_kernel_arg(kernel_apply_boundary_conditions_local, 7, sizeof(cl_mem), &dev_U);
      ezcl_set_kernel_arg(kernel_apply_boundary_conditions_local, 8, sizeof(cl_mem), &dev_V);
      ezcl_enqueue_ndrange_kernel(command_queue, kernel_apply_boundary_conditions_local,   1, NULL, &global_work_size, &local_work_size, NULL);
    
        /*
        __kernel void copy_state_data_cl(
                         const int    isize,         // 0
                __global      state_t *H,            // 1
                __global      state_t *U,            // 2
                __global      state_t *V,            // 3
                __global      state_t *H_new,        // 4
                __global      state_t *U_new,        // 5
                __global      state_t *V_new)        // 6
        */

      ezcl_set_kernel_arg(kernel_copy_state_data, 0, sizeof(cl_int), (void *)&ncells);
      ezcl_set_kernel_arg(kernel_copy_state_data, 1, sizeof(cl_mem), (void *)&dev_H);
      ezcl_set_kernel_arg(kernel_copy_state_data, 2, sizeof(cl_mem), (void *)&dev_U);
      ezcl_set_kernel_arg(kernel_copy_state_data, 3, sizeof(cl_mem), (void *)&dev_V);
      ezcl_set_kernel_arg(kernel_copy_state_data, 4, sizeof(cl_mem), (void *)&dev_H_new);
      ezcl_set_kernel_arg(kernel_copy_state_data, 5, sizeof(cl_mem), (void *)&dev_U_new);
      ezcl_set_kernel_arg(kernel_copy_state_data, 6, sizeof(cl_mem), (void *)&dev_V_new);

      //ezcl_enqueue_ndrange_kernel(command_queue, kernel_copy_state_data,   1, NULL, &global_work_size, &local_work_size, &copy_state_data_event);
      ezcl_enqueue_ndrange_kernel(command_queue, kernel_copy_state_data,   1, NULL, &global_work_size, &local_work_size, NULL);

      dev_H = (cl_mem)gpu_state_memory.memory_replace(dev_H, dev_H_new);
      dev_U = (cl_mem)gpu_state_memory.memory_replace(dev_U, dev_U_new);
      dev_V = (cl_mem)gpu_state_memory.memory_replace(dev_V, dev_V_new);

      L7_Dev_Update(dev_H, L7_STATE_T, mesh->cell_handle);
      L7_Dev_Update(dev_U, L7_STATE_T, mesh->cell_handle);
      L7_Dev_Update(dev_V, L7_STATE_T, mesh->cell_handle);

      dev_H_new = (cl_mem)gpu_state_memory.memory_malloc(ncells_ghost, sizeof(cl_state_t), const_cast<char *>("dev_H_new"), DEVICE_REGULAR_MEMORY);
      dev_U_new = (cl_mem)gpu_state_memory.memory_malloc(ncells_ghost, sizeof(cl_state_t), const_cast<char *>("dev_U_new"), DEVICE_REGULAR_MEMORY);
      dev_V_new = (cl_mem)gpu_state_memory.memory_malloc(ncells_ghost, sizeof(cl_state_t), const_cast<char *>("dev_V_new"), DEVICE_REGULAR_MEMORY);

      ezcl_set_kernel_arg(kernel_apply_boundary_conditions_ghost, 0, sizeof(cl_int), &ncells);
      ezcl_set_kernel_arg(kernel_apply_boundary_conditions_ghost, 1, sizeof(cl_mem), &dev_celltype);
      ezcl_set_kernel_arg(kernel_apply_boundary_conditions_ghost, 2, sizeof(cl_mem), &dev_nlft);
      ezcl_set_kernel_arg(kernel_apply_boundary_conditions_ghost, 3, sizeof(cl_mem), &dev_nrht);
      ezcl_set_kernel_arg(kernel_apply_boundary_conditions_ghost, 4, sizeof(cl_mem), &dev_ntop);
      ezcl_set_kernel_arg(kernel_apply_boundary_conditions_ghost, 5, sizeof(cl_mem), &dev_nbot);
      ezcl_set_kernel_arg(kernel_apply_boundary_conditions_ghost, 6, sizeof(cl_mem), &dev_H);
      ezcl_set_kernel_arg(kernel_apply_boundary_conditions_ghost, 7, sizeof(cl_mem), &dev_U);
      ezcl_set_kernel_arg(kernel_apply_boundary_conditions_ghost, 8, sizeof(cl_mem), &dev_V);
      ezcl_enqueue_ndrange_kernel(command_queue, kernel_apply_boundary_conditions_ghost,   1, NULL, &global_work_size, &local_work_size, NULL);
   } else {
      ezcl_set_kernel_arg(kernel_apply_boundary_conditions, 0, sizeof(cl_int), &ncells);
      ezcl_set_kernel_arg(kernel_apply_boundary_conditions, 1, sizeof(cl_mem), &dev_celltype);
      ezcl_set_kernel_arg(kernel_apply_boundary_conditions, 2, sizeof(cl_mem), &dev_nlft);
      ezcl_set_kernel_arg(kernel_apply_boundary_conditions, 3, sizeof(cl_mem), &dev_nrht);
      ezcl_set_kernel_arg(kernel_apply_boundary_conditions, 4, sizeof(cl_mem), &dev_ntop);
      ezcl_set_kernel_arg(kernel_apply_boundary_conditions, 5, sizeof(cl_mem), &dev_nbot);
      ezcl_set_kernel_arg(kernel_apply_boundary_conditions, 6, sizeof(cl_mem), &dev_H);
      ezcl_set_kernel_arg(kernel_apply_boundary_conditions, 7, sizeof(cl_mem), &dev_U);
      ezcl_set_kernel_arg(kernel_apply_boundary_conditions, 8, sizeof(cl_mem), &dev_V);
      ezcl_enqueue_ndrange_kernel(command_queue, kernel_apply_boundary_conditions,   1, NULL, &global_work_size, &local_work_size, NULL);
   }
#else
   ezcl_set_kernel_arg(kernel_apply_boundary_conditions, 0, sizeof(cl_int), &ncells);
   ezcl_set_kernel_arg(kernel_apply_boundary_conditions, 1, sizeof(cl_mem), (void *)&dev_celltype);
   ezcl_set_kernel_arg(kernel_apply_boundary_conditions, 2, sizeof(cl_mem), (void *)&dev_nlft);
   ezcl_set_kernel_arg(kernel_apply_boundary_conditions, 3, sizeof(cl_mem), (void *)&dev_nrht);
   ezcl_set_kernel_arg(kernel_apply_boundary_conditions, 4, sizeof(cl_mem), (void *)&dev_ntop);
   ezcl_set_kernel_arg(kernel_apply_boundary_conditions, 5, sizeof(cl_mem), (void *)&dev_nbot);
   ezcl_set_kernel_arg(kernel_apply_boundary_conditions, 6, sizeof(cl_mem), (void *)&dev_H);
   ezcl_set_kernel_arg(kernel_apply_boundary_conditions, 7, sizeof(cl_mem), (void *)&dev_U);
   ezcl_set_kernel_arg(kernel_apply_boundary_conditions, 8, sizeof(cl_mem), (void *)&dev_V);

    //gpu_timers[STATE_TIMER_WRITE] += (long)(cpu_timer_stop(tstart_cpu_part)*1.0e9);
    //cpu_timer_start(&tstart_cpu_part);

   ezcl_enqueue_ndrange_kernel(command_queue, kernel_apply_boundary_conditions,   1, NULL, &global_work_size, &local_work_size, NULL);
#endif

    //gpu_timers[STATE_TIMER_FINITE_DIFFERENCE_PART2] += (long)(cpu_timer_stop(tstart_cpu_part)*1.0e9);
    //cpu_timer_start(&tstart_cpu_part);


   cl_mem &dev_map_xface2cell_lower = mesh->dev_map_xface2cell_lower;
   cl_mem &dev_map_xface2cell_upper = mesh->dev_map_xface2cell_upper;
   cl_mem &dev_map_xcell2face_left1 = mesh->dev_map_xcell2face_left1;
   cl_mem &dev_map_xcell2face_left2 = mesh->dev_map_xcell2face_left2;
   cl_mem &dev_map_xcell2face_right1 = mesh->dev_map_xcell2face_right1;
   cl_mem &dev_map_xcell2face_right2 = mesh->dev_map_xcell2face_right2;
   cl_mem &dev_map_yface2cell_lower = mesh->dev_map_yface2cell_lower;
   cl_mem &dev_map_yface2cell_upper = mesh->dev_map_yface2cell_upper;
   cl_mem &dev_map_ycell2face_bot1 = mesh->dev_map_ycell2face_bot1;
   cl_mem &dev_map_ycell2face_bot2 = mesh->dev_map_ycell2face_bot2;
   cl_mem &dev_map_ycell2face_top1 = mesh->dev_map_ycell2face_top1;
   cl_mem &dev_map_ycell2face_top2 = mesh->dev_map_ycell2face_top2;
    cl_mem dev_xface_level = mesh->dev_xface_level;
    cl_mem dev_xface_i = mesh->dev_xface_i;
    cl_mem dev_xface_j = mesh->dev_xface_j;
    cl_mem dev_ixmin_level = mesh->dev_ixmin_level;
    cl_mem dev_ixmax_level = mesh->dev_ixmax_level;
    cl_mem dev_jxmin_level = mesh->dev_jxmin_level;
    cl_mem dev_jxmax_level = mesh->dev_jxmax_level;
    cl_mem dev_yface_level = mesh->dev_yface_level;
    cl_mem dev_yface_i = mesh->dev_yface_i;
    cl_mem dev_yface_j = mesh->dev_yface_j;
    cl_mem dev_iymin_level = mesh->dev_iymin_level;
    cl_mem dev_iymax_level = mesh->dev_iymax_level;
    cl_mem dev_jymin_level = mesh->dev_jymin_level;
    cl_mem dev_jymax_level = mesh->dev_jymax_level;
   assert(dev_map_xface2cell_lower);
   assert(dev_map_xface2cell_upper);
   assert(dev_map_xcell2face_left1);
   assert(dev_map_xcell2face_left2);
   assert(dev_map_xcell2face_right1);
   assert(dev_map_xcell2face_right2);
   assert(dev_map_yface2cell_lower);
   assert(dev_map_yface2cell_upper);
   assert(dev_map_ycell2face_bot1);
   assert(dev_map_ycell2face_bot2);
   assert(dev_map_ycell2face_top1);
   assert(dev_map_ycell2face_top2);
   assert(dev_xface_level);
   assert(dev_xface_i);
   assert(dev_xface_j);
   assert(dev_ixmin_level);
   assert(dev_ixmax_level);
   assert(dev_jxmin_level);
   assert(dev_jxmax_level);
   assert(dev_yface_level);
   assert(dev_yface_i);
   assert(dev_yface_j);
   assert(dev_iymin_level);
   assert(dev_iymax_level);
   assert(dev_jymin_level);
   assert(dev_jymax_level);

    size_t mem_request;
    //ezcl_enqueue_read_buffer(command_queue, dev_nface,     CL_TRUE, 0, sizeof(cl_int), &mem_requestx, NULL);
    //ezcl_enqueue_read_buffer(command_queue, dev_nface,     CL_TRUE, 1*sizeof(cl_int), sizeof(cl_int), &mem_requesty, NULL);
    //printf("\nMem requests %d and %d\n", mem_requestx, mem_requesty);
    //printf("\n%d\n", ncells);
    //mem_requestx = mesh->pxfaceCnt;
    //mem_requesty = mesh->pyfaceCnt;
    mem_request = ncells;
    
    //gpu_faces_setup_phantom(mem_requestx, mem_requesty);
    gpu_faces_setup_phantom(mem_request);

   cl_event calc_finite_difference_in_place_cell_event, calc_finite_difference_in_place_fill_new_event, calc_finite_difference_in_place_fixup_event;

   //size_t local_face_work = 128;
   //size_t global_face_work = ((pcellCnt+local_face_work - 1) /local_face_work) * local_face_work;
   //printf("\nglobal face work %d\n", global_face_work);

   struct timespec tstart_cpu;
   cpu_timer_start(&tstart_cpu);

    /*
__kernel void calc_finite_difference_in_place_cell_comps_cl (
                        const int       ncells,                     // 0 Number of cells (not including phantom)
            __global    const int       *nfaces,                    // 1 Number of x faces
                        const int       levmx,                      // 2 Maximum level
            __global    const state_t   *H,                         // 3
            __global    const state_t   *U,                         // 4
            __global    const state_t   *V,                         // 5
            __global    const uchar_t   *level,                     // 6 Array of level information
                        const real_t    deltaT,                     // 7 Size of time step
            __global    const real_t    *lev_dx,                    // 8
            __global    const real_t    *lev_dy,                    // 9
            __local           state4_t  *tile,                      // 10 Tile size in state4_t
            __local           int8      *itile,                     // 11 Tile size in int8
            __local           int8      *xface,                     // 12 xFace size in int8
            __local           int8      *yface,                     // 13 yFace size in int8 
            __global    const int       *map_xface2cell_lower,      // 14 A face's left cell 
            __global    const int       *map_xface2cell_upper,      // 15 A face's left cell 
            __global    const int       *map_yface2cell_lower,      // 16 A face's below cell 
            __global    const int       *map_yface2cell_upper,      // 17 A face's above cell 
            __global    const int       *map_xcell2face_left1,      // 18 A cell's left primary face 
            __global    const int       *map_xcell2face_right1,     // 19 A cell's right primary face 
            __global    const int       *map_ycell2face_bot1,       // 20 A cell's bot primary face 
            __global    const int       *map_ycell2face_top1,       // 21 A cell's top primary face 
            __global          state_t   *Hxfluxplus,                // 22
            __global          state_t   *Hxfluxminus,               // 23
            __global          state_t   *Uxfluxplus,                // 24
            __global          state_t   *Uxfluxminus,               // 25
            __global          state_t   *Vxfluxplus,                // 26
            __global          state_t   *Vxfluxminus,               // 27
            __global          state_t   *Hyfluxplus,                // 28
            __global          state_t   *Hyfluxminus,               // 29
            __global          state_t   *Uyfluxplus,                // 30
            __global          state_t   *Uyfluxminus,               // 31
            __global          state_t   *Vyfluxplus,                // 32
            __global          state_t   *Vyfluxminus,               // 33
            __global          state_t   *wplusx_H,                  // 34
            __global          state_t   *wminusx_H,                 // 35
            __global          state_t   *wplusx_U,                  // 36
            __global          state_t   *wminusx_U,                 // 37
            __global          state_t   *wplusy_H,                  // 38
            __global          state_t   *wminusy_H,                 // 39
            __global          state_t   *wplusy_V,                  // 40
            __global          state_t   *wminusy_V) {               // 41
            */

   ezcl_set_kernel_arg(kernel_calc_finite_difference_in_place_cell_comps, 0, sizeof(cl_int), (void *)&ncells); 
   ezcl_set_kernel_arg(kernel_calc_finite_difference_in_place_cell_comps, 1, sizeof(cl_mem), (void *)&dev_nface); 
   ezcl_set_kernel_arg(kernel_calc_finite_difference_in_place_cell_comps, 2, sizeof(cl_int), (void *)&levmx); 
   ezcl_set_kernel_arg(kernel_calc_finite_difference_in_place_cell_comps, 3, sizeof(cl_mem), (void *)&dev_H); 
   ezcl_set_kernel_arg(kernel_calc_finite_difference_in_place_cell_comps, 4, sizeof(cl_mem), (void *)&dev_U); 
   ezcl_set_kernel_arg(kernel_calc_finite_difference_in_place_cell_comps, 5, sizeof(cl_mem), (void *)&dev_V); 
   ezcl_set_kernel_arg(kernel_calc_finite_difference_in_place_cell_comps, 6, sizeof(cl_mem), (void *)&dev_level); 
   ezcl_set_kernel_arg(kernel_calc_finite_difference_in_place_cell_comps, 7, sizeof(cl_real_t), (void *)&deltaT_local); 
   ezcl_set_kernel_arg(kernel_calc_finite_difference_in_place_cell_comps, 8, sizeof(cl_mem), (void *)&dev_levdx); 
   ezcl_set_kernel_arg(kernel_calc_finite_difference_in_place_cell_comps, 9, sizeof(cl_mem), (void *)&dev_levdy); 
   ezcl_set_kernel_arg(kernel_calc_finite_difference_in_place_cell_comps, 10, sizeof(cl_state4_t), NULL); 
   ezcl_set_kernel_arg(kernel_calc_finite_difference_in_place_cell_comps, 11, sizeof(cl_int8), NULL); 
   ezcl_set_kernel_arg(kernel_calc_finite_difference_in_place_cell_comps, 12, sizeof(cl_int8), NULL); 
   ezcl_set_kernel_arg(kernel_calc_finite_difference_in_place_cell_comps, 13, sizeof(cl_int8), NULL); 
   ezcl_set_kernel_arg(kernel_calc_finite_difference_in_place_cell_comps, 14, sizeof(cl_mem), (void *)&dev_map_xface2cell_lower); 
   ezcl_set_kernel_arg(kernel_calc_finite_difference_in_place_cell_comps, 15, sizeof(cl_mem), (void *)&dev_map_xface2cell_upper); 
   ezcl_set_kernel_arg(kernel_calc_finite_difference_in_place_cell_comps, 16, sizeof(cl_mem), (void *)&dev_map_yface2cell_lower); 
   ezcl_set_kernel_arg(kernel_calc_finite_difference_in_place_cell_comps, 17, sizeof(cl_mem), (void *)&dev_map_yface2cell_upper); 
   ezcl_set_kernel_arg(kernel_calc_finite_difference_in_place_cell_comps, 18, sizeof(cl_mem), (void *)&dev_map_xcell2face_left1); 
   ezcl_set_kernel_arg(kernel_calc_finite_difference_in_place_cell_comps, 19, sizeof(cl_mem), (void *)&dev_map_xcell2face_right1); 
   ezcl_set_kernel_arg(kernel_calc_finite_difference_in_place_cell_comps, 20, sizeof(cl_mem), (void *)&dev_map_ycell2face_bot1); 
   ezcl_set_kernel_arg(kernel_calc_finite_difference_in_place_cell_comps, 21, sizeof(cl_mem), (void *)&dev_map_ycell2face_top1); 
   ezcl_set_kernel_arg(kernel_calc_finite_difference_in_place_cell_comps, 22, sizeof(cl_mem), (void *)&dev_Hxfluxplus); 
   ezcl_set_kernel_arg(kernel_calc_finite_difference_in_place_cell_comps, 23, sizeof(cl_mem), (void *)&dev_Hxfluxminus); 
   ezcl_set_kernel_arg(kernel_calc_finite_difference_in_place_cell_comps, 24, sizeof(cl_mem), (void *)&dev_Uxfluxplus); 
   ezcl_set_kernel_arg(kernel_calc_finite_difference_in_place_cell_comps, 25, sizeof(cl_mem), (void *)&dev_Uxfluxminus); 
   ezcl_set_kernel_arg(kernel_calc_finite_difference_in_place_cell_comps, 26, sizeof(cl_mem), (void *)&dev_Vxfluxplus); 
   ezcl_set_kernel_arg(kernel_calc_finite_difference_in_place_cell_comps, 27, sizeof(cl_mem), (void *)&dev_Vxfluxminus); 
   ezcl_set_kernel_arg(kernel_calc_finite_difference_in_place_cell_comps, 28, sizeof(cl_mem), (void *)&dev_Hyfluxplus); 
   ezcl_set_kernel_arg(kernel_calc_finite_difference_in_place_cell_comps, 29, sizeof(cl_mem), (void *)&dev_Hyfluxminus); 
   ezcl_set_kernel_arg(kernel_calc_finite_difference_in_place_cell_comps, 30, sizeof(cl_mem), (void *)&dev_Uyfluxplus); 
   ezcl_set_kernel_arg(kernel_calc_finite_difference_in_place_cell_comps, 31, sizeof(cl_mem), (void *)&dev_Uyfluxminus); 
   ezcl_set_kernel_arg(kernel_calc_finite_difference_in_place_cell_comps, 32, sizeof(cl_mem), (void *)&dev_Vyfluxplus); 
   ezcl_set_kernel_arg(kernel_calc_finite_difference_in_place_cell_comps, 33, sizeof(cl_mem), (void *)&dev_Vyfluxminus); 
   ezcl_set_kernel_arg(kernel_calc_finite_difference_in_place_cell_comps, 34, sizeof(cl_mem), (void *)&dev_Wplusx_H); 
   ezcl_set_kernel_arg(kernel_calc_finite_difference_in_place_cell_comps, 35, sizeof(cl_mem), (void *)&dev_Wminusx_H); 
   ezcl_set_kernel_arg(kernel_calc_finite_difference_in_place_cell_comps, 36, sizeof(cl_mem), (void *)&dev_Wplusx_U); 
   ezcl_set_kernel_arg(kernel_calc_finite_difference_in_place_cell_comps, 37, sizeof(cl_mem), (void *)&dev_Wminusx_U); 
   ezcl_set_kernel_arg(kernel_calc_finite_difference_in_place_cell_comps, 38, sizeof(cl_mem), (void *)&dev_Wplusy_H); 
   ezcl_set_kernel_arg(kernel_calc_finite_difference_in_place_cell_comps, 39, sizeof(cl_mem), (void *)&dev_Wminusy_H); 
   ezcl_set_kernel_arg(kernel_calc_finite_difference_in_place_cell_comps, 40, sizeof(cl_mem), (void *)&dev_Wplusy_V); 
   ezcl_set_kernel_arg(kernel_calc_finite_difference_in_place_cell_comps, 41, sizeof(cl_mem), (void *)&dev_Wminusy_V); 

    //gpu_timers[STATE_TIMER_WRITE] += (long)(cpu_timer_stop(tstart_cpu_part)*1.0e9);
    //cpu_timer_start(&tstart_cpu_part);

   ezcl_enqueue_ndrange_kernel(command_queue, kernel_calc_finite_difference_in_place_cell_comps, 1, NULL, &global_work_size, &local_work_size, &calc_finite_difference_in_place_cell_event);

   ezcl_wait_for_events(1, &calc_finite_difference_in_place_cell_event);
   ezcl_event_release(calc_finite_difference_in_place_cell_event);
   
    //gpu_timers[STATE_TIMER_FINITE_DIFFERENCE_PART3] += (long)(cpu_timer_stop(tstart_cpu_part)*1.0e9);
    //cpu_timer_start(&tstart_cpu_part);

   /*
__kernel void calc_finite_difference_in_place_fixup_cl(
                        const int        nxfixup,                   // 0
                        const int        nyfixup,                   // 1
            __global    const int       *xrecvCIdx,                 // 2
            __global    const int       *xplusCell2Idx,             // 3
            __global    const int       *xminusCell2Idx,            // 4
            __global    const int       *xsendIdx1,                 // 5
            __global    const int       *xsendIdx2,                 // 6
            __global    const int       *yrecvCIdx,                 // 7
            __global    const int       *yplusCell2Idx,             // 8
            __global    const int       *yminusCell2Idx,            // 9
            __global    const int       *ysendIdx1,                 // 10
            __global    const int       *ysendIdx2,                 // 11
            __global    const int       *map_xface2cell_lower,      // 12
            __global    const int       *map_xface2cell_upper,      // 13
            __global    const int       *map_yface2cell_lower,      // 14 
            __global    const int       *map_yface2cell_upper,      // 15 
            __global          state_t   *Hxfluxplus,                // 16
            __global          state_t   *Hxfluxminus,               // 17
            __global          state_t   *Uxfluxplus,                // 18
            __global          state_t   *Uxfluxminus,               // 19
            __global          state_t   *Vxfluxplus,                // 20
            __global          state_t   *Vxfluxminus,               // 21
            __global          state_t   *Hyfluxplus,                // 22
            __global          state_t   *Hyfluxminus,               // 23
            __global          state_t   *Uyfluxplus,                // 24
            __global          state_t   *Uyfluxminus,               // 25
            __global          state_t   *Vyfluxplus,                // 26
            __global          state_t   *Vyfluxminus,               // 27
            __global          state_t   *wplusx_H,                  // 28
            __global          state_t   *wminusx_H,                 // 29
            __global          state_t   *wplusx_U,                  // 30
            __global          state_t   *wminusx_U,                 // 31
            __global          state_t   *wplusy_H,                  // 32
            __global          state_t   *wminusy_H,                 // 33
            __global          state_t   *wplusy_V,                  // 34
            __global          state_t   *wminusy_V) {               // 35
    */

    ezcl_set_kernel_arg(kernel_calc_finite_difference_in_place_fixup, 0, sizeof(cl_int), (void *)&nxfixup); 
    ezcl_set_kernel_arg(kernel_calc_finite_difference_in_place_fixup, 1, sizeof(cl_int), (void *)&nyfixup); 
    ezcl_set_kernel_arg(kernel_calc_finite_difference_in_place_fixup, 2, sizeof(cl_mem), (void *)&dev_xrecvCIdx);
    ezcl_set_kernel_arg(kernel_calc_finite_difference_in_place_fixup, 3, sizeof(cl_mem), (void *)&dev_xplusCell2Idx);
    ezcl_set_kernel_arg(kernel_calc_finite_difference_in_place_fixup, 4, sizeof(cl_mem), (void *)&dev_xminusCell2Idx);
    ezcl_set_kernel_arg(kernel_calc_finite_difference_in_place_fixup, 5, sizeof(cl_mem), (void *)&dev_xsendIdx1);
    ezcl_set_kernel_arg(kernel_calc_finite_difference_in_place_fixup, 6, sizeof(cl_mem), (void *)&dev_xsendIdx2);
    ezcl_set_kernel_arg(kernel_calc_finite_difference_in_place_fixup, 7, sizeof(cl_mem), (void *)&dev_yrecvCIdx);
    ezcl_set_kernel_arg(kernel_calc_finite_difference_in_place_fixup, 8, sizeof(cl_mem), (void *)&dev_yplusCell2Idx);
    ezcl_set_kernel_arg(kernel_calc_finite_difference_in_place_fixup, 9, sizeof(cl_mem), (void *)&dev_yminusCell2Idx);
    ezcl_set_kernel_arg(kernel_calc_finite_difference_in_place_fixup, 10, sizeof(cl_mem), (void *)&dev_ysendIdx1);
    ezcl_set_kernel_arg(kernel_calc_finite_difference_in_place_fixup, 11, sizeof(cl_mem), (void *)&dev_ysendIdx2);
    ezcl_set_kernel_arg(kernel_calc_finite_difference_in_place_fixup, 12, sizeof(cl_mem), (void *)&dev_map_xface2cell_lower); 
    ezcl_set_kernel_arg(kernel_calc_finite_difference_in_place_fixup, 13, sizeof(cl_mem), (void *)&dev_map_xface2cell_upper); 
    ezcl_set_kernel_arg(kernel_calc_finite_difference_in_place_fixup, 14, sizeof(cl_mem), (void *)&dev_map_yface2cell_lower); 
    ezcl_set_kernel_arg(kernel_calc_finite_difference_in_place_fixup, 15, sizeof(cl_mem), (void *)&dev_map_yface2cell_upper); 
    ezcl_set_kernel_arg(kernel_calc_finite_difference_in_place_fixup, 16, sizeof(cl_mem), (void *)&dev_Hxfluxplus); 
    ezcl_set_kernel_arg(kernel_calc_finite_difference_in_place_fixup, 17, sizeof(cl_mem), (void *)&dev_Hxfluxminus); 
    ezcl_set_kernel_arg(kernel_calc_finite_difference_in_place_fixup, 18, sizeof(cl_mem), (void *)&dev_Uxfluxplus); 
    ezcl_set_kernel_arg(kernel_calc_finite_difference_in_place_fixup, 19, sizeof(cl_mem), (void *)&dev_Uxfluxminus); 
    ezcl_set_kernel_arg(kernel_calc_finite_difference_in_place_fixup, 20, sizeof(cl_mem), (void *)&dev_Vxfluxplus); 
    ezcl_set_kernel_arg(kernel_calc_finite_difference_in_place_fixup, 21, sizeof(cl_mem), (void *)&dev_Vxfluxminus); 
    ezcl_set_kernel_arg(kernel_calc_finite_difference_in_place_fixup, 22, sizeof(cl_mem), (void *)&dev_Hyfluxplus); 
    ezcl_set_kernel_arg(kernel_calc_finite_difference_in_place_fixup, 23, sizeof(cl_mem), (void *)&dev_Hyfluxminus); 
    ezcl_set_kernel_arg(kernel_calc_finite_difference_in_place_fixup, 24, sizeof(cl_mem), (void *)&dev_Uyfluxplus); 
    ezcl_set_kernel_arg(kernel_calc_finite_difference_in_place_fixup, 25, sizeof(cl_mem), (void *)&dev_Uyfluxminus); 
    ezcl_set_kernel_arg(kernel_calc_finite_difference_in_place_fixup, 26, sizeof(cl_mem), (void *)&dev_Vyfluxplus); 
    ezcl_set_kernel_arg(kernel_calc_finite_difference_in_place_fixup, 27, sizeof(cl_mem), (void *)&dev_Vyfluxminus); 
    ezcl_set_kernel_arg(kernel_calc_finite_difference_in_place_fixup, 28, sizeof(cl_mem), (void *)&dev_Wplusx_H); 
    ezcl_set_kernel_arg(kernel_calc_finite_difference_in_place_fixup, 29, sizeof(cl_mem), (void *)&dev_Wminusx_H); 
    ezcl_set_kernel_arg(kernel_calc_finite_difference_in_place_fixup, 30, sizeof(cl_mem), (void *)&dev_Wplusx_U); 
    ezcl_set_kernel_arg(kernel_calc_finite_difference_in_place_fixup, 31, sizeof(cl_mem), (void *)&dev_Wminusx_U); 
    ezcl_set_kernel_arg(kernel_calc_finite_difference_in_place_fixup, 32, sizeof(cl_mem), (void *)&dev_Wplusy_H); 
    ezcl_set_kernel_arg(kernel_calc_finite_difference_in_place_fixup, 33, sizeof(cl_mem), (void *)&dev_Wminusy_H); 
    ezcl_set_kernel_arg(kernel_calc_finite_difference_in_place_fixup, 34, sizeof(cl_mem), (void *)&dev_Wplusy_V); 
    ezcl_set_kernel_arg(kernel_calc_finite_difference_in_place_fixup, 35, sizeof(cl_mem), (void *)&dev_Wminusy_V); 

    //gpu_timers[STATE_TIMER_WRITE] += (long)(cpu_timer_stop(tstart_cpu_part)*1.0e9);
    //cpu_timer_start(&tstart_cpu_part);

    ezcl_enqueue_ndrange_kernel(command_queue, kernel_calc_finite_difference_in_place_fixup, 1, NULL, &global_work_size, &local_work_size, &calc_finite_difference_in_place_fixup_event);

    ezcl_wait_for_events(1, &calc_finite_difference_in_place_fixup_event);
    ezcl_event_release(calc_finite_difference_in_place_fixup_event);

    //gpu_timers[STATE_TIMER_FINITE_DIFFERENCE_PART4] += (long)(cpu_timer_stop(tstart_cpu_part)*1.0e9);
    //cpu_timer_start(&tstart_cpu_part);

  /*  
   vector<real_t>H_loc(ncells, 0);
   ezcl_enqueue_read_buffer(command_queue, dev_Wplusx_H,     CL_TRUE, 0, ncells*sizeof(cl_real_t), &H_loc[0], NULL);
   vector<real_t>U_loc(ncells, 0);
   ezcl_enqueue_read_buffer(command_queue, dev_Wminusx_H,     CL_TRUE, 0, ncells*sizeof(cl_real_t), &U_loc[0], NULL);
   vector<real_t>V_loc(ncells, 0);
   ezcl_enqueue_read_buffer(command_queue, dev_Wplusy_H,     CL_TRUE, 0, ncells*sizeof(cl_real_t), &V_loc[0], NULL);
   vector<real_t>Z_loc(ncells, 0);
   ezcl_enqueue_read_buffer(command_queue, dev_Wminusy_H,     CL_TRUE, 0, ncells*sizeof(cl_real_t), &Z_loc[0], NULL);
   printf("\n");
   for (int jello = 0; jello < ncells; jello++) { printf("%d) %f | %f | %f | %f\n", jello, H_loc[jello], U_loc[jello], V_loc[jello], Z_loc[jello]); }
*/


/*
__kernel void calc_finite_difference_in_place_fill_new_cl(
                        const int       ncells,                     // 0 Number of cells (not including phantom)
                        const real_t    deltaT,                     // 1 Size of time step
            __global    const real_t    *lev_dx,                    // 2
            __global    const real_t    *lev_dy,                    // 3
            __global          state_t   *Hxfluxplus,                // 4
            __global          state_t   *Hxfluxminus,               // 5
            __global          state_t   *Uxfluxplus,                // 6
            __global          state_t   *Uxfluxminus,               // 7
            __global          state_t   *Vxfluxplus,                // 8
            __global          state_t   *Vxfluxminus,               // 9
            __global          state_t   *Hyfluxplus,                // 10
            __global          state_t   *Hyfluxminus,               // 11
            __global          state_t   *Uyfluxplus,                // 12
            __global          state_t   *Uyfluxminus,               // 13
            __global          state_t   *Vyfluxplus,                // 14
            __global          state_t   *Vyfluxminus,               // 15
            __global          state_t   *wplusx_H,                  // 16
            __global          state_t   *wminusx_H,                 // 17
            __global          state_t   *wplusx_U,                  // 18
            __global          state_t   *wminusx_U,                 // 19
            __global          state_t   *wplusy_H,                  // 20
            __global          state_t   *wminusy_H,                 // 21
            __global          state_t   *wplusy_V,                  // 22
            __global          state_t   *wminusy_V,                 // 23
            __global          state_t   *level,                     // 24
            __global          state_t   *H,                         // 25
            __global          state_t   *U,                         // 26
            __global          state_t   *V,                         // 27
            __global          state_t   *H_new,                     // 28
            __global          state_t   *U_new,                     // 29
            __global          state_t   *V_new) {                   // 30
*/
  

    ezcl_set_kernel_arg(kernel_calc_finite_difference_in_place_fill_new, 0, sizeof(cl_int), (void *)&ncells); 
    ezcl_set_kernel_arg(kernel_calc_finite_difference_in_place_fill_new, 1, sizeof(cl_real_t), (void *)&deltaT_local); 
    ezcl_set_kernel_arg(kernel_calc_finite_difference_in_place_fill_new, 2, sizeof(cl_mem), (void *)&dev_levdx); 
    ezcl_set_kernel_arg(kernel_calc_finite_difference_in_place_fill_new, 3, sizeof(cl_mem), (void *)&dev_levdy); 
    ezcl_set_kernel_arg(kernel_calc_finite_difference_in_place_fill_new, 4, sizeof(cl_mem), (void *)&dev_Hxfluxplus); 
    ezcl_set_kernel_arg(kernel_calc_finite_difference_in_place_fill_new, 5, sizeof(cl_mem), (void *)&dev_Hxfluxminus); 
    ezcl_set_kernel_arg(kernel_calc_finite_difference_in_place_fill_new, 6, sizeof(cl_mem), (void *)&dev_Uxfluxplus); 
    ezcl_set_kernel_arg(kernel_calc_finite_difference_in_place_fill_new, 7, sizeof(cl_mem), (void *)&dev_Uxfluxminus); 
    ezcl_set_kernel_arg(kernel_calc_finite_difference_in_place_fill_new, 8, sizeof(cl_mem), (void *)&dev_Vxfluxplus); 
    ezcl_set_kernel_arg(kernel_calc_finite_difference_in_place_fill_new, 9, sizeof(cl_mem), (void *)&dev_Vxfluxminus); 
    ezcl_set_kernel_arg(kernel_calc_finite_difference_in_place_fill_new, 10, sizeof(cl_mem), (void *)&dev_Hyfluxplus); 
    ezcl_set_kernel_arg(kernel_calc_finite_difference_in_place_fill_new, 11, sizeof(cl_mem), (void *)&dev_Hyfluxminus); 
    ezcl_set_kernel_arg(kernel_calc_finite_difference_in_place_fill_new, 12, sizeof(cl_mem), (void *)&dev_Uyfluxplus); 
    ezcl_set_kernel_arg(kernel_calc_finite_difference_in_place_fill_new, 13, sizeof(cl_mem), (void *)&dev_Uyfluxminus); 
    ezcl_set_kernel_arg(kernel_calc_finite_difference_in_place_fill_new, 14, sizeof(cl_mem), (void *)&dev_Vyfluxplus); 
    ezcl_set_kernel_arg(kernel_calc_finite_difference_in_place_fill_new, 15, sizeof(cl_mem), (void *)&dev_Vyfluxminus); 
    ezcl_set_kernel_arg(kernel_calc_finite_difference_in_place_fill_new, 16, sizeof(cl_mem), (void *)&dev_Wplusx_H); 
    ezcl_set_kernel_arg(kernel_calc_finite_difference_in_place_fill_new, 17, sizeof(cl_mem), (void *)&dev_Wminusx_H); 
    ezcl_set_kernel_arg(kernel_calc_finite_difference_in_place_fill_new, 18, sizeof(cl_mem), (void *)&dev_Wplusx_U); 
    ezcl_set_kernel_arg(kernel_calc_finite_difference_in_place_fill_new, 19, sizeof(cl_mem), (void *)&dev_Wminusx_U); 
    ezcl_set_kernel_arg(kernel_calc_finite_difference_in_place_fill_new, 20, sizeof(cl_mem), (void *)&dev_Wplusy_H); 
    ezcl_set_kernel_arg(kernel_calc_finite_difference_in_place_fill_new, 21, sizeof(cl_mem), (void *)&dev_Wminusy_H); 
    ezcl_set_kernel_arg(kernel_calc_finite_difference_in_place_fill_new, 22, sizeof(cl_mem), (void *)&dev_Wplusy_V); 
    ezcl_set_kernel_arg(kernel_calc_finite_difference_in_place_fill_new, 23, sizeof(cl_mem), (void *)&dev_Wminusy_V); 
    ezcl_set_kernel_arg(kernel_calc_finite_difference_in_place_fill_new, 24, sizeof(cl_mem), (void *)&dev_level); 
    ezcl_set_kernel_arg(kernel_calc_finite_difference_in_place_fill_new, 25, sizeof(cl_mem), (void *)&dev_map_xface2cell_lower); 
    ezcl_set_kernel_arg(kernel_calc_finite_difference_in_place_fill_new, 26, sizeof(cl_mem), (void *)&dev_map_xface2cell_upper); 
    ezcl_set_kernel_arg(kernel_calc_finite_difference_in_place_fill_new, 27, sizeof(cl_mem), (void *)&dev_map_yface2cell_lower); 
    ezcl_set_kernel_arg(kernel_calc_finite_difference_in_place_fill_new, 28, sizeof(cl_mem), (void *)&dev_map_yface2cell_upper); 
    ezcl_set_kernel_arg(kernel_calc_finite_difference_in_place_fill_new, 29, sizeof(cl_mem), (void *)&dev_map_xcell2face_left1); 
    ezcl_set_kernel_arg(kernel_calc_finite_difference_in_place_fill_new, 30, sizeof(cl_mem), (void *)&dev_map_xcell2face_right1); 
    ezcl_set_kernel_arg(kernel_calc_finite_difference_in_place_fill_new, 31, sizeof(cl_mem), (void *)&dev_map_ycell2face_bot1); 
    ezcl_set_kernel_arg(kernel_calc_finite_difference_in_place_fill_new, 32, sizeof(cl_mem), (void *)&dev_map_ycell2face_top1); 
    ezcl_set_kernel_arg(kernel_calc_finite_difference_in_place_fill_new, 33, sizeof(cl_mem), (void *)&dev_H); 
    ezcl_set_kernel_arg(kernel_calc_finite_difference_in_place_fill_new, 34, sizeof(cl_mem), (void *)&dev_U); 
    ezcl_set_kernel_arg(kernel_calc_finite_difference_in_place_fill_new, 35, sizeof(cl_mem), (void *)&dev_V); 
    ezcl_set_kernel_arg(kernel_calc_finite_difference_in_place_fill_new, 36, sizeof(cl_mem), (void *)&dev_H_new); 
    ezcl_set_kernel_arg(kernel_calc_finite_difference_in_place_fill_new, 37, sizeof(cl_mem), (void *)&dev_U_new); 
    ezcl_set_kernel_arg(kernel_calc_finite_difference_in_place_fill_new, 38, sizeof(cl_mem), (void *)&dev_V_new); 

    ezcl_enqueue_ndrange_kernel(command_queue, kernel_calc_finite_difference_in_place_fill_new, 1, NULL, &global_work_size, &local_work_size, &calc_finite_difference_in_place_fill_new_event);

    ezcl_wait_for_events(1, &calc_finite_difference_in_place_fill_new_event);
    ezcl_event_release(calc_finite_difference_in_place_fill_new_event);

    //gpu_timers[STATE_TIMER_FINITE_DIFFERENCE_PART5] += (long)(cpu_timer_stop(tstart_cpu_part)*1.0e9);
    //cpu_timer_start(&tstart_cpu_part);
/*
   vector<real_t>H_loc(ncells);
   ezcl_enqueue_read_buffer(command_queue, dev_H_new,     CL_TRUE, 0, ncells*sizeof(cl_real_t), &H_loc[0], NULL);
   vector<real_t>U_loc(ncells);
   ezcl_enqueue_read_buffer(command_queue, dev_U_new,     CL_TRUE, 0, ncells*sizeof(cl_real_t), &U_loc[0], NULL);
   vector<real_t>V_loc(ncells);
   ezcl_enqueue_read_buffer(command_queue, dev_V_new,     CL_TRUE, 0, ncells*sizeof(cl_real_t), &V_loc[0], NULL);
   printf("\n");
   for (int jello = 0; jello < ncells; jello++) { printf("%d) %f | %f | %f\n", jello, H_loc[jello], U_loc[jello], V_loc[jello]); }
*/

   /* 
   int nxface;
   ezcl_enqueue_read_buffer(command_queue, dev_nface,     CL_TRUE, 0, sizeof(cl_int), &nxface, NULL);
   vector<int>lefty(nxface);
   ezcl_enqueue_read_buffer(command_queue, dev_map_xface2cell_lower,     CL_TRUE, 0, nxface*sizeof(cl_int), &lefty[0], NULL);
   vector<int>righty(nxface);
   ezcl_enqueue_read_buffer(command_queue, dev_map_xface2cell_upper,     CL_TRUE, 0, nxface*sizeof(cl_int), &righty[0], NULL);
   vector<state_t>Hx_loc(nxface);
   ezcl_enqueue_read_buffer(command_queue, dev_HxFlux,     CL_TRUE, 0, nxface*sizeof(cl_real_t), &Hx_loc[0], NULL);
   vector<state_t>Ux_loc(nxface);
   ezcl_enqueue_read_buffer(command_queue, dev_UxFlux,     CL_TRUE, 0, nxface*sizeof(cl_real_t), &Ux_loc[0], NULL);
   vector<state_t>Vx_loc(nxface);
   ezcl_enqueue_read_buffer(command_queue, dev_VxFlux,     CL_TRUE, 0, nxface*sizeof(cl_real_t), &Vx_loc[0], NULL);
   printf("\n");
   for (int jello = 0; jello < nxface; jello++) { printf("%d) %d %d\n\t%d) %f | %f | %f\n", jello, lefty[jello], righty[jello], jello, Hx_loc[jello], Ux_loc[jello], Vx_loc[jello]); }
   */

     
   dev_H = (cl_mem)gpu_state_memory.memory_replace(dev_H, dev_H_new);
   dev_U = (cl_mem)gpu_state_memory.memory_replace(dev_U, dev_U_new);
   dev_V = (cl_mem)gpu_state_memory.memory_replace(dev_V, dev_V_new);

   gpu_timers[STATE_TIMER_FINITE_DIFFERENCE] += (long)(cpu_timer_stop(tstart_cpu)*1.0e9);
    
   gpu_faces_delete_phantom();
   mesh->gpu_wbidirmap_delete();
}

void State::gpu_calc_finite_difference_via_face_in_place(double deltaT)
{
   //struct timespec tstart_cpu_part;
   //cpu_timer_start(&tstart_cpu_part);

   cl_command_queue command_queue = ezcl_get_command_queue();


   //cl_mem dev_ptr = NULL;

   size_t &ncells    = mesh->ncells;
   size_t &ncells_ghost = mesh->ncells_ghost;
   if (ncells_ghost < ncells) ncells_ghost = ncells;
   int &levmx           = mesh->levmx;
   real_t deltaT_local = deltaT;
   cl_mem &dev_nface    = mesh->dev_nface;
   cl_mem &dev_celltype = mesh->dev_celltype;
   cl_mem &dev_nlft     = mesh->dev_nlft;
   cl_mem &dev_nrht     = mesh->dev_nrht;
   cl_mem &dev_nbot     = mesh->dev_nbot;
   cl_mem &dev_ntop     = mesh->dev_ntop;
   cl_mem &dev_level    = mesh->dev_level;
   cl_mem &dev_levdx    = mesh->dev_levdx;
   cl_mem &dev_levdy    = mesh->dev_levdy;
   cl_mem &dev_xrecvIdx = mesh->dev_xrecvIdx;
   cl_mem &dev_xsendIdx1 = mesh->dev_xsendIdx1;
   cl_mem &dev_xsendIdx2 = mesh->dev_xsendIdx2;
   cl_mem &dev_yrecvIdx = mesh->dev_yrecvIdx;
   cl_mem &dev_ysendIdx1 = mesh->dev_ysendIdx1;
   cl_mem &dev_ysendIdx2 = mesh->dev_ysendIdx2;
   int &nxfixup = mesh->nxfixup;
   int &nyfixup = mesh->nyfixup;

   assert(dev_nface);
   assert(dev_H);
   assert(dev_U);
   assert(dev_V);
   assert(dev_nlft);
   assert(dev_nrht);
   assert(dev_nbot);
   assert(dev_ntop);
   assert(dev_level);
   assert(dev_levdx);
   assert(dev_levdy);

   mesh->gpu_wbidirmap_setup();
    mesh->gpu_calc_face_list_wbidirmap_phantom(gpu_state_memory, deltaT_local);
    gpu_memory_reset_ptrs();

    //gpu_timers[STATE_TIMER_FINITE_DIFFERENCE_PART1] += (long)(cpu_timer_stop(tstart_cpu_part)*1.0e9);
    //cpu_timer_start(&tstart_cpu_part);

   cl_mem dev_H_new = (cl_mem)gpu_state_memory.memory_malloc(ncells, sizeof(cl_state_t), const_cast<char *>("dev_H_new"), DEVICE_REGULAR_MEMORY);
   cl_mem dev_U_new = (cl_mem)gpu_state_memory.memory_malloc(ncells, sizeof(cl_state_t), const_cast<char *>("dev_U_new"), DEVICE_REGULAR_MEMORY);
   cl_mem dev_V_new = (cl_mem)gpu_state_memory.memory_malloc(ncells, sizeof(cl_state_t), const_cast<char *>("dev_V_new"), DEVICE_REGULAR_MEMORY);
 
    //gpu_timers[STATE_TIMER_WRITE] += (long)(cpu_timer_stop(tstart_cpu_part)*1.0e9);
    //cpu_timer_start(&tstart_cpu_part);

   size_t local_work_size = 256;
   size_t global_work_size = ((ncells+local_work_size - 1) /local_work_size) * local_work_size;
   //size_t global_work_size = MAX(MAX(mesh->nxface, mesh->nyface), ncells);

#ifdef HAVE_MPI
   if (mesh->numpe > 1) {
      ezcl_set_kernel_arg(kernel_apply_boundary_conditions_local, 0, sizeof(cl_int), &ncells);
      ezcl_set_kernel_arg(kernel_apply_boundary_conditions_local, 1, sizeof(cl_mem), &dev_celltype);
      ezcl_set_kernel_arg(kernel_apply_boundary_conditions_local, 2, sizeof(cl_mem), &dev_nlft);
      ezcl_set_kernel_arg(kernel_apply_boundary_conditions_local, 3, sizeof(cl_mem), &dev_nrht);
      ezcl_set_kernel_arg(kernel_apply_boundary_conditions_local, 4, sizeof(cl_mem), &dev_ntop);
      ezcl_set_kernel_arg(kernel_apply_boundary_conditions_local, 5, sizeof(cl_mem), &dev_nbot);
      ezcl_set_kernel_arg(kernel_apply_boundary_conditions_local, 6, sizeof(cl_mem), &dev_H);
      ezcl_set_kernel_arg(kernel_apply_boundary_conditions_local, 7, sizeof(cl_mem), &dev_U);
      ezcl_set_kernel_arg(kernel_apply_boundary_conditions_local, 8, sizeof(cl_mem), &dev_V);
      ezcl_enqueue_ndrange_kernel(command_queue, kernel_apply_boundary_conditions_local,   1, NULL, &global_work_size, &local_work_size, NULL);
    
        /*
        __kernel void copy_state_data_cl(
                         const int    isize,         // 0
                __global      state_t *H,            // 1
                __global      state_t *U,            // 2
                __global      state_t *V,            // 3
                __global      state_t *H_new,        // 4
                __global      state_t *U_new,        // 5
                __global      state_t *V_new)        // 6
        */

      ezcl_set_kernel_arg(kernel_copy_state_data, 0, sizeof(cl_int), (void *)&ncells);
      ezcl_set_kernel_arg(kernel_copy_state_data, 1, sizeof(cl_mem), (void *)&dev_H);
      ezcl_set_kernel_arg(kernel_copy_state_data, 2, sizeof(cl_mem), (void *)&dev_U);
      ezcl_set_kernel_arg(kernel_copy_state_data, 3, sizeof(cl_mem), (void *)&dev_V);
      ezcl_set_kernel_arg(kernel_copy_state_data, 4, sizeof(cl_mem), (void *)&dev_H_new);
      ezcl_set_kernel_arg(kernel_copy_state_data, 5, sizeof(cl_mem), (void *)&dev_U_new);
      ezcl_set_kernel_arg(kernel_copy_state_data, 6, sizeof(cl_mem), (void *)&dev_V_new);

      //ezcl_enqueue_ndrange_kernel(command_queue, kernel_copy_state_data,   1, NULL, &global_work_size, &local_work_size, &copy_state_data_event);
      ezcl_enqueue_ndrange_kernel(command_queue, kernel_copy_state_data,   1, NULL, &global_work_size, &local_work_size, NULL);

      dev_H = (cl_mem)gpu_state_memory.memory_replace(dev_H, dev_H_new);
      dev_U = (cl_mem)gpu_state_memory.memory_replace(dev_U, dev_U_new);
      dev_V = (cl_mem)gpu_state_memory.memory_replace(dev_V, dev_V_new);

      L7_Dev_Update(dev_H, L7_STATE_T, mesh->cell_handle);
      L7_Dev_Update(dev_U, L7_STATE_T, mesh->cell_handle);
      L7_Dev_Update(dev_V, L7_STATE_T, mesh->cell_handle);

      dev_H_new = (cl_mem)gpu_state_memory.memory_malloc(ncells_ghost, sizeof(cl_state_t), const_cast<char *>("dev_H_new"), DEVICE_REGULAR_MEMORY);
      dev_U_new = (cl_mem)gpu_state_memory.memory_malloc(ncells_ghost, sizeof(cl_state_t), const_cast<char *>("dev_U_new"), DEVICE_REGULAR_MEMORY);
      dev_V_new = (cl_mem)gpu_state_memory.memory_malloc(ncells_ghost, sizeof(cl_state_t), const_cast<char *>("dev_V_new"), DEVICE_REGULAR_MEMORY);

      ezcl_set_kernel_arg(kernel_apply_boundary_conditions_ghost, 0, sizeof(cl_int), &ncells);
      ezcl_set_kernel_arg(kernel_apply_boundary_conditions_ghost, 1, sizeof(cl_mem), &dev_celltype);
      ezcl_set_kernel_arg(kernel_apply_boundary_conditions_ghost, 2, sizeof(cl_mem), &dev_nlft);
      ezcl_set_kernel_arg(kernel_apply_boundary_conditions_ghost, 3, sizeof(cl_mem), &dev_nrht);
      ezcl_set_kernel_arg(kernel_apply_boundary_conditions_ghost, 4, sizeof(cl_mem), &dev_ntop);
      ezcl_set_kernel_arg(kernel_apply_boundary_conditions_ghost, 5, sizeof(cl_mem), &dev_nbot);
      ezcl_set_kernel_arg(kernel_apply_boundary_conditions_ghost, 6, sizeof(cl_mem), &dev_H);
      ezcl_set_kernel_arg(kernel_apply_boundary_conditions_ghost, 7, sizeof(cl_mem), &dev_U);
      ezcl_set_kernel_arg(kernel_apply_boundary_conditions_ghost, 8, sizeof(cl_mem), &dev_V);
      ezcl_enqueue_ndrange_kernel(command_queue, kernel_apply_boundary_conditions_ghost,   1, NULL, &global_work_size, &local_work_size, NULL);
   } else {
      ezcl_set_kernel_arg(kernel_apply_boundary_conditions, 0, sizeof(cl_int), &ncells);
      ezcl_set_kernel_arg(kernel_apply_boundary_conditions, 1, sizeof(cl_mem), &dev_celltype);
      ezcl_set_kernel_arg(kernel_apply_boundary_conditions, 2, sizeof(cl_mem), &dev_nlft);
      ezcl_set_kernel_arg(kernel_apply_boundary_conditions, 3, sizeof(cl_mem), &dev_nrht);
      ezcl_set_kernel_arg(kernel_apply_boundary_conditions, 4, sizeof(cl_mem), &dev_ntop);
      ezcl_set_kernel_arg(kernel_apply_boundary_conditions, 5, sizeof(cl_mem), &dev_nbot);
      ezcl_set_kernel_arg(kernel_apply_boundary_conditions, 6, sizeof(cl_mem), &dev_H);
      ezcl_set_kernel_arg(kernel_apply_boundary_conditions, 7, sizeof(cl_mem), &dev_U);
      ezcl_set_kernel_arg(kernel_apply_boundary_conditions, 8, sizeof(cl_mem), &dev_V);
      ezcl_enqueue_ndrange_kernel(command_queue, kernel_apply_boundary_conditions,   1, NULL, &global_work_size, &local_work_size, NULL);
   }
#else
   ezcl_set_kernel_arg(kernel_apply_boundary_conditions, 0, sizeof(cl_int), &ncells);
   ezcl_set_kernel_arg(kernel_apply_boundary_conditions, 1, sizeof(cl_mem), (void *)&dev_celltype);
   ezcl_set_kernel_arg(kernel_apply_boundary_conditions, 2, sizeof(cl_mem), (void *)&dev_nlft);
   ezcl_set_kernel_arg(kernel_apply_boundary_conditions, 3, sizeof(cl_mem), (void *)&dev_nrht);
   ezcl_set_kernel_arg(kernel_apply_boundary_conditions, 4, sizeof(cl_mem), (void *)&dev_ntop);
   ezcl_set_kernel_arg(kernel_apply_boundary_conditions, 5, sizeof(cl_mem), (void *)&dev_nbot);
   ezcl_set_kernel_arg(kernel_apply_boundary_conditions, 6, sizeof(cl_mem), (void *)&dev_H);
   ezcl_set_kernel_arg(kernel_apply_boundary_conditions, 7, sizeof(cl_mem), (void *)&dev_U);
   ezcl_set_kernel_arg(kernel_apply_boundary_conditions, 8, sizeof(cl_mem), (void *)&dev_V);

    //gpu_timers[STATE_TIMER_WRITE] += (long)(cpu_timer_stop(tstart_cpu_part)*1.0e9);
    //cpu_timer_start(&tstart_cpu_part);

   ezcl_enqueue_ndrange_kernel(command_queue, kernel_apply_boundary_conditions,   1, NULL, &global_work_size, &local_work_size, NULL);
#endif

    //gpu_timers[STATE_TIMER_FINITE_DIFFERENCE_PART2] += (long)(cpu_timer_stop(tstart_cpu_part)*1.0e9);
    //cpu_timer_start(&tstart_cpu_part);

   cl_mem &dev_map_xface2cell_lower = mesh->dev_map_xface2cell_lower;
   cl_mem &dev_map_xface2cell_upper = mesh->dev_map_xface2cell_upper;
   cl_mem &dev_map_xcell2face_left1 = mesh->dev_map_xcell2face_left1;
   cl_mem &dev_map_xcell2face_left2 = mesh->dev_map_xcell2face_left2;
   cl_mem &dev_map_xcell2face_right1 = mesh->dev_map_xcell2face_right1;
   cl_mem &dev_map_xcell2face_right2 = mesh->dev_map_xcell2face_right2;
   cl_mem &dev_map_yface2cell_lower = mesh->dev_map_yface2cell_lower;
   cl_mem &dev_map_yface2cell_upper = mesh->dev_map_yface2cell_upper;
   cl_mem &dev_map_ycell2face_bot1 = mesh->dev_map_ycell2face_bot1;
   cl_mem &dev_map_ycell2face_bot2 = mesh->dev_map_ycell2face_bot2;
   cl_mem &dev_map_ycell2face_top1 = mesh->dev_map_ycell2face_top1;
   cl_mem &dev_map_ycell2face_top2 = mesh->dev_map_ycell2face_top2;
    cl_mem dev_xface_level = mesh->dev_xface_level;
    cl_mem dev_xface_i = mesh->dev_xface_i;
    cl_mem dev_xface_j = mesh->dev_xface_j;
    cl_mem dev_ixmin_level = mesh->dev_ixmin_level;
    cl_mem dev_ixmax_level = mesh->dev_ixmax_level;
    cl_mem dev_jxmin_level = mesh->dev_jxmin_level;
    cl_mem dev_jxmax_level = mesh->dev_jxmax_level;
    cl_mem dev_yface_level = mesh->dev_yface_level;
    cl_mem dev_yface_i = mesh->dev_yface_i;
    cl_mem dev_yface_j = mesh->dev_yface_j;
    cl_mem dev_iymin_level = mesh->dev_iymin_level;
    cl_mem dev_iymax_level = mesh->dev_iymax_level;
    cl_mem dev_jymin_level = mesh->dev_jymin_level;
    cl_mem dev_jymax_level = mesh->dev_jymax_level;
   assert(dev_map_xface2cell_lower);
   assert(dev_map_xface2cell_upper);
   assert(dev_map_xcell2face_left1);
   assert(dev_map_xcell2face_left2);
   assert(dev_map_xcell2face_right1);
   assert(dev_map_xcell2face_right2);
   assert(dev_map_yface2cell_lower);
   assert(dev_map_yface2cell_upper);
   assert(dev_map_ycell2face_bot1);
   assert(dev_map_ycell2face_bot2);
   assert(dev_map_ycell2face_top1);
   assert(dev_map_ycell2face_top2);
   assert(dev_xface_level);
   assert(dev_xface_i);
   assert(dev_xface_j);
   assert(dev_ixmin_level);
   assert(dev_ixmax_level);
   assert(dev_jxmin_level);
   assert(dev_jxmax_level);
   assert(dev_yface_level);
   assert(dev_yface_i);
   assert(dev_yface_j);
   assert(dev_iymin_level);
   assert(dev_iymax_level);
   assert(dev_jymin_level);
   assert(dev_jymax_level);

    size_t mem_requestx, mem_requesty;
    //size_t nfacex, nfacey;
    //ezcl_enqueue_read_buffer(command_queue, dev_nface,     CL_TRUE, 0, sizeof(cl_int), &nfacex, NULL);
    //ezcl_enqueue_read_buffer(command_queue, dev_nface,     CL_TRUE, 1*sizeof(cl_int), sizeof(cl_int), &nfacey, NULL);
    //printf("\nMem requests %d and %d\n", nfacex, nfacey);
    mem_requestx = mesh->pxfaceCnt;
    mem_requesty = mesh->pyfaceCnt;
    //printf("\n%d %d\n", mem_requestx, mem_requesty);
    
    gpu_faces_setup(mem_requestx, mem_requesty);

   cl_event calc_finite_difference_in_place_face_event, calc_finite_difference_in_place_fill_new_event, calc_finite_difference_in_place_fixup_event;

   //size_t local_face_work = CL_DEVICE_MAX_WORK_GROUP_SIZE;
   //size_t global_face_work = ((pcellCnt+local_face_work - 1) /local_face_work) * local_face_work;
   //printf("\nglobal face work %d\n", global_face_work);

   struct timespec tstart_cpu;
   cpu_timer_start(&tstart_cpu);

   /* 
   vector<int>H_loc(p);
   //ezcl_enqueue_read_buffer(command_queue, dev_map_xface2cell_lower,     CL_TRUE, 0, mem_requestx*sizeof(cl_int), &H_loc[0], NULL);
   ezcl_enqueue_read_buffer(command_queue, dev_map_xface2cell_lower,     CL_TRUE, 0, mem_requestx*sizeof(cl_int), &H_loc[0], NULL);
   vector<int>U_loc(mem_requestx);
   //ezcl_enqueue_read_buffer(command_queue, dev_map_xface2cell_upper,     CL_TRUE, 0, mem_requestx*sizeof(cl_int), &U_loc[0], NULL);
   ezcl_enqueue_read_buffer(command_queue, dev_map_xface2cell_upper,     CL_TRUE, 0, mem_requestx*sizeof(cl_int), &U_loc[0], NULL);
   vector<int>V_loc(mem_requesty);
   //ezcl_enqueue_read_buffer(command_queue, dev_map_yface2cell_lower,     CL_TRUE, 0, mem_requesty*sizeof(cl_int), &V_loc[0], NULL);
   ezcl_enqueue_read_buffer(command_queue, dev_map_yface2cell_lower,     CL_TRUE, 0, mem_requesty*sizeof(cl_int), &V_loc[0], NULL);
   vector<int>Z_loc(mem_requesty);
   //ezcl_enqueue_read_buffer(command_queue, dev_map_yface2cell_upper,     CL_TRUE, 0, mem_requesty*sizeof(cl_int), &Z_loc[0], NULL);
   ezcl_enqueue_read_buffer(command_queue, dev_map_yface2cell_upper,     CL_TRUE, 0, mem_requesty*sizeof(cl_int), &Z_loc[0], NULL);
   printf("\n");
   for (int jello = 0; jello < mem_requestx; jello++) { printf("%d) %d | %d | %d | %d\n", jello, H_loc[jello], U_loc[jello], V_loc[jello], Z_loc[jello]); }
*/
/*
   vector<real_t>H_loc(pcellCnt);
   ezcl_enqueue_read_buffer(command_queue, dev_H,     CL_TRUE, 0, pcellCnt*sizeof(cl_real_t), &H_loc[0], NULL);
   vector<real_t>U_loc(pcellCnt);
   ezcl_enqueue_read_buffer(command_queue, dev_U,     CL_TRUE, 0, pcellCnt*sizeof(cl_real_t), &U_loc[0], NULL);
   vector<real_t>V_loc(pcellCnt);
   ezcl_enqueue_read_buffer(command_queue, dev_V,     CL_TRUE, 0, pcellCnt*sizeof(cl_real_t), &V_loc[0], NULL);
   printf("\n");
   for (int jello = 0; jello < pcellCnt; jello++) { printf("%d) %f | %f | %f\n", jello, H_loc[jello], U_loc[jello], V_loc[jello]); }
*/

    /*
__kernel void calc_finite_difference_via_face_in_place_face_comps_cl(
                        const int       ncells,                     // 0 Number of cells (not including phantom)
            __global    const int       *nfaces,                    // 1 Number of faces
                        const int       levmx,                      // 2 Maximum level
            __global    const state_t   *H,                         // 3
            __global    const state_t   *U,                         // 4
            __global    const state_t   *V,                         // 5
            __global    const uchar_t   *level,                     // 6 Array of level information
                        const real_t    deltaT,                     // 7 Size of time step
            __global    const real_t    *lev_dx,                    // 8
            __global    const real_t    *lev_dy,                    // 9
            __global    const int       *map_xface2cell_lower,      // 10 A face's left cell 
            __global    const int       *map_xface2cell_upper,      // 11 A face's left cell 
            __global    const int       *map_yface2cell_lower,      // 12 A face's below cell 
            __global    const int       *map_yface2cell_upper,      // 13 A face's above cell 
            __global    const int       *map_xcell2face_left1,      // 14 A cell's left primary face 
            __global    const int       *map_xcell2face_right1,     // 15 A cell's right primary face 
            __global    const int       *map_ycell2face_bot1,       // 16 A cell's bot primary face 
            __global    const int       *map_ycell2face_top1,       // 17 A cell's top primary face 
            __global          state_t   *HxFlux,                    // 18
            __global          state_t   *UxFlux,                    // 19
            __global          state_t   *VxFlux,                    // 20
            __global          state_t   *HyFlux,                    // 21
            __global          state_t   *UyFlux,                    // 22
            __global          state_t   *VyFlux,                    // 23
            __global          state_t   *Wx_H,                      // 24
            __global          state_t   *Wx_U,                      // 25
            __global          state_t   *Wy_H,                      // 26
            __global          state_t   *Wy_V) {                    // 27

            */

   ezcl_set_kernel_arg(kernel_calc_finite_difference_via_face_in_place_face_comps, 0, sizeof(cl_int), (void *)&ncells); 
   ezcl_set_kernel_arg(kernel_calc_finite_difference_via_face_in_place_face_comps, 1, sizeof(cl_mem), (void *)&dev_nface); 
   ezcl_set_kernel_arg(kernel_calc_finite_difference_via_face_in_place_face_comps, 2, sizeof(cl_int), (void *)&levmx); 
   ezcl_set_kernel_arg(kernel_calc_finite_difference_via_face_in_place_face_comps, 3, sizeof(cl_mem), (void *)&dev_H); 
   ezcl_set_kernel_arg(kernel_calc_finite_difference_via_face_in_place_face_comps, 4, sizeof(cl_mem), (void *)&dev_U); 
   ezcl_set_kernel_arg(kernel_calc_finite_difference_via_face_in_place_face_comps, 5, sizeof(cl_mem), (void *)&dev_V); 
   ezcl_set_kernel_arg(kernel_calc_finite_difference_via_face_in_place_face_comps, 6, sizeof(cl_mem), (void *)&dev_level); 
   ezcl_set_kernel_arg(kernel_calc_finite_difference_via_face_in_place_face_comps, 7, sizeof(cl_real_t), (void *)&deltaT_local); 
   ezcl_set_kernel_arg(kernel_calc_finite_difference_via_face_in_place_face_comps, 8, sizeof(cl_mem), (void *)&dev_levdx); 
   ezcl_set_kernel_arg(kernel_calc_finite_difference_via_face_in_place_face_comps, 9, sizeof(cl_mem), (void *)&dev_levdy); 
   ezcl_set_kernel_arg(kernel_calc_finite_difference_via_face_in_place_face_comps, 10, sizeof(cl_mem), (void *)&dev_map_xface2cell_lower); 
   ezcl_set_kernel_arg(kernel_calc_finite_difference_via_face_in_place_face_comps, 11, sizeof(cl_mem), (void *)&dev_map_xface2cell_upper); 
   ezcl_set_kernel_arg(kernel_calc_finite_difference_via_face_in_place_face_comps, 12, sizeof(cl_mem), (void *)&dev_map_yface2cell_lower); 
   ezcl_set_kernel_arg(kernel_calc_finite_difference_via_face_in_place_face_comps, 13, sizeof(cl_mem), (void *)&dev_map_yface2cell_upper); 
   ezcl_set_kernel_arg(kernel_calc_finite_difference_via_face_in_place_face_comps, 14, sizeof(cl_mem), (void *)&dev_map_xcell2face_left1); 
   ezcl_set_kernel_arg(kernel_calc_finite_difference_via_face_in_place_face_comps, 15, sizeof(cl_mem), (void *)&dev_map_xcell2face_right1); 
   ezcl_set_kernel_arg(kernel_calc_finite_difference_via_face_in_place_face_comps, 16, sizeof(cl_mem), (void *)&dev_map_ycell2face_bot1); 
   ezcl_set_kernel_arg(kernel_calc_finite_difference_via_face_in_place_face_comps, 17, sizeof(cl_mem), (void *)&dev_map_ycell2face_top1); 
   ezcl_set_kernel_arg(kernel_calc_finite_difference_via_face_in_place_face_comps, 18, sizeof(cl_mem), (void *)&dev_HxFlux); 
   ezcl_set_kernel_arg(kernel_calc_finite_difference_via_face_in_place_face_comps, 19, sizeof(cl_mem), (void *)&dev_UxFlux); 
   ezcl_set_kernel_arg(kernel_calc_finite_difference_via_face_in_place_face_comps, 20, sizeof(cl_mem), (void *)&dev_VxFlux); 
   ezcl_set_kernel_arg(kernel_calc_finite_difference_via_face_in_place_face_comps, 21, sizeof(cl_mem), (void *)&dev_HyFlux); 
   ezcl_set_kernel_arg(kernel_calc_finite_difference_via_face_in_place_face_comps, 22, sizeof(cl_mem), (void *)&dev_UyFlux); 
   ezcl_set_kernel_arg(kernel_calc_finite_difference_via_face_in_place_face_comps, 23, sizeof(cl_mem), (void *)&dev_VyFlux); 
   ezcl_set_kernel_arg(kernel_calc_finite_difference_via_face_in_place_face_comps, 24, sizeof(cl_mem), (void *)&dev_Wx_H); 
   ezcl_set_kernel_arg(kernel_calc_finite_difference_via_face_in_place_face_comps, 25, sizeof(cl_mem), (void *)&dev_Wx_U); 
   ezcl_set_kernel_arg(kernel_calc_finite_difference_via_face_in_place_face_comps, 26, sizeof(cl_mem), (void *)&dev_Wy_H); 
   ezcl_set_kernel_arg(kernel_calc_finite_difference_via_face_in_place_face_comps, 27, sizeof(cl_mem), (void *)&dev_Wy_V); 

    //gpu_timers[STATE_TIMER_WRITE] += (long)(cpu_timer_stop(tstart_cpu_part)*1.0e9);
    //cpu_timer_start(&tstart_cpu_part);

   ezcl_enqueue_ndrange_kernel(command_queue, kernel_calc_finite_difference_via_face_in_place_face_comps, 1, NULL, &global_work_size, &local_work_size, &calc_finite_difference_in_place_face_event);

   ezcl_wait_for_events(1, &calc_finite_difference_in_place_face_event);
   ezcl_event_release(calc_finite_difference_in_place_face_event);
   
    //gpu_timers[STATE_TIMER_FINITE_DIFFERENCE_PART3] += (long)(cpu_timer_stop(tstart_cpu_part)*1.0e9);
    //cpu_timer_start(&tstart_cpu_part);

   /*
__kernel void calc_finite_difference_via_face_in_place_fixup_cl(
                        const int        nxfixup,                   // 0
                        const int        nyfixup,                   // 1
            __global    const int       *xrecvIdx,                  // 2
            __global    const int       *xsendIdx1,                 // 3
            __global    const int       *xsendIdx2,                 // 4
            __global    const int       *yrecvIdx,                  // 5
            __global    const int       *ysendIdx1,                 // 6
            __global    const int       *ysendIdx2,                 // 7
            __global    const int       *map_xface2cell_lower,      // 8
            __global    const int       *map_xface2cell_upper,      // 9
            __global    const int       *map_yface2cell_lower,      // 10 
            __global    const int       *map_yface2cell_upper,      // 11 
            __global          state_t   *HxFlux,                    // 12
            __global          state_t   *UxFlux,                    // 13
            __global          state_t   *VxFlux,                    // 14
            __global          state_t   *HyFlux,                    // 15
            __global          state_t   *UyFlux,                    // 16
            __global          state_t   *VyFlux,                    // 17
            __global          state_t   *Wx_H,                      // 18
            __global          state_t   *Wx_U,                      // 19
            __global          state_t   *Wy_H,                      // 20
            __global          state_t   *Wy_V) {                    // 21
    */

    ezcl_set_kernel_arg(kernel_calc_finite_difference_via_face_in_place_fixup, 0, sizeof(cl_int), (void *)&nxfixup); 
    ezcl_set_kernel_arg(kernel_calc_finite_difference_via_face_in_place_fixup, 1, sizeof(cl_int), (void *)&nyfixup); 
    ezcl_set_kernel_arg(kernel_calc_finite_difference_via_face_in_place_fixup, 2, sizeof(cl_mem), (void *)&dev_xrecvIdx);
    ezcl_set_kernel_arg(kernel_calc_finite_difference_via_face_in_place_fixup, 3, sizeof(cl_mem), (void *)&dev_xsendIdx1);
    ezcl_set_kernel_arg(kernel_calc_finite_difference_via_face_in_place_fixup, 4, sizeof(cl_mem), (void *)&dev_xsendIdx2);
    ezcl_set_kernel_arg(kernel_calc_finite_difference_via_face_in_place_fixup, 5, sizeof(cl_mem), (void *)&dev_yrecvIdx);
    ezcl_set_kernel_arg(kernel_calc_finite_difference_via_face_in_place_fixup, 6, sizeof(cl_mem), (void *)&dev_ysendIdx1);
    ezcl_set_kernel_arg(kernel_calc_finite_difference_via_face_in_place_fixup, 7, sizeof(cl_mem), (void *)&dev_ysendIdx2);
    ezcl_set_kernel_arg(kernel_calc_finite_difference_via_face_in_place_fixup, 8, sizeof(cl_mem), (void *)&dev_map_xface2cell_lower); 
    ezcl_set_kernel_arg(kernel_calc_finite_difference_via_face_in_place_fixup, 9, sizeof(cl_mem), (void *)&dev_map_xface2cell_upper); 
    ezcl_set_kernel_arg(kernel_calc_finite_difference_via_face_in_place_fixup, 10, sizeof(cl_mem), (void *)&dev_map_yface2cell_lower); 
    ezcl_set_kernel_arg(kernel_calc_finite_difference_via_face_in_place_fixup, 11, sizeof(cl_mem), (void *)&dev_map_yface2cell_upper); 
    ezcl_set_kernel_arg(kernel_calc_finite_difference_via_face_in_place_fixup, 12, sizeof(cl_mem), (void *)&dev_HxFlux); 
    ezcl_set_kernel_arg(kernel_calc_finite_difference_via_face_in_place_fixup, 13, sizeof(cl_mem), (void *)&dev_UxFlux); 
    ezcl_set_kernel_arg(kernel_calc_finite_difference_via_face_in_place_fixup, 14, sizeof(cl_mem), (void *)&dev_VxFlux); 
    ezcl_set_kernel_arg(kernel_calc_finite_difference_via_face_in_place_fixup, 15, sizeof(cl_mem), (void *)&dev_HyFlux); 
    ezcl_set_kernel_arg(kernel_calc_finite_difference_via_face_in_place_fixup, 16, sizeof(cl_mem), (void *)&dev_UyFlux); 
    ezcl_set_kernel_arg(kernel_calc_finite_difference_via_face_in_place_fixup, 17, sizeof(cl_mem), (void *)&dev_VyFlux); 
    ezcl_set_kernel_arg(kernel_calc_finite_difference_via_face_in_place_fixup, 18, sizeof(cl_mem), (void *)&dev_Wx_H); 
    ezcl_set_kernel_arg(kernel_calc_finite_difference_via_face_in_place_fixup, 19, sizeof(cl_mem), (void *)&dev_Wx_U); 
    ezcl_set_kernel_arg(kernel_calc_finite_difference_via_face_in_place_fixup, 20, sizeof(cl_mem), (void *)&dev_Wy_H); 
    ezcl_set_kernel_arg(kernel_calc_finite_difference_via_face_in_place_fixup, 21, sizeof(cl_mem), (void *)&dev_Wy_V); 

    //gpu_timers[STATE_TIMER_WRITE] += (long)(cpu_timer_stop(tstart_cpu_part)*1.0e9);
    //cpu_timer_start(&tstart_cpu_part);

    ezcl_enqueue_ndrange_kernel(command_queue, kernel_calc_finite_difference_via_face_in_place_fixup, 1, NULL, &global_work_size, &local_work_size, &calc_finite_difference_in_place_fixup_event);

    ezcl_wait_for_events(1, &calc_finite_difference_in_place_fixup_event);
    ezcl_event_release(calc_finite_difference_in_place_fixup_event);

    //gpu_timers[STATE_TIMER_FINITE_DIFFERENCE_PART4] += (long)(cpu_timer_stop(tstart_cpu_part)*1.0e9);
    //cpu_timer_start(&tstart_cpu_part);

/*
   vector<real_t>H_loc(mem_requestx);
   ezcl_enqueue_read_buffer(command_queue, dev_HxFlux,     CL_TRUE, 0, mem_requestx*sizeof(cl_real_t), &H_loc[0], NULL);
   vector<real_t>U_loc(mem_requestx);
   ezcl_enqueue_read_buffer(command_queue, dev_UxFlux,     CL_TRUE, 0, mem_requestx*sizeof(cl_real_t), &U_loc[0], NULL);
   vector<real_t>V_loc(mem_requestx);
   ezcl_enqueue_read_buffer(command_queue, dev_VxFlux,     CL_TRUE, 0, mem_requestx*sizeof(cl_real_t), &V_loc[0], NULL);
   printf("\n");
   for (int jello = 0; jello < mem_requestx; jello++) { printf("%d) %f | %f | %f\n", jello, H_loc[jello], U_loc[jello], V_loc[jello]); }
*/


/*
__kernel void calc_finite_difference_via_face_in_place_fill_new_cl(
                        const int       ncells,                     // 0 Number of cells (not including phantom)
                        const real_t    deltaT,                     // 1 Size of time step
            __global    const real_t    *lev_dx,                    // 2
            __global    const real_t    *lev_dy,                    // 3
            __global    const int       *level,                     // 4
            __global    const int       *map_xcell2face_left1,      // 5 A cell's left primary face 
            __global    const int       *map_xcell2face_right1,     // 6 A cell's right primary face 
            __global    const int       *map_ycell2face_bot1,       // 7 A cell's bot primary face 
            __global    const int       *map_ycell2face_top1,       // 8 A cell's top primary face 
            __global          state_t   *HxFlux,                    // 9
            __global          state_t   *UxFlux,                    // 10
            __global          state_t   *VxFlux,                    // 11
            __global          state_t   *HyFlux,                    // 12
            __global          state_t   *UyFlux,                    // 13
            __global          state_t   *VyFlux,                    // 14
            __global          state_t   *Wx_H,                      // 15
            __global          state_t   *Wx_U,                      // 16
            __global          state_t   *Wy_H,                      // 17
            __global          state_t   *Wy_V,                      // 18
            __global    const state_t   *H,                         // 19
            __global    const state_t   *U,                         // 20
            __global    const state_t   *V,                         // 21
            __global          state_t   *H_new,                     // 22
            __global          state_t   *U_new,                     // 23
            __global          state_t   *V_new) {                   // 24
*/
  

    ezcl_set_kernel_arg(kernel_calc_finite_difference_via_face_in_place_fill_new, 0, sizeof(cl_int), (void *)&ncells); 
    ezcl_set_kernel_arg(kernel_calc_finite_difference_via_face_in_place_fill_new, 1, sizeof(cl_real_t), (void *)&deltaT_local); 
    ezcl_set_kernel_arg(kernel_calc_finite_difference_via_face_in_place_fill_new, 2, sizeof(cl_mem), (void *)&dev_levdx); 
    ezcl_set_kernel_arg(kernel_calc_finite_difference_via_face_in_place_fill_new, 3, sizeof(cl_mem), (void *)&dev_levdy); 
    ezcl_set_kernel_arg(kernel_calc_finite_difference_via_face_in_place_fill_new, 4, sizeof(cl_mem), (void *)&dev_level); 
    ezcl_set_kernel_arg(kernel_calc_finite_difference_via_face_in_place_fill_new, 5, sizeof(cl_mem), (void *)&dev_map_xcell2face_left1); 
    ezcl_set_kernel_arg(kernel_calc_finite_difference_via_face_in_place_fill_new, 6, sizeof(cl_mem), (void *)&dev_map_xcell2face_right1); 
    ezcl_set_kernel_arg(kernel_calc_finite_difference_via_face_in_place_fill_new, 7, sizeof(cl_mem), (void *)&dev_map_ycell2face_bot1); 
    ezcl_set_kernel_arg(kernel_calc_finite_difference_via_face_in_place_fill_new, 8, sizeof(cl_mem), (void *)&dev_map_ycell2face_top1); 
    ezcl_set_kernel_arg(kernel_calc_finite_difference_via_face_in_place_fill_new, 9, sizeof(cl_mem), (void *)&dev_HxFlux); 
    ezcl_set_kernel_arg(kernel_calc_finite_difference_via_face_in_place_fill_new, 10, sizeof(cl_mem), (void *)&dev_UxFlux); 
    ezcl_set_kernel_arg(kernel_calc_finite_difference_via_face_in_place_fill_new, 11, sizeof(cl_mem), (void *)&dev_VxFlux); 
    ezcl_set_kernel_arg(kernel_calc_finite_difference_via_face_in_place_fill_new, 12, sizeof(cl_mem), (void *)&dev_HyFlux); 
    ezcl_set_kernel_arg(kernel_calc_finite_difference_via_face_in_place_fill_new, 13, sizeof(cl_mem), (void *)&dev_UyFlux); 
    ezcl_set_kernel_arg(kernel_calc_finite_difference_via_face_in_place_fill_new, 14, sizeof(cl_mem), (void *)&dev_VyFlux); 
    ezcl_set_kernel_arg(kernel_calc_finite_difference_via_face_in_place_fill_new, 15, sizeof(cl_mem), (void *)&dev_Wx_H); 
    ezcl_set_kernel_arg(kernel_calc_finite_difference_via_face_in_place_fill_new, 16, sizeof(cl_mem), (void *)&dev_Wx_U); 
    ezcl_set_kernel_arg(kernel_calc_finite_difference_via_face_in_place_fill_new, 17, sizeof(cl_mem), (void *)&dev_Wy_H); 
    ezcl_set_kernel_arg(kernel_calc_finite_difference_via_face_in_place_fill_new, 18, sizeof(cl_mem), (void *)&dev_Wy_V); 
    ezcl_set_kernel_arg(kernel_calc_finite_difference_via_face_in_place_fill_new, 19, sizeof(cl_mem), (void *)&dev_H); 
    ezcl_set_kernel_arg(kernel_calc_finite_difference_via_face_in_place_fill_new, 20, sizeof(cl_mem), (void *)&dev_U); 
    ezcl_set_kernel_arg(kernel_calc_finite_difference_via_face_in_place_fill_new, 21, sizeof(cl_mem), (void *)&dev_V); 
    ezcl_set_kernel_arg(kernel_calc_finite_difference_via_face_in_place_fill_new, 22, sizeof(cl_mem), (void *)&dev_H_new); 
    ezcl_set_kernel_arg(kernel_calc_finite_difference_via_face_in_place_fill_new, 23, sizeof(cl_mem), (void *)&dev_U_new); 
    ezcl_set_kernel_arg(kernel_calc_finite_difference_via_face_in_place_fill_new, 24, sizeof(cl_mem), (void *)&dev_V_new); 

    //gpu_timers[STATE_TIMER_WRITE] += (long)(cpu_timer_stop(tstart_cpu_part)*1.0e9);
    //cpu_timer_start(&tstart_cpu_part);

    ezcl_enqueue_ndrange_kernel(command_queue, kernel_calc_finite_difference_via_face_in_place_fill_new, 1, NULL, &global_work_size, &local_work_size, &calc_finite_difference_in_place_fill_new_event);

    ezcl_wait_for_events(1, &calc_finite_difference_in_place_fill_new_event);
    ezcl_event_release(calc_finite_difference_in_place_fill_new_event);

    //gpu_timers[STATE_TIMER_FINITE_DIFFERENCE_PART5] += (long)(cpu_timer_stop(tstart_cpu_part)*1.0e9);
    //cpu_timer_start(&tstart_cpu_part);


   /* 
   int nxface;
   ezcl_enqueue_read_buffer(command_queue, dev_nface,     CL_TRUE, 0, sizeof(cl_int), &nxface, NULL);
   vector<int>lefty(nxface);
   ezcl_enqueue_read_buffer(command_queue, dev_map_xface2cell_lower,     CL_TRUE, 0, nxface*sizeof(cl_int), &lefty[0], NULL);
   vector<int>righty(nxface);
   ezcl_enqueue_read_buffer(command_queue, dev_map_xface2cell_upper,     CL_TRUE, 0, nxface*sizeof(cl_int), &righty[0], NULL);
   vector<state_t>Hx_loc(nxface);
   ezcl_enqueue_read_buffer(command_queue, dev_HxFlux,     CL_TRUE, 0, nxface*sizeof(cl_real_t), &Hx_loc[0], NULL);
   vector<state_t>Ux_loc(nxface);
   ezcl_enqueue_read_buffer(command_queue, dev_UxFlux,     CL_TRUE, 0, nxface*sizeof(cl_real_t), &Ux_loc[0], NULL);
   vector<state_t>Vx_loc(nxface);
   ezcl_enqueue_read_buffer(command_queue, dev_VxFlux,     CL_TRUE, 0, nxface*sizeof(cl_real_t), &Vx_loc[0], NULL);
   printf("\n");
   for (int jello = 0; jello < nxface; jello++) { printf("%d) %d %d\n\t%d) %f | %f | %f\n", jello, lefty[jello], righty[jello], jello, Hx_loc[jello], Ux_loc[jello], Vx_loc[jello]); }
   */

     
   dev_H = (cl_mem)gpu_state_memory.memory_replace(dev_H, dev_H_new);
   dev_U = (cl_mem)gpu_state_memory.memory_replace(dev_U, dev_U_new);
   dev_V = (cl_mem)gpu_state_memory.memory_replace(dev_V, dev_V_new);

   gpu_timers[STATE_TIMER_FINITE_DIFFERENCE] += (long)(cpu_timer_stop(tstart_cpu)*1.0e9);
   
    
   gpu_faces_delete();
   mesh->gpu_wbidirmap_only_essentials();
   mesh->gpu_wbidirmap_delete();
}

void State::gpu_reggrid_setup(size_t mem_request)
{
    size_t level_size = 3;
    dev_H_reg_lev = ezcl_malloc(NULL, const_cast<char *>("dev_H_reg_lev"), &mem_request, sizeof(cl_real_t), CL_MEM_READ_WRITE, 0);
    dev_U_reg_lev = ezcl_malloc(NULL, const_cast<char *>("dev_U_reg_lev"), &mem_request, sizeof(cl_real_t), CL_MEM_READ_WRITE, 0);
    dev_V_reg_lev = ezcl_malloc(NULL, const_cast<char *>("dev_V_reg_lev"), &mem_request, sizeof(cl_real_t), CL_MEM_READ_WRITE, 0);
    dev_lev_jregmin = ezcl_malloc(NULL, const_cast<char *>("dev_lev_jregmin"), &level_size, sizeof(cl_int), CL_MEM_READ_WRITE, 0);
    dev_lev_iregmin = ezcl_malloc(NULL, const_cast<char *>("dev_lev_iregmin"), &level_size, sizeof(cl_int), CL_MEM_READ_WRITE, 0);
    dev_lev_jregsize = ezcl_malloc(NULL, const_cast<char *>("dev_lev_jregsize"), &level_size, sizeof(cl_int), CL_MEM_READ_WRITE, 0);
    dev_lev_iregsize = ezcl_malloc(NULL, const_cast<char *>("dev_lev_iregsize"), &level_size, sizeof(cl_int), CL_MEM_READ_WRITE, 0);
    dev_reg_start = ezcl_malloc(NULL, const_cast<char *>("dev_reg_start"), &level_size, sizeof(cl_int), CL_MEM_READ_WRITE, 0);
}

void State::gpu_reggrid_delete()
{
    ezcl_device_memory_delete(dev_H_reg_lev);
    ezcl_device_memory_delete(dev_U_reg_lev);
    ezcl_device_memory_delete(dev_V_reg_lev);
    ezcl_device_memory_delete(dev_lev_jregmin);
    ezcl_device_memory_delete(dev_lev_iregmin);
    ezcl_device_memory_delete(dev_lev_jregsize);
    ezcl_device_memory_delete(dev_lev_iregsize);
    ezcl_device_memory_delete(dev_reg_start);
}

void State::gpu_calc_finite_difference_regular_cells(double deltaT)
{

   //struct timespec tstart_cpu_part;
   //cpu_timer_start(&tstart_cpu_part);

   cl_command_queue command_queue = ezcl_get_command_queue();


   //cl_mem dev_ptr = NULL;

   size_t &ncells    = mesh->ncells;
   size_t &ncells_ghost = mesh->ncells_ghost;
   if (ncells_ghost < ncells) ncells_ghost = ncells;
   int &levmx           = mesh->levmx;
   real_t deltaT_local = deltaT;

   mesh->gpu_wbidirmap_setup();
    mesh->gpu_calc_face_list_wbidirmap_phantom(gpu_state_memory, deltaT_local);
    gpu_memory_reset_ptrs();

   mesh->generate_regular_cell_meshes(gpu_state_memory);

   cl_mem &dev_nface    = mesh->dev_nface;
   cl_mem &dev_celltype = mesh->dev_celltype;
   cl_mem &dev_nlft     = mesh->dev_nlft;
   cl_mem &dev_nrht     = mesh->dev_nrht;
   cl_mem &dev_nbot     = mesh->dev_nbot;
   cl_mem &dev_ntop     = mesh->dev_ntop;
   cl_mem &dev_level    = mesh->dev_level;
   cl_mem &dev_levdx    = mesh->dev_levdx;
   cl_mem &dev_levdy    = mesh->dev_levdy;
   cl_mem &dev_j        = mesh->dev_j;
   cl_mem &dev_i        = mesh->dev_i;
   cl_mem &dev_map_xface2cell_lower = mesh->dev_map_xface2cell_lower;
   cl_mem &dev_map_xface2cell_upper = mesh->dev_map_xface2cell_upper;
   cl_mem &dev_map_yface2cell_lower = mesh->dev_map_yface2cell_lower;
   cl_mem &dev_map_yface2cell_upper = mesh->dev_map_yface2cell_upper;
   cl_mem &dev_xrecvCIdx = mesh->dev_xrecvCIdx;
   cl_mem &dev_xminusCell2Idx = mesh->dev_xminusCell2Idx;
   cl_mem &dev_xplusCell2Idx = mesh->dev_xplusCell2Idx;
   cl_mem &dev_xsendIdx1 = mesh->dev_xsendIdx1;
   cl_mem &dev_xsendIdx2 = mesh->dev_xsendIdx2;
   cl_mem &dev_yrecvCIdx = mesh->dev_yrecvCIdx;
   cl_mem &dev_yminusCell2Idx = mesh->dev_yminusCell2Idx;
   cl_mem &dev_yplusCell2Idx = mesh->dev_yplusCell2Idx;
   cl_mem &dev_ysendIdx1 = mesh->dev_ysendIdx1;
   cl_mem &dev_ysendIdx2 = mesh->dev_ysendIdx2;
   int &nxfixup = mesh->nxfixup;
   int &nyfixup = mesh->nyfixup;

   assert(dev_nface);
   assert(dev_H);
   assert(dev_U);
   assert(dev_V);
   assert(dev_nlft);
   assert(dev_nrht);
   assert(dev_nbot);
   assert(dev_ntop);
   assert(dev_level);
   assert(dev_levdx);
   assert(dev_levdy);

   cl_mem dev_H_new = (cl_mem)gpu_state_memory.memory_malloc(ncells, sizeof(cl_state_t), const_cast<char *>("dev_H_new"), DEVICE_REGULAR_MEMORY);
   cl_mem dev_U_new = (cl_mem)gpu_state_memory.memory_malloc(ncells, sizeof(cl_state_t), const_cast<char *>("dev_U_new"), DEVICE_REGULAR_MEMORY);
   cl_mem dev_V_new = (cl_mem)gpu_state_memory.memory_malloc(ncells, sizeof(cl_state_t), const_cast<char *>("dev_V_new"), DEVICE_REGULAR_MEMORY);
 
   size_t local_work_size = 128;
   size_t global_work_size = ((ncells+local_work_size - 1) /local_work_size) * local_work_size;
    //gpu_timers[STATE_TIMER_WRITE] += (long)(cpu_timer_stop(tstart_cpu_part)*1.0e9);
    //cpu_timer_start(&tstart_cpu_part);

   //size_t global_work_size = MAX(MAX(mesh->nxface, mesh->nyface), ncells);

#ifdef HAVE_MPI
   if (mesh->numpe > 1) {
      ezcl_set_kernel_arg(kernel_apply_boundary_conditions_local, 0, sizeof(cl_int), &ncells);
      ezcl_set_kernel_arg(kernel_apply_boundary_conditions_local, 1, sizeof(cl_mem), &dev_celltype);
      ezcl_set_kernel_arg(kernel_apply_boundary_conditions_local, 2, sizeof(cl_mem), &dev_nlft);
      ezcl_set_kernel_arg(kernel_apply_boundary_conditions_local, 3, sizeof(cl_mem), &dev_nrht);
      ezcl_set_kernel_arg(kernel_apply_boundary_conditions_local, 4, sizeof(cl_mem), &dev_ntop);
      ezcl_set_kernel_arg(kernel_apply_boundary_conditions_local, 5, sizeof(cl_mem), &dev_nbot);
      ezcl_set_kernel_arg(kernel_apply_boundary_conditions_local, 6, sizeof(cl_mem), &dev_H);
      ezcl_set_kernel_arg(kernel_apply_boundary_conditions_local, 7, sizeof(cl_mem), &dev_U);
      ezcl_set_kernel_arg(kernel_apply_boundary_conditions_local, 8, sizeof(cl_mem), &dev_V);
      ezcl_enqueue_ndrange_kernel(command_queue, kernel_apply_boundary_conditions_local,   1, NULL, &global_work_size, &local_work_size, NULL);
    
        /*
        __kernel void copy_state_data_cl(
                         const int    isize,         // 0
                __global      state_t *H,            // 1
                __global      state_t *U,            // 2
                __global      state_t *V,            // 3
                __global      state_t *H_new,        // 4
                __global      state_t *U_new,        // 5
                __global      state_t *V_new)        // 6
        */

      ezcl_set_kernel_arg(kernel_copy_state_data, 0, sizeof(cl_int), (void *)&ncells);
      ezcl_set_kernel_arg(kernel_copy_state_data, 1, sizeof(cl_mem), (void *)&dev_H);
      ezcl_set_kernel_arg(kernel_copy_state_data, 2, sizeof(cl_mem), (void *)&dev_U);
      ezcl_set_kernel_arg(kernel_copy_state_data, 3, sizeof(cl_mem), (void *)&dev_V);
      ezcl_set_kernel_arg(kernel_copy_state_data, 4, sizeof(cl_mem), (void *)&dev_H_new);
      ezcl_set_kernel_arg(kernel_copy_state_data, 5, sizeof(cl_mem), (void *)&dev_U_new);
      ezcl_set_kernel_arg(kernel_copy_state_data, 6, sizeof(cl_mem), (void *)&dev_V_new);

      //ezcl_enqueue_ndrange_kernel(command_queue, kernel_copy_state_data,   1, NULL, &global_work_size, &local_work_size, &copy_state_data_event);
      ezcl_enqueue_ndrange_kernel(command_queue, kernel_copy_state_data,   1, NULL, &global_work_size, &local_work_size, NULL);

      dev_H = (cl_mem)gpu_state_memory.memory_replace(dev_H, dev_H_new);
      dev_U = (cl_mem)gpu_state_memory.memory_replace(dev_U, dev_U_new);
      dev_V = (cl_mem)gpu_state_memory.memory_replace(dev_V, dev_V_new);

      L7_Dev_Update(dev_H, L7_STATE_T, mesh->cell_handle);
      L7_Dev_Update(dev_U, L7_STATE_T, mesh->cell_handle);
      L7_Dev_Update(dev_V, L7_STATE_T, mesh->cell_handle);

      dev_H_new = (cl_mem)gpu_state_memory.memory_malloc(ncells_ghost, sizeof(cl_state_t), const_cast<char *>("dev_H_new"), DEVICE_REGULAR_MEMORY);
      dev_U_new = (cl_mem)gpu_state_memory.memory_malloc(ncells_ghost, sizeof(cl_state_t), const_cast<char *>("dev_U_new"), DEVICE_REGULAR_MEMORY);
      dev_V_new = (cl_mem)gpu_state_memory.memory_malloc(ncells_ghost, sizeof(cl_state_t), const_cast<char *>("dev_V_new"), DEVICE_REGULAR_MEMORY);

      ezcl_set_kernel_arg(kernel_apply_boundary_conditions_ghost, 0, sizeof(cl_int), &ncells);
      ezcl_set_kernel_arg(kernel_apply_boundary_conditions_ghost, 1, sizeof(cl_mem), &dev_celltype);
      ezcl_set_kernel_arg(kernel_apply_boundary_conditions_ghost, 2, sizeof(cl_mem), &dev_nlft);
      ezcl_set_kernel_arg(kernel_apply_boundary_conditions_ghost, 3, sizeof(cl_mem), &dev_nrht);
      ezcl_set_kernel_arg(kernel_apply_boundary_conditions_ghost, 4, sizeof(cl_mem), &dev_ntop);
      ezcl_set_kernel_arg(kernel_apply_boundary_conditions_ghost, 5, sizeof(cl_mem), &dev_nbot);
      ezcl_set_kernel_arg(kernel_apply_boundary_conditions_ghost, 6, sizeof(cl_mem), &dev_H);
      ezcl_set_kernel_arg(kernel_apply_boundary_conditions_ghost, 7, sizeof(cl_mem), &dev_U);
      ezcl_set_kernel_arg(kernel_apply_boundary_conditions_ghost, 8, sizeof(cl_mem), &dev_V);
      ezcl_enqueue_ndrange_kernel(command_queue, kernel_apply_boundary_conditions_ghost,   1, NULL, &global_work_size, &local_work_size, NULL);
   } else {
      ezcl_set_kernel_arg(kernel_apply_boundary_conditions, 0, sizeof(cl_int), &ncells);
      ezcl_set_kernel_arg(kernel_apply_boundary_conditions, 1, sizeof(cl_mem), &dev_celltype);
      ezcl_set_kernel_arg(kernel_apply_boundary_conditions, 2, sizeof(cl_mem), &dev_nlft);
      ezcl_set_kernel_arg(kernel_apply_boundary_conditions, 3, sizeof(cl_mem), &dev_nrht);
      ezcl_set_kernel_arg(kernel_apply_boundary_conditions, 4, sizeof(cl_mem), &dev_ntop);
      ezcl_set_kernel_arg(kernel_apply_boundary_conditions, 5, sizeof(cl_mem), &dev_nbot);
      ezcl_set_kernel_arg(kernel_apply_boundary_conditions, 6, sizeof(cl_mem), &dev_H);
      ezcl_set_kernel_arg(kernel_apply_boundary_conditions, 7, sizeof(cl_mem), &dev_U);
      ezcl_set_kernel_arg(kernel_apply_boundary_conditions, 8, sizeof(cl_mem), &dev_V);
      ezcl_enqueue_ndrange_kernel(command_queue, kernel_apply_boundary_conditions,   1, NULL, &global_work_size, &local_work_size, NULL);
   }
#else
   ezcl_set_kernel_arg(kernel_apply_boundary_conditions, 0, sizeof(cl_int), &ncells);
   ezcl_set_kernel_arg(kernel_apply_boundary_conditions, 1, sizeof(cl_mem), (void *)&dev_celltype);
   ezcl_set_kernel_arg(kernel_apply_boundary_conditions, 2, sizeof(cl_mem), (void *)&dev_nlft);
   ezcl_set_kernel_arg(kernel_apply_boundary_conditions, 3, sizeof(cl_mem), (void *)&dev_nrht);
   ezcl_set_kernel_arg(kernel_apply_boundary_conditions, 4, sizeof(cl_mem), (void *)&dev_ntop);
   ezcl_set_kernel_arg(kernel_apply_boundary_conditions, 5, sizeof(cl_mem), (void *)&dev_nbot);
   ezcl_set_kernel_arg(kernel_apply_boundary_conditions, 6, sizeof(cl_mem), (void *)&dev_H);
   ezcl_set_kernel_arg(kernel_apply_boundary_conditions, 7, sizeof(cl_mem), (void *)&dev_U);
   ezcl_set_kernel_arg(kernel_apply_boundary_conditions, 8, sizeof(cl_mem), (void *)&dev_V);

    //gpu_timers[STATE_TIMER_WRITE] += (long)(cpu_timer_stop(tstart_cpu_part)*1.0e9);
    //cpu_timer_start(&tstart_cpu_part);

   ezcl_enqueue_ndrange_kernel(command_queue, kernel_apply_boundary_conditions,   1, NULL, &global_work_size, &local_work_size, NULL);
#endif

   size_t mem_request;
   mem_request = ncells;

   gpu_faces_setup_phantom(mem_request);

   // Get total number of elements for 1D conversion from pstate
   int total_size = 0;
   for (int lev=0; lev < levmx+1; lev++) {
      total_size +=  mesh->lev_jregsize[lev] * mesh->lev_iregsize[lev];
   //    printf("%d\n", total_size);
   }
   size_t malloc_size = total_size;

   // Allocate reg_lev arrays with that total size
   gpu_reggrid_setup(malloc_size);

   // Copy regular mesh arrays
   ezcl_enqueue_write_buffer(command_queue, dev_lev_jregmin, CL_TRUE, 0, 3*sizeof(cl_int), &mesh->lev_jregmin[0], NULL);
   ezcl_enqueue_write_buffer(command_queue, dev_lev_iregmin, CL_TRUE, 0, 3*sizeof(cl_int), &mesh->lev_iregmin[0], NULL);
   ezcl_enqueue_write_buffer(command_queue, dev_lev_jregsize, CL_TRUE, 0, 3*sizeof(cl_int), &mesh->lev_jregsize[0], NULL);
   ezcl_enqueue_write_buffer(command_queue, dev_lev_iregsize, CL_TRUE, 0, 3*sizeof(cl_int), &mesh->lev_iregsize[0], NULL);
   //ezcl_enqueue_write_buffer(command_queue, dev_j, CL_TRUE, 0, ncells*sizeof(cl_int), &mesh->j[0], NULL);
   //ezcl_enqueue_write_buffer(command_queue, dev_i, CL_TRUE, 0, ncells*sizeof(cl_int), &mesh->i[0], NULL);
   //ezcl_enqueue_write_buffer(command_queue, dev_level, CL_TRUE, 0, ncells*sizeof(cl_int), &mesh->level[0], NULL);
   //ezcl_enqueue_write_buffer(command_queue, dev_levdx, CL_TRUE, 0, 3*sizeof(cl_int), &mesh->lev_deltax[0], NULL);
   //ezcl_enqueue_write_buffer(command_queue, dev_levdy, CL_TRUE, 0, 3*sizeof(cl_int), &mesh->lev_deltay[0], NULL);

   vector<int>reg_start(levmx+1);

   // dev_H_reg_lev = ezcl_malloc(NULL, const_cast<char *>("dev_H_reg_lev"), &malloc_size, sizeof(cl_real_t), CL_MEM_READ_WRITE, 0);
   // dev_U_reg_lev = ezcl_malloc(NULL, const_cast<char *>("dev_U_reg_lev"), &malloc_size, sizeof(cl_real_t), CL_MEM_READ_WRITE, 0);
   // dev_V_reg_lev = ezcl_malloc(NULL, const_cast<char *>("dev_V_reg_lev"), &malloc_size, sizeof(cl_real_t), CL_MEM_READ_WRITE, 0);
   vector<real_t>H_reg(total_size);
   vector<real_t>U_reg(total_size);
   vector<real_t>V_reg(total_size);

   //printf("Malloc %d\n", malloc_size);
   // Get values for reg_start and also copy pstate to reg_levs
   reg_start[0] = 0;
   int startCnt = 0;
   for (int lev=0; lev < levmx+1; lev++) {
      int llsize = mesh->lev_jregsize[lev] * mesh->lev_iregsize[lev];
      //printf("startCnt %d llsize %d\n", startCnt, llsize);
      for (int copy = 0; copy < llsize; copy++) {
        H_reg[startCnt+copy] = mesh->meshes[lev].pstate[0][copy/mesh->lev_iregsize[lev]][copy%mesh->lev_iregsize[lev]];
        U_reg[startCnt+copy] = mesh->meshes[lev].pstate[1][copy/mesh->lev_iregsize[lev]][copy%mesh->lev_iregsize[lev]];
        V_reg[startCnt+copy] = mesh->meshes[lev].pstate[2][copy/mesh->lev_iregsize[lev]][copy%mesh->lev_iregsize[lev]];
      }
      startCnt += llsize;
      if (lev < levmx) reg_start[lev+1] = startCnt;
   }


      ezcl_enqueue_write_buffer(command_queue, dev_H_reg_lev, CL_TRUE, 0*sizeof(cl_real_t), total_size*sizeof(cl_real_t), &H_reg[0], NULL);
      ezcl_enqueue_write_buffer(command_queue, dev_U_reg_lev, CL_TRUE, 0*sizeof(cl_real_t), total_size*sizeof(cl_real_t), &U_reg[0], NULL);
      ezcl_enqueue_write_buffer(command_queue, dev_V_reg_lev, CL_TRUE, 0*sizeof(cl_real_t), total_size*sizeof(cl_real_t), &V_reg[0], NULL);

   // Copy up reg_start
   ezcl_enqueue_write_buffer(command_queue, dev_reg_start, CL_TRUE, 0, 3*sizeof(cl_int), &reg_start[0], NULL);

   cl_event calc_finite_difference_regular_cells_comps_event, calc_finite_difference_regular_cells_fill_event, calc_finite_difference_in_place_fixup_event;

   //for (int ic = 0; ic < ncells; ic++) {
       //int ll = mesh->level[ic];
       //int jj = mesh->j[ic] - mesh->lev_jregmin[ll];
       //int ii = mesh->i[ic] - mesh->lev_iregmin[ll];
       //printf("%d) %d %d %d\n", ic, ll, jj, ii);
       //printf("%f ", V_reg[reg_start[ll]+jj*mesh->lev_iregsize[ll]+ii]);
       //printf("%f ", V_reg[reg_start[ll]+jj*mesh->lev_iregsize[ll]+ii+1]);
       //printf("%f ", V_reg[reg_start[ll]+jj*mesh->lev_iregsize[ll]+ii-1]);
       //printf("%f ", V_reg[reg_start[ll]+(jj+1)*mesh->lev_iregsize[ll]+ii]);
       //printf("%f\n", V_reg[reg_start[ll]+(jj-1)*mesh->lev_iregsize[ll]+ii]);
        //printf("%d) ll %d jj %d ii %d startIdx %d\n", i, mesh->level[i], mesh->j[i] - mesh->lev_jregmin[mesh->level[i]], mesh->i[i] - mesh->lev_iregmin[mesh->level[i]], reg_start[mesh->level[i]]);
    //printf("%d) %f\n", i, mesh->meshes[ll].pstate[0][jj][ii]);
   //}

   struct timespec tstart_cpu;
   cpu_timer_start(&tstart_cpu);

   /*
__kernel void calc_finite_difference_regular_cells(
                        const int       ncells,                     // 0  Total number of cells.
                        const real_t    deltaT,                     // 1 Size of time step
            __global    const real_t    *lev_dx,                    // 2
            __global    const real_t    *lev_dy,                    // 3
            __global    const state_t   *Hxfluxplus,                // 4
            __global    const state_t   *Hxfluxminus,               // 5
            __global    const state_t   *Uxfluxplus,                // 6
            __global    const state_t   *Uxfluxminus,               // 7
            __global    const state_t   *Vxfluxplus,                // 8
            __global    const state_t   *Vxfluxminus,               // 9
            __global    const state_t   *Hyfluxplus,                // 10
            __global    const state_t   *Hyfluxminus,               // 11
            __global    const state_t   *Uyfluxplus,                // 12
            __global    const state_t   *Uyfluxminus,               // 13
            __global    const state_t   *Vyfluxplus,                // 14
            __global          state_t   *Vyfluxminus,               // 15
            __global    const state_t   *wplusx_H,                  // 16
            __global    const state_t   *wminusx_H,                 // 17
            __global    const state_t   *wplusx_U,                  // 18
            __global    const state_t   *wminusx_U,                 // 19
            __global    const state_t   *wplusy_H,                  // 20
            __global    const state_t   *wminusy_H,                 // 21
            __global    const state_t   *wplusy_V,                  // 22
            __global    const state_t   *wminusy_V,                 // 23
            __global    const int       *level,                     // 24
            __global    const int       *reg_start,                 // 25
            __global    const real_t    *H_reg_lev,                 // 26
            __global    const real_t    *U_reg_lev,                 // 27
            __global    const real_t    *V_reg_lev,                 // 28
            __global    const int       *j,                         // 29
            __global    const int       *i,                         // 30
            __global    const int       *lev_jregmin,               // 31
            __global    const int       *lev_iregmin,               // 32
            __global    const int       *lev_jregsize,              // 33
            __global    const int       *lev_iregsize)              // 34
{
    */


   ezcl_set_kernel_arg(kernel_calc_finite_difference_regular_cells_comps, 0, sizeof(cl_int), (void *)&ncells); 
   ezcl_set_kernel_arg(kernel_calc_finite_difference_regular_cells_comps, 1, sizeof(cl_real_t), (void *)&deltaT_local); 
   ezcl_set_kernel_arg(kernel_calc_finite_difference_regular_cells_comps, 2, sizeof(cl_mem), (void *)&dev_levdx); 
   ezcl_set_kernel_arg(kernel_calc_finite_difference_regular_cells_comps, 3, sizeof(cl_mem), (void *)&dev_levdy); 
   ezcl_set_kernel_arg(kernel_calc_finite_difference_regular_cells_comps, 4, sizeof(cl_mem), (void *)&dev_Hxfluxplus); 
   ezcl_set_kernel_arg(kernel_calc_finite_difference_regular_cells_comps, 5, sizeof(cl_mem), (void *)&dev_Hxfluxminus); 
   ezcl_set_kernel_arg(kernel_calc_finite_difference_regular_cells_comps, 6, sizeof(cl_mem), (void *)&dev_Uxfluxplus); 
   ezcl_set_kernel_arg(kernel_calc_finite_difference_regular_cells_comps, 7, sizeof(cl_mem), (void *)&dev_Uxfluxminus); 
   ezcl_set_kernel_arg(kernel_calc_finite_difference_regular_cells_comps, 8, sizeof(cl_mem), (void *)&dev_Vxfluxplus); 
   ezcl_set_kernel_arg(kernel_calc_finite_difference_regular_cells_comps, 9, sizeof(cl_mem), (void *)&dev_Vxfluxminus); 
   ezcl_set_kernel_arg(kernel_calc_finite_difference_regular_cells_comps, 10, sizeof(cl_mem), (void *)&dev_Hyfluxplus); 
   ezcl_set_kernel_arg(kernel_calc_finite_difference_regular_cells_comps, 11, sizeof(cl_mem), (void *)&dev_Hyfluxminus); 
   ezcl_set_kernel_arg(kernel_calc_finite_difference_regular_cells_comps, 12, sizeof(cl_mem), (void *)&dev_Uyfluxplus); 
   ezcl_set_kernel_arg(kernel_calc_finite_difference_regular_cells_comps, 13, sizeof(cl_mem), (void *)&dev_Uyfluxminus); 
   ezcl_set_kernel_arg(kernel_calc_finite_difference_regular_cells_comps, 14, sizeof(cl_mem), (void *)&dev_Vyfluxplus); 
   ezcl_set_kernel_arg(kernel_calc_finite_difference_regular_cells_comps, 15, sizeof(cl_mem), (void *)&dev_Vyfluxminus); 
   ezcl_set_kernel_arg(kernel_calc_finite_difference_regular_cells_comps, 16, sizeof(cl_mem), (void *)&dev_Wplusx_H); 
   ezcl_set_kernel_arg(kernel_calc_finite_difference_regular_cells_comps, 17, sizeof(cl_mem), (void *)&dev_Wminusx_H); 
   ezcl_set_kernel_arg(kernel_calc_finite_difference_regular_cells_comps, 18, sizeof(cl_mem), (void *)&dev_Wplusx_U); 
   ezcl_set_kernel_arg(kernel_calc_finite_difference_regular_cells_comps, 19, sizeof(cl_mem), (void *)&dev_Wminusx_U); 
   ezcl_set_kernel_arg(kernel_calc_finite_difference_regular_cells_comps, 20, sizeof(cl_mem), (void *)&dev_Wplusy_H); 
   ezcl_set_kernel_arg(kernel_calc_finite_difference_regular_cells_comps, 21, sizeof(cl_mem), (void *)&dev_Wminusy_H); 
   ezcl_set_kernel_arg(kernel_calc_finite_difference_regular_cells_comps, 22, sizeof(cl_mem), (void *)&dev_Wplusy_V); 
   ezcl_set_kernel_arg(kernel_calc_finite_difference_regular_cells_comps, 23, sizeof(cl_mem), (void *)&dev_Wminusy_V); 
   ezcl_set_kernel_arg(kernel_calc_finite_difference_regular_cells_comps, 24, sizeof(cl_mem), (void *)&dev_level); 
   ezcl_set_kernel_arg(kernel_calc_finite_difference_regular_cells_comps, 25, sizeof(cl_mem), (void *)&dev_reg_start); 
   ezcl_set_kernel_arg(kernel_calc_finite_difference_regular_cells_comps, 26, sizeof(cl_mem), (void *)&dev_H_reg_lev); 
   ezcl_set_kernel_arg(kernel_calc_finite_difference_regular_cells_comps, 27, sizeof(cl_mem), (void *)&dev_U_reg_lev); 
   ezcl_set_kernel_arg(kernel_calc_finite_difference_regular_cells_comps, 28, sizeof(cl_mem), (void *)&dev_V_reg_lev); 
   ezcl_set_kernel_arg(kernel_calc_finite_difference_regular_cells_comps, 29, sizeof(cl_mem), (void *)&dev_j); 
   ezcl_set_kernel_arg(kernel_calc_finite_difference_regular_cells_comps, 30, sizeof(cl_mem), (void *)&dev_i); 
   ezcl_set_kernel_arg(kernel_calc_finite_difference_regular_cells_comps, 31, sizeof(cl_mem), (void *)&dev_lev_jregmin); 
   ezcl_set_kernel_arg(kernel_calc_finite_difference_regular_cells_comps, 32, sizeof(cl_mem), (void *)&dev_lev_iregmin); 
   ezcl_set_kernel_arg(kernel_calc_finite_difference_regular_cells_comps, 33, sizeof(cl_mem), (void *)&dev_lev_jregsize); 
   ezcl_set_kernel_arg(kernel_calc_finite_difference_regular_cells_comps, 34, sizeof(cl_mem), (void *)&dev_lev_iregsize); 

   ezcl_enqueue_ndrange_kernel(command_queue, kernel_calc_finite_difference_regular_cells_comps, 1, NULL, &global_work_size, &local_work_size, &calc_finite_difference_regular_cells_comps_event);

   ezcl_wait_for_events(1, &calc_finite_difference_regular_cells_comps_event);
   ezcl_event_release(calc_finite_difference_regular_cells_comps_event);
/*
   vector<real_t>H_loc(ncells);
   ezcl_enqueue_read_buffer(command_queue, dev_Hyfluxminus,     CL_TRUE, 0, ncells*sizeof(cl_real_t), &H_loc[0], NULL);
   vector<real_t>U_loc(ncells);
   ezcl_enqueue_read_buffer(command_queue, dev_Uyfluxminus,     CL_TRUE, 0, ncells*sizeof(cl_real_t), &U_loc[0], NULL);
   vector<real_t>V_loc(ncells);
   ezcl_enqueue_read_buffer(command_queue, dev_Vyfluxminus,     CL_TRUE, 0, ncells*sizeof(cl_real_t), &V_loc[0], NULL);
   printf("\n");
   for (int jello = 0; jello < ncells; jello++) { printf("%d) %f | %f | %f\n", jello, H_loc[jello], U_loc[jello], V_loc[jello]); }
*/

/*
   vector<real_t>H_loc(total_size);
   ezcl_enqueue_read_buffer(command_queue, dev_H_reg_lev,     CL_TRUE, 0, total_size*sizeof(cl_real_t), &H_loc[0], NULL);
   printf("\n");
   for (int jello = 0; jello < total_size; jello++) { printf("%d) %f\n", jello, H_loc[jello]); }
*/
   /*
__kernel void calc_finite_difference_in_place_fixup_cl(
                        const int        nxfixup,                   // 0
                        const int        nyfixup,                   // 1
            __global    const int       *xrecvCIdx,                 // 2
            __global    const int       *xplusCell2Idx,             // 3
            __global    const int       *xminusCell2Idx,            // 4
            __global    const int       *xsendIdx1,                 // 5
            __global    const int       *xsendIdx2,                 // 6
            __global    const int       *yrecvCIdx,                 // 7
            __global    const int       *yplusCell2Idx,             // 8
            __global    const int       *yminusCell2Idx,            // 9
            __global    const int       *ysendIdx1,                 // 10
            __global    const int       *ysendIdx2,                 // 11
            __global    const int       *map_xface2cell_lower,      // 12
            __global    const int       *map_xface2cell_upper,      // 13
            __global    const int       *map_yface2cell_lower,      // 14 
            __global    const int       *map_yface2cell_upper,      // 15 
            __global          state_t   *Hxfluxplus,                // 16
            __global          state_t   *Hxfluxminus,               // 17
            __global          state_t   *Uxfluxplus,                // 18
            __global          state_t   *Uxfluxminus,               // 19
            __global          state_t   *Vxfluxplus,                // 20
            __global          state_t   *Vxfluxminus,               // 21
            __global          state_t   *Hyfluxplus,                // 22
            __global          state_t   *Hyfluxminus,               // 23
            __global          state_t   *Uyfluxplus,                // 24
            __global          state_t   *Uyfluxminus,               // 25
            __global          state_t   *Vyfluxplus,                // 26
            __global          state_t   *Vyfluxminus,               // 27
            __global          state_t   *wplusx_H,                  // 28
            __global          state_t   *wminusx_H,                 // 29
            __global          state_t   *wplusx_U,                  // 30
            __global          state_t   *wminusx_U,                 // 31
            __global          state_t   *wplusy_H,                  // 32
            __global          state_t   *wminusy_H,                 // 33
            __global          state_t   *wplusy_V,                  // 34
            __global          state_t   *wminusy_V) {               // 35
    */

    ezcl_set_kernel_arg(kernel_calc_finite_difference_in_place_fixup, 0, sizeof(cl_int), (void *)&nxfixup); 
    ezcl_set_kernel_arg(kernel_calc_finite_difference_in_place_fixup, 1, sizeof(cl_int), (void *)&nyfixup); 
    ezcl_set_kernel_arg(kernel_calc_finite_difference_in_place_fixup, 2, sizeof(cl_mem), (void *)&dev_xrecvCIdx);
    ezcl_set_kernel_arg(kernel_calc_finite_difference_in_place_fixup, 3, sizeof(cl_mem), (void *)&dev_xplusCell2Idx);
    ezcl_set_kernel_arg(kernel_calc_finite_difference_in_place_fixup, 4, sizeof(cl_mem), (void *)&dev_xminusCell2Idx);
    ezcl_set_kernel_arg(kernel_calc_finite_difference_in_place_fixup, 5, sizeof(cl_mem), (void *)&dev_xsendIdx1);
    ezcl_set_kernel_arg(kernel_calc_finite_difference_in_place_fixup, 6, sizeof(cl_mem), (void *)&dev_xsendIdx2);
    ezcl_set_kernel_arg(kernel_calc_finite_difference_in_place_fixup, 7, sizeof(cl_mem), (void *)&dev_yrecvCIdx);
    ezcl_set_kernel_arg(kernel_calc_finite_difference_in_place_fixup, 8, sizeof(cl_mem), (void *)&dev_yplusCell2Idx);
    ezcl_set_kernel_arg(kernel_calc_finite_difference_in_place_fixup, 9, sizeof(cl_mem), (void *)&dev_yminusCell2Idx);
    ezcl_set_kernel_arg(kernel_calc_finite_difference_in_place_fixup, 10, sizeof(cl_mem), (void *)&dev_ysendIdx1);
    ezcl_set_kernel_arg(kernel_calc_finite_difference_in_place_fixup, 11, sizeof(cl_mem), (void *)&dev_ysendIdx2);
    ezcl_set_kernel_arg(kernel_calc_finite_difference_in_place_fixup, 12, sizeof(cl_mem), (void *)&dev_map_xface2cell_lower); 
    ezcl_set_kernel_arg(kernel_calc_finite_difference_in_place_fixup, 13, sizeof(cl_mem), (void *)&dev_map_xface2cell_upper); 
    ezcl_set_kernel_arg(kernel_calc_finite_difference_in_place_fixup, 14, sizeof(cl_mem), (void *)&dev_map_yface2cell_lower); 
    ezcl_set_kernel_arg(kernel_calc_finite_difference_in_place_fixup, 15, sizeof(cl_mem), (void *)&dev_map_yface2cell_upper); 
    ezcl_set_kernel_arg(kernel_calc_finite_difference_in_place_fixup, 16, sizeof(cl_mem), (void *)&dev_Hxfluxplus); 
    ezcl_set_kernel_arg(kernel_calc_finite_difference_in_place_fixup, 17, sizeof(cl_mem), (void *)&dev_Hxfluxminus); 
    ezcl_set_kernel_arg(kernel_calc_finite_difference_in_place_fixup, 18, sizeof(cl_mem), (void *)&dev_Uxfluxplus); 
    ezcl_set_kernel_arg(kernel_calc_finite_difference_in_place_fixup, 19, sizeof(cl_mem), (void *)&dev_Uxfluxminus); 
    ezcl_set_kernel_arg(kernel_calc_finite_difference_in_place_fixup, 20, sizeof(cl_mem), (void *)&dev_Vxfluxplus); 
    ezcl_set_kernel_arg(kernel_calc_finite_difference_in_place_fixup, 21, sizeof(cl_mem), (void *)&dev_Vxfluxminus); 
    ezcl_set_kernel_arg(kernel_calc_finite_difference_in_place_fixup, 22, sizeof(cl_mem), (void *)&dev_Hyfluxplus); 
    ezcl_set_kernel_arg(kernel_calc_finite_difference_in_place_fixup, 23, sizeof(cl_mem), (void *)&dev_Hyfluxminus); 
    ezcl_set_kernel_arg(kernel_calc_finite_difference_in_place_fixup, 24, sizeof(cl_mem), (void *)&dev_Uyfluxplus); 
    ezcl_set_kernel_arg(kernel_calc_finite_difference_in_place_fixup, 25, sizeof(cl_mem), (void *)&dev_Uyfluxminus); 
    ezcl_set_kernel_arg(kernel_calc_finite_difference_in_place_fixup, 26, sizeof(cl_mem), (void *)&dev_Vyfluxplus); 
    ezcl_set_kernel_arg(kernel_calc_finite_difference_in_place_fixup, 27, sizeof(cl_mem), (void *)&dev_Vyfluxminus); 
    ezcl_set_kernel_arg(kernel_calc_finite_difference_in_place_fixup, 28, sizeof(cl_mem), (void *)&dev_Wplusx_H); 
    ezcl_set_kernel_arg(kernel_calc_finite_difference_in_place_fixup, 29, sizeof(cl_mem), (void *)&dev_Wminusx_H); 
    ezcl_set_kernel_arg(kernel_calc_finite_difference_in_place_fixup, 30, sizeof(cl_mem), (void *)&dev_Wplusx_U); 
    ezcl_set_kernel_arg(kernel_calc_finite_difference_in_place_fixup, 31, sizeof(cl_mem), (void *)&dev_Wminusx_U); 
    ezcl_set_kernel_arg(kernel_calc_finite_difference_in_place_fixup, 32, sizeof(cl_mem), (void *)&dev_Wplusy_H); 
    ezcl_set_kernel_arg(kernel_calc_finite_difference_in_place_fixup, 33, sizeof(cl_mem), (void *)&dev_Wminusy_H); 
    ezcl_set_kernel_arg(kernel_calc_finite_difference_in_place_fixup, 34, sizeof(cl_mem), (void *)&dev_Wplusy_V); 
    ezcl_set_kernel_arg(kernel_calc_finite_difference_in_place_fixup, 35, sizeof(cl_mem), (void *)&dev_Wminusy_V); 

    //gpu_timers[STATE_TIMER_WRITE] += (long)(cpu_timer_stop(tstart_cpu_part)*1.0e9);
    //cpu_timer_start(&tstart_cpu_part);

    ezcl_enqueue_ndrange_kernel(command_queue, kernel_calc_finite_difference_in_place_fixup, 1, NULL, &global_work_size, &local_work_size, &calc_finite_difference_in_place_fixup_event);

    ezcl_wait_for_events(1, &calc_finite_difference_in_place_fixup_event);
    ezcl_event_release(calc_finite_difference_in_place_fixup_event);

    /*
   vector<real_t>H_loc(ncells);
   ezcl_enqueue_read_buffer(command_queue, dev_H,     CL_TRUE, 0, ncells*sizeof(cl_real_t), &H_loc[0], NULL);
   vector<real_t>U_loc(ncells);
   ezcl_enqueue_read_buffer(command_queue, dev_U,     CL_TRUE, 0, ncells*sizeof(cl_real_t), &U_loc[0], NULL);
   vector<real_t>V_loc(ncells);
   ezcl_enqueue_read_buffer(command_queue, dev_V,     CL_TRUE, 0, ncells*sizeof(cl_real_t), &V_loc[0], NULL);
   printf("\n");
   for (int jello = 0; jello < ncells; jello++) { printf("%d) %f | %f | %f\n", jello, H_loc[jello], U_loc[jello], V_loc[jello]); }
   */
/*
__kernel void calc_finite_difference_regular_cells_fill_cl(
                              int       ncells,                     // 0  Total number of cells.
                              real_t    deltaT,                     // 1 Size of time step
            __global    const real_t    *lev_dx,                    // 2
            __global    const real_t    *lev_dy,                    // 3
            __global          state_t   *Hxfluxplus,                // 4
            __global          state_t   *Hxfluxminus,               // 5
            __global          state_t   *Uxfluxplus,                // 6
            __global          state_t   *Uxfluxminus,               // 7
            __global          state_t   *Vxfluxplus,                // 8
            __global          state_t   *Vxfluxminus,               // 9
            __global          state_t   *Hyfluxplus,                // 10
            __global          state_t   *Hyfluxminus,               // 11
            __global          state_t   *Uyfluxplus,                // 12
            __global          state_t   *Uyfluxminus,               // 13
            __global          state_t   *Vyfluxplus,                // 14
            __global          state_t   *Vyfluxminus,               // 15
            __global          state_t   *wplusx_H,                  // 16
            __global          state_t   *wminusx_H,                 // 17
            __global          state_t   *wplusx_U,                  // 18
            __global          state_t   *wminusx_U,                 // 19
            __global          state_t   *wplusy_H,                  // 20
            __global          state_t   *wminusy_H,                 // 21
            __global          state_t   *wplusy_V,                  // 22
            __global          state_t   *wminusy_V,                 // 23
            __global    const int       *level,                     // 24
            __global    const int       *reg_start,                 // 25
            __global    const real_t    *H_reg_lev,                 // 26
            __global    const real_t    *U_reg_lev,                 // 27
            __global    const real_t    *V_reg_lev,                 // 28
            __global    const int       *j,                         // 29
            __global    const int       *i,                         // 30
            __global    const int       *lev_jregmin,               // 31
            __global    const int       *lev_iregmin,               // 32
            __global    const int       *lev_jregsize,              // 33
            __global    const int       *lev_iregsize,              // 34
            __global    const real_t    *H_state_new,               // 35
            __global    const real_t    *U_state_new,               // 36
            __global    const real_t    *V_state_new)               // 37
*/

   ezcl_set_kernel_arg(kernel_calc_finite_difference_regular_cells_fill, 0, sizeof(cl_int), (void *)&ncells); 
   ezcl_set_kernel_arg(kernel_calc_finite_difference_regular_cells_fill, 1, sizeof(cl_real_t), (void *)&deltaT_local); 
   ezcl_set_kernel_arg(kernel_calc_finite_difference_regular_cells_fill, 2, sizeof(cl_mem), (void *)&dev_levdx); 
   ezcl_set_kernel_arg(kernel_calc_finite_difference_regular_cells_fill, 3, sizeof(cl_mem), (void *)&dev_levdy); 
   ezcl_set_kernel_arg(kernel_calc_finite_difference_regular_cells_fill, 4, sizeof(cl_mem), (void *)&dev_Hxfluxplus); 
   ezcl_set_kernel_arg(kernel_calc_finite_difference_regular_cells_fill, 5, sizeof(cl_mem), (void *)&dev_Hxfluxminus); 
   ezcl_set_kernel_arg(kernel_calc_finite_difference_regular_cells_fill, 6, sizeof(cl_mem), (void *)&dev_Uxfluxplus); 
   ezcl_set_kernel_arg(kernel_calc_finite_difference_regular_cells_fill, 7, sizeof(cl_mem), (void *)&dev_Uxfluxminus); 
   ezcl_set_kernel_arg(kernel_calc_finite_difference_regular_cells_fill, 8, sizeof(cl_mem), (void *)&dev_Vxfluxplus); 
   ezcl_set_kernel_arg(kernel_calc_finite_difference_regular_cells_fill, 9, sizeof(cl_mem), (void *)&dev_Vxfluxminus); 
   ezcl_set_kernel_arg(kernel_calc_finite_difference_regular_cells_fill, 10, sizeof(cl_mem), (void *)&dev_Hyfluxplus); 
   ezcl_set_kernel_arg(kernel_calc_finite_difference_regular_cells_fill, 11, sizeof(cl_mem), (void *)&dev_Hyfluxminus); 
   ezcl_set_kernel_arg(kernel_calc_finite_difference_regular_cells_fill, 12, sizeof(cl_mem), (void *)&dev_Uyfluxplus); 
   ezcl_set_kernel_arg(kernel_calc_finite_difference_regular_cells_fill, 13, sizeof(cl_mem), (void *)&dev_Uyfluxminus); 
   ezcl_set_kernel_arg(kernel_calc_finite_difference_regular_cells_fill, 14, sizeof(cl_mem), (void *)&dev_Vyfluxplus); 
   ezcl_set_kernel_arg(kernel_calc_finite_difference_regular_cells_fill, 15, sizeof(cl_mem), (void *)&dev_Vyfluxminus); 
   ezcl_set_kernel_arg(kernel_calc_finite_difference_regular_cells_fill, 16, sizeof(cl_mem), (void *)&dev_Wplusx_H); 
   ezcl_set_kernel_arg(kernel_calc_finite_difference_regular_cells_fill, 17, sizeof(cl_mem), (void *)&dev_Wminusx_H); 
   ezcl_set_kernel_arg(kernel_calc_finite_difference_regular_cells_fill, 18, sizeof(cl_mem), (void *)&dev_Wplusx_U); 
   ezcl_set_kernel_arg(kernel_calc_finite_difference_regular_cells_fill, 19, sizeof(cl_mem), (void *)&dev_Wminusx_U); 
   ezcl_set_kernel_arg(kernel_calc_finite_difference_regular_cells_fill, 20, sizeof(cl_mem), (void *)&dev_Wplusy_H); 
   ezcl_set_kernel_arg(kernel_calc_finite_difference_regular_cells_fill, 21, sizeof(cl_mem), (void *)&dev_Wminusy_H); 
   ezcl_set_kernel_arg(kernel_calc_finite_difference_regular_cells_fill, 22, sizeof(cl_mem), (void *)&dev_Wplusy_V); 
   ezcl_set_kernel_arg(kernel_calc_finite_difference_regular_cells_fill, 23, sizeof(cl_mem), (void *)&dev_Wminusy_V); 
   ezcl_set_kernel_arg(kernel_calc_finite_difference_regular_cells_fill, 24, sizeof(cl_mem), (void *)&dev_level); 
   ezcl_set_kernel_arg(kernel_calc_finite_difference_regular_cells_fill, 25, sizeof(cl_mem), (void *)&dev_reg_start); 
   ezcl_set_kernel_arg(kernel_calc_finite_difference_regular_cells_fill, 26, sizeof(cl_mem), (void *)&dev_H_reg_lev); 
   ezcl_set_kernel_arg(kernel_calc_finite_difference_regular_cells_fill, 27, sizeof(cl_mem), (void *)&dev_U_reg_lev); 
   ezcl_set_kernel_arg(kernel_calc_finite_difference_regular_cells_fill, 28, sizeof(cl_mem), (void *)&dev_V_reg_lev); 
   ezcl_set_kernel_arg(kernel_calc_finite_difference_regular_cells_fill, 29, sizeof(cl_mem), (void *)&dev_j); 
   ezcl_set_kernel_arg(kernel_calc_finite_difference_regular_cells_fill, 30, sizeof(cl_mem), (void *)&dev_i); 
   ezcl_set_kernel_arg(kernel_calc_finite_difference_regular_cells_fill, 31, sizeof(cl_mem), (void *)&dev_lev_jregmin); 
   ezcl_set_kernel_arg(kernel_calc_finite_difference_regular_cells_fill, 32, sizeof(cl_mem), (void *)&dev_lev_iregmin); 
   ezcl_set_kernel_arg(kernel_calc_finite_difference_regular_cells_fill, 33, sizeof(cl_mem), (void *)&dev_lev_jregsize); 
   ezcl_set_kernel_arg(kernel_calc_finite_difference_regular_cells_fill, 34, sizeof(cl_mem), (void *)&dev_lev_iregsize); 
   ezcl_set_kernel_arg(kernel_calc_finite_difference_regular_cells_fill, 35, sizeof(cl_mem), (void *)&dev_H_new); 
   ezcl_set_kernel_arg(kernel_calc_finite_difference_regular_cells_fill, 36, sizeof(cl_mem), (void *)&dev_U_new); 
   ezcl_set_kernel_arg(kernel_calc_finite_difference_regular_cells_fill, 37, sizeof(cl_mem), (void *)&dev_V_new); 

   ezcl_enqueue_ndrange_kernel(command_queue, kernel_calc_finite_difference_regular_cells_fill, 1, NULL, &global_work_size, &local_work_size, &calc_finite_difference_regular_cells_fill_event);

   ezcl_wait_for_events(1, &calc_finite_difference_regular_cells_fill_event);
   ezcl_event_release(calc_finite_difference_regular_cells_fill_event);

   gpu_timers[STATE_TIMER_FINITE_DIFFERENCE] += (long)(cpu_timer_stop(tstart_cpu)*1.0e9);

   dev_H = (cl_mem)gpu_state_memory.memory_replace(dev_H, dev_H_new);
   dev_U = (cl_mem)gpu_state_memory.memory_replace(dev_U, dev_U_new);
   dev_V = (cl_mem)gpu_state_memory.memory_replace(dev_V, dev_V_new);



   gpu_faces_delete_phantom();
   gpu_reggrid_delete();
   mesh->destroy_regular_cell_meshes(state_memory);
   mesh->gpu_wbidirmap_delete();
}

void State::gpu_calc_finite_difference_regular_cells_by_faces(double deltaT)
{

   //struct timespec tstart_cpu_part;
   //cpu_timer_start(&tstart_cpu_part);

   cl_command_queue command_queue = ezcl_get_command_queue();


   //cl_mem dev_ptr = NULL;

   size_t &ncells    = mesh->ncells;
   size_t &ncells_ghost = mesh->ncells_ghost;
   if (ncells_ghost < ncells) ncells_ghost = ncells;
   int &levmx           = mesh->levmx;
   real_t deltaT_local = deltaT;

   mesh->gpu_wbidirmap_setup();
    mesh->gpu_calc_face_list_wbidirmap_phantom(gpu_state_memory, deltaT_local);
    gpu_memory_reset_ptrs();

   mesh->generate_regular_cell_meshes(gpu_state_memory);

   cl_mem &dev_nface    = mesh->dev_nface;
   cl_mem &dev_celltype = mesh->dev_celltype;
   cl_mem &dev_nlft     = mesh->dev_nlft;
   cl_mem &dev_nrht     = mesh->dev_nrht;
   cl_mem &dev_nbot     = mesh->dev_nbot;
   cl_mem &dev_ntop     = mesh->dev_ntop;
   cl_mem &dev_level    = mesh->dev_level;
   cl_mem &dev_levdx    = mesh->dev_levdx;
   cl_mem &dev_levdy    = mesh->dev_levdy;
   cl_mem &dev_j        = mesh->dev_j;
   cl_mem &dev_i        = mesh->dev_i;
   cl_mem &dev_xface_level    = mesh->dev_xface_level;
   cl_mem &dev_xface_j    = mesh->dev_xface_j;
   cl_mem &dev_xface_i    = mesh->dev_xface_i;
   cl_mem &dev_yface_level    = mesh->dev_yface_level;
   cl_mem &dev_yface_j    = mesh->dev_yface_j;
   cl_mem &dev_yface_i    = mesh->dev_yface_i;
   cl_mem &dev_map_xface2cell_lower = mesh->dev_map_xface2cell_lower;
   cl_mem &dev_map_xface2cell_upper = mesh->dev_map_xface2cell_upper;
   cl_mem &dev_map_yface2cell_lower = mesh->dev_map_yface2cell_lower;
   cl_mem &dev_map_yface2cell_upper = mesh->dev_map_yface2cell_upper;
   cl_mem &dev_map_xcell2face_left1 = mesh->dev_map_xcell2face_left1;
   cl_mem &dev_map_xcell2face_right1 = mesh->dev_map_xcell2face_right1;
   cl_mem &dev_map_ycell2face_bot1 = mesh->dev_map_ycell2face_bot1;
   cl_mem &dev_map_ycell2face_top1 = mesh->dev_map_ycell2face_top1;
   cl_mem &dev_xrecvIdx = mesh->dev_xrecvIdx;
   cl_mem &dev_xsendIdx1 = mesh->dev_xsendIdx1;
   cl_mem &dev_xsendIdx2 = mesh->dev_xsendIdx2;
   cl_mem &dev_yrecvIdx = mesh->dev_yrecvIdx;
   cl_mem &dev_ysendIdx1 = mesh->dev_ysendIdx1;
   cl_mem &dev_ysendIdx2 = mesh->dev_ysendIdx2;
   int &nxfixup = mesh->nxfixup;
   int &nyfixup = mesh->nyfixup;

   assert(dev_nface);
   assert(dev_H);
   assert(dev_U);
   assert(dev_V);
   assert(dev_nlft);
   assert(dev_nrht);
   assert(dev_nbot);
   assert(dev_ntop);
   assert(dev_level);
   assert(dev_levdx);
   assert(dev_levdy);

   cl_mem dev_H_new = (cl_mem)gpu_state_memory.memory_malloc(ncells, sizeof(cl_state_t), const_cast<char *>("dev_H_new"), DEVICE_REGULAR_MEMORY);
   cl_mem dev_U_new = (cl_mem)gpu_state_memory.memory_malloc(ncells, sizeof(cl_state_t), const_cast<char *>("dev_U_new"), DEVICE_REGULAR_MEMORY);
   cl_mem dev_V_new = (cl_mem)gpu_state_memory.memory_malloc(ncells, sizeof(cl_state_t), const_cast<char *>("dev_V_new"), DEVICE_REGULAR_MEMORY);
 
   size_t local_work_size = 128;
   size_t global_work_size = ((ncells+local_work_size - 1) /local_work_size) * local_work_size;
    //gpu_timers[STATE_TIMER_WRITE] += (long)(cpu_timer_stop(tstart_cpu_part)*1.0e9);
    //cpu_timer_start(&tstart_cpu_part);

   //size_t global_work_size = MAX(MAX(mesh->nxface, mesh->nyface), ncells);

#ifdef HAVE_MPI
   if (mesh->numpe > 1) {
      ezcl_set_kernel_arg(kernel_apply_boundary_conditions_local, 0, sizeof(cl_int), &ncells);
      ezcl_set_kernel_arg(kernel_apply_boundary_conditions_local, 1, sizeof(cl_mem), &dev_celltype);
      ezcl_set_kernel_arg(kernel_apply_boundary_conditions_local, 2, sizeof(cl_mem), &dev_nlft);
      ezcl_set_kernel_arg(kernel_apply_boundary_conditions_local, 3, sizeof(cl_mem), &dev_nrht);
      ezcl_set_kernel_arg(kernel_apply_boundary_conditions_local, 4, sizeof(cl_mem), &dev_ntop);
      ezcl_set_kernel_arg(kernel_apply_boundary_conditions_local, 5, sizeof(cl_mem), &dev_nbot);
      ezcl_set_kernel_arg(kernel_apply_boundary_conditions_local, 6, sizeof(cl_mem), &dev_H);
      ezcl_set_kernel_arg(kernel_apply_boundary_conditions_local, 7, sizeof(cl_mem), &dev_U);
      ezcl_set_kernel_arg(kernel_apply_boundary_conditions_local, 8, sizeof(cl_mem), &dev_V);
      ezcl_enqueue_ndrange_kernel(command_queue, kernel_apply_boundary_conditions_local,   1, NULL, &global_work_size, &local_work_size, NULL);
    
        /*
        __kernel void copy_state_data_cl(
                         const int    isize,         // 0
                __global      state_t *H,            // 1
                __global      state_t *U,            // 2
                __global      state_t *V,            // 3
                __global      state_t *H_new,        // 4
                __global      state_t *U_new,        // 5
                __global      state_t *V_new)        // 6
        */

      ezcl_set_kernel_arg(kernel_copy_state_data, 0, sizeof(cl_int), (void *)&ncells);
      ezcl_set_kernel_arg(kernel_copy_state_data, 1, sizeof(cl_mem), (void *)&dev_H);
      ezcl_set_kernel_arg(kernel_copy_state_data, 2, sizeof(cl_mem), (void *)&dev_U);
      ezcl_set_kernel_arg(kernel_copy_state_data, 3, sizeof(cl_mem), (void *)&dev_V);
      ezcl_set_kernel_arg(kernel_copy_state_data, 4, sizeof(cl_mem), (void *)&dev_H_new);
      ezcl_set_kernel_arg(kernel_copy_state_data, 5, sizeof(cl_mem), (void *)&dev_U_new);
      ezcl_set_kernel_arg(kernel_copy_state_data, 6, sizeof(cl_mem), (void *)&dev_V_new);

      //ezcl_enqueue_ndrange_kernel(command_queue, kernel_copy_state_data,   1, NULL, &global_work_size, &local_work_size, &copy_state_data_event);
      ezcl_enqueue_ndrange_kernel(command_queue, kernel_copy_state_data,   1, NULL, &global_work_size, &local_work_size, NULL);

      dev_H = (cl_mem)gpu_state_memory.memory_replace(dev_H, dev_H_new);
      dev_U = (cl_mem)gpu_state_memory.memory_replace(dev_U, dev_U_new);
      dev_V = (cl_mem)gpu_state_memory.memory_replace(dev_V, dev_V_new);

      L7_Dev_Update(dev_H, L7_STATE_T, mesh->cell_handle);
      L7_Dev_Update(dev_U, L7_STATE_T, mesh->cell_handle);
      L7_Dev_Update(dev_V, L7_STATE_T, mesh->cell_handle);

      dev_H_new = (cl_mem)gpu_state_memory.memory_malloc(ncells_ghost, sizeof(cl_state_t), const_cast<char *>("dev_H_new"), DEVICE_REGULAR_MEMORY);
      dev_U_new = (cl_mem)gpu_state_memory.memory_malloc(ncells_ghost, sizeof(cl_state_t), const_cast<char *>("dev_U_new"), DEVICE_REGULAR_MEMORY);
      dev_V_new = (cl_mem)gpu_state_memory.memory_malloc(ncells_ghost, sizeof(cl_state_t), const_cast<char *>("dev_V_new"), DEVICE_REGULAR_MEMORY);

      ezcl_set_kernel_arg(kernel_apply_boundary_conditions_ghost, 0, sizeof(cl_int), &ncells);
      ezcl_set_kernel_arg(kernel_apply_boundary_conditions_ghost, 1, sizeof(cl_mem), &dev_celltype);
      ezcl_set_kernel_arg(kernel_apply_boundary_conditions_ghost, 2, sizeof(cl_mem), &dev_nlft);
      ezcl_set_kernel_arg(kernel_apply_boundary_conditions_ghost, 3, sizeof(cl_mem), &dev_nrht);
      ezcl_set_kernel_arg(kernel_apply_boundary_conditions_ghost, 4, sizeof(cl_mem), &dev_ntop);
      ezcl_set_kernel_arg(kernel_apply_boundary_conditions_ghost, 5, sizeof(cl_mem), &dev_nbot);
      ezcl_set_kernel_arg(kernel_apply_boundary_conditions_ghost, 6, sizeof(cl_mem), &dev_H);
      ezcl_set_kernel_arg(kernel_apply_boundary_conditions_ghost, 7, sizeof(cl_mem), &dev_U);
      ezcl_set_kernel_arg(kernel_apply_boundary_conditions_ghost, 8, sizeof(cl_mem), &dev_V);
      ezcl_enqueue_ndrange_kernel(command_queue, kernel_apply_boundary_conditions_ghost,   1, NULL, &global_work_size, &local_work_size, NULL);
   } else {
      ezcl_set_kernel_arg(kernel_apply_boundary_conditions, 0, sizeof(cl_int), &ncells);
      ezcl_set_kernel_arg(kernel_apply_boundary_conditions, 1, sizeof(cl_mem), &dev_celltype);
      ezcl_set_kernel_arg(kernel_apply_boundary_conditions, 2, sizeof(cl_mem), &dev_nlft);
      ezcl_set_kernel_arg(kernel_apply_boundary_conditions, 3, sizeof(cl_mem), &dev_nrht);
      ezcl_set_kernel_arg(kernel_apply_boundary_conditions, 4, sizeof(cl_mem), &dev_ntop);
      ezcl_set_kernel_arg(kernel_apply_boundary_conditions, 5, sizeof(cl_mem), &dev_nbot);
      ezcl_set_kernel_arg(kernel_apply_boundary_conditions, 6, sizeof(cl_mem), &dev_H);
      ezcl_set_kernel_arg(kernel_apply_boundary_conditions, 7, sizeof(cl_mem), &dev_U);
      ezcl_set_kernel_arg(kernel_apply_boundary_conditions, 8, sizeof(cl_mem), &dev_V);
      ezcl_enqueue_ndrange_kernel(command_queue, kernel_apply_boundary_conditions,   1, NULL, &global_work_size, &local_work_size, NULL);
   }
#else
   ezcl_set_kernel_arg(kernel_apply_boundary_conditions, 0, sizeof(cl_int), &ncells);
   ezcl_set_kernel_arg(kernel_apply_boundary_conditions, 1, sizeof(cl_mem), (void *)&dev_celltype);
   ezcl_set_kernel_arg(kernel_apply_boundary_conditions, 2, sizeof(cl_mem), (void *)&dev_nlft);
   ezcl_set_kernel_arg(kernel_apply_boundary_conditions, 3, sizeof(cl_mem), (void *)&dev_nrht);
   ezcl_set_kernel_arg(kernel_apply_boundary_conditions, 4, sizeof(cl_mem), (void *)&dev_ntop);
   ezcl_set_kernel_arg(kernel_apply_boundary_conditions, 5, sizeof(cl_mem), (void *)&dev_nbot);
   ezcl_set_kernel_arg(kernel_apply_boundary_conditions, 6, sizeof(cl_mem), (void *)&dev_H);
   ezcl_set_kernel_arg(kernel_apply_boundary_conditions, 7, sizeof(cl_mem), (void *)&dev_U);
   ezcl_set_kernel_arg(kernel_apply_boundary_conditions, 8, sizeof(cl_mem), (void *)&dev_V);

    //gpu_timers[STATE_TIMER_WRITE] += (long)(cpu_timer_stop(tstart_cpu_part)*1.0e9);
    //cpu_timer_start(&tstart_cpu_part);

   ezcl_enqueue_ndrange_kernel(command_queue, kernel_apply_boundary_conditions,   1, NULL, &global_work_size, &local_work_size, NULL);
#endif

   size_t mem_requestx, mem_requesty;

   mem_requestx = mesh->pxfaceCnt;
   mem_requesty = mesh->pyfaceCnt;
   gpu_faces_setup(mem_requestx, mem_requesty);

   //size_t nfacex, nfacey;
   //ezcl_enqueue_read_buffer(command_queue, dev_nface,     CL_TRUE, 0, sizeof(cl_int), &nfacex, NULL);
   //ezcl_enqueue_read_buffer(command_queue, dev_nface,     CL_TRUE, 1*sizeof(cl_int), sizeof(cl_int), &nfacey, NULL);
   //printf("\nMem requests %d and %d\n", nfacex, nfacey);

   // Get total number of elements for 1D conversion from pstate
   int total_size = 0;
   for (int lev=0; lev < levmx+1; lev++) {
      total_size +=  mesh->lev_jregsize[lev] * mesh->lev_iregsize[lev];
   //    printf("%d\n", total_size);
   }
   size_t malloc_size = total_size;

   // Allocate reg_lev arrays with that total size
   gpu_reggrid_setup(malloc_size);

   // Copy regular mesh arrays
   ezcl_enqueue_write_buffer(command_queue, dev_lev_jregmin, CL_TRUE, 0, 3*sizeof(cl_int), &mesh->lev_jregmin[0], NULL);
   ezcl_enqueue_write_buffer(command_queue, dev_lev_iregmin, CL_TRUE, 0, 3*sizeof(cl_int), &mesh->lev_iregmin[0], NULL);
   ezcl_enqueue_write_buffer(command_queue, dev_lev_jregsize, CL_TRUE, 0, 3*sizeof(cl_int), &mesh->lev_jregsize[0], NULL);
   ezcl_enqueue_write_buffer(command_queue, dev_lev_iregsize, CL_TRUE, 0, 3*sizeof(cl_int), &mesh->lev_iregsize[0], NULL);
   //ezcl_enqueue_write_buffer(command_queue, dev_j, CL_TRUE, 0, ncells*sizeof(cl_int), &mesh->j[0], NULL);
   //ezcl_enqueue_write_buffer(command_queue, dev_i, CL_TRUE, 0, ncells*sizeof(cl_int), &mesh->i[0], NULL);
   //ezcl_enqueue_write_buffer(command_queue, dev_level, CL_TRUE, 0, ncells*sizeof(cl_int), &mesh->level[0], NULL);
   //ezcl_enqueue_write_buffer(command_queue, dev_levdx, CL_TRUE, 0, 3*sizeof(cl_int), &mesh->lev_deltax[0], NULL);
   //ezcl_enqueue_write_buffer(command_queue, dev_levdy, CL_TRUE, 0, 3*sizeof(cl_int), &mesh->lev_deltay[0], NULL);

   vector<int>reg_start(levmx+1);

   // dev_H_reg_lev = ezcl_malloc(NULL, const_cast<char *>("dev_H_reg_lev"), &malloc_size, sizeof(cl_real_t), CL_MEM_READ_WRITE, 0);
   // dev_U_reg_lev = ezcl_malloc(NULL, const_cast<char *>("dev_U_reg_lev"), &malloc_size, sizeof(cl_real_t), CL_MEM_READ_WRITE, 0);
   // dev_V_reg_lev = ezcl_malloc(NULL, const_cast<char *>("dev_V_reg_lev"), &malloc_size, sizeof(cl_real_t), CL_MEM_READ_WRITE, 0);
   vector<real_t>H_reg(total_size);
   vector<real_t>U_reg(total_size);
   vector<real_t>V_reg(total_size);

   //printf("Malloc %d\n", malloc_size);
   // Get values for reg_start and also copy pstate to reg_levs
   reg_start[0] = 0;
   int startCnt = 0;
   for (int lev=0; lev < levmx+1; lev++) {
      int llsize = mesh->lev_jregsize[lev] * mesh->lev_iregsize[lev];
      //printf("startCnt %d llsize %d\n", startCnt, llsize);
      for (int copy = 0; copy < llsize; copy++) {
        H_reg[startCnt+copy] = mesh->meshes[lev].pstate[0][copy/mesh->lev_iregsize[lev]][copy%mesh->lev_iregsize[lev]];
        U_reg[startCnt+copy] = mesh->meshes[lev].pstate[1][copy/mesh->lev_iregsize[lev]][copy%mesh->lev_iregsize[lev]];
        V_reg[startCnt+copy] = mesh->meshes[lev].pstate[2][copy/mesh->lev_iregsize[lev]][copy%mesh->lev_iregsize[lev]];
      }
      startCnt += llsize;
      if (lev < levmx) reg_start[lev+1] = startCnt;
   }


      ezcl_enqueue_write_buffer(command_queue, dev_H_reg_lev, CL_TRUE, 0*sizeof(cl_real_t), total_size*sizeof(cl_real_t), &H_reg[0], NULL);
      ezcl_enqueue_write_buffer(command_queue, dev_U_reg_lev, CL_TRUE, 0*sizeof(cl_real_t), total_size*sizeof(cl_real_t), &U_reg[0], NULL);
      ezcl_enqueue_write_buffer(command_queue, dev_V_reg_lev, CL_TRUE, 0*sizeof(cl_real_t), total_size*sizeof(cl_real_t), &V_reg[0], NULL);

   // Copy up reg_start
   ezcl_enqueue_write_buffer(command_queue, dev_reg_start, CL_TRUE, 0, 3*sizeof(cl_int), &reg_start[0], NULL);

   cl_event calc_finite_difference_regular_cells_comps_event, calc_finite_difference_regular_cells_fill_event, calc_finite_difference_in_place_fixup_event;

   //for (int ic = 0; ic < ncells; ic++) {
       //int ll = mesh->level[ic];
       //int jj = mesh->j[ic] - mesh->lev_jregmin[ll];
       //int ii = mesh->i[ic] - mesh->lev_iregmin[ll];
       //printf("%d) %d %d %d\n", ic, ll, jj, ii);
       //printf("%f ", V_reg[reg_start[ll]+jj*mesh->lev_iregsize[ll]+ii]);
       //printf("%f ", V_reg[reg_start[ll]+jj*mesh->lev_iregsize[ll]+ii+1]);
       //printf("%f ", V_reg[reg_start[ll]+jj*mesh->lev_iregsize[ll]+ii-1]);
       //printf("%f ", V_reg[reg_start[ll]+(jj+1)*mesh->lev_iregsize[ll]+ii]);
       //printf("%f\n", V_reg[reg_start[ll]+(jj-1)*mesh->lev_iregsize[ll]+ii]);
        //printf("%d) ll %d jj %d ii %d startIdx %d\n", i, mesh->level[i], mesh->j[i] - mesh->lev_jregmin[mesh->level[i]], mesh->i[i] - mesh->lev_iregmin[mesh->level[i]], reg_start[mesh->level[i]]);
    //printf("%d) %f\n", i, mesh->meshes[ll].pstate[0][jj][ii]);
   //}
   for (int i = 0; i < mesh->pyfaceCnt; i++) {
//        printf("%d %d %d\n", mesh->yface_level[i], mesh->yface_i[i], mesh->yface_j[i]);
   }


   struct timespec tstart_cpu;
   cpu_timer_start(&tstart_cpu);

   /*
__kernel void calc_finite_difference_regular_cells_face_comps_cl(
                              int       ncells,                     // 0  Total number of cells.
            __global    const int       *nfaces,                    // 1 Number of faces
                              real_t    deltaT,                     // 2 Size of time step
            __global    const real_t    *lev_dx,                    // 3
            __global    const real_t    *lev_dy,                    // 4
            __global          state_t   *HxFlux,                    // 5
            __global          state_t   *UxFlux,                    // 6
            __global          state_t   *VxFlux,                    // 7
            __global          state_t   *HyFlux,                    // 8
            __global          state_t   *UyFlux,                    // 9
            __global          state_t   *VyFlux,                    // 10
            __global          state_t   *Wx_H,                      // 11
            __global          state_t   *Wx_U,                      // 12
            __global          state_t   *Wy_H,                      // 13
            __global          state_t   *Wy_V,                      // 14
            __global    const int       *reg_start,                 // 15
            __global    const real_t    *H_reg_lev,                 // 16
            __global    const real_t    *U_reg_lev,                 // 17
            __global    const real_t    *V_reg_lev,                 // 18
            __global    const int       *xface_level,               // 19
            __global    const int       *xface_j,                   // 20
            __global    const int       *xface_i,                   // 21
            __global    const int       *yface_level,               // 22
            __global    const int       *yface_j,                   // 23
            __global    const int       *yface_i,                   // 24
            __global    const int       *lev_jregmin,               // 25
            __global    const int       *lev_iregmin,               // 26
            __global    const int       *lev_jregsize,              // 27
            __global    const int       *lev_iregsize)              // 28
    */


   ezcl_set_kernel_arg(kernel_calc_finite_difference_regular_cells_face_comps, 0, sizeof(cl_int), (void *)&ncells); 
   ezcl_set_kernel_arg(kernel_calc_finite_difference_regular_cells_face_comps, 1, sizeof(cl_mem), (void *)&dev_nface); 
   ezcl_set_kernel_arg(kernel_calc_finite_difference_regular_cells_face_comps, 2, sizeof(cl_real_t), (void *)&deltaT_local); 
   ezcl_set_kernel_arg(kernel_calc_finite_difference_regular_cells_face_comps, 3, sizeof(cl_mem), (void *)&dev_levdx); 
   ezcl_set_kernel_arg(kernel_calc_finite_difference_regular_cells_face_comps, 4, sizeof(cl_mem), (void *)&dev_levdy); 
   ezcl_set_kernel_arg(kernel_calc_finite_difference_regular_cells_face_comps, 5, sizeof(cl_mem), (void *)&dev_HxFlux); 
   ezcl_set_kernel_arg(kernel_calc_finite_difference_regular_cells_face_comps, 6, sizeof(cl_mem), (void *)&dev_UxFlux); 
   ezcl_set_kernel_arg(kernel_calc_finite_difference_regular_cells_face_comps, 7, sizeof(cl_mem), (void *)&dev_VxFlux); 
   ezcl_set_kernel_arg(kernel_calc_finite_difference_regular_cells_face_comps, 8, sizeof(cl_mem), (void *)&dev_HyFlux); 
   ezcl_set_kernel_arg(kernel_calc_finite_difference_regular_cells_face_comps, 9, sizeof(cl_mem), (void *)&dev_UyFlux); 
   ezcl_set_kernel_arg(kernel_calc_finite_difference_regular_cells_face_comps, 10, sizeof(cl_mem), (void *)&dev_VyFlux); 
   ezcl_set_kernel_arg(kernel_calc_finite_difference_regular_cells_face_comps, 11, sizeof(cl_mem), (void *)&dev_Wx_H); 
   ezcl_set_kernel_arg(kernel_calc_finite_difference_regular_cells_face_comps, 12, sizeof(cl_mem), (void *)&dev_Wx_U); 
   ezcl_set_kernel_arg(kernel_calc_finite_difference_regular_cells_face_comps, 13, sizeof(cl_mem), (void *)&dev_Wy_H); 
   ezcl_set_kernel_arg(kernel_calc_finite_difference_regular_cells_face_comps, 14, sizeof(cl_mem), (void *)&dev_Wy_V); 
   ezcl_set_kernel_arg(kernel_calc_finite_difference_regular_cells_face_comps, 15, sizeof(cl_mem), (void *)&dev_reg_start); 
   ezcl_set_kernel_arg(kernel_calc_finite_difference_regular_cells_face_comps, 16, sizeof(cl_mem), (void *)&dev_H_reg_lev); 
   ezcl_set_kernel_arg(kernel_calc_finite_difference_regular_cells_face_comps, 17, sizeof(cl_mem), (void *)&dev_U_reg_lev); 
   ezcl_set_kernel_arg(kernel_calc_finite_difference_regular_cells_face_comps, 18, sizeof(cl_mem), (void *)&dev_V_reg_lev); 
   ezcl_set_kernel_arg(kernel_calc_finite_difference_regular_cells_face_comps, 19, sizeof(cl_mem), (void *)&dev_xface_level); 
   ezcl_set_kernel_arg(kernel_calc_finite_difference_regular_cells_face_comps, 20, sizeof(cl_mem), (void *)&dev_xface_j); 
   ezcl_set_kernel_arg(kernel_calc_finite_difference_regular_cells_face_comps, 21, sizeof(cl_mem), (void *)&dev_xface_i); 
   ezcl_set_kernel_arg(kernel_calc_finite_difference_regular_cells_face_comps, 22, sizeof(cl_mem), (void *)&dev_yface_level); 
   ezcl_set_kernel_arg(kernel_calc_finite_difference_regular_cells_face_comps, 23, sizeof(cl_mem), (void *)&dev_yface_j); 
   ezcl_set_kernel_arg(kernel_calc_finite_difference_regular_cells_face_comps, 24, sizeof(cl_mem), (void *)&dev_yface_i); 
   ezcl_set_kernel_arg(kernel_calc_finite_difference_regular_cells_face_comps, 25, sizeof(cl_mem), (void *)&dev_lev_jregmin); 
   ezcl_set_kernel_arg(kernel_calc_finite_difference_regular_cells_face_comps, 26, sizeof(cl_mem), (void *)&dev_lev_iregmin); 
   ezcl_set_kernel_arg(kernel_calc_finite_difference_regular_cells_face_comps, 27, sizeof(cl_mem), (void *)&dev_lev_jregsize); 
   ezcl_set_kernel_arg(kernel_calc_finite_difference_regular_cells_face_comps, 28, sizeof(cl_mem), (void *)&dev_lev_iregsize); 

   ezcl_enqueue_ndrange_kernel(command_queue, kernel_calc_finite_difference_regular_cells_face_comps, 1, NULL, &global_work_size, &local_work_size, &calc_finite_difference_regular_cells_comps_event);

   ezcl_wait_for_events(1, &calc_finite_difference_regular_cells_comps_event);
   ezcl_event_release(calc_finite_difference_regular_cells_comps_event);
/*
   vector<real_t>H_loc(ncells);
   ezcl_enqueue_read_buffer(command_queue, dev_Hyfluxminus,     CL_TRUE, 0, ncells*sizeof(cl_real_t), &H_loc[0], NULL);
   vector<real_t>U_loc(ncells);
   ezcl_enqueue_read_buffer(command_queue, dev_Uyfluxminus,     CL_TRUE, 0, ncells*sizeof(cl_real_t), &U_loc[0], NULL);
   vector<real_t>V_loc(ncells);
   ezcl_enqueue_read_buffer(command_queue, dev_Vyfluxminus,     CL_TRUE, 0, ncells*sizeof(cl_real_t), &V_loc[0], NULL);
   printf("\n");
   for (int jello = 0; jello < ncells; jello++) { printf("%d) %f | %f | %f\n", jello, H_loc[jello], U_loc[jello], V_loc[jello]); }
*/

/*
   vector<real_t>H_loc(total_size);
   ezcl_enqueue_read_buffer(command_queue, dev_H_reg_lev,     CL_TRUE, 0, total_size*sizeof(cl_real_t), &H_loc[0], NULL);
   printf("\n");
   for (int jello = 0; jello < total_size; jello++) { printf("%d) %f\n", jello, H_loc[jello]); }
*/

   /*
__kernel void calc_finite_difference_via_face_in_place_fixup_cl(
                        const int        nxfixup,                   // 0
                        const int        nyfixup,                   // 1
            __global    const int       *xrecvIdx,                  // 2
            __global    const int       *xsendIdx1,                 // 3
            __global    const int       *xsendIdx2,                 // 4
            __global    const int       *yrecvIdx,                  // 5
            __global    const int       *ysendIdx1,                 // 6
            __global    const int       *ysendIdx2,                 // 7
            __global    const int       *map_xface2cell_lower,      // 8
            __global    const int       *map_xface2cell_upper,      // 9
            __global    const int       *map_yface2cell_lower,      // 10 
            __global    const int       *map_yface2cell_upper,      // 11 
            __global          state_t   *HxFlux,                    // 12
            __global          state_t   *UxFlux,                    // 13
            __global          state_t   *VxFlux,                    // 14
            __global          state_t   *HyFlux,                    // 15
            __global          state_t   *UyFlux,                    // 16
            __global          state_t   *VyFlux,                    // 17
            __global          state_t   *Wx_H,                      // 18
            __global          state_t   *Wx_U,                      // 19
            __global          state_t   *Wy_H,                      // 20
            __global          state_t   *Wy_V) {                    // 21
    */

    ezcl_set_kernel_arg(kernel_calc_finite_difference_via_face_in_place_fixup, 0, sizeof(cl_int), (void *)&nxfixup); 
    ezcl_set_kernel_arg(kernel_calc_finite_difference_via_face_in_place_fixup, 1, sizeof(cl_int), (void *)&nyfixup); 
    ezcl_set_kernel_arg(kernel_calc_finite_difference_via_face_in_place_fixup, 2, sizeof(cl_mem), (void *)&dev_xrecvIdx);
    ezcl_set_kernel_arg(kernel_calc_finite_difference_via_face_in_place_fixup, 3, sizeof(cl_mem), (void *)&dev_xsendIdx1);
    ezcl_set_kernel_arg(kernel_calc_finite_difference_via_face_in_place_fixup, 4, sizeof(cl_mem), (void *)&dev_xsendIdx2);
    ezcl_set_kernel_arg(kernel_calc_finite_difference_via_face_in_place_fixup, 5, sizeof(cl_mem), (void *)&dev_yrecvIdx);
    ezcl_set_kernel_arg(kernel_calc_finite_difference_via_face_in_place_fixup, 6, sizeof(cl_mem), (void *)&dev_ysendIdx1);
    ezcl_set_kernel_arg(kernel_calc_finite_difference_via_face_in_place_fixup, 7, sizeof(cl_mem), (void *)&dev_ysendIdx2);
    ezcl_set_kernel_arg(kernel_calc_finite_difference_via_face_in_place_fixup, 8, sizeof(cl_mem), (void *)&dev_map_xface2cell_lower); 
    ezcl_set_kernel_arg(kernel_calc_finite_difference_via_face_in_place_fixup, 9, sizeof(cl_mem), (void *)&dev_map_xface2cell_upper); 
    ezcl_set_kernel_arg(kernel_calc_finite_difference_via_face_in_place_fixup, 10, sizeof(cl_mem), (void *)&dev_map_yface2cell_lower); 
    ezcl_set_kernel_arg(kernel_calc_finite_difference_via_face_in_place_fixup, 11, sizeof(cl_mem), (void *)&dev_map_yface2cell_upper); 
    ezcl_set_kernel_arg(kernel_calc_finite_difference_via_face_in_place_fixup, 12, sizeof(cl_mem), (void *)&dev_HxFlux); 
    ezcl_set_kernel_arg(kernel_calc_finite_difference_via_face_in_place_fixup, 13, sizeof(cl_mem), (void *)&dev_UxFlux); 
    ezcl_set_kernel_arg(kernel_calc_finite_difference_via_face_in_place_fixup, 14, sizeof(cl_mem), (void *)&dev_VxFlux); 
    ezcl_set_kernel_arg(kernel_calc_finite_difference_via_face_in_place_fixup, 15, sizeof(cl_mem), (void *)&dev_HyFlux); 
    ezcl_set_kernel_arg(kernel_calc_finite_difference_via_face_in_place_fixup, 16, sizeof(cl_mem), (void *)&dev_UyFlux); 
    ezcl_set_kernel_arg(kernel_calc_finite_difference_via_face_in_place_fixup, 17, sizeof(cl_mem), (void *)&dev_VyFlux); 
    ezcl_set_kernel_arg(kernel_calc_finite_difference_via_face_in_place_fixup, 18, sizeof(cl_mem), (void *)&dev_Wx_H); 
    ezcl_set_kernel_arg(kernel_calc_finite_difference_via_face_in_place_fixup, 19, sizeof(cl_mem), (void *)&dev_Wx_U); 
    ezcl_set_kernel_arg(kernel_calc_finite_difference_via_face_in_place_fixup, 20, sizeof(cl_mem), (void *)&dev_Wy_H); 
    ezcl_set_kernel_arg(kernel_calc_finite_difference_via_face_in_place_fixup, 21, sizeof(cl_mem), (void *)&dev_Wy_V); 

    //gpu_timers[STATE_TIMER_WRITE] += (long)(cpu_timer_stop(tstart_cpu_part)*1.0e9);
    //cpu_timer_start(&tstart_cpu_part);

    ezcl_enqueue_ndrange_kernel(command_queue, kernel_calc_finite_difference_via_face_in_place_fixup, 1, NULL, &global_work_size, &local_work_size, &calc_finite_difference_in_place_fixup_event);

    ezcl_wait_for_events(1, &calc_finite_difference_in_place_fixup_event);
    ezcl_event_release(calc_finite_difference_in_place_fixup_event);

    /*
   vector<real_t>H_loc(ncells);
   ezcl_enqueue_read_buffer(command_queue, dev_H,     CL_TRUE, 0, ncells*sizeof(cl_real_t), &H_loc[0], NULL);
   vector<real_t>U_loc(ncells);
   ezcl_enqueue_read_buffer(command_queue, dev_U,     CL_TRUE, 0, ncells*sizeof(cl_real_t), &U_loc[0], NULL);
   vector<real_t>V_loc(ncells);
   ezcl_enqueue_read_buffer(command_queue, dev_V,     CL_TRUE, 0, ncells*sizeof(cl_real_t), &V_loc[0], NULL);
   printf("\n");
   for (int jello = 0; jello < ncells; jello++) { printf("%d) %f | %f | %f\n", jello, H_loc[jello], U_loc[jello], V_loc[jello]); }
   */
/*
__kernel void calc_finite_difference_regular_cells_fill_cl(
                              int       ncells,                     // 0  Total number of cells.
                              real_t    deltaT,                     // 1 Size of time step
            __global    const real_t    *lev_dx,                    // 2
            __global    const real_t    *lev_dy,                    // 3
            __global    const state_t   *Hxfluxplus,                // 4
            __global    const state_t   *Uxfluxplus,                // 5
            __global    const state_t   *Vxfluxplus,                // 6
            __global    const state_t   *Hyfluxplus,                // 7
            __global    const state_t   *Uyfluxplus,                // 8
            __global    const state_t   *Vyfluxplus,                // 9
            __global    const state_t   *wplusx_H,                  // 10
            __global    const state_t   *wplusx_U,                  // 11
            __global    const state_t   *wplusy_H,                  // 12
            __global    const state_t   *wplusy_V,                  // 13
            __global    const int       *level,                     // 14
            __global    const int       *map_xcell2face_left1,      // 15 A cell's left primary face 
            __global    const int       *map_xcell2face_right1,     // 16 A cell's right primary face 
            __global    const int       *map_ycell2face_bot1,       // 17 A cell's bot primary face 
            __global    const int       *map_ycell2face_top1,       // 18 A cell's top primary face 
            __global    const int       *reg_start,                 // 19
            __global    const real_t    *H_reg_lev,                 // 20
            __global    const real_t    *U_reg_lev,                 // 21
            __global    const real_t    *V_reg_lev,                 // 22
            __global    const int       *j,                         // 23
            __global    const int       *i,                         // 24
            __global    const int       *lev_jregmin,               // 25
            __global    const int       *lev_iregmin,               // 26
            __global    const int       *lev_jregsize,              // 27
            __global    const int       *lev_iregsize,              // 28
            __global          real_t    *H_state_new,               // 29
            __global          real_t    *U_state_new,               // 30
            __global          real_t    *V_state_new)               // 31
*/

   ezcl_set_kernel_arg(kernel_calc_finite_difference_regular_cells_by_faces_fill, 0, sizeof(cl_int), (void *)&ncells); 
   ezcl_set_kernel_arg(kernel_calc_finite_difference_regular_cells_by_faces_fill, 1, sizeof(cl_real_t), (void *)&deltaT_local); 
   ezcl_set_kernel_arg(kernel_calc_finite_difference_regular_cells_by_faces_fill, 2, sizeof(cl_mem), (void *)&dev_levdx); 
   ezcl_set_kernel_arg(kernel_calc_finite_difference_regular_cells_by_faces_fill, 3, sizeof(cl_mem), (void *)&dev_levdy); 
   ezcl_set_kernel_arg(kernel_calc_finite_difference_regular_cells_by_faces_fill, 4, sizeof(cl_mem), (void *)&dev_HxFlux); 
   ezcl_set_kernel_arg(kernel_calc_finite_difference_regular_cells_by_faces_fill, 5, sizeof(cl_mem), (void *)&dev_UxFlux); 
   ezcl_set_kernel_arg(kernel_calc_finite_difference_regular_cells_by_faces_fill, 6, sizeof(cl_mem), (void *)&dev_VxFlux); 
   ezcl_set_kernel_arg(kernel_calc_finite_difference_regular_cells_by_faces_fill, 7, sizeof(cl_mem), (void *)&dev_HyFlux); 
   ezcl_set_kernel_arg(kernel_calc_finite_difference_regular_cells_by_faces_fill, 8, sizeof(cl_mem), (void *)&dev_UyFlux); 
   ezcl_set_kernel_arg(kernel_calc_finite_difference_regular_cells_by_faces_fill, 9, sizeof(cl_mem), (void *)&dev_VyFlux); 
   ezcl_set_kernel_arg(kernel_calc_finite_difference_regular_cells_by_faces_fill, 10, sizeof(cl_mem), (void *)&dev_Wx_H); 
   ezcl_set_kernel_arg(kernel_calc_finite_difference_regular_cells_by_faces_fill, 11, sizeof(cl_mem), (void *)&dev_Wx_U); 
   ezcl_set_kernel_arg(kernel_calc_finite_difference_regular_cells_by_faces_fill, 12, sizeof(cl_mem), (void *)&dev_Wy_H); 
   ezcl_set_kernel_arg(kernel_calc_finite_difference_regular_cells_by_faces_fill, 13, sizeof(cl_mem), (void *)&dev_Wy_V); 
   ezcl_set_kernel_arg(kernel_calc_finite_difference_regular_cells_by_faces_fill, 14, sizeof(cl_mem), (void *)&dev_level); 
   ezcl_set_kernel_arg(kernel_calc_finite_difference_regular_cells_by_faces_fill, 15, sizeof(cl_mem), (void *)&dev_map_xcell2face_left1); 
   ezcl_set_kernel_arg(kernel_calc_finite_difference_regular_cells_by_faces_fill, 16, sizeof(cl_mem), (void *)&dev_map_xcell2face_right1); 
   ezcl_set_kernel_arg(kernel_calc_finite_difference_regular_cells_by_faces_fill, 17, sizeof(cl_mem), (void *)&dev_map_ycell2face_bot1); 
   ezcl_set_kernel_arg(kernel_calc_finite_difference_regular_cells_by_faces_fill, 18, sizeof(cl_mem), (void *)&dev_map_ycell2face_top1); 
   ezcl_set_kernel_arg(kernel_calc_finite_difference_regular_cells_by_faces_fill, 19, sizeof(cl_mem), (void *)&dev_reg_start); 
   ezcl_set_kernel_arg(kernel_calc_finite_difference_regular_cells_by_faces_fill, 20, sizeof(cl_mem), (void *)&dev_H_reg_lev); 
   ezcl_set_kernel_arg(kernel_calc_finite_difference_regular_cells_by_faces_fill, 21, sizeof(cl_mem), (void *)&dev_U_reg_lev); 
   ezcl_set_kernel_arg(kernel_calc_finite_difference_regular_cells_by_faces_fill, 22, sizeof(cl_mem), (void *)&dev_V_reg_lev); 
   ezcl_set_kernel_arg(kernel_calc_finite_difference_regular_cells_by_faces_fill, 23, sizeof(cl_mem), (void *)&dev_j); 
   ezcl_set_kernel_arg(kernel_calc_finite_difference_regular_cells_by_faces_fill, 24, sizeof(cl_mem), (void *)&dev_i); 
   ezcl_set_kernel_arg(kernel_calc_finite_difference_regular_cells_by_faces_fill, 25, sizeof(cl_mem), (void *)&dev_lev_jregmin); 
   ezcl_set_kernel_arg(kernel_calc_finite_difference_regular_cells_by_faces_fill, 26, sizeof(cl_mem), (void *)&dev_lev_iregmin); 
   ezcl_set_kernel_arg(kernel_calc_finite_difference_regular_cells_by_faces_fill, 27, sizeof(cl_mem), (void *)&dev_lev_jregsize); 
   ezcl_set_kernel_arg(kernel_calc_finite_difference_regular_cells_by_faces_fill, 28, sizeof(cl_mem), (void *)&dev_lev_iregsize); 
   ezcl_set_kernel_arg(kernel_calc_finite_difference_regular_cells_by_faces_fill, 29, sizeof(cl_mem), (void *)&dev_H_new); 
   ezcl_set_kernel_arg(kernel_calc_finite_difference_regular_cells_by_faces_fill, 30, sizeof(cl_mem), (void *)&dev_U_new); 
   ezcl_set_kernel_arg(kernel_calc_finite_difference_regular_cells_by_faces_fill, 31, sizeof(cl_mem), (void *)&dev_V_new); 

   ezcl_enqueue_ndrange_kernel(command_queue, kernel_calc_finite_difference_regular_cells_by_faces_fill, 1, NULL, &global_work_size, &local_work_size, &calc_finite_difference_regular_cells_fill_event);

   ezcl_wait_for_events(1, &calc_finite_difference_regular_cells_fill_event);
   ezcl_event_release(calc_finite_difference_regular_cells_fill_event);

   gpu_timers[STATE_TIMER_FINITE_DIFFERENCE] += (long)(cpu_timer_stop(tstart_cpu)*1.0e9);

   dev_H = (cl_mem)gpu_state_memory.memory_replace(dev_H, dev_H_new);
   dev_U = (cl_mem)gpu_state_memory.memory_replace(dev_U, dev_U_new);
   dev_V = (cl_mem)gpu_state_memory.memory_replace(dev_V, dev_V_new);



   gpu_faces_delete();
   gpu_reggrid_delete();
   mesh->destroy_regular_cell_meshes(state_memory);
   mesh->gpu_wbidirmap_delete();
}
#endif

void State::symmetry_check(const char *string, vector<int> sym_index, double eps,
                           SIGN_RULE sign_rule, int &flag)
{
   size_t &ncells = mesh->ncells;

   double xsign = 1.0, ysign = 1.0;

   if (sign_rule == DIAG_RULE || sign_rule == X_RULE) {
      xsign = -1.0;
   }

   if (sign_rule == DIAG_RULE || sign_rule == Y_RULE) {
      ysign = -1.0;
   }

   for (uint ic=0; ic<ncells; ic++) {
      /*  Symmetrical check */
      if (fabs(H[ic] - H[sym_index[ic]]) > eps) {
         printf("%s ic %d sym %d H[ic] %lf Hsym %lf diff %lf\n",
                string,ic,sym_index[ic],H[ic],H[sym_index[ic]],fabs(H[ic]-H[sym_index[ic]]));
         flag++;
      }
      if (fabs(U[ic] - xsign*U[sym_index[ic]]) > eps) {
         printf("%s ic %d sym %d U[ic] %lf Usym %lf diff %lf\n",
                string,ic,sym_index[ic],U[ic],U[sym_index[ic]],fabs(U[ic]-xsign*U[sym_index[ic]]));
         flag++;
      }
      if (fabs(V[ic] - ysign*V[sym_index[ic]]) > eps) {
         printf("%s ic %d sym %d V[ic] %lf Vsym %lf diff %lf\n",
                string,ic,sym_index[ic],V[ic],V[sym_index[ic]],fabs(V[ic]-ysign*V[sym_index[ic]]));
         flag++;
      }
   }

}

size_t State::calc_refine_potential(vector<char_t> &mpot,int &icount, int &jcount)
{
   
  struct timespec tstart_cpu;
#ifdef _OPENMP
#pragma omp parallel 
{
#endif

  struct timespec tstart_lev2;

#ifdef _OPENMP
#pragma omp master
{
#endif
   cpu_timer_start(&tstart_cpu);
   if (TIMING_LEVEL >= 2) cpu_timer_start(&tstart_lev2);
#ifdef _OPENMP
}
#endif

   int *nlft, *nrht, *nbot, *ntop;
   uchar_t *level;
   
   nlft  = mesh->nlft;
   nrht  = mesh->nrht;
   nbot  = mesh->nbot;
   ntop  = mesh->ntop;
   level = mesh->level;

#ifdef _OPENMP
#pragma omp master
   {
#endif
   icount=0;
   jcount=0;
#ifdef _OPENMP
   }
#pragma omp barrier
#endif

   // We need to update the ghost regions and boundary regions for the state
   // variables since they were changed in the finite difference routine. We
   // want to use the updated values for refinement decisions
   apply_boundary_conditions();

#ifdef _OPENMP
#pragma omp barrier
#endif
/*****HIGH LEVEL OMP******/

   int lowerBound, upperBound;
   //mesh->set_bounds(ncells);
   mesh->get_bounds(lowerBound,upperBound);
#ifndef _OPENMP
#ifdef _OPENMP_SIMD
#pragma omp simd
#endif
#endif
   for (int ic=lowerBound; ic<upperBound; ic++) {

      if (mesh->celltype[ic] != REAL_CELL) continue;

      state_t Hic = H[ic];
      //state_t Uic = U[ic];
      //state_t Vic = V[ic];

      int nl = nlft[ic];
      state_t Hl = H[nl];
      //state_t Ul = U[nl];
      //state_t Vl = V[nl];

      if (level[nl] > level[ic]){
         int nlt = ntop[nl];
         Hl = REFINE_HALF * (Hl + H[nlt]);
      }

      int nr = nrht[ic];
      state_t Hr = H[nr];
      //state_t Ur = U[nr];
      //state_t Vr = V[nr];

      if (level[nr] > level[ic]){
         int nrt = ntop[nr];
         Hr = REFINE_HALF * (Hr + H[nrt]);
      }

      int nb = nbot[ic];
      state_t Hb = H[nb];
      //state_t Ub = U[nb];
      //state_t Vb = V[nb];

      if (level[nb] > level[ic]){
         int nbr = nrht[nb];
         Hb = REFINE_HALF * (Hb + H[nbr]);
      }

      int nt = ntop[ic];
      state_t Ht = H[nt];
      //state_t Ut = U[nt];
      //state_t Vt = V[nt];

      if (level[nt] > level[ic]){
         int ntr = nrht[nt];
         Ht = REFINE_HALF * (Ht + H[ntr]);
      }

      state_t duplus1; //, duplus2;
      state_t duhalf1; //, duhalf2;
      state_t duminus1; //, duminus2;

      duplus1 = Hr-Hic;
      //duplus2 = Ur-Uic;
      duhalf1 = Hic-Hl;
      //duhalf2 = Uic-Ul;

      real_t qmax = REFINE_NEG_THOUSAND;

      state_t qpot = max(fabs(duplus1/Hic), fabs(duhalf1/Hic));
      if (qpot > qmax) qmax = qpot;

      duminus1 = Hic-Hl;
      //duminus2 = Uic-Ul;
      duhalf1 = Hr-Hic;
      //duhalf2 = Ur-Uic;

      qpot = max(fabs(duminus1/Hic), fabs(duhalf1/Hic));
      if (qpot > qmax) qmax = qpot;

      duplus1 = Ht-Hic;
      //duplus2 = Vt-Vic;
      duhalf1 = Hic-Hb;
      //duhalf2 = Vic-Vb;

      qpot = max(fabs(duplus1/Hic), fabs(duhalf1/Hic));
      if (qpot > qmax) qmax = qpot;

      duminus1 = Hic-Hb;
      //duminus2 = Vic-Vb;
      duhalf1 = Ht-Hic;
      //duhalf2 = Vt-Vic;

      qpot = max(fabs(duminus1/Hic), fabs(duhalf1/Hic));
      if (qpot > qmax) qmax = qpot;

      mpot[ic]=0;
      if (qmax > REFINE_GRADIENT && (int)level[ic] < mesh->levmx) {
         mpot[ic]=1;
      } else if (qmax < COARSEN_GRADIENT && level[ic] > 0) {
         mpot[ic] = -1;
      }
      //if (mpot[ic]) printf("DEBUG cpu cell is %d mpot %d\n",ic,mpot[ic]);
   }

#ifdef _OPENMP
#pragma omp master
    {
#endif
   if (TIMING_LEVEL >= 2) {
      cpu_timers[STATE_TIMER_CALC_MPOT] += cpu_timer_stop(tstart_lev2);
   }
#ifdef _OPENMP
}
#endif

#ifdef _OPENMP
}
#pragma omp barrier
#endif
   int newcount = mesh->refine_smooth(mpot, icount, jcount);
   //printf("DEBUG -- after refine smooth in file %s line %d icount %d jcount %d newcount %d\n",__FILE__,__LINE__,icount,jcount,newcount);

   cpu_timers[STATE_TIMER_REFINE_POTENTIAL] += cpu_timer_stop(tstart_cpu);

   return(newcount);
}

#ifdef HAVE_OPENCL
size_t State::gpu_calc_refine_potential(int &icount, int &jcount)
{
   struct timespec tstart_cpu;
   cpu_timer_start(&tstart_cpu);

   struct timespec tstart_lev2;
   if (TIMING_LEVEL >= 2) cpu_timer_start(&tstart_lev2);

   cl_command_queue command_queue = ezcl_get_command_queue();

   size_t &ncells       = mesh->ncells;
   int &levmx           = mesh->levmx;
   cl_mem &dev_nlft     = mesh->dev_nlft;
   cl_mem &dev_nrht     = mesh->dev_nrht;
   cl_mem &dev_nbot     = mesh->dev_nbot;
   cl_mem &dev_ntop     = mesh->dev_ntop;
   //cl_mem &dev_mpot     = mesh->dev_mpot;
   cl_mem &dev_i        = mesh->dev_i;
   cl_mem &dev_j        = mesh->dev_j;
   cl_mem &dev_level    = mesh->dev_level;
   cl_mem &dev_celltype = mesh->dev_celltype;
   cl_mem &dev_levdx    = mesh->dev_levdx;
   cl_mem &dev_levdy    = mesh->dev_levdy;

   assert(dev_H);
   assert(dev_U);
   assert(dev_V);
   assert(dev_nlft);
   assert(dev_nrht);
   assert(dev_nbot);
   assert(dev_ntop);
   assert(dev_i);
   assert(dev_j);
   assert(dev_level);
   //assert(dev_mpot);
   //assert(dev_ioffset);
   assert(dev_levdx);
   assert(dev_levdy);

   icount = 0;
   jcount = 0;

   size_t local_work_size = 128;
   size_t global_work_size = ((ncells+local_work_size - 1) /local_work_size) * local_work_size;
   size_t block_size = global_work_size/local_work_size;

#ifdef HAVE_MPI
   //size_t nghost_local = mesh->ncells_ghost - ncells;

   if (mesh->numpe > 1) {
      ezcl_set_kernel_arg(kernel_apply_boundary_conditions_local, 0, sizeof(cl_int), &ncells);
      ezcl_set_kernel_arg(kernel_apply_boundary_conditions_local, 1, sizeof(cl_mem), &dev_celltype);
      ezcl_set_kernel_arg(kernel_apply_boundary_conditions_local, 2, sizeof(cl_mem), &dev_nlft);
      ezcl_set_kernel_arg(kernel_apply_boundary_conditions_local, 3, sizeof(cl_mem), &dev_nrht);
      ezcl_set_kernel_arg(kernel_apply_boundary_conditions_local, 4, sizeof(cl_mem), &dev_ntop);
      ezcl_set_kernel_arg(kernel_apply_boundary_conditions_local, 5, sizeof(cl_mem), &dev_nbot);
      ezcl_set_kernel_arg(kernel_apply_boundary_conditions_local, 6, sizeof(cl_mem), &dev_H);
      ezcl_set_kernel_arg(kernel_apply_boundary_conditions_local, 7, sizeof(cl_mem), &dev_U);
      ezcl_set_kernel_arg(kernel_apply_boundary_conditions_local, 8, sizeof(cl_mem), &dev_V);
      ezcl_enqueue_ndrange_kernel(command_queue, kernel_apply_boundary_conditions_local,   1, NULL, &global_work_size, &local_work_size, NULL);

      L7_Dev_Update(dev_H, L7_STATE_T, mesh->cell_handle);
      L7_Dev_Update(dev_U, L7_STATE_T, mesh->cell_handle);
      L7_Dev_Update(dev_V, L7_STATE_T, mesh->cell_handle);

      ezcl_set_kernel_arg(kernel_apply_boundary_conditions_ghost, 0, sizeof(cl_int), &ncells);
      ezcl_set_kernel_arg(kernel_apply_boundary_conditions_ghost, 1, sizeof(cl_mem), &dev_celltype);
      ezcl_set_kernel_arg(kernel_apply_boundary_conditions_ghost, 2, sizeof(cl_mem), &dev_nlft);
      ezcl_set_kernel_arg(kernel_apply_boundary_conditions_ghost, 3, sizeof(cl_mem), &dev_nrht);
      ezcl_set_kernel_arg(kernel_apply_boundary_conditions_ghost, 4, sizeof(cl_mem), &dev_ntop);
      ezcl_set_kernel_arg(kernel_apply_boundary_conditions_ghost, 5, sizeof(cl_mem), &dev_nbot);
      ezcl_set_kernel_arg(kernel_apply_boundary_conditions_ghost, 6, sizeof(cl_mem), &dev_H);
      ezcl_set_kernel_arg(kernel_apply_boundary_conditions_ghost, 7, sizeof(cl_mem), &dev_U);
      ezcl_set_kernel_arg(kernel_apply_boundary_conditions_ghost, 8, sizeof(cl_mem), &dev_V);
      ezcl_enqueue_ndrange_kernel(command_queue, kernel_apply_boundary_conditions_ghost,   1, NULL, &global_work_size, &local_work_size, NULL);
   } else {
      ezcl_set_kernel_arg(kernel_apply_boundary_conditions, 0, sizeof(cl_int), &ncells);
      ezcl_set_kernel_arg(kernel_apply_boundary_conditions, 1, sizeof(cl_mem), &dev_celltype);
      ezcl_set_kernel_arg(kernel_apply_boundary_conditions, 2, sizeof(cl_mem), &dev_nlft);
      ezcl_set_kernel_arg(kernel_apply_boundary_conditions, 3, sizeof(cl_mem), &dev_nrht);
      ezcl_set_kernel_arg(kernel_apply_boundary_conditions, 4, sizeof(cl_mem), &dev_ntop);
      ezcl_set_kernel_arg(kernel_apply_boundary_conditions, 5, sizeof(cl_mem), &dev_nbot);
      ezcl_set_kernel_arg(kernel_apply_boundary_conditions, 6, sizeof(cl_mem), &dev_H);
      ezcl_set_kernel_arg(kernel_apply_boundary_conditions, 7, sizeof(cl_mem), &dev_U);
      ezcl_set_kernel_arg(kernel_apply_boundary_conditions, 8, sizeof(cl_mem), &dev_V);
      ezcl_enqueue_ndrange_kernel(command_queue, kernel_apply_boundary_conditions,   1, NULL, &global_work_size, &local_work_size, NULL);
   }
#else
   ezcl_set_kernel_arg(kernel_apply_boundary_conditions, 0, sizeof(cl_int), &ncells);
   ezcl_set_kernel_arg(kernel_apply_boundary_conditions, 1, sizeof(cl_mem), &dev_celltype);
   ezcl_set_kernel_arg(kernel_apply_boundary_conditions, 2, sizeof(cl_mem), &dev_nlft);
   ezcl_set_kernel_arg(kernel_apply_boundary_conditions, 3, sizeof(cl_mem), &dev_nrht);
   ezcl_set_kernel_arg(kernel_apply_boundary_conditions, 4, sizeof(cl_mem), &dev_ntop);
   ezcl_set_kernel_arg(kernel_apply_boundary_conditions, 5, sizeof(cl_mem), &dev_nbot);
   ezcl_set_kernel_arg(kernel_apply_boundary_conditions, 6, sizeof(cl_mem), &dev_H);
   ezcl_set_kernel_arg(kernel_apply_boundary_conditions, 7, sizeof(cl_mem), &dev_U);
   ezcl_set_kernel_arg(kernel_apply_boundary_conditions, 8, sizeof(cl_mem), &dev_V);
   ezcl_enqueue_ndrange_kernel(command_queue, kernel_apply_boundary_conditions,   1, NULL, &global_work_size, &local_work_size, NULL);
#endif

#ifdef BOUNDS_CHECK
      {
         vector<int> nlft_tmp(mesh->ncells_ghost);
         vector<int> nrht_tmp(mesh->ncells_ghost);
         vector<int> nbot_tmp(mesh->ncells_ghost);
         vector<int> ntop_tmp(mesh->ncells_ghost);
         vector<uchar_t> level_tmp(mesh->ncells_ghost);
         vector<state_t> H_tmp(mesh->ncells_ghost);
         ezcl_enqueue_read_buffer(command_queue, dev_nlft,  CL_FALSE, 0, mesh->ncells_ghost*sizeof(cl_int), &nlft_tmp[0],  NULL);
         ezcl_enqueue_read_buffer(command_queue, dev_nrht,  CL_FALSE, 0, mesh->ncells_ghost*sizeof(cl_int), &nrht_tmp[0],  NULL);
         ezcl_enqueue_read_buffer(command_queue, dev_nbot,  CL_FALSE, 0, mesh->ncells_ghost*sizeof(cl_int), &nbot_tmp[0],  NULL);
         ezcl_enqueue_read_buffer(command_queue, dev_ntop,  CL_TRUE,  0, mesh->ncells_ghost*sizeof(cl_int), &ntop_tmp[0],  NULL);
         ezcl_enqueue_read_buffer(command_queue, dev_level, CL_TRUE,  0, mesh->ncells_ghost*sizeof(cl_uchar_t), &level_tmp[0], NULL);
         ezcl_enqueue_read_buffer(command_queue, dev_H,     CL_TRUE,  0, mesh->ncells_ghost*sizeof(cl_int), &H_tmp[0],     NULL);
         for (uint ic=0; ic<ncells; ic++){
            int nl = nlft_tmp[ic];
            if (nl<0 || nl>= (int)mesh->ncells_ghost) printf("%d: Warning at line %d cell %d nlft %d\n",mesh->mype,__LINE__,ic,nl);
            if (level_tmp[nl] > level_tmp[ic]){
               int ntl = ntop_tmp[nl];
               if (ntl<0 || ntl>= (int)mesh->ncells_ghost) printf("%d: Warning at line %d cell %d global %d nlft %d ntop of nlft %d\n",mesh->mype,__LINE__,ic,ic+mesh->noffset,nl,ntl);
            }
            int nr = nrht_tmp[ic];
            if (nr<0 || nr>= (int)mesh->ncells_ghost) printf("%d: Warning at line %d cell %d nrht %d\n",mesh->mype,__LINE__,ic,nr);
            if (level_tmp[nr] > level_tmp[ic]){
               int ntr = ntop_tmp[nr];
               if (ntr<0 || ntr>= (int)mesh->ncells_ghost) printf("%d: Warning at line %d cell %d ntop of nrht %d\n",mesh->mype,__LINE__,ic,ntr);
            }
            int nb = nbot_tmp[ic];
            if (nb<0 || nb>= (int)mesh->ncells_ghost) printf("%d: Warning at line %d cell %d nbot %d\n",mesh->mype,__LINE__,ic,nb);
            if (level_tmp[nb] > level_tmp[ic]){
               int nrb = nrht_tmp[nb];
               if (nrb<0 || nrb>= (int)mesh->ncells_ghost) printf("%d: Warning at line %d cell %d nrht of nbot %d\n",mesh->mype,__LINE__,ic,nrb);
            }
            int nt = ntop_tmp[ic];
            if (nt<0 || nt>= (int)mesh->ncells_ghost) printf("%d: Warning at line %d cell %d ntop %d\n",mesh->mype,__LINE__,ic,nt);
            if (level_tmp[nt] > level_tmp[ic]){
               int nrt = nrht_tmp[nt];
               if (nrt<0 || nrt>= (int)mesh->ncells_ghost) printf("%d: Warning at line %d cell %d nrht of ntop %d\n",mesh->mype,__LINE__,ic,nrt);
            }
         }
         for (uint ic=0; ic<mesh->ncells_ghost; ic++){
            if (H_tmp[ic] < 1.0) printf("%d: Warning at line %d cell %d H %lf\n",mesh->mype,__LINE__,ic,H_tmp[ic]);
         }
      }
#endif

   size_t result_size = 1;
   cl_mem dev_result     = ezcl_malloc(NULL, const_cast<char *>("dev_result"),     &result_size,        sizeof(cl_int2), CL_MEM_READ_WRITE, 0);
   cl_mem dev_redscratch = ezcl_malloc(NULL, const_cast<char *>("dev_redscratch"), &block_size,         sizeof(cl_int2), CL_MEM_READ_WRITE, 0);

   dev_mpot              = ezcl_malloc(NULL, const_cast<char *>("dev_mpot"),       &mesh->ncells_ghost, sizeof(cl_char_t),  CL_MEM_READ_WRITE, 0);

     /*
     __kernel void refine_potential
              const int      ncells,     // 0  Total number of cells.
              const int      levmx,      // 1  Maximum level
     __global       state_t *H,          // 2
     __global       state_t *U,          // 3
     __global       state_t *V,          // 4
     __global const int     *nlft,       // 5  Array of left neighbors.
     __global const int     *nrht,       // 6  Array of right neighbors.
     __global const int     *ntop,       // 7  Array of bottom neighbors.
     __global const int     *nbot,       // 8  Array of top neighbors.
     __global const uchar_t *level,      // 9  Array of level information.
     __global const char_t  *celltype,   // 10  Array of celltype information.
     __global       char_t  *mpot,       // 11  Array of mesh potential information.
     __global       int2    *redscratch, // 12
     __global const real_t  *lev_dx,     // 13
     __global const real_t  *lev_dy,     // 14
     __global       int2    *result,     // 15
     __local        state_t *tile,       // 16  Tile size in real4.
     __local        int8    *itile)      // 17  Tile size in int8.
     */

   ezcl_set_kernel_arg(kernel_refine_potential, 0, sizeof(cl_int),  (void *)&ncells);
   ezcl_set_kernel_arg(kernel_refine_potential, 1, sizeof(cl_int),  (void *)&levmx);
   ezcl_set_kernel_arg(kernel_refine_potential, 2, sizeof(cl_mem),  (void *)&dev_H);
   ezcl_set_kernel_arg(kernel_refine_potential, 3, sizeof(cl_mem),  (void *)&dev_U);
   ezcl_set_kernel_arg(kernel_refine_potential, 4, sizeof(cl_mem),  (void *)&dev_V);
   ezcl_set_kernel_arg(kernel_refine_potential, 5, sizeof(cl_mem),  (void *)&dev_nlft);
   ezcl_set_kernel_arg(kernel_refine_potential, 6, sizeof(cl_mem),  (void *)&dev_nrht);
   ezcl_set_kernel_arg(kernel_refine_potential, 7, sizeof(cl_mem),  (void *)&dev_ntop);
   ezcl_set_kernel_arg(kernel_refine_potential, 8, sizeof(cl_mem),  (void *)&dev_nbot);
   ezcl_set_kernel_arg(kernel_refine_potential, 9, sizeof(cl_mem),  (void *)&dev_i);
   ezcl_set_kernel_arg(kernel_refine_potential,10, sizeof(cl_mem),  (void *)&dev_j);
   ezcl_set_kernel_arg(kernel_refine_potential,11, sizeof(cl_mem),  (void *)&dev_level);
   ezcl_set_kernel_arg(kernel_refine_potential,12, sizeof(cl_mem),  (void *)&dev_celltype);
   ezcl_set_kernel_arg(kernel_refine_potential,13, sizeof(cl_mem),  (void *)&dev_levdx);
   ezcl_set_kernel_arg(kernel_refine_potential,14, sizeof(cl_mem),  (void *)&dev_levdy);
   ezcl_set_kernel_arg(kernel_refine_potential,15, sizeof(cl_mem),  (void *)&dev_mpot);
   ezcl_set_kernel_arg(kernel_refine_potential,16, sizeof(cl_mem),  (void *)&dev_redscratch);
   ezcl_set_kernel_arg(kernel_refine_potential,17, sizeof(cl_mem),  (void *)&dev_result);
   ezcl_set_kernel_arg(kernel_refine_potential,18, local_work_size*sizeof(cl_state_t),    NULL);
   ezcl_set_kernel_arg(kernel_refine_potential,19, local_work_size*sizeof(cl_int8),    NULL);

   ezcl_enqueue_ndrange_kernel(command_queue, kernel_refine_potential, 1, NULL, &global_work_size, &local_work_size, NULL);

   mesh->gpu_rezone_count2(block_size, local_work_size, dev_redscratch, dev_result);

   int count[2] = {0, 0};
   ezcl_enqueue_read_buffer(command_queue, dev_result, CL_TRUE, 0, sizeof(cl_int2), count, NULL);
   icount  = count[0];
   jcount  = count[1];
   //size_t result = ncells + icount - jcount;

   //int mpot_check[ncells];
   //ezcl_enqueue_read_buffer(command_queue, dev_mpot, CL_TRUE, 0, ncells*sizeof(cl_char_t), mpot_check, NULL);
   //for (int ic=0; ic<ncells; ic++){
   //   if (mpot_check[ic]) printf("DEBUG -- cell %d mpot %d\n",ic,mpot_check[ic]);
   //}

   //printf("result = %lu after first refine potential icount %d jcount %d\n",result, icount, jcount);
//   int which_smooth = 1;

   ezcl_device_memory_delete(dev_redscratch);
   ezcl_device_memory_delete(dev_result);

   if (TIMING_LEVEL >= 2) {
      gpu_timers[STATE_TIMER_CALC_MPOT] += (long)(cpu_timer_stop(tstart_lev2)*1.0e9);
   }

   int my_result = mesh->gpu_refine_smooth(dev_mpot, icount, jcount);
   //printf("DEBUG gpu calc refine potential %d icount %d jcount %d\n",my_result,icount,jcount);

   gpu_timers[STATE_TIMER_REFINE_POTENTIAL] += (long)(cpu_timer_stop(tstart_cpu)*1.0e9);

   return((size_t)my_result);
}
#endif

double State::mass_sum(int enhanced_precision_sum)
{
   size_t &ncells = mesh->ncells;
   char_t *celltype = mesh->celltype;
   uchar_t *level    = mesh->level;

#ifdef HAVE_MPI
   //int &mype = mesh->mype;
#endif

   struct timespec tstart_cpu;
   cpu_timer_start(&tstart_cpu);

   double summer = 0.0;
   double total_sum = 0.0;

   if (enhanced_precision_sum == SUM_KAHAN) {
      //printf("DEBUG -- kahan_sum\n");
      double corrected_next_term, new_sum;
      struct esum_type local;
#ifdef HAVE_MPI
      struct esum_type global;
#endif

      local.sum = 0.0;
      local.correction = 0.0;
      int ic;
      for (ic = 0; ic < (int)ncells; ic++) {
         if (celltype[ic] == REAL_CELL) {
            //  Exclude boundary cells.
            corrected_next_term= H[ic]*mesh->lev_deltax[level[ic]]*mesh->lev_deltay[level[ic]] + local.correction;
            new_sum            = local.sum + local.correction;
            local.correction   = corrected_next_term - (new_sum - local.sum);
            local.sum          = new_sum;
         }
      }

#ifdef HAVE_MPI
      if (mesh->parallel) {
         MPI_Allreduce(&local, &global, 1, MPI_TWO_DOUBLES, KNUTH_SUM, MPI_COMM_WORLD);
         total_sum = global.sum + global.correction;
      } else {
         total_sum = local.sum + local.correction;
      }

//if(mype == 0) printf("MYPE %d: Line %d Iteration %d \t local_sum = %12.6lg, global_sum = %12.6lg\n", mype, __LINE__, mesh->m_ncycle, local.sum, global.sum);

#else
      total_sum = local.sum + local.correction;
#endif

   } else if (enhanced_precision_sum == SUM_REGULAR) {
      //printf("DEBUG -- regular_sum\n");
      for (uint ic=0; ic < ncells; ic++){
         if (celltype[ic] == REAL_CELL) {
            summer += H[ic]*mesh->lev_deltax[level[ic]]*mesh->lev_deltay[level[ic]];
         }
      }
#ifdef HAVE_MPI
      if (mesh->parallel) {
         MPI_Allreduce(&summer, &total_sum, 1, MPI_DOUBLE, MPI_SUM, MPI_COMM_WORLD);
      } else {
         total_sum = summer;
      }
#else
      total_sum = summer;
#endif
   }

   cpu_timers[STATE_TIMER_MASS_SUM] += cpu_timer_stop(tstart_cpu);

   return(total_sum);
}

#ifdef HAVE_OPENCL
double State::gpu_mass_sum(int enhanced_precision_sum)
{
   struct timespec tstart_cpu;
   cpu_timer_start(&tstart_cpu);

   cl_command_queue command_queue = ezcl_get_command_queue();

   size_t &ncells       = mesh->ncells;
   cl_mem &dev_levdx    = mesh->dev_levdx;
   cl_mem &dev_levdy    = mesh->dev_levdy;
   cl_mem &dev_celltype = mesh->dev_celltype;
   cl_mem &dev_level    = mesh->dev_level;

   assert(dev_H);
   assert(dev_level);
   assert(dev_levdx);
   assert(dev_levdy);
   assert(dev_celltype);

   size_t one = 1;
   cl_mem dev_mass_sum, dev_redscratch;
   double gpu_mass_sum_total;

   size_t local_work_size = 128;
   size_t global_work_size = ((ncells+local_work_size - 1) /local_work_size) * local_work_size;
   size_t block_size     = global_work_size/local_work_size;

   if (enhanced_precision_sum) {
      dev_mass_sum = ezcl_malloc(NULL, const_cast<char *>("dev_mass_sum"), &one,    sizeof(cl_real2_t), CL_MEM_READ_WRITE, 0);
      dev_redscratch = ezcl_malloc(NULL, const_cast<char *>("dev_redscratch"), &block_size, sizeof(cl_real2_t), CL_MEM_READ_WRITE, 0);

        /*
        __kernel void reduce_sum_cl(
                         const int      isize,      // 0
                __global       state_t *array,      // 1   Array to be reduced.
                __global       uchar_t *level,      // 2
                __global       int     *levdx,      // 3
                __global       int     *levdy,      // 4
                __global       char_t   *celltype,   // 5
                __global       real_t  *redscratch, // 6   Final result of operation.
                __local        real_t  *tile)       // 7
        */
      ezcl_set_kernel_arg(kernel_reduce_epsum_mass_stage1of2, 0, sizeof(cl_int), (void *)&ncells);
      ezcl_set_kernel_arg(kernel_reduce_epsum_mass_stage1of2, 1, sizeof(cl_mem), (void *)&dev_H);
      ezcl_set_kernel_arg(kernel_reduce_epsum_mass_stage1of2, 2, sizeof(cl_mem), (void *)&dev_level);
      ezcl_set_kernel_arg(kernel_reduce_epsum_mass_stage1of2, 3, sizeof(cl_mem), (void *)&dev_levdx);
      ezcl_set_kernel_arg(kernel_reduce_epsum_mass_stage1of2, 4, sizeof(cl_mem), (void *)&dev_levdy);
      ezcl_set_kernel_arg(kernel_reduce_epsum_mass_stage1of2, 5, sizeof(cl_mem), (void *)&dev_celltype);
      ezcl_set_kernel_arg(kernel_reduce_epsum_mass_stage1of2, 6, sizeof(cl_mem), (void *)&dev_mass_sum);
      ezcl_set_kernel_arg(kernel_reduce_epsum_mass_stage1of2, 7, sizeof(cl_mem), (void *)&dev_redscratch);
      ezcl_set_kernel_arg(kernel_reduce_epsum_mass_stage1of2, 8, local_work_size*sizeof(cl_real2_t), NULL);

      ezcl_enqueue_ndrange_kernel(command_queue, kernel_reduce_epsum_mass_stage1of2, 1, NULL, &global_work_size, &local_work_size, NULL);

      if (block_size > 1) {
           /*
           __kernel void reduce_sum_cl(
                            const int      isize,      // 0
                   __global       int     *redscratch, // 1   Array to be reduced.
                   __local        real_t  *tile)       // 2
           */

         ezcl_set_kernel_arg(kernel_reduce_epsum_mass_stage2of2, 0, sizeof(cl_int), (void *)&block_size);
         ezcl_set_kernel_arg(kernel_reduce_epsum_mass_stage2of2, 1, sizeof(cl_mem), (void *)&dev_mass_sum);
         ezcl_set_kernel_arg(kernel_reduce_epsum_mass_stage2of2, 2, sizeof(cl_mem), (void *)&dev_redscratch);
         ezcl_set_kernel_arg(kernel_reduce_epsum_mass_stage2of2, 3, local_work_size*sizeof(cl_real2_t), NULL);

         ezcl_enqueue_ndrange_kernel(command_queue, kernel_reduce_epsum_mass_stage2of2, 1, NULL, &local_work_size, &local_work_size, NULL);
      }

      struct esum_type local, global;
      real2_t mass_sum;

      ezcl_enqueue_read_buffer(command_queue, dev_mass_sum, CL_TRUE, 0, 1*sizeof(cl_real2_t), &mass_sum, NULL);

      local.sum = mass_sum.s0;
      local.correction = mass_sum.s1;
      global.sum = local.sum;
      global.correction = local.correction;
#ifdef HAVE_MPI
      MPI_Allreduce(&local, &global, 1, MPI_TWO_DOUBLES, KNUTH_SUM, MPI_COMM_WORLD);
#endif
      gpu_mass_sum_total = global.sum + global.correction;
   } else {
      dev_mass_sum = ezcl_malloc(NULL, const_cast<char *>("dev_mass_sum"), &one,    sizeof(cl_real_t), CL_MEM_READ_WRITE, 0);
      dev_redscratch = ezcl_malloc(NULL, const_cast<char *>("dev_redscratch"), &block_size, sizeof(cl_real_t), CL_MEM_READ_WRITE, 0);

        /*
        __kernel void reduce_sum_cl(
                         const int      isize,      // 0
                __global       state_t *array,      // 1   Array to be reduced.
                __global       uchar_t *level,      // 2
                __global       int     *levdx,      // 3
                __global       int     *levdy,      // 4
                __global       char_t   *celltype,   // 5
                __global       real_t  *redscratch, // 6   Final result of operation.
                __local        real_t  *tile)       // 7
        */
      ezcl_set_kernel_arg(kernel_reduce_sum_mass_stage1of2, 0, sizeof(cl_int), (void *)&ncells);
      ezcl_set_kernel_arg(kernel_reduce_sum_mass_stage1of2, 1, sizeof(cl_mem), (void *)&dev_H);
      ezcl_set_kernel_arg(kernel_reduce_sum_mass_stage1of2, 2, sizeof(cl_mem), (void *)&dev_level);
      ezcl_set_kernel_arg(kernel_reduce_sum_mass_stage1of2, 3, sizeof(cl_mem), (void *)&dev_levdx);
      ezcl_set_kernel_arg(kernel_reduce_sum_mass_stage1of2, 4, sizeof(cl_mem), (void *)&dev_levdy);
      ezcl_set_kernel_arg(kernel_reduce_sum_mass_stage1of2, 5, sizeof(cl_mem), (void *)&dev_celltype);
      ezcl_set_kernel_arg(kernel_reduce_sum_mass_stage1of2, 6, sizeof(cl_mem), (void *)&dev_mass_sum);
      ezcl_set_kernel_arg(kernel_reduce_sum_mass_stage1of2, 7, sizeof(cl_mem), (void *)&dev_redscratch);
      ezcl_set_kernel_arg(kernel_reduce_sum_mass_stage1of2, 8, local_work_size*sizeof(cl_real_t), NULL);

      ezcl_enqueue_ndrange_kernel(command_queue, kernel_reduce_sum_mass_stage1of2, 1, NULL, &global_work_size, &local_work_size, NULL);

      if (block_size > 1) {
           /*
           __kernel void reduce_sum_cl(
                            const int     isize,      // 0
                   __global       int    *redscratch, // 1   Array to be reduced.
                   __local        real_t  *tile)       // 2
           */

         ezcl_set_kernel_arg(kernel_reduce_sum_mass_stage2of2, 0, sizeof(cl_int), (void *)&block_size);
         ezcl_set_kernel_arg(kernel_reduce_sum_mass_stage2of2, 1, sizeof(cl_mem), (void *)&dev_mass_sum);
         ezcl_set_kernel_arg(kernel_reduce_sum_mass_stage2of2, 2, sizeof(cl_mem), (void *)&dev_redscratch);
         ezcl_set_kernel_arg(kernel_reduce_sum_mass_stage2of2, 3, local_work_size*sizeof(cl_real_t), NULL);

         ezcl_enqueue_ndrange_kernel(command_queue, kernel_reduce_sum_mass_stage2of2, 1, NULL, &local_work_size, &local_work_size, NULL);
      }

      double local_sum, global_sum;
      real_t mass_sum;

      ezcl_enqueue_read_buffer(command_queue, dev_mass_sum, CL_TRUE, 0, 1*sizeof(cl_real_t), &mass_sum, NULL);
      
      local_sum = mass_sum;
      global_sum = local_sum;
#ifdef HAVE_MPI
      MPI_Allreduce(&local_sum, &global_sum, 1, MPI_DOUBLE, MPI_SUM, MPI_COMM_WORLD);
#endif
      gpu_mass_sum_total = global_sum;
   }

   ezcl_device_memory_delete(dev_redscratch);
   ezcl_device_memory_delete(dev_mass_sum);

   gpu_timers[STATE_TIMER_MASS_SUM] += (long)(cpu_timer_stop(tstart_cpu)*1.0e9);

   return(gpu_mass_sum_total);
}
#endif

#ifdef HAVE_OPENCL
void State::allocate_device_memory(size_t ncells)
{
   dev_H = (cl_mem)gpu_state_memory.memory_malloc(ncells, sizeof(cl_state_t), const_cast<char *>("dev_H"), DEVICE_REGULAR_MEMORY);
   dev_U = (cl_mem)gpu_state_memory.memory_malloc(ncells, sizeof(cl_state_t), const_cast<char *>("dev_U"), DEVICE_REGULAR_MEMORY);
   dev_V = (cl_mem)gpu_state_memory.memory_malloc(ncells, sizeof(cl_state_t), const_cast<char *>("dev_V"), DEVICE_REGULAR_MEMORY);
   //dev_Hx = (cl_mem)gpu_state_memory.memory_malloc(mesh->nxface, sizeof(cl_state_t), const_cast<char *>("dev_Hx"), DEVICE_REGULAR_MEMORY);
   //dev_Ux = (cl_mem)gpu_state_memory.memory_malloc(mesh->nxface, sizeof(cl_state_t), const_cast<char *>("dev_Ux"), DEVICE_REGULAR_MEMORY);
   //dev_Vx = (cl_mem)gpu_state_memory.memory_malloc(mesh->nxface, sizeof(cl_state_t), const_cast<char *>("dev_Vx"), DEVICE_REGULAR_MEMORY);
   //dev_Hy = (cl_mem)gpu_state_memory.memory_malloc(mesh->nyface, sizeof(cl_state_t), const_cast<char *>("dev_Hy"), DEVICE_REGULAR_MEMORY);
   //dev_Uy = (cl_mem)gpu_state_memory.memory_malloc(mesh->nyface, sizeof(cl_state_t), const_cast<char *>("dev_Uy"), DEVICE_REGULAR_MEMORY);
   //dev_Vy = (cl_mem)gpu_state_memory.memory_malloc(mesh->nyface, sizeof(cl_state_t), const_cast<char *>("dev_Vy"), DEVICE_REGULAR_MEMORY);
}
#endif

void State::resize_old_device_memory(size_t ncells)
{
#ifdef HAVE_OPENCL
   //gpu_state_memory.memory_delete(dev_H);
   //gpu_state_memory.memory_delete(dev_U);
   //gpu_state_memory.memory_delete(dev_V);
   //gpu_state_memory.memory_delete(dev_Hx);
   //gpu_state_memory.memory_delete(dev_Ux);
   //gpu_state_memory.memory_delete(dev_Vx);
   //gpu_state_memory.memory_delete(dev_Hy);
   //gpu_state_memory.memory_delete(dev_Uy);
   //gpu_state_memory.memory_delete(dev_Vy);
   dev_H = (cl_mem)gpu_state_memory.memory_malloc(ncells, sizeof(cl_state_t), const_cast<char *>("dev_H"), DEVICE_REGULAR_MEMORY);
   dev_U = (cl_mem)gpu_state_memory.memory_malloc(ncells, sizeof(cl_state_t), const_cast<char *>("dev_U"), DEVICE_REGULAR_MEMORY);
   dev_V = (cl_mem)gpu_state_memory.memory_malloc(ncells, sizeof(cl_state_t), const_cast<char *>("dev_V"), DEVICE_REGULAR_MEMORY);
   //dev_Hx = (cl_mem)gpu_state_memory.memory_malloc(mesh->nxface, sizeof(cl_state_t), const_cast<char *>("dev_Hx"), DEVICE_REGULAR_MEMORY);
   //dev_Ux = (cl_mem)gpu_state_memory.memory_malloc(mesh->nxface, sizeof(cl_state_t), const_cast<char *>("dev_Ux"), DEVICE_REGULAR_MEMORY);
   //dev_Vx = (cl_mem)gpu_state_memory.memory_malloc(mesh->nxface, sizeof(cl_state_t), const_cast<char *>("dev_Vx"), DEVICE_REGULAR_MEMORY);
   //dev_Hy = (cl_mem)gpu_state_memory.memory_malloc(mesh->nyface, sizeof(cl_state_t), const_cast<char *>("dev_Hy"), DEVICE_REGULAR_MEMORY);
   //dev_Uy = (cl_mem)gpu_state_memory.memory_malloc(mesh->nyface, sizeof(cl_state_t), const_cast<char *>("dev_Uy"), DEVICE_REGULAR_MEMORY);
   //dev_Vy = (cl_mem)gpu_state_memory.memory_malloc(mesh->nyface, sizeof(cl_state_t), const_cast<char *>("dev_Vy"), DEVICE_REGULAR_MEMORY);
#else
   // Just to block compiler warnings
   //if (1 == 2) printf("DEBUG -- ncells is %ld\n",ncells);
#endif
}

#ifdef HAVE_MPI
void State::do_load_balance_local(size_t &numcells){
   mesh->do_load_balance_local(numcells, NULL, state_memory);
   memory_reset_ptrs();
}
#endif
#ifdef HAVE_OPENCL
#ifdef HAVE_MPI
void State::gpu_do_load_balance_local(size_t &numcells){
   if (mesh->gpu_do_load_balance_local(numcells, NULL, gpu_state_memory) ){
      //gpu_state_memory.memory_report();
      dev_H = (cl_mem)gpu_state_memory.get_memory_ptr("dev_H");
      dev_U = (cl_mem)gpu_state_memory.get_memory_ptr("dev_U");
      dev_V = (cl_mem)gpu_state_memory.get_memory_ptr("dev_V");
/*
      if (dev_H == NULL){
         dev_H = (cl_mem)gpu_state_memory.get_memory_ptr("dev_H_new");
         dev_U = (cl_mem)gpu_state_memory.get_memory_ptr("dev_U_new");
         dev_V = (cl_mem)gpu_state_memory.get_memory_ptr("dev_V_new");
      }
      printf("DEBUG memory for proc %d dev_H is %p dev_U is %p dev_V is %p\n",mesh->mype,dev_H,dev_U,dev_V);
*/
   }
}
#endif
#endif

static double reference_time = 0.0;

void State::output_timing_info(int do_cpu_calc, int do_gpu_calc, double total_elapsed_time)
{
   int parallel = mesh->parallel;

   double cpu_time_compute = 0.0;
   double gpu_time_compute = 0.0;

   double cpu_elapsed_time = 0.0;
   double gpu_elapsed_time = 0.0;

   double cpu_mesh_time = 0.0;
   double gpu_mesh_time = 0.0;

   if (do_cpu_calc) {
      cpu_time_compute = get_cpu_timer(STATE_TIMER_SET_TIMESTEP) +
                         get_cpu_timer(STATE_TIMER_FINITE_DIFFERENCE) +
                         get_cpu_timer(STATE_TIMER_REFINE_POTENTIAL) +
                         get_cpu_timer(STATE_TIMER_REZONE_ALL) +
                         mesh->get_cpu_timer(MESH_TIMER_CALC_NEIGHBORS) +
                         mesh->get_cpu_timer(MESH_TIMER_LOAD_BALANCE) +
                         get_cpu_timer(STATE_TIMER_MASS_SUM) +
                         mesh->get_cpu_timer(MESH_TIMER_CALC_SPATIAL_COORDINATES) +
                         mesh->get_cpu_timer(MESH_TIMER_PARTITION);
      cpu_elapsed_time = cpu_time_compute;
      cpu_mesh_time = mesh->get_cpu_timer(MESH_TIMER_CALC_NEIGHBORS) +
                      get_cpu_timer(STATE_TIMER_REZONE_ALL) +
                      mesh->get_cpu_timer(MESH_TIMER_REFINE_SMOOTH) +
                      mesh->get_cpu_timer(MESH_TIMER_LOAD_BALANCE);
   }
   if (do_gpu_calc) {
      gpu_time_compute = get_gpu_timer(STATE_TIMER_APPLY_BCS) +
                         get_gpu_timer(STATE_TIMER_SET_TIMESTEP) +
                         get_gpu_timer(STATE_TIMER_FINITE_DIFFERENCE) +
                         get_gpu_timer(STATE_TIMER_REFINE_POTENTIAL) +
                         get_gpu_timer(STATE_TIMER_REZONE_ALL) +
                         mesh->get_gpu_timer(MESH_TIMER_CALC_NEIGHBORS) +
                         mesh->get_gpu_timer(MESH_TIMER_LOAD_BALANCE) +
                         get_gpu_timer(STATE_TIMER_MASS_SUM) +
                         mesh->get_gpu_timer(MESH_TIMER_CALC_SPATIAL_COORDINATES) +
                         mesh->get_gpu_timer(MESH_TIMER_COUNT_BCS);
      gpu_elapsed_time = get_gpu_timer(STATE_TIMER_WRITE) + gpu_time_compute + get_gpu_timer(STATE_TIMER_READ);
      gpu_mesh_time = mesh->get_gpu_timer(MESH_TIMER_CALC_NEIGHBORS) +
                      get_gpu_timer(STATE_TIMER_REZONE_ALL) +
                      mesh->get_gpu_timer(MESH_TIMER_REFINE_SMOOTH) +
                      mesh->get_gpu_timer(MESH_TIMER_LOAD_BALANCE);
   }

   if (! parallel && do_cpu_calc) reference_time = cpu_elapsed_time;

   double speedup_ratio = 0.0;
   if (reference_time > 0.0){
      if (do_cpu_calc && parallel) speedup_ratio = reference_time/cpu_elapsed_time;
      if (do_gpu_calc) speedup_ratio = reference_time/gpu_elapsed_time;
   }

   if (do_cpu_calc) {
      output_timer_block(MESH_DEVICE_CPU, cpu_elapsed_time, cpu_mesh_time, cpu_time_compute, total_elapsed_time, speedup_ratio);
   }
   if (do_gpu_calc) {
      output_timer_block(MESH_DEVICE_GPU, gpu_elapsed_time, gpu_mesh_time, gpu_time_compute, total_elapsed_time, speedup_ratio);
   }
}

void State::output_timer_block(mesh_device_types device_type, double elapsed_time,
   double mesh_time, double compute_time, double total_elapsed_time, double speedup_ratio)
{
   int mype  = mesh->mype;
   int parallel = mesh->parallel;

   int rank = mype;
   if (! parallel) {
      // We need to get rank info for check routines
#ifdef HAVE_MPI
      MPI_Comm_rank(MPI_COMM_WORLD, &rank);
#endif
   }

   if (! parallel && rank) return;

   char device_string[10];
   if (device_type == MESH_DEVICE_CPU) {
      sprintf(device_string,"CPU");
   } else {
      sprintf(device_string,"GPU");
   }

   if (rank == 0) {
      printf("\n");
      printf("~~~~~~~~~~~~~~~~ Device timing information ~~~~~~~~~~~~~~~~~~\n");
   }

   if (rank == 0 && parallel) {
      printf("\n%3s: Parallel timings\n\n",device_string);
   }

   if (device_type == MESH_DEVICE_GPU) {
      mesh->parallel_output("GPU: Write to device          time was",  get_gpu_timer(STATE_TIMER_WRITE), 0, "s");
      mesh->parallel_output("GPU: Read from device         time was",  get_gpu_timer(STATE_TIMER_READ),  0, "s");
   }

   const char *device_compute_string[2] = {
      "CPU: Device compute           time was",
      "GPU: Device compute           time was"
   };
   mesh->parallel_output(device_compute_string[device_type], compute_time, 0, "s");

   timer_output(STATE_TIMER_SET_TIMESTEP,                  device_type, 1);
   timer_output(STATE_TIMER_FINITE_DIFFERENCE,             device_type, 1);
   timer_output(STATE_TIMER_FINITE_DIFFERENCE_PART1,       device_type, 2);
   timer_output(STATE_TIMER_FINITE_DIFFERENCE_PART2,       device_type, 2);
   timer_output(STATE_TIMER_FINITE_DIFFERENCE_PART3,       device_type, 2);
   timer_output(STATE_TIMER_FINITE_DIFFERENCE_PART4,       device_type, 2);
   timer_output(STATE_TIMER_FINITE_DIFFERENCE_PART5,       device_type, 2);
   timer_output(STATE_TIMER_FINITE_DIFFERENCE_PART6,       device_type, 2);
   mesh->timer_output(MESH_TIMER_BIDIR,                    device_type, 3);
   mesh->timer_output(MESH_TIMER_BIDIRPART1,               device_type, 3);
   mesh->timer_output(MESH_TIMER_BIDIRPART2,               device_type, 3);
   mesh->timer_output(MESH_TIMER_BIDIRPART3,               device_type, 3);
   mesh->timer_output(MESH_TIMER_BIDIRPART4,               device_type, 3);
   mesh->timer_output(MESH_TIMER_BIDIRPART5,               device_type, 3);
   mesh->timer_output(MESH_TIMER_BIDIRPART6,               device_type, 3);
   mesh->timer_output(MESH_TIMER_BIDIRPART7,               device_type, 3);
   mesh->timer_output(MESH_TIMER_BIDIRPART8,               device_type, 3);
   mesh->timer_output(MESH_TIMER_BIDIRPART9,               device_type, 3);
   mesh->timer_output(MESH_TIMER_BIDIRPART10,              device_type, 3);
   mesh->timer_output(MESH_TIMER_BIDIRPART11,              device_type, 3);
   mesh->timer_output(MESH_TIMER_BIDIRPART12,              device_type, 3);
   timer_output(STATE_TIMER_REFINE_POTENTIAL,              device_type, 1);
   timer_output(STATE_TIMER_CALC_MPOT,                     device_type, 2);
   mesh->timer_output(MESH_TIMER_REFINE_SMOOTH,            device_type, 2);
   timer_output(STATE_TIMER_REZONE_ALL,                    device_type, 1);
   mesh->timer_output(MESH_TIMER_PARTITION,                device_type, 1);
   mesh->timer_output(MESH_TIMER_CALC_NEIGHBORS,           device_type, 1);
   if (mesh->get_calc_neighbor_type() == HASH_TABLE) {
      mesh->timer_output(MESH_TIMER_HASH_SETUP,            device_type, 2);
      mesh->timer_output(MESH_TIMER_HASH_QUERY,            device_type, 2);
      if (parallel) {
         mesh->timer_output(MESH_TIMER_FIND_BOUNDARY,      device_type, 2);
         mesh->timer_output(MESH_TIMER_PUSH_SETUP,         device_type, 2);
         mesh->timer_output(MESH_TIMER_PUSH_BOUNDARY,      device_type, 2);
         mesh->timer_output(MESH_TIMER_LOCAL_LIST,         device_type, 2);
         mesh->timer_output(MESH_TIMER_LAYER1,             device_type, 2);
         mesh->timer_output(MESH_TIMER_LAYER2,             device_type, 2);
         mesh->timer_output(MESH_TIMER_LAYER_LIST,         device_type, 2);
         mesh->timer_output(MESH_TIMER_COPY_MESH_DATA,     device_type, 2);
         mesh->timer_output(MESH_TIMER_FILL_MESH_GHOST,    device_type, 2);
         mesh->timer_output(MESH_TIMER_FILL_NEIGH_GHOST,   device_type, 2);
         mesh->timer_output(MESH_TIMER_SET_CORNER_NEIGH,   device_type, 2);
         mesh->timer_output(MESH_TIMER_NEIGH_ADJUST,       device_type, 2);
         mesh->timer_output(MESH_TIMER_SETUP_COMM,         device_type, 2);
      }
   } else {
      mesh->timer_output(MESH_TIMER_KDTREE_SETUP,          device_type, 2);
      mesh->timer_output(MESH_TIMER_KDTREE_QUERY,          device_type, 2);
   }
   timer_output(STATE_TIMER_MASS_SUM,                      device_type, 1);
   if (parallel) {
      mesh->timer_output(MESH_TIMER_LOAD_BALANCE,          device_type, 1);
   }
   mesh->timer_output(MESH_TIMER_CALC_SPATIAL_COORDINATES, device_type, 1);
   if (! mesh->have_boundary) {
      mesh->timer_output(MESH_TIMER_COUNT_BCS,             device_type, 1);
   }
   if (rank == 0) printf("=============================================================\n");

   const char *profile_string[2] = {
      "Profiling: Total CPU          time was",
      "Profiling: Total GPU          time was"
   };
   mesh->parallel_output(profile_string[device_type], elapsed_time, 0, "s");
   if (elapsed_time > 600.0){
      mesh->parallel_output("                                  or  ", elapsed_time/60.0, 0, "min");
   }

   if (rank == 0) printf("-------------------------------------------------------------\n");
   mesh->parallel_output("Mesh Ops (Neigh+rezone+smooth+balance) ",mesh_time, 0, "s");
   mesh->parallel_output("Mesh Ops Percentage                    ",mesh_time/elapsed_time*100.0, 0, "percent");
   if (rank == 0) printf("=============================================================\n");

   mesh->parallel_output("Profiling: Total              time was",total_elapsed_time, 0, "s");
   if (elapsed_time > 600.0){
      mesh->parallel_output("                                  or  ",total_elapsed_time/60.0, 0, "min");
   }

   if (speedup_ratio > 0.0) {
      mesh->parallel_output("Parallel Speed-up:                    ",speedup_ratio, 0, "Reference Serial CPU");
   }

   if (rank == 0) printf("=============================================================\n");
}

void State::timer_output(state_timer_category category, mesh_device_types device_type, int timer_level)
{
   int mype = mesh->mype;

   double local_time = 0.0;
   if (device_type == MESH_DEVICE_CPU){
      local_time = get_cpu_timer(category);
   } else {
      local_time = get_gpu_timer(category);
   }

   char string[80] = "/0";

   if (mype == 0) {
      const char *blank="          ";

      const char *device_string[2] = {
         "CPU",
         "GPU"
      };

      sprintf(string,"%3s: %.*s%-30.30s\t", device_string[device_type],
         2*timer_level, blank, state_timer_descriptor[category]);
   }

   mesh->parallel_output(string, local_time, timer_level, "s");
}

#ifdef HAVE_OPENCL
void State::compare_state_gpu_global_to_cpu_global(const char* string, int cycle, uint ncells)
{
   cl_command_queue command_queue = ezcl_get_command_queue();

   vector<state_t>H_check(ncells);
   vector<state_t>U_check(ncells);
   vector<state_t>V_check(ncells);
   ezcl_enqueue_read_buffer(command_queue, dev_H, CL_FALSE, 0, ncells*sizeof(cl_state_t), &H_check[0], NULL);
   ezcl_enqueue_read_buffer(command_queue, dev_U, CL_FALSE, 0, ncells*sizeof(cl_state_t), &U_check[0], NULL);
   ezcl_enqueue_read_buffer(command_queue, dev_V, CL_TRUE,  0, ncells*sizeof(cl_state_t), &V_check[0], NULL);
   for (uint ic = 0; ic < ncells; ic++){
      if (fabs(H[ic]-H_check[ic]) > STATE_EPS) printf("DEBUG %s at cycle %d H & H_check %d %lf %lf\n",string,cycle,ic,H[ic],H_check[ic]);
      if (fabs(U[ic]-U_check[ic]) > STATE_EPS) printf("DEBUG %s at cycle %d U & U_check %d %lf %lf\n",string,cycle,ic,U[ic],U_check[ic]);
      if (fabs(V[ic]-V_check[ic]) > STATE_EPS) printf("DEBUG %s at cycle %d V & V_check %d %lf %lf\n",string,cycle,ic,V[ic],V_check[ic]);
   }
}
#endif

void State::compare_state_cpu_local_to_cpu_global(State *state_global, const char* string, int cycle, uint ncells, uint ncells_global, int *nsizes, int *ndispl)
{
   state_t *H_global = state_global->H;
   state_t *U_global = state_global->U;
   state_t *V_global = state_global->V;

   vector<state_t>H_check(ncells_global);
   vector<state_t>U_check(ncells_global);
   vector<state_t>V_check(ncells_global);
#ifdef HAVE_MPI
   MPI_Allgatherv(&H[0], ncells, MPI_STATE_T, &H_check[0], &nsizes[0], &ndispl[0], MPI_STATE_T, MPI_COMM_WORLD);
   MPI_Allgatherv(&U[0], ncells, MPI_STATE_T, &U_check[0], &nsizes[0], &ndispl[0], MPI_STATE_T, MPI_COMM_WORLD);
   MPI_Allgatherv(&V[0], ncells, MPI_STATE_T, &V_check[0], &nsizes[0], &ndispl[0], MPI_STATE_T, MPI_COMM_WORLD);
#else
   // Just to block compiler warnings
   //if (1 == 2) printf("DEBUG -- ncells %u nsizes %d ndispl %d\n",ncells, nsizes[0],ndispl[0]);
#endif

   if (mesh->mype == 0) {
      for (uint ic = 0; ic < ncells_global; ic++){
         if (fabs(H_global[ic]-H_check[ic]) > STATE_EPS) printf("DEBUG %s at cycle %d H & H_check %d %lf %lf\n",string,cycle,ic,H_global[ic],H_check[ic]);
         if (fabs(U_global[ic]-U_check[ic]) > STATE_EPS) printf("DEBUG %s at cycle %d U & U_check %d %lf %lf\n",string,cycle,ic,U_global[ic],U_check[ic]);
         if (fabs(V_global[ic]-V_check[ic]) > STATE_EPS) printf("DEBUG %s at cycle %d V & V_check %d %lf %lf\n",string,cycle,ic,V_global[ic],V_check[ic]);
      }
   }
}

#ifdef HAVE_OPENCL
void State::compare_state_all_to_gpu_local(State *state_global, uint ncells, uint ncells_global, int mype, int ncycle, int *nsizes, int *ndispl)
{
#ifdef HAVE_MPI
   cl_command_queue command_queue = ezcl_get_command_queue();

   state_t *H_global = state_global->H;
   state_t *U_global = state_global->U;
   state_t *V_global = state_global->V;
   cl_mem &dev_H_global = state_global->dev_H;
   cl_mem &dev_U_global = state_global->dev_U;
   cl_mem &dev_V_global = state_global->dev_V;

   // Need to compare dev_H to H, etc
   vector<state_t>H_save(ncells);
   vector<state_t>U_save(ncells);
   vector<state_t>V_save(ncells);
   ezcl_enqueue_read_buffer(command_queue, dev_H, CL_FALSE, 0, ncells*sizeof(cl_state_t), &H_save[0], NULL);
   ezcl_enqueue_read_buffer(command_queue, dev_U, CL_FALSE, 0, ncells*sizeof(cl_state_t), &U_save[0], NULL);
   ezcl_enqueue_read_buffer(command_queue, dev_V, CL_TRUE,  0, ncells*sizeof(cl_state_t), &V_save[0], NULL);
   for (uint ic = 0; ic < ncells; ic++){
      if (fabs(H[ic]-H_save[ic]) > STATE_EPS) printf("%d: DEBUG finite_difference 1 at cycle %d H & H_save %d %lf %lf \n",mype,ncycle,ic,H[ic],H_save[ic]);
      if (fabs(U[ic]-U_save[ic]) > STATE_EPS) printf("%d: DEBUG finite_difference 1 at cycle %d U & U_save %d %lf %lf \n",mype,ncycle,ic,U[ic],U_save[ic]);
      if (fabs(V[ic]-V_save[ic]) > STATE_EPS) printf("%d: DEBUG finite_difference 1 at cycle %d V & V_save %d %lf %lf \n",mype,ncycle,ic,V[ic],V_save[ic]);
   }

   // And compare dev_H gathered to H_global, etc
   vector<state_t>H_save_global(ncells_global);
   vector<state_t>U_save_global(ncells_global);
   vector<state_t>V_save_global(ncells_global);
   MPI_Allgatherv(&H_save[0], nsizes[mype], MPI_STATE_T, &H_save_global[0], &nsizes[0], &ndispl[0], MPI_STATE_T, MPI_COMM_WORLD);
   MPI_Allgatherv(&U_save[0], nsizes[mype], MPI_STATE_T, &U_save_global[0], &nsizes[0], &ndispl[0], MPI_STATE_T, MPI_COMM_WORLD);
   MPI_Allgatherv(&V_save[0], nsizes[mype], MPI_STATE_T, &V_save_global[0], &nsizes[0], &ndispl[0], MPI_STATE_T, MPI_COMM_WORLD);
   if (mype == 0) {
      for (uint ic = 0; ic < ncells_global; ic++){
         if (fabs(H_global[ic]-H_save_global[ic]) > STATE_EPS) printf("%d: DEBUG finite_difference 2 at cycle %d H_global & H_save_global %d %lf %lf \n",mype,ncycle,ic,H_global[ic],H_save_global[ic]);
         if (fabs(U_global[ic]-U_save_global[ic]) > STATE_EPS) printf("%d: DEBUG finite_difference 2 at cycle %d U_global & U_save_global %d %lf %lf \n",mype,ncycle,ic,U_global[ic],U_save_global[ic]);
         if (fabs(V_global[ic]-V_save_global[ic]) > STATE_EPS) printf("%d: DEBUG finite_difference 2 at cycle %d V_global & V_save_global %d %lf %lf \n",mype,ncycle,ic,V_global[ic],V_save_global[ic]);
      }
   }

   // And compare H gathered to H_global, etc
   MPI_Allgatherv(&H[0], nsizes[mype], MPI_STATE_T, &H_save_global[0], &nsizes[0], &ndispl[0], MPI_STATE_T, MPI_COMM_WORLD);
   MPI_Allgatherv(&U[0], nsizes[mype], MPI_STATE_T, &U_save_global[0], &nsizes[0], &ndispl[0], MPI_STATE_T, MPI_COMM_WORLD);
   MPI_Allgatherv(&V[0], nsizes[mype], MPI_STATE_T, &V_save_global[0], &nsizes[0], &ndispl[0], MPI_STATE_T, MPI_COMM_WORLD);
   if (mype == 0) {
      for (uint ic = 0; ic < ncells_global; ic++){
         if (fabs(H_global[ic]-H_save_global[ic]) > STATE_EPS) printf("DEBUG finite_difference 3 at cycle %d H_global & H_save_global %d %lf %lf \n",ncycle,ic,H_global[ic],H_save_global[ic]);
         if (fabs(U_global[ic]-U_save_global[ic]) > STATE_EPS) printf("DEBUG finite_difference 3 at cycle %d U_global & U_save_global %d %lf %lf \n",ncycle,ic,U_global[ic],U_save_global[ic]);
         if (fabs(V_global[ic]-V_save_global[ic]) > STATE_EPS) printf("DEBUG finite_difference 3 at cycle %d V_global & V_save_global %d %lf %lf \n",ncycle,ic,V_global[ic],V_save_global[ic]);
      }
   }

   // Now the global dev_H_global to H_global, etc
   ezcl_enqueue_read_buffer(command_queue, dev_H_global, CL_FALSE, 0, ncells_global*sizeof(cl_state_t), &H_save_global[0], NULL);
   ezcl_enqueue_read_buffer(command_queue, dev_U_global, CL_FALSE, 0, ncells_global*sizeof(cl_state_t), &U_save_global[0], NULL);
   ezcl_enqueue_read_buffer(command_queue, dev_V_global, CL_TRUE,  0, ncells_global*sizeof(cl_state_t), &V_save_global[0], NULL);
   if (mype == 0) {
      for (uint ic = 0; ic < ncells_global; ic++){
         if (fabs(H_global[ic]-H_save_global[ic]) > STATE_EPS) printf("%d: DEBUG finite_difference 4 at cycle %d H_global & H_save_global %d %lf %lf \n",mype,ncycle,ic,H_global[ic],H_save_global[ic]);
         if (fabs(U_global[ic]-U_save_global[ic]) > STATE_EPS) printf("%d: DEBUG finite_difference 4 at cycle %d U_global & U_save_global %d %lf %lf \n",mype,ncycle,ic,U_global[ic],U_save_global[ic]);
         if (fabs(V_global[ic]-V_save_global[ic]) > STATE_EPS) printf("%d: DEBUG finite_difference 4 at cycle %d V_global & V_save_global %d %lf %lf \n",mype,ncycle,ic,V_global[ic],V_save_global[ic]);
      }
   }
#else
   // Just to get rid of compiler warnings
   //if (1 == 2) printf("%d: DEBUG -- ncells %d ncells_global %d ncycle %d nsizes[0] %d ndispl %d state_global %p\n",
   //   mype,ncells,ncells_global,ncycle,nsizes[0],ndispl[0],state_global);
#endif
}
#endif

void State::print_object_info(void)
{
   printf(" ---- State object info -----\n");

#ifdef HAVE_OPENCL
   int num_elements, elsize;

   num_elements = ezcl_get_device_mem_nelements(dev_H);
   elsize = ezcl_get_device_mem_elsize(dev_H);
   printf("dev_H       ptr : %p nelements %d elsize %d\n",dev_H,num_elements,elsize);
   num_elements = ezcl_get_device_mem_nelements(dev_U);
   elsize = ezcl_get_device_mem_elsize(dev_U);
   printf("dev_U       ptr : %p nelements %d elsize %d\n",dev_U,num_elements,elsize);
   num_elements = ezcl_get_device_mem_nelements(dev_V);
   elsize = ezcl_get_device_mem_elsize(dev_V);
   printf("dev_V       ptr : %p nelements %d elsize %d\n",dev_V,num_elements,elsize);
   num_elements = ezcl_get_device_mem_nelements(dev_mpot);
   elsize = ezcl_get_device_mem_elsize(dev_mpot);
   printf("dev_mpot    ptr : %p nelements %d elsize %d\n",dev_mpot,num_elements,elsize);
   //num_elements = ezcl_get_device_mem_nelements(dev_ioffset);
   //elsize = ezcl_get_device_mem_elsize(dev_ioffset);
   //printf("dev_ioffset ptr : %p nelements %d elsize %d\n",dev_ioffset,num_elements,elsize);
#endif
   state_memory.memory_report();
   //printf("vector H    ptr : %p nelements %ld elsize %ld\n",&H[0],H.size(),sizeof(H[0]));
   //printf("vector U    ptr : %p nelements %ld elsize %ld\n",&U[0],U.size(),sizeof(U[0]));
   //printf("vector V    ptr : %p nelements %ld elsize %ld\n",&V[0],V.size(),sizeof(V[0]));
}

void State::print(void)
{  //printf("size is %lu %lu %lu %lu %lu\n",index.size(), i.size(), level.size(), nlft.size(), x.size());

   if (mesh->fp == NULL) {
      char filename[10];
      sprintf(filename,"out%1d",mesh->mype);
      mesh->fp=fopen(filename,"w");
   }

   if (mesh->mesh_memory.get_memory_size(mesh->nlft) >= mesh->ncells_ghost){
      fprintf(mesh->fp,"%d:   index global  i     j     lev   nlft  nrht  nbot  ntop \n",mesh->mype);
      for (uint ic=0; ic<mesh->ncells; ic++) {
         fprintf(mesh->fp,"%d: %6d  %6d %4d  %4d   %4d  %4d  %4d  %4d  %4d \n", mesh->mype,ic, ic+mesh->noffset,mesh->i[ic], mesh->j[ic], mesh->level[ic], mesh->nlft[ic], mesh->nrht[ic], mesh->nbot[ic], mesh->ntop[ic]);
      }
      for (uint ic=mesh->ncells; ic<mesh->ncells_ghost; ic++) {
         fprintf(mesh->fp,"%d: %6d  %6d %4d  %4d   %4d  %4d  %4d  %4d  %4d \n", mesh->mype,ic, ic+mesh->noffset,mesh->i[ic], mesh->j[ic], mesh->level[ic], mesh->nlft[ic], mesh->nrht[ic], mesh->nbot[ic], mesh->ntop[ic]);
      }
   } else {
      fprintf(mesh->fp,"%d:  index     H        U         V      i     j     lev\n",mesh->mype);
      for (uint ic=0; ic<mesh->ncells_ghost; ic++) {
         fprintf(mesh->fp,"%d: %6d %lf %lf %lf %4d  %4d   %4d  \n", mesh->mype,ic, H[ic], U[ic], V[ic], mesh->i[ic], mesh->j[ic], mesh->level[ic]);
      }
   }
   fflush(mesh->fp);
}

void State::print_data_dump(int ncycle)

{  //printf("size is %lu %lu %lu %lu %lu\n",index.size(), i.size(), level.size(), nlft.size(), x.size());

   char filename[20];
   sprintf(filename,"out%1d.%05d",mesh->mype,ncycle);
   FILE *fp=fopen(filename,"w");

   uint nlength = (mesh->mesh_memory.get_memory_size(mesh->nlft) >= mesh->ncells_ghost) ? mesh->ncells_ghost : mesh->ncells;

fprintf(fp,"%d:\tindex\tglobal\ti\tj\tlev\tnlft\tnrht\tnbot\tntop\tH\tHhex\tU\tUhex\tV\tVhex\n",mesh->mype);
   for (uint ic=0; ic < nlength; ic++) {
      fprintf(fp,"%d:\t%d\t%d\t%d\t%d\t%d\t%d\t%d\t%d\t%d", mesh->mype,ic,
                ic+mesh->noffset,mesh->i[ic], mesh->j[ic], mesh->level[ic], mesh->nlft[ic],
                mesh->nrht[ic], mesh->nbot[ic], mesh->ntop[ic]);
 
      fprintf(fp,"\t%lf\t", H[ic]);
 
      doubleToHex(fp, H[ic]);
 
      fprintf(fp,"\t%lf\t",U[ic]);
 
      doubleToHex(fp, U[ic]);
 
      fprintf(fp,"\t%lf\t",V[ic]);
 
      doubleToHex(fp, V[ic]);
 
      fprintf(fp,"\n");

	}
   fclose(fp);
}

size_t State::get_checkpoint_size(void)
{
#ifdef FULL_PRECISION
   size_t nsize = mesh->ncells*3*sizeof(double);
#else
   size_t nsize = mesh->ncells*3*sizeof(float);
#endif
   nsize += num_int_vals*sizeof(int);
   nsize += mesh->get_checkpoint_size();
   return(nsize);
}

void State::store_checkpoint(Crux *crux)
{
   // Store mesh data first
   mesh->store_checkpoint(crux);

   size_t save_size = state_memory.get_memory_capacity(H);
   state_memory.set_restart_length(H,mesh->ncells);
   state_memory.set_restart_length(U,mesh->ncells);
   state_memory.set_restart_length(V,mesh->ncells);
   crux->store_MallocPlus(state_memory);
   state_memory.set_restart_length(H,save_size);
   state_memory.set_restart_length(U,save_size);
   state_memory.set_restart_length(V,save_size);
}

void State::restore_checkpoint(Crux *crux)
{
   // Restore mesh data first
   mesh->restore_checkpoint(crux);

   // Clear memory for restoring data into
   int_vals[0] = 0;

   // allocate is a state method

   state_memory.memory_delete(H);
   state_memory.memory_delete(U);
   state_memory.memory_delete(V);
   allocate(mesh->ncells);
   memory_reset_ptrs();

   // Restore memory database
   crux->restore_MallocPlus(state_memory);

   // Check version number
   if (int_vals[ 0] != CRUX_STATE_VERSION) {
      printf("CRUX version mismatch for state data, version on file is %d, version in code is %d\n",
         int_vals[0], CRUX_STATE_VERSION);
      exit(0);
   }

#ifdef DEBUG_RESTORE_VALS
   if (DEBUG_RESTORE_VALS) {
      printf("\n");
      printf("       === Restored state cpu timers ===\n");
      for (int i = 0; i < STATE_TIMER_SIZE; i++){
         printf("       %-30s %lg\n",state_timer_descriptor[i], cpu_timers[i]);
      }
      printf("       === Restored state cpu timers ===\n");
      printf("\n");
   }
#endif

#ifdef DEBUG_RESTORED_VALS
   if (DEBUG_RESTORED_VALS) {
      printf("\n");
      printf("       === Restored state gpu timers ===\n");
      for (int i = 0; i < STATE_TIMER_SIZE; i++){
         printf("       %-30s %lld\n",state_timer_descriptor[i], gpu_timers[i]);
      }
      printf("       === Restored state gpu_timers ===\n");
      printf("\n");
   }
#endif

   memory_reset_ptrs();
}

// Added overloaded print to get mesh information to print in each cycle
// Brian Atkinson (5-29-14)
void State::print(int iteration, double simTime, double initial_mass, double iteration_mass, double mass_diff_percentage)
{  //printf("size is %lu %lu %lu %lu %lu\n",index.size(), i.size(), level.size(), nlft.size(), x.size());

      char filename[40];
      sprintf(filename,"iteration%d",iteration);
      mesh->fp=fopen(filename,"w");

      if(iteration_mass == 0.0){
         fprintf(mesh->fp,"Iteration = %d\t\tSimuation Time = %lf\n", iteration, simTime);
         fprintf(mesh->fp,"mesh->ncells = %lu\t\tmesh->ncells_ghost = %lu\n", mesh->ncells, mesh->ncells_ghost);
         fprintf(mesh->fp,"Initial Mass: %14.12lg\t\tSimulation Time: %lf\n", initial_mass, simTime);
      }
      else{
         double mass_diff = iteration_mass - initial_mass;
         fprintf(mesh->fp,"Iteration = %d\t\tSimuation Time = %lf\n", iteration, simTime);
         fprintf(mesh->fp,"mesh->ncells = %lu\t\tmesh->ncells_ghost = %lu\n", mesh->ncells, mesh->ncells_ghost);
         fprintf(mesh->fp,"Initial Mass: %14.12lg\t\tIteration Mass: %14.12lg\n", initial_mass, iteration_mass);
         fprintf(mesh->fp,"Mass Difference: %12.6lg\t\tMass Difference Percentage: %12.6lg%%\n", mass_diff, mass_diff_percentage);
      }

   if (mesh->mesh_memory.get_memory_size(mesh->nlft) >= mesh->ncells_ghost){
      fprintf(mesh->fp,"%d:   index global  i     j     lev   nlft  nrht  nbot  ntop \n",mesh->mype);
      for (uint ic=0; ic<mesh->ncells; ic++) {
         fprintf(mesh->fp,"%d: %6d  %6d %4d  %4d   %4d  %4d  %4d  %4d  %4d \n", mesh->mype,ic, ic+mesh->noffset,mesh->i[ic], mesh->j[ic], mesh->level[ic], mesh->nlft[ic], mesh->nrht[ic], mesh->nbot[ic], mesh->ntop[ic]);
      }
      for (uint ic=mesh->ncells; ic<mesh->ncells_ghost; ic++) {
         fprintf(mesh->fp,"%d: %6d  %6d %4d  %4d   %4d  %4d  %4d  %4d  %4d \n", mesh->mype,ic, ic+mesh->noffset,mesh->i[ic], mesh->j[ic], mesh->level[ic], mesh->nlft[ic], mesh->nrht[ic], mesh->nbot[ic], mesh->ntop[ic]);
      }
   } else {
      fprintf(mesh->fp,"%d:  index     H        U         V      i     j     lev\n",mesh->mype);
      for (uint ic=0; ic<mesh->ncells_ghost; ic++) {
         fprintf(mesh->fp,"%d: %6d %lf %lf %lf %4d  %4d   %4d  \n", mesh->mype,ic, H[ic], U[ic], V[ic], mesh->i[ic], mesh->j[ic], mesh->level[ic]);
      }
   }
}

void State::print_local(int ncycle)
{  //printf("size is %lu %lu %lu %lu %lu\n",index.size(), i.size(), level.size(), nlft.size(), x.size());

   if (mesh->fp == NULL) {
      char filename[10];
      sprintf(filename,"out%1d",mesh->mype);
      mesh->fp=fopen(filename,"w");
   }

   fprintf(mesh->fp,"DEBUG in print_local ncycle is %d\n",ncycle);
   if (mesh->nlft != NULL){
      fprintf(mesh->fp,"%d:  index     H        U         V      i     j     lev   nlft   nrht   nbot   ntop\n",mesh->mype);
      uint state_size = state_memory.get_memory_size(H);
      for (uint ic=0; ic<mesh->ncells_ghost; ic++) {
         if (ic >= state_size){
            fprintf(mesh->fp,"%d: %6d                              %4d  %4d   %4d  %4d  %4d  %4d  %4d\n", mesh->mype,ic, mesh->i[ic], mesh->j[ic], mesh->level[ic], mesh->nlft[ic], mesh->nrht[ic], mesh->nbot[ic], mesh->ntop[ic]);
         } else {
            fprintf(mesh->fp,"%d: %6d %lf %lf %lf %4d  %4d   %4d  %4d  %4d  %4d  %4d\n", mesh->mype,ic, H[ic], U[ic], V[ic], mesh->i[ic], mesh->j[ic], mesh->level[ic], mesh->nlft[ic], mesh->nrht[ic], mesh->nbot[ic], mesh->ntop[ic]);
         }
      }
   } else {
      fprintf(mesh->fp,"%d:  index     H        U         V      i     j     lev\n",mesh->mype);
      for (uint ic=0; ic<mesh->ncells_ghost; ic++) {
         fprintf(mesh->fp,"%d: %6d %lf %lf %lf %4d  %4d   %4d\n", mesh->mype,ic, H[ic], U[ic], V[ic], mesh->i[ic], mesh->j[ic], mesh->level[ic]);
      }
   }
}

void State::print_failure_log(int iteration, double simTime, double initial_mass, double iteration_mass, double mass_diff_percentage, bool got_nan){
   char filename[] = {"failure.log"};
   mesh->fp=fopen(filename,"w");

   double mass_diff = iteration_mass - initial_mass;
   if(got_nan){
      fprintf(mesh->fp,"Failed because of nan for H_sum was equal to NAN\n");
   }
   else{
      fprintf(mesh->fp,"Failed because mass difference is outside of accepted percentage\n");
   }
   fprintf(mesh->fp,"Iteration = %d\t\tSimuation Time = %lf\n", iteration, simTime);
   fprintf(mesh->fp,"mesh->ncells = %lu\t\tmesh->ncells_ghost = %lu\n", mesh->ncells, mesh->ncells_ghost);
   fprintf(mesh->fp,"Initial Mass: %14.12lg\t\tIteration Mass: %14.12lg\n", initial_mass, iteration_mass);
   fprintf(mesh->fp,"Mass Difference: %12.6lg\t\tMass Difference Percentage: %12.6lg%%\n", mass_diff, mass_diff_percentage);

   if (mesh->mesh_memory.get_memory_size(mesh->nlft) >= mesh->ncells_ghost){
      fprintf(mesh->fp,"%d:   index global  i     j     lev   nlft  nrht  nbot  ntop \n",mesh->mype);
      for (uint ic=0; ic<mesh->ncells; ic++) {
         fprintf(mesh->fp,"%d: %6d  %6d %4d  %4d   %4d  %4d  %4d  %4d  %4d \n", mesh->mype,ic, ic+mesh->noffset,mesh->i[ic], mesh->j[ic], mesh->level[ic], mesh->nlft[ic], mesh->nrht[ic], mesh->nbot[ic], mesh->ntop[ic]);
      }
      for (uint ic=mesh->ncells; ic<mesh->ncells_ghost; ic++) {
         fprintf(mesh->fp,"%d: %6d  %6d %4d  %4d   %4d  %4d  %4d  %4d  %4d \n", mesh->mype,ic, ic+mesh->noffset,mesh->i[ic], mesh->j[ic], mesh->level[ic], mesh->nlft[ic], mesh->nrht[ic], mesh->nbot[ic], mesh->ntop[ic]);
      }
   } else {
      fprintf(mesh->fp,"%d:  index     H        U         V      i     j     lev\n",mesh->mype);
      for (uint ic=0; ic<mesh->ncells_ghost; ic++) {
         fprintf(mesh->fp,"%d: %6d %lf %lf %lf %4d  %4d   %4d  \n", mesh->mype,ic, H[ic], U[ic], V[ic], mesh->i[ic], mesh->j[ic], mesh->level[ic]);
      }
   }
}

void State::print_rollback_log(int iteration, double simTime, double initial_mass, double iteration_mass, double mass_diff_percentage, int backup_attempt, int num_of_attempts, int error_status){
   char filename[40];
   sprintf(filename, "rollback%d.log",backup_attempt);
   mesh->fp=fopen(filename,"w");

   double mass_diff = iteration_mass - initial_mass;
   if(error_status == STATUS_NAN){
      fprintf(mesh->fp,"Rolling back because of nan for H_sum was equal to NAN\n");
   }
   else{
      fprintf(mesh->fp,"Rolling back because mass difference is outside of accepted percentage\n");
   }
   fprintf(mesh->fp,"Rollback attempt %d of %d ---> Number of attempts left:%d\n", backup_attempt, num_of_attempts, num_of_attempts - backup_attempt);
   fprintf(mesh->fp,"Iteration = %d\t\tSimuation Time = %lf\n", iteration, simTime);
   fprintf(mesh->fp,"mesh->ncells = %lu\t\tmesh->ncells_ghost = %lu\n", mesh->ncells, mesh->ncells_ghost);
   fprintf(mesh->fp,"Initial Mass: %14.12lg\t\tIteration Mass: %14.12lg\n", initial_mass, iteration_mass);
   fprintf(mesh->fp,"Mass Difference: %12.6lg\t\tMass Difference Percentage: %12.6lg%%\n", mass_diff, mass_diff_percentage);

   if (mesh->mesh_memory.get_memory_size(mesh->nlft) >= mesh->ncells_ghost){
      fprintf(mesh->fp,"%d:   index global  i     j     lev   nlft  nrht  nbot  ntop \n",mesh->mype);
      for (uint ic=0; ic<mesh->ncells; ic++) {
         fprintf(mesh->fp,"%d: %6d  %6d %4d  %4d   %4d  %4d  %4d  %4d  %4d \n", mesh->mype,ic, ic+mesh->noffset,mesh->i[ic], mesh->j[ic], mesh->level[ic], mesh->nlft[ic], mesh->nrht[ic], mesh->nbot[ic], mesh->ntop[ic]);
      }
      for (uint ic=mesh->ncells; ic<mesh->ncells_ghost; ic++) {
         fprintf(mesh->fp,"%d: %6d  %6d %4d  %4d   %4d  %4d  %4d  %4d  %4d \n", mesh->mype,ic, ic+mesh->noffset,mesh->i[ic], mesh->j[ic], mesh->level[ic], mesh->nlft[ic], mesh->nrht[ic], mesh->nbot[ic], mesh->ntop[ic]);
      }
   } else {
      fprintf(mesh->fp,"%d:  index     H        U         V      i     j     lev\n",mesh->mype);
      for (uint ic=0; ic<mesh->ncells_ghost; ic++) {
         fprintf(mesh->fp,"%d: %6d %lf %lf %lf %4d  %4d   %4d  \n", mesh->mype,ic, H[ic], U[ic], V[ic], mesh->i[ic], mesh->j[ic], mesh->level[ic]);
      }
   }
}

Mesh_CLAMR::Mesh_CLAMR(int nx, int ny, int levmx_in, int ndim_in, double deltax_in, double deltay_in, int boundary, int parallel_in, int do_gpu_calc) : Mesh(nx,ny,levmx_in,ndim_in,deltax_in,deltay_in,boundary,parallel_in,do_gpu_calc){};


void Mesh_CLAMR::interpolate(int scheme, int index, int cell_lower, int cell_upper, double deltaT, MallocPlus &state_memory){
   switch(scheme){
      case 0: // fine cell, x-direction, right cell more refined
      interpolate_fine_x(0,index,cell_lower,cell_upper,deltaT,state_memory);
      break;

      case 1: // fine cell, x-direction, left cell more refined
      interpolate_fine_x(1,index+2,cell_lower,cell_upper,deltaT,state_memory);
      break;

      case 2: // fine cell, y-direction, top cell more refined
      interpolate_fine_y(0,index,cell_lower,cell_upper,deltaT,state_memory);
      break;

      case 3: // fine cell, y-direction, bottom cell more refined
      interpolate_fine_y(1,index+2,cell_lower,cell_upper,deltaT,state_memory);
      break;

      case 4: // course cell, x-direction, right cell more refined
      interpolate_course_x(0,index+2,cell_lower,cell_upper,deltaT,state_memory);
      break;

      case 5: // course cell, x-direction, left cell more refined
      interpolate_course_x(1,index,cell_lower,cell_upper,deltaT,state_memory);
      break;

      case 6: // course cell, y-direction, top cell more refined
      interpolate_course_y(0,index+2,cell_lower,cell_upper,deltaT,state_memory);
      break;

      case 7: // course cell, y-direction, bottom cell more refined
      interpolate_course_x(1,index,cell_lower,cell_upper,deltaT,state_memory);
      break;
   } 
}

void Mesh_CLAMR::interpolate_fine_x(int scheme, int index, int cell_lower, int cell_upper, double deltaT, MallocPlus &state_memory){
   state_t* H = (state_t *)state_memory.get_memory_ptr("H");
   state_t* U = (state_t *)state_memory.get_memory_ptr("U");
   state_t* V = (state_t *)state_memory.get_memory_ptr("V");

   real_t dx_lower = lev_deltax[level[cell_lower]];
   real_t dx_upper = lev_deltax[level[cell_upper]];
   real_t FA_lower = dx_lower;
   real_t FA_upper = dx_upper;
   real_t FA_lolim = FA_lower*min(ONE, FA_upper/FA_lower);
   real_t FA_uplim = FA_upper*min(ONE, FA_lower/FA_upper);
   real_t CV_lower = SQ(dx_lower);
   real_t CV_upper = SQ(dx_upper);
   real_t CV_lolim = CV_lower*min(HALF, CV_upper/CV_lower);
   real_t CV_uplim = CV_upper*min(HALF, CV_lower/CV_upper);
   real_t g     = 9.80;   // gravitational constant
   real_t ghalf = 0.5*g;
   switch(scheme){
      case 0: // H,U,V interpolation, right cell more refined
         H[index] = (2*(dx_lower*H[cell_upper]+dx_lower*H[cell_lower])/(dx_lower+dx_upper)
         - deltaT*((FA_uplim*HXFLUX(cell_upper)-FA_lolim*HXFLUX(cell_lower))/(CV_uplim+CV_lolim)
         - (HXFLUX(cell_upper)-HXFLUX(cell_lower))/dx_upper) - H[cell_upper]);
         U[index] = (2*(dx_lower*U[cell_upper]+dx_lower*U[cell_lower])/(dx_lower+dx_upper)
         - deltaT*((FA_uplim*UXFLUX(cell_upper)-FA_lolim*UXFLUX(cell_lower))/(CV_uplim+CV_lolim)
         - (UXFLUX(cell_upper)-UXFLUX(cell_lower))/dx_upper) - U[cell_upper]);
         V[index] = (2*(dx_lower*V[cell_upper]+dx_lower*V[cell_lower])/(dx_lower+dx_upper)
         - deltaT*((FA_uplim*UVFLUX(cell_upper)-FA_lolim*UVFLUX(cell_lower))/(CV_uplim+CV_lolim)
         - (UVFLUX(cell_upper)-UVFLUX(cell_lower))/dx_upper) - V[cell_upper]);
      break;   
      case 1: // H,U,V interpolation, left cell more refined
         H[index] = (2*(dx_lower*H[cell_upper]+dx_lower*H[cell_lower])/(dx_lower+dx_upper)
         - deltaT*((FA_uplim*HXFLUX(cell_upper)-FA_lolim*HXFLUX(cell_lower))/(CV_uplim+CV_lolim)
         - (HXFLUX(cell_upper)-HXFLUX(cell_lower))/dx_upper) - H[cell_lower]);
         U[index] = (2*(dx_lower*U[cell_upper]+dx_lower*U[cell_lower])/(dx_lower+dx_upper)
         - deltaT*((FA_uplim*UXFLUX(cell_upper)-FA_lolim*UXFLUX(cell_lower))/(CV_uplim+CV_lolim)
         - (UXFLUX(cell_upper)-UXFLUX(cell_lower))/dx_upper) - U[cell_lower]);
         V[index] = (2*(dx_lower*V[cell_upper]+dx_lower*V[cell_lower])/(dx_lower+dx_upper)
         - deltaT*((FA_uplim*UVFLUX(cell_upper)-FA_lolim*UVFLUX(cell_lower))/(CV_uplim+CV_lolim)
         - (UVFLUX(cell_upper)-UVFLUX(cell_lower))/dx_upper) - V[cell_lower]);
      break;
   }
   //printf("DEBUG MESH: ID %d) LOWER:  %d, UPPER: %d, POS: lft\n",index,cell_lower,cell_upper);
   //printf("            H:  %f, U: %f, V: %f\n", H[index],U[index],V[index]);
} 

void Mesh_CLAMR::interpolate_fine_y(int scheme, int index, int cell_lower, int cell_upper, double deltaT, MallocPlus &state_memory){
   state_t* H = (state_t *)state_memory.get_memory_ptr("H");
   state_t* U = (state_t *)state_memory.get_memory_ptr("U");
   state_t* V = (state_t *)state_memory.get_memory_ptr("V");

   real_t dy_lower = lev_deltay[level[cell_lower]];
   real_t dy_upper = lev_deltay[level[cell_upper]];
   real_t FA_lower = dy_lower;
   real_t FA_upper = dy_upper;
   real_t FA_lolim = FA_lower*min(ONE, FA_upper/FA_lower);
   real_t FA_uplim = FA_upper*min(ONE, FA_lower/FA_upper);
   real_t CV_lower = SQ(dy_lower);
   real_t CV_upper = SQ(dy_upper);
   real_t CV_lolim = CV_lower*min(HALF, CV_upper/CV_lower);
   real_t CV_uplim = CV_upper*min(HALF, CV_lower/CV_upper);
   real_t g     = 9.80;   // gravitational constant
   real_t ghalf = 0.5*g;
   switch(scheme){
      case 0: // H,U,V interpolation, top cell more refined 
         H[index] = (2*(dy_lower*H[cell_upper]+dy_lower*H[cell_lower])/(dy_lower+dy_upper)
         - deltaT*((FA_uplim*HYFLUX(cell_upper)-FA_lolim*HYFLUX(cell_lower))/(CV_uplim+CV_lolim)
         - (HYFLUX(cell_upper)-HYFLUX(cell_lower))/dy_upper) - H[cell_upper]);
         U[index] = (2*(dy_lower*U[cell_upper]+dy_lower*U[cell_lower])/(dy_lower+dy_upper)
         - deltaT*((FA_uplim*UVFLUX(cell_upper)-FA_lolim*UVFLUX(cell_lower))/(CV_uplim+CV_lolim)
         - (UVFLUX(cell_upper)-UVFLUX(cell_lower))/dy_upper) - U[cell_upper]);
         V[index] = (2*(dy_lower*V[cell_upper]+dy_lower*V[cell_lower])/(dy_lower+dy_upper)
         - deltaT*((FA_uplim*VYFLUX(cell_upper)-FA_lolim*VYFLUX(cell_lower))/(CV_uplim+CV_lolim)
         - (VYFLUX(cell_upper)-VYFLUX(cell_lower))/dy_upper) - V[cell_upper]);
   break;
      case 1: // H,U,V interpolation, bottom cell more refined
         H[index] = (2*(dy_lower*H[cell_upper]+dy_lower*H[cell_lower])/(dy_lower+dy_upper)
         - deltaT*((FA_uplim*HYFLUX(cell_upper)-FA_lolim*HYFLUX(cell_lower))/(CV_uplim+CV_lolim)
         - (HYFLUX(cell_upper)-HYFLUX(cell_lower))/dy_upper) - H[cell_lower]);
         U[index] = (2*(dy_lower*U[cell_upper]+dy_lower*U[cell_lower])/(dy_lower+dy_upper)
         - deltaT*((FA_uplim*UVFLUX(cell_upper)-FA_lolim*UVFLUX(cell_lower))/(CV_uplim+CV_lolim)
         - (UVFLUX(cell_upper)-UVFLUX(cell_lower))/dy_upper) - U[cell_lower]);
         V[index] = (2*(dy_lower*V[cell_upper]+dy_lower*V[cell_lower])/(dy_lower+dy_upper)
         - deltaT*((FA_uplim*VYFLUX(cell_upper)-FA_lolim*VYFLUX(cell_lower))/(CV_uplim+CV_lolim)
         - (VYFLUX(cell_upper)-VYFLUX(cell_lower))/dy_upper) - V[cell_lower]);
      break;
   }
   //printf("DEBUG MESH: ID %d) LOWER:  %d, UPPER: %d, POS: lft\n",index,cell_lower,cell_upper);
   //printf("            H:  %f, U: %f, V: %f\n", H[index],U[index],V[index]);
}

void Mesh_CLAMR::interpolate_course_x(int scheme, int index, int cell_lower, int cell_upper, double deltaT, MallocPlus &state_memory){
   state_t* H = (state_t *)state_memory.get_memory_ptr("H");
   state_t* U = (state_t *)state_memory.get_memory_ptr("U");
   state_t* V = (state_t *)state_memory.get_memory_ptr("V");

   int cell_course, cell_bot, cell_top;
   real_t dx_lower_bot, dx_lower_top, dx_upper_bot, dx_upper_top;
   real_t FA_lower_bot, FA_lower_top, FA_upper_bot, FA_upper_top;
   real_t FA_lolim_bot, FA_lolim_top, FA_uplim_bot, FA_uplim_top;
   real_t CV_lower_bot, CV_lower_top, CV_upper_bot, CV_upper_top;
   real_t CV_lolim_bot, CV_lolim_top, CV_uplim_bot, CV_uplim_top;
   real_t Hx_bot, Hx_top, Ux_bot, Ux_top, Vx_bot, Vx_top;
   real_t g = 9.80;   // gravitational constant
   real_t ghalf = 0.5*g;
   switch(scheme){
      case 0: // H,U,V interpolation, right cell more refined  
         cell_course = cell_lower;
         cell_bot = nrht[cell_course];
         cell_top = ntop[cell_bot];
         dx_lower_bot = lev_deltax[level[cell_course]];
         dx_lower_top = lev_deltax[level[cell_course]];
         dx_upper_bot = lev_deltax[level[cell_bot]];
         dx_upper_top = lev_deltax[level[cell_top]];
         FA_lower_bot = dx_lower_bot;
         FA_lower_top = dx_lower_top;
         FA_upper_bot = dx_upper_bot;
         FA_upper_top = dx_upper_top;
         FA_lolim_bot = FA_lower_bot*min(ONE, FA_upper_bot/FA_lower_bot);
         FA_lolim_top = FA_lower_top*min(ONE, FA_upper_top/FA_lower_top);
         FA_uplim_bot = FA_upper_bot*min(ONE, FA_lower_bot/FA_upper_bot);
         FA_uplim_top = FA_upper_top*min(ONE, FA_lower_top/FA_upper_top);
         CV_lower_bot = SQ(dx_lower_bot);
         CV_lower_top = SQ(dx_lower_top);
         CV_upper_bot = SQ(dx_upper_bot);
         CV_upper_top = SQ(dx_upper_top);
         CV_lolim_bot = CV_lower_bot*min(HALF, CV_upper_bot/CV_lower_bot);
         CV_lolim_top = CV_lower_top*min(HALF, CV_upper_top/CV_lower_top);
         CV_uplim_bot = CV_upper_bot*min(HALF, CV_lower_bot/CV_upper_bot);
         CV_uplim_top = CV_upper_top*min(HALF, CV_lower_top/CV_upper_top);

         Hx_bot = (dx_lower_bot*H[cell_bot]+dx_upper_bot*H[cell_course])/(dx_lower_bot+dx_upper_bot) -
                 HALF*deltaT*( (FA_uplim_bot*HXFLUX(cell_bot))-(FA_lolim_bot*HXFLUX(cell_course)) )/
                 (CV_uplim_bot+CV_lolim_bot);
         Hx_top = (dx_lower_top*H[cell_top]+dx_upper_top*H[cell_course])/(dx_lower_top+dx_upper_top) -
                 HALF*deltaT*( (FA_uplim_top*HXFLUX(cell_top))-(FA_lolim_top*HXFLUX(cell_course)) )/
                 (CV_uplim_bot+CV_lolim_bot);
         Ux_bot = (dx_lower_bot*U[cell_bot]+dx_upper_bot*U[cell_course])/(dx_lower_bot+dx_upper_bot) -
                 HALF*deltaT*( (FA_uplim_bot*UXFLUX(cell_bot))-(FA_lolim_bot*UXFLUX(cell_course)) )/
                 (CV_uplim_bot+CV_lolim_bot);
         Ux_top = (dx_lower_top*U[cell_top]+dx_upper_top*U[cell_course])/(dx_lower_top+dx_upper_top) -
                 HALF*deltaT*( (FA_uplim_top*UXFLUX(cell_top))-(FA_lolim_top*UXFLUX(cell_course)) )/
                 (CV_uplim_top+CV_lolim_top);
         Vx_bot = (dx_lower_bot*V[cell_bot]+dx_upper_bot*V[cell_course])/(dx_lower_bot+dx_upper_bot) -
                 HALF*deltaT*( (FA_uplim_bot*UVFLUX(cell_bot))-(FA_lolim_bot*UVFLUX(cell_course)) )/
                 (CV_uplim_bot+CV_lolim_bot);
         Vx_top = (dx_lower_top*V[cell_top]+dx_upper_top*V[cell_course])/(dx_lower_top+dx_upper_top) -
                 HALF*deltaT*( (FA_uplim_top*UVFLUX(cell_top))-(FA_lolim_top*UVFLUX(cell_course)) )/
                 (CV_uplim_top+CV_lolim_top);

         H[index] = (2*(Hx_bot+Hx_top) + deltaT*(HXFLUX(cell_upper)-HXFLUX(cell_lower))/dx_lower_bot - H[cell_upper]);
         U[index] = (2*(Ux_bot+Ux_top) + deltaT*(UXFLUX(cell_upper)-UXFLUX(cell_lower))/dx_lower_bot - U[cell_upper]);
         V[index] = (2*(Vx_bot+Vx_top) + deltaT*(UVFLUX(cell_upper)-UVFLUX(cell_lower))/dx_lower_bot - V[cell_upper]);
      break;
      case 1: // H,U,V interpolation, left cell more refined  
         cell_course = cell_upper;
         cell_bot = nlft[cell_course];
         cell_top = ntop[cell_bot];         
         dx_lower_bot = lev_deltax[level[cell_bot]];
         dx_lower_top = lev_deltax[level[cell_top]];
         dx_upper_bot = lev_deltax[level[cell_course]];
         dx_upper_top = lev_deltax[level[cell_course]];
         FA_lower_bot = dx_lower_bot;
         FA_lower_top = dx_lower_top;
         FA_upper_bot = dx_upper_bot;
         FA_upper_top = dx_upper_top;
         FA_lolim_bot = FA_lower_bot*min(ONE, FA_upper_bot/FA_lower_bot);
         FA_lolim_top = FA_lower_top*min(ONE, FA_upper_top/FA_lower_top);
         FA_uplim_bot = FA_upper_bot*min(ONE, FA_lower_bot/FA_upper_bot);
         FA_uplim_top = FA_upper_top*min(ONE, FA_lower_top/FA_upper_top);
         CV_lower_bot = SQ(dx_lower_bot);
         CV_lower_top = SQ(dx_lower_top);
         CV_upper_bot = SQ(dx_upper_bot);
         CV_upper_top = SQ(dx_upper_top);
         CV_lolim_bot = CV_lower_bot*min(HALF, CV_upper_bot/CV_lower_bot);
         CV_lolim_top = CV_lower_top*min(HALF, CV_upper_top/CV_lower_top);
         CV_uplim_bot = CV_upper_bot*min(HALF, CV_lower_bot/CV_upper_bot);
         CV_uplim_top = CV_upper_top*min(HALF, CV_lower_top/CV_upper_top);

         Hx_bot = (dx_lower_bot*H[cell_bot]+dx_upper_bot*H[cell_course])/(dx_lower_bot+dx_upper_bot) -
                 HALF*deltaT*( (FA_uplim_bot*HXFLUX(cell_bot))-(FA_lolim_bot*HXFLUX(cell_course)) )/
                 (CV_uplim_bot+CV_lolim_bot);
         Hx_top = (dx_lower_top*H[cell_top]+dx_upper_top*H[cell_course])/(dx_lower_top+dx_upper_top) -
                 HALF*deltaT*( (FA_uplim_top*HXFLUX(cell_top))-(FA_lolim_top*HXFLUX(cell_course)) )/
                 (CV_uplim_bot+CV_lolim_bot);
         Ux_bot = (dx_lower_bot*U[cell_bot]+dx_upper_bot*U[cell_course])/(dx_lower_bot+dx_upper_bot) -
                 HALF*deltaT*( (FA_uplim_bot*UXFLUX(cell_bot))-(FA_lolim_bot*UXFLUX(cell_course)) )/
                 (CV_uplim_bot+CV_lolim_bot);
         Ux_top = (dx_lower_top*U[cell_top]+dx_upper_top*U[cell_course])/(dx_lower_top+dx_upper_top) -
                 HALF*deltaT*( (FA_uplim_top*UXFLUX(cell_top))-(FA_lolim_top*UXFLUX(cell_course)) )/
                 (CV_uplim_top+CV_lolim_top);
         Vx_bot = (dx_lower_bot*V[cell_bot]+dx_upper_bot*V[cell_course])/(dx_lower_bot+dx_upper_bot) -
                 HALF*deltaT*( (FA_uplim_bot*UVFLUX(cell_bot))-(FA_lolim_bot*UVFLUX(cell_course)) )/
                 (CV_uplim_bot+CV_lolim_bot);
         Vx_top = (dx_lower_top*V[cell_top]+dx_upper_top*V[cell_course])/(dx_lower_top+dx_upper_top) -
                 HALF*deltaT*( (FA_uplim_top*UVFLUX(cell_top))-(FA_lolim_top*UVFLUX(cell_course)) )/
                 (CV_uplim_top+CV_lolim_top);

         H[index] = (2*(Hx_bot+Hx_top) + deltaT*(HXFLUX(cell_upper)-HXFLUX(cell_lower))/dx_upper_bot - H[cell_lower]);
         U[index] = (2*(Ux_bot+Ux_top) + deltaT*(UXFLUX(cell_upper)-UXFLUX(cell_lower))/dx_upper_bot - U[cell_lower]);
         V[index] = (2*(Vx_bot+Vx_top) + deltaT*(UVFLUX(cell_upper)-UVFLUX(cell_lower))/dx_upper_bot - V[cell_lower]);
      break;
   }
   //printf("DEBUG MESH: ID %d) LOWER:  %d, UPPER: %d, POS: lft\n",index,cell_lower,cell_upper);
   //printf("            H:  %f, U: %f, V: %f\n", H[index],U[index],V[index]);
}

void Mesh_CLAMR::interpolate_course_y(int scheme, int index, int cell_lower, int cell_upper , double deltaT, MallocPlus &state_memory){
   state_t* H = (state_t *)state_memory.get_memory_ptr("H");
   state_t* U = (state_t *)state_memory.get_memory_ptr("U");
   state_t* V = (state_t *)state_memory.get_memory_ptr("V");

   int cell_course, cell_left, cell_right;
   real_t dy_lower_left, dy_lower_right, dy_upper_left, dy_upper_right;
   real_t FA_lower_left, FA_lower_right, FA_upper_left, FA_upper_right;
   real_t FA_lolim_left, FA_lolim_right, FA_uplim_left, FA_uplim_right;
   real_t CV_lower_left, CV_lower_right, CV_upper_left, CV_upper_right;
   real_t CV_lolim_left, CV_lolim_right, CV_uplim_left, CV_uplim_right;
   real_t Hy_left, Hy_right, Uy_left, Uy_right, Vy_left, Vy_right;
   real_t g = 9.80;   // gravitational constant
   real_t ghalf = 0.5*g;

   switch(scheme){
      case 0: // H,U,V interpolation, top cell more refined  
         cell_course = cell_lower;
         cell_left = ntop[cell_course];
         cell_right = nrht[cell_left];
         dy_lower_left = lev_deltay[level[cell_course]];
         dy_lower_right = lev_deltay[level[cell_course]];
         dy_upper_left = lev_deltay[level[cell_left]];
         dy_upper_right = lev_deltay[level[cell_right]];
         FA_lower_left = dy_lower_left;
         FA_lower_right = dy_lower_right;
         FA_upper_left = dy_upper_left;
         FA_upper_right = dy_upper_right;
         FA_lolim_left = FA_lower_left*min(ONE, FA_upper_left/FA_lower_left);
         FA_lolim_right = FA_lower_right*min(ONE, FA_upper_right/FA_lower_right);
         FA_uplim_left = FA_upper_left*min(ONE, FA_lower_left/FA_upper_left);
         FA_uplim_right = FA_upper_right*min(ONE, FA_lower_right/FA_upper_right);
         CV_lower_left = SQ(dy_lower_left);
         CV_lower_right = SQ(dy_lower_right);
         CV_upper_left = SQ(dy_upper_left);
         CV_upper_right = SQ(dy_upper_right);
         CV_lolim_left = CV_lower_left*min(HALF, CV_upper_left/CV_lower_left);
         CV_lolim_right = CV_lower_right*min(HALF, CV_upper_right/CV_lower_right);
         CV_uplim_left = CV_upper_left*min(HALF, CV_lower_left/CV_upper_left);
         CV_uplim_right = CV_upper_right*min(HALF, CV_lower_right/CV_upper_right);

         Hy_left = (dy_lower_left*H[cell_left]+dy_upper_left*H[cell_course])/(dy_lower_left+dy_upper_left) -
                 HALF*deltaT*( (FA_uplim_left*HYFLUX(cell_left))-(FA_lolim_left*HYFLUX(cell_course)) )/
                 (CV_uplim_left+CV_lolim_left);
         Hy_right = (dy_lower_right*H[cell_right]+dy_upper_right*H[cell_course])/(dy_lower_right+dy_upper_right) -
                 HALF*deltaT*( (FA_uplim_right*HYFLUX(cell_right))-(FA_lolim_right*HYFLUX(cell_course)) )/
                 (CV_uplim_left+CV_lolim_left);
         Uy_left = (dy_lower_left*U[cell_left]+dy_upper_left*U[cell_course])/(dy_lower_left+dy_upper_left) -
                 HALF*deltaT*( (FA_uplim_left*UVFLUX(cell_left))-(FA_lolim_left*UVFLUX(cell_course)) )/
                 (CV_uplim_left+CV_lolim_left);
         Uy_right = (dy_lower_right*U[cell_right]+dy_upper_right*U[cell_course])/(dy_lower_right+dy_upper_right) -
                 HALF*deltaT*( (FA_uplim_right*UVFLUX(cell_right))-(FA_lolim_right*UVFLUX(cell_course)) )/
                 (CV_uplim_right+CV_lolim_right);
         Vy_left = (dy_lower_left*V[cell_left]+dy_upper_left*V[cell_course])/(dy_lower_left+dy_upper_left) -
                 HALF*deltaT*( (FA_uplim_left*VYFLUX(cell_left))-(FA_lolim_left*VYFLUX(cell_course)) )/
                 (CV_uplim_left+CV_lolim_left);
         Vy_right = (dy_lower_right*V[cell_right]+dy_upper_right*V[cell_course])/(dy_lower_right+dy_upper_right) -
                 HALF*deltaT*( (FA_uplim_right*VYFLUX(cell_right))-(FA_lolim_right*VYFLUX(cell_course)) )/
                 (CV_uplim_right+CV_lolim_right);

         H[index] = (2*(Hy_left+Hy_right) + deltaT*(HYFLUX(cell_upper)-HYFLUX(cell_lower))/dy_lower_left - H[cell_upper]);
         U[index] = (2*(Uy_left+Uy_right) + deltaT*(UVFLUX(cell_upper)-UVFLUX(cell_lower))/dy_lower_left - U[cell_upper]);
         V[index] = (2*(Vy_left+Vy_right) + deltaT*(VYFLUX(cell_upper)-VYFLUX(cell_lower))/dy_lower_left - V[cell_upper]);
      break;
      case 1: // H,U,V interpolation, bottom cell more refined  
         cell_course = cell_upper;
         cell_left = nbot[cell_course];
         cell_right = nrht[cell_left];
         dy_lower_left = lev_deltay[level[cell_left]];
         dy_lower_right = lev_deltay[level[cell_right]];
         dy_upper_left = lev_deltay[level[cell_course]];
         dy_upper_right = lev_deltay[level[cell_course]];
         FA_lower_left = dy_lower_left;
         FA_lower_right = dy_lower_right;
         FA_upper_left = dy_upper_left;
         FA_upper_right = dy_upper_right;
         FA_lolim_left = FA_lower_left*min(ONE, FA_upper_left/FA_lower_left);
         FA_lolim_right = FA_lower_right*min(ONE, FA_upper_right/FA_lower_right);
         FA_uplim_left = FA_upper_left*min(ONE, FA_lower_left/FA_upper_left);
         FA_uplim_right = FA_upper_right*min(ONE, FA_lower_right/FA_upper_right);
         CV_lower_left = SQ(dy_lower_left);
         CV_lower_right = SQ(dy_lower_right);
         CV_upper_left = SQ(dy_upper_left);
         CV_upper_right = SQ(dy_upper_right);
         CV_lolim_left = CV_lower_left*min(HALF, CV_upper_left/CV_lower_left);
         CV_lolim_right = CV_lower_right*min(HALF, CV_upper_right/CV_lower_right);
         CV_uplim_left = CV_upper_left*min(HALF, CV_lower_left/CV_upper_left);
         CV_uplim_right = CV_upper_right*min(HALF, CV_lower_right/CV_upper_right);

         Hy_left = (dy_lower_left*H[cell_left]+dy_upper_left*H[cell_course])/(dy_lower_left+dy_upper_left) -
                 HALF*deltaT*( (FA_uplim_left*HYFLUX(cell_left))-(FA_lolim_left*HYFLUX(cell_course)) )/
                 (CV_uplim_left+CV_lolim_left);
         Hy_right = (dy_lower_right*H[cell_right]+dy_upper_right*H[cell_course])/(dy_lower_right+dy_upper_right) -
                 HALF*deltaT*( (FA_uplim_right*HYFLUX(cell_right))-(FA_lolim_right*HYFLUX(cell_course)) )/
                 (CV_uplim_left+CV_lolim_left);
         Uy_left = (dy_lower_left*U[cell_left]+dy_upper_left*U[cell_course])/(dy_lower_left+dy_upper_left) -
                 HALF*deltaT*( (FA_uplim_left*UVFLUX(cell_left))-(FA_lolim_left*UVFLUX(cell_course)) )/
                 (CV_uplim_left+CV_lolim_left);
         Uy_right = (dy_lower_right*U[cell_right]+dy_upper_right*U[cell_course])/(dy_lower_right+dy_upper_right) -
                 HALF*deltaT*( (FA_uplim_right*UVFLUX(cell_right))-(FA_lolim_right*UVFLUX(cell_course)) )/
                 (CV_uplim_right+CV_lolim_right);
         Vy_left = (dy_lower_left*V[cell_left]+dy_upper_left*V[cell_course])/(dy_lower_left+dy_upper_left) -
                 HALF*deltaT*( (FA_uplim_left*VYFLUX(cell_left))-(FA_lolim_left*VYFLUX(cell_course)) )/
                 (CV_uplim_left+CV_lolim_left);
         Vy_right = (dy_lower_right*V[cell_right]+dy_upper_right*V[cell_course])/(dy_lower_right+dy_upper_right) -
                 HALF*deltaT*( (FA_uplim_right*VYFLUX(cell_right))-(FA_lolim_right*VYFLUX(cell_course)) )/
                 (CV_uplim_right+CV_lolim_right);

         H[index] = (2*(Hy_left+Hy_right) + deltaT*(HYFLUX(cell_upper)-HYFLUX(cell_lower))/dy_upper_left - H[cell_lower]);
         U[index] = (2*(Uy_left+Uy_right) + deltaT*(UVFLUX(cell_upper)-UVFLUX(cell_lower))/dy_upper_left - U[cell_lower]);
         V[index] = (2*(Vy_left+Vy_right) + deltaT*(VYFLUX(cell_upper)-VYFLUX(cell_lower))/dy_upper_left - V[cell_lower]);
      break;
   }
   //printf("DEBUG MESH: ID %d) LOWER:  %d, UPPER: %d, POS: lft\n",index,cell_lower,cell_upper);
   //printf("            H:  %f, U: %f, V: %f\n", H[index],U[index],V[index]);
}
