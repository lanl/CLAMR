#include "mesh.h"
#include "l7.h"
#include "mpi.h"
#include "genmalloc.h"

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

void set_mesh_values(int *i, int *j, int *level);
void set_xface_check(int *xface_i_check, int *xface_j_check, int *xface_level_check);
void set_yface_check(int *yface_i_check, int *yface_j_check, int *yface_level_check);
void set_map_xface2cell_check(int *map_xface2cell_lower_check, int *map_xface2cell_upper_check);
void set_map_yface2cell_check(int *map_yface2cell_lower_check, int *map_yface2cell_upper_check);
int **set_xface_flag_check(int jsize, int isize, int lev);
int **set_yface_flag_check(int jsize, int isize, int lev);
int **set_zone_flag_check(int jsize, int isize, int lev);
int **set_zone_cell_check(int jsize, int isize, int lev);
int *set_map_xcell2face_left1_check(int ncells);
int *set_map_xcell2face_left2_check(int ncells);
int *set_map_xcell2face_right1_check(int ncells);
int *set_map_xcell2face_right2_check(int ncells);
int *set_map_ycell2face_bot1_check(int ncells);
int *set_map_ycell2face_bot2_check(int ncells);
int *set_map_ycell2face_top1_check(int ncells);
int *set_map_ycell2face_top2_check(int ncells);

int main (int argc, char **argv)
{

   int nx = 4;
   int ny = 4;
   int levmx = 2;
   int ndim = 2;
   int boundary = 1;
   int parallel_in = 0;
   int do_gpu_calc = 0;
   double circ_radius = 6.0;
   circ_radius *= (double)nx / 128.0; // scaling radius for problem size

   double deltax_in = 1.0;
   double deltay_in = 1.0;
   mesh = new Mesh(nx, ny, levmx, ndim, deltax_in, deltay_in, boundary, parallel_in, do_gpu_calc);

   mesh->init(nx, ny, circ_radius, initial_order, do_gpu_calc);

   mesh->ncells = 112;
   mesh->i = (int *)mesh->mesh_memory.memory_realloc(mesh->ncells, mesh->i);
   mesh->j = (int *)mesh->mesh_memory.memory_realloc(mesh->ncells, mesh->j);
   mesh->level = (int *)mesh->mesh_memory.memory_realloc(mesh->ncells, mesh->level);
   
   mesh->index.resize(mesh->ncells);
   for (uint ic = 0; ic < mesh->ncells; ic++){
      mesh->index[ic] = ic;
   }

   int *i = mesh->i;
   int *j = mesh->j;
   int *level = mesh->level;

   set_mesh_values(i, j, level);

   vector<double> H(mesh->ncells);

   mesh->calc_neighbors(mesh->ncells);
   mesh->calc_celltype(mesh->ncells);
   //mesh->calc_spatial_coordinates(0);

   int nxface_check = 88;
   vector<int>xface_i_check(nxface_check);
   vector<int>xface_j_check(nxface_check);
   vector<int>xface_level_check(nxface_check);

   set_xface_check(&xface_i_check[0], &xface_j_check[0], &xface_level_check[0]);

   int nyface_check = 88;
   vector<int>yface_i_check(nyface_check);
   vector<int>yface_j_check(nyface_check);
   vector<int>yface_level_check(nyface_check);

   set_yface_check(&yface_i_check[0], &yface_j_check[0], &yface_level_check[0]);

   vector<int>map_xface2cell_lower_check(nxface_check);
   vector<int>map_xface2cell_upper_check(nxface_check);

   set_map_xface2cell_check(&map_xface2cell_lower_check[0], &map_xface2cell_upper_check[0]);

   vector<int>map_yface2cell_lower_check(nyface_check);
   vector<int>map_yface2cell_upper_check(nyface_check);

   set_map_yface2cell_check(&map_yface2cell_lower_check[0], &map_yface2cell_upper_check[0]);

   mesh->calc_face_list();

   int icount_err = 0;

   if (mesh->nxface != nxface_check){
      printf("Error -- nxface does not match, nxface %d nxface_check %d\n",
         mesh->nxface, nxface_check);
      icount_err++;
   //} else {
   //   printf("Passed -- nxface check\n");
   }

   if (mesh->nyface != nyface_check){
      printf("Error -- nyface does not match, nyface %d nyface_check %d\n",
         mesh->nyface, nyface_check);
      icount_err++;
   //} else {
   //   printf("Passed -- nyface check\n");
   }

   for (int iface = 0; iface < mesh->nxface; iface++){
      if (mesh->xface_i[iface] != xface_i_check[iface]){
         printf("Error -- xface_i does not match for face %d xface_i %d xface_i_check %d\n",
            iface, mesh->xface_i[iface], xface_i_check[iface]);
         icount_err++;
      }
      if (mesh->xface_j[iface] != xface_j_check[iface]){
         printf("Error -- xface_j does not match for face %d xface_j %d xface_j_check %d\n",
            iface, mesh->xface_j[iface], xface_j_check[iface]);
         icount_err++;
      }
      if (mesh->xface_level[iface] != xface_level_check[iface]){
         printf("Error -- xface_level does not match for face %d xface_level %d xface_level_check %d\n",
            iface, mesh->xface_level[iface], xface_level_check[iface]);
         icount_err++;
      }
   }

   for (int iface = 0; iface < mesh->nyface; iface++){
      if (mesh->yface_i[iface] != yface_i_check[iface]){
         printf("Error -- yface_i does not match for face %d yface_i %d yface_i_check %d\n",
            iface, mesh->yface_i[iface], yface_i_check[iface]);
         icount_err++;
      }
      if (mesh->yface_j[iface] != yface_j_check[iface]){
         printf("Error -- yface_j does not match for face %d yface_j %d yface_j_check %d\n",
            iface, mesh->yface_j[iface], yface_j_check[iface]);
         icount_err++;
      }
      if (mesh->yface_level[iface] != yface_level_check[iface]){
         printf("Error -- yface_level does not match for face %d yface_level %d yface_level_check %d\n",
            iface, mesh->yface_level[iface], yface_level_check[iface]);
         icount_err++;
      }
   }


   //if (mype == 0){
      if (icount_err){
          printf("  Error with calc face\n");
      }
      else{
          printf("  PASSED calc face\n");
      }
   //}

   mesh->calc_face_list_wmap();

   icount_err = 0;

   if (mesh->nxface != nxface_check){
      printf("Error -- nxface does not match, nxface %d nxface_check %d\n",
         mesh->nxface, nxface_check);
      icount_err++;
   //} else {
   //   printf("Passed -- nxface check\n");
   }

   if (mesh->nyface != nyface_check){
      printf("Error -- nyface does not match, nyface %d nyface_check %d\n",
         mesh->nyface, nyface_check);
      icount_err++;
   //} else {
   //   printf("Passed -- nyface check\n");
   }

   for (int iface = 0; iface < mesh->nxface; iface++){
      if (mesh->xface_i[iface] != xface_i_check[iface]){
         printf("Error -- xface_i does not match for face %d xface_i %d xface_i_check %d\n",
            iface, mesh->xface_i[iface], xface_i_check[iface]);
         icount_err++;
      }
      if (mesh->xface_j[iface] != xface_j_check[iface]){
         printf("Error -- xface_j does not match for face %d xface_j %d xface_j_check %d\n",
            iface, mesh->xface_j[iface], xface_j_check[iface]);
         icount_err++;
      }
      if (mesh->xface_level[iface] != xface_level_check[iface]){
         printf("Error -- xface_level does not match for face %d xface_level %d xface_level_check %d\n",
            iface, mesh->xface_level[iface], xface_level_check[iface]);
         icount_err++;
      }
      if (mesh->map_xface2cell_lower[iface] != map_xface2cell_lower_check[iface]){
         printf("Error -- map_xface2cell_lower does not match for face %d map_xface2cell_lower %d map_xface2cell_lower_check %d\n",
            iface, mesh->map_xface2cell_lower[iface], map_xface2cell_lower_check[iface]);
         icount_err++;
      }
      if (mesh->map_xface2cell_upper[iface] != map_xface2cell_upper_check[iface]){
         printf("Error -- map_xface2cell_upper does not match for face %d map_xface2cell_upper %d map_xface2cell_upper_check %d\n",
            iface, mesh->map_xface2cell_upper[iface], map_xface2cell_upper_check[iface]);
         icount_err++;
      }
   }

   for (int iface = 0; iface < mesh->nyface; iface++){
      if (mesh->yface_i[iface] != yface_i_check[iface]){
         printf("Error -- yface_i does not match for face %d yface_i %d yface_i_check %d\n",
            iface, mesh->yface_i[iface], yface_i_check[iface]);
         icount_err++;
      }
      if (mesh->yface_j[iface] != yface_j_check[iface]){
         printf("Error -- yface_j does not match for face %d yface_j %d yface_j_check %d\n",
            iface, mesh->yface_j[iface], yface_j_check[iface]);
         icount_err++;
      }
      if (mesh->yface_level[iface] != yface_level_check[iface]){
         printf("Error -- yface_level does not match for face %d yface_level %d yface_level_check %d\n",
            iface, mesh->yface_level[iface], yface_level_check[iface]);
         icount_err++;
      }
      if (mesh->map_yface2cell_lower[iface] != map_yface2cell_lower_check[iface]){
         printf("Error -- map_yface2cell_lower does not match for face %d map_yface2cell_lower %d map_yface2cell_lower_check %d\n",
            iface, mesh->map_yface2cell_lower[iface], map_yface2cell_lower_check[iface]);
         icount_err++;
      }
      if (mesh->map_yface2cell_upper[iface] != map_yface2cell_upper_check[iface]){
         printf("Error -- map_yface2cell_upper does not match for face %d map_yface2cell_upper %d map_yface2cell_upper_check %d\n",
            iface, mesh->map_yface2cell_upper[iface], map_yface2cell_upper_check[iface]);
         icount_err++;
      }
   }

   //if (mype == 0){
      if (icount_err){
          printf("  Error with calc face with map\n");
      } else{
          printf("  PASSED calc face with map\n");
      }
   //}

   mesh->calc_face_list_wbidirmap();

   int *map_xcell2face_left1_check = set_map_xcell2face_left1_check(mesh->ncells);
   int *map_xcell2face_left2_check = set_map_xcell2face_left2_check(mesh->ncells);
   int *map_xcell2face_right1_check = set_map_xcell2face_right1_check(mesh->ncells);
   int *map_xcell2face_right2_check = set_map_xcell2face_right2_check(mesh->ncells);

   int *map_ycell2face_bot1_check = set_map_ycell2face_bot1_check(mesh->ncells);
   int *map_ycell2face_bot2_check = set_map_ycell2face_bot2_check(mesh->ncells);
   int *map_ycell2face_top1_check = set_map_ycell2face_top1_check(mesh->ncells);
   int *map_ycell2face_top2_check = set_map_ycell2face_top2_check(mesh->ncells);

   icount_err = 0;

   if (mesh->nxface != nxface_check){
      printf("Error -- nxface does not match, nxface %d nxface_check %d\n",
         mesh->nxface, nxface_check);
      icount_err++;
   }

   if (mesh->nyface != nyface_check){
      printf("Error -- nyface does not match, nyface %d nyface_check %d\n",
         mesh->nyface, nyface_check);
      icount_err++;
   }

   for (int iface = 0; iface < mesh->nxface; iface++){
      if (mesh->xface_i[iface] != xface_i_check[iface]){
         printf("Error -- xface_i does not match for face %d xface_i %d xface_i_check %d\n",
            iface, mesh->xface_i[iface], xface_i_check[iface]);
         icount_err++;
      }
      if (mesh->xface_j[iface] != xface_j_check[iface]){
         printf("Error -- xface_j does not match for face %d xface_j %d xface_j_check %d\n",
            iface, mesh->xface_j[iface], xface_j_check[iface]);
         icount_err++;
      }
      if (mesh->xface_level[iface] != xface_level_check[iface]){
         printf("Error -- xface_level does not match for face %d xface_level %d xface_level_check %d\n",
            iface, mesh->xface_level[iface], xface_level_check[iface]);
         icount_err++;
      }
      if (mesh->map_xface2cell_lower[iface] != map_xface2cell_lower_check[iface]){
         printf("Error -- map_xface2cell_lower does not match for face %d map_xface2cell_lower %d map_xface2cell_lower_check %d\n",
            iface, mesh->map_xface2cell_lower[iface], map_xface2cell_lower_check[iface]);
         icount_err++;
      }
      if (mesh->map_xface2cell_upper[iface] != map_xface2cell_upper_check[iface]){
         printf("Error -- map_xface2cell_upper does not match for face %d map_xface2cell_upper %d map_xface2cell_upper_check %d\n",
            iface, mesh->map_xface2cell_upper[iface], map_xface2cell_upper_check[iface]);
         icount_err++;
      }
   }

   for (int iface = 0; iface < mesh->nyface; iface++){
      if (mesh->yface_i[iface] != yface_i_check[iface]){
         printf("Error -- yface_i does not match for face %d yface_i %d yface_i_check %d\n",
            iface, mesh->yface_i[iface], yface_i_check[iface]);
         icount_err++;
      }
      if (mesh->yface_j[iface] != yface_j_check[iface]){
         printf("Error -- yface_j does not match for face %d yface_j %d yface_j_check %d\n",
            iface, mesh->yface_j[iface], yface_j_check[iface]);
         icount_err++;
      }
      if (mesh->yface_level[iface] != yface_level_check[iface]){
         printf("Error -- yface_level does not match for face %d yface_level %d yface_level_check %d\n",
            iface, mesh->yface_level[iface], yface_level_check[iface]);
         icount_err++;
      }
      if (mesh->map_yface2cell_lower[iface] != map_yface2cell_lower_check[iface]){
         printf("Error -- map_yface2cell_lower does not match for face %d map_yface2cell_lower %d map_yface2cell_lower_check %d\n",
            iface, mesh->map_yface2cell_lower[iface], map_yface2cell_lower_check[iface]);
         icount_err++;
      }
      if (mesh->map_yface2cell_upper[iface] != map_yface2cell_upper_check[iface]){
         printf("Error -- map_yface2cell_upper does not match for face %d map_yface2cell_upper %d map_yface2cell_upper_check %d\n",
            iface, mesh->map_yface2cell_upper[iface], map_yface2cell_upper_check[iface]);
         icount_err++;
      }
   }

   for (uint iz = 0; iz < mesh->ncells; iz++){
      if (mesh->map_xcell2face_left1[iz] != map_xcell2face_left1_check[iz]) {
         printf("Error -- map_xcell2face_left1 does not match for cell %d map_xcell2face_left1 %d check %d\n",
            iz,mesh->map_xcell2face_left1[iz], map_xcell2face_left1_check[iz]);
         //printf("   map_xcell2face_left1_check[%d] = %d;\n", iz,mesh->map_xcell2face_left1[iz]);
      }
      if (mesh->map_xcell2face_left2[iz] != map_xcell2face_left2_check[iz]) {
         printf("Error -- map_xcell2face_left2 does not match for cell %d map_xcell2face_left2 %d check %d\n",
            iz,mesh->map_xcell2face_left2[iz], map_xcell2face_left2_check[iz]);
         //printf("   map_xcell2face_left2_check[%d] = %d;\n", iz,mesh->map_xcell2face_left2[iz]);
      }
      if (mesh->map_xcell2face_right1[iz] != map_xcell2face_right1_check[iz]) {
         printf("Error -- map_xcell2face_right1 does not match for cell %d map_xcell2face_right1 %d check %d\n",
            iz,mesh->map_xcell2face_right1[iz], map_xcell2face_right1_check[iz]);
         //printf("   map_xcell2face_right1_check[%d] = %d;\n", iz,mesh->map_xcell2face_right1[iz]);
      }
      if (mesh->map_xcell2face_right2[iz] != map_xcell2face_right2_check[iz]) {
         printf("Error -- map_xcell2face_right2 does not match for cell %d map_xcell2face_right2 %d check %d\n",
            iz,mesh->map_xcell2face_right2[iz], map_xcell2face_right2_check[iz]);
         //printf("   map_xcell2face_right2_check[%d] = %d;\n", iz,mesh->map_xcell2face_right2[iz]);
      }

      if (mesh->map_ycell2face_bot1[iz] != map_ycell2face_bot1_check[iz]) {
         printf("Error -- map_ycell2face_bot1 does not match for cell %d map_ycell2face_bot1 %d check %d\n",
            iz,mesh->map_ycell2face_bot1[iz], map_ycell2face_bot1_check[iz]);
         //printf("   map_ycell2face_bot1_check[%d] = %d;\n", iz,mesh->map_ycell2face_bot1[iz]);
      }
      if (mesh->map_ycell2face_bot2[iz] != map_ycell2face_bot2_check[iz]) {
         printf("Error -- map_ycell2face_bot2 does not match for cell %d map_ycell2face_bot2 %d check %d\n",
            iz,mesh->map_ycell2face_bot2[iz], map_ycell2face_bot2_check[iz]);
         //printf("   map_ycell2face_bot2_check[%d] = %d;\n", iz,mesh->map_ycell2face_bot2[iz]);
      }
      if (mesh->map_ycell2face_top1[iz] != map_ycell2face_top1_check[iz]) {
         printf("Error -- map_ycell2face_top1 does not match for cell %d map_ycell2face_top1 %d check %d\n",
            iz,mesh->map_ycell2face_top1[iz], map_ycell2face_top1_check[iz]);
         //printf("   map_ycell2face_top1_check[%d] = %d;\n", iz,mesh->map_ycell2face_top1[iz]);
      }
      if (mesh->map_ycell2face_top2[iz] != map_ycell2face_top2_check[iz]) {
         printf("Error -- map_ycell2face_top2 does not match for cell %d map_ycell2face_top2 %d check %d\n",
            iz,mesh->map_ycell2face_top2[iz], map_ycell2face_top2_check[iz]);
         //printf("   map_ycell2face_top2_check[%d] = %d;\n", iz,mesh->map_ycell2face_top2[iz]);
      }
   }

   //if (mype == 0){
      if (icount_err){
          printf("  Error with calc face with bi-directional map\n");
      } else{
          printf("  PASSED calc face with bi-directional map\n");
      }
   //}

   icount_err = 0;

   for (int lev = 1; lev < levmx+1; lev++){

      int isize, jsize;

      if (lev == 0) {
         isize = 0;
         jsize = 0;
      } else if (lev == 1) {
         isize = 7;
         jsize = 8;
      } else if (lev == 2) {
         isize = 9;
         jsize = 8;
      }

      if (isize == 0 || jsize ==0) continue;

      int **xface_flag_check = set_xface_flag_check(jsize, isize, lev);

      int **xface_flag = mesh->get_xface_flag(lev);

      for (int jj = 0; jj < jsize; jj++){
         for (int ii = 0; ii < isize; ii++){
            if (xface_flag[jj][ii] != xface_flag_check[jj][ii]){
               printf("Error -- xface_flag does not match for j %d i %d flag %d check %d\n",
                  jj, ii, xface_flag[jj][ii], xface_flag_check[jj][ii]);
               icount_err++;
            }
         }
      }

      genmatrixfree((void **)xface_flag_check);
      genmatrixfree((void **)xface_flag);
   }

   //if (mype == 0){
      if (icount_err){
          printf("  Error with get_xface_flag\n");
      } else{
          printf("  PASSED get_xface_flag\n");
      }
   //}

   icount_err = 0;

   for (int lev = 1; lev < levmx+1; lev++){

      int isize, jsize;

      if (lev == 0) {
         isize = 0;
         jsize = 0;
      } else if (lev == 1) {
         isize = 8;
         jsize = 7;
      } else if (lev == 2) {
         isize = 8;
         jsize = 9;
      }

      if (isize == 0 || jsize ==0) continue;

      int **yface_flag_check = set_yface_flag_check(jsize, isize, lev);

      int **yface_flag = mesh->get_yface_flag(lev);

      for (int jj = 0; jj < jsize; jj++){
         for (int ii = 0; ii < isize; ii++){
            if (yface_flag[jj][ii] != yface_flag_check[jj][ii]){
               printf("Error -- yface_flag does not match for j %d i %d flag %d check %d\n",
                  jj, ii, yface_flag[jj][ii], yface_flag_check[jj][ii]);
               icount_err++;
            }
         }
      }

      genmatrixfree((void **)yface_flag_check);
      genmatrixfree((void **)yface_flag);
   }

   //if (mype == 0){
      if (icount_err){
          printf("  Error with get_yface_flag\n");
      } else{
          printf("  PASSED get_yface_flag\n");
      }
   //}

   icount_err = 0;

   for (int lev = 1; lev < levmx+1; lev++){

      int isize, jsize;

      int **zone_flag;
      int **zone_cell;

      if (lev == 1){
         isize = 10;
         jsize = 10;
      } else if (lev == 2){
         isize = 12;
         jsize = 12;
      }

      mesh->get_flat_grid(lev, &zone_flag, &zone_cell);

      int **zone_flag_check = set_zone_flag_check(jsize, isize, lev);
      int **zone_cell_check = set_zone_cell_check(jsize, isize, lev);

      for (int jj = 0; jj < jsize; jj++){
         for (int ii = 0; ii < isize; ii++){
            if (zone_flag[jj][ii] != zone_flag_check[jj][ii]){
               printf("Error -- zone_flag does not match for j %d i %d flag %d check %d\n",
                  jj, ii, zone_flag[jj][ii], zone_flag_check[jj][ii]);
               icount_err++;
            }
            if (zone_cell[jj][ii] != zone_cell_check[jj][ii]){
               printf("Error -- zone_cell does not match for j %d i %d cell %d check %d\n",
                  jj, ii, zone_cell[jj][ii], zone_cell_check[jj][ii]);
               icount_err++;
            }
         }
      }

      genmatrixfree((void **)zone_flag_check);
      genmatrixfree((void **)zone_flag);
      genmatrixfree((void **)zone_cell_check);
      genmatrixfree((void **)zone_cell);
   }

   //if (mype == 0){
      if (icount_err){
          printf("  Error with get_flat_grid\n");
      } else{
          printf("  PASSED get_flat_grid\n");
      }
   //}


   mesh->mesh_memory.memory_delete(mesh->i);
   mesh->mesh_memory.memory_delete(mesh->j);
   mesh->mesh_memory.memory_delete(mesh->level);
   mesh->mesh_memory.memory_delete(mesh->celltype);
   mesh->mesh_memory.memory_delete(mesh->nlft);
   mesh->mesh_memory.memory_delete(mesh->nrht);
   mesh->mesh_memory.memory_delete(mesh->nbot);
   mesh->mesh_memory.memory_delete(mesh->ntop);

   mesh->mesh_memory.memory_report();

   mesh->terminate();
   delete mesh;

   //L7_Terminate();

   exit(0);
}

void set_mesh_values(int *i, int *j, int *level)
{
   i[0] = 1;
   i[1] = 1;
   i[2] = 0;
   i[3] = 1;
   i[4] = 1;
   i[5] = 2;
   i[6] = 2;
   i[7] = 3;
   i[8] = 3;
   i[9] = 4;
   i[10] = 8;
   i[11] = 8;
   i[12] = 9;
   i[13] = 9;
   i[14] = 10;
   i[15] = 10;
   i[16] = 11;
   i[17] = 11;
   i[18] = 11;
   i[19] = 10;
   i[20] = 10;
   i[21] = 11;
   i[22] = 5;
   i[23] = 4;
   i[24] = 4;
   i[25] = 5;
   i[26] = 4;
   i[27] = 5;
   i[28] = 6;
   i[29] = 7;
   i[30] = 6;
   i[31] = 6;
   i[32] = 7;
   i[33] = 7;
   i[34] = 4;
   i[35] = 5;
   i[36] = 4;
   i[37] = 10;
   i[38] = 10;
   i[39] = 9;
   i[40] = 9;
   i[41] = 8;
   i[42] = 8;
   i[43] = 15;
   i[44] = 14;
   i[45] = 14;
   i[46] = 15;
   i[47] = 7;
   i[48] = 13;
   i[49] = 13;
   i[50] = 12;
   i[51] = 12;
   i[52] = 12;
   i[53] = 13;
   i[54] = 13;
   i[55] = 12;
   i[56] = 12;
   i[57] = 13;
   i[58] = 13;
   i[59] = 12;
   i[60] = 12;
   i[61] = 12;
   i[62] = 13;
   i[63] = 13;
   i[64] = 7;
   i[65] = 15;
   i[66] = 14;
   i[67] = 14;
   i[68] = 15;
   i[69] = 8;
   i[70] = 8;
   i[71] = 9;
   i[72] = 9;
   i[73] = 10;
   i[74] = 10;
   i[75] = 4;
   i[76] = 5;
   i[77] = 4;
   i[78] = 7;
   i[79] = 7;
   i[80] = 6;
   i[81] = 6;
   i[82] = 7;
   i[83] = 6;
   i[84] = 5;
   i[85] = 4;
   i[86] = 5;
   i[87] = 4;
   i[88] = 4;
   i[89] = 5;
   i[90] = 11;
   i[91] = 10;
   i[92] = 10;
   i[93] = 11;
   i[94] = 11;
   i[95] = 11;
   i[96] = 10;
   i[97] = 10;
   i[98] = 9;
   i[99] = 9;
   i[100] = 8;
   i[101] = 8;
   i[102] = 4;
   i[103] = 3;
   i[104] = 3;
   i[105] = 2;
   i[106] = 2;
   i[107] = 1;
   i[108] = 1;
   i[109] = 0;
   i[110] = 1;
   i[111] = 1;
   
   j[0] = 0;
   j[1] = 1;
   j[2] = 1;
   j[3] = 4;
   j[4] = 5;
   j[5] = 4;
   j[6] = 5;
   j[7] = 5;
   j[8] = 4;
   j[9] = 4;
   j[10] = 10;
   j[11] = 11;
   j[12] = 11;
   j[13] = 10;
   j[14] = 10;
   j[15] = 11;
   j[16] = 11;
   j[17] = 10;
   j[18] = 9;
   j[19] = 9;
   j[20] = 8;
   j[21] = 8;
   j[22] = 3;
   j[23] = 3;
   j[24] = 2;
   j[25] = 2;
   j[26] = 1;
   j[27] = 1;
   j[28] = 1;
   j[29] = 1;
   j[30] = 2;
   j[31] = 3;
   j[32] = 3;
   j[33] = 2;
   j[34] = 0;
   j[35] = 1;
   j[36] = 1;
   j[37] = 4;
   j[38] = 5;
   j[39] = 5;
   j[40] = 4;
   j[41] = 4;
   j[42] = 5;
   j[43] = 11;
   j[44] = 11;
   j[45] = 10;
   j[46] = 10;
   j[47] = 4;
   j[48] = 9;
   j[49] = 8;
   j[50] = 8;
   j[51] = 9;
   j[52] = 10;
   j[53] = 10;
   j[54] = 11;
   j[55] = 11;
   j[56] = 12;
   j[57] = 12;
   j[58] = 13;
   j[59] = 13;
   j[60] = 14;
   j[61] = 15;
   j[62] = 15;
   j[63] = 14;
   j[64] = 7;
   j[65] = 13;
   j[66] = 13;
   j[67] = 12;
   j[68] = 12;
   j[69] = 6;
   j[70] = 7;
   j[71] = 7;
   j[72] = 6;
   j[73] = 6;
   j[74] = 7;
   j[75] = 4;
   j[76] = 4;
   j[77] = 5;
   j[78] = 9;
   j[79] = 8;
   j[80] = 8;
   j[81] = 9;
   j[82] = 10;
   j[83] = 10;
   j[84] = 10;
   j[85] = 10;
   j[86] = 9;
   j[87] = 9;
   j[88] = 8;
   j[89] = 8;
   j[90] = 15;
   j[91] = 15;
   j[92] = 14;
   j[93] = 14;
   j[94] = 13;
   j[95] = 12;
   j[96] = 12;
   j[97] = 13;
   j[98] = 13;
   j[99] = 12;
   j[100] = 12;
   j[101] = 13;
   j[102] = 7;
   j[103] = 7;
   j[104] = 6;
   j[105] = 6;
   j[106] = 7;
   j[107] = 6;
   j[108] = 7;
   j[109] = 4;
   j[110] = 4;
   j[111] = 5;

   level[0] = 0;
   level[1] = 0;
   level[2] = 0;
   level[3] = 1;
   level[4] = 1;
   level[5] = 1;
   level[6] = 1;
   level[7] = 1;
   level[8] = 1;
   level[9] = 1;
   level[10] = 2;
   level[11] = 2;
   level[12] = 2;
   level[13] = 2;
   level[14] = 2;
   level[15] = 2;
   level[16] = 2;
   level[17] = 2;
   level[18] = 2;
   level[19] = 2;
   level[20] = 2;
   level[21] = 2;
   level[22] = 1;
   level[23] = 1;
   level[24] = 1;
   level[25] = 1;
   level[26] = 1;
   level[27] = 1;
   level[28] = 1;
   level[29] = 1;
   level[30] = 1;
   level[31] = 1;
   level[32] = 1;
   level[33] = 1;
   level[34] = 0;
   level[35] = 0;
   level[36] = 0;
   level[37] = 1;
   level[38] = 1;
   level[39] = 1;
   level[40] = 1;
   level[41] = 1;
   level[42] = 1;
   level[43] = 2;
   level[44] = 2;
   level[45] = 2;
   level[46] = 2;
   level[47] = 1;
   level[48] = 2;
   level[49] = 2;
   level[50] = 2;
   level[51] = 2;
   level[52] = 2;
   level[53] = 2;
   level[54] = 2;
   level[55] = 2;
   level[56] = 2;
   level[57] = 2;
   level[58] = 2;
   level[59] = 2;
   level[60] = 2;
   level[61] = 2;
   level[62] = 2;
   level[63] = 2;
   level[64] = 1;
   level[65] = 2;
   level[66] = 2;
   level[67] = 2;
   level[68] = 2;
   level[69] = 1;
   level[70] = 1;
   level[71] = 1;
   level[72] = 1;
   level[73] = 1;
   level[74] = 1;
   level[75] = 0;
   level[76] = 0;
   level[77] = 0;
   level[78] = 1;
   level[79] = 1;
   level[80] = 1;
   level[81] = 1;
   level[82] = 1;
   level[83] = 1;
   level[84] = 1;
   level[85] = 1;
   level[86] = 1;
   level[87] = 1;
   level[88] = 1;
   level[89] = 1;
   level[90] = 2;
   level[91] = 2;
   level[92] = 2;
   level[93] = 2;
   level[94] = 2;
   level[95] = 2;
   level[96] = 2;
   level[97] = 2;
   level[98] = 2;
   level[99] = 2;
   level[100] = 2;
   level[101] = 2;
   level[102] = 1;
   level[103] = 1;
   level[104] = 1;
   level[105] = 1;
   level[106] = 1;
   level[107] = 1;
   level[108] = 1;
   level[109] = 0;
   level[110] = 0;
   level[111] = 0;
}

void set_xface_check(int *xface_i_check, int *xface_j_check, int *xface_level_check)
{
   xface_i_check[0] = 1;
   xface_i_check[1] = 1;
   xface_i_check[2] = 0;
   xface_i_check[3] = 0;
   xface_i_check[4] = 0;
   xface_i_check[5] = 0;
   xface_i_check[6] = 1;
   xface_i_check[7] = 2;
   xface_i_check[8] = 2;
   xface_i_check[9] = 1;
   xface_i_check[10] = 1;
   xface_i_check[11] = 2;
   xface_i_check[12] = 2;
   xface_i_check[13] = 3;
   xface_i_check[14] = 3;
   xface_i_check[15] = 4;
   xface_i_check[16] = 4;
   xface_i_check[17] = 4;
   xface_i_check[18] = 3;
   xface_i_check[19] = 3;
   xface_i_check[20] = 4;
   xface_i_check[21] = 3;
   xface_i_check[22] = 2;
   xface_i_check[23] = 2;
   xface_i_check[24] = 3;
   xface_i_check[25] = 4;
   xface_i_check[26] = 4;
   xface_i_check[27] = 5;
   xface_i_check[28] = 5;
   xface_i_check[29] = 6;
   xface_i_check[30] = 6;
   xface_i_check[31] = 8;
   xface_i_check[32] = 7;
   xface_i_check[33] = 7;
   xface_i_check[34] = 8;
   xface_i_check[35] = 5;
   xface_i_check[36] = 6;
   xface_i_check[37] = 6;
   xface_i_check[38] = 5;
   xface_i_check[39] = 5;
   xface_i_check[40] = 5;
   xface_i_check[41] = 6;
   xface_i_check[42] = 6;
   xface_i_check[43] = 5;
   xface_i_check[44] = 5;
   xface_i_check[45] = 6;
   xface_i_check[46] = 6;
   xface_i_check[47] = 5;
   xface_i_check[48] = 5;
   xface_i_check[49] = 5;
   xface_i_check[50] = 6;
   xface_i_check[51] = 6;
   xface_i_check[52] = 5;
   xface_i_check[53] = 8;
   xface_i_check[54] = 7;
   xface_i_check[55] = 7;
   xface_i_check[56] = 8;
   xface_i_check[57] = 6;
   xface_i_check[58] = 6;
   xface_i_check[59] = 5;
   xface_i_check[60] = 5;
   xface_i_check[61] = 4;
   xface_i_check[62] = 4;
   xface_i_check[63] = 3;
   xface_i_check[64] = 2;
   xface_i_check[65] = 2;
   xface_i_check[66] = 3;
   xface_i_check[67] = 4;
   xface_i_check[68] = 3;
   xface_i_check[69] = 3;
   xface_i_check[70] = 4;
   xface_i_check[71] = 4;
   xface_i_check[72] = 4;
   xface_i_check[73] = 3;
   xface_i_check[74] = 3;
   xface_i_check[75] = 2;
   xface_i_check[76] = 2;
   xface_i_check[77] = 1;
   xface_i_check[78] = 1;
   xface_i_check[79] = 2;
   xface_i_check[80] = 2;
   xface_i_check[81] = 1;
   xface_i_check[82] = 0;
   xface_i_check[83] = 0;
   xface_i_check[84] = 0;
   xface_i_check[85] = 0;
   xface_i_check[86] = 1;
   xface_i_check[87] = 1;

   xface_j_check[0] = 0;
   xface_j_check[1] = 1;
   xface_j_check[2] = 2;
   xface_j_check[3] = 3;
   xface_j_check[4] = 2;
   xface_j_check[5] = 3;
   xface_j_check[6] = 2;
   xface_j_check[7] = 0;
   xface_j_check[8] = 1;
   xface_j_check[9] = 2;
   xface_j_check[10] = 3;
   xface_j_check[11] = 3;
   xface_j_check[12] = 2;
   xface_j_check[13] = 2;
   xface_j_check[14] = 3;
   xface_j_check[15] = 3;
   xface_j_check[16] = 2;
   xface_j_check[17] = 1;
   xface_j_check[18] = 1;
   xface_j_check[19] = 0;
   xface_j_check[20] = 0;
   xface_j_check[21] = 1;
   xface_j_check[22] = 1;
   xface_j_check[23] = 0;
   xface_j_check[24] = 0;
   xface_j_check[25] = 0;
   xface_j_check[26] = 1;
   xface_j_check[27] = 1;
   xface_j_check[28] = 0;
   xface_j_check[29] = 2;
   xface_j_check[30] = 3;
   xface_j_check[31] = 3;
   xface_j_check[32] = 3;
   xface_j_check[33] = 2;
   xface_j_check[34] = 2;
   xface_j_check[35] = 2;
   xface_j_check[36] = 1;
   xface_j_check[37] = 0;
   xface_j_check[38] = 0;
   xface_j_check[39] = 1;
   xface_j_check[40] = 2;
   xface_j_check[41] = 2;
   xface_j_check[42] = 3;
   xface_j_check[43] = 3;
   xface_j_check[44] = 4;
   xface_j_check[45] = 4;
   xface_j_check[46] = 5;
   xface_j_check[47] = 5;
   xface_j_check[48] = 6;
   xface_j_check[49] = 7;
   xface_j_check[50] = 7;
   xface_j_check[51] = 6;
   xface_j_check[52] = 5;
   xface_j_check[53] = 5;
   xface_j_check[54] = 5;
   xface_j_check[55] = 4;
   xface_j_check[56] = 4;
   xface_j_check[57] = 4;
   xface_j_check[58] = 5;
   xface_j_check[59] = 7;
   xface_j_check[60] = 6;
   xface_j_check[61] = 6;
   xface_j_check[62] = 7;
   xface_j_check[63] = 7;
   xface_j_check[64] = 7;
   xface_j_check[65] = 6;
   xface_j_check[66] = 6;
   xface_j_check[67] = 7;
   xface_j_check[68] = 7;
   xface_j_check[69] = 6;
   xface_j_check[70] = 6;
   xface_j_check[71] = 5;
   xface_j_check[72] = 4;
   xface_j_check[73] = 4;
   xface_j_check[74] = 5;
   xface_j_check[75] = 5;
   xface_j_check[76] = 4;
   xface_j_check[77] = 4;
   xface_j_check[78] = 5;
   xface_j_check[79] = 6;
   xface_j_check[80] = 7;
   xface_j_check[81] = 5;
   xface_j_check[82] = 4;
   xface_j_check[83] = 5;
   xface_j_check[84] = 4;
   xface_j_check[85] = 5;
   xface_j_check[86] = 6;
   xface_j_check[87] = 7;

   xface_level_check[0] = 1;
   xface_level_check[1] = 1;
   xface_level_check[2] = 1;
   xface_level_check[3] = 1;
   xface_level_check[4] = 2;
   xface_level_check[5] = 2;
   xface_level_check[6] = 1;
   xface_level_check[7] = 2;
   xface_level_check[8] = 2;
   xface_level_check[9] = 2;
   xface_level_check[10] = 2;
   xface_level_check[11] = 2;
   xface_level_check[12] = 2;
   xface_level_check[13] = 2;
   xface_level_check[14] = 2;
   xface_level_check[15] = 2;
   xface_level_check[16] = 2;
   xface_level_check[17] = 2;
   xface_level_check[18] = 2;
   xface_level_check[19] = 2;
   xface_level_check[20] = 2;
   xface_level_check[21] = 1;
   xface_level_check[22] = 1;
   xface_level_check[23] = 1;
   xface_level_check[24] = 1;
   xface_level_check[25] = 1;
   xface_level_check[26] = 1;
   xface_level_check[27] = 1;
   xface_level_check[28] = 1;
   xface_level_check[29] = 1;
   xface_level_check[30] = 1;
   xface_level_check[31] = 2;
   xface_level_check[32] = 2;
   xface_level_check[33] = 2;
   xface_level_check[34] = 2;
   xface_level_check[35] = 1;
   xface_level_check[36] = 2;
   xface_level_check[37] = 2;
   xface_level_check[38] = 2;
   xface_level_check[39] = 2;
   xface_level_check[40] = 2;
   xface_level_check[41] = 2;
   xface_level_check[42] = 2;
   xface_level_check[43] = 2;
   xface_level_check[44] = 2;
   xface_level_check[45] = 2;
   xface_level_check[46] = 2;
   xface_level_check[47] = 2;
   xface_level_check[48] = 2;
   xface_level_check[49] = 2;
   xface_level_check[50] = 2;
   xface_level_check[51] = 2;
   xface_level_check[52] = 1;
   xface_level_check[53] = 2;
   xface_level_check[54] = 2;
   xface_level_check[55] = 2;
   xface_level_check[56] = 2;
   xface_level_check[57] = 1;
   xface_level_check[58] = 1;
   xface_level_check[59] = 1;
   xface_level_check[60] = 1;
   xface_level_check[61] = 1;
   xface_level_check[62] = 1;
   xface_level_check[63] = 1;
   xface_level_check[64] = 1;
   xface_level_check[65] = 1;
   xface_level_check[66] = 1;
   xface_level_check[67] = 2;
   xface_level_check[68] = 2;
   xface_level_check[69] = 2;
   xface_level_check[70] = 2;
   xface_level_check[71] = 2;
   xface_level_check[72] = 2;
   xface_level_check[73] = 2;
   xface_level_check[74] = 2;
   xface_level_check[75] = 2;
   xface_level_check[76] = 2;
   xface_level_check[77] = 2;
   xface_level_check[78] = 2;
   xface_level_check[79] = 2;
   xface_level_check[80] = 2;
   xface_level_check[81] = 1;
   xface_level_check[82] = 2;
   xface_level_check[83] = 2;
   xface_level_check[84] = 1;
   xface_level_check[85] = 1;
   xface_level_check[86] = 1;
   xface_level_check[87] = 1;
}

void set_yface_check(int *yface_i_check, int *yface_j_check, int *yface_level_check)
{
   yface_i_check[0] = 0;
   yface_i_check[1] = 1;
   yface_i_check[2] = 0;
   yface_i_check[3] = 0;
   yface_i_check[4] = 1;
   yface_i_check[5] = 1;
   yface_i_check[6] = 0;
   yface_i_check[7] = 1;
   yface_i_check[8] = 0;
   yface_i_check[9] = 0;
   yface_i_check[10] = 1;
   yface_i_check[11] = 1;
   yface_i_check[12] = 2;
   yface_i_check[13] = 2;
   yface_i_check[14] = 3;
   yface_i_check[15] = 3;
   yface_i_check[16] = 3;
   yface_i_check[17] = 2;
   yface_i_check[18] = 2;
   yface_i_check[19] = 3;
   yface_i_check[20] = 2;
   yface_i_check[21] = 3;
   yface_i_check[22] = 2;
   yface_i_check[23] = 2;
   yface_i_check[24] = 3;
   yface_i_check[25] = 4;
   yface_i_check[26] = 4;
   yface_i_check[27] = 5;
   yface_i_check[28] = 5;
   yface_i_check[29] = 5;
   yface_i_check[30] = 6;
   yface_i_check[31] = 7;
   yface_i_check[32] = 7;
   yface_i_check[33] = 7;
   yface_i_check[34] = 6;
   yface_i_check[35] = 6;
   yface_i_check[36] = 7;
   yface_i_check[37] = 6;
   yface_i_check[38] = 6;
   yface_i_check[39] = 7;
   yface_i_check[40] = 6;
   yface_i_check[41] = 7;
   yface_i_check[42] = 5;
   yface_i_check[43] = 5;
   yface_i_check[44] = 4;
   yface_i_check[45] = 4;
   yface_i_check[46] = 4;
   yface_i_check[47] = 5;
   yface_i_check[48] = 5;
   yface_i_check[49] = 4;
   yface_i_check[50] = 4;
   yface_i_check[51] = 5;
   yface_i_check[52] = 5;
   yface_i_check[53] = 4;
   yface_i_check[54] = 4;
   yface_i_check[55] = 4;
   yface_i_check[56] = 5;
   yface_i_check[57] = 5;
   yface_i_check[58] = 5;
   yface_i_check[59] = 7;
   yface_i_check[60] = 6;
   yface_i_check[61] = 6;
   yface_i_check[62] = 7;
   yface_i_check[63] = 6;
   yface_i_check[64] = 6;
   yface_i_check[65] = 7;
   yface_i_check[66] = 7;
   yface_i_check[67] = 5;
   yface_i_check[68] = 4;
   yface_i_check[69] = 2;
   yface_i_check[70] = 3;
   yface_i_check[71] = 3;
   yface_i_check[72] = 2;
   yface_i_check[73] = 2;
   yface_i_check[74] = 3;
   yface_i_check[75] = 3;
   yface_i_check[76] = 3;
   yface_i_check[77] = 2;
   yface_i_check[78] = 2;
   yface_i_check[79] = 1;
   yface_i_check[80] = 1;
   yface_i_check[81] = 0;
   yface_i_check[82] = 0;
   yface_i_check[83] = 2;
   yface_i_check[84] = 1;
   yface_i_check[85] = 1;
   yface_i_check[86] = 0;
   yface_i_check[87] = 0;

   yface_j_check[0] = 1;
   yface_j_check[1] = 1;
   yface_j_check[2] = 2;
   yface_j_check[3] = 3;
   yface_j_check[4] = 3;
   yface_j_check[5] = 2;
   yface_j_check[6] = 2;
   yface_j_check[7] = 2;
   yface_j_check[8] = 3;
   yface_j_check[9] = 4;
   yface_j_check[10] = 4;
   yface_j_check[11] = 3;
   yface_j_check[12] = 3;
   yface_j_check[13] = 4;
   yface_j_check[14] = 4;
   yface_j_check[15] = 3;
   yface_j_check[16] = 2;
   yface_j_check[17] = 2;
   yface_j_check[18] = 1;
   yface_j_check[19] = 1;
   yface_j_check[20] = 0;
   yface_j_check[21] = 0;
   yface_j_check[22] = 1;
   yface_j_check[23] = 0;
   yface_j_check[24] = 0;
   yface_j_check[25] = 0;
   yface_j_check[26] = 0;
   yface_j_check[27] = 0;
   yface_j_check[28] = 1;
   yface_j_check[29] = 0;
   yface_j_check[30] = 1;
   yface_j_check[31] = 1;
   yface_j_check[32] = 3;
   yface_j_check[33] = 2;
   yface_j_check[34] = 2;
   yface_j_check[35] = 3;
   yface_j_check[36] = 4;
   yface_j_check[37] = 4;
   yface_j_check[38] = 3;
   yface_j_check[39] = 3;
   yface_j_check[40] = 2;
   yface_j_check[41] = 2;
   yface_j_check[42] = 2;
   yface_j_check[43] = 1;
   yface_j_check[44] = 1;
   yface_j_check[45] = 2;
   yface_j_check[46] = 3;
   yface_j_check[47] = 3;
   yface_j_check[48] = 4;
   yface_j_check[49] = 4;
   yface_j_check[50] = 5;
   yface_j_check[51] = 5;
   yface_j_check[52] = 6;
   yface_j_check[53] = 6;
   yface_j_check[54] = 7;
   yface_j_check[55] = 8;
   yface_j_check[56] = 8;
   yface_j_check[57] = 7;
   yface_j_check[58] = 5;
   yface_j_check[59] = 6;
   yface_j_check[60] = 6;
   yface_j_check[61] = 5;
   yface_j_check[62] = 5;
   yface_j_check[63] = 4;
   yface_j_check[64] = 5;
   yface_j_check[65] = 5;
   yface_j_check[66] = 4;
   yface_j_check[67] = 6;
   yface_j_check[68] = 6;
   yface_j_check[69] = 6;
   yface_j_check[70] = 6;
   yface_j_check[71] = 8;
   yface_j_check[72] = 8;
   yface_j_check[73] = 7;
   yface_j_check[74] = 7;
   yface_j_check[75] = 6;
   yface_j_check[76] = 5;
   yface_j_check[77] = 5;
   yface_j_check[78] = 6;
   yface_j_check[79] = 6;
   yface_j_check[80] = 5;
   yface_j_check[81] = 5;
   yface_j_check[82] = 6;
   yface_j_check[83] = 5;
   yface_j_check[84] = 5;
   yface_j_check[85] = 4;
   yface_j_check[86] = 4;
   yface_j_check[87] = 5;

   yface_level_check[0] = 1;
   yface_level_check[1] = 1;
   yface_level_check[2] = 1;
   yface_level_check[3] = 1;
   yface_level_check[4] = 1;
   yface_level_check[5] = 1;
   yface_level_check[6] = 2;
   yface_level_check[7] = 2;
   yface_level_check[8] = 2;
   yface_level_check[9] = 2;
   yface_level_check[10] = 2;
   yface_level_check[11] = 2;
   yface_level_check[12] = 2;
   yface_level_check[13] = 2;
   yface_level_check[14] = 2;
   yface_level_check[15] = 2;
   yface_level_check[16] = 2;
   yface_level_check[17] = 2;
   yface_level_check[18] = 2;
   yface_level_check[19] = 2;
   yface_level_check[20] = 2;
   yface_level_check[21] = 2;
   yface_level_check[22] = 1;
   yface_level_check[23] = 1;
   yface_level_check[24] = 1;
   yface_level_check[25] = 1;
   yface_level_check[26] = 2;
   yface_level_check[27] = 2;
   yface_level_check[28] = 1;
   yface_level_check[29] = 1;
   yface_level_check[30] = 1;
   yface_level_check[31] = 1;
   yface_level_check[32] = 1;
   yface_level_check[33] = 1;
   yface_level_check[34] = 1;
   yface_level_check[35] = 1;
   yface_level_check[36] = 2;
   yface_level_check[37] = 2;
   yface_level_check[38] = 2;
   yface_level_check[39] = 2;
   yface_level_check[40] = 2;
   yface_level_check[41] = 2;
   yface_level_check[42] = 2;
   yface_level_check[43] = 2;
   yface_level_check[44] = 2;
   yface_level_check[45] = 2;
   yface_level_check[46] = 2;
   yface_level_check[47] = 2;
   yface_level_check[48] = 2;
   yface_level_check[49] = 2;
   yface_level_check[50] = 2;
   yface_level_check[51] = 2;
   yface_level_check[52] = 2;
   yface_level_check[53] = 2;
   yface_level_check[54] = 2;
   yface_level_check[55] = 2;
   yface_level_check[56] = 2;
   yface_level_check[57] = 2;
   yface_level_check[58] = 1;
   yface_level_check[59] = 2;
   yface_level_check[60] = 2;
   yface_level_check[61] = 2;
   yface_level_check[62] = 2;
   yface_level_check[63] = 1;
   yface_level_check[64] = 1;
   yface_level_check[65] = 1;
   yface_level_check[66] = 1;
   yface_level_check[67] = 1;
   yface_level_check[68] = 1;
   yface_level_check[69] = 1;
   yface_level_check[70] = 1;
   yface_level_check[71] = 2;
   yface_level_check[72] = 2;
   yface_level_check[73] = 2;
   yface_level_check[74] = 2;
   yface_level_check[75] = 2;
   yface_level_check[76] = 2;
   yface_level_check[77] = 2;
   yface_level_check[78] = 2;
   yface_level_check[79] = 2;
   yface_level_check[80] = 2;
   yface_level_check[81] = 2;
   yface_level_check[82] = 2;
   yface_level_check[83] = 1;
   yface_level_check[84] = 1;
   yface_level_check[85] = 1;
   yface_level_check[86] = 1;
   yface_level_check[87] = 1;
}

void set_map_xface2cell_check(int *map_xface2cell_lower_check, int *map_xface2cell_upper_check)
{
   map_xface2cell_lower_check[0] = 1;
   map_xface2cell_lower_check[1] = 1;
   map_xface2cell_lower_check[2] = 5;
   map_xface2cell_lower_check[3] = 6;
   map_xface2cell_lower_check[4] = 7;
   map_xface2cell_lower_check[5] = 7;
   map_xface2cell_lower_check[6] = 8;
   map_xface2cell_lower_check[7] = 9;
   map_xface2cell_lower_check[8] = 9;
   map_xface2cell_lower_check[9] = 10;
   map_xface2cell_lower_check[10] = 11;
   map_xface2cell_lower_check[11] = 12;
   map_xface2cell_lower_check[12] = 13;
   map_xface2cell_lower_check[13] = 14;
   map_xface2cell_lower_check[14] = 15;
   map_xface2cell_lower_check[15] = 16;
   map_xface2cell_lower_check[16] = 17;
   map_xface2cell_lower_check[17] = 18;
   map_xface2cell_lower_check[18] = 19;
   map_xface2cell_lower_check[19] = 20;
   map_xface2cell_lower_check[20] = 21;
   map_xface2cell_lower_check[21] = 22;
   map_xface2cell_lower_check[22] = 23;
   map_xface2cell_lower_check[23] = 24;
   map_xface2cell_lower_check[24] = 25;
   map_xface2cell_lower_check[25] = 30;
   map_xface2cell_lower_check[26] = 31;
   map_xface2cell_lower_check[27] = 32;
   map_xface2cell_lower_check[28] = 33;
   map_xface2cell_lower_check[29] = 41;
   map_xface2cell_lower_check[30] = 42;
   map_xface2cell_lower_check[31] = 43;
   map_xface2cell_lower_check[32] = 44;
   map_xface2cell_lower_check[33] = 45;
   map_xface2cell_lower_check[34] = 46;
   map_xface2cell_lower_check[35] = 47;
   map_xface2cell_lower_check[36] = 48;
   map_xface2cell_lower_check[37] = 49;
   map_xface2cell_lower_check[38] = 50;
   map_xface2cell_lower_check[39] = 51;
   map_xface2cell_lower_check[40] = 52;
   map_xface2cell_lower_check[41] = 53;
   map_xface2cell_lower_check[42] = 54;
   map_xface2cell_lower_check[43] = 55;
   map_xface2cell_lower_check[44] = 56;
   map_xface2cell_lower_check[45] = 57;
   map_xface2cell_lower_check[46] = 58;
   map_xface2cell_lower_check[47] = 59;
   map_xface2cell_lower_check[48] = 60;
   map_xface2cell_lower_check[49] = 61;
   map_xface2cell_lower_check[50] = 62;
   map_xface2cell_lower_check[51] = 63;
   map_xface2cell_lower_check[52] = 64;
   map_xface2cell_lower_check[53] = 65;
   map_xface2cell_lower_check[54] = 66;
   map_xface2cell_lower_check[55] = 67;
   map_xface2cell_lower_check[56] = 68;
   map_xface2cell_lower_check[57] = 69;
   map_xface2cell_lower_check[58] = 70;
   map_xface2cell_lower_check[59] = 78;
   map_xface2cell_lower_check[60] = 79;
   map_xface2cell_lower_check[61] = 80;
   map_xface2cell_lower_check[62] = 81;
   map_xface2cell_lower_check[63] = 86;
   map_xface2cell_lower_check[64] = 87;
   map_xface2cell_lower_check[65] = 88;
   map_xface2cell_lower_check[66] = 89;
   map_xface2cell_lower_check[67] = 90;
   map_xface2cell_lower_check[68] = 91;
   map_xface2cell_lower_check[69] = 92;
   map_xface2cell_lower_check[70] = 93;
   map_xface2cell_lower_check[71] = 94;
   map_xface2cell_lower_check[72] = 95;
   map_xface2cell_lower_check[73] = 96;
   map_xface2cell_lower_check[74] = 97;
   map_xface2cell_lower_check[75] = 98;
   map_xface2cell_lower_check[76] = 99;
   map_xface2cell_lower_check[77] = 100;
   map_xface2cell_lower_check[78] = 101;
   map_xface2cell_lower_check[79] = 102;
   map_xface2cell_lower_check[80] = 102;
   map_xface2cell_lower_check[81] = 103;
   map_xface2cell_lower_check[82] = 104;
   map_xface2cell_lower_check[83] = 104;
   map_xface2cell_lower_check[84] = 105;
   map_xface2cell_lower_check[85] = 106;
   map_xface2cell_lower_check[86] = 110;
   map_xface2cell_lower_check[87] = 110;

   map_xface2cell_upper_check[0] = 24;
   map_xface2cell_upper_check[1] = 23;
   map_xface2cell_upper_check[2] = 8;
   map_xface2cell_upper_check[3] = 7;
   map_xface2cell_upper_check[4] = 10;
   map_xface2cell_upper_check[5] = 11;
   map_xface2cell_upper_check[6] = 9;
   map_xface2cell_upper_check[7] = 20;
   map_xface2cell_upper_check[8] = 19;
   map_xface2cell_upper_check[9] = 13;
   map_xface2cell_upper_check[10] = 12;
   map_xface2cell_upper_check[11] = 15;
   map_xface2cell_upper_check[12] = 14;
   map_xface2cell_upper_check[13] = 17;
   map_xface2cell_upper_check[14] = 16;
   map_xface2cell_upper_check[15] = 55;
   map_xface2cell_upper_check[16] = 52;
   map_xface2cell_upper_check[17] = 51;
   map_xface2cell_upper_check[18] = 18;
   map_xface2cell_upper_check[19] = 21;
   map_xface2cell_upper_check[20] = 50;
   map_xface2cell_upper_check[21] = 31;
   map_xface2cell_upper_check[22] = 22;
   map_xface2cell_upper_check[23] = 25;
   map_xface2cell_upper_check[24] = 30;
   map_xface2cell_upper_check[25] = 33;
   map_xface2cell_upper_check[26] = 32;
   map_xface2cell_upper_check[27] = 36;
   map_xface2cell_upper_check[28] = 36;
   map_xface2cell_upper_check[29] = 40;
   map_xface2cell_upper_check[30] = 39;
   map_xface2cell_upper_check[31] = 42;
   map_xface2cell_upper_check[32] = 43;
   map_xface2cell_upper_check[33] = 46;
   map_xface2cell_upper_check[34] = 42;
   map_xface2cell_upper_check[35] = 41;
   map_xface2cell_upper_check[36] = 47;
   map_xface2cell_upper_check[37] = 47;
   map_xface2cell_upper_check[38] = 49;
   map_xface2cell_upper_check[39] = 48;
   map_xface2cell_upper_check[40] = 53;
   map_xface2cell_upper_check[41] = 45;
   map_xface2cell_upper_check[42] = 44;
   map_xface2cell_upper_check[43] = 54;
   map_xface2cell_upper_check[44] = 57;
   map_xface2cell_upper_check[45] = 67;
   map_xface2cell_upper_check[46] = 66;
   map_xface2cell_upper_check[47] = 58;
   map_xface2cell_upper_check[48] = 63;
   map_xface2cell_upper_check[49] = 62;
   map_xface2cell_upper_check[50] = 64;
   map_xface2cell_upper_check[51] = 64;
   map_xface2cell_upper_check[52] = 70;
   map_xface2cell_upper_check[53] = 69;
   map_xface2cell_upper_check[54] = 65;
   map_xface2cell_upper_check[55] = 68;
   map_xface2cell_upper_check[56] = 69;
   map_xface2cell_upper_check[57] = 72;
   map_xface2cell_upper_check[58] = 71;
   map_xface2cell_upper_check[59] = 75;
   map_xface2cell_upper_check[60] = 75;
   map_xface2cell_upper_check[61] = 79;
   map_xface2cell_upper_check[62] = 78;
   map_xface2cell_upper_check[63] = 81;
   map_xface2cell_upper_check[64] = 86;
   map_xface2cell_upper_check[65] = 89;
   map_xface2cell_upper_check[66] = 80;
   map_xface2cell_upper_check[67] = 61;
   map_xface2cell_upper_check[68] = 90;
   map_xface2cell_upper_check[69] = 93;
   map_xface2cell_upper_check[70] = 60;
   map_xface2cell_upper_check[71] = 59;
   map_xface2cell_upper_check[72] = 56;
   map_xface2cell_upper_check[73] = 95;
   map_xface2cell_upper_check[74] = 94;
   map_xface2cell_upper_check[75] = 97;
   map_xface2cell_upper_check[76] = 96;
   map_xface2cell_upper_check[77] = 99;
   map_xface2cell_upper_check[78] = 98;
   map_xface2cell_upper_check[79] = 92;
   map_xface2cell_upper_check[80] = 91;
   map_xface2cell_upper_check[81] = 102;
   map_xface2cell_upper_check[82] = 100;
   map_xface2cell_upper_check[83] = 101;
   map_xface2cell_upper_check[84] = 104;
   map_xface2cell_upper_check[85] = 103;
   map_xface2cell_upper_check[86] = 88;
   map_xface2cell_upper_check[87] = 87;
}

void set_map_yface2cell_check(int *map_yface2cell_lower_check, int *map_yface2cell_upper_check)
{
   map_yface2cell_lower_check[0] = 1;
   map_yface2cell_lower_check[1] = 1;
   map_yface2cell_lower_check[2] = 5;
   map_yface2cell_lower_check[3] = 6;
   map_yface2cell_lower_check[4] = 7;
   map_yface2cell_lower_check[5] = 8;
   map_yface2cell_lower_check[6] = 9;
   map_yface2cell_lower_check[7] = 9;
   map_yface2cell_lower_check[8] = 10;
   map_yface2cell_lower_check[9] = 11;
   map_yface2cell_lower_check[10] = 12;
   map_yface2cell_lower_check[11] = 13;
   map_yface2cell_lower_check[12] = 14;
   map_yface2cell_lower_check[13] = 15;
   map_yface2cell_lower_check[14] = 16;
   map_yface2cell_lower_check[15] = 17;
   map_yface2cell_lower_check[16] = 18;
   map_yface2cell_lower_check[17] = 19;
   map_yface2cell_lower_check[18] = 20;
   map_yface2cell_lower_check[19] = 21;
   map_yface2cell_lower_check[20] = 22;
   map_yface2cell_lower_check[21] = 22;
   map_yface2cell_lower_check[22] = 23;
   map_yface2cell_lower_check[23] = 24;
   map_yface2cell_lower_check[24] = 25;
   map_yface2cell_lower_check[25] = 30;
   map_yface2cell_lower_check[26] = 31;
   map_yface2cell_lower_check[27] = 31;
   map_yface2cell_lower_check[28] = 32;
   map_yface2cell_lower_check[29] = 33;
   map_yface2cell_lower_check[30] = 36;
   map_yface2cell_lower_check[31] = 36;
   map_yface2cell_lower_check[32] = 39;
   map_yface2cell_lower_check[33] = 40;
   map_yface2cell_lower_check[34] = 41;
   map_yface2cell_lower_check[35] = 42;
   map_yface2cell_lower_check[36] = 43;
   map_yface2cell_lower_check[37] = 44;
   map_yface2cell_lower_check[38] = 45;
   map_yface2cell_lower_check[39] = 46;
   map_yface2cell_lower_check[40] = 47;
   map_yface2cell_lower_check[41] = 47;
   map_yface2cell_lower_check[42] = 48;
   map_yface2cell_lower_check[43] = 49;
   map_yface2cell_lower_check[44] = 50;
   map_yface2cell_lower_check[45] = 51;
   map_yface2cell_lower_check[46] = 52;
   map_yface2cell_lower_check[47] = 53;
   map_yface2cell_lower_check[48] = 54;
   map_yface2cell_lower_check[49] = 55;
   map_yface2cell_lower_check[50] = 56;
   map_yface2cell_lower_check[51] = 57;
   map_yface2cell_lower_check[52] = 58;
   map_yface2cell_lower_check[53] = 59;
   map_yface2cell_lower_check[54] = 60;
   map_yface2cell_lower_check[55] = 61;
   map_yface2cell_lower_check[56] = 62;
   map_yface2cell_lower_check[57] = 63;
   map_yface2cell_lower_check[58] = 64;
   map_yface2cell_lower_check[59] = 65;
   map_yface2cell_lower_check[60] = 66;
   map_yface2cell_lower_check[61] = 67;
   map_yface2cell_lower_check[62] = 68;
   map_yface2cell_lower_check[63] = 69;
   map_yface2cell_lower_check[64] = 70;
   map_yface2cell_lower_check[65] = 71;
   map_yface2cell_lower_check[66] = 72;
   map_yface2cell_lower_check[67] = 79;
   map_yface2cell_lower_check[68] = 80;
   map_yface2cell_lower_check[69] = 88;
   map_yface2cell_lower_check[70] = 89;
   map_yface2cell_lower_check[71] = 90;
   map_yface2cell_lower_check[72] = 91;
   map_yface2cell_lower_check[73] = 92;
   map_yface2cell_lower_check[74] = 93;
   map_yface2cell_lower_check[75] = 94;
   map_yface2cell_lower_check[76] = 95;
   map_yface2cell_lower_check[77] = 96;
   map_yface2cell_lower_check[78] = 97;
   map_yface2cell_lower_check[79] = 98;
   map_yface2cell_lower_check[80] = 99;
   map_yface2cell_lower_check[81] = 100;
   map_yface2cell_lower_check[82] = 101;
   map_yface2cell_lower_check[83] = 102;
   map_yface2cell_lower_check[84] = 103;
   map_yface2cell_lower_check[85] = 104;
   map_yface2cell_lower_check[86] = 105;
   map_yface2cell_lower_check[87] = 106;

   map_yface2cell_upper_check[0] = 5;
   map_yface2cell_upper_check[1] = 8;
   map_yface2cell_upper_check[2] = 6;
   map_yface2cell_upper_check[3] = 105;
   map_yface2cell_upper_check[4] = 104;
   map_yface2cell_upper_check[5] = 7;
   map_yface2cell_upper_check[6] = 10;
   map_yface2cell_upper_check[7] = 13;
   map_yface2cell_upper_check[8] = 11;
   map_yface2cell_upper_check[9] = 100;
   map_yface2cell_upper_check[10] = 99;
   map_yface2cell_upper_check[11] = 12;
   map_yface2cell_upper_check[12] = 15;
   map_yface2cell_upper_check[13] = 96;
   map_yface2cell_upper_check[14] = 95;
   map_yface2cell_upper_check[15] = 16;
   map_yface2cell_upper_check[16] = 17;
   map_yface2cell_upper_check[17] = 14;
   map_yface2cell_upper_check[18] = 19;
   map_yface2cell_upper_check[19] = 18;
   map_yface2cell_upper_check[20] = 20;
   map_yface2cell_upper_check[21] = 21;
   map_yface2cell_upper_check[22] = 9;
   map_yface2cell_upper_check[23] = 23;
   map_yface2cell_upper_check[24] = 22;
   map_yface2cell_upper_check[25] = 31;
   map_yface2cell_upper_check[26] = 50;
   map_yface2cell_upper_check[27] = 49;
   map_yface2cell_upper_check[28] = 47;
   map_yface2cell_upper_check[29] = 32;
   map_yface2cell_upper_check[30] = 41;
   map_yface2cell_upper_check[31] = 40;
   map_yface2cell_upper_check[32] = 72;
   map_yface2cell_upper_check[33] = 39;
   map_yface2cell_upper_check[34] = 42;
   map_yface2cell_upper_check[35] = 69;
   map_yface2cell_upper_check[36] = 68;
   map_yface2cell_upper_check[37] = 67;
   map_yface2cell_upper_check[38] = 44;
   map_yface2cell_upper_check[39] = 43;
   map_yface2cell_upper_check[40] = 45;
   map_yface2cell_upper_check[41] = 46;
   map_yface2cell_upper_check[42] = 53;
   map_yface2cell_upper_check[43] = 48;
   map_yface2cell_upper_check[44] = 51;
   map_yface2cell_upper_check[45] = 52;
   map_yface2cell_upper_check[46] = 55;
   map_yface2cell_upper_check[47] = 54;
   map_yface2cell_upper_check[48] = 57;
   map_yface2cell_upper_check[49] = 56;
   map_yface2cell_upper_check[50] = 59;
   map_yface2cell_upper_check[51] = 58;
   map_yface2cell_upper_check[52] = 63;
   map_yface2cell_upper_check[53] = 60;
   map_yface2cell_upper_check[54] = 61;
   map_yface2cell_upper_check[55] = 80;
   map_yface2cell_upper_check[56] = 80;
   map_yface2cell_upper_check[57] = 62;
   map_yface2cell_upper_check[58] = 79;
   map_yface2cell_upper_check[59] = 64;
   map_yface2cell_upper_check[60] = 64;
   map_yface2cell_upper_check[61] = 66;
   map_yface2cell_upper_check[62] = 65;
   map_yface2cell_upper_check[63] = 70;
   map_yface2cell_upper_check[64] = 75;
   map_yface2cell_upper_check[65] = 75;
   map_yface2cell_upper_check[66] = 71;
   map_yface2cell_upper_check[67] = 78;
   map_yface2cell_upper_check[68] = 81;
   map_yface2cell_upper_check[69] = 87;
   map_yface2cell_upper_check[70] = 86;
   map_yface2cell_upper_check[71] = 89;
   map_yface2cell_upper_check[72] = 89;
   map_yface2cell_upper_check[73] = 91;
   map_yface2cell_upper_check[74] = 90;
   map_yface2cell_upper_check[75] = 93;
   map_yface2cell_upper_check[76] = 94;
   map_yface2cell_upper_check[77] = 97;
   map_yface2cell_upper_check[78] = 92;
   map_yface2cell_upper_check[79] = 102;
   map_yface2cell_upper_check[80] = 98;
   map_yface2cell_upper_check[81] = 101;
   map_yface2cell_upper_check[82] = 102;
   map_yface2cell_upper_check[83] = 88;
   map_yface2cell_upper_check[84] = 110;
   map_yface2cell_upper_check[85] = 103;
   map_yface2cell_upper_check[86] = 106;
   map_yface2cell_upper_check[87] = 110;
}

int **set_xface_flag_check(int jsize, int isize, int lev)
{
   int **xface_flag_check = (int **)genmatrix(jsize, isize, sizeof(int));

   for (int jj = 0; jj < jsize; jj++){
      for (int ii = 0; ii < isize; ii++){
         xface_flag_check[jj][ii] = -1;
      }
   }

   if (lev == 1){
      xface_flag_check[7][1] = 1;
      xface_flag_check[7][2] = 1;
      xface_flag_check[7][3] = 1;
      xface_flag_check[7][4] = 1;
      xface_flag_check[7][5] = 1;
      xface_flag_check[6][1] = 1;
      xface_flag_check[6][2] = 1;
      xface_flag_check[6][3] = 1;
      xface_flag_check[6][4] = 1;
      xface_flag_check[6][5] = 1;
      xface_flag_check[5][0] = 1;
      xface_flag_check[5][1] = 1;
      xface_flag_check[5][5] = 1;
      xface_flag_check[5][6] = 1;
      xface_flag_check[4][0] = 1;
      xface_flag_check[4][6] = 1;
      xface_flag_check[3][0] = 1;
      xface_flag_check[3][6] = 1;
      xface_flag_check[2][0] = 1;
      xface_flag_check[2][1] = 1;
      xface_flag_check[2][5] = 1;
      xface_flag_check[2][6] = 1;
      xface_flag_check[1][1] = 1;
      xface_flag_check[1][2] = 1;
      xface_flag_check[1][3] = 1;
      xface_flag_check[1][4] = 1;
      xface_flag_check[1][5] = 1;
      xface_flag_check[0][1] = 1;
      xface_flag_check[0][2] = 1;
      xface_flag_check[0][3] = 1;
      xface_flag_check[0][4] = 1;
      xface_flag_check[0][5] = 1;
   }
   if (lev == 2){
      xface_flag_check[7][2] = 1;
      xface_flag_check[7][3] = 1;
      xface_flag_check[7][4] = 1;
      xface_flag_check[7][5] = 1;
      xface_flag_check[7][6] = 1;
      xface_flag_check[6][2] = 1;
      xface_flag_check[6][3] = 1;
      xface_flag_check[6][4] = 1;
      xface_flag_check[6][5] = 1;
      xface_flag_check[6][6] = 1;
      xface_flag_check[5][0] = 1;
      xface_flag_check[5][1] = 1;
      xface_flag_check[5][2] = 1;
      xface_flag_check[5][3] = 1;
      xface_flag_check[5][4] = 1;
      xface_flag_check[5][5] = 1;
      xface_flag_check[5][6] = 1;
      xface_flag_check[5][7] = 1;
      xface_flag_check[5][8] = 1;
      xface_flag_check[4][0] = 1;
      xface_flag_check[4][1] = 1;
      xface_flag_check[4][2] = 1;
      xface_flag_check[4][3] = 1;
      xface_flag_check[4][4] = 1;
      xface_flag_check[4][5] = 1;
      xface_flag_check[4][6] = 1;
      xface_flag_check[4][7] = 1;
      xface_flag_check[4][8] = 1;
      xface_flag_check[3][0] = 1;
      xface_flag_check[3][1] = 1;
      xface_flag_check[3][2] = 1;
      xface_flag_check[3][3] = 1;
      xface_flag_check[3][4] = 1;
      xface_flag_check[3][5] = 1;
      xface_flag_check[3][6] = 1;
      xface_flag_check[3][7] = 1;
      xface_flag_check[3][8] = 1;
      xface_flag_check[2][0] = 1;
      xface_flag_check[2][1] = 1;
      xface_flag_check[2][2] = 1;
      xface_flag_check[2][3] = 1;
      xface_flag_check[2][4] = 1;
      xface_flag_check[2][5] = 1;
      xface_flag_check[2][6] = 1;
      xface_flag_check[2][7] = 1;
      xface_flag_check[2][8] = 1;
      xface_flag_check[1][2] = 1;
      xface_flag_check[1][3] = 1;
      xface_flag_check[1][4] = 1;
      xface_flag_check[1][5] = 1;
      xface_flag_check[1][6] = 1;
      xface_flag_check[0][2] = 1;
      xface_flag_check[0][3] = 1;
      xface_flag_check[0][4] = 1;
      xface_flag_check[0][5] = 1;
      xface_flag_check[0][6] = 1;
   }

   return(xface_flag_check);
}

int **set_yface_flag_check(int jsize, int isize, int lev)
{
   int **yface_flag_check = (int **)genmatrix(jsize, isize, sizeof(int));

   for (int jj = 0; jj < jsize; jj++){
      for (int ii = 0; ii < isize; ii++){
         yface_flag_check[jj][ii] = -1;
      }
   }

   if (lev == 1){
      yface_flag_check[6][2] = 1;
      yface_flag_check[6][3] = 1;
      yface_flag_check[6][4] = 1;
      yface_flag_check[6][5] = 1;
      yface_flag_check[5][0] = 1;
      yface_flag_check[5][1] = 1;
      yface_flag_check[5][2] = 1;
      yface_flag_check[5][5] = 1;
      yface_flag_check[5][6] = 1;
      yface_flag_check[5][7] = 1;
      yface_flag_check[4][0] = 1;
      yface_flag_check[4][1] = 1;
      yface_flag_check[4][6] = 1;
      yface_flag_check[4][7] = 1;
      yface_flag_check[3][0] = 1;
      yface_flag_check[3][1] = 1;
      yface_flag_check[3][6] = 1;
      yface_flag_check[3][7] = 1;
      yface_flag_check[2][0] = 1;
      yface_flag_check[2][1] = 1;
      yface_flag_check[2][6] = 1;
      yface_flag_check[2][7] = 1;
      yface_flag_check[1][0] = 1;
      yface_flag_check[1][1] = 1;
      yface_flag_check[1][2] = 1;
      yface_flag_check[1][5] = 1;
      yface_flag_check[1][6] = 1;
      yface_flag_check[1][7] = 1;
      yface_flag_check[0][2] = 1;
      yface_flag_check[0][3] = 1;
      yface_flag_check[0][4] = 1;
      yface_flag_check[0][5] = 1;
   } else if (lev == 2) {
      yface_flag_check[8][2] = 1;
      yface_flag_check[8][3] = 1;
      yface_flag_check[8][4] = 1;
      yface_flag_check[8][5] = 1;
      yface_flag_check[7][2] = 1;
      yface_flag_check[7][3] = 1;
      yface_flag_check[7][4] = 1;
      yface_flag_check[7][5] = 1;
      yface_flag_check[6][0] = 1;
      yface_flag_check[6][1] = 1;
      yface_flag_check[6][2] = 1;
      yface_flag_check[6][3] = 1;
      yface_flag_check[6][4] = 1;
      yface_flag_check[6][5] = 1;
      yface_flag_check[6][6] = 1;
      yface_flag_check[6][7] = 1;
      yface_flag_check[5][0] = 1;
      yface_flag_check[5][1] = 1;
      yface_flag_check[5][2] = 1;
      yface_flag_check[5][3] = 1;
      yface_flag_check[5][4] = 1;
      yface_flag_check[5][5] = 1;
      yface_flag_check[5][6] = 1;
      yface_flag_check[5][7] = 1;
      yface_flag_check[4][0] = 1;
      yface_flag_check[4][1] = 1;
      yface_flag_check[4][2] = 1;
      yface_flag_check[4][3] = 1;
      yface_flag_check[4][4] = 1;
      yface_flag_check[4][5] = 1;
      yface_flag_check[4][6] = 1;
      yface_flag_check[4][7] = 1;
      yface_flag_check[3][0] = 1;
      yface_flag_check[3][1] = 1;
      yface_flag_check[3][2] = 1;
      yface_flag_check[3][3] = 1;
      yface_flag_check[3][4] = 1;
      yface_flag_check[3][5] = 1;
      yface_flag_check[3][6] = 1;
      yface_flag_check[3][7] = 1;
      yface_flag_check[2][0] = 1;
      yface_flag_check[2][1] = 1;
      yface_flag_check[2][2] = 1;
      yface_flag_check[2][3] = 1;
      yface_flag_check[2][4] = 1;
      yface_flag_check[2][5] = 1;
      yface_flag_check[2][6] = 1;
      yface_flag_check[2][7] = 1;
      yface_flag_check[1][2] = 1;
      yface_flag_check[1][3] = 1;
      yface_flag_check[1][4] = 1;
      yface_flag_check[1][5] = 1;
      yface_flag_check[0][2] = 1;
      yface_flag_check[0][3] = 1;
      yface_flag_check[0][4] = 1;
      yface_flag_check[0][5] = 1;
   }

   return(yface_flag_check);
}

int **set_zone_flag_check(int jsize, int isize, int lev)
{
   int **zone_flag_check = (int **)genmatrix(jsize, isize, sizeof(int));

   for (int jj = 0; jj < jsize; jj++){
      for (int ii = 0; ii < isize; ii++){
         zone_flag_check[jj][ii] = -1;
      }
   }

   if (lev == 1){
      zone_flag_check[8][2] = 1;
      zone_flag_check[8][3] = 1;
      zone_flag_check[8][4] = 1;
      zone_flag_check[8][5] = 1;
      zone_flag_check[8][6] = 1;
      zone_flag_check[8][7] = 1;
      zone_flag_check[7][1] = 1;
      zone_flag_check[7][2] = 1;
      zone_flag_check[7][3] = 1;
      zone_flag_check[7][4] = 1;
      zone_flag_check[7][5] = 1;
      zone_flag_check[7][6] = 1;
      zone_flag_check[7][7] = 1;
      zone_flag_check[7][8] = 1;
      zone_flag_check[6][1] = 1;
      zone_flag_check[6][2] = 1;
      zone_flag_check[6][3] = 1;
      zone_flag_check[6][6] = 1;
      zone_flag_check[6][7] = 1;
      zone_flag_check[6][8] = 1;
      zone_flag_check[5][1] = 1;
      zone_flag_check[5][2] = 1;
      zone_flag_check[5][7] = 1;
      zone_flag_check[5][8] = 1;
      zone_flag_check[4][1] = 1;
      zone_flag_check[4][2] = 1;
      zone_flag_check[4][7] = 1;
      zone_flag_check[4][8] = 1;
      zone_flag_check[3][1] = 1;
      zone_flag_check[3][2] = 1;
      zone_flag_check[3][3] = 1;
      zone_flag_check[3][6] = 1;
      zone_flag_check[3][7] = 1;
      zone_flag_check[3][8] = 1;
      zone_flag_check[2][1] = 1;
      zone_flag_check[2][2] = 1;
      zone_flag_check[2][3] = 1;
      zone_flag_check[2][4] = 1;
      zone_flag_check[2][5] = 1;
      zone_flag_check[2][6] = 1;
      zone_flag_check[2][7] = 1;
      zone_flag_check[2][8] = 1;
      zone_flag_check[1][2] = 1;
      zone_flag_check[1][3] = 1;
      zone_flag_check[1][4] = 1;
      zone_flag_check[1][5] = 1;
      zone_flag_check[1][6] = 1;
      zone_flag_check[1][7] = 1;
   } else if (lev == 2){
      zone_flag_check[10][4] = 1;
      zone_flag_check[10][5] = 1;
      zone_flag_check[10][6] = 1;
      zone_flag_check[10][7] = 1;
      zone_flag_check[9][3] = 1;
      zone_flag_check[9][4] = 1;
      zone_flag_check[9][5] = 1;
      zone_flag_check[9][6] = 1;
      zone_flag_check[9][7] = 1;
      zone_flag_check[9][8] = 1;
      zone_flag_check[8][2] = 1;
      zone_flag_check[8][3] = 1;
      zone_flag_check[8][4] = 1;
      zone_flag_check[8][5] = 1;
      zone_flag_check[8][6] = 1;
      zone_flag_check[8][7] = 1;
      zone_flag_check[8][8] = 1;
      zone_flag_check[8][9] = 1;
      zone_flag_check[7][1] = 1;
      zone_flag_check[7][2] = 1;
      zone_flag_check[7][3] = 1;
      zone_flag_check[7][4] = 1;
      zone_flag_check[7][5] = 1;
      zone_flag_check[7][6] = 1;
      zone_flag_check[7][7] = 1;
      zone_flag_check[7][8] = 1;
      zone_flag_check[7][9] = 1;
      zone_flag_check[7][10] = 1;
      zone_flag_check[6][1] = 1;
      zone_flag_check[6][2] = 1;
      zone_flag_check[6][3] = 1;
      zone_flag_check[6][4] = 1;
      zone_flag_check[6][5] = 1;
      zone_flag_check[6][6] = 1;
      zone_flag_check[6][7] = 1;
      zone_flag_check[6][8] = 1;
      zone_flag_check[6][9] = 1;
      zone_flag_check[6][10] = 1;
      zone_flag_check[5][1] = 1;
      zone_flag_check[5][2] = 1;
      zone_flag_check[5][3] = 1;
      zone_flag_check[5][4] = 1;
      zone_flag_check[5][5] = 1;
      zone_flag_check[5][6] = 1;
      zone_flag_check[5][7] = 1;
      zone_flag_check[5][8] = 1;
      zone_flag_check[5][9] = 1;
      zone_flag_check[5][10] = 1;
      zone_flag_check[4][1] = 1;
      zone_flag_check[4][2] = 1;
      zone_flag_check[4][3] = 1;
      zone_flag_check[4][4] = 1;
      zone_flag_check[4][5] = 1;
      zone_flag_check[4][6] = 1;
      zone_flag_check[4][7] = 1;
      zone_flag_check[4][8] = 1;
      zone_flag_check[4][9] = 1;
      zone_flag_check[4][10] = 1;
      zone_flag_check[3][2] = 1;
      zone_flag_check[3][3] = 1;
      zone_flag_check[3][4] = 1;
      zone_flag_check[3][5] = 1;
      zone_flag_check[3][6] = 1;
      zone_flag_check[3][7] = 1;
      zone_flag_check[3][8] = 1;
      zone_flag_check[3][9] = 1;
      zone_flag_check[2][3] = 1;
      zone_flag_check[2][4] = 1;
      zone_flag_check[2][5] = 1;
      zone_flag_check[2][6] = 1;
      zone_flag_check[2][7] = 1;
      zone_flag_check[2][8] = 1;
      zone_flag_check[1][4] = 1;
      zone_flag_check[1][5] = 1;
      zone_flag_check[1][6] = 1;
      zone_flag_check[1][7] = 1;
   }

   return(zone_flag_check);
}

int **set_zone_cell_check(int jsize, int isize, int lev)
{
   int **zone_cell_check = (int **)genmatrix(jsize, isize, sizeof(int));

   for (int jj = 0; jj < jsize; jj++){
      for (int ii = 0; ii < isize; ii++){
         zone_cell_check[jj][ii] = -1;
      }
   }

   if (lev == 1){
      zone_cell_check[9][3] = 85;
      zone_cell_check[9][4] = 84;
      zone_cell_check[9][5] = 83;
      zone_cell_check[9][6] = 82;
      zone_cell_check[8][1] = 110;
      zone_cell_check[8][2] = 110;
      zone_cell_check[8][3] = 87;
      zone_cell_check[8][4] = 86;
      zone_cell_check[8][5] = 81;
      zone_cell_check[8][6] = 78;
      zone_cell_check[8][7] = 75;
      zone_cell_check[8][8] = 75;
      zone_cell_check[7][1] = 110;
      zone_cell_check[7][2] = 110;
      zone_cell_check[7][3] = 88;
      zone_cell_check[7][4] = 89;
      zone_cell_check[7][5] = 80;
      zone_cell_check[7][6] = 79;
      zone_cell_check[7][7] = 75;
      zone_cell_check[7][8] = 75;
      zone_cell_check[6][0] = 108;
      zone_cell_check[6][1] = 106;
      zone_cell_check[6][2] = 103;
      zone_cell_check[6][3] = 102;
      zone_cell_check[6][4] = 91;
      zone_cell_check[6][5] = 61;
      zone_cell_check[6][6] = 64;
      zone_cell_check[6][7] = 70;
      zone_cell_check[6][8] = 71;
      zone_cell_check[6][9] = 74;
      zone_cell_check[5][0] = 107;
      zone_cell_check[5][1] = 105;
      zone_cell_check[5][2] = 104;
      zone_cell_check[5][3] = 101;
      zone_cell_check[5][6] = 66;
      zone_cell_check[5][7] = 69;
      zone_cell_check[5][8] = 72;
      zone_cell_check[5][9] = 73;
      zone_cell_check[4][0] = 4;
      zone_cell_check[4][1] = 6;
      zone_cell_check[4][2] = 7;
      zone_cell_check[4][3] = 10;
      zone_cell_check[4][6] = 45;
      zone_cell_check[4][7] = 42;
      zone_cell_check[4][8] = 39;
      zone_cell_check[4][9] = 38;
      zone_cell_check[3][0] = 3;
      zone_cell_check[3][1] = 5;
      zone_cell_check[3][2] = 8;
      zone_cell_check[3][3] = 9;
      zone_cell_check[3][4] = 20;
      zone_cell_check[3][5] = 50;
      zone_cell_check[3][6] = 47;
      zone_cell_check[3][7] = 41;
      zone_cell_check[3][8] = 40;
      zone_cell_check[3][9] = 37;
      zone_cell_check[2][1] = 1;
      zone_cell_check[2][2] = 1;
      zone_cell_check[2][3] = 23;
      zone_cell_check[2][4] = 22;
      zone_cell_check[2][5] = 31;
      zone_cell_check[2][6] = 32;
      zone_cell_check[2][7] = 36;
      zone_cell_check[2][8] = 36;
      zone_cell_check[1][1] = 1;
      zone_cell_check[1][2] = 1;
      zone_cell_check[1][3] = 24;
      zone_cell_check[1][4] = 25;
      zone_cell_check[1][5] = 30;
      zone_cell_check[1][6] = 33;
      zone_cell_check[1][7] = 36;
      zone_cell_check[1][8] = 36;
      zone_cell_check[0][3] = 26;
      zone_cell_check[0][4] = 27;
      zone_cell_check[0][5] = 28;
      zone_cell_check[0][6] = 29;
   } else if (lev == 2) {
      zone_cell_check[11][4] = 89;
      zone_cell_check[11][5] = 89;
      zone_cell_check[11][6] = 80;
      zone_cell_check[11][7] = 80;
      zone_cell_check[10][4] = 89;
      zone_cell_check[10][5] = 89;
      zone_cell_check[10][6] = 80;
      zone_cell_check[10][7] = 80;
      zone_cell_check[9][2] = 102;
      zone_cell_check[9][3] = 102;
      zone_cell_check[9][4] = 91;
      zone_cell_check[9][5] = 90;
      zone_cell_check[9][6] = 61;
      zone_cell_check[9][7] = 62;
      zone_cell_check[9][8] = 64;
      zone_cell_check[9][9] = 64;
      zone_cell_check[8][2] = 102;
      zone_cell_check[8][3] = 102;
      zone_cell_check[8][4] = 92;
      zone_cell_check[8][5] = 93;
      zone_cell_check[8][6] = 60;
      zone_cell_check[8][7] = 63;
      zone_cell_check[8][8] = 64;
      zone_cell_check[8][9] = 64;
      zone_cell_check[7][0] = 104;
      zone_cell_check[7][1] = 104;
      zone_cell_check[7][2] = 101;
      zone_cell_check[7][3] = 98;
      zone_cell_check[7][4] = 97;
      zone_cell_check[7][5] = 94;
      zone_cell_check[7][6] = 59;
      zone_cell_check[7][7] = 58;
      zone_cell_check[7][8] = 66;
      zone_cell_check[7][9] = 65;
      zone_cell_check[7][10] = 69;
      zone_cell_check[7][11] = 69;
      zone_cell_check[6][0] = 104;
      zone_cell_check[6][1] = 104;
      zone_cell_check[6][2] = 100;
      zone_cell_check[6][3] = 99;
      zone_cell_check[6][4] = 96;
      zone_cell_check[6][5] = 95;
      zone_cell_check[6][6] = 56;
      zone_cell_check[6][7] = 57;
      zone_cell_check[6][8] = 67;
      zone_cell_check[6][9] = 68;
      zone_cell_check[6][10] = 69;
      zone_cell_check[6][11] = 69;
      zone_cell_check[5][0] = 7;
      zone_cell_check[5][1] = 7;
      zone_cell_check[5][2] = 11;
      zone_cell_check[5][3] = 12;
      zone_cell_check[5][4] = 15;
      zone_cell_check[5][5] = 16;
      zone_cell_check[5][6] = 55;
      zone_cell_check[5][7] = 54;
      zone_cell_check[5][8] = 44;
      zone_cell_check[5][9] = 43;
      zone_cell_check[5][10] = 42;
      zone_cell_check[5][11] = 42;
      zone_cell_check[4][0] = 7;
      zone_cell_check[4][1] = 7;
      zone_cell_check[4][2] = 10;
      zone_cell_check[4][3] = 13;
      zone_cell_check[4][4] = 14;
      zone_cell_check[4][5] = 17;
      zone_cell_check[4][6] = 52;
      zone_cell_check[4][7] = 53;
      zone_cell_check[4][8] = 45;
      zone_cell_check[4][9] = 46;
      zone_cell_check[4][10] = 42;
      zone_cell_check[4][11] = 42;
      zone_cell_check[3][2] = 9;
      zone_cell_check[3][3] = 9;
      zone_cell_check[3][4] = 19;
      zone_cell_check[3][5] = 18;
      zone_cell_check[3][6] = 51;
      zone_cell_check[3][7] = 48;
      zone_cell_check[3][8] = 47;
      zone_cell_check[3][9] = 47;
      zone_cell_check[2][2] = 9;
      zone_cell_check[2][3] = 9;
      zone_cell_check[2][4] = 20;
      zone_cell_check[2][5] = 21;
      zone_cell_check[2][6] = 50;
      zone_cell_check[2][7] = 49;
      zone_cell_check[2][8] = 47;
      zone_cell_check[2][9] = 47;
      zone_cell_check[1][4] = 22;
      zone_cell_check[1][5] = 22;
      zone_cell_check[1][6] = 31;
      zone_cell_check[1][7] = 31;
      zone_cell_check[0][4] = 22;
      zone_cell_check[0][5] = 22;
      zone_cell_check[0][6] = 31;
      zone_cell_check[0][7] = 31;
   }

   return(zone_cell_check);
}

int *set_map_xcell2face_left1_check(int ncells)
{
   int *map_xcell2face_left1_check = (int *)malloc(ncells*sizeof(int));

   map_xcell2face_left1_check[0] = -1;
   map_xcell2face_left1_check[1] = -1;
   map_xcell2face_left1_check[2] = -1;
   map_xcell2face_left1_check[3] = -1;
   map_xcell2face_left1_check[4] = -1;
   map_xcell2face_left1_check[5] = -1;
   map_xcell2face_left1_check[6] = -1;
   map_xcell2face_left1_check[7] = 3;
   map_xcell2face_left1_check[8] = 2;
   map_xcell2face_left1_check[9] = 6;
   map_xcell2face_left1_check[10] = 4;
   map_xcell2face_left1_check[11] = 5;
   map_xcell2face_left1_check[12] = 10;
   map_xcell2face_left1_check[13] = 9;
   map_xcell2face_left1_check[14] = 12;
   map_xcell2face_left1_check[15] = 11;
   map_xcell2face_left1_check[16] = 14;
   map_xcell2face_left1_check[17] = 13;
   map_xcell2face_left1_check[18] = 18;
   map_xcell2face_left1_check[19] = 8;
   map_xcell2face_left1_check[20] = 7;
   map_xcell2face_left1_check[21] = 19;
   map_xcell2face_left1_check[22] = 22;
   map_xcell2face_left1_check[23] = 1;
   map_xcell2face_left1_check[25] = 23;
   map_xcell2face_left1_check[26] = -1;
   map_xcell2face_left1_check[27] = -1;
   map_xcell2face_left1_check[28] = -1;
   map_xcell2face_left1_check[29] = -1;
   map_xcell2face_left1_check[30] = 24;
   map_xcell2face_left1_check[31] = 21;
   map_xcell2face_left1_check[32] = 26;
   map_xcell2face_left1_check[33] = 25;
   map_xcell2face_left1_check[34] = -1;
   map_xcell2face_left1_check[35] = -1;
   map_xcell2face_left1_check[36] = 28;
   map_xcell2face_left1_check[37] = -1;
   map_xcell2face_left1_check[38] = -1;
   map_xcell2face_left1_check[39] = 30;
   map_xcell2face_left1_check[40] = 29;
   map_xcell2face_left1_check[41] = 35;
   map_xcell2face_left1_check[42] = 34;
   map_xcell2face_left1_check[43] = 32;
   map_xcell2face_left1_check[44] = 42;
   map_xcell2face_left1_check[45] = 41;
   map_xcell2face_left1_check[46] = 33;
   map_xcell2face_left1_check[47] = 37;
   map_xcell2face_left1_check[48] = 39;
   map_xcell2face_left1_check[49] = 38;
   map_xcell2face_left1_check[50] = 20;
   map_xcell2face_left1_check[51] = 17;
   map_xcell2face_left1_check[52] = 16;
   map_xcell2face_left1_check[53] = 40;
   map_xcell2face_left1_check[54] = 43;
   map_xcell2face_left1_check[55] = 15;
   map_xcell2face_left1_check[56] = 72;
   map_xcell2face_left1_check[57] = 44;
   map_xcell2face_left1_check[58] = 47;
   map_xcell2face_left1_check[59] = 71;
   map_xcell2face_left1_check[60] = 70;
   map_xcell2face_left1_check[61] = 67;
   map_xcell2face_left1_check[62] = 49;
   map_xcell2face_left1_check[63] = 48;
   map_xcell2face_left1_check[64] = 51;
   map_xcell2face_left1_check[65] = 54;
   map_xcell2face_left1_check[66] = 46;
   map_xcell2face_left1_check[67] = 45;
   map_xcell2face_left1_check[68] = 55;
   map_xcell2face_left1_check[69] = 56;
   map_xcell2face_left1_check[70] = 52;
   map_xcell2face_left1_check[71] = 58;
   map_xcell2face_left1_check[72] = 57;
   map_xcell2face_left1_check[73] = -1;
   map_xcell2face_left1_check[74] = -1;
   map_xcell2face_left1_check[75] = 60;
   map_xcell2face_left1_check[76] = -1;
   map_xcell2face_left1_check[77] = -1;
   map_xcell2face_left1_check[78] = 62;
   map_xcell2face_left1_check[79] = 61;
   map_xcell2face_left1_check[80] = 66;
   map_xcell2face_left1_check[81] = 63;
   map_xcell2face_left1_check[82] = -1;
   map_xcell2face_left1_check[83] = -1;
   map_xcell2face_left1_check[84] = -1;
   map_xcell2face_left1_check[85] = -1;
   map_xcell2face_left1_check[86] = 64;
   map_xcell2face_left1_check[87] = 87;
   map_xcell2face_left1_check[88] = 86;
   map_xcell2face_left1_check[89] = 65;
   map_xcell2face_left1_check[90] = 68;
   map_xcell2face_left1_check[91] = 80;
   map_xcell2face_left1_check[92] = 79;
   map_xcell2face_left1_check[93] = 69;
   map_xcell2face_left1_check[94] = 74;
   map_xcell2face_left1_check[95] = 73;
   map_xcell2face_left1_check[96] = 76;
   map_xcell2face_left1_check[97] = 75;
   map_xcell2face_left1_check[98] = 78;
   map_xcell2face_left1_check[99] = 77;
   map_xcell2face_left1_check[100] = 82;
   map_xcell2face_left1_check[101] = 83;
   map_xcell2face_left1_check[102] = 81;
   map_xcell2face_left1_check[103] = 85;
   map_xcell2face_left1_check[104] = 84;
   map_xcell2face_left1_check[105] = -1;
   map_xcell2face_left1_check[106] = -1;
   map_xcell2face_left1_check[107] = -1;
   map_xcell2face_left1_check[108] = -1;
   map_xcell2face_left1_check[109] = -1;
   map_xcell2face_left1_check[110] = -1;
   map_xcell2face_left1_check[111] = -1;

   return(map_xcell2face_left1_check);
}

int *set_map_xcell2face_left2_check(int ncells)
{
   int *map_xcell2face_left2_check = (int *)malloc(ncells*sizeof(int));

   map_xcell2face_left2_check[0] = -1;
   map_xcell2face_left2_check[1] = -1;
   map_xcell2face_left2_check[2] = -1;
   map_xcell2face_left2_check[3] = -1;
   map_xcell2face_left2_check[4] = -1;
   map_xcell2face_left2_check[5] = -1;
   map_xcell2face_left2_check[6] = -1;
   map_xcell2face_left2_check[7] = -1;
   map_xcell2face_left2_check[8] = -1;
   map_xcell2face_left2_check[9] = -1;
   map_xcell2face_left2_check[10] = -1;
   map_xcell2face_left2_check[11] = -1;
   map_xcell2face_left2_check[12] = -1;
   map_xcell2face_left2_check[13] = -1;
   map_xcell2face_left2_check[14] = -1;
   map_xcell2face_left2_check[15] = -1;
   map_xcell2face_left2_check[16] = -1;
   map_xcell2face_left2_check[17] = -1;
   map_xcell2face_left2_check[18] = -1;
   map_xcell2face_left2_check[19] = -1;
   map_xcell2face_left2_check[20] = -1;
   map_xcell2face_left2_check[21] = -1;
   map_xcell2face_left2_check[22] = -1;
   map_xcell2face_left2_check[23] = -1;
   map_xcell2face_left2_check[24] = -1;
   map_xcell2face_left2_check[25] = -1;
   map_xcell2face_left2_check[26] = -1;
   map_xcell2face_left2_check[27] = -1;
   map_xcell2face_left2_check[28] = -1;
   map_xcell2face_left2_check[29] = -1;
   map_xcell2face_left2_check[30] = -1;
   map_xcell2face_left2_check[31] = -1;
   map_xcell2face_left2_check[32] = -1;
   map_xcell2face_left2_check[33] = -1;
   map_xcell2face_left2_check[34] = -1;
   map_xcell2face_left2_check[35] = -1;
   map_xcell2face_left2_check[36] = 27;
   map_xcell2face_left2_check[37] = -1;
   map_xcell2face_left2_check[38] = -1;
   map_xcell2face_left2_check[39] = -1;
   map_xcell2face_left2_check[40] = -1;
   map_xcell2face_left2_check[41] = -1;
   map_xcell2face_left2_check[42] = 31;
   map_xcell2face_left2_check[43] = -1;
   map_xcell2face_left2_check[44] = -1;
   map_xcell2face_left2_check[45] = -1;
   map_xcell2face_left2_check[46] = -1;
   map_xcell2face_left2_check[47] = 36;
   map_xcell2face_left2_check[48] = -1;
   map_xcell2face_left2_check[49] = -1;
   map_xcell2face_left2_check[50] = -1;
   map_xcell2face_left2_check[51] = -1;
   map_xcell2face_left2_check[52] = -1;
   map_xcell2face_left2_check[53] = -1;
   map_xcell2face_left2_check[54] = -1;
   map_xcell2face_left2_check[55] = -1;
   map_xcell2face_left2_check[56] = -1;
   map_xcell2face_left2_check[57] = -1;
   map_xcell2face_left2_check[58] = -1;
   map_xcell2face_left2_check[59] = -1;
   map_xcell2face_left2_check[60] = -1;
   map_xcell2face_left2_check[61] = -1;
   map_xcell2face_left2_check[62] = -1;
   map_xcell2face_left2_check[63] = -1;
   map_xcell2face_left2_check[64] = 50;
   map_xcell2face_left2_check[65] = -1;
   map_xcell2face_left2_check[66] = -1;
   map_xcell2face_left2_check[67] = -1;
   map_xcell2face_left2_check[68] = -1;
   map_xcell2face_left2_check[69] = 53;
   map_xcell2face_left2_check[70] = -1;
   map_xcell2face_left2_check[71] = -1;
   map_xcell2face_left2_check[72] = -1;
   map_xcell2face_left2_check[73] = -1;
   map_xcell2face_left2_check[74] = -1;
   map_xcell2face_left2_check[75] = 59;
   map_xcell2face_left2_check[76] = -1;
   map_xcell2face_left2_check[77] = -1;
   map_xcell2face_left2_check[78] = -1;
   map_xcell2face_left2_check[79] = -1;
   map_xcell2face_left2_check[80] = -1;
   map_xcell2face_left2_check[81] = -1;
   map_xcell2face_left2_check[82] = -1;
   map_xcell2face_left2_check[83] = -1;
   map_xcell2face_left2_check[84] = -1;
   map_xcell2face_left2_check[85] = -1;
   map_xcell2face_left2_check[86] = -1;
   map_xcell2face_left2_check[87] = -1;
   map_xcell2face_left2_check[88] = -1;
   map_xcell2face_left2_check[89] = -1;
   map_xcell2face_left2_check[90] = -1;
   map_xcell2face_left2_check[91] = -1;
   map_xcell2face_left2_check[92] = -1;
   map_xcell2face_left2_check[93] = -1;
   map_xcell2face_left2_check[94] = -1;
   map_xcell2face_left2_check[95] = -1;
   map_xcell2face_left2_check[96] = -1;
   map_xcell2face_left2_check[97] = -1;
   map_xcell2face_left2_check[98] = -1;
   map_xcell2face_left2_check[99] = -1;
   map_xcell2face_left2_check[100] = -1;
   map_xcell2face_left2_check[101] = -1;
   map_xcell2face_left2_check[102] = -1;
   map_xcell2face_left2_check[103] = -1;
   map_xcell2face_left2_check[104] = -1;
   map_xcell2face_left2_check[105] = -1;
   map_xcell2face_left2_check[106] = -1;
   map_xcell2face_left2_check[107] = -1;
   map_xcell2face_left2_check[108] = -1;
   map_xcell2face_left2_check[109] = -1;
   map_xcell2face_left2_check[110] = -1;
   map_xcell2face_left2_check[111] = -1;

   return(map_xcell2face_left2_check);
}

int *set_map_xcell2face_right1_check(int ncells)
{
   int *map_xcell2face_right1_check = (int *)malloc(ncells*sizeof(int));

   map_xcell2face_right1_check[0] = -1;
   map_xcell2face_right1_check[2] = -1;
   map_xcell2face_right1_check[3] = -1;
   map_xcell2face_right1_check[4] = -1;
   map_xcell2face_right1_check[5] = 2;
   map_xcell2face_right1_check[6] = 3;
   map_xcell2face_right1_check[7] = 4;
   map_xcell2face_right1_check[8] = 6;
   map_xcell2face_right1_check[9] = 7;
   map_xcell2face_right1_check[10] = 9;
   map_xcell2face_right1_check[11] = 10;
   map_xcell2face_right1_check[12] = 11;
   map_xcell2face_right1_check[13] = 12;
   map_xcell2face_right1_check[14] = 13;
   map_xcell2face_right1_check[15] = 14;
   map_xcell2face_right1_check[16] = 15;
   map_xcell2face_right1_check[17] = 16;
   map_xcell2face_right1_check[18] = 17;
   map_xcell2face_right1_check[19] = 18;
   map_xcell2face_right1_check[20] = 19;
   map_xcell2face_right1_check[21] = 20;
   map_xcell2face_right1_check[22] = 21;
   map_xcell2face_right1_check[23] = 22;
   map_xcell2face_right1_check[24] = 23;
   map_xcell2face_right1_check[25] = 24;
   map_xcell2face_right1_check[26] = -1;
   map_xcell2face_right1_check[27] = -1;
   map_xcell2face_right1_check[28] = -1;
   map_xcell2face_right1_check[29] = -1;
   map_xcell2face_right1_check[30] = 25;
   map_xcell2face_right1_check[31] = 26;
   map_xcell2face_right1_check[32] = 27;
   map_xcell2face_right1_check[33] = 28;
   map_xcell2face_right1_check[34] = -1;
   map_xcell2face_right1_check[35] = -1;
   map_xcell2face_right1_check[36] = -1;
   map_xcell2face_right1_check[37] = -1;
   map_xcell2face_right1_check[38] = -1;
   map_xcell2face_right1_check[39] = -1;
   map_xcell2face_right1_check[40] = -1;
   map_xcell2face_right1_check[41] = 29;
   map_xcell2face_right1_check[42] = 30;
   map_xcell2face_right1_check[43] = 31;
   map_xcell2face_right1_check[44] = 32;
   map_xcell2face_right1_check[45] = 33;
   map_xcell2face_right1_check[46] = 34;
   map_xcell2face_right1_check[47] = 35;
   map_xcell2face_right1_check[48] = 36;
   map_xcell2face_right1_check[49] = 37;
   map_xcell2face_right1_check[50] = 38;
   map_xcell2face_right1_check[51] = 39;
   map_xcell2face_right1_check[52] = 40;
   map_xcell2face_right1_check[53] = 41;
   map_xcell2face_right1_check[54] = 42;
   map_xcell2face_right1_check[55] = 43;
   map_xcell2face_right1_check[56] = 44;
   map_xcell2face_right1_check[57] = 45;
   map_xcell2face_right1_check[58] = 46;
   map_xcell2face_right1_check[59] = 47;
   map_xcell2face_right1_check[60] = 48;
   map_xcell2face_right1_check[61] = 49;
   map_xcell2face_right1_check[62] = 50;
   map_xcell2face_right1_check[63] = 51;
   map_xcell2face_right1_check[64] = 52;
   map_xcell2face_right1_check[65] = 53;
   map_xcell2face_right1_check[66] = 54;
   map_xcell2face_right1_check[67] = 55;
   map_xcell2face_right1_check[68] = 56;
   map_xcell2face_right1_check[69] = 57;
   map_xcell2face_right1_check[70] = 58;
   map_xcell2face_right1_check[71] = -1;
   map_xcell2face_right1_check[72] = -1;
   map_xcell2face_right1_check[73] = -1;
   map_xcell2face_right1_check[74] = -1;
   map_xcell2face_right1_check[75] = -1;
   map_xcell2face_right1_check[76] = -1;
   map_xcell2face_right1_check[77] = -1;
   map_xcell2face_right1_check[78] = 59;
   map_xcell2face_right1_check[79] = 60;
   map_xcell2face_right1_check[80] = 61;
   map_xcell2face_right1_check[81] = 62;
   map_xcell2face_right1_check[82] = -1;
   map_xcell2face_right1_check[83] = -1;
   map_xcell2face_right1_check[84] = -1;
   map_xcell2face_right1_check[85] = -1;
   map_xcell2face_right1_check[86] = 63;
   map_xcell2face_right1_check[87] = 64;
   map_xcell2face_right1_check[88] = 65;
   map_xcell2face_right1_check[89] = 66;
   map_xcell2face_right1_check[90] = 67;
   map_xcell2face_right1_check[91] = 68;
   map_xcell2face_right1_check[92] = 69;
   map_xcell2face_right1_check[93] = 70;
   map_xcell2face_right1_check[94] = 71;
   map_xcell2face_right1_check[95] = 72;
   map_xcell2face_right1_check[96] = 73;
   map_xcell2face_right1_check[97] = 74;
   map_xcell2face_right1_check[98] = 75;
   map_xcell2face_right1_check[99] = 76;
   map_xcell2face_right1_check[100] = 77;
   map_xcell2face_right1_check[101] = 78;
   map_xcell2face_right1_check[102] = 79;
   map_xcell2face_right1_check[103] = 81;
   map_xcell2face_right1_check[104] = 82;
   map_xcell2face_right1_check[105] = 84;
   map_xcell2face_right1_check[106] = 85;
   map_xcell2face_right1_check[107] = -1;
   map_xcell2face_right1_check[108] = -1;
   map_xcell2face_right1_check[109] = -1;
   map_xcell2face_right1_check[110] = 86;
   map_xcell2face_right1_check[111] = -1;

   return(map_xcell2face_right1_check);
}

int *set_map_xcell2face_right2_check(int ncells)
{
   int *map_xcell2face_right2_check = (int *)malloc(ncells*sizeof(int));

   map_xcell2face_right2_check[0] = -1;
   map_xcell2face_right2_check[1] = 1;
   map_xcell2face_right2_check[2] = -1;
   map_xcell2face_right2_check[3] = -1;
   map_xcell2face_right2_check[4] = -1;
   map_xcell2face_right2_check[5] = -1;
   map_xcell2face_right2_check[6] = -1;
   map_xcell2face_right2_check[7] = 5;
   map_xcell2face_right2_check[8] = -1;
   map_xcell2face_right2_check[9] = 8;
   map_xcell2face_right2_check[10] = -1;
   map_xcell2face_right2_check[11] = -1;
   map_xcell2face_right2_check[12] = -1;
   map_xcell2face_right2_check[13] = -1;
   map_xcell2face_right2_check[14] = -1;
   map_xcell2face_right2_check[15] = -1;
   map_xcell2face_right2_check[16] = -1;
   map_xcell2face_right2_check[17] = -1;
   map_xcell2face_right2_check[18] = -1;
   map_xcell2face_right2_check[19] = -1;
   map_xcell2face_right2_check[20] = -1;
   map_xcell2face_right2_check[21] = -1;
   map_xcell2face_right2_check[22] = -1;
   map_xcell2face_right2_check[23] = -1;
   map_xcell2face_right2_check[24] = -1;
   map_xcell2face_right2_check[25] = -1;
   map_xcell2face_right2_check[26] = -1;
   map_xcell2face_right2_check[27] = -1;
   map_xcell2face_right2_check[28] = -1;
   map_xcell2face_right2_check[29] = -1;
   map_xcell2face_right2_check[30] = -1;
   map_xcell2face_right2_check[31] = -1;
   map_xcell2face_right2_check[32] = -1;
   map_xcell2face_right2_check[33] = -1;
   map_xcell2face_right2_check[34] = -1;
   map_xcell2face_right2_check[35] = -1;
   map_xcell2face_right2_check[36] = -1;
   map_xcell2face_right2_check[37] = -1;
   map_xcell2face_right2_check[38] = -1;
   map_xcell2face_right2_check[39] = -1;
   map_xcell2face_right2_check[40] = -1;
   map_xcell2face_right2_check[41] = -1;
   map_xcell2face_right2_check[42] = -1;
   map_xcell2face_right2_check[43] = -1;
   map_xcell2face_right2_check[44] = -1;
   map_xcell2face_right2_check[45] = -1;
   map_xcell2face_right2_check[46] = -1;
   map_xcell2face_right2_check[47] = -1;
   map_xcell2face_right2_check[48] = -1;
   map_xcell2face_right2_check[49] = -1;
   map_xcell2face_right2_check[50] = -1;
   map_xcell2face_right2_check[51] = -1;
   map_xcell2face_right2_check[52] = -1;
   map_xcell2face_right2_check[53] = -1;
   map_xcell2face_right2_check[54] = -1;
   map_xcell2face_right2_check[55] = -1;
   map_xcell2face_right2_check[56] = -1;
   map_xcell2face_right2_check[57] = -1;
   map_xcell2face_right2_check[58] = -1;
   map_xcell2face_right2_check[59] = -1;
   map_xcell2face_right2_check[60] = -1;
   map_xcell2face_right2_check[61] = -1;
   map_xcell2face_right2_check[62] = -1;
   map_xcell2face_right2_check[63] = -1;
   map_xcell2face_right2_check[64] = -1;
   map_xcell2face_right2_check[65] = -1;
   map_xcell2face_right2_check[66] = -1;
   map_xcell2face_right2_check[67] = -1;
   map_xcell2face_right2_check[68] = -1;
   map_xcell2face_right2_check[69] = -1;
   map_xcell2face_right2_check[70] = -1;
   map_xcell2face_right2_check[71] = -1;
   map_xcell2face_right2_check[72] = -1;
   map_xcell2face_right2_check[73] = -1;
   map_xcell2face_right2_check[74] = -1;
   map_xcell2face_right2_check[75] = -1;
   map_xcell2face_right2_check[76] = -1;
   map_xcell2face_right2_check[77] = -1;
   map_xcell2face_right2_check[78] = -1;
   map_xcell2face_right2_check[79] = -1;
   map_xcell2face_right2_check[80] = -1;
   map_xcell2face_right2_check[81] = -1;
   map_xcell2face_right2_check[82] = -1;
   map_xcell2face_right2_check[83] = -1;
   map_xcell2face_right2_check[84] = -1;
   map_xcell2face_right2_check[85] = -1;
   map_xcell2face_right2_check[86] = -1;
   map_xcell2face_right2_check[87] = -1;
   map_xcell2face_right2_check[88] = -1;
   map_xcell2face_right2_check[89] = -1;
   map_xcell2face_right2_check[90] = -1;
   map_xcell2face_right2_check[91] = -1;
   map_xcell2face_right2_check[92] = -1;
   map_xcell2face_right2_check[93] = -1;
   map_xcell2face_right2_check[94] = -1;
   map_xcell2face_right2_check[95] = -1;
   map_xcell2face_right2_check[96] = -1;
   map_xcell2face_right2_check[97] = -1;
   map_xcell2face_right2_check[98] = -1;
   map_xcell2face_right2_check[99] = -1;
   map_xcell2face_right2_check[100] = -1;
   map_xcell2face_right2_check[101] = -1;
   map_xcell2face_right2_check[102] = 80;
   map_xcell2face_right2_check[103] = -1;
   map_xcell2face_right2_check[104] = 83;
   map_xcell2face_right2_check[105] = -1;
   map_xcell2face_right2_check[106] = -1;
   map_xcell2face_right2_check[107] = -1;
   map_xcell2face_right2_check[108] = -1;
   map_xcell2face_right2_check[109] = -1;
   map_xcell2face_right2_check[110] = 87;
   map_xcell2face_right2_check[111] = -1;

   return(map_xcell2face_right2_check);
}

int *set_map_ycell2face_bot1_check(int ncells)
{
   int *map_ycell2face_bot1_check = (int *)malloc(ncells*sizeof(int));

   map_ycell2face_bot1_check[0] = -1;
   map_ycell2face_bot1_check[1] = -1;
   map_ycell2face_bot1_check[2] = -1;
   map_ycell2face_bot1_check[3] = -1;
   map_ycell2face_bot1_check[4] = -1;
   map_ycell2face_bot1_check[6] = 2;
   map_ycell2face_bot1_check[7] = 5;
   map_ycell2face_bot1_check[8] = 1;
   map_ycell2face_bot1_check[9] = 22;
   map_ycell2face_bot1_check[10] = 6;
   map_ycell2face_bot1_check[11] = 8;
   map_ycell2face_bot1_check[12] = 11;
   map_ycell2face_bot1_check[13] = 7;
   map_ycell2face_bot1_check[14] = 17;
   map_ycell2face_bot1_check[15] = 12;
   map_ycell2face_bot1_check[16] = 15;
   map_ycell2face_bot1_check[17] = 16;
   map_ycell2face_bot1_check[18] = 19;
   map_ycell2face_bot1_check[19] = 18;
   map_ycell2face_bot1_check[20] = 20;
   map_ycell2face_bot1_check[21] = 21;
   map_ycell2face_bot1_check[22] = 24;
   map_ycell2face_bot1_check[23] = 23;
   map_ycell2face_bot1_check[24] = -1;
   map_ycell2face_bot1_check[25] = -1;
   map_ycell2face_bot1_check[26] = -1;
   map_ycell2face_bot1_check[27] = -1;
   map_ycell2face_bot1_check[28] = -1;
   map_ycell2face_bot1_check[29] = -1;
   map_ycell2face_bot1_check[30] = -1;
   map_ycell2face_bot1_check[31] = 25;
   map_ycell2face_bot1_check[32] = 29;
   map_ycell2face_bot1_check[33] = -1;
   map_ycell2face_bot1_check[34] = -1;
   map_ycell2face_bot1_check[35] = -1;
   map_ycell2face_bot1_check[36] = -1;
   map_ycell2face_bot1_check[37] = -1;
   map_ycell2face_bot1_check[38] = -1;
   map_ycell2face_bot1_check[39] = 33;
   map_ycell2face_bot1_check[40] = 31;
   map_ycell2face_bot1_check[41] = 30;
   map_ycell2face_bot1_check[42] = 34;
   map_ycell2face_bot1_check[43] = 39;
   map_ycell2face_bot1_check[44] = 38;
   map_ycell2face_bot1_check[45] = 40;
   map_ycell2face_bot1_check[46] = 41;
   map_ycell2face_bot1_check[47] = 28;
   map_ycell2face_bot1_check[48] = 43;
   map_ycell2face_bot1_check[49] = 27;
   map_ycell2face_bot1_check[50] = 26;
   map_ycell2face_bot1_check[51] = 44;
   map_ycell2face_bot1_check[52] = 45;
   map_ycell2face_bot1_check[53] = 42;
   map_ycell2face_bot1_check[54] = 47;
   map_ycell2face_bot1_check[55] = 46;
   map_ycell2face_bot1_check[56] = 49;
   map_ycell2face_bot1_check[57] = 48;
   map_ycell2face_bot1_check[58] = 51;
   map_ycell2face_bot1_check[59] = 50;
   map_ycell2face_bot1_check[60] = 53;
   map_ycell2face_bot1_check[61] = 54;
   map_ycell2face_bot1_check[62] = 57;
   map_ycell2face_bot1_check[63] = 52;
   map_ycell2face_bot1_check[64] = 60;
   map_ycell2face_bot1_check[65] = 62;
   map_ycell2face_bot1_check[66] = 61;
   map_ycell2face_bot1_check[67] = 37;
   map_ycell2face_bot1_check[68] = 36;
   map_ycell2face_bot1_check[69] = 35;
   map_ycell2face_bot1_check[70] = 63;
   map_ycell2face_bot1_check[71] = 66;
   map_ycell2face_bot1_check[72] = 32;
   map_ycell2face_bot1_check[73] = -1;
   map_ycell2face_bot1_check[74] = -1;
   map_ycell2face_bot1_check[75] = 64;
   map_ycell2face_bot1_check[76] = -1;
   map_ycell2face_bot1_check[77] = -1;
   map_ycell2face_bot1_check[78] = 67;
   map_ycell2face_bot1_check[79] = 58;
   map_ycell2face_bot1_check[80] = 55;
   map_ycell2face_bot1_check[81] = 68;
   map_ycell2face_bot1_check[82] = -1;
   map_ycell2face_bot1_check[83] = -1;
   map_ycell2face_bot1_check[84] = -1;
   map_ycell2face_bot1_check[85] = -1;
   map_ycell2face_bot1_check[86] = 70;
   map_ycell2face_bot1_check[87] = 69;
   map_ycell2face_bot1_check[88] = 83;
   map_ycell2face_bot1_check[89] = 72;
   map_ycell2face_bot1_check[90] = 74;
   map_ycell2face_bot1_check[91] = 73;
   map_ycell2face_bot1_check[92] = 78;
   map_ycell2face_bot1_check[93] = 75;
   map_ycell2face_bot1_check[94] = 76;
   map_ycell2face_bot1_check[95] = 14;
   map_ycell2face_bot1_check[96] = 13;
   map_ycell2face_bot1_check[97] = 77;
   map_ycell2face_bot1_check[98] = 80;
   map_ycell2face_bot1_check[99] = 10;
   map_ycell2face_bot1_check[100] = 9;
   map_ycell2face_bot1_check[101] = 81;
   map_ycell2face_bot1_check[102] = 82;
   map_ycell2face_bot1_check[103] = 85;
   map_ycell2face_bot1_check[104] = 4;
   map_ycell2face_bot1_check[105] = 3;
   map_ycell2face_bot1_check[106] = 86;
   map_ycell2face_bot1_check[107] = -1;
   map_ycell2face_bot1_check[108] = -1;
   map_ycell2face_bot1_check[109] = -1;
   map_ycell2face_bot1_check[110] = 87;
   map_ycell2face_bot1_check[111] = -1;

   return(map_ycell2face_bot1_check);
}

int *set_map_ycell2face_bot2_check(int ncells)
{
   int *map_ycell2face_bot2_check = (int *)malloc(ncells*sizeof(int));

   map_ycell2face_bot2_check[0] = -1;
   map_ycell2face_bot2_check[1] = -1;
   map_ycell2face_bot2_check[2] = -1;
   map_ycell2face_bot2_check[3] = -1;
   map_ycell2face_bot2_check[4] = -1;
   map_ycell2face_bot2_check[5] = -1;
   map_ycell2face_bot2_check[6] = -1;
   map_ycell2face_bot2_check[7] = -1;
   map_ycell2face_bot2_check[8] = -1;
   map_ycell2face_bot2_check[9] = -1;
   map_ycell2face_bot2_check[10] = -1;
   map_ycell2face_bot2_check[11] = -1;
   map_ycell2face_bot2_check[12] = -1;
   map_ycell2face_bot2_check[13] = -1;
   map_ycell2face_bot2_check[14] = -1;
   map_ycell2face_bot2_check[15] = -1;
   map_ycell2face_bot2_check[16] = -1;
   map_ycell2face_bot2_check[17] = -1;
   map_ycell2face_bot2_check[18] = -1;
   map_ycell2face_bot2_check[19] = -1;
   map_ycell2face_bot2_check[20] = -1;
   map_ycell2face_bot2_check[21] = -1;
   map_ycell2face_bot2_check[22] = -1;
   map_ycell2face_bot2_check[23] = -1;
   map_ycell2face_bot2_check[24] = -1;
   map_ycell2face_bot2_check[25] = -1;
   map_ycell2face_bot2_check[26] = -1;
   map_ycell2face_bot2_check[27] = -1;
   map_ycell2face_bot2_check[28] = -1;
   map_ycell2face_bot2_check[29] = -1;
   map_ycell2face_bot2_check[30] = -1;
   map_ycell2face_bot2_check[31] = -1;
   map_ycell2face_bot2_check[32] = -1;
   map_ycell2face_bot2_check[33] = -1;
   map_ycell2face_bot2_check[34] = -1;
   map_ycell2face_bot2_check[35] = -1;
   map_ycell2face_bot2_check[36] = -1;
   map_ycell2face_bot2_check[37] = -1;
   map_ycell2face_bot2_check[38] = -1;
   map_ycell2face_bot2_check[39] = -1;
   map_ycell2face_bot2_check[40] = -1;
   map_ycell2face_bot2_check[41] = -1;
   map_ycell2face_bot2_check[42] = -1;
   map_ycell2face_bot2_check[43] = -1;
   map_ycell2face_bot2_check[44] = -1;
   map_ycell2face_bot2_check[45] = -1;
   map_ycell2face_bot2_check[46] = -1;
   map_ycell2face_bot2_check[47] = -1;
   map_ycell2face_bot2_check[48] = -1;
   map_ycell2face_bot2_check[49] = -1;
   map_ycell2face_bot2_check[50] = -1;
   map_ycell2face_bot2_check[51] = -1;
   map_ycell2face_bot2_check[52] = -1;
   map_ycell2face_bot2_check[53] = -1;
   map_ycell2face_bot2_check[54] = -1;
   map_ycell2face_bot2_check[55] = -1;
   map_ycell2face_bot2_check[56] = -1;
   map_ycell2face_bot2_check[57] = -1;
   map_ycell2face_bot2_check[58] = -1;
   map_ycell2face_bot2_check[59] = -1;
   map_ycell2face_bot2_check[60] = -1;
   map_ycell2face_bot2_check[61] = -1;
   map_ycell2face_bot2_check[62] = -1;
   map_ycell2face_bot2_check[63] = -1;
   map_ycell2face_bot2_check[64] = 59;
   map_ycell2face_bot2_check[65] = -1;
   map_ycell2face_bot2_check[66] = -1;
   map_ycell2face_bot2_check[67] = -1;
   map_ycell2face_bot2_check[68] = -1;
   map_ycell2face_bot2_check[69] = -1;
   map_ycell2face_bot2_check[70] = -1;
   map_ycell2face_bot2_check[71] = -1;
   map_ycell2face_bot2_check[72] = -1;
   map_ycell2face_bot2_check[73] = -1;
   map_ycell2face_bot2_check[74] = -1;
   map_ycell2face_bot2_check[75] = 65;
   map_ycell2face_bot2_check[76] = -1;
   map_ycell2face_bot2_check[77] = -1;
   map_ycell2face_bot2_check[78] = -1;
   map_ycell2face_bot2_check[79] = -1;
   map_ycell2face_bot2_check[80] = 56;
   map_ycell2face_bot2_check[81] = -1;
   map_ycell2face_bot2_check[82] = -1;
   map_ycell2face_bot2_check[83] = -1;
   map_ycell2face_bot2_check[84] = -1;
   map_ycell2face_bot2_check[85] = -1;
   map_ycell2face_bot2_check[86] = -1;
   map_ycell2face_bot2_check[87] = -1;
   map_ycell2face_bot2_check[88] = -1;
   map_ycell2face_bot2_check[89] = 71;
   map_ycell2face_bot2_check[90] = -1;
   map_ycell2face_bot2_check[91] = -1;
   map_ycell2face_bot2_check[92] = -1;
   map_ycell2face_bot2_check[93] = -1;
   map_ycell2face_bot2_check[94] = -1;
   map_ycell2face_bot2_check[95] = -1;
   map_ycell2face_bot2_check[96] = -1;
   map_ycell2face_bot2_check[97] = -1;
   map_ycell2face_bot2_check[98] = -1;
   map_ycell2face_bot2_check[99] = -1;
   map_ycell2face_bot2_check[100] = -1;
   map_ycell2face_bot2_check[101] = -1;
   map_ycell2face_bot2_check[102] = 79;
   map_ycell2face_bot2_check[103] = -1;
   map_ycell2face_bot2_check[104] = -1;
   map_ycell2face_bot2_check[105] = -1;
   map_ycell2face_bot2_check[106] = -1;
   map_ycell2face_bot2_check[107] = -1;
   map_ycell2face_bot2_check[108] = -1;
   map_ycell2face_bot2_check[109] = -1;
   map_ycell2face_bot2_check[110] = 84;
   map_ycell2face_bot2_check[111] = -1;

   return(map_ycell2face_bot2_check);
}

int *set_map_ycell2face_top1_check(int ncells)
{
   int *map_ycell2face_top1_check = (int *)malloc(ncells*sizeof(int));

   map_ycell2face_top1_check[0] = -1;
   map_ycell2face_top1_check[2] = -1;
   map_ycell2face_top1_check[3] = -1;
   map_ycell2face_top1_check[4] = -1;
   map_ycell2face_top1_check[5] = 2;
   map_ycell2face_top1_check[6] = 3;
   map_ycell2face_top1_check[7] = 4;
   map_ycell2face_top1_check[8] = 5;
   map_ycell2face_top1_check[9] = 6;
   map_ycell2face_top1_check[10] = 8;
   map_ycell2face_top1_check[11] = 9;
   map_ycell2face_top1_check[12] = 10;
   map_ycell2face_top1_check[13] = 11;
   map_ycell2face_top1_check[14] = 12;
   map_ycell2face_top1_check[15] = 13;
   map_ycell2face_top1_check[16] = 14;
   map_ycell2face_top1_check[17] = 15;
   map_ycell2face_top1_check[18] = 16;
   map_ycell2face_top1_check[19] = 17;
   map_ycell2face_top1_check[20] = 18;
   map_ycell2face_top1_check[21] = 19;
   map_ycell2face_top1_check[22] = 20;
   map_ycell2face_top1_check[23] = 22;
   map_ycell2face_top1_check[24] = 23;
   map_ycell2face_top1_check[25] = 24;
   map_ycell2face_top1_check[26] = -1;
   map_ycell2face_top1_check[27] = -1;
   map_ycell2face_top1_check[28] = -1;
   map_ycell2face_top1_check[29] = -1;
   map_ycell2face_top1_check[30] = 25;
   map_ycell2face_top1_check[31] = 26;
   map_ycell2face_top1_check[32] = 28;
   map_ycell2face_top1_check[33] = 29;
   map_ycell2face_top1_check[34] = -1;
   map_ycell2face_top1_check[35] = -1;
   map_ycell2face_top1_check[36] = 30;
   map_ycell2face_top1_check[37] = -1;
   map_ycell2face_top1_check[38] = -1;
   map_ycell2face_top1_check[39] = 32;
   map_ycell2face_top1_check[40] = 33;
   map_ycell2face_top1_check[41] = 34;
   map_ycell2face_top1_check[42] = 35;
   map_ycell2face_top1_check[43] = 36;
   map_ycell2face_top1_check[44] = 37;
   map_ycell2face_top1_check[45] = 38;
   map_ycell2face_top1_check[46] = 39;
   map_ycell2face_top1_check[47] = 40;
   map_ycell2face_top1_check[48] = 42;
   map_ycell2face_top1_check[49] = 43;
   map_ycell2face_top1_check[50] = 44;
   map_ycell2face_top1_check[51] = 45;
   map_ycell2face_top1_check[52] = 46;
   map_ycell2face_top1_check[53] = 47;
   map_ycell2face_top1_check[54] = 48;
   map_ycell2face_top1_check[55] = 49;
   map_ycell2face_top1_check[56] = 50;
   map_ycell2face_top1_check[57] = 51;
   map_ycell2face_top1_check[58] = 52;
   map_ycell2face_top1_check[59] = 53;
   map_ycell2face_top1_check[60] = 54;
   map_ycell2face_top1_check[61] = 55;
   map_ycell2face_top1_check[62] = 56;
   map_ycell2face_top1_check[63] = 57;
   map_ycell2face_top1_check[64] = 58;
   map_ycell2face_top1_check[65] = 59;
   map_ycell2face_top1_check[66] = 60;
   map_ycell2face_top1_check[67] = 61;
   map_ycell2face_top1_check[68] = 62;
   map_ycell2face_top1_check[69] = 63;
   map_ycell2face_top1_check[70] = 64;
   map_ycell2face_top1_check[71] = 65;
   map_ycell2face_top1_check[72] = 66;
   map_ycell2face_top1_check[73] = -1;
   map_ycell2face_top1_check[74] = -1;
   map_ycell2face_top1_check[75] = -1;
   map_ycell2face_top1_check[76] = -1;
   map_ycell2face_top1_check[77] = -1;
   map_ycell2face_top1_check[78] = -1;
   map_ycell2face_top1_check[79] = 67;
   map_ycell2face_top1_check[80] = 68;
   map_ycell2face_top1_check[81] = -1;
   map_ycell2face_top1_check[82] = -1;
   map_ycell2face_top1_check[83] = -1;
   map_ycell2face_top1_check[84] = -1;
   map_ycell2face_top1_check[85] = -1;
   map_ycell2face_top1_check[86] = -1;
   map_ycell2face_top1_check[87] = -1;
   map_ycell2face_top1_check[88] = 69;
   map_ycell2face_top1_check[89] = 70;
   map_ycell2face_top1_check[90] = 71;
   map_ycell2face_top1_check[91] = 72;
   map_ycell2face_top1_check[92] = 73;
   map_ycell2face_top1_check[93] = 74;
   map_ycell2face_top1_check[94] = 75;
   map_ycell2face_top1_check[95] = 76;
   map_ycell2face_top1_check[96] = 77;
   map_ycell2face_top1_check[97] = 78;
   map_ycell2face_top1_check[98] = 79;
   map_ycell2face_top1_check[99] = 80;
   map_ycell2face_top1_check[100] = 81;
   map_ycell2face_top1_check[101] = 82;
   map_ycell2face_top1_check[102] = 83;
   map_ycell2face_top1_check[103] = 84;
   map_ycell2face_top1_check[104] = 85;
   map_ycell2face_top1_check[105] = 86;
   map_ycell2face_top1_check[106] = 87;
   map_ycell2face_top1_check[107] = -1;
   map_ycell2face_top1_check[108] = -1;
   map_ycell2face_top1_check[109] = -1;
   map_ycell2face_top1_check[110] = -1;
   map_ycell2face_top1_check[111] = -1;

   return(map_ycell2face_top1_check);
}

int *set_map_ycell2face_top2_check(int ncells)
{
   int *map_ycell2face_top2_check = (int *)malloc(ncells*sizeof(int));

   map_ycell2face_top2_check[0] = -1;
   map_ycell2face_top2_check[1] = 1;
   map_ycell2face_top2_check[2] = -1;
   map_ycell2face_top2_check[3] = -1;
   map_ycell2face_top2_check[4] = -1;
   map_ycell2face_top2_check[5] = -1;
   map_ycell2face_top2_check[6] = -1;
   map_ycell2face_top2_check[7] = -1;
   map_ycell2face_top2_check[8] = -1;
   map_ycell2face_top2_check[9] = 7;
   map_ycell2face_top2_check[10] = -1;
   map_ycell2face_top2_check[11] = -1;
   map_ycell2face_top2_check[12] = -1;
   map_ycell2face_top2_check[13] = -1;
   map_ycell2face_top2_check[14] = -1;
   map_ycell2face_top2_check[15] = -1;
   map_ycell2face_top2_check[16] = -1;
   map_ycell2face_top2_check[17] = -1;
   map_ycell2face_top2_check[18] = -1;
   map_ycell2face_top2_check[19] = -1;
   map_ycell2face_top2_check[20] = -1;
   map_ycell2face_top2_check[21] = -1;
   map_ycell2face_top2_check[22] = 21;
   map_ycell2face_top2_check[23] = -1;
   map_ycell2face_top2_check[24] = -1;
   map_ycell2face_top2_check[25] = -1;
   map_ycell2face_top2_check[26] = -1;
   map_ycell2face_top2_check[27] = -1;
   map_ycell2face_top2_check[28] = -1;
   map_ycell2face_top2_check[29] = -1;
   map_ycell2face_top2_check[30] = -1;
   map_ycell2face_top2_check[31] = 27;
   map_ycell2face_top2_check[32] = -1;
   map_ycell2face_top2_check[33] = -1;
   map_ycell2face_top2_check[34] = -1;
   map_ycell2face_top2_check[35] = -1;
   map_ycell2face_top2_check[36] = 31;
   map_ycell2face_top2_check[37] = -1;
   map_ycell2face_top2_check[38] = -1;
   map_ycell2face_top2_check[39] = -1;
   map_ycell2face_top2_check[40] = -1;
   map_ycell2face_top2_check[41] = -1;
   map_ycell2face_top2_check[42] = -1;
   map_ycell2face_top2_check[43] = -1;
   map_ycell2face_top2_check[44] = -1;
   map_ycell2face_top2_check[45] = -1;
   map_ycell2face_top2_check[46] = -1;
   map_ycell2face_top2_check[47] = 41;
   map_ycell2face_top2_check[48] = -1;
   map_ycell2face_top2_check[49] = -1;
   map_ycell2face_top2_check[50] = -1;
   map_ycell2face_top2_check[51] = -1;
   map_ycell2face_top2_check[52] = -1;
   map_ycell2face_top2_check[53] = -1;
   map_ycell2face_top2_check[54] = -1;
   map_ycell2face_top2_check[55] = -1;
   map_ycell2face_top2_check[56] = -1;
   map_ycell2face_top2_check[57] = -1;
   map_ycell2face_top2_check[58] = -1;
   map_ycell2face_top2_check[59] = -1;
   map_ycell2face_top2_check[60] = -1;
   map_ycell2face_top2_check[61] = -1;
   map_ycell2face_top2_check[62] = -1;
   map_ycell2face_top2_check[63] = -1;
   map_ycell2face_top2_check[64] = -1;
   map_ycell2face_top2_check[65] = -1;
   map_ycell2face_top2_check[66] = -1;
   map_ycell2face_top2_check[67] = -1;
   map_ycell2face_top2_check[68] = -1;
   map_ycell2face_top2_check[69] = -1;
   map_ycell2face_top2_check[70] = -1;
   map_ycell2face_top2_check[71] = -1;
   map_ycell2face_top2_check[72] = -1;
   map_ycell2face_top2_check[73] = -1;
   map_ycell2face_top2_check[74] = -1;
   map_ycell2face_top2_check[75] = -1;
   map_ycell2face_top2_check[76] = -1;
   map_ycell2face_top2_check[77] = -1;
   map_ycell2face_top2_check[78] = -1;
   map_ycell2face_top2_check[79] = -1;
   map_ycell2face_top2_check[80] = -1;
   map_ycell2face_top2_check[81] = -1;
   map_ycell2face_top2_check[82] = -1;
   map_ycell2face_top2_check[83] = -1;
   map_ycell2face_top2_check[84] = -1;
   map_ycell2face_top2_check[85] = -1;
   map_ycell2face_top2_check[86] = -1;
   map_ycell2face_top2_check[87] = -1;
   map_ycell2face_top2_check[88] = -1;
   map_ycell2face_top2_check[89] = -1;
   map_ycell2face_top2_check[90] = -1;
   map_ycell2face_top2_check[91] = -1;
   map_ycell2face_top2_check[92] = -1;
   map_ycell2face_top2_check[93] = -1;
   map_ycell2face_top2_check[94] = -1;
   map_ycell2face_top2_check[95] = -1;
   map_ycell2face_top2_check[96] = -1;
   map_ycell2face_top2_check[97] = -1;
   map_ycell2face_top2_check[98] = -1;
   map_ycell2face_top2_check[99] = -1;
   map_ycell2face_top2_check[100] = -1;
   map_ycell2face_top2_check[101] = -1;
   map_ycell2face_top2_check[102] = -1;
   map_ycell2face_top2_check[103] = -1;
   map_ycell2face_top2_check[104] = -1;
   map_ycell2face_top2_check[105] = -1;
   map_ycell2face_top2_check[106] = -1;
   map_ycell2face_top2_check[107] = -1;
   map_ycell2face_top2_check[108] = -1;
   map_ycell2face_top2_check[109] = -1;
   map_ycell2face_top2_check[110] = -1;
   map_ycell2face_top2_check[111] = -1;

   return(map_ycell2face_top2_check);
}
