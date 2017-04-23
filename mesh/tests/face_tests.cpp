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

   int nxface_check = 110;
   vector<int>xface_i_check(nxface_check);
   vector<int>xface_j_check(nxface_check);
   vector<int>xface_level_check(nxface_check);

   set_xface_check(&xface_i_check[0], &xface_j_check[0], &xface_level_check[0]);

   int nyface_check = 110;
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
   }

   if (mesh->nyface != nyface_check){
      printf("Error -- nyface does not match, nyface %d nyface_check %d\n",
         mesh->nyface, nyface_check);
      icount_err++;
   }

   for (int iface = 0; iface < mesh->nxface; iface++){
      //printf("   xface_i_check[%d] = %d;\n",iface, mesh->xface_i[iface]);
      //printf("   xface_j_check[%d] = %d;\n",iface, mesh->xface_j[iface]);
      //printf("   xface_level_check[%d] = %d;\n",iface, mesh->xface_level[iface]);
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
      //printf("   yface_i_check[%d] = %d;\n",iface, mesh->yface_i[iface]);
      //printf("   yface_j_check[%d] = %d;\n",iface, mesh->yface_j[iface]);
      //printf("   yface_level_check[%d] = %d;\n",iface, mesh->yface_level[iface]);
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
      //printf("   map_xface2cell_lower_check[%d] = %d;\n",iface, mesh->map_xface2cell_lower[iface]);
      //printf("   map_xface2cell_upper_check[%d] = %d;\n",iface, mesh->map_xface2cell_upper[iface]);
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
      //printf("   map_yface2cell_lower_check[%d] = %d;\n",iface, mesh->map_yface2cell_lower[iface]);
      //printf("   map_yface2cell_upper_check[%d] = %d;\n",iface, mesh->map_yface2cell_upper[iface]);
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
      //printf("   map_xcell2face_left1_check[%d] = %d;\n", iz,mesh->map_xcell2face_left1[iz]);
      if (mesh->map_xcell2face_left1[iz] != map_xcell2face_left1_check[iz]) {
         printf("Error -- map_xcell2face_left1 does not match for cell %d map_xcell2face_left1 %d check %d\n",
            iz,mesh->map_xcell2face_left1[iz], map_xcell2face_left1_check[iz]);
      }
      //printf("   map_xcell2face_left2_check[%d] = %d;\n", iz,mesh->map_xcell2face_left2[iz]);
      if (mesh->map_xcell2face_left2[iz] != map_xcell2face_left2_check[iz]) {
         printf("Error -- map_xcell2face_left2 does not match for cell %d map_xcell2face_left2 %d check %d\n",
            iz,mesh->map_xcell2face_left2[iz], map_xcell2face_left2_check[iz]);
      }
      //printf("   map_xcell2face_right1_check[%d] = %d;\n", iz,mesh->map_xcell2face_right1[iz]);
      if (mesh->map_xcell2face_right1[iz] != map_xcell2face_right1_check[iz]) {
         printf("Error -- map_xcell2face_right1 does not match for cell %d map_xcell2face_right1 %d check %d\n",
            iz,mesh->map_xcell2face_right1[iz], map_xcell2face_right1_check[iz]);
      }
      //printf("   map_xcell2face_right2_check[%d] = %d;\n", iz,mesh->map_xcell2face_right2[iz]);
      if (mesh->map_xcell2face_right2[iz] != map_xcell2face_right2_check[iz]) {
         printf("Error -- map_xcell2face_right2 does not match for cell %d map_xcell2face_right2 %d check %d\n",
            iz,mesh->map_xcell2face_right2[iz], map_xcell2face_right2_check[iz]);
      }

      //printf("   map_ycell2face_bot1_check[%d] = %d;\n", iz,mesh->map_ycell2face_bot1[iz]);
      if (mesh->map_ycell2face_bot1[iz] != map_ycell2face_bot1_check[iz]) {
         printf("Error -- map_ycell2face_bot1 does not match for cell %d map_ycell2face_bot1 %d check %d\n",
            iz,mesh->map_ycell2face_bot1[iz], map_ycell2face_bot1_check[iz]);
      }
      //printf("   map_ycell2face_bot2_check[%d] = %d;\n", iz,mesh->map_ycell2face_bot2[iz]);
      if (mesh->map_ycell2face_bot2[iz] != map_ycell2face_bot2_check[iz]) {
         printf("Error -- map_ycell2face_bot2 does not match for cell %d map_ycell2face_bot2 %d check %d\n",
            iz,mesh->map_ycell2face_bot2[iz], map_ycell2face_bot2_check[iz]);
      }
      //printf("   map_ycell2face_top1_check[%d] = %d;\n", iz,mesh->map_ycell2face_top1[iz]);
      if (mesh->map_ycell2face_top1[iz] != map_ycell2face_top1_check[iz]) {
         printf("Error -- map_ycell2face_top1 does not match for cell %d map_ycell2face_top1 %d check %d\n",
            iz,mesh->map_ycell2face_top1[iz], map_ycell2face_top1_check[iz]);
      }
      //printf("   map_ycell2face_top2_check[%d] = %d;\n", iz,mesh->map_ycell2face_top2[iz]);
      if (mesh->map_ycell2face_top2[iz] != map_ycell2face_top2_check[iz]) {
         printf("Error -- map_ycell2face_top2 does not match for cell %d map_ycell2face_top2 %d check %d\n",
            iz,mesh->map_ycell2face_top2[iz], map_ycell2face_top2_check[iz]);
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

   for (int lev = 0; lev < levmx+1; lev++){

      int isize, jsize;

      if (lev == 0) {
         isize = 0;
         jsize = 0;
      } else if (lev == 1) {
         isize = 9;
         jsize = 10;
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

   for (int lev = 0; lev < levmx+1; lev++){

      int isize, jsize;

      if (lev == 0) {
         isize = 0;
         jsize = 0;
      } else if (lev == 1) {
         isize = 10;
         jsize = 9;
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
         isize = 12;
         jsize = 12;
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
   xface_i_check[0] = 2;
   xface_i_check[1] = 2;
   xface_i_check[2] = 2;
   xface_i_check[3] = 0;
   xface_i_check[4] = 0;
   xface_i_check[5] = 0;
   xface_i_check[6] = 1;
   xface_i_check[7] = 1;
   xface_i_check[8] = 0;
   xface_i_check[9] = 0;
   xface_i_check[10] = 2;
   xface_i_check[11] = 2;
   xface_i_check[12] = 2;
   xface_i_check[13] = 1;
   xface_i_check[14] = 1;
   xface_i_check[15] = 2;
   xface_i_check[16] = 2;
   xface_i_check[17] = 3;
   xface_i_check[18] = 3;
   xface_i_check[19] = 4;
   xface_i_check[20] = 4;
   xface_i_check[21] = 4;
   xface_i_check[22] = 3;
   xface_i_check[23] = 3;
   xface_i_check[24] = 4;
   xface_i_check[25] = 4;
   xface_i_check[26] = 3;
   xface_i_check[27] = 3;
   xface_i_check[28] = 4;
   xface_i_check[29] = 3;
   xface_i_check[30] = 4;
   xface_i_check[31] = 5;
   xface_i_check[32] = 6;
   xface_i_check[33] = 5;
   xface_i_check[34] = 5;
   xface_i_check[35] = 6;
   xface_i_check[36] = 6;
   xface_i_check[37] = 4;
   xface_i_check[38] = 8;
   xface_i_check[39] = 8;
   xface_i_check[40] = 7;
   xface_i_check[41] = 7;
   xface_i_check[42] = 8;
   xface_i_check[43] = 7;
   xface_i_check[44] = 7;
   xface_i_check[45] = 8;
   xface_i_check[46] = 6;
   xface_i_check[47] = 6;
   xface_i_check[48] = 6;
   xface_i_check[49] = 5;
   xface_i_check[50] = 5;
   xface_i_check[51] = 5;
   xface_i_check[52] = 6;
   xface_i_check[53] = 6;
   xface_i_check[54] = 5;
   xface_i_check[55] = 5;
   xface_i_check[56] = 6;
   xface_i_check[57] = 6;
   xface_i_check[58] = 5;
   xface_i_check[59] = 5;
   xface_i_check[60] = 5;
   xface_i_check[61] = 6;
   xface_i_check[62] = 6;
   xface_i_check[63] = 6;
   xface_i_check[64] = 8;
   xface_i_check[65] = 7;
   xface_i_check[66] = 7;
   xface_i_check[67] = 8;
   xface_i_check[68] = 7;
   xface_i_check[69] = 7;
   xface_i_check[70] = 8;
   xface_i_check[71] = 8;
   xface_i_check[72] = 4;
   xface_i_check[73] = 6;
   xface_i_check[74] = 6;
   xface_i_check[75] = 5;
   xface_i_check[76] = 5;
   xface_i_check[77] = 6;
   xface_i_check[78] = 5;
   xface_i_check[79] = 4;
   xface_i_check[80] = 3;
   xface_i_check[81] = 4;
   xface_i_check[82] = 3;
   xface_i_check[83] = 3;
   xface_i_check[84] = 4;
   xface_i_check[85] = 4;
   xface_i_check[86] = 3;
   xface_i_check[87] = 3;
   xface_i_check[88] = 4;
   xface_i_check[89] = 4;
   xface_i_check[90] = 4;
   xface_i_check[91] = 3;
   xface_i_check[92] = 3;
   xface_i_check[93] = 2;
   xface_i_check[94] = 2;
   xface_i_check[95] = 1;
   xface_i_check[96] = 1;
   xface_i_check[97] = 2;
   xface_i_check[98] = 2;
   xface_i_check[99] = 2;
   xface_i_check[100] = 0;
   xface_i_check[101] = 0;
   xface_i_check[102] = 1;
   xface_i_check[103] = 1;
   xface_i_check[104] = 0;
   xface_i_check[105] = 0;
   xface_i_check[106] = 0;
   xface_i_check[107] = 2;
   xface_i_check[108] = 2;
   xface_i_check[109] = 2;

   xface_j_check[0] = 0;
   xface_j_check[1] = 1;
   xface_j_check[2] = 2;
   xface_j_check[3] = 0;
   xface_j_check[4] = 3;
   xface_j_check[5] = 4;
   xface_j_check[6] = 3;
   xface_j_check[7] = 4;
   xface_j_check[8] = 2;
   xface_j_check[9] = 3;
   xface_j_check[10] = 3;
   xface_j_check[11] = 0;
   xface_j_check[12] = 1;
   xface_j_check[13] = 2;
   xface_j_check[14] = 3;
   xface_j_check[15] = 3;
   xface_j_check[16] = 2;
   xface_j_check[17] = 2;
   xface_j_check[18] = 3;
   xface_j_check[19] = 3;
   xface_j_check[20] = 2;
   xface_j_check[21] = 1;
   xface_j_check[22] = 1;
   xface_j_check[23] = 0;
   xface_j_check[24] = 0;
   xface_j_check[25] = 2;
   xface_j_check[26] = 2;
   xface_j_check[27] = 1;
   xface_j_check[28] = 1;
   xface_j_check[29] = 0;
   xface_j_check[30] = 0;
   xface_j_check[31] = 0;
   xface_j_check[32] = 0;
   xface_j_check[33] = 1;
   xface_j_check[34] = 2;
   xface_j_check[35] = 2;
   xface_j_check[36] = 1;
   xface_j_check[37] = 0;
   xface_j_check[38] = 4;
   xface_j_check[39] = 3;
   xface_j_check[40] = 3;
   xface_j_check[41] = 4;
   xface_j_check[42] = 3;
   xface_j_check[43] = 3;
   xface_j_check[44] = 2;
   xface_j_check[45] = 2;
   xface_j_check[46] = 3;
   xface_j_check[47] = 1;
   xface_j_check[48] = 0;
   xface_j_check[49] = 0;
   xface_j_check[50] = 1;
   xface_j_check[51] = 2;
   xface_j_check[52] = 2;
   xface_j_check[53] = 3;
   xface_j_check[54] = 3;
   xface_j_check[55] = 4;
   xface_j_check[56] = 4;
   xface_j_check[57] = 5;
   xface_j_check[58] = 5;
   xface_j_check[59] = 6;
   xface_j_check[60] = 7;
   xface_j_check[61] = 7;
   xface_j_check[62] = 6;
   xface_j_check[63] = 6;
   xface_j_check[64] = 5;
   xface_j_check[65] = 5;
   xface_j_check[66] = 4;
   xface_j_check[67] = 4;
   xface_j_check[68] = 5;
   xface_j_check[69] = 6;
   xface_j_check[70] = 6;
   xface_j_check[71] = 5;
   xface_j_check[72] = 3;
   xface_j_check[73] = 8;
   xface_j_check[74] = 7;
   xface_j_check[75] = 7;
   xface_j_check[76] = 8;
   xface_j_check[77] = 9;
   xface_j_check[78] = 9;
   xface_j_check[79] = 9;
   xface_j_check[80] = 9;
   xface_j_check[81] = 8;
   xface_j_check[82] = 8;
   xface_j_check[83] = 7;
   xface_j_check[84] = 7;
   xface_j_check[85] = 7;
   xface_j_check[86] = 7;
   xface_j_check[87] = 6;
   xface_j_check[88] = 6;
   xface_j_check[89] = 5;
   xface_j_check[90] = 4;
   xface_j_check[91] = 4;
   xface_j_check[92] = 5;
   xface_j_check[93] = 5;
   xface_j_check[94] = 4;
   xface_j_check[95] = 4;
   xface_j_check[96] = 5;
   xface_j_check[97] = 6;
   xface_j_check[98] = 7;
   xface_j_check[99] = 6;
   xface_j_check[100] = 4;
   xface_j_check[101] = 5;
   xface_j_check[102] = 5;
   xface_j_check[103] = 6;
   xface_j_check[104] = 5;
   xface_j_check[105] = 6;
   xface_j_check[106] = 3;
   xface_j_check[107] = 7;
   xface_j_check[108] = 8;
   xface_j_check[109] = 9;

   xface_level_check[0] = 1;
   xface_level_check[1] = 1;
   xface_level_check[2] = 1;
   xface_level_check[3] = 0;
   xface_level_check[4] = 1;
   xface_level_check[5] = 1;
   xface_level_check[6] = 1;
   xface_level_check[7] = 1;
   xface_level_check[8] = 2;
   xface_level_check[9] = 2;
   xface_level_check[10] = 1;
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
   xface_level_check[21] = 2;
   xface_level_check[22] = 2;
   xface_level_check[23] = 2;
   xface_level_check[24] = 2;
   xface_level_check[25] = 1;
   xface_level_check[26] = 1;
   xface_level_check[27] = 1;
   xface_level_check[28] = 1;
   xface_level_check[29] = 1;
   xface_level_check[30] = 1;
   xface_level_check[31] = 1;
   xface_level_check[32] = 1;
   xface_level_check[33] = 1;
   xface_level_check[34] = 1;
   xface_level_check[35] = 1;
   xface_level_check[36] = 1;
   xface_level_check[37] = 0;
   xface_level_check[38] = 1;
   xface_level_check[39] = 1;
   xface_level_check[40] = 1;
   xface_level_check[41] = 1;
   xface_level_check[42] = 2;
   xface_level_check[43] = 2;
   xface_level_check[44] = 2;
   xface_level_check[45] = 2;
   xface_level_check[46] = 1;
   xface_level_check[47] = 2;
   xface_level_check[48] = 2;
   xface_level_check[49] = 2;
   xface_level_check[50] = 2;
   xface_level_check[51] = 2;
   xface_level_check[52] = 2;
   xface_level_check[53] = 2;
   xface_level_check[54] = 2;
   xface_level_check[55] = 2;
   xface_level_check[56] = 2;
   xface_level_check[57] = 2;
   xface_level_check[58] = 2;
   xface_level_check[59] = 2;
   xface_level_check[60] = 2;
   xface_level_check[61] = 2;
   xface_level_check[62] = 2;
   xface_level_check[63] = 1;
   xface_level_check[64] = 2;
   xface_level_check[65] = 2;
   xface_level_check[66] = 2;
   xface_level_check[67] = 2;
   xface_level_check[68] = 1;
   xface_level_check[69] = 1;
   xface_level_check[70] = 1;
   xface_level_check[71] = 1;
   xface_level_check[72] = 0;
   xface_level_check[73] = 1;
   xface_level_check[74] = 1;
   xface_level_check[75] = 1;
   xface_level_check[76] = 1;
   xface_level_check[77] = 1;
   xface_level_check[78] = 1;
   xface_level_check[79] = 1;
   xface_level_check[80] = 1;
   xface_level_check[81] = 1;
   xface_level_check[82] = 1;
   xface_level_check[83] = 1;
   xface_level_check[84] = 1;
   xface_level_check[85] = 2;
   xface_level_check[86] = 2;
   xface_level_check[87] = 2;
   xface_level_check[88] = 2;
   xface_level_check[89] = 2;
   xface_level_check[90] = 2;
   xface_level_check[91] = 2;
   xface_level_check[92] = 2;
   xface_level_check[93] = 2;
   xface_level_check[94] = 2;
   xface_level_check[95] = 2;
   xface_level_check[96] = 2;
   xface_level_check[97] = 2;
   xface_level_check[98] = 2;
   xface_level_check[99] = 1;
   xface_level_check[100] = 2;
   xface_level_check[101] = 2;
   xface_level_check[102] = 1;
   xface_level_check[103] = 1;
   xface_level_check[104] = 1;
   xface_level_check[105] = 1;
   xface_level_check[106] = 0;
   xface_level_check[107] = 1;
   xface_level_check[108] = 1;
   xface_level_check[109] = 1;
}

void set_yface_check(int *yface_i_check, int *yface_j_check, int *yface_level_check)
{
   yface_i_check[0] = 0;
   yface_i_check[1] = 1;
   yface_i_check[2] = 2;
   yface_i_check[3] = 0;
   yface_i_check[4] = 0;
   yface_i_check[5] = 0;
   yface_i_check[6] = 1;
   yface_i_check[7] = 1;
   yface_i_check[8] = 2;
   yface_i_check[9] = 2;
   yface_i_check[10] = 0;
   yface_i_check[11] = 1;
   yface_i_check[12] = 0;
   yface_i_check[13] = 0;
   yface_i_check[14] = 1;
   yface_i_check[15] = 1;
   yface_i_check[16] = 2;
   yface_i_check[17] = 2;
   yface_i_check[18] = 3;
   yface_i_check[19] = 3;
   yface_i_check[20] = 3;
   yface_i_check[21] = 2;
   yface_i_check[22] = 2;
   yface_i_check[23] = 3;
   yface_i_check[24] = 2;
   yface_i_check[25] = 3;
   yface_i_check[26] = 3;
   yface_i_check[27] = 3;
   yface_i_check[28] = 4;
   yface_i_check[29] = 3;
   yface_i_check[30] = 4;
   yface_i_check[31] = 5;
   yface_i_check[32] = 6;
   yface_i_check[33] = 5;
   yface_i_check[34] = 4;
   yface_i_check[35] = 5;
   yface_i_check[36] = 6;
   yface_i_check[37] = 6;
   yface_i_check[38] = 3;
   yface_i_check[39] = 9;
   yface_i_check[40] = 7;
   yface_i_check[41] = 8;
   yface_i_check[42] = 9;
   yface_i_check[43] = 9;
   yface_i_check[44] = 8;
   yface_i_check[45] = 8;
   yface_i_check[46] = 7;
   yface_i_check[47] = 7;
   yface_i_check[48] = 7;
   yface_i_check[49] = 6;
   yface_i_check[50] = 6;
   yface_i_check[51] = 7;
   yface_i_check[52] = 6;
   yface_i_check[53] = 7;
   yface_i_check[54] = 5;
   yface_i_check[55] = 5;
   yface_i_check[56] = 4;
   yface_i_check[57] = 4;
   yface_i_check[58] = 4;
   yface_i_check[59] = 5;
   yface_i_check[60] = 5;
   yface_i_check[61] = 4;
   yface_i_check[62] = 4;
   yface_i_check[63] = 5;
   yface_i_check[64] = 5;
   yface_i_check[65] = 4;
   yface_i_check[66] = 4;
   yface_i_check[67] = 4;
   yface_i_check[68] = 5;
   yface_i_check[69] = 5;
   yface_i_check[70] = 6;
   yface_i_check[71] = 7;
   yface_i_check[72] = 6;
   yface_i_check[73] = 6;
   yface_i_check[74] = 7;
   yface_i_check[75] = 7;
   yface_i_check[76] = 7;
   yface_i_check[77] = 8;
   yface_i_check[78] = 8;
   yface_i_check[79] = 9;
   yface_i_check[80] = 9;
   yface_i_check[81] = 3;
   yface_i_check[82] = 6;
   yface_i_check[83] = 6;
   yface_i_check[84] = 5;
   yface_i_check[85] = 5;
   yface_i_check[86] = 4;
   yface_i_check[87] = 3;
   yface_i_check[88] = 3;
   yface_i_check[89] = 4;
   yface_i_check[90] = 3;
   yface_i_check[91] = 2;
   yface_i_check[92] = 2;
   yface_i_check[93] = 3;
   yface_i_check[94] = 3;
   yface_i_check[95] = 3;
   yface_i_check[96] = 2;
   yface_i_check[97] = 2;
   yface_i_check[98] = 1;
   yface_i_check[99] = 1;
   yface_i_check[100] = 0;
   yface_i_check[101] = 0;
   yface_i_check[102] = 3;
   yface_i_check[103] = 2;
   yface_i_check[104] = 2;
   yface_i_check[105] = 1;
   yface_i_check[106] = 1;
   yface_i_check[107] = 0;
   yface_i_check[108] = 0;
   yface_i_check[109] = 0;

   yface_j_check[0] = 0;
   yface_j_check[1] = 2;
   yface_j_check[2] = 2;
   yface_j_check[3] = 2;
   yface_j_check[4] = 3;
   yface_j_check[5] = 4;
   yface_j_check[6] = 3;
   yface_j_check[7] = 4;
   yface_j_check[8] = 4;
   yface_j_check[9] = 3;
   yface_j_check[10] = 2;
   yface_j_check[11] = 2;
   yface_j_check[12] = 3;
   yface_j_check[13] = 4;
   yface_j_check[14] = 4;
   yface_j_check[15] = 3;
   yface_j_check[16] = 3;
   yface_j_check[17] = 4;
   yface_j_check[18] = 4;
   yface_j_check[19] = 3;
   yface_j_check[20] = 2;
   yface_j_check[21] = 2;
   yface_j_check[22] = 1;
   yface_j_check[23] = 1;
   yface_j_check[24] = 0;
   yface_j_check[25] = 0;
   yface_j_check[26] = 2;
   yface_j_check[27] = 1;
   yface_j_check[28] = 1;
   yface_j_check[29] = 0;
   yface_j_check[30] = 0;
   yface_j_check[31] = 0;
   yface_j_check[32] = 0;
   yface_j_check[33] = 1;
   yface_j_check[34] = 0;
   yface_j_check[35] = 0;
   yface_j_check[36] = 2;
   yface_j_check[37] = 1;
   yface_j_check[38] = 0;
   yface_j_check[39] = 2;
   yface_j_check[40] = 2;
   yface_j_check[41] = 2;
   yface_j_check[42] = 3;
   yface_j_check[43] = 4;
   yface_j_check[44] = 4;
   yface_j_check[45] = 3;
   yface_j_check[46] = 3;
   yface_j_check[47] = 4;
   yface_j_check[48] = 4;
   yface_j_check[49] = 4;
   yface_j_check[50] = 3;
   yface_j_check[51] = 3;
   yface_j_check[52] = 2;
   yface_j_check[53] = 2;
   yface_j_check[54] = 2;
   yface_j_check[55] = 1;
   yface_j_check[56] = 1;
   yface_j_check[57] = 2;
   yface_j_check[58] = 3;
   yface_j_check[59] = 3;
   yface_j_check[60] = 4;
   yface_j_check[61] = 4;
   yface_j_check[62] = 5;
   yface_j_check[63] = 5;
   yface_j_check[64] = 6;
   yface_j_check[65] = 6;
   yface_j_check[66] = 7;
   yface_j_check[67] = 8;
   yface_j_check[68] = 8;
   yface_j_check[69] = 7;
   yface_j_check[70] = 6;
   yface_j_check[71] = 6;
   yface_j_check[72] = 6;
   yface_j_check[73] = 5;
   yface_j_check[74] = 5;
   yface_j_check[75] = 5;
   yface_j_check[76] = 6;
   yface_j_check[77] = 6;
   yface_j_check[78] = 5;
   yface_j_check[79] = 5;
   yface_j_check[80] = 6;
   yface_j_check[81] = 4;
   yface_j_check[82] = 8;
   yface_j_check[83] = 7;
   yface_j_check[84] = 7;
   yface_j_check[85] = 8;
   yface_j_check[86] = 8;
   yface_j_check[87] = 8;
   yface_j_check[88] = 7;
   yface_j_check[89] = 7;
   yface_j_check[90] = 8;
   yface_j_check[91] = 8;
   yface_j_check[92] = 7;
   yface_j_check[93] = 7;
   yface_j_check[94] = 6;
   yface_j_check[95] = 5;
   yface_j_check[96] = 5;
   yface_j_check[97] = 6;
   yface_j_check[98] = 6;
   yface_j_check[99] = 5;
   yface_j_check[100] = 5;
   yface_j_check[101] = 6;
   yface_j_check[102] = 6;
   yface_j_check[103] = 6;
   yface_j_check[104] = 5;
   yface_j_check[105] = 5;
   yface_j_check[106] = 6;
   yface_j_check[107] = 5;
   yface_j_check[108] = 6;
   yface_j_check[109] = 4;

   yface_level_check[0] = 0;
   yface_level_check[1] = 1;
   yface_level_check[2] = 1;
   yface_level_check[3] = 1;
   yface_level_check[4] = 1;
   yface_level_check[5] = 1;
   yface_level_check[6] = 1;
   yface_level_check[7] = 1;
   yface_level_check[8] = 1;
   yface_level_check[9] = 1;
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
   yface_level_check[22] = 2;
   yface_level_check[23] = 2;
   yface_level_check[24] = 2;
   yface_level_check[25] = 2;
   yface_level_check[26] = 1;
   yface_level_check[27] = 1;
   yface_level_check[28] = 1;
   yface_level_check[29] = 1;
   yface_level_check[30] = 1;
   yface_level_check[31] = 1;
   yface_level_check[32] = 1;
   yface_level_check[33] = 1;
   yface_level_check[34] = 2;
   yface_level_check[35] = 2;
   yface_level_check[36] = 1;
   yface_level_check[37] = 1;
   yface_level_check[38] = 0;
   yface_level_check[39] = 1;
   yface_level_check[40] = 1;
   yface_level_check[41] = 1;
   yface_level_check[42] = 1;
   yface_level_check[43] = 1;
   yface_level_check[44] = 1;
   yface_level_check[45] = 1;
   yface_level_check[46] = 1;
   yface_level_check[47] = 1;
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
   yface_level_check[58] = 2;
   yface_level_check[59] = 2;
   yface_level_check[60] = 2;
   yface_level_check[61] = 2;
   yface_level_check[62] = 2;
   yface_level_check[63] = 2;
   yface_level_check[64] = 2;
   yface_level_check[65] = 2;
   yface_level_check[66] = 2;
   yface_level_check[67] = 2;
   yface_level_check[68] = 2;
   yface_level_check[69] = 2;
   yface_level_check[70] = 1;
   yface_level_check[71] = 2;
   yface_level_check[72] = 2;
   yface_level_check[73] = 2;
   yface_level_check[74] = 2;
   yface_level_check[75] = 1;
   yface_level_check[76] = 1;
   yface_level_check[77] = 1;
   yface_level_check[78] = 1;
   yface_level_check[79] = 1;
   yface_level_check[80] = 1;
   yface_level_check[81] = 0;
   yface_level_check[82] = 1;
   yface_level_check[83] = 1;
   yface_level_check[84] = 1;
   yface_level_check[85] = 1;
   yface_level_check[86] = 1;
   yface_level_check[87] = 1;
   yface_level_check[88] = 1;
   yface_level_check[89] = 1;
   yface_level_check[90] = 2;
   yface_level_check[91] = 2;
   yface_level_check[92] = 2;
   yface_level_check[93] = 2;
   yface_level_check[94] = 2;
   yface_level_check[95] = 2;
   yface_level_check[96] = 2;
   yface_level_check[97] = 2;
   yface_level_check[98] = 2;
   yface_level_check[99] = 2;
   yface_level_check[100] = 2;
   yface_level_check[101] = 2;
   yface_level_check[102] = 1;
   yface_level_check[103] = 1;
   yface_level_check[104] = 1;
   yface_level_check[105] = 1;
   yface_level_check[106] = 1;
   yface_level_check[107] = 1;
   yface_level_check[108] = 1;
   yface_level_check[109] = 0;
}

void set_map_xface2cell_check(int *map_xface2cell_lower_check, int *map_xface2cell_upper_check)
{
   map_xface2cell_lower_check[0] = 0;
   map_xface2cell_lower_check[1] = 1;
   map_xface2cell_lower_check[2] = 1;
   map_xface2cell_lower_check[3] = 2;
   map_xface2cell_lower_check[4] = 3;
   map_xface2cell_lower_check[5] = 4;
   map_xface2cell_lower_check[6] = 5;
   map_xface2cell_lower_check[7] = 6;
   map_xface2cell_lower_check[8] = 7;
   map_xface2cell_lower_check[9] = 7;
   map_xface2cell_lower_check[10] = 8;
   map_xface2cell_lower_check[11] = 9;
   map_xface2cell_lower_check[12] = 9;
   map_xface2cell_lower_check[13] = 10;
   map_xface2cell_lower_check[14] = 11;
   map_xface2cell_lower_check[15] = 12;
   map_xface2cell_lower_check[16] = 13;
   map_xface2cell_lower_check[17] = 14;
   map_xface2cell_lower_check[18] = 15;
   map_xface2cell_lower_check[19] = 16;
   map_xface2cell_lower_check[20] = 17;
   map_xface2cell_lower_check[21] = 18;
   map_xface2cell_lower_check[22] = 19;
   map_xface2cell_lower_check[23] = 20;
   map_xface2cell_lower_check[24] = 21;
   map_xface2cell_lower_check[25] = 22;
   map_xface2cell_lower_check[26] = 23;
   map_xface2cell_lower_check[27] = 24;
   map_xface2cell_lower_check[28] = 25;
   map_xface2cell_lower_check[29] = 26;
   map_xface2cell_lower_check[30] = 27;
   map_xface2cell_lower_check[31] = 28;
   map_xface2cell_lower_check[32] = 29;
   map_xface2cell_lower_check[33] = 30;
   map_xface2cell_lower_check[34] = 31;
   map_xface2cell_lower_check[35] = 32;
   map_xface2cell_lower_check[36] = 33;
   map_xface2cell_lower_check[37] = 36;
   map_xface2cell_lower_check[38] = 39;
   map_xface2cell_lower_check[39] = 40;
   map_xface2cell_lower_check[40] = 41;
   map_xface2cell_lower_check[41] = 42;
   map_xface2cell_lower_check[42] = 43;
   map_xface2cell_lower_check[43] = 44;
   map_xface2cell_lower_check[44] = 45;
   map_xface2cell_lower_check[45] = 46;
   map_xface2cell_lower_check[46] = 47;
   map_xface2cell_lower_check[47] = 48;
   map_xface2cell_lower_check[48] = 49;
   map_xface2cell_lower_check[49] = 50;
   map_xface2cell_lower_check[50] = 51;
   map_xface2cell_lower_check[51] = 52;
   map_xface2cell_lower_check[52] = 53;
   map_xface2cell_lower_check[53] = 54;
   map_xface2cell_lower_check[54] = 55;
   map_xface2cell_lower_check[55] = 56;
   map_xface2cell_lower_check[56] = 57;
   map_xface2cell_lower_check[57] = 58;
   map_xface2cell_lower_check[58] = 59;
   map_xface2cell_lower_check[59] = 60;
   map_xface2cell_lower_check[60] = 61;
   map_xface2cell_lower_check[61] = 62;
   map_xface2cell_lower_check[62] = 63;
   map_xface2cell_lower_check[63] = 64;
   map_xface2cell_lower_check[64] = 65;
   map_xface2cell_lower_check[65] = 66;
   map_xface2cell_lower_check[66] = 67;
   map_xface2cell_lower_check[67] = 68;
   map_xface2cell_lower_check[68] = 69;
   map_xface2cell_lower_check[69] = 70;
   map_xface2cell_lower_check[70] = 71;
   map_xface2cell_lower_check[71] = 72;
   map_xface2cell_lower_check[72] = 75;
   map_xface2cell_lower_check[73] = 78;
   map_xface2cell_lower_check[74] = 79;
   map_xface2cell_lower_check[75] = 80;
   map_xface2cell_lower_check[76] = 81;
   map_xface2cell_lower_check[77] = 82;
   map_xface2cell_lower_check[78] = 83;
   map_xface2cell_lower_check[79] = 84;
   map_xface2cell_lower_check[80] = 85;
   map_xface2cell_lower_check[81] = 86;
   map_xface2cell_lower_check[82] = 87;
   map_xface2cell_lower_check[83] = 88;
   map_xface2cell_lower_check[84] = 89;
   map_xface2cell_lower_check[85] = 90;
   map_xface2cell_lower_check[86] = 91;
   map_xface2cell_lower_check[87] = 92;
   map_xface2cell_lower_check[88] = 93;
   map_xface2cell_lower_check[89] = 94;
   map_xface2cell_lower_check[90] = 95;
   map_xface2cell_lower_check[91] = 96;
   map_xface2cell_lower_check[92] = 97;
   map_xface2cell_lower_check[93] = 98;
   map_xface2cell_lower_check[94] = 99;
   map_xface2cell_lower_check[95] = 100;
   map_xface2cell_lower_check[96] = 101;
   map_xface2cell_lower_check[97] = 102;
   map_xface2cell_lower_check[98] = 102;
   map_xface2cell_lower_check[99] = 103;
   map_xface2cell_lower_check[100] = 104;
   map_xface2cell_lower_check[101] = 104;
   map_xface2cell_lower_check[102] = 105;
   map_xface2cell_lower_check[103] = 106;
   map_xface2cell_lower_check[104] = 107;
   map_xface2cell_lower_check[105] = 108;
   map_xface2cell_lower_check[106] = 109;
   map_xface2cell_lower_check[107] = 110;
   map_xface2cell_lower_check[108] = 110;
   map_xface2cell_lower_check[109] = 111;

   map_xface2cell_upper_check[0] = 26;
   map_xface2cell_upper_check[1] = 24;
   map_xface2cell_upper_check[2] = 23;
   map_xface2cell_upper_check[3] = 1;
   map_xface2cell_upper_check[4] = 5;
   map_xface2cell_upper_check[5] = 6;
   map_xface2cell_upper_check[6] = 8;
   map_xface2cell_upper_check[7] = 7;
   map_xface2cell_upper_check[8] = 10;
   map_xface2cell_upper_check[9] = 11;
   map_xface2cell_upper_check[10] = 9;
   map_xface2cell_upper_check[11] = 20;
   map_xface2cell_upper_check[12] = 19;
   map_xface2cell_upper_check[13] = 13;
   map_xface2cell_upper_check[14] = 12;
   map_xface2cell_upper_check[15] = 15;
   map_xface2cell_upper_check[16] = 14;
   map_xface2cell_upper_check[17] = 17;
   map_xface2cell_upper_check[18] = 16;
   map_xface2cell_upper_check[19] = 55;
   map_xface2cell_upper_check[20] = 52;
   map_xface2cell_upper_check[21] = 51;
   map_xface2cell_upper_check[22] = 18;
   map_xface2cell_upper_check[23] = 21;
   map_xface2cell_upper_check[24] = 50;
   map_xface2cell_upper_check[25] = 31;
   map_xface2cell_upper_check[26] = 22;
   map_xface2cell_upper_check[27] = 25;
   map_xface2cell_upper_check[28] = 30;
   map_xface2cell_upper_check[29] = 27;
   map_xface2cell_upper_check[30] = 28;
   map_xface2cell_upper_check[31] = 29;
   map_xface2cell_upper_check[32] = 34;
   map_xface2cell_upper_check[33] = 33;
   map_xface2cell_upper_check[34] = 32;
   map_xface2cell_upper_check[35] = 36;
   map_xface2cell_upper_check[36] = 36;
   map_xface2cell_upper_check[37] = 35;
   map_xface2cell_upper_check[38] = 38;
   map_xface2cell_upper_check[39] = 37;
   map_xface2cell_upper_check[40] = 40;
   map_xface2cell_upper_check[41] = 39;
   map_xface2cell_upper_check[42] = 42;
   map_xface2cell_upper_check[43] = 43;
   map_xface2cell_upper_check[44] = 46;
   map_xface2cell_upper_check[45] = 42;
   map_xface2cell_upper_check[46] = 41;
   map_xface2cell_upper_check[47] = 47;
   map_xface2cell_upper_check[48] = 47;
   map_xface2cell_upper_check[49] = 49;
   map_xface2cell_upper_check[50] = 48;
   map_xface2cell_upper_check[51] = 53;
   map_xface2cell_upper_check[52] = 45;
   map_xface2cell_upper_check[53] = 44;
   map_xface2cell_upper_check[54] = 54;
   map_xface2cell_upper_check[55] = 57;
   map_xface2cell_upper_check[56] = 67;
   map_xface2cell_upper_check[57] = 66;
   map_xface2cell_upper_check[58] = 58;
   map_xface2cell_upper_check[59] = 63;
   map_xface2cell_upper_check[60] = 62;
   map_xface2cell_upper_check[61] = 64;
   map_xface2cell_upper_check[62] = 64;
   map_xface2cell_upper_check[63] = 70;
   map_xface2cell_upper_check[64] = 69;
   map_xface2cell_upper_check[65] = 65;
   map_xface2cell_upper_check[66] = 68;
   map_xface2cell_upper_check[67] = 69;
   map_xface2cell_upper_check[68] = 72;
   map_xface2cell_upper_check[69] = 71;
   map_xface2cell_upper_check[70] = 74;
   map_xface2cell_upper_check[71] = 73;
   map_xface2cell_upper_check[72] = 76;
   map_xface2cell_upper_check[73] = 75;
   map_xface2cell_upper_check[74] = 75;
   map_xface2cell_upper_check[75] = 79;
   map_xface2cell_upper_check[76] = 78;
   map_xface2cell_upper_check[77] = 77;
   map_xface2cell_upper_check[78] = 82;
   map_xface2cell_upper_check[79] = 83;
   map_xface2cell_upper_check[80] = 84;
   map_xface2cell_upper_check[81] = 81;
   map_xface2cell_upper_check[82] = 86;
   map_xface2cell_upper_check[83] = 89;
   map_xface2cell_upper_check[84] = 80;
   map_xface2cell_upper_check[85] = 61;
   map_xface2cell_upper_check[86] = 90;
   map_xface2cell_upper_check[87] = 93;
   map_xface2cell_upper_check[88] = 60;
   map_xface2cell_upper_check[89] = 59;
   map_xface2cell_upper_check[90] = 56;
   map_xface2cell_upper_check[91] = 95;
   map_xface2cell_upper_check[92] = 94;
   map_xface2cell_upper_check[93] = 97;
   map_xface2cell_upper_check[94] = 96;
   map_xface2cell_upper_check[95] = 99;
   map_xface2cell_upper_check[96] = 98;
   map_xface2cell_upper_check[97] = 92;
   map_xface2cell_upper_check[98] = 91;
   map_xface2cell_upper_check[99] = 102;
   map_xface2cell_upper_check[100] = 100;
   map_xface2cell_upper_check[101] = 101;
   map_xface2cell_upper_check[102] = 104;
   map_xface2cell_upper_check[103] = 103;
   map_xface2cell_upper_check[104] = 105;
   map_xface2cell_upper_check[105] = 106;
   map_xface2cell_upper_check[106] = 110;
   map_xface2cell_upper_check[107] = 88;
   map_xface2cell_upper_check[108] = 87;
   map_xface2cell_upper_check[109] = 85;
}

void set_map_yface2cell_check(int *map_yface2cell_lower_check, int *map_yface2cell_upper_check)
{
   map_yface2cell_lower_check[0] = 0;
   map_yface2cell_lower_check[1] = 1;
   map_yface2cell_lower_check[2] = 1;
   map_yface2cell_lower_check[3] = 2;
   map_yface2cell_lower_check[4] = 3;
   map_yface2cell_lower_check[5] = 4;
   map_yface2cell_lower_check[6] = 5;
   map_yface2cell_lower_check[7] = 6;
   map_yface2cell_lower_check[8] = 7;
   map_yface2cell_lower_check[9] = 8;
   map_yface2cell_lower_check[10] = 9;
   map_yface2cell_lower_check[11] = 9;
   map_yface2cell_lower_check[12] = 10;
   map_yface2cell_lower_check[13] = 11;
   map_yface2cell_lower_check[14] = 12;
   map_yface2cell_lower_check[15] = 13;
   map_yface2cell_lower_check[16] = 14;
   map_yface2cell_lower_check[17] = 15;
   map_yface2cell_lower_check[18] = 16;
   map_yface2cell_lower_check[19] = 17;
   map_yface2cell_lower_check[20] = 18;
   map_yface2cell_lower_check[21] = 19;
   map_yface2cell_lower_check[22] = 20;
   map_yface2cell_lower_check[23] = 21;
   map_yface2cell_lower_check[24] = 22;
   map_yface2cell_lower_check[25] = 22;
   map_yface2cell_lower_check[26] = 23;
   map_yface2cell_lower_check[27] = 24;
   map_yface2cell_lower_check[28] = 25;
   map_yface2cell_lower_check[29] = 26;
   map_yface2cell_lower_check[30] = 27;
   map_yface2cell_lower_check[31] = 28;
   map_yface2cell_lower_check[32] = 29;
   map_yface2cell_lower_check[33] = 30;
   map_yface2cell_lower_check[34] = 31;
   map_yface2cell_lower_check[35] = 31;
   map_yface2cell_lower_check[36] = 32;
   map_yface2cell_lower_check[37] = 33;
   map_yface2cell_lower_check[38] = 34;
   map_yface2cell_lower_check[39] = 35;
   map_yface2cell_lower_check[40] = 36;
   map_yface2cell_lower_check[41] = 36;
   map_yface2cell_lower_check[42] = 37;
   map_yface2cell_lower_check[43] = 38;
   map_yface2cell_lower_check[44] = 39;
   map_yface2cell_lower_check[45] = 40;
   map_yface2cell_lower_check[46] = 41;
   map_yface2cell_lower_check[47] = 42;
   map_yface2cell_lower_check[48] = 43;
   map_yface2cell_lower_check[49] = 44;
   map_yface2cell_lower_check[50] = 45;
   map_yface2cell_lower_check[51] = 46;
   map_yface2cell_lower_check[52] = 47;
   map_yface2cell_lower_check[53] = 47;
   map_yface2cell_lower_check[54] = 48;
   map_yface2cell_lower_check[55] = 49;
   map_yface2cell_lower_check[56] = 50;
   map_yface2cell_lower_check[57] = 51;
   map_yface2cell_lower_check[58] = 52;
   map_yface2cell_lower_check[59] = 53;
   map_yface2cell_lower_check[60] = 54;
   map_yface2cell_lower_check[61] = 55;
   map_yface2cell_lower_check[62] = 56;
   map_yface2cell_lower_check[63] = 57;
   map_yface2cell_lower_check[64] = 58;
   map_yface2cell_lower_check[65] = 59;
   map_yface2cell_lower_check[66] = 60;
   map_yface2cell_lower_check[67] = 61;
   map_yface2cell_lower_check[68] = 62;
   map_yface2cell_lower_check[69] = 63;
   map_yface2cell_lower_check[70] = 64;
   map_yface2cell_lower_check[71] = 65;
   map_yface2cell_lower_check[72] = 66;
   map_yface2cell_lower_check[73] = 67;
   map_yface2cell_lower_check[74] = 68;
   map_yface2cell_lower_check[75] = 69;
   map_yface2cell_lower_check[76] = 70;
   map_yface2cell_lower_check[77] = 71;
   map_yface2cell_lower_check[78] = 72;
   map_yface2cell_lower_check[79] = 73;
   map_yface2cell_lower_check[80] = 74;
   map_yface2cell_lower_check[81] = 75;
   map_yface2cell_lower_check[82] = 78;
   map_yface2cell_lower_check[83] = 79;
   map_yface2cell_lower_check[84] = 80;
   map_yface2cell_lower_check[85] = 81;
   map_yface2cell_lower_check[86] = 86;
   map_yface2cell_lower_check[87] = 87;
   map_yface2cell_lower_check[88] = 88;
   map_yface2cell_lower_check[89] = 89;
   map_yface2cell_lower_check[90] = 90;
   map_yface2cell_lower_check[91] = 91;
   map_yface2cell_lower_check[92] = 92;
   map_yface2cell_lower_check[93] = 93;
   map_yface2cell_lower_check[94] = 94;
   map_yface2cell_lower_check[95] = 95;
   map_yface2cell_lower_check[96] = 96;
   map_yface2cell_lower_check[97] = 97;
   map_yface2cell_lower_check[98] = 98;
   map_yface2cell_lower_check[99] = 99;
   map_yface2cell_lower_check[100] = 100;
   map_yface2cell_lower_check[101] = 101;
   map_yface2cell_lower_check[102] = 102;
   map_yface2cell_lower_check[103] = 103;
   map_yface2cell_lower_check[104] = 104;
   map_yface2cell_lower_check[105] = 105;
   map_yface2cell_lower_check[106] = 106;
   map_yface2cell_lower_check[107] = 107;
   map_yface2cell_lower_check[108] = 108;
   map_yface2cell_lower_check[109] = 110;

   map_yface2cell_upper_check[0] = 1;
   map_yface2cell_upper_check[1] = 5;
   map_yface2cell_upper_check[2] = 8;
   map_yface2cell_upper_check[3] = 3;
   map_yface2cell_upper_check[4] = 4;
   map_yface2cell_upper_check[5] = 107;
   map_yface2cell_upper_check[6] = 6;
   map_yface2cell_upper_check[7] = 105;
   map_yface2cell_upper_check[8] = 104;
   map_yface2cell_upper_check[9] = 7;
   map_yface2cell_upper_check[10] = 10;
   map_yface2cell_upper_check[11] = 13;
   map_yface2cell_upper_check[12] = 11;
   map_yface2cell_upper_check[13] = 100;
   map_yface2cell_upper_check[14] = 99;
   map_yface2cell_upper_check[15] = 12;
   map_yface2cell_upper_check[16] = 15;
   map_yface2cell_upper_check[17] = 96;
   map_yface2cell_upper_check[18] = 95;
   map_yface2cell_upper_check[19] = 16;
   map_yface2cell_upper_check[20] = 17;
   map_yface2cell_upper_check[21] = 14;
   map_yface2cell_upper_check[22] = 19;
   map_yface2cell_upper_check[23] = 18;
   map_yface2cell_upper_check[24] = 20;
   map_yface2cell_upper_check[25] = 21;
   map_yface2cell_upper_check[26] = 9;
   map_yface2cell_upper_check[27] = 23;
   map_yface2cell_upper_check[28] = 22;
   map_yface2cell_upper_check[29] = 24;
   map_yface2cell_upper_check[30] = 25;
   map_yface2cell_upper_check[31] = 30;
   map_yface2cell_upper_check[32] = 33;
   map_yface2cell_upper_check[33] = 31;
   map_yface2cell_upper_check[34] = 50;
   map_yface2cell_upper_check[35] = 49;
   map_yface2cell_upper_check[36] = 47;
   map_yface2cell_upper_check[37] = 32;
   map_yface2cell_upper_check[38] = 36;
   map_yface2cell_upper_check[39] = 37;
   map_yface2cell_upper_check[40] = 41;
   map_yface2cell_upper_check[41] = 40;
   map_yface2cell_upper_check[42] = 38;
   map_yface2cell_upper_check[43] = 73;
   map_yface2cell_upper_check[44] = 72;
   map_yface2cell_upper_check[45] = 39;
   map_yface2cell_upper_check[46] = 42;
   map_yface2cell_upper_check[47] = 69;
   map_yface2cell_upper_check[48] = 68;
   map_yface2cell_upper_check[49] = 67;
   map_yface2cell_upper_check[50] = 44;
   map_yface2cell_upper_check[51] = 43;
   map_yface2cell_upper_check[52] = 45;
   map_yface2cell_upper_check[53] = 46;
   map_yface2cell_upper_check[54] = 53;
   map_yface2cell_upper_check[55] = 48;
   map_yface2cell_upper_check[56] = 51;
   map_yface2cell_upper_check[57] = 52;
   map_yface2cell_upper_check[58] = 55;
   map_yface2cell_upper_check[59] = 54;
   map_yface2cell_upper_check[60] = 57;
   map_yface2cell_upper_check[61] = 56;
   map_yface2cell_upper_check[62] = 59;
   map_yface2cell_upper_check[63] = 58;
   map_yface2cell_upper_check[64] = 63;
   map_yface2cell_upper_check[65] = 60;
   map_yface2cell_upper_check[66] = 61;
   map_yface2cell_upper_check[67] = 80;
   map_yface2cell_upper_check[68] = 80;
   map_yface2cell_upper_check[69] = 62;
   map_yface2cell_upper_check[70] = 79;
   map_yface2cell_upper_check[71] = 64;
   map_yface2cell_upper_check[72] = 64;
   map_yface2cell_upper_check[73] = 66;
   map_yface2cell_upper_check[74] = 65;
   map_yface2cell_upper_check[75] = 70;
   map_yface2cell_upper_check[76] = 75;
   map_yface2cell_upper_check[77] = 75;
   map_yface2cell_upper_check[78] = 71;
   map_yface2cell_upper_check[79] = 74;
   map_yface2cell_upper_check[80] = 76;
   map_yface2cell_upper_check[81] = 77;
   map_yface2cell_upper_check[82] = 82;
   map_yface2cell_upper_check[83] = 78;
   map_yface2cell_upper_check[84] = 81;
   map_yface2cell_upper_check[85] = 83;
   map_yface2cell_upper_check[86] = 84;
   map_yface2cell_upper_check[87] = 85;
   map_yface2cell_upper_check[88] = 87;
   map_yface2cell_upper_check[89] = 86;
   map_yface2cell_upper_check[90] = 89;
   map_yface2cell_upper_check[91] = 89;
   map_yface2cell_upper_check[92] = 91;
   map_yface2cell_upper_check[93] = 90;
   map_yface2cell_upper_check[94] = 93;
   map_yface2cell_upper_check[95] = 94;
   map_yface2cell_upper_check[96] = 97;
   map_yface2cell_upper_check[97] = 92;
   map_yface2cell_upper_check[98] = 102;
   map_yface2cell_upper_check[99] = 98;
   map_yface2cell_upper_check[100] = 101;
   map_yface2cell_upper_check[101] = 102;
   map_yface2cell_upper_check[102] = 88;
   map_yface2cell_upper_check[103] = 110;
   map_yface2cell_upper_check[104] = 103;
   map_yface2cell_upper_check[105] = 106;
   map_yface2cell_upper_check[106] = 110;
   map_yface2cell_upper_check[107] = 108;
   map_yface2cell_upper_check[108] = 109;
   map_yface2cell_upper_check[109] = 111;
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
      xface_flag_check[9][2] = 1;
      xface_flag_check[9][3] = 1;
      xface_flag_check[9][4] = 1;
      xface_flag_check[9][5] = 1;
      xface_flag_check[9][6] = 1;
      xface_flag_check[8][2] = 1;
      xface_flag_check[8][3] = 1;
      xface_flag_check[8][4] = 1;
      xface_flag_check[8][5] = 1;
      xface_flag_check[8][6] = 1;
      xface_flag_check[7][2] = 1;
      xface_flag_check[7][3] = 1;
      xface_flag_check[7][4] = 1;
      xface_flag_check[7][5] = 1;
      xface_flag_check[7][6] = 1;
      xface_flag_check[6][0] = 1;
      xface_flag_check[6][1] = 1;
      xface_flag_check[6][2] = 1;
      xface_flag_check[6][6] = 1;
      xface_flag_check[6][7] = 1;
      xface_flag_check[6][8] = 1;
      xface_flag_check[5][0] = 1;
      xface_flag_check[5][1] = 1;
      xface_flag_check[5][7] = 1;
      xface_flag_check[5][8] = 1;
      xface_flag_check[4][0] = 1;
      xface_flag_check[4][1] = 1;
      xface_flag_check[4][7] = 1;
      xface_flag_check[4][8] = 1;
      xface_flag_check[3][0] = 1;
      xface_flag_check[3][1] = 1;
      xface_flag_check[3][2] = 1;
      xface_flag_check[3][6] = 1;
      xface_flag_check[3][7] = 1;
      xface_flag_check[3][8] = 1;
      xface_flag_check[2][2] = 1;
      xface_flag_check[2][3] = 1;
      xface_flag_check[2][4] = 1;
      xface_flag_check[2][5] = 1;
      xface_flag_check[2][6] = 1;
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
   } else if (lev == 2){
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
      yface_flag_check[8][3] = 1;
      yface_flag_check[8][4] = 1;
      yface_flag_check[8][5] = 1;
      yface_flag_check[8][6] = 1;
      yface_flag_check[7][3] = 1;
      yface_flag_check[7][4] = 1;
      yface_flag_check[7][5] = 1;
      yface_flag_check[7][6] = 1;
      yface_flag_check[6][0] = 1;
      yface_flag_check[6][1] = 1;
      yface_flag_check[6][2] = 1;
      yface_flag_check[6][3] = 1;
      yface_flag_check[6][6] = 1;
      yface_flag_check[6][7] = 1;
      yface_flag_check[6][8] = 1;
      yface_flag_check[6][9] = 1;
      yface_flag_check[5][0] = 1;
      yface_flag_check[5][1] = 1;
      yface_flag_check[5][2] = 1;
      yface_flag_check[5][7] = 1;
      yface_flag_check[5][8] = 1;
      yface_flag_check[5][9] = 1;
      yface_flag_check[4][0] = 1;
      yface_flag_check[4][1] = 1;
      yface_flag_check[4][2] = 1;
      yface_flag_check[4][7] = 1;
      yface_flag_check[4][8] = 1;
      yface_flag_check[4][9] = 1;
      yface_flag_check[3][0] = 1;
      yface_flag_check[3][1] = 1;
      yface_flag_check[3][2] = 1;
      yface_flag_check[3][7] = 1;
      yface_flag_check[3][8] = 1;
      yface_flag_check[3][9] = 1;
      yface_flag_check[2][0] = 1;
      yface_flag_check[2][1] = 1;
      yface_flag_check[2][2] = 1;
      yface_flag_check[2][3] = 1;
      yface_flag_check[2][6] = 1;
      yface_flag_check[2][7] = 1;
      yface_flag_check[2][8] = 1;
      yface_flag_check[2][9] = 1;
      yface_flag_check[1][3] = 1;
      yface_flag_check[1][4] = 1;
      yface_flag_check[1][5] = 1;
      yface_flag_check[1][6] = 1;
      yface_flag_check[0][3] = 1;
      yface_flag_check[0][4] = 1;
      yface_flag_check[0][5] = 1;
      yface_flag_check[0][6] = 1;
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
      zone_flag_check[10][3] = 1;
      zone_flag_check[10][4] = 1;
      zone_flag_check[10][5] = 1;
      zone_flag_check[10][6] = 1;
      zone_flag_check[10][7] = 1;
      zone_flag_check[10][8] = 1;
      zone_flag_check[9][3] = 1;
      zone_flag_check[9][4] = 1;
      zone_flag_check[9][5] = 1;
      zone_flag_check[9][6] = 1;
      zone_flag_check[9][7] = 1;
      zone_flag_check[9][8] = 1;
      zone_flag_check[8][1] = 1;
      zone_flag_check[8][2] = 1;
      zone_flag_check[8][3] = 1;
      zone_flag_check[8][4] = 1;
      zone_flag_check[8][5] = 1;
      zone_flag_check[8][6] = 1;
      zone_flag_check[8][7] = 1;
      zone_flag_check[8][8] = 1;
      zone_flag_check[8][9] = 1;
      zone_flag_check[8][10] = 1;
      zone_flag_check[7][1] = 1;
      zone_flag_check[7][2] = 1;
      zone_flag_check[7][3] = 1;
      zone_flag_check[7][4] = 1;
      zone_flag_check[7][7] = 1;
      zone_flag_check[7][8] = 1;
      zone_flag_check[7][9] = 1;
      zone_flag_check[7][10] = 1;
      zone_flag_check[6][1] = 1;
      zone_flag_check[6][2] = 1;
      zone_flag_check[6][3] = 1;
      zone_flag_check[6][8] = 1;
      zone_flag_check[6][9] = 1;
      zone_flag_check[6][10] = 1;
      zone_flag_check[5][1] = 1;
      zone_flag_check[5][2] = 1;
      zone_flag_check[5][3] = 1;
      zone_flag_check[5][8] = 1;
      zone_flag_check[5][9] = 1;
      zone_flag_check[5][10] = 1;
      zone_flag_check[4][1] = 1;
      zone_flag_check[4][2] = 1;
      zone_flag_check[4][3] = 1;
      zone_flag_check[4][4] = 1;
      zone_flag_check[4][7] = 1;
      zone_flag_check[4][8] = 1;
      zone_flag_check[4][9] = 1;
      zone_flag_check[4][10] = 1;
      zone_flag_check[3][1] = 1;
      zone_flag_check[3][2] = 1;
      zone_flag_check[3][3] = 1;
      zone_flag_check[3][4] = 1;
      zone_flag_check[3][5] = 1;
      zone_flag_check[3][6] = 1;
      zone_flag_check[3][7] = 1;
      zone_flag_check[3][8] = 1;
      zone_flag_check[3][9] = 1;
      zone_flag_check[3][10] = 1;
      zone_flag_check[2][3] = 1;
      zone_flag_check[2][4] = 1;
      zone_flag_check[2][5] = 1;
      zone_flag_check[2][6] = 1;
      zone_flag_check[2][7] = 1;
      zone_flag_check[2][8] = 1;
      zone_flag_check[1][3] = 1;
      zone_flag_check[1][4] = 1;
      zone_flag_check[1][5] = 1;
      zone_flag_check[1][6] = 1;
      zone_flag_check[1][7] = 1;
      zone_flag_check[1][8] = 1;
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
      zone_cell_check[11][4] = 85;
      zone_cell_check[11][5] = 84;
      zone_cell_check[11][6] = 83;
      zone_cell_check[11][7] = 82;
      zone_cell_check[10][2] = 111;
      zone_cell_check[10][3] = 111;
      zone_cell_check[10][4] = 85;
      zone_cell_check[10][5] = 84;
      zone_cell_check[10][6] = 83;
      zone_cell_check[10][7] = 82;
      zone_cell_check[10][8] = 77;
      zone_cell_check[10][9] = 77;
      zone_cell_check[9][1] = 109;
      zone_cell_check[9][2] = 110;
      zone_cell_check[9][3] = 110;
      zone_cell_check[9][4] = 87;
      zone_cell_check[9][5] = 86;
      zone_cell_check[9][6] = 81;
      zone_cell_check[9][7] = 78;
      zone_cell_check[9][8] = 75;
      zone_cell_check[9][9] = 75;
      zone_cell_check[9][10] = 76;
      zone_cell_check[8][1] = 109;
      zone_cell_check[8][2] = 110;
      zone_cell_check[8][3] = 110;
      zone_cell_check[8][4] = 88;
      zone_cell_check[8][5] = 89;
      zone_cell_check[8][6] = 80;
      zone_cell_check[8][7] = 79;
      zone_cell_check[8][8] = 75;
      zone_cell_check[8][9] = 75;
      zone_cell_check[8][10] = 76;
      zone_cell_check[7][0] = 108;
      zone_cell_check[7][1] = 108;
      zone_cell_check[7][2] = 106;
      zone_cell_check[7][3] = 103;
      zone_cell_check[7][4] = 102;
      zone_cell_check[7][5] = 91;
      zone_cell_check[7][6] = 61;
      zone_cell_check[7][7] = 64;
      zone_cell_check[7][8] = 70;
      zone_cell_check[7][9] = 71;
      zone_cell_check[7][10] = 74;
      zone_cell_check[7][11] = 74;
      zone_cell_check[6][0] = 107;
      zone_cell_check[6][1] = 107;
      zone_cell_check[6][2] = 105;
      zone_cell_check[6][3] = 104;
      zone_cell_check[6][4] = 101;
      zone_cell_check[6][7] = 66;
      zone_cell_check[6][8] = 69;
      zone_cell_check[6][9] = 72;
      zone_cell_check[6][10] = 73;
      zone_cell_check[6][11] = 73;
      zone_cell_check[5][0] = 4;
      zone_cell_check[5][1] = 4;
      zone_cell_check[5][2] = 6;
      zone_cell_check[5][3] = 7;
      zone_cell_check[5][4] = 10;
      zone_cell_check[5][7] = 45;
      zone_cell_check[5][8] = 42;
      zone_cell_check[5][9] = 39;
      zone_cell_check[5][10] = 38;
      zone_cell_check[5][11] = 38;
      zone_cell_check[4][0] = 3;
      zone_cell_check[4][1] = 3;
      zone_cell_check[4][2] = 5;
      zone_cell_check[4][3] = 8;
      zone_cell_check[4][4] = 9;
      zone_cell_check[4][5] = 20;
      zone_cell_check[4][6] = 50;
      zone_cell_check[4][7] = 47;
      zone_cell_check[4][8] = 41;
      zone_cell_check[4][9] = 40;
      zone_cell_check[4][10] = 37;
      zone_cell_check[4][11] = 37;
      zone_cell_check[3][1] = 2;
      zone_cell_check[3][2] = 1;
      zone_cell_check[3][3] = 1;
      zone_cell_check[3][4] = 23;
      zone_cell_check[3][5] = 22;
      zone_cell_check[3][6] = 31;
      zone_cell_check[3][7] = 32;
      zone_cell_check[3][8] = 36;
      zone_cell_check[3][9] = 36;
      zone_cell_check[3][10] = 35;
      zone_cell_check[2][1] = 2;
      zone_cell_check[2][2] = 1;
      zone_cell_check[2][3] = 1;
      zone_cell_check[2][4] = 24;
      zone_cell_check[2][5] = 25;
      zone_cell_check[2][6] = 30;
      zone_cell_check[2][7] = 33;
      zone_cell_check[2][8] = 36;
      zone_cell_check[2][9] = 36;
      zone_cell_check[2][10] = 35;
      zone_cell_check[1][2] = 0;
      zone_cell_check[1][3] = 0;
      zone_cell_check[1][4] = 26;
      zone_cell_check[1][5] = 27;
      zone_cell_check[1][6] = 28;
      zone_cell_check[1][7] = 29;
      zone_cell_check[1][8] = 34;
      zone_cell_check[1][9] = 34;
      zone_cell_check[0][4] = 26;
      zone_cell_check[0][5] = 27;
      zone_cell_check[0][6] = 28;
      zone_cell_check[0][7] = 29;
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
   map_xcell2face_left1_check[1] = 3;
   map_xcell2face_left1_check[2] = -1;
   map_xcell2face_left1_check[3] = -1;
   map_xcell2face_left1_check[4] = -1;
   map_xcell2face_left1_check[5] = 4;
   map_xcell2face_left1_check[6] = 5;
   map_xcell2face_left1_check[7] = 7;
   map_xcell2face_left1_check[8] = 6;
   map_xcell2face_left1_check[9] = 10;
   map_xcell2face_left1_check[10] = 8;
   map_xcell2face_left1_check[11] = 9;
   map_xcell2face_left1_check[12] = 14;
   map_xcell2face_left1_check[13] = 13;
   map_xcell2face_left1_check[14] = 16;
   map_xcell2face_left1_check[15] = 15;
   map_xcell2face_left1_check[16] = 18;
   map_xcell2face_left1_check[17] = 17;
   map_xcell2face_left1_check[18] = 22;
   map_xcell2face_left1_check[19] = 12;
   map_xcell2face_left1_check[20] = 11;
   map_xcell2face_left1_check[21] = 23;
   map_xcell2face_left1_check[22] = 26;
   map_xcell2face_left1_check[23] = 2;
   map_xcell2face_left1_check[24] = 1;
   map_xcell2face_left1_check[25] = 27;
   map_xcell2face_left1_check[26] = -1;
   map_xcell2face_left1_check[27] = 29;
   map_xcell2face_left1_check[28] = 30;
   map_xcell2face_left1_check[29] = 31;
   map_xcell2face_left1_check[30] = 28;
   map_xcell2face_left1_check[31] = 25;
   map_xcell2face_left1_check[32] = 34;
   map_xcell2face_left1_check[33] = 33;
   map_xcell2face_left1_check[34] = 32;
   map_xcell2face_left1_check[35] = 37;
   map_xcell2face_left1_check[36] = 36;
   map_xcell2face_left1_check[37] = 39;
   map_xcell2face_left1_check[38] = 38;
   map_xcell2face_left1_check[39] = 41;
   map_xcell2face_left1_check[40] = 40;
   map_xcell2face_left1_check[41] = 46;
   map_xcell2face_left1_check[42] = 45;
   map_xcell2face_left1_check[43] = 43;
   map_xcell2face_left1_check[44] = 53;
   map_xcell2face_left1_check[45] = 52;
   map_xcell2face_left1_check[46] = 44;
   map_xcell2face_left1_check[47] = 48;
   map_xcell2face_left1_check[48] = 50;
   map_xcell2face_left1_check[49] = 49;
   map_xcell2face_left1_check[50] = 24;
   map_xcell2face_left1_check[51] = 21;
   map_xcell2face_left1_check[52] = 20;
   map_xcell2face_left1_check[53] = 51;
   map_xcell2face_left1_check[54] = 54;
   map_xcell2face_left1_check[55] = 19;
   map_xcell2face_left1_check[56] = 90;
   map_xcell2face_left1_check[57] = 55;
   map_xcell2face_left1_check[58] = 58;
   map_xcell2face_left1_check[59] = 89;
   map_xcell2face_left1_check[60] = 88;
   map_xcell2face_left1_check[61] = 85;
   map_xcell2face_left1_check[62] = 60;
   map_xcell2face_left1_check[63] = 59;
   map_xcell2face_left1_check[64] = 62;
   map_xcell2face_left1_check[65] = 65;
   map_xcell2face_left1_check[66] = 57;
   map_xcell2face_left1_check[67] = 56;
   map_xcell2face_left1_check[68] = 66;
   map_xcell2face_left1_check[69] = 67;
   map_xcell2face_left1_check[70] = 63;
   map_xcell2face_left1_check[71] = 69;
   map_xcell2face_left1_check[72] = 68;
   map_xcell2face_left1_check[73] = 71;
   map_xcell2face_left1_check[74] = 70;
   map_xcell2face_left1_check[75] = 74;
   map_xcell2face_left1_check[76] = 72;
   map_xcell2face_left1_check[77] = 77;
   map_xcell2face_left1_check[78] = 76;
   map_xcell2face_left1_check[79] = 75;
   map_xcell2face_left1_check[80] = 84;
   map_xcell2face_left1_check[81] = 81;
   map_xcell2face_left1_check[82] = 78;
   map_xcell2face_left1_check[83] = 79;
   map_xcell2face_left1_check[84] = 80;
   map_xcell2face_left1_check[85] = 109;
   map_xcell2face_left1_check[86] = 82;
   map_xcell2face_left1_check[87] = 108;
   map_xcell2face_left1_check[88] = 107;
   map_xcell2face_left1_check[89] = 83;
   map_xcell2face_left1_check[90] = 86;
   map_xcell2face_left1_check[91] = 98;
   map_xcell2face_left1_check[92] = 97;
   map_xcell2face_left1_check[93] = 87;
   map_xcell2face_left1_check[94] = 92;
   map_xcell2face_left1_check[95] = 91;
   map_xcell2face_left1_check[96] = 94;
   map_xcell2face_left1_check[97] = 93;
   map_xcell2face_left1_check[98] = 96;
   map_xcell2face_left1_check[99] = 95;
   map_xcell2face_left1_check[100] = 100;
   map_xcell2face_left1_check[101] = 101;
   map_xcell2face_left1_check[102] = 99;
   map_xcell2face_left1_check[103] = 103;
   map_xcell2face_left1_check[104] = 102;
   map_xcell2face_left1_check[105] = 104;
   map_xcell2face_left1_check[106] = 105;
   map_xcell2face_left1_check[107] = -1;
   map_xcell2face_left1_check[108] = -1;
   map_xcell2face_left1_check[109] = -1;
   map_xcell2face_left1_check[110] = 106;
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
   map_xcell2face_left2_check[34] = 36;
   map_xcell2face_left2_check[35] = -1;
   map_xcell2face_left2_check[36] = 35;
   map_xcell2face_left2_check[37] = -1;
   map_xcell2face_left2_check[38] = -1;
   map_xcell2face_left2_check[39] = -1;
   map_xcell2face_left2_check[40] = -1;
   map_xcell2face_left2_check[41] = -1;
   map_xcell2face_left2_check[42] = 42;
   map_xcell2face_left2_check[43] = -1;
   map_xcell2face_left2_check[44] = -1;
   map_xcell2face_left2_check[45] = -1;
   map_xcell2face_left2_check[46] = -1;
   map_xcell2face_left2_check[47] = 47;
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
   map_xcell2face_left2_check[64] = 61;
   map_xcell2face_left2_check[65] = -1;
   map_xcell2face_left2_check[66] = -1;
   map_xcell2face_left2_check[67] = -1;
   map_xcell2face_left2_check[68] = -1;
   map_xcell2face_left2_check[69] = 64;
   map_xcell2face_left2_check[70] = -1;
   map_xcell2face_left2_check[71] = -1;
   map_xcell2face_left2_check[72] = -1;
   map_xcell2face_left2_check[73] = -1;
   map_xcell2face_left2_check[74] = -1;
   map_xcell2face_left2_check[75] = 73;
   map_xcell2face_left2_check[76] = -1;
   map_xcell2face_left2_check[77] = 77;
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

   map_xcell2face_right1_check[0] = 0;
   map_xcell2face_right1_check[1] = 1;
   map_xcell2face_right1_check[2] = 3;
   map_xcell2face_right1_check[3] = 4;
   map_xcell2face_right1_check[4] = 5;
   map_xcell2face_right1_check[5] = 6;
   map_xcell2face_right1_check[6] = 7;
   map_xcell2face_right1_check[7] = 8;
   map_xcell2face_right1_check[8] = 10;
   map_xcell2face_right1_check[9] = 11;
   map_xcell2face_right1_check[10] = 13;
   map_xcell2face_right1_check[11] = 14;
   map_xcell2face_right1_check[12] = 15;
   map_xcell2face_right1_check[13] = 16;
   map_xcell2face_right1_check[14] = 17;
   map_xcell2face_right1_check[15] = 18;
   map_xcell2face_right1_check[16] = 19;
   map_xcell2face_right1_check[17] = 20;
   map_xcell2face_right1_check[18] = 21;
   map_xcell2face_right1_check[19] = 22;
   map_xcell2face_right1_check[20] = 23;
   map_xcell2face_right1_check[21] = 24;
   map_xcell2face_right1_check[22] = 25;
   map_xcell2face_right1_check[23] = 26;
   map_xcell2face_right1_check[24] = 27;
   map_xcell2face_right1_check[25] = 28;
   map_xcell2face_right1_check[26] = 29;
   map_xcell2face_right1_check[27] = 30;
   map_xcell2face_right1_check[28] = 31;
   map_xcell2face_right1_check[29] = 32;
   map_xcell2face_right1_check[30] = 33;
   map_xcell2face_right1_check[31] = 34;
   map_xcell2face_right1_check[32] = 35;
   map_xcell2face_right1_check[33] = 36;
   map_xcell2face_right1_check[34] = -1;
   map_xcell2face_right1_check[35] = -1;
   map_xcell2face_right1_check[36] = 37;
   map_xcell2face_right1_check[37] = -1;
   map_xcell2face_right1_check[38] = -1;
   map_xcell2face_right1_check[39] = 38;
   map_xcell2face_right1_check[40] = 39;
   map_xcell2face_right1_check[41] = 40;
   map_xcell2face_right1_check[42] = 41;
   map_xcell2face_right1_check[43] = 42;
   map_xcell2face_right1_check[44] = 43;
   map_xcell2face_right1_check[45] = 44;
   map_xcell2face_right1_check[46] = 45;
   map_xcell2face_right1_check[47] = 46;
   map_xcell2face_right1_check[48] = 47;
   map_xcell2face_right1_check[49] = 48;
   map_xcell2face_right1_check[50] = 49;
   map_xcell2face_right1_check[51] = 50;
   map_xcell2face_right1_check[52] = 51;
   map_xcell2face_right1_check[53] = 52;
   map_xcell2face_right1_check[54] = 53;
   map_xcell2face_right1_check[55] = 54;
   map_xcell2face_right1_check[56] = 55;
   map_xcell2face_right1_check[57] = 56;
   map_xcell2face_right1_check[58] = 57;
   map_xcell2face_right1_check[59] = 58;
   map_xcell2face_right1_check[60] = 59;
   map_xcell2face_right1_check[61] = 60;
   map_xcell2face_right1_check[62] = 61;
   map_xcell2face_right1_check[63] = 62;
   map_xcell2face_right1_check[64] = 63;
   map_xcell2face_right1_check[65] = 64;
   map_xcell2face_right1_check[66] = 65;
   map_xcell2face_right1_check[67] = 66;
   map_xcell2face_right1_check[68] = 67;
   map_xcell2face_right1_check[69] = 68;
   map_xcell2face_right1_check[70] = 69;
   map_xcell2face_right1_check[71] = 70;
   map_xcell2face_right1_check[72] = 71;
   map_xcell2face_right1_check[73] = -1;
   map_xcell2face_right1_check[74] = -1;
   map_xcell2face_right1_check[75] = 72;
   map_xcell2face_right1_check[76] = -1;
   map_xcell2face_right1_check[77] = -1;
   map_xcell2face_right1_check[78] = 73;
   map_xcell2face_right1_check[79] = 74;
   map_xcell2face_right1_check[80] = 75;
   map_xcell2face_right1_check[81] = 76;
   map_xcell2face_right1_check[82] = 77;
   map_xcell2face_right1_check[83] = 78;
   map_xcell2face_right1_check[84] = 79;
   map_xcell2face_right1_check[85] = 80;
   map_xcell2face_right1_check[86] = 81;
   map_xcell2face_right1_check[87] = 82;
   map_xcell2face_right1_check[88] = 83;
   map_xcell2face_right1_check[89] = 84;
   map_xcell2face_right1_check[90] = 85;
   map_xcell2face_right1_check[91] = 86;
   map_xcell2face_right1_check[92] = 87;
   map_xcell2face_right1_check[93] = 88;
   map_xcell2face_right1_check[94] = 89;
   map_xcell2face_right1_check[95] = 90;
   map_xcell2face_right1_check[96] = 91;
   map_xcell2face_right1_check[97] = 92;
   map_xcell2face_right1_check[98] = 93;
   map_xcell2face_right1_check[99] = 94;
   map_xcell2face_right1_check[100] = 95;
   map_xcell2face_right1_check[101] = 96;
   map_xcell2face_right1_check[102] = 97;
   map_xcell2face_right1_check[103] = 99;
   map_xcell2face_right1_check[104] = 100;
   map_xcell2face_right1_check[105] = 102;
   map_xcell2face_right1_check[106] = 103;
   map_xcell2face_right1_check[107] = 104;
   map_xcell2face_right1_check[108] = 105;
   map_xcell2face_right1_check[109] = 106;
   map_xcell2face_right1_check[110] = 107;
   map_xcell2face_right1_check[111] = 109;

   return(map_xcell2face_right1_check);
}

int *set_map_xcell2face_right2_check(int ncells)
{
   int *map_xcell2face_right2_check = (int *)malloc(ncells*sizeof(int));

   map_xcell2face_right2_check[0] = -1;
   map_xcell2face_right2_check[1] = 2;
   map_xcell2face_right2_check[2] = -1;
   map_xcell2face_right2_check[3] = -1;
   map_xcell2face_right2_check[4] = -1;
   map_xcell2face_right2_check[5] = -1;
   map_xcell2face_right2_check[6] = -1;
   map_xcell2face_right2_check[7] = 9;
   map_xcell2face_right2_check[8] = -1;
   map_xcell2face_right2_check[9] = 12;
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
   map_xcell2face_right2_check[102] = 98;
   map_xcell2face_right2_check[103] = -1;
   map_xcell2face_right2_check[104] = 101;
   map_xcell2face_right2_check[105] = -1;
   map_xcell2face_right2_check[106] = -1;
   map_xcell2face_right2_check[107] = -1;
   map_xcell2face_right2_check[108] = -1;
   map_xcell2face_right2_check[109] = -1;
   map_xcell2face_right2_check[110] = 108;
   map_xcell2face_right2_check[111] = -1;

   return(map_xcell2face_right2_check);
}

int *set_map_ycell2face_bot1_check(int ncells)
{
   int *map_ycell2face_bot1_check = (int *)malloc(ncells*sizeof(int));

   map_ycell2face_bot1_check[0] = -1;
   map_ycell2face_bot1_check[1] = 0;
   map_ycell2face_bot1_check[2] = -1;
   map_ycell2face_bot1_check[3] = -1;
   map_ycell2face_bot1_check[4] = 4;
   map_ycell2face_bot1_check[5] = 1;
   map_ycell2face_bot1_check[6] = 6;
   map_ycell2face_bot1_check[7] = 9;
   map_ycell2face_bot1_check[8] = 2;
   map_ycell2face_bot1_check[9] = 26;
   map_ycell2face_bot1_check[10] = 10;
   map_ycell2face_bot1_check[11] = 12;
   map_ycell2face_bot1_check[12] = 15;
   map_ycell2face_bot1_check[13] = 11;
   map_ycell2face_bot1_check[14] = 21;
   map_ycell2face_bot1_check[15] = 16;
   map_ycell2face_bot1_check[16] = 19;
   map_ycell2face_bot1_check[17] = 20;
   map_ycell2face_bot1_check[18] = 23;
   map_ycell2face_bot1_check[19] = 22;
   map_ycell2face_bot1_check[20] = 24;
   map_ycell2face_bot1_check[21] = 25;
   map_ycell2face_bot1_check[22] = 28;
   map_ycell2face_bot1_check[23] = 27;
   map_ycell2face_bot1_check[24] = 29;
   map_ycell2face_bot1_check[25] = 30;
   map_ycell2face_bot1_check[26] = -1;
   map_ycell2face_bot1_check[27] = -1;
   map_ycell2face_bot1_check[28] = -1;
   map_ycell2face_bot1_check[29] = -1;
   map_ycell2face_bot1_check[30] = 31;
   map_ycell2face_bot1_check[31] = 33;
   map_ycell2face_bot1_check[32] = 37;
   map_ycell2face_bot1_check[33] = 32;
   map_ycell2face_bot1_check[34] = -1;
   map_ycell2face_bot1_check[35] = -1;
   map_ycell2face_bot1_check[36] = 38;
   map_ycell2face_bot1_check[37] = 39;
   map_ycell2face_bot1_check[38] = 42;
   map_ycell2face_bot1_check[39] = 45;
   map_ycell2face_bot1_check[40] = 41;
   map_ycell2face_bot1_check[41] = 40;
   map_ycell2face_bot1_check[42] = 46;
   map_ycell2face_bot1_check[43] = 51;
   map_ycell2face_bot1_check[44] = 50;
   map_ycell2face_bot1_check[45] = 52;
   map_ycell2face_bot1_check[46] = 53;
   map_ycell2face_bot1_check[47] = 36;
   map_ycell2face_bot1_check[48] = 55;
   map_ycell2face_bot1_check[49] = 35;
   map_ycell2face_bot1_check[50] = 34;
   map_ycell2face_bot1_check[51] = 56;
   map_ycell2face_bot1_check[52] = 57;
   map_ycell2face_bot1_check[53] = 54;
   map_ycell2face_bot1_check[54] = 59;
   map_ycell2face_bot1_check[55] = 58;
   map_ycell2face_bot1_check[56] = 61;
   map_ycell2face_bot1_check[57] = 60;
   map_ycell2face_bot1_check[58] = 63;
   map_ycell2face_bot1_check[59] = 62;
   map_ycell2face_bot1_check[60] = 65;
   map_ycell2face_bot1_check[61] = 66;
   map_ycell2face_bot1_check[62] = 69;
   map_ycell2face_bot1_check[63] = 64;
   map_ycell2face_bot1_check[64] = 72;
   map_ycell2face_bot1_check[65] = 74;
   map_ycell2face_bot1_check[66] = 73;
   map_ycell2face_bot1_check[67] = 49;
   map_ycell2face_bot1_check[68] = 48;
   map_ycell2face_bot1_check[69] = 47;
   map_ycell2face_bot1_check[70] = 75;
   map_ycell2face_bot1_check[71] = 78;
   map_ycell2face_bot1_check[72] = 44;
   map_ycell2face_bot1_check[73] = 43;
   map_ycell2face_bot1_check[74] = 79;
   map_ycell2face_bot1_check[75] = 76;
   map_ycell2face_bot1_check[76] = 80;
   map_ycell2face_bot1_check[77] = 81;
   map_ycell2face_bot1_check[78] = 83;
   map_ycell2face_bot1_check[79] = 70;
   map_ycell2face_bot1_check[80] = 67;
   map_ycell2face_bot1_check[81] = 84;
   map_ycell2face_bot1_check[82] = 82;
   map_ycell2face_bot1_check[83] = 85;
   map_ycell2face_bot1_check[84] = 86;
   map_ycell2face_bot1_check[85] = 87;
   map_ycell2face_bot1_check[86] = 89;
   map_ycell2face_bot1_check[87] = 88;
   map_ycell2face_bot1_check[88] = 102;
   map_ycell2face_bot1_check[89] = 91;
   map_ycell2face_bot1_check[90] = 93;
   map_ycell2face_bot1_check[91] = 92;
   map_ycell2face_bot1_check[92] = 97;
   map_ycell2face_bot1_check[93] = 94;
   map_ycell2face_bot1_check[94] = 95;
   map_ycell2face_bot1_check[95] = 18;
   map_ycell2face_bot1_check[96] = 17;
   map_ycell2face_bot1_check[97] = 96;
   map_ycell2face_bot1_check[98] = 99;
   map_ycell2face_bot1_check[99] = 14;
   map_ycell2face_bot1_check[100] = 13;
   map_ycell2face_bot1_check[101] = 100;
   map_ycell2face_bot1_check[102] = 101;
   map_ycell2face_bot1_check[103] = 104;
   map_ycell2face_bot1_check[104] = 8;
   map_ycell2face_bot1_check[105] = 7;
   map_ycell2face_bot1_check[106] = 105;
   map_ycell2face_bot1_check[107] = 5;
   map_ycell2face_bot1_check[108] = 107;
   map_ycell2face_bot1_check[109] = 108;
   map_ycell2face_bot1_check[110] = 106;
   map_ycell2face_bot1_check[111] = 109;

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
   map_ycell2face_bot2_check[64] = 71;
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
   map_ycell2face_bot2_check[75] = 77;
   map_ycell2face_bot2_check[76] = 80;
   map_ycell2face_bot2_check[77] = -1;
   map_ycell2face_bot2_check[78] = -1;
   map_ycell2face_bot2_check[79] = -1;
   map_ycell2face_bot2_check[80] = 68;
   map_ycell2face_bot2_check[81] = -1;
   map_ycell2face_bot2_check[82] = -1;
   map_ycell2face_bot2_check[83] = -1;
   map_ycell2face_bot2_check[84] = -1;
   map_ycell2face_bot2_check[85] = -1;
   map_ycell2face_bot2_check[86] = -1;
   map_ycell2face_bot2_check[87] = -1;
   map_ycell2face_bot2_check[88] = -1;
   map_ycell2face_bot2_check[89] = 90;
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
   map_ycell2face_bot2_check[102] = 98;
   map_ycell2face_bot2_check[103] = -1;
   map_ycell2face_bot2_check[104] = -1;
   map_ycell2face_bot2_check[105] = -1;
   map_ycell2face_bot2_check[106] = -1;
   map_ycell2face_bot2_check[107] = -1;
   map_ycell2face_bot2_check[108] = -1;
   map_ycell2face_bot2_check[109] = 106;
   map_ycell2face_bot2_check[110] = 103;
   map_ycell2face_bot2_check[111] = -1;

   return(map_ycell2face_bot2_check);
}

int *set_map_ycell2face_top1_check(int ncells)
{
   int *map_ycell2face_top1_check = (int *)malloc(ncells*sizeof(int));

   map_ycell2face_top1_check[0] = 0;
   map_ycell2face_top1_check[1] = 1;
   map_ycell2face_top1_check[2] = 3;
   map_ycell2face_top1_check[3] = 4;
   map_ycell2face_top1_check[4] = 5;
   map_ycell2face_top1_check[5] = 6;
   map_ycell2face_top1_check[6] = 7;
   map_ycell2face_top1_check[7] = 8;
   map_ycell2face_top1_check[8] = 9;
   map_ycell2face_top1_check[9] = 10;
   map_ycell2face_top1_check[10] = 12;
   map_ycell2face_top1_check[11] = 13;
   map_ycell2face_top1_check[12] = 14;
   map_ycell2face_top1_check[13] = 15;
   map_ycell2face_top1_check[14] = 16;
   map_ycell2face_top1_check[15] = 17;
   map_ycell2face_top1_check[16] = 18;
   map_ycell2face_top1_check[17] = 19;
   map_ycell2face_top1_check[18] = 20;
   map_ycell2face_top1_check[19] = 21;
   map_ycell2face_top1_check[20] = 22;
   map_ycell2face_top1_check[21] = 23;
   map_ycell2face_top1_check[22] = 24;
   map_ycell2face_top1_check[23] = 26;
   map_ycell2face_top1_check[24] = 27;
   map_ycell2face_top1_check[25] = 28;
   map_ycell2face_top1_check[26] = 29;
   map_ycell2face_top1_check[27] = 30;
   map_ycell2face_top1_check[28] = 31;
   map_ycell2face_top1_check[29] = 32;
   map_ycell2face_top1_check[30] = 33;
   map_ycell2face_top1_check[31] = 34;
   map_ycell2face_top1_check[32] = 36;
   map_ycell2face_top1_check[33] = 37;
   map_ycell2face_top1_check[34] = 38;
   map_ycell2face_top1_check[35] = 39;
   map_ycell2face_top1_check[36] = 40;
   map_ycell2face_top1_check[37] = 42;
   map_ycell2face_top1_check[38] = 43;
   map_ycell2face_top1_check[39] = 44;
   map_ycell2face_top1_check[40] = 45;
   map_ycell2face_top1_check[41] = 46;
   map_ycell2face_top1_check[42] = 47;
   map_ycell2face_top1_check[43] = 48;
   map_ycell2face_top1_check[44] = 49;
   map_ycell2face_top1_check[45] = 50;
   map_ycell2face_top1_check[46] = 51;
   map_ycell2face_top1_check[47] = 52;
   map_ycell2face_top1_check[48] = 54;
   map_ycell2face_top1_check[49] = 55;
   map_ycell2face_top1_check[50] = 56;
   map_ycell2face_top1_check[51] = 57;
   map_ycell2face_top1_check[52] = 58;
   map_ycell2face_top1_check[53] = 59;
   map_ycell2face_top1_check[54] = 60;
   map_ycell2face_top1_check[55] = 61;
   map_ycell2face_top1_check[56] = 62;
   map_ycell2face_top1_check[57] = 63;
   map_ycell2face_top1_check[58] = 64;
   map_ycell2face_top1_check[59] = 65;
   map_ycell2face_top1_check[60] = 66;
   map_ycell2face_top1_check[61] = 67;
   map_ycell2face_top1_check[62] = 68;
   map_ycell2face_top1_check[63] = 69;
   map_ycell2face_top1_check[64] = 70;
   map_ycell2face_top1_check[65] = 71;
   map_ycell2face_top1_check[66] = 72;
   map_ycell2face_top1_check[67] = 73;
   map_ycell2face_top1_check[68] = 74;
   map_ycell2face_top1_check[69] = 75;
   map_ycell2face_top1_check[70] = 76;
   map_ycell2face_top1_check[71] = 77;
   map_ycell2face_top1_check[72] = 78;
   map_ycell2face_top1_check[73] = 79;
   map_ycell2face_top1_check[74] = 80;
   map_ycell2face_top1_check[75] = 81;
   map_ycell2face_top1_check[76] = -1;
   map_ycell2face_top1_check[77] = -1;
   map_ycell2face_top1_check[78] = 82;
   map_ycell2face_top1_check[79] = 83;
   map_ycell2face_top1_check[80] = 84;
   map_ycell2face_top1_check[81] = 85;
   map_ycell2face_top1_check[82] = -1;
   map_ycell2face_top1_check[83] = -1;
   map_ycell2face_top1_check[84] = -1;
   map_ycell2face_top1_check[85] = -1;
   map_ycell2face_top1_check[86] = 86;
   map_ycell2face_top1_check[87] = 87;
   map_ycell2face_top1_check[88] = 88;
   map_ycell2face_top1_check[89] = 89;
   map_ycell2face_top1_check[90] = 90;
   map_ycell2face_top1_check[91] = 91;
   map_ycell2face_top1_check[92] = 92;
   map_ycell2face_top1_check[93] = 93;
   map_ycell2face_top1_check[94] = 94;
   map_ycell2face_top1_check[95] = 95;
   map_ycell2face_top1_check[96] = 96;
   map_ycell2face_top1_check[97] = 97;
   map_ycell2face_top1_check[98] = 98;
   map_ycell2face_top1_check[99] = 99;
   map_ycell2face_top1_check[100] = 100;
   map_ycell2face_top1_check[101] = 101;
   map_ycell2face_top1_check[102] = 102;
   map_ycell2face_top1_check[103] = 103;
   map_ycell2face_top1_check[104] = 104;
   map_ycell2face_top1_check[105] = 105;
   map_ycell2face_top1_check[106] = 106;
   map_ycell2face_top1_check[107] = 107;
   map_ycell2face_top1_check[108] = 108;
   map_ycell2face_top1_check[109] = -1;
   map_ycell2face_top1_check[110] = 109;
   map_ycell2face_top1_check[111] = -1;

   return(map_ycell2face_top1_check);
}

int *set_map_ycell2face_top2_check(int ncells)
{
   int *map_ycell2face_top2_check = (int *)malloc(ncells*sizeof(int));

   map_ycell2face_top2_check[0] = -1;
   map_ycell2face_top2_check[1] = 2;
   map_ycell2face_top2_check[2] = -1;
   map_ycell2face_top2_check[3] = -1;
   map_ycell2face_top2_check[4] = -1;
   map_ycell2face_top2_check[5] = -1;
   map_ycell2face_top2_check[6] = -1;
   map_ycell2face_top2_check[7] = -1;
   map_ycell2face_top2_check[8] = -1;
   map_ycell2face_top2_check[9] = 11;
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
   map_ycell2face_top2_check[22] = 25;
   map_ycell2face_top2_check[23] = -1;
   map_ycell2face_top2_check[24] = -1;
   map_ycell2face_top2_check[25] = -1;
   map_ycell2face_top2_check[26] = -1;
   map_ycell2face_top2_check[27] = -1;
   map_ycell2face_top2_check[28] = -1;
   map_ycell2face_top2_check[29] = -1;
   map_ycell2face_top2_check[30] = -1;
   map_ycell2face_top2_check[31] = 35;
   map_ycell2face_top2_check[32] = -1;
   map_ycell2face_top2_check[33] = -1;
   map_ycell2face_top2_check[34] = -1;
   map_ycell2face_top2_check[35] = -1;
   map_ycell2face_top2_check[36] = 41;
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
   map_ycell2face_top2_check[47] = 53;
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
