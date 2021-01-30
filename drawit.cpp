#include <stdio.h>
#include <stdlib.h>
#include <string.h>
#include "graphics/display.h"
#include "graphics/graphics.h"

int main(int argc, char *argv[])
{
   char buffer[81];
   FILE *fp = fopen("data.txt", "r");

   int ncells;
   float xmin =  10000000.0;
   float xmax = -10000000.0;
   float ymin =  10000000.0;
   float ymax = -10000000.0;
   int ic = 0;
   while (fgets(buffer, sizeof(buffer), fp) !=NULL){ 
      ic++;
   }
   ncells = ic;
   rewind(fp);

   int *proc = (int *)malloc(sizeof(int)*ncells);
   int *cellnum = (int *)malloc(sizeof(int)*ncells);
   int *i = (int *)malloc(sizeof(int)*ncells);
   int *j = (int *)malloc(sizeof(int)*ncells);
   int *level = (int *)malloc(sizeof(int)*ncells);
   float *x = (float *)malloc(sizeof(float)*ncells);
   float *y = (float *)malloc(sizeof(float)*ncells);
   float *dx = (float *)malloc(sizeof(float)*ncells);
   float *dy = (float *)malloc(sizeof(float)*ncells);
   float *data = (float *)malloc(sizeof(float)*ncells);
   ic = 0;
   for (int ic = 0; ic < ncells; ic++){
      fgets(buffer, sizeof(buffer), fp);
//    printf("DEBUG buffer is %s\n",buffer);
      sscanf(buffer, "%d: ic %d i %d j %d level %d", &proc[ic], &cellnum[ic], &i[ic], &j[ic], &level[ic]);
      dx[ic]= 1.0/(level[ic]+1);
      x[ic] = (float)i[ic]*dx[ic];
      dy[ic]= 1.0/(level[ic]+1);
      y[ic] = (float)j[ic]*dy[ic];
      data[ic] = (float)proc[ic];
//    printf("DEBUG -- %d: ic %d i %d j %d level %d dx %f x %f dy %f y %f\n",proc[ic],cellnum[ic],i[ic],j[ic],level[ic],dx[ic],x[ic],dy[ic],y[ic]);
      if (x[ic] < xmin) xmin = x[ic];
      if (x[ic]+dx[ic] > xmax) xmax = x[ic]+dx[ic];
      if (y[ic] < ymin) ymin = y[ic];
      if (y[ic]+dy[ic] > ymax) ymax = y[ic]+dy[ic];
   }
   xmin -= 1.0;
   ymin -= 1.0;

   bool do_display_graphics = true;
   set_display_mysize(ncells);
   set_display_cell_coordinates_float(x, dx, y, dy);

   set_display_cell_proc(proc);
   set_display_cell_data_float(data);
   set_display_cellnumber_data(cellnum);
   set_display_indexing_data(i, j, level);
   set_display_autoscale();

   set_display_window(xmin, xmax,
                      ymin, ymax);
   bool outline = true;
   set_display_outline((int)outline);
   int view_mode=1;
   set_display_viewmode(view_mode);

   set_graphics_outline(outline);
   set_graphics_window(xmin, xmax,
                       ymin, ymax);
   set_graphics_mysize(ncells);
   set_graphics_cell_coordinates_float(x, dx, y, dy);
   set_graphics_cell_proc(proc);
   set_graphics_cell_data_float(data);
   set_graphics_viewmode(view_mode);

   init_graphics_output();
   set_graphics_cell_proc(proc);
   //write_graphics_info(0,0,0.0,0,0);

   init_display(&argc, argv, "Shallow Water");
   draw_scene();
   start_main_loop();
}
