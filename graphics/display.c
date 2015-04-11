/*
 *  Copyright (c) 2011, Los Alamos National Security, LLC.
 *  All rights Reserved.
 *
 *  Copyright 2011. Los Alamos National Security, LLC. This software was produced 
 *  under U.S. Government contract DE-AC52-06NA25396 for Los Alamos National 
 *  Laboratory (LANL), which is operated by Los Alamos National Security, LLC 
 *  for the U.S. Department of Energy. The U.S. Government has rights to use, 
 *  reproduce, and distribute this software.  NEITHER THE GOVERNMENT NOR LOS 
 *  ALAMOS NATIONAL SECURITY, LLC MAKES ANY WARRANTY, EXPRESS OR IMPLIED, OR 
 *  ASSUMES ANY LIABILITY FOR THE USE OF THIS SOFTWARE.  If software is modified
 *  to produce derivative works, such modified software should be clearly marked,
 *  so as not to confuse it with the version available from LANL.
 *
 *  Additionally, redistribution and use in source and binary forms, with or without
 *  modification, are permitted provided that the following conditions are met:
 *     * Redistributions of source code must retain the above copyright
 *       notice, this list of conditions and the following disclaimer.
 *     * Redistributions in binary form must reproduce the above copyright
 *       notice, this list of conditions and the following disclaimer in the
 *       documentation and/or other materials provided with the distribution.
 *     * Neither the name of the Los Alamos National Security, LLC, Los Alamos 
 *       National Laboratory, LANL, the U.S. Government, nor the names of its 
 *       contributors may be used to endorse or promote products derived from 
 *       this software without specific prior written permission.
 *  
 *  THIS SOFTWARE IS PROVIDED BY THE LOS ALAMOS NATIONAL SECURITY, LLC AND 
 *  CONTRIBUTORS "AS IS" AND ANY EXPRESS OR IMPLIED WARRANTIES, INCLUDING, BUT 
 *  NOT LIMITED TO, THE IMPLIED WARRANTIES OF MERCHANTABILITY AND FITNESS FOR
 *  A PARTICULAR PURPOSE ARE DISCLAIMED. IN NO EVENT SHALL LOS ALAMOS NATIONAL
 *  SECURITY, LLC OR CONTRIBUTORS BE LIABLE FOR ANY DIRECT, INDIRECT, INCIDENTAL,
 *  SPECIAL, EXEMPLARY, OR CONSEQUENTIAL DAMAGES (INCLUDING, BUT NOT LIMITED TO,
 *  PROCUREMENT OF SUBSTITUTE GOODS OR SERVICES; LOSS OF USE, DATA, OR PROFITS;
 *  OR BUSINESS INTERRUPTION) HOWEVER CAUSED AND ON ANY THEORY OF LIABILITY,
 *  WHETHER IN CONTRACT, STRICT LIABILITY, OR TORT (INCLUDING NEGLIGENCE OR OTHERWISE)
 *  ARISING IN ANY WAY OUT OF THE USE OF THIS SOFTWARE, EVEN IF ADVISED OF THE
 *  POSSIBILITY OF SUCH DAMAGE.
 *  
 *  CLAMR -- LA-CC-11-094
 *  This research code is being developed as part of the 
 *  2011 X Division Summer Workshop for the express purpose
 *  of a collaborative code for development of ideas in
 *  the implementation of AMR codes for Exascale platforms
 *  
 *  AMR implementation of the Wave code previously developed
 *  as a demonstration code for regular grids on Exascale platforms
 *  as part of the Supercomputing Challenge and Los Alamos 
 *  National Laboratory
 *  
 *  Authors: Bob Robey       XCP-2   brobey@lanl.gov
 *           Neal Davis              davis68@lanl.gov, davis68@illinois.edu
 *           David Nicholaeff        dnic@lanl.gov, mtrxknight@aol.com
 *           Dennis Trujillo         dptrujillo@lanl.gov, dptru10@gmail.com
 * 
 */
#include <stdio.h>
#include <stdlib.h>
#include <unistd.h>
#include <math.h>
#include <string.h>
#include "display.h"                                                                                                   
#ifdef HAVE_CONFIG_H
#include "config.h"
#endif

#if defined(MINIMUM_PRECISION)
   typedef float real_t;
#elif defined(MIXED_PRECISION)
   typedef double real_t;
#elif defined(FULL_PRECISION)
   typedef double real_t;
#endif

#define ESCAPE 27
#ifdef HAVE_OPENGL
#define NCOLORS 1000
#else
#ifdef HAVE_MPE
#define NCOLORS 256
#else
#define NCOLORS 1000
#endif
#endif

#define WINSIZE 800

int DrawString(float x, float y, float z, char* string);
void InitGL(int width, int height);
void DrawSquares(void);
void DrawBoxes(void);
void SelectionRec(void);
void mouseClick(int button, int state, int x, int y);
void mouseDrag(int x, int y);
void keyPressed(unsigned char key, int x, int y);
void Scale();
void mpe_main_loop(void);
void display_get_event(void);

struct ColorTable {
   float Red;
   float Green;
   float Blue;
};

static int autoscale = 0;
static double display_circle_radius=-1.0;

#ifndef HAVE_MPI
#if defined(HAVE_MPE)
static int MPI_COMM_WORLD=0;
#endif
#endif

#ifdef HAVE_OPENGL
static struct ColorTable Rainbow[NCOLORS];
static int window;
#endif
#ifdef HAVE_MPE
static MPE_Color *color_array;
static int ncolors=256;
static MPE_XGraph window;
static double xconv = 0.0;
static double yconv = 0.0;
static XFontStruct *font_info;
void (*idlefunction)(void);
#endif

#ifdef HAVE_OPENGL
static real_t xstart, ystart, xend, yend;
#endif
enum mode_choice {EYE, MOVE, DRAW};
static int mode = MOVE;

#if defined(HAVE_MPE)
static int height;
#endif
static int width;
static float display_xmin=0.0, display_xmax=0.0, display_ymin=0.0, display_ymax=0.0;

#ifdef HAVE_OPENGL
static GLfloat xrot = 0.0, yrot = 0.0, xloc = 0.0, zloc = 0.0;
#else
#ifdef HAVE_MPE
static double xrot = 0.0, yrot = 0.0, xloc = 0.0, zloc = 0.0;
#else
static double xrot = 0.0, yrot = 0.0, xloc = 0.0, zloc = 0.0;
#endif
#endif

static int display_outline;
static int display_view_mode = 0;
static int display_mysize    = 0;

enum spatial_data_type {SPATIAL_DOUBLE, SPATIAL_FLOAT};
static int spatial_type = SPATIAL_FLOAT;

static double *x_double=NULL, *y_double=NULL, *dx_double=NULL, *dy_double=NULL;
static float *x_float=NULL, *y_float=NULL, *dx_float=NULL, *dy_float=NULL;

enum plot_data_type {DATA_DOUBLE, DATA_FLOAT};
static int data_type = DATA_FLOAT;
static double *data_double=NULL;
static float *data_float=NULL;
static int *display_proc=NULL;
#if defined(HAVE_MPI) || defined(HAVE_OPENGL)
static int rank = 0;
#endif

int DrawString(float x, float y, float z, char* string) {
#ifdef HAVE_OPENGL
   char *c;
   glColor3f(0.0f, 0.0f, 0.0f);
   glRasterPos3f(x, y, z);
   for(c = string; *c != '\0'; c++) {
      glutBitmapCharacter(GLUT_BITMAP_TIMES_ROMAN_10, *c);
      //glutBitmapCharacter(GLUT_BITMAP_HELVETICA_18, *c);
   }
#endif
#ifdef HAVE_MPE
   int xloc = (int)((x-display_xmin)*xconv);
   int yloc = (int)((display_ymax-y)*yconv);
   //MPE_Draw_string(window, xloc, yloc, MPE_BLACK, string);
#endif
#if ! defined(HAVE_OPENGL) && ! defined(HAVE_MPE)
   // Just to get rid of compiler warnings
   if (1 == 2) printf("DEBUG -- x %f y %f z %f string %s\n",x,y,z,string);
#endif
   return 1;
}

#ifdef HAVE_OPENGL
void InitGL(int width, int height) {
   glClearColor(1.0f, 1.0f, 1.0f, 1.0f);
   glDepthFunc(GL_LESS);
   glShadeModel(GL_SMOOTH);

   glMatrixMode(GL_PROJECTION);
   glLoadIdentity();

   // Just to get rid of compiler warnings
   if (1 == 2) printf("DEBUG -- width %d height %d\n",width,height);
}
#endif
void init_display(int *argc, char **argv, const char *title){

#if ! defined(HAVE_OPENGL) && ! defined(HAVE_MPE)
   // Just to get rid of compiler warnings
   if (1 == 2) printf("DEBUG -- argc %d argv %s title %s\n",*argc,argv[0],title);
#endif

#ifdef HAVE_MPI
   int mpi_init_flag=0;
   MPI_Initialized(&mpi_init_flag);
   if (mpi_init_flag){
      MPI_Comm_rank(MPI_COMM_WORLD, &rank);
   } else {
      MPI_Init(argc, &argv);
   }
#endif

   width = (WINSIZE / (display_ymax - display_ymin)) * (display_xmax - display_xmin);
#ifdef HAVE_OPENGL
   glutInit(argc, argv);
   glutInitDisplayMode(GLUT_RGBA | GLUT_DOUBLE);
   if (rank == 0) {
      glutInitWindowSize(width, WINSIZE);
      glutInitWindowPosition(20, 20);
   } else {
      glutInitWindowSize(1,1);
      glutInitWindowPosition(5, 5);
   }

   window = glutCreateWindow(title);
#endif

#ifdef HAVE_MPE
   char fontname[30];
   char *displayname=NULL;

   /* Open the graphics display */

   if (*argc > 2 && strcmp( argv[1], "-display" ) == 0) {
      displayname = (char *)malloc( strlen( argv[2] ) + 1 );
      strcpy(displayname, argv[2]);
   }

   if (displayname == NULL){
      displayname = getenv("DISPLAY");
      if (displayname == NULL){
         displayname = (char *)malloc((size_t)28);
         displayname = ":0";
      }
   }

   if (MPE_Open_graphics( &window, MPI_COMM_WORLD, displayname,
                         -1, -1, WINSIZE, WINSIZE, 0) != MPE_SUCCESS){
      printf("Problem opening graphics window\n");
      exit(-1);
   }

   xconv = (double)WINSIZE/ (display_xmax-display_xmin);
   yconv = (double)WINSIZE/(display_ymax-display_ymin);

   if (rank == 0){
      XSelectInput( window->xwin->disp, window->xwin->win, MPE_XEVT_IDLE_MASK |
                  ButtonPressMask   |
                  ButtonReleaseMask |
                  Button1MotionMask |
                  KeyPressMask      |
                  ExposureMask      |
                  StructureNotifyMask);
   }

   MPE_Update(window);

   // if we want to change fonts, we can load a new one with the following
   strcpy(fontname,"fixed");
   font_info = XLoadQueryFont( window->xwin->disp, fontname );
#endif

   Scale();

#ifdef HAVE_OPENGL
   glutDisplayFunc(&draw_scene);

   glutKeyboardFunc(&keyPressed);
   glutMouseFunc(&mouseClick);
   glutMotionFunc(&mouseDrag);
   InitGL(width, WINSIZE);
#endif
}

void set_idle_function(void (*function)(void)){
#ifdef HAVE_OPENGL
   glutIdleFunc(function);
#endif
#ifdef HAVE_MPE
   idlefunction = function;
#endif

#if ! defined(HAVE_OPENGL) && ! defined(HAVE_MPE)
   // Just to get rid of compiler warnings
   if (1 == 2) printf("DEBUG -- function %p\n",function);
#endif
}

void start_main_loop(void){
#ifdef HAVE_OPENGL
   glutMainLoop();
#endif
#ifdef HAVE_MPE
   mpe_main_loop();
#endif
}
   
void set_display_window(float display_xmin_in, float display_xmax_in,
                        float display_ymin_in, float display_ymax_in){
   display_xmin = display_xmin_in;
   display_xmax = display_xmax_in;
   display_ymin = display_ymin_in;
   display_ymax = display_ymax_in;
}
void set_display_cell_data_double(double *data_in){
   data_type = DATA_DOUBLE;
   data_double = data_in;
}
void set_display_cell_data_float(float *data_in){
   data_type = DATA_FLOAT;
   data_float = data_in;
}
void set_display_cell_proc(int *display_proc_in){
   display_proc = display_proc_in;
}
void set_display_cell_coordinates_double(double *x_in, double *dx_in, double *y_in, double *dy_in){
   spatial_type = SPATIAL_DOUBLE;
   x_double = x_in;
   dx_double = dx_in;
   y_double = y_in;
   dy_double = dy_in;
}
void set_display_cell_coordinates_float(float *x_in, float *dx_in, float *y_in, float *dy_in){
   spatial_type = SPATIAL_FLOAT;
   x_float = x_in;
   dx_float = dx_in;
   y_float = y_in;
   dy_float = dy_in;
}
void free_display(void){
#ifdef HAVE_OPENGL
   glutDestroyWindow(window);
#endif
}
void DisplayState(void) {
   double scaleMax = 25.0, scaleMin = 0.0;
   int i;
#if defined(HAVE_OPENGL) || defined(HAVE_MPE)
   int color;
#endif
   //vector<real_t> &H = state->H;

   if (autoscale) {
      scaleMax=-1.0e30;
      scaleMin=1.0e30;
      if (data_type == DATA_DOUBLE){
         for(i = 0; i<display_mysize; i++) {
            if (data_double[i] > scaleMax) scaleMax = data_double[i];
            if (data_double[i] < scaleMin) scaleMin = data_double[i];
         }
      } else {
         for(i = 0; i<display_mysize; i++) {
            if (data_float[i] > scaleMax) scaleMax = data_float[i];
            if (data_float[i] < scaleMin) scaleMin = data_float[i];
         }
      }
   }

#if defined(HAVE_OPENGL) || defined(HAVE_MPE)
   double step = NCOLORS/(scaleMax - scaleMin);
#endif
  
#ifdef HAVE_OPENGL
   //set up 2D viewing
   gluOrtho2D(display_xmin, display_xmax, display_ymin, display_ymax);
#endif
  
#ifdef HAVE_OPENGL
   for(i = 0; i < display_mysize; i++) {
      /*Draw filled cells with color set by state value*/
      glPolygonMode(GL_FRONT_AND_BACK, GL_FILL);
      glBegin(GL_QUADS);

      if (data_type == DATA_DOUBLE){
         color = (int)(data_double[i]-scaleMin)*step;
      } else {
         color = (int)(data_float[i]-scaleMin)*step;
      }
      if (color < 0) {
         color=0;
      }
      if (color >= NCOLORS) color = NCOLORS-1;
   
      glColor3f(Rainbow[color].Red, Rainbow[color].Green, Rainbow[color].Blue);
   
      if (spatial_type == SPATIAL_DOUBLE){
         glVertex2f(x_double[i],              y_double[i]);
         glVertex2f(x_double[i]+dx_double[i], y_double[i]);
         glVertex2f(x_double[i]+dx_double[i], y_double[i]+dy_double[i]);
         glVertex2f(x_double[i],              y_double[i]+dy_double[i]);
      } else {
         glVertex2f(x_float[i],             y_float[i]);
         glVertex2f(x_float[i]+dx_float[i], y_float[i]);
         glVertex2f(x_float[i]+dx_float[i], y_float[i]+dy_float[i]);
         glVertex2f(x_float[i],             y_float[i]+dy_float[i]);
      }
      glEnd();
   
      /*Draw cells again as outline*/
      if (display_outline) {
         glColor3f(0.0f, 0.0f, 0.0f);
         glPolygonMode(GL_FRONT_AND_BACK, GL_LINE);
         glBegin(GL_QUADS);
         if (spatial_type == SPATIAL_DOUBLE){
            glVertex2f(x_double[i],              y_double[i]);
            glVertex2f(x_double[i]+dx_double[i], y_double[i]);
            glVertex2f(x_double[i]+dx_double[i], y_double[i]+dy_double[i]);
            glVertex2f(x_double[i],              y_double[i]+dy_double[i]);
         } else {
            glVertex2f(x_float[i],             y_float[i]);
            glVertex2f(x_float[i]+dx_float[i], y_float[i]);
            glVertex2f(x_float[i]+dx_float[i], y_float[i]+dy_float[i]);
            glVertex2f(x_float[i],             y_float[i]+dy_float[i]);
         }
         glEnd();
      }
   }
#endif
#ifdef HAVE_MPE
   int xloc, xwid, yloc, ywid;
   int xloc1, xloc2, yloc1, yloc2;
   if (spatial_type == SPATIAL_DOUBLE){
      for(i = 0; i < display_mysize; i++) {
         if (data_type == DATA_DOUBLE){
            color = (int)(data_double[i]-scaleMin)*step;
         } else {
            color = (int)(data_float[i]-scaleMin)*step;
         }
         color = NCOLORS-color;
         if (color < 0) {
            color=0;
         }
         if (color >= NCOLORS) color = NCOLORS-1;

         //printf("DEBUG mesh -- xconv is %lf xmin %lf x[%d]=%lf\n",xconv,display_xmin,i,x[i]);
         //printf("DEBUG mesh -- yconv is %lf ymin %lf ymax %lf y[%d]=%lf\n",yconv,display_ymin,display_ymax,i,y[i]);
         xloc = (int)((x_double[i]-display_xmin)*xconv);
         xwid = (int)((x_double[i]+dx_double[i]-display_xmin)*xconv-xloc);
         yloc = (int)((display_ymax-(y_double[i]+dy_double[i]))*yconv);
         ywid = (int)((display_ymax-y_double[i])*yconv);
         ywid -= yloc;
         //printf("xloc is %d xwid %d yloc is %d ywid %d color is  %d step is %d\n",xloc,xwid,yloc,ywid,color,step);
         MPE_Fill_rectangle(window, xloc, yloc, xwid, ywid, color_array[color]);

         if (display_outline) {
            xloc1 = (int)((x_double[i]-display_xmin)*xconv);
            xloc2 = (int)((x_double[i]+dx_double[i]-display_xmin)*xconv);
            yloc1 = (int)((display_ymax-y_double[i])*yconv);
            yloc2 = (int)((display_ymax-(y_double[i]+dy_double[i]))*yconv);
            //printf("xloc1 is %d xloc2 %d yloc1 is %d yloc2 %d\n",xloc1,xloc2,yloc1,yloc2);
            MPE_Draw_line(window, xloc1, yloc2, xloc2, yloc2, MPE_BLACK);
            MPE_Draw_line(window, xloc1, yloc1, xloc2, yloc1, MPE_BLACK);
            MPE_Draw_line(window, xloc1, yloc1, xloc1, yloc2, MPE_BLACK);
            MPE_Draw_line(window, xloc2, yloc1, xloc2, yloc2, MPE_BLACK);
         }
      }
   } else {
      for(i = 0; i < display_mysize; i++) {
         if (data_type == DATA_DOUBLE){
            color = (int)(data_double[i]-scaleMin)*step;
         } else {
            color = (int)(data_float[i]-scaleMin)*step;
         }
         color = NCOLORS-color;
         if (color < 0) {
            color=0;
         }
         if (color >= NCOLORS) color = NCOLORS-1;

         //printf("DEBUG mesh -- xconv is %lf xmin %lf x[%d]=%lf\n",xconv,display_xmin,i,x[i]);
         //printf("DEBUG mesh -- yconv is %lf ymin %lf ymax %lf y[%d]=%lf\n",yconv,display_ymin,display_ymax,i,y[i]);
         xloc = (int)((x_float[i]-display_xmin)*xconv);
         xwid = (int)((x_float[i]+dx_float[i]-display_xmin)*xconv-xloc);
         yloc = (int)((display_ymax-(y_float[i]+dy_float[i]))*yconv);
         ywid = (int)((display_ymax-y_float[i])*yconv);
         ywid -= yloc;
         //printf("xloc is %d xwid %d yloc is %d ywid %d color is  %d step is %d\n",xloc,xwid,yloc,ywid,color,step);
         MPE_Fill_rectangle(window, xloc, yloc, xwid, ywid, color_array[color]);

         if (display_outline) {
            xloc1 = (int)((x_float[i]-display_xmin)*xconv);
            xloc2 = (int)((x_float[i]+dx_float[i]-display_xmin)*xconv);
            yloc1 = (int)((display_ymax-y_float[i])*yconv);
            yloc2 = (int)((display_ymax-(y_float[i]+dy_float[i]))*yconv);
            //printf("xloc1 is %d xloc2 %d yloc1 is %d yloc2 %d\n",xloc1,xloc2,yloc1,yloc2);
            MPE_Draw_line(window, xloc1, yloc2, xloc2, yloc2, MPE_BLACK);
            MPE_Draw_line(window, xloc1, yloc1, xloc2, yloc1, MPE_BLACK);
            MPE_Draw_line(window, xloc1, yloc1, xloc1, yloc2, MPE_BLACK);
            MPE_Draw_line(window, xloc2, yloc1, xloc2, yloc2, MPE_BLACK);
         }
      }
   }
   MPE_Update(window);
#endif

  
#ifdef HAVE_OPENGL
   glColor3f(0.0f, 0.0f, 0.0f);

   /* Draw circle for initial conditions */
   if (display_circle_radius > 0.0) {
      double radians;
      glBegin(GL_LINE_LOOP);
      for (double angle = 0.0; angle <= 360.0; angle+=1.0) {
         radians=angle*3.14159/180.0;
         glVertex2f(sin(radians) * display_circle_radius,cos(radians) * display_circle_radius);
      }
      glEnd();
   }
#endif
}

void DrawSquares(void) {
   if (display_proc == NULL) return;
#if defined(HAVE_OPENGL) || defined(HAVE_MPE)
   int i, color;
   int step = NCOLORS/(display_proc[display_mysize-1]+1);
#endif
   //step utilizes whole range of color, assumes last element of proc is last processor

#ifdef HAVE_OPENGL
   gluOrtho2D(display_xmin, display_xmax, display_ymin, display_ymax);
   //set up 2D viewing

   for(i = 0; i < display_mysize; i++) {
      /*Draw filled cells with color set by processor it is assigned to*/
      glPolygonMode(GL_FRONT_AND_BACK, GL_FILL);
      glBegin(GL_QUADS);
         color = display_proc[i]*step;
         glColor3f(Rainbow[color].Red, Rainbow[color].Green, Rainbow[color].Blue);

         if (spatial_type == SPATIAL_DOUBLE){
            glVertex2f(x_double[i],              y_double[i]);
            glVertex2f(x_double[i]+dx_double[i], y_double[i]);
            glVertex2f(x_double[i]+dx_double[i], y_double[i]+dy_double[i]);
            glVertex2f(x_double[i],              y_double[i]+dy_double[i]);
         } else {
            glVertex2f(x_float[i],             y_float[i]);
            glVertex2f(x_float[i]+dx_float[i], y_float[i]);
            glVertex2f(x_float[i]+dx_float[i], y_float[i]+dy_float[i]);
            glVertex2f(x_float[i],             y_float[i]+dy_float[i]);
         }
      glEnd();

      /*Draw cells again as outline*/
      glColor3f(0.0f, 0.0f, 0.0f);
      glPolygonMode(GL_FRONT_AND_BACK, GL_LINE);
      glBegin(GL_QUADS);
         if (spatial_type == SPATIAL_DOUBLE) {
            glVertex2f(x_double[i],              y_double[i]);
            glVertex2f(x_double[i]+dx_double[i], y_double[i]);
            glVertex2f(x_double[i]+dx_double[i], y_double[i]+dy_double[i]);
            glVertex2f(x_double[i],              y_double[i]+dy_double[i]);
         } else {
            glVertex2f(x_float[i],             y_float[i]);
            glVertex2f(x_float[i]+dx_float[i], y_float[i]);
            glVertex2f(x_float[i]+dx_float[i], y_float[i]+dy_float[i]);
            glVertex2f(x_float[i],             y_float[i]+dy_float[i]);
         }
      glEnd();
   }
#endif
#ifdef HAVE_MPE
   int xloc, xwid, yloc, ywid;
   int xloc1, xloc2, yloc1, yloc2;
   if (spatial_type == SPATIAL_DOUBLE){
      for(i = 0; i < display_mysize; i++) {
         //printf("DEBUG mesh -- xconv is %lf xmin %lf x[%d]=%lf\n",xconv,display_xmin,i,x[i]);
         xloc = (int)((x_double[i]-display_xmin)*xconv);
         xwid = (int)((x_double[i]+dx_double[i]-display_xmin)*xconv-xloc);
         yloc = (int)((display_ymax-(y_double[i]+dy_double[i]))*yconv);
         ywid = (int)((display_ymax-y_double[i])*yconv);
         ywid -= yloc;
         color = display_proc[i]*step;
         //printf("xloc is %d xwid %d yloc is %d ywid %d color is  %d step is %d\n",xloc,xwid,yloc,ywid,color,step);
         MPE_Fill_rectangle(window, xloc, yloc, xwid, ywid, color_array[color]);

         xloc1 = (int)((x_double[i]-display_xmin)*xconv);
         xloc2 = (int)((x_double[i]+dx_double[i]-display_xmin)*xconv);
         yloc1 = (int)((display_ymax-y_double[i])*yconv);
         yloc2 = (int)((display_ymax-(y_double[i]+dy_double[i]))*yconv);
         //printf("xloc1 is %d xloc2 %d yloc1 is %d yloc2 %d\n",xloc1,xloc2,yloc1,yloc2);
         MPE_Draw_line(window, xloc1, yloc2, xloc2, yloc2, MPE_BLACK);
         MPE_Draw_line(window, xloc1, yloc1, xloc2, yloc1, MPE_BLACK);
         MPE_Draw_line(window, xloc1, yloc1, xloc1, yloc2, MPE_BLACK);
         MPE_Draw_line(window, xloc2, yloc1, xloc2, yloc2, MPE_BLACK);

         //MPE_Fill_rectangle(window, xloc, yloc, xwid, ywid, MPE_RED);
      }
   } else {
      for(i = 0; i < display_mysize; i++) {
         //printf("DEBUG mesh -- xconv is %lf xmin %lf x[%d]=%lf\n",xconv,display_xmin,i,x[i]);
         xloc = (int)((x_float[i]-display_xmin)*xconv);
         xwid = (int)((x_float[i]+dx_float[i]-display_xmin)*xconv-xloc);
         yloc = (int)((display_ymax-(y_float[i]+dy_float[i]))*yconv);
         ywid = (int)((display_ymax-y_float[i])*yconv);
         ywid -= yloc;
         color = display_proc[i]*step;
         //printf("xloc is %d xwid %d yloc is %d ywid %d color is  %d step is %d\n",xloc,xwid,yloc,ywid,color,step);
         MPE_Fill_rectangle(window, xloc, yloc, xwid, ywid, color_array[color]);

         xloc1 = (int)((x_float[i]-display_xmin)*xconv);
         xloc2 = (int)((x_float[i]+dx_float[i]-display_xmin)*xconv);
         yloc1 = (int)((display_ymax-y_float[i])*yconv);
         yloc2 = (int)((display_ymax-(y_float[i]+dy_float[i]))*yconv);
         //printf("xloc1 is %d xloc2 %d yloc1 is %d yloc2 %d\n",xloc1,xloc2,yloc1,yloc2);
         MPE_Draw_line(window, xloc1, yloc2, xloc2, yloc2, MPE_BLACK);
         MPE_Draw_line(window, xloc1, yloc1, xloc2, yloc1, MPE_BLACK);
         MPE_Draw_line(window, xloc1, yloc1, xloc1, yloc2, MPE_BLACK);
         MPE_Draw_line(window, xloc2, yloc1, xloc2, yloc2, MPE_BLACK);

         //MPE_Fill_rectangle(window, xloc, yloc, xwid, ywid, MPE_RED);
      }
   }
#endif

#ifdef HAVE_OPENGL
   /* Draw circle for initial conditions */
   if (display_circle_radius > 0.0) {
      double radians;
      glBegin(GL_LINE_LOOP);
      for (double angle = 0.0; angle <= 360.0; angle+=1.0) {
         radians=angle*3.14159/180.0;
         glVertex2f(sin(radians) * display_circle_radius,cos(radians) * display_circle_radius);
      }
      glEnd();
   }
#endif
  
#ifdef HAVE_OPENGL
   /*Trace order of cells with line going from center to center*/
   glBegin(GL_LINE_STRIP);
      if (spatial_type == SPATIAL_DOUBLE) {
         for(i = 0; i < display_mysize; i++) {
            glVertex2f(x_double[i]+dx_double[i]/2, y_double[i]+dy_double[i]/2);
         }
      } else {
         for(i = 0; i < display_mysize; i++) {
            glVertex2f(x_float[i]+dx_float[i]/2, y_float[i]+dy_float[i]/2);
         }
      }
   glEnd();
#endif
#ifdef HAVE_MPE
   int xloc_old, yloc_old;
   if (spatial_type == SPATIAL_DOUBLE){
      xloc_old = (int)((x_double[0]+0.5*dx_double[0]-display_xmin)*xconv);
      yloc_old = (int)((y_double[0]+0.5*dy_double[0]-display_ymin)*yconv);
      for(i = 1; i < display_mysize; i++) {
         xloc = (int)((x_double[i]+0.5*dx_double[i]-display_xmin)*xconv);
         yloc = (int)((y_double[i]+0.5*dy_double[i]-display_ymin)*yconv);
         MPE_Draw_line(window,xloc_old,yloc_old,xloc,yloc,MPE_BLACK);
         xloc_old = xloc;
         yloc_old = yloc;
      }
   } else {
      xloc_old = (int)((x_float[0]+0.5*dx_float[0]-display_xmin)*xconv);
      yloc_old = (int)((y_float[0]+0.5*dy_float[0]-display_ymin)*yconv);
      for(i = 1; i < display_mysize; i++) {
         xloc = (int)((x_float[i]+0.5*dx_float[i]-display_xmin)*xconv);
         yloc = (int)((y_float[i]+0.5*dy_float[i]-display_ymin)*yconv);
         MPE_Draw_line(window,xloc_old,yloc_old,xloc,yloc,MPE_BLACK);
         xloc_old = xloc;
         yloc_old = yloc;
      }
   }
#ifdef HAVE_MPI
   MPI_Barrier(MPI_COMM_WORLD);
#endif
   MPE_Update(window);
#endif
}

void DrawBoxes(void) {

#ifdef HAVE_OPENGL
   int i, color, step = NCOLORS/(display_proc[display_mysize-1]+1);
   //step utilizes whole range of color, assumes last element of proc is last processor

   glFrustum(display_xmin-1, display_xmax, display_ymin-1, display_ymax, 5.0, 10.0);
   glTranslatef(0.0f, 0.0f, -6.0f); //must move into screen to draw

   if (spatial_type == SPATIAL_DOUBLE){
      for(i = 0; i < display_mysize; i++) {
         /*Draw filled cells with color set by processor*/
         glPolygonMode(GL_FRONT_AND_BACK, GL_FILL);
         glBegin(GL_QUADS);
            color = display_proc[i]*step;
            glColor3f(Rainbow[color].Red, Rainbow[color].Green, Rainbow[color].Blue);

            /*Front Face*/
            glVertex3f(x_double[i],              y_double[i],              0.0f);
            glVertex3f(x_double[i]+dx_double[i], y_double[i],              0.0f);
            glVertex3f(x_double[i]+dx_double[i], y_double[i]+dy_double[i], 0.0f);
            glVertex3f(x_double[i],              y_double[i]+dy_double[i], 0.0f);
            /*Right Face*/
            glVertex3f(x_double[i]+dx_double[i], y_double[i],              0.0f);
            glVertex3f(x_double[i]+dx_double[i], y_double[i]+dy_double[i], 0.0f);
            glVertex3f(x_double[i]+dx_double[i], y_double[i]+dy_double[i],-1.0f);
            glVertex3f(x_double[i]+dx_double[i], y_double[i],             -1.0f);
            /*Back Face*/
            glVertex3f(x_double[i]+dx_double[i], y_double[i],             -1.0f);
            glVertex3f(x_double[i]+dx_double[i], y_double[i]+dy_double[i],-1.0f);
            glVertex3f(x_double[i],              y_double[i]+dy_double[i],-1.0f);
            glVertex3f(x_double[i],              y_double[i],             -1.0f);
            /*Left Face*/
            glVertex3f(x_double[i],              y_double[i],              0.0f);
            glVertex3f(x_double[i],              y_double[i],             -1.0f);
            glVertex3f(x_double[i],              y_double[i]+dy_double[i],-1.0f);
            glVertex3f(x_double[i],              y_double[i]+dy_double[i], 0.0f);
            /*Bottom*/
            glVertex3f(x_double[i],              y_double[i],              0.0f);
            glVertex3f(x_double[i],              y_double[i],             -1.0f);
            glVertex3f(x_double[i]+dx_double[i], y_double[i],             -1.0f);
            glVertex3f(x_double[i]+dx_double[i], y_double[i],              0.0f);
            /*Top*/
            glVertex3f(x_double[i],              y_double[i]+dy_double[i], 0.0f);
            glVertex3f(x_double[i],              y_double[i]+dy_double[i],-1.0f);
            glVertex3f(x_double[i]+dx_double[i], y_double[i]+dy_double[i],-1.0f);
            glVertex3f(x_double[i]+dx_double[i], y_double[i]+dy_double[i], 0.0f);
         glEnd();
      }
   } else {
      for(i = 0; i < display_mysize; i++) {
         /*Draw filled cells with color set by processor*/
         glPolygonMode(GL_FRONT_AND_BACK, GL_FILL);
         glBegin(GL_QUADS);
            color = display_proc[i]*step;
            glColor3f(Rainbow[color].Red, Rainbow[color].Green, Rainbow[color].Blue);

            /*Front Face*/
            glVertex3f(x_float[i],             y_float[i],             0.0f);
            glVertex3f(x_float[i]+dx_float[i], y_float[i],             0.0f);
            glVertex3f(x_float[i]+dx_float[i], y_float[i]+dy_float[i], 0.0f);
            glVertex3f(x_float[i],             y_float[i]+dy_float[i], 0.0f);
            /*Right Face*/
            glVertex3f(x_float[i]+dx_float[i], y_float[i],             0.0f);
            glVertex3f(x_float[i]+dx_float[i], y_float[i]+dy_float[i], 0.0f);
            glVertex3f(x_float[i]+dx_float[i], y_float[i]+dy_float[i],-1.0f);
            glVertex3f(x_float[i]+dx_float[i], y_float[i],            -1.0f);
            /*Back Face*/
            glVertex3f(x_float[i]+dx_float[i], y_float[i],            -1.0f);
            glVertex3f(x_float[i]+dx_float[i], y_float[i]+dy_float[i],-1.0f);
            glVertex3f(x_float[i],             y_float[i]+dy_float[i],-1.0f);
            glVertex3f(x_float[i],             y_float[i],            -1.0f);
            /*Left Face*/
            glVertex3f(x_float[i],             y_float[i],             0.0f);
            glVertex3f(x_float[i],             y_float[i],            -1.0f);
            glVertex3f(x_float[i],             y_float[i]+dy_float[i],-1.0f);
            glVertex3f(x_float[i],             y_float[i]+dy_float[i], 0.0f);
            /*Bottom*/
            glVertex3f(x_float[i],             y_float[i],             0.0f);
            glVertex3f(x_float[i],             y_float[i],            -1.0f);
            glVertex3f(x_float[i]+dx_float[i], y_float[i],            -1.0f);
            glVertex3f(x_float[i]+dx_float[i], y_float[i],             0.0f);
            /*Top*/
            glVertex3f(x_float[i],             y_float[i]+dy_float[i], 0.0f);
            glVertex3f(x_float[i],             y_float[i]+dy_float[i],-1.0f);
            glVertex3f(x_float[i]+dx_float[i], y_float[i]+dy_float[i],-1.0f);
            glVertex3f(x_float[i]+dx_float[i], y_float[i]+dy_float[i], 0.0f);
         glEnd();
      }
   }
#endif
}



void set_display_viewmode(int display_view_mode_in){
   display_view_mode = display_view_mode_in;
}
void set_display_mysize(int display_mysize_in){
   display_mysize = display_mysize_in;
}
void set_display_circle_radius(double display_circle_radius_in){
   display_circle_radius = display_circle_radius_in;
}
void set_display_outline(int display_outline_in){
   display_outline = display_outline_in;
}

void draw_scene(void) {
#ifdef HAVE_OPENGL
   if (rank) return;

   glClear(GL_COLOR_BUFFER_BIT | GL_DEPTH_BUFFER_BIT);
   glLoadIdentity();
#endif
#ifdef HAVE_MPE
   if (rank == 0) {
      MPE_Fill_rectangle(window, 0, 0, WINSIZE, WINSIZE, MPE_WHITE);
   }
   MPE_Update(window);
#ifdef HAVE_MPI
   MPI_Barrier(MPI_COMM_WORLD);
#endif
#endif


   if (display_view_mode == 0) {
      DrawSquares();
   } else {
      DisplayState();
   }

   if (display_mysize <=500) {
      char c[10];
      if (data_type == DATA_DOUBLE){
         for(int i = 0; i < display_mysize; i++) {
            sprintf(c, "%d", i);
            DrawString(x_double[i]+0.5*dx_double[i], y_double[i]+0.5*dy_double[i], 0.0, c);
         }
      } else {
         for(int i = 0; i < display_mysize; i++) {
            sprintf(c, "%d", i);
            DrawString(x_float[i]+0.5*dx_float[i], y_float[i]+0.5*dy_float[i], 0.0, c);
         }
      }
   }

#ifdef HAVE_OPENGL
   if(mode == DRAW) {
      SelectionRec();
   }

   glLoadIdentity();

   glutSwapBuffers();
#endif

#ifdef HAVE_MPE
#ifdef HAVE_MPI
   MPI_Barrier(MPI_COMM_WORLD);
#endif
   MPE_Update(window);
#endif

   //if (display_mysize <= 500) sleep(1);
}


#ifdef HAVE_OPENGL
void SelectionRec(void) {
   glColor3f(0.0f, 0.0f, 0.0f);
   glLineWidth(2.0f);
   glBegin(GL_QUADS);
      glVertex2f(xstart, ystart);
      glVertex2f(xstart, yend);
      glVertex2f(xend, yend);
      glVertex2f(xend, ystart);
   glEnd();
   glLineWidth(1.0f);
}
#endif

void mouseClick(int button, int state, int x, int y) {
#ifdef HAVE_OPENGL
   if(state == GLUT_DOWN) {
      mode = EYE;
      xstart = (display_xmax-display_ymin)*(x/width)+display_ymin;
      ystart = (display_ymax-display_ymin)*((float)(WINSIZE-y)/WINSIZE)+display_ymin;
   }
   if(state == GLUT_UP) {
      glutPostRedisplay();
      mode = DRAW;
      xend = (display_xmax-display_ymin)*(x/width)+display_ymin;
      yend = (display_xmax-display_ymin)*((float)(WINSIZE-y)/WINSIZE)+display_ymin;
   }
   // Just to get rid of compiler warnings
   if (1 == 2) printf("DEBUG -- button %d\n",button);
#else
   // Just to get rid of compiler warnings
   if (1 == 2) printf("DEBUG -- x %d y %d button %d state %d\n",x,y,button,state);
#endif
}
void mouseDrag(int x, int y) {
#ifdef HAVE_OPENGL
   glutPostRedisplay();
   mode = DRAW;
   xend = (display_xmax-display_ymin)*(x/width)+display_ymin;
   yend = (display_xmax-display_ymin)*((float)(WINSIZE-y)/WINSIZE)+display_ymin;
#else
   // Just to get rid of compiler warnings
   if (1 == 2) printf("DEBUG -- x %d y %d\n",x,y);
#endif
}

void keyPressed(unsigned char key, int x, int y) {

    //usleep(100);

    if(key == ESCAPE) {
       free_display();
       exit(0);
    }
    if(key == 'm')   { mode = MOVE; }
    if(key == 'e')   { mode = EYE; }
    if(mode == EYE) {
      if(key == 'o') { xrot -= 5.0; }
      if(key == 'l') { xrot += 5.0; }
      if(key == 'k') { yrot -= 5.0; }
      if(key == ';') { yrot += 5.0; }
    }
    if(mode == MOVE) {
      if(key == 'o') { zloc += 1.0; }
      if(key == 'l') { zloc -= 1.0; }
      if(key == 'k') { xloc += 1.0; }
      if(key == ';') { xloc -= 1.0; }
    }

   // Just to get rid of compiler warnings
   if (1 == 2) printf("DEBUG -- x %d y %d\n",x,y);
}
void Scale() {
#ifdef HAVE_OPENGL
   int i;
   double r;
   for (i=0, r=0.0; i<200; i++, r+=.005) {
         Rainbow[i].Red = 0.0;
         Rainbow[i].Green = r;
         Rainbow[i].Blue = 1.0;
   }
   for (i=0, r=1.0; i<200; i++, r-=.005) {
         Rainbow[200+i].Red = 0.0;
         Rainbow[200+i].Green = 1.0;
         Rainbow[200+i].Blue = r;
   }
   for (i=0, r=0.0; i<200; i++, r+=.005) {
         Rainbow[400+i].Red = r;
         Rainbow[400+i].Green = 1.0;
         Rainbow[400+i].Blue = 0.0;
   }
   for (i=0, r=1.0; i<200; i++, r-=.005) {
         Rainbow[600+i].Red = 1.0;
         Rainbow[600+i].Green = r;
         Rainbow[600+i].Blue = 0.0;
   }
   for (i=0, r=0.0; i<200; i++, r+=.005) {
         Rainbow[800+i].Red = 1.0;
         Rainbow[800+i].Green = 0.0;
         Rainbow[800+i].Blue = r;
   }
#endif

#ifdef HAVE_MPE
   color_array = (MPE_Color *) malloc(sizeof(MPE_Color)*ncolors);

   int ierr = MPE_Make_color_array(window, ncolors, color_array);
   if (ierr && rank == 0) printf("Error(Make_color_array): ierr is %d\n",ierr);

#endif
}

void mpe_main_loop(void)
{
#ifdef HAVE_MPE
   while (1) {
      display_get_event();
      idlefunction();
   }
#endif
}

/********************************************************************************/
void display_get_event(void)
/********************************************************************************/
{

#ifdef HAVE_MPE
   //double xmid, ymid;
   //int xrel, yrel;
   XEvent event;
   long EventMask;
   int EventFlag;
   char keys[20];
   int numChar;
   KeySym keysym;
   XComposeStatus compose;
   int rank=0;
   unsigned char key;
   int button, xloc, yloc, special_event;
   double xcor, ycor;


#ifdef HAVE_MPI
   MPI_Comm_rank(MPI_COMM_WORLD, &rank);
#endif

   if (window->Cookie != MPE_G_COOKIE) {
      printf("Handle argument is incorrect or corrupted\n" );
   }

   key = '\0';
   button = -1;
   xcor = 0.0;
   ycor = 0.0;
   special_event = 0;

   EventMask= ButtonPressMask | ButtonReleaseMask | KeyPressMask | ExposureMask | StructureNotifyMask ;


   EventFlag=0;
   //if (rank ==0){
      if(XCheckWindowEvent(window->xwin->disp,window->xwin->win,EventMask,&event)){
         EventFlag=1;
      }
   //}
#ifdef HAVE_MPI
   MPI_Bcast(&EventFlag, 1, MPI_INT, 0, MPI_COMM_WORLD);
#endif
     
   if (EventFlag){
#ifdef HAVE_MPI
      MPI_Bcast(&event, (int)sizeof(XEvent), MPI_BYTE, 0, MPI_COMM_WORLD);
#endif

      //printf("Event type is %d\n",event.type);
      int state;
      switch(event.type){
         case ButtonPress:
            //printf("Button Press is %d\n",event.xbutton.button);
            button=event.xbutton.button;
            xloc=event.xbutton.x;
            yloc=event.xbutton.y;
            xcor=display_xmin+(double)xloc/xconv;
            ycor=display_ymin+(double)(height-yloc -1)/yconv;

            state = 1; // Button Down
            mouseClick((int)event.xbutton.button, state, xloc, yloc);

            //printf("DEBUG graphics -- button is %d loc is %d %d cor is %lf %lf\n",
            //                       button,    xloc,yloc,xcor,ycor);
            break;
         case KeyPress:
            if (rank == 0) {
               numChar=XLookupString((XKeyEvent *) &event, keys, 20, &keysym, &compose);
               printf("%c:%c:%d:\n",keys[0],keys[1],(int)keysym);
               if(((keysym>=XK_KP_Space)&&(keysym<=XK_KP_9))
                  ||((keysym>XK_space)&&(keysym<XK_asciitilde))){
                  key=keys[0];
               }
            }
#ifdef HAVE_MPI
            MPI_Bcast(&key,           1, MPI_CHAR,   0, MPI_COMM_WORLD);
#endif
     
            xloc=event.xkey.x;
            yloc=event.xkey.y;
     
            xcor=display_xmin+(double)xloc/xconv;
            ycor=display_ymin+(double)(height-yloc-1)/yconv;

            keyPressed(key, xloc, yloc);
     
            //printf("DEBUG graphics -- key is '%c' loc is %d %d cor is %lf %lf\n",
            //                         key,      xloc,yloc,xcor,ycor);
            break;
         case Expose:
            if (event.xexpose.count == 0){
               special_event=1;
            }
            //printf("DEBUG graphics -- special event is \n",*special_event);
            break;
         case ConfigureNotify:
            /* Window has been resized */
/*
            width=event.xconfigure.width; //-SCALEWIDTH;
            double ratio = (double)event.xconfigure.height/(double)height;
            height *= ratio;
     
            conHeight *= ratio;
            conWidth  = width - SCALEWIDTH;
     
            graHeight = height - conHeight - 2*TEXTHEIGHT;
            if (graHeight < 0) graHeight  = 0;
            graWidth  = width - SCALEWIDTH;
     
            xconv = (double)conWidth/ (display_xmax-xmin);
            yconv = (double)conHeight/(display_ymax-display_ymin);
*/

            //printf("DEBUG graphics -- window has been resized to %d %d\n",width, height);
            break;
      } // switch
     
      //printf("DEBUG -- button is %d key is '%c' loc is %d %d cor is %lf %lf\n",
      //                button,     key,      xloc,yloc,xcor,ycor);
   }
     
   //MPE_Update(window);
   //usleep(300000);
     
#endif

   return;
}
