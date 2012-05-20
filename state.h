/*
 *  Copyright (c) 2011-2012, Los Alamos National Security, LLC.
 *  All rights Reserved.
 *
 *  Copyright 2011-2012. Los Alamos National Security, LLC. This software was produced 
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
#ifndef STATE_H_
#define STATE_H_

#include "mesh.h"
#include "ezcl/ezcl.h"
#include "l7/l7.h"

extern "C" void do_calc(void);

enum SIGN_RULE {
   DIAG_RULE,
   X_RULE,
   Y_RULE,
};

using namespace std;

class State {
   
public:
   vector<real> H;
   vector<real> U;
   vector<real> V;

   cl_mem dev_H;
   cl_mem dev_U;
   cl_mem dev_V;

   cl_mem dev_mass_sum;
   cl_mem dev_deltaT;

   cl_event apply_BCs_event;

   double   cpu_time_apply_BCs,
            cpu_time_set_timestep,
            cpu_time_finite_difference,
            cpu_time_refine_potential,
            cpu_time_rezone_all,
            cpu_time_mass_sum;

   long     gpu_time_apply_BCs,
            gpu_time_set_timestep,
            gpu_time_finite_difference,
            gpu_time_refine_potential,
            gpu_time_rezone_all,
            gpu_time_mass_sum,
            gpu_time_read,
            gpu_time_write;

   // constructor -- allocates state arrays to size ncells
   State(int ncells, cl_context context);

   State(const State&); // To block copy constructor so copies are not made inadvertently

   void init(int ncells, cl_context context, int do_gpu_calc);

   double get_cpu_time_apply_BCs(void)         {return(cpu_time_apply_BCs);};
   double get_cpu_time_set_timestep(void)      {return(cpu_time_set_timestep);};
   double get_cpu_time_finite_difference(void) {return(cpu_time_finite_difference);};
   double get_cpu_time_refine_potential(void)  {return(cpu_time_refine_potential);};
   double get_cpu_time_rezone_all(void)        {return(cpu_time_rezone_all);};
   double get_cpu_time_mass_sum(void)          {return(cpu_time_mass_sum);};

   long get_gpu_time_apply_BCs(void)         {return(gpu_time_apply_BCs);};
   long get_gpu_time_set_timestep(void)      {return(gpu_time_set_timestep);};
   long get_gpu_time_finite_difference(void) {return(gpu_time_finite_difference);};
   long get_gpu_time_refine_potential(void)  {return(gpu_time_refine_potential);};
   long get_gpu_time_rezone_all(void)        {return(gpu_time_rezone_all);};
   long get_gpu_time_mass_sum(void)          {return(gpu_time_mass_sum);};
   long get_gpu_time_read(void)              {return(gpu_time_read);};
   long get_gpu_time_write(void)             {return(gpu_time_write);};

   void allocate_device_memory(size_t ncells);
   void resize_old_device_memory(size_t ncells);

   void add_boundary_cells(Mesh *mesh);
   void apply_boundary_conditions(Mesh *mesh);
   void remove_boundary_cells(Mesh *mesh);
   double set_timestep(Mesh *mesh, double g, double sigma);
   double gpu_set_timestep(cl_command_queue command_queue, Mesh *mesh, double sigma);
   
   void fill_circle(Mesh *mesh, double circ_radius, double fill_value, double background);
   void state_reorder(vector<int> iorder);
   void rezone_all(Mesh *mesh, vector<int> mpot, int add_ncells);
   void gpu_rezone_all(cl_command_queue command_queue, Mesh *mesh, size_t &ncells, size_t new_ncells, size_t old_ncells, bool localStencil, 
      cl_mem dev_mpot, cl_mem dev_ioffset);
   void gpu_rezone_all_local(cl_command_queue command_queue, Mesh *mesh, size_t &ncells, size_t new_ncells, size_t old_ncells, bool localStencil, 
      cl_mem dev_mpot, cl_mem dev_ioffset);
   void calc_refine_potential(Mesh *mesh, vector<int> &mpot, int &icount, int &jcount);
   void gpu_calc_refine_potential(cl_command_queue command_queue, Mesh *mesh, cl_mem dev_mpot, cl_mem dev_result, cl_mem dev_ioffset);
   void gpu_calc_refine_potential_local(cl_command_queue command_queue, Mesh *mesh, cl_mem dev_mpot, cl_mem dev_result, cl_mem dev_ioffset);
   
   void calc_finite_difference(Mesh *mesh, double deltaT);
   void calc_finite_difference_local(Mesh *mesh, double deltaT);
   void gpu_calc_finite_difference(cl_command_queue, Mesh *mesh, double deltaT);
   void gpu_calc_finite_difference_local(cl_command_queue, Mesh *mesh, double deltaT);

   void symmetry_check(Mesh *mesh, const char *string, vector<int> sym_index, double eps, 
                       SIGN_RULE sign_rule, int &flag);
   double mass_sum(Mesh *mesh, bool enhanced_precision_sum);
   double gpu_mass_sum(cl_command_queue command_queue, Mesh *mesh, bool enhanced_precision_sum);
   double gpu_mass_sum_local(cl_command_queue command_queue, Mesh *mesh, bool enhanced_precision_sum);

   void output_timing_info(Mesh *mesh, int do_cpu_calc, int do_gpu_calc, double elapsed_time);
   void compare_state_gpu_global_to_cpu_global(cl_command_queue command_queue, const char* string, int cycle, uint ncells);
   void compare_state_cpu_local_to_cpu_global(State *state_global, const char* string, int cycle, uint ncells, uint ncells_global, int *nsizes, int *ndispl);
private:
   void parallel_timer_output(int numpe, int mype, const char *string, double local_time);

};

#endif // ifndef STATE_H_

