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
 *  
 *  This file and the associated header is based on a file from the capablanca
 *  project available under the MIT open-source license.  As author of that code,
 *  I, Neal Davis, permit repurposing and redistribution for CLAMR under the New
 *  BSD License used above.
 *      http://code.google.com/p/capablanca/
 */
#include <fstream>
#include <iostream>
#include <cmath>
#include <cstdio>
#include <cstdlib>
#include <cstring>

#include "partition.h"
#include "mesh.h"

#ifdef HAVE_CONFIG_H
#include "config.h"
#endif

#define OUTPUT_INTERVAL 20
#define COARSE_GRID_RES 128
#define MAX_TIME_STEP 5000

using namespace std;

//  Global variables.
char progName[12];      //  Program name.
char progVers[8];       //  Program version.

//  External global variables.
extern bool verbose,
            localStencil,
            outline,
            enhanced_precision_sum;
extern int  outputInterval,
            tmax,
            levmx,
            nx,
            ny,
            niter,
            measure_type,
            calc_neighbor_type,
            initial_order,
            cycle_reorder;

//extern int  do_cpu_calc,
//            do_gpu_calc;

void outputHelp()
{   cout << "CLAMR is an experimental adaptive mesh refinement code for the GPU." << endl
         << "Version is " << PACKAGE_VERSION << endl << endl
         << "Usage:  " << progName << " [options]..." << endl
         << "  -c                turn on CPU profiling;" << endl
         << "  -g                turn on GPU profiling;" << endl
         << "  -h                display this help message;" << endl
         << "  -i <I>            specify I steps between output files;" << endl
         << "  -l <l>            max number of levels;" << endl
         << "  -m <m>            specify partition measure type;" << endl
         << "      \"with_duplicates\"" << endl
         << "      \"without_duplicates\"" << endl
         << "  -N <n>            specify calc neighbor type;" << endl
         << "      \"hash_table\"" << endl
         << "      \"kdtree\"" << endl
         << "  -n <N>            specify coarse grid resolution of NxN;" << endl
         << "  -o                turn off outlines;" << endl
         << "  -P <P>            specify initial order P;" << endl
         << "      \"original_order\"" << endl
         << "      \"hilbert_sort\"" << endl
         << "      \"hilbert_partition\"" << endl
         << "      \"z_order\"" << endl
         << "  -p <p>            specify ordering P every cycle;" << endl
         << "      \"original_order\"" << endl
         << "      \"hilbert_sort\"" << endl
         << "      \"hilbert_partition\"" << endl
         << "      \"local_hilbert\"" << endl
         << "      \"local_fixed\"" << endl
         << "      \"z_order\"" << endl
         << "  -r                regular sum instead of enhanced precision sum (Kahan sum);" << endl
         << "  -s <s>            specify space-filling curve method S;" << endl
         << "  -T                execute with TVD;" << endl
         << "  -t <t>            specify T time steps to run;" << endl
         << "  -V                use verbose output;" << endl
         << "  -v                display version information." << endl; }

void outputVersion()
{   cout << progName << " " << progVers << endl; }

/*  parseInput(const int argc, char** argv)
 *  
 *  Interpret the command line input.
 */
void parseInput(const int argc, char** argv)
{   strcpy(progName, "clamr");
    strcpy(progVers, PACKAGE_VERSION);
    
    //	Reconstruct command line argument as a string.
    char progCL[256];       //  Complete program command line.
    strcpy(progCL, argv[0]);
    for (int i = 1; i < argc; i++)
    {   strcat(progCL, " ");
        strcat(progCL, argv[i]); }
    
    //  Set variables to defaults, which may be overridden by CLI.
    verbose            = false;
    localStencil       = true;
    outline            = true;
    outputInterval     = OUTPUT_INTERVAL;
    nx                 = COARSE_GRID_RES;
    ny                 = COARSE_GRID_RES;
    niter              = MAX_TIME_STEP;
    measure_type       = CVALUE;
    calc_neighbor_type = HASH_TABLE;
    initial_order      = HILBERT_SORT;
    cycle_reorder      = ORIGINAL_ORDER;
    levmx              = 1;
    enhanced_precision_sum = true;
    
    char   *val;
    if (argc > 1)
    {   int i = 1;
        val = strtok(argv[i++], " ,.-");
        while (val != NULL)
        {   switch (val[0])
            {   case 'c':   //  Turn on CPU profiling.
                    //do_cpu_calc = 1;
                    break;
                    
                case 'g':   //  Turn on GPU profiling.
                    //do_gpu_calc = 1;
                    break;
                    
                case 'h':   //  Output help.
                    outputHelp();
                    cout.flush();
                    exit(EXIT_SUCCESS);
                    break;
                    
                case 'i':   //  Output interval specified.
                    val = strtok(argv[i++], " ,.-");
                    outputInterval = atoi(val);
                    break;
                    
                case 'l':   //  max level specified.
                    val = strtok(argv[i++], " ,");
                    levmx = atoi(val);
                    break;
                    
                case 'm':   //  partition measure specified.
                    val = strtok(argv[i++], " ,");
                    if (! strcmp(val,"with_duplicates") ) {
                       measure_type = WITH_DUPLICATES;
                    } else if (! strcmp(val,"without_duplicates") ) {
                       measure_type = WITHOUT_DUPLICATES;
                    } else if (! strcmp(val,"cvalue") ) {
                       measure_type = CVALUE;
                    }
                    break;
                    
                case 'N':   //  calc neighbor type specified.
                    val = strtok(argv[i++], " ,");
                    if (! strcmp(val,"hash_table") ) {
                       calc_neighbor_type = HASH_TABLE;
                    } else if (! strcmp(val,"kdtree") ) {
                       calc_neighbor_type = KDTREE;
                    }
                    break;
                    
                case 'n':   //  Domain grid resolution specified.
                    val = strtok(argv[i++], " ,");
                    nx = atoi(val);
                    ny = nx;
                    break;
                    
                case 'o':   //  Turn off outlines on mesh drawing.
                    outline = false;
                    break;
                    
                case 'P':   //  Initial order specified.
                    val = strtok(argv[i++], " ,");
                    if (! strcmp(val,"original_order") ) {
                       initial_order = ORIGINAL_ORDER;
                    } else if (! strcmp(val,"hilbert_sort") ) {
                       initial_order = HILBERT_SORT;
                    } else if (! strcmp(val,"hilbert_partition") ) {
                       initial_order = HILBERT_PARTITION;
                    } else if (! strcmp(val,"z_order") ) {
                       initial_order = ZORDER;
                    }
                    break;
                    
                case 'p':   //  Initial order specified.
                    val = strtok(argv[i++], " ,");
                    if (! strcmp(val,"original_order") ) {
                       cycle_reorder = ORIGINAL_ORDER;
                       localStencil = false;
                    } else if (! strcmp(val,"hilbert_sort") ) {
                       cycle_reorder = HILBERT_SORT;
                       localStencil = false;
                    } else if (! strcmp(val,"hilbert_partition") ) {
                       cycle_reorder = HILBERT_PARTITION;
                       localStencil = false;
                    } else if (! strcmp(val,"local_hilbert") ) {
                       cycle_reorder = ORIGINAL_ORDER;
                       localStencil = true;
                    } else if (! strcmp(val,"local_fixed") ) {
                       cycle_reorder = ORIGINAL_ORDER;
                       localStencil = false;
                    } else if (! strcmp(val,"z_order") ) {
                       cycle_reorder = ZORDER;
                       localStencil = false;
                    }
                    break;
                    
                case 'r':   //  Regular sum instead of enhanced precision sum.
                    enhanced_precision_sum = false;
                    break;
                    
                case 's':   //  Space-filling curve method specified (default HILBERT_SORT).
                //  Add different problem setups such as sloped wave in x, y and diagonal directions to help check algorithm
                    //  HILBERT_SORT
                    break;
                    
                case 'T':   //  TVD inclusion specified.
                    break;
                    
                case 't':   //  Number of time steps specified.
                    val = strtok(argv[i++], " ,.-");
                    niter = atoi(val);
                    break;
                    
                case 'V':   //  Verbose output desired.
                    verbose = true;
                    break;
                    
                case 'v':   //  Version.
                    outputVersion();
                    cout.flush();
                    exit(EXIT_SUCCESS);
                    break;
                    
                default:    //  Unknown parameter encountered.
                    cout << "⚠ Unknown input parameter " << val << endl;
                    outputHelp();
                    cout.flush();
                    exit(EXIT_FAILURE);
                    break; }
            
            val = strtok(argv[i++], " ,.-"); } } }
