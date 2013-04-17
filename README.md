The CLAMR code is a cell-based adaptive mesh refinement (AMR) mini-app developed
as a testbed for hybrid algorithm development using MPI and OpenCL GPU code. 

The CLAMR code is open-sourced under its LANL copyright
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
 *  See LICENSE file for full copyright 
 *  
 */

Contributions to CLAMR are welcomed as long as they do not substantially change
the nature of the code.

To build the CLAMR executables. CLAMR now used cmake for its builds

// In-tree build
cmake . 
// out-of-tree build
cmake <path-to-source>

// Optimized build (part optimized RelWithDebInfo is default)
cmake -DCMAKE_BUILD_TYPE=release .
// Graphics options (OpenGL is the default)
cmake -DGRAPHICS_TYPE=[None|OpenGL|MPE] <path-to-source>

There are two real-time graphics packages. The default is OpenGL. An alternative
real-time graphics package uses MPE. To use this package, configure with 
"cmake -DGRAPHICS_TYPE=MPE". The OpenGL option is automatically turned off when
selecting MPE. MPE is part of the MPICH package from Argonne National Laboratory.
It usually is not on a system and will need to be installed. A setup script and
version 1.9.1 of the MPE package are available in the download directory.

make

Seven executables are currently built. The first four are "standalone" versions which
run one implementation of the routines

clamr: Calls the MPI/GPU versions of each call. Option to check the results against
other verisons

clamr_gpuonly: Calls the GPU versions of each call. Option to check the results
against the cpu calls

clamr_cpuonly: Calls the CPU versions of each call. Option to check the results
against the gpu calls

clamr_mpionly: Calls the MPI versions of each call. Option to check the results
against the cpu calls

Check versions. These versions run multiple implementations and check correctness of
the implementations

clamr_gpucheck: Calls the GPU and CPU versions of each call and checks the results
against each other

clamr_mpicheck: Calls the CPU and MPI/CPU versions of each call and checks the results
against each other

clamr_checkall: Calls the GPU, CPU, MPI and GPU/MPI versions of each call and checks
the results against each other.

More executables are planned

Currently the executables run only on NVIDIA GPUs. Fixing the kernels to run on
ATI GPUs and MIC is of great interest

The numerical algorithm still does not handle "dry" conditions properly and will
crash

Many other limitations exist -- coarsening has not been implemented and boundary
conditions need some more work

Current performance shows about a 30x speedup on the GPU versus the CPU using NVIDIA
Tesla 2090s and Intel CPUs

See the PAPERS file for a list of publications related to the CLAMR code (Papers.bib for 
bibtex format)
