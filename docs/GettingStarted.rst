===============
Getting Started
===============

--------------
Building CLAMR
--------------

------------
Requirements
------------

MPI

OpenCL capable device -- Nvidia GPU

OpenMP capable compiler

Optional
   MPE graphics (for alternate real-time graphics)
      MPE is part of the MPICH package from Argonne National Laboratory.
      It usually is not on a system and will need to be installed. A setup script and
      version 1.9.1 of the MPE package are available in the download directory.

      Set MPEHOME to location of MPE install or MPE_INCLUDE_DIR and MPE_LIBRARIES to 
      the include path and the location of the MPE libraries

   ImageMagick (for output of graphics post-processing files)

   mplayer (mencode utility to create movies)

--------
Checkout
--------

To checkout the code

% git clone git@github.com:losalamos/CLAMR.git

---------
Configure
---------

*CLAMR* uses cmake so the configure step is 

    cmake .                         // in-tree build
    cmake <path-to-src>             // out-of-tree build

    cmake -DCMAKE_BUILD_TYPE=RelWithDebInfo <path-to-src> (currently the default)
    cmake -DCMAKE_BUILD_TYPE=Release <path-to-src>
    cmake -DCMAKE_BUILD_TYPE=Debug <path-to-src>

There are two real-time graphics packages. The default is no real-time graphics. The
recommended, primary version is OpenGL. An alternative real-time graphics package
uses MPE. To use MPE, configure with "cmake -DGRAPHICS_TYPE=MPE". The OpenGL
option is automatically turned off when selecting MPE. 

    cmake -DGRAPHICS_TYPE=None   <path-to-src>  // Turn off real-time graphics (currently the default)
    cmake -DGRAPHICS_TYPE=OpenGL <path-to-src>  // Primary real-time graphics
    cmake -DGRAPHICS_TYPE=MPE    <path-to-src>  // Alternate real-time graphics

    cmake -DPRECISION_TYPE=full_precision     // full double precision
    cmake -DPRECISION_TYPE=mixed_precision    // intermediates double, arrays single
    cmake -DPRECISION_TYPE=minimum_precision  // all single precision

    cmake -DMIC_NATIVE=yes   // For MIC Knights Corner accelerator cross-compile

-----
Build
-----

make

Seven executables are currently built. The first four are "standalone" versions which
run one implementation of the routines

 * clamr: Runs the MPI/GPU versions on multiple GPUs.

 * clamr_gpuonly: Runs the GPU version on a single GPU.

 * clamr_cpuonly: Runs a standard CPU version on a single CPU.

 * clamr_mpionly: Runs a MPI parallel version on multiple shared or distribured
   memory CPUs.

 * clamr_openmponly: Runs the OpenMP version on multiple shared-memory CPUs.

 * clamr_mpiopenmponly: Runs the MPI/OpenMP version on hybrid shared/distributed
   systems. This is the version used on the Intel MIC Knight's Corner accelerator.

Check versions. These versions run multiple implementations and check correctness of
the implementations

 * clamr_gpucheck: Calls the GPU and CPU versions of each call and checks the results
   against each other

 * clamr_mpicheck: Calls the CPU and MPI/CPU versions of each call and checks the results
   against each other

 * clamr_checkall: Calls the GPU, CPU, MPI and GPU/MPI versions of each call and checks
   the results against each other.

-------------
Special Cases
-------------

^^^^^^^^^^^^^^^^
Darwin MIC-NODE:
^^^^^^^^^^^^^^^^
ssh darwin
salloc -p knc-mic
ssh [NODE]
// where [NODE] is $SLURM_NODELIST
// Load modules for build
module load compilers/intel/15.0.0
module load mpi/intelmpi-4.1.3.045-mic
module load cmake/2.8.11.1

/// These can be setup in the .login for the Intel MICs on darwin with
module load cmake/2.8.11.1    compilers/intel/15.0.0    mpi/intelmpi-4.1.3.045-mic
//The following source command can also be put in the .login with the first version for
// csh and the second for sh. When these are needed and what they do is still being studied
source /projects/opt/intel/compilers/composer_xe_2015.0.090/bin/compilervars.csh intel64
source /projects/opt/intel/compilers/composer_xe_2015.0.090/bin/compilervars.sh intel64

// For building on the host, change the module load to drop the -mic on the mpi
 module load cmake/2.8.11.1    compilers/intel/15.0.0    mpi/intelmpi-4.1.3.045

// Build CLAMR on MIC:
cd CLAMR
cmake -DMIC_NATIVE=yes
make
// To execute:
ssh [NODE]-mic0

// Execute clamr_openmponly
export OMP_NUM_THREADS=120
In the CLAMR directory
cd CLAMR
./clamr_openmponly


-----
Notes
-----

Currently the executables run on NVIDIA GPUs and MICs.

The numerical algorithm still does not handle "dry" conditions properly and will
crash under some resolutions

Current performance shows about a 30x speedup on the GPU versus the CPU using NVIDIA
Tesla 2090s and Intel CPUs

See the PAPERS file for a list of publications related to the CLAMR code (Papers.bib for 
bibtex format)
