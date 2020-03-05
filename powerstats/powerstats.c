/*
 *  Copyright (c) 2011-2019, Triad National Security, LLC.
 *  All rights Reserved.
 *
 *  CLAMR -- LA-CC-11-094
 *
 *  Copyright 2011-2019. Triad National Security, LLC. This software was produced 
 *  under U.S. Government contract 89233218CNA000001 for Los Alamos National 
 *  Laboratory (LANL), which is operated by Triad National Security, LLC 
 *  for the U.S. Department of Energy. The U.S. Government has rights to use, 
 *  reproduce, and distribute this software.  NEITHER THE GOVERNMENT NOR
 *  TRIAD NATIONAL SECURITY, LLC MAKES ANY WARRANTY, EXPRESS OR IMPLIED, OR 
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
 *     * Neither the name of the Triad National Security, LLC, Los Alamos 
 *       National Laboratory, LANL, the U.S. Government, nor the names of its 
 *       contributors may be used to endorse or promote products derived from 
 *       this software without specific prior written permission.
 *  
 *  THIS SOFTWARE IS PROVIDED BY THE TRIAD NATIONAL SECURITY, LLC AND 
 *  CONTRIBUTORS "AS IS" AND ANY EXPRESS OR IMPLIED WARRANTIES, INCLUDING, BUT 
 *  NOT LIMITED TO, THE IMPLIED WARRANTIES OF MERCHANTABILITY AND FITNESS FOR
 *  A PARTICULAR PURPOSE ARE DISCLAIMED. IN NO EVENT SHALL TRIAD NATIONAL
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
#include <stdlib.h>
#include <stdio.h>
#include <unistd.h>
#include <string.h>

#ifdef HAVE_CONFIG_H
#include "config.h"
#endif

#include "powerstats.h"

#if defined(HAVE_RAPL_PERF) || defined(HAVE_RAPL_SYSFS) || defined(HAVE_RAPL_MSR)
#define HAVE_RAPL 1
int detect_cpu(void);
int detect_packages(void);
#endif
#ifdef HAVE_RAPL_PERF
int rapl_perf(int core);
#endif
#ifdef HAVE_RAPL_SYSFS
int rapl_sysfs(int core);
#endif
#ifdef HAVE_RAPL_MSR
int rapl_msr(int core);
#endif
int (*rapl_read)(int core) = NULL;

void extract_power_data();

#ifdef HAVE_POWER_GADGET
   #include <EnergyLib.h>
#endif
int rapl_perf(int core);

void powerstats_init(void){
#ifdef HAVE_RAPL
#ifdef HAVE_RAPL_MSR
   detect_cpu();
#endif
   detect_packages();

#ifdef HAVE_RAPL_PERF
   if (rapl_perf(0) != -1){
      rapl_read = &rapl_perf;
   }
#endif

#ifdef HAVE_RAPL_SYSFS
   if (rapl_sysfs(0) != -1){
      printf("Found working RAPL_SYSFS\n");
      rapl_read = &rapl_sysfs;
   } else {
      printf("RAPL_SYSFS fails\n");
   }
#endif

#ifdef HAVE_RAPL_MSR
   if (rapl_msr(0) != -1){
      printf("Found working RAPL_MSR\n");
      rapl_read = &rapl_msr;
   } else {
      printf("RAPL_MSR fails\n");
   }
#endif

#endif

#ifdef HAVE_POWER_GADGET
   IntelEnergyLibInitialize();
   //StartLog((char *)"PowerGadgetLog.csv"); // causes a sample to be read
   ReadSample();
   ReadSample();

   extract_power_data();
#endif
}

void powerstats_sample(void){
#ifdef HAVE_POWER_GADGET
      ReadSample();
      extract_power_data();
#endif
}

void powerstats_finalize(void){
#ifdef HAVE_POWER_GADGET
      printf("Closing power gadget log file\n");
      //StopLog();// causes a sample to be read
      ReadSample();

      extract_power_data();
      printf("\n");
      printf("Processor      Energy(mWh) = %10.5f\n",  power_gadget_processor_mWh_sum);
      printf("IA             Energy(mWh) = %10.5f\n",  power_gadget_ia_mWh_sum);
      printf("DRAM           Energy(mWh) = %10.5f\n",  power_gadget_dram_mWh_sum);
      printf("Processor      Power (W)   = %10.5f\n",  power_gadget_processor_power_W_sum/power_gadget_time_secs_sum);
      printf("IA             Power (W)   = %10.5f\n",  power_gadget_ia_power_W_sum/power_gadget_time_secs_sum);
      printf("DRAM           Power (W)   = %10.5f\n",  power_gadget_dram_power_W_sum/power_gadget_time_secs_sum);
      printf("Average Frequency          = %10.5f\n",  power_gadget_avg_frequency/power_gadget_time_secs_sum);
      printf("Average Temperature (C)    = %10.5f\n",  power_gadget_avg_temperature/power_gadget_time_secs_sum);
      printf("Time Expended (secs)       = %10.5f\n",  power_gadget_time_secs_sum);

         double time_expended;
#endif
}

void extract_power_data(void)
{
#ifdef HAVE_POWER_GADGET
   int numMsrs = 0;
   GetNumMsrs(&numMsrs);

   for (int j = 0; j < numMsrs; j++) {
      int funcID;
      char szName[1024];
      GetMsrFunc(j, &funcID);
      GetMsrName(j, szName);
      int nData;
      double data[3];
      GetPowerData(0, j, data, &nData);
      double time_expended;

      // Frequency
      if (funcID == MSR_FUNC_FREQ) {
         GetTimeInterval(&time_expended);

         power_gadget_avg_frequency += data[0]*time_expended;;
         sample_count ++;
         if (POWER_GADGET_VERBOSE == 1) {
            printf("%-14s             = %4.0f\n",  szName, data[0]);
         }

         power_gadget_time_secs_sum += time_expended;
      }

      // Power
      else if (funcID == MSR_FUNC_POWER) {
         if (POWER_GADGET_VERBOSE == 1) {
            printf("%-14s Power (W)   = %10.5f\n", szName, data[0]);
            printf("%-14s Energy(J)   = %10.5f\n",  szName, data[1]);
            printf("%-14s Energy(mWh) = %10.5f\n",  szName, data[2]);
         }
         if (! strcmp(szName,"Processor")){
            power_gadget_processor_power_W_sum += data[0];
            power_gadget_processor_mWh_sum += data[2];
         } else if (! strcmp(szName,"IA")){
            power_gadget_ia_power_W_sum += data[0];
            power_gadget_ia_mWh_sum += data[2];
         } else if (! strcmp(szName,"DRAM")){
            power_gadget_dram_power_W_sum += data[0];
            power_gadget_dram_mWh_sum += data[2];
         }
      }

      // Temperature
      else if (funcID == MSR_FUNC_TEMP) {
         power_gadget_avg_temperature += data[0]*time_expended;
         if (POWER_GADGET_VERBOSE == 1) {
            printf("%-14s Temp (C)    = %4.0f\n",  szName, data[0]);
         }
      }
   }
   if (POWER_GADGET_VERBOSE == 1) {
      printf("\n");
   }
#endif
}
