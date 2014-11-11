#include <stdio.h>
#include "MallocPlus.h"

int main(int argc, char **argv)
{
   printf("\n\t\tRunning the MallocPlus tests\n\n");

   MallocPlus *my_mem = new MallocPlus();

   my_mem->memory_report();

   double *X;
   X = (double *)my_mem->memory_malloc(10, sizeof(double), "X");
   my_mem->memory_report();
   X = (double *)my_mem->memory_delete(X);

   my_mem->memory_report();

   delete my_mem;


   printf("\n\t\tFinished the MallocPlus tests\n\n");
}
