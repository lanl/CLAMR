#include <stdio.h>
#include "MallocPlus.h"

int main(int argc, char **argv)
{
   printf("\n\t\tRunning the MallocPlus tests\n\n");

   MallocPlus *my_mem = new MallocPlus();

   double *X;
   X = (double *)my_mem->memory_malloc(10, sizeof(double), "X");
   X = (double *)my_mem->memory_delete(X);


   delete my_mem;


   printf("\n\t\tFinished the MallocPlus tests\n\n");
}
