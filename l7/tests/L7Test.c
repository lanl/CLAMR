#include "stdio.h"
#include "stdlib.h"
#include "l7.h"
#include "unistd.h"
#include <mpi.h>

void broadcast_test();
void reduction_test();
void update_test();

int nchars;
char buf[20];
int main(int argc, char *argv[])
{
    int mype =0,
        numpes=0,
        ierr;

   long long location;
   double xarray[10];

   int do_quo_setup=0;
   int lttrace_on=0;
   ierr = L7_Init(&mype, &numpes, &argc, argv, do_quo_setup, lttrace_on);

   if (mype == 0)
      printf("\n\t\tStarting the L7 tests\n\n");

   if (mype == 0){
      if (ierr != L7_OK){
          printf("  Error with L7_Init\n");
      }
      else{
          printf("  PASSED L7_Init\n");
      }
   }

   broadcast_test();

   reduction_test();

   update_test();

#ifdef XXX
   location = L7_Address(xarray);
   if (mype == 0){
      if (location != xarray){
          printf("  Error with L7_Address\n");
      }
      else{
          printf("  PASSED L7_Address\n");
      }
   }
#endif

   if (mype == 0)
      printf("\n\t\tRunning the L7_Terminate test\n\n");

   ierr = L7_Terminate();
   if (mype == 0){
      if (ierr != L7_OK){
          printf("  Error with L7_Terminate\n");
      }
      else{
          printf("  PASSED L7_Terminate\n");
      }
   }

   if (mype == 0)
      printf("\n\t\tFinished the L7 tests\n\n");

   exit(0);
}
