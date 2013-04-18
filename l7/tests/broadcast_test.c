#include <stdio.h>
#include <stdlib.h>
#include <string.h>
#include "l7.h"

void broadcast_test()
{
   int mype,
       numpes,
       ierr;
   
   char c[10];
   int i;
   long ilong;
   long long ilonglong;
   float xfloat;
   double xdouble;

   mype  =L7_Get_Rank();
   numpes=L7_Get_Numpes();

   if (mype == 0)
      printf("\n\t\tRunning the Broadcast tests\n\n");
   
   strcpy(c,"\0");
   if (mype == numpes-1) strcpy(c,"test\0");
   ierr = L7_Broadcast(c, 5, L7_CHAR, numpes-1 );
   if (mype ==0) {
       if (ierr != L7_OK || strncmp(c,"test",4) ){
        printf("  Error with L7_Broadcast of char type\n");
       }
       else{
        printf("  PASSED L7_Broadcast of char type\n");
       }
   }

   i=0;
   if (mype == numpes-1) i = 1;
   ierr = L7_Broadcast(&i, 1, L7_INT, numpes-1 );
   if (mype == 0){
      if (ierr != L7_OK || i != 1){
          printf("  Error with L7_Broadcast of int type\n");
       }
       else{
         printf("  PASSED L7_Broadcast of int type\n");
      }
   }
   
   ilong=0L;
   if (mype == numpes-1) ilong = 1L;
   ierr = L7_Broadcast(&ilong, 1, L7_LONG, numpes-1 );
   if (mype == 0){
      if (ierr != L7_OK || ilong != 1L){
         printf("  Error with L7_Broadcast of long type\n");
       }
       else{
          printf("  PASSED L7_Broadcast of long type\n");
       }
   }
   
   ilonglong=0L;
   if (mype == numpes-1) ilonglong = 1L;
   ierr = L7_Broadcast(&ilonglong, 1, L7_LONG_LONG_INT, numpes-1 );
   if (mype == 0){
      if (ierr != L7_OK || ilonglong != 1L){
         printf("  Error with L7_Broadcast of long long type\n");
       }
       else{
          printf("  PASSED L7_Broadcast of long long type\n");
       }
   }

   xfloat=0.0;
   if (mype == numpes-1) xfloat = 1.0;
   ierr = L7_Broadcast(&xfloat, 1, L7_FLOAT, numpes-1 );
   if (mype == 0){
      if (ierr != L7_OK || xfloat != 1.0){
         printf("  Error with L7_Broadcast of float type\n");
       }
       else{
          printf("  PASSED L7_Broadcast of float type\n");
       }
   }

   xdouble=0.0;
   if (mype == numpes-1) xdouble = 1.0;
   ierr = L7_Broadcast(&xdouble, 1, L7_DOUBLE, numpes-1 );
   if (mype == 0){
      if (ierr != L7_OK || xdouble != 1.0){
         printf("  Error with L7_Broadcast of double type\n");
       }
       else{
          printf("  PASSED L7_Broadcast of double type\n");
       }
   }

}
