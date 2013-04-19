#include <stdio.h>
#include "l7.h"

void reduction_test()
{
   int mype,
       numpes,
       ierr;
   
   int i, iout, ivalout;
   long ilong, ilongout;
   long long ilonglong, ilonglongout;
   float xfloat, xfloatout;
   double xdouble, xdoubleout;

   mype  =L7_Get_Rank();
   numpes=L7_Get_Numpes();

   if (mype == 0)
      printf("\n\t\tRunning the Reduction tests\n\n");

   i = 1;
   L7_Sum(&i, 1, L7_INT, &iout);
   if (mype == 0){
      if (iout != numpes){
         printf("  Error with L7_Sum of int type\n");
       }
       else{
          printf("  PASSED L7_Sum of int type\n");
       }
   }
   
   ilong = 1L;
   L7_Sum(&ilong, 1, L7_LONG, &ilongout);
   if (mype == 0){
      if (ilongout != numpes){
         printf("  Error with L7_Sum of long type\n");
       }
       else{
          printf("  PASSED L7_Sum of long type\n");
       }
   }
   
   ilonglong = 1L;
   L7_Sum(&ilonglong, 1, L7_LONG_LONG_INT, &ilonglongout);
   if (mype == 0){
      if (ilonglongout != numpes){
         printf("  Error with L7_Sum of long long type\n");
       }
       else{
          printf("  PASSED L7_Sum of long long type\n");
       }
   }
   
   xfloat = 1.0;
   L7_Sum(&xfloat, 1, L7_FLOAT, &xfloatout);
   if (mype == 0){
      if (xfloatout != (float)numpes){
         printf("  Error with L7_Sum of float type\n");
       }
       else{
          printf("  PASSED L7_Sum of float type\n");
       }
   }

   xdouble = 1.0;
   L7_Sum(&xdouble, 1, L7_DOUBLE, &xdoubleout);
   if (mype == 0){
      if (xdoubleout != (double)numpes){
         printf("  Error with L7_Sum of double type\n");
       }
       else{
          printf("  PASSED L7_Sum of double type\n");
       }
   }
   
   i = 1;
   if (mype == numpes-1) i = 5;
   ierr = L7_Max(&i, 1, L7_INT, &iout);
   if (mype == 0){
      if (ierr != L7_OK || iout != 5){
         printf("  Error with L7_Max of int type\n");
       }
       else{
          printf("  PASSED L7_Max of int type\n");
       }
   }

   ilong = 1L;
   if (mype == numpes-1) ilong = 5L;
   ierr = L7_Max(&ilong, 1, L7_LONG, &ilongout);
   if (mype == 0){
      if (ierr != L7_OK || ilongout != 5L){
         printf("  Error with L7_Max of long type\n");
       }
       else{
          printf("  PASSED L7_Max of long type\n");
       }
   }

   ilonglong = 1L;
   if (mype == numpes-1) ilonglong = 5L;
   ierr = L7_Max(&ilonglong, 1, L7_LONG_LONG_INT, &ilonglongout);
   if (mype == 0){
      if (ierr != L7_OK || ilonglongout != 5L){
         printf("  Error with L7_Max of long long type\n");
       }
       else{
          printf("  PASSED L7_Max of long long type\n");
       }
   }

   xfloat = 1.0;
   if (mype == numpes-1) xfloat = 5.0;
   ierr = L7_Max(&xfloat, 1, L7_FLOAT, &xfloatout);
   if (mype == 0){
      if (ierr != L7_OK || xfloatout != 5.0){
         printf("  Error with L7_Max of float type\n");
       }
       else{
          printf("  PASSED L7_Max of float type\n");
       }
   }

   xdouble = 1.0;
   if (mype == numpes-1) xdouble = 5.0;
   ierr = L7_Max(&xdouble, 1, L7_DOUBLE, &xdoubleout);
   if (mype == 0){
      if (ierr != L7_OK || xdoubleout != 5.0){
         printf("  Error with L7_Max of double type\n");
       }
       else{
          printf("  PASSED L7_Max of double type\n");
       }
   }

   i = 10;
   if (mype == numpes-1) i = 5;
   ierr = L7_Min(&i, 1, L7_INT, &iout);
   if (mype == 0){
      if (ierr != L7_OK || iout != 5){
         printf("  Error with L7_Min of int type\n");
       }
       else{
          printf("  PASSED L7_Min of int type\n");
       }
   }

   ilong = 10L;
   if (mype == numpes-1) ilong = 5L;
   ierr = L7_Min(&ilong, 1, L7_LONG, &ilongout);
   if (mype == 0){
      if (ierr != L7_OK || ilongout != 5L){
         printf("  Error with L7_Min of long type\n");
       }
       else{
          printf("  PASSED L7_Min of long type\n");
       }
   }

   ilonglong = 10L;
   if (mype == numpes-1) ilonglong = 5L;
   ierr = L7_Min(&ilonglong, 1, L7_LONG_LONG_INT, &ilonglongout);
   if (mype == 0){
      if (ierr != L7_OK || ilonglongout != 5L){
         printf("  Error with L7_Min of long long type\n");
       }
       else{
          printf("  PASSED L7_Min of long long type\n");
       }
   }

   xfloat = 10.0;
   if (mype == numpes-1) xfloat = 5.0;
   ierr = L7_Min(&xfloat, 1, L7_FLOAT, &xfloatout);
   if (mype == 0){
      if (ierr != L7_OK || xfloatout != 5.0){
         printf("  Error with L7_Min of float type\n");
       }
       else{
          printf("  PASSED L7_Min of float type\n");
       }
   }

   xdouble = 10.0;
   if (mype == numpes-1) xdouble = 5.0;
   ierr = L7_Min(&xdouble, 1, L7_DOUBLE, &xdoubleout);
   if (mype == 0){
      if (ierr != L7_OK || xdoubleout != 5.0){
         printf("  Error with L7_Min of double type\n");
       }
       else{
          printf("  PASSED L7_Min of double type\n");
       }
   }

   i = 0;
   if (mype == numpes-1) i = 1;
   ierr = L7_Any(&i, 1, L7_INT, &iout);
   if (mype == 0){
      if (ierr != L7_OK || iout != 1){
         printf("  Error with L7_Any of int type\n");
       }
       else{
          printf("  PASSED L7_Any of int type\n");
       }
   }

   ilong = 0L;
   if (mype == numpes-1) ilong = 1L;
   ierr = L7_Any(&ilong, 1, L7_LONG, &ilongout);
   if (mype == 0){
      if (ierr != L7_OK || ilongout != 1L){
         printf("  Error with L7_Any of long type\n");
       }
       else{
          printf("  PASSED L7_Any of long type\n");
       }
   }

   ilonglong = 0L;
   if (mype == numpes-1) ilonglong = 1L;
   ierr = L7_Any(&ilonglong, 1, L7_LONG_LONG_INT, &ilonglongout);
   if (mype == 0){
      if (ierr != L7_OK || ilonglongout != 1L){
         printf("  Error with L7_Any of long long type\n");
       }
       else{
          printf("  PASSED L7_Any of long long type\n");
       }
   }

   i = 1;
   if (mype == numpes-1) i = 0;
   ierr = L7_All(&i, 1, L7_INT, &iout);
   if (mype == 0){
      if (ierr != L7_OK || iout != 0){
         printf("  Error with L7_All of int type\n");
       }
       else{
          printf("  PASSED L7_All of int type\n");
       }
   }

   ilong = 1L;
   if (mype == numpes-1) ilong = 0L;
   ierr = L7_All(&ilong, 1, L7_LONG, &ilongout);
   if (mype == 0){
      if (ierr != L7_OK || ilongout != 0L){
         printf("  Error with L7_All of long type\n");
       }
       else{
          printf("  PASSED L7_All of long type\n");
       }
   }

   ilonglong = 1L;
   if (mype == numpes-1) ilonglong = 0L;
   ierr = L7_All(&ilonglong, 1, L7_LONG_LONG_INT, &ilonglongout);
   if (mype == 0){
      if (ierr != L7_OK || ilonglongout != 0L){
         printf("  Error with L7_All of long long type\n");
       }
       else{
          printf("  PASSED L7_All of long long type\n");
       }
   }

   i = 0;
   if (mype == numpes-1) i = 1;
   ierr = L7_MaxLoc(&i, 1, L7_INT, &iout);
   if (mype == 0){
      if (ierr != L7_OK || iout != numpes-1){
         printf("  Error with L7_MaxLoc of int type\n");
       }
       else{
          printf("  PASSED L7_MaxLoc of int type\n");
       }
   }

   ilong = 0L;
   if (mype == numpes-1) ilong = 1L;
   ierr = L7_MaxLoc(&i, 1, L7_LONG, &iout);
   if (mype == 0){
      if (ierr != L7_OK || iout != numpes-1){
         printf("  Error with L7_MaxLoc of long type\n");
       }
       else{
          printf("  PASSED L7_MaxLoc of long type\n");
       }
   }

   xfloat = 0.0;
   if (mype == numpes-1) xfloat = 1.0;
   ierr = L7_MaxLoc(&xfloat, 1, L7_FLOAT, &iout);
   if (mype == 0){
      if (ierr != L7_OK || iout != numpes-1){
         printf("  Error with L7_MaxLoc of float type\n");
       }
       else{
          printf("  PASSED L7_MaxLoc of float type\n");
       }
   }


   xdouble = 0.0;
   if (mype == numpes-1) xdouble = 1.0;
   ierr = L7_MaxLoc(&xdouble, 1, L7_DOUBLE, &iout);
   if (mype == 0){
      if (ierr != L7_OK || iout != numpes-1){
         printf("  Error with L7_MaxLoc of double type\n");
       }
       else{
          printf("  PASSED L7_MaxLoc of double type\n");
       }
   }

   
   i = 1;
   if (mype == numpes-1) i = 0;
   ierr = L7_MinLoc(&i, 1, L7_INT, &iout);
   if (mype == 0){
      if (ierr != L7_OK || iout != numpes-1){
         printf("  Error with L7_MinLoc of int type\n");
       }
       else{
          printf("  PASSED L7_MinLoc of int type\n");
       }
   }

   ilong = 1L;
   if (mype == numpes-1) ilong = 0L;
   ierr = L7_MinLoc(&ilong, 1, L7_LONG, &iout);
   if (mype == 0){
      if (ierr != L7_OK || iout != numpes-1){
         printf("  Error with L7_MinLoc of long type\n");
       }
       else{
          printf("  PASSED L7_MinLoc of long type\n");
       }
   }

   xfloat = 1.0;
   if (mype == numpes-1) xfloat = 0.0;
   ierr = L7_MinLoc(&xfloat, 1, L7_FLOAT, &iout);
   if (mype == 0){
      if (ierr != L7_OK || iout != numpes-1){
         printf("  Error with L7_MinLoc of float type\n");
       }
       else{
          printf("  PASSED L7_MinLoc of float type\n");
       }
   }


   xdouble = 1.0;
   if (mype == numpes-1) xdouble = 0.0;
   ierr = L7_MinLoc(&xdouble, 1, L7_DOUBLE, &iout);
   if (mype == 0){
      if (ierr != L7_OK || iout != numpes-1){
         printf("  Error with L7_MinLoc of double type\n");
       }
       else{
          printf("  PASSED L7_MinLoc of double type\n");
       }
   }

   
   i = 0;
   if (mype == numpes-1) i = 1;
   ierr = L7_MaxValLoc(&i, 1, L7_INT, &ivalout, &iout);
   if (mype == 0){
      if (ierr != L7_OK || iout != numpes-1 || ivalout != 1){
         printf("  Error with L7_MaxValLoc of int type\n");
       }
       else{
          printf("  PASSED L7_MaxValLoc of int type\n");
       }
   }

   ilong = 0L;
   if (mype == numpes-1) ilong = 1L;
   ierr = L7_MaxValLoc(&ilong, 1, L7_LONG, &ilongout, &iout);
   if (mype == 0){
      if (ierr != L7_OK || iout != numpes-1 || ilongout != 1){
         printf("  Error with L7_MaxValLoc of long type\n");
       }
       else{
          printf("  PASSED L7_MaxValLoc of long type\n");
       }
   }

   xfloat = 0.0;
   if (mype == numpes-1) xfloat = 1.0;
   ierr = L7_MaxValLoc(&xfloat, 1, L7_FLOAT, &xfloatout, &iout);
   if (mype == 0){
      if (ierr != L7_OK || iout != numpes-1 || xfloatout != 1.0){
         printf("  Error with L7_MaxValLoc of float type\n");
       }
       else{
          printf("  PASSED L7_MaxValLoc of float type\n");
       }
   }

   xdouble = 0.0;
   if (mype == numpes-1) xdouble = 1.0;
   ierr = L7_MaxValLoc(&xdouble, 1, L7_DOUBLE, &xdoubleout, &iout);
   if (mype == 0){
      if (ierr != L7_OK || iout != numpes-1 || xdoubleout != 1.0){
         printf("  Error with L7_MaxValLoc of double type\n");
       }
       else{
          printf("  PASSED L7_MaxValLoc of double type\n");
       }
   }

   i = 1;
   if (mype == numpes-1) i = 0;
   ierr = L7_MinValLoc(&i, 1, L7_INT, &ivalout, &iout);
   if (mype == 0){
      if (ierr != L7_OK || iout != numpes-1 || ivalout != 0){
         printf("  Error with L7_MinValLoc of int type\n");
       }
       else{
          printf("  PASSED L7_MinValLoc of int type\n");
       }
   }

   ilong = 1L;
   if (mype == numpes-1) ilong = 0L;
   ierr = L7_MinValLoc(&ilong, 1, L7_LONG, &ilongout, &iout);
   if (mype == 0){
      if (ierr != L7_OK || iout != numpes-1 || ilongout != 0){
         printf("  Error with L7_MinValLoc of long type\n");
       }
       else{
          printf("  PASSED L7_MinValLoc of long type\n");
       }
   }

   xfloat = 1.0;
   if (mype == numpes-1) xfloat = 0.0;
   ierr = L7_MinValLoc(&xfloat, 1, L7_FLOAT, &xfloatout, &iout);
   if (mype == 0){
      if (ierr != L7_OK || iout != numpes-1 || xfloatout != 0.0){
         printf("  Error with L7_MinValLoc of float type\n");
       }
       else{
          printf("  PASSED L7_MinValLoc of float type\n");
       }
   }

   xdouble = 1.0;
   if (mype == numpes-1) xdouble = 0.0;
   ierr = L7_MinValLoc(&xdouble, 1, L7_DOUBLE, &xdoubleout, &iout);
   if (mype == 0){
      if (ierr != L7_OK || iout != numpes-1 || xdoubleout != 0.0){
         printf("  Error with L7_MinValLoc of double type\n");
       }
       else{
          printf("  PASSED L7_MinValLoc of double type\n");
       }
   }

}
