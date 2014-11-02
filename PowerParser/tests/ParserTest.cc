#include <cstring>
#include "PowerParser.hh"
//#include "PowerParser/Comm.hh"
//#include "PowerParser/Parse.hh"
#include "genmalloc/genmalloc.h"

using Comm_ns::Comm;
using Support_ns::Parse;

using namespace std;

static Comm  *comm;
static Parse *parse;

int main(int argc, char **argv)
{
   printf("\n\t\tRunning the Parser read tests\n\n");

   // Setup parser's comm routines
   comm = new Comm;
   comm->initialize();
   int mype = comm->getProcRank();
   int npes = comm->getNumProcs();

   // Create the parser
   parse = new Parse(*comm);

   // Process input file
   string filename("parsetest.in",15);
   parse->parse_file(filename);

   parse->compile_buffer();
   parse->echo_input_start();

   string sline = "";
   int icount = 0;

/*
   while (parse->get_ssfout_line(sline) ) {
      printf("%s:%d: %s\n",&filename[0],icount++,&sline[0]);
   }
*/

   bool skip = false;
   vector<int> size;

   vector<string> sinput(1);
   string test_name("",40);
   string cname_stringinput("string_input");
   parse->get_char(cname_stringinput, sinput, size, false, skip);
   test_name = sinput[0];

   if(test_name.compare("parsetest") == 0 ){
      printf("  PASSED: Read string %s\n",test_name.c_str() );
   } else {
      printf("  FAILED: Read string %s should be %s\n",test_name.c_str(), "parsetest");
   }

   int intvalue = -1;
   string cname_intinput("int_input");
   parse->get_int(cname_intinput, &intvalue, size, skip);

   if (intvalue == 10){
      printf("  PASSED: Read int %d\n",intvalue);
   } else {
      printf("  FAILED: Read int %d should be %d\n",intvalue, 10);
   }
   
   double doublevalue = -1;
   string cname_doubleinput("double_input");
   parse->get_real(cname_doubleinput, &doublevalue, size, skip);

   if (doublevalue == 5.5){
      printf("  PASSED: Read real %lf\n",doublevalue);
   } else {
      printf("  FAILED: Read real %lf should be %lf\n",doublevalue, 5.5);
   }

   doublevalue = -1;
   cname_doubleinput = "dx";
   parse->get_real(cname_doubleinput, &doublevalue, size, skip);

   if (doublevalue == 1.0){
      printf("  PASSED: Read real %lf\n",doublevalue);
   } else {
      printf("  FAILED: Read real %lf should be %lf\n",doublevalue, 1.0);
   }

   doublevalue = -1;
   cname_doubleinput = "xreal";
   parse->get_real(cname_doubleinput, &doublevalue, size, skip);

   if (doublevalue == 1.0){
      printf("  PASSED: Read real %lf\n",doublevalue);
   } else {
      printf("  FAILED: Read real %lf should be %lf\n",doublevalue, 1.0);
   }

   intvalue = -1;
   cname_intinput = "ivalue";
   parse->get_int(cname_intinput, &intvalue, size, skip);

   if (intvalue == 100){
      printf("  PASSED: Read int %d\n",intvalue);
   } else {
      printf("  FAILED: Read int %d should be %d\n",intvalue, 100);
   }
   
   cname_stringinput = "string";
   parse->get_char(cname_stringinput, sinput, size, false, skip);
   test_name = sinput[0];

   if(test_name.compare("malarky") == 0 ){
      printf("  PASSED: Read string %s\n",test_name.c_str() );
   } else {
      printf("  FAILED: Read string %s should be %s\n",test_name.c_str(), "malarky");
   }

   bool boolvalue = false;
   string cname_boolinput("doThing");
   parse->get_bool(cname_boolinput, &boolvalue, size, skip);

   if (boolvalue == true){
      printf("  PASSED: Read boolean %s\n",boolvalue ? "true" : "false");
   } else {
      printf("  FAILED: Read bool %s should be %s\n",boolvalue ? "true" : "false", "true");
   }

   boolvalue = false;
   cname_boolinput = "doThing1";
   parse->get_bool(cname_boolinput, &boolvalue, size, skip);

   if (boolvalue == true){
      printf("  PASSED: Read boolean %s\n",boolvalue ? "true" : "false");
   } else {
      printf("  FAILED: Read bool %s should be %s\n",boolvalue ? "true" : "false", "true");
   }
   
   boolvalue = false;
   cname_boolinput = "doThing2";
   parse->get_bool(cname_boolinput, &boolvalue, size, skip);

   if (boolvalue == true){
      printf("  PASSED: Read boolean %s\n",boolvalue ? "true" : "false");
   } else {
      printf("  FAILED: Read bool %s should be %s\n",boolvalue ? "true" : "false", "true");
   }
   
   boolvalue = true;
   cname_boolinput = "doThing3";
   parse->get_bool(cname_boolinput, &boolvalue, size, skip);

   if (boolvalue == false){
      printf("  PASSED: Read boolean %s\n",boolvalue ? "true" : "false");
   } else {
      printf("  FAILED: Read bool %s should be %s\n",boolvalue ? "true" : "false", "false");
   }
   
   boolvalue = true;
   cname_boolinput = "doThing4";
   parse->get_bool(cname_boolinput, &boolvalue, size, skip);

   if (boolvalue == false){
      printf("  PASSED: Read boolean %s\n",boolvalue ? "true" : "false");
   } else {
      printf("  FAILED: Read bool %s should be %s\n",boolvalue ? "true" : "false", "false");
   }
   
   double doublearray[6] = {-1.0, -1.0, -1.0, -1.0, -1.0, -1.0};
   cname_doubleinput = "array1d";
   size.push_back(6);
   parse->get_real(cname_doubleinput, doublearray, size, skip);

   double dcheckarray[6] = {1.0, 2.3, -5.6, 7.1e19, 3.0, -3.4e-23};

   bool iflag = false;
   for (int i = 0; i < 6; i++){
      if (doublearray[i] != dcheckarray[i]) iflag = true;
   }

   if (! iflag){
      printf("  PASSED: Read real[0] %lg real[1] %lg real[2] %lg real[3] %lg real[4] %lg real[5] %lg\n",
             doublearray[0], doublearray[1], doublearray[2], doublearray[3], doublearray[4], doublearray[5]);
   } else {
      printf("  FAILED: Read real[0] %lg real[1] %lg real[2] %lg real[3] %lg real[4] %lg real[5] %lg\n",
             doublearray[0], doublearray[1], doublearray[2], doublearray[3], doublearray[4], doublearray[5]);
      printf("     should be real[0] %lg real[1] %lg real[2] %lg real[3] %lg real[4] %lg real[5] %lg\n",
             dcheckarray[0], dcheckarray[1], dcheckarray[2], dcheckarray[3], dcheckarray[4], dcheckarray[5]);
   }

   for (int i = 0; i < 6; i++){
      doublearray[i] = -1.0;
   }
   cname_doubleinput = "array1d_1";
   parse->get_real(cname_doubleinput, doublearray, size, skip);

   for (int i = 0; i < 6; i++){
      if (doublearray[i] != dcheckarray[i]) iflag = true;
   }

   if (! iflag){
      printf("  PASSED: Read real[0] %lg real[1] %lg real[2] %lg real[3] %lg real[4] %lg real[5] %lg\n",
             doublearray[0], doublearray[1], doublearray[2], doublearray[3], doublearray[4], doublearray[5]);
   } else {
      printf("  FAILED: Read real[0] %lg real[1] %lg real[2] %lg real[3] %lg real[4] %lg real[5] %lg\n",
             doublearray[0], doublearray[1], doublearray[2], doublearray[3], doublearray[4], doublearray[5]);
      printf("     should be real[0] %lg real[1] %lg real[2] %lg real[3] %lg real[4] %lg real[5] %lg\n",
             dcheckarray[0], dcheckarray[1], dcheckarray[2], dcheckarray[3], dcheckarray[4], dcheckarray[5]);
   }

   for (int i = 0; i < 6; i++){
      doublearray[i] = -1.0;
   }
   cname_doubleinput = "array1d_2";
   parse->get_real(cname_doubleinput, doublearray, size, skip);

   for (int i = 0; i < 6; i++){
      if (doublearray[i] != dcheckarray[i]) iflag = true;
   }

   if (! iflag){
      printf("  PASSED: Read real[0] %lg real[1] %lg real[2] %lg real[3] %lg real[4] %lg real[5] %lg\n",
             doublearray[0], doublearray[1], doublearray[2], doublearray[3], doublearray[4], doublearray[5]);
   } else {
      printf("  FAILED: Read real[0] %lg real[1] %lg real[2] %lg real[3] %lg real[4] %lg real[5] %lg\n",
             doublearray[0], doublearray[1], doublearray[2], doublearray[3], doublearray[4], doublearray[5]);
      printf("     should be real[0] %lg real[1] %lg real[2] %lg real[3] %lg real[4] %lg real[5] %lg\n",
             dcheckarray[0], dcheckarray[1], dcheckarray[2], dcheckarray[3], dcheckarray[4], dcheckarray[5]);
   }

   for (int i = 0; i < 6; i++){
      doublearray[i] = -1.0;
   }
   cname_doubleinput = "array1d_3";
   parse->get_real(cname_doubleinput, doublearray, size, skip);

   for (int i = 0; i < 6; i++){
      if (doublearray[i] != dcheckarray[i]) iflag = true;
   }

   if (! iflag){
      printf("  PASSED: Read real[0] %lg real[1] %lg real[2] %lg real[3] %lg real[4] %lg real[5] %lg\n",
             doublearray[0], doublearray[1], doublearray[2], doublearray[3], doublearray[4], doublearray[5]);
   } else {
      printf("  FAILED: Read real[0] %lg real[1] %lg real[2] %lg real[3] %lg real[4] %lg real[5] %lg\n",
             doublearray[0], doublearray[1], doublearray[2], doublearray[3], doublearray[4], doublearray[5]);
      printf("     should be real[0] %lg real[1] %lg real[2] %lg real[3] %lg real[4] %lg real[5] %lg\n",
             dcheckarray[0], dcheckarray[1], dcheckarray[2], dcheckarray[3], dcheckarray[4], dcheckarray[5]);
   }

   double **doublearray2d = (double **)genmatrix(2, 3, sizeof(double));
   cname_doubleinput = "array2d";

   size.clear();
   size.push_back(3);
   size.push_back(2);
   parse->get_real(cname_doubleinput, &doublearray2d[0][0], size, skip);

   double **dcheckarray2d = (double **)genmatrix(2, 3, sizeof(double));

   dcheckarray2d[0][0] = 1.0;
   dcheckarray2d[0][1] = 2.3;
   dcheckarray2d[0][2] = -5.6;
   dcheckarray2d[1][0] = 7.1e19;
   dcheckarray2d[1][1] = 3.0;
   dcheckarray2d[1][2] = -3.4e-23;
   
   for (int j = 0; j < size[1]; j++){
      for (int i = 0; i < size[0]; i++){
         if (doublearray2d[j][i] != dcheckarray2d[j][i]) iflag = true;
      }
   }

   if (! iflag){
      printf("  PASSED: Read real[0][0] %lg real[0][1] %lg real[0][2] %lg real[1][0] %lg real[1][1] %lg real[1][2] %lg\n",
             doublearray2d[0][0], doublearray2d[0][1], doublearray2d[0][2],
             doublearray2d[1][0], doublearray2d[1][1], doublearray2d[1][2]);
   } else {
      printf("  FAILED: Read real[0][0] %lg real[0][1] %lg real[0][2] %lg real[1][0] %lg real[1][1] %lg real[1][2] %lg\n",
             doublearray2d[0][0], doublearray2d[0][1], doublearray2d[0][2],
             doublearray2d[1][0], doublearray2d[1][1], doublearray2d[1][2]);
      printf("     should be real[0][0] %lg real[0][1] %lg real[0][2] %lg real[1][0] %lg real[1][1] %lg real[1][2] %lg\n",
             dcheckarray2d[0][0], dcheckarray2d[0][1], dcheckarray2d[0][2],
             dcheckarray2d[1][0], dcheckarray2d[1][1], dcheckarray2d[1][2]);
   }

   genmatrixfree((void **)doublearray2d);
   genmatrixfree((void **)dcheckarray2d);

   doublearray2d = (double **)genmatrix(3, 2, sizeof(double));
   cname_doubleinput = "array2d_1";

   size.clear();
   size.push_back(2);
   size.push_back(3);
   parse->get_real(cname_doubleinput, &doublearray2d[0][0], size, skip);

   dcheckarray2d = (double **)genmatrix(3, 2, sizeof(double));

   dcheckarray2d[0][0] = 1.0;
   dcheckarray2d[0][1] = 2.3;
   dcheckarray2d[1][0] = -5.6;
   dcheckarray2d[1][1] = 7.1e19;
   dcheckarray2d[2][0] = 3.0;
   dcheckarray2d[2][1] = -3.4e-23;
   
   for (int j = 0; j < size[1]; j++){
      for (int i = 0; i < size[0]; i++){
         if (doublearray2d[j][i] != dcheckarray2d[j][i]) iflag = true;
      }
   }

   if (! iflag){
      printf("  PASSED: Read real[0][0] %lg real[0][1] %lg real[1][0] %lg real[1][1] %lg real[2][0] %lg real[2][1] %lg\n",
             doublearray2d[0][0], doublearray2d[0][1],
             doublearray2d[1][0], doublearray2d[1][1],
             doublearray2d[2][0], doublearray2d[2][1]);
   } else {
      printf("  FAILED: Read real[0][0] %lg real[0][1] %lg real[1][0] %lg real[1][1] %lg real[2][0] %lg real[2][1] %lg\n",
             doublearray2d[0][0], doublearray2d[0][1],
             doublearray2d[1][0], doublearray2d[1][1],
             doublearray2d[2][0], doublearray2d[2][1]);
      printf("     should be real[0][0] %lg real[0][1] %lg real[1][0] %lg real[1][1] %lg real[2][0] %lg real[2][1] %lg\n",
             dcheckarray2d[0][0], dcheckarray2d[0][1],
             dcheckarray2d[1][0], dcheckarray2d[1][1],
             dcheckarray2d[2][0], dcheckarray2d[2][1]);
   }

   genmatrixfree((void **)doublearray2d);
   genmatrixfree((void **)dcheckarray2d);

   printf("\n\t\tFinished the Parser read tests\n\n");
   
}


