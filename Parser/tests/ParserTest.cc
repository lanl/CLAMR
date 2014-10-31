#include <cstring>
#include "Parser/Comm.hh"
#include "Parser/Parse.hh"

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

   printf("\n\t\tFinished the Parser read tests\n\n");
   
}


