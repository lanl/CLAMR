#include <stdio.h>
#include <assert.h>
#include <string.h>
#include "MallocPlus.h"

#define VERBOSE 1

bool check_memory_dict(MallocPlus *my_mem,
                       double *mem_ptr_check,
                       const char *var_name_check,
                       int elemsize_check,
                       size_t nelems_check,
                       size_t ncapacity_check,
                       int flags_check)
{
   int elemsize          =           my_mem->get_memory_elemsize(mem_ptr_check);
   size_t nelems         =           my_mem->get_memory_size(mem_ptr_check);
   size_t ncapacity      =           my_mem->get_memory_capacity(mem_ptr_check);
   int flags             =           my_mem->get_memory_flags(mem_ptr_check);
   double *mem_ptr       = (double *)my_mem->get_memory_ptr(var_name_check);
   const char *var_name  =   (char *)my_mem->get_memory_name(mem_ptr_check);

   bool rflag = true;
   if (elemsize  == elemsize_check  &&
       nelems    == nelems_check    &&
       ncapacity == ncapacity_check &&
       flags     == flags_check     &&
       mem_ptr   == mem_ptr_check   &&
       !strcmp(var_name, var_name_check) ) {
      rflag = true;
      if (VERBOSE) {
         printf("PASSED: elemsize %d nelems is %lu capacity %lu flags %d\n",
                elemsize, nelems, ncapacity, flags);
         printf("        mem_ptr  = %p mem_ptr_check  = %p\n", mem_ptr, mem_ptr_check);
         printf("        var_name = %s var_name_check = %s\n", var_name, var_name_check);
      }
   } else {
      rflag = false;
      if (VERBOSE) {
         printf(" Failed: elemsize %d nelems is %lu capacity %lu flags %d\n",
                elemsize, nelems, ncapacity, flags);
         printf("        mem_ptr  = %p mem_ptr_check  = %p\n", mem_ptr, mem_ptr_check);
         printf("        var_name = %s var_name_check = %s\n", var_name, var_name_check);
      }
   }

   assert(elemsize  == elemsize_check);
   assert(nelems    == nelems_check);
   assert(ncapacity == ncapacity_check);
   assert(flags     == flags_check);
   assert(mem_ptr   == mem_ptr_check);
   assert(!strcmp(var_name, var_name_check));

   return(rflag);
}

bool check_memory_dict(MallocPlus *my_mem,
                       float *mem_ptr_check,
                       const char *var_name_check,
                       int elemsize_check,
                       size_t nelems_check,
                       size_t ncapacity_check,
                       int flags_check)
{
   int elemsize          =           my_mem->get_memory_elemsize(mem_ptr_check);
   size_t nelems         =           my_mem->get_memory_size(mem_ptr_check);
   size_t ncapacity      =           my_mem->get_memory_capacity(mem_ptr_check);
   int flags             =           my_mem->get_memory_flags(mem_ptr_check);
   float *mem_ptr        =  (float *)my_mem->get_memory_ptr(var_name_check);
   const char *var_name  =   (char *)my_mem->get_memory_name(mem_ptr_check);

   bool rflag = true;
   if (elemsize  == elemsize_check  &&
       nelems    == nelems_check    &&
       ncapacity == ncapacity_check &&
       flags     == flags_check     &&
       mem_ptr   == mem_ptr_check   &&
       !strcmp(var_name, var_name_check) ) {
      rflag = true;
      if (VERBOSE) {
         printf("PASSED: elemsize %d nelems is %lu capacity %lu flags %d\n",
                elemsize, nelems, ncapacity, flags);
         printf("        mem_ptr  = %p mem_ptr_check  = %p\n", mem_ptr, mem_ptr_check);
         printf("        var_name = %s var_name_check = %s\n", var_name, var_name_check);
      }
   } else {
      rflag = false;
      if (VERBOSE) {
         printf(" Failed: elemsize %d nelems is %lu capacity %lu flags %d\n",
                elemsize, nelems, ncapacity, flags);
         printf("        mem_ptr  = %p mem_ptr_check  = %p\n", mem_ptr, mem_ptr_check);
         printf("        var_name = %s var_name_check = %s\n", var_name, var_name_check);
      }
   }

   assert(elemsize  == elemsize_check);
   assert(nelems    == nelems_check);
   assert(ncapacity == ncapacity_check);
   assert(flags     == flags_check);
   assert(mem_ptr   == mem_ptr_check);
   assert(!strcmp(var_name, var_name_check));

   return(rflag);
}

bool check_memory_dict(MallocPlus *my_mem,
                       int *mem_ptr_check,
                       const char *var_name_check,
                       int elemsize_check,
                       size_t nelems_check,
                       size_t ncapacity_check,
                       int flags_check)
{
   int elemsize          =           my_mem->get_memory_elemsize(mem_ptr_check);
   size_t nelems         =           my_mem->get_memory_size(mem_ptr_check);
   size_t ncapacity      =           my_mem->get_memory_capacity(mem_ptr_check);
   int flags             =           my_mem->get_memory_flags(mem_ptr_check);
   int *mem_ptr          =    (int *)my_mem->get_memory_ptr(var_name_check);
   const char *var_name  =   (char *)my_mem->get_memory_name(mem_ptr_check);

   bool rflag = true;
   if (elemsize  == elemsize_check  &&
       nelems    == nelems_check    &&
       ncapacity == ncapacity_check &&
       flags     == flags_check     &&
       mem_ptr   == mem_ptr_check   &&
       !strcmp(var_name, var_name_check) ) {
      rflag = true;
      if (VERBOSE) {
         printf("PASSED: elemsize %d nelems is %lu capacity %lu flags %d\n",
                elemsize, nelems, ncapacity, flags);
         printf("        mem_ptr  = %p mem_ptr_check  = %p\n", mem_ptr, mem_ptr_check);
         printf("        var_name = %s var_name_check = %s\n", var_name, var_name_check);
      }
   } else {
      rflag = false;
      if (VERBOSE) {
         printf(" Failed: elemsize %d nelems is %lu capacity %lu flags %d\n",
                elemsize, nelems, ncapacity, flags);
         printf("        mem_ptr  = %p mem_ptr_check  = %p\n", mem_ptr, mem_ptr_check);
         printf("        var_name = %s var_name_check = %s\n", var_name, var_name_check);
      }
   }

   assert(elemsize  == elemsize_check);
   assert(nelems    == nelems_check);
   assert(ncapacity == ncapacity_check);
   assert(flags     == flags_check);
   assert(mem_ptr   == mem_ptr_check);
   assert(!strcmp(var_name, var_name_check));

   return(rflag);
}


int main(int argc, char **argv)
{
   printf("\n\t\tRunning the MallocPlus tests\n\n");

   printf("\n   Starting test of simple malloc and delete by pointer\n");
   {
      // This tests a simple malloc and delete by pointer for memory errors.
      // Uses:
      //   new MallocPlus
      //   memory_malloc
      //   memory_delete
      //   delete my_mem
      double *X;
      MallocPlus *my_mem = new MallocPlus();
      X = (double *)my_mem->memory_malloc(10, sizeof(double), "X");
      check_memory_dict(my_mem, X, "X", sizeof(double), 10, 10, HOST_REGULAR_MEMORY);
      X = (double *)my_mem->memory_delete(X);
      delete my_mem;
   }
   printf("   Ending simple malloc test\n\n");

   printf("\n   Starting test of simple malloc and delete by name\n");
   {
      // This tests a simple malloc and delete by name for memory errors.
      // Uses:
      //   new MallocPlus
      //   memory_malloc
      //   memory_delete
      //   delete my_mem
      double *X;
      MallocPlus *my_mem = new MallocPlus();
      X = (double *)my_mem->memory_malloc(10, sizeof(double), "X");
      check_memory_dict(my_mem, X, "X", sizeof(double), 10, 10, HOST_REGULAR_MEMORY);
      X = (double *)my_mem->memory_delete("X");
      delete my_mem;
   }
   printf("   Ending simple malloc test\n\n");

   printf("\n   Starting test of managed memory malloc and delete\n");
   {
      // This tests a managed memory malloc and delete for memory errors
      // Uses:
      //   new MallocPlus
      //   memory_malloc with flags
      //   memory_realloc
      //   memory_delete
      //   delete my_mem
      double *X;
      int flags = HOST_MANAGED_MEMORY;
      MallocPlus *my_mem = new MallocPlus();
      X = (double *)my_mem->memory_malloc(10, sizeof(double), "X", flags);
      check_memory_dict(my_mem, X, "X", sizeof(double), 10, 20, HOST_MANAGED_MEMORY);
      X = (double *)my_mem->memory_realloc(30, "X");
      check_memory_dict(my_mem, X, "X", sizeof(double), 30, 60, HOST_MANAGED_MEMORY);
      X = (double *)my_mem->memory_realloc(100, X);
      check_memory_dict(my_mem, X, "X", sizeof(double), 100, 200, HOST_MANAGED_MEMORY);
      X = (double *)my_mem->memory_delete(X);
      delete my_mem;
   }
   printf("   Ending managed memory malloc\n\n");

   printf("\n   Starting realloc of named variable\n");
   {
      // This tests a named realloc for memory errors
      // Uses
      //   new MallocPlus
      //   memory_malloc
      //   memory_realloc
      //   memory_delete
      //   delete my_mem
      double *X;
      MallocPlus *my_mem = new MallocPlus();
      X = (double *)my_mem->memory_malloc(10, sizeof(double), "X");
      X = (double *)my_mem->memory_realloc(20, "X");
      check_memory_dict(my_mem, X, "X", sizeof(double), 20, 20, HOST_REGULAR_MEMORY);
      X = (double *)my_mem->memory_delete(X);
      delete my_mem;
   }
   printf("   Ending named realloc test\n\n");

   printf("\n   Starting pointer realloc\n");
   {
      // This tests the pointer realloc for memory errors
      // Uses
      //   new MallocPlus
      //   memory_malloc
      //   memory_realloc
      //   memory_delete
      //   delete my_mem
      double *X;
      MallocPlus *my_mem = new MallocPlus();
      X = (double *)my_mem->memory_malloc(10, sizeof(double), "X");
      X = (double *)my_mem->memory_realloc(20, X);
      check_memory_dict(my_mem, X, "X", sizeof(double), 20, 20, HOST_REGULAR_MEMORY);
      X = (double *)my_mem->memory_delete(X);
      delete my_mem;
   }
   printf("   Ending pointer realloc test\n\n");

   printf("\n   Starting pointer realloc_all\n");
   {
      // This tests the pointer realloc for memory errors
      // Uses
      //   new MallocPlus
      //   memory_malloc
      //   memory_realloc_all
      //   get_memory_ptr
      //   memory_delete
      //   delete my_mem
      double *X, *Y;
      MallocPlus *my_mem = new MallocPlus();
      X = (double *)my_mem->memory_malloc(10, sizeof(double), "X");
      Y = (double *)my_mem->memory_malloc(20, sizeof(double), "Y");
      my_mem->memory_realloc_all(30);
      // Must get the memory pointers since they have changed!
      X = (double *)my_mem->get_memory_ptr("X");
      Y = (double *)my_mem->get_memory_ptr("Y");
      check_memory_dict(my_mem, X, "X", sizeof(double), 30, 30, HOST_REGULAR_MEMORY);
      check_memory_dict(my_mem, Y, "Y", sizeof(double), 30, 30, HOST_REGULAR_MEMORY);
      X = (double *)my_mem->memory_delete(X);
      Y = (double *)my_mem->memory_delete(Y);
      delete my_mem;
   }
   printf("   Ending pointer realloc_all test\n\n");

   printf("\n   Starting test of managed memory request all\n");
   {
      // This tests a managed memory malloc and delete for memory errors
      // Uses:
      //   new MallocPlus
      //   memory_malloc with flags
      //   get_memory_ptr
      //   memory_request
      //   memory_delete
      //   delete my_mem
      double *X, *Y;
      int flags = HOST_MANAGED_MEMORY;
      MallocPlus *my_mem = new MallocPlus();
      X = (double *)my_mem->memory_malloc(10, sizeof(double), "X", flags);
      Y = (double *)my_mem->memory_malloc(20, sizeof(double), "Y", flags);
      my_mem->memory_request_all(50);
      X = (double *)my_mem->get_memory_ptr("X");
      Y = (double *)my_mem->get_memory_ptr("Y");
      check_memory_dict(my_mem, X, "X", sizeof(double), 10, 50, HOST_MANAGED_MEMORY);
      check_memory_dict(my_mem, Y, "Y", sizeof(double), 20, 50, HOST_MANAGED_MEMORY);
      X = (double *)my_mem->memory_request(200, X);
      Y = (double *)my_mem->memory_request(200, "Y");
      check_memory_dict(my_mem, X, "X", sizeof(double), 10, 200, HOST_MANAGED_MEMORY);
      check_memory_dict(my_mem, Y, "Y", sizeof(double), 20, 200, HOST_MANAGED_MEMORY);
      X = (double *)my_mem->memory_delete(X);
      Y = (double *)my_mem->memory_delete(Y);
      delete my_mem;
   }
   printf("   Ending managed memory malloc\n\n");

   printf("\n   Starting test of simple malloc and duplicate\n");
   {
      // This tests a simple malloc and duplicate.
      // Uses:
      //   new MallocPlus
      //   memory_malloc
      //   memory_duplicate
      //   delete my_mem
      double *X, *X_new;
      MallocPlus *my_mem = new MallocPlus();
      X = (double *)my_mem->memory_malloc(10, sizeof(double), "X");
      check_memory_dict(my_mem, X, "X", sizeof(double), 10, 10, HOST_REGULAR_MEMORY);
      X_new = (double *)my_mem->memory_duplicate(X, "X_new");
      check_memory_dict(my_mem, X_new, "X_new", sizeof(double), 10, 10, HOST_REGULAR_MEMORY);
      X = (double *)my_mem->memory_delete(X);
      X_new = (double *)my_mem->memory_delete(X_new);
      delete my_mem;
   }
   printf("   Ending simple malloc and duplicate\n\n");

   printf("\n   Starting pointer swap\n");
   {
      // This tests the pointer realloc for memory errors
      // Uses
      //   new MallocPlus
      //   memory_malloc
      //   memory_swap
      //   memory_delete
      //   delete my_mem
      double *X, *Y;
      MallocPlus *my_mem = new MallocPlus();
      X = (double *)my_mem->memory_malloc(10, sizeof(double), "X");
      Y = (double *)my_mem->memory_malloc(20, sizeof(double), "Y");
      my_mem->memory_swap(&X, &Y);
      check_memory_dict(my_mem, X, "X", sizeof(double), 20, 20, HOST_REGULAR_MEMORY);
      check_memory_dict(my_mem, Y, "Y", sizeof(double), 10, 10, HOST_REGULAR_MEMORY);
      X = (double *)my_mem->memory_delete(X);
      Y = (double *)my_mem->memory_delete(Y);
      delete my_mem;
   }
   printf("   Ending pointer swap test\n\n");

   printf("\n   Starting float pointer swap\n");
   {
      // This tests the pointer realloc for memory errors
      // Uses
      //   new MallocPlus
      //   memory_malloc
      //   memory_swap
      //   memory_delete
      //   delete my_mem
      float *X, *Y;
      MallocPlus *my_mem = new MallocPlus();
      X = (float *)my_mem->memory_malloc(10, sizeof(float), "X");
      Y = (float *)my_mem->memory_malloc(20, sizeof(float), "Y");
      my_mem->memory_swap(&X, &Y);
      check_memory_dict(my_mem, X, "X", sizeof(float), 20, 20, HOST_REGULAR_MEMORY);
      check_memory_dict(my_mem, Y, "Y", sizeof(float), 10, 10, HOST_REGULAR_MEMORY);
      X = (float *)my_mem->memory_delete(X);
      Y = (float *)my_mem->memory_delete(Y);
      delete my_mem;
   }
   printf("   Ending float pointer swap test\n\n");

   printf("\n   Starting int pointer swap\n");
   {
      // This tests the pointer realloc for memory errors
      // Uses
      //   new MallocPlus
      //   memory_malloc
      //   memory_realloc
      //   memory_delete
      //   delete my_mem
      int *X, *Y;
      MallocPlus *my_mem = new MallocPlus();
      X = (int *)my_mem->memory_malloc(10, sizeof(int), "X");
      Y = (int *)my_mem->memory_malloc(20, sizeof(int), "Y");
      my_mem->memory_swap(&X, &Y);
      check_memory_dict(my_mem, X, "X", sizeof(int), 20, 20, HOST_REGULAR_MEMORY);
      check_memory_dict(my_mem, Y, "Y", sizeof(int), 10, 10, HOST_REGULAR_MEMORY);
      X = (int *)my_mem->memory_delete(X);
      Y = (int *)my_mem->memory_delete(Y);
      delete my_mem;
   }
   printf("   Ending int pointer swap test\n\n");

   printf("\n   Starting pointer replace\n");
   {
      // This tests the pointer replace
      // Uses
      //   new MallocPlus
      //   memory_malloc
      //   memory_replace
      //   memory_delete
      //   delete my_mem
      double *X, *X_new;
      MallocPlus *my_mem = new MallocPlus();
      X = (double *)my_mem->memory_malloc(10, sizeof(double), "X");
      X_new = (double *)my_mem->memory_malloc(20, sizeof(double), "Y");
      X = (double *)my_mem->memory_replace(X, X_new);
      check_memory_dict(my_mem, X, "X", sizeof(double), 20, 20, HOST_REGULAR_MEMORY);
      X = (double *)my_mem->memory_delete(X);
      delete my_mem;
   }
   printf("   Ending pointer replace test\n\n");

   printf("\n   Starting test of simple malloc, add, remove, and delete\n");
   {
      // This tests a simple malloc and delete by pointer for memory errors.
      // Uses:
      //   new MallocPlus
      //   memory_malloc
      //   memory_add
      //   memory_remove
      //   delete my_mem
      double *X;
      MallocPlus *my_mem = new MallocPlus();
      X = (double *)malloc(10*sizeof(double));
      my_mem->memory_add(X, 10, sizeof(double), "X");
      check_memory_dict(my_mem, X, "X", sizeof(double), 10, 10, HOST_REGULAR_MEMORY);
      my_mem->memory_remove(X);
      free(X);
      delete my_mem;
   }
   printf("   Ending simple malloc and and remove test\n\n");

   printf("\n   Starting test of simple malloc, add, remove, and delete by name\n");
   {
      // This tests a simple malloc and delete by pointer for memory errors.
      // Uses:
      //   new MallocPlus
      //   memory_malloc
      //   memory_add
      //   memory_remove
      //   delete my_mem
      double *X;
      MallocPlus *my_mem = new MallocPlus();
      X = (double *)malloc(10*sizeof(double));
      my_mem->memory_add(X, 10, sizeof(double), "X");
      check_memory_dict(my_mem, X, "X", sizeof(double), 10, 10, HOST_REGULAR_MEMORY);
      my_mem->memory_remove("X");
      free(X);
      delete my_mem;
   }
   printf("   Ending simple malloc and and remove test\n\n");

   printf("\n   Starting test of simple malloc set attributes\n");
   {
      // This tests a simple malloc set attributes.
      // Uses:
      //   new MallocPlus
      //   memory_malloc
      //   set_memory_attribute
      //   check_memory_attribute
      //   clear_memory_attribute
      //   memory_delete
      //   delete my_mem
      double *X;
      MallocPlus *my_mem = new MallocPlus();
      X = (double *)my_mem->memory_malloc(10, sizeof(double), "X");
      my_mem->set_memory_attribute(X,HOST_MANAGED_MEMORY);
      check_memory_dict(my_mem, X, "X", sizeof(double), 10, 10, HOST_MANAGED_MEMORY);
      bool flag_set = my_mem->check_memory_attribute(X, HOST_MANAGED_MEMORY);
      assert(flag_set);
      my_mem->clear_memory_attribute(X, HOST_MANAGED_MEMORY);
      check_memory_dict(my_mem, X, "X", sizeof(double), 10, 10, HOST_REGULAR_MEMORY);
      X = (double *)my_mem->memory_delete(X);
      delete my_mem;
   }
   printf("   Ending simple malloc test\n\n");

   printf("\n   Starting check of delete_all\n");
   {
      // This tests the pointer realloc for memory errors
      // Uses
      //   new MallocPlus
      //   memory_malloc
      //   memory_realloc_all
      //   get_memory_ptr
      //   memory_delete_all
      //   delete my_mem
      double *X, *Y;
      MallocPlus *my_mem = new MallocPlus();
      X = (double *)my_mem->memory_malloc(10, sizeof(double), "X");
      Y = (double *)my_mem->memory_malloc(20, sizeof(double), "Y");
      my_mem->memory_realloc_all(30);
      // Must get the memory pointers since they have changed!
      X = (double *)my_mem->get_memory_ptr("X");
      Y = (double *)my_mem->get_memory_ptr("Y");
      check_memory_dict(my_mem, X, "X", sizeof(double), 30, 30, HOST_REGULAR_MEMORY);
      check_memory_dict(my_mem, Y, "Y", sizeof(double), 30, 30, HOST_REGULAR_MEMORY);
      my_mem->memory_delete_all();
      X = NULL;
      Y = NULL;
      delete my_mem;
   }
   printf("   Ending check of delete_all test\n\n");

   printf("\n   Starting test of reorder of double array\n");
   {
      // This tests a reorder of an array of doubles.
      // Uses:
      //   new MallocPlus
      //   memory_malloc
      //   memory_reorder
      //   memory_delete
      //   delete my_mem
      double *X;
      MallocPlus *my_mem = new MallocPlus();
      X = (double *)my_mem->memory_malloc(10, sizeof(double), "X");
      int *iorder = (int *)malloc(10*sizeof(int));
      for (int i = 0; i < 10; i++){
         X[i] = (double)(9-i);
         iorder[i]=(9-i);
      }
      check_memory_dict(my_mem, X, "X", sizeof(double), 10, 10, HOST_REGULAR_MEMORY);
      X = my_mem->memory_reorder(X, iorder);
      check_memory_dict(my_mem, X, "X", sizeof(double), 10, 10, HOST_REGULAR_MEMORY);
      for (int i = 0; i < 10; i++){
         assert(i == (int)X[i]);
      }
      X = (double *)my_mem->memory_delete(X);
      delete my_mem;
      free(iorder);
   }
   printf("   Ending reorder test\n\n");

   printf("\n   Starting test a reorder of an array of floats\n");
   {
      // This tests a reorder of an array of floats.
      // Uses:
      //   new MallocPlus
      //   memory_malloc
      //   memory_reorder
      //   memory_delete
      //   delete my_mem
      float *X;
      MallocPlus *my_mem = new MallocPlus();
      X = (float *)my_mem->memory_malloc(10, sizeof(float), "X");
      int *iorder = (int *)malloc(10*sizeof(int));
      for (int i = 0; i < 10; i++){
         X[i] = (float)(9-i);
         iorder[i]=(9-i);
      }
      check_memory_dict(my_mem, X, "X", sizeof(float), 10, 10, HOST_REGULAR_MEMORY);
      X = my_mem->memory_reorder(X, iorder);
      check_memory_dict(my_mem, X, "X", sizeof(float), 10, 10, HOST_REGULAR_MEMORY);
      for (int i = 0; i < 10; i++){
         assert(i == (int)X[i]);
      }
      X = (float *)my_mem->memory_delete(X);
      delete my_mem;
      free(iorder);
   }
   printf("   Ending reorder test of an array of floats\n\n");

   printf("\n   Starting test a reorder of an array of ints\n");
   {
      // This tests a reorder of an array of ints.
      // Uses:
      //   new MallocPlus
      //   memory_malloc
      //   memory_reorder
      //   memory_delete
      //   delete my_mem
      int *X;
      MallocPlus *my_mem = new MallocPlus();
      X = (int *)my_mem->memory_malloc(10, sizeof(int), "X");
      int *iorder = (int *)malloc(10*sizeof(int));
      for (int i = 0; i < 10; i++){
         X[i] = 9-i;
         iorder[i]=(9-i);
      }
      check_memory_dict(my_mem, X, "X", sizeof(int), 10, 10, HOST_REGULAR_MEMORY);
      X = my_mem->memory_reorder(X, iorder);
      check_memory_dict(my_mem, X, "X", sizeof(int), 10, 10, HOST_REGULAR_MEMORY);
      for (int i = 0; i < 10; i++){
         assert(i == X[i]);
      }
      X = (int *)my_mem->memory_delete(X);
      delete my_mem;
      free(iorder);
   }
   printf("   Ending reorder test of an array of ints\n\n");

   printf("\n   Starting test of reorder_all of double arrays\n");
   {
      // This tests a reorder_all of two arrays of doubles.
      // Uses:
      //   new MallocPlus
      //   memory_malloc
      //   memory_reorder_all
      //   memory_delete
      //   delete my_mem
      double *X, *Y;
      MallocPlus *my_mem = new MallocPlus();
      X = (double *)my_mem->memory_malloc(10, sizeof(double), "X");
      Y = (double *)my_mem->memory_malloc(10, sizeof(double), "Y");
      int *iorder = (int *)malloc(10*sizeof(int));
      for (int i = 0; i < 10; i++){
         X[i] = (double)(9-i);
         Y[i] = (double)(2*(9-i));
         iorder[i]=(9-i);
      }
      check_memory_dict(my_mem, X, "X", sizeof(double), 10, 10, HOST_REGULAR_MEMORY);
      check_memory_dict(my_mem, Y, "Y", sizeof(double), 10, 10, HOST_REGULAR_MEMORY);
      my_mem->memory_reorder_all(iorder);
      X = (double *)my_mem->get_memory_ptr("X");
      Y = (double *)my_mem->get_memory_ptr("Y");
      check_memory_dict(my_mem, X, "X", sizeof(double), 10, 10, HOST_REGULAR_MEMORY);
      check_memory_dict(my_mem, Y, "Y", sizeof(double), 10, 10, HOST_REGULAR_MEMORY);
      for (int i = 0; i < 10; i++){
         assert(i == (int)X[i]);
         assert(i == (int)(Y[i]/2.0));
      }
      my_mem->memory_delete_all();
      delete my_mem;
      free(iorder);
   }
   printf("   Ending reorder test\n\n");

   printf("\n   Starting test a reorder of an index array of ints\n");
   {
      // This tests a reorder of an index array of ints.
      // Uses:
      //   new MallocPlus
      //   memory_malloc
      //   memory_reorder_indexarray
      //   memory_delete
      //   delete my_mem
      int *X;
      MallocPlus *my_mem = new MallocPlus();
      X = (int *)my_mem->memory_malloc(10, sizeof(int), "X");
      int *iorder = (int *)malloc(10*sizeof(int));
      int *inv_iorder = (int *)malloc(10*sizeof(int));
      for (int i = 0; i < 10; i++){
         X[i] = 9-i;
         iorder[i]=(9-i);
         inv_iorder[i]=i;
      }
      check_memory_dict(my_mem, X, "X", sizeof(int), 10, 10, HOST_REGULAR_MEMORY);
      X = my_mem->memory_reorder_indexarray(X, iorder, inv_iorder);
      check_memory_dict(my_mem, X, "X", sizeof(int), 10, 10, HOST_REGULAR_MEMORY);
      for (int i = 0; i < 10; i++){
         assert(i == X[i]);
      }
      X = (int *)my_mem->memory_delete(X);
      delete my_mem;
      free(iorder);
      free(inv_iorder);
   }
   printf("   Ending reorder test of an index array of ints\n\n");

   printf("\n   Starting begin, end loop\n");
   {
      // This tests the begin, end loop
      // Uses
      //   new MallocPlus
      //   memory_malloc
      //   memory_begin
      //   memory_next
      //   get_memory_elemsize
      //   memory_delete
      //   delete my_mem
      double *X, *Y;
      MallocPlus *my_mem = new MallocPlus();
      X = (double *)my_mem->memory_malloc(10, sizeof(double), "X");
      Y = (double *)my_mem->memory_malloc(20, sizeof(double), "Y");
      for (void *mem_ptr = my_mem->memory_begin(); mem_ptr != NULL; mem_ptr = my_mem->memory_next() ){
         int elsize = my_mem->get_memory_elemsize(mem_ptr);
         assert(elsize == sizeof(double) );
      }
      check_memory_dict(my_mem, X, "X", sizeof(double), 10, 10, HOST_REGULAR_MEMORY);
      check_memory_dict(my_mem, Y, "Y", sizeof(double), 20, 20, HOST_REGULAR_MEMORY);
      X = (double *)my_mem->memory_delete(X);
      Y = (double *)my_mem->memory_delete(Y);
      delete my_mem;
   }
   printf("   Ending begin, end loop test\n\n");

   printf("\n\t\tFinished the MallocPlus tests\n\n");
}

// x void *memory_malloc(size_t nelem, size_t elsize, const char *name, int flags=0);
// x void *memory_duplicate(void *malloc_mem_ptr, const char *addname);
// x void *memory_realloc(size_t nelem, void *malloc_mem_ptr);
// x void *memory_realloc(size_t nelem, const char *name);
// x void *memory_request(size_t new_capacity, void *malloc_mem_ptr);
// x void *memory_request(size_t new_capacity, const char *name);
// x void memory_realloc_all(size_t nelem);
// x void memory_request_all(size_t new_capacity);
// x void *memory_replace(void *malloc_mem_ptr_old, void * const malloc_mem_ptr_new);
// x void memory_swap(int **malloc_mem_ptr_old, int **malloc_mem_ptr_new);
// x void memory_swap(float **malloc_mem_ptr_old, float **malloc_mem_ptr_new);
// x void memory_swap(double **malloc_mem_ptr_old, double **malloc_mem_ptr_new);
// x void *memory_add(void *malloc_mem_ptr, size_t nelem, size_t elsize,
//      const char *name, int flags=0);
// x double *memory_reorder(double *malloc_mem_ptr, int *iorder);
// x float *memory_reorder(float *malloc_mem_ptr, int *iorder);
// x int *memory_reorder(int *malloc_mem_ptr, int *iorder);
// x int *memory_reorder_indexarray(int *malloc_mem_ptr, int *iorder, int *inv_iorder);
// x void memory_reorder_all(int *iorder);
//   void memory_report(void);
// x void *memory_delete(void *malloc_mem_ptr);
// x void *memory_delete(const char *name);
// x void memory_delete_all(void);
// x void memory_remove(void *malloc_mem_ptr);
// x void memory_remove(const char *name);
//   void *memory_begin(void);
//   void *memory_next(void);
//   malloc_plus_memory_entry *memory_entry_begin(void);
//   malloc_plus_memory_entry *memory_entry_next(void);
//   malloc_plus_memory_entry *memory_entry_end(void);
// x size_t get_memory_size(void *malloc_mem_ptr);
// x size_t get_memory_capacity(void *malloc_mem_ptr);
// x int get_memory_elemsize(void *malloc_mem_ptr);
// x const char *get_memory_name(void *malloc_mem_ptr);
// x void *get_memory_ptr(const char *name);
// x void set_memory_attribute(void *malloc_mem_ptr, int attribute);
// x void clear_memory_attribute(void *malloc_mem_ptr, int attribute);
// x int  get_memory_flags(void *malloc_mem_ptr);
// x bool  check_memory_attribute(void *malloc_mem_ptr, int attribute);
