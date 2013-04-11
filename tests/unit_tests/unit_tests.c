
#include <stdio.h>
#include <string.h>
#include <CUnit/Basic.h>
#include "../../lib/global_defs.h"
#include "../../lib/material_properties.h"
#include "../../lib/output.h"
#include "../../lib/parallel_related.h"



int layers_r(struct All_variables *E,float a)
{
	return 0;
}

FILE* output_open(char* refstate_file, char* a)
{
	return 0;
}

void parallel_process_termination()
{

}


/* The suite initialization function.
 * Opens the temporary file used by the tests.
 * Returns zero on success, non-zero otherwise.
 */
int init_suite1(void)
{
	return 0;
}

/* The suite cleanup function.
 * Closes the temporary file used by the tests.
 * Returns zero on success, non-zero otherwise.
 */
int clean_suite1(void)
{
    return 0;
}

/* Simple test of idxNZ().
 *
 *
 */
void testidxNz(void)
{
      CU_ASSERT(0 == idxNz(0, 5));
      CU_ASSERT(-1 == idxNz(10,0));
      CU_ASSERT(1 == idxNz(1, 10));
      CU_ASSERT(1 == idxNz(11, 10));
      CU_ASSERT(10 == idxNz(10,10));
}

void testidxTemp(void)
{
      CU_ASSERT(1 == idxTemp(0.0, 15.0, 10));
      CU_ASSERT(-1 == idxTemp(10.0, 0, 5));
      CU_ASSERT(-1 == idxTemp(10.0, 5.0, 0));
      CU_ASSERT(4 == idxTemp(10.0, 1.0, 5));
      CU_ASSERT(6 == idxTemp(10.0, 2.0, 20));
      CU_ASSERT(9 == idxTemp(20.0, 2.0, 10));
}

/*
void testCompressorString(void)
{
	char compressor_string[128];
	get_compressor_string(1,compressor_string);
    CU_ASSERT(0 == strcmp(compressor_string,""));

	get_compressor_string(0,compressor_string);
    CU_ASSERT(0 == strcmp(compressor_string," compressor=\"vtkZLibDataCompressor\" byte_order=\"LittleEndian\""));

}*/


/* The main() function for setting up and running the tests.
 * Returns a CUE_SUCCESS on successful running, another
 * CUnit error code on failure.
 */
int main()
{
   CU_pSuite pSuite = NULL;

   /* initialize the CUnit test registry */
   if (CUE_SUCCESS != CU_initialize_registry())
      return CU_get_error();

   /* add a suite to the registry */
   pSuite = CU_add_suite("Suite_1", init_suite1, clean_suite1);
   if (NULL == pSuite) {
      CU_cleanup_registry();
      return CU_get_error();
   }

   /* add the tests to the suite */
   /* NOTE - ORDER IS IMPORTANT - MUST TEST fread() AFTER fprintf() */
   if (NULL == CU_add_test(pSuite, "test of idxNz()", testidxNz) ||
		   NULL == CU_add_test(pSuite, "test of idxTemp()", testidxTemp)  )
   {
      CU_cleanup_registry();
      return CU_get_error();
   }

   /* Run all tests using the CUnit Basic interface */
   CU_basic_set_mode(CU_BRM_VERBOSE);
   CU_basic_run_tests();
   CU_cleanup_registry();
   return CU_get_error();
}


