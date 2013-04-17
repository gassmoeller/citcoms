
#include <stdio.h>
#include <string.h>
#include <CUnit/Basic.h>
#include "../../lib/global_defs.h"
#include "../../lib/material_properties.h"
#include "../../lib/output.h"
#include "../../lib/parallel_related.h"
#include "../../lib/material_properties_perplex.h"


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

void testReadFieldOrder(void)
{

    struct IDs* read_field_order(char* bufline, int nfields, struct IDs *);

    //int *testfieldids;
    int i;
    char *testline = malloc (255*sizeof(char));

    snprintf(testline,255,"rho,kg/m3      alpha,1/K      cp,J/K/kg      vp,km/s        vs,km/s");

    struct IDs testfieldids;

    testfieldids = *(read_field_order(testline,5,&(testfieldids)));

    CU_ASSERT(0 == testfieldids.rho);
    CU_ASSERT(1 == testfieldids.alpha);
    CU_ASSERT(2 == testfieldids.cp);
    CU_ASSERT(3 == testfieldids.vp);
    CU_ASSERT(4 == testfieldids.vs);

    snprintf(testline,255,"vs,km/s  alpha,1/K   rho,kg/m3       cp,J/K/kg      vp,km/s       ");
    testfieldids = *(read_field_order(testline,5,&(testfieldids)));

    CU_ASSERT(2 == testfieldids.rho);
    CU_ASSERT(1 == testfieldids.alpha);
    CU_ASSERT(3 == testfieldids.cp);
    CU_ASSERT(4 == testfieldids.vp);
    CU_ASSERT(0 == testfieldids.vs);

    free(testline);
}

void testReadPerplexHeader(void)
{

    void read_perplex_header (FILE *perplex_file,
            struct table_properties *perplex_table);

    struct table_properties test_props;

    FILE *test_file = fopen("pyrolite.tab","r");
    read_perplex_header (test_file, &test_props);
    fclose(test_file);


    CU_ASSERT(5 == test_props.ninput_fields);
    CU_ASSERT(361 == test_props.TP[0].ndeps);
    CU_ASSERT(361 == test_props.TP[test_props.field_ids.T].ndeps);

    CU_ASSERT(262 == test_props.TP[1].ndeps);
    CU_ASSERT(fabs(test_props.TP[1].delta- 4999.9900000050002) <= 1e-14);

    test_file = fopen("pyrolite_reversed.tab","r");
    read_perplex_header (test_file, &test_props);
    fclose(test_file);

    CU_ASSERT(5 == test_props.ninput_fields);
    CU_ASSERT(361 == test_props.TP[1].ndeps);
    CU_ASSERT(361 == test_props.TP[test_props.field_ids.T].ndeps);

    CU_ASSERT(262 == test_props.TP[0].ndeps);

    CU_ASSERT(fabs(test_props.TP[0].delta- 4999.9900000050002) <= 1e-14);
}

void maketesttable(struct table_properties *table)
{
    table->ninput_fields = 3;
    table->field_ids.P = 1;
    table->field_ids.T = 0;
    table->field_ids.rho = 0;
    table->field_ids.alpha = 1;
    table->field_ids.cp = 2;

    int i;
    for (i=0;i<2;i++)
    {
    table->TP[i].ndeps = 1001;
    table->TP[i].start = 0;
    table->TP[i].delta = 1;
    table->TP[i].end = 1000;
    }
}

void maketestdata(struct table_properties *table, double ***data)
{
    int i;
    int j;
    for (i=0;i<table->TP[0].ndeps;i++)
        for (j=0;j<table->TP[1].ndeps;j++)
        {
            data[i][j][table->field_ids.rho] = 3000;
            data[i][j][table->field_ids.alpha] = 1e-4;
            data[i][j][table->field_ids.cp] = 1000;
        }
}

void testCalcAdiabaticConditionsConst(void)
{
    double Padi[2];
    double Tadi[2];

    struct table_properties table;
    maketesttable(&table);

    double ***data;
    int i,j;
    data = (double ***) malloc(table.TP[0].ndeps * sizeof(double **));
    for (i = 0;i<table.TP[0].ndeps;i++)
        {
        data[i] = (double **) malloc(table.TP[1].ndeps * sizeof(double *));
        for (j=0;j<table.TP[1].ndeps;j++)
            data[i][j] = (double *) malloc(table.ninput_fields * sizeof(double));
        }

    maketestdata(&table,data);

    Padi[1] = 0.0;
    Tadi[1] = 500.0;

    calc_adiabatic_conditions(&table,(const double ***)data,0,1000.0,10.0,&Padi,&Tadi);

    CU_ASSERT(fabs(300   - Padi[0]) < 1e-5);
    CU_ASSERT(fabs(500.500250 - Tadi[0]) < 1e-5);

}


//TODO: MAke Initial_temperature.c linkable to remove this code copy
double add_layer_point(double zInterface, double r1, double half_width, double maxdT)
{
    double trans;
    double dT = 0.0;
    if (r1 <= zInterface + 3 * half_width){
        trans = 0.5 * (1.0 + tanh((zInterface-r1) / half_width));
        dT += trans * maxdT;
    }
    return dT;
}

void testAddLayerPoint(void)
{
    double testdT;

    testdT = add_layer_point(0.5,0.4,0.01,1.0);
    fprintf(stderr,"%f\n",testdT);
    CU_ASSERT(fabs(1.0 - testdT) <= 1e-7);

    testdT = add_layer_point(0.5,0.5,0.01,1.0);
    fprintf(stderr,"%f\n",testdT);
    CU_ASSERT(fabs(0.5 - testdT) <= 1e-7);

    testdT = add_layer_point(0.5,0.6,0.01,1.0);
    fprintf(stderr,"%f\n",testdT);
    CU_ASSERT(fabs(0.0 - testdT) <= 1e-7);

    testdT = add_layer_point(0.5,0.51,0.01,1.0);
    fprintf(stderr,"%f\n",testdT);
    CU_ASSERT(fabs(0.119203 - testdT) <= 1e-7);

    testdT = add_layer_point(0.5,0.49,0.01,1.0);
    fprintf(stderr,"%f\n",testdT);
    CU_ASSERT(fabs(0.880797 - testdT) <= 1e-7);

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
   if (NULL == CU_add_test(pSuite, "test of idxNz()", testidxNz) ||
           NULL == CU_add_test(pSuite, "test of idxTemp()", testidxTemp) ||
           NULL == CU_add_test(pSuite, "test of read_field_order()", testReadFieldOrder) ||
           NULL == CU_add_test(pSuite, "test of calc_adiabatic_conditions()", testCalcAdiabaticConditionsConst) ||
           NULL == CU_add_test(pSuite, "test of add_layer_point()", testAddLayerPoint) ||
           NULL == CU_add_test(pSuite, "test of read_perplex_header()", testReadPerplexHeader)
           )
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


