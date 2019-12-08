/****************************************************************************
 * boundaryTest.c                                                           *
 *                                                                          *
 * Tests the convergence of the Jacobi algorithm versus the Gauss-Seidel    *
 *   algorithm                                                              *
 *                                                                          *
 *    11/26/04 WQ -- file created                                           *
 ****************************************************************************/

#include "gaussS.h"
#include "const.h"

#include <stdio.h> //for printf
#include <stdlib.h>
#include <time.h>
#include <math.h>
#include <string.h>
#include <sys/time.h> //for gettimeofday

/*
 * define data structures
 */
struct Config 
{
    double tolerance;
    unsigned int max_iteration;
    double max_element_value;
    unsigned int max_testing;
    
} g_myConfig;

/*
 * define constants for this file
 */
//#define MATRIX_DIM 100
//#define MATRIX_DEBUG 0
int MATRIX_DEBUG;
const char* g_csFileConfig = "config.ini";

boolean GetConfig();
void MakeInitialMatrix_a(double a[MATRIX_DIM][MATRIX_DIM], double val);
void MakeInitialMatrix_b(double b[MATRIX_DIM], double val);
void MakeInitialMatrix_x(double x[MATRIX_DIM], double val);
void PrintSolution(double x[MATRIX_DIM]);
int  GetTimeDiff(struct timeval tvStart, struct timeval tvEnd) {
    if (tvStart.tv_usec > tvEnd.tv_usec) {
        tvEnd.tv_usec += 1000000;
        tvEnd.tv_sec--;
    }
    return ((tvEnd.tv_sec - tvStart.tv_sec) * 1000000 
         + (tvEnd.tv_usec - tvStart.tv_usec));
}

int main(int argc, char* argv[])
{
    double a[MATRIX_DIM][MATRIX_DIM];
    double x[MATRIX_DIM];
    double b[MATRIX_DIM];

    unsigned int result;

    struct timeval tvStart, tvEnd;
    struct timezone tzp;
    
    int n_uSecLapsedJ ;

    MATRIX_DEBUG = 0;
    if ( argc == 2 ) {
        MATRIX_DEBUG = atoi(argv[1]); 
    }

    /* get configuration from config file */
    if ( GetConfig() == FALSE ) 
        printf("Failed to get settings from config file, use default"
               " values for calculation.\n"); 

    /* display welcome message */
    if ( MATRIX_DEBUG > 0 ) {
        printf("Current settings: \n"
                 "Matrix Dimension=%d\n"
                 "TOLERANCE=%f\n"
                 "MAX_ITERATION=%d\n"
                 "MAX_TESTING=%d\n"
                 "OVERFLOW=%f\n", 
                 MATRIX_DIM,
                 g_myConfig.tolerance, 
                 g_myConfig.max_iteration,
                 g_myConfig.max_testing,
                 g_myConfig.max_element_value
                );
    }

    /* initialize the random number generator */
    srand((unsigned)time(NULL));

    /* initialize matrix, set up boundary condition */
    MakeInitialMatrix_a(a, 2.0);
    MakeInitialMatrix_b(b, 2.0);
    MakeInitialMatrix_x(x, 1.0);

    /* test the matrix with the Jacobi algorithm */
    gettimeofday(&tvStart, &tzp); 
    result = gaussS((double *)a,
              (double *)b,
              (double *)x,
              MATRIX_DIM, 
              g_myConfig.max_iteration, 
              g_myConfig.max_element_value, 
              g_myConfig.tolerance);
    gettimeofday(&tvEnd, &tzp); 
    n_uSecLapsedJ = GetTimeDiff(tvStart, tvEnd);
    if ( MATRIX_DEBUG > 0 ) {
        printf("Gauss Seidel Solution: \n");
        PrintSolution(x);
    }

    /* display the results */
    if ( MATRIX_DEBUG >= 0 ) {
        printf("Results:\n");
        printf("Gauss Seidel: Iterations=%d Time Cost=%d\n",
            result, n_uSecLapsedJ);
    }

    return EXIT_SUCCESS;
}

boolean GetConfig()
{
    FILE * input_file;
    char sTag[256];
    float fValue;

    /* initialize */
    g_myConfig.tolerance = 0.001;
    g_myConfig.max_iteration = 1000;
    g_myConfig.max_element_value= 10000;
    g_myConfig.max_testing = 10;

    if ((input_file = fopen(g_csFileConfig, "r")) == NULL) 
    {
        fprintf(stderr, "Cannot open %s\n", "output_file");
        return FALSE;
    }

    while ( fscanf(input_file, "%s%f", sTag, &fValue) != EOF )
    {
        if ( strcmp(sTag, "TOLERANCE") == 0 )
            g_myConfig.tolerance = fValue;
        else if ( strcmp(sTag, "MAX_ITERATION") == 0 )
            g_myConfig.max_iteration = (int)fValue;
        else if ( strcmp(sTag, "MAX_ELEMENT_VALUE") == 0 )
            g_myConfig.max_element_value= fValue;
        else if ( strcmp(sTag, "MAX_TESTING") == 0 )
            g_myConfig.max_testing = (int)fValue;
    }
    fclose(input_file);
    return TRUE;
}

void MakeInitialMatrix_a(double a[MATRIX_DIM][MATRIX_DIM], double val)
{
    if ( MATRIX_DEBUG > 0 ) 
        printf("\nPrint matrix a[][]\n");

    for (int i = 0; i < MATRIX_DIM; i++)
        for (int j = 0; j < MATRIX_DIM; j++)
                a[i][j] = val;

    if ( MATRIX_DEBUG > 0 ) {
        for (int i = 0; i < MATRIX_DIM; i++) {
            for (int j = 0; j < MATRIX_DIM; j++) {
                printf("%.2f\t", a[i][j]);
                
            }
            printf("\n");
        }
    }
}

void MakeInitialMatrix_b(double b[MATRIX_DIM], double val)
{
    if ( MATRIX_DEBUG > 0 ) 
        printf("\nPrint matrix b[]\n");

    for (int i = 0; i < MATRIX_DIM; i++)
        b[i] = val;

    if ( MATRIX_DEBUG > 0 ) {
        for (int i = 0; i < MATRIX_DIM; i++) {
            printf("%.2f\t", b[i]);
        }
        printf("\n");
    }
}

void MakeInitialMatrix_x(double x[MATRIX_DIM], double val)
{
    if ( MATRIX_DEBUG > 0 ) 
        printf("\nPrint matrix x[]\n");

    for (int i = 0; i < MATRIX_DIM; i++)
        x[i] = val;

    if ( MATRIX_DEBUG > 0 ) {
        for (int i = 0; i < MATRIX_DIM; i++) {
            printf("%.2f\t", x[i]);
        }
        printf("\n");
    }
}

void PrintSolution(double x[MATRIX_DIM])
{
    if ( MATRIX_DEBUG < 0 ) return;

    for (int i = 0; i<MATRIX_DIM; i++){
        printf("%.3f\t", x[i]);
    }
    printf("\n");
}