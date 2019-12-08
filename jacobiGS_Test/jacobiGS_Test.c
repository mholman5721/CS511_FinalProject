/****************************************************************************
 * boundaryTest.c                                                           *
 *                                                                          *
 * Tests the convergence of the Jacobi algorithm versus the Gauss-Seidel    *
 *   algorithm                                                              *
 *                                                                          *
 *    11/26/04 WQ -- file created                                           *
 ****************************************************************************/

#include "jacobi.h"
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
    
    double up;
    double down;
    double right;
    double left;
} g_myConfig;

/*
 * define constants for this file
 */
//#define MATRIX_DIM 100
//#define MATRIX_DEBUG 0
int MATRIX_DEBUG;
const char* g_csFileConfig = "config.ini";

boolean GetConfig();
void MakeInitialMatrix(double x[MATRIX_DIM+2][MATRIX_DIM+2]);
void PrintSolution(double x[MATRIX_DIM+2][MATRIX_DIM+2]);
void GetCalcOrder(int o[MATRIX_DIM*MATRIX_DIM]);
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
    double solution_vector[MATRIX_DIM+2][MATRIX_DIM+2];
    unsigned int calc_order_vector[MATRIX_DIM][MATRIX_DIM];
    unsigned int jacobi_convergences = 0, jacobi_iterations = 0;
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

    GetCalcOrder((int*)calc_order_vector);

    /* initialize the random number generator */
    srand((unsigned)time(NULL));

    g_myConfig.up = 10;
    g_myConfig.down = 10;
    g_myConfig.right = 10;
    g_myConfig.left = 10;

    /* initialize matrix, set up boundary condition */
    MakeInitialMatrix(solution_vector);

    /* test the matrix with the Jacobi algorithm */
    gettimeofday(&tvStart, &tzp); 
    result = JacobiBC((double *)solution_vector,
              (int*)calc_order_vector, MATRIX_DIM+2, 
              g_myConfig.max_iteration, 
              g_myConfig.max_element_value, 
              g_myConfig.tolerance);
    gettimeofday(&tvEnd, &tzp); 
    n_uSecLapsedJ = GetTimeDiff(tvStart, tvEnd);
    if ( MATRIX_DEBUG > 0 ) {
        printf("Jacobi Solution: \n");
        PrintSolution(solution_vector);
    }

    /* display the results */
    if ( MATRIX_DEBUG >= 0 ) {
        printf("Results:\n");
        printf("JacobiBC: Iterations=%d Time Cost=%d\n",
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

void MakeInitialMatrix(double x[MATRIX_DIM+2][MATRIX_DIM+2])
{
    unsigned int i, j;
    unsigned int modelDimension;
    modelDimension = MATRIX_DIM + 2;
    if ( MATRIX_DEBUG > 0 ) 
        printf("\nPrint matrix x[][]\n");

    g_myConfig.right = 10;
    for (i=0; i<modelDimension; i++) {
        x[0][i] = g_myConfig.up;	
        x[modelDimension-1][i] = g_myConfig.down;	
        x[i][0] = g_myConfig.left;	
        x[i][modelDimension-1] = g_myConfig.right;	
    }
    for (i=1; i<modelDimension-1; i++)
        for (j=1; j<modelDimension-1; j++)
                x[i][j] = 1.0;

    if ( MATRIX_DEBUG > 0 ) {
        for (i=0; i<modelDimension; i++) {
            for (j=0; j<modelDimension; j++)
            if ( i==0 || i==(modelDimension-1)
                || j==0 || j==(modelDimension-1))
                printf("%.2f\t", x[i][j]);
            else
                printf("%.3f\t", x[i][j]);
            printf("\n");
        }
    }
}

void PrintSolution(double x[MATRIX_DIM+2][MATRIX_DIM+2])
{
    if ( MATRIX_DEBUG < 0 ) return;
    unsigned int i, j;
    unsigned int modelDimension;
    modelDimension = MATRIX_DIM+2;

    for (i=0; i<modelDimension; i++)
    {
        for (j=0; j<modelDimension; j++) {
            if ( i==0 || i==(modelDimension-1)
                || j==0 || j==(modelDimension-1))
                printf("%.2f\t", x[i][j]);
            else
                printf("%.3f\t", x[i][j]);
        }
        printf("\n");
    }
}
void GetCalcOrder(int o[MATRIX_DIM*MATRIX_DIM])
{
    int i, radius;
    int ki, kj;
    int modelDim;
    int m, low, high;
    modelDim = MATRIX_DIM;

    radius = (modelDim + 1) / 2; 
    i = 0;
    for ( m = 0; m < radius; m++ ) 
    {
        for (ki = 0; ki < modelDim; ki++) 
            for ( kj=0; kj < modelDim; kj++)
            {
                low = m;
                high = modelDim - 1 - m;
                if(((ki==low || ki==high) && (kj >= low && kj <= high)) 
                || ((kj==low || kj==high) && (ki >= low && ki <= high))  
                  ) 
                {
                    o[i] = ki*modelDim+kj;
                    o[i] += ki * 2 + modelDim + 2 + 1;
                    i++;
                }
            }
    } 

    if ( MATRIX_DEBUG > 0 ) 
    {
        printf("\nCalcOrder{");
        for (i=0; i<MATRIX_DIM*MATRIX_DIM; i++)
            printf("%d,", o[i]);
        printf("}\n");
    }
}
