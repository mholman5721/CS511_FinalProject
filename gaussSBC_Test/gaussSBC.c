/****************************************************************************
 * jacobi.c                                                                 *
 *                                                                          *
 * Solves a system of linear equations using Jacobi iteration               *
 *                                                                          *
 *    11/26/04 WQ -- file created                                           *
 ****************************************************************************/

#include "gaussSBC.h"
#include <math.h>
#include <string.h> /* for memcpy */
#include <stdio.h>
#include <stdlib.h>
#include <omp.h>

/* update x after calculation on each ring, calculated from outside to center */
int GaussSBC(
                                 double * x, 
                                 int * o,
                                 const unsigned int matrix_size, 
                                 const unsigned int max_iterations,
                                 const double max_element_val,
                                 const double tolerance)
{
    double error = 0.0;
    int nModelDim;
    int nthreads;
    int iteration[NUM_THREADS];
    for(int i = 0; i < NUM_THREADS; i++){
        iteration[i] = 0;
    }
    boolean done, bOverflow;
    
    nModelDim = matrix_size-2;

    /* verify inputs */
    if ( matrix_size == 0 || tolerance <= 0 || max_iterations == 0 )
        return 0; /* invalid inputs */ 

    /* create a dynamic temp array */
    double * last_x;
    last_x = (double * )malloc(matrix_size * matrix_size * sizeof(double));

    /* initialize the last iteration matrix */
    (void)memcpy(last_x, x, sizeof(double) * matrix_size *  matrix_size);

    /* calculate iterations */
    done = FALSE;
    bOverflow = FALSE;

    /* ring is determined by following method      */
    /* outmost ring ... 4 * ( nModelDim - 1 )      */
    /* next ring ...  4 * ( nModelDim - 2 - 1)     */
    /* until nModelDim - 2*i == 2 or 1(inner most) */
    /* update last_x[i] when ring level is changed */
    #pragma omp parallel num_threads(NUM_THREADS) 
    {
        double errVal = 0.0;
        int m, i;
        int nRingCnt = 0;
        int nRingLevel = 0;
        int ID = omp_get_thread_num();
        int nthrds = omp_get_num_threads();
        if(ID == 0) nthreads = nthrds;

        do 
        {

            int nOA_size = nModelDim  * nModelDim;

            /* determine if we're done */
            errVal = 0;

            /* initialize */
            nRingLevel = nModelDim; 
            nRingCnt = 0;

            #pragma omp for schedule(auto)
            for (m=0; m<nOA_size; m++){
                /* get next one from calc order vector */
                i = o[m];
                x[i] = 0.25 * (x[i-1] + x[i+1] + x[i-matrix_size] + x[i+matrix_size]);

                nRingCnt++;

                /* determine error before overwrite last_x */
                if ( fabs(x[i] - last_x[i])/fabs(x[i]) > errVal ) {
                    errVal = fabs(x[i] - last_x[i])/fabs(x[i]);
                }

                /* if any entry is greater than ELEMENT_MAX, consider it an overflow and abort */
                if ((x[i] >= max_element_val) || (x[i] <= -max_element_val)){
                    printf("OVERFLOW! %lf\n", x[i]);
                    bOverflow = TRUE;
                    //break;
                }

                if ( nRingCnt == 4 * (nRingLevel - 1) ) { 
                    #pragma omp critical
                    {
                        /* update last_x */ 
                        memcpy(last_x, x, sizeof(double) * matrix_size *  matrix_size);
                    }   
                    /* reset nRingCnt to zero */
                    nRingCnt = 0;

                    /* go to next level */
                    nRingLevel = nRingLevel - 2;
                }
            }
            #pragma omp critical
            {
                error += errVal;
                /* increment the iteration counter */
                iteration[ID]++;

                /* copy next_iteration to last_iteration */
                memcpy(last_x, x, sizeof(double) * matrix_size *  matrix_size);
            }
            
            /* we are done if the iteration count exceeds the maximum number of iterations or the calculation converge */
            if (iteration[ID] > max_iterations || errVal < tolerance || bOverflow) {
                done = TRUE;
            }
        } while (!done);
    }

    free(last_x);
  
    if ( bOverflow ) 
        iteration[0] = -iteration[0];

    printf("nthreads: %d\n", nthreads);

    /* success if iteration between 0 and max_iterations*/
    int totalIterations = 0;
    for(int i = 0; i < NUM_THREADS; i++){
        totalIterations += iteration[i];
    }
    return totalIterations;
}//GaussSBC