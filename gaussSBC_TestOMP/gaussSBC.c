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
    int iteration;
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
        int nOA_size = nModelDim  * nModelDim;
        if(ID == 0) nthreads = nthrds;

        /* create a dynamic temp array */
        double * calc_x;
        calc_x = (double * )malloc(matrix_size * matrix_size * sizeof(double));

        /* initialize the last iteration matrix */
        (void)memcpy(calc_x, x, sizeof(double) * matrix_size *  matrix_size);

        do 
        {
            /* determine if we're done */
            errVal = 0;
            error = 0;

            /* initialize */
            nRingLevel = nModelDim; 
            nRingCnt = 0;

            for (m=ID; m < nOA_size; m = m + nthrds){
                /* get next one from calc order vector */
                i = o[m];
                calc_x[i] = 0.25 * (x[i-1] + x[i+1] + x[i-matrix_size] + x[i+matrix_size]);

                nRingCnt++;

                /* determine error before overwrite last_x */
                if ( fabs(calc_x[i] - last_x[i])/fabs(calc_x[i]) > errVal ) {
                    errVal = fabs(calc_x[i] - last_x[i])/fabs(calc_x[i]);
                }

                /* if any entry is greater than ELEMENT_MAX, consider it an overflow and abort */
                if ((calc_x[i] >= max_element_val) || (calc_x[i] <= -max_element_val)){
                    #pragma omp critical
                    {
                        printf("OVERFLOW! %lf\n", calc_x[i]);
                        bOverflow = TRUE;
                    }
                    //break;
                }
                
                if ( nRingCnt == 4 * (nRingLevel - 1) ) { 
                    #pragma omp critical
                    {   
                        /* update calc_x */ 
                        memcpy(x, calc_x, sizeof(double) * matrix_size *  matrix_size); 

                        /* update last_x */ 
                        memcpy(last_x, x, sizeof(double) * matrix_size *  matrix_size); 

                        /* reset nRingCnt to zero */
                        nRingCnt = 0;

                        /* go to next level */
                        nRingLevel = nRingLevel - 2;
                    }
                }
            }
            #pragma omp critical
            {
                /* update error value over all threads */
                if(errVal > error){
                    error = errVal;
                }
                
                /* increment the iteration counter */
                iteration++;

                /* update calc_x */ 
                memcpy(x, calc_x, sizeof(double) * matrix_size *  matrix_size);

                /* copy next_iteration to last_iteration */
                memcpy(last_x, x, sizeof(double) * matrix_size *  matrix_size);

                /* we are done if the iteration count exceeds the maximum number of iterations or the calculation converge */
                if (iteration > max_iterations || error < tolerance || bOverflow) {
                    done = TRUE;
                }
            }
            //printf("Thread: %d\n", ID);
        } while (!done);
    }

    free(last_x);

    printf("nthreads: %d\n", nthreads);

    /* success if iteration between 0 and max_iterations*/
    if ( bOverflow ) 
        iteration = -iteration;
        
    return iteration;
}//GaussSBC