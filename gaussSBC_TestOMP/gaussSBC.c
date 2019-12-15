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
    double errVal[NUM_THREADS];
    int nModelDim = matrix_size-2;
    int matrix_size_2 = matrix_size * matrix_size;
    int nthreads;
    int iteration;

    /* verify inputs */
    if ( matrix_size == 0 || tolerance <= 0 || max_iterations == 0 )
        return 0; /* invalid inputs */ 

    /* create a dynamic temp array */
    double * last_x;
    last_x = (double * )malloc(sizeof(double) * matrix_size_2);

    /* initialize the last iteration matrix */
    (void)memcpy(last_x, x, sizeof(double) * matrix_size_2);

    /* ring is determined by following method      */
    /* outmost ring ... 4 * ( nModelDim - 1 )      */
    /* next ring ...  4 * ( nModelDim - 2 - 1)     */
    /* until nModelDim - 2*i == 2 or 1(inner most) */
    /* update last_x[i] when ring level is changed */
    #pragma omp parallel num_threads(NUM_THREADS) 
    {
        boolean done, bOverflow; 
        int m, i;
        int nRingCnt = 0;
        int nRingLevel = 0;
        int ID = omp_get_thread_num();
        int nthrds = omp_get_num_threads();
        int nOA_size = nModelDim  * nModelDim;
        if(ID == 0) nthreads = nthrds;

        for(int b = 0; b < NUM_THREADS; b++){
            errVal[b] = 0.0;
        }

        /* calculate iterations */
        done = FALSE;
        bOverflow = FALSE;

        /* create a dynamic temp array */
        double * calc_x;
        calc_x = (double * )malloc(sizeof(double) * matrix_size_2);

        /* initialize the last iteration matrix */
        (void)memcpy(calc_x, x, sizeof(double) * matrix_size_2);

        /* create a dynamic temp array */
        double * calcLast_x;
        calcLast_x = (double * )malloc(sizeof(double) * matrix_size_2);

        /* initialize the last iteration matrix */
        (void)memcpy(calcLast_x, x, sizeof(double) * matrix_size_2);

        //int chunkSize = nOA_size / NUM_THREADS;

        #pragma omp barrier

        do 
        {
            /* determine if we're done */
            errVal[ID] = 0.0;

            /* initialize */
            nRingLevel = nModelDim; 
            nRingCnt = 0;

            for (m = ID; m < nOA_size; m = m + nthrds){
            //for (m = (ID * chunkSize); m < (ID * chunkSize) + chunkSize; m++){
                /* get next one from calc order vector */
                i = o[m];
                calc_x[i] = 0.25 * (calc_x[i-1] + calc_x[i+1] + calc_x[i-matrix_size] + calc_x[i+matrix_size]);

                nRingCnt++;

                /* determine error before overwrite last_x */
                if ( fabs(calc_x[i] - calcLast_x[i])/fabs(calc_x[i]) > errVal[ID] ) {
                    errVal[ID] = fabs(calc_x[i] - calcLast_x[i])/fabs(calc_x[i]);
                }

                /* if any entry is greater than ELEMENT_MAX, consider it an overflow and abort */
                if ((calc_x[i] >= max_element_val) || (calc_x[i] <= -max_element_val)){
                    #pragma omp critical
                    {
                        printf("OVERFLOW! %lf\n", calc_x[i]);
                        bOverflow = TRUE;
                    }
                    break;
                }

                if ( (nRingCnt * NUM_THREADS) == 4 * (nRingLevel - 1) ) { 
                    #pragma omp critical
                    {   
                        /* update calc_x */ 
                        for(int n = 0; n < matrix_size_2; n++){
                            if(x[n] < calc_x[n]) {
                                x[n] = calc_x[n];
                            }
                        }

                        /* update calcLast_x */ 
                        for(int n = 0; n < matrix_size_2; n++){
                            if(calcLast_x[n] < x[n]) {
                                calcLast_x[n] = x[n];
                            }
                        }

                        /* update last_x */ 
                        for(int n = 0; n < matrix_size_2; n++){
                            if(last_x[n] < x[n]) {
                                last_x[n] = x[n];
                            }
                        }

                        /* reset nRingCnt to zero */
                        nRingCnt = 0;

                        /* go to next level */
                        nRingLevel = nRingLevel - 2;
                    }
                } else if (nModelDim % 2 != 0 && i == o[nOA_size-1]){
                    printf("HERE\n");
                    #pragma omp critical
                    {   
                        /* update calc_x */ 
                        for(int n = 0; n < matrix_size_2; n++){
                            if(x[n] < calc_x[n]) {
                                x[n] = calc_x[n];
                            }
                        }

                        /* update calcLast_x */ 
                        for(int n = 0; n < matrix_size_2; n++){
                            if(calcLast_x[n] < x[n]) {
                                calcLast_x[n] = x[n];
                            }
                        }

                        /* update last_x */ 
                        for(int n = 0; n < matrix_size_2; n++){
                            if(last_x[n] < x[n]) {
                                last_x[n] = x[n];
                            }
                        }

                        /* reset nRingCnt to zero */
                        nRingCnt = 0;

                        /* go to next level */
                        nRingLevel = nRingLevel - 2;
                    }
                }
            }
            
            #pragma omp critical
            {                
                /* increment the iteration counter */
                iteration++;

                /* update calc_x */ 
                for(int n = 0; n < matrix_size_2; n++){
                    if(x[n] < calc_x[n]) {
                        x[n] = calc_x[n];
                    }
                }

                /* update calcLast_x */ 
                for(int n = 0; n < matrix_size_2; n++){
                    if(calcLast_x[n] < x[n]) {
                        calcLast_x[n] = x[n];
                    }
                }

                /* update last_x */ 
                for(int n = 0; n < matrix_size_2; n++){
                    if(last_x[n] < x[n]) {
                        last_x[n] = x[n];
                    }
                }

                int errCount = 0;
                for(int n = 0; n < NUM_THREADS; n++){
                    if(errVal[n] < tolerance){
                        errCount++;
                    }
                }

                /* we are done if the iteration count exceeds the maximum number of iterations or the calculation converge */
                if (iteration > max_iterations){ 
                    printf("ITERATIONS: %d / %d\n", iteration, max_iterations);
                    done = TRUE;
                } else if (bOverflow) {
                    printf("OVERFLOW: TRUE\n");
                    iteration = -iteration;
                    done = TRUE;
                } else if (errCount == NUM_THREADS){
                    printf("ERROR [ %d ]: %lf / %lf\n", ID, errVal, tolerance);
                    done = TRUE;
                } 
            }
            //printf("Thread: %d\n", ID);
        } while (!done);
    }

    free(last_x);

    printf("nthreads: %d\n", nthreads);

    /* success if iteration between 0 and max_iterations*/
    //if ( bOverflow ) 
    //    iteration = -iteration;

    return iteration;
}//GaussSBC

//printf("errVal [ %d ] = %lf\n", ID, errVal);

/*for(int n = 0; n < matrix_size_2; n++){
    if(calc_x[n] > x[n]) {
        x[n] = calc_x[n];
    }
}*/

//memcpy(last_x, x, sizeof(double) * matrix_size_2);