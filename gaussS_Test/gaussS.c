/****************************************************************************
 * jacobi.c                                                                 *
 *                                                                          *
 * Solves a system of linear equations using Jacobi iteration               *
 *                                                                          *
 *    11/26/04 WQ -- file created                                           *
 ****************************************************************************/

#include "gaussS.h"
#include <math.h>
#include <string.h> /* for memcpy */
#include <stdio.h>
#include <stdlib.h>

/*  
 * Jacobi()
 * Solves a system of linear equations in the form a*x=b using Jacobi iteration
 *
 * Parameters:
 *   a, b - matrices containing the system of linear equations
 *   x - destination matrix used to store the solution to the system
 *   matrix_size - matrix size, a would be a m * m two dimension array while b and x are one 
 *       dimension array
 *   max_iterations
 *   tolerance -  
 * Return Value:
 *   If successful, the number of iterations required to reach the requested
 *     precision; 0 if the solution does not converge in limited iterations
 */
int Jacobi(double * a,
                                 double * b,
                                 double * x, 
                                 const unsigned int matrix_size, 
                                 const unsigned int max_iterations,
                                 const double max_element_val,
                                 const double tolerance)
{
    int iteration = 0;
    unsigned int i, j;
    double sigma;
    double error;
    boolean   done, bOverflow;

    /* verify inputs */
    if ( matrix_size == 0 || tolerance <= 0 || max_iterations == 0 )
        return 0; /* invalid inputs */ 

    /* create a dynamic temp array */
    double * last_x;
    last_x = (double * )malloc(matrix_size * sizeof(double));

    /* initialize the last iteration matrix (x[n-1]) */
    for (i=0; i<matrix_size; i++)
        last_x[i] = 1.0; /* use 1 instead of 0 to avoid multiplying
                                    by 0 */

    /* calculate iterations */
    done = FALSE;
    bOverflow = FALSE;
    do 
    {
        /* calculate the next iteration */
        for (i=0; i<matrix_size; i++)
        {
            sigma = 0.0;
            for (j=0; j<matrix_size; j++)
            {
                if (j != i)
                    sigma += a[i*matrix_size+j] * last_x[j];
            }
            x[i] = (b[i] - sigma) / a[i*matrix_size+i];
        }

        /* determine if we're done */
        error = 0;
        for (i=0; i<matrix_size; i++)
        {
            if ( fabs(x[i] - last_x[i])/fabs(x[i]) > error ) 
                error = fabs(x[i] - last_x[i])/fabs(x[i]);
        }

        /* increment the iteration counter */
        iteration++;

        /* if any entry is greater than ELEMENT_MAX, consider it an overflow
           and abort */
        for (i=0; i<matrix_size; i++)
            if ((x[i] >= max_element_val) || (x[i] <= -max_element_val))
            {
                bOverflow = TRUE;
                break;
            }

        /* we are done if the iteration count exceeds the maximum number of iterations
           or the calculation converge */
        if ( iteration > max_iterations 
          || error < tolerance
          || bOverflow) 
        {
            done = TRUE;
        }

        /* copy next_iteration to last_iteration */
        (void)memcpy(last_x, x, sizeof(double) * matrix_size);
    } while (!done);

    free(last_x);
  
    if ( bOverflow ) 
        iteration = -iteration;

    /* success if iteration between 0 and max_iterations*/
    return iteration;
}//jacobi

/* update x after calculation on each ring, calculated from outside to center */
int JacobiBC(
                                 double * x, 
                                 int * o,
                                 const unsigned int matrix_size, 
                                 const unsigned int max_iterations,
                                 const double max_element_val,
                                 const double tolerance)
{
    int iteration = 0;
    unsigned int i, j;
    double error;
    boolean   done, bOverflow;
    int nModelDim, nRingLevel, nRingCnt;
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
    do 
    {
        /* initialize */
        nRingLevel = nModelDim; 
        nRingCnt = 0;

        /* calculate the next iteration */
        int m;
        /* determine if we're done */
        error = 0;

        int nOA_size = nModelDim  * nModelDim;
        for (m=0; m<nOA_size; m++)
        {
            /* get next one from calc order vector */
            i = o[m];
            x[i] = 0.25 * ( last_x[i-1] + last_x[i+1]
                //+ last_x[i-nModelDim] + last_x[i+nModelDim]);
                + last_x[i-matrix_size] + last_x[i+matrix_size]);

            nRingCnt++;
            //printf("%d--%d  ", nRingLevel, i);

            /* determine error before overwrite last_x */
            if ( fabs(x[i] - last_x[i])/fabs(x[i]) > error ) {
                error = fabs(x[i] - last_x[i])/fabs(x[i]);
                //printf("%.5f %0.5f--", x[m], last_x[m]);
            }
            //printf("++%.3f ", x[i]);

            /* if any entry is greater than ELEMENT_MAX, consider it an overflow
               and abort */
            if ((x[i] >= max_element_val) || (x[i] <= -max_element_val))
            {
                bOverflow = TRUE;
                break;
            }

            if ( nRingCnt == 4 * (nRingLevel - 1) ) 
            { 
                /* update last_x */ 
    		memcpy(last_x, x, sizeof(double) * matrix_size *  matrix_size);
                
                /* reset nRingCnt to zero */
                nRingCnt = 0;

                /* go to next level */
                nRingLevel = nRingLevel - 2;
                //printf("\n");
            }
        }

        /* increment the iteration counter */
        iteration++;

        /* we are done if the iteration count exceeds the maximum number of iterations
           or the calculation converge */
        if ( iteration > max_iterations 
          || error < tolerance
          || bOverflow) 
        {
            done = TRUE;
        }

        /* copy next_iteration to last_iteration */
    	memcpy(last_x, x, sizeof(double) * matrix_size *  matrix_size);
    } while (!done);

    free(last_x);
  
    if ( bOverflow ) 
        iteration = -iteration;

    /* success if iteration between 0 and max_iterations*/
    return iteration;
}//JacobiBC