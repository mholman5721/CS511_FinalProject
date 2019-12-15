/****************************************************************************
 * jacobiBC.c                                                               *
 *                                                                          *
 * Boundary condition - ringed jacobi with no parallelization               *
 *                                                                          *
 *    11/26/04 WQ -- file created                                           *
 ****************************************************************************/

#include "jacobi.h"
#include <math.h>
#include <string.h> /* for memcpy */
#include <stdio.h>


/*  
 * JacobiBC()
 * Given a relaxation mesh. Iterate until the boundary has been distributed through the mesh
 *
 * Parameters:
 *   x, o - x contains the boundary and relaxation mesh (solution) values
 *   matrix_size - matrix size, a would be a m * m two dimension array while b and x are one 
 *       dimension array
 *   max_iterations
 *   tolerance -  
 * Return Value:
 *   If successful, the number of iterations required to reach the requested
 *     precision; 0 if the solution does not converge in limited iterations
 */

/* update x after calculation on each ring, calculated from outside to center */
/* overwrite tmp array element by element rather than using memcpy */
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
	  int m,mm;
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

		  //memcpy(last_x, x, sizeof(double) * matrix_size *  matrix_size);
		  for (mm=1; mm<=nRingCnt; mm++) {
		    last_x[o[m-nRingCnt+mm]] = x[o[m-nRingCnt+mm]]; 
		  }

		  /* reset nRingCnt to zero */
		  nRingCnt = 0;

		  /* go to next level */
		  nRingLevel = nRingLevel - 2;
		  //printf("\n");
		}
	      else if ( nModelDim % 2 != 0 && m == nOA_size - 1)
		last_x[o[m]] = x[o[m]]; 

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
	  //memcpy(last_x, x, sizeof(double) * matrix_size *  matrix_size);
	} while (!done);

      free(last_x);

      if ( bOverflow )
	iteration = -iteration;

      /* success if iteration between 0 and max_iterations*/
      return iteration;
}//JacobiIterative
