/****************************************************************************
 * jacobiBC2.c                                                                 *
 *                                                                          *
 * Solves a system of linear equations using Jacobi iteration               *
 *                                                                          *
 *    11/26/04 WQ -- file created                                           *
 ****************************************************************************/

#include "jacobi.h"
#include <math.h>
#include <string.h> /* for memcpy */
#include <stdio.h>


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

/* update x after calculation on each ring, calculated from outside to center */
/* overwrite tmp array element by element rather than using memcpy */
int JacobiIterative(
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
  int nModelDim;
  int index; 
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

	  /* determine if we're done */
	  error = 0;

	  for (i = 1; i <= nModelDim; i++){
	    for(j = 1; j <= nModelDim; j++){
	    
	      /* get next one from calc order vector */
	      index = (i * (matrix_size) + j); 
	      x[index] = 0.25 * ( last_x[index-1] + last_x[index+1]
			      //+ last_x[index-nModelDim] + last_x[index+nModelDim]);
			      + last_x[index-matrix_size] + last_x[index+matrix_size]);

	      /* determine error before overwrite last_x */
	      if ( fabs(x[index] - last_x[index])/fabs(x[index]) > error ) {
		error = fabs(x[index] - last_x[index])/fabs(x[index]);
		//printf("%.5f %0.5f--", x[m], last_x[m]);
	      }
	      //printf("++%.3f ", x[i]);

	      /* if any entry is greater than ELEMENT_MAX, consider it an overflow
		 and abort */
	      if ((x[index] >= max_element_val) || (x[index] <= -max_element_val))
		{
		  bOverflow = TRUE;
		  break;
		}

	      last_x[index] = x[index]; 
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
	  //memcpy(last_x, x, sizeof(double) * matrix_size *  matrix_size);
	} while (!done);

      free(last_x);

      if ( bOverflow )
	iteration = -iteration;

      /* success if iteration between 0 and max_iterations*/
      return iteration;
}//JacobiIterative
