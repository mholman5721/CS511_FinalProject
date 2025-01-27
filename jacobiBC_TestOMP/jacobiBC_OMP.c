/****************************************************************************
 * jacobiBC_OMP.c                                                           *
 *                                                                          *
 * A parallel version of the ringed-boundary jacobi using OpenMP            *
 *                                                                          *
 *    11/26/04 WQ -- file created                                           *
 ****************************************************************************/

#include "jacobi.h"
#include <math.h>
#include <string.h> /* for memcpy */
#include <stdio.h>
#include <omp.h>

/*  
 * Jacobi()
 * Solves a system of linear equations in the form a*x=b using Jacobi iteration
 *
 * Parameters:
 *   x, o - x is a matrix containing the boundary condition information and the relaxation mesh (solutions)
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

int new_min(int num1, int num2) { return (num1 > num2 ) ? num2 : num1; }

int JacobiBC_OMP(
	      double * x,
	      int * o,
	      const unsigned int matrix_size,
	      const unsigned int max_iterations,
	      const double max_element_val,
	      const double tolerance)
{
  int iteration = 0;
  int i, j, k;
  double error[NUM_THREADS];
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


  // calculate the number of rings
  int ring_num = 0;
  //create a variable to simulate the countdown
  //of rings like regular JacobiBC to evaluate
  //how many rings we have 
  int temp = nModelDim;
  
  while(temp > 0){
    ring_num++;
    temp -= 2;
  }

  //create an array and store the amount of work to be done
  //in each ring. (how many values in each)
  
  int ring_work[ring_num];
  temp = nModelDim; 
  
  for(i = 0; i < ring_num; i++){  
    //small fix for the last value if matrix size is odd
    if(nModelDim % 2 != 0 && temp == 1)
      ring_work[i] = 1;
    else
      ring_work[i] = 4 * (temp - 1); 
  
    temp -= 2;
  
  }

  int starting_index = 0; 

  
#pragma omp parallel num_threads(NUM_THREADS) shared(last_x, starting_index, done)
  {
    do
      {
	
	/* determine if we're done */
	error[omp_get_thread_num()] = 0;
	int my_rank = omp_get_thread_num();
	
	for (j = 0; j < ring_num;)
	  {
	    
	    // we have to determine where we begin inside the total array
	    // we do this by summing all the work done previously (if any)
	    
            #pragma omp single
	    {
	      if(j == 0)
		starting_index = 0;
	      else{
		starting_index += ring_work[j - 1];
	      }
	    }
	    
	    #pragma omp barrier

	    if(done == true){
	      break;
	    }
	   
	    int THREAD_COUNT = omp_get_num_threads();
	    int my_first_index, my_last_index;
	    bool my_update = true; 
	    // if there is more work than threads (read indices) have some threads just
	    // skip the whole thing
	    
	    if(ring_work[j] - 1 >= my_rank){
	      
	      // calculate the starting and ending indices on the given ring
	      // using the amount of work, the rank, and the starting index
	      
	      // make sure the amount of working threads (which we use to divide the work)
	      // is equal to or less than our thread count (so if there's 2 indices we want two
	      // threads to do the work not four and if there's 10 indices we want four not ten
	      int working_threads = new_min(ring_work[j], THREAD_COUNT);
	      int local_m = ceil(ring_work[j] / working_threads);
	      my_first_index = (my_rank * local_m) + starting_index;
	      // take the minimum between our assigned work and the end index of our current ring
	      // we do this because due to the threads that might divide the work unevenly our ceil
	      // may assign work beyond our current ring
	      my_last_index = new_min((((my_rank + 1) * local_m) + starting_index), (starting_index + ring_work[j]));   

	      for(int m = my_first_index; m < my_last_index; m++){
		
		/* get next one from calc order vector */
		i = o[m];
		x[i] = 0.25 * ( last_x[i-1] + last_x[i+1]
				+ last_x[i-matrix_size] + last_x[i+matrix_size]);
		
		/* determine error before overwrite last_x */
		if ( fabs(x[i] - last_x[i])/fabs(x[i]) > error[my_rank] )
		  error[my_rank] = fabs(x[i] - last_x[i])/fabs(x[i]);
		
		/* if any entry is greater than ELEMENT_MAX, consider it an overflow
		   and abort */
		if ((x[i] >= max_element_val) || (x[i] <= -max_element_val))
		  {
		    bOverflow = TRUE;
		    my_update = false;
		    done = true;
		    break;
		  }
	      }
	    }
	    
	    #pragma omp barrier 
	    
	    //update the indices if they were a part of the work. we can't just wrap the if statement
	    //around all the work because then the threads who may not be doing work won't get to access it and we hang
	    if(ring_work[j] >= my_rank && my_update == true){
	      for(int m = my_first_index; m < my_last_index; m++)
		last_x[o[m]] = x[o[m]];
	    }
	  
       
	    if(done == true)
	       break;

	    #pragma omp single
	    j++;
	     
	  }
	

	/* increment the iteration counter but only once each time */
	#pragma omp single
	iteration++;
	
	/* we are done if the iteration count exceeds the maximum number of iterations
	   or the calculation converge */
	if ( iteration > max_iterations
	     || error[my_rank] < tolerance
	     || bOverflow)
	  {
	    done = TRUE;
	  }

	//printf("at end\n"); 

      #pragma omp barrier 
	
      } while (!done);
  }
  free(last_x);
  
  if ( bOverflow )
    iteration = -iteration;
  
  /* success if iteration between 0 and max_iterations*/
  return iteration;
}//JacobiIterative
