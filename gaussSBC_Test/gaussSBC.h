/****************************************************************************
 * jacobi.h                                                                 *
 *                                                                          *
 * Header file for jacobi.c                                                 *
 *                                                                          *
 *    11/26/04 WQ -- file created                                           *
 ****************************************************************************/

#ifndef __JACOBI_H_11262004__
#define __JACOBI_H_11262004__

#include "const.h"

int GaussSBC(
                                 double * x, 
                                 int * o,
                                 const unsigned int matrix_size,
                                 const unsigned int max_iterations,
                                 const double max_element_val,
                                 const double tolerance);

#endif /* __JACOBI_H_11262004__ */

