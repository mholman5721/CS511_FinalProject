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

/* solves a system of linear equations using Jacobi iteration */
int Jacobi(double *a,
                                 double * b,
                                 double * x, 
                                 const unsigned int matrix_size,
                                 const unsigned int max_iterations,
                                 const double max_element_val,
                                 const double tolerance);

int JacobiC(
                                 double * x, 
                                 const unsigned int matrix_size,
                                 const unsigned int max_iterations,
                                 const double max_element_val,
                                 const double tolerance);
int JacobiR(
                                 double * x, 
                                 const unsigned int matrix_size,
                                 const unsigned int max_iterations,
                                 const double max_element_val,
                                 const double tolerance);

int JacobiBC(
                                 double * x, 
                                 int * o,
                                 const unsigned int matrix_size,
                                 const unsigned int max_iterations,
                                 const double max_element_val,
                                 const double tolerance);

#endif /* __JACOBI_H_11262004__ */

