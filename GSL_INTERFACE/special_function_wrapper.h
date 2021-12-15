/* C Wrapper to GNU GSL special functions for use in
   FORTRAN programs 
   
   Created by: Mark A. Herndon
   Lehigh University, Department of Mechanical Engineering
   and Mechanics */
 
#include <gsl/gsl_sf_bessel.h>

double bessel_j0_wrapper(double *x);

double bessel_j1_wrapper(double *x);

double bessel_jn_wrapper(int *n, double *x);

double bessel_y0_wrapper(double *x);

double bessel_y1_wrapper(double *x);

double bessel_yn_wrapper(int *n, double *x);

double bessel_i0_wrapper(double *x);

double bessel_i1_wrapper(double *x);

double bessel_in_wrapper(int *n, double *x);

double bessel_k0_wrapper(double *x);

double bessel_k1_wrapper(double *x);

double bessel_kn_wrapper(int *n, double *x);
