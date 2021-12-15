/* C Wrapper to GNU GSL special functions for use in
   FORTRAN programs 
   
   Created by: Mark A. Herndon
   Lehigh University, Department of Mechanical Engineering
   and Mechanics */
 
#include <gsl/gsl_sf_bessel.h>
#include "special_function_wrapper.h"

double bessel_j0_wrapper(double *x){
	return gsl_sf_bessel_J0(*x);
};

double bessel_j1_wrapper(double *x){
	return gsl_sf_bessel_J1(*x);
};

double bessel_jn_wrapper(int *n, double *x){
	return gsl_sf_bessel_Jn(*n, *x);
};

double bessel_y0_wrapper(double *x){
	return gsl_sf_bessel_Y0(*x);
};

double bessel_y1_wrapper(double *x){
	return gsl_sf_bessel_Y1(*x);
};

double bessel_yn_wrapper(int *n, double *x){
	return gsl_sf_bessel_Yn(*n, *x);
};

double bessel_i0_wrapper(double *x){
	return gsl_sf_bessel_I0(*x);
};

double bessel_i1_wrapper(double *x){
	return gsl_sf_bessel_I1(*x);
};

double bessel_in_wrapper(int *n, double *x){
	return gsl_sf_bessel_In(*n, *x);
};
double bessel_k0_wrapper(double *x){
	return gsl_sf_bessel_K0(*x);
};

double bessel_k1_wrapper(double *x){
	return gsl_sf_bessel_K1(*x);
};

double bessel_kn_wrapper(int *n, double *x){
	return gsl_sf_bessel_Kn(*n, *x);
};
