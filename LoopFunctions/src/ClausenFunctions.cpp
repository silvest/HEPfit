/* 
 * File:   ClausenFunctions.cpp
 * Author: mishima
 */

#include <cstdlib>
#include <cmath>
#include <gsl/gsl_monte.h>
#include <gsl/gsl_monte_vegas.h>
#include <gsl/gsl_sf.h>
#include "ClausenFunctions.h"


ClausenFunctions::ClausenFunctions() {
}

//ClausenFunctions::ClausenFunctions(const ClausenFunctions& orig) {
//}

ClausenFunctions::~ClausenFunctions() {
}


////////////////////////////////////////////////////////////////////////

double ClausenFunctions::Cl2(const double phi) const {
    return ( gsl_sf_clausen(phi) );    
}

double ClausenFunctions::Cl3(const double phi) const {
    struct ClausenFunctions::my_f_params params = { phi }; 
    
    // Randam number generator
    const gsl_rng_type *T;
    gsl_rng *r;
    gsl_rng_env_setup();
    T = gsl_rng_default;
    r = gsl_rng_alloc(T);        
    // gsl_momte_vegas  
    size_t calls = 50000, itmx = 5, dim = 1;
    double xl[1] = {0}, xu[1] = {1};
    double res, err;
    gsl_monte_function func = {ClausenFunctions::integrand_for_Li3_imag, dim, &params};
    gsl_monte_vegas_state *st = gsl_monte_vegas_alloc(dim);
    st->iterations = itmx;
    gsl_monte_vegas_integrate(&func, xl, xu, dim, calls/50, r, st, &res, &err);
    do {
        gsl_monte_vegas_integrate(&func, xl, xu, dim, calls/5, r, st, &res, &err);
    } while (fabs(st->chisq - 1.0) > 0.5);	
    gsl_monte_vegas_free(st);
    gsl_rng_free(r);        

    return res;    
}

double ClausenFunctions::integrand_for_Li3_imag(double *k, size_t dim, void *params) {
    struct my_f_params* fp = (struct my_f_params*)params;
    double phi = fp->phi;
    if (1.0 - 2.0*k[0]*cos(phi) + k[0]*k[0]==0.0 || k[0]==0.0) 
        throw "Error in ClausenFunctions::integrand_for_Li3_imag()";
    double log_k0 = log(k[0]);
    return ( log_k0*log_k0*sin(phi)/(1.0 - 2.0*k[0]*cos(phi) + k[0]*k[0])/2.0 );
}



