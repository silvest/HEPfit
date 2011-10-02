/* 
 * File:   Polylogarithms.cpp
 * Author: mishima
 */

#include <cstdlib>
#include <cmath>
#include <gsl/gsl_monte.h>
#include <gsl/gsl_monte_vegas.h>
#include <gsl/gsl_sf.h>
#include "Polylogarithms.h"


Polylogarithms::Polylogarithms() {
    B0[0] = 1.0; B0[1] = -1.0/2.0; B0[2] = 1.0/6.0; B0[3] = 0.0; 
    B0[4] = -1.0/30.0; B0[5] = 0.0; B0[6] = 1.0/42.0; B0[7] = 0.0;
    B0[8] = -1.0/30.0; B0[9] = 0.0; B0[10] = 5.0/66.0; B0[11] = 0.0; 
    B0[12] = -691.0/2730.0; B0[13] = 0.0; B0[14] = 7.0/6.0; B0[15] = 0.0;
    B0[16] = -3617.0/510.0; B0[17] = 0.0; B0[18] = 43867.0/798.0;
}

//Polylogarithms::Polylogarithms(const Polylogarithms& orig) {
//}

Polylogarithms::~Polylogarithms() {
}


////////////////////////////////////////////////////////////////////////

double Polylogarithms::integrand_for_Li3(double *k, size_t dim, void *params) {
        return ( gsl_sf_dilog(k[0])/k[0] );
}

double Polylogarithms::Li3(const double x) const {
    double Li3;
    if (x < 0.0) {
        Li3 = -gsl_sf_fermi_dirac_2(log(-x));
    } else if (x == 0.0) {
        Li3 = 0.0;    
    } else if (x > 0.0 && x < 0.5) {
        // Randam number generator
        const gsl_rng_type *T;
        gsl_rng *r;
        gsl_rng_env_setup();
        T = gsl_rng_default;
        r = gsl_rng_alloc(T);        
        // gsl_momte_vegas  
        size_t calls = 50000, itmx = 5, dim = 1;
        double xl[1] = {0}, xu[1] = {x};
        double res, err;
        gsl_monte_function func = {Polylogarithms::integrand_for_Li3, dim, 0};
        gsl_monte_vegas_state *st = gsl_monte_vegas_alloc(dim);
        st->iterations = itmx;
        gsl_monte_vegas_integrate(&func, xl, xu, dim, calls/50, r, st, &res, &err);
        do {
            gsl_monte_vegas_integrate(&func, xl, xu, dim, calls/5, r, st, &res, &err);
        } while (fabs(st->chisq - 1.0) > 0.5);	
        gsl_monte_vegas_free(st);
        gsl_rng_free(r);        
        //
        Li3 = res;
    } else if (x == 0.5) {
        double log2 = log(2.0);
        double zeta3 = gsl_sf_zeta_int(3);
        Li3 = (4.0*pow(log2, 3.0) - 2.0*M_PI*M_PI*log2 + 21.0*zeta3)/24.0;
    } else if (x > 0.5 && x < 1.0) {
        double log_x = log(x);
        double S12 = 0.0, lfactorial = 1.0;
        for (int l=0; l<19; l++) {
            if (l!=0) lfactorial *= (double)l;
            S12 += 0.5 * B0[l]/((double)l+2.0)/lfactorial 
                   * pow(-log_x, (double)l+2.0);
        }
        Li3 = - S12 - log_x*gsl_sf_dilog(1.0-x) - 0.5*log_x*log_x*log(1.0-x)
              + gsl_sf_zeta_int(2)*log_x + gsl_sf_zeta_int(3);
    } else if (x == 1.0) {
        Li3 = gsl_sf_zeta_int(3);
    } else {
        throw "x is out of range in Polylogarithms::Li3()";
    }
    return (Li3);
}







