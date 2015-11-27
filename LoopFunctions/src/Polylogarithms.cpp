/* 
 * Copyright (C) 2012 HEPfit Collaboration
 * All rights reserved.
 *
 * For the licensing terms see doc/COPYING.
 */

#include <cstdlib>
#include <stdexcept>
#include <cmath>
#include <gsl/gsl_complex.h>
#include <gsl/gsl_sf.h>
#include "Polylogarithms.h"


Polylogarithms::Polylogarithms() 
{
}

////////////////////////////////////////////////////////////////////////

gslpp::complex Polylogarithms::Li2(const double x) const 
{
    gsl_sf_result re, im;
    gsl_sf_complex_dilog_xy_e(x, 0.0, &re, &im);
    return gslpp::complex(re.val, im.val, false);
}

gslpp::complex Polylogarithms::Li2(const gslpp::complex z) const 
{
    gsl_sf_result re, im;
    gsl_sf_complex_dilog_xy_e(z.real(), z.imag(), &re, &im);
    return gslpp::complex(re.val, im.val, false);
}

double Polylogarithms::Li3(const double x) const 
{
    double Li3 = 0.0;
    if (x < 0.0)
        Li3 = -gsl_sf_fermi_dirac_2(log(-x));
    else if (x == 0.0)
        Li3 = 0.0;    
    else if (x > 0.0 && x < 0.5) {
        double log_1mx = log(1.0 - x);
        double lfactorial = 1.0, kfactorial = 1.0;
        for (int l=0; l<19; l++) {
            if (l!=0) lfactorial *= (double)l;
            kfactorial = 1.0;
            for (int k=0; k<19; k++) {
                if (k!=0) kfactorial *= (double)k;            
                Li3 += B[l]*B[k]/((double)l+1.0)/((double)l+(double)k+1.0)
                       /lfactorial/kfactorial 
                       * pow(-log_1mx, (double)l+(double)k+1.0);
            }
        }
    } else if (x == 0.5) {
        double log2 = log(2.0);
        double zeta3 = gsl_sf_zeta_int(3);
        Li3 = (4.0*pow(log2, 3.0) - 2.0*M_PI*M_PI*log2 + 21.0*zeta3)/24.0;
    } else if (x > 0.5 && x < 1.0) {
        double log_x = log(x);
        double S12 = 0.0, lfactorial = 1.0;
        for (int l=0; l<19; l++) {
            if (l!=0) lfactorial *= (double)l;
            S12 += 0.5 * B[l]/((double)l+2.0)/lfactorial 
                   * pow(-log_x, (double)l+2.0);
        }
        Li3 = - S12 - log_x*gsl_sf_dilog(1.0-x) - 0.5*log_x*log_x*log(1.0-x)
              + gsl_sf_zeta_int(2)*log_x + gsl_sf_zeta_int(3);
    } else if (x == 1.0)
        Li3 = gsl_sf_zeta_int(3);
    else
        throw std::runtime_error("Polylogarithms::Li3(): x is out of range!");
    return (Li3);
}



