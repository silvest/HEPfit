/* 
 * Copyright (C) 2012 HEPfit Collaboration
 *
 *
 * For the licensing terms see doc/COPYING.
 */

#include <cstdlib>
#include <stdexcept>
#include <cmath>
#include <gsl/gsl_sf.h>
#include "ClausenFunctions.h"


ClausenFunctions::ClausenFunctions() 
{
}

////////////////////////////////////////////////////////////////////////

double ClausenFunctions::Cl2(const double phi) const 
{
    return ( gsl_sf_clausen(phi) );    
}

double ClausenFunctions::Cl3(const double phi) const 
{
    phi = remainder(phi, 2.0 * M_PI);   // reduce phi to the range [-pi, pi]
    if (phi < 0.0) phi = -phi; // Clausen function of third order is even, so we can use the positive value
    
    if (phi==0.0) return ( gsl_sf_zeta_int(3) );
    
    double TMP = 0.0, l_double = 0.0, lfactorial = 1.0, sign = 1.0;
    for (int l=2; l<19; l++) {
        l_double = (double)l;
        lfactorial *= l_double;
        if (l%4) { sign = - 1.0; } else { sign = 1.0; }
        TMP += B[l]*pow(phi,l_double)/l_double/(l_double + 1.0)/(l_double + 2.0)
               /lfactorial * sign;
    }
    return ( gsl_sf_zeta_int(3) - phi*phi*(3.0/4.0 - log(phi)/2.0 - TMP) );
}


