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
    auto myphi = remainder(phi, 2.0 * M_PI);   // reduce phi to the range [-pi, pi]
    if (myphi < 0.0) myphi = -myphi; // Clausen function of third order is even, so we can use the positive value

    if (myphi==0.0) return ( gsl_sf_zeta_int(3) );

    double TMP = 0.0, l_double = 0.0, lfactorial = 1.0, sign = 1.0;
    for (int l=2; l<19; l++) {
        l_double = (double)l;
        lfactorial *= l_double;
        if (l%4) { sign = - 1.0; } else { sign = 1.0; }
        TMP += B[l]*pow(myphi,l_double)/l_double/(l_double + 1.0)/(l_double + 2.0)
               /lfactorial * sign;
    }
    return ( gsl_sf_zeta_int(3) - myphi*myphi*(3.0/4.0 - log(myphi)/2.0 - TMP) );
}


