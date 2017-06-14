/* 
 * Copyright (C) 2017 HEPfit Collaboration
 *
 *
 * For the licensing terms see doc/COPYING.
 */

/* 
 * File:   gslpp_special_functions.h
 * Author: marco
 *
 * Created on 14 marzo 2017, 17.56
 */

#ifndef GSLPP_SPECIAL_FUNCTIONS_H
#define GSLPP_SPECIAL_FUNCTIONS_H

#include <gsl/gsl_sf_zeta.h>
#include <gsl/gsl_sf_dilog.h>
#include <gsl/gsl_sf_clausen.h>

namespace gslpp_special_functions { 
    double zeta(int i) { return gsl_sf_zeta_int(i); }
    double clausen(double x) { return gsl_sf_clausen(x); }
    double dilog(double x) { return gsl_sf_dilog(x); }
}

#endif /* GSLPP_SPECIAL_FUNCTIONS_H */
