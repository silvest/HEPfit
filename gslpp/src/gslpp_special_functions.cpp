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

#include <gsl/gsl_sf_zeta.h>
#include <gsl/gsl_sf_dilog.h>
#include <gsl/gsl_sf_clausen.h>
#include "gslpp_special_functions.h"

double gslpp_special_functions::zeta(int i) {
    return gsl_sf_zeta_int(i);
}

double gslpp_special_functions::clausen(double x) {
    return gsl_sf_clausen(x);
}

double gslpp_special_functions::dilog(double x) {
    return gsl_sf_dilog(x);
}
