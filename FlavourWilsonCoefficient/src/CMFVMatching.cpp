/* 
 * Copyright (C) 2016 HEPfit Collaboration
 *
 * For the licensing terms see doc/COPYING.
 */

#include "CMFVMatching.h"
#include "CMFV.h"
#include <stdexcept>

CMFVMatching::CMFVMatching(const CMFV & CMFV_i) :
    StandardModelMatching(CMFV_i),
    myCMFV(CMFV_i) {};

double CMFVMatching::S0(double x1, double x2) const {
    if (x1 > .02 && x2 > .02)
        return (myCMFV.getFtt() + StandardModelMatching::S0(x1, x2));
    return StandardModelMatching::S0(x1,x2);
}