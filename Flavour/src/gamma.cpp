/* 
 * Copyright (C) 2018 HEPfit Collaboration
 *
 *
 * For the licensing terms see doc/COPYING.
 */

#include "gamma.h"
#include "StandardModel.h"

CKMGamma::CKMGamma(const StandardModel& SM_i) : ThObservable(SM_i) 
{}

double CKMGamma::computeThValue() 
{
    return(SM.getCKM().computeGamma()/M_PI*180.);
}

CKMGamma_rad::CKMGamma_rad(const StandardModel& SM_i) : ThObservable(SM_i) 
{}

double CKMGamma_rad::computeThValue() 
{
    return(SM.getCKM().computeGamma());
}