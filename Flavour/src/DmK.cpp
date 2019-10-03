/* 
 * Copyright (C) 2012 HEPfit Collaboration
 *
 *
 * For the licensing terms see doc/COPYING.
 */

#include "DmK.h"
#include "StandardModel.h"
#include "std_make_vector.h"

DmK::DmK(const StandardModel& SM_i) : ThObservable(SM_i), AmpDK2(SM_i) 
{
    setParametersForObservable(make_vector<std::string>() << "DmkSM");
};

double DmK::computeThValue() 
{
    return(SM.getCDMK()* (2.*AmpDMKNP(FULLNLO).real() + SM.getOptionalParameter("DmkSM")));
}
