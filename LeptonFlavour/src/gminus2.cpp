/*
 * Copyright (C) 2015 HEPfit Collaboration
 * All rights reserved.
 *
 * For the licensing terms see doc/COPYING.
 */

#include "gminus2.h"
#include "StandardModel.h"

//gminus2::gminus2(const StandardModel& SM_i): ThObservable(SM_i)
//{
//};
//
//double gminus2::computeThValue()
//{
//    return 0.0;
//}

//gminus2_mu::gminus2_mu(const StandardModel& SM_i)
//: gminus2(SM_i), mySM(SM_i)
//{}

gminus2_mu::gminus2_mu(const StandardModel& SM_i)
: ThObservable(SM_i), mySM(SM_i)
{}

double gminus2_mu::computeThValue()
{
    gslpp::vector<gslpp::complex> ** allcoeff_gminus2mu = mySM.getMyLeptonFlavour()->ComputeCoeffgminus2mu();

    return ((*(allcoeff_gminus2mu[LO]))(0)+(*(allcoeff_gminus2mu[LO]))(1)).abs();
}
