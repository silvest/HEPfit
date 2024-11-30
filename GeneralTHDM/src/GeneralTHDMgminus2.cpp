/*
 * Copyright (C) 2017 HEPfit Collaboration
 *
 *
 * For the licensing terms see doc/COPYING.
 */

#include "GeneralTHDMgminus2.h"
#include "LeptonFlavour.h"

GeneralTHDMgminus2_mu::GeneralTHDMgminus2_mu(const StandardModel& SM_i)
: ThObservable(SM_i)
{}

double GeneralTHDMgminus2_mu::computeThValue()
{

    gslpp::vector<gslpp::complex> ** allcoeff_gminus2mu = SM.getMyLeptonFlavour()->ComputeCoeffgminus2mu();
    
    //Temporary implementation. This returns directly the coefficients that are 
    //evaluated to the \Delta a_\mu contributions from 1502.04199
    return ((*(allcoeff_gminus2mu[NLO]))(0)+(*(allcoeff_gminus2mu[NLO]))(1)).real();
}
