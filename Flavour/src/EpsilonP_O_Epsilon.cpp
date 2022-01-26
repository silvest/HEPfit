/* 
 * Copyright (C) 2012 HEPfit Collaboration
 *
 *
 * For the licensing terms see doc/COPYING.
 */

#include "EpsilonP_O_Epsilon.h"
#include "StandardModel.h"
#include "std_make_vector.h"

EpsilonP_O_Epsilon::EpsilonP_O_Epsilon(const StandardModel& SM_i)
: ThObservable(SM_i), AmpDS1(SM_i)
{
    setParametersForObservable(make_vector<std::string>() << "ReA0_Kd" << "ReA2_Kd" << "Omega_eta_etap" << "MP0" << "FP0");
};

double EpsilonP_O_Epsilon::computeThValue()
{
    return (M_SQRT1_2 * (SM.getOptionalParameter("ReA2_Kd") / SM.getOptionalParameter("ReA0_Kd") / SM.getOptionalParameter("ReA0_Kd")) * ((1 / SM.getOptionalParameter("ReA2_Kd") / SM.getOptionalParameter("ReA0_Kd")) * AmpDS1pp2(NLO).imag() - (1 - SM.getOptionalParameter("Omega_eta_etap")) * AmpDS1pp0(NLO).imag()));
}
