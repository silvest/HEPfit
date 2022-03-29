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
    setParametersForObservable(make_vector<std::string>() << "Br_Ks_P0P0" << "Br_Ks_PpPm" << "Br_Kp_P0Pp" << "Omega_eta_etap" << "Delta_0" << "Delta_2" << "EpsK" << "phiEpsK"
                                                          << "Zqq00" << "Zqq11" << "Zqq12" << "Zqq13" << "Zqq14" << "Zqq21" << "Zqq22" << "Zqq23" << "Zqq24"
                                                          << "Zqq31" << "Zqq32" << "Zqq33" << "Zqq34" << "Zqq41" << "Zqq42" << "Zqq43" << "Zqq44" << "Zqq55"
                                                          << "Zqq56" << "Zqq65" << "Zqq66");
}

double EpsilonP_O_Epsilon::computeThValue()
{
  //Evaluate Re(eps'/eps) as defined in ArXiv:2004.09440
  double phase = -sin(((SM.getOptionalParameter("Delta_2")-SM.getOptionalParameter("Delta_0"))-SM.getOptionalParameter("phiEpsK"))*M_PI/180.);
  return M_SQRT1_2 * phase * (getReA2()/getReA0()) * ( (AmpDS1pp2(NLO).imag() / getReA2()) - ( (1.-SM.getOptionalParameter("Omega_eta_etap")) * (AmpDS1pp0(NLO).imag() / getReA0()) ) )/SM.getOptionalParameter("EpsK");
}
