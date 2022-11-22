/*
 * Copyright (C) 2012 HEPfit Collaboration
 *
 *
 * For the licensing terms see doc/COPYING.
 */

#include "EpsilonP_O_Epsilon.h"
#include "StandardModel.h"
#include "std_make_vector.h"

EpsilonP_O_Epsilon::EpsilonP_O_Epsilon(const StandardModel& SM_i, unsigned int part_i)
: ThObservable(SM_i), AmpDS1(SM_i)
{
    setParametersForObservable(make_vector<std::string>() << "Br_Ks_P0P0" << "Br_Ks_PpPm" << "Br_Kp_P0Pp" << "Omega_eta_etap" << "Delta_0" << "Delta_2" << "EpsK" << "phiEpsK"
                                                          << "Zqq00" << "Zqq11" << "Zqq12" << "Zqq13" << "Zqq14" << "Zqq21" << "Zqq22" << "Zqq23" << "Zqq24"
                                                          << "Zqq31" << "Zqq32" << "Zqq33" << "Zqq34" << "Zqq41" << "Zqq42" << "Zqq43" << "Zqq44" << "Zqq55"
                                                          << "Zqq56" << "Zqq65" << "Zqq66");
    part = part_i;
}

double EpsilonP_O_Epsilon::computeThValue()
{
  double phase = -sin(((SM.getOptionalParameter("Delta_2")-SM.getOptionalParameter("Delta_0"))-SM.getOptionalParameter("phiEpsK"))*M_PI/180.);
  switch (part) {
      case 0:
        // return ReA0 using exp info
        return getReA0();
      case 1:
        // return ReA2 using exp info
        return getReA2();
      case 2:
        // return ReA0 using lattice info
        return AmpDS1pp0pureLAT(NLO).real();
      case 3:
        // return ReA2 using lattice info
        return AmpDS1pp2(NLO).real();
      case 4:
        // eps'/eps: state-of-the-art from lattice + exp measurement of ReA0,2
        ReA0 = getReA0();
        ImA0 = AmpDS1pp0(NLO).imag();
        ReA2 = getReA2();
        ImA2 = AmpDS1pp2(NLO).imag();
        break;
      case 5:
        // eps'/eps: prediction solely based on lattice results
        ReA0 = AmpDS1pp0pureLAT(NLO).real();
        ImA0 = AmpDS1pp0pureLAT(NLO).imag();
        ReA2 = AmpDS1pp2(NLO).real();
        ImA2 = AmpDS1pp2(NLO).imag();
        break;
      case 6:
          // return Im A0 using exp info
          return AmpDS1pp0(NLO).imag();
      case 7:
          // return Im A0 using pure lattice
          return AmpDS1pp0pureLAT(NLO).imag();
      case 8:
          // return Im A2
          return AmpDS1pp2(NLO).imag();
      default:
          throw ("incorrect value of part in EpsilonP_O_Epsilon");
  }
  //Evaluate Re(eps'/eps) as defined in ArXiv:2004.09440
  return M_SQRT1_2 * phase * (ReA2/ReA0) * ( (ImA2/ReA2) - ((1.-SM.getOptionalParameter("Omega_eta_etap")) * (ImA0/ReA0)))/SM.getOptionalParameter("EpsK");
}
