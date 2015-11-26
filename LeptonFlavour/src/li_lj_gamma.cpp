/*
 * Copyright (C) 2012 HEPfit Collaboration
 * All rights reserved.
 *
 * For the licensing terms see doc/COPYING.
 */

#include "li_lj_gamma.h"
#include "StandardModel.h"

//li_lj_gamma::li_lj_gamma(const StandardModel& SM_i): ThObservable(SM_i)
//{
//};
//
//double li_lj_gamma::computeThValue()
//{
//    return 0.0;
//}

mu_e_gamma::mu_e_gamma(const StandardModel& SM_i)
: ThObservable(SM_i), mySM(SM_i)
{}

double mu_e_gamma::computeThValue()
{
    double alph = mySM.getAle();
    double mMU = mySM.getLeptons(StandardModel::MU).getMass();
    gslpp::vector<gslpp::complex> ** allcoeff_me = mySM.getMyLeptonFlavour()->ComputeCoeffli_lj_gamma(1);
    double BR_mu_e_gamma = (alph*pow(mMU,5) * ((*(allcoeff_me[LO])) * (*(allcoeff_me[LO])).conjugate()).abs() / (2.99598e-19) );
//    double BR_mu_e_gamma = log10((alph*pow(mMU,5.0) * ((*(allcoeff_me[LO])) * (*(allcoeff_me[LO])).conjugate()).abs() / (3.0e-19) )); /*This is the log!*/
    return BR_mu_e_gamma;
}

tau_mu_gamma::tau_mu_gamma(const StandardModel& SM_i)
: ThObservable(SM_i), mySM(SM_i)
{}

double tau_mu_gamma::computeThValue()
{
    double alph = mySM.getAle();
    double mTAU = mySM.getLeptons(StandardModel::TAU).getMass();
    gslpp::vector<gslpp::complex> ** allcoeff_tm = mySM.getMyLeptonFlavour()->ComputeCoeffli_lj_gamma(2);
    double BR_tau_mu_gamma = (alph*pow(mTAU,5) * ((*(allcoeff_tm[LO])) * (*(allcoeff_tm[LO])).conjugate()).abs() / (2.26735e-12) );
    return BR_tau_mu_gamma;
}

tau_e_gamma::tau_e_gamma(const StandardModel& SM_i)
: ThObservable(SM_i), mySM(SM_i)
{}

double tau_e_gamma::computeThValue()
{
    double alph = mySM.getAle();
    double mTAU = mySM.getLeptons(StandardModel::TAU).getMass();
    gslpp::vector<gslpp::complex> ** allcoeff_te = mySM.getMyLeptonFlavour()->ComputeCoeffli_lj_gamma(3);
    double BR_tau_e_gamma = (alph*pow(mTAU,5) * ((*(allcoeff_te[LO])) * (*(allcoeff_te[LO])).conjugate()).abs() / (2.26735e-12) );
    return BR_tau_e_gamma;
}
