/*
 * Copyright (C) 2012 HEPfit Collaboration
 * All rights reserved.
 *
 * For the licensing terms see doc/COPYING.
 */

#include "li_lj_gamma.h"
#include "StandardModel.h"
// #include "SUSYMatching.h"

li_lj_gamma::li_lj_gamma(const StandardModel& SM_i): ThObservable(SM_i), mySM(SM_i){
//li_lj_gamma::li_lj_gamma(const StandardModel& SM_i, int obsFlag): ThObservable(SM_i), mySM(SM_i){
//    if (obsFlag > 0 and obsFlag < 4) obs = obsFlag;
//    else throw std::runtime_error("li_to_lj can only be"
//            "1 (mu to electron) or"
//            "2 (tau to mu) or"
//            "3 (tau to electron)");
};

void li_lj_gamma::computeParameters()
{

    double alph = mySM.getAle();
    double mMU = mySM.getLeptons(StandardModel::MU).getMass();
    double mTAU = mySM.getLeptons(StandardModel::TAU).getMass();

        gslpp::vector<complex> ** allcoeff_me = mySM.getMyLeptonFlavour()->ComputeCoeffli_lj_gamma();
        BR_mu_e_gamma = (alph*pow(mMU,5.0) * ((*(allcoeff_me[LO])) * (*(allcoeff_me[LO])).conjugate()).abs() / (3.0e-19) );
//        BR_tau_e_gamma = (alph*pow(mTAU,5.0) * ((*(allcoeff_te[LO])) * (*(allcoeff_te[LO])).conjugate()).abs() / (3.0e-19) );
//        BR_tau_mu_gamma = (alph*pow(mTAU,5.0) * ((*(allcoeff_tm[LO])) * (*(allcoeff_tm[LO])).conjugate()).abs() / (3.0e-19) );

//        return (alph*ml * ((*(allcoeff[LO])) * (*(allcoeff[LO])).conjugate()).abs());

//        gslpp::vector<complex> C7m2e(2, 0.);
//        C7m2e.assign(0, mySUSYMatching.C7_Lepton(1)(0));
//        C7m2e.assign(1, mySUSYMatching.C7_Lepton(1)(1));
////        C7(obs).assign( 0, mySUSYMatching.C7_Lepton(obs));
////        C7(obs).assign( 1, mySUSYMatching.C7_Lepton(obs));
//        return (alph*pow(mMU,5.0) * ( C7m2e(0).abs2() + C7m2e(1).abs2() ));

//    throw std::runtime_error("EWPO::computeThValue(): Observable type not defined. Can be only any of (1) or (2) or (3)");
//    return (EXIT_FAILURE);

}

double li_lj_gamma::computeThValue()
{
    return 0;
}

mu_e_gamma::mu_e_gamma(const StandardModel& SM_i)
: li_lj_gamma(SM_i)
{}

double mu_e_gamma::computeThValue()
{
    computeParameters();
    return BR_mu_e_gamma;
}

tau_e_gamma::tau_e_gamma(const StandardModel& SM_i)
: li_lj_gamma(SM_i)
{}

double tau_e_gamma::computeThValue()
{
//    computeParameters();
    return 0.0;
}

tau_mu_gamma::tau_mu_gamma(const StandardModel& SM_i)
: li_lj_gamma(SM_i)
{}

double tau_mu_gamma::computeThValue()
{
//    computeParameters();
    return 0.0;
}

