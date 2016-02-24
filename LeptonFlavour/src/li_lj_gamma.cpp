/*
 * Copyright (C) 2012 HEPfit Collaboration
 * All rights reserved.
 *
 * For the licensing terms see doc/COPYING.
 */

#include "li_lj_gamma.h"
#include "StandardModel.h"

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

Robs_mu_e_gamma::Robs_mu_e_gamma(const StandardModel& SM_i)
: ThObservable(SM_i), mySM(SM_i)
{}

double Robs_mu_e_gamma::computeThValue()
{
    double alph = mySM.getAle();
    double mMU = mySM.getLeptons(StandardModel::MU).getMass();
    gslpp::vector<gslpp::complex> ** allcoeff_me = mySM.getMyLeptonFlavour()->ComputeCoeffli_lj_gamma(1);
    double BR_mu_e_gamma = (alph*pow(mMU,5) * ((*(allcoeff_me[LO])) * (*(allcoeff_me[LO])).conjugate()).abs() / (2.99598e-19) );
    double BR_mu_e_gamma_obs=0.57e-12;  //from Phys.Rev.Lett. 110 (2013) 201801 (Adam:2013mnn)
    double BR_mu_e_gamma_exp=0.77e-12;  //from Phys.Rev.Lett. 110 (2013) 201801 (Adam:2013mnn)
    double ntos=1.64485362695147;
    
    return (1+(BR_mu_e_gamma-BR_mu_e_gamma_obs)/BR_mu_e_gamma_exp)*ntos;
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

Robs_tau_mu_gamma::Robs_tau_mu_gamma(const StandardModel& SM_i)
: ThObservable(SM_i), mySM(SM_i)
{}

double Robs_tau_mu_gamma::computeThValue()
{
    double alph = mySM.getAle();
    double mTAU = mySM.getLeptons(StandardModel::TAU).getMass();
    gslpp::vector<gslpp::complex> ** allcoeff_tm = mySM.getMyLeptonFlavour()->ComputeCoeffli_lj_gamma(2);
    double BR_tau_mu_gamma = (alph*pow(mTAU,5) * ((*(allcoeff_tm[LO])) * (*(allcoeff_tm[LO])).conjugate()).abs() / (2.26735e-12) );
    double BR_tau_mu_gamma_obs=4.4e-8;    //from Phys.Rev.Lett. 104 (2010) 021802 (Aubert:2009ag)
    double BR_tau_mu_gamma_exp=8.2e-8;    //from Phys.Rev.Lett. 104 (2010) 021802 (Aubert:2009ag)
    double ntos=1.64485362695147;
    
    return (1+(BR_tau_mu_gamma-BR_tau_mu_gamma_obs)/BR_tau_mu_gamma_exp)*ntos;
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

Robs_tau_e_gamma::Robs_tau_e_gamma(const StandardModel& SM_i)
: ThObservable(SM_i), mySM(SM_i)
{}

double Robs_tau_e_gamma::computeThValue()
{
    double alph = mySM.getAle();
    double mTAU = mySM.getLeptons(StandardModel::TAU).getMass();
    gslpp::vector<gslpp::complex> ** allcoeff_te = mySM.getMyLeptonFlavour()->ComputeCoeffli_lj_gamma(3);
    double BR_tau_e_gamma = (alph*pow(mTAU,5) * ((*(allcoeff_te[LO])) * (*(allcoeff_te[LO])).conjugate()).abs() / (2.26735e-12) );
    double BR_tau_e_gamma_obs=3.3e-8;    //from Phys.Rev.Lett. 104 (2010) 021802 (Aubert:2009ag)
    double BR_tau_e_gamma_exp=9.8e-8;    //from Phys.Rev.Lett. 104 (2010) 021802 (Aubert:2009ag)
    double ntos=1.64485362695147;
    
    return (1+(BR_tau_e_gamma-BR_tau_e_gamma_obs)/BR_tau_e_gamma_exp)*ntos;
}
