/*
 * Copyright (C) 2012 HEPfit Collaboration
 *
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
    return BR_mu_e_gamma;
}

log_meg::log_meg(const StandardModel& SM_i)
: ThObservable(SM_i), mySM(SM_i)
{}

double log_meg::computeThValue()
{
    double alph = mySM.getAle();
    double mMU = mySM.getLeptons(StandardModel::MU).getMass();
    gslpp::vector<gslpp::complex> ** allcoeff_me = mySM.getMyLeptonFlavour()->ComputeCoeffli_lj_gamma(1);
    double BR_mu_e_gamma = (alph*pow(mMU,5) * ((*(allcoeff_me[LO])) * (*(allcoeff_me[LO])).conjugate()).abs() / (2.99598e-19) );
    return log10(BR_mu_e_gamma);
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
    double BR_mu_e_gamma_obs=4.2e-13;  //from 1605.05081v1 (TheMEG:2016wtm)
    double BR_mu_e_gamma_exp=5.3e-13;  //from 1605.05081v1 (TheMEG:2016wtm)
    double ntos=1.64485362695147; //90% CL
    
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

log_tmg::log_tmg(const StandardModel& SM_i)
: ThObservable(SM_i), mySM(SM_i)
{}

double log_tmg::computeThValue()
{
    double alph = mySM.getAle();
    double mTAU = mySM.getLeptons(StandardModel::TAU).getMass();
    gslpp::vector<gslpp::complex> ** allcoeff_tm = mySM.getMyLeptonFlavour()->ComputeCoeffli_lj_gamma(2);
    double BR_tau_mu_gamma = (alph*pow(mTAU,5) * ((*(allcoeff_tm[LO])) * (*(allcoeff_tm[LO])).conjugate()).abs() / (2.26735e-12) );
    return log10(BR_tau_mu_gamma);
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
    double ntos=1.64485362695147; //90% CL
    
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

log_teg::log_teg(const StandardModel& SM_i)
: ThObservable(SM_i), mySM(SM_i)
{}

double log_teg::computeThValue()
{
    double alph = mySM.getAle();
    double mTAU = mySM.getLeptons(StandardModel::TAU).getMass();
    gslpp::vector<gslpp::complex> ** allcoeff_te = mySM.getMyLeptonFlavour()->ComputeCoeffli_lj_gamma(3);
    double BR_tau_e_gamma = (alph*pow(mTAU,5) * ((*(allcoeff_te[LO])) * (*(allcoeff_te[LO])).conjugate()).abs() / (2.26735e-12) );
    return log10(BR_tau_e_gamma);
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
    double ntos=1.64485362695147; //90% CL
    
    return (1+(BR_tau_e_gamma-BR_tau_e_gamma_obs)/BR_tau_e_gamma_exp)*ntos;
}
