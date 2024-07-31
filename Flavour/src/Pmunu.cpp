/*
 * Copyright (C) 2024 HEPfit Collaboration
 *
 *
 * For the licensing terms see doc/COPYING.
 */



#include "Pmunu.h"
#include "StandardModel.h"
#include "std_make_vector.h"



Pmunu::Pmunu(const StandardModel& SM_i)
: ThObservable(SM_i)
{
    

    setParametersForObservable(make_vector<std::string>() << "fK" << "delPl2em" << "fKofP");
    SM.initializeMeson(QCD::P_P);
    
};


double Pmunu::computeThValue()
{
    //The WC are written in the LEFT basis of arxiv:1709.04486 the expressions can be found in arxiv:1706.00410 and arxiv:1605.07114
    gslpp::vector<gslpp::complex> ** allcoeff_P_mu_nu = SM.getFlavour().ComputeCoeffdiujlknu(0,0,1,SM.getQuarks(QCD::UP).getMass_scale());

    double mP = SM.getMesons(QCD::P_P).getMass();
    
    double mu = SM.getQuarks(QCD::UP).getMass();
    double md = SM.getQuarks(QCD::DOWN).getMass();
    double mmu = SM.getLeptons(StandardModel::MU).getMass();
    
    double fP = SM.getOptionalParameter("fK") / SM.getOptionalParameter("fKofP"); // this is the pion decay constant in the SU(2)-symmetric limit 
    double delPl2em = SM.getOptionalParameter("delPl2em"); //this is the non-perturbative EM + strong isospin violation correction to the decay rate
    
    double amp2 = ((*(allcoeff_P_mu_nu[LO]))(0) -(*(allcoeff_P_mu_nu[LO]))(1) + mP * mP / (md+mu) / mmu * ((*(allcoeff_P_mu_nu[LO]))(2) - (*(allcoeff_P_mu_nu[LO]))(3))).abs2();
    
    //double WCnum = 1.0;
    //double WCden = 1.0;
   
    //return 1.;
    return 1. / (64. * M_PI) * mmu * mmu * mP * (1-mmu*mmu/mP/mP) * (1-mmu*mmu/mP/mP) * fP * fP * (1+delPl2em) * amp2 / HCUT;

}