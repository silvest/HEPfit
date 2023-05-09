/*
 * Copyright (C) 2023 HEPfit Collaboration
 *
 *
 * For the licensing terms see doc/COPYING.
 */



#include "Kmunu_o_Pmunu.h"
#include "StandardModel.h"
#include "std_make_vector.h"



Kmunu_o_Pmunu::Kmunu_o_Pmunu(const StandardModel& SM_i)
: ThObservable(SM_i)
{
    

    setParametersForObservable(make_vector<std::string>() << "fKofP" << "delKl2Pl2em");
    SM.initializeMeson(QCD::P_P);
    SM.initializeMeson(QCD::K_P);
    
};


double Kmunu_o_Pmunu::computeThValue()
{
    //The WC are written in the LEFT basis of arxiv:1709.04486 the expressions can be found in arxiv:1706.00410 and arxiv:1605.07114 in a similar basis
    gslpp::vector<gslpp::complex> ** allcoeff_K_mu_nu = SM.getFlavour().ComputeCoeffsleptonnu(QCD::K_P, QCD::MU);
    gslpp::vector<gslpp::complex> ** allcoeff_P_mu_nu = SM.getFlavour().ComputeCoeffuleptonnu(QCD::P_P, QCD::MU);


    double mK = SM.getMesons(QCD::K_P).getMass();
    double mP = SM.getMesons(QCD::P_P).getMass();
    
    
    double mu = SM.getQuarks(QCD::UP).getMass();
    double md = SM.getQuarks(QCD::DOWN).getMass();
    double ms = SM.getQuarks(QCD::STRANGE).getMass();
    double mmu = SM.getLeptons(StandardModel::MU).getMass();

    

    
    //double delKl2Pl2em = -0.0070;
    //double fKofP = 1.1932;
    
    
    double fKofP = SM.getOptionalParameter("fKofP");
    double delKl2Pl2em = SM.getOptionalParameter("delKl2Pl2em");
    
    
    double WCnum = ((*(allcoeff_K_mu_nu[LO]))(0) -(*(allcoeff_K_mu_nu[LO]))(1) + mK * mK / (ms+mu) / mmu * ((*(allcoeff_K_mu_nu[LO]))(2) - (*(allcoeff_K_mu_nu[LO]))(3))).abs2();
    double WCden = ((*(allcoeff_P_mu_nu[LO]))(0) -(*(allcoeff_P_mu_nu[LO]))(1) + mP * mP / (md+mu) / mmu * ((*(allcoeff_P_mu_nu[LO]))(2) - (*(allcoeff_P_mu_nu[LO]))(3))).abs2();
    
    //double WCnum = 1.0;
    //double WCden = 1.0;
   
    //return 1.;
    return mK/mP * (1-mmu*mmu/mK/mK) * (1-mmu*mmu/mK/mK)/(1-mmu*mmu/mP/mP)/(1-mmu*mmu/mP/mP)*fKofP*fKofP*(1+delKl2Pl2em)*WCnum/WCden;

}
