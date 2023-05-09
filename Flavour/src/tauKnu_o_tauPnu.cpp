/* 
 * Copyright (C) 2023 HEPfit Collaboration
 *
 *
 * For the licensing terms see doc/COPYING.
 */


#include "tauKnu_o_tauPnu.h"
#include "StandardModel.h"
#include "std_make_vector.h"




tauKnu_o_tauPnu::tauKnu_o_tauPnu(const StandardModel& SM_i)
: ThObservable(SM_i)
{
    

    setParametersForObservable(make_vector<std::string>() << "fKofP" << "deltauK2tauP2em");
    SM.initializeMeson(QCD::P_P);
    SM.initializeMeson(QCD::K_P);
    
};


double tauKnu_o_tauPnu::computeThValue()
{
    //The WC are written in the LEFT basis of arxiv:1709.04486 the expressions can be found in arxiv:1706.00410 and arxiv:1605.07114 in a similar basis
    gslpp::vector<gslpp::complex> ** allcoeff_tau_K_nu = SM.getFlavour().ComputeCoeffsleptonnu(QCD::K_P, QCD::TAU);
    gslpp::vector<gslpp::complex> ** allcoeff_tau_P_nu = SM.getFlavour().ComputeCoeffuleptonnu(QCD::P_P, QCD::TAU);


    double mK = SM.getMesons(QCD::K_P).getMass();
    double mP = SM.getMesons(QCD::P_P).getMass();
    
    
    double mu = SM.getQuarks(QCD::UP).getMass();
    double md = SM.getQuarks(QCD::DOWN).getMass();
    double ms = SM.getQuarks(QCD::STRANGE).getMass();
    double mtau = SM.getLeptons(StandardModel::TAU).getMass();

    

    
    //double  deltauK2tauP2em= 0.0004;
    //double fKofP = 1.1932;
    
    
    double fKofP = SM.getOptionalParameter("fKofP");
    double deltauK2tauP2em = SM.getOptionalParameter("deltauK2tauP2em");
    
    
    double WCnum = ((*(allcoeff_tau_K_nu[LO]))(0) -(*(allcoeff_tau_K_nu[LO]))(1) + mK * mK / (ms+mu) / mtau * ((*(allcoeff_tau_K_nu[LO]))(2) - (*(allcoeff_tau_K_nu[LO]))(3))).abs2();
    double WCden = ((*(allcoeff_tau_P_nu[LO]))(0) -(*(allcoeff_tau_P_nu[LO]))(1) + mP * mP / (md+mu) / mtau * ((*(allcoeff_tau_P_nu[LO]))(2) - (*(allcoeff_tau_P_nu[LO]))(3))).abs2();
    
    //double WCnum = 1.0;
    //double WCden = 1.0;
   
    //return 1.;
    return (1-mK*mK/mtau/mtau) * (1-mK*mK/mtau/mtau)/(1-mP*mP/mtau/mtau)/(1-mP*mP/mtau/mtau)*fKofP*fKofP*(1+deltauK2tauP2em)*WCnum/WCden;

}