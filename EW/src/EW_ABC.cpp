/* 
 * Copyright (C) 2012 SUSYfit Collaboration
 * All rights reserved.
 *
 * For the licensing terms see doc/COPYING.
 */

#include "EW_ABC.h"


double EW_ABC::Mw(const double eps1, const double eps2, const double eps3) const {
    return ( SM.getMz()/sqrt(2.0)
             *sqrt(1.0 + sqrt(1.0 - 4.0*M_PI*SM.alphaMz()
                                    /(sqrt(2.0)*SM.getGF()*SM.getMz()*SM.getMz()
                                      *(1.0 - Delta_rW(eps1,eps2,eps3))))) );
}


double EW_ABC::Gamma_l(StandardModel::lepton l, const double eps1, const double eps3) const {
    double factor = SM.getGF()*SM.getMz()*SM.getMz()*SM.getMz()/6.0/M_PI/sqrt(2.0);
    return ( factor*gAl(l,eps1)*gAl(l,eps1)
             *(1.0 + gVl_over_gAl(l,eps1,eps3)*gVl_over_gAl(l,eps1,eps3))
             *(1.0 + 3.0/4.0*SM.alphaMz()/M_PI) );
}


double EW_ABC::Gamma_q(StandardModel::quark q, const double eps1, const double eps3) const {
    if (q==StandardModel::BOTTOM || q==StandardModel::TOP)
        throw std::runtime_error("Error in EW_ABC::Gamma_q()");  
    double a = SM.Als(SM.getMz(), FULLNNLO)/M_PI;
    double RQCD = 1.0 + 1.2*a - 1.1*a*a - 13.0*a*a*a;
    double beta = 1.0;
    double Nc = 3.0; 
    double factor = SM.getGF()*SM.getMz()*SM.getMz()*SM.getMz()/6.0/M_PI/sqrt(2.0);
    
    return ( factor*beta*((3.0 - beta*beta)/2.0*gVq(q,eps1,eps3)*gVq(q,eps1,eps3)
                          + beta*beta*gAq(q,eps1)*gAq(q,eps1))
             *Nc*RQCD*(1.0 + SM.alphaMz()/12.0/M_PI) );
}


double EW_ABC::Gamma_b(const double eps1, const double eps3, const double epsb) const {
    double a = SM.Als(SM.getMz(), FULLNNLO)/M_PI;
    double RQCD = 1.0 + 1.2*a - 1.1*a*a - 13.0*a*a*a;
    double mb = 4.7;
    double beta = sqrt(1.0 - 4.0*mb*mb/SM.getMz()/SM.getMz());
    double Nc = 3.0; 
    double factor = SM.getGF()*SM.getMz()*SM.getMz()*SM.getMz()/6.0/M_PI/sqrt(2.0);
    
    return ( factor*beta*((3.0 - beta*beta)/2.0*gVb(eps1,eps3,epsb)*gVb(eps1,eps3,epsb)
                          + beta*beta*gAb(eps1,epsb)*gAb(eps1,epsb))
             *Nc*RQCD*(1.0 + SM.alphaMz()/12.0/M_PI) );
}

    
double EW_ABC::Gamma_had(const double eps1, const double eps3, const double epsb) const {
    return ( Gamma_q(SM.UP,eps1,eps3) + Gamma_q(SM.DOWN,eps1,eps3) 
             + Gamma_q(SM.CHARM,eps1,eps3) + Gamma_q(SM.STRANGE,eps1,eps3)  
             + Gamma_b(eps1,eps3,epsb) );
}


double EW_ABC::GammaZ(const double eps1, const double eps3, const double epsb) const {
    return ( Gamma_l(SM.NEUTRINO_1,eps1,eps3) + Gamma_l(SM.NEUTRINO_2,eps1,eps3) 
             + Gamma_l(SM.NEUTRINO_3,eps1,eps3) + Gamma_l(SM.ELECTRON,eps1,eps3)
             + Gamma_l(SM.MU,eps1,eps3) + Gamma_l(SM.TAU,eps1,eps3) 
             + Gamma_had(eps1,eps3,epsb) );
}

   
double EW_ABC::R_l(const double eps1, const double eps3, const double epsb) const {
    return ( Gamma_had(eps1,eps3,epsb)/Gamma_l(SM.ELECTRON,eps1,eps3) );
}


double EW_ABC::R_c(const double eps1, const double eps3, const double epsb) const {
    return ( Gamma_q(SM.CHARM,eps1,eps3)/Gamma_had(eps1,eps3,epsb) );    
}


double EW_ABC::R_b(const double eps1, const double eps3, const double epsb) const {
    return ( Gamma_b(eps1,eps3,epsb)/Gamma_had(eps1,eps3,epsb) );        
}


double EW_ABC::sigma0_had(const double eps1, const double eps3, const double epsb) const {
    return ( 12.0*M_PI/SM.getMz()/SM.getMz()
             *Gamma_l(SM.ELECTRON,eps1,eps3)*Gamma_had(eps1,eps3,epsb)
             /GammaZ(eps1,eps3,epsb)/GammaZ(eps1,eps3,epsb) );
}


double EW_ABC::A_l(StandardModel::lepton l, const double eps1, const double eps3) const {
    double x = gVl_over_gAl(l,eps1,eps3);
    return ( 2.0*x/(1.0 + x*x) );    
}


double EW_ABC::A_q(StandardModel::quark q, const double eps1, const double eps3) const {
    if (q==StandardModel::BOTTOM || q==StandardModel::TOP)
        throw std::runtime_error("Error in EW_ABC::A_q()");  
    double x = gVq_over_gAq(q,eps1,eps3);
    return ( 2.0*x/(1.0 + x*x) );    
}


double EW_ABC::A_b(const double eps1, const double eps3, const double epsb) const {
    double x = gVb_over_gAb(eps1,eps3,epsb);
    return ( 2.0*x/(1.0 + x*x) );
}


double EW_ABC::AFB_l(StandardModel::lepton l, const double eps1, const double eps3) const {
    double x = gVl_over_gAl(l,eps1,eps3);
    return ( 3.0*x*x/(1.0 + x*x)/(1.0 + x*x) );
}

double EW_ABC::AFB_c(const double eps1, const double eps3) const {
    double x = gVq_over_gAq(SM.CHARM,eps1,eps3);
    return ( 3.0*x*x/(1.0 + x*x)/(1.0 + x*x) );    
}


double EW_ABC::AFB_b(const double eps1, const double eps3, const double epsb) const {
    double x = gVb_over_gAb(eps1,eps3,epsb);
    return ( 3.0*x*x/(1.0 + x*x)/(1.0 + x*x) );    
}

    
double EW_ABC::sin2thetaEff(const double eps1, const double eps3) const {
    double x = gVl_over_gAl(SM.ELECTRON,eps1,eps3);
    return ( (1.0 - x)/4.0 );
}

    
double EW_ABC::gVl(StandardModel::lepton l, const double eps1, const double eps3) const {
    return ( gVl_over_gAl(l,eps1,eps3)*gAl(l,eps1) );
}


double EW_ABC::gAl(StandardModel::lepton l, const double eps1) const {
    return ( SM.getLeptons(l).getIsospin()*(1.0 + Delta_rho(eps1)/2.0) );
}


double EW_ABC::gVl_over_gAl(StandardModel::lepton l, const double eps1, const double eps3) const {
    double Ql = SM.getLeptons(l).getCharge();
    return ( 1.0 - 4.0*fabs(Ql)*(1.0 + Delta_kappa(eps1,eps3))*SM.s02() );
}


double EW_ABC::gVq(StandardModel::quark q, const double eps1, const double eps3) const {
    if (q==StandardModel::BOTTOM || q==StandardModel::TOP)
        throw std::runtime_error("Error in EW_ABC::gVq()");    
    return ( gVq_over_gAq(q,eps1,eps3)*gAq(q,eps1) );
}


double EW_ABC::gAq(StandardModel::quark q, const double eps1) const {
    if (q==StandardModel::BOTTOM || q==StandardModel::TOP)
        throw std::runtime_error("Error in EW_ABC::gAq()");  
    return ( SM.getQuarks(q).getIsospin()*(1.0 + Delta_rho(eps1)/2.0) );
}


double EW_ABC::gVq_over_gAq(StandardModel::quark q, const double eps1, const double eps3) const {
    if (q==StandardModel::BOTTOM || q==StandardModel::TOP)
        throw std::runtime_error("Error in EW_ABC::gVq_over_gAq()");  
    double Qq = SM.getQuarks(q).getCharge();
    return ( 1.0 - 4.0*fabs(Qq)*(1.0 + Delta_kappa(eps1,eps3))*SM.s02() );
}


double EW_ABC::gVb(const double eps1, const double eps3, const double epsb) const {
    return ( gVb_over_gAb(eps1,eps3,epsb)*gAb(eps1,epsb) );
}


double EW_ABC::gAb(const double eps1, const double epsb) const {
    return ( - (1.0 + Delta_rho(eps1)/2.0)/2.0*(1.0 + epsb) );
}


double EW_ABC::gVb_over_gAb(const double eps1, const double eps3, const double epsb) const {
    double sin2theff = (1.0 + Delta_kappa(eps1,eps3))*SM.s02();
    return ( (1.0 - 4.0/3.0*sin2theff + epsb)/(1.0 + epsb) );
}


double EW_ABC::Delta_rW(const double eps1, const double eps2, const double eps3) const {
    return ( (eps2 - SM.c02()*Delta_rho(eps1) + 2.0*SM.s02()*Delta_kappa(eps1,eps3))
              *(SM.c02() - SM.s02())/SM.s02() );
}


double EW_ABC::Delta_rho(const double eps1) const {
    return eps1;
}


double EW_ABC::Delta_kappa(const double eps1, const double eps3) const {
    return ( (eps3 - SM.c02()*Delta_rho(eps1))/(SM.c02() - SM.s02()) );
}




