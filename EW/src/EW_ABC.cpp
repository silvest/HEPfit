/* 
 * Copyright (C) 2012 SusyFit Collaboration
 * All rights reserved.
 *
 * For the licensing terms see doc/COPYING.
 */

#include "EW_ABC.h"


double EW_ABC::Mw(const double eps1, const double eps2, const double eps3) const {
    return myEWepsilons.Mw(eps1, eps2, eps3);
}


double EW_ABC::Gamma_l(StandardModel::lepton l, const double eps1, const double eps3) const {
    double Qf = SM.getLeptons(l).getCharge();
    double RQED = 1.0 + 3.0*SM.alphaMz()/4.0/M_PI*Qf*Qf;
    double mf = SM.getLeptons(l).getMass();
    double beta = sqrt(1.0 - 4.0*mf*mf/SM.getMz()/SM.getMz());
    double factor = SM.getGF()*SM.getMz()*SM.getMz()*SM.getMz()/6.0/M_PI/sqrt(2.0);
    return ( factor*beta*( (3.0 - beta*beta)/2.0*gVl(l,eps1,eps3).abs2()
                           + beta*beta*gAl(l,eps1).abs2() )*RQED );
}


double EW_ABC::Gamma_q(StandardModel::quark q, const double eps1, const double eps3) const {
    if (q==StandardModel::BOTTOM || q==StandardModel::TOP)
        throw std::runtime_error("Error in EW_ABC::Gamma_q()");  
    double Qf = SM.getQuarks(q).getCharge();
    double RQED = 1.0 + 3.0*SM.alphaMz()/4.0/M_PI*Qf*Qf;
    double a = SM.Als(SM.getMz(), FULLNNLO)/M_PI;
    double RQCD = 1.0 + 1.2*a - 1.1*a*a - 13.0*a*a*a;
    double mf = 0.0;
    if (q==StandardModel::CHARM)
        mf = 1.67; // pole mass (PDG2012)
    else
        mf = 0.0;
        //mf = SM.Mrun(SM.getMz(),SM.getQuarks(q).getMass_scale(),SM.getQuarks(q).getMass(),FULLNNLO);
    double beta = sqrt(1.0 - 4.0*mf*mf/SM.getMz()/SM.getMz());
    double Nc = 3.0; 
    double factor = SM.getGF()*SM.getMz()*SM.getMz()*SM.getMz()/6.0/M_PI/sqrt(2.0);
    
    return ( factor*beta
             *((3.0 - beta*beta)/2.0*gVq(q,eps1,eps3).abs2()
               + beta*beta*gAq(q,eps1).abs2())*Nc*RQCD*RQED );
}


double EW_ABC::Gamma_b(const double eps1, const double eps3, const double epsb) const {
    double RQED = 1.0 + SM.alphaMz()/12.0/M_PI;
    double a = SM.Als(SM.getMz(), FULLNNLO)/M_PI;    
    double RQCD = 1.0 + 1.2*a - 1.1*a*a - 13.0*a*a*a;
    //double mb = 4.78; // pole mass (PDG2012)
    double mb = 4.7; // used by Altarelli et al (1998)
    double beta = sqrt(1.0 - 4.0*mb*mb/SM.getMz()/SM.getMz());
    double Nc = 3.0; 
    double factor = SM.getGF()*SM.getMz()*SM.getMz()*SM.getMz()/6.0/M_PI/sqrt(2.0);
    
    return ( factor*beta
             *((3.0 - beta*beta)/2.0*gVb(eps1,eps3,epsb).abs2()
               + beta*beta*gAb(eps1,epsb).abs2())*Nc*RQCD*RQED );
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
    double x = gVl_over_gAl(l,eps1,eps3).real();
    return ( 2.0*x/(1.0 + x*x) );    
}


double EW_ABC::A_q(StandardModel::quark q, const double eps1, const double eps3) const {
    if (q==StandardModel::BOTTOM || q==StandardModel::TOP)
        throw std::runtime_error("Error in EW_ABC::A_q()");  
    double x = gVq_over_gAq(q,eps1,eps3).real();
    return ( 2.0*x/(1.0 + x*x) );    
}


double EW_ABC::A_b(const double eps1, const double eps3, const double epsb) const {
    double x = gVb_over_gAb(eps1,eps3,epsb).real();
    return ( 2.0*x/(1.0 + x*x) );
}


double EW_ABC::AFB_l(StandardModel::lepton l, const double eps1, const double eps3) const {
    double x = gVl_over_gAl(l,eps1,eps3).real();
    return ( 3.0*x*x/(1.0 + x*x)/(1.0 + x*x) );
}

double EW_ABC::AFB_c(const double eps1, const double eps3) const {
    double x = gVl_over_gAl(SM.ELECTRON,eps1,eps3).real();
    double xc = gVq_over_gAq(SM.CHARM,eps1,eps3).real();
    return ( 3.0*x*xc/(1.0 + x*x)/(1.0 + xc*xc) );    
}


double EW_ABC::AFB_b(const double eps1, const double eps3, const double epsb) const {
    double x = gVl_over_gAl(SM.ELECTRON,eps1,eps3).real();
    double xb = gVb_over_gAb(eps1,eps3,epsb).real();
    return ( 3.0*x*xb/(1.0 + x*x)/(1.0 + xb*xb) );    
}

    
double EW_ABC::sin2thetaEff(const double eps1, const double eps3) const {
    double x = gVl_over_gAl(SM.ELECTRON,eps1,eps3).real();
    return ( (1.0 - x)/4.0 );
}

    
complex EW_ABC::gVl(StandardModel::lepton l, const double eps1, const double eps3) const {
    return myEWepsilons.gVl(l, eps1, eps3);
}


complex EW_ABC::gAl(StandardModel::lepton l, const double eps1) const {
    return myEWepsilons.gAl(l, eps1);
}


complex EW_ABC::gVl_over_gAl(StandardModel::lepton l, const double eps1, const double eps3) const {
    return ( gVl(l,eps1,eps3)/gAl(l,eps1) );
}


complex EW_ABC::gVq(StandardModel::quark q, const double eps1, const double eps3) const {
    if (q==StandardModel::BOTTOM || q==StandardModel::TOP)
        throw std::runtime_error("Error in EW_ABC::gVq()");    
    return myEWepsilons.gVq(q, eps1, eps3);
}


complex EW_ABC::gAq(StandardModel::quark q, const double eps1) const {
    if (q==StandardModel::BOTTOM || q==StandardModel::TOP)
        throw std::runtime_error("Error in EW_ABC::gAq()");  
    return myEWepsilons.gAq(q, eps1);
}


complex EW_ABC::gVq_over_gAq(StandardModel::quark q, const double eps1, const double eps3) const {
    if (q==StandardModel::BOTTOM || q==StandardModel::TOP)
        throw std::runtime_error("Error in EW_ABC::gVq_over_gAq()");  
    return ( gVq(q,eps1,eps3)/gAq(q,eps1) );
}


complex EW_ABC::gVb(const double eps1, const double eps3, const double epsb) const {
    return myEWepsilons.gVb(eps1, eps3, epsb);
}


complex EW_ABC::gAb(const double eps1, const double epsb) const {
    return myEWepsilons.gAb(eps1, epsb);
}


complex EW_ABC::gVb_over_gAb(const double eps1, const double eps3, const double epsb) const {
    return ( gVb(eps1,eps3,epsb)/gAb(eps1,epsb) );
}

