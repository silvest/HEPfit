/* 
 * Copyright (C) 2012-2013 SusyFit Collaboration
 * All rights reserved.
 *
 * For the licensing terms see doc/COPYING.
 */

#include <NPEpsilons.h>
#include <NPEpsilons_pureNP.h>
#include "EW_ABC.h"


EW_ABC::EW_ABC(const StandardModel& SM_i)
: SM(SM_i), myEWNPEpsilons(SM_i)
{
    if (SM.ModelName().compare("NPEpsilons") != 0)
        throw std::runtime_error("EW_ABC::EW_ABC() is not applicable for the model " + SM.ModelName());
}


double EW_ABC::Mw(const bool bAlternative) const
{
    if (!bAlternative)
        return myEWNPEpsilons.Mw(eps1(), eps2(), eps3());
    else {
        double delta_alpha = (SM.alphaMz() - 1.0/128.90)/SM.getAle();
        double cW2_Born = 0.768905*(1.0 - 0.40*delta_alpha);
        double cW2 = cW2_Born*(1.0 + 1.43*eps1() - 1.00*eps2() - 0.86*eps3());
        return ( sqrt(cW2)*SM.getMz() );
    }
}


double EW_ABC::Gamma_l(StandardModel::lepton l) const 
{
    double Qf = SM.getLeptons(l).getCharge();
    double RQED = 1.0 + 3.0*SM.alphaMz()/4.0/M_PI*Qf*Qf;
    double mf = SM.getLeptons(l).getMass();
    double beta = sqrt(1.0 - 4.0*mf*mf/SM.getMz()/SM.getMz());
    double factor = SM.getGF()*SM.getMz()*SM.getMz()*SM.getMz()/6.0/M_PI/sqrt(2.0);
    return ( factor*beta*( (3.0 - beta*beta)/2.0*gVl(l).abs2()
                           + beta*beta*gAl(l).abs2() )*RQED );
}


double EW_ABC::Gamma_q(StandardModel::quark q) const 
{
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
             *((3.0 - beta*beta)/2.0*gVq(q).abs2()
               + beta*beta*gAq(q).abs2())*Nc*RQCD*RQED );
}


double EW_ABC::Gamma_b() const
{
    double RQED = 1.0 + SM.alphaMz()/12.0/M_PI;
    double a = SM.Als(SM.getMz(), FULLNNLO)/M_PI;    
    double RQCD = 1.0 + 1.2*a - 1.1*a*a - 13.0*a*a*a;
    //double mb = 4.78; // pole mass (PDG2012)
    double mb = 4.7; // used by Altarelli et al (1998)
    double beta = sqrt(1.0 - 4.0*mb*mb/SM.getMz()/SM.getMz());
    double Nc = 3.0; 
    double factor = SM.getGF()*SM.getMz()*SM.getMz()*SM.getMz()/6.0/M_PI/sqrt(2.0);
    
    return ( factor*beta
             *((3.0 - beta*beta)/2.0*gVb().abs2()
               + beta*beta*gAb().abs2())*Nc*RQCD*RQED );
}

    
double EW_ABC::Gamma_had() const
{
    return ( Gamma_q(SM.UP) + Gamma_q(SM.DOWN) 
             + Gamma_q(SM.CHARM) + Gamma_q(SM.STRANGE)  
             + Gamma_b() );
}


double EW_ABC::GammaZ(const bool bAlternative) const
{
    if (!bAlternative)
        return ( Gamma_l(SM.NEUTRINO_1) + Gamma_l(SM.NEUTRINO_2)
                 + Gamma_l(SM.NEUTRINO_3) + Gamma_l(SM.ELECTRON)
                 + Gamma_l(SM.MU) + Gamma_l(SM.TAU)
                 + Gamma_had() );
    else {
        double delta_als = (SM.Als(SM.getMz(),FULLNNLO) - 0.119)/M_PI;
        double delta_alpha = (SM.alphaMz() - 1.0/128.90)/SM.getAle();
        double Gamma_T0 = 2.48946*(1.0 + 0.73*delta_als - 0.35*delta_alpha);
        return ( Gamma_T0*(1.0 + 1.35*eps1() - 0.46*eps3() + 0.35*epsb()) );
    }
}

   
double EW_ABC::R_l(const bool bAlternative) const
{
    if (!bAlternative)
        return ( Gamma_had()/Gamma_l(SM.ELECTRON) );
    else {
        double delta_als = (SM.Als(SM.getMz(),FULLNNLO) - 0.119)/M_PI;
        double delta_alpha = (SM.alphaMz() - 1.0/128.90)/SM.getAle();
        double R_0 = 20.8228*(1.0 + 1.05*delta_als - 0.28*delta_alpha);
        return ( R_0*(1.0 + 0.28*eps1() - 0.36*eps3() + 0.50*epsb()) );
    }
}


double EW_ABC::R_c() const 
{
    return ( Gamma_q(SM.CHARM)/Gamma_had() );    
}


double EW_ABC::R_b(const bool bAlternative) const
{
    if (!bAlternative)
        return ( Gamma_b()/Gamma_had() );
    else {
        double R_b0 = 0.2182355;
        return ( R_b0*(1.0 - 0.06*eps1() + 0.07*eps3() + 1.79*epsb()) );
    }
}


double EW_ABC::sigma0_had(const bool bAlternative) const
{
    if (!bAlternative)
        return ( 12.0*M_PI/SM.getMz()/SM.getMz()
                 *Gamma_l(SM.ELECTRON)*Gamma_had()
                 /GammaZ(false)/GammaZ(false) );
    else {
        double delta_als = (SM.Als(SM.getMz(),FULLNNLO) - 0.119)/M_PI;
        double delta_alpha = (SM.alphaMz() - 1.0/128.90)/SM.getAle();
        double sigma_h0 = 41.420*(1.0 - 0.41*delta_als + 0.03*delta_alpha);
        return ( sigma_h0*(1.0 - 0.03*eps1() + 0.04*eps3() - 0.20*epsb()) );
    }
}


double EW_ABC::A_l(StandardModel::lepton l, const bool bAlternative) const
{
    if (!bAlternative) {
        double x = gVl_over_gAl(l).real();
        return ( 2.0*x/(1.0 + x*x) );
    } else {
        double delta_alpha = (SM.alphaMz() - 1.0/128.90)/SM.getAle();
        double x0 = 0.075619 - 1.32*delta_alpha;
        double x = x0*(1.0 + 17.6*eps1() - 22.9*eps3());
        return ( 2.0*x/(1.0 + x*x) );
    }
}


double EW_ABC::A_q(StandardModel::quark q) const 
{
    if (q==StandardModel::BOTTOM || q==StandardModel::TOP)
        throw std::runtime_error("Error in EW_ABC::A_q()");  
    double x = gVq_over_gAq(q).real();
    return ( 2.0*x/(1.0 + x*x) );    
}


double EW_ABC::A_b() const 
{
    double x = gVb_over_gAb().real();
    return ( 2.0*x/(1.0 + x*x) );
}


double EW_ABC::AFB_l(StandardModel::lepton l, const bool bAlternative) const
{
    if (!bAlternative) {
        double x = gVl_over_gAl(l).real();
        return ( 3.0*x*x/(1.0 + x*x)/(1.0 + x*x) );
    } else {
        double delta_als = (SM.Als(SM.getMz(),FULLNNLO) - 0.119)/M_PI;
        double AFB_l_Born = 0.01696*(1.0 - 34.0*delta_als);
        return ( AFB_l_Born*(1.0 + 34.72*eps1() - 45.15*eps3()) );
    }
}

double EW_ABC::AFB_c() const 
{
    double x = gVl_over_gAl(SM.ELECTRON).real();
    double xc = gVq_over_gAq(SM.CHARM).real();
    return ( 3.0*x*xc/(1.0 + x*x)/(1.0 + xc*xc) );    
}


double EW_ABC::AFB_b() const
{
    double x = gVl_over_gAl(SM.ELECTRON).real();
    double xb = gVb_over_gAb().real();
    return ( 3.0*x*xb/(1.0 + x*x)/(1.0 + xb*xb) );    
}

    
double EW_ABC::sin2thetaEff(const bool bAlternative) const
{
    if (!bAlternative) {
        double x = gVl_over_gAl(SM.ELECTRON).real();
        return ( (1.0 - x)/4.0 );
   } else {
        double delta_als = (SM.Als(SM.getMz(),FULLNNLO) - 0.119)/M_PI;
        double x0 = 0.075619 - 1.32*delta_als;
        double x = x0*(1.0 + 17.6*eps1() - 22.9*eps3());
        return ( (1.0 - x)/4.0 );
   }
}


////////////////////////////////////////////////////////////////////////

double EW_ABC::eps1() const
{
    if (SM.ModelName().compare("NPEpsilons") == 0)
        return (static_cast<const NPEpsilons*> (&SM))->epsilon1();
    else if (SM.ModelName().compare("NPEpsilons_pureNP") == 0)
        return (static_cast<const NPEpsilons_pureNP*> (&SM))->epsilon1();
    else
        throw std::runtime_error("EW_ABC::eps1() is not applicable to " + SM.ModelName());
}


double EW_ABC::eps2() const
{
    if (SM.ModelName().compare("NPEpsilons") == 0)
        return (static_cast<const NPEpsilons*> (&SM))->epsilon2();
    else if (SM.ModelName().compare("NPEpsilons_pureNP") == 0)
        return (static_cast<const NPEpsilons_pureNP*> (&SM))->epsilon2();
    else
        throw std::runtime_error("EW_ABC::eps2() is not applicable to " + SM.ModelName());
}


double EW_ABC::eps3() const
{
    if (SM.ModelName().compare("NPEpsilons") == 0)
        return (static_cast<const NPEpsilons*> (&SM))->epsilon3();
    else if (SM.ModelName().compare("NPEpsilons_pureNP") == 0)
        return (static_cast<const NPEpsilons_pureNP*> (&SM))->epsilon3();
    else
        throw std::runtime_error("EW_ABC::eps3() is not applicable to " + SM.ModelName());
}


double EW_ABC::epsb() const
{
    if (SM.ModelName().compare("NPEpsilons") == 0)
        return (static_cast<const NPEpsilons*> (&SM))->epsilonb();
    else if (SM.ModelName().compare("NPEpsilons_pureNP") == 0)
        return (static_cast<const NPEpsilons_pureNP*> (&SM))->epsilonb();
    else
        throw std::runtime_error("EW_ABC::epsb() is not applicable to " + SM.ModelName());
}


////////////////////////////////////////////////////////////////////////
    
complex EW_ABC::gVl(StandardModel::lepton l) const 
{
    return myEWNPEpsilons.gVl(l, eps1(), eps3());
}


complex EW_ABC::gAl(StandardModel::lepton l) const 
{
    return myEWNPEpsilons.gAl(l, eps1());
}


complex EW_ABC::gVl_over_gAl(StandardModel::lepton l) const
{
    return ( gVl(l)/gAl(l) );
}


complex EW_ABC::gVq(StandardModel::quark q) const 
{
    if (q==StandardModel::BOTTOM || q==StandardModel::TOP)
        throw std::runtime_error("Error in EW_ABC::gVq()");    
    return myEWNPEpsilons.gVq(q, eps1(), eps3());
}


complex EW_ABC::gAq(StandardModel::quark q) const 
{
    if (q==StandardModel::BOTTOM || q==StandardModel::TOP)
        throw std::runtime_error("Error in EW_ABC::gAq()");  
    return myEWNPEpsilons.gAq(q, eps1());
}


complex EW_ABC::gVq_over_gAq(StandardModel::quark q) const
{
    if (q==StandardModel::BOTTOM || q==StandardModel::TOP)
        throw std::runtime_error("Error in EW_ABC::gVq_over_gAq()");  
    return ( gVq(q)/gAq(q) );
}


complex EW_ABC::gVb() const 
{
    return myEWNPEpsilons.gVb(eps1(), eps3(), epsb());
}


complex EW_ABC::gAb() const 
{
    return myEWNPEpsilons.gAb(eps1(), epsb());
}


complex EW_ABC::gVb_over_gAb() const
{
    return ( gVb()/gAb() );
}

