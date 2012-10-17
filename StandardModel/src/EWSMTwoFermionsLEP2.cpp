/* 
 * File:   EWSMTwoFermionsLEP2.cpp
 * Author: mishima
 */

#include <stdexcept>
#include <cmath>
#include "EWSMTwoFermionsLEP2.h"
#include "EWSMThreeLoopEW.h"


EWSMTwoFermionsLEP2::EWSMTwoFermionsLEP2(const StandardModel& SM_i, 
                                         const bool bDebug_i) : bDebug(bDebug_i), 
        SM(SM_i), myOneLoopEW(SM_i) {
    bWeak = true;
    bWeakBox = true;
    bISR = true; 
    bQEDFSR = true; 
    bQCDFSR = true; 
}


//////////////////////////////////////////////////////////////////////// 

void EWSMTwoFermionsLEP2::setFlag(const std::string str, const bool flag) {
    if (str=="Weak")
        bWeak = flag;
    else if (str=="WeakBox")
        bWeakBox = flag;
    else if (str=="ISR")
        bISR = flag;
    else if (str=="QEDFSR")
        bQEDFSR = flag;
    else if (str=="QCDFSR")
        bQCDFSR = flag;
    else
        throw std::runtime_error("Error in EWSMTwoFermionsLEP2::setFlag()");             
}
    
//////////////////////////////////////////////////////////////////////// 

double EWSMTwoFermionsLEP2::sigma_l(const StandardModel::lepton l, 
                                    const double s, const double Mw, 
                                    const double GammaZ) const {
    double mf = ml(l);
    double I3f = SM.getLeptons(l).getIsospin();
    double Qf = SM.getLeptons(l).getCharge();

    return sigma(s, Mw, GammaZ, I3f, Qf, 0.0, mf, 1.0);
}


double EWSMTwoFermionsLEP2::sigma_q(const StandardModel::quark q, 
                                    const double s, const double Mw, 
                                    const double GammaZ) const {
    double mf = mq(q, sqrt(s));
    double I3f = SM.getQuarks(q).getIsospin();
    double Qf = SM.getQuarks(q).getCharge();
    double mfp;
    if (q==SM.TOP)
        throw std::runtime_error("Error in EWSMTwoFermionsLEP2::sigma_q()");
    else if (q==SM.BOTTOM) 
        mfp = mq(SM.TOP, sqrt(s));
    else 
        mfp = 0.0;

    return sigma(s, Mw, GammaZ, I3f, Qf, mfp, mf, 3.0);
}


double EWSMTwoFermionsLEP2::AFB_l(const StandardModel::lepton l, 
                                  const double s, const double Mw, 
                                  const double GammaZ) const {
    double mf = ml(l);
    double I3f = SM.getLeptons(l).getIsospin();
    double Qf = SM.getLeptons(l).getCharge();

    return AFB(s, Mw, GammaZ, I3f, Qf, 0.0, mf);
}


double EWSMTwoFermionsLEP2::AFB_q(const StandardModel::quark q, 
                                  const double s, const double Mw, 
                                  const double GammaZ) const {
    double mf = mq(q, sqrt(s));
    double I3f = SM.getQuarks(q).getIsospin();
    double Qf = SM.getQuarks(q).getCharge();
    double mfp;
    if (q==SM.TOP)
        throw std::runtime_error("Error in EWSMTwoFermionsLEP2::AFB_q()");
    else if (q==SM.BOTTOM) 
        mfp = mq(SM.TOP, sqrt(s));
    else 
        mfp = 0.0;

    return AFB(s, Mw, GammaZ, I3f, Qf, mfp, mf);
}


////////////////////////////////////////////////////////////////////////

double EWSMTwoFermionsLEP2::ml(const StandardModel::lepton l) const {
    return SM.getLeptons(l).getMass();
}


double EWSMTwoFermionsLEP2::mq(const StandardModel::quark q, const double mu, 
                               const orders order) const {
    switch(q) {
        case StandardModel::UP:
        case StandardModel::DOWN:
        case StandardModel::STRANGE:
            if (bDebug)
                return SM.getQuarks(q).getMass(); // for debug
            else
                return SM.Mrun(mu, SM.getQuarks(q).getMass_scale(), 
                        SM.getQuarks(q).getMass(), order);                    
        case StandardModel::CHARM:
        case StandardModel::BOTTOM:
            if (bDebug)
                return SM.getQuarks(q).getMass(); // for debug
            else
                return SM.Mrun(mu, SM.getQuarks(q).getMass(), order);
        case StandardModel::TOP:
            return SM.getMtpole(); // the pole mass
        default:
            throw std::runtime_error("Error in EWSMTwoFermionsLEP2::mq()"); 
    }
}


////////////////////////////////////////////////////////////////////////

complex EWSMTwoFermionsLEP2::Vpol_inv(const double s) const {
    complex V_inv;
    if (bDebug)
        V_inv = 1.0/complex(1.0715119759, -0.0186242179, false); // for debug
    else
        V_inv = SM.getAle()/SM.alphaMz(); //!!TEST

    return V_inv;    
}


complex EWSMTwoFermionsLEP2::chi_Z(const double s, const double Mw, 
                                   const double GammaZ) const {
    double Mz = SM.getMz();
    complex denom = complex(s - Mz*Mz, GammaZ/Mz*s, false);
    double prefactor = SM.getGF()*Mz*Mz/(sqrt(2.0)*8.0*M_PI*SM.getAle());
    
    return ( prefactor*s/denom );
}

    
complex EWSMTwoFermionsLEP2::rho_ef(const double s, const double Mw, const double I3f, 
                                    const double Qf, const double mfp) const {
    complex G = complex(1.0, 0.0, false);
    
    return G;       
}

complex EWSMTwoFermionsLEP2::kappa_e(const double s, const double Mw, const double I3f, 
                                     const double Qf) const {
    complex G = complex(1.0, 0.0, false);
    
    return G;       
}


complex EWSMTwoFermionsLEP2::kappa_f(const double s, const double Mw, const double I3f, 
                                     const double Qf, const double mfp) const {
    complex G = complex(1.0, 0.0, false);
    
    return G;       
}


complex EWSMTwoFermionsLEP2::kappa_ef(const double s, const double Mw, const double I3f, 
                                      const double Qf, const double mfp) const {
    complex G = complex(1.0, 0.0, false);
    
    return G;       
}


complex EWSMTwoFermionsLEP2::I2e(const double s, const double Mw, const double I3f, 
                                 const double Qf) const {
    double Mz = SM.getMz(), sW2 = 1.0 - Mw*Mw/(Mz*Mz);
    double alpha = SM.getAle()/Vpol_inv(s).real();
    double ReKappa_e = 1.0;
    return ( 35.0*alpha*alpha/18.0*( 1.0 - 8.0/3.0*ReKappa_e*sW2 ) );
}


complex EWSMTwoFermionsLEP2::I2f(const double s, const double Mw, const double I3f, 
                                 const double Qf, const double mfp) const {
    double Mz = SM.getMz(), sW2 = 1.0 - Mw*Mw/(Mz*Mz);
    double alpha = SM.getAle()/Vpol_inv(s).real();
    double ReKappa_f = 1.0;
    return ( 35.0*alpha*alpha/18.0*( 1.0 - 8.0/3.0*ReKappa_f*sW2 ) );
}


complex EWSMTwoFermionsLEP2::G_e(const double s, const double Mw, const double I3f, 
                                 const double Qf) const {
    double Mz = SM.getMz(), sW2 = 1.0 - Mw*Mw/(Mz*Mz);
    return ( 1.0 - 4.0*( kappa_e(s,Mw,I3f,Qf)*sW2 + I2e(s,Mw,I3f,Qf)) ); 
}


complex EWSMTwoFermionsLEP2::G_f(const double s, const double Mw, const double I3f, 
                                 const double Qf, const double mfp) const {
    double Mz = SM.getMz(), sW2 = 1.0 - Mw*Mw/(Mz*Mz);
    return ( 1.0 - 4.0*fabs(Qf)*( kappa_f(s,Mw,I3f,Qf,mfp)*sW2 + I2f(s,Mw,I3f,Qf,mfp)) );     
}


complex EWSMTwoFermionsLEP2::G_ef(const double s, const double Mw, const double I3f, 
                                  const double Qf, const double mfp) const {
    double Mz = SM.getMz(), sW2 = 1.0 - Mw*Mw/(Mz*Mz);
    return ( - 1.0 + G_e(s,Mw,I3f,Qf) + G_f(s,Mw,I3f,Qf,mfp) 
             + 16.0*fabs(Qf)
               *( kappa_ef(s,Mw,I3f,Qf,mfp)*sW2*sW2 
                 + (kappa_e(s,Mw,I3f,Qf)*I2e(s,Mw,I3f,Qf) 
                    + kappa_f(s,Mw,I3f,Qf,mfp)*I2f(s,Mw,I3f,Qf,mfp))*sW2 ) );
}

   
double EWSMTwoFermionsLEP2::G_1(const double s, const double Mw, 
                                const double GammaZ, const double I3f, 
                                const double Qf, const double mfp) const {
    complex Vpol = 1.0/Vpol_inv(s);
    complex rhoef, Ge, Gf, Gef;
    if (bWeak) {
        rhoef = rho_ef(s, Mw, I3f, Qf, mfp);
        Ge = G_e(s, Mw, I3f, Qf);
        Gf = G_f(s, Mw, I3f, Qf, mfp);
        Gef = G_ef(s, Mw, I3f, Qf, mfp);
    } else {
        double Mz = SM.getMz(), sW2 = 1.0 - Mw*Mw/(Mz*Mz);
        rhoef = 1.0;
        Ge = 1.0 - 4.0*sW2;
        Gf = 1.0 - 4.0*fabs(Qf)*sW2;
        Gef = - 1.0 + Ge + Gf + 16.0*fabs(Qf)*sW2*sW2;
    }
    
    return ( Qf*Qf*Vpol.abs2() 
             + 2.0*fabs(Qf)*(Vpol.conjugate()*rhoef*Gef*chi_Z(s,Mw,GammaZ)).real() 
             + rhoef.abs2()*(Gef.abs2() + Gf.abs2() + Ge.abs2() + 1.0) 
               *chi_Z(s,Mw,GammaZ).abs2() );
}


double EWSMTwoFermionsLEP2::G_2(const double s, const double Mw, 
                                const double GammaZ, const double I3f, 
                                const double Qf, const double mfp) const {
    complex Vpol = 1.0/Vpol_inv(s);
    complex rhoef, Ge, Gf, Gef;
    if (bWeak) {
        rhoef = rho_ef(s, Mw, I3f, Qf, mfp);
        Ge = G_e(s, Mw, I3f, Qf);
        Gf = G_f(s, Mw, I3f, Qf, mfp);
        Gef = G_ef(s, Mw, I3f, Qf, mfp);
    } else {
        double Mz = SM.getMz(), sW2 = 1.0 - Mw*Mw/(Mz*Mz);
        rhoef = 1.0;
        Ge = 1.0 - 4.0*sW2;
        Gf = 1.0 - 4.0*fabs(Qf)*sW2;
        Gef = - 1.0 + Ge + Gf + 16.0*fabs(Qf)*sW2*sW2;
    }   

    return ( Qf*Qf*Vpol.abs2() 
             + 2.0*fabs(Qf)*(Vpol.conjugate()*rhoef*Gef*chi_Z(s,Mw,GammaZ)).real() 
             + rhoef.abs2()*(Gef.abs2() + Gf.abs2())*chi_Z(s,Mw,GammaZ).abs2() ); 
}


double EWSMTwoFermionsLEP2::G_3(const double s, const double Mw, 
                                const double GammaZ, const double I3f, 
                                const double Qf, const double mfp) const {
    complex Vpol = 1.0/Vpol_inv(s);
    complex rhoef, Ge, Gf, Gef;
    if (bWeak) {
        rhoef = rho_ef(s, Mw, I3f, Qf, mfp);
        Ge = G_e(s, Mw, I3f, Qf);
        Gf = G_f(s, Mw, I3f, Qf, mfp);
        Gef = G_ef(s, Mw, I3f, Qf, mfp);
    } else {
        double Mz = SM.getMz(), sW2 = 1.0 - Mw*Mw/(Mz*Mz);
        rhoef = 1.0;
        Ge = 1.0 - 4.0*sW2;
        Gf = 1.0 - 4.0*fabs(Qf)*sW2;
        Gef = - 1.0 + Ge + Gf + 16.0*fabs(Qf)*sW2*sW2;
    }  

    return ( 2.0*fabs(Qf)*(Vpol.conjugate()*chi_Z(s,Mw,GammaZ)).real() 
             + 2.0*(Ge*Gf.conjugate() + Gef).real()*chi_Z(s,Mw,GammaZ).abs2() );
}


double EWSMTwoFermionsLEP2::sigma(const double s, const double Mw, 
                                  const double GammaZ, const double I3f, 
                                  const double Qf, const double mfp,
                                  const double mf, const double Ncf) const {
    double betaf = sqrt(1.0 - 4.0*mf*mf/s);
    double G1 = G_1(s, Mw, GammaZ, I3f, Qf, mfp);
    double G2 = G_2(s, Mw, GammaZ, I3f, Qf, mfp);

    return ( 4.0*M_PI*SM.getAle()*SM.getAle()/(3.0*s)*Ncf*betaf
             *( G1 + 2.0*mf*mf/s*G2) );
}


double EWSMTwoFermionsLEP2::AFB(const double s, const double Mw, 
                                const double GammaZ, const double I3f, 
                                const double Qf, const double mfp,
                                const double mf) const {
    double betaf = sqrt(1.0 - 4.0*mf*mf/s);
    double G1 = G_1(s, Mw, GammaZ, I3f, Qf, mfp);
    double G2 = G_2(s, Mw, GammaZ, I3f, Qf, mfp);
    double G3 = G_3(s, Mw, GammaZ, I3f, Qf, mfp);

    return ( 3.0/4.0*betaf*G3/( G1 + 2.0*mf*mf/s*G2) );
}


////////////////////////////////////////////////////////////////////////




