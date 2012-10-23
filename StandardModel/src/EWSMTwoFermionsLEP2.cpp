/* 
 * File:   EWSMTwoFermionsLEP2.cpp
 * Author: mishima
 */

#include <stdexcept>
#include <cmath>
#include "EWSMTwoFermionsLEP2.h"
#include "EWSMThreeLoopEW.h"


EWSMTwoFermionsLEP2::EWSMTwoFermionsLEP2(const StandardModel& SM_i, 
                                         const bool bKeepNonUnitary_i) : SM(SM_i), 
        myCache(SM_i), myOneLoopEW(myCache) {
    bDebug = SM_i.isBDebug();
    bKeepNonUnitary = bKeepNonUnitary_i;
}


//////////////////////////////////////////////////////////////////////// 

double EWSMTwoFermionsLEP2::G_1(const double s, const double t, 
                                const double Mw, const double GammaZ, 
                                const double I3f, const double Qf, 
                                const double mf, const double mfp, 
                                const bool bWeak, const bool bWWbox, 
                                const bool bZZbox) const {
    complex Vpol = V_pol(s);
    complex rhoef, Ge, Gf, Gef;
    if (bWeak) {
        rhoef = rho_ef(s, t, Mw, I3f, Qf, mf, mfp, bWWbox, bZZbox);
        Ge = G_e(s, t, Mw, I3f, Qf, mf, mfp, bWWbox, bZZbox);
        Gf = G_f(s, t, Mw, I3f, Qf, mf, mfp, bWWbox, bZZbox);
        Gef = G_ef(s, t, Mw, I3f, Qf, mf, mfp, bWWbox, bZZbox);
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


double EWSMTwoFermionsLEP2::G_2(const double s, const double t, 
                                const double Mw, const double GammaZ, 
                                const double I3f, const double Qf, 
                                const double mf, const double mfp, 
                                const bool bWeak, const bool bWWbox, 
                                const bool bZZbox) const {
    complex Vpol = V_pol(s);
    complex rhoef, Ge, Gf, Gef;
    if (bWeak) {
        rhoef = rho_ef(s, t, Mw, I3f, Qf, mf, mfp, bWWbox, bZZbox);
        Ge = G_e(s, t, Mw, I3f, Qf, mf, mfp, bWWbox, bZZbox);
        Gf = G_f(s, t, Mw, I3f, Qf, mf, mfp, bWWbox, bZZbox);
        Gef = G_ef(s, t, Mw, I3f, Qf, mf, mfp, bWWbox, bZZbox);
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


double EWSMTwoFermionsLEP2::G_3(const double s, const double t, 
                                const double Mw, const double GammaZ, 
                                const double I3f, const double Qf, 
                                const double mf, const double mfp, 
                                const bool bWeak, const bool bWWbox, 
                                const bool bZZbox) const {
    complex Vpol = V_pol(s);
    complex rhoef, Ge, Gf, Gef;
    if (bWeak) {
        rhoef = rho_ef(s, t, Mw, I3f, Qf, mf, mfp, bWWbox, bZZbox);
        Ge = G_e(s, t, Mw, I3f, Qf, mf, mfp, bWWbox, bZZbox);
        Gf = G_f(s, t, Mw, I3f, Qf, mf, mfp, bWWbox, bZZbox);
        Gef = G_ef(s, t, Mw, I3f, Qf, mf, mfp, bWWbox, bZZbox);
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


double EWSMTwoFermionsLEP2::G_1_noBox(const double s, const double Mw, 
                                      const double GammaZ, const double I3f, 
                                      const double Qf, const double mf, 
                                      const double mfp, const bool bWeak) const {
    bool bWWbox = false, bZZbox = false;
    double t = 0.0; //dummy
    return ( G_1(s, t, Mw, GammaZ, I3f, Qf, mf, mfp, bWeak, bWWbox, bZZbox) );
}


double EWSMTwoFermionsLEP2::G_2_noBox(const double s, const double Mw, 
                                      const double GammaZ, const double I3f, 
                                      const double Qf, const double mf, 
                                      const double mfp, const bool bWeak) const {
    bool bWWbox = false, bZZbox = false;
    double t = 0.0; //dummy
    return ( G_2(s, t, Mw, GammaZ, I3f, Qf, mf, mfp, bWeak, bWWbox, bZZbox) );
}


double EWSMTwoFermionsLEP2::G_3_noBox(const double s, const double Mw, 
                                      const double GammaZ, const double I3f, 
                                      const double Qf, const double mf, 
                                      const double mfp, const bool bWeak) const {
    bool bWWbox = false, bZZbox = false;
    double t = 0.0; //dummy
    return ( G_3(s, t, Mw, GammaZ, I3f, Qf, mf, mfp, bWeak, bWWbox, bZZbox) );
}


double EWSMTwoFermionsLEP2::G_1_box(const double s, const double t, 
                                    const double Mw, const double GammaZ, 
                                    const double I3f, const double Qf, 
                                    const double mf, const double mfp) const {
    //double u = 2.0*mf*mf - s - t;
    double u = - s - t;
        
    double Mz = SM.getMz(), sW2 = 1.0 - Mw*Mw/(Mz*Mz); 
    double mu = Mw; // renormalization scale
    complex Vpol = V_pol(s);
    complex D_rho_ef = Delta_rho_ef_WW_hat(s, t, u, Mw, I3f)
                       + Delta_rho_ef_ZZ(mu, s, t, u, Mw, I3f, Qf);
    complex D_kappa_e = Delta_kappa_ef_WW_hat(s, t, u, Mw, I3f)
                        + Delta_kappa_e_ZZ(mu, s, t, u, Mw, I3f, Qf);
    complex D_kappa_f = Delta_kappa_e_WW_hat(s, t, u, Mw, I3f)
                        + Delta_kappa_f_ZZ(mu, s, t, u, Mw, I3f, Qf);
    complex D_kappa_ef = Delta_kappa_f_WW_hat(s, t, u, Mw, I3f)
                         + Delta_kappa_ef_ZZ(mu, s, t, u, Mw, I3f, Qf);
    if (I3f==SM.getQuarks(SM.BOTTOM).getIsospin() 
            && Qf==SM.getQuarks(SM.BOTTOM).getCharge() 
            && mfp!=0.0) {
        D_rho_ef += Delta_rho_ef_WW_TOP_hat(s, t, u, Mw);
        D_kappa_e += Delta_kappa_e_WW_TOP_hat(s, t, u, Mw);
        D_kappa_f += Delta_kappa_f_WW_TOP_hat(s, t, u, Mw);
        D_kappa_ef += Delta_kappa_ef_WW_TOP_hat(s, t, u, Mw);
    }
    
    return ( 2.0*fabs(Qf)
             *( Vpol.conjugate()*(D_rho_ef - 4.0*sW2*D_kappa_e 
                                  - 4.0*fabs(Qf)*sW2*D_kappa_f 
                                  + 16.0*fabs(Qf)*sW2*sW2*D_kappa_ef)
                *chi_Z(s,Mw,GammaZ) ).real()
             + 8.0*( D_rho_ef - 2.0*sW2*D_kappa_e - 2.0*fabs(Qf)*sW2*D_kappa_f 
                     + 4.0*fabs(Qf)*sW2*sW2*D_kappa_ef  ).real()
               *chi_Z(s,Mw,GammaZ).abs2() );
}


double EWSMTwoFermionsLEP2::G_2_box(const double s, const double t, 
                                    const double Mw, const double GammaZ, 
                                    const double I3f, const double Qf, 
                                    const double mf, const double mfp) const {
    //double u = 2.0*mf*mf - s - t;
    double u = - s - t;
        
    double Mz = SM.getMz(), sW2 = 1.0 - Mw*Mw/(Mz*Mz); 
    double mu = Mw; // renormalization scale
    complex Vpol = V_pol(s);
    complex D_rho_ef = Delta_rho_ef_WW_hat(s, t, u, Mw, I3f)
                       + Delta_rho_ef_ZZ(mu, s, t, u, Mw, I3f, Qf);
    complex D_kappa_e = Delta_kappa_ef_WW_hat(s, t, u, Mw, I3f)
                        + Delta_kappa_e_ZZ(mu, s, t, u, Mw, I3f, Qf);
    complex D_kappa_f = Delta_kappa_e_WW_hat(s, t, u, Mw, I3f)
                        + Delta_kappa_f_ZZ(mu, s, t, u, Mw, I3f, Qf);
    complex D_kappa_ef = Delta_kappa_f_WW_hat(s, t, u, Mw, I3f)
                         + Delta_kappa_ef_ZZ(mu, s, t, u, Mw, I3f, Qf);
    if (I3f==SM.getQuarks(SM.BOTTOM).getIsospin() 
            && Qf==SM.getQuarks(SM.BOTTOM).getCharge() 
            && mfp!=0.0) {
        D_rho_ef += Delta_rho_ef_WW_TOP_hat(s, t, u, Mw);
        D_kappa_e += Delta_kappa_e_WW_TOP_hat(s, t, u, Mw);
        D_kappa_f += Delta_kappa_f_WW_TOP_hat(s, t, u, Mw);
        D_kappa_ef += Delta_kappa_ef_WW_TOP_hat(s, t, u, Mw);
    }
        
    return ( 2.0*fabs(Qf)
             *( Vpol.conjugate()*(D_rho_ef - 4.0*sW2*D_kappa_e 
                                  - 4.0*fabs(Qf)*sW2*D_kappa_f 
                                  + 16.0*fabs(Qf)*sW2*sW2*D_kappa_ef)
                *chi_Z(s,Mw,GammaZ) ).real()
             + 4.0*( D_rho_ef - 2.0*sW2*D_kappa_e - 4.0*fabs(Qf)*sW2*D_kappa_f 
                     + 8.0*fabs(Qf)*sW2*sW2*D_kappa_ef  ).real()
               *chi_Z(s,Mw,GammaZ).abs2() );
}


double EWSMTwoFermionsLEP2::G_3_box(const double s, const double t, 
                                    const double Mw, const double GammaZ, 
                                    const double I3f, const double Qf, 
                                    const double mf, const double mfp) const {
    //double u = 2.0*mf*mf - s - t;
    double u = - s - t;

    double Mz = SM.getMz(), sW2 = 1.0 - Mw*Mw/(Mz*Mz); 
    double mu = Mw; // renormalization scale
    complex Vpol = V_pol(s);
    complex D_rho_ef = Delta_rho_ef_WW_hat(s, t, u, Mw, I3f)
                       + Delta_rho_ef_ZZ(mu, s, t, u, Mw, I3f, Qf);
    complex D_kappa_e = Delta_kappa_ef_WW_hat(s, t, u, Mw, I3f)
                        + Delta_kappa_e_ZZ(mu, s, t, u, Mw, I3f, Qf);
    complex D_kappa_f = Delta_kappa_e_WW_hat(s, t, u, Mw, I3f)
                        + Delta_kappa_f_ZZ(mu, s, t, u, Mw, I3f, Qf);
    complex D_kappa_ef = Delta_kappa_f_WW_hat(s, t, u, Mw, I3f)
                         + Delta_kappa_ef_ZZ(mu, s, t, u, Mw, I3f, Qf);
    if (I3f==SM.getQuarks(SM.BOTTOM).getIsospin() 
            && Qf==SM.getQuarks(SM.BOTTOM).getCharge() 
            && mfp!=0.0) {
        D_rho_ef += Delta_rho_ef_WW_TOP_hat(s, t, u, Mw);
        D_kappa_e += Delta_kappa_e_WW_TOP_hat(s, t, u, Mw);
        D_kappa_f += Delta_kappa_f_WW_TOP_hat(s, t, u, Mw);
        D_kappa_ef += Delta_kappa_ef_WW_TOP_hat(s, t, u, Mw);
    }
        
    return ( 2.0*fabs(Qf)
             *( Vpol.conjugate()*D_rho_ef*chi_Z(s,Mw,GammaZ) ).real()
             + 8.0*( D_rho_ef - 3.0*sW2*D_kappa_e - 3.0*fabs(Qf)*sW2*D_kappa_f 
                     + 4.0*fabs(Qf)*sW2*sW2*D_kappa_ef  ).real()
               *chi_Z(s,Mw,GammaZ).abs2() );
}


//////////////////////////////////////////////////////////////////////// 

complex EWSMTwoFermionsLEP2::V_pol(const double s) const {
    complex V;
    if (bDebug)
        V = complex(1.0715119759, -0.0186242179, false); // for debug
    else {
        V = SM.ale_OS(sqrt(s), FULLNLO)/SM.getAle() + myOneLoopEW.DeltaAlpha_t(s);
        //V = SM.ale_OS(sqrt(s), FULLNLO)/SM.getAle();
        //V = complex(1.0715119759, -0.0186242179, false); //!!TEST
    }
    return V;    
}


complex EWSMTwoFermionsLEP2::chi_Z(const double s, const double Mw, 
                                   const double GammaZ) const {
    double Mz = SM.getMz();
    complex denom = complex(s - Mz*Mz, GammaZ/Mz*s, false);
    double prefactor = SM.getGF()*Mz*Mz/(sqrt(2.0)*8.0*M_PI*SM.getAle());
    
    return ( prefactor*s/denom );
}


complex EWSMTwoFermionsLEP2::G_e(const double s, const double t, 
                                 const double Mw, const double I3f, 
                                 const double Qf, const double mf, 
                                 const double mfp, const bool bWWbox, 
                                 const bool bZZbox) const {
    double Mz = SM.getMz(), sW2 = 1.0 - Mw*Mw/(Mz*Mz);
    return ( 1.0 - 4.0*( kappa_e(s,t,Mw,I3f,Qf,mf,mfp,bWWbox,bZZbox)*sW2 
             + I2e(s,Mw)) ); 
}


complex EWSMTwoFermionsLEP2::G_f(const double s, const double t, 
                                 const double Mw, const double I3f, 
                                 const double Qf, const double mf, 
                                 const double mfp, const bool bWWbox, 
                                 const bool bZZbox) const {
    double Mz = SM.getMz(), sW2 = 1.0 - Mw*Mw/(Mz*Mz);
    return ( 1.0 - 4.0*fabs(Qf)*( kappa_f(s,t,Mw,I3f,Qf,mf,mfp,bWWbox,bZZbox)*sW2 
             + I2f(s,Mw)) );     
}


complex EWSMTwoFermionsLEP2::G_ef(const double s, const double t, 
                                  const double Mw, const double I3f, 
                                  const double Qf, const double mf, 
                                  const double mfp, const bool bWWbox, 
                                  const bool bZZbox) const {
    double Mz = SM.getMz(), sW2 = 1.0 - Mw*Mw/(Mz*Mz);
    return ( - 1.0 + G_e(s,t,Mw,I3f,Qf,mf,mfp,bWWbox,bZZbox) 
             + G_f(s,t,Mw,I3f,Qf,mf,mfp,bWWbox,bZZbox) 
             + 16.0*fabs(Qf)
               *( kappa_ef(s,t,Mw,I3f,Qf,mf,mfp,bWWbox,bZZbox)*sW2*sW2 
                 + (kappa_e(s,t,Mw,I3f,Qf,mf,mfp,bWWbox,bZZbox)*I2e(s,Mw) 
                    + kappa_f(s,t,Mw,I3f,Qf,mf,mfp,bWWbox,bZZbox)*I2f(s,Mw))*sW2 ) );
}


complex EWSMTwoFermionsLEP2::rho_ef(const double s, const double t, 
                                    const double Mw, const double I3f, 
                                    const double Qf, const double mf, 
                                    const double mfp, const bool bWWbox, 
                                    const bool bZZbox) const {
    double u = 2.0*mf*mf - s - t;
    double Mz = SM.getMz(), cW2 = Mw*Mw/(Mz*Mz), sW2 = 1.0 - cW2;
    double ve = - 0.5 + 2.0*sW2, ae = -0.5;
    double vf = I3f - 2.0*Qf*sW2, af = I3f;
    double mu = Mw; // renormalization scale

    complex rhoef = 1.0;
    rhoef += SM.getAle()/4.0/M_PI/sW2
             *( - DeltaRhobarZ(mu, Mw) + D_Z_hat(s, Mw) 
                + 5.0/3.0*PV.B0(mu, s, Mw, Mw) - 9.0*cW2/4.0/sW2*log(cW2) 
                - 6.0 + 5.0*cW2/8.0*(1.0 + cW2) 
                + (3.0*ve*ve + ae*ae + 3.0*vf*vf + af*af)/4.0/cW2*F_za_0(s, Mw)
                + 2.0*F_W_0_hat(s, Mw) );
    
    // WW box    
    if (bWWbox)
        rhoef += Delta_rho_ef_WW_hat(s, t, u, Mw, I3f);
    
    // ZZ box
    if (bZZbox)
        rhoef += Delta_rho_ef_ZZ(mu, s, t, u, Mw, I3f, Qf);
    
    // Top quark contribution for e^+ e^- -> b bbar
    if (I3f==SM.getQuarks(SM.BOTTOM).getIsospin() 
            && Qf==SM.getQuarks(SM.BOTTOM).getCharge() 
            && mfp!=0.0)
        rhoef += Delta_rho_ef_TOP(s, t, u, Mw, bWWbox);
    
    return rhoef;
}

complex EWSMTwoFermionsLEP2::kappa_e(const double s, const double t, 
                                     const double Mw, const double I3f, 
                                     const double Qf, const double mf, 
                                     const double mfp, const bool bWWbox, 
                                     const bool bZZbox) const {
    double u = 2.0*mf*mf - s - t;
    double Mz = SM.getMz(), cW2 = Mw*Mw/(Mz*Mz), sW2 = 1.0 - cW2;
    double ve = - 0.5 + 2.0*sW2, ae = -0.5, sigmae = ve + ae;
    double vfa = 0.5 - 2.0*fabs(Qf)*sW2;
    double mu = Mw; // renormalization scale
    
    double Qfp;
    if (Qf==-1.0)
        Qfp = 0.0;
    else if (Qf==0.0)
        Qfp = -1.0;
    else if (Qf==2.0/3.0)
        Qfp = -1.0/3.0;
    else if (Qf==-1.0/3.0)
        Qfp = 2.0/3.0;
    else 
        throw std::runtime_error("Error in EWSMTwoFermionsLEP2::kappa_e()");
    
    complex kappae = 1.0;
    kappae += SM.getAle()/4.0/M_PI/sW2    
              *( - cW2/sW2*DeltaRhobar(mu, Mw) + Pibar_Zgamma_hat(s, Mw)
                 - PV.B0(mu, s, Mw, Mw)/6.0 - 1.0/9.0 
                 - ve*sigmae/2.0/cW2*F_za_0(s, Mw) - F_W_0_hat(s, Mw)
                 + ( Mz*Mz/s - 1.0 )
                   *( fabs(Qf)*vfa*F_za_0(s, Mw)
                      + cW2*(F_Wn_0_hat(s, Mw) - fabs(Qfp)*F_Wa_0(s, Mw)) ) );

    // WW box    
    if (bWWbox)
        kappae += Delta_kappa_e_WW_hat(s, t, u, Mw, I3f);
    
    // ZZ box
    if (bZZbox)
        kappae += Delta_kappa_e_ZZ(mu, s, t, u, Mw, I3f, Qf);
    
    // Top quark contribution for e^+ e^- -> b bbar
    if (I3f==SM.getQuarks(SM.BOTTOM).getIsospin() 
            && Qf==SM.getQuarks(SM.BOTTOM).getCharge() 
            && mfp!=0.0)
        kappae += Delta_kappa_e_TOP(s, t, u, Mw, bWWbox);

    return kappae;
}


complex EWSMTwoFermionsLEP2::kappa_f(const double s, const double t, 
                                     const double Mw, const double I3f, 
                                     const double Qf, const double mf, 
                                     const double mfp, const bool bWWbox, 
                                     const bool bZZbox) const {
    double u = 2.0*mf*mf - s - t;
    double Mz = SM.getMz(), cW2 = Mw*Mw/(Mz*Mz), sW2 = 1.0 - cW2;
    double vea = 0.5 - 2.0*sW2;
    double vf = I3f - 2.0*Qf*sW2, af = I3f, sigmaf = vf + af;
    double mu = Mw; // renormalization scale

    complex kappaf = 1.0;
    kappaf += SM.getAle()/4.0/M_PI/sW2    
              *( - cW2/sW2*DeltaRhobar(mu, Mw) + Pibar_Zgamma_hat(s, Mw)
                 - PV.B0(mu, s, Mw, Mw)/6.0 - 1.0/9.0 
                 - vf*sigmaf/2.0/cW2*F_za_0(s, Mw) - F_W_0_hat(s, Mw)
                 + ( Mz*Mz/s - 1.0 )
                   *( vea*F_za_0(s, Mw) + cW2*F_Wn_0_hat(s, Mw) ) );

    // WW box    
    if (bWWbox)
        kappaf += Delta_kappa_f_WW_hat(s, t, u, Mw, I3f);
    
    // ZZ box
    if (bZZbox)
        kappaf += Delta_kappa_f_ZZ(mu, s, t, u, Mw, I3f, Qf);
    
    // Top quark contribution for e^+ e^- -> b bbar
    if (I3f==SM.getQuarks(SM.BOTTOM).getIsospin() 
            && Qf==SM.getQuarks(SM.BOTTOM).getCharge() 
            && mfp!=0.0)
        kappaf += Delta_kappa_f_TOP(s, t, u, Mw, bWWbox);
    
    return kappaf;
}


complex EWSMTwoFermionsLEP2::kappa_ef(const double s, const double t, 
                                      const double Mw, const double I3f, 
                                      const double Qf, const double mf, 
                                      const double mfp, const bool bWWbox, 
                                      const bool bZZbox) const {
    double u = 2.0*mf*mf - s - t;
    double Mz = SM.getMz(), cW2 = Mw*Mw/(Mz*Mz), sW2 = 1.0 - cW2;
    double ve = - 0.5 + 2.0*sW2, ae = -0.5, deltae = ve - ae;
    double vf = I3f - 2.0*Qf*sW2, af = I3f, deltaf = vf - af;
    double mu = Mw; // renormalization scale

    complex kappaef = 1.0;
    kappaef += SM.getAle()/4.0/M_PI/sW2    
              *( - 2.0*cW2/sW2*DeltaRhobar(mu, Mw) + 2.0*Pibar_Zgamma_hat(s, Mw)
                 - PV.B0(mu, s, Mw, Mw)/3.0 - 2.0/9.0 
                 - ((deltae*deltae + deltaf*deltaf)/sW2*(Mw*Mw/s - 1.0)
                    + 3.0*ve*ve + ae*ae + 3.0*vf*vf + af*af)*F_za_0(s, Mw)/4.0/cW2
                 - 2.0*F_W_0_hat(s, Mw)
                 + (Mz*Mz/s - 1.0)*cW2*(2.0/3.0 + Pibar_gg_bos_hat(s, Mw)) );

    // WW box    
    if (bWWbox) 
        kappaef += Delta_kappa_ef_WW_hat(s, t, u, Mw, I3f);
    
    // ZZ box
    if (bZZbox)
        kappaef += Delta_kappa_ef_ZZ(mu, s, t, u, Mw, I3f, Qf);
    
    // Top quark contribution for e^+ e^- -> b bbar
    if (I3f==SM.getQuarks(SM.BOTTOM).getIsospin() 
            && Qf==SM.getQuarks(SM.BOTTOM).getCharge() 
            && mfp!=0.0)
        kappaef += Delta_kappa_ef_TOP(s, t, u, Mw, bWWbox);
    
    return kappaef;
}


complex EWSMTwoFermionsLEP2::Delta_rho_ef_TOP(const double s, const double t, 
                                              const double u, const double Mw, 
                                              const bool bWWbox) const {
    double Mz = SM.getMz(), cW2 = Mw*Mw/(Mz*Mz), sW2 = 1.0 - cW2;
    double Mt = SM.getMtpole();
    double mu = Mw; // renormalization scale
    
    complex Bww = complex(0.0, 0.0, false);
    if (bWWbox)
        Bww = Delta_rho_ef_WW_TOP_hat(s, t, u, Mw);
        
    return ( SM.getAle()/4.0/M_PI/sW2
             *( F_W_t_hat(s, Mw) 
                - Mt*Mt/4.0/Mw/Mw*(PV.B0(mu, s, Mw, Mw) + 1.0) ) + Bww ); 
}


complex EWSMTwoFermionsLEP2::Delta_kappa_e_TOP(const double s, const double t, 
                                               const double u, const double Mw, 
                                               const bool bWWbox) const {
    double Mz = SM.getMz(), cW2 = Mw*Mw/(Mz*Mz), sW2 = 1.0 - cW2;
    double Qfp = SM.getQuarks(SM.TOP).getCharge();
    
    complex Bww = complex(0.0, 0.0, false);
    if (bWWbox)
        Bww = Delta_kappa_e_WW_TOP_hat(s, t, u, Mw);

    return ( SM.getAle()/4.0/M_PI/sW2    
             *( Mz*Mz/s - 1.0 )*cW2
             *( F_Wn_t_hat(s, Mw) - fabs(Qfp)*F_Wa_t(s, Mw) ) + Bww );
}


complex EWSMTwoFermionsLEP2::Delta_kappa_f_TOP(const double s, const double t, 
                                               const double u, const double Mw, 
                                               const bool bWWbox) const {
    double Mz = SM.getMz(), cW2 = Mw*Mw/(Mz*Mz), sW2 = 1.0 - cW2;
    double Mt = SM.getMtpole();
    double mu = Mw; // renormalization scale

    complex Bww = complex(0.0, 0.0, false);
    if (bWWbox)
        Bww = Delta_kappa_f_WW_TOP_hat(s, t, u, Mw);
    
    return ( SM.getAle()/4.0/M_PI/sW2    
             *( - F_W_t_hat(s, Mw) 
                + Mt*Mt/4.0/Mw/Mw*( PV.B0(mu, s, Mw, Mw) + 1.0) ) + Bww );
}


complex EWSMTwoFermionsLEP2::Delta_kappa_ef_TOP(const double s, const double t, 
                                                const double u, const double Mw, 
                                                const bool bWWbox) const {
    double Mz = SM.getMz(), cW2 = Mw*Mw/(Mz*Mz), sW2 = 1.0 - cW2;
    double Mt = SM.getMtpole();
    double mu = Mw; // renormalization scale

    complex Bww = complex(0.0, 0.0, false);
    if (bWWbox)
        Bww = Delta_kappa_ef_WW_TOP_hat(s, t, u, Mw);
    
    return ( SM.getAle()/4.0/M_PI/sW2    
             *( - F_W_t_hat(s, Mw) 
                + Mt*Mt/4.0/Mw/Mw*( PV.B0(mu, s, Mw, Mw) + 1.0) ) + Bww );
}
    

complex EWSMTwoFermionsLEP2::Delta_rho_ef_WW_hat(const double s, const double t, 
                                                 const double u, const double Mw,
                                                 const double I3f) const {
    double Mz = SM.getMz(), cW2 = Mw*Mw/(Mz*Mz), sW2 = 1.0 - cW2;

    if (I3f==SM.getLeptons(SM.ELECTRON).getIsospin() 
            || I3f==SM.getQuarks(SM.DOWN).getIsospin()) 
        return ( SM.getAle()/4.0/M_PI/sW2
                 *( - cW2*(Mz*Mz - s)*B_WW_d_0_hat(s, t, u, Mw) ) );
    else if (I3f==SM.getLeptons(SM.NEUTRINO_1).getIsospin() 
                || I3f==SM.getQuarks(SM.UP).getIsospin()) 
        return ( SM.getAle()/4.0/M_PI/sW2
                 *( cW2*(Mz*Mz - s)*B_WW_c_0_hat(s, t, u, Mw) ) );
    else 
        throw std::runtime_error("Error in EWSMTwoFermionsLEP2::Delta_rho_ef_WW_hat()");
}


complex EWSMTwoFermionsLEP2::Delta_kappa_e_WW_hat(const double s, const double t, 
                                                  const double u, const double Mw,
                                                  const double I3f) const {
    double Mz = SM.getMz(), cW2 = Mw*Mw/(Mz*Mz), sW2 = 1.0 - cW2;

    if (I3f==SM.getLeptons(SM.ELECTRON).getIsospin() 
            || I3f==SM.getQuarks(SM.DOWN).getIsospin()) 
        return ( SM.getAle()/4.0/M_PI/sW2
                 *( cW2*(Mz*Mz - s)*B_WW_d_0_hat(s, t, u, Mw) ) );
    else if (I3f==SM.getLeptons(SM.NEUTRINO_1).getIsospin() 
                || I3f==SM.getQuarks(SM.UP).getIsospin()) 
        return ( SM.getAle()/4.0/M_PI/sW2
                 *( - cW2*(Mz*Mz - s)*B_WW_c_0_hat(s, t, u, Mw) ) );
    else 
        throw std::runtime_error("Error in EWSMTwoFermionsLEP2::Delta_kappa_e_WW_hat()");    
}


complex EWSMTwoFermionsLEP2::Delta_kappa_f_WW_hat(const double s, const double t, 
                                                  const double u, const double Mw,
                                                  const double I3f) const {
    double Mz = SM.getMz(), cW2 = Mw*Mw/(Mz*Mz), sW2 = 1.0 - cW2;
       
    if (I3f==SM.getLeptons(SM.ELECTRON).getIsospin() 
            || I3f==SM.getQuarks(SM.DOWN).getIsospin()) 
        return ( SM.getAle()/4.0/M_PI/sW2
                 *( cW2*(Mz*Mz - s)*B_WW_d_0_hat(s, t, u, Mw) ) );
    else if (I3f==SM.getLeptons(SM.NEUTRINO_1).getIsospin() 
                || I3f==SM.getQuarks(SM.UP).getIsospin())  
        return ( SM.getAle()/4.0/M_PI/sW2
                 *( - cW2*(Mz*Mz - s)*B_WW_c_0_hat(s, t, u, Mw) ) );
    else 
        throw std::runtime_error("Error in EWSMTwoFermionsLEP2::Delta_kappa_f_WW_hat()");
}


complex EWSMTwoFermionsLEP2::Delta_kappa_ef_WW_hat(const double s, const double t, 
                                                   const double u, const double Mw,
                                                   const double I3f) const {
    double Mz = SM.getMz(), cW2 = Mw*Mw/(Mz*Mz), sW2 = 1.0 - cW2;

    if (I3f==SM.getLeptons(SM.ELECTRON).getIsospin() 
            || I3f==SM.getQuarks(SM.DOWN).getIsospin()) 
        return ( SM.getAle()/4.0/M_PI/sW2
                 *( cW2*(Mz*Mz - s)*B_WW_d_0_hat(s, t, u, Mw) ) );
    else if (I3f==SM.getLeptons(SM.NEUTRINO_1).getIsospin() 
                || I3f==SM.getQuarks(SM.UP).getIsospin()) 
        return ( SM.getAle()/4.0/M_PI/sW2
                 *( - cW2*(Mz*Mz - s)*B_WW_c_0_hat(s, t, u, Mw) ) );
    else 
        throw std::runtime_error("Error in EWSMTwoFermionsLEP2::Delta_kappa_ef_WW_hat()");
}


complex EWSMTwoFermionsLEP2::Delta_rho_ef_WW_TOP_hat(const double s, const double t, 
                                                     const double u, const double Mw) const {
    double Mz = SM.getMz(), cW2 = Mw*Mw/(Mz*Mz), sW2 = 1.0 - cW2;
    return ( SM.getAle()/4.0/M_PI/sW2
             *( - cW2*(Mz*Mz - s)*Delta_B_WW_d_hat(s, t, u, Mw) ) );    
}


complex EWSMTwoFermionsLEP2::Delta_kappa_e_WW_TOP_hat(const double s, const double t, 
                                                      const double u, const double Mw) const {
    double Mz = SM.getMz(), cW2 = Mw*Mw/(Mz*Mz), sW2 = 1.0 - cW2;
    return ( SM.getAle()/4.0/M_PI/sW2
             *( cW2*(Mz*Mz - s)*Delta_B_WW_d_hat(s, t, u, Mw) ) );
}


complex EWSMTwoFermionsLEP2::Delta_kappa_f_WW_TOP_hat(const double s, const double t, 
                                                      const double u, const double Mw) const {
    double Mz = SM.getMz(), cW2 = Mw*Mw/(Mz*Mz), sW2 = 1.0 - cW2;
    return ( SM.getAle()/4.0/M_PI/sW2
             *( cW2*(Mz*Mz - s)*Delta_B_WW_d_hat(s, t, u, Mw) ) );    
}


complex EWSMTwoFermionsLEP2::Delta_kappa_ef_WW_TOP_hat(const double s, const double t, 
                                                       const double u, const double Mw) const {
    double Mz = SM.getMz(), cW2 = Mw*Mw/(Mz*Mz), sW2 = 1.0 - cW2;
    return ( SM.getAle()/4.0/M_PI/sW2
             *( cW2*(Mz*Mz - s)*Delta_B_WW_d_hat(s, t, u, Mw) ) );    
}


complex EWSMTwoFermionsLEP2::Delta_rho_ef_ZZ(const double mu, const double s, 
                                             const double t, const double u, 
                                             const double Mw, const double I3f, 
                                             const double Qf) const {
    double Mz = SM.getMz(), cW2 = Mw*Mw/(Mz*Mz), sW2 = 1.0 - cW2;
    double I3e = -0.5;
    double ve = - 0.5 + 2.0*sW2, ae = -0.5;
    double vf = I3f - 2.0*Qf*sW2, af = I3f;

    return ( SM.getAle()/4.0/M_PI/sW2*(s - Mz*Mz)/2.0/cW2
             *( (4.0*I3e*I3f*(ve*ve + ae*ae)*(vf*vf + af*af) + ve*vf)
                  *B_ZZ_0(mu, s, t, u)
                + (4.0*I3e*I3f*(ve*ve + ae*ae)*(vf*vf + af*af) - ve*vf)
                  *B_ZZ_0(mu, s, u, t) ) );
}
    

complex EWSMTwoFermionsLEP2::Delta_kappa_e_ZZ(const double mu, const double s, 
                                              const double t, const double u, 
                                              const double Mw, const double I3f,
                                              const double Qf) const {
    double Mz = SM.getMz(), cW2 = Mw*Mw/(Mz*Mz), sW2 = 1.0 - cW2;
    double ve = - 0.5 + 2.0*sW2, ae = -0.5, deltae = ve - ae;
    double vf = I3f - 2.0*Qf*sW2, af = I3f, sigmaf = vf + af, deltaf = vf - af;
    
    return ( - SM.getAle()/4.0/M_PI/sW2*(s - Mz*Mz)/2.0/cW2*deltae*I3f
             *( deltaf*deltaf*B_ZZ_0(mu, s, t, u)
                + sigmaf*sigmaf*B_ZZ_0(mu, s, u, t) ) 
             - Delta_rho_ef_ZZ(mu,s,t,u,Mw,I3f,Qf) );    
}


complex EWSMTwoFermionsLEP2::Delta_kappa_f_ZZ(const double mu, const double s, 
                                              const double t, const double u, 
                                              const double Mw, const double I3f,
                                              const double Qf) const {
    double Mz = SM.getMz(), cW2 = Mw*Mw/(Mz*Mz), sW2 = 1.0 - cW2;
    double I3e = -0.5;
    double ve = - 0.5 + 2.0*sW2, ae = -0.5, sigmae = ve + ae, deltae = ve - ae;
    double vf = I3f - 2.0*Qf*sW2, af = I3f, deltaf = vf - af;
    
    return ( - SM.getAle()/4.0/M_PI/sW2*(s - Mz*Mz)/2.0/cW2*deltaf*I3e
             *( deltae*deltae*B_ZZ_0(mu, s, t, u)
                + sigmae*sigmae*B_ZZ_0(mu, s, u, t) ) 
             - Delta_rho_ef_ZZ(mu,s,t,u,Mw,I3f,Qf) );  
}


complex EWSMTwoFermionsLEP2::Delta_kappa_ef_ZZ(const double mu, const double s, 
                                               const double t, const double u, 
                                               const double Mw, const double I3f,
                                               const double Qf) const {
    double Mz = SM.getMz(), cW2 = Mw*Mw/(Mz*Mz), sW2 = 1.0 - cW2;
    double ve = - 0.5 + 2.0*sW2, ae = -0.5, deltae = ve - ae;
    double vf = I3f - 2.0*Qf*sW2, af = I3f, deltaf = vf - af;
    
    return ( SM.getAle()/4.0/M_PI/sW2*(s - Mz*Mz)/2.0/cW2*deltae*deltaf
             *B_ZZ_0(mu, s, t, u)
             - Delta_rho_ef_ZZ(mu,s,t,u,Mw,I3f,Qf) );      
}


////////////////////////////////////////////////////////////////////////

complex EWSMTwoFermionsLEP2::I2e(const double s, const double Mw) const {
    double Mz = SM.getMz(), sW2 = 1.0 - Mw*Mw/(Mz*Mz);
    double alpha = SM.getAle()*V_pol(s).real();
    double ReKappa_e = 1.0;
    return ( 35.0*alpha*alpha/18.0*( 1.0 - 8.0/3.0*ReKappa_e*sW2 ) );
}


complex EWSMTwoFermionsLEP2::I2f(const double s, const double Mw) const {
    double Mz = SM.getMz(), sW2 = 1.0 - Mw*Mw/(Mz*Mz);
    double alpha = SM.getAle()*V_pol(s).real();
    double ReKappa_f = 1.0;
    return ( 35.0*alpha*alpha/18.0*( 1.0 - 8.0/3.0*ReKappa_f*sW2 ) );
}


complex EWSMTwoFermionsLEP2::DeltaRhobar(const double mu, const double Mw) const {
    return myOneLoopEW.DeltaRhobar(mu,Mw);
}


complex EWSMTwoFermionsLEP2::DeltaRhobarZ(const double mu, const double Mw) const {
    return ( myOneLoopEW.DeltaRhobar(mu,Mw) + myOneLoopEW.DeltaRhobarW(mu,Mw) );
}


complex EWSMTwoFermionsLEP2::D_Z(const double mu, const double s, 
                                 const double Mw) const {
    double Mz = SM.getMz(), cW2 = Mw*Mw/(Mz*Mz);
    complex D_Z_bos = (myOneLoopEW.SigmaZZ_bos(mu,s,Mw) 
                       - myOneLoopEW.SigmaZZ_bos(mu,Mz*Mz,Mw))/cW2/(Mz*Mz - s);
    complex D_Z_fer = (myOneLoopEW.SigmaZZ_fer(mu,s,Mw) 
                       - myOneLoopEW.SigmaZZ_fer(mu,Mz*Mz,Mw))/cW2/(Mz*Mz - s);
    return ( D_Z_bos + D_Z_fer );    
}


complex EWSMTwoFermionsLEP2::Pibar_Zgamma(const double mu, const double s, 
                                          const double Mw) const {
    complex Pibar_Zgamma_bos = myOneLoopEW.PiZgamma_bos(mu,s,Mw);
    complex Pibar_Zgamma_fer = myOneLoopEW.PiZgamma_fer(mu,s,Mw);
    return ( Pibar_Zgamma_bos + Pibar_Zgamma_fer );    
}


complex EWSMTwoFermionsLEP2::Pibar_gg_bos(const double mu, const double s, 
                                          const double Mw) const {
    complex Pibar_gg_bos = myOneLoopEW.PiGammaGamma_bos(mu,s,Mw);
    return Pibar_gg_bos;  
}


complex EWSMTwoFermionsLEP2::F_za_0(const double s, const double Mw) const {
    return myOneLoopEW.FZa_0(s, Mw);
}


complex EWSMTwoFermionsLEP2::F_Wa_0(const double s, const double Mw) const {
    return myOneLoopEW.FWa_0(s, Mw);  
} 


complex EWSMTwoFermionsLEP2::F_Wa_t(const double s, const double Mw) const {
    return myOneLoopEW.FWa_t(s, Mw);    
}


complex EWSMTwoFermionsLEP2::F_Wn_0(const double s, const double Mw) const {
    return myOneLoopEW.FWn_0(s, Mw);  
}


complex EWSMTwoFermionsLEP2::F_Wn_t(const double s, const double Mw) const {
    return myOneLoopEW.FWn_t(s, Mw);
}

    
complex EWSMTwoFermionsLEP2::F_W_0(const double s, const double Mw) const {
    double Mz = SM.getMz(), cW2 = Mw*Mw/(Mz*Mz);
    return ( cW2*F_Wn_0(s, Mw) - F_Wa_0(s, Mw)/2.0 );
}


complex EWSMTwoFermionsLEP2::F_W_t(const double s, const double Mw) const {
    double Mz = SM.getMz(), cW2 = Mw*Mw/(Mz*Mz);
    return ( cW2*F_Wn_t(s, Mw) - F_Wa_t(s, Mw)/2.0 
             - myOneLoopEW.FbarWa_t(s, Mw)/2.0 );        
}


complex EWSMTwoFermionsLEP2::B_WW_d_0(const double mu, const double s, 
                                      const double t, const double u, 
                                      const double Mw) const {
    double s2 = s*s, t2 = t*t, u2 = u*u, Mw2 = Mw*Mw;
    return ( ( - t*(1.0 + t2/u2) - 4.0*Mw2*t2/u2 + 2.0*Mw2*Mw2/u*(1.0 + 2.0*s/u) )
              *PV.D0(s, t, Mw, 0.0, Mw, 0.0)
             - 2.0*( 2.0 + 2.0*s/u + s2/u2 - 2.0*Mw2*s/u2 )*PV.C0(s, Mw, 0.0, Mw)
             + 2.0*( 2.0 + 3.0*s/u + s2/u2 + 2.0*Mw2*t/u2 )*PV.C0(t, 0.0, Mw, 0.0)
             + ( - 2.0/u - 5.0/3.0/Mw2 - s/12.0/Mw2/Mw2 )*PV.B0(mu, s, Mw, Mw)
             + 2.0/u*PV.B0(mu, t, 0.0, 0.0) - 1.0/6.0/Mw2/Mw2*PV.A0(mu, Mw)
             + 1.0/3.0/Mw2 - s/18.0/Mw2/Mw2 );
}


complex EWSMTwoFermionsLEP2::B_WW_d(const double mu, const double s, 
                                    const double t, const double u, 
                                    const double Mw) const {
    double Mt = SM.getMtpole(), Mt2 = Mt*Mt;
    double s2 = s*s, t2 = t*t, u2 = u*u, Mw2 = Mw*Mw;
    return ( ( - t*(1.0 + t2/u2) - 4.0*Mw2*t2/u2 + 2.0*Mw2*Mw2/u*(1.0 + 2.0*s/u) 
               + Mt2*(2.0 + 3.0*s/u + 2.0*s2/u2 - 2.0*Mw2/u*(1.0 + 2.0*s/u) ) 
               + Mt2*Mt2*s/u2 )*PV.D0(s, t, Mw, 0.0, Mw, Mt)
             + ( - 2.0 - 2.0*s/u - s2/u2 + 2.0*Mw2*s/u2 
                 + Mt2/2.0/Mw2*(4.0 - Mw2/s*(1.0 + 2.0*s2/u2)) 
                 - Mt2*Mt2/2.0/Mw2/Mw2*(1.0 - 2.0*Mw2/s) 
                 - Mt2*Mt2*Mt2/2.0/Mw2/Mw2/s )*PV.C0(s, Mw, Mt, Mw)
             - ( 2.0 + 2.0*s/u + s2/u2 - 2.0*Mw2*s/u2 + Mt2*s/u2)*PV.C0(s, Mw, 0.0, Mw)
             + ( 2.0 + 3.0*s/u + s2/u2 + 2.0*Mw2*t/u2 - Mt2*t/u2)
               *(PV.C0(t, 0.0, Mw, Mt) + PV.C0(t, Mt, Mw, 0.0))
             + ( - 2.0/u - 5.0/3.0/Mw2 - s/12.0/Mw2/Mw2 
                 - Mt2/4.0/Mw2/s*(2.0 - s/Mw2) + Mt2*Mt2/2.0/Mw2/Mw2/s )        
               *PV.B0(mu, s, Mw, Mw)
             + 2.0/u*PV.B0(mu, t, Mt, 0.0) + Mt2/2.0/Mw2/Mw2/s*PV.A0(mu,Mt)
             - 1.0/6.0/Mw2/Mw2*(1.0 + 3.0*Mt2/s)*PV.A0(mu, Mw)
             + 1.0/3.0/Mw2*(1.0 + 3.0*Mt2/4.0/Mw2 - s/6.0/Mw2) );
}


complex EWSMTwoFermionsLEP2::Delta_B_WW_d(const double mu, const double s, 
                                          const double t, const double u, 
                                          const double Mw) const {
    return ( B_WW_d(mu, s, t, u, Mw) - B_WW_d_0(mu, s, t, u, Mw) );
}


complex EWSMTwoFermionsLEP2::B_WW_c_0(const double mu, const double s, 
                                      const double t, const double u, 
                                      const double Mw) const {
    double Mw2 = Mw*Mw;
    return ( 2.0*u*PV.D0(s, u, Mw, 0.0, Mw, 0.0)
             + 4.0*PV.C0(s, Mw, 0.0, Mw) 
             + (20.0 + s/Mw2)/12.0/Mw2*PV.B0(mu, s, Mw, Mw)
             + PV.A0(mu, Mw)/6.0/Mw2/Mw2 - (1.0 - s/6.0/Mw2)/3.0/Mw2 );
}


complex EWSMTwoFermionsLEP2::B_ZZ_0(const double mu, const double s,
                                    const double t, const double u) const {
    double Mz = SM.getMz(), Mz2 = Mz*Mz;
    double t2 = t*t, u2 = u*u;
    return ( 2.0*u*PV.D0(s, u, Mz, 0.0, Mz, 0.0)
             + ( - 2.0*u - t*(3.0 + t/u)*(3.0 + t/u) 
                 + 2.0*(Mz2 - s)*(1.0 + 3.0*t/u - Mz2/u*(1.0 + 2.0*t/u)) )
               *PV.D0(s, t, Mz, 0.0, Mz, 0.0)
             + 2.0*( 3.0 + 4.0*t/u + t2/u2 - 2.0*s*(s - Mz2)/u2 )
               *PV.C0(s, Mz, 0.0, Mz)
             - 2.0*t/u*( 3.0 + t/u + 2.0*(s - Mz2)/u )*PV.C0(t, 0.0, Mz, 0.0)
             - 2.0/u*( PV.B0(mu, s, Mz, Mz) - PV.B0(mu, t, 0.0, 0.0) ) );
}


////////////////////////////////////////////////////////////////////////

complex EWSMTwoFermionsLEP2::Pibar_Zgamma_hat(const double s, const double Mw) const {
    double mu = Mw;

    complex add = complex(0.0, 0.0, false);
    if (!bKeepNonUnitary) {
        double Mz = SM.getMz(), cW2 = Mw*Mw/(Mz*Mz);
        double Rw = Mw*Mw/s;
        add = - cW2*( (4.0/3.0/Rw + 1.0/12.0/Rw/Rw)*PV.B0(mu,s,Mw,Mw)
                      + 1.0/18.0/Rw/Rw - 13.0/18.0/Rw );
    }
    
    return ( Pibar_Zgamma(mu, s, Mw) + add );
}


complex EWSMTwoFermionsLEP2::Pibar_gg_bos_hat(const double s, const double Mw) const {
    double mu = Mw;
    
    complex add = complex(0.0, 0.0, false);
    if (!bKeepNonUnitary) {
        double Mz = SM.getMz(), cW2 = Mw*Mw/(Mz*Mz);
        double Rw = Mw*Mw/s;
        add = s/(s - Mz*Mz)
                *( 1.0/12.0/Rw/cW2*PV.B0(mu,s,Mw,Mw) + 1.0/18.0/Rw/cW2 )
              - s/(s - Mz*Mz)
                *( (4.0/3.0/Rw + 1.0/12.0/Rw/Rw)*PV.B0(mu,s,Mw,Mw)
                   + 1.0/18.0/Rw/Rw - 13.0/18.0/Rw );
    }
    
    return ( Pibar_gg_bos(mu, s, Mw) + add );
}


complex EWSMTwoFermionsLEP2::D_Z_hat(const double s, const double Mw) const {
    double mu = Mw;
    
    complex add = complex(0.0, 0.0, false);
    if (!bKeepNonUnitary) {
        double Mz = SM.getMz(), cW2 = Mw*Mw/(Mz*Mz);
        double Rz = Mz*Mz/s;
        add = ( (1.0/12.0/cW2 + 4.0/3.0)/Rz + 1.0/12.0/cW2/Rz/Rz )*PV.B0(mu,s,Mw,Mw)
              + ( (1.0/cW2 - 13.0)/Rz + 1.0/cW2/Rz/Rz )/18.0;
    }
  
    return ( D_Z(mu, s, Mw) + add );
}


complex EWSMTwoFermionsLEP2::F_Wn_0_hat(const double s, const double Mw) const {
    complex add = complex(0.0, 0.0, false);
    if (!bKeepNonUnitary) {
        double Mz = SM.getMz(), cW2 = Mw*Mw/(Mz*Mz);
        double Rw = Mw*Mw/s;
        double mu = Mw;
        add = - s/(s-Mz*Mz)*( - 1.0/12.0/Rw/cW2 + 3.0/2.0/Rw + 1.0/12.0/Rw/Rw )
                *PV.B0(mu,s,Mw,Mw) 
              - s/(s-Mz*Mz)*( -1.0/18.0/Rw/cW2 - 11.0/18.0/Rw + 1.0/18.0/Rw/Rw );
    }
  
    return ( F_Wn_0(s, Mw) + add );
}   


complex EWSMTwoFermionsLEP2::F_Wn_t_hat(const double s, const double Mw) const {
    complex add = complex(0.0, 0.0, false);
    if (!bKeepNonUnitary) {
        double Mz = SM.getMz(), Rw = Mw*Mw/s, Mt = SM.getMtpole();
        double mu = Mw;
        add = Mt*Mt*s/4.0/Rw/Mw/Mw/(s-Mz*Mz)*(PV.B0(mu,s,Mw,Mw) + 1.0);
    }
  
    return ( F_Wn_t(s, Mw) + add );
}   


complex EWSMTwoFermionsLEP2::F_W_0_hat(const double s, const double Mw) const {
    double Mz = SM.getMz(), cW2 = Mw*Mw/(Mz*Mz);
    complex F_W_0 = cW2*F_Wn_0(s, Mw) - F_Wa_0(s, Mw)/2.0;

    complex add = complex(0.0, 0.0, false);
    if (!bKeepNonUnitary) {
        double Rw = Mw*Mw/s;
        double mu = Mw;
        add = cW2*( - 3.0/2.0/Rw - 1.0/12.0/Rw/Rw )*PV.B0(mu,s,Mw,Mw) 
              + cW2*( 11.0/18.0/Rw - 1.0/18.0/Rw/Rw );
    }
    
    return ( F_W_0 + add );
}


complex EWSMTwoFermionsLEP2::F_W_t_hat(const double s, const double Mw) const {
    double Mz = SM.getMz(), cW2 = Mw*Mw/(Mz*Mz);
    complex F_W_t = cW2*F_Wn_t(s, Mw) - F_Wa_t(s, Mw)/2.0 
                    - myOneLoopEW.FbarWa_t(s, Mw)/2.0;

    complex add = complex(0.0, 0.0, false);
    if (!bKeepNonUnitary) {
        double Mt = SM.getMtpole();
        double Rw = Mw*Mw/s;
        double mu = Mw;
        add = Mt*Mt/4.0/Rw/Mz/Mz*(PV.B0(mu,s,Mw,Mw) + 1.0);
    }
    
    return ( F_W_t + add );
}

complex EWSMTwoFermionsLEP2::B_WW_d_0_hat(const double s, const double t, 
                                          const double u, const double Mw) const {
    double mu = Mw;
    
    complex add = complex(0.0, 0.0, false);
    if (!bKeepNonUnitary) {
        double Mz = SM.getMz(), cW2 = Mw*Mw/(Mz*Mz), Rw = Mw*Mw/s;
        add = (5.0/3.0 - 1.0/12.0/cW2 + 1.0/12.0/Rw)/Rw/(s - Mz*Mz)*PV.B0(mu,s,Mw,Mw)
               - (1.0/2.0 + 1.0/18.0/cW2 - 1.0/18.0/Rw)/Rw/(s - Mz*Mz);
    }
  
    return ( B_WW_d_0(mu, s, t, u, Mw) + add );
}

complex EWSMTwoFermionsLEP2::B_WW_d_0_hat_TEST(const double s, const double t, 
                                               const double u, const double Mw) const {
    double mu = Mw;
    
    if (!bKeepNonUnitary) {
        double s2 = s*s, t2 = t*t, u2 = u*u, Mw2 = Mw*Mw;
        return ( ( - t*(1.0 + t2/u2) - 4.0*Mw2*t2/u2 + 2.0*Mw2*Mw2/u*(1.0 + 2.0*s/u) )
                  *PV.D0(s, t, Mw, 0.0, Mw, 0.0)
                 - 2.0*( 2.0 + 2.0*s/u + s2/u2 - 2.0*Mw2*s/u2 )*PV.C0(s, Mw, 0.0, Mw)
                 + 2.0*( 2.0 + 3.0*s/u + s2/u2 + 2.0*Mw2*t/u2 )*PV.C0(t, 0.0, Mw, 0.0)
                 - 2.0/u*( PV.B0(mu, s, Mw, Mw) - PV.B0(mu, t, 0.0, 0.0) ) );
    } else {
        return B_WW_d_0(mu, s, t, u, Mw);
    }
}


complex EWSMTwoFermionsLEP2::Delta_B_WW_d_hat(const double s, const double t, 
                                              const double u, const double Mw) const {
    double mu = Mw;
    
    complex add = complex(0.0, 0.0, false);
    if (!bKeepNonUnitary) {
        double Mz = SM.getMz(), Rw = Mw*Mw/s, Mt = SM.getMtpole();
        add = - Mt*Mt/4.0/Rw/Mw/Mw/(s - Mz*Mz)*( PV.B0(mu, s, Mw, Mw) + 1.0 );
    }
  
    return ( Delta_B_WW_d(mu, s, t, u, Mw) + add );
}


complex EWSMTwoFermionsLEP2::B_WW_c_0_hat(const double s, const double t, 
                                          const double u, const double Mw) const {
    double mu = Mw;
    
    complex add = complex(0.0, 0.0, false);
    if (!bKeepNonUnitary) {
        double Mz = SM.getMz(), cW2 = Mw*Mw/(Mz*Mz), Rw = Mw*Mw/s;
        add = - (5.0/3.0 - 1.0/12.0/cW2 + 1.0/12.0/Rw)/Rw/(s - Mz*Mz)*PV.B0(mu,s,Mw,Mw)
              + (1.0/2.0 + 1.0/18.0/cW2 - 1.0/18.0/Rw)/Rw/(s - Mz*Mz);
    }
  
    return ( B_WW_c_0(mu, s, t, u, Mw) + add );
}



