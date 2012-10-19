/* 
 * File:   EWSMTwoFermionsLEP2.cpp
 * Author: mishima
 */

#include <stdexcept>
#include <cmath>
#include "EWSMTwoFermionsLEP2.h"
#include "EWSMThreeLoopEW.h"


EWSMTwoFermionsLEP2::EWSMTwoFermionsLEP2(const StandardModel& SM_i) : SM(SM_i), 
        myCache(SM_i), myOneLoopEW(myCache) {
    bDebug = SM_i.isBDebug();
}


//////////////////////////////////////////////////////////////////////// 

complex EWSMTwoFermionsLEP2::Vpol_inv(const double s) const {
    complex V_inv;
    if (bDebug)
        V_inv = 1.0/complex(1.0715119759, -0.0186242179, false); // for debug
    else {
        //!!!!!!
        V_inv = SM.getAle()/SM.alphaMz(); //!!TEST
        //V_inv = 1.0 - SM.DeltaAlphaLepton(s) - SM.getDAle5Mz(); // !!TEST
        //V_inv = 1.0 - myOneLoopEW.DeltaAlpha_l(s) - SM.getDAle5Mz()
        //        - myOneLoopEW.DeltaAlpha_t(s);
        //V_inv = 1.0/complex(1.0715119759, 0.0, false); //!!TEST
        //V_inv = 1.0/complex(1.0715119759, -0.0186242179, false); //!!TEST
        //!!!!!!
    }
    return V_inv;    
}


complex EWSMTwoFermionsLEP2::chi_Z(const double s, const double Mw, 
                                   const double GammaZ) const {
    double Mz = SM.getMz();
    complex denom = complex(s - Mz*Mz, GammaZ/Mz*s, false);
    double prefactor = SM.getGF()*Mz*Mz/(sqrt(2.0)*8.0*M_PI*SM.getAle());
    
    return ( prefactor*s/denom );
}

    
complex EWSMTwoFermionsLEP2::rho_ef(const double s, const double Mw, 
                                    const double I3f, const double Qf, 
                                    const double mfp) const {
    double Mz = SM.getMz(), cW2 = Mw*Mw/(Mz*Mz), sW2 = 1.0 - cW2;
    double ve = - 0.5 + 2.0*sW2, ae = -0.5;
    double vf = I3f - 2.0*Qf*sW2, af = I3f;
    double mu = Mw; // renormalization scale
    
    complex FW_hat;
    if (mfp==SM.getMtpole())
        FW_hat = F_W_t_hat(mu, s, Mw);
    else 
        FW_hat = F_W_0_hat(mu, s, Mw);
    
    complex rhoef = 1.0;
    rhoef += SM.getAle()/4.0/M_PI/sW2
             *( - DeltaRhobarZ(mu, Mw) + D_Z_hat(mu, s, Mw) 
                + 5.0/3.0*PV.B0(mu, s, Mw, Mw) - 9.0*cW2/4.0/sW2*log(cW2) - 6.0
                + 5.0*cW2/8.0*(1.0 + cW2) 
                + (3.0*ve*ve + ae*ae + 3.0*vf*vf + af*af)/4.0/cW2*F_za_0(s, Mw)
                + F_W_0_hat(mu, s, Mw) + FW_hat
                - mfp*mfp/4.0/Mw/Mw*(PV.B0(mu, s, Mw, Mw) - 1.0) );
    return rhoef;       
}

complex EWSMTwoFermionsLEP2::kappa_e(const double s, const double Mw, 
                                     const double I3f, const double Qf, 
                                     const double mfp) const {
    double Mz = SM.getMz(), cW2 = Mw*Mw/(Mz*Mz), sW2 = 1.0 - cW2;
    double ve = - 0.5 + 2.0*sW2, ae = -0.5, sigmae = ve + ae;
    double vfa = 0.5 - 2.0*fabs(Qf)*sW2;
    double mu = Mw; // renormalization scale
    
    complex FWn_hat, FWa;
    if (mfp==SM.getMtpole()) {
        FWn_hat = F_Wn_t_hat(mu, s, Mw);
        FWa = F_Wa_t(s, Mw);
    } else {
        FWn_hat = F_Wn_0_hat(mu, s, Mw);     
        FWa = F_Wa_0(s, Mw);
    }
    
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
              *( - cW2/sW2*DeltaRhobar(mu, Mw) + Pibar_Zgamma_hat(mu, s, Mw)
                 - PV.B0(mu, s, Mw, Mw)/6.0 - 1.0/9.0 
                 - ve*sigmae/2.0/cW2*F_za_0(s, Mw) - F_W_0_hat(mu, s, Mw)
                 + ( Mz*Mz/s - 1.0 )
                   *( fabs(Qf)*vfa*F_za_0(s, Mw)
                      + cW2*(FWn_hat - fabs(Qfp)*FWa) ) );
    return kappae;
}


complex EWSMTwoFermionsLEP2::kappa_f(const double s, const double Mw, 
                                     const double I3f, const double Qf, 
                                     const double mfp) const {
    double Mz = SM.getMz(), cW2 = Mw*Mw/(Mz*Mz), sW2 = 1.0 - cW2;
    double vea = 0.5 - 2.0*sW2;
    double vf = I3f - 2.0*Qf*sW2, af = I3f, sigmaf = vf + af;
    double mu = Mw; // renormalization scale
    
    complex FW_hat;
    if (mfp==SM.getMtpole())
        FW_hat = F_W_t_hat(mu, s, Mw);
    else
        FW_hat = F_W_0_hat(mu, s, Mw); 

    complex kappaf = 1.0;
    kappaf += SM.getAle()/4.0/M_PI/sW2    
              *( - cW2/sW2*DeltaRhobar(mu, Mw) + Pibar_Zgamma_hat(mu, s, Mw)
                 - PV.B0(mu, s, Mw, Mw)/6.0 - 1.0/9.0 
                 - vf*sigmaf/2.0/cW2*F_za_0(s, Mw) - FW_hat
                 + ( Mz*Mz/s - 1.0 )
                   *( vea*F_za_0(s, Mw) + cW2*F_Wn_0_hat(mu, s, Mw))
                 + mfp*mfp/4.0/Mw/Mw*( PV.B0(mu, s, Mw, Mw) + 1.0) );
    return kappaf;
}


complex EWSMTwoFermionsLEP2::kappa_ef(const double s, const double Mw, 
                                      const double I3f, const double Qf, 
                                      const double mfp) const {
    double Mz = SM.getMz(), cW2 = Mw*Mw/(Mz*Mz), sW2 = 1.0 - cW2;
    double ve = - 0.5 + 2.0*sW2, ae = -0.5, deltae = ve - ae;
    double vf = I3f - 2.0*Qf*sW2, af = I3f, deltaf = vf - af;
    double mu = Mw; // renormalization scale
    
    complex FW_hat;
    if (mfp==SM.getMtpole())
        FW_hat = F_W_t_hat(mu, s, Mw);
    else 
        FW_hat = F_W_0_hat(mu, s, Mw);      

    complex kappaef = 1.0;
    kappaef += SM.getAle()/4.0/M_PI/sW2    
              *( - 2.0*cW2/sW2*DeltaRhobar(mu, Mw) + 2.0*Pibar_Zgamma_hat(mu, s, Mw)
                 - PV.B0(mu, s, Mw, Mw)/3.0 - 2.0/9.0 
                 - ((deltae*deltae + deltaf*deltaf)/sW2*(Mw*Mw/s - 1.0)
                    + 3.0*ve*ve + ae*ae + 3.0*vf*vf + af*af)*F_za_0(s, Mw)/4.0/cW2
                 - F_W_0_hat(mu, s, Mw) - FW_hat
                 + (Mz*Mz/s - 1.0)*cW2*(2.0/3.0 + Pibar_gg_bos_hat(mu, s, Mw))
                 + mfp*mfp/4.0/Mw/Mw*( PV.B0(mu, s, Mw, Mw) + 1.0) );
    return kappaef;
}


complex EWSMTwoFermionsLEP2::G_e(const double s, const double Mw, const double I3f, 
                                 const double Qf, const double mfp) const {
    double Mz = SM.getMz(), sW2 = 1.0 - Mw*Mw/(Mz*Mz);
    return ( 1.0 - 4.0*( kappa_e(s,Mw,I3f,Qf,mfp)*sW2 + I2e(s,Mw)) ); 
}


complex EWSMTwoFermionsLEP2::G_f(const double s, const double Mw, const double I3f, 
                                 const double Qf, const double mfp) const {
    double Mz = SM.getMz(), sW2 = 1.0 - Mw*Mw/(Mz*Mz);
    return ( 1.0 - 4.0*fabs(Qf)*( kappa_f(s,Mw,I3f,Qf,mfp)*sW2 + I2f(s,Mw)) );     
}


complex EWSMTwoFermionsLEP2::G_ef(const double s, const double Mw, const double I3f, 
                                  const double Qf, const double mfp) const {
    double Mz = SM.getMz(), sW2 = 1.0 - Mw*Mw/(Mz*Mz);
    return ( - 1.0 + G_e(s,Mw,I3f,Qf,mfp) + G_f(s,Mw,I3f,Qf,mfp) 
             + 16.0*fabs(Qf)
               *( kappa_ef(s,Mw,I3f,Qf,mfp)*sW2*sW2 
                 + (kappa_e(s,Mw,I3f,Qf,mfp)*I2e(s,Mw) 
                    + kappa_f(s,Mw,I3f,Qf,mfp)*I2f(s,Mw))*sW2 ) );
}

   
double EWSMTwoFermionsLEP2::G_1(const double s, const double Mw, 
                                const double GammaZ, const double I3f, 
                                const double Qf, const double mfp,
                                const bool bWeak) const {
    complex Vpol = 1.0/Vpol_inv(s);
    complex rhoef, Ge, Gf, Gef;
    if (bWeak) {
        rhoef = rho_ef(s, Mw, I3f, Qf, mfp);
        Ge = G_e(s, Mw, I3f, Qf, mfp);
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
                                const double Qf, const double mfp,
                                const bool bWeak) const {
    complex Vpol = 1.0/Vpol_inv(s);
    complex rhoef, Ge, Gf, Gef;
    if (bWeak) {
        rhoef = rho_ef(s, Mw, I3f, Qf, mfp);
        Ge = G_e(s, Mw, I3f, Qf, mfp);
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
                                const double Qf, const double mfp,
                                const bool bWeak) const {
    complex Vpol = 1.0/Vpol_inv(s);
    complex rhoef, Ge, Gf, Gef;
    if (bWeak) {
        rhoef = rho_ef(s, Mw, I3f, Qf, mfp);
        Ge = G_e(s, Mw, I3f, Qf, mfp);
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


////////////////////////////////////////////////////////////////////////

complex EWSMTwoFermionsLEP2::I2e(const double s, const double Mw) const {
    double Mz = SM.getMz(), sW2 = 1.0 - Mw*Mw/(Mz*Mz);
    double alpha = SM.getAle()/Vpol_inv(s).real();
    double ReKappa_e = 1.0;
    return ( 35.0*alpha*alpha/18.0*( 1.0 - 8.0/3.0*ReKappa_e*sW2 ) );
}


complex EWSMTwoFermionsLEP2::I2f(const double s, const double Mw) const {
    double Mz = SM.getMz(), sW2 = 1.0 - Mw*Mw/(Mz*Mz);
    double alpha = SM.getAle()/Vpol_inv(s).real();
    double ReKappa_f = 1.0;
    return ( 35.0*alpha*alpha/18.0*( 1.0 - 8.0/3.0*ReKappa_f*sW2 ) );
}


complex EWSMTwoFermionsLEP2::DeltaRhobar(const double mu, const double Mw) const {
    return myOneLoopEW.DeltaRhobar(mu,Mw);
}


complex EWSMTwoFermionsLEP2::DeltaRhobarZ(const double mu, const double Mw) const {
    return ( myOneLoopEW.DeltaRhobar(mu,Mw) + myOneLoopEW.DeltaRhobarW(mu,Mw) );
}


complex EWSMTwoFermionsLEP2::D_Z_hat(const double mu, const double s, 
                                     const double Mw) const {
    double Mz = SM.getMz(), cW2 = Mw*Mw/(Mz*Mz);
    double Rz = Mz*Mz/s;
    complex D_Z_bos = (myOneLoopEW.SigmaZZ_bos(mu,s,Mw) 
                       - myOneLoopEW.SigmaZZ_bos(mu,Mz*Mz,Mw))/cW2/(Mz*Mz - s);
    complex D_Z_fer = (myOneLoopEW.SigmaZZ_fer(mu,s,Mw) 
                       - myOneLoopEW.SigmaZZ_fer(mu,Mz*Mz,Mw))/cW2/(Mz*Mz - s);
    complex D_Z_bos_hat = D_Z_bos 
                          + ( (1.0/12.0/cW2 + 4.0/3.0)/Rz + 1.0/12.0/cW2/Rz/Rz )
                             *PV.B0(mu,s,Mw,Mw)
                          + ( (1.0/cW2 - 13.0)/Rz + 1.0/cW2/Rz/Rz )/18.0;
    return ( D_Z_bos_hat + D_Z_fer );
}


complex EWSMTwoFermionsLEP2::Pibar_Zgamma_hat(const double mu, const double s, 
                                              const double Mw) const {
    double Mz = SM.getMz(), cW2 = Mw*Mw/(Mz*Mz);
    double Rw = Mw*Mw/s;
    complex Pibar_Zgamma_bos = myOneLoopEW.PiZgamma_bos(mu,s,Mw);
    complex Pibar_Zgamma_fer = myOneLoopEW.PiZgamma_fer(mu,s,Mw);
    complex Pibar_Zgamma_bos_hat = Pibar_Zgamma_bos 
                                   - cW2*(4.0*Rw + 17.0/3.0)*PV.B0(mu,s,Mw,Mw);
    return ( Pibar_Zgamma_bos_hat + Pibar_Zgamma_fer );
}


complex EWSMTwoFermionsLEP2::Pibar_gg_bos_hat(const double mu, const double s, 
                                              const double Mw) const {
    double Rw = Mw*Mw/s;
    complex Pibar_gg_bos = myOneLoopEW.PiGammaGamma_bos(mu,s,Mw);
    complex Pibar_gg_bos_hat = Pibar_gg_bos 
                               - (4.0*Rw + 17.0/3.0)*PV.B0(mu,s,Mw,Mw);
    return Pibar_gg_bos_hat;    
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


complex EWSMTwoFermionsLEP2::F_Wn_0_hat(const double mu, const double s, 
                                        const double Mw) const {
    double Rw = Mw*Mw/s;
    complex Fwn0 = myOneLoopEW.FWn_0(s, Mw);  
    complex Fwn0_hat = Fwn0 - (3.0/2.0/Rw + 1.0/12.0/Rw/Rw)*PV.B0(mu,s,Mw,Mw) 
                       + 11.0/18.0/Rw - 1.0/18.0/Rw/Rw;
    return Fwn0_hat;
}   


complex EWSMTwoFermionsLEP2::F_Wn_t_hat(const double mu, const double s, 
                                        const double Mw) const {
    double Rw = Mw*Mw/s, wt = SM.getMtpole()*SM.getMtpole()/Mw/Mw;
    complex Fwnt = myOneLoopEW.FWn_t(s, Mw);
    complex Fwnt_hat = Fwnt + wt/4.0/Rw*(PV.B0(mu,s,Mw,Mw) + 1.0);
    return Fwnt_hat;
}   


complex EWSMTwoFermionsLEP2::F_W_0_hat(const double mu, const double s, 
                                       const double Mw) const {
    double Mz = SM.getMz(), cW2 = Mw*Mw/(Mz*Mz);
    return ( cW2*F_Wn_0_hat(mu, s, Mw) - F_Wa_0(s, Mw)/2.0 );
}


complex EWSMTwoFermionsLEP2::F_W_t_hat(const double mu, const double s, 
                                       const double Mw) const {
    double Mz = SM.getMz(), cW2 = Mw*Mw/(Mz*Mz);
    return ( cW2*F_Wn_t_hat(mu, s, Mw) - F_Wa_t(s, Mw)/2.0 
             - myOneLoopEW.FbarWa_t(s, Mw)/2.0 );
}



