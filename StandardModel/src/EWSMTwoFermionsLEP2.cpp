/* 
 * File:   EWSMTwoFermionsLEP2.cpp
 * Author: mishima
 */

#include <cmath>
#include <stdexcept>
#include "EWSMTwoFermionsLEP2.h"


EWSMTwoFermionsLEP2::EWSMTwoFermionsLEP2(const StandardModel& SM_i) : SM(SM_i), 
        myOneLoopEW_HV(SM_i) {
}


////////////////////////////////////////////////////////////////////////

double EWSMTwoFermionsLEP2::sigma_l(const StandardModel::lepton l, const double s,
                                    const double Mw, const double GammaZ, 
                                    const bool bDP, const bool bQED) const {
    double mf = SM.getLeptons(l).getMass();
    double Qf = SM.getLeptons(l).getCharge();
    double I3f = SM.getLeptons(l).getIsospin();
    
    //-----------------------------------
    // Test for renormalization scale dependence of the renormalized self-energies
    std::cout << "TEST1 (mu=Mz)      : " 
              << Sigma_hat_ZZ(SM.getMz(),s,Mw) << " "
              << Sigma_hat_gZ(SM.getMz(),s,Mw) << " "           
              << Sigma_hat_gg(SM.getMz(),s,Mw) << std::endl;
    std::cout << "TEST2 (mu=sqrt(s)) : " 
              << Sigma_hat_ZZ(sqrt(s),s,Mw) << " "
              << Sigma_hat_gZ(sqrt(s),s,Mw) << " "           
              << Sigma_hat_gg(sqrt(s),s,Mw) << std::endl;
    //-----------------------------------
    
    return sigma_f(s, Mw, GammaZ, mf, Qf, I3f, 1.0, bDP, bQED);
}


double EWSMTwoFermionsLEP2::sigma_q(const StandardModel::quark q, const double s,
                                    const double Mw, const double GammaZ, 
                                    const bool bDP, const bool bQED) const {
    double mf = myOneLoopEW_HV.mq(q, s); // m_q(s)
    double Qf = SM.getQuarks(q).getCharge();
    double I3f = SM.getQuarks(q).getIsospin();
    
    return sigma_f(s, Mw, GammaZ, mf, Qf, I3f, 3.0, bDP, bQED);
}


double EWSMTwoFermionsLEP2::sigma_f(const double s, const double Mw, const double GammaZ, 
                                    const double mf, const double Qf,  const double I3f, 
                                    const double Ncf,
                                    const bool bDP, const bool bQED) const {
    double betaf = sqrt(1.0 - 4.0*mf*mf/s);
    double Qe = SM.getLeptons(SM.ELECTRON).getCharge();
    double I3e = SM.getLeptons(SM.ELECTRON).getIsospin();
    double Mz = SM.getMz();
    double cW2 = Mw*Mw/Mz/Mz, cW = sqrt(cW2);
    double sW2 = 1.0 - cW2, sW = sqrt(sW2);
    double ve = (I3e - 2.0*Qe*sW2)/(2.0*sW*cW);
    double ae = I3e/(2.0*sW*cW);
    double vf = (I3f - 2.0*Qf*sW2)/(2.0*sW*cW);
    double af = I3f/(2.0*sW*cW);
    double Qe2 = Qe*Qe, Qf2 = Qf*Qf, betaf2 = betaf*betaf;
    double ve2 = ve*ve, ae2 = ae*ae, vf2 = vf*vf, af2 = af*af;

    // Tree-level or dressed gauge-boson propagators
    complex chiG, chiZ, chiGZ; 
    double mu = sqrt(s); // renormalization scale (Check scale dependence!!)
    //mu = Mz; // test for scale dependence

    if (bDP) {
        chiG = chi_gamma(mu,s,Mw);
        chiZ = chi_Z(mu,s,Mw);
        chiGZ = chi_gammaZ(mu,s,Mw);
    } else {
        chiG = complex(1.0, 0.0, false);
        complex denom = complex(s - Mz*Mz, Mz*GammaZ, false);
        chiZ = s/denom;
        chiGZ = complex(0.0, 0.0, false);
    }    
    
    double G1 = Qe2*Qf2*chiG.abs2()
                + 2.0*ve*vf*Qe*Qf*(chiZ*chiG.conjugate()).real()
                + (ve2 + ae2)*(vf2 + betaf2*af2)*chiZ.abs2()
                + (Qe2*(vf2 + af2) + 2.0*ve*vf*Qe*Qf + (ve2 + ae2)*Qf2)*chiGZ.abs2()
                + 2.0*(vf*Qe2*Qf + ve*Qe*Qf2)*(chiG*chiGZ.conjugate()).real()
                + 2.0*(ve*(vf2 + af2)*Qe + (ve2 + ae2)*vf*Qf)*(chiZ*chiGZ.conjugate()).real();
    double G2 = Qe2*Qf2*chiG.abs2()
                + 2.0*ve*vf*Qe*Qf*(chiZ*chiG.conjugate()).real()
                + (ve2 + ae2)*vf2*chiZ.abs2();
    
    // QED corrections
    if (bQED) {
        G1 += Qf*Qf*C11V(s,mf,Qf)*chiG.abs2()
              + 2.0*Qe*Qf*( (ve*vf*C12V(s,GammaZ,mf,Qf) + ae*af*C12A(s,mf,Qf))
                            *chiG*chiGZ.conjugate() ).real()
              + ( (ve2 + ae2)*(vf2 + af2)*C22V(s,GammaZ,mf,Qf) 
                  + 4.0*ve*ae*vf*af*C22A(s,mf,Qf) )*chiZ.abs2();
    }
    
    return ( 4.0*M_PI*SM.getAle()*SM.getAle()/(3.0*s)*Ncf*betaf*(G1 + 2.0*mf*mf/s*G2) );
}


//////////////////////////////////////////////////////////////////////// 
// Renormalized self-energies

complex EWSMTwoFermionsLEP2::Sigma_hat_ZZ(const double mu, const double s, 
                                          const double Mw) const {
    double Mw2 = Mw*Mw;
    double Mz = SM.getMz(), Mz2 = Mz*Mz;
    double cW2 = Mw2/Mz2, cW = sqrt(cW2);
    double sW2 = 1.0 - cW2, sW = sqrt(sW2);
    
    // Bosonic contributions to self-energies
    complex Sigma_WW_Mw2 = myOneLoopEW_HV.SigmaWW_bos(mu, Mw2, Mw);
    complex Sigma_ZZ_s   = myOneLoopEW_HV.SigmaZZ_bos(mu, s, Mw);
    complex Sigma_ZZ_Mz2 = myOneLoopEW_HV.SigmaZZ_bos(mu, Mz2, Mw);
    complex Sigma_Zg_0   = myOneLoopEW_HV.SigmaZgamma_bos(mu, 0.0, Mw);
    complex Pi_gg_0      = myOneLoopEW_HV.PiGammaGamma_bos(mu, 0.0, Mw);
    
    //-- TEST (use the self-energies in Hollik's paper) --
    Sigma_WW_Mw2 = myOneLoopEW_HV.SigmaWW_bos_Hollik(mu, Mw2, Mw); // for test
    Sigma_ZZ_s   = myOneLoopEW_HV.SigmaZZ_bos_Hollik(mu, s, Mw); // for test
    Sigma_ZZ_Mz2 = myOneLoopEW_HV.SigmaZZ_bos_Hollik(mu, Mz2, Mw); // for test
    Sigma_Zg_0   = myOneLoopEW_HV.SigmaZgamma_bos_Hollik(mu, 0.0, Mw); // for test
    Pi_gg_0      = myOneLoopEW_HV.PiGammaGamma_bos_Hollik(mu, 0.0, Mw); // for test
    //--------------------
    
    // Fermionic contributions to self-energies
    double muForMq = s; // renormalization scale for the running quark mass
    Sigma_WW_Mw2 += myOneLoopEW_HV.SigmaWW_fer(mu, muForMq, Mw2);
    Sigma_ZZ_s   += myOneLoopEW_HV.SigmaZZ_fer(mu, muForMq, s, Mw);
    Sigma_ZZ_Mz2 += myOneLoopEW_HV.SigmaZZ_fer(mu, muForMq, Mz2, Mw);
    Sigma_Zg_0   += myOneLoopEW_HV.SigmaZgamma_fer(mu, muForMq, 0.0, Mw); 
    Pi_gg_0      += myOneLoopEW_HV.PiGammaGamma_fer(mu, muForMq, 0.0);
    
    // Refactoring
    Sigma_WW_Mw2 *= SM.getAle()/4.0/M_PI/sW2;
    Sigma_ZZ_s   *= SM.getAle()/4.0/M_PI/sW2/cW2;
    Sigma_ZZ_Mz2 *= SM.getAle()/4.0/M_PI/sW2/cW2;
    Sigma_Zg_0   *= - SM.getAle()/4.0/M_PI/sW/cW;
    Pi_gg_0      *= SM.getAle()/4.0/M_PI;
    
    // Counter terms for the mass renormalization
    double deltaMw2 = Sigma_WW_Mw2.real();
    double deltaMz2 = Sigma_ZZ_Mz2.real();
    
    // Counter terms for the wave-function renormalization
    complex deltaZz = - Pi_gg_0 + (cW2-sW2)/sW2*(deltaMz2/Mz2 - deltaMw2/Mw2 
                                                 + 2.0*sW/cW*Sigma_Zg_0/Mz2);
    
    return ( Sigma_ZZ_s + (s-Mz2)*deltaZz - deltaMz2 );
}


complex EWSMTwoFermionsLEP2::Sigma_hat_gZ(const double mu, const double s,
                                          const double Mw) const {
    double Mw2 = Mw*Mw;
    double Mz = SM.getMz(), Mz2 = Mz*Mz;
    double cW2 = Mw2/Mz2, cW = sqrt(cW2);
    double sW2 = 1.0 - cW2, sW = sqrt(sW2);

    // Bosonic contributions to self-energies
    complex Sigma_WW_Mw2 = myOneLoopEW_HV.SigmaWW_bos(mu, Mw2, Mw);
    complex Sigma_ZZ_Mz2 = myOneLoopEW_HV.SigmaZZ_bos(mu, Mz2, Mw);
    complex Sigma_Zg_s   = myOneLoopEW_HV.SigmaZgamma_bos(mu, s, Mw);
    complex Sigma_Zg_0   = myOneLoopEW_HV.SigmaZgamma_bos(mu, 0.0, Mw);
    
    //-- TEST (use the self-energies in Hollik's paper) --
    Sigma_WW_Mw2 = myOneLoopEW_HV.SigmaWW_bos_Hollik(mu, Mw2, Mw); // for test
    Sigma_ZZ_Mz2 = myOneLoopEW_HV.SigmaZZ_bos_Hollik(mu, Mz2, Mw); // for test
    Sigma_Zg_s   = myOneLoopEW_HV.SigmaZgamma_bos_Hollik(mu, s, Mw); // for test
    Sigma_Zg_0   = myOneLoopEW_HV.SigmaZgamma_bos_Hollik(mu, 0.0, Mw); // for test
    //--------------------
    
    // Fermionic contributions to self-energies
    double muForMq = s; // renormalization scale for the running quark mass
    Sigma_WW_Mw2 += myOneLoopEW_HV.SigmaWW_fer(mu, muForMq, Mw2);
    Sigma_ZZ_Mz2 += myOneLoopEW_HV.SigmaZZ_fer(mu, muForMq, Mz2, Mw);
    Sigma_Zg_s   += myOneLoopEW_HV.SigmaZgamma_fer(mu, muForMq, s, Mw); 
    Sigma_Zg_0   += myOneLoopEW_HV.SigmaZgamma_fer(mu, muForMq, 0.0, Mw); 

    // Refactoring
    Sigma_WW_Mw2 *= SM.getAle()/4.0/M_PI/sW2;
    Sigma_ZZ_Mz2 *= SM.getAle()/4.0/M_PI/sW2/cW2;
    Sigma_Zg_s   *= - SM.getAle()/4.0/M_PI/sW/cW;
    Sigma_Zg_0   *= - SM.getAle()/4.0/M_PI/sW/cW;
    
    // Counter terms for the mass renormalization
    double deltaMw2 = Sigma_WW_Mw2.real();
    double deltaMz2 = Sigma_ZZ_Mz2.real();
    
    return ( Sigma_Zg_s - Sigma_Zg_0 
             + s*( cW/sW*(deltaMz2/Mz2 - deltaMw2/Mw2) + 2.0*Sigma_Zg_0/Mz2 ) );
}


complex EWSMTwoFermionsLEP2::Sigma_hat_gg(const double mu, const double s,
                                          const double Mw) const {
    // Bosonic contributions to self-energies
    complex Sigma_gg_s = myOneLoopEW_HV.SigmaGammaGamma_bos(mu, s, Mw);
    complex Pi_gg_0    = myOneLoopEW_HV.PiGammaGamma_bos(mu, 0.0, Mw);
    
    //-- TEST (use the self-energies in Hollik's paper) --
    Sigma_gg_s  = myOneLoopEW_HV.SigmaGammaGamma_bos_Hollik(mu, s, Mw); // for test
    Pi_gg_0  = myOneLoopEW_HV.PiGammaGamma_bos_Hollik(mu, 0.0, Mw); // for test
    //--------------------
    
    // Fermionic contributions to self-energies
    double muForMq = s; // renormalization scale for the running quark mass
    Sigma_gg_s += myOneLoopEW_HV.SigmaGammaGamma_fer(mu, muForMq, s);
    Pi_gg_0    += myOneLoopEW_HV.PiGammaGamma_fer(mu, muForMq, 0.0);
    
    // Refactoring
    Sigma_gg_s *= SM.getAle()/4.0/M_PI;
    Pi_gg_0    *= SM.getAle()/4.0/M_PI;
    
    return ( Sigma_gg_s - s*Pi_gg_0 );
}


//////////////////////////////////////////////////////////////////////// 
// Dressed gauge-boson propagators

complex EWSMTwoFermionsLEP2::chi_Z(const double mu, const double s, 
                                   const double Mw) const {
    double Mz = SM.getMz();

    complex D_Z = 1.0/(s - Mz*Mz + Sigma_hat_ZZ(mu,s,Mw));
    return (s*D_Z);
}


complex EWSMTwoFermionsLEP2::chi_gamma(const double mu, const double s, 
                                       const double Mw) const {
    complex D_gamma = 1.0/(s + Sigma_hat_gg(mu,s,Mw));
    return (s*D_gamma);
}


complex EWSMTwoFermionsLEP2::chi_gammaZ(const double mu, const double s, 
                                        const double Mw) const {
    // O(alpha) approximation
    return ( Sigma_hat_gZ(mu,s,Mw)/s*chi_Z(mu,s,Mw) );
}


//////////////////////////////////////////////////////////////////////// 
// QED corrections    

double EWSMTwoFermionsLEP2::delta() const {
    return ( 1.0 - 0.85*0.85 ); // sqrt{s'} > 0.85*sqrt{s}
    
}


double EWSMTwoFermionsLEP2::Bf(const double s, const double mf) const {
    return ( log(s/mf/mf) - 1.0 );
    
}


double EWSMTwoFermionsLEP2::gamma_delta(const double s, const double mf, const double Qf) const {
    double me = SM.getLeptons(SM.ELECTRON).getMass();

    return ( 2.0*SM.getAle()/M_PI
             *( Bf(s,me) + Qf*Qf*Bf(s,mf) )*log(delta()) );
}


complex EWSMTwoFermionsLEP2::gamma_delta_int(const double s, const double GammaZ, 
                                             const double mf, const double Qf) const {
    double me = SM.getLeptons(SM.ELECTRON).getMass();
    double Mz = SM.getMz();
    complex M2 = complex(Mz*Mz, -Mz*GammaZ, false);
    double d = delta();

    return ( 2.0*SM.getAle()/M_PI
             *( Bf(s,me)*log(d*(s-M2)/(s-s*d-M2)) + Qf*Qf*Bf(s,mf)*log(d) ) );
}


double EWSMTwoFermionsLEP2::gamma_delta_res(const double s, const double GammaZ, 
                                            const double mf, const double Qf) const {
    double me = SM.getLeptons(SM.ELECTRON).getMass();
    double Mz = SM.getMz();
    complex M2 = complex(Mz*Mz, -Mz*GammaZ, false);
    double d = delta();

    return ( 2.0*SM.getAle()/M_PI
             *( Bf(s,me)*log( (d*(s-M2)/(s-s*d-M2)).abs() ) 
                + Qf*Qf*Bf(s,mf)*log(d) ) );
}


double EWSMTwoFermionsLEP2::gamma_tail(const double s, const double GammaZ) const {
    double me = SM.getLeptons(SM.ELECTRON).getMass();
    double Mz = SM.getMz();
    double d = delta();
    
    return ( 2.0*SM.getAle()/M_PI
             *Bf(s,me)*(s-Mz*Mz)/Mz/GammaZ
             *( atan((Mz*Mz-s+s*d)/Mz/GammaZ) - atan((Mz*Mz-s)/Mz/GammaZ) ) );
}


double EWSMTwoFermionsLEP2::gamma_fin(const double s, const double mf, const double Qf) const {
    double me = SM.getLeptons(SM.ELECTRON).getMass();

    return ( 3.0*SM.getAle()/2.0/M_PI * (Bf(s,me) + Qf*Qf*Bf(s,mf))
             + SM.getAle()/M_PI * (1 + Qf*Qf)*(M_PI*M_PI/3.0 - 1.0/2.0) );
}


double EWSMTwoFermionsLEP2::C11V(const double s, const double mf, const double Qf) const {
    return ( gamma_delta(s,mf,Qf) + gamma_fin(s,mf,Qf) );
}  


double EWSMTwoFermionsLEP2::C11A(const double s, const double mf, const double Qf) const {
    return 0.0;
}


complex EWSMTwoFermionsLEP2::C12V(const double s, const double GammaZ, 
                                  const double mf, const double Qf) const {
    return ( gamma_delta_int(s,GammaZ,mf,Qf).conjugate() + gamma_fin(s,mf,Qf) );
}


complex EWSMTwoFermionsLEP2::C12A(const double s, const double mf, const double Qf) const {
    return complex(0.0, 0.0, false);    
}


double EWSMTwoFermionsLEP2::C22V(const double s, const double GammaZ, 
                                 const double mf, const double Qf) const {
    return ( gamma_delta_res(s,GammaZ,mf,Qf) + gamma_tail(s,GammaZ) + gamma_fin(s,mf,Qf) );
}


double EWSMTwoFermionsLEP2::C22A(const double s, const double mf, const double Qf) const {
    return 0.0;   
}





