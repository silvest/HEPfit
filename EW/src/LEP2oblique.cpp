/* 
 * File:   LEP2oblique.cpp
 * Author: mishima
 */

#include "LEP2oblique.h"


LEP2oblique::LEP2oblique(const EW& EW_i) : myEW(EW_i) {
}


////////////////////////////////////////////////////////////////////////

double LEP2oblique::sigma_l_LEP2_NP(const StandardModel::lepton l, 
                                    const double s) const {
    double Ncf = 1.0;
    double betaf = sqrt(1.0 - 4.0*ml(l)*ml(l)/s);
    double alpha0 = myEW.getSM().alphaMz();
    
    return ( 4.0*M_PI*alpha0*alpha0/(3.0*s)*Ncf*betaf*G1_l_NP(l, s) );
}


double LEP2oblique::sigma_q_LEP2_NP(const StandardModel::quark q, 
                                    const double s) const {
    double Ncf = 3.0;
    double m_q = mq(q,sqrt(s));
    double betaf = sqrt(1.0 - 4.0*m_q*m_q/s);
    double alpha0 = myEW.getSM().alphaMz();
    
    return ( 4.0*M_PI*alpha0*alpha0/(3.0*s)*Ncf*betaf*G1_q_NP(q, s) );
}


double LEP2oblique::AFB_l_LEP2_NP(const StandardModel::lepton l, 
                                  const double s) const {
    double mf2 = ml(l)*ml(l);
    double betaf = sqrt(1.0 - 4.0*mf2/s);
    double G1SM0 = G1_l_SM0(l, s), G2SM0 = G2_l_SM0(l, s), G3SM0 = G3_l_SM0(l, s);
    double AFB_Born0 = 3.0/4.0*betaf*G3SM0/(G1SM0 + 2.0*mf2/s*G2SM0);
    
    return ( - AFB_Born0*G1_l_NP(l, s)/(G1SM0 + 2.0*mf2/s*G2SM0)
             + AFB_Born0*G3_l_NP(l, s)/G3SM0 );
}


double LEP2oblique::AFB_q_LEP2_NP(const StandardModel::quark q, 
                                  const double s) const {
    double m_q = mq(q,sqrt(s)), mf2 = m_q*m_q;
    double betaf = sqrt(1.0 - 4.0*mf2/s);
    double G1SM0 = G1_q_SM0(q, s), G2SM0 = G2_q_SM0(q, s), G3SM0 = G3_q_SM0(q, s);
    double AFB_Born0 = 3.0/4.0*betaf*G3SM0/(G1SM0 + 2.0*mf2/s*G2SM0);
    
    return ( - AFB_Born0*G1_q_NP(q, s)/(G1SM0 + 2.0*mf2/s*G2SM0)
             + AFB_Born0*G3_q_NP(q, s)/G3SM0 );  
}


double LEP2oblique::R_q_LEP2_NP(const StandardModel::quark q, 
                                const double s) const {
    double sigma_q_SM0 = sigma_q_LEP2_SM0(q, s);
    double sigma_had_SM0 = sigma_q_LEP2_SM0(StandardModel::UP, s)
                         + sigma_q_LEP2_SM0(StandardModel::DOWN, s)
                         + sigma_q_LEP2_SM0(StandardModel::CHARM, s)
                         + sigma_q_LEP2_SM0(StandardModel::STRANGE, s)
                         + sigma_q_LEP2_SM0(StandardModel::BOTTOM, s);
    double sigma_q_NP = sigma_q_LEP2_NP(q, s);
    double sigma_had_NP = sigma_q_LEP2_NP(StandardModel::UP, s)
                        + sigma_q_LEP2_NP(StandardModel::DOWN, s)
                        + sigma_q_LEP2_NP(StandardModel::CHARM, s)
                        + sigma_q_LEP2_NP(StandardModel::STRANGE, s)
                        + sigma_q_LEP2_NP(StandardModel::BOTTOM, s);
    
    return ( - sigma_q_SM0/(sigma_had_SM0*sigma_had_SM0)*sigma_had_NP
             + sigma_q_NP/sigma_had_SM0 );
}


////////////////////////////////////////////////////////////////////////
   
double LEP2oblique::DeltaEpsilon_1() const {
    double c0 = sqrt(myEW.c02()), s0 = sqrt(myEW.s02());
    return ( myEW.That() - myEW.W() + 2.0*s0/c0*myEW.X() - s0*s0/c0/c0*myEW.Y() );
}


double LEP2oblique::DeltaEpsilon_2() const {
    double c0 = sqrt(myEW.c02()), s0 = sqrt(myEW.s02());
    return ( myEW.Uhat() - myEW.V() - myEW.W() + 2.0*s0/c0*myEW.X() );
}


double LEP2oblique::DeltaEpsilon_3() const {
    double c0 = sqrt(myEW.c02()), s0 = sqrt(myEW.s02());
    return ( myEW.Shat() - myEW.W() + myEW.X()/s0/c0 - myEW.Y() );
}


double LEP2oblique::epsilonZZ() const {
    double c0 = sqrt(myEW.c02()), s0 = sqrt(myEW.s02());
    return ( myEW.c02()*myEW.W() - 2.0*s0*c0*myEW.X() + myEW.s02()*myEW.Y() );
}


double LEP2oblique::epsilonGammaGamma() const {
    double c0 = sqrt(myEW.c02()), s0 = sqrt(myEW.s02());
    return ( myEW.s02()*myEW.W() + 2.0*s0*c0*myEW.X() + myEW.c02()*myEW.Y() );
}


double LEP2oblique::epsilonGammaZ() const {
    double c0 = sqrt(myEW.c02()), s0 = sqrt(myEW.s02());
    return ( (myEW.c02() - myEW.s02())*myEW.X() + s0*c0*(myEW.W() - myEW.Y()) );    
}


double LEP2oblique::vl(const StandardModel::lepton l) const {
    double c0 = sqrt(myEW.c02()), s0 = sqrt(myEW.s02());
    return ( - (myEW.getSM().getLeptons(l).getIsospin() 
                - 2.0*myEW.Ql(l)*myEW.s02())/(2.0*s0*c0) );
}


double LEP2oblique::vq(const StandardModel::quark q) const {
    double c0 = sqrt(myEW.c02()), s0 = sqrt(myEW.s02());
    return ( - (myEW.getSM().getQuarks(q).getIsospin() 
                - 2.0*myEW.Qq(q)*myEW.s02())/(2.0*s0*c0) );
}


double LEP2oblique::al(const StandardModel::lepton l) const {
    double c0 = sqrt(myEW.c02()), s0 = sqrt(myEW.s02());
    return ( - myEW.getSM().getLeptons(l).getIsospin()/(2.0*s0*c0) ); 
}


double LEP2oblique::aq(const StandardModel::quark q) const {
    double c0 = sqrt(myEW.c02()), s0 = sqrt(myEW.s02());
    return ( - myEW.getSM().getQuarks(q).getIsospin()/(2.0*s0*c0) );     
}


double LEP2oblique::G1_NP(const double s, const double Qf, 
                          const double vf, const double af) const {
    double epsilonbarGamma = - s/(myEW.Mw0()*myEW.Mw0())*epsilonGammaGamma();
    double epsilonbarZ = - s/(myEW.Mw0()*myEW.Mw0())*epsilonZZ();
    double epsilonbarGammaZ = s/(myEW.Mw0()*myEW.Mw0())*epsilonGammaZ();
    double Qe = myEW.Ql(StandardModel::ELECTRON);
    double ve = vl(StandardModel::ELECTRON);
    double ae = al(StandardModel::ELECTRON);
    double Qe2 = Qe*Qe, Qf2 = Qf*Qf;
    double ve2 = ve*ve, vf2 = vf*vf, ae2 = ae*ae, af2 = af*af;
    double Mz = myEW.getSM().getMz();
    double GammaZ0 = 7.0*myEW.getSM().alphaMz()*Mz/(16.0*myEW.s02()*myEW.c02());
    complex denom = complex(s - Mz*Mz, Mz*GammaZ0, false);
    double Zprop = (1.0/denom).real();
    
    return ( 2.0*Qe2*Qf2*epsilonbarGamma
             + 2.0*ve*vf*Qe*Qf*(epsilonbarZ + s*epsilonbarGamma*Zprop)
             + 2.0*(ve2 + ae2)*(vf2 + af2)*s*epsilonbarZ*Zprop
             + 2.0*(vf*Qe2*Qf + ve*Qe*Qf2)*epsilonbarGammaZ
             + 2.0*(ve*(vf2 + af2)*Qe + (ve2 + ae2)*vf*Qf)*s*epsilonbarGammaZ*Zprop );    
}


double LEP2oblique::G1_l_NP(const StandardModel::lepton l, const double s) const {
    double Qf = myEW.Ql(l), vf = vl(l), af = al(l);    
    return ( G1_NP(s, Qf, vf, af) );
}


double LEP2oblique::G1_q_NP(const StandardModel::quark q, const double s) const {
    double Qf = myEW.Qq(q), vf = vq(q), af = aq(q);    
    return ( G1_NP(s, Qf, vf, af) );    
}


double LEP2oblique::G3_NP(const double s, const double Qf, 
                          const double vf, const double af) const {
    double epsilonbarGamma = - s/(myEW.Mw0()*myEW.Mw0())*epsilonGammaGamma();
    double epsilonbarZ = - s/(myEW.Mw0()*myEW.Mw0())*epsilonZZ();
    double epsilonbarGammaZ = s/(myEW.Mw0()*myEW.Mw0())*epsilonGammaZ();
    double Qe = myEW.Ql(StandardModel::ELECTRON);
    double ve = vl(StandardModel::ELECTRON);
    double ae = al(StandardModel::ELECTRON);
    double Mz = myEW.getSM().getMz();
    double GammaZ0 = 7.0*myEW.getSM().alphaMz()*Mz/(16.0*myEW.s02()*myEW.c02());
    complex denom = complex(s - Mz*Mz, Mz*GammaZ0, false);
    double Zprop = (1.0/denom).real();
    
    return ( 2.0*ae*af*Qe*Qf*(epsilonbarZ + s*epsilonbarGamma*Zprop)
             + 8.0*ve*ae*vf*af*s*epsilonbarZ*Zprop
             + 4.0*(ae*vf*af*Qe + ve*ae*af*Qf)*s*epsilonbarGammaZ*Zprop );   
}


double LEP2oblique::G3_l_NP(const StandardModel::lepton l, const double s) const {
    double Qf = myEW.Ql(l), vf = vl(l), af = al(l);    
    return ( G3_NP(s, Qf, vf, af) );
}


double LEP2oblique::G3_q_NP(const StandardModel::quark q, const double s) const {
    double Qf = myEW.Qq(q), vf = vq(q), af = aq(q);    
    return ( G3_NP(s, Qf, vf, af) ); 
}


double LEP2oblique::G1_SM0(const double s, const double Qf, 
                           const double vf, const double af) const {
    double Qe = myEW.Ql(StandardModel::ELECTRON);
    double ve = vl(StandardModel::ELECTRON);
    double ae = al(StandardModel::ELECTRON);
    double Qe2 = Qe*Qe, Qf2 = Qf*Qf;
    double ve2 = ve*ve, vf2 = vf*vf, ae2 = ae*ae, af2 = af*af;
    double Mz = myEW.getSM().getMz();
    double GammaZ0 = 7.0*myEW.getSM().alphaMz()*Mz/(16.0*myEW.s02()*myEW.c02());
    complex denom = complex(s - Mz*Mz, Mz*GammaZ0, false);
    complex chiZ = s/denom;
    
    return ( Qe2*Qf2 + 2.0*ve*vf*Qe*Qf*chiZ.real()
             + (ve2 + ae2)*(vf2 + af2)*chiZ.abs2() );
}    


double LEP2oblique::G1_l_SM0(const StandardModel::lepton l, const double s) const {
    double Qf = myEW.Ql(l), vf = vl(l), af = al(l);    
    return ( G1_SM0(s, Qf, vf, af) );
}


double LEP2oblique::G1_q_SM0(const StandardModel::quark q, const double s) const {
    double Qf = myEW.Qq(q), vf = vq(q), af = aq(q);    
    return ( G1_SM0(s, Qf, vf, af) ); 
}


double LEP2oblique::G2_SM0(const double s, const double Qf, 
                           const double vf, const double af) const {
    double Qe = myEW.Ql(StandardModel::ELECTRON);
    double ve = vl(StandardModel::ELECTRON);
    double ae = al(StandardModel::ELECTRON);
    double Qe2 = Qe*Qe, Qf2 = Qf*Qf;
    double ve2 = ve*ve, vf2 = vf*vf, ae2 = ae*ae;
    double Mz = myEW.getSM().getMz();
    double GammaZ0 = 7.0*myEW.getSM().alphaMz()*Mz/(16.0*myEW.s02()*myEW.c02());
    complex denom = complex(s - Mz*Mz, Mz*GammaZ0, false);
    complex chiZ = s/denom;
    
    return ( Qe2*Qf2 + 2.0*ve*vf*Qe*Qf*chiZ.real()
             + (ve2 + ae2)*vf2*chiZ.abs2() );
}


double LEP2oblique::G2_l_SM0(const StandardModel::lepton l, const double s) const {
    double Qf = myEW.Ql(l), vf = vl(l), af = al(l);    
    return ( G2_SM0(s, Qf, vf, af) );    
}


double LEP2oblique::G2_q_SM0(const StandardModel::quark q, const double s) const {
    double Qf = myEW.Qq(q), vf = vq(q), af = aq(q);    
    return ( G2_SM0(s, Qf, vf, af) ); 
    
}


double LEP2oblique::G3_SM0(const double s, const double Qf, 
                           const double vf, const double af) const {
    double Qe = myEW.Ql(StandardModel::ELECTRON);
    double ve = vl(StandardModel::ELECTRON);
    double ae = al(StandardModel::ELECTRON);
    double Mz = myEW.getSM().getMz();
    double GammaZ0 = 7.0*myEW.getSM().alphaMz()*Mz/(16.0*myEW.s02()*myEW.c02());
    complex denom = complex(s - Mz*Mz, Mz*GammaZ0, false);
    complex chiZ = s/denom;

    return ( 2.0*ae*af*Qe*Qf*chiZ.real() + 4.0*ve*ae*vf*af*chiZ.abs2() );
}


double LEP2oblique::G3_l_SM0(const StandardModel::lepton l, const double s) const {
    double Qf = myEW.Ql(l), vf = vl(l), af = al(l);    
    return ( G3_SM0(s, Qf, vf, af) );
}


double LEP2oblique::G3_q_SM0(const StandardModel::quark q, const double s) const {
    double Qf = myEW.Qq(q), vf = vq(q), af = aq(q);    
    return ( G3_SM0(s, Qf, vf, af) );    
}


double LEP2oblique::sigma_l_LEP2_SM0(const StandardModel::lepton l, 
                                     const double s) const {
    double Ncf = 1.0;
    double m_l = ml(l), mf2 = m_l*m_l;
    double betaf = sqrt(1.0 - 4.0*mf2/s);
    double alpha0 = myEW.getSM().alphaMz();
    double G1SM0 = G1_l_SM0(l, s), G2SM0 = G2_l_SM0(l, s);

    return ( 4.0*M_PI*alpha0*alpha0/(3.0*s)*Ncf*betaf
             *(G1SM0 + 2.0*mf2/s*G2SM0) );    
}


double LEP2oblique::sigma_q_LEP2_SM0(const StandardModel::quark q, 
                                     const double s) const {
    double Ncf = 3.0;
    double m_q = mq(q,sqrt(s)), mf2 = m_q*m_q;
    double betaf = sqrt(1.0 - 4.0*mf2/s);
    double alpha0 = myEW.getSM().alphaMz();
    double G1SM0 = G1_q_SM0(q, s), G2SM0 = G2_q_SM0(q, s);

    return ( 4.0*M_PI*alpha0*alpha0/(3.0*s)*Ncf*betaf
             *(G1SM0 + 2.0*mf2/s*G2SM0) );
}


