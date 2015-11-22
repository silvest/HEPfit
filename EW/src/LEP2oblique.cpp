/* 
 * Copyright (C) 2012 HEPfit Collaboration
 * All rights reserved.
 *
 * For the licensing terms see doc/COPYING.
 */

#include "LEP2oblique.h"


LEP2oblique::LEP2oblique(const StandardModel& SM_i) 
: SM(SM_i) 
{
}


////////////////////////////////////////////////////////////////////////

double LEP2oblique::sigma_l_LEP2_NP(const StandardModel::lepton l, 
                                    const double s, const double ml,
                                    const double ObParam_i[]) const 
{
    double alpha0 = alpha_at_s(s);
    double Ncf = 1.0;
    double betaf = sqrt(1.0 - 4.0*ml*ml/s);
    double G1NP = G1_l_NP(l, s, alpha0, ObParam_i);
    
    return ( 4.0*M_PI*alpha0*alpha0/(3.0*s)*Ncf*betaf*G1NP );
}


double LEP2oblique::sigma_q_LEP2_NP(const QCD::quark q, 
                                    const double s, const double mq,
                                    const double ObParam_i[]) const
{
    double alpha0 = alpha_at_s(s);
    double Ncf = 3.0;
    double betaf = sqrt(1.0 - 4.0*mq*mq/s);
    double G1NP = G1_q_NP(q, s, alpha0, ObParam_i);
    
    return ( 4.0*M_PI*alpha0*alpha0/(3.0*s)*Ncf*betaf*G1NP );
}


double LEP2oblique::AFB_l_LEP2_NP(const StandardModel::lepton l, 
                                  const double s, const double ml,
                                  const double ObParam_i[]) const
{
    double alpha0 = alpha_at_s(s);
    double mf2 = ml*ml;
    double betaf = sqrt(1.0 - 4.0*mf2/s);
    double G1SM0 = G1_l_SM0(l, s, alpha0);
    double G2SM0 = G2_l_SM0(l, s, alpha0);
    double G3SM0 = G3_l_SM0(l, s, alpha0);
    double AFB_Born0 = 3.0/4.0*betaf*G3SM0/(G1SM0 + 2.0*mf2/s*G2SM0);
    double G1NP = G1_l_NP(l, s, alpha0, ObParam_i);
    double G3NP = G3_l_NP(l, s, alpha0, ObParam_i);
    
    return ( - AFB_Born0*G1NP/(G1SM0 + 2.0*mf2/s*G2SM0)
             + AFB_Born0*G3NP/G3SM0 );
}


double LEP2oblique::AFB_q_LEP2_NP(const QCD::quark q, 
                                  const double s, const double mq,
                                  const double ObParam_i[]) const
{
    double alpha0 = alpha_at_s(s);
    double mf2 = mq*mq;
    double betaf = sqrt(1.0 - 4.0*mf2/s);
    double G1SM0 = G1_q_SM0(q, s, alpha0);
    double G2SM0 = G2_q_SM0(q, s, alpha0);
    double G3SM0 = G3_q_SM0(q, s, alpha0);
    double AFB_Born0 = 3.0/4.0*betaf*G3SM0/(G1SM0 + 2.0*mf2/s*G2SM0);
    double G1NP = G1_q_NP(q, s, alpha0, ObParam_i);
    double G3NP = G3_q_NP(q, s, alpha0, ObParam_i);
    
    return ( - AFB_Born0*G1NP/(G1SM0 + 2.0*mf2/s*G2SM0)
             + AFB_Born0*G3NP/G3SM0 );  
}


double LEP2oblique::R_q_LEP2_NP(const QCD::quark q, 
                                const double s, const double mq,
                                const double ObParam_i[]) const
{
    double alpha0 = alpha_at_s(s);
    double sigma_q_SM0 = sigma_q_LEP2_SM0(q, s, alpha0, mq);
    double sigma_had_SM0 = sigma_q_LEP2_SM0(QCD::UP, s, alpha0, mq)
                         + sigma_q_LEP2_SM0(QCD::DOWN, s, alpha0, mq)
                         + sigma_q_LEP2_SM0(QCD::CHARM, s, alpha0, mq)
                         + sigma_q_LEP2_SM0(QCD::STRANGE, s, alpha0, mq)
                         + sigma_q_LEP2_SM0(QCD::BOTTOM, s, alpha0, mq);
    double sigma_q_NP = sigma_q_LEP2_NP(q, s, mq, ObParam_i);
    double sigma_had_NP = sigma_q_LEP2_NP(QCD::UP, s, mq, ObParam_i)
                        + sigma_q_LEP2_NP(QCD::DOWN, s, mq, ObParam_i)
                        + sigma_q_LEP2_NP(QCD::CHARM, s, mq, ObParam_i)
                        + sigma_q_LEP2_NP(QCD::STRANGE, s, mq, ObParam_i)
                        + sigma_q_LEP2_NP(QCD::BOTTOM, s, mq, ObParam_i);
    
    return ( - sigma_q_SM0/(sigma_had_SM0*sigma_had_SM0)*sigma_had_NP
             + sigma_q_NP/sigma_had_SM0 );
}


////////////////////////////////////////////////////////////////////////
   
double LEP2oblique::DeltaEpsilon_1(const double alpha0, 
                                   const double ObParam_i[]) const
{
    double c0 = sqrt(c02(alpha0)), s0 = sqrt(s02(alpha0));
    return ( ObParam_i[That] - ObParam_i[W] + 2.0*s0/c0*ObParam_i[X] 
             - s0*s0/c0/c0*ObParam_i[Y] );
}


double LEP2oblique::DeltaEpsilon_2(const double alpha0, 
                                   const double ObParam_i[]) const
{
    double c0 = sqrt(c02(alpha0)), s0 = sqrt(s02(alpha0));
    return ( ObParam_i[Uhat] - ObParam_i[V] - ObParam_i[W] 
             + 2.0*s0/c0*ObParam_i[X] );
}


double LEP2oblique::DeltaEpsilon_3(const double alpha0, 
                                   const double ObParam_i[]) const
{
    double c0 = sqrt(c02(alpha0)), s0 = sqrt(s02(alpha0));
    return ( ObParam_i[Shat] - ObParam_i[W] + ObParam_i[X]/s0/c0 
             - ObParam_i[Y] );
}


double LEP2oblique::epsilonZZ(const double alpha0, 
                              const double ObParam_i[]) const 
{
    double c0 = sqrt(c02(alpha0)), s0 = sqrt(s02(alpha0));
    return ( c02(alpha0)*ObParam_i[W] - 2.0*s0*c0*ObParam_i[X] 
             + s02(alpha0)*ObParam_i[Y] );
}


double LEP2oblique::epsilonGammaGamma(const double alpha0, 
                                      const double ObParam_i[]) const 
{
    double c0 = sqrt(c02(alpha0)), s0 = sqrt(s02(alpha0));
    return ( s02(alpha0)*ObParam_i[W] + 2.0*s0*c0*ObParam_i[X] 
             + c02(alpha0)*ObParam_i[Y] );
}


double LEP2oblique::epsilonGammaZ(const double alpha0, 
                                  const double ObParam_i[]) const 
{
    double c0 = sqrt(c02(alpha0)), s0 = sqrt(s02(alpha0));
    return ( (c02(alpha0) - s02(alpha0))*ObParam_i[X] + s0*c0*(ObParam_i[W] 
             - ObParam_i[Y]) );    
}


double LEP2oblique::vl(const StandardModel::lepton l, const double alpha0) const 
{
    double c0 = sqrt(c02(alpha0)), s0 = sqrt(s02(alpha0));
    double Q = SM.getLeptons(l).getCharge();
    return ( - (SM.getLeptons(l).getIsospin() 
                - 2.0*Q*s02(alpha0))/(2.0*s0*c0) );
}


double LEP2oblique::vq(const QCD::quark q, const double alpha0) const 
{
    double c0 = sqrt(c02(alpha0)), s0 = sqrt(s02(alpha0));
    double Q = SM.getQuarks(q).getCharge();
    return ( - (SM.getQuarks(q).getIsospin() 
                - 2.0*Q*s02(alpha0))/(2.0*s0*c0) );
}


double LEP2oblique::al(const StandardModel::lepton l, const double alpha0) const
{
    double c0 = sqrt(c02(alpha0)), s0 = sqrt(s02(alpha0));
    return ( - SM.getLeptons(l).getIsospin()/(2.0*s0*c0) ); 
}


double LEP2oblique::aq(const QCD::quark q, const double alpha0) const 
{
    double c0 = sqrt(c02(alpha0)), s0 = sqrt(s02(alpha0));
    return ( - SM.getQuarks(q).getIsospin()/(2.0*s0*c0) );     
}


double LEP2oblique::G1_NP(const double s, const double alpha0, const double Qf, 
                          const double vf, const double af, 
                          const double ObParam_i[]) const 
{
    double c0 = sqrt(c02(alpha0)), s0 = sqrt(s02(alpha0));
    double Qe = SM.getLeptons(StandardModel::ELECTRON).getCharge();
    double ve = vl(StandardModel::ELECTRON, alpha0);
    double ae = al(StandardModel::ELECTRON, alpha0);
    double Qe2 = Qe*Qe, Qf2 = Qf*Qf;
    double ve2 = ve*ve, vf2 = vf*vf, ae2 = ae*ae, af2 = af*af;
    double Mz = SM.getMz();
    double GammaZ0 = 7.0*alpha0*Mz/(16.0*s02(alpha0)*c02(alpha0));
    gslpp::complex denom = gslpp::complex(s - Mz*Mz, Mz*GammaZ0, false);
    double Zprop = (1.0/denom).real();

    double epsilonbarGamma = - s/(Mw0(alpha0)*Mw0(alpha0))
                               *epsilonGammaGamma(alpha0, ObParam_i);
    double epsilonbarZ = s*Zprop*DeltaEpsilon_1(alpha0, ObParam_i)
                         - s/(Mw0(alpha0)*Mw0(alpha0))*epsilonZZ(alpha0, ObParam_i);
    double epsilonbarGammaZ = c0/s0*s*Zprop
                              *( DeltaEpsilon_1(alpha0, ObParam_i) 
                                 - DeltaEpsilon_2(alpha0, ObParam_i) )
                              - s0/c0*s*Zprop*DeltaEpsilon_3(alpha0, ObParam_i)
                              + s/(Mw0(alpha0)*Mw0(alpha0))
                                *epsilonGammaZ(alpha0, ObParam_i);
    
    return ( 2.0*Qe2*Qf2*epsilonbarGamma
             + 2.0*ve*vf*Qe*Qf*(epsilonbarZ + s*epsilonbarGamma*Zprop)
             + 2.0*(ve2 + ae2)*(vf2 + af2)*s*epsilonbarZ*Zprop
             + 2.0*(vf*Qe2*Qf + ve*Qe*Qf2)*epsilonbarGammaZ
             + 2.0*(ve*(vf2 + af2)*Qe + (ve2 + ae2)*vf*Qf)*s*epsilonbarGammaZ*Zprop );    
}


double LEP2oblique::G1_l_NP(const StandardModel::lepton l, 
                            const double s, const double alpha0, 
                            const double ObParam_i[]) const 
{
    double Qf = SM.getLeptons(l).getCharge();
    double vf = vl(l, alpha0), af = al(l, alpha0);
    return ( G1_NP(s, alpha0, Qf, vf, af, ObParam_i) );
}


double LEP2oblique::G1_q_NP(const QCD::quark q, 
                            const double s, const double alpha0, 
                            const double ObParam_i[]) const
{
    double Qf = SM.getQuarks(q).getCharge();
    double vf = vq(q, alpha0), af = aq(q, alpha0);    
    return ( G1_NP(s, alpha0, Qf, vf, af, ObParam_i) );    
}


double LEP2oblique::G3_NP(const double s, const double alpha0, const double Qf, 
                          const double vf, const double af, 
                          const double ObParam_i[]) const
{
    double c0 = sqrt(c02(alpha0)), s0 = sqrt(s02(alpha0));
    double Qe = SM.getLeptons(StandardModel::ELECTRON).getCharge();
    double ve = vl(StandardModel::ELECTRON, alpha0);
    double ae = al(StandardModel::ELECTRON, alpha0);
    double Mz = SM.getMz();
    double GammaZ0 = 7.0*alpha0*Mz/(16.0*s02(alpha0)*c02(alpha0));
    gslpp::complex denom = gslpp::complex(s - Mz*Mz, Mz*GammaZ0, false);
    double Zprop = (1.0/denom).real();
    
    double epsilonbarGamma = - s/(Mw0(alpha0)*Mw0(alpha0))
                               *epsilonGammaGamma(alpha0, ObParam_i);
    double epsilonbarZ = s*Zprop*DeltaEpsilon_1(alpha0, ObParam_i)
                         - s/(Mw0(alpha0)*Mw0(alpha0))*epsilonZZ(alpha0, ObParam_i);
    double epsilonbarGammaZ = c0/s0*s*Zprop*( DeltaEpsilon_1(alpha0, ObParam_i) 
                                              - DeltaEpsilon_2(alpha0, ObParam_i) )
                              - s0/c0*s*Zprop*DeltaEpsilon_3(alpha0, ObParam_i)
                              + s/(Mw0(alpha0)*Mw0(alpha0))
                                *epsilonGammaZ(alpha0, ObParam_i);
    
    return ( 2.0*ae*af*Qe*Qf*(epsilonbarZ + s*epsilonbarGamma*Zprop)
             + 8.0*ve*ae*vf*af*s*epsilonbarZ*Zprop
             + 4.0*(ae*vf*af*Qe + ve*ae*af*Qf)*s*epsilonbarGammaZ*Zprop );   
}


double LEP2oblique::G3_l_NP(const StandardModel::lepton l, 
                            const double s, const double alpha0, 
                            const double ObParam_i[]) const
{
    double Qf = SM.getLeptons(l).getCharge();
    double vf = vl(l, alpha0), af = al(l, alpha0);    
    return ( G3_NP(s, alpha0, Qf, vf, af, ObParam_i) );
}


double LEP2oblique::G3_q_NP(const QCD::quark q, 
                            const double s, const double alpha0, 
                            const double ObParam_i[]) const
{
    double Qf = SM.getQuarks(q).getCharge();
    double vf = vq(q, alpha0), af = aq(q, alpha0);    
    return ( G3_NP(s, alpha0, Qf, vf, af, ObParam_i) ); 
}


double LEP2oblique::G1_SM0(const double s, const double alpha0, const double Qf, 
                           const double vf, const double af) const
{
    double Qe = SM.getLeptons(StandardModel::ELECTRON).getCharge();
    double ve = vl(StandardModel::ELECTRON, alpha0);
    double ae = al(StandardModel::ELECTRON, alpha0);
    double Qe2 = Qe*Qe, Qf2 = Qf*Qf;
    double ve2 = ve*ve, vf2 = vf*vf, ae2 = ae*ae, af2 = af*af;
    double Mz = SM.getMz();
    double GammaZ0 = 7.0*alpha0*Mz/(16.0*s02(alpha0)*c02(alpha0));
    gslpp::complex denom = gslpp::complex(s - Mz*Mz, Mz*GammaZ0, false);
    gslpp::complex chiZ = s/denom;
    
    return ( Qe2*Qf2 + 2.0*ve*vf*Qe*Qf*chiZ.real()
             + (ve2 + ae2)*(vf2 + af2)*chiZ.abs2() );
}    


double LEP2oblique::G1_l_SM0(const StandardModel::lepton l, 
                             const double s, const double alpha0) const 
{
    double Qf = SM.getLeptons(l).getCharge();
    double vf = vl(l, alpha0), af = al(l, alpha0);
    return ( G1_SM0(s, alpha0, Qf, vf, af) );
}


double LEP2oblique::G1_q_SM0(const QCD::quark q, 
                             const double s, const double alpha0) const 
{
    double Qf = SM.getQuarks(q).getCharge();
    double vf = vq(q, alpha0), af = aq(q, alpha0);
    return ( G1_SM0(s, alpha0, Qf, vf, af) ); 
}


double LEP2oblique::G2_SM0(const double s, const double alpha0, const double Qf, 
                           const double vf, const double af) const 
{
    double Qe = SM.getLeptons(StandardModel::ELECTRON).getCharge();
    double ve = vl(StandardModel::ELECTRON, alpha0);
    double ae = al(StandardModel::ELECTRON, alpha0);
    double Qe2 = Qe*Qe, Qf2 = Qf*Qf;
    double ve2 = ve*ve, vf2 = vf*vf, ae2 = ae*ae;
    double Mz = SM.getMz();
    double GammaZ0 = 7.0*alpha0*Mz/(16.0*s02(alpha0)*c02(alpha0));
    gslpp::complex denom = gslpp::complex(s - Mz*Mz, Mz*GammaZ0, false);
    gslpp::complex chiZ = s/denom;
    
    return ( Qe2*Qf2 + 2.0*ve*vf*Qe*Qf*chiZ.real()
             + (ve2 + ae2)*vf2*chiZ.abs2() );
}


double LEP2oblique::G2_l_SM0(const StandardModel::lepton l, 
                             const double s, const double alpha0) const
{
    double Qf = SM.getLeptons(l).getCharge();
    double vf = vl(l, alpha0), af = al(l, alpha0);
    return ( G2_SM0(s, alpha0, Qf, vf, af) );    
}


double LEP2oblique::G2_q_SM0(const QCD::quark q, 
                             const double s, const double alpha0) const
{
    double Qf = SM.getQuarks(q).getCharge();
    double vf = vq(q, alpha0), af = aq(q, alpha0);
    return ( G2_SM0(s, alpha0, Qf, vf, af) ); 
    
}


double LEP2oblique::G3_SM0(const double s, const double alpha0, const double Qf, 
                           const double vf, const double af) const 
{
    double Qe = SM.getLeptons(StandardModel::ELECTRON).getCharge();
    double ve = vl(StandardModel::ELECTRON, alpha0);
    double ae = al(StandardModel::ELECTRON, alpha0);
    double Mz = SM.getMz();
    double GammaZ0 = 7.0*alpha0*Mz/(16.0*s02(alpha0)*c02(alpha0));
    gslpp::complex denom = gslpp::complex(s - Mz*Mz, Mz*GammaZ0, false);
    gslpp::complex chiZ = s/denom;

    return ( 2.0*ae*af*Qe*Qf*chiZ.real() + 4.0*ve*ae*vf*af*chiZ.abs2() );
}


double LEP2oblique::G3_l_SM0(const StandardModel::lepton l, 
                             const double s, const double alpha0) const 
{
    double Qf = SM.getLeptons(l).getCharge();
    double vf = vl(l, alpha0), af = al(l, alpha0);
    return ( G3_SM0(s, alpha0, Qf, vf, af) );
}


double LEP2oblique::G3_q_SM0(const QCD::quark q, 
                             const double s, const double alpha0) const 
{
    double Qf = SM.getQuarks(q).getCharge();
    double vf = vq(q, alpha0), af = aq(q, alpha0);
    return ( G3_SM0(s, alpha0, Qf, vf, af) );    
}


double LEP2oblique::sigma_l_LEP2_SM0(const StandardModel::lepton l, 
                                     const double s, const double alpha0, 
                                     const double ml) const 
{
    double Ncf = 1.0;
    double mf2 = ml*ml;
    double betaf = sqrt(1.0 - 4.0*mf2/s);
    double G1SM0 = G1_l_SM0(l, s, alpha0), G2SM0 = G2_l_SM0(l, s, alpha0);

    return ( 4.0*M_PI*alpha0*alpha0/(3.0*s)*Ncf*betaf
             *(G1SM0 + 2.0*mf2/s*G2SM0) );    
}


double LEP2oblique::sigma_q_LEP2_SM0(const QCD::quark q, 
                                     const double s, const double alpha0, 
                                     const double mq) const 
{
    double Ncf = 3.0;
    double mf2 = mq*mq;
    double betaf = sqrt(1.0 - 4.0*mf2/s);
    double G1SM0 = G1_q_SM0(q, s, alpha0), G2SM0 = G2_q_SM0(q, s, alpha0);

    return ( 4.0*M_PI*alpha0*alpha0/(3.0*s)*Ncf*betaf
             *(G1SM0 + 2.0*mf2/s*G2SM0) );
}


