/* 
 * Copyright (C) 2012 HEPfit Collaboration
 * All rights reserved.
 *
 * For the licensing terms see doc/COPYING.
 */

#include "EWSMTwoFermionsLEP2.h"
#include "LEP2TwoFermions.h"


LEP2TwoFermions::LEP2TwoFermions(const StandardModel& SM_i) 
: SM(SM_i) 
{
}


//////////////////////////////////////////////////////////////////////// 

double LEP2TwoFermions::dsigma_l(const StandardModel::lepton l, const double mf, 
                                 const double s, const double cosTheta,
                                 const double Mw, const double GammaZ, 
                                 const bool bWeak) const 
{
    double I3f = SM.getLeptons(l).getIsospin();
    double Qf = SM.getLeptons(l).getCharge();

    return ( dsigma(s, cosTheta, Mw, GammaZ, I3f, Qf, mf, 0.0, 1.0, bWeak) );
}


double LEP2TwoFermions::dsigma_q(const QCD::quark q, const double mf, 
                                 const double s, const double cosTheta,
                                 const double Mw, const double GammaZ, 
                                 const bool bWeak) const 
{
    double I3f = SM.getQuarks(q).getIsospin();
    double Qf = SM.getQuarks(q).getCharge();
    double mfp;
    if (q==SM.TOP)
        throw std::runtime_error("Error in LEP2TwoFermions::sigma_q()");
    else if (q==SM.BOTTOM)
        mfp = SM.getMtpole();
    else
        mfp = 0.0;
    
    return ( dsigma(s, cosTheta, Mw, GammaZ, I3f, Qf, mf, mfp, 3.0, bWeak) );
}


double LEP2TwoFermions::dsigma_l_box(const StandardModel::lepton l, const double mf, 
                                     const double s, const double cosTheta,
                                     const double Mw, const double GammaZ) const 
{
    double I3f = SM.getLeptons(l).getIsospin();
    double Qf = SM.getLeptons(l).getCharge();

    return ( dsigma_box(s, cosTheta, Mw, GammaZ, I3f, Qf, mf, 0.0, 1.0) );
}


double LEP2TwoFermions::dsigma_q_box(const QCD::quark q, const double mf,
                                     const double s, const double cosTheta,
                                     const double Mw, const double GammaZ) const 
{
    double I3f = SM.getQuarks(q).getIsospin();
    double Qf = SM.getQuarks(q).getCharge();
    double mfp;
    if (q==SM.TOP)
        throw std::runtime_error("Error in LEP2TwoFermions::sigma_q()");
    else if (q==SM.BOTTOM)
        mfp = SM.getMtpole();
    else
        mfp = 0.0;
    
    return ( dsigma_box(s, cosTheta, Mw, GammaZ, I3f, Qf, mf, mfp, 3.0) );
}


double LEP2TwoFermions::sigma_l(const StandardModel::lepton l, const double mf, 
                                const double s, const double Mw, 
                                const double GammaZ, const bool bWeak) const 
{
    double I3f = SM.getLeptons(l).getIsospin();
    double Qf = SM.getLeptons(l).getCharge();

    return ( sigma(s, Mw, GammaZ, I3f, Qf, mf, 0.0, 1.0, bWeak) );
}


double LEP2TwoFermions::sigma_q(const QCD::quark q, const double mf, 
                                const double s, const double Mw, 
                                const double GammaZ, const bool bWeak) const 
{
    double I3f = SM.getQuarks(q).getIsospin();
    double Qf = SM.getQuarks(q).getCharge();
    double mfp;
    if (q==SM.TOP)
        throw std::runtime_error("Error in LEP2TwoFermions::sigma_q()");
    else if (q==SM.BOTTOM)
        mfp = SM.getMtpole();
    else
        mfp = 0.0;
    
    return ( sigma(s, Mw, GammaZ, I3f, Qf, mf, mfp, 3.0, bWeak) );
}


double LEP2TwoFermions::AFB_l(const StandardModel::lepton l, const double mf, 
                              const double s, const double Mw, 
                              const double GammaZ, const bool bWeak) const 
{
    double I3f = SM.getLeptons(l).getIsospin();
    double Qf = SM.getLeptons(l).getCharge();

    return AFB(s, Mw, GammaZ, I3f, Qf, mf, 0.0, bWeak);
}


double LEP2TwoFermions::AFB_q(const QCD::quark q, const double mf, 
                              const double s, const double Mw, 
                              const double GammaZ, const bool bWeak) const 
{
    double I3f = SM.getQuarks(q).getIsospin();
    double Qf = SM.getQuarks(q).getCharge();
    double mfp;
    if (q==SM.TOP)
        throw std::runtime_error("Error in LEP2TwoFermions::AFB_q()");
    else if (q==SM.BOTTOM)
        mfp = SM.getMtpole();
    else
        mfp = 0.0;
    
    return ( AFB(s, Mw, GammaZ, I3f, Qf, mf, mfp, bWeak) );
}


double LEP2TwoFermions::QCD_FSR_forSigma(const double s) const 
{
    return ( 1.0 + SM.Als(sqrt(s), FULLNLO)/M_PI );
}
    

double LEP2TwoFermions::QCD_FSR_forAFB(const QCD::quark q, 
                                       const double mf, const double s) const 
{
    return ( 1.0 - SM.Als(sqrt(s), FULLNLO)/M_PI*(1.0 - 16.0/3.0*mf/sqrt(s)) );
}


double LEP2TwoFermions::QED_FSR_forSigma(const double s, const double Qf) const
{
    //double alpha = SM.getAle();
    double alpha = alpha_at_s(s);
    
    return ( 1.0 + 3.0*alpha/(4.0*M_PI)*Qf*Qf );
}


double LEP2TwoFermions::H_ISR(const double x, const double s) const 
{
    double me = SM.getLeptons(SM.ELECTRON).getMass();
    double alphaOverPi = SM.getAle()/M_PI; // alpha(0)/Pi
    double L = log(s/(me*me));
    double beta = 2.0*alphaOverPi*(L - 1.0);
    double deltaVS_1 = 3.0/2.0*L + M_PI*M_PI/3.0 - 2.0;
    double deltaH_1 = - (2.0 - x)*(L - 1.0);   
    
    return ( beta*pow(x, beta-1.0)*(1.0 + alphaOverPi*deltaVS_1) 
             + alphaOverPi*deltaH_1 );
}


double LEP2TwoFermions::H_ISR_FB(const double x, const double s) const 
{
    double me = SM.getLeptons(SM.ELECTRON).getMass();
    double alphaOverPi = SM.getAle()/M_PI; // alpha(0)/Pi
    double L = log(s/(me*me));
    double beta = 2.0*alphaOverPi*(L - 1.0);
    double deltaVS_1 = 3.0/2.0*L + M_PI*M_PI/3.0 - 2.0;
    double tmp = (1.0-x)/(1.0-x/2.0)/(1.0-x/2.0);
    double deltaH_FB_1 = (1.0+(1.0-x)*(1.0-x))/x*tmp*(L - 1.0 - log(tmp))
                         - 2.0/x*(L - 1.0);
    
    return ( beta*pow(x, beta-1.0)*(1.0 + alphaOverPi*deltaVS_1) 
             + alphaOverPi*deltaH_FB_1 );
}


double LEP2TwoFermions::G_3prime_l(const StandardModel::lepton l, 
                                   const double mf, const double s,
                                   const double Mw, const double GammaZ, 
                                   const bool bWeak) const 
{
    double betaf = sqrt(1.0 - 4.0*mf*mf/s);
    double I3f = SM.getLeptons(l).getIsospin();
    double Qf = SM.getLeptons(l).getCharge();
    double G3 = SM.getMyTwoFermionsLEP2()->G_3_noBox(s, Mw, GammaZ, I3f, Qf, mf, 0.0, bWeak);    
    
    return ( betaf*betaf*G3 );    
}


double LEP2TwoFermions::G_3prime_q(const QCD::quark q, 
                                   const double mf, const double s,
                                   const double Mw, const double GammaZ, 
                                   const bool bWeak) const 
{
    double betaf = sqrt(1.0 - 4.0*mf*mf/s);
    double I3f = SM.getQuarks(q).getIsospin();
    double Qf = SM.getQuarks(q).getCharge();
    double mfp;
    if (q==SM.TOP)
        throw std::runtime_error("Error in LEP2TwoFermions::G_3prime_q()");
    else if (q==SM.BOTTOM)
        mfp = SM.getMtpole();
    else
        mfp = 0.0;
    double G3 = SM.getMyTwoFermionsLEP2()->G_3_noBox(s, Mw, GammaZ, I3f, Qf, mf, mfp, bWeak);
    
    return ( betaf*betaf*G3 );    
}


////////////////////////////////////////////////////////////////////////     

double LEP2TwoFermions::alpha_at_s(const double s) const 
{
    double alpha;
    //alpha = SM.getAle()/complex(1.0715119759, -0.0186242179, false).real(); // for debug, s=(200GeV)^2
    alpha = SM.ale_OS(sqrt(s), FULLNLO);

    return alpha;
}


double LEP2TwoFermions::dsigma(const double s, const double cosTheta, 
                               const double Mw, const double GammaZ, 
                               const double I3f, const double Qf, 
                               const double mf, const double mfp, 
                               const double Ncf, const bool bWeak) const 
{
    double betaf = sqrt(1.0 - 4.0*mf*mf/s);
    double G1 = SM.getMyTwoFermionsLEP2()->G_1_noBox(s, Mw, GammaZ, I3f, Qf, mf, mfp, bWeak);
    double G2 = SM.getMyTwoFermionsLEP2()->G_2_noBox(s, Mw, GammaZ, I3f, Qf, mf, mfp, bWeak);
    double G3 = SM.getMyTwoFermionsLEP2()->G_3_noBox(s, Mw, GammaZ, I3f, Qf, mf, mfp, bWeak);

    return ( M_PI*SM.getAle()*SM.getAle()/s*Ncf*betaf
             *( G1*(1.0 + cosTheta*cosTheta)/2.0 
                + 2.0*mf*mf/s*G2*(1.0 - cosTheta*cosTheta)
                + betaf*G3*cosTheta ) );
}


double LEP2TwoFermions::dsigma_box(const double s, const double cosTheta, 
                                   const double Mw, const double GammaZ, 
                                   const double I3f, const double Qf, 
                                   const double mf, const double mfp, 
                                   const double Ncf) const
{
    double betaf = sqrt(1.0 - 4.0*mf*mf/s);
    
    //double t = mf*mf - s/2.0*(1.0 - betaf*cosTheta);
    double t = - s/2.0*(1.0 - betaf*cosTheta);
    
    double G1 = SM.getMyTwoFermionsLEP2()->G_1_box(s, t, Mw, GammaZ, I3f, Qf, mf, mfp);
    double G2 = SM.getMyTwoFermionsLEP2()->G_2_box(s, t, Mw, GammaZ, I3f, Qf, mf, mfp);
    double G3 = SM.getMyTwoFermionsLEP2()->G_3_box(s, t, Mw, GammaZ, I3f, Qf, mf, mfp);

    return ( M_PI*SM.getAle()*SM.getAle()/s*Ncf*betaf
             *( G1*(1.0 + cosTheta*cosTheta)/2.0 
                + 2.0*mf*mf/s*G2*(1.0 - cosTheta*cosTheta)
                + betaf*G3*cosTheta ) );
}


double LEP2TwoFermions::sigma(const double s, const double Mw, 
                              const double GammaZ, const double I3f, 
                              const double Qf, const double mf,
                              const double mfp, const double Ncf, 
                              const bool bWeak) const 
{
    double betaf = sqrt(1.0 - 4.0*mf*mf/s);
    double G1 = SM.getMyTwoFermionsLEP2()->G_1_noBox(s, Mw, GammaZ, I3f, Qf, mf, mfp, bWeak);
    double G2 = SM.getMyTwoFermionsLEP2()->G_2_noBox(s, Mw, GammaZ, I3f, Qf, mf, mfp, bWeak);

    return ( 4.0*M_PI*SM.getAle()*SM.getAle()/(3.0*s)*Ncf*betaf
             *(G1 + 2.0*mf*mf/s*G2) );
}


double LEP2TwoFermions::AFB(const double s, const double Mw, 
                            const double GammaZ, const double I3f, 
                            const double Qf, const double mf,
                            const double mfp, const bool bWeak) const 
{
    double betaf = sqrt(1.0 - 4.0*mf*mf/s);
    double G1 = SM.getMyTwoFermionsLEP2()->G_1_noBox(s, Mw, GammaZ, I3f, Qf, mf, mfp, bWeak);
    double G2 = SM.getMyTwoFermionsLEP2()->G_2_noBox(s, Mw, GammaZ, I3f, Qf, mf, mfp, bWeak);
    double G3 = SM.getMyTwoFermionsLEP2()->G_3_noBox(s, Mw, GammaZ, I3f, Qf, mf, mfp, bWeak);

    return ( 3.0/4.0*betaf*G3/(G1 + 2.0*mf*mf/s*G2) );
}

