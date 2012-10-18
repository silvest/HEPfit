/* 
 * File:   LEP2TwoFermions.cpp
 * Author: mishima
 */

#include "LEP2TwoFermions.h"


LEP2TwoFermions::LEP2TwoFermions(const StandardModel& SM_i) : EW(SM_i) {
}


//////////////////////////////////////////////////////////////////////// 

double LEP2TwoFermions::sigma_l(const StandardModel::lepton l, 
                                const double s, const double Mw, 
                                const double GammaZ, const bool bRCs[]) const {
    double mf = ml(l);
    double I3f = SM.getLeptons(l).getIsospin();
    double Qf = SM.getLeptons(l).getCharge();

    return sigma(s, Mw, GammaZ, I3f, Qf, 0.0, mf, 1.0, bRCs);
}


double LEP2TwoFermions::sigma_q(const StandardModel::quark q, 
                                const double s, const double Mw, 
                                const double GammaZ, const bool bRCs[]) const {
    double mf = mq(q, sqrt(s));
    double I3f = SM.getQuarks(q).getIsospin();
    double Qf = SM.getQuarks(q).getCharge();
    double mfp;
    if (q==SM.TOP)
        throw std::runtime_error("Error in LEP2TwoFermions::sigma_q()");
    else if (q==SM.BOTTOM) 
        mfp = mq(SM.TOP, sqrt(s));
    else 
        mfp = 0.0;

    // Final-state QCD radiations
    double QCD_FSR = 1.0;
    if(bRCs[QCDFSR])
        QCD_FSR += SM.Als(sqrt(s), FULLNLO)/M_PI;
    
    return ( sigma(s, Mw, GammaZ, I3f, Qf, mfp, mf, 3.0, bRCs)*QCD_FSR );
}


double LEP2TwoFermions::AFB_l(const StandardModel::lepton l, 
                              const double s, const double Mw, 
                              const double GammaZ, const bool bRCs[]) const {
    double mf = ml(l);
    double I3f = SM.getLeptons(l).getIsospin();
    double Qf = SM.getLeptons(l).getCharge();

    return AFB(s, Mw, GammaZ, I3f, Qf, 0.0, mf, bRCs);
}


double LEP2TwoFermions::AFB_q(const StandardModel::quark q, 
                              const double s, const double Mw, 
                              const double GammaZ, const bool bRCs[]) const {
    double mf = mq(q, sqrt(s));
    double I3f = SM.getQuarks(q).getIsospin();
    double Qf = SM.getQuarks(q).getCharge();
    double mfp;
    if (q==SM.TOP)
        throw std::runtime_error("Error in LEP2TwoFermions::AFB_q()");
    else if (q==SM.BOTTOM) 
        mfp = mq(SM.TOP, sqrt(s));
    else 
        mfp = 0.0;

    // Final-state QCD radiations
    double QCD_FSR = 1.0;
    if(bRCs[QCDFSR])
        QCD_FSR -= SM.Als(sqrt(s), FULLNLO)/M_PI*(1.0 - 16.0/3.0*mf/sqrt(s));
    
    return ( AFB(s, Mw, GammaZ, I3f, Qf, mfp, mf, bRCs)*QCD_FSR );
}


double LEP2TwoFermions::H_ISR(const double x, const double s) const {
    double me = SM.getLeptons(SM.ELECTRON).getMass();
    double alphaOverPi = SM.getAle()/M_PI; // alpha(0)/Pi
    double L = log(s/(me*me));
    double beta = 2.0*alphaOverPi*(L - 1.0);
    double deltaVS_1 = 3.0/2.0*L + M_PI*M_PI/3.0 - 2.0;
    double deltaH_1 = - (2.0 - x)*(L - 1.0);   
    
    return ( beta*pow(x, beta-1.0)*(1.0 + alphaOverPi*deltaVS_1) 
             + alphaOverPi*deltaH_1 );
}


double LEP2TwoFermions::H_ISR_FB(const double x, const double s) const {
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


double LEP2TwoFermions::G_3prime_l(const StandardModel::lepton l, const double s, 
                                   const double Mw, const double GammaZ, 
                                   const bool bRCs[]) const {
    double mf = ml(l);
    double betaf = sqrt(1.0 - 4.0*mf*mf/s);
    double I3f = SM.getLeptons(l).getIsospin();
    double Qf = SM.getLeptons(l).getCharge();
    double G3 = SM.G_3(s, Mw, GammaZ, I3f, Qf, 0.0, bRCs[Weak]);    
    
    return ( betaf*betaf*G3 );    
}


double LEP2TwoFermions::G_3prime_q(const StandardModel::quark q, const double s, 
                                   const double Mw, const double GammaZ, 
                                   const bool bRCs[]) const {
    double mf = mq(q, sqrt(s));
    double betaf = sqrt(1.0 - 4.0*mf*mf/s);
    double I3f = SM.getQuarks(q).getIsospin();
    double Qf = SM.getQuarks(q).getCharge();
    double mfp;
    if (q==SM.TOP)
        throw std::runtime_error("Error in LEP2TwoFermions::G_3prime_q()");
    else if (q==SM.BOTTOM) 
        mfp = mq(SM.TOP, sqrt(s));
    else 
        mfp = 0.0;
    double G3 = SM.G_3(s, Mw, GammaZ, I3f, Qf, mfp, bRCs[Weak]);    
    
    return ( betaf*betaf*G3 );    
}


////////////////////////////////////////////////////////////////////////     

double LEP2TwoFermions::alpha_at_s(const double s) const {
    double alpha;
    if(bDebug)
        alpha = SM.getAle()/complex(1.0715119759, -0.0186242179, false).real(); // for debug, s=(200GeV)^2
    else
        alpha = SM.alphaMz(); //!!TEST    

    return alpha;
}


double LEP2TwoFermions::ml(const StandardModel::lepton l) const {
    return SM.getLeptons(l).getMass();
}


double LEP2TwoFermions::mq(const StandardModel::quark q, const double mu, 
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
            throw std::runtime_error("Error in LEP2TwoFermions::mq()"); 
    }
}


double LEP2TwoFermions::sigma(const double s, const double Mw, 
                              const double GammaZ, const double I3f, 
                              const double Qf, const double mfp,
                              const double mf, const double Ncf, 
                              const bool bRCs[]) const {
    double betaf = sqrt(1.0 - 4.0*mf*mf/s);
    double G1 = SM.G_1(s, Mw, GammaZ, I3f, Qf, mfp, bRCs[Weak]);
    double G2 = SM.G_2(s, Mw, GammaZ, I3f, Qf, mfp, bRCs[Weak]);

    // Final-state QED radiations    
    double QED_FSR = 1.0;
    if(bRCs[QEDFSR])
        QED_FSR += 3.0*SM.getAle()/(4.0*M_PI)*Qf*Qf;
    
    return ( 4.0*M_PI*SM.getAle()*SM.getAle()/(3.0*s)*Ncf*betaf
             *( G1 + 2.0*mf*mf/s*G2)*QED_FSR );
}


double LEP2TwoFermions::AFB(const double s, const double Mw, 
                            const double GammaZ, const double I3f, 
                            const double Qf, const double mfp,
                            const double mf, const bool bRCs[]) const {
    double betaf = sqrt(1.0 - 4.0*mf*mf/s);
    double G1 = SM.G_1(s, Mw, GammaZ, I3f, Qf, mfp, bRCs[Weak]);
    double G2 = SM.G_2(s, Mw, GammaZ, I3f, Qf, mfp, bRCs[Weak]);
    double G3 = SM.G_3(s, Mw, GammaZ, I3f, Qf, mfp, bRCs[Weak]);

    return ( 3.0/4.0*betaf*G3/( G1 + 2.0*mf*mf/s*G2) );
}




