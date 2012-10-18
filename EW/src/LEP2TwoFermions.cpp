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

    return sigma(s, Mw, GammaZ, I3f, Qf, mfp, mf, 3.0, bRCs);
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

    return AFB(s, Mw, GammaZ, I3f, Qf, mfp, mf, bRCs);
}


////////////////////////////////////////////////////////////////////////     

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


////////////////////////////////////////////////////////////////////////     

double LEP2TwoFermions::sigma(const double s, const double Mw, 
                              const double GammaZ, const double I3f, 
                              const double Qf, const double mfp,
                              const double mf, const double Ncf, 
                              const bool bRCs[]) const {
    double betaf = sqrt(1.0 - 4.0*mf*mf/s);
    double G1 = SM.G_1(s, Mw, GammaZ, I3f, Qf, mfp, bRCs[Weak]);
    double G2 = SM.G_2(s, Mw, GammaZ, I3f, Qf, mfp, bRCs[Weak]);

    double FSR = 1.0;
    if(bRCs[QEDFSR]) {
        double alpha;
        if(bDebug)
            alpha = SM.getAle()/complex(1.0715119759, -0.0186242179, false).real(); // for debug
        else
            alpha = SM.alphaMz(); //!!TEST
        FSR += 3.0*alpha/(4.0*M_PI)*Qf*Qf;
    }
    
    return ( 4.0*M_PI*SM.getAle()*SM.getAle()/(3.0*s)*Ncf*betaf
             *( G1 + 2.0*mf*mf/s*G2)*FSR );
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




