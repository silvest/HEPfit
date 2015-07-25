/* 
 * Copyright (C) 2012 HEPfit Collaboration
 * All rights reserved.
 *
 * For the licensing terms see doc/COPYING.
 */

#ifndef LEP2THOBSERVABLE_H
#define	LEP2THOBSERVABLE_H

#include <stdexcept>
#include <cstring>
#include <ThObservable.h>
#include <NPSTUVWXY.h>
#include <StandardModel.h>
#include "LEP2TwoFermions.h"
#include "LEP2oblique.h"
#include "LEP2test.h"

// TEST: use the Zfitter outputs defined in LEP2test class for SM predictions
//#define LEP2TEST

/**
 * @class LEP2ThObservable
 * @ingroup EW 
 * @brief A class for the LEP-II observables. 
 * @author HEPfit Collaboration
 * @copyright GNU General Public License
 * @details  
 */
class LEP2ThObservable : public ThObservable  {
public:

    // Radiative Corrections for the LEP-II observables
    enum LEP2RCs {Weak=0, WeakBox, ISR, QEDFSR, QCDFSR, NUMofLEP2RCs};
        
    /**
     * @brief LEP2ThObservable constructor
     * @param[in] SM_i a reference to an object of type StandardModel
     * @param[in] sqrt_s_i the CM energy of the e^+ e^- pair
     * @param[in] bSigmaForAFB_i true for the denominator of A_FB
     * @param[in] bSigmaForR_i true for the denominator of R_b or R_c
     */
    LEP2ThObservable(const StandardModel& SM_i, const double sqrt_s_i, 
                     const bool bSigmaForAFB_i=false,
                     const bool bSigmaForR_i=false) 
            : ThObservable(SM_i),
            myTwoFermions(SM_i), myLEP2oblique(SM_i), 
            sqrt_s(sqrt_s_i), s(sqrt_s_i*sqrt_s_i), bSigmaForAFB(bSigmaForAFB_i),
            bSigmaForR(bSigmaForR_i)
    {
        flag[Weak] = true;
        flag[WeakBox] = true;
        flag[ISR] = true;
        flag[QEDFSR] = true;
        flag[QCDFSR] = true;
        
        // test
        //flag[Weak] = false;
        //flag[WeakBox] = false;
        //flag[ISR] = false;
        //flag[QEDFSR] = false;
        //flag[QCDFSR] = false; 
    }

    bool checkLEP2NP() const
    {
        std::string Model = SM.ModelName();
        if ( (Model.compare("NPSTU") == 0
                || Model.compare("NPSTUVWXY") == 0
                || Model.compare("NPHiggs") == 0
                || Model.compare("THDM") == 0)
                && ((dynamic_cast<const NPSTUVWXY*> (&SM))->obliqueShat()!=0.0
                || (dynamic_cast<const NPSTUVWXY*> (&SM))->obliqueThat()!=0.0
                || (dynamic_cast<const NPSTUVWXY*> (&SM))->obliqueUhat()!=0.0
                || (dynamic_cast<const NPSTUVWXY*> (&SM))->obliqueV()!=0.0
                || (dynamic_cast<const NPSTUVWXY*> (&SM))->obliqueW()!=0.0
                || (dynamic_cast<const NPSTUVWXY*> (&SM))->obliqueX()!=0.0
                || (dynamic_cast<const NPSTUVWXY*> (&SM))->obliqueY()!=0.0) )
            return true;
        else
            return false;
    }    

    /**
     * @brief set a flag to control radiative corrections
     * @param[in] str "Weak", "WeakBox", "ISR", "QEDFSR" or "QCDFSR"
     * @param[in] flag_i boolean variable
     */
    void setFlag(const std::string str, const bool flag_i) 
    {
        if (str=="Weak")
            flag[Weak] = flag_i;
        else if (str=="WeakBox")
            flag[WeakBox] = flag_i;
        else if (str=="ISR")
            flag[ISR] = flag_i;
        else if (str=="QEDFSR")
            flag[QEDFSR] = flag_i;
        else if (str=="QCDFSR")
            flag[QCDFSR] = flag_i;
        else
            throw std::runtime_error("Error in LEP2ThObservable::setFlag()");
    } 

    void setFlags(const bool flag_i[])
    {
        setFlag("Weak", flag_i[Weak]);
        setFlag("WeakBox", flag_i[WeakBox]);
        setFlag("ISR", flag_i[ISR]);
        setFlag("QEDFSR", flag_i[QEDFSR]);
        setFlag("QCDFSR", flag_i[QCDFSR]);        
    }
    
    bool checkSMparams(const double s, const double Mw, const double GammaZ) const 
    {
        // 23 SM parameters in checkSMparams() + s, Mw, GammaZ + 5 booleans
        //bool bCache = true;
        //bCache &= SM.checkSMparams(SMparams_cache);
        //
        //if (SMparams_cache[StandardModel::NumSMParamsForEWPO] != s) {
        //    SMparams_cache[StandardModel::NumSMParamsForEWPO] = s;
        //    bCache &= false;
        //}
        //if (SMparams_cache[StandardModel::NumSMParamsForEWPO+1] != Mw) {
        //    SMparams_cache[StandardModel::NumSMParamsForEWPO+1] = Mw;
        //    bCache &= false;
        //}
        //if (SMparams_cache[StandardModel::NumSMParamsForEWPO+2] != GammaZ) {
        //    SMparams_cache[StandardModel::NumSMParamsForEWPO+2] = GammaZ;
        //    bCache &= false;
        //}
        //for (int i=0; i<NUMofLEP2RCs; i++) {
        //    if (flag_cache[i] != flag[i]) {
        //        flag_cache[i] = flag[i];
        //        bCache &= false;
        //    }
        //}
        //
        //return bCache;

        return false;
    }
    
    
protected:
    // These variables have to be initialized in child classes. 
    StandardModel::lepton l_flavor;
    QCD::quark q_flavor;
    
    const LEP2TwoFermions myTwoFermions;
    const LEP2oblique myLEP2oblique;
    const LEP2test myTEST;
    const double sqrt_s, s;

    bool flag[NUMofLEP2RCs]; 
    bool bSigmaForAFB;
    bool bSigmaForR;
    
    double Mw, GammaZ;
    
    // caches for the SM prediction
    mutable double SMparams_cache[StandardModel::NumSMParamsForEWPO+3];
    mutable double SMresult_cache; 
    mutable bool flag_cache[NUMofLEP2RCs];
    mutable double ml_cache, mq_cache, mqForHad_cache[6];
    
    
    void SetObParam(LEP2oblique::Oblique ob, double ObParam_i[]) const 
    {
        for (int i=0; i<7; i++) {
            if (i==ob) ObParam_i[ob] = 1.0;
            else ObParam_i[i] = 0.0;
        }
    }
    
    double m_l(const StandardModel::lepton l) const 
    {
        return SM.getLeptons(l).getMass();
    }

    double m_q(const QCD::quark q, const double mu, const orders order=FULLNLO) const 
    {
        switch(q) {
            case QCD::UP:
            case QCD::DOWN:
            case QCD::STRANGE:
                return SM.Mrun(mu, SM.getQuarks(q).getMass_scale(), 
                                         SM.getQuarks(q).getMass(), order);
            case QCD::CHARM:
            case QCD::BOTTOM:
                return SM.Mrun(mu, SM.getQuarks(q).getMass(), order);
            case QCD::TOP:
                return SM.getMtpole(); // the pole mass
            default:
                throw std::runtime_error("Error in LEP2ThObservable::mq()"); 
        }
    }
    
    double sigma_NoISR_l() 
    {
        double sigma = myTwoFermions.sigma_l(l_flavor, ml_cache, s, Mw, GammaZ, flag[Weak]);

        if (!bSigmaForAFB && flag[QEDFSR])
            sigma *= myTwoFermions.QED_FSR_forSigma(s, SM.getLeptons(l_flavor).getCharge());

        return sigma;
    }   
    
    double sigma_NoISR_q() 
    {
        double sigma = myTwoFermions.sigma_q(q_flavor, mq_cache, s, Mw, GammaZ, flag[Weak]);
    
        if (!bSigmaForAFB && flag[QEDFSR])
            sigma *= myTwoFermions.QED_FSR_forSigma(s, SM.getQuarks(q_flavor).getCharge());
        
        if (!bSigmaForAFB && flag[QCDFSR])
            sigma *= myTwoFermions.QCD_FSR_forSigma(s);

        return sigma;
    }      
    
    double AFB_NoISR_l() 
    {
        double AFB = myTwoFermions.AFB_l(l_flavor, ml_cache, s, Mw, GammaZ, flag[Weak]);

        return AFB;
    }
    
    double AFB_NoISR_q() 
    {
        double AFB = myTwoFermions.AFB_q(q_flavor, mq_cache, s, Mw, GammaZ, flag[Weak]);
        
        if (flag[QCDFSR])
            AFB *= myTwoFermions.QCD_FSR_forAFB(q_flavor, mq_cache, s);
        
        return AFB;
    }
    
    double Integrand_sigmaWithISR_l(double x) 
    {
        double sprime = (1.0 - x)*s;
        double sigma = myTwoFermions.sigma_l(l_flavor, ml_cache, sprime, Mw, GammaZ, 
                                             flag[Weak]);
        double H = myTwoFermions.H_ISR(x, s);

        if (!bSigmaForAFB && flag[QEDFSR])
            sigma *= myTwoFermions.QED_FSR_forSigma(sprime, SM.getLeptons(l_flavor).getCharge());

        return ( H*sigma );
    }    
    
    double Integrand_sigmaWithISR_q(double x) 
    {
        double sprime = (1.0 - x)*s;
        double sigma = myTwoFermions.sigma_q(q_flavor, mq_cache, sprime, Mw, GammaZ, 
                                             flag[Weak]);
        double H = myTwoFermions.H_ISR(x, s);
    
        if (!bSigmaForAFB && flag[QEDFSR])
            sigma *= myTwoFermions.QED_FSR_forSigma(sprime, 
                                                    SM.getQuarks(q_flavor).getCharge());
        
        if (!bSigmaForAFB && flag[QCDFSR])
            sigma *= myTwoFermions.QCD_FSR_forSigma(sprime);

        return ( H*sigma );
    }    
    
    double Integrand_dsigmaBox_l(double cosTheta) 
    {
        return ( myTwoFermions.dsigma_l_box(l_flavor, ml_cache, s, cosTheta, Mw, GammaZ) );
    }       
    
    double Integrand_dsigmaBox_q(double cosTheta) 
    {
        return ( myTwoFermions.dsigma_q_box(q_flavor, mq_cache, s, cosTheta, Mw, GammaZ) );
    }       
    
    double Integrand_AFBnumeratorWithISR_l(double x)
    {
        double sprime = (1.0 - x)*s;
        double Ncf = 1.0;
        double G3prime = myTwoFermions.G_3prime_l(l_flavor, ml_cache, sprime, Mw, GammaZ, 
                                                  flag[Weak]);
        double H = myTwoFermions.H_ISR_FB(x, s);

        return ( M_PI*SM.getAle()*SM.getAle()*Ncf*H*G3prime/sprime );
    }
    
    double Integrand_AFBnumeratorWithISR_q(double x) 
    {
        double sprime = (1.0 - x)*s;
        double Ncf = 3.0;
        double G3prime = myTwoFermions.G_3prime_q(q_flavor, mq_cache, sprime, Mw, GammaZ, 
                                                  flag[Weak]);
        double H = myTwoFermions.H_ISR_FB(x, s);
        
        if (flag[QCDFSR])
            G3prime *= myTwoFermions.QCD_FSR_forAFB(q_flavor, mq_cache, sprime);
        
        return ( M_PI*SM.getAle()*SM.getAle()*Ncf*H*G3prime/sprime );
    }
    
    
private:   
      
};

#endif	/* LEP2THOBSERVABLE_H */

