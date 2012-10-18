/* 
 * File:   LEP2ThObservable.h
 * Author: mishima
 */

#ifndef LEP2THOBSERVABLE_H
#define	LEP2THOBSERVABLE_H

#include <stdexcept>
#include <cstring>
#include <ThObservable.h>
#include "EWSM.h"
#include "EW.h"
#include "LEP2TwoFermions.h"
#include "LEP2oblique.h"


class LEP2ThObservable : public ThObservable  {
public:
    
    /**
     * @brief LEP2ThObservable constructor
     * @param[in] EW_i an object of EW class
     * @param[in] sqrt_s_i the CM energy of the e^+ e^- pair
     */
    LEP2ThObservable(const EW& EW_i, const double sqrt_s_i) : ThObservable(EW_i), 
            myEW(EW_i), myTwoFermions(EW_i.getSM()), myLEP2oblique(EW_i), 
            sqrt_s(sqrt_s_i), s(sqrt_s_i*sqrt_s_i) {
        bRCs[LEP2TwoFermions::Weak] = true;
        bRCs[LEP2TwoFermions::WeakBox] = true;
        bRCs[LEP2TwoFermions::ISR] = true;
        bRCs[LEP2TwoFermions::QEDFSR] = true;
        bRCs[LEP2TwoFermions::QCDFSR] = true;
    }

    /**
     * @brief set a flag to control radiative corrections
     * @param[in] str "Weak", "WeakBox", "ISR", "QEDFSR" or "QCDFSR"
     * @param[in] flag boolean variable
     */
    void setFlag(const std::string str, const bool flag) {
        if (str=="Weak")
            bRCs[LEP2TwoFermions::Weak] = flag;
        else if (str=="WeakBox")
            bRCs[LEP2TwoFermions::WeakBox] = flag;
        else if (str=="ISR")
            bRCs[LEP2TwoFermions::ISR] = flag;
        else if (str=="QEDFSR")
            bRCs[LEP2TwoFermions::QEDFSR] = flag;
        else if (str=="QCDFSR")
            bRCs[LEP2TwoFermions::QCDFSR] = flag;
        else
            throw std::runtime_error("Error in LEP2ThObservable::setFlag()");
    }

    
    bool checkSMparams(const double s, const double Mw, const double GammaZ) const {
        // 21 SM parameters in checkSMparams() + s, Mw, GammaZ + 5 booleans
        bool bCache = true;
        bCache &= myEW.getSM().getEWSM()->checkSMparams(SMparams_cache);
        
        if (SMparams_cache[EWSM::NumSMParams] != s) { 
            SMparams_cache[EWSM::NumSMParams] = s;
            bCache &= false;
        }    
        if (SMparams_cache[EWSM::NumSMParams+1] != Mw) { 
            SMparams_cache[EWSM::NumSMParams+1] = Mw;
            bCache &= false;
        }    
        if (SMparams_cache[EWSM::NumSMParams+2] != GammaZ) { 
            SMparams_cache[EWSM::NumSMParams+2] = GammaZ;
            bCache &= false;
        }    
        for (int i=0; i<LEP2TwoFermions::NUMofLEP2RCs; i++) {
            if (bRCs_cache[i] != bRCs[i]) { 
                bRCs_cache[i] = bRCs[i];
                bCache &= false;
            }    
        }

        return bCache;
    }
    
protected:
    bool bRCs[LEP2TwoFermions::NUMofLEP2RCs]; // flags for radiative corrections

    // used in integrands
    StandardModel::lepton l_flavor;
    StandardModel::quark q_flavor;
    
    const EW& myEW;
    const LEP2TwoFermions myTwoFermions;
    const LEP2oblique myLEP2oblique;
    const double sqrt_s, s;
    
    // caches for the SM prediction
    mutable double SMparams_cache[EWSM::NumSMParams+3];
    mutable double SMresult_cache; 
    mutable bool bRCs_cache[LEP2TwoFermions::NUMofLEP2RCs];

protected:
    
    double IntegrandISR_AFB_l(double x) {
        double sprime = (1.0 - x)*s;
        double Mw = SM.Mw(); 
        double GammaZ = myEW.Gamma_Z();
        double Ncf = 1.0;
        
        double G3prime = myTwoFermions.G_3prime_l(l_flavor, sprime, Mw, GammaZ, bRCs);
        double H = myTwoFermions.H_ISR_FB(x, s);
        
        return ( Ncf*H*G3prime/sprime );
    }
    
    double IntegrandISR_AFB_q(double x) {
        double sprime = (1.0 - x)*s;
        double Mw = SM.Mw(); 
        double GammaZ = myEW.Gamma_Z();
        double Ncf = 3.0;
        
        double G3prime = myTwoFermions.G_3prime_q(q_flavor, sprime, Mw, GammaZ, bRCs);
        double H = myTwoFermions.H_ISR_FB(x, s);
        
        return ( Ncf*H*G3prime/sprime );
    }

    double IntegrandISR_sigma_l(double x) {
        double sprime = (1.0 - x)*s;
        double Mw = SM.Mw(); 
        double GammaZ = myEW.Gamma_Z();
        
        bool bRCs_noQCDFSR[5];
        bRCs_noQCDFSR[LEP2TwoFermions::Weak] = bRCs[LEP2TwoFermions::Weak];
        bRCs_noQCDFSR[LEP2TwoFermions::WeakBox] = bRCs[LEP2TwoFermions::WeakBox];
        bRCs_noQCDFSR[LEP2TwoFermions::ISR] = bRCs[LEP2TwoFermions::ISR];
        bRCs_noQCDFSR[LEP2TwoFermions::QEDFSR] = bRCs[LEP2TwoFermions::QEDFSR];
        bRCs_noQCDFSR[LEP2TwoFermions::QCDFSR] = false;
        
        double sigma = myTwoFermions.sigma_l(l_flavor, sprime, Mw, GammaZ, bRCs_noQCDFSR);
        double H = myTwoFermions.H_ISR(x, s);
        
        return ( H*sigma );
    }    
    
    double IntegrandISR_sigma_q(double x) {
        double sprime = (1.0 - x)*s;
        double Mw = SM.Mw(); 
        double GammaZ = myEW.Gamma_Z();
        
        bool bRCs_noQCDFSR[5];
        bRCs_noQCDFSR[LEP2TwoFermions::Weak] = bRCs[LEP2TwoFermions::Weak];
        bRCs_noQCDFSR[LEP2TwoFermions::WeakBox] = bRCs[LEP2TwoFermions::WeakBox];
        bRCs_noQCDFSR[LEP2TwoFermions::ISR] = bRCs[LEP2TwoFermions::ISR];
        bRCs_noQCDFSR[LEP2TwoFermions::QEDFSR] = bRCs[LEP2TwoFermions::QEDFSR];
        bRCs_noQCDFSR[LEP2TwoFermions::QCDFSR] = false;
        
        double sigma = myTwoFermions.sigma_q(q_flavor, sprime, Mw, GammaZ, bRCs_noQCDFSR);
        double H = myTwoFermions.H_ISR(x, s);
        
        return ( H*sigma );
    }    
    
private:    
    
};

#endif	/* LEP2THOBSERVABLE_H */

