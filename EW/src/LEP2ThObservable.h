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

    // Radiative Corrections for the LEP-II observables
    enum LEP2RCs {Weak=0, WeakBox, ISR, QEDFSR, QCDFSR, NUMofLEP2RCs};
        
    /**
     * @brief LEP2ThObservable constructor
     * @param[in] EW_i an object of EW class
     * @param[in] sqrt_s_i the CM energy of the e^+ e^- pair
     * @param[in] bSigmaForAFB_i true for the denominator of A_FB
     * @param[in] bSigmaForR_i true for the denominator of R_b or R_c
     */
    LEP2ThObservable(const EW& EW_i, const double sqrt_s_i, 
                     const bool bSigmaForAFB_i=false,
                     const bool bSigmaForR_i=false) : ThObservable(EW_i), 
            myEW(EW_i), myTwoFermions(EW_i.getSM()), myLEP2oblique(EW_i), 
            sqrt_s(sqrt_s_i), s(sqrt_s_i*sqrt_s_i), bSigmaForAFB(bSigmaForAFB_i),
            bSigmaForR(bSigmaForR_i){
        flag[Weak] = true;
        flag[WeakBox] = true;
        flag[ISR] = true;
        flag[QEDFSR] = true;
        flag[QCDFSR] = true;
    }

    /**
     * @brief set a flag to control radiative corrections
     * @param[in] str "Weak", "WeakBox", "ISR", "QEDFSR" or "QCDFSR"
     * @param[in] flag_i boolean variable
     */
    void setFlag(const std::string str, const bool flag_i) {
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

    void setFlags(const bool flag_i[]) {
        setFlag("Weak", flag_i[Weak]);
        setFlag("WeakBox", flag_i[WeakBox]);
        setFlag("ISR", flag_i[ISR]);
        setFlag("QEDFSR", flag_i[QEDFSR]);
        setFlag("QCDFSR", flag_i[QCDFSR]);        
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
        for (int i=0; i<NUMofLEP2RCs; i++) {
            if (flag_cache[i] != flag[i]) { 
                flag_cache[i] = flag[i];
                bCache &= false;
            }    
        }

        return bCache;
    }
    
    
protected:
    // These variables have to be initialized in child classes. 
    StandardModel::lepton l_flavor;
    StandardModel::quark q_flavor;
    
    const EW& myEW;
    const LEP2TwoFermions myTwoFermions;
    const LEP2oblique myLEP2oblique;
    const double sqrt_s, s;

    bool flag[NUMofLEP2RCs]; 
    bool bSigmaForAFB;
    bool bSigmaForR;
    
    // caches for the SM prediction
    mutable double SMparams_cache[EWSM::NumSMParams+3];
    mutable double SMresult_cache; 
    mutable bool flag_cache[NUMofLEP2RCs];


    double sigma_NoISR_l() {
        double Mw = SM.Mw(); 
        double GammaZ = myEW.Gamma_Z();
        double sigma = myTwoFermions.sigma_l(l_flavor, s, Mw, GammaZ, flag[Weak]);

        if (!bSigmaForAFB && flag[QEDFSR])
            sigma *= myTwoFermions.QED_FSR_forSigma(s, SM.getLeptons(l_flavor).getCharge());

        return sigma;
    }   
    
    double sigma_NoISR_q() {
        double Mw = SM.Mw(); 
        double GammaZ = myEW.Gamma_Z();
        double sigma = myTwoFermions.sigma_q(q_flavor, s, Mw, GammaZ, flag[Weak]);
    
        if (!bSigmaForAFB && flag[QEDFSR])
            sigma *= myTwoFermions.QED_FSR_forSigma(s, SM.getQuarks(q_flavor).getCharge());
        
        if (!bSigmaForAFB && flag[QCDFSR])
            sigma *= myTwoFermions.QCD_FSR_forSigma(s);

        return sigma;
    }      
    
    double AFB_NoISR_l() {
        double Mw = SM.Mw(); 
        double GammaZ = myEW.Gamma_Z();
        double AFB = myTwoFermions.AFB_l(l_flavor, s, Mw, GammaZ, flag[Weak]);

        return AFB;
    }
    
    double AFB_NoISR_q() {
        double Mw = SM.Mw(); 
        double GammaZ = myEW.Gamma_Z();
        double AFB = myTwoFermions.AFB_q(q_flavor, s, Mw, GammaZ, flag[Weak]);
        
        if (flag[QCDFSR])
            AFB *= myTwoFermions.QCD_FSR_forAFB(q_flavor, s);
        
        return AFB;
    }
    
    double Integrand_sigmaWithISR_l(double x) {
        double sprime = (1.0 - x)*s;
        double Mw = SM.Mw(); 
        double GammaZ = myEW.Gamma_Z();
        double sigma = myTwoFermions.sigma_l(l_flavor, sprime, Mw, GammaZ, 
                                             flag[Weak]);
        double H = myTwoFermions.H_ISR(x, s);

        if (!bSigmaForAFB && flag[QEDFSR])
            sigma *= myTwoFermions.QED_FSR_forSigma(sprime, SM.getLeptons(l_flavor).getCharge());

        return ( H*sigma );
    }    
    
    double Integrand_sigmaWithISR_q(double x) {
        double sprime = (1.0 - x)*s;
        double Mw = SM.Mw(); 
        double GammaZ = myEW.Gamma_Z();
        double sigma = myTwoFermions.sigma_q(q_flavor, sprime, Mw, GammaZ, 
                                             flag[Weak]);
        double H = myTwoFermions.H_ISR(x, s);
    
        if (!bSigmaForAFB && flag[QEDFSR])
            sigma *= myTwoFermions.QED_FSR_forSigma(sprime, 
                                                    SM.getQuarks(q_flavor).getCharge());
        
        if (!bSigmaForAFB && flag[QCDFSR])
            sigma *= myTwoFermions.QCD_FSR_forSigma(sprime);

        return ( H*sigma );
    }    
    
    double Integrand_dsigmaBox_l(double cosTheta) {
        double Mw = SM.Mw(); 
        double GammaZ = myEW.Gamma_Z();
        return ( myTwoFermions.dsigma_l_box(l_flavor, s, cosTheta, Mw, GammaZ) );
    }       
    
    double Integrand_dsigmaBox_q(double cosTheta) {
        double Mw = SM.Mw(); 
        double GammaZ = myEW.Gamma_Z();
        return ( myTwoFermions.dsigma_q_box(q_flavor, s, cosTheta, Mw, GammaZ) );
    }       
    
    double Integrand_AFBnumeratorWithISR_l(double x) {
        double sprime = (1.0 - x)*s;
        double Mw = SM.Mw(); 
        double GammaZ = myEW.Gamma_Z();
        double Ncf = 1.0;
        double G3prime = myTwoFermions.G_3prime_l(l_flavor, sprime, Mw, GammaZ, 
                                                  flag[Weak]);
        double H = myTwoFermions.H_ISR_FB(x, s);

        return ( M_PI*SM.getAle()*SM.getAle()*Ncf*H*G3prime/sprime );
    }
    
    double Integrand_AFBnumeratorWithISR_q(double x) {
        double sprime = (1.0 - x)*s;
        double Mw = SM.Mw(); 
        double GammaZ = myEW.Gamma_Z();
        double Ncf = 3.0;
        double G3prime = myTwoFermions.G_3prime_q(q_flavor, sprime, Mw, GammaZ, 
                                                  flag[Weak]);
        double H = myTwoFermions.H_ISR_FB(x, s);
        
        if (flag[QCDFSR])
            G3prime *= myTwoFermions.QCD_FSR_forAFB(q_flavor, sprime);
        
        return ( M_PI*SM.getAle()*SM.getAle()*Ncf*H*G3prime/sprime );
    }
    
    
private:   
      
};

#endif	/* LEP2THOBSERVABLE_H */

