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
#include "LEP2oblique.h"


class LEP2ThObservable : public ThObservable  {
public:
    
    /**
     * @brief LEP2ThObservable constructor
     * @param[in] EW_i an object of EW class
     * @param[in] sqrt_s_i the CM energy of the e^+ e^- pair
     */
    LEP2ThObservable(const EW& EW_i, const double sqrt_s_i) : ThObservable(EW_i), 
            myEW(EW_i), myLEP2oblique(EW_i), sqrt_s(sqrt_s_i) {
        bRCs[0] = true;
        bRCs[1] = true;
        bRCs[2] = true; 
        bRCs[3] = true; 
        bRCs[4] = true; 
    }

    /**
     * @brief set a flag to control radiative corrections
     * @param[in] str "Weak", "WeakBox", "ISR", "QEDFSR" or "QCDFSR"
     * @param[in] flag boolean variable
     */
    void setFlag(const std::string str, const bool flag) {
        if (str=="Weak")
            bRCs[0] = flag;
        else if (str=="WeakBox")
            bRCs[1] = flag;
        else if (str=="ISR")
            bRCs[2] = flag;
        else if (str=="QEDFSR")
            bRCs[3] = flag;
        else if (str=="QCDFSR")
            bRCs[4] = flag;
        else
            throw std::runtime_error("Error in LEP2ThObservable::setFlag()");             
    }
    
protected:
    bool bRCs[5]; // flags for radiative corrections

    const EW& myEW;
    const LEP2oblique myLEP2oblique;
    const double sqrt_s;
    
    // caches for the SM prediction
    mutable double SMparams_cache[EWSM::NumSMParams+3];
    mutable double SMresult_cache; 
    mutable bool bRCs_cache[5];

private:

};

#endif	/* LEP2THOBSERVABLE_H */

