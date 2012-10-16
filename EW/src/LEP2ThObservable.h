/* 
 * File:   LEP2ThObservable.h
 * Author: mishima
 */

#ifndef LEP2THOBSERVABLE_H
#define	LEP2THOBSERVABLE_H

#include <map>
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
    }

    /**
     * @brief set flags for radiative corrections
     */
    void setFlags(const std::string name, const bool& value) {
        if (name.compare("Weak") == 0)
            this->Flags["Weak"] = value;
    }
    
protected:
    const EW& myEW;
    const LEP2oblique myLEP2oblique;
    const double sqrt_s;
    std::map<std::string, bool> Flags;
    
    // caches for the SM prediction
    mutable double SMparams_cache[EWSM::NumSMParams+3];
    mutable bool   bool_cache[3];
    mutable double SMresult_cache;    

private:

};

#endif	/* LEP2THOBSERVABLE_H */

