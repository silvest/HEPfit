/* 
 * Copyright (C) 2012 SUSYfit Collaboration
 * All rights reserved.
 *
 * For the licensing terms see doc/COPYING.
 */

#ifndef LEP2SIGMATAU_H
#define	LEP2SIGMATAU_H

#include <ThObservable.h>
#include "EWSM.h"
#include "EW.h"
#include "LEP2oblique.h"


class LEP2sigmaTau : public ThObservable {
public:

    /**
     * @brief LEP2sigmaTau constructor
     * @param[in] EW_i an object of EW class
     * @param[in] sqrt_s_i the CM energy of the e^+ e^- pair
     */
    LEP2sigmaTau(const EW& EW_i, const double sqrt_s_i);

    /**
     * @return the cross section for e^+ e^- -> tau^+ tau^- at sqrt_s in pb
     */
    double getThValue();

    /**
     * @brief set flags for radiative corrections
     * @param[in] bDP_i with/without dressed gauge-boson propagators
     * @param[in] bWEAK_i with/without weak corrections
     * @param[in] bQED_i with/without QED corrections
     */
    void setFlags(const bool bDP_i, const bool bWEAK_i, const bool bQED_i) {
        this->bDP = bDP_i;
        this->bWEAK = bWEAK_i;
        this->bQED = bQED_i;
    }
    

private:
    const EW& myEW;
    const LEP2oblique myLEP2oblique;
    const double sqrt_s;
    bool bDP, bWEAK, bQED;        
    
    // caches for the SM prediction
    mutable double SMparams_cache[EWSM::NumSMParams+3];
    mutable bool   bool_cache[3];
    mutable double SMresult_cache;
};

#endif	/* LEP2SIGMATAU_H */

