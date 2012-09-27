/* 
 * File:   LEP2AFBtau.h
 * Author: mishima
 */

#ifndef LEP2AFBTAU_H
#define	LEP2AFBTAU_H

#include <ThObservable.h>
#include "EW.h"


class LEP2AFBtau : public ThObservable {
public:

    /**
     * @brief LEP2AFBtau constructor
     * @param[in] EW_i an object of EW class
     * @param[in] sqrt_s_i the CM energy of the e^+ e^- pair
     */
    LEP2AFBtau(const EW& EW_i, const double sqrt_s_i);

    /**
     * @return the forward-backward asymmetry for e^+ e^- -> tau^+ tau^- at sqrt_s
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
    const double sqrt_s;
    bool bDP, bWEAK, bQED;
};

#endif	/* LEP2AFBTAU_H */

