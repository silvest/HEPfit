/* 
 * File:   LEP2sigmaMu.h
 * Author: mishima
 */

#ifndef LEP2SIGMAMU_H
#define	LEP2SIGMAMU_H

#include <ThObservable.h>
#include "EW.h"


class LEP2sigmaMu : public ThObservable {
public:

    /**
     * @brief LEP2sigmaMu constructor
     * @param[in] EW_i an object of EW class
     * @param[in] sqrt_s_i the CM energy of the e^+ e^- pair
     */
    LEP2sigmaMu(const EW& EW_i, const double sqrt_s_i);

    /**
     * @return the cross section for e^+ e^- -> mu^+ mu^- at sqrt_s in pb
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
    double sqrt_s;
    bool bDP, bWEAK, bQED;
};

#endif	/* LEP2SIGMAMU_H */

