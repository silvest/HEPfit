/* 
 * Copyright (C) 2012 SUSYfit Collaboration
 * All rights reserved.
 *
 * For the licensing terms see doc/COPYING.
 */

#ifndef LEP2SIGMAMU_H
#define	LEP2SIGMAMU_H

#include "LEP2ThObservable.h"


class LEP2sigmaMu : public LEP2ThObservable {
public:

    /**
     * @brief LEP2sigmaMu constructor
     * @param[in] EW_i an object of EW class
     * @param[in] sqrt_s_i the CM energy of the e^+ e^- pair
     * @param[in] bSigmaForAFB_i true for the denominator of A_FB
     */
    LEP2sigmaMu(const EW& EW_i, const double sqrt_s_i, 
                const bool bSigmaForAFB_i=false) : LEP2ThObservable(EW_i, sqrt_s_i, bSigmaForAFB_i) {
        l_flavor = StandardModel::MU;
    }

    /**
     * @return the cross section for e^+ e^- -> mu^+ mu^- at sqrt_s in pb
     */
    double getThValue();

private:
    
};

#endif	/* LEP2SIGMAMU_H */

