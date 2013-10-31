/* 
 * Copyright (C) 2012-2013 SusyFit Collaboration
 * All rights reserved.
 *
 * For the licensing terms see doc/COPYING.
 */

#ifndef LEP2SIGMAMU_H
#define	LEP2SIGMAMU_H

#include "LEP2ThObservable.h"

/**
 * @class LEP2sigmaMu
 * @ingroup EW
 * @brief A class for the cross section of @f$e^+e^-\to \mu^+\mu^-@f$ above the @f$Z@f$ pole.
 * @author SusyFit Collaboration
 * @copyright GNU General Public License
 * @details 
 */
class LEP2sigmaMu : public LEP2ThObservable {
public:

    /**
     * @brief LEP2sigmaMu constructor
     * @param[in] EW_i an object of EW class
     * @param[in] sqrt_s_i the CM energy of the e^+ e^- pair
     * @param[in] bSigmaForAFB_i true for the denominator of A_FB
     */
    LEP2sigmaMu(const EW& EW_i, const double sqrt_s_i, 
                const bool bSigmaForAFB_i=false) 
    : LEP2ThObservable(EW_i, sqrt_s_i, bSigmaForAFB_i) 
    {
        l_flavor = StandardModel::MU;
    }

    /**
     * @return the cross section for e^+ e^- -> mu^+ mu^- at sqrt_s in pb
     */
    double computeThValue();

private:
    
};

#endif	/* LEP2SIGMAMU_H */

