/* 
 * Copyright (C) 2012-2013 SusyFit Collaboration
 * All rights reserved.
 *
 * For the licensing terms see doc/COPYING.
 */

#ifndef LEP2AFBTAU_H
#define	LEP2AFBTAU_H

#include "LEP2ThObservable.h"
#include "LEP2sigmaTau.h"

/**
 * @class LEP2AFBtau
 * @ingroup EW
 * @brief A class for the forward-backward asymmetry of @f$e^+e^-\to \mu^+\mu^-@f$ 
 * above the @f$Z@f$ pole.
 * @author SusyFit Collaboration
 * @copyright GNU General Public License
 * @details 
 */
class LEP2AFBtau : public LEP2ThObservable {
public:

    /**
     * @brief LEP2AFBtau constructor
     * @param[in] EW_i an object of EW class
     * @param[in] sqrt_s_i the CM energy of the e^+ e^- pair
     */
    LEP2AFBtau(const EW& EW_i, const double sqrt_s_i) 
    : LEP2ThObservable(EW_i, sqrt_s_i), myLEP2sigmaTau(EW_i, sqrt_s_i, true)
    {
        l_flavor = StandardModel::TAU;
    }

    /**
     * @return the forward-backward asymmetry for e^+ e^- -> tau^+ tau^- at sqrt_s
     */
    double computeThValue();

private:
    LEP2sigmaTau myLEP2sigmaTau;
    
};

#endif	/* LEP2AFBTAU_H */

