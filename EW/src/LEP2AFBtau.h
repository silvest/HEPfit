/* 
 * Copyright (C) 2012 HEPfit Collaboration
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
 * @author HEPfit Collaboration
 * @copyright GNU General Public License
 * @details 
 */
class LEP2AFBtau : public LEP2ThObservable {
public:

    /**
     * @brief LEP2AFBtau constructor
     * @param[in] SM_i a reference to an object of type StandardModel
     * @param[in] sqrt_s_i the CM energy of the e^+ e^- pair
     */
    LEP2AFBtau(const StandardModel& SM_i, const double sqrt_s_i) 
    : LEP2ThObservable(SM_i, sqrt_s_i), myLEP2sigmaTau(SM_i, sqrt_s_i, true)
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

