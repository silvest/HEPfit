/* 
 * Copyright (C) 2012 HEPfit Collaboration
 *
 *
 * For the licensing terms see doc/COPYING.
 */

#ifndef LEP2AFBCHARM_H
#define	LEP2AFBCHARM_H

#include "LEP2ThObservable.h"

/**
 * @class LEP2AFBcharm
 * @ingroup EW
 * @brief A class for the forward-backward asymmetry of @f$e^+e^-\to c\bar{c}@f$ 
 * above the @f$Z@f$ pole.
 * @author HEPfit Collaboration
 * @copyright GNU General Public License
 * @details 
 */
class LEP2AFBcharm : public LEP2ThObservable {
public:

    /**
     * @brief LEP2AFBcharm constructor
     * @param[in] SM_i a reference to an object of type StandardModel
     * @param[in] sqrt_s_i the CM energy of the e^+ e^- pair
     */
    LEP2AFBcharm(const StandardModel& SM_i, const double sqrt_s_i) 
    : LEP2ThObservable(SM_i, sqrt_s_i)
    {
    }

    /**
     * @return the forward-backward asymmetry for e^+ e^- -> c cbar at sqrt_s
     */
    double computeThValue();

private:
    
};

#endif	/* LEP2AFBCHARM_H */

