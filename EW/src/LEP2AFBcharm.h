/* 
 * Copyright (C) 2012 SusyFit Collaboration
 * All rights reserved.
 *
 * For the licensing terms see doc/COPYING.
 */

#ifndef LEP2AFBCHARM_H
#define	LEP2AFBCHARM_H

#include "LEP2ThObservable.h"
#include "LEP2sigmaCharm.h"

/**
 * @class LEP2AFBcharm
 * @ingroup EW
 * @brief A class for the forward-backward asymmetry of @f$e^+e^-\to c\bar{c}@f$ 
 * above the @f$Z@f$ pole.
 * @author SusyFit Collaboration
 * @copyright GNU General Public License
 * @details 
 */
class LEP2AFBcharm : public LEP2ThObservable {
public:

    /**
     * @brief LEP2AFBcharm constructor
     * @param[in] EW_i an object of EW class
     * @param[in] sqrt_s_i the CM energy of the e^+ e^- pair
     */
    LEP2AFBcharm(const StandardModel& SM_i, const double sqrt_s_i) 
    : LEP2ThObservable(SM_i, sqrt_s_i), myLEP2sigmaCharm(SM_i, sqrt_s_i, true) 
    {
        q_flavor = QCD::CHARM;
    }

    /**
     * @return the forward-backward asymmetry for e^+ e^- -> c cbar at sqrt_s
     */
    double computeThValue();

private:
    LEP2sigmaCharm myLEP2sigmaCharm;
    
};

#endif	/* LEP2AFBCHARM_H */

