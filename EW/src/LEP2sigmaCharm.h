/* 
 * Copyright (C) 2012 HEPfit Collaboration
 *
 *
 * For the licensing terms see doc/COPYING.
 */

#ifndef LEP2SIGMACHARM_H
#define	LEP2SIGMACHARM_H

#include "LEP2ThObservable.h"

/**
 * @class LEP2sigmaCharm
 * @ingroup EW
 * @brief A class for the cross section of @f$e^+e^-\to c\bar{c}@f$ above the @f$Z@f$ pole.
 * @author HEPfit Collaboration
 * @copyright GNU General Public License
 * @details 
 */
class LEP2sigmaCharm : public LEP2ThObservable {
public:

    /**
     * @brief LEP2sigmaCharm constructor
     * @param[in] SM_i a reference to an object of type StandardModel
     * @param[in] sqrt_s_i the CM energy of the e^+ e^- pair
     * @param[in] bSigmaForAFB_i true for the denominator of A_FB
     */
    LEP2sigmaCharm(const StandardModel& SM_i, const double sqrt_s_i) 
    : LEP2ThObservable(SM_i, sqrt_s_i) 
    {

    }

    /**
     * @return the cross section for e^+ e^- -> c cbar at sqrt_s in pb
     */
    double computeThValue();

private:
    
};

#endif	/* LEP2SIGMACHARM_H */

