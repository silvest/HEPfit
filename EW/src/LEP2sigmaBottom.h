/* 
 * Copyright (C) 2012 HEPfit Collaboration
 * All rights reserved.
 *
 * For the licensing terms see doc/COPYING.
 */

#ifndef LEP2SIGMABOTTOM_H
#define	LEP2SIGMABOTTOM_H

#include "LEP2ThObservable.h"

/**
 * @class LEP2sigmaBottom
 * @ingroup EW
 * @brief A class for the cross section of @f$e^+e^-\to b\bar{b}@f$ above the @f$Z@f$ pole.
 * @author HEPfit Collaboration
 * @copyright GNU General Public License
 * @details 
 */
class LEP2sigmaBottom : public LEP2ThObservable {
public:

    /**
     * @brief LEP2sigmaBottom constructor
     * @param[in] SM_i a reference to an object of type StandardModel
     * @param[in] sqrt_s_i the CM energy of the e^+ e^- pair
     * @param[in] bSigmaForAFB_i true for the denominator of A_FB
     */
    LEP2sigmaBottom(const StandardModel& SM_i, const double sqrt_s_i, 
                    const bool bSigmaForAFB_i=false) 
    : LEP2ThObservable(SM_i, sqrt_s_i, bSigmaForAFB_i)
    {
        q_flavor = QCD::BOTTOM;
    }

    /**
     * @return the cross section for e^+ e^- -> b bbar at sqrt_s in pb
     */
    double computeThValue();

private:
    
};

#endif	/* LEP2SIGMABOTTOM_H */

