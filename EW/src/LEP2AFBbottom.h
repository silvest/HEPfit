/* 
 * Copyright (C) 2012 HEPfit Collaboration
 * All rights reserved.
 *
 * For the licensing terms see doc/COPYING.
 */

#ifndef LEP2AFBBOTTOM_H
#define	LEP2AFBBOTTOM_H

#include "LEP2ThObservable.h"
#include "LEP2sigmaBottom.h"

/**
 * @class LEP2AFBbottom
 * @ingroup EW
 * @brief A class for the forward-backward asymmetry of @f$e^+e^-\to b\bar{b}@f$ 
 * above the @f$Z@f$ pole.
 * @author HEPfit Collaboration
 * @copyright GNU General Public License
 * @details 
 */
class LEP2AFBbottom : public LEP2ThObservable {
public:

    /**
     * @brief LEP2AFBbottom constructor
     * @param[in] SM_i a reference to an object of type StandardModel
     * @param[in] sqrt_s_i the CM energy of the e^+ e^- pair
     */
    LEP2AFBbottom(const StandardModel& SM_i, const double sqrt_s_i) 
    : LEP2ThObservable(SM_i, sqrt_s_i), myLEP2sigmaBottom(SM_i, sqrt_s_i, true) 
    {
        q_flavor = QCD::BOTTOM;
    }

    /**
     * @return the forward-backward asymmetry for e^+ e^- -> b bbar at sqrt_s
     */
    double computeThValue();

private:
    LEP2sigmaBottom myLEP2sigmaBottom;
     
};

#endif	/* LEP2AFBBOTTOM_H */

