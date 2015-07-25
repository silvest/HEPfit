/* 
 * Copyright (C) 2012 HEPfit Collaboration
 * All rights reserved.
 *
 * For the licensing terms see doc/COPYING.
 */

#ifndef LEP2AFBMU_H
#define	LEP2AFBMU_H

#include "LEP2ThObservable.h"
#include "LEP2sigmaMu.h"

/**
 * @class LEP2AFBmu
 * @ingroup EW
 * @brief A class for the forward-backward asymmetry of @f$e^+e^-\to \mu^+\mu^-@f$ 
 * above the @f$Z@f$ pole.
 * @author HEPfit Collaboration
 * @copyright GNU General Public License
 * @details 
 */
class LEP2AFBmu : public LEP2ThObservable {
public:

    /**
     * @brief LEP2AFBmu constructor
     * @param[in] SM_i a reference to an object of type StandardModel
     * @param[in] sqrt_s_i the CM energy of the e^+ e^- pair
     */
    LEP2AFBmu(const StandardModel& SM_i, const double sqrt_s_i) 
    : LEP2ThObservable(SM_i, sqrt_s_i), myLEP2sigmaMu(SM_i, sqrt_s_i, true) 
    {
        l_flavor = StandardModel::MU;
    }
    
    /**
     * @return the forward-backward asymmetry for e^+ e^- -> mu^+ mu^- at sqrt_s
     */
    double computeThValue();

private:
    LEP2sigmaMu myLEP2sigmaMu;
    
};

#endif	/* LEP2AFBMU_H */

