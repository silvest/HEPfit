/* 
 * Copyright (C) 2012-2013 SusyFit Collaboration
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
 * @author SusyFit Collaboration
 * @copyright GNU General Public License
 * @details 
 */
class LEP2AFBmu : public LEP2ThObservable {
public:

    /**
     * @brief LEP2AFBmu constructor
     * @param[in] EW_i an object of EW class
     * @param[in] sqrt_s_i the CM energy of the e^+ e^- pair
     */
    LEP2AFBmu(const EW& EW_i, const double sqrt_s_i) 
    : LEP2ThObservable(EW_i, sqrt_s_i), myLEP2sigmaMu(EW_i, sqrt_s_i, true) 
    {
        l_flavor = StandardModel::MU;
    }
    
    /**
     * @return the forward-backward asymmetry for e^+ e^- -> mu^+ mu^- at sqrt_s
     */
    double getThValue();

private:
    LEP2sigmaMu myLEP2sigmaMu;
    
};

#endif	/* LEP2AFBMU_H */

