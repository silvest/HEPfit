/* 
 * Copyright (C) 2012 HEPfit Collaboration
 *
 *
 * For the licensing terms see doc/COPYING.
 */

#ifndef LEP2RCHARM_H
#define	LEP2RCHARM_H

#include "LEP2ThObservable.h"

/**
 * @class LEP2Rcharm
 * @ingroup EW
 * @brief A class for @f$R_c^0@f$ above the @f$Z@f$ pole.
 * @author HEPfit Collaboration
 * @copyright GNU General Public License
 * @details 
 */
class LEP2Rcharm : public LEP2ThObservable {
public:

    /**
     * @brief LEP2Rcharm constructor
     * @param[in] SM_i a reference to an object of type StandardModel
     * @param[in] sqrt_s_i the CM energy of the e^+ e^- pair
     */
    LEP2Rcharm(const StandardModel& SM_i, const double sqrt_s_i) 
    : LEP2ThObservable(SM_i, sqrt_s_i)
    {

    }

    /**
     * @return the ratio of the c-cbar cross section to the hadronic cross section at sqrt_s
     */
    double computeThValue();

private:
    
};

#endif	/* LEP2RCHARM_H */

