/* 
 * Copyright (C) 2012 HEPfit Collaboration
 *
 *
 * For the licensing terms see doc/COPYING.
 */

#ifndef ALPHA_2A_H
#define	ALPHA_2A_H

#include "ThObservable.h"

/**
* @class Alpha_2a
* @ingroup Flavour
* @brief A class for the CKM angle @f$\alpha @f$ in degrees. 
* @author HEPfit Collaboration
* @copyright GNU General Public License
* @details This class is used to compute the CKM angle @f$\alpha @f$ in degrees,
* extracted as @f$\pi + 1/2 arg(%AmpDB2) - \gamma @f$. The results is given in the
 * range [0, 180] deg.
*/
class Alpha_2a : public ThObservable, AmpDB2 {
public:
    /**
    * @brief Constructor.
    * @param[in] SM_i a reference to an object of type StandardModel
    */
    Alpha_2a(const StandardModel& SM_i) : ThObservable(SM_i), AmpDB2(SM_i) {};

    /**
     * @return The CKM angle @f$\alpha @f$ in degrees
     */
    double computeThValue();
};

#endif	/* ALPHA_2A_H */

