/* 
 * Copyright (C) 2012 HEPfit Collaboration
 *
 *
 * For the licensing terms see doc/COPYING.
 */

#ifndef ALPHA_H
#define	ALPHA_H

#include "ThObservable.h"
#include "AmpDB2.h"

/**
* @class Alpha
* @ingroup Flavour
* @brief A class for the CKM angle @f$\alpha @f$ in degrees. 
* @author HEPfit Collaboration
* @copyright GNU General Public License
* @details This class is used to compute the CKM angle @f$\alpha @f$ in degrees,
* extracted as @f$\pi + 1/2 arg(%AmpDB2) - \gamma @f$.
*/
class Alpha : public ThObservable, AmpDB2 {
public:
    /**
    * @brief Constructor.
    * @param[in] SM_i a reference to an object of type StandardModel
    */
    Alpha(const StandardModel& SM_i);

    /**
     * @return The CKM angle @f$\alpha @f$ in degrees
     */
    double computeThValue();
};

#endif	/* ALPHA_H */
