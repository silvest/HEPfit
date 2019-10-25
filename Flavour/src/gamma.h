/* 
 * Copyright (C) 2012 HEPfit Collaboration
 *
 *
 * For the licensing terms see doc/COPYING.
 */

#ifndef GAMMA_H
#define	GAMMA_H

#include "ThObservable.h"


/**
* @class CKMGamma
* @ingroup Flavour
* @brief A class for the CKM angle @f$\gamma @f$ in degrees. 
* @author HEPfit Collaboration
* @copyright GNU General Public License
* @details This class is used to compute the CKM angle @f$\gamma @f$ in degrees.
*/
class CKMGamma : public ThObservable {
public:
    /**
    * @brief Constructor.
    * @param[in] SM_i a reference to an object of type StandardModel
    */
    CKMGamma(const StandardModel& SM_i);

    /**
     * @return The CKM angle @f$\gamma @f$ in degrees
     */
    double computeThValue();
};

#endif	/* GAMMA_H */
