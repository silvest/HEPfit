/* 
 * Copyright (C) 2012-2013 SusyFit Collaboration
 * All rights reserved.
 *
 * For the licensing terms see doc/COPYING.
 */

#ifndef RCHARM_H
#define	RCHARM_H

#include <stdexcept>
#include <ThObservable.h>
#include "EW.h"

/**
 * @class Rcharm
 * @ingroup EW 
 * @brief A class for @f$R_c^0=\Gamma(Z\to c\bar{c})/\Gamma(Z\to {\rm hadrons})@f$ 
 * at the @f$Z@f$ pole.
 * @author SusyFit Collaboration
 * @copyright GNU General Public License
 * @details This class is used to compute the ratio of the @f$Z\to c\bar{c}@f$ 
 * width to the @f$Z@f$-boson hadronic width at the @f$Z@f$ pole. 
 */
class Rcharm : public ThObservable {
public:

    /**
     * @brief Constructor.
     * @param[in] EW_i A reference to an object of EW class, which is the base 
     * class of the electroweak precision observables.
     */
    Rcharm(const EW& EW_i) 
    : ThObservable(EW_i), myEW(EW_i)
    {
    };

    /**
     * @brief The theory prediction for \f$R_c^0\f$ at the \f$Z\f$ pole.
     * @return The ratio of the @f$Z\to c\bar{c}@f$ width to the @f$Z@f$-boson 
     * hadronic width at the @f$Z@f$ pole. 
     */
    double computeThValue();

    
private:

    /**
     * A reference to an object of EW class, which is the base class of the 
     * electroweak precision observables.
     */
    const EW& myEW;///< A reference to an object of the EW class.
};

#endif	/* RCHARM_H */

