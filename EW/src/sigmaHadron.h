/* 
 * Copyright (C) 2012-2013 SusyFit Collaboration
 * All rights reserved.
 *
 * For the licensing terms see doc/COPYING.
 */

#ifndef SIGMAHADRON_H
#define	SIGMAHADRON_H

#include <stdexcept>
#include <ThObservable.h>
#include "EW.h"

/**
 * @class sigmaHadron 
 * @ingroup EW 
 * @brief A class for the hadronic cross section of the @f$Z@f$ boson at the 
 * @f$Z@f$ pole
 * @author SusyFit Collaboration
 * @copyright GNU General Public License
 * @details This class is used to compute the hadronic cross section of the 
 * @f$Z@f$ boson at the @f$Z@f$ pole. 
 */
class sigmaHadron : public ThObservable {
public:

    /**
     * A constructor.
     * @param[in] EW_i A reference to an object of EW class, which is the base 
     * class of the electroweak precision observables.
     */
    sigmaHadron(const EW& EW_i) 
    : ThObservable(EW_i), myEW(EW_i), myEWTYPE(EW_i.getEWTYPE()) 
    {
    };

    /**
     * @return The hadronic cross section of the @f$Z@f$ boson at the @f$Z@f$ pole in nb. 
     */
    double getThValue();


private:

    /**
     * A reference to an object of EW class, which is the base class of the 
     * electroweak precision observables.
     */
    const EW& myEW;

    /**
     * An enumerator controlling the formulae used in the computation.
     */
    const EW::EWTYPE myEWTYPE;
};

#endif	/* SIGMAHADRON_H */

