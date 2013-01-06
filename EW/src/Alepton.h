/* 
 * Copyright (C) 2012 SusyFit Collaboration
 * All rights reserved.
 *
 * For the licensing terms see doc/COPYING.
 */

#ifndef ALEPTON_H
#define	ALEPTON_H

#include <stdexcept>
#include <ThObservable.h>
#include "EW.h"

/**
 * @class Alepton
 * @ingroup EW 
 * @brief A class for the left-right asymmetry of @f$Z\to \ell\bar{\ell}@f$ at 
 * the @f$Z@f$ pole.
 * @author SusyFit Collaboration
 * @copyright GNU General Public License
 * @details This class is used to compute the left-right asymmetry of a 
 * @f$Z\to \ell\bar{ell}@f$ channel at the @f$Z@f$ pole. 
 */
class Alepton : public ThObservable {
public:

    /**
     * A constructor.
     * @param[in] EW_i A reference to an object of EW class, which is the base class of 
     * the electroweak precision observables.
     */
    Alepton(const EW& EW_i) 
    : ThObservable(EW_i), myEW(EW_i), myEWTYPE(EW_i.getEWTYPE()) 
    {
    };

    /**
     * @return The left-right asymmetry of a @f$Z\to \ell\bar{\ell}@f$ channel 
     * at the @f$Z@f$ pole, where @f$\ell@f$ denotes a charged-lepton.
     */
    double getThValue();

    
private:

    /**
     * A reference to an object of EW class, which is the base class of the electroweak 
     * precision observables.
     */
    const EW& myEW;

    /**
     * An enumerator controlling the formulae used in the computation.
     */
    const EW::EWTYPE myEWTYPE;
};


#endif	/* ALEPTON_H */

