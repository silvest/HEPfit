/* 
 * Copyright (C) 2012-2013 SusyFit Collaboration
 * All rights reserved.
 *
 * For the licensing terms see doc/COPYING.
 */

#ifndef AFBLEPTON_H
#define	AFBLEPTON_H

#include <stdexcept>
#include <ThObservable.h>
#include "EW.h"

/**
 * @class AFBlepton
 * @ingroup EW 
 * @brief A class for the forward-backward asymmetry of @f$Z\to \ell\bar{\ell}@f$ 
 * at the @f$Z@f$ pole.
 * @author SusyFit Collaboration
 * @copyright GNU General Public License
 * @details This class is used to compute the forward-backward asymmetry of 
 * a @f$Z\to \ell\bar{\ell}@f$ channel at the @f$Z@f$ pole. 
 */
class AFBlepton : public ThObservable {
public:
 
    /**
     * A constructor.
     * @param[in] EW_i A reference to an object of EW class, which is the base class of 
     * the electroweak precision observables.
     */
    AFBlepton(const EW& EW_i) 
    : ThObservable(EW_i), myEW(EW_i)
    {
    };

    /**
     * @return The forward-backward asymmetry of a @f$Z\to \ell\bar{\ell}@f$ 
     * channel at the @f$Z@f$ pole, where @f$\ell@f$ denotes a charged-lepton.
     */
    double computeThValue();

    
private:

    /**
     * A reference to an object of EW class, which is the base class of the electroweak 
     * precision observables.
     */
    const EW& myEW;
};

#endif	/* AFBLEPTON_H */

