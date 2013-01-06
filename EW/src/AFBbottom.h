/* 
 * Copyright (C) 2012 SusyFit Collaboration
 * All rights reserved.
 *
 * For the licensing terms see doc/COPYING.
 */

#ifndef AFBBOTTOM_H
#define	AFBBOTTOM_H

#include <stdexcept>
#include <ThObservable.h>
#include "EW.h"

/**
 * @class AFBbottom 
 * @ingroup EW
 * @brief A class for the forward-backward asymmetry of @f$Z\to b\bar{b}@f$ 
 * at the @f$Z@f$ pole.
 * @author SusyFit Collaboration
 * @copyright GNU General Public License
 * @details This class is used to compute the forward-backward asymmetry of 
 * the @f$Z\to b\bar{b}@f$ channel at the @f$Z@f$ pole. 
 */
class AFBbottom : public ThObservable {
public:

    /**
     * A constructor.
     * @param[in] EW_i A reference to an object of EW class, which is the base 
     * class of the electroweak precision observables.
     */
    AFBbottom(const EW& EW_i) 
    : ThObservable(EW_i), myEW(EW_i), myEWTYPE(EW_i.getEWTYPE()) 
    {
    };

    /**
     * @return The forward-backward asymmetry of @f$Z\to b\bar{b}@f$ at the 
     * @f$Z@f$ pole.
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

#endif	/* AFBBOTTOM_H */
