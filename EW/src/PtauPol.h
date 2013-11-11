/* 
 * Copyright (C) 2012-2013 SusyFit Collaboration
 * All rights reserved.
 *
 * For the licensing terms see doc/COPYING.
 */

#ifndef PTAUPOL_H
#define	PTAUPOL_H

#include <stdexcept>
#include <ThObservable.h>
#include "EW.h"

/**
 * @class PtauPol
 * @ingroup EW 
 * @brief A class for the longitudinal polarization of @f$Z\to \tau\bar{\tau}@f$ 
 * at the @f$Z@f$ pole.
 * @author SusyFit Collaboration
 * @copyright GNU General Public License
 * @details This class is used to compute the longitudinal polarization of 
 * @f$Z\to \tau\bar{\tau}@f$ at the @f$Z@f$ pole.
 */
class PtauPol : public ThObservable {
public:

    /**
     * A constructor.
     * @param[in] EW_i A reference to an object of EW class, which is the base 
     * class of the electroweak precision observables.
     */
    PtauPol(const EW& EW_i) 
    : ThObservable(EW_i), myEW(EW_i)
    {
    };

    /**
     * @return The longitudinal polarization of @f$Z\to \tau\bar{\tau}@f$ at the 
     * @f$Z@f$ pole.
     */
    double computeThValue();

    
private:

    /**
     * A reference to an object of EW class, which is the base class of the 
     * electroweak precision observables.
     */
    const EW& myEW;
};

#endif	/* PTAUPOL_H */

