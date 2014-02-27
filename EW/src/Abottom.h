/* 
 * Copyright (C) 2012 SusyFit Collaboration
 * All rights reserved.
 *
 * For the licensing terms see doc/COPYING.
 */

#ifndef ABOTTOM_H
#define	ABOTTOM_H

#include <stdexcept>
#include <ThObservable.h>
#include "EW.h"

/**
 * @addtogroup EW
 * @brief A module for electroweak precision observables.
 * @details 
 * @{
 */

/**
 * @class Abottom 
 * @brief An observable class for the left-right asymmetry in
 * @f$e^+e^-\to Z\to b\bar{b}@f$ at the @f$Z@f$ pole.
 * @author SusyFit Collaboration
 * @copyright GNU General Public License
 * @details This class is used to compute the left-right asymmetry for the
 * @f$e^+e^-\to Z\to b\bar{b}@f$ process at the @f$Z@f$ pole.
 * 
 * @sa EW_NPZff::Abottom() and the detailed description of EW class
 * for the inclusion of new physics contribution
 *
 */
class Abottom : public ThObservable {
public:

    /**
     * @brief Constructor.
     * @param[in] EW_i a reference to an object of type EW
     */
    Abottom(const EW& EW_i) 
    : ThObservable(EW_i), myEW(EW_i)
    {
    };

    /**
     * @brief The left-right asymmetry for @f$e^+e^-\to Z\to b\bar{b}@f$
     * at the @f$Z@f$ pole, @f$\mathcal{A}_b@f$.
     * @return @f$\mathcal{A}_b@f$
     */
    double computeThValue();

    
private:
    const EW& myEW;///< A reference to an object of type EW.

};

/** 
 * @}
 */

#endif	/* ABOTTOM_H */

