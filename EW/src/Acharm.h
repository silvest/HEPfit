/* 
 * Copyright (C) 2012-2014 SusyFit Collaboration
 * All rights reserved.
 *
 * For the licensing terms see doc/COPYING.
 */

#ifndef ACHARM_H
#define	ACHARM_H

#include <stdexcept>
#include <ThObservable.h>
#include "EW.h"

/**
 * @class Acharm
 * @ingroup EW 
 * @brief An observable class for the left-right asymmetry in
 * @f$e^+e^-\to Z\to c\bar{c}@f$ at the @f$Z@f$ pole.
 * @author SusyFit Collaboration
 * @copyright GNU General Public License
 * @details This class is used to compute the left-right asymmetry for the
 * @f$e^+e^-\to Z\to c\bar{c}@f$ process at the @f$Z@f$ pole.
 *
 * @sa EW_NPZff::Acharm() and the detailed description of EW class
 * for the inclusion of new physics contribution
 *
 */
class Acharm : public ThObservable {
public:

    /**
     * @brief Constructor.
     * @param[in] EW_i a reference to an object of type EW
     */
    Acharm(const EW& EW_i) 
    : ThObservable(EW_i), myEW(EW_i)
    {
    };

    /**
     * @brief The left-right asymmetry for @f$e^+e^-\to Z\to c\bar{c}@f$
     * at the @f$Z@f$ pole, @f$\mathcal{A}_c@f$.
     * @return @f$\mathcal{A}_c@f$
     */
    double computeThValue();

    
private:
    const EW& myEW;///< A reference to an object of type EW.

};

#endif	/* ACHARM_H */

