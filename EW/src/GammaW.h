/* 
 * Copyright (C) 2012-2014 SusyFit Collaboration
 * All rights reserved.
 *
 * For the licensing terms see doc/COPYING.
 */

#ifndef GAMMAW_H
#define	GAMMAW_H

#include <stdexcept>
#include <ThObservable.h>
#include "EW.h"

/**
 * @class GammaW 
 * @ingroup EW 
 * @brief An observable class for the total decay width of the @f$W@f$ boson.
 * @author SusyFit Collaboration
 * @copyright GNU General Public License
 * @details This class is used to compute the total decay width of the @f$W@f$
 * boson, @f$\Gamma_W@f$.
 */
class GammaW : public ThObservable {
public:
    
    /**
     * @brief Constructor.
     * @param[in] EW_i a reference to an object of type EW
     */
    GammaW(const EW& EW_i) 
    : ThObservable(EW_i)
    {
    };

    /**
     * @brief The total decay width of the @f$W@f$ boson, @f$\Gamma_W@f$,
     * in units of GeV.
     * @return @f$\Gamma_W@f$ in units of GeV
     */
    double computeThValue();

};

#endif	/* GAMMAW_H */

