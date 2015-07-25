/* 
 * Copyright (C) 2012 HEPfit Collaboration
 * All rights reserved.
 *
 * For the licensing terms see doc/COPYING.
 */

#ifndef GAMMAW_H
#define	GAMMAW_H

#include <stdexcept>
#include <ThObservable.h>

/**
 * @class GammaW 
 * @ingroup EW 
 * @brief An observable class for the total decay width of the @f$W@f$ boson.
 * @author HEPfit Collaboration
 * @copyright GNU General Public License
 * @details This class is used to compute the total decay width of the @f$W@f$
 * boson, @f$\Gamma_W@f$.
 */
class GammaW : public ThObservable {
public:
    
    /**
     * @brief Constructor.
     * @param[in] SM_i a reference to an object of type StandardModel
     */
    GammaW(const StandardModel& SM_i) 
    : ThObservable(SM_i)
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

