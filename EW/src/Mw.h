/* 
 * Copyright (C) 2012 HEPfit Collaboration
 * All rights reserved.
 *
 * For the licensing terms see doc/COPYING.
 */

#ifndef MW_H
#define	MW_H

#include <stdexcept>
#include <ThObservable.h>

/**
 * @class Mw 
 * @ingroup EW 
 * @brief An observable class for the @f$W@f$-boson mass.
 * @author HEPfit Collaboration
 * @copyright GNU General Public License
 * @details This class is used to compute the @f$W@f$-boson mass.
 */
class Mw : public ThObservable {
public:

    /**
     * @brief Constructor.
     * @param[in] SM_i a reference to an object of type StandardModel
     */
    Mw(const StandardModel& SM_i) 
    : ThObservable(SM_i)
    {
    };

    /**
     * @brief The @f$W@f$-boson mass, @f$M_W@f$, in units of GeV.
     * @return @f$M_W@f$ in units of GeV
     */
    double computeThValue();

};

#endif	/* MW_H */

