/* 
 * Copyright (C) 2015 HEPfit Collaboration
 * All rights reserved.
 *
 * For the licensing terms see doc/COPYING.
 */

#ifndef GLOBALMINIMUM_H
#define	GLOBALMINIMUM_H

#include <stdexcept>
#include "ThObservable.h"
#include "THDM.h"

/**
 * @class globalminimum 
 * @ingroup THDM 
 * @brief An observable class to prevent the occurence of a second deeper vacuum.
 * @author HEPfit Collaboration
 * @copyright GNU General Public License
 * @details This class is used to assure that our vacuum is the global minimum of the potential.
 * Modified discriminant from @cite Barroso:2013awa. (Scaled by a factor @f$10^{-7}@f$ and some always positive factors removed.)
 */
class globalminimum : public ThObservable {
public:

    /**
     * @brief Constructor.
     * @param[in] ?
     */
    globalminimum(const StandardModel& SM_i) 
    : ThObservable(SM_i), myTHDM(static_cast<const THDM*> (&SM_i))
    {
    };

    /**
     * @brief A condition to assure that our vacuum is global.
     * @return 
     */
    double computeThValue();

    private:
        const THDM * myTHDM;
};

#endif	/* GLOBALMINIMUM_H */

