/* 
 * Copyright (C) 2015 SusyFit Collaboration
 * All rights reserved.
 *
 * For the licensing terms see doc/COPYING.
 */

#ifndef LAMBDA3_H
#define	LAMBDA3_H

#include <stdexcept>
#include <ThObservable.h>
#include "THDM.h"

/**
 * @class lambda3
 * @ingroup THDM 
 * @brief An observable class for the quartic Higgs potential coupling @f$\lambda_3@f$.
 * @author SusyFit Collaboration
 * @copyright GNU General Public License
 * @details This class is used to compute the quartic Higgs potential coupling @f$\lambda_3@f$.
 */
class lambda3 : public ThObservable {
public:

    /**
     * @brief Constructor.
     * @param[in] ?
     */
    lambda3(const StandardModel& SM_i) 
    : ThObservable(SM_i), myTHDM(static_cast<const THDM*> (&SM_i))
    {
    };

    /**
     * @brief The quartic coupling @f$\lambda_3@f$.
     * @return @f$\lambda_3@f$
     */
    double computeThValue();

    private:
        const THDM * myTHDM;
};

#endif	/* LAMBDA3_H */


