/* 
 * Copyright (C) 2015 SusyFit Collaboration
 * All rights reserved.
 *
 * For the licensing terms see doc/COPYING.
 */

#ifndef LAMBDA4_H
#define	LAMBDA4_H

#include <stdexcept>
#include <ThObservable.h>
#include "THDM.h"

/**
 * @class lambda4
 * @ingroup THDM 
 * @brief An observable class for the quartic Higgs potential coupling @f$\lambda_4@f$.
 * @author SusyFit Collaboration
 * @copyright GNU General Public License
 * @details This class is used to compute the quartic Higgs potential coupling @f$\lambda_4@f$.
 */
class lambda4 : public ThObservable {
public:

    /**
     * @brief Constructor.
     * @param[in] ?
     */
    lambda4(const StandardModel& SM_i) 
    : ThObservable(SM_i), myTHDM(static_cast<const THDM*> (&SM_i))
    {
    };

    /**
     * @brief The quartic coupling @f$\lambda_4@f$.
     * @return @f$\lambda_4@f$
     */
    double computeThValue();

    private:
        const THDM * myTHDM;
};

#endif	/* LAMBDA4_H */
