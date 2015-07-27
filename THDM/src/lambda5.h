/* 
 * Copyright (C) 2015 SusyFit Collaboration
 * All rights reserved.
 *
 * For the licensing terms see doc/COPYING.
 */

#ifndef LAMBDA5_H
#define	LAMBDA5_H

#include <stdexcept>
#include <ThObservable.h>
#include "THDM.h"

/**
 * @class lambda5
 * @ingroup THDM 
 * @brief An observable class for the quartic Higgs potential coupling @f$\lambda_5@f$.
 * @author SusyFit Collaboration
 * @copyright GNU General Public License
 * @details This class is used to compute the quartic Higgs potential coupling @f$\lambda_5@f$.
 */
class lambda5 : public ThObservable {
public:

    /**
     * @brief Constructor.
     * @param[in] ?
     */
    lambda5(const StandardModel& SM_i) 
    : ThObservable(SM_i), myTHDM(static_cast<const THDM*> (&SM_i))
    {
    };

    /**
     * @brief The quartic coupling @f$\lambda_5@f$.
     * @return @f$\lambda_5@f$
     */
    double computeThValue();

    private:
        const THDM * myTHDM;
};

#endif	/* LAMBDA5_H */



