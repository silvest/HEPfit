/*
 * Copyright (C) 2013 HEPfit Collaboration
 * All rights reserved.
 *
 * For the licensing terms see doc/COPYING.
 */

#ifndef MSDOWN_H
#define	MSDOWN_H

#include <ThObservable.h>


#include "SUSY.h"

/**
 * @class Msdown
 * @ingroup SUSY
 * @brief A class for the down-type-squark masses.
 * @author HEPfit Collaboration
 * @copyright GNU General Public License
 * @details
 */
class Msdown : public ThObservable {
public:

    Msdown(const StandardModel& SM_i, const int ind)
    : ThObservable(SM_i), index(ind), mySUSY(static_cast<const SUSY*> (&SM_i))
    {
        if (mySUSY->isModelSUSY() == false)
            throw std::runtime_error("\nERROR: The sdown mass spectrum can only be computed in a SUSY model. Please check your observables list.\n");
    };

    double computeThValue()
    {
        return (sqrt(mySUSY->getMsd2()(index)));
    };

private:
    const int index;
    const SUSY * mySUSY;
};

#endif	/* MSDOWN_H */

