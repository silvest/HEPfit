/*
 * Copyright (C) 2013 HEPfit Collaboration
 * All rights reserved.
 *
 * For the licensing terms see doc/COPYING.
 */

#ifndef MSUP_H
#define	MSUP_H

#include <ThObservable.h>


#include "SUSY.h"

/**
 * @class Msup
 * @ingroup SUSY
 * @brief A class for the up-type-squark masses.
 * @author HEPfit Collaboration
 * @copyright GNU General Public License
 * @details
 */
class Msup : public ThObservable {
public:

    Msup(const StandardModel& SM_i, const int ind)
    : ThObservable(SM_i), index(ind), mySUSY(static_cast<const SUSY*> (&SM_i))
    {
        if (mySUSY->isModelSUSY() == false)
            throw std::runtime_error("\nERROR: The msup mass spectrum can only be computed in a SUSY model. Please check your observables list.\n");
    };

    double computeThValue()
    {
        return (sqrt(mySUSY->getMsu2()(index)));
    };

private:
    const int index;
    const SUSY * mySUSY;
};

#endif	/* MSUP_H */

