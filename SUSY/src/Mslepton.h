/*
 * Copyright (C) 2015 HEPfit Collaboration
 * All rights reserved.
 *
 * For the licensing terms see doc/COPYING.
 */

#ifndef MSLEPTON_H
#define	MSLEPTON_H

#include <ThObservable.h>


#include "SUSY.h"

/**
 * @class Mslepton
 * @ingroup SUSY
 * @brief A class for the charged slepton masses.
 * @author HEPfit Collaboration
 * @copyright GNU General Public License
 * @details
 */
class Mslepton : public ThObservable {
public:

    Mslepton(const StandardModel& SM_i, const int ind)
    : ThObservable(SM_i), index(ind), mySUSY(static_cast<const SUSY*> (&SM_i))
    {
        if (mySUSY->isModelSUSY() == false)
            throw std::runtime_error("\nERROR: The mslepton mass spectrum can only be computed in a SUSY model. Please check your observables list.\n");
    };

    double computeThValue()
    {
        return (sqrt(mySUSY->getMse2()(index)));
    };

private:
    const int index;
    const SUSY * mySUSY;
};

#endif	/* MSLEPTON_H */
