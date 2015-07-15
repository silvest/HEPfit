/*
 * Copyright (C) 2013 HEPfit Collaboration
 * All rights reserved.
 *
 * For the licensing terms see doc/COPYING.
 */

#ifndef MCHARGINO_H
#define	MCHARGINO_H

#include <ThObservable.h>
#include "SUSY.h"

/**
 * @class Mchargino
 * @ingroup SUSY
 * @brief A class for the chargino masses.
 * @author HEPfit Collaboration
 * @copyright GNU General Public License
 * @details
 */
class Mchargino : public ThObservable {
public:

    Mchargino(const StandardModel& SM_i, const int ind)
    : ThObservable(SM_i), index(ind), mySUSY(static_cast<const SUSY*> (&SM_i))
    {
        if (mySUSY->isModelSUSY() == false)
            throw std::runtime_error("\nERROR: The chargino mass spectrum can only be computed in a SUSY model. Please check your observables list.\n");
    };

    double computeThValue()
    {
        return (mySUSY->getMch()(index));
    };

private:
    const int index;
    const SUSY * mySUSY;
};

#endif	/* MCHARGINO_H */

