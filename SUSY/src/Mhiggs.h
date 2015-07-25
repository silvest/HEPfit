/*
 * Copyright (C) 2013 HEPfit Collaboration
 * All rights reserved.
 *
 * For the licensing terms see doc/COPYING.
 */

#ifndef MHIGGS_H
#define	MHIGGS_H

#include <stdexcept>
#include <ThObservable.h>
#include "SUSY.h"

/**
 * @class Mhiggs
 * @ingroup SUSY
 * @brief A class for the Higgs masses.
 * @author HEPfit Collaboration
 * @copyright GNU General Public License
 * @details
 */
class Mhiggs : public ThObservable {
public:

    Mhiggs(const StandardModel& SM_i, const int ind)
    : ThObservable(SM_i), index(ind), mySUSY(static_cast<const SUSY*> (&SM_i))
    {
        if (mySUSY->isModelSUSY() == false)
            throw std::runtime_error("\nERROR: The Higgs mass spectrum can only be computed in a SUSY model. Please check your observables list.\n");
    };

    double computeThValue()
    {
        switch(index) {
            case 0:
                return mySUSY->getMHl();
            case 1:
                return mySUSY->getMHh();
            case 2:
                return mySUSY->getMHa();
            case 3:
                return mySUSY->getMHp();
            default:
                throw std::runtime_error("Mhiggs::computeThValue(): undefined index");
        }
    };

private:
    const int index;
    const SUSY * mySUSY;
};



#endif	/* MHIGGS_H */

