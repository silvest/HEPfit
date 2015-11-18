/*
 * Copyright (C) 2015 HEPfit Collaboration
 * All rights reserved.
 *
 * For the licensing terms see doc/COPYING.
 */

#ifndef MSNEUTRINO_H
#define	MSNEUTRINO_H

#include <ThObservable.h>


#include "SUSY.h"

/**
 * @class Msneutrino
 * @ingroup SUSY
 * @brief A class for the sneutrino masses.
 * @author HEPfit Collaboration
 * @copyright GNU General Public License
 * @details
 */
class Msneutrino : public ThObservable {
public:

    Msneutrino(const StandardModel& SM_i, const int ind)
    : ThObservable(SM_i), index(ind), mySUSY(static_cast<const SUSY*> (&SM_i))
    {
        if (mySUSY->isModelSUSY() == false)
            throw std::runtime_error("\nERROR: The msneutrino mass spectrum can only be computed in a SUSY model. Please check your observables list.\n");
    };

    double computeThValue()
    {
        return (sqrt(mySUSY->getMsn2()(index)));
    };

private:
    const int index;
    const SUSY * mySUSY;
};

#endif	/* MSNEUTRINO_H */
