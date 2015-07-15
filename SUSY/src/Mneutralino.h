/*
 * Copyright (C) 2013 HEPfit Collaboration
 * All rights reserved.
 *
 * For the licensing terms see doc/COPYING.
 */

#ifndef MNEUTRALINO_H
#define	MNEUTRALINO_H

#include <ThObservable.h>


#include "SUSY.h"

/**
 * @class Mneutralino
 * @ingroup SUSY
 * @brief A class for the neutralino masses.
 * @author HEPfit Collaboration
 * @copyright GNU General Public License
 * @details
 */
class Mneutralino : public ThObservable {
public:

    Mneutralino(const StandardModel& SM_i, const int ind)
    : ThObservable(SM_i), index(ind), mySUSY(static_cast<const SUSY*> (&SM_i))
    {
        if (mySUSY->isModelSUSY() == false)
            throw std::runtime_error("\nERROR: The neutralino mass spectrum can only be computed in a SUSY model. Please check your observables list.\n");
    };

    double computeThValue()
    {
        return (mySUSY->getMneu()(index));
    };
    
private:
    const int index;
    const SUSY * mySUSY;

};

#endif	/* MNEUTRALINO_H */

