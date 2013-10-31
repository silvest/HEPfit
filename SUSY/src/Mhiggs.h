/*
 * Copyright (C) 2013 SusyFit Collaboration
 * All rights reserved.
 *
 * For the licensing terms see doc/COPYING.
 */

#ifndef MHIGGS_H
#define	MHIGGS_H

#include <stdexcept>
#include <ThObservable.h>
#include <ThObsType.h>
#include "SUSY.h"

/**
 * @class Mhiggs
 * @ingroup SUSY
 * @brief A class for the Higgs masses.
 * @author SusyFit Collaboration
 * @copyright GNU General Public License
 * @details
 */
class Mhiggs : public ThObservable {
public:

    Mhiggs(const ThObsType& ObsType, const int ind)
    : ThObservable(ObsType), index(ind)
    {
    };

    double computeThValue()
    {
        switch(index) {
            case 0:
                return (static_cast<const SUSY*> (&SM))->getMHl();
            case 1:
                return (static_cast<const SUSY*> (&SM))->getMHh();
            case 2:
                return (static_cast<const SUSY*> (&SM))->getMHa();
            case 3:
                return (static_cast<const SUSY*> (&SM))->getMHp();
            default:
                throw std::runtime_error("Mhiggs::computeThValue(): undefined index");
        }
    };

private:
    const int index;
};



#endif	/* MHIGGS_H */

