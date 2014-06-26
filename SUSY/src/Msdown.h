/*
 * Copyright (C) 2013 SusyFit Collaboration
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
 * @author SusyFit Collaboration
 * @copyright GNU General Public License
 * @details
 */
class Msdown : public ThObservable {
public:

    Msdown(const ThObsType& ObsType, const int ind)
    : ThObservable(ObsType), index(ind)
    {
    };

    double computeThValue()
    {
        return (sqrt((static_cast<const SUSY*> (&SM))->getMsd2()(index)));
    };

private:
    const int index;
    
};

#endif	/* MSDOWN_H */

