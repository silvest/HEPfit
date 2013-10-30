/*
 * Copyright (C) 2013 SusyFit Collaboration
 * All rights reserved.
 *
 * For the licensing terms see doc/COPYING.
 */

#ifndef MSUP_H
#define	MSUP_H

#include <ThObservable.h>
#include <ThObsType.h>
#include "SUSY.h"

/**
 * @class Msup
 * @ingroup SUSY
 * @brief A class for the up-type-squark masses.
 * @author SusyFit Collaboration
 * @copyright GNU General Public License
 * @details
 */
class Msup : public ThObservable {
public:

    Msup(const ThObsType& ObsType, const int ind)
    : ThObservable(ObsType), index(ind)
    {
    };

    double computeThValue()
    {
        return (sqrt((static_cast<const SUSY*> (&SM))->getMsu2()(index)));
    };

private:
    const int index;

};

#endif	/* MSUP_H */

