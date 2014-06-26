/*
 * Copyright (C) 2013 SusyFit Collaboration
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
 * @author SusyFit Collaboration
 * @copyright GNU General Public License
 * @details
 */
class Mneutralino : public ThObservable {
public:

    Mneutralino(const ThObsType& ObsType, const int ind)
    : ThObservable(ObsType), index(ind)
    {
    };

    double computeThValue()
    {
        return ((static_cast<const SUSY*> (&SM))->getMneu()(index));
    };
    
private:
    const int index;

};

#endif	/* MNEUTRALINO_H */

