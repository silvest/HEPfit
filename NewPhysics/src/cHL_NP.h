/*
 * Copyright (C) 2013 SusyFit Collaboration
 * All rights reserved.
 *
 * For the licensing terms see doc/COPYING.
 */

#ifndef CHL_NP_H
#define	CHL_NP_H

#include <ThObservable.h>
#include <ThObsType.h>
#include "NPEffective.h"


class cHL_NP : public ThObservable {
public:

    cHL_NP(const ThObsType& ObsType)
    : ThObservable(ObsType)
    {
    };

    double computeThValue()
    {
        double cHL1 = (static_cast<const NPEffective*> (&SM))->getCHL1();
        double cHL2 = (static_cast<const NPEffective*> (&SM))->getCHL2();
        double cHL3 = (static_cast<const NPEffective*> (&SM))->getCHL3();
        if ( (cHL1 == cHL2) && (cHL2 == cHL3) )
            return cHL1;
        else
            throw std::runtime_error("cHL_NP::computeThValue(): No lepton-flavor universality!");
    };

private:

};

#endif	/* CHL_NP_H */

