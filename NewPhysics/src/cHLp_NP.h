/*
 * Copyright (C) 2013 SusyFit Collaboration
 * All rights reserved.
 *
 * For the licensing terms see doc/COPYING.
 */

#ifndef CHLP_NP_H
#define	CHLP_NP_H

#include <ThObservable.h>
#include <ThObsType.h>
#include "NPEffective.h"


class cHLp_NP : public ThObservable {
public:

    cHLp_NP(const ThObsType& ObsType)
    : ThObservable(ObsType)
    {
    };

    double computeThValue()
    {
        double cHL1p = (static_cast<const NPEffective*> (&SM))->getCHL1p();
        double cHL2p = (static_cast<const NPEffective*> (&SM))->getCHL2p();
        double cHL3p = (static_cast<const NPEffective*> (&SM))->getCHL3p();
        if ( (cHL1p == cHL2p) && (cHL2p == cHL3p) ) 
            return cHL1p;
        else
            throw std::runtime_error("cHLp_NP::computeThValue(): No lepton-flavor universality!");
    };

private:

};

#endif	/* CHLP_NP_H */

