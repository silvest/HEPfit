/*
 * Copyright (C) 2013 SusyFit Collaboration
 * All rights reserved.
 *
 * For the licensing terms see doc/COPYING.
 */

#ifndef CHQP_NP_H
#define	CHQP_NP_H

#include <ThObservable.h>
#include <ThObsType.h>
#include "NPEffective.h"


class cHQp_NP : public ThObservable {
public:

    cHQp_NP(const ThObsType& ObsType)
    : ThObservable(ObsType)
    {
    };

    double getThValue()
    {
        double cHQ1p = (static_cast<const NPEffective*> (&SM))->getCHQ1p();
        double cHQ2p = (static_cast<const NPEffective*> (&SM))->getCHQ2p();
        double cHQ3p = (static_cast<const NPEffective*> (&SM))->getCHQ3p();
        if ( (cHQ1p == cHQ2p) && (cHQ2p == cHQ3p) )
            return cHQ1p;
        else
            throw std::runtime_error("cHQp_NP::getThValue(): No quark-flavor universality!");
    };

private:

};

#endif	/* CHQP_NP_H */

