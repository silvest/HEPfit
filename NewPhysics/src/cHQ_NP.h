/*
 * Copyright (C) 2013 SusyFit Collaboration
 * All rights reserved.
 *
 * For the licensing terms see doc/COPYING.
 */

#ifndef CHQ_NP_H
#define	CHQ_NP_H

#include <ThObservable.h>
#include <ThObsType.h>
#include "NPEffective.h"


class cHQ_NP : public ThObservable {
public:

    cHQ_NP(const ThObsType& ObsType)
    : ThObservable(ObsType)
    {
    };

    double computeThValue()
    {
        double cHQ1 = (static_cast<const NPEffective*> (&SM))->getCHQ1();
        double cHQ2 = (static_cast<const NPEffective*> (&SM))->getCHQ2();
        double cHQ3 = (static_cast<const NPEffective*> (&SM))->getCHQ3();
        if ( (cHQ1 == cHQ2) && (cHQ2 == cHQ3) )
            return cHQ1;
        else
            throw std::runtime_error("cHQ_NP::computeThValue(): No quark-flavor universality!");
    };

private:

};

#endif	/* CHQ_NP_H */

