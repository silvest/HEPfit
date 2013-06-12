/*
 * Copyright (C) 2013 SusyFit Collaboration
 * All rights reserved.
 *
 * For the licensing terms see doc/COPYING.
 */

#ifndef CHE_NP_H
#define	CHE_NP_H

#include <ThObservable.h>
#include <ThObsType.h>
#include "NPEffective.h"


class cHE_NP : public ThObservable {
public:

    cHE_NP(const ThObsType& ObsType)
    : ThObservable(ObsType)
    {
    };

    double getThValue()
    {
        double cHE1 = (static_cast<const NPEffective*> (&SM))->getCHE1();
        double cHE2 = (static_cast<const NPEffective*> (&SM))->getCHE2();
        double cHE3 = (static_cast<const NPEffective*> (&SM))->getCHE3();
        if ( (cHE1 == cHE2) && (cHE2 == cHE3) )
            return cHE1;
        else
            throw std::runtime_error("cHE_NP::getThValue(): No lepton-flavor universality!");
    };

private:

};

#endif	/* CHE_NP_H */

