/*
 * Copyright (C) 2013 SusyFit Collaboration
 * All rights reserved.
 *
 * For the licensing terms see doc/COPYING.
 */

#ifndef CHQ2PMINUSCHQ2_NP_H
#define	CHQ2PMINUSCHQ2_NP_H

#include <ThObservable.h>
#include <ThObsType.h>
#include "NPEffective.h"

class cHQ2pMINUScHQ2_NP : public ThObservable {
public:

    cHQ2pMINUScHQ2_NP(const ThObsType& ObsType)
    : ThObservable(ObsType)
    {
    };

    double computeThValue()
    {
        return ( (static_cast<const NPEffective*> (&SM))->getCHQ2p()
                  - (static_cast<const NPEffective*> (&SM))->getCHQ2() );
    };

private:

};

#endif	/* CHQ2PMINUSCHQ2_NP_H */

