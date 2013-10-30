/*
 * Copyright (C) 2013 SusyFit Collaboration
 * All rights reserved.
 *
 * For the licensing terms see doc/COPYING.
 */

#ifndef CHQ1PPLUSCHQ2P_NP_H
#define	CHQ1PPLUSCHQ2P_NP_H

#include <ThObservable.h>
#include <ThObsType.h>
#include "NPEffective.h"

class cHQ1pPLUScHQ2p_NP : public ThObservable {
public:

    cHQ1pPLUScHQ2p_NP(const ThObsType& ObsType)
    : ThObservable(ObsType)
    {
    };

    double computeThValue()
    {
        return ( (static_cast<const NPEffective*> (&SM))->getCHQ1p()
                  + (static_cast<const NPEffective*> (&SM))->getCHQ2p() );
    };

private:

};


#endif	/* CHQ1PPLUSCHQ2P_NP_H */

