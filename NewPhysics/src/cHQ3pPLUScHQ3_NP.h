/*
 * Copyright (C) 2013 SusyFit Collaboration
 * All rights reserved.
 *
 * For the licensing terms see doc/COPYING.
 */

#ifndef CHQ3PPLUSCHQ3_NP_H
#define	CHQ3PPLUSCHQ3_NP_H

#include <ThObservable.h>
#include <ThObsType.h>
#include "NPEffective.h"

class cHQ3pPLUScHQ3_NP : public ThObservable {
public:

    cHQ3pPLUScHQ3_NP(const ThObsType& ObsType)
    : ThObservable(ObsType)
    {
    };

    double computeThValue()
    {
        return ( (static_cast<const NPEffective*> (&SM))->getCHQ3p()
                  + (static_cast<const NPEffective*> (&SM))->getCHQ3() );
    };

private:

};

#endif	/* CHQ3PPLUSCHQ3_NP_H */

