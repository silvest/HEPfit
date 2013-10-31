/*
 * Copyright (C) 2013 SusyFit Collaboration
 * All rights reserved.
 *
 * For the licensing terms see doc/COPYING.
 */

#ifndef CHD3_NP_H
#define	CHD3_NP_H

#include <ThObservable.h>
#include <ThObsType.h>
#include "NPEffective.h"


class cHD3_NP : public ThObservable {
public:

    cHD3_NP(const ThObsType& ObsType)
    : ThObservable(ObsType)
    {
    };

    double computeThValue()
    {
        return (static_cast<const NPEffective*> (&SM))->getCHD3();
    };

private:

};

#endif	/* CHD3_NP_H */

