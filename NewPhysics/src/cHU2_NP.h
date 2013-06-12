/*
 * Copyright (C) 2013 SusyFit Collaboration
 * All rights reserved.
 *
 * For the licensing terms see doc/COPYING.
 */

#ifndef CHU2_NP_H
#define	CHU2_NP_H

#include <ThObservable.h>
#include <ThObsType.h>
#include "NPEffective.h"


class cHU2_NP : public ThObservable {
public:

    cHU2_NP(const ThObsType& ObsType)
    : ThObservable(ObsType)
    {
    };

    double getThValue()
    {
        return (static_cast<const NPEffective*> (&SM))->getCHU2();
    };

private:

};

#endif	/* CHU2_NP_H */

