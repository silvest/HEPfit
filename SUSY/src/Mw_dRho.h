/*
 * Copyright (C) 2013 SusyFit Collaboration
 * All rights reserved.
 *
 * For the licensing terms see doc/COPYING.
 */

#ifndef MW_DRHO_H
#define	MW_DRHO_H

#include <ThObservable.h>
#include <ThObsType.h>
#include "SUSY.h"

class Mw_dRho : public ThObservable {
public:

    Mw_dRho(const ThObsType& ObsType)
    : ThObservable(ObsType)
    {
    };

    double getThValue()
    {
        return (static_cast<const SUSY*> (&SM))->Mw_dRho();
    };

private:

};

#endif	/* MW_DRHO_H */

