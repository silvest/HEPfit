/*
 * Copyright (C) 2013 SusyFit Collaboration
 * All rights reserved.
 *
 * For the licensing terms see doc/COPYING.
 */

#ifndef C_AE_NP_H
#define	C_AE_NP_H

#include <ThObservable.h>
#include <ThObsType.h>
#include "NPEffective.h"

class c_Ae_NP : public ThObservable {
public:

    c_Ae_NP(const ThObsType& ObsType)
    : ThObservable(ObsType)
    {
    };

    double computeThValue()
    {
        double delGVe = SM.deltaGVl(SM.ELECTRON);
        double delGAe = SM.deltaGAl(SM.ELECTRON);

        double gVe = SM.StandardModel::gVl(SM.ELECTRON).real();
        double gAe = SM.StandardModel::gAl(SM.ELECTRON).real();

        double Lam = (static_cast<const NPEffective*> (&SM))->getLambdaNP();

        return ( (gAe*delGVe - gVe*delGAe)/2.0*Lam*Lam/SM.v()/SM.v() );
    };

private:

};

#endif	/* C_AE_NP_H */

