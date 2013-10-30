/*
 * Copyright (C) 2013 SusyFit Collaboration
 * All rights reserved.
 *
 * For the licensing terms see doc/COPYING.
 */

#ifndef C_GAMMAZ_UDS_NP_H
#define	C_GAMMAZ_UDS_NP_H

#include <ThObservable.h>
#include <ThObsType.h>
#include "NPEffective.h"

class c_GammaZ_uds_NP : public ThObservable {
public:

    c_GammaZ_uds_NP(const ThObsType& ObsType)
    : ThObservable(ObsType)
    {
    };

    double computeThValue()
    {
        double delGVu = SM.deltaGVq(SM.UP);
        double delGVd = SM.deltaGVq(SM.DOWN);
        double delGVs = SM.deltaGVq(SM.STRANGE);

        double delGAu = SM.deltaGAq(SM.UP);
        double delGAd = SM.deltaGAq(SM.DOWN);
        double delGAs = SM.deltaGAq(SM.STRANGE);

        double gVu = SM.StandardModel::gVq(SM.UP).real();
        double gVd = SM.StandardModel::gVq(SM.DOWN).real();
        double gVs = SM.StandardModel::gVq(SM.STRANGE).real();

        double gAu = SM.StandardModel::gAq(SM.UP).real();
        double gAd = SM.StandardModel::gAq(SM.DOWN).real();
        double gAs = SM.StandardModel::gAq(SM.STRANGE).real();

        double Lam = (static_cast<const NPEffective*> (&SM))->getLambdaNP();

        return ( (gVu*delGVu + gAu*delGAu + gVd*delGVd + gAd*delGAd
                     + gVs*delGVs + gAs*delGAs )/2.0
                  *Lam*Lam/SM.v()/SM.v() );
    };

private:

};

#endif	/* C_GAMMAZ_UDS_NP_H */

