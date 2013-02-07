/* 
 * Copyright (C) 2012 SusyFit Collaboration
 * All rights reserved.
 *
 * For the licensing terms see doc/COPYING.
 */

#ifndef GAMMA_H
#define	GAMMA_H

#include <ThObservable.h>
#include <ThObsType.h>

class CKMGamma : public ThObservable {
public:
    CKMGamma(const ThObsType& ObsType) : ThObservable(ObsType) {};

    double getThValue() { 
        return(SM.getCKM().getGamma()/M_PI*180.);
    };
};

#endif	/* GAMMA_H */
