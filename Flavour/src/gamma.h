/* 
 * Copyright (C) 2012 HEPfit Collaboration
 * All rights reserved.
 *
 * For the licensing terms see doc/COPYING.
 */

#ifndef GAMMA_H
#define	GAMMA_H

#include <ThObservable.h>



class CKMGamma : public ThObservable {
public:
    CKMGamma(const StandardModel& SM_i) : ThObservable(SM_i) {};

    double computeThValue() { 
        return(SM.computeGamma()/M_PI*180.);
    };
};

#endif	/* GAMMA_H */
