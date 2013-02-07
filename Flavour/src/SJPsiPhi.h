/* 
 * Copyright (C) 2012 SusyFit Collaboration
 * All rights reserved.
 *
 * For the licensing terms see doc/COPYING.
 */

#ifndef SJPSIPHI_H
#define	SJPSIPHI_H


#include <ThObservable.h>
#include <Flavour.h>
#include <AmpDB2.h>

class SJPsiPhi : public ThObservable, AmpDB2 {
public:
    SJPsiPhi(Flavour& Flavour) : ThObservable(Flavour), AmpDB2(Flavour) {};

    double getThValue() {
        return sin(AmpBs(NLO).arg());
    };
};


#endif	/* SJPSIPHI_H */

