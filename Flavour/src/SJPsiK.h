/* 
 * Copyright (C) 2012 SusyFit Collaboration
 * All rights reserved.
 *
 * For the licensing terms see doc/COPYING.
 */

#ifndef SJPSIK_H
#define	SJPSIK_H

#include <ThObservable.h>
#include <Flavour.h>
#include <AmpDB2.h>

class SJPsiK : public ThObservable, AmpDB2 {
public:
    SJPsiK(Flavour& Flavour) : ThObservable(Flavour), AmpDB2(Flavour) {};

    double getThValue() {
        return sin(-AmpBd(NLO).arg());
    };
};


#endif	/* SJPSIK_H */

