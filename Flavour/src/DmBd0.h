/* 
 * Copyright (C) 2012 SUSYfit Collaboration
 * All rights reserved.
 *
 * For the licensing terms see doc/COPYING.
 */

#ifndef DMBD0_H
#define	DMBD0_H

#include <ThObservable.h>
#include <Flavour.h>
#include <AmpDB2.h>

class DmBd0 : public ThObservable, AmpDB2 {
public:
    DmBd0(Flavour& Flavour) : ThObservable(Flavour), AmpDB2(Flavour) {};

    double getThValue() {
        return(2.*AmpBd(LO).abs());
    };
};

#endif	/* DMBD0_H */
