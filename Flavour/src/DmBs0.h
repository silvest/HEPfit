/* 
 * Copyright (C) 2012 SUSYfit Collaboration
 * All rights reserved.
 *
 * For the licensing terms see doc/COPYING.
 */

#ifndef DMBS0_H
#define	DMBS0_H

#include <ThObservable.h>
#include <Flavour.h>
#include <AmpDB2.h>

class DmBs0 : public ThObservable, AmpDB2 {
public:
    DmBs0(Flavour& Flavour) : ThObservable(Flavour), AmpDB2(Flavour) {};

    double getThValue() {
        return(2.*AmpBs(LO).abs());
    };
};



#endif	/* DMBS0_H */

