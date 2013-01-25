/* 
 * Copyright (C) 2013 SusyFit Collaboration
 * All rights reserved.
 *
 * For the licensing terms see doc/COPYING.
 */

#ifndef DELTAGRB_H
#define	DELTAGRB_H

#include <ThObservable.h>
#include <ThObsType.h>
#include "NPZbbbar.h"

class deltaGRb : public ThObservable {
public:

    deltaGRb(const ThObsType& ObsType) 
    : ThObservable(ObsType) 
    {
    };
    
    double getThValue()
    {
        return ( (SM.deltaGVb() - SM.deltaGAb())/2.0 );
    };

private:

};

#endif	/* DELTAGRB_H */

