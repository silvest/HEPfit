/* 
 * Copyright (C) 2013 SusyFit Collaboration
 * All rights reserved.
 *
 * For the licensing terms see doc/COPYING.
 */

#ifndef DELTAGVB_H
#define	DELTAGVB_H

#include <ThObservable.h>
#include <ThObsType.h>
#include "NPZbbbar.h"

class deltaGVb : public ThObservable {
public:

    deltaGVb(const ThObsType& ObsType) 
    : ThObservable(ObsType) 
    {
    };
    
    double getThValue()
    {
        return SM.deltaGVb();
    };
    
private:

};

#endif	/* DELTAGVB_H */

