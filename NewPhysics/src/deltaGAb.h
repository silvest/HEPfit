/* 
 * Copyright (C) 2013 SusyFit Collaboration
 * All rights reserved.
 *
 * For the licensing terms see doc/COPYING.
 */

#ifndef DELTAGAB_H
#define	DELTAGAB_H

#include <ThObservable.h>
#include <ThObsType.h>
#include "NPZbbbar.h"

class deltaGAb : public ThObservable {
public:

    deltaGAb(const ThObsType& ObsType) 
    : ThObservable(ObsType) 
    {
    };
    
    double getThValue()
    {
        return SM.deltaGAb();
    };
    
private:

};

#endif	/* DELTAGAB_H */

