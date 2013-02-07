/* 
 * Copyright (C) 2013 SusyFit Collaboration
 * All rights reserved.
 *
 * For the licensing terms see doc/COPYING.
 */

#ifndef MZ_H
#define	MZ_H

#include <ThObservable.h>
#include <ThObsType.h>
#include "QCD.h"

class mZ : public ThObservable {
public:

    mZ(const ThObsType& ObsType) 
    : ThObservable(ObsType) 
    {
    };
    
    double getThValue()
    {
        return SM.getMz();
    };

private:

};

#endif	/* MZ_H */

