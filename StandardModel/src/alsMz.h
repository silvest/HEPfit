/* 
 * Copyright (C) 2013 SusyFit Collaboration
 * All rights reserved.
 *
 * For the licensing terms see doc/COPYING.
 */

#ifndef ALSMZ_H
#define	ALSMZ_H

#include <ThObservable.h>
#include <ThObsType.h>
#include "QCD.h"

class alsMz : public ThObservable {
public:

    alsMz(const ThObsType& ObsType) 
    : ThObservable(ObsType) 
    {
    };
    
    double getThValue()
    {
        return SM.getAlsMz();
    };

private:

};

#endif	/* ALSMZ_H */

