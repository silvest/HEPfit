/* 
 * Copyright (C) 2013 SusyFit Collaboration
 * All rights reserved.
 *
 * For the licensing terms see doc/COPYING.
 */

#ifndef MHL_H
#define	MHL_H

#include <ThObservable.h>
#include <ThObsType.h>
#include "StandardModel.h"

class mHl : public ThObservable {
public:

    mHl(const ThObsType& ObsType) 
    : ThObservable(ObsType) 
    {
    };
    
    double getThValue()
    {
        return SM.getMHl();
    };

private:

};

#endif	/* MHL_H */

