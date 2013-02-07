/* 
 * Copyright (C) 2012 SusyFit Collaboration
 * All rights reserved.
 *
 * For the licensing terms see doc/COPYING.
 */

#ifndef MU2_H
#define	MU2_H

#include <ThObservable.h>
#include <ThObsType.h>
#include "SUSY.h"

class Msu2 : public ThObservable {
public:
    Msu2(const ThObsType& ObsType) : ThObservable(ObsType) {};
    double getThValue(){
        
        return (sqrt((static_cast<const SUSY*> (&SM))->getMsu2()(1)));
    };
private:

};

#endif	/* MU2_H */

