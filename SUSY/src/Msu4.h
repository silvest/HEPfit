/* 
 * Copyright (C) 2012 SusyFit Collaboration
 * All rights reserved.
 *
 * For the licensing terms see doc/COPYING.
 */

#ifndef MU4_H
#define	MU4_H

#include <ThObservable.h>
#include <ThObsType.h>
#include "SUSY.h"

class Msu4 : public ThObservable {
public:
    Msu4(const ThObsType& ObsType) : ThObservable(ObsType) {};
    double getThValue(){
        
        return (sqrt((static_cast<const SUSY*> (&SM))->getMsu2()(3)));
    };
private:

};

#endif	/* MU4_H */

