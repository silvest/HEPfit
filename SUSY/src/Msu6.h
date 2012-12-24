/* 
 * Copyright (C) 2012 SusyFit Collaboration
 * All rights reserved.
 *
 * For the licensing terms see doc/COPYING.
 */

#ifndef MU6_H
#define	MU6_H

#include <ThObservable.h>
#include <ThObsType.h>
#include "SUSY.h"

class Msu6 : public ThObservable {
public:
    Msu6(const ThObsType& ObsType) : ThObservable(ObsType) {};
    double getThValue(){
        
        return (sqrt((static_cast<const SUSY*> (&SM))->getMsu2()(5)));
    };
private:

};

#endif	/* MU6_H */

