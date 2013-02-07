/* 
 * Copyright (C) 2012 SusyFit Collaboration
 * All rights reserved.
 *
 * For the licensing terms see doc/COPYING.
 */

#ifndef MU1_H
#define	MU1_H

#include <ThObservable.h>
#include <ThObsType.h>
#include "SUSY.h"

class Msu1 : public ThObservable {
public:
    Msu1(const ThObsType& ObsType) : ThObservable(ObsType) {};
    double getThValue(){
        
        return (sqrt((static_cast<const SUSY*> (&SM))->getMsu2()(0)));
    };
private:

};


#endif	/* MU1_H */

