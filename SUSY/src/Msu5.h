/* 
 * Copyright (C) 2012 SusyFit Collaboration
 * All rights reserved.
 *
 * For the licensing terms see doc/COPYING.
 */

#ifndef MU5_H
#define	MU5_H

#include <ThObservable.h>
#include <ThObsType.h>
#include "SUSY.h"

class Msu5 : public ThObservable {
public:
    Msu5(const ThObsType& ObsType) : ThObservable(ObsType) {};
    double getThValue(){
        
        return (sqrt((static_cast<const SUSY*> (&SM))->getMsu2()(4)));
    };
private:

};

#endif	/* MU5_H */

