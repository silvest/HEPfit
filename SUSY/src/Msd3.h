/* 
 * Copyright (C) 2012 SusyFit Collaboration
 * All rights reserved.
 *
 * For the licensing terms see doc/COPYING.
 */

#ifndef MD3_H
#define	MD3_H

#include <ThObservable.h>
#include <ThObsType.h>
#include "SUSY.h"

class Msd3 : public ThObservable {
public:
    Msd3(const ThObsType& ObsType) : ThObservable(ObsType) {};
    double getThValue(){
        
        return (sqrt((static_cast<const SUSY*> (&SM))->getMsd2()(2)));
    };
private:

};

#endif	/* MD3_H */

