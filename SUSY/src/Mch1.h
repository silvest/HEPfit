/* 
 * Copyright (C) 2012 SusyFit Collaboration
 * All rights reserved.
 *
 * For the licensing terms see doc/COPYING.
 */

#ifndef MCH1_H
#define	MCH1_H

#include <ThObservable.h>
#include <ThObsType.h>
#include "SUSY.h"

class Mch1 : public ThObservable {
public:
    Mch1(const ThObsType& ObsType) : ThObservable(ObsType) {};
    double getThValue(){        
 
        return ((static_cast<const SUSY*> (&SM))->getMch()(0));
    };
private:

};


#endif	/* MCH1_H */

