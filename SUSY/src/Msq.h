/* 
 * Copyright (C) 2012 SusyFit Collaboration
 * All rights reserved.
 *
 * For the licensing terms see doc/COPYING.
 */

#ifndef MSQ_H
#define	MSQ_H

#include <ThObservable.h>
#include <ThObsType.h>
#include "SUSY.h"

class Msq : public ThObservable {
public:
    Msq(const ThObsType& ObsType) : ThObservable(ObsType) {};
    double getThValue(){
        double res = (static_cast<const SUSY*> (&SM))->getMsd2()(0);
        return ((static_cast<const SUSY*> (&SM))->getMsd2()(0));
    };
private:

};

#endif	/* MSQ_H */

