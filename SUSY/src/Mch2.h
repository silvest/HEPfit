/* 
 * Copyright (C) 2012 SusyFit Collaboration
 * All rights reserved.
 *
 * For the licensing terms see doc/COPYING.
 */

#ifndef MCH2_H
#define	MCH2_H

#include <ThObservable.h>
#include <ThObsType.h>
#include "SUSY.h"

#include <ThObservable.h>
#include <ThObsType.h>
#include "SUSY.h"

class Mch2 : public ThObservable {
public:
    Mch2(const ThObsType& ObsType) : ThObservable(ObsType) {};
    double getThValue(){
 
        return ((static_cast<const SUSY*> (&SM))->getMch()(1));
    };
private:

};

#endif	/* MCH2_H */

