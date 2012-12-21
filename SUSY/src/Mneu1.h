/* 
 * Copyright (C) 2012 SUSYfit Collaboration
 * All rights reserved.
 *
 * For the licensing terms see doc/COPYING.
 */

#ifndef MNEU1_H
#define	MNEU1_H

#include <ThObservable.h>
#include <ThObsType.h>
#include "SUSY.h"

class Mneu1 : public ThObservable {
public:
    Mneu1(const ThObsType& ObsType) : ThObservable(ObsType) {};
    double getThValue(){
        
        return ((static_cast<const SUSY*> (&SM))->getMneu()(0));
    };
private:

};

#endif	/* MNEU1_H */

