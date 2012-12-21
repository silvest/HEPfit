/* 
 * Copyright (C) 2012 SUSYfit Collaboration
 * All rights reserved.
 *
 * For the licensing terms see doc/COPYING.
 */

#ifndef MD4_H
#define	MD4_H

#include <ThObservable.h>
#include <ThObsType.h>
#include "SUSY.h"

class Msd4 : public ThObservable {
public:
    Msd4(const ThObsType& ObsType) : ThObservable(ObsType) {};
    double getThValue(){
        
        return (sqrt((static_cast<const SUSY*> (&SM))->getMsd2()(3)));
    };
private:

};

#endif	/* MD4_H */

