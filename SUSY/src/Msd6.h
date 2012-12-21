/* 
 * Copyright (C) 2012 SUSYfit Collaboration
 * All rights reserved.
 *
 * For the licensing terms see doc/COPYING.
 */

#ifndef MD6_H
#define	MD6_H

#include <ThObservable.h>
#include <ThObsType.h>
#include "SUSY.h"

class Msd6 : public ThObservable {
public:
    Msd6(const ThObsType& ObsType) : ThObservable(ObsType) {};
    double getThValue(){
        
        return (sqrt((static_cast<const SUSY*> (&SM))->getMsd2()(5)));
    };
private:

};

#endif	/* MD6_H */

