/* 
 * Copyright (C) 2012 SusyFit Collaboration
 * All rights reserved.
 *
 * For the licensing terms see doc/COPYING.
 */

#ifndef MD1_H
#define	MD1_H

#include <ThObservable.h>
#include <ThObsType.h>
#include "SUSY.h"

class Msd1 : public ThObservable {
public:
    Msd1(const ThObsType& ObsType) : ThObservable(ObsType) {};
    double getThValue(){
        
        return (sqrt((static_cast<const SUSY*> (&SM))->getMsd2()(0)));
    };
private:

};


#endif	/* MD1_H */

