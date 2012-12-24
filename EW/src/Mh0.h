/* 
 * Copyright (C) 2012 SusyFit Collaboration
 * All rights reserved.
 *
 * For the licensing terms see doc/COPYING.
 */

#ifndef MH0_H
#define	MH0_H

#include <ThObservable.h>
#include "EW.h"

class Mh0 : public ThObservable {
public: 
    Mh0(const EW& EW_i) : ThObservable(EW_i){ };
    double getThValue(){
        return SM.getMHl();
    }
private:

};

#endif	/* MH0_H */

