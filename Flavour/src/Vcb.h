/* 
 * Copyright (C) 2012 SusyFit Collaboration
 * All rights reserved.
 *
 * For the licensing terms see doc/COPYING.
 */

#ifndef VCB_H
#define	VCB_H

#include <ThObservable.h>

class Vcb : public ThObservable {
public:
    Vcb(const StandardModel& SM_i) : ThObservable(SM_i) {};

    double computeThValue();
};

#endif	/* VCB_H */
