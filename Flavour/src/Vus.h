/* 
 * Copyright (C) 2012 SusyFit Collaboration
 * All rights reserved.
 *
 * For the licensing terms see doc/COPYING.
 */

#ifndef VUS_H
#define	VUS_H

#include <ThObservable.h>

class Vus : public ThObservable {
public:
    Vus(const StandardModel& SM_i) : ThObservable(SM_i) {};

    double computeThValue();
};

#endif	/* VUS_H */
