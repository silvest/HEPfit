/* 
 * Copyright (C) 2012 SusyFit Collaboration
 * All rights reserved.
 *
 * For the licensing terms see doc/COPYING.
 */

#ifndef VUB_H
#define	VUB_H

#include <ThObservable.h>

class Vub : public ThObservable {
public:
    Vub(const StandardModel& SM_i) : ThObservable(SM_i) {};

    double computeThValue();
};

#endif	/* VUB_H */
