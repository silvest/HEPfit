/* 
 * Copyright (C) 2012 HEPfit Collaboration
 * All rights reserved.
 *
 * For the licensing terms see doc/COPYING.
 */

#ifndef ALPHA_2A_H
#define	ALPHA_2A_H

#include <ThObservable.h>



class Alpha_2a : public ThObservable {
public:
    Alpha_2a(const StandardModel& SM_i) : ThObservable(SM_i) {};

    double computeThValue();
};

#endif	/* ALPHA_2A_H */

