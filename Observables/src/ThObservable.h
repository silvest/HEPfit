/* 
 * Copyright (C) 2012 SusyFit Collaboration
 * All rights reserved.
 *
 * For the licensing terms see doc/COPYING.
 */

#ifndef THOBSERVABLE_H
#define	THOBSERVABLE_H


#include <StandardModel.h>
#include "ThObsType.h"

class ThObservable {
public:
    ThObservable(const ThObsType& ObsType_i);
    ThObservable(const ThObservable& orig);
    virtual ~ThObservable();
    virtual double getThValue() = 0;

    static const double GeVminus2_to_nb;// conversion factor from GeV^{-2} to nb

protected:
    const ThObsType& ObsType;
    const StandardModel& SM;
};

#endif	/* THOBSERVABLE_H */
