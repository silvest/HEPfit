/* 
 * File:   ThObservable.h
 * Author: silvest
 *
 * Created on March 29, 2011, 12:15 PM
 */

#ifndef THOBSERVABLE_H
#define	THOBSERVABLE_H


#include <StandardModel.h>
#include "ThObsType.h"

class ThObservable {
public:
    ThObservable(const ThObsType& ObsType_i);
    ThObservable(const StandardModel& SM_i);
    ThObservable(const ThObservable& orig);
    virtual ~ThObservable();
    virtual double getThValue() = 0;

    static double const GeVminus2_to_nb;// conversion factor from GeV^{-2} to nb

protected:
    const ThObsType& ObsType;
    const StandardModel& SM;
};

#endif	/* THOBSERVABLE_H */
