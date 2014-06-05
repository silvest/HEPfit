/* 
 * File:   ParamObs.h
 * Author: marco
 *
 * Created on May 16, 2014, 5:48 PM
 */

#ifndef PARAMOBS_H
#define	PARAMOBS_H

#include "ThObservable.h"

class ParamObs : public ThObservable {
public:
    ParamObs(const StandardModel& SM, const std::string name);

    double computeThValue();
    
private:
    const double& param;
};

#endif	/* PARAMOBS_H */

