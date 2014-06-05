/* 
 * File:   ParamObs.cpp
 * Author: marco
 * 
 * Created on May 16, 2014, 5:48 PM
 */

#include "ParamObs.h"

ParamObs::ParamObs(const StandardModel& SM_i, std::string name) 
: ThObservable(SM_i), param(SM_i.getModelParam(name))
{
}

double ParamObs::computeThValue()
{
    return param;
}
