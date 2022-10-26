/* 
 * Copyright (C) 2022 HEPfit Collaboration
 *
 *
 * For the licensing terms see doc/COPYING.
 */

#include "BBd.h"
#include "StandardModel.h"

BBd::BBd(const StandardModel& SM_i) : ThObservable(SM_i) 
{
  SM.initializeBParameter("BBd");
};

double BBd::computeThValue() 
{
  return(SM.getBBd().getBpars()(0));
}
