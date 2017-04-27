/*
 * Copyright (C) 2012 SusyFit Collaboration
 *
 *
 * For the licensing terms see doc/COPYING.
 */

#include "alpha_s.h"
#include "StandardModel.h"

alpha_s::alpha_s(const StandardModel& SM_i, orders order)
: ThObservable(SM_i)
{
    if (order >= LO and order <= FULLNNLO) this->order = order;
    else throw std::runtime_error("orders can only be LO, NLO, FULLNLO, NNLO or FULLNNLO");
};

double alpha_s::computeThValue()
{
    return SM.Als(getBinMin(), order);   
}