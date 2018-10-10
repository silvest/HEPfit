/*
 * Copyright (C) 2016 HEPfit Collaboration
 *
 * For the licensing terms see doc/COPYING.
 */

#ifndef BDMUMUOBSMUMU_H
#define BDMUMUOBSMUMU_H

class StandardModel;
#include "ThObservable.h"

class BdmumuOBsmumu : public ThObservable {
public:
    BdmumuOBsmumu(const StandardModel& SM_i);

    double computeThValue();
};

#endif /* BDMUMUOBSMUMU_H */

