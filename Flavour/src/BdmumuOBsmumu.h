/*
 * Copyright (C) 2016 SusyFit Collaboration
 *
 * For the licensing terms see doc/COPYING.
 */

#ifndef BDMUMUOBSMUMU_H
#define BDMUMUOBSMUMU_H

#include "ThObservable.h"
#include "Flavour.h"
#include "StandardModel.h"

class BdmumuOBsmumu : public ThObservable {
public:
    BdmumuOBsmumu(const StandardModel& SM_i);

    double computeThValue();
};

#endif /* BDMUMUOBSMUMU_H */

