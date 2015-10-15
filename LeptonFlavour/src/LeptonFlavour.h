/*
 * Copyright (C) 2013 HEPfit Collaboration
 * All rights reserved.
 *
 * For the licensing terms see doc/COPYING.
 */

#ifndef LEPTONFLAVOUR_H
#define	LEPTONFLAVOUR_H

#include <StandardModel.h>
#include "HeffDLij.h"
#include "HeffDLi3j.h"

using namespace gslpp;

class LeptonFlavour {
public:

    LeptonFlavour(const StandardModel& SM_i) :
            HDLij(SM_i),HDLi3j(SM_i)//, HDS1(SM_i), HDB1(SM_i)
    {};

    const HeffDLij& getHDLij() const {
        return HDLij;
    }

    vector<complex>** ComputeCoeffli_lj_gamma(int li_lj) {
        return HDLij.ComputeCoeffDLij(li_lj);
    }

    const HeffDLi3j& getHDLi3j() const {
        return HDLi3j;
    }

    vector<complex>** ComputeCoeffli_3lj(int li_lj) {
        return HDLi3j.ComputeCoeffDLi3j(li_lj);
    }

private:
    HeffDLij HDLij;
    HeffDLi3j HDLi3j;
};

#endif	/* LEPTONFLAVOUR_H */

