/*
 * Copyright (C) 2013 HEPfit Collaboration
 * All rights reserved.
 *
 * For the licensing terms see doc/COPYING.
 */

#ifndef LEPTONFLAVOUR_H
#define	LEPTONFLAVOUR_H

#include <StandardModel.h>
#include "HeffDL1.h"

using namespace gslpp;

class LeptonFlavour {
public:

    LeptonFlavour(const StandardModel& SM_i) :
            HDL1(SM_i)//, HDS1(SM_i), HDB1(SM_i)
    {};

    const HeffDL1& getHDL1() const {
        return HDL1;
    }
    
    vector<complex>** ComputeCoeffli_lj_gamma() {
        return HDL1.ComputeCoeffDL1();
    }
    
private:
    HeffDL1 HDL1;
};

#endif	/* LEPTONFLAVOUR_H */

