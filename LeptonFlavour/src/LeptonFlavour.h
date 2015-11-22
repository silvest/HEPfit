/*
 * Copyright (C) 2013 HEPfit Collaboration
 * All rights reserved.
 *
 * For the licensing terms see doc/COPYING.
 */

#ifndef LEPTONFLAVOUR_H
#define	LEPTONFLAVOUR_H

#include "StandardModel.h"
#include "HeffDLij.h"
#include "HeffDLi3j.h"
#include "Heffmueconv.h"
#include "Heffgminus2.h"

class LeptonFlavour {
public:

    LeptonFlavour(const StandardModel& SM_i) :
            HDLij(SM_i),HDLi3j(SM_i),Hmueconv(SM_i),Hgminus2(SM_i)
    {};

    const HeffDLij& getHDLij() const {
        return HDLij;
    }

    gslpp::vector<gslpp::complex>** ComputeCoeffli_lj_gamma(int li_lj) {
        return HDLij.ComputeCoeffDLij(li_lj);
    }

    const HeffDLi3j& getHDLi3j() const {
        return HDLi3j;
    }

    gslpp::vector<gslpp::complex>** ComputeCoeffli_3lj(int li_lj) {
        return HDLi3j.ComputeCoeffDLi3j(li_lj);
    }

    const Heffmueconv& getHmueconv() const {
        return Hmueconv;
    }

    gslpp::vector<gslpp::complex>** ComputeCoeffmueconversion() {
        return Hmueconv.ComputeCoeffmueconv();
    }

    const Heffgminus2& getHgminus2() const {
        return Hgminus2;
    }

    gslpp::vector<gslpp::complex>** ComputeCoeffgminus2mu() {
        return Hgminus2.ComputeCoeffgm2mu();
    }

private:
    HeffDLij HDLij;
    HeffDLi3j HDLi3j;
    Heffmueconv Hmueconv;
    Heffgminus2 Hgminus2;
};

#endif	/* LEPTONFLAVOUR_H */
