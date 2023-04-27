/*
 * Copyright (C) 2023 HEPfit Collaboration
 * 
 * 
 * For the licensing terms see doc/COPYING.
 */


#include "DGamma.h"
#include "StandardModel.h"
#include "AmpDB2.h"
#include "std_make_vector.h"

DGamma_d::DGamma_d(const StandardModel& SM_i)
: ThObservable(SM_i), AmpDB2(SM_i) {
    setParametersForObservable(make_vector<std::string>() << "DmBd_ex" );
}

DGamma_d::~DGamma_d() {
}

double DGamma_d::computeThValue() {
    return -Gamma12oM12_Bd(FULLNLO).real() * SM.getOptionalParameter("DmBd_ex");
}

DGamma_s::DGamma_s(const StandardModel& SM_i)
: ThObservable(SM_i), AmpDB2(SM_i) {
    setParametersForObservable(make_vector<std::string>() << "DmBs_ex" );
}

DGamma_s::~DGamma_s() {
}

double DGamma_s::computeThValue() {
    return -Gamma12oM12_Bs(FULLNLO).real() * SM.getOptionalParameter("DmBs_ex");
}

