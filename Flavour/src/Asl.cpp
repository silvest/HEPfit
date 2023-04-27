/* 
 * Copyright (C) 2023 HEPfit Collaboration
 *
 *
 * For the licensing terms see doc/COPYING.
 */

#include "Asl.h"
#include "StandardModel.h"
#include "AmpDB2.h"
#include "std_make_vector.h"

Asl_d::Asl_d(const StandardModel& SM_i)
: ThObservable(SM_i), AmpDB2(SM_i) {
}

Asl_d::~Asl_d() {
}

double Asl_d::computeThValue() {
    return Gamma12oM12_Bd(FULLNLO).imag();
}

Asl_s::Asl_s(const StandardModel& SM_i)
: ThObservable(SM_i), AmpDB2(SM_i) {
}

Asl_s::~Asl_s() {
}

double Asl_s::computeThValue() {
    return Gamma12oM12_Bs(FULLNLO).imag();
}