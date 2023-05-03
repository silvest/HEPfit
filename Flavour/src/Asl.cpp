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
    return Gamma12overM12_Bd(FULLNNLO).imag();
}

Asl_s::Asl_s(const StandardModel& SM_i)
: ThObservable(SM_i), AmpDB2(SM_i) {
}

Asl_s::~Asl_s() {
}

double Asl_s::computeThValue() {
    return Gamma12overM12_Bs(FULLNNLO).imag();
}

Asl_d_NLO::Asl_d_NLO(const StandardModel& SM_i)
: ThObservable(SM_i), AmpDB2(SM_i) {
}

Asl_d_NLO::~Asl_d_NLO() {
}

double Asl_d_NLO::computeThValue() {
    return Gamma12overM12_Bd(FULLNLO).imag();
}

Asl_s_NLO::Asl_s_NLO(const StandardModel& SM_i)
: ThObservable(SM_i), AmpDB2(SM_i) {
}

Asl_s_NLO::~Asl_s_NLO() {
}

double Asl_s_NLO::computeThValue() {
    return Gamma12overM12_Bs(FULLNLO).imag();
}

Asl_d_NLO1::Asl_d_NLO1(const StandardModel& SM_i)
: ThObservable(SM_i), AmpDB2(SM_i) {
}

Asl_d_NLO1::~Asl_d_NLO1() {
}

double Asl_d_NLO1::computeThValue() {
    return Gamma12overM12_BdFULLNLO1().imag();
}

Asl_s_NLO1::Asl_s_NLO1(const StandardModel& SM_i)
: ThObservable(SM_i), AmpDB2(SM_i) {
}

Asl_s_NLO1::~Asl_s_NLO1() {
}

double Asl_s_NLO1::computeThValue() {
    return Gamma12overM12_BsFULLNLO1().imag();
}
