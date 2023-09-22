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


Asl_d_pole::Asl_d_pole(const StandardModel& SM_i)
: ThObservable(SM_i), AmpDB2(SM_i) {
}

Asl_d_pole::~Asl_d_pole() {
}

double Asl_d_pole::computeThValue() {
    return -getGamma21overM21_Bd(FULLNNLO, pole).imag();
}

Asl_s_pole::Asl_s_pole(const StandardModel& SM_i)
: ThObservable(SM_i), AmpDB2(SM_i) {
}

Asl_s_pole::~Asl_s_pole() {
}

double Asl_s_pole::computeThValue() {
    return -getGamma21overM21_Bs(FULLNNLO, pole).imag();
}

Asl_d_MSbar::Asl_d_MSbar(const StandardModel& SM_i)
: ThObservable(SM_i), AmpDB2(SM_i) {
}

Asl_d_MSbar::~Asl_d_MSbar() {
}

double Asl_d_MSbar::computeThValue() {
    return -getGamma21overM21_Bd(FULLNNLO, MSbar).imag();
}

Asl_s_MSbar::Asl_s_MSbar(const StandardModel& SM_i)
: ThObservable(SM_i), AmpDB2(SM_i) {
}

Asl_s_MSbar::~Asl_s_MSbar() {
}

double Asl_s_MSbar::computeThValue() {
    return -getGamma21overM21_Bs(FULLNNLO, MSbar).imag();
}

Asl_d_PS::Asl_d_PS(const StandardModel& SM_i)
: ThObservable(SM_i), AmpDB2(SM_i) {
}

Asl_d_PS::~Asl_d_PS() {
}

double Asl_d_PS::computeThValue() {
    return -getGamma21overM21_Bd(FULLNNLO, PS).imag();
}

Asl_s_PS::Asl_s_PS(const StandardModel& SM_i)
: ThObservable(SM_i), AmpDB2(SM_i) {
}

Asl_s_PS::~Asl_s_PS() {
}

double Asl_s_PS::computeThValue() {
    return -getGamma21overM21_Bs(FULLNNLO, PS).imag();
}

Asl_d_MSbar_NLO::Asl_d_MSbar_NLO(const StandardModel& SM_i)
: ThObservable(SM_i), AmpDB2(SM_i) {
}

Asl_d_MSbar_NLO::~Asl_d_MSbar_NLO() {
}

double Asl_d_MSbar_NLO::computeThValue() {
    return -getGamma21overM21_Bd(FULLNLO, MSbar).imag();
}

Asl_s_MSbar_NLO::Asl_s_MSbar_NLO(const StandardModel& SM_i)
: ThObservable(SM_i), AmpDB2(SM_i) {
}

Asl_s_MSbar_NLO::~Asl_s_MSbar_NLO() {
}

double Asl_s_MSbar_NLO::computeThValue() {
    return -getGamma21overM21_Bs(FULLNLO, MSbar).imag();
}

Asl_d_MSbar_NLO_tradBasis::Asl_d_MSbar_NLO_tradBasis(const StandardModel& SM_i)
: ThObservable(SM_i), AmpDB2(SM_i, true) {
}

Asl_d_MSbar_NLO_tradBasis::~Asl_d_MSbar_NLO_tradBasis() {
}

double Asl_d_MSbar_NLO_tradBasis::computeThValue() {
    return -Gamma21overM21_BdFULLNLO_tradBasis().imag();
}

Asl_s_MSbar_NLO_tradBasis::Asl_s_MSbar_NLO_tradBasis(const StandardModel& SM_i)
: ThObservable(SM_i), AmpDB2(SM_i, true) {
}

Asl_s_MSbar_NLO_tradBasis::~Asl_s_MSbar_NLO_tradBasis() {
}

double Asl_s_MSbar_NLO_tradBasis::computeThValue() {
    return -Gamma21overM21_BsFULLNLO_tradBasis().imag();
}
