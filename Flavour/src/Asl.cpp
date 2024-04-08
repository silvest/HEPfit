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

Asl_s_pole_NLO::Asl_s_pole_NLO(const StandardModel& SM_i)
: ThObservable(SM_i), AmpDB2(SM_i) {
}

Asl_s_pole_NLO::~Asl_s_pole_NLO() {
}

double Asl_s_pole_NLO::computeThValue() {
    return -getGamma21overM21_Bs(FULLNLO, pole).imag();
}

Asl_s_pole_LO::Asl_s_pole_LO(const StandardModel& SM_i)
: ThObservable(SM_i), AmpDB2(SM_i) {
}

Asl_s_pole_LO::~Asl_s_pole_LO() {
}

double Asl_s_pole_LO::computeThValue() {
    return -getGamma21overM21_Bs(LO, pole).imag();
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

Asl_s_MSbar_NLO::Asl_s_MSbar_NLO(const StandardModel& SM_i)
: ThObservable(SM_i), AmpDB2(SM_i) {
}

Asl_s_MSbar_NLO::~Asl_s_MSbar_NLO() {
}

double Asl_s_MSbar_NLO::computeThValue() {
    return -getGamma21overM21_Bs(FULLNLO, MSbar).imag();
}

Asl_s_MSbar_LO::Asl_s_MSbar_LO(const StandardModel& SM_i)
: ThObservable(SM_i), AmpDB2(SM_i) {
}

Asl_s_MSbar_LO::~Asl_s_MSbar_LO() {
}

double Asl_s_MSbar_LO::computeThValue() {
    return -getGamma21overM21_Bs(LO, MSbar).imag();
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

Asl_s_PS_NLO::Asl_s_PS_NLO(const StandardModel& SM_i)
: ThObservable(SM_i), AmpDB2(SM_i) {
}

Asl_s_PS_NLO::~Asl_s_PS_NLO() {
}

double Asl_s_PS_NLO::computeThValue() {
    return -getGamma21overM21_Bs(FULLNLO, PS).imag();
}

Asl_s_PS_LO::Asl_s_PS_LO(const StandardModel& SM_i)
: ThObservable(SM_i), AmpDB2(SM_i) {
}

Asl_s_PS_LO::~Asl_s_PS_LO() {
}

double Asl_s_PS_LO::computeThValue() {
    return -getGamma21overM21_Bs(LO, PS).imag();
}


Asl_s_pole_fixmub::Asl_s_pole_fixmub(const StandardModel& SM_i)
: ThObservable(SM_i), AmpDB2(SM_i, true) {
}

Asl_s_pole_fixmub::~Asl_s_pole_fixmub() {
}

double Asl_s_pole_fixmub::computeThValue() {
    return -getGamma21overM21_Bs(FULLNNLO, pole).imag();
}

Asl_s_MSbar_fixmub::Asl_s_MSbar_fixmub(const StandardModel& SM_i)
: ThObservable(SM_i), AmpDB2(SM_i, true) {
}

Asl_s_MSbar_fixmub::~Asl_s_MSbar_fixmub() {
}

double Asl_s_MSbar_fixmub::computeThValue() {
    return -getGamma21overM21_Bs(FULLNNLO, MSbar).imag();
}

Asl_s_PS_fixmub::Asl_s_PS_fixmub(const StandardModel& SM_i)
: ThObservable(SM_i), AmpDB2(SM_i, true) {
}

Asl_s_PS_fixmub::~Asl_s_PS_fixmub() {
}

double Asl_s_PS_fixmub::computeThValue() {
    return -getGamma21overM21_Bs(FULLNNLO, PS).imag();
}

Asl_d_only1overmb::Asl_d_only1overmb(const StandardModel& SM_i)
: ThObservable(SM_i), AmpDB2(SM_i) {
}

Asl_d_only1overmb::~Asl_d_only1overmb() {
}

double Asl_d_only1overmb::computeThValue() {
    return -getGamma21overM21_Bd(FULLNNLO, only1overmb).imag();
}

Asl_s_only1overmb::Asl_s_only1overmb(const StandardModel& SM_i)
: ThObservable(SM_i), AmpDB2(SM_i) {
}

Asl_s_only1overmb::~Asl_s_only1overmb() {
}

double Asl_s_only1overmb::computeThValue() {
    return -getGamma21overM21_Bs(FULLNNLO, only1overmb).imag();
}

Asl_d_MSbar_NLO_tradBasis::Asl_d_MSbar_NLO_tradBasis(const StandardModel& SM_i)
: ThObservable(SM_i), AmpDB2(SM_i) {
}

Asl_d_MSbar_NLO_tradBasis::~Asl_d_MSbar_NLO_tradBasis() {
}

double Asl_d_MSbar_NLO_tradBasis::computeThValue() {
    return -getGamma21overM21_Bd_tradBasis(FULLNLO).imag();
}

Asl_s_MSbar_NLO_tradBasis::Asl_s_MSbar_NLO_tradBasis(const StandardModel& SM_i)
: ThObservable(SM_i), AmpDB2(SM_i) {
}

Asl_s_MSbar_NLO_tradBasis::~Asl_s_MSbar_NLO_tradBasis() {
}

double Asl_s_MSbar_NLO_tradBasis::computeThValue() {
    return -getGamma21overM21_Bs_tradBasis(FULLNLO).imag();
}

Asl_s_MSbar_NLO_RI::Asl_s_MSbar_NLO_RI(const StandardModel& SM_i)
: ThObservable(SM_i), AmpDB2(SM_i, false, true) {
}

Asl_s_MSbar_NLO_RI::~Asl_s_MSbar_NLO_RI() {
}

double Asl_s_MSbar_NLO_RI::computeThValue() {
    return -getGamma21overM21_Bs(FULLNLO, MSbar).imag();
}

Asl_s_MSbar_NLO_RI_tradBasis::Asl_s_MSbar_NLO_RI_tradBasis(const StandardModel& SM_i)
: ThObservable(SM_i), AmpDB2(SM_i, false, true) {
}

Asl_s_MSbar_NLO_RI_tradBasis::~Asl_s_MSbar_NLO_RI_tradBasis() {
}

double Asl_s_MSbar_NLO_RI_tradBasis::computeThValue() {
    return -getGamma21overM21_Bs_tradBasis(FULLNLO).imag();
}

Asl_s_PS_NLO_RI::Asl_s_PS_NLO_RI(const StandardModel& SM_i)
: ThObservable(SM_i), AmpDB2(SM_i, false, true) {
}

Asl_s_PS_NLO_RI::~Asl_s_PS_NLO_RI() {
}

double Asl_s_PS_NLO_RI::computeThValue() {
    return -getGamma21overM21_Bs(FULLNLO, PS).imag();
}

Asl_s_PS_NLO_RI_tradBasis::Asl_s_PS_NLO_RI_tradBasis(const StandardModel& SM_i)
: ThObservable(SM_i), AmpDB2(SM_i, false, true) {
}

Asl_s_PS_NLO_RI_tradBasis::~Asl_s_PS_NLO_RI_tradBasis() {
}

double Asl_s_PS_NLO_RI_tradBasis::computeThValue() {
    return -getGamma21overM21_Bs_tradBasis(FULLNLO).imag();
}

Asl_s_MSbar_partialNNLO::Asl_s_MSbar_partialNNLO(const StandardModel& SM_i)
: ThObservable(SM_i), AmpDB2(SM_i) {
}

Asl_s_MSbar_partialNNLO::~Asl_s_MSbar_partialNNLO() {
}

double Asl_s_MSbar_partialNNLO::computeThValue() {
    return -getGamma21overM21_Bs(FULLNNLO, MSbar_partialNNLO).imag();
}

Asl_s_PS_partialNNLO::Asl_s_PS_partialNNLO(const StandardModel& SM_i)
: ThObservable(SM_i), AmpDB2(SM_i) {
}

Asl_s_PS_partialNNLO::~Asl_s_PS_partialNNLO() {
}

double Asl_s_PS_partialNNLO::computeThValue() {
    return -getGamma21overM21_Bs(FULLNNLO, PS_partialNNLO).imag();
}

Asl_s_MSbar_partialN3LO::Asl_s_MSbar_partialN3LO(const StandardModel& SM_i)
: ThObservable(SM_i), AmpDB2(SM_i) {
}

Asl_s_MSbar_partialN3LO::~Asl_s_MSbar_partialN3LO() {
}

double Asl_s_MSbar_partialN3LO::computeThValue() {
    return -getGamma21overM21_Bs(FULLNNNLO, MSbar_partialN3LO).imag();
}

Asl_s_PS_partialN3LO::Asl_s_PS_partialN3LO(const StandardModel& SM_i)
: ThObservable(SM_i), AmpDB2(SM_i) {
}

Asl_s_PS_partialN3LO::~Asl_s_PS_partialN3LO() {
}

double Asl_s_PS_partialN3LO::computeThValue() {
    return -getGamma21overM21_Bs(FULLNNNLO, PS_partialN3LO).imag();
}