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


DGamma_d_pole::DGamma_d_pole(const StandardModel& SM_i)
: ThObservable(SM_i), AmpDB2(SM_i) {
    setParametersForObservable(make_vector<std::string>() << "DmBd_ex" );
}

DGamma_d_pole::~DGamma_d_pole() {
}

double DGamma_d_pole::computeThValue() {
    return -getGamma21overM21_Bd(FULLNNLO, pole).real() * SM.getOptionalParameter("DmBd_ex");
}

DGamma_s_pole::DGamma_s_pole(const StandardModel& SM_i)
: ThObservable(SM_i), AmpDB2(SM_i) {
    setParametersForObservable(make_vector<std::string>() << "DmBs_ex" );
}

DGamma_s_pole::~DGamma_s_pole() {
}

double DGamma_s_pole::computeThValue() {
    return -getGamma21overM21_Bs(FULLNNLO, pole).real() * SM.getOptionalParameter("DmBs_ex");
}

DGamma_s_pole_NLO::DGamma_s_pole_NLO(const StandardModel& SM_i)
: ThObservable(SM_i), AmpDB2(SM_i) {
    setParametersForObservable(make_vector<std::string>() << "DmBs_ex" );
}

DGamma_s_pole_NLO::~DGamma_s_pole_NLO() {
}

double DGamma_s_pole_NLO::computeThValue() {
    return -getGamma21overM21_Bs(FULLNLO, pole).real() * SM.getOptionalParameter("DmBs_ex");
}

DGamma_s_pole_LO::DGamma_s_pole_LO(const StandardModel& SM_i)
: ThObservable(SM_i), AmpDB2(SM_i) {
    setParametersForObservable(make_vector<std::string>() << "DmBs_ex" );
}

DGamma_s_pole_LO::~DGamma_s_pole_LO() {
}

double DGamma_s_pole_LO::computeThValue() {
    return -getGamma21overM21_Bs(LO, pole).real() * SM.getOptionalParameter("DmBs_ex");
}

DGamma_d_MSbar::DGamma_d_MSbar(const StandardModel& SM_i)
: ThObservable(SM_i), AmpDB2(SM_i) {
    setParametersForObservable(make_vector<std::string>() << "DmBd_ex" );
}

DGamma_d_MSbar::~DGamma_d_MSbar() {
}

double DGamma_d_MSbar::computeThValue() {
    return -getGamma21overM21_Bd(FULLNNLO, MSbar).real() * SM.getOptionalParameter("DmBd_ex");
}

DGamma_s_MSbar::DGamma_s_MSbar(const StandardModel& SM_i)
: ThObservable(SM_i), AmpDB2(SM_i) {
    setParametersForObservable(make_vector<std::string>() << "DmBs_ex" );
}

DGamma_s_MSbar::~DGamma_s_MSbar() {
}

double DGamma_s_MSbar::computeThValue() {
    return -getGamma21overM21_Bs(FULLNNLO, MSbar).real() * SM.getOptionalParameter("DmBs_ex");
}

DGamma_d_PS::DGamma_d_PS(const StandardModel& SM_i)
: ThObservable(SM_i), AmpDB2(SM_i) {
    setParametersForObservable(make_vector<std::string>() << "DmBd_ex" );
}

DGamma_d_PS::~DGamma_d_PS() {
}

double DGamma_d_PS::computeThValue() {
    return -getGamma21overM21_Bd(FULLNNLO, PS).real() * SM.getOptionalParameter("DmBd_ex");
}

DGamma_s_PS::DGamma_s_PS(const StandardModel& SM_i)
: ThObservable(SM_i), AmpDB2(SM_i) {
    setParametersForObservable(make_vector<std::string>() << "DmBs_ex" );
}

DGamma_s_PS::~DGamma_s_PS() {
}

double DGamma_s_PS::computeThValue() {
    return -getGamma21overM21_Bs(FULLNNLO, PS).real() * SM.getOptionalParameter("DmBs_ex");
}

////////////

DGamma_d_MSbar_NLO::DGamma_d_MSbar_NLO(const StandardModel& SM_i)
: ThObservable(SM_i), AmpDB2(SM_i) {
    setParametersForObservable(make_vector<std::string>() << "DmBd_ex" );
}

DGamma_d_MSbar_NLO::~DGamma_d_MSbar_NLO() {
}

double DGamma_d_MSbar_NLO::computeThValue() {
    return -getGamma21overM21_Bd(FULLNLO, MSbar).real() * SM.getOptionalParameter("DmBd_ex");
}

DGamma_s_MSbar_NLO::DGamma_s_MSbar_NLO(const StandardModel& SM_i)
: ThObservable(SM_i), AmpDB2(SM_i) {
    setParametersForObservable(make_vector<std::string>() << "DmBs_ex" );
}

DGamma_s_MSbar_NLO::~DGamma_s_MSbar_NLO() {
}

double DGamma_s_MSbar_NLO::computeThValue() {
    return -getGamma21overM21_Bs(FULLNLO, MSbar).real() * SM.getOptionalParameter("DmBs_ex");
}

DGamma_d_PS_NLO::DGamma_d_PS_NLO(const StandardModel& SM_i)
: ThObservable(SM_i), AmpDB2(SM_i) {
    setParametersForObservable(make_vector<std::string>() << "DmBd_ex" );
}

DGamma_d_PS_NLO::~DGamma_d_PS_NLO() {
}

double DGamma_d_PS_NLO::computeThValue() {
    return -getGamma21overM21_Bd(FULLNLO, PS).real() * SM.getOptionalParameter("DmBd_ex");
}

DGamma_s_PS_NLO::DGamma_s_PS_NLO(const StandardModel& SM_i)
: ThObservable(SM_i), AmpDB2(SM_i) {
    setParametersForObservable(make_vector<std::string>() << "DmBs_ex" );
}

DGamma_s_PS_NLO::~DGamma_s_PS_NLO() {
}

double DGamma_s_PS_NLO::computeThValue() {
    return -getGamma21overM21_Bs(FULLNLO, PS).real() * SM.getOptionalParameter("DmBs_ex");
}


DGamma_d_MSbar_LO::DGamma_d_MSbar_LO(const StandardModel& SM_i)
: ThObservable(SM_i), AmpDB2(SM_i) {
    setParametersForObservable(make_vector<std::string>() << "DmBd_ex" );
}

DGamma_d_MSbar_LO::~DGamma_d_MSbar_LO() {
}

double DGamma_d_MSbar_LO::computeThValue() {
    return -getGamma21overM21_Bd(LO, MSbar).real() * SM.getOptionalParameter("DmBd_ex");
}

DGamma_s_MSbar_LO::DGamma_s_MSbar_LO(const StandardModel& SM_i)
: ThObservable(SM_i), AmpDB2(SM_i) {
    setParametersForObservable(make_vector<std::string>() << "DmBs_ex" );
}

DGamma_s_MSbar_LO::~DGamma_s_MSbar_LO() {
}

double DGamma_s_MSbar_LO::computeThValue() {
    return -getGamma21overM21_Bs(LO, MSbar).real() * SM.getOptionalParameter("DmBs_ex");
}

DGamma_d_PS_LO::DGamma_d_PS_LO(const StandardModel& SM_i)
: ThObservable(SM_i), AmpDB2(SM_i) {
    setParametersForObservable(make_vector<std::string>() << "DmBd_ex" );
}

DGamma_d_PS_LO::~DGamma_d_PS_LO() {
}

double DGamma_d_PS_LO::computeThValue() {
    return -getGamma21overM21_Bd(LO, PS).real() * SM.getOptionalParameter("DmBd_ex");
}

DGamma_s_PS_LO::DGamma_s_PS_LO(const StandardModel& SM_i)
: ThObservable(SM_i), AmpDB2(SM_i) {
    setParametersForObservable(make_vector<std::string>() << "DmBs_ex" );
}

DGamma_s_PS_LO::~DGamma_s_PS_LO() {
}

double DGamma_s_PS_LO::computeThValue() {
    return -getGamma21overM21_Bs(LO, PS).real() * SM.getOptionalParameter("DmBs_ex");
}

DGamma_d_only1overmb::DGamma_d_only1overmb(const StandardModel& SM_i)
: ThObservable(SM_i), AmpDB2(SM_i) {
    setParametersForObservable(make_vector<std::string>() << "DmBd_ex" );
}

DGamma_d_only1overmb::~DGamma_d_only1overmb() {
}

double DGamma_d_only1overmb::computeThValue() {
    return -getGamma21overM21_Bd(FULLNNLO, only1overmb).real() * SM.getOptionalParameter("DmBd_ex");
}

DGamma_s_only1overmb::DGamma_s_only1overmb(const StandardModel& SM_i)
: ThObservable(SM_i), AmpDB2(SM_i) {
    setParametersForObservable(make_vector<std::string>() << "DmBs_ex" );
}

DGamma_s_only1overmb::~DGamma_s_only1overmb() {
}

double DGamma_s_only1overmb::computeThValue() {
    return -getGamma21overM21_Bs(FULLNNLO, only1overmb).real() * SM.getOptionalParameter("DmBs_ex");
}

DGamma_d_NLO_tradBasis::DGamma_d_NLO_tradBasis(const StandardModel& SM_i)
: ThObservable(SM_i), AmpDB2(SM_i) {
    setParametersForObservable(make_vector<std::string>() << "DmBd_ex" );
}

DGamma_d_NLO_tradBasis::~DGamma_d_NLO_tradBasis() {
}

double DGamma_d_NLO_tradBasis::computeThValue() {
    return -getGamma21overM21_Bd_tradBasis(FULLNLO).real() * SM.getOptionalParameter("DmBd_ex");
}

DGamma_s_NLO_tradBasis::DGamma_s_NLO_tradBasis(const StandardModel& SM_i)
: ThObservable(SM_i), AmpDB2(SM_i) {
    setParametersForObservable(make_vector<std::string>() << "DmBs_ex" );
}

DGamma_s_NLO_tradBasis::~DGamma_s_NLO_tradBasis() {
}

double DGamma_s_NLO_tradBasis::computeThValue() {
    return -getGamma21overM21_Bs_tradBasis(FULLNLO).real() * SM.getOptionalParameter("DmBs_ex");
}

DGamma_s_LO_tradBasis::DGamma_s_LO_tradBasis(const StandardModel& SM_i)
: ThObservable(SM_i), AmpDB2(SM_i) {
    setParametersForObservable(make_vector<std::string>() << "DmBs_ex" );
}

DGamma_s_LO_tradBasis::~DGamma_s_LO_tradBasis() {
}

double DGamma_s_LO_tradBasis::computeThValue() {
    return -getGamma21overM21_Bs_tradBasis(LO).real() * SM.getOptionalParameter("DmBs_ex");
}
