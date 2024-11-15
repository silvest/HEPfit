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
: ThObservable(SM_i) {
    setParametersForObservable(make_vector<std::string>() << "DmBd_ex" );
            SM.getFlavour().getDB2(0);
}

DGamma_d_pole::~DGamma_d_pole() {
}

double DGamma_d_pole::computeThValue() {
    return -(SM.getFlavour().getDB2(0).getGamma21overM21(FULLNNLO, SM.getFlavour().getDB2(0).pole) / gslpp::complex(SM.getCBd(), -2. * SM.getPhiBd(), true)).real() * SM.getOptionalParameter("DmBd_ex");
}

DGamma_s_pole::DGamma_s_pole(const StandardModel& SM_i)
: ThObservable(SM_i) {
    setParametersForObservable(make_vector<std::string>() << "DmBs_ex" );
            SM.getFlavour().getDB2(1);
}

DGamma_s_pole::~DGamma_s_pole() {
}

double DGamma_s_pole::computeThValue() {
    return -(SM.getFlavour().getDB2(1).getGamma21overM21(FULLNNLO, SM.getFlavour().getDB2(1).pole) / gslpp::complex(SM.getCBs(), 2. * SM.getPhiBs(), true)).real() * SM.getOptionalParameter("DmBs_ex");
}

DGamma_s_pole_NLO::DGamma_s_pole_NLO(const StandardModel& SM_i)
: ThObservable(SM_i) {
    setParametersForObservable(make_vector<std::string>() << "DmBs_ex" );
            SM.getFlavour().getDB2(1);
}

DGamma_s_pole_NLO::~DGamma_s_pole_NLO() {
}

double DGamma_s_pole_NLO::computeThValue() {
    return -(SM.getFlavour().getDB2(1).getGamma21overM21(FULLNLO, SM.getFlavour().getDB2(1).pole) / gslpp::complex(SM.getCBs(), 2. * SM.getPhiBs(), true)).real() * SM.getOptionalParameter("DmBs_ex");
}

DGamma_s_pole_LO::DGamma_s_pole_LO(const StandardModel& SM_i)
: ThObservable(SM_i) {
    setParametersForObservable(make_vector<std::string>() << "DmBs_ex" );
            SM.getFlavour().getDB2(1);
}

DGamma_s_pole_LO::~DGamma_s_pole_LO() {
}

double DGamma_s_pole_LO::computeThValue() {
    return -(SM.getFlavour().getDB2(1).getGamma21overM21(LO, SM.getFlavour().getDB2(1).pole) / gslpp::complex(SM.getCBs(), 2. * SM.getPhiBs(), true)).real() * SM.getOptionalParameter("DmBs_ex");
}

DGamma_d_MSbar::DGamma_d_MSbar(const StandardModel& SM_i)
: ThObservable(SM_i) {
    setParametersForObservable(make_vector<std::string>() << "DmBd_ex" );
            SM.getFlavour().getDB2(0);
}

DGamma_d_MSbar::~DGamma_d_MSbar() {
}

double DGamma_d_MSbar::computeThValue() {
    return -(SM.getFlavour().getDB2(0).getGamma21overM21(FULLNNLO, SM.getFlavour().getDB2(0).MSbar) / gslpp::complex(SM.getCBd(), -2. * SM.getPhiBd(), true)).real() * SM.getOptionalParameter("DmBd_ex");
}

DGamma_s_MSbar::DGamma_s_MSbar(const StandardModel& SM_i)
: ThObservable(SM_i) {
    setParametersForObservable(make_vector<std::string>() << "DmBs_ex" );
            SM.getFlavour().getDB2(1);
}

DGamma_s_MSbar::~DGamma_s_MSbar() {
}

double DGamma_s_MSbar::computeThValue() {
    return -(SM.getFlavour().getDB2(1).getGamma21overM21(FULLNNLO, SM.getFlavour().getDB2(1).MSbar) / gslpp::complex(SM.getCBs(), 2. * SM.getPhiBs(), true)).real() * SM.getOptionalParameter("DmBs_ex");
}

DGamma_s_MSbar_NLO::DGamma_s_MSbar_NLO(const StandardModel& SM_i)
: ThObservable(SM_i) {
    setParametersForObservable(make_vector<std::string>() << "DmBs_ex" );
            SM.getFlavour().getDB2(1);
}

DGamma_s_MSbar_NLO::~DGamma_s_MSbar_NLO() {
}

double DGamma_s_MSbar_NLO::computeThValue() {
    return -(SM.getFlavour().getDB2(1).getGamma21overM21(FULLNLO, SM.getFlavour().getDB2(1).MSbar) / gslpp::complex(SM.getCBs(), 2. * SM.getPhiBs(), true)).real() * SM.getOptionalParameter("DmBs_ex");
}

DGamma_s_MSbar_LO::DGamma_s_MSbar_LO(const StandardModel& SM_i)
: ThObservable(SM_i) {
    setParametersForObservable(make_vector<std::string>() << "DmBs_ex" );
            SM.getFlavour().getDB2(1);
}

DGamma_s_MSbar_LO::~DGamma_s_MSbar_LO() {
}

double DGamma_s_MSbar_LO::computeThValue() {
    return -(SM.getFlavour().getDB2(1).getGamma21overM21(LO, SM.getFlavour().getDB2(1).MSbar) / gslpp::complex(SM.getCBs(), 2. * SM.getPhiBs(), true)).real() * SM.getOptionalParameter("DmBs_ex");
}

DGamma_d_PS::DGamma_d_PS(const StandardModel& SM_i)
: ThObservable(SM_i) {
    setParametersForObservable(make_vector<std::string>() << "DmBd_ex" );
            SM.getFlavour().getDB2(0);
}

DGamma_d_PS::~DGamma_d_PS() {
}

double DGamma_d_PS::computeThValue() {
    return -(SM.getFlavour().getDB2(0).getGamma21overM21(FULLNNLO, SM.getFlavour().getDB2(0).PS) / gslpp::complex(SM.getCBd(), -2. * SM.getPhiBd(), true)).real() * SM.getOptionalParameter("DmBd_ex");
}

DGamma_s_PS::DGamma_s_PS(const StandardModel& SM_i)
: ThObservable(SM_i) {
    setParametersForObservable(make_vector<std::string>() << "DmBs_ex" );
            SM.getFlavour().getDB2(1);
}

DGamma_s_PS::~DGamma_s_PS() {
}

double DGamma_s_PS::computeThValue() {
    return -(SM.getFlavour().getDB2(1).getGamma21overM21(FULLNNLO, SM.getFlavour().getDB2(1).PS) / gslpp::complex(SM.getCBs(), 2. * SM.getPhiBs(), true)).real() * SM.getOptionalParameter("DmBs_ex");
}

DGamma_s_PS_NLO::DGamma_s_PS_NLO(const StandardModel& SM_i)
: ThObservable(SM_i) {
    setParametersForObservable(make_vector<std::string>() << "DmBs_ex" );
            SM.getFlavour().getDB2(1);
}

DGamma_s_PS_NLO::~DGamma_s_PS_NLO() {
}

double DGamma_s_PS_NLO::computeThValue() {
    return -(SM.getFlavour().getDB2(1).getGamma21overM21(FULLNLO, SM.getFlavour().getDB2(1).PS) / gslpp::complex(SM.getCBs(), 2. * SM.getPhiBs(), true)).real() * SM.getOptionalParameter("DmBs_ex");
}

DGamma_s_PS_LO::DGamma_s_PS_LO(const StandardModel& SM_i)
: ThObservable(SM_i) {
    setParametersForObservable(make_vector<std::string>() << "DmBs_ex" );
            SM.getFlavour().getDB2(1);
}

DGamma_s_PS_LO::~DGamma_s_PS_LO() {
}

double DGamma_s_PS_LO::computeThValue() {
    return -(SM.getFlavour().getDB2(1).getGamma21overM21(LO, SM.getFlavour().getDB2(1).PS) / gslpp::complex(SM.getCBs(), 2. * SM.getPhiBs(), true)).real() * SM.getOptionalParameter("DmBs_ex");
}


DGamma_s_pole_fixmub::DGamma_s_pole_fixmub(const StandardModel& SM_i)
: ThObservable(SM_i){
    setParametersForObservable(make_vector<std::string>() << "DmBs_ex" );
            SM.getFlavour().getDB2(1);
}

DGamma_s_pole_fixmub::~DGamma_s_pole_fixmub() {
}

double DGamma_s_pole_fixmub::computeThValue() {
    return -(SM.getFlavour().getDB2(1, true).getGamma21overM21(FULLNNLO, SM.getFlavour().getDB2(1,true).pole) / gslpp::complex(SM.getCBs(), 2. * SM.getPhiBs(), true)).real() * SM.getOptionalParameter("DmBs_ex");
}

DGamma_s_MSbar_fixmub::DGamma_s_MSbar_fixmub(const StandardModel& SM_i)
: ThObservable(SM_i) {
    setParametersForObservable(make_vector<std::string>() << "DmBs_ex" );
            SM.getFlavour().getDB2(1);
}

DGamma_s_MSbar_fixmub::~DGamma_s_MSbar_fixmub() {
}

double DGamma_s_MSbar_fixmub::computeThValue() {
    return -(SM.getFlavour().getDB2(1,true).getGamma21overM21(FULLNNLO, SM.getFlavour().getDB2(1,true).MSbar) / gslpp::complex(SM.getCBs(), 2. * SM.getPhiBs(), true)).real() * SM.getOptionalParameter("DmBs_ex");
}

DGamma_s_PS_fixmub::DGamma_s_PS_fixmub(const StandardModel& SM_i)
: ThObservable(SM_i) {
    setParametersForObservable(make_vector<std::string>() << "DmBs_ex" );
            SM.getFlavour().getDB2(1);
}

DGamma_s_PS_fixmub::~DGamma_s_PS_fixmub() {
}

double DGamma_s_PS_fixmub::computeThValue() {
    return -(SM.getFlavour().getDB2(1,true).getGamma21overM21(FULLNNLO, SM.getFlavour().getDB2(1,true).PS) / gslpp::complex(SM.getCBs(), 2. * SM.getPhiBs(), true)).real() * SM.getOptionalParameter("DmBs_ex");
}

DGamma_d_only1overmb::DGamma_d_only1overmb(const StandardModel& SM_i)
: ThObservable(SM_i) {
    setParametersForObservable(make_vector<std::string>() << "DmBd_ex" );
            SM.getFlavour().getDB2(0);
}

DGamma_d_only1overmb::~DGamma_d_only1overmb() {
}

double DGamma_d_only1overmb::computeThValue() {
    return -(SM.getFlavour().getDB2(0).getGamma21overM21(FULLNNLO, SM.getFlavour().getDB2(0).only1overmb) / gslpp::complex(SM.getCBd(), -2. * SM.getPhiBd(), true)).real() * SM.getOptionalParameter("DmBd_ex");
}

DGamma_s_only1overmb::DGamma_s_only1overmb(const StandardModel& SM_i)
: ThObservable(SM_i) {
    setParametersForObservable(make_vector<std::string>() << "DmBs_ex" );
            SM.getFlavour().getDB2(1);
}

DGamma_s_only1overmb::~DGamma_s_only1overmb() {
}

double DGamma_s_only1overmb::computeThValue() {
    return -(SM.getFlavour().getDB2(1).getGamma21overM21(FULLNNLO, SM.getFlavour().getDB2(1).only1overmb) / gslpp::complex(SM.getCBs(), 2. * SM.getPhiBs(), true)).real() * SM.getOptionalParameter("DmBs_ex");
}

DGamma_d_NLO_tradBasis::DGamma_d_NLO_tradBasis(const StandardModel& SM_i)
: ThObservable(SM_i) {
    setParametersForObservable(make_vector<std::string>() << "DmBd_ex" );
            SM.getFlavour().getDB2(0);
}

DGamma_d_NLO_tradBasis::~DGamma_d_NLO_tradBasis() {
}

double DGamma_d_NLO_tradBasis::computeThValue() {
    return -(SM.getFlavour().getDB2(0).getGamma21overM21_tradBasis(FULLNLO) / gslpp::complex(SM.getCBd(), -2. * SM.getPhiBd(), true)).real() * SM.getOptionalParameter("DmBd_ex");
}

DGamma_s_NLO_tradBasis::DGamma_s_NLO_tradBasis(const StandardModel& SM_i)
: ThObservable(SM_i) {
    setParametersForObservable(make_vector<std::string>() << "DmBs_ex" );
            SM.getFlavour().getDB2(1);
}

DGamma_s_NLO_tradBasis::~DGamma_s_NLO_tradBasis() {
}

double DGamma_s_NLO_tradBasis::computeThValue() {
    return -(SM.getFlavour().getDB2(1).getGamma21overM21_tradBasis(FULLNLO) / gslpp::complex(SM.getCBs(), 2. * SM.getPhiBs(), true)).real() * SM.getOptionalParameter("DmBs_ex");
}

DGamma_d_LO_tradBasis::DGamma_d_LO_tradBasis(const StandardModel& SM_i)
: ThObservable(SM_i) {
    setParametersForObservable(make_vector<std::string>() << "DmBd_ex" );
            SM.getFlavour().getDB2(0);
}

DGamma_d_LO_tradBasis::~DGamma_d_LO_tradBasis() {
}

double DGamma_d_LO_tradBasis::computeThValue() {
    return -(SM.getFlavour().getDB2(0).getGamma21overM21_tradBasis(LO) / gslpp::complex(SM.getCBd(), -2. * SM.getPhiBd(), true)).real() * SM.getOptionalParameter("DmBd_ex");
}

DGamma_s_LO_tradBasis::DGamma_s_LO_tradBasis(const StandardModel& SM_i)
: ThObservable(SM_i) {
    setParametersForObservable(make_vector<std::string>() << "DmBs_ex" );
            SM.getFlavour().getDB2(1);
}

DGamma_s_LO_tradBasis::~DGamma_s_LO_tradBasis() {
}

double DGamma_s_LO_tradBasis::computeThValue() {
    return -(SM.getFlavour().getDB2(1).getGamma21overM21_tradBasis(LO) / gslpp::complex(SM.getCBs(), 2. * SM.getPhiBs(), true)).real() * SM.getOptionalParameter("DmBs_ex");
}

DGamma_s_MSbar_NLO_RI::DGamma_s_MSbar_NLO_RI(const StandardModel& SM_i)
: ThObservable(SM_i) {
    setParametersForObservable(make_vector<std::string>() << "DmBs_ex" );
            SM.getFlavour().getDB2(1);
}

DGamma_s_MSbar_NLO_RI::~DGamma_s_MSbar_NLO_RI() {
}

double DGamma_s_MSbar_NLO_RI::computeThValue() {
    return -(SM.getFlavour().getDB2(1,false, true).getGamma21overM21(FULLNLO, SM.getFlavour().getDB2(1,false,true).MSbar) / gslpp::complex(SM.getCBs(), 2. * SM.getPhiBs(), true)).real() * SM.getOptionalParameter("DmBs_ex");
}

DGamma_s_MSbar_NLO_RI_tradBasis::DGamma_s_MSbar_NLO_RI_tradBasis(const StandardModel& SM_i)
: ThObservable(SM_i) {
    setParametersForObservable(make_vector<std::string>() << "DmBs_ex" );
            SM.getFlavour().getDB2(1);
}

DGamma_s_MSbar_NLO_RI_tradBasis::~DGamma_s_MSbar_NLO_RI_tradBasis() {
}

double DGamma_s_MSbar_NLO_RI_tradBasis::computeThValue() {
    return -(SM.getFlavour().getDB2(1,false,true).getGamma21overM21_tradBasis(FULLNLO) / gslpp::complex(SM.getCBs(), 2. * SM.getPhiBs(), true)).real() * SM.getOptionalParameter("DmBs_ex");
}

DGamma_s_PS_NLO_RI::DGamma_s_PS_NLO_RI(const StandardModel& SM_i)
: ThObservable(SM_i) {
    setParametersForObservable(make_vector<std::string>() << "DmBs_ex" );
            SM.getFlavour().getDB2(1);
}

DGamma_s_PS_NLO_RI::~DGamma_s_PS_NLO_RI() {
}

double DGamma_s_PS_NLO_RI::computeThValue() {
    return -(SM.getFlavour().getDB2(1,false,true).getGamma21overM21(FULLNLO, SM.getFlavour().getDB2(1,false,true).PS) / gslpp::complex(SM.getCBs(), 2. * SM.getPhiBs(), true)).real() * SM.getOptionalParameter("DmBs_ex");
}

DGamma_s_MSbar_partialNNLO::DGamma_s_MSbar_partialNNLO(const StandardModel& SM_i)
: ThObservable(SM_i) {
    setParametersForObservable(make_vector<std::string>() << "DmBs_ex" );
            SM.getFlavour().getDB2(1);
}

DGamma_s_MSbar_partialNNLO::~DGamma_s_MSbar_partialNNLO() {
}

double DGamma_s_MSbar_partialNNLO::computeThValue() {
    return -(SM.getFlavour().getDB2(1).getGamma21overM21(FULLNNLO, SM.getFlavour().getDB2(1).MSbar_partialNNLO) / gslpp::complex(SM.getCBs(), 2. * SM.getPhiBs(), true)).real() * SM.getOptionalParameter("DmBs_ex");
}

DGamma_s_PS_partialNNLO::DGamma_s_PS_partialNNLO(const StandardModel& SM_i)
: ThObservable(SM_i) {
    setParametersForObservable(make_vector<std::string>() << "DmBs_ex" );
            SM.getFlavour().getDB2(1);
}

DGamma_s_PS_partialNNLO::~DGamma_s_PS_partialNNLO() {
}

double DGamma_s_PS_partialNNLO::computeThValue() {
    return -(SM.getFlavour().getDB2(1).getGamma21overM21(FULLNNLO, SM.getFlavour().getDB2(1).PS_partialNNLO) / gslpp::complex(SM.getCBs(), 2. * SM.getPhiBs(), true)).real() * SM.getOptionalParameter("DmBs_ex");
}

DGamma_s_MSbar_partialN3LO::DGamma_s_MSbar_partialN3LO(const StandardModel& SM_i)
: ThObservable(SM_i) {
    setParametersForObservable(make_vector<std::string>() << "DmBs_ex" );
            SM.getFlavour().getDB2(1);
}

DGamma_s_MSbar_partialN3LO::~DGamma_s_MSbar_partialN3LO() {
}

double DGamma_s_MSbar_partialN3LO::computeThValue() {
    return -(SM.getFlavour().getDB2(1).getGamma21overM21(FULLNNNLO, SM.getFlavour().getDB2(1).MSbar_partialN3LO) / gslpp::complex(SM.getCBs(), 2. * SM.getPhiBs(), true)).real() * SM.getOptionalParameter("DmBs_ex");
}

DGamma_s_PS_partialN3LO::DGamma_s_PS_partialN3LO(const StandardModel& SM_i)
: ThObservable(SM_i) {
    setParametersForObservable(make_vector<std::string>() << "DmBs_ex" );
            SM.getFlavour().getDB2(1);
}

DGamma_s_PS_partialN3LO::~DGamma_s_PS_partialN3LO() {
}

double DGamma_s_PS_partialN3LO::computeThValue() {
    return -(SM.getFlavour().getDB2(1).getGamma21overM21(FULLNNNLO, SM.getFlavour().getDB2(1).PS_partialN3LO) / gslpp::complex(SM.getCBs(), 2. * SM.getPhiBs(), true)).real() * SM.getOptionalParameter("DmBs_ex");
}