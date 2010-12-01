/* 
 * File:   StandardModel.cpp
 * Author: silvest
 * 
 * Created on November 30, 2010, 1:27 PM
 */

#include "StandardModel.h"

StandardModel::StandardModel(const gslpp::matrix<gslpp::complex>& VCKM_i,
        double mu_i, double md_i, double mc_i, double ms_i, double mt_i,
        double mb_i, const gslpp::matrix<gslpp::complex>& UPMNS_i, double me_i,
        double mmu_i, double mtau_i, double mnu1_i, double mnu2_i,
        double mnu3_i) : UPMNS(UPMNS_i),VCKM(VCKM_i) {
    mu = mu_i;
    md = md_i;
    mc = mc_i;
    ms = ms_i;
    mt = mt_i;
    mb = mb_i;
    me = me_i;
    mmu = mmu_i;
    mtau = mtau_i;
    mnu1 = mnu1_i;
    mnu2 = mnu2_i;
    mnu3 = mnu3_i;
}

StandardModel::StandardModel(const StandardModel& orig) :
UPMNS(orig.GetUPMNS()),VCKM(orig.GetVCKM()) {
    mu = orig.GetMu();
    md = orig.GetMd();
    mc = orig.GetMc();
    ms = orig.GetMs();
    mt = orig.GetMt();
    mb = orig.GetMb();
    me = orig.GetMe();
    mmu = orig.GetMmu();
    mtau = orig.GetMtau();
    mnu1 = orig.GetMnu1();
    mnu2 = orig.GetMnu2();
    mnu3 = orig.GetMnu3();
}

StandardModel::~StandardModel() {
}
