/* 
 * File:   StandardModel.cpp
 * Author: silvest
 * 
 * Created on November 30, 2010, 1:27 PM
 */

#include "StandardModel.h"
#include <math.h>

StandardModel::StandardModel(const gslpp::matrix<gslpp::complex>& VCKM_i,
        double mu_i, double md_i, double mc_i, double ms_i, double mt_i,
        double mb_i, const gslpp::matrix<gslpp::complex>& UPMNS_i, double me_i,
        double mmu_i, double mtau_i, double mnu1_i, double mnu2_i,
        double mnu3_i) : UPMNS(UPMNS_i), VCKM(VCKM_i),
        Yd(3,3,0.), Yu(3,3,0.), Ye(3,3,0.), Yn(3,3,0.) {
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
    Yd.assign(0,0,md/v()*sqrt(2.));
    Yd.assign(1,1,ms/v()*sqrt(2.));
    Yd.assign(2,2,mb/v()*sqrt(2.));
    Yu.assign(0,0,mu/v()*sqrt(2.));
    Yu.assign(1,1,mc/v()*sqrt(2.));
    Yu.assign(2,2,mt/v()*sqrt(2.));
    Yu = Yu*VCKM;
    GF = 1.;
}

StandardModel::StandardModel(const StandardModel& orig) :
        UPMNS(orig.getUPMNS()), VCKM(orig.getVCKM()), Yd(3,3,0.), Yu(3,3,0.),
        Ye(3,3,0.), Yn(3,3,0.) {
    StandardModel(orig.getVCKM(), orig.getMu(), orig.getMd(), orig.getMc(),
            orig.getMs(), orig.getMt(), orig.getMb(), orig.getUPMNS(),
            orig.getMe(), orig.getMmu(), orig.getMtau(), orig.getMnu1(),
            orig.getMnu2(), orig.getMnu3());
}

StandardModel::~StandardModel() {
}

double StandardModel::v() {
    return 1./sqrt(sqrt(2.)*GF);
}