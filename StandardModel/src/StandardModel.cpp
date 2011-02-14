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

double StandardModel::mW() const {
    // Eqs. (6), (7) and (9) in hep-ph/0311148
    // applicable for 100 GeV <= mHl <= 1 TeV

    const double Mw0 = 80.3800;
    const double c1 = 0.05253;
    const double c2 = 0.010345;
    const double c3 = 0.001021;
    const double c4 = -0.000070;
    const double c5 = 1.077;
    const double c6 = 0.5270;
    const double c7 = 0.0698;
    const double c8 = 0.004055;
    const double c9 = 0.000110;
    const double c10 = 0.0716;
    const double c11 = 115.0;

    // mt, mZ, deltaalphae and alphasmZ have to be varied within their
    // combined 2 sigma region around their central values adopted below.
    double dH = log(mHl/100.0);
    double dh = pow((mHl/100.0), 2.0);
    double dt = pow((mt/174.3), 2.0) - 1.0;
    double dZ = mZ/91.1875 - 1.0;
    double dalphae = deltaalphae/0.05907 - 1.0;
    double dalphas = alphasmZ/0.119 - 1.0;

    double mw = Mw0 - c1*dH - c2*dH*dH + c3*pow(dH, 4.0)
                + c4*(dh - 1.0) - c5*dalphae + c6*dt - c7*dt*dt
                - c8*dH*dt + c9*dh*dt - c10*dalphas + c11*dZ;
    return mw;
}


