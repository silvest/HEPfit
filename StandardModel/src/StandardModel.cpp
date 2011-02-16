/* 
 * File:   StandardModel.cpp
 * Author: silvest
 * 
 * Created on November 30, 2010, 1:27 PM
 */

#include "StandardModel.h"
#include <boost/assign/list_of.hpp> // for 'map_list_of()'
#include <math.h>

/*
const std::map<std::string,std::string[] > Deps =
    boost::assign::map_list_of
        ("mZ",{"mW"})
        ("mt",{"mW"})
        ("dAle5Mz",{"mW"})
        ("alsMz",{"mW"})
        ("mHl",{"mW"})
        ("GF",{"v"});
*/

StandardModel::StandardModel(const gslpp::matrix<gslpp::complex>& VCKM_i,
        double mu_i, double md_i, double mc_i, double ms_i, double mt_i,
        double mb_i, const gslpp::matrix<gslpp::complex>& UPMNS_i, double me_i,
        double mmu_i, double mtau_i, double mnu1_i, double mnu2_i,
        double mnu3_i, double GF_i, double alsMz_i, double ale_i, double mZ_i,
        double dAle5Mz_i, double mHl_i) : UPMNS(UPMNS_i), VCKM(VCKM_i),
        Yd(3,3,0.), Yu(3,3,0.), Ye(3,3,0.), Yn(3,3,0.) {
    GF = GF_i;
    alsMz = alsMz_i;
    ale = ale_i;
    mZ = mZ_i;
    dAle5Mz = dAle5Mz_i;
    mHl = mHl_i;
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
    Ye.assign(0,0,me/v()*sqrt(2.));
    Ye.assign(1,1,mmu/v()*sqrt(2.));
    Ye.assign(2,2,mtau/v()*sqrt(2.));
    Yn.assign(0,0,mnu1/v()*sqrt(2.));
    Yn.assign(1,1,mnu2/v()*sqrt(2.));
    Yn.assign(2,2,mnu3/v()*sqrt(2.));
    Yn = Yn*UPMNS.hconjugate();
}

StandardModel::StandardModel(Parameters& Par): UPMNS(3,3,0.), VCKM(3,3,0.),
        Yd(3,3,0.), Yu(3,3,0.), Ye(3,3,0.), Yn(3,3,0.) {
    Par.Get("mu",mu);
    Par.Get("md",md);
    Par.Get("mc",mc);
    Par.Get("ms",ms);
    Par.Get("mt",mt);
    Par.Get("mb",mb);
    Par.Get("me",me);
    Par.Get("mmu",mmu);
    Par.Get("mtau",mtau);
    Par.Get("mnu1",mnu1);
    Par.Get("mnu2",mnu2);
    Par.Get("mnu3",mnu3);
    Yd.assign(0,0,md/v()*sqrt(2.));
    Yd.assign(1,1,ms/v()*sqrt(2.));
    Yd.assign(2,2,mb/v()*sqrt(2.));
    Yu.assign(0,0,mu/v()*sqrt(2.));
    Yu.assign(1,1,mc/v()*sqrt(2.));
    Yu.assign(2,2,mt/v()*sqrt(2.));
    Par.Get("VCKM",VCKM);
    Par.Get("UPMNS",UPMNS);
    Yu = Yu*VCKM;
    Par.Get("GF", GF);
    Par.Get("alsMz", alsMz);
    Par.Get("ale", ale);
    Par.Get("mZ", mZ);
    Par.Get("dAle5Mz", dAle5Mz);
    Par.Get("mHl", mHl);
}


StandardModel::StandardModel(const StandardModel& orig) :
        UPMNS(orig.getUPMNS()), VCKM(orig.getVCKM()), Yd(3,3,0.), Yu(3,3,0.),
        Ye(3,3,0.), Yn(3,3,0.) {
    StandardModel(orig.getVCKM(), orig.getMu(), orig.getMd(), orig.getMc(),
            orig.getMs(), orig.getMt(), orig.getMb(), orig.getUPMNS(),
            orig.getMe(), orig.getMmu(), orig.getMtau(), orig.getMnu1(),
            orig.getMnu2(), orig.getMnu3(), orig.getGF(), orig.getAlsMz(),
            orig.getAle(), orig.getMZ(), orig.getDAle5Mz(), orig.getMHl() );
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

    // mt, mZ, dAle5Mz and alsMz have to be varied within their combined 
    // 2 sigma region around their central values (year 2003) adopted below.
    double dH = log(mHl/100.0);
    double dh = pow((mHl/100.0), 2.0);
    double dt = pow((mt/174.3), 2.0) - 1.0;
    double dZ = mZ/91.1875 - 1.0;
    double dalphae = dAle5Mz/0.05907 - 1.0;
    double dalphas = alsMz/0.119 - 1.0;

    double mw = Mw0 - c1*dH - c2*dH*dH + c3*pow(dH, 4.0)
                + c4*(dh - 1.0) - c5*dalphae + c6*dt - c7*dt*dt
                - c8*dH*dt + c9*dh*dt - c10*dalphas + c11*dZ;
    return mw;
}

double StandardModel::sin2thw() const {
    // Effective leptonic weak mixing angle
    // hep-ph/0407317, hep-ph/0608099
    // applicable for 10 GeV <= mHl <= 1 TeV

    const double s0 = 0.2312527;
    const double d1 = 4.729*0.0001;
    const double d2 = 2.07*0.00001;
    const double d3 = 3.85*0.000001;
    const double d4 = -1.85*0.000001;
    const double d5 = 0.0207;
    const double d6 = -0.002851;
    const double d7 = 1.82*0.0001;
    const double d8 = -9.74*0.000001;
    const double d9 = 3.98*0.0001;
    const double d10 = -0.655;

    // mt, mZ, dAle5Mz and alsMz have to be varied within their combined
    // 2 sigma region around their central values (year 2003) adopted below.
    double L_H = log(mHl/100.0);
    double Delta_H = mHl/100.0;
    double Delta_alphae = dAle5Mz/0.05907 - 1.0;
    double Delta_t = pow((mt/178.0), 2.0) - 1.0;
    double Delta_alphas = alsMz/0.117 - 1.0;
    double Delta_Z = mZ/91.1876 - 1.0;

    double sin2t = s0 + d1*L_H + d2*L_H*L_H + d3*L_H*L_H*L_H*L_H
                 + d4*(Delta_H*Delta_H - 1.0) + d5*Delta_alphae + d6*Delta_t
                 + d7*Delta_t*Delta_t + d8*Delta_t*(Delta_H - 1.0)
                 + d9*Delta_alphas + d10*Delta_Z;
    return sin2t;
}

    double StandardModel::sin2thwb() const{
    //Effective mixing angle for b quarks
    //http://arXiv.org/abs/0811.1364v2
    // applicable for 10 GeV <= mHl <= 1 TeV
        double const s0=0.2327580;
        double const d1=4.749*0.0001;  
        double const d2=2.03*0.00001;  
        double const d3=3.94*0.000001;
        double const d4=-1.84*0.000001;
        double const d5=2.08*0.01;
        double const d6=-9.93*0.0001;
        double const d7=7.08*0.00001;
        double const d8=-7.61*0.000001;
        double const d9=4.03*0.0001;
        double const d10=0.661;
        
        double deltaa=dAle5Mz/0.05907-1.0;
        double deltas=alsMz/0.117-1.0;
        double  Lh=log(mHl/100.0);
        double deltH=mHl/100.0;
        double deltt=(mt/178)*(mt/178)-1.0;
        double deltz=mZ/91.1876-1.0;
        double deltw=mW()/80.404-1.0;
        double sin2b=s0 + d1*Lh + d2*Lh*Lh + d3*pow(Lh,4.0)+d4*(deltH*deltH-1.0)+d5*deltaa
                     + d6*deltt + d7*deltt*deltt + d8*deltt*(deltH-1.0) + d9*deltas + d10*deltz  ;
        return sin2b;
     }

double StandardModel::GammaW() const {
    return 0.0;
}

double StandardModel::GammaZ() const {
    return 0.0;
}

double StandardModel::sigma_had() const {
    return 0.0;
}

double StandardModel::R_l() const {
    return 0.0;
}

double StandardModel::R_c() const {
    return 0.0;
}

double StandardModel::R_b() const {
    return 0.0;
}

double StandardModel::AFB_l() const {
    return 0.0;
}

double StandardModel::AFB_c() const {
    return 0.0;
}

double StandardModel::AFB_b() const {
    return 0.0;
}

double StandardModel::A_l() const {
    return 0.0;
}

double StandardModel::A_c() const {
    return 0.0;
}

double StandardModel::A_b() const {
    return 0.0;
}

double StandardModel::S() const {
    return 0.0;
}

double StandardModel::T() const {
    return 0.0;
}

double StandardModel::U() const {
    return 0.0;
}




