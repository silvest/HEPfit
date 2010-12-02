/* 
 * File:   MFV.cpp
 * Author: silvest
 * 
 * Created on September 24, 2010, 10:53 AM
 */

#include "MFV.h"

MFV::MFV(double mQtilde_i, double mUtilde_i, double mDtilde_i,
        double Au_i, double Ad_i, double mLtilde_i, double mEtilde_i,
        double mNtilde_i, double Ae_i, double An_i, double m1_i,
        double m2_i, double m3_i, double muH_i, double tanb_i, double mHp_i,
        const gslpp::matrix<gslpp::complex>& VCKM_i,
        double mu_i, double md_i, double mc_i, double ms_i, double mt_i,
        double mb_i, const gslpp::matrix<gslpp::complex>& UPMNS_i,
        double me_i, double mmu_i, double mtau_i, double mnu1_i,
        double mnu2_i, double mnu3_i) : StandardModel(VCKM_i, mu_i, md_i, mc_i,
        ms_i, mt_i, mb_i, UPMNS_i, me_i, mmu_i, mtau_i,mnu1_i, mnu2_i, mnu3_i)
{
    mQtilde = mQtilde_i;
    mUtilde = mUtilde_i;
    mDtilde = mDtilde_i,
    Au = Au_i;
    Ad = Ad_i;
    mLtilde = mLtilde_i;
    mEtilde = mEtilde_i;
    mNtilde = mNtilde_i;
    Ae = Ae_i;
    An = An_i;
    m1 = m1_i;
    m2 = m2_i;
    m3 = m3_i;
    muH = muH_i;
    tanb = tanb_i;
    mHp = mHp_i;
}

MFV::MFV(double mQtilde_i, double mUtilde_i, double mDtilde_i,
        double Au_i, double Ad_i, double mLtilde_i, double mEtilde_i,
        double mNtilde_i, double Ae_i, double An_i, double m1_i,
        double m2_i, double m3_i, double muH_i, double tanb_i, double mHp_i,
        const StandardModel& SM_i) : StandardModel(SM_i)
{
    mQtilde = mQtilde_i;
    mUtilde = mUtilde_i;
    mDtilde = mDtilde_i,
    Au = Au_i;
    Ad = Ad_i;
    mLtilde = mLtilde_i;
    mEtilde = mEtilde_i;
    mNtilde = mNtilde_i;
    Ae = Ae_i;
    An = An_i;
    m1 = m1_i;
    m2 = m2_i;
    m3 = m3_i;
    muH = muH_i;
    tanb = tanb_i;
    mHp = mHp_i;
}

MFV::MFV(const MFV& orig) : StandardModel(orig.GetVCKM(), orig.GetMu(),
        orig.GetMd(), orig.GetMc(), orig.GetMs(), orig.GetMt(), orig.GetMb(),
        orig.GetUPMNS(), orig.GetMe(), orig.GetMmu(), orig.GetMtau(),
        orig.GetMnu1(), orig.GetMnu2(), orig.GetMnu3())
{
    mQtilde = orig.GetMQtilde();
    mUtilde = orig.GetMUtilde();
    mDtilde = orig.GetMDtilde(),
    Au = orig.GetAu();
    Ad = orig.GetAd();
    mLtilde = orig.GetMLtilde();
    mEtilde = orig.GetMEtilde();
    mNtilde = orig.GetMNtilde();
    Ae = orig.GetAe();
    An = orig.GetAn();
    m1 = orig.GetM1();
    m2 = orig.GetM2();
    m3 = orig.GetM3();
    muH = orig.GetMuH();
    tanb = orig.GetTanb();
    mHp = orig.GetMHp();
}

MFV::~MFV() {
}
