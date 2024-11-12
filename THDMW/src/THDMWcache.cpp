/* 
 * Copyright (C) 2017 HEPfit Collaboration
 * All rights reserved.
 *
 * For the licensing terms see doc/COPYING.
 */

#include "THDMWcache.h"
#include <fstream>
#include "gslpp.h"
#include <sstream>
#include <string>

THDMWcache::THDMWcache(const StandardModel& SM_i)
        :unitarityeigenvalues(11, 0.),
        NLOunitarityeigenvalues(11, 0.),
        myTHDMW(static_cast<const THDMW*> (&SM_i)),
        PV(false),
        ATLAS8_gg_phi_tt(53, 2, 0.),
        ATLAS8_gg_phi_tt_e(53, 2, 0.),
        CMS8_pp_H_hh_bbbb(167, 2, 0.),
        CMS8_pp_H_hh_bbbb_e(167, 2, 0.),
        CMS8_bb_phi_bb(81, 2, 0.),
        CMS8_bb_phi_bb_e(81, 2, 0.),
        Dummy(167, 2, 0.),
        ATLAS13_bb_phi_tt(61,2,0.),
        ATLAS13_tt_phi_tt(61,2,0.),
        ATLAS13_pp_H_hh_bbbb(271,2,0.),
        ATLAS13_bb_phi_tt_e(61,2,0.),
        ATLAS13_tt_phi_tt_e(61,2,0.),
        ATLAS13_pp_H_hh_bbbb_e(271,2,0.),
        CMS13_pp_phi_bb(66,2,0.),
        CMS8_pp_phi_bb(88,2,0.),
        CMS13_pp_H_hh_bbbb(95,2,0.),
        CMS13_ggF_H_hh_bbbb(226,2,0.),
        CMS13_pp_phi_bb_e(66,2,0.),
        CMS13_pp_H_hh_bbbb_e(95,2,0.),
        CMS13_ggF_H_hh_bbbb_e(226,2,0.),
        CMS13_pp_R_gg(241,2,0.),
        
        ATLAS8_pp_Hpm_tb(41,2,0.),
        ATLAS8_pp_Hpm_tb_e(41,2,0.),
        CMS8_pp_Hp_tb(43,2,0.),
        CMS8_pp_Hp_tb_e(43,2,0.),
        CMS13_bb_H_bb(101,2,0.),
        ATLAS13_pp_Hp_tb(181,2,0.),
//        ATLAS13_pp_Hp_tb1(71,2,0.),
//        ATLAS13_pp_Hp_tb2(181,2,0.),
//        ATLAS13_pp_Hp_tb1_e(71,2,0.),
//        ATLAS13_pp_Hp_tb2_e(181,2,0.),
        ATLAS13_pp_Gkk_tt(131,2,0.),
        ATLAS13_pp_SS_jjjj(126,2,0.),
        MadGraph_pp_Sr_tt(22800,5,0.),
        MadGraph_pp_Srtt_tttt(22800,5,0.),
        MadGraph_pp_Sr_jj(2940,5,0.),
        MadGraph_pp_SrSr_jjjj(4200,5,0.),
        MadGraph_pp_Stb_tbtb(4332,4,0.),
        MadGraph_pp_Sitt_tttt(9360,4,0.),
        MadGraph_pp_Srbb_bbbb(15960,5,0.),
        MadGraph_pp_Srbb_bbbb_8TeV(15960,5,0.),
        MadGraph_pp_Sibb_bbbb(8892,4,0.),
        MadGraph_pp_Sibb_bbbb_8TeV(8892,4,0.),
        MadGraph_pp_Sr_bb(15960,5,0.),
        MadGraph_pp_Sr_bb_8TeV(15960,5,0.),
        MadGraph_pp_Si_bb(8892,4,0.),
        MadGraph_pp_Si_bb_8TeV(8892,4,0.),
        arraybsgamma(1111, 3, 0.),
        betaeigenvalues(11, 0.)
        //myTHDMW(static_cast<const THDMW*> (&SM_i))
        


{
    myRunnerTHDMW=new RunnerTHDMW(SM_i);
    read();  
}

THDMWcache::~THDMWcache()
{
  delete myRunnerTHDMW;
}
//
///////////////////////////////////////////////////////////////////////////////////////////////////
//
int THDMWcache::CacheCheck(const gslpp::complex cache[][CacheSize], 
                          const int NumPar, const double params[]) const {
    bool bCache;
    for(int i=0; i<CacheSize; i++) {
        bCache = true;
        for(int j=0; j<NumPar; j++)
            bCache &= (params[j] == cache[j][i].real());
        if (bCache) return i;
    }
    return -1;
}

int THDMWcache::CacheCheckReal(const double cache[][CacheSize], 
                          const int NumPar, const double params[]) const {
    bool bCache;
    for(int i=0; i<CacheSize; i++) {
        bCache = true;
        for(int j=0; j<NumPar; j++)
            bCache &= (params[j] == cache[j][i]);
        if (bCache) return i;
    }
    return -1;
}

void THDMWcache::CacheShift(gslpp::complex cache[][CacheSize], const int NumPar, 
                           const double params[], const gslpp::complex newResult) const {
    // shift old parameters and result
    for(int i=CacheSize-1; i>0; i--)
        for(int j=0; j<NumPar+1; j++)
            cache[j][i] = cache[j][i-1];

    // store new parameters and result
    for(int j=0; j<NumPar; j++) {
        cache[j][0] = gslpp::complex(params[j], 0.0, false);
        cache[NumPar][0] = newResult;
    }
}

void THDMWcache::CacheShiftReal(double cache[][CacheSize], const int NumPar,
                           const double params[], const double newResult) const {
    // shift old parameters and result
    for(int i=CacheSize-1; i>0; i--)
        for(int j=0; j<NumPar+1; j++)
            cache[j][i] = cache[j][i-1];

    // store new parameters and result
    for(int j=0; j<NumPar; j++) {
        cache[j][0] = params[j];
        cache[NumPar][0] = newResult;
    }
}

////////////////////////////////////////////////////
//--------- Passarino Veltman Functions for THDMW ---------//
////////////////////////////////////////////////////

gslpp::complex THDMWcache::A0_MZ2_mSp2(const double MZ2, const double mSp2) const {
    int NumPar = 2;
    double params[] = {MZ2, mSp2};

    int i = CacheCheck(A0_MZ2_mSp2_cache, NumPar, params);
    if (i>=0) {
        return ( A0_MZ2_mSp2_cache[NumPar][i] );
    } else {
        gslpp::complex newResult = PV.A0(MZ2, mSp2); 
        CacheShift(A0_MZ2_mSp2_cache, NumPar, params, newResult);
        return newResult;
    } 
}

gslpp::complex THDMWcache::A0_MZ2_mSr2(const double MZ2, const double mSr2) const {
    int NumPar = 2;
    double params[] = {MZ2, mSr2};

    int i = CacheCheck(A0_MZ2_mSr2_cache, NumPar, params);
    if (i>=0) {
        return ( A0_MZ2_mSr2_cache[NumPar][i] );
    } else {
        gslpp::complex newResult = PV.A0(MZ2, mSr2);
        CacheShift(A0_MZ2_mSr2_cache, NumPar, params, newResult);
        return newResult;
    } 
}

gslpp::complex THDMWcache::A0_MZ2_mSi2(const double MZ2, const double mSi2) const {
    int NumPar = 2;
    double params[] = {MZ2, mSi2};

    int i = CacheCheck(A0_MZ2_mSi2_cache, NumPar, params);
    if (i>=0) {
        return ( A0_MZ2_mSi2_cache[NumPar][i] );
    } else {
        gslpp::complex newResult = PV.A0(MZ2, mSi2);
        CacheShift(A0_MZ2_mSi2_cache, NumPar, params, newResult);
        return newResult;
    } 
}

gslpp::complex THDMWcache::B0_MZ2_0_mSp2_mSp2(const double MZ2, const double mSp2) const {
    int NumPar = 2;
    double params[] = {MZ2, mSp2};

    int i = CacheCheck(B0_MZ2_0_mSp2_mSp2_cache, NumPar, params);
    if (i>=0) {
        return ( B0_MZ2_0_mSp2_mSp2_cache[NumPar][i] );
    } else {
        gslpp::complex newResult = PV.B0(MZ2,0. ,mSp2 , mSp2);
        CacheShift(B0_MZ2_0_mSp2_mSp2_cache, NumPar, params, newResult);
        return newResult;
    } 
}
 /*gslpp::complex THDMWcache::B00_MZ2_0_mSr2_mSp2(const double MZ2, const double mSr2, const double mSp2) const {
    int NumPar = 3;
    double params[] = {MZ2, mSr2, mSp2};

    int i = CacheCheck(B00_MZ2_0_mSr2_mSp2_cache, NumPar, params);
    if (i>=0) {
        return ( B00_MZ2_0_mSr2_mSp2_cache[NumPar][i] );
    } else {
        gslpp::complex newResult = PV.B00(MZ2, 0. , mSr2, mSp2);
        CacheShift(B00_MZ2_0_mSr2_mSp2_cache, NumPar, params, newResult);
        return newResult;
    } 
}

gslpp::complex THDMWcache::B00_MZ2_0_mSi2_mSp2(const double MZ2, const double mSi2, const double mSp2) const {
    int NumPar = 3;
    double params[] = {MZ2, mSi2, mSp2};

    int i = CacheCheck(B00_MZ2_0_mSi2_mSp2_cache, NumPar, params);
    if (i>=0) {
        return ( B00_MZ2_0_mSi2_mSp2_cache[NumPar][i] );
    } else {
        gslpp::complex newResult = PV.B00(MZ2, 0. , mSi2, mSp2);
        CacheShift(B00_MZ2_0_mSi2_mSp2_cache, NumPar, params, newResult);
        return newResult;
    } 
}

gslpp::complex THDMWcache::B00_MZ2_0_mSp2_mSp2(const double MZ2, const double mSp2) const {
    int NumPar = 2;
    double params[] = {MZ2, mSp2};

    int i = CacheCheck(B00_MZ2_0_mSp2_mSp2_cache, NumPar, params);
    if (i>=0) {
        return ( B00_MZ2_0_mSp2_mSp2_cache[NumPar][i] );
    } else {
        gslpp::complex newResult = PV.B00(MZ2, 0. , mSp2, mSp2);
        CacheShift(B00_MZ2_0_mSp2_mSp2_cache, NumPar, params, newResult);
        return newResult;
    } 
}*/

gslpp::complex THDMWcache::B00_MZ2_MZ2_mSr2_mSp2(const double MZ2, const double mSr2, const double mSp2) const {
    int NumPar = 3;
    double params[] = {MZ2, mSr2, mSp2};

    int i = CacheCheck(B00_MZ2_MZ2_mSr2_mSp2_cache, NumPar, params);
    if (i>=0) {
        return ( B00_MZ2_MZ2_mSr2_mSp2_cache[NumPar][i] );
    } else {
        gslpp::complex newResult = PV.B00(MZ2, MZ2 , mSr2, mSp2);
        CacheShift(B00_MZ2_MZ2_mSr2_mSp2_cache, NumPar, params, newResult);
        return newResult;
    } 
}

gslpp::complex THDMWcache::B00_MZ2_MZ2_mSi2_mSp2(const double MZ2, const double mSi2, const double mSp2) const {
    int NumPar = 3;
    double params[] = {MZ2, mSi2, mSp2};

    int i = CacheCheck(B00_MZ2_MZ2_mSi2_mSp2_cache, NumPar, params);
    if (i>=0) {
        return ( B00_MZ2_MZ2_mSi2_mSp2_cache[NumPar][i] );
    } else {
        gslpp::complex newResult = PV.B00(MZ2, MZ2 , mSi2, mSp2);
        CacheShift(B00_MZ2_MZ2_mSi2_mSp2_cache, NumPar, params, newResult);
        return newResult;
    } 
}

gslpp::complex THDMWcache::B00_MZ2_MZ2_mSr2_mSi2(const double MZ2, const double mSr2, const double mSi2) const {
    int NumPar = 3;
    double params[] = {MZ2, mSr2, mSi2};

    int i = CacheCheck(B00_MZ2_MZ2_mSr2_mSi2_cache, NumPar, params);
    if (i>=0) {
        return ( B00_MZ2_MZ2_mSr2_mSi2_cache[NumPar][i] );
    } else {
        gslpp::complex newResult = PV.B00(MZ2, MZ2 , mSr2, mSi2);
        CacheShift(B00_MZ2_MZ2_mSr2_mSi2_cache, NumPar, params, newResult);
        return newResult;
    } 
}



gslpp::complex THDMWcache::B00_MZ2_MZ2_mSp2_mSp2(const double MZ2, const double mSp2) const {
    int NumPar = 2;
    double params[] = {MZ2, mSp2};

    int i = CacheCheck(B00_MZ2_MZ2_mSp2_mSp2_cache, NumPar, params);
    if (i>=0) {
        return ( B00_MZ2_MZ2_mSp2_mSp2_cache[NumPar][i] );
    } else {
        gslpp::complex newResult = PV.B00(MZ2, MZ2 , mSp2, mSp2);
        CacheShift(B00_MZ2_MZ2_mSp2_mSp2_cache, NumPar, params, newResult);
        return newResult;
    } 
}











////////////////////////////////////////////////////
//--------- End Passarino Veltman Functions for THDMW ---------//
////////////////////////////////////////////////////




//double THDMWcache::ghHpHm(const double mHp2, const double tanb, const double m12_2, const double bma, const double mHl2, const double vev) const {
//    int NumPar = 6;
//    double params[] = {mHp2, tanb, m12_2, bma, mHl2, vev};
//
//    int i = CacheCheckReal(ghHpHm_cache, NumPar, params);
//    if (i>=0) {
//        return ( ghHpHm_cache[NumPar][i] );
//    } else {
//        double newResult = ((cos(bma)*mHl2*(tanb*tanb-1.0))/tanb 
//                                    -(mHl2+2.0*mHp2)*sin(bma) 
//                                    +(m12_2*(cos(bma)*(1.0-tanb*tanb)+2.0*sin(bma)*tanb)*(1.0+tanb*tanb))/(tanb*tanb))/vev;
//        CacheShiftReal(ghHpHm_cache, NumPar, params, newResult);
//        return newResult;
//    }
//}
//
//double THDMWcache::g_HH_HpHm(const double mHp2, const double mHh2, const double tanb, const double m12_2, const double bma, const double vev) const {
//    int NumPar = 6;
//    double params[] = {mHp2, mHh2, tanb, m12_2, bma, vev};
//
//    int i = CacheCheckReal(g_HH_HpHm_cache, NumPar, params);
//    if (i>=0) {
//        return ( g_HH_HpHm_cache[NumPar][i] );
//    } else {
//        double newResult = (cos(bma)*(mHh2-2.0*mHp2)
//                                    +((m12_2-mHh2*tanb+m12_2*tanb*tanb)
//                                      *(2.0*cos(bma)*tanb+sin(bma)*(tanb*tanb-1.0)))/(tanb*tanb))/vev;
//        CacheShiftReal(g_HH_HpHm_cache, NumPar, params, newResult);
//        return newResult;
//    }
//}

gslpp::complex THDMWcache::I_h_U(const double mHl2, const double Mu, const double Mc, const double Mt) const {
    int NumPar = 4;
    double params[] = {mHl2, Mu, Mc, Mt};

    int i = CacheCheck(I_h_U_cache, NumPar, params);
    if (i>=0) {
        return ( I_h_U_cache[NumPar][i] );
    } else {
    	double TAUu=4.0*Mu*Mu/mHl2;
    	double TAUc=4.0*Mc*Mc/mHl2;
    	double TAUt=4.0*Mt*Mt/mHl2;
        gslpp::complex newResult = -(8./3.)*(TAUu*(1.0+(1.0-TAUu)*f_func(TAUu))
                         +TAUc*(1.0+(1.0-TAUc)*f_func(TAUc))+TAUt*(1.0+(1.0-TAUt)*f_func(TAUt)));
        CacheShift(I_h_U_cache, NumPar, params, newResult);
        return newResult;
    }
}

gslpp::complex THDMWcache::I_HH_U(const double mHh2, const double Mc, const double Mt) const {
    int NumPar = 3;
    double params[] = {mHh2, Mc, Mt};

    int i = CacheCheck(I_HH_U_cache, NumPar, params);
    if (i>=0) {
        return ( I_HH_U_cache[NumPar][i] );
    } else {
    	double TAUc=4.0*Mc*Mc/mHh2;
    	double TAUt=4.0*Mt*Mt/mHh2;
        gslpp::complex newResult = -(8./3.)*(TAUc*(1.0+(1.0-TAUc)*f_func(TAUc))
                      +TAUt*(1.0+(1.0-TAUt)*f_func(TAUt)));
        CacheShift(I_HH_U_cache, NumPar, params, newResult);
        return newResult;
    }
}

gslpp::complex THDMWcache::I_A_U(const double mA2, const double Mc, const double Mt) const {
    int NumPar = 3;
    double params[] = {mA2, Mc, Mt};

    int i = CacheCheck(I_A_U_cache, NumPar, params);
    if (i>=0) {
        return ( I_A_U_cache[NumPar][i] );
    } else {
    	double TAUc=4.0*Mc*Mc/mA2;
    	double TAUt=4.0*Mt*Mt/mA2;
        gslpp::complex newResult = -(8./3.)*(TAUc*f_func(TAUc)+TAUt*f_func(TAUt));
        CacheShift(I_A_U_cache, NumPar, params, newResult);
        return newResult;
    }
}

gslpp::complex THDMWcache::I_h_D(const double mHl2, const double Md, const double Ms, const double Mb) const {
    int NumPar = 4;
    double params[] = {mHl2, Md, Ms, Mb};

    int i = CacheCheck(I_h_D_cache, NumPar, params);
    if (i>=0) {
        return ( I_h_D_cache[NumPar][i] );
    } else {
    	double TAUd=4.0*Md*Md/mHl2;
    	double TAUs=4.0*Ms*Ms/mHl2;
    	double TAUb=4.0*Mb*Mb/mHl2;
        gslpp::complex newResult = -(2./3.)*(TAUd*(1.0+(1.0-TAUd)*f_func(TAUd))
                         +TAUs*(1.0+(1.0-TAUs)*f_func(TAUs))+TAUb*(1.0+(1.0-TAUb)*f_func(TAUb)));
        CacheShift(I_h_D_cache, NumPar, params, newResult);
        return newResult;
    }
}

gslpp::complex THDMWcache::I_HH_D(const double mHh2, const double Ms, const double Mb) const {
    int NumPar = 3;
    double params[] = {mHh2, Ms, Mb};

    int i = CacheCheck(I_HH_D_cache, NumPar, params);
    if (i>=0) {
        return ( I_HH_D_cache[NumPar][i] );
    } else {
    	double TAUs=4.0*Ms*Ms/mHh2;
    	double TAUb=4.0*Mb*Mb/mHh2;
        gslpp::complex newResult = -(2./3.)*(TAUs*(1.0+(1.0-TAUs)*f_func(TAUs))
                      +TAUb*(1.0+(1.0-TAUb)*f_func(TAUb)));
        CacheShift(I_HH_D_cache, NumPar, params, newResult);
        return newResult;
    }
}

gslpp::complex THDMWcache::I_A_D(const double mA2, const double Ms, const double Mb) const {
    int NumPar = 3;
    double params[] = {mA2, Ms, Mb};

    int i = CacheCheck(I_A_D_cache, NumPar, params);
    if (i>=0) {
        return ( I_A_D_cache[NumPar][i] );
    } else {
    	double TAUs=4.0*Ms*Ms/mA2;
    	double TAUb=4.0*Mb*Mb/mA2;
        gslpp::complex newResult = -(2./3.)*(TAUs*f_func(TAUs)+TAUb*f_func(TAUb));
        CacheShift(I_A_D_cache, NumPar, params, newResult);
        return newResult;
    }
}

gslpp::complex THDMWcache::I_h_L(const double mHl2, const double Me, const double Mmu, const double Mtau) const {
    int NumPar = 4;
    double params[] = {mHl2, Me, Mmu, Mtau};

    int i = CacheCheck(I_h_L_cache, NumPar, params);
    if (i>=0) {
        return ( I_h_L_cache[NumPar][i] );
    } else {
    	double TAUe=4.0*Me*Me/mHl2;
    	double TAUmu=4.0*Mmu*Mmu/mHl2;
    	double TAUtau=4.0*Mtau*Mtau/mHl2;
        gslpp::complex newResult = -2.0*(TAUe*(1.0+(1.0-TAUe)*f_func(TAUe))
                         +TAUmu*(1.0+(1.0-TAUmu)*f_func(TAUmu))
                         +TAUtau*(1.0+(1.0-TAUtau)*f_func(TAUtau)));
        CacheShift(I_h_L_cache, NumPar, params, newResult);
        return newResult;
    }
}

gslpp::complex THDMWcache::I_HH_L(const double mHh2, const double Mmu, const double Mtau) const {
    int NumPar = 3;
    double params[] = {mHh2, Mmu, Mtau};

    int i = CacheCheck(I_HH_L_cache, NumPar, params);
    if (i>=0) {
        return ( I_HH_L_cache[NumPar][i] );
    } else {
    	double TAUmu=4.0*Mmu*Mmu/mHh2;
    	double TAUtau=4.0*Mtau*Mtau/mHh2;
        gslpp::complex newResult = -2.0*(TAUmu*(1.0+(1.0-TAUmu)*f_func(TAUmu))+
                           TAUtau*(1.0+(1.0-TAUtau)*f_func(TAUtau)));
        CacheShift(I_HH_L_cache, NumPar, params, newResult);
        return newResult;
    }
}

gslpp::complex THDMWcache::I_A_L(const double mA2, const double Mmu, const double Mtau) const {
    int NumPar = 3;
    double params[] = {mA2, Mmu, Mtau};

    int i = CacheCheck(I_A_L_cache, NumPar, params);
    if (i>=0) {
        return ( I_A_L_cache[NumPar][i] );
    } else {
    	double TAUmu=4.0*Mmu*Mmu/mA2;
    	double TAUtau=4.0*Mtau*Mtau/mA2;
        gslpp::complex newResult = -2.0*(TAUmu*f_func(TAUmu)+TAUtau*f_func(TAUtau));
        CacheShift(I_A_L_cache, NumPar, params, newResult);
        return newResult;
    }
}

gslpp::complex THDMWcache::I_H_W(const double mH, const double MW) const {
    int NumPar = 2;
    double params[] = {mH, MW};

    int i = CacheCheck(I_H_W_cache, NumPar, params);
    if (i>=0) {
        return ( I_H_W_cache[NumPar][i] );
    } else {
        double TAUw=4.0*MW*MW/(mH*mH);
        gslpp::complex newResult = 2.0 + 3.0*TAUw + 3.0*TAUw*(2.0-TAUw)*f_func(TAUw);
        CacheShift(I_H_W_cache, NumPar, params, newResult);
        return newResult;
    }
}

gslpp::complex THDMWcache::I_H_Hp(const double mHp2, const double mH) const {
    int NumPar = 2;
    double params[] = {mHp2, mH};

    int i = CacheCheck(I_H_Hp_cache, NumPar, params);
    if (i>=0) {
        return ( I_H_Hp_cache[NumPar][i] );
    } else {
        double TAUhp=4.0*mHp2/(mH*mH);
        gslpp::complex newResult = -TAUhp*(1.0-TAUhp*f_func(TAUhp));
        CacheShift(I_H_Hp_cache, NumPar, params, newResult);
        return newResult;
    }
}

gslpp::complex THDMWcache::A_h_U(const double mHl2, const double cW2, const double Mu, const double Mc, const double Mt, const double MZ) const {
    int NumPar = 6;
    double params[] = {mHl2, cW2, Mu, Mc, Mt, MZ};

    int i = CacheCheck(A_h_U_cache, NumPar, params);
    if (i>=0) {
        return ( A_h_U_cache[NumPar][i] );
    } else {
    	double TAUu=4.0*Mu*Mu/mHl2;
    	double TAUc=4.0*Mc*Mc/mHl2;
    	double TAUt=4.0*Mt*Mt/mHl2;
    	double LAMu=4.0*Mu*Mu/(MZ*MZ);
    	double LAMc=4.0*Mc*Mc/(MZ*MZ);
    	double LAMt=4.0*Mt*Mt/(MZ*MZ);
	double sW2=1.0-cW2;
        gslpp::complex newResult = -4.0*(1.0/2.0-4.0/3.0*sW2)*(Int1(TAUu,LAMu)+Int1(TAUc,LAMc)
                           +Int1(TAUt,LAMt)-Int2(TAUu,LAMu)-Int2(TAUc,LAMc)-Int2(TAUt,LAMt));
        CacheShift(A_h_U_cache, NumPar, params, newResult);
        return newResult;
    }
}

gslpp::complex THDMWcache::A_HH_U(const double mHh2, const double cW2, const double Mc, const double Mt, const double MZ) const {
    int NumPar = 5;
    double params[] = {mHh2, cW2, Mc, Mt, MZ};

    int i = CacheCheck(A_HH_U_cache, NumPar, params);
    if (i>=0) {
        return ( A_HH_U_cache[NumPar][i] );
    } else {
    	double TAUc=4.0*Mc*Mc/mHh2;
    	double TAUt=4.0*Mt*Mt/mHh2;
    	double LAMc=4.0*Mc*Mc/(MZ*MZ);
    	double LAMt=4.0*Mt*Mt/(MZ*MZ);
	double sW2=1.0-cW2;
        gslpp::complex newResult = -4.0*(1.0/2.0-4.0/3.0*sW2)*(Int1(TAUc,LAMc)-Int2(TAUc,LAMc)
                                         +Int1(TAUt,LAMt)-Int2(TAUt,LAMt));
        CacheShift(A_HH_U_cache, NumPar, params, newResult);
        return newResult;
    }
}

gslpp::complex THDMWcache::A_A_U(const double mA2, const double cW2, const double Mc, const double Mt, const double MZ) const {
    int NumPar = 5;
    double params[] = {mA2, cW2, Mc, Mt, MZ};

    int i = CacheCheck(A_A_U_cache, NumPar, params);
    if (i>=0) {
        return ( A_A_U_cache[NumPar][i] );
    } else {
    	double TAUc=4.0*Mc*Mc/mA2;
    	double TAUt=4.0*Mt*Mt/mA2;
    	double LAMc=4.0*Mc*Mc/(MZ*MZ);
    	double LAMt=4.0*Mt*Mt/(MZ*MZ);
	double sW2=1.0-cW2;
        gslpp::complex newResult = -4.0*(1.0/2.0-4.0/3.0*sW2)*(-Int2(TAUc,LAMc)-Int2(TAUt,LAMt))/sqrt(sW2*cW2);
        CacheShift(A_A_U_cache, NumPar, params, newResult);
        return newResult;
    }
}

gslpp::complex THDMWcache::A_h_D(const double mHl2, const double cW2, const double Md, const double Ms, const double Mb, const double MZ) const {
    int NumPar = 6;
    double params[] = {mHl2, cW2, Md, Ms, Mb, MZ};

    int i = CacheCheck(A_h_D_cache, NumPar, params);
    if (i>=0) {
        return ( A_h_D_cache[NumPar][i] );
    } else {
    	double TAUd=4.0*Md*Md/mHl2;
    	double TAUs=4.0*Ms*Ms/mHl2;
    	double TAUb=4.0*Mb*Mb/mHl2;
    	double LAMd=4.0*Md*Md/(MZ*MZ);
    	double LAMs=4.0*Ms*Ms/(MZ*MZ);
	double LAMb=4.0*Mb*Mb/(MZ*MZ);
	double sW2=1.0-cW2;
        gslpp::complex newResult = 2.0*(-1.0/2.0+2.0/3.0*sW2)*(Int1(TAUd,LAMd)+Int1(TAUs,LAMs)
                           +Int1(TAUb,LAMb)-Int2(TAUd,LAMd)-Int2(TAUs,LAMs)-Int2(TAUb,LAMb));
        CacheShift(A_h_D_cache, NumPar, params, newResult);
        return newResult;
    }
}

gslpp::complex THDMWcache::A_HH_D(const double mHh2, const double cW2, const double Ms, const double Mb, const double MZ) const {
    int NumPar = 5;
    double params[] = {mHh2, cW2, Ms, Mb, MZ};

    int i = CacheCheck(A_HH_D_cache, NumPar, params);
    if (i>=0) {
        return ( A_HH_D_cache[NumPar][i] );
    } else {
    	double TAUs=4.0*Ms*Ms/mHh2;
    	double TAUb=4.0*Mb*Mb/mHh2;
    	double LAMs=4.0*Ms*Ms/(MZ*MZ);
	double LAMb=4.0*Mb*Mb/(MZ*MZ);
	double sW2=1.0-cW2;
        gslpp::complex newResult = 2.0*(-1.0/2.0+2.0/3.0*sW2)*(Int1(TAUs,LAMs)-Int2(TAUs,LAMs)
                                          +Int1(TAUb,LAMb)-Int2(TAUb,LAMb));
        CacheShift(A_HH_D_cache, NumPar, params, newResult);
        return newResult;
    }
}

gslpp::complex THDMWcache::A_A_D(const double mA2, const double cW2, const double Ms, const double Mb, const double MZ) const {
    int NumPar = 5;
    double params[] = {mA2, cW2, Ms, Mb, MZ};

    int i = CacheCheck(A_A_D_cache, NumPar, params);
    if (i>=0) {
        return ( A_A_D_cache[NumPar][i] );
    } else {
    	double TAUs=4.0*Ms*Ms/mA2;
    	double TAUb=4.0*Mb*Mb/mA2;
    	double LAMs=4.0*Ms*Ms/(MZ*MZ);
	double LAMb=4.0*Mb*Mb/(MZ*MZ);
	double sW2=1.0-cW2;
        gslpp::complex newResult = 2.0*(-1.0/2.0+2.0/3.0*sW2)*(-Int2(TAUs,LAMs)-Int2(TAUb,LAMb))/sqrt(sW2*cW2);
        CacheShift(A_A_D_cache, NumPar, params, newResult);
        return newResult;
    }
}

gslpp::complex THDMWcache::A_h_L(const double mHl2, const double cW2, const double Me, const double Mmu, const double Mtau, const double MZ) const {
    int NumPar = 6;
    double params[] = {mHl2, cW2, Me, Mmu, Mtau, MZ};

    int i = CacheCheck(A_h_L_cache, NumPar, params);
    if (i>=0) {
        return ( A_h_L_cache[NumPar][i] );
    } else {
    	double TAUe=4.0*Me*Me/mHl2;
    	double TAUmu=4.0*Mmu*Mmu/mHl2;
    	double TAUtau=4.0*Mtau*Mtau/mHl2;
    	double LAMe=4.0*Me*Me/(MZ*MZ);
    	double LAMmu=4.0*Mmu*Mmu/(MZ*MZ);
	double LAMtau=4.0*Mtau*Mtau/(MZ*MZ);
	double sW2=1.0-cW2;
        gslpp::complex newResult = 2.0*(-1.0/2.0+2.0*sW2)*(Int1(TAUe,LAMe)+Int1(TAUmu,LAMmu)
                            +Int1(TAUtau,LAMtau)-Int2(TAUe,LAMe)-Int2(TAUmu,LAMmu)
                            -Int2(TAUtau,LAMtau));
        CacheShift(A_h_L_cache, NumPar, params, newResult);
        return newResult;
    }
}

gslpp::complex THDMWcache::A_HH_L(const double mHh2, const double cW2, const double Mmu, const double Mtau, const double MZ) const {
    int NumPar = 5;
    double params[] = {mHh2, cW2, Mmu, Mtau, MZ};

    int i = CacheCheck(A_HH_L_cache, NumPar, params);
    if (i>=0) {
        return ( A_HH_L_cache[NumPar][i] );
    } else {
    	double TAUmu=4.0*Mmu*Mmu/mHh2;
    	double TAUtau=4.0*Mtau*Mtau/mHh2;
    	double LAMmu=4.0*Mmu*Mmu/(MZ*MZ);
	double LAMtau=4.0*Mtau*Mtau/(MZ*MZ);
	double sW2=1.0-cW2;
        gslpp::complex newResult = 2.0*(-1.0/2.0+2.0*sW2)*(Int1(TAUmu,LAMmu)-Int2(TAUmu,LAMmu)
                                      +Int1(TAUtau,LAMtau)-Int2(TAUtau,LAMtau));
        CacheShift(A_HH_L_cache, NumPar, params, newResult);
        return newResult;
    }
}

gslpp::complex THDMWcache::A_A_L(const double mA2, const double cW2, const double Mmu, const double Mtau, const double MZ) const {
    int NumPar = 5;
    double params[] = {mA2, cW2, Mmu, Mtau, MZ};

    int i = CacheCheck(A_A_L_cache, NumPar, params);
    if (i>=0) {
        return ( A_A_L_cache[NumPar][i] );
    } else {
    	double TAUmu=4.0*Mmu*Mmu/mA2;
    	double TAUtau=4.0*Mtau*Mtau/mA2;
    	double LAMmu=4.0*Mmu*Mmu/(MZ*MZ);
	double LAMtau=4.0*Mtau*Mtau/(MZ*MZ);
	double sW2=1.0-cW2;
        gslpp::complex newResult = 2.0*(-1.0/2.0+2.0*sW2)*(-Int2(TAUmu,LAMmu)-Int2(TAUtau,LAMtau))/sqrt(sW2*cW2);
        CacheShift(A_A_L_cache, NumPar, params, newResult);
        return newResult;
    }
}

gslpp::complex THDMWcache::A_H_W(const double mH, const double cW2, const double MW, const double MZ) const {
    int NumPar = 4;
    double params[] = {mH, cW2, MW, MZ};

    int i = CacheCheck(A_H_W_cache, NumPar, params);
    if (i>=0) {
        return ( A_H_W_cache[NumPar][i] );
    } else {
        double TAUw=4.0*MW*MW/(mH*mH);
        double LAMw=4.0*MW*MW/(MZ*MZ);
	double sW2=1.0-cW2;
        gslpp::complex newResult = -sqrt(cW2/sW2)*(4.0*(3.0-sW2/cW2)*Int2(TAUw,LAMw)
                            +((1.0+2.0/TAUw)*sW2/cW2-(5.0+2.0/TAUw))*Int1(TAUw,LAMw));
        CacheShift(A_H_W_cache, NumPar, params, newResult);
        return newResult;
    }
}

gslpp::complex THDMWcache::A_H_Hp(const double mHp2, const double mH, const double cW2, const double MZ) const {
    int NumPar = 4;
    double params[] = {mHp2, mH, cW2, MZ};

    int i = CacheCheck(A_H_Hp_cache, NumPar, params);
    if (i>=0) {
        return ( A_H_Hp_cache[NumPar][i] );
    } else {
        double TAUhp=4.0*mHp2/(mH*mH);
        double LAMhp=4.0*mHp2/(MZ*MZ);
	double sW2=1.0-cW2;
        gslpp::complex newResult = (1.0-2.0*sW2)/sqrt(cW2*sW2)*Int1(TAUhp,LAMhp);
        CacheShift(A_H_Hp_cache, NumPar, params, newResult);
        return newResult;
    }
}

gslpp::complex THDMWcache::f_func(const double x) const{
    if(x<1) {
    gslpp::complex z = -gslpp::complex::i()*M_PI;
    return -pow(log((1.0+sqrt(1.0-x))/(1.0-sqrt(1.0-x)))+z,2)/4.0;
    }
    else {
        return pow(asin(sqrt(1.0/x)),2);
    }
}

gslpp::complex THDMWcache::g_func(const double x) const{
    if(x<1) {
    gslpp::complex z = -gslpp::complex::i()*M_PI;
    gslpp::complex gs1 = sqrt(1.0-x)*(log((1.0+sqrt(1.0-x))/(1.0-sqrt(1.0-x)))+z)/2.0;
    return gs1;
    }
    else {
        gslpp::complex gg1 = sqrt(x-1.0)*asin(sqrt(1.0/x));
        return gg1;
    }
}

gslpp::complex THDMWcache::Int1(const double tau, const double lambda) const{
    return tau*lambda/(tau-lambda)/2.0+tau*tau*lambda*lambda/((tau-lambda)
           *(tau-lambda))/2.0*(f_func(tau)-f_func(lambda))+tau*tau*lambda/((tau-lambda)
           *(tau-lambda))*(g_func(tau)-g_func(lambda));
}

gslpp::complex THDMWcache::Int2(const double tau, const double lambda) const{
    return -tau*lambda/(tau-lambda)/2.0*(f_func(tau)-f_func(lambda));
}

void THDMWcache::computeSignalStrengthQuantities()
{
    double Mt = myTHDMW->getQuarks(QCD::TOP).getMass();
    double Mb = myTHDMW->getQuarks(QCD::BOTTOM).getMass();
    double Mc = myTHDMW->getQuarks(QCD::CHARM).getMass();
    double Ms = myTHDMW->getQuarks(QCD::STRANGE).getMass();
    double Mu = myTHDMW->getQuarks(QCD::UP).getMass();
    double Md = myTHDMW->getQuarks(QCD::DOWN).getMass();
    double Mtau = myTHDMW->getLeptons(StandardModel::TAU).getMass();
    double Mmu = myTHDMW->getLeptons(StandardModel::MU).getMass();
    double Me = myTHDMW->getLeptons(StandardModel::ELECTRON).getMass();
    double MW = myTHDMW->Mw();
    double cW2 = myTHDMW->c02();
    double sW2=1.0-cW2;

    double BrSM_htobb = 5.77e-1;
    double BrSM_htotautau = 6.32e-2;
    double BrSM_htogaga = 2.28e-3;
    double BrSM_htoWW = 2.15e-1;
    double BrSM_htoZZ = 2.64e-2;
    double BrSM_htogg = 8.57e-2;
    double BrSM_htoZga = 1.54e-3;
    double BrSM_htocc = 2.91e-2;

    if( THDMWmodel == "ManoharWise" || THDMWmodel == "custodialMW" ) {
        rh_QuQu = 1.0;
        rh_VV = 1.0;
        rh_QdQd = 1.0;
        rh_ll = 1.0;

        //gluon coupling
        gslpp::complex fermU = I_h_U(mhsq,Mu,Mc,Mt);
        gslpp::complex fermD = I_h_D(mhsq,Md,Ms,Mb);
        double ch_p=-nu1*vev*vev/(4.0*mSpsq);//Victor Miralles, non-custodial
        double ch_r=-(nu1+nu2+2.0*nu3)*vev*vev/(4.0*mSRsq);//Victor Miralles, non-custodial
        double ch_i=-(nu1+nu2-2.0*nu3)*vev*vev/(4.0*mSIsq);//Victor Miralles, non-custodial
        gslpp::complex I_h_Sp = 4.5*ch_p*I_H_Hp(mSpsq,sqrt(mhsq));   //Factor 3 to normalize Higgs Hunters Guide to 1606.01298
        gslpp::complex I_h_SR = 2.25*ch_r*I_H_Hp(mSRsq,sqrt(mhsq));   //Factor 3 to normalize Higgs Hunters Guide to 1606.01298
        gslpp::complex I_h_SI = 2.25*ch_i*I_H_Hp(mSIsq,sqrt(mhsq));   //Factor 3 to normalize Higgs Hunters Guide to 1606.01298
        double ABSggMW=(-9.0/16.0*(fermU+4.0*fermD)+I_h_Sp+I_h_SR+I_h_SI).abs2();   //Factor 9/16 and 4 to normalize Higgs Hunters Guide to 1606.01298
        double ABSggSM=(-9.0/16.0*(fermU+4.0*fermD)).abs2();   //Factor 9/16 and 4 to normalize Higgs Hunters Guide to 1606.01298
        rh_gg=ABSggMW/ABSggSM;
        //std::cout<<"I_h_Sp*2.0 = "<<I_h_Sp*2.0<<std::endl;
        //std::cout<<"-9.0/16.0*fermU = "<<-9.0/16.0*fermU<<std::endl;
        //photon coupling
        gslpp::complex fermL = I_h_L(mhsq,Me,Mmu,Mtau);
        gslpp::complex I_hSM_W = I_H_W(sqrt(mhsq),MW);
        gslpp::complex I_h_S = -8.0*ch_p*I_H_Hp(mSpsq,sqrt(mhsq));//Factor of 1/3 cancels with the normalization factor 3 between Higgs Hunters Guide and 1606.01298
        double ABSgagaMW=(fermU+fermD+fermL+I_hSM_W+I_h_S).abs2();//-8 times (31) in 1606.01298
        double ABSgagaSM=(fermU+fermD+fermL+I_hSM_W).abs2();//-8 times (31) in 1606.01298
        rh_gaga=ABSgagaMW/ABSgagaSM;

        //Z photon coupling
        gslpp::complex A_h_Ux = A_h_U(mhsq,cW2,Mu,Mc,Mt,MZ);
        gslpp::complex A_h_Dx = A_h_D(mhsq,cW2,Md,Ms,Mb,MZ);
        gslpp::complex A_h_Lx  = A_h_L(mhsq,cW2,Me,Mmu,Mtau,MZ);
        gslpp::complex A_h_F = (A_h_Ux+A_h_Dx+A_h_Lx)/sqrt(sW2*cW2);
        gslpp::complex A_hSM_W = A_H_W(sqrt(mhsq),cW2,MW,MZ);
        gslpp::complex A_h_S = -8.0*ch_p*A_H_Hp(mSpsq,sqrt(mhsq),cW2,MZ);//Just analagous to the diphoton loop
        double ABSZgaMW=(A_h_F+A_hSM_W+A_h_S).abs2();
        double ABSZgaSM=(A_h_F+A_hSM_W).abs2();
        rh_Zga=ABSZgaMW/ABSZgaSM;

    }
    else if( THDMWmodel == "custodial1" ) {
        rh_QuQu = cosa*cosa/(sinb*sinb);
        rh_VV = sin(bma)*sin(bma);
        rh_QdQd = rh_QuQu;
        rh_ll = rh_QuQu;

        //gluon coupling
        gslpp::complex fermU = I_h_U(mhsq,Mu,Mc,Mt);
        gslpp::complex fermD = I_h_D(mhsq,Md,Ms,Mb);
        double ch_p=(-nu1*sina*cosb+omega1*cosa*sinb+kappa1*(cosa*cosb-sina*sinb))*vev*vev/(4.0*mSpsq);
        double ch_r=(-(nu1+2.0*nu2)*sina*cosb+(omega1+2.0*omega2)*cosa*sinb+(kappa1+2.0*kappa2)*(cosa*cosb-sina*sinb))*vev*vev/(4.0*mSRsq);
        double ch_i=ch_p;
        gslpp::complex I_h_Sp = 4.5*ch_p*I_H_Hp(mSpsq,sqrt(mhsq));   //Factor 3 to normalize Higgs Hunters Guide to 1606.01298
        gslpp::complex I_h_SR = 2.25*ch_r*I_H_Hp(mSRsq,sqrt(mhsq));   //Factor 3 to normalize Higgs Hunters Guide to 1606.01298
        gslpp::complex I_h_SI = 2.25*ch_i*I_H_Hp(mSIsq,sqrt(mhsq));   //Factor 3 to normalize Higgs Hunters Guide to 1606.01298
        double ABSggTHDMW=(-9.0/16.0*(cosa/sinb)*(fermU+4.0*fermD)+I_h_Sp+I_h_SR+I_h_SI).abs2();   //Factor 9/16 to normalize Higgs Hunters Guide to 1606.01298
        double ABSggSM=(-9.0/16.0*(fermU+4.0*fermD)).abs2();   //Factor 9/16 to normalize Higgs Hunters Guide to 1606.01298
        rh_gg=ABSggTHDMW/ABSggSM;

        //photon coupling
        double ghHpHm = vev*vev/mAsq * (-lambda1*sina*sinb*sinb*cosb+lambda2*cosa*sinb*cosb*cosb
                                        +lambda3*(cosa*sinb*sinb*sinb-sina*cosb*cosb*cosb)
                                        -2.0*lambda4*(cosa*cosb-sina*sinb)*sinb*cosb);
        gslpp::complex fermL = I_h_L(mhsq,Me,Mmu,Mtau);
        gslpp::complex I_hSM_W = I_H_W(sqrt(mhsq),MW);
        gslpp::complex I_h_Hp = -0.5*ghHpHm*I_H_Hp(mSpsq,sqrt(mhsq));
        gslpp::complex I_h_S = -8.0*ch_p*I_H_Hp(mSpsq,sqrt(mhsq));//Factor of 1/3 cancels with the normalization factor 3 between Higgs Hunters Guide and 1606.01298
        double ABSgagaTHDMW=((cosa/sinb)*(fermU+fermD+fermL)+sin(bma)*I_hSM_W+I_h_Hp+I_h_S).abs2();//-8 times (31) in 1606.01298
        double ABSgagaSM=(fermU+fermD+fermL+I_hSM_W).abs2();//-8 times (31) in 1606.01298
        rh_gaga=ABSgagaTHDMW/ABSgagaSM;

        //Z photon coupling
        gslpp::complex A_h_Ux = A_h_U(mhsq,cW2,Mu,Mc,Mt,MZ);
        gslpp::complex A_h_Dx = A_h_D(mhsq,cW2,Md,Ms,Mb,MZ);
        gslpp::complex A_h_Lx  = A_h_L(mhsq,cW2,Me,Mmu,Mtau,MZ);
        gslpp::complex A_h_F = cosa/sinb*(A_h_Ux+A_h_Dx+A_h_Lx)/sqrt(sW2*cW2);
        gslpp::complex A_hSM_W = A_H_W(sqrt(mhsq),cW2,MW,MZ);
        gslpp::complex A_h_Hp = -0.5*ghHpHm*A_H_Hp(mSpsq,sqrt(mhsq),cW2,MZ);
        gslpp::complex A_h_S = -8.0*ch_p*A_H_Hp(mSpsq,sqrt(mhsq),cW2,MZ);//Just analagous to the diphoton loop
        double ABSZgaTHDMW=(A_h_F+sin(bma)*A_hSM_W+A_h_Hp+A_h_S).abs2();
        double ABSZgaSM=(A_h_F+A_hSM_W).abs2();
        rh_Zga=ABSZgaTHDMW/ABSZgaSM;

    }
    else {
        throw std::runtime_error("THDMWmodel can only be \"ManoharWise\", \"custodialMW\" or \"custodial1\".");
    }

    sumModBRs = rh_QdQd*BrSM_htobb + rh_VV*(BrSM_htoWW+BrSM_htoZZ) + rh_ll*BrSM_htotautau +
          rh_gaga*BrSM_htogaga + rh_gg*BrSM_htogg + rh_Zga*BrSM_htoZga + rh_QuQu*BrSM_htocc;

    Gamma_h = sumModBRs*myTHDMW->computeGammaHTotal();
    
    THDM_BR_h_bb = rh_QdQd*BrSM_htobb/sumModBRs;
    THDM_BR_h_gaga = rh_gaga*BrSM_htogaga/sumModBRs;
    THDM_BR_h_tautau = rh_ll*BrSM_htotautau/sumModBRs;
    THDM_BR_h_WW = rh_VV*BrSM_htoWW/sumModBRs;
    THDM_BR_h_ZZ = rh_VV*BrSM_htoZZ/sumModBRs;
}

void THDMWcache::runTHDMWparameters()
{

    std::string RGEorder=myTHDMW->getRGEorderflag();
    //flag will be used to transport information about model and RGEorder to the Runner:
    //flag=0 for LO, 1 for approxNLO (and 2 for NLO - not implemented yet)
    int flag;
    if( RGEorder == "LO" ) flag=0;
    else if( RGEorder == "approxNLO" ) flag=1;
//    else if( RGEorder == "NLO" ) flag=2;
    else {
        throw std::runtime_error("RGEorder can be only any of \"LO\", \"approxNLO\" or \"NLO\"");
    }

    if( THDMWmodel == "custodial1")
    {
        double lambda1_at_MZ=lambda1;
        double lambda2_at_MZ=lambda2;
        double lambda3_at_MZ=lambda3;
        double lambda4_at_MZ=lambda4;
        double mu1_at_MZ=mu1;
        double mu3_at_MZ=mu3;
        double mu4_at_MZ=mu4;
        double nu1_at_MZ=nu1;
        double omega1_at_MZ=omega1;
        double kappa1_at_MZ=kappa1;
        double nu2_at_MZ=nu2;
        double omega2_at_MZ=omega2;
        double kappa2_at_MZ=kappa2;
        double nu4_at_MZ=nu4;
        double omega4_at_MZ=omega4;
        double NLOuniscale=myTHDMW->getNLOuniscaleTHDMW();

        if(fabs(Q_THDMW-log10(MZ))<0.005)   //at MZ scale
        {
            Q_cutoff=log10(MZ);

            lambda1_at_Q = lambda1_at_MZ;
            lambda2_at_Q = lambda2_at_MZ;
            lambda3_at_Q = lambda3_at_MZ;
            lambda4_at_Q = lambda4_at_MZ;
            mu1_at_Q = mu1_at_MZ;
            mu3_at_Q = mu3_at_MZ;
            mu4_at_Q = mu4_at_MZ;
            nu1_at_Q = nu1_at_MZ;
            omega1_at_Q = omega1_at_MZ;
            kappa1_at_Q = kappa1_at_MZ;
            nu2_at_Q = nu2_at_MZ;
            omega2_at_Q = omega2_at_MZ;
            kappa2_at_Q = kappa2_at_MZ;
            nu4_at_Q = nu4_at_MZ;
            omega4_at_Q = omega4_at_MZ;
        }
        else   //at some other scale
        {
            double InitVals[15];
            InitVals[0]=lambda1_at_MZ;
            InitVals[1]=lambda2_at_MZ;
            InitVals[2]=lambda3_at_MZ;
            InitVals[3]=lambda4_at_MZ;
            InitVals[4]=mu1_at_MZ;
            InitVals[5]=mu3_at_MZ;
            InitVals[6]=mu4_at_MZ;
            InitVals[7]=nu1_at_MZ;
            InitVals[8]=omega1_at_MZ;
            InitVals[9]=kappa1_at_MZ;
            InitVals[10]=nu2_at_MZ;
            InitVals[11]=omega2_at_MZ;
            InitVals[12]=kappa2_at_MZ;
            InitVals[13]=nu4_at_MZ;
            InitVals[14]=omega4_at_MZ;

            Q_cutoff=myRunnerTHDMW->RGERunnerTHDMW(InitVals, 15, log10(MZ), Q_THDMW, flag, RpepsTHDMW, NLOuniscale);  //Running up to Q_cutoff<=Q_THDM

            lambda1_at_Q = InitVals[0];
            lambda2_at_Q = InitVals[1];
            lambda3_at_Q = InitVals[2];
            lambda4_at_Q = InitVals[3];
            mu1_at_Q=InitVals[4];
            mu3_at_Q=InitVals[5];
            mu4_at_Q = InitVals[6];
            nu1_at_Q = InitVals[7];
            omega1_at_Q = InitVals[8];
            kappa1_at_Q = InitVals[9];
            nu2_at_Q = InitVals[10];
            omega2_at_Q = InitVals[11];
            kappa2_at_Q = InitVals[12];
            nu4_at_Q = InitVals[13];
            omega4_at_Q = InitVals[14];
        }
    }//End custodial1 case
    else if( THDMWmodel == "ManoharWise")
    {
        double lambda1_at_MZ=lambda1;
        double nu1_at_MZ=nu1;
        double nu2_at_MZ=nu2;
        double nu3_at_MZ=nu3;
        double nu4_at_MZ=nu4;
        double nu5_at_MZ=nu5;
        double mu1_at_MZ=mu1;
        double mu2_at_MZ=mu2;
        double mu3_at_MZ=mu3;
        double mu4_at_MZ=mu4;
        double mu5_at_MZ=mu5;
        double mu6_at_MZ=mu6;
        double NLOuniscale=myTHDMW->getNLOuniscaleTHDMW();

        if(fabs(Q_THDMW-log10(MZ))<0.005)   //at MZ scale
        {
            Q_cutoff=log10(MZ);

            lambda1_at_Q = lambda1_at_MZ;
            nu1_at_Q = nu1_at_MZ;
            nu2_at_Q = nu2_at_MZ;
            nu3_at_Q = nu3_at_MZ;
            nu4_at_Q = nu4_at_MZ;
            nu5_at_Q = nu5_at_MZ;
            mu1_at_Q = mu1_at_MZ;
            mu2_at_Q = mu2_at_MZ;
            mu3_at_Q = mu3_at_MZ;
            mu4_at_Q = mu4_at_MZ;
            mu5_at_Q = mu5_at_MZ;
            mu6_at_Q = mu6_at_MZ;
        }
        else   //at some other scale
        {
            double InitVals[12];
            InitVals[0]=lambda1_at_MZ;
            InitVals[1]=nu1_at_MZ;
            InitVals[2]=nu2_at_MZ;
            InitVals[3]=nu3_at_MZ;
            InitVals[4]=nu4_at_MZ;
            InitVals[5]=nu5_at_MZ;
            InitVals[6]=mu1_at_MZ;
            InitVals[7]=mu2_at_MZ;
            InitVals[8]=mu3_at_MZ;
            InitVals[9]=mu4_at_MZ;
            InitVals[10]=mu5_at_MZ;
            InitVals[11]=mu6_at_MZ;

            Q_cutoff=myRunnerTHDMW->RGERunnerMW(InitVals, 12, log10(MZ), Q_THDMW, flag, RpepsTHDMW, NLOuniscale);  //Running up to Q_cutoff<=Q_THDM

            lambda1_at_Q = InitVals[0];
            nu1_at_Q = InitVals[1];
            nu2_at_Q = InitVals[2];
            nu3_at_Q = InitVals[3];
            nu4_at_Q = InitVals[4];
            nu5_at_Q = InitVals[5];
            mu1_at_Q=InitVals[6];
            mu2_at_Q=InitVals[7];
            mu3_at_Q=InitVals[8];
            mu4_at_Q = InitVals[9];
            mu5_at_Q=InitVals[10];
            mu6_at_Q=InitVals[11];
        }
    }//End ManoharWise case
    else if( THDMWmodel == "custodialMW")
    {
        double lambda1_at_MZ=lambda1;
        double nu1_at_MZ=nu1;
        double nu2_at_MZ=nu2;
        double nu4_at_MZ=nu4;
        double mu1_at_MZ=mu1;
        double mu3_at_MZ=mu3;
        double mu4_at_MZ=mu4;
        double NLOuniscale=myTHDMW->getNLOuniscaleTHDMW();

        if(fabs(Q_THDMW-log10(MZ))<0.005)   //at MZ scale
        {
            Q_cutoff=log10(MZ);

            lambda1_at_Q = lambda1_at_MZ;
            nu1_at_Q = nu1_at_MZ;
            nu2_at_Q = nu2_at_MZ;
            nu4_at_Q = nu4_at_MZ;
            mu1_at_Q = mu1_at_MZ;
            mu3_at_Q = mu3_at_MZ;
            mu4_at_Q = mu4_at_MZ;
        }
        else   //at some other scale
        {
            double InitVals[12];
            InitVals[0]=lambda1_at_MZ;
            InitVals[1]=nu1_at_MZ;
            InitVals[2]=nu2_at_MZ;
            InitVals[3]=nu4_at_MZ;
            InitVals[4]=mu1_at_MZ;
            InitVals[5]=mu3_at_MZ;
            InitVals[6]=mu4_at_MZ;

            Q_cutoff=myRunnerTHDMW->RGERunnercustodialMW(InitVals, 7, log10(MZ), Q_THDMW, flag, RpepsTHDMW, NLOuniscale);  //Running up to Q_cutoff<=Q_THDM

            lambda1_at_Q = InitVals[0];
            nu1_at_Q = InitVals[1];
            nu2_at_Q = InitVals[2];
            nu4_at_Q = InitVals[3];
            mu1_at_Q=InitVals[4];
            mu3_at_Q=InitVals[5];
            mu4_at_Q = InitVals[6];
        }
    }//End custodialMW case
}

void THDMWcache::computeUnitarity()
{
    if( THDMWmodel != "custodial1" && THDMWmodel != "ManoharWise" && THDMWmodel != "custodialMW")
    {
        throw std::runtime_error("THDMW unitarity constraints are only implemented for the \"custodial1\", the \"ManoharWise\" and the \"custodialMW\" model.");
    }

    if( THDMWmodel == "custodial1")
    {
        double pi=M_PI;
        gslpp::matrix<gslpp::complex> Smatrix1(4,4,0.), Smatrix2(4,4,0.);
        gslpp::matrix<gslpp::complex> Sbmatrix1(4,4,0.), Sbmatrix2(4,4,0.);
        gslpp::matrix<gslpp::complex> Seigenvectors1(4,4,0.), Seigenvectors2(4,4,0.);
        gslpp::matrix<gslpp::complex> Seigenvectors1T(4,4,0.), Seigenvectors2T(4,4,0.);
        gslpp::vector<double> Seigenvalues1(4,0.), Seigenvalues2(4,0.);
        gslpp::vector<gslpp::complex> Sbeigenvalues1(4,0.), Sbeigenvalues2(4,0.);

        /*
        *******   LO part   *************
        */

        // Definition of the blocks of the S-matrix
        Smatrix1.assign(0,0, 3.0*lambda1/(16.0*pi));
        Smatrix1.assign(0,1, (2.0*lambda3+lambda4)/(16.0*pi));
        Smatrix1.assign(1,0, Smatrix1(0,1));
        Smatrix1.assign(0,3, (2.0*nu1+nu2)/(8.0*sqrt(2.0)*pi));
        Smatrix1.assign(3,0, Smatrix1(0,3));
        Smatrix1.assign(1,1, 3.0*lambda2/(16.0*pi));
        Smatrix1.assign(1,3, (2.0*omega1+omega2)/(8.0*sqrt(2.0)*pi));
        Smatrix1.assign(3,1, Smatrix1(1,3));
        Smatrix1.assign(2,2, (lambda3+5.0*lambda4)/(16.0*pi));
        Smatrix1.assign(2,3, (4.0*kappa1+2.0*kappa2)/(16.0*pi));
        Smatrix1.assign(3,2, Smatrix1(2,3));
        Smatrix1.assign(3,3, (26.0*mu1+17.0*mu3+13.0*mu4)/(32.0*pi));

        Smatrix2.assign(0,0, lambda1/(16.0*pi));
        Smatrix2.assign(0,1, lambda4/(16.0*pi));
        Smatrix2.assign(1,0, Smatrix2(0,1));
        Smatrix2.assign(0,3, nu2/(8.0*sqrt(2.0)*pi));
        Smatrix2.assign(3,0, Smatrix2(0,3));
        Smatrix2.assign(1,1, lambda2/(16.0*pi));
        Smatrix2.assign(1,3, omega2/(8.0*sqrt(2.0)*pi));
        Smatrix2.assign(3,1, Smatrix2(1,3));
        Smatrix2.assign(2,2, (lambda3+lambda4)/(16.0*pi));
        Smatrix2.assign(2,3, kappa2/(8.0*pi));
        Smatrix2.assign(3,2, Smatrix2(2,3));
        Smatrix2.assign(3,3, (14.0*mu1+3.0*mu3+27.0*mu4)/(96.0*pi));

        Smatrix1.eigensystem(Seigenvectors1, Seigenvalues1);
        Smatrix2.eigensystem(Seigenvectors2, Seigenvalues2);

        for (int i=0; i < 4; i++) {
            unitarityeigenvalues.assign(i, Seigenvalues1(i));
            unitarityeigenvalues.assign(4+i, Seigenvalues2(i));
        }
        unitarityeigenvalues.assign(8, (lambda3-lambda4)/(16.0*pi));
        unitarityeigenvalues.assign(9, sqrt(15.0)*nu4/(16.0*pi));
        unitarityeigenvalues.assign(10, sqrt(15.0)*omega4/(16.0*pi));

        /*
        *******   NLO part   *************
        */

        double blambda1=(12.0*lambda1*lambda1 + 4.0*lambda3*lambda3 + 4.0*lambda3*lambda4 + 4.0*lambda4*lambda4 
                         + 8.0*nu1*nu1 + 8.0*nu1*nu2 + 8.0*nu2*nu2)/(16.0*pi*pi);
        double blambda2=(12.0*lambda2*lambda2 + 4.0*lambda3*lambda3 + 4.0*lambda3*lambda4 + 4.0*lambda4*lambda4
                         + 8.0*omega1*omega1 + 8.0*omega1*omega2 + 8.0*omega2*omega2)/(16.0*pi*pi);
        double blambda3=(4.0*lambda3*lambda3 + 4.0*lambda4*lambda4 + (lambda1+lambda2)*(6.0*lambda3+2.0*lambda4) 
                         + 8.0*kappa2*kappa2 + 8.0*nu1*omega1 + 4.0*nu2*omega1 + 4.0*nu1*omega2)/(16.0*pi*pi);
        double blambda4=(lambda1*lambda4 + lambda2*lambda4 + 4.0*lambda3*lambda4 + 6.0*lambda4*lambda4
                         + 4.0*kappa1*kappa1 + 4.0*kappa1*kappa2 + 2.0*kappa2*kappa2 + 2.0*nu2*omega2)/(8.0*pi*pi);
        double bmu1=(11.0*mu1*mu1 + 3.0*mu1*mu4 + mu1*(2.0*mu1+6.0*mu3+3.0*mu4)
                   + 3.0*nu4*nu4 + 3.0*omega4*omega4)/(16.0*pi*pi);
        double bmu3=(18.0*kappa1*kappa1 + 18.0*kappa1*kappa2 + 134.0*mu1*mu1 + 6.0*mu1*(39.0*mu3 + 22.0*mu4)
                   + 3.0*(30.0*mu3*mu3 + 39.0*mu3*mu4 + 9.0*mu4*mu4 
                          + 3.0*nu1*nu1 + 3.0*nu1*nu2 - 5.0*nu4*nu4
                          + 3.0*omega1*omega1 + 3.0*omega1*omega2 - 5.0*omega4*omega4))/(72.0*pi*pi);
        double bmu4=(18.0*kappa2*kappa2 + 4.0*mu1*mu1 + 156.0*mu1*mu4 + 54.0*mu3*mu4 + 144.0*mu4*mu4
                   + 9.0*nu2*nu2 + 6.0*nu4*nu4 + 9.0*omega2*omega2 + 6.0*omega4*omega4)/(144.0*pi*pi);
        double bnu1=(6.0*kappa1*kappa1 + 6.0*kappa2*kappa2 + 18.0*lambda1*nu1
                   + 78.0*mu1*nu1 + 51.0*mu3*nu1 + 39.0*mu4*nu1 + 6.0*nu1*nu1
                   + 6.0*lambda1*nu2 + 32.0*mu1*nu2 + 24.0*mu3*nu2 + 6.0*mu4*nu2
                   + 6.0*nu2*nu2 + 10.0*nu4*nu4
                   + 12.0*lambda3*omega1 + 6.0*lambda4*omega1 + 6.0*lambda3*omega2)/(48.0*pi*pi);
        double bomega1=(6.0*kappa1*kappa1 + 6.0*kappa2*kappa2 
                   + 12.0*lambda3*nu1 + 6.0*lambda4*nu1 + 6.0*lambda3*nu2
                   + 18.0*lambda2*omega1 + 78.0*mu1*omega1 + 51.0*mu3*omega1 + 39.0*mu4*omega1 + 6.0*omega1*omega1
                   + 6.0*lambda2*omega2 + 32.0*mu1*omega2 + 24.0*mu3*omega2 + 6.0*mu4*omega2 + 6.0*omega2*omega2
                   + 10.0*omega4*omega4)/(48.0*pi*pi);
        double bkappa1=(6.0*kappa1*(2.0*lambda3 + 10.0*lambda4 + 18.0*mu1 + 17.0*mu3 + 13.0*mu4 + 2.0*nu1 + 2.0*omega1)
                   + kappa2*(24.0*lambda4 + 64.0*mu1 + 48.0*mu3 + 24.0*mu4 + 9.0*nu2 + 9.0*omega2)
                   + 20.0*nu4*omega4)/(96.0*pi*pi);
        double bnu2=(4.0*kappa1*kappa2 + 6.0*kappa2*kappa2 + 2.0*lambda1*nu2 + ((14.0*mu1)/3.0 + mu3 + 9.0*mu4)*nu2 
                    + 4.0*nu1*nu2 + 6.0*nu2*nu2 + (25.0*nu4*nu4)/3.0 + 2.0*lambda4*omega2)/(16.0*pi*pi);
        double bomega2=(4.0*kappa1*kappa2 + 6.0*kappa2*kappa2 + 2.0*lambda4*nu2 + 2.0*lambda2*omega2 
                    + ((14.0*mu1)/3.0 + mu3 + 9.0*mu4)*omega2 + 4.0*omega1*omega2 + 6.0*omega2*omega2 
                    + (25.0*omega4*omega4)/3.0)/(16.0*pi*pi);
        double bkappa2=(kappa2*(6.0*lambda3 + 6.0*lambda4 + 14.0*mu1 + 3.0*mu3 + 27.0*mu4
                         + 6.0*nu1 + 12.0*nu2 + 6.0*omega1 + 12.0*omega2)
                    + 6.0*kappa1*(nu2 + omega2) + 42.0*nu4*omega4)/(48.0*pi*pi);
        double bnu4=(11.0*mu1*nu4 + 3.0*mu3*nu4 + 9.0*mu4*nu4 + 3.0*nu1*nu4 + 9.0*nu2*nu4 
                    + 3.0*kappa1*omega4 + 9.0*kappa2*omega4)/(16.0*pi*pi);
        double bomega4=(3.0*kappa1*nu4 + 9.0*kappa2*nu4 
                    + (11.0*mu1 + 3.0*(mu3 + 3.0*mu4 + omega1 + 3.0*omega2))*omega4)/(16.0*pi*pi);

        Sbmatrix1.assign(0,0, 3.0*blambda1/(16.0*pi));
        Sbmatrix1.assign(0,1, (2.0*blambda3+blambda4)/(16.0*pi));
        Sbmatrix1.assign(1,0, Sbmatrix1(0,1));
        Sbmatrix1.assign(0,3, (2.0*bnu1+bnu2)/(8.0*sqrt(2.0)*pi));
        Sbmatrix1.assign(3,0, Sbmatrix1(0,3));
        Sbmatrix1.assign(1,1, 3.0*blambda2/(16.0*pi));
        Sbmatrix1.assign(1,3, (2.0*bomega1+bomega2)/(8.0*sqrt(2.0)*pi));
        Sbmatrix1.assign(3,1, Sbmatrix1(1,3));
        Sbmatrix1.assign(2,2, (blambda3+5.0*blambda4)/(16.0*pi));
        Sbmatrix1.assign(2,3, (4.0*bkappa1+2.0*bkappa2)/(16.0*pi));
        Sbmatrix1.assign(3,2, Sbmatrix1(2,3));
        Sbmatrix1.assign(3,3, (26.0*bmu1+17.0*bmu3+13.0*bmu4)/(32.0*pi));

        Sbmatrix2.assign(0,0, blambda1/(16.0*pi));
        Sbmatrix2.assign(0,1, blambda4/(16.0*pi));
        Sbmatrix2.assign(1,0, Sbmatrix2(0,1));
        Sbmatrix2.assign(0,3, bnu2/(8.0*sqrt(2.0)*pi));
        Sbmatrix2.assign(3,0, Sbmatrix2(0,3));
        Sbmatrix2.assign(1,1, blambda2/(16.0*pi));
        Sbmatrix2.assign(1,3, bomega2/(8.0*sqrt(2.0)*pi));
        Sbmatrix2.assign(3,1, Sbmatrix2(1,3));
        Sbmatrix2.assign(2,2, (blambda3+blambda4)/(16.0*pi));
        Sbmatrix2.assign(2,3, bkappa2/(8.0*pi));
        Sbmatrix2.assign(3,2, Sbmatrix2(2,3));
        Sbmatrix2.assign(3,3, (14.0*bmu1+3.0*bmu3+27.0*bmu4)/(96.0*pi));

        Seigenvectors1T=Seigenvectors1.hconjugate();
        Seigenvectors2T=Seigenvectors2.hconjugate();

        for (int i=0; i < 4; i++) {
            for (int k=0; k < 4; k++) {
                for (int l=0; l < 4; l++) {
                    Sbeigenvalues1.assign(i, Sbeigenvalues1(i) + Seigenvectors1T(i,k) * Sbmatrix1(k,l) * Seigenvectors1(l,i) );
                    Sbeigenvalues2.assign(i, Sbeigenvalues2(i) + Seigenvectors2T(i,k) * Sbmatrix2(k,l) * Seigenvectors2(l,i) );
                }                
            }
            betaeigenvalues.assign(i, -1.5 * Sbeigenvalues1(i));
            betaeigenvalues.assign(i+4, -1.5 * Sbeigenvalues2(i));
        }

        betaeigenvalues.assign(8, -1.5 * (blambda3-blambda4)/(16.0*pi));
        betaeigenvalues.assign(9, -1.5 * sqrt(15.0)*bnu4/(16.0*pi));
        betaeigenvalues.assign(10, -1.5 * sqrt(15.0)*bomega4/(16.0*pi));

        for (int i=0; i < 11; i++) {
            NLOunitarityeigenvalues.assign(i, -(gslpp::complex::i()-1.0/pi)*unitarityeigenvalues(i)*unitarityeigenvalues(i) + betaeigenvalues(i) );
        }
    }//End of the custodial1 case

    if( THDMWmodel == "ManoharWise")
    {
        double pi=M_PI;

        /*
        *******   LO part   *************
        */

        // Eigenvalues of the S-matrix, calculated by Li Cheng and Victor Miralles
        double muA = 4.0*mu1+4.0*mu2+8.5*mu3+5.0*mu4+1.5*mu5+2.5*mu6;
        double muB = (4.0*mu1+4.0*mu2+1.5*mu3+12.0*mu4+1.5*mu5-0.5*mu6)/3.0;
        double muC = (-0.5*mu1-0.5*mu2+1.5*mu3+1.5*mu4+12.0*mu5+4.0*mu6)/3.0;
        double MA1 = 3.0*lambda1 + muA - sqrt(9.0*lambda1*lambda1-6.0*lambda1*muA+muA*muA+32.0*nu1*nu1+32.0*nu1*nu2+8.0*nu2*nu2);
        double MA2 = 3.0*lambda1 + muA + sqrt(9.0*lambda1*lambda1-6.0*lambda1*muA+muA*muA+32.0*nu1*nu1+32.0*nu1*nu2+8.0*nu2*nu2);
        double MB1 = lambda1 + muB - sqrt(lambda1*lambda1-2.0*lambda1*muB+muB*muB+8.0*nu2*nu2);
        double MB2 = lambda1 + muB + sqrt(lambda1*lambda1-2.0*lambda1*muB+muB*muB+8.0*nu2*nu2);
        double MC1 = lambda1 + muC - sqrt(lambda1*lambda1-2.0*lambda1*muC+muC*muC+32.0*nu3*nu3);
        double MC2 = lambda1 + muC + sqrt(lambda1*lambda1-2.0*lambda1*muC+muC*muC+32.0*nu3*nu3);
        unitarityeigenvalues.assign(0, MA1/(32.0*pi));
        unitarityeigenvalues.assign(1, MA2/(32.0*pi));
        unitarityeigenvalues.assign(2, MB1/(32.0*pi));
        unitarityeigenvalues.assign(3, MB2/(32.0*pi));
        unitarityeigenvalues.assign(4, MC1/(32.0*pi));
        unitarityeigenvalues.assign(5, MC2/(32.0*pi));
        unitarityeigenvalues.assign(6, lambda1/(16.0*pi));
        unitarityeigenvalues.assign(7, sqrt(15.0)*(nu4+nu5)/(64.0*pi));
        
        /*
        *******   NLO part   *************
        */
        
            //beta_lambda1
        double betalambda1 = (12.0*lambda1*lambda1 + 8.0*nu1*nu1 + 8.0*nu1*nu2 + 4.0*nu2*nu2 + 16.0*nu3*nu3)/(16.0*pi*pi);
        //beta_nu1
        double betanu1 = (2.0*nu1*nu1 + nu2*nu2 + 4.0*nu3*nu3 + 2.0*lambda1*(3.0*nu1+nu2)
               + (7.0*nu4*nu4 - 4.0*nu4*nu5 + 7.0*nu5*nu5)/3.0
               + nu1*(8.0*mu1 + 8.0*mu2 + 17.0*mu3 + 10.0*mu4 + 3.0*mu5 + 5.0*mu6)
               + nu2*(8.0*mu1 + 8.0*mu2 + 24.0*mu3 + 3.0*mu4
                      + 3.0*mu5 + 8.0*mu6)/3.0)/(16.0*pi*pi);
        //beta_nu2
        double betanu2 = (2.0*nu2*nu2 + 4.0*nu1*nu2 + 16.0*nu3*nu3 + 2.0*lambda1*nu2
               + (4.0*nu4*nu4 + 17.0*nu4*nu5 + 4.0*nu5*nu5)/3.0
               + nu2*(8.0*mu1 + 8.0*mu2 + 3.0*mu3 + 24.0*mu4
                      + 3.0*mu5 - mu6)/3.0)/(16.0*pi*pi);
        //beta_nu3
        double betanu3 = (2.0*nu3*(lambda1 + 2.0*nu1 + 3.0*nu2)
               + (17.0*nu4*nu4 + 16.0*nu4*nu5 + 17.0*nu5*nu5)/12.0
               + nu3*(-mu1 - mu2 + 3.0*mu3 + 3.0*mu4
                      + 24.0*mu5 + 8.0*mu6)/3.0)/(16.0*pi*pi);
        //beta_nu4
        double betanu4 = (8.0*nu3*nu4 + 2.0*nu3*nu5
               + nu5*(2.0*nu2 - mu2 + 2.0*mu4 + 4.0*mu5 + mu6)
               + nu4*(3.0*nu1 + 2.0*nu2 + 6.0*mu1 + 2.0*mu2 + 3.0*mu3
                      + 2.0*mu4 + mu5 + mu6))/(16.0*pi*pi);
        //beta_nu5
        double betanu5 = (2.0*nu3*nu4 + 8.0*nu3*nu5
               + nu4*(2.0*nu2 - mu1 + 2.0*mu4 + 4.0*mu5 + mu6)
               + nu5*(3.0*nu1 + 2.0*nu2 + 6.0*mu1 + 2.0*mu2 + 3.0*mu3
                      + 2.0*mu4 + mu5 + mu6))/(16.0*pi*pi);
        //beta_mu1
        double betamu1 = (3.0*nu4*nu4 + 7.0*mu1*mu1
               + mu1*(6.0*mu2 + 6.0*mu3 + 4.0*mu4 - mu5 - 2.0*mu6)
               + mu2*(4.0*mu4 - mu5)
               - 2.0*mu4*mu6 + 2.0*mu5*mu6 + mu6*mu6)/(16.0*pi*pi);
        //beta_mu2
        double betamu2 = (3.0*nu5*nu5 + 7.0*mu2*mu2
               + mu2*(6.0*mu1 + 6.0*mu3 + 4.0*mu4 - mu5 - 2.0*mu6)
               + mu1*(4.0*mu4 - mu5)
               - 2.0*mu4*mu6 + 2.0*mu5*mu6 + mu6*mu6)/(16.0*pi*pi);
        //beta_mu3
        double betamu3 = (20.0*mu3*mu3
               + mu3*(288.0*mu1 + 288.0*mu2 + 360.0*mu4 + 108.0*mu5 + 180.0*mu6)/18.0
               + (36.0*nu1*nu1 + 36.0*nu1*nu2 - 24.0*nu4*nu4 - 12.0*nu4*nu5
                  - 24.0*nu5*nu5 + 62.0*mu1*mu1 + 64.0*mu1*mu2 + 62.0*mu2*mu2
                  + (96.0*mu4 + 18.0*mu5 + 58.0*mu6)*(mu1 + mu2)
                  + 54.0*mu4*mu4 + 36.0*mu4*mu5 + 132.0*mu4*mu6 + 18.0*mu5*mu5
                  + 18.0*mu5*mu6 + 29.0*mu6*mu6)/18.0)/(16.0*pi*pi);
        //beta_mu4
        double betamu4 = (nu2*nu2 - (nu4*nu4 - 4.0*nu4*nu5 + nu5*nu5)/3.0 + 10.0*mu4*mu4 /*mu4??*/
               + mu5*(mu1 + mu2 + mu6)
               + mu4*(4.0*(4.0*mu1 + 4.0*mu2 + mu6)/3.0 + 2.0*mu5 + 6.0*mu4)
               + 4.0*mu5*mu5
               + (mu1*mu1 + mu2*mu2 - 4.0*mu6*(mu1+mu2) - 2.0*mu6*mu6)/9.0
               + 26.0/9.0*mu1*mu2)/(16.0*pi*pi);
        //beta_mu5
        double betamu5 = (4.0*nu3*nu3 - (nu4*nu4 - 4.0*nu4*nu5 + nu5*nu5)/3.0
               + mu5*((mu1 + mu2 + 19.0*mu6)/3.0 + 8.0*mu4 + 6.0*mu3)
               + 2.0*mu4*mu6 + 8.0*mu5*mu5
               + (mu1*mu1 + mu2*mu2 - 4.0*mu6*(mu1+mu2) + 7.0*mu6*mu6)/9.0
               - 10.0/9.0*mu1*mu2)/(16.0*pi*pi);
        //beta_mu6
        double betamu6 = (0.5*mu6*mu6 + 3.0*nu4*nu4 + 3.0*nu5*nu5
               - 2.0*(mu1*mu1 + mu2*mu2) + 6.0*mu5*(mu1 + mu2)
               + 7.0*mu6*(mu1 + mu2 + mu3))/(16.0*pi*pi);
        
        
        double betamuA = 4.0*betamu1+4.0*betamu2+8.5*betamu3+5.0*betamu4+1.5*betamu5+2.5*betamu6;
        double betamuB = (4.0*betamu1+4.0*betamu2+1.5*betamu3+12.0*betamu4+1.5*betamu5-0.5*betamu6)/3.0;
        double betamuC = (-0.5*betamu1-0.5*betamu2+1.5*betamu3+1.5*betamu4+12.0*betamu5+4.0*betamu6)/3.0;
        double betaMA1 = 3.0*betalambda1 + betamuA 
                         - sqrt(9.0*betalambda1*betalambda1-6.0*betalambda1*betamuA+betamuA*betamuA
                                +32.0*betanu1*betanu1+32.0*betanu1*betanu2+8.0*betanu2*betanu2);
        double betaMA2 = 3.0*betalambda1 + betamuA 
                         + sqrt(9.0*betalambda1*betalambda1-6.0*betalambda1*betamuA+betamuA*betamuA
                                +32.0*betanu1*betanu1+32.0*betanu1*betanu2+8.0*betanu2*betanu2);
        double betaMB1 = betalambda1 + betamuB - sqrt(betalambda1*betalambda1-2.0*betalambda1*betamuB+betamuB*betamuB+8.0*betanu2*betanu2);
        double betaMB2 = betalambda1 + betamuB + sqrt(betalambda1*betalambda1-2.0*betalambda1*betamuB+betamuB*betamuB+8.0*betanu2*betanu2);
        double betaMC1 = betalambda1 + betamuC - sqrt(betalambda1*betalambda1-2.0*betalambda1*betamuC+betamuC*betamuC+32.0*betanu3*betanu3);
        double betaMC2 = betalambda1 + betamuC + sqrt(betalambda1*betalambda1-2.0*betalambda1*betamuC+betamuC*betamuC+32.0*betanu3*betanu3);

        betaeigenvalues.assign(0, -1.5 * betaMA1/(32.0*pi));
        betaeigenvalues.assign(1, -1.5 * betaMA2/(32.0*pi));
        betaeigenvalues.assign(2, -1.5 * betaMB1/(32.0*pi));
        betaeigenvalues.assign(3, -1.5 * betaMB2/(32.0*pi));
        betaeigenvalues.assign(4, -1.5 * betaMC1/(32.0*pi));
        betaeigenvalues.assign(5, -1.5 * betaMC2/(32.0*pi));
        betaeigenvalues.assign(6, -1.5 * betalambda1/(16.0*pi));
        betaeigenvalues.assign(7, -1.5 * sqrt(15.0)*(betanu4+betanu5)/(64.0*pi));
        
        
        for (int i=0; i < 8; i++) {
            NLOunitarityeigenvalues.assign(i, -(gslpp::complex::i()-1.0/pi)*unitarityeigenvalues(i)*unitarityeigenvalues(i) + betaeigenvalues(i) );
        }
        
    }//End of the ManoharWise case
    if( THDMWmodel == "custodialMW")
    {
        double pi=M_PI;
        gslpp::matrix<gslpp::complex> Smatrix1(2,2,0.), Smatrix2(2,2,0.);
        gslpp::matrix<gslpp::complex> Seigenvectors1(2,2,0.), Seigenvectors2(2,2,0.);
        gslpp::vector<double> Seigenvalues1(2,0.), Seigenvalues2(2,0.);

        /*
        *******   LO part   *************
        */

        // Definition of the blocks of the S-matrix, taken from 1303.4848
        Smatrix1.assign(0,0, 3.0*lambda1/(16.0*pi));
        Smatrix1.assign(0,1, (2.0*nu1+nu2)/(8.0*sqrt(2.0)*pi));
        Smatrix1.assign(1,0, Smatrix1(0,1));
        Smatrix1.assign(1,1, (26.0*mu1+17.0*mu3+13.0*mu4)/(32.0*pi));

        Smatrix2.assign(0,0, lambda1/(16.0*pi));
        Smatrix2.assign(0,1, nu2/(8.0*sqrt(2.0)*pi));
        Smatrix2.assign(1,0, Smatrix2(0,1));
        Smatrix2.assign(1,1, (14.0*mu1+3.0*mu3+27.0*mu4)/(96.0*pi));

        Smatrix1.eigensystem(Seigenvectors1, Seigenvalues1);
        Smatrix2.eigensystem(Seigenvectors2, Seigenvalues2);

        for (int i=0; i < 2; i++) {
            unitarityeigenvalues.assign(i, Seigenvalues1(i));
            unitarityeigenvalues.assign(2+i, Seigenvalues2(i));
        }
        unitarityeigenvalues.assign(4, sqrt(15.0)*nu4/(16.0*pi)); //non-custodial limit from 1606.01298
    }//End of the custodialMW case
}

// Direct Searches



gslpp::matrix<double> THDMWcache::readTable(std::string filename, int rowN, int colN){

    std::ifstream INfile;
    std::string lineTab;
    INfile.open( filename.c_str() );
    if(INfile.fail()){
        std::cout<<"error: in THDMWcache, table doesn't exist!"<< filename <<std::endl;
    }

    gslpp::matrix<double> arrayTab(rowN, colN, 0.);
    int a =0;
    int b=0;
    double v;

    while(INfile.good()){
        while(getline(INfile, lineTab)){
            if( lineTab[0]=='#' )continue;
            else{
            std::istringstream streamTab(lineTab);
            b=0;
            while(streamTab >>v){
                arrayTab.assign(a,b,v);
                b++;
            }
            a++;
            }
        }
    }

    INfile.close();
    
    return arrayTab;
}

//1D interpolation

double THDMWcache::interpolate(gslpp::matrix<double> arrayTab, double x){

    int rowN=arrayTab.size_i();
    
    double xmin = arrayTab(0,0);
    double xmax = arrayTab(rowN-1,0);
    double interval = arrayTab(1,0)-arrayTab(0,0);
    int Nintervals = (x-xmin)/interval;
    double y = 0.0;
       
    if(x<xmin){
//        std::cout<<"warning: your table parameter value is smaller than the minimum allowed value"<<std::endl;
        return 0.;
    }
    else if(x>xmax){
//        std::cout<<"warning: your table parameter value is greater than the maximum allowed value"<<std::endl;
        return 0.;
    }
    else{
        y =(arrayTab(Nintervals+1,1)-arrayTab(Nintervals,1))/(arrayTab(Nintervals+1,0)
                   -arrayTab(Nintervals,0))*(x-arrayTab(Nintervals,0))+arrayTab(Nintervals,1);
        return y;
    }
}


//3D interpolation

double THDMWcache::interpolate3D(gslpp::matrix<double> arrayTab, double x, double y, double z){

    int rowN=arrayTab.size_i();
    double xmin = arrayTab(0,0);
    double xmax = arrayTab(rowN-1,0);
    double ymin = arrayTab(0,1);
    double ymax = arrayTab(rowN-1,1);
    double zmin = arrayTab(0,2);
    double zmax = arrayTab(rowN-1,2);
    double intervalx = arrayTab(1,0)-arrayTab(0,0);
    int iy=1;
    do iy++;
    while(arrayTab(iy,1)-arrayTab(iy-1,1)==0&&iy<6000000);
    double intervaly = arrayTab(iy,1)-arrayTab(iy-1,1);
    int iz=1;
    do iz++;
    while(arrayTab(iz,2)-arrayTab(iz-1,2)==0&&iz<6000000);
    double intervalz = arrayTab(iz,2)-arrayTab(iz-1,2);
    int Nintervalsx = (x-xmin)/intervalx;
    int Nintervalsy = (y-ymin)/intervaly;
    int Nintervalsz = (z-zmin)/intervalz;  
    int MaxNintervalx = round((xmax-xmin)/intervalx);
    int MaxNintervaly = round((ymax-ymin)/intervaly);
    int MaxNintervalz = round((zmax-zmin)/intervalz);
    //std::cout<<"MaxNintervalx="<<MaxNintervalx<<std::endl;
    //std::cout<<"MaxNintervaly="<<MaxNintervaly<<std::endl;
    //std::cout<<"MaxNintervalz="<<MaxNintervalz<<std::endl;
    //std::cout<<"imax="<<iz*Nintervalsz+iy*Nintervalsy+Nintervalsx<<std::endl;
    //std::cout<<"imax+1="<<(iz)*(Nintervalsz+1)+(iy)*(Nintervalsy+1)+Nintervalsx+1<<std::endl;
    //std::cout<<"Nintervalx="<<Nintervalsx<<std::endl;
    //std::cout<<"Nintervaly="<<Nintervalsy<<std::endl;
    //std::cout<<"Nintervalz="<<Nintervalsz<<std::endl;
    if(x<xmin||Nintervalsx>MaxNintervalx||y<ymin||Nintervalsy>MaxNintervaly||z<zmin||Nintervalsz>MaxNintervalz){
        //std::cout<<"warning: the parameter point lies outside the table"<<std::endl;
        //std::cout<<"x="<<x<<std::endl;
        //std::cout<<"y="<<y<<std::endl;
        //std::cout<<"z="<<z<<std::endl;
        return 0.;
    }
    else{
    
    double x1=arrayTab(iz*Nintervalsz+iy*Nintervalsy+Nintervalsx,0);
    double x2=arrayTab(iz*(Nintervalsz)+iy*(Nintervalsy)+Nintervalsx+1,0);
    double y1=arrayTab(iz*Nintervalsz+iy*Nintervalsy+Nintervalsx,1);
    double y2=arrayTab(iz*(Nintervalsz)+iy*(Nintervalsy+1)+Nintervalsx,1);
    double z1=arrayTab(iz*Nintervalsz+iy*Nintervalsy+Nintervalsx,2);
    double z2=arrayTab(iz*(Nintervalsz+1)+iy*(Nintervalsy)+Nintervalsx,2);

    return (arrayTab(iz*Nintervalsz+iy*Nintervalsy+Nintervalsx,3) * (x2-x) * (y2-y) * (z2-z)
            +arrayTab(iz*Nintervalsz+iy*Nintervalsy+Nintervalsx+1,3) * (x-x1) * (y2-y) * (z2-z)
            +arrayTab(iz*Nintervalsz+iy*(Nintervalsy+1)+Nintervalsx,3) * (x2-x) * (y-y1) * (z2-z) 
            +arrayTab(iz*(Nintervalsz+1)+iy*Nintervalsy+Nintervalsx,3) * (x2-x) * (y2-y) * (z-z1) 
            +arrayTab(iz*Nintervalsz+iy*(Nintervalsy+1)+Nintervalsx+1,3) * (x-x1) * (y-y1) * (z2-z)
            +arrayTab(iz*(Nintervalsz+1)+iy*Nintervalsy+Nintervalsx+1,3) * (x-x1) * (y2-y) * (z-z1)
            +arrayTab(iz*(Nintervalsz+1)+iy*(Nintervalsy+1)+Nintervalsx,3) * (x2-x) * (y-y1) * (z-z1)
            +arrayTab(iz*(Nintervalsz+1)+iy*(Nintervalsy+1)+Nintervalsx+1,3) * (x-x1) * (y-y1) * (z-z1))/((x2-x1)*(y2-y1)*(z2-z1));
    
    }
}



//4D interpolation

double THDMWcache::interpolate4D(gslpp::matrix<double> arrayTab, double x, double y, double z, double v){

    int rowN=arrayTab.size_i();

    double xmin = arrayTab(0,0);
    double xmax = arrayTab(rowN-1,0);
    double ymin = arrayTab(0,1);
    double ymax = arrayTab(rowN-1,1);
    double zmin = arrayTab(0,2);
    double zmax = arrayTab(rowN-1,2);
    double vmin = arrayTab(0,3);
    double vmax = arrayTab(rowN-1,3);
    double intervalx = arrayTab(1,0)-arrayTab(0,0);
    int iy=1;
    do iy++;
    while(arrayTab(iy,1)-arrayTab(iy-1,1)==0&&iy<6000000);
    double intervaly = arrayTab(iy,1)-arrayTab(iy-1,1);
    int iz=1;
    do iz++;
    while(arrayTab(iz,2)-arrayTab(iz-1,2)==0&&iz<6000000);
    double intervalz = arrayTab(iz,2)-arrayTab(iz-1,2);
    int iv=1;
    do iv++;
    while(arrayTab(iv,3)-arrayTab(iv-1,3)==0&&iv<6000000);
    double intervalv = arrayTab(iv,3)-arrayTab(iv-1,3);
    int Nintervalsx = (x-xmin)/intervalx;
    int Nintervalsy = (y-ymin)/intervaly;
    int Nintervalsz = (z-zmin)/intervalz;
    int Nintervalsv = (v-vmin)/intervalv;  
    int MaxNintervalx = round((xmax-xmin)/intervalx);
    int MaxNintervaly = round((ymax-ymin)/intervaly);
    int MaxNintervalz = round((zmax-zmin)/intervalz);
    int MaxNintervalv = round((vmax-vmin)/intervalv);
    //std::cout<<"xmin="<<xmin<<std::endl;
    //std::cout<<"xmax="<<xmax<<std::endl;
    //std::cout<<"ymin="<<ymin<<std::endl;
    //std::cout<<"ymax="<<ymax<<std::endl;
    //std::cout<<"zmin="<<zmin<<std::endl;
    //std::cout<<"zmax="<<zmax<<std::endl;
    //std::cout<<"vmin="<<vmin<<std::endl;
    //std::cout<<"vmax="<<vmax<<std::endl;
    //std::cout<<"intervalx="<<intervalx<<std::endl;
    //std::cout<<"intervaly="<<intervaly<<std::endl;
    //std::cout<<"intervalz="<<intervalz<<std::endl;
    //std::cout<<"intervalv="<<intervalv<<std::endl;
    //std::cout<<"Nintervalsx="<<Nintervalsx<<std::endl;
    //std::cout<<"Nintervalsy="<<Nintervalsy<<std::endl;
    //std::cout<<"Nintervalsz="<<Nintervalsz<<std::endl;
    //std::cout<<"Nintervalsv="<<Nintervalsv<<std::endl;
    
    if(x<xmin||Nintervalsx>MaxNintervalx||y<ymin||Nintervalsy>MaxNintervaly||z<zmin||Nintervalsz>MaxNintervalz||v<vmin||Nintervalsv>MaxNintervalv){
        //std::cout<<"warning: the parameter point lies outside the table"<<std::endl;
        //std::cout<<"x="<<x<<std::endl;
        //std::cout<<"y="<<y<<std::endl;
        //std::cout<<"z="<<z<<std::endl;
        //std::cout<<"v="<<v<<std::endl;
        return 0.;
    }
    else{
    double x1=arrayTab(iv*Nintervalsv+iz*Nintervalsz+iy*Nintervalsy+Nintervalsx,0);
    double x2=arrayTab(iv*(Nintervalsv)+iz*(Nintervalsz)+iy*(Nintervalsy)+Nintervalsx+1,0);
    double y1=arrayTab(iv*Nintervalsv+iz*Nintervalsz+iy*Nintervalsy+Nintervalsx,1);
    double y2=arrayTab(iv*(Nintervalsv)+iz*(Nintervalsz)+iy*(Nintervalsy+1)+Nintervalsx,1);
    double z1=arrayTab(iv*Nintervalsv+iz*Nintervalsz+iy*Nintervalsy+Nintervalsx,2);
    double z2=arrayTab(iv*(Nintervalsv)+iz*(Nintervalsz+1)+iy*(Nintervalsy)+Nintervalsx,2);
    double v1=arrayTab(iv*Nintervalsv+iz*Nintervalsz+iy*Nintervalsy+Nintervalsx,3);
    double v2=arrayTab(iv*(Nintervalsv+1)+iz*(Nintervalsz)+iy*(Nintervalsy)+Nintervalsx,3);
    //std::cout<<"Nmax="<<iv*(Nintervalsv+1)+iz*(Nintervalsz+1)+iy*(Nintervalsy+1)+Nintervalsx+1<<std::endl;
    //std::cout<<"x1="<<x1<<std::endl;
    //std::cout<<"x2="<<x2<<std::endl;
    //std::cout<<"y1="<<y1<<std::endl;
    //std::cout<<"y2="<<y2<<std::endl;
    //std::cout<<"z1="<<z1<<std::endl;
    //std::cout<<"z2="<<z2<<std::endl;
    //std::cout<<"v1="<<v1<<std::endl;
    //std::cout<<"v2="<<v2<<std::endl;
    /*std::cout<<"Interpolation= "<<(arrayTab(iv*Nintervalsv+iz*Nintervalsz+iy*Nintervalsy+Nintervalsx,4) * (x2-x) * (y2-y) * (z2-z) * (v2-v)
            +arrayTab(iv*Nintervalsv+iz*Nintervalsz+iy*Nintervalsy+Nintervalsx+1,4) * (x-x1) * (y2-y) * (z2-z) * (v2-v)
            +arrayTab(iv*Nintervalsv+iz*Nintervalsz+iy*(Nintervalsy+1)+Nintervalsx,4) * (x2-x) * (y-y1) * (z2-z) * (v2-v)
            +arrayTab(iv*Nintervalsv+iz*(Nintervalsz+1)+iy*Nintervalsy+Nintervalsx,4) * (x2-x) * (y2-y) * (z-z1) * (v2-v)
            +arrayTab(iv*(Nintervalsv+1)+iz*Nintervalsz+iy*Nintervalsy+Nintervalsx,4) * (x2-x) * (y2-y) * (z2-z) * (v-v1)
            +arrayTab(iv*Nintervalsv+iz*Nintervalsz+iy*(Nintervalsy+1)+Nintervalsx+1,4) * (x-x1) * (y-y1) * (z2-z) * (v2-v)
            +arrayTab(iv*Nintervalsv+iz*(Nintervalsz+1)+iy*Nintervalsy+Nintervalsx+1,4) * (x-x1) * (y2-y) * (z-z1) * (v2-v)
            +arrayTab(iv*(Nintervalsv+1)+iz*Nintervalsz+iy*Nintervalsy+Nintervalsx+1,4) * (x-x1) * (y2-y) * (z2-z) * (v-v1)
            +arrayTab(iv*Nintervalsv+iz*(Nintervalsz+1)+iy*(Nintervalsy+1)+Nintervalsx,4) * (x2-x) * (y-y1) * (z-z1) * (v2-v)
            +arrayTab(iv*(Nintervalsv+1)+iz*Nintervalsz+iy*(Nintervalsy+1)+Nintervalsx,4) * (x2-x) * (y-y1) * (z2-z) * (v-v1)
            +arrayTab(iv*(Nintervalsv+1)+iz*(Nintervalsz+1)+iy*Nintervalsy+Nintervalsx,4) * (x2-x) * (y2-y) * (z-z1) * (v-v1)
            +arrayTab(iv*Nintervalsv+iz*(Nintervalsz+1)+iy*(Nintervalsy+1)+Nintervalsx+1,4) * (x-x1) * (y-y1) * (z-z1) * (v2-v)
            +arrayTab(iv*(Nintervalsv+1)+iz*Nintervalsz+iy*(Nintervalsy+1)+Nintervalsx+1,4) * (x-x1) * (y-y1) * (z2-z) * (v-v1)
            +arrayTab(iv*(Nintervalsv+1)+iz*(Nintervalsz+1)+iy*Nintervalsy+Nintervalsx+1,4) * (x-x1) * (y2-y) * (z-z1) * (v-v1)
            +arrayTab(iv*(Nintervalsv+1)+iz*(Nintervalsz+1)+iy*(Nintervalsy+1)+Nintervalsx,4) * (x2-x) * (y-y1) * (z-z1) * (v-v1)
            +arrayTab(iv*(Nintervalsv+1)+iz*(Nintervalsz+1)+iy*(Nintervalsy+1)+Nintervalsx+1,4) * (x-x1) * (y-y1) * (z-z1) * (v-v1))/((x2-x1)*(y2-y1)*(z2-z1)*(v2-v1))<<std::endl;*/
    return (arrayTab(iv*Nintervalsv+iz*Nintervalsz+iy*Nintervalsy+Nintervalsx,4) * (x2-x) * (y2-y) * (z2-z) * (v2-v)
            +arrayTab(iv*Nintervalsv+iz*Nintervalsz+iy*Nintervalsy+Nintervalsx+1,4) * (x-x1) * (y2-y) * (z2-z) * (v2-v)
            +arrayTab(iv*Nintervalsv+iz*Nintervalsz+iy*(Nintervalsy+1)+Nintervalsx,4) * (x2-x) * (y-y1) * (z2-z) * (v2-v)
            +arrayTab(iv*Nintervalsv+iz*(Nintervalsz+1)+iy*Nintervalsy+Nintervalsx,4) * (x2-x) * (y2-y) * (z-z1) * (v2-v)
            +arrayTab(iv*(Nintervalsv+1)+iz*Nintervalsz+iy*Nintervalsy+Nintervalsx,4) * (x2-x) * (y2-y) * (z2-z) * (v-v1)
            +arrayTab(iv*Nintervalsv+iz*Nintervalsz+iy*(Nintervalsy+1)+Nintervalsx+1,4) * (x-x1) * (y-y1) * (z2-z) * (v2-v)
            +arrayTab(iv*Nintervalsv+iz*(Nintervalsz+1)+iy*Nintervalsy+Nintervalsx+1,4) * (x-x1) * (y2-y) * (z-z1) * (v2-v)
            +arrayTab(iv*(Nintervalsv+1)+iz*Nintervalsz+iy*Nintervalsy+Nintervalsx+1,4) * (x-x1) * (y2-y) * (z2-z) * (v-v1)
            +arrayTab(iv*Nintervalsv+iz*(Nintervalsz+1)+iy*(Nintervalsy+1)+Nintervalsx,4) * (x2-x) * (y-y1) * (z-z1) * (v2-v)
            +arrayTab(iv*(Nintervalsv+1)+iz*Nintervalsz+iy*(Nintervalsy+1)+Nintervalsx,4) * (x2-x) * (y-y1) * (z2-z) * (v-v1)
            +arrayTab(iv*(Nintervalsv+1)+iz*(Nintervalsz+1)+iy*Nintervalsy+Nintervalsx,4) * (x2-x) * (y2-y) * (z-z1) * (v-v1)
            +arrayTab(iv*Nintervalsv+iz*(Nintervalsz+1)+iy*(Nintervalsy+1)+Nintervalsx+1,4) * (x-x1) * (y-y1) * (z-z1) * (v2-v)
            +arrayTab(iv*(Nintervalsv+1)+iz*Nintervalsz+iy*(Nintervalsy+1)+Nintervalsx+1,4) * (x-x1) * (y-y1) * (z2-z) * (v-v1)
            +arrayTab(iv*(Nintervalsv+1)+iz*(Nintervalsz+1)+iy*Nintervalsy+Nintervalsx+1,4) * (x-x1) * (y2-y) * (z-z1) * (v-v1)
            +arrayTab(iv*(Nintervalsv+1)+iz*(Nintervalsz+1)+iy*(Nintervalsy+1)+Nintervalsx,4) * (x2-x) * (y-y1) * (z-z1) * (v-v1)
            +arrayTab(iv*(Nintervalsv+1)+iz*(Nintervalsz+1)+iy*(Nintervalsy+1)+Nintervalsx+1,4) * (x-x1) * (y-y1) * (z-z1) * (v-v1))/((x2-x1)*(y2-y1)*(z2-z1)*(v2-v1));
    }
}

/*
//4D log interpolation

double THDMWcache::loginterpolate4D(gslpp::matrix<double> arrayTab, double x, double y, double z, double v){

    int rowN=arrayTab.size_i();

    double xmin = arrayTab(0,0);
    double xmax = arrayTab(rowN-1,0);
    double ymin = arrayTab(0,1);
    double ymax = arrayTab(rowN-1,1);
    double zmin = arrayTab(0,2);
    double zmax = arrayTab(rowN-1,2);
    double vmin = arrayTab(0,3);
    double vmax = arrayTab(rowN-1,3);
    double intervalx = arrayTab(1,0)-arrayTab(0,0);
    int iy=1;
    do iy++;
    while(arrayTab(iy,1)-arrayTab(iy-1,1)==0&&iy<6000000);
    double intervaly = arrayTab(iy,1)-arrayTab(iy-1,1);
    int iz=1;
    do iz++;
    while(arrayTab(iz,2)-arrayTab(iz-1,2)==0&&iz<6000000);
    double intervalz = arrayTab(iz,2)-arrayTab(iz-1,2);
    int iv=1;
    do iv++;
    while(arrayTab(iv,3)-arrayTab(iv-1,3)==0&&iv<6000000);
    double intervalv = arrayTab(iv,3)-arrayTab(iv-1,3);
    int Nintervalsx = (x-xmin)/intervalx;
    int Nintervalsy = (y-ymin)/intervaly;
    int Nintervalsz = (z-zmin)/intervalz;
    int Nintervalsv = (v-vmin)/intervalv;       
    
    if(x<xmin||x>xmax||y<ymin||y>ymax||z<zmin||z>zmax||v<vmin||v>vmax){
        std::cout<<"warning: the parameter point lies outside the table"<<std::endl;
        return 0.;
    }
    else{
    double x1=arrayTab(iv*Nintervalsv+iz*Nintervalsz+iy*Nintervalsy+Nintervalsx,0);
    double x2=arrayTab(iv*(Nintervalsv)+iz*(Nintervalsz)+iy*(Nintervalsy)+Nintervalsx+1,0);
    double y1=arrayTab(iv*Nintervalsv+iz*Nintervalsz+iy*Nintervalsy+Nintervalsx,1);
    double y2=arrayTab(iv*(Nintervalsv)+iz*(Nintervalsz)+iy*(Nintervalsy+1)+Nintervalsx,1);
    double z1=arrayTab(iv*Nintervalsv+iz*Nintervalsz+iy*Nintervalsy+Nintervalsx,2);
    double z2=arrayTab(iv*(Nintervalsv)+iz*(Nintervalsz+1)+iy*(Nintervalsy)+Nintervalsx,2);
    double v1=arrayTab(iv*Nintervalsv+iz*Nintervalsz+iy*Nintervalsy+Nintervalsx,3);
    double v2=arrayTab(iv*(Nintervalsv+1)+iz*(Nintervalsz)+iy*(Nintervalsy)+Nintervalsx,3);
    double N1=1e-15;
    double N2=1e-15;
    double N3=1e-15;
    double N4=1e-15;
    double N5=1e-15;
    double N6=1e-15;
    double N7=1e-15;
    double N8=1e-15;
    double N9=1e-15;
    double N10=1e-15;
    double N11=1e-15;
    double N12=1e-15;
    double N13=1e-15;
    double N14=1e-15;
    double N15=1e-15;
    double N16=1e-15;
    
    if(sqrt(pow(arrayTab(iv*Nintervalsv+iz*Nintervalsz+iy*Nintervalsy+Nintervalsx,4),2))>1e-15) N1=arrayTab(iv*Nintervalsv+iz*Nintervalsz+iy*Nintervalsy+Nintervalsx,4);
    if(sqrt(pow(arrayTab(iv*Nintervalsv+iz*Nintervalsz+iy*Nintervalsy+Nintervalsx+1,4),2)>1e-15) N2=arrayTab(iv*Nintervalsv+iz*Nintervalsz+iy*Nintervalsy+Nintervalsx+1,4);
    if(sqrt(pow(arrayTab(iv*Nintervalsv+iz*Nintervalsz+iy*(Nintervalsy+1)+Nintervalsx,4),2))>1e-15) N3=arrayTab(iv*Nintervalsv+iz*Nintervalsz+iy*(Nintervalsy+1)+Nintervalsx,4);
    if(sqrt(pow(arrayTab(iv*Nintervalsv+iz*(Nintervalsz+1)+iy*Nintervalsy+Nintervalsx,4),2))>1e-15) N4=arrayTab(iv*Nintervalsv+iz*(Nintervalsz+1)+iy*Nintervalsy+Nintervalsx,4);
    if(sqrt(pow(arrayTab(iv*(Nintervalsv+1)+iz*Nintervalsz+iy*Nintervalsy+Nintervalsx,4),2))>1e-15) N5=arrayTab(iv*(Nintervalsv+1)+iz*Nintervalsz+iy*Nintervalsy+Nintervalsx,4);
    if(sqrt(pow(arrayTab(iv*Nintervalsv+iz*Nintervalsz+iy*(Nintervalsy+1)+Nintervalsx+1,4),2))>1e-15) N6=arrayTab(iv*Nintervalsv+iz*Nintervalsz+iy*(Nintervalsy+1)+Nintervalsx+1,4);
    if(sqrt(pow(arrayTab(iv*Nintervalsv+iz*(Nintervalsz+1)+iy*Nintervalsy+Nintervalsx+1,4),2))>1e-15) N7=arrayTab(iv*Nintervalsv+iz*(Nintervalsz+1)+iy*Nintervalsy+Nintervalsx+1,4);
    if(sqrt(pow(arrayTab(iv*(Nintervalsv+1)+iz*Nintervalsz+iy*Nintervalsy+Nintervalsx+1,4),2))>1e-15) N8=arrayTab(iv*(Nintervalsv+1)+iz*Nintervalsz+iy*Nintervalsy+Nintervalsx+1,4);
    if(sqrt(pow(arrayTab(iv*Nintervalsv+1+iz*(Nintervalsz+1)+iy*(Nintervalsy+1)+Nintervalsx,4),2))>1e-15) N9=arrayTab(iv*Nintervalsv+1+iz*(Nintervalsz+1)+iy*(Nintervalsy+1)+Nintervalsx,4);
    if(sqrt(pow(arrayTab(iv*(Nintervalsv+1)+iz*Nintervalsz+iy*(Nintervalsy+1)+Nintervalsx,4),2))>1e-15) N10=arrayTab(iv*(Nintervalsv+1)+iz*Nintervalsz+iy*(Nintervalsy+1)+Nintervalsx,4);
    if(sqrt(pow(arrayTab(iv*(Nintervalsv+1)+iz*(Nintervalsz+1)+iy*Nintervalsy+Nintervalsx,4),2))>1e-15) N11=arrayTab(iv*(Nintervalsv+1)+iz*(Nintervalsz+1)+iy*Nintervalsy+Nintervalsx,4);
    if(sqrt(pow(arrayTab(iv*Nintervalsv+iz*(Nintervalsz+1)+iy*(Nintervalsy+1)+Nintervalsx+1,4),2))>1e-15) N12=arrayTab(iv*Nintervalsv+iz*(Nintervalsz+1)+iy*(Nintervalsy+1)+Nintervalsx+1,4);
    if(sqrt(pow(arrayTab(iv*(Nintervalsv+1)+iz*Nintervalsz+iy*(Nintervalsy+1)+Nintervalsx+1,4),2))>1e-15) N13=arrayTab(iv*(Nintervalsv+1)+iz*Nintervalsz+iy*(Nintervalsy+1)+Nintervalsx+1,4);
    if(sqrt(pow(arrayTab(iv*(Nintervalsv+1)+iz*(Nintervalsz+1)+iy*Nintervalsy+Nintervalsx+1,4),2))>1e-15) N14=arrayTab(iv*(Nintervalsv+1)+iz*(Nintervalsz+1)+iy*Nintervalsy+Nintervalsx+1,4);
    if(sqrt(pow(arrayTab(iv*(Nintervalsv+1)+iz*(Nintervalsz+1)+iy*(Nintervalsy+1)+Nintervalsx,4),2))>1e-15) N15=arrayTab(iv*(Nintervalsv+1)+iz*(Nintervalsz+1)+iy*(Nintervalsy+1)+Nintervalsx,4);
    if(sqrt(pow(arrayTab(iv*(Nintervalsv+1)+iz*(Nintervalsz+1)+iy*(Nintervalsy+1)+Nintervalsx+1,4),2))>1e-15) N16=arrayTab(iv*(Nintervalsv+1)+iz*(Nintervalsz+1)+iy*(Nintervalsy+1)+Nintervalsx+1,4);
    
    
    
    
    return (log10(N1) * (x2-x) * (y2-y) * (z2-z) * (v2-v)
            +log10(N2) * (x-x1) * (y2-y) * (z2-z) * (v2-v)
            +log10(N3) * (x2-x) * (y-y1) * (z2-z) * (v2-v)
            +log10(N4) * (x2-x) * (y2-y) * (z-z1) * (v2-v)
            +log10(N5) * (x2-x) * (y2-y) * (z2-z) * (v-v1)
            +log10(N6) * (x-x1) * (y-y1) * (z2-z) * (v2-v)
            +log10(N7) * (x-x1) * (y2-y) * (z-z1) * (v2-v)
            +log10(N8) * (x-x1) * (y2-y) * (z2-z) * (v-v1)
            +log10(N9) * (x2-x) * (y-y1) * (z-z1) * (v2-v)
            +log10(N10) * (x2-x) * (y-y1) * (z2-z) * (v-v1)
            +log10(N11) * (x2-x) * (y2-y) * (z-z1) * (v-v1)
            +log10(N12) * (x-x1) * (y-y1) * (z-z1) * (v2-v)
            +log10(N13) * (x-x1) * (y-y1) * (z2-z) * (v-v1)
            +log10(N14) * (x-x1) * (y2-y) * (z-z1) * (v-v1)
            +log10(N15) * (x2-x) * (y-y1) * (z-z1) * (v-v1)
            +log10(N16) * (x-x1) * (y-y1) * (z-z1) * (v-v1))/((x2-x1)*(y2-y1)*(z2-z1)*(v2-v1));
    }
}
*/

/*
double THDMWcache::ip_cs_ppto2Sto4t_13(double etaD, double etaU, double THDMW_nu4, double mSR){
    int NumPar = 4;
    double params[] = {etaD, etaU, THDMW_nu4, mSR};

    int i = CacheCheck(ip_cs_ppto2Sto4t_13_cache, NumPar, params);
    if (i>=0) {
        return ( ip_cs_ppto2Sto4t_13_cache[NumPar][i] );
    } else {
        double newResult = 0.0;
        if (mSR>=500. && mSR <=1500.) {
            newResult = interpolate4D(log_cs_ggHp_8, etaD, etaU, THDMW_nu4, mSR);
        }
        CacheShift(ip_cs_ppto2Sto4t_13_cache, NumPar, params, newResult);
        return newResult;
    }
}*/



//double THDMWcache::ip_cs_ppto2Sto4t_13(double THDMW_nu1, double THDMW_nu2, double THDMW_nu4, double THDMW_mS2){
//    int NumPar = 4;
//    double params[] = {THDMW_nu1, THDMW_nu2, THDMW_nu4, THDMW_mS2};

////    int i = CacheCheck(ip_cs_ppto2Sto4t_13_cache, NumPar, params);
////    if (i>=0) {
////        return ( ip_cs_ppto2Sto4t_13_cache[NumPar][i] );
////    } else {
//        double newResult = 0.0;
//        if (THDMW_mS2>=250000. && THDMW_mS2 <=1000000.) {
//            newResult = interpolate4D(log_cs_ggH_13, THDMW_nu1, THDMW_nu2, THDMW_nu4, sqrt(THDMW_mS2));
//        }
////        CacheShift(ip_cs_ppto2Sto4t_13_cache, NumPar, params, newResult);
//        return newResult;
//   }

void THDMWcache::read(){



    std::stringstream ex0,ex1,ex2,ex3;
    std::stringstream ex1e,ex2e,ex3e;
//    std::stringstream ex14ep2,ex14em2;
    std::stringstream ex4,ex5,ex6,ex7,ex8;
    std::stringstream ex4e,ex5e,ex6e,ex7e,ex8e;
    std::stringstream ex9,ex10,ex13;
    std::stringstream ex9e,ex10e,ex13e;
    std::stringstream ex14, ex15,ex16,ex17,ex18,ex19;
    std::stringstream th1,th2,th3,th4,th5,th6,th7,th8,th9,th10,th11,th12,th13,th14;

    std::stringstream bsg1;

    std::cout<<"reading tables"<<std::endl;

//    std::cout << "HEPFITTABS = " << getenv("HEPFITPATH") << std::endl;
    std::stringstream path;
    path << getenv("HEPFITTABS") << "/THDM/tabs/";
//    path << "/Users/victormirallesaznar/tabs/";
//    std::cout << path.str() << std::endl;
    std::string tablepath=path.str();
//    std::cout << tablepath << std::endl;



//// THIS IS FOR THE FUTURE IMPLEMENTATION INTO HEADERS:
//    std::cout<<"br_tt="<<br_tt<<std::endl;
//    double brtt1[4][2];
//    brtt1[0][1]=1;
//        gslpp::matrix<double> brtt1(19861,2,0.);
//    std::stringstream br1x;
//    br1x << "log_cs_ggH_13.h";
//      //brtt1(2)=(3.,4.);
//      brtt1=readTable(br1x.str(),20,2);
//    std::cout<<"brtt1="<<bla1<<std::endl;



    ex0 << tablepath << "150304114.dat";//Dummy will be deleted by Scientific Linux  
    Dummy = readTable(ex0.str(),167,2); 
    ex1 << tablepath << "150304114.dat";              
    CMS8_pp_H_hh_bbbb = readTable(ex1.str(),167,2);   
    ex1e << tablepath << "150304114_e.dat";          
    CMS8_pp_H_hh_bbbb_e = readTable(ex1e.str(),167,2); 
    ex2 << tablepath << "150608329.dat";        
    CMS8_bb_phi_bb = readTable(ex2.str(),81,2); 
    ex2e << tablepath << "150608329_e.dat";      
    CMS8_bb_phi_bb_e = readTable(ex2e.str(),81,2);
    ex3 << tablepath << "150507018.dat";         
    ATLAS8_gg_phi_tt = readTable(ex3.str(),53,2); 
    ex3e << tablepath << "150507018_e.dat";       
    ATLAS8_gg_phi_tt_e = readTable(ex3e.str(),53,2);

//    ex14ep1 << tablepath << "150602301_ep1.dat";
//    CMS_ggF_phi_gaga_ep1 = readTable(ex14ep1.str(),141,2);
    //CHANGE THIS DEFINITION!
//    ex14ep2 << tablepath << "150602301_e.dat";
//    CMS_ggF_phi_gaga_ep2 = readTable(ex14ep2.str(),141,2);
//    ex14em1 << tablepath << "150602301_em1.dat";
//    CMS_ggF_phi_gaga_em1 = readTable(ex14em1.str(),141,2);
    //CHANGE THIS DEFINITION!
//    ex14em2 << tablepath << "150602301_e.dat";
//    CMS_ggF_phi_gaga_em2 = readTable(ex14em2.str(),141,2);


    
    ex4 << tablepath << "ATLAS-CONF-2016-104_b.dat";    
    ATLAS13_bb_phi_tt = readTable(ex4.str(),61,2);      
    ex4e << tablepath << "ATLAS-CONF-2016-104_b_e.dat"; 
    ATLAS13_bb_phi_tt_e = readTable(ex4e.str(),61,2);   
    ex5 << tablepath << "180711883.dat";    
    ATLAS13_tt_phi_tt = readTable(ex5.str(),61,2);      
    ex5e << tablepath << "ATLAS-CONF-2016-104_a_e.dat";
    ATLAS13_tt_phi_tt_e = readTable(ex5e.str(),61,2);  
    ex6 << tablepath << "ATLAS-CONF-2016-049.dat";     
    ATLAS13_pp_H_hh_bbbb = readTable(ex6.str(),271,2); 
    ex6e << tablepath << "ATLAS-CONF-2016-049_e.dat";  
    ATLAS13_pp_H_hh_bbbb_e = readTable(ex6e.str(),271,2); 
    ex7 << tablepath << "CMS-PAS-HIG-16-025.dat";
    CMS13_pp_phi_bb = readTable(ex7.str(),66,2); 
    ex7e << tablepath << "CMS-PAS-HIG-16-025_e.dat"; 
    CMS13_pp_phi_bb_e = readTable(ex7e.str(),66,2);  
    ex8 << tablepath << "180603548.dat";
    CMS13_pp_H_hh_bbbb = readTable(ex8.str(),95,2);
    ex8e << tablepath << "180603548_e.dat";
    CMS13_pp_H_hh_bbbb_e = readTable(ex8e.str(),95,2);



    ex9 << tablepath << "151203704.dat";
    ATLAS8_pp_Hpm_tb = readTable(ex9.str(),41,2); 
    ex9e << tablepath << "151203704_e.dat";       
    ATLAS8_pp_Hpm_tb_e = readTable(ex9e.str(),41,2);
    ex10 << tablepath << "150807774_b.dat";
    CMS8_pp_Hp_tb = readTable(ex10.str(),43,2);  
    ex10e << tablepath << "150807774_b_e.dat";   
    CMS8_pp_Hp_tb_e = readTable(ex10e.str(),43,2);
    ex17 << tablepath << "210210076.dat";
    ATLAS13_pp_Hp_tb = readTable(ex17.str(),181,2); 
    ex18 << tablepath << "180512191.dat";
    CMS13_bb_H_bb = readTable(ex18.str(),101,2); 
//    ex11 << tablepath << "ATLAS-CONF-2016-089.dat";
//    ATLAS13_pp_Hp_tb1 = readTable(ex11.str(),71,2);        
//    ex11e << tablepath << "ATLAS-CONF-2016-089_e.dat";     
//    ATLAS13_pp_Hp_tb1_e = readTable(ex11e.str(),71,2);     
//    ex12 << tablepath << "ATLAS-CONF-2016-104_c.dat"; 
//    ATLAS13_pp_Hp_tb2 = readTable(ex12.str(),181,2);  
//    ex12e << tablepath << "ATLAS-CONF-2016-104_c_e.dat"; 
//    ATLAS13_pp_Hp_tb2_e = readTable(ex12e.str(),181,2);  
    ex13 << tablepath << "171004960.dat";
    CMS13_ggF_H_hh_bbbb = readTable(ex13.str(),226,2);
    ex13e << tablepath << "171004960_e.dat";
    CMS13_ggF_H_hh_bbbb_e = readTable(ex13e.str(),226,2);
    ex14 << tablepath << "180410823_b.dat";
    ATLAS13_pp_Gkk_tt = readTable(ex14.str(),131,2);
    ex15 << tablepath << "CMS-CR-2018-204.dat";
    CMS13_pp_R_gg = readTable(ex15.str(),241,2);
    ex16 << tablepath << "171007171.dat";
    ATLAS13_pp_SS_jjjj = readTable(ex16.str(),126,2);
    
    ex19 << tablepath << "180206149.dat";
    CMS8_pp_phi_bb = readTable(ex19.str(),88,2); 

    th1 << tablepath << "Generated_data_S2t_Fixed_Steps.dat";
    MadGraph_pp_Sr_tt = readTable(th1.str(),22800,5);
    
    th2 << tablepath << "Generated_data_Stt_tttt_Fixed_Steps.dat";
    MadGraph_pp_Srtt_tttt = readTable(th2.str(),22800,5);
    
    th3 << tablepath << "Generated_data_S_jj_Fixed_Steps.dat";
    MadGraph_pp_Sr_jj = readTable(th3.str(),2940,5);
    
    th4 << tablepath << "Generated_data_SS_jjjj_Fixed_Steps.dat";
    MadGraph_pp_SrSr_jjjj = readTable(th4.str(),4200,5);

    th5 << tablepath << "Generated_data_Stb_tbtb_Fixed_Steps.dat";
    MadGraph_pp_Stb_tbtb = readTable(th5.str(),4332,4);
    
    th6 << tablepath << "Generated_data_Soddtt_tttt_Fixed_Steps.dat";
    MadGraph_pp_Sitt_tttt = readTable(th6.str(),9360,4);
    
    th7 << tablepath << "Generated_data_Srbb_bbbb_Fixed_Steps.dat";
    MadGraph_pp_Srbb_bbbb = readTable(th7.str(),15960,5);
    
    th8 << tablepath << "Generated_data_Sibb_bbbb_Fixed_Steps.dat";
    MadGraph_pp_Sibb_bbbb = readTable(th8.str(),8892,4);
    
    th9 << tablepath << "Generated_data_Sr_bb_Fixed_Steps.dat";
    MadGraph_pp_Sr_bb = readTable(th9.str(),15960,5);
    
    th10 << tablepath << "Generated_data_Sr_bb_8TeV_Fixed_Steps.dat";
    MadGraph_pp_Sr_bb_8TeV = readTable(th10.str(),15960,5);
    
    th11 << tablepath << "Generated_data_Si_bb_Fixed_Steps.dat";
    MadGraph_pp_Si_bb = readTable(th11.str(),8892,4);
    
    th12 << tablepath << "Generated_data_Si_bb_8TeV_Fixed_Steps.dat";
    MadGraph_pp_Si_bb_8TeV = readTable(th12.str(),8892,4);
    
    th13 << tablepath << "Generated_data_Srbb_bbbb_8TeV_Fixed_Steps.dat";
    MadGraph_pp_Srbb_bbbb_8TeV = readTable(th13.str(),15960,5);
    
    th14 << tablepath << "Generated_data_Sibb_bbbb_8TeV_Fixed_Steps.dat";
    MadGraph_pp_Sibb_bbbb_8TeV = readTable(th14.str(),8892,4);
    
    
    bsg1 << tablepath << "bsgammatable.dat";
    arraybsgamma = readTable(bsg1.str(),1111,3);
}    

double THDMWcache::ip_ex_pp_phi_hh_bbbb_CMS8(double mass){
    int NumPar = 1;
    double params[] = {mass};

    int i = CacheCheckReal(ip_ex_pp_phi_hh_bbbb_CMS8_cache, NumPar, params);
    if (i>=0) {
        return ( ip_ex_pp_phi_hh_bbbb_CMS8_cache[NumPar][i] );
    } else {
        double newResult = interpolate(CMS8_pp_H_hh_bbbb,mass);
        CacheShiftReal(ip_ex_pp_phi_hh_bbbb_CMS8_cache, NumPar, params, newResult);
        return newResult;
    }
}

double THDMWcache::ip_ex_pp_phi_hh_bbbb_CMS8_e(double mass){
    int NumPar = 1;
    double params[] = {mass};

    int i = CacheCheckReal(ip_ex_pp_phi_hh_bbbb_CMS8_cache_e, NumPar, params);
    if (i>=0) {
        return ( ip_ex_pp_phi_hh_bbbb_CMS8_cache_e[NumPar][i] );
    } else {
        double newResult = interpolate(CMS8_pp_H_hh_bbbb_e,mass);
        CacheShiftReal(ip_ex_pp_phi_hh_bbbb_CMS8_cache_e, NumPar, params, newResult);
        return newResult;
    }
}
        
double THDMWcache::ip_ex_bb_phi_bb_CMS8(double mass){
    int NumPar = 1;
    double params[] = {mass};

    int i = CacheCheckReal(ip_ex_bb_phi_bb_CMS8_cache, NumPar, params);
    if (i>=0) {
        return ( ip_ex_bb_phi_bb_CMS8_cache[NumPar][i] );
    } else {
        double newResult = interpolate (CMS8_bb_phi_bb,mass);
        CacheShiftReal(ip_ex_bb_phi_bb_CMS8_cache, NumPar, params, newResult);
        return newResult;
    }
}



double THDMWcache::ip_ex_bb_phi_bb_CMS8_e(double mass){
    int NumPar = 1;
    double params[] = {mass};

    int i = CacheCheckReal(ip_ex_bb_phi_bb_CMS8_cache_e, NumPar, params);
    if (i>=0) {
        return ( ip_ex_bb_phi_bb_CMS8_cache_e[NumPar][i] );
    } else {
        double newResult = interpolate (CMS8_bb_phi_bb_e,mass);
        CacheShiftReal(ip_ex_bb_phi_bb_CMS8_cache_e, NumPar, params, newResult);
        return newResult;
    }
}   

double THDMWcache::ip_ex_gg_phi_tt_ATLAS8(double mass){
    int NumPar = 1;
    double params[] = {mass};

    int i = CacheCheckReal(ip_ex_gg_phi_tt_ATLAS8_cache, NumPar, params);
    if (i>=0) {
        return ( ip_ex_gg_phi_tt_ATLAS8_cache[NumPar][i] );
    } else {
        double newResult = interpolate (ATLAS8_gg_phi_tt,mass);
        CacheShiftReal(ip_ex_gg_phi_tt_ATLAS8_cache, NumPar, params, newResult);
        return newResult;
    }
}



double THDMWcache::ip_ex_gg_phi_tt_ATLAS8_e(double mass){
    int NumPar = 1;
    double params[] = {mass};

    int i = CacheCheckReal(ip_ex_gg_phi_tt_ATLAS8_cache_e, NumPar, params);
    if (i>=0) {
        return ( ip_ex_gg_phi_tt_ATLAS8_cache_e[NumPar][i] );
    } else {
        double newResult = interpolate (ATLAS8_gg_phi_tt_e,mass);
        CacheShiftReal(ip_ex_gg_phi_tt_ATLAS8_cache_e, NumPar, params, newResult);
        return newResult;
    }
}

double THDMWcache::ip_ex_bb_phi_tt_ATLAS13(double mass){
    int NumPar = 1;
    double params[] = {mass};

    int i = CacheCheckReal(ip_ex_bb_phi_tt_ATLAS13_cache, NumPar, params);
    if (i>=0) {
        return ( ip_ex_bb_phi_tt_ATLAS13_cache[NumPar][i] );
    } else {
        double newResult = interpolate (ATLAS13_bb_phi_tt,mass);
        CacheShiftReal(ip_ex_bb_phi_tt_ATLAS13_cache, NumPar, params, newResult);
        return newResult;
    }
}



double THDMWcache::ip_ex_bb_phi_tt_ATLAS13_e(double mass){
    int NumPar = 1;
    double params[] = {mass};

    int i = CacheCheckReal(ip_ex_bb_phi_tt_ATLAS13_cache_e, NumPar, params);
    if (i>=0) {
        return ( ip_ex_bb_phi_tt_ATLAS13_cache_e[NumPar][i] );
    } else {
        double newResult = interpolate (ATLAS13_bb_phi_tt_e,mass);
        CacheShiftReal(ip_ex_bb_phi_tt_ATLAS13_cache_e, NumPar, params, newResult);
        return newResult;
    }
}



double THDMWcache::ip_ex_tt_phi_tt_ATLAS13(double mass){
    int NumPar = 1;
    double params[] = {mass};

    int i = CacheCheckReal(ip_ex_tt_phi_tt_ATLAS13_cache, NumPar, params);
    if (i>=0) {
        return(ip_ex_tt_phi_tt_ATLAS13_cache[NumPar][i] );
    } else {
        double newResult = interpolate (ATLAS13_tt_phi_tt,mass);
        CacheShiftReal(ip_ex_tt_phi_tt_ATLAS13_cache, NumPar, params, newResult);
        return newResult;
    }
}



double THDMWcache::ip_ex_tt_phi_tt_ATLAS13_e(double mass){
    int NumPar = 1;
    double params[] = {mass};

    int i = CacheCheckReal(ip_ex_tt_phi_tt_ATLAS13_cache_e, NumPar, params);
    if (i>=0) {
        return(ip_ex_tt_phi_tt_ATLAS13_cache_e[NumPar][i] );
    } else {
        double newResult = interpolate (ATLAS13_tt_phi_tt_e,mass);
        CacheShiftReal(ip_ex_tt_phi_tt_ATLAS13_cache_e, NumPar, params, newResult);
        return newResult;
    }
}

double THDMWcache::ip_ex_pp_H_hh_bbbb_ATLAS13(double mass){
    int NumPar = 1;
    double params[] = {mass};

    int i = CacheCheckReal(ip_ex_pp_H_hh_bbbb_ATLAS13_cache, NumPar, params);
    if (i>=0) {
        return(ip_ex_pp_H_hh_bbbb_ATLAS13_cache[NumPar][i] );
    } else {
        double newResult = interpolate (ATLAS13_pp_H_hh_bbbb,mass);
        CacheShiftReal(ip_ex_pp_H_hh_bbbb_ATLAS13_cache, NumPar, params, newResult);
        return newResult;
    }
}



double THDMWcache::ip_ex_pp_H_hh_bbbb_ATLAS13_e(double mass){
    int NumPar = 1;
    double params[] = {mass};

    int i = CacheCheckReal(ip_ex_pp_H_hh_bbbb_ATLAS13_cache_e, NumPar, params);
    if (i>=0) {
        return(ip_ex_pp_H_hh_bbbb_ATLAS13_cache_e[NumPar][i] );
    } else {
        double newResult = interpolate (ATLAS13_pp_H_hh_bbbb_e,mass);
        CacheShiftReal(ip_ex_pp_H_hh_bbbb_ATLAS13_cache_e, NumPar, params, newResult);
        return newResult;
    }
}


double THDMWcache::ip_ex_pp_phi_bb_CMS13(double mass){
    int NumPar = 1;
    double params[] = {mass};

    int i = CacheCheckReal(ip_ex_pp_phi_bb_CMS13_cache, NumPar, params);
    if (i>=0) {
        return(ip_ex_pp_phi_bb_CMS13_cache[NumPar][i] );
    } else {
        double newResult = interpolate (CMS13_pp_phi_bb,mass);
        CacheShiftReal(ip_ex_pp_phi_bb_CMS13_cache, NumPar, params, newResult);
        return newResult;
    }
}

double THDMWcache::ip_ex_pp_phi_bb_CMS8(double mass){
    int NumPar = 1;
    double params[] = {mass};

    int i = CacheCheckReal(ip_ex_pp_phi_bb_CMS8_cache, NumPar, params);
    if (i>=0) {
        return(ip_ex_pp_phi_bb_CMS8_cache[NumPar][i] );
    } else {
        double newResult = interpolate (CMS8_pp_phi_bb,mass);
        CacheShiftReal(ip_ex_pp_phi_bb_CMS8_cache, NumPar, params, newResult);
        return newResult;
    }
}





double THDMWcache::ip_ex_pp_phi_bb_CMS13_e(double mass){
    int NumPar = 1;
    double params[] = {mass};

    int i = CacheCheckReal(ip_ex_pp_phi_bb_CMS13_cache_e, NumPar, params);
    if (i>=0) {
        return(ip_ex_pp_phi_bb_CMS13_cache_e[NumPar][i] );
    } else {
        double newResult = interpolate (CMS13_pp_phi_bb_e,mass);
        CacheShiftReal(ip_ex_pp_phi_bb_CMS13_cache_e, NumPar, params, newResult);
        return newResult;
    }
}

double THDMWcache::ip_ex_pp_H_hh_bbbb_CMS13(double mass){
    int NumPar = 1;
    double params[] = {mass};

    int i = CacheCheckReal(ip_ex_pp_H_hh_bbbb_CMS13_cache, NumPar, params);
    if (i>=0) {
        return(ip_ex_pp_H_hh_bbbb_CMS13_cache[NumPar][i] );
    } else {
        double newResult = interpolate (CMS13_pp_H_hh_bbbb,mass);
        CacheShiftReal(ip_ex_pp_H_hh_bbbb_CMS13_cache, NumPar, params, newResult);
        return newResult;
    }
}



double THDMWcache::ip_ex_pp_H_hh_bbbb_CMS13_e(double mass){
    int NumPar = 1;
    double params[] = {mass};

    int i = CacheCheckReal(ip_ex_pp_H_hh_bbbb_CMS13_cache_e, NumPar, params);
    if (i>=0) {
        return(ip_ex_pp_H_hh_bbbb_CMS13_cache_e[NumPar][i] );
    } else {
        double newResult = interpolate (CMS13_pp_H_hh_bbbb_e,mass);
        CacheShiftReal(ip_ex_pp_H_hh_bbbb_CMS13_cache_e, NumPar, params, newResult);
        return newResult;
    }
}


double THDMWcache::ip_ex_pp_Hpm_tb_ATLAS8(double mass){
    int NumPar = 1;
    double params[] = {mass};

    int i = CacheCheckReal(ip_ex_pp_Hpm_tb_ATLAS8_cache, NumPar, params);
    if (i>=0) {
        return(ip_ex_pp_Hpm_tb_ATLAS8_cache[NumPar][i] );
    } else {
        double newResult = interpolate (ATLAS8_pp_Hpm_tb,mass);
        CacheShiftReal(ip_ex_pp_Hpm_tb_ATLAS8_cache, NumPar, params, newResult);
        return newResult;
    }
}



double THDMWcache::ip_ex_pp_Hpm_tb_ATLAS8_e(double mass){
    int NumPar = 1;
    double params[] = {mass};

    int i = CacheCheckReal(ip_ex_pp_Hpm_tb_ATLAS8_cache_e, NumPar, params);
    if (i>=0) {
        return(ip_ex_pp_Hpm_tb_ATLAS8_cache_e[NumPar][i] );
    } else {
        double newResult = interpolate (ATLAS8_pp_Hpm_tb_e,mass);
        CacheShiftReal(ip_ex_pp_Hpm_tb_ATLAS8_cache_e, NumPar, params, newResult);
        return newResult;
    }
}



double THDMWcache::ip_ex_pp_Hp_tb_CMS8(double mass){
    int NumPar = 1;
    double params[] = {mass};

    int i = CacheCheckReal(ip_ex_pp_Hp_tb_CMS8_cache, NumPar, params);
    if (i>=0) {
        return(ip_ex_pp_Hp_tb_CMS8_cache[NumPar][i] );
    } else {
        double newResult = interpolate (CMS8_pp_Hp_tb,mass);
        CacheShiftReal(ip_ex_pp_Hp_tb_CMS8_cache, NumPar, params, newResult);
        return newResult;
    }
}



double THDMWcache::ip_ex_pp_Hp_tb_CMS8_e(double mass){
    int NumPar = 1;
    double params[] = {mass};

    int i = CacheCheckReal(ip_ex_pp_Hp_tb_CMS8_cache_e, NumPar, params);
    if (i>=0) {
        return(ip_ex_pp_Hp_tb_CMS8_cache_e[NumPar][i] );
    } else {
        double newResult = interpolate (CMS8_pp_Hp_tb_e,mass);
        CacheShiftReal(ip_ex_pp_Hp_tb_CMS8_cache_e, NumPar, params, newResult);
        return newResult;
    }
}

double THDMWcache::ip_ex_pp_Hp_tb_ATLAS13(double mass){
    int NumPar = 1;
    double params[] = {mass};

    int i = CacheCheckReal(ip_ex_pp_Hp_tb_ATLAS13_cache, NumPar, params);
    if (i>=0) {
        return(ip_ex_pp_Hp_tb_ATLAS13_cache[NumPar][i] );
    } else {
        double newResult = interpolate (ATLAS13_pp_Hp_tb,mass);
        CacheShiftReal(ip_ex_pp_Hp_tb_ATLAS13_cache, NumPar, params, newResult);
        return newResult;
    }
}




double THDMWcache::ip_ex_ggF_H_hh_bbbb_CMS13(double mass){
    int NumPar = 1;
    double params[] = {mass};

    int i = CacheCheckReal(ip_ex_ggF_H_hh_bbbb_CMS13_cache, NumPar, params);
    if (i>=0) {
        return(ip_ex_ggF_H_hh_bbbb_CMS13_cache[NumPar][i] );
    } else {
        double newResult = interpolate (CMS13_ggF_H_hh_bbbb,mass);
        CacheShiftReal(ip_ex_ggF_H_hh_bbbb_CMS13_cache, NumPar, params, newResult);
        return newResult;
    }
}



double THDMWcache::ip_ex_ggF_H_hh_bbbb_CMS13_e(double mass){
    int NumPar = 1;
    double params[] = {mass};

    int i = CacheCheckReal(ip_ex_ggF_H_hh_bbbb_CMS13_cache_e, NumPar, params);
    if (i>=0) {
        return(ip_ex_ggF_H_hh_bbbb_CMS13_cache_e[NumPar][i] );
    } else {
        double newResult = interpolate (CMS13_ggF_H_hh_bbbb_e,mass);
        CacheShiftReal(ip_ex_ggF_H_hh_bbbb_CMS13_cache_e, NumPar, params, newResult);
        return newResult;
    }
}

double THDMWcache::ip_ex_pp_Gkk_tt_ATLAS13(double mass){
    int NumPar = 1;
    double params[] = {mass};

    int i = CacheCheckReal(ip_ex_pp_Gkk_tt_ATLAS13_cache, NumPar, params);
    if (i>=0) {
        return(ip_ex_pp_Gkk_tt_ATLAS13_cache[NumPar][i] );
    } else {
        double newResult = interpolate(ATLAS13_pp_Gkk_tt,mass);
        CacheShiftReal(ip_ex_pp_Gkk_tt_ATLAS13_cache, NumPar, params, newResult);
        return newResult;
    }
}

double THDMWcache::ip_ex_pp_R_gg_CMS13(double mass){
    int NumPar = 1;
    double params[] = {mass};

    int i = CacheCheckReal(ip_ex_pp_R_gg_CMS13_cache, NumPar, params);
    if (i>=0) {
        return(ip_ex_pp_R_gg_CMS13_cache[NumPar][i] );
    } else {
        double newResult = interpolate(CMS13_pp_R_gg,mass);
        CacheShiftReal(ip_ex_pp_R_gg_CMS13_cache, NumPar, params, newResult);
        return newResult;
    }
}

double THDMWcache::ip_ex_pp_SS_jjjj_ATLAS13(double mass){
    int NumPar = 1;
    double params[] = {mass};

    int i = CacheCheckReal(ip_ex_pp_SS_jjjj_ATLAS13_cache, NumPar, params);
    if (i>=0) {
        return(ip_ex_pp_SS_jjjj_ATLAS13_cache[NumPar][i] );
    } else {
        double newResult = interpolate(ATLAS13_pp_SS_jjjj,mass);
        CacheShiftReal(ip_ex_pp_SS_jjjj_ATLAS13_cache, NumPar, params, newResult);
        return newResult;
    }
}

double THDMWcache::ip_ex_bb_H_bb_CMS13(double mass){
    int NumPar = 1;
    double params[] = {mass};

    int i = CacheCheckReal(ip_ex_bb_H_bb_CMS13_cache, NumPar, params);
    if (i>=0) {
        return ( ip_ex_bb_H_bb_CMS13_cache[NumPar][i] );
    } else {
        double newResult = interpolate(CMS13_bb_H_bb,mass);
        CacheShiftReal(ip_ex_bb_H_bb_CMS13_cache, NumPar, params, newResult);
        return newResult;
    }
}





double THDMWcache::ip_th_pp_Sr_tt(double etaD, double etaU, double Lambda4, double mSr){
    int NumPar = 4;
    double params[] = {etaD, etaU, Lambda4, mSr};

    int i = CacheCheckReal(ip_th_pp_Sr_tt_cache, NumPar, params);
    if (i>=0) {
        return(ip_th_pp_Sr_tt_cache[NumPar][i] );
    } else {
        double newResult = interpolate4D (MadGraph_pp_Sr_tt,etaD,etaU,Lambda4,mSr);
        CacheShiftReal(ip_th_pp_Sr_tt_cache, NumPar, params, newResult);
        return newResult;
    }
}

double THDMWcache::ip_th_pp_Srtt_tttt(double etaD, double etaU, double Lambda4, double mSr){
    int NumPar = 4;
    double params[] = {etaD, etaU, Lambda4, mSr};

    int i = CacheCheckReal(ip_th_pp_Srtt_tttt_cache, NumPar, params);
    if (i>=0) {
        return(ip_th_pp_Srtt_tttt_cache[NumPar][i] );
    } else {
        double newResult = interpolate4D (MadGraph_pp_Srtt_tttt,etaD,etaU,Lambda4,mSr);
        CacheShiftReal(ip_th_pp_Srtt_tttt_cache, NumPar, params, newResult);
        return newResult;
    }
}

double THDMWcache::ip_th_pp_Sr_jj(double etaD, double etaU, double Lambda4, double mSr){
    int NumPar = 4;
    double params[] = {etaD, etaU, Lambda4, mSr};

    int i = CacheCheckReal(ip_th_pp_Sr_jj_cache, NumPar, params);
    if (i>=0) {
        return(ip_th_pp_Sr_jj_cache[NumPar][i] );
    } else {
        double newResult = interpolate4D (MadGraph_pp_Sr_jj,etaD,etaU,Lambda4,mSr);
        CacheShiftReal(ip_th_pp_Sr_jj_cache, NumPar, params, newResult);
        return newResult;
    }
}

double THDMWcache::ip_th_pp_SrSr_jjjj(double etaD, double etaU, double Lambda4, double mSr){
    int NumPar = 4;
    double params[] = {etaD, etaU, Lambda4, mSr};

    int i = CacheCheckReal(ip_th_pp_SrSr_jjjj_cache, NumPar, params);
    if (i>=0) {
        return(ip_th_pp_SrSr_jjjj_cache[NumPar][i] );
    } else {
        double newResult = interpolate4D (MadGraph_pp_SrSr_jjjj,etaD,etaU,Lambda4,mSr);
        CacheShiftReal(ip_th_pp_SrSr_jjjj_cache, NumPar, params, newResult);
        return newResult;
    }
}


double THDMWcache::ip_th_pp_Stb_tbtb(double etaD, double etaU, double mS){
    int NumPar = 3;
    double params[] = {etaD, etaU, mS};

    int i = CacheCheckReal(ip_th_pp_Stb_tbtb_cache, NumPar, params);
    if (i>=0) {
        return(ip_th_pp_Stb_tbtb_cache[NumPar][i] );
    } else {
        double newResult = interpolate3D (MadGraph_pp_Stb_tbtb,etaD,etaU,mS);
        CacheShiftReal(ip_th_pp_Stb_tbtb_cache, NumPar, params, newResult);
        return newResult;
    }
}


double THDMWcache::ip_th_pp_Sitt_tttt(double etaD, double etaU, double mS){
    int NumPar = 3;
    double params[] = {etaD, etaU, mS};

    int i = CacheCheckReal(ip_th_pp_Sitt_tttt_cache, NumPar, params);
    if (i>=0) {
        return(ip_th_pp_Sitt_tttt_cache[NumPar][i] );
    } else {
        double newResult = interpolate3D (MadGraph_pp_Sitt_tttt,etaD,etaU,mS);
        CacheShiftReal(ip_th_pp_Sitt_tttt_cache, NumPar, params, newResult);
        return newResult;
    }
}


double THDMWcache::ip_th_pp_Srbb_bbbb(double etaD, double etaU, double Lambda4, double mSr){
    int NumPar = 4;
    double params[] = {etaD, etaU, Lambda4, mSr};

    int i = CacheCheckReal(ip_th_pp_Srbb_bbbb_cache, NumPar, params);
    if (i>=0) {
        return(ip_th_pp_Srbb_bbbb_cache[NumPar][i] );
    } else {
        double newResult = interpolate4D (MadGraph_pp_Srbb_bbbb,etaD,etaU,Lambda4,mSr);
        CacheShiftReal(ip_th_pp_Srbb_bbbb_cache, NumPar, params, newResult);
        return newResult;
    }
}


double THDMWcache::ip_th_pp_Srbb_bbbb_8TeV(double etaD, double etaU, double Lambda4, double mSr){
    int NumPar = 4;
    double params[] = {etaD, etaU, Lambda4, mSr};

    int i = CacheCheckReal(ip_th_pp_Srbb_bbbb_8TeV_cache, NumPar, params);
    if (i>=0) {
        return(ip_th_pp_Srbb_bbbb_8TeV_cache[NumPar][i] );
    } else {
        double newResult = interpolate4D (MadGraph_pp_Srbb_bbbb_8TeV,etaD,etaU,Lambda4,mSr);
        CacheShiftReal(ip_th_pp_Srbb_bbbb_8TeV_cache, NumPar, params, newResult);
        return newResult;
    }
}


double THDMWcache::ip_th_pp_Sibb_bbbb(double etaD, double etaU, double mS){
    int NumPar = 3;
    double params[] = {etaD, etaU, mS};

    int i = CacheCheckReal(ip_th_pp_Sibb_bbbb_cache, NumPar, params);
    if (i>=0) {
        return(ip_th_pp_Sibb_bbbb_cache[NumPar][i] );
    } else {
        double newResult = interpolate3D (MadGraph_pp_Sibb_bbbb,etaD,etaU,mS);
        //std::cout<<"check"<<std::endl;
        CacheShiftReal(ip_th_pp_Sibb_bbbb_cache, NumPar, params, newResult);
        return newResult;
    }
}



double THDMWcache::ip_th_pp_Sibb_bbbb_8TeV(double etaD, double etaU, double mS){
    int NumPar = 3;
    double params[] = {etaD, etaU, mS};

    int i = CacheCheckReal(ip_th_pp_Sibb_bbbb_8TeV_cache, NumPar, params);
    if (i>=0) {
        return(ip_th_pp_Sibb_bbbb_8TeV_cache[NumPar][i] );
    } else {
        double newResult = interpolate3D (MadGraph_pp_Sibb_bbbb_8TeV,etaD,etaU,mS);
        //std::cout<<"check"<<std::endl;
        CacheShiftReal(ip_th_pp_Sibb_bbbb_8TeV_cache, NumPar, params, newResult);
        return newResult;
    }
}





double THDMWcache::ip_th_pp_Sr_bb(double etaD, double etaU, double Lambda4, double mSr){
    int NumPar = 4;
    double params[] = {etaD, etaU, Lambda4, mSr};

    int i = CacheCheckReal(ip_th_pp_Sr_bb_cache, NumPar, params);
    if (i>=0) {
        return(ip_th_pp_Sr_bb_cache[NumPar][i] );
    } else {
        double newResult = interpolate4D (MadGraph_pp_Sr_bb,etaD,etaU,Lambda4,mSr);
        CacheShiftReal(ip_th_pp_Sr_bb_cache, NumPar, params, newResult);
        return newResult;
    }
}



double THDMWcache::ip_th_pp_Sr_bb_8TeV(double etaD, double etaU, double Lambda4, double mSr){
    int NumPar = 4;
    double params[] = {etaD, etaU, Lambda4, mSr};

    int i = CacheCheckReal(ip_th_pp_Sr_bb_8TeV_cache, NumPar, params);
    if (i>=0) {
        return(ip_th_pp_Sr_bb_8TeV_cache[NumPar][i] );
    } else {
        double newResult = interpolate4D (MadGraph_pp_Sr_bb_8TeV,etaD,etaU,Lambda4,mSr);
        CacheShiftReal(ip_th_pp_Sr_bb_8TeV_cache, NumPar, params, newResult);
        return newResult;
    }
}


double THDMWcache::ip_th_pp_Si_bb(double etaD, double etaU, double mS){
    int NumPar = 3;
    double params[] = {etaD, etaU, mS};

    int i = CacheCheckReal(ip_th_pp_Si_bb_cache, NumPar, params);
    if (i>=0) {
        return(ip_th_pp_Si_bb_cache[NumPar][i] );
    } else {
        double newResult = interpolate3D (MadGraph_pp_Si_bb,etaD,etaU,mS);
        //std::cout<<"check"<<std::endl;
        CacheShiftReal(ip_th_pp_Si_bb_cache, NumPar, params, newResult);
        return newResult;
    }
}




double THDMWcache::ip_th_pp_Si_bb_8TeV(double etaD, double etaU, double mS){
    int NumPar = 3;
    double params[] = {etaD, etaU, mS};

    int i = CacheCheckReal(ip_th_pp_Si_bb_8TeV_cache, NumPar, params);
    if (i>=0) {
        return(ip_th_pp_Si_bb_8TeV_cache[NumPar][i] );
    } else {
        double newResult = interpolate3D (MadGraph_pp_Si_bb_8TeV,etaD,etaU,mS);
        //std::cout<<"check"<<std::endl;
        CacheShiftReal(ip_th_pp_Si_bb_8TeV_cache, NumPar, params, newResult);
        return newResult;
    }
}




/*
double THDMWcache::logip_th_pp_SrSr_jjjj(double etaD, double etaU, double Lambda4, double mSr){
    int NumPar = 4;
    double params[] = {etaD, etaU, Lambda4, mSr};

    int i = CacheCheckReal(logip_th_pp_SrSr_jjjj_cache, NumPar, params);
    if (i>=0) {
        return(logip_th_pp_SrSr_jjjj_cache[NumPar][i] );
    } else {
        double newResult = loginterpolate4D (MadGraph_pp_SrSr_jjjj,etaD,etaU,Lambda4,mSr);
        CacheShiftReal(logip_th_pp_SrSr_jjjj_cache, NumPar, params, newResult);
        return newResult;
    }
}
*/



void THDMWcache::computeHHlimits()
{
    double mSr=sqrt(mSRsq);
    double mSp=sqrt(mSpsq);
    double mSi=sqrt(mSIsq);
    double MW=myTHDMW->Mw();
    //double mW2=myTHDMW->Mw();
    double SqrtEtaU=copysign(sqrt(sqrt(pow(etaU,2))),etaU);
    double SqrtEtaD=copysign(sqrt(sqrt(pow(etaD,2))),etaD);
    double nu45=(nu4+nu5)/2;
    //EtaU and EtaD in Sqrt Units!!!
    THoEX_pp_Sr_tt=0.;
    THoEX_pp_Srtt_tttt=0.;
    THoEX_pp_Sr_jj=0.;
    THoEX_pp_SrSr_jjjj=0.;
    THoEX_pp_Stb_tbtb=0.;
    THoEX_pp_Sitt_tttt=0.;
    THoEX_pp_Srbb_bbbb=0.;
    THoEX_pp_Srbb_bbbb_8TeV=0.;
    THoEX_pp_Sibb_bbbb=0.;
    THoEX_pp_Sibb_bbbb_8TeV=0.;
    THoEX_pp_Sr_bb=0.;
    THoEX_pp_Sr_bb_8TeV=0.;
    THoEX_pp_Si_bb=0.;
    THoEX_pp_Si_bb_8TeV=0.;
    
    pp_Sr_tt_TH13 = 1.0e-15;
    pp_Srtt_tttt_TH13 = 1.0e-15;
    pp_Sr_jj_TH13=1.0e-15;
    pp_SrSr_jjjj_TH13=1.0e-15;
    pp_Stb_tbtb_TH13=1.0e-15;
    pp_Srtt_tttt_TH13 = 1.0e-15;
    pp_Sitt_tttt_TH13 = 1.0e-15;
    pp_Srbb_bbbb_TH13= 1.0e-15;
    pp_Srbb_bbbb_TH8= 1.0e-15;
    pp_Sibb_bbbb_TH13= 1.0e-15;
    pp_Sibb_bbbb_TH8= 1.0e-15;
    pp_Sr_bb_TH13= 1.0e-15;
    pp_Sr_bb_TH8= 1.0e-15;
    pp_Si_bb_TH13= 1.0e-15;
    pp_Si_bb_TH8= 1.0e-15;
    //logpp_SrSr_jjjj_TH13=-15;
    
    
    
    //std::cout<<"mSr="<<mSr<<std::endl;
    //std::cout<<"nu45="<<nu45<<std::endl;
    //std::cout<<"etaU="<<etaU<<std::endl;
    //std::cout<<"etaD="<<etaD<<std::endl;
    //std::cout<<"MW="<<MW<<std::endl;
    //std::cout<<"MZ="<<MZ<<std::endl;
    
    if(mSRsq > 1.6001e5 && mSRsq<2.2499e6 && etaD*etaD<399.99 && etaU*etaU<56.2499 && nu45*nu45<168.999 &&  mSRsq<(mSpsq+MW*MW) && mSRsq<=(mSIsq+MZ*MZ)) pp_Sr_tt_TH13=ip_th_pp_Sr_tt(SqrtEtaD,SqrtEtaU,nu45,mSr);
    if(mSRsq > 1.6001e5 && mSRsq<2.2499e6 && etaD*etaD<399.99 && etaU*etaU<56.2499 && nu45*nu45<168.999 &&  mSRsq<(mSpsq+MW*MW) && mSRsq<=(mSIsq+MZ*MZ)) pp_Srtt_tttt_TH13=ip_th_pp_Srtt_tttt(SqrtEtaD,SqrtEtaU,nu45,mSr);
    if(mSRsq > 1.6001e5 && mSRsq<2.2499e6 && etaD*etaD<399.99 && etaU*etaU<3.9999 && nu45*nu45<168.999 &&  mSRsq<(mSpsq+MW*MW) && mSRsq<=(mSIsq+MZ*MZ)) pp_Sr_jj_TH13=ip_th_pp_Sr_jj(SqrtEtaD,SqrtEtaU,nu45,mSr);
    if(mSRsq > 1.6001e5 && mSRsq<2.2499e6 && etaD*etaD<399.99 && etaU*etaU<3.9999 && nu45*nu45<168.999 &&  mSRsq<(mSpsq+MW*MW) && mSRsq<=(mSIsq+MZ*MZ)) pp_SrSr_jjjj_TH13=ip_th_pp_SrSr_jjjj(SqrtEtaD,SqrtEtaU,nu45,mSr);
    if(mSpsq > 1.6001e5 && mSpsq<2.2499e6 && etaD*etaD<399.99 && etaU*etaU<56.2499 && mSpsq<(mSRsq+MW*MW) && mSp<=(mSIsq+MW*MW)) pp_Stb_tbtb_TH13=ip_th_pp_Stb_tbtb(SqrtEtaD,SqrtEtaU,mSp);
    if(mSIsq > 1.6001e5 && mSIsq<2.2499e6 && etaD*etaD<399.99 && etaU*etaU<56.2499 && mSIsq<(mSpsq+MW*MW) && mSIsq<=(mSRsq+MZ*MZ)) pp_Sitt_tttt_TH13=ip_th_pp_Sitt_tttt(SqrtEtaD,SqrtEtaU,mSi);
    if(mSRsq > 1.6001e5 && mSRsq<2.2499e6 && etaD*etaD<399.99 && etaU*etaU<3.9999 && nu45*nu45<168.999 &&  mSRsq<(mSpsq+MW*MW) && mSRsq<=(mSIsq+MZ*MZ)) pp_Srbb_bbbb_TH13=ip_th_pp_Srbb_bbbb(SqrtEtaD,SqrtEtaU,nu45,mSr);
    if(mSRsq > 1.6001e5 && mSRsq<2.2499e6 && etaD*etaD<399.99 && etaU*etaU<3.9999 && nu45*nu45<168.999 &&  mSRsq<(mSpsq+MW*MW) && mSRsq<=(mSIsq+MZ*MZ)) pp_Srbb_bbbb_TH8=ip_th_pp_Srbb_bbbb_8TeV(SqrtEtaD,SqrtEtaU,nu45,mSr);
    if(mSIsq > 1.6001e5 && mSIsq<2.2499e6 && etaD*etaD<399.99 && etaU*etaU<3.999 &&  mSIsq<(mSpsq+MW*MW) && mSIsq<=(mSRsq+MZ*MZ)) pp_Sibb_bbbb_TH13=ip_th_pp_Sibb_bbbb(SqrtEtaD,SqrtEtaU,mSi);
    if(mSIsq > 1.6001e5 && mSIsq<2.2499e6 && etaD*etaD<399.99 && etaU*etaU<3.999 &&  mSIsq<(mSpsq+MW*MW) && mSIsq<=(mSRsq+MZ*MZ)) pp_Sibb_bbbb_TH8=ip_th_pp_Sibb_bbbb_8TeV(SqrtEtaD,SqrtEtaU,mSi);
    if(mSRsq > 1.6001e5 && mSRsq<2.2499e6 && etaD*etaD<399.99 && etaU*etaU<3.9999 && nu45*nu45<168.999 &&  mSRsq<(mSpsq+MW*MW) && mSRsq<=(mSIsq+MZ*MZ)) pp_Sr_bb_TH13=ip_th_pp_Sr_bb(SqrtEtaD,SqrtEtaU,nu45,mSr);
    if(mSRsq > 1.6001e5 && mSRsq<2.2499e6 && etaD*etaD<399.99 && etaU*etaU<3.9999 && nu45*nu45<168.999 &&  mSRsq<(mSpsq+MW*MW) && mSRsq<=(mSIsq+MZ*MZ)) pp_Sr_bb_TH8 =ip_th_pp_Sr_bb_8TeV(SqrtEtaD,SqrtEtaU,nu45,mSr);
    if(mSIsq > 1.6001e5 && mSIsq<2.2499e6 && etaD*etaD<399.99 && etaU*etaU<3.999 &&  mSIsq<(mSpsq+MW*MW) && mSIsq<=(mSRsq+MZ*MZ)) pp_Si_bb_TH13=ip_th_pp_Si_bb(SqrtEtaD,SqrtEtaU,mSi);
    if(mSIsq > 1.6001e5 && mSIsq<2.2499e6 && etaD*etaD<399.99 && etaU*etaU<3.999 &&  mSIsq<(mSpsq+MW*MW) && mSIsq<=(mSRsq+MZ*MZ)) pp_Si_bb_TH8=ip_th_pp_Si_bb_8TeV(SqrtEtaD,SqrtEtaU,mSi);



    //std::cout<<"Sibb2="<< interpolate3D(MadGraph_pp_Sibb_bbbb,4.47213,1.4141,1450)<<std::endl;
    //std::cout<<"Sibb="<< interpolate3D(MadGraph_pp_Sibb_bbbb,-4.47214,-1.4142,350)<<std::endl;
    //std::cout<<"Sitt="<< interpolate3D(MadGraph_pp_Sitt_tttt,SqrtEtaD,SqrtEtaU,mSi)<<std::endl;
    //std::cout<<"Sitt="<<ip_th_pp_Sitt_tttt(SqrtEtaD,0,mSi)<<std::endl;
    //if(mSr>= 400 && mSr<=1500) logpp_SrSr_jjjj_TH13=logip_th_pp_SrSr_jjjj(SqrtEtaD,SqrtEtaU,nu45,mSr);

    //std::cout<<"pp_Stb_tbtb_TH13 first="<<ip_th_pp_Stb_tbtb(SqrtEtaD,SqrtEtaU,mSp)<<std::endl;
    //std::cout<<"pp_Stb_tbtb_TH13 second="<<interpolate3D(MadGraph_pp_Stb_tbtb,SqrtEtaD,SqrtEtaU,mSp)<<std::endl;
    if(mSr>= 400 && mSr<=1500) THoEX_pp_Sr_tt=pp_Sr_tt_TH13/ip_ex_pp_Gkk_tt_ATLAS13(mSr);
    if(mSr>= 400 && mSr<=1000) THoEX_pp_Srtt_tttt=pp_Srtt_tttt_TH13/ip_ex_tt_phi_tt_ATLAS13(mSr);
    if(mSr>= 600 && mSr<=1500) THoEX_pp_Sr_jj=pp_Sr_jj_TH13/ip_ex_pp_R_gg_CMS13(mSr);
    if(mSr>= 500 && mSr<=1500) THoEX_pp_SrSr_jjjj=pp_SrSr_jjjj_TH13/ip_ex_pp_SS_jjjj_ATLAS13(mSr);
    if(mSp>= 400 && mSp<=1500) THoEX_pp_Stb_tbtb=pp_Stb_tbtb_TH13/ip_ex_pp_Hp_tb_ATLAS13(mSp);
    if(mSi>= 400 && mSi<=1000) THoEX_pp_Sitt_tttt=pp_Sitt_tttt_TH13/ip_ex_tt_phi_tt_ATLAS13(mSi);
    if(mSr>= 400 && mSr<=1300) THoEX_pp_Srbb_bbbb=pp_Srbb_bbbb_TH13/ip_ex_bb_H_bb_CMS13(mSr);
    if(mSr>= 400 && mSr<=900) THoEX_pp_Srbb_bbbb_8TeV=pp_Srbb_bbbb_TH8/ip_ex_bb_phi_bb_CMS8(mSr);
    if(mSi>= 400 && mSi<=1300) THoEX_pp_Sibb_bbbb=pp_Sibb_bbbb_TH13/ip_ex_bb_H_bb_CMS13(mSi);
    if(mSi>= 400 && mSi<=900) THoEX_pp_Sibb_bbbb_8TeV=pp_Sibb_bbbb_TH8/ip_ex_bb_phi_bb_CMS8(mSi);
    if(mSr>= 550 && mSr<=1200) THoEX_pp_Sr_bb=pp_Sr_bb_TH13/ip_ex_pp_phi_bb_CMS13(mSr);
    if(mSr>= 400 && mSr<=1200) THoEX_pp_Sr_bb_8TeV=pp_Sr_bb_TH8/ip_ex_pp_phi_bb_CMS8(mSr);
    if(mSi>= 550 && mSi<=1200) THoEX_pp_Si_bb=pp_Si_bb_TH13/ip_ex_pp_phi_bb_CMS13(mSi);
    if(mSr>= 400 && mSr<=1200) THoEX_pp_Si_bb_8TeV=pp_Si_bb_TH8/ip_ex_pp_phi_bb_CMS8(mSr);
    //std::cout<<"ip_ex_pp_phi_bb_CMS13(mSi)="<< ip_ex_pp_phi_bb_CMS13(550) <<std::endl;
    //std::cout<<"pp_Si_bb_TH13="<< ip_th_pp_Si_bb(4.469,0,550) <<std::endl;
    //std::cout<<"pp_Si_bb_TH13="<< ip_th_pp_Si_bb(4.469,1.41,550) <<std::endl;
}



double THDMWcache::setOtherParameters()
{
    double sin2b=2.0*sinb*cosb;
    double cos2b=cosb*cosb-sinb*sinb;
    double tan2b=sin2b/cos2b;
    double cot2b=1.0/tan2b;
    double sin2a=2.0*sina*cosa;
    double cos2a=cosa*cosa-sina*sina;
    double tan2a=sin2a/cos2a;
    double cot2a=1.0/tan2a;
    double lambda345=lambda3+lambda4+lambda5;

    m11sq = vev*vev*(lambda2*sinb*sinb*tanb/(cot2a-2.0*cot2b)
                     +(lambda1*(cosb*cosb - (4.0*cosb*cosb-3.0)*cosb*tan2a/sinb)
                       -lambda345*(sinb*sinb + cos2b*tan2a*tanb))/(4.0*cot2b*tan2a-2.0));

    m22sq = vev*vev*(-lambda1*cosb*cosb*cosb/sinb/(cot2a-2.0*cot2b)
                     +(lambda2*(sinb*sinb + (4.0*sinb*sinb-3.0)*tanb*tan2a)
                       -lambda345*(cosb*cosb + cos2b*tan2a*cosb/sinb))/(4.0*cot2b*tan2a-2.0));

    m12sq = vev*vev*(-lambda345*sin2b
                     +2.0*(lambda1*cosb*cosb - lambda2*sinb*sinb)*tan2a/(4.0*tan2a/tan2b-2.0));
    
    mHsq = vev*vev*(lambda1*cosa*cosa*cosb*cosb + lambda2*sina*sina*sinb*sinb
                    +lambda345*sin2a*cosb*sinb
                    +sin(bma)*sin(bma)*(lambda345 + (lambda2 - lambda1/(tanb*tanb))*tan2a*tanb)/(1.0 - 2.0*cot2b*tan2a));

    mAsq = vev*vev*(lambda3+lambda4 + tan2a*(-lambda1*cosb/sinb + lambda2*tanb + 2.0*lambda5*cot2b))/(1.0 - 2.0*cot2b*tan2a);

    mhsq = vev*vev*(lambda1*sina*sina*cosb*cosb + lambda2*cosa*cosa*sinb*sinb
                    -lambda345*sin2a*cosb*sinb
                    +cos(bma)*cos(bma)*(lambda345 + (lambda2 - lambda1/(tanb*tanb))*tan2a*tanb)/(1.0 - 2.0*cot2b*tan2a));

    mSRsq = mSsq + vev*vev*((nu1+nu2+2.0*nu3)*cosb*cosb + (omega1+omega2+2.0*omega3)*sinb*sinb
                            +(kappa1+kappa2+kappa3)*sin2b)/4.0;

    mSIsq = mSsq + vev*vev*((nu1+nu2-2.0*nu3)*cosb*cosb + (omega1+omega2-2.0*omega3)*sinb*sinb
                            +(kappa1+kappa2-kappa3)*sin2b)/4.0;

    if( THDMWmodel == "custodial1" ) {
        mHpsq = mAsq;
        mSpsq = mSIsq;
    }
    else if( THDMWmodel == "ManoharWise" ) {
        mhsq = vev*vev*lambda1;
        mSRsq = mSsq + vev*vev*(nu1+nu2+2.0*nu3)/4.0;
        mSIsq = mSsq + vev*vev*(nu1+nu2-2.0*nu3)/4.0;
        mSpsq = mSsq + vev*vev*nu1/4.0;
    }
    else if( THDMWmodel == "custodialMW" ) {
        mhsq = vev*vev*lambda1;
        mSRsq = mSsq + vev*vev*(nu1+2.0*nu2)/4.0;
        mSIsq = mSsq + vev*vev*nu1/4.0;
        mSpsq = mSIsq;
    }
    else {
        mHpsq = vev*vev*(lambda345 + tan2a*(-lambda1*cosb/sinb + lambda2*tanb + (lambda4+lambda5)*cot2b))/(1.0 - 2.0*cot2b*tan2a);
        mSpsq = mSsq + vev*vev*(nu1*cosb*cosb + omega1*sinb*sinb + kappa1*sin2b)/4.0;
    }

    if(mhsq < 0 || mHsq < 0 || mAsq < 0 || mSRsq < 0 || mSIsq < 0 || mHpsq < 0 || mSpsq < 0)
    {
        return std::numeric_limits<double>::quiet_NaN();
    }
    else
    {
        return 1.;
    }
}

void THDMWcache::updateCache()
{
    THDMWmodel=myTHDMW->getModelTypeTHDMWflag();
    Q_THDMW=myTHDMW->getQ_THDMW();
    MZ=myTHDMW->getMz();
    vev=myTHDMW->v();
    tanb=myTHDMW->getTHDMW_tanb();
    sinb=myTHDMW->getTHDMW_sinb();
    cosb=myTHDMW->getTHDMW_cosb();
    bma=myTHDMW->getTHDMW_bma();
    sina=myTHDMW->getTHDMW_sina();
    cosa=myTHDMW->getTHDMW_cosa();
    lambda1=myTHDMW->getTHDMW_lambda1();
    lambda2=myTHDMW->getTHDMW_lambda2();
    lambda3=myTHDMW->getTHDMW_lambda3();
    lambda4=myTHDMW->getTHDMW_lambda4();
    lambda5=myTHDMW->getTHDMW_lambda5();
    mSsq=myTHDMW->getTHDMW_mS2();
    mu1=myTHDMW->getTHDMW_mu1();
    mu2=myTHDMW->getTHDMW_mu2();
    mu3=myTHDMW->getTHDMW_mu3();
    mu4=myTHDMW->getTHDMW_mu4();
    mu5=myTHDMW->getTHDMW_mu5();
    mu6=myTHDMW->getTHDMW_mu6();
    nu1=myTHDMW->getTHDMW_nu1();
    nu2=myTHDMW->getTHDMW_nu2();
    nu3=myTHDMW->getTHDMW_nu3();
    nu4=myTHDMW->getTHDMW_nu4();
    nu5=myTHDMW->getTHDMW_nu5();
    omega1=myTHDMW->getTHDMW_omega1();
    omega2=myTHDMW->getTHDMW_omega2();
    omega3=myTHDMW->getTHDMW_omega3();
    omega4=myTHDMW->getTHDMW_omega4();
    kappa1=myTHDMW->getTHDMW_kappa1();
    kappa2=myTHDMW->getTHDMW_kappa2();
    kappa3=myTHDMW->getTHDMW_kappa3();
    etaU=myTHDMW->getTHDMW_etaU();
    etaD=myTHDMW->getTHDMW_etaD();
    rho_b=myTHDMW->getTHDMW_rho_b();
    S_b=myTHDMW->getTHDMW_S_b();
    RpepsTHDMW=myTHDMW->getRpepsTHDMW();


    setOtherParameters();
    computeHHlimits();
    runTHDMWparameters();
    computeUnitarity();
    computeSignalStrengthQuantities();
}
