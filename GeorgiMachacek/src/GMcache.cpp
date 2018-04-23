/* 
 * Copyright (C) 2017 HEPfit Collaboration
 * All rights reserved.
 *
 * For the licensing terms see doc/COPYING.
 */

#include "gslpp.h"
#include "GMcache.h"
#include <fstream>
#include <sstream>
#include <string>

GMcache::GMcache(const StandardModel& SM_i)
    :br_tt(19961, 2, 0.),
    br_bb(19961, 2, 0.),
    br_tautau(19961, 2, 0.),
    br_cc(19961, 2, 0.),
    br_mumu(19961, 2, 0.),
    br_ZZ(19961, 2, 0.),
    br_WW(19961, 2, 0.),
    GammaHtot_SM(19961, 2, 0.),
    log_cs_ggH_8(199, 2, 0.),
    log_cs_VBF_8(199, 2, 0.),
    log_cs_WH_8(199, 2, 0.),
    log_cs_ZH_8(199, 2, 0.),
    log_cs_ggH_13(199, 2, 0.),
    log_cs_VBF_13(199, 2, 0.),
    log_cs_WH_13(199, 2, 0.),
    log_cs_ZH_13(199, 2, 0.),
    log_cs_ttH_8(199, 2, 0.),
    log_cs_ttH_13(199, 2, 0.),
    log_cs_bbH_8(199, 2, 0.),
    log_cs_bbH_13(199, 2, 0.),
    log_cs_ggA_8(199, 2, 0.),
    log_cs_ttA_8(199, 2, 0.),
    log_cs_bbA_8(199, 2, 0.),
    log_cs_ggA_13(199, 2, 0.),
    log_cs_ttA_13(199, 2, 0.),
    log_cs_bbA_13(199, 2, 0.),
    log_cs_ggHp_8(744, 3, 0.),
    log_cs_ggHp_13(1104, 3, 0.),
    ATLAS8_gg_phi_tautau(92, 2, 0.),
    ATLAS8_bb_phi_tautau(92, 2, 0.),
    ATLAS8_gg_phi_tautau_e(92, 2, 0.),
    ATLAS8_bb_phi_tautau_e(92, 2, 0.),
    CMS13_pp_H_hh_bbbb(95,2,0.),
    CMS13_pp_H_hh_bbbb_e(95,2,0.),
    unitarityeigenvalues(17, 0.),
    myGM(static_cast<const GeorgiMachacek*> (&SM_i))
{
  read();
//    myRunnerGM=new RunnerGM(SM_i);
}

GMcache::~GMcache()
{
//  delete myRunnerGM;
}
//
///////////////////////////////////////////////////////////////////////////////////////////////////
//
int GMcache::CacheCheck(const gslpp::complex cache[][CacheSize], 
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

int GMcache::CacheCheckReal(const double cache[][CacheSize], 
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

void GMcache::CacheShift(gslpp::complex cache[][CacheSize], const int NumPar, 
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

void GMcache::CacheShiftReal(double cache[][CacheSize], const int NumPar,
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

void GMcache::read(){

    std::stringstream br1,br2,br3,br4,br5,br6,br7;
    std::stringstream dw1;
    std::stringstream cs1,cs2,cs3,cs4,cs5,cs6,cs7,cs8,cs9;
    std::stringstream cs11,cs12,cs13,cs14,cs15,cs16,cs17,cs18,cs19;
    std::stringstream cs20,cs21;
    std::stringstream ex6,ex7,ex53;
    std::stringstream ex6e,ex7e,ex53e;

    std::cout<<"reading tables"<<std::endl;

    std::stringstream path;
    path << getenv("HEPFITTABS") << "/THDM/tabs/";
    std::string tablepath=path.str();

    br1 << tablepath << "br1.dat";
    br_tt = readTable(br1.str(),19961,2);
    br2 << tablepath << "br2.dat";
    br_bb = readTable(br2.str(),19961,2);
    br3 << tablepath << "br3.dat";
    br_tautau = readTable(br3.str(),19961,2); 
    br4 << tablepath << "br4.dat";
    br_cc = readTable(br4.str(),19961,2);
    br5 << tablepath << "br5.dat";
    br_mumu = readTable(br5.str(),19961,2);
    br6 << tablepath << "br6.dat";
    br_ZZ = readTable(br6.str(),19961,2);
    br7 << tablepath << "br7.dat";
    br_WW = readTable(br7.str(),19961,2);
    dw1 << tablepath << "dw1.dat";
    GammaHtot_SM = readTable(dw1.str(),19961,2);
    cs1 << tablepath << "log_cs_ggH_8.dat";
    log_cs_ggH_8 = readTable(cs1.str(),199,2);
    cs11 << tablepath << "log_cs_ggH_13.dat";
    log_cs_ggH_13 = readTable(cs11.str(),199,2);
    cs2 << tablepath << "log_cs_VBF_8.dat";
    log_cs_VBF_8 = readTable(cs2.str(),199,2);
    cs12 << tablepath << "log_cs_VBF_13.dat";
    log_cs_VBF_13 = readTable(cs12.str(),199,2);
    cs3 << tablepath << "log_cs_WH_8.dat";
    log_cs_WH_8 = readTable(cs3.str(),199,2);
    cs13 << tablepath << "log_cs_WH_13.dat";
    log_cs_WH_13 = readTable(cs13.str(),199,2);
    cs4 << tablepath << "log_cs_ZH_8.dat";
    log_cs_ZH_8 = readTable(cs4.str(),199,2);
    cs14 << tablepath << "log_cs_ZH_13.dat";
    log_cs_ZH_13 = readTable(cs14.str(),199,2);
    cs5 << tablepath << "log_cs_ttH_8.dat";
    log_cs_ttH_8 = readTable(cs5.str(),199,2);
    cs15 << tablepath << "log_cs_ttH_13.dat";
    log_cs_ttH_13 = readTable(cs15.str(),199,2);
    cs6 << tablepath << "log_cs_bbH_8.dat";
    log_cs_bbH_8 = readTable(cs6.str(),199,2);
    cs16 << tablepath << "log_cs_bbH_13.dat";
    log_cs_bbH_13 = readTable(cs16.str(),199,2);
    cs7 << tablepath << "log_cs_ggA_8.dat";
    log_cs_ggA_8 = readTable(cs7.str(),199,2);
    cs17 << tablepath << "log_cs_ggA_13.dat";
    log_cs_ggA_13 = readTable(cs17.str(),199,2);
    cs8 << tablepath << "log_cs_ttA_8.dat";
    log_cs_ttA_8 = readTable(cs8.str(),199,2);
    cs18 << tablepath << "log_cs_ttA_13.dat";
    log_cs_ttA_13 = readTable(cs18.str(),199,2);
    cs9 << tablepath << "log_cs_bbA_8.dat";
    log_cs_bbA_8 = readTable(cs9.str(),199,2);
    cs19 << tablepath << "log_cs_bbA_13.dat";
    log_cs_bbA_13 = readTable(cs19.str(),199,2);
    cs20 << tablepath << "log_cs_ggHp_8.dat";
    log_cs_ggHp_8 = readTable(cs20.str(),744,3);
    cs21 << tablepath << "log_cs_ggHp_13.dat";
    log_cs_ggHp_13 = readTable(cs21.str(),1104,3);

    ex6 << tablepath << "14096064_a.dat";
    ATLAS8_gg_phi_tautau = readTable(ex6.str(),92,2);
    ex6e << tablepath << "14096064_a_e.dat";
    ATLAS8_gg_phi_tautau_e = readTable(ex6e.str(),92,2);
    ex7 << tablepath << "14096064_b.dat";
    ATLAS8_bb_phi_tautau = readTable(ex7.str(),92,2);
    ex7e << tablepath << "14096064_b_e.dat";
    ATLAS8_bb_phi_tautau_e = readTable(ex7e.str(),92,2);

    ex53 << tablepath << "CMS-PAS-HIG-17-009.dat";
    CMS13_pp_H_hh_bbbb = readTable(ex53.str(),95,2);
    ex53e << tablepath << "CMS-PAS-HIG-17-009_e.dat";
    CMS13_pp_H_hh_bbbb_e = readTable(ex53e.str(),95,2);

}    
    


double GMcache::ip_Br_HPtott(double mass){
    int NumPar = 1;
    double params[] = {mass};

    int i = CacheCheckReal(ip_Br_HPtott_cache, NumPar, params);
    if (i>=0) {
        return ( ip_Br_HPtott_cache[NumPar][i] );
    } else {
        double newResult = pow(10.0,interpolate(br_tt,mass));
        CacheShiftReal(ip_Br_HPtott_cache, NumPar, params, newResult);
        return newResult;
    }
}



double GMcache::ip_Br_HPtobb(double mass){
    int NumPar = 1;
    double params[] = {mass};

    int i = CacheCheckReal(ip_Br_HPtobb_cache, NumPar, params);
    if (i>=0) {
        return ( ip_Br_HPtobb_cache[NumPar][i] );
    } else {
        double newResult = pow(10.0,interpolate(br_bb,mass));
        CacheShiftReal(ip_Br_HPtobb_cache, NumPar, params, newResult);
        return newResult;
    }
}



double GMcache::ip_Br_HPtotautau(double mass){
    int NumPar = 1;
    double params[] = {mass};

    int i = CacheCheckReal(ip_Br_HPtotautau_cache, NumPar, params);
    if (i>=0) {
        return ( ip_Br_HPtotautau_cache[NumPar][i] );
    } else {
        double newResult = pow(10.0,interpolate(br_tautau,mass));
        CacheShiftReal(ip_Br_HPtotautau_cache, NumPar, params, newResult);
        return newResult;
    }
}



double GMcache::ip_Br_HPtocc(double mass){
    int NumPar = 1;
    double params[] = {mass};

    int i = CacheCheckReal(ip_Br_HPtocc_cache, NumPar, params);
    if (i>=0) {
        return ( ip_Br_HPtocc_cache[NumPar][i] );
    } else {
        double newResult = pow(10.0,interpolate(br_cc,mass));
        CacheShiftReal(ip_Br_HPtocc_cache, NumPar, params, newResult);
        return newResult;
    }
}



double GMcache::ip_Br_HPtomumu(double mass){
    int NumPar = 1;
    double params[] = {mass};

    int i = CacheCheckReal(ip_Br_HPtomumu_cache, NumPar, params);
    if (i>=0) {
        return ( ip_Br_HPtomumu_cache[NumPar][i] );
    } else {
        double newResult = pow(10.0,interpolate(br_mumu,mass));
        CacheShiftReal(ip_Br_HPtomumu_cache, NumPar, params, newResult);
        return newResult;
    }
}



double GMcache::ip_Br_HPtoZZ(double mass){
    int NumPar = 1;
    double params[] = {mass};

    int i = CacheCheckReal(ip_Br_HPtoZZ_cache, NumPar, params);
    if (i>=0) {
        return ( ip_Br_HPtoZZ_cache[NumPar][i] );
    } else {
        double newResult = pow(10.0,interpolate(br_ZZ,mass));
        CacheShiftReal(ip_Br_HPtoZZ_cache, NumPar, params, newResult);
        return newResult;
    }
}



double GMcache::ip_Br_HPtoWW(double mass){
    int NumPar = 1;
    double params[] = {mass};

    int i = CacheCheckReal(ip_Br_HPtoWW_cache, NumPar, params);
    if (i>=0) {
        return ( ip_Br_HPtoWW_cache[NumPar][i] );
    } else {
        double newResult = pow(10.0,interpolate(br_WW,mass));
        CacheShiftReal(ip_Br_HPtoWW_cache, NumPar, params, newResult);
        return newResult;
    }
}



double GMcache::ip_GammaHPtotSM(double mass){
    int NumPar = 1;
    double params[] = {mass};

    int i = CacheCheckReal(ip_GammaHPtotSM_cache, NumPar, params);
    if (i>=0) {
        return ( ip_GammaHPtotSM_cache[NumPar][i] );
    } else {
        double newResult = pow(10.0,interpolate(GammaHtot_SM,mass));
        CacheShiftReal(ip_GammaHPtotSM_cache, NumPar, params, newResult);
        return newResult;
    }
}



double GMcache::ip_cs_ggtoH_8(double mass){
    int NumPar = 1;
    double params[] = {mass};

    int i = CacheCheckReal(ip_cs_ggtoH_8_cache, NumPar, params);
    if (i>=0) {
        return ( ip_cs_ggtoH_8_cache[NumPar][i] );
    } else {
        double newResult = 0.0;
        if (mass>=20. && mass <=2000.) {
            newResult = pow(10.0,interpolate(log_cs_ggH_8,mass));
        }
        CacheShiftReal(ip_cs_ggtoH_8_cache, NumPar, params, newResult);
        return newResult;
    }
}



double GMcache::ip_cs_ggtoH_13(double mass){
    int NumPar = 1;
    double params[] = {mass};

    int i = CacheCheckReal(ip_cs_ggtoH_13_cache, NumPar, params);
    if (i>=0) {
        return ( ip_cs_ggtoH_13_cache[NumPar][i] );
    } else {
        double newResult = 0.0;
        if (mass>=20. && mass <=2000.) {
            newResult = pow(10.0,interpolate (log_cs_ggH_13,mass));
        }
        CacheShiftReal(ip_cs_ggtoH_13_cache, NumPar, params, newResult);
        return newResult;
    }
}



double GMcache::ip_cs_VBFtoH_8(double mass){
    int NumPar = 1;
    double params[] = {mass};

    int i = CacheCheckReal(ip_cs_VBFtoH_8_cache, NumPar, params);
    if (i>=0) {
        return ( ip_cs_VBFtoH_8_cache[NumPar][i] );
    } else {
        double newResult = 0.0;
        if (mass>=20. && mass <=2000.) {
            newResult = pow(10.0,interpolate (log_cs_VBF_8,mass));
        }
        CacheShiftReal(ip_cs_VBFtoH_8_cache, NumPar, params, newResult);
        return newResult;
    }
}



double GMcache::ip_cs_VBFtoH_13(double mass){
    int NumPar = 1;
    double params[] = {mass};

    int i = CacheCheckReal(ip_cs_VBFtoH_13_cache, NumPar, params);
    if (i>=0) {
        return ( ip_cs_VBFtoH_13_cache[NumPar][i] );
    } else {
        double newResult = 0.0;
        if (mass>=20. && mass <=2000.) {
            newResult = pow(10.0,interpolate (log_cs_VBF_13,mass));
        }
        CacheShiftReal(ip_cs_VBFtoH_13_cache, NumPar, params, newResult);
        return newResult;
    }
}



double GMcache::ip_cs_WtoWH_8(double mass){
    int NumPar = 1;
    double params[] = {mass};

    int i = CacheCheckReal(ip_cs_WtoWH_8_cache, NumPar, params);
    if (i>=0) {
        return ( ip_cs_WtoWH_8_cache[NumPar][i] );
    } else {
        double newResult = 0.0;
        if (mass>=20. && mass <=2000.) {
            newResult = pow(10.0,interpolate (log_cs_WH_8,mass));
        }
        CacheShiftReal(ip_cs_WtoWH_8_cache, NumPar, params, newResult);
        return newResult;
    }
}



double GMcache::ip_cs_WtoWH_13(double mass){
    int NumPar = 1;
    double params[] = {mass};

    int i = CacheCheckReal(ip_cs_WtoWH_13_cache, NumPar, params);
    if (i>=0) {
        return ( ip_cs_WtoWH_13_cache[NumPar][i] );
    } else {
        double newResult = 0.0;
        if (mass>=20. && mass <=2000.) {
            newResult = pow(10.0,interpolate (log_cs_WH_13,mass));
        }
        CacheShiftReal(ip_cs_WtoWH_13_cache, NumPar, params, newResult);
        return newResult;
    }
}



double GMcache::ip_cs_ZtoZH_8(double mass){
    int NumPar = 1;
    double params[] = {mass};

    int i = CacheCheckReal(ip_cs_ZtoZH_8_cache, NumPar, params);
    if (i>=0) {
        return ( ip_cs_ZtoZH_8_cache[NumPar][i] );
    } else {
        double newResult = 0.0;
        if (mass>=20. && mass <=2000.) {
            newResult = pow(10.0,interpolate (log_cs_ZH_8,mass));
        }
        CacheShiftReal(ip_cs_ZtoZH_8_cache, NumPar, params, newResult);
        return newResult;
    }
}



double GMcache::ip_cs_ZtoZH_13(double mass){
    int NumPar = 1;
    double params[] = {mass};

    int i = CacheCheckReal(ip_cs_ZtoZH_13_cache, NumPar, params);
    if (i>=0) {
        return ( ip_cs_ZtoZH_13_cache[NumPar][i] );
    } else {
        double newResult = 0.0;
        if (mass>=20. && mass <=2000.) {
            newResult = pow(10.0,interpolate (log_cs_ZH_13,mass));
        }
        CacheShiftReal(ip_cs_ZtoZH_13_cache, NumPar, params, newResult);
        return newResult;
    }
}



double GMcache::ip_cs_pptottH_8(double mass){
    int NumPar = 1;
    double params[] = {mass};

    int i = CacheCheckReal(ip_cs_pptottH_8_cache, NumPar, params);
    if (i>=0) {
        return ( ip_cs_pptottH_8_cache[NumPar][i] );
    } else {
        double newResult = 0.0;
        if (mass>=20. && mass <=2000.) {
            newResult = pow(10.0,interpolate (log_cs_ttH_8,mass));
        }
        CacheShiftReal(ip_cs_pptottH_8_cache, NumPar, params, newResult);
        return newResult;
    }
}



double GMcache::ip_cs_pptottH_13(double mass){
    int NumPar = 1;
    double params[] = {mass};

    int i = CacheCheckReal(ip_cs_pptottH_13_cache, NumPar, params);
    if (i>=0) {
        return ( ip_cs_pptottH_13_cache[NumPar][i] );
    } else {
        double newResult = 0.0;
        if (mass>=20. && mass <=2000.) {
            newResult = pow(10.0,interpolate (log_cs_ttH_13,mass));
        }
        CacheShiftReal(ip_cs_pptottH_13_cache, NumPar, params, newResult);
        return newResult;
    }
}



double GMcache::ip_cs_pptobbH_8(double mass){
    int NumPar = 1;
    double params[] = {mass};

    int i = CacheCheckReal(ip_cs_pptobbH_8_cache, NumPar, params);
    if (i>=0) {
        return ( ip_cs_pptobbH_8_cache[NumPar][i] );
    } else {
        double newResult = 0.0;
        if (mass>=20. && mass <=2000.) {
            newResult = pow(10.0,interpolate (log_cs_bbH_8,mass));
        }
        CacheShiftReal(ip_cs_pptobbH_8_cache, NumPar, params, newResult);
        return newResult;
    }
}



double GMcache::ip_cs_pptobbH_13(double mass){
    int NumPar = 1;
    double params[] = {mass};

    int i = CacheCheckReal(ip_cs_pptobbH_13_cache, NumPar, params);
    if (i>=0) {
        return ( ip_cs_pptobbH_13_cache[NumPar][i] );
    } else {
        double newResult = 0.0;
        if (mass>=20. && mass <=2000.) {
            newResult = pow(10.0,interpolate (log_cs_bbH_13,mass));
        }
        CacheShiftReal(ip_cs_pptobbH_13_cache, NumPar, params, newResult);
        return newResult;
    }
}



double GMcache::ip_cs_ggtoA_8(double mass){
    int NumPar = 1;
    double params[] = {mass};

    int i = CacheCheckReal(ip_cs_ggtoA_8_cache, NumPar, params);
    if (i>=0) {
        return ( ip_cs_ggtoA_8_cache[NumPar][i] );
    } else {
        double newResult = 0.0;
        if (mass>=20. && mass <=2000.) {
            newResult = pow(10.0,interpolate (log_cs_ggA_8,mass));
        }
        CacheShiftReal(ip_cs_ggtoA_8_cache, NumPar, params, newResult);
        return newResult;
    }
}



double GMcache::ip_cs_ggtoA_13(double mass){
    int NumPar = 1;
    double params[] = {mass};

    int i = CacheCheckReal(ip_cs_ggtoA_13_cache, NumPar, params);
    if (i>=0) {
        return ( ip_cs_ggtoA_13_cache[NumPar][i] );
    } else {
        double newResult = 0.0;
        if (mass>=20. && mass <=2000.) {
            newResult = pow(10.0,interpolate (log_cs_ggA_13,mass));
        }
        CacheShiftReal(ip_cs_ggtoA_13_cache, NumPar, params, newResult);
        return newResult;
    }
}



double GMcache::ip_cs_pptottA_8(double mass){
    int NumPar = 1;
    double params[] = {mass};

    int i = CacheCheckReal(ip_cs_pptottA_8_cache, NumPar, params);
    if (i>=0) {
        return ( ip_cs_pptottA_8_cache[NumPar][i] );
    } else {
        double newResult = 0.0;
        if (mass>=20. && mass <=2000.) {
            newResult = pow(10.0,interpolate (log_cs_ttA_8,mass));
        }
        CacheShiftReal(ip_cs_pptottA_8_cache, NumPar, params, newResult);
        return newResult;
    }
}



double GMcache::ip_cs_pptottA_13(double mass){
    int NumPar = 1;
    double params[] = {mass};

    int i = CacheCheckReal(ip_cs_pptottA_13_cache, NumPar, params);
    if (i>=0) {
        return ( ip_cs_pptottA_13_cache[NumPar][i] );
    } else {
        double newResult = 0.0;
        if (mass>=20. && mass <=2000.) {
            newResult = pow(10.0,interpolate (log_cs_ttA_13,mass));
        }
        CacheShiftReal(ip_cs_pptottA_13_cache, NumPar, params, newResult);
        return newResult;
    }
}



double GMcache::ip_cs_pptobbA_8(double mass){
    int NumPar = 1;
    double params[] = {mass};

    int i = CacheCheckReal(ip_cs_pptobbA_8_cache, NumPar, params);
    if (i>=0) {
        return ( ip_cs_pptobbA_8_cache[NumPar][i] );
    } else {
        double newResult = 0.0;
        if (mass>=20. && mass <=2000.) {
            newResult = pow(10.0,interpolate (log_cs_bbA_8,mass));
        }
        CacheShiftReal(ip_cs_pptobbA_8_cache, NumPar, params, newResult);
        return newResult;
    }
}



double GMcache::ip_cs_pptobbA_13(double mass){
    int NumPar = 1;
    double params[] = {mass};

    int i = CacheCheckReal(ip_cs_pptobbA_13_cache, NumPar, params);
    if (i>=0) {
        return ( ip_cs_pptobbA_13_cache[NumPar][i] );
    } else {
        double newResult = 0.0;
        if (mass>=20. && mass <=2000.) {
            newResult = pow(10.0,interpolate (log_cs_bbA_13,mass));
        }
        CacheShiftReal(ip_cs_pptobbA_13_cache, NumPar, params, newResult);
        return newResult;
    }
}



double GMcache::ip_cs_ggtoHp_8(double mHp, double logtb){
    int NumPar = 2;
    double params[] = {mHp, logtb};

    int i = CacheCheckReal(ip_cs_ggtoHp_8_cache, NumPar, params);
    if (i>=0) {
        return ( ip_cs_ggtoHp_8_cache[NumPar][i] );
    } else {
        double newResult = 0.0;
        if (mHp>=180. && mHp <=1400. && logtb>=-1. && logtb<=1.75) {
            newResult = pow(10.0,interpolate2D(log_cs_ggHp_8, logtb, mHp));
        }
        CacheShiftReal(ip_cs_ggtoHp_8_cache, NumPar, params, newResult);
        return newResult;
    }
}



double GMcache::ip_cs_ggtoHp_13(double mHp, double logtb){
    int NumPar = 2;
    double params[] = {mHp, logtb};

    int i = CacheCheckReal(ip_cs_ggtoHp_13_cache, NumPar, params);
    if (i>=0) {
        return ( ip_cs_ggtoHp_13_cache[NumPar][i] );
    } else {
        double newResult = 0.0;
        if (mHp>=180. && mHp <=2000. && logtb>=-1. && logtb<=1.75) {
            newResult = pow(10.0,interpolate2D(log_cs_ggHp_13, logtb, mHp));
        }
        CacheShiftReal(ip_cs_ggtoHp_13_cache, NumPar, params, newResult);
        return newResult;
    }
}




double GMcache::ip_ex_gg_phi_tautau_ATLAS8(double mass){
    int NumPar = 1;
    double params[] = {mass};

    int i = CacheCheckReal(ip_ex_gg_phi_tautau_ATLAS8_cache, NumPar, params);
    if (i>=0) {
        return ( ip_ex_gg_phi_tautau_ATLAS8_cache[NumPar][i] );
    } else {
        double newResult = interpolate(ATLAS8_gg_phi_tautau,mass);
        CacheShiftReal(ip_ex_gg_phi_tautau_ATLAS8_cache, NumPar, params, newResult);
        return newResult;
    }
}



double GMcache::ip_ex_gg_phi_tautau_ATLAS8_e(double mass){
    int NumPar = 1;
    double params[] = {mass};

    int i = CacheCheckReal(ip_ex_gg_phi_tautau_ATLAS8_cache_e, NumPar, params);
    if (i>=0) {
        return ( ip_ex_gg_phi_tautau_ATLAS8_cache_e[NumPar][i] );
    } else {
        double newResult = interpolate(ATLAS8_gg_phi_tautau_e,mass);
        CacheShiftReal(ip_ex_gg_phi_tautau_ATLAS8_cache_e, NumPar, params, newResult);
        return newResult;
    }
}



double GMcache::ip_ex_bb_phi_tautau_ATLAS8(double mass){
    int NumPar = 1;
    double params[] = {mass};

    int i = CacheCheckReal(ip_ex_bb_phi_tautau_ATLAS8_cache, NumPar, params);
    if (i>=0) {
        return ( ip_ex_bb_phi_tautau_ATLAS8_cache[NumPar][i] );
    } else {
        double newResult = interpolate(ATLAS8_bb_phi_tautau,mass);
        CacheShiftReal(ip_ex_bb_phi_tautau_ATLAS8_cache, NumPar, params, newResult);
        return newResult;
    }
}



double GMcache::ip_ex_bb_phi_tautau_ATLAS8_e(double mass){
    int NumPar = 1;
    double params[] = {mass};

    int i = CacheCheckReal(ip_ex_bb_phi_tautau_ATLAS8_cache_e, NumPar, params);
    if (i>=0) {
        return ( ip_ex_bb_phi_tautau_ATLAS8_cache_e[NumPar][i] );
    } else {
        double newResult = interpolate(ATLAS8_bb_phi_tautau_e,mass);
        CacheShiftReal(ip_ex_bb_phi_tautau_ATLAS8_cache_e, NumPar, params, newResult);
        return newResult;
    }
}



double GMcache::ip_ex_pp_H_hh_bbbb_CMS13(double mass){
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



double GMcache::ip_ex_pp_H_hh_bbbb_CMS13_e(double mass){
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



gslpp::matrix<double> GMcache::readTable(std::string filename, int rowN, int colN){

    std::ifstream INfile;
    std::string lineTab;
    INfile.open( filename.c_str() );
    if(INfile.fail()){
        std::cout<<"error: in GMcache, table doesn't exist!"<<std::endl;
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

double GMcache::interpolate(gslpp::matrix<double> arrayTab, double x){

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

//2D interpolation

double GMcache::interpolate2D(gslpp::matrix<double> arrayTab, double x, double y){

    int rowN=arrayTab.size_i();

    double xmin = arrayTab(0,0);
    double xmax = arrayTab(rowN-1,0);
    double ymin = arrayTab(0,1);
    double ymax = arrayTab(rowN-1,1);
    double intervalx = arrayTab(1,0)-arrayTab(0,0);
    int i=1;
    do i++;
    while(arrayTab(i,1)-arrayTab(i-1,1)==0&&i<30000);
    double intervaly = arrayTab(i,1)-arrayTab(i-1,1);
    int Nintervalsx = (x-xmin)/intervalx;
    int Nintervalsy = (y-ymin)/intervaly;
    if(x<xmin||x>xmax||y<ymin||y>ymax){
//        std::cout<<"warning: the parameter point lies outside the table"<<std::endl;
        return 0.;
    }
    else{
    double x1=arrayTab(i*Nintervalsy+Nintervalsx,0);
    double x2=arrayTab(i*Nintervalsy+Nintervalsx+1,0);
    double y1=arrayTab(i*Nintervalsy+Nintervalsx,1);
    double y2=arrayTab(i*(Nintervalsy+1)+Nintervalsx,1);
    return (arrayTab(i*Nintervalsy+Nintervalsx,2) * (x2-x) * (y2-y)
            +arrayTab(i*Nintervalsy+Nintervalsx+1,2) * (x-x1) * (y2-y)
            +arrayTab(i*(Nintervalsy+1)+Nintervalsx,2) * (x2-x) * (y-y1)
            +arrayTab(i*(Nintervalsy+1)+Nintervalsx+1,2) * (x-x1) * (y-y1))
           /((x2-x1)*(y2-y1));
    }
}


gslpp::complex GMcache::I_h_U(const double mHl2, const double Mu, const double Mc, const double Mt) const {
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

gslpp::complex GMcache::I_HH_U(const double mHh2, const double Mc, const double Mt) const {
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

gslpp::complex GMcache::I_A_U(const double mA2, const double Mc, const double Mt) const {
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

gslpp::complex GMcache::I_h_D(const double mHl2, const double Md, const double Ms, const double Mb) const {
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

gslpp::complex GMcache::I_HH_D(const double mHh2, const double Ms, const double Mb) const {
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

gslpp::complex GMcache::I_A_D(const double mA2, const double Ms, const double Mb) const {
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

gslpp::complex GMcache::I_h_L(const double mHl2, const double Me, const double Mmu, const double Mtau) const {
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

gslpp::complex GMcache::I_HH_L(const double mHh2, const double Mmu, const double Mtau) const {
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

gslpp::complex GMcache::I_A_L(const double mA2, const double Mmu, const double Mtau) const {
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

gslpp::complex GMcache::I_H_W(const double mH, const double MW) const {
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

gslpp::complex GMcache::I_H_Hp(const double mHp2, const double mH) const {
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

gslpp::complex GMcache::A_h_U(const double mHl2, const double cW2, const double Mu, const double Mc, const double Mt, const double MZ) const {
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

gslpp::complex GMcache::A_HH_U(const double mHh2, const double cW2, const double Mc, const double Mt, const double MZ) const {
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

gslpp::complex GMcache::A_A_U(const double mA2, const double cW2, const double Mc, const double Mt, const double MZ) const {
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

gslpp::complex GMcache::A_h_D(const double mHl2, const double cW2, const double Md, const double Ms, const double Mb, const double MZ) const {
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

gslpp::complex GMcache::A_HH_D(const double mHh2, const double cW2, const double Ms, const double Mb, const double MZ) const {
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

gslpp::complex GMcache::A_A_D(const double mA2, const double cW2, const double Ms, const double Mb, const double MZ) const {
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

gslpp::complex GMcache::A_h_L(const double mHl2, const double cW2, const double Me, const double Mmu, const double Mtau, const double MZ) const {
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

gslpp::complex GMcache::A_HH_L(const double mHh2, const double cW2, const double Mmu, const double Mtau, const double MZ) const {
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

gslpp::complex GMcache::A_A_L(const double mA2, const double cW2, const double Mmu, const double Mtau, const double MZ) const {
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

gslpp::complex GMcache::A_H_W(const double mH, const double cW2, const double MW, const double MZ) const {
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

gslpp::complex GMcache::A_H_Hp(const double mHp2, const double mH, const double cW2, const double MZ) const {
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

gslpp::complex GMcache::f_func(const double x) const{
    if(x<1) {
    gslpp::complex z = -gslpp::complex::i()*M_PI;
    return -pow(log((1.0+sqrt(1.0-x))/(1.0-sqrt(1.0-x)))+z,2)/4.0;
    }
    else {
        return pow(asin(sqrt(1.0/x)),2);
    }
}

gslpp::complex GMcache::g_func(const double x) const{
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

gslpp::complex GMcache::Int1(const double tau, const double lambda) const{
    return tau*lambda/(tau-lambda)/2.0+tau*tau*lambda*lambda/((tau-lambda)
           *(tau-lambda))/2.0*(f_func(tau)-f_func(lambda))+tau*tau*lambda/((tau-lambda)
           *(tau-lambda))*(g_func(tau)-g_func(lambda));
}

gslpp::complex GMcache::Int2(const double tau, const double lambda) const{
    return -tau*lambda/(tau-lambda)/2.0*(f_func(tau)-f_func(lambda));
}

void GMcache::computeSignalStrengthQuantities()
{
    double Mt = myGM->getQuarks(QCD::TOP).getMass();
    double Mb = myGM->getQuarks(QCD::BOTTOM).getMass();
    double Mc = myGM->getQuarks(QCD::CHARM).getMass();
    double Ms = myGM->getQuarks(QCD::STRANGE).getMass();
    double Mu = myGM->getQuarks(QCD::UP).getMass();
    double Md = myGM->getQuarks(QCD::DOWN).getMass();
    double Mtau = myGM->getLeptons(StandardModel::TAU).getMass();
    double Mmu = myGM->getLeptons(StandardModel::MU).getMass();
    double Me = myGM->getLeptons(StandardModel::ELECTRON).getMass();

    double BrSM_htobb = myGM->computeBrHtobb();
    double BrSM_htotautau = myGM->computeBrHtotautau();
    double BrSM_htogaga = myGM->computeBrHtogaga();
    double BrSM_htoWW = myGM->computeBrHtoWW();
    double BrSM_htoZZ = myGM->computeBrHtoZZ();
    double BrSM_htogg = myGM->computeBrHtogg();
    double BrSM_htoZga = myGM->computeBrHtoZga();
    double BrSM_htocc = myGM->computeBrHtocc();

//    depending on which Higgs is heavier

    double gh_H3=0., gh_H5=0.;

    if(mH1sq>=mHl2)
    {
        rh_ff = cosa*cosa/(sinb*sinb);
        rh_VV = cosa*cosa*sinb*sinb - sqrt(32.0/3.0)*cosa*sina*cosb*sinb + 8.0/3.0*sina*sina*cosb*cosb;
        gh_H3 = ghH3pH3m;
        gh_H5 = ghH5pH5m;
    }
    else
    {
        rh_ff = sina*sina/(sinb*sinb);
        rh_VV = sina*sina*sinb*sinb + sqrt(32.0/3.0)*sina*cosa*cosb*sinb + 8.0/3.0*cosa*cosa*cosb*cosb;
        gh_H3 = gHH3pH3m;
        gh_H5 = gHH5pH5m;
    }

    //rh_gaga formula = abs(I_h_F+I_h_W+I_h_Hp)^2 / abs(I_hSM_F+I_hSM_W)^2

    gslpp::complex I_hSM_F = I_h_U(mHl2,Mu,Mc,Mt)+I_h_D(mHl2,Md,Ms,Mb)+I_h_L(mHl2,Me,Mmu,Mtau);
    gslpp::complex I_hSM_W = I_H_W(mHl2,MW);
    gslpp::complex I_h_F = sqrt(rh_ff)*I_hSM_F;
    gslpp::complex I_h_W = sqrt(rh_VV)*I_hSM_W;
    gslpp::complex I_h_Hp = 0.5*vev*( gh_H3*I_H_Hp(mAsq,mHl2)/mAsq + 5.0*gh_H5*I_H_Hp(mH5sq,mHl2)/mH5sq );

    double ABSgagaGM=(I_h_F+I_h_W+I_h_Hp).abs2();
    double ABSgagaSM=(I_hSM_F+I_hSM_W).abs2();
    rh_gaga=ABSgagaGM/ABSgagaSM;

    //rh_Zga formula = abs(A_h_F+A_h_W+A_h_Hp)^2 / abs(A_hSM_F+A_hSM_W)^2

    gslpp::complex A_hSM_F = A_h_U(mHl2,cW2,Mu,Mc,Mt,MZ)+A_h_D(mHl2,cW2,Md,Ms,Mb,MZ)+A_h_L(mHl2,cW2,Me,Mmu,Mtau,MZ);
    gslpp::complex A_hSM_W = A_H_W(mHl2,cW2,MW,MZ);
    gslpp::complex A_h_F = sqrt(rh_ff)*A_hSM_F;
    gslpp::complex A_h_W = sqrt(rh_VV)*A_hSM_W;
    gslpp::complex A_h_Hp = 0.5*vev*( gh_H3*A_H_Hp(mAsq,mHl2,cW2,MZ)/mAsq + 5.0*gh_H5*A_H_Hp(mH5sq,mHl2,cW2,MZ)/mH5sq );

    double ABSZgaGM=(A_h_F+A_h_W+A_h_Hp).abs2();
    double ABSZgaSM=(A_hSM_F+A_hSM_W).abs2();
    rh_Zga=ABSZgaGM/ABSZgaSM;

    rh_gg=rh_ff;

    sumModBRs = rh_ff*(BrSM_htobb+BrSM_htotautau+BrSM_htocc) + rh_VV*(BrSM_htoWW+BrSM_htoZZ) 
                + rh_gaga*BrSM_htogaga + rh_gg*BrSM_htogg + rh_Zga*BrSM_htoZga;

    Gamma_h = sumModBRs*myGM->computeGammaHTotal();

    GM_BR_h_bb = rh_ff*BrSM_htobb/sumModBRs;
//    GM_BR_h_gaga = rh_gaga*BrSM_htogaga/sumModBRs;
//    GM_BR_h_tautau = rh_ff*BrSM_htotautau/sumModBRs;
//    GM_BR_h_WW = rh_VV*BrSM_htoWW/sumModBRs;
//    GM_BR_h_ZZ = rh_VV*BrSM_htoZZ/sumModBRs;
}

int GMcache::HSTheta (const double x) const{
    if(x<0) return 0.0;
    else return 1.0;
}

double GMcache::KaellenFunction(const double a2, const double b2, const double c2) const{
    int NumPar = 3;
    double params[] = {a2, b2, c2};

    int i = CacheCheckReal(KaellenFunction_cache, NumPar, params);
    if (i>=0) {
        return ( KaellenFunction_cache[NumPar][i] );
    }
    else {
        double newResult = 0.0;
        double x = (a2-b2-c2)*(a2-b2-c2)-4.0*b2*c2;
        if(x>0) newResult = sqrt(std::fabs(x/a2))/2.0;        
        CacheShiftReal(KaellenFunction_cache, NumPar, params, newResult);
        return newResult;
    }
}

void GMcache::computeOtherHiggsProperties()
{
    double Als=myGM->getAlsMz();
    double sW2=1.0-cW2;
    double mHh=sqrt(mH1sq);
    double mA=sqrt(mAsq);
    double mH5=sqrt(mH5sq);
    double Mt = myGM->getQuarks(QCD::TOP).getMass();
    double Mb = myGM->getQuarks(QCD::BOTTOM).getMass();
    double Mc = myGM->getQuarks(QCD::CHARM).getMass();
    double Ms = myGM->getQuarks(QCD::STRANGE).getMass();
//    double Mu = myGM->getQuarks(QCD::UP).getMass();
//    double Md = myGM->getQuarks(QCD::DOWN).getMass();
    double Mtau = myGM->getLeptons(StandardModel::TAU).getMass();
    double Mmu = myGM->getLeptons(StandardModel::MU).getMass();
//    double Me = myGM->getLeptons(StandardModel::ELECTRON).getMass();

//    /* rX_ii is the squared ratio between the GM vertex coupling of the neutral heavy Higgs X to
//     * the particle i and the corresponding coupling of the SM Higgs boson.*/
    double rHH_ff=0.0;
    double rA_ff=0.0;
    double rH5_ff=0.0;
    double rHH_VV=0.0;
    double rA_VV=0.0;
    double rH5_VV=0.0;

    double gHH_H3=0., gHH_H5=0.;

    if(mH1sq>=mHl2)
    {
        rHH_ff = sina*sina/(sinb*sinb);
        rHH_VV = sina*sina*sinb*sinb + sqrt(32.0/3.0)*sina*cosa*cosb*sinb + 8.0/3.0*cosa*cosa*cosb*cosb;
        gHH_H3 = gHH3pH3m;
        gHH_H5 = gHH5pH5m;
    }
    else
    {
        rHH_ff = cosa*cosa/(sinb*sinb);
        rHH_VV = cosa*cosa*sinb*sinb - sqrt(32.0/3.0)*cosa*sina*cosb*sinb + 8.0/3.0*sina*sina*cosb*cosb;
        gHH_H3 = ghH3pH3m;
        gHH_H5 = ghH5pH5m;
    }

    gslpp::complex I_HH_F=0.0;//It depends on the modelType
    gslpp::complex I_HH_Ux=I_HH_U(mH1sq,Mc,Mt);
    gslpp::complex I_HH_Dx=I_HH_D(mH1sq,Ms,Mb);
    gslpp::complex I_HH_Lx=I_HH_L(mH1sq,Mmu,Mtau);
    gslpp::complex I_HH_W=sqrt(rHH_VV)*I_H_W(mHh,MW);
    gslpp::complex I_HH_Hp = 0.5*vev*( gHH_H3*I_H_Hp(mAsq,mH1sq)/mAsq + 5.0*gHH_H5*I_H_Hp(mH5sq,mH1sq)/mH5sq );

    //    /*A_HH_F, A_HH_W and A_HH_Hp are needed for Gamma_HZga;
//     * their expressions can be found in "The Higgs Hunter's Guide", Appendix C, C.12*/
    gslpp::complex A_HH_F = 0.0;//It depends on the modelType
    gslpp::complex A_HH_Ux = A_HH_U(mH1sq,cW2,Mc,Mt,MZ);
    gslpp::complex A_HH_Dx = A_HH_D(mH1sq,cW2,Ms,Mb,MZ);
    gslpp::complex A_HH_Lx = A_HH_L(mH1sq,cW2,Mmu,Mtau,MZ);
//    /*A_HH_W expression can be found in "The Higgs Hunter's Guide", Appendix C, C.13*/
    gslpp::complex A_HH_W = sqrt(rHH_VV)*A_H_W(mHh,cW2,MW,MZ);
//    /*A_HH_Hp expression can be found in "The Higgs Hunter's Guide", Appendix C, C.14*/
    gslpp::complex A_HH_Hp = 0.5*vev*( gHH_H3*A_H_Hp(mAsq,mH1sq,cW2,MZ)/mAsq + 5.0*gHH_H5*A_H_Hp(mH5sq,mH1sq,cW2,MZ)/mH5sq );

    I_HH_F=sqrt(rHH_ff)*(I_HH_Ux+I_HH_Dx+I_HH_Lx);
    A_HH_F=sqrt(rHH_ff)*(A_HH_Ux+A_HH_Dx+A_HH_Lx)/sqrt(sW2*cW2);
    /*Gamma_Hgaga expression can be found in arXiv:1412.7387, eq. (65)*/
    double Gamma_Hgaga=Ale*Ale*mHh*mHh*mHh/(256.0*M_PI*M_PI*M_PI*vev*vev)
                *(I_HH_F+I_HH_W+I_HH_Hp).abs2();
//
//    /*Gamma_HZga expression can be found in arXiv:0902.4665v3, Appendix A, A.9*/
    double Gamma_HZga=HSTheta(mHh-MZ)*GF*Ale*Ale*mHh*mHh*mHh/(sqrt(2.0)*64.0*M_PI*M_PI*M_PI)
               *(1.0-MZ*MZ/(mH1sq))*(1.0-MZ*MZ/(mH1sq))*(1.0-MZ*MZ/(mH1sq))
               *(A_HH_F+A_HH_W+A_HH_Hp).abs2();

    rHH_gg=rHH_ff;
    /*Gamma_Hgg expression can be found in arXiv:0902.4665v3, Appendix A, A.10 or in the Higgs Hunter's Guide (2.30); relative coupling see above*/
    double Gamma_Hgg=rHH_gg*GF*Als*Als*mHh*mHh*mHh/(sqrt(2.0)*16.0*M_PI*M_PI*M_PI)
                     *(9.0/4.0)*(I_HH_Ux/4.0+I_HH_Dx).abs2();

    double Gamma_Agaga=0.0;
    double Gamma_H5gaga=0.0;
    double Gamma_Agg=0.0;
    double Gamma_H5gg=0.0;
    double Gamma_AZga=0.0;
    double Gamma_H5Zga=0.0;

    SigmaggF_H8=ip_cs_ggtoH_8(mHh)*rHH_ff;
    SigmabbF_H8=ip_cs_pptobbH_8(mHh)*rHH_ff;
//    SigmaVBF_H8=ip_cs_VBFtoH_8(mHh)*rHH_VV;
//    double SigmattF_H8=ip_cs_pptottH_8(mHh)*rHH_QuQu;
//    double SigmaVH_H8=(ip_cs_WtoWH_8(mHh)+ip_cs_ZtoZH_8(mHh))*rHH_VV;
//    SigmaTotSM_H8 = 1.0e-15;
//    if (mHh>=20. && mHh <=2000.) {
//            SigmaTotSM_H8=ip_cs_ggtoH_8(mHh)+ip_cs_VBFtoH_8(mHh)+ip_cs_WtoWH_8(mHh)+ip_cs_ZtoZH_8(mHh)+ip_cs_pptottH_8(mHh)+ip_cs_pptobbH_8(mHh);
//    }
//    SigmaSumH8 = SigmaggF_H8 + SigmaVBF_H8 + SigmaVH_H8 + SigmattF_H8 + SigmabbF_H8;
//
    SigmaggF_H13=ip_cs_ggtoH_13(mHh)*rHH_ff;
    SigmabbF_H13=ip_cs_pptobbH_13(mHh)*rHH_ff;
    SigmaVBF_H13=ip_cs_VBFtoH_13(mHh)*rHH_VV;
    SigmattF_H13=ip_cs_pptottH_13(mHh)*rHH_ff;
    SigmaVH_H13=(ip_cs_WtoWH_13(mHh)+ip_cs_ZtoZH_13(mHh))*rHH_VV;
////    double SigmaTotSM_H13 = 1.0e-15;
////    if (mHh>=20. && mHh <=2000.) {
////            SigmaTotSM_H13=ip_cs_ggtoH_13(mHh)+ip_cs_VBFtoH_13(mHh)+ip_cs_WtoWH_13(mHh)+ip_cs_ZtoZH_13(mHh)+ip_cs_pptottH_13(mHh)+ip_cs_pptobbH_13(mHh);
////    }
    SigmaSumH13 = SigmaggF_H13 + SigmaVBF_H13 + SigmaVH_H13 + SigmattF_H13 + SigmabbF_H13;

    double BrSM_Htott=ip_Br_HPtott(mHh);
    double BrSM_Atott=ip_Br_HPtott(mA);
    double BrSM_H5tott=ip_Br_HPtott(mH5);

    double BrSM_Htocc=ip_Br_HPtocc(mHh);
    double BrSM_Atocc=ip_Br_HPtocc(mA);
    double BrSM_H5tocc=ip_Br_HPtocc(mH5);

    double BrSM_Htobb=ip_Br_HPtobb(mHh);
    double BrSM_Atobb=ip_Br_HPtobb(mA);
    double BrSM_H5tobb=ip_Br_HPtobb(mH5);

    double BrSM_Htotautau=ip_Br_HPtotautau(mHh);
    double BrSM_Atotautau=ip_Br_HPtotautau(mA);
    double BrSM_H5totautau=ip_Br_HPtotautau(mH5);

    double BrSM_Htomumu=ip_Br_HPtomumu(mHh);
    double BrSM_Atomumu=ip_Br_HPtomumu(mA);
    double BrSM_H5tomumu=ip_Br_HPtomumu(mH5);

    double BrSM_HtoWW =ip_Br_HPtoWW(mHh);
    double BrSM_AtoWW =ip_Br_HPtoWW(mA);
    double BrSM_H5toWW =ip_Br_HPtoWW(mH5);

    double BrSM_HtoZZ =ip_Br_HPtoZZ(mHh);
    double BrSM_AtoZZ =ip_Br_HPtoZZ(mA);
    double BrSM_H5toZZ =ip_Br_HPtoZZ(mH5);

    double GammaHtotSM=ip_GammaHPtotSM(mHh);
    double GammaAtotSM=ip_GammaHPtotSM(mA);
    double GammaH5totSM=ip_GammaHPtotSM(mH5);

    
    
//    ghHH = cosa*sinb*(4.0-6.0*cosa*cosa-2.0*sqrt(6.0)*sina*cosa*tanb)/vev * M1sq
//    gHHH = (2.0*cosa*cosa*sinb*(-3.0*sina+sqrt(6.0)*cosa*tanb))/vev * M1sq
//    ghH3H3 = -(2.0*sqrt(2.0/3.0)*sina)/(vev*cosb) * M1sq
//    gHH3H3 = (2.0*sqrt(2.0/3.0)*cosa)/(vev*cosb) * M1sq
//    ghH3pH3m = ghH3H3;
//    ghH5H5 = 2.0*sinb*(2.0*cosa+sqrt(6.0)*sina*tanb)/vev * (M1sq-2.0*mAsq)
//    gHH5H5 = 2.0*sinb*(2.0*sina-sqrt(6.0)*cosa*tanb)/vev * (M1sq-2.0*mAsq)
//    ghH5pH5m = ghH5H5;
//    gHH5pH5m = gHH5H5;
//    ghH5ppH5mm = ghH5H5;
//    gHH5ppH5mm = gHH5H5;
//    gH3H3H5 = 2.0*((mH5sq-2.0*mAsq)*cosb-(2.0*M1sq-3.0*mAsq+mH5sq)/cosb)/(sqrt(3.0)*vev);
//    gH5H5H5 = 2.0*((3.0*mH5sq-4.0*M2sq)/cosb+(6.0*M1sq-9.0*mAsq)*sinb*tanb)/(sqrt(3.0)*vev);
//    gH3H3pH5m = gslpp::complex::i()*(4.0*M1sq-4.0*mAsq+mH5sq+(2.0*mAsq-mH5sq)*(cosb*cosb-sinb*sinb))/(2.0*vev*cosb);
//    gH5H3pH3m = (4.0*M1sq-4.0*mAsq+mH5sq+(2.0*mAsq-mH5sq)*(cosb*cosb-sinb*sinb))/(2.0*sqrt(3.0)*vev*cosb);
//    gH5H5pH5m = 0.5*gH5H5H5;
//    gH5H5ppH5mm = -gH5H5H5;
//    gH5ppH3mH3m = sqrt(6.0)*gH5H3pH3m;
//    gH5ppH5mH5m = gH5H5pH5m;

    
    //Following partial widths stem from the GMcalc manual 1412.7387
    //lambda^{1/2}(x,y) is equivalent to 2*KaellenFunction(1,x,y)

    // H3 -> f f' decays
    
    rA_ff = 1.0/(tanb*tanb);

    // H1 -> V H2 decays

    double GammaHAZ       = HSTheta(mHh-sqrt(mAsq)-MZ)
                            *KaellenFunction(mH1sq,MZ*MZ,mAsq)*KaellenFunction(mH1sq,MZ*MZ,mAsq)*KaellenFunction(mH1sq,MZ*MZ,mAsq)
                            /(2.0*M_PI*MZ*MZ);

    double GammaHHpW      = HSTheta(mHh-sqrt(mAsq)-MW)
                            *KaellenFunction(mH1sq,MW*MW,mAsq)*KaellenFunction(mH1sq,MW*MW,mAsq)*KaellenFunction(mH1sq,MW*MW,mAsq)
                            /(2.0*M_PI*MW*MW);

    double GammaAhZ=0.;
    double GammaAHZ=0.;
    double GammaAHpW=0.;
    double GammaAH5Z=0.;
    double GammaAH5pW=0.;
    double GammaHphW=0.;
    double GammaHpHW=0.;
    double GammaHpZW=0.;
    double GammaHpH5pZ=0.;
    double GammaHpH5ppW   = 0.0;
    double GammaH5AZ=0.;
    double GammaH5HpW=0.;
    double GammaH5pZW=0.;
    double GammaH5pAW=0.;
    double GammaH5pHpZ=0.;
    double GammaH5ppWW=0.;

    // H1 -> H2 H3 decays

    double GammaHhh       = HSTheta(mHh-2.0*sqrt(mHl2))*fabs(ghhH)*fabs(ghhH)
                            *KaellenFunction(1.0,mHl2/mH1sq,mHl2/mH1sq)/(16.0*mHh*M_PI);

    double GammaHHpHm     = HSTheta(mHh-2.0*mA)*fabs(gHH3pH3m)*fabs(gHH3pH3m)
                            *KaellenFunction(1.0,mAsq/mH1sq,mAsq/mH1sq)/(8.0*mHh*M_PI);

    double GammaHAA       = HSTheta(mHh-2.0*mA)*fabs(gHH3H3)*fabs(gHH3H3)
                            *KaellenFunction(1.0,mAsq/mH1sq,mAsq/mH1sq)/(16.0*mHh*M_PI);
    
    double GammaHH5H5     = HSTheta(mHh-2.0*mH5)*fabs(gHH5H5)*fabs(gHH5H5)
                            *KaellenFunction(1.0,mH5sq/mH1sq,mH5sq/mH1sq)/(16.0*mHh*M_PI);

    double GammaHHp5H5m   = HSTheta(mHh-2.0*mH5)*fabs(gHH5pH5m)*fabs(gHH5pH5m)
                            *KaellenFunction(1.0,mH5sq/mH1sq,mH5sq/mH1sq)/(8.0*mHh*M_PI);

    double GammaHH5ppH5mm = HSTheta(mHh-2.0*mH5)*fabs(gHH5ppH5mm)*fabs(gHH5ppH5mm)
                            *KaellenFunction(1.0,mH5sq/mH1sq,mH5sq/mH1sq)/(8.0*mHh*M_PI);

    double GammaH5hh      = 0.;
    double GammaH5HH      = 0.;
    double GammaH5HpHm    = 0.;
    double GammaH5AA      = 0.;
    double GammaH5ppHpHp  = 0.;

    GammaH1tot=  (rHH_ff*(BrSM_Htott+BrSM_Htocc+BrSM_Htobb+BrSM_Htotautau+BrSM_Htomumu)
                  +rHH_VV*(BrSM_HtoWW+BrSM_HtoZZ))*GammaHtotSM
                  +Gamma_Hgg+Gamma_Hgaga+Gamma_HZga
                  +GammaHAZ+GammaHHpW
                  +GammaHhh+GammaHHpHm+GammaHAA
                  +GammaHH5H5+GammaHHp5H5m+GammaHH5ppH5mm;

    GammaH3tot=  (rA_ff*(BrSM_Atott+BrSM_Atocc+BrSM_Atobb+BrSM_Atotautau+BrSM_Atomumu)
                  +rA_VV*(BrSM_AtoWW+BrSM_AtoZZ))*GammaAtotSM
                  +Gamma_Agg+Gamma_Agaga+Gamma_AZga
                  +GammaAhZ+GammaAHZ+GammaAHpW
                  +GammaAH5Z+GammaAH5pW;

    GammaH3ptot=  0.0 /*fermionic and vector boson terms go here*/
                  +GammaHphW+GammaHpHW+GammaHpZW
                  +GammaHpH5pZ+GammaHpH5ppW;

    GammaH5tot=  (rH5_ff*(BrSM_H5tott+BrSM_H5tocc+BrSM_H5tobb+BrSM_H5totautau+BrSM_H5tomumu)
                  +rH5_VV*(BrSM_H5toWW+BrSM_H5toZZ))*GammaH5totSM
                  +Gamma_H5gg+Gamma_H5gaga+Gamma_H5Zga
                  +GammaH5AZ+GammaH5HpW
                  +GammaH5hh+GammaH5HH+GammaH5HpHm+GammaH5AA;

    GammaH5ptot=  0.0 /*fermionic and vector boson terms go here*/
                  +GammaH5pZW+GammaH5pAW+GammaH5pHpZ;

    GammaH5pptot= 0.0 /*fermionic and vector boson terms go here*/
                  +GammaH5ppWW+GammaH5ppHpHp;

//    Br_Htott=BrSM_Htott*rHH_QuQu*GammaHtotSM/GammaHtot;
//    Br_Htobb=BrSM_Htobb*rHH_QdQd*GammaHtotSM/GammaHtot;
    Br_Htotautau=BrSM_Htotautau*rHH_ff*GammaHtotSM/GammaH1tot;
//    Br_HtoWW=BrSM_HtoWW*rHH_VV*GammaHtotSM/GammaHtot;
//    Br_HtoZZ=BrSM_HtoZZ*rHH_VV*GammaHtotSM/GammaHtot;
//    Br_Htogaga=Gamma_Hgaga/GammaHtot;
//    Br_HtoZga=Gamma_HZga/GammaHtot;
    Br_Htohh=GammaHhh/GammaH1tot;
//    Br_HtoAA=GammaHAA/GammaHtot;
//    Br_HtoHpHm=GammaHHpHm/GammaHtot;
//    Br_HtoAZ=GammaHAZ/GammaHtot;
//    Br_HtoHpW=GammaHHpW/GammaHtot;

}    

void GMcache::computeDirectSearchQuantities()
{
    computeOtherHiggsProperties();
    double mHh=sqrt(mH1sq);

    ggF_H_tautau_TH8=SigmaggF_H8*Br_Htotautau;
    bbF_H_tautau_TH8=SigmabbF_H8*Br_Htotautau;
    pp_H_hh_bbbb_TH13=SigmaSumH13*Br_Htohh*GM_BR_h_bb*GM_BR_h_bb;

    THoEX_ggF_H_tautau_ATLAS8=0.0;
    R_ggF_H_tautau_ATLAS8=0.0;
    THoEX_pp_H_hh_bbbb_CMS13=0.0;
    R_pp_H_hh_bbbb_CMS13=0.0;

    //95% to 1 sigma conversion factor, roughly sqrt(3.84)
    double nftos=1.95996398454;

    if(mHh>=65.0 && mHh<90.0)
    {
    }
    else if(mHh>=90.0 && mHh<100.0)
    {
        THoEX_ggF_H_tautau_ATLAS8=ggF_H_tautau_TH8/ip_ex_gg_phi_tautau_ATLAS8(mHh);
        R_ggF_H_tautau_ATLAS8=(1+(ggF_H_tautau_TH8-ip_ex_gg_phi_tautau_ATLAS8(mHh))/ip_ex_gg_phi_tautau_ATLAS8_e(mHh) ) * nftos;
        THoEX_bbF_H_tautau_ATLAS8=bbF_H_tautau_TH8/ip_ex_bb_phi_tautau_ATLAS8(mHh);
        R_bbF_H_tautau_ATLAS8=(1+(bbF_H_tautau_TH8-ip_ex_bb_phi_tautau_ATLAS8(mHh))/ip_ex_bb_phi_tautau_ATLAS8_e(mHh) ) * nftos;
    }
    else if(mHh>=100.0 && mHh<130.0)
    {
        THoEX_ggF_H_tautau_ATLAS8=ggF_H_tautau_TH8/ip_ex_gg_phi_tautau_ATLAS8(mHh);
        R_ggF_H_tautau_ATLAS8=(1+(ggF_H_tautau_TH8-ip_ex_gg_phi_tautau_ATLAS8(mHh))/ip_ex_gg_phi_tautau_ATLAS8_e(mHh) ) * nftos;
        THoEX_bbF_H_tautau_ATLAS8=bbF_H_tautau_TH8/ip_ex_bb_phi_tautau_ATLAS8(mHh);
        R_bbF_H_tautau_ATLAS8=(1+(bbF_H_tautau_TH8-ip_ex_bb_phi_tautau_ATLAS8(mHh))/ip_ex_bb_phi_tautau_ATLAS8_e(mHh) ) * nftos;
    }
    else if(mHh>=130.0 && mHh<140.0)
    {
        THoEX_ggF_H_tautau_ATLAS8=ggF_H_tautau_TH8/ip_ex_gg_phi_tautau_ATLAS8(mHh);
        R_ggF_H_tautau_ATLAS8=(1+(ggF_H_tautau_TH8-ip_ex_gg_phi_tautau_ATLAS8(mHh))/ip_ex_gg_phi_tautau_ATLAS8_e(mHh) ) * nftos;
        THoEX_bbF_H_tautau_ATLAS8=bbF_H_tautau_TH8/ip_ex_bb_phi_tautau_ATLAS8(mHh);
        R_bbF_H_tautau_ATLAS8=(1+(bbF_H_tautau_TH8-ip_ex_bb_phi_tautau_ATLAS8(mHh))/ip_ex_bb_phi_tautau_ATLAS8_e(mHh) ) * nftos;
    }
    else if(mHh>=140.0 && mHh<145.0)
    {
        THoEX_ggF_H_tautau_ATLAS8=ggF_H_tautau_TH8/ip_ex_gg_phi_tautau_ATLAS8(mHh);
        R_ggF_H_tautau_ATLAS8=(1+(ggF_H_tautau_TH8-ip_ex_gg_phi_tautau_ATLAS8(mHh))/ip_ex_gg_phi_tautau_ATLAS8_e(mHh) ) * nftos;
        THoEX_bbF_H_tautau_ATLAS8=bbF_H_tautau_TH8/ip_ex_bb_phi_tautau_ATLAS8(mHh);
        R_bbF_H_tautau_ATLAS8=(1+(bbF_H_tautau_TH8-ip_ex_bb_phi_tautau_ATLAS8(mHh))/ip_ex_bb_phi_tautau_ATLAS8_e(mHh) ) * nftos;
    }
    else if(mHh>=145.0 && mHh<150.0)
    {
        THoEX_ggF_H_tautau_ATLAS8=ggF_H_tautau_TH8/ip_ex_gg_phi_tautau_ATLAS8(mHh);
        R_ggF_H_tautau_ATLAS8=(1+(ggF_H_tautau_TH8-ip_ex_gg_phi_tautau_ATLAS8(mHh))/ip_ex_gg_phi_tautau_ATLAS8_e(mHh) ) * nftos;
        THoEX_bbF_H_tautau_ATLAS8=bbF_H_tautau_TH8/ip_ex_bb_phi_tautau_ATLAS8(mHh);
        R_bbF_H_tautau_ATLAS8=(1+(bbF_H_tautau_TH8-ip_ex_bb_phi_tautau_ATLAS8(mHh))/ip_ex_bb_phi_tautau_ATLAS8_e(mHh) ) * nftos;
    }
    else if(mHh>=150.0 && mHh<175.0)
    {
        THoEX_ggF_H_tautau_ATLAS8=ggF_H_tautau_TH8/ip_ex_gg_phi_tautau_ATLAS8(mHh);
        R_ggF_H_tautau_ATLAS8=(1+(ggF_H_tautau_TH8-ip_ex_gg_phi_tautau_ATLAS8(mHh))/ip_ex_gg_phi_tautau_ATLAS8_e(mHh) ) * nftos;
        THoEX_bbF_H_tautau_ATLAS8=bbF_H_tautau_TH8/ip_ex_bb_phi_tautau_ATLAS8(mHh);
        R_bbF_H_tautau_ATLAS8=(1+(bbF_H_tautau_TH8-ip_ex_bb_phi_tautau_ATLAS8(mHh))/ip_ex_bb_phi_tautau_ATLAS8_e(mHh) ) * nftos;
    }
    else if(mHh>=175.0 && mHh<200.0)
    {
        THoEX_ggF_H_tautau_ATLAS8=ggF_H_tautau_TH8/ip_ex_gg_phi_tautau_ATLAS8(mHh);
        R_ggF_H_tautau_ATLAS8=(1+(ggF_H_tautau_TH8-ip_ex_gg_phi_tautau_ATLAS8(mHh))/ip_ex_gg_phi_tautau_ATLAS8_e(mHh) ) * nftos;
        THoEX_bbF_H_tautau_ATLAS8=bbF_H_tautau_TH8/ip_ex_bb_phi_tautau_ATLAS8(mHh);
        R_bbF_H_tautau_ATLAS8=(1+(bbF_H_tautau_TH8-ip_ex_bb_phi_tautau_ATLAS8(mHh))/ip_ex_bb_phi_tautau_ATLAS8_e(mHh) ) * nftos;
    }
    else if(mHh>=200.0 && mHh<220.0)
    {
        THoEX_ggF_H_tautau_ATLAS8=ggF_H_tautau_TH8/ip_ex_gg_phi_tautau_ATLAS8(mHh);
        R_ggF_H_tautau_ATLAS8=(1+(ggF_H_tautau_TH8-ip_ex_gg_phi_tautau_ATLAS8(mHh))/ip_ex_gg_phi_tautau_ATLAS8_e(mHh) ) * nftos;
        THoEX_bbF_H_tautau_ATLAS8=bbF_H_tautau_TH8/ip_ex_bb_phi_tautau_ATLAS8(mHh);
        R_bbF_H_tautau_ATLAS8=(1+(bbF_H_tautau_TH8-ip_ex_bb_phi_tautau_ATLAS8(mHh))/ip_ex_bb_phi_tautau_ATLAS8_e(mHh) ) * nftos;
    }
    else if(mHh>=220.0 && mHh<250.0)
    {
        THoEX_ggF_H_tautau_ATLAS8=ggF_H_tautau_TH8/ip_ex_gg_phi_tautau_ATLAS8(mHh);
        R_ggF_H_tautau_ATLAS8=(1+(ggF_H_tautau_TH8-ip_ex_gg_phi_tautau_ATLAS8(mHh))/ip_ex_gg_phi_tautau_ATLAS8_e(mHh) ) * nftos;
        THoEX_bbF_H_tautau_ATLAS8=bbF_H_tautau_TH8/ip_ex_bb_phi_tautau_ATLAS8(mHh);
        R_bbF_H_tautau_ATLAS8=(1+(bbF_H_tautau_TH8-ip_ex_bb_phi_tautau_ATLAS8(mHh))/ip_ex_bb_phi_tautau_ATLAS8_e(mHh) ) * nftos;
    }
    else if(mHh>=250.0 && mHh<260.0)
    {
        THoEX_ggF_H_tautau_ATLAS8=ggF_H_tautau_TH8/ip_ex_gg_phi_tautau_ATLAS8(mHh);
        R_ggF_H_tautau_ATLAS8=(1+(ggF_H_tautau_TH8-ip_ex_gg_phi_tautau_ATLAS8(mHh))/ip_ex_gg_phi_tautau_ATLAS8_e(mHh) ) * nftos;
        THoEX_bbF_H_tautau_ATLAS8=bbF_H_tautau_TH8/ip_ex_bb_phi_tautau_ATLAS8(mHh);
        R_bbF_H_tautau_ATLAS8=(1+(bbF_H_tautau_TH8-ip_ex_bb_phi_tautau_ATLAS8(mHh))/ip_ex_bb_phi_tautau_ATLAS8_e(mHh) ) * nftos;
    }
    else if(mHh>=260.0 && mHh<270.0)
    {
        THoEX_ggF_H_tautau_ATLAS8=ggF_H_tautau_TH8/ip_ex_gg_phi_tautau_ATLAS8(mHh);
        R_ggF_H_tautau_ATLAS8=(1+(ggF_H_tautau_TH8-ip_ex_gg_phi_tautau_ATLAS8(mHh))/ip_ex_gg_phi_tautau_ATLAS8_e(mHh) ) * nftos;
        THoEX_bbF_H_tautau_ATLAS8=bbF_H_tautau_TH8/ip_ex_bb_phi_tautau_ATLAS8(mHh);
        R_bbF_H_tautau_ATLAS8=(1+(bbF_H_tautau_TH8-ip_ex_bb_phi_tautau_ATLAS8(mHh))/ip_ex_bb_phi_tautau_ATLAS8_e(mHh) ) * nftos;
        THoEX_pp_H_hh_bbbb_CMS13=pp_H_hh_bbbb_TH13/ip_ex_pp_H_hh_bbbb_CMS13(mHh);
        R_pp_H_hh_bbbb_CMS13=(1+(pp_H_hh_bbbb_TH13-ip_ex_pp_H_hh_bbbb_CMS13(mHh))/ip_ex_pp_H_hh_bbbb_CMS13_e(mHh) ) * nftos;
    }
    else if(mHh>=270.0 && mHh<275.0)
    {
        THoEX_ggF_H_tautau_ATLAS8=ggF_H_tautau_TH8/ip_ex_gg_phi_tautau_ATLAS8(mHh);
        R_ggF_H_tautau_ATLAS8=(1+(ggF_H_tautau_TH8-ip_ex_gg_phi_tautau_ATLAS8(mHh))/ip_ex_gg_phi_tautau_ATLAS8_e(mHh) ) * nftos;
        THoEX_bbF_H_tautau_ATLAS8=bbF_H_tautau_TH8/ip_ex_bb_phi_tautau_ATLAS8(mHh);
        R_bbF_H_tautau_ATLAS8=(1+(bbF_H_tautau_TH8-ip_ex_bb_phi_tautau_ATLAS8(mHh))/ip_ex_bb_phi_tautau_ATLAS8_e(mHh) ) * nftos;
        THoEX_pp_H_hh_bbbb_CMS13=pp_H_hh_bbbb_TH13/ip_ex_pp_H_hh_bbbb_CMS13(mHh);
        R_pp_H_hh_bbbb_CMS13=(1+(pp_H_hh_bbbb_TH13-ip_ex_pp_H_hh_bbbb_CMS13(mHh))/ip_ex_pp_H_hh_bbbb_CMS13_e(mHh) ) * nftos;
    }
    else if(mHh>=275.0 && mHh<300.0)
    {
        THoEX_ggF_H_tautau_ATLAS8=ggF_H_tautau_TH8/ip_ex_gg_phi_tautau_ATLAS8(mHh);
        R_ggF_H_tautau_ATLAS8=(1+(ggF_H_tautau_TH8-ip_ex_gg_phi_tautau_ATLAS8(mHh))/ip_ex_gg_phi_tautau_ATLAS8_e(mHh) ) * nftos;
        THoEX_bbF_H_tautau_ATLAS8=bbF_H_tautau_TH8/ip_ex_bb_phi_tautau_ATLAS8(mHh);
        R_bbF_H_tautau_ATLAS8=(1+(bbF_H_tautau_TH8-ip_ex_bb_phi_tautau_ATLAS8(mHh))/ip_ex_bb_phi_tautau_ATLAS8_e(mHh) ) * nftos;
        THoEX_pp_H_hh_bbbb_CMS13=pp_H_hh_bbbb_TH13/ip_ex_pp_H_hh_bbbb_CMS13(mHh);
        R_pp_H_hh_bbbb_CMS13=(1+(pp_H_hh_bbbb_TH13-ip_ex_pp_H_hh_bbbb_CMS13(mHh))/ip_ex_pp_H_hh_bbbb_CMS13_e(mHh) ) * nftos;
    }
    else if(mHh>=300.0 && mHh<350.0)
    {
        THoEX_ggF_H_tautau_ATLAS8=ggF_H_tautau_TH8/ip_ex_gg_phi_tautau_ATLAS8(mHh);
        R_ggF_H_tautau_ATLAS8=(1+(ggF_H_tautau_TH8-ip_ex_gg_phi_tautau_ATLAS8(mHh))/ip_ex_gg_phi_tautau_ATLAS8_e(mHh) ) * nftos;
        THoEX_bbF_H_tautau_ATLAS8=bbF_H_tautau_TH8/ip_ex_bb_phi_tautau_ATLAS8(mHh);
        R_bbF_H_tautau_ATLAS8=(1+(bbF_H_tautau_TH8-ip_ex_bb_phi_tautau_ATLAS8(mHh))/ip_ex_bb_phi_tautau_ATLAS8_e(mHh) ) * nftos;
        THoEX_pp_H_hh_bbbb_CMS13=pp_H_hh_bbbb_TH13/ip_ex_pp_H_hh_bbbb_CMS13(mHh);
        R_pp_H_hh_bbbb_CMS13=(1+(pp_H_hh_bbbb_TH13-ip_ex_pp_H_hh_bbbb_CMS13(mHh))/ip_ex_pp_H_hh_bbbb_CMS13_e(mHh) ) * nftos;
    }
    else if(mHh>=350.0 && mHh<400.0)
    {
        THoEX_ggF_H_tautau_ATLAS8=ggF_H_tautau_TH8/ip_ex_gg_phi_tautau_ATLAS8(mHh);
        R_ggF_H_tautau_ATLAS8=(1+(ggF_H_tautau_TH8-ip_ex_gg_phi_tautau_ATLAS8(mHh))/ip_ex_gg_phi_tautau_ATLAS8_e(mHh) ) * nftos;
        THoEX_bbF_H_tautau_ATLAS8=bbF_H_tautau_TH8/ip_ex_bb_phi_tautau_ATLAS8(mHh);
        R_bbF_H_tautau_ATLAS8=(1+(bbF_H_tautau_TH8-ip_ex_bb_phi_tautau_ATLAS8(mHh))/ip_ex_bb_phi_tautau_ATLAS8_e(mHh) ) * nftos;
        THoEX_pp_H_hh_bbbb_CMS13=pp_H_hh_bbbb_TH13/ip_ex_pp_H_hh_bbbb_CMS13(mHh);
        R_pp_H_hh_bbbb_CMS13=(1+(pp_H_hh_bbbb_TH13-ip_ex_pp_H_hh_bbbb_CMS13(mHh))/ip_ex_pp_H_hh_bbbb_CMS13_e(mHh) ) * nftos;
    }
    else if(mHh>=400.0 && mHh<500.0)
    {
        THoEX_ggF_H_tautau_ATLAS8=ggF_H_tautau_TH8/ip_ex_gg_phi_tautau_ATLAS8(mHh);
        R_ggF_H_tautau_ATLAS8=(1+(ggF_H_tautau_TH8-ip_ex_gg_phi_tautau_ATLAS8(mHh))/ip_ex_gg_phi_tautau_ATLAS8_e(mHh) ) * nftos;
        THoEX_bbF_H_tautau_ATLAS8=bbF_H_tautau_TH8/ip_ex_bb_phi_tautau_ATLAS8(mHh);
        R_bbF_H_tautau_ATLAS8=(1+(bbF_H_tautau_TH8-ip_ex_bb_phi_tautau_ATLAS8(mHh))/ip_ex_bb_phi_tautau_ATLAS8_e(mHh) ) * nftos;
        THoEX_pp_H_hh_bbbb_CMS13=pp_H_hh_bbbb_TH13/ip_ex_pp_H_hh_bbbb_CMS13(mHh);
        R_pp_H_hh_bbbb_CMS13=(1+(pp_H_hh_bbbb_TH13-ip_ex_pp_H_hh_bbbb_CMS13(mHh))/ip_ex_pp_H_hh_bbbb_CMS13_e(mHh) ) * nftos;
    }
    else if(mHh>=500.0 && mHh<550.0)
    {
        THoEX_ggF_H_tautau_ATLAS8=ggF_H_tautau_TH8/ip_ex_gg_phi_tautau_ATLAS8(mHh);
        R_ggF_H_tautau_ATLAS8=(1+(ggF_H_tautau_TH8-ip_ex_gg_phi_tautau_ATLAS8(mHh))/ip_ex_gg_phi_tautau_ATLAS8_e(mHh) ) * nftos;
        THoEX_bbF_H_tautau_ATLAS8=bbF_H_tautau_TH8/ip_ex_bb_phi_tautau_ATLAS8(mHh);
        R_bbF_H_tautau_ATLAS8=(1+(bbF_H_tautau_TH8-ip_ex_bb_phi_tautau_ATLAS8(mHh))/ip_ex_bb_phi_tautau_ATLAS8_e(mHh) ) * nftos;
        THoEX_pp_H_hh_bbbb_CMS13=pp_H_hh_bbbb_TH13/ip_ex_pp_H_hh_bbbb_CMS13(mHh);
        R_pp_H_hh_bbbb_CMS13=(1+(pp_H_hh_bbbb_TH13-ip_ex_pp_H_hh_bbbb_CMS13(mHh))/ip_ex_pp_H_hh_bbbb_CMS13_e(mHh) ) * nftos;
    }
    else if(mHh>=550.0 && mHh<600.0)
    {
        THoEX_ggF_H_tautau_ATLAS8=ggF_H_tautau_TH8/ip_ex_gg_phi_tautau_ATLAS8(mHh);
        R_ggF_H_tautau_ATLAS8=(1+(ggF_H_tautau_TH8-ip_ex_gg_phi_tautau_ATLAS8(mHh))/ip_ex_gg_phi_tautau_ATLAS8_e(mHh) ) * nftos;
        THoEX_bbF_H_tautau_ATLAS8=bbF_H_tautau_TH8/ip_ex_bb_phi_tautau_ATLAS8(mHh);
        R_bbF_H_tautau_ATLAS8=(1+(bbF_H_tautau_TH8-ip_ex_bb_phi_tautau_ATLAS8(mHh))/ip_ex_bb_phi_tautau_ATLAS8_e(mHh) ) * nftos;
        THoEX_pp_H_hh_bbbb_CMS13=pp_H_hh_bbbb_TH13/ip_ex_pp_H_hh_bbbb_CMS13(mHh);
        R_pp_H_hh_bbbb_CMS13=(1+(pp_H_hh_bbbb_TH13-ip_ex_pp_H_hh_bbbb_CMS13(mHh))/ip_ex_pp_H_hh_bbbb_CMS13_e(mHh) ) * nftos;
    }
    else if(mHh>=600.0 && mHh<650.0)
    {
        THoEX_ggF_H_tautau_ATLAS8=ggF_H_tautau_TH8/ip_ex_gg_phi_tautau_ATLAS8(mHh);
        R_ggF_H_tautau_ATLAS8=(1+(ggF_H_tautau_TH8-ip_ex_gg_phi_tautau_ATLAS8(mHh))/ip_ex_gg_phi_tautau_ATLAS8_e(mHh) ) * nftos;
        THoEX_bbF_H_tautau_ATLAS8=bbF_H_tautau_TH8/ip_ex_bb_phi_tautau_ATLAS8(mHh);
        R_bbF_H_tautau_ATLAS8=(1+(bbF_H_tautau_TH8-ip_ex_bb_phi_tautau_ATLAS8(mHh))/ip_ex_bb_phi_tautau_ATLAS8_e(mHh) ) * nftos;
        THoEX_pp_H_hh_bbbb_CMS13=pp_H_hh_bbbb_TH13/ip_ex_pp_H_hh_bbbb_CMS13(mHh);
        R_pp_H_hh_bbbb_CMS13=(1+(pp_H_hh_bbbb_TH13-ip_ex_pp_H_hh_bbbb_CMS13(mHh))/ip_ex_pp_H_hh_bbbb_CMS13_e(mHh) ) * nftos;
    }
    else if(mHh>=650.0 && mHh<760.0)
    {
        THoEX_ggF_H_tautau_ATLAS8=ggF_H_tautau_TH8/ip_ex_gg_phi_tautau_ATLAS8(mHh);
        R_ggF_H_tautau_ATLAS8=(1+(ggF_H_tautau_TH8-ip_ex_gg_phi_tautau_ATLAS8(mHh))/ip_ex_gg_phi_tautau_ATLAS8_e(mHh) ) * nftos;
        THoEX_bbF_H_tautau_ATLAS8=bbF_H_tautau_TH8/ip_ex_bb_phi_tautau_ATLAS8(mHh);
        R_bbF_H_tautau_ATLAS8=(1+(bbF_H_tautau_TH8-ip_ex_bb_phi_tautau_ATLAS8(mHh))/ip_ex_bb_phi_tautau_ATLAS8_e(mHh) ) * nftos;
        THoEX_pp_H_hh_bbbb_CMS13=pp_H_hh_bbbb_TH13/ip_ex_pp_H_hh_bbbb_CMS13(mHh);
        R_pp_H_hh_bbbb_CMS13=(1+(pp_H_hh_bbbb_TH13-ip_ex_pp_H_hh_bbbb_CMS13(mHh))/ip_ex_pp_H_hh_bbbb_CMS13_e(mHh) ) * nftos;
    }
    else if(mHh>=760.0 && mHh<850.0)
    {
        THoEX_ggF_H_tautau_ATLAS8=ggF_H_tautau_TH8/ip_ex_gg_phi_tautau_ATLAS8(mHh);
        R_ggF_H_tautau_ATLAS8=(1+(ggF_H_tautau_TH8-ip_ex_gg_phi_tautau_ATLAS8(mHh))/ip_ex_gg_phi_tautau_ATLAS8_e(mHh) ) * nftos;
        THoEX_bbF_H_tautau_ATLAS8=bbF_H_tautau_TH8/ip_ex_bb_phi_tautau_ATLAS8(mHh);
        R_bbF_H_tautau_ATLAS8=(1+(bbF_H_tautau_TH8-ip_ex_bb_phi_tautau_ATLAS8(mHh))/ip_ex_bb_phi_tautau_ATLAS8_e(mHh) ) * nftos;
        THoEX_pp_H_hh_bbbb_CMS13=pp_H_hh_bbbb_TH13/ip_ex_pp_H_hh_bbbb_CMS13(mHh);
        R_pp_H_hh_bbbb_CMS13=(1+(pp_H_hh_bbbb_TH13-ip_ex_pp_H_hh_bbbb_CMS13(mHh))/ip_ex_pp_H_hh_bbbb_CMS13_e(mHh) ) * nftos;
    }
    else if(mHh>=850.0 && mHh<900.0)
    {
        THoEX_ggF_H_tautau_ATLAS8=ggF_H_tautau_TH8/ip_ex_gg_phi_tautau_ATLAS8(mHh);
        R_ggF_H_tautau_ATLAS8=(1+(ggF_H_tautau_TH8-ip_ex_gg_phi_tautau_ATLAS8(mHh))/ip_ex_gg_phi_tautau_ATLAS8_e(mHh) ) * nftos;
        THoEX_bbF_H_tautau_ATLAS8=bbF_H_tautau_TH8/ip_ex_bb_phi_tautau_ATLAS8(mHh);
        R_bbF_H_tautau_ATLAS8=(1+(bbF_H_tautau_TH8-ip_ex_bb_phi_tautau_ATLAS8(mHh))/ip_ex_bb_phi_tautau_ATLAS8_e(mHh) ) * nftos;
        THoEX_pp_H_hh_bbbb_CMS13=pp_H_hh_bbbb_TH13/ip_ex_pp_H_hh_bbbb_CMS13(mHh);
        R_pp_H_hh_bbbb_CMS13=(1+(pp_H_hh_bbbb_TH13-ip_ex_pp_H_hh_bbbb_CMS13(mHh))/ip_ex_pp_H_hh_bbbb_CMS13_e(mHh) ) * nftos;
    }
    else if(mHh>=900.0 && mHh<1000.0)
    {
        THoEX_ggF_H_tautau_ATLAS8=ggF_H_tautau_TH8/ip_ex_gg_phi_tautau_ATLAS8(mHh);
        R_ggF_H_tautau_ATLAS8=(1+(ggF_H_tautau_TH8-ip_ex_gg_phi_tautau_ATLAS8(mHh))/ip_ex_gg_phi_tautau_ATLAS8_e(mHh) ) * nftos;
        THoEX_bbF_H_tautau_ATLAS8=bbF_H_tautau_TH8/ip_ex_bb_phi_tautau_ATLAS8(mHh);
        R_bbF_H_tautau_ATLAS8=(1+(bbF_H_tautau_TH8-ip_ex_bb_phi_tautau_ATLAS8(mHh))/ip_ex_bb_phi_tautau_ATLAS8_e(mHh) ) * nftos;
        THoEX_pp_H_hh_bbbb_CMS13=pp_H_hh_bbbb_TH13/ip_ex_pp_H_hh_bbbb_CMS13(mHh);
        R_pp_H_hh_bbbb_CMS13=(1+(pp_H_hh_bbbb_TH13-ip_ex_pp_H_hh_bbbb_CMS13(mHh))/ip_ex_pp_H_hh_bbbb_CMS13_e(mHh) ) * nftos;
    }
    else if(mHh>=1000.0 && mHh<1100.0)
    {
        THoEX_pp_H_hh_bbbb_CMS13=pp_H_hh_bbbb_TH13/ip_ex_pp_H_hh_bbbb_CMS13(mHh);
        R_pp_H_hh_bbbb_CMS13=(1+(pp_H_hh_bbbb_TH13-ip_ex_pp_H_hh_bbbb_CMS13(mHh))/ip_ex_pp_H_hh_bbbb_CMS13_e(mHh) ) * nftos;
    }
    else if(mHh>=1100.0 && mHh<1200.0)
    {
        THoEX_pp_H_hh_bbbb_CMS13=pp_H_hh_bbbb_TH13/ip_ex_pp_H_hh_bbbb_CMS13(mHh);
        R_pp_H_hh_bbbb_CMS13=(1+(pp_H_hh_bbbb_TH13-ip_ex_pp_H_hh_bbbb_CMS13(mHh))/ip_ex_pp_H_hh_bbbb_CMS13_e(mHh) ) * nftos;
    }

}    


//void GMcache::runGMparameters()
//{
////    vev=myGM->v();
////    double cosb=myGM->getcosb();
//
//    std::string RGEorder=myGM->getRGEorderflag();
//    //flag will be used to transport information about model and RGEorder to the Runner:
//    //flag=0 for LO, 1 for approxNLO (and 2 for NLO - not implemented yet)
//    int flag;
//    if( RGEorder == "LO" ) flag=0;
//    else if( RGEorder == "approxNLO" ) flag=1;
////    else if( RGEorder == "NLO" ) flag=2;
//    else {
//        throw std::runtime_error("RGEorder can be only any of \"LO\", \"approxNLO\" or \"NLO\"");
//    }
//
//    double lambda1_at_MZ=lambda1;
//    double lambda2_at_MZ=lambda2;
//    double lambda3_at_MZ=lambda3;
//    double lambda4_at_MZ=lambda4;
//    double mu1_at_MZ=mu1;
//    double mu3_at_MZ=mu3;
//    double mu4_at_MZ=mu4;
//    double nu1_at_MZ=nu1;
//    double omega1_at_MZ=omega1;
//    double kappa1_at_MZ=kappa1;
//    double nu2_at_MZ=nu2;
//    double omega2_at_MZ=omega2;
//    double kappa2_at_MZ=kappa2;
//    double nu4_at_MZ=nu4;
//    double omega4_at_MZ=omega4;
//    double NLOuniscale=myGM->getNLOuniscaleGM();
//
//    if(fabs(Q_GM-log10(MZ))<0.005)   //at MZ scale
//    {
//        Q_cutoff=log10(MZ);
//
//        lambda1_at_Q = lambda1_at_MZ;
//        lambda2_at_Q = lambda2_at_MZ;
//        lambda3_at_Q = lambda3_at_MZ;
//        lambda4_at_Q = lambda4_at_MZ;
//        mu1_at_Q = mu1_at_MZ;
//        mu3_at_Q = mu3_at_MZ;
//        mu4_at_Q = mu4_at_MZ;
//        nu1_at_Q = nu1_at_MZ;
//        omega1_at_Q = omega1_at_MZ;
//        kappa1_at_Q = kappa1_at_MZ;
//        nu2_at_Q = nu2_at_MZ;
//        omega2_at_Q = omega2_at_MZ;
//        kappa2_at_Q = kappa2_at_MZ;
//        nu4_at_Q = nu4_at_MZ;
//        omega4_at_Q = omega4_at_MZ;
//    }
//    else   //at some other scale
//    {
//        double InitVals[15];
//        InitVals[0]=lambda1_at_MZ;
//        InitVals[1]=lambda2_at_MZ;
//        InitVals[2]=lambda3_at_MZ;
//        InitVals[3]=lambda4_at_MZ;
//        InitVals[4]=mu1_at_MZ;
//        InitVals[5]=mu3_at_MZ;
//        InitVals[6]=mu4_at_MZ;
//        InitVals[7]=nu1_at_MZ;
//        InitVals[8]=omega1_at_MZ;
//        InitVals[9]=kappa1_at_MZ;
//        InitVals[10]=nu2_at_MZ;
//        InitVals[11]=omega2_at_MZ;
//        InitVals[12]=kappa2_at_MZ;
//        InitVals[13]=nu4_at_MZ;
//        InitVals[14]=omega4_at_MZ;
//
//        Q_cutoff=myRunnerGM->RGERunnerGM(InitVals, 15, log10(MZ), Q_GM, flag, RpepsGM, NLOuniscale);  //Running up to Q_cutoff<=Q_GM
//
//        lambda1_at_Q = InitVals[0];
//        lambda2_at_Q = InitVals[1];
//        lambda3_at_Q = InitVals[2];
//        lambda4_at_Q = InitVals[3];
//        mu1_at_Q=InitVals[4];
//        mu3_at_Q=InitVals[5];
//        mu4_at_Q = InitVals[6];
//        nu1_at_Q = InitVals[7];
//        omega1_at_Q = InitVals[8];
//        kappa1_at_Q = InitVals[9];
//        nu2_at_Q = InitVals[10];
//        omega2_at_Q = InitVals[11];
//        kappa2_at_Q = InitVals[12];
//        nu4_at_Q = InitVals[13];
//        omega4_at_Q = InitVals[14];
//    }
//
//}

void GMcache::computeUnitarity()
{
//    std::string ModelType=myGM->getModelTypeGMflag();
//    if( ModelType != "custodial1" )
//    {
//        throw std::runtime_error("GM unitarity constraints are only implemented for the \"custodial1\" model.");
//    }
//
    double pi=M_PI;
//    /*
//    *******   LO part   *************
//    */

    // Taken from 0712.4053v2, eq. (3.96) to (3.107)
    // (all unitarityeigenvalues should be smaller than 0.5 in magnitude)
    // Note that when going from their convention to ours, one has to swap lambda3 and lambda4.
    unitarityeigenvalues.assign(0, (6.0*lambda1+11.0*lambda2+7.0*lambda3+sqrt((6.0*lambda1-11.0*lambda2-7.0*lambda3)*(6.0*lambda1-11.0*lambda2-7.0*lambda3)+36.0*lambda4*lambda4))/(8.0*pi));
    unitarityeigenvalues.assign(1, (6.0*lambda1+11.0*lambda2+7.0*lambda3-sqrt((6.0*lambda1-11.0*lambda2-7.0*lambda3)*(6.0*lambda1-11.0*lambda2-7.0*lambda3)+36.0*lambda4*lambda4))/(8.0*pi));
    unitarityeigenvalues.assign(2, (2.0*lambda1+2.0*lambda2-lambda3+sqrt((2.0*lambda1-2.0*lambda2+lambda3)*(2.0*lambda1-2.0*lambda2+lambda3)+lambda5*lambda5))/(8.0*pi));
    unitarityeigenvalues.assign(3, (2.0*lambda1+2.0*lambda2-lambda3-sqrt((2.0*lambda1-2.0*lambda2+lambda3)*(2.0*lambda1-2.0*lambda2+lambda3)+lambda5*lambda5))/(8.0*pi));
    unitarityeigenvalues.assign(4, (2.0*lambda1+2.0*lambda2+sqrt((2.0*lambda1-2.0*lambda2)*(2.0*lambda1-2.0*lambda2)+lambda5*lambda5))/(8.0*pi));
    unitarityeigenvalues.assign(5, (2.0*lambda1+2.0*lambda2-sqrt((2.0*lambda1-2.0*lambda2)*(2.0*lambda1-2.0*lambda2)+lambda5*lambda5))/(8.0*pi));
    unitarityeigenvalues.assign(6, (4.0*lambda1+2.0*lambda2-lambda3+sqrt((4.0*lambda1-2.0*lambda2+lambda3)*(4.0*lambda1-2.0*lambda2+lambda3)+2.0*lambda5*lambda5))/(8.0*pi));
    unitarityeigenvalues.assign(7, (4.0*lambda1+2.0*lambda2-lambda3-sqrt((4.0*lambda1-2.0*lambda2+lambda3)*(4.0*lambda1-2.0*lambda2+lambda3)+2.0*lambda5*lambda5))/(8.0*pi));
    unitarityeigenvalues.assign(8, (6.0*lambda2+7.0*lambda3+sqrt(4.0*lambda2*lambda2+4.0*lambda2*lambda3+17.0*lambda3*lambda3))/(8.0*pi));
    unitarityeigenvalues.assign(9, (6.0*lambda2+7.0*lambda3-sqrt(4.0*lambda2*lambda2+4.0*lambda2*lambda3+17.0*lambda3*lambda3))/(8.0*pi));
    unitarityeigenvalues.assign(10, (lambda2+2.0*lambda3)/(2.0*pi));
    unitarityeigenvalues.assign(11, (2.0*lambda2+lambda3)/(4.0*pi));
    unitarityeigenvalues.assign(12, (4.0*lambda4+lambda5)/(16.0*pi));
    unitarityeigenvalues.assign(13, (2.0*lambda4-lambda5)/(8.0*pi));
    unitarityeigenvalues.assign(14, (lambda4+lambda5)/(4.0*pi));
    unitarityeigenvalues.assign(15, (2.0*lambda2+(2.0*sqrt(2.0))*lambda3)/(4.0*pi));
    unitarityeigenvalues.assign(16, (2.0*lambda2+(2.0*sqrt(2.0))*lambda3)/(4.0*pi));

    /*
    *******   NLO part   *************
    */

    // beta functions taken from 
//    double blambda1=(12.0*lambda1*lambda1 + 4.0*lambda3*lambda3 + 4.0*lambda3*lambda4 + 4.0*lambda4*lambda4 
//                     + 8.0*nu1*nu1 + 8.0*nu1*nu2 + 8.0*nu2*nu2)/(16.0*pi*pi);
//    double blambda2=(12.0*lambda2*lambda2 + 4.0*lambda3*lambda3 + 4.0*lambda3*lambda4 + 4.0*lambda4*lambda4
//                     + 8.0*omega1*omega1 + 8.0*omega1*omega2 + 8.0*omega2*omega2)/(16.0*pi*pi);
//    double blambda3=(4.0*lambda3*lambda3 + 4.0*lambda4*lambda4 + (lambda1+lambda2)*(6.0*lambda3+2.0*lambda4) 
//                     + 8.0*kappa2*kappa2 + 8.0*nu1*omega1 + 4.0*nu2*omega1 + 4.0*nu1*omega2)/(16.0*pi*pi);
//    double blambda4=(lambda1*lambda4 + lambda2*lambda4 + 4.0*lambda3*lambda4 + 6.0*lambda4*lambda4
//                     + 4.0*kappa1*kappa1 + 4.0*kappa1*kappa2 + 2.0*kappa2*kappa2 + 2.0*nu2*omega2)/(8.0*pi*pi);
//    double bmu1=(11.0*mu1*mu1 + 3.0*mu1*mu4 + mu1*(2.0*mu1+6.0*mu3+3.0*mu4)
//               + 3.0*nu4*nu4 + 3.0*omega4*omega4)/(16.0*pi*pi);
//    double bmu3=(18.0*kappa1*kappa1 + 18.0*kappa1*kappa2 + 134.0*mu1*mu1 + 6.0*mu1*(39.0*mu3 + 22.0*mu4)
//               + 3.0*(30.0*mu3*mu3 + 39.0*mu3*mu4 + 9.0*mu4*mu4 
//                      + 3.0*nu1*nu1 + 3.0*nu1*nu2 - 5.0*nu4*nu4
//                      + 3.0*omega1*omega1 + 3.0*omega1*omega2 - 5.0*omega4*omega4))/(72.0*pi*pi);
//    double bmu4=(18.0*kappa2*kappa2 + 4.0*mu1*mu1 + 156.0*mu1*mu4 + 54.0*mu3*mu4 + 144.0*mu4*mu4
//               + 9.0*nu2*nu2 + 6.0*nu4*nu4 + 9.0*omega2*omega2 + 6.0*omega4*omega4)/(144.0*pi*pi);
//    double bnu1=(6.0*kappa1*kappa1 + 6.0*kappa2*kappa2 + 18.0*lambda1*nu1
//               + 78.0*mu1*nu1 + 51.0*mu3*nu1 + 39.0*mu4*nu1 + 6.0*nu1*nu1
//               + 6.0*lambda1*nu2 + 32.0*mu1*nu2 + 24.0*mu3*nu2 + 6.0*mu4*nu2
//               + 6.0*nu2*nu2 + 10.0*nu4*nu4
//               + 12.0*lambda3*omega1 + 6.0*lambda4*omega1 + 6.0*lambda3*omega2)/(48.0*pi*pi);
//    double bomega1=(6.0*kappa1*kappa1 + 6.0*kappa2*kappa2 
//               + 12.0*lambda3*nu1 + 6.0*lambda4*nu1 + 6.0*lambda3*nu2
//               + 18.0*lambda2*omega1 + 78.0*mu1*omega1 + 51.0*mu3*omega1 + 39.0*mu4*omega1 + 6.0*omega1*omega1
//               + 6.0*lambda2*omega2 + 32.0*mu1*omega2 + 24.0*mu3*omega2 + 6.0*mu4*omega2 + 6.0*omega2*omega2
//               + 10.0*omega4*omega4)/(48.0*pi*pi);
//    double bkappa1=(6.0*kappa1*(2.0*lambda3 + 10.0*lambda4 + 18.0*mu1 + 17.0*mu3 + 13.0*mu4 + 2.0*nu1 + 2.0*omega1)
//               + kappa2*(24.0*lambda4 + 64.0*mu1 + 48.0*mu3 + 24.0*mu4 + 9.0*nu2 + 9.0*omega2)
//               + 20.0*nu4*omega4)/(96.0*pi*pi);
//    double bnu2=(4.0*kappa1*kappa2 + 6.0*kappa2*kappa2 + 2.0*lambda1*nu2 + ((14.0*mu1)/3.0 + mu3 + 9.0*mu4)*nu2 
//                + 4.0*nu1*nu2 + 6.0*nu2*nu2 + (25.0*nu4*nu4)/3.0 + 2.0*lambda4*omega2)/(16.0*pi*pi);
//    double bomega2=(4.0*kappa1*kappa2 + 6.0*kappa2*kappa2 + 2.0*lambda4*nu2 + 2.0*lambda2*omega2 
//                + ((14.0*mu1)/3.0 + mu3 + 9.0*mu4)*omega2 + 4.0*omega1*omega2 + 6.0*omega2*omega2 
//                + (25.0*omega4*omega4)/3.0)/(16.0*pi*pi);
//    double bkappa2=(kappa2*(6.0*lambda3 + 6.0*lambda4 + 14.0*mu1 + 3.0*mu3 + 27.0*mu4
//                     + 6.0*nu1 + 12.0*nu2 + 6.0*omega1 + 12.0*omega2)
//                + 6.0*kappa1*(nu2 + omega2) + 42.0*nu4*omega4)/(48.0*pi*pi);
//    double bnu4=(11.0*mu1*nu4 + 3.0*mu3*nu4 + 9.0*mu4*nu4 + 3.0*nu1*nu4 + 9.0*nu2*nu4 
//                + 3.0*kappa1*omega4 + 9.0*kappa2*omega4)/(16.0*pi*pi);
//    double bomega4=(3.0*kappa1*nu4 + 9.0*kappa2*nu4 
//                + (11.0*mu1 + 3.0*(mu3 + 3.0*mu4 + omega1 + 3.0*omega2))*omega4)/(16.0*pi*pi);
//
//    Sbmatrix1.assign(0,0, 3.0*blambda1/(16.0*pi));
//    Sbmatrix1.assign(0,1, (2.0*blambda3+blambda4)/(16.0*pi));
//    Sbmatrix1.assign(1,0, Sbmatrix1(0,1));
//    Sbmatrix1.assign(0,3, (2.0*bnu1+bnu2)/(8.0*sqrt(2.0)*pi));
//    Sbmatrix1.assign(3,0, Sbmatrix1(0,3));
//    Sbmatrix1.assign(1,1, 3.0*blambda2/(16.0*pi));
//    Sbmatrix1.assign(1,3, (2.0*bomega1+bomega2)/(8.0*sqrt(2.0)*pi));
//    Sbmatrix1.assign(3,1, Sbmatrix1(1,3));
//    Sbmatrix1.assign(2,2, (blambda3+5.0*blambda4)/(16.0*pi));
//    Sbmatrix1.assign(2,3, (4.0*bkappa1+2.0*bkappa2)/(16.0*pi));
//    Sbmatrix1.assign(3,2, Sbmatrix1(2,3));
//    Sbmatrix1.assign(3,3, (26.0*bmu1+17.0*bmu3+13.0*bmu4)/(32.0*pi));
//
//    Sbmatrix2.assign(0,0, blambda1/(16.0*pi));
//    Sbmatrix2.assign(0,1, blambda4/(16.0*pi));
//    Sbmatrix2.assign(1,0, Sbmatrix2(0,1));
//    Sbmatrix2.assign(0,3, bnu2/(8.0*sqrt(2.0)*pi));
//    Sbmatrix2.assign(3,0, Sbmatrix2(0,3));
//    Sbmatrix2.assign(1,1, blambda2/(16.0*pi));
//    Sbmatrix2.assign(1,3, bomega2/(8.0*sqrt(2.0)*pi));
//    Sbmatrix2.assign(3,1, Sbmatrix2(1,3));
//    Sbmatrix2.assign(2,2, (blambda3+blambda4)/(16.0*pi));
//    Sbmatrix2.assign(2,3, bkappa2/(8.0*pi));
//    Sbmatrix2.assign(3,2, Sbmatrix2(2,3));
//    Sbmatrix2.assign(3,3, (14.0*bmu1+3.0*bmu3+27.0*bmu4)/(96.0*pi));
//
//    Seigenvectors1T=Seigenvectors1.hconjugate();
//    Seigenvectors2T=Seigenvectors2.hconjugate();
//
//    for (int i=0; i < 4; i++) {
//        for (int k=0; k < 4; k++) {
//            for (int l=0; l < 4; l++) {
//                Sbeigenvalues1.assign(i, Sbeigenvalues1(i) + Seigenvectors1T(i,k) * Sbmatrix1(k,l) * Seigenvectors1(l,i) );
//                Sbeigenvalues2.assign(i, Sbeigenvalues2(i) + Seigenvectors2T(i,k) * Sbmatrix2(k,l) * Seigenvectors2(l,i) );
//            }                
//        }
//        betaeigenvalues.assign(i, -1.5 * Sbeigenvalues1(i));
//        betaeigenvalues.assign(i+4, -1.5 * Sbeigenvalues2(i));
//    }
//
//    betaeigenvalues.assign(8, -1.5 * (blambda3-blambda4)/(16.0*pi));
//    betaeigenvalues.assign(9, -1.5 * sqrt(15.0)*bnu4/(16.0*pi));
//    betaeigenvalues.assign(10, -1.5 * sqrt(15.0)*bomega4/(16.0*pi));
//
//    for (int i=0; i < 11; i++) {
//        NLOunitarityeigenvalues.assign(i, -(gslpp::complex::i()-1.0/pi)*unitarityeigenvalues(i)*unitarityeigenvalues(i) + betaeigenvalues(i) );
//    }
}

double GMcache::cW2GM(const double c02) const{
    return c02;
}


double GMcache::updateCache()
{
//    GMmodel=myGM->getModelTypeGMflag();
    Q_GM=myGM->getQ_GM();
    vev=myGM->v();
    mHl=myGM->getMHl();
    mHl2=mHl*mHl;
    GF=1/(sqrt(2.0)*vev*vev);
    Ale=myGM->getAle();
    MZ=myGM->getMz();
    MW=myGM->Mw();
    cW2=cW2GM(myGM->c02()); /*This might have to be replaced by the GM corrected value.*/
    tanb=myGM->gettanb();
    sinb=myGM->getsinb();
    cosb=myGM->getcosb();
    sina=myGM->getsina();
    cosa=myGM->getcosa();
    mH1sq=myGM->getinputmHh2();
    mAsq=myGM->getmAsq();
    mH5sq=myGM->getmH5sq();
    Mu1=myGM->getMu1();
    Mu2=myGM->getMu2();
    M1sq=-vev*Mu1/(sqrt(2.0)*cosb);
    M2sq=-sqrt(18.0)*cosb*vev*Mu2;

    vPhi=vev*sinb;
    vDelta=vev*cosb/sqrt(8.0);
    double mHlsq, mHhsq;
    if(mH1sq>=mHl2)
    {
        mHlsq=mHl2;
        mHhsq=mH1sq;
    }
    else
    {
        mHlsq=mH1sq;
        mHhsq=mHl2;
    }

//    double cos2b=cosb*cosb-sinb*sinb;
    lambda1 = (mHlsq*cosa*cosa+mHhsq*sina*sina)/(8.0*vev*vev*sinb*sinb);
    lambda2 = (M2sq+2.0*(mAsq-M1sq)*sinb*sinb
               +2.0/3.0*(mHlsq*sina*sina+mHhsq*cosa*cosa-mH5sq))/(2.0*vev*vev*cosb*cosb);
    lambda3 = (mH5sq-M2sq+(2.0*M1sq-3.0*mAsq)*sinb*sinb)/(vev*vev*cosb*cosb);
    lambda4 = (mAsq-0.5*M1sq+(mHhsq-mHlsq)*sina*cosa/(sqrt(6.0)*sinb*cosb))/(vev*vev);
    lambda5 = 2.0*(M1sq-mAsq)/(vev*vev);

    //triple scalar couplings
    ghhh = -(2.0*sina*sina*sinb*(3.0*cosa+sqrt(6.0)*sina*tanb))/vev * M1sq
           +(sqrt(2.0/3.0)*sina*sina*sina)/(vev*cosb) * M2sq
           +(-3.0*cosa*cosa*cosa/sinb+2.0*sqrt(6.0)*sina*sina*sina/cosb)/vev * mHlsq;
    ghhH = sina*sinb*(4.0-6.0*sina*sina+2.0*sqrt(6.0)*sina*cosa*tanb)/vev * M1sq
            -(sqrt(2.0/3.0)*cosa*sina*sina)/(vev*cosb) * M2sq
            -sina*(cosa*cosa/sinb+2.0*sqrt(2.0/3.0)*sina*cosa/cosb)/vev * (2.0*mHlsq+mHhsq);
    ghHH = cosa*sinb*(4.0-6.0*cosa*cosa-2.0*sqrt(6.0)*sina*cosa*tanb)/vev * M1sq
            +(sqrt(2.0/3.0)*cosa*cosa*sina)/(vev*cosb) * M2sq
            +cosa*(-sina*sina/sinb+2.0*sqrt(2.0/3.0)*sina*cosa/cosb)/vev * (mHlsq+2.0*mHhsq);
    gHHH = (2.0*cosa*cosa*sinb*(-3.0*sina+sqrt(6.0)*cosa*tanb))/vev * M1sq
           -(sqrt(2.0/3.0)*cosa*cosa*cosa)/(vev*cosb) * M2sq
           -(3.0*sina*sina*sina/sinb+2.0*sqrt(6.0)*cosa*cosa*cosa/cosb)/vev * mHhsq;
    ghH3H3 = -(2.0*sqrt(2.0/3.0)*sina)/(vev*cosb) * M1sq
             +(4.0*sqrt(2.0/3.0)*sina*cosb-2.0*cosa*sinb)/vev * mAsq
             +(2.0*sqrt(2.0/3.0)*sina*sinb*tanb-cosa*cosb/tanb)/vev * mHlsq;
    gHH3H3 = (2.0*sqrt(2.0/3.0)*cosa)/(vev*cosb) * M1sq
             -(4.0*sqrt(2.0/3.0)*cosa*cosb+2.0*sina*sinb)/vev * mAsq
             -(2.0*sqrt(2.0/3.0)*cosa*sinb*tanb+sina*cosb/tanb)/vev * mHhsq;
    ghH3pH3m = ghH3H3;
    gHH3pH3m = gHH3H3;
    ghH5H5 = 2.0*sinb*(2.0*cosa+sqrt(6.0)*sina*tanb)/vev * (M1sq-2.0*mAsq)
             +(2.0*sqrt(2.0/3.0)*sina)/(vev*cosb) * (mHlsq + 2.0*mH5sq - M2sq)
             +(2.0*cosa*sinb)/vev * mAsq;
    gHH5H5 = 2.0*sinb*(2.0*sina-sqrt(6.0)*cosa*tanb)/vev * (M1sq-2.0*mAsq)
             -(2.0*sqrt(2.0/3.0)*cosa)/(vev*cosb) * (mHhsq + 2.0*mH5sq - M2sq)
             +(2.0*sina*sinb)/vev * mAsq;
    ghH5pH5m = ghH5H5;
    gHH5pH5m = gHH5H5;
    ghH5ppH5mm = ghH5H5;
    gHH5ppH5mm = gHH5H5;
    gH3H3H5 = 2.0*((mH5sq-2.0*mAsq)*cosb-(2.0*M1sq-3.0*mAsq+mH5sq)/cosb)/(sqrt(3.0)*vev);
    gH5H5H5 = 2.0*((3.0*mH5sq-4.0*M2sq)/cosb+(6.0*M1sq-9.0*mAsq)*sinb*tanb)/(sqrt(3.0)*vev);
    gH3H3pH5m = gslpp::complex::i()*(4.0*M1sq-4.0*mAsq+mH5sq+(2.0*mAsq-mH5sq)*(cosb*cosb-sinb*sinb))/(2.0*vev*cosb);
    gH5H3pH3m = (4.0*M1sq-4.0*mAsq+mH5sq+(2.0*mAsq-mH5sq)*(cosb*cosb-sinb*sinb))/(2.0*sqrt(3.0)*vev*cosb);
    gH5H5pH5m = 0.5*gH5H5H5;
    gH5H5ppH5mm = -gH5H5H5;
    gH5ppH3mH3m = sqrt(6.0)*gH5H3pH3m;
    gH5ppH5mH5m = gH5H5pH5m;

    //    runGMparameters();
    computeUnitarity();
    computeSignalStrengthQuantities();
    computeDirectSearchQuantities();
    return mHl2;
}
