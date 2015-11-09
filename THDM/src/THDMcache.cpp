/* 
 * Copyright (C) 2012 SusyFit Collaboration
 * All rights reserved.
 *
 * For the licensing terms see doc/COPYING.
 */

#include "THDMcache.h"
#include <fstream>
#include "gslpp.h"
#include <sstream>
#include <string>

THDMcache::THDMcache()
        
        : PV(true),
        
        array1(19861, 2, 0.),
        array2(19861, 2, 0.),
        array3(19861, 2, 0.),
        array4(19861, 2, 0.),
        array5(19861, 2, 0.),
        array6(19861, 2, 0.),
        array7(19861, 2, 0.),
        array11(1971, 2, 0.),
        array12(1971, 2, 0.),
        array13(1971, 2, 0.),
        array14(1971, 2, 0.),
        array15(1971, 2, 0.),
        array16(9851, 2, 0.),
        array17(9851, 2, 0.),
        array18(19861, 2, 0.),
        array19(186, 2, 0.),
        array20(186, 2, 0.),
        array21(186, 2, 0.),
        array22(186, 2, 0.),
        array26(992, 2, 0.),
        array27(185, 2, 0.),
        array29(986, 2, 0.),
        array31(985, 2, 0.),
        array32(496, 2, 0.),
        array33(496, 2, 0.),
        array34(991, 2, 0.),
        array35(496, 2, 0.),
        array36(496, 2, 0.),
        array37(100, 2, 0.),
        array38(100, 2, 0.),
        array39(986, 2, 0.),
        array44(986, 2, 0.),
        array45(986, 2, 0.),
        arrayX_bb(199, 2, 0.),
        arrayX_tt(198, 2, 0.) 
{
  read();
//  std::cout<<"read"<<std::endl;
}


/////////////////////////////////////////////////////////////////////////////////////////////////

int THDMcache::CacheCheck(const gslpp::complex cache[][CacheSize], 
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

void THDMcache::CacheShift(gslpp::complex cache[][CacheSize], const int NumPar, 
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

////////////////////////////////////////////////////////////////////////////////
/*One-loop functions*/


gslpp::complex THDMcache::B0_Mz_Mz2_Mz_mH(const double Mz, const double mH) const {
    int NumPar = 2;
    double params[] = {Mz, mH};

    int i = CacheCheck(B0_Mz_Mz2_Mz_mH_cache, NumPar, params);
    if (i>=0) {
        return ( B0_Mz_Mz2_Mz_mH_cache[NumPar][i] );
    } else {
        gslpp::complex newResult = PV.B0(Mz*Mz, Mz*Mz, Mz*Mz, mH*mH);
        CacheShift(B0_Mz_Mz2_Mz_mH_cache, NumPar, params, newResult);
        return newResult;
    } 
}

gslpp::complex THDMcache::B0_Mz_0_Mz_mH(const double Mz, const double mH) const {
    int NumPar = 2;
    double params[] = {Mz, mH};

    int i = CacheCheck(B0_Mz_0_Mz_mH_cache, NumPar, params);
    if (i>=0) {
        return ( B0_Mz_0_Mz_mH_cache[NumPar][i] );
    } else {
        gslpp::complex newResult = PV.B0(Mz*Mz, 0., Mz*Mz, mH*mH);
        CacheShift(B0_Mz_0_Mz_mH_cache, NumPar, params, newResult);
        return newResult;
    } 
}

gslpp::complex THDMcache::B0_Mz_Mz2_Mz_mh(const double Mz, const double mh) const {
    int NumPar = 2;
    double params[] = {Mz, mh};

    int i = CacheCheck(B0_Mz_Mz2_Mz_mh_cache, NumPar, params);
    if (i>=0) {
        return ( B0_Mz_Mz2_Mz_mh_cache[NumPar][i] );
    } else {
        gslpp::complex newResult = PV.B0(Mz*Mz, Mz*Mz, Mz*Mz, mh*mh);
        CacheShift(B0_Mz_Mz2_Mz_mh_cache, NumPar, params, newResult);
        return newResult;
    } 
}

gslpp::complex THDMcache::B0_Mz_0_Mz_mh(const double Mz, const double mh) const {
    int NumPar = 2;
    double params[] = {Mz, mh};

    int i = CacheCheck(B0_Mz_0_Mz_mh_cache, NumPar, params);
    if (i>=0) {
        return ( B0_Mz_0_Mz_mh_cache[NumPar][i] );
    } else {
        gslpp::complex newResult = PV.B0(Mz*Mz, 0., Mz*Mz, mh*mh);
        CacheShift(B0_Mz_0_Mz_mh_cache, NumPar, params, newResult);
        return newResult;
    } 
}

gslpp::complex THDMcache::B0_Mz_Mw2_Mw_mH(const double Mz, const double Mw, const double mH) const {
    int NumPar = 3;
    double params[] = {Mz, Mw, mH};

    int i = CacheCheck(B0_Mz_Mw2_Mw_mH_cache, NumPar, params);
    if (i>=0) {
        return ( B0_Mz_Mw2_Mw_mH_cache[NumPar][i] );
    } else {
        gslpp::complex newResult = PV.B0(Mz*Mz, Mw*Mw, Mw*Mz, mH*mH);
        CacheShift(B0_Mz_Mw2_Mw_mH_cache, NumPar, params, newResult);
        return newResult;
    } 
}

gslpp::complex THDMcache::B0_Mz_0_Mw_mH(const double Mz, const double Mw, const double mH) const {
    int NumPar = 3;
    double params[] = {Mz, Mw, mH};

    int i = CacheCheck(B0_Mz_0_Mw_mH_cache, NumPar, params);
    if (i>=0) {
        return ( B0_Mz_0_Mw_mH_cache[NumPar][i] );
    } else {
        gslpp::complex newResult = PV.B0(Mz*Mz, 0., Mw*Mw, mH*mH);
        CacheShift(B0_Mz_0_Mw_mH_cache, NumPar, params, newResult);
        return newResult;
    } 
}

gslpp::complex THDMcache::B0_Mz_Mw2_Mw_mh(const double Mz, const double Mw, const double mh) const {
    int NumPar = 3;
    double params[] = {Mz, Mw, mh};

    int i = CacheCheck(B0_Mz_Mw2_Mw_mh_cache, NumPar, params);
    if (i>=0) {
        return ( B0_Mz_Mw2_Mw_mh_cache[NumPar][i] );
    } else {
        gslpp::complex newResult = PV.B0(Mz*Mz, Mw*Mw, Mw*Mw, mh*mh);
        CacheShift(B0_Mz_Mw2_Mw_mh_cache, NumPar, params, newResult);
        return newResult;
    } 
}

gslpp::complex THDMcache::B0_Mz_0_Mw_mh(const double Mz, const double Mw, const double mh) const {
    int NumPar = 3;
    double params[] = {Mz, Mw, mh};

    int i = CacheCheck(B0_Mz_0_Mw_mh_cache, NumPar, params);
    if (i>=0) {
        return ( B0_Mz_0_Mw_mh_cache[NumPar][i] );
    } else {
        gslpp::complex newResult = PV.B0(Mz*Mz, 0., Mw*Mw, mh*mh);
        CacheShift(B0_Mz_0_Mw_mh_cache, NumPar, params, newResult);
        return newResult;
    } 
}

gslpp::complex THDMcache::B00_Mz_Mz2_mH_mA(const double Mz, const double mH, const double mA) const {
    int NumPar = 3;
    double params[] = {Mz, mH, mA};

    int i = CacheCheck(B00_Mz_Mz2_mH_mA_cache, NumPar, params);
    if (i>=0) {
        return ( B00_Mz_Mz2_mH_mA_cache[NumPar][i] );
    } else {
        gslpp::complex newResult = PV.B00(Mz*Mz, Mz*Mz, mH*mH, mA*mA);
        CacheShift(B00_Mz_Mz2_mH_mA_cache, NumPar, params, newResult);
        return newResult;
    } 
}

gslpp::complex THDMcache::B00_Mz_0_mH_mA(const double Mz, const double mH, const double mA) const {
    int NumPar = 3;
    double params[] = {Mz, mH, mA};

    int i = CacheCheck(B00_Mz_0_mH_mA_cache, NumPar, params);
    if (i>=0) {
        return ( B00_Mz_0_mH_mA_cache[NumPar][i] );
    } else {
        gslpp::complex newResult = PV.B00(Mz*Mz, 0., mH*mH, mA*mA);
        CacheShift(B00_Mz_0_mH_mA_cache, NumPar, params, newResult);
        return newResult;
    } 
}

gslpp::complex THDMcache::B00_Mz_Mz2_mHp_mHp(const double Mz, const double mHp) const {
    int NumPar = 2;
    double params[] = {Mz, mHp};

    int i = CacheCheck(B00_Mz_Mz2_mHp_mHp_cache, NumPar, params);
    if (i>=0) {
        return ( B00_Mz_Mz2_mHp_mHp_cache[NumPar][i] );
    } else {
        gslpp::complex newResult = PV.B00(Mz*Mz, Mz*Mz, mHp*mHp, mHp*mHp);
        CacheShift(B00_Mz_Mz2_mHp_mHp_cache, NumPar, params, newResult);
        return newResult;
    } 
}

gslpp::complex THDMcache::B00_Mz_0_mHp_mHp(const double Mz, const double mHp) const {
    int NumPar = 2;
    double params[] = {Mz, mHp};

    int i = CacheCheck(B00_Mz_0_mHp_mHp_cache, NumPar, params);
    if (i>=0) {
        return ( B00_Mz_0_mHp_mHp_cache[NumPar][i] );
    } else {
        gslpp::complex newResult = PV.B00(Mz*Mz, 0., mHp*mHp, mHp*mHp);
        CacheShift(B00_Mz_0_mHp_mHp_cache, NumPar, params, newResult);
        return newResult;
    } 
}

gslpp::complex THDMcache::B00_Mz_Mz2_mh_mA(const double Mz, const double mh, const double mA) const {
    int NumPar = 3;
    double params[] = {Mz, mh, mA};

    int i = CacheCheck(B00_Mz_Mz2_mh_mA_cache, NumPar, params);
    if (i>=0) {
        return ( B00_Mz_Mz2_mh_mA_cache[NumPar][i] );
    } else {
        gslpp::complex newResult = PV.B00(Mz*Mz, Mz*Mz, mh*mh, mA*mA);
        CacheShift(B00_Mz_Mz2_mh_mA_cache, NumPar, params, newResult);
        return newResult;
    } 
}

gslpp::complex THDMcache::B00_Mz_0_mh_mA(const double Mz, const double mh, const double mA) const {
    int NumPar = 3;
    double params[] = {Mz, mh, mA};

    int i = CacheCheck(B00_Mz_0_mh_mA_cache, NumPar, params);
    if (i>=0) {
        return ( B00_Mz_0_mh_mA_cache[NumPar][i] );
    } else {
        gslpp::complex newResult = PV.B00(Mz*Mz, 0., mh*mh, mA*mA);
        CacheShift(B00_Mz_0_mh_mA_cache, NumPar, params, newResult);
        return newResult;
    } 
}

gslpp::complex THDMcache::B00_Mz_Mz2_Mz_mH(const double Mz, const double mH) const {
    int NumPar = 2;
    double params[] = {Mz, mH};

    int i = CacheCheck(B00_Mz_Mz2_Mz_mH_cache, NumPar, params);
    if (i>=0) {
        return ( B00_Mz_Mz2_Mz_mH_cache[NumPar][i] );
    } else {
        gslpp::complex newResult = PV.B00(Mz*Mz, Mz*Mz, Mz*Mz, mH*mH);
        CacheShift(B00_Mz_Mz2_Mz_mH_cache, NumPar, params, newResult);
        return newResult;
    } 
}

gslpp::complex THDMcache::B00_Mz_0_Mz_mH(const double Mz, const double mH) const {
    int NumPar = 2;
    double params[] = {Mz, mH};

    int i = CacheCheck(B00_Mz_0_Mz_mH_cache, NumPar, params);
    if (i>=0) {
        return ( B00_Mz_0_Mz_mH_cache[NumPar][i] );
    } else {
        gslpp::complex newResult = PV.B00(Mz*Mz, 0., Mz*Mz, mH*mH);
        CacheShift(B00_Mz_0_Mz_mH_cache, NumPar, params, newResult);
        return newResult;
    } 
}

gslpp::complex THDMcache::B00_Mz_Mz2_Mz_mh(const double Mz, const double mh) const {
    int NumPar = 2;
    double params[] = {Mz, mh};

    int i = CacheCheck(B00_Mz_Mz2_Mz_mh_cache, NumPar, params);
    if (i>=0) {
        return ( B00_Mz_Mz2_Mz_mh_cache[NumPar][i] );
    } else {
        gslpp::complex newResult = PV.B00(Mz*Mz, Mz*Mz, Mz*Mz, mh*mh);
        CacheShift(B00_Mz_Mz2_Mz_mh_cache, NumPar, params, newResult);
        return newResult;
    } 
}

gslpp::complex THDMcache::B00_Mz_0_Mz_mh(const double Mz, const double mh) const {
    int NumPar = 2;
    double params[] = {Mz, mh};

    int i = CacheCheck(B00_Mz_0_Mz_mh_cache, NumPar, params);
    if (i>=0) {
        return ( B00_Mz_0_Mz_mh_cache[NumPar][i] );
    } else {
        gslpp::complex newResult = PV.B00(Mz*Mz, 0., Mz*Mz, mh*mh);
        CacheShift(B00_Mz_0_Mz_mh_cache, NumPar, params, newResult);
        return newResult;
    } 
}

gslpp::complex THDMcache::B00_Mz_Mw2_mA_mHp(const double Mz, const double Mw, const double mA, const double mHp) const {
    int NumPar = 4;
    double params[] = {Mz,Mw, mA, mHp};

    int i = CacheCheck(B00_Mz_Mw2_mA_mHp_cache, NumPar, params);
    if (i>=0) {
        return ( B00_Mz_Mw2_mA_mHp_cache[NumPar][i] );
    } else {
        gslpp::complex newResult = PV.B00(Mz*Mz, Mw*Mw, mA*mA, mHp*mHp);
        CacheShift(B00_Mz_Mw2_mA_mHp_cache, NumPar, params, newResult);
        return newResult;
    } 
}

gslpp::complex THDMcache::B00_Mz_0_mA_mHp(const double Mz, const double mA, const double mHp) const {
    int NumPar = 3;
    double params[] = {Mz, mA, mHp};

    int i = CacheCheck(B00_Mz_0_mA_mHp_cache, NumPar, params);
    if (i>=0) {
        return ( B00_Mz_0_mA_mHp_cache[NumPar][i] );
    } else {
        gslpp::complex newResult = PV.B00(Mz*Mz, 0., mA*mA, mHp*mHp);
        CacheShift(B00_Mz_0_mA_mHp_cache, NumPar, params, newResult);
        return newResult;
    } 
}

gslpp::complex THDMcache::B00_Mz_Mw2_mHp_mHp(const double Mz, const double Mw, const double mHp) const {
    int NumPar = 3;
    double params[] = {Mz,Mw, mHp};

    int i = CacheCheck(B00_Mz_Mw2_mHp_mHp_cache, NumPar, params);
    if (i>=0) {
        return ( B00_Mz_Mw2_mHp_mHp_cache[NumPar][i] );
    } else {
        gslpp::complex newResult = PV.B00(Mz*Mz, Mw*Mw, mHp*mHp, mHp*mHp);
        CacheShift(B00_Mz_Mw2_mHp_mHp_cache, NumPar, params, newResult);
        return newResult;
    } 
}

gslpp::complex THDMcache::B00_Mz_Mw2_Mw_mH(const double Mz, const double Mw, const double mH) const {
    int NumPar = 3;
    double params[] = {Mz,Mw, mH};

    int i = CacheCheck(B00_Mz_Mw2_Mw_mH_cache, NumPar, params);
    if (i>=0) {
        return ( B00_Mz_Mw2_Mw_mH_cache[NumPar][i] );
    } else {
        gslpp::complex newResult = PV.B00(Mz*Mz, Mw*Mw, Mw*Mw, mH*mH);
        CacheShift(B00_Mz_Mw2_Mw_mH_cache, NumPar, params, newResult);
        return newResult;
    } 
}

gslpp::complex THDMcache::B00_Mz_0_Mw_mH(const double Mz, const double Mw, const double mH) const {
    int NumPar = 3;
    double params[] = {Mz,Mw, mH};

    int i = CacheCheck(B00_Mz_0_Mw_mH_cache, NumPar, params);
    if (i>=0) {
        return ( B00_Mz_0_Mw_mH_cache[NumPar][i] );
    } else {
        gslpp::complex newResult = PV.B00(Mz*Mz, Mw*Mw, Mw*Mw, mH*mH);
        CacheShift(B00_Mz_0_Mw_mH_cache, NumPar, params, newResult);
        return newResult;
    } 
}


gslpp::complex THDMcache::B00_Mz_Mw2_Mw_mh(const double Mz, const double Mw, const double mh) const {
    int NumPar = 3;
    double params[] = {Mz,Mw, mh};

    int i = CacheCheck(B00_Mz_Mw2_Mw_mh_cache, NumPar, params);
    if (i>=0) {
        return ( B00_Mz_Mw2_Mw_mh_cache[NumPar][i] );
    } else {
        gslpp::complex newResult = PV.B00(Mz*Mz, Mw*Mw, Mw*Mw, mh*mh);
        CacheShift(B00_Mz_Mw2_Mw_mh_cache, NumPar, params, newResult);
        return newResult;
    } 
}

gslpp::complex THDMcache::B00_Mz_0_Mw_mh(const double Mz, const double Mw, const double mh) const {
    int NumPar = 3;
    double params[] = {Mz,Mw, mh};

    int i = CacheCheck(B00_Mz_0_Mw_mh_cache, NumPar, params);
    if (i>=0) {
        return ( B00_Mz_0_Mw_mh_cache[NumPar][i] );
    } else {
        gslpp::complex newResult = PV.B00(Mz*Mz, 0., Mw*Mw, mh*mh);
        CacheShift(B00_Mz_0_Mw_mh_cache, NumPar, params, newResult);
        return newResult;
    } 
}

gslpp::complex THDMcache::B00_Mz_Mw2_mH_mHp(const double Mz, const double Mw, const double mH, const double mHp) const {
    int NumPar = 4;
    double params[] = {Mz,Mw, mH, mHp};

    int i = CacheCheck(B00_Mz_Mw2_mH_mHp_cache, NumPar, params);
    if (i>=0) {
        return ( B00_Mz_Mw2_mH_mHp_cache[NumPar][i] );
    } else {
        gslpp::complex newResult = PV.B00(Mz*Mz, Mw*Mw, mH*mH, mHp*mHp);
        CacheShift(B00_Mz_Mw2_mH_mHp_cache, NumPar, params, newResult);
        return newResult;
    } 
}

gslpp::complex THDMcache::B00_Mz_0_mH_mHp(const double Mz, const double mH, const double mHp) const {
    int NumPar = 3;
    double params[] = {Mz, mH, mHp};

    int i = CacheCheck(B00_Mz_0_mH_mHp_cache, NumPar, params);
    if (i>=0) {
        return ( B00_Mz_0_mH_mHp_cache[NumPar][i] );
    } else {
        gslpp::complex newResult = PV.B00(Mz*Mz, 0., mH*mH, mHp*mHp);
        CacheShift(B00_Mz_0_mH_mHp_cache, NumPar, params, newResult);
        return newResult;
    } 
}

gslpp::complex THDMcache::B00_Mz_Mw2_mh_mHp(const double Mz, const double Mw, const double mh, const double mHp) const {
    int NumPar = 4;
    double params[] = {Mz,Mw, mh, mHp};

    int i = CacheCheck(B00_Mz_Mw2_mh_mHp_cache, NumPar, params);
    if (i>=0) {
        return ( B00_Mz_Mw2_mh_mHp_cache[NumPar][i] );
    } else {
        gslpp::complex newResult = PV.B00(Mz*Mz, Mw*Mw, mh*mh, mHp*mHp);
        CacheShift(B00_Mz_Mw2_mh_mHp_cache, NumPar, params, newResult);
        return newResult;
    } 
}

gslpp::complex THDMcache::B00_Mz_0_mh_mHp(const double Mz, const double mh, const double mHp) const {
    int NumPar = 3;
    double params[] = {Mz, mh, mHp};

    int i = CacheCheck(B00_Mz_0_mh_mHp_cache, NumPar, params);
    if (i>=0) {
        return ( B00_Mz_0_mh_mHp_cache[NumPar][i] );
    } else {
        gslpp::complex newResult = PV.B00(Mz*Mz, 0., mh*mh, mHp*mHp);
        CacheShift(B00_Mz_0_mh_mHp_cache, NumPar, params, newResult);
        return newResult;
    } 
}


//HEAVY


void THDMcache::read(){
    

    std::string tablepath="/Users/Roma1/Desktop/HEPfit/THDM/tabs/";
    std::stringstream br1,br2,br3,br4,br5,br6,br7;
    std::stringstream pc1,pc2,pc3,pc4,pc5;
    std::stringstream dw1;
    std::stringstream cs1,cs2,cs3,cs4,cs5,cs6;
    std::stringstream ex1,ex2,ex3,ex4,ex5,ex6,ex7,ex8,ex9,ex10,ex11,ex12,ex13,ex14,ex15,ex16;

    std::cout<<"reading tables"<<std::endl;

    br1 << tablepath << "GridSM1.dat";
    array1 = readTable(br1.str(),19861);
    br2 << tablepath << "GridSM2.dat";
    array2 = readTable(br2.str(),19861);
    br3 << tablepath << "GridSM3.dat";
    array3 = readTable(br3.str(),19861); 
    br4 << tablepath << "GridSM4.dat";
    array4 = readTable(br4.str(),19861);
    br5 << tablepath << "GridSM5.dat";
    array5 = readTable(br5.str(),19861);
    br6 << tablepath << "GridSM6.dat";
    array6 = readTable(br6.str(),19861);
    br7 << tablepath << "GridSM7.dat";
    array7 = readTable(br7.str(),19861);
    pc1 << tablepath << "GridSM11.dat";
    array11 = readTable(pc1.str(),1971);
    pc2 << tablepath << "GridSM12.dat";
    array12 = readTable(pc2.str(),1971);
    pc3 << tablepath << "GridSM13.dat";
    array13 = readTable(pc3.str(),1971);
    pc4 << tablepath << "GridSM14.dat";
    array14 = readTable(pc4.str(),1971);
    pc5 << tablepath << "GridSM15.dat";
    array15 = readTable(pc5.str(),1971);
    ex1 << tablepath << "GridSM16.dat";
    array16 = readTable(ex1.str(),9851);
    ex2 << tablepath << "GridSM17.dat";
    array17 = readTable(ex2.str(),9851);
    dw1 << tablepath << "GridSM18.dat";
    array18 = readTable(dw1.str(),19861);
    cs1 << tablepath << "GridSM19.dat";
    array19 = readTable(cs1.str(),186);
    cs2 << tablepath << "GridSM20.dat";
    array20 = readTable(cs2.str(),186);
    cs3 << tablepath << "GridSM21.dat";
    array21 = readTable(cs3.str(),186);
    cs4 << tablepath << "GridSM22.dat";
    array22 = readTable(cs4.str(),186);
    cs5 << tablepath << "GridSM26.dat";
    array26 = readTable(cs5.str(),992);
    cs6 << tablepath << "GridSM27.dat";
    array27 = readTable(cs6.str(),185);
    ex3 << tablepath << "GridSM29.dat";
    array29 = readTable(ex3.str(),986);
    ex4 << tablepath << "GridSM31.dat";
    array31 = readTable(ex4.str(),985);
    ex5 << tablepath << "GridSM32.dat";
    array32 = readTable(ex5.str(),496);
    ex6 << tablepath << "GridSM33.dat";
    array33 = readTable(ex6.str(),496);
    ex7 << tablepath << "GridSM34.dat";
    array34 = readTable(ex7.str(),991);
    ex8 << tablepath << "GridSM35.dat";
    array35 = readTable(ex8.str(),496);
    ex9 << tablepath << "GridSM36.dat";
    array36 = readTable(ex9.str(),496);
    ex10 << tablepath << "GridSM37.dat";
    array37 = readTable(ex10.str(),100);
    ex11 << tablepath << "GridSM38.dat";
    array38 = readTable(ex11.str(),100);
    ex12 << tablepath << "GridSM39.dat";
    array39 = readTable(ex12.str(),986);
    ex13 << tablepath << "GridSM44.dat";
    array44 = readTable(ex13.str(),986);
    ex14 << tablepath << "GridSM45.dat";
    array45 = readTable(ex14.str(),986);
    ex15 << tablepath << "Grid_X_bb.dat";
    arrayX_bb = readTable(ex15.str(),199);
    ex16 << tablepath << "Grid_X_tt.dat";
    arrayX_tt = readTable(ex16.str(),198);
}    
    
    
    


double THDMcache::Br_HPtott(double mass){
    int NumPar = 1;
    double params[] = {mass};
    return pow(10.0,interpolate(array1,mass));  
}



double THDMcache::Br_HPtobb(double mass){
    int NumPar = 1;
    double params[] = {mass};
    return pow(10.0,interpolate(array2,mass));  
}



double THDMcache::Br_HPtotautau(double mass){
    int NumPar = 1;
    double params[] = {mass};
    return pow(10.0,interpolate(array3,mass));  
}



double THDMcache::Br_HPtocc(double mass){
    int NumPar = 1;
    double params[] = {mass};
    return pow(10.0,interpolate(array4,mass));  
}



double THDMcache::Br_HPtomumu(double mass){
    int NumPar = 1;
    double params[] = {mass};
    return pow(10.0,interpolate(array5,mass));  
}



double THDMcache::Br_HPtoZZ(double mass){
    int NumPar = 1;
    double params[] = {mass};
    return pow(10.0,interpolate(array6,mass));  
}



double THDMcache::Br_HPtoWW(double mass){
    int NumPar = 1;
    double params[] = {mass};
    return pow(10.0,interpolate(array7,mass));  
}



double THDMcache::pc_ggFtoHP(double mass){
    int NumPar = 1;
    double params[] = {mass};
    return interpolate(array11,mass);  
}



double THDMcache::pc_VBFtoHP(double mass){
    int NumPar = 1;
    double params[] = {mass};
    return interpolate(array12,mass);  
}



double THDMcache::pc_WHP_HP(double mass){
    int NumPar = 1;
    double params[] = {mass};
    return interpolate(array13,mass);  
}


double THDMcache::pc_ZHP_HP(double mass){
    int NumPar = 1;
    double params[] = {mass};
    return interpolate(array14,mass);  
}



double THDMcache::pc_ttFtoHP(double mass){
    int NumPar = 1;
    double params[] = {mass};
    return interpolate(array15,mass);  
}



double THDMcache::ex_HP_ZZ(double mass){
    int NumPar = 1;
    double params[] = {mass};
    return interpolate(array16,mass);  
}



double THDMcache::ex_A_tautau(double mass){
    int NumPar = 1;
    double params[] = {mass};
    return interpolate(array17,mass);  
}



double THDMcache::GammaHPtotSM(double mass){
    int NumPar = 1;
    double params[] = {mass};
    return pow(10.0,interpolate(array18,mass));  
}



double THDMcache::cs_ggFtoHP(double mass){
    int NumPar = 1;
    double params[] = {mass}; 
    return pow(10.0,interpolate (array19,log10(mass)));  
}


double THDMcache::cs_ggHP_tt(double mass){
    int NumPar = 1;
    double params[] = {mass};
    return pow(10.0,interpolate (array20,log10(mass)));  
}



double THDMcache::cs_ggHP_bb(double mass){
    int NumPar = 1;
    double params[] = {mass};
    return pow(10.0,interpolate (array21,log10(mass)));  
}



double THDMcache::cs_ggFtoA(double mass){
    int NumPar = 1;
    double params[] = {mass};
    return pow(10.0,interpolate (array22,log10(mass)));  
}



double THDMcache::cs_ggFtoHPbbtotautaubb(double mass){
    int NumPar = 1;
    double params[] = {mass};
    return interpolate(array26,mass);  
}



double THDMcache::cs_bbFtoHP(double mass){
    int NumPar = 1;
    double params[] = {mass};
    return pow(10.0,interpolate (array27,log10(mass)));  
}



double THDMcache::ex_ggF_A_hZ(double mass){
    int NumPar = 1;
    double params[] = {mass};
    return interpolate(array29,mass); 
}


double THDMcache::ex_H_WW(double mass){
    int NumPar = 1;
    double params[] = {mass};
    return interpolate(array31,mass); 
}



double THDMcache::ex_H_hh_gagabb(double mass){
    int NumPar = 1;
    double params[] = {mass};
    return interpolate(array32,mass);  
}



double THDMcache::ex_H_hh_bbbb(double mass){
    int NumPar = 1;
    double params[] = {mass};
    return interpolate(array33,mass);  
}



double THDMcache::ex_H_gaga(double mass){
    int NumPar = 1;
    double params[] = {mass};
    return interpolate(array34,mass);  
}



double THDMcache::ex_ggF_H_tautau(double mass){
    int NumPar = 1;
    double params[] = {mass};
    return interpolate(array35,mass);    
}



double THDMcache::ex_bbF_H_tautau(double mass){
    int NumPar = 1;
    double params[] = {mass};
    return interpolate(array36,mass);  
}



double THDMcache::ex_ggF_H_WW(double mass){
    int NumPar = 1;
    double params[] = {mass};
    return interpolate(array37,mass);  
}



double THDMcache::ex_VBF_H_WW(double mass){
    int NumPar = 1;
    double params[] = {mass};
    return interpolate(array38,mass);  
}



double THDMcache::ex_H_hh(double mass){
    int NumPar = 1;
    double params[] = {mass}; 
    return interpolate(array39,mass);  
}



double THDMcache::ex_ggF_A_hZ_tautauZ(double mass){
    int NumPar = 1;
    double params[] = {mass};
    return interpolate(array44,mass);  
}



double THDMcache::ex_ggF_A_hZ_bbZ(double mass){
    int NumPar = 1;
    double params[] = {mass}; 
    return interpolate(array45,mass);
}



double THDMcache::ex_H_bb(double mass){
    int NumPar = 1;
    double params[] = {mass}; 
    return interpolate (arrayX_bb,mass);  
}



double THDMcache::ex_H_tt(double mass){
    int NumPar = 1;
    double params[] = {mass};
    return interpolate (arrayX_tt,mass);  
}



gslpp::matrix<double> THDMcache::readTable(std::string filename, int rowN){

    std::ifstream INfile;
    std::string lineTab;
    INfile.open( filename.c_str() );
    if(INfile.fail()){
        std::cout<<"error: in HEAVY"<<" doesn't exist!"<<std::endl;
    }

    gslpp::matrix<double> arrayTab(rowN, 2, 0.);
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

//Old Interpolate function

//double THDMcache::OLDinterpolate (gslpp::matrix<double> arrayTab, double mass){
//
//    int rowN=arrayTab.size_i();
//
//    double y = 0.0;
//    if(mass<arrayTab(0,0)){
//        std::cout<<"error: your mass value is smaller than the minimum mass value"<<std::endl;
//    }
//    else if(mass>arrayTab(rowN-1,0)){
//        std::cout<<"error: your mass value is greater than the maximum mass value"<<std::endl;
//    }
//    else{
//        for(int i=0; i<rowN; i++){
//            if( (mass>=arrayTab(i,0)) && (mass<=arrayTab(i+1,0)) ){
//                y =(arrayTab(i+1,1)-arrayTab(i,1))/(arrayTab(i+1,0)
//                   -arrayTab(i,0))*(mass-arrayTab(i,0))+arrayTab(i,1);
//                 break;       
//        }
//            
//        }
//        return y;
//    }
//    
//}

//New Interpolate function

double THDMcache::interpolate (gslpp::matrix<double> arrayTab, double mass){

    int rowN=arrayTab.size_i();
    
    double Mmin = arrayTab(0,0);
    double Mmax = arrayTab(rowN-1,0);
    double interval = arrayTab(1,0)-arrayTab(0,0);
    int Nintervals = (mass-Mmin)/interval;
    double y = 0.0;
       
    if(mass<Mmin){
        std::cout<<"error: your mass value is smaller than the minimum mass value"<<std::endl;
    }
    else if(mass>Mmax){
        std::cout<<"error: your mass value is greater than the maximum mass value"<<std::endl;
    }
    else{
        
        y =(arrayTab(Nintervals+1,1)-arrayTab(Nintervals,1))/(arrayTab(Nintervals+1,0)
                   -arrayTab(Nintervals,0))*(mass-arrayTab(Nintervals,0))+arrayTab(Nintervals,1);
        
        return y;
    }
}