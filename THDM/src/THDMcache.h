/* 
 * Copyright (C) 2012 SusyFit Collaboration
 * All rights reserved.
 *
 * For the licensing terms see doc/COPYING.
 */

#ifndef THDMCACHE_H
#define	THDMCACHE_H

#include <cmath>
#include <PVfunctions.h>
//#include "THDM.h"

#include <stdexcept>
#include <gslpp.h>


class THDMcache {
    
public:
    
    /**
     * @brief THDMcache constructor
     * @param[in] SM_i reference to a StandardModel object
     * @param[in] THDM_i reference to a THDM object
     */
    THDMcache();
    
    static const int CacheSize = 5;

    
    int CacheCheck(const gslpp::complex cache[][CacheSize], 
                   const int NumPar, const double params[]) const;
    
    void CacheShift(gslpp::complex cache[][CacheSize], const int NumPar, 
                    const double params[], const gslpp::complex newResult) const; 

////////////////////////////////////////////////////////////////////////////////
    
    /**
     * @return an object of PVfunctions class
     */
    const PVfunctions getPV() const {
        return PV;
    }
    
    /**
     * @return a reference to the THDM object in THDMcache class
     */
//   const StandardModel& GetTHDM() const {
//        return myTHDM;
//    }

    
    gslpp::matrix<double> readTable(std::string filename, int rowN);
//    double OLDinterpolate (gslpp::matrix<double> arrayTab, double mass);
    double interpolate (gslpp::matrix<double> arrayTab, double mass);
    
    void read();
    
    gslpp::matrix<double> array1;
    gslpp::matrix<double> array2;
    gslpp::matrix<double> array3;
    gslpp::matrix<double> array4;
    gslpp::matrix<double> array5;
    gslpp::matrix<double> array6;
    gslpp::matrix<double> array7;
    gslpp::matrix<double> array11;
    gslpp::matrix<double> array12;
    gslpp::matrix<double> array13;
    gslpp::matrix<double> array14;
    gslpp::matrix<double> array15;
    gslpp::matrix<double> array16;
    gslpp::matrix<double> array17;
    gslpp::matrix<double> array18;
    gslpp::matrix<double> array19;
    gslpp::matrix<double> array20;
    gslpp::matrix<double> array21;
    gslpp::matrix<double> array22;
    gslpp::matrix<double> array26;
    gslpp::matrix<double> array27;
    gslpp::matrix<double> array29;
    gslpp::matrix<double> array31;
    gslpp::matrix<double> array32;
    gslpp::matrix<double> array33;
    gslpp::matrix<double> array34;
    gslpp::matrix<double> array35;
    gslpp::matrix<double> array36;
    gslpp::matrix<double> array37;
    gslpp::matrix<double> array38;
    gslpp::matrix<double> array39;
    gslpp::matrix<double> array44;
    gslpp::matrix<double> array45;
    gslpp::matrix<double> arrayX_bb;
    gslpp::matrix<double> arrayX_tt;
    
    double Br_HPtott(double mass);
    double Br_HPtobb(double mass);
    double Br_HPtotautau(double mass);
    double Br_HPtocc(double mass);
    double Br_HPtomumu(double mass);
    double Br_HPtoZZ(double mass);
    double Br_HPtoWW(double mass);
    double pc_ggFtoHP(double mass);
    double pc_VBFtoHP(double mass);
    double pc_WHP_HP(double mass);
    double pc_ZHP_HP(double mass);
    double pc_ttFtoHP(double mass);
    double ex_HP_ZZ(double mass);
    double ex_A_tautau(double mass);
    double GammaHPtotSM(double mass);
    double cs_ggFtoHP(double mass);
    double cs_ggHP_tt(double mass);
    double cs_ggHP_bb(double mass);
    double cs_ggFtoA(double mass);
    double cs_ggFtoHPbbtotautaubb(double mass);
    double cs_bbFtoHP(double mass);
    double ex_ggF_A_hZ(double mass);
    double ex_H_WW(double mass);
    double ex_H_hh_gagabb(double mass);
    double ex_H_hh_bbbb(double mass);
    double ex_H_gaga(double mass);
    double ex_ggF_H_tautau(double mass);
    double ex_bbF_H_tautau(double mass);
    double ex_ggF_H_WW(double mass);
    double ex_VBF_H_WW(double mass);
    double ex_H_hh(double mass);
    double ex_ggF_A_hZ_tautauZ(double mass);
    double ex_ggF_A_hZ_bbZ(double mass);
    double ex_H_bb(double mass);
    double ex_H_tt(double mass);
    
    /*One-loop functions*/
    
    gslpp::complex B0_Mz_Mz2_Mz_mH(const double Mz, const double mH) const;
    gslpp::complex B0_Mz_0_Mz_mH(const double Mz, const double mH) const;
    gslpp::complex B0_Mz_Mz2_Mz_mh(const double Mz, const double mh) const;
    gslpp::complex B0_Mz_0_Mz_mh(const double Mz, const double mh) const;
    gslpp::complex B0_Mz_Mw2_Mw_mH(const double Mz, const double Mw, const double mH) const;
    gslpp::complex B0_Mz_0_Mw_mH(const double Mz, const double Mw, const double mH) const;
    gslpp::complex B0_Mz_Mw2_Mw_mh(const double Mz, const double Mw, const double mh) const;
    gslpp::complex B0_Mz_0_Mw_mh(const double Mz, const double Mw, const double mh) const;
    
    gslpp::complex B00_Mz_Mz2_mH_mA(const double Mz, const double mH, const double mA) const;
    gslpp::complex B00_Mz_0_mH_mA(const double Mz, const double mH, const double mA) const;
    gslpp::complex B00_Mz_Mz2_mHp_mHp(const double Mz, const double mHp) const;
    gslpp::complex B00_Mz_0_mHp_mHp(const double Mz, const double mHp) const;
    gslpp::complex B00_Mz_Mz2_mh_mA(const double Mz, const double mh, const double mA) const;
    gslpp::complex B00_Mz_0_mh_mA(const double Mz, const double mh, const double mA) const;
    gslpp::complex B00_Mz_Mz2_Mz_mH(const double Mz, const double mH) const;
    gslpp::complex B00_Mz_0_Mz_mH(const double Mz, const double mH) const;
    gslpp::complex B00_Mz_Mz2_Mz_mh(const double Mz, const double mh) const;
    gslpp::complex B00_Mz_0_Mz_mh(const double Mz, const double mh) const;
    gslpp::complex B00_Mz_Mw2_mA_mHp(const double Mz, const double Mw, const double mA, const double mHp) const;
    gslpp::complex B00_Mz_0_mA_mHp(const double Mz, const double mA, const double mHp) const;
    gslpp::complex B00_Mz_Mw2_mHp_mHp(const double Mz, const double Mw, const double mHp) const;
    gslpp::complex B00_Mz_Mw2_Mw_mH(const double Mz, const double Mw, const double mH) const;
    gslpp::complex B00_Mz_0_Mw_mH(const double Mz, const double Mw, const double mH) const;
    gslpp::complex B00_Mz_Mw2_Mw_mh(const double Mz, const double Mw, const double mh) const;
    gslpp::complex B00_Mz_0_Mw_mh(const double Mz, const double Mw, const double mh) const;
    gslpp::complex B00_Mz_Mw2_mH_mHp(const double Mz, const double Mw, const double mH, const double mHp) const;
    gslpp::complex B00_Mz_0_mH_mHp(const double Mz, const double mH, const double mHp) const;
    gslpp::complex B00_Mz_Mw2_mh_mHp(const double Mz, const double Mw, const double mh, const double mHp) const;
    gslpp::complex B00_Mz_0_mh_mHp(const double Mz, const double mh, const double mHp) const;
    
 ///////////////////////////////////////////////////////////////////////////////   
    
private:

    const PVfunctions PV;
    //const THDM myTHDM;
    
    
    ////////////////////////////////////////////////////////////////////////////
    //Caches
    
    /*One-loop functions*/
    mutable gslpp::complex B0_Mz_Mz2_Mz_mH_cache[3][CacheSize];
    mutable gslpp::complex B0_Mz_0_Mz_mH_cache[3][CacheSize];
    mutable gslpp::complex B0_Mz_Mz2_Mz_mh_cache[3][CacheSize];
    mutable gslpp::complex B0_Mz_0_Mz_mh_cache[3][CacheSize];
    mutable gslpp::complex B0_Mz_Mw2_Mw_mH_cache[4][CacheSize];
    mutable gslpp::complex B0_Mz_0_Mw_mH_cache[4][CacheSize];
    mutable gslpp::complex B0_Mz_Mw2_Mw_mh_cache[4][CacheSize];
    mutable gslpp::complex B0_Mz_0_Mw_mh_cache[4][CacheSize];
    
    mutable gslpp::complex B00_Mz_Mz2_mH_mA_cache[4][CacheSize];
    mutable gslpp::complex B00_Mz_0_mH_mA_cache[4][CacheSize];
    mutable gslpp::complex B00_Mz_Mz2_mHp_mHp_cache[3][CacheSize];
    mutable gslpp::complex B00_Mz_0_mHp_mHp_cache[3][CacheSize];
    mutable gslpp::complex B00_Mz_Mz2_mh_mA_cache[4][CacheSize];
    mutable gslpp::complex B00_Mz_0_mh_mA_cache[4][CacheSize];
    mutable gslpp::complex B00_Mz_Mz2_Mz_mH_cache[3][CacheSize];
    mutable gslpp::complex B00_Mz_0_Mz_mH_cache[3][CacheSize];
    mutable gslpp::complex B00_Mz_Mz2_Mz_mh_cache[3][CacheSize];
    mutable gslpp::complex B00_Mz_0_Mz_mh_cache[3][CacheSize];
    mutable gslpp::complex B00_Mz_Mw2_mA_mHp_cache[5][CacheSize];
    mutable gslpp::complex B00_Mz_0_mA_mHp_cache[4][CacheSize];
    mutable gslpp::complex B00_Mz_Mw2_mHp_mHp_cache[4][CacheSize];
    mutable gslpp::complex B00_Mz_Mw2_Mw_mH_cache[4][CacheSize];
    mutable gslpp::complex B00_Mz_0_Mw_mH_cache[4][CacheSize];
    mutable gslpp::complex B00_Mz_Mw2_Mw_mh_cache[4][CacheSize];
    mutable gslpp::complex B00_Mz_0_Mw_mh_cache[4][CacheSize];
    mutable gslpp::complex B00_Mz_Mw2_mH_mHp_cache[5][CacheSize];
    mutable gslpp::complex B00_Mz_0_mH_mHp_cache[4][CacheSize];
    mutable gslpp::complex B00_Mz_Mw2_mh_mHp_cache[5][CacheSize];
    mutable gslpp::complex B00_Mz_0_mh_mHp_cache[4][CacheSize];

};

#endif	/* THDMCACHE_H */
