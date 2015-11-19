/* 
 * Copyright (C) 2012 HEPfit Collaboration
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
    
    int CacheCheckReal(const double cache[][CacheSize],
                   const int NumPar, const double params[]) const;
    
    void CacheShift(gslpp::complex cache[][CacheSize], const int NumPar,
                    const double params[], const gslpp::complex newResult) const; 

    void CacheShiftReal(double cache[][CacheSize], const int NumPar,
                    const double params[], const double newResult) const; 

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

    
    gslpp::matrix<double> readTable(std::string filename, int rowN, int colN);

    double interpolate (gslpp::matrix<double> arrayTab, double mass);
    double interpolate2D (gslpp::matrix<double> arrayTab, double x, double y);
    
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
    gslpp::matrix<double> array18;
    gslpp::matrix<double> array19;
    gslpp::matrix<double> array20;
    gslpp::matrix<double> array21;
    gslpp::matrix<double> array22;
    gslpp::matrix<double> array27;
    gslpp::matrix<double> array29;
    gslpp::matrix<double> array32;
    gslpp::matrix<double> array33;
    gslpp::matrix<double> array34;
    gslpp::matrix<double> array35;
    gslpp::matrix<double> array36;
    gslpp::matrix<double> array44;
    gslpp::matrix<double> array45;
    gslpp::matrix<double> bbF_phi_bb_CMS;
    gslpp::matrix<double> pp_phi_tt_ATLAS;
    gslpp::matrix<double> ggF_phi_tautau_CMS;
    gslpp::matrix<double> bbF_phi_tautau_CMS;
    gslpp::matrix<double> ggF_phi_gaga_CMS;
    gslpp::matrix<double> ggF_H_WW_ATLAS;
    gslpp::matrix<double> VBF_H_WW_ATLAS;
    gslpp::matrix<double> ggF_H_hh_ATLAS;
    gslpp::matrix<double> ggF_H_hh_bbtautau_CMS;
    gslpp::matrix<double> ggF_A_hZ_tautaull_CMS;
    gslpp::matrix<double> arraybsgamma;

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

    double GammaHPtotSM(double mass);

    double cs_ggFtoHP(double mass);
    double cs_ggHP_tt(double mass);
    double cs_ggHP_bb(double mass);
    double cs_ggFtoA(double mass);
    double cs_bbFtoHP(double mass);

    double ex_pp_H_ZZ_CMS(double mass);            //
    double ex_ggF_A_hZ_bbll_CMS(double mass);      //
    double ex_pp_phi_hh_gagabb_CMS(double mass);   //
    double ex_pp_phi_hh_bbbb_CMS(double mass);     //
    double ex_ggF_phi_gaga_ATLAS(double mass);     //
    double ex_ggF_phi_tautau_ATLAS(double mass);   //
    double ex_bbF_phi_tautau_ATLAS(double mass);   //
    double ex_ggF_A_hZ_tautauZ_ATLAS(double mass); //
    double ex_ggF_A_hZ_bbZ_ATLAS(double mass);     //
    double ex_bbF_phi_bb_CMS(double mass);         //CMS-HIG-14-017, Figure 6
    double ex_pp_phi_tt_ATLAS(double mass);        //arXiv: 1505.07018, Figure 11d
    double ex_ggF_phi_tautau_CMS(double mass);     //CMS-PAS-HIG-14-029, Figure 10-b
    double ex_bbF_phi_tautau_CMS(double mass);     //CMS-PAS-HIG-14-029, Figure 10-b
    double ex_ggF_phi_gaga_CMS(double mass);       //CMS-HIG-14-006, Figure 7a
    double ex_ggF_H_WW_ATLAS(double mass);         //arXiv:1509.00389, Figure 13, left
    double ex_VBF_H_WW_ATLAS(double mass);         //arXiv:1509.00389, Figure 13, right
    double ex_ggF_H_hh_ATLAS(double mass);         //arXiv:1509.04670, Figure 6
    double ex_ggF_H_hh_bbtautau_CMS(double mass);  //arXiv:1510.01181, Figure 8, bottom right
    double ex_ggF_A_hZ_tautaull_CMS(double mass);  //arXiv:1510.01181, Figure 10 left

    double ex_bsgamma(double logtb, double logmHp);

    /*One-loop functions*/

    double A0_MZ2_mHp2(const double MZ2, const double mHp2) const;
    double A0_MZ2_MW2(const double MZ2, const double MW2) const;

    gslpp::complex B0_MZ2_0_MW2_mHh2(const double MZ2, const double MW2, const double mHh2) const;
    gslpp::complex B0_MZ2_0_MW2_mHl2(const double MZ2, const double MW2, const double mHl2) const;
    gslpp::complex B0_MZ2_0_MZ2_mHh2(const double MZ2, const double mHh2) const;
    gslpp::complex B0_MZ2_0_MZ2_mHl2(const double MZ2, const double mHl2) const;
    gslpp::complex B0_MZ2_MW2_MW2_mHh2(const double MZ2, const double MW2, const double mHh2) const;
    gslpp::complex B0_MZ2_MW2_MW2_mHl2(const double MZ2, const double MW2, const double mHl2) const;
    gslpp::complex B0_MZ2_MZ2_MZ2_mHh2(const double MZ2, const double mHh2) const;
    gslpp::complex B0_MZ2_MZ2_MZ2_mHl2(const double MZ2, const double mHl2) const;

    gslpp::complex B00_MZ2_0_mA2_mHp2(const double MZ2, const double mA2, const double mHp2) const;
    gslpp::complex B00_MZ2_0_mHh2_mA2(const double MZ2, const double mHh2, const double mA2) const;
    gslpp::complex B00_MZ2_0_mHh2_mHp2(const double MZ2, const double mHh2, const double mHp2) const;
    gslpp::complex B00_MZ2_0_mHh2_MW2(const double MZ2, const double mHh2, const double MW2) const; //new
    gslpp::complex B00_MZ2_0_mHl2_mA2(const double MZ2, const double mHl2, const double mA2) const;
    gslpp::complex B00_MZ2_0_mHl2_mHp2(const double MZ2, const double mHl2, const double mHp2) const;
    gslpp::complex B00_MZ2_0_mHl2_MW2(const double MZ2, const double mHl2, const double MW2) const; //new
    gslpp::complex B00_MZ2_0_mHp2_mHp2(const double MZ2, const double mHp2) const;
    gslpp::complex B00_MZ2_0_MW2_mHh2(const double MZ2, const double MW2, const double mHh2) const;
    gslpp::complex B00_MZ2_0_MW2_mHl2(const double MZ2, const double MW2, const double mHl2) const;
    gslpp::complex B00_MZ2_0_MW2_MZ2(const double MZ2, const double MW2) const; //new
    gslpp::complex B00_MZ2_0_MZ2_mHh2(const double MZ2, const double mHh2) const;
    gslpp::complex B00_MZ2_0_MZ2_mHl2(const double MZ2, const double mHl2) const;
    gslpp::complex B00_MZ2_0_MZ2_MW2(const double MZ2, const double MW2) const; //new
    gslpp::complex B00_MZ2_MW2_mA2_mHp2(const double MZ2, const double MW2, const double mA2, const double mHp2) const;
    gslpp::complex B00_MZ2_MW2_mHh2_mHp2(const double MZ2, const double MW2, const double mHh2, const double mHp2) const;
    gslpp::complex B00_MZ2_MW2_mHh2_MW2(const double MZ2, const double MW2, const double mHh2) const;   //new
    gslpp::complex B00_MZ2_MW2_mHl2_mHp2(const double MZ2, const double MW2, const double mHl2, const double mHp2) const;
    gslpp::complex B00_MZ2_MW2_mHl2_MW2(const double MZ2, const double MW2, const double mHl2) const;   //new
    gslpp::complex B00_MZ2_MW2_mHp2_mHp2(const double MZ2, const double MW2, const double mHp2) const;
    gslpp::complex B00_MZ2_MW2_MW2_mHh2(const double MZ2, const double MW2, const double mHh2) const;
    gslpp::complex B00_MZ2_MW2_MW2_mHl2(const double MZ2, const double MW2, const double mHl2) const;
    gslpp::complex B00_MZ2_MW2_MW2_MZ2(const double MZ2, const double MW2) const;   //new
    gslpp::complex B00_MZ2_MW2_MZ2_MW2(const double MZ2, const double MW2) const;   //new
    gslpp::complex B00_MZ2_MZ2_mHh2_mA2(const double MZ2, const double mHh2, const double mA2) const;
    gslpp::complex B00_MZ2_MZ2_mHh2_MZ2(const double MZ2, const double mHh2) const; //new
    gslpp::complex B00_MZ2_MZ2_mHl2_mA2(const double MZ2, const double mHl2, const double mA2) const;
    gslpp::complex B00_MZ2_MZ2_mHl2_MZ2(const double MZ2, const double mHl2) const; //new
    gslpp::complex B00_MZ2_MZ2_mHp2_mHp2(const double MZ2, const double mHp2) const;
    gslpp::complex B00_MZ2_MZ2_MZ2_mHh2(const double MZ2, const double mHh2) const;
    gslpp::complex B00_MZ2_MZ2_MZ2_mHl2(const double MZ2, const double mHl2) const;

    gslpp::complex B00p_0_mHp2_mHp2(const double MZ2, const double mHp2) const; /*check whether this is the correct PV function!*/

 ///////////////////////////////////////////////////////////////////////////////

private:

    const PVfunctions PV;
    //const THDM myTHDM;

    ////////////////////////////////////////////////////////////////////////////
    //Caches
    
    /*One-loop functions*/

    mutable double A0_MZ2_MHp2_cache[3][CacheSize];
    mutable double A0_MZ2_MW2_cache[3][CacheSize];

    mutable gslpp::complex B0_MZ2_0_MW2_mHh2_cache[4][CacheSize];
    mutable gslpp::complex B0_MZ2_0_MW2_mHl2_cache[4][CacheSize];
    mutable gslpp::complex B0_MZ2_0_MZ2_mHh2_cache[3][CacheSize];
    mutable gslpp::complex B0_MZ2_0_MZ2_mHl2_cache[3][CacheSize];
    mutable gslpp::complex B0_MZ2_MW2_MW2_mHh2_cache[4][CacheSize];
    mutable gslpp::complex B0_MZ2_MW2_MW2_mHl2_cache[4][CacheSize];
    mutable gslpp::complex B0_MZ2_MZ2_MZ2_mHh2_cache[3][CacheSize];
    mutable gslpp::complex B0_MZ2_MZ2_MZ2_mHl2_cache[3][CacheSize];

    mutable gslpp::complex B00_MZ2_0_mA2_mHp2_cache[4][CacheSize];
    mutable gslpp::complex B00_MZ2_0_mHh2_mA2_cache[4][CacheSize];
    mutable gslpp::complex B00_MZ2_0_mHh2_mHp2_cache[4][CacheSize];
    mutable gslpp::complex B00_MZ2_0_mHh2_MW2_cache[4][CacheSize]; //new
    mutable gslpp::complex B00_MZ2_0_mHl2_mA2_cache[4][CacheSize];
    mutable gslpp::complex B00_MZ2_0_mHl2_mHp2_cache[4][CacheSize];
    mutable gslpp::complex B00_MZ2_0_mHl2_MW2_cache[4][CacheSize]; //new
    mutable gslpp::complex B00_MZ2_0_mHp2_mHp2_cache[3][CacheSize];
    mutable gslpp::complex B00_MZ2_0_MW2_mHh2_cache[4][CacheSize];
    mutable gslpp::complex B00_MZ2_0_MW2_mHl2_cache[4][CacheSize];
    mutable gslpp::complex B00_MZ2_0_MW2_MZ2_cache[3][CacheSize]; //new
    mutable gslpp::complex B00_MZ2_0_MZ2_mHh2_cache[3][CacheSize];
    mutable gslpp::complex B00_MZ2_0_MZ2_mHl2_cache[3][CacheSize];
    mutable gslpp::complex B00_MZ2_0_MZ2_MW2_cache[3][CacheSize]; //new
    mutable gslpp::complex B00_MZ2_MW2_mA2_mHp2_cache[5][CacheSize];
    mutable gslpp::complex B00_MZ2_MW2_mHh2_mHp2_cache[5][CacheSize];
    mutable gslpp::complex B00_MZ2_MW2_mHh2_MW2_cache[4][CacheSize]; //new
    mutable gslpp::complex B00_MZ2_MW2_mHl2_mHp2_cache[5][CacheSize];
    mutable gslpp::complex B00_MZ2_MW2_mHl2_MW2_cache[4][CacheSize]; //new
    mutable gslpp::complex B00_MZ2_MW2_mHp2_mHp2_cache[4][CacheSize];
    mutable gslpp::complex B00_MZ2_MW2_MW2_mHh2_cache[4][CacheSize];
    mutable gslpp::complex B00_MZ2_MW2_MW2_mHl2_cache[4][CacheSize];
    mutable gslpp::complex B00_MZ2_MW2_MW2_MZ2_cache[3][CacheSize]; //new
    mutable gslpp::complex B00_MZ2_MW2_MZ2_MW2_cache[3][CacheSize]; //new
    mutable gslpp::complex B00_MZ2_MZ2_mHh2_mA2_cache[4][CacheSize];
    mutable gslpp::complex B00_MZ2_MZ2_mHh2_MZ2_cache[3][CacheSize]; //new
    mutable gslpp::complex B00_MZ2_MZ2_mHl2_mA2_cache[4][CacheSize];
    mutable gslpp::complex B00_MZ2_MZ2_mHl2_MZ2_cache[3][CacheSize]; //new
    mutable gslpp::complex B00_MZ2_MZ2_mHp2_mHp2_cache[3][CacheSize];
    mutable gslpp::complex B00_MZ2_MZ2_MZ2_mHh2_cache[3][CacheSize];
    mutable gslpp::complex B00_MZ2_MZ2_MZ2_mHl2_cache[3][CacheSize];

    mutable gslpp::complex B00p_0_mHp2_mHp2_cache[3][CacheSize]; /*check whether this is the correct PV function!*/
};

#endif	/* THDMCACHE_H */
