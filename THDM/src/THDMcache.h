/* 
 * Copyright (C) 2012 SUSYfit Collaboration
 * All rights reserved.
 *
 * For the licensing terms see doc/COPYING.
 */

#ifndef THDMCACHE_H
#define	THDMCACHE_H

#include <cmath>
#include <PVfunctions.h>
//#include "THDM.h"

using namespace gslpp;


class THDMcache {
    
public:
    
    /**
     * @brief THDMcache constructor
     * @param[in] SM_i reference to a StandardModel object
     * @param[in] THDM_i reference to a THDM object
     */
    THDMcache();
    
    static const int CacheSize = 5;

    
    int CacheCheck(const complex cache[][CacheSize], 
                   const int NumPar, const double params[]) const;
    
    void CacheShift(complex cache[][CacheSize], const int NumPar, 
                    const double params[], const complex newResult) const; 

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
    
    
    
    /*One-loop functions*/
    
    complex B0_Mz_Mz2_Mz_mH(const double Mz, const double mH) const;
    complex B0_Mz_0_Mz_mH(const double Mz, const double mH) const;
    complex B0_Mz_Mz2_Mz_mh(const double Mz, const double mh) const;
    complex B0_Mz_0_Mz_mh(const double Mz, const double mh) const;
    complex B0_Mz_Mw2_Mw_mH(const double Mz, const double Mw, const double mH) const;
    complex B0_Mz_0_Mw_mH(const double Mz, const double Mw, const double mH) const;
    complex B0_Mz_Mw2_Mw_mh(const double Mz, const double Mw, const double mh) const;
    complex B0_Mz_0_Mw_mh(const double Mz, const double Mw, const double mh) const;
    
    complex B22_Mz_Mz2_mH_mA(const double Mz, const double mH, const double mA) const;
    complex B22_Mz_0_mH_mA(const double Mz, const double mH, const double mA) const;
    complex B22_Mz_Mz2_mHp_mHp(const double Mz, const double mHp) const;
    complex B22_Mz_0_mHp_mHp(const double Mz, const double mHp) const;
    complex B22_Mz_Mz2_mh_mA(const double Mz, const double mh, const double mA) const;
    complex B22_Mz_0_mh_mA(const double Mz, const double mh, const double mA) const;
    complex B22_Mz_Mz2_Mz_mH(const double Mz, const double mH) const;
    complex B22_Mz_0_Mz_mH(const double Mz, const double mH) const;
    complex B22_Mz_Mz2_Mz_mh(const double Mz, const double mh) const;
    complex B22_Mz_0_Mz_mh(const double Mz, const double mh) const;
    complex B22_Mz_Mw2_mA_mHp(const double Mz, const double Mw, const double mA, const double mHp) const;
    complex B22_Mz_0_mA_mHp(const double Mz, const double mA, const double mHp) const;
    complex B22_Mz_Mw2_mHp_mHp(const double Mz, const double Mw, const double mHp) const;
    complex B22_Mz_Mw2_Mw_mH(const double Mz, const double Mw, const double mH) const;
    complex B22_Mz_0_Mw_mH(const double Mz, const double Mw, const double mH) const;
    complex B22_Mz_Mw2_Mw_mh(const double Mz, const double Mw, const double mh) const;
    complex B22_Mz_0_Mw_mh(const double Mz, const double Mw, const double mh) const;
    complex B22_Mz_Mw2_mH_mHp(const double Mz, const double Mw, const double mH, const double mHp) const; 
    complex B22_Mz_0_mH_mHp(const double Mz, const double mH, const double mHp) const;
    complex B22_Mz_Mw2_mh_mHp(const double Mz, const double Mw, const double mh, const double mHp) const;
    complex B22_Mz_0_mh_mHp(const double Mz, const double mh, const double mHp) const;
    
 ///////////////////////////////////////////////////////////////////////////////   
    
private:

    const PVfunctions PV;
    //const THDM myTHDM;
    
    
    ////////////////////////////////////////////////////////////////////////////
    //Caches
    
    /*One-loop functions*/
    mutable complex B0_Mz_Mz2_Mz_mH_cache[3][CacheSize];
    mutable complex B0_Mz_0_Mz_mH_cache[3][CacheSize];
    mutable complex B0_Mz_Mz2_Mz_mh_cache[3][CacheSize];
    mutable complex B0_Mz_0_Mz_mh_cache[3][CacheSize];
    mutable complex B0_Mz_Mw2_Mw_mH_cache[4][CacheSize];
    mutable complex B0_Mz_0_Mw_mH_cache[4][CacheSize];
    mutable complex B0_Mz_Mw2_Mw_mh_cache[4][CacheSize];
    mutable complex B0_Mz_0_Mw_mh_cache[4][CacheSize];
    
    mutable complex B22_Mz_Mz2_mH_mA_cache[4][CacheSize];
    mutable complex B22_Mz_0_mH_mA_cache[4][CacheSize];
    mutable complex B22_Mz_Mz2_mHp_mHp_cache[3][CacheSize];
    mutable complex B22_Mz_0_mHp_mHp_cache[3][CacheSize];
    mutable complex B22_Mz_Mz2_mh_mA_cache[4][CacheSize];
    mutable complex B22_Mz_0_mh_mA_cache[4][CacheSize];
    mutable complex B22_Mz_Mz2_Mz_mH_cache[3][CacheSize];
    mutable complex B22_Mz_0_Mz_mH_cache[3][CacheSize];
    mutable complex B22_Mz_Mz2_Mz_mh_cache[3][CacheSize];
    mutable complex B22_Mz_0_Mz_mh_cache[3][CacheSize];
    mutable complex B22_Mz_Mw2_mA_mHp_cache[5][CacheSize];
    mutable complex B22_Mz_0_mA_mHp_cache[4][CacheSize];
    mutable complex B22_Mz_Mw2_mHp_mHp_cache[4][CacheSize];
    mutable complex B22_Mz_Mw2_Mw_mH_cache[4][CacheSize];
    mutable complex B22_Mz_0_Mw_mH_cache[4][CacheSize];
    mutable complex B22_Mz_Mw2_Mw_mh_cache[4][CacheSize];
    mutable complex B22_Mz_0_Mw_mh_cache[4][CacheSize];
    mutable complex B22_Mz_Mw2_mH_mHp_cache[5][CacheSize]; 
    mutable complex B22_Mz_0_mH_mHp_cache[4][CacheSize];
    mutable complex B22_Mz_Mw2_mh_mHp_cache[5][CacheSize];
    mutable complex B22_Mz_0_mh_mHp_cache[4][CacheSize];

};

#endif	/* THDMCACHE_H */

