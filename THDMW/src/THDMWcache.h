/* 
 * Copyright (C) 2017 HEPfit Collaboration
 * All rights reserved.
 *
 * For the licensing terms see doc/COPYING.
 */

#ifndef THDMWCACHE_H
#define	THDMWCACHE_H

#include "THDMW.h"
#include "RunnerTHDMW.h"

/**
 * @class THDMWcache
 * @ingroup THDMW
 * @brief A class for the caching of some %THDMW objects.
 * @author HEPfit Collaboration
 * @copyright GNU General Public License
 */
class THDMWcache {
    
public:

    /**
     * @brief THDMWcache constructor.
     * @details Reads all the tables values and stores them in the memory.
     */
    THDMWcache(const StandardModel& SM_i);

    /**
     * @brief THDMWcache destructor.
     */
    ~THDMWcache();

    void updateCache();
    double setOtherParameters();

    double Q_cutoff;
//    double g1_at_Q;
//    double g2_at_Q;
//    double g3_at_Q;
//    double Ytop_at_Q;
//    double Ybottom1_at_Q;
//    double Ybottom2_at_Q;
//    double Ytau1_at_Q;
//    double Ytau2_at_Q;
    double lambda1_at_Q;
    double lambda2_at_Q;
    double lambda3_at_Q;
    double lambda4_at_Q;
    double mu1_at_Q;
    double mu3_at_Q;
    double mu4_at_Q;
    double nu1_at_Q;
    double omega1_at_Q;
    double kappa1_at_Q;
    double nu2_at_Q;
    double omega2_at_Q;
    double kappa2_at_Q;
    double nu4_at_Q;
    double omega4_at_Q;
    double m12sq;
    double m11sq;
    double m22sq;
    double mhsq;
    double mHsq;
    double mAsq;
    double mSRsq;
    double mSIsq;
    double mHpsq;
    double mSpsq;

    double RpepsTHDMW;
    gslpp::vector<gslpp::complex> unitarityeigenvalues;
    gslpp::vector<gslpp::complex> NLOunitarityeigenvalues;

    double rh_QuQu, rh_VV, rh_gg, rh_QdQd, rh_ll, rh_gaga, rh_Zga;
    double sumModBRs, Gamma_h, THDM_BR_h_bb, THDM_BR_h_gaga, THDM_BR_h_tautau, THDM_BR_h_WW, THDM_BR_h_ZZ;
    
protected:

private:

    const THDMW * myTHDMW;
    RunnerTHDMW * myRunnerTHDMW;

    /**
     * @brief Cache size.
     * @details Determines the size of the cache. If it is set to 5, the cache will remember the last five function calls and store their results.
     */
    static const int CacheSize = 5;

    /**
     * @brief Check whether for the latest set of parameters a value is in the cache.
     * @details Takes a complex value.
     */
    int CacheCheck(const gslpp::complex cache[][CacheSize],
                   const int NumPar, const double params[]) const;

    /**
     * @brief Check whether for the latest set of parameters a value is in the cache.
     * @details Takes a real value.
     */
//    int CacheCheckReal(const double cache[][CacheSize],
//                   const int NumPar, const double params[]) const;

    /**
     * @brief Adds a new result and its parameters into the cache.
     * @details The new values are added on top. The oldest set on the stack is deleted. Takes a complex value.
     */
    void CacheShift(gslpp::complex cache[][CacheSize], const int NumPar,
                    const double params[], const gslpp::complex newResult) const; 

    /**
     * @brief Adds a new result and its parameters into the cache.
     * @details The new values are added on top. The oldest set on the stack is deleted. Takes a real value.
     */
//    void CacheShiftReal(double cache[][CacheSize], const int NumPar,
//                    const double params[], const double newResult) const; 

    gslpp::complex I_h_U(const double mHl2, const double Mu, const double Mc, const double Mt) const;
    gslpp::complex I_HH_U(const double mHh2, const double Mc, const double Mt) const;
    gslpp::complex I_A_U(const double mA2, const double Mc, const double Mt) const;
    gslpp::complex I_h_D(const double mHl2, const double Md, const double Ms, const double Mb) const;
    gslpp::complex I_HH_D(const double mHh2, const double Ms, const double Mb) const;
    gslpp::complex I_A_D(const double mA2, const double Ms, const double Mb) const;
    gslpp::complex I_h_L(const double mHl2, const double Me, const double Mmu, const double Mtau) const;
    gslpp::complex I_HH_L(const double mHh2, const double Mmu, const double Mtau) const;
    gslpp::complex I_A_L(const double mA2, const double Mmu, const double Mtau) const;
    gslpp::complex I_H_W(const double mH, const double MW) const;
    gslpp::complex I_H_Hp(const double mHp2, const double mH) const;

    gslpp::complex A_h_U(const double mHl2, const double cW2, const double Mu, const double Mc, const double Mt, const double MZ) const;
    gslpp::complex A_HH_U(const double mHh2, const double cW2, const double Mc, const double Mt, const double MZ) const;
    gslpp::complex A_A_U(const double mA2, const double cW2, const double Mc, const double Mt, const double MZ) const;
    gslpp::complex A_h_D(const double mHl2, const double cW2, const double Md, const double Ms, const double Mb, const double MZ) const;
    gslpp::complex A_HH_D(const double mHh2, const double cW2, const double Ms, const double Mb, const double MZ) const;
    gslpp::complex A_A_D(const double mA2, const double cW2, const double Ms, const double Mb, const double MZ) const;
    gslpp::complex A_h_L(const double mHl2, const double cW2, const double Me, const double Mmu, const double Mtau, const double MZ) const;
    gslpp::complex A_HH_L(const double mHh2, const double cW2, const double Mmu, const double Mtau, const double MZ) const;
    gslpp::complex A_A_L(const double mA2, const double cW2, const double Mmu, const double Mtau, const double MZ) const;
    gslpp::complex A_H_W(const double mH, const double cW2, const double MW, const double MZ) const;
    gslpp::complex A_H_Hp(const double mHp2, const double mH, const double cW2, const double MZ) const;

    mutable gslpp::complex I_h_U_cache[5][CacheSize];
    mutable gslpp::complex I_HH_U_cache[4][CacheSize];
    mutable gslpp::complex I_A_U_cache[4][CacheSize];
    mutable gslpp::complex I_h_D_cache[5][CacheSize];
    mutable gslpp::complex I_HH_D_cache[4][CacheSize];
    mutable gslpp::complex I_A_D_cache[4][CacheSize];
    mutable gslpp::complex I_h_L_cache[5][CacheSize];
    mutable gslpp::complex I_HH_L_cache[4][CacheSize];
    mutable gslpp::complex I_A_L_cache[4][CacheSize];
    mutable gslpp::complex I_H_W_cache[3][CacheSize];
    mutable gslpp::complex I_H_Hp_cache[3][CacheSize];

    mutable gslpp::complex A_h_U_cache[7][CacheSize];
    mutable gslpp::complex A_HH_U_cache[6][CacheSize];
    mutable gslpp::complex A_A_U_cache[6][CacheSize];
    mutable gslpp::complex A_h_D_cache[7][CacheSize];
    mutable gslpp::complex A_HH_D_cache[6][CacheSize];
    mutable gslpp::complex A_A_D_cache[6][CacheSize];
    mutable gslpp::complex A_h_L_cache[7][CacheSize];
    mutable gslpp::complex A_HH_L_cache[6][CacheSize];
    mutable gslpp::complex A_A_L_cache[6][CacheSize];
    mutable gslpp::complex A_H_W_cache[5][CacheSize];
    mutable gslpp::complex A_H_Hp_cache[5][CacheSize];

    gslpp::complex f_func(const double x) const;
    gslpp::complex g_func(const double x) const;

    gslpp::complex Int1(const double tau, const double lambda) const;
    gslpp::complex Int2(const double tau, const double lambda) const;

    void runTHDMWparameters();

    void computeUnitarity();
    gslpp::vector<gslpp::complex> betaeigenvalues;

    void computeSignalStrengthQuantities();

    std::string THDMWmodel;
    double Q_THDMW;
    double MZ;
    double vev;
    double tanb;
    double sinb;
    double cosb;
    double bma;
    double sina;
    double cosa;
    double lambda1;
    double lambda2;
    double lambda3;
    double lambda4;
    double lambda5;
    double mSsq;
    double mu1;
    double mu2;
    double mu3;
    double mu4;
    double mu5;
    double mu6;
    double nu1;
    double omega1;
    double kappa1;
    double nu2;
    double omega2;
    double kappa2;
    double nu3;
    double omega3;
    double kappa3;
    double nu4;
    double omega4;
};

#endif	/* THDMWCACHE_H */
