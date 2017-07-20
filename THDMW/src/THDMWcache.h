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
    
    /**
     * @brief Cache size.
     * @details Determines the size of the cache. If it is set to 5, the cache will remember the last five function calls and store their results.
     */
//    static const int CacheSize = 5;

    /**
     * @brief Check whether for the latest set of parameters a value is in the cache.
     * @details Takes a complex value.
     */
//    int CacheCheck(const gslpp::complex cache[][CacheSize],
//                   const int NumPar, const double params[]) const;

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
//    void CacheShift(gslpp::complex cache[][CacheSize], const int NumPar,
//                    const double params[], const gslpp::complex newResult) const; 

    /**
     * @brief Adds a new result and its parameters into the cache.
     * @details The new values are added on top. The oldest set on the stack is deleted. Takes a real value.
     */
//    void CacheShiftReal(double cache[][CacheSize], const int NumPar,
//                    const double params[], const double newResult) const; 

    void updateCache();
    void setOtherParameters();

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

    double RpepsTHDMW;
    gslpp::vector<gslpp::complex> unitarityeigenvalues;
    gslpp::vector<gslpp::complex> NLOunitarityeigenvalues;

private:

    const THDMW * myTHDMW;
    RunnerTHDMW * myRunnerTHDMW;

    void runTHDMWparameters();

    void computeUnitarity();
    gslpp::vector<gslpp::complex> betaeigenvalues;

    double Q_THDMW;
    double MZ;
    double lambda1;
    double lambda2;
    double lambda3;
    double lambda4;
    double mu1;
    double mu3;
    double mu4;
    double nu1;
    double omega1;
    double kappa1;
    double nu2;
    double omega2;
    double kappa2;
    double nu4;
    double omega4;
    double mHpsq;
    double mAsq;
    double mhsq;
    double mHsq;
    double m12sq;
};

#endif	/* THDMWCACHE_H */
