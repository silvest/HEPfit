/* 
 * Copyright (C) 2016 HEPfit Collaboration
 * All rights reserved.
 *
 * For the licensing terms see doc/COPYING.
 */

#ifndef GENERALTHDMSTABILITY_H
#define GENERALTHDMSTABILITY_H

#include "ThObservable.h"
#include <gslpp.h>

class GeneralTHDM;
class GeneralTHDMcache;

/**
 * This computes the bounded from below conditions according to [2203.11462]
 */
class stability_GTHDM {
public:
    
    /**
     * @brief Constructor.
     */
    stability_GTHDM(const StandardModel& SM_i);
    
    /**
     * @brief Destructor.
     */
    virtual ~stability_GTHDM();
    
    /**
     * @brief Getter function for the stability conditions.
     */
    gslpp::vector<double> getStability();
    
    
    /**
    * @brief Computes the eigenvalues of the @f$\Lambda^{\mu}_{\nu}$ matrix of [hep-ph/0609018]
    */
    bool CalcStabeigen(gslpp::matrix<gslpp::complex>& Stabeigvec_i, gslpp::vector<double>& Stabeigval_i);
    
    /**
     * @brief Getter function for the vacuum stability condition.
     */
    double getVacuumStability();
    
    
private:
    const GeneralTHDM& myGTHDM;
    
    gslpp::matrix<gslpp::complex> Lambmat, Lambeigvec;
    gslpp::vector<double> vecMinus1, vecStability, Lambeigval;

};










/**
 * 
 */
class stability1_GTHDM: public ThObservable {
public:

    /**
     * @brief stability1_GTHDM constructor.
     */
    stability1_GTHDM(const StandardModel& SM_i);

    /**
     * @return
     */
    double computeThValue();

private:
    stability_GTHDM mystability_GTHDM;
};


/**
 * 
 */
class stability2_GTHDM: public ThObservable {
public:

    /**
     * @brief stability2_GTHDM constructor.
     */
    stability2_GTHDM(const StandardModel& SM_i);

    /**
     * @return
     */
    double computeThValue();

private:
    stability_GTHDM mystability_GTHDM;
};


/**
 * 
 */
class stability3_GTHDM: public ThObservable {
public:

    /**
     * @brief stability3_GTHDM constructor.
     */
    stability3_GTHDM(const StandardModel& SM_i);

    /**
     * @return
     */
    double computeThValue();

private:
    stability_GTHDM mystability_GTHDM;
};


/**
 * 
 */
class stability4_GTHDM: public ThObservable {
public:

    /**
     * @brief stability4_GTHDM constructor.
     */
    stability4_GTHDM(const StandardModel& SM_i);

    /**
     * @return
     */
    double computeThValue();

private:
    stability_GTHDM mystability_GTHDM;
};



/**
 * 
 */
class stability5_GTHDM: public ThObservable {
public:

    /**
     * @brief stability5_GTHDM constructor.
     */
    stability5_GTHDM(const StandardModel& SM_i);

    /**
     * @return
     */
    double computeThValue();

private:
    stability_GTHDM mystability_GTHDM;
};



/**
 * 
 */
class stability6_GTHDM: public ThObservable {
public:

    /**
     * @brief stability6_GTHDM constructor.
     */
    stability6_GTHDM(const StandardModel& SM_i);

    /**
     * @return
     */
    double computeThValue();

private:
    stability_GTHDM mystability_GTHDM;
};



/**
 * 
 */
class stability7_GTHDM: public ThObservable {
public:

    /**
     * @brief stability7_GTHDM constructor.
     */
    stability7_GTHDM(const StandardModel& SM_i);

    /**
     * @return
     */
    double computeThValue();

private:
    stability_GTHDM mystability_GTHDM;
};



/**
 * 
 */
class vacuumstability_GTHDM: public ThObservable {
public:

    /**
     * @brief vacuumstability_GTHDM constructor.
     */
    vacuumstability_GTHDM(const StandardModel& SM_i);

    /**
     * @return
     */
    double computeThValue();

private:
    stability_GTHDM mystability_GTHDM;
};





#endif /* GENERALTHDMSTABILITY_H */

