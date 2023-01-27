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
    * @brief Computes the eigenvalues of the @f$\Lambda^{\mu}_{\nu}$ matrix of [hep-ph/0609018]
    */
    bool CalcStabeigen(gslpp::matrix<gslpp::complex>& Stabeigvec_i, gslpp::vector<gslpp::complex>& Stabeigval_i);
    
    
    
    
    
    
    /**
     * @brief Getter function for the stability conditions.
     */
    gslpp::vector<double> getStability();
    
    
    
    
    
private:
    const GeneralTHDM& myGTHDM;
    
    gslpp::matrix<double> LambmatE;
    gslpp::matrix<gslpp::complex> Lambeigvec;
    gslpp::vector<double> vecMinus1, vecStability;
    gslpp::vector<gslpp::complex> Lambeigval;

};










/**
 * 
 */
class bounded_from_below_GTHDM: public ThObservable {
public:

    /**
     * @brief bounded_from_below_GTHDM constructor.
     */
    bounded_from_below_GTHDM(const StandardModel& SM_i);

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
class vacuum_stability_GTHDM: public ThObservable {
public:

    /**
     * @brief vacuum_stability_GTHDM constructor.
     */
    vacuum_stability_GTHDM(const StandardModel& SM_i);

    /**
     * @return
     */
    double computeThValue();

private:
    stability_GTHDM mystability_GTHDM;
};





#endif /* GENERALTHDMSTABILITY_H */

