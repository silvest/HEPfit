/* 
 * Copyright (C) 2016 HEPfit Collaboration
 * All rights reserved.
 *
 * For the licensing terms see doc/COPYING.
 */

#ifndef GENERALTHDMSTABILITY_H
#define GENERALTHDMSTABILITY_H

#include "ThObservable.h"
#include "GeneralTHDM.h"
#include "GeneralTHDMcache.h"
#include <gslpp.h>

/**
 * 
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
     * @brief Destructor.
     */
    gslpp::vector<double> getStability();
    
private:
    const GeneralTHDM& myGTHDM;
    gslpp::vector<double> vecMinus1, vecStability;
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

#endif /* GENERALTHDMSTABILITY_H */

