/* 
 * Copyright (C) 2016 HEPfit Collaboration
 * All rights reserved.
 *
 * For the licensing terms see doc/COPYING.
 */

#ifndef GENERALTHDMEFFECTIVEPOT_H
#define GENERALTHDMEFFECTIVEPOT_H

#include "ThObservable.h"
#include "gslpp.h"

class GeneralTHDM;
class GeneralTHDMcache;

/**
 * 
 */
class EffectivePot_GTHDM {
public:

    /**
     * @brief EffectivePot_GTHDM constructor.
     */
    EffectivePot_GTHDM(const StandardModel& SM_i);

    double potentialfunction(const double *Svec);
    
    gslpp::matrix<double> Mneutral_2(const double S1, const double S2, const double S3);
    
    const double* potentialminimizer(double S1_start, double S2_start, double S3_start);

private:
    
//    const GeneralTHDM& myGTHDM;
    gslpp::matrix<double> mat_neutral;
};

/**
 * 
 */
class EffectivePotMin1_GTHDM: public ThObservable {
public:

    /**
     * @brief EffectivePotMin1_GTHDM constructor.
     */
    EffectivePotMin1_GTHDM(const StandardModel& SM_i);

    /**
     * @return
     */
    double computeThValue();

private:
    EffectivePot_GTHDM myEffectivePot_GTHDM;
};

/**
 * 
 */
class EffectivePotMin2_GTHDM: public ThObservable {
public:

    /**
     * @brief EffectivePotMin2_GTHDM constructor.
     */
    EffectivePotMin2_GTHDM(const StandardModel& SM_i);

    /**
     * @return
     */
    double computeThValue();

private:
    EffectivePot_GTHDM myEffectivePot_GTHDM;
};

#endif /* GENERALTHDMEFFECTIVEPOT_H */

