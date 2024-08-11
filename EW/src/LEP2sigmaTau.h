/* 
 * Copyright (C) 2012 HEPfit Collaboration
 *
 *
 * For the licensing terms see doc/COPYING.
 */

#ifndef LEP2SIGMATAU_H
#define	LEP2SIGMATAU_H

#include "LEP2ThObservable.h"

/**
 * @class LEP2sigmaTau
 * @ingroup EW
 * @brief A class for the cross section of @f$e^+e^-\to \tau^+\tau^-@f$ above the @f$Z@f$ pole.
 * @author HEPfit Collaboration
 * @copyright GNU General Public License
 * @details 
 */
class LEP2sigmaTau : public LEP2ThObservable {
public:

    /**
     * @brief LEP2sigmaTau constructor
     * @param[in] SM_i a reference to an object of type StandardModel
     * @param[in] sqrt_s_i the CM energy of the e^+ e^- pair
     * @param[in] bSigmaForAFB_i true for the denominator of A_FB
     */
    LEP2sigmaTau(const StandardModel& SM_i, const double sqrt_s_i) 
    : LEP2ThObservable(SM_i, sqrt_s_i) 
    {

    }

    /**
     * @return the cross section for e^+ e^- -> tau^+ tau^- at sqrt_s in pb
     */
    double computeThValue();

private:

};


/**
 * @class LEP2dsigmadcosTau
 * @ingroup EW
 * @brief A class for the cross section of @f$e^+e^-\to \tau^+\tau^-@f$ above the @f$Z@f$ pole.
 * @author HEPfit Collaboration
 * @copyright GNU General Public License
 * @details 
 */
class LEP2dsigmadcosTau : public LEP2ThDiffObservable {
public:

    /**
     * @brief LEP2dsigmadcosTau constructor
     * @param[in] SM_i a reference to an object of type StandardModel
     * @param[in] sqrt_s_i the CM energy of the e^+ e^- pair
     * @param[in] cos_i the polar angle of the final state particle wrt e^-
     */
    LEP2dsigmadcosTau(const StandardModel& SM_i, const double sqrt_s_i, const double cos_i) 
    : LEP2ThDiffObservable(SM_i, sqrt_s_i, cos_i) 
    {

    }

    /**
     * @return the cross section for e^+ e^- -> tau^+ tau^- at sqrt_s in pb
     */
    double computeThValue();

private:

};

#endif	/* LEP2SIGMATAU_H */

