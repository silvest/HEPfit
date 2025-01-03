/* 
 * Copyright (C) 2012 HEPfit Collaboration
 *
 *
 * For the licensing terms see doc/COPYING.
 */

#ifndef LEP2SIGMABOTTOM_H
#define	LEP2SIGMABOTTOM_H

#include "LEP2ThObservable.h"

/**
 * @class LEP2sigmaBottom
 * @ingroup EW
 * @brief A class for the cross section of @f$e^+e^-\to b\bar{b}@f$ above the @f$Z@f$ pole.
 * @author HEPfit Collaboration
 * @copyright GNU General Public License
 * @details 
 */
class LEP2sigmaBottom : public LEP2ThObservable {
public:

    /**
     * @brief LEP2sigmaBottom constructor
     * @param[in] SM_i a reference to an object of type StandardModel
     * @param[in] sqrt_s_i the CM energy of the e^+ e^- pair
     * @param[in] bSigmaForAFB_i true for the denominator of A_FB
     */
    LEP2sigmaBottom(const StandardModel& SM_i, const double sqrt_s_i) 
    : LEP2ThObservable(SM_i, sqrt_s_i)
    {
    }

    /**
     * @return the cross section for e^+ e^- -> b bbar at sqrt_s in pb
     */
    double computeThValue();

private:
    
};


/**
 * @class eeffsigmaBottom
 * @ingroup EW
 * @brief A class for the cross section of @f$e^+e^-\to b\bar{b}@f$ above the @f$Z@f$ pole.
 * @author HEPfit Collaboration
 * @copyright GNU General Public License
 * @details 
 */
class eeffsigmaBottom : public eeffThObservable {
public:

    /**
     * @brief eeffsigmaBottom constructor
     * @param[in] SM_i a reference to an object of type StandardModel
     * @param[in] pol_e_i, pol_p_i the electron and positron polarizations
     * @param[in] sqrt_s_i the CM energy of the e^+ e^- pair
     */
    eeffsigmaBottom(const StandardModel& SM_i, const double pol_e_i, const double pol_p_i, const double sqrt_s_i) 
    : eeffThObservable(SM_i, pol_e_i, pol_p_i, sqrt_s_i)
    {
    }

    /**
     * @return the cross section for e^+ e^- -> b bbar at sqrt_s in pb
     */
    double computeThValue();

private:
    
};

#endif	/* LEP2SIGMABOTTOM_H */

