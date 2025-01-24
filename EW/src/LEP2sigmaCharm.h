/* 
 * Copyright (C) 2012 HEPfit Collaboration
 *
 *
 * For the licensing terms see doc/COPYING.
 */

#ifndef LEP2SIGMACHARM_H
#define	LEP2SIGMACHARM_H

#include "LEP2ThObservable.h"

/**
 * @class LEP2sigmaCharm
 * @ingroup EW
 * @brief A class for the cross section of @f$e^+e^-\to c\bar{c}@f$ above the @f$Z@f$ pole.
 * @author HEPfit Collaboration
 * @copyright GNU General Public License
 * @details 
 */
class LEP2sigmaCharm : public LEP2ThObservable {
public:

    /**
     * @brief LEP2sigmaCharm constructor
     * @param[in] SM_i a reference to an object of type StandardModel
     * @param[in] sqrt_s_i the CM energy of the e^+ e^- pair
     * @param[in] bSigmaForAFB_i true for the denominator of A_FB
     */
    LEP2sigmaCharm(const StandardModel& SM_i, const double sqrt_s_i) 
    : LEP2ThObservable(SM_i, sqrt_s_i) 
    {

    }

    /**
     * @return the cross section for e^+ e^- -> c cbar at sqrt_s in pb
     */
    double computeThValue();

private:
    
};


/**
 * @class eeffsigmaCharm
 * @ingroup EW
 * @brief A class for the cross section of @f$e^+e^-\to c\bar{c}@f$ off the @f$Z@f$ pole.
 * @author HEPfit Collaboration
 * @copyright GNU General Public License
 * @details 
 */
class eeffsigmaCharm : public eeffThObservable {
public:

    /**
     * @brief eeffsigmaCharm constructor
     * @param[in] SM_i a reference to an object of type StandardModel
     * @param[in] pol_e_i, pol_p_i the electron and positron polarizations
     * @param[in] sqrt_s_i the CM energy of the e^+ e^- pair
     */
    eeffsigmaCharm(const StandardModel& SM_i, const double pol_e_i, const double pol_p_i, const double sqrt_s_i) 
    : eeffThObservable(SM_i, pol_e_i, pol_p_i, sqrt_s_i) 
    {

    }

    /**
     * @return the cross section for e^+ e^- -> c cbar at sqrt_s in pb
     */
    double computeThValue();

private:
    
};


/**
 * @class eeffRcharm
 * @ingroup EW
 * @brief A class for the ratio of cross sections of @f$R_c=\sigma(e^+e^-\to c\bar{c})/\sigma(e^+e^-\to hadrons)@f$ off the @f$Z@f$ pole.
 * @author HEPfit Collaboration
 * @copyright GNU General Public License
 * @details 
 */
class eeffRcharm : public eeffThObservable {
public:

    /**
     * @brief eeffRcharm constructor
     * @param[in] SM_i a reference to an object of type StandardModel
     * @param[in] pol_e_i, pol_p_i the electron and positron polarizations
     * @param[in] sqrt_s_i the CM energy of the e^+ e^- pair
     */
    eeffRcharm(const StandardModel& SM_i, const double pol_e_i, const double pol_p_i, const double sqrt_s_i) 
    : eeffThObservable(SM_i, pol_e_i, pol_p_i, sqrt_s_i)
    {
    }

    /**
     * @return the ratio @f$R_c=\sigma(e^+e^-\to c\bar{c})/\sigma(e^+e^-\to hadrons)@f$ at sqrt_s
     */
    double computeThValue();

private:
    
};

#endif	/* LEP2SIGMACHARM_H */

