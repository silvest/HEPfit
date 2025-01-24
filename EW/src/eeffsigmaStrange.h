/* 
 * Copyright (C) 2012 HEPfit Collaboration
 *
 *
 * For the licensing terms see doc/COPYING.
 */

#ifndef EEFFSIGMASTRANGE_H
#define	EEFFSIGMASTRANGE_H

#include "LEP2ThObservable.h"


/**
 * @class eeffsigmaStrange
 * @ingroup EW
 * @brief A class for the cross section of @f$e^+e^-\to s\bar{s}@f$ off the @f$Z@f$ pole.
 * @author HEPfit Collaboration
 * @copyright GNU General Public License
 * @details 
 */
class eeffsigmaStrange : public eeffThObservable {
public:

    /**
     * @brief eeffsigmaStrange constructor
     * @param[in] SM_i a reference to an object of type StandardModel
     * @param[in] pol_e_i, pol_p_i the electron and positron polarizations
     * @param[in] sqrt_s_i the CM energy of the e^+ e^- pair
     */
    eeffsigmaStrange(const StandardModel& SM_i, const double pol_e_i, const double pol_p_i, const double sqrt_s_i) 
    : eeffThObservable(SM_i, pol_e_i, pol_p_i, sqrt_s_i)
    {
    }

    /**
     * @return the cross section for e^+ e^- -> s sbar at sqrt_s in pb
     */
    double computeThValue();

private:
    
};


/**
 * @class eeffRstrange
 * @ingroup EW
 * @brief A class for the ratio of cross sections of @f$R_s=\sigma(e^+e^-\to s\bar{s})/\sigma(e^+e^-\to hadrons)@f$ off the @f$Z@f$ pole.
 * @author HEPfit Collaboration
 * @copyright GNU General Public License
 * @details 
 */
class eeffRstrange : public eeffThObservable {
public:

    /**
     * @brief eeffRstrange constructor
     * @param[in] SM_i a reference to an object of type StandardModel
     * @param[in] pol_e_i, pol_p_i the electron and positron polarizations
     * @param[in] sqrt_s_i the CM energy of the e^+ e^- pair
     */
    eeffRstrange(const StandardModel& SM_i, const double pol_e_i, const double pol_p_i, const double sqrt_s_i) 
    : eeffThObservable(SM_i, pol_e_i, pol_p_i, sqrt_s_i)
    {
    }

    /**
     * @return the ratio @f$R_b=\sigma(e^+e^-\to s\bar{s})/\sigma(e^+e^-\to hadrons)@f$ at sqrt_s
     */
    double computeThValue();

private:
    
};

#endif	/* EEFFSIGMASTRANGE_H */

