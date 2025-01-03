/* 
 * Copyright (C) 2012 HEPfit Collaboration
 *
 *
 * For the licensing terms see doc/COPYING.
 */

#ifndef LEP2AFBBOTTOM_H
#define	LEP2AFBBOTTOM_H

#include "LEP2ThObservable.h"

/**
 * @class LEP2AFBbottom
 * @ingroup EW
 * @brief A class for the forward-backward asymmetry of @f$e^+e^-\to b\bar{b}@f$ 
 * above the @f$Z@f$ pole.
 * @author HEPfit Collaboration
 * @copyright GNU General Public License
 * @details 
 */
class LEP2AFBbottom : public LEP2ThObservable {
public:

    /**
     * @brief LEP2AFBbottom constructor
     * @param[in] SM_i a reference to an object of type StandardModel
     * @param[in] sqrt_s_i the CM energy of the e^+ e^- pair
     */
    LEP2AFBbottom(const StandardModel& SM_i, const double sqrt_s_i) 
    : LEP2ThObservable(SM_i, sqrt_s_i) 
    {
    }

    /**
     * @return the forward-backward asymmetry for e^+ e^- -> b bbar at sqrt_s
     */
    double computeThValue();

private:
     
};


/**
 * @class eeffAFBbottom
 * @ingroup EW
 * @brief A class for the forward-backward asymmetry of @f$e^+e^-\to b\bar{b}@f$ 
 * above the @f$Z@f$ pole.
 * @author HEPfit Collaboration
 * @copyright GNU General Public License
 * @details 
 */
class eeffAFBbottom : public eeffThObservable {
public:

    /**
     * @brief eeffAFBbottom constructor
     * @param[in] pol_e_i, pol_p_i the electron and positron polarizations
     * @param[in] sqrt_s_i the CM energy of the e^+ e^- pair
     */
    eeffAFBbottom(const StandardModel& SM_i, const double pol_e_i, const double pol_p_i, const double sqrt_s_i) 
    : eeffThObservable(SM_i, pol_e_i, pol_p_i, sqrt_s_i)
    {
    }

    /**
     * @return the forward-backward asymmetry for e^+ e^- -> b bbar at sqrt_s
     */
    double computeThValue();

private:
     
};

#endif	/* LEP2AFBBOTTOM_H */

