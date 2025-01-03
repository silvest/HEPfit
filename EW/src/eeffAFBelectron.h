/* 
 * Copyright (C) 2012 HEPfit Collaboration
 *
 *
 * For the licensing terms see doc/COPYING.
 */

#ifndef EEFFAFBELECTRON_H
#define	EEFFAFBELECTRON_H

#include "LEP2ThObservable.h"

/**
 * @class eeffAFBelectron
 * @ingroup EW
 * @brief A class for the forward-backward asymmetry of @f$e^+e^-\to e^+e^-@f$ 
 * above the @f$Z@f$ pole.
 * @author HEPfit Collaboration
 * @copyright GNU General Public License
 * @details 
 */
class eeffAFBelectron : public eeffThObservable {
public:

    /**
     * @brief eeffAFBelectron constructor
     * @param[in] pol_e_i, pol_p_i the electron and positron polarizations
     * @param[in] sqrt_s_i the CM energy of the e^+ e^- pair
     */
    eeffAFBelectron(const StandardModel& SM_i, const double pol_e_i, const double pol_p_i, const double sqrt_s_i) 
    : eeffThObservable(SM_i, pol_e_i, pol_p_i, sqrt_s_i)
    {
    }
    
    /**
     * @return the forward-backward asymmetry for e^+ e^- -> e^+ e^- at sqrt_s
     */
    double computeThValue();

private:
    
};

#endif	/* EEFFAFBELECTRON_H */

