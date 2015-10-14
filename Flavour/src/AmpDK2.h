/* 
 * Copyright (C) 2012 HEPfit Collaboration
 * All rights reserved.
 *
 * For the licensing terms see doc/COPYING.
 */

#ifndef AMPDK2_H
#define	AMPDK2_H

#include <gslpp.h>
#include "Flavour.h"
#include <StandardModel.h>

/**
 * @class AmpDK2
 * @ingroup Flavour
 * @brief A class for calculating the amplitudes contributing to @f$\epsilon_K@f$
 * and @f$\Delta M_K@f$.
 * @author HEPfit Collaboration
 * @copyright GNU General Public License
 * @details This class is used to compute the theoretical value of
 * amplitudes contributing to @f$\epsilon_K@f$ and @f$\Delta M_K@f$. The
 * hadronic matrix elements are defined for the operators @f$O_1, \ldots, O_5@f$
 * in the chiral limit as can be found in ...
 *
 */
class AmpDK2 {
public:
    /**
     * 
     * @brief compute the amplitude for kaon oscillations
     * @param[in] SM_i
     */
    AmpDK2(const StandardModel& SM_i);

protected:
    gslpp::complex AmpDK(orders order);
    gslpp::complex AmpMK(orders order);
    
private:
    
    const StandardModel& mySM;

};

#endif	/* AMPDK2_H */


