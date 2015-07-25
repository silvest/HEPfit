/* 
 * Copyright (C) 2012 HEPfit Collaboration
 * All rights reserved.
 *
 * For the licensing terms see doc/COPYING.
 */

#ifndef DMK_H
#define	DMK_H

#include <ThObservable.h>
#include "Flavour.h"
#include "AmpDK2.h"

/**
 * @class DMK
 * @ingroup Flavour
 * @brief A class for @f$\Delta M_K@f$, the difference in mass
 * of between the heavy and light mass eigenstates of the neutral
 * Kaon system.
 * @author HEPfit Collaboration
 * @copyright GNU General Public License
 * @details This class is used to compute the theoretical value of
 * @f$\Delta M_K@f$. While the short distance contributions to this 
 * parameter is calculable both in the SM and in NP models, the long distance
 * contribution can only be estimated and might contribute to as much as half the
 * value of the parameter. This long distance contribution is put in as an
 * input parameter in this code through SM.getDMK().
 */
class DmK : public ThObservable, AmpDK2 {
public:
    
    /**
     * constructor
     * @param Flavour
     */
    
    DmK(const StandardModel& SM_i) : ThObservable(SM_i), AmpDK2(SM_i) {};
    
    /**
     * 
     * @return theoretical value of @f$\Delta M_K@f$ 
     */
    virtual double computeThValue();
};

#endif	/* DMK_H */
