/*
 * Copyright (C) 2012 HEPfit Collaboration
 * All rights reserved.
 *
 * For the licensing terms see doc/COPYING.
 */

#ifndef DMBD_H
#define	DMBD_H

#include <ThObservable.h>
#include "Flavour.h"
#include "AmpDB2.h"

/**
 * @class DmBd
 * @ingroup Flavour
 * @brief A class for @f$\Delta M_{B_d}@f$, the difference in mass
 * of between the heavy and light mass eigenstates of the neutral
 * @f$B_d@f$ system.
 * @author HEPfit Collaboration
 * @copyright GNU General Public License
 * @details This class is used to compute the theoretical value of
 * @f$\Delta M_{B_d}@f$. The short distance contributions to this
 * parameter is calculable both in the SM and in NP models. In the SM
 * the short distance contribution is dominated by the top loop in the
 * box diagram.
 */
class DmBd : public ThObservable, AmpDB2 {
public:
    
    /**
     * constructor
     * @param Flavour
     */
    
    DmBd(const StandardModel& SM_i) : ThObservable(SM_i), AmpDB2(SM_i) {};
    
    /**
     *
     * @return theoretical value of @f$\Delta M_{B_d}@f$
     */
    
    double computeThValue();
};

#endif	/* DMBD_H */
