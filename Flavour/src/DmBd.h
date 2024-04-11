/*
 * Copyright (C) 2012 HEPfit Collaboration
 *
 *
 * For the licensing terms see doc/COPYING.
 */

#ifndef DMBD_H
#define	DMBD_H

#include "ThObservable.h"

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
class DmBd : public ThObservable {
public:
    
    /**
    * @brief Constructor.
    * @param[in] SM_i a reference to an object of type StandardModel
    */
    DmBd(const StandardModel& SM_i) : ThObservable(SM_i){};
    
    /**
     *
     * @return theoretical value of @f$\Delta M_{B_d}@f$
     */
    
    double computeThValue();
};

#endif	/* DMBD_H */
