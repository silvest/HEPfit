/* 
 * Copyright (C) 2012 HEPfit Collaboration
 *
 *
 * For the licensing terms see doc/COPYING.
 */

#ifndef DMK_H
#define DMK_H

#include "ThObservable.h"
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
 * @f$\Delta M_K@f$. The Standard Model contribution can only be estimated 
 * non-perturbatively and therefore it is implemented as an
 * input parameter in this code through SM.getDMK().
 */
class DmK : public ThObservable, AmpDK2 {
public:

    /**
    * @brief Constructor.
    * @param[in] SM_i a reference to an object of type StandardModel
    */
    DmK(const StandardModel& SM_i);

    /**
     * 
     * @return theoretical value of @f$\Delta M_K@f$ 
     */
    virtual double computeThValue();
};

#endif /* DMK_H */
