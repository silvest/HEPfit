/* 
 * Copyright (C) 2015 SusyFit Collaboration
 * All rights reserved.
 *
 * For the licensing terms see doc/COPYING.
 */

#ifndef EWPO_H
#define	EWPO_H

#include <stdexcept>
#include <ThObservable.h>
#include "THDM.h"

/**
 * @class EWPO
 * @ingroup THDM 
 * @brief An observable class to calculate the electroweak precision observables.
 * @author SusyFit Collaboration
 * @copyright GNU General Public License
 * @details An observable class to calculate ??? in context of the 2HDM.
 */
class EWPO : public ThObservable {
public:

    /**
     * @brief Constructor.
     * @param[in] ?
     */
    EWPO(const StandardModel& SM_i, int obsFlag);

    /**
     * @brief EWPO.
     * @return 
     */
    double computeThValue();

    private:
        const THDM * myTHDM;
        int obs;
};

#endif	/* EWPO_H */
