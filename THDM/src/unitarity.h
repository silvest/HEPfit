/* 
 * Copyright (C) 2015 SusyFit Collaboration
 * All rights reserved.
 *
 * For the licensing terms see doc/COPYING.
 */

#ifndef UNITARITY_H
#define	UNITARITY_H

#include <stdexcept>
#include <ThObservable.h>
#include "THDM.h"

/**
 * @class unitarity
 * @ingroup THDM 
 * @brief An observable class for the requirement of tree level perturbative 
 * unitarity.
 * @author SusyFit Collaboration
 * @copyright GNU General Public License
 * @details This class is used to require the unitarity for all the tree level 
 * scalar-scalar scattering amplitudes.
 */
class unitarity : public ThObservable {
public:
    /**
     * @brief Constructor.
     * @param[in] ?
     */
   unitarity(const StandardModel& SM_i, int obsFlag);
     
    /**
     * @brief The unitarity conditions for all the tree level scalar-scalar 
     * scattering amplitudes.
     * @return
     */
    double computeThValue();

    
    private:
        const THDM * myTHDM;
        int obs;
};

#endif	/* UNITARITY_H */

