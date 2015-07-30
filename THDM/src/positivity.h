/* 
 * Copyright (C) 2015 SusyFit Collaboration
 * All rights reserved.
 *
 * For the licensing terms see doc/COPYING.
 */

#ifndef POSITIVITY_H
#define	POSITIVITY_H

#include <stdexcept>
#include <ThObservable.h>
#include "THDM.h"
#include "lambda1.h"
#include "lambda2.h"

/**
 * @class positivity
 * @ingroup THDM 
 * @brief An observable class for the positivity conditions of the Higgs potential.
 * @author SusyFit Collaboration
 * @copyright GNU General Public License
 * @details This class is used to compute the positivity conditions of the Higgs potential.
 */
class positivity : public ThObservable {
public:
    /**
     * @brief Constructor.
     * @param[in] ?
     */
   positivity(const StandardModel& SM_i, int obsFlag);
     
   ~positivity();

   /**
     * @brief Two positivity conditions of the Higgs potential.
     * @return
     */
    double computeThValue();

    
    private:
        const THDM * myTHDM;
        int obs;
        lambda1 * mylambda1;
        lambda2 * mylambda2;
};

#endif	/* POSITIVITY_H */



