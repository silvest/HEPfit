/* 
 * Copyright (C) 2015 HEPfit Collaboration
 * All rights reserved.
 *
 * For the licensing terms see doc/COPYING.
 */

#ifndef POSITIVITY_H
#define	POSITIVITY_H

#include <stdexcept>
#include "ThObservable.h"
#include "THDM.h"
#include "THDMquantities.h"

/**
 * @class positivity
 * @ingroup THDM 
 * @brief An observable class for the positivity conditions of the Higgs potential.
 * @author HEPfit Collaboration
 * @copyright GNU General Public License
 * @details This class is used to compute the positivity conditions of the Higgs potential.
 * Formulae can be found in @cite Deshpande:1977rw.
 */
class positivity : public ThObservable {
public:
    /**
     * @brief Constructor.
     * @param[in] ?
     */
   positivity(const StandardModel& SM_i);
     
   ~positivity();

   /**
     * @brief Two positivity conditions of the Higgs potential.
     * @return
     */
    double computeThValue();
    const THDM * myTHDM;
    lambda1 * mylambda1;
    lambda2 * mylambda2;
};

class positivity1: public positivity {
public:

    /**
     * @brief Constructor.
     */
    positivity1(const StandardModel& SM_i);

    /**
     * @return positivity1
     */
    double computeThValue();
};

class positivity2: public positivity {
public:

    /**
     * @brief Constructor.
     */
    positivity2(const StandardModel& SM_i);

    /**
     * @return positivity2
     */
    double computeThValue();
};

#endif	/* POSITIVITY_H */
