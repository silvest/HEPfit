/* 
 * Copyright (C) 2015 HEPfit Collaboration
 *
 *
 * For the licensing terms see doc/COPYING.
 */

#ifndef POSITIVITY_H
#define	POSITIVITY_H

#include "ThObservable.h"


class THDM;
class lambda1;
class lambda2;

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

/**
 * @class positivity1
 * @ingroup THDM 
 * @brief Controls that the scalar %THDM potential is bounded from below.
 * @details @f$\lambda_3>-\sqrt{\lambda_1 \lambda_2}@f$.
 */
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

/**
 * @class positivity2
 * @ingroup THDM 
 * @brief Controls that the scalar %THDM potential is bounded from below.
 * @details @f$\lambda_3+\lambda_4-|\lambda_5|>-\sqrt{\lambda_1 \lambda_2}@f$.
 */
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
