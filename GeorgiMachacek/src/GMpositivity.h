/*
 * Copyright (C) 2017 HEPfit Collaboration
 *
 *
 * For the licensing terms see doc/COPYING.
 */

#ifndef GMPOSITIVITY_H
#define	GMPOSITIVITY_H

#include "ThObservable.h"
#include "GeorgiMachacek.h"

/**
 * @class GMpositivity
 * @ingroup GeorgiMachacek 
 * @brief Base class for "boundedness-from-below" constraints
 * @author HEPfit Collaboration
 * @copyright GNU General Public License
 * @details Formulae taken from @cite Chiang/Kuo/Yamada.
 */

/**
 * @class GMpositivity1
 * @brief Controls that the scalar %THDMW potential is bounded from below.
 * @details @f$\lambda_2+\lambda_3>0@f$.
 * @author HEPfit Collaboration
 * @copyright GNU General Public License
 */
class GMpositivity1: public ThObservable {
public:

    /**
     * @brief GMpositivity1 constructor.
     */
    GMpositivity1(const StandardModel& SM_i);

    /**
     * @return @f$\lambda_2+\lambda_3@f$
     */
    double computeThValue();
private:
    const GeorgiMachacek& myGM;
};

#endif	/* GMPOSITIVITY_H */
