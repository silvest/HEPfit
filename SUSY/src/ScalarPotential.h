/* 
 * Copyright (C) 2016 HEPfit Collaboration
 *
 *
 * For the licensing terms see doc/COPYING.
 */

#ifndef SCALARPOTENTIAL_H
#define SCALARPOTENTIAL_H

#include "SUSY.h"

/**
 * @class SUSYScalarPotential
 * @ingroup SUSY
 * @brief SUSYScalarPotential.
 */
class SUSYScalarPotential//: public ThObservable
{
public:

    /**
     * @brief SUSYScalarPotential constructor.
     */
    SUSYScalarPotential(const StandardModel& SM_i);

    /**
     * @return @f$SUSYScalarPotential@f$
     */
    gslpp::vector<double> coefficients();

    /**
     * @return @f$V_{scalar}@f$
     */
    double potential(gslpp::vector<double> coefficients, double field1, double field2, double field3);

    /**
     * @return @f$dV_{scalar}@f$
     */
    gslpp::vector<double> potentialderivative(gslpp::vector<double> coefficients, double field1, double field2, double field3);

    /**
     * @return @f$d^2V_{scalar}@f$
     */
    gslpp::vector<double> secondpotentialderivative(gslpp::vector<double> coefficients, double field1, double field2, double field3);

private:
    const SUSY& mySUSY;
};

#endif /* SCALARPOTENTIAL_H */

