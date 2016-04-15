/* 
 * Copyright (C) 2016 HEPfit Collaboration
 *
 *
 * For the licensing terms see doc/COPYING.
 */

#ifndef SFERMIONDELTAS_H
#define	SFERMIONDELTAS_H

//#include <stdexcept>
#include "ThObservable.h"
#include "SUSY.h"

/**
 * @class deltaRL_12_u
 * @ingroup LeptonFlavour
 * @brief deltaRL_12_u.
 */
class deltaRL_12_u: public ThObservable {
public:

    /**
     * @brief deltaRL_12_u constructor.
     */
    deltaRL_12_u(const StandardModel& SM_i);

    /**
     * @return @f$\delta_RL12u@f$
     */
    double computeThValue();

private:
    const SUSY& mySUSY;
};

#endif	/* SFERMIONDELTAS_H */
