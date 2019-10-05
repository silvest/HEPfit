/* 
 * Copyright (C) 2012 HEPfit Collaboration
 *
 *
 * For the licensing terms see doc/COPYING.
 */

#ifndef M12D_H
#define	M12D_H

#include "ThObservable.h"
#include "AmpDD2.h"
#include "std_make_vector.h"

class M12D : public ThObservable, AmpDD2 {
/**
 * @brief a class for the absolute value of the complex amplitude of D^{0} oscillations
 * @param Flavour an object of Flavour class
 */
public:
    /**
     * @brief M12D constructor
     * @param Flavour an object of Flavour class
     */
    M12D(const StandardModel& SM_i) : ThObservable(SM_i), AmpDD2(SM_i) 
    {
        setParametersForObservable(make_vector<std::string>() << "SM_M12D");
    };
    /**
     * @brief a method returning the absolute value of the complex amplitude for 
     * the absorptive part of the\f$ | \Delta C = 2 | \f$ mixing
     * @return the absolute value of the complex amplitude of \f$ D^{0} \f$ oscillations
     */
    double computeThValue() {
        return((AmpDD(FULLNLO) + SM.getOptionalParameter("SM_M12D")).abs());
         
    };
};

#endif	/* MD12_H */

