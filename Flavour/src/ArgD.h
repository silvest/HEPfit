/* 
 * Copyright (C) 2012 HEPfit Collaboration
 *
 *
 * For the licensing terms see doc/COPYING.
 */

#ifndef ARGD_H
#define	ARGD_H

#include "ThObservable.h"
#include "gslpp.h"
#include "AmpDD2.h"
#include "std_make_vector.h"

class ArgD : public ThObservable, AmpDD2 {
/**
 * @brief a class for the argument of the complex amplitude of D^{0} oscillations
 * @param Flavour an object of Flavour class
 */
public:
    /**
     * @brief ArgD constructor
     * @param Flavour an object of Flavour class
     */
    ArgD(const StandardModel& SM_i) : ThObservable(SM_i), AmpDD2(SM_i) 
    {
        setParametersForObservable(make_vector<std::string>() << "SM_M12D");
    };
    /**
     * @brief a method returning the argument of the complex amplitude for 
     * the absorptive part of the \f$ | \Delta C = 2 | \f$ mixing
     * @return the argument of the complex amplitude of \f$ D^{0} \f$ oscillations 
     */
    double computeThValue() {   
       return(AmpDD(NLO).arg() / M_PI * 180. );
    };
};

#endif	/* ARGD_H */
