/* 
 * Copyright (C) 2012 HEPfit Collaboration
 * All rights reserved.
 *
 * For the licensing terms see doc/COPYING.
 */

#ifndef AMPDD2_H
#define	AMPDD2_H

#include <gslpp_complex.h>
#include "Flavour.h"
#include <StandardModel.h>

class AmpDD2 {
/**
 * @brief \f$ | \Delta C = 2 | \f$ Amplitude Class
 * @details This class is related to the calculation of the \f$ D^{0}-\bar{D}^{0}\f$
 * mixing amplitude in relation to the absorptive part of the mixing;  
 * it is written to consider also the contribution coming from supersymmetric theories
 * in the squark-gluino-quark sector in the mass insertion approximation.
 * @param Flavour an object of Flavour class
 */
public:
/**
 * @brief AmpDD2 constructor  
 * @param Flavour an object of Flavour class
 */
    AmpDD2(const StandardModel& SM_i);
protected:
/**
 * @brief a method returning the complex amplitude for the \f$ D^{0}-\bar{D}^{0}\f$ mixing
 * @param an enum "orders" for the order of perturbation theory of the evolutor 
 * @return the complex \f$ | \Delta C = 2 | \f$ amplitude 
 */
    gslpp::complex AmpDD(orders order);
    
private:
    
    const StandardModel& mySM;
};

#endif	/* AMPDD2_H */

