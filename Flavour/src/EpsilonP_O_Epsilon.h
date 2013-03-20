/* 
 * Copyright (C) 2012 SusyFit Collaboration
 * All rights reserved.
 *
 * For the licensing terms see doc/COPYING.
 */

#ifndef EPSILONP_O_EPSILON_H
#define	EPSILONP_O_EPSILON_H

#include <ThObservable.h>
#include "Flavour.h"
#include "AmpDS1.h"

using namespace gslpp;

class EpsilonP_O_Epsilon : public ThObservable, AmpDS1 {
public:   
    /**
     * constructor
     * @param Flavour
     */
    EpsilonP_O_Epsilon(Flavour& Flavour): ThObservable(Flavour), AmpDS1(Flavour) {};
    
    /**
     * 
     * @return theoretical value of |\f$ \epsilon ' / \epsilon \f$| 
     */
    double getThValue();
    
private:
    
};

#endif	/* EPSILONP_O_EPSILON_H */
