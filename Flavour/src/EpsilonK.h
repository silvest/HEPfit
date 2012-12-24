/* 
 * Copyright (C) 2012 SusyFit Collaboration
 * All rights reserved.
 *
 * For the licensing terms see doc/COPYING.
 */

#ifndef EPSILONK_H
#define	EPSILONK_H

#include <ThObservable.h>
#include "Flavour.h"
#include "AmpDK2.h"

using namespace gslpp;

class EpsilonK : public ThObservable, AmpDK2 {
public:   
    /**
     * constructor
     * @param Flavour
     */
    EpsilonK(Flavour& Flavour): ThObservable(Flavour), AmpDK2(Flavour) {};
    
    /**
     * 
     * @return theoretical value of Epsilon_K 
     */
    double getThValue();
    
};

#endif	/* EPSILONK_H */



