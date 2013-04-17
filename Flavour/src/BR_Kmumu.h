/* 
 * Copyright (C) 2012 SusyFit Collaboration
 * All rights reserved.
 *
 * For the licensing terms see doc/COPYING.
 */

#ifndef BR_KMUMU_H
#define	BR_KMUMU_H

#include <ThObservable.h>
#include "Flavour.h"
#include "CPenguinBoxMu.h"

using namespace gslpp;

class BR_Kmumu : public ThObservable {
public:   
    /**
     * constructor
     * @param Flavour
     */
    BR_Kmumu(Flavour& Flavour): ThObservable(Flavour), myFlavour(Flavour), CPB(SM){};
    
    /**
     * 
     * @return theoretical value of |\f$ BR(K_L \rightarrow \mu \bar{\mu}) \f$|, 
     * for example see hep-ph/0603079 section 2.3
     */
    double getThValue();
    
    
protected:
    
    /**
     * 
     * @param order
     * @param order_ew
     * @return the short distance contribution to the 
     * |\f$ BR(K_L \rightarrow \mu \bar{\mu}) \f$|, 
     */
    complex BRKmumu(orders order);
    
private:
    Flavour& myFlavour;
    CPenguinBoxMu CPB;
};

#endif	/* BR_KMUMU_H */