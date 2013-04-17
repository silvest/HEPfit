/* 
 * Copyright (C) 2012 SusyFit Collaboration
 * All rights reserved.
 *
 * For the licensing terms see doc/COPYING.
 */

#ifndef BR_BDMUMU_H
#define	BR_BDMUMU_H

#include <ThObservable.h>
#include "Flavour.h"

using namespace gslpp;

class BR_Bdmumu : public ThObservable {
public:   
    /**
     * constructor
     * @param Flavour
     */
    BR_Bdmumu(Flavour& Flavour): ThObservable(Flavour), myFlavour(Flavour){};
    
    /**
     * 
     * @return theoretical value of |\f$ BR(B_d \rightarrow \mu \bar{\mu}) \f$|
     */
    double getThValue();
    
    
protected:
    
    /**
     * 
     * @param order
     * @param order_ew
     * @return the short distance contribution to the 
     * |\f$ BR(B_d \rightarrow \mu \bar{\mu}) \f$|
     */
    complex BRBdmumu(orders order);
    
private:
    Flavour& myFlavour;
};

#endif	/* BR_BDMUMU_H */