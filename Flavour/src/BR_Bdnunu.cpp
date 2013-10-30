/* 
 * Copyright (C) 2012 SusyFit Collaboration
 * All rights reserved.
 *
 * For the licensing terms see doc/COPYING.
 */

#ifndef BR_BDNUNU_H
#define	BR_BDNUNU_H

#include <ThObservable.h>
#include "Flavour.h"

using namespace gslpp;

class BR_Bdnunu : public ThObservable {
public:   
    /**
     * constructor
     * @param Flavour
     */
    BR_Bdnunu(Flavour& Flavour): ThObservable(Flavour), myFlavour(Flavour){};
    
    /**
     * 
     * @return theoretical value of |\f$ BR(B_s \rightarrow \nu \bar{\nu}) \f$|
     */
    double computeThValue();
    
    
protected:
    
    /**
     * 
     * @param order
     * @param order_ew
     * @return the short distance contribution to the 
     * |\f$ BR(B_s \rightarrow \nu \bar{\nu}) \f$|
     */
    complex BRBdnunu(orders order);
    
private:
    Flavour& myFlavour;
};

#endif	/* BR_BDNUNU_H */