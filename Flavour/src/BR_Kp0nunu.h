/* 
 * Copyright (C) 2012 SusyFit Collaboration
 * All rights reserved.
 *
 * For the licensing terms see doc/COPYING.
 */

#ifndef BR_KP0NUNU_H
#define	BR_KP0NUNU_H

#include <ThObservable.h>
#include "Flavour.h"

using namespace gslpp;

class BR_Kp0nunu : public ThObservable {
public:   
    /**
     * constructor
     * @param Flavour
     */
    BR_Kp0nunu(Flavour& Flavour): ThObservable(Flavour), myFlavour(Flavour){};
    
    /**
     * 
     * @return theoretical value of |\f$ BR(K_L \rightarrow \pi^0 \nu \bar{\nu}) \f$|, 
     * for example see hep-ph/0603079 section 2.3
     */
    double getThValue();
    
    
protected:
    
    /**
     * 
     * @param order
     * @param order_ew
     * @return the short distance contribution to the 
     * |\f$ BR(K_L \rightarrow \pi^0 \nu \bar{\nu}) \f$|, for example
     * see hep-ph/0603079 section 2.3
     */
    complex BRKp0nunu(orders order, orders_ew order_ew);
    
private:
    Flavour& myFlavour;
};

#endif	/* BR_KP0NUNU_H */
