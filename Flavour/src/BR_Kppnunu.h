/* 
 * Copyright (C) 2012 SusyFit Collaboration
 * All rights reserved.
 *
 * For the licensing terms see doc/COPYING.
 */

#ifndef BR_KPPNUNU_H
#define	BR_KPPNUNU_H

#include <ThObservable.h>
#include "Flavour.h"
#include "Charm_Kpnunu.h"

using namespace gslpp;

class BR_Kppnunu : public ThObservable {
public:   
    /**
     * constructor
     * @param Flavour
     */
    BR_Kppnunu(Flavour& Flavour): ThObservable(Flavour), myFlavour(Flavour), CKpnunu(SM) {};
    
    /**
     * 
     * @return theoretical value of |\f$ BR(K^{+} \rightarrow \pi^+ \nu \bar{\nu}) \f$|, 
     * for example see hep-ph/0603079 section 2.3
     */
    double computeThValue();
    
    
protected:
    
    /**
     * 
     * @param order
     * @param order_ew
     * @return the short distance contribution to the 
     * |\f$ BR(K^{+} \rightarrow \pi^{+} \nu \bar{\nu}) \f$|, for example
     * see hep-ph/0603079 section 2.3
     */
    complex BRKppnunu(orders order, orders_ew order_ew);
    
    /**
     * 
     * @param order
     * @return \f$ P_{C} \f$ defined for exmple in hep-ph/0603079 
     */
    complex P_C(orders order);
    
private:
    Flavour& myFlavour;
    Charm_Kpnunu CKpnunu;
};

#endif	/* BR_KPPNUNU_H */
