/* 
 * Copyright (C) 2012 SusyFit Collaboration
 * All rights reserved.
 *
 * For the licensing terms see doc/COPYING.
 */

#ifndef BR_BSNUNU_H
#define	BR_BSNUNU_H

#include <ThObservable.h>
#include "Flavour.h"

using namespace gslpp;

class BR_Bsnunu : public ThObservable {
public:   
    /**
     * constructor
     * @param Flavour
     */
    BR_Bsnunu(Flavour& Flavour): ThObservable(Flavour), myFlavour(Flavour){};
    
    /**
     * 
     * @brief hep-ph/9512380v2
     * @return theoretical value of |\f$ BR(B_s \rightarrow \nu \bar{\nu}) \f$|
     */
    double getThValue();
    
    
protected:
    
    /**
     * 
     * @param order
     * @param order_ew
     * @return the short distance contribution to the 
     * |\f$ BR(B_s \rightarrow \nu \bar{\nu}) \f$|
     */
    complex BRBsnunu(orders order);
    
private:
    Flavour& myFlavour;
};

#endif	/* BR_BSNUNU_H */