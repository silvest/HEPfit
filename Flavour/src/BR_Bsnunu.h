/* 
 * Copyright (C) 2012 HEPfit Collaboration
 * All rights reserved.
 *
 * For the licensing terms see doc/COPYING.
 */

#ifndef BR_BSNUNU_H
#define	BR_BSNUNU_H

#include <ThObservable.h>
#include "Flavour.h"

class BR_Bsnunu : public ThObservable {
public:   
    /**
     * constructor
     * @param Flavour
     */
    BR_Bsnunu(StandardModel& SM_i): ThObservable(SM_i), mySM(SM_i){};
    
    /**
     * 
     * @brief hep-ph/9512380v2
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
    gslpp::complex BRBsnunu(orders order);
    
private:
    
    StandardModel& mySM;
};

#endif	/* BR_BSNUNU_H */