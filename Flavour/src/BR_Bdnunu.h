/* 
 * Copyright (C) 2012 HEPfit Collaboration
 * All rights reserved.
 *
 * For the licensing terms see doc/COPYING.
 */

#ifndef BR_BDNUNU_H
#define	BR_BDNUNU_H

#include <ThObservable.h>
#include "Flavour.h"
#include <StandardModel.h>

class BR_Bdnunu : public ThObservable {
public:   
    /**
     * constructor
     * @param Flavour
     */
    BR_Bdnunu(StandardModel& SM_i): ThObservable(SM_i)//, mySM(SM_i)
    {};
    
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
    gslpp::complex BRBdnunu(orders order);
    
private:
    
};

#endif	/* BR_BDNUNU_H */