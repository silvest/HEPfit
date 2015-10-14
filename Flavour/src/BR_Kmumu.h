/* 
 * Copyright (C) 2012 HEPfit Collaboration
 * All rights reserved.
 *
 * For the licensing terms see doc/COPYING.
 */

#ifndef BR_KMUMU_H
#define	BR_KMUMU_H

#include <ThObservable.h>
#include "Flavour.h"
#include "CPenguinBoxMu.h"
#include <StandardModel.h>

class BR_Kmumu : public ThObservable {
public:   
    /**
     * constructor
     * @param Flavour
     */
    BR_Kmumu(StandardModel& SM_i): ThObservable(SM_i), mySM(SM_i), CPB(SM){};
    
    /**
     * 
     * @return theoretical value of |\f$ BR(K_L \rightarrow \mu \bar{\mu}) \f$|, 
     * for example see hep-ph/0603079 section 2.3
     */
    double computeThValue();
    
    
protected:
    
    /**
     * 
     * @param order
     * @param order_ew
     * @return the short distance contribution to the 
     * |\f$ BR(K_L \rightarrow \mu \bar{\mu}) \f$|, 
     */
    gslpp::complex BRKmumu(orders order);
    
private:
    
    StandardModel& mySM;
    CPenguinBoxMu CPB;
};

#endif	/* BR_KMUMU_H */