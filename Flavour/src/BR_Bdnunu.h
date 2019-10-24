/* 
 * Copyright (C) 2012 HEPfit Collaboration
 *
 *
 * For the licensing terms see doc/COPYING.
 */

#ifndef BR_BDNUNU_H
#define	BR_BDNUNU_H

#include "ThObservable.h"
#include "OrderScheme.h"
#include "gslpp.h"
class StandardModel;

class BR_Bdnunu : public ThObservable {
public:   
    /**
    * @brief Constructor.
    * @param[in] SM_i a reference to an object of type StandardModel
    */
    BR_Bdnunu(StandardModel& SM_i): ThObservable(SM_i)//, mySM(SM_i)
    {};
    
    /**
     * 
     * @return theoretical value of |\f$ BR(B_d \rightarrow \nu \bar{\nu}) \f$|
     */
    double computeThValue();
    
    
protected:
    
    /**
    * @brief the short distance contribution to the 
    * |\f$ BR(B_d \rightarrow \nu \bar{\nu}) \f$|
    * @param[in] order the %QCD order of the computation
    */
    gslpp::complex BRBdnunu(orders order);
    
private:
    
};

#endif	/* BR_BDNUNU_H */