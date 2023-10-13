/* 
 * Copyright (C) 2012 HEPfit Collaboration
 *
 *
 * For the licensing terms see doc/COPYING.
 */

#ifndef NPBR_KPNUNU_H
#define NPBR_KPNUNU_H

class StandardModel;
#include "ThObservable.h"
#include "OrderScheme.h"
#include "gslpp.h"
#include "Charm_Kpnunu.h"
#include "StandardModel.h"
#include "StandardModelMatching.h"
#include "EpsilonK.h"
#include "BR_Kppnunu.h"
#include "BR_Kp0nunu.h"
#include "WilsonCoefficient.h"


class NPBR_kp0nunu : public ThObservable {
public:
    /**
     * constructor
     * @param Flavour
     */
    NPBR_kp0nunu(const StandardModel& SM_i);
    
    /**
     * 
     * @return theoretical value of \f$ BR(K_L \rightarrow \pi^{0} \nu \bar{\nu}) \f$, 
     * for example see hep-ph/0603079 section 2.3
     */
    double computeThValue();
    
    /**
     * 
     * @return The full contribution w/o approximation at highest order. 
     * Here Xt is defined by: Xt_SM * R_K * exp(-i*Theta_K).  
     * (see hep-ph/9607447v1 eq. 9 for the SM formula) 
     */
    double BRKp0nunu_NP(orders order, orders_qed order_qed);
    
    
    

private:
    const StandardModel& mySM;
    BR_Kp0nunu myBR_Kp0nunu;

};


class NPBR_kppnunu : public ThObservable {
public:
    /**
     * constructor
     * @param Flavour
     */
    NPBR_kppnunu(const StandardModel& SM_i);
    
    /**
     * 
     * @return theoretical value of |\f$ BR(K^{+} \rightarrow \pi^+ \nu \bar{\nu}) \f$| with new physics, 
     * See arxiv: 0705.2025v2 eq. 7 for SM formula
     */
    double computeThValue();

    
    /**
     * 
     * @param order
     * @param order_qed
     * @return the short distance contribution to the 
     * |\f$ BR(K^{+} \rightarrow \pi^{+} \nu \bar{\nu}) with New Physics contribution:
     * Xt is defined by: Xt_SM * R_K * exp(-i*Theta_K). \f$|
     */
    double BRKppnunu_NP(orders order, orders_qed order_qed);
    

private:
    const StandardModel& mySM;
    BR_Kppnunu myBR_Kppnunu;
};



#endif /* NPBR_KPNUNU_H */

