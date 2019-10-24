/* 
 * Copyright (C) 2012 HEPfit Collaboration
 *
 *
 * For the licensing terms see doc/COPYING.
 */

#ifndef BR_KPPNUNU_H
#define	BR_KPPNUNU_H

class StandardModel;
#include "ThObservable.h"
#include "OrderScheme.h"
#include "gslpp.h"
#include "Charm_Kpnunu.h"


/**
 * @class BR_Kppnunu
 * @ingroup Flavour
 * @brief A class for the branching ratio of \f$K^+\to\pi^+\nu\bar{\nu}\f$
 * @author HEPfit Collaboration
 * @copyright GNU General Public License
 * @details This class is used to compute the theoretical value of
 * the branching ratio of \f$K^+\to\pi^+\nu\bar{\nu}\f$.
 * 
 * 
 *
 * @anchor BR_Kppnunu
 * <h3>%Model parameters</h3>
 *
 * The model parameters of %BR_Kppnunu are summarized below:
 * <table class="model">
 * <tr>
 *   <th>Label</th>
 *   <th>LaTeX symbol</th>
 *   <th>Description</th>
 * </tr>
 * <tr>
 *   <td class="mod_name">%Br_Kp_P0enu</td>
 *   <td class="mod_symb">@f$\mathrm{BR}(K^+\to\pi^0e^+\nu)@f$</td>
 *   <td class="mod_desc">The experimental value for the branching ratio of \f$K^+\to\pi^0e^+\nu\f$.</td>
 * </tr>
 * <tr>
 *   <td class="mod_name">%DeltaP_cu</td>
 *   <td class="mod_symb">@f$@f$</td>
 *   <td class="mod_desc">The long-distance correction to the charm contribution of \f$K^+\to\pi^+\nu\bar{\nu}\f$.</td>
 * </tr>
 * <tr>
 *   <td class="mod_name">%IB_Kp</td>
 *   <td class="mod_symb">@f$@f$</td>
 *   <td class="mod_desc">The isospin breaking corrections between @f$K^+\to\pi^+ \nu\bar{\nu}@f$ and \f$K^+\to\pi^0 e^+\nu\f$.</td>
 * </tr>
 * </table>
 * 
 */
class BR_Kppnunu : public ThObservable {
public:   
    /**
     * constructor
     * @param Flavour
     */
    BR_Kppnunu(StandardModel& SM_i);
    
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
     * @param order_qed
     * @return the short distance contribution to the 
     * |\f$ BR(K^{+} \rightarrow \pi^{+} \nu \bar{\nu}) \f$|, for example
     * see hep-ph/0603079 section 2.3
     */
    gslpp::complex BRKppnunu(orders order, orders_qed order_qed);
    
    /**
     * 
     * @param order
     * @return \f$ P_{C} \f$ defined for exmple in hep-ph/0603079 
     */
    gslpp::complex P_C(orders order);
    
private:
    
    StandardModel& mySM;
    Charm_Kpnunu CKpnunu;
};

#endif	/* BR_KPPNUNU_H */
