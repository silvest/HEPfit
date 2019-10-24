/* 
 * Copyright (C) 2012 HEPfit Collaboration
 *
 *
 * For the licensing terms see doc/COPYING.
 */

#ifndef BR_KP0NUNU_H
#define	BR_KP0NUNU_H

class StandardModel;
#include "ThObservable.h"
#include "OrderScheme.h"
#include "gslpp.h"



/**
 * @class BR_Kp0nunu
 * @ingroup Flavour
 * @brief A class for the branching ratio of \f$K_L\to\pi^0\nu\bar{\nu}\f$
 * @author HEPfit Collaboration
 * @copyright GNU General Public License
 * @details This class is used to compute the theoretical value of
 * the branching ratio of \f$K_L\to\pi^0\nu\bar{\nu}\f$.
 * 
 * 
 *
 * @anchor BR_Kp0nunu
 * <h3>%Model parameters</h3>
 *
 * The model parameters of %BR_Kp0nunu are summarized below:
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
 *   <td class="mod_name">%IB_Kl</td>
 *   <td class="mod_symb">@f$@f$</td>
 *   <td class="mod_desc">the isospin breaking corrections between @f$K_L\to\pi^0\nu\bar{\nu}@f$ and \f$K^+\to\pi^0 e^+\nu\f$.</td>
 * </tr>
 * </table>
 * 
 */
class BR_Kp0nunu : public ThObservable {
public:   
    /**
     * constructor
     * @param Flavour
     */
    BR_Kp0nunu(StandardModel& SM_i);
    
    /**
     * 
     * @return theoretical value of \f$ BR(K_L \rightarrow \pi^{0} \nu \bar{\nu}) \f$, 
     * for example see hep-ph/0603079 section 2.3
     */
    double computeThValue();
    
    
protected:
    
    /**
     * 
     * @param order
     * @param order_qed
     * @return the short distance contribution to the 
     * \f$ BR(K_{L} \rightarrow \pi^{0} \nu \bar{\nu}) \f$, for example
     * see hep-ph/0603079 section 2.3
     */
    gslpp::complex BRKp0nunu(orders order, orders_qed order_qed);
    
private:
    
    StandardModel& mySM;
};

#endif	/* BR_KP0NUNU_H */
