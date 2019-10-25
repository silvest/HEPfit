/* 
 * Copyright (C) 2012 HEPfit Collaboration
 *
 *
 * For the licensing terms see doc/COPYING.
 */

#ifndef BR_KMUMU_H
#define	BR_KMUMU_H

#include "ThObservable.h"
#include "OrderScheme.h"
#include "gslpp.h"
#include "CPenguinBoxMu.h"
class StandardModel;


/**
 * @class BR_Kmumu
 * @ingroup Flavour
 * @brief A class for the branching ratio of \f$K\to\mu^+\mu^-\f$
 * @author HEPfit Collaboration
 * @copyright GNU General Public License
 * @details This class is used to compute the theoretical value of
 * the branching ratio of \f$K\to\mu^+\mu^-\f$.
 * 
 * 
 *
 * @anchor BR_Kmumu
 * <h3>%Model parameters</h3>
 *
 * The model parameters of %BR_Kmumu are summarized below:
 * <table class="model">
 * <tr>
 *   <th>Label</th>
 *   <th>LaTeX symbol</th>
 *   <th>Description</th>
 * </tr>
 * <tr>
 *   <td class="mod_name">%Br_Kp_munu</td>
 *   <td class="mod_symb">@f$\mathrm{BR}(K^+\to\mu^+\nu)@f$</td>
 *   <td class="mod_desc">The experimental value for the branching ratio of \f$K^+\to\mu^+\nu\f$.</td>
 * </tr>
 * <tr>
 *   <td class="mod_name">%DeltaP_cu</td>
 *   <td class="mod_symb">@f$@f$</td>
 *   <td class="mod_desc">The long-distance correction to the charm contribution of \f$K^+\to\pi^+\nu\bar{\nu}\f$.</td>
 * </tr>
 * </table>
 * 
 */
class BR_Kmumu : public ThObservable {
public:   
    /**
     * constructor
     * @param Flavour
     */
    BR_Kmumu(StandardModel& SM_i);
    
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
     * @param order_qed
     * @return the short distance contribution to the 
     * |\f$ BR(K_L \rightarrow \mu \bar{\mu}) \f$|, 
     */
    gslpp::complex BRKmumu(orders order);
    
private:
    
    StandardModel& mySM;
    CPenguinBoxMu CPB;
};

#endif	/* BR_KMUMU_H */