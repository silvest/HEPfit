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
 *  <tr>
 *   <td class="mod_name">%PhSp_KP</td>
 *   <td class="mod_symb">@f$@f$</td>
 *   <td class="mod_desc">The phase space integral of  @f$K^+\to\pi^+ \nu\bar{\nu}@f$ .</td>
 * </tr>
 * <tr>
 *   <td class="mod_name">%DeltaP_cu</td>
 *   <td class="mod_symb">@f$@f$</td>
 *   <td class="mod_desc">The long-distance correction to the charm contribution of \f$K^+\to\pi^+\nu\bar{\nu}\f$.</td>
 * </tr>
 * <tr>
 *   <td class="mod_name">%IB_Kp</td>
 *   <td class="mod_symb">@f$@f$</td>
 *   <td class="mod_desc"> Defined as @f$ f_+^{K^+\pi^+}/f_+^{K^0\pi^+}@f$</td>
 * </tr>
 * <tr>
 *   <td class="mod_name">%Vus_fpK0Pip</td>
 *   <td class="mod_symb">@f$@f$</td>
 *   <td class="mod_desc">Defined as the norm of @f$ f_+^{K^0\pi^+} * V_{us}@f$</td>
 * </tr>
 * <tr>
 *   <td class="mod_name">%Delta_EM</td>
 *   <td class="mod_symb">@f$@f$</td>
 *   <td class="mod_desc"> Electromagnetic corrections to the decay.</td>
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
    BR_Kppnunu(const StandardModel& SM_i);
    
    /**
     * 
     * @return theoretical value of |\f$ BR(K^{+} \rightarrow \pi^+ \nu \bar{\nu}) \f$|, 
     * See arxiv: 0705.2025v2 eq. 7
     */
    double computeThValue();
    
    
    /**
     * 
     * @return The prefactor of the branching ratio. The inclusion of the factor c0 (prefactor of the Wilson Coefficients) is done in this function. See arxiv: 0705.2025v2  eq. 52
     */
    double k_plus();
   
    /**
     * 
     * @param order
     * @param order_qed
     * @return the short distance contribution to the 
     * |\f$ BR(K^{+} \rightarrow \pi^{+} \nu \bar{\nu}) \f$|, for example
     * see 0805.4119 eq 4. Warning: There is no division by lambda^10 here
     */
    double BRKppnunu(orders order, orders_qed order_qed);
    
     /**
     * 
     * @return the long distance contribution of the charm: c0 * \lambda_c * \delta P_c * \lambda^4 
      *        where co is the prefactor of the Wilson Coefficients:  4. * GF / sqrt(2.) * alphaMz / 2. / M_PI / sW2_ND
     */
    gslpp::complex LongDistance();
    
private:
    
    const StandardModel& mySM;
};

#endif	/* BR_KPPNUNU_H */
