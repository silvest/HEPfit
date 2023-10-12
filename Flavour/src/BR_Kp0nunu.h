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
#include "Charm_Kpnunu.h"
#include "StandardModelMatching.h"
#include "EpsilonK.h"


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
  * The model parameters of %BR_Kppnunu are summarized below:
 * <table class="model">
 * <tr>
 *   <th>Label</th>
 *   <th>LaTeX symbol</th>
 *   <th>Description</th>
 * </tr>
 *  <tr>
 *   <td class="mod_name">%PhSp_KL</td>
 *   <td class="mod_symb">@f$@f$</td>
 *   <td class="mod_desc">The phase space integral of the decay.</td>
 * </tr>
 * <tr>
 *   <td class="mod_name">%DeltaP_cu</td>
 *   <td class="mod_symb">@f$@f$</td>
 *   <td class="mod_desc">The long-distance correction to the charm contribution.</td>
 * </tr>
 * <tr>
 *   <td class="mod_name">%IB_Kp</td>
 *   <td class="mod_symb">@f$@f$</td>
 *   <td class="mod_desc"> Defined as @f$ f_+^{K^+\pi^+}/f_+^{K^0\pi^+}@f$</td>
 * </tr>
 * <tr>
 *   <td class="mod_name">%IB_KL</td>
 *   <td class="mod_symb">@f$@f$</td>
 *   <td class="mod_desc">Defined as @f$ (f_+^{K^+\pi^0}/f_+^{K^0\pi^+})*(f_+^{K^0\pi^0}/f_+^{K^+\pi^+}) @f$</td>
 * </tr>
 * <tr>
 *   <td class="mod_name">%IB_K0p</td>
 *   <td class="mod_symb">@f$@f$</td>
 *   <td class="mod_desc">Defined as @f$ f_+^{K^+\pi^0}/f_+^{K^0\pi^+}@f$</td>
 * </tr>
 * <tr>
 *   <td class="mod_name">%Vus_fpK0Pip</td>
 *   <td class="mod_symb">@f$@f$</td>
 *   <td class="mod_desc">Defined as the norm of @f$ f_+^{K^0\pi^+} * V_{us}@f$</td>
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
    BR_Kp0nunu(const StandardModel& SM_i);
    
    /**
     * 
     * @return theoretical value of \f$ BR(K_L \rightarrow \pi^{0} \nu \bar{\nu}) \f$, 
     * for example see hep-ph/0603079 section 2.3
     */
    double computeThValue();
    
    
    /**
     * 
     * @return The prefactor of the branching ratio. The inclusion of the factor c0 (prefactor of the Wilson Coefficients) is done in this function. See hep-ph/9607447v1 
     */
    double k_zero();
    
    /**
     * 
     * @param order
     * @param order_qed
     * @return the short distance contribution to the 
     * \f$ BR(K_{L} \rightarrow \pi^{0} \nu \bar{\nu}) \f$,
     * See hep-ph/9607447v1 for the non approximate formula
     */
    double BRKp0nunu(orders order, orders_qed order_qed);
    
    /**
     * 
     * @return the long distance contribution of the charm: c0 * \lambda_c * \delta P_c * \lambda^4 
      *        where co is the prefactor of the Wilson Coefficients:  4. * GF / sqrt(2.) * alphaMz / 2. / M_PI / sW2_ND
     */
    gslpp::complex LongDistance();
    

private:
    
    const StandardModel& mySM;
    
};

#endif	/* BR_KP0NUNU_H */
