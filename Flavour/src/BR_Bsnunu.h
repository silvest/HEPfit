/* 
 * Copyright (C) 2012 HEPfit Collaboration
 *
 *
 * For the licensing terms see doc/COPYING.
 */

#ifndef BR_BSNUNU_H
#define	BR_BSNUNU_H

#include "ThObservable.h"
#include "OrderScheme.h"
#include "gslpp.h"


/**
 * @class BR_Bsnunu
 * @ingroup Flavour
 * @brief A class for the branching ratio of \f$B_s\to \nu\nu\f$
 * @author HEPfit Collaboration
 * @copyright GNU General Public License
 * @details This class is used to compute the theoretical value of
 * the branching ratio of \f$B_s\to \nu\nu\f$.
 * 
 * 
 *
 * @anchor BR_Bsnunu
 * <h3>%Model parameters</h3>
 *
 * The model parameters of %BR_Bsnunu are summarized below:
 * <table class="model">
 * <tr>
 *   <th>Label</th>
 *   <th>LaTeX symbol</th>
 *   <th>Description</th>
 * </tr>
 * <tr>
 *   <td class="mod_name">%Br_B_Xcenu</td>
 *   <td class="mod_symb">@f$\mathrm{BR}(B\to X_ce\nu)@f$</td>
 *   <td class="mod_desc">The experimental value for the branching ratio of \f$B\to X_c e\nu\f$.</td>
 * </tr>
 * </table>
 * 
 */
class BR_Bsnunu : public ThObservable {
public:   
    /**
    * @brief Constructor.
    * @param[in] SM_i a reference to an object of type StandardModel
    */
    BR_Bsnunu(StandardModel& SM_i);
    
    /**
     * 
     * @brief hep-ph/9512380v2
     * @return theoretical value of |\f$ BR(B_s \rightarrow \nu \bar{\nu}) \f$|
     */
    double computeThValue();
    
    
protected:
    
    /**
    * @brief the short distance contribution to the 
    * |\f$ BR(B_d \rightarrow \nu \bar{\nu}) \f$|
    * @param[in] order the %QCD order of the computation
    */
    gslpp::complex BRBsnunu(orders order);
    
private:
    
    StandardModel& mySM;
};

#endif	/* BR_BSNUNU_H */