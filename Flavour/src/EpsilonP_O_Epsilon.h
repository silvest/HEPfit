/* 
 * Copyright (C) 2012 HEPfit Collaboration
 *
 *
 * For the licensing terms see doc/COPYING.
 */

#ifndef EPSILONP_O_EPSILON_H
#define	EPSILONP_O_EPSILON_H

#include "ThObservable.h"
#include "AmpDS1.h"

/**
 * @class EpsilonP_O_Epsilon
 * @ingroup Flavour
 * @brief A class for @f$|\epsilon'_K/\epsilon_K|@f$ that parametrizes
 * direct CPV in the Kaon sector
 * @author HEPfit Collaboration
 * @copyright GNU General Public License
 * @details This class is used to compute the theoretical value of
 * @f$|\epsilon'_K/\epsilon_K|@f$. This parameter gets contributions both
 * from the SM and many NP models.
 * 
 * 
 *
 * @anchor EpsilonP_O_EpsilonParameters
 * <h3>%Model parameters</h3>
 *
 * The model parameters of %EpsilonP_O_Epsilon are summarized below:
 * <table class="model">
 * <tr>
 *   <td class="mod_name">%ReA0_Kd</td>
 *   <td class="mod_symb">@f${\cal Re}(A_0(K\to\pi\pi))@f$</td>
 *   <td class="mod_desc">The experimental value of the real part of the amplitude for \f$K^0\to\pi\pi\f$ with \f$\Delta I=0\f$.</td>
 * </tr>
 * <tr>
 *   <td class="mod_name">%ReA2_Kd</td>
 *   <td class="mod_symb">@f${\cal Re}(A_2(K\to\pi\pi))@f$</td>
 *   <td class="mod_desc">the experimental value of the real part of the amplitude for \f$K^0\to\pi\pi\f$ with \f$\Delta I=2\f$.</td>
 * </tr>
 * <tr>
 *   <td class="mod_name">%Omega_eta_etap</td>
 *   <td class="mod_symb">@f$\Omega_{\eta/\eta'}@f$</td>
 *   <td class="mod_desc">The isospin breaking contribution in \f$K^0\to\pi\pi\f$.</td>
 * </tr>
 * </table>
 * 
 */
class EpsilonP_O_Epsilon : public ThObservable, AmpDS1 {
public:   
    /**
     * constructor
     * @param Flavour
     */
    EpsilonP_O_Epsilon(const StandardModel& SM_i);
    
    /**
     * 
     * @return theoretical value of @f$|\epsilon'_K/\epsilon_K|@f$
     */
    double computeThValue();
    
private:
    
};

#endif	/* EPSILONP_O_EPSILON_H */
