/* 
 * Copyright (C) 2022 HEPfit Collaboration
 *
 *
 * For the licensing terms see doc/COPYING.
 */

#ifndef EPSILONP_O_EPSILON_TH_H
#define	EPSILONP_O_EPSILON_TH_H

#include "ThObservable.h"
#include "AmpDS1.h"
#include "AmpDK2.h"

/**
 * @class EpsilonP_O_Epsilon_TH
 * @ingroup Flavour
 * @brief A class for @f$|\epsilon'_K/\epsilon_K|@f$ that parametrizes
 * direct CPV in the Kaon sector. In this implementation, both 
 * @f$|\epsilon'_K|@f$ and @f$|\epsilon_K|@f$ are computed
 * using their theoretical expression.
 * @author HEPfit Collaboration
 * @copyright GNU General Public License
 * @details This class is used to compute the theoretical value of
 * @f$|\epsilon'_K/\epsilon_K|@f$, using theory for numerator and denominator. 
 * This parameter gets contributions both
 * from the SM and many NP models.
 * 
 * 
 *
 * @anchor EpsilonP_O_Epsilon_THParameters
 * <h3>%Model parameters</h3>
 *
 * The model parameters of %EpsilonP_O_Epsilon_TH are summarized below:
 * <table class="model">
 * <tr>
 *   <th>Label</th>
 *   <th>LaTeX symbol</th>
 *   <th>Description</th>
 * </tr>
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
class EpsilonP_O_Epsilon_TH : public ThObservable, AmpDS1, AmpDK2 {
public:   
    /**
     * constructor
     * @param SM_i a reference to an object of class StandardModel, or of one of its extensions 
     */
    EpsilonP_O_Epsilon_TH(const StandardModel& SM_i);
    
    /**
     * 
     * @return theoretical value of @f$|\epsilon'_K/\epsilon_K|@f$
     */
    double computeThValue();
    
private:
    
};

#endif	/* EPSILONP_O_EPSILON_TH_H */

