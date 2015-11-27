/*
 * Copyright (C) 2014 HEPfit Collaboration
 * All rights reserved.
 *
 * For the licensing terms see doc/COPYING.
 */

#ifndef NPEFFECTIVEGIMR_H
#define	NPEFFECTIVEGIMR_H

#include "gslpp.h"
#include "NPbase.h"
#include <string.h>
#include <stdexcept>


/**
 * @class NPEffectiveGIMR
 * @brief A model class for new physics in the form of the
 * dimension-six effective Lagrangian.
 * @ingroup NewPhysics
 * @author HEPfit Collaboration
 * @copyright GNU General Public License
 * @details This is a Model class containing parameters and functions
 * associated with the general dimension-six effective Lagrangian.
 * (Use the model name "NPEffectiveGIMR_LFU_QFU" to asumme lepton and quark flavour universality)
 *
 * In this class we consider the dimension-six effective Lagrangian
 * @f[
 * \mathcal{L}_\mathrm{eff}
 * = \mathcal{L}_\mathrm{SM}
 *   + \sum_i \frac{C_i}{\Lambda^2} \mathcal{O}_i
 * @f]
 * as written in the basis of \cite Grzadkowski:2010es.
 *
 * @anchor NPEffectiveGIMRInitialization
 * <h3>Initialization</h3>
 *
 * After creating an instance of the current class with the constructor
 * NPEffectiveGIMR(), it is required to call the initialization method
 * InitializeModel().
 * In the Monte Carlo run, the constructor as well as the initialization
 * method are called in InputParser::ReadParameters().
 *
 *
 * @anchor NPEffectiveGIMRParameters
 * <h3>%Model parameters</h3>
 *
 * The model parameters of %NPEffectiveGIMR are summarized below:
 * <table class="model">
 * <tr>
 *   <th>Label</th>
 *   <th>LaTeX symbol</th>
 *   <th>Description</th>
 * </tr>
 * <tr>
 *   <td class="mod_name">%CW </td>
 *   <td class="mod_symb">\f$C_{W} \f$</td>
 *   <td class="mod_desc">The coefficient of the operator
 *   \f${\cal O}_{W}=\varepsilon_{abc}W_{\mu}^{a\nu} W_{\nu}^{b\rho}W_{\rho}^{b\mu}\f$. </td>
 * </tr>
 * <tr>
 *   <td class="mod_name">%CHG </td>
 *   <td class="mod_symb">\f$C_{HG} \f$</td>
 *   <td class="mod_desc">The coefficient of the operator
 *   \f${\cal O}_{HG}=\big(H^\dagger H\big)G_{\mu\nu}^A G^{A\mu\nu}\f$. </td>
 * </tr>
 * <tr>
 *   <td class="mod_name">%CHW </td>
 *   <td class="mod_symb">\f$C_{HW} \f$</td>
 *   <td class="mod_desc">The coefficient of the operator
 *   \f${\cal O}_{HW}=\big(H^\dagger H\big)W_{\mu\nu}^a W^{a\mu\nu}\f$. </td>
 * </tr>
 * <tr>
 *   <td class="mod_name">%CHB </td>
 *   <td class="mod_symb">\f$C_{HB} \f$</td>
 *   <td class="mod_desc">The coefficient of the operator
 *   \f${\cal O}_{HB}=\big(H^\dagger H\big)B_{\mu\nu} B^{\mu\nu}\f$. </td>
 * </tr>
 * <tr>
 *   <td class="mod_name">%CWB </td>
 *   <td class="mod_symb">\f$C_{WB} \f$</td>
 *   <td class="mod_desc">The coefficient of the operator
 *   \f${\cal O}_{HWB}=\big(H^\dagger\tau^a H\big)W_{\mu\nu}^a B^{\mu\nu}\f$. </td>
 * </tr>
 * <tr>
 *   <td class="mod_name">%CHD </td>
 *   <td class="mod_symb">\f$C_{HD}\f$</td>
 *   <td class="mod_desc">The coefficient of the operator
 *   \f${\cal O}_{HD}=\big|H^\dagger D_\mu H\big|^2\f$. </td>
 * </tr>
 * <tr>
 *   <td class="mod_name">%CHbox </td>
 *   <td class="mod_symb">\f$C_{H\Box}\f$</td>
 *   <td class="mod_desc">The coefficient of the operator
 *   \f${\cal O}_{H\Box}=\big(H^\dagger H\big)\Box\big(H^\dagger H\big)\f$. </td>
 * </tr>
 * <tr>
 *   <td class="mod_name">%CH </td>
 *   <td class="mod_symb">\f$C_{H}\f$</td>
 *   <td class="mod_desc">The coefficient of the operator
 *   \f${\cal O}_{H}=\big(H^\dagger H\big)^3\f$. </td>
 * </tr>
 * <tr>
 *   <td class="mod_name">%CHL1_kk, CHL1_klr, CHL1_kli </td>
 *   <td class="mod_symb">\f$ (C_{HL}^{(1)})_{kk}, \mbox{Re}\big[(C_{HL}^{(1)})_{kl}\big], \mbox{Im}\big[(C_{HL}^{(1)})_{kl}\big] \f$</td>
 *   <td class="mod_desc">The real and imaginary parts of the coefficient of the operator
 *   \f$({\cal O}_{HL}^{(1)})_{ij} =i\big(H^\dagger D_\mu H\big)
 *   \big(\overline{L^i}\,\gamma^\mu L^j\big)\f$, for \f$i,j=1,2,3\f$. </td>
 * </tr>
 * <tr>
 *   <td class="mod_name">%CHL3_kk, CHL3_klr, CHL3_kli </td>
 *   <td class="mod_symb">\f$ (C_{HL}^{(3)})_{kk}, \mbox{Re}\big[(C_{HL}^{(3)})_{kl}\big], \mbox{Im}\big[(C_{HL}^{(3)})_{kl}\big] \f$</td>
 *   <td class="mod_desc">The real and imaginary parts of the coefficient of the operator
 *   \f$({\cal O}_{HL}^{(3)})_{ij} =i\big(H^\dagger \tau^a D_\mu H\big)
 *   \big(\overline{L^i}\,\gamma^\mu \tau^a L^j\big)\f$, for \f$i,j=1,2,3\f$. </td>
 * </tr>
 * <tr>
 *   <td class="mod_name">%CHe_kk, CHe_klr, CHe_kli </td>
 *   <td class="mod_symb">\f$ (C_{HE})_{kk}, \mbox{Re}\big[(C_{HE})_{kl}\big], \mbox{Im}\big[(C_{HE})_{kl}\big] \f$</td>
 *   <td class="mod_desc">The real and imaginary parts of the coefficient of the operator
 *   \f$({\cal O}_{HE})_{ij} =i\big(H^\dagger D_\mu H\big)
 *   \big(\overline{E^i}\,\gamma^\mu E^j\big)\f$, for \f$i,j=1,2,3\f$. </td>
 * </tr>
 * <tr>
 *   <td class="mod_name">%CHQ1_kk, CHQ1_klr, CHQ1_kli </td>
 *   <td class="mod_symb">\f$ (C_{HQ}^{(1)})_{kk}, \mbox{Re}\big[(C_{HQ}^{(1)})_{kl}\big], \mbox{Im}\big[(C_{HQ}^{(1)})_{kl}\big] \f$</td>
 *   <td class="mod_desc">The real and imaginary parts of the coefficient of the operator
 *   \f$({\cal O}_{HQ}^{(1)})_{ij} =i\big(H^\dagger D_\mu H\big)
 *   \big(\overline{Q^i}\,\gamma^\mu Q^j\big)\f$, for \f$i,j=1,2,3\f$. </td>
 * </tr>
 * <tr>
 *   <td class="mod_name">%CHQ3_kk, CHQ3_klr, CHQ3_kli </td>
 *   <td class="mod_symb">\f$ (C_{HQ}^{(3)})_{kk}, \mbox{Re}\big[(C_{HQ}^{(3)})_{kl}\big], \mbox{Im}\big[(C_{HQ}^{(3)})_{kl}\big] \f$</td>
 *   <td class="mod_desc">The real and imaginary parts of the coefficient of the operator
 *   \f$({\cal O}_{HQ}^{(3)})_{ij} =i\big(H^\dagger \tau^a D_\mu H\big)
 *   \big(\overline{Q^i}\,\gamma^\mu \tau^a Q^j\big)\f$, for \f$i,j=1,2,3\f$. </td>
 * </tr>
 * <tr>
 *   <td class="mod_name">%CHu_kk, CHu_klr, CHu_kli </td>
 *   <td class="mod_symb">\f$ (C_{HU})_{kk}, \mbox{Re}\big[(C_{HU})_{kl}\big], \mbox{Im}\big[(C_{HU})_{kl}\big] \f$</td>
 *   <td class="mod_desc">The real and imaginary parts of the coefficient of the operator
 *   \f$({\cal O}_{HU})_{ij} =i\big(H^\dagger D_\mu H\big)
 *   \big(\overline{U^i}\,\gamma^\mu U^j\big)\f$, for \f$i,j=1,2,3\f$. </td>
 * </tr>
 * <tr>
 *   <td class="mod_name">%CHd_kk, CHd_klr, CHd_kli </td>
 *   <td class="mod_symb">\f$ (C_{HD})_{kk}, \mbox{Re}\big[(C_{HD})_{kl}\big], \mbox{Im}\big[(C_{HD})_{kl}\big] \f$</td>
 *   <td class="mod_desc">The real and imaginary parts of the coefficient of the operator
 *   \f$({\cal O}_{HD})_{ij} =i\big(H^\dagger D_\mu H\big)
 *   \big(\overline{D^i}\,\gamma^\mu D^j\big)\f$, for \f$i,j=1,2,3\f$. </td>
 * </tr>
  * <tr>
 *   <td class="mod_name">%CHud_klr, CHud_kli </td>
 *   <td class="mod_symb">\f$\mbox{Re}\big[(C_{HUD})_{kl}\big], \mbox{Im}\big[(C_{HUD})_{kl}\big] \f$</td>
 *   <td class="mod_desc">The real and imaginary parts of the coefficient of the operator
 *   \f$({\cal O}_{HUD})_{ij} =i\big(\widetilde{H}^\dagger D_\mu H\big)
 *   \big(\overline{U^i}\,\gamma^\mu D^j\big)\f$, for \f$i,j=1,2,3\f$. </td>
 * </tr>
 * <tr>
 *   <td class="mod_name">%CeH_klr, CeH_kli </td>
 *   <td class="mod_symb">\f$\mbox{Re}\big[(C_{EH})_{kl}\big], \mbox{Im}\big[(C_{EH})_{kl}\big] \f$</td>
 *   <td class="mod_desc">The real and imaginary parts of the coefficient of the operator
 *   \f$({\cal O}_{EH})_{ij} =\big(H^\dagger H\big)
 *   \big(\overline{L^i}\,H E^j\big)\f$, for \f$i,j=1,2,3\f$. </td>
 * </tr>
 * <tr>
 *   <td class="mod_name">%CuH_klr, CuH_kli </td>
 *   <td class="mod_symb">\f$\mbox{Re}\big[(C_{UH})_{kl}\big], \mbox{Im}\big[(C_{UH})_{kl}\big] \f$</td>
 *   <td class="mod_desc">The real and imaginary parts of the coefficient of the operator
 *   \f$({\cal O}_{UH})_{ij} =\big(H^\dagger H\big)
 *   \big(\overline{Q^i}\,\widetilde{H} U^j\big)\f$, for \f$i,j=1,2,3\f$. </td>
 * </tr>
 * <tr>
 *   <td class="mod_name">%CdH_klr, CdH_kli </td>
 *   <td class="mod_symb">\f$\mbox{Re}\big[(C_{DH})_{kl}\big], \mbox{Im}\big[(C_{DH})_{kl}\big] \f$</td>
 *   <td class="mod_desc">The real and imaginary parts of the coefficient of the operator
 *   \f$({\cal O}_{DH})_{ij} =\big(H^\dagger H\big)
 *   \big(\overline{Q^i}\,H D^j\big)\f$, for \f$i,j=1,2,3\f$. </td>
 * </tr>
 * <tr>
 *   <td class="mod_name">%CLL_1221, CLL_2112 </td>
 *   <td class="mod_symb">\f$(C_{LL})_{1221,2112}\f$</td>
 *   <td class="mod_desc">The coefficient of the operator
 *   \f$({\cal O}_{LL})_{ijkl}=\big(\overline{L^i}\,\gamma^\mu L^j\big)
 *   \big(\overline{L^k}\,\gamma_\mu L^l\big)\f$, for \f$ijkl=1221,2112\f$. </td>
 * </tr>
 * <tr>
 *   <td class="mod_name">%CLQ1 </td>
 *   <td class="mod_symb">\f$C_{LQ}^{(1)}\f$</td>
 *   <td class="mod_desc">The coefficient of the operator
 *   \f$({\cal O}_{LQ}^{(1)})_{ijkl}=\big(\overline{L^i}\,\gamma^\mu L^j\big)
 *   \big(\overline{Q^k}\,\gamma_\mu Q^l\big)\f$. </td>
 * </tr>
 * <tr>
 *   <td class="mod_name">%CLQ3 </td>
 *   <td class="mod_symb">\f$C_{LQ}^{(3)}\f$</td>
 *   <td class="mod_desc">The coefficient of the operator
 *   \f$({\cal O}_{LQ}^{(3)})_{ijkl}=\big(\overline{L^i}\,\gamma^\mu \tau_a L^j\big)
 *   \big(\overline{Q^k}\,\gamma_\mu \tau_a Q^l\big)\f$. </td>
 * </tr>
 * <tr>
 *   <td class="mod_name">%Cee </td>
 *   <td class="mod_symb">\f$C_{EE}\f$</td>
 *   <td class="mod_desc">The coefficient of the operator
 *   \f$({\cal O}_{EE})_{ijkl}=\big(\overline{E^i}\,\gamma^\mu E^j\big)
 *   \big(\overline{E^k}\,\gamma_\mu E^l\big)\f$. </td>
 * </tr>
 * <tr>
 *   <td class="mod_name">%Ceu </td>
 *   <td class="mod_symb">\f$C_{EU}\f$</td>
 *   <td class="mod_desc">The coefficient of the operator
 *   \f$({\cal O}_{EU})_{ijkl}=\big(\overline{E^i}\,\gamma^\mu E^j\big)
 *   \big(\overline{U^k}\,\gamma_\mu U^l\big)\f$. </td>
 * </tr>
 * <tr>
 *   <td class="mod_name">%Ced </td>
 *   <td class="mod_symb">\f$C_{ED}\f$</td>
 *   <td class="mod_desc">The coefficient of the operator
 *   \f$({\cal O}_{ED})_{ijkl}=\big(\overline{E^i}\,\gamma^\mu E^j\big)
 *   \big(\overline{D^k}\,\gamma_\mu D^l\big)\f$. </td>
 * </tr>
 * <tr>
 *   <td class="mod_name">%CLe </td>
 *   <td class="mod_symb">\f$C_{LE}\f$</td>
 *   <td class="mod_desc">The coefficient of the operator
 *   \f$({\cal O}_{LE})_{ijkl}=\big(\overline{L^i}\,\gamma^\mu L^j\big)
 *   \big(\overline{E^k}\,\gamma_\mu E^l\big)\f$. </td>
 * </tr>
 * <tr>
 *   <td class="mod_name">%CLu </td>
 *   <td class="mod_symb">\f$C_{LU}\f$</td>
 *   <td class="mod_desc">The coefficient of the operator
 *   \f$({\cal O}_{LU})_{ijkl}=\big(\overline{L^i}\,\gamma^\mu L^j\big)
 *   \big(\overline{U^k}\,\gamma_\mu U^l\big)\f$. </td>
 * </tr>
 * <tr>
 *   <td class="mod_name">%CLd </td>
 *   <td class="mod_symb">\f$C_{LD}\f$</td>
 *   <td class="mod_desc">The coefficient of the operator
 *   \f$({\cal O}_{LD})_{ijkl}=\big(\overline{L^i}\,\gamma^\mu L^j\big)
 *   \big(\overline{D^k}\,\gamma_\mu D^l\big)\f$. </td>
 * </tr>
 * <tr>
 *   <td class="mod_name">%CQe </td>
 *   <td class="mod_symb">\f$C_{QE}\f$</td>
 *   <td class="mod_desc">The coefficient of the operator
 *   \f$({\cal O}_{QE})_{ijkl}=\big(\overline{Q^i}\,\gamma^\mu Q^j\big)
 *   \big(\overline{E^k}\,\gamma_\mu E^l\big)\f$. </td>
 * </tr>
 * <tr>
 *   <td class="mod_name">%Lambda_NP </td>
 *   <td class="mod_symb">\f$\Lambda \f$</td>
 *   <td class="mod_desc">The new physics scale. </td>
 * </tr>
 * <tr>
 *   <td class="mod_name">%eVBFE_i </td>
 *   <td class="mod_symb">\f$\varepsilon_{VBF}^i(E)\f$</td>
 *   <td class="mod_desc">The theoretical uncertainty in the coefficient multiplying
 *  the effective coupling \f$g_i\f$ in the VBF production cross section at 
 *  Tevatron (\f$E=2\f$) or the LHC (\f$E=78\f$).
 * \f$(g_i=g_{HZZ}^{(1,2,3)}, g_{HZA}^{(1,2)}, g_{HAA}, g_{HWW}^{(1,2,3)}, g_{Hgg}, g_{HZuu,HZdd}^{L,R}, 
 * g_{HWud}^{L}, g_{Zuu,Zdd}^{L,R}, g_{Wud}^{L})\f$ </td>
 * </tr>
 * <tr>
 *   <td class="mod_name">%eWHE_i </td>
 *   <td class="mod_symb">\f$\varepsilon_{WH}^i(E)\f$</td>
 *   <td class="mod_desc">The theoretical uncertainty in the coefficient multiplying
 *  the effective coupling \f$g_i\f$ in the WH production cross section at 
 *  Tevatron (\f$E=2\f$) or the LHC (\f$E=78\f$). 
 * \f$(g_i= g_{HWW}^{(1,2,3)}, g_{HWud}^{L}, g_{Wud}^{L})\f$ </td>
 * </tr>
 * <tr>
 *   <td class="mod_name">%eZHE_i </td>
 *   <td class="mod_symb">\f$\varepsilon_{ZH}^i(E)\f$</td>
 *   <td class="mod_desc">The theoretical uncertainty in the coefficient multiplying
 *  the effective coupling \f$g_i\f$ in the ZH production cross section at 
 *  Tevatron (\f$E=2\f$) or the LHC (\f$E=78\f$). 
 * \f$(g_i=g_{HZZ}^{(1,2,3)}, g_{HZA}^{(1,2)}, g_{HZuu,HZdd}^{L,R}, 
 * g_{Zuu,Zdd}^{L,R})\f$ </td>
 * </tr>
 * <tr>
 *   <td class="mod_name">%ettHE_i </td>
 *   <td class="mod_symb">\f$\varepsilon_{ttH}^i(E)\f$</td>
 *   <td class="mod_desc">The theoretical uncertainty in the coefficient multiplying
 *  the effective coupling \f$g_i\f$ in the ttH production cross section at 
 *  Tevatron (\f$E=2\f$) or the LHC (\f$E=78\f$). 
 * \f$(g_i= g_{Htt}, g_{Hgg})\f$ </td>
 * </tr>
 * </table>
 *
 * Alternatively, when using the model name "NPEffectiveGIMR_LFU_QFU", where lepton
 * and quark flavour universality are assumed, the parameters to be used as inputs
 * for the dimension six coefficients are the following:
 * <table class="model">
 * <tr>
 *   <th>Label</th>
 *   <th>LaTeX symbol</th>
 *   <th>Description</th>
 * </tr>
 * <tr>
 *   <td class="mod_name">%CW </td>
 *   <td class="mod_symb">\f$C_{W} \f$</td>
 *   <td class="mod_desc">The coefficient of the operator
 *   \f${\cal O}_{W}=\varepsilon_{abc}W_{\mu}^{a\nu} W_{\nu}^{b\rho}W_{\rho}^{b\mu}\f$. </td>
 * </tr>
 * <tr>
 *   <td class="mod_name">%CHG </td>
 *   <td class="mod_symb">\f$C_{HG} \f$</td>
 *   <td class="mod_desc">The coefficient of the operator
 *   \f${\cal O}_{HG}=\big(H^\dagger H\big)G_{\mu\nu}^A G^{A\mu\nu}\f$. </td>
 * </tr>
 * <tr>
 *   <td class="mod_name">%CHW </td>
 *   <td class="mod_symb">\f$C_{HW} \f$</td>
 *   <td class="mod_desc">The coefficient of the operator
 *   \f${\cal O}_{HW}=\big(H^\dagger H\big)W_{\mu\nu}^a W^{a\mu\nu}\f$. </td>
 * </tr>
 * <tr>
 *   <td class="mod_name">%CHB </td>
 *   <td class="mod_symb">\f$C_{HB} \f$</td>
 *   <td class="mod_desc">The coefficient of the operator
 *   \f${\cal O}_{HB}=\big(H^\dagger H\big)B_{\mu\nu} B^{\mu\nu}\f$. </td>
 * </tr>
 * <tr>
 *   <td class="mod_name">%CWB </td>
 *   <td class="mod_symb">\f$C_{WB} \f$</td>
 *   <td class="mod_desc">The coefficient of the operator
 *   \f${\cal O}_{HWB}=\big(H^\dagger\tau^a H\big)W_{\mu\nu}^a B^{\mu\nu}\f$. </td>
 * </tr>
 * <tr>
 *   <td class="mod_name">%CHD </td>
 *   <td class="mod_symb">\f$C_{HD}\f$</td>
 *   <td class="mod_desc">The coefficient of the operator
 *   \f${\cal O}_{HD}=\big|H^\dagger D_\mu H\big|^2\f$. </td>
 * </tr>
 * <tr>
 *   <td class="mod_name">%CHbox </td>
 *   <td class="mod_symb">\f$C_{H\Box}\f$</td>
 *   <td class="mod_desc">The coefficient of the operator
 *   \f${\cal O}_{H\Box}=\big(H^\dagger H\big)\Box\big(H^\dagger H\big)\f$. </td>
 * </tr>
 * <tr>
 *   <td class="mod_name">%CH </td>
 *   <td class="mod_symb">\f$C_{H}\f$</td>
 *   <td class="mod_desc">The coefficient of the operator
 *   \f${\cal O}_{H}=\big(H^\dagger H\big)^3\f$. </td>
 * </tr>
 * <tr>
 *   <td class="mod_name">%CHL1</td>
 *   <td class="mod_symb">\f$ (C_{HL}^{(1)})_{ii} \f$</td>
 *   <td class="mod_desc">The coefficient of the operator
 *   \f$({\cal O}_{HL}^{(1)})_{ii} =i\big(H^\dagger D_\mu H\big)
 *   \big(\overline{L^i}\,\gamma^\mu L^i\big)\f$ (flavor universal). </td>
 * </tr>
 * <tr>
 *   <td class="mod_name">%CHL3 </td>
 *   <td class="mod_symb">\f$ (C_{HL}^{(3)})_{ii} \f$</td>
 *   <td class="mod_desc">The coefficient of the operator
 *   \f$({\cal O}_{HL}^{(3)})_{ii} =i\big(H^\dagger \tau^a D_\mu H\big)
 *   \big(\overline{L^i}\,\gamma^\mu \tau^a L^i\big)\f$ (flavor universal). </td>
 * </tr>
 * <tr>
 *   <td class="mod_name">%CHe </td>
 *   <td class="mod_symb">\f$ (C_{HE})_{ii} \f$</td>
 *   <td class="mod_desc">The coefficient of the operator
 *   \f$({\cal O}_{HE})_{ij} =i\big(H^\dagger D_\mu H\big)
 *   \big(\overline{E^i}\,\gamma^\mu E^i\big)\f$ (flavor universal). </td>
 * </tr>
 * <tr>
 *   <td class="mod_name">%CHQ1 </td>
 *   <td class="mod_symb">\f$ (C_{HQ}^{(1)})_{ii} \f$</td>
 *   <td class="mod_desc">The coefficient of the operator
 *   \f$({\cal O}_{HQ}^{(1)})_{ii} =i\big(H^\dagger D_\mu H\big)
 *   \big(\overline{Q^i}\,\gamma^\mu Q^i\big)\f$ (flavor universal). </td>
 * </tr>
 * <tr>
 *   <td class="mod_name">%CHQ3 </td>
 *   <td class="mod_symb">\f$ (C_{HQ}^{(3)})_{ii}\f$</td>
 *   <td class="mod_desc">The coefficient of the operator
 *   \f$({\cal O}_{HQ}^{(3)})_{ii} =i\big(H^\dagger \tau^a D_\mu H\big)
 *   \big(\overline{Q^i}\,\gamma^\mu \tau^a Q^i\big)\f$ (flavor universal). </td>
 * </tr>
 * <tr>
 *   <td class="mod_name">%CHu </td>
 *   <td class="mod_symb">\f$ (C_{HU})_{ii} \f$</td>
 *   <td class="mod_desc">The coefficient of the operator
 *   \f$({\cal O}_{HU})_{ii} =i\big(H^\dagger D_\mu H\big)
 *   \big(\overline{U^i}\,\gamma^\mu U^i\big)\f$ (flavor universal). </td>
 * </tr>
 * <tr>
 *   <td class="mod_name">%CHd </td>
 *   <td class="mod_symb">\f$ (C_{HD})_{ii} \f$</td>
 *   <td class="mod_desc">The coefficient of the operator
 *   \f$({\cal O}_{HD})_{ii} =i\big(H^\dagger D_\mu H\big)
 *   \big(\overline{D^i}\,\gamma^\mu D^i\big)\f$ (flavor universal). </td>
 * </tr>
  * <tr>
 *   <td class="mod_name">%CHud_r, CHud_i </td>
 *   <td class="mod_symb">\f$\mbox{Re}\big[(C_{HUD})_{ii}\big], \mbox{Im}\big[(C_{HUD})_{ii}\big] \f$</td>
 *   <td class="mod_desc">The real and imaginary parts of the coefficient of the operator
 *   \f$({\cal O}_{HUD})_{ii} =i\big(\widetilde{H}^\dagger D_\mu H\big)
 *   \big(\overline{U^i}\,\gamma^\mu D^i\big)\f$ (flavor universal). </td>
 * </tr>
 * <tr>
 *   <td class="mod_name">%CeH_r, CeH_i </td>
 *   <td class="mod_symb">\f$\mbox{Re}\big[(C_{EH})_{ii}\big], \mbox{Im}\big[(C_{EH})_{ii}\big] \f$</td>
 *   <td class="mod_desc">The real and imaginary parts of the coefficient of the operator
 *   \f$({\cal O}_{EH})_{ii} =\big(H^\dagger H\big)
 *   \big(\overline{L^i}\,H E^i\big)\f$ (flavor universal). </td>
 * </tr>
 * <tr>
 *   <td class="mod_name">%CuH_r, CuH_i </td>
 *   <td class="mod_symb">\f$\mbox{Re}\big[(C_{UH})_{ii}\big], \mbox{Im}\big[(C_{UH})_{ii}\big] \f$</td>
 *   <td class="mod_desc">The real and imaginary parts of the coefficient of the operator
 *   \f$({\cal O}_{UH})_{ii} =\big(H^\dagger H\big)
 *   \big(\overline{Q^i}\,\widetilde{H} U^i\big)\f$ (flavor universal). </td>
 * </tr>
 * <tr>
 *   <td class="mod_name">%CdH_r, CdH_i </td>
 *   <td class="mod_symb">\f$\mbox{Re}\big[(C_{DH})_{ii}\big], \mbox{Im}\big[(C_{DH})_{ii}\big] \f$</td>
 *   <td class="mod_desc">The real and imaginary parts of the coefficient of the operator
 *   \f$({\cal O}_{DH})_{ii} =\big(H^\dagger H\big)
 *   \big(\overline{Q^i}\,H D^i\big)\f$ (flavor universal). </td>
 * </tr>
 * <tr>
 *   <td class="mod_name">%CLL </td>
 *   <td class="mod_symb">\f$(C_{LL})_{1221,2112}\f$</td>
 *   <td class="mod_desc">The coefficient of the operator
 *   \f$({\cal O}_{LL})_{ijkl}=\big(\overline{L^i}\,\gamma^\mu L^j\big)
 *   \big(\overline{L^k}\,\gamma_\mu L^l\big)\f$, for \f$ijkl=1221,2112\f$. </td>
 * </tr>
 * <tr>
 *   <td class="mod_name">%CLQ1 </td>
 *   <td class="mod_symb">\f$C_{LQ}^{(1)}\f$</td>
 *   <td class="mod_desc">The coefficient of the operator
 *   \f$({\cal O}_{LQ}^{(1)})_{ijkl}=\big(\overline{L^i}\,\gamma^\mu L^j\big)
 *   \big(\overline{Q^k}\,\gamma_\mu Q^l\big)\f$. </td>
 * </tr>
 * <tr>
 *   <td class="mod_name">%CLQ3 </td>
 *   <td class="mod_symb">\f$C_{LQ}^{(3)}\f$</td>
 *   <td class="mod_desc">The coefficient of the operator
 *   \f$({\cal O}_{LQ}^{(3)})_{ijkl}=\big(\overline{L^i}\,\gamma^\mu \tau_a L^j\big)
 *   \big(\overline{Q^k}\,\gamma_\mu \tau_a Q^l\big)\f$. </td>
 * </tr>
 * <tr>
 *   <td class="mod_name">%Cee </td>
 *   <td class="mod_symb">\f$C_{EE}\f$</td>
 *   <td class="mod_desc">The coefficient of the operator
 *   \f$({\cal O}_{EE})_{ijkl}=\big(\overline{E^i}\,\gamma^\mu E^j\big)
 *   \big(\overline{E^k}\,\gamma_\mu E^l\big)\f$. </td>
 * </tr>
 * <tr>
 *   <td class="mod_name">%Ceu </td>
 *   <td class="mod_symb">\f$C_{EU}\f$</td>
 *   <td class="mod_desc">The coefficient of the operator
 *   \f$({\cal O}_{EU})_{ijkl}=\big(\overline{E^i}\,\gamma^\mu E^j\big)
 *   \big(\overline{U^k}\,\gamma_\mu U^l\big)\f$. </td>
 * </tr>
 * <tr>
 *   <td class="mod_name">%Ced </td>
 *   <td class="mod_symb">\f$C_{ED}\f$</td>
 *   <td class="mod_desc">The coefficient of the operator
 *   \f$({\cal O}_{ED})_{ijkl}=\big(\overline{E^i}\,\gamma^\mu E^j\big)
 *   \big(\overline{D^k}\,\gamma_\mu D^l\big)\f$. </td>
 * </tr>
 * <tr>
 *   <td class="mod_name">%CLe </td>
 *   <td class="mod_symb">\f$C_{LE}\f$</td>
 *   <td class="mod_desc">The coefficient of the operator
 *   \f$({\cal O}_{LE})_{ijkl}=\big(\overline{L^i}\,\gamma^\mu L^j\big)
 *   \big(\overline{E^k}\,\gamma_\mu E^l\big)\f$. </td>
 * </tr>
 * <tr>
 *   <td class="mod_name">%CLu </td>
 *   <td class="mod_symb">\f$C_{LU}\f$</td>
 *   <td class="mod_desc">The coefficient of the operator
 *   \f$({\cal O}_{LU})_{ijkl}=\big(\overline{L^i}\,\gamma^\mu L^j\big)
 *   \big(\overline{U^k}\,\gamma_\mu U^l\big)\f$. </td>
 * </tr>
 * <tr>
 *   <td class="mod_name">%CLd </td>
 *   <td class="mod_symb">\f$C_{LD}\f$</td>
 *   <td class="mod_desc">The coefficient of the operator
 *   \f$({\cal O}_{LD})_{ijkl}=\big(\overline{L^i}\,\gamma^\mu L^j\big)
 *   \big(\overline{D^k}\,\gamma_\mu D^l\big)\f$. </td>
 * </tr>
 * <tr>
 *   <td class="mod_name">%CQe </td>
 *   <td class="mod_symb">\f$C_{QE}\f$</td>
 *   <td class="mod_desc">The coefficient of the operator
 *   \f$({\cal O}_{QE})_{ijkl}=\big(\overline{Q^i}\,\gamma^\mu Q^j\big)
 *   \big(\overline{E^k}\,\gamma_\mu E^l\big)\f$. </td>
 * </tr>
 * <tr>
 *   <td class="mod_name">%Lambda_NP </td>
 *   <td class="mod_symb">\f$\Lambda \f$</td>
 *   <td class="mod_desc">The new physics scale. </td>
 * </tr>
 * </table>
 * (The parameters associated to the theoretical uncertainties, \f$\varepsilon_{X}^i(E)\f$,
 * are the same for both "NPEffectiveGIMR" and "NPEffectiveGIMR_LFU_QFU".)
 * 
 * Finally, if the flag %MwInput (see below) is set to TRUE, one must also specify the input value for
 * the W mass via the following parameter (Warning: The W width is not implemented in this case):
 * 
 * <table class="model">
 * <tr>
 *   <th>Label</th>
 *   <th>LaTeX symbol</th>
 *   <th>Description</th>
 * </tr>
 * <tr>
 *   <td class="mod_name">%MwInput </td>
 *   <td class="mod_symb">\f$M_{W} \f$</td>
 *   <td class="mod_desc">The input value for the W mass in GeV.</td>
 * </tr>
 * </table>
 *
 * @anchor NPEffectiveGIMRFlags
 * <h3>%Model flags</h3>
 *
 * The Flags of NPEffectiveGIMR are summarized below:
 * <table class="model">
 * <tr>
 *   <th>Label</th>
 *   <th>Value</th>
 *   <th>Description</th>
 * </tr>
 * <tr>
 *   <td class="mod_name">%MwInput</td>
 *   <td class="mod_valu">TRUE&nbsp;/&nbsp;<b>FALSE</b></td>
 *   <td class="mod_desc">This flag is set to TRUE if the W mass is taken as an input parameter.
 *   (Warning: The W width is not implemented in this case.) The default value is FALSE. </td>
 * </tr>
 * <tr>
 *   <td class="mod_name">%QuadraticTerms</td>
 *   <td class="mod_valu">TRUE&nbsp;/&nbsp;<b>FALSE</b></td>
 *   <td class="mod_desc">This flag is set to TRUE if the quadratic terms in Higgs cross sections and widths are switched on.
 *   The default value is FALSE; new physics contributions are linearized.</td>
 * </tr>
 * </table>
 *
 *
 * @anchor NPEffectiveGIMRFunctions
 * <h3>Important member functions</h3>
 *
 * See the base classes of the current class.
 *
 */
class NPEffectiveGIMR : public NPbase {
public:

    /**
     *　@brief The number of the model parameters in %NPEffectiveGIMR. 
     */
    static const int NNPEffectiveGIMRVars = 210;

    /**
     * @brief A string array containing the labels of the model parameters in
     * %NPEffectiveGIMR.
     */
    static const std::string NPEffectiveGIMRVars[NNPEffectiveGIMRVars];

    /**
     *　@brief The number of the model parameters in %NPEffectiveGIMR
     * with lepton and quark flavour universalities.
     */
    static const int NNPEffectiveGIMRVars_LFU_QFU = 114;

    /**
     * @brief A string array containing the labels of the model parameters in
     * %NPEffectiveGIMR with lepton and quark flavour universalities.
     */
    static const std::string NPEffectiveGIMRVars_LFU_QFU[NNPEffectiveGIMRVars_LFU_QFU];

    /**
     * @brief Constructor.
     * @param[in] FlagLeptonUniversal_in an internal boolean flag that is true if assuming lepton flavour universality
     * @param[in] FlagQuarkUniversal_in an internal boolean flag that is true if assuming quark flavour universality
     */
    NPEffectiveGIMR(const bool FlagLeptonUniversal_in = false, const bool FlagQuarkUniversal_in = false);

    /**
     * @brief The post-update method for %NPEffectiveGIMR.
     * @details This method runs all the procedures that are need to be executed
     * after the model is successfully updated.
     * @return a boolean that is true if the execution is successful
     */
    virtual bool PostUpdate();

    /**
     * @brief A method to check if all the mandatory parameters for %NPEffectiveGIMR
     * have been provided in model initialization.
     * @param[in] DPars a map of the parameters that are being updated in the Monte Carlo run
     * (including parameters that are varied and those that are held constant)
     * @return a boolean that is true if the execution is successful
     */
    virtual bool CheckParameters(const std::map<std::string, double>& DPars);

    /**
     * @brief A method to set a flag of %NPEffectiveGIMR.
     * @param[in] name name of a model flag
     * @param[in] value the boolean to be assigned to the flag specified by name
     * @return a boolean that is true if the execution is successful
     */
    virtual bool setFlag(const std::string name, const bool value);


    ////////////////////////////////////////////////////////////////////////

    /**
     * @brief New physics contribution to the Fermi constant.
     * @details The new physics contribution @f$\Delta G@f$ is defined as
     * @f[
     * G_\mu = G_{\mu,\mathrm{SM}}(1+\Delta G)\,,
     * @f]
     * where @f$G_\mu@f$ is the experimental value measured through muon decays, 
     * and @f$G_{\mu,\mathrm{SM}}@f$ is the Fermi constant in the SM.
     * @return @f$\Delta G@f$
     */
    virtual double DeltaGF() const;

    /**
     * @brief The oblique parameter \f$S\f$.
     * @return the value of @f$S@f$
     */
    virtual double obliqueS() const;

    /**
     * @brief The oblique parameter \f$T\f$.
     * @return the value of @f$T@f$
     */
    virtual double obliqueT() const;

    /**
     * @brief The oblique parameter \f$U\f$.
     * @return the value of @f$U@f$
     */
    virtual double obliqueU() const;

    /**
     * @brief The mass of the @f$W@f$ boson, @f$M_W@f$.
     * @return @f$M_W@f$ in GeV
     */
    virtual double Mw() const;

    /**
     * @brief The total width of the @f$W@f$ boson, @f$\Gamma_W@f$.
     * @return @f$\Gamma_W@f$ in GeV
     */
    virtual double GammaW() const;

    /**
     * @brief New physics contribution to the neutral-current vector coupling @f$g_V^f@f$.
     * @param[in] f a lepton or quark
     * @return @f$\delta g_V^f@f$
     */
    virtual double deltaGV_f(const Particle p) const;

    /**
     * @brief New physics contribution to the neutral-current axial-vector coupling @f$g_A^f@f$.
     * @param[in] f a lepton or quark
     * @return @f$\delta g_A^f@f$
     */
    virtual double deltaGA_f(const Particle p) const;

    /**
     * @brief New physics contribution to the neutral-current left-handed coupling @f$g_L^f@f$.
     * @param[in] f a lepton or quark
     * @return @f$\delta g_L^f@f$
     */
        // no generation mixing
    double deltaGL_f(const Particle p) const;

    /**
     * @brief New physics contribution to the neutral-current right-handed coupling @f$g_R^f@f$.
     * @param[in] f a lepton or quark
     * @return @f$\delta g_R^f@f$
     */
        // no generation mixing
    double deltaGR_f(const Particle p) const;


        ////////////////////////////////////////////////////////////////////////
 
    /**
     * @brief The total cross section for @f$e^+e^- \to f\bar{f} @f$ [CHECK!]
     * @param[in] p a lepton or quark
     * @param[in] sqrt_s the center-of-mass energy in ???
     * @return @f$\sigma(e^+e^- \to \mathrm{hadrons})@f$ in ???
     */
    virtual double sigma_eeTOffbar(const Particle p, const double sqrt_s) const;
 
    /**
     * @brief The forward cross section for @f$e^+e^- \to f\bar{f} @f$ [CHECK!]
     * @param[in] p a lepton or quark
     * @param[in] sqrt_s the center-of-mass energy in ???
     * @return @f$\sigma_F(e^+e^- \to f\bar{f})@f$ in ???
     */
    virtual gslpp::complex sigma_eeTOffbarF(const Particle p, const double sqrt_s) const;
    /**
     * @brief The backward cross section for @f$e^+e^- \to f\bar{f} @f$ [CHECK!]
     * @param[in] p a lepton or quark
     * @param[in] sqrt_s the center-of-mass energy in ???
     * @return @f$\sigma_B(e^+e^- \to f\bar{f})@f$ in ???
     */
    virtual gslpp::complex sigma_eeTOffbarB(const Particle p, const double sqrt_s) const;
 
    /**
     * @brief The total cross section for @f$e^+e^- \to \mu^+ \mu^- @f$ [CHECK!]
     * @param[in] sqrt_s the center-of-mass energy in ???
     * @return @f$\sigma(e^+e^- \to \mu^+ \mu^-)@f$ in ???
     */
    double sigma_eeTOmumu(const double sqrt_s) const;
    
    /**
     * @brief The total cross section for @f$e^+e^- \to \mathrm{hadrons} @f$ [CHECK!]
     * @param[in] sqrt_s the center-of-mass energy in ???
     * @return @f$\sigma(e^+e^- \to \mathrm{hadrons})@f$ in ???
     */
    double sigma_eeTOqq(const double sqrt_s) const;

    /**
     * @brief The forward-backward asymmetry for @f$e^+e^- \to \mu^+ \mu^- @f$ [CHECK!]
     * @param[in] sqrt_s the center-of-mass energy in ???
     * @return @f$A_{FB}^{\mu}@f$
     */    
    double AFB_mu(const double sqrt_s) const;

    ////////////////////////////////////////////////////////////////////////
    
    ////////////////////////////////////////////////////////////////////////

    /**
     * @brief New physics contribution to the charged current coupling @f$W_\mu \bar{f_L}\gamma^mu f_L@f$.
     * @param[in] pbar a lepton or quark
     * @param[in] p a lepton or quark
     * @return @f$\delta g_{Wff}^{L}@f$
     */
    // no generation mixing
    gslpp::complex deltaGL_Wff(const Particle pbar, const Particle p) const;
    /**
     * @brief New physics contribution to the charged current coupling @f$W_\mu \bar{f_R}\gamma^mu f_R@f$.
     * @param[in] pbar a lepton or quark
     * @param[in] p a lepton or quark
     * @return @f$\delta g_{Wff}^{R}@f$
     */
    // no generation mixing
    gslpp::complex deltaGR_Wff(const Particle pbar, const Particle p) const;

    /**
     * @brief The new physics contribution to the coupling of the effective interaction @f$H G_{\mu\nu}^AG^{A \mu\nu}@f$.
     * @return @f$\delta g_{HGG}@f$
     */
    double deltaG_hgg() const;
    /**
     * @brief The new physics contribution to the coupling of the effective interaction @f$H W_{\mu\nu}^\dagger W^{\mu\nu}@f$.
     * @return @f$\delta g_{HWW}^{(1)}@f$
     */
    double deltaG1_hWW() const;
    /**
     * @brief The new physics contribution to the coupling of the effective interaction @f$H W_{\nu}^\dagger \partial^\mu W^{\mu\nu}@f$.
     * @return @f$\delta g_{HWW}^{(2)}@f$
     */
    double deltaG2_hWW() const;
    /**
     * @brief The new physics contribution to the coupling of the effective interaction @f$H W_{\mu}^\dagger W^{\mu}@f$.
     * @return @f$\delta g_{HWW}^{(3)}@f$
     */
    double deltaG3_hWW() const;
    /**
     * @brief The new physics contribution to the coupling of the effective interaction @f$H Z_{\mu\nu} Z^{\mu\nu}@f$.
     * @return @f$\delta g_{HZZ}^{(1)}@f$
     */
    double deltaG1_hZZ() const;
    /**
     * @brief The new physics contribution to the coupling of the effective interaction @f$H Z_{\nu} \partial^\mu Z^{\mu\nu}@f$.
     * @return @f$\delta g_{HZZ}^{(2)}@f$
     */
    double deltaG2_hZZ() const;
    /**
     * @brief The new physics contribution to the coupling of the effective interaction @f$H Z_{\mu} Z^{\mu}@f$.
     * @return @f$\delta g_{HZZ}^{(3)}@f$
     */
    double deltaG3_hZZ() const;
    /**
     * @brief The new physics contribution to the coupling of the effective interaction @f$H Z_{\mu\nu} F^{\mu\nu}@f$.
     * @return @f$\delta g_{HZA}^{(1)}@f$
     */
    double deltaG1_hZA() const;
    /**
     * @brief The new physics contribution to the coupling of the effective interaction @f$H Z_{\nu} \partial^\mu F^{\mu\nu}@f$.
     * @return @f$\delta g_{HZA}^{(2)}@f$
     */
    double deltaG2_hZA() const;
    /**
     * @brief The new physics contribution to the coupling of the effective interaction @f$H F_{\mu\nu} F^{\mu\nu}@f$.
     * @return @f$\delta g_{HAA}@f$
     */
    double deltaG_hAA() const;

    /**
     * @brief The new physics contribution to the coupling of the effective interaction @f$H f\bar{f}@f$.
     * @param[in] p a lepton or quark
     * @return @f$\delta g_{Hff}@f$
     */
    // no generation mixing
    gslpp::complex deltaG_hff(const Particle p) const;

    /**
     * @brief The new physics contribution to the coupling of the effective interaction @f$H W_\mu \bar{f_L}\gamma^mu f_L@f$.
     * @param[in] pbar a lepton or quark
     * @param[in] p a lepton or quark
     * @return @f$\delta g_{WffH}^{L}@f$
     */
    // no generation mixing
    gslpp::complex deltaGL_Wffh(const Particle pbar, const Particle p) const;
    
    /**
     * @brief The new physics contribution to the coupling of the effective interaction @f$H W_\mu \bar{f_R}\gamma^mu f_R@f$.
     * @param[in] pbar a lepton or quark
     * @param[in] p a lepton or quark
     * @return @f$\delta g_{WffH}^{R}@f$
     */
    // no generation mixing
    gslpp::complex deltaGR_Wffh(const Particle pbar, const Particle p) const;
    
    /**
     * @brief The new physics contribution to the coupling of the effective interaction @f$H Z_\mu \bar{f_L}\gamma^mu f_L@f$.
     * @param[in] p a lepton or quark
     * @return @f$\delta g_{ZffH}^{L}@f$
     */
    // no generation mixing
    double deltaGL_Zffh(const Particle p) const;
    
    /**
     * @brief The new physics contribution to the coupling of the effective interaction @f$H Z_\mu \bar{f_R}\gamma^mu f_R@f$.
     * @param[in] p a lepton or quark
     * @return @f$\delta g_{ZffH}^{R}@f$
     */
    // no generation mixing 
    double deltaGR_Zffh(const Particle p) const;


    ////////////////////////////////////////////////////////////////////////

    /**
     * @brief Loop function entering in the calculation of the effective @f$Hgg@f$ and @f$H\gamma\gamma@f$ couplings.
     * @param[in] @f$\tau=4 M^2/m_h^2@f$, with @f$M@f$ the mass of the particle in the loop.
     * @return @f$f(\tau)@f$
     */
    gslpp::complex f_triangle(const double tau) const;
    /**
     * @brief Fermionic loop function entering in the calculation of the effective @f$Hgg@f$ and @f$H\gamma\gamma@f$ couplings.
     * @details
     * @f$A^H_f(\tau)=2\tau [1+(1-\tau)f(\tau)]@f$
     * @param[in] @f$\tau=4 M^2/m_h^2@f$, with @f$M@f$ the mass of the fermion in the loop.
     * @return @f$A^H_f(\tau)@f$
     */
    gslpp::complex AH_f(const double tau) const;

    /**
     * @brief The ratio @f$\mu_{ggH}@f$ between the gluon-gluon fusion Higgs
     * production cross-section in the current model and in the Standard Model.
     * @param[in] sqrt_s the center-of-mass energy in TeV
     * @return @f$\mu_{ggH}@f$
     */
    virtual double muggH(const double sqrt_s) const;
    /**
     * @brief The ratio @f$\mu_{VBF}@f$ between the vector-boson fusion Higgs
     * production cross-section in the current model and in the Standard Model. 
     * @param[in] sqrt_s the center-of-mass energy in TeV
     * @return @f$\mu_{VBF}@f$
     */
    virtual double muVBF(const double sqrt_s) const;
    /**
     * @brief The ratio @f$\mu_{WH}@f$ between the W-Higgs associated production
     * cross-section in the current model and in the Standard Model. 
     * @param[in] sqrt_s the center-of-mass energy in TeV
     * @return @f$\mu_{WH}@f$
     */
    virtual double muWH(const double sqrt_s) const;
    /**
     * @brief The ratio @f$\mu_{ZH}@f$ between the Z-Higgs associated production
     * cross-section in the current model and in the Standard Model.
     * @param[in] sqrt_s the center-of-mass energy in TeV
     * @return @f$\mu_{ZH}@f$
     */
    virtual double muZH(const double sqrt_s) const;
    /**
     * @brief The ratio @f$\mu_{eeZH}@f$ between the 
     * @f$e^{+}e^{-}\to ZH@f$ associated production
     * cross-section in the current model and in the Standard Model.
     * @param[in] sqrt_s the center-of-mass energy in TeV
     * @return @f$\mu_{eeZH}@f$
     */
    virtual double mueeZH(const double sqrt_s) const;
    /**
     * @brief The ratio @f$\mu_{VH}@f$ between the WH+ZH associated production
     * cross-section in the current model and in the Standard Model.
     * @param[in] sqrt_s the center-of-mass energy in TeV
     * @return @f$\mu_{VH}@f$
     */
    virtual double muVH(const double sqrt_s) const;
    /**
     * @brief The ratio @f$\mu_{VBF+VH}@f$ between the sum of VBF and WH+ZH associated production
     * cross-section in the current model and in the Standard Model.
     * @param[in] sqrt_s the center-of-mass energy in TeV
     * @return @f$\mu_{VBF+VH}@f$
     */
    virtual double muVBFpVH(const double sqrt_s) const;
    /**
     * @brief The ratio @f$\mu_{ttH}@f$ between the t-tbar-Higgs associated 
     * production cross-section in the current model and in the Standard Model.
     * @param[in] sqrt_s the center-of-mass energy in TeV
     * @return @f$\mu_{ttH}@f$
     */
    virtual double muttH(const double sqrt_s) const;
    /**
     * @brief The ratio @f$\mu_{ggH+ttH}@f$ between the sum of gluon-gluon fusion
     * and t-tbar-Higgs associated 
     * production cross-section in the current model and in the Standard Model.
     * @param[in] sqrt_s the center-of-mass energy in TeV
     * @return @f$\mu_{ggH+ttH}@f$
     */
    virtual double muggHpttH(const double sqrt_s) const;
    /**
     * @brief The ratio of the Br@f$(H\to gg)@f$ in the current model
     * and in the Standard Model.
     * @return Br@f$(H\to gg)@f$/Br@f$(H\to gg)_{\mathrm{SM}}@f$
     */
    virtual double BrHggRatio() const;
    /**
     * @brief The ratio of the Br@f$(H\to WW)@f$ in the current model
     * and in the Standard Model.
     * @return Br@f$(H\to WW)@f$/Br@f$(H\to WW)_{\mathrm{SM}}@f$
     */
    virtual double BrHWWRatio() const;
    /**
     * @brief The ratio of the Br@f$(H\to ZZ)@f$ in the current model
     * and in the Standard Model.
     * @return Br@f$(H\to ZZ)@f$/Br@f$(H\to ZZ)_{\mathrm{SM}}@f$
     */
    virtual double BrHZZRatio() const;
    /**
     * @brief The ratio of the Br@f$(H\to Z\gamma)@f$ in the current model
     * and in the Standard Model.
     * @return Br@f$(H\to Z\gamma)@f$/Br@f$(H\to Z\gamma)_{\mathrm{SM}}@f$
     */
    virtual double BrHZgaRatio() const;
    /**
     * @brief The ratio of the Br@f$(H\to \gamma\gamma)@f$ in the current model
     * and in the Standard Model.
     * @return Br@f$(H\to \gamma\gamma)@f$/Br@f$(H\to \gamma\gamma)_{\mathrm{SM}}@f$
     */
    virtual double BrHgagaRatio() const;
    /**
     * @brief The ratio of the Br@f$(H\to \tau^+\tau^-)@f$ in the current model
     * and in the Standard Model.
     * @return Br@f$(H\to \tau^+\tau^-)@f$/Br@f$(H\to \tau^+\tau^-)_{\mathrm{SM}}@f$
     */
    virtual double BrHtautauRatio() const;
    /**
     * @brief The ratio of the Br@f$(H\to c\bar{c})@f$ in the current model
     * and in the Standard Model.
     * @return Br@f$(H\to c\bar{c})@f$/Br@f$(H\to c\bar{c})_{\mathrm{SM}}@f$
     */
    virtual double BrHccRatio() const;
    /**
     * @brief The ratio of the Br@f$(H\to b\bar{b})@f$ in the current model
     * and in the Standard Model.
     * @return Br@f$(H\to b\bar{b})@f$/Br@f$(H\to b\bar{b})_{\mathrm{SM}}@f$
     */
    virtual double BrHbbRatio() const;
    /**
     * @brief The ratio of the @f$\Gamma(H)@f$ in the current model
     * and in the Standard Model.
     * @return @f$\Gamma(H)@f$/@f$\Gamma(H)_{\mathrm{SM}}@f$
     */
    virtual double computeGammaTotalRatio() const;
    /**
     * @brief The new physics contribution to the ratio of the @f$\Gamma(H)@f$ 
     * in the current model and in the Standard Model. Only terms that are
     * linear in the effective Lagrangian coefficients.
     * @return @f$\delta \Gamma(H)@f$/@f$\Gamma(H)_{\mathrm{SM}}@f$
     */
    virtual double deltaGammaTotalRatio1() const;
    /**
     * @brief The new physics contribution to the ratio of the @f$\Gamma(H)@f$ 
     * in the current model and in the Standard Model. Only terms that are
     * quadratic in the effective Lagrangian coefficients.
     * @return @f$\delta \Gamma(H)@f$/@f$\Gamma(H)_{\mathrm{SM}}@f$
     */
    virtual double deltaGammaTotalRatio2() const;
    /**
     * @brief The ratio of the @f$\Gamma(H\to gg)@f$ in the current model
     * and in the Standard Model.
     * @return @f$\Gamma(H\to gg)@f$/@f$\Gamma(H\to gg)_{\mathrm{SM}}@f$
     */
    double GammaHggRatio() const;
    /**
     * @brief The ratio of the @f$\Gamma(H\to WW)@f$ in the current model
     * and in the Standard Model.
     * @return @f$\Gamma(H\to WW)@f$/@f$\Gamma(H\to WW)_{\mathrm{SM}}@f$
     */
    double GammaHWWRatio() const;
    /**
     * @brief The ratio of the @f$\Gamma(H\to ZZ)@f$ in the current model
     * and in the Standard Model.
     * @return @f$\Gamma(H\to ZZ)@f$/@f$\Gamma(H\to ZZ)_{\mathrm{SM}}@f$
     */
    double GammaHZZRatio() const;
    /**
     * @brief The ratio of the @f$\Gamma(H\to Z\gamma)@f$ in the current model
     * and in the Standard Model.
     * @return @f$\Gamma(H\to Z\gamma)@f$/@f$\Gamma(H\to Z\gamma)_{\mathrm{SM}}@f$
     */
    double GammaHZgaRatio() const;
    /**
     * @brief The ratio of the @f$\Gamma(H\to \gamma\gamma)@f$ in the current model
     * and in the Standard Model.
     * @return @f$\Gamma(H\to \gamma\gamma)@f$/@f$\Gamma(H\to \gamma\gamma)_{\mathrm{SM}}@f$
     */
    double GammaHgagaRatio() const;
    /**
     * @brief The ratio of the @f$\Gamma(H\to \tau\tau)@f$ in the current model
     * and in the Standard Model.
     * @return @f$\Gamma(H\to \tau\tau)@f$/@f$\Gamma(H\to \tau\tau)_{\mathrm{SM}}@f$
     */
    double GammaHtautauRatio() const;
    /**
     * @brief The ratio of the @f$\Gamma(H\to cc)@f$ in the current model
     * and in the Standard Model.
     * @return @f$\Gamma(H\to cc)@f$/@f$\Gamma(H\to cc)_{\mathrm{SM}}@f$
     */
    double GammaHccRatio() const;
    /**
     * @brief The ratio of the @f$\Gamma(H\to bb)@f$ in the current model
     * and in the Standard Model.
     * @return @f$\Gamma(H\to bb)@f$/@f$\Gamma(H\to bb)_{\mathrm{SM}}@f$
     */
    double GammaHbbRatio() const;
    /**
     * @brief The new physics contribution to the ratio of the @f$\Gamma(H\to gg)@f$
     * in the current model and in the Standard Model. Only terms that are
     * linear in the effective Lagrangian coefficients.
     * @return @f$\delta \Gamma(H\to gg)@f$/@f$\Gamma(H\to gg)_{\mathrm{SM}}@f$
     */
    double deltaGammaHggRatio1() const;
    /**
     * @brief The new physics contribution to the ratio of the @f$\Gamma(H\to WW)@f$
     * in the current model and in the Standard Model. (Only terms that are
     * linear in the effective Lagrangian coefficients.)
     * @return @f$\delta \Gamma(H\to WW)@f$/@f$\Gamma(H\to WW)_{\mathrm{SM}}@f$
     */
    double deltaGammaHWWRatio1() const;
    /**
     * @brief The new physics contribution to the ratio of the @f$\Gamma(H\to ZZ)@f$
     * in the current model and in the Standard Model. (Only terms that are
     * linear in the effective Lagrangian coefficients.)
     * @return @f$\delta \Gamma(H\to ZZ)@f$/@f$\Gamma(H\to ZZ)_{\mathrm{SM}}@f$
     */
    double deltaGammaHZZRatio1() const;
    /**
     * @brief The new physics contribution to the ratio of the @f$\Gamma(H\to Z\gamma)@f$
     * in the current model and in the Standard Model. (Only terms that are
     * linear in the effective Lagrangian coefficients.)
     * @return @f$\delta \Gamma(H\to Z\gamma)@f$/@f$\Gamma(H\to Z\gamma)_{\mathrm{SM}}@f$
     */
    double deltaGammaHZgaRatio1() const;
    /**
     * @brief The new physics contribution to the ratio of the @f$\Gamma(H\to \gamma\gamma)@f$
     * in the current model and in the Standard Model. (Only terms that are
     * linear in the effective Lagrangian coefficients.)
     * @return @f$\delta \Gamma(H\to \gamma\gamma)@f$/@f$\Gamma(H\to \gamma\gamma)_{\mathrm{SM}}@f$
     */
    double deltaGammaHgagaRatio1() const;
    /**
     * @brief The new physics contribution to the ratio of the @f$\Gamma(H\to \tau\tau)@f$
     * in the current model and in the Standard Model. (Only terms that are
     * linear in the effective Lagrangian coefficients.)
     * @return @f$\delta \Gamma(H\to \tau\tau)@f$/@f$\Gamma(H\to \tau\tau)_{\mathrm{SM}}@f$
     */
    double deltaGammaHtautauRatio1() const;
    /**
     * @brief The new physics contribution to the ratio of the @f$\Gamma(H\to cc)@f$
     * in the current model and in the Standard Model. (Only terms that are
     * linear in the effective Lagrangian coefficients.)
     * @return @f$\delta \Gamma(H\to cc)@f$/@f$\Gamma(H\to cc)_{\mathrm{SM}}@f$
     */
    double deltaGammaHccRatio1() const;
    /**
     * @brief The new physics contribution to the ratio of the @f$\Gamma(H\to bb)@f$
     * in the current model and in the Standard Model. (Only terms that are
     * linear in the effective Lagrangian coefficients.)
     * @return @f$\delta \Gamma(H\to bb)@f$/@f$\Gamma(H\to bb)_{\mathrm{SM}}@f$
     */
    double deltaGammaHbbRatio1() const;
    /**
     * @brief The new physics contribution to the ratio of the @f$\Gamma(H\to gg)@f$
     * in the current model and in the Standard Model. (Only terms that are
     * quadratic in the effective Lagrangian coefficients.)
     * @return @f$\delta \Gamma(H\to gg)@f$/@f$\Gamma(H\to gg)_{\mathrm{SM}}@f$
     */   
    double deltaGammaHggRatio2() const;
    /**
     * @brief The new physics contribution to the ratio of the @f$\Gamma(H\to WW)@f$
     * in the current model and in the Standard Model. (Only terms that are
     * quadratic in the effective Lagrangian coefficients.)
     * @return @f$\delta \Gamma(H\to WW)@f$/@f$\Gamma(H\to WW)_{\mathrm{SM}}@f$
     */
    double deltaGammaHWWRatio2() const;
    /**
     * @brief The new physics contribution to the ratio of the @f$\Gamma(H\to ZZ)@f$
     * in the current model and in the Standard Model. (Only terms that are
     * quadratic in the effective Lagrangian coefficients.)
     * @return @f$\delta \Gamma(H\to ZZ)@f$/@f$\Gamma(H\to ZZ)_{\mathrm{SM}}@f$
     */
    double deltaGammaHZZRatio2() const;
    /**
     * @brief The new physics contribution to the ratio of the @f$\Gamma(H\to Z\gamma)@f$
     * in the current model and in the Standard Model. (Only terms that are
     * quadratic in the effective Lagrangian coefficients.)
     * @return @f$\delta \Gamma(H\to Z\gamma)@f$/@f$\Gamma(H\to Z\gamma)_{\mathrm{SM}}@f$
     */
    double deltaGammaHZgaRatio2() const;
    /**
     * @brief The new physics contribution to the ratio of the @f$\Gamma(H\to \gamma\gamma)@f$
     * in the current model and in the Standard Model. (Only terms that are
     * quadratic in the effective Lagrangian coefficients.)
     * @return @f$\delta \Gamma(H\to \gamma\gamma)@f$/@f$\Gamma(H\to \gamma\gamma)_{\mathrm{SM}}@f$
     */
    double deltaGammaHgagaRatio2() const;
    /**
     * @brief The new physics contribution to the ratio of the @f$\Gamma(H\to \tau\tau)@f$
     * in the current model and in the Standard Model. (Only terms that are
     * quadratic in the effective Lagrangian coefficients.)
     * @return @f$\delta \Gamma(H\to \tau\tau)@f$/@f$\Gamma(H\to \tau\tau)_{\mathrm{SM}}@f$
     */
    double deltaGammaHtautauRatio2() const;
    /**
     * @brief The new physics contribution to the ratio of the @f$\Gamma(H\to cc)@f$
     * in the current model and in the Standard Model. (Only terms that are
     * quadratic in the effective Lagrangian coefficients.)
     * @return @f$\delta \Gamma(H\to cc)@f$/@f$\Gamma(H\to cc)_{\mathrm{SM}}@f$
     */
    double deltaGammaHccRatio2() const;
    /**
     * @brief The new physics contribution to the ratio of the @f$\Gamma(H\to bb)@f$
     * in the current model and in the Standard Model. (Only terms that are
     * quadratic in the effective Lagrangian coefficients.)
     * @return @f$\delta \Gamma(H\to bb)@f$/@f$\Gamma(H\to bb)_{\mathrm{SM}}@f$
     */
    double deltaGammaHbbRatio2() const;
    
    ////////////////////////////////////////////////////////////////////////
protected:

    /**
     * @brief @copybrief Model::setParameter()
     * @copydetails Model::setParameter()
     */
    virtual void setParameter(const std::string name, const double& value);

    double CW; ///< The dimension-6 operator coefficient \f$C_{W}\f$.
    double CHG; ///< The dimension-6 operator coefficient \f$C_{HG}\f$.
    double CHW; ///< The dimension-6 operator coefficient \f$C_{HW}\f$.
    double CHB; ///< The dimension-6 operator coefficient \f$C_{HB}\f$.
    double CHWB; ///< The dimension-6 operator coefficient \f$C_{HWB}\f$.
    double CHD; ///< The dimension-6 operator coefficient \f$C_{HD}\f$.
    double CHbox; ///< The dimension-6 operator coefficient \f$C_{H\Box}\f$.
    double CH; ///< The dimension-6 operator coefficient \f$C_{H}\f$.
    double CHL1_11; ///< The dimension-6 operator coefficient \f$(C_{HL}^{(1)})_{11}\f$.
    double CHL1_12r; ///< The dimension-6 operator coefficient \f$(C_{HL}^{(1)})_{12}\f$ (real part).
    double CHL1_13r; ///< The dimension-6 operator coefficient \f$(C_{HL}^{(1)})_{13}\f$ (real part).
    double CHL1_22; ///< The dimension-6 operator coefficient \f$(C_{HL}^{(1)})_{22}\f$.
    double CHL1_23r; ///< The dimension-6 operator coefficient \f$(C_{HL}^{(1)})_{23}\f$ (real part).
    double CHL1_33; ///< The dimension-6 operator coefficient \f$(C_{HL}^{(1)})_{33}\f$.
    double CHL1_12i; ///< The dimension-6 operator coefficient \f$(C_{HL}^{(1)})_{12}\f$ (imaginary part). 
    double CHL1_13i; ///< The dimension-6 operator coefficient \f$(C_{HL}^{(1)})_{13}\f$ (imaginary part).
    double CHL1_23i; ///< The dimension-6 operator coefficient \f$(C_{HL}^{(1)})_{23}\f$ (imaginary part).
    double CHL3_11; ///< The dimension-6 operator coefficient \f$(C_{HL}^{(3)})_{11}\f$.
    double CHL3_12r; ///< The dimension-6 operator coefficient \f$(C_{HL}^{(3)})_{12}\f$ (real part).
    double CHL3_13r; ///< The dimension-6 operator coefficient \f$(C_{HL}^{(3)})_{13}\f$ (real part).
    double CHL3_22; ///< The dimension-6 operator coefficient \f$(C_{HL}^{(3)})_{22}\f$.
    double CHL3_23r; ///< The dimension-6 operator coefficient \f$(C_{HL}^{(3)})_{23}\f$ (real part).
    double CHL3_33; ///< The dimension-6 operator coefficient \f$(C_{HL}^{(3)})_{33}\f$.
    double CHL3_12i; ///< The dimension-6 operator coefficient \f$(C_{HL}^{(3)})_{12}\f$ (real part).
    double CHL3_13i; ///< The dimension-6 operator coefficient \f$(C_{HL}^{(3)})_{13}\f$ (real part).
    double CHL3_23i; ///< The dimension-6 operator coefficient \f$(C_{HL}^{(3)})_{23}\f$ (real part).
    double CHe_11; ///< The dimension-6 operator coefficient \f$(C_{HE})_{11}\f$.
    double CHe_12r; ///< The dimension-6 operator coefficient \f$(C_{HE})_{12}\f$ (real part).
    double CHe_13r; ///< The dimension-6 operator coefficient \f$(C_{HE})_{13}\f$ (real part).
    double CHe_22; ///< The dimension-6 operator coefficient \f$(C_{HE})_{22}\f$.
    double CHe_23r; ///< The dimension-6 operator coefficient \f$(C_{HE})_{23}\f$ (real part).
    double CHe_33; ///< The dimension-6 operator coefficient \f$(C_{HE})_{33}\f$.
    double CHe_12i; ///< The dimension-6 operator coefficient \f$(C_{HE})_{12}\f$ (imaginary part).
    double CHe_13i; ///< The dimension-6 operator coefficient \f$(C_{HE})_{13}\f$ (imaginary part).
    double CHe_23i; ///< The dimension-6 operator coefficient \f$(C_{HE})_{23}\f$ (imaginary part).
    double CHQ1_11; ///< The dimension-6 operator coefficient \f$(C_{HQ}^{(1)})_{11}\f$.
    double CHQ1_12r; ///< The dimension-6 operator coefficient \f$(C_{HQ}^{(1)})_{12}\f$ (real part).
    double CHQ1_13r; ///< The dimension-6 operator coefficient \f$(C_{HQ}^{(1)})_{13}\f$ (real part).
    double CHQ1_22; ///< The dimension-6 operator coefficient \f$(C_{HQ}^{(1)})_{22}\f$.
    double CHQ1_23r; ///< The dimension-6 operator coefficient \f$(C_{HQ}^{(1)})_{23}\f$ (real part).
    double CHQ1_33; ///< The dimension-6 operator coefficient \f$(C_{HQ}^{(1)})_{33}\f$.
    double CHQ1_12i; ///< The dimension-6 operator coefficient \f$(C_{HQ}^{(1)})_{12}\f$ (imaginary part).
    double CHQ1_13i; ///< The dimension-6 operator coefficient \f$(C_{HQ}^{(1)})_{13}\f$ (imaginary part).
    double CHQ1_23i; ///< The dimension-6 operator coefficient \f$(C_{HQ}^{(1)})_{23}\f$ (imaginary part).
    double CHQ3_11; ///< The dimension-6 operator coefficient \f$(C_{HQ}^{(3)})_{11}\f$.
    double CHQ3_12r; ///< The dimension-6 operator coefficient \f$(C_{HQ}^{(3)})_{12}\f$ (real part).
    double CHQ3_13r; ///< The dimension-6 operator coefficient \f$(C_{HQ}^{(3)})_{13}\f$ (real part).
    double CHQ3_22; ///< The dimension-6 operator coefficient \f$(C_{HQ}^{(3)})_{22}\f$.
    double CHQ3_23r; ///< The dimension-6 operator coefficient \f$(C_{HQ}^{(3)})_{23}\f$ (real part).
    double CHQ3_33; ///< The dimension-6 operator coefficient \f$(C_{HQ}^{(3)})_{33}\f$.
    double CHQ3_12i; ///< The dimension-6 operator coefficient \f$(C_{HQ}^{(3)})_{12}\f$ (imaginary part).
    double CHQ3_13i; ///< The dimension-6 operator coefficient \f$(C_{HQ}^{(3)})_{13}\f$ (imaginary part).
    double CHQ3_23i; ///< The dimension-6 operator coefficient \f$(C_{HQ}^{(3)})_{23}\f$ (imaginary part).
    double CHu_11; ///< The dimension-6 operator coefficient \f$(C_{HU})_{11}\f$.
    double CHu_12r; ///< The dimension-6 operator coefficient \f$(C_{HU})_{12}\f$ (real part).
    double CHu_13r; ///< The dimension-6 operator coefficient \f$(C_{HU})_{13}\f$ (real part).
    double CHu_22; ///< The dimension-6 operator coefficient \f$(C_{HU})_{22}\f$.
    double CHu_23r; ///< The dimension-6 operator coefficient \f$(C_{HU})_{23}\f$ (real part).
    double CHu_33; ///< The dimension-6 operator coefficient \f$(C_{HU})_{33}\f$.
    double CHu_12i; ///< The dimension-6 operator coefficient \f$(C_{HU})_{12}\f$ (imaginary part).
    double CHu_13i; ///< The dimension-6 operator coefficient \f$(C_{HU})_{13}\f$ (imaginary part).
    double CHu_23i; ///< The dimension-6 operator coefficient \f$(C_{HU})_{23}\f$ (imaginary part).
    double CHd_11; ///< The dimension-6 operator coefficient \f$(C_{HD})_{11}\f$.
    double CHd_12r; ///< The dimension-6 operator coefficient \f$(C_{HD})_{12}\f$ (real part).
    double CHd_13r; ///< The dimension-6 operator coefficient \f$(C_{HD})_{13}\f$ (real part).
    double CHd_22; ///< The dimension-6 operator coefficient \f$(C_{HD})_{22}\f$.
    double CHd_23r; ///< The dimension-6 operator coefficient \f$(C_{HD})_{23}\f$ (real part).
    double CHd_33; ///< The dimension-6 operator coefficient \f$(C_{HD})_{33}\f$.
    double CHd_12i; ///< The dimension-6 operator coefficient \f$(C_{HD})_{12}\f$ (imaginary part).
    double CHd_13i; ///< The dimension-6 operator coefficient \f$(C_{HD})_{13}\f$ (imaginary part).
    double CHd_23i; ///< The dimension-6 operator coefficient \f$(C_{HD})_{23}\f$ (imaginary part).
    double CHud_11r; ///< The dimension-6 operator coefficient \f$(C_{HUD})_{11}\f$ (real part).
    double CHud_12r; ///< The dimension-6 operator coefficient \f$(C_{HUD})_{12}\f$ (real part).
    double CHud_13r; ///< The dimension-6 operator coefficient \f$(C_{HUD})_{13}\f$ (real part).
    double CHud_22r; ///< The dimension-6 operator coefficient \f$(C_{HUD})_{22}\f$ (real part).
    double CHud_23r; ///< The dimension-6 operator coefficient \f$(C_{HUD})_{23}\f$ (real part).
    double CHud_33r; ///< The dimension-6 operator coefficient \f$(C_{HUD})_{33}\f$ (real part).
    double CHud_11i; ///< The dimension-6 operator coefficient \f$(C_{HUD})_{11}\f$ (imaginary part).
    double CHud_12i; ///< The dimension-6 operator coefficient \f$(C_{HUD})_{12}\f$ (imaginary part).
    double CHud_13i; ///< The dimension-6 operator coefficient \f$(C_{HUD})_{13}\f$ (imaginary part).
    double CHud_22i; ///< The dimension-6 operator coefficient \f$(C_{HUD})_{22}\f$ (imaginary part).
    double CHud_23i; ///< The dimension-6 operator coefficient \f$(C_{HUD})_{23}\f$ (imaginary part).
    double CHud_33i; ///< The dimension-6 operator coefficient \f$(C_{HUD})_{33}\f$ (imaginary part).
    double CeH_11r; ///< The dimension-6 operator coefficient \f$(C_{EH})_{11}\f$ (real part).
    double CeH_12r; ///< The dimension-6 operator coefficient \f$(C_{EH})_{12}\f$ (real part).
    double CeH_13r; ///< The dimension-6 operator coefficient \f$(C_{EH})_{13}\f$ (real part).
    double CeH_22r; ///< The dimension-6 operator coefficient \f$(C_{EH})_{22}\f$ (real part).
    double CeH_23r; ///< The dimension-6 operator coefficient \f$(C_{EH})_{23}\f$ (real part).
    double CeH_33r; ///< The dimension-6 operator coefficient \f$(C_{EH})_{33}\f$ (real part).
    double CeH_11i; ///< The dimension-6 operator coefficient \f$(C_{EH})_{11}\f$ (imaginary part).
    double CeH_12i; ///< The dimension-6 operator coefficient \f$(C_{EH})_{12}\f$ (imaginary part).
    double CeH_13i; ///< The dimension-6 operator coefficient \f$(C_{EH})_{13}\f$ (imaginary part).
    double CeH_22i; ///< The dimension-6 operator coefficient \f$(C_{EH})_{22}\f$ (imaginary part).
    double CeH_23i; ///< The dimension-6 operator coefficient \f$(C_{EH})_{23}\f$ (imaginary part).
    double CeH_33i; ///< The dimension-6 operator coefficient \f$(C_{EH})_{33}\f$ (imaginary part).
    double CuH_11r; ///< The dimension-6 operator coefficient \f$(C_{UH})_{11}\f$ (real part).
    double CuH_12r; ///< The dimension-6 operator coefficient \f$(C_{UH})_{12}\f$ (real part).
    double CuH_13r; ///< The dimension-6 operator coefficient \f$(C_{UH})_{13}\f$ (real part).
    double CuH_22r; ///< The dimension-6 operator coefficient \f$(C_{UH})_{22}\f$ (real part).
    double CuH_23r; ///< The dimension-6 operator coefficient \f$(C_{UH})_{23}\f$ (real part).
    double CuH_33r; ///< The dimension-6 operator coefficient \f$(C_{UH})_{33}\f$ (real part).
    double CuH_11i; ///< The dimension-6 operator coefficient \f$(C_{UH})_{11}\f$ (imaginary part).
    double CuH_12i; ///< The dimension-6 operator coefficient \f$(C_{UH})_{12}\f$ (imaginary part).
    double CuH_13i; ///< The dimension-6 operator coefficient \f$(C_{UH})_{13}\f$ (imaginary part).
    double CuH_22i; ///< The dimension-6 operator coefficient \f$(C_{UH})_{22}\f$ (imaginary part).
    double CuH_23i; ///< The dimension-6 operator coefficient \f$(C_{UH})_{23}\f$ (imaginary part).
    double CuH_33i; ///< The dimension-6 operator coefficient \f$(C_{UH})_{33}\f$ (imaginary part).
    double CdH_11r; ///< The dimension-6 operator coefficient \f$(C_{DH})_{11}\f$ (real part).
    double CdH_12r; ///< The dimension-6 operator coefficient \f$(C_{DH})_{12}\f$ (real part).
    double CdH_13r; ///< The dimension-6 operator coefficient \f$(C_{DH})_{13}\f$ (real part).
    double CdH_22r; ///< The dimension-6 operator coefficient \f$(C_{DH})_{22}\f$ (real part).
    double CdH_23r; ///< The dimension-6 operator coefficient \f$(C_{DH})_{23}\f$ (real part).
    double CdH_33r; ///< The dimension-6 operator coefficient \f$(C_{DH})_{33}\f$ (real part).
    double CdH_11i; ///< The dimension-6 operator coefficient \f$(C_{DH})_{11}\f$ (imaginary part).
    double CdH_12i; ///< The dimension-6 operator coefficient \f$(C_{DH})_{12}\f$ (imaginary part).
    double CdH_13i; ///< The dimension-6 operator coefficient \f$(C_{DH})_{13}\f$ (imaginary part).
    double CdH_22i; ///< The dimension-6 operator coefficient \f$(C_{DH})_{22}\f$ (imaginary part).
    double CdH_23i; ///< The dimension-6 operator coefficient \f$(C_{DH})_{23}\f$ (imaginary part).
    double CdH_33i; ///< The dimension-6 operator coefficient \f$(C_{DH})_{33}\f$ (imaginary part).
    double CLL_1221; ///< The dimension-6 operator coefficient \f$(C_{LL})_{1221}\f$.
    double CLL_2112; ///< The dimension-6 operator coefficient \f$(C_{LL})_{2112}\f$.
    double CLQ1; ///< The dimension-6 (four-fermion) operator coefficient \f$C_{LQ}^{(1)}\f$.
    double CLQ3; ///< The dimension-6 (four-fermion) operator coefficient \f$C_{LQ}^{(3)}\f$.
    double Cee; ///< The dimension-6 (four-fermion) operator coefficient \f$C_{EE}\f$.
    double Ceu; ///< The dimension-6 (four-fermion) operator coefficient \f$C_{EU}\f$.
    double Ced; ///< The dimension-6 (four-fermion) operator coefficient \f$C_{ED}\f$.
    double CLe; ///< The dimension-6 (four-fermion) operator coefficient \f$C_{LE}\f$.
    double CLu; ///< The dimension-6 (four-fermion) operator coefficient \f$C_{LU}\f$.
    double CLd; ///< The dimension-6 (four-fermion) operator coefficient \f$C_{LD}\f$.
    double CQe; ///< The dimension-6 (four-fermion) operator coefficient \f$C_{QE}\f$.
    double Lambda_NP; ///< The new physics scale [GeV].
// The error in the parameters multiplying the dimension-6 operator coefficients in the production cross sections.
    double eVBF2_HZZ1;///< Theoretical uncertainty in the (linear) new physics contribution from \f$g_{HZZ}^{(1)}\f$ to VBF production at Tevatron (1.96 TeV).
    double eVBF2_HZZ2;///< Theoretical uncertainty in the (linear) new physics contribution from \f$g_{HZZ}^{(2)}\f$ to VBF production at Tevatron (1.96 TeV).
    double eVBF2_HZZ3;///< Theoretical uncertainty in the (linear) new physics contribution from \f$g_{HZZ}^{(3)}\f$ to VBF production at Tevatron (1.96 TeV).
    double eVBF2_HZA1;///< Theoretical uncertainty in the (linear) new physics contribution from \f$g_{HZA}^{(1)}\f$ to VBF production at Tevatron (1.96 TeV).
    double eVBF2_HZA2;///< Theoretical uncertainty in the (linear) new physics contribution from \f$g_{HZA}^{(2)}\f$ to VBF production at Tevatron (1.96 TeV).
    double eVBF2_HAA;///< Theoretical uncertainty in the (linear) new physics contribution from \f$g_{HAA}\f$ to VBF production at Tevatron (1.96 TeV).
    double eVBF2_HWW1;///< Theoretical uncertainty in the (linear) new physics contribution from \f$g_{HWW}^{(1)}\f$ to VBF production at Tevatron (1.96 TeV).
    double eVBF2_HWW2;///< Theoretical uncertainty in the (linear) new physics contribution from \f$g_{HWW}^{(2)}\f$ to VBF production at Tevatron (1.96 TeV).
    double eVBF2_HWW3;///< Theoretical uncertainty in the (linear) new physics contribution from \f$g_{HWW}^{(3)}\f$ to VBF production at Tevatron (1.96 TeV).
    double eVBF2_Hgg;///< Theoretical uncertainty in the (linear) new physics contribution from \f$g_{Hgg}\f$ to VBF production at Tevatron (1.96 TeV).
    double eVBF2_HZuL;///< Theoretical uncertainty in the (linear) new physics contribution from \f$g_{HZuu}^{L}\f$ to VBF production at Tevatron (1.96 TeV).
    double eVBF2_HZuR;///< Theoretical uncertainty in the (linear) new physics contribution from \f$g_{HZuu}^{R}\f$ to VBF production at Tevatron (1.96 TeV).
    double eVBF2_HZdL;///< Theoretical uncertainty in the (linear) new physics contribution from \f$g_{HZdd}^{L}\f$ to VBF production at Tevatron (1.96 TeV).
    double eVBF2_HZdR;///< Theoretical uncertainty in the (linear) new physics contribution from \f$g_{HZdd}^{R}\f$ to VBF production at Tevatron (1.96 TeV).
    double eVBF2_HWud;///< Theoretical uncertainty in the (linear) new physics contribution from \f$g_{HWud}^{L}\f$ to VBF production at Tevatron (1.96 TeV).
    double eVBF2_ZuL;///< Theoretical uncertainty in the (linear) new physics contribution from \f$g_{Zuu}^{L}\f$ to VBF production at Tevatron (1.96 TeV).
    double eVBF2_ZuR;///< Theoretical uncertainty in the (linear) new physics contribution from \f$g_{Zuu}^{R}\f$ to VBF production at Tevatron (1.96 TeV).
    double eVBF2_ZdL;///< Theoretical uncertainty in the (linear) new physics contribution from \f$g_{Zdd}^{L}\f$ to VBF production at Tevatron (1.96 TeV).
    double eVBF2_ZdR;///< Theoretical uncertainty in the (linear) new physics contribution from \f$g_{Zdd}^{R}\f$ to VBF production at Tevatron (1.96 TeV).
    double eVBF2_Wud;///< Theoretical uncertainty in the (linear) new physics contribution from \f$g_{Wud}^{L}\f$ to VBF production at Tevatron (1.96 TeV).
    double eVBF78_HZZ1;///< Theoretical uncertainty in the (linear) new physics contribution from \f$g_{HZZ}^{(1)}\f$ to VBF production at the LHC (7 & 8 TeV).
    double eVBF78_HZZ2;///< Theoretical uncertainty in the (linear) new physics contribution from \f$g_{HZZ}^{(2)}\f$ to VBF production at the LHC (7 & 8 TeV).
    double eVBF78_HZZ3;///< Theoretical uncertainty in the (linear) new physics contribution from \f$g_{HZZ}^{(3)}\f$ to VBF production at the LHC (7 & 8 TeV).
    double eVBF78_HZA1;///< Theoretical uncertainty in the (linear) new physics contribution from \f$g_{HZA}^{(1)}\f$ to VBF production at the LHC (7 & 8 TeV).
    double eVBF78_HZA2;///< Theoretical uncertainty in the (linear) new physics contribution from \f$g_{HZA}^{(2)}\f$ to VBF production at the LHC (7 & 8 TeV).
    double eVBF78_HAA;///< Theoretical uncertainty in the (linear) new physics contribution from \f$g_{HAA}\f$ to VBF production at the LHC (7 & 8 TeV).
    double eVBF78_HWW1;///< Theoretical uncertainty in the (linear) new physics contribution from \f$g_{HWW}^{(1)}\f$ to VBF production at the LHC (7 & 8 TeV).
    double eVBF78_HWW2;///< Theoretical uncertainty in the (linear) new physics contribution from \f$g_{HWW}^{(2)}\f$ to VBF production at the LHC (7 & 8 TeV).
    double eVBF78_HWW3;///< Theoretical uncertainty in the (linear) new physics contribution from \f$g_{HWW}^{(3)}\f$ to VBF production at the LHC (7 & 8 TeV).
    double eVBF78_Hgg;///< Theoretical uncertainty in the (linear) new physics contribution from \f$g_{Hgg}\f$ to VBF production at the LHC (7 & 8 TeV).
    double eVBF78_HZuL;///< Theoretical uncertainty in the (linear) new physics contribution from \f$g_{HZuu}^{L}\f$ to VBF production at the LHC (7 & 8 TeV).
    double eVBF78_HZuR;///< Theoretical uncertainty in the (linear) new physics contribution from \f$g_{HZuu}^{R}\f$ to VBF production at the LHC (7 & 8 TeV).
    double eVBF78_HZdL;///< Theoretical uncertainty in the (linear) new physics contribution from \f$g_{HZdd}^{L}\f$ to VBF production at the LHC (7 & 8 TeV).
    double eVBF78_HZdR;///< Theoretical uncertainty in the (linear) new physics contribution from \f$g_{HZdd}^{R}\f$ to VBF production at the LHC (7 & 8 TeV).
    double eVBF78_HWud;///< Theoretical uncertainty in the (linear) new physics contribution from \f$g_{HWud}^{L}\f$ to VBF production at the LHC (7 & 8 TeV).
    double eVBF78_ZuL;///< Theoretical uncertainty in the (linear) new physics contribution from \f$g_{Zuu}^{L}\f$ to VBF production at the LHC (7 & 8 TeV).
    double eVBF78_ZuR;///< Theoretical uncertainty in the (linear) new physics contribution from \f$g_{Zuu}^{R}\f$ to VBF production at the LHC (7 & 8 TeV).
    double eVBF78_ZdL;///< Theoretical uncertainty in the (linear) new physics contribution from \f$g_{Zdd}^{L}\f$ to VBF production at the LHC (7 & 8 TeV).
    double eVBF78_ZdR;///< Theoretical uncertainty in the (linear) new physics contribution from \f$g_{Zdd}^{R}\f$ to VBF production at the LHC (7 & 8 TeV).
    double eVBF78_Wud;///< Theoretical uncertainty in the (linear) new physics contribution from \f$g_{Wud}^{L}\f$ to VBF production at the LHC (7 & 8 TeV).
    
    double eWH2_HWW1;///< Theoretical uncertainty in the (linear) new physics contribution from \f$g_{HWW}^{(1)}\f$ to WH production at Tevatron (1.96 TeV).
    double eWH2_HWW2;///< Theoretical uncertainty in the (linear) new physics contribution from \f$g_{HWW}^{(2)}\f$ to WH production at Tevatron (1.96 TeV).
    double eWH2_HWW3;///< Theoretical uncertainty in the (linear) new physics contribution from \f$g_{HWW}^{(3)}\f$ to WH production at Tevatron (1.96 TeV).
    double eWH2_HWud;///< Theoretical uncertainty in the (linear) new physics contribution from \f$g_{HWud}^{L}\f$ to WH production at Tevatron (1.96 TeV).
    double eWH2_Wud;///< Theoretical uncertainty in the (linear) new physics contribution from \f$g_{Wud}^{L}\f$ to WH production at Tevatron (1.96 TeV).
    double eWH78_HWW1;///< Theoretical uncertainty in the (linear) new physics contribution from \f$g_{HWW}^{(1)}\f$ to WH production at the LHC (7 & 8 TeV).
    double eWH78_HWW2;///< Theoretical uncertainty in the (linear) new physics contribution from \f$g_{HWW}^{(2)}\f$ to WH production at the LHC (7 & 8 TeV).
    double eWH78_HWW3;///< Theoretical uncertainty in the (linear) new physics contribution from \f$g_{HWW}^{(3)}\f$ to WH production at the LHC (7 & 8 TeV).
    double eWH78_HWud;///< Theoretical uncertainty in the (linear) new physics contribution from \f$g_{HWud}^{L}\f$ to WH production at the LHC (7 & 8 TeV).
    double eWH78_Wud;///< Theoretical uncertainty in the (linear) new physics contribution from \f$g_{Wud}^{L}\f$ to WH production at the LHC (7 & 8 TeV).
    
    double eZH2_HZZ1;///< Theoretical uncertainty in the (linear) new physics contribution from \f$g_{HZZ}^{(1)}\f$ to ZH production at Tevatron (1.96 TeV).
    double eZH2_HZZ2;///< Theoretical uncertainty in the (linear) new physics contribution from \f$g_{HZZ}^{(2)}\f$ to ZH production at Tevatron (1.96 TeV).
    double eZH2_HZZ3;///< Theoretical uncertainty in the (linear) new physics contribution from \f$g_{HZZ}^{(3)}\f$ to ZH production at Tevatron (1.96 TeV).
    double eZH2_HZA1;///< Theoretical uncertainty in the (linear) new physics contribution from \f$g_{HZA}^{(1)}\f$ to ZH production at Tevatron (1.96 TeV).
    double eZH2_HZA2;///< Theoretical uncertainty in the (linear) new physics contribution from \f$g_{HZA}^{(2)}\f$ to ZH production at Tevatron (1.96 TeV).
    double eZH2_HZuL;///< Theoretical uncertainty in the (linear) new physics contribution from \f$g_{HZuu}^{L}\f$ to ZH production at Tevatron (1.96 TeV).
    double eZH2_HZuR;///< Theoretical uncertainty in the (linear) new physics contribution from \f$g_{HZuu}^{R}\f$ to ZH production at Tevatron (1.96 TeV).
    double eZH2_HZdL;///< Theoretical uncertainty in the (linear) new physics contribution from \f$g_{HZdd}^{L}\f$ to ZH production at Tevatron (1.96 TeV).
    double eZH2_HZdR;///< Theoretical uncertainty in the (linear) new physics contribution from \f$g_{HZdd}^{R}\f$ to ZH production at Tevatron (1.96 TeV).
    double eZH2_ZuL;///< Theoretical uncertainty in the (linear) new physics contribution from \f$g_{Zuu}^{L}\f$ to ZH production at Tevatron (1.96 TeV).
    double eZH2_ZuR;///< Theoretical uncertainty in the (linear) new physics contribution from \f$g_{Zuu}^{R}\f$ to ZH production at Tevatron (1.96 TeV).
    double eZH2_ZdL;///< Theoretical uncertainty in the (linear) new physics contribution from \f$g_{Zdd}^{L}\f$ to ZH production at Tevatron (1.96 TeV).
    double eZH2_ZdR;///< Theoretical uncertainty in the (linear) new physics contribution from \f$g_{Zdd}^{R}\f$ to ZH production at Tevatron (1.96 TeV).
    double eZH78_HZZ1;///< Theoretical uncertainty in the (linear) new physics contribution from \f$g_{HZZ}^{(1)}\f$ to ZH production at the LHC (7 & 8 TeV).
    double eZH78_HZZ2;///< Theoretical uncertainty in the (linear) new physics contribution from \f$g_{HZZ}^{(2)}\f$ to ZH production at the LHC (7 & 8 TeV).
    double eZH78_HZZ3;///< Theoretical uncertainty in the (linear) new physics contribution from \f$g_{HZZ}^{(3)}\f$ to ZH production at the LHC (7 & 8 TeV).
    double eZH78_HZA1;///< Theoretical uncertainty in the (linear) new physics contribution from \f$g_{HZA}^{(1)}\f$ to ZH production at the LHC (7 & 8 TeV).
    double eZH78_HZA2;///< Theoretical uncertainty in the (linear) new physics contribution from \f$g_{HZA}^{(2)}\f$ to ZH production at the LHC (7 & 8 TeV).
    double eZH78_HZuL;///< Theoretical uncertainty in the (linear) new physics contribution from \f$g_{HZuu}^{L}\f$ to ZH production at the LHC (7 & 8 TeV).
    double eZH78_HZuR;///< Theoretical uncertainty in the (linear) new physics contribution from \f$g_{HZuu}^{R}\f$ to ZH production at the LHC (7 & 8 TeV).
    double eZH78_HZdL;///< Theoretical uncertainty in the (linear) new physics contribution from \f$g_{HZdd}^{L}\f$ to ZH production at the LHC (7 & 8 TeV).
    double eZH78_HZdR;///< Theoretical uncertainty in the (linear) new physics contribution from \f$g_{HZdd}^{R}\f$ to ZH production at the LHC (7 & 8 TeV).
    double eZH78_ZuL;///< Theoretical uncertainty in the (linear) new physics contribution from \f$g_{Zuu}^{L}\f$ to ZH production at the LHC (7 & 8 TeV).
    double eZH78_ZuR;///< Theoretical uncertainty in the (linear) new physics contribution from \f$g_{Zuu}^{R}\f$ to ZH production at the LHC (7 & 8 TeV).
    double eZH78_ZdL;///< Theoretical uncertainty in the (linear) new physics contribution from \f$g_{Zdd}^{L}\f$ to ZH production at the LHC (7 & 8 TeV).
    double eZH78_ZdR;///< Theoretical uncertainty in the (linear) new physics contribution from \f$g_{Zdd}^{R}\f$ to ZH production at the LHC (7 & 8 TeV).
    
    double ettH2_Htt;///< Theoretical uncertainty in the (linear) new physics contribution from \f$g_{Htt}\f$ to ttH production at Tevatron (1.96 TeV).
    double ettH2_Hgg;///< Theoretical uncertainty in the (linear) new physics contribution from \f$g_{Hgg}\f$ to ttH production at Tevatron (1.96 TeV).
    double ettH78_Htt;///< Theoretical uncertainty in the (linear) new physics contribution from \f$g_{Htt}\f$ to ttH production at the LHC (7 & 8 TeV).
    double ettH78_Hgg;///< Theoretical uncertainty in the (linear) new physics contribution from \f$g_{Hgg}\f$ to ttH production at the LHC (7 & 8 TeV).
    
    double MwInput;///< The input value for the \f$W\f$-boson mass if FlagMwInput is true.

    double LambdaNP2;///< The square of the new physics scale [GeV\f$^2\f$].
    double v2_over_LambdaNP2;///< The ratio between the EW vev and the new physics scale, squared \f$v^2/\Lambda^2\f$.
    double cW_tree;///< The tree level values for the cosine of the weak angle.
    double sW_tree;///< The tree level values for the sine of the weak angle.
    double cW2_tree;///< The sqaure of the tree level values for the cosine of the weak angle.
    double sW2_tree;///< The sqaure of the tree level values for the sine of the weak angle.
    double delta_ZZ;///< Combination of dimension 6 coefficients modifying the \f$Z_\mu\f$ canonical field definition.
    double delta_AA;///< Combination of dimension 6 coefficients modifying the \f$A_\mu\f$ canonical field definition.
    double delta_AZ;///< Combination of dimension 6 coefficients modifying the \f$A_\mu\f$ canonical field definition.
    double delta_h;///< Combinations of dimension 6 coefficients modifying the \f$H\f$ canonical field definition.

    /**
     * @brief The diagonal entry of the dimension-6 operator coefficient \f$C_{HL,HQ}^{(1)}\f$ corresponding to particle F.
     * @param[in] F a lepton or quark
     * @return @f$(\f$C_{HF}^{(1)})_{FF}\f$@f$
     */
    double CHF1_diag(const Particle F) const;
    
    /**
     * @brief The diagonal entry of the dimension-6 operator coefficient \f$C_{HL,HQ}^{(3)}\f$ corresponding to particle F.
     * @param[in] F a lepton or quark
     * @return @f$(\f$C_{HF}^{(3)})_{FF}\f$@f$
     */
    double CHF3_diag(const Particle F) const;

    /**
     * @brief The diagonal entry of the dimension-6 operator coefficient \f$C_{HE,HU,HD}\f$ corresponding to particle f.
     * @param[in] f a lepton or quark
     * @return @f$(\f$C_{Hf})_{ff}\f$@f$
     */
    double CHf_diag(const Particle f) const;
    
    /**
     * @brief The diagonal entry of the dimension-6 operator coefficient \f$C_{HUD}\f$ corresponding to particle f.
     * @param[in] u a quark
     * @return @f$(\f$C_{HUD})_{ud}\f$@f$
     */
    gslpp::complex CHud_diag(const Particle u) const;
    
    /**
     * @brief The diagonal entry of the dimension-6 operator coefficient \f$C_{EH,UH,DH}\f$ corresponding to particle f.
     * @param[in] f a lepton or quark
     * @return @f$(\f$C_{fH})_{ff}\f$@f$
     */
    gslpp::complex CfH_diag(const Particle f) const;


    ////////////////////////////////////////////////////////////////////////
private:

    bool FlagMwInput; ///< A boolean flag that is true if the W mass is taken as an input parameter. (Warning: The W width is not implemented in this case.)
    bool FlagQuadraticTerms; ///< A boolean flag that is true if the quadratic terms in cross sections and widths are switched on.

    /**
     * @brief An internal boolean flag that is true if assuming lepton flavour
     * universality.
     */
    const bool FlagLeptonUniversal;

    /**
     * @brief An internal boolean flag that is true if assuming quark flavour
     * universality.
     */
    const bool FlagQuarkUniversal;

};

#endif	/* NPEFFECTIVEGIMR_H */

