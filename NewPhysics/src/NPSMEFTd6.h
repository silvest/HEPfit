/*
 * Copyright (C) 2014 HEPfit Collaboration
 *
 *
 * For the licensing terms see doc/COPYING.
 */

#ifndef NPSMEFTD6_H
#define	NPSMEFTD6_H

#include "gslpp.h"
#include "NPbase.h"
#include "NPSMEFTd6Matching.h"
#include <string.h>
#include <stdexcept>


/**
 * @class NPSMEFTd6
 * @brief A model class for new physics in the form of the
 * dimension-six effective Lagrangian.
 * @ingroup NewPhysics
 * @author HEPfit Collaboration
 * @copyright GNU General Public License
 * @details This is a Model class containing parameters and functions
 * associated with the general dimension-six effective Lagrangian.
 * (Use the model name "NPSMEFTd6_LFU_QFU" to asumme lepton and quark flavour universality)
 *
 * In this class we consider the dimension-six effective Lagrangian
 * @f[
 * \mathcal{L}_\mathrm{eff}
 * = \mathcal{L}_\mathrm{SM}
 *   + \sum_i \frac{C_i}{\Lambda^2} \mathcal{O}_i.
 * @f]
 * The implementation is written in the basis of \cite Grzadkowski:2010es.
 * For convenience, the parameterization also includes operators appearing in
 * other common bases. In particular, the complete set of parameters containts 2
 * redundancies, given by coefficients \f$C_{DHB} \f$ and \f$C_{DHW} \f$,
 * which correspond to operators not included in the basis of \cite Grzadkowski:2010es.
 * For meaningful physical results one must make sure to include only
 * a complete set of interactions in a given analysis.
 *
 * @anchor NPSMEFTd6Initialization
 * <h3>Initialization</h3>
 *
 * After creating an instance of the current class with the constructor
 * NPSMEFTd6(), it is required to call the initialization method
 * InitializeModel().
 * In the Monte Carlo run, the constructor as well as the initialization
 * method are called in InputParser::ReadParameters().
 *
 *
 * @anchor NPSMEFTd6Parameters
 * <h3>%Model parameters</h3>
 *
 * The model parameters of %NPSMEFTd6 are summarized below:
 * <table class="model">
 * <tr>
 *   <th>Label</th>
 *   <th>LaTeX symbol</th>
 *   <th>Description</th>
 * </tr>
 * <tr>
 *   <td class="mod_name">%CG </td>
 *   <td class="mod_symb">\f$C_{G} \f$</td>
 *   <td class="mod_desc">The coefficient of the operator
 *   \f${\cal O}_{G}=f_{ABC}G_{\mu}^{A\nu} G_{\nu}^{B\rho}W_{\rho}^{C\mu}\f$. </td>
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
 *   <td class="mod_name">%CDHB </td>
 *   <td class="mod_symb">\f$C_{DHB} \f$</td>
 *   <td class="mod_desc">The coefficient of the operator
 *   \f${\cal O}_{DHB}=i\big(D^\mu H^\dagger D^\nu H\big) B_{\mu\nu}\f$. </td>
 * </tr>
 * <tr>
 *   <td class="mod_name">%CDHW </td>
 *   <td class="mod_symb">\f$C_{DHW}\f$</td>
 *   <td class="mod_desc">The coefficient of the operator
 *   \f${\cal O}_{DHW}=i\big(D^\mu H^\dagger \tau^a D^\nu H\big) W_{\mu\nu}^a\f$. </td>
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
 *   \f$({\cal O}_{HL}^{(1)})_{ij} =i\big(H^\dagger \overset{\leftrightarrow}{D}_\mu H\big)
 *   \big(\overline{L^i}\,\gamma^\mu L^j\big)\f$, for \f$i,j=1,2,3\f$. </td>
 * </tr>
 * <tr>
 *   <td class="mod_name">%CHL3_kk, CHL3_klr, CHL3_kli </td>
 *   <td class="mod_symb">\f$ (C_{HL}^{(3)})_{kk}, \mbox{Re}\big[(C_{HL}^{(3)})_{kl}\big], \mbox{Im}\big[(C_{HL}^{(3)})_{kl}\big] \f$</td>
 *   <td class="mod_desc">The real and imaginary parts of the coefficient of the operator
 *   \f$({\cal O}_{HL}^{(3)})_{ij} =i\big(H^\dagger \overset{\leftrightarrow}{D^a_\mu} H\big)
 *   \big(\overline{L^i}\,\gamma^\mu \tau^a L^j\big)\f$, for \f$i,j=1,2,3\f$. </td>
 * </tr>
 * <tr>
 *   <td class="mod_name">%CHe_kk, CHe_klr, CHe_kli </td>
 *   <td class="mod_symb">\f$ (C_{He})_{kk}, \mbox{Re}\big[(C_{He})_{kl}\big], \mbox{Im}\big[(C_{He})_{kl}\big] \f$</td>
 *   <td class="mod_desc">The real and imaginary parts of the coefficient of the operator
 *   \f$({\cal O}_{He})_{ij} =i\big(H^\dagger \overset{\leftrightarrow}{D}_\mu H\big)
 *   \big(\overline{E^i}\,\gamma^\mu E^j\big)\f$, for \f$i,j=1,2,3\f$. </td>
 * </tr>
 * <tr>
 *   <td class="mod_name">%CHQ1_kk, CHQ1_klr, CHQ1_kli </td>
 *   <td class="mod_symb">\f$ (C_{HQ}^{(1)})_{kk}, \mbox{Re}\big[(C_{HQ}^{(1)})_{kl}\big], \mbox{Im}\big[(C_{HQ}^{(1)})_{kl}\big] \f$</td>
 *   <td class="mod_desc">The real and imaginary parts of the coefficient of the operator
 *   \f$({\cal O}_{HQ}^{(1)})_{ij} =i\big(H^\dagger \overset{\leftrightarrow}{D}_\mu H\big)
 *   \big(\overline{Q^i}\,\gamma^\mu Q^j\big)\f$, for \f$i,j=1,2,3\f$. </td>
 * </tr>
 * <tr>
 *   <td class="mod_name">%CHQ3_kk, CHQ3_klr, CHQ3_kli </td>
 *   <td class="mod_symb">\f$ (C_{HQ}^{(3)})_{kk}, \mbox{Re}\big[(C_{HQ}^{(3)})_{kl}\big], \mbox{Im}\big[(C_{HQ}^{(3)})_{kl}\big] \f$</td>
 *   <td class="mod_desc">The real and imaginary parts of the coefficient of the operator
 *   \f$({\cal O}_{HQ}^{(3)})_{ij} =i\big(H^\dagger \overset{\leftrightarrow}{D^a_\mu} H\big)
 *   \big(\overline{Q^i}\,\gamma^\mu \tau^a Q^j\big)\f$, for \f$i,j=1,2,3\f$. </td>
 * </tr>
 * <tr>
 *   <td class="mod_name">%CHu_kk, CHu_klr, CHu_kli </td>
 *   <td class="mod_symb">\f$ (C_{Hu})_{kk}, \mbox{Re}\big[(C_{Hu})_{kl}\big], \mbox{Im}\big[(C_{Hu})_{kl}\big] \f$</td>
 *   <td class="mod_desc">The real and imaginary parts of the coefficient of the operator
 *   \f$({\cal O}_{Hu})_{ij} =i\big(H^\dagger \overset{\leftrightarrow}{D}_\mu H\big)
 *   \big(\overline{U^i}\,\gamma^\mu U^j\big)\f$, for \f$i,j=1,2,3\f$. </td>
 * </tr>
 * <tr>
 *   <td class="mod_name">%CHd_kk, CHd_klr, CHd_kli </td>
 *   <td class="mod_symb">\f$ (C_{Hd})_{kk}, \mbox{Re}\big[(C_{Hd})_{kl}\big], \mbox{Im}\big[(C_{Hd})_{kl}\big] \f$</td>
 *   <td class="mod_desc">The real and imaginary parts of the coefficient of the operator
 *   \f$({\cal O}_{Hd})_{ij} =i\big(H^\dagger \overset{\leftrightarrow}{D}_\mu H\big)
 *   \big(\overline{D^i}\,\gamma^\mu D^j\big)\f$, for \f$i,j=1,2,3\f$. </td>
 * </tr>
  * <tr>
 *   <td class="mod_name">%CHud_klr, CHud_kli </td>
 *   <td class="mod_symb">\f$\mbox{Re}\big[(C_{Hud})_{kl}\big], \mbox{Im}\big[(C_{Hud})_{kl}\big] \f$</td>
 *   <td class="mod_desc">The real and imaginary parts of the coefficient of the operator
 *   \f$({\cal O}_{Hud})_{ij} =i\big(\widetilde{H}^\dagger D_\mu H\big)
 *   \big(\overline{U^i}\,\gamma^\mu D^j\big)\f$, for \f$i,j=1,2,3\f$. </td>
 * </tr>
 * <tr>
 *   <td class="mod_name">%CeH_klr, CeH_kli </td>
 *   <td class="mod_symb">\f$\mbox{Re}\big[(C_{eH})_{kl}\big], \mbox{Im}\big[(C_{eH})_{kl}\big] \f$</td>
 *   <td class="mod_desc">The real and imaginary parts of the coefficient of the operator
 *   \f$({\cal O}_{eH})_{ij} =\big(H^\dagger H\big)
 *   \big(\overline{L^i}\,H E^j\big)\f$, for \f$i,j=1,2,3\f$. </td>
 * </tr>
 * <tr>
 *   <td class="mod_name">%CuH_klr, CuH_kli </td>
 *   <td class="mod_symb">\f$\mbox{Re}\big[(C_{uH})_{kl}\big], \mbox{Im}\big[(C_{uH})_{kl}\big] \f$</td>
 *   <td class="mod_desc">The real and imaginary parts of the coefficient of the operator
 *   \f$({\cal O}_{uH})_{ij} =\big(H^\dagger H\big)
 *   \big(\overline{Q^i}\,\widetilde{H} U^j\big)\f$, for \f$i,j=1,2,3\f$. </td>
 * </tr>
 * <tr>
 *   <td class="mod_name">%CdH_klr, CdH_kli </td>
 *   <td class="mod_symb">\f$\mbox{Re}\big[(C_{dH})_{kl}\big], \mbox{Im}\big[(C_{dH})_{kl}\big] \f$</td>
 *   <td class="mod_desc">The real and imaginary parts of the coefficient of the operator
 *   \f$({\cal O}_{dH})_{ij} =\big(H^\dagger H\big)
 *   \big(\overline{Q^i}\,H D^j\big)\f$, for \f$i,j=1,2,3\f$. </td>
 * </tr>
 * <tr>
 *   <td class="mod_name">%CuG_klr, CuG_kli </td>
 *   <td class="mod_symb">\f$\mbox{Re}\big[(C_{uG})_{kl}\big], \mbox{Im}\big[(C_{uG})_{kl}\big] \f$</td>
 *   <td class="mod_desc">The real and imaginary parts of the coefficient of the operator
 *   \f$({\cal O}_{uG})_{ij} =\big(\overline{Q^i}\sigma^{\mu\nu} T_A U^j\big)\widetilde{H} G_{\mu\nu}^A\f$, for \f$i,j=1,2,3\f$. </td>
 * </tr>
 * <tr>
 *   <td class="mod_name">%CuW_klr, CuW_kli </td>
 *   <td class="mod_symb">\f$\mbox{Re}\big[(C_{uW})_{kl}\big], \mbox{Im}\big[(C_{uW})_{kl}\big] \f$</td>
 *   <td class="mod_desc">The real and imaginary parts of the coefficient of the operator
 *   \f$({\cal O}_{uW})_{ij} =\big(\overline{Q^i}\sigma^{\mu\nu} \tau_a U^j\big)\widetilde{H} W_{\mu\nu}^a\f$, for \f$i,j=1,2,3\f$. </td>
 * </tr>
 * <tr>
 *   <td class="mod_name">%CuB_klr, CuB_kli </td>
 *   <td class="mod_symb">\f$\mbox{Re}\big[(C_{uB})_{kl}\big], \mbox{Im}\big[(C_{uB})_{kl}\big] \f$</td>
 *   <td class="mod_desc">The real and imaginary parts of the coefficient of the operator
 *   \f$({\cal O}_{uB})_{ij} =\big(\overline{Q^i}\sigma^{\mu\nu} U^j\big)\widetilde{H} B_{\mu\nu}\f$, for \f$i,j=1,2,3\f$. </td>
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
 *   <td class="mod_name">%eXint </td>
 *   <td class="mod_symb">\f$\varepsilon_{X}^{int}\f$</td>
 *   <td class="mod_desc">The relative intrinsic theoretical uncertainty for the process X.
 *    (Only for Higgs observables and assumed to be constant in the energy.)</td>
 * </tr>
 * <tr>
 *   <td class="mod_name">%eXpar </td>
 *   <td class="mod_symb">\f$\varepsilon_{X}^{par}\f$</td>
 *   <td class="mod_desc">The relative parametric theoretical uncertainty for the process X.
 *    (Only for Higgs observables and assumed to be constant in the energy.)</td>
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
 * Where the hermitian derivatives are defined as 
 * @f[
 * H^\dagger i \overset{\leftrightarrow}{D}_\mu H\equiv H^\dagger i(D_\mu - \overset{\leftarrow}{D}_\mu)H
 * @f]
 * and
 * @f[
 * H^\dagger i \overset{\leftrightarrow}{D^a_\mu} H\equiv H^\dagger i (\tau^a D_\mu - \overset{\leftarrow}{D}_\mu \tau^a)H.  
 * @f]
 *
 * Alternatively, when using the model name "NPSMEFTd6_LFU_QFU", where lepton
 * and quark flavour universality are assumed 
 * (except for \f$(C_{fH})_{ij},~f=e,u,d\f$ which are assumed to be diagonal only), 
 * the parameters to be used as inputs
 * for the dimension six coefficients are the following:
 * <table class="model">
 * <tr>
 *   <th>Label</th>
 *   <th>LaTeX symbol</th>
 *   <th>Description</th>
 * </tr>
 * <tr>
 *   <td class="mod_name">%CG </td>
 *   <td class="mod_symb">\f$C_{G} \f$</td>
 *   <td class="mod_desc">The coefficient of the operator
 *   \f${\cal O}_{G}=f_{ABC}G_{\mu}^{A\nu} G_{\nu}^{B\rho}W_{\rho}^{C\mu}\f$. </td>
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
 *   <td class="mod_name">%CDHB </td>
 *   <td class="mod_symb">\f$C_{DHB} \f$</td>
 *   <td class="mod_desc">The coefficient of the operator
 *   \f${\cal O}_{DHB}=i\big(D^\mu H^\dagger D^\nu H\big) B_{\mu\nu}\f$. </td>
 * </tr>
 * <tr>
 *   <td class="mod_name">%CDHW </td>
 *   <td class="mod_symb">\f$C_{DHW}\f$</td>
 *   <td class="mod_desc">The coefficient of the operator
 *   \f${\cal O}_{DHW}=i\big(D^\mu H^\dagger \tau^a D^\nu H\big) W_{\mu\nu}^a\f$. </td>
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
 *   \f$({\cal O}_{HL}^{(1)})_{ii} =i\big(H^\dagger \overset{\leftrightarrow}{D}_\mu H\big)
 *   \big(\overline{L^i}\,\gamma^\mu L^i\big)\f$ (flavor universal). </td>
 * </tr>
 * <tr>
 *   <td class="mod_name">%CHL3 </td>
 *   <td class="mod_symb">\f$ (C_{HL}^{(3)})_{ii} \f$</td>
 *   <td class="mod_desc">The coefficient of the operator
 *   \f$({\cal O}_{HL}^{(3)})_{ii} =i\big(H^\dagger \overset{\leftrightarrow}{D^a_\mu} H\big)
 *   \big(\overline{L^i}\,\gamma^\mu \tau^a L^i\big)\f$ (flavor universal). </td>
 * </tr>
 * <tr>
 *   <td class="mod_name">%CHe </td>
 *   <td class="mod_symb">\f$ (C_{He})_{ii} \f$</td>
 *   <td class="mod_desc">The coefficient of the operator
 *   \f$({\cal O}_{He})_{ij} =i\big(H^\dagger \overset{\leftrightarrow}{D}_\mu H\big)
 *   \big(\overline{E^i}\,\gamma^\mu E^i\big)\f$ (flavor universal). </td>
 * </tr>
 * <tr>
 *   <td class="mod_name">%CHQ1 </td>
 *   <td class="mod_symb">\f$ (C_{HQ}^{(1)})_{ii} \f$</td>
 *   <td class="mod_desc">The coefficient of the operator
 *   \f$({\cal O}_{HQ}^{(1)})_{ii} =i\big(H^\dagger \overset{\leftrightarrow}{D}_\mu H\big)
 *   \big(\overline{Q^i}\,\gamma^\mu Q^i\big)\f$ (flavor universal). </td>
 * </tr>
 * <tr>
 *   <td class="mod_name">%CHQ3 </td>
 *   <td class="mod_symb">\f$ (C_{HQ}^{(3)})_{ii}\f$</td>
 *   <td class="mod_desc">The coefficient of the operator
 *   \f$({\cal O}_{HQ}^{(3)})_{ii} =i\big(H^\dagger \overset{\leftrightarrow}{D^a_\mu} H\big)
 *   \big(\overline{Q^i}\,\gamma^\mu \tau^a Q^i\big)\f$ (flavor universal). </td>
 * </tr>
 * <tr>
 *   <td class="mod_name">%CHu </td>
 *   <td class="mod_symb">\f$ (C_{Hu})_{ii} \f$</td>
 *   <td class="mod_desc">The coefficient of the operator
 *   \f$({\cal O}_{Hu})_{ii} =i\big(H^\dagger \overset{\leftrightarrow}{D}_\mu H\big)
 *   \big(\overline{U^i}\,\gamma^\mu U^i\big)\f$ (flavor universal). </td>
 * </tr>
 * <tr>
 *   <td class="mod_name">%CHd </td>
 *   <td class="mod_symb">\f$ (C_{Hd})_{ii} \f$</td>
 *   <td class="mod_desc">The coefficient of the operator
 *   \f$({\cal O}_{Hd})_{ii} =i\big(H^\dagger \overset{\leftrightarrow}{D}_\mu H\big)
 *   \big(\overline{D^i}\,\gamma^\mu D^i\big)\f$ (flavor universal). </td>
 * </tr>
  * <tr>
 *   <td class="mod_name">%CHud_r, CHud_i </td>
 *   <td class="mod_symb">\f$\mbox{Re}\big[(C_{Hud})_{ii}\big], \mbox{Im}\big[(C_{Hud})_{ii}\big] \f$</td>
 *   <td class="mod_desc">The real and imaginary parts of the coefficient of the operator
 *   \f$({\cal O}_{Hud})_{ii} =i\big(\widetilde{H}^\dagger D_\mu H\big)
 *   \big(\overline{U^i}\,\gamma^\mu D^i\big)\f$ (flavor universal). </td>
 * </tr>
 * <tr>
 *   <td class="mod_name">%CeH_jjr, CeH_jji </td>
 *   <td class="mod_symb">\f$\mbox{Re}\big[(C_{eH})_{jj}\big], \mbox{Im}\big[(C_{eH})_{jj}\big] \f$</td>
 *   <td class="mod_desc">The real and imaginary parts of the coefficient of the operator
 *   \f$({\cal O}_{eH})_{jj} =\big(H^\dagger H\big)
 *   \big(\overline{L^j}\,H E^j\big)\f$ (flavor universal). </td>
 * </tr>
 * <tr>
 *   <td class="mod_name">%CuH_jjr, CuH_jji </td>
 *   <td class="mod_symb">\f$\mbox{Re}\big[(C_{uH})_{jj}\big], \mbox{Im}\big[(C_{uH})_{jj}\big] \f$</td>
 *   <td class="mod_desc">The real and imaginary parts of the coefficient of the operator
 *   \f$({\cal O}_{uH})_{jj} =\big(H^\dagger H\big)
 *   \big(\overline{Q^j}\,\widetilde{H} U^j\big)\f$ (flavor universal). </td>
 * </tr>
 * <tr>
 *   <td class="mod_name">%CdH_jjr, CdH_jji </td>
 *   <td class="mod_symb">\f$\mbox{Re}\big[(C_{dH})_{jj}\big], \mbox{Im}\big[(C_{dH})_{jj}\big] \f$</td>
 *   <td class="mod_desc">The real and imaginary parts of the coefficient of the operator
 *   \f$({\cal O}_{dH})_{jj} =\big(H^\dagger H\big)
 *   \big(\overline{Q^j}\,H D^j\big)\f$ (flavor universal). </td>
 * </tr>
 * <tr>
 *   <td class="mod_name">%CuG_klr, CuG_kli </td>
 *   <td class="mod_symb">\f$\mbox{Re}\big[(C_{uG})_{kl}\big], \mbox{Im}\big[(C_{uG})_{kl}\big] \f$</td>
 *   <td class="mod_desc">The real and imaginary parts of the coefficient of the operator
 *   \f$({\cal O}_{uG})_{ij} =\big(\overline{Q^i}\sigma^{\mu\nu} T_A U^j\big)\widetilde{H} G_{\mu\nu}^A\f$, for \f$i,j=1,2,3\f$. </td>
 * </tr>
 * <tr>
 *   <td class="mod_name">%CuW_klr, CuW_kli </td>
 *   <td class="mod_symb">\f$\mbox{Re}\big[(C_{uW})_{kl}\big], \mbox{Im}\big[(C_{uW})_{kl}\big] \f$</td>
 *   <td class="mod_desc">The real and imaginary parts of the coefficient of the operator
 *   \f$({\cal O}_{uW})_{ij} =\big(\overline{Q^i}\sigma^{\mu\nu} \tau_a U^j\big)\widetilde{H} W_{\mu\nu}^a\f$, for \f$i,j=1,2,3\f$. </td>
 * </tr>
 * <tr>
 *   <td class="mod_name">%CuB_klr, CuB_kli </td>
 *   <td class="mod_symb">\f$\mbox{Re}\big[(C_{uB})_{kl}\big], \mbox{Im}\big[(C_{uB})_{kl}\big] \f$</td>
 *   <td class="mod_desc">The real and imaginary parts of the coefficient of the operator
 *   \f$({\cal O}_{uB})_{ij} =\big(\overline{Q^i}\sigma^{\mu\nu} U^j\big)\widetilde{H} B_{\mu\nu}\f$, for \f$i,j=1,2,3\f$. </td>
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
 * (The parameters associated to the theoretical uncertainties:
 * \f$\varepsilon_{X}^{int}\f$, \f$\varepsilon_{X}^{par}\f$ and
 * \f$\varepsilon_{X}^i(E)\f$,
 * are the same for both "NPSMEFTd6" and "NPSMEFTd6_LFU_QFU".)
 *
 * @anchor NPSMEFTd6Flags
 * <h3>%Model flags</h3>
 *
 * The Flags of NPSMEFTd6 are summarized below:
 * <table class="model">
 * <tr>
 *   <th>Label</th>
 *   <th>Value</th>
 *   <th>Description</th>
 * </tr>
 * <tr>
 *   <td class="mod_name">%QuadraticTerms</td>
 *   <td class="mod_valu">TRUE&nbsp;/&nbsp;<b>FALSE</b></td>
 *   <td class="mod_desc">This flag is set to TRUE if the quadratic terms in Higgs cross sections and widths are switched on.
 *   The default value is FALSE; new physics contributions are linearized.</td>
 * </tr>
 * <tr>
 *   <td class="mod_name">%RotateCHWCHB</td>
 *   <td class="mod_valu">TRUE&nbsp;/&nbsp;<b>FALSE</b></td>
 *   <td class="mod_desc">This flag is set to TRUE if using {sW2*CHW+cW2*CHB, -cW2*CHW+sW2*CHB}
 *   instead of {CHW, CHB} as floating parameters.
 *   The default value is FALSE.</td>
 * </tr>
 * </table>
 *
 *
 * @anchor NPSMEFTd6Functions
 * <h3>Important member functions</h3>
 *
 * See the base classes of the current class.
 *
 */
class NPSMEFTd6 : public NPbase {
public:

    /**
     *　@brief The number of the model parameters in %NPSMEFTd6. 
     */
    static const int NNPSMEFTd6Vars = 393;

    /**
     * @brief A string array containing the labels of the model parameters in
     * %NPSMEFTd6 if the model flag FlagRotateCHWCHB=false.
     */
    static const std::string NPSMEFTd6Vars[NNPSMEFTd6Vars];
      
    /**
     * @brief A string array containing the labels of the model parameters in
     * %NPSMEFTd6 if the model flag FlagRotateCHWCHB=true.
     */
    static const std::string NPSMEFTd6VarsRot[NNPSMEFTd6Vars];

    /**
     *　@brief The number of the model parameters in %NPSMEFTd6
     * with lepton and quark flavour universalities.
     */
    static const int NNPSMEFTd6Vars_LFU_QFU = 197;

    /**
     * @brief A string array containing the labels of the model parameters in
     * %NPSMEFTd6 with lepton and quark flavour universalities
     * if the model flag FlagRotateCHWCHB=false.
     */
    static const std::string NPSMEFTd6Vars_LFU_QFU[NNPSMEFTd6Vars_LFU_QFU];
      
    /**
     * @brief A string array containing the labels of the model parameters in
     * %NPSMEFTd6 with lepton and quark flavour universalities
     * if the model flag FlagRotateCHWCHB=true.
     */
    static const std::string NPSMEFTd6VarsRot_LFU_QFU[NNPSMEFTd6Vars_LFU_QFU];

    /**
     * @brief Constructor.
     * @param[in] FlagLeptonUniversal_in an internal boolean flag that is true if assuming lepton flavour universality
     * @param[in] FlagQuarkUniversal_in an internal boolean flag that is true if assuming quark flavour universality
     */
    NPSMEFTd6(const bool FlagLeptonUniversal_in = false, const bool FlagQuarkUniversal_in = false);

    /**
     * @brief The post-update method for %NPSMEFTd6.
     * @details This method runs all the procedures that are need to be executed
     * after the model is successfully updated.
     * @return a boolean that is true if the execution is successful
     */
    virtual bool PostUpdate();

    /**
     * @brief A method to check if all the mandatory parameters for %NPSMEFTd6
     * have been provided in model initialization.
     * @param[in] DPars a map of the parameters that are being updated in the Monte Carlo run
     * (including parameters that are varied and those that are held constant)
     * @return a boolean that is true if the execution is successful
     */
    virtual bool CheckParameters(const std::map<std::string, double>& DPars);

    /**
     * @brief A method to set a flag of %NPSMEFTd6.
     * @param[in] name name of a model flag
     * @param[in] value the boolean to be assigned to the flag specified by name
     * @return a boolean that is true if the execution is successful
     */
    virtual bool setFlag(const std::string name, const bool value);
    
    /**
     * @brief A method to get the Matching object for this model.
     * @return The matching object for this model
     */
    virtual NPSMEFTd6Matching& getMatching() const
    {
        return NPSMEFTd6M.getObj();
    }


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
    
    virtual double deltaMz2() const;
    
    /**
     * @brief The mass of the @f$W@f$ boson, @f$M_W@f$.
     * @return @f$M_W@f$ in GeV
     */
    virtual double Mw() const;

    /**
     * @brief A partial decay width of the @f$W@f$ boson decay into a SM fermion pair.
     * @param[in] fi a lepton or quark
     * @param[in] fj a lepton or quark
     * @return @f$\Gamma^W_{ij}@f$
     */
    virtual double GammaW(const Particle fi, const Particle fj) const;
    
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
 
    double CLL_mu() const;
    double CLL_tau() const;
    double CLL_up() const;
    double CLL_down() const;
    double CLL_charm() const;
    double CLL_strange() const;
    double CLL_bottom() const;
    
    double CLR_mu() const;
    double CLR_tau() const;
    double CLR_up() const;
    double CLR_down() const;
    double CLR_charm() const;
    double CLR_strange() const;
    double CLR_bottom() const;
    
    double CRL_mu() const;
    double CRL_tau() const;
    double CRL_up() const;
    double CRL_down() const;
    double CRL_charm() const;
    double CRL_strange() const;
    double CRL_bottom() const;
    
    double CRR_mu() const;
    double CRR_tau() const;
    double CRR_up() const;
    double CRR_down() const;
    double CRR_charm() const;
    double CRR_strange() const;
    double CRR_bottom() const;

    ////////////////////////////////////////////////////////////////////////
    
    ////////////////////////////////////////////////////////////////////////

    /**
     * @brief New physics contribution to the charged current coupling @f$W_\mu \bar{f_L}\gamma^mu f_L@f$.
     * @param[in] pbar a lepton or quark
     * @param[in] p a lepton or quark
     * @return @f$\delta g_{Wff}^{L}@f$
     */
    // no generation mixing
    virtual gslpp::complex deltaGL_Wff(const Particle pbar, const Particle p) const;
    /**
     * @brief New physics contribution to the charged current coupling @f$W_\mu \bar{f_R}\gamma^mu f_R@f$.
     * @param[in] pbar a lepton or quark
     * @param[in] p a lepton or quark
     * @return @f$\delta g_{Wff}^{R}@f$
     */
    // no generation mixing
    virtual gslpp::complex deltaGR_Wff(const Particle pbar, const Particle p) const;

    /**
     * @brief The new physics contribution to the coupling of the effective interaction @f$H G_{\mu\nu}^AG^{A \mu\nu}@f$.
     * @return @f$\delta g_{HGG}@f$
     */
    virtual double deltaG_hgg() const;
    /**
     * @brief The full new physics contribution to the coupling of the effective interaction @f$H G_{\mu\nu}^AG^{A \mu\nu}@f$,
     * including new local terms and modifications on the SM-loops. Normalized to the SM value.
     * @return @f$\delta g_{HGG}/g_{HGG}^SM}@f$
     */
    virtual double deltaG_hggRatio() const;
    /**
     * @brief The new physics contribution to the coupling of the effective interaction @f$H W_{\mu\nu}^\dagger W^{\mu\nu}@f$.
     * @return @f$\delta g_{HWW}^{(1)}@f$
     */
    virtual double deltaG1_hWW() const;
    /**
     * @brief The new physics contribution to the coupling of the effective interaction @f$H W_{\nu}^\dagger \partial^\mu W^{\mu\nu}@f$.
     * @return @f$\delta g_{HWW}^{(2)}@f$
     */
    virtual double deltaG2_hWW() const;
    /**
     * @brief The new physics contribution to the coupling of the effective interaction @f$H W_{\mu}^\dagger W^{\mu}@f$.
     * @return @f$\delta g_{HWW}^{(3)}@f$
     */
    virtual double deltaG3_hWW() const;
    /**
     * @brief The new physics contribution to the coupling of the effective interaction @f$H Z_{\mu\nu} Z^{\mu\nu}@f$.
     * @return @f$\delta g_{HZZ}^{(1)}@f$
     */
    virtual double deltaG1_hZZ() const;
    /**
     * @brief The new physics contribution to the coupling of the effective interaction @f$H Z_{\nu} \partial^\mu Z^{\mu\nu}@f$.
     * @return @f$\delta g_{HZZ}^{(2)}@f$
     */
    virtual double deltaG2_hZZ() const;
    /**
     * @brief The new physics contribution to the coupling of the effective interaction @f$H Z_{\mu} Z^{\mu}@f$.
     * @return @f$\delta g_{HZZ}^{(3)}@f$
     */
    virtual double deltaG3_hZZ() const;
    /**
     * @brief The new physics contribution to the coupling of the effective interaction @f$H Z_{\mu\nu} F^{\mu\nu}@f$.
     * @return @f$\delta g_{HZA}^{(1)}@f$
     */
    virtual double deltaG1_hZA() const;
    /**
     * @brief The full new physics contribution to the coupling of the effective interaction @f$H Z_{\mu\nu} F^{A \mu\nu}@f$,
     * including new local terms and modifications on the SM-loops. Normalized to the SM value.
     * @return @f$\delta g_{HZA}^{(1)}/g_{HZA}^{(1),SM}@f$
     */
    virtual double deltaG1_hZARatio() const;
    /**
     * @brief The new physics contribution to the coupling of the effective interaction @f$H Z_{\nu} \partial^\mu F^{\mu\nu}@f$.
     * @return @f$\delta g_{HZA}^{(2)}@f$
     */
    virtual double deltaG2_hZA() const;
    /**
     * @brief The new physics contribution to the coupling of the effective interaction @f$H F_{\mu\nu} F^{\mu\nu}@f$.
     * @return @f$\delta g_{HAA}@f$
     */
    virtual double deltaG_hAA() const;
    /**
     * @brief The full new physics contribution to the coupling of the effective interaction @f$H F_{\mu\nu} F^{\mu\nu}@f$,
     * including new local terms and modifications on the SM-loops. Normalized to the SM value.
     * @return @f$\delta g_{HAA}/g_{HAA}^SM}@f$
     */
    virtual double deltaG_hAARatio() const;
    /**
     * @brief The new physics contribution to the coupling of the effective interaction @f$H f\bar{f}@f$.
     * @param[in] p a lepton or quark
     * @return @f$\delta g_{Hff}@f$
     */
    // no generation mixing
    virtual gslpp::complex deltaG_hff(const Particle p) const;

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
      
    /**
     * @brief The new physics contribution to the coupling of the effective interaction @f$H G_{\mu\nu} \bar{f}\sigmma^{\mu\nu} f@f$.
     * @param[in] p a lepton or quark
     * @return @f$\delta g_{hGff}@f$
     */
    // no generation mixing
    gslpp::complex deltaG_hGff(const Particle p) const;
      
    /**
     * @brief The new physics contribution to the coupling of the effective interaction @f$H Z_{\mu\nu} \bar{f}\sigmma^{\mu\nu} f@f$.
     * @param[in] p a lepton or quark
     * @return @f$\delta g_{hZff}@f$
     */
    // no generation mixing
    gslpp::complex deltaG_hZff(const Particle p) const;
      
    /**
     * @brief The new physics contribution to the coupling of the effective interaction @f$H A_{\mu\nu} \bar{f}\sigmma^{\mu\nu} f@f$.
     * @param[in] p a lepton or quark
     * @return @f$\delta g_{hAff}@f$
     */
    // no generation mixing
    gslpp::complex deltaG_hAff(const Particle p) const;
      
    /**
     * @brief The new physics contribution to the coupling of the effective interaction @f$G_{\mu\nu} \bar{f}\sigmma^{\mu\nu} f@f$.
     * @param[in] p a lepton or quark
     * @return @f$\delta g_{Gff}@f$
     */
    // no generation mixing
    gslpp::complex deltaG_Gff(const Particle p) const;
      
    /**
     * @brief The new physics contribution to the coupling of the effective interaction @f$Z_{\mu\nu} \bar{f}\sigmma^{\mu\nu} f@f$.
     * @param[in] p a lepton or quark
     * @return @f$\delta g_{Zff}@f$
     */
    // no generation mixing
    gslpp::complex deltaG_Zff(const Particle p) const;
      
    /**
     * @brief The new physics contribution to the coupling of the effective interaction @f$A_{\mu\nu} \bar{f}\sigmma^{\mu\nu} f@f$.
     * @param[in] p a lepton or quark
     * @return @f$\delta g_{Aff}@f$
     */
    // no generation mixing
    gslpp::complex deltaG_Aff(const Particle p) const;
      
    /**
     * @brief The new physics contribution to the coupling of the effective interaction @f$f_{ABC} G_{\mu\nu}^A G_{\nu\rho}^B G_{\rho\mu}^C@f$.
     * @return @f$\delta g_{3G}@f$
     */
    double deltag3G() const;


    ////////////////////////////////////////////////////////////////////////

    /**
     * @brief Loop function entering in the calculation of the effective @f$Hgg@f$ and @f$H\gamma\gamma@f$ couplings.
     * @param[in] @f$\tau=4 M^2/m_h^2@f$, with @f$M@f$ the mass of the particle in the loop.
     * @return @f$f(\tau)@f$
     */
    gslpp::complex f_triangle(const double tau) const;
    /**
     * @brief Loop function entering in the calculation of the effective @f$HZ\gamma@f$ coupling.
     * @param[in] @f$\tau=4 M^2/m_h^2@f$, with @f$M@f$ the mass of the particle in the loop.
     * @return @f$g(\tau)@f$
     */
    gslpp::complex g_triangle(const double tau) const;
    /**
     * @brief Loop function entering in the calculation of the effective @f$HZ\gamma@f$ coupling.
     * @param[in] @f$\tau=4 M^2/m_h^2@f$, @f$\lambda=4 M^2/m_Z^2@f$, with @f$M@f$ the mass of the particle in the loop.
     * @return @f$I_1(\tau,\lambda)@f$
     */
    gslpp::complex I_triangle_1(const double tau, const double lambda) const;
    /**
     * @brief Loop function entering in the calculation of the effective @f$HZ\gamma@f$ coupling.
     * @param[in] @f$\tau=4 M^2/m_h^2@f$, @f$\lambda=4 M^2/m_Z^2@f$, with @f$M@f$ the mass of the particle in the loop.
     * @return @f$I_2(\tau,\lambda)@f$
     */
    gslpp::complex I_triangle_2(const double tau, const double lambda) const;
    /**
     * @brief Fermionic loop function entering in the calculation of the effective @f$Hgg@f$ and @f$H\gamma\gamma@f$ couplings.
     * @details
     * @f$A^H_f(\tau)=2\tau [1+(1-\tau)f(\tau)]@f$
     * @param[in] @f$\tau=4 M^2/m_h^2@f$, with @f$M@f$ the mass of the fermion in the loop.
     * @return @f$A^H_f(\tau)@f$
     */
    gslpp::complex AH_f(const double tau) const;
    
    /**
     * @brief W loop function entering in the calculation of the effective @f$H\gamma\gamma@f$ coupling.
     * @details
     * @f$A^H_W(\tau)=-[2+3\tau + 3\tau*(2-\tau) f(\tau)]@f$
     * @param[in] @f$\tau=4 M^2/m_h^2@f$, with @f$M@f$ the mass of the fermion in the loop.
     * @return @f$A^H_W(\tau)@f$
     */
    gslpp::complex AH_W(const double tau) const;
    
    /**
     * @brief Fermionic loop function entering in the calculation of the effective @f$HZ\gamma@f$ coupling.
     * @details
     * @param[in] @f$\tau=4 M^2/m_h^2@f$, @f$\lambda=4 M^2/m_Z^2@f$, with @f$M@f$ the mass of the fermion in the loop.
     * @return @f$A^{HZ\gamma}_f(\tau,\lambda)@f$
     */
    gslpp::complex AHZga_f(const double tau, const double lambda) const;
    
    /**
     * @brief W loop function entering in the calculation of the effective @f$HZ\gamma@f$ coupling.
     * @details
     * @param[in] @f$\tau=4 M^2/m_h^2@f$, @f$\lambda=4 M^2/m_Z^2@f$, with @f$M@f$ the mass of the fermion in the loop.
     * @return @f$A^{HZ\gamma}_W(\tau,\lambda)@f$
     */
    gslpp::complex AHZga_W(const double tau, const double lambda) const;

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
     * @brief The ratio @f$\mu_{eeWBF}@f$ between the 
     * @f$ e^{+}e^{-}\to \nu\bar{\nu} H @f$ production
     * cross-section in the current model and in the Standard Model.
     * @param[in] sqrt_s the center-of-mass energy in TeV
     * @return @f$\mu_{eeWBF}@f$
     */
    virtual double mueeWBF(const double sqrt_s) const;
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
     * @brief The ratio @f$\mu_{eettH}@f$ between the 
     * @f$ e^{+}e^{-}\to t\bar{t} H @f$ production
     * cross-section in the current model and in the Standard Model.
     * @param[in] sqrt_s the center-of-mass energy in TeV
     * @return @f$\mu_{eettH}@f$
     */
    virtual double mueettH(const double sqrt_s) const;
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
     * @brief The ratio of the Br@f$(H\to \mu^+\mu^-)@f$ in the current model
     * and in the Standard Model.
     * @return Br@f$(H\to \mu^+\mu^-)@f$/Br@f$(H\to \mu^+\mu^-)_{\mathrm{SM}}@f$
     */
    virtual double BrHmumuRatio() const;
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
     * @brief The ratio of the @f$\Gamma(H\to \mu\mu)@f$ in the current model
     * and in the Standard Model.
     * @return @f$\Gamma(H\to \mu\mu)@f$/@f$\Gamma(H\to \mu\mu)_{\mathrm{SM}}@f$
     */
    double GammaHmumuRatio() const;
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
     * @brief The new physics contribution to the ratio of the @f$\Gamma(H\to \mu\mu)@f$
     * in the current model and in the Standard Model. (Only terms that are
     * linear in the effective Lagrangian coefficients.)
     * @return @f$\delta \Gamma(H\to \mu\mu)@f$/@f$\Gamma(H\to \mu\mu)_{\mathrm{SM}}@f$
     */
    double deltaGammaHmumuRatio1() const;
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
     * @brief The new physics contribution to the ratio of the @f$\Gamma(H\to \mu\mu)@f$
     * in the current model and in the Standard Model. (Only terms that are
     * quadratic in the effective Lagrangian coefficients.)
     * @return @f$\delta \Gamma(H\to \mu\mu)@f$/@f$\Gamma(H\to \mu\mu)_{\mathrm{SM}}@f$
     */
    double deltaGammaHmumuRatio2() const;
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
      
    /**
     * @brief The new physics contribution to the anomalous triple gauge coupling @f$g_{1,Z}@f$.
     * @return @f$\delta g_{1,Z}@f$
     */
    virtual double deltag1ZNP() const;
      
    /**
     * @brief The new physics contribution to the anomalous triple gauge coupling @f$\kappa_{\gamma}@f$.
     * @return @f$\delta \kappa_{\gamma}@f$
     */
    virtual double deltaKgammaNP() const;
      
    /**
     * @brief The new physics contribution to the anomalous triple gauge coupling @f$\lambda_{Z}@f$.
     * @return @f$\lambda_{Z}@f$
     */
    virtual double lambdaZNP() const;
    
    ////////////////////////////////////////////////////////////////////////
protected:

    /**
     * @brief @copybrief Model::setParameter()
     * @copydetails Model::setParameter()
     */
    virtual void setParameter(const std::string name, const double& value);
    mutable Matching<NPSMEFTd6Matching,NPSMEFTd6> NPSMEFTd6M;

    double CG; ///< The dimension-6 operator coefficient \f$C_{G}\f$.
    double CW; ///< The dimension-6 operator coefficient \f$C_{W}\f$.
    double CHG; ///< The dimension-6 operator coefficient \f$C_{HG}\f$.
    double CHW; ///< The dimension-6 operator coefficient \f$C_{HW}\f$.
    double CHB; ///< The dimension-6 operator coefficient \f$C_{HB}\f$.
    double CHWHB_gaga; ///< The combination of dimension-6 operator coefficients entering in \f$\delta_{AA}\f$: \f$s_W^2 C_{HW} + c_W^2 C_{HW}\f$.
    double CHWHB_gagaorth; ///< The combination of dimension-6 operator coefficients \f$-c_W^2 C_{HW} + s_W^2 C_{HW}\f$.
    double CDHB; ///< The dimension-6 operator coefficient \f$C_{DHB}\f$.
    double CDHW; ///< The dimension-6 operator coefficient \f$C_{DHW}\f$.
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
    double CHe_11; ///< The dimension-6 operator coefficient \f$(C_{He})_{11}\f$.
    double CHe_12r; ///< The dimension-6 operator coefficient \f$(C_{He})_{12}\f$ (real part).
    double CHe_13r; ///< The dimension-6 operator coefficient \f$(C_{He})_{13}\f$ (real part).
    double CHe_22; ///< The dimension-6 operator coefficient \f$(C_{He})_{22}\f$.
    double CHe_23r; ///< The dimension-6 operator coefficient \f$(C_{He})_{23}\f$ (real part).
    double CHe_33; ///< The dimension-6 operator coefficient \f$(C_{He})_{33}\f$.
    double CHe_12i; ///< The dimension-6 operator coefficient \f$(C_{He})_{12}\f$ (imaginary part).
    double CHe_13i; ///< The dimension-6 operator coefficient \f$(C_{He})_{13}\f$ (imaginary part).
    double CHe_23i; ///< The dimension-6 operator coefficient \f$(C_{He})_{23}\f$ (imaginary part).
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
    double CHu_11; ///< The dimension-6 operator coefficient \f$(C_{Hu})_{11}\f$.
    double CHu_12r; ///< The dimension-6 operator coefficient \f$(C_{Hu})_{12}\f$ (real part).
    double CHu_13r; ///< The dimension-6 operator coefficient \f$(C_{Hu})_{13}\f$ (real part).
    double CHu_22; ///< The dimension-6 operator coefficient \f$(C_{Hu})_{22}\f$.
    double CHu_23r; ///< The dimension-6 operator coefficient \f$(C_{Hu})_{23}\f$ (real part).
    double CHu_33; ///< The dimension-6 operator coefficient \f$(C_{Hu})_{33}\f$.
    double CHu_12i; ///< The dimension-6 operator coefficient \f$(C_{Hu})_{12}\f$ (imaginary part).
    double CHu_13i; ///< The dimension-6 operator coefficient \f$(C_{Hu})_{13}\f$ (imaginary part).
    double CHu_23i; ///< The dimension-6 operator coefficient \f$(C_{Hu})_{23}\f$ (imaginary part).
    double CHd_11; ///< The dimension-6 operator coefficient \f$(C_{Hd})_{11}\f$.
    double CHd_12r; ///< The dimension-6 operator coefficient \f$(C_{Hd})_{12}\f$ (real part).
    double CHd_13r; ///< The dimension-6 operator coefficient \f$(C_{Hd})_{13}\f$ (real part).
    double CHd_22; ///< The dimension-6 operator coefficient \f$(C_{Hd})_{22}\f$.
    double CHd_23r; ///< The dimension-6 operator coefficient \f$(C_{Hd})_{23}\f$ (real part).
    double CHd_33; ///< The dimension-6 operator coefficient \f$(C_{Hd})_{33}\f$.
    double CHd_12i; ///< The dimension-6 operator coefficient \f$(C_{Hd})_{12}\f$ (imaginary part).
    double CHd_13i; ///< The dimension-6 operator coefficient \f$(C_{Hd})_{13}\f$ (imaginary part).
    double CHd_23i; ///< The dimension-6 operator coefficient \f$(C_{Hd})_{23}\f$ (imaginary part).
    double CHud_11r; ///< The dimension-6 operator coefficient \f$(C_{Hud})_{11}\f$ (real part).
    double CHud_12r; ///< The dimension-6 operator coefficient \f$(C_{Hud})_{12}\f$ (real part).
    double CHud_13r; ///< The dimension-6 operator coefficient \f$(C_{Hud})_{13}\f$ (real part).
    double CHud_22r; ///< The dimension-6 operator coefficient \f$(C_{Hud})_{22}\f$ (real part).
    double CHud_23r; ///< The dimension-6 operator coefficient \f$(C_{Hud})_{23}\f$ (real part).
    double CHud_33r; ///< The dimension-6 operator coefficient \f$(C_{Hud})_{33}\f$ (real part).
    double CHud_11i; ///< The dimension-6 operator coefficient \f$(C_{Hud})_{11}\f$ (imaginary part).
    double CHud_12i; ///< The dimension-6 operator coefficient \f$(C_{Hud})_{12}\f$ (imaginary part).
    double CHud_13i; ///< The dimension-6 operator coefficient \f$(C_{Hud})_{13}\f$ (imaginary part).
    double CHud_22i; ///< The dimension-6 operator coefficient \f$(C_{Hud})_{22}\f$ (imaginary part).
    double CHud_23i; ///< The dimension-6 operator coefficient \f$(C_{Hud})_{23}\f$ (imaginary part).
    double CHud_33i; ///< The dimension-6 operator coefficient \f$(C_{Hud})_{33}\f$ (imaginary part).
    double CeH_11r; ///< The dimension-6 operator coefficient \f$(C_{eH})_{11}\f$ (real part).
    double CeH_12r; ///< The dimension-6 operator coefficient \f$(C_{eH})_{12}\f$ (real part).
    double CeH_13r; ///< The dimension-6 operator coefficient \f$(C_{eH})_{13}\f$ (real part).
    double CeH_22r; ///< The dimension-6 operator coefficient \f$(C_{eH})_{22}\f$ (real part).
    double CeH_23r; ///< The dimension-6 operator coefficient \f$(C_{eH})_{23}\f$ (real part).
    double CeH_33r; ///< The dimension-6 operator coefficient \f$(C_{eH})_{33}\f$ (real part).
    double CeH_11i; ///< The dimension-6 operator coefficient \f$(C_{eH})_{11}\f$ (imaginary part).
    double CeH_12i; ///< The dimension-6 operator coefficient \f$(C_{eH})_{12}\f$ (imaginary part).
    double CeH_13i; ///< The dimension-6 operator coefficient \f$(C_{eH})_{13}\f$ (imaginary part).
    double CeH_22i; ///< The dimension-6 operator coefficient \f$(C_{eH})_{22}\f$ (imaginary part).
    double CeH_23i; ///< The dimension-6 operator coefficient \f$(C_{eH})_{23}\f$ (imaginary part).
    double CeH_33i; ///< The dimension-6 operator coefficient \f$(C_{eH})_{33}\f$ (imaginary part).
    double CuH_11r; ///< The dimension-6 operator coefficient \f$(C_{uH})_{11}\f$ (real part).
    double CuH_12r; ///< The dimension-6 operator coefficient \f$(C_{uH})_{12}\f$ (real part).
    double CuH_13r; ///< The dimension-6 operator coefficient \f$(C_{uH})_{13}\f$ (real part).
    double CuH_22r; ///< The dimension-6 operator coefficient \f$(C_{uH})_{22}\f$ (real part).
    double CuH_23r; ///< The dimension-6 operator coefficient \f$(C_{uH})_{23}\f$ (real part).
    double CuH_33r; ///< The dimension-6 operator coefficient \f$(C_{uH})_{33}\f$ (real part).
    double CuH_11i; ///< The dimension-6 operator coefficient \f$(C_{uH})_{11}\f$ (imaginary part).
    double CuH_12i; ///< The dimension-6 operator coefficient \f$(C_{uH})_{12}\f$ (imaginary part).
    double CuH_13i; ///< The dimension-6 operator coefficient \f$(C_{uH})_{13}\f$ (imaginary part).
    double CuH_22i; ///< The dimension-6 operator coefficient \f$(C_{uH})_{22}\f$ (imaginary part).
    double CuH_23i; ///< The dimension-6 operator coefficient \f$(C_{uH})_{23}\f$ (imaginary part).
    double CuH_33i; ///< The dimension-6 operator coefficient \f$(C_{uH})_{33}\f$ (imaginary part).
    double CdH_11r; ///< The dimension-6 operator coefficient \f$(C_{dH})_{11}\f$ (real part).
    double CdH_12r; ///< The dimension-6 operator coefficient \f$(C_{dH})_{12}\f$ (real part).
    double CdH_13r; ///< The dimension-6 operator coefficient \f$(C_{dH})_{13}\f$ (real part).
    double CdH_22r; ///< The dimension-6 operator coefficient \f$(C_{dH})_{22}\f$ (real part).
    double CdH_23r; ///< The dimension-6 operator coefficient \f$(C_{dH})_{23}\f$ (real part).
    double CdH_33r; ///< The dimension-6 operator coefficient \f$(C_{dH})_{33}\f$ (real part).
    double CdH_11i; ///< The dimension-6 operator coefficient \f$(C_{dH})_{11}\f$ (imaginary part).
    double CdH_12i; ///< The dimension-6 operator coefficient \f$(C_{dH})_{12}\f$ (imaginary part).
    double CdH_13i; ///< The dimension-6 operator coefficient \f$(C_{dH})_{13}\f$ (imaginary part).
    double CdH_22i; ///< The dimension-6 operator coefficient \f$(C_{dH})_{22}\f$ (imaginary part).
    double CdH_23i; ///< The dimension-6 operator coefficient \f$(C_{dH})_{23}\f$ (imaginary part).
    double CdH_33i; ///< The dimension-6 operator coefficient \f$(C_{dH})_{33}\f$ (imaginary part).
    double CuG_11r; ///< The dimension-6 operator coefficient \f$(C_{uG})_{11}\f$ (real part).
    double CuG_12r; ///< The dimension-6 operator coefficient \f$(C_{uG})_{12}\f$ (real part).
    double CuG_13r; ///< The dimension-6 operator coefficient \f$(C_{uG})_{13}\f$ (real part).
    double CuG_22r; ///< The dimension-6 operator coefficient \f$(C_{uG})_{22}\f$ (real part).
    double CuG_23r; ///< The dimension-6 operator coefficient \f$(C_{uG})_{23}\f$ (real part).
    double CuG_33r; ///< The dimension-6 operator coefficient \f$(C_{uG})_{33}\f$ (real part).
    double CuG_11i; ///< The dimension-6 operator coefficient \f$(C_{uG})_{11}\f$ (imaginary part).
    double CuG_12i; ///< The dimension-6 operator coefficient \f$(C_{uG})_{12}\f$ (imaginary part).
    double CuG_13i; ///< The dimension-6 operator coefficient \f$(C_{uG})_{13}\f$ (imaginary part).
    double CuG_22i; ///< The dimension-6 operator coefficient \f$(C_{uG})_{22}\f$ (imaginary part).
    double CuG_23i; ///< The dimension-6 operator coefficient \f$(C_{uG})_{23}\f$ (imaginary part).
    double CuG_33i; ///< The dimension-6 operator coefficient \f$(C_{uG})_{33}\f$ (imaginary part).
    double CuW_11r; ///< The dimension-6 operator coefficient \f$(C_{uW})_{11}\f$ (real part).
    double CuW_12r; ///< The dimension-6 operator coefficient \f$(C_{uW})_{12}\f$ (real part).
    double CuW_13r; ///< The dimension-6 operator coefficient \f$(C_{uW})_{13}\f$ (real part).
    double CuW_22r; ///< The dimension-6 operator coefficient \f$(C_{uW})_{22}\f$ (real part).
    double CuW_23r; ///< The dimension-6 operator coefficient \f$(C_{uW})_{23}\f$ (real part).
    double CuW_33r; ///< The dimension-6 operator coefficient \f$(C_{uW})_{33}\f$ (real part).
    double CuW_11i; ///< The dimension-6 operator coefficient \f$(C_{uW})_{11}\f$ (imaginary part).
    double CuW_12i; ///< The dimension-6 operator coefficient \f$(C_{uW})_{12}\f$ (imaginary part).
    double CuW_13i; ///< The dimension-6 operator coefficient \f$(C_{uW})_{13}\f$ (imaginary part).
    double CuW_22i; ///< The dimension-6 operator coefficient \f$(C_{uW})_{22}\f$ (imaginary part).
    double CuW_23i; ///< The dimension-6 operator coefficient \f$(C_{uW})_{23}\f$ (imaginary part).
    double CuW_33i; ///< The dimension-6 operator coefficient \f$(C_{uW})_{33}\f$ (imaginary part).
    double CuB_11r; ///< The dimension-6 operator coefficient \f$(C_{uB})_{11}\f$ (real part).
    double CuB_12r; ///< The dimension-6 operator coefficient \f$(C_{uB})_{12}\f$ (real part).
    double CuB_13r; ///< The dimension-6 operator coefficient \f$(C_{uB})_{13}\f$ (real part).
    double CuB_22r; ///< The dimension-6 operator coefficient \f$(C_{uB})_{22}\f$ (real part).
    double CuB_23r; ///< The dimension-6 operator coefficient \f$(C_{uB})_{23}\f$ (real part).
    double CuB_33r; ///< The dimension-6 operator coefficient \f$(C_{uB})_{33}\f$ (real part).
    double CuB_11i; ///< The dimension-6 operator coefficient \f$(C_{uB})_{11}\f$ (imaginary part).
    double CuB_12i; ///< The dimension-6 operator coefficient \f$(C_{uB})_{12}\f$ (imaginary part).
    double CuB_13i; ///< The dimension-6 operator coefficient \f$(C_{uB})_{13}\f$ (imaginary part).
    double CuB_22i; ///< The dimension-6 operator coefficient \f$(C_{uB})_{22}\f$ (imaginary part).
    double CuB_23i; ///< The dimension-6 operator coefficient \f$(C_{uB})_{23}\f$ (imaginary part).
    double CuB_33i; ///< The dimension-6 operator coefficient \f$(C_{uB})_{33}\f$ (imaginary part).
    double CLL_1111;
    double CLL_1221, CLL_2112, CLL_1122, CLL_2211;
    double CLL_1331, CLL_3113, CLL_1133, CLL_3311;
    double CLQ1_1111;
    double CLQ1_1122, CLQ1_2211, CLQ1_1221, CLQ1_2112;
    double CLQ1_1133, CLQ1_3311, CLQ1_1331, CLQ1_3113;
    double CLQ1_1123, CLQ1_2223, CLQ1_3323;
    double CLQ1_1132, CLQ1_2232, CLQ1_3332;
    double CLQ3_1111;
    double CLQ3_1122, CLQ3_2211, CLQ3_1221, CLQ3_2112;
    double CLQ3_1133, CLQ3_3311, CLQ3_1331, CLQ3_3113;
    double CLQ3_1123, CLQ3_2223, CLQ3_3323;
    double CLQ3_1132, CLQ3_2232, CLQ3_3332;
    double Cee_1111;
    double Cee_1122, Cee_2211;
    double Cee_1133, Cee_3311;
    double Ceu_1111;
    double Ceu_1122, Ceu_2211;
    double Ceu_1133, Ceu_3311;
    double Ced_1111;
    double Ced_1122, Ced_2211;
    double Ced_1133, Ced_3311;
    double Ced_1123, Ced_2223, Ced_3323;
    double Ced_1132, Ced_2232, Ced_3332;
    double CLe_1111;
    double CLe_1122, CLe_2211;
    double CLe_1133, CLe_3311;
    double CLu_1111;
    double CLu_1122, CLu_2211;
    double CLu_1133, CLu_3311;
    double CLd_1111;
    double CLd_1122, CLd_2211;
    double CLd_1133, CLd_3311;
    double CLd_1123, CLd_2223, CLd_3323;
    double CLd_1132, CLd_2232, CLd_3332;
    double CQe_1111;
    double CQe_1122, CQe_2211;
    double CQe_1133, CQe_3311;
    double CQe_2311, CQe_2322, CQe_2333;
    double CQe_3211, CQe_3222, CQe_3233;
    double Lambda_NP; ///< The new physics scale [GeV].
// The intrinsic and parametric theory relative errors in the Higgs observables.
    double eggFint; ///< Intrinsic relative theoretical error in ggF production. (Assumed to be constant in energy.)
    double eggFpar; ///< Parametric relative theoretical error in ggF production. (Assumed to be constant in energy.)
    double ettHint; ///< Intrinsic relative theoretical error in ttH production. (Assumed to be constant in energy.)
    double ettHpar; ///< Parametric relative theoretical error in ttH production. (Assumed to be constant in energy.)
    double eVBFint; ///< Intrinsic relative theoretical error in VBF production. (Assumed to be constant in energy.)
    double eVBFpar; ///< Parametric relative theoretical error in VBF production. (Assumed to be constant in energy.)
    double eWHint; ///< Intrinsic relative theoretical error in WH production. (Assumed to be constant in energy.)
    double eWHpar; ///< Parametric relative theoretical error in WH production. (Assumed to be constant in energy.)
    double eZHint; ///< Intrinsic relative theoretical error in ZH production. (Assumed to be constant in energy.)
    double eZHpar; ///< Parametric relative theoretical error in ZH production. (Assumed to be constant in energy.)
    double eeeWBFint; ///< Intrinsic relative theoretical error in \f$e^+ e^- \to H \bar{\nu} \nu\f$. (Assumed to be constant in energy.)
    double eeeWBFpar; ///< Parametric relative theoretical error in \f$e^+ e^- \to H \bar{\nu} \nu\f$. (Assumed to be constant in energy.)
    double eeeZHint; ///< Intrinsic relative theoretical error in \f$e^+ e^- \to Z H\f$. (Assumed to be constant in energy.)
    double eeeZHpar; ///< Parametric relative theoretical error in \f$e^+ e^- \to Z H\f$. (Assumed to be constant in energy.)
    double eeettHint; ///< Intrinsic relative theoretical error in \f$e^+ e^- \to t \bar{t} H\f$. (Assumed to be constant in energy.)
    double eeettHpar; ///< Parametric relative theoretical error in \f$e^+ e^- \to t \bar{t} H\f$. (Assumed to be constant in energy.)
    double eHggint; ///< Intrinsic relative theoretical error in \f$H \to g g\f$.
    double eHggpar; ///< Parametric relative theoretical error in \f$H \to g g\f$.
    double eHWWint; ///< Intrinsic relative theoretical error in \f$H \to W W\f$.
    double eHWWpar; ///< Parametric relative theoretical error in \f$H \to W W\f$.
    double eHZZint; ///< Intrinsic relative theoretical error in \f$H \to Z Z\f$.
    double eHZZpar; ///< Parametric relative theoretical error in \f$H \to Z Z\f$.
    double eHZgaint; ///< Intrinsic relative theoretical error in \f$H \to Z \gamma\f$.
    double eHZgapar; ///< Parametric relative theoretical error in \f$H \to Z \gamma\f$.
    double eHgagaint; ///< Intrinsic relative theoretical error in \f$H \to \gamma\gamma\f$.
    double eHgagapar; ///< Parametric relative theoretical error in \f$H \to \gamma\gamma\f$.
    double eHmumuint; ///< Intrinsic relative theoretical error in \f$H \to \mu^+ \mu^-\f$.
    double eHmumupar; ///< Parametric relative theoretical error in \f$H \to \mu^+ \mu^-\f$.
    double eHtautauint; ///< Intrinsic relative theoretical error in \f$H \to \tau^+ \tau^-\f$.
    double eHtautaupar; ///< Parametric relative theoretical error in \f$H \to \tau^+ \tau^-\f$.
    double eHccint; ///< Intrinsic relative theoretical error in \f$H \to c\bar{c}\f$.
    double eHccpar; ///< Parametric relative theoretical error in \f$H \to c\bar{c}\f$.
    double eHbbint; ///< Intrinsic relative theoretical error in \f$H \to b\bar{b}\f$.
    double eHbbpar; ///< Parametric relative theoretical error in \f$H \to b\bar{b}\f$.
// The error in the parameters multiplying the dimension-6 operator coefficients in the production cross sections.
    double eVBF_2_Hbox;///< Theoretical uncertainty in the (linear) new physics contribution from \f$C_{H\Box}\f$ to VBF production at Tevatron (1.96 TeV).
    double eVBF_2_HQ1_11;///< Theoretical uncertainty in the (linear) new physics contribution from \f$(C_{HQ}^{(1)})_{11}\f$ to VBF production at Tevatron (1.96 TeV).
    double eVBF_2_Hu_11;///< Theoretical uncertainty in the (linear) new physics contribution from \f$(C_{Hu})_{11}\f$ to VBF production at Tevatron (1.96 TeV).
    double eVBF_2_Hd_11;///< Theoretical uncertainty in the (linear) new physics contribution from \f$(C_{Hd})_{11}\f$ to VBF production at Tevatron (1.96 TeV).
    double eVBF_2_HQ3_11;///< Theoretical uncertainty in the (linear) new physics contribution from \f$(C_{HQ}^{(3)})_{11}\f$ to VBF production at Tevatron (1.96 TeV).
    double eVBF_2_HD;///< Theoretical uncertainty in the (linear) new physics contribution from \f$C_{HD}\f$ to VBF production at Tevatron (1.96 TeV).
    double eVBF_2_HB;///< Theoretical uncertainty in the (linear) new physics contribution from \f$C_{HB}\f$ to VBF production at Tevatron (1.96 TeV).
    double eVBF_2_HW;///< Theoretical uncertainty in the (linear) new physics contribution from \f$C_{HW}\f$ to VBF production at Tevatron (1.96 TeV).
    double eVBF_2_HWB;///< Theoretical uncertainty in the (linear) new physics contribution from \f$C_{HWB}\f$ to VBF production at Tevatron (1.96 TeV).
    double eVBF_2_HG;///< Theoretical uncertainty in the (linear) new physics contribution from \f$C_{HG}\f$ to VBF production at Tevatron (1.96 TeV).
    double eVBF_2_DHB;///< Theoretical uncertainty in the (linear) new physics contribution from \f$C_{DHB}\f$ to VBF production at Tevatron (1.96 TeV).
    double eVBF_2_DHW;///< Theoretical uncertainty in the (linear) new physics contribution from \f$C_{DHW}\f$ to VBF production at Tevatron (1.96 TeV).
    double eVBF_2_DeltaGF;///< Theoretical uncertainty in the (linear) new physics contribution from \f$\delta_{G_F}\f$ to VBF production at Tevatron (1.96 TeV).
      
    double eVBF_78_Hbox;///< Theoretical uncertainty in the (linear) new physics contribution from \f$C_{H\Box}\f$ to VBF production at Tevatron (7 & 8 TeV).
    double eVBF_78_HQ1_11;///< Theoretical uncertainty in the (linear) new physics contribution from \f$(C_{HQ}^{(1)})_{11}\f$ to VBF production at Tevatron (7 & 8 TeV).
    double eVBF_78_Hu_11;///< Theoretical uncertainty in the (linear) new physics contribution from \f$(C_{Hu})_{11}\f$ to VBF production at Tevatron (7 & 8 TeV).
    double eVBF_78_Hd_11;///< Theoretical uncertainty in the (linear) new physics contribution from \f$(C_{Hd})_{11}\f$ to VBF production at Tevatron (7 & 8 TeV).
    double eVBF_78_HQ3_11;///< Theoretical uncertainty in the (linear) new physics contribution from \f$(C_{HQ}^{(3)})_{11}\f$ to VBF production at Tevatron (7 & 8 TeV).
    double eVBF_78_HD;///< Theoretical uncertainty in the (linear) new physics contribution from \f$C_{HD}\f$ to VBF production at Tevatron (7 & 8 TeV).
    double eVBF_78_HB;///< Theoretical uncertainty in the (linear) new physics contribution from \f$C_{HB}\f$ to VBF production at Tevatron (7 & 8 TeV).
    double eVBF_78_HW;///< Theoretical uncertainty in the (linear) new physics contribution from \f$C_{HW}\f$ to VBF production at Tevatron (7 & 8 TeV).
    double eVBF_78_HWB;///< Theoretical uncertainty in the (linear) new physics contribution from \f$C_{HWB}\f$ to VBF production at Tevatron (7 & 8 TeV).
    double eVBF_78_HG;///< Theoretical uncertainty in the (linear) new physics contribution from \f$C_{HG}\f$ to VBF production at Tevatron (7 & 8 TeV).
    double eVBF_78_DHB;///< Theoretical uncertainty in the (linear) new physics contribution from \f$C_{DHB}\f$ to VBF production at Tevatron (7 & 8 TeV).
    double eVBF_78_DHW;///< Theoretical uncertainty in the (linear) new physics contribution from \f$C_{DHW}\f$ to VBF production at Tevatron (7 & 8 TeV).
    double eVBF_78_DeltaGF;///< Theoretical uncertainty in the (linear) new physics contribution from \f$\delta_{G_F}\f$ to VBF production at Tevatron (7 & 8 TeV).
      
    double eVBF_1314_Hbox;///< Theoretical uncertainty in the (linear) new physics contribution from \f$C_{H\Box}\f$ to VBF production at Tevatron (13 & 14 TeV).
    double eVBF_1314_HQ1_11;///< Theoretical uncertainty in the (linear) new physics contribution from \f$(C_{HQ}^{(1)})_{11}\f$ to VBF production at Tevatron (13 & 14 TeV).
    double eVBF_1314_Hu_11;///< Theoretical uncertainty in the (linear) new physics contribution from \f$(C_{Hu})_{11}\f$ to VBF production at Tevatron (13 & 14 TeV).
    double eVBF_1314_Hd_11;///< Theoretical uncertainty in the (linear) new physics contribution from \f$(C_{Hd})_{11}\f$ to VBF production at Tevatron (13 & 14 TeV).
    double eVBF_1314_HQ3_11;///< Theoretical uncertainty in the (linear) new physics contribution from \f$(C_{HQ}^{(3)})_{11}\f$ to VBF production at Tevatron (13 & 14 TeV).
    double eVBF_1314_HD;///< Theoretical uncertainty in the (linear) new physics contribution from \f$C_{HD}\f$ to VBF production at Tevatron (13 & 14 TeV).
    double eVBF_1314_HB;///< Theoretical uncertainty in the (linear) new physics contribution from \f$C_{HB}\f$ to VBF production at Tevatron (13 & 14 TeV).
    double eVBF_1314_HW;///< Theoretical uncertainty in the (linear) new physics contribution from \f$C_{HW}\f$ to VBF production at Tevatron (13 & 14 TeV).
    double eVBF_1314_HWB;///< Theoretical uncertainty in the (linear) new physics contribution from \f$C_{HWB}\f$ to VBF production at Tevatron (13 & 14 TeV).
    double eVBF_1314_HG;///< Theoretical uncertainty in the (linear) new physics contribution from \f$C_{HG}\f$ to VBF production at Tevatron (13 & 14 TeV).
    double eVBF_1314_DHB;///< Theoretical uncertainty in the (linear) new physics contribution from \f$C_{DHB}\f$ to VBF production at Tevatron (13 & 14 TeV).
    double eVBF_1314_DHW;///< Theoretical uncertainty in the (linear) new physics contribution from \f$C_{DHW}\f$ to VBF production at Tevatron (13 & 14 TeV).
    double eVBF_1314_DeltaGF;///< Theoretical uncertainty in the (linear) new physics contribution from \f$\delta_{G_F}\f$ to VBF production at Tevatron (13 & 14 TeV).

    double eWH_2_Hbox;///< Theoretical uncertainty in the (linear) new physics contribution from \f$C_{H\Box}\f$ to WH production at Tevatron (1.96 TeV).
    double eWH_2_HQ3_11;///< Theoretical uncertainty in the (linear) new physics contribution from \f$(C_{HQ}^{(3)})_{11}\f$ to WH production at Tevatron (1.96 TeV).
    double eWH_2_HD;///< Theoretical uncertainty in the (linear) new physics contribution from \f$C_{HD}\f$ to WH production at Tevatron (1.96 TeV).
    double eWH_2_HW;///< Theoretical uncertainty in the (linear) new physics contribution from \f$C_{HW}\f$ to WH production at Tevatron (1.96 TeV).
    double eWH_2_HWB;///< Theoretical uncertainty in the (linear) new physics contribution from \f$C_{HWB}\f$ to WH production at Tevatron (1.96 TeV).
    double eWH_2_DHW;///< Theoretical uncertainty in the (linear) new physics contribution from \f$C_{DHW}\f$ to WH production at the LHC (1.96 TeV).
    double eWH_2_DeltaGF;///< Theoretical uncertainty in the (linear) new physics contribution from \f$\delta_{G_F}\f$ to WH production at the LHC (1.96 TeV).
      
    double eWH_78_Hbox;///< Theoretical uncertainty in the (linear) new physics contribution from \f$C_{H\Box}\f$ to WH production at Tevatron (7 & 8 TeV).
    double eWH_78_HQ3_11;///< Theoretical uncertainty in the (linear) new physics contribution from \f$(C_{HQ}^{(3)})_{11}\f$ to WH production at Tevatron (7 & 8 TeV).
    double eWH_78_HD;///< Theoretical uncertainty in the (linear) new physics contribution from \f$C_{HD}\f$ to WH production at Tevatron (7 & 8 TeV).
    double eWH_78_HW;///< Theoretical uncertainty in the (linear) new physics contribution from \f$C_{HW}\f$ to WH production at Tevatron (7 & 8 TeV).
    double eWH_78_HWB;///< Theoretical uncertainty in the (linear) new physics contribution from \f$C_{HWB}\f$ to WH production at Tevatron (7 & 8 TeV).
    double eWH_78_DHW;///< Theoretical uncertainty in the (linear) new physics contribution from \f$C_{DHW}\f$ to WH production at the LHC (7 & 8 TeV).
    double eWH_78_DeltaGF;///< Theoretical uncertainty in the (linear) new physics contribution from \f$\delta_{G_F}\f$ to WH production at the LHC (7 & 8 TeV).
      
    double eWH_1314_Hbox;///< Theoretical uncertainty in the (linear) new physics contribution from \f$C_{H\Box}\f$ to WH production at Tevatron (13 & 14 TeV).
    double eWH_1314_HQ3_11;///< Theoretical uncertainty in the (linear) new physics contribution from \f$(C_{HQ}^{(3)})_{11}\f$ to WH production at Tevatron (13 & 14 TeV).
    double eWH_1314_HD;///< Theoretical uncertainty in the (linear) new physics contribution from \f$C_{HD}\f$ to WH production at Tevatron (13 & 14 TeV).
    double eWH_1314_HW;///< Theoretical uncertainty in the (linear) new physics contribution from \f$C_{HW}\f$ to WH production at Tevatron (13 & 14 TeV).
    double eWH_1314_HWB;///< Theoretical uncertainty in the (linear) new physics contribution from \f$C_{HWB}\f$ to WH production at Tevatron (13 & 14 TeV).
    double eWH_1314_DHW;///< Theoretical uncertainty in the (linear) new physics contribution from \f$C_{DHW}\f$ to WH production at the LHC (13 & 14 TeV).
    double eWH_1314_DeltaGF;///< Theoretical uncertainty in the (linear) new physics contribution from \f$\delta_{G_F}\f$ to WH production at the LHC (13 & 14 TeV).
      
    double eZH_2_Hbox;///< Theoretical uncertainty in the (linear) new physics contribution from \f$C_{H\Box}\f$ to ZH production at Tevatron (1.96 TeV).
    double eZH_2_HQ1_11;///< Theoretical uncertainty in the (linear) new physics contribution from \f$(C_{HQ}^{(1)})_{11}\f$ to ZH production at Tevatron (1.96 TeV).
    double eZH_2_Hu_11;///< Theoretical uncertainty in the (linear) new physics contribution from \f$(C_{Hu})_{11}\f$ to ZH production at Tevatron (1.96 TeV).
    double eZH_2_Hd_11;///< Theoretical uncertainty in the (linear) new physics contribution from \f$(C_{Hd})_{11}\f$ to ZH production at Tevatron (1.96 TeV).
    double eZH_2_HQ3_11;///< Theoretical uncertainty in the (linear) new physics contribution from \f$(C_{HQ}^{(3)})_{11}\f$ to ZH production at Tevatron (1.96 TeV).
    double eZH_2_HD;///< Theoretical uncertainty in the (linear) new physics contribution from \f$C_{HD}\f$ to ZH production at Tevatron (1.96 TeV).
    double eZH_2_HB;///< Theoretical uncertainty in the (linear) new physics contribution from \f$C_{HB}\f$ to ZH production at Tevatron (1.96 TeV).
    double eZH_2_HW;///< Theoretical uncertainty in the (linear) new physics contribution from \f$C_{HW}\f$ to ZH production at Tevatron (1.96 TeV).
    double eZH_2_HWB;///< Theoretical uncertainty in the (linear) new physics contribution from \f$C_{HWB}\f$ to ZH production at Tevatron (1.96 TeV).
    double eZH_2_DHB;///< Theoretical uncertainty in the (linear) new physics contribution from \f$C_{DHB}\f$ to ZH production at Tevatron (1.96 TeV).
    double eZH_2_DHW;///< Theoretical uncertainty in the (linear) new physics contribution from \f$C_{DHW}\f$ to ZH production at Tevatron (1.96 TeV).
    double eZH_2_DeltaGF;///< Theoretical uncertainty in the (linear) new physics contribution from \f$\delta_{G_F}\f$ to ZH production at Tevatron (1.96 TeV).
      
    double eZH_78_Hbox;///< Theoretical uncertainty in the (linear) new physics contribution from \f$C_{H\Box}\f$ to ZH production at Tevatron (7 & 8 TeV).
    double eZH_78_HQ1_11;///< Theoretical uncertainty in the (linear) new physics contribution from \f$(C_{HQ}^{(1)})_{11}\f$ to ZH production at Tevatron (7 & 8 TeV).
    double eZH_78_Hu_11;///< Theoretical uncertainty in the (linear) new physics contribution from \f$(C_{Hu})_{11}\f$ to ZH production at Tevatron (7 & 8 TeV).
    double eZH_78_Hd_11;///< Theoretical uncertainty in the (linear) new physics contribution from \f$(C_{Hd})_{11}\f$ to ZH production at Tevatron (7 & 8 TeV).
    double eZH_78_HQ3_11;///< Theoretical uncertainty in the (linear) new physics contribution from \f$(C_{HQ}^{(3)})_{11}\f$ to ZH production at Tevatron (7 & 8 TeV).
    double eZH_78_HD;///< Theoretical uncertainty in the (linear) new physics contribution from \f$C_{HD}\f$ to ZH production at Tevatron (7 & 8 TeV).
    double eZH_78_HB;///< Theoretical uncertainty in the (linear) new physics contribution from \f$C_{HB}\f$ to ZH production at Tevatron (7 & 8 TeV).
    double eZH_78_HW;///< Theoretical uncertainty in the (linear) new physics contribution from \f$C_{HW}\f$ to ZH production at Tevatron (7 & 8 TeV).
    double eZH_78_HWB;///< Theoretical uncertainty in the (linear) new physics contribution from \f$C_{HWB}\f$ to ZH production at Tevatron (7 & 8 TeV).
    double eZH_78_DHB;///< Theoretical uncertainty in the (linear) new physics contribution from \f$C_{DHB}\f$ to ZH production at Tevatron (7 & 8 TeV).
    double eZH_78_DHW;///< Theoretical uncertainty in the (linear) new physics contribution from \f$C_{DHW}\f$ to ZH production at Tevatron (7 & 8 TeV).
    double eZH_78_DeltaGF;///< Theoretical uncertainty in the (linear) new physics contribution from \f$\delta_{G_F}\f$ to ZH production at Tevatron (7 & 8 TeV).
      
    double eZH_1314_Hbox;///< Theoretical uncertainty in the (linear) new physics contribution from \f$C_{H\Box}\f$ to ZH production at Tevatron (13 & 14 TeV).
    double eZH_1314_HQ1_11;///< Theoretical uncertainty in the (linear) new physics contribution from \f$(C_{HQ}^{(1)})_{11}\f$ to ZH production at Tevatron (13 & 14 TeV).
    double eZH_1314_Hu_11;///< Theoretical uncertainty in the (linear) new physics contribution from \f$(C_{Hu})_{11}\f$ to ZH production at Tevatron (13 & 14 TeV).
    double eZH_1314_Hd_11;///< Theoretical uncertainty in the (linear) new physics contribution from \f$(C_{Hd})_{11}\f$ to ZH production at Tevatron (13 & 14 TeV).
    double eZH_1314_HQ3_11;///< Theoretical uncertainty in the (linear) new physics contribution from \f$(C_{HQ}^{(3)})_{11}\f$ to ZH production at Tevatron (13 & 14 TeV).
    double eZH_1314_HD;///< Theoretical uncertainty in the (linear) new physics contribution from \f$C_{HD}\f$ to ZH production at Tevatron (13 & 14 TeV).
    double eZH_1314_HB;///< Theoretical uncertainty in the (linear) new physics contribution from \f$C_{HB}\f$ to ZH production at Tevatron (13 & 14 TeV).
    double eZH_1314_HW;///< Theoretical uncertainty in the (linear) new physics contribution from \f$C_{HW}\f$ to ZH production at Tevatron (13 & 14 TeV).
    double eZH_1314_HWB;///< Theoretical uncertainty in the (linear) new physics contribution from \f$C_{HWB}\f$ to ZH production at Tevatron (13 & 14 TeV).
    double eZH_1314_DHB;///< Theoretical uncertainty in the (linear) new physics contribution from \f$C_{DHB}\f$ to ZH production at Tevatron (13 & 14 TeV).
    double eZH_1314_DHW;///< Theoretical uncertainty in the (linear) new physics contribution from \f$C_{DHW}\f$ to ZH production at Tevatron (13 & 14 TeV).
    double eZH_1314_DeltaGF;///< Theoretical uncertainty in the (linear) new physics contribution from \f$\delta_{G_F}\f$ to ZH production at Tevatron (13 & 14 TeV).

    double ettH_2_HG;///< Theoretical uncertainty in the (linear) new physics contribution from \f$C_{HG}\f$ to ttH production at Tevatron (1.96 TeV).
    double ettH_2_G;///< Theoretical uncertainty in the (linear) new physics contribution from \f$C_{G}\f$ to ttH production at Tevatron (1.96 TeV).
    double ettH_2_uG_33r;///< Theoretical uncertainty in the (linear) new physics contribution from \f$(C_{uG})_{33}\f$ to ttH production at the LHC (1.96 TeV).
    double ettH_2_DeltagHt;///< Theoretical uncertainty in the (linear) new physics contribution from \f$\delta g_{Htt}\f$ to ttH production at the LHC (1.96 TeV).
      
    double ettH_78_HG;///< Theoretical uncertainty in the (linear) new physics contribution from \f$C_{HG}\f$ to ttH production at Tevatron (7 & 8 TeV).
    double ettH_78_G;///< Theoretical uncertainty in the (linear) new physics contribution from \f$C_{G}\f$ to ttH production at Tevatron (7 & 8 TeV).
    double ettH_78_uG_33r;///< Theoretical uncertainty in the (linear) new physics contribution from \f$(C_{uG})_{33}\f$ to ttH production at the LHC (7 & 8 TeV).
    double ettH_78_DeltagHt;///< Theoretical uncertainty in the (linear) new physics contribution from \f$\delta g_{Htt}\f$ to ttH production at the LHC (7 & 8 TeV).
      
    double ettH_1314_HG;///< Theoretical uncertainty in the (linear) new physics contribution from \f$C_{HG}\f$ to ttH production at Tevatron (13 & 14 TeV).
    double ettH_1314_G;///< Theoretical uncertainty in the (linear) new physics contribution from \f$C_{G}\f$ to ttH production at Tevatron (13 & 14 TeV).
    double ettH_1314_uG_33r;///< Theoretical uncertainty in the (linear) new physics contribution from \f$(C_{uG})_{33}\f$ to ttH production at the LHC (13 & 14 TeV).
    double ettH_1314_DeltagHt;///< Theoretical uncertainty in the (linear) new physics contribution from \f$\delta g_{Htt}\f$ to ttH production at the LHC (13 & 14 TeV).

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
     * @brief The diagonal entry of the dimension-6 operator coefficient \f$C_{Hud}\f$ corresponding to particle f.
     * @param[in] u a quark
     * @return @f$(\f$C_{Hud})_{ud}\f$@f$
     */
    gslpp::complex CHud_diag(const Particle u) const;
    
    /**
     * @brief The diagonal entry of the dimension-6 operator coefficient \f$C_{EH,UH,DH}\f$ corresponding to particle f.
     * @param[in] f a lepton or quark
     * @return @f$(\f$C_{fH})_{ff}\f$@f$
     */
    gslpp::complex CfH_diag(const Particle f) const;
      
    /**
     * @brief The diagonal entry of the dimension-6 operator coefficient \f$C_{UG,DG}\f$ corresponding to particle f.
     * @param[in] f a lepton or quark
     * @return @f$(\f$C_{fG})_{ff}\f$@f$
     */
    gslpp::complex CfG_diag(const Particle f) const;
      
    /**
     * @brief The diagonal entry of the dimension-6 operator coefficient \f$C_{EW,UW,DW}\f$ corresponding to particle f.
     * @param[in] f a lepton or quark
     * @return @f$(\f$C_{fW})_{ff}\f$@f$
     */
    gslpp::complex CfW_diag(const Particle f) const;
      
    /**
     * @brief The diagonal entry of the dimension-6 operator coefficient \f$C_{EB,UB,DB}\f$ corresponding to particle f.
     * @param[in] f a lepton or quark
     * @return @f$(\f$C_{fB})_{ff}\f$@f$
     */
    gslpp::complex CfB_diag(const Particle f) const;


    ////////////////////////////////////////////////////////////////////////
private:

    bool FlagQuadraticTerms; ///< A boolean flag that is true if the quadratic terms in cross sections and widths are switched on.
    bool FlagRotateCHWCHB; ///< A boolean flag that is true if we use as parameters CHWHB_gaga and CHWHB_gagaorth instead of CHW and CHB.

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

#endif	/* NPSMEFTD6_H */

