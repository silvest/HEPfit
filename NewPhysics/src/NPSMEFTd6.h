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
#include <gsl/gsl_integration.h>


/**
 * @class NPSMEFTd6
 * @brief A model class for new physics in the form of the
 * dimension-six effective Lagrangian.
 * @ingroup NewPhysics
 * @author HEPfit Collaboration
 * @copyright GNU General Public License
 * @details This is a Model class containing parameters and functions
 * associated with the general dimension-six effective Lagrangian.
 * (Use the model name "NPSMEFTd6_LFU_QFU" to assume lepton and quark flavour universality)
 *
 * In this class we consider the dimension-six effective Lagrangian
 * @f[
 * \mathcal{L}_\mathrm{eff}
 * = \mathcal{L}_\mathrm{SM}
 *   + \sum_i \frac{C_i}{\Lambda^2} \mathcal{O}_i.
 * @f]
 * The implementation is written in the basis of \cite Grzadkowski:2010es.
 * For convenience, the parameterization also includes operators appearing in
 * other common bases. In particular, the complete set of parameters contains 4
 * redundancies, given by the coefficients \f$C_{2B,2W,DHB,DHW,DB,DW} \f$,
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
 *   <td class="mod_name">%C2B </td>
 *   <td class="mod_symb">\f$C_{2B} \f$</td>
 *   <td class="mod_desc">The coefficient of the operator
 *   \f${\cal O}_{2B}=\frac 12 (\partial_\rho B_{\mu\nu})^2\f$. 
 *  (Implemented via EOM.)</td>
 * </tr>
 * <tr>
 *   <td class="mod_name">%C2W </td>
 *   <td class="mod_symb">\f$C_{2W} \f$</td>
 *   <td class="mod_desc">The coefficient of the operator
 *   \f${\cal O}_{2W}=\frac 12 (D_\rho W_{\mu\nu}^{a})^2\f$. 
 *  (Implemented via EOM.)</td>
 * </tr>
 * <tr>
 *   <td class="mod_name">%C2BS </td>
 *   <td class="mod_symb">\f$C_{2B}^{SILH} \f$</td>
 *   <td class="mod_desc">The coefficient of the SILH operator
 *   \f${\cal O}_{2B}^{SILH}=\frac 12 (\partial^\mu B_{\mu\nu})(\partial_\rho B^{\rho\nu})\f$. 
 *   (Implemented via EOM.) </td>
 * </tr>
 * <tr>
 *   <td class="mod_name">%C2WS </td>
 *   <td class="mod_symb">\f$C_{2W}^{SILH} \f$</td>
 *   <td class="mod_desc">The coefficient of the operator
 *   \f${\cal O}_{2W}^{SILH}=\frac 12 (D_\mu W^{a~\!\mu\nu})(D^\rho W_{\rho\nu}^{a})\f$. 
 *   (Implemented via EOM.) </td>
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
 *   \f${\cal O}_{DHW}=i\big(D^\mu H^\dagger \sigma^a D^\nu H\big) W_{\mu\nu}^a\f$. </td>
 * </tr>
 * <tr>
 *   <td class="mod_name">%CDB </td>
 *   <td class="mod_symb">\f$C_{DB} \f$</td>
 *   <td class="mod_desc">The coefficient of the operator
 *   \f${\cal O}_{DB}=\frac{i}{2}\big(H^\dagger \overset{\leftrightarrow}{D}^\mu H\big) \partial^\nu B_{\mu\nu}\f$. 
 *   (Implemented via EOM.)</td>
 * </tr>
 * <tr>
 *   <td class="mod_name">%CDW </td>
 *   <td class="mod_symb">\f$C_{DW}\f$</td>
 *   <td class="mod_desc">The coefficient of the operator
 *   \f${\cal O}_{DW}=\frac{i}{2}\big(H^\dagger \overset{\leftrightarrow}{D}^{a~\!\mu} H\big) D^\nu W_{\mu\nu}^a\f$.
 *  (Implemented via EOM.) </td>
 * </tr>
 * <tr>
 *   <td class="mod_name">%CWB </td>
 *   <td class="mod_symb">\f$C_{WB} \f$</td>
 *   <td class="mod_desc">The coefficient of the operator
 *   \f${\cal O}_{HWB}=\big(H^\dagger\sigma^a H\big)W_{\mu\nu}^a B^{\mu\nu}\f$. </td>
 * </tr>
 * <tr>
 *   <td class="mod_name">%CHD </td>
 *   <td class="mod_symb">\f$C_{HD}\f$</td>
 *   <td class="mod_desc">The coefficient of the operator
 *   \f${\cal O}_{HD}=\big|H^\dagger D_\mu H\big|^2\f$. </td>
 * </tr>
 * <tr>
 *   <td class="mod_name">%CT </td>
 *   <td class="mod_symb">\f$C_{T}\f$</td>
 *   <td class="mod_desc">The coefficient of the operator
 *   \f${\cal O}_{T}=\frac{1}{2} \big(H^\dagger \overset{\leftrightarrow}{D}_\mu H\big)^2\f$. </td>
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
 *   \big(\overline{L^i}\,\gamma^\mu \sigma^a L^j\big)\f$, for \f$i,j=1,2,3\f$. </td>
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
 *   \big(\overline{Q^i}\,\gamma^\mu \sigma^a Q^j\big)\f$, for \f$i,j=1,2,3\f$. </td>
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
 *   \f$({\cal O}_{uW})_{ij} =\big(\overline{Q^i}\sigma^{\mu\nu} \sigma_a U^j\big)\widetilde{H} W_{\mu\nu}^a\f$, for \f$i,j=1,2,3\f$. </td>
 * </tr>
 * <tr>
 *   <td class="mod_name">%CuB_klr, CuB_kli </td>
 *   <td class="mod_symb">\f$\mbox{Re}\big[(C_{uB})_{kl}\big], \mbox{Im}\big[(C_{uB})_{kl}\big] \f$</td>
 *   <td class="mod_desc">The real and imaginary parts of the coefficient of the operator
 *   \f$({\cal O}_{uB})_{ij} =\big(\overline{Q^i}\sigma^{\mu\nu} U^j\big)\widetilde{H} B_{\mu\nu}\f$, for \f$i,j=1,2,3\f$. </td>
 * </tr>
 * <tr>
 *   <td class="mod_name">%CdG_klr, CdG_kli </td>
 *   <td class="mod_symb">\f$\mbox{Re}\big[(C_{dG})_{kl}\big], \mbox{Im}\big[(C_{dG})_{kl}\big] \f$</td>
 *   <td class="mod_desc">The real and imaginary parts of the coefficient of the operator
 *   \f$({\cal O}_{dG})_{ij} =\big(\overline{Q^i}\sigma^{\mu\nu} T_A D^j\big)H G_{\mu\nu}^A\f$, for \f$i,j=1,2,3\f$. </td>
 * </tr>
 * <tr>
 *   <td class="mod_name">%CdW_klr, CdW_kli </td>
 *   <td class="mod_symb">\f$\mbox{Re}\big[(C_{dW})_{kl}\big], \mbox{Im}\big[(C_{dW})_{kl}\big] \f$</td>
 *   <td class="mod_desc">The real and imaginary parts of the coefficient of the operator
 *   \f$({\cal O}_{dW})_{ij} =\big(\overline{Q^i}\sigma^{\mu\nu} \sigma_a D^j\big)H W_{\mu\nu}^a\f$, for \f$i,j=1,2,3\f$. </td>
 * </tr>
 * <tr>
 *   <td class="mod_name">%CdB_klr, CdB_kli </td>
 *   <td class="mod_symb">\f$\mbox{Re}\big[(C_{dB})_{kl}\big], \mbox{Im}\big[(C_{dB})_{kl}\big] \f$</td>
 *   <td class="mod_desc">The real and imaginary parts of the coefficient of the operator
 *   \f$({\cal O}_{dB})_{ij} =\big(\overline{Q^i}\sigma^{\mu\nu} D^j\big)H B_{\mu\nu}\f$, for \f$i,j=1,2,3\f$. </td>
 * </tr>
 * <tr>
 *   <td class="mod_name">%CeW_klr, CeW_kli </td>
 *   <td class="mod_symb">\f$\mbox{Re}\big[(C_{eW})_{kl}\big], \mbox{Im}\big[(C_{eW})_{kl}\big] \f$</td>
 *   <td class="mod_desc">The real and imaginary parts of the coefficient of the operator
 *   \f$({\cal O}_{eW})_{ij} =\big(\overline{L^i}\sigma^{\mu\nu} \sigma_a E^j\big)H W_{\mu\nu}^a\f$, for \f$i,j=1,2,3\f$. </td>
 * </tr>
 * <tr>
 *   <td class="mod_name">%CeB_klr, CeB_kli </td>
 *   <td class="mod_symb">\f$\mbox{Re}\big[(C_{eB})_{kl}\big], \mbox{Im}\big[(C_{eB})_{kl}\big] \f$</td>
 *   <td class="mod_desc">The real and imaginary parts of the coefficient of the operator
 *   \f$({\cal O}_{eB})_{ij} =\big(\overline{L^i}\sigma^{\mu\nu} E^j\big)H B_{\mu\nu}\f$, for \f$i,j=1,2,3\f$. </td>
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
 *   \f$({\cal O}_{LQ}^{(3)})_{ijkl}=\big(\overline{L^i}\,\gamma^\mu \sigma_a L^j\big)
 *   \big(\overline{Q^k}\,\gamma_\mu \sigma_a Q^l\big)\f$. </td>
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
 *   <td class="mod_name">%CQQ1_1133, CQQ1_1331 </td>
 *   <td class="mod_symb">\f$C_{QQ}^{(1)}\f$</td>
 *   <td class="mod_desc">The coefficient of the operator
 *   \f$({\cal O}_{QQ}^{(1)})_{ijkl}=\big(\overline{Q^i}\,\gamma^\mu Q^j\big)
 *   \big(\overline{Q^k}\,\gamma_\mu Q^l\big)\f$. </td>
 * </tr>
 * <tr>
 *   <td class="mod_name">%CQQ3_1133, CQQ3_1331 </td>
 *   <td class="mod_symb">\f$C_{QQ}^{(3)}\f$</td>
 *   <td class="mod_desc">The coefficient of the operator
 *   \f$({\cal O}_{QQ}^{(3)})_{ijkl}=\big(\overline{Q^i}\,\gamma^\mu \sigma_a Q^j\big)
 *   \big(\overline{Q^k}\,\gamma_\mu \sigma_a Q^l\big)\f$. </td>
 * </tr>
 * <tr>
 *   <td class="mod_name">%Cuu_1133, Cuu_1331 </td>
 *   <td class="mod_symb">\f$C_{uu}\f$</td>
 *   <td class="mod_desc">The coefficient of the operator
 *   \f$({\cal O}_{uu})_{ijkl}=\big(\overline{u^i}\,\gamma^\mu u^j\big)
 *   \big(\overline{u^k}\,\gamma_\mu u^l\big)\f$. </td>
 * </tr>
 * <tr>
 *   <td class="mod_name">%Cud1_3311 </td>
 *   <td class="mod_symb">\f$C_{ud}^{(1)}\f$</td>
 *   <td class="mod_desc">The coefficient of the operator
 *   \f$({\cal O}_{ud}^{(1)})_{ijkl}=\big(\overline{u^i}\,\gamma^\mu u^j\big)
 *   \big(\overline{d^k}\,\gamma_\mu d^l\big)\f$. </td>
 * </tr>
 * <tr>
 *   <td class="mod_name">%Cud8_3311 </td>
 *   <td class="mod_symb">\f$C_{ud}^{(8)}\f$</td>
 *   <td class="mod_desc">The coefficient of the operator
 *   \f$({\cal O}_{ud}^{(8)})_{ijkl}=\big(\overline{u^i}\,\gamma^\mu T^A u^j\big)
 *   \big(\overline{d^k}\,\gamma_\mu T^A d^l\big)\f$. </td>
 * </tr>
 * <tr>
 *   <td class="mod_name">%CQu1_1133, CQu1_3311 </td>
 *   <td class="mod_symb">\f$C_{Qu}^{(1)}\f$</td>
 *   <td class="mod_desc">The coefficient of the operator
 *   \f$({\cal O}_{Qu}^{(1)})_{ijkl}=\big(\overline{Q^i}\,\gamma^\mu Q^j\big)
 *   \big(\overline{u^k}\,\gamma_\mu u^l\big)\f$. </td>
 * </tr>
 * <tr>
 *   <td class="mod_name">%CQu8_1133, CQu8_3311 </td>
 *   <td class="mod_symb">\f$C_{Qu}^{(8)}\f$</td>
 *   <td class="mod_desc">The coefficient of the operator
 *   \f$({\cal O}_{Qu}^{(8)})_{ijkl}=\big(\overline{Q^i}\,\gamma^\mu T^A Q^j\big)
 *   \big(\overline{u^k}\,\gamma_\mu T^A u^l\big)\f$. </td>
 * </tr>
 * <tr>
 *   <td class="mod_name">%CQd1_3311 </td>
 *   <td class="mod_symb">\f$C_{Qd}^{(1)}\f$</td>
 *   <td class="mod_desc">The coefficient of the operator
 *   \f$({\cal O}_{Qd}^{(1)})_{ijkl}=\big(\overline{Q^i}\,\gamma^\mu Q^j\big)
 *   \big(\overline{d^k}\,\gamma_\mu d^l\big)\f$. </td>
 * </tr>
 * <tr>
 *   <td class="mod_name">%CQd8_3311 </td>
 *   <td class="mod_symb">\f$C_{Qd}^{(8)}\f$</td>
 *   <td class="mod_desc">The coefficient of the operator
 *   \f$({\cal O}_{Qd}^{(8)})_{ijkl}=\big(\overline{Q^i}\,\gamma^\mu T^A Q^j\big)
 *   \big(\overline{d^k}\,\gamma_\mu T^A d^l\big)\f$. </td>
 * </tr> 
 * <tr>
 *   <td class="mod_name">%Lambda_NP </td>
 *   <td class="mod_symb">\f$\Lambda \f$</td>
 *   <td class="mod_desc">The new physics scale. </td>
 * </tr>
 * <tr>
 *   <td class="mod_name">%BrHinv</td>
 *   <td class="mod_symb">Br@f$(H\to invisible)@f$</td>
 *   <td class="mod_desc">The branching ratio of invisible Higgs decays. Only the absolute value of this parameter is considered.(Not part of the EFT. Only for tests.)</td>
 * </tr>
 * <tr>
 *   <td class="mod_name">%BrHexo</td>
 *   <td class="mod_symb">Br@f$(H\to exotic)@f$</td>
 *   <td class="mod_desc">The branching ratio of exotic Higgs decays. Only the absolute value of this parameter is considered. (Not part of the EFT. Only for tests.)</td>
 * </tr>
 * <tr>
 *   <td class="mod_name">%dg1Z</td>
 *   <td class="mod_symb">@f$\delta g_{1Z}@f$</td>
 *   <td class="mod_desc">Independent contribution to aTGC. (extra contribution to the one from the EFT. Only for tests.)</td>
 * </tr>
  * <tr>
 *   <td class="mod_name">%dKappaga</td>
 *   <td class="mod_symb">@f$\delta \kappa_{\gamma}@f$</td>
 *   <td class="mod_desc">Independent contribution to aTGC. (extra contribution to the one from the EFT. Only for tests.)</td>
 * </tr>
  * <tr>
 *   <td class="mod_name">%lambZ</td>
 *   <td class="mod_symb">@f$\lambda_{Z}@f$</td>
 *   <td class="mod_desc">Independent contribution to aTGC. (extra contribution to the one from the EFT. Only for tests.)</td>
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
 * H^\dagger i \overset{\leftrightarrow}{D^a_\mu} H\equiv H^\dagger i (\sigma^a D_\mu - \overset{\leftarrow}{D}_\mu \sigma^a)H.  
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
 *   <td class="mod_name">%C2B </td>
 *   <td class="mod_symb">\f$C_{2B} \f$</td>
 *   <td class="mod_desc">The coefficient of the operator
 *   \f${\cal O}_{2B}=\frac 12 (\partial_\rho B_{\mu\nu})^2\f$. 
 *   (Implemented via EOM.) </td>
 * </tr>
 * <tr>
 *   <td class="mod_name">%C2W </td>
 *   <td class="mod_symb">\f$C_{2W} \f$</td>
 *   <td class="mod_desc">The coefficient of the operator
 *   \f${\cal O}_{2W}=\frac 12 (D_\rho W_{\mu\nu}^{a})^2\f$. 
 *   (Implemented via EOM.) </td>
 * </tr>
 * <tr>
 *   <td class="mod_name">%C2BS </td>
 *   <td class="mod_symb">\f$C_{2B}^{SILH} \f$</td>
 *   <td class="mod_desc">The coefficient of the SILH operator
 *   \f${\cal O}_{2B}^{SILH}=\frac 12 (\partial^\mu B_{\mu\nu})(\partial_\rho B^{\rho\nu})\f$. 
 *   (Implemented via EOM.) </td>
 * </tr>
 * <tr>
 *   <td class="mod_name">%C2WS </td>
 *   <td class="mod_symb">\f$C_{2W}^{SILH} \f$</td>
 *   <td class="mod_desc">The coefficient of the operator
 *   \f${\cal O}_{2W}^{SILH}=\frac 12 (D_\mu W^{a~\!\mu\nu})(D^\rho W_{\rho\nu}^{a})\f$. 
 *   (Implemented via EOM.) </td>
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
 *   \f${\cal O}_{DHW}=i\big(D^\mu H^\dagger \sigma^a D^\nu H\big) W_{\mu\nu}^a\f$. </td>
 * </tr>
 * <tr>
 *   <td class="mod_name">%CDB </td>
 *   <td class="mod_symb">\f$C_{DB} \f$</td>
 *   <td class="mod_desc">The coefficient of the operator
 *   \f${\cal O}_{DB}=\frac{i}{2}\big(H^\dagger \overset{\leftrightarrow}{D}^\mu H\big) \partial^\nu B_{\mu\nu}\f$. 
 *   (Implemented via EOM.)</td>
 * </tr>
 * <tr>
 *   <td class="mod_name">%CDW </td>
 *   <td class="mod_symb">\f$C_{DW}\f$</td>
 *   <td class="mod_desc">The coefficient of the operator
 *   \f${\cal O}_{DW}=\frac{i}{2}\big(H^\dagger \overset{\leftrightarrow}{D}^{a~\!\mu} H\big) D^\nu W_{\mu\nu}^a\f$.
 *  (Implemented via EOM.) </td>
 * </tr>
 * <tr>
 *   <td class="mod_name">%CWB </td>
 *   <td class="mod_symb">\f$C_{WB} \f$</td>
 *   <td class="mod_desc">The coefficient of the operator
 *   \f${\cal O}_{HWB}=\big(H^\dagger\sigma^a H\big)W_{\mu\nu}^a B^{\mu\nu}\f$. </td>
 * </tr>
 * <tr>
 *   <td class="mod_name">%CHD </td>
 *   <td class="mod_symb">\f$C_{HD}\f$</td>
 *   <td class="mod_desc">The coefficient of the operator
 *   \f${\cal O}_{HD}=\big|H^\dagger D_\mu H\big|^2\f$. </td>
 * </tr>
 * <tr>
 *   <td class="mod_name">%CT </td>
 *   <td class="mod_symb">\f$C_{T}\f$</td>
 *   <td class="mod_desc">The coefficient of the operator
 *   \f${\cal O}_{T}=\frac{1}{2} \big(H^\dagger \overset{\leftrightarrow}{D}_\mu H\big)^2\f$. </td>
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
 *   \big(\overline{L^i}\,\gamma^\mu \sigma^a L^i\big)\f$ (flavor universal). </td>
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
 *   \big(\overline{Q^i}\,\gamma^\mu \sigma^a Q^i\big)\f$ (flavor universal). </td>
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
 *   <td class="mod_name">%CuG_jjr, CuG_jji </td>
 *   <td class="mod_symb">\f$\mbox{Re}\big[(C_{uG})_{jj}\big], \mbox{Im}\big[(C_{uG})_{jj}\big] \f$</td>
 *   <td class="mod_desc">The real and imaginary parts of the coefficient of the operator
 *   \f$({\cal O}_{uG})_{jj} =\big(\overline{Q^j}\sigma^{\mu\nu} T_A U^j\big)\widetilde{H} G_{\mu\nu}^A\f$ (flavor universal). </td>
 * </tr>
 * <tr>
 *   <td class="mod_name">%CuW_jjr, CuW_jji </td>
 *   <td class="mod_symb">\f$\mbox{Re}\big[(C_{uW})_{jj}\big], \mbox{Im}\big[(C_{uW})_{jj}\big] \f$</td>
 *   <td class="mod_desc">The real and imaginary parts of the coefficient of the operator
 *   \f$({\cal O}_{uW})_{jj} =\big(\overline{Q^j}\sigma^{\mu\nu} \sigma_a U^j\big)\widetilde{H} W_{\mu\nu}^a\f$ (flavor universal). </td>
 * </tr>
 * <tr>
 *   <td class="mod_name">%CuB_jjr, CuB_jji </td>
 *   <td class="mod_symb">\f$\mbox{Re}\big[(C_{uB})_{jj}\big], \mbox{Im}\big[(C_{uB})_{jj}\big] \f$</td>
 *   <td class="mod_desc">The real and imaginary parts of the coefficient of the operator
 *   \f$({\cal O}_{uB})_{jj} =\big(\overline{Q^j}\sigma^{\mu\nu} U^j\big)\widetilde{H} B_{\mu\nu}\f$ (flavor universal). </td>
 * </tr>
 * <tr>
 *   <td class="mod_name">%CdG_jjr, CdG_jji </td>
 *   <td class="mod_symb">\f$\mbox{Re}\big[(C_{dG})_{jj}\big], \mbox{Im}\big[(C_{dG})_{jj}\big] \f$</td>
 *   <td class="mod_desc">The real and imaginary parts of the coefficient of the operator
 *   \f$({\cal O}_{dG})_{jj} =\big(\overline{Q^j}\sigma^{\mu\nu} T_A D^j\big)H G_{\mu\nu}^A\f$ (flavor universal). </td>
 * </tr>
 * <tr>
 *   <td class="mod_name">%CdW_jjr, CdW_jji </td>
 *   <td class="mod_symb">\f$\mbox{Re}\big[(C_{dW})_{jj}\big], \mbox{Im}\big[(C_{dW})_{jj}\big] \f$</td>
 *   <td class="mod_desc">The real and imaginary parts of the coefficient of the operator
 *   \f$({\cal O}_{dW})_{jj} =\big(\overline{Q^j}\sigma^{\mu\nu} \sigma_a D^j\big)H W_{\mu\nu}^a\f$ (flavor universal). </td>
 * </tr>
 * <tr>
 *   <td class="mod_name">%CdB_jjr, CdB_jji </td>
 *   <td class="mod_symb">\f$\mbox{Re}\big[(C_{dB})_{jj}\big], \mbox{Im}\big[(C_{dB})_{jj}\big] \f$</td>
 *   <td class="mod_desc">The real and imaginary parts of the coefficient of the operator
 *   \f$({\cal O}_{dB})_{jj} =\big(\overline{Q^j}\sigma^{\mu\nu} D^j\big)H B_{\mu\nu}\f$ (flavor universal). </td>
 * </tr>
 * <tr>
 *   <td class="mod_name">%CeW_jjr, CeW_jji </td>
 *   <td class="mod_symb">\f$\mbox{Re}\big[(C_{eW})_{jj}\big], \mbox{Im}\big[(C_{eW})_{jj}\big] \f$</td>
 *   <td class="mod_desc">The real and imaginary parts of the coefficient of the operator
 *   \f$({\cal O}_{eW})_{jj} =\big(\overline{L^j}\sigma^{\mu\nu} \sigma_a E^j\big)H W_{\mu\nu}^a\f$ (flavor universal). </td>
 * </tr>
 * <tr>
 *   <td class="mod_name">%CeB_jjr, CeB_jji </td>
 *   <td class="mod_symb">\f$\mbox{Re}\big[(C_{eB})_{jj}\big], \mbox{Im}\big[(C_{eB})_{jj}\big] \f$</td>
 *   <td class="mod_desc">The real and imaginary parts of the coefficient of the operator
 *   \f$({\cal O}_{eB})_{jj} =\big(\overline{L^j}\sigma^{\mu\nu} E^j\big)H B_{\mu\nu}\f$ (flavor universal). </td>
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
 *   \f$({\cal O}_{LQ}^{(3)})_{ijkl}=\big(\overline{L^i}\,\gamma^\mu \sigma_a L^j\big)
 *   \big(\overline{Q^k}\,\gamma_\mu \sigma_a Q^l\big)\f$. </td>
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
 *   <td class="mod_name">%CQQ1 </td>
 *   <td class="mod_symb">\f$C_{QQ}^{(1)}\f$</td>
 *   <td class="mod_desc">The coefficient of the operator
 *   \f$({\cal O}_{QQ}^{(1)})_{ijkl}=\big(\overline{Q^i}\,\gamma^\mu Q^j\big)
 *   \big(\overline{Q^k}\,\gamma_\mu Q^l\big)\f$. </td>
 * </tr>
 * <tr>
 *   <td class="mod_name">%CQQ3 </td>
 *   <td class="mod_symb">\f$C_{QQ}^{(3)}\f$</td>
 *   <td class="mod_desc">The coefficient of the operator
 *   \f$({\cal O}_{QQ}^{(3)})_{ijkl}=\big(\overline{Q^i}\,\gamma^\mu \sigma_a Q^j\big)
 *   \big(\overline{Q^k}\,\gamma_\mu \sigma_a Q^l\big)\f$. </td>
 * </tr>
 * <tr>
 *   <td class="mod_name">%Cuu </td>
 *   <td class="mod_symb">\f$C_{uu}\f$</td>
 *   <td class="mod_desc">The coefficient of the operator
 *   \f$({\cal O}_{uu})_{ijkl}=\big(\overline{u^i}\,\gamma^\mu u^j\big)
 *   \big(\overline{u^k}\,\gamma_\mu u^l\big)\f$. </td>
 * </tr>
 * <tr>
 *   <td class="mod_name">%Cud1 </td>
 *   <td class="mod_symb">\f$C_{ud}^{(1)}\f$</td>
 *   <td class="mod_desc">The coefficient of the operator
 *   \f$({\cal O}_{ud}^{(1)})_{ijkl}=\big(\overline{u^i}\,\gamma^\mu u^j\big)
 *   \big(\overline{d^k}\,\gamma_\mu d^l\big)\f$. </td>
 * </tr>
 * <tr>
 *   <td class="mod_name">%Cud8 </td>
 *   <td class="mod_symb">\f$C_{ud}^{(8)}\f$</td>
 *   <td class="mod_desc">The coefficient of the operator
 *   \f$({\cal O}_{ud}^{(8)})_{ijkl}=\big(\overline{u^i}\,\gamma^\mu T^A u^j\big)
 *   \big(\overline{d^k}\,\gamma_\mu T^A d^l\big)\f$. </td>
 * </tr>
 * <tr>
 *   <td class="mod_name">%CQu1 </td>
 *   <td class="mod_symb">\f$C_{Qu}^{(1)}\f$</td>
 *   <td class="mod_desc">The coefficient of the operator
 *   \f$({\cal O}_{Qu}^{(1)})_{ijkl}=\big(\overline{Q^i}\,\gamma^\mu Q^j\big)
 *   \big(\overline{u^k}\,\gamma_\mu u^l\big)\f$. </td>
 * </tr>
 * <tr>
 *   <td class="mod_name">%CQu8 </td>
 *   <td class="mod_symb">\f$C_{Qu}^{(8)}\f$</td>
 *   <td class="mod_desc">The coefficient of the operator
 *   \f$({\cal O}_{Qu}^{(8)})_{ijkl}=\big(\overline{Q^i}\,\gamma^\mu T^A Q^j\big)
 *   \big(\overline{u^k}\,\gamma_\mu T^A u^l\big)\f$. </td>
 * </tr>
 * <tr>
 *   <td class="mod_name">%CQd1 </td>
 *   <td class="mod_symb">\f$C_{Qd}^{(1)}\f$</td>
 *   <td class="mod_desc">The coefficient of the operator
 *   \f$({\cal O}_{Qd}^{(1)})_{ijkl}=\big(\overline{Q^i}\,\gamma^\mu Q^j\big)
 *   \big(\overline{d^k}\,\gamma_\mu d^l\big)\f$. </td>
 * </tr>
 * <tr>
 *   <td class="mod_name">%CQd8 </td>
 *   <td class="mod_symb">\f$C_{Qd}^{(8)}\f$</td>
 *   <td class="mod_desc">The coefficient of the operator
 *   \f$({\cal O}_{Qd}^{(8)})_{ijkl}=\big(\overline{Q^i}\,\gamma^\mu T^A Q^j\big)
 *   \big(\overline{d^k}\,\gamma_\mu T^A d^l\big)\f$. </td>
 * </tr>
 * <tr>
 *   <td class="mod_name">%Lambda_NP </td>
 *   <td class="mod_symb">\f$\Lambda \f$</td>
 *   <td class="mod_desc">The new physics scale. </td>
 * </tr>
 * <tr>
 *   <td class="mod_name">%BrHinv</td>
 *   <td class="mod_symb">Br@f$(H\to invisible)@f$</td>
 *   <td class="mod_desc">The branching ratio of invisible Higgs decays. (Not part of the EFT. Only for tests.)</td>
 * </tr>
 * <tr>
 *   <td class="mod_name">%BrHexo</td>
 *   <td class="mod_symb">Br@f$(H\to exotic)@f$</td>
 *   <td class="mod_desc">The branching ratio of exotic Higgs decays. (Not part of the EFT. Only for tests.)</td>
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
 * <tr>
 *   <td class="mod_name">%PartialQFU</td>
 *   <td class="mod_valu">TRUE&nbsp;/&nbsp;<b>FALSE</b></td>
 *   <td class="mod_desc">This flag is set to TRUE if using CHQ1_11=CHQ1_22, 
 *   CHQ3_11=CHQ3_22, CHU_11=CHU_22, CHD_11=CHD_22, CHud_11=CHud_22.}
 *   Only applies in the Non QFU case. In that case only the (1,1) component is taken into account.
 *   The default value is FALSE.</td>
 * </tr>
 * <tr>
 *   <td class="mod_name">%FlavU3OfX</td>
 *   <td class="mod_valu">TRUE&nbsp;/&nbsp;<b>FALSE</b></td>
 *   <td class="mod_desc">This flag is set to TRUE if using @f$U(3)^5@f$ flavour symmetry relations
 *   in the coefficients of the operators @f$O_{fH}@f$ and @f$O_{fV}@f$. If TRUE, the operator coefficient is proportional
 *   to the corresponding Yukawa matrix (diagonal), with the proportionality coefficient given by 
 *   the Model parameter corresponding to the coefficient of third family.
 *   (Implemented only for the real and diagonal elements of the @f$O_{fH}@f$ and @f$O_{fV}@f$ operators.)
 *   The default value is FALSE.</td>
 * </tr>
 * <tr>
 *   <td class="mod_name">%FlagUnivOfX</td>
 *   <td class="mod_valu">TRUE&nbsp;/&nbsp;<b>FALSE</b></td>
 *   <td class="mod_desc">This flag is set to TRUE if using @f$U(3)^5@f$ flavour symmetry relations
 *   in the coefficients of the operators @f$O_{fH}@f$ and @f$O_{fV}@f$ plus they are the same for all fermions. 
 *   If TRUE, all the operator coefficients are proportional
 *   to the corresponding Yukawa matrix (diagonal), with the proportionality coefficient given by 
 *   the Model parameter corresponding to the coefficients of third family for @f$O_{uH}@f$ and @f$O_{uV}@f$, respectively.
 *   (Implemented only for the real and diagonal elements of the @f$O_{fH}@f$ and @f$O_{fV}@f$ operators.)
 *   The default value is FALSE.</td>
 * </tr>
 * <tr>
 *   <td class="mod_name">%HiggsSM</td>
 *   <td class="mod_valu">TRUE&nbsp;/&nbsp;<b>FALSE</b></td>
 *   <td class="mod_desc">This flag is set to TRUE if including dependence on small variations of the SM parameters (dependence is linearized). 
 *   Available only in selected Higgs observables.
 *   The default value is FALSE.</td>
 * </tr>
 * <tr>
 *   <td class="mod_name">%LoopHd6</td>
 *   <td class="mod_valu">TRUE&nbsp;/&nbsp;<b>FALSE</b></td>
 *   <td class="mod_desc">This flag is set to TRUE if including modifications in the SM loops in Higgs observables due to the dim 6 interactions.
 *   The default value is FALSE.</td>
 * </tr>
 * <tr>
 *   <td class="mod_name">%LoopH3d6Quad</td>
 *   <td class="mod_valu">TRUE&nbsp;/&nbsp;<b>FALSE</b></td>
 *   <td class="mod_desc">This flag is set to TRUE if including quadratic modifications in the SM loops in Higgs observables due 
 *   to the dim 6 interactions that contribute to the trilinear Higgs coupling. 
 *   Works independently of the flag QuadraticTerms (the quadratic contributions are also added if the latter is true).
 *   The default value is FALSE.</td>
 * </tr>
 * <tr>
 *   <td class="mod_name">%RGEciLLA</td>
 *   <td class="mod_valu">TRUE&nbsp;/&nbsp;<b>FALSE</b></td>
 *   <td class="mod_desc">This flag is set to TRUE if including log-enhanced 1-loop corrections proportional to the dim-6 Wilson coefficients via renormalization group effects. 
 *   Only valid working exactly inside the Warsaw basis. Logs evaluated at the top mass scale. (NOT ACTIVE YET. PLACEHOLDER.)
 *   The default value is FALSE.</td>
 * </tr>
 * 
 * 
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
    static const int NNPSMEFTd6Vars = 522;

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
    static const int NNPSMEFTd6Vars_LFU_QFU = 271;

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
    * @brief Return Lambda_NP
    * @return \f$ Lambda_NP \f$ 
    */
    double getLambda_NP() const{
        return Lambda_NP;
    }

    /**
    * @brief Return NP coeff CQe_2322
    * @return \f$ CQe_2322 \f$ 
    */
    double getCQe_2311() const{
        return CQe_2311;
    }
    
    /**
    * @brief Return NP coeff CQe_2322
    * @return \f$ CQe_2311 \f$ 
    */
    double getCQe_2322() const{
        return CQe_2322;
    }

    /**
    * @brief Return NP coeff CLQ1_1123
    * @return \f$ CLQ1_1123 \f$ 
    */
    double getCLQ1_1123() const{
        return CLQ1_1123;
    }
    
    /**
    * @brief Return NP coeff CLQ1_2223
    * @return \f$ CLQ1_2223 \f$ 
    */
    double getCLQ1_2223() const{
        return CLQ1_2223;
    }    

    /**
    * @brief Return NP coeff CLQ3_1123
    * @return \f$ CLQ3_1123 \f$ 
    */
    double getCLQ3_1123() const{
        return CLQ3_1123;
    }
    
    /**
    * @brief Return NP coeff CLQ3_2223
    * @return \f$ CLQ3_2223 \f$ 
    */
    double getCLQ3_2223() const{
        return CLQ3_2223;
    }
    
    /**
    * @brief Return NP coeff CLd_1123
    * @return \f$ CLd_1123 \f$ 
    */
    double getCLd_1123() const{
        return CLd_1123;
    }
    
    /**
    * @brief Return NP coeff CLd_2223
    * @return \f$ CLd_2223 \f$ 
    */
    double getCLd_2223() const{
        return CLd_2223;
    }
    
    /**
    * @brief Return NP coeff Ced_1123
    * @return \f$ Ced_1123 \f$ 
    */
    double getCed_1123() const{
        return Ced_1123;
    }
    
    /**
    * @brief Return NP coeff Ced_2223
    * @return \f$ Ced_2223 \f$ 
    */
    double getCed_2223() const{
        return Ced_2223;
    }
    
    /**
    * @brief Return NP coeff CHL1_11
    * @return \f$ CHL1_11 \f$ 
    */
    double getCHL1_11() const{
        return CHL1_11;
    }
    
    /**
    * @brief Return NP coeff CHL1_22
    * @return \f$ CHL1_22 \f$ 
    */
    double getCHL1_22() const{
        return CHL1_22;
    }

    /**
    * @brief Return NP coeff CHL3_11
    * @return \f$ CHL3_11 \f$ 
    */
    double getCHL3_11() const{
        return CHL3_11;
    }
    
    /**
    * @brief Return NP coeff CHL3_22
    * @return \f$ CHL3_22 \f$ 
    */
    double getCHL3_22() const{
        return CHL3_22;
    }
    
    /**
    * @brief Return NP coeff CHe_11
    * @return \f$ CHe_11 \f$ 
    */
    double getCHe_11() const{
        return CHe_11;
    }
        
    /**
    * @brief Return NP coeff CHe_22
    * @return \f$ CHe_22 \f$ 
    */
    double getCHe_22() const{
        return CHe_22;
    }
    
    /**
    * @brief Return NP coeff Ceu_1133
    * @return \f$ Ceu_1133 \f$ 
    */
    double getCeu_1133() const{
        return Ceu_1133;
    }
        
    /**
    * @brief Return NP coeff Ceu_2233
    * @return \f$ Ceu_2233 \f$ 
    */
    double getCeu_2233() const{
        return Ceu_2233;
    }
    
    /**
    * @brief Return NP coeff CLu_1133
    * @return \f$ CLu_1133 \f$ 
    */
    double getCLu_1133() const{
        return CLu_1133;
    }
        
    /**
    * @brief Return NP coeff CLu_2233
    * @return \f$ CLu_2233 \f$ 
    */
    double getCLu_2233() const{
        return CLu_2233;
    }
    
    /**
    * @brief Return NP coeff CLedq_11
    * @return \f$ CLedq_11 \f$ 
    */
    double getCLedQ_11() const{
        return CLedQ_11;
    }

    /**
    * @brief Return NP coeff CLedq_22
    * @return \f$ CLedq_22 \f$ 
    */
    double getCLedQ_22() const{
        return CLedQ_22;
    }
    
        /**
    * @brief Return NP coeff CpLedq_11
    * @return \f$ CpLedq_11 \f$ 
    */
    double getCpLedQ_11() const{
        return CpLedQ_11;
    }
    
    /**
    * @brief Return NP coeff CpLedq_22
    * @return \f$ CpLedq \f$ 
    */
    double getCpLedQ_22() const{
        return CpLedQ_22;
    }
    
    ////////////////////////////////////////////////////////////////////////
    
    /**
     * @brief A function to apply the 1st leading log corrections to the Wilson coefficients, according to the d6 SMEFT RGEs
     * @details This method runs inside the PostUpdate() method, after the definition of the internal Ci coefficients,
     * adding the log corrections according to the d6 SMEFT RGEs. It also sets the values of the corresponding anomalous dimensions.
     */
    virtual bool RGd6SMEFTlogs();

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
     * (Simplified implementation. Contribution only from @f$O_{HWB}@f$.)
     * @return the value of @f$S@f$
     */
    virtual double obliqueS() const;

    /**
     * @brief The oblique parameter \f$T\f$.
     * (Simplified implementation. Contribution only from @f$O_{HD}@f$.)
     * @return the value of @f$T@f$
     */
    virtual double obliqueT() const;

    /**
     * @brief The oblique parameter \f$U\f$.
     * @return the value of @f$U@f$
     */
    virtual double obliqueU() const;
    
    /**
     * @brief The oblique parameter \f$W\f$.
     * (Simplified implementation. Contribution only from @f$O_{2W}@f$.)
     * @return the value of @f$W@f$
     */
    virtual double obliqueW() const;

    /**
     * @brief The oblique parameter \f$Y\f$.
     * (Simplified implementation. Contribution only from @f$O_{2B}@f$.)
     * @return the value of @f$Y@f$
     */
    virtual double obliqueY() const;
    
    
    /////////////////////////////// Deviations in the experimental values of the SM input parameters /////////////////////////////////////////
    
    /**
     * @brief The relative correction to the mass of the @f$Z@f$ boson, @f$\delta M_Z/M_Z@f$, with respect to ref. point used in the SM calculation of Higgs observables.
     * @return @f$\delta M_Z/M_Z@f$
     */
    virtual double deltaMz() const;
    
    /**
     * @brief The relative correction to the mass of the @f$Z@f$ boson squared, @f$(\delta M_Z/M_Z)^2@f$, with respect to ref. point used in the SM calculation of Higgs observables.
     * @return @f$(\delta M_Z/M_Z)^2@f$
     */
    virtual double deltaMz2() const;
        
    /**
     * @brief The relative correction to the mass of the @f$H@f$ boson, @f$\delta M_H/M_H@f$, with respect to ref. point used in the SM calculation of Higgs observables.
     * @return @f$\delta M_H/M_H@f$
     */
    virtual double deltaMh() const;
    
    /**
     * @brief The relative correction to the mass of the @f$H@f$ boson squared, @f$(\delta M_H/M_H)^2@f$, with respect to ref. point used in the SM calculation of Higgs observables.
     * @return @f$(\delta M_H/M_H)^2@f$
     */
    virtual double deltaMh2() const;
    
    /**
     * @brief The relative correction to the mass of the @f$t@f$ quark, @f$\delta m_t/m_t@f$, with respect to ref. point used in the SM calculation of Higgs observables.
     * @return @f$\delta m_t/m_t@f$
     */
    virtual double deltamt() const;
    
    /**
     * @brief The relative correction to the mass of the @f$t@f$ quark squared, @f$(\delta m_t/m_t)^2@f$, with respect to ref. point used in the SM calculation of Higgs observables.
     * @return @f$(\delta m_t/m_t)^2@f$
     */
    virtual double deltamt2() const;
    
    /**
     * @brief The relative correction to the mass of the @f$b@f$ quark, @f$\delta m_b/m_b@f$, with respect to ref. point used in the SM calculation of Higgs observables.
     * @return @f$\delta m_b/m_b@f$
     */
    virtual double deltamb() const;
    
    /**
     * @brief The relative correction to the mass of the @f$b@f$ quark squared, @f$(\delta m_b/m_b)^2@f$, with respect to ref. point used in the SM calculation of Higgs observables.
     * @return @f$(\delta m_b/m_b)^2@f$
     */
    virtual double deltamb2() const;
    
    /**
     * @brief The relative correction to the mass of the @f$c@f$ quark, @f$\delta m_c/m_c@f$, with respect to ref. point used in the SM calculation of Higgs observables.
     * @return @f$\delta m_c/m_c@f$
     */
    virtual double deltamc() const;
    
    /**
     * @brief The relative correction to the mass of the @f$c@f$ quark squared, @f$(\delta m_c/m_c)^2@f$, with respect to ref. point used in the SM calculation of Higgs observables.
     * @return @f$(\delta m_c/m_c)^2@f$
     */
    virtual double deltamc2() const;
   
    /**
     * @brief The relative correction to the mass of the @f$\tau@f$ lepton, @f$\delta m_\tau/m_\tau@f$, with respect to ref. point used in the SM calculation of Higgs observables.
     * @return @f$\delta m_\tau/m_\tau@f$
     */
    virtual double deltamtau() const;
    
    /**
     * @brief The relative correction to the mass of the @f$\tau@f$ lepton squared, @f$(\delta m_\tau/m_\tau)^2@f$, with respect to ref. point used in the SM calculation of Higgs observables.
     * @return @f$(\delta m_\tau/m_\tau)^2@f$
     */
    virtual double deltamtau2() const;
    
    /**
     * @brief The relative correction to the muon decay constant, @f$\delta G_\mu/G_\mu@f$, with respect to ref. point used in the SM calculation of Higgs observables.
     * @return @f$\delta G_\mu/G_\mu@f$
     */
    virtual double deltaGmu() const;
    
    /**
     * @brief The relative correction to the muon decay constant, @f$(\delta G_\mu/G_\mu)^2@f$, with respect to ref. point used in the SM calculation of Higgs observables.
     * @return @f$(\delta G_\mu/G_\mu)^2@f$
     */
    virtual double deltaGmu2() const;
    
    /**
     * @brief The relative correction to the electromagnetic constant at the Z pole, @f$\delta \alpha(M_Z^2)/\alpha(M_Z^2)@f$, with respect to ref. point used in the SM calculation of Higgs observables.
     * @return @f$\delta \alpha(M_Z^2)/\alpha(M_Z^2)@f$
     */
    virtual double deltaaMZ() const;
    
    /**
     * @brief The relative correction to the electromagnetic constant at the Z pole, @f$(\delta \alpha(M_Z^2)/\alpha(M_Z^2))^2@f$, with respect to ref. point used in the SM calculation of Higgs observables.
     * @return @f$(\delta \alpha(M_Z^2)/\alpha(M_Z^2))^2@f$
     */
    virtual double deltaaMZ2() const;
    
    /**
     * @brief The relative correction to the electromagnetic constant at zero momentum, @f$\delta \alpha(0)/\alpha(0)@f$, with respect to ref. point used in the SM calculation of Higgs observables.
     * @return @f$\delta \alpha(0)/\alpha(0)@f$
     */
    virtual double deltaa0() const;
    
    /**
     * @brief The relative correction to the electromagnetic constant at zero momentum, @f$(\delta \alpha(0)/\alpha(0))^2@f$, with respect to ref. point used in the SM calculation of Higgs observables.
     * @return @f$(\delta \alpha(0)/\alpha(0))^2@f$
     */
    virtual double deltaa02() const;
    
    /**
     * @brief The relative correction to the strong coupling constant at the Z pole, @f$\delta \alpha_S(M_Z^2)/\alpha_S(M_Z^2)@f$, with respect to ref. point used in the SM calculation of Higgs observables.
     * @return @f$\delta \alpha_S(M_Z^2)/\alpha_S(M_Z^2)@f$
     */
    virtual double deltaaSMZ() const;
    
    /**
     * @brief The relative correction to the strong coupling constant at the Z pole, @f$(\delta \alpha_S(M_Z^2)/\alpha_S(M_Z^2))^2@f$, with respect to ref. point used in the SM calculation of Higgs observables.
     * @return @f$(\delta \alpha_S(M_Z^2)/\alpha_S(M_Z^2))^2@f$
     */
    virtual double deltaaSMZ2() const;
    
    ////////////////////////////////////////////////////////////////////////
    
    /**
     * @brief The mass of the @f$W@f$ boson, @f$M_W@f$.
     * @return @f$M_W@f$ in GeV
     */
    virtual double Mw() const;
    
    /**
     * @brief The relative NP corrections to the mass of the @f$W@f$ boson, @f$\delta M_W/M_W@f$.
     * @return @f$\delta M_W/M_W@f$
     */
    virtual double deltaMwd6() const;
    
    /**
     * @brief The relative NP corrections to the mass of the @f$W@f$ boson squared, @f$(\delta M_W/M_W)^2@f$.
     * @return @f$(\delta M_W/M_W)^2@f$
     */
    virtual double deltaMwd62() const;
    
    /**
     * @brief The new physics contribution to the decay width of the @f$W@f$ boson into a given fermion pair, @f$\delta \Gamma_Z^{f}@f$.
     * @param[in] fi a lepton or quark
     * @param[in] fj a lepton or quark
     * @return @f$\delta \Gamma_W^{ff}@f$ in GeV
     */
    virtual double deltaGamma_Wff(const Particle fi, const Particle fj) const;

    /**
     * @brief A partial decay width of the @f$W@f$ boson decay into a SM fermion pair.
     * @param[in] fi a lepton or quark
     * @param[in] fj a lepton or quark
     * @return @f$\Gamma^W_{ij}@f$
     */
    virtual double GammaW(const Particle fi, const Particle fj) const;
    
    /**
     * @brief The new physics contribution to the total decay width of the @f$W@f$ boson, @f$\delta \Gamma_W@f$.
     * @return @f$\delta \Gamma_W@f$ in GeV
     */
    virtual double deltaGamma_W() const;
    
    /**
     * @brief The total width of the @f$W@f$ boson, @f$\Gamma_W@f$.
     * @return @f$\Gamma_W@f$ in GeV
     */
    virtual double GammaW() const;
    
    /**
     * @brief The relative NP corrections to the width of the @f$W@f$ boson, @f$\delta \Gamma_W/\Gamma_W@f$.
     * @return @f$\delta \Gamma_W/\Gamma_W@f$
     */
    virtual double deltaGwd6() const;
    
    /**
     * @brief The relative NP corrections to the width of the @f$W@f$ boson squared, @f$(\delta \Gamma_W/\Gamma_W)^2@f$.
     * @return @f$(\delta \Gamma_W/\Gamma_W)^2@f$
     */
    virtual double deltaGwd62() const;
    
    /**
     * @brief The relative NP corrections to the width of the @f$Z@f$ boson, @f$\delta \Gamma_Z/\Gamma_Z@f$.
     * @return @f$\delta \Gamma_Z/\Gamma_Z@f$
     */
    virtual double deltaGzd6() const;
    
    /**
     * @brief The relative NP corrections to the width of the @f$Z@f$ boson squared, @f$(\delta \Gamma_Z/\Gamma_Z)^2@f$.
     * @return @f$(\delta \Gamma_Z/\Gamma_Z)^2@f$
     */
    virtual double deltaGzd62() const;

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
        
    /**
     * @brief The branching ratio of the @f$W@f$ boson decaying into a SM fermion pair, @f$Br(W\to f_i f_j)@f$.
     * @return @f$Br(W\to f_i f_j)@f$ in GeV
     */
    virtual double BrW(const Particle fi, const Particle fj) const;
    
    /**
     * @brief The lepton universality ratio @f$R_{W,l_i/l_j)=\Gamma(W\to l_i \nu_i)/\Gamma(W\to l_j \nu_j)@f$.
     * @return @f$R_{W,l_i/l_j)@f$ in GeV
     */
    virtual double RWlilj(const Particle li, const Particle lj) const;
    
    /**
     * @brief The lepton universality ratio @f$R_{Z,l_i/l_j)=\Gamma(Z\to l_i^+ l_i^-)/\Gamma(Z\to l_j^+ l_j^-)@f$.
     * @return @f$R_{Z,l_i/l_j)@f$ in GeV
     */
    virtual double RZlilj(const Particle li, const Particle lj) const;


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
     * @brief The new physics contribution to the Higgs self-coupling @f$ H H H@f$. Normalized to the SM value.
     * @return @f$\delta g_{HHH}/g_{HHH}^SM}@f$
     */
    virtual double deltaG_hhhRatio() const;

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
     * @brief The ratio @f$\mu_{ggHH}@f$ between the gluon-gluon fusion di-Higgs
     * production cross-section in the current model and in the Standard Model.
     * (From arXiv: 1502.00539 [hpe-ph].)
     * @param[in] sqrt_s the center-of-mass energy in TeV
     * @return @f$\mu_{ggHH}@f$
     */
    virtual double muggHH(const double sqrt_s) const;
    /**
     * @brief The ratio @f$\mu_{VBF}@f$ between the vector-boson fusion Higgs
     * production cross-section in the current model and in the Standard Model. 
     * @param[in] sqrt_s the center-of-mass energy in TeV
     * @return @f$\mu_{VBF}@f$
     */
    virtual double muVBF(const double sqrt_s) const;
    /**
     * @brief The ratio @f$\mu_{VBF+\gamma}@f$ between the vector-boson fusion Higgs
     * production cross-section in association with a hard photon in the current model
     * and in the Standard Model. 
     * @param[in] sqrt_s the center-of-mass energy in TeV
     * @return @f$\mu_{VBF+\gamma}@f$
     */
    virtual double muVBFgamma(const double sqrt_s) const;
    /**
     * @brief The ratio @f$\mu_{eeWBF}@f$ between the 
     * @f$ e^{+}e^{-}\to \nu\bar{\nu} H @f$ production
     * cross-section in the current model and in the Standard Model.
     * @param[in] sqrt_s the center-of-mass energy in TeV
     * @return @f$\mu_{eeWBF}@f$
     */
    virtual double mueeWBF(const double sqrt_s) const;
    /**
     * @brief The ratio @f$\mu_{eeWBF}@f$ between the 
     * @f$ e^{+}e^{-}\to \nu\bar{\nu} H @f$ production
     * cross-section in the current model and in the Standard Model.
     * @param[in] sqrt_s the center-of-mass energy in TeV, Pol_em and Pol_ep
     * are the polarization of electrons and positrons, respectively
     * @return @f$\mu_{eeWBF}@f$
     */
    virtual double mueeWBFPol(const double sqrt_s, const double Pol_em, const double Pol_ep) const;
    /**
     * @brief The ratio @f$\mu_{e^+e^- \to H\nu\bar{\nu}}@f$ between the 
     * @f$ e^+e^- \to H\nu\bar{\nu} @f$ associated production
     * cross-section in the current model and in the Standard Model.
     * @param[in] sqrt_s the center-of-mass energy in TeV
     * @return @f$\mu_{e^+e^- \to H\nu\bar{\nu}}@f$
     */
    virtual double mueeHvv(const double sqrt_s) const;
    /**
     * @brief The ratio @f$\mu_{e^+e^- \to H\nu\bar{\nu}}@f$ between the 
     * @f$ e^+e^- \to H\nu\bar{\nu} @f$ associated production
     * cross-section in the current model and in the Standard Model.
     * @param[in] sqrt_s the center-of-mass energy in TeV, Pol_em and Pol_ep
     * are the polarization of electrons and positrons, respectively 
     * @return @f$\mu_{e^+e^- \to H\nu\bar{\nu}}@f$
     */
    virtual double mueeHvvPol(const double sqrt_s, const double Pol_em, const double Pol_ep) const;
    /**
     * @brief The ratio @f$\mu_{eeZBF}@f$ between the 
     * @f$ e^{+}e^{-}\to e^{+}e^{-} H @f$ production
     * cross-section in the current model and in the Standard Model.
     * @param[in] sqrt_s the center-of-mass energy in TeV
     * @return @f$\mu_{eeZBF}@f$
     */
    virtual double mueeZBF(const double sqrt_s) const;
    /**
     * @brief The ratio @f$\mu_{eeZBF}@f$ between the 
     * @f$ e^{+}e^{-}\to e^{+}e^{-} H @f$ production
     * cross-section in the current model and in the Standard Model.
     * @param[in] sqrt_s the center-of-mass energy in TeV, Pol_em and Pol_ep
     * are the polarization of electrons and positrons, respectively 
     * @return @f$\mu_{eeZBF}@f$
     */
    virtual double mueeZBFPol(const double sqrt_s, const double Pol_em, const double Pol_ep) const;
    /**
     * @brief The ratio @f$\mu_{epWBF}@f$ between the 
     * @f$ e^{-} p\to \nu j H @f$ production
     * cross-section in the current model and in the Standard Model.
     * @param[in] sqrt_s the center-of-mass energy in TeV
     * @return @f$\mu_{epWBF}@f$
     */
    virtual double muepWBF(const double sqrt_s) const;
    /**
     * @brief The ratio @f$\mu_{epZBF}@f$ between the 
     * @f$ e^{-} p\to e^{-} j H @f$ production
     * cross-section in the current model and in the Standard Model.
     * @param[in] sqrt_s the center-of-mass energy in TeV
     * @return @f$\mu_{epZBF}@f$
     */
    virtual double muepZBF(const double sqrt_s) const;
    /**
     * @brief The ratio @f$\mu_{WH}@f$ between the W-Higgs associated production
     * cross-section in the current model and in the Standard Model. 
     * @param[in] sqrt_s the center-of-mass energy in TeV
     * @return @f$\mu_{WH}@f$
     */
    virtual double muWH(const double sqrt_s) const;
    /**
     * @brief The ratio @f$\mu_{WH}@f$ between the W-Higgs associated production
     * cross-section in the current model and in the Standard Model, with @f$p_{T,H}>250@f$ GeV.
     * @param[in] sqrt_s the center-of-mass energy in TeV
     * @return @f$\mu_{WH}@f$
     */
    virtual double muWHpT250(const double sqrt_s) const;
    /**
     * @brief The ratio @f$\mu_{ZH}@f$ between the Z-Higgs associated production
     * cross-section in the current model and in the Standard Model.
     * @param[in] sqrt_s the center-of-mass energy in TeV
     * @return @f$\mu_{ZH}@f$
     */
    virtual double muZH(const double sqrt_s) const;
    /**
     * @brief The ratio @f$\mu_{ZH}@f$ between the Z-Higgs associated production
     * cross-section in the current model and in the Standard Model, with @f$p_{T,H}>250@f$ GeV.
     * @param[in] sqrt_s the center-of-mass energy in TeV
     * @return @f$\mu_{ZH}@f$
     */
    virtual double muZHpT250(const double sqrt_s) const;
    /**
     * @brief The ratio @f$\mu_{eeZH}@f$ between the 
     * @f$e^{+}e^{-}\to ZH@f$ associated production
     * cross-section in the current model and in the Standard Model.
     * @param[in] sqrt_s the center-of-mass energy in TeV
     * @return @f$\mu_{eeZH}@f$
     */
    virtual double mueeZH(const double sqrt_s) const;
    /**
     * @brief The ratio @f$\mu_{eeZH, Z \to e^+ e^-, \mu^+ \mu^-}@f$ between the 
     * @f$ e^{+}e^{-}\to ZH, Z \to e^+ e^-, \mu^+ \mu^- @f$ associated production
     * cross-section in the current model and in the Standard Model.
     * @param[in] sqrt_s the center-of-mass energy in TeV
     * @return @f$\mu_{eeZH, Z \to e^+ e^-, \mu^+ \mu^-}@f$
     */
    virtual double mueeZllH(const double sqrt_s) const;
    /**
     * @brief The ratio @f$\mu_{eeZH, Z \to q \bar{q}}@f$ between the 
     * @f$ e^{+}e^{-}\to ZH, Z \to q \bar{q} @f$ associated production
     * cross-section in the current model and in the Standard Model.
     * @param[in] sqrt_s the center-of-mass energy in TeV
     * @return @f$\mu_{eeZH, Z \to q \bar{q}}@f$
     */
    virtual double mueeZqqH(const double sqrt_s) const;
    /**
     * @brief The ratio @f$\mu_{eeZH}@f$ between the 
     * @f$ e^{+}e^{-}\to ZH @f$ associated production
     * cross-section in the current model and in the Standard Model.
     * @param[in] sqrt_s the center-of-mass energy in TeV, Pol_em and Pol_ep
     * are the polarization of electrons and positrons, respectively 
     * @return @f$\mu_{eeZH}@f$
     */
    virtual double mueeZHPol(const double sqrt_s, const double Pol_em, const double Pol_ep) const;
    /**
     * @brief The ratio @f$\mu_{eeZH, Z \to e^+ e^-, \mu^+ \mu^-}@f$ between the 
     * @f$ e^{+}e^{-}\to ZH, Z \to e^+ e^-, \mu^+ \mu^- @f$ associated production
     * cross-section in the current model and in the Standard Model.
     * @param[in] sqrt_s the center-of-mass energy in TeV, Pol_em and Pol_ep
     * are the polarization of electrons and positrons, respectively 
     * @return @f$\mu_{eeZH, Z \to e^+ e^-, \mu^+ \mu^-}@f$
     */
    virtual double mueeZllHPol(const double sqrt_s, const double Pol_em, const double Pol_ep) const;    
    /**
     * @brief The ratio @f$\mu_{eeZH, Z \to q \bar{q}}@f$ between the 
     * @f$ e^{+}e^{-}\to ZH, Z \to q \bar{q} @f$ associated production
     * cross-section in the current model and in the Standard Model.
     * @param[in] sqrt_s the center-of-mass energy in TeV, Pol_em and Pol_ep
     * are the polarization of electrons and positrons, respectively 
     * @return @f$\mu_{eeZH, Z \to q \bar{q}}@f$
     */
    virtual double mueeZqqHPol(const double sqrt_s, const double Pol_em, const double Pol_ep) const;
    /**
     * @brief the angular parameter @f$a@f$ from 
     * @f$\mu_{e^+e^- \to ZH}@f$ (arXiv:1708.09079 [hep-ph]).
     * @param[in] sqrt_s the center-of-mass energy in TeV, Pol_em and Pol_ep
     * are the polarization of electrons and positrons, respectively 
     * @return @f$a_{eeZH}@f$
     */
    virtual double aPskPol(const double sqrt_s, const double Pol_em, const double Pol_ep) const;
    /**
     * @brief the angular parameter @f$b@f$ from 
     * @f$\mu_{e^+e^- \to ZH}@f$ (arXiv:1708.09079 [hep-ph]).
     * @param[in] sqrt_s the center-of-mass energy in TeV, Pol_em and Pol_ep
     * are the polarization of electrons and positrons, respectively 
     * @return @f$b_{eeZH}@f$
     */
    virtual double bPskPol(const double sqrt_s, const double Pol_em, const double Pol_ep) const;
    /**
     * @brief The ratio @f$\mu_{VH}@f$ between the WH+ZH associated production
     * cross-section in the current model and in the Standard Model.
     * @param[in] sqrt_s the center-of-mass energy in TeV
     * @return @f$\mu_{VH}@f$
     */
    virtual double muVH(const double sqrt_s) const;
    /**
     * @brief The ratio @f$\mu_{VH}@f$ between the WH+ZH associated production
     * cross-section in the current model and in the Standard Model, with @f$p_{T,H}>250@f$ GeV.
     * @param[in] sqrt_s the center-of-mass energy in TeV
     * @return @f$\mu_{VH}@f$
     */
    virtual double muVHpT250(const double sqrt_s) const;
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
     * @brief The ratio @f$\mu_{tHq}@f$ between the t-q-Higgs associated 
     * production cross-section in the current model and in the Standard Model.
     * @param[in] sqrt_s the center-of-mass energy in TeV
     * @return @f$\mu_{tHq}@f$
     */
    virtual double mutHq(const double sqrt_s) const;
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
     * @brief The ratio @f$\mu_{eettH}@f$ between the 
     * @f$ e^{+}e^{-}\to t\bar{t} H @f$ production
     * cross-section in the current model and in the Standard Model.
     * @param[in] sqrt_s the center-of-mass energy in TeV, Pol_em and Pol_ep
     * are the polarization of electrons and positrons, respectively
     * @return @f$\mu_{eettH}@f$
     */
    virtual double mueettHPol(const double sqrt_s, const double Pol_em, const double Pol_ep) const;
     /**
     * @brief The ratio @f$\mu_{\mu\mu H}@f$ between the @f$\sigma(\mu \mu \to H)}@f$
     * production cross-section in the current model and in the Standard Model.
     * @param[in] sqrt_s the center-of-mass energy in TeV
     * @return @f$\mu_{\mu\mu H}@f$
     */
    virtual double mummH(const double sqrt_s) const;
    
    ////////////HIGGS DECAY WIDTHS AND BRANCHING RATIOS/////////////
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
     * linear in the effective Lagrangian coefficients. Neglecting SM theory errors.
     * @return @f$\delta \Gamma(H)@f$/@f$\Gamma(H)_{\mathrm{SM}}@f$
     */
    virtual double deltaGammaTotalRatio1noError() const;
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
     * @brief The new physics contribution to the ratio of the @f$\Gamma(H\to gg)@f$
     * in the current model and in the Standard Model. Only terms that are
     * linear in the effective Lagrangian coefficients.
     * @return @f$\delta \Gamma(H\to gg)@f$/@f$\Gamma(H\to gg)_{\mathrm{SM}}@f$
     */
    double deltaGammaHggRatio1() const;
    /**
     * @brief The new physics contribution to the ratio of the @f$\Gamma(H\to gg)@f$
     * in the current model and in the Standard Model. (Only terms that are
     * quadratic in the effective Lagrangian coefficients.)
     * @return @f$\delta \Gamma(H\to gg)@f$/@f$\Gamma(H\to gg)_{\mathrm{SM}}@f$
     */   
    double deltaGammaHggRatio2() const;
    /**
     * @brief The ratio of the Br@f$(H\to gg)@f$ in the current model
     * and in the Standard Model.
     * @return Br@f$(H\to gg)@f$/Br@f$(H\to gg)_{\mathrm{SM}}@f$
     */
    virtual double BrHggRatio() const;
    
    /**
     * @brief The ratio of the @f$\Gamma(H\to WW)@f$ in the current model
     * and in the Standard Model.
     * @return @f$\Gamma(H\to WW)@f$/@f$\Gamma(H\to WW)_{\mathrm{SM}}@f$
     */
    double GammaHWWRatio() const;
    /**
     * @brief The new physics contribution to the ratio of the @f$\Gamma(H\to WW)@f$
     * in the current model and in the Standard Model. (Only terms that are
     * linear in the effective Lagrangian coefficients.)
     * @return @f$\delta \Gamma(H\to WW)@f$/@f$\Gamma(H\to WW)_{\mathrm{SM}}@f$
     */
    double deltaGammaHWWRatio1() const;
    /**
     * @brief The new physics contribution to the ratio of the @f$\Gamma(H\to WW)@f$
     * in the current model and in the Standard Model. (Only terms that are
     * quadratic in the effective Lagrangian coefficients.)
     * @return @f$\delta \Gamma(H\to WW)@f$/@f$\Gamma(H\to WW)_{\mathrm{SM}}@f$
     */
    double deltaGammaHWWRatio2() const;
    /**
     * @brief The ratio of the Br@f$(H\to WW)@f$ in the current model
     * and in the Standard Model.
     * @return Br@f$(H\to WW)@f$/Br@f$(H\to WW)_{\mathrm{SM}}@f$
     */
    virtual double BrHWWRatio() const;
    
    /**
     * @brief The ratio of the @f$\Gamma(H\to W l\nu)@f$ (@f$l=e,\mu @f$) in the current model
     * and in the Standard Model.
     * @return @f$\Gamma(H\to Wl\nu)@f$/@f$\Gamma(H\to Wl\nu)_{\mathrm{SM}}@f$
     */
    
    double GammaHWlvRatio() const;
    /**
     * @brief The new physics contribution to the ratio of the @f$\Gamma(H\to Wl\nu)@f$ (@f$l=e,\mu @f$)
     * in the current model and in the Standard Model. (Only terms that are
     * linear in the effective Lagrangian coefficients.)
     * @return @f$\delta \Gamma(H\to Wl\nu)@f$/@f$\Gamma(H\to Wl\nu)_{\mathrm{SM}}@f$
     */
    double deltaGammaHWlvRatio1() const;
    /**
     * @brief The new physics contribution to the ratio of the @f$\Gamma(H\to Wl\nu)@f$ (@f$l=e,\mu @f$)
     * in the current model and in the Standard Model. (Only terms that are
     * quadratic in the effective Lagrangian coefficients.)
     * @return @f$\delta \Gamma(H\to Wl\nu)@f$/@f$\Gamma(H\to Wl\nu)_{\mathrm{SM}}@f$
     */
    double deltaGammaHWlvRatio2() const;
    /**
     * @brief The ratio of the Br@f$(H\to W l\nu)@f$ (@f$l=e,\mu @f$) in the current model
     * and in the Standard Model.
     * @return Br@f$(H\to Wl\nu)@f$/Br@f$(H\to Wl\nu)_{\mathrm{SM}}@f$
     */
    virtual double BrHWlvRatio() const;
    
    /**
     * @brief The ratio of the @f$\Gamma(H\to WW^*\to l\nu l\nu)@f$ (@f$l=e,\mu @f$) in the current model
     * and in the Standard Model.
     * @return @f$\Gamma(H\to WW^*\to l\nu l\nu)@f$/@f$\Gamma(H\to WW^*\to l\nu l\nu)_{\mathrm{SM}}@f$
     */
    double GammaHWW2l2vRatio() const;
    /**
     * @brief The new physics contribution to the ratio of the @f$\Gamma(H\to WW^*\to l\nu l\nu)@f$ (@f$l=e,\mu @f$)
     * in the current model and in the Standard Model. (Only terms that are
     * linear in the effective Lagrangian coefficients.)
     * @return @f$\delta \Gamma(H\to WW^*\to l\nu l\nu)@f$/@f$\Gamma(H\to WW^*\to l\nu l\nu)_{\mathrm{SM}}@f$
     */
    double deltaGammaHWW2l2vRatio1() const;   
    /**
     * @brief The new physics contribution to the ratio of the @f$\Gamma(H\to WW^*\to l\nu l\nu)@f$ (@f$l=e,\mu @f$)
     * in the current model and in the Standard Model. (Only terms that are
     * quadratic in the effective Lagrangian coefficients.)
     * @return @f$\delta \Gamma(H\to WW^*\to l\nu l\nu)@f$/@f$\Gamma(H\to WW^*\to l\nu l\nu)_{\mathrm{SM}}@f$
     */
    double deltaGammaHWW2l2vRatio2() const;
    /**
     * @brief The ratio of the Br@f$(H\to WW^*\to l\nu l\nu)@f$ (@f$l=e,\mu @f$) in the current model
     * and in the Standard Model.
     * @return Br@f$(H\to WW^*\to l\nu l\nu)@f$/Br@f$(H\to WW^*\to l\nu l\nu)_{\mathrm{SM}}@f$
     */
    virtual double BrHWW2l2vRatio() const;
    
    /**
     * @brief The ratio of the @f$\Gamma(H\to W j j)@f$ in the current model
     * and in the Standard Model.
     * @return @f$\Gamma(H\to W j j)@f$/@f$\Gamma(H\to W j j)_{\mathrm{SM}}@f$
     */
    double GammaHWjjRatio() const;
    /**
     * @brief The new physics contribution to the ratio of the @f$\Gamma(H\to W j j)@f$
     * in the current model and in the Standard Model. (Only terms that are
     * linear in the effective Lagrangian coefficients.)
     * @return @f$\delta \Gamma(H\to W j j)@f$/@f$\Gamma(H\to W j j)_{\mathrm{SM}}@f$
     */
    double deltaGammaHWjjRatio1() const;
    /**
     * @brief The new physics contribution to the ratio of the @f$\Gamma(H\to W j j)@f$
     * in the current model and in the Standard Model. (Only terms that are
     * quadratic in the effective Lagrangian coefficients.)
     * @return @f$\delta \Gamma(H\to W j j)@f$/@f$\Gamma(H\to W j j)_{\mathrm{SM}}@f$
     */
    double deltaGammaHWjjRatio2() const;
    /**
     * @brief The ratio of the Br@f$(H\to W j j)@f$ in the current model
     * and in the Standard Model.
     * @return Br@f$(H\to W j j)@f$/Br@f$(H\to W j j)_{\mathrm{SM}}@f$
     */
    virtual double BrHWjjRatio() const;
    
    /**
     * @brief The ratio of the @f$\Gamma(H\to WW^*\to 4j)@f$ in the current model
     * and in the Standard Model.
     * @return @f$\Gamma(H\to WW^*\to 4j)@f$/@f$\Gamma(H\to WW^*\to 4j)_{\mathrm{SM}}@f$
     */
    double GammaHWW4jRatio() const;
    /**
     * @brief The new physics contribution to the ratio of the @f$\Gamma(H\to WW^*\to 4j)@f$
     * in the current model and in the Standard Model. (Only terms that are
     * linear in the effective Lagrangian coefficients.)
     * @return @f$\delta \Gamma(H\to WW^*\to 4j)@f$/@f$\Gamma(H\to WW^*\to 4j)_{\mathrm{SM}}@f$
     */
    double deltaGammaHWW4jRatio1() const;
    /**
     * @brief The new physics contribution to the ratio of the @f$\Gamma(H\to WW^*\to 4j)@f$
     * in the current model and in the Standard Model. (Only terms that are
     * quadratic in the effective Lagrangian coefficients.)
     * @return @f$\delta \Gamma(H\to WW^*\to 4j)@f$/@f$\Gamma(H\to WW^*\to 4j)_{\mathrm{SM}}@f$
     */
    double deltaGammaHWW4jRatio2() const;
    /**
     * @brief The ratio of the Br@f$(H\to WW^*\to 4j)@f$ in the current model
     * and in the Standard Model.
     * @return Br@f$(H\to WW^*\to 4j)@f$/Br@f$(H\to WW^*\to 4j)_{\mathrm{SM}}@f$
     */
    virtual double BrHWW4jRatio() const;
    
    /**
     * @brief The ratio of the @f$\Gamma(H\to W f f)@f$, with @f$f@f$ any fermion, in the current model
     * and in the Standard Model.
     * @return @f$\Gamma(H\to W f f)@f$/@f$\Gamma(H\to W f f)_{\mathrm{SM}}@f$
     */
    double GammaHWffRatio() const;
    /**
     * @brief The new physics contribution to the ratio of the @f$\Gamma(H\to W f f)@f$, with @f$f@f$ any fermion,
     * in the current model and in the Standard Model. (Only terms that are
     * linear in the effective Lagrangian coefficients.)
     * @return @f$\delta \Gamma(H\to W f f)@f$/@f$\Gamma(H\to W f f)_{\mathrm{SM}}@f$
     */
    double deltaGammaHWffRatio1() const;
    /**
     * @brief The new physics contribution to the ratio of the @f$\Gamma(H\to W f f)@f$, with @f$f@f$ any fermion,
     * in the current model and in the Standard Model. (Only terms that are
     * quadratic in the effective Lagrangian coefficients.)
     * @return @f$\delta \Gamma(H\to W f f)@f$/@f$\Gamma(H\to W f f)_{\mathrm{SM}}@f$
     */
    double deltaGammaHWffRatio2() const;
    /**
     * @brief The ratio of the Br@f$(H\to W f f)@f$, with @f$f@f$ any fermion, in the current model
     * and in the Standard Model.
     * @return Br@f$(H\to W f f)@f$/Br@f$(H\to W f f)_{\mathrm{SM}}@f$
     */
    virtual double BrHWffRatio() const;
    
    /**
     * @brief The ratio of the @f$\Gamma(H\to WW^*\to 4f)@f$, with @f$f@f$ any fermion, in the current model
     * and in the Standard Model.
     * @return @f$\Gamma(H\to WW^*\to 4f)@f$/@f$\Gamma(H\to WW^*\to 4f)_{\mathrm{SM}}@f$
     */
    double GammaHWW4fRatio() const;
    /**
     * @brief The new physics contribution to the ratio of the @f$\Gamma(H\to WW^*\to 4f)@f$, with @f$f@f$ any fermion,
     * in the current model and in the Standard Model. (Only terms that are
     * linear in the effective Lagrangian coefficients.)
     * @return @f$\delta \Gamma(H\to WW^*\to 4f)@f$/@f$\Gamma(H\to WW^*\to 4f)_{\mathrm{SM}}@f$
     */
    double deltaGammaHWW4fRatio1() const;
    /**
     * @brief The new physics contribution to the ratio of the @f$\Gamma(H\to WW^*\to 4f)@f$, with @f$f@f$ any fermion,
     * in the current model and in the Standard Model. (Only terms that are
     * quadratic in the effective Lagrangian coefficients.)
     * @return @f$\delta \Gamma(H\to WW^*\to 4f)@f$/@f$\Gamma(H\to WW^*\to 4f)_{\mathrm{SM}}@f$
     */
    double deltaGammaHWW4fRatio2() const; 
    /**
     * @brief The ratio of the Br@f$(H\to WW^*\to 4f)@f$, with @f$f@f$ any fermion, in the current model
     * and in the Standard Model.
     * @return Br@f$(H\to WW^*\to 4f)@f$/Br@f$(H\to WW^*\to 4f)_{\mathrm{SM}}@f$
     */
    virtual double BrHWW4fRatio() const;
    
    /**
     * @brief The ratio of the @f$\Gamma(H\to ZZ)@f$ in the current model
     * and in the Standard Model.
     * @return @f$\Gamma(H\to ZZ)@f$/@f$\Gamma(H\to ZZ)_{\mathrm{SM}}@f$
     */
    double GammaHZZRatio() const;
    /**
     * @brief The new physics contribution to the ratio of the @f$\Gamma(H\to ZZ)@f$
     * in the current model and in the Standard Model. (Only terms that are
     * linear in the effective Lagrangian coefficients.)
     * @return @f$\delta \Gamma(H\to ZZ)@f$/@f$\Gamma(H\to ZZ)_{\mathrm{SM}}@f$
     */
    double deltaGammaHZZRatio1() const;    
    /**
     * @brief The new physics contribution to the ratio of the @f$\Gamma(H\to ZZ)@f$
     * in the current model and in the Standard Model. (Only terms that are
     * quadratic in the effective Lagrangian coefficients.)
     * @return @f$\delta \Gamma(H\to ZZ)@f$/@f$\Gamma(H\to ZZ)_{\mathrm{SM}}@f$
     */
    double deltaGammaHZZRatio2() const;
    /**
     * @brief The ratio of the Br@f$(H\to ZZ)@f$ in the current model
     * and in the Standard Model.
     * @return Br@f$(H\to ZZ)@f$/Br@f$(H\to ZZ)_{\mathrm{SM}}@f$
     */
    virtual double BrHZZRatio() const;
    
    /**
     * @brief The ratio of the @f$\Gamma(H\to Zll)@f$ (@f$l=e,\mu @f$) in the current model
     * and in the Standard Model.
     * @return @f$\Gamma(H\to Zll)@f$/@f$\Gamma(H\to Zll)_{\mathrm{SM}}@f$
     */
    double GammaHZllRatio() const;
    /**
     * @brief The new physics contribution to the ratio of the @f$\Gamma(H\to Zll)@f$ (@f$l=e,\mu @f$)
     * in the current model and in the Standard Model. (Only terms that are
     * linear in the effective Lagrangian coefficients.)
     * @return @f$\delta \Gamma(H\to Zll)@f$/@f$\Gamma(H\to Zll)_{\mathrm{SM}}@f$
     */
    double deltaGammaHZllRatio1() const;
    /**
     * @brief The new physics contribution to the ratio of the @f$\Gamma(H\to Zll)@f$ (@f$l=e,\mu @f$)
     * in the current model and in the Standard Model. (Only terms that are
     * quadratic in the effective Lagrangian coefficients.)
     * @return @f$\delta \Gamma(H\to Zll)@f$/@f$\Gamma(H\to Zll)_{\mathrm{SM}}@f$
     */
    double deltaGammaHZllRatio2() const;
    /**
     * @brief The ratio of the Br@f$(H\to Zll)@f$ (@f$l=e,\mu @f$) in the current model
     * and in the Standard Model.
     * @return Br@f$(H\to Zll)@f$/Br@f$(H\to Zll)_{\mathrm{SM}}@f$
     */
    virtual double BrHZllRatio() const;
    
    /**
     * @brief The ratio of the @f$\Gamma(H\to Zee)@f$ in the current model
     * and in the Standard Model.
     * @return @f$\Gamma(H\to Zee)@f$/@f$\Gamma(H\to Zee)_{\mathrm{SM}}@f$
     */
    double GammaHZeeRatio() const;
    /**
     * @brief The new physics contribution to the ratio of the @f$\Gamma(H\to Zee)@f$
     * in the current model and in the Standard Model. (Only terms that are
     * linear in the effective Lagrangian coefficients.)
     * @return @f$\delta \Gamma(H\to Zee)@f$/@f$\Gamma(H\to Zee)_{\mathrm{SM}}@f$
     */    
    double deltaGammaHZeeRatio1() const;
    /**
     * @brief The new physics contribution to the ratio of the @f$\Gamma(H\to Zee)@f$
     * in the current model and in the Standard Model. (Only terms that are
     * quadratic in the effective Lagrangian coefficients.)
     * @return @f$\delta \Gamma(H\to Zee)@f$/@f$\Gamma(H\to Zee)_{\mathrm{SM}}@f$
     */
    double deltaGammaHZeeRatio2() const;
    
    /**
     * @brief The ratio of the @f$\Gamma(H\to Z\mu\mu)@f$ in the current model
     * and in the Standard Model.
     * @return @f$\Gamma(H\to Z\mu\mu)@f$/@f$\Gamma(H\to Z\mu\mu)_{\mathrm{SM}}@f$
     */
    double GammaHZmumuRatio() const;    
    /**
     * @brief The new physics contribution to the ratio of the @f$\Gamma(H\to Z\mu\mu)@f$
     * in the current model and in the Standard Model. (Only terms that are
     * linear in the effective Lagrangian coefficients.)
     * @return @f$\delta \Gamma(H\to Z\mu\mu)@f$/@f$\Gamma(H\to Z\mu\mu)_{\mathrm{SM}}@f$
     */
    double deltaGammaHZmumuRatio1() const;
    /**
     * @brief The new physics contribution to the ratio of the @f$\Gamma(H\to Z\mu\mu)@f$
     * in the current model and in the Standard Model. (Only terms that are
     * quadratic in the effective Lagrangian coefficients.)
     * @return @f$\delta \Gamma(H\to Z\mu\mu)@f$/@f$\Gamma(H\to Z\mu\mu)_{\mathrm{SM}}@f$
     */
    double deltaGammaHZmumuRatio2() const;
    
    /**
     * @brief The ratio of the @f$\Gamma(H\to ZZ* \to 4l)@f$ (@f$l=e,\mu @f$) in the current model
     * and in the Standard Model.
     * @return @f$\Gamma(H\to ZZ* \to 4l)@f$/@f$\Gamma(H\to ZZ* \to 4l)_{\mathrm{SM}}@f$
     */
    double GammaHZZ4lRatio() const;    
    /**
     * @brief The new physics contribution to the ratio of the @f$\Gamma(H\to ZZ* \to 4l)@f$ (@f$l=e,\mu @f$)
     * in the current model and in the Standard Model. (Only terms that are
     * linear in the effective Lagrangian coefficients.)
     * @return @f$\delta \Gamma(H\to ZZ* \to 4l)@f$/@f$\Gamma(H\to ZZ* \to 4l)_{\mathrm{SM}}@f$
     */
    double deltaGammaHZZ4lRatio1() const;  
    /**
     * @brief The new physics contribution to the ratio of the @f$\Gamma(H\to ZZ* \to 4l)@f$ (@f$l=e,\mu @f$)
     * in the current model and in the Standard Model. (Only terms that are
     * quadratic in the effective Lagrangian coefficients.)
     * @return @f$\delta \Gamma(H\to ZZ* \to 4l)@f$/@f$\Gamma(H\to ZZ* \to 4l)_{\mathrm{SM}}@f$
     */
    double deltaGammaHZZ4lRatio2() const;
    /**
     * @brief The ratio of the Br@f$(H\to ZZ* \to 4l)@f$ (@f$l=e,\mu @f$) in the current model
     * and in the Standard Model.
     * @return Br@f$(H\to ZZ* \to 4l)@f$/Br@f$(H\to ZZ* \to 4l)_{\mathrm{SM}}@f$
     */
    virtual double BrHZZ4lRatio() const;
    
    /**
     * @brief The ratio of the @f$\Gamma(H\to ZZ* \to 4e)@f$ in the current model
     * and in the Standard Model.
     * @return @f$\Gamma(H\to ZZ* \to 4e)@f$/@f$\Gamma(H\to ZZ* \to 4e)_{\mathrm{SM}}@f$
     */
    double GammaHZZ4eRatio() const;
    /**
     * @brief The new physics contribution to the ratio of the @f$\Gamma(H\to ZZ* \to 4e)@f$
     * in the current model and in the Standard Model. (Only terms that are
     * linear in the effective Lagrangian coefficients.)
     * @return @f$\delta \Gamma(H\to ZZ* \to 4e)@f$/@f$\Gamma(H\to ZZ* \to 4e)_{\mathrm{SM}}@f$
     */
    double deltaGammaHZZ4eRatio1() const;
    /**
     * @brief The new physics contribution to the ratio of the @f$\Gamma(H\to ZZ* \to 4e)@f$
     * in the current model and in the Standard Model. (Only terms that are
     * quadratic in the effective Lagrangian coefficients.)
     * @return @f$\delta \Gamma(H\to ZZ* \to 4e)@f$/@f$\Gamma(H\to ZZ* \to 4e)_{\mathrm{SM}}@f$
     */
    double deltaGammaHZZ4eRatio2() const;
    /**
     * @brief The ratio of the Br@f$(H\to ZZ* \to 4e)@f$ in the current model
     * and in the Standard Model.
     * @return Br@f$(H\to ZZ* \to 4e)@f$/Br@f$(H\to ZZ* \to 4e)_{\mathrm{SM}}@f$
     */
    virtual double BrHZZ4eRatio() const;
    
    /**
     * @brief The ratio of the @f$\Gamma(H\to ZZ* \to 2e2\mu)@f$ in the current model
     * and in the Standard Model.
     * @return @f$\Gamma(H\to ZZ* \to 2e2\mu)@f$/@f$\Gamma(H\to ZZ* \to 2e2\mu)_{\mathrm{SM}}@f$
     */
    double GammaHZZ2e2muRatio() const;
    /**
     * @brief The new physics contribution to the ratio of the @f$\Gamma(H\to ZZ* \to 2e2\mu)@f$
     * in the current model and in the Standard Model. (Only terms that are
     * linear in the effective Lagrangian coefficients.)
     * @return @f$\delta \Gamma(H\to ZZ* \to 2e2\mu)@f$/@f$\Gamma(H\to ZZ* \to 2e2\mu)_{\mathrm{SM}}@f$
     */
    double deltaGammaHZZ2e2muRatio1() const;
    /**
     * @brief The new physics contribution to the ratio of the @f$\Gamma(H\to ZZ* \to 2e2\mu)@f$
     * in the current model and in the Standard Model. (Only terms that are
     * quadratic in the effective Lagrangian coefficients.)
     * @return @f$\delta \Gamma(H\to ZZ* \to 2e2\mu)@f$/@f$\Gamma(H\to ZZ* \to 2e2\mu)_{\mathrm{SM}}@f$
     */
    double deltaGammaHZZ2e2muRatio2() const;
    /**
     * @brief The ratio of the Br@f$(H\to ZZ* \to 2e 2\mu)@f$ in the current model
     * and in the Standard Model.
     * @return Br@f$(H\to ZZ* \to 2e 2\mu)@f$/Br@f$(H\to ZZ* \to 2e 2\mu)_{\mathrm{SM}}@f$
     */
    virtual double BrHZZ2e2muRatio() const;
       
    /**
     * @brief The ratio of the @f$\Gamma(H\to ZZ* \to 4\mu)@f$ in the current model
     * and in the Standard Model.
     * @return @f$\Gamma(H\to ZZ* \to 4\mu)@f$/@f$\Gamma(H\to ZZ* \to 4\mu)_{\mathrm{SM}}@f$
     */
    double GammaHZZ4muRatio() const;    
    /**
     * @brief The new physics contribution to the ratio of the @f$\Gamma(H\to ZZ* \to 4\mu)@f$
     * in the current model and in the Standard Model. (Only terms that are
     * linear in the effective Lagrangian coefficients.)
     * @return @f$\delta \Gamma(H\to ZZ* \to 4\mu)@f$/@f$\Gamma(H\to ZZ* \to 4\mu)_{\mathrm{SM}}@f$
     */
    double deltaGammaHZZ4muRatio1() const;
    /**
     * @brief The new physics contribution to the ratio of the @f$\Gamma(H\to ZZ* \to 4\mu)@f$
     * in the current model and in the Standard Model. (Only terms that are
     * quadratic in the effective Lagrangian coefficients.)
     * @return @f$\delta \Gamma(H\to ZZ* \to 4\mu)@f$/@f$\Gamma(H\to ZZ* \to 4\mu)_{\mathrm{SM}}@f$
     */
    double deltaGammaHZZ4muRatio2() const;
    /**
     * @brief The ratio of the Br@f$(H\to ZZ* \to 4\mu)@f$ in the current model
     * and in the Standard Model.
     * @return Br@f$(H\to ZZ* \to 4\mu)@f$/Br@f$(H\to ZZ* \to 4\mu)_{\mathrm{SM}}@f$
     */
    virtual double BrHZZ4muRatio() const;
    
    /**
     * @brief The ratio of the @f$\Gamma(H\to Z\nu\nu)@f$ in the current model
     * and in the Standard Model.
     * @return @f$\Gamma(H\to Z\nu\nu)@f$/@f$\Gamma(H\to Z\nu\nu)_{\mathrm{SM}}@f$
     */
    double GammaHZvvRatio() const;    
    /**
     * @brief The new physics contribution to the ratio of the @f$\Gamma(H\to Z\nu\nu)@f$
     * in the current model and in the Standard Model. (Only terms that are
     * linear in the effective Lagrangian coefficients.)
     * @return @f$\delta \Gamma(H\to Z\nu\nu)@f$/@f$\Gamma(H\to Z\nu\nu)_{\mathrm{SM}}@f$
     */
    double deltaGammaHZvvRatio1() const;
    /**
     * @brief The new physics contribution to the ratio of the @f$\Gamma(H\to Z\nu\nu)@f$
     * in the current model and in the Standard Model. (Only terms that are
     * quadratic in the effective Lagrangian coefficients.)
     * @return @f$\delta \Gamma(H\to Z\nu\nu)@f$/@f$\Gamma(H\to Z\nu\nu)_{\mathrm{SM}}@f$
     */
    double deltaGammaHZvvRatio2() const;
    /**
     * @brief The ratio of the Br@f$(H\to Z\nu\nu)@f$ in the current model
     * and in the Standard Model.
     * @return Br@f$(H\to Z\nu\nu)@f$/Br@f$(H\to Z\nu\nu)_{\mathrm{SM}}@f$
     */
    virtual double BrHZvvRatio() const;
    
    /**
     * @brief The ratio of the @f$\Gamma(H\to ZZ* \to 4\nu)@f$ in the current model
     * and in the Standard Model.
     * @return @f$\Gamma(H\to ZZ* \to 4\nu)@f$/@f$\Gamma(H\to ZZ* \to 4\nu)_{\mathrm{SM}}@f$
     */
    double GammaHZZ4vRatio() const;    
    /**
     * @brief The new physics contribution to the ratio of the @f$\Gamma(H\to ZZ* \to 4\nu)@f$
     * in the current model and in the Standard Model. (Only terms that are
     * linear in the effective Lagrangian coefficients.)
     * @return @f$\delta \Gamma(H\to ZZ* \to 4\nu)@f$/@f$\Gamma(H\to ZZ* \to 4\nu)_{\mathrm{SM}}@f$
     */
    double deltaGammaHZZ4vRatio1() const;
    /**
     * @brief The new physics contribution to the ratio of the @f$\Gamma(H\to ZZ* \to 4\nu)@f$
     * in the current model and in the Standard Model. (Only terms that are
     * quadratic in the effective Lagrangian coefficients.)
     * @return @f$\delta \Gamma(H\to ZZ* \to 4\nu)@f$/@f$\Gamma(H\to ZZ* \to 4\nu)_{\mathrm{SM}}@f$
     */
    double deltaGammaHZZ4vRatio2() const;
    /**
     * @brief The ratio of the Br@f$(H\to ZZ* \to 4\nu)@f$ in the current model
     * and in the Standard Model.
     * @return Br@f$(H\to ZZ* \to 4\nu)@f$/Br@f$(H\to ZZ* \to 4\nu)_{\mathrm{SM}}@f$
     */
    virtual double BrHZZ4vRatio() const;
    
    /**
     * @brief The ratio of the @f$\Gamma(H\to Zu u)@f$ (@f$u=u,c @f$) in the current model
     * and in the Standard Model.
     * @return @f$\Gamma(H\to Zu u)@f$/@f$\Gamma(H\to Zu u)_{\mathrm{SM}}@f$
     */
    double GammaHZuuRatio() const;    
    /**
     * @brief The new physics contribution to the ratio of the @f$\Gamma(H\to Z u u)@f$ (@f$u=u,c @f$)
     * in the current model and in the Standard Model. (Only terms that are
     * linear in the effective Lagrangian coefficients.)
     * @return @f$\delta \Gamma(H\to Z u u)@f$/@f$\Gamma(H\to Z u u)_{\mathrm{SM}}@f$
     */
    double deltaGammaHZuuRatio1() const;
    /**
     * @brief The new physics contribution to the ratio of the @f$\Gamma(H\to Z u u)@f$ (@f$u=u,c @f$)
     * in the current model and in the Standard Model. (Only terms that are
     * quadratic in the effective Lagrangian coefficients.)
     * @return @f$\delta \Gamma(H\to Z u u)@f$/@f$\Gamma(H\to Z u u)_{\mathrm{SM}}@f$
     */
    double deltaGammaHZuuRatio2() const;
    /**
     * @brief The ratio of the Br@f$(H\to Z u u)@f$ (@f$u=u,c @f$) in the current model
     * and in the Standard Model.
     * @return Br@f$(H\to Z u u)@f$/Br@f$(H\to Z u u)_{\mathrm{SM}}@f$
     */
    virtual double BrHZuuRatio() const;
    
    /**
     * @brief The ratio of the @f$\Gamma(H\to ZZ* \to 4 u)@f$ (@f$u=u,c @f$) in the current model
     * and in the Standard Model.
     * @return @f$\Gamma(H\to ZZ* \to 4u)@f$/@f$\Gamma(H\to ZZ* \to 4u)_{\mathrm{SM}}@f$
     */
    double GammaHZZ4uRatio() const
    {
        return 1.0;
    };
    /**
     * @brief The new physics contribution to the ratio of the @f$\Gamma(H\to ZZ* \to 4 u)@f$ (@f$u=u,c @f$)
     * in the current model and in the Standard Model. (Only terms that are
     * linear in the effective Lagrangian coefficients.)
     * @return @f$\delta \Gamma(H\to ZZ* \to 4 u)@f$/@f$\Gamma(H\to ZZ* \to 4 u)_{\mathrm{SM}}@f$
     */
    double deltaGammaHZZ4uRatio1() const
    {
        return 0.0;
    };
    /**
     * @brief The new physics contribution to the ratio of the @f$\Gamma(H\to ZZ* \to 4 u)@f$ (@f$u=u,c @f$)
     * in the current model and in the Standard Model. (Only terms that are
     * quadratic in the effective Lagrangian coefficients.)
     * @return @f$\delta \Gamma(H\to ZZ* \to 4 u)@f$/@f$\Gamma(H\to ZZ* \to 4 u)_{\mathrm{SM}}@f$
     */
    double deltaGammaHZZ4uRatio2() const
    {
        return 0.0;
    };
    /**
     * @brief The ratio of the Br@f$(H\to ZZ* \to 4 u)@f$ (@f$u=u,c @f$) in the current model
     * and in the Standard Model.
     * @return Br@f$(H\to ZZ* \to 4 u)@f$/Br@f$(H\to ZZ* \to 4 u)_{\mathrm{SM}}@f$
     */
    virtual double BrHZZ4uRatio() const;
    
    /**
     * @brief The ratio of the @f$\Gamma(H\to Zd d)@f$ (@f$d=d,s,b @f$) in the current model
     * and in the Standard Model.
     * @return @f$\Gamma(H\to Zd d)@f$/@f$\Gamma(H\to Zd d)_{\mathrm{SM}}@f$
     */
    double GammaHZddRatio() const;
    /**
     * @brief The new physics contribution to the ratio of the @f$\Gamma(H\to Z d d)@f$ (@f$d=d,s,b @f$)
     * in the current model and in the Standard Model. (Only terms that are
     * linear in the effective Lagrangian coefficients.)
     * @return @f$\delta \Gamma(H\to Z d d)@f$/@f$\Gamma(H\to Z d d)_{\mathrm{SM}}@f$
     */
    double deltaGammaHZddRatio1() const;
    /**
     * @brief The new physics contribution to the ratio of the @f$\Gamma(H\to Z d d)@f$ (@f$d=d,s,b @f$)
     * in the current model and in the Standard Model. (Only terms that are
     * quadratic in the effective Lagrangian coefficients.)
     * @return @f$\delta \Gamma(H\to Z d d)@f$/@f$\Gamma(H\to Z d d)_{\mathrm{SM}}@f$
     */
    double deltaGammaHZddRatio2() const;
    /**
     * @brief The ratio of the Br@f$(H\to Z d d)@f$ (@f$d=d,s,b @f$) in the current model
     * and in the Standard Model.
     * @return Br@f$(H\to Z d d)@f$/Br@f$(H\to Z d d)_{\mathrm{SM}}@f$
     */
    virtual double BrHZddRatio() const;
    
    /**
     * @brief The ratio of the @f$\Gamma(H\to ZZ* \to 4 d)@f$ (@f$d=d,s,b @f$) in the current model
     * and in the Standard Model.
     * @return @f$\Gamma(H\to ZZ* \to 4 d)@f$/@f$\Gamma(H\to ZZ* \to 4 d)_{\mathrm{SM}}@f$
     */
    double GammaHZZ4dRatio() const
    {
        return 1.0;
    };
    /**
     * @brief The new physics contribution to the ratio of the @f$\Gamma(H\to ZZ* \to 4 d)@f$ (@f$d=d,s,b @f$)
     * in the current model and in the Standard Model. (Only terms that are
     * linear in the effective Lagrangian coefficients.)
     * @return @f$\delta \Gamma(H\to ZZ* \to 4 d)@f$/@f$\Gamma(H\to ZZ* \to 4 d)_{\mathrm{SM}}@f$
     */
    double deltaGammaHZZ4dRatio1() const
    {
        return 0.0;
    };
    /**
     * @brief The new physics contribution to the ratio of the @f$\Gamma(H\to ZZ* \to 4 d)@f$ (@f$d=d,s,b @f$)
     * in the current model and in the Standard Model. (Only terms that are
     * quadratic in the effective Lagrangian coefficients.)
     * @return @f$\delta \Gamma(H\to ZZ* \to 4 d)@f$/@f$\Gamma(H\to ZZ* \to 4 d)_{\mathrm{SM}}@f$
     */
    double deltaGammaHZZ4dRatio2() const
    {
        return 0.0;
    };
    /**
     * @brief The ratio of the Br@f$(H\to ZZ* \to 4 d)@f$ (@f$d=d,s,b @f$) in the current model
     * and in the Standard Model.
     * @return Br@f$(H\to ZZ* \to 4 d)@f$/Br@f$(H\to ZZ* \to 4 d)_{\mathrm{SM}}@f$
     */
    virtual double BrHZZ4dRatio() const;
    
    /**
     * @brief The ratio of the @f$\Gamma(H\to Zff)@f$, with @f$f@f$ any fermion, in the current model
     * and in the Standard Model.
     * @return @f$\Gamma(H\to Zff)@f$/@f$\Gamma(H\to Zff)_{\mathrm{SM}}@f$
     */
    double GammaHZffRatio() const;
    /**
     * @brief The new physics contribution to the ratio of the @f$\Gamma(H\to Z ff)@f$, with @f$f@f$ any fermion,
     * in the current model and in the Standard Model. (Only terms that are
     * linear in the effective Lagrangian coefficients.)
     * @return @f$\delta \Gamma(H\to Z f f)@f$/@f$\Gamma(H\to Z f f)_{\mathrm{SM}}@f$
     */
    double deltaGammaHZffRatio1() const;
    /**
     * @brief The new physics contribution to the ratio of the @f$\Gamma(H\to Z ff)@f$, with @f$f@f$ any fermion,
     * in the current model and in the Standard Model. (Only terms that are
     * quadratic in the effective Lagrangian coefficients.)
     * @return @f$\delta \Gamma(H\to Z f f)@f$/@f$\Gamma(H\to Z f f)_{\mathrm{SM}}@f$
     */
    double deltaGammaHZffRatio2() const;
    /**
     * @brief The ratio of the Br@f$(H\to Zff)@f$, with @f$f@f$ any fermion, in the current model
     * and in the Standard Model.
     * @return Br@f$(H\to Zff)@f$/Br@f$(H\to Zff)_{\mathrm{SM}}@f$
     */
    virtual double BrHZffRatio() const;
    
    /**
     * @brief The ratio of the @f$\Gamma(H\to ZZ* \to 4f)@f$, with @f$f@f$ any fermion, in the current model
     * and in the Standard Model.
     * @return @f$\Gamma(H\to ZZ* \to 4f)@f$/@f$\Gamma(H\to ZZ* \to 4f)_{\mathrm{SM}}@f$
     */
    double GammaHZZ4fRatio() const;
    /**
     * @brief The new physics contribution to the ratio of the @f$\Gamma(H\to ZZ* \to 4f)@f$, with @f$f@f$ any fermion,
     * in the current model and in the Standard Model. (Only terms that are
     * linear in the effective Lagrangian coefficients.)
     * @return @f$\delta \Gamma(H\to ZZ* \to 4f)@f$/@f$\Gamma(H\to ZZ* \to 4f)_{\mathrm{SM}}@f$
     */
    double deltaGammaHZZ4fRatio1() const;
    /**
     * @brief The new physics contribution to the ratio of the @f$\Gamma(H\to ZZ* \to 4f)@f$, with @f$f@f$ any fermion,
     * in the current model and in the Standard Model. (Only terms that are
     * quadratic in the effective Lagrangian coefficients.)
     * @return @f$\delta \Gamma(H\to ZZ* \to 4f)@f$/@f$\Gamma(H\to ZZ* \to 4f)_{\mathrm{SM}}@f$
     */
    double deltaGammaHZZ4fRatio2() const; 
    /**
     * @brief The ratio of the Br@f$(H\to ZZ* \to 4f)@f$, with @f$f@f$ any fermion, in the current model
     * and in the Standard Model.
     * @return Br@f$(H\to ZZ* \to 4f)@f$/Br@f$(H\to ZZ* \to 4f)_{\mathrm{SM}}@f$
     */
    virtual double BrHZZ4fRatio() const;
    
    /**
     * @brief The ratio of the @f$\Gamma(H\to Z\gamma)@f$ in the current model
     * and in the Standard Model.
     * @return @f$\Gamma(H\to Z\gamma)@f$/@f$\Gamma(H\to Z\gamma)_{\mathrm{SM}}@f$
     */
    double GammaHZgaRatio() const;
    /**
     * @brief The new physics contribution to the ratio of the @f$\Gamma(H\to Z\gamma)@f$
     * in the current model and in the Standard Model. (Only terms that are
     * linear in the effective Lagrangian coefficients.)
     * @return @f$\delta \Gamma(H\to Z\gamma)@f$/@f$\Gamma(H\to Z\gamma)_{\mathrm{SM}}@f$
     */
    double deltaGammaHZgaRatio1() const;
    /**
     * @brief The new physics contribution to the ratio of the @f$\Gamma(H\to Z\gamma)@f$
     * in the current model and in the Standard Model. (Only terms that are
     * quadratic in the effective Lagrangian coefficients.)
     * @return @f$\delta \Gamma(H\to Z\gamma)@f$/@f$\Gamma(H\to Z\gamma)_{\mathrm{SM}}@f$
     */
    double deltaGammaHZgaRatio2() const;
    /**
     * @brief The ratio of the Br@f$(H\to Z\gamma)@f$ in the current model
     * and in the Standard Model.
     * @return Br@f$(H\to Z\gamma)@f$/Br@f$(H\to Z\gamma)_{\mathrm{SM}}@f$
     */
    virtual double BrHZgaRatio() const;
    /**
     * @brief The ratio of the Br@f$(H\to Z\gamma\to ll\gamma)@f$ (@f$l=e,\mu @f$) in the current model
     * and in the Standard Model.
     * @return Br@f$(H\to Z\gamma\to ll\gamma)@f$/Br@f$(H\to Z\gamma\to ll\gamma)_{\mathrm{SM}}@f$
     */
    virtual double BrHZgallRatio() const;
    /**
     * @brief The ratio of the Br@f$(H\to Z\gamma\to ee\gamma)@f$ in the current model
     * and in the Standard Model.
     * @return Br@f$(H\to Z\gamma\to ee\gamma)@f$/Br@f$(H\to Z\gamma\to ee\gamma)_{\mathrm{SM}}@f$
     */
    virtual double BrHZgaeeRatio() const;
    /**
     * @brief The ratio of the Br@f$(H\to Z\gamma\to \mu\mu\gamma)@f$ in the current model
     * and in the Standard Model.
     * @return Br@f$(H\to Z\gamma\to \mu\mu\gamma)@f$/Br@f$(H\to Z\gamma\to \mu\mu\gamma)_{\mathrm{SM}}@f$
     */
    virtual double BrHZgamumuRatio() const;
    
    /**
     * @brief The ratio of the @f$\Gamma(H\to \gamma\gamma)@f$ in the current model
     * and in the Standard Model.
     * @return @f$\Gamma(H\to \gamma\gamma)@f$/@f$\Gamma(H\to \gamma\gamma)_{\mathrm{SM}}@f$
     */
    double GammaHgagaRatio() const;
    /**
     * @brief The new physics contribution to the ratio of the @f$\Gamma(H\to \gamma\gamma)@f$
     * in the current model and in the Standard Model. (Only terms that are
     * linear in the effective Lagrangian coefficients.)
     * @return @f$\delta \Gamma(H\to \gamma\gamma)@f$/@f$\Gamma(H\to \gamma\gamma)_{\mathrm{SM}}@f$
     */
    double deltaGammaHgagaRatio1() const;
    /**
     * @brief The new physics contribution to the ratio of the @f$\Gamma(H\to \gamma\gamma)@f$
     * in the current model and in the Standard Model. (Only terms that are
     * quadratic in the effective Lagrangian coefficients.)
     * @return @f$\delta \Gamma(H\to \gamma\gamma)@f$/@f$\Gamma(H\to \gamma\gamma)_{\mathrm{SM}}@f$
     */
    double deltaGammaHgagaRatio2() const;
    /**
     * @brief The ratio of the Br@f$(H\to \gamma\gamma)@f$ in the current model
     * and in the Standard Model.
     * @return Br@f$(H\to \gamma\gamma)@f$/Br@f$(H\to \gamma\gamma)_{\mathrm{SM}}@f$
     */
    virtual double BrHgagaRatio() const;
    
    /**
     * @brief The ratio of the @f$\Gamma(H\to \mu\mu)@f$ in the current model
     * and in the Standard Model.
     * @return @f$\Gamma(H\to \mu\mu)@f$/@f$\Gamma(H\to \mu\mu)_{\mathrm{SM}}@f$
     */
    double GammaHmumuRatio() const;
    /**
     * @brief The new physics contribution to the ratio of the @f$\Gamma(H\to \mu\mu)@f$
     * in the current model and in the Standard Model. (Only terms that are
     * linear in the effective Lagrangian coefficients.)
     * @return @f$\delta \Gamma(H\to \mu\mu)@f$/@f$\Gamma(H\to \mu\mu)_{\mathrm{SM}}@f$
     */
    double deltaGammaHmumuRatio1() const;
    /**
     * @brief The new physics contribution to the ratio of the @f$\Gamma(H\to \mu\mu)@f$
     * in the current model and in the Standard Model. (Only terms that are
     * quadratic in the effective Lagrangian coefficients.)
     * @return @f$\delta \Gamma(H\to \mu\mu)@f$/@f$\Gamma(H\to \mu\mu)_{\mathrm{SM}}@f$
     */
    double deltaGammaHmumuRatio2() const;
    /**
     * @brief The ratio of the Br@f$(H\to \mu^+\mu^-)@f$ in the current model
     * and in the Standard Model.
     * @return Br@f$(H\to \mu^+\mu^-)@f$/Br@f$(H\to \mu^+\mu^-)_{\mathrm{SM}}@f$
     */
    virtual double BrHmumuRatio() const;
    
    /**
     * @brief The ratio of the @f$\Gamma(H\to \tau\tau)@f$ in the current model
     * and in the Standard Model.
     * @return @f$\Gamma(H\to \tau\tau)@f$/@f$\Gamma(H\to \tau\tau)_{\mathrm{SM}}@f$
     */
    double GammaHtautauRatio() const;
    /**
     * @brief The new physics contribution to the ratio of the @f$\Gamma(H\to \tau\tau)@f$
     * in the current model and in the Standard Model. (Only terms that are
     * linear in the effective Lagrangian coefficients.)
     * @return @f$\delta \Gamma(H\to \tau\tau)@f$/@f$\Gamma(H\to \tau\tau)_{\mathrm{SM}}@f$
     */
    double deltaGammaHtautauRatio1() const;
    /**
     * @brief The new physics contribution to the ratio of the @f$\Gamma(H\to \tau\tau)@f$
     * in the current model and in the Standard Model. (Only terms that are
     * quadratic in the effective Lagrangian coefficients.)
     * @return @f$\delta \Gamma(H\to \tau\tau)@f$/@f$\Gamma(H\to \tau\tau)_{\mathrm{SM}}@f$
     */
    double deltaGammaHtautauRatio2() const;
    /**
     * @brief The ratio of the Br@f$(H\to \tau^+\tau^-)@f$ in the current model
     * and in the Standard Model.
     * @return Br@f$(H\to \tau^+\tau^-)@f$/Br@f$(H\to \tau^+\tau^-)_{\mathrm{SM}}@f$
     */
    virtual double BrHtautauRatio() const;
    
    /**
     * @brief The ratio of the @f$\Gamma(H\to cc)@f$ in the current model
     * and in the Standard Model.
     * @return @f$\Gamma(H\to cc)@f$/@f$\Gamma(H\to cc)_{\mathrm{SM}}@f$
     */
    double GammaHccRatio() const;
    /**
     * @brief The new physics contribution to the ratio of the @f$\Gamma(H\to cc)@f$
     * in the current model and in the Standard Model. (Only terms that are
     * linear in the effective Lagrangian coefficients.)
     * @return @f$\delta \Gamma(H\to cc)@f$/@f$\Gamma(H\to cc)_{\mathrm{SM}}@f$
     */
    double deltaGammaHccRatio1() const;
    /**
     * @brief The new physics contribution to the ratio of the @f$\Gamma(H\to cc)@f$
     * in the current model and in the Standard Model. (Only terms that are
     * quadratic in the effective Lagrangian coefficients.)
     * @return @f$\delta \Gamma(H\to cc)@f$/@f$\Gamma(H\to cc)_{\mathrm{SM}}@f$
     */
    double deltaGammaHccRatio2() const;
    /**
     * @brief The ratio of the Br@f$(H\to c\bar{c})@f$ in the current model
     * and in the Standard Model.
     * @return Br@f$(H\to c\bar{c})@f$/Br@f$(H\to c\bar{c})_{\mathrm{SM}}@f$
     */
    virtual double BrHccRatio() const;
    
    /**
     * @brief The ratio of the @f$\Gamma(H\to bb)@f$ in the current model
     * and in the Standard Model.
     * @return @f$\Gamma(H\to bb)@f$/@f$\Gamma(H\to bb)_{\mathrm{SM}}@f$
     */
    double GammaHbbRatio() const;
    /**
     * @brief The new physics contribution to the ratio of the @f$\Gamma(H\to bb)@f$
     * in the current model and in the Standard Model. (Only terms that are
     * linear in the effective Lagrangian coefficients.)
     * @return @f$\delta \Gamma(H\to bb)@f$/@f$\Gamma(H\to bb)_{\mathrm{SM}}@f$
     */
    double deltaGammaHbbRatio1() const;
    /**
     * @brief The new physics contribution to the ratio of the @f$\Gamma(H\to bb)@f$
     * in the current model and in the Standard Model. (Only terms that are
     * quadratic in the effective Lagrangian coefficients.)
     * @return @f$\delta \Gamma(H\to bb)@f$/@f$\Gamma(H\to bb)_{\mathrm{SM}}@f$
     */
    double deltaGammaHbbRatio2() const;
    /**
     * @brief The ratio of the Br@f$(H\to b\bar{b})@f$ in the current model
     * and in the Standard Model.
     * @return Br@f$(H\to b\bar{b})@f$/Br@f$(H\to b\bar{b})_{\mathrm{SM}}@f$
     */
    virtual double BrHbbRatio() const;
    
    // General 4 fermion final states
    
    /**
     * @brief The ratio of the @f$\Gamma(H\to 2L2L)@f$ (@f$L=e,\mu,\tau@f$) in the current model
     * and in the Standard Model.
     * @return @f$\Gamma(H\to 2L2L)@f$/@f$\Gamma(H\to 2L2L)_{\mathrm{SM}}@f$
     */
    double GammaH2L2LRatio() const;    
    /**
     * @brief The new physics contribution to the ratio of the @f$\Gamma(H\to 2L2L)@f$ (@f$L=e,\mu,\tau@f$)
     * in the current model and in the Standard Model. (Only terms that are
     * linear in the effective Lagrangian coefficients.)
     * @return @f$\delta \Gamma(H\to 2L2L)@f$/@f$\Gamma(H\to 2L2L)_{\mathrm{SM}}@f$
     */
    double deltaGammaH2L2LRatio1() const;
    /**
     * @brief The new physics contribution to the ratio of the @f$\Gamma(H\to 2L2L)@f$ (@f$L=e,\mu,\tau@f$)
     * in the current model and in the Standard Model. (Only terms that are
     * quadratic in the effective Lagrangian coefficients.)
     * @return @f$\delta \Gamma(H\to 2L2L)@f$/@f$\Gamma(H\to 2L2L)_{\mathrm{SM}}@f$
     */
    double deltaGammaH2L2LRatio2() const;
    /**
     * @brief The ratio of the Br@f$(H\to 2L2L)@f$ (@f$L=e,\mu,\tau@f$) in the current model
     * and in the Standard Model.
     * @return Br@f$(H\to 2L2L)@f$/Br@f$(H\to 2L2L)_{\mathrm{SM}}@f$
     */
    virtual double BrH2L2LRatio() const;
    
    /**
     * @brief The ratio of the @f$\Gamma(H\to 2L2L)@f$ (@f$L=e,\mu@f$) in the current model
     * and in the Standard Model.
     * @return @f$\Gamma(H\to 2L2L)@f$/@f$\Gamma(H\to 2L2L)_{\mathrm{SM}}@f$
     */
    double GammaH2L2L2Ratio() const;    
    /**
     * @brief The new physics contribution to the ratio of the @f$\Gamma(H\to 2L2L)@f$ (@f$L=e,\mu@f$)
     * in the current model and in the Standard Model. (Only terms that are
     * linear in the effective Lagrangian coefficients.)
     * @return @f$\delta \Gamma(H\to 2L2L)@f$/@f$\Gamma(H\to 2L2L)_{\mathrm{SM}}@f$
     */
    double deltaGammaH2L2L2Ratio1() const;
    /**
     * @brief The new physics contribution to the ratio of the @f$\Gamma(H\to 2L2L)@f$ (@f$L=e,\mu@f$)
     * in the current model and in the Standard Model. (Only terms that are
     * quadratic in the effective Lagrangian coefficients.)
     * @return @f$\delta \Gamma(H\to 2L2L)@f$/@f$\Gamma(H\to 2L2L)_{\mathrm{SM}}@f$
     */
    double deltaGammaH2L2L2Ratio2() const;
    /**
     * @brief The ratio of the Br@f$(H\to 2L2L)@f$ (@f$L=e,\mu@f$) in the current model
     * and in the Standard Model.
     * @return Br@f$(H\to 2L2L)@f$/Br@f$(H\to 2L2L)_{\mathrm{SM}}@f$
     */
    virtual double BrH2L2L2Ratio() const;
    
    /**
     * @brief The ratio of the @f$\Gamma(H\to 2v2v)@f$ in the current model
     * and in the Standard Model.
     * @return @f$\Gamma(H\to 2v2v)@f$/@f$\Gamma(H\to 2v2v)_{\mathrm{SM}}@f$
     */
    double GammaH2v2vRatio() const;    
    /**
     * @brief The new physics contribution to the ratio of the @f$\Gamma(H\to 2v2v)@f$
     * in the current model and in the Standard Model. (Only terms that are
     * linear in the effective Lagrangian coefficients.)
     * @return @f$\delta \Gamma(H\to 2v2v)@f$/@f$\Gamma(H\to 2v2v)_{\mathrm{SM}}@f$
     */
    double deltaGammaH2v2vRatio1() const;
    /**
     * @brief The new physics contribution to the ratio of the @f$\Gamma(H\to 2v2v)@f$
     * in the current model and in the Standard Model. (Only terms that are
     * quadratic in the effective Lagrangian coefficients.)
     * @return @f$\delta \Gamma(H\to 2v2v)@f$/@f$\Gamma(H\to 2v2v)_{\mathrm{SM}}@f$
     */
    double deltaGammaH2v2vRatio2() const;
    /**
     * @brief The ratio of the Br@f$(H\to 2v2v)@f$ in the current model
     * and in the Standard Model.
     * @return Br@f$(H\to 2v2v)@f$/Br@f$(H\to 2v2v)_{\mathrm{SM}}@f$
     */
    virtual double BrH2v2vRatio() const;
    
    /**
     * @brief The ratio of the @f$\Gamma(H\to 2L2v)@f$ (@f$L=e,\mu,\tau@f$) in the current model
     * and in the Standard Model.
     * @return @f$\Gamma(H\to 2L2v)@f$/@f$\Gamma(H\to 2L2v)_{\mathrm{SM}}@f$
     */
    double GammaH2L2vRatio() const;    
    /**
     * @brief The new physics contribution to the ratio of the @f$\Gamma(H\to 2L2v)@f$ (@f$L=e,\mu,\tau@f$)
     * in the current model and in the Standard Model. (Only terms that are
     * linear in the effective Lagrangian coefficients.)
     * @return @f$\delta \Gamma(H\to 2L2v)@f$/@f$\Gamma(H\to 2L2v)_{\mathrm{SM}}@f$
     */
    double deltaGammaH2L2vRatio1() const;
    /**
     * @brief The new physics contribution to the ratio of the @f$\Gamma(H\to 2L2v)@f$ (@f$L=e,\mu,\tau@f$)
     * in the current model and in the Standard Model. (Only terms that are
     * quadratic in the effective Lagrangian coefficients.)
     * @return @f$\delta \Gamma(H\to 2L2v)@f$/@f$\Gamma(H\to 2L2v)_{\mathrm{SM}}@f$
     */
    double deltaGammaH2L2vRatio2() const;
    /**
     * @brief The ratio of the Br@f$(H\to 2L2v)@f$ (@f$L=e,\mu,\tau@f$) in the current model
     * and in the Standard Model.
     * @return Br@f$(H\to 2L2v)@f$/Br@f$(H\to 2L2v)_{\mathrm{SM}}@f$
     */
    virtual double BrH2L2vRatio() const;
    
    
    /**
     * @brief The ratio of the @f$\Gamma(H\to 2L2v)@f$ (@f$L=e,\mu@f$) in the current model
     * and in the Standard Model.
     * @return @f$\Gamma(H\to 2L2v)@f$/@f$\Gamma(H\to 2L2v)_{\mathrm{SM}}@f$
     */
    double GammaH2L2v2Ratio() const;    
    /**
     * @brief The new physics contribution to the ratio of the @f$\Gamma(H\to 2L2v)@f$ (@f$L=e,\mu@f$)
     * in the current model and in the Standard Model. (Only terms that are
     * linear in the effective Lagrangian coefficients.)
     * @return @f$\delta \Gamma(H\to 2L2v)@f$/@f$\Gamma(H\to 2L2v)_{\mathrm{SM}}@f$
     */
    double deltaGammaH2L2v2Ratio1() const;
    /**
     * @brief The new physics contribution to the ratio of the @f$\Gamma(H\to 2L2v)@f$ (@f$L=e,\mu@f$)
     * in the current model and in the Standard Model. (Only terms that are
     * quadratic in the effective Lagrangian coefficients.)
     * @return @f$\delta \Gamma(H\to 2L2v)@f$/@f$\Gamma(H\to 2L2v)_{\mathrm{SM}}@f$
     */
    double deltaGammaH2L2v2Ratio2() const;
    /**
     * @brief The ratio of the Br@f$(H\to 2L2v)@f$ (@f$L=e,\mu@f$) in the current model
     * and in the Standard Model.
     * @return Br@f$(H\to 2L2v)@f$/Br@f$(H\to 2L2v)_{\mathrm{SM}}@f$
     */
    virtual double BrH2L2v2Ratio() const;
    
    
    /**
     * @brief The ratio of the @f$\Gamma(H\to 2u2u)@f$ in the current model
     * and in the Standard Model.
     * @return @f$\Gamma(H\to 2u2u)@f$/@f$\Gamma(H\to 2u2u)_{\mathrm{SM}}@f$
     */
    double GammaH2u2uRatio() const;    
    /**
     * @brief The new physics contribution to the ratio of the @f$\Gamma(H\to 2u2u)@f$
     * in the current model and in the Standard Model. (Only terms that are
     * linear in the effective Lagrangian coefficients.)
     * @return @f$\delta \Gamma(H\to 2u2u)@f$/@f$\Gamma(H\to 2u2u)_{\mathrm{SM}}@f$
     */
    double deltaGammaH2u2uRatio1() const;
    /**
     * @brief The new physics contribution to the ratio of the @f$\Gamma(H\to 2u2u)@f$
     * in the current model and in the Standard Model. (Only terms that are
     * quadratic in the effective Lagrangian coefficients.)
     * @return @f$\delta \Gamma(H\to 2u2u)@f$/@f$\Gamma(H\to 2u2u)_{\mathrm{SM}}@f$
     */
    double deltaGammaH2u2uRatio2() const;
    /**
     * @brief The ratio of the Br@f$(H\to 2u2u)@f$ in the current model
     * and in the Standard Model.
     * @return Br@f$(H\to 2u2u)@f$/Br@f$(H\to 2u2u)_{\mathrm{SM}}@f$
     */
    virtual double BrH2u2uRatio() const;
    
    /**
     * @brief The ratio of the @f$\Gamma(H\to 2d2d)@f$ in the current model
     * and in the Standard Model.
     * @return @f$\Gamma(H\to 2d2d)@f$/@f$\Gamma(H\to 2d2d)_{\mathrm{SM}}@f$
     */
    double GammaH2d2dRatio() const;    
    /**
     * @brief The new physics contribution to the ratio of the @f$\Gamma(H\to 2d2d)@f$
     * in the current model and in the Standard Model. (Only terms that are
     * linear in the effective Lagrangian coefficients.)
     * @return @f$\delta \Gamma(H\to 2d2d)@f$/@f$\Gamma(H\to 2d2d)_{\mathrm{SM}}@f$
     */
    double deltaGammaH2d2dRatio1() const;
    /**
     * @brief The new physics contribution to the ratio of the @f$\Gamma(H\to 2d2d)@f$
     * in the current model and in the Standard Model. (Only terms that are
     * quadratic in the effective Lagrangian coefficients.)
     * @return @f$\delta \Gamma(H\to 2d2d)@f$/@f$\Gamma(H\to 2d2d)_{\mathrm{SM}}@f$
     */
    double deltaGammaH2d2dRatio2() const;
    /**
     * @brief The ratio of the Br@f$(H\to 2d2d)@f$ in the current model
     * and in the Standard Model.
     * @return Br@f$(H\to 2d2d)@f$/Br@f$(H\to 2d2d)_{\mathrm{SM}}@f$
     */
    virtual double BrH2d2dRatio() const;
    
    /**
     * @brief The ratio of the @f$\Gamma(H\to 2u2d)@f$ in the current model
     * and in the Standard Model.
     * @return @f$\Gamma(H\to 2u2d)@f$/@f$\Gamma(H\to 2u2d)_{\mathrm{SM}}@f$
     */
    double GammaH2u2dRatio() const;    
    /**
     * @brief The new physics contribution to the ratio of the @f$\Gamma(H\to 2u2d)@f$
     * in the current model and in the Standard Model. (Only terms that are
     * linear in the effective Lagrangian coefficients.)
     * @return @f$\delta \Gamma(H\to 2u2d)@f$/@f$\Gamma(H\to 2u2d)_{\mathrm{SM}}@f$
     */
    double deltaGammaH2u2dRatio1() const;
    /**
     * @brief The new physics contribution to the ratio of the @f$\Gamma(H\to 2u2d)@f$
     * in the current model and in the Standard Model. (Only terms that are
     * quadratic in the effective Lagrangian coefficients.)
     * @return @f$\delta \Gamma(H\to 2u2d)@f$/@f$\Gamma(H\to 2u2d)_{\mathrm{SM}}@f$
     */
    double deltaGammaH2u2dRatio2() const;
    /**
     * @brief The ratio of the Br@f$(H\to 2u2d)@f$ in the current model
     * and in the Standard Model.
     * @return Br@f$(H\to 2u2d)@f$/Br@f$(H\to 2u2d)_{\mathrm{SM}}@f$
     */
    virtual double BrH2u2dRatio() const;
    
    /**
     * @brief The ratio of the @f$\Gamma(H\to 2L2u)@f$ (@f$L=e,\mu,\tau@f$) in the current model
     * and in the Standard Model.
     * @return @f$\Gamma(H\to 2L2u)@f$/@f$\Gamma(H\to 2L2u)_{\mathrm{SM}}@f$
     */
    double GammaH2L2uRatio() const;    
    /**
     * @brief The new physics contribution to the ratio of the @f$\Gamma(H\to 2L2u)@f$ (@f$L=e,\mu,\tau@f$)
     * in the current model and in the Standard Model. (Only terms that are
     * linear in the effective Lagrangian coefficients.)
     * @return @f$\delta \Gamma(H\to 2L2u)@f$/@f$\Gamma(H\to 2L2u)_{\mathrm{SM}}@f$
     */
    double deltaGammaH2L2uRatio1() const;
    /**
     * @brief The new physics contribution to the ratio of the @f$\Gamma(H\to 2L2u)@f$ (@f$L=e,\mu,\tau@f$)
     * in the current model and in the Standard Model. (Only terms that are
     * quadratic in the effective Lagrangian coefficients.)
     * @return @f$\delta \Gamma(H\to 2L2u)@f$/@f$\Gamma(H\to 2L2u)_{\mathrm{SM}}@f$
     */
    double deltaGammaH2L2uRatio2() const;
    /**
     * @brief The ratio of the Br@f$(H\to 2L2u)@f$ (@f$L=e,\mu,\tau@f$) in the current model
     * and in the Standard Model.
     * @return Br@f$(H\to 2L2u)@f$/Br@f$(H\to 2L2u)_{\mathrm{SM}}@f$
     */
    virtual double BrH2L2uRatio() const;
    
    /**
     * @brief The ratio of the @f$\Gamma(H\to 2L2d)@f$ (@f$L=e,\mu,\tau@f$) in the current model
     * and in the Standard Model.
     * @return @f$\Gamma(H\to 2L2d)@f$/@f$\Gamma(H\to 2L2d)_{\mathrm{SM}}@f$
     */
    double GammaH2L2dRatio() const;    
    /**
     * @brief The new physics contribution to the ratio of the @f$\Gamma(H\to 2L2d)@f$ (@f$L=e,\mu,\tau@f$)
     * in the current model and in the Standard Model. (Only terms that are
     * linear in the effective Lagrangian coefficients.)
     * @return @f$\delta \Gamma(H\to 2L2d)@f$/@f$\Gamma(H\to 2L2d)_{\mathrm{SM}}@f$
     */
    double deltaGammaH2L2dRatio1() const;
    /**
     * @brief The new physics contribution to the ratio of the @f$\Gamma(H\to 2L2d)@f$ (@f$L=e,\mu,\tau@f$)
     * in the current model and in the Standard Model. (Only terms that are
     * quadratic in the effective Lagrangian coefficients.)
     * @return @f$\delta \Gamma(H\to 2L2d)@f$/@f$\Gamma(H\to 2L2d)_{\mathrm{SM}}@f$
     */
    double deltaGammaH2L2dRatio2() const;
    /**
     * @brief The ratio of the Br@f$(H\to 2L2d)@f$ (@f$L=e,\mu,\tau@f$) in the current model
     * and in the Standard Model.
     * @return Br@f$(H\to 2L2d)@f$/Br@f$(H\to 2L2d)_{\mathrm{SM}}@f$
     */
    virtual double BrH2L2dRatio() const;
    
    /**
     * @brief The ratio of the @f$\Gamma(H\to 2v2u)@f$ in the current model
     * and in the Standard Model.
     * @return @f$\Gamma(H\to 2v2u)@f$/@f$\Gamma(H\to 2v2u)_{\mathrm{SM}}@f$
     */
    double GammaH2v2uRatio() const;    
    /**
     * @brief The new physics contribution to the ratio of the @f$\Gamma(H\to 2v2u)@f$
     * in the current model and in the Standard Model. (Only terms that are
     * linear in the effective Lagrangian coefficients.)
     * @return @f$\delta \Gamma(H\to 2v2u)@f$/@f$\Gamma(H\to 2v2u)_{\mathrm{SM}}@f$
     */
    double deltaGammaH2v2uRatio1() const;
    /**
     * @brief The new physics contribution to the ratio of the @f$\Gamma(H\to 2v2u)@f$
     * in the current model and in the Standard Model. (Only terms that are
     * quadratic in the effective Lagrangian coefficients.)
     * @return @f$\delta \Gamma(H\to 2v2u)@f$/@f$\Gamma(H\to 2v2u)_{\mathrm{SM}}@f$
     */
    double deltaGammaH2v2uRatio2() const;
    /**
     * @brief The ratio of the Br@f$(H\to 2v2u)@f$ in the current model
     * and in the Standard Model.
     * @return Br@f$(H\to 2v2u)@f$/Br@f$(H\to 2v2u)_{\mathrm{SM}}@f$
     */
    virtual double BrH2v2uRatio() const;
    
    /**
     * @brief The ratio of the @f$\Gamma(H\to 2v2d)@f$ in the current model
     * and in the Standard Model.
     * @return @f$\Gamma(H\to 2v2d)@f$/@f$\Gamma(H\to 2v2d)_{\mathrm{SM}}@f$
     */
    double GammaH2v2dRatio() const;    
    /**
     * @brief The new physics contribution to the ratio of the @f$\Gamma(H\to 2v2d)@f$
     * in the current model and in the Standard Model. (Only terms that are
     * linear in the effective Lagrangian coefficients.)
     * @return @f$\delta \Gamma(H\to 2v2d)@f$/@f$\Gamma(H\to 2v2d)_{\mathrm{SM}}@f$
     */
    double deltaGammaH2v2dRatio1() const;
    /**
     * @brief The new physics contribution to the ratio of the @f$\Gamma(H\to 2v2d)@f$
     * in the current model and in the Standard Model. (Only terms that are
     * quadratic in the effective Lagrangian coefficients.)
     * @return @f$\delta \Gamma(H\to 2v2d)@f$/@f$\Gamma(H\to 2v2d)_{\mathrm{SM}}@f$
     */
    double deltaGammaH2v2dRatio2() const;
    /**
     * @brief The ratio of the Br@f$(H\to 2v2d)@f$ in the current model
     * and in the Standard Model.
     * @return Br@f$(H\to 2v2d)@f$/Br@f$(H\to 2v2d)_{\mathrm{SM}}@f$
     */
    virtual double BrH2v2dRatio() const;
    
    /**
     * @brief The ratio of the @f$\Gamma(H\to 4L)@f$ (@f$L=e,\mu,\tau@f$) in the current model
     * and in the Standard Model.
     * @return @f$\Gamma(H\to 4L)@f$/@f$\Gamma(H\to 4L)_{\mathrm{SM}}@f$
     */
    double GammaH4LRatio() const;    
    /**
     * @brief The new physics contribution to the ratio of the @f$\Gamma(H\to 4L)@f$ (@f$L=e,\mu,\tau@f$)
     * in the current model and in the Standard Model. (Only terms that are
     * linear in the effective Lagrangian coefficients.)
     * @return @f$\delta \Gamma(H\to 4L)@f$/@f$\Gamma(H\to 4L)_{\mathrm{SM}}@f$
     */
    double deltaGammaH4LRatio1() const;
    /**
     * @brief The new physics contribution to the ratio of the @f$\Gamma(H\to 4L)@f$ (@f$L=e,\mu,\tau@f$)
     * in the current model and in the Standard Model. (Only terms that are
     * quadratic in the effective Lagrangian coefficients.)
     * @return @f$\delta \Gamma(H\to 4L)@f$/@f$\Gamma(H\to 4L)_{\mathrm{SM}}@f$
     */
    double deltaGammaH4LRatio2() const;
    /**
     * @brief The ratio of the Br@f$(H\to 4L)@f$ (@f$L=e,\mu,\tau@f$) in the current model
     * and in the Standard Model.
     * @return Br@f$(H\to 4L)@f$/Br@f$(H\to 4L)_{\mathrm{SM}}@f$
     */
    virtual double BrH4LRatio() const;
    
    
    /**
     * @brief The ratio of the @f$\Gamma(H\to 4L)@f$ (@f$L=e,\mu@f$) in the current model
     * and in the Standard Model.
     * @return @f$\Gamma(H\to 4L)@f$/@f$\Gamma(H\to 4L)_{\mathrm{SM}}@f$
     */
    double GammaH4L2Ratio() const;    
    /**
     * @brief The new physics contribution to the ratio of the @f$\Gamma(H\to 4L)@f$ (@f$L=e,\mu@f$)
     * in the current model and in the Standard Model. (Only terms that are
     * linear in the effective Lagrangian coefficients.)
     * @return @f$\delta \Gamma(H\to 4L)@f$/@f$\Gamma(H\to 4L)_{\mathrm{SM}}@f$
     */
    double deltaGammaH4L2Ratio1() const;
    /**
     * @brief The new physics contribution to the ratio of the @f$\Gamma(H\to 4L)@f$ (@f$L=e,\mu@f$)
     * in the current model and in the Standard Model. (Only terms that are
     * quadratic in the effective Lagrangian coefficients.)
     * @return @f$\delta \Gamma(H\to 4L)@f$/@f$\Gamma(H\to 4L)_{\mathrm{SM}}@f$
     */
    double deltaGammaH4L2Ratio2() const;
    /**
     * @brief The ratio of the Br@f$(H\to 4L)@f$ (@f$L=e,\mu@f$) in the current model
     * and in the Standard Model.
     * @return Br@f$(H\to 4L)@f$/Br@f$(H\to 4L)_{\mathrm{SM}}@f$
     */
    virtual double BrH4L2Ratio() const;
    
    
    /**
     * @brief The ratio of the @f$\Gamma(H\to 4v)@f$ in the current model
     * and in the Standard Model.
     * @return @f$\Gamma(H\to 4v)@f$/@f$\Gamma(H\to 4v)_{\mathrm{SM}}@f$
     */
    double GammaH4vRatio() const;    
    /**
     * @brief The new physics contribution to the ratio of the @f$\Gamma(H\to 4v)@f$
     * in the current model and in the Standard Model. (Only terms that are
     * linear in the effective Lagrangian coefficients.)
     * @return @f$\delta \Gamma(H\to 4v)@f$/@f$\Gamma(H\to 4v)_{\mathrm{SM}}@f$
     */
    double deltaGammaH4vRatio1() const;
    /**
     * @brief The new physics contribution to the ratio of the @f$\Gamma(H\to 4v)@f$
     * in the current model and in the Standard Model. (Only terms that are
     * quadratic in the effective Lagrangian coefficients.)
     * @return @f$\delta \Gamma(H\to 4v)@f$/@f$\Gamma(H\to 4v)_{\mathrm{SM}}@f$
     */
    double deltaGammaH4vRatio2() const;
    /**
     * @brief The ratio of the Br@f$(H\to 4v)@f$ in the current model
     * and in the Standard Model.
     * @return Br@f$(H\to 4v)@f$/Br@f$(H\to 4v)_{\mathrm{SM}}@f$
     */
    virtual double BrH4vRatio() const;
    
    /**
     * @brief The ratio of the @f$\Gamma(H\to 4u)@f$ in the current model
     * and in the Standard Model.
     * @return @f$\Gamma(H\to 4u)@f$/@f$\Gamma(H\to 4u)_{\mathrm{SM}}@f$
     */
    double GammaH4uRatio() const;    
    /**
     * @brief The new physics contribution to the ratio of the @f$\Gamma(H\to 4u)@f$
     * in the current model and in the Standard Model. (Only terms that are
     * linear in the effective Lagrangian coefficients.)
     * @return @f$\delta \Gamma(H\to 4u)@f$/@f$\Gamma(H\to 4u)_{\mathrm{SM}}@f$
     */
    double deltaGammaH4uRatio1() const;
    /**
     * @brief The new physics contribution to the ratio of the @f$\Gamma(H\to 4u)@f$
     * in the current model and in the Standard Model. (Only terms that are
     * quadratic in the effective Lagrangian coefficients.)
     * @return @f$\delta \Gamma(H\to 4u)@f$/@f$\Gamma(H\to 4u)_{\mathrm{SM}}@f$
     */
    double deltaGammaH4uRatio2() const;
    /**
     * @brief The ratio of the Br@f$(H\to 4u)@f$ in the current model
     * and in the Standard Model.
     * @return Br@f$(H\to 4u)@f$/Br@f$(H\to 4u)_{\mathrm{SM}}@f$
     */
    virtual double BrH4uRatio() const;
    
    /**
     * @brief The ratio of the @f$\Gamma(H\to 4d)@f$ in the current model
     * and in the Standard Model.
     * @return @f$\Gamma(H\to 4d)@f$/@f$\Gamma(H\to 4d)_{\mathrm{SM}}@f$
     */
    double GammaH4dRatio() const;    
    /**
     * @brief The new physics contribution to the ratio of the @f$\Gamma(H\to 4d)@f$
     * in the current model and in the Standard Model. (Only terms that are
     * linear in the effective Lagrangian coefficients.)
     * @return @f$\delta \Gamma(H\to 4d)@f$/@f$\Gamma(H\to 4d)_{\mathrm{SM}}@f$
     */
    double deltaGammaH4dRatio1() const;
    /**
     * @brief The new physics contribution to the ratio of the @f$\Gamma(H\to 4d)@f$
     * in the current model and in the Standard Model. (Only terms that are
     * quadratic in the effective Lagrangian coefficients.)
     * @return @f$\delta \Gamma(H\to 4d)@f$/@f$\Gamma(H\to 4d)_{\mathrm{SM}}@f$
     */
    double deltaGammaH4dRatio2() const;
    /**
     * @brief The ratio of the Br@f$(H\to 4d)@f$ in the current model
     * and in the Standard Model.
     * @return Br@f$(H\to 4d)@f$/Br@f$(H\to 4d)_{\mathrm{SM}}@f$
     */
    virtual double BrH4dRatio() const;
    
    /**
     * @brief The ratio of the @f$\Gamma(H\to LvvL)@f$ (@f$L=e,\mu,\tau@f$) in the current model
     * and in the Standard Model.
     * @return @f$\Gamma(H\to LvvL)@f$/@f$\Gamma(H\to LvvL)_{\mathrm{SM}}@f$
     */
    double GammaHLvvLRatio() const;    
    /**
     * @brief The new physics contribution to the ratio of the @f$\Gamma(H\to LvvL)@f$ (@f$L=e,\mu,\tau@f$)
     * in the current model and in the Standard Model. (Only terms that are
     * linear in the effective Lagrangian coefficients.)
     * @return @f$\delta \Gamma(H\to LvvL)@f$/@f$\Gamma(H\to LvvL)_{\mathrm{SM}}@f$
     */
    double deltaGammaHLvvLRatio1() const;
    /**
     * @brief The new physics contribution to the ratio of the @f$\Gamma(H\to LvvL)@f$ (@f$L=e,\mu,\tau@f$)
     * in the current model and in the Standard Model. (Only terms that are
     * quadratic in the effective Lagrangian coefficients.)
     * @return @f$\delta \Gamma(H\to LvvL)@f$/@f$\Gamma(H\to LvvL)_{\mathrm{SM}}@f$
     */
    double deltaGammaHLvvLRatio2() const;
    /**
     * @brief The ratio of the Br@f$(H\to LvvL)@f$ (@f$L=e,\mu,\tau@f$) in the current model
     * and in the Standard Model.
     * @return Br@f$(H\to LvvL)@f$/Br@f$(H\to LvvL)_{\mathrm{SM}}@f$
     */
    virtual double BrHLvvLRatio() const;
    
    /**
     * @brief The ratio of the @f$\Gamma(H\to LvvL)@f$ (@f$L=e,\mu@f$) in the current model
     * and in the Standard Model.
     * @return @f$\Gamma(H\to LvvL)@f$/@f$\Gamma(H\to LvvL)_{\mathrm{SM}}@f$
     */
    double GammaHLvvL2Ratio() const;    
    /**
     * @brief The new physics contribution to the ratio of the @f$\Gamma(H\to LvvL)@f$ (@f$L=e,\mu@f$)
     * in the current model and in the Standard Model. (Only terms that are
     * linear in the effective Lagrangian coefficients.)
     * @return @f$\delta \Gamma(H\to LvvL)@f$/@f$\Gamma(H\to LvvL)_{\mathrm{SM}}@f$
     */
    double deltaGammaHLvvL2Ratio1() const;
    /**
     * @brief The new physics contribution to the ratio of the @f$\Gamma(H\to LvvL)@f$ (@f$L=e,\mu@f$)
     * in the current model and in the Standard Model. (Only terms that are
     * quadratic in the effective Lagrangian coefficients.)
     * @return @f$\delta \Gamma(H\to LvvL)@f$/@f$\Gamma(H\to LvvL)_{\mathrm{SM}}@f$
     */
    double deltaGammaHLvvL2Ratio2() const;
    /**
     * @brief The ratio of the Br@f$(H\to LvvL)@f$ (@f$L=e,\mu@f$) in the current model
     * and in the Standard Model.
     * @return Br@f$(H\to LvvL)@f$/Br@f$(H\to LvvL)_{\mathrm{SM}}@f$
     */
    virtual double BrHLvvL2Ratio() const;
    
    
    /**
     * @brief The ratio of the @f$\Gamma(H\to uddu)@f$ in the current model
     * and in the Standard Model.
     * @return @f$\Gamma(H\to uddu)@f$/@f$\Gamma(H\to uddu)_{\mathrm{SM}}@f$
     */
    double GammaHudduRatio() const;    
    /**
     * @brief The new physics contribution to the ratio of the @f$\Gamma(H\to uddu)@f$
     * in the current model and in the Standard Model. (Only terms that are
     * linear in the effective Lagrangian coefficients.)
     * @return @f$\delta \Gamma(H\to uddu)@f$/@f$\Gamma(H\to uddu)_{\mathrm{SM}}@f$
     */
    double deltaGammaHudduRatio1() const;
    /**
     * @brief The new physics contribution to the ratio of the @f$\Gamma(H\to uddu)@f$
     * in the current model and in the Standard Model. (Only terms that are
     * quadratic in the effective Lagrangian coefficients.)
     * @return @f$\delta \Gamma(H\to uddu)@f$/@f$\Gamma(H\to uddu)_{\mathrm{SM}}@f$
     */
    double deltaGammaHudduRatio2() const;
    /**
     * @brief The ratio of the Br@f$(H\to uddu)@f$ in the current model
     * and in the Standard Model.
     * @return Br@f$(H\to uddu)@f$/Br@f$(H\to uddu)_{\mathrm{SM}}@f$
     */
    virtual double BrHudduRatio() const;
    
    /**
     * @brief The ratio of the @f$\Gamma(H\to Lvud)@f$ (@f$L=e,\mu,\tau@f$) in the current model
     * and in the Standard Model.
     * @return @f$\Gamma(H\to Lvud)@f$/@f$\Gamma(H\to Lvud)_{\mathrm{SM}}@f$
     */
    double GammaHLvudRatio() const;    
    /**
     * @brief The new physics contribution to the ratio of the @f$\Gamma(H\to Lvud)@f$ (@f$L=e,\mu,\tau@f$)
     * in the current model and in the Standard Model. (Only terms that are
     * linear in the effective Lagrangian coefficients.)
     * @return @f$\delta \Gamma(H\to Lvud)@f$/@f$\Gamma(H\to Lvud)_{\mathrm{SM}}@f$
     */
    double deltaGammaHLvudRatio1() const;
    /**
     * @brief The new physics contribution to the ratio of the @f$\Gamma(H\to Lvud)@f$ (@f$L=e,\mu,\tau@f$)
     * in the current model and in the Standard Model. (Only terms that are
     * quadratic in the effective Lagrangian coefficients.)
     * @return @f$\delta \Gamma(H\to Lvud)@f$/@f$\Gamma(H\to Lvud)_{\mathrm{SM}}@f$
     */
    double deltaGammaHLvudRatio2() const;
    /**
     * @brief The ratio of the Br@f$(H\to Lvud)@f$ (@f$L=e,\mu,\tau@f$) in the current model
     * and in the Standard Model.
     * @return Br@f$(H\to Lvud)@f$/Br@f$(H\to Lvud)_{\mathrm{SM}}@f$
     */
    virtual double BrHLvudRatio() const;
    
    /**
     * @brief The ratio of the @f$\Gamma(H\to 2ud)@f$ in the current model
     * and in the Standard Model.
     * @return @f$\Gamma(H\to 2ud)@f$/@f$\Gamma(H\to 2ud)_{\mathrm{SM}}@f$
     */
    double GammaH2udRatio() const;    
    /**
     * @brief The new physics contribution to the ratio of the @f$\Gamma(H\to 2ud)@f$
     * in the current model and in the Standard Model. (Only terms that are
     * linear in the effective Lagrangian coefficients.)
     * @return @f$\delta \Gamma(H\to 2ud)@f$/@f$\Gamma(H\to 2ud)_{\mathrm{SM}}@f$
     */
    double deltaGammaH2udRatio1() const;
    /**
     * @brief The new physics contribution to the ratio of the @f$\Gamma(H\to 2ud)@f$
     * in the current model and in the Standard Model. (Only terms that are
     * quadratic in the effective Lagrangian coefficients.)
     * @return @f$\delta \Gamma(H\to 2ud)@f$/@f$\Gamma(H\to 2ud)_{\mathrm{SM}}@f$
     */
    double deltaGammaH2udRatio2() const;
    /**
     * @brief The ratio of the Br@f$(H\to 2ud)@f$ in the current model
     * and in the Standard Model.
     * @return Br@f$(H\to 2ud)@f$/Br@f$(H\to 2ud)_{\mathrm{SM}}@f$
     */
    virtual double BrH2udRatio() const;
    
    /**
     * @brief The ratio of the @f$\Gamma(H\to 2Lv)@f$ (@f$L=e,\mu,\tau@f$) in the current model
     * and in the Standard Model.
     * @return @f$\Gamma(H\to 2Lv)@f$/@f$\Gamma(H\to 2Lv)_{\mathrm{SM}}@f$
     */
    double GammaH2LvRatio() const;    
    /**
     * @brief The new physics contribution to the ratio of the @f$\Gamma(H\to 2Lv)@f$ (@f$L=e,\mu,\tau@f$)
     * in the current model and in the Standard Model. (Only terms that are
     * linear in the effective Lagrangian coefficients.)
     * @return @f$\delta \Gamma(H\to 2Lv)@f$/@f$\Gamma(H\to 2Lv)_{\mathrm{SM}}@f$
     */
    double deltaGammaH2LvRatio1() const;
    /**
     * @brief The new physics contribution to the ratio of the @f$\Gamma(H\to 2Lv)@f$ (@f$L=e,\mu,\tau@f$)
     * in the current model and in the Standard Model. (Only terms that are
     * quadratic in the effective Lagrangian coefficients.)
     * @return @f$\delta \Gamma(H\to 2Lv)@f$/@f$\Gamma(H\to 2Lv)_{\mathrm{SM}}@f$
     */
    double deltaGammaH2LvRatio2() const;
    /**
     * @brief The ratio of the Br@f$(H\to 2Lv)@f$ (@f$L=e,\mu,\tau@f$) in the current model
     * and in the Standard Model.
     * @return Br@f$(H\to 2Lv)@f$/Br@f$(H\to 2Lv)_{\mathrm{SM}}@f$
     */
    virtual double BrH2LvRatio() const;
    
    /**
     * @brief The ratio of the @f$\Gamma(H\to 2Lv)@f$ (@f$L=e,\mu@f$) in the current model
     * and in the Standard Model.
     * @return @f$\Gamma(H\to 2Lv)@f$/@f$\Gamma(H\to 2Lv)_{\mathrm{SM}}@f$
     */
    double GammaH2Lv2Ratio() const;    
    /**
     * @brief The new physics contribution to the ratio of the @f$\Gamma(H\to 2Lv)@f$ (@f$L=e,\mu@f$)
     * in the current model and in the Standard Model. (Only terms that are
     * linear in the effective Lagrangian coefficients.)
     * @return @f$\delta \Gamma(H\to 2Lv)@f$/@f$\Gamma(H\to 2Lv)_{\mathrm{SM}}@f$
     */
    double deltaGammaH2Lv2Ratio1() const;
    /**
     * @brief The new physics contribution to the ratio of the @f$\Gamma(H\to 2Lv)@f$ (@f$L=e,\mu@f$)
     * in the current model and in the Standard Model. (Only terms that are
     * quadratic in the effective Lagrangian coefficients.)
     * @return @f$\delta \Gamma(H\to 2Lv)@f$/@f$\Gamma(H\to 2Lv)_{\mathrm{SM}}@f$
     */
    double deltaGammaH2Lv2Ratio2() const;
    /**
     * @brief The ratio of the Br@f$(H\to 2Lv)@f$ (@f$L=e,\mu@f$) in the current model
     * and in the Standard Model.
     * @return Br@f$(H\to 2Lv)@f$/Br@f$(H\to 2Lv)_{\mathrm{SM}}@f$
     */
    virtual double BrH2Lv2Ratio() const;
    
    /**
     * @brief The ratio of the @f$\Gamma(H\to 4f)@f$ in the current model
     * and in the Standard Model.
     * @return @f$\Gamma(H\to 4f)@f$/@f$\Gamma(H\to 4f)_{\mathrm{SM}}@f$
     */
    double GammaH4fRatio() const;    
    /**
     * @brief The new physics contribution to the ratio of the @f$\Gamma(H\to 4f)@f$
     * in the current model and in the Standard Model. (Only terms that are
     * linear in the effective Lagrangian coefficients.)
     * @return @f$\delta \Gamma(H\to 4f)@f$/@f$\Gamma(H\to 4f)_{\mathrm{SM}}@f$
     */
    double deltaGammaH4fRatio1() const;
    /**
     * @brief The new physics contribution to the ratio of the @f$\Gamma(H\to 4f)@f$
     * in the current model and in the Standard Model. (Only terms that are
     * quadratic in the effective Lagrangian coefficients.)
     * @return @f$\delta \Gamma(H\to 4f)@f$/@f$\Gamma(H\to 4f)_{\mathrm{SM}}@f$
     */
    double deltaGammaH4fRatio2() const;
    /**
     * @brief The ratio of the Br@f$(H\to 4f)@f$ in the current model
     * and in the Standard Model.
     * @return Br@f$(H\to 4f)@f$/Br@f$(H\to 4f)_{\mathrm{SM}}@f$
     */
    virtual double BrH4fRatio() const;
    
    // DECAYS INVOLVING ONLY ELECTRONS, MUONS OR NEUTRINOS IN THE FINAL STATES 
    
    /**
     * @brief The ratio of the @f$\Gamma(H\to 4l)@f$ (@f$l=e,\mu@f$) in the current model
     * and in the Standard Model.
     * @return @f$\Gamma(H\to 4l)@f$/@f$\Gamma(H\to 4l)_{\mathrm{SM}}@f$
     */
    double GammaH4lRatio() const;    
    /**
     * @brief The new physics contribution to the ratio of the @f$\Gamma(H\to 4l)@f$ (@f$l=e,\mu@f$)
     * in the current model and in the Standard Model. (Only terms that are
     * linear in the effective Lagrangian coefficients.)
     * @return @f$\delta \Gamma(H\to 4l)@f$/@f$\Gamma(H\to 4l)_{\mathrm{SM}}@f$
     */
    double deltaGammaH4lRatio1() const;
    /**
     * @brief The new physics contribution to the ratio of the @f$\Gamma(H\to 4l)@f$ (@f$l=e,\mu@f$)
     * in the current model and in the Standard Model. (Only terms that are
     * quadratic in the effective Lagrangian coefficients.)
     * @return @f$\delta \Gamma(H\to 4l)@f$/@f$\Gamma(H\to 4l)_{\mathrm{SM}}@f$
     */
    double deltaGammaH4lRatio2() const;
    /**
     * @brief The ratio of the Br@f$(H\to 4l)@f$ (@f$l=e,\mu@f$) in the current model
     * and in the Standard Model.
     * @return Br@f$(H\to 4l)@f$/Br@f$(H\to 4l)_{\mathrm{SM}}@f$
     */
    virtual double BrH4lRatio() const;
    
    /**
     * @brief The ratio of the @f$\Gamma(H\to 2l2v)@f$ (@f$l=e,\mu@f$) in the current model
     * and in the Standard Model.
     * @return @f$\Gamma(H\to 2l2v)@f$/@f$\Gamma(H\to 2l2v)_{\mathrm{SM}}@f$
     */
    double GammaH2l2vRatio() const;    
    /**
     * @brief The new physics contribution to the ratio of the @f$\Gamma(H\to 2l2v)@f$ (@f$l=e,\mu@f$)
     * in the current model and in the Standard Model. (Only terms that are
     * linear in the effective Lagrangian coefficients.)
     * @return @f$\delta \Gamma(H\to 2l2v)@f$/@f$\Gamma(H\to 2l2v)_{\mathrm{SM}}@f$
     */
    double deltaGammaH2l2vRatio1() const;
    /**
     * @brief The new physics contribution to the ratio of the @f$\Gamma(H\to 2l2v)@f$ (@f$l=e,\mu@f$)
     * in the current model and in the Standard Model. (Only terms that are
     * quadratic in the effective Lagrangian coefficients.)
     * @return @f$\delta \Gamma(H\to 2l2v)@f$/@f$\Gamma(H\to 2l2v)_{\mathrm{SM}}@f$
     */
    double deltaGammaH2l2vRatio2() const;
    /**
     * @brief The ratio of the Br@f$(H\to 2l2v)@f$ (@f$l=e,\mu@f$) in the current model
     * and in the Standard Model.
     * @return Br@f$(H\to 2l2v)@f$/Br@f$(H\to 2l2v)_{\mathrm{SM}}@f$
     */
    virtual double BrH2l2vRatio() const;
    
    /**
     * @brief The ratio of the @f$\Gamma(H\to lvvl)@f$ (@f$l=e,\mu@f$) in the current model
     * and in the Standard Model.
     * @return @f$\Gamma(H\to lvvl)@f$/@f$\Gamma(H\to lvvl)_{\mathrm{SM}}@f$
     */
    double GammaHlvvlRatio() const;    
    /**
     * @brief The new physics contribution to the ratio of the @f$\Gamma(H\to lvvl)@f$ (@f$l=e,\mu@f$)
     * in the current model and in the Standard Model. (Only terms that are
     * linear in the effective Lagrangian coefficients.)
     * @return @f$\delta \Gamma(H\to lvvl)@f$/@f$\Gamma(H\to lvvl)_{\mathrm{SM}}@f$
     */
    double deltaGammaHlvvlRatio1() const;
    /**
     * @brief The new physics contribution to the ratio of the @f$\Gamma(H\to lvvl)@f$ (@f$l=e,\mu@f$)
     * in the current model and in the Standard Model. (Only terms that are
     * quadratic in the effective Lagrangian coefficients.)
     * @return @f$\delta \Gamma(H\to lvvl)@f$/@f$\Gamma(H\to lvvl)_{\mathrm{SM}}@f$
     */
    double deltaGammaHlvvlRatio2() const;
    /**
     * @brief The ratio of the Br@f$(H\to lvvl)@f$ (@f$l=e,\mu@f$) in the current model
     * and in the Standard Model.
     * @return Br@f$(H\to lvvl)@f$/Br@f$(H\to lvvl)_{\mathrm{SM}}@f$
     */
    virtual double BrHlvvlRatio() const;
    
    ///////////////////////OTHER HIGGS BRANCHING RATIOS/////////////////////////    

    /**
     * @brief The branching ratio of the of the Higgs into exotic particles.
     * @return Br@f$(H\to exotic)@f$
     */
    virtual double Br_H_exo() const;    
    /**
     * @brief The branching ratio of the of the Higgs into invisible particles.
     * @return Br@f$(H\to invisible)@f$
     */
    virtual double Br_H_inv() const;    
    /**
     * @brief The branching ratio of the of the Higgs into invisible particles 
     * (only invisible new particles).
     * @return Br@f$(H\to invisible,NP)@f$
     */
    virtual double Br_H_inv_NP() const;   
    /**
     * @brief The ratio of the Br@f$(H\to visible)@f$ in the current model
     * and in the Standard Model.
     * @return Br@f$(H\to visible)@f$/Br@f$(H\to visible)_{\mathrm{SM}}@f$
     */
    virtual double BrHvisRatio() const;
    /**
     * @brief The ratio of the Br@f$(H\to invisible)@f$ in the current model
     * and in the Standard Model.
     * @return Br@f$(H\to invisible)@f$/Br@f$(H\to ZZ \to invisible)_{\mathrm{SM}}@f$
     */
    virtual double BrHtoinvRatio() const;

    ///////////////////////SPECIAL OBSERVABLES/////////////////////////
    
    /**
     * @brief The ratio @f$\sigma(ttH)/\sigma(ttZ)@f$ 
     * in the @f$H,Z\to b\bar{b}@f$ channel in the current model and in the Standard Model.
     * @param[in] sqrt_s the center-of-mass energy in TeV
     * @return @f$\sigma(ttH)/\sigma(ttZ)@f$ normalized to the SM
     */
    virtual double muttHZbbboost(const double sqrt_s) const;

    /**
     * @brief The ratio @f$\mu_{ggH,\gamma\gamma}@f$ between the gluon-gluon fusion Higgs
     * production cross-section with subsequent decay into 2 photons in the
     * current model and in the Standard Model.
     * @param[in] sqrt_s the center-of-mass energy in TeV
     * @return @f$\mu_{ggH,\gamma\gamma}@f$
     */    
    virtual double muggHgaga(const double sqrt_s) const;    
    /**
     * @brief The ratio @f$\mu_{VBF,\gamma\gamma}@f$ between the VBF Higgs
     * production cross-section with subsequent decay into 2 photons in the
     * current model and in the Standard Model.
     * @param[in] sqrt_s the center-of-mass energy in TeV
     * @return @f$\mu_{VBF,\gamma\gamma}@f$
     */    
    virtual double muVBFHgaga(const double sqrt_s) const;
    /**
     * @brief The ratio @f$\mu_{ZH,\gamma\gamma}@f$ between the ZH
     * production cross-section with subsequent decay into 2 photons in the
     * current model and in the Standard Model.
     * @param[in] sqrt_s the center-of-mass energy in TeV
     * @return @f$\mu_{ZH,\gamma\gamma}@f$
     */
    virtual double muZHgaga(const double sqrt_s) const;
    /**
     * @brief The ratio @f$\mu_{WH,\gamma\gamma}@f$ between the WH
     * production cross-section with subsequent decay into 2 photons in the
     * current model and in the Standard Model.
     * @param[in] sqrt_s the center-of-mass energy in TeV
     * @return @f$\mu_{WH,\gamma\gamma}@f$
     */
    virtual double muWHgaga(const double sqrt_s) const;
    /**
     * @brief The ratio @f$\mu_{VH,\gamma\gamma}@f$ between the VH
     * production cross-section with subsequent decay into 2 photons in the
     * current model and in the Standard Model.
     * @param[in] sqrt_s the center-of-mass energy in TeV
     * @return @f$\mu_{VH,\gamma\gamma}@f$
     */
    virtual double muVHgaga(const double sqrt_s) const;
    /**
     * @brief The ratio @f$\mu_{ttH,\gamma\gamma}@f$ between the ttH
     * production cross-section with subsequent decay into 2 photons in the
     * current model and in the Standard Model.
     * @param[in] sqrt_s the center-of-mass energy in TeV
     * @return @f$\mu_{ttH,\gamma\gamma}@f$
     */
    virtual double muttHgaga(const double sqrt_s) const;  
    /**
     * @brief The ratio @f$\mu_{ggH,Z\gamma}@f$ between the gluon-gluon fusion Higgs
     * production cross-section with subsequent decay into @f$Z \gamma@f$ in the
     * current model and in the Standard Model.
     * @param[in] sqrt_s the center-of-mass energy in TeV
     * @return @f$\mu_{ggH,Z\gamma}@f$
     */
    virtual double muggHZga(const double sqrt_s) const;
    /**
     * @brief The ratio @f$\mu_{VBF,Z\gamma}@f$ between the VBF Higgs
     * production cross-section with subsequent decay into @f$Z \gamma@f$ in the
     * current model and in the Standard Model.
     * @param[in] sqrt_s the center-of-mass energy in TeV
     * @return @f$\mu_{VBF,Z\gamma}@f$
     */
    virtual double muVBFHZga(const double sqrt_s) const;    
    /**
     * @brief The ratio @f$\mu_{ZH,Z\gamma}@f$ between the ZH
     * production cross-section with subsequent decay into @f$Z \gamma@f$ in the
     * current model and in the Standard Model.
     * @param[in] sqrt_s the center-of-mass energy in TeV
     * @return @f$\mu_{ZH,Z\gamma}@f$
     */
    virtual double muZHZga(const double sqrt_s) const;
    /**
     * @brief The ratio @f$\mu_{WH,Z\gamma}@f$ between the WH
     * production cross-section with subsequent decay into @f$Z \gamma@f$ in the
     * current model and in the Standard Model.
     * @param[in] sqrt_s the center-of-mass energy in TeV
     * @return @f$\mu_{WH,Z\gamma}@f$
     */
    virtual double muWHZga(const double sqrt_s) const;    
    /**
     * @brief The ratio @f$\mu_{VH,Z\gamma}@f$ between the VH
     * production cross-section with subsequent decay into @f$Z \gamma@f$ in the
     * current model and in the Standard Model.
     * @param[in] sqrt_s the center-of-mass energy in TeV
     * @return @f$\mu_{VH,Z\gamma}@f$
     */
    virtual double muVHZga(const double sqrt_s) const;
    /**
     * @brief The ratio @f$\mu_{ttH,Z\gamma}@f$ between the ttH
     * production cross-section with subsequent decay into @f$Z \gamma@f$ in the
     * current model and in the Standard Model.
     * @param[in] sqrt_s the center-of-mass energy in TeV
     * @return @f$\mu_{ttH,Z\gamma}@f$
     */
    virtual double muttHZga(const double sqrt_s) const;    
    /**
     * @brief The ratio @f$\mu_{ggH,ZZ}@f$ between the gluon-gluon fusion Higgs
     * production cross-section with subsequent decay into @f$Z Z^*@f$ in the
     * current model and in the Standard Model.
     * @param[in] sqrt_s the center-of-mass energy in TeV
     * @return @f$\mu_{ggH,ZZ}@f$
     */
    virtual double muggHZZ(const double sqrt_s) const;
    /**
     * @brief The ratio @f$\mu_{VBF,ZZ}@f$ between the VBF Higgs
     * production cross-section with subsequent decay into @f$Z Z^*@f$ in the
     * current model and in the Standard Model.
     * @param[in] sqrt_s the center-of-mass energy in TeV
     * @return @f$\mu_{VBF,ZZ}@f$
     */
    virtual double muVBFHZZ(const double sqrt_s) const;    
    /**
     * @brief The ratio @f$\mu_{ZH,ZZ}@f$ between the ZH
     * production cross-section with subsequent decay into @f$Z Z^*@f$ in the
     * current model and in the Standard Model.
     * @param[in] sqrt_s the center-of-mass energy in TeV
     * @return @f$\mu_{ZH,ZZ}@f$
     */
    virtual double muZHZZ(const double sqrt_s) const;
    /**
     * @brief The ratio @f$\mu_{WH,ZZ}@f$ between the WH
     * production cross-section with subsequent decay into @f$Z Z^*@f$ in the
     * current model and in the Standard Model.
     * @param[in] sqrt_s the center-of-mass energy in TeV
     * @return @f$\mu_{WH,ZZ}@f$
     */
    virtual double muWHZZ(const double sqrt_s) const;    
    /**
     * @brief The ratio @f$\mu_{VH,ZZ}@f$ between the VH
     * production cross-section with subsequent decay into @f$Z Z^*@f$ in the
     * current model and in the Standard Model.
     * @param[in] sqrt_s the center-of-mass energy in TeV
     * @return @f$\mu_{VH,ZZ}@f$
     */
    virtual double muVHZZ(const double sqrt_s) const;
    /**
     * @brief The ratio @f$\mu_{ttH,ZZ}@f$ between the ttH
     * production cross-section with subsequent decay into @f$Z Z^*@f$ in the
     * current model and in the Standard Model.
     * @param[in] sqrt_s the center-of-mass energy in TeV
     * @return @f$\mu_{ttH,ZZ}@f$
     */
    virtual double muttHZZ(const double sqrt_s) const;
        
    /**
     * @brief The ratio @f$\mu_{ggH,ZZ\to 4l}@f$ between the gluon-gluon fusion Higgs
     * production cross-section with subsequent decay into @f$Z Z^*\to 4l@f$ in the
     * current model and in the Standard Model.
     * @param[in] sqrt_s the center-of-mass energy in TeV
     * @return @f$\mu_{ggH,ZZ\to 4l}@f$
     */
    virtual double muggHZZ4l(const double sqrt_s) const;
    /**
     * @brief The ratio @f$\mu_{VBF,ZZ\to 4l}@f$ between the VBF Higgs
     * production cross-section with subsequent decay into @f$Z Z^*\to 4l@f$ in the
     * current model and in the Standard Model.
     * @param[in] sqrt_s the center-of-mass energy in TeV
     * @return @f$\mu_{VBF,ZZ\to 4l}@f$
     */
    virtual double muVBFHZZ4l(const double sqrt_s) const;    
    /**
     * @brief The ratio @f$\mu_{ZH,ZZ\to 4l}@f$ between the ZH
     * production cross-section with subsequent decay into @f$Z Z^*\to 4l@f$ in the
     * current model and in the Standard Model.
     * @param[in] sqrt_s the center-of-mass energy in TeV
     * @return @f$\mu_{ZH,ZZ\to 4l}@f$
     */
    virtual double muZHZZ4l(const double sqrt_s) const;
    /**
     * @brief The ratio @f$\mu_{WH,ZZ\to 4l}@f$ between the WH
     * production cross-section with subsequent decay into @f$Z Z^*\to 4l@f$ in the
     * current model and in the Standard Model.
     * @param[in] sqrt_s the center-of-mass energy in TeV
     * @return @f$\mu_{WH,ZZ\to 4l}@f$
     */
    virtual double muWHZZ4l(const double sqrt_s) const;    
    /**
     * @brief The ratio @f$\mu_{VH,ZZ\to 4l}@f$ between the VH
     * production cross-section with subsequent decay into @f$Z Z^*\to 4l@f$ in the
     * current model and in the Standard Model.
     * @param[in] sqrt_s the center-of-mass energy in TeV
     * @return @f$\mu_{VH,ZZ\to 4l}@f$
     */
    virtual double muVHZZ4l(const double sqrt_s) const;
    /**
     * @brief The ratio @f$\mu_{ttH,ZZ\to 4l}@f$ between the ttH
     * production cross-section with subsequent decay into @f$Z Z^*\to 4l@f$ in the
     * current model and in the Standard Model.
     * @param[in] sqrt_s the center-of-mass energy in TeV
     * @return @f$\mu_{ttH,ZZ\to 4l}@f$
     */
    virtual double muttHZZ4l(const double sqrt_s) const;
       
    /**
     * @brief The ratio @f$\mu_{ggH,WW}@f$ between the gluon-gluon fusion Higgs
     * production cross-section with subsequent decay into @f$W W^*@f$ in the
     * current model and in the Standard Model.
     * @param[in] sqrt_s the center-of-mass energy in TeV
     * @return @f$\mu_{ggH,WW}@f$
     */
    virtual double muggHWW(const double sqrt_s) const;
    /**
     * @brief The ratio @f$\mu_{VBF,WW}@f$ between the VBF Higgs
     * production cross-section with subsequent decay into @f$W W^*@f$ in the
     * current model and in the Standard Model.
     * @param[in] sqrt_s the center-of-mass energy in TeV
     * @return @f$\mu_{VBF,WW}@f$
     */
    virtual double muVBFHWW(const double sqrt_s) const;   
    /**
     * @brief The ratio @f$\mu_{ZH,WW}@f$ between the ZH
     * production cross-section with subsequent decay into @f$W W^*@f$ in the
     * current model and in the Standard Model.
     * @param[in] sqrt_s the center-of-mass energy in TeV
     * @return @f$\mu_{ZH,WW}@f$
     */
    virtual double muZHWW(const double sqrt_s) const;
    /**
     * @brief The ratio @f$\mu_{WH,WW}@f$ between the WH
     * production cross-section with subsequent decay into @f$W W^*@f$ in the
     * current model and in the Standard Model.
     * @param[in] sqrt_s the center-of-mass energy in TeV
     * @return @f$\mu_{WH,WW}@f$
     */
    virtual double muWHWW(const double sqrt_s) const;    
    /**
     * @brief The ratio @f$\mu_{VH,WW}@f$ between the VH
     * production cross-section with subsequent decay into @f$W W^*@f$ in the
     * current model and in the Standard Model.
     * @param[in] sqrt_s the center-of-mass energy in TeV
     * @return @f$\mu_{VH,WW}@f$
     */
    virtual double muVHWW(const double sqrt_s) const;
    /**
     * @brief The ratio @f$\mu_{ttH,WW}@f$ between the ttH
     * production cross-section with subsequent decay into @f$W W^*@f$ in the
     * current model and in the Standard Model.
     * @param[in] sqrt_s the center-of-mass energy in TeV
     * @return @f$\mu_{ttH,WW}@f$
     */
    virtual double muttHWW(const double sqrt_s) const;   
    /**
     * @brief The ratio @f$\mu_{ggH,WW\to 2l2\nu}@f$ between the gluon-gluon fusion Higgs
     * production cross-section with subsequent decay into @f$W W^*\to 2l2\nu@f$ in the
     * current model and in the Standard Model.
     * @param[in] sqrt_s the center-of-mass energy in TeV
     * @return @f$\mu_{ggH,WW\to 2l2\nu}@f$
     */
    virtual double muggHWW2l2v(const double sqrt_s) const;
    /**
     * @brief The ratio @f$\mu_{VBF,WW\to 2l2\nu}@f$ between the VBF Higgs
     * production cross-section with subsequent decay into @f$W W^*\to 2l2\nu@f$ in the
     * current model and in the Standard Model.
     * @param[in] sqrt_s the center-of-mass energy in TeV
     * @return @f$\mu_{VBF,WW\to 2l2\nu}@f$
     */
    virtual double muVBFHWW2l2v(const double sqrt_s) const;   
    /**
     * @brief The ratio @f$\mu_{ZH,WW\to 2l2\nu}@f$ between the ZH
     * production cross-section with subsequent decay into @f$W W^*\to 2l2\nu@f$ in the
     * current model and in the Standard Model.
     * @param[in] sqrt_s the center-of-mass energy in TeV
     * @return @f$\mu_{ZH,WW\to 2l2\nu}@f$
     */
    virtual double muZHWW2l2v(const double sqrt_s) const;
    /**
     * @brief The ratio @f$\mu_{WH,WW\to 2l2\nu}@f$ between the WH
     * production cross-section with subsequent decay into @f$W W^*\to 2l2\nu@f$ in the
     * current model and in the Standard Model.
     * @param[in] sqrt_s the center-of-mass energy in TeV
     * @return @f$\mu_{WH,WW\to 2l2\nu}@f$
     */
    virtual double muWHWW2l2v(const double sqrt_s) const;    
    /**
     * @brief The ratio @f$\mu_{VH,WW\to 2l2\nu}@f$ between the VH
     * production cross-section with subsequent decay into @f$W W^*\to 2l2\nu@f$ in the
     * current model and in the Standard Model.
     * @param[in] sqrt_s the center-of-mass energy in TeV
     * @return @f$\mu_{VH,WW\to 2l2\nu}@f$
     */
    virtual double muVHWW2l2v(const double sqrt_s) const;
    /**
     * @brief The ratio @f$\mu_{ttH,WW\to 2l2\nu}@f$ between the ttH
     * production cross-section with subsequent decay into @f$W W^*\to 2l2\nu@f$ in the
     * current model and in the Standard Model.
     * @param[in] sqrt_s the center-of-mass energy in TeV
     * @return @f$\mu_{ttH,WW\to 2l2\nu}@f$
     */
    virtual double muttHWW2l2v(const double sqrt_s) const;  
    /**
     * @brief The ratio @f$\mu_{ggH,\mu\mu}@f$ between the gluon-gluon fusion Higgs
     * production cross-section with subsequent decay into @f$\mu\mu@f$ in the
     * current model and in the Standard Model.
     * @param[in] sqrt_s the center-of-mass energy in TeV
     * @return @f$\mu_{ggH,\mu\mu}@f$
     */
    virtual double muggHmumu(const double sqrt_s) const;
    /**
     * @brief The ratio @f$\mu_{VBF,\mu\mu}@f$ between the VBF Higgs
     * production cross-section with subsequent decay into @f$\mu\mu@f$ in the
     * current model and in the Standard Model.
     * @param[in] sqrt_s the center-of-mass energy in TeV
     * @return @f$\mu_{VBF,\mu\mu}@f$
     */
    virtual double muVBFHmumu(const double sqrt_s) const;    
    /**
     * @brief The ratio @f$\mu_{ZH,\mu\mu}@f$ between the ZH
     * production cross-section with subsequent decay into @f$\mu\mu@f$ in the
     * current model and in the Standard Model.
     * @param[in] sqrt_s the center-of-mass energy in TeV
     * @return @f$\mu_{ZH,\mu\mu}@f$
     */
    virtual double muZHmumu(const double sqrt_s) const;
    /**
     * @brief The ratio @f$\mu_{WH,\mu\mu}@f$ between the WH
     * production cross-section with subsequent decay into @f$\mu\mu@f$ in the
     * current model and in the Standard Model.
     * @param[in] sqrt_s the center-of-mass energy in TeV
     * @return @f$\mu_{WH,\mu\mu}@f$
     */
    virtual double muWHmumu(const double sqrt_s) const;    
    /**
     * @brief The ratio @f$\mu_{VH,\mu\mu}@f$ between the VH
     * production cross-section with subsequent decay into @f$\mu\mu@f$ in the
     * current model and in the Standard Model.
     * @param[in] sqrt_s the center-of-mass energy in TeV
     * @return @f$\mu_{VH,\mu\mu}@f$
     */
    virtual double muVHmumu(const double sqrt_s) const;
    /**
     * @brief The ratio @f$\mu_{ttH,\mu\mu}@f$ between the ttH
     * production cross-section with subsequent decay into @f$\mu\mu@f$ in the
     * current model and in the Standard Model.
     * @param[in] sqrt_s the center-of-mass energy in TeV
     * @return @f$\mu_{ttH,\mu\mu}@f$
     */
    virtual double muttHmumu(const double sqrt_s) const;
    /**
     * @brief The ratio @f$\mu_{ggH,\tau\tau}@f$ between the gluon-gluon fusion Higgs
     * production cross-section with subsequent decay into @f$\tau\tau@f$ in the
     * current model and in the Standard Model.
     * @param[in] sqrt_s the center-of-mass energy in TeV
     * @return @f$\mu_{ggH,\tau\tau}@f$
     */
    virtual double muggHtautau(const double sqrt_s) const;
    /**
     * @brief The ratio @f$\mu_{VBF,\tau\tau}@f$ between the VBF Higgs
     * production cross-section with subsequent decay into @f$\tau\tau@f$ in the
     * current model and in the Standard Model.
     * @param[in] sqrt_s the center-of-mass energy in TeV
     * @return @f$\mu_{VBF,\tau\tau}@f$
     */
    virtual double muVBFHtautau(const double sqrt_s) const;    
    /**
     * @brief The ratio @f$\mu_{ZH,\tau\tau}@f$ between the ZH
     * production cross-section with subsequent decay into @f$\tau\tau@f$ in the
     * current model and in the Standard Model.
     * @param[in] sqrt_s the center-of-mass energy in TeV
     * @return @f$\mu_{ZH,\tau\tau}@f$
     */
    virtual double muZHtautau(const double sqrt_s) const;
    /**
     * @brief The ratio @f$\mu_{WH,\tau\tau}@f$ between the WH
     * production cross-section with subsequent decay into @f$\tau\tau@f$ in the
     * current model and in the Standard Model.
     * @param[in] sqrt_s the center-of-mass energy in TeV
     * @return @f$\mu_{WH,\tau\tau}@f$
     */
    virtual double muWHtautau(const double sqrt_s) const;    
    /**
     * @brief The ratio @f$\mu_{VH,\tau\tau}@f$ between the VH
     * production cross-section with subsequent decay into @f$\tau\tau@f$ in the
     * current model and in the Standard Model.
     * @param[in] sqrt_s the center-of-mass energy in TeV
     * @return @f$\mu_{VH,\tau\tau}@f$
     */
    virtual double muVHtautau(const double sqrt_s) const;
    /**
     * @brief The ratio @f$\mu_{ttH,\tau\tau}@f$ between the ttH
     * production cross-section with subsequent decay into @f$\tau\tau@f$ in the
     * current model and in the Standard Model.
     * @param[in] sqrt_s the center-of-mass energy in TeV
     * @return @f$\mu_{ttH,\tau\tau}@f$
     */
    virtual double muttHtautau(const double sqrt_s) const;
    /**
     * @brief The ratio @f$\mu_{ggH,bb}@f$ between the gluon-gluon fusion Higgs
     * production cross-section with subsequent decay into @f$bb@f$ in the
     * current model and in the Standard Model.
     * @param[in] sqrt_s the center-of-mass energy in TeV
     * @return @f$\mu_{ggH,bb}@f$
     */
    virtual double muggHbb(const double sqrt_s) const;
    /**
     * @brief The ratio @f$\mu_{VBF,bb}@f$ between the VBF Higgs
     * production cross-section with subsequent decay into @f$bb@f$ in the
     * current model and in the Standard Model.
     * @param[in] sqrt_s the center-of-mass energy in TeV
     * @return @f$\mu_{VBF,bb}@f$
     */
    virtual double muVBFHbb(const double sqrt_s) const;    
    /**
     * @brief The ratio @f$\mu_{ZH,bb}@f$ between the ZH
     * production cross-section with subsequent decay into @f$bb@f$ in the
     * current model and in the Standard Model.
     * @param[in] sqrt_s the center-of-mass energy in TeV
     * @return @f$\mu_{ZH,bb}@f$
     */
    virtual double muZHbb(const double sqrt_s) const;
    /**
     * @brief The ratio @f$\mu_{WH,bb}@f$ between the WH
     * production cross-section with subsequent decay into @f$bb@f$ in the
     * current model and in the Standard Model.
     * @param[in] sqrt_s the center-of-mass energy in TeV
     * @return @f$\mu_{WH,bb}@f$
     */
    virtual double muWHbb(const double sqrt_s) const;
    /**
     * @brief The ratio @f$\mu_{VH,bb}@f$ between the VH
     * production cross-section with subsequent decay into @f$bb@f$ in the
     * current model and in the Standard Model.
     * @param[in] sqrt_s the center-of-mass energy in TeV
     * @return @f$\mu_{VH,bb}@f$
     */
    virtual double muVHbb(const double sqrt_s) const;
    /**
     * @brief The ratio @f$\mu_{ttH,bb}@f$ between the ttH
     * production cross-section with subsequent decay into @f$bb@f$ in the
     * current model and in the Standard Model.
     * @param[in] sqrt_s the center-of-mass energy in TeV
     * @return @f$\mu_{ttH,bb}@f$
     */
    virtual double muttHbb(const double sqrt_s) const;
    
////////////////////////////////////////////////////////////////////////////////////////////
//-----------------------------------------------------------------------------------------
//-- Special Hadron collider signal strengths with separate full TH unc U(prod x decay) ---
//-----------------------------------------------------------------------------------------
//////////////////////////////////////////////////////////////////////////////////////////// 
    
    /**
     * @brief The ratio @f$\mu_{ggH,\gamma\gamma}@f$ between the gluon-gluon fusion Higgs
     * production cross-section with subsequent decay into 2 photons in the
     * current model and in the Standard Model.
     * @param[in] sqrt_s the center-of-mass energy in TeV
     * @return @f$\mu_{ggH,\gamma\gamma}@f$
     */    
    virtual double muTHUggHgaga(const double sqrt_s) const;    
    /**
     * @brief The ratio @f$\mu_{VBF,\gamma\gamma}@f$ between the VBF Higgs
     * production cross-section with subsequent decay into 2 photons in the
     * current model and in the Standard Model.
     * @param[in] sqrt_s the center-of-mass energy in TeV
     * @return @f$\mu_{VBF,\gamma\gamma}@f$
     */    
    virtual double muTHUVBFHgaga(const double sqrt_s) const;
    /**
     * @brief The ratio @f$\mu_{ZH,\gamma\gamma}@f$ between the ZH
     * production cross-section with subsequent decay into 2 photons in the
     * current model and in the Standard Model.
     * @param[in] sqrt_s the center-of-mass energy in TeV
     * @return @f$\mu_{ZH,\gamma\gamma}@f$
     */
    virtual double muTHUZHgaga(const double sqrt_s) const;
    /**
     * @brief The ratio @f$\mu_{WH,\gamma\gamma}@f$ between the WH
     * production cross-section with subsequent decay into 2 photons in the
     * current model and in the Standard Model.
     * @param[in] sqrt_s the center-of-mass energy in TeV
     * @return @f$\mu_{WH,\gamma\gamma}@f$
     */
    virtual double muTHUWHgaga(const double sqrt_s) const;
    /**
     * @brief The ratio @f$\mu_{VH,\gamma\gamma}@f$ between the VH
     * production cross-section with subsequent decay into 2 photons in the
     * current model and in the Standard Model.
     * @param[in] sqrt_s the center-of-mass energy in TeV
     * @return @f$\mu_{VH,\gamma\gamma}@f$
     */
    virtual double muTHUVHgaga(const double sqrt_s) const;
    /**
     * @brief The ratio @f$\mu_{ttH,\gamma\gamma}@f$ between the ttH
     * production cross-section with subsequent decay into 2 photons in the
     * current model and in the Standard Model.
     * @param[in] sqrt_s the center-of-mass energy in TeV
     * @return @f$\mu_{ttH,\gamma\gamma}@f$
     */
    virtual double muTHUttHgaga(const double sqrt_s) const;
    /**
     * @brief The ratio @f$\mu_{ggH,Z\gamma}@f$ between the gluon-gluon fusion Higgs
     * production cross-section with subsequent decay into @f$Z \gamma@f$ in the
     * current model and in the Standard Model.
     * @param[in] sqrt_s the center-of-mass energy in TeV
     * @return @f$\mu_{ggH,Z\gamma}@f$
     */
    virtual double muTHUggHZga(const double sqrt_s) const;
    /**
     * @brief The ratio @f$\mu_{VBF,Z\gamma}@f$ between the VBF Higgs
     * production cross-section with subsequent decay into @f$Z \gamma@f$ in the
     * current model and in the Standard Model.
     * @param[in] sqrt_s the center-of-mass energy in TeV
     * @return @f$\mu_{VBF,Z\gamma}@f$
     */
    virtual double muTHUVBFHZga(const double sqrt_s) const;    
    /**
     * @brief The ratio @f$\mu_{ZH,Z\gamma}@f$ between the ZH
     * production cross-section with subsequent decay into @f$Z \gamma@f$ in the
     * current model and in the Standard Model.
     * @param[in] sqrt_s the center-of-mass energy in TeV
     * @return @f$\mu_{ZH,Z\gamma}@f$
     */
    virtual double muTHUZHZga(const double sqrt_s) const;
    /**
     * @brief The ratio @f$\mu_{WH,Z\gamma}@f$ between the WH
     * production cross-section with subsequent decay into @f$Z \gamma@f$ in the
     * current model and in the Standard Model.
     * @param[in] sqrt_s the center-of-mass energy in TeV
     * @return @f$\mu_{WH,Z\gamma}@f$
     */
    virtual double muTHUWHZga(const double sqrt_s) const;    
    /**
     * @brief The ratio @f$\mu_{VH,Z\gamma}@f$ between the VH
     * production cross-section with subsequent decay into @f$Z \gamma@f$ in the
     * current model and in the Standard Model.
     * @param[in] sqrt_s the center-of-mass energy in TeV
     * @return @f$\mu_{VH,Z\gamma}@f$
     */
    virtual double muTHUVHZga(const double sqrt_s) const;
    /**
     * @brief The ratio @f$\mu_{ttH,Z\gamma}@f$ between the ttH
     * production cross-section with subsequent decay into @f$Z \gamma@f$ in the
     * current model and in the Standard Model.
     * @param[in] sqrt_s the center-of-mass energy in TeV
     * @return @f$\mu_{ttH,Z\gamma}@f$
     */
    virtual double muTHUttHZga(const double sqrt_s) const; 
    /**
     * @brief The ratio @f$\mu_{ggH,ZZ}@f$ between the gluon-gluon fusion Higgs
     * production cross-section with subsequent decay into @f$Z Z^*@f$ in the
     * current model and in the Standard Model.
     * @param[in] sqrt_s the center-of-mass energy in TeV
     * @return @f$\mu_{ggH,ZZ}@f$
     */
    virtual double muTHUggHZZ(const double sqrt_s) const;
    /**
     * @brief The ratio @f$\mu_{VBF,ZZ}@f$ between the VBF Higgs
     * production cross-section with subsequent decay into @f$Z Z^*@f$ in the
     * current model and in the Standard Model.
     * @param[in] sqrt_s the center-of-mass energy in TeV
     * @return @f$\mu_{VBF,ZZ}@f$
     */
    virtual double muTHUVBFHZZ(const double sqrt_s) const;
    /**
     * @brief The ratio @f$\mu_{ZH,ZZ}@f$ between the ZH
     * production cross-section with subsequent decay into @f$Z Z^*@f$ in the
     * current model and in the Standard Model.
     * @param[in] sqrt_s the center-of-mass energy in TeV
     * @return @f$\mu_{ZH,ZZ}@f$
     */
    virtual double muTHUZHZZ(const double sqrt_s) const;
    /**
     * @brief The ratio @f$\mu_{WH,ZZ}@f$ between the WH
     * production cross-section with subsequent decay into @f$Z Z^*@f$ in the
     * current model and in the Standard Model.
     * @param[in] sqrt_s the center-of-mass energy in TeV
     * @return @f$\mu_{WH,ZZ}@f$
     */
    virtual double muTHUWHZZ(const double sqrt_s) const;
    /**
     * @brief The ratio @f$\mu_{VH,ZZ}@f$ between the VH
     * production cross-section with subsequent decay into @f$Z Z^*@f$ in the
     * current model and in the Standard Model.
     * @param[in] sqrt_s the center-of-mass energy in TeV
     * @return @f$\mu_{VH,ZZ}@f$
     */
    virtual double muTHUVHZZ(const double sqrt_s) const;
    /**
     * @brief The ratio @f$\mu_{ttH,ZZ}@f$ between the ttH
     * production cross-section with subsequent decay into @f$Z Z^*@f$ in the
     * current model and in the Standard Model.
     * @param[in] sqrt_s the center-of-mass energy in TeV
     * @return @f$\mu_{ttH,ZZ}@f$
     */
    virtual double muTHUttHZZ(const double sqrt_s) const;
    
    /**
     * @brief The ratio @f$\mu_{ggH,ZZ\to 4l}@f$ between the gluon-gluon fusion Higgs
     * production cross-section with subsequent decay into @f$Z Z^*\to 4l@f$ in the
     * current model and in the Standard Model.
     * @param[in] sqrt_s the center-of-mass energy in TeV
     * @return @f$\mu_{ggH,ZZ\to 4l}@f$
     */
    virtual double muTHUggHZZ4l(const double sqrt_s) const;
    /**
     * @brief The ratio @f$\mu_{VBF,ZZ\to 4l}@f$ between the VBF Higgs
     * production cross-section with subsequent decay into @f$Z Z^*\to 4l@f$ in the
     * current model and in the Standard Model.
     * @param[in] sqrt_s the center-of-mass energy in TeV
     * @return @f$\mu_{VBF,ZZ\to 4l}@f$
     */
    virtual double muTHUVBFHZZ4l(const double sqrt_s) const;    
    /**
     * @brief The ratio @f$\mu_{ZH,ZZ\to 4l}@f$ between the ZH
     * production cross-section with subsequent decay into @f$Z Z^*\to 4l@f$ in the
     * current model and in the Standard Model.
     * @param[in] sqrt_s the center-of-mass energy in TeV
     * @return @f$\mu_{ZH,ZZ\to 4l}@f$
     */
    virtual double muTHUZHZZ4l(const double sqrt_s) const;
    /**
     * @brief The ratio @f$\mu_{WH,ZZ\to 4l}@f$ between the WH
     * production cross-section with subsequent decay into @f$Z Z^*\to 4l@f$ in the
     * current model and in the Standard Model.
     * @param[in] sqrt_s the center-of-mass energy in TeV
     * @return @f$\mu_{WH,ZZ\to 4l}@f$
     */
    virtual double muTHUWHZZ4l(const double sqrt_s) const;    
    /**
     * @brief The ratio @f$\mu_{VH,ZZ\to 4l}@f$ between the VH
     * production cross-section with subsequent decay into @f$Z Z^*\to 4l@f$ in the
     * current model and in the Standard Model.
     * @param[in] sqrt_s the center-of-mass energy in TeV
     * @return @f$\mu_{VH,ZZ\to 4l}@f$
     */
    virtual double muTHUVHZZ4l(const double sqrt_s) const;
    /**
     * @brief The ratio @f$\mu_{ttH,ZZ\to 4l}@f$ between the ttH
     * production cross-section with subsequent decay into @f$Z Z^*\to 4l@f$ in the
     * current model and in the Standard Model.
     * @param[in] sqrt_s the center-of-mass energy in TeV
     * @return @f$\mu_{ttH,ZZ\to 4l}@f$
     */
    virtual double muTHUttHZZ4l(const double sqrt_s) const;
    
    /**
     * @brief The ratio @f$\mu_{ggH,WW}@f$ between the gluon-gluon fusion Higgs
     * production cross-section with subsequent decay into @f$W W^*@f$ in the
     * current model and in the Standard Model.
     * @param[in] sqrt_s the center-of-mass energy in TeV
     * @return @f$\mu_{ggH,WW}@f$
     */
    virtual double muTHUggHWW(const double sqrt_s) const;
    /**
     * @brief The ratio @f$\mu_{VBF,WW}@f$ between the VBF Higgs
     * production cross-section with subsequent decay into @f$W W^*@f$ in the
     * current model and in the Standard Model.
     * @param[in] sqrt_s the center-of-mass energy in TeV
     * @return @f$\mu_{VBF,WW}@f$
     */
    virtual double muTHUVBFHWW(const double sqrt_s) const;
    /**
     * @brief The ratio @f$\mu_{ZH,WW}@f$ between the ZH
     * production cross-section with subsequent decay into @f$W W^*@f$ in the
     * current model and in the Standard Model.
     * @param[in] sqrt_s the center-of-mass energy in TeV
     * @return @f$\mu_{ZH,WW}@f$
     */
    virtual double muTHUZHWW(const double sqrt_s) const;
    /**
     * @brief The ratio @f$\mu_{WH,WW}@f$ between the WH
     * production cross-section with subsequent decay into @f$W W^*@f$ in the
     * current model and in the Standard Model.
     * @param[in] sqrt_s the center-of-mass energy in TeV
     * @return @f$\mu_{WH,WW}@f$
     */
    virtual double muTHUWHWW(const double sqrt_s) const; 
    /**
     * @brief The ratio @f$\mu_{VH,WW}@f$ between the VH
     * production cross-section with subsequent decay into @f$W W^*@f$ in the
     * current model and in the Standard Model.
     * @param[in] sqrt_s the center-of-mass energy in TeV
     * @return @f$\mu_{VH,WW}@f$
     */
    virtual double muTHUVHWW(const double sqrt_s) const;
    /**
     * @brief The ratio @f$\mu_{ttH,WW}@f$ between the ttH
     * production cross-section with subsequent decay into @f$W W^*@f$ in the
     * current model and in the Standard Model.
     * @param[in] sqrt_s the center-of-mass energy in TeV
     * @return @f$\mu_{ttH,WW}@f$
     */
    virtual double muTHUttHWW(const double sqrt_s) const;
    /**
     * @brief The ratio @f$\mu_{ggH,WW\to 2l2\nu}@f$ between the gluon-gluon fusion Higgs
     * production cross-section with subsequent decay into @f$W W^*\to 2l2\nu@f$ in the
     * current model and in the Standard Model.
     * @param[in] sqrt_s the center-of-mass energy in TeV
     * @return @f$\mu_{ggH,WW\to 2l2\nu}@f$
     */
    virtual double muTHUggHWW2l2v(const double sqrt_s) const;
    /**
     * @brief The ratio @f$\mu_{VBF,WW\to 2l2\nu}@f$ between the VBF Higgs
     * production cross-section with subsequent decay into @f$W W^*\to 2l2\nu@f$ in the
     * current model and in the Standard Model.
     * @param[in] sqrt_s the center-of-mass energy in TeV
     * @return @f$\mu_{VBF,WW\to 2l2\nu}@f$
     */
    virtual double muTHUVBFHWW2l2v(const double sqrt_s) const;   
    /**
     * @brief The ratio @f$\mu_{ZH,WW\to 2l2\nu}@f$ between the ZH
     * production cross-section with subsequent decay into @f$W W^*\to 2l2\nu@f$ in the
     * current model and in the Standard Model.
     * @param[in] sqrt_s the center-of-mass energy in TeV
     * @return @f$\mu_{ZH,WW\to 2l2\nu}@f$
     */
    virtual double muTHUZHWW2l2v(const double sqrt_s) const;
    /**
     * @brief The ratio @f$\mu_{WH,WW\to 2l2\nu}@f$ between the WH
     * production cross-section with subsequent decay into @f$W W^*\to 2l2\nu@f$ in the
     * current model and in the Standard Model.
     * @param[in] sqrt_s the center-of-mass energy in TeV
     * @return @f$\mu_{WH,WW\to 2l2\nu}@f$
     */
    virtual double muTHUWHWW2l2v(const double sqrt_s) const;    
    /**
     * @brief The ratio @f$\mu_{VH,WW\to 2l2\nu}@f$ between the VH
     * production cross-section with subsequent decay into @f$W W^*\to 2l2\nu@f$ in the
     * current model and in the Standard Model.
     * @param[in] sqrt_s the center-of-mass energy in TeV
     * @return @f$\mu_{VH,WW\to 2l2\nu}@f$
     */
    virtual double muTHUVHWW2l2v(const double sqrt_s) const;
    /**
     * @brief The ratio @f$\mu_{ttH,WW\to 2l2\nu}@f$ between the ttH
     * production cross-section with subsequent decay into @f$W W^*\to 2l2\nu@f$ in the
     * current model and in the Standard Model.
     * @param[in] sqrt_s the center-of-mass energy in TeV
     * @return @f$\mu_{ttH,WW\to 2l2\nu}@f$
     */
    virtual double muTHUttHWW2l2v(const double sqrt_s) const;
    /**
     * @brief The ratio @f$\mu_{ggH,\mu\mu}@f$ between the gluon-gluon fusion Higgs
     * production cross-section with subsequent decay into @f$\mu\mu@f$ in the
     * current model and in the Standard Model.
     * @param[in] sqrt_s the center-of-mass energy in TeV
     * @return @f$\mu_{ggH,\mu\mu}@f$
     */
    virtual double muTHUggHmumu(const double sqrt_s) const;
    /**
     * @brief The ratio @f$\mu_{VBF,\mu\mu}@f$ between the VBF Higgs
     * production cross-section with subsequent decay into @f$\mu\mu@f$ in the
     * current model and in the Standard Model.
     * @param[in] sqrt_s the center-of-mass energy in TeV
     * @return @f$\mu_{VBF,\mu\mu}@f$
     */
    virtual double muTHUVBFHmumu(const double sqrt_s) const;    
    /**
     * @brief The ratio @f$\mu_{ZH,\mu\mu}@f$ between the ZH
     * production cross-section with subsequent decay into @f$\mu\mu@f$ in the
     * current model and in the Standard Model.
     * @param[in] sqrt_s the center-of-mass energy in TeV
     * @return @f$\mu_{ZH,\mu\mu}@f$
     */
    virtual double muTHUZHmumu(const double sqrt_s) const;
    /**
     * @brief The ratio @f$\mu_{WH,\mu\mu}@f$ between the WH
     * production cross-section with subsequent decay into @f$\mu\mu@f$ in the
     * current model and in the Standard Model.
     * @param[in] sqrt_s the center-of-mass energy in TeV
     * @return @f$\mu_{WH,\mu\mu}@f$
     */
    virtual double muTHUWHmumu(const double sqrt_s) const;    
    /**
     * @brief The ratio @f$\mu_{VH,\mu\mu}@f$ between the VH
     * production cross-section with subsequent decay into @f$\mu\mu@f$ in the
     * current model and in the Standard Model.
     * @param[in] sqrt_s the center-of-mass energy in TeV
     * @return @f$\mu_{VH,\mu\mu}@f$
     */
    virtual double muTHUVHmumu(const double sqrt_s) const;
    /**
     * @brief The ratio @f$\mu_{ttH,\mu\mu}@f$ between the ttH
     * production cross-section with subsequent decay into @f$\mu\mu@f$ in the
     * current model and in the Standard Model.
     * @param[in] sqrt_s the center-of-mass energy in TeV
     * @return @f$\mu_{ttH,\mu\mu}@f$
     */
    virtual double muTHUttHmumu(const double sqrt_s) const;
    /**
     * @brief The ratio @f$\mu_{ggH,\tau\tau}@f$ between the gluon-gluon fusion Higgs
     * production cross-section with subsequent decay into @f$\tau\tau@f$ in the
     * current model and in the Standard Model.
     * @param[in] sqrt_s the center-of-mass energy in TeV
     * @return @f$\mu_{ggH,\tau\tau}@f$
     */
    virtual double muTHUggHtautau(const double sqrt_s) const;
    /**
     * @brief The ratio @f$\mu_{VBF,\tau\tau}@f$ between the VBF Higgs
     * production cross-section with subsequent decay into @f$\tau\tau@f$ in the
     * current model and in the Standard Model.
     * @param[in] sqrt_s the center-of-mass energy in TeV
     * @return @f$\mu_{VBF,\tau\tau}@f$
     */
    virtual double muTHUVBFHtautau(const double sqrt_s) const;
    /**
     * @brief The ratio @f$\mu_{ZH,\tau\tau}@f$ between the ZH
     * production cross-section with subsequent decay into @f$\tau\tau@f$ in the
     * current model and in the Standard Model.
     * @param[in] sqrt_s the center-of-mass energy in TeV
     * @return @f$\mu_{ZH,\tau\tau}@f$
     */
    virtual double muTHUZHtautau(const double sqrt_s) const;
    /**
     * @brief The ratio @f$\mu_{WH,\tau\tau}@f$ between the WH
     * production cross-section with subsequent decay into @f$\tau\tau@f$ in the
     * current model and in the Standard Model.
     * @param[in] sqrt_s the center-of-mass energy in TeV
     * @return @f$\mu_{WH,\tau\tau}@f$
     */
    virtual double muTHUWHtautau(const double sqrt_s) const;
    /**
     * @brief The ratio @f$\mu_{VH,\tau\tau}@f$ between the VH
     * production cross-section with subsequent decay into @f$\tau\tau@f$ in the
     * current model and in the Standard Model.
     * @param[in] sqrt_s the center-of-mass energy in TeV
     * @return @f$\mu_{VH,\tau\tau}@f$
     */
    virtual double muTHUVHtautau(const double sqrt_s) const;
    /**
     * @brief The ratio @f$\mu_{ttH,\tau\tau}@f$ between the ttH
     * production cross-section with subsequent decay into @f$\tau\tau@f$ in the
     * current model and in the Standard Model.
     * @param[in] sqrt_s the center-of-mass energy in TeV
     * @return @f$\mu_{ttH,\tau\tau}@f$
     */
    virtual double muTHUttHtautau(const double sqrt_s) const;
    /**
     * @brief The ratio @f$\mu_{ggH,bb}@f$ between the gluon-gluon fusion Higgs
     * production cross-section with subsequent decay into @f$bb@f$ in the
     * current model and in the Standard Model.
     * @param[in] sqrt_s the center-of-mass energy in TeV
     * @return @f$\mu_{ggH,bb}@f$
     */
    virtual double muTHUggHbb(const double sqrt_s) const;
    /**
     * @brief The ratio @f$\mu_{VBF,bb}@f$ between the VBF Higgs
     * production cross-section with subsequent decay into @f$bb@f$ in the
     * current model and in the Standard Model.
     * @param[in] sqrt_s the center-of-mass energy in TeV
     * @return @f$\mu_{VBF,bb}@f$
     */
    virtual double muTHUVBFHbb(const double sqrt_s) const;
    /**
     * @brief The ratio @f$\mu_{ZH,bb}@f$ between the ZH
     * production cross-section with subsequent decay into @f$bb@f$ in the
     * current model and in the Standard Model.
     * @param[in] sqrt_s the center-of-mass energy in TeV
     * @return @f$\mu_{ZH,bb}@f$
     */
    virtual double muTHUZHbb(const double sqrt_s) const;
    /**
     * @brief The ratio @f$\mu_{WH,bb}@f$ between the WH
     * production cross-section with subsequent decay into @f$bb@f$ in the
     * current model and in the Standard Model.
     * @param[in] sqrt_s the center-of-mass energy in TeV
     * @return @f$\mu_{WH,bb}@f$
     */
    virtual double muTHUWHbb(const double sqrt_s) const;
    /**
     * @brief The ratio @f$\mu_{VH,bb}@f$ between the VH
     * production cross-section with subsequent decay into @f$bb@f$ in the
     * current model and in the Standard Model.
     * @param[in] sqrt_s the center-of-mass energy in TeV
     * @return @f$\mu_{VH,bb}@f$
     */
    virtual double muTHUVHbb(const double sqrt_s) const;
    /**
     * @brief The ratio @f$\mu_{ttH,bb}@f$ between the ttH
     * production cross-section with subsequent decay into @f$bb@f$ in the
     * current model and in the Standard Model.
     * @param[in] sqrt_s the center-of-mass energy in TeV
     * @return @f$\mu_{ttH,bb}@f$
     */
    virtual double muTHUttHbb(const double sqrt_s) const;
    
    /**
     * @brief The ratio @f$\mu_{VBF}@f$ between the VBF
     * production cross-section in the
     * current model and in the Standard Model, multiplied by the 
     * total (SM+new physics) invisible decay branching ratio.
     * @param[in] sqrt_s the center-of-mass energy in TeV
     * @return @f$\mu_{VBF}BR_{inv}@f$
     */    
    virtual double muTHUVBFBRinv(const double sqrt_s) const;
    /**
     * @brief The ratio @f$\mu_{VBF,inv}@f$ between the VBF
     * production cross-section with subsequent decay into invisible states in the
     * current model and in the Standard Model.
     * @param[in] sqrt_s the center-of-mass energy in TeV
     * @return @f$\mu_{VBF,inv}@f$
     */     
    virtual double muTHUVBFHinv(const double sqrt_s) const;
    /**
     * @brief The ratio @f$\mu_{VH}@f$ between the VH
     * production cross-section in the
     * current model and in the Standard Model, multiplied by the 
     * total (SM+new physics) invisible decay branching ratio.
     * @param[in] sqrt_s the center-of-mass energy in TeV
     * @return @f$\mu_{VH}BR_{inv}@f$
     */     
    virtual double muTHUVHBRinv(const double sqrt_s) const;
    /**
     * @brief The ratio @f$\mu_{VH,inv}@f$ between the VH
     * production cross-section with subsequent decay into invisible states in the
     * current model and in the Standard Model.
     * @param[in] sqrt_s the center-of-mass energy in TeV
     * @return @f$\mu_{VH,inv}@f$
     */    
    virtual double muTHUVHinv(const double sqrt_s) const;
    
    /**
     * @brief The ratio @f$\mu_{ggH,ZZ\to 4\mu}@f$ between the gluon-gluon fusion Higgs
     * production cross-section with subsequent decay into @f$Z Z^*\to 4\mu@f$ in the
     * current model and in the Standard Model.
     * @param[in] sqrt_s the center-of-mass energy in TeV
     * @return @f$\mu_{ggH,ZZ\to 4\mu}@f$
     */    
    virtual double muTHUggHZZ4mu(const double sqrt_s) const;    
    /**
     * @brief The ratio @f$\mu_{ggH,Z\gamma\to \gamma 2\mu}@f$ between the gluon-gluon fusion Higgs
     * production cross-section with subsequent decay into @f$Z \gamma\to \gamma 2\mu@f$ in the
     * current model and in the Standard Model.
     * @param[in] sqrt_s the center-of-mass energy in TeV
     * @return @f$\mu_{ggH,Z\gamma\to \gamma 2\mu}@f$
     */    
    virtual double muTHUggHZgamumu(const double sqrt_s) const;
      
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
      
    /**
     * @brief The new physics contribution to the effective anomalous triple 
     * gauge coupling @f$g_{1,Z}^{Eff}@f$ from arXiv: 1708.09079 [hep-ph].
     * @return @f$\delta g_{1,Z}@f$
     */
    virtual double deltag1ZNPEff() const;
      
    /**
     * @brief The new physics contribution to the effective anomalous triple 
     * gauge coupling @f$\kappa_{\gamma}^{Eff}@f$ from arXiv: 1708.09079 [hep-ph].
     * @return @f$\delta \kappa_{\gamma}@f$
     */
    virtual double deltaKgammaNPEff() const;
    
    ////////////////////////////////////////////////////////////////////////
    
    /**
     * @brief The new physics contribution to the cross section in pb for @f$e^+ e^- \to W^+ W^- \to 4f @f$, 
     * with @f$ 4f = 0 (jjjj), 1 (e v jj), 2 (mu v jj), 3 (tau v jj),
     * 4 (e v e v), 5 (mu v mu v), 6 (tau v tau v),
     * 7 (e v mu v), 8 (e v tau v), 9 (mu v tau v), 10 (l v jj), 11 (l v l v) @f$
     * the different fermion final states for C.O.M. energies in 188-208 GeV.
     * From arXiv: 1606.06693 [hep-ph].
     * @return @f$\delta sigma@f$ [pb]
     */
    virtual double deltaxseeWW4fLEP2(const double sqrt_s, const int fstate) const;
    
    /**
     * @brief The cross section in pb for @f$e^+ e^- \to W^+ W^- \to 4f @f$, 
     * with @f$ 4f = 0 (jjjj), 1 (e v jj), 2 (mu v jj), 3 (tau v jj),
     * 4 (e v e v), 5 (mu v mu v), 6 (tau v tau v),
     * 7 (e v mu v), 8 (e v tau v), 9 (mu v tau v), 10 (l v jj), 11 (l v l v) @f$
     * the different fermion final states for C.O.M. energies in 188-208 GeV.
     * From arXiv: 1606.06693 [hep-ph].
     * @return @f$sigma@f$ [pb]
     */
    virtual double xseeWW4fLEP2(const double sqrt_s, const int fstate) const;

    /**
     * @brief The new physics contribution to the total cross section in pb for @f$e^+ e^- \to W^+ W^-@f$, 
     * summing over all final states for C.O.M. energies in 188-208 GeV.
     * From arXiv: 1606.06693 [hep-ph].
     * @return @f$\delta sigma@f$ [pb]
     */
    virtual double deltaxseeWWtotLEP2(const double sqrt_s) const;

    
    /**
     * @brief The total cross section in pb for @f$e^+ e^- \to W^+ W^-@f$, 
     * summing over all final states for C.O.M. energies in 188-208 GeV.
     * From arXiv: 1606.06693 [hep-ph].
     * @return @f$sigma@f$ [pb]
     */
    virtual double xseeWWtotLEP2(const double sqrt_s) const;
    
    /**
     * @brief The new physics contribution to the differential cross section in pb for @f$e^+ e^- \to W^+ W^- \to lv jj @f$, 
     * with @f$ l= e,\mu @f$ for the 4 @f$ cos{\theta}@f$ bins defined in arXiv: 1606.06693 [hep-ph].
     * for the C.O.M. energies of 182.6 and 205.9 GeV.
     * From arXiv: 1606.06693 [hep-ph].
     * @return @f$\delta d\sigma/d\cos{\theta}@f$ [pb]
     */
    virtual double deltadxsdcoseeWWlvjjLEP2(const double sqrt_s, const int bin) const;
    
    /**
     * @brief The differential cross section in pb for @f$e^+ e^- \to W^+ W^- \to lv jj @f$, 
     * with @f$ l= e,\mu @f$ for the 4 @f$ cos{\theta}@f$ bins defined in arXiv: 1606.06693 [hep-ph].
     * for the C.O.M. energies of 182.6 and 205.9 GeV.
     * From arXiv: 1606.06693 [hep-ph].
     * @return @f$d\sigma/d\cos{\theta}@f$ [pb]
     */
    virtual double dxsdcoseeWWlvjjLEP2(const double sqrt_s, const int bin) const;
    
    /**
     * @brief The differential distribution for @f$e^+ e^- \to W^+ W^- \to jj \ell \nu@f$, 
     * with @f$\ell= e, \mu@f$, as a function of the @f$W@f$ polar angle.
     * @return @f$d\sigma/d\cos{\theta}@f$
     */
    virtual double dxseeWWdcos(const double sqrt_s, const double cos) const;
    
    /**
     * @brief The integral of differential distribution for @f$e^+ e^- \to W^+ W^- \to jj \ell \nu@f$, 
     * with @f$\ell= e, \mu@f$ in a given bin of the @f$W@f$ polar angle.
     * @return @f$\int_{\cos{\theta_1}}^{\cos{\theta_2}} d\sigma/d\cos{\theta}@f$
     */
    virtual double dxseeWWdcosBin(const double sqrt_s, const double cos1, const double cos2) const;
    
    /**
     * @brief Total @f$e^+ e^- \to W^+ W^- \to jj \ell \nu@f$ cross section in pb, 
     * with @f$\ell= e, \mu@f$.
     * @return @f$\sigma(e^+ e^- \to W^+ W^- \to jj \ell \nu) @f$
     */
    virtual double xseeWW(const double sqrt_s) const;
    
    /**
     * @brief The ratio @f$\mu_{eeWW}@f$ between the 
     * @f$ e^{+}e^{-}\to W^{+}W^{-} @f$ production
     * cross-section in the current model and in the Standard Model.
     * @param[in] sqrt_s the center-of-mass energy in TeV
     * @return @f$\mu_{eeWW}@f$
     */
    virtual double mueeWW(const double sqrt_s) const;
    
    /**
     * @brief The ratio @f$\mu_{eeWW}@f$ between the 
     * @f$ e^{+}e^{-}\to W^{+}W^{-} @f$ production
     * cross-section in the current model and in the Standard Model.
     * @param[in] sqrt_s the center-of-mass energy in TeV, Pol_em and Pol_ep
     * are the polarization of electrons and positrons, respectively 
     * @return @f$\mu_{eeWW}@f$
     */
    virtual double mueeWWPol(const double sqrt_s, const double Pol_em, const double Pol_ep) const;
    
    ////////////////////////////////////////////////////////////////////////
    
    //----- High Energy diboson observables at hadron colliders

    /**
     * @brief The direction constrained by @f$ p p \to Z H@f$ in the boosted regime, @f$g_p^Z@f$.
     * From arXiv:1807.01796 and the contribution to FCC CDR Vol 1. Implemented only in NPSMEFTd6 class.
     * @return @f$g_p^Z@f$
     */
    virtual double ppZHprobe(const double sqrt_s) const;
    
    /**
     * @brief The number of events in  @f$ p p \to WZ@f$
     * in a given @f$p_{TV}@f$ bin, normalized to the SM prediction.
     * From arXiv: 1712.01310 [hep-ph] and private communication.
     * Implemented only in NPSMEFTd6 class.
     * @return @f$N_{ev}^{p_{TV}}/N_{ev,SM}^{p_{TV}}@f$
     */
    virtual double mupTVppWZ(const double sqrt_s, const double pTV1, const double pTV2) const;
    
    
    ////////////////////////////////////////////////////////////////////////
    
    //----- Simplified Template Cross Sections Bins
    
    //----- Stage 0
    
    /**
     * @brief The STXS0 bin @f$pp \to H qq@f$.
     * @param[in] sqrt_s the center-of-mass energy in TeV
     */
    virtual double STXS0_qqH(const double sqrt_s) const;
    
    
    //----- Stage 1
    
    /**
     * @brief The STXS bin @f$gg \to H@f$.
     * @param[in] sqrt_s the center-of-mass energy in TeV
     */
    virtual double STXS_ggH_VBFtopo_j3v(const double sqrt_s) const;
    
    
    /**
     * @brief The STXS bin @f$gg \to H@f$.
     * @param[in] sqrt_s the center-of-mass energy in TeV
     */
    virtual double STXS_ggH_VBFtopo_j3(const double sqrt_s) const;
    
    /**
     * @brief The STXS bin @f$gg \to H@f$.
     * @param[in] sqrt_s the center-of-mass energy in TeV
     */
    virtual double STXS_ggH0j(const double sqrt_s) const;
    
    /**
     * @brief The STXS bin @f$gg \to H@f$.
     * @param[in] sqrt_s the center-of-mass energy in TeV
     */
    virtual double STXS_ggH1j_pTH_0_60(const double sqrt_s) const;
    
    /**
     * @brief The STXS bin @f$gg \to H@f$.
     * @param[in] sqrt_s the center-of-mass energy in TeV
     */
    virtual double STXS_ggH1j_pTH_60_120(const double sqrt_s) const;
    
    /**
     * @brief The STXS bin @f$gg \to H@f$.
     * @param[in] sqrt_s the center-of-mass energy in TeV
     */
    virtual double STXS_ggH1j_pTH_120_200(const double sqrt_s) const;
    
    /**
     * @brief The STXS bin @f$gg \to H@f$.
     * @param[in] sqrt_s the center-of-mass energy in TeV
     */
    virtual double STXS_ggH1j_pTH_200(const double sqrt_s) const;
    
    /**
     * @brief The STXS bin @f$gg \to H@f$.
     * @param[in] sqrt_s the center-of-mass energy in TeV
     */
    virtual double STXS_ggH2j_pTH_0_200(const double sqrt_s) const;
    
    /**
     * @brief The STXS bin @f$gg \to H@f$.
     * @param[in] sqrt_s the center-of-mass energy in TeV
     */
    virtual double STXS_ggH2j_pTH_0_60(const double sqrt_s) const;
    
    
    /**
     * @brief The STXS bin @f$gg \to H@f$.
     * @param[in] sqrt_s the center-of-mass energy in TeV
     */
    virtual double STXS_ggH2j_pTH_60_120(const double sqrt_s) const;
    
        
    /**
     * @brief The STXS bin @f$gg \to H@f$.
     * @param[in] sqrt_s the center-of-mass energy in TeV
     */
    virtual double STXS_ggH2j_pTH_120_200(const double sqrt_s) const;
    
    
    /**
     * @brief The STXS bin @f$gg \to H@f$.
     * @param[in] sqrt_s the center-of-mass energy in TeV
     */
    virtual double STXS_ggH2j_pTH_200(const double sqrt_s) const;
    
    /**
     * @brief The STXS bin @f$qq \to H qq@f$.
     * @param[in] sqrt_s the center-of-mass energy in TeV
     */
    virtual double STXS_qqHqq_VBFtopo_Rest(const double sqrt_s) const;
    
    
    /**
     * @brief The STXS bin @f$qq \to H qq@f$.
     * @param[in] sqrt_s the center-of-mass energy in TeV
     */
    virtual double STXS_qqHqq_VBFtopo_j3v(const double sqrt_s) const;
    
    
    /**
     * @brief The STXS bin @f$qq \to H qq@f$.
     * @param[in] sqrt_s the center-of-mass energy in TeV
     */
    virtual double STXS_qqHqq_VBFtopo_j3(const double sqrt_s) const;
    
    /**
     * @brief The STXS bin @f$qq \to H qq@f$.
     * @param[in] sqrt_s the center-of-mass energy in TeV
     */
    virtual double STXS_qqHqq_nonVHtopo(const double sqrt_s) const;
    
    /**
     * @brief The STXS bin @f$qq \to H qq@f$.
     * @param[in] sqrt_s the center-of-mass energy in TeV
     */
    virtual double STXS_qqHqq_VHtopo(const double sqrt_s) const;
    
    /**
     * @brief The STXS bin @f$qq \to H qq@f$.
     * @param[in] sqrt_s the center-of-mass energy in TeV
     */
    virtual double STXS_qqHqq_Rest(const double sqrt_s) const;
    
    /**
     * @brief The STXS bin @f$qq \to H qq@f$.
     * @param[in] sqrt_s the center-of-mass energy in TeV
     */
    virtual double STXS_qqHqq_pTj_200(const double sqrt_s) const;
    
    /**
     * @brief The STXS bin @f$qq \to H \ell \nu@f$.
     * @param[in] sqrt_s the center-of-mass energy in TeV
     */
    virtual double STXS_qqHlv_pTV_0_250(const double sqrt_s) const;
    
 
    /**
     * @brief The STXS bin @f$qq \to H \ell \nu@f$.
     * @param[in] sqrt_s the center-of-mass energy in TeV
     */
    virtual double STXS_qqHlv_pTV_0_150(const double sqrt_s) const;
    
    
    /**
     * @brief The STXS bin @f$qq \to H \ell \nu@f$.
     * @param[in] sqrt_s the center-of-mass energy in TeV
     */
    virtual double STXS_qqHlv_pTV_150_250_0j(const double sqrt_s) const;
    
    
    /**
     * @brief The STXS bin @f$qq \to H \ell \nu@f$.
     * @param[in] sqrt_s the center-of-mass energy in TeV
     */
    virtual double STXS_qqHlv_pTV_150_250_1j(const double sqrt_s) const;
    
    
    /**
     * @brief The STXS bin @f$qq \to H \ell \nu@f$.
     * @param[in] sqrt_s the center-of-mass energy in TeV
     */
    virtual double STXS_qqHlv_pTV_250(const double sqrt_s) const;
    
    /**
     * @brief The STXS bin @f$qq \to H \ell \ell@f$.
     * @param[in] sqrt_s the center-of-mass energy in TeV
     */
    virtual double STXS_qqHll_pTV_0_150(const double sqrt_s) const;
    
    /**
     * @brief The STXS bin @f$qq \to H \ell \ell@f$.
     * @param[in] sqrt_s the center-of-mass energy in TeV
     */
    virtual double STXS_qqHll_pTV_150_250(const double sqrt_s) const;
    
    
    /**
     * @brief The STXS bin @f$qq \to H \ell \ell@f$.
     * @param[in] sqrt_s the center-of-mass energy in TeV
     */
    virtual double STXS_qqHll_pTV_150_250_0j(const double sqrt_s) const;
    
    
    /**
     * @brief The STXS bin @f$qq \to H \ell \ell@f$.
     * @param[in] sqrt_s the center-of-mass energy in TeV
     */
    virtual double STXS_qqHll_pTV_150_250_1j(const double sqrt_s) const;
    
    /**
     * @brief The STXS bin @f$qq \to H \ell \ell@f$.
     * @param[in] sqrt_s the center-of-mass energy in TeV
     */
    virtual double STXS_qqHll_pTV_250(const double sqrt_s) const;
    
    
    /**
     * @brief The STXS bin @f$ ttH + tH @f$.
     * @param[in] sqrt_s the center-of-mass energy in TeV
     */
    virtual double STXS_ttHtH(const double sqrt_s) const;
    
    /**
     * @brief The STXS bin @f$ qq \to WH \to H qq @f$.
     * @param[in] sqrt_s the center-of-mass energy in TeV
     */
    virtual double STXS_WHqqHqq_VBFtopo_j3v(const double sqrt_s) const;
    
    /**
     * @brief The STXS bin @f$ qq \to WH \to H qq @f$.
     * @param[in] sqrt_s the center-of-mass energy in TeV
     */
    virtual double STXS_WHqqHqq_VBFtopo_j3(const double sqrt_s) const;
    
    /**
     * @brief The STXS bin @f$ qq \to WH \to H qq @f$.
     * @param[in] sqrt_s the center-of-mass energy in TeV
     */
    virtual double STXS_WHqqHqq_VH2j(const double sqrt_s) const;
    
    /**
     * @brief The STXS bin @f$ qq \to WH \to H qq @f$.
     * @param[in] sqrt_s the center-of-mass energy in TeV
     */
    virtual double STXS_WHqqHqq_Rest(const double sqrt_s) const;
    
    /**
     * @brief The STXS bin @f$ qq \to WH \to H qq @f$.
     * @param[in] sqrt_s the center-of-mass energy in TeV
     */
    virtual double STXS_WHqqHqq_pTj1_200(const double sqrt_s) const;
    
    /**
     * @brief The STXS bin @f$ qq \to ZH \to H qq @f$.
     * @param[in] sqrt_s the center-of-mass energy in TeV
     */
    virtual double STXS_ZHqqHqq_VBFtopo_j3v(const double sqrt_s) const;
    
    /**
     * @brief The STXS bin @f$ qq \to ZH \to H qq @f$.
     * @param[in] sqrt_s the center-of-mass energy in TeV
     */
    virtual double STXS_ZHqqHqq_VBFtopo_j3(const double sqrt_s) const;
    
    /**
     * @brief The STXS bin @f$ qq \to ZH \to H qq @f$.
     * @param[in] sqrt_s the center-of-mass energy in TeV
     */
    virtual double STXS_ZHqqHqq_VH2j(const double sqrt_s) const;
    
    /**
     * @brief The STXS bin @f$ qq \to ZH \to H qq @f$.
     * @param[in] sqrt_s the center-of-mass energy in TeV
     */
    virtual double STXS_ZHqqHqq_Rest(const double sqrt_s) const;
    
    /**
     * @brief The STXS bin @f$ qq \to ZH \to H qq @f$.
     * @param[in] sqrt_s the center-of-mass energy in TeV
     */
    virtual double STXS_ZHqqHqq_pTj1_200(const double sqrt_s) const;
    
    
    ////////////////////////////////////////////////////////////////////////
    
    /**
     * @brief The effective coupling @f$\kappa_{\mu,eff}=\sqrt{\Gamma_{H\mu\mu}/\Gamma_{H\mu\mu}^{SM}}@f$.
     * @return @f$\kappa_{\mu,eff}@f$
     */
    virtual double kappamueff() const;
    
    /**
     * @brief The effective coupling @f$\kappa_{\tau,eff}=\sqrt{\Gamma_{H\tau\tau}/\Gamma_{H\tau\tau}^{SM}}@f$.
     * @return @f$\kappa_{\tau,eff}@f$
     */
    virtual double kappataueff() const;
    
    /**
     * @brief The effective coupling @f$\kappa_{c,eff}=\sqrt{\Gamma_{Hcc}/\Gamma_{Hcc}^{SM}}@f$.
     * @return @f$\kappa_{c,eff}@f$
     */
    virtual double kappaceff() const;
    
    /**
     * @brief The effective coupling @f$\kappa_{b,eff}=\sqrt{\Gamma_{Hbb}/\Gamma_{Hbb}^{SM}}@f$.
     * @return @f$\kappa_{b,eff}@f$
     */
    virtual double kappabeff() const;
    
    /**
     * @brief The effective coupling @f$\kappa_{G,eff}=\sqrt{\Gamma_{HGG}/\Gamma_{HGG}^{SM}}@f$.
     * @return @f$\kappa_{G,eff}@f$
     */
    virtual double kappaGeff() const;
    
    /**
     * @brief The effective coupling @f$\kappa_{Z,eff}=\sqrt{\Gamma_{HZZ}/\Gamma_{HZZ}^{SM}}@f$.
     * @return @f$\kappa_{Z,eff}@f$
     */
    virtual double kappaZeff() const;
    
    /**
     * @brief The effective coupling @f$\kappa_{W,eff}=\sqrt{\Gamma_{HWW}/\Gamma_{HWW}^{SM}}@f$.
     * @return @f$\kappa_{W,eff}@f$
     */
    virtual double kappaWeff() const;
    
    /**
     * @brief The effective coupling @f$\kappa_{A,eff}=\sqrt{\Gamma_{HAA}/\Gamma_{HAA}^{SM}}@f$.
     * @return @f$\kappa_{A,eff}@f$
     */
    virtual double kappaAeff() const;
    
    /**
     * @brief The effective coupling @f$\kappa_{ZA,eff}=\sqrt{\Gamma_{HZA}/\Gamma_{HZA}^{SM}}@f$.
     * @return @f$\kappa_{ZA,eff}@f$
     */
    virtual double kappaZAeff() const;
    
    /////////////Basic interactions of the so-called Higgs basis////////////////
    
    /**
     * @brief The Higgs-basis coupling @f$\delta y_t@f$.
     * (See LHCHXSWG-INT-2015-001 document.)
     * Note that the Lagrangian definition of the Higgs-basis parameters coincides with the one of 
     * some of the @f$g_i, \delta g_i@f$ couplings defined above.
     * In the Higgs basis, however, one uses the freedom to perform certain field redefinitions and 
     * operations to demand that the mass eigenstate Lagrangian has specific features. (See pag. 5,6 in the reference.)
     * Therefore, the actual expression in terms of dim 6 coefficients may differ from the one for @f$g_i, \delta g_i@f$.
     * @return @f$\delta y_t@f$
     */
    virtual double deltayt_HB() const;
    
    /**
     * @brief The Higgs-basis coupling @f$\delta y_b@f$.
     * (See LHCHXSWG-INT-2015-001 document.)
     * Note that the Lagrangian definition of the Higgs-basis parameters coincides with the one of 
     * some of the @f$g_i, \delta g_i@f$ couplings defined above.
     * In the Higgs basis, however, one uses the freedom to perform certain field redefinitions and 
     * operations to demand that the mass eigenstate Lagrangian has specific features. (See pag. 5,6 in the reference.)
     * Therefore, the actual expression in terms of dim 6 coefficients may differ from the one for @f$g_i, \delta g_i@f$.
     * @return @f$\delta y_b@f$
     */
    virtual double deltayb_HB() const;
    
    /**
     * @brief The Higgs-basis coupling @f$\delta y_\tau@f$.
     * (See LHCHXSWG-INT-2015-001 document.)
     * Note that the Lagrangian definition of the Higgs-basis parameters coincides with the one of 
     * some of the @f$g_i, \delta g_i@f$ couplings defined above.
     * In the Higgs basis, however, one uses the freedom to perform certain field redefinitions and 
     * operations to demand that the mass eigenstate Lagrangian has specific features. (See pag. 5,6 in the reference.)
     * Therefore, the actual expression in terms of dim 6 coefficients may differ from the one for @f$g_i, \delta g_i@f$.
     * @return @f$\delta y_\tau@f$
     */
    virtual double deltaytau_HB() const;
    
    /**
     * @brief The Higgs-basis coupling @f$\delta y_c@f$.
     * (See LHCHXSWG-INT-2015-001 document.)
     * Note that the Lagrangian definition of the Higgs-basis parameters coincides with the one of 
     * some of the @f$g_i, \delta g_i@f$ couplings defined above.
     * In the Higgs basis, however, one uses the freedom to perform certain field redefinitions and 
     * operations to demand that the mass eigenstate Lagrangian has specific features. (See pag. 5,6 in the reference.)
     * Therefore, the actual expression in terms of dim 6 coefficients may differ from the one for @f$g_i, \delta g_i@f$.
     * @return @f$\delta y_c@f$
     */
    virtual double deltayc_HB() const;
    
    /**
     * @brief The Higgs-basis coupling @f$\delta y_\mu@f$.
     * (See LHCHXSWG-INT-2015-001 document.) 
     * Note that the Lagrangian definition of the Higgs-basis parameters coincides with the one of 
     * some of the @f$g_i, \delta g_i@f$ couplings defined above.
     * In the Higgs basis, however, one uses the freedom to perform certain field redefinitions and 
     * operations to demand that the mass eigenstate Lagrangian has specific features. (See pag. 5,6 in the reference.)
     * Therefore, the actual expression in terms of dim 6 coefficients may differ from the one for @f$g_i, \delta g_i@f$.
     * @return @f$\delta y_\mu@f$
     */
    virtual double deltaymu_HB() const;
    
    /**
     * @brief The Higgs-basis coupling @f$\delta c_z@f$.
     * (See LHCHXSWG-INT-2015-001 document.)
     * Note that the Lagrangian definition of the Higgs-basis parameters coincides with the one of 
     * some of the @f$g_i, \delta g_i@f$ couplings defined above.
     * In the Higgs basis, however, one uses the freedom to perform certain field redefinitions and 
     * operations to demand that the mass eigenstate Lagrangian has specific features. (See pag. 5,6 in the reference.)
     * Therefore, the actual expression in terms of dim 6 coefficients may differ from the one for @f$g_i, \delta g_i@f$.
     * @return @f$\delta c_z@f$
     */
    virtual double deltacZ_HB() const;
    
    /**
     * @brief The Higgs-basis coupling @f$c_{z\Box}@f$.
     * (See LHCHXSWG-INT-2015-001 document.)
     * Note that the Lagrangian definition of the Higgs-basis parameters coincides with the one of 
     * some of the @f$g_i, \delta g_i@f$ couplings defined above.
     * In the Higgs basis, however, one uses the freedom to perform certain field redefinitions and 
     * operations to demand that the mass eigenstate Lagrangian has specific features. (See pag. 5,6 in the reference.)
     * Therefore, the actual expression in terms of dim 6 coefficients may differ from the one for @f$g_i, \delta g_i@f$.
     * @return @f$c_{z\Box}@f$
     */
    virtual double cZBox_HB() const;
    
    /**
     * @brief The Higgs-basis coupling @f$c_{zz}@f$.
     * (See LHCHXSWG-INT-2015-001 document.)
     * Note that the Lagrangian definition of the Higgs-basis parameters coincides with the one of 
     * some of the @f$g_i, \delta g_i@f$ couplings defined above.
     * In the Higgs basis, however, one uses the freedom to perform certain field redefinitions and 
     * operations to demand that the mass eigenstate Lagrangian has specific features. (See pag. 5,6 in the reference.)
     * Therefore, the actual expression in terms of dim 6 coefficients may differ from the one for @f$g_i, \delta g_i@f$.
     * @return @f$c_{zz}@f$
     */
    virtual double cZZ_HB() const;
    
    /**
     * @brief The Higgs-basis coupling @f$c_{z\gamma}@f$.
     * (See LHCHXSWG-INT-2015-001 document.)
     * Note that the Lagrangian definition of the Higgs-basis parameters coincides with the one of 
     * some of the @f$g_i, \delta g_i@f$ couplings defined above.
     * In the Higgs basis, however, one uses the freedom to perform certain field redefinitions and 
     * operations to demand that the mass eigenstate Lagrangian has specific features. (See pag. 5,6 in the reference.)
     * Therefore, the actual expression in terms of dim 6 coefficients may differ from the one for @f$g_i, \delta g_i@f$.
     * @return @f$c_{z\gamma}@f$
     */
    virtual double cZga_HB() const;
    
    /**
     * @brief The Higgs-basis coupling @f$c_{\gamma\gamma}@f$.
     * (See LHCHXSWG-INT-2015-001 document.)
     * Note that the Lagrangian definition of the Higgs-basis parameters coincides with the one of 
     * some of the @f$g_i, \delta g_i@f$ couplings defined above.
     * In the Higgs basis, however, one uses the freedom to perform certain field redefinitions and 
     * operations to demand that the mass eigenstate Lagrangian has specific features. (See pag. 5,6 in the reference.)
     * Therefore, the actual expression in terms of dim 6 coefficients may differ from the one for @f$g_i, \delta g_i@f$.
     * @return @f$c_{\gamma\gamma}@f$
     */
    virtual double cgaga_HB() const;
    
    /**
     * @brief The Higgs-basis coupling @f$c_{gg}@f$.
     * (See LHCHXSWG-INT-2015-001 document.)
     * Note that the Lagrangian definition of the Higgs-basis parameters coincides with the one of 
     * some of the @f$g_i, \delta g_i@f$ couplings defined above.
     * In the Higgs basis, however, one uses the freedom to perform certain field redefinitions and 
     * operations to demand that the mass eigenstate Lagrangian has specific features. (See pag. 5,6 in the reference.)
     * Therefore, the actual expression in terms of dim 6 coefficients may differ from the one for @f$g_i, \delta g_i@f$.
     * @return @f$c_{gg}@f$
     */
    virtual double cgg_HB() const;
    
    /**
     * @brief The effective Higgs-basis coupling @f$c_{gg}^{Eff}@f$. (Similar to cgg_HB but including modifications of SM loops.)
     * (See arXiv: 1505.00046 [hep-ph] document.)
     * Note that the Lagrangian definition of the Higgs-basis parameters coincides with the one of 
     * some of the @f$g_i, \delta g_i@f$ couplings defined above.
     * In the Higgs basis, however, one uses the freedom to perform certain field redefinitions and 
     * operations to demand that the mass eigenstate Lagrangian has specific features. (See pag. 5,6 in the reference.)
     * Therefore, the actual expression in terms of dim 6 coefficients may differ from the one for @f$g_i, \delta g_i@f$.
     * @return @f$c_{gg}^{Eff}@f$
     */
    virtual double cggEff_HB() const;
    
    /**
     * @brief The Higgs-basis coupling @f$\lambda_{z}@f$.
     * (See LHCHXSWG-INT-2015-001 document.)
     * Note that the Lagrangian definition of the Higgs-basis parameters coincides with the one of 
     * some of the @f$g_i, \delta g_i@f$ couplings defined above.
     * In the Higgs basis, however, one uses the freedom to perform certain field redefinitions and 
     * operations to demand that the mass eigenstate Lagrangian has specific features. (See pag. 5,6 in the reference.)
     * Therefore, the actual expression in terms of dim 6 coefficients may differ from the one for @f$g_i, \delta g_i@f$.
     * @return @f$\lambda_{z}@f$
     */
    virtual double lambz_HB() const;
    
    /////////////Combinations of Warsaw basis coefficients constrained by EWPO////////////////
    
    /**
     * @brief Combination of coefficients of the Warsaw basis constrained by EWPO
     * @f$(\hat{C}_{HL}^{(1)})_{11}@f$.
     * @return @f$(\hat{C}_{HL}^{(1)})_{11}@f$
     */
    virtual double CEWHL111() const;
    
    /**
     * @brief Combination of coefficients of the Warsaw basis constrained by EWPO
     * @f$(\hat{C}_{HL}^{(1)})_{22}@f$.
     * @return @f$(\hat{C}_{HL}^{(1)})_{22}@f$
     */
    virtual double CEWHL122() const;
    
    /**
     * @brief Combination of coefficients of the Warsaw basis constrained by EWPO
     * @f$(\hat{C}_{HL}^{(1)})_{33}@f$.
     * @return @f$(\hat{C}_{HL}^{(1)})_{33}@f$
     */
    virtual double CEWHL133() const;
    
    /**
     * @brief Combination of coefficients of the Warsaw basis constrained by EWPO
     * @f$(\hat{C}_{HL}^{(3)})_{11}@f$.
     * @return @f$(\hat{C}_{HL}^{(3)})_{11}@f$
     */
    virtual double CEWHL311() const;
    
    /**
     * @brief Combination of coefficients of the Warsaw basis constrained by EWPO
     * @f$(\hat{C}_{HL}^{(3)})_{22}@f$.
     * @return @f$(\hat{C}_{HL}^{(3)})_{22}@f$
     */
    virtual double CEWHL322() const;
    
    /**
     * @brief Combination of coefficients of the Warsaw basis constrained by EWPO
     * @f$(\hat{C}_{HL}^{(3)})_{33}@f$.
     * @return @f$(\hat{C}_{HL}^{(3)})_{33}@f$
     */
    virtual double CEWHL333() const;
    
    /**
     * @brief Combination of coefficients of the Warsaw basis constrained by EWPO
     * @f$(\hat{C}_{HQ}^{(1)})_{11}@f$.
     * @return @f$(\hat{C}_{HQ}^{(1)})_{11}@f$
     */
    virtual double CEWHQ111() const;
    
    /**
     * @brief Combination of coefficients of the Warsaw basis constrained by EWPO
     * @f$(\hat{C}_{HQ}^{(1)})_{22}@f$.
     * @return @f$(\hat{C}_{HQ}^{(1)})_{22}@f$
     */
    virtual double CEWHQ122() const;
    
    /**
     * @brief Combination of coefficients of the Warsaw basis constrained by EWPO
     * @f$(\hat{C}_{HQ}^{(1)})_{33}@f$.
     * @return @f$(\hat{C}_{HQ}^{(1)})_{33}@f$
     */
    virtual double CEWHQ133() const;
    
    /**
     * @brief Combination of coefficients of the Warsaw basis constrained by EWPO
     * @f$(\hat{C}_{HQ}^{(3)})_{11}@f$.
     * @return @f$(\hat{C}_{HQ}^{(3)})_{11}@f$
     */
    virtual double CEWHQ311() const;
    
    /**
     * @brief Combination of coefficients of the Warsaw basis constrained by EWPO
     * @f$(\hat{C}_{HQ}^{(3)})_{22}@f$.
     * @return @f$(\hat{C}_{HQ}^{(3)})_{22}@f$
     */
    virtual double CEWHQ322() const;
    
    /**
     * @brief Combination of coefficients of the Warsaw basis constrained by EWPO
     * @f$(\hat{C}_{HQ}^{(3)})_{33}@f$.
     * @return @f$(\hat{C}_{HQ}^{(3)})_{33}@f$
     */
    virtual double CEWHQ333() const;
    
    
    /**
     * @brief Combination of coefficients of the Warsaw basis constrained by EWPO
     * @f$(\hat{C}_{HQ}^{(d)})_{33}@f$.
     * @return @f$(\hat{C}_{HQ}^{(d)})_{33}@f$
     */
    virtual double CEWHQd33() const;
    
        
    /**
     * @brief Combination of coefficients of the Warsaw basis constrained by EWPO
     * @f$(\hat{C}_{He})_{11}@f$.
     * @return @f$(\hat{C}_{He})_{11}@f$
     */
    virtual double CEWHe11() const;
    
    /**
     * @brief Combination of coefficients of the Warsaw basis constrained by EWPO
     * @f$(\hat{C}_{He})_{22}@f$.
     * @return @f$(\hat{C}_{He})_{22}@f$
     */
    virtual double CEWHe22() const;
    
    /**
     * @brief Combination of coefficients of the Warsaw basis constrained by EWPO
     * @f$(\hat{C}_{He})_{33}@f$.
     * @return @f$(\hat{C}_{He})_{33}@f$
     */
    virtual double CEWHe33() const;
    
    
    /**
     * @brief Combination of coefficients of the Warsaw basis constrained by EWPO
     * @f$(\hat{C}_{Hu})_{11}@f$.
     * @return @f$(\hat{C}_{Hu})_{11}@f$
     */
    virtual double CEWHu11() const;
    
    /**
     * @brief Combination of coefficients of the Warsaw basis constrained by EWPO
     * @f$(\hat{C}_{Hu})_{22}@f$.
     * @return @f$(\hat{C}_{Hu})_{22}@f$
     */
    virtual double CEWHu22() const;
    
    /**
     * @brief Combination of coefficients of the Warsaw basis constrained by EWPO
     * @f$(\hat{C}_{Hu})_{33}@f$.
     * @return @f$(\hat{C}_{Hu})_{33}@f$
     */
    virtual double CEWHu33() const;
    

    /**
     * @brief Combination of coefficients of the Warsaw basis constrained by EWPO
     * @f$(\hat{C}_{Hd})_{11}@f$.
     * @return @f$(\hat{C}_{Hd})_{11}@f$
     */
    virtual double CEWHd11() const;
    
    /**
     * @brief Combination of coefficients of the Warsaw basis constrained by EWPO
     * @f$(\hat{C}_{Hd})_{22}@f$.
     * @return @f$(\hat{C}_{Hd})_{22}@f$
     */
    virtual double CEWHd22() const;
    
    /**
     * @brief Combination of coefficients of the Warsaw basis constrained by EWPO
     * @f$(\hat{C}_{Hd})_{33}@f$.
     * @return @f$(\hat{C}_{Hd})_{33}@f$
     */
    virtual double CEWHd33() const;

    
    /////////////Auxiliary observables////////////////
    
    /**
     * @brief Auxiliary observable AuxObs_NP1
     * (See code for details.)
     * @return AuxObs_NP1
     */
    virtual double AuxObs_NP1() const;
    
    /**
     * @brief Auxiliary observable AuxObs_NP2
     * (See code for details.)
     * @return AuxObs_NP2
     */
    virtual double AuxObs_NP2() const;
    
    /**
     * @brief Auxiliary observable AuxObs_NP3
     * (See code for details.)
     * @return AuxObs_NP3
     */
    virtual double AuxObs_NP3() const;
    
    /**
     * @brief Auxiliary observable AuxObs_NP4
     * (See code for details.)
     * @return AuxObs_NP4
     */
    virtual double AuxObs_NP4() const;
    
    /**
     * @brief Auxiliary observable AuxObs_NP5
     * (See code for details.)
     * @return AuxObs_NP5
     */
    virtual double AuxObs_NP5() const;
    
    /**
     * @brief Auxiliary observable AuxObs_NP6
     * (See code for details.)
     * @return AuxObs_NP6
     */
    virtual double AuxObs_NP6() const;
    
    /**
     * @brief Auxiliary observable AuxObs_NP7
     * (See code for details.)
     * @return AuxObs_NP7
     */
    virtual double AuxObs_NP7() const;
    
    /**
     * @brief Auxiliary observable AuxObs_NP8
     * (See code for details.)
     * @return AuxObs_NP8
     */
    virtual double AuxObs_NP8() const;
    
    /**
     * @brief Auxiliary observable AuxObs_NP9
     * (See code for details.)
     * @return AuxObs_NP9
     */
    virtual double AuxObs_NP9() const;
    
    /**
     * @brief Auxiliary observable AuxObs_NP10
     * (See code for details.)
     * @return AuxObs_NP10
     */
    virtual double AuxObs_NP10() const;
    
    /**
     * @brief Auxiliary observable AuxObs_NP11
     * (See code for details.)
     * @return AuxObs_NP11
     */
    virtual double AuxObs_NP11() const;
        
    /**
     * @brief Auxiliary observable AuxObs_NP12
     * (See code for details.)
     * @return AuxObs_NP12
     */
    virtual double AuxObs_NP12() const;
    
    /**
     * @brief Auxiliary observable AuxObs_NP13
     * @return AuxObs_NP13
     */
    virtual double AuxObs_NP13() const;
    
    /**
     * @brief Auxiliary observable AuxObs_NP14
     * @return AuxObs_NP14
     */
    virtual double AuxObs_NP14() const;
    
    /**
     * @brief Auxiliary observable AuxObs_NP15
     * @return AuxObs_NP15
     */
    virtual double AuxObs_NP15() const;
    
    /**
     * @brief Auxiliary observable AuxObs_NP16
     * @return AuxObs_NP16
     */
    virtual double AuxObs_NP16() const;
    
    /**
     * @brief Auxiliary observable AuxObs_NP17
     * @return AuxObs_NP17
     */
    virtual double AuxObs_NP17() const;
    
    /**
     * @brief Auxiliary observable AuxObs_NP18
     * @return AuxObs_NP18
     */
    virtual double AuxObs_NP18() const;
    
    /**
     * @brief Auxiliary observable AuxObs_NP19
     * @return AuxObs_NP19
     */
    virtual double AuxObs_NP19() const;
    
    /**
     * @brief Auxiliary observable AuxObs_NP20
     * @return AuxObs_NP20
     */
    virtual double AuxObs_NP20() const;
    
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
    double C2B; ///< The dimension-6 operator coefficient \f$C_{2W}\f$.
    double C2W; ///< The dimension-6 operator coefficient \f$C_{2B}\f$.
    double C2BS; ///< The dimension-6 operator coefficient \f$C_{2W}^{SILH}\f$.
    double C2WS; ///< The dimension-6 operator coefficient \f$C_{2B}^{SILH}\f$.
    double CHG; ///< The dimension-6 operator coefficient \f$C_{HG}\f$.
    double CHW; ///< The dimension-6 operator coefficient \f$C_{HW}\f$.
    double CHB; ///< The dimension-6 operator coefficient \f$C_{HB}\f$.
    double CHWHB_gaga; ///< The combination of dimension-6 operator coefficients entering in \f$\delta_{AA}\f$: \f$s_W^2 C_{HW} + c_W^2 C_{HW}\f$.
    double CHWHB_gagaorth; ///< The combination of dimension-6 operator coefficients \f$-c_W^2 C_{HW} + s_W^2 C_{HW}\f$.
    double CDHB; ///< The dimension-6 operator coefficient \f$C_{DHB}\f$.
    double CDHW; ///< The dimension-6 operator coefficient \f$C_{DHW}\f$.
    double CDB; ///< The dimension-6 operator coefficient \f$C_{DB}\f$.
    double CDW; ///< The dimension-6 operator coefficient \f$C_{DW}\f$.
    double CHWB; ///< The dimension-6 operator coefficient \f$C_{HWB}\f$.
    double CHD; ///< The dimension-6 operator coefficient \f$C_{HD}\f$.
    double CT; ///< The dimension-6 operator coefficient \f$C_{T}\f$.
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
    double CdG_11r; ///< The dimension-6 operator coefficient \f$(C_{dG})_{11}\f$ (real part).
    double CdG_12r; ///< The dimension-6 operator coefficient \f$(C_{dG})_{12}\f$ (real part).
    double CdG_13r; ///< The dimension-6 operator coefficient \f$(C_{dG})_{13}\f$ (real part).
    double CdG_22r; ///< The dimension-6 operator coefficient \f$(C_{dG})_{22}\f$ (real part).
    double CdG_23r; ///< The dimension-6 operator coefficient \f$(C_{dG})_{23}\f$ (real part).
    double CdG_33r; ///< The dimension-6 operator coefficient \f$(C_{dG})_{33}\f$ (real part).
    double CdG_11i; ///< The dimension-6 operator coefficient \f$(C_{dG})_{11}\f$ (imaginary part).
    double CdG_12i; ///< The dimension-6 operator coefficient \f$(C_{dG})_{12}\f$ (imaginary part).
    double CdG_13i; ///< The dimension-6 operator coefficient \f$(C_{dG})_{13}\f$ (imaginary part).
    double CdG_22i; ///< The dimension-6 operator coefficient \f$(C_{dG})_{22}\f$ (imaginary part).
    double CdG_23i; ///< The dimension-6 operator coefficient \f$(C_{dG})_{23}\f$ (imaginary part).
    double CdG_33i; ///< The dimension-6 operator coefficient \f$(C_{dG})_{33}\f$ (imaginary part).
    double CdW_11r; ///< The dimension-6 operator coefficient \f$(C_{dW})_{11}\f$ (real part).
    double CdW_12r; ///< The dimension-6 operator coefficient \f$(C_{dW})_{12}\f$ (real part).
    double CdW_13r; ///< The dimension-6 operator coefficient \f$(C_{dW})_{13}\f$ (real part).
    double CdW_22r; ///< The dimension-6 operator coefficient \f$(C_{dW})_{22}\f$ (real part).
    double CdW_23r; ///< The dimension-6 operator coefficient \f$(C_{dW})_{23}\f$ (real part).
    double CdW_33r; ///< The dimension-6 operator coefficient \f$(C_{dW})_{33}\f$ (real part).
    double CdW_11i; ///< The dimension-6 operator coefficient \f$(C_{dW})_{11}\f$ (imaginary part).
    double CdW_12i; ///< The dimension-6 operator coefficient \f$(C_{dW})_{12}\f$ (imaginary part).
    double CdW_13i; ///< The dimension-6 operator coefficient \f$(C_{dW})_{13}\f$ (imaginary part).
    double CdW_22i; ///< The dimension-6 operator coefficient \f$(C_{dW})_{22}\f$ (imaginary part).
    double CdW_23i; ///< The dimension-6 operator coefficient \f$(C_{dW})_{23}\f$ (imaginary part).
    double CdW_33i; ///< The dimension-6 operator coefficient \f$(C_{dW})_{33}\f$ (imaginary part).
    double CdB_11r; ///< The dimension-6 operator coefficient \f$(C_{dB})_{11}\f$ (real part).
    double CdB_12r; ///< The dimension-6 operator coefficient \f$(C_{dB})_{12}\f$ (real part).
    double CdB_13r; ///< The dimension-6 operator coefficient \f$(C_{dB})_{13}\f$ (real part).
    double CdB_22r; ///< The dimension-6 operator coefficient \f$(C_{dB})_{22}\f$ (real part).
    double CdB_23r; ///< The dimension-6 operator coefficient \f$(C_{dB})_{23}\f$ (real part).
    double CdB_33r; ///< The dimension-6 operator coefficient \f$(C_{dB})_{33}\f$ (real part).
    double CdB_11i; ///< The dimension-6 operator coefficient \f$(C_{dB})_{11}\f$ (imaginary part).
    double CdB_12i; ///< The dimension-6 operator coefficient \f$(C_{dB})_{12}\f$ (imaginary part).
    double CdB_13i; ///< The dimension-6 operator coefficient \f$(C_{dB})_{13}\f$ (imaginary part).
    double CdB_22i; ///< The dimension-6 operator coefficient \f$(C_{dB})_{22}\f$ (imaginary part).
    double CdB_23i; ///< The dimension-6 operator coefficient \f$(C_{dB})_{23}\f$ (imaginary part).
    double CdB_33i; ///< The dimension-6 operator coefficient \f$(C_{dB})_{33}\f$ (imaginary part).
    double CeW_11r; ///< The dimension-6 operator coefficient \f$(C_{eW})_{11}\f$ (real part).
    double CeW_12r; ///< The dimension-6 operator coefficient \f$(C_{eW})_{12}\f$ (real part).
    double CeW_13r; ///< The dimension-6 operator coefficient \f$(C_{eW})_{13}\f$ (real part).
    double CeW_22r; ///< The dimension-6 operator coefficient \f$(C_{eW})_{22}\f$ (real part).
    double CeW_23r; ///< The dimension-6 operator coefficient \f$(C_{eW})_{23}\f$ (real part).
    double CeW_33r; ///< The dimension-6 operator coefficient \f$(C_{eW})_{33}\f$ (real part).
    double CeW_11i; ///< The dimension-6 operator coefficient \f$(C_{eW})_{11}\f$ (imaginary part).
    double CeW_12i; ///< The dimension-6 operator coefficient \f$(C_{eW})_{12}\f$ (imaginary part).
    double CeW_13i; ///< The dimension-6 operator coefficient \f$(C_{eW})_{13}\f$ (imaginary part).
    double CeW_22i; ///< The dimension-6 operator coefficient \f$(C_{eW})_{22}\f$ (imaginary part).
    double CeW_23i; ///< The dimension-6 operator coefficient \f$(C_{eW})_{23}\f$ (imaginary part).
    double CeW_33i; ///< The dimension-6 operator coefficient \f$(C_{eW})_{33}\f$ (imaginary part).
    double CeB_11r; ///< The dimension-6 operator coefficient \f$(C_{eB})_{11}\f$ (real part).
    double CeB_12r; ///< The dimension-6 operator coefficient \f$(C_{eB})_{12}\f$ (real part).
    double CeB_13r; ///< The dimension-6 operator coefficient \f$(C_{eB})_{13}\f$ (real part).
    double CeB_22r; ///< The dimension-6 operator coefficient \f$(C_{eB})_{22}\f$ (real part).
    double CeB_23r; ///< The dimension-6 operator coefficient \f$(C_{eB})_{23}\f$ (real part).
    double CeB_33r; ///< The dimension-6 operator coefficient \f$(C_{eB})_{33}\f$ (real part).
    double CeB_11i; ///< The dimension-6 operator coefficient \f$(C_{eB})_{11}\f$ (imaginary part).
    double CeB_12i; ///< The dimension-6 operator coefficient \f$(C_{eB})_{12}\f$ (imaginary part).
    double CeB_13i; ///< The dimension-6 operator coefficient \f$(C_{eB})_{13}\f$ (imaginary part).
    double CeB_22i; ///< The dimension-6 operator coefficient \f$(C_{eB})_{22}\f$ (imaginary part).
    double CeB_23i; ///< The dimension-6 operator coefficient \f$(C_{eB})_{23}\f$ (imaginary part).
    double CeB_33i; ///< The dimension-6 operator coefficient \f$(C_{eB})_{33}\f$ (imaginary part).
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
    double Ceu_2233;
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
    double CLu_2233;
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
    double CLedQ_11, CLedQ_22, CpLedQ_11, CpLedQ_22;
    double CQQ1_1133, CQQ1_1331, CQQ3_1133, CQQ3_1331;
    double Cuu_1133, Cuu_1331, Cud1_3311, Cud8_3311;
    double CQu1_1133, CQu1_3311, CQu8_1133, CQu8_3311;
    double CQd1_3311, CQd8_3311; 
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
    double eepWBFint; ///< Intrinsic relative theoretical error in \f$e^- p \to H e^- j\f$ via WBF. (Assumed to be constant in energy.)
    double eepWBFpar; ///< Parametric relative theoretical error in \f$e^- p \to H e^- j\f$ via WBF. (Assumed to be constant in energy.)    
    double eepZBFint; ///< Intrinsic relative theoretical error in \f$e^- p \to H e^- j\f$ via ZBF. (Assumed to be constant in energy.)
    double eepZBFpar; ///< Parametric relative theoretical error in \f$e^- p \to H e^- j\f$ via ZBF. (Assumed to be constant in energy.)
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
    
//  Intrinsic errors in \f$ee \to WW\f$ observables: error in total and differential cross section
    double eeeWWint,edeeWWdcint; ///< Intrinsic relative theoretical error in \f$e^+ e^- \to W^+ W^-\f$: total cross section and \f$d/d\cos{\theta}\f$ distribution.
    
    double eggFHgaga,eggFHZga,eggFHZZ,eggFHWW,eggFHtautau,eggFHbb,eggFHmumu; ///< Total relative theoretical error in \f$gg \to H \to X\f$.   
    double eVBFHgaga,eVBFHZga,eVBFHZZ,eVBFHWW,eVBFHtautau,eVBFHbb,eVBFHmumu; ///< Total relative theoretical error in \f$pp \to Hjj (VBF) \to X jj\f$.   
    double eWHgaga,eWHZga,eWHZZ,eWHWW,eWHtautau,eWHbb,eWHmumu; ///< Total relative theoretical error in \f$pp \to WH \to W X\f$.    
    double eZHgaga,eZHZga,eZHZZ,eZHWW,eZHtautau,eZHbb,eZHmumu; ///< Total relative theoretical error in \f$pp \to ZH \to Z X\f$.
    double ettHgaga,ettHZga,ettHZZ,ettHWW,ettHtautau,ettHbb,ettHmumu; ///< Total relative theoretical error in \f$pp \to ttH \to tt X\f$.
    double eVBFHinv,eVHinv; ///< Total relative theoretical error in \f$pp \to X H \to X + invisible\f$.
    
    double eHwidth; ///< Total relative theoretical error in the Higgs width
    
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
        
    double BrHinv; ///< The branching ratio of invisible Higgs decays.
    double BrHexo; ///< The branching ratio of exotic (not invisible) Higgs decays.
    
    double dg1Z; ///< Independent contribution to aTGC.
    double dKappaga; ///< Independent contribution to aTGC.    
    double lambZ; ///< Independent contribution to aTGC.

//  ----------------------------------------------------------------------------
//  Internal parameters    
//  ----------------------------------------------------------------------------
    
    double LambdaNP2;///< The square of the new physics scale [GeV\f$^2\f$].
    
    // Internal values for some of the dimension-6 coefficients (to allow changes of bases in Post-Update)
    // Define also the corresponding anomalous dimension, gADX, for each parameter
    double CiHL1_11;
    double CiHL1_22;
    double CiHL1_33;
    double CiHL3_11;
    double CiHL3_22;
    double CiHL3_33;
    
    double gADHL1_11;
    double gADHL1_22;
    double gADHL1_33;
    double gADHL3_11;
    double gADHL3_22;
    double gADHL3_33;
    
    double CiHQ1_11;
    double CiHQ1_22;
    double CiHQ1_33;
    double CiHQ3_11;
    double CiHQ3_22;
    double CiHQ3_33;
    
    double gADHQ1_11;
    double gADHQ1_22;
    double gADHQ1_33;
    double gADHQ3_11;
    double gADHQ3_22;
    double gADHQ3_33;
    
    double CiHe_11;
    double CiHe_22;
    double CiHe_33;
    
    double gADHe_11;
    double gADHe_22;
    double gADHe_33;
    
    double CiHu_11;
    double CiHu_22;
    double CiHu_33;
    
    double gADHu_11;
    double gADHu_22;
    double gADHu_33;
    
    double CiHd_11;
    double CiHd_22;
    double CiHd_33;
    
    double gADHd_11;
    double gADHd_22;
    double gADHd_33;
    
    double CiW;
    double CiG;   
    
    double gADW;
    double gADG;

    double CiHG;   
    double CiHW;
    double CiHB; 
    double CiHWB;
    double CiDHB; 
    double CiDHW;
    
    double gADHG;   
    double gADHW;
    double gADHB; 
    double gADHWB;
    double gADDHB; 
    double gADDHW;
    
    double CiHbox;
    double CiHD;
    double CiH;
    
    double gADHbox;
    double gADHD;
    double gADH;
    
    double CieH_11r;
    double CieH_22r;
    double CieH_33r;
    
    double gADeH_11r;
    double gADeH_22r;
    double gADeH_33r;
    
    double CiuH_11r;
    double CiuH_22r;
    double CiuH_33r;
    
    double gADuH_11r;
    double gADuH_22r;
    double gADuH_33r;
    
    double CidH_11r;
    double CidH_22r;
    double CidH_33r;
    
    double gADdH_11r;
    double gADdH_22r;
    double gADdH_33r;
    
    double CiuG_11r;
    double CiuG_22r;
    double CiuG_33r;
    
    double gADuG_11r;
    double gADuG_22r;
    double gADuG_33r;
    
    double CiuW_11r;
    double CiuW_22r;
    double CiuW_33r;
    
    double gADuW_11r;
    double gADuW_22r;
    double gADuW_33r;
    
    double CiuB_11r;
    double CiuB_22r;
    double CiuB_33r;
    
    double gADuB_11r;
    double gADuB_22r;
    double gADuB_33r;
    
    double CiLL_1221;
    double CiLL_2112;   
    
    double gADLL_1221;
    
    double v2;///< The square of the EW vev.
    double v2_over_LambdaNP2;///< The ratio between the EW vev and the new physics scale, squared \f$v^2/\Lambda^2\f$.
    double aleMz;///< The em constant at Mz. 
    double eeMz;///< The em coupling at Mz.
    double eeMz2;///< The em coupling squared (at Mz).
    double cW_tree;///< The tree level values for the cosine of the weak angle.
    double sW_tree;///< The tree level values for the sine of the weak angle.
    double cW2_tree;///< The square of the tree level values for the cosine of the weak angle.
    double sW2_tree;///< The square of the tree level values for the sine of the weak angle.
    
    double g1_tree;///< The tree level value of the \f$U(1)_Y\f$ gauge coupling contant (at the \f$Z\f$ pole).
    double g2_tree;///< The tree level value of the \f$SU(2)_L\f$ gauge coupling contant (at the \f$Z\f$ pole).
    double g3_tree;///< The tree level value of the \f$SU(3)_c\f$ gauge coupling contant (at the \f$Z\f$ pole).
    
    double gZvL;///< The tree level value of the \f$Z\bar{\nu}\nu\f$ couplings in the SM.
    double gZlL,gZlR;///< The tree level value of the \f$Z\ell^+\ell^-\f$ couplings in the SM.
    double gZuL,gZuR;///< The tree level value of the \f$Z\bar{u}u\f$ couplings in the SM.
    double gZdL,gZdR;///< The tree level value of the \f$Z\bar{d}d\f$ couplings in the SM.
    
    double UevL;///< The tree level value of the \f$W^-\bar{\ell}\nu\f$ couplings in the SM. (Neglecting PMNS effects.)
    double VudL;///< The tree level value of the \f$W^+\bar{u}d\f$ couplings in the SM. (Neglecting CKM effects.)
        
    double delta_ZZ;///< Combination of dimension 6 coefficients modifying the \f$Z_\mu\f$ canonical field definition.
    double delta_AA;///< Combination of dimension 6 coefficients modifying the \f$A_\mu\f$ canonical field definition.
    double delta_AZ;///< Combination of dimension 6 coefficients modifying the \f$A_\mu\f$ canonical field definition.
    double delta_h;///< Combinations of dimension 6 coefficients modifying the \f$H\f$ canonical field definition.
    
    double lambdaH_tree;///< The SM tree level value of the scalar quartic coupling in the potential.
    
    double dZH;///< Higgs self-coupling contribution to the universal resummed Higgs wave function renormalization
    
    double cHSM;///< Parameter to control the inclusion of modifications of SM parameters in selected Higgs processes.
    
    double cLHd6;///< Parameter to control the inclusion of modifications of SM loops in Higgs processes due to dim 6 interactions.
    
    double cLH3d62;///< Parameter to control the inclusion of modifications of SM loops in Higgs processes due to dim 6 interactions modifying the Higgs trilinear coupling (Quadratic terms).

    double cRGE;///< Parameter to control the inclusion of log-enhanced contributions via RG effects. If activated then it takes the value multiplying the anomalous dimension: \f$-\log(\Lambda/\mu)/16 \pi^2 \Lambda^2\f$
    
    double Yuke,Yukmu,Yuktau;///< SM lepton Yukawas
    double Yuku,Yukc,Yukt;///< SM u-quark Yukawas
    double Yukd,Yuks,Yukb;///< SM d-quark Yukawas
    
    // STXS ci: denote them as ai to separate from the normal conventions 
    double aiG,ai3G,ai2G;
    double aiT,aiH,aiWW,aiB,aiHW,aiHB;
    double aiA;
    double aiHQ,aipHQ,aiHL,aipHL,aiHu,aiHd,aiHe;
    double aiu,aiuG;
    
    // Cache variables
    double dGammaHTotR1,dGammaHTotR2,GammaHTotR;///< NP contributions and Total to Higgs width ratio with SM
    
    double delta_GF;///< The dimension 6 correction to the Fermi constant, as extracted from muon decay
    
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
    bool FlagPartialQFU; ///< A boolean flag that is true if assuming partial quark flavour universality between the 1st and 2nd family in the CHF operators.
    bool FlagFlavU3OfX; ///< A boolean flag that is true if assuming U(3)^5 symmetry in the CfH and CfV operator coefficients.
    bool FlagUnivOfX; ///< A boolean flag that is true if assuming U(3)^5 symmetry in the CfH and CfV operator coefficients and all proportional to the same coefficient (CuH_33 and CuV_33 respectively).
    bool FlagHiggsSM; ///< A boolean flag that is true if including dependence on small variations of the SM parameters (dependence is linearized). Available only in selected Higgs observables. 
    bool FlagLoopHd6; ///< A boolean flag that is true if including modifications in the SM loops in Higgs observables due to the dim 6 interactions.
    bool FlagLoopH3d6Quad; ///< A boolean flag that is true if including quadratic modifications in the SM loops in Higgs observables due to the dim 6 interactions that contribute to the trilinear Higgs coupling.
    bool FlagRGEciLLA; ///< A flag that is TRUE if including log-enhanced 1-loop corrections propotional to the dim-6 Wilson coefficients. See comment in documentation above.
    
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
    
    gsl_integration_cquad_workspace * w_WW;/**< Gsl integral variable */

};

#endif	/* NPSMEFTD6_H */
