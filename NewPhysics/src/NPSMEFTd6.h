/*
 * Copyright (C) 2014 HEPfit Collaboration
 *
 *
 * For the licensing terms see doc/COPYING.
 */

#ifndef NPSMEFTD6_H
#define NPSMEFTD6_H

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
 *   <td class="mod_name">%CQQ1_1133, CQQ1_1331, CQQ1_2233, CQQ1_2332, CQQ1_3333 </td>
 *   <td class="mod_symb">\f$C_{QQ}^{(1)}\f$</td>
 *   <td class="mod_desc">The coefficient of the operator
 *   \f$({\cal O}_{QQ}^{(1)})_{ijkl}=\big(\overline{Q^i}\,\gamma^\mu Q^j\big)
 *   \big(\overline{Q^k}\,\gamma_\mu Q^l\big)\f$. </td>
 * </tr>
 * <tr>
 *   <td class="mod_name">%CQQ3_1133, CQQ3_1331, CQQ3_2233, CQQ3_2332, CQQ3_3333 </td>
 *   <td class="mod_symb">\f$C_{QQ}^{(3)}\f$</td>
 *   <td class="mod_desc">The coefficient of the operator
 *   \f$({\cal O}_{QQ}^{(3)})_{ijkl}=\big(\overline{Q^i}\,\gamma^\mu \sigma_a Q^j\big)
 *   \big(\overline{Q^k}\,\gamma_\mu \sigma_a Q^l\big)\f$. </td>
 * </tr>
 * <tr>
 *   <td class="mod_name">%Cuu_1133, Cuu_1331, Cuu_2233, Cuu_2332, Cuu_3333 </td>
 *   <td class="mod_symb">\f$C_{uu}\f$</td>
 *   <td class="mod_desc">The coefficient of the operator
 *   \f$({\cal O}_{uu})_{ijkl}=\big(\overline{u^i}\,\gamma^\mu u^j\big)
 *   \big(\overline{u^k}\,\gamma_\mu u^l\big)\f$. </td>
 * </tr>
 * <tr>
 *   <td class="mod_name">%Cud1_3311, Cud1_3322, Cud1_3333 </td>
 *   <td class="mod_symb">\f$C_{ud}^{(1)}\f$</td>
 *   <td class="mod_desc">The coefficient of the operator
 *   \f$({\cal O}_{ud}^{(1)})_{ijkl}=\big(\overline{u^i}\,\gamma^\mu u^j\big)
 *   \big(\overline{d^k}\,\gamma_\mu d^l\big)\f$. </td>
 * </tr>
 * <tr>
 *   <td class="mod_name">%Cud8_3311, Cud8_3322, Cud8_3333 </td>
 *   <td class="mod_symb">\f$C_{ud}^{(8)}\f$</td>
 *   <td class="mod_desc">The coefficient of the operator
 *   \f$({\cal O}_{ud}^{(8)})_{ijkl}=\big(\overline{u^i}\,\gamma^\mu T^A u^j\big)
 *   \big(\overline{d^k}\,\gamma_\mu T^A d^l\big)\f$. </td>
 * </tr>
 * <tr>
 *   <td class="mod_name">%CQu1_1133, CQu1_3311, CQu1_2233, CQu1_3322, CQu1_3333 </td>
 *   <td class="mod_symb">\f$C_{Qu}^{(1)}\f$</td>
 *   <td class="mod_desc">The coefficient of the operator
 *   \f$({\cal O}_{Qu}^{(1)})_{ijkl}=\big(\overline{Q^i}\,\gamma^\mu Q^j\big)
 *   \big(\overline{u^k}\,\gamma_\mu u^l\big)\f$. </td>
 * </tr>
 * <tr>
 *   <td class="mod_name">%CQu8_1133, CQu8_3311, CQu8_2233, CQu8_3322, CQu8_3333 </td>
 *   <td class="mod_symb">\f$C_{Qu}^{(8)}\f$</td>
 *   <td class="mod_desc">The coefficient of the operator
 *   \f$({\cal O}_{Qu}^{(8)})_{ijkl}=\big(\overline{Q^i}\,\gamma^\mu T^A Q^j\big)
 *   \big(\overline{u^k}\,\gamma_\mu T^A u^l\big)\f$. </td>
 * </tr>
 * <tr>
 *   <td class="mod_name">%CQd1_3311, CQd1_3322, CQd1_3333 </td>
 *   <td class="mod_symb">\f$C_{Qd}^{(1)}\f$</td>
 *   <td class="mod_desc">The coefficient of the operator
 *   \f$({\cal O}_{Qd}^{(1)})_{ijkl}=\big(\overline{Q^i}\,\gamma^\mu Q^j\big)
 *   \big(\overline{d^k}\,\gamma_\mu d^l\big)\f$. </td>
 * </tr>
 * <tr>
 *   <td class="mod_name">%CQd8_3311, CQd8_3322, CQd8_3333 </td>
 *   <td class="mod_symb">\f$C_{Qd}^{(8)}\f$</td>
 *   <td class="mod_desc">The coefficient of the operator
 *   \f$({\cal O}_{Qd}^{(8)})_{ijkl}=\big(\overline{Q^i}\,\gamma^\mu T^A Q^j\big)
 *   \big(\overline{d^k}\,\gamma_\mu T^A d^l\big)\f$. </td>
 * </tr> 
 * <tr>
 *   <td class="mod_name">%CQuQd1_3333 </td>
 *   <td class="mod_symb">\f$C_{QuQd}^{(1)}\f$</td>
 *   <td class="mod_desc">The coefficient of the operator
 *   \f$({\cal O}_{QuQd}^{(1)})_{ijkl}=\big(\overline{Q^i}\, u^j\big) i\sigma_2
 *   \big(\overline{Q^k}^T\, d^l\big)\f$. </td>
 * </tr>
 * <tr>
 *   <td class="mod_name">%CQuQd8_3333 </td>
 *   <td class="mod_symb">\f$C_{QuQd}^{(8)}\f$</td>
 *   <td class="mod_desc">The coefficient of the operator
 *   \f$({\cal O}_{QuQd}^{(8)})_{ijkl}=\big(\overline{Q^i}\, T^A u^j\big) i\sigma_2
 *   \big(\overline{Q^k}^T\, T^A d^l\big)\f$. </td>
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
 *   <td class="mod_name">%CQuQd1 </td>
 *   <td class="mod_symb">\f$C_{QuQd}^{(1)}\f$</td>
 *   <td class="mod_desc">The coefficient of the operator
 *   \f$({\cal O}_{QuQd}^{(1)})_{ijkl}=\big(\overline{Q^i}\, u^j\big) i\sigma_2
 *   \big(\overline{Q^k}^T\, d^l\big)\f$. </td>
 * </tr>
 * <tr>
 *   <td class="mod_name">%CQuQd8 </td>
 *   <td class="mod_symb">\f$C_{QuQd}^{(8)}\f$</td>
 *   <td class="mod_desc">The coefficient of the operator
 *   \f$({\cal O}_{QuQd}^{(8)})_{ijkl}=\big(\overline{Q^i}\, T^A u^j\big) i\sigma_2
 *   \big(\overline{Q^k}^T\, T^A d^l\big)\f$. </td>
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
 * <tr>
 *   <td class="mod_name">%MWinput</td>
 *   <td class="mod_valu">TRUE&nbsp;/&nbsp;<b>FALSE</b></td>
 *   <td class="mod_desc">This auxiliary flag is used for setting the W mass as a SM input, 
 *   instead of the electromagnetic constant parameter %dAle5Mz. 
 *   The default value is FALSE. This flag must be used together with the flag of the same name in StandardModel</td>
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
    static const int NNPSMEFTd6Vars = 557;

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
    static const int NNPSMEFTd6Vars_LFU_QFU = 292; //AG:modified,hat

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
    virtual bool CheckParameters(const std::map<std::string, double >& DPars);

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
    virtual NPSMEFTd6Matching& getMatching() const {
        return NPSMEFTd6M.getObj();
    }

    ////////////////////////////////////////////////////////////////////////

    /**
     * @brief Return Lambda_NP
     * @return \f$ Lambda_NP \f$ 
     */
    const double getLambda_NP() const {
        return Lambda_NP;
    }

    /**
     * @brief Return NP coeff CQe_2322
     * @return \f$ CQe_2322 \f$ 
     */
    const double getCQe_2311() const {
        return CQe_2311;
    }

    /**
     * @brief Return NP coeff CQe_2322
     * @return \f$ CQe_2311 \f$ 
     */
    const double getCQe_2322() const {
        return CQe_2322;
    }

    /**
     * @brief Return NP coeff CLQ1_1123
     * @return \f$ CLQ1_1123 \f$ 
     */
    const double getCLQ1_1123() const {
        return CLQ1_1123;
    }

    /**
     * @brief Return NP coeff CLQ1_2223
     * @return \f$ CLQ1_2223 \f$ 
     */
    const double getCLQ1_2223() const {
        return CLQ1_2223;
    }

    /**
     * @brief Return NP coeff CLQ1_3323
     * @return \f$ CLQ1_3323 \f$ 
     */
    const double getCLQ1_3323() const {
        return CLQ1_3323;
    }

    /**
     * @brief Return NP coeff CLQ3_1123
     * @return \f$ CLQ3_1123 \f$ 
     */
    const double getCLQ3_1123() const {
        return CLQ3_1123;
    }

    /**
     * @brief Return NP coeff CLQ3_2223
     * @return \f$ CLQ3_2223 \f$ 
     */
    const double getCLQ3_2223() const {
        return CLQ3_2223;
    }

    /**
     * @brief Return NP coeff CLQ3_3323
     * @return \f$ CLQ3_3323 \f$ 
     */
    const double getCLQ3_3323() const {
        return CLQ3_3323;
    }

    /**
     * @brief Return NP coeff CLd_1123
     * @return \f$ CLd_1123 \f$ 
     */
    const double getCLd_1123() const {
        return CLd_1123;
    }

    /**
     * @brief Return NP coeff CLd_2223
     * @return \f$ CLd_2223 \f$ 
     */
    const double getCLd_2223() const {
        return CLd_2223;
    }

    /**
     * @brief Return NP coeff CLd_3323
     * @return \f$ CLd_3323 \f$ 
     */
    const double getCLd_3323() const {
        return CLd_3323;
    }

    /**
     * @brief Return NP coeff Ced_1123
     * @return \f$ Ced_1123 \f$ 
     */
    const double getCed_1123() const {
        return Ced_1123;
    }

    /**
     * @brief Return NP coeff Ced_2223
     * @return \f$ Ced_2223 \f$ 
     */
    const double getCed_2223() const {
        return Ced_2223;
    }

    /**
     * @brief Return NP coeff CHL1_11
     * @return \f$ CHL1_11 \f$ 
     */
    const double getCHL1_11() const {
        return CHL1_11;
    }

    /**
     * @brief Return NP coeff CHL1_22
     * @return \f$ CHL1_22 \f$ 
     */
    const double getCHL1_22() const {
        return CHL1_22;
    }

    /**
     * @brief Return NP coeff CHL3_11
     * @return \f$ CHL3_11 \f$ 
     */
    const double getCHL3_11() const {
        return CHL3_11;
    }

    /**
     * @brief Return NP coeff CHL3_22
     * @return \f$ CHL3_22 \f$ 
     */
    const double getCHL3_22() const {
        return CHL3_22;
    }

    /**
     * @brief Return NP coeff CHe_11
     * @return \f$ CHe_11 \f$ 
     */
    const double getCHe_11() const {
        return CHe_11;
    }

    /**
     * @brief Return NP coeff CHe_22
     * @return \f$ CHe_22 \f$ 
     */
    const double getCHe_22() const {
        return CHe_22;
    }

    /**
     * @brief Return NP coeff Ceu_1133
     * @return \f$ Ceu_1133 \f$ 
     */
    const double getCeu_1133() const {
        return Ceu_1133;
    }

    /**
     * @brief Return NP coeff Ceu_2233
     * @return \f$ Ceu_2233 \f$ 
     */
    const double getCeu_2233() const {
        return Ceu_2233;
    }

    /**
     * @brief Return NP coeff CLu_1133
     * @return \f$ CLu_1133 \f$ 
     */
    const double getCLu_1133() const {
        return CLu_1133;
    }

    /**
     * @brief Return NP coeff CLu_2233
     * @return \f$ CLu_2233 \f$ 
     */
    const double getCLu_2233() const {
        return CLu_2233;
    }

    /**
     * @brief Return NP coeff CLedq_11
     * @return \f$ CLedq_11 \f$ 
     */
    const double getCLedQ_11() const {
        return CLedQ_11;
    }

    /**
     * @brief Return NP coeff CLedq_22
     * @return \f$ CLedq_22 \f$ 
     */
    const double getCLedQ_22() const {
        return CLedQ_22;
    }

    /**
     * @brief Return NP coeff CpLedq_11
     * @return \f$ CpLedq_11 \f$ 
     */
    const double getCpLedQ_11() const {
        return CpLedQ_11;
    }

    /**
     * @brief Return NP coeff CpLedq_22
     * @return \f$ CpLedq \f$ 
     */
    const double getCpLedQ_22() const {
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

    //  Functions used to compute the 1-loop dependence of single Higgs observables 
    //  on the Higgs self-coupling
        
    /**
     * @brief The coefficient of the 1-loop linear term in the Higgs selfcoupling
     */
    virtual const double deltaH3L1(double C1) const;
    
    /**
     * @brief The coefficient of the 1-loop quadratic term in the Higgs selfcoupling
     */
    virtual const double deltaH3L2(double C1) const;
    
    /////////////// QUADRATIC CORRECTIONS ON SIGMA AND SIGMA X BR FROM H3 ////////
    
    /**
     * @brief Quadratic contribution from the Higgs self-couplings modifications to the signal strength for an observable @f$\sigma@f$ in the current model.
     * @return @f$\delta^{(2)}_{h^3}@f$
     */
    virtual const double delta2sH3(const double C1) const;
    
    /**
     * @brief Quadratic contribution from the Higgs self-couplings modifications to the signal strength for @f$\sigma \times BR(H\to xx)@f$ in the current model.
     * @return @f$\delta^{(2)}_{h^3}@f$
     */
    virtual const double delta2sBRH3(const double C1prod, const double C1Hxx) const;
        
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
    virtual const double DeltaGF() const;

    /**
     * @brief The oblique parameter \f$S\f$.
     * (Simplified implementation. Contribution only from @f$O_{HWB}@f$.)
     * @return the value of @f$S@f$
     */
    virtual const double obliqueS() const;

    /**
     * @brief The oblique parameter \f$T\f$.
     * (Simplified implementation. Contribution only from @f$O_{HD}@f$.)
     * @return the value of @f$T@f$
     */
    virtual const double obliqueT() const;

    /**
     * @brief The oblique parameter \f$U\f$.
     * @return the value of @f$U@f$
     */
    virtual const double obliqueU() const;

    /**
     * @brief The oblique parameter \f$W\f$.
     * (Simplified implementation. Contribution only from @f$O_{2W}@f$.)
     * @return the value of @f$W@f$
     */
    virtual const double obliqueW() const;

    /**
     * @brief The oblique parameter \f$Y\f$.
     * (Simplified implementation. Contribution only from @f$O_{2B}@f$.)
     * @return the value of @f$Y@f$
     */
    virtual const double obliqueY() const;


    /////////////////////////////// Deviations in the experimental values of the SM input parameters /////////////////////////////////////////

    /**
     * @brief The relative correction to the mass of the @f$Z@f$ boson, @f$\delta M_Z/M_Z@f$, with respect to ref. point used in the SM calculation of Higgs observables.
     * @return @f$\delta M_Z/M_Z@f$
     */
    virtual const double deltaMz() const;

    /**
     * @brief The relative correction to the mass of the @f$Z@f$ boson squared, @f$(\delta M_Z/M_Z)^2@f$, with respect to ref. point used in the SM calculation of Higgs observables.
     * @return @f$(\delta M_Z/M_Z)^2@f$
     */
    virtual const double deltaMz2() const;

    /**
     * @brief The relative correction to the mass of the @f$H@f$ boson, @f$\delta M_H/M_H@f$, with respect to ref. point used in the SM calculation of Higgs observables.
     * @return @f$\delta M_H/M_H@f$
     */
    virtual const double deltaMh() const;

    /**
     * @brief The relative correction to the mass of the @f$H@f$ boson squared, @f$(\delta M_H/M_H)^2@f$, with respect to ref. point used in the SM calculation of Higgs observables.
     * @return @f$(\delta M_H/M_H)^2@f$
     */
    virtual const double deltaMh2() const;

    /**
     * @brief The relative correction to the mass of the @f$t@f$ quark, @f$\delta m_t/m_t@f$, with respect to ref. point used in the SM calculation of Higgs observables.
     * @return @f$\delta m_t/m_t@f$
     */
    virtual const double deltamt() const;

    /**
     * @brief The relative correction to the mass of the @f$t@f$ quark squared, @f$(\delta m_t/m_t)^2@f$, with respect to ref. point used in the SM calculation of Higgs observables.
     * @return @f$(\delta m_t/m_t)^2@f$
     */
    virtual const double deltamt2() const;

    /**
     * @brief The relative correction to the mass of the @f$b@f$ quark, @f$\delta m_b/m_b@f$, with respect to ref. point used in the SM calculation of Higgs observables.
     * @return @f$\delta m_b/m_b@f$
     */
    virtual const double deltamb() const;

    /**
     * @brief The relative correction to the mass of the @f$b@f$ quark squared, @f$(\delta m_b/m_b)^2@f$, with respect to ref. point used in the SM calculation of Higgs observables.
     * @return @f$(\delta m_b/m_b)^2@f$
     */
    virtual const double deltamb2() const;

    /**
     * @brief The relative correction to the mass of the @f$c@f$ quark, @f$\delta m_c/m_c@f$, with respect to ref. point used in the SM calculation of Higgs observables.
     * @return @f$\delta m_c/m_c@f$
     */
    virtual const double deltamc() const;

    /**
     * @brief The relative correction to the mass of the @f$c@f$ quark squared, @f$(\delta m_c/m_c)^2@f$, with respect to ref. point used in the SM calculation of Higgs observables.
     * @return @f$(\delta m_c/m_c)^2@f$
     */
    virtual const double deltamc2() const;

    /**
     * @brief The relative correction to the mass of the @f$\tau@f$ lepton, @f$\delta m_\tau/m_\tau@f$, with respect to ref. point used in the SM calculation of Higgs observables.
     * @return @f$\delta m_\tau/m_\tau@f$
     */
    virtual const double deltamtau() const;

    /**
     * @brief The relative correction to the mass of the @f$\tau@f$ lepton squared, @f$(\delta m_\tau/m_\tau)^2@f$, with respect to ref. point used in the SM calculation of Higgs observables.
     * @return @f$(\delta m_\tau/m_\tau)^2@f$
     */
    virtual const double deltamtau2() const;

    /**
     * @brief The relative correction to the muon decay constant, @f$\delta G_\mu/G_\mu@f$, with respect to ref. point used in the SM calculation of Higgs observables.
     * @return @f$\delta G_\mu/G_\mu@f$
     */
    virtual const double deltaGmu() const;

    /**
     * @brief The relative correction to the muon decay constant, @f$(\delta G_\mu/G_\mu)^2@f$, with respect to ref. point used in the SM calculation of Higgs observables.
     * @return @f$(\delta G_\mu/G_\mu)^2@f$
     */
    virtual const double deltaGmu2() const;

    /**
     * @brief The relative correction to the electromagnetic constant at the Z pole, @f$\delta \alpha(M_Z^2)/\alpha(M_Z^2)@f$, with respect to ref. point used in the SM calculation of Higgs observables.
     * @return @f$\delta \alpha(M_Z^2)/\alpha(M_Z^2)@f$
     */
    virtual const double deltaaMZ() const;

    /**
     * @brief The relative correction to the electromagnetic constant at the Z pole, @f$(\delta \alpha(M_Z^2)/\alpha(M_Z^2))^2@f$, with respect to ref. point used in the SM calculation of Higgs observables.
     * @return @f$(\delta \alpha(M_Z^2)/\alpha(M_Z^2))^2@f$
     */
    virtual const double deltaaMZ2() const;

    /**
     * @brief The relative correction to the electromagnetic constant at zero momentum, @f$\delta \alpha(0)/\alpha(0)@f$, with respect to ref. point used in the SM calculation of Higgs observables.
     * @return @f$\delta \alpha(0)/\alpha(0)@f$
     */
    virtual const double deltaa0() const;

    /**
     * @brief The relative correction to the electromagnetic constant at zero momentum, @f$(\delta \alpha(0)/\alpha(0))^2@f$, with respect to ref. point used in the SM calculation of Higgs observables.
     * @return @f$(\delta \alpha(0)/\alpha(0))^2@f$
     */
    virtual const double deltaa02() const;

    /**
     * @brief The relative correction to the strong coupling constant at the Z pole, @f$\delta \alpha_S(M_Z^2)/\alpha_S(M_Z^2)@f$, with respect to ref. point used in the SM calculation of Higgs observables.
     * @return @f$\delta \alpha_S(M_Z^2)/\alpha_S(M_Z^2)@f$
     */
    virtual const double deltaaSMZ() const;

    /**
     * @brief The relative correction to the strong coupling constant at the Z pole, @f$(\delta \alpha_S(M_Z^2)/\alpha_S(M_Z^2))^2@f$, with respect to ref. point used in the SM calculation of Higgs observables.
     * @return @f$(\delta \alpha_S(M_Z^2)/\alpha_S(M_Z^2))^2@f$
     */
    virtual const double deltaaSMZ2() const;

    /**
     * @brief The relative correction to the mass of the @f$W@f$ boson, @f$\delta M_W/M_W@f$, with respect to ref. point used in the SM calculation of Higgs observables.
     * @return @f$\delta M_W/M_W@f$
     */
    virtual const double deltaMw() const;

    /**
     * @brief The relative correction to the mass of the @f$W@f$ boson squared, @f$(\delta M_W/M_W)^2@f$, with respect to ref. point used in the SM calculation of Higgs observables.
     * @return @f$(\delta M_W/M_W)^2@f$
     */
    virtual const double deltaMw2() const;

    ////////////////////////////////////////////////////////////////////////

    /**
     * @brief The electromagnetic coupling at the @f$Z@f$-mass scale.
     * @details  
     * \f$\newcommand{\nc}{\newcommand}\f$
     * \f$\nc{\gwt}{\widetilde{g}_W}\f$
     * \f$\nc{\gbt}{\widetilde{g}_1}\f$
     * \f$\nc{\gbL}{\Delta g_1^{(1)}}\f$
     * \f$\nc{\gwL}{\Delta g_W^{(1)}}\f$
     * \f$\nc{\gbQ}{\Delta g_1^{(2)}}\f$
     * \f$\nc{\gwQ}{\Delta g_W^{(2)}}\f$
     * \f$\nc{\vt}{\widetilde{v}}\f$
     * \f$\nc{\dGf}{\delta_{G_F}}\f$
     * \f$\nc{\swt}{s_W}\f$
     * \f$\nc{\cwt}{c_W}\f$
     * \f[
     * \alpha (M_Z) = \alpha (M_Z)^{SM} + \Delta \alpha (Mz)^{(1)} + \Delta \alpha (Mz)^{(2)}
     * \f]
     * where,
     * \f[
     * \Delta \alpha^{(1)} = 2\swt^2\frac{\gwL}{\gwt} + 2\cwt^2\frac{\gbL}{\gbt} - 4\swt\cwt C_{\varphi WB}\frac{\vt^2}{\Lambda^2}
     * \f]
     * \f{eqnarray*}{
     * \Delta \alpha^{(2)} &=& 2\swt^2\frac{\gwQ}{\gwt} + 2\cwt^2\frac{\gbQ}{\gbt}
     * + \swt^2(\swt^2-3\cwt^2)\left(\frac{\gwL}{\gwt}\right)^2 + \cwt^2(\cwt^2-3\swt^2)\left(\frac{\gbL}{\gbt}\right)^2 
     * + 8 \swt^2\cwt^2 \frac{\gwL}{\gwt}\frac{\gbL}{\gbt} \\
     * &+& 2 \swt\cwt\left((\cwt^2-3\swt^2)\frac{\gwL}{\gwt} 
     * + (\swt^2-3\cwt^2)\frac{\gbL}{\gbt}\right) \hat{C}_{\varphi WB} \vt^2 
     * - 2 \swt \cwt \left( -2\swt\cwt \hat{C}_{\varphi WB}\vt^2 
     * + (\hat{C}_{\varphi W}+\hat{C}_{\varphi B})\vt^2 + \dGf^{(1)} \right) \hat{C}_{\varphi WB}\vt^2
     * \f}
     * @return @f$\alpha(M_Z^2)@f$, including @f$\mathcal{O}(\Lambda^{-2})@f$ and @f$\mathcal{O}(\Lambda^{-4})@f$ NP contributions.
     */
    virtual const double alphaMz() const; //AG:modified

    /**
     * @brief The mass of the @f$W@f$ boson, @f$M_W@f$.
     * @details
     * \f[
     * M_W = M_W^{SM}( 1 + \Delta M_W^{(1)} + \Delta M_W^{(2)})
     * \f]
     * @return @f$M_W@f$ in GeV, including @f$\mathcal{O}(\Lambda^{-2})@f$ and @f$\mathcal{O}(\Lambda^{-4})@f$ NP contributions.
     */
    virtual const double Mw() const; //AG:modified

    /**
     * @brief The relative NP corrections to the mass of the @f$W@f$ boson, @f$\delta M_W/M_W@f$.
     * @return @f$\delta M_W/M_W@f$
     */
    virtual const double deltaMwd6() const;

    /**
     * @brief The relative NP corrections to the mass of the @f$W@f$ boson squared, @f$(\delta M_W/M_W)^2@f$.
     * @return @f$(\delta M_W/M_W)^2@f$
     */
    virtual const double deltaMwd62() const;

    /**
     * @brief The relative \f$\mathcal{O}(\Lambda^{-4})\f$ NP corrections to the mass of the @f$W@f$ boson, @f$\Delta M_W^{(2)}@f$.
     * @details
     * \f$\nc{\dGf}{\delta_{G_F}}\f$
     * \f[
     * \Delta M_W^{(2)} =  \frac{\gwQ}{\gwt} + \frac{\dGf^{(2)}}{2} + \frac{\gwL \dGf^{(1)}}{2\gwt} - \frac{(\dGf^{(1)})^2}{8}
     * \f]
     * @return @f$\Delta M_W^{(2)}@f$
     */
    virtual const double deltaMwd6_2() const; //AG:added

    virtual const double deltaGamma_Wff_2(const Particle fi, const Particle fj) const; //AG:added

    /**
     * @brief The new physics contribution to the decay width of the @f$W@f$ boson into a given fermion pair, @f$\delta \Gamma_Z^{f}@f$.
     * @param[in] fi a lepton or quark
     * @param[in] fj a lepton or quark
     * @return @f$\delta \Gamma_W^{ff}@f$ in GeV
     */
    virtual const double deltaGamma_Wff(const Particle fi, const Particle fj) const;

    /**
     * @brief A partial decay width of the @f$W@f$ boson decay into a SM fermion pair.
     * @param[in] fi a lepton or quark
     * @param[in] fj a lepton or quark
     * @return @f$\Gamma^W_{ij}@f$
     */
    virtual const double GammaW(const Particle fi, const Particle fj) const;

    virtual const double deltaGamma_W_2() const; //AG:added

    /**
     * @brief The new physics contribution to the total decay width of the @f$W@f$ boson, @f$\delta \Gamma_W@f$.
     * @return @f$\delta \Gamma_W@f$ in GeV
     */
    virtual const double deltaGamma_W() const;

    /**
     * @brief The total width of the @f$W@f$ boson, @f$\Gamma_W@f$.
     * @return @f$\Gamma_W@f$ in GeV
     */
    virtual const double GammaW() const;

    /**
     * @brief The relative NP corrections to the width of the @f$W@f$ boson, @f$\delta \Gamma_W/\Gamma_W@f$.
     * @return @f$\delta \Gamma_W/\Gamma_W@f$
     */
    virtual const double deltaGwd6() const;

    /**
     * @brief The relative NP corrections to the width of the @f$W@f$ boson squared, @f$(\delta \Gamma_W/\Gamma_W)^2@f$.
     * @return @f$(\delta \Gamma_W/\Gamma_W)^2@f$
     */
    virtual const double deltaGwd62() const;

    /**
     * @brief The relative NP corrections to the width of the @f$Z@f$ boson, @f$\delta \Gamma_Z/\Gamma_Z@f$.
     * @return @f$\delta \Gamma_Z/\Gamma_Z@f$
     */
    virtual const double deltaGzd6() const;

    /**
     * @brief The relative NP corrections to the width of the @f$Z@f$ boson squared, @f$(\delta \Gamma_Z/\Gamma_Z)^2@f$.
     * @return @f$(\delta \Gamma_Z/\Gamma_Z)^2@f$
     */
    virtual const double deltaGzd62() const;

    /**
     * @brief The \f$\mathcal{O}(\Lambda^{-4})\f$ new physics contribution to the neutral-current vector coupling @f$g_V^f@f$.
     * @details
     * \f[
     * \Delta g_{V,f}^{(2)} = \Delta g_{L,f}^{(2)} + \Delta g_{R,f}^{(2)}
     * \f]
     * @param[in] f a lepton or quark
     * @return @f$\Delta g_{V,f}^{(2)}@f$
     */
    virtual const double deltaGV_f_2(const Particle p) const; //AG:added

    /**
     * @brief New physics contribution to the neutral-current vector coupling @f$g_V^f@f$.
     * @param[in] f a lepton or quark
     * @return @f$\delta g_V^f@f$
     */
    virtual const double deltaGV_f(const Particle p) const;

    /**
     * @brief The \f$\mathcal{O}(\Lambda^{-4})\f$ new physics contribution to the neutral-current vector coupling @f$g_A^f@f$.
     * @details
     * \f[
     * \Delta g_{A,f}^{(2)} = \Delta g_{L,f}^{(2)} - \Delta g_{R,f}^{(2)}
     * \f]
     * @param[in] f a lepton or quark
     * @return @f$\Delta g_{A,f}^{(2)}@f$
     */
    virtual const double deltaGA_f_2(const Particle p) const; //AG:added

    /**
     * @brief New physics contribution to the neutral-current axial-vector coupling @f$g_A^f@f$.
     * @param[in] f a lepton or quark
     * @return @f$\delta g_A^f@f$
     */
    virtual const double deltaGA_f(const Particle p) const;

    /**
     * @brief The \f$\mathcal{O}(\Lambda^{-4})\f$ new physics contribution to the left-handed coupling @f$g_L^f@f$.
     * @details
     * \f$\newcommand{\nc}{\newcommand}\f$
     * \f$\nc{\vt}{\widetilde{v}}\f$
     * \f$\nc{\xWZb}{\mathbb{\bar{X}}_{(WZ)}}\f$
     * \f$\nc{\xBZb}{\mathbb{\bar{X}}_{(BZ)}}\f$
     * \f$\nc{\xWZt}{\mathbb{\widetilde{X}}_{(WZ)}}\f$
     * \f$\nc{\xBZt}{\mathbb{\widetilde{X}}_{(BZ)}}\f$
     * \f$\nc{\xWZL}{\Delta \mathbb{X}_{WZ}^{(1)}}\f$
     * \f$\nc{\xBZL}{\Delta \mathbb{X}_{BZ}^{(1)}}\f$
     * \f$\nc{\xWZQ}{\Delta \mathbb{X}_{WZ}^{(2)}}\f$
     * \f$\nc{\xBZQ}{\Delta \mathbb{X}_{BZ}^{(2)}}\f$
     * \f$\nc{\gwt}{\widetilde{g}_W}\f$
     * \f$\nc{\gbt}{\widetilde{g}_1}\f$
     * \f$\nc{\gbL}{\Delta g_1^{(1)}}\f$
     * \f$\nc{\gwL}{\Delta g_W^{(1)}}\f$
     * \f$\nc{\gbQ}{\Delta g_1^{(2)}}\f$
     * \f$\nc{\gwQ}{\Delta g_W^{(2)}}\f$
     * \f{eqnarray*}{
     * \Delta g_L^{(2)} &=& \frac{1}{\sqrt{\gbt^2+\gwt^2}} \left[ (Q-I_3) ~\left(\gbt ~\xBZQ + \gbL \xBZL + \xBZt \gbQ \right) 
     * + I_3 ~\left(\gwt \xWZQ + \gwL \xWZL + \xWZt \gwQ \right) \right] \\
     * &+& \frac{(\hat{C}_{Hf1} - 2I_3\hat{C}_{Hf3}) \vt^2}{2} ~\frac{\xBZt \gbL + \gbt \xBZL  - \xWZt \gwL - \gwt \xWZL }{\sqrt{\gbt^2+\gwt^2}}
     * \f}
     * @param[in] f a lepton or quark
     * @return @f$\Delta g_{L,f}^{(2)}@f$
     */
    const double deltaGL_f_2(const Particle p) const; //AG:added    

    /**
     * @brief New physics contribution to the neutral-current left-handed coupling @f$g_L^f@f$.
     * @param[in] f a lepton or quark
     * @return @f$\delta g_L^f@f$
     */
    // no generation mixing
    const double deltaGL_f(const Particle p) const;

    /**
     * @brief The \f$\mathcal{O}(\Lambda^{-4})\f$ new physics contribution to the right-handed coupling @f$g_R^f@f$.
     * @details
     * \f[
     * \Delta g_R^{(2)} = \frac{Q}{\sqrt{\gbt^2+\gwt^2}} \left( \gbt \xBZQ + \gbL \xBZL + \xBZt \gbQ \right) 
     * + \frac{\hat{C}_{Hf} \vt^2}{2} ~\frac{\xBZt \gbL + \gbt \xBZL  - \xWZt \gwL - \gwt \xWZL }{\sqrt{\gbt^2+\gwt^2}}
     * \f]
     * @param[in] f a lepton or quark
     * @return @f$\Delta g_{R,f}^{(2)}@f$
     */
    const double deltaGR_f_2(const Particle p) const; //AG:added    

    /**
     * @brief New physics contribution to the neutral-current right-handed coupling @f$g_R^f@f$.
     * @param[in] f a lepton or quark
     * @return @f$\delta g_R^f@f$
     */
    // no generation mixing
    const double deltaGR_f(const Particle p) const;
    
    /**
     * @brief The branching ratio of the @f$W@f$ boson decaying into a SM fermion pair, @f$Br(W\to f_i f_j)@f$.
     * @details
     * \f[ BrW_{fi,fj} = BrW_{fi,fj}^{SM} + \Delta BrW_{fi,fj}^{(1)} + \Delta BrW_{fi,fj}^{(2)} \f]
     * where,
     * \f[
     * \Delta BrW_{fi,fj}^{(1)} = \frac{\Gamma_{W, f_i f_j}^{SM}}{\Gamma_W^{SM}} \left( \frac{\Delta \Gamma_{W, f_i f_j}^{(1)}}{\Gamma_{W, f_i f_j}^{SM}} - \frac{\Delta \Gamma_{W}^{(1)}}{\Gamma_{W}^{SM}} \right)
     * \f]
     * \f[
     * \Delta BrW_{fi,fj}^{(2)} = \frac{\Gamma_{W,f_i f_j}^{SM}}{\Gamma_W^{SM}} \left( \frac{\Delta \Gamma_{W, f_i f_j}^{(2)}}{\Gamma_{W, f_i f_j}^{SM}}
     * - \frac{\Delta \Gamma_{W}^{(2)}}{\Gamma_{W}^{SM}}
     * + \frac{(\Delta \Gamma_{W}^{(1)})^2}{(\Gamma_{W}^{SM})^2} 
     * + \frac{\Delta \Gamma_{W,f_i f_j}^{(1)} \Delta \Gamma_{W}^{(1)}}{\Gamma_{W, f_i f_j}^{SM} \Gamma_{W}^{SM}} \right)
     * \f]
     * @return @f$Br(W\to f_i f_j)@f$ in GeV, including @f$\mathcal{O}(\Lambda^{-2})@f$ and @f$\mathcal{O}(\Lambda^{-4})@f$ NP contributions.
     */
    virtual const double BrW(const Particle fi, const Particle fj) const; //AG:modified

    /**
     * @brief The lepton universality ratio @f$R_{W,l_i/l_j)=\Gamma(W\to l_i \nu_i)/\Gamma(W\to l_j \nu_j)@f$.
     * @details
     * \f[ RWc = RWc^{SM} + \Delta RWc^{(1)} + \Delta RWc^{(2)} \f]
     * where,
     * \f[
     * \Delta RWc^{(2)} = \frac{\Delta \Gamma_{W,cx}^{(2)}}{\Gamma_{W,had}^{SM}} 
     * + \frac{\Gamma_{W,cx}^{SM} ~\Delta \Gamma_{W,had}^{(1)}}{\Gamma_{W,had}^{SM}} 
     * + \frac{\Gamma_{W,cx}^{SM} ~(\Delta \Gamma_{W,had}^{(1)})^2}{(\Gamma_{W,had}^{SM})^3} 
     * - \frac{\Delta \Gamma_{W,cx}^{(1)} ~\Delta \Gamma_{W,had}^{(1)}}{(\Gamma_{W,had}^{SM})^2}
     * \f]
     * @return @f$R_{W,l_i/l_j)@f$ in GeV
     */
    virtual const double RWlilj(const Particle li, const Particle lj) const;

    /**
     * @brief The ratio @f$ R_{W,c}=\Gamma(W\to c + X)/\Gamma(W\to had) @f$.
     * @return \f$R_{W,c}\f$ in GeV, including @f$\mathcal{O}(\Lambda^{-2})@f$ and @f$\mathcal{O}(\Lambda^{-4})@f$ NP contributions.
     */
    virtual const double RWc() const; //AG:modified

    /**
     * @brief The lepton universality ratio @f$R_{Z,l_i/l_j)=\Gamma(Z\to l_i^+ l_i^-)/\Gamma(Z\to l_j^+ l_j^-)@f$.
     * @return @f$R_{Z,l_i/l_j)@f$ in GeV
     */
    virtual const double RZlilj(const Particle li, const Particle lj) const;

    ////////////////////////////////////////////////////////////////////////
    
    // LEP2 definitions for 2 to 2 fermion processes

    const double CeeLL_e() const;
    const double CeeLL_mu() const;
    const double CeeLL_tau() const;
    
    const double CeeLL_up() const;
    const double CeeLL_charm() const;
    const double CeeLL_top() const;
    
    const double CeeLL_down() const;
    const double CeeLL_strange() const;
    const double CeeLL_bottom() const;


    const double CeeLR_e() const;
    const double CeeLR_mu() const;
    const double CeeLR_tau() const;
    
    const double CeeLR_up() const;
    const double CeeLR_charm() const;
    const double CeeLR_top() const;
    
    const double CeeLR_down() const;
    const double CeeLR_strange() const;
    const double CeeLR_bottom() const;

    const double CeeRL_e() const;
    const double CeeRL_mu() const;
    const double CeeRL_tau() const;
    
    const double CeeRL_up() const;
    const double CeeRL_charm() const;
    const double CeeRL_top() const;
    
    const double CeeRL_down() const;
    const double CeeRL_strange() const;
    const double CeeRL_bottom() const;

    const double CeeRR_e() const;
    const double CeeRR_mu() const;
    const double CeeRR_tau() const;
    
    const double CeeRR_up() const;
    const double CeeRR_charm() const;
    const double CeeRR_top() const;
    
    const double CeeRR_down() const;
    const double CeeRR_strange() const;
    const double CeeRR_bottom() const;
    
//  Absolute corrections to the helicity amplitudes squared
  
    const double deltaMLR2_f(const Particle f, const double s) const;
    const double deltaMRL2_f(const Particle f, const double s) const;

    const double deltaMLL2_f(const Particle f, const double s, const double t) const;       
    const double deltaMRR2_f(const Particle f, const double s, const double t) const;

//  Some simple functions for cos \theta integrals   
    const double tovers2(const double cosmin, const double cosmax) const;
    const double uovers2(const double cosmin, const double cosmax) const;

//  Absolute corrections to the differential cross section        
    virtual const double delta_Dsigma_f(const Particle f, const double pol_e, const double pol_p, const double s, const double cos) const;
    
//  Absolute corrections to the differential cross section integrated in [cos \theta_{min},cos \theta_{max}] 
//  Valid for f=/=e
    virtual const double delta_sigma_f(const Particle f, const double pol_e, const double pol_p, const double s, const double cosmin, const double cosmax) const;
    
    virtual const double delta_sigma_had(const double pol_e, const double pol_p, const double s, const double cosmin, const double cosmax) const;
    
//  Total cross sections  (full acceptance)
    virtual const double delta_sigmaTot_f(const Particle f, const double pol_e, const double pol_p, const double s) const;
    
//  Forward-Backward asymmetry (full acceptance). Valid for f!=e
    virtual const double delta_AFB_f(const Particle f, const double pol_e, const double pol_p, const double s) const;
    
//  Expressions for f=e   
    
//  Corrections to t-channel LR and RL amplitudes
    const double deltaMLR2t_e(const double t) const;      
    const double deltaMRL2t_e(const double t) const;    

//  Integrals of the SM squared amplitudes x (t/s)^2, (s/t)^2, (u/s)^2 in [t0, t1]
    virtual const double intMeeLR2SMts2(const double s, const double t0, const double t1) const;
    virtual const double intMeeLRtilde2SMst2(const double s, const double t0, const double t1) const;
    virtual const double intMeeLL2SMus2(const double s, const double t0, const double t1) const;
    virtual const double intMeeRR2SMus2(const double s, const double t0, const double t1) const;
    
//  Integrals of the corrections to the squared amplitudes x (t/s)^2, (s/t)^2, (u/s)^2 in [t0, t1]
    virtual const double intDMLL2eus2(const double s, const double t0, const double t1) const;
    virtual const double intDMRR2eus2(const double s, const double t0, const double t1) const;
    virtual const double intDMLR2ets2(const double s, const double t0, const double t1) const;
    virtual const double intDMRL2ets2(const double s, const double t0, const double t1) const;
    virtual const double intDMLR2etildest2(const double s, const double t0, const double t1) const;
    virtual const double intDMRL2etildest2(const double s, const double t0, const double t1) const;

//  SM cross section integrated in [cos \theta_{min},cos \theta_{max}] 
    virtual const double sigmaSM_ee(const double pol_e, const double pol_p, const double s, const double cosmin, const double cosmax) const;
//  Absolute corrections to the differential cross section integrated in [cos \theta_{min},cos \theta_{max}] 
    virtual const double delta_sigma_ee(const double pol_e, const double pol_p, const double s, const double cosmin, const double cosmax) const;
//  Absolute corrections to the total cross section 
    virtual const double delta_sigmaTot_ee(const double pol_e, const double pol_p, const double s) const;
//  Absolute corrections to the FB asymmetry 
    virtual const double delta_AFB_ee(const double pol_e, const double pol_p, const double s) const;

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
    virtual const double deltaG_hgg() const;
    /**
     * @brief The full new physics contribution to the coupling of the effective interaction @f$H G_{\mu\nu}^AG^{A \mu\nu}@f$,
     * including new local terms and modifications on the SM-loops. Normalized to the SM value.
     * @return @f$\delta g_{HGG}/g_{HGG}^SM}@f$
     */
    virtual const double deltaG_hggRatio() const;
    /**
     * @brief The new physics contribution to the coupling of the effective interaction @f$H W_{\mu\nu}^\dagger W^{\mu\nu}@f$.
     * @return @f$\delta g_{HWW}^{(1)}@f$
     */
    virtual const double deltaG1_hWW() const;
    /**
     * @brief The new physics contribution to the coupling of the effective interaction @f$H W_{\nu}^\dagger \partial^\mu W^{\mu\nu}@f$.
     * @return @f$\delta g_{HWW}^{(2)}@f$
     */
    virtual const double deltaG2_hWW() const;
    /**
     * @brief The new physics contribution to the coupling of the effective interaction @f$H W_{\mu}^\dagger W^{\mu}@f$.
     * @return @f$\delta g_{HWW}^{(3)}@f$
     */
    virtual const double deltaG3_hWW() const;
    /**
     * @brief The new physics contribution to the coupling of the effective interaction @f$H Z_{\mu\nu} Z^{\mu\nu}@f$.
     * @return @f$\delta g_{HZZ}^{(1)}@f$
     */
    virtual const double deltaG1_hZZ() const;
    /**
     * @brief The new physics contribution to the coupling of the effective interaction @f$H Z_{\nu} \partial^\mu Z^{\mu\nu}@f$.
     * @return @f$\delta g_{HZZ}^{(2)}@f$
     */
    virtual const double deltaG2_hZZ() const;
    /**
     * @brief The new physics contribution to the coupling of the effective interaction @f$H Z_{\mu} Z^{\mu}@f$.
     * @return @f$\delta g_{HZZ}^{(3)}@f$
     */
    virtual const double deltaG3_hZZ() const;
    /**
     * @brief The new physics contribution to the coupling of the effective interaction @f$H Z_{\mu\nu} F^{\mu\nu}@f$.
     * @return @f$\delta g_{HZA}^{(1)}@f$
     */
    virtual const double deltaG1_hZA() const;
    /**
     * @brief The full new physics contribution to the coupling of the effective interaction @f$H Z_{\mu\nu} F^{A \mu\nu}@f$,
     * including new local terms and modifications on the SM-loops. Normalized to the SM value.
     * @return @f$\delta g_{HZA}^{(1)}/g_{HZA}^{(1),SM}@f$
     */
    virtual const double deltaG1_hZARatio() const;
    /**
     * @brief The new physics contribution to the coupling of the effective interaction @f$H Z_{\nu} \partial^\mu F^{\mu\nu}@f$.
     * @return @f$\delta g_{HZA}^{(2)}@f$
     */
    virtual const double deltaG2_hZA() const;
    /**
     * @brief The new physics contribution to the coupling of the effective interaction @f$H F_{\mu\nu} F^{\mu\nu}@f$.
     * @return @f$\delta g_{HAA}@f$
     */
    virtual const double deltaG_hAA() const;
    /**
     * @brief The full new physics contribution to the coupling of the effective interaction @f$H F_{\mu\nu} F^{\mu\nu}@f$,
     * including new local terms and modifications on the SM-loops. Normalized to the SM value.
     * @return @f$\delta g_{HAA}/g_{HAA}^SM}@f$
     */
    virtual const double deltaG_hAARatio() const;
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
    virtual const double deltaG_hhhRatio() const;

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
    const double deltaGL_Zffh(const Particle p) const;

    /**
     * @brief The new physics contribution to the coupling of the effective interaction @f$H Z_\mu \bar{f_R}\gamma^mu f_R@f$.
     * @param[in] p a lepton or quark
     * @return @f$\delta g_{ZffH}^{R}@f$
     */
    // no generation mixing 
    const double deltaGR_Zffh(const Particle p) const;

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
    const double deltag3G() const;


    ////////////////////////////////////////////////////////////////////////

    /**
     * @brief Loop function entering in the calculation of the effective @f$Hgg@f$ and @f$H\gamma\gamma@f$ couplings.
     * @param[in] @f$\tau=4 M^2/m_h^2@f$, with @f$M@f$ the mass of the particle in the loop.
     * @return @f$f(\tau)@f$
     */
    gslpp::complex f_triangle(double tau) const;
    /**
     * @brief Loop function entering in the calculation of the effective @f$HZ\gamma@f$ coupling.
     * @param[in] @f$\tau=4 M^2/m_h^2@f$, with @f$M@f$ the mass of the particle in the loop.
     * @return @f$g(\tau)@f$
     */
    gslpp::complex g_triangle(double tau) const;
    /**
     * @brief Loop function entering in the calculation of the effective @f$HZ\gamma@f$ coupling.
     * @param[in] @f$\tau=4 M^2/m_h^2@f$, @f$\lambda=4 M^2/m_Z^2@f$, with @f$M@f$ the mass of the particle in the loop.
     * @return @f$I_1(\tau,\lambda)@f$
     */
    gslpp::complex I_triangle_1(double tau, double lambda) const;
    /**
     * @brief Loop function entering in the calculation of the effective @f$HZ\gamma@f$ coupling.
     * @param[in] @f$\tau=4 M^2/m_h^2@f$, @f$\lambda=4 M^2/m_Z^2@f$, with @f$M@f$ the mass of the particle in the loop.
     * @return @f$I_2(\tau,\lambda)@f$
     */
    gslpp::complex I_triangle_2(double tau, double lambda) const;
    /**
     * @brief Fermionic loop function entering in the calculation of the effective @f$Hgg@f$ and @f$H\gamma\gamma@f$ couplings.
     * @details
     * @f$A^H_f(\tau)=2\tau [1+(1-\tau)f(\tau)]@f$
     * @param[in] @f$\tau=4 M^2/m_h^2@f$, with @f$M@f$ the mass of the fermion in the loop.
     * @return @f$A^H_f(\tau)@f$
     */
    gslpp::complex AH_f(double tau) const;

    /**
     * @brief W loop function entering in the calculation of the effective @f$H\gamma\gamma@f$ coupling.
     * @details
     * @f$A^H_W(\tau)=-[2+3\tau + 3\tau*(2-\tau) f(\tau)]@f$
     * @param[in] @f$\tau=4 M^2/m_h^2@f$, with @f$M@f$ the mass of the fermion in the loop.
     * @return @f$A^H_W(\tau)@f$
     */
    gslpp::complex AH_W(double tau) const;

    /**
     * @brief Fermionic loop function entering in the calculation of the effective @f$HZ\gamma@f$ coupling.
     * @details
     * @param[in] @f$\tau=4 M^2/m_h^2@f$, @f$\lambda=4 M^2/m_Z^2@f$, with @f$M@f$ the mass of the fermion in the loop.
     * @return @f$A^{HZ\gamma}_f(\tau,\lambda)@f$
     */
    gslpp::complex AHZga_f(double tau, double lambda) const;

    /**
     * @brief W loop function entering in the calculation of the effective @f$HZ\gamma@f$ coupling.
     * @details
     * @param[in] @f$\tau=4 M^2/m_h^2@f$, @f$\lambda=4 M^2/m_Z^2@f$, with @f$M@f$ the mass of the fermion in the loop.
     * @return @f$A^{HZ\gamma}_W(\tau,\lambda)@f$
     */
    gslpp::complex AHZga_W(double tau, double lambda) const;

    /**
     * @brief The ratio @f$\mu_{ggH}@f$ between the gluon-gluon fusion Higgs
     * production cross-section in the current model and in the Standard Model.
     * @param[in] sqrt_s the center-of-mass energy in TeV
     * @return @f$\mu_{ggH}@f$
     */
    virtual const double muggH(double sqrt_s) const;
    /**
     * @brief The ratio @f$\mu_{ggHH}@f$ between the gluon-gluon fusion di-Higgs
     * production cross-section in the current model and in the Standard Model.
     * (From arXiv: 1502.00539 [hpe-ph].)
     * @param[in] sqrt_s the center-of-mass energy in TeV
     * @return @f$\mu_{ggHH}@f$
     */
    virtual const double muggHH(double sqrt_s) const;
    /**
     * @brief The ratio @f$\mu_{VBF}@f$ between the vector-boson fusion Higgs
     * production cross-section in the current model and in the Standard Model. 
     * @param[in] sqrt_s the center-of-mass energy in TeV
     * @return @f$\mu_{VBF}@f$
     */
    virtual const double muVBF(double sqrt_s) const;
    /**
     * @brief The ratio @f$\mu_{VBF+\gamma}@f$ between the vector-boson fusion Higgs
     * production cross-section in association with a hard photon in the current model
     * and in the Standard Model. 
     * @param[in] sqrt_s the center-of-mass energy in TeV
     * @return @f$\mu_{VBF+\gamma}@f$
     */
    virtual const double muVBFgamma(double sqrt_s) const;
    /**
     * @brief The ratio @f$\mu_{eeWBF}@f$ between the 
     * @f$ e^{+}e^{-}\to \nu\bar{\nu} H @f$ production
     * cross-section in the current model and in the Standard Model.
     * @param[in] sqrt_s the center-of-mass energy in TeV, Pol_em and Pol_ep
     * are the polarization of electrons and positrons, respectively
     * @return @f$\mu_{eeWBF}@f$
     */
    virtual const double mueeWBF(double sqrt_s, double Pol_em, double Pol_ep) const;
    /**
     * @brief The ratio @f$\mu_{e^+e^- \to H\nu\bar{\nu}}@f$ between the 
     * @f$ e^+e^- \to H\nu\bar{\nu} @f$ associated production
     * cross-section in the current model and in the Standard Model.
     * @param[in] sqrt_s the center-of-mass energy in TeV, Pol_em and Pol_ep
     * are the polarization of electrons and positrons, respectively
     * @return @f$\mu_{e^+e^- \to H\nu\bar{\nu}}@f$
     */
    virtual const double mueeHvv(double sqrt_s, double Pol_em, double Pol_ep) const;
    /**
     * @brief The ratio @f$\mu_{e^+e^- \to H\nu\bar{\nu}}@f$ between the 
     * @f$ e^+e^- \to H\nu\bar{\nu} @f$ associated production
     * cross-section in the current model and in the Standard Model.
     * @param[in] sqrt_s the center-of-mass energy in TeV, Pol_em and Pol_ep
     * are the polarization of electrons and positrons, respectively 
     * @return @f$\mu_{e^+e^- \to H\nu\bar{\nu}}@f$
     */
    virtual const double mueeHvvPol(double sqrt_s, double Pol_em, double Pol_ep) const;
    /**
     * @brief The ratio @f$\mu_{eeZBF}@f$ between the 
     * @f$ e^{+}e^{-}\to e^{+}e^{-} H @f$ production
     * cross-section in the current model and in the Standard Model.
     * @param[in] sqrt_s the center-of-mass energy in TeV, Pol_em and Pol_ep
     * are the polarization of electrons and positrons, respectively
     * @return @f$\mu_{eeZBF}@f$
     */
    virtual const double mueeZBF(double sqrt_s, double Pol_em, double Pol_ep) const;
    /**
     * @brief The ratio @f$\mu_{eeZBF}@f$ between the 
     * @f$ e^{+}e^{-}\to e^{+}e^{-} H @f$ production
     * cross-section in the current model and in the Standard Model.
     * @param[in] sqrt_s the center-of-mass energy in TeV, Pol_em and Pol_ep
     * are the polarization of electrons and positrons, respectively 
     * @return @f$\mu_{eeZBF}@f$
     */
    virtual const double mueeZBFPol(double sqrt_s, double Pol_em, double Pol_ep) const;
    /**
     * @brief The ratio @f$\mu_{epWBF}@f$ between the 
     * @f$ e^{-} p\to \nu j H @f$ production
     * cross-section in the current model and in the Standard Model.
     * @param[in] sqrt_s the center-of-mass energy in TeV
     * @return @f$\mu_{epWBF}@f$
     */
    virtual const double muepWBF(double sqrt_s) const;
    /**
     * @brief The ratio @f$\mu_{epZBF}@f$ between the 
     * @f$ e^{-} p\to e^{-} j H @f$ production
     * cross-section in the current model and in the Standard Model.
     * @param[in] sqrt_s the center-of-mass energy in TeV
     * @return @f$\mu_{epZBF}@f$
     */
    virtual const double muepZBF(double sqrt_s) const;
    /**
     * @brief The ratio @f$\mu_{WH}@f$ between the W-Higgs associated production
     * cross-section in the current model and in the Standard Model. 
     * @param[in] sqrt_s the center-of-mass energy in TeV
     * @return @f$\mu_{WH}@f$
     */
    virtual const double muWH(double sqrt_s) const;
    /**
     * @brief The ratio @f$\mu_{WH}@f$ between the W-Higgs associated production
     * cross-section in the current model and in the Standard Model, with @f$p_{T,H}>250@f$ GeV.
     * @param[in] sqrt_s the center-of-mass energy in TeV
     * @return @f$\mu_{WH}@f$
     */
    virtual const double muWHpT250(double sqrt_s) const;
    /**
     * @brief The ratio @f$\mu_{ZH}@f$ between the Z-Higgs associated production
     * cross-section in the current model and in the Standard Model.
     * @param[in] sqrt_s the center-of-mass energy in TeV
     * @return @f$\mu_{ZH}@f$
     */
    virtual const double muZH(double sqrt_s) const;
    /**
     * @brief The ratio @f$\mu_{ZH}@f$ between the Z-Higgs associated production
     * cross-section in the current model and in the Standard Model, with @f$p_{T,H}>250@f$ GeV.
     * @param[in] sqrt_s the center-of-mass energy in TeV
     * @return @f$\mu_{ZH}@f$
     */
    virtual const double muZHpT250(double sqrt_s) const;
    /**
     * @brief The ratio @f$\mu_{eeZH}@f$ between the 
     * @f$e^{+}e^{-}\to ZH@f$ associated production
     * cross-section in the current model and in the Standard Model.
     * @param[in] sqrt_s the center-of-mass energy in TeV, Pol_em and Pol_ep
     * are the polarization of electrons and positrons, respectively 
     * @return @f$\mu_{eeZH}@f$
     */
    virtual const double mueeZH(double sqrt_s, const double Pol_em, const double Pol_ep) const;
    /**
     * @brief The ratio @f$\mu_{eeZH, Z \to e^+ e^-, \mu^+ \mu^-}@f$ between the 
     * @f$ e^{+}e^{-}\to ZH, Z \to e^+ e^-, \mu^+ \mu^- @f$ associated production
     * cross-section in the current model and in the Standard Model.
     * @param[in] sqrt_s the center-of-mass energy in TeV, Pol_em and Pol_ep
     * are the polarization of electrons and positrons, respectively
     * @return @f$\mu_{eeZH, Z \to e^+ e^-, \mu^+ \mu^-}@f$
     */
    virtual const double mueeZllH(double sqrt_s, double Pol_em, double Pol_ep) const;
    /**
     * @brief The ratio @f$\mu_{eeZH, Z \to q \bar{q}}@f$ between the 
     * @f$ e^{+}e^{-}\to ZH, Z \to q \bar{q} @f$ associated production
     * cross-section in the current model and in the Standard Model.
     * @param[in] sqrt_s the center-of-mass energy in TeV, Pol_em and Pol_ep
     * are the polarization of electrons and positrons, respectively
     * @return @f$\mu_{eeZH, Z \to q \bar{q}}@f$
     */
    virtual const double mueeZqqH(double sqrt_s, double Pol_em, double Pol_ep) const;
    /**
     * @brief The ratio @f$\mu_{eeZH}@f$ between the 
     * @f$ e^{+}e^{-}\to ZH @f$ associated production
     * cross-section in the current model and in the Standard Model.
     * @param[in] sqrt_s the center-of-mass energy in TeV, Pol_em and Pol_ep
     * are the polarization of electrons and positrons, respectively 
     * @return @f$\mu_{eeZH}@f$
     */
    virtual const double mueeZHPol(double sqrt_s, double Pol_em, double Pol_ep) const;
    /**
     * @brief The ratio @f$\mu_{eeZH, Z \to e^+ e^-, \mu^+ \mu^-}@f$ between the 
     * @f$ e^{+}e^{-}\to ZH, Z \to e^+ e^-, \mu^+ \mu^- @f$ associated production
     * cross-section in the current model and in the Standard Model.
     * @param[in] sqrt_s the center-of-mass energy in TeV, Pol_em and Pol_ep
     * are the polarization of electrons and positrons, respectively 
     * @return @f$\mu_{eeZH, Z \to e^+ e^-, \mu^+ \mu^-}@f$
     */
    virtual const double mueeZllHPol(double sqrt_s, double Pol_em, double Pol_ep) const;
    /**
     * @brief The ratio @f$\mu_{eeZH, Z \to q \bar{q}}@f$ between the 
     * @f$ e^{+}e^{-}\to ZH, Z \to q \bar{q} @f$ associated production
     * cross-section in the current model and in the Standard Model.
     * @param[in] sqrt_s the center-of-mass energy in TeV, Pol_em and Pol_ep
     * are the polarization of electrons and positrons, respectively 
     * @return @f$\mu_{eeZH, Z \to q \bar{q}}@f$
     */
    virtual const double mueeZqqHPol(double sqrt_s, double Pol_em, double Pol_ep) const;
    /**
     * @brief the angular parameter @f$a@f$ from 
     * @f$\mu_{e^+e^- \to ZH}@f$ (arXiv:1708.09079 [hep-ph]).
     * @param[in] sqrt_s the center-of-mass energy in TeV, Pol_em and Pol_ep
     * are the polarization of electrons and positrons, respectively 
     * @return @f$a_{eeZH}@f$
     */
    virtual const double aPskPol(double sqrt_s, double Pol_em, double Pol_ep) const;
    /**
     * @brief the angular parameter @f$b@f$ from 
     * @f$\mu_{e^+e^- \to ZH}@f$ (arXiv:1708.09079 [hep-ph]).
     * @param[in] sqrt_s the center-of-mass energy in TeV, Pol_em and Pol_ep
     * are the polarization of electrons and positrons, respectively 
     * @return @f$b_{eeZH}@f$
     */
    virtual const double bPskPol(double sqrt_s, double Pol_em, double Pol_ep) const;
    /**
     * @brief The ratio @f$\mu_{VH}@f$ between the WH+ZH associated production
     * cross-section in the current model and in the Standard Model.
     * @param[in] sqrt_s the center-of-mass energy in TeV
     * @return @f$\mu_{VH}@f$
     */
    virtual const double muVH(double sqrt_s) const;
    /**
     * @brief The ratio @f$\mu_{VH}@f$ between the WH+ZH associated production
     * cross-section in the current model and in the Standard Model, with @f$p_{T,H}>250@f$ GeV.
     * @param[in] sqrt_s the center-of-mass energy in TeV
     * @return @f$\mu_{VH}@f$
     */
    virtual const double muVHpT250(double sqrt_s) const;
    /**
     * @brief The ratio @f$\mu_{VBF+VH}@f$ between the sum of VBF and WH+ZH associated production
     * cross-section in the current model and in the Standard Model.
     * @param[in] sqrt_s the center-of-mass energy in TeV
     * @return @f$\mu_{VBF+VH}@f$
     */
    virtual const double muVBFpVH(double sqrt_s) const;
    /**
     * @brief The ratio @f$\mu_{ttH}@f$ between the t-tbar-Higgs associated 
     * production cross-section in the current model and in the Standard Model.
     * @param[in] sqrt_s the center-of-mass energy in TeV
     * @return @f$\mu_{ttH}@f$
     */
    virtual const double muttH(double sqrt_s) const;
    /**
     * @brief The ratio @f$\mu_{tHq}@f$ between the t-q-Higgs associated 
     * production cross-section in the current model and in the Standard Model.
     * @param[in] sqrt_s the center-of-mass energy in TeV
     * @return @f$\mu_{tHq}@f$
     */
    virtual const double mutHq(double sqrt_s) const;
    /**
     * @brief The ratio @f$\mu_{ggH+ttH}@f$ between the sum of gluon-gluon fusion
     * and t-tbar-Higgs associated 
     * production cross-section in the current model and in the Standard Model.
     * @param[in] sqrt_s the center-of-mass energy in TeV
     * @return @f$\mu_{ggH+ttH}@f$
     */
    virtual const double muggHpttH(double sqrt_s) const;
    /**
     * @brief The ratio @f$\mu_{eettH}@f$ between the 
     * @f$ e^{+}e^{-}\to t\bar{t} H @f$ production
     * cross-section in the current model and in the Standard Model.
     * @param[in] sqrt_s the center-of-mass energy in TeV, Pol_em and Pol_ep
     * are the polarization of electrons and positrons, respectively
     * @return @f$\mu_{eettH}@f$
     */
    virtual const double mueettH(double sqrt_s, double Pol_em, double Pol_ep) const;
    /**
     * @brief The ratio @f$\mu_{eettH}@f$ between the 
     * @f$ e^{+}e^{-}\to t\bar{t} H @f$ production
     * cross-section in the current model and in the Standard Model.
     * @param[in] sqrt_s the center-of-mass energy in TeV, Pol_em and Pol_ep
     * are the polarization of electrons and positrons, respectively
     * @return @f$\mu_{eettH}@f$
     */
    virtual const double mueettHPol(double sqrt_s, double Pol_em, double Pol_ep) const;
    /**
     * @brief The ratio @f$\mu_{\mu\mu H}@f$ between the @f$\sigma(\mu \mu \to H)}@f$
     * production cross-section in the current model and in the Standard Model.
     * @param[in] sqrt_s the center-of-mass energy in TeV
     * @return @f$\mu_{\mu\mu H}@f$
     */
    virtual const double mummH(double sqrt_s) const;
    /**
     * @brief The ratio @f$\mu_{\mu\mu H}@f$ between the @f$\sigma(\mu \mu \to H)}@f$
     * production cross-section in the current model and in the Standard Model, 
     * in the narrow width approximation.
     * @param[in] sqrt_s the center-of-mass energy in TeV
     * @return @f$\mu_{\mu\mu H}@f$
     */
    virtual const double mummHNWA(double sqrt_s) const;
    /**
     * @brief The ratio @f$\mu_{\mu\mu ZH}@f$ between the @f$\sigma(\mu \mu \to Z H)}@f$
     * production cross-section in the current model and in the Standard Model.
     * @param[in] sqrt_s the center-of-mass energy in TeV
     * @return @f$\mu_{\mu\mu ZH}@f$
     */
    virtual const double mummZH(double sqrt_s) const;
    /**
     * @brief The ratio @f$\mu_{\mu\mu H\nu\nu}@f$ between the @f$\sigma(\mu \mu \to H \nu \nu)}@f$
     * production cross-section in the current model and in the Standard Model.
     * @param[in] sqrt_s the center-of-mass energy in TeV
     * @return @f$\mu_{\mu\mu H\nu\nu}@f$
     */
    virtual const double mummHvv(double sqrt_s) const;
    /**
     * @brief The ratio @f$\mu_{\mu\mu H\mu\mu}@f$ between the @f$\sigma(\mu \mu \to H \mu \mu)}@f$
     * production cross-section in the current model and in the Standard Model.
     * @param[in] sqrt_s the center-of-mass energy in TeV
     * @return @f$\mu_{\mu\mu H\mu\mu}@f$
     */
    virtual const double mummHmm(double sqrt_s) const;
    /**
     * @brief The ratio @f$\mu_{\mu\mu ttH}@f$ between the @f$\sigma(\mu \mu \to t\bar{t} H )}@f$
     * production cross-section in the current model and in the Standard Model.
     * @param[in] sqrt_s the center-of-mass energy in TeV
     * @return @f$\mu_{\mu\mu ttH}@f$
     */
    virtual const double mummttH(double sqrt_s) const;
    
    //AG:added, begin
    virtual const double delta_muggH_1(const double sqrt_s) const;    ///<The SMEFT linear correction to the ratio @f$\mu_{ggH}@f$ between the gluon-gluon fusion Higgs production cross-section in the current model and in the Standard Model.
    virtual const double delta_muggH_2(const double sqrt_s) const
    {
        return 0.0;
    }	///<The SMEFT quadratic correction to the ratio @f$\mu_{ggH}@f$ between the gluon-gluon fusion Higgs production cross-section in the current model and in the Standard Model.
    
    virtual const double delta_muVBF_1(const double sqrt_s) const;	///<The SMEFT linear correction to the ratio @f$\mu_{VBF}@f$ between the vector-boson fusion Higgs production cross-section in the current model and in the Standard Model.
    virtual const double delta_muVBF_2(const double sqrt_s) const
    {
        return 0.0;
    }	///<The SMEFT quadratic correction to the ratio @f$\mu_{VBF}@f$ between the vector-boson fusion Higgs production cross-section in the current model and in the Standard Model.
    
    virtual const double delta_muWH_1(const double sqrt_s) const;	///<The SMEFT linear correction to the ratio @f$\mu_{WH}@f$ between the W-Higgs associated production cross-section in the current model and in the Standard Model.
    virtual const double delta_muWH_2(const double sqrt_s) const
    {
        return 0.0;
    }	///<The SMEFT quadratic correction to the ratio @f$\mu_{WH}@f$ between the W-Higgs associated production cross-section in the current model and in the Standard Model.
    
    virtual const double delta_muZH_1(const double sqrt_s) const;	///<The SMEFT linear correction to the ratio @f$\mu_{ZH}@f$ between the Z-Higgs associated production cross-section in the current model and in the Standard Model.
    virtual const double delta_muZH_2(const double sqrt_s) const
    {
        return 0.0;
    }	///<The SMEFT quadratic correction to the ratio @f$\mu_{ZH}@f$ between the Z-Higgs associated production cross-section in the current model and in the Standard Model.
    
    virtual const double delta_muVH_1(const double sqrt_s) const;	///<The SMEFT linear correction to the ratio @f$\mu_{VH}@f$ between the Z-Higgs and W-Higgs associated production cross-section in the current model and in the Standard Model.
    virtual const double delta_muVH_2(const double sqrt_s) const
    {
        return 0.0;
    }	///<The SMEFT quadratic correction to the ratio @f$\mu_{VH}@f$ between the Z-Higgs and W-Higgs associated production cross-section in the current model and in the Standard Model.
    
    virtual const double delta_muttH_1(const double sqrt_s) const;	///<The SMEFT linear correction to the ratio @f$\mu_{ttH}@f$ between the t-tbar-Higgs associated production cross-section in the current model and in the Standard Model.
    virtual const double delta_muttH_2(const double sqrt_s) const
    {
        return 0.0;
    }	///<The SMEFT quadratic correction to the ratio @f$\mu_{ttH}@f$ between the t-tbar-Higgs associated production cross-section in the current model and in the Standard Model.
    
    virtual const double delta_mutH_1(const double sqrt_s) const
    {
        return 0.0;
    }	///<The SMEFT linear correction to the ratio @f$\mu_{tH}@f$ between the t-Higgs associated production cross-section in the current model and in the Standard Model.
    
    virtual const double delta_mutH_2(const double sqrt_s) const
    {
        return 0.0;
    }	///<The SMEFT quadratic correction to the ratio @f$\mu_{tH}@f$ between the t-Higgs associated production cross-section in the current model and in the Standard Model.
    
    virtual const double delta_mubbH_1(const double sqrt_s) const
    {
        return 0.0;
    }	///<The SMEFT linear correction to the ratio @f$\mu_{bbH}@f$ between the bbH production cross-section in the current model and in the Standard Model.
    virtual const double delta_mubbH_2(const double sqrt_s) const
    {
        return 0.0;
    }	///<The SMEFT quadratic correction to the ratio @f$\mu_{bbH}@f$ between the bbH production cross-section in the current model and in the Standard Model.
    //AG:end

    ////////////HIGGS DECAY WIDTHS AND BRANCHING RATIOS/////////////
    /**
     * @brief The ratio of the @f$\Gamma(H)@f$ in the current model
     * and in the Standard Model.
     * @return @f$\Gamma(H)@f$/@f$\Gamma(H)_{\mathrm{SM}}@f$
     */
    virtual const double computeGammaTotalRatio() const;
    /**
     * @brief The new physics contribution to the ratio of the @f$\Gamma(H)@f$ 
     * in the current model and in the Standard Model. Only terms that are
     * linear in the effective Lagrangian coefficients.
     * @return @f$\delta \Gamma(H)@f$/@f$\Gamma(H)_{\mathrm{SM}}@f$
     */
    virtual const double deltaGammaTotalRatio1() const;
    /**
     * @brief The new physics contribution to the ratio of the @f$\Gamma(H)@f$ 
     * in the current model and in the Standard Model. Only terms that are
     * linear in the effective Lagrangian coefficients. Neglecting SM theory errors.
     * @return @f$\delta \Gamma(H)@f$/@f$\Gamma(H)_{\mathrm{SM}}@f$
     */
    virtual const double deltaGammaTotalRatio1noError() const;
    /**
     * @brief The new physics contribution to the ratio of the @f$\Gamma(H)@f$ 
     * in the current model and in the Standard Model. Only terms that are
     * quadratic in the effective Lagrangian coefficients.
     * @return @f$\delta \Gamma(H)@f$/@f$\Gamma(H)_{\mathrm{SM}}@f$
     */
    virtual const double deltaGammaTotalRatio2() const;

    /**
     * @brief The ratio of the @f$\Gamma(H\to gg)@f$ in the current model
     * and in the Standard Model.
     * @return @f$\Gamma(H\to gg)@f$/@f$\Gamma(H\to gg)_{\mathrm{SM}}@f$
     */
    const double GammaHggRatio() const;
    /**
     * @brief The new physics contribution to the ratio of the @f$\Gamma(H\to gg)@f$
     * in the current model and in the Standard Model. Only terms that are
     * linear in the effective Lagrangian coefficients.
     * @return @f$\delta \Gamma(H\to gg)@f$/@f$\Gamma(H\to gg)_{\mathrm{SM}}@f$
     */
    const double deltaGammaHggRatio1() const;
    /**
     * @brief The new physics contribution to the ratio of the @f$\Gamma(H\to gg)@f$
     * in the current model and in the Standard Model. (Only terms that are
     * quadratic in the effective Lagrangian coefficients.)
     * @return @f$\delta \Gamma(H\to gg)@f$/@f$\Gamma(H\to gg)_{\mathrm{SM}}@f$
     */
    const double deltaGammaHggRatio2() const;
    /**
     * @brief The ratio of the Br@f$(H\to gg)@f$ in the current model
     * and in the Standard Model.
     * @return Br@f$(H\to gg)@f$/Br@f$(H\to gg)_{\mathrm{SM}}@f$
     */
    virtual const double BrHggRatio() const;

    /**
     * @brief The ratio of the @f$\Gamma(H\to WW)@f$ in the current model
     * and in the Standard Model.
     * @return @f$\Gamma(H\to WW)@f$/@f$\Gamma(H\to WW)_{\mathrm{SM}}@f$
     */
    const double GammaHWWRatio() const;
    /**
     * @brief The new physics contribution to the ratio of the @f$\Gamma(H\to WW)@f$
     * in the current model and in the Standard Model. (Only terms that are
     * linear in the effective Lagrangian coefficients.)
     * @return @f$\delta \Gamma(H\to WW)@f$/@f$\Gamma(H\to WW)_{\mathrm{SM}}@f$
     */
    const double deltaGammaHWWRatio1() const;
    /**
     * @brief The new physics contribution to the ratio of the @f$\Gamma(H\to WW)@f$
     * in the current model and in the Standard Model. (Only terms that are
     * quadratic in the effective Lagrangian coefficients.)
     * @return @f$\delta \Gamma(H\to WW)@f$/@f$\Gamma(H\to WW)_{\mathrm{SM}}@f$
     */
    const double deltaGammaHWWRatio2() const;
    /**
     * @brief The ratio of the Br@f$(H\to WW)@f$ in the current model
     * and in the Standard Model.
     * @return Br@f$(H\to WW)@f$/Br@f$(H\to WW)_{\mathrm{SM}}@f$
     */
    virtual const double BrHWWRatio() const;

    /**
     * @brief The ratio of the @f$\Gamma(H\to WW^*\to 4f)@f$, with @f$f@f$ any fermion, in the current model
     * and in the Standard Model.
     * @return @f$\Gamma(H\to WW^*\to 4f)@f$/@f$\Gamma(H\to WW^*\to 4f)_{\mathrm{SM}}@f$
     */
    const double GammaHWW4fRatio() const;
    /**
     * @brief The new physics contribution to the ratio of the @f$\Gamma(H\to WW^*\to 4f)@f$, with @f$f@f$ any fermion,
     * in the current model and in the Standard Model. (Only terms that are
     * linear in the effective Lagrangian coefficients.)
     * @return @f$\delta \Gamma(H\to WW^*\to 4f)@f$/@f$\Gamma(H\to WW^*\to 4f)_{\mathrm{SM}}@f$
     */
    const double deltaGammaHWW4fRatio1() const;
    /**
     * @brief The new physics contribution to the ratio of the @f$\Gamma(H\to WW^*\to 4f)@f$, with @f$f@f$ any fermion,
     * in the current model and in the Standard Model. (Only terms that are
     * quadratic in the effective Lagrangian coefficients.)
     * @return @f$\delta \Gamma(H\to WW^*\to 4f)@f$/@f$\Gamma(H\to WW^*\to 4f)_{\mathrm{SM}}@f$
     */
    const double deltaGammaHWW4fRatio2() const;
    /**
     * @brief The ratio of the Br@f$(H\to WW^*\to 4f)@f$, with @f$f@f$ any fermion, in the current model
     * and in the Standard Model.
     * @return Br@f$(H\to WW^*\to 4f)@f$/Br@f$(H\to WW^*\to 4f)_{\mathrm{SM}}@f$
     */
    virtual const double BrHWW4fRatio() const;

    /**
     * @brief The ratio of the @f$\Gamma(H\to ZZ)@f$ in the current model
     * and in the Standard Model.
     * @return @f$\Gamma(H\to ZZ)@f$/@f$\Gamma(H\to ZZ)_{\mathrm{SM}}@f$
     */
    const double GammaHZZRatio() const;
    /**
     * @brief The new physics contribution to the ratio of the @f$\Gamma(H\to ZZ)@f$
     * in the current model and in the Standard Model. (Only terms that are
     * linear in the effective Lagrangian coefficients.)
     * @return @f$\delta \Gamma(H\to ZZ)@f$/@f$\Gamma(H\to ZZ)_{\mathrm{SM}}@f$
     */
    const double deltaGammaHZZRatio1() const;
    /**
     * @brief The new physics contribution to the ratio of the @f$\Gamma(H\to ZZ)@f$
     * in the current model and in the Standard Model. (Only terms that are
     * quadratic in the effective Lagrangian coefficients.)
     * @return @f$\delta \Gamma(H\to ZZ)@f$/@f$\Gamma(H\to ZZ)_{\mathrm{SM}}@f$
     */
    const double deltaGammaHZZRatio2() const;
    /**
     * @brief The ratio of the Br@f$(H\to ZZ)@f$ in the current model
     * and in the Standard Model.
     * @return Br@f$(H\to ZZ)@f$/Br@f$(H\to ZZ)_{\mathrm{SM}}@f$
     */
    virtual const double BrHZZRatio() const;

    /**
     * @brief The ratio of the @f$\Gamma(H\to ZZ* \to 4f)@f$, with @f$f@f$ any fermion, in the current model
     * and in the Standard Model.
     * @return @f$\Gamma(H\to ZZ* \to 4f)@f$/@f$\Gamma(H\to ZZ* \to 4f)_{\mathrm{SM}}@f$
     */
    const double GammaHZZ4fRatio() const;
    /**
     * @brief The new physics contribution to the ratio of the @f$\Gamma(H\to ZZ* \to 4f)@f$, with @f$f@f$ any fermion,
     * in the current model and in the Standard Model. (Only terms that are
     * linear in the effective Lagrangian coefficients.)
     * @return @f$\delta \Gamma(H\to ZZ* \to 4f)@f$/@f$\Gamma(H\to ZZ* \to 4f)_{\mathrm{SM}}@f$
     */
    const double deltaGammaHZZ4fRatio1() const;
    /**
     * @brief The new physics contribution to the ratio of the @f$\Gamma(H\to ZZ* \to 4f)@f$, with @f$f@f$ any fermion,
     * in the current model and in the Standard Model. (Only terms that are
     * quadratic in the effective Lagrangian coefficients.)
     * @return @f$\delta \Gamma(H\to ZZ* \to 4f)@f$/@f$\Gamma(H\to ZZ* \to 4f)_{\mathrm{SM}}@f$
     */
    const double deltaGammaHZZ4fRatio2() const;
    /**
     * @brief The ratio of the Br@f$(H\to ZZ* \to 4f)@f$, with @f$f@f$ any fermion, in the current model
     * and in the Standard Model.
     * @return Br@f$(H\to ZZ* \to 4f)@f$/Br@f$(H\to ZZ* \to 4f)_{\mathrm{SM}}@f$
     */
    virtual const double BrHZZ4fRatio() const;

    /**
     * @brief The ratio of the Br@f$(H\to VV)@f$ in the current model
     * and in the Standard Model.
     * @return Br@f$(H\to VV)@f$/Br@f$(H\to VV)_{\mathrm{SM}}@f$
     */
    virtual const double BrHVVRatio() const;

    /**
     * @brief The ratio of the @f$\Gamma(H\to Z\gamma)@f$ in the current model
     * and in the Standard Model.
     * @return @f$\Gamma(H\to Z\gamma)@f$/@f$\Gamma(H\to Z\gamma)_{\mathrm{SM}}@f$
     */
    const double GammaHZgaRatio() const;
    /**
     * @brief The new physics contribution to the ratio of the @f$\Gamma(H\to Z\gamma)@f$
     * in the current model and in the Standard Model. (Only terms that are
     * linear in the effective Lagrangian coefficients.)
     * @return @f$\delta \Gamma(H\to Z\gamma)@f$/@f$\Gamma(H\to Z\gamma)_{\mathrm{SM}}@f$
     */
    const double deltaGammaHZgaRatio1() const;
    /**
     * @brief The new physics contribution to the ratio of the @f$\Gamma(H\to Z\gamma)@f$
     * in the current model and in the Standard Model. (Only terms that are
     * quadratic in the effective Lagrangian coefficients.)
     * @return @f$\delta \Gamma(H\to Z\gamma)@f$/@f$\Gamma(H\to Z\gamma)_{\mathrm{SM}}@f$
     */
    const double deltaGammaHZgaRatio2() const;
    /**
     * @brief The ratio of the Br@f$(H\to Z\gamma)@f$ in the current model
     * and in the Standard Model.
     * @return Br@f$(H\to Z\gamma)@f$/Br@f$(H\to Z\gamma)_{\mathrm{SM}}@f$
     */
    virtual const double BrHZgaRatio() const;
    /**
     * @brief The ratio of the Br@f$(H\to Z\gamma\to ll\gamma)@f$ (@f$l=e,\mu @f$) in the current model
     * and in the Standard Model.
     * @return Br@f$(H\to Z\gamma\to ll\gamma)@f$/Br@f$(H\to Z\gamma\to ll\gamma)_{\mathrm{SM}}@f$
     */
    virtual const double BrHZgallRatio() const;
    /**
     * @brief The ratio of the Br@f$(H\to Z\gamma\to ee\gamma)@f$ in the current model
     * and in the Standard Model.
     * @return Br@f$(H\to Z\gamma\to ee\gamma)@f$/Br@f$(H\to Z\gamma\to ee\gamma)_{\mathrm{SM}}@f$
     */
    virtual const double BrHZgaeeRatio() const;
    /**
     * @brief The ratio of the Br@f$(H\to Z\gamma\to \mu\mu\gamma)@f$ in the current model
     * and in the Standard Model.
     * @return Br@f$(H\to Z\gamma\to \mu\mu\gamma)@f$/Br@f$(H\to Z\gamma\to \mu\mu\gamma)_{\mathrm{SM}}@f$
     */
    virtual const double BrHZgamumuRatio() const;

    /**
     * @brief The ratio of the @f$\Gamma(H\to \gamma\gamma)@f$ in the current model
     * and in the Standard Model.
     * @return @f$\Gamma(H\to \gamma\gamma)@f$/@f$\Gamma(H\to \gamma\gamma)_{\mathrm{SM}}@f$
     */
    const double GammaHgagaRatio() const;
    /**
     * @brief The new physics contribution to the ratio of the @f$\Gamma(H\to \gamma\gamma)@f$
     * in the current model and in the Standard Model. (Only terms that are
     * linear in the effective Lagrangian coefficients.)
     * @return @f$\delta \Gamma(H\to \gamma\gamma)@f$/@f$\Gamma(H\to \gamma\gamma)_{\mathrm{SM}}@f$
     */
    const double deltaGammaHgagaRatio1() const;
    /**
     * @brief The new physics contribution to the ratio of the @f$\Gamma(H\to \gamma\gamma)@f$
     * in the current model and in the Standard Model. (Only terms that are
     * quadratic in the effective Lagrangian coefficients.)
     * @return @f$\delta \Gamma(H\to \gamma\gamma)@f$/@f$\Gamma(H\to \gamma\gamma)_{\mathrm{SM}}@f$
     */
    const double deltaGammaHgagaRatio2() const;
    /**
     * @brief The ratio of the Br@f$(H\to \gamma\gamma)@f$ in the current model
     * and in the Standard Model.
     * @return Br@f$(H\to \gamma\gamma)@f$/Br@f$(H\to \gamma\gamma)_{\mathrm{SM}}@f$
     */
    virtual const double BrHgagaRatio() const;

    /////////////////////// HIGGS TO 2 FERMION DECAYS /////////////////////////

    /**
     * @brief The ratio of the @f$\Gamma(H\to \mu\mu)@f$ in the current model
     * and in the Standard Model.
     * @return @f$\Gamma(H\to \mu\mu)@f$/@f$\Gamma(H\to \mu\mu)_{\mathrm{SM}}@f$
     */
    const double GammaHmumuRatio() const;
    /**
     * @brief The new physics contribution to the ratio of the @f$\Gamma(H\to \mu\mu)@f$
     * in the current model and in the Standard Model. (Only terms that are
     * linear in the effective Lagrangian coefficients.)
     * @return @f$\delta \Gamma(H\to \mu\mu)@f$/@f$\Gamma(H\to \mu\mu)_{\mathrm{SM}}@f$
     */
    const double deltaGammaHmumuRatio1() const;
    /**
     * @brief The new physics contribution to the ratio of the @f$\Gamma(H\to \mu\mu)@f$
     * in the current model and in the Standard Model. (Only terms that are
     * quadratic in the effective Lagrangian coefficients.)
     * @return @f$\delta \Gamma(H\to \mu\mu)@f$/@f$\Gamma(H\to \mu\mu)_{\mathrm{SM}}@f$
     */
    const double deltaGammaHmumuRatio2() const;
    /**
     * @brief The ratio of the Br@f$(H\to \mu^+\mu^-)@f$ in the current model
     * and in the Standard Model.
     * @return Br@f$(H\to \mu^+\mu^-)@f$/Br@f$(H\to \mu^+\mu^-)_{\mathrm{SM}}@f$
     */
    virtual const double BrHmumuRatio() const;

    /**
     * @brief The ratio of the @f$\Gamma(H\to \tau\tau)@f$ in the current model
     * and in the Standard Model.
     * @return @f$\Gamma(H\to \tau\tau)@f$/@f$\Gamma(H\to \tau\tau)_{\mathrm{SM}}@f$
     */
    const double GammaHtautauRatio() const;
    /**
     * @brief The new physics contribution to the ratio of the @f$\Gamma(H\to \tau\tau)@f$
     * in the current model and in the Standard Model. (Only terms that are
     * linear in the effective Lagrangian coefficients.)
     * @return @f$\delta \Gamma(H\to \tau\tau)@f$/@f$\Gamma(H\to \tau\tau)_{\mathrm{SM}}@f$
     */
    const double deltaGammaHtautauRatio1() const;
    /**
     * @brief The new physics contribution to the ratio of the @f$\Gamma(H\to \tau\tau)@f$
     * in the current model and in the Standard Model. (Only terms that are
     * quadratic in the effective Lagrangian coefficients.)
     * @return @f$\delta \Gamma(H\to \tau\tau)@f$/@f$\Gamma(H\to \tau\tau)_{\mathrm{SM}}@f$
     */
    const double deltaGammaHtautauRatio2() const;
    /**
     * @brief The ratio of the Br@f$(H\to \tau^+\tau^-)@f$ in the current model
     * and in the Standard Model.
     * @return Br@f$(H\to \tau^+\tau^-)@f$/Br@f$(H\to \tau^+\tau^-)_{\mathrm{SM}}@f$
     */
    virtual const double BrHtautauRatio() const;

    /**
     * @brief The ratio of the @f$\Gamma(H\to cc)@f$ in the current model
     * and in the Standard Model.
     * @return @f$\Gamma(H\to cc)@f$/@f$\Gamma(H\to cc)_{\mathrm{SM}}@f$
     */
    const double GammaHccRatio() const;
    /**
     * @brief The new physics contribution to the ratio of the @f$\Gamma(H\to cc)@f$
     * in the current model and in the Standard Model. (Only terms that are
     * linear in the effective Lagrangian coefficients.)
     * @return @f$\delta \Gamma(H\to cc)@f$/@f$\Gamma(H\to cc)_{\mathrm{SM}}@f$
     */
    const double deltaGammaHccRatio1() const;
    /**
     * @brief The new physics contribution to the ratio of the @f$\Gamma(H\to cc)@f$
     * in the current model and in the Standard Model. (Only terms that are
     * quadratic in the effective Lagrangian coefficients.)
     * @return @f$\delta \Gamma(H\to cc)@f$/@f$\Gamma(H\to cc)_{\mathrm{SM}}@f$
     */
    const double deltaGammaHccRatio2() const;
    /**
     * @brief The ratio of the Br@f$(H\to c\bar{c})@f$ in the current model
     * and in the Standard Model.
     * @return Br@f$(H\to c\bar{c})@f$/Br@f$(H\to c\bar{c})_{\mathrm{SM}}@f$
     */
    virtual const double BrHccRatio() const;

    /**
     * @brief The ratio of the @f$\Gamma(H\to bb)@f$ in the current model
     * and in the Standard Model.
     * @return @f$\Gamma(H\to bb)@f$/@f$\Gamma(H\to bb)_{\mathrm{SM}}@f$
     */
    const double GammaHbbRatio() const;
    /**
     * @brief The new physics contribution to the ratio of the @f$\Gamma(H\to bb)@f$
     * in the current model and in the Standard Model. (Only terms that are
     * linear in the effective Lagrangian coefficients.)
     * @return @f$\delta \Gamma(H\to bb)@f$/@f$\Gamma(H\to bb)_{\mathrm{SM}}@f$
     */
    const double deltaGammaHbbRatio1() const;
    /**
     * @brief The new physics contribution to the ratio of the @f$\Gamma(H\to bb)@f$
     * in the current model and in the Standard Model. (Only terms that are
     * quadratic in the effective Lagrangian coefficients.)
     * @return @f$\delta \Gamma(H\to bb)@f$/@f$\Gamma(H\to bb)_{\mathrm{SM}}@f$
     */
    const double deltaGammaHbbRatio2() const;
    /**
     * @brief The ratio of the Br@f$(H\to b\bar{b})@f$ in the current model
     * and in the Standard Model.
     * @return Br@f$(H\to b\bar{b})@f$/Br@f$(H\to b\bar{b})_{\mathrm{SM}}@f$
     */
    virtual const double BrHbbRatio() const;

    /////////////////////// HIGGS TO 4 FERMION DECAYS /////////////////////////

    /**
     * @brief The ratio of the @f$\Gamma(H\to 2L2L')@f$ (@f$L,L'=e,\mu,\tau@f$) in the current model
     * and in the Standard Model.
     * @return @f$\Gamma(H\to 2L2L')@f$/@f$\Gamma(H\to 2L2L')_{\mathrm{SM}}@f$
     */
    const double GammaH2L2LRatio() const;
    /**
     * @brief The new physics contribution to the ratio of the @f$\Gamma(H\to 2L2L')@f$ (@f$L,L'=e,\mu,\tau@f$)
     * in the current model and in the Standard Model. (Only terms that are
     * linear in the effective Lagrangian coefficients.)
     * @return @f$\delta \Gamma(H\to 2L2L')@f$/@f$\Gamma(H\to 2L2L')_{\mathrm{SM}}@f$
     */
    const double deltaGammaH2L2LRatio1() const;
    /**
     * @brief The new physics contribution to the ratio of the @f$\Gamma(H\to 2L2L')@f$ (@f$L,L'=e,\mu,\tau@f$)
     * in the current model and in the Standard Model. (Only terms that are
     * quadratic in the effective Lagrangian coefficients.)
     * @return @f$\delta \Gamma(H\to 2L2L')@f$/@f$\Gamma(H\to 2L2L')_{\mathrm{SM}}@f$
     */
    const double deltaGammaH2L2LRatio2() const;
    /**
     * @brief The ratio of the Br@f$(H\to 2L2L')@f$ (@f$L,L'=e,\mu,\tau@f$) in the current model
     * and in the Standard Model.
     * @return Br@f$(H\to 2L2L')@f$/Br@f$(H\to 2L2L')_{\mathrm{SM}}@f$
     */
    virtual const double BrH2L2LRatio() const;

    /**
     * @brief The ratio of the @f$\Gamma(H\to 2e 2\mu)@f$ in the current model
     * and in the Standard Model.
     * @return @f$\Gamma(H\to 2e 2\mu)@f$/@f$\Gamma(H\to 2e 2\mu)_{\mathrm{SM}}@f$
     */
    const double GammaH2e2muRatio() const;
    /**
     * @brief The new physics contribution to the ratio of the @f$\Gamma(H\to 2e 2\mu)@f$
     * in the current model and in the Standard Model. (Only terms that are
     * linear in the effective Lagrangian coefficients.)
     * @return @f$\delta \Gamma(H\to 2e 2\mu)@f$/@f$\Gamma(H\to 2e 2\mu)_{\mathrm{SM}}@f$
     */
    const double deltaGammaH2e2muRatio1() const;
    /**
     * @brief The new physics contribution to the ratio of the @f$\Gamma(H\to 2e 2\mu)@f$
     * in the current model and in the Standard Model. (Only terms that are
     * quadratic in the effective Lagrangian coefficients.)
     * @return @f$\delta \Gamma(H\to 2e 2\mu)@f$/@f$\Gamma(H\to 2e 2\mu)_{\mathrm{SM}}@f$
     */
    const double deltaGammaH2e2muRatio2() const;
    /**
     * @brief The ratio of the Br@f$(H\to 2e 2\mu)@f$ in the current model
     * and in the Standard Model.
     * @return Br@f$(H\to 2L2L)@f$/Br@f$(H\to 2e 2\mu)_{\mathrm{SM}}@f$
     */
    virtual const double BrH2e2muRatio() const;

    /**
     * @brief The ratio of the @f$\Gamma(H\to 2v2v)@f$ in the current model
     * and in the Standard Model.
     * @return @f$\Gamma(H\to 2v2v)@f$/@f$\Gamma(H\to 2v2v)_{\mathrm{SM}}@f$
     */
    const double GammaH2v2vRatio() const;
    /**
     * @brief The new physics contribution to the ratio of the @f$\Gamma(H\to 2v2v)@f$
     * in the current model and in the Standard Model. (Only terms that are
     * linear in the effective Lagrangian coefficients.)
     * @return @f$\delta \Gamma(H\to 2v2v)@f$/@f$\Gamma(H\to 2v2v)_{\mathrm{SM}}@f$
     */
    const double deltaGammaH2v2vRatio1() const;
    /**
     * @brief The new physics contribution to the ratio of the @f$\Gamma(H\to 2v2v)@f$
     * in the current model and in the Standard Model. (Only terms that are
     * quadratic in the effective Lagrangian coefficients.)
     * @return @f$\delta \Gamma(H\to 2v2v)@f$/@f$\Gamma(H\to 2v2v)_{\mathrm{SM}}@f$
     */
    const double deltaGammaH2v2vRatio2() const;
    /**
     * @brief The ratio of the Br@f$(H\to 2v2v)@f$ in the current model
     * and in the Standard Model.
     * @return Br@f$(H\to 2v2v)@f$/Br@f$(H\to 2v2v)_{\mathrm{SM}}@f$
     */
    virtual const double BrH2v2vRatio() const;

    /**
     * @brief The ratio of the @f$\Gamma(H\to 2L2v)@f$ (@f$L=e,\mu,\tau@f$) in the current model
     * and in the Standard Model.
     * @return @f$\Gamma(H\to 2L2v)@f$/@f$\Gamma(H\to 2L2v)_{\mathrm{SM}}@f$
     */
    const double GammaH2L2vRatio() const;
    /**
     * @brief The new physics contribution to the ratio of the @f$\Gamma(H\to 2L2v)@f$ (@f$L=e,\mu,\tau@f$)
     * in the current model and in the Standard Model. (Only terms that are
     * linear in the effective Lagrangian coefficients.)
     * @return @f$\delta \Gamma(H\to 2L2v)@f$/@f$\Gamma(H\to 2L2v)_{\mathrm{SM}}@f$
     */
    const double deltaGammaH2L2vRatio1() const;
    /**
     * @brief The new physics contribution to the ratio of the @f$\Gamma(H\to 2L2v)@f$ (@f$L=e,\mu,\tau@f$)
     * in the current model and in the Standard Model. (Only terms that are
     * quadratic in the effective Lagrangian coefficients.)
     * @return @f$\delta \Gamma(H\to 2L2v)@f$/@f$\Gamma(H\to 2L2v)_{\mathrm{SM}}@f$
     */
    const double deltaGammaH2L2vRatio2() const;
    /**
     * @brief The ratio of the Br@f$(H\to 2L2v)@f$ (@f$L=e,\mu,\tau@f$) in the current model
     * and in the Standard Model.
     * @return Br@f$(H\to 2L2v)@f$/Br@f$(H\to 2L2v)_{\mathrm{SM}}@f$
     */
    virtual const double BrH2L2vRatio() const;


    /**
     * @brief The ratio of the @f$\Gamma(H\to 2L2v)@f$ (@f$L=e,\mu@f$) in the current model
     * and in the Standard Model.
     * @return @f$\Gamma(H\to 2L2v)@f$/@f$\Gamma(H\to 2L2v)_{\mathrm{SM}}@f$
     */
    const double GammaH2L2v2Ratio() const;
    /**
     * @brief The new physics contribution to the ratio of the @f$\Gamma(H\to 2L2v)@f$ (@f$L=e,\mu@f$)
     * in the current model and in the Standard Model. (Only terms that are
     * linear in the effective Lagrangian coefficients.)
     * @return @f$\delta \Gamma(H\to 2L2v)@f$/@f$\Gamma(H\to 2L2v)_{\mathrm{SM}}@f$
     */
    const double deltaGammaH2L2v2Ratio1() const;
    /**
     * @brief The new physics contribution to the ratio of the @f$\Gamma(H\to 2L2v)@f$ (@f$L=e,\mu@f$)
     * in the current model and in the Standard Model. (Only terms that are
     * quadratic in the effective Lagrangian coefficients.)
     * @return @f$\delta \Gamma(H\to 2L2v)@f$/@f$\Gamma(H\to 2L2v)_{\mathrm{SM}}@f$
     */
    const double deltaGammaH2L2v2Ratio2() const;
    /**
     * @brief The ratio of the Br@f$(H\to 2L2v)@f$ (@f$L=e,\mu@f$) in the current model
     * and in the Standard Model.
     * @return Br@f$(H\to 2L2v)@f$/Br@f$(H\to 2L2v)_{\mathrm{SM}}@f$
     */
    virtual const double BrH2L2v2Ratio() const;


    /**
     * @brief The ratio of the @f$\Gamma(H\to 2e2v)@f$ in the current model
     * and in the Standard Model.
     * @return @f$\Gamma(H\to 2e2v)@f$/@f$\Gamma(H\to 2e2v)_{\mathrm{SM}}@f$
     */
    const double GammaH2e2vRatio() const;
    /**
     * @brief The new physics contribution to the ratio of the @f$\Gamma(H\to 2e2v)@f$
     * in the current model and in the Standard Model. (Only terms that are
     * linear in the effective Lagrangian coefficients.)
     * @return @f$\delta \Gamma(H\to 2e2v)@f$/@f$\Gamma(H\to 2e2v)_{\mathrm{SM}}@f$
     */
    const double deltaGammaH2e2vRatio1() const;
    /**
     * @brief The new physics contribution to the ratio of the @f$\Gamma(H\to 2e2v)@f$
     * in the current model and in the Standard Model. (Only terms that are
     * quadratic in the effective Lagrangian coefficients.)
     * @return @f$\delta \Gamma(H\to 2e2v)@f$/@f$\Gamma(H\to 2e2v)_{\mathrm{SM}}@f$
     */
    const double deltaGammaH2e2vRatio2() const;
    /**
     * @brief The ratio of the Br@f$(H\to 2e2v)@f$ in the current model
     * and in the Standard Model.
     * @return Br@f$(H\to 2e2v)@f$/Br@f$(H\to 2e2v)_{\mathrm{SM}}@f$
     */
    virtual const double BrH2e2vRatio() const;


    /**
     * @brief The ratio of the @f$\Gamma(H\to 2\mu 2v)@f$ in the current model
     * and in the Standard Model.
     * @return @f$\Gamma(H\to 2\mu 2v)@f$/@f$\Gamma(H\to 2\mu 2v)_{\mathrm{SM}}@f$
     */
    const double GammaH2mu2vRatio() const;
    /**
     * @brief The new physics contribution to the ratio of the @f$\Gamma(H\to 2\mu 2v)@f$
     * in the current model and in the Standard Model. (Only terms that are
     * linear in the effective Lagrangian coefficients.)
     * @return @f$\delta \Gamma(H\to 2\mu 2v)@f$/@f$\Gamma(H\to 2\mu 2v)_{\mathrm{SM}}@f$
     */
    const double deltaGammaH2mu2vRatio1() const;
    /**
     * @brief The new physics contribution to the ratio of the @f$\Gamma(H\to 2\mu 2v)@f$
     * in the current model and in the Standard Model. (Only terms that are
     * quadratic in the effective Lagrangian coefficients.)
     * @return @f$\delta \Gamma(H\to 2\mu 2v)@f$/@f$\Gamma(H\to 2\mu 2v)_{\mathrm{SM}}@f$
     */
    const double deltaGammaH2mu2vRatio2() const;
    /**
     * @brief The ratio of the Br@f$(H\to 2\mu 2v)@f$ in the current model
     * and in the Standard Model.
     * @return Br@f$(H\to 2\mu 2v)@f$/Br@f$(H\to 2\mu 2v)_{\mathrm{SM}}@f$
     */
    virtual const double BrH2mu2vRatio() const;


    /**
     * @brief The ratio of the @f$\Gamma(H\to 2u2u)@f$ in the current model
     * and in the Standard Model.
     * @return @f$\Gamma(H\to 2u2u)@f$/@f$\Gamma(H\to 2u2u)_{\mathrm{SM}}@f$
     */
    const double GammaH2u2uRatio() const;
    /**
     * @brief The new physics contribution to the ratio of the @f$\Gamma(H\to 2u2u)@f$
     * in the current model and in the Standard Model. (Only terms that are
     * linear in the effective Lagrangian coefficients.)
     * @return @f$\delta \Gamma(H\to 2u2u)@f$/@f$\Gamma(H\to 2u2u)_{\mathrm{SM}}@f$
     */
    const double deltaGammaH2u2uRatio1() const;
    /**
     * @brief The new physics contribution to the ratio of the @f$\Gamma(H\to 2u2u)@f$
     * in the current model and in the Standard Model. (Only terms that are
     * quadratic in the effective Lagrangian coefficients.)
     * @return @f$\delta \Gamma(H\to 2u2u)@f$/@f$\Gamma(H\to 2u2u)_{\mathrm{SM}}@f$
     */
    const double deltaGammaH2u2uRatio2() const;
    /**
     * @brief The ratio of the Br@f$(H\to 2u2u)@f$ in the current model
     * and in the Standard Model.
     * @return Br@f$(H\to 2u2u)@f$/Br@f$(H\to 2u2u)_{\mathrm{SM}}@f$
     */
    virtual const double BrH2u2uRatio() const;

    /**
     * @brief The ratio of the @f$\Gamma(H\to 2d2d)@f$ in the current model
     * and in the Standard Model.
     * @return @f$\Gamma(H\to 2d2d)@f$/@f$\Gamma(H\to 2d2d)_{\mathrm{SM}}@f$
     */
    const double GammaH2d2dRatio() const;
    /**
     * @brief The new physics contribution to the ratio of the @f$\Gamma(H\to 2d2d)@f$
     * in the current model and in the Standard Model. (Only terms that are
     * linear in the effective Lagrangian coefficients.)
     * @return @f$\delta \Gamma(H\to 2d2d)@f$/@f$\Gamma(H\to 2d2d)_{\mathrm{SM}}@f$
     */
    const double deltaGammaH2d2dRatio1() const;
    /**
     * @brief The new physics contribution to the ratio of the @f$\Gamma(H\to 2d2d)@f$
     * in the current model and in the Standard Model. (Only terms that are
     * quadratic in the effective Lagrangian coefficients.)
     * @return @f$\delta \Gamma(H\to 2d2d)@f$/@f$\Gamma(H\to 2d2d)_{\mathrm{SM}}@f$
     */
    const double deltaGammaH2d2dRatio2() const;
    /**
     * @brief The ratio of the Br@f$(H\to 2d2d)@f$ in the current model
     * and in the Standard Model.
     * @return Br@f$(H\to 2d2d)@f$/Br@f$(H\to 2d2d)_{\mathrm{SM}}@f$
     */
    virtual const double BrH2d2dRatio() const;

    /**
     * @brief The ratio of the @f$\Gamma(H\to 2u2d)@f$ in the current model
     * and in the Standard Model.
     * @return @f$\Gamma(H\to 2u2d)@f$/@f$\Gamma(H\to 2u2d)_{\mathrm{SM}}@f$
     */
    const double GammaH2u2dRatio() const;
    /**
     * @brief The new physics contribution to the ratio of the @f$\Gamma(H\to 2u2d)@f$
     * in the current model and in the Standard Model. (Only terms that are
     * linear in the effective Lagrangian coefficients.)
     * @return @f$\delta \Gamma(H\to 2u2d)@f$/@f$\Gamma(H\to 2u2d)_{\mathrm{SM}}@f$
     */
    const double deltaGammaH2u2dRatio1() const;
    /**
     * @brief The new physics contribution to the ratio of the @f$\Gamma(H\to 2u2d)@f$
     * in the current model and in the Standard Model. (Only terms that are
     * quadratic in the effective Lagrangian coefficients.)
     * @return @f$\delta \Gamma(H\to 2u2d)@f$/@f$\Gamma(H\to 2u2d)_{\mathrm{SM}}@f$
     */
    const double deltaGammaH2u2dRatio2() const;
    /**
     * @brief The ratio of the Br@f$(H\to 2u2d)@f$ in the current model
     * and in the Standard Model.
     * @return Br@f$(H\to 2u2d)@f$/Br@f$(H\to 2u2d)_{\mathrm{SM}}@f$
     */
    virtual const double BrH2u2dRatio() const;

    /**
     * @brief The ratio of the @f$\Gamma(H\to 2L2u)@f$ (@f$L=e,\mu,\tau@f$) in the current model
     * and in the Standard Model.
     * @return @f$\Gamma(H\to 2L2u)@f$/@f$\Gamma(H\to 2L2u)_{\mathrm{SM}}@f$
     */
    const double GammaH2L2uRatio() const;
    /**
     * @brief The new physics contribution to the ratio of the @f$\Gamma(H\to 2L2u)@f$ (@f$L=e,\mu,\tau@f$)
     * in the current model and in the Standard Model. (Only terms that are
     * linear in the effective Lagrangian coefficients.)
     * @return @f$\delta \Gamma(H\to 2L2u)@f$/@f$\Gamma(H\to 2L2u)_{\mathrm{SM}}@f$
     */
    const double deltaGammaH2L2uRatio1() const;
    /**
     * @brief The new physics contribution to the ratio of the @f$\Gamma(H\to 2L2u)@f$ (@f$L=e,\mu,\tau@f$)
     * in the current model and in the Standard Model. (Only terms that are
     * quadratic in the effective Lagrangian coefficients.)
     * @return @f$\delta \Gamma(H\to 2L2u)@f$/@f$\Gamma(H\to 2L2u)_{\mathrm{SM}}@f$
     */
    const double deltaGammaH2L2uRatio2() const;
    /**
     * @brief The ratio of the Br@f$(H\to 2L2u)@f$ (@f$L=e,\mu,\tau@f$) in the current model
     * and in the Standard Model.
     * @return Br@f$(H\to 2L2u)@f$/Br@f$(H\to 2L2u)_{\mathrm{SM}}@f$
     */
    virtual const double BrH2L2uRatio() const;

    /**
     * @brief The ratio of the @f$\Gamma(H\to 2L2d)@f$ (@f$L=e,\mu,\tau@f$) in the current model
     * and in the Standard Model.
     * @return @f$\Gamma(H\to 2L2d)@f$/@f$\Gamma(H\to 2L2d)_{\mathrm{SM}}@f$
     */
    const double GammaH2L2dRatio() const;
    /**
     * @brief The new physics contribution to the ratio of the @f$\Gamma(H\to 2L2d)@f$ (@f$L=e,\mu,\tau@f$)
     * in the current model and in the Standard Model. (Only terms that are
     * linear in the effective Lagrangian coefficients.)
     * @return @f$\delta \Gamma(H\to 2L2d)@f$/@f$\Gamma(H\to 2L2d)_{\mathrm{SM}}@f$
     */
    const double deltaGammaH2L2dRatio1() const;
    /**
     * @brief The new physics contribution to the ratio of the @f$\Gamma(H\to 2L2d)@f$ (@f$L=e,\mu,\tau@f$)
     * in the current model and in the Standard Model. (Only terms that are
     * quadratic in the effective Lagrangian coefficients.)
     * @return @f$\delta \Gamma(H\to 2L2d)@f$/@f$\Gamma(H\to 2L2d)_{\mathrm{SM}}@f$
     */
    const double deltaGammaH2L2dRatio2() const;
    /**
     * @brief The ratio of the Br@f$(H\to 2L2d)@f$ (@f$L=e,\mu,\tau@f$) in the current model
     * and in the Standard Model.
     * @return Br@f$(H\to 2L2d)@f$/Br@f$(H\to 2L2d)_{\mathrm{SM}}@f$
     */
    virtual const double BrH2L2dRatio() const;

    /**
     * @brief The ratio of the @f$\Gamma(H\to 2v2u)@f$ in the current model
     * and in the Standard Model.
     * @return @f$\Gamma(H\to 2v2u)@f$/@f$\Gamma(H\to 2v2u)_{\mathrm{SM}}@f$
     */
    const double GammaH2v2uRatio() const;
    /**
     * @brief The new physics contribution to the ratio of the @f$\Gamma(H\to 2v2u)@f$
     * in the current model and in the Standard Model. (Only terms that are
     * linear in the effective Lagrangian coefficients.)
     * @return @f$\delta \Gamma(H\to 2v2u)@f$/@f$\Gamma(H\to 2v2u)_{\mathrm{SM}}@f$
     */
    const double deltaGammaH2v2uRatio1() const;
    /**
     * @brief The new physics contribution to the ratio of the @f$\Gamma(H\to 2v2u)@f$
     * in the current model and in the Standard Model. (Only terms that are
     * quadratic in the effective Lagrangian coefficients.)
     * @return @f$\delta \Gamma(H\to 2v2u)@f$/@f$\Gamma(H\to 2v2u)_{\mathrm{SM}}@f$
     */
    const double deltaGammaH2v2uRatio2() const;
    /**
     * @brief The ratio of the Br@f$(H\to 2v2u)@f$ in the current model
     * and in the Standard Model.
     * @return Br@f$(H\to 2v2u)@f$/Br@f$(H\to 2v2u)_{\mathrm{SM}}@f$
     */
    virtual const double BrH2v2uRatio() const;

    /**
     * @brief The ratio of the @f$\Gamma(H\to 2v2d)@f$ in the current model
     * and in the Standard Model.
     * @return @f$\Gamma(H\to 2v2d)@f$/@f$\Gamma(H\to 2v2d)_{\mathrm{SM}}@f$
     */
    const double GammaH2v2dRatio() const;
    /**
     * @brief The new physics contribution to the ratio of the @f$\Gamma(H\to 2v2d)@f$
     * in the current model and in the Standard Model. (Only terms that are
     * linear in the effective Lagrangian coefficients.)
     * @return @f$\delta \Gamma(H\to 2v2d)@f$/@f$\Gamma(H\to 2v2d)_{\mathrm{SM}}@f$
     */
    const double deltaGammaH2v2dRatio1() const;
    /**
     * @brief The new physics contribution to the ratio of the @f$\Gamma(H\to 2v2d)@f$
     * in the current model and in the Standard Model. (Only terms that are
     * quadratic in the effective Lagrangian coefficients.)
     * @return @f$\delta \Gamma(H\to 2v2d)@f$/@f$\Gamma(H\to 2v2d)_{\mathrm{SM}}@f$
     */
    const double deltaGammaH2v2dRatio2() const;
    /**
     * @brief The ratio of the Br@f$(H\to 2v2d)@f$ in the current model
     * and in the Standard Model.
     * @return Br@f$(H\to 2v2d)@f$/Br@f$(H\to 2v2d)_{\mathrm{SM}}@f$
     */
    virtual const double BrH2v2dRatio() const;

    /**
     * @brief The ratio of the @f$\Gamma(H\to 4L)@f$ (@f$L=e,\mu,\tau@f$) in the current model
     * and in the Standard Model.
     * @return @f$\Gamma(H\to 4L)@f$/@f$\Gamma(H\to 4L)_{\mathrm{SM}}@f$
     */
    const double GammaH4LRatio() const;
    /**
     * @brief The new physics contribution to the ratio of the @f$\Gamma(H\to 4L)@f$ (@f$L=e,\mu,\tau@f$)
     * in the current model and in the Standard Model. (Only terms that are
     * linear in the effective Lagrangian coefficients.)
     * @return @f$\delta \Gamma(H\to 4L)@f$/@f$\Gamma(H\to 4L)_{\mathrm{SM}}@f$
     */
    const double deltaGammaH4LRatio1() const;
    /**
     * @brief The new physics contribution to the ratio of the @f$\Gamma(H\to 4L)@f$ (@f$L=e,\mu,\tau@f$)
     * in the current model and in the Standard Model. (Only terms that are
     * quadratic in the effective Lagrangian coefficients.)
     * @return @f$\delta \Gamma(H\to 4L)@f$/@f$\Gamma(H\to 4L)_{\mathrm{SM}}@f$
     */
    const double deltaGammaH4LRatio2() const;
    /**
     * @brief The ratio of the Br@f$(H\to 4L)@f$ (@f$L=e,\mu,\tau@f$) in the current model
     * and in the Standard Model.
     * @return Br@f$(H\to 4L)@f$/Br@f$(H\to 4L)_{\mathrm{SM}}@f$
     */
    virtual const double BrH4LRatio() const;

    /**
     * @brief The ratio of the @f$\Gamma(H\to 4L)@f$ (@f$L=e,\mu@f$) in the current model
     * and in the Standard Model.
     * @return @f$\Gamma(H\to 4L)@f$/@f$\Gamma(H\to 4L)_{\mathrm{SM}}@f$
     */
    const double GammaH4L2Ratio() const;
    /**
     * @brief The new physics contribution to the ratio of the @f$\Gamma(H\to 4L)@f$ (@f$L=e,\mu@f$)
     * in the current model and in the Standard Model. (Only terms that are
     * linear in the effective Lagrangian coefficients.)
     * @return @f$\delta \Gamma(H\to 4L)@f$/@f$\Gamma(H\to 4L)_{\mathrm{SM}}@f$
     */
    const double deltaGammaH4L2Ratio1() const;
    /**
     * @brief The new physics contribution to the ratio of the @f$\Gamma(H\to 4L)@f$ (@f$L=e,\mu@f$)
     * in the current model and in the Standard Model. (Only terms that are
     * quadratic in the effective Lagrangian coefficients.)
     * @return @f$\delta \Gamma(H\to 4L)@f$/@f$\Gamma(H\to 4L)_{\mathrm{SM}}@f$
     */
    const double deltaGammaH4L2Ratio2() const;
    /**
     * @brief The ratio of the Br@f$(H\to 4L)@f$ (@f$L=e,\mu@f$) in the current model
     * and in the Standard Model.
     * @return Br@f$(H\to 4L)@f$/Br@f$(H\to 4L)_{\mathrm{SM}}@f$
     */
    virtual const double BrH4L2Ratio() const;


    /**
     * @brief The ratio of the @f$\Gamma(H\to 4e)@f$ in the current model
     * and in the Standard Model.
     * @return @f$\Gamma(H\to 4e)@f$/@f$\Gamma(H\to 4e)_{\mathrm{SM}}@f$
     */
    const double GammaH4eRatio() const;
    /**
     * @brief The new physics contribution to the ratio of the @f$\Gamma(H\to 4e)@f$
     * in the current model and in the Standard Model. (Only terms that are
     * linear in the effective Lagrangian coefficients.)
     * @return @f$\delta \Gamma(H\to 4e)@f$/@f$\Gamma(H\to 4e)_{\mathrm{SM}}@f$
     */
    const double deltaGammaH4eRatio1() const;
    /**
     * @brief The new physics contribution to the ratio of the @f$\Gamma(H\to 4e)@f$
     * in the current model and in the Standard Model. (Only terms that are
     * quadratic in the effective Lagrangian coefficients.)
     * @return @f$\delta \Gamma(H\to 4e)@f$/@f$\Gamma(H\to 4e)_{\mathrm{SM}}@f$
     */
    const double deltaGammaH4eRatio2() const;
    /**
     * @brief The ratio of the Br@f$(H\to 4e)@f$ in the current model
     * and in the Standard Model.
     * @return Br@f$(H\to 4e)@f$/Br@f$(H\to 4e)_{\mathrm{SM}}@f$
     */
    virtual const double BrH4eRatio() const;


    /**
     * @brief The ratio of the @f$\Gamma(H\to 4\mu)@f$ in the current model
     * and in the Standard Model.
     * @return @f$\Gamma(H\to 4\mu)@f$/@f$\Gamma(H\to 4\mu)_{\mathrm{SM}}@f$
     */
    const double GammaH4muRatio() const;
    /**
     * @brief The new physics contribution to the ratio of the @f$\Gamma(H\to 4\mu)@f$
     * in the current model and in the Standard Model. (Only terms that are
     * linear in the effective Lagrangian coefficients.)
     * @return @f$\delta \Gamma(H\to 4\mu)@f$/@f$\Gamma(H\to 4\mu)_{\mathrm{SM}}@f$
     */
    const double deltaGammaH4muRatio1() const;
    /**
     * @brief The new physics contribution to the ratio of the @f$\Gamma(H\to 4\mu)@f$
     * in the current model and in the Standard Model. (Only terms that are
     * quadratic in the effective Lagrangian coefficients.)
     * @return @f$\delta \Gamma(H\to 4\mu)@f$/@f$\Gamma(H\to 4\mu)_{\mathrm{SM}}@f$
     */
    const double deltaGammaH4muRatio2() const;
    /**
     * @brief The ratio of the Br@f$(H\to 4\mu)@f$ in the current model
     * and in the Standard Model.
     * @return Br@f$(H\to 4\mu)@f$/Br@f$(H\to 4\mu)_{\mathrm{SM}}@f$
     */
    virtual const double BrH4muRatio() const;


    /**
     * @brief The ratio of the @f$\Gamma(H\to 4v)@f$ in the current model
     * and in the Standard Model.
     * @return @f$\Gamma(H\to 4v)@f$/@f$\Gamma(H\to 4v)_{\mathrm{SM}}@f$
     */
    const double GammaH4vRatio() const;
    /**
     * @brief The new physics contribution to the ratio of the @f$\Gamma(H\to 4v)@f$
     * in the current model and in the Standard Model. (Only terms that are
     * linear in the effective Lagrangian coefficients.)
     * @return @f$\delta \Gamma(H\to 4v)@f$/@f$\Gamma(H\to 4v)_{\mathrm{SM}}@f$
     */
    const double deltaGammaH4vRatio1() const;
    /**
     * @brief The new physics contribution to the ratio of the @f$\Gamma(H\to 4v)@f$
     * in the current model and in the Standard Model. (Only terms that are
     * quadratic in the effective Lagrangian coefficients.)
     * @return @f$\delta \Gamma(H\to 4v)@f$/@f$\Gamma(H\to 4v)_{\mathrm{SM}}@f$
     */
    const double deltaGammaH4vRatio2() const;
    /**
     * @brief The ratio of the Br@f$(H\to 4v)@f$ in the current model
     * and in the Standard Model.
     * @return Br@f$(H\to 4v)@f$/Br@f$(H\to 4v)_{\mathrm{SM}}@f$
     */
    virtual const double BrH4vRatio() const;

    /**
     * @brief The ratio of the @f$\Gamma(H\to 4u)@f$ in the current model
     * and in the Standard Model.
     * @return @f$\Gamma(H\to 4u)@f$/@f$\Gamma(H\to 4u)_{\mathrm{SM}}@f$
     */
    const double GammaH4uRatio() const;
    /**
     * @brief The new physics contribution to the ratio of the @f$\Gamma(H\to 4u)@f$
     * in the current model and in the Standard Model. (Only terms that are
     * linear in the effective Lagrangian coefficients.)
     * @return @f$\delta \Gamma(H\to 4u)@f$/@f$\Gamma(H\to 4u)_{\mathrm{SM}}@f$
     */
    const double deltaGammaH4uRatio1() const;
    /**
     * @brief The new physics contribution to the ratio of the @f$\Gamma(H\to 4u)@f$
     * in the current model and in the Standard Model. (Only terms that are
     * quadratic in the effective Lagrangian coefficients.)
     * @return @f$\delta \Gamma(H\to 4u)@f$/@f$\Gamma(H\to 4u)_{\mathrm{SM}}@f$
     */
    const double deltaGammaH4uRatio2() const;
    /**
     * @brief The ratio of the Br@f$(H\to 4u)@f$ in the current model
     * and in the Standard Model.
     * @return Br@f$(H\to 4u)@f$/Br@f$(H\to 4u)_{\mathrm{SM}}@f$
     */
    virtual const double BrH4uRatio() const;

    /**
     * @brief The ratio of the @f$\Gamma(H\to 4d)@f$ in the current model
     * and in the Standard Model.
     * @return @f$\Gamma(H\to 4d)@f$/@f$\Gamma(H\to 4d)_{\mathrm{SM}}@f$
     */
    const double GammaH4dRatio() const;
    /**
     * @brief The new physics contribution to the ratio of the @f$\Gamma(H\to 4d)@f$
     * in the current model and in the Standard Model. (Only terms that are
     * linear in the effective Lagrangian coefficients.)
     * @return @f$\delta \Gamma(H\to 4d)@f$/@f$\Gamma(H\to 4d)_{\mathrm{SM}}@f$
     */
    const double deltaGammaH4dRatio1() const;
    /**
     * @brief The new physics contribution to the ratio of the @f$\Gamma(H\to 4d)@f$
     * in the current model and in the Standard Model. (Only terms that are
     * quadratic in the effective Lagrangian coefficients.)
     * @return @f$\delta \Gamma(H\to 4d)@f$/@f$\Gamma(H\to 4d)_{\mathrm{SM}}@f$
     */
    const double deltaGammaH4dRatio2() const;
    /**
     * @brief The ratio of the Br@f$(H\to 4d)@f$ in the current model
     * and in the Standard Model.
     * @return Br@f$(H\to 4d)@f$/Br@f$(H\to 4d)_{\mathrm{SM}}@f$
     */
    virtual const double BrH4dRatio() const;

    /**
     * @brief The ratio of the @f$\Gamma(H\to LvvL)@f$ (@f$L=e,\mu,\tau@f$) in the current model
     * and in the Standard Model.
     * @return @f$\Gamma(H\to LvvL)@f$/@f$\Gamma(H\to LvvL)_{\mathrm{SM}}@f$
     */
    const double GammaHLvvLRatio() const;
    /**
     * @brief The new physics contribution to the ratio of the @f$\Gamma(H\to LvvL)@f$ (@f$L=e,\mu,\tau@f$)
     * in the current model and in the Standard Model. (Only terms that are
     * linear in the effective Lagrangian coefficients.)
     * @return @f$\delta \Gamma(H\to LvvL)@f$/@f$\Gamma(H\to LvvL)_{\mathrm{SM}}@f$
     */
    const double deltaGammaHLvvLRatio1() const;
    /**
     * @brief The new physics contribution to the ratio of the @f$\Gamma(H\to LvvL)@f$ (@f$L=e,\mu,\tau@f$)
     * in the current model and in the Standard Model. (Only terms that are
     * quadratic in the effective Lagrangian coefficients.)
     * @return @f$\delta \Gamma(H\to LvvL)@f$/@f$\Gamma(H\to LvvL)_{\mathrm{SM}}@f$
     */
    const double deltaGammaHLvvLRatio2() const;
    /**
     * @brief The ratio of the Br@f$(H\to LvvL)@f$ (@f$L=e,\mu,\tau@f$) in the current model
     * and in the Standard Model.
     * @return Br@f$(H\to LvvL)@f$/Br@f$(H\to LvvL)_{\mathrm{SM}}@f$
     */
    virtual const double BrHLvvLRatio() const;

    /**
     * @brief The ratio of the @f$\Gamma(H\to e\nu \mu\nu)@f$ in the current model
     * and in the Standard Model.
     * @return @f$\Gamma(H\to e\nu \mu\nu)@f$/@f$\Gamma(H\to e\nu \mu\nu)_{\mathrm{SM}}@f$
     */
    const double GammaHevmuvRatio() const;
    /**
     * @brief The new physics contribution to the ratio of the @f$\Gamma(H\to e\nu \mu\nu)@f$
     * in the current model and in the Standard Model. (Only terms that are
     * linear in the effective Lagrangian coefficients.)
     * @return @f$\delta \Gamma(H\to e\nu \mu\nu)@f$/@f$\Gamma(H\to e\nu \mu\nu)_{\mathrm{SM}}@f$
     */
    const double deltaGammaHevmuvRatio1() const;
    /**
     * @brief The new physics contribution to the ratio of the @f$\Gamma(H\to e\nu \mu\nu)@f$
     * in the current model and in the Standard Model. (Only terms that are
     * quadratic in the effective Lagrangian coefficients.)
     * @return @f$\delta \Gamma(H\to e\nu \mu\nu)@f$/@f$\Gamma(H\to e\nu \mu\nu)_{\mathrm{SM}}@f$
     */
    const double deltaGammaHevmuvRatio2() const;
    /**
     * @brief The ratio of the Br@f$(H\to e\nu \mu\nu)@f$ in the current model
     * and in the Standard Model.
     * @return Br@f$(H\to e\nu \mu\nu)@f$/Br@f$(H\to e\nu \mu\nu)_{\mathrm{SM}}@f$
     */
    virtual const double BrHevmuvRatio() const;


    /**
     * @brief The ratio of the @f$\Gamma(H\to uddu)@f$ in the current model
     * and in the Standard Model.
     * @return @f$\Gamma(H\to uddu)@f$/@f$\Gamma(H\to uddu)_{\mathrm{SM}}@f$
     */
    const double GammaHudduRatio() const;
    /**
     * @brief The new physics contribution to the ratio of the @f$\Gamma(H\to uddu)@f$
     * in the current model and in the Standard Model. (Only terms that are
     * linear in the effective Lagrangian coefficients.)
     * @return @f$\delta \Gamma(H\to uddu)@f$/@f$\Gamma(H\to uddu)_{\mathrm{SM}}@f$
     */
    const double deltaGammaHudduRatio1() const;
    /**
     * @brief The new physics contribution to the ratio of the @f$\Gamma(H\to uddu)@f$
     * in the current model and in the Standard Model. (Only terms that are
     * quadratic in the effective Lagrangian coefficients.)
     * @return @f$\delta \Gamma(H\to uddu)@f$/@f$\Gamma(H\to uddu)_{\mathrm{SM}}@f$
     */
    const double deltaGammaHudduRatio2() const;
    /**
     * @brief The ratio of the Br@f$(H\to uddu)@f$ in the current model
     * and in the Standard Model.
     * @return Br@f$(H\to uddu)@f$/Br@f$(H\to uddu)_{\mathrm{SM}}@f$
     */
    virtual const double BrHudduRatio() const;

    /**
     * @brief The ratio of the @f$\Gamma(H\to Lvud)@f$ (@f$L=e,\mu,\tau@f$) in the current model
     * and in the Standard Model.
     * @return @f$\Gamma(H\to Lvud)@f$/@f$\Gamma(H\to Lvud)_{\mathrm{SM}}@f$
     */
    const double GammaHLvudRatio() const;
    /**
     * @brief The new physics contribution to the ratio of the @f$\Gamma(H\to Lvud)@f$ (@f$L=e,\mu,\tau@f$)
     * in the current model and in the Standard Model. (Only terms that are
     * linear in the effective Lagrangian coefficients.)
     * @return @f$\delta \Gamma(H\to Lvud)@f$/@f$\Gamma(H\to Lvud)_{\mathrm{SM}}@f$
     */
    const double deltaGammaHLvudRatio1() const;
    /**
     * @brief The new physics contribution to the ratio of the @f$\Gamma(H\to Lvud)@f$ (@f$L=e,\mu,\tau@f$)
     * in the current model and in the Standard Model. (Only terms that are
     * quadratic in the effective Lagrangian coefficients.)
     * @return @f$\delta \Gamma(H\to Lvud)@f$/@f$\Gamma(H\to Lvud)_{\mathrm{SM}}@f$
     */
    const double deltaGammaHLvudRatio2() const;
    /**
     * @brief The ratio of the Br@f$(H\to Lvud)@f$ (@f$L=e,\mu,\tau@f$) in the current model
     * and in the Standard Model.
     * @return Br@f$(H\to Lvud)@f$/Br@f$(H\to Lvud)_{\mathrm{SM}}@f$
     */
    virtual const double BrHLvudRatio() const;

    /**
     * @brief The ratio of the @f$\Gamma(H\to 2ud)@f$ in the current model
     * and in the Standard Model.
     * @return @f$\Gamma(H\to 2ud)@f$/@f$\Gamma(H\to 2ud)_{\mathrm{SM}}@f$
     */
    const double GammaH2udRatio() const;
    /**
     * @brief The new physics contribution to the ratio of the @f$\Gamma(H\to 2ud)@f$
     * in the current model and in the Standard Model. (Only terms that are
     * linear in the effective Lagrangian coefficients.)
     * @return @f$\delta \Gamma(H\to 2ud)@f$/@f$\Gamma(H\to 2ud)_{\mathrm{SM}}@f$
     */
    const double deltaGammaH2udRatio1() const;
    /**
     * @brief The new physics contribution to the ratio of the @f$\Gamma(H\to 2ud)@f$
     * in the current model and in the Standard Model. (Only terms that are
     * quadratic in the effective Lagrangian coefficients.)
     * @return @f$\delta \Gamma(H\to 2ud)@f$/@f$\Gamma(H\to 2ud)_{\mathrm{SM}}@f$
     */
    const double deltaGammaH2udRatio2() const;
    /**
     * @brief The ratio of the Br@f$(H\to 2ud)@f$ in the current model
     * and in the Standard Model.
     * @return Br@f$(H\to 2ud)@f$/Br@f$(H\to 2ud)_{\mathrm{SM}}@f$
     */
    virtual const double BrH2udRatio() const;

    /**
     * @brief The ratio of the @f$\Gamma(H\to 2Lv)@f$ (@f$L=e,\mu,\tau@f$) in the current model
     * and in the Standard Model.
     * @return @f$\Gamma(H\to 2Lv)@f$/@f$\Gamma(H\to 2Lv)_{\mathrm{SM}}@f$
     */
    const double GammaH2LvRatio() const;
    /**
     * @brief The new physics contribution to the ratio of the @f$\Gamma(H\to 2Lv)@f$ (@f$L=e,\mu,\tau@f$)
     * in the current model and in the Standard Model. (Only terms that are
     * linear in the effective Lagrangian coefficients.)
     * @return @f$\delta \Gamma(H\to 2Lv)@f$/@f$\Gamma(H\to 2Lv)_{\mathrm{SM}}@f$
     */
    const double deltaGammaH2LvRatio1() const;
    /**
     * @brief The new physics contribution to the ratio of the @f$\Gamma(H\to 2Lv)@f$ (@f$L=e,\mu,\tau@f$)
     * in the current model and in the Standard Model. (Only terms that are
     * quadratic in the effective Lagrangian coefficients.)
     * @return @f$\delta \Gamma(H\to 2Lv)@f$/@f$\Gamma(H\to 2Lv)_{\mathrm{SM}}@f$
     */
    const double deltaGammaH2LvRatio2() const;
    /**
     * @brief The ratio of the Br@f$(H\to 2Lv)@f$ (@f$L=e,\mu,\tau@f$) in the current model
     * and in the Standard Model.
     * @return Br@f$(H\to 2Lv)@f$/Br@f$(H\to 2Lv)_{\mathrm{SM}}@f$
     */
    virtual const double BrH2LvRatio() const;

    /**
     * @brief The ratio of the @f$\Gamma(H\to 2Lv)@f$ (@f$L=e,\mu@f$) in the current model
     * and in the Standard Model.
     * @return @f$\Gamma(H\to 2Lv)@f$/@f$\Gamma(H\to 2Lv)_{\mathrm{SM}}@f$
     */
    const double GammaH2Lv2Ratio() const;
    /**
     * @brief The new physics contribution to the ratio of the @f$\Gamma(H\to 2Lv)@f$ (@f$L=e,\mu@f$)
     * in the current model and in the Standard Model. (Only terms that are
     * linear in the effective Lagrangian coefficients.)
     * @return @f$\delta \Gamma(H\to 2Lv)@f$/@f$\Gamma(H\to 2Lv)_{\mathrm{SM}}@f$
     */
    const double deltaGammaH2Lv2Ratio1() const;
    /**
     * @brief The new physics contribution to the ratio of the @f$\Gamma(H\to 2Lv)@f$ (@f$L=e,\mu@f$)
     * in the current model and in the Standard Model. (Only terms that are
     * quadratic in the effective Lagrangian coefficients.)
     * @return @f$\delta \Gamma(H\to 2Lv)@f$/@f$\Gamma(H\to 2Lv)_{\mathrm{SM}}@f$
     */
    const double deltaGammaH2Lv2Ratio2() const;
    /**
     * @brief The ratio of the Br@f$(H\to 2Lv)@f$ (@f$L=e,\mu@f$) in the current model
     * and in the Standard Model.
     * @return Br@f$(H\to 2Lv)@f$/Br@f$(H\to 2Lv)_{\mathrm{SM}}@f$
     */
    virtual const double BrH2Lv2Ratio() const;


    /**
     * @brief The ratio of the @f$\Gamma(H\to 2ev)@f$ in the current model
     * and in the Standard Model.
     * @return @f$\Gamma(H\to 2ev)@f$/@f$\Gamma(H\to 2ev)_{\mathrm{SM}}@f$
     */
    const double GammaH2evRatio() const;
    /**
     * @brief The new physics contribution to the ratio of the @f$\Gamma(H\to 2ev)@f$
     * in the current model and in the Standard Model. (Only terms that are
     * linear in the effective Lagrangian coefficients.)
     * @return @f$\delta \Gamma(H\to 2ev)@f$/@f$\Gamma(H\to 2ev)_{\mathrm{SM}}@f$
     */
    const double deltaGammaH2evRatio1() const;
    /**
     * @brief The new physics contribution to the ratio of the @f$\Gamma(H\to 2ev)@f$
     * in the current model and in the Standard Model. (Only terms that are
     * quadratic in the effective Lagrangian coefficients.)
     * @return @f$\delta \Gamma(H\to 2ev)@f$/@f$\Gamma(H\to 2ev)_{\mathrm{SM}}@f$
     */
    const double deltaGammaH2evRatio2() const;
    /**
     * @brief The ratio of the Br@f$(H\to 2ev)@f$ in the current model
     * and in the Standard Model.
     * @return Br@f$(H\to 2ev)@f$/Br@f$(H\to 2ev)_{\mathrm{SM}}@f$
     */
    virtual const double BrH2evRatio() const;


    /**
     * @brief The ratio of the @f$\Gamma(H\to 2\mu v)@f$ in the current model
     * and in the Standard Model.
     * @return @f$\Gamma(H\to 2\mu v)@f$/@f$\Gamma(H\to 2\mu v)_{\mathrm{SM}}@f$
     */
    const double GammaH2muvRatio() const;
    /**
     * @brief The new physics contribution to the ratio of the @f$\Gamma(H\to 2\mu v)@f$
     * in the current model and in the Standard Model. (Only terms that are
     * linear in the effective Lagrangian coefficients.)
     * @return @f$\delta \Gamma(H\to 2\mu v)@f$/@f$\Gamma(H\to 2\mu v)_{\mathrm{SM}}@f$
     */
    const double deltaGammaH2muvRatio1() const;
    /**
     * @brief The new physics contribution to the ratio of the @f$\Gamma(H\to 2\mu v)@f$
     * in the current model and in the Standard Model. (Only terms that are
     * quadratic in the effective Lagrangian coefficients.)
     * @return @f$\delta \Gamma(H\to 2\mu v)@f$/@f$\Gamma(H\to 2\mu v)_{\mathrm{SM}}@f$
     */
    const double deltaGammaH2muvRatio2() const;
    /**
     * @brief The ratio of the Br@f$(H\to 2ev)@f$ in the current model
     * and in the Standard Model.
     * @return Br@f$(H\to 2\mu v)@f$/Br@f$(H\to 2\mu v)_{\mathrm{SM}}@f$
     */
    virtual const double BrH2muvRatio() const;


    /**
     * @brief The ratio of the @f$\Gamma(H\to 4f)@f$ in the current model
     * and in the Standard Model.
     * @return @f$\Gamma(H\to 4f)@f$/@f$\Gamma(H\to 4f)_{\mathrm{SM}}@f$
     */
    const double GammaH4fRatio() const;
    /**
     * @brief The new physics contribution to the ratio of the @f$\Gamma(H\to 4f)@f$
     * in the current model and in the Standard Model. (Only terms that are
     * linear in the effective Lagrangian coefficients.)
     * @return @f$\delta \Gamma(H\to 4f)@f$/@f$\Gamma(H\to 4f)_{\mathrm{SM}}@f$
     */
    const double deltaGammaH4fRatio1() const;
    /**
     * @brief The new physics contribution to the ratio of the @f$\Gamma(H\to 4f)@f$
     * in the current model and in the Standard Model. (Only terms that are
     * quadratic in the effective Lagrangian coefficients.)
     * @return @f$\delta \Gamma(H\to 4f)@f$/@f$\Gamma(H\to 4f)_{\mathrm{SM}}@f$
     */
    const double deltaGammaH4fRatio2() const;
    /**
     * @brief The ratio of the Br@f$(H\to 4f)@f$ in the current model
     * and in the Standard Model.
     * @return Br@f$(H\to 4f)@f$/Br@f$(H\to 4f)_{\mathrm{SM}}@f$
     */
    virtual const double BrH4fRatio() const;

    // DECAYS INVOLVING ONLY ELECTRONS, MUONS OR NEUTRINOS IN THE FINAL STATES 

    /**
     * @brief The ratio of the @f$\Gamma(H\to 4l)@f$ (@f$l=e,\mu@f$) in the current model
     * and in the Standard Model.
     * @return @f$\Gamma(H\to 4l)@f$/@f$\Gamma(H\to 4l)_{\mathrm{SM}}@f$
     */
    const double GammaH4lRatio() const;
    /**
     * @brief The new physics contribution to the ratio of the @f$\Gamma(H\to 4l)@f$ (@f$l=e,\mu@f$)
     * in the current model and in the Standard Model. (Only terms that are
     * linear in the effective Lagrangian coefficients.)
     * @return @f$\delta \Gamma(H\to 4l)@f$/@f$\Gamma(H\to 4l)_{\mathrm{SM}}@f$
     */
    const double deltaGammaH4lRatio1() const;
    /**
     * @brief The new physics contribution to the ratio of the @f$\Gamma(H\to 4l)@f$ (@f$l=e,\mu@f$)
     * in the current model and in the Standard Model. (Only terms that are
     * quadratic in the effective Lagrangian coefficients.)
     * @return @f$\delta \Gamma(H\to 4l)@f$/@f$\Gamma(H\to 4l)_{\mathrm{SM}}@f$
     */
    const double deltaGammaH4lRatio2() const;
    /**
     * @brief The ratio of the Br@f$(H\to 4l)@f$ (@f$l=e,\mu@f$) in the current model
     * and in the Standard Model.
     * @return Br@f$(H\to 4l)@f$/Br@f$(H\to 4l)_{\mathrm{SM}}@f$
     */
    virtual const double BrH4lRatio() const;

    /**
     * @brief The ratio of the @f$\Gamma(H\to 2l2v)@f$ (@f$l=e,\mu@f$) in the current model
     * and in the Standard Model.
     * @return @f$\Gamma(H\to 2l2v)@f$/@f$\Gamma(H\to 2l2v)_{\mathrm{SM}}@f$
     */
    const double GammaH2l2vRatio() const;
    /**
     * @brief The new physics contribution to the ratio of the @f$\Gamma(H\to 2l2v)@f$ (@f$l=e,\mu@f$)
     * in the current model and in the Standard Model. (Only terms that are
     * linear in the effective Lagrangian coefficients.)
     * @return @f$\delta \Gamma(H\to 2l2v)@f$/@f$\Gamma(H\to 2l2v)_{\mathrm{SM}}@f$
     */
    const double deltaGammaH2l2vRatio1() const;
    /**
     * @brief The new physics contribution to the ratio of the @f$\Gamma(H\to 2l2v)@f$ (@f$l=e,\mu@f$)
     * in the current model and in the Standard Model. (Only terms that are
     * quadratic in the effective Lagrangian coefficients.)
     * @return @f$\delta \Gamma(H\to 2l2v)@f$/@f$\Gamma(H\to 2l2v)_{\mathrm{SM}}@f$
     */
    const double deltaGammaH2l2vRatio2() const;
    /**
     * @brief The ratio of the Br@f$(H\to 2l2v)@f$ (@f$l=e,\mu@f$) in the current model
     * and in the Standard Model.
     * @return Br@f$(H\to 2l2v)@f$/Br@f$(H\to 2l2v)_{\mathrm{SM}}@f$
     */
    virtual const double BrH2l2vRatio() const;

    ///////////////////////OTHER DEDICATED (SEMI-)LEPTONIC 4 FERMION DECAYS/////////////////////////

    /**
     * @brief The ratio of the @f$\Gamma(H\to l l j j)@f$ (@f$l=e,\mu@f,~~j\not=b$) in the current model
     * and in the Standard Model.
     * @return @f$\Gamma(H\to l l j j)@f$/@f$\Gamma(H\to l l j j)_{\mathrm{SM}}@f$
     */
    const double GammaHlljjRatio() const;
    /**
     * @brief The new physics contribution to the ratio of the @f$\Gamma(H\to l l j j)@f$ (@f$l=e,\mu,~~j\not=b@f$)
     * in the current model and in the Standard Model. (Only terms that are
     * linear in the effective Lagrangian coefficients.)
     * @return @f$\delta \Gamma(H\to l l j j)@f$/@f$\Gamma(H\to l l j j)_{\mathrm{SM}}@f$
     */
    const double deltaGammaHlljjRatio1() const;
    /**
     * @brief The new physics contribution to the ratio of the @f$\Gamma(H\to l l j j)@f$ (@f$l=e,\mu,~~j\not=b@f$)
     * in the current model and in the Standard Model. (Only terms that are
     * quadratic in the effective Lagrangian coefficients.)
     * @return @f$\delta \Gamma(H\to l l j j)@f$/@f$\Gamma(H\to l l j j)_{\mathrm{SM}}@f$
     */
    const double deltaGammaHlljjRatio2() const;
    /**
     * @brief The ratio of the Br@f$(H\to l l j j)@f$ (@f$l=e,\mu,~~j\not=b@f$) in the current model
     * and in the Standard Model.
     * @return Br@f$(H\to l l j j)@f$/Br@f$(H\to l l j j)_{\mathrm{SM}}@f$
     */
    virtual const double BrHlljjRatio() const;


    /**
     * @brief The ratio of the @f$\Gamma(H\to l \nu j j)@f$ (@f$l=e,\mu@f,~~j\not=b$) in the current model
     * and in the Standard Model.
     * @return @f$\Gamma(H\to l \nu j j)@f$/@f$\Gamma(H\to l \nu j j)_{\mathrm{SM}}@f$
     */
    const double GammaHlvjjRatio() const;
    /**
     * @brief The new physics contribution to the ratio of the @f$\Gamma(H\to l \nu j j)@f$ (@f$l=e,\mu,~~j\not=b@f$)
     * in the current model and in the Standard Model. (Only terms that are
     * linear in the effective Lagrangian coefficients.)
     * @return @f$\delta \Gamma(H\to l \nu j j)@f$/@f$\Gamma(H\to l \nu j j)_{\mathrm{SM}}@f$
     */
    const double deltaGammaHlvjjRatio1() const;
    /**
     * @brief The new physics contribution to the ratio of the @f$\Gamma(H\to l \nu j j)@f$ (@f$l=e,\mu,~~j\not=b@f$)
     * in the current model and in the Standard Model. (Only terms that are
     * quadratic in the effective Lagrangian coefficients.)
     * @return @f$\delta \Gamma(H\to l \nu j j)@f$/@f$\Gamma(H\to l \nu j j)_{\mathrm{SM}}@f$
     */
    const double deltaGammaHlvjjRatio2() const;
    /**
     * @brief The ratio of the Br@f$(H\to l \nu j j)@f$ (@f$l=e,\mu,~~j\not=b@f$) in the current model
     * and in the Standard Model.
     * @return Br@f$(H\to l \nu j j)@f$/Br@f$(H\to l \nu j j)_{\mathrm{SM}}@f$
     */
    virtual const double BrHlvjjRatio() const;


    /**
     * @brief The ratio of the @f$\Gamma(H\to l \nu l \nu, l \nu j j)@f$ (@f$l=e,\mu,~~j\not=b@f$) in the current model
     * and in the Standard Model.
     * @return @f$\Gamma(H\to l \nu l \nu, l \nu j j)@f$/@f$\Gamma(H\to l \nu l \nu, l \nu j j)_{\mathrm{SM}}@f$
     */
    const double GammaHlv_lvorjjRatio() const;
    /**
     * @brief The new physics contribution to the ratio of the @f$\Gamma(H\to l \nu l \nu, l \nu j j)@f$ (@f$l=e,\mu,~~j\not=b@f$)
     * in the current model and in the Standard Model. (Only terms that are
     * linear in the effective Lagrangian coefficients.)
     * @return @f$\delta \Gamma(H\to l \nu l \nu, l \nu j j)@f$/@f$\Gamma(H\to l \nu l \nu, l \nu j j)_{\mathrm{SM}}@f$
     */
    const double deltaGammaHlv_lvorjjRatio1() const;
    /**
     * @brief The new physics contribution to the ratio of the @f$\Gamma(H\to l \nu l \nu, l \nu j j)@f$ (@f$l=e,\mu,~~j\not=b@f$)
     * in the current model and in the Standard Model. (Only terms that are
     * quadratic in the effective Lagrangian coefficients.)
     * @return @f$\delta \Gamma(H\to l \nu l \nu, l \nu j j)@f$/@f$\Gamma(H\to l \nu l \nu, l \nu j j)_{\mathrm{SM}}@f$
     */
    const double deltaGammaHlv_lvorjjRatio2() const;
    /**
     * @brief The ratio of the Br@f$(H\to l \nu l \nu, l \nu j j)@f$ (@f$l=e,\mu,~~j\not=b@f$) in the current model
     * and in the Standard Model.
     * @return Br@f$(H\to l \nu l \nu, l \nu j j)@f$/Br@f$(H\to l \nu l \nu, l \nu j j)_{\mathrm{SM}}@f$
     */
    virtual const double BrHlv_lvorjjRatio() const;


    /**
     * @brief The ratio of the @f$\Gamma(H\to l l \nu\nu, l l j j)@f$ (@f$l=e,\mu,~~j\not=b@f$) in the current model
     * and in the Standard Model.
     * @return @f$\Gamma(H\to l l \nu\nu, l l j j)@f$/@f$\Gamma(H\to l l \nu\nu, l l j j)_{\mathrm{SM}}@f$
     */
    const double GammaHll_vvorjjRatio() const;
    /**
     * @brief The new physics contribution to the ratio of the @f$\Gamma(H\to l l \nu\nu, l l j j)@f$ (@f$l=e,\mu,~~j\not=b@f$)
     * in the current model and in the Standard Model. (Only terms that are
     * linear in the effective Lagrangian coefficients.)
     * @return @f$\delta \Gamma(H\to l l \nu\nu, l l j j)@f$/@f$\Gamma(H\to l l \nu\nu, l l j j)_{\mathrm{SM}}@f$
     */
    const double deltaGammaHll_vvorjjRatio1() const;
    /**
     * @brief The new physics contribution to the ratio of the @f$\Gamma(H\to l l \nu\nu, l l j j)@f$ (@f$l=e,\mu,~~j\not=b@f$)
     * in the current model and in the Standard Model. (Only terms that are
     * quadratic in the effective Lagrangian coefficients.)
     * @return @f$\delta \Gamma(H\to l l \nu\nu, l l j j)@f$/@f$\Gamma(H\to l l \nu\nu, l l j j)_{\mathrm{SM}}@f$
     */
    const double deltaGammaHll_vvorjjRatio2() const;
    /**
     * @brief The ratio of the Br@f$(H\to l l \nu\nu, l l j j)@f$ (@f$l=e,\mu,~~j\not=b@f$) in the current model
     * and in the Standard Model.
     * @return Br@f$(H\to l l \nu\nu, l l j j)@f$/Br@f$(H\to l l \nu\nu, l l j j)_{\mathrm{SM}}@f$
     */
    virtual const double BrHll_vvorjjRatio() const;

    ///////////////////////OTHER HIGGS BRANCHING RATIOS/////////////////////////    

    /**
     * @brief The branching ratio of the of the Higgs into exotic particles.
     * @return Br@f$(H\to exotic)@f$
     */
    virtual const double Br_H_exo() const;
    /**
     * @brief The branching ratio of the of the Higgs into invisible particles.
     * @return Br@f$(H\to invisible)@f$
     */
    virtual const double Br_H_inv() const;
    /**
     * @brief The branching ratio of the of the Higgs into invisible particles 
     * (only invisible new particles).
     * @return Br@f$(H\to invisible,NP)@f$
     */
    virtual const double Br_H_inv_NP() const;
    /**
     * @brief The ratio of the Br@f$(H\to visible)@f$ in the current model
     * and in the Standard Model.
     * @return Br@f$(H\to visible)@f$/Br@f$(H\to visible)_{\mathrm{SM}}@f$
     */
    virtual const double BrHvisRatio() const;
    /**
     * @brief The ratio of the Br@f$(H\to invisible)@f$ in the current model
     * and in the Standard Model.
     * @return Br@f$(H\to invisible)@f$/Br@f$(H\to ZZ \to invisible)_{\mathrm{SM}}@f$
     */
    virtual const double BrHtoinvRatio() const;
    
    ///////////////////////SPECIAL OBSERVABLES/////////////////////////

    /**
     * @brief The ratio @f$\sigma(ttH)/\sigma(ttZ)@f$ 
     * in the @f$H,Z\to b\bar{b}@f$ channel in the current model and in the Standard Model.
     * @param[in] sqrt_s the center-of-mass energy in TeV
     * @return @f$\sigma(ttH)/\sigma(ttZ)@f$ normalized to the SM
     */
    virtual const double muttHZbbboost(double sqrt_s) const;
    
    /**
     * @brief The ratio @f$\sigma(ttH)/\sigma(ttZ)@f$ 
     * in the @f$H\to b\bar{b}@f$, @f$Z\to e^+e^-@f$ channel channel in the current model and in the Standard Model.
     * @param[in] sqrt_s the center-of-mass energy in TeV
     * @return @f$\sigma(ttH)/\sigma(ttZ)@f$ normalized to the SM
     */
    virtual const double muttHgagaZeeboost(const double sqrt_s) const;

    /**
     * @brief The ratio @f$\mu_{ggH,\gamma\gamma}@f$ between the gluon-gluon fusion Higgs
     * production cross-section with subsequent decay into 2 photons in the
     * current model and in the Standard Model.
     * @param[in] sqrt_s the center-of-mass energy in TeV
     * @return @f$\mu_{ggH,\gamma\gamma}@f$
     */
    virtual const double muggHgaga(double sqrt_s) const;
    /**
     * @brief The ratio @f$\mu_{VBF,\gamma\gamma}@f$ between the VBF Higgs
     * production cross-section with subsequent decay into 2 photons in the
     * current model and in the Standard Model.
     * @param[in] sqrt_s the center-of-mass energy in TeV
     * @return @f$\mu_{VBF,\gamma\gamma}@f$
     */
    virtual const double muVBFHgaga(double sqrt_s) const;
    /**
     * @brief The ratio @f$\mu_{ZH,\gamma\gamma}@f$ between the ZH
     * production cross-section with subsequent decay into 2 photons in the
     * current model and in the Standard Model.
     * @param[in] sqrt_s the center-of-mass energy in TeV
     * @return @f$\mu_{ZH,\gamma\gamma}@f$
     */
    virtual const double muZHgaga(double sqrt_s) const;
    /**
     * @brief The ratio @f$\mu_{WH,\gamma\gamma}@f$ between the WH
     * production cross-section with subsequent decay into 2 photons in the
     * current model and in the Standard Model.
     * @param[in] sqrt_s the center-of-mass energy in TeV
     * @return @f$\mu_{WH,\gamma\gamma}@f$
     */
    virtual const double muWHgaga(double sqrt_s) const;
    /**
     * @brief The ratio @f$\mu_{VH,\gamma\gamma}@f$ between the VH
     * production cross-section with subsequent decay into 2 photons in the
     * current model and in the Standard Model.
     * @param[in] sqrt_s the center-of-mass energy in TeV
     * @return @f$\mu_{VH,\gamma\gamma}@f$
     */
    virtual const double muVHgaga(double sqrt_s) const;
    /**
     * @brief The ratio @f$\mu_{ttH,\gamma\gamma}@f$ between the ttH
     * production cross-section with subsequent decay into 2 photons in the
     * current model and in the Standard Model.
     * @param[in] sqrt_s the center-of-mass energy in TeV
     * @return @f$\mu_{ttH,\gamma\gamma}@f$
     */
    virtual const double muttHgaga(double sqrt_s) const;
    /**
     * @brief The ratio @f$\mu_{ggH,Z\gamma}@f$ between the gluon-gluon fusion Higgs
     * production cross-section with subsequent decay into @f$Z \gamma@f$ in the
     * current model and in the Standard Model.
     * @param[in] sqrt_s the center-of-mass energy in TeV
     * @return @f$\mu_{ggH,Z\gamma}@f$
     */
    virtual const double muggHZga(double sqrt_s) const;
    /**
     * @brief The ratio @f$\mu_{VBF,Z\gamma}@f$ between the VBF Higgs
     * production cross-section with subsequent decay into @f$Z \gamma@f$ in the
     * current model and in the Standard Model.
     * @param[in] sqrt_s the center-of-mass energy in TeV
     * @return @f$\mu_{VBF,Z\gamma}@f$
     */
    virtual const double muVBFHZga(double sqrt_s) const;
    /**
     * @brief The ratio @f$\mu_{ZH,Z\gamma}@f$ between the ZH
     * production cross-section with subsequent decay into @f$Z \gamma@f$ in the
     * current model and in the Standard Model.
     * @param[in] sqrt_s the center-of-mass energy in TeV
     * @return @f$\mu_{ZH,Z\gamma}@f$
     */
    virtual const double muZHZga(double sqrt_s) const;
    /**
     * @brief The ratio @f$\mu_{WH,Z\gamma}@f$ between the WH
     * production cross-section with subsequent decay into @f$Z \gamma@f$ in the
     * current model and in the Standard Model.
     * @param[in] sqrt_s the center-of-mass energy in TeV
     * @return @f$\mu_{WH,Z\gamma}@f$
     */
    virtual const double muWHZga(double sqrt_s) const;
    /**
     * @brief The ratio @f$\mu_{VH,Z\gamma}@f$ between the VH
     * production cross-section with subsequent decay into @f$Z \gamma@f$ in the
     * current model and in the Standard Model.
     * @param[in] sqrt_s the center-of-mass energy in TeV
     * @return @f$\mu_{VH,Z\gamma}@f$
     */
    virtual const double muVHZga(double sqrt_s) const;
    /**
     * @brief The ratio @f$\mu_{ttH,Z\gamma}@f$ between the ttH
     * production cross-section with subsequent decay into @f$Z \gamma@f$ in the
     * current model and in the Standard Model.
     * @param[in] sqrt_s the center-of-mass energy in TeV
     * @return @f$\mu_{ttH,Z\gamma}@f$
     */
    virtual const double muttHZga(double sqrt_s) const;
    /**
     * @brief The ratio @f$\mu_{ggH,ZZ}@f$ between the gluon-gluon fusion Higgs
     * production cross-section with subsequent decay into @f$Z Z^*@f$ in the
     * current model and in the Standard Model.
     * @param[in] sqrt_s the center-of-mass energy in TeV
     * @return @f$\mu_{ggH,ZZ}@f$
     */
    virtual const double muggHZZ(double sqrt_s) const;
    /**
     * @brief The ratio @f$\mu_{VBF,ZZ}@f$ between the VBF Higgs
     * production cross-section with subsequent decay into @f$Z Z^*@f$ in the
     * current model and in the Standard Model.
     * @param[in] sqrt_s the center-of-mass energy in TeV
     * @return @f$\mu_{VBF,ZZ}@f$
     */
    virtual const double muVBFHZZ(double sqrt_s) const;
    /**
     * @brief The ratio @f$\mu_{ZH,ZZ}@f$ between the ZH
     * production cross-section with subsequent decay into @f$Z Z^*@f$ in the
     * current model and in the Standard Model.
     * @param[in] sqrt_s the center-of-mass energy in TeV
     * @return @f$\mu_{ZH,ZZ}@f$
     */
    virtual const double muZHZZ(double sqrt_s) const;
    /**
     * @brief The ratio @f$\mu_{WH,ZZ}@f$ between the WH
     * production cross-section with subsequent decay into @f$Z Z^*@f$ in the
     * current model and in the Standard Model.
     * @param[in] sqrt_s the center-of-mass energy in TeV
     * @return @f$\mu_{WH,ZZ}@f$
     */
    virtual const double muWHZZ(double sqrt_s) const;
    /**
     * @brief The ratio @f$\mu_{VH,ZZ}@f$ between the VH
     * production cross-section with subsequent decay into @f$Z Z^*@f$ in the
     * current model and in the Standard Model.
     * @param[in] sqrt_s the center-of-mass energy in TeV
     * @return @f$\mu_{VH,ZZ}@f$
     */
    virtual const double muVHZZ(double sqrt_s) const;
    /**
     * @brief The ratio @f$\mu_{ttH,ZZ}@f$ between the ttH
     * production cross-section with subsequent decay into @f$Z Z^*@f$ in the
     * current model and in the Standard Model.
     * @param[in] sqrt_s the center-of-mass energy in TeV
     * @return @f$\mu_{ttH,ZZ}@f$
     */
    virtual const double muttHZZ(double sqrt_s) const;

    /**
     * @brief The ratio @f$\mu_{ggH,ZZ\to 4l}@f$ between the gluon-gluon fusion Higgs
     * production cross-section with subsequent decay into @f$Z Z^*\to 4l@f$ in the
     * current model and in the Standard Model.
     * @param[in] sqrt_s the center-of-mass energy in TeV
     * @return @f$\mu_{ggH,ZZ\to 4l}@f$
     */
    virtual const double muggHZZ4l(double sqrt_s) const;
    /**
     * @brief The ratio @f$\mu_{VBF,ZZ\to 4l}@f$ between the VBF Higgs
     * production cross-section with subsequent decay into @f$Z Z^*\to 4l@f$ in the
     * current model and in the Standard Model.
     * @param[in] sqrt_s the center-of-mass energy in TeV
     * @return @f$\mu_{VBF,ZZ\to 4l}@f$
     */
    virtual const double muVBFHZZ4l(double sqrt_s) const;
    /**
     * @brief The ratio @f$\mu_{ZH,ZZ\to 4l}@f$ between the ZH
     * production cross-section with subsequent decay into @f$Z Z^*\to 4l@f$ in the
     * current model and in the Standard Model.
     * @param[in] sqrt_s the center-of-mass energy in TeV
     * @return @f$\mu_{ZH,ZZ\to 4l}@f$
     */
    virtual const double muZHZZ4l(double sqrt_s) const;
    /**
     * @brief The ratio @f$\mu_{WH,ZZ\to 4l}@f$ between the WH
     * production cross-section with subsequent decay into @f$Z Z^*\to 4l@f$ in the
     * current model and in the Standard Model.
     * @param[in] sqrt_s the center-of-mass energy in TeV
     * @return @f$\mu_{WH,ZZ\to 4l}@f$
     */
    virtual const double muWHZZ4l(double sqrt_s) const;
    /**
     * @brief The ratio @f$\mu_{VH,ZZ\to 4l}@f$ between the VH
     * production cross-section with subsequent decay into @f$Z Z^*\to 4l@f$ in the
     * current model and in the Standard Model.
     * @param[in] sqrt_s the center-of-mass energy in TeV
     * @return @f$\mu_{VH,ZZ\to 4l}@f$
     */
    virtual const double muVHZZ4l(double sqrt_s) const;
    /**
     * @brief The ratio @f$\mu_{ttH,ZZ\to 4l}@f$ between the ttH
     * production cross-section with subsequent decay into @f$Z Z^*\to 4l@f$ in the
     * current model and in the Standard Model.
     * @param[in] sqrt_s the center-of-mass energy in TeV
     * @return @f$\mu_{ttH,ZZ\to 4l}@f$
     */
    virtual const double muttHZZ4l(double sqrt_s) const;

    /**
     * @brief The ratio @f$\mu_{ggH,WW}@f$ between the gluon-gluon fusion Higgs
     * production cross-section with subsequent decay into @f$W W^*@f$ in the
     * current model and in the Standard Model.
     * @param[in] sqrt_s the center-of-mass energy in TeV
     * @return @f$\mu_{ggH,WW}@f$
     */
    virtual const double muggHWW(double sqrt_s) const;
    /**
     * @brief The ratio @f$\mu_{VBF,WW}@f$ between the VBF Higgs
     * production cross-section with subsequent decay into @f$W W^*@f$ in the
     * current model and in the Standard Model.
     * @param[in] sqrt_s the center-of-mass energy in TeV
     * @return @f$\mu_{VBF,WW}@f$
     */
    virtual const double muVBFHWW(double sqrt_s) const;
    /**
     * @brief The ratio @f$\mu_{ZH,WW}@f$ between the ZH
     * production cross-section with subsequent decay into @f$W W^*@f$ in the
     * current model and in the Standard Model.
     * @param[in] sqrt_s the center-of-mass energy in TeV
     * @return @f$\mu_{ZH,WW}@f$
     */
    virtual const double muZHWW(double sqrt_s) const;
    /**
     * @brief The ratio @f$\mu_{WH,WW}@f$ between the WH
     * production cross-section with subsequent decay into @f$W W^*@f$ in the
     * current model and in the Standard Model.
     * @param[in] sqrt_s the center-of-mass energy in TeV
     * @return @f$\mu_{WH,WW}@f$
     */
    virtual const double muWHWW(double sqrt_s) const;
    /**
     * @brief The ratio @f$\mu_{VH,WW}@f$ between the VH
     * production cross-section with subsequent decay into @f$W W^*@f$ in the
     * current model and in the Standard Model.
     * @param[in] sqrt_s the center-of-mass energy in TeV
     * @return @f$\mu_{VH,WW}@f$
     */
    virtual const double muVHWW(double sqrt_s) const;
    /**
     * @brief The ratio @f$\mu_{ttH,WW}@f$ between the ttH
     * production cross-section with subsequent decay into @f$W W^*@f$ in the
     * current model and in the Standard Model.
     * @param[in] sqrt_s the center-of-mass energy in TeV
     * @return @f$\mu_{ttH,WW}@f$
     */
    virtual const double muttHWW(double sqrt_s) const;
    /**
     * @brief The ratio @f$\mu_{ggH,WW\to 2l2\nu}@f$ between the gluon-gluon fusion Higgs
     * production cross-section with subsequent decay into @f$W W^*\to 2l2\nu@f$ in the
     * current model and in the Standard Model.
     * @param[in] sqrt_s the center-of-mass energy in TeV
     * @return @f$\mu_{ggH,WW\to 2l2\nu}@f$
     */
    virtual const double muggHWW2l2v(double sqrt_s) const;
    /**
     * @brief The ratio @f$\mu_{VBF,WW\to 2l2\nu}@f$ between the VBF Higgs
     * production cross-section with subsequent decay into @f$W W^*\to 2l2\nu@f$ in the
     * current model and in the Standard Model.
     * @param[in] sqrt_s the center-of-mass energy in TeV
     * @return @f$\mu_{VBF,WW\to 2l2\nu}@f$
     */
    virtual const double muVBFHWW2l2v(double sqrt_s) const;
    /**
     * @brief The ratio @f$\mu_{ZH,WW\to 2l2\nu}@f$ between the ZH
     * production cross-section with subsequent decay into @f$W W^*\to 2l2\nu@f$ in the
     * current model and in the Standard Model.
     * @param[in] sqrt_s the center-of-mass energy in TeV
     * @return @f$\mu_{ZH,WW\to 2l2\nu}@f$
     */
    virtual const double muZHWW2l2v(double sqrt_s) const;
    /**
     * @brief The ratio @f$\mu_{WH,WW\to 2l2\nu}@f$ between the WH
     * production cross-section with subsequent decay into @f$W W^*\to 2l2\nu@f$ in the
     * current model and in the Standard Model.
     * @param[in] sqrt_s the center-of-mass energy in TeV
     * @return @f$\mu_{WH,WW\to 2l2\nu}@f$
     */
    virtual const double muWHWW2l2v(double sqrt_s) const;
    /**
     * @brief The ratio @f$\mu_{VH,WW\to 2l2\nu}@f$ between the VH
     * production cross-section with subsequent decay into @f$W W^*\to 2l2\nu@f$ in the
     * current model and in the Standard Model.
     * @param[in] sqrt_s the center-of-mass energy in TeV
     * @return @f$\mu_{VH,WW\to 2l2\nu}@f$
     */
    virtual const double muVHWW2l2v(double sqrt_s) const;
    /**
     * @brief The ratio @f$\mu_{ttH,WW\to 2l2\nu}@f$ between the ttH
     * production cross-section with subsequent decay into @f$W W^*\to 2l2\nu@f$ in the
     * current model and in the Standard Model.
     * @param[in] sqrt_s the center-of-mass energy in TeV
     * @return @f$\mu_{ttH,WW\to 2l2\nu}@f$
     */
    virtual const double muttHWW2l2v(double sqrt_s) const;
    /**
     * @brief The ratio @f$\mu_{ggH,\mu\mu}@f$ between the gluon-gluon fusion Higgs
     * production cross-section with subsequent decay into @f$\mu\mu@f$ in the
     * current model and in the Standard Model.
     * @param[in] sqrt_s the center-of-mass energy in TeV
     * @return @f$\mu_{ggH,\mu\mu}@f$
     */
    virtual const double muggHmumu(double sqrt_s) const;
    /**
     * @brief The ratio @f$\mu_{VBF,\mu\mu}@f$ between the VBF Higgs
     * production cross-section with subsequent decay into @f$\mu\mu@f$ in the
     * current model and in the Standard Model.
     * @param[in] sqrt_s the center-of-mass energy in TeV
     * @return @f$\mu_{VBF,\mu\mu}@f$
     */
    virtual const double muVBFHmumu(double sqrt_s) const;
    /**
     * @brief The ratio @f$\mu_{ZH,\mu\mu}@f$ between the ZH
     * production cross-section with subsequent decay into @f$\mu\mu@f$ in the
     * current model and in the Standard Model.
     * @param[in] sqrt_s the center-of-mass energy in TeV
     * @return @f$\mu_{ZH,\mu\mu}@f$
     */
    virtual const double muZHmumu(double sqrt_s) const;
    /**
     * @brief The ratio @f$\mu_{WH,\mu\mu}@f$ between the WH
     * production cross-section with subsequent decay into @f$\mu\mu@f$ in the
     * current model and in the Standard Model.
     * @param[in] sqrt_s the center-of-mass energy in TeV
     * @return @f$\mu_{WH,\mu\mu}@f$
     */
    virtual const double muWHmumu(double sqrt_s) const;
    /**
     * @brief The ratio @f$\mu_{VH,\mu\mu}@f$ between the VH
     * production cross-section with subsequent decay into @f$\mu\mu@f$ in the
     * current model and in the Standard Model.
     * @param[in] sqrt_s the center-of-mass energy in TeV
     * @return @f$\mu_{VH,\mu\mu}@f$
     */
    virtual const double muVHmumu(double sqrt_s) const;
    /**
     * @brief The ratio @f$\mu_{ttH,\mu\mu}@f$ between the ttH
     * production cross-section with subsequent decay into @f$\mu\mu@f$ in the
     * current model and in the Standard Model.
     * @param[in] sqrt_s the center-of-mass energy in TeV
     * @return @f$\mu_{ttH,\mu\mu}@f$
     */
    virtual const double muttHmumu(double sqrt_s) const;
    /**
     * @brief The ratio @f$\mu_{ggH,\tau\tau}@f$ between the gluon-gluon fusion Higgs
     * production cross-section with subsequent decay into @f$\tau\tau@f$ in the
     * current model and in the Standard Model.
     * @param[in] sqrt_s the center-of-mass energy in TeV
     * @return @f$\mu_{ggH,\tau\tau}@f$
     */
    virtual const double muggHtautau(double sqrt_s) const;
    /**
     * @brief The ratio @f$\mu_{VBF,\tau\tau}@f$ between the VBF Higgs
     * production cross-section with subsequent decay into @f$\tau\tau@f$ in the
     * current model and in the Standard Model.
     * @param[in] sqrt_s the center-of-mass energy in TeV
     * @return @f$\mu_{VBF,\tau\tau}@f$
     */
    virtual const double muVBFHtautau(double sqrt_s) const;
    /**
     * @brief The ratio @f$\mu_{ZH,\tau\tau}@f$ between the ZH
     * production cross-section with subsequent decay into @f$\tau\tau@f$ in the
     * current model and in the Standard Model.
     * @param[in] sqrt_s the center-of-mass energy in TeV
     * @return @f$\mu_{ZH,\tau\tau}@f$
     */
    virtual const double muZHtautau(double sqrt_s) const;
    /**
     * @brief The ratio @f$\mu_{WH,\tau\tau}@f$ between the WH
     * production cross-section with subsequent decay into @f$\tau\tau@f$ in the
     * current model and in the Standard Model.
     * @param[in] sqrt_s the center-of-mass energy in TeV
     * @return @f$\mu_{WH,\tau\tau}@f$
     */
    virtual const double muWHtautau(double sqrt_s) const;
    /**
     * @brief The ratio @f$\mu_{VH,\tau\tau}@f$ between the VH
     * production cross-section with subsequent decay into @f$\tau\tau@f$ in the
     * current model and in the Standard Model.
     * @param[in] sqrt_s the center-of-mass energy in TeV
     * @return @f$\mu_{VH,\tau\tau}@f$
     */
    virtual const double muVHtautau(double sqrt_s) const;
    /**
     * @brief The ratio @f$\mu_{ttH,\tau\tau}@f$ between the ttH
     * production cross-section with subsequent decay into @f$\tau\tau@f$ in the
     * current model and in the Standard Model.
     * @param[in] sqrt_s the center-of-mass energy in TeV
     * @return @f$\mu_{ttH,\tau\tau}@f$
     */
    virtual const double muttHtautau(double sqrt_s) const;
    /**
     * @brief The ratio @f$\mu_{ggH,bb}@f$ between the gluon-gluon fusion Higgs
     * production cross-section with subsequent decay into @f$bb@f$ in the
     * current model and in the Standard Model.
     * @param[in] sqrt_s the center-of-mass energy in TeV
     * @return @f$\mu_{ggH,bb}@f$
     */
    virtual const double muggHbb(double sqrt_s) const;
    /**
     * @brief The ratio @f$\mu_{VBF,bb}@f$ between the VBF Higgs
     * production cross-section with subsequent decay into @f$bb@f$ in the
     * current model and in the Standard Model.
     * @param[in] sqrt_s the center-of-mass energy in TeV
     * @return @f$\mu_{VBF,bb}@f$
     */
    virtual const double muVBFHbb(double sqrt_s) const;
    /**
     * @brief The ratio @f$\mu_{ZH,bb}@f$ between the ZH
     * production cross-section with subsequent decay into @f$bb@f$ in the
     * current model and in the Standard Model.
     * @param[in] sqrt_s the center-of-mass energy in TeV
     * @return @f$\mu_{ZH,bb}@f$
     */
    virtual const double muZHbb(double sqrt_s) const;
    /**
     * @brief The ratio @f$\mu_{WH,bb}@f$ between the WH
     * production cross-section with subsequent decay into @f$bb@f$ in the
     * current model and in the Standard Model.
     * @param[in] sqrt_s the center-of-mass energy in TeV
     * @return @f$\mu_{WH,bb}@f$
     */
    virtual const double muWHbb(double sqrt_s) const;
    /**
     * @brief The ratio @f$\mu_{VH,bb}@f$ between the VH
     * production cross-section with subsequent decay into @f$bb@f$ in the
     * current model and in the Standard Model.
     * @param[in] sqrt_s the center-of-mass energy in TeV
     * @return @f$\mu_{VH,bb}@f$
     */
    virtual const double muVHbb(double sqrt_s) const;
    /**
     * @brief The ratio @f$\mu_{ttH,bb}@f$ between the ttH
     * production cross-section with subsequent decay into @f$bb@f$ in the
     * current model and in the Standard Model.
     * @param[in] sqrt_s the center-of-mass energy in TeV
     * @return @f$\mu_{ttH,bb}@f$
     */
    virtual const double muttHbb(double sqrt_s) const;

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
    virtual const double muTHUggHgaga(double sqrt_s) const;
    /**
     * @brief The ratio @f$\mu_{VBF,\gamma\gamma}@f$ between the VBF Higgs
     * production cross-section with subsequent decay into 2 photons in the
     * current model and in the Standard Model.
     * @param[in] sqrt_s the center-of-mass energy in TeV
     * @return @f$\mu_{VBF,\gamma\gamma}@f$
     */
    virtual const double muTHUVBFHgaga(double sqrt_s) const;
    /**
     * @brief The ratio @f$\mu_{ZH,\gamma\gamma}@f$ between the ZH
     * production cross-section with subsequent decay into 2 photons in the
     * current model and in the Standard Model.
     * @param[in] sqrt_s the center-of-mass energy in TeV
     * @return @f$\mu_{ZH,\gamma\gamma}@f$
     */
    virtual const double muTHUZHgaga(double sqrt_s) const;
    /**
     * @brief The ratio @f$\mu_{WH,\gamma\gamma}@f$ between the WH
     * production cross-section with subsequent decay into 2 photons in the
     * current model and in the Standard Model.
     * @param[in] sqrt_s the center-of-mass energy in TeV
     * @return @f$\mu_{WH,\gamma\gamma}@f$
     */
    virtual const double muTHUWHgaga(double sqrt_s) const;
    /**
     * @brief The ratio @f$\mu_{VH,\gamma\gamma}@f$ between the VH
     * production cross-section with subsequent decay into 2 photons in the
     * current model and in the Standard Model.
     * @param[in] sqrt_s the center-of-mass energy in TeV
     * @return @f$\mu_{VH,\gamma\gamma}@f$
     */
    virtual const double muTHUVHgaga(double sqrt_s) const;
    /**
     * @brief The ratio @f$\mu_{ttH,\gamma\gamma}@f$ between the ttH
     * production cross-section with subsequent decay into 2 photons in the
     * current model and in the Standard Model.
     * @param[in] sqrt_s the center-of-mass energy in TeV
     * @return @f$\mu_{ttH,\gamma\gamma}@f$
     */
    virtual const double muTHUttHgaga(double sqrt_s) const;
    /**
     * @brief The ratio @f$\mu_{ggH,Z\gamma}@f$ between the gluon-gluon fusion Higgs
     * production cross-section with subsequent decay into @f$Z \gamma@f$ in the
     * current model and in the Standard Model.
     * @param[in] sqrt_s the center-of-mass energy in TeV
     * @return @f$\mu_{ggH,Z\gamma}@f$
     */
    virtual const double muTHUggHZga(double sqrt_s) const;
    /**
     * @brief The ratio @f$\mu_{VBF,Z\gamma}@f$ between the VBF Higgs
     * production cross-section with subsequent decay into @f$Z \gamma@f$ in the
     * current model and in the Standard Model.
     * @param[in] sqrt_s the center-of-mass energy in TeV
     * @return @f$\mu_{VBF,Z\gamma}@f$
     */
    virtual const double muTHUVBFHZga(double sqrt_s) const;
    /**
     * @brief The ratio @f$\mu_{ZH,Z\gamma}@f$ between the ZH
     * production cross-section with subsequent decay into @f$Z \gamma@f$ in the
     * current model and in the Standard Model.
     * @param[in] sqrt_s the center-of-mass energy in TeV
     * @return @f$\mu_{ZH,Z\gamma}@f$
     */
    virtual const double muTHUZHZga(double sqrt_s) const;
    /**
     * @brief The ratio @f$\mu_{WH,Z\gamma}@f$ between the WH
     * production cross-section with subsequent decay into @f$Z \gamma@f$ in the
     * current model and in the Standard Model.
     * @param[in] sqrt_s the center-of-mass energy in TeV
     * @return @f$\mu_{WH,Z\gamma}@f$
     */
    virtual const double muTHUWHZga(double sqrt_s) const;
    /**
     * @brief The ratio @f$\mu_{VH,Z\gamma}@f$ between the VH
     * production cross-section with subsequent decay into @f$Z \gamma@f$ in the
     * current model and in the Standard Model.
     * @param[in] sqrt_s the center-of-mass energy in TeV
     * @return @f$\mu_{VH,Z\gamma}@f$
     */
    virtual const double muTHUVHZga(double sqrt_s) const;
    /**
     * @brief The ratio @f$\mu_{ttH,Z\gamma}@f$ between the ttH
     * production cross-section with subsequent decay into @f$Z \gamma@f$ in the
     * current model and in the Standard Model.
     * @param[in] sqrt_s the center-of-mass energy in TeV
     * @return @f$\mu_{ttH,Z\gamma}@f$
     */
    virtual const double muTHUttHZga(double sqrt_s) const;
    /**
     * @brief The ratio @f$\mu_{ggH,ZZ}@f$ between the gluon-gluon fusion Higgs
     * production cross-section with subsequent decay into @f$Z Z^*@f$ in the
     * current model and in the Standard Model.
     * @param[in] sqrt_s the center-of-mass energy in TeV
     * @return @f$\mu_{ggH,ZZ}@f$
     */
    virtual const double muTHUggHZZ(double sqrt_s) const;
    /**
     * @brief The ratio @f$\mu_{VBF,ZZ}@f$ between the VBF Higgs
     * production cross-section with subsequent decay into @f$Z Z^*@f$ in the
     * current model and in the Standard Model.
     * @param[in] sqrt_s the center-of-mass energy in TeV
     * @return @f$\mu_{VBF,ZZ}@f$
     */
    virtual const double muTHUVBFHZZ(double sqrt_s) const;
    /**
     * @brief The ratio @f$\mu_{ZH,ZZ}@f$ between the ZH
     * production cross-section with subsequent decay into @f$Z Z^*@f$ in the
     * current model and in the Standard Model.
     * @param[in] sqrt_s the center-of-mass energy in TeV
     * @return @f$\mu_{ZH,ZZ}@f$
     */
    virtual const double muTHUZHZZ(double sqrt_s) const;
    /**
     * @brief The ratio @f$\mu_{WH,ZZ}@f$ between the WH
     * production cross-section with subsequent decay into @f$Z Z^*@f$ in the
     * current model and in the Standard Model.
     * @param[in] sqrt_s the center-of-mass energy in TeV
     * @return @f$\mu_{WH,ZZ}@f$
     */
    virtual const double muTHUWHZZ(double sqrt_s) const;
    /**
     * @brief The ratio @f$\mu_{VH,ZZ}@f$ between the VH
     * production cross-section with subsequent decay into @f$Z Z^*@f$ in the
     * current model and in the Standard Model.
     * @param[in] sqrt_s the center-of-mass energy in TeV
     * @return @f$\mu_{VH,ZZ}@f$
     */
    virtual const double muTHUVHZZ(double sqrt_s) const;
    /**
     * @brief The ratio @f$\mu_{ttH,ZZ}@f$ between the ttH
     * production cross-section with subsequent decay into @f$Z Z^*@f$ in the
     * current model and in the Standard Model.
     * @param[in] sqrt_s the center-of-mass energy in TeV
     * @return @f$\mu_{ttH,ZZ}@f$
     */
    virtual const double muTHUttHZZ(double sqrt_s) const;

    /**
     * @brief The ratio @f$\mu_{ggH,ZZ\to 4l}@f$ between the gluon-gluon fusion Higgs
     * production cross-section with subsequent decay into @f$Z Z^*\to 4l@f$ in the
     * current model and in the Standard Model.
     * @param[in] sqrt_s the center-of-mass energy in TeV
     * @return @f$\mu_{ggH,ZZ\to 4l}@f$
     */
    virtual const double muTHUggHZZ4l(double sqrt_s) const;
    /**
     * @brief The ratio @f$\mu_{VBF,ZZ\to 4l}@f$ between the VBF Higgs
     * production cross-section with subsequent decay into @f$Z Z^*\to 4l@f$ in the
     * current model and in the Standard Model.
     * @param[in] sqrt_s the center-of-mass energy in TeV
     * @return @f$\mu_{VBF,ZZ\to 4l}@f$
     */
    virtual const double muTHUVBFHZZ4l(double sqrt_s) const;
    /**
     * @brief The ratio @f$\mu_{ZH,ZZ\to 4l}@f$ between the ZH
     * production cross-section with subsequent decay into @f$Z Z^*\to 4l@f$ in the
     * current model and in the Standard Model.
     * @param[in] sqrt_s the center-of-mass energy in TeV
     * @return @f$\mu_{ZH,ZZ\to 4l}@f$
     */
    virtual const double muTHUZHZZ4l(double sqrt_s) const;
    /**
     * @brief The ratio @f$\mu_{WH,ZZ\to 4l}@f$ between the WH
     * production cross-section with subsequent decay into @f$Z Z^*\to 4l@f$ in the
     * current model and in the Standard Model.
     * @param[in] sqrt_s the center-of-mass energy in TeV
     * @return @f$\mu_{WH,ZZ\to 4l}@f$
     */
    virtual const double muTHUWHZZ4l(double sqrt_s) const;
    /**
     * @brief The ratio @f$\mu_{VH,ZZ\to 4l}@f$ between the VH
     * production cross-section with subsequent decay into @f$Z Z^*\to 4l@f$ in the
     * current model and in the Standard Model.
     * @param[in] sqrt_s the center-of-mass energy in TeV
     * @return @f$\mu_{VH,ZZ\to 4l}@f$
     */
    virtual const double muTHUVHZZ4l(double sqrt_s) const;
    /**
     * @brief The ratio @f$\mu_{ttH,ZZ\to 4l}@f$ between the ttH
     * production cross-section with subsequent decay into @f$Z Z^*\to 4l@f$ in the
     * current model and in the Standard Model.
     * @param[in] sqrt_s the center-of-mass energy in TeV
     * @return @f$\mu_{ttH,ZZ\to 4l}@f$
     */
    virtual const double muTHUttHZZ4l(double sqrt_s) const;

    /**
     * @brief The ratio @f$\mu_{ggH,WW}@f$ between the gluon-gluon fusion Higgs
     * production cross-section with subsequent decay into @f$W W^*@f$ in the
     * current model and in the Standard Model.
     * @param[in] sqrt_s the center-of-mass energy in TeV
     * @return @f$\mu_{ggH,WW}@f$
     */
    virtual const double muTHUggHWW(double sqrt_s) const;
    /**
     * @brief The ratio @f$\mu_{VBF,WW}@f$ between the VBF Higgs
     * production cross-section with subsequent decay into @f$W W^*@f$ in the
     * current model and in the Standard Model.
     * @param[in] sqrt_s the center-of-mass energy in TeV
     * @return @f$\mu_{VBF,WW}@f$
     */
    virtual const double muTHUVBFHWW(double sqrt_s) const;
    /**
     * @brief The ratio @f$\mu_{ZH,WW}@f$ between the ZH
     * production cross-section with subsequent decay into @f$W W^*@f$ in the
     * current model and in the Standard Model.
     * @param[in] sqrt_s the center-of-mass energy in TeV
     * @return @f$\mu_{ZH,WW}@f$
     */
    virtual const double muTHUZHWW(double sqrt_s) const;
    /**
     * @brief The ratio @f$\mu_{WH,WW}@f$ between the WH
     * production cross-section with subsequent decay into @f$W W^*@f$ in the
     * current model and in the Standard Model.
     * @param[in] sqrt_s the center-of-mass energy in TeV
     * @return @f$\mu_{WH,WW}@f$
     */
    virtual const double muTHUWHWW(double sqrt_s) const;
    /**
     * @brief The ratio @f$\mu_{VH,WW}@f$ between the VH
     * production cross-section with subsequent decay into @f$W W^*@f$ in the
     * current model and in the Standard Model.
     * @param[in] sqrt_s the center-of-mass energy in TeV
     * @return @f$\mu_{VH,WW}@f$
     */
    virtual const double muTHUVHWW(double sqrt_s) const;
    /**
     * @brief The ratio @f$\mu_{ttH,WW}@f$ between the ttH
     * production cross-section with subsequent decay into @f$W W^*@f$ in the
     * current model and in the Standard Model.
     * @param[in] sqrt_s the center-of-mass energy in TeV
     * @return @f$\mu_{ttH,WW}@f$
     */
    virtual const double muTHUttHWW(double sqrt_s) const;
    /**
     * @brief The ratio @f$\mu_{ggH,WW\to 2l2\nu}@f$ between the gluon-gluon fusion Higgs
     * production cross-section with subsequent decay into @f$W W^*\to 2l2\nu@f$ in the
     * current model and in the Standard Model.
     * @param[in] sqrt_s the center-of-mass energy in TeV
     * @return @f$\mu_{ggH,WW\to 2l2\nu}@f$
     */
    virtual const double muTHUggHWW2l2v(double sqrt_s) const;
    /**
     * @brief The ratio @f$\mu_{VBF,WW\to 2l2\nu}@f$ between the VBF Higgs
     * production cross-section with subsequent decay into @f$W W^*\to 2l2\nu@f$ in the
     * current model and in the Standard Model.
     * @param[in] sqrt_s the center-of-mass energy in TeV
     * @return @f$\mu_{VBF,WW\to 2l2\nu}@f$
     */
    virtual const double muTHUVBFHWW2l2v(double sqrt_s) const;
    /**
     * @brief The ratio @f$\mu_{ZH,WW\to 2l2\nu}@f$ between the ZH
     * production cross-section with subsequent decay into @f$W W^*\to 2l2\nu@f$ in the
     * current model and in the Standard Model.
     * @param[in] sqrt_s the center-of-mass energy in TeV
     * @return @f$\mu_{ZH,WW\to 2l2\nu}@f$
     */
    virtual const double muTHUZHWW2l2v(double sqrt_s) const;
    /**
     * @brief The ratio @f$\mu_{WH,WW\to 2l2\nu}@f$ between the WH
     * production cross-section with subsequent decay into @f$W W^*\to 2l2\nu@f$ in the
     * current model and in the Standard Model.
     * @param[in] sqrt_s the center-of-mass energy in TeV
     * @return @f$\mu_{WH,WW\to 2l2\nu}@f$
     */
    virtual const double muTHUWHWW2l2v(double sqrt_s) const;
    /**
     * @brief The ratio @f$\mu_{VH,WW\to 2l2\nu}@f$ between the VH
     * production cross-section with subsequent decay into @f$W W^*\to 2l2\nu@f$ in the
     * current model and in the Standard Model.
     * @param[in] sqrt_s the center-of-mass energy in TeV
     * @return @f$\mu_{VH,WW\to 2l2\nu}@f$
     */
    virtual const double muTHUVHWW2l2v(double sqrt_s) const;
    /**
     * @brief The ratio @f$\mu_{ttH,WW\to 2l2\nu}@f$ between the ttH
     * production cross-section with subsequent decay into @f$W W^*\to 2l2\nu@f$ in the
     * current model and in the Standard Model.
     * @param[in] sqrt_s the center-of-mass energy in TeV
     * @return @f$\mu_{ttH,WW\to 2l2\nu}@f$
     */
    virtual const double muTHUttHWW2l2v(double sqrt_s) const;
    /**
     * @brief The ratio @f$\mu_{ggH,\mu\mu}@f$ between the gluon-gluon fusion Higgs
     * production cross-section with subsequent decay into @f$\mu\mu@f$ in the
     * current model and in the Standard Model.
     * @param[in] sqrt_s the center-of-mass energy in TeV
     * @return @f$\mu_{ggH,\mu\mu}@f$
     */
    virtual const double muTHUggHmumu(double sqrt_s) const;
    /**
     * @brief The ratio @f$\mu_{VBF,\mu\mu}@f$ between the VBF Higgs
     * production cross-section with subsequent decay into @f$\mu\mu@f$ in the
     * current model and in the Standard Model.
     * @param[in] sqrt_s the center-of-mass energy in TeV
     * @return @f$\mu_{VBF,\mu\mu}@f$
     */
    virtual const double muTHUVBFHmumu(double sqrt_s) const;
    /**
     * @brief The ratio @f$\mu_{ZH,\mu\mu}@f$ between the ZH
     * production cross-section with subsequent decay into @f$\mu\mu@f$ in the
     * current model and in the Standard Model.
     * @param[in] sqrt_s the center-of-mass energy in TeV
     * @return @f$\mu_{ZH,\mu\mu}@f$
     */
    virtual const double muTHUZHmumu(double sqrt_s) const;
    /**
     * @brief The ratio @f$\mu_{WH,\mu\mu}@f$ between the WH
     * production cross-section with subsequent decay into @f$\mu\mu@f$ in the
     * current model and in the Standard Model.
     * @param[in] sqrt_s the center-of-mass energy in TeV
     * @return @f$\mu_{WH,\mu\mu}@f$
     */
    virtual const double muTHUWHmumu(double sqrt_s) const;
    /**
     * @brief The ratio @f$\mu_{VH,\mu\mu}@f$ between the VH
     * production cross-section with subsequent decay into @f$\mu\mu@f$ in the
     * current model and in the Standard Model.
     * @param[in] sqrt_s the center-of-mass energy in TeV
     * @return @f$\mu_{VH,\mu\mu}@f$
     */
    virtual const double muTHUVHmumu(double sqrt_s) const;
    /**
     * @brief The ratio @f$\mu_{ttH,\mu\mu}@f$ between the ttH
     * production cross-section with subsequent decay into @f$\mu\mu@f$ in the
     * current model and in the Standard Model.
     * @param[in] sqrt_s the center-of-mass energy in TeV
     * @return @f$\mu_{ttH,\mu\mu}@f$
     */
    virtual const double muTHUttHmumu(double sqrt_s) const;
    /**
     * @brief The ratio @f$\mu_{ggH,\tau\tau}@f$ between the gluon-gluon fusion Higgs
     * production cross-section with subsequent decay into @f$\tau\tau@f$ in the
     * current model and in the Standard Model.
     * @param[in] sqrt_s the center-of-mass energy in TeV
     * @return @f$\mu_{ggH,\tau\tau}@f$
     */
    virtual const double muTHUggHtautau(double sqrt_s) const;
    /**
     * @brief The ratio @f$\mu_{VBF,\tau\tau}@f$ between the VBF Higgs
     * production cross-section with subsequent decay into @f$\tau\tau@f$ in the
     * current model and in the Standard Model.
     * @param[in] sqrt_s the center-of-mass energy in TeV
     * @return @f$\mu_{VBF,\tau\tau}@f$
     */
    virtual const double muTHUVBFHtautau(double sqrt_s) const;
    /**
     * @brief The ratio @f$\mu_{ZH,\tau\tau}@f$ between the ZH
     * production cross-section with subsequent decay into @f$\tau\tau@f$ in the
     * current model and in the Standard Model.
     * @param[in] sqrt_s the center-of-mass energy in TeV
     * @return @f$\mu_{ZH,\tau\tau}@f$
     */
    virtual const double muTHUZHtautau(double sqrt_s) const;
    /**
     * @brief The ratio @f$\mu_{WH,\tau\tau}@f$ between the WH
     * production cross-section with subsequent decay into @f$\tau\tau@f$ in the
     * current model and in the Standard Model.
     * @param[in] sqrt_s the center-of-mass energy in TeV
     * @return @f$\mu_{WH,\tau\tau}@f$
     */
    virtual const double muTHUWHtautau(double sqrt_s) const;
    /**
     * @brief The ratio @f$\mu_{VH,\tau\tau}@f$ between the VH
     * production cross-section with subsequent decay into @f$\tau\tau@f$ in the
     * current model and in the Standard Model.
     * @param[in] sqrt_s the center-of-mass energy in TeV
     * @return @f$\mu_{VH,\tau\tau}@f$
     */
    virtual const double muTHUVHtautau(double sqrt_s) const;
    /**
     * @brief The ratio @f$\mu_{ttH,\tau\tau}@f$ between the ttH
     * production cross-section with subsequent decay into @f$\tau\tau@f$ in the
     * current model and in the Standard Model.
     * @param[in] sqrt_s the center-of-mass energy in TeV
     * @return @f$\mu_{ttH,\tau\tau}@f$
     */
    virtual const double muTHUttHtautau(double sqrt_s) const;
    /**
     * @brief The ratio @f$\mu_{ggH,bb}@f$ between the gluon-gluon fusion Higgs
     * production cross-section with subsequent decay into @f$bb@f$ in the
     * current model and in the Standard Model.
     * @param[in] sqrt_s the center-of-mass energy in TeV
     * @return @f$\mu_{ggH,bb}@f$
     */
    virtual const double muTHUggHbb(double sqrt_s) const;
    /**
     * @brief The ratio @f$\mu_{VBF,bb}@f$ between the VBF Higgs
     * production cross-section with subsequent decay into @f$bb@f$ in the
     * current model and in the Standard Model.
     * @param[in] sqrt_s the center-of-mass energy in TeV
     * @return @f$\mu_{VBF,bb}@f$
     */
    virtual const double muTHUVBFHbb(double sqrt_s) const;
    /**
     * @brief The ratio @f$\mu_{ZH,bb}@f$ between the ZH
     * production cross-section with subsequent decay into @f$bb@f$ in the
     * current model and in the Standard Model.
     * @param[in] sqrt_s the center-of-mass energy in TeV
     * @return @f$\mu_{ZH,bb}@f$
     */
    virtual const double muTHUZHbb(double sqrt_s) const;
    /**
     * @brief The ratio @f$\mu_{WH,bb}@f$ between the WH
     * production cross-section with subsequent decay into @f$bb@f$ in the
     * current model and in the Standard Model.
     * @param[in] sqrt_s the center-of-mass energy in TeV
     * @return @f$\mu_{WH,bb}@f$
     */
    virtual const double muTHUWHbb(double sqrt_s) const;
    /**
     * @brief The ratio @f$\mu_{VH,bb}@f$ between the VH
     * production cross-section with subsequent decay into @f$bb@f$ in the
     * current model and in the Standard Model.
     * @param[in] sqrt_s the center-of-mass energy in TeV
     * @return @f$\mu_{VH,bb}@f$
     */
    virtual const double muTHUVHbb(double sqrt_s) const;
    /**
     * @brief The ratio @f$\mu_{ttH,bb}@f$ between the ttH
     * production cross-section with subsequent decay into @f$bb@f$ in the
     * current model and in the Standard Model.
     * @param[in] sqrt_s the center-of-mass energy in TeV
     * @return @f$\mu_{ttH,bb}@f$
     */
    virtual const double muTHUttHbb(double sqrt_s) const;

    /**
     * @brief The ratio @f$\mu_{VBF}@f$ between the VBF
     * production cross-section in the
     * current model and in the Standard Model, multiplied by the 
     * total (SM+new physics) invisible decay branching ratio.
     * @param[in] sqrt_s the center-of-mass energy in TeV
     * @return @f$\mu_{VBF}BR_{inv}@f$
     */
    virtual const double muTHUVBFBRinv(double sqrt_s) const;
    /**
     * @brief The ratio @f$\mu_{VBF,inv}@f$ between the VBF
     * production cross-section with subsequent decay into invisible states in the
     * current model and in the Standard Model.
     * @param[in] sqrt_s the center-of-mass energy in TeV
     * @return @f$\mu_{VBF,inv}@f$
     */
    virtual const double muTHUVBFHinv(double sqrt_s) const;
    /**
     * @brief The ratio @f$\mu_{VH}@f$ between the VH
     * production cross-section in the
     * current model and in the Standard Model, multiplied by the 
     * total (SM+new physics) invisible decay branching ratio.
     * @param[in] sqrt_s the center-of-mass energy in TeV
     * @return @f$\mu_{VH}BR_{inv}@f$
     */
    virtual const double muTHUVHBRinv(double sqrt_s) const;
    /**
     * @brief The ratio @f$\mu_{VH,inv}@f$ between the VH
     * production cross-section with subsequent decay into invisible states in the
     * current model and in the Standard Model.
     * @param[in] sqrt_s the center-of-mass energy in TeV
     * @return @f$\mu_{VH,inv}@f$
     */
    virtual const double muTHUVHinv(double sqrt_s) const;

    /**
     * @brief The ratio @f$\mu_{ggH,ZZ\to 4\mu}@f$ between the gluon-gluon fusion Higgs
     * production cross-section with subsequent decay into @f$Z Z^*\to 4\mu@f$ in the
     * current model and in the Standard Model.
     * @param[in] sqrt_s the center-of-mass energy in TeV
     * @return @f$\mu_{ggH,ZZ\to 4\mu}@f$
     */
    virtual const double muTHUggHZZ4mu(double sqrt_s) const;
    /**
     * @brief The ratio @f$\mu_{ggH,Z\gamma\to \gamma 2\mu}@f$ between the gluon-gluon fusion Higgs
     * production cross-section with subsequent decay into @f$Z \gamma\to \gamma 2\mu@f$ in the
     * current model and in the Standard Model.
     * @param[in] sqrt_s the center-of-mass energy in TeV
     * @return @f$\mu_{ggH,Z\gamma\to \gamma 2\mu}@f$
     */
    virtual const double muTHUggHZgamumu(double sqrt_s) const;

    ////////////////////////////////////////////////////////////////////////


    /**
     * @brief The new physics contribution to the anomalous triple gauge coupling @f$g_{1,Z}@f$.
     * @return @f$\delta g_{1,Z}@f$
     */
    virtual const double deltag1ZNP(const double mu) const;

    /**
     * @brief The new physics contribution to the anomalous triple gauge coupling @f$\kappa_{Z}@f$.
     * @return @f$\delta \kappa_{Z}@f$
     */
    virtual const double deltaKZNP(const double mu) const;

    /**
     * @brief The new physics contribution to the anomalous triple gauge coupling @f$g_{1,\gamma}@f$.
     * @return @f$\delta g_{1,\gamma}@f$
     */
    virtual const double deltag1gaNP(const double mu) const;

    /**
     * @brief The new physics contribution to the anomalous triple gauge coupling @f$\kappa_{\gamma}@f$.
     * @return @f$\delta \kappa_{\gamma}@f$
     */
    virtual const double deltaKgammaNP(const double mu) const;

    /**
     * @brief The new physics contribution to the anomalous triple gauge coupling @f$\lambda_{Z}@f$.
     * @return @f$\lambda_{Z}@f$
     */
    virtual const double lambdaZNP(const double mu) const;

    ////////////////////////////////////////////////////////////////////////

    /**
     * @brief The new physics contribution to the effective anomalous triple 
     * gauge coupling @f$g_{1,Z}^{Eff}@f$ from arXiv: 1708.09079 [hep-ph].
     * @return @f$\delta g_{1,Z}@f$
     */
    virtual const double deltag1ZNPEff() const;

    /**
     * @brief The new physics contribution to the effective anomalous triple 
     * gauge coupling @f$\kappa_{\gamma}^{Eff}@f$ from arXiv: 1708.09079 [hep-ph].
     * @return @f$\delta \kappa_{\gamma}@f$
     */
    virtual const double deltaKgammaNPEff() const;

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
    virtual const double deltaxseeWW4fLEP2(double sqrt_s, const int fstate) const;

    /**
     * @brief The cross section in pb for @f$e^+ e^- \to W^+ W^- \to 4f @f$, 
     * with @f$ 4f = 0 (jjjj), 1 (e v jj), 2 (mu v jj), 3 (tau v jj),
     * 4 (e v e v), 5 (mu v mu v), 6 (tau v tau v),
     * 7 (e v mu v), 8 (e v tau v), 9 (mu v tau v), 10 (l v jj), 11 (l v l v) @f$
     * the different fermion final states for C.O.M. energies in 188-208 GeV.
     * From arXiv: 1606.06693 [hep-ph].
     * @return @f$sigma@f$ [pb]
     */
    virtual const double xseeWW4fLEP2(double sqrt_s, const int fstate) const;

    /**
     * @brief The new physics contribution to the total cross section in pb for @f$e^+ e^- \to W^+ W^-@f$, 
     * summing over all final states for C.O.M. energies in 188-208 GeV.
     * From arXiv: 1606.06693 [hep-ph].
     * @return @f$\delta sigma@f$ [pb]
     */
    virtual const double deltaxseeWWtotLEP2(double sqrt_s) const;


    /**
     * @brief The total cross section in pb for @f$e^+ e^- \to W^+ W^-@f$, 
     * summing over all final states for C.O.M. energies in 188-208 GeV.
     * From arXiv: 1606.06693 [hep-ph].
     * @return @f$sigma@f$ [pb]
     */
    virtual const double xseeWWtotLEP2(double sqrt_s) const;

    /**
     * @brief The new physics contribution to the differential cross section in pb for @f$e^+ e^- \to W^+ W^- \to lv jj @f$, 
     * with @f$ l= e,\mu @f$ for the 4 @f$ cos{\theta}@f$ bins defined in arXiv: 1606.06693 [hep-ph].
     * for the C.O.M. energies of 182.6 and 205.9 GeV.
     * From arXiv: 1606.06693 [hep-ph].
     * @return @f$\delta d\sigma/d\cos{\theta}@f$ [pb]
     */
    virtual const double deltadxsdcoseeWWlvjjLEP2(double sqrt_s, const int bin) const;

    /**
     * @brief The differential cross section in pb for @f$e^+ e^- \to W^+ W^- \to lv jj @f$, 
     * with @f$ l= e,\mu @f$ for the 4 @f$ cos{\theta}@f$ bins defined in arXiv: 1606.06693 [hep-ph].
     * for the C.O.M. energies of 182.6 and 205.9 GeV.
     * From arXiv: 1606.06693 [hep-ph].
     * @return @f$d\sigma/d\cos{\theta}@f$ [pb]
     */
    virtual const double dxsdcoseeWWlvjjLEP2(double sqrt_s, const int bin) const;

    /**
     * @brief The differential distribution for @f$e^+ e^- \to W^+ W^- \to jj \ell \nu@f$, 
     * with @f$\ell= e, \mu@f$, as a function of the @f$W@f$ polar angle.
     * @return @f$d\sigma/d\cos{\theta}@f$
     */
    virtual const double dxseeWWdcos(double sqrt_s, double cos) const;

    /**
     * @brief The integral of differential distribution for @f$e^+ e^- \to W^+ W^- \to jj \ell \nu@f$, 
     * with @f$\ell= e, \mu@f$ in a given bin of the @f$W@f$ polar angle.
     * @return @f$\int_{\cos{\theta_1}}^{\cos{\theta_2}} d\sigma/d\cos{\theta}@f$
     */
    virtual const double dxseeWWdcosBin(double sqrt_s, double cos1, double cos2) const;

    /**
     * @brief Total @f$e^+ e^- \to W^+ W^- \to jj \ell \nu@f$ cross section in pb, 
     * with @f$\ell= e, \mu@f$.
     * @return @f$\sigma(e^+ e^- \to W^+ W^- \to jj \ell \nu) @f$
     */
    virtual const double xseeWW(double sqrt_s) const;

    /**
     * @brief The ratio @f$\mu_{eeWW}@f$ between the 
     * @f$ e^{+}e^{-}\to W^{+}W^{-} @f$ production
     * cross-section in the current model and in the Standard Model.
     * @param[in] sqrt_s the center-of-mass energy in TeV, Pol_em and Pol_ep
     * are the polarization of electrons and positrons, respectively 
     * @return @f$\mu_{eeWW}@f$
     */
    virtual const double mueeWW(double sqrt_s, double Pol_em, double Pol_ep) const;

    /**
     * @brief The ratio @f$\mu_{eeWW}@f$ between the 
     * @f$ e^{+}e^{-}\to W^{+}W^{-} @f$ production
     * cross-section in the current model and in the Standard Model.
     * @param[in] sqrt_s the center-of-mass energy in TeV, Pol_em and Pol_ep
     * are the polarization of electrons and positrons, respectively 
     * @return @f$\mu_{eeWW}@f$
     */
    virtual const double mueeWWPol(double sqrt_s, double Pol_em, double Pol_ep) const;

    ////////////////////////////////////////////////////////////////////////

    //----- High Energy diboson observables at hadron colliders

    /**
     * @brief The direction constrained by @f$ p p \to Z H@f$ in the boosted regime, @f$g_p^Z@f$.
     * From arXiv:1807.01796 and the contribution to FCC CDR Vol 1. Implemented only in NPSMEFTd6 class.
     * @return @f$g_p^Z@f$
     */
    virtual const double ppZHprobe(double sqrt_s) const;

    /**
     * @brief The number of events in  @f$ p p \to WZ@f$
     * in a given @f$p_{TV}@f$ bin, normalized to the SM prediction.
     * From arXiv: 1712.01310 [hep-ph] and private communication.
     * Implemented only in NPSMEFTd6 class.
     * @return @f$N_{ev}^{p_{TV}}/N_{ev,SM}^{p_{TV}}@f$
     */
    virtual const double mupTVppWZ(double sqrt_s, double pTV1, double pTV2) const;


    ////////////////////////////////////////////////////////////////////////

    //----- Simplified Template Cross Sections Bins

    //----- Stage 0

    /**
     * @brief The STXS0 bin @f$pp \to H qq@f$.
     * @param[in] sqrt_s the center-of-mass energy in TeV
     */
    virtual const double STXS0_qqH(double sqrt_s) const;


    //----- Stage 1

    /**
     * @brief The STXS bin @f$gg \to H@f$.
     * @param[in] sqrt_s the center-of-mass energy in TeV
     */
    virtual const double STXS_ggH_VBFtopo_j3v(double sqrt_s) const;


    /**
     * @brief The STXS bin @f$gg \to H@f$.
     * @param[in] sqrt_s the center-of-mass energy in TeV
     */
    virtual const double STXS_ggH_VBFtopo_j3(double sqrt_s) const;

    /**
     * @brief The STXS bin @f$gg \to H@f$.
     * @param[in] sqrt_s the center-of-mass energy in TeV
     */
    virtual const double STXS_ggH0j(double sqrt_s) const;

    /**
     * @brief The STXS bin @f$gg \to H@f$.
     * @param[in] sqrt_s the center-of-mass energy in TeV
     */
    virtual const double STXS_ggH1j_pTH_0_60(double sqrt_s) const;

    /**
     * @brief The STXS bin @f$gg \to H@f$.
     * @param[in] sqrt_s the center-of-mass energy in TeV
     */
    virtual const double STXS_ggH1j_pTH_60_120(double sqrt_s) const;

    /**
     * @brief The STXS bin @f$gg \to H@f$.
     * @param[in] sqrt_s the center-of-mass energy in TeV
     */
    virtual const double STXS_ggH1j_pTH_120_200(double sqrt_s) const;

    /**
     * @brief The STXS bin @f$gg \to H@f$.
     * @param[in] sqrt_s the center-of-mass energy in TeV
     */
    virtual const double STXS_ggH1j_pTH_200(double sqrt_s) const;

    /**
     * @brief The STXS bin @f$gg \to H@f$.
     * @param[in] sqrt_s the center-of-mass energy in TeV
     */
    virtual const double STXS_ggH2j_pTH_0_200(double sqrt_s) const;

    /**
     * @brief The STXS bin @f$gg \to H@f$.
     * @param[in] sqrt_s the center-of-mass energy in TeV
     */
    virtual const double STXS_ggH2j_pTH_0_60(double sqrt_s) const;


    /**
     * @brief The STXS bin @f$gg \to H@f$.
     * @param[in] sqrt_s the center-of-mass energy in TeV
     */
    virtual const double STXS_ggH2j_pTH_60_120(double sqrt_s) const;


    /**
     * @brief The STXS bin @f$gg \to H@f$.
     * @param[in] sqrt_s the center-of-mass energy in TeV
     */
    virtual const double STXS_ggH2j_pTH_120_200(double sqrt_s) const;


    /**
     * @brief The STXS bin @f$gg \to H@f$.
     * @param[in] sqrt_s the center-of-mass energy in TeV
     */
    virtual const double STXS_ggH2j_pTH_200(double sqrt_s) const;

    /**
     * @brief The STXS bin @f$qq \to H qq@f$.
     * @param[in] sqrt_s the center-of-mass energy in TeV
     */
    virtual const double STXS_qqHqq_VBFtopo_Rest(double sqrt_s) const;


    /**
     * @brief The STXS bin @f$qq \to H qq@f$.
     * @param[in] sqrt_s the center-of-mass energy in TeV
     */
    virtual const double STXS_qqHqq_VBFtopo_j3v(double sqrt_s) const;


    /**
     * @brief The STXS bin @f$qq \to H qq@f$.
     * @param[in] sqrt_s the center-of-mass energy in TeV
     */
    virtual const double STXS_qqHqq_VBFtopo_j3(double sqrt_s) const;

    /**
     * @brief The STXS bin @f$qq \to H qq@f$.
     * @param[in] sqrt_s the center-of-mass energy in TeV
     */
    virtual const double STXS_qqHqq_nonVHtopo(double sqrt_s) const;

    /**
     * @brief The STXS bin @f$qq \to H qq@f$.
     * @param[in] sqrt_s the center-of-mass energy in TeV
     */
    virtual const double STXS_qqHqq_VHtopo(double sqrt_s) const;

    /**
     * @brief The STXS bin @f$qq \to H qq@f$.
     * @param[in] sqrt_s the center-of-mass energy in TeV
     */
    virtual const double STXS_qqHqq_Rest(double sqrt_s) const;

    /**
     * @brief The STXS bin @f$qq \to H qq@f$.
     * @param[in] sqrt_s the center-of-mass energy in TeV
     */
    virtual const double STXS_qqHqq_pTj_200(double sqrt_s) const;

    /**
     * @brief The STXS bin @f$qq \to H \ell \nu@f$.
     * @param[in] sqrt_s the center-of-mass energy in TeV
     */
    virtual const double STXS_qqHlv_pTV_0_250(double sqrt_s) const;


    /**
     * @brief The STXS bin @f$qq \to H \ell \nu@f$.
     * @param[in] sqrt_s the center-of-mass energy in TeV
     */
    virtual const double STXS_qqHlv_pTV_0_150(double sqrt_s) const;


    /**
     * @brief The STXS bin @f$qq \to H \ell \nu@f$.
     * @param[in] sqrt_s the center-of-mass energy in TeV
     */
    virtual const double STXS_qqHlv_pTV_150_250_0j(double sqrt_s) const;


    /**
     * @brief The STXS bin @f$qq \to H \ell \nu@f$.
     * @param[in] sqrt_s the center-of-mass energy in TeV
     */
    virtual const double STXS_qqHlv_pTV_150_250_1j(double sqrt_s) const;


    /**
     * @brief The STXS bin @f$qq \to H \ell \nu@f$.
     * @param[in] sqrt_s the center-of-mass energy in TeV
     */
    virtual const double STXS_qqHlv_pTV_250(double sqrt_s) const;

    /**
     * @brief The STXS bin @f$qq \to H \ell \ell@f$.
     * @param[in] sqrt_s the center-of-mass energy in TeV
     */
    virtual const double STXS_qqHll_pTV_0_150(double sqrt_s) const;

    /**
     * @brief The STXS bin @f$qq \to H \ell \ell@f$.
     * @param[in] sqrt_s the center-of-mass energy in TeV
     */
    virtual const double STXS_qqHll_pTV_150_250(double sqrt_s) const;


    /**
     * @brief The STXS bin @f$qq \to H \ell \ell@f$.
     * @param[in] sqrt_s the center-of-mass energy in TeV
     */
    virtual const double STXS_qqHll_pTV_150_250_0j(double sqrt_s) const;


    /**
     * @brief The STXS bin @f$qq \to H \ell \ell@f$.
     * @param[in] sqrt_s the center-of-mass energy in TeV
     */
    virtual const double STXS_qqHll_pTV_150_250_1j(double sqrt_s) const;

    /**
     * @brief The STXS bin @f$qq \to H \ell \ell@f$.
     * @param[in] sqrt_s the center-of-mass energy in TeV
     */
    virtual const double STXS_qqHll_pTV_250(double sqrt_s) const;


    /**
     * @brief The STXS bin @f$ ttH + tH @f$.
     * @param[in] sqrt_s the center-of-mass energy in TeV
     */
    virtual const double STXS_ttHtH(double sqrt_s) const;

    /**
     * @brief The STXS bin @f$ qq \to WH \to H qq @f$.
     * @param[in] sqrt_s the center-of-mass energy in TeV
     */
    virtual const double STXS_WHqqHqq_VBFtopo_j3v(double sqrt_s) const;

    /**
     * @brief The STXS bin @f$ qq \to WH \to H qq @f$.
     * @param[in] sqrt_s the center-of-mass energy in TeV
     */
    virtual const double STXS_WHqqHqq_VBFtopo_j3(double sqrt_s) const;

    /**
     * @brief The STXS bin @f$ qq \to WH \to H qq @f$.
     * @param[in] sqrt_s the center-of-mass energy in TeV
     */
    virtual const double STXS_WHqqHqq_VH2j(double sqrt_s) const;

    /**
     * @brief The STXS bin @f$ qq \to WH \to H qq @f$.
     * @param[in] sqrt_s the center-of-mass energy in TeV
     */
    virtual const double STXS_WHqqHqq_Rest(double sqrt_s) const;

    /**
     * @brief The STXS bin @f$ qq \to WH \to H qq @f$.
     * @param[in] sqrt_s the center-of-mass energy in TeV
     */
    virtual const double STXS_WHqqHqq_pTj1_200(double sqrt_s) const;

    /**
     * @brief The STXS bin @f$ qq \to ZH \to H qq @f$.
     * @param[in] sqrt_s the center-of-mass energy in TeV
     */
    virtual const double STXS_ZHqqHqq_VBFtopo_j3v(double sqrt_s) const;

    /**
     * @brief The STXS bin @f$ qq \to ZH \to H qq @f$.
     * @param[in] sqrt_s the center-of-mass energy in TeV
     */
    virtual const double STXS_ZHqqHqq_VBFtopo_j3(double sqrt_s) const;

    /**
     * @brief The STXS bin @f$ qq \to ZH \to H qq @f$.
     * @param[in] sqrt_s the center-of-mass energy in TeV
     */
    virtual const double STXS_ZHqqHqq_VH2j(double sqrt_s) const;

    /**
     * @brief The STXS bin @f$ qq \to ZH \to H qq @f$.
     * @param[in] sqrt_s the center-of-mass energy in TeV
     */
    virtual const double STXS_ZHqqHqq_Rest(double sqrt_s) const;

    /**
     * @brief The STXS bin @f$ qq \to ZH \to H qq @f$.
     * @param[in] sqrt_s the center-of-mass energy in TeV
     */
    virtual const double STXS_ZHqqHqq_pTj1_200(double sqrt_s) const;


    //----- Stage 1.2 
    // From ATLAS-CONF-2020-053
    // Expressions valid in the {G_F, M_Z, M_W} scheme

    /**
     * @brief The STXS BR @f$ H \to 4l @f$, @f$l=e,\mu@f$.
     */
    virtual const double STXS12_BrH4lRatio() const;

    /**
     * @brief The STXS BR @f$ H \to e\nu \mu\nu @f$.
     */
    virtual const double STXS12_BrHevmuvRatio() const;

    /**
     * @brief The STXS BR @f$ H \to \gamma \gamma @f$.
     */
    virtual const double STXS12_BrHgagaRatio() const;

    /**
     * @brief The STXS BR @f$ H \to bb @f$.
     */
    virtual const double STXS12_BrHbbRatio() const;

    /**
     * @brief The STXS bin @f$gg \to H@f$, @f$N_j\leq 1,~200<p_{TH} [GeV]<300@f$.
     * @param[in] sqrt_s the center-of-mass energy in TeV
     */
    virtual const double STXS12_ggH_pTH200_300_Nj01(double sqrt_s) const;

    /**
     * @brief The STXS bin @f$gg \to H@f$, @f$N_j\leq 1,~300<p_{TH} [GeV]<450@f$.
     * @param[in] sqrt_s the center-of-mass energy in TeV
     */
    virtual const double STXS12_ggH_pTH300_450_Nj01(double sqrt_s) const;

    /**
     * @brief The STXS bin @f$gg \to H@f$, @f$N_j\leq 1,~450<p_{TH} [GeV]<650@f$.
     * @param[in] sqrt_s the center-of-mass energy in TeV
     */
    virtual const double STXS12_ggH_pTH450_650_Nj01(double sqrt_s) const;

    /**
     * @brief The STXS bin @f$gg \to H@f$, @f$N_j\leq 1,650<p_{TH} [GeV]@f$.
     * @param[in] sqrt_s the center-of-mass energy in TeV
     */
    virtual const double STXS12_ggH_pTH650_Inf_Nj01(double sqrt_s) const;

    /**
     * @brief The STXS bin @f$gg \to H@f$, @f$N_j=0,~p_{TH} [GeV]<10@f$.
     * @param[in] sqrt_s the center-of-mass energy in TeV
     */
    virtual const double STXS12_ggH_pTH0_10_Nj0(double sqrt_s) const;

    /**
     * @brief The STXS bin @f$gg \to H@f$, @f$N_j=0,~10<p_{TH} [GeV]@f$.
     * @param[in] sqrt_s the center-of-mass energy in TeV
     */
    virtual const double STXS12_ggH_pTH10_Inf_Nj0(double sqrt_s) const;

    /**
     * @brief The STXS bin @f$gg \to H@f$, @f$N_j = 1,~p_{TH} [GeV]<60@f$.
     * @param[in] sqrt_s the center-of-mass energy in TeV
     */
    virtual const double STXS12_ggH_pTH0_60_Nj1(double sqrt_s) const;

    /**
     * @brief The STXS bin @f$gg \to H@f$, @f$N_j = 1,~60<p_{TH} [GeV]<120@f$.
     * @param[in] sqrt_s the center-of-mass energy in TeV
     */
    virtual const double STXS12_ggH_pTH60_120_Nj1(double sqrt_s) const;

    /**
     * @brief The STXS bin @f$gg \to H@f$, @f$N_j = 1,~120<p_{TH} [GeV]<200@f$.
     * @param[in] sqrt_s the center-of-mass energy in TeV
     */
    virtual const double STXS12_ggH_pTH120_200_Nj1(double sqrt_s) const;

    /**
     * @brief The STXS bin @f$gg \to H@f$, @f$N_j \geq 2,~m_{jj}[GeV]<350,~p_{TH} [GeV]<60@f$.
     * @param[in] sqrt_s the center-of-mass energy in TeV
     */
    virtual const double STXS12_ggH_mjj0_350_pTH0_60_Nj2(double sqrt_s) const;

    /**
     * @brief The STXS bin @f$gg \to H@f$, @f$N_j \geq 2,~m_{jj}[GeV]<350,~60<p_{TH} [GeV]<120@f$.
     * @param[in] sqrt_s the center-of-mass energy in TeV
     */
    virtual const double STXS12_ggH_mjj0_350_pTH60_120_Nj2(double sqrt_s) const;

    /**
     * @brief The STXS bin @f$gg \to H@f$, @f$N_j \geq 2,~m_{jj}[GeV]<350,~120<p_{TH} [GeV]<200@f$.
     * @param[in] sqrt_s the center-of-mass energy in TeV
     */
    virtual const double STXS12_ggH_mjj0_350_pTH120_200_Nj2(double sqrt_s) const;

    /**
     * @brief The STXS bin @f$gg \to H@f$, @f$N_j \geq 2,~350<m_{jj}[GeV]<700,~p_{TH} [GeV]<200,~p_{THjj}[GeV]<25@f$.
     * @param[in] sqrt_s the center-of-mass energy in TeV
     */
    virtual const double STXS12_ggH_mjj350_700_pTH0_200_ptHjj0_25_Nj2(double sqrt_s) const;

    /**
     * @brief The STXS bin @f$gg \to H@f$, @f$N_j \geq 2,~350<m_{jj}[GeV]<700,~p_{TH} [GeV]<200,~25<p_{THjj}[GeV]@f$.
     * @param[in] sqrt_s the center-of-mass energy in TeV
     */
    virtual const double STXS12_ggH_mjj350_700_pTH0_200_ptHjj25_Inf_Nj2(double sqrt_s) const;

    /**
     * @brief The STXS bin @f$gg \to H@f$, @f$N_j \geq 2,~700<m_{jj}[GeV],~p_{TH} [GeV]<200,~p_{THjj}[GeV]<25@f$.
     * @param[in] sqrt_s the center-of-mass energy in TeV
     */
    virtual const double STXS12_ggH_mjj700_Inf_pTH0_200_ptHjj0_25_Nj2(double sqrt_s) const;

    /**
     * @brief The STXS bin @f$gg \to H@f$, @f$N_j \geq 2,~700<m_{jj}[GeV],~p_{TH} [GeV]<200,~25<p_{THjj}[GeV]@f$.
     * @param[in] sqrt_s the center-of-mass energy in TeV
     */
    virtual const double STXS12_ggH_mjj700_Inf_pTH0_200_ptHjj25_Inf_Nj2(double sqrt_s) const;

    /**
     * @brief The STXS bin @f$gg \to H\ell\ell@f$, @f$p_{TV}[GeV]<75@f$.
     * @param[in] sqrt_s the center-of-mass energy in TeV
     */
    virtual const double STXS12_ggHll_pTV0_75(double sqrt_s) const;

    /**
     * @brief The STXS bin @f$gg \to H\ell\ell@f$, @f$75<p_{TV}[GeV]<150@f$.
     * @param[in] sqrt_s the center-of-mass energy in TeV
     */
    virtual const double STXS12_ggHll_pTV75_150(double sqrt_s) const;

    /**
     * @brief The STXS bin @f$gg \to H\ell\ell@f$, @f$N_j = 0,~150<p_{TV}[GeV]<250@f$.
     * @param[in] sqrt_s the center-of-mass energy in TeV
     */
    virtual const double STXS12_ggHll_pTV150_250_Nj0(double sqrt_s) const;

    /**
     * @brief The STXS bin @f$gg \to H\ell\ell@f$, @f$N_j = 1,~150<p_{TV}[GeV]<250@f$.
     * @param[in] sqrt_s the center-of-mass energy in TeV
     */
    virtual const double STXS12_ggHll_pTV150_250_Nj1(double sqrt_s) const;

    /**
     * @brief The STXS bin @f$gg \to H\ell\ell@f$, @f$250 < p_{TV}[GeV]@f$.
     * @param[in] sqrt_s the center-of-mass energy in TeV
     */
    virtual const double STXS12_ggHll_pTV250_Inf(double sqrt_s) const;


    /**
     * @brief The STXS bin @f$qq \to Hqq@f$, @f$N_j = 0@f$.
     * @param[in] sqrt_s the center-of-mass energy in TeV
     */
    virtual const double STXS12_qqHqq_Nj0(double sqrt_s) const;

    /**
     * @brief The STXS bin @f$qq \to Hqq@f$, @f$N_j = 1@f$.
     * @param[in] sqrt_s the center-of-mass energy in TeV
     */
    virtual const double STXS12_qqHqq_Nj1(double sqrt_s) const;

    /**
     * @brief The STXS bin @f$qq \to Hqq@f$, @f$N_j \geq 2,~m_{jj}[GeV]<60@f$.
     * @param[in] sqrt_s the center-of-mass energy in TeV
     */
    virtual const double STXS12_qqHqq_mjj0_60_Nj2(double sqrt_s) const;

    /**
     * @brief The STXS bin @f$qq \to Hqq@f$, @f$N_j \geq 2,~60<m_{jj}[GeV]<120@f$.
     * @param[in] sqrt_s the center-of-mass energy in TeV
     */
    virtual const double STXS12_qqHqq_mjj60_120_Nj2(double sqrt_s) const;

    /**
     * @brief The STXS bin @f$qq \to Hqq@f$, @f$N_j \geq 2,~120<m_{jj}[GeV]<350@f$.
     * @param[in] sqrt_s the center-of-mass energy in TeV
     */
    virtual const double STXS12_qqHqq_mjj120_350_Nj2(double sqrt_s) const;

    /**
     * @brief The STXS bin @f$qq \to Hqq@f$, @f$N_j \geq 2,~350<m_{jj}[GeV],~200<p_{TH}[GeV]@f$.
     * @param[in] sqrt_s the center-of-mass energy in TeV
     */
    virtual const double STXS12_qqHqq_mjj350_Inf_pTH200_Inf_Nj2(double sqrt_s) const;

    /**
     * @brief The STXS bin @f$qq \to Hqq@f$, @f$N_j \geq 2,~350<m_{jj}[GeV]<700,~p_{TH}[GeV]<200,~p_{THjj}[GeV]<25@f$.
     * @param[in] sqrt_s the center-of-mass energy in TeV
     */
    virtual const double STXS12_qqHqq_mjj350_700_pTH0_200_pTHjj0_25_Nj2(double sqrt_s) const;

    /**
     * @brief The STXS bin @f$qq \to Hqq@f$, @f$N_j \geq 2,~350<m_{jj}[GeV]<700,~p_{TH}[GeV]<200,~25<p_{THjj}[GeV]@f$.
     * @param[in] sqrt_s the center-of-mass energy in TeV
     */
    virtual const double STXS12_qqHqq_mjj350_700_pTH0_200_pTHjj25_Inf_Nj2(double sqrt_s) const;

    /**
     * @brief The STXS bin @f$qq \to Hqq@f$, @f$N_j \geq 2,~700<m_{jj}[GeV],~p_{TH}[GeV]<200,~p_{THjj}[GeV]<25@f$.
     * @param[in] sqrt_s the center-of-mass energy in TeV
     */
    virtual const double STXS12_qqHqq_mjj700_Inf_pTH0_200_pTHjj0_25_Nj2(double sqrt_s) const;

    /**
     * @brief The STXS bin @f$qq \to Hqq@f$, @f$N_j \geq 2,~700<m_{jj}[GeV],~p_{TH}[GeV]<200,~25<p_{THjj}[GeV]@f$.
     * @param[in] sqrt_s the center-of-mass energy in TeV
     */
    virtual const double STXS12_qqHqq_mjj700_Inf_pTH0_200_pTHjj25_Inf_Nj2(double sqrt_s) const;

    /**
     * @brief The STXS bin @f$qq \to H\ell\nu@f$, @f$p_{TV}[GeV]<75@f$.
     * @param[in] sqrt_s the center-of-mass energy in TeV
     */
    virtual const double STXS12_qqHlv_pTV0_75(double sqrt_s) const;

    /**
     * @brief The STXS bin @f$qq \to H\ell\nu@f$, @f$75<p_{TV}[GeV]<150@f$.
     * @param[in] sqrt_s the center-of-mass energy in TeV
     */
    virtual const double STXS12_qqHlv_pTV75_150(double sqrt_s) const;

    /**
     * @brief The STXS bin @f$qq \to H\ell\nu@f$, @f$N_j = 0,~150<p_{TV}[GeV]<250@f$.
     * @param[in] sqrt_s the center-of-mass energy in TeV
     */
    virtual const double STXS12_qqHlv_pTV150_250_Nj0(double sqrt_s) const;

    /**
     * @brief The STXS bin @f$qq \to H\ell\nu@f$, @f$N_j \geq 1,~150<p_{TV}[GeV]<250@f$.
     * @param[in] sqrt_s the center-of-mass energy in TeV
     */
    virtual const double STXS12_qqHlv_pTV150_250_Nj1(double sqrt_s) const;

    /**
     * @brief The STXS bin @f$qq \to H\ell\nu@f$, @f$250<p_{TV}[GeV]@f$.
     * @param[in] sqrt_s the center-of-mass energy in TeV
     */
    virtual const double STXS12_qqHlv_pTV250_Inf(double sqrt_s) const;

    /**
     * @brief The STXS bin @f$qq \to H\ell\ell@f$, @f$p_{TV}[GeV]<75@f$.
     * @param[in] sqrt_s the center-of-mass energy in TeV
     */
    virtual const double STXS12_qqHll_pTV0_75(double sqrt_s) const;

    /**
     * @brief The STXS bin @f$qq \to H\ell\ell@f$, @f$75<p_{TV}[GeV]<150@f$.
     * @param[in] sqrt_s the center-of-mass energy in TeV
     */
    virtual const double STXS12_qqHll_pTV75_150(double sqrt_s) const;

    /**
     * @brief The STXS bin @f$qq \to H\ell\ell@f$, @f$N_j = 0,~150<p_{TV}[GeV]<250@f$.
     * @param[in] sqrt_s the center-of-mass energy in TeV
     */
    virtual const double STXS12_qqHll_pTV150_250_Nj0(double sqrt_s) const;

    /**
     * @brief The STXS bin @f$qq \to H\ell\ell@f$, @f$N_j \geq 1,~150<p_{TV}[GeV]<250@f$.
     * @param[in] sqrt_s the center-of-mass energy in TeV
     */
    virtual const double STXS12_qqHll_pTV150_250_Nj1(double sqrt_s) const;

    /**
     * @brief The STXS bin @f$qq \to H\ell\ell@f$, @f$250<p_{TV}[GeV]@f$.
     * @param[in] sqrt_s the center-of-mass energy in TeV
     */
    virtual const double STXS12_qqHll_pTV250_Inf(double sqrt_s) const;

    /**
     * @brief The STXS bin @f$pp \to ttH@f$, @f$p_{TH}[GeV]<60@f$.
     * @param[in] sqrt_s the center-of-mass energy in TeV
     */
    virtual const double STXS12_ttH_pTH0_60(double sqrt_s) const;

    /**
     * @brief The STXS bin @f$pp \to ttH@f$, @f$60<p_{TH}[GeV]<120@f$.
     * @param[in] sqrt_s the center-of-mass energy in TeV
     */
    virtual const double STXS12_ttH_pTH60_120(double sqrt_s) const;

    /**
     * @brief The STXS bin @f$pp \to ttH@f$, @f$120<p_{TH}[GeV]<200@f$.
     * @param[in] sqrt_s the center-of-mass energy in TeV
     */
    virtual const double STXS12_ttH_pTH120_200(double sqrt_s) const;

    /**
     * @brief The STXS bin @f$pp \to ttH@f$, @f$200<p_{TH}[GeV]<300@f$.
     * @param[in] sqrt_s the center-of-mass energy in TeV
     */
    virtual const double STXS12_ttH_pTH200_300(double sqrt_s) const;

    /**
     * @brief The STXS bin @f$pp \to ttH@f$, @f$300<p_{TH}[GeV]@f$.
     * @param[in] sqrt_s the center-of-mass energy in TeV
     */
    virtual const double STXS12_ttH_pTH300_Inf(double sqrt_s) const;

    /**
     * @brief The STXS bin @f$pp \to tH@f$.
     * @param[in] sqrt_s the center-of-mass energy in TeV
     */
    virtual const double STXS12_tH(double sqrt_s) const;

    ////////////////////////////////////////////////////////////////////////

    /**
     * @brief The effective coupling @f$\kappa_{\mu,eff}=\sqrt{\Gamma_{H\mu\mu}/\Gamma_{H\mu\mu}^{SM}}@f$.
     * @return @f$\kappa_{\mu,eff}@f$
     */
    virtual const double kappamueff() const;

    /**
     * @brief The effective coupling @f$\kappa_{\tau,eff}=\sqrt{\Gamma_{H\tau\tau}/\Gamma_{H\tau\tau}^{SM}}@f$.
     * @return @f$\kappa_{\tau,eff}@f$
     */
    virtual const double kappataueff() const;

    /**
     * @brief The effective coupling @f$\kappa_{c,eff}=\sqrt{\Gamma_{Hcc}/\Gamma_{Hcc}^{SM}}@f$.
     * @return @f$\kappa_{c,eff}@f$
     */
    virtual const double kappaceff() const;

    /**
     * @brief The effective coupling @f$\kappa_{b,eff}=\sqrt{\Gamma_{Hbb}/\Gamma_{Hbb}^{SM}}@f$.
     * @return @f$\kappa_{b,eff}@f$
     */
    virtual const double kappabeff() const;

    /**
     * @brief The effective coupling @f$\kappa_{G,eff}=\sqrt{\Gamma_{HGG}/\Gamma_{HGG}^{SM}}@f$.
     * @return @f$\kappa_{G,eff}@f$
     */
    virtual const double kappaGeff() const;

    /**
     * @brief The effective coupling @f$\kappa_{Z,eff}=\sqrt{\Gamma_{HZZ}/\Gamma_{HZZ}^{SM}}@f$.
     * @return @f$\kappa_{Z,eff}@f$
     */
    virtual const double kappaZeff() const;

    /**
     * @brief The effective coupling @f$\kappa_{W,eff}=\sqrt{\Gamma_{HWW}/\Gamma_{HWW}^{SM}}@f$.
     * @return @f$\kappa_{W,eff}@f$
     */
    virtual const double kappaWeff() const;

    /**
     * @brief The effective coupling @f$\kappa_{A,eff}=\sqrt{\Gamma_{HAA}/\Gamma_{HAA}^{SM}}@f$.
     * @return @f$\kappa_{A,eff}@f$
     */
    virtual const double kappaAeff() const;

    /**
     * @brief The effective coupling @f$\kappa_{ZA,eff}=\sqrt{\Gamma_{HZA}/\Gamma_{HZA}^{SM}}@f$.
     * @return @f$\kappa_{ZA,eff}@f$
     */
    virtual const double kappaZAeff() const;

    /////////////Basic interactions of the so-called Higgs basis////////////////

    /**
     * @brief The Higgs-basis coupling @f$\delta y_t@f$.
     * (See LHCHXSWG-INT-2015-001 document.)
     * Note that the Lagrangian definition of the Higgs-basis parameters coincides with the one of 
     * some of the @f$g_i, \delta g_i@f$ couplings defined above.
     * In the Higgs basis, however, one uses the freedom to perform certain field redefinitions and 
     * operations to demand that the mass eigenstate Lagrangian has specific features. (See pag. 5,6 in the reference.)
     * Therefore, the actual expression in terms of dim 6 coefficients may differ from the one for @f$g_i, \delta g_i@f$.
     * @param[in] mu the RG scale to be used in the evaluation
     * @return @f$\delta y_t@f$
     */
    virtual const double deltayt_HB(const double mu) const;

    /**
     * @brief The Higgs-basis coupling @f$\delta y_b@f$.
     * (See LHCHXSWG-INT-2015-001 document.)
     * Note that the Lagrangian definition of the Higgs-basis parameters coincides with the one of 
     * some of the @f$g_i, \delta g_i@f$ couplings defined above.
     * In the Higgs basis, however, one uses the freedom to perform certain field redefinitions and 
     * operations to demand that the mass eigenstate Lagrangian has specific features. (See pag. 5,6 in the reference.)
     * Therefore, the actual expression in terms of dim 6 coefficients may differ from the one for @f$g_i, \delta g_i@f$.
     * @param[in] mu the RG scale to be used in the evaluation
     * @return @f$\delta y_b@f$
     */
    virtual const double deltayb_HB(const double mu) const;

    /**
     * @brief The Higgs-basis coupling @f$\delta y_\tau@f$.
     * (See LHCHXSWG-INT-2015-001 document.)
     * Note that the Lagrangian definition of the Higgs-basis parameters coincides with the one of 
     * some of the @f$g_i, \delta g_i@f$ couplings defined above.
     * In the Higgs basis, however, one uses the freedom to perform certain field redefinitions and 
     * operations to demand that the mass eigenstate Lagrangian has specific features. (See pag. 5,6 in the reference.)
     * Therefore, the actual expression in terms of dim 6 coefficients may differ from the one for @f$g_i, \delta g_i@f$.
     * @param[in] mu the RG scale to be used in the evaluation
     * @return @f$\delta y_\tau@f$
     */
    virtual const double deltaytau_HB(const double mu) const;

    /**
     * @brief The Higgs-basis coupling @f$\delta y_c@f$.
     * (See LHCHXSWG-INT-2015-001 document.)
     * Note that the Lagrangian definition of the Higgs-basis parameters coincides with the one of 
     * some of the @f$g_i, \delta g_i@f$ couplings defined above.
     * In the Higgs basis, however, one uses the freedom to perform certain field redefinitions and 
     * operations to demand that the mass eigenstate Lagrangian has specific features. (See pag. 5,6 in the reference.)
     * Therefore, the actual expression in terms of dim 6 coefficients may differ from the one for @f$g_i, \delta g_i@f$.
     * @param[in] mu the RG scale to be used in the evaluation
     * @return @f$\delta y_c@f$
     */
    virtual const double deltayc_HB(const double mu) const;
    

    /**
     * @brief The Higgs-basis coupling @f$\delta y_s@f$.
     * (See LHCHXSWG-INT-2015-001 document.)
     * Note that the Lagrangian definition of the Higgs-basis parameters coincides with the one of 
     * some of the @f$g_i, \delta g_i@f$ couplings defined above.
     * In the Higgs basis, however, one uses the freedom to perform certain field redefinitions and 
     * operations to demand that the mass eigenstate Lagrangian has specific features. (See pag. 5,6 in the reference.)
     * Therefore, the actual expression in terms of dim 6 coefficients may differ from the one for @f$g_i, \delta g_i@f$.
     * @param[in] mu the RG scale to be used in the evaluation
     * @return @f$\delta y_s@f$
     */
    virtual const double deltays_HB(const double mu) const;

    /**
     * @brief The Higgs-basis coupling @f$\delta y_\mu@f$.
     * (See LHCHXSWG-INT-2015-001 document.) 
     * Note that the Lagrangian definition of the Higgs-basis parameters coincides with the one of 
     * some of the @f$g_i, \delta g_i@f$ couplings defined above.
     * In the Higgs basis, however, one uses the freedom to perform certain field redefinitions and 
     * operations to demand that the mass eigenstate Lagrangian has specific features. (See pag. 5,6 in the reference.)
     * Therefore, the actual expression in terms of dim 6 coefficients may differ from the one for @f$g_i, \delta g_i@f$.
     * @param[in] mu the RG scale to be used in the evaluation
     * @return @f$\delta y_\mu@f$
     */
    virtual const double deltaymu_HB(const double mu) const;

    /**
     * @brief The Higgs-basis coupling @f$\delta c_z@f$.
     * (See LHCHXSWG-INT-2015-001 document.)
     * Note that the Lagrangian definition of the Higgs-basis parameters coincides with the one of 
     * some of the @f$g_i, \delta g_i@f$ couplings defined above.
     * In the Higgs basis, however, one uses the freedom to perform certain field redefinitions and 
     * operations to demand that the mass eigenstate Lagrangian has specific features. (See pag. 5,6 in the reference.)
     * Therefore, the actual expression in terms of dim 6 coefficients may differ from the one for @f$g_i, \delta g_i@f$.
     * @param[in] mu the RG scale to be used in the evaluation
     * @return @f$\delta c_z@f$
     */
    virtual const double deltacZ_HB(const double mu) const;

    /**
     * @brief The Higgs-basis coupling @f$c_{z\Box}@f$.
     * (See LHCHXSWG-INT-2015-001 document.)
     * Note that the Lagrangian definition of the Higgs-basis parameters coincides with the one of 
     * some of the @f$g_i, \delta g_i@f$ couplings defined above.
     * In the Higgs basis, however, one uses the freedom to perform certain field redefinitions and 
     * operations to demand that the mass eigenstate Lagrangian has specific features. (See pag. 5,6 in the reference.)
     * Therefore, the actual expression in terms of dim 6 coefficients may differ from the one for @f$g_i, \delta g_i@f$.
     * @param[in] mu the RG scale to be used in the evaluation
     * @return @f$c_{z\Box}@f$
     */
    virtual const double cZBox_HB(const double mu) const;

    /**
     * @brief The Higgs-basis coupling @f$c_{zz}@f$.
     * (See LHCHXSWG-INT-2015-001 document.)
     * Note that the Lagrangian definition of the Higgs-basis parameters coincides with the one of 
     * some of the @f$g_i, \delta g_i@f$ couplings defined above.
     * In the Higgs basis, however, one uses the freedom to perform certain field redefinitions and 
     * operations to demand that the mass eigenstate Lagrangian has specific features. (See pag. 5,6 in the reference.)
     * Therefore, the actual expression in terms of dim 6 coefficients may differ from the one for @f$g_i, \delta g_i@f$.
     * @param[in] mu the RG scale to be used in the evaluation
     * @return @f$c_{zz}@f$
     */
    virtual const double cZZ_HB(const double mu) const;

    /**
     * @brief The Higgs-basis coupling @f$c_{z\gamma}@f$.
     * (See LHCHXSWG-INT-2015-001 document.)
     * Note that the Lagrangian definition of the Higgs-basis parameters coincides with the one of 
     * some of the @f$g_i, \delta g_i@f$ couplings defined above.
     * In the Higgs basis, however, one uses the freedom to perform certain field redefinitions and 
     * operations to demand that the mass eigenstate Lagrangian has specific features. (See pag. 5,6 in the reference.)
     * Therefore, the actual expression in terms of dim 6 coefficients may differ from the one for @f$g_i, \delta g_i@f$.
     * @param[in] mu the RG scale to be used in the evaluation
     * @return @f$c_{z\gamma}@f$
     */
    virtual const double cZga_HB(const double mu) const;

    /**
     * @brief The Higgs-basis coupling @f$c_{\gamma\gamma}@f$.
     * (See LHCHXSWG-INT-2015-001 document.)
     * Note that the Lagrangian definition of the Higgs-basis parameters coincides with the one of 
     * some of the @f$g_i, \delta g_i@f$ couplings defined above.
     * In the Higgs basis, however, one uses the freedom to perform certain field redefinitions and 
     * operations to demand that the mass eigenstate Lagrangian has specific features. (See pag. 5,6 in the reference.)
     * Therefore, the actual expression in terms of dim 6 coefficients may differ from the one for @f$g_i, \delta g_i@f$.
     * @param[in] mu the RG scale to be used in the evaluation
     * @return @f$c_{\gamma\gamma}@f$
     */
    virtual const double cgaga_HB(const double mu) const;

    /**
     * @brief The Higgs-basis coupling @f$c_{gg}@f$.
     * (See LHCHXSWG-INT-2015-001 document.)
     * Note that the Lagrangian definition of the Higgs-basis parameters coincides with the one of 
     * some of the @f$g_i, \delta g_i@f$ couplings defined above.
     * In the Higgs basis, however, one uses the freedom to perform certain field redefinitions and 
     * operations to demand that the mass eigenstate Lagrangian has specific features. (See pag. 5,6 in the reference.)
     * Therefore, the actual expression in terms of dim 6 coefficients may differ from the one for @f$g_i, \delta g_i@f$.
     * @param[in] mu the RG scale to be used in the evaluation
     * @return @f$c_{gg}@f$
     */
    virtual const double cgg_HB(const double mu) const;

    /**
     * @brief The effective Higgs-basis coupling @f$c_{gg}^{Eff}@f$. (Similar to cgg_HB but including modifications of SM loops.)
     * (See arXiv: 1505.00046 [hep-ph] document.)
     * Note that the Lagrangian definition of the Higgs-basis parameters coincides with the one of 
     * some of the @f$g_i, \delta g_i@f$ couplings defined above.
     * In the Higgs basis, however, one uses the freedom to perform certain field redefinitions and 
     * operations to demand that the mass eigenstate Lagrangian has specific features. (See pag. 5,6 in the reference.)
     * Therefore, the actual expression in terms of dim 6 coefficients may differ from the one for @f$g_i, \delta g_i@f$.
     * @param[in] mu the RG scale to be used in the evaluation
     * @return @f$c_{gg}^{Eff}@f$
     */
    virtual const double cggEff_HB(const double mu) const;

    /**
     * @brief The Higgs-basis coupling @f$\lambda_{z}@f$.
     * (See LHCHXSWG-INT-2015-001 document.)
     * Note that the Lagrangian definition of the Higgs-basis parameters coincides with the one of 
     * some of the @f$g_i, \delta g_i@f$ couplings defined above.
     * In the Higgs basis, however, one uses the freedom to perform certain field redefinitions and 
     * operations to demand that the mass eigenstate Lagrangian has specific features. (See pag. 5,6 in the reference.)
     * Therefore, the actual expression in terms of dim 6 coefficients may differ from the one for @f$g_i, \delta g_i@f$.
     * @param[in] mu the RG scale to be used in the evaluation
     * @return @f$\lambda_{z}@f$
     */
    virtual const double lambz_HB(const double mu) const;

    /////////////Combinations of Warsaw basis coefficients constrained by EWPO////////////////

    /**
     * @brief Combination of coefficients of the Warsaw basis constrained by EWPO
     * @f$(\hat{C}_{HL}^{(1)})_{11}@f$.
     * @return @f$(\hat{C}_{HL}^{(1)})_{11}@f$
     */
    virtual const double CEWHL111() const;

    /**
     * @brief Combination of coefficients of the Warsaw basis constrained by EWPO
     * @f$(\hat{C}_{HL}^{(1)})_{22}@f$.
     * @return @f$(\hat{C}_{HL}^{(1)})_{22}@f$
     */
    virtual const double CEWHL122() const;

    /**
     * @brief Combination of coefficients of the Warsaw basis constrained by EWPO
     * @f$(\hat{C}_{HL}^{(1)})_{33}@f$.
     * @return @f$(\hat{C}_{HL}^{(1)})_{33}@f$
     */
    virtual const double CEWHL133() const;

    /**
     * @brief Combination of coefficients of the Warsaw basis constrained by EWPO
     * @f$(\hat{C}_{HL}^{(3)})_{11}@f$.
     * @return @f$(\hat{C}_{HL}^{(3)})_{11}@f$
     */
    virtual const double CEWHL311() const;

    /**
     * @brief Combination of coefficients of the Warsaw basis constrained by EWPO
     * @f$(\hat{C}_{HL}^{(3)})_{22}@f$.
     * @return @f$(\hat{C}_{HL}^{(3)})_{22}@f$
     */
    virtual const double CEWHL322() const;

    /**
     * @brief Combination of coefficients of the Warsaw basis constrained by EWPO
     * @f$(\hat{C}_{HL}^{(3)})_{33}@f$.
     * @return @f$(\hat{C}_{HL}^{(3)})_{33}@f$
     */
    virtual const double CEWHL333() const;

    /**
     * @brief Combination of coefficients of the Warsaw basis constrained by EWPO
     * @f$(\hat{C}_{HQ}^{(1)})_{11}@f$.
     * @return @f$(\hat{C}_{HQ}^{(1)})_{11}@f$
     */
    virtual const double CEWHQ111() const;

    /**
     * @brief Combination of coefficients of the Warsaw basis constrained by EWPO
     * @f$(\hat{C}_{HQ}^{(1)})_{22}@f$.
     * @return @f$(\hat{C}_{HQ}^{(1)})_{22}@f$
     */
    virtual const double CEWHQ122() const;

    /**
     * @brief Combination of coefficients of the Warsaw basis constrained by EWPO
     * @f$(\hat{C}_{HQ}^{(1)})_{33}@f$.
     * @return @f$(\hat{C}_{HQ}^{(1)})_{33}@f$
     */
    virtual const double CEWHQ133() const;

    /**
     * @brief Combination of coefficients of the Warsaw basis constrained by EWPO
     * @f$(\hat{C}_{HQ}^{(3)})_{11}@f$.
     * @return @f$(\hat{C}_{HQ}^{(3)})_{11}@f$
     */
    virtual const double CEWHQ311() const;

    /**
     * @brief Combination of coefficients of the Warsaw basis constrained by EWPO
     * @f$(\hat{C}_{HQ}^{(3)})_{22}@f$.
     * @return @f$(\hat{C}_{HQ}^{(3)})_{22}@f$
     */
    virtual const double CEWHQ322() const;

    /**
     * @brief Combination of coefficients of the Warsaw basis constrained by EWPO
     * @f$(\hat{C}_{HQ}^{(3)})_{33}@f$.
     * @return @f$(\hat{C}_{HQ}^{(3)})_{33}@f$
     */
    virtual const double CEWHQ333() const;


    /**
     * @brief Combination of coefficients of the Warsaw basis constrained by EWPO
     * @f$(\hat{C}_{HQ}^{(d)})_{33}@f$.
     * @return @f$(\hat{C}_{HQ}^{(d)})_{33}@f$
     */
    virtual const double CEWHQd33() const;
    
    /**
     * @brief Combination of coefficients of the Warsaw basis not constrained by EWPO (at LO)
     * @f$(\hat{C}_{HQ}^{(u)})_{33}@f$.
     * @return @f$(\hat{C}_{HQ}^{(u)})_{33}@f$
     */
    virtual const double CEWHQu33() const;


    /**
     * @brief Combination of coefficients of the Warsaw basis constrained by EWPO
     * @f$(\hat{C}_{He})_{11}@f$.
     * @return @f$(\hat{C}_{He})_{11}@f$
     */
    virtual const double CEWHe11() const;

    /**
     * @brief Combination of coefficients of the Warsaw basis constrained by EWPO
     * @f$(\hat{C}_{He})_{22}@f$.
     * @return @f$(\hat{C}_{He})_{22}@f$
     */
    virtual const double CEWHe22() const;

    /**
     * @brief Combination of coefficients of the Warsaw basis constrained by EWPO
     * @f$(\hat{C}_{He})_{33}@f$.
     * @return @f$(\hat{C}_{He})_{33}@f$
     */
    virtual const double CEWHe33() const;


    /**
     * @brief Combination of coefficients of the Warsaw basis constrained by EWPO
     * @f$(\hat{C}_{Hu})_{11}@f$.
     * @return @f$(\hat{C}_{Hu})_{11}@f$
     */
    virtual const double CEWHu11() const;

    /**
     * @brief Combination of coefficients of the Warsaw basis constrained by EWPO
     * @f$(\hat{C}_{Hu})_{22}@f$.
     * @return @f$(\hat{C}_{Hu})_{22}@f$
     */
    virtual const double CEWHu22() const;

    /**
     * @brief Combination of coefficients of the Warsaw basis constrained by EWPO
     * @f$(\hat{C}_{Hu})_{33}@f$.
     * @return @f$(\hat{C}_{Hu})_{33}@f$
     */
    virtual const double CEWHu33() const;


    /**
     * @brief Combination of coefficients of the Warsaw basis constrained by EWPO
     * @f$(\hat{C}_{Hd})_{11}@f$.
     * @return @f$(\hat{C}_{Hd})_{11}@f$
     */
    virtual const double CEWHd11() const;

    /**
     * @brief Combination of coefficients of the Warsaw basis constrained by EWPO
     * @f$(\hat{C}_{Hd})_{22}@f$.
     * @return @f$(\hat{C}_{Hd})_{22}@f$
     */
    virtual const double CEWHd22() const;

    /**
     * @brief Combination of coefficients of the Warsaw basis constrained by EWPO
     * @f$(\hat{C}_{Hd})_{33}@f$.
     * @return @f$(\hat{C}_{Hd})_{33}@f$
     */
    virtual const double CEWHd33() const;
    
    ///////////Collider observables: LHC dilepton events////////////////////////
    
    /**
     * @brief Number of di-electron events at the LHC at 13 TeV
     * @return NevLHCppee13
     */
    virtual const double NevLHCppee13(const int i_bin) const;
    
    /**
     * @brief Number of di-muon events at the LHC at 13 TeV
     * @return NevLHCppmumu13
     */
    virtual const double NevLHCppmumu13(const int i_bin) const;
    
    /**
     * @brief Number of di-tau events at the LHC at 13 TeV
     * @return NevLHCpptautau13
     */
    virtual const double NevLHCpptautau13(const int i_bin) const;
    
    ///////////Collider observables: LHC mono-lepton events////////////////////////
    
    /**
     * @brief Number of mono-electron events at the LHC at 13 TeV
     * @return NevLHCppenu13
     */
    virtual const double NevLHCppenu13(const int i_bin) const;
    
    /**
     * @brief Number of mono-muon events at the LHC at 13 TeV
     * @return NevLHCppmunu13
     */
    virtual const double NevLHCppmunu13(const int i_bin) const;
    
    /**
     * @brief Number of mono-tau events at the LHC at 13 TeV
     * @return NevLHCpptaunu13
     */
    virtual const double NevLHCpptaunu13(const int i_bin) const;

    /////////////Auxiliary observables////////////////

    /**
     * @brief Auxiliary observable AuxObs_NP1
     * (See code for details.)
     * @return AuxObs_NP1
     */
    virtual const double AuxObs_NP1() const;

    /**
     * @brief Auxiliary observable AuxObs_NP2
     * (See code for details.)
     * @return AuxObs_NP2
     */
    virtual const double AuxObs_NP2() const;

    /**
     * @brief Auxiliary observable AuxObs_NP3
     * (See code for details.)
     * @return AuxObs_NP3
     */
    virtual const double AuxObs_NP3() const;

    /**
     * @brief Auxiliary observable AuxObs_NP4
     * (See code for details.)
     * @return AuxObs_NP4
     */
    virtual const double AuxObs_NP4() const;

    /**
     * @brief Auxiliary observable AuxObs_NP5
     * (See code for details.)
     * @return AuxObs_NP5
     */
    virtual const double AuxObs_NP5() const;

    /**
     * @brief Auxiliary observable AuxObs_NP6
     * (See code for details.)
     * @return AuxObs_NP6
     */
    virtual const double AuxObs_NP6() const;

    /**
     * @brief Auxiliary observable AuxObs_NP7
     * (See code for details.)
     * @return AuxObs_NP7
     */
    virtual const double AuxObs_NP7() const;

    /**
     * @brief Auxiliary observable AuxObs_NP8
     * (See code for details.)
     * @return AuxObs_NP8
     */
    virtual const double AuxObs_NP8() const;

    /**
     * @brief Auxiliary observable AuxObs_NP9
     * (See code for details.)
     * @return AuxObs_NP9
     */
    virtual const double AuxObs_NP9() const;

    /**
     * @brief Auxiliary observable AuxObs_NP10
     * (See code for details.)
     * @return AuxObs_NP10
     */
    virtual const double AuxObs_NP10() const;

    /**
     * @brief Auxiliary observable AuxObs_NP11
     * (See code for details.)
     * @return AuxObs_NP11
     */
    virtual const double AuxObs_NP11() const;

    /**
     * @brief Auxiliary observable AuxObs_NP12
     * (See code for details.)
     * @return AuxObs_NP12
     */
    virtual const double AuxObs_NP12() const;

    /**
     * @brief Auxiliary observable AuxObs_NP13
     * @return AuxObs_NP13
     */
    virtual const double AuxObs_NP13() const;

    /**
     * @brief Auxiliary observable AuxObs_NP14
     * @return AuxObs_NP14
     */
    virtual const double AuxObs_NP14() const;

    /**
     * @brief Auxiliary observable AuxObs_NP15
     * @return AuxObs_NP15
     */
    virtual const double AuxObs_NP15() const;

    /**
     * @brief Auxiliary observable AuxObs_NP16
     * @return AuxObs_NP16
     */
    virtual const double AuxObs_NP16() const;

    /**
     * @brief Auxiliary observable AuxObs_NP17
     * @return AuxObs_NP17
     */
    virtual const double AuxObs_NP17() const;

    /**
     * @brief Auxiliary observable AuxObs_NP18
     * @return AuxObs_NP18
     */
    virtual const double AuxObs_NP18() const;

    /**
     * @brief Auxiliary observable AuxObs_NP19
     * @return AuxObs_NP19
     */
    virtual const double AuxObs_NP19() const;

    /**
     * @brief Auxiliary observable AuxObs_NP20
     * @return AuxObs_NP20
     */
    virtual const double AuxObs_NP20() const;

    /**
     * @brief Auxiliary observable AuxObs_NP21
     * (See code for details.)
     * @return AuxObs_NP21
     */
    virtual const double AuxObs_NP21() const;

    /**
     * @brief Auxiliary observable AuxObs_NP22
     * (See code for details.)
     * @return AuxObs_NP22
     */
    virtual const double AuxObs_NP22() const;

    /**
     * @brief Auxiliary observable AuxObs_NP23
     * @return AuxObs_NP23
     */
    virtual const double AuxObs_NP23() const;

    /**
     * @brief Auxiliary observable AuxObs_NP24
     * @return AuxObs_NP24
     */
    virtual const double AuxObs_NP24() const;

    /**
     * @brief Auxiliary observable AuxObs_NP25
     * @return AuxObs_NP25
     */
    virtual const double AuxObs_NP25() const;

    /**
     * @brief Auxiliary observable AuxObs_NP26
     * @return AuxObs_NP26
     */
    virtual const double AuxObs_NP26() const;

    /**
     * @brief Auxiliary observable AuxObs_NP27
     * @return AuxObs_NP27
     */
    virtual const double AuxObs_NP27() const;

    /**
     * @brief Auxiliary observable AuxObs_NP28
     * @return AuxObs_NP28
     */
    virtual const double AuxObs_NP28() const;

    /**
     * @brief Auxiliary observable AuxObs_NP29
     * @return AuxObs_NP29
     */
    virtual const double AuxObs_NP29() const;

    /**
     * @brief Auxiliary observable AuxObs_NP30
     * @return AuxObs_NP20
     */
    virtual const double AuxObs_NP30() const;

    //AG:begin
    virtual int OutputOrder() const; ///<Type of contributions to be included in the EWPOs. Takes a numerica values depending on the choice.
    bool hatCis() const; ///<If True, explicitly defines the 8 'hat' coefficients in the EWPOs (Z-couplings, dGf, W-width)              
    bool flagCHWpCHB() const; ///< If True, uses the coefficient CHWpCHW instead of the sum CiHW+CiHB.
    //AG:end

    ////////////////////////////////////////////////////////////////////////
protected:

    /**
     * @brief @copybrief Model::setParameter()
     * @copydetails Model::setParameter()
     */
    virtual void setParameter(const std::string name, const double & value);
    mutable Matching<NPSMEFTd6Matching, NPSMEFTd6> NPSMEFTd6M;

    //AGhat:start
    double CHL1hat;
    double CHL3hat;
    double CHQ1hat;
    double CHQ3hat;
    double CHdhat;
    double CHuhat;
    double CHehat;
    double CLLhat;
    double CHWpCHB;
    //AGhat:end
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
    double CQQ1_1133, CQQ1_1331, CQQ1_2233, CQQ1_2332, CQQ1_3333, CQQ3_1133, CQQ3_1331, CQQ3_2233, CQQ3_2332, CQQ3_3333;
    double Cuu_1133, Cuu_1331, Cuu_2233, Cuu_2332, Cuu_3333, Cud1_3311, Cud1_3322, Cud1_3333, Cud8_3311, Cud8_3322, Cud8_3333;
    double CQu1_1133, CQu1_3311, CQu1_2233, CQu1_3322, CQu1_3333, CQu8_1133, CQu8_3311, CQu8_2233, CQu8_3322, CQu8_3333;
    double CQd1_3311, CQd1_3322, CQd1_3333, CQd8_3311, CQd8_3322, CQd8_3333;
    double CQuQd1_3333, CQuQd8_3333;
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
    double eeeWWint, edeeWWdcint; ///< Intrinsic relative theoretical error in \f$e^+ e^- \to W^+ W^-\f$: total cross section and \f$d/d\cos{\theta}\f$ distribution.

    double eggFHgaga, eggFHZga, eggFHZZ, eggFHWW, eggFHtautau, eggFHbb, eggFHmumu; ///< Total relative theoretical error in \f$gg \to H \to X\f$.   
    double eVBFHgaga, eVBFHZga, eVBFHZZ, eVBFHWW, eVBFHtautau, eVBFHbb, eVBFHmumu; ///< Total relative theoretical error in \f$pp \to Hjj (VBF) \to X jj\f$.   
    double eWHgaga, eWHZga, eWHZZ, eWHWW, eWHtautau, eWHbb, eWHmumu; ///< Total relative theoretical error in \f$pp \to WH \to W X\f$.    
    double eZHgaga, eZHZga, eZHZZ, eZHWW, eZHtautau, eZHbb, eZHmumu; ///< Total relative theoretical error in \f$pp \to ZH \to Z X\f$.
    double ettHgaga, ettHZga, ettHZZ, ettHWW, ettHtautau, ettHbb, ettHmumu; ///< Total relative theoretical error in \f$pp \to ttH \to tt X\f$.
    double eVBFHinv, eVHinv; ///< Total relative theoretical error in \f$pp \to X H \to X + invisible\f$.

    double eHwidth; ///< Total relative theoretical error in the Higgs width

    double nuisP1, nuisP2, nuisP3, nuisP4, nuisP5, nuisP6, nuisP7, nuisP8, nuisP9, nuisP10; ///< Nuisance parameters to be used in observables

    // The error in the parameters multiplying the dimension-6 operator coefficients in the production cross sections.
    double eVBF_2_Hbox; ///< Theoretical uncertainty in the (linear) new physics contribution from \f$C_{H\Box}\f$ to VBF production at Tevatron (1.96 TeV).
    double eVBF_2_HQ1_11; ///< Theoretical uncertainty in the (linear) new physics contribution from \f$(C_{HQ}^{(1)})_{11}\f$ to VBF production at Tevatron (1.96 TeV).
    double eVBF_2_Hu_11; ///< Theoretical uncertainty in the (linear) new physics contribution from \f$(C_{Hu})_{11}\f$ to VBF production at Tevatron (1.96 TeV).
    double eVBF_2_Hd_11; ///< Theoretical uncertainty in the (linear) new physics contribution from \f$(C_{Hd})_{11}\f$ to VBF production at Tevatron (1.96 TeV).
    double eVBF_2_HQ3_11; ///< Theoretical uncertainty in the (linear) new physics contribution from \f$(C_{HQ}^{(3)})_{11}\f$ to VBF production at Tevatron (1.96 TeV).
    double eVBF_2_HD; ///< Theoretical uncertainty in the (linear) new physics contribution from \f$C_{HD}\f$ to VBF production at Tevatron (1.96 TeV).
    double eVBF_2_HB; ///< Theoretical uncertainty in the (linear) new physics contribution from \f$C_{HB}\f$ to VBF production at Tevatron (1.96 TeV).
    double eVBF_2_HW; ///< Theoretical uncertainty in the (linear) new physics contribution from \f$C_{HW}\f$ to VBF production at Tevatron (1.96 TeV).
    double eVBF_2_HWB; ///< Theoretical uncertainty in the (linear) new physics contribution from \f$C_{HWB}\f$ to VBF production at Tevatron (1.96 TeV).
    double eVBF_2_HG; ///< Theoretical uncertainty in the (linear) new physics contribution from \f$C_{HG}\f$ to VBF production at Tevatron (1.96 TeV).
    double eVBF_2_DHB; ///< Theoretical uncertainty in the (linear) new physics contribution from \f$C_{DHB}\f$ to VBF production at Tevatron (1.96 TeV).
    double eVBF_2_DHW; ///< Theoretical uncertainty in the (linear) new physics contribution from \f$C_{DHW}\f$ to VBF production at Tevatron (1.96 TeV).
    double eVBF_2_DeltaGF; ///< Theoretical uncertainty in the (linear) new physics contribution from \f$\delta_{G_F}\f$ to VBF production at Tevatron (1.96 TeV).

    double eVBF_78_Hbox; ///< Theoretical uncertainty in the (linear) new physics contribution from \f$C_{H\Box}\f$ to VBF production at Tevatron (7 & 8 TeV).
    double eVBF_78_HQ1_11; ///< Theoretical uncertainty in the (linear) new physics contribution from \f$(C_{HQ}^{(1)})_{11}\f$ to VBF production at Tevatron (7 & 8 TeV).
    double eVBF_78_Hu_11; ///< Theoretical uncertainty in the (linear) new physics contribution from \f$(C_{Hu})_{11}\f$ to VBF production at Tevatron (7 & 8 TeV).
    double eVBF_78_Hd_11; ///< Theoretical uncertainty in the (linear) new physics contribution from \f$(C_{Hd})_{11}\f$ to VBF production at Tevatron (7 & 8 TeV).
    double eVBF_78_HQ3_11; ///< Theoretical uncertainty in the (linear) new physics contribution from \f$(C_{HQ}^{(3)})_{11}\f$ to VBF production at Tevatron (7 & 8 TeV).
    double eVBF_78_HD; ///< Theoretical uncertainty in the (linear) new physics contribution from \f$C_{HD}\f$ to VBF production at Tevatron (7 & 8 TeV).
    double eVBF_78_HB; ///< Theoretical uncertainty in the (linear) new physics contribution from \f$C_{HB}\f$ to VBF production at Tevatron (7 & 8 TeV).
    double eVBF_78_HW; ///< Theoretical uncertainty in the (linear) new physics contribution from \f$C_{HW}\f$ to VBF production at Tevatron (7 & 8 TeV).
    double eVBF_78_HWB; ///< Theoretical uncertainty in the (linear) new physics contribution from \f$C_{HWB}\f$ to VBF production at Tevatron (7 & 8 TeV).
    double eVBF_78_HG; ///< Theoretical uncertainty in the (linear) new physics contribution from \f$C_{HG}\f$ to VBF production at Tevatron (7 & 8 TeV).
    double eVBF_78_DHB; ///< Theoretical uncertainty in the (linear) new physics contribution from \f$C_{DHB}\f$ to VBF production at Tevatron (7 & 8 TeV).
    double eVBF_78_DHW; ///< Theoretical uncertainty in the (linear) new physics contribution from \f$C_{DHW}\f$ to VBF production at Tevatron (7 & 8 TeV).
    double eVBF_78_DeltaGF; ///< Theoretical uncertainty in the (linear) new physics contribution from \f$\delta_{G_F}\f$ to VBF production at Tevatron (7 & 8 TeV).

    double eVBF_1314_Hbox; ///< Theoretical uncertainty in the (linear) new physics contribution from \f$C_{H\Box}\f$ to VBF production at Tevatron (13 & 14 TeV).
    double eVBF_1314_HQ1_11; ///< Theoretical uncertainty in the (linear) new physics contribution from \f$(C_{HQ}^{(1)})_{11}\f$ to VBF production at Tevatron (13 & 14 TeV).
    double eVBF_1314_Hu_11; ///< Theoretical uncertainty in the (linear) new physics contribution from \f$(C_{Hu})_{11}\f$ to VBF production at Tevatron (13 & 14 TeV).
    double eVBF_1314_Hd_11; ///< Theoretical uncertainty in the (linear) new physics contribution from \f$(C_{Hd})_{11}\f$ to VBF production at Tevatron (13 & 14 TeV).
    double eVBF_1314_HQ3_11; ///< Theoretical uncertainty in the (linear) new physics contribution from \f$(C_{HQ}^{(3)})_{11}\f$ to VBF production at Tevatron (13 & 14 TeV).
    double eVBF_1314_HD; ///< Theoretical uncertainty in the (linear) new physics contribution from \f$C_{HD}\f$ to VBF production at Tevatron (13 & 14 TeV).
    double eVBF_1314_HB; ///< Theoretical uncertainty in the (linear) new physics contribution from \f$C_{HB}\f$ to VBF production at Tevatron (13 & 14 TeV).
    double eVBF_1314_HW; ///< Theoretical uncertainty in the (linear) new physics contribution from \f$C_{HW}\f$ to VBF production at Tevatron (13 & 14 TeV).
    double eVBF_1314_HWB; ///< Theoretical uncertainty in the (linear) new physics contribution from \f$C_{HWB}\f$ to VBF production at Tevatron (13 & 14 TeV).
    double eVBF_1314_HG; ///< Theoretical uncertainty in the (linear) new physics contribution from \f$C_{HG}\f$ to VBF production at Tevatron (13 & 14 TeV).
    double eVBF_1314_DHB; ///< Theoretical uncertainty in the (linear) new physics contribution from \f$C_{DHB}\f$ to VBF production at Tevatron (13 & 14 TeV).
    double eVBF_1314_DHW; ///< Theoretical uncertainty in the (linear) new physics contribution from \f$C_{DHW}\f$ to VBF production at Tevatron (13 & 14 TeV).
    double eVBF_1314_DeltaGF; ///< Theoretical uncertainty in the (linear) new physics contribution from \f$\delta_{G_F}\f$ to VBF production at Tevatron (13 & 14 TeV).

    double eWH_2_Hbox; ///< Theoretical uncertainty in the (linear) new physics contribution from \f$C_{H\Box}\f$ to WH production at Tevatron (1.96 TeV).
    double eWH_2_HQ3_11; ///< Theoretical uncertainty in the (linear) new physics contribution from \f$(C_{HQ}^{(3)})_{11}\f$ to WH production at Tevatron (1.96 TeV).
    double eWH_2_HD; ///< Theoretical uncertainty in the (linear) new physics contribution from \f$C_{HD}\f$ to WH production at Tevatron (1.96 TeV).
    double eWH_2_HW; ///< Theoretical uncertainty in the (linear) new physics contribution from \f$C_{HW}\f$ to WH production at Tevatron (1.96 TeV).
    double eWH_2_HWB; ///< Theoretical uncertainty in the (linear) new physics contribution from \f$C_{HWB}\f$ to WH production at Tevatron (1.96 TeV).
    double eWH_2_DHW; ///< Theoretical uncertainty in the (linear) new physics contribution from \f$C_{DHW}\f$ to WH production at the LHC (1.96 TeV).
    double eWH_2_DeltaGF; ///< Theoretical uncertainty in the (linear) new physics contribution from \f$\delta_{G_F}\f$ to WH production at the LHC (1.96 TeV).

    double eWH_78_Hbox; ///< Theoretical uncertainty in the (linear) new physics contribution from \f$C_{H\Box}\f$ to WH production at Tevatron (7 & 8 TeV).
    double eWH_78_HQ3_11; ///< Theoretical uncertainty in the (linear) new physics contribution from \f$(C_{HQ}^{(3)})_{11}\f$ to WH production at Tevatron (7 & 8 TeV).
    double eWH_78_HD; ///< Theoretical uncertainty in the (linear) new physics contribution from \f$C_{HD}\f$ to WH production at Tevatron (7 & 8 TeV).
    double eWH_78_HW; ///< Theoretical uncertainty in the (linear) new physics contribution from \f$C_{HW}\f$ to WH production at Tevatron (7 & 8 TeV).
    double eWH_78_HWB; ///< Theoretical uncertainty in the (linear) new physics contribution from \f$C_{HWB}\f$ to WH production at Tevatron (7 & 8 TeV).
    double eWH_78_DHW; ///< Theoretical uncertainty in the (linear) new physics contribution from \f$C_{DHW}\f$ to WH production at the LHC (7 & 8 TeV).
    double eWH_78_DeltaGF; ///< Theoretical uncertainty in the (linear) new physics contribution from \f$\delta_{G_F}\f$ to WH production at the LHC (7 & 8 TeV).

    double eWH_1314_Hbox; ///< Theoretical uncertainty in the (linear) new physics contribution from \f$C_{H\Box}\f$ to WH production at Tevatron (13 & 14 TeV).
    double eWH_1314_HQ3_11; ///< Theoretical uncertainty in the (linear) new physics contribution from \f$(C_{HQ}^{(3)})_{11}\f$ to WH production at Tevatron (13 & 14 TeV).
    double eWH_1314_HD; ///< Theoretical uncertainty in the (linear) new physics contribution from \f$C_{HD}\f$ to WH production at Tevatron (13 & 14 TeV).
    double eWH_1314_HW; ///< Theoretical uncertainty in the (linear) new physics contribution from \f$C_{HW}\f$ to WH production at Tevatron (13 & 14 TeV).
    double eWH_1314_HWB; ///< Theoretical uncertainty in the (linear) new physics contribution from \f$C_{HWB}\f$ to WH production at Tevatron (13 & 14 TeV).
    double eWH_1314_DHW; ///< Theoretical uncertainty in the (linear) new physics contribution from \f$C_{DHW}\f$ to WH production at the LHC (13 & 14 TeV).
    double eWH_1314_DeltaGF; ///< Theoretical uncertainty in the (linear) new physics contribution from \f$\delta_{G_F}\f$ to WH production at the LHC (13 & 14 TeV).

    double eZH_2_Hbox; ///< Theoretical uncertainty in the (linear) new physics contribution from \f$C_{H\Box}\f$ to ZH production at Tevatron (1.96 TeV).
    double eZH_2_HQ1_11; ///< Theoretical uncertainty in the (linear) new physics contribution from \f$(C_{HQ}^{(1)})_{11}\f$ to ZH production at Tevatron (1.96 TeV).
    double eZH_2_Hu_11; ///< Theoretical uncertainty in the (linear) new physics contribution from \f$(C_{Hu})_{11}\f$ to ZH production at Tevatron (1.96 TeV).
    double eZH_2_Hd_11; ///< Theoretical uncertainty in the (linear) new physics contribution from \f$(C_{Hd})_{11}\f$ to ZH production at Tevatron (1.96 TeV).
    double eZH_2_HQ3_11; ///< Theoretical uncertainty in the (linear) new physics contribution from \f$(C_{HQ}^{(3)})_{11}\f$ to ZH production at Tevatron (1.96 TeV).
    double eZH_2_HD; ///< Theoretical uncertainty in the (linear) new physics contribution from \f$C_{HD}\f$ to ZH production at Tevatron (1.96 TeV).
    double eZH_2_HB; ///< Theoretical uncertainty in the (linear) new physics contribution from \f$C_{HB}\f$ to ZH production at Tevatron (1.96 TeV).
    double eZH_2_HW; ///< Theoretical uncertainty in the (linear) new physics contribution from \f$C_{HW}\f$ to ZH production at Tevatron (1.96 TeV).
    double eZH_2_HWB; ///< Theoretical uncertainty in the (linear) new physics contribution from \f$C_{HWB}\f$ to ZH production at Tevatron (1.96 TeV).
    double eZH_2_DHB; ///< Theoretical uncertainty in the (linear) new physics contribution from \f$C_{DHB}\f$ to ZH production at Tevatron (1.96 TeV).
    double eZH_2_DHW; ///< Theoretical uncertainty in the (linear) new physics contribution from \f$C_{DHW}\f$ to ZH production at Tevatron (1.96 TeV).
    double eZH_2_DeltaGF; ///< Theoretical uncertainty in the (linear) new physics contribution from \f$\delta_{G_F}\f$ to ZH production at Tevatron (1.96 TeV).

    double eZH_78_Hbox; ///< Theoretical uncertainty in the (linear) new physics contribution from \f$C_{H\Box}\f$ to ZH production at Tevatron (7 & 8 TeV).
    double eZH_78_HQ1_11; ///< Theoretical uncertainty in the (linear) new physics contribution from \f$(C_{HQ}^{(1)})_{11}\f$ to ZH production at Tevatron (7 & 8 TeV).
    double eZH_78_Hu_11; ///< Theoretical uncertainty in the (linear) new physics contribution from \f$(C_{Hu})_{11}\f$ to ZH production at Tevatron (7 & 8 TeV).
    double eZH_78_Hd_11; ///< Theoretical uncertainty in the (linear) new physics contribution from \f$(C_{Hd})_{11}\f$ to ZH production at Tevatron (7 & 8 TeV).
    double eZH_78_HQ3_11; ///< Theoretical uncertainty in the (linear) new physics contribution from \f$(C_{HQ}^{(3)})_{11}\f$ to ZH production at Tevatron (7 & 8 TeV).
    double eZH_78_HD; ///< Theoretical uncertainty in the (linear) new physics contribution from \f$C_{HD}\f$ to ZH production at Tevatron (7 & 8 TeV).
    double eZH_78_HB; ///< Theoretical uncertainty in the (linear) new physics contribution from \f$C_{HB}\f$ to ZH production at Tevatron (7 & 8 TeV).
    double eZH_78_HW; ///< Theoretical uncertainty in the (linear) new physics contribution from \f$C_{HW}\f$ to ZH production at Tevatron (7 & 8 TeV).
    double eZH_78_HWB; ///< Theoretical uncertainty in the (linear) new physics contribution from \f$C_{HWB}\f$ to ZH production at Tevatron (7 & 8 TeV).
    double eZH_78_DHB; ///< Theoretical uncertainty in the (linear) new physics contribution from \f$C_{DHB}\f$ to ZH production at Tevatron (7 & 8 TeV).
    double eZH_78_DHW; ///< Theoretical uncertainty in the (linear) new physics contribution from \f$C_{DHW}\f$ to ZH production at Tevatron (7 & 8 TeV).
    double eZH_78_DeltaGF; ///< Theoretical uncertainty in the (linear) new physics contribution from \f$\delta_{G_F}\f$ to ZH production at Tevatron (7 & 8 TeV).

    double eZH_1314_Hbox; ///< Theoretical uncertainty in the (linear) new physics contribution from \f$C_{H\Box}\f$ to ZH production at Tevatron (13 & 14 TeV).
    double eZH_1314_HQ1_11; ///< Theoretical uncertainty in the (linear) new physics contribution from \f$(C_{HQ}^{(1)})_{11}\f$ to ZH production at Tevatron (13 & 14 TeV).
    double eZH_1314_Hu_11; ///< Theoretical uncertainty in the (linear) new physics contribution from \f$(C_{Hu})_{11}\f$ to ZH production at Tevatron (13 & 14 TeV).
    double eZH_1314_Hd_11; ///< Theoretical uncertainty in the (linear) new physics contribution from \f$(C_{Hd})_{11}\f$ to ZH production at Tevatron (13 & 14 TeV).
    double eZH_1314_HQ3_11; ///< Theoretical uncertainty in the (linear) new physics contribution from \f$(C_{HQ}^{(3)})_{11}\f$ to ZH production at Tevatron (13 & 14 TeV).
    double eZH_1314_HD; ///< Theoretical uncertainty in the (linear) new physics contribution from \f$C_{HD}\f$ to ZH production at Tevatron (13 & 14 TeV).
    double eZH_1314_HB; ///< Theoretical uncertainty in the (linear) new physics contribution from \f$C_{HB}\f$ to ZH production at Tevatron (13 & 14 TeV).
    double eZH_1314_HW; ///< Theoretical uncertainty in the (linear) new physics contribution from \f$C_{HW}\f$ to ZH production at Tevatron (13 & 14 TeV).
    double eZH_1314_HWB; ///< Theoretical uncertainty in the (linear) new physics contribution from \f$C_{HWB}\f$ to ZH production at Tevatron (13 & 14 TeV).
    double eZH_1314_DHB; ///< Theoretical uncertainty in the (linear) new physics contribution from \f$C_{DHB}\f$ to ZH production at Tevatron (13 & 14 TeV).
    double eZH_1314_DHW; ///< Theoretical uncertainty in the (linear) new physics contribution from \f$C_{DHW}\f$ to ZH production at Tevatron (13 & 14 TeV).
    double eZH_1314_DeltaGF; ///< Theoretical uncertainty in the (linear) new physics contribution from \f$\delta_{G_F}\f$ to ZH production at Tevatron (13 & 14 TeV).

    double ettH_2_HG; ///< Theoretical uncertainty in the (linear) new physics contribution from \f$C_{HG}\f$ to ttH production at Tevatron (1.96 TeV).
    double ettH_2_G; ///< Theoretical uncertainty in the (linear) new physics contribution from \f$C_{G}\f$ to ttH production at Tevatron (1.96 TeV).
    double ettH_2_uG_33r; ///< Theoretical uncertainty in the (linear) new physics contribution from \f$(C_{uG})_{33}\f$ to ttH production at the LHC (1.96 TeV).
    double ettH_2_DeltagHt; ///< Theoretical uncertainty in the (linear) new physics contribution from \f$\delta g_{Htt}\f$ to ttH production at the LHC (1.96 TeV).

    double ettH_78_HG; ///< Theoretical uncertainty in the (linear) new physics contribution from \f$C_{HG}\f$ to ttH production at Tevatron (7 & 8 TeV).
    double ettH_78_G; ///< Theoretical uncertainty in the (linear) new physics contribution from \f$C_{G}\f$ to ttH production at Tevatron (7 & 8 TeV).
    double ettH_78_uG_33r; ///< Theoretical uncertainty in the (linear) new physics contribution from \f$(C_{uG})_{33}\f$ to ttH production at the LHC (7 & 8 TeV).
    double ettH_78_DeltagHt; ///< Theoretical uncertainty in the (linear) new physics contribution from \f$\delta g_{Htt}\f$ to ttH production at the LHC (7 & 8 TeV).

    double ettH_1314_HG; ///< Theoretical uncertainty in the (linear) new physics contribution from \f$C_{HG}\f$ to ttH production at Tevatron (13 & 14 TeV).
    double ettH_1314_G; ///< Theoretical uncertainty in the (linear) new physics contribution from \f$C_{G}\f$ to ttH production at Tevatron (13 & 14 TeV).
    double ettH_1314_uG_33r; ///< Theoretical uncertainty in the (linear) new physics contribution from \f$(C_{uG})_{33}\f$ to ttH production at the LHC (13 & 14 TeV).
    double ettH_1314_DeltagHt; ///< Theoretical uncertainty in the (linear) new physics contribution from \f$\delta g_{Htt}\f$ to ttH production at the LHC (13 & 14 TeV).

    double BrHinv; ///< The branching ratio of invisible Higgs decays.
    double BrHexo; ///< The branching ratio of exotic (not invisible) Higgs decays.

    double dg1Z; ///< Independent contribution to aTGC.
    double dKappaga; ///< Independent contribution to aTGC.    
    double lambZ; ///< Independent contribution to aTGC.

    //  ----------------------------------------------------------------------------
    //  Internal parameters    
    //  ----------------------------------------------------------------------------

    double LambdaNP2; ///< The square of the new physics scale [GeV\f$^2\f$].

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

    double v2; ///< The square of the EW vev.
    double v2_over_LambdaNP2; ///< The ratio between the EW vev and the new physics scale, squared \f$v^2/\Lambda^2\f$.
    double aleMz; ///< The em constant at Mz. 
    double eeMz; ///< The em coupling at Mz.
    double eeMz2; ///< The em coupling squared (at Mz).
    double cW_tree; ///< The tree level values for the cosine of the weak angle.
    double sW_tree; ///< The tree level values for the sine of the weak angle.
    double cW2_tree; ///< The square of the tree level values for the cosine of the weak angle.
    double sW2_tree; ///< The square of the tree level values for the sine of the weak angle.

    double g1_tree; ///< The tree level value of the \f$U(1)_Y\f$ gauge coupling contant (at the \f$Z\f$ pole).
    double g2_tree; ///< The tree level value of the \f$SU(2)_L\f$ gauge coupling contant (at the \f$Z\f$ pole).
    double g3_tree; ///< The tree level value of the \f$SU(3)_c\f$ gauge coupling contant (at the \f$Z\f$ pole).

    double Mw_tree; ///< The tree level value of the \f$W\f$ boson mass.

    double lambdaH_tree; ///< The SM tree level value of the scalar quartic coupling in the potential.

    double gZvL; ///< The tree level value of the \f$Z\bar{\nu}\nu\f$ couplings in the SM.
    double gZlL, gZlR; ///< The tree level value of the \f$Z\ell^+\ell^-\f$ couplings in the SM.
    double gZuL, gZuR; ///< The tree level value of the \f$Z\bar{u}u\f$ couplings in the SM.
    double gZdL, gZdR; ///< The tree level value of the \f$Z\bar{d}d\f$ couplings in the SM.

    double UevL; ///< The tree level value of the \f$W^-\bar{\ell}\nu\f$ couplings in the SM. (Neglecting PMNS effects.)
    double VudL; ///< The tree level value of the \f$W^+\bar{u}d\f$ couplings in the SM. (Neglecting CKM effects.)

    double delta_ZZ; ///< Combination of dimension 6 coefficients modifying the \f$Z_\mu\f$ canonical field definition.
    double delta_AA; ///< Combination of dimension 6 coefficients modifying the \f$A_\mu\f$ canonical field definition.
    double delta_AZ; ///< Combination of dimension 6 coefficients modifying the \f$A_\mu\f$ canonical field definition.

    double delta_Z; ///< Combination of dimension 6 coefficients modifying the \f$Z_\mu\f$ canonical field definition for EWPO.
    double delta_A; ///< Combination of dimension 6 coefficients modifying the \f$A_\mu\f$ canonical field definition for EWPO.
    double delta_ZA; ///< Combination of dimension 6 coefficients modifying the \f$A_\mu\f$ canonical field definition for EWPO.

    double delta_h; ///< Combinations of dimension 6 coefficients modifying the \f$H\f$ canonical field definition.

    double dZH, dZH1, dZH2; ///< Higgs self-coupling contribution to the universal resummed Higgs wave function renormalization and combinations entering in the linear and quadratic pieces

    double cHSM; ///< Parameter to control the inclusion of modifications of SM parameters in selected Higgs processes.

    double cLHd6; ///< Parameter to control the inclusion of modifications of SM loops in Higgs processes due to dim 6 interactions.

    double cLH3d62; ///< Parameter to control the inclusion of modifications of SM loops in Higgs processes due to dim 6 interactions modifying the Higgs trilinear coupling (Quadratic terms).

    double cRGE; ///< Parameter to control the inclusion of log-enhanced contributions via RG effects. If activated then it takes the value multiplying the anomalous dimension: \f$-\log(\Lambda/\mu)/16 \pi^2 \Lambda^2\f$
    double cRGEon; ///< Another parameter to control the inclusion of log-enhanced contributions via RG effects. It multiplies the logarithmic contribution to some observables, so it takes values 0/1 to deactivate/activate the log terms.

    double cAsch, cWsch; ///< Parameters to control the SM EW input scheme: Alpha or MW.

    double Yuke, Yukmu, Yuktau; ///< SM lepton Yukawas
    double Yuku, Yukc, Yukt; ///< SM u-quark Yukawas
    double Yukd, Yuks, Yukb; ///< SM d-quark Yukawas

    // STXS ci: denote them as ai to separate from the normal conventions 
    double aiG, ai3G, ai2G;
    double aiT, aiH, aiWW, aiB, aiHW, aiHB;
    double aiA;
    double aiHQ, aipHQ, aiHL, aipHL, aiHu, aiHd, aiHe;
    double aiu, aiuG;

    // Cache variables
    double dGammaHTotR1, dGammaHTotR2, GammaHTotR; ///< NP contributions and Total to Higgs width ratio with SM

    double delta_GF; ///< The dimension 6 correction to the Fermi constant, as extracted from muon decay

    double delta_MZ; ///< The dimension 6 correction to Z mass Lagrangian parameter
    double delta_MW; ///< The dimension 6 correction to W mass Lagrangian parameter
    double delta_v; ///< The dimension 6 correction to the vev, as extracted from GF
    double delta_e; ///< The dimension 6 correction to the electric constant parameter
    double delta_sW2; ///< The dimension 6 correction to the weak mixing angle
    double delta_UgNC; ///< The dimension 6 universal correction to neutral current EW couplings
    double delta_QgNC; ///< The dimension 6 charge correction to neutral current EW couplings
    double delta_UgCC; ///< The dimension 6 universal correction to charged current EW couplings
    
    double delta_em; ///< The relative dimension 6 correction to the QED interaction vertex
    
    double C1Htotal; ///< The C1 coefficient controlling the H^3 corrections to the total Higgs width from the Higgs trilinear coupling.

    //AG:begin
    /**
     * @brief The \f$\mathcal{O}(\Lambda^{-4})\f$ dimension 6 correction to the Fermi constant.
     * @details
     * \f[
     * \dGf^{(2)} = \left( (C_{\varphi l}^{(3)11})^2 +(C_{\varphi l}^{(3)22})^2 +3C_{\varphi l}^{(3)11} C_{\varphi l}^{(3)22} \right) \frac{\vt^4 }{\Lambda^4}
     * - \left(C_{\varphi l}^{(3)11} + C_{\varphi l}^{(3)22}\right) \left( C_{ll}^{1221}+C_{ll}^{2112}\right) \frac{\vt^4 }{\Lambda^4}
     * + \frac{1}{4} \left( C_{ll}^{1221}+C_{ll}^{2112}\right)^2 \frac{\vt^4 }{\Lambda^4}
     * \f]
     */
    double delta_GF_2;

    /**
     * @brief The \f$\mathcal{O}(\Lambda^{-2})\f$ dimension 6 correction to the electromagnetic coupling.
     * @details
     * \f$\nc{\daEW}{\delta_{\alpha}}\f$
     * \f$\nc{\dMZZ}{\delta_{M_Z^2}}\f$
     * \f$\nc{\dGf}{\delta_{G_F}}\f$
     * \f[
     * \daEW^{(1)} = -2 s_W c_W \frac{C_{\varphi WB}\vt^2}{\Lambda^2}
     * \f]
     */
    double delta_ale;

    /**
     * @brief The \f$\mathcal{O}(\Lambda^{-4})\f$ dimension 6 correction to the electromagnetic coupling.
     * @details
     * \f[
     * \daEW^{(2)} = -2 s_W c_W \frac{C_{\varphi WB} \vt^2}{\Lambda^2} \left( \frac{\left(C_{\varphi W}+C_{\varphi B}\right)\vt^2}{\Lambda^2} 
     * + \frac{1}{4} \frac{C_{\varphi D}\vt^2}{\Lambda^2} + 3\dGf^{(1)}\right)
     * \f]
     */
    double delta_ale_2;

    /**
     * @brief The \f$\mathcal{O}(\Lambda^{-2})\f$ dimension 6 correction to the Z-boson mass squared.
     * @details
     * \f[
     * \dMZZ^{(1)} = 2 s_W c_W \frac{C_{\varphi WB}\vt^2}{\Lambda^2} + \frac{1}{2} \frac{C_{\varphi D}\vt^2}{\Lambda^2}
     * \f]
     */
    double delta_Mz2;

    /**
     * @brief The \f$\mathcal{O}(\Lambda^{-4})\f$ dimension 6 correction to the Z-boson mass squared.
     * @details
     * \f{eqnarray*}{
     * \dMZZ^{(2)} &=& \text{cAsch} \left[ \left(3 s_W c_W C_{\varphi WB}\frac{vt^2}{\Lambda^2} + \frac{1}{2} C_{\varphi D}\frac{\vt^2}{\Lambda^2}\right)\dGf^{(1)}
     * + 2 s_W c_W C_{\varphi WB}\left(C_{\varphi W}+C_{\varphi B}+\frac{3}{4}C_{\varphi D}\right)\frac{\vt^4}{\Lambda^4} + (1+4s_W^2 c_W^2)C_{\varphi WB}^2 \frac{\vt^4}{\Lambda^4} \right]\\
     * &+& \text{cWsch} \left[ \left(2s_W c_W C_{\varphi WB}\frac{\vt^2}{\Lambda^2} + \frac{1}{2} C_{\varphi D}\frac{\vt^2}{\Lambda^2} \right) \dGf^{(1)}
     * + (1+2 c_W^2 - 4 c_W^4) C_{\varphi WB}^2 \frac{\vt^4}{\Lambda^4} 
     * + 2 s_W c_W C_{\varphi WB} \left(C_{\varphi W} + C_{\varphi B} + \frac{1}{2}C_{\varphi D}\right) \frac{\vt^4}{\Lambda^4}
     * + \frac{(1-2 c_W^2)c_W}{2s_W} C_{\varphi WB} C_{\varphi D} \frac{\vt^4}{\Lambda^4} \right] 
     * \f}
     */
    double delta_Mz2_2;

    /**
     * @brief The \f$\mathcal{O}(\Lambda^{-2})\f$ dimension 6 correction to the \f$U(1)_Y\f$ gauge coupling, 
     * for the Alpha-Scheme (cAsch=1, cWsch=0) or the Mw-Scheme (cAsch=0, cWsch=1).
     * @details
     * \f$\nc{\daEW}{\delta_{\alpha}}\f$
     * \f$\nc{\dMZZ}{\delta_{M_Z^2}}\f$
     * \f$\nc{\dGf}{\delta_{G_F}}\f$
     * \f[
     * \gbL =  \text{cAsch} \left[\gbt \frac{c_W^2\daEW^{(1)}-s_W^2\left(\dMZZ^{(1)}+\dGf^{(1)}\right)}{2(-1+2 s_W^2)} \right] 
     * + \text{cWsch} \left[\gbt \frac{-\dMZZ^{(1)}-s_W^2\dGf^{(1)}}{2 s_W^2}\right]
     * \f]
     */
    double delta_g1;

    /**
     * @brief The \f$\mathcal{O}(\Lambda^{-4})\f$ dimension 6 correction to the \f$U(1)_Y\f$ gauge coupling.
     * @details
     * Alpha-Scheme:
     * \f{eqnarray*}{
     * \gbQ &=&  \gbt \frac{4 (-1+2 s_W^2)^2 \left( c_W^2 \daEW^{(2)} - s_W^2\big(\dGf^{(2)} + \dMZZ^{(2)}\big) \right)
     * + \big( -3+12s_W^2-19s_W^4+10s_W^6 \big) \left(\daEW^{(1)}\right)^2 
     * + s_W^4\big(-7+10s_W^2\big) \left((\dMZZ^{(1)})^2+(\dGf^{(1)})^2\right)}{8(-1+2s_W^2)^3} \\
     * &+& \gbt \frac{ s_W^2 (3-5s_W^2+2s_W^4) \left( \daEW^{(1)}\dMZZ^{(1)} + \daEW^{(1)}\dGf^{(1)}\right)
     * +  s_W^2 (-2+s_W^2+2s_W^4) ~ \dMZZ^{(1)} ~ \dGf^{(1)}}{4(-1+2s_W^2)^3}  
     * \f}
     * Mw-Scheme:
     * \f[
     * \gbQ = \gbt \left[ -\frac{\dMZZ^{(2)}}{2s_W^2} - \frac{\dGf^{(2)}}{2} 
     * - \frac{ (1-4 s_W^2)\left(\dMZZ^{(1)}\right)^2}{8s_W^4} + \frac{3\left(\dGf^{(1)}\right)^2}{8} + \frac{\dMZZ^{(1)} \dGf^{(1)}}{4s_W^2}\right]
     * \f]
     */
    double delta_g1_2;

    /**
     * @brief The \f$\mathcal{O}(\Lambda^{-2})\f$ dimension 6 correction to the \f$SU(2)_L\f$ gauge coupling, 
     * for the Alpha-Scheme (cAsch=1, cWsch=0) or the Mw-Scheme (cAsch=0, cWsch=1).
     * @details
     * \f[
     * \gwL =  \text{cAsch}\left[ \gwt\left( -\frac{\daEW^{(1)}}{2} - \frac{c_W^2(-\daEW^{(1)}+\dGf^{(1)}+\dMZZ^{(1)})}{2(c_W^2-s_W^2)} \right) \right]
     * + \text{cWsch} \left[ \gwt \left(-\frac{\dGf^{(1)}}{2}\right) \right]
     * \f]
     */
    double delta_g2;

    /**
     * @brief The \f$\mathcal{O}(\Lambda^{-4})\f$ dimension 6 correction to the \f$SU(2)_L\f$ gauge coupling.
     * @details
     * Alpha-Scheme:
     * \f{eqnarray*}{
     * \gwQ &=& \gwt \frac{ 4 (-1+2s_W^2)^2 \left( -s_W^2 \daEW^{(2)} + c_W^2\big(\dGf^{(2)} + \dMZZ^{(2)}\big) \right) 
     * + s_W^2 \big(4-11s_W^2+10s_W^4\big) \left(\daEW^{(1)}\right)^2 
     * + c_W^4 \big(-3+10s_W^2\big) \left((\dMZZ^{(1)})^2+(\dGf^{(1)})^2\right) }{8(-1+2s_W^2)^3}  \\
     * &+& \gwt \frac{  s_W^2 \big(-1-s_W^2+2s_W^4\big) \left( \daEW^{(1)}\dMZZ^{(1)} + \daEW^{(1)}\dGf^{(1)} \right)
     * + \big(-1+6s_W^2-7s_W^4+2s_W^6\big) ~ \dMZZ^{(1)} ~ \dGf^{(1)} }{4(-1+2s_W^2)^3}
     * \f}
     * Mw-Scheme:
     * \f[
     * \gwQ = \gwt \left[-\frac{\dGf^{(2)}}{2} + \frac{3\left(\dGf^{(1)}\right)^2}{8}\right]
     * \f]
     */
    double delta_g2_2;

    /**
     *@brief The tree level component of the matrix that transform the gauge field \f$W_{\mu}^3\f$ into \f$Z_{\mu}\f$.
     *@details 
     * \f[
     * \mathbb{\widetilde{X}}_{(WZ)} = c_W
     * \f]
     * where, \f$c_W\f$ is the tree level value for the cosine of the weak angle. 
     */
    double xWZ_tree;

    /**
     *@brief The tree level component of the matrix that transform the gauge field \f$B_{\mu}\f$ into \f$Z_{\mu}\f$.
     *@details 
     * \f[
     * \mathbb{\widetilde{X}}_{(BZ)} = - s_W
     * \f]
     * where, \f$s_W\f$ is the tree level value for the sine of the weak angle. 
     */
    double xBZ_tree;

    /**
     *@brief The \f$\mathcal{O}(\Lambda^{-2})\f$ dimension 6 correction to the component of the matrix that transform the gauge field \f$W_{\mu}^3\f$ into \f$Z_{\mu}\f$.
     *@details 
     * \f[
     * \Delta \mathbb{X}_{WZ}^{(1)} = s_W^2 c_W \left(\frac{\gwL}{\gwt}-\frac{\gbL}{\gbt}\right) + s_W^3 \frac{C_{\varphi WB}\vt^2}{\Lambda^2 }
     * \f]
     * where, \f$s_W\f$ and \f$c_W\f$ are the tree level values for the sine and cosine of the weak angle. 
     * Whereas \f$\gwt\f$ and \f$\gbt\f$ are tree level values of the \f$SU(2)_{L}\f$ and \f$U(1)_Y\f$ gauge couplings at the Z-pole, 
     * and \f$\Delta g^{(1)}\f$ their \f$\mathcal{O}(\Lambda^{-2})\f$ new physics corrections.
     */
    double delta_xWZ;

    /**
     *@brief The \f$\mathcal{O}(\Lambda^{-2})\f$ dimension 6 correction to the component of the matrix that transform the gauge field \f$B_{\mu}\f$ into \f$Z_{\mu}\f$.
     *@details 
     * \f[
     * \Delta \mathbb{X}_{BZ}^{(1)} = c_W^2 s_W \left(\frac{\gwL}{\gwt}-\frac{\gbL}{\gbt}\right) - c_W^3 \frac{C_{\varphi WB}\vt^2}{\Lambda^2 } 
     * \f]
     * where, \f$s_W\f$ and \f$c_W\f$ are the tree level values for the sine and cosine of the weak angle. 
     * Whereas \f$\gwt\f$ and \f$\gbt\f$ are tree level values of the \f$SU(2)_{L}\f$ and \f$U(1)_Y\f$ gauge couplings at the Z-pole, 
     * and \f$\Delta g^{(1)}\f$ their \f$\mathcal{O}(\Lambda^{-2})\f$ new physics corrections.
     */
    double delta_xBZ;

    /**
     *@brief The \f$\mathcal{O}(\Lambda^{-4})\f$ dimension 6 correction to the component of the matrix that transform the gauge field \f$W_{\mu}^3\f$ into \f$Z_{\mu}\f$.
     *@details 
     * \f{eqnarray*}{
     * \Delta \mathbb{X}_{WZ}^{(2)} &=& -\frac{3 s_W^2 c_W^3}{2}\left(\frac{\gwL}{\gwt}\right)^2 
     * + \frac{s_W^2 c_W(2 s_W^2- c_W^2)}{2}\left(\frac{\gbL}{\gbt}\right)^2
     * - s_W^2 c_W(s_W^2-2 c_W^2)\frac{\gwL}{\gwt}\frac{\gbL}{\gbt} 
     * + s_W^2 c_W\left(\frac{\gwQ}{\gwt}-\frac{\gbQ}{\gbt}\right) \\
     * &-& 3 s_W^3 c_W^2 \left(\frac{\gwL}{\gwt}-\frac{\gbL}{\gbt}\right)\frac{C_{\varphi WB}\vt^2}{\Lambda^2 }
     * -c_W(s_W^4-3s_W^2 c_W^2-c_W^4)\frac{C_{\varphi WB}\vt^4}{\Lambda^4 }
     * + s_W^3\left(C_{\varphi WB}(C_{\varphi W}+C_{\varphi B})\frac{\vt^4}{\Lambda^4} + \delta_{G_F}^{(1)} \frac{C_{\varphi WB} \vt^2}{\Lambda^2}  \right)
     * \f}
     * where, \f$s_W\f$ and \f$c_W\f$ are the tree level values for the sine and cosine of the weak angle. 
     * Whereas \f$\gwt\f$ and \f$\gbt\f$ are tree level values of the \f$SU(2)_{L}\f$ and \f$U(1)_Y\f$ gauge couplings at the Z-pole;
     * \f$\Delta g^{(1)}\f$ and \f$\Delta g^{(2)}\f$ their respective \f$\mathcal{O}(\Lambda^{-2})\f$ and \f$\mathcal{O}(\Lambda^{-4})\f$ new physics corrections.
     */
    double delta_xWZ_2;

    /**
     *@brief The \f$\mathcal{O}(\Lambda^{-4})\f$ dimension 6 correction to the component of the matrix that transform the gauge field \f$B_{\mu}\f$ into \f$Z_{\mu}\f$.
     *@details 
     * \f{eqnarray*}{
     * \Delta \mathbb{X}_{BZ}^{(2)} &=& \frac{3 s_W^3 c_W^2}{2}\left(\frac{\gbL}{\gbt^2}\right)^2 
     * + \frac{s_W c_W^2(s_W^2-2c_W^2)}{2}\left(\frac{\gwL}{\gwt}\right)^2 
     * - s_W c_W(2s_W^2-c_W^2)\frac{\gwL}{\gwt}\frac{\gbL}{\gbt}
     * + s_W c_W^2\left(\frac{\gwQ}{\gwt}-\frac{\gbQ}{\gbt}\right) \\
     * &-& 3s_W^2 c_W^3\left(\frac{\gwL}{\gwt}-\frac{\gbL}{\gbt}\right)\frac{C_{\varphi WB}\vt^2}{\Lambda^2 }
     * - \frac{1}{2}s_W(s_W^4+3s_W^2 c_W^2- c_W^4)\frac{C_{\varphi WB}^2\vt^4}{\Lambda^4 }
     * + c_W^3\left(C_{\varphi WB}(C_{\varphi W}+C_{\varphi B})\frac{\vt^4}{\Lambda^4} + \delta_{G_F}^{(1)} \frac{C_{\varphi WB} \vt^2}{\Lambda^2}  \right)
     * \f}
     * where, \f$s_W\f$ and \f$c_W\f$ are the tree level values for the sine and cosine of the weak angle. 
     * Whereas \f$\gwt\f$ and \f$\gbt\f$ are tree level values of the \f$SU(2)_{L}\f$ and \f$U(1)_Y\f$ gauge couplings at the Z-pole;
     * \f$\Delta g^{(1)}\f$ and \f$\Delta g^{(2)}\f$ their respective \f$\mathcal{O}(\Lambda^{-2})\f$ and \f$\mathcal{O}(\Lambda^{-4})\f$ new physics corrections.
     */
    double delta_xBZ_2;
    //AG:end

    /**
     * @brief The diagonal entry of the dimension-6 operator coefficient \f$C_{HL,HQ}^{(1)}\f$ corresponding to particle F.
     * @param[in] F a lepton or quark
     * @return @f$(\f$C_{HF}^{(1)})_{FF}\f$@f$
     */
    const double CHF1_diag(const Particle F) const;

    /**
     * @brief The diagonal entry of the dimension-6 operator coefficient \f$C_{HL,HQ}^{(3)}\f$ corresponding to particle F.
     * @param[in] F a lepton or quark
     * @return @f$(\f$C_{HF}^{(3)})_{FF}\f$@f$
     */
    const double CHF3_diag(const Particle F) const;

    /**
     * @brief The diagonal entry of the dimension-6 operator coefficient \f$C_{HE,HU,HD}\f$ corresponding to particle f.
     * @param[in] f a lepton or quark
     * @return @f$(\f$C_{Hf})_{ff}\f$@f$
     */
    const double CHf_diag(const Particle f) const;

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
    bool FlagMWinput; ///< A boolean for the model flag %MWinput. 

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

    gsl_integration_cquad_workspace * w_WW; /**< Gsl integral variable */

};

#endif /* NPSMEFTD6_H */
