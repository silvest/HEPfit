
/* 
 * File:   NPSMEFTd6General.h
 * Author: miralles
 *
 * Created on 12 September 2022, 14:17
 */

#ifndef NPSMEFTD6GENERAL_H
#define NPSMEFTD6GENERAL_H

#include "gslpp.h"
#include "NPbase.h"
#include "NPSMEFTd6GeneralMatching.h"
#include <string.h>
#include <stdexcept>
#include <gsl/gsl_integration.h>
#include <RGESolver.h>

/**
 * @class NPSMEFTd6General
 * @brief A model class for new physics in the form of the
 * dimension-six effective Lagrangian.
 * @ingroup NewPhysics
 * @author HEPfit Collaboration
 * @copyright GNU General Public License
 * @details This is a Model class containing parameters and functions
 * associated with the general dimension-six effective Lagrangian.
 * (Use the model "NPSMEFTd6_LFU_QFU" to assume lepton and quark flavour universality)
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
 * @anchor NPSMEFTd6GeneralInitialization
 * <h3>Initialization</h3>
 *
 * After creating an instance of the current class with the constructor
 * NPSMEFTd6General(), it is required to call the initialization method
 * InitializeModel().
 * In the Monte Carlo run, the constructor as well as the initialization
 * method are called in InputParser::ReadParameters().
 *
 *
 * @anchor NPSMEFTd6GeneralParameters
 * <h3>%Model parameters</h3>
 *
 * The model parameters of %NPSMEFTd6General are summarized below:
 * <table class="model">
 * <tr>
 *   <th>Label</th>
 *   <th>LaTeX symbol</th>
 *   <th>Description</th>
 * </tr>
 * <tr>
 *   <td class="mod_name">%g1_LNP </td>
 *   <td class="mod_symb">\f$g_1(\Lambda_\mathrm{NP})\f$</td>
 *   <td class="mod_desc">The hypercharge coupling constant at the New Physics scale. </td>
 * </tr>
 * <tr>
 *   <td class="mod_name">%g2_LNP </td>
 *   <td class="mod_symb">\f$g_2(\Lambda_\mathrm{NP})\f$</td>
 *   <td class="mod_desc">The SU(2)@f$_L@f$ coupling constant at the New Physics scale. </td>
 * </tr>
 * <tr>
 *   <td class="mod_name">%g1_LNP </td>
 *   <td class="mod_symb">\f$g_3(\Lambda_\mathrm{NP})\f$</td>
 *   <td class="mod_desc">The QCD coupling constant at the New Physics scale. </td>
 * </tr>
 * <tr>
 *   <td class="mod_name">%lambdaH_LNP </td>
 *   <td class="mod_symb">\f$\lambda_H(\Lambda_\mathrm{NP})\f$</td>
 *   <td class="mod_desc">The coefficient of the @f$(H^\dagger H)^2@f$ term at the New Physics scale. </td>
 * </tr>
 * <tr>
 *   <td class="mod_name">%muH_LNP </td>
 *   <td class="mod_symb">\f$\mu_H(\Lambda_\mathrm{NP})\f$</td>
 *   <td class="mod_desc">The square root of the coefficient of the @f$(H^\dagger H)@f$ term at the New Physics scale. </td>
 * </tr>
 * <tr>
 *   <td class="mod_name">%Ye_LNP </td>
 *   <td class="mod_symb">\f$Y_e(\Lambda_\mathrm{NP})\f$</td>
 *   <td class="mod_desc">The electron Yukawa coupling at the New Physics scale. </td>
 * </tr>
 * <tr>
 *   <td class="mod_name">%Ymu_LNP </td>
 *   <td class="mod_symb">\f$Y_\mu(\Lambda_\mathrm{NP})\f$</td>
 *   <td class="mod_desc">The muon Yukawa coupling at the New Physics scale. </td>
 * </tr>
 * <tr>
 *   <td class="mod_name">%Ytau_LNP </td>
 *   <td class="mod_symb">\f$Y_\tau(\Lambda_\mathrm{NP})\f$</td>
 *   <td class="mod_desc">The tau Yukawa coupling at the New Physics scale. </td>
 * </tr>
 * <tr>
 *   <td class="mod_name">%Yu_LNP </td>
 *   <td class="mod_symb">\f$Y_u(\Lambda_\mathrm{NP})\f$</td>
 *   <td class="mod_desc">The up quark Yukawa coupling at the New Physics scale in the up basis. </td>
 * </tr>
 * <tr>
 *   <td class="mod_name">%Yc_LNP </td>
 *   <td class="mod_symb">\f$Y_c(\Lambda_\mathrm{NP})\f$</td>
 *   <td class="mod_desc">The charm quark Yukawa coupling at the New Physics scale in the up basis. </td>
 * </tr>
 * <tr>
 *   <td class="mod_name">%Yt_LNP </td>
 *   <td class="mod_symb">\f$Y_t(\Lambda_\mathrm{NP})\f$</td>
 *   <td class="mod_desc">The top quark Yukawa coupling at the New Physics scale in the up basis. </td>
 * </tr>
 * <tr>
 *   <td class="mod_name">%Yd_LNP </td>
 *   <td class="mod_symb">\f$Y_d(\Lambda_\mathrm{NP})\f$</td>
 *   <td class="mod_desc">The down quark Yukawa coupling at the New Physics scale in the down basis. </td>
 * </tr>
 * <tr>
 *   <td class="mod_name">%Ys_LNP </td>
 *   <td class="mod_symb">\f$Y_s(\Lambda_\mathrm{NP})\f$</td>
 *   <td class="mod_desc">The strange quark Yukawa coupling at the New Physics scale in the down basis. </td>
 * </tr>
 * <tr>
 *   <td class="mod_name">%Yb_LNP </td>
 *   <td class="mod_symb">\f$Y_b(\Lambda_\mathrm{NP})\f$</td>
 *   <td class="mod_desc">The bottom quark Yukawa coupling at the New Physics scale in the down basis. </td>
 * </tr>
 * <tr>
 *   <td class="mod_name">%s12CKM_LNP </td>
 *   <td class="mod_symb">\f$s_{12}(\Lambda_\mathrm{NP})\f$</td>
 *   <td class="mod_desc">The sine of the @f$\theta_{12}@f$ angle parameterizing the relative orientation of up and down Yukawa couplings at the New Physics scale. </td>
 * </tr>
 * <tr>
 *   <td class="mod_name">%s13CKM_LNP </td>
 *   <td class="mod_symb">\f$s_{13}(\Lambda_\mathrm{NP})\f$</td>
 *   <td class="mod_desc">The sine of the @f$\theta_{13}@f$ angle parameterizing the relative orientation of up and down Yukawa couplings at the New Physics scale. </td>
 * </tr>
 * <tr>
 *   <td class="mod_name">%s23CKM_LNP </td>
 *   <td class="mod_symb">\f$s_{23}(\Lambda_\mathrm{NP})\f$</td>
 *   <td class="mod_desc">The sine of the @f$\theta_{23}@f$ angle parameterizing the relative orientation of up and down Yukawa couplings at the New Physics scale. </td>
 * </tr>
 * <tr>
 *   <td class="mod_name">%dCKM_LNP </td>
 *   <td class="mod_symb">\f$\delta(\Lambda_\mathrm{NP})\f$</td>
 *   <td class="mod_desc">The CP-violating phase of the matrix describing the relative orientation of up and down Yukawa couplings at the New Physics scale. </td>
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
 *   <td class="mod_name">%CHWB </td>
 *   <td class="mod_symb">\f$C_{HWB} \f$</td>
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
 *   <td class="mod_name">%CHl1_kk, CHl1_klr, CHl1_kli </td>
 *   <td class="mod_symb">\f$ (C_{HL}^{(1)})_{kk}, \mbox{Re}\big[(C_{HL}^{(1)})_{kl}\big], \mbox{Im}\big[(C_{HL}^{(1)})_{kl}\big] \f$</td>
 *   <td class="mod_desc">The real and imaginary parts of the coefficient of the operator
 *   \f$({\cal O}_{HL}^{(1)})_{ij} =i\big(H^\dagger \overset{\leftrightarrow}{D}_\mu H\big)
 *   \big(\overline{L^i}\,\gamma^\mu L^j\big)\f$, for \f$i,j=1,2,3\f$. </td>
 * </tr>
 * <tr>
 *   <td class="mod_name">%CHl3_kk, CHl3_klr, CHl3_kli </td>
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
 *   <td class="mod_name">%CHq1_kk, CHq1_klr, CHq1_kli </td>
 *   <td class="mod_symb">\f$ (C_{HQ}^{(1)})_{kk}, \mbox{Re}\big[(C_{HQ}^{(1)})_{kl}\big], \mbox{Im}\big[(C_{HQ}^{(1)})_{kl}\big] \f$</td>
 *   <td class="mod_desc">The real and imaginary parts of the coefficient of the operator
 *   \f$({\cal O}_{HQ}^{(1)})_{ij} =i\big(H^\dagger \overset{\leftrightarrow}{D}_\mu H\big)
 *   \big(\overline{Q^i}\,\gamma^\mu Q^j\big)\f$, for \f$i,j=1,2,3\f$. </td>
 * </tr>
 * <tr>
 *   <td class="mod_name">%CHq3_kk, CHq3_klr, CHq3_kli </td>
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
 *   <td class="mod_name">%Cll_1221, Cll_2112 </td>
 *   <td class="mod_symb">\f$(C_{LL})_{1221,2112}\f$</td>
 *   <td class="mod_desc">The coefficient of the operator
 *   \f$({\cal O}_{LL})_{ijkl}=\big(\overline{L^i}\,\gamma^\mu L^j\big)
 *   \big(\overline{L^k}\,\gamma_\mu L^l\big)\f$, for \f$ijkl=1221,2112\f$. </td>
 * </tr>
 * <tr>
 *   <td class="mod_name">%Clq1 </td>
 *   <td class="mod_symb">\f$C_{LQ}^{(1)}\f$</td>
 *   <td class="mod_desc">The coefficient of the operator
 *   \f$({\cal O}_{LQ}^{(1)})_{ijkl}=\big(\overline{L^i}\,\gamma^\mu L^j\big)
 *   \big(\overline{Q^k}\,\gamma_\mu Q^l\big)\f$. </td>
 * </tr>
 * <tr>
 *   <td class="mod_name">%Clq3 </td>
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
 *   <td class="mod_name">%Cle </td>
 *   <td class="mod_symb">\f$C_{LE}\f$</td>
 *   <td class="mod_desc">The coefficient of the operator
 *   \f$({\cal O}_{LE})_{ijkl}=\big(\overline{L^i}\,\gamma^\mu L^j\big)
 *   \big(\overline{E^k}\,\gamma_\mu E^l\big)\f$. </td>
 * </tr>
 * <tr>
 *   <td class="mod_name">%Clu </td>
 *   <td class="mod_symb">\f$C_{LU}\f$</td>
 *   <td class="mod_desc">The coefficient of the operator
 *   \f$({\cal O}_{LU})_{ijkl}=\big(\overline{L^i}\,\gamma^\mu L^j\big)
 *   \big(\overline{U^k}\,\gamma_\mu U^l\big)\f$. </td>
 * </tr>
 * <tr>
 *   <td class="mod_name">%Cld </td>
 *   <td class="mod_symb">\f$C_{LD}\f$</td>
 *   <td class="mod_desc">The coefficient of the operator
 *   \f$({\cal O}_{LD})_{ijkl}=\big(\overline{L^i}\,\gamma^\mu L^j\big)
 *   \big(\overline{D^k}\,\gamma_\mu D^l\big)\f$. </td>
 * </tr>
 * <tr>
 *   <td class="mod_name">%Cqe </td>
 *   <td class="mod_symb">\f$C_{QE}\f$</td>
 *   <td class="mod_desc">The coefficient of the operator
 *   \f$({\cal O}_{QE})_{ijkl}=\big(\overline{Q^i}\,\gamma^\mu Q^j\big)
 *   \big(\overline{E^k}\,\gamma_\mu E^l\big)\f$. </td>
 * </tr> 
 * <tr>
 *   <td class="mod_name">%Cqq1_1133, Cqq1_1331, Cqq1_3333 </td>
 *   <td class="mod_symb">\f$C_{QQ}^{(1)}\f$</td>
 *   <td class="mod_desc">The coefficient of the operator
 *   \f$({\cal O}_{QQ}^{(1)})_{ijkl}=\big(\overline{Q^i}\,\gamma^\mu Q^j\big)
 *   \big(\overline{Q^k}\,\gamma_\mu Q^l\big)\f$. </td>
 * </tr>
 * <tr>
 *   <td class="mod_name">%Cqq3_1133, Cqq3_1331, Cqq3_3333 </td>
 *   <td class="mod_symb">\f$C_{QQ}^{(3)}\f$</td>
 *   <td class="mod_desc">The coefficient of the operator
 *   \f$({\cal O}_{QQ}^{(3)})_{ijkl}=\big(\overline{Q^i}\,\gamma^\mu \sigma_a Q^j\big)
 *   \big(\overline{Q^k}\,\gamma_\mu \sigma_a Q^l\big)\f$. </td>
 * </tr>
 * <tr>
 *   <td class="mod_name">%Cuu_1133, Cuu_1331, Cuu_3333 </td>
 *   <td class="mod_symb">\f$C_{uu}\f$</td>
 *   <td class="mod_desc">The coefficient of the operator
 *   \f$({\cal O}_{uu})_{ijkl}=\big(\overline{u^i}\,\gamma^\mu u^j\big)
 *   \big(\overline{u^k}\,\gamma_\mu u^l\big)\f$. </td>
 * </tr>
 * <tr>
 *   <td class="mod_name">%Cud1_3311, Cud1_3333 </td>
 *   <td class="mod_symb">\f$C_{ud}^{(1)}\f$</td>
 *   <td class="mod_desc">The coefficient of the operator
 *   \f$({\cal O}_{ud}^{(1)})_{ijkl}=\big(\overline{u^i}\,\gamma^\mu u^j\big)
 *   \big(\overline{d^k}\,\gamma_\mu d^l\big)\f$. </td>
 * </tr>
 * <tr>
 *   <td class="mod_name">%Cud8_3311, Cud8_3333 </td>
 *   <td class="mod_symb">\f$C_{ud}^{(8)}\f$</td>
 *   <td class="mod_desc">The coefficient of the operator
 *   \f$({\cal O}_{ud}^{(8)})_{ijkl}=\big(\overline{u^i}\,\gamma^\mu T^A u^j\big)
 *   \big(\overline{d^k}\,\gamma_\mu T^A d^l\big)\f$. </td>
 * </tr>
 * <tr>
 *   <td class="mod_name">%Cqu1_1133, Cqu1_3311, Cqu1_3333 </td>
 *   <td class="mod_symb">\f$C_{Qu}^{(1)}\f$</td>
 *   <td class="mod_desc">The coefficient of the operator
 *   \f$({\cal O}_{Qu}^{(1)})_{ijkl}=\big(\overline{Q^i}\,\gamma^\mu Q^j\big)
 *   \big(\overline{u^k}\,\gamma_\mu u^l\big)\f$. </td>
 * </tr>
 * <tr>
 *   <td class="mod_name">%Cqu8_1133, Cqu8_3311, Cqu8_3333 </td>
 *   <td class="mod_symb">\f$C_{Qu}^{(8)}\f$</td>
 *   <td class="mod_desc">The coefficient of the operator
 *   \f$({\cal O}_{Qu}^{(8)})_{ijkl}=\big(\overline{Q^i}\,\gamma^\mu T^A Q^j\big)
 *   \big(\overline{u^k}\,\gamma_\mu T^A u^l\big)\f$. </td>
 * </tr>
 * <tr>
 *   <td class="mod_name">%Cqd1_3311, Cqd1_3333 </td>
 *   <td class="mod_symb">\f$C_{Qd}^{(1)}\f$</td>
 *   <td class="mod_desc">The coefficient of the operator
 *   \f$({\cal O}_{Qd}^{(1)})_{ijkl}=\big(\overline{Q^i}\,\gamma^\mu Q^j\big)
 *   \big(\overline{d^k}\,\gamma_\mu d^l\big)\f$. </td>
 * </tr>
 * <tr>
 *   <td class="mod_name">%Cqd8_3311, Cqd8_3333 </td>
 *   <td class="mod_symb">\f$C_{Qd}^{(8)}\f$</td>
 *   <td class="mod_desc">The coefficient of the operator
 *   \f$({\cal O}_{Qd}^{(8)})_{ijkl}=\big(\overline{Q^i}\,\gamma^\mu T^A Q^j\big)
 *   \big(\overline{d^k}\,\gamma_\mu T^A d^l\big)\f$. </td>
 * </tr> 
 * <tr>
 *   <td class="mod_name">%Cquqd1_3333 </td>
 *   <td class="mod_symb">\f$C_{QuQd}^{(1)}\f$</td>
 *   <td class="mod_desc">The coefficient of the operator
 *   \f$({\cal O}_{QuQd}^{(1)})_{ijkl}=\big(\overline{Q^i}\, u^j\big) i\sigma_2
 *   \big(\overline{Q^k}^T\, d^l\big)\f$. </td>
 * </tr>
 * <tr>
 *   <td class="mod_name">%Cquqd8_3333 </td>
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
 * @anchor NPSMEFTd6GeneralFlags
 * <h3>%Model flags</h3>
 *
 * The Flags of NPSMEFTd6General are summarized below:
 * <table class="model">
 * <tr>
 *   <th>Label</th>
 *   <th>Value</th>
 *   <th>Description</th>
 * </tr>
 * <tr>
 *   <td class="mod_name">%MWinput</td>
 *   <td class="mod_valu">TRUE&nbsp;/&nbsp;<b>FALSE</b></td>
 *   <td class="mod_desc">This auxiliary flag is used for setting the W mass as a SM input, 
 *   instead of the electromagnetic constant parameter %dAle5Mz. 
 *   The default value is FALSE. This flag will override the flag of the same name in StandardModel</td>
 * </tr> 
 * <tr>
 *   <td class="mod_name">%QuadraticTerms</td>
 *   <td class="mod_valu">TRUE&nbsp;/&nbsp;<b>FALSE</b></td>
 *   <td class="mod_desc">This flag is set to TRUE if the quadratic terms in Higgs cross sections and widths are switched on.
 *   The default value is FALSE; new physics contributions are linearized.</td>
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
 *   <td class="mod_name">%RGEci</td>
 *   <td class="mod_valu">TRUE&nbsp;/&nbsp;<b>FALSE</b></td>
 *   <td class="mod_desc">This flag is set to TRUE if including the SMEFT RGE effects.
 *   The default value is TRUE.</td>
 * </tr>
 * <tr>
 *   <td class="mod_name">%finiteNLO</td>
 *   <td class="mod_valu">TRUE&nbsp;/&nbsp;<b>FALSE</b></td>
 *   <td class="mod_desc">This flag is set to TRUE if including SMEFT finite NLO terms (where available).
 *   The default value is FALSE.</td>
 * </tr>
 * <tr>
 *   <td class="mod_name">%matchLEFT</td>
 *   <td class="mod_valu">TRUE&nbsp;/&nbsp;<b>FALSE</b></td>
 *   <td class="mod_desc">This flag is set to TRUE if the matching between SMEFT and LEFT needs to be computed.
 *   The default value is TRUE.</td>
 * </tr>  
 * </table>
 *
 *
 * @anchor NPSMEFTd6GeneralFunctions
 * <h3>Important member functions</h3>
 *
 * See the base classes of the current class.
 *
 */
class NPSMEFTd6General : public NPbase {
public:

    /**
     *　@brief The number of the model parameters in %NPSMEFTd6General (including the 18 parameters needed for the SM and 79 auxiliary parameters). 
     */
    //static const int NNPSMEFTd6GeneralVars = 2708 + 18;
    static const int NNPSMEFTd6GeneralVars = 2708-208 + 79 ;

    
    /**
     * @brief A string array containing the labels of the model parameters in
     * %NPSMEFTd6General.
     */
    static const std::string NPSMEFTd6GeneralVars[NNPSMEFTd6GeneralVars];

    /**
     * @brief Constructor.
     */
    NPSMEFTd6General();
    
    
    /**
     * @brief A method to initialize the model parameters.
     * @param[in] DPars a map of the parameters that are being updated in the Monte Carlo run
     * (including parameters that are varied and those that are held constant)
     * @return a boolean that is true if the execution is successful
     */
    virtual bool Init(const std::map<std::string, double>& DPars);
    

    /**
     * @brief The pre-update method for %NPSMEFTd6General.
     * @details This method runs all the procedures that are need to be executed
     * before the model is successfully updated.
     * @return a boolean that is true if the execution is successful
     */
    virtual bool PreUpdate();

    /**
     * @brief The post-update method for %NPSMEFTd6General.
     * @details This method runs all the procedures that are need to be executed
     * after the model is successfully updated.
     * @return a boolean that is true if the execution is successful
     */
    virtual bool PostUpdate();

    
    //I think that for this model, with so many parameters, it's better to include only those we're really interested in
    //Since all the parameters are initialized to zero we can just skip the CheckParameters for this model
    /**
     * @brief A method to check if all the mandatory parameters for %NPSMEFTd6General
     * have been provided in model initialization.
     * @param[in] DPars a map of the parameters that are being updated in the Monte Carlo run
     * (including parameters that are varied and those that are held constant)
     * @return a boolean that is true if the execution is successful
     */
    //virtual bool CheckParameters(const std::map<std::string, double>& DPars);

    /**
     * @brief A method to set a flag of %NPSMEFTd6General.
     * @param[in] name name of a model flag
     * @param[in] value the boolean to be assigned to the flag specified by name
     * @return a boolean that is true if the execution is successful
     */
    virtual bool setFlag(const std::string name, const bool value);

    /**
     * @brief A method to set a flag of %NPSMEFTd6General.
     * @param[in] name name of a model flag
     * @param[in] value the value to be assigned to the flag specified by name
     * @return a boolean that is true if the execution is successful
     */
    virtual bool setFlagStr(const std::string name, const std::string value);

    /**
     * @brief A method to get the Matching object for this model.
     * @return The matching object for this model
     */
    virtual NPSMEFTd6GeneralMatching& getMatching() const {
        return NPSMEFTd6GM.getObj();
    }

    ////////////////////////////////////////////////////////////////////////

    /**
     * @brief Return Lambda_NP
     * @return \f$ Lambda_NP \f$ 
     */
    double getLambda_NP() const {
        return Lambda_NP;
    }
    
    
    /**
     * @brief Return CG_LNP
     * @return \f$ CG_LNP \f$ 
     */
    double getCG_LNP() const {
        return CG_LNP;
    }
    
    //////////////////////////////////////////////////////////////////////// 
    // Functions to facilitate the change of scheme (at tree level)
    
    /**
     * @brief Difference in prediction in @f$\alpha@f$ scheme and W mass scheme, 
     * computed from observable in @f$\alpha@f$ scheme. 
     * Difference at tree level.
     * @details The difference @f$\Delta O@f$ is defined as
     * @f[
     * \Delta O = O_{\alpha} - O_{W} \,,
     * @f]
     * where @f$O@f$ is a given observable in the SMEFT.
     * @param[in] dOSMdalpha the derivative with respect to @f$\alpha@f$ of 
     * the SM tree-level prediction of @f$O@f$ (in the @f$\alpha@f$ scheme)
     * @param[in] mu the RG scale associated to the observable
     * @return The tree level difference between schemes @f$\Delta O@f$
     */
    virtual const double  DeltaOalphtoW(const double dOSMdalpha, const double mu) const;
    
    /**
     * @brief Difference in prediction in @f$\alpha@f$ scheme and W mass scheme, 
     * computed from observable in W mass scheme. 
     * Difference at tree level.
     * @details The difference @f$\Delta O@f$ is defined as
     * @f[
     * \Delta O = O_{\alpha} - O_{W} \,,
     * @f]
     * where @f$O@f$ is a given observable in the SMEFT.
     * @param[in] dOSMdMW the derivative with respect to @f$M_W@f$ of 
     * the SM tree-level prediction of @f$O@f$ (in the @f$M_W@f$ scheme)
     * @param[in] mu the RG scale associated to the observable
     * @return The tree level difference between schemes @f$\Delta O@f$
     */
    virtual const double  DeltaOWtoalph(const double dOSMdMW, const double mu) const;

    ////////////////////////////////////////////////////////////////////////

    //test function to test speed of code
    //double test_direct() const;

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
    virtual const double  DeltaGF() const;
    
    ////////////////////////////////////////////////////////////////////////
    // Functions to compute quantities involved in indirect corrections, depending on RG scale
    
    /**
     * @brief Correction to photon WF.
     * @details 
     * @param[in] mu the RG scale associated to the observable where the function is used
     * @return @f$\delta A(\mu)@f$
     */
    virtual const double  del_A_mu(const double mu) const;
    
    /**
     * @brief Correction to Z WF.
     * @details 
     * @param[in] mu the RG scale associated to the observable where the function is used
     * @return @f$\delta Z(\mu)@f$
     */
    virtual const double  del_Z_mu(const double mu) const;
    
    /**
     * @brief Correction to Z-A mixing.
     * @details 
     * @param[in] mu the RG scale associated to the observable where the function is used
     * @return @f$\delta ZA(\mu)@f$
     */
    virtual const double  del_ZA_mu(const double mu) const;
    
    
    /**
     * @brief Correction to electric charge.
     * @details 
     * @param[in] mu the RG scale associated to the observable where the function is used
     * @return @f$\delta e(\mu)@f$
     */
    virtual const double  del_e_mu(const double mu) const;
    
    /**
     * @brief Correction to (sin squared of) weak mixing angle.
     * @details 
     * @param[in] mu the RG scale associated to the observable where the function is used
     * @return @f$\delta s_W^2(\mu)@f$
     */
    virtual const double  del_sW2_mu(const double mu) const;
    
    /**
     * @brief Universal indirect correction to EW neutral currents.
     * @details 
     * @param[in] mu the RG scale associated to the observable where the function is used
     * @return @f$\delta g_{NC}^{U}(\mu)@f$
     */
    virtual const double  delU_gNC(const double mu) const;
    
    /**
     * @brief Separate, charge-proportional, indirect correction to EW neutral currents.
     * @details 
     * @param[in] mu the RG scale associated to the observable where the function is used
     * @return @f$\delta g_{NC}^{Q}(\mu)@f$
     */
    virtual const double  delQ_gNC(const double mu) const;
    
    /**
     * @brief Universal indirect correction to EW charged currents.
     * @details 
     * @param[in] mu the RG scale associated to the observable where the function is used
     * @return @f$\delta g_{CC}^{U}(\mu)@f$
     */
    virtual const double  delU_gCC(const double mu) const;
    
    ///////////////////////////// Oblique parameters ////////////////////////////

    /**
     * @brief The oblique parameter \f$S\f$.
     * (Simplified implementation. Contribution only from @f$O_{HWB}@f$.)
     * @return the value of @f$S@f$
     */
    virtual const double  obliqueS() const;

    /**
     * @brief The oblique parameter \f$T\f$.
     * (Simplified implementation. Contribution only from @f$O_{HD}@f$.)
     * @return the value of @f$T@f$
     */
    virtual const double  obliqueT() const;

    /**
     * @brief The oblique parameter \f$U\f$.
     * @return the value of @f$U@f$
     */
    virtual const double  obliqueU() const;

    /**
     * @brief The oblique parameter \f$W\f$.
     * (Simplified implementation. Contribution only from @f$O_{2W}@f$.)
     * @return the value of @f$W@f$
     */
    virtual const double  obliqueW() const;

    /**
     * @brief The oblique parameter \f$Y\f$.
     * (Simplified implementation. Contribution only from @f$O_{2B}@f$.)
     * @return the value of @f$Y@f$
     */
    virtual const double  obliqueY() const;


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
     * @return @f$\alpha(M_Z^2)@f$
     */
    virtual const double  alphaMz() const;
    
    /////////////////////// EWPO: W mass and widths /////////////////////////

    /**
     * @brief The mass of the @f$W@f$ boson, @f$M_W@f$.
     * @return @f$M_W@f$ in GeV
     */
    virtual const double  Mw() const;

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
    virtual const double  GammaW(const Particle fi, const Particle fj) const;

    /**
     * @brief The new physics contribution to the total decay width of the @f$W@f$ boson, @f$\delta \Gamma_W@f$.
     * @return @f$\delta \Gamma_W@f$ in GeV
     */
    virtual const double deltaGamma_W() const;

    /**
     * @brief The total width of the @f$W@f$ boson, @f$\Gamma_W@f$.
     * @return @f$\Gamma_W@f$ in GeV
     */
    virtual const double  GammaW() const;

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
    
    /////////////////////// EWPO: Z-pole ///////////////////////// 
    
    
    /**
     * @brief The new physics contribution to the decay width of the @f$Z@f$ boson into a given fermion pair, @f$\delta \Gamma_Z^{f}@f$.
     * @param[in] f a lepton or quark
     * @return @f$\delta \Gamma_Z^{f}@f$ in GeV
     */
    virtual const double deltaGamma_Zf(const Particle f) const;
    
    /**
     * @brief The decay width of the @f$Z@f$ boson into a given fermion pair, @f$\Gamma_Z^{f}@f$.
     * @details
     * \f[
     * \Gamma_Z^{f} = \Gamma_{Z,f}^{SM} + \Delta \Gamma_{Z,f}^{(1)} + \Delta \Gamma_{Z,f}^{(2)}
     * \f]
     * @param[in] f a lepton or quark
     * @return @f$\Gamma_Z^{f}@f$ in GeV, including SM plus @f$\mathcal{O}(\Lambda^{-2})@f$ and @f$\mathcal{O}(\Lambda^{-4})@f$ NP contributions
     *
     * @attention This function is applicable only to the NP model classes that
     * are inherited from NPbase.
     */
    virtual const double Gamma_Zf(const Particle f) const;
    
    /**
     * @brief The new physics contribution to the total decay width of the @f$Z@f$ boson, @f$\delta \Gamma_Z@f$.
     * @return @f$\delta \Gamma_Z@f$ in GeV
     */
    virtual const double deltaGamma_Z() const;

    /**
     * @brief The total decay width of the @f$Z@f$ boson, @f$\Gamma_Z@f$.
     * @details
     * \f[
     * \Gamma_Z = \Gamma_Z^{SM} + \Delta \Gamma_Z^{(1)} + \Delta \Gamma_Z^{(2)}
     * \f]
     * @return @f$\Gamma_Z@f$ in GeV, including SM plus @f$\mathcal{O}(\Lambda^{-2})@f$ and @f$\mathcal{O}(\Lambda^{-4})@f$ NP contributions
     *
     * @attention This function is applicable only to the NP model classes that
     * are inherited from NPbase.
     */
    virtual const double Gamma_Z() const; 

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
     * @brief The branching ratio of the @f$W@f$ boson decaying into a SM fermion pair, @f$Br(W\to f_i f_j)@f$.
     * @return @f$Br(W\to f_i f_j)@f$ in GeV
     */
    virtual const double  BrW(const Particle fi, const Particle fj) const;

    /**
     * @brief The lepton universality ratio @f$R_{W,l_i/l_j)=\Gamma(W\to l_i \nu_i)/\Gamma(W\to l_j \nu_j)@f$.
     * @return @f$R_{W,l_i/l_j)@f$ in GeV
     */
    virtual const double  RWlilj(const Particle li, const Particle lj) const;

    /**
     * @brief The ratio @f$R_{W,c)=\Gamma(W\to c + X)/\Gamma(W\to had)@f$.
     * @return @f$R_{W,c)@f$ in GeV
     */
    virtual const double  RWc() const;

    /**
     * @brief The lepton universality ratio @f$R_{Z,l_i/l_j)=\Gamma(Z\to l_i^+ l_i^-)/\Gamma(Z\to l_j^+ l_j^-)@f$.
     * @return @f$R_{Z,l_i/l_j)@f$ in GeV
     */
    virtual const double  RZlilj(const Particle li, const Particle lj) const;
    
    /**
     * @brief The new physics contribution to the cross section for the process @f$e^+ e^-\to Z\to \mathrm{hadrons}@f$
     * at the @f$Z@f$ pole, @f$\delta \sigma_h^0@f$.
     * @return @f$\delta \sigma_h^0@f$ in GeV@f$^{-2}@f$
     */
    virtual const double deltaSigmaHadron() const;
    
    
    /**
     * @brief The cross section for the process @f$e^+ e^-\to Z\to \mathrm{hadrons}@f$
     * at the @f$Z@f$ pole, @f$\sigma_h^0@f$.
     * @details
     * \f[
     * \sigma_h = \sigma_h^{SM} + \Delta \sigma_h^{(1)} + \Delta \sigma_h^{(2)}
     * \f]
     * @return @f$\sigma_h^0@f$ in GeV@f$^{-2}@f$, including SM plus @f$\mathcal{O}(\Lambda^{-2})@f$ and @f$\mathcal{O}(\Lambda^{-4})@f$ NP contributions
     *
     * @attention This function is applicable only to the NP model classes that
     * are inherited from NPbase.
     */
    virtual const double sigma0_had() const;                             //AG:modified

    /**
     * @brief The new physics contribution to the left-right asymmetry in @f$e^+e^-\to Z\to f \bar{f}@f$ at the
     * @f$Z@f$-pole, @f$\delta \mathcal{A}_f@f$. 
     * @param[in] f a lepton or quark
     * @return @f$\delta \mathcal{A}_f@f$
     */
    virtual const double deltaA_f(const Particle f) const;
    
    /**
     * @brief The left-right asymmetry in @f$e^+e^-\to Z\to f \bar{f}@f$ at the
     * @f$Z@f$-pole, @f$\mathcal{A}_f@f$.
     * @details
     * \f[
     * \mathcal{A}_f = \mathcal{A}_f^{SM} + \Delta \mathcal{A}_f^{(1)} + \Delta \mathcal{A}_f^{(2)}
     * \f]
     * @param[in] f a lepton or quark
     * @return @f$\mathcal{A}_f@f$, including SM plus @f$\mathcal{O}(\Lambda^{-2})@f$ and @f$\mathcal{O}(\Lambda^{-4})@f$ NP contributions
     *
     * @attention This function is applicable only to the NP model classes that
     * are inherited from NPbase.
     */
    virtual const double A_f(const Particle f) const;
    
    /**
     * @brief The new physics contribution to the forward-backward asymmetry in @f$e^+e^-\to Z\to f \bar{f}@f$ at the
     * @f$Z@f$-pole, @f$\delta A^f_{FB}@f$. 
     * @param[in] f a lepton or quark
     * @return @f$\delta A^f_{FB}@f$
     */
    virtual const double deltaAFB(const Particle f) const;
    
    /**
     * @brief The forward-backward asymmetry in @f$e^+e^-\to Z\to f \bar{f}@f$ at the
     * @f$Z@f$-pole, @f$A^f_{FB}@f$.
     * @details
     * \f[
     * A_{FB,f} = A_{FB,f}^{SM} + \Delta A_{FB,f}^{(1)} + \Delta A_{FB,f}^{(2)}
     * \f]
     * @param[in] f a lepton or quark
     * @return @f$A^f_{FB}@f$, including SM plus @f$\mathcal{O}(\Lambda^{-2})@f$ and @f$\mathcal{O}(\Lambda^{-4})@f$ NP contributions
     */
    virtual const double AFB(const Particle f) const;
    
    /**
     * @brief The new physics contribution to the ratio @f$R_\ell^0=\Gamma_{\mathrm{had}}/\Gamma_\ell@f$,
     * @f$R_q^0=\Gamma_q/\Gamma_{\mathrm{had}}@f$ and @f$R_\nu^0=\Gamma_\nu/\Gamma_{\mathrm{had}}@f$, 
     * for charged leptons, quarks and neutrinos, respectively.
     * @param f a lepton or quark
     * @return @f$\delta R_f^0@f$
     */
    virtual const double deltaR0_f(const Particle f) const;  
    
    /**
     * @brief The ratio @f$R_\ell^0=\Gamma_{\mathrm{had}}/\Gamma_\ell@f$,
     * @f$R_q^0=\Gamma_q/\Gamma_{\mathrm{had}}@f$ and @f$R_\nu^0=\Gamma_\nu/\Gamma_{\mathrm{had}}@f$, 
     * for charged leptons, quarks and neutrinos, respectively. 
     * @details
     * \f[
     * R^0_f = R_f^{SM} + \Delta R_f^{(1)} + \Delta R_f^{(2)}
     * \f]
     * @param[in] f a lepton or quark
     * @return @f$R_f^0@f$, including SM plus @f$\mathcal{O}(\Lambda^{-2})@f$ and @f$\mathcal{O}(\Lambda^{-4})@f$ NP contributions
     */
    virtual const double R0_f(const Particle f) const; 
    
    ////////////////////////////////////////////////////////////////////////
    
    // LEP2 definitions for 2 to 2 fermion processes

    const double CeeLL_e(const double mu) const;
    const double CeeLL_mu(const double mu) const;
    const double CeeLL_tau(const double mu) const;
    
    const double CeeLL_up(const double mu) const;
    const double CeeLL_charm(const double mu) const;
    const double CeeLL_top(const double mu) const;
    
    const double CeeLL_down(const double mu) const;
    const double CeeLL_strange(const double mu) const;
    const double CeeLL_bottom(const double mu) const;


    const double CeeLR_e(const double mu) const;
    const double CeeLR_mu(const double mu) const;
    const double CeeLR_tau(const double mu) const;
    
    const double CeeLR_up(const double mu) const;
    const double CeeLR_charm(const double mu) const;
    const double CeeLR_top(const double mu) const;
    
    const double CeeLR_down(const double mu) const;
    const double CeeLR_strange(const double mu) const;
    const double CeeLR_bottom(const double mu) const;

    const double CeeRL_e(const double mu) const;
    const double CeeRL_mu(const double mu) const;
    const double CeeRL_tau(const double mu) const;
    
    const double CeeRL_up(const double mu) const;
    const double CeeRL_charm(const double mu) const;
    const double CeeRL_top(const double mu) const;
    
    const double CeeRL_down(const double mu) const;
    const double CeeRL_strange(const double mu) const;
    const double CeeRL_bottom(const double mu) const;

    const double CeeRR_e(const double mu) const;
    const double CeeRR_mu(const double mu) const;
    const double CeeRR_tau(const double mu) const;
    
    const double CeeRR_up(const double mu) const;
    const double CeeRR_charm(const double mu) const;
    const double CeeRR_top(const double mu) const;
    
    const double CeeRR_down(const double mu) const;
    const double CeeRR_strange(const double mu) const;
    const double CeeRR_bottom(const double mu) const;
    
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
    const double deltaMLR2t_e(const double s, const double t) const;      
    const double deltaMRL2t_e(const double s, const double t) const;    

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
    // Top Wilson coefficients in the notation of LHC Top WG arXiv: 1802.07237

    virtual const double cHQpOO(const double mu) const;

    virtual const double cHQmOO(const double mu) const;

    virtual const double cHtOO(const double mu) const;
    
    virtual const double ctWOO(const double mu) const;

    virtual const double ctZOO(const double mu) const;
    
    virtual const double ImctWOO(const double mu) const;

    virtual const double ImctZOO(const double mu) const;
    
    virtual const double ctHOO(const double mu) const;

    virtual const double clQpOO(const double mu) const;

    virtual const double clQmOO(const double mu) const;

    virtual const double cltOO(const double mu) const;

    virtual const double cQeOO(const double mu) const;

    virtual const double cetOO(const double mu) const;
    
    ////////////////////////////////////////////////////////////////////////     
    // EW low-energy observables: Muon g-2
       
    /**
     * @brief The computation of the anomalous magnetic moment of the muon @f$a_\mu=(g_\mu-2)/2@f$
     * @return @f$a_\mu=(g_\mu-2)/2@f$
     */
      virtual const double delta_amuon() const;
      
    // EW low-energy observables: Parity violation

    /**
     * @brief The computation of the electron's weak charge
     * @param[in] q2 the @f$Q^2@f$ at which the weak charge is measured
     * @param[in] y
     * @return @f$Q_{w}(e)@f$
     */
      virtual const double delta_Qwemoller(const double q2, const double y) const;


    /**
     * @brief The computation of the parity violating asymmetry in Moller scattering
     * @param[in] q2 the @f$Q^2@f$ of the process
     * @param[in] y
     * @return @f$A_{LR}@f$
     */
      virtual const double delta_alrmoller(const double q2, const double y) const;


    /**
     * @brief The computation of the proton weak charge: Qwp
     * @return @f$Q_{W}(p)@f$
     */
      virtual const double delta_Qwp() const;

      
    /**
     * @brief The computation of the neutron weak charge: Qwn
     * @return @f$Q_{W}(n)@f$
     */
      virtual const double delta_Qwn() const;
      
    // EW low-energy observables: neutrino scattering
      
    /**
     * @brief The computation of the correction to the effective neutrino nucleon LH coupling: delta_gLnuN2
     * @return @f$\Delta g_L^2(\nu N)@f$
     */
      virtual const double delta_gLnuN2() const;
      
    /**
     * @brief The computation of the correction to the effective neutrino nucleon RH coupling: delta_gRnuN2
     * @return @f$\Delta g_R^2(\nu N)@f$
     */
      virtual const double delta_gRnuN2() const;
      
    /**
     * @brief The computation of the correction to the effective (muon) neutrino-electron vector coupling: delta_gVnue
     * @details
     * @return @f$\Delta g_V^{\nu_\mu e}@f$
     */
      virtual const double delta_gVnue() const;
      
    /**
     * @brief The computation of the correction to the effective (muon) neutrino-electron vector coupling: delta_gAnue
     * @details
     * @return @f$\Delta g_A^{\nu_\mu e}@f$
     */
      virtual const double delta_gAnue() const;
      
    ////////////////////////////////////////////////////////////////////////     
    // Lepton decays
      
    // Lepton Flavor universality tests in Tau decays

    /**
     * @brief The computation of the correction to the LFU ratio @f$g_\mu/ g_e @f$
     * @details
     * @return @f$\delta g_\mu/ g_e @f$
     */
      virtual const double delta_TauLFU_gmuge() const;


    /**
     * @brief The computation of the correction to the LFU ratio @f$g_\tau/ g_\mu @f$
     * @details
     * @return @f$\delta g_\tau/ g_\mu @f$
     */
      virtual const double delta_TauLFU_gtaugmu() const;


    /**
     * @brief The computation of the correction to the LFU ratio @f$g_\tau/ g_e @f$
     * @details
     * @return @f$\delta g_\tau/ g_e @f$
     */
      virtual const double delta_TauLFU_gtauge() const;
      
           
    /**
     * @brief The computation of the correction to the LFU ratio @f$\left(g_\tau/ g_\mu\right)_\pi @f$
     * @details
     * @return @f$\delta \left(g_\tau/ g_\mu\right)_\pi @f$
     */
      virtual const double delta_TauLFU_gtaugmuPi() const;
      
    /**
     * @brief The computation of the correction to the LFU ratio @f$\left(g_\tau/ g_\mu\right)_K @f$
     * @details
     * @return @f$\delta\left(g_\tau/ g_\mu\right)_K @f$
     */
      virtual const double delta_TauLFU_gtaugmuK() const;
      
    ////////////////////////////////////////////////////////////////////////
    // Modifications of EW/Higgs couplings
    
    // Modifications of Zff couplings
      
    /**
     * @brief New physics contribution to the neutral-current vector coupling @f$g_V^f@f$.
     * @param[in] f a lepton or quark
     * @return @f$\delta g_V^f@f$
     */
    virtual const double deltaGV_f(const Particle p) const;

    /**
     * @brief New physics contribution to the neutral-current axial-vector coupling @f$g_A^f@f$.
     * @param[in] f a lepton or quark
     * @return @f$\delta g_A^f@f$
     */
    virtual const double deltaGA_f(const Particle p) const;

    /**
     * @brief New physics contribution to the neutral-current left-handed coupling @f$g_L^f@f$.
     * @param[in] f a lepton or quark
     * @return @f$\delta g_L^f@f$
     */
    // no generation mixing
    const double deltaGL_f(const Particle p) const;

    /**
     * @brief New physics contribution to the neutral-current right-handed coupling @f$g_R^f@f$.
     * @param[in] f a lepton or quark
     * @return @f$\delta g_R^f@f$
     */
    // no generation mixing
    const double deltaGR_f(const Particle p) const;
    
    /**
     * @brief New physics contribution to the neutral-current left-handed coupling @f$g_L^f@f$.
     * @param[in] f a lepton or quark
     * @param[in] mu the RG scale to be used in the evaluation
     * @return @f$\delta g_L^f(\mu)@f$
     */
    // no generation mixing
    const double deltaGL_f_mu(const Particle p, const double mu) const;

    /**
     * @brief New physics contribution to the neutral-current right-handed coupling @f$g_R^f@f$.
     * @param[in] f a lepton or quark
     * @param[in] mu the RG scale to be used in the evaluation
     * @return @f$\delta g_R^f(\mu)@f$
     */
    // no generation mixing
    const double deltaGR_f_mu(const Particle p, const double mu) const;
    
    // Modifications of Wff couplings

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
     * @brief New physics contribution to the charged current coupling @f$W_\mu \bar{f_L}\gamma^mu f_L@f$.
     * @param[in] pbar a lepton or quark
     * @param[in] p a lepton or quark
     * @param[in] mu the RG scale to be used in the evaluation (for those models where it is available)
     * @return @f$\delta g_{Wff}^{L}@f$
     */
    // no generation mixing
    virtual gslpp::complex deltaGL_Wff_mu(const Particle pbar, const Particle p, const double mu) const;
    /**
     * @brief New physics contribution to the charged current coupling @f$W_\mu \bar{f_R}\gamma^mu f_R@f$.
     * @param[in] pbar a lepton or quark
     * @param[in] p a lepton or quark
     * @param[in] mu the RG scale to be used in the evaluation (for those models where it is available)
     * @return @f$\delta g_{Wff}^{R}@f$
     */
    // no generation mixing
    virtual gslpp::complex deltaGR_Wff_mu(const Particle pbar, const Particle p, const double mu) const;
    
    // Modifications of Higgs couplings

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
     * @brief The new physics contribution to the coupling of the effective interaction @f$H G_{\mu\nu} \bar{f}\sigma^{\mu\nu} f@f$.
     * @param[in] p a lepton or quark
     * @return @f$\delta g_{hGff}@f$
     */
    // no generation mixing
    gslpp::complex deltaG_hGff(const Particle p) const;

    /**
     * @brief The new physics contribution to the coupling of the effective interaction @f$H Z_{\mu\nu} \bar{f}\sigma^{\mu\nu} f@f$.
     * @param[in] p a lepton or quark
     * @return @f$\delta g_{hZff}@f$
     */
    // no generation mixing
    gslpp::complex deltaG_hZff(const Particle p) const;

    /**
     * @brief The new physics contribution to the coupling of the effective interaction @f$H A_{\mu\nu} \bar{f}\sigma^{\mu\nu} f@f$.
     * @param[in] p a lepton or quark
     * @return @f$\delta g_{hAff}@f$
     */
    // no generation mixing
    gslpp::complex deltaG_hAff(const Particle p) const;

    /**
     * @brief The new physics contribution to the coupling of the effective interaction @f$G_{\mu\nu}^A \bar{f}\sigma^{\mu\nu} T_A f@f$.
     * @param[in] p a lepton or quark
     * @return @f$\delta g_{Gff}@f$
     */
    // no generation mixing
    gslpp::complex deltaG_Gff(const Particle p) const;

    /**
     * @brief The new physics contribution to the coupling of the effective interaction @f$Z_{\mu\nu} \bar{f}\sigma^{\mu\nu} f@f$.
     * @param[in] p a lepton or quark
     * @return @f$\delta g_{Zff}@f$
     */
    // no generation mixing
    gslpp::complex deltaG_Zff(const Particle p) const;

    /**
     * @brief The new physics contribution to the coupling of the effective interaction @f$A_{\mu\nu} \bar{f}\sigma^{\mu\nu} f@f$.
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
     * @brief The ratio @f$\mu_{ggH}@f$ between the gluon-gluon fusion Higgs
     * production cross-section in the current model and in the Standard Model.
     * @param[in] sqrt_s the center-of-mass energy in TeV
     * @return @f$\mu_{ggH}@f$
     */
    virtual const double muggH(const double sqrt_s) const;
    /**
     * @brief The ratio @f$\mu_{ggHH}@f$ between the gluon-gluon fusion di-Higgs
     * production cross-section in the current model and in the Standard Model.
     * (From arXiv: 1502.00539 [hpe-ph].)
     * @param[in] sqrt_s the center-of-mass energy in TeV
     * @return @f$\mu_{ggHH}@f$
     */
    virtual const double muggHH(const double sqrt_s) const;
    /**
     * @brief The ratio @f$\mu_{VBF}@f$ between the vector-boson fusion Higgs
     * production cross-section in the current model and in the Standard Model. 
     * @param[in] sqrt_s the center-of-mass energy in TeV
     * @return @f$\mu_{VBF}@f$
     */
    virtual const double muVBF(const double sqrt_s) const;
    /**
     * @brief The ratio @f$\mu_{VBF+\gamma}@f$ between the vector-boson fusion Higgs
     * production cross-section in association with a hard photon in the current model
     * and in the Standard Model. 
     * @param[in] sqrt_s the center-of-mass energy in TeV
     * @return @f$\mu_{VBF+\gamma}@f$
     */
    virtual const double muVBFgamma(const double sqrt_s) const;
    /**
     * @brief The ratio @f$\mu_{eeWBF}@f$ between the 
     * @f$ e^{+}e^{-}\to \nu\bar{\nu} H @f$ production
     * cross-section in the current model and in the Standard Model.
     * @param[in] sqrt_s the center-of-mass energy in TeV
     * @return @f$\mu_{eeWBF}@f$
     */
    virtual const double mueeWBF(const double sqrt_s) const;
    /**
     * @brief The ratio @f$\mu_{eeWBF}@f$ between the 
     * @f$ e^{+}e^{-}\to \nu\bar{\nu} H @f$ production
     * cross-section in the current model and in the Standard Model.
     * @param[in] sqrt_s the center-of-mass energy in TeV, Pol_em and Pol_ep
     * are the polarization of electrons and positrons, respectively
     * @return @f$\mu_{eeWBF}@f$
     */
    virtual const double mueeWBFPol(const double sqrt_s, const double Pol_em, const double Pol_ep) const;
    /**
     * @brief The ratio @f$\mu_{e^+e^- \to H\nu\bar{\nu}}@f$ between the 
     * @f$ e^+e^- \to H\nu\bar{\nu} @f$ associated production
     * cross-section in the current model and in the Standard Model.
     * @param[in] sqrt_s the center-of-mass energy in TeV
     * @return @f$\mu_{e^+e^- \to H\nu\bar{\nu}}@f$
     */
    virtual const double mueeHvv(const double sqrt_s) const;
    /**
     * @brief The ratio @f$\mu_{e^+e^- \to H\nu\bar{\nu}}@f$ between the 
     * @f$ e^+e^- \to H\nu\bar{\nu} @f$ associated production
     * cross-section in the current model and in the Standard Model.
     * @param[in] sqrt_s the center-of-mass energy in TeV, Pol_em and Pol_ep
     * are the polarization of electrons and positrons, respectively 
     * @return @f$\mu_{e^+e^- \to H\nu\bar{\nu}}@f$
     */
    virtual const double mueeHvvPol(const double sqrt_s, const double Pol_em, const double Pol_ep) const;
    /**
     * @brief The ratio @f$\mu_{eeZBF}@f$ between the 
     * @f$ e^{+}e^{-}\to e^{+}e^{-} H @f$ production
     * cross-section in the current model and in the Standard Model.
     * @param[in] sqrt_s the center-of-mass energy in TeV
     * @return @f$\mu_{eeZBF}@f$
     */
    virtual const double mueeZBF(const double sqrt_s) const;
    /**
     * @brief The ratio @f$\mu_{eeZBF}@f$ between the 
     * @f$ e^{+}e^{-}\to e^{+}e^{-} H @f$ production
     * cross-section in the current model and in the Standard Model.
     * @param[in] sqrt_s the center-of-mass energy in TeV, Pol_em and Pol_ep
     * are the polarization of electrons and positrons, respectively 
     * @return @f$\mu_{eeZBF}@f$
     */
    virtual const double mueeZBFPol(const double sqrt_s, const double Pol_em, const double Pol_ep) const;
    /**
     * @brief The ratio @f$\mu_{epWBF}@f$ between the 
     * @f$ e^{-} p\to \nu j H @f$ production
     * cross-section in the current model and in the Standard Model.
     * @param[in] sqrt_s the center-of-mass energy in TeV
     * @return @f$\mu_{epWBF}@f$
     */
    virtual const double muepWBF(const double sqrt_s) const;
    /**
     * @brief The ratio @f$\mu_{epZBF}@f$ between the 
     * @f$ e^{-} p\to e^{-} j H @f$ production
     * cross-section in the current model and in the Standard Model.
     * @param[in] sqrt_s the center-of-mass energy in TeV
     * @return @f$\mu_{epZBF}@f$
     */
    virtual const double muepZBF(const double sqrt_s) const;
    /**
     * @brief The ratio @f$\mu_{WH}@f$ between the W-Higgs associated production
     * cross-section in the current model and in the Standard Model. 
     * @param[in] sqrt_s the center-of-mass energy in TeV
     * @return @f$\mu_{WH}@f$
     */
    virtual const double muWH(const double sqrt_s) const;
    /**
     * @brief The ratio @f$\mu_{WH}@f$ between the W-Higgs associated production
     * cross-section in the current model and in the Standard Model, with @f$p_{T,H}>250@f$ GeV.
     * @param[in] sqrt_s the center-of-mass energy in TeV
     * @return @f$\mu_{WH}@f$
     */
    virtual const double muWHpT250(const double sqrt_s) const;
    /**
     * @brief The ratio @f$\mu_{ZH}@f$ between the Z-Higgs associated production
     * cross-section in the current model and in the Standard Model.
     * @param[in] sqrt_s the center-of-mass energy in TeV
     * @return @f$\mu_{ZH}@f$
     */
    virtual const double muZH(const double sqrt_s) const;
    /**
     * @brief The ratio @f$\mu_{ZH}@f$ between the Z-Higgs associated production
     * cross-section in the current model and in the Standard Model, with @f$p_{T,H}>250@f$ GeV.
     * @param[in] sqrt_s the center-of-mass energy in TeV
     * @return @f$\mu_{ZH}@f$
     */
    virtual const double muZHpT250(const double sqrt_s) const;
    
    /**
     * @brief The ratio @f$\mu_{eeZH}@f$ between the 
     * @f$ e^{+}e^{-}\to ZH @f$ associated production
     * cross-section in the current model and in the Standard Model.
     * @param[in] sqrt_s the center-of-mass energy in TeV, Pol_em and Pol_ep
     * are the polarization of electrons and positrons, respectively 
     * @return @f$\mu_{eeZH}@f$
     */
    virtual const double mueeZHGen(const double sqrt_s, const double Pol_em, const double Pol_ep) const;
    /**
     * @brief The ratio @f$\mu_{eeZH}@f$ between the 
     * @f$e^{+}e^{-}\to ZH@f$ associated production
     * cross-section in the current model and in the Standard Model.
     * @param[in] sqrt_s the center-of-mass energy in TeV, Pol_em and Pol_ep
     * are the polarization of electrons and positrons, respectively 
     * @return @f$\mu_{eeZH}@f$
     */
    virtual const double mueeZH(const double sqrt_s, const double Pol_em, const double Pol_ep) const;
    /**
     * @brief The ratio @f$\mu_{eeZH, Z \to e^+ e^-, \mu^+ \mu^-}@f$ between the 
     * @f$ e^{+}e^{-}\to ZH, Z \to e^+ e^-, \mu^+ \mu^- @f$ associated production
     * cross-section in the current model and in the Standard Model.
     * @param[in] sqrt_s the center-of-mass energy in TeV
     * @return @f$\mu_{eeZH, Z \to e^+ e^-, \mu^+ \mu^-}@f$
     */
    virtual const double mueeZllH(const double sqrt_s) const;
    /**
     * @brief The ratio @f$\mu_{eeZH, Z \to q \bar{q}}@f$ between the 
     * @f$ e^{+}e^{-}\to ZH, Z \to q \bar{q} @f$ associated production
     * cross-section in the current model and in the Standard Model.
     * @param[in] sqrt_s the center-of-mass energy in TeV
     * @return @f$\mu_{eeZH, Z \to q \bar{q}}@f$
     */
    virtual const double mueeZqqH(const double sqrt_s) const;
    /**
     * @brief The ratio @f$\mu_{eeZH}@f$ between the 
     * @f$ e^{+}e^{-}\to ZH @f$ associated production
     * cross-section in the current model and in the Standard Model.
     * @param[in] sqrt_s the center-of-mass energy in TeV, Pol_em and Pol_ep
     * are the polarization of electrons and positrons, respectively 
     * @return @f$\mu_{eeZH}@f$
     */
    virtual const double mueeZHPol(const double sqrt_s, const double Pol_em, const double Pol_ep) const;
    /**
     * @brief The ratio @f$\mu_{eeZH, Z \to e^+ e^-, \mu^+ \mu^-}@f$ between the 
     * @f$ e^{+}e^{-}\to ZH, Z \to e^+ e^-, \mu^+ \mu^- @f$ associated production
     * cross-section in the current model and in the Standard Model.
     * @param[in] sqrt_s the center-of-mass energy in TeV, Pol_em and Pol_ep
     * are the polarization of electrons and positrons, respectively 
     * @return @f$\mu_{eeZH, Z \to e^+ e^-, \mu^+ \mu^-}@f$
     */
    virtual const double mueeZllHPol(const double sqrt_s, const double Pol_em, const double Pol_ep) const;
    /**
     * @brief The ratio @f$\mu_{eeZH, Z \to q \bar{q}}@f$ between the 
     * @f$ e^{+}e^{-}\to ZH, Z \to q \bar{q} @f$ associated production
     * cross-section in the current model and in the Standard Model.
     * @param[in] sqrt_s the center-of-mass energy in TeV, Pol_em and Pol_ep
     * are the polarization of electrons and positrons, respectively 
     * @return @f$\mu_{eeZH, Z \to q \bar{q}}@f$
     */
    virtual const double mueeZqqHPol(const double sqrt_s, const double Pol_em, const double Pol_ep) const;
    /**
     * @brief the angular parameter @f$a@f$ from 
     * @f$\mu_{e^+e^- \to ZH}@f$ (arXiv:1708.09079 [hep-ph]).
     * @param[in] sqrt_s the center-of-mass energy in TeV, Pol_em and Pol_ep
     * are the polarization of electrons and positrons, respectively 
     * @return @f$a_{eeZH}@f$
     */
    virtual const double aPskPol(const double sqrt_s, const double Pol_em, const double Pol_ep) const;
    /**
     * @brief the angular parameter @f$b@f$ from 
     * @f$\mu_{e^+e^- \to ZH}@f$ (arXiv:1708.09079 [hep-ph]).
     * @param[in] sqrt_s the center-of-mass energy in TeV, Pol_em and Pol_ep
     * are the polarization of electrons and positrons, respectively 
     * @return @f$b_{eeZH}@f$
     */
    virtual const double bPskPol(const double sqrt_s, const double Pol_em, const double Pol_ep) const;
    /**
     * @brief The ratio @f$\mu_{VH}@f$ between the WH+ZH associated production
     * cross-section in the current model and in the Standard Model.
     * @param[in] sqrt_s the center-of-mass energy in TeV
     * @return @f$\mu_{VH}@f$
     */
    virtual const double muVH(const double sqrt_s) const;
    /**
     * @brief The ratio @f$\mu_{VH}@f$ between the WH+ZH associated production
     * cross-section in the current model and in the Standard Model, with @f$p_{T,H}>250@f$ GeV.
     * @param[in] sqrt_s the center-of-mass energy in TeV
     * @return @f$\mu_{VH}@f$
     */
    virtual const double muVHpT250(const double sqrt_s) const;
    /**
     * @brief The ratio @f$\mu_{VBF+VH}@f$ between the sum of VBF and WH+ZH associated production
     * cross-section in the current model and in the Standard Model.
     * @param[in] sqrt_s the center-of-mass energy in TeV
     * @return @f$\mu_{VBF+VH}@f$
     */
    virtual const double muVBFpVH(const double sqrt_s) const;
    /**
     * @brief The ratio @f$\mu_{ttH}@f$ between the t-tbar-Higgs associated 
     * production cross-section in the current model and in the Standard Model.
     * @param[in] sqrt_s the center-of-mass energy in TeV
     * @return @f$\mu_{ttH}@f$
     */
    virtual const double muttH(const double sqrt_s) const;
    /**
     * @brief The ratio @f$\mu_{tH}@f$ between the t-Higgs associated 
     * production cross-section in the current model and in the Standard Model.
     * @param[in] sqrt_s the center-of-mass energy in TeV
     * @return @f$\mu_{tH}@f$
     */
    virtual const double mutH(const double sqrt_s) const;                             //AG:added
    /**
     * @brief The ratio @f$\mu_{bbH}@f$ between the bbH
     *  production cross-section in the current model and in the Standard Model.
     * @param[in] sqrt_s the center-of-mass energy in TeV
     * @return @f$\mu_{bbH}@f$
     */
    virtual const double mubbH(const double sqrt_s) const
    {
        return 1.0;
    };                             //AG:added
    
    /**
     * @brief The ratio @f$\mu_{tHq}@f$ between the t-q-Higgs associated 
     * production cross-section in the current model and in the Standard Model.
     * @param[in] sqrt_s the center-of-mass energy in TeV
     * @return @f$\mu_{tHq}@f$
     */
    virtual const double mutHq(const double sqrt_s) const;
    /**
     * @brief The ratio @f$\mu_{ggH+ttH}@f$ between the sum of gluon-gluon fusion
     * and t-tbar-Higgs associated 
     * production cross-section in the current model and in the Standard Model.
     * @param[in] sqrt_s the center-of-mass energy in TeV
     * @return @f$\mu_{ggH+ttH}@f$
     */
    virtual const double muggHpttH(const double sqrt_s) const;
    /**
     * @brief The ratio @f$\mu_{eettH}@f$ between the 
     * @f$ e^{+}e^{-}\to t\bar{t} H @f$ production
     * cross-section in the current model and in the Standard Model.
     * @param[in] sqrt_s the center-of-mass energy in TeV
     * @return @f$\mu_{eettH}@f$
     */
    virtual const double mueettH(const double sqrt_s) const;
    /**
     * @brief The ratio @f$\mu_{eettH}@f$ between the 
     * @f$ e^{+}e^{-}\to t\bar{t} H @f$ production
     * cross-section in the current model and in the Standard Model.
     * @param[in] sqrt_s the center-of-mass energy in TeV, Pol_em and Pol_ep
     * are the polarization of electrons and positrons, respectively
     * @return @f$\mu_{eettH}@f$
     */
    virtual const double mueettHPol(const double sqrt_s, const double Pol_em, const double Pol_ep) const;
    /**
     * @brief The ratio @f$\mu_{\mu\mu H}@f$ between the @f$\sigma(\mu \mu \to H)}@f$
     * production cross-section in the current model and in the Standard Model.
     * @param[in] sqrt_s the center-of-mass energy in TeV
     * @return @f$\mu_{\mu\mu H}@f$
     */
    virtual const double mummH(const double sqrt_s) const;
    /**
     * @brief The ratio @f$\mu_{\mu\mu H}@f$ between the @f$\sigma(\mu \mu \to H)}@f$
     * production cross-section in the current model and in the Standard Model, 
     * in the narrow width approximation.
     * @param[in] sqrt_s the center-of-mass energy in TeV
     * @return @f$\mu_{\mu\mu H}@f$
     */
    virtual const double mummHNWA(const double sqrt_s) const;
    /**
     * @brief The ratio @f$\mu_{\mu\mu ZH}@f$ between the @f$\sigma(\mu \mu \to Z H)}@f$
     * production cross-section in the current model and in the Standard Model.
     * @param[in] sqrt_s the center-of-mass energy in TeV
     * @return @f$\mu_{\mu\mu ZH}@f$
     */
    virtual const double mummZH(const double sqrt_s) const;
    /**
     * @brief The ratio @f$\mu_{\mu\mu H\nu\nu}@f$ between the @f$\sigma(\mu \mu \to H \nu \nu)}@f$
     * production cross-section in the current model and in the Standard Model.
     * @param[in] sqrt_s the center-of-mass energy in TeV
     * @return @f$\mu_{\mu\mu H\nu\nu}@f$
     */
    virtual const double mummHvv(const double sqrt_s) const;
    /**
     * @brief The ratio @f$\mu_{\mu\mu H\mu\mu}@f$ between the @f$\sigma(\mu \mu \to H \mu \mu)}@f$
     * production cross-section in the current model and in the Standard Model.
     * @param[in] sqrt_s the center-of-mass energy in TeV
     * @return @f$\mu_{\mu\mu H\mu\mu}@f$
     */
    virtual const double mummHmm(const double sqrt_s) const;
    /**
     * @brief The ratio @f$\mu_{\mu\mu ttH}@f$ between the @f$\sigma(\mu \mu \to t\bar{t} H )}@f$
     * production cross-section in the current model and in the Standard Model.
     * @param[in] sqrt_s the center-of-mass energy in TeV
     * @return @f$\mu_{\mu\mu ttH}@f$
     */
    virtual const double mummttH(const double sqrt_s) const;
    
    //AG:added, begin
    virtual const double delta_muggH_1(const double sqrt_s) const;    ///<The SMEFT linear correction to the ratio @f$\mu_{ggH}@f$ between the gluon-gluon fusion Higgs production cross-section in the current model and in the Standard Model.
    virtual const double delta_muggH_2(const double sqrt_s) const;    ///<The SMEFT quadratic correction to the ratio @f$\mu_{ggH}@f$ between the gluon-gluon fusion Higgs production cross-section in the current model and in the Standard Model.
    
    virtual const double delta_muVBF_1(const double sqrt_s) const;	///<The SMEFT linear correction to the ratio @f$\mu_{VBF}@f$ between the vector-boson fusion Higgs production cross-section in the current model and in the Standard Model.
    virtual const double delta_muVBF_2(const double sqrt_s) const;	///<The SMEFT quadratic correction to the ratio @f$\mu_{VBF}@f$ between the vector-boson fusion Higgs production cross-section in the current model and in the Standard Model.
    
    virtual const double delta_muWH_1(const double sqrt_s) const;	///<The SMEFT linear correction to the ratio @f$\mu_{WH}@f$ between the W-Higgs associated production cross-section in the current model and in the Standard Model.
    virtual const double delta_muWH_2(const double sqrt_s) const;	///<The SMEFT quadratic correction to the ratio @f$\mu_{WH}@f$ between the W-Higgs associated production cross-section in the current model and in the Standard Model.
    
    virtual const double delta_muZH_1(const double sqrt_s) const;	///<The SMEFT linear correction to the ratio @f$\mu_{ZH}@f$ between the Z-Higgs associated production cross-section in the current model and in the Standard Model.
    virtual const double delta_muZH_2(const double sqrt_s) const;	///<The SMEFT quadratic correction to the ratio @f$\mu_{ZH}@f$ between the Z-Higgs associated production cross-section in the current model and in the Standard Model.
    
    virtual const double delta_muVH_1(const double sqrt_s) const;	///<The SMEFT linear correction to the ratio @f$\mu_{VH}@f$ between the Z-Higgs and W-Higgs associated production cross-section in the current model and in the Standard Model.
    virtual const double delta_muVH_2(const double sqrt_s) const;	///<The SMEFT quadratic correction to the ratio @f$\mu_{VH}@f$ between the Z-Higgs and W-Higgs associated production cross-section in the current model and in the Standard Model.
    
    virtual const double delta_muttH_1(const double sqrt_s) const;	///<The SMEFT linear correction to the ratio @f$\mu_{ttH}@f$ between the t-tbar-Higgs associated production cross-section in the current model and in the Standard Model.
    virtual const double delta_muttH_2(const double sqrt_s) const;	///<The SMEFT quadratic correction to the ratio @f$\mu_{ttH}@f$ between the t-tbar-Higgs associated production cross-section in the current model and in the Standard Model.
    
    virtual const double delta_mutH_1(const double sqrt_s) const;	///<The SMEFT linear correction to the ratio @f$\mu_{tH}@f$ between the t-Higgs associated production cross-section in the current model and in the Standard Model.
    virtual const double delta_mutH_2(const double sqrt_s) const;	///<The SMEFT quadratic correction to the ratio @f$\mu_{tH}@f$ between the t-Higgs associated production cross-section in the current model and in the Standard Model.
    
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
    const double   deltaGammaHggRatio1() const;
    /**
     * @brief The new physics contribution to the ratio of the @f$\Gamma(H\to gg)@f$
     * in the current model and in the Standard Model. (Only terms that are
     * quadratic in the effective Lagrangian coefficients.)
     * @return @f$\delta \Gamma(H\to gg)@f$/@f$\Gamma(H\to gg)_{\mathrm{SM}}@f$
     */
    const double   deltaGammaHggRatio2() const;
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
    const double  GammaHWWRatio() const;
    /**
     * @brief The new physics contribution to the ratio of the @f$\Gamma(H\to WW)@f$
     * in the current model and in the Standard Model. (Only terms that are
     * linear in the effective Lagrangian coefficients.)
     * @return @f$\delta \Gamma(H\to WW)@f$/@f$\Gamma(H\to WW)_{\mathrm{SM}}@f$
     */
    const double   deltaGammaHWWRatio1() const;
    /**
     * @brief The new physics contribution to the ratio of the @f$\Gamma(H\to WW)@f$
     * in the current model and in the Standard Model. (Only terms that are
     * quadratic in the effective Lagrangian coefficients.)
     * @return @f$\delta \Gamma(H\to WW)@f$/@f$\Gamma(H\to WW)_{\mathrm{SM}}@f$
     */
    const double   deltaGammaHWWRatio2() const;
    /**
     * @brief The ratio of the Br@f$(H\to WW)@f$ in the current model
     * and in the Standard Model.
     * @return Br@f$(H\to WW)@f$/Br@f$(H\to WW)_{\mathrm{SM}}@f$
     */
    virtual const double BrHWWRatio() const;

    /**
     * @brief The ratio of the @f$\Gamma(H\to W l\nu)@f$ (@f$l=e,\mu @f$) in the current model
     * and in the Standard Model.
     * @return @f$\Gamma(H\to Wl\nu)@f$/@f$\Gamma(H\to Wl\nu)_{\mathrm{SM}}@f$
     */

    const double  GammaHWlvRatio() const;
    /**
     * @brief The new physics contribution to the ratio of the @f$\Gamma(H\to Wl\nu)@f$ (@f$l=e,\mu @f$)
     * in the current model and in the Standard Model. (Only terms that are
     * linear in the effective Lagrangian coefficients.)
     * @return @f$\delta \Gamma(H\to Wl\nu)@f$/@f$\Gamma(H\to Wl\nu)_{\mathrm{SM}}@f$
     */
    const double deltaGammaHWlvRatio1() const;
    /**
     * @brief The new physics contribution to the ratio of the @f$\Gamma(H\to Wl\nu)@f$ (@f$l=e,\mu @f$)
     * in the current model and in the Standard Model. (Only terms that are
     * quadratic in the effective Lagrangian coefficients.)
     * @return @f$\delta \Gamma(H\to Wl\nu)@f$/@f$\Gamma(H\to Wl\nu)_{\mathrm{SM}}@f$
     */
    const double deltaGammaHWlvRatio2() const;
    /**
     * @brief The ratio of the Br@f$(H\to W l\nu)@f$ (@f$l=e,\mu @f$) in the current model
     * and in the Standard Model.
     * @return Br@f$(H\to Wl\nu)@f$/Br@f$(H\to Wl\nu)_{\mathrm{SM}}@f$
     */
    virtual const double BrHWlvRatio() const;

    /**
     * @brief The ratio of the @f$\Gamma(H\to WW^*\to l\nu l\nu)@f$ (@f$l=e,\mu @f$) in the current model
     * and in the Standard Model.
     * @return @f$\Gamma(H\to WW^*\to l\nu l\nu)@f$/@f$\Gamma(H\to WW^*\to l\nu l\nu)_{\mathrm{SM}}@f$
     */
    const double  GammaHWW2l2vRatio() const;
    /**
     * @brief The new physics contribution to the ratio of the @f$\Gamma(H\to WW^*\to l\nu l\nu)@f$ (@f$l=e,\mu @f$)
     * in the current model and in the Standard Model. (Only terms that are
     * linear in the effective Lagrangian coefficients.)
     * @return @f$\delta \Gamma(H\to WW^*\to l\nu l\nu)@f$/@f$\Gamma(H\to WW^*\to l\nu l\nu)_{\mathrm{SM}}@f$
     */
    const double deltaGammaHWW2l2vRatio1() const;
    /**
     * @brief The new physics contribution to the ratio of the @f$\Gamma(H\to WW^*\to l\nu l\nu)@f$ (@f$l=e,\mu @f$)
     * in the current model and in the Standard Model. (Only terms that are
     * quadratic in the effective Lagrangian coefficients.)
     * @return @f$\delta \Gamma(H\to WW^*\to l\nu l\nu)@f$/@f$\Gamma(H\to WW^*\to l\nu l\nu)_{\mathrm{SM}}@f$
     */
    const double deltaGammaHWW2l2vRatio2() const;
    /**
     * @brief The ratio of the Br@f$(H\to WW^*\to l\nu l\nu)@f$ (@f$l=e,\mu @f$) in the current model
     * and in the Standard Model.
     * @return Br@f$(H\to WW^*\to l\nu l\nu)@f$/Br@f$(H\to WW^*\to l\nu l\nu)_{\mathrm{SM}}@f$
     */
    virtual const double BrHWW2l2vRatio() const;

    /**
     * @brief The ratio of the @f$\Gamma(H\to W j j)@f$ in the current model
     * and in the Standard Model.
     * @return @f$\Gamma(H\to W j j)@f$/@f$\Gamma(H\to W j j)_{\mathrm{SM}}@f$
     */
    const double  GammaHWjjRatio() const;
    /**
     * @brief The new physics contribution to the ratio of the @f$\Gamma(H\to W j j)@f$
     * in the current model and in the Standard Model. (Only terms that are
     * linear in the effective Lagrangian coefficients.)
     * @return @f$\delta \Gamma(H\to W j j)@f$/@f$\Gamma(H\to W j j)_{\mathrm{SM}}@f$
     */
    const double deltaGammaHWjjRatio1() const;
    /**
     * @brief The new physics contribution to the ratio of the @f$\Gamma(H\to W j j)@f$
     * in the current model and in the Standard Model. (Only terms that are
     * quadratic in the effective Lagrangian coefficients.)
     * @return @f$\delta \Gamma(H\to W j j)@f$/@f$\Gamma(H\to W j j)_{\mathrm{SM}}@f$
     */
    const double deltaGammaHWjjRatio2() const;
    /**
     * @brief The ratio of the Br@f$(H\to W j j)@f$ in the current model
     * and in the Standard Model.
     * @return Br@f$(H\to W j j)@f$/Br@f$(H\to W j j)_{\mathrm{SM}}@f$
     */
    virtual const double BrHWjjRatio() const;

    /**
     * @brief The ratio of the @f$\Gamma(H\to WW^*\to 4j)@f$ in the current model
     * and in the Standard Model.
     * @return @f$\Gamma(H\to WW^*\to 4j)@f$/@f$\Gamma(H\to WW^*\to 4j)_{\mathrm{SM}}@f$
     */
    const double  GammaHWW4jRatio() const;
    /**
     * @brief The new physics contribution to the ratio of the @f$\Gamma(H\to WW^*\to 4j)@f$
     * in the current model and in the Standard Model. (Only terms that are
     * linear in the effective Lagrangian coefficients.)
     * @return @f$\delta \Gamma(H\to WW^*\to 4j)@f$/@f$\Gamma(H\to WW^*\to 4j)_{\mathrm{SM}}@f$
     */
    const double deltaGammaHWW4jRatio1() const;
    /**
     * @brief The new physics contribution to the ratio of the @f$\Gamma(H\to WW^*\to 4j)@f$
     * in the current model and in the Standard Model. (Only terms that are
     * quadratic in the effective Lagrangian coefficients.)
     * @return @f$\delta \Gamma(H\to WW^*\to 4j)@f$/@f$\Gamma(H\to WW^*\to 4j)_{\mathrm{SM}}@f$
     */
    const double deltaGammaHWW4jRatio2() const;
    /**
     * @brief The ratio of the Br@f$(H\to WW^*\to 4j)@f$ in the current model
     * and in the Standard Model.
     * @return Br@f$(H\to WW^*\to 4j)@f$/Br@f$(H\to WW^*\to 4j)_{\mathrm{SM}}@f$
     */
    virtual const double BrHWW4jRatio() const;

    /**
     * @brief The ratio of the @f$\Gamma(H\to W f f)@f$, with @f$f@f$ any fermion, in the current model
     * and in the Standard Model.
     * @return @f$\Gamma(H\to W f f)@f$/@f$\Gamma(H\to W f f)_{\mathrm{SM}}@f$
     */
    const double  GammaHWffRatio() const;
    /**
     * @brief The new physics contribution to the ratio of the @f$\Gamma(H\to W f f)@f$, with @f$f@f$ any fermion,
     * in the current model and in the Standard Model. (Only terms that are
     * linear in the effective Lagrangian coefficients.)
     * @return @f$\delta \Gamma(H\to W f f)@f$/@f$\Gamma(H\to W f f)_{\mathrm{SM}}@f$
     */
    const double deltaGammaHWffRatio1() const;
    /**
     * @brief The new physics contribution to the ratio of the @f$\Gamma(H\to W f f)@f$, with @f$f@f$ any fermion,
     * in the current model and in the Standard Model. (Only terms that are
     * quadratic in the effective Lagrangian coefficients.)
     * @return @f$\delta \Gamma(H\to W f f)@f$/@f$\Gamma(H\to W f f)_{\mathrm{SM}}@f$
     */
    const double deltaGammaHWffRatio2() const;
    /**
     * @brief The ratio of the Br@f$(H\to W f f)@f$, with @f$f@f$ any fermion, in the current model
     * and in the Standard Model.
     * @return Br@f$(H\to W f f)@f$/Br@f$(H\to W f f)_{\mathrm{SM}}@f$
     */
    virtual const double BrHWffRatio() const;

    /**
     * @brief The ratio of the @f$\Gamma(H\to WW^*\to 4f)@f$, with @f$f@f$ any fermion, in the current model
     * and in the Standard Model.
     * @return @f$\Gamma(H\to WW^*\to 4f)@f$/@f$\Gamma(H\to WW^*\to 4f)_{\mathrm{SM}}@f$
     */
    const double  GammaHWW4fRatio() const;
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
    const double  GammaHZZRatio() const;
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
     * @brief The ratio of the @f$\Gamma(H\to Zll)@f$ (@f$l=e,\mu @f$) in the current model
     * and in the Standard Model.
     * @return @f$\Gamma(H\to Zll)@f$/@f$\Gamma(H\to Zll)_{\mathrm{SM}}@f$
     */
    const double  GammaHZllRatio() const;
    /**
     * @brief The new physics contribution to the ratio of the @f$\Gamma(H\to Zll)@f$ (@f$l=e,\mu @f$)
     * in the current model and in the Standard Model. (Only terms that are
     * linear in the effective Lagrangian coefficients.)
     * @return @f$\delta \Gamma(H\to Zll)@f$/@f$\Gamma(H\to Zll)_{\mathrm{SM}}@f$
     */
    const double deltaGammaHZllRatio1() const;
    /**
     * @brief The new physics contribution to the ratio of the @f$\Gamma(H\to Zll)@f$ (@f$l=e,\mu @f$)
     * in the current model and in the Standard Model. (Only terms that are
     * quadratic in the effective Lagrangian coefficients.)
     * @return @f$\delta \Gamma(H\to Zll)@f$/@f$\Gamma(H\to Zll)_{\mathrm{SM}}@f$
     */
    const double deltaGammaHZllRatio2() const;
    /**
     * @brief The ratio of the Br@f$(H\to Zll)@f$ (@f$l=e,\mu @f$) in the current model
     * and in the Standard Model.
     * @return Br@f$(H\to Zll)@f$/Br@f$(H\to Zll)_{\mathrm{SM}}@f$
     */
    virtual const double BrHZllRatio() const;

    /**
     * @brief The ratio of the @f$\Gamma(H\to Zee)@f$ in the current model
     * and in the Standard Model.
     * @return @f$\Gamma(H\to Zee)@f$/@f$\Gamma(H\to Zee)_{\mathrm{SM}}@f$
     */
    const double  GammaHZeeRatio() const;
    /**
     * @brief The new physics contribution to the ratio of the @f$\Gamma(H\to Zee)@f$
     * in the current model and in the Standard Model. (Only terms that are
     * linear in the effective Lagrangian coefficients.)
     * @return @f$\delta \Gamma(H\to Zee)@f$/@f$\Gamma(H\to Zee)_{\mathrm{SM}}@f$
     */
    const double deltaGammaHZeeRatio1() const;
    /**
     * @brief The new physics contribution to the ratio of the @f$\Gamma(H\to Zee)@f$
     * in the current model and in the Standard Model. (Only terms that are
     * quadratic in the effective Lagrangian coefficients.)
     * @return @f$\delta \Gamma(H\to Zee)@f$/@f$\Gamma(H\to Zee)_{\mathrm{SM}}@f$
     */
    const double deltaGammaHZeeRatio2() const;

    /**
     * @brief The ratio of the @f$\Gamma(H\to Z\mu\mu)@f$ in the current model
     * and in the Standard Model.
     * @return @f$\Gamma(H\to Z\mu\mu)@f$/@f$\Gamma(H\to Z\mu\mu)_{\mathrm{SM}}@f$
     */
    const double  GammaHZmumuRatio() const;
    /**
     * @brief The new physics contribution to the ratio of the @f$\Gamma(H\to Z\mu\mu)@f$
     * in the current model and in the Standard Model. (Only terms that are
     * linear in the effective Lagrangian coefficients.)
     * @return @f$\delta \Gamma(H\to Z\mu\mu)@f$/@f$\Gamma(H\to Z\mu\mu)_{\mathrm{SM}}@f$
     */
    const double deltaGammaHZmumuRatio1() const;
    /**
     * @brief The new physics contribution to the ratio of the @f$\Gamma(H\to Z\mu\mu)@f$
     * in the current model and in the Standard Model. (Only terms that are
     * quadratic in the effective Lagrangian coefficients.)
     * @return @f$\delta \Gamma(H\to Z\mu\mu)@f$/@f$\Gamma(H\to Z\mu\mu)_{\mathrm{SM}}@f$
     */
    const double deltaGammaHZmumuRatio2() const;

    /**
     * @brief The ratio of the @f$\Gamma(H\to ZZ* \to 4l)@f$ (@f$l=e,\mu @f$) in the current model
     * and in the Standard Model.
     * @return @f$\Gamma(H\to ZZ* \to 4l)@f$/@f$\Gamma(H\to ZZ* \to 4l)_{\mathrm{SM}}@f$
     */
    const double  GammaHZZ4lRatio() const;
    /**
     * @brief The new physics contribution to the ratio of the @f$\Gamma(H\to ZZ* \to 4l)@f$ (@f$l=e,\mu @f$)
     * in the current model and in the Standard Model. (Only terms that are
     * linear in the effective Lagrangian coefficients.)
     * @return @f$\delta \Gamma(H\to ZZ* \to 4l)@f$/@f$\Gamma(H\to ZZ* \to 4l)_{\mathrm{SM}}@f$
     */
    const double deltaGammaHZZ4lRatio1() const;
    /**
     * @brief The new physics contribution to the ratio of the @f$\Gamma(H\to ZZ* \to 4l)@f$ (@f$l=e,\mu @f$)
     * in the current model and in the Standard Model. (Only terms that are
     * quadratic in the effective Lagrangian coefficients.)
     * @return @f$\delta \Gamma(H\to ZZ* \to 4l)@f$/@f$\Gamma(H\to ZZ* \to 4l)_{\mathrm{SM}}@f$
     */
    const double deltaGammaHZZ4lRatio2() const;
    /**
     * @brief The ratio of the Br@f$(H\to ZZ* \to 4l)@f$ (@f$l=e,\mu @f$) in the current model
     * and in the Standard Model.
     * @return Br@f$(H\to ZZ* \to 4l)@f$/Br@f$(H\to ZZ* \to 4l)_{\mathrm{SM}}@f$
     */
    virtual const double BrHZZ4lRatio() const;

    /**
     * @brief The ratio of the @f$\Gamma(H\to ZZ* \to 4e)@f$ in the current model
     * and in the Standard Model.
     * @return @f$\Gamma(H\to ZZ* \to 4e)@f$/@f$\Gamma(H\to ZZ* \to 4e)_{\mathrm{SM}}@f$
     */
    const double  GammaHZZ4eRatio() const;
    /**
     * @brief The new physics contribution to the ratio of the @f$\Gamma(H\to ZZ* \to 4e)@f$
     * in the current model and in the Standard Model. (Only terms that are
     * linear in the effective Lagrangian coefficients.)
     * @return @f$\delta \Gamma(H\to ZZ* \to 4e)@f$/@f$\Gamma(H\to ZZ* \to 4e)_{\mathrm{SM}}@f$
     */
    const double deltaGammaHZZ4eRatio1() const;
    /**
     * @brief The new physics contribution to the ratio of the @f$\Gamma(H\to ZZ* \to 4e)@f$
     * in the current model and in the Standard Model. (Only terms that are
     * quadratic in the effective Lagrangian coefficients.)
     * @return @f$\delta \Gamma(H\to ZZ* \to 4e)@f$/@f$\Gamma(H\to ZZ* \to 4e)_{\mathrm{SM}}@f$
     */
    const double deltaGammaHZZ4eRatio2() const;
    /**
     * @brief The ratio of the Br@f$(H\to ZZ* \to 4e)@f$ in the current model
     * and in the Standard Model.
     * @return Br@f$(H\to ZZ* \to 4e)@f$/Br@f$(H\to ZZ* \to 4e)_{\mathrm{SM}}@f$
     */
    virtual const double BrHZZ4eRatio() const;

    /**
     * @brief The ratio of the @f$\Gamma(H\to ZZ* \to 2e2\mu)@f$ in the current model
     * and in the Standard Model.
     * @return @f$\Gamma(H\to ZZ* \to 2e2\mu)@f$/@f$\Gamma(H\to ZZ* \to 2e2\mu)_{\mathrm{SM}}@f$
     */
    const double  GammaHZZ2e2muRatio() const;
    /**
     * @brief The new physics contribution to the ratio of the @f$\Gamma(H\to ZZ* \to 2e2\mu)@f$
     * in the current model and in the Standard Model. (Only terms that are
     * linear in the effective Lagrangian coefficients.)
     * @return @f$\delta \Gamma(H\to ZZ* \to 2e2\mu)@f$/@f$\Gamma(H\to ZZ* \to 2e2\mu)_{\mathrm{SM}}@f$
     */
    const double deltaGammaHZZ2e2muRatio1() const;
    /**
     * @brief The new physics contribution to the ratio of the @f$\Gamma(H\to ZZ* \to 2e2\mu)@f$
     * in the current model and in the Standard Model. (Only terms that are
     * quadratic in the effective Lagrangian coefficients.)
     * @return @f$\delta \Gamma(H\to ZZ* \to 2e2\mu)@f$/@f$\Gamma(H\to ZZ* \to 2e2\mu)_{\mathrm{SM}}@f$
     */
    const double deltaGammaHZZ2e2muRatio2() const;
    /**
     * @brief The ratio of the Br@f$(H\to ZZ* \to 2e 2\mu)@f$ in the current model
     * and in the Standard Model.
     * @return Br@f$(H\to ZZ* \to 2e 2\mu)@f$/Br@f$(H\to ZZ* \to 2e 2\mu)_{\mathrm{SM}}@f$
     */
    virtual const double BrHZZ2e2muRatio() const;

    /**
     * @brief The ratio of the @f$\Gamma(H\to ZZ* \to 4\mu)@f$ in the current model
     * and in the Standard Model.
     * @return @f$\Gamma(H\to ZZ* \to 4\mu)@f$/@f$\Gamma(H\to ZZ* \to 4\mu)_{\mathrm{SM}}@f$
     */
    const double  GammaHZZ4muRatio() const;
    /**
     * @brief The new physics contribution to the ratio of the @f$\Gamma(H\to ZZ* \to 4\mu)@f$
     * in the current model and in the Standard Model. (Only terms that are
     * linear in the effective Lagrangian coefficients.)
     * @return @f$\delta \Gamma(H\to ZZ* \to 4\mu)@f$/@f$\Gamma(H\to ZZ* \to 4\mu)_{\mathrm{SM}}@f$
     */
    const double deltaGammaHZZ4muRatio1() const;
    /**
     * @brief The new physics contribution to the ratio of the @f$\Gamma(H\to ZZ* \to 4\mu)@f$
     * in the current model and in the Standard Model. (Only terms that are
     * quadratic in the effective Lagrangian coefficients.)
     * @return @f$\delta \Gamma(H\to ZZ* \to 4\mu)@f$/@f$\Gamma(H\to ZZ* \to 4\mu)_{\mathrm{SM}}@f$
     */
    const double deltaGammaHZZ4muRatio2() const;
    /**
     * @brief The ratio of the Br@f$(H\to ZZ* \to 4\mu)@f$ in the current model
     * and in the Standard Model.
     * @return Br@f$(H\to ZZ* \to 4\mu)@f$/Br@f$(H\to ZZ* \to 4\mu)_{\mathrm{SM}}@f$
     */
    virtual const double BrHZZ4muRatio() const;

    /**
     * @brief The ratio of the @f$\Gamma(H\to Z\nu\nu)@f$ in the current model
     * and in the Standard Model.
     * @return @f$\Gamma(H\to Z\nu\nu)@f$/@f$\Gamma(H\to Z\nu\nu)_{\mathrm{SM}}@f$
     */
    const double  GammaHZvvRatio() const;
    /**
     * @brief The new physics contribution to the ratio of the @f$\Gamma(H\to Z\nu\nu)@f$
     * in the current model and in the Standard Model. (Only terms that are
     * linear in the effective Lagrangian coefficients.)
     * @return @f$\delta \Gamma(H\to Z\nu\nu)@f$/@f$\Gamma(H\to Z\nu\nu)_{\mathrm{SM}}@f$
     */
    const double deltaGammaHZvvRatio1() const;
    /**
     * @brief The new physics contribution to the ratio of the @f$\Gamma(H\to Z\nu\nu)@f$
     * in the current model and in the Standard Model. (Only terms that are
     * quadratic in the effective Lagrangian coefficients.)
     * @return @f$\delta \Gamma(H\to Z\nu\nu)@f$/@f$\Gamma(H\to Z\nu\nu)_{\mathrm{SM}}@f$
     */
    const double deltaGammaHZvvRatio2() const;
    /**
     * @brief The ratio of the Br@f$(H\to Z\nu\nu)@f$ in the current model
     * and in the Standard Model.
     * @return Br@f$(H\to Z\nu\nu)@f$/Br@f$(H\to Z\nu\nu)_{\mathrm{SM}}@f$
     */
    virtual const double BrHZvvRatio() const;

    /**
     * @brief The ratio of the @f$\Gamma(H\to ZZ* \to 4\nu)@f$ in the current model
     * and in the Standard Model.
     * @return @f$\Gamma(H\to ZZ* \to 4\nu)@f$/@f$\Gamma(H\to ZZ* \to 4\nu)_{\mathrm{SM}}@f$
     */
    const double  GammaHZZ4vRatio() const;
    /**
     * @brief The new physics contribution to the ratio of the @f$\Gamma(H\to ZZ* \to 4\nu)@f$
     * in the current model and in the Standard Model. (Only terms that are
     * linear in the effective Lagrangian coefficients.)
     * @return @f$\delta \Gamma(H\to ZZ* \to 4\nu)@f$/@f$\Gamma(H\to ZZ* \to 4\nu)_{\mathrm{SM}}@f$
     */
    const double deltaGammaHZZ4vRatio1() const;
    /**
     * @brief The new physics contribution to the ratio of the @f$\Gamma(H\to ZZ* \to 4\nu)@f$
     * in the current model and in the Standard Model. (Only terms that are
     * quadratic in the effective Lagrangian coefficients.)
     * @return @f$\delta \Gamma(H\to ZZ* \to 4\nu)@f$/@f$\Gamma(H\to ZZ* \to 4\nu)_{\mathrm{SM}}@f$
     */
    const double deltaGammaHZZ4vRatio2() const;
    /**
     * @brief The ratio of the Br@f$(H\to ZZ* \to 4\nu)@f$ in the current model
     * and in the Standard Model.
     * @return Br@f$(H\to ZZ* \to 4\nu)@f$/Br@f$(H\to ZZ* \to 4\nu)_{\mathrm{SM}}@f$
     */
    virtual const double BrHZZ4vRatio() const;

    /**
     * @brief The ratio of the @f$\Gamma(H\to Zu u)@f$ (@f$u=u,c @f$) in the current model
     * and in the Standard Model.
     * @return @f$\Gamma(H\to Zu u)@f$/@f$\Gamma(H\to Zu u)_{\mathrm{SM}}@f$
     */
    const double  GammaHZuuRatio() const;
    /**
     * @brief The new physics contribution to the ratio of the @f$\Gamma(H\to Z u u)@f$ (@f$u=u,c @f$)
     * in the current model and in the Standard Model. (Only terms that are
     * linear in the effective Lagrangian coefficients.)
     * @return @f$\delta \Gamma(H\to Z u u)@f$/@f$\Gamma(H\to Z u u)_{\mathrm{SM}}@f$
     */
    const double deltaGammaHZuuRatio1() const;
    /**
     * @brief The new physics contribution to the ratio of the @f$\Gamma(H\to Z u u)@f$ (@f$u=u,c @f$)
     * in the current model and in the Standard Model. (Only terms that are
     * quadratic in the effective Lagrangian coefficients.)
     * @return @f$\delta \Gamma(H\to Z u u)@f$/@f$\Gamma(H\to Z u u)_{\mathrm{SM}}@f$
     */
    const double deltaGammaHZuuRatio2() const;
    /**
     * @brief The ratio of the Br@f$(H\to Z u u)@f$ (@f$u=u,c @f$) in the current model
     * and in the Standard Model.
     * @return Br@f$(H\to Z u u)@f$/Br@f$(H\to Z u u)_{\mathrm{SM}}@f$
     */
    virtual const double BrHZuuRatio() const;

    /**
     * @brief The ratio of the @f$\Gamma(H\to ZZ* \to 4 u)@f$ (@f$u=u,c @f$) in the current model
     * and in the Standard Model.
     * @return @f$\Gamma(H\to ZZ* \to 4u)@f$/@f$\Gamma(H\to ZZ* \to 4u)_{\mathrm{SM}}@f$
     */
    const double  GammaHZZ4uRatio() const {
        return 1.0;
    };

    /**
     * @brief The new physics contribution to the ratio of the @f$\Gamma(H\to ZZ* \to 4 u)@f$ (@f$u=u,c @f$)
     * in the current model and in the Standard Model. (Only terms that are
     * linear in the effective Lagrangian coefficients.)
     * @return @f$\delta \Gamma(H\to ZZ* \to 4 u)@f$/@f$\Gamma(H\to ZZ* \to 4 u)_{\mathrm{SM}}@f$
     */
    const double deltaGammaHZZ4uRatio1() const {
        return 0.0;
    };

    /**
     * @brief The new physics contribution to the ratio of the @f$\Gamma(H\to ZZ* \to 4 u)@f$ (@f$u=u,c @f$)
     * in the current model and in the Standard Model. (Only terms that are
     * quadratic in the effective Lagrangian coefficients.)
     * @return @f$\delta \Gamma(H\to ZZ* \to 4 u)@f$/@f$\Gamma(H\to ZZ* \to 4 u)_{\mathrm{SM}}@f$
     */
    const double deltaGammaHZZ4uRatio2() const {
        return 0.0;
    };
    /**
     * @brief The ratio of the Br@f$(H\to ZZ* \to 4 u)@f$ (@f$u=u,c @f$) in the current model
     * and in the Standard Model.
     * @return Br@f$(H\to ZZ* \to 4 u)@f$/Br@f$(H\to ZZ* \to 4 u)_{\mathrm{SM}}@f$
     */
    virtual const double BrHZZ4uRatio() const;

    /**
     * @brief The ratio of the @f$\Gamma(H\to Zd d)@f$ (@f$d=d,s,b @f$) in the current model
     * and in the Standard Model.
     * @return @f$\Gamma(H\to Zd d)@f$/@f$\Gamma(H\to Zd d)_{\mathrm{SM}}@f$
     */
    const double  GammaHZddRatio() const;
    /**
     * @brief The new physics contribution to the ratio of the @f$\Gamma(H\to Z d d)@f$ (@f$d=d,s,b @f$)
     * in the current model and in the Standard Model. (Only terms that are
     * linear in the effective Lagrangian coefficients.)
     * @return @f$\delta \Gamma(H\to Z d d)@f$/@f$\Gamma(H\to Z d d)_{\mathrm{SM}}@f$
     */
    const double deltaGammaHZddRatio1() const;
    /**
     * @brief The new physics contribution to the ratio of the @f$\Gamma(H\to Z d d)@f$ (@f$d=d,s,b @f$)
     * in the current model and in the Standard Model. (Only terms that are
     * quadratic in the effective Lagrangian coefficients.)
     * @return @f$\delta \Gamma(H\to Z d d)@f$/@f$\Gamma(H\to Z d d)_{\mathrm{SM}}@f$
     */
    const double deltaGammaHZddRatio2() const;
    /**
     * @brief The ratio of the Br@f$(H\to Z d d)@f$ (@f$d=d,s,b @f$) in the current model
     * and in the Standard Model.
     * @return Br@f$(H\to Z d d)@f$/Br@f$(H\to Z d d)_{\mathrm{SM}}@f$
     */
    virtual const double BrHZddRatio() const;

    /**
     * @brief The ratio of the @f$\Gamma(H\to ZZ* \to 4 d)@f$ (@f$d=d,s,b @f$) in the current model
     * and in the Standard Model.
     * @return @f$\Gamma(H\to ZZ* \to 4 d)@f$/@f$\Gamma(H\to ZZ* \to 4 d)_{\mathrm{SM}}@f$
     */
    const double  GammaHZZ4dRatio() const {
        return 1.0;
    };

    /**
     * @brief The new physics contribution to the ratio of the @f$\Gamma(H\to ZZ* \to 4 d)@f$ (@f$d=d,s,b @f$)
     * in the current model and in the Standard Model. (Only terms that are
     * linear in the effective Lagrangian coefficients.)
     * @return @f$\delta \Gamma(H\to ZZ* \to 4 d)@f$/@f$\Gamma(H\to ZZ* \to 4 d)_{\mathrm{SM}}@f$
     */
    const double deltaGammaHZZ4dRatio1() const {
        return 0.0;
    };

    /**
     * @brief The new physics contribution to the ratio of the @f$\Gamma(H\to ZZ* \to 4 d)@f$ (@f$d=d,s,b @f$)
     * in the current model and in the Standard Model. (Only terms that are
     * quadratic in the effective Lagrangian coefficients.)
     * @return @f$\delta \Gamma(H\to ZZ* \to 4 d)@f$/@f$\Gamma(H\to ZZ* \to 4 d)_{\mathrm{SM}}@f$
     */
    const double deltaGammaHZZ4dRatio2() const {
        return 0.0;
    };
    /**
     * @brief The ratio of the Br@f$(H\to ZZ* \to 4 d)@f$ (@f$d=d,s,b @f$) in the current model
     * and in the Standard Model.
     * @return Br@f$(H\to ZZ* \to 4 d)@f$/Br@f$(H\to ZZ* \to 4 d)_{\mathrm{SM}}@f$
     */
    virtual const double BrHZZ4dRatio() const;

    /**
     * @brief The ratio of the @f$\Gamma(H\to Zff)@f$, with @f$f@f$ any fermion, in the current model
     * and in the Standard Model.
     * @return @f$\Gamma(H\to Zff)@f$/@f$\Gamma(H\to Zff)_{\mathrm{SM}}@f$
     */
    const double  GammaHZffRatio() const;
    /**
     * @brief The new physics contribution to the ratio of the @f$\Gamma(H\to Z ff)@f$, with @f$f@f$ any fermion,
     * in the current model and in the Standard Model. (Only terms that are
     * linear in the effective Lagrangian coefficients.)
     * @return @f$\delta \Gamma(H\to Z f f)@f$/@f$\Gamma(H\to Z f f)_{\mathrm{SM}}@f$
     */
    const double deltaGammaHZffRatio1() const;
    /**
     * @brief The new physics contribution to the ratio of the @f$\Gamma(H\to Z ff)@f$, with @f$f@f$ any fermion,
     * in the current model and in the Standard Model. (Only terms that are
     * quadratic in the effective Lagrangian coefficients.)
     * @return @f$\delta \Gamma(H\to Z f f)@f$/@f$\Gamma(H\to Z f f)_{\mathrm{SM}}@f$
     */
    const double deltaGammaHZffRatio2() const;
    /**
     * @brief The ratio of the Br@f$(H\to Zff)@f$, with @f$f@f$ any fermion, in the current model
     * and in the Standard Model.
     * @return Br@f$(H\to Zff)@f$/Br@f$(H\to Zff)_{\mathrm{SM}}@f$
     */
    virtual const double BrHZffRatio() const;

    /**
     * @brief The ratio of the @f$\Gamma(H\to ZZ* \to 4f)@f$, with @f$f@f$ any fermion, in the current model
     * and in the Standard Model.
     * @return @f$\Gamma(H\to ZZ* \to 4f)@f$/@f$\Gamma(H\to ZZ* \to 4f)_{\mathrm{SM}}@f$
     */
    const double  GammaHZZ4fRatio() const;
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
    const double  GammaHZgaRatio() const;
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
    const double  GammaHgagaRatio() const;
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

    /**
     * @brief The ratio of the @f$\Gamma(H\to ff)@f$ in the current model
     * and in the Standard Model.
     * @return @f$\Gamma(H\to ff)@f$/@f$\Gamma(H\to ff)_{\mathrm{SM}}@f$
     */
    const double deltaGammaHffRatio1(const double mf, const double CifH) const;          //AG:added
    /**
     * @brief The \f$\mathcal{O}(\Lambda^{-4})\f$ new physics contribution to the
     * ratio of the @f$\Gamma(H\to ff)@f$ in the current model
     * and in the Standard Model at order Lambd.
     * @return @f$\Gamma(H\to ff)@f$/@f$\Gamma(H\to ff)_{\mathrm{SM}}@f$
     */
    const double deltaGammaHffRatio2(const double mf, const double CifH) const;         //AG:added
    
    /**
     * @brief The ratio of the @f$\Gamma(H\to \mu\mu)@f$ in the current model
     * and in the Standard Model.
     * @return @f$\Gamma(H\to \mu\mu)@f$/@f$\Gamma(H\to \mu\mu)_{\mathrm{SM}}@f$
     */
    const double  GammaHmumuRatio() const;
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
    const double  GammaHtautauRatio() const;
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
    const double  GammaHccRatio() const;
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
     * @brief The ratio of the @f$\Gamma(H\to ss)@f$ in the current model
     * and in the Standard Model.
     * @return @f$\Gamma(H\to ss)@f$/@f$\Gamma(H\to ss)_{\mathrm{SM}}@f$
     */
    const double  GammaHssRatio() const;
    /**
     * @brief The new physics contribution to the ratio of the @f$\Gamma(H\to ss)@f$
     * in the current model and in the Standard Model. (Only terms that are
     * linear in the effective Lagrangian coefficients.)
     * @return @f$\delta \Gamma(H\to ss)@f$/@f$\Gamma(H\to ss)_{\mathrm{SM}}@f$
     */
    const double deltaGammaHssRatio1() const;
    /**
     * @brief The new physics contribution to the ratio of the @f$\Gamma(H\to ss)@f$
     * in the current model and in the Standard Model. (Only terms that are
     * quadratic in the effective Lagrangian coefficients.)
     * @return @f$\delta \Gamma(H\to ss)@f$/@f$\Gamma(H\to ss)_{\mathrm{SM}}@f$
     */
    const double deltaGammaHssRatio2() const;
    /**
     * @brief The ratio of the Br@f$(H\to s\bar{s})@f$ in the current model
     * and in the Standard Model.
     * @return Br@f$(H\to s\bar{s})@f$/Br@f$(H\to s\bar{s})_{\mathrm{SM}}@f$
     */
    virtual const double BrHssRatio() const;
    /**
     * @brief The ratio of the @f$\Gamma(H\to bb)@f$ in the current model
     * and in the Standard Model.
     * @return @f$\Gamma(H\to bb)@f$/@f$\Gamma(H\to bb)_{\mathrm{SM}}@f$
     */
    const double  GammaHbbRatio() const;
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

    // General 4 fermion final states

    /**
     * @brief The ratio of the @f$\Gamma(H\to 2L2L')@f$ (@f$L,L'=e,\mu,\tau@f$) in the current model
     * and in the Standard Model.
     * @return @f$\Gamma(H\to 2L2L')@f$/@f$\Gamma(H\to 2L2L')_{\mathrm{SM}}@f$
     */
    const double  GammaH2L2LRatio() const;
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
    const double  GammaH2e2muRatio() const;
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
    const double  GammaH2v2vRatio() const;
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
    const double  GammaH2L2vRatio() const;
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
    const double  GammaH2L2v2Ratio() const;
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
    const double  GammaH2e2vRatio() const;
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
    const double  GammaH2mu2vRatio() const;
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
    const double  GammaH2u2uRatio() const;
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
    const double  GammaH2d2dRatio() const;
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
    const double  GammaH2u2dRatio() const;
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
    const double  GammaH2L2uRatio() const;
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
    const double  GammaH2L2dRatio() const;
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
    const double  GammaH2v2uRatio() const;
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
    const double  GammaH2v2dRatio() const;
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
    const double  GammaH4LRatio() const;
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
    const double  GammaH4L2Ratio() const;
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
    const double  GammaH4eRatio() const;
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
    const double  GammaH4muRatio() const;
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
    const double  GammaH4vRatio() const;
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
    const double  GammaH4uRatio() const;
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
    const double  GammaH4dRatio() const;
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
    const double  GammaHLvvLRatio() const;
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
    const double  GammaHevmuvRatio() const;
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
    const double  GammaHudduRatio() const;
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
    const double  GammaHLvudRatio() const;
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
    const double  GammaH2udRatio() const;
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
    const double  GammaH2LvRatio() const;
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
    const double  GammaH2Lv2Ratio() const;
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
    const double  GammaH2evRatio() const;
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
    const double  GammaH2muvRatio() const;
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
    const double  GammaH4fRatio() const;
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
    const double  GammaH4lRatio() const;
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
    const double  GammaH2l2vRatio() const;
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
    const double  GammaHlljjRatio() const;
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
    const double  GammaHlvjjRatio() const;
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
    const double  GammaHlv_lvorjjRatio() const;
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
    const double  GammaHll_vvorjjRatio() const;
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
    virtual const double muttHZbbboost(const double sqrt_s) const;

    /**
     * @brief The ratio @f$\mu_{ggH,\gamma\gamma}@f$ between the gluon-gluon fusion Higgs
     * production cross-section with subsequent decay into 2 photons in the
     * current model and in the Standard Model.
     * @param[in] sqrt_s the center-of-mass energy in TeV
     * @return @f$\mu_{ggH,\gamma\gamma}@f$
     */
    virtual const double muggHgaga(const double sqrt_s) const;
    /**
     * @brief The ratio @f$\mu_{VBF,\gamma\gamma}@f$ between the VBF Higgs
     * production cross-section with subsequent decay into 2 photons in the
     * current model and in the Standard Model.
     * @param[in] sqrt_s the center-of-mass energy in TeV
     * @return @f$\mu_{VBF,\gamma\gamma}@f$
     */
    virtual const double muVBFHgaga(const double sqrt_s) const;
    /**
     * @brief The ratio @f$\mu_{ZH,\gamma\gamma}@f$ between the ZH
     * production cross-section with subsequent decay into 2 photons in the
     * current model and in the Standard Model.
     * @param[in] sqrt_s the center-of-mass energy in TeV
     * @return @f$\mu_{ZH,\gamma\gamma}@f$
     */
    virtual const double muZHgaga(const double sqrt_s) const;
    /**
     * @brief The ratio @f$\mu_{WH,\gamma\gamma}@f$ between the WH
     * production cross-section with subsequent decay into 2 photons in the
     * current model and in the Standard Model.
     * @param[in] sqrt_s the center-of-mass energy in TeV
     * @return @f$\mu_{WH,\gamma\gamma}@f$
     */
    virtual const double muWHgaga(const double sqrt_s) const;
    /**
     * @brief The ratio @f$\mu_{VH,\gamma\gamma}@f$ between the VH
     * production cross-section with subsequent decay into 2 photons in the
     * current model and in the Standard Model.
     * @param[in] sqrt_s the center-of-mass energy in TeV
     * @return @f$\mu_{VH,\gamma\gamma}@f$
     */
    virtual const double muVHgaga(const double sqrt_s) const;
    /**
     * @brief The ratio @f$\mu_{ttH,\gamma\gamma}@f$ between the ttH
     * production cross-section with subsequent decay into 2 photons in the
     * current model and in the Standard Model.
     * @param[in] sqrt_s the center-of-mass energy in TeV
     * @return @f$\mu_{ttH,\gamma\gamma}@f$
     */
    virtual const double muttHgaga(const double sqrt_s) const;
    /**
     * @brief The ratio @f$\mu_{ggH,Z\gamma}@f$ between the gluon-gluon fusion Higgs
     * production cross-section with subsequent decay into @f$Z \gamma@f$ in the
     * current model and in the Standard Model.
     * @param[in] sqrt_s the center-of-mass energy in TeV
     * @return @f$\mu_{ggH,Z\gamma}@f$
     */
    virtual const double muggHZga(const double sqrt_s) const;
    /**
     * @brief The ratio @f$\mu_{VBF,Z\gamma}@f$ between the VBF Higgs
     * production cross-section with subsequent decay into @f$Z \gamma@f$ in the
     * current model and in the Standard Model.
     * @param[in] sqrt_s the center-of-mass energy in TeV
     * @return @f$\mu_{VBF,Z\gamma}@f$
     */
    virtual const double muVBFHZga(const double sqrt_s) const;
    /**
     * @brief The ratio @f$\mu_{ZH,Z\gamma}@f$ between the ZH
     * production cross-section with subsequent decay into @f$Z \gamma@f$ in the
     * current model and in the Standard Model.
     * @param[in] sqrt_s the center-of-mass energy in TeV
     * @return @f$\mu_{ZH,Z\gamma}@f$
     */
    virtual const double muZHZga(const double sqrt_s) const;
    /**
     * @brief The ratio @f$\mu_{WH,Z\gamma}@f$ between the WH
     * production cross-section with subsequent decay into @f$Z \gamma@f$ in the
     * current model and in the Standard Model.
     * @param[in] sqrt_s the center-of-mass energy in TeV
     * @return @f$\mu_{WH,Z\gamma}@f$
     */
    virtual const double muWHZga(const double sqrt_s) const;
    /**
     * @brief The ratio @f$\mu_{VH,Z\gamma}@f$ between the VH
     * production cross-section with subsequent decay into @f$Z \gamma@f$ in the
     * current model and in the Standard Model.
     * @param[in] sqrt_s the center-of-mass energy in TeV
     * @return @f$\mu_{VH,Z\gamma}@f$
     */
    virtual const double muVHZga(const double sqrt_s) const;
    /**
     * @brief The ratio @f$\mu_{ttH,Z\gamma}@f$ between the ttH
     * production cross-section with subsequent decay into @f$Z \gamma@f$ in the
     * current model and in the Standard Model.
     * @param[in] sqrt_s the center-of-mass energy in TeV
     * @return @f$\mu_{ttH,Z\gamma}@f$
     */
    virtual const double muttHZga(const double sqrt_s) const;
    /**
     * @brief The ratio @f$\mu_{ggH,ZZ}@f$ between the gluon-gluon fusion Higgs
     * production cross-section with subsequent decay into @f$Z Z^*@f$ in the
     * current model and in the Standard Model.
     * @param[in] sqrt_s the center-of-mass energy in TeV
     * @return @f$\mu_{ggH,ZZ}@f$
     */
    virtual const double muggHZZ(const double sqrt_s) const;
    /**
     * @brief The ratio @f$\mu_{VBF,ZZ}@f$ between the VBF Higgs
     * production cross-section with subsequent decay into @f$Z Z^*@f$ in the
     * current model and in the Standard Model.
     * @param[in] sqrt_s the center-of-mass energy in TeV
     * @return @f$\mu_{VBF,ZZ}@f$
     */
    virtual const double muVBFHZZ(const double sqrt_s) const;
    /**
     * @brief The ratio @f$\mu_{ZH,ZZ}@f$ between the ZH
     * production cross-section with subsequent decay into @f$Z Z^*@f$ in the
     * current model and in the Standard Model.
     * @param[in] sqrt_s the center-of-mass energy in TeV
     * @return @f$\mu_{ZH,ZZ}@f$
     */
    virtual const double muZHZZ(const double sqrt_s) const;
    /**
     * @brief The ratio @f$\mu_{WH,ZZ}@f$ between the WH
     * production cross-section with subsequent decay into @f$Z Z^*@f$ in the
     * current model and in the Standard Model.
     * @param[in] sqrt_s the center-of-mass energy in TeV
     * @return @f$\mu_{WH,ZZ}@f$
     */
    virtual const double muWHZZ(const double sqrt_s) const;
    /**
     * @brief The ratio @f$\mu_{VH,ZZ}@f$ between the VH
     * production cross-section with subsequent decay into @f$Z Z^*@f$ in the
     * current model and in the Standard Model.
     * @param[in] sqrt_s the center-of-mass energy in TeV
     * @return @f$\mu_{VH,ZZ}@f$
     */
    virtual const double muVHZZ(const double sqrt_s) const;
    /**
     * @brief The ratio @f$\mu_{ttH,ZZ}@f$ between the ttH
     * production cross-section with subsequent decay into @f$Z Z^*@f$ in the
     * current model and in the Standard Model.
     * @param[in] sqrt_s the center-of-mass energy in TeV
     * @return @f$\mu_{ttH,ZZ}@f$
     */
    virtual const double muttHZZ(const double sqrt_s) const;

    /**
     * @brief The ratio @f$\mu_{ggH,ZZ\to 4l}@f$ between the gluon-gluon fusion Higgs
     * production cross-section with subsequent decay into @f$Z Z^*\to 4l@f$ in the
     * current model and in the Standard Model.
     * @param[in] sqrt_s the center-of-mass energy in TeV
     * @return @f$\mu_{ggH,ZZ\to 4l}@f$
     */
    virtual const double muggHZZ4l(const double sqrt_s) const;
    /**
     * @brief The ratio @f$\mu_{VBF,ZZ\to 4l}@f$ between the VBF Higgs
     * production cross-section with subsequent decay into @f$Z Z^*\to 4l@f$ in the
     * current model and in the Standard Model.
     * @param[in] sqrt_s the center-of-mass energy in TeV
     * @return @f$\mu_{VBF,ZZ\to 4l}@f$
     */
    virtual const double muVBFHZZ4l(const double sqrt_s) const;
    /**
     * @brief The ratio @f$\mu_{ZH,ZZ\to 4l}@f$ between the ZH
     * production cross-section with subsequent decay into @f$Z Z^*\to 4l@f$ in the
     * current model and in the Standard Model.
     * @param[in] sqrt_s the center-of-mass energy in TeV
     * @return @f$\mu_{ZH,ZZ\to 4l}@f$
     */
    virtual const double muZHZZ4l(const double sqrt_s) const;
    /**
     * @brief The ratio @f$\mu_{WH,ZZ\to 4l}@f$ between the WH
     * production cross-section with subsequent decay into @f$Z Z^*\to 4l@f$ in the
     * current model and in the Standard Model.
     * @param[in] sqrt_s the center-of-mass energy in TeV
     * @return @f$\mu_{WH,ZZ\to 4l}@f$
     */
    virtual const double muWHZZ4l(const double sqrt_s) const;
    /**
     * @brief The ratio @f$\mu_{VH,ZZ\to 4l}@f$ between the VH
     * production cross-section with subsequent decay into @f$Z Z^*\to 4l@f$ in the
     * current model and in the Standard Model.
     * @param[in] sqrt_s the center-of-mass energy in TeV
     * @return @f$\mu_{VH,ZZ\to 4l}@f$
     */
    virtual const double muVHZZ4l(const double sqrt_s) const;
    /**
     * @brief The ratio @f$\mu_{ttH,ZZ\to 4l}@f$ between the ttH
     * production cross-section with subsequent decay into @f$Z Z^*\to 4l@f$ in the
     * current model and in the Standard Model.
     * @param[in] sqrt_s the center-of-mass energy in TeV
     * @return @f$\mu_{ttH,ZZ\to 4l}@f$
     */
    virtual const double muttHZZ4l(const double sqrt_s) const;

    /**
     * @brief The ratio @f$\mu_{ggH,WW}@f$ between the gluon-gluon fusion Higgs
     * production cross-section with subsequent decay into @f$W W^*@f$ in the
     * current model and in the Standard Model.
     * @param[in] sqrt_s the center-of-mass energy in TeV
     * @return @f$\mu_{ggH,WW}@f$
     */
    virtual const double muggHWW(const double sqrt_s) const;
    /**
     * @brief The ratio @f$\mu_{VBF,WW}@f$ between the VBF Higgs
     * production cross-section with subsequent decay into @f$W W^*@f$ in the
     * current model and in the Standard Model.
     * @param[in] sqrt_s the center-of-mass energy in TeV
     * @return @f$\mu_{VBF,WW}@f$
     */
    virtual const double muVBFHWW(const double sqrt_s) const;
    /**
     * @brief The ratio @f$\mu_{ZH,WW}@f$ between the ZH
     * production cross-section with subsequent decay into @f$W W^*@f$ in the
     * current model and in the Standard Model.
     * @param[in] sqrt_s the center-of-mass energy in TeV
     * @return @f$\mu_{ZH,WW}@f$
     */
    virtual const double muZHWW(const double sqrt_s) const;
    /**
     * @brief The ratio @f$\mu_{WH,WW}@f$ between the WH
     * production cross-section with subsequent decay into @f$W W^*@f$ in the
     * current model and in the Standard Model.
     * @param[in] sqrt_s the center-of-mass energy in TeV
     * @return @f$\mu_{WH,WW}@f$
     */
    virtual const double muWHWW(const double sqrt_s) const;
    /**
     * @brief The ratio @f$\mu_{VH,WW}@f$ between the VH
     * production cross-section with subsequent decay into @f$W W^*@f$ in the
     * current model and in the Standard Model.
     * @param[in] sqrt_s the center-of-mass energy in TeV
     * @return @f$\mu_{VH,WW}@f$
     */
    virtual const double muVHWW(const double sqrt_s) const;
    /**
     * @brief The ratio @f$\mu_{ttH,WW}@f$ between the ttH
     * production cross-section with subsequent decay into @f$W W^*@f$ in the
     * current model and in the Standard Model.
     * @param[in] sqrt_s the center-of-mass energy in TeV
     * @return @f$\mu_{ttH,WW}@f$
     */
    virtual const double muttHWW(const double sqrt_s) const;
    /**
     * @brief The ratio @f$\mu_{ggH,WW\to 2l2\nu}@f$ between the gluon-gluon fusion Higgs
     * production cross-section with subsequent decay into @f$W W^*\to 2l2\nu@f$ in the
     * current model and in the Standard Model.
     * @param[in] sqrt_s the center-of-mass energy in TeV
     * @return @f$\mu_{ggH,WW\to 2l2\nu}@f$
     */
    virtual const double muggHWW2l2v(const double sqrt_s) const;
    /**
     * @brief The ratio @f$\mu_{VBF,WW\to 2l2\nu}@f$ between the VBF Higgs
     * production cross-section with subsequent decay into @f$W W^*\to 2l2\nu@f$ in the
     * current model and in the Standard Model.
     * @param[in] sqrt_s the center-of-mass energy in TeV
     * @return @f$\mu_{VBF,WW\to 2l2\nu}@f$
     */
    virtual const double muVBFHWW2l2v(const double sqrt_s) const;
    /**
     * @brief The ratio @f$\mu_{ZH,WW\to 2l2\nu}@f$ between the ZH
     * production cross-section with subsequent decay into @f$W W^*\to 2l2\nu@f$ in the
     * current model and in the Standard Model.
     * @param[in] sqrt_s the center-of-mass energy in TeV
     * @return @f$\mu_{ZH,WW\to 2l2\nu}@f$
     */
    virtual const double muZHWW2l2v(const double sqrt_s) const;
    /**
     * @brief The ratio @f$\mu_{WH,WW\to 2l2\nu}@f$ between the WH
     * production cross-section with subsequent decay into @f$W W^*\to 2l2\nu@f$ in the
     * current model and in the Standard Model.
     * @param[in] sqrt_s the center-of-mass energy in TeV
     * @return @f$\mu_{WH,WW\to 2l2\nu}@f$
     */
    virtual const double muWHWW2l2v(const double sqrt_s) const;
    /**
     * @brief The ratio @f$\mu_{VH,WW\to 2l2\nu}@f$ between the VH
     * production cross-section with subsequent decay into @f$W W^*\to 2l2\nu@f$ in the
     * current model and in the Standard Model.
     * @param[in] sqrt_s the center-of-mass energy in TeV
     * @return @f$\mu_{VH,WW\to 2l2\nu}@f$
     */
    virtual const double muVHWW2l2v(const double sqrt_s) const;
    /**
     * @brief The ratio @f$\mu_{ttH,WW\to 2l2\nu}@f$ between the ttH
     * production cross-section with subsequent decay into @f$W W^*\to 2l2\nu@f$ in the
     * current model and in the Standard Model.
     * @param[in] sqrt_s the center-of-mass energy in TeV
     * @return @f$\mu_{ttH,WW\to 2l2\nu}@f$
     */
    virtual const double muttHWW2l2v(const double sqrt_s) const;
    /**
     * @brief The ratio @f$\mu_{ggH,\mu\mu}@f$ between the gluon-gluon fusion Higgs
     * production cross-section with subsequent decay into @f$\mu\mu@f$ in the
     * current model and in the Standard Model.
     * @param[in] sqrt_s the center-of-mass energy in TeV
     * @return @f$\mu_{ggH,\mu\mu}@f$
     */
    virtual const double muggHmumu(const double sqrt_s) const;
    /**
     * @brief The ratio @f$\mu_{VBF,\mu\mu}@f$ between the VBF Higgs
     * production cross-section with subsequent decay into @f$\mu\mu@f$ in the
     * current model and in the Standard Model.
     * @param[in] sqrt_s the center-of-mass energy in TeV
     * @return @f$\mu_{VBF,\mu\mu}@f$
     */
    virtual const double muVBFHmumu(const double sqrt_s) const;
    /**
     * @brief The ratio @f$\mu_{ZH,\mu\mu}@f$ between the ZH
     * production cross-section with subsequent decay into @f$\mu\mu@f$ in the
     * current model and in the Standard Model.
     * @param[in] sqrt_s the center-of-mass energy in TeV
     * @return @f$\mu_{ZH,\mu\mu}@f$
     */
    virtual const double muZHmumu(const double sqrt_s) const;
    /**
     * @brief The ratio @f$\mu_{WH,\mu\mu}@f$ between the WH
     * production cross-section with subsequent decay into @f$\mu\mu@f$ in the
     * current model and in the Standard Model.
     * @param[in] sqrt_s the center-of-mass energy in TeV
     * @return @f$\mu_{WH,\mu\mu}@f$
     */
    virtual const double muWHmumu(const double sqrt_s) const;
    /**
     * @brief The ratio @f$\mu_{VH,\mu\mu}@f$ between the VH
     * production cross-section with subsequent decay into @f$\mu\mu@f$ in the
     * current model and in the Standard Model.
     * @param[in] sqrt_s the center-of-mass energy in TeV
     * @return @f$\mu_{VH,\mu\mu}@f$
     */
    virtual const double muVHmumu(const double sqrt_s) const;
    /**
     * @brief The ratio @f$\mu_{ttH,\mu\mu}@f$ between the ttH
     * production cross-section with subsequent decay into @f$\mu\mu@f$ in the
     * current model and in the Standard Model.
     * @param[in] sqrt_s the center-of-mass energy in TeV
     * @return @f$\mu_{ttH,\mu\mu}@f$
     */
    virtual const double muttHmumu(const double sqrt_s) const;
    /**
     * @brief The ratio @f$\mu_{ggH,\tau\tau}@f$ between the gluon-gluon fusion Higgs
     * production cross-section with subsequent decay into @f$\tau\tau@f$ in the
     * current model and in the Standard Model.
     * @param[in] sqrt_s the center-of-mass energy in TeV
     * @return @f$\mu_{ggH,\tau\tau}@f$
     */
    virtual const double muggHtautau(const double sqrt_s) const;
    /**
     * @brief The ratio @f$\mu_{VBF,\tau\tau}@f$ between the VBF Higgs
     * production cross-section with subsequent decay into @f$\tau\tau@f$ in the
     * current model and in the Standard Model.
     * @param[in] sqrt_s the center-of-mass energy in TeV
     * @return @f$\mu_{VBF,\tau\tau}@f$
     */
    virtual const double muVBFHtautau(const double sqrt_s) const;
    /**
     * @brief The ratio @f$\mu_{ZH,\tau\tau}@f$ between the ZH
     * production cross-section with subsequent decay into @f$\tau\tau@f$ in the
     * current model and in the Standard Model.
     * @param[in] sqrt_s the center-of-mass energy in TeV
     * @return @f$\mu_{ZH,\tau\tau}@f$
     */
    virtual const double muZHtautau(const double sqrt_s) const;
    /**
     * @brief The ratio @f$\mu_{WH,\tau\tau}@f$ between the WH
     * production cross-section with subsequent decay into @f$\tau\tau@f$ in the
     * current model and in the Standard Model.
     * @param[in] sqrt_s the center-of-mass energy in TeV
     * @return @f$\mu_{WH,\tau\tau}@f$
     */
    virtual const double muWHtautau(const double sqrt_s) const;
    /**
     * @brief The ratio @f$\mu_{VH,\tau\tau}@f$ between the VH
     * production cross-section with subsequent decay into @f$\tau\tau@f$ in the
     * current model and in the Standard Model.
     * @param[in] sqrt_s the center-of-mass energy in TeV
     * @return @f$\mu_{VH,\tau\tau}@f$
     */
    virtual const double muVHtautau(const double sqrt_s) const;
    /**
     * @brief The ratio @f$\mu_{ttH,\tau\tau}@f$ between the ttH
     * production cross-section with subsequent decay into @f$\tau\tau@f$ in the
     * current model and in the Standard Model.
     * @param[in] sqrt_s the center-of-mass energy in TeV
     * @return @f$\mu_{ttH,\tau\tau}@f$
     */
    virtual const double muttHtautau(const double sqrt_s) const;
    /**
     * @brief The ratio @f$\mu_{ggH,bb}@f$ between the gluon-gluon fusion Higgs
     * production cross-section with subsequent decay into @f$bb@f$ in the
     * current model and in the Standard Model.
     * @param[in] sqrt_s the center-of-mass energy in TeV
     * @return @f$\mu_{ggH,bb}@f$
     */
    virtual const double muggHbb(const double sqrt_s) const;
    /**
     * @brief The ratio @f$\mu_{VBF,bb}@f$ between the VBF Higgs
     * production cross-section with subsequent decay into @f$bb@f$ in the
     * current model and in the Standard Model.
     * @param[in] sqrt_s the center-of-mass energy in TeV
     * @return @f$\mu_{VBF,bb}@f$
     */
    virtual const double muVBFHbb(const double sqrt_s) const;
    /**
     * @brief The ratio @f$\mu_{ZH,bb}@f$ between the ZH
     * production cross-section with subsequent decay into @f$bb@f$ in the
     * current model and in the Standard Model.
     * @param[in] sqrt_s the center-of-mass energy in TeV
     * @return @f$\mu_{ZH,bb}@f$
     */
    virtual const double muZHbb(const double sqrt_s) const;
    /**
     * @brief The ratio @f$\mu_{WH,bb}@f$ between the WH
     * production cross-section with subsequent decay into @f$bb@f$ in the
     * current model and in the Standard Model.
     * @param[in] sqrt_s the center-of-mass energy in TeV
     * @return @f$\mu_{WH,bb}@f$
     */
    virtual const double muWHbb(const double sqrt_s) const;
    /**
     * @brief The ratio @f$\mu_{VH,bb}@f$ between the VH
     * production cross-section with subsequent decay into @f$bb@f$ in the
     * current model and in the Standard Model.
     * @param[in] sqrt_s the center-of-mass energy in TeV
     * @return @f$\mu_{VH,bb}@f$
     */
    virtual const double muVHbb(const double sqrt_s) const;
    /**
     * @brief The ratio @f$\mu_{ttH,bb}@f$ between the ttH
     * production cross-section with subsequent decay into @f$bb@f$ in the
     * current model and in the Standard Model.
     * @param[in] sqrt_s the center-of-mass energy in TeV
     * @return @f$\mu_{ttH,bb}@f$
     */
    virtual const double muttHbb(const double sqrt_s) const;

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
    virtual const double muTHUggHgaga(const double sqrt_s) const;
    /**
     * @brief The ratio @f$\mu_{VBF,\gamma\gamma}@f$ between the VBF Higgs
     * production cross-section with subsequent decay into 2 photons in the
     * current model and in the Standard Model.
     * @param[in] sqrt_s the center-of-mass energy in TeV
     * @return @f$\mu_{VBF,\gamma\gamma}@f$
     */
    virtual const double muTHUVBFHgaga(const double sqrt_s) const;
    /**
     * @brief The ratio @f$\mu_{ZH,\gamma\gamma}@f$ between the ZH
     * production cross-section with subsequent decay into 2 photons in the
     * current model and in the Standard Model.
     * @param[in] sqrt_s the center-of-mass energy in TeV
     * @return @f$\mu_{ZH,\gamma\gamma}@f$
     */
    virtual const double muTHUZHgaga(const double sqrt_s) const;
    /**
     * @brief The ratio @f$\mu_{WH,\gamma\gamma}@f$ between the WH
     * production cross-section with subsequent decay into 2 photons in the
     * current model and in the Standard Model.
     * @param[in] sqrt_s the center-of-mass energy in TeV
     * @return @f$\mu_{WH,\gamma\gamma}@f$
     */
    virtual const double muTHUWHgaga(const double sqrt_s) const;
    /**
     * @brief The ratio @f$\mu_{VH,\gamma\gamma}@f$ between the VH
     * production cross-section with subsequent decay into 2 photons in the
     * current model and in the Standard Model.
     * @param[in] sqrt_s the center-of-mass energy in TeV
     * @return @f$\mu_{VH,\gamma\gamma}@f$
     */
    virtual const double muTHUVHgaga(const double sqrt_s) const;
    /**
     * @brief The ratio @f$\mu_{ttH,\gamma\gamma}@f$ between the ttH
     * production cross-section with subsequent decay into 2 photons in the
     * current model and in the Standard Model.
     * @param[in] sqrt_s the center-of-mass energy in TeV
     * @return @f$\mu_{ttH,\gamma\gamma}@f$
     */
    virtual const double muTHUttHgaga(const double sqrt_s) const;
    /**
     * @brief The ratio @f$\mu_{ggH,Z\gamma}@f$ between the gluon-gluon fusion Higgs
     * production cross-section with subsequent decay into @f$Z \gamma@f$ in the
     * current model and in the Standard Model.
     * @param[in] sqrt_s the center-of-mass energy in TeV
     * @return @f$\mu_{ggH,Z\gamma}@f$
     */
    virtual const double muTHUggHZga(const double sqrt_s) const;
    /**
     * @brief The ratio @f$\mu_{VBF,Z\gamma}@f$ between the VBF Higgs
     * production cross-section with subsequent decay into @f$Z \gamma@f$ in the
     * current model and in the Standard Model.
     * @param[in] sqrt_s the center-of-mass energy in TeV
     * @return @f$\mu_{VBF,Z\gamma}@f$
     */
    virtual const double muTHUVBFHZga(const double sqrt_s) const;
    /**
     * @brief The ratio @f$\mu_{ZH,Z\gamma}@f$ between the ZH
     * production cross-section with subsequent decay into @f$Z \gamma@f$ in the
     * current model and in the Standard Model.
     * @param[in] sqrt_s the center-of-mass energy in TeV
     * @return @f$\mu_{ZH,Z\gamma}@f$
     */
    virtual const double muTHUZHZga(const double sqrt_s) const;
    /**
     * @brief The ratio @f$\mu_{WH,Z\gamma}@f$ between the WH
     * production cross-section with subsequent decay into @f$Z \gamma@f$ in the
     * current model and in the Standard Model.
     * @param[in] sqrt_s the center-of-mass energy in TeV
     * @return @f$\mu_{WH,Z\gamma}@f$
     */
    virtual const double muTHUWHZga(const double sqrt_s) const;
    /**
     * @brief The ratio @f$\mu_{VH,Z\gamma}@f$ between the VH
     * production cross-section with subsequent decay into @f$Z \gamma@f$ in the
     * current model and in the Standard Model.
     * @param[in] sqrt_s the center-of-mass energy in TeV
     * @return @f$\mu_{VH,Z\gamma}@f$
     */
    virtual const double muTHUVHZga(const double sqrt_s) const;
    /**
     * @brief The ratio @f$\mu_{ttH,Z\gamma}@f$ between the ttH
     * production cross-section with subsequent decay into @f$Z \gamma@f$ in the
     * current model and in the Standard Model.
     * @param[in] sqrt_s the center-of-mass energy in TeV
     * @return @f$\mu_{ttH,Z\gamma}@f$
     */
    virtual const double muTHUttHZga(const double sqrt_s) const;
    /**
     * @brief The ratio @f$\mu_{ggH,ZZ}@f$ between the gluon-gluon fusion Higgs
     * production cross-section with subsequent decay into @f$Z Z^*@f$ in the
     * current model and in the Standard Model.
     * @param[in] sqrt_s the center-of-mass energy in TeV
     * @return @f$\mu_{ggH,ZZ}@f$
     */
    virtual const double muTHUggHZZ(const double sqrt_s) const;
    /**
     * @brief The ratio @f$\mu_{VBF,ZZ}@f$ between the VBF Higgs
     * production cross-section with subsequent decay into @f$Z Z^*@f$ in the
     * current model and in the Standard Model.
     * @param[in] sqrt_s the center-of-mass energy in TeV
     * @return @f$\mu_{VBF,ZZ}@f$
     */
    virtual const double muTHUVBFHZZ(const double sqrt_s) const;
    /**
     * @brief The ratio @f$\mu_{ZH,ZZ}@f$ between the ZH
     * production cross-section with subsequent decay into @f$Z Z^*@f$ in the
     * current model and in the Standard Model.
     * @param[in] sqrt_s the center-of-mass energy in TeV
     * @return @f$\mu_{ZH,ZZ}@f$
     */
    virtual const double muTHUZHZZ(const double sqrt_s) const;
    /**
     * @brief The ratio @f$\mu_{WH,ZZ}@f$ between the WH
     * production cross-section with subsequent decay into @f$Z Z^*@f$ in the
     * current model and in the Standard Model.
     * @param[in] sqrt_s the center-of-mass energy in TeV
     * @return @f$\mu_{WH,ZZ}@f$
     */
    virtual const double muTHUWHZZ(const double sqrt_s) const;
    /**
     * @brief The ratio @f$\mu_{VH,ZZ}@f$ between the VH
     * production cross-section with subsequent decay into @f$Z Z^*@f$ in the
     * current model and in the Standard Model.
     * @param[in] sqrt_s the center-of-mass energy in TeV
     * @return @f$\mu_{VH,ZZ}@f$
     */
    virtual const double muTHUVHZZ(const double sqrt_s) const;
    /**
     * @brief The ratio @f$\mu_{ttH,ZZ}@f$ between the ttH
     * production cross-section with subsequent decay into @f$Z Z^*@f$ in the
     * current model and in the Standard Model.
     * @param[in] sqrt_s the center-of-mass energy in TeV
     * @return @f$\mu_{ttH,ZZ}@f$
     */
    virtual const double muTHUttHZZ(const double sqrt_s) const;

    /**
     * @brief The ratio @f$\mu_{ggH,ZZ\to 4l}@f$ between the gluon-gluon fusion Higgs
     * production cross-section with subsequent decay into @f$Z Z^*\to 4l@f$ in the
     * current model and in the Standard Model.
     * @param[in] sqrt_s the center-of-mass energy in TeV
     * @return @f$\mu_{ggH,ZZ\to 4l}@f$
     */
    virtual const double muTHUggHZZ4l(const double sqrt_s) const;
    /**
     * @brief The ratio @f$\mu_{VBF,ZZ\to 4l}@f$ between the VBF Higgs
     * production cross-section with subsequent decay into @f$Z Z^*\to 4l@f$ in the
     * current model and in the Standard Model.
     * @param[in] sqrt_s the center-of-mass energy in TeV
     * @return @f$\mu_{VBF,ZZ\to 4l}@f$
     */
    virtual const double muTHUVBFHZZ4l(const double sqrt_s) const;
    /**
     * @brief The ratio @f$\mu_{ZH,ZZ\to 4l}@f$ between the ZH
     * production cross-section with subsequent decay into @f$Z Z^*\to 4l@f$ in the
     * current model and in the Standard Model.
     * @param[in] sqrt_s the center-of-mass energy in TeV
     * @return @f$\mu_{ZH,ZZ\to 4l}@f$
     */
    virtual const double muTHUZHZZ4l(const double sqrt_s) const;
    /**
     * @brief The ratio @f$\mu_{WH,ZZ\to 4l}@f$ between the WH
     * production cross-section with subsequent decay into @f$Z Z^*\to 4l@f$ in the
     * current model and in the Standard Model.
     * @param[in] sqrt_s the center-of-mass energy in TeV
     * @return @f$\mu_{WH,ZZ\to 4l}@f$
     */
    virtual const double muTHUWHZZ4l(const double sqrt_s) const;
    /**
     * @brief The ratio @f$\mu_{VH,ZZ\to 4l}@f$ between the VH
     * production cross-section with subsequent decay into @f$Z Z^*\to 4l@f$ in the
     * current model and in the Standard Model.
     * @param[in] sqrt_s the center-of-mass energy in TeV
     * @return @f$\mu_{VH,ZZ\to 4l}@f$
     */
    virtual const double muTHUVHZZ4l(const double sqrt_s) const;
    /**
     * @brief The ratio @f$\mu_{ttH,ZZ\to 4l}@f$ between the ttH
     * production cross-section with subsequent decay into @f$Z Z^*\to 4l@f$ in the
     * current model and in the Standard Model.
     * @param[in] sqrt_s the center-of-mass energy in TeV
     * @return @f$\mu_{ttH,ZZ\to 4l}@f$
     */
    virtual const double muTHUttHZZ4l(const double sqrt_s) const;

    /**
     * @brief The ratio @f$\mu_{ggH,WW}@f$ between the gluon-gluon fusion Higgs
     * production cross-section with subsequent decay into @f$W W^*@f$ in the
     * current model and in the Standard Model.
     * @param[in] sqrt_s the center-of-mass energy in TeV
     * @return @f$\mu_{ggH,WW}@f$
     */
    virtual const double muTHUggHWW(const double sqrt_s) const;
    /**
     * @brief The ratio @f$\mu_{VBF,WW}@f$ between the VBF Higgs
     * production cross-section with subsequent decay into @f$W W^*@f$ in the
     * current model and in the Standard Model.
     * @param[in] sqrt_s the center-of-mass energy in TeV
     * @return @f$\mu_{VBF,WW}@f$
     */
    virtual const double muTHUVBFHWW(const double sqrt_s) const;
    /**
     * @brief The ratio @f$\mu_{ZH,WW}@f$ between the ZH
     * production cross-section with subsequent decay into @f$W W^*@f$ in the
     * current model and in the Standard Model.
     * @param[in] sqrt_s the center-of-mass energy in TeV
     * @return @f$\mu_{ZH,WW}@f$
     */
    virtual const double muTHUZHWW(const double sqrt_s) const;
    /**
     * @brief The ratio @f$\mu_{WH,WW}@f$ between the WH
     * production cross-section with subsequent decay into @f$W W^*@f$ in the
     * current model and in the Standard Model.
     * @param[in] sqrt_s the center-of-mass energy in TeV
     * @return @f$\mu_{WH,WW}@f$
     */
    virtual const double muTHUWHWW(const double sqrt_s) const;
    /**
     * @brief The ratio @f$\mu_{VH,WW}@f$ between the VH
     * production cross-section with subsequent decay into @f$W W^*@f$ in the
     * current model and in the Standard Model.
     * @param[in] sqrt_s the center-of-mass energy in TeV
     * @return @f$\mu_{VH,WW}@f$
     */
    virtual const double muTHUVHWW(const double sqrt_s) const;
    /**
     * @brief The ratio @f$\mu_{ttH,WW}@f$ between the ttH
     * production cross-section with subsequent decay into @f$W W^*@f$ in the
     * current model and in the Standard Model.
     * @param[in] sqrt_s the center-of-mass energy in TeV
     * @return @f$\mu_{ttH,WW}@f$
     */
    virtual const double muTHUttHWW(const double sqrt_s) const;
    /**
     * @brief The ratio @f$\mu_{ggH,WW\to 2l2\nu}@f$ between the gluon-gluon fusion Higgs
     * production cross-section with subsequent decay into @f$W W^*\to 2l2\nu@f$ in the
     * current model and in the Standard Model.
     * @param[in] sqrt_s the center-of-mass energy in TeV
     * @return @f$\mu_{ggH,WW\to 2l2\nu}@f$
     */
    virtual const double muTHUggHWW2l2v(const double sqrt_s) const;
    /**
     * @brief The ratio @f$\mu_{VBF,WW\to 2l2\nu}@f$ between the VBF Higgs
     * production cross-section with subsequent decay into @f$W W^*\to 2l2\nu@f$ in the
     * current model and in the Standard Model.
     * @param[in] sqrt_s the center-of-mass energy in TeV
     * @return @f$\mu_{VBF,WW\to 2l2\nu}@f$
     */
    virtual const double muTHUVBFHWW2l2v(const double sqrt_s) const;
    /**
     * @brief The ratio @f$\mu_{ZH,WW\to 2l2\nu}@f$ between the ZH
     * production cross-section with subsequent decay into @f$W W^*\to 2l2\nu@f$ in the
     * current model and in the Standard Model.
     * @param[in] sqrt_s the center-of-mass energy in TeV
     * @return @f$\mu_{ZH,WW\to 2l2\nu}@f$
     */
    virtual const double muTHUZHWW2l2v(const double sqrt_s) const;
    /**
     * @brief The ratio @f$\mu_{WH,WW\to 2l2\nu}@f$ between the WH
     * production cross-section with subsequent decay into @f$W W^*\to 2l2\nu@f$ in the
     * current model and in the Standard Model.
     * @param[in] sqrt_s the center-of-mass energy in TeV
     * @return @f$\mu_{WH,WW\to 2l2\nu}@f$
     */
    virtual const double muTHUWHWW2l2v(const double sqrt_s) const;
    /**
     * @brief The ratio @f$\mu_{VH,WW\to 2l2\nu}@f$ between the VH
     * production cross-section with subsequent decay into @f$W W^*\to 2l2\nu@f$ in the
     * current model and in the Standard Model.
     * @param[in] sqrt_s the center-of-mass energy in TeV
     * @return @f$\mu_{VH,WW\to 2l2\nu}@f$
     */
    virtual const double muTHUVHWW2l2v(const double sqrt_s) const;
    /**
     * @brief The ratio @f$\mu_{ttH,WW\to 2l2\nu}@f$ between the ttH
     * production cross-section with subsequent decay into @f$W W^*\to 2l2\nu@f$ in the
     * current model and in the Standard Model.
     * @param[in] sqrt_s the center-of-mass energy in TeV
     * @return @f$\mu_{ttH,WW\to 2l2\nu}@f$
     */
    virtual const double muTHUttHWW2l2v(const double sqrt_s) const;
    /**
     * @brief The ratio @f$\mu_{ggH,\mu\mu}@f$ between the gluon-gluon fusion Higgs
     * production cross-section with subsequent decay into @f$\mu\mu@f$ in the
     * current model and in the Standard Model.
     * @param[in] sqrt_s the center-of-mass energy in TeV
     * @return @f$\mu_{ggH,\mu\mu}@f$
     */
    virtual const double muTHUggHmumu(const double sqrt_s) const;
    /**
     * @brief The ratio @f$\mu_{VBF,\mu\mu}@f$ between the VBF Higgs
     * production cross-section with subsequent decay into @f$\mu\mu@f$ in the
     * current model and in the Standard Model.
     * @param[in] sqrt_s the center-of-mass energy in TeV
     * @return @f$\mu_{VBF,\mu\mu}@f$
     */
    virtual const double muTHUVBFHmumu(const double sqrt_s) const;
    /**
     * @brief The ratio @f$\mu_{ZH,\mu\mu}@f$ between the ZH
     * production cross-section with subsequent decay into @f$\mu\mu@f$ in the
     * current model and in the Standard Model.
     * @param[in] sqrt_s the center-of-mass energy in TeV
     * @return @f$\mu_{ZH,\mu\mu}@f$
     */
    virtual const double muTHUZHmumu(const double sqrt_s) const;
    /**
     * @brief The ratio @f$\mu_{WH,\mu\mu}@f$ between the WH
     * production cross-section with subsequent decay into @f$\mu\mu@f$ in the
     * current model and in the Standard Model.
     * @param[in] sqrt_s the center-of-mass energy in TeV
     * @return @f$\mu_{WH,\mu\mu}@f$
     */
    virtual const double muTHUWHmumu(const double sqrt_s) const;
    /**
     * @brief The ratio @f$\mu_{VH,\mu\mu}@f$ between the VH
     * production cross-section with subsequent decay into @f$\mu\mu@f$ in the
     * current model and in the Standard Model.
     * @param[in] sqrt_s the center-of-mass energy in TeV
     * @return @f$\mu_{VH,\mu\mu}@f$
     */
    virtual const double muTHUVHmumu(const double sqrt_s) const;
    /**
     * @brief The ratio @f$\mu_{ttH,\mu\mu}@f$ between the ttH
     * production cross-section with subsequent decay into @f$\mu\mu@f$ in the
     * current model and in the Standard Model.
     * @param[in] sqrt_s the center-of-mass energy in TeV
     * @return @f$\mu_{ttH,\mu\mu}@f$
     */
    virtual const double muTHUttHmumu(const double sqrt_s) const;
    /**
     * @brief The ratio @f$\mu_{ggH,\tau\tau}@f$ between the gluon-gluon fusion Higgs
     * production cross-section with subsequent decay into @f$\tau\tau@f$ in the
     * current model and in the Standard Model.
     * @param[in] sqrt_s the center-of-mass energy in TeV
     * @return @f$\mu_{ggH,\tau\tau}@f$
     */
    virtual const double muTHUggHtautau(const double sqrt_s) const;
    /**
     * @brief The ratio @f$\mu_{VBF,\tau\tau}@f$ between the VBF Higgs
     * production cross-section with subsequent decay into @f$\tau\tau@f$ in the
     * current model and in the Standard Model.
     * @param[in] sqrt_s the center-of-mass energy in TeV
     * @return @f$\mu_{VBF,\tau\tau}@f$
     */
    virtual const double muTHUVBFHtautau(const double sqrt_s) const;
    /**
     * @brief The ratio @f$\mu_{ZH,\tau\tau}@f$ between the ZH
     * production cross-section with subsequent decay into @f$\tau\tau@f$ in the
     * current model and in the Standard Model.
     * @param[in] sqrt_s the center-of-mass energy in TeV
     * @return @f$\mu_{ZH,\tau\tau}@f$
     */
    virtual const double muTHUZHtautau(const double sqrt_s) const;
    /**
     * @brief The ratio @f$\mu_{WH,\tau\tau}@f$ between the WH
     * production cross-section with subsequent decay into @f$\tau\tau@f$ in the
     * current model and in the Standard Model.
     * @param[in] sqrt_s the center-of-mass energy in TeV
     * @return @f$\mu_{WH,\tau\tau}@f$
     */
    virtual const double muTHUWHtautau(const double sqrt_s) const;
    /**
     * @brief The ratio @f$\mu_{VH,\tau\tau}@f$ between the VH
     * production cross-section with subsequent decay into @f$\tau\tau@f$ in the
     * current model and in the Standard Model.
     * @param[in] sqrt_s the center-of-mass energy in TeV
     * @return @f$\mu_{VH,\tau\tau}@f$
     */
    virtual const double muTHUVHtautau(const double sqrt_s) const;
    /**
     * @brief The ratio @f$\mu_{ttH,\tau\tau}@f$ between the ttH
     * production cross-section with subsequent decay into @f$\tau\tau@f$ in the
     * current model and in the Standard Model.
     * @param[in] sqrt_s the center-of-mass energy in TeV
     * @return @f$\mu_{ttH,\tau\tau}@f$
     */
    virtual const double muTHUttHtautau(const double sqrt_s) const;
    /**
     * @brief The ratio @f$\mu_{ggH,bb}@f$ between the gluon-gluon fusion Higgs
     * production cross-section with subsequent decay into @f$bb@f$ in the
     * current model and in the Standard Model.
     * @param[in] sqrt_s the center-of-mass energy in TeV
     * @return @f$\mu_{ggH,bb}@f$
     */
    virtual const double muTHUggHbb(const double sqrt_s) const;
    /**
     * @brief The ratio @f$\mu_{VBF,bb}@f$ between the VBF Higgs
     * production cross-section with subsequent decay into @f$bb@f$ in the
     * current model and in the Standard Model.
     * @param[in] sqrt_s the center-of-mass energy in TeV
     * @return @f$\mu_{VBF,bb}@f$
     */
    virtual const double muTHUVBFHbb(const double sqrt_s) const;
    /**
     * @brief The ratio @f$\mu_{ZH,bb}@f$ between the ZH
     * production cross-section with subsequent decay into @f$bb@f$ in the
     * current model and in the Standard Model.
     * @param[in] sqrt_s the center-of-mass energy in TeV
     * @return @f$\mu_{ZH,bb}@f$
     */
    virtual const double muTHUZHbb(const double sqrt_s) const;
    /**
     * @brief The ratio @f$\mu_{WH,bb}@f$ between the WH
     * production cross-section with subsequent decay into @f$bb@f$ in the
     * current model and in the Standard Model.
     * @param[in] sqrt_s the center-of-mass energy in TeV
     * @return @f$\mu_{WH,bb}@f$
     */
    virtual const double muTHUWHbb(const double sqrt_s) const;
    /**
     * @brief The ratio @f$\mu_{VH,bb}@f$ between the VH
     * production cross-section with subsequent decay into @f$bb@f$ in the
     * current model and in the Standard Model.
     * @param[in] sqrt_s the center-of-mass energy in TeV
     * @return @f$\mu_{VH,bb}@f$
     */
    virtual const double muTHUVHbb(const double sqrt_s) const;
    /**
     * @brief The ratio @f$\mu_{ttH,bb}@f$ between the ttH
     * production cross-section with subsequent decay into @f$bb@f$ in the
     * current model and in the Standard Model.
     * @param[in] sqrt_s the center-of-mass energy in TeV
     * @return @f$\mu_{ttH,bb}@f$
     */
    virtual const double muTHUttHbb(const double sqrt_s) const;

    /**
     * @brief The ratio @f$\mu_{VBF}@f$ between the VBF
     * production cross-section in the
     * current model and in the Standard Model, multiplied by the 
     * total (SM+new physics) invisible decay branching ratio.
     * @param[in] sqrt_s the center-of-mass energy in TeV
     * @return @f$\mu_{VBF}BR_{inv}@f$
     */
    virtual const double muTHUVBFBRinv(const double sqrt_s) const;
    /**
     * @brief The ratio @f$\mu_{VBF,inv}@f$ between the VBF
     * production cross-section with subsequent decay into invisible states in the
     * current model and in the Standard Model.
     * @param[in] sqrt_s the center-of-mass energy in TeV
     * @return @f$\mu_{VBF,inv}@f$
     */
    virtual const double muTHUVBFHinv(const double sqrt_s) const;
    /**
     * @brief The ratio @f$\mu_{VH}@f$ between the VH
     * production cross-section in the
     * current model and in the Standard Model, multiplied by the 
     * total (SM+new physics) invisible decay branching ratio.
     * @param[in] sqrt_s the center-of-mass energy in TeV
     * @return @f$\mu_{VH}BR_{inv}@f$
     */
    virtual const double muTHUVHBRinv(const double sqrt_s) const;
    /**
     * @brief The ratio @f$\mu_{VH,inv}@f$ between the VH
     * production cross-section with subsequent decay into invisible states in the
     * current model and in the Standard Model.
     * @param[in] sqrt_s the center-of-mass energy in TeV
     * @return @f$\mu_{VH,inv}@f$
     */
    virtual const double muTHUVHinv(const double sqrt_s) const;

    /**
     * @brief The ratio @f$\mu_{ggH,ZZ\to 4\mu}@f$ between the gluon-gluon fusion Higgs
     * production cross-section with subsequent decay into @f$Z Z^*\to 4\mu@f$ in the
     * current model and in the Standard Model.
     * @param[in] sqrt_s the center-of-mass energy in TeV
     * @return @f$\mu_{ggH,ZZ\to 4\mu}@f$
     */
    virtual const double muTHUggHZZ4mu(const double sqrt_s) const;
    /**
     * @brief The ratio @f$\mu_{ggH,Z\gamma\to \gamma 2\mu}@f$ between the gluon-gluon fusion Higgs
     * production cross-section with subsequent decay into @f$Z \gamma\to \gamma 2\mu@f$ in the
     * current model and in the Standard Model.
     * @param[in] sqrt_s the center-of-mass energy in TeV
     * @return @f$\mu_{ggH,Z\gamma\to \gamma 2\mu}@f$
     */
    virtual const double muTHUggHZgamumu(const double sqrt_s) const;

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
    virtual const double deltaxseeWW4fLEP2(const double sqrt_s, const int fstate) const;

    /**
     * @brief The cross section in pb for @f$e^+ e^- \to W^+ W^- \to 4f @f$, 
     * with @f$ 4f = 0 (jjjj), 1 (e v jj), 2 (mu v jj), 3 (tau v jj),
     * 4 (e v e v), 5 (mu v mu v), 6 (tau v tau v),
     * 7 (e v mu v), 8 (e v tau v), 9 (mu v tau v), 10 (l v jj), 11 (l v l v) @f$
     * the different fermion final states for C.O.M. energies in 188-208 GeV.
     * From arXiv: 1606.06693 [hep-ph].
     * @return @f$sigma@f$ [pb]
     */
    virtual const double xseeWW4fLEP2(const double sqrt_s, const int fstate) const;

    /**
     * @brief The new physics contribution to the total cross section in pb for @f$e^+ e^- \to W^+ W^-@f$, 
     * summing over all final states for C.O.M. energies in 188-208 GeV.
     * From arXiv: 1606.06693 [hep-ph].
     * @return @f$\delta sigma@f$ [pb]
     */
    virtual const double deltaxseeWWtotLEP2(const double sqrt_s) const;


    /**
     * @brief The total cross section in pb for @f$e^+ e^- \to W^+ W^-@f$, 
     * summing over all final states for C.O.M. energies in 188-208 GeV.
     * From arXiv: 1606.06693 [hep-ph].
     * @return @f$sigma@f$ [pb]
     */
    virtual const double xseeWWtotLEP2(const double sqrt_s) const;

    /**
     * @brief The new physics contribution to the differential cross section in pb for @f$e^+ e^- \to W^+ W^- \to lv jj @f$, 
     * with @f$ l= e,\mu @f$ for the 4 @f$ cos{\theta}@f$ bins defined in arXiv: 1606.06693 [hep-ph].
     * for the C.O.M. energies of 182.6 and 205.9 GeV.
     * From arXiv: 1606.06693 [hep-ph].
     * @return @f$\delta d\sigma/d\cos{\theta}@f$ [pb]
     */
    virtual const double deltadxsdcoseeWWlvjjLEP2(const double sqrt_s, const int bin) const;

    /**
     * @brief The differential cross section in pb for @f$e^+ e^- \to W^+ W^- \to lv jj @f$, 
     * with @f$ l= e,\mu @f$ for the 4 @f$ cos{\theta}@f$ bins defined in arXiv: 1606.06693 [hep-ph].
     * for the C.O.M. energies of 182.6 and 205.9 GeV.
     * From arXiv: 1606.06693 [hep-ph].
     * @return @f$d\sigma/d\cos{\theta}@f$ [pb]
     */
    virtual const double dxsdcoseeWWlvjjLEP2(const double sqrt_s, const int bin) const;

    /**
     * @brief The differential distribution for @f$e^+ e^- \to W^+ W^- \to jj \ell \nu@f$, 
     * with @f$\ell= e, \mu@f$, as a function of the @f$W@f$ polar angle.
     * @return @f$d\sigma/d\cos{\theta}@f$
     */
    virtual const double dxseeWWdcos(const double sqrt_s, const double cos) const;

    /**
     * @brief The integral of differential distribution for @f$e^+ e^- \to W^+ W^- \to jj \ell \nu@f$, 
     * with @f$\ell= e, \mu@f$ in a given bin of the @f$W@f$ polar angle.
     * @return @f$\int_{\cos{\theta_1}}^{\cos{\theta_2}} d\sigma/d\cos{\theta}@f$
     */
    virtual const double dxseeWWdcosBin(const double sqrt_s, const double cos1, const double cos2) const;

    /**
     * @brief Total @f$e^+ e^- \to W^+ W^- \to jj \ell \nu@f$ cross section in pb, 
     * with @f$\ell= e, \mu@f$.
     * @return @f$\sigma(e^+ e^- \to W^+ W^- \to jj \ell \nu) @f$
     */
    virtual const double xseeWW(const double sqrt_s) const;

    /**
     * @brief The ratio @f$\mu_{eeWW}@f$ between the 
     * @f$ e^{+}e^{-}\to W^{+}W^{-} @f$ production
     * cross-section in the current model and in the Standard Model.
     * @param[in] sqrt_s the center-of-mass energy in TeV
     * @return @f$\mu_{eeWW}@f$
     */
    virtual const double mueeWW(const double sqrt_s) const;

    /**
     * @brief The ratio @f$\mu_{eeWW}@f$ between the 
     * @f$ e^{+}e^{-}\to W^{+}W^{-} @f$ production
     * cross-section in the current model and in the Standard Model.
     * @param[in] sqrt_s the center-of-mass energy in TeV, Pol_em and Pol_ep
     * are the polarization of electrons and positrons, respectively 
     * @return @f$\mu_{eeWW}@f$
     */
    virtual const double mueeWWPol(const double sqrt_s, const double Pol_em, const double Pol_ep) const;

    ////////////////////////////////////////////////////////////////////////

    //----- High Energy diboson observables at hadron colliders

    /**
     * @brief The direction constrained by @f$ p p \to Z H@f$ in the boosted regime, @f$g_p^Z@f$.
     * From arXiv:1807.01796 and the contribution to FCC CDR Vol 1. Implemented only in NPSMEFTd6General class.
     * @return @f$g_p^Z@f$
     */
    virtual const double ppZHprobe(const double sqrt_s) const;

    /**
     * @brief The number of events in  @f$ p p \to WZ@f$
     * in a given @f$p_{TV}@f$ bin, normalized to the SM prediction.
     * From arXiv: 1712.01310 [hep-ph] and private communication.
     * Implemented only in NPSMEFTd6General class.
     * @return @f$N_{ev}^{p_{TV}}/N_{ev,SM}^{p_{TV}}@f$
     */
    virtual const double mupTVppWZ(const double sqrt_s, const double pTV1, const double pTV2) const;


    ////////////////////////////////////////////////////////////////////////

    //----- Simplified Template Cross Sections Bins

    //----- Stage 0

    /**
     * @brief The STXS0 bin @f$pp \to H qq@f$.
     * @param[in] sqrt_s the center-of-mass energy in TeV
     */
    virtual const double STXS0_qqH(const double sqrt_s) const;


    //----- Stage 1

    /**
     * @brief The STXS bin @f$gg \to H@f$.
     * @param[in] sqrt_s the center-of-mass energy in TeV
     */
    virtual const double STXS_ggH_VBFtopo_j3v(const double sqrt_s) const;


    /**
     * @brief The STXS bin @f$gg \to H@f$.
     * @param[in] sqrt_s the center-of-mass energy in TeV
     */
    virtual const double STXS_ggH_VBFtopo_j3(const double sqrt_s) const;

    /**
     * @brief The STXS bin @f$gg \to H@f$.
     * @param[in] sqrt_s the center-of-mass energy in TeV
     */
    virtual const double STXS_ggH0j(const double sqrt_s) const;

    /**
     * @brief The STXS bin @f$gg \to H@f$.
     * @param[in] sqrt_s the center-of-mass energy in TeV
     */
    virtual const double STXS_ggH1j_pTH_0_60(const double sqrt_s) const;

    /**
     * @brief The STXS bin @f$gg \to H@f$.
     * @param[in] sqrt_s the center-of-mass energy in TeV
     */
    virtual const double STXS_ggH1j_pTH_60_120(const double sqrt_s) const;

    /**
     * @brief The STXS bin @f$gg \to H@f$.
     * @param[in] sqrt_s the center-of-mass energy in TeV
     */
    virtual const double STXS_ggH1j_pTH_120_200(const double sqrt_s) const;

    /**
     * @brief The STXS bin @f$gg \to H@f$.
     * @param[in] sqrt_s the center-of-mass energy in TeV
     */
    virtual const double STXS_ggH1j_pTH_200(const double sqrt_s) const;

    /**
     * @brief The STXS bin @f$gg \to H@f$.
     * @param[in] sqrt_s the center-of-mass energy in TeV
     */
    virtual const double STXS_ggH2j_pTH_0_200(const double sqrt_s) const;

    /**
     * @brief The STXS bin @f$gg \to H@f$.
     * @param[in] sqrt_s the center-of-mass energy in TeV
     */
    virtual const double STXS_ggH2j_pTH_0_60(const double sqrt_s) const;


    /**
     * @brief The STXS bin @f$gg \to H@f$.
     * @param[in] sqrt_s the center-of-mass energy in TeV
     */
    virtual const double STXS_ggH2j_pTH_60_120(const double sqrt_s) const;


    /**
     * @brief The STXS bin @f$gg \to H@f$.
     * @param[in] sqrt_s the center-of-mass energy in TeV
     */
    virtual const double STXS_ggH2j_pTH_120_200(const double sqrt_s) const;


    /**
     * @brief The STXS bin @f$gg \to H@f$.
     * @param[in] sqrt_s the center-of-mass energy in TeV
     */
    virtual const double STXS_ggH2j_pTH_200(const double sqrt_s) const;

    /**
     * @brief The STXS bin @f$qq \to H qq@f$.
     * @param[in] sqrt_s the center-of-mass energy in TeV
     */
    virtual const double STXS_qqHqq_VBFtopo_Rest(const double sqrt_s) const;


    /**
     * @brief The STXS bin @f$qq \to H qq@f$.
     * @param[in] sqrt_s the center-of-mass energy in TeV
     */
    virtual const double STXS_qqHqq_VBFtopo_j3v(const double sqrt_s) const;


    /**
     * @brief The STXS bin @f$qq \to H qq@f$.
     * @param[in] sqrt_s the center-of-mass energy in TeV
     */
    virtual const double STXS_qqHqq_VBFtopo_j3(const double sqrt_s) const;

    /**
     * @brief The STXS bin @f$qq \to H qq@f$.
     * @param[in] sqrt_s the center-of-mass energy in TeV
     */
    virtual const double STXS_qqHqq_nonVHtopo(const double sqrt_s) const;

    /**
     * @brief The STXS bin @f$qq \to H qq@f$.
     * @param[in] sqrt_s the center-of-mass energy in TeV
     */
    virtual const double STXS_qqHqq_VHtopo(const double sqrt_s) const;

    /**
     * @brief The STXS bin @f$qq \to H qq@f$.
     * @param[in] sqrt_s the center-of-mass energy in TeV
     */
    virtual const double STXS_qqHqq_Rest(const double sqrt_s) const;

    /**
     * @brief The STXS bin @f$qq \to H qq@f$.
     * @param[in] sqrt_s the center-of-mass energy in TeV
     */
    virtual const double STXS_qqHqq_pTj_200(const double sqrt_s) const;

    /**
     * @brief The STXS bin @f$qq \to H \ell \nu@f$.
     * @param[in] sqrt_s the center-of-mass energy in TeV
     */
    virtual const double STXS_qqHlv_pTV_0_250(const double sqrt_s) const;


    /**
     * @brief The STXS bin @f$qq \to H \ell \nu@f$.
     * @param[in] sqrt_s the center-of-mass energy in TeV
     */
    virtual const double STXS_qqHlv_pTV_0_150(const double sqrt_s) const;


    /**
     * @brief The STXS bin @f$qq \to H \ell \nu@f$.
     * @param[in] sqrt_s the center-of-mass energy in TeV
     */
    virtual const double STXS_qqHlv_pTV_150_250_0j(const double sqrt_s) const;


    /**
     * @brief The STXS bin @f$qq \to H \ell \nu@f$.
     * @param[in] sqrt_s the center-of-mass energy in TeV
     */
    virtual const double STXS_qqHlv_pTV_150_250_1j(const double sqrt_s) const;


    /**
     * @brief The STXS bin @f$qq \to H \ell \nu@f$.
     * @param[in] sqrt_s the center-of-mass energy in TeV
     */
    virtual const double STXS_qqHlv_pTV_250(const double sqrt_s) const;

    /**
     * @brief The STXS bin @f$qq \to H \ell \ell@f$.
     * @param[in] sqrt_s the center-of-mass energy in TeV
     */
    virtual const double STXS_qqHll_pTV_0_150(const double sqrt_s) const;

    /**
     * @brief The STXS bin @f$qq \to H \ell \ell@f$.
     * @param[in] sqrt_s the center-of-mass energy in TeV
     */
    virtual const double STXS_qqHll_pTV_150_250(const double sqrt_s) const;


    /**
     * @brief The STXS bin @f$qq \to H \ell \ell@f$.
     * @param[in] sqrt_s the center-of-mass energy in TeV
     */
    virtual const double STXS_qqHll_pTV_150_250_0j(const double sqrt_s) const;


    /**
     * @brief The STXS bin @f$qq \to H \ell \ell@f$.
     * @param[in] sqrt_s the center-of-mass energy in TeV
     */
    virtual const double STXS_qqHll_pTV_150_250_1j(const double sqrt_s) const;

    /**
     * @brief The STXS bin @f$qq \to H \ell \ell@f$.
     * @param[in] sqrt_s the center-of-mass energy in TeV
     */
    virtual const double STXS_qqHll_pTV_250(const double sqrt_s) const;


    /**
     * @brief The STXS bin @f$ ttH + tH @f$.
     * @param[in] sqrt_s the center-of-mass energy in TeV
     */
    virtual const double STXS_ttHtH(const double sqrt_s) const;

    /**
     * @brief The STXS bin @f$ qq \to WH \to H qq @f$.
     * @param[in] sqrt_s the center-of-mass energy in TeV
     */
    virtual const double STXS_WHqqHqq_VBFtopo_j3v(const double sqrt_s) const;

    /**
     * @brief The STXS bin @f$ qq \to WH \to H qq @f$.
     * @param[in] sqrt_s the center-of-mass energy in TeV
     */
    virtual const double STXS_WHqqHqq_VBFtopo_j3(const double sqrt_s) const;

    /**
     * @brief The STXS bin @f$ qq \to WH \to H qq @f$.
     * @param[in] sqrt_s the center-of-mass energy in TeV
     */
    virtual const double STXS_WHqqHqq_VH2j(const double sqrt_s) const;

    /**
     * @brief The STXS bin @f$ qq \to WH \to H qq @f$.
     * @param[in] sqrt_s the center-of-mass energy in TeV
     */
    virtual const double STXS_WHqqHqq_Rest(const double sqrt_s) const;

    /**
     * @brief The STXS bin @f$ qq \to WH \to H qq @f$.
     * @param[in] sqrt_s the center-of-mass energy in TeV
     */
    virtual const double STXS_WHqqHqq_pTj1_200(const double sqrt_s) const;

    /**
     * @brief The STXS bin @f$ qq \to ZH \to H qq @f$.
     * @param[in] sqrt_s the center-of-mass energy in TeV
     */
    virtual const double STXS_ZHqqHqq_VBFtopo_j3v(const double sqrt_s) const;

    /**
     * @brief The STXS bin @f$ qq \to ZH \to H qq @f$.
     * @param[in] sqrt_s the center-of-mass energy in TeV
     */
    virtual const double STXS_ZHqqHqq_VBFtopo_j3(const double sqrt_s) const;

    /**
     * @brief The STXS bin @f$ qq \to ZH \to H qq @f$.
     * @param[in] sqrt_s the center-of-mass energy in TeV
     */
    virtual const double STXS_ZHqqHqq_VH2j(const double sqrt_s) const;

    /**
     * @brief The STXS bin @f$ qq \to ZH \to H qq @f$.
     * @param[in] sqrt_s the center-of-mass energy in TeV
     */
    virtual const double STXS_ZHqqHqq_Rest(const double sqrt_s) const;

    /**
     * @brief The STXS bin @f$ qq \to ZH \to H qq @f$.
     * @param[in] sqrt_s the center-of-mass energy in TeV
     */
    virtual const double STXS_ZHqqHqq_pTj1_200(const double sqrt_s) const;


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
    virtual const double STXS12_ggH_pTH200_300_Nj01(const double sqrt_s) const;

    /**
     * @brief The STXS bin @f$gg \to H@f$, @f$N_j\leq 1,~300<p_{TH} [GeV]<450@f$.
     * @param[in] sqrt_s the center-of-mass energy in TeV
     */
    virtual const double STXS12_ggH_pTH300_450_Nj01(const double sqrt_s) const;

    /**
     * @brief The STXS bin @f$gg \to H@f$, @f$N_j\leq 1,~450<p_{TH} [GeV]<650@f$.
     * @param[in] sqrt_s the center-of-mass energy in TeV
     */
    virtual const double STXS12_ggH_pTH450_650_Nj01(const double sqrt_s) const;

    /**
     * @brief The STXS bin @f$gg \to H@f$, @f$N_j\leq 1,650<p_{TH} [GeV]@f$.
     * @param[in] sqrt_s the center-of-mass energy in TeV
     */
    virtual const double STXS12_ggH_pTH650_Inf_Nj01(const double sqrt_s) const;

    /**
     * @brief The STXS bin @f$gg \to H@f$, @f$N_j=0,~p_{TH} [GeV]<10@f$.
     * @param[in] sqrt_s the center-of-mass energy in TeV
     */
    virtual const double STXS12_ggH_pTH0_10_Nj0(const double sqrt_s) const;
    
    
    /**
     * @brief The STXS bin @f$gg \to H@f$, @f$N_j=0,~p_{TH} [GeV]<10@f$.
     * @param[in] sqrt_s the center-of-mass energy in TeV
     */
    virtual const double STXS12_ggH_pTH10_200_Nj0(const double sqrt_s) const;

    /**
     * @brief The STXS bin @f$gg \to H@f$, @f$N_j=0,~10<p_{TH} [GeV]@f$.
     * @param[in] sqrt_s the center-of-mass energy in TeV
     */
    virtual const double STXS12_ggH_pTH10_Inf_Nj0(const double sqrt_s) const;
    
    /**
     * @brief The STXS bin @f$gg \to H@f$, @f$N_j = 1,~p_{TH} [GeV]<60@f$.
     * @param[in] sqrt_s the center-of-mass energy in TeV
     */
    virtual const double STXS12_ggH_pTH0_60_Nj1(const double sqrt_s) const;
    
    
    /**
     * @brief The STXS bin @f$gg \to H@f$, @f$N_j \geq 1,~m_{jj}[GeV]<350,~p_{TH} [GeV]<60@f$.
     * @param[in] sqrt_s the center-of-mass energy in TeV
     */
    virtual const double STXS12_ggH_mjj0_350_pTH0_60_Nj1(const double sqrt_s) const;

    /**
     * @brief The STXS bin @f$gg \to H@f$, @f$N_j = 1,~60<p_{TH} [GeV]<120@f$.
     * @param[in] sqrt_s the center-of-mass energy in TeV
     */
    virtual const double STXS12_ggH_pTH60_120_Nj1(const double sqrt_s) const;

    /**
     * @brief The STXS bin @f$gg \to H@f$, @f$N_j = 1,~120<p_{TH} [GeV]<200@f$.
     * @param[in] sqrt_s the center-of-mass energy in TeV
     */
    virtual const double STXS12_ggH_pTH120_200_Nj1(const double sqrt_s) const;

    /**
     * @brief The STXS bin @f$gg \to H@f$, @f$N_j \geq 2,~m_{jj}[GeV]<350,~p_{TH} [GeV]<60@f$.
     * @param[in] sqrt_s the center-of-mass energy in TeV
     */
    virtual const double STXS12_ggH_mjj0_350_pTH0_60_Nj2(const double sqrt_s) const;

    /**
     * @brief The STXS bin @f$gg \to H@f$, @f$N_j \geq 2,~m_{jj}[GeV]<350,~60<p_{TH} [GeV]<120@f$.
     * @param[in] sqrt_s the center-of-mass energy in TeV
     */
    virtual const double STXS12_ggH_mjj0_350_pTH60_120_Nj2(const double sqrt_s) const;

    /**
     * @brief The STXS bin @f$gg \to H@f$, @f$N_j \geq 2,~m_{jj}[GeV]<350,~120<p_{TH} [GeV]<200@f$.
     * @param[in] sqrt_s the center-of-mass energy in TeV
     */
    virtual const double STXS12_ggH_mjj0_350_pTH120_200_Nj2(const double sqrt_s) const;

    /**
     * @brief The STXS bin @f$gg \to H@f$, @f$N_j \geq 2,~350<m_{jj}[GeV]<700,~p_{TH} [GeV]<200,~p_{THjj}[GeV]<25@f$.
     * @param[in] sqrt_s the center-of-mass energy in TeV
     */
    virtual const double STXS12_ggH_mjj350_700_pTH0_200_ptHjj0_25_Nj2(const double sqrt_s) const;

    /**
     * @brief The STXS bin @f$gg \to H@f$, @f$N_j \geq 2,~350<m_{jj}[GeV]<700,~p_{TH} [GeV]<200,~25<p_{THjj}[GeV]@f$.
     * @param[in] sqrt_s the center-of-mass energy in TeV
     */
    virtual const double STXS12_ggH_mjj350_700_pTH0_200_ptHjj25_Inf_Nj2(const double sqrt_s) const;

    /**
     * @brief The STXS bin @f$gg \to H@f$, @f$N_j \geq 2,~700<m_{jj}[GeV],~p_{TH} [GeV]<200,~p_{THjj}[GeV]<25@f$.
     * @param[in] sqrt_s the center-of-mass energy in TeV
     */
    virtual const double STXS12_ggH_mjj700_Inf_pTH0_200_ptHjj0_25_Nj2(const double sqrt_s) const;

    /**
     * @brief The STXS bin @f$gg \to H@f$, @f$N_j \geq 2,~700<m_{jj}[GeV],~p_{TH} [GeV]<200,~25<p_{THjj}[GeV]@f$.
     * @param[in] sqrt_s the center-of-mass energy in TeV
     */
    virtual const double STXS12_ggH_mjj700_Inf_pTH0_200_ptHjj25_Inf_Nj2(const double sqrt_s) const;

    /**
     * @brief The STXS bin @f$gg \to H@f$, @f$N_j \geq 2,~350<m_{jj}[GeV]<700,~<p_{TH} [GeV]<200@f$.
     * @param[in] sqrt_s the center-of-mass energy in TeV
     */
    virtual const double STXS12_ggH_mjj350_700_pTH0_200_Nj2(const double sqrt_s) const;           //AG:added
    
    /**
     * @brief The STXS bin @f$gg \to H@f$, @f$N_j \geq 2,~700<m_{jj}[GeV],~p_{TH} [GeV]<200@f$.
     * @param[in] sqrt_s the center-of-mass energy in TeV
     */
    virtual const double STXS12_ggH_mjj700_Inf_pTH0_200_Nj2(const double sqrt_s) const;           //AG:added
    
    
    /**
     * @brief The STXS bin @f$gg \to H@f$, @f$,200<~p_{TH} [GeV]<300@f$.
     * @param[in] sqrt_s the center-of-mass energy in TeV
     */
    virtual const double STXS12_ggH_pTH200_300(const double sqrt_s) const;
    
    
    /**
     * @brief The STXS bin @f$gg \to H@f$, @f$,300<~p_{TH} [GeV]<450@f$.
     * @param[in] sqrt_s the center-of-mass energy in TeV
     */
    virtual const double STXS12_ggH_pTH300_450(const double sqrt_s) const;
    
    
    /**
     * @brief The STXS bin @f$gg \to H@f$, @f$,450<~p_{TH} [GeV]@f$.
     * @param[in] sqrt_s the center-of-mass energy in TeV
     */
    virtual const double STXS12_ggH_pTH450_Inf(const double sqrt_s) const;
    
    
    /**
     * @brief The STXS bin @f$gg \to H@f$, @f$450<~p_{TH} [GeV]<650@f$.
     * @param[in] sqrt_s the center-of-mass energy in TeV
     */
    virtual const double STXS12_ggH_pTH450_650(const double sqrt_s) const;
    
    
    /**
     * @brief The STXS bin @f$gg \to H@f$, @f$~p_{TH} [GeV]>650@f$.
     * @param[in] sqrt_s the center-of-mass energy in TeV
     */
    virtual const double STXS12_ggH_pTH650_Inf(const double sqrt_s) const;
    
    
    /**
     * @brief The STXS bin @f$gg \to H\ell\ell@f$, @f$p_{TV}[GeV]<75@f$.
     * @param[in] sqrt_s the center-of-mass energy in TeV
     */
    virtual const double STXS12_ggHll_pTV0_75(const double sqrt_s) const;

    /**
     * @brief The STXS bin @f$gg \to H\ell\ell@f$, @f$75<p_{TV}[GeV]<150@f$.
     * @param[in] sqrt_s the center-of-mass energy in TeV
     */
    virtual const double STXS12_ggHll_pTV75_150(const double sqrt_s) const;

    /**
     * @brief The STXS bin @f$gg \to H\ell\ell@f$, @f$N_j = 0,~150<p_{TV}[GeV]<250@f$.
     * @param[in] sqrt_s the center-of-mass energy in TeV
     */
    virtual const double STXS12_ggHll_pTV150_250_Nj0(const double sqrt_s) const;

    /**
     * @brief The STXS bin @f$gg \to H\ell\ell@f$, @f$N_j = 1,~150<p_{TV}[GeV]<250@f$.
     * @param[in] sqrt_s the center-of-mass energy in TeV
     */
    virtual const double STXS12_ggHll_pTV150_250_Nj1(const double sqrt_s) const;

    /**
     * @brief The STXS bin @f$gg \to H\ell\ell@f$, @f$250 < p_{TV}[GeV]@f$.
     * @param[in] sqrt_s the center-of-mass energy in TeV
     */
    virtual const double STXS12_ggHll_pTV250_Inf(const double sqrt_s) const;


    
    /**
     * @brief The STXS bin @f$qq \to Hqq@f$, @f$N_j = 0,1@f$ VH-veto Ref. 2402.05742.
     * @param[in] sqrt_s the center-of-mass energy in TeV
     */
    virtual const double STXS12_qqHqq_VH_veto_Nj01(const double sqrt_s) const;
    
    
    /**
     * @brief The STXS bin @f$qq \to Hqq@f$, @f$N_j \ge 2@f$ VH-had Ref. 2402.05742.
     * @param[in] sqrt_s the center-of-mass energy in TeV
     */
    virtual const double STXS12_qqHqq_VH_had_Nj2(const double sqrt_s) const;
    
    
    /**
     * @brief The STXS bin @f$qq \to Hqq@f$, @f$N_j = 0@f$.
     * @param[in] sqrt_s the center-of-mass energy in TeV
     */
    virtual const double STXS12_qqHqq_Nj0(const double sqrt_s) const;

    /**
     * @brief The STXS bin @f$qq \to Hqq@f$, @f$N_j = 1@f$.
     * @param[in] sqrt_s the center-of-mass energy in TeV
     */
    virtual const double STXS12_qqHqq_Nj1(const double sqrt_s) const;

    /**
     * @brief The STXS bin @f$qq \to Hqq@f$, @f$N_j \geq 2,~m_{jj}[GeV]<60@f$.
     * @param[in] sqrt_s the center-of-mass energy in TeV
     */
    virtual const double STXS12_qqHqq_mjj0_60_Nj2(const double sqrt_s) const;

    /**
     * @brief The STXS bin @f$qq \to Hqq@f$, @f$N_j \geq 2,~60<m_{jj}[GeV]<120@f$.
     * @param[in] sqrt_s the center-of-mass energy in TeV
     */
    virtual const double STXS12_qqHqq_mjj60_120_Nj2(const double sqrt_s) const;

    /**
     * @brief The STXS bin @f$qq \to Hqq@f$, @f$N_j \geq 2,~120<m_{jj}[GeV]<350@f$.
     * @param[in] sqrt_s the center-of-mass energy in TeV
     */
    virtual const double STXS12_qqHqq_mjj120_350_Nj2(const double sqrt_s) const;

    /**
     * @brief The STXS bin @f$qq \to Hqq@f$, @f$N_j \geq 2,~350<m_{jj}[GeV],~200<p_{TH}[GeV]@f$.
     * @param[in] sqrt_s the center-of-mass energy in TeV
     */
    virtual const double STXS12_qqHqq_mjj350_Inf_pTH200_Inf_Nj2(const double sqrt_s) const;

    /**
     * @brief The STXS bin @f$qq \to Hqq@f$, @f$N_j \geq 2,~350<m_{jj}[GeV]<700,~p_{TH}[GeV]<200,~p_{THjj}[GeV]<25@f$.
     * @param[in] sqrt_s the center-of-mass energy in TeV
     */
    virtual const double STXS12_qqHqq_mjj350_700_pTH0_200_pTHjj0_25_Nj2(const double sqrt_s) const;

    /**
     * @brief The STXS bin @f$qq \to Hqq@f$, @f$N_j \geq 2,~350<m_{jj}[GeV]<700,~p_{TH}[GeV]<200,~25<p_{THjj}[GeV]@f$.
     * @param[in] sqrt_s the center-of-mass energy in TeV
     */
    virtual const double STXS12_qqHqq_mjj350_700_pTH0_200_pTHjj25_Inf_Nj2(const double sqrt_s) const;

    /**
     * @brief The STXS bin @f$qq \to Hqq@f$, @f$N_j \geq 2,~700<m_{jj}[GeV],~p_{TH}[GeV]<200,~p_{THjj}[GeV]<25@f$.
     * @param[in] sqrt_s the center-of-mass energy in TeV
     */
    virtual const double STXS12_qqHqq_mjj700_Inf_pTH0_200_pTHjj0_25_Nj2(const double sqrt_s) const;

    /**
     * @brief The STXS bin @f$qq \to Hqq@f$, @f$N_j \geq 2,~700<m_{jj}[GeV],~p_{TH}[GeV]<200,~25<p_{THjj}[GeV]@f$.
     * @param[in] sqrt_s the center-of-mass energy in TeV
     */
    virtual const double STXS12_qqHqq_mjj700_Inf_pTH0_200_pTHjj25_Inf_Nj2(const double sqrt_s) const;

    /**
     * @brief The STXS bin @f$qq \to Hqq@f$, @f$N_j \geq 2,~350<m_{jj}[GeV]<700,~p_{TH}[GeV]<200@f$.
     * @param[in] sqrt_s the center-of-mass energy in TeV
     */
    virtual const double STXS12_qqHqq_mjj350_700_pTH0_200_Nj2(const double sqrt_s) const;     //AG:added
    
    /**
     * @brief The STXS bin @f$qq \to Hqq@f$, @f$N_j \geq 2,~700<m_{jj}[GeV]<1000,~p_{TH}[GeV]<200@f$.
     * @param[in] sqrt_s the center-of-mass energy in TeV
     */
    virtual const double STXS12_qqHqq_mjj700_1000_pTH0_200_Nj2(const double sqrt_s) const;    //AG:added
    
    /**
     * @brief The STXS bin @f$qq \to Hqq@f$, @f$N_j \geq 2,~1000<m_{jj}[GeV]<1500,~p_{TH}[GeV]<200@f$.
     * @param[in] sqrt_s the center-of-mass energy in TeV
     */
    virtual const double STXS12_qqHqq_mjj1000_1500_pTH0_200_Nj2(const double sqrt_s) const;   //AG:added
    
    /**
     * @brief The STXS bin @f$qq \to Hqq@f$, @f$N_j \geq 2,~1500<m_{jj}[GeV],~p_{TH}[GeV]<200@f$.
     * @param[in] sqrt_s the center-of-mass energy in TeV
     */
    virtual const double STXS12_qqHqq_mjj1500_Inf_pTH0_200_Nj2(const double sqrt_s) const;    //AG:added
    
    /**
     * @brief The STXS bin @f$qq \to Hqq@f$, @f$N_j \geq 2,~350<m_{jj}[GeV]<1000,~p_{TH}[GeV]>200@f$.
     * @param[in] sqrt_s the center-of-mass energy in TeV
     */
    virtual const double STXS12_qqHqq_mjj350_1000_pTH200_Inf_Nj2(const double sqrt_s) const;    //AG:added
    
    /**
     * @brief The STXS bin @f$qq \to Hqq@f$, @f$N_j \geq 2,~1000<m_{jj}[GeV],~p_{TH}[GeV]>200@f$.
     * @param[in] sqrt_s the center-of-mass energy in TeV
     */
    virtual const double STXS12_qqHqq_mjj1000_Inf_pTH200_Inf_Nj2(const double sqrt_s) const;    //AG:added
    
    /**
     * @brief The STXS bin @f$qq \to H\ell\nu@f$, @f$p_{TV}[GeV]<75@f$.
     * @param[in] sqrt_s the center-of-mass energy in TeV
     */
    virtual const double STXS12_qqHlv_pTV0_75(const double sqrt_s) const;

    /**
     * @brief The STXS bin @f$qq \to H\ell\nu@f$, @f$75<p_{TV}[GeV]<150@f$.
     * @param[in] sqrt_s the center-of-mass energy in TeV
     */
    virtual const double STXS12_qqHlv_pTV75_150(const double sqrt_s) const;

    /**
     * @brief The STXS bin @f$qq \to H\ell\nu@f$, @f$N_j = 0,~150<p_{TV}[GeV]<250@f$.
     * @param[in] sqrt_s the center-of-mass energy in TeV
     */
    virtual const double STXS12_qqHlv_pTV150_250_Nj0(const double sqrt_s) const;

    /**
     * @brief The STXS bin @f$qq \to H\ell\nu@f$, @f$N_j \geq 1,~150<p_{TV}[GeV]<250@f$.
     * @param[in] sqrt_s the center-of-mass energy in TeV
     */
    virtual const double STXS12_qqHlv_pTV150_250_Nj1(const double sqrt_s) const;

    /**
     * @brief The STXS bin @f$qq \to H\ell\nu@f$, @f$250<p_{TV}[GeV]@f$.
     * @param[in] sqrt_s the center-of-mass energy in TeV
     */
    virtual const double STXS12_qqHlv_pTV250_Inf(const double sqrt_s) const;

    /**
     * @brief The STXS bin @f$qq \to H\ell\nu@f$, @f$0<p_{TV}<150[GeV]@f$.
     * @param[in] sqrt_s the center-of-mass energy in TeV
     */
    virtual const double STXS12_qqHlv_pTV0_150(const double sqrt_s) const;            //AG:added
    
    /**
     * @brief The STXS bin @f$qq \to H\ell\nu@f$, @f$250<p_{TV}<400[GeV]@f$.
     * @param[in] sqrt_s the center-of-mass energy in TeV
     */
    virtual const double STXS12_qqHlv_pTV250_400(const double sqrt_s) const;          //AG:added
    
    /**
     * @brief The STXS bin @f$qq \to H\ell\nu@f$, @f$400<p_{TV}[GeV]@f$.
     * @param[in] sqrt_s the center-of-mass energy in TeV
     */
    virtual const double STXS12_qqHlv_pTV400_Inf(const double sqrt_s) const;          //AG:added    
    
    /**
     * @brief The STXS bin @f$qq \to H\ell\ell@f$, @f$p_{TV}[GeV]<75@f$.
     * @param[in] sqrt_s the center-of-mass energy in TeV
     */
    virtual const double STXS12_qqHll_pTV0_75(const double sqrt_s) const;

    /**
     * @brief The STXS bin @f$qq \to H\ell\ell@f$, @f$75<p_{TV}[GeV]<150@f$.
     * @param[in] sqrt_s the center-of-mass energy in TeV
     */
    virtual const double STXS12_qqHll_pTV75_150(const double sqrt_s) const;

    /**
     * @brief The STXS bin @f$qq \to H\ell\ell@f$, @f$N_j = 0,~150<p_{TV}[GeV]<250@f$.
     * @param[in] sqrt_s the center-of-mass energy in TeV
     */
    virtual const double STXS12_qqHll_pTV150_250_Nj0(const double sqrt_s) const;

    /**
     * @brief The STXS bin @f$qq \to H\ell\ell@f$, @f$N_j \geq 1,~150<p_{TV}[GeV]<250@f$.
     * @param[in] sqrt_s the center-of-mass energy in TeV
     */
    virtual const double STXS12_qqHll_pTV150_250_Nj1(const double sqrt_s) const;

    /**
     * @brief The STXS bin @f$qq \to H\ell\ell@f$, @f$250<p_{TV}[GeV]@f$.
     * @param[in] sqrt_s the center-of-mass energy in TeV
     */
    virtual const double STXS12_qqHll_pTV250_Inf(const double sqrt_s) const;

    /**
     * @brief The STXS bin @f$qq \to H\ell\ell@f$, @f$0<p_{TV}<150[GeV]@f$.
     * @param[in] sqrt_s the center-of-mass energy in TeV
     */
    virtual const double STXS12_qqHll_pTV0_150(const double sqrt_s) const;            //AG:added
    
    /**
     * @brief The STXS bin @f$qq \to H\ell\ell@f$, @f$250<p_{TV}<400[GeV]@f$.
     * @param[in] sqrt_s the center-of-mass energy in TeV
     */
    virtual const double STXS12_qqHll_pTV250_400(const double sqrt_s) const;          //AG:added
    
    /**
     * @brief The STXS bin @f$qq \to H\ell\ell@f$, @f$400<p_{TV}[GeV]@f$.
     * @param[in] sqrt_s the center-of-mass energy in TeV
     */
    virtual const double STXS12_qqHll_pTV400_Inf(const double sqrt_s) const;          //AG:added
    
    /**
     * @brief The STXS bin @f$pp \to ttH@f$, @f$p_{TH}[GeV]<60@f$.
     * @param[in] sqrt_s the center-of-mass energy in TeV
     */
    virtual const double STXS12_ttH_pTH0_60(const double sqrt_s) const;

    /**
     * @brief The STXS bin @f$pp \to ttH@f$, @f$60<p_{TH}[GeV]<120@f$.
     * @param[in] sqrt_s the center-of-mass energy in TeV
     */
    virtual const double STXS12_ttH_pTH60_120(const double sqrt_s) const;

    /**
     * @brief The STXS bin @f$pp \to ttH@f$, @f$120<p_{TH}[GeV]<200@f$.
     * @param[in] sqrt_s the center-of-mass energy in TeV
     */
    virtual const double STXS12_ttH_pTH120_200(const double sqrt_s) const;

    /**
     * @brief The STXS bin @f$pp \to ttH@f$, @f$200<p_{TH}[GeV]<300@f$.
     * @param[in] sqrt_s the center-of-mass energy in TeV
     */
    virtual const double STXS12_ttH_pTH200_300(const double sqrt_s) const;

    /**
     * @brief The STXS bin @f$pp \to ttH@f$, @f$300<p_{TH}[GeV]@f$.
     * @param[in] sqrt_s the center-of-mass energy in TeV
     */
    virtual const double STXS12_ttH_pTH300_Inf(const double sqrt_s) const;

    /**
     * @brief The STXS bin @f$pp \to ttH@f$, @f$300<p_{TH}[GeV]<450@f$.
     * @param[in] sqrt_s the center-of-mass energy in TeV
     */
    virtual const double STXS12_ttH_pTH300_450(const double sqrt_s) const;            //AG:added
    
    /**
     * @brief The STXS bin @f$pp \to ttH@f$, @f$450<p_{TH}[GeV]@f$.
     * @param[in] sqrt_s the center-of-mass energy in TeV
     */
    virtual const double STXS12_ttH_pTH450_Inf(const double sqrt_s) const;            //AG:added   
    
    /**
     * @brief The STXS bin @f$pp \to tH@f$.
     * @param[in] sqrt_s the center-of-mass energy in TeV
     */
    virtual const double STXS12_tH(const double sqrt_s) const;

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
     * @brief The effective coupling @f$\kappa_{s,eff}=\sqrt{\Gamma_{Hss}/\Gamma_{Hss}^{SM}}@f$.
     * @return @f$\kappa_{s,eff}@f$
     */
    virtual const double kappaseff() const;

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
     * @return @f$\delta y_t@f$
     */
    virtual const double deltayt_HB() const;

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
    virtual const double deltayb_HB() const;

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
    virtual const double deltaytau_HB() const;

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
    virtual const double deltayc_HB() const;
    

    /**
     * @brief The Higgs-basis coupling @f$\delta y_s@f$.
     * (See LHCHXSWG-INT-2015-001 document.)
     * Note that the Lagrangian definition of the Higgs-basis parameters coincides with the one of 
     * some of the @f$g_i, \delta g_i@f$ couplings defined above.
     * In the Higgs basis, however, one uses the freedom to perform certain field redefinitions and 
     * operations to demand that the mass eigenstate Lagrangian has specific features. (See pag. 5,6 in the reference.)
     * Therefore, the actual expression in terms of dim 6 coefficients may differ from the one for @f$g_i, \delta g_i@f$.
     * @return @f$\delta y_s@f$
     */
    virtual const double deltays_HB() const;

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
    virtual const double deltaymu_HB() const;

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
    virtual const double deltacZ_HB() const;

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
    virtual const double cZBox_HB() const;

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
    virtual const double cZZ_HB() const;

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
    virtual const double cZga_HB() const;

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
    virtual const double cgaga_HB() const;

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
    virtual const double cgg_HB() const;

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
    virtual const double cggEff_HB() const;

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
    virtual const double lambz_HB() const;

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

    const RGESolver& getSMEFTEvol() const {
        return SMEFTEvolEW;
    }

    gslpp::matrix<gslpp::complex> getYd() const {
        return Yd;
    }

    gslpp::matrix<gslpp::complex> getYe() const {
        return Ye;
    }

    gslpp::matrix<gslpp::complex> getYu() const {
        return Yu;
    }

    inline double getSMEFTCoeffEW(const std::string name) const
    {
        return SMEFTEvolEW.GetCoefficient(name);
    }
    
    inline double getSMEFTCoeffEW(const std::string name, int i, int j) const
    {
        return SMEFTEvolEW.GetCoefficient(name, i, j);
    }
    
    inline double getSMEFTCoeffEW(const std::string name, int i, int j, int k, int l) const
    {
        return SMEFTEvolEW.GetCoefficient(name, i, j, k, l);
    }
        
    void printNonVanishingSMEFTCoeffEW() const;
    
    RGESolver getSMEFTEvolEW() const
    {
        return SMEFTEvolEW;
    }

    
    // Functions to select the operators at different scales
    double getSMEFTCoeff(const std::string name, const double mu) const;
    
    double getSMEFTCoeff(const std::string name, int i, int j, const double mu) const;
    
    double getSMEFTCoeff(const std::string name, int i, int j, int k, int l, const double mu) const;
    
    
    const double getDelta_v() const
    {
        return delta_v;
    }

    const double getCW_tree() const
    {
        return cW_tree;
    }

    const double getSW_tree() const
    {
        return sW_tree;
    }

    const double getEeMz() const
    {
        return eeMz;
    }

    const double getGZdL() const
    {
        return gZdL;
    }

    const double getGZdR() const
    {
        return gZdR;
    }

    const double getGZlL() const
    {
        return gZlL;
    }

    const double getGZlR() const
    {
        return gZlR;
    }

    const double getGZuL() const
    {
        return gZuL;
    }

    const double getGZuR() const
    {
        return gZuR;
    }

    const double getGZvL() const
    {
        return gZvL;
    }

    const double getG1_tree() const
    {
        return g1_tree;
    }

    const double getG2_tree() const
    {
        return g2_tree;
    }

    const double getG3_tree() const
    {
        return g3_tree;
    }

    const double getDelta_Mz2() const
    {
        return delta_Mz2;
    }

    const double getDelta_ale() const
    {
        return delta_ale;
    }

    const double getDelta_g1() const
    {
        return delta_g1;
    }

    const double getDelta_g2() const
    {
        return delta_g2;
    }

    const double getDelta_xBZ() const
    {
        return delta_xBZ;
    }

    const double getDelta_xWZ() const
    {
        return delta_xWZ;
    }

    const double getXBZ_tree() const
    {
        return xBZ_tree;
    }

    const double getXWZ_tree() const
    {
        return xWZ_tree;
    }

    const double getMw_tree() const
    {
        return Mw_tree;
    }
    
    /**
     * getter for the rotation matrix for left-handed down-type quarks
     * @return a reference to the rotation matrix for left-handed down-type quarks
     */
    const gslpp::matrix<gslpp::complex>& getVdL() const {
        return VdL;
    }

    /**
     * getter for the Hermitean conjugate of the rotation matrix for left-handed down-type quarks
     * @return a reference to the Hermitean conjugate of the rotation matrix for left-handed down-type quarks
     */
    const gslpp::matrix<gslpp::complex>& getVdLd() const {
        return VdLd;
    }

    /**
     * getter for the rotation matrix for right-handed down-type quarks
     * @return a reference to the rotation matrix for right-handed down-type quarks
     */
    const gslpp::matrix<gslpp::complex>& getVdR() const {
        return VdR;
    }

    /**
     * getter for the Hermitean conjugate of the rotation matrix for right-handed down-type quarks
     * @return a reference to the Hermitean conjugate of the rotation matrix for right-handed down-type quarks
     */
    const gslpp::matrix<gslpp::complex>& getVdRd() const {
        return VdRd;
    }

    /**
     * getter for the rotation matrix for left-handed charged leptons
     * @return a reference to the rotation matrix for left-handed charged leptons
     */
    const gslpp::matrix<gslpp::complex>& getVeL() const {
        return VeL;
    }

    /**
     * getter for the Hermitean conjugate of the rotation matrix for left-handed charged leptons
     * @return a reference to the Hermitean conjugate of the rotation matrix for left-handed charged leptons
     */
    const gslpp::matrix<gslpp::complex>& getVeLd() const {
        return VeLd;
    }

    /**
     * getter for the rotation matrix for right-handed charged leptons
     * @return a reference to the rotation matrix for right-handed charged leptons
     */
    const gslpp::matrix<gslpp::complex>& getVeR() const {
        return VeR;
    }

    /**
     * getter for the Hermitean conjugate of the rotation matrix for right-handed charged leptons
     * @return a reference to the Hermitean conjugate of the rotation matrix for right-handed dcharged leptons
     */
    const gslpp::matrix<gslpp::complex>& getVeRd() const {
        return VeRd;
    }

    /**
     * getter for the rotation matrix for left-handed up-type quarks
     * @return a reference to the rotation matrix for left-handed up-type quarks
     */
    const gslpp::matrix<gslpp::complex>& getVuL() const {
        return VuL;
    }

    /**
     * getter for the Hermitean conjugate of the rotation matrix for left-handed up-type quarks
     * @return a reference to the Hermitean conjugate of the rotation matrix for left-handed up-type quarks
     */
    const gslpp::matrix<gslpp::complex>& getVuLd() const {
        return VuLd;
    }

    /**
     * getter for the rotation matrix for right-handed up-type quarks
     * @return a reference to the rotation matrix for right-handed up-type quarks
     */
    const gslpp::matrix<gslpp::complex>& getVuR() const {
        return VuR;
    }

    /**
     * getter for the Hermitean conjugate of the rotation matrix for right-handed up-type quarks
     * @return a reference to the Hermitean conjugate of the rotation matrix for right-handed up-type quarks
     */
    const gslpp::matrix<gslpp::complex>& getVuRd() const {
        return VuRd;
    }
   
    ////////////////////////////////////////////////////////////////////////
protected:

    //Old definition we won't use it
    //virtual void ComputeQuarkMassFromMCParameters(quark q);

    /**
     * @brief @copybrief Model::setParameter()
     * @copydetails Model::setParameter()
     */
    virtual void setParameter(const std::string name, const double& value);
    
    
    /**
     * @brief This function computes the SM parameters needed for the evolutor, neglecting any NP contribution.
     * 
     */
    void ChangeToEvolutorsBasisPureSM();
    
    
    /**
     * @brief This function computes the SM parameters needed for the evolutor, including a LO approximation to NP 
     * contributions where the Wilson coefficients are evaluated at the UV.
     * 
     */
    void ChangeToEvolutorsBasisSMEFTtoSM();
    
    
    /**
     * @brief The method to compute the %Yukawas matrix.
     */
    //virtual void computeYukawas();   TO BE FIXED ASAP!!!!!
    
    
    /**
     * @brief An auxiliary method to compute quark masses and CKM parameters from Yukawa couplings
     */
    void computeQuarkMassesAndCKMFromYukawas();
    
    
    /**
     * @brief An auxiliary method to set the WC on the evolutor
     */
    void setSMEFTEvolWC(RGESolver& RGevol);
    
    
    /**
     * @brief An auxiliary method to get the WC from the evolutor
     */
//    void getWCFromEvolutor();
    
    mutable Matching<NPSMEFTd6GeneralMatching, NPSMEFTd6General> NPSMEFTd6GM;

    std::string SMEFTBasisFlag;

    //For the SM parameters used in the evolutor we start at the EW scale
    double g1_LEW = 0;
    double g2_LEW = 0;
    double g3_LEW = 0;
    double lambdaH_LEW = 0;
    double mH2_LEW = 0;
    double me_LEW = 0;
    double mmu_LEW = 0;
    double mtau_LEW = 0;
    double mu_LEW = 0;
    double mc_LEW = 0;
    double mt_LEW = 0;
    double md_LEW = 0;
    double ms_LEW = 0;
    double mb_LEW = 0;
    double s12CKM_LEW = 0.;
    double s13CKM_LEW = 0.;
    double s23CKM_LEW = 0.;
    double dCKM_LEW = 0.;

    //We won't need the SM at the NP scale anymore
    //double g1_LNP = 0;
    //double g2_LNP = 0;
    //double g3_LNP = 0;
    //double lambdaH_LNP = 0;
    //double muH_LNP = 0;
    //double Ye_LNP = 0;
    //double Ymu_LNP = 0;
    //double Ytau_LNP = 0;
    //double Yu_LNP = 0;
    //double Yc_LNP = 0;
    //double Yt_LNP = 0;
    //double Yd_LNP = 0;
    //double Ys_LNP = 0;
    //double Yb_LNP = 0;
    //double s12CKM_LNP = 0.;
    //double s13CKM_LNP = 0.;
    //double s23CKM_LNP = 0.;
    //double dCKM_LNP = 0.;
    double CG_LNP = 0.; ///< The dimension-6 operator coefficient \f$C_{G}\f$.
    double CW_LNP = 0.; ///< The dimension-6 operator coefficient \f$C_{W}\f$.
    double CHG_LNP = 0.; ///< The dimension-6 operator coefficient \f$C_{HG}\f$.
    double CHW_LNP = 0.; ///< The dimension-6 operator coefficient \f$C_{HW}\f$.
    double CHB_LNP = 0.; ///< The dimension-6 operator coefficient \f$C_{HB}\f$.
    double CDB_LNP = 0; ///< The dimension-6 operator coefficient \f$C_{DB}\f$.
    double CDW_LNP = 0; ///< The dimension-6 operator coefficient \f$C_{DW}\f$.
    double CHWB_LNP = 0.; ///< The dimension-6 operator coefficient \f$C_{HWB}\f$.
    double CHD_LNP = 0.; ///< The dimension-6 operator coefficient \f$C_{HD}\f$.
    double CHbox_LNP = 0.; ///< The dimension-6 operator coefficient \f$C_{H\Box}\f$.
    double CH_LNP = 0.; ///< The dimension-6 operator coefficient \f$C_{H}\f$.
    double CGtilde_LNP = 0.; ///< The dimension-6 operator coefficient \f$C_{\tilde{G}}\f$.
    double CWtilde_LNP = 0.; ///< The dimension-6 operator coefficient \f$C_{\tilde{W}}\f$.
    double CHGtilde_LNP = 0.; ///< The dimension-6 operator coefficient \f$C_{H\tilde{G}}\f$.
    double CHWtilde_LNP = 0.; ///< The dimension-6 operator coefficient \f$C_{H\tilde{W}}\f$.
    double CHBtilde_LNP = 0.; ///< The dimension-6 operator coefficient \f$C_{H\tilde{B}}\f$.
    double CHWtildeB_LNP = 0.; ///< The dimension-6 operator coefficient \f$C_{H\tilde{W}B}\f$.
    double CHl1_11r_LNP = 0., CHl1_12r_LNP = 0., CHl1_13r_LNP = 0., CHl1_22r_LNP = 0., CHl1_23r_LNP = 0., CHl1_33r_LNP = 0.; ///< The dimension-6 operator coefficient \f$(C_{HL}^{(1)})_{ij}\f$ (Real part and pure real operator).
    double CHl1_12i_LNP = 0., CHl1_13i_LNP = 0., CHl1_23i_LNP = 0.; ///< The dimension-6 operator coefficient \f$(C_{HL}^{(1)})_{ij}\f$ (Imaginary part).
    double CHl3_11r_LNP = 0., CHl3_12r_LNP = 0., CHl3_13r_LNP = 0., CHl3_22r_LNP = 0., CHl3_23r_LNP = 0., CHl3_33r_LNP = 0.; ///< The dimension-6 operator coefficient \f$(C_{HL}^{(3)})_{ij}\f$ (Real part and pure real operator).
    double CHl3_12i_LNP = 0., CHl3_13i_LNP = 0., CHl3_23i_LNP = 0.; ///< The dimension-6 operator coefficient \f$(C_{HL}^{(3)})_{ij}\f$ (Imaginary part).
    double CHe_11r_LNP = 0., CHe_12r_LNP = 0., CHe_13r_LNP = 0., CHe_22r_LNP = 0., CHe_23r_LNP = 0., CHe_33r_LNP = 0.; ///< The dimension-6 operator coefficient \f$(C_{He})_{ij}\f$ (Real part and pure real operator).
    double CHe_12i_LNP = 0., CHe_13i_LNP = 0., CHe_23i_LNP = 0.; ///< The dimension-6 operator coefficient \f$(C_{He})_{ij}\f$ (Imaginary part).
    double CHq1_11r_LNP = 0., CHq1_12r_LNP = 0., CHq1_13r_LNP = 0., CHq1_22r_LNP = 0., CHq1_23r_LNP = 0., CHq1_33r_LNP = 0.; ///< The dimension-6 operator coefficient \f$(C_{HQ}^{(1)})_{ij}\f$ (Real part and pure real operator).
    double CHq1_12i_LNP = 0., CHq1_13i_LNP = 0., CHq1_23i_LNP = 0.; ///< The dimension-6 operator coefficient \f$(C_{HQ}^{(1)})_{ij}\f$ (Imaginary part).
    double CHq3_11r_LNP = 0., CHq3_12r_LNP = 0., CHq3_13r_LNP = 0., CHq3_22r_LNP = 0., CHq3_23r_LNP = 0., CHq3_33r_LNP = 0.; ///< The dimension-6 operator coefficient \f$(C_{HQ}^{(3)})_{ij}\f$ (Real part and pure real operator).
    double CHq3_12i_LNP = 0., CHq3_13i_LNP = 0., CHq3_23i_LNP = 0.; ///< The dimension-6 operator coefficient \f$(C_{HQ}^{(3)})_{ij}\f$ (Imaginary part).
    double CHu_11r_LNP = 0., CHu_12r_LNP = 0., CHu_13r_LNP = 0., CHu_22r_LNP = 0., CHu_23r_LNP = 0., CHu_33r_LNP = 0.; ///< The dimension-6 operator coefficient \f$(C_{Hu})_{ij}\f$ (Real part and pure real operator).
    double CHu_12i_LNP = 0., CHu_13i_LNP = 0., CHu_23i_LNP = 0.; ///< The dimension-6 operator coefficient \f$(C_{Hu})_{ij}\f$ (Imaginary part).
    double CHd_11r_LNP = 0., CHd_12r_LNP = 0., CHd_13r_LNP = 0., CHd_22r_LNP = 0., CHd_23r_LNP = 0., CHd_33r_LNP = 0.; ///< The dimension-6 operator coefficient \f$(C_{Hd})_{ij}\f$ (Real part and pure real operator).
    double CHd_12i_LNP = 0., CHd_13i_LNP = 0., CHd_23i_LNP = 0.; ///< The dimension-6 operator coefficient \f$(C_{Hd})_{ij}\f$ (Imaginary part).
    double CHud_11r_LNP = 0., CHud_12r_LNP = 0., CHud_13r_LNP = 0., CHud_21r_LNP = 0., CHud_22r_LNP = 0., CHud_23r_LNP = 0., CHud_31r_LNP = 0., CHud_32r_LNP = 0., CHud_33r_LNP = 0.; ///< The dimension-6 operator coefficient \f$(C_{Hud})_{ij}\f$ (Real part and pure real operator).
    double CHud_11i_LNP = 0., CHud_12i_LNP = 0., CHud_13i_LNP = 0., CHud_21i_LNP = 0., CHud_22i_LNP = 0., CHud_23i_LNP = 0., CHud_31i_LNP = 0., CHud_32i_LNP = 0., CHud_33i_LNP = 0.; ///< The dimension-6 operator coefficient \f$(C_{Hud})_{ij}\f$ (Imaginary part).
    double CeH_11r_LNP = 0., CeH_12r_LNP = 0., CeH_13r_LNP = 0., CeH_21r_LNP = 0., CeH_22r_LNP = 0., CeH_23r_LNP = 0., CeH_31r_LNP = 0., CeH_32r_LNP = 0., CeH_33r_LNP = 0.; ///< The dimension-6 operator coefficient \f$(C_{eH})_{ij}\f$ (Real part and pure real operator).
    double CeH_11i_LNP = 0., CeH_12i_LNP = 0., CeH_13i_LNP = 0., CeH_21i_LNP = 0., CeH_22i_LNP = 0., CeH_23i_LNP = 0., CeH_31i_LNP = 0., CeH_32i_LNP = 0., CeH_33i_LNP = 0.; ///< The dimension-6 operator coefficient \f$(C_{eH})_{ij}\f$ (Imaginary part).
    double CuH_11r_LNP = 0., CuH_12r_LNP = 0., CuH_13r_LNP = 0., CuH_21r_LNP = 0., CuH_22r_LNP = 0., CuH_23r_LNP = 0., CuH_31r_LNP = 0., CuH_32r_LNP = 0., CuH_33r_LNP = 0.; ///< The dimension-6 operator coefficient \f$(C_{uH})_{ij}\f$ (Real part and pure real operator).
    double CuH_11i_LNP = 0., CuH_12i_LNP = 0., CuH_13i_LNP = 0., CuH_21i_LNP = 0., CuH_22i_LNP = 0., CuH_23i_LNP = 0., CuH_31i_LNP = 0., CuH_32i_LNP = 0., CuH_33i_LNP = 0.; ///< The dimension-6 operator coefficient \f$(C_{uH})_{ij}\f$ (Imaginary part).
    double CdH_11r_LNP = 0., CdH_12r_LNP = 0., CdH_13r_LNP = 0., CdH_21r_LNP = 0., CdH_22r_LNP = 0., CdH_23r_LNP = 0., CdH_31r_LNP = 0., CdH_32r_LNP = 0., CdH_33r_LNP = 0.; ///< The dimension-6 operator coefficient \f$(C_{dH})_{ij}\f$ (Real part and pure real operator).
    double CdH_11i_LNP = 0., CdH_12i_LNP = 0., CdH_13i_LNP = 0., CdH_21i_LNP = 0., CdH_22i_LNP = 0., CdH_23i_LNP = 0., CdH_31i_LNP = 0., CdH_32i_LNP = 0., CdH_33i_LNP = 0.; ///< The dimension-6 operator coefficient \f$(C_{dH})_{ij}\f$ (Imaginary part).
    double CuG_11r_LNP = 0., CuG_12r_LNP = 0., CuG_13r_LNP = 0., CuG_21r_LNP = 0., CuG_22r_LNP = 0., CuG_23r_LNP = 0., CuG_31r_LNP = 0., CuG_32r_LNP = 0., CuG_33r_LNP = 0.; ///< The dimension-6 operator coefficient \f$(C_{uG})_{ij}\f$ (Real part and pure real operator).
    double CuG_11i_LNP = 0., CuG_12i_LNP = 0., CuG_13i_LNP = 0., CuG_21i_LNP = 0., CuG_22i_LNP = 0., CuG_23i_LNP = 0., CuG_31i_LNP = 0., CuG_32i_LNP = 0., CuG_33i_LNP = 0.; ///< The dimension-6 operator coefficient \f$(C_{uG})_{ij}\f$ (Imaginary part).
    double CuW_11r_LNP = 0., CuW_12r_LNP = 0., CuW_13r_LNP = 0., CuW_21r_LNP = 0., CuW_22r_LNP = 0., CuW_23r_LNP = 0., CuW_31r_LNP = 0., CuW_32r_LNP = 0., CuW_33r_LNP = 0.; ///< The dimension-6 operator coefficient \f$(C_{uW})_{ij}\f$ (Real part and pure real operator).
    double CuW_11i_LNP = 0., CuW_12i_LNP = 0., CuW_13i_LNP = 0., CuW_21i_LNP = 0., CuW_22i_LNP = 0., CuW_23i_LNP = 0., CuW_31i_LNP = 0., CuW_32i_LNP = 0., CuW_33i_LNP = 0.; ///< The dimension-6 operator coefficient \f$(C_{uW})_{ij}\f$ (Imaginary part).
    double CuB_11r_LNP = 0., CuB_12r_LNP = 0., CuB_13r_LNP = 0., CuB_21r_LNP = 0., CuB_22r_LNP = 0., CuB_23r_LNP = 0., CuB_31r_LNP = 0., CuB_32r_LNP = 0., CuB_33r_LNP = 0.; ///< The dimension-6 operator coefficient \f$(C_{uB})_{ij}\f$ (Real part and pure real operator).
    double CuB_11i_LNP = 0., CuB_12i_LNP = 0., CuB_13i_LNP = 0., CuB_21i_LNP = 0., CuB_22i_LNP = 0., CuB_23i_LNP = 0., CuB_31i_LNP = 0., CuB_32i_LNP = 0., CuB_33i_LNP = 0.; ///< The dimension-6 operator coefficient \f$(C_{uB})_{ij}\f$ (Imaginary part).
    double CdG_11r_LNP = 0., CdG_12r_LNP = 0., CdG_13r_LNP = 0., CdG_21r_LNP = 0., CdG_22r_LNP = 0., CdG_23r_LNP = 0., CdG_31r_LNP = 0., CdG_32r_LNP = 0., CdG_33r_LNP = 0.; ///< The dimension-6 operator coefficient \f$(C_{dG})_{ij}\f$ (Real part and pure real operator).
    double CdG_11i_LNP = 0., CdG_12i_LNP = 0., CdG_13i_LNP = 0., CdG_21i_LNP = 0., CdG_22i_LNP = 0., CdG_23i_LNP = 0., CdG_31i_LNP = 0., CdG_32i_LNP = 0., CdG_33i_LNP = 0.; ///< The dimension-6 operator coefficient \f$(C_{dG})_{ij}\f$ (Imaginary part).
    double CdW_11r_LNP = 0., CdW_12r_LNP = 0., CdW_13r_LNP = 0., CdW_21r_LNP = 0., CdW_22r_LNP = 0., CdW_23r_LNP = 0., CdW_31r_LNP = 0., CdW_32r_LNP = 0., CdW_33r_LNP = 0.; ///< The dimension-6 operator coefficient \f$(C_{dW})_{ij}\f$ (Real part and pure real operator).
    double CdW_11i_LNP = 0., CdW_12i_LNP = 0., CdW_13i_LNP = 0., CdW_21i_LNP = 0., CdW_22i_LNP = 0., CdW_23i_LNP = 0., CdW_31i_LNP = 0., CdW_32i_LNP = 0., CdW_33i_LNP = 0.; ///< The dimension-6 operator coefficient \f$(C_{dW})_{ij}\f$ (Imaginary part).
    double CdB_11r_LNP = 0., CdB_12r_LNP = 0., CdB_13r_LNP = 0., CdB_21r_LNP = 0., CdB_22r_LNP = 0., CdB_23r_LNP = 0., CdB_31r_LNP = 0., CdB_32r_LNP = 0., CdB_33r_LNP = 0.; ///< The dimension-6 operator coefficient \f$(C_{dB})_{ij}\f$ (Real part and pure real operator).
    double CdB_11i_LNP = 0., CdB_12i_LNP = 0., CdB_13i_LNP = 0., CdB_21i_LNP = 0., CdB_22i_LNP = 0., CdB_23i_LNP = 0., CdB_31i_LNP = 0., CdB_32i_LNP = 0., CdB_33i_LNP = 0.; ///< The dimension-6 operator coefficient \f$(C_{dB})_{ij}\f$ (Imaginary part).
    double CeW_11r_LNP = 0., CeW_12r_LNP = 0., CeW_13r_LNP = 0., CeW_21r_LNP = 0., CeW_22r_LNP = 0., CeW_23r_LNP = 0., CeW_31r_LNP = 0., CeW_32r_LNP = 0., CeW_33r_LNP = 0.; ///< The dimension-6 operator coefficient \f$(C_{eW})_{ij}\f$ (Real part and pure real operator).
    double CeW_11i_LNP = 0., CeW_12i_LNP = 0., CeW_13i_LNP = 0., CeW_21i_LNP = 0., CeW_22i_LNP = 0., CeW_23i_LNP = 0., CeW_31i_LNP = 0., CeW_32i_LNP = 0., CeW_33i_LNP = 0.; ///< The dimension-6 operator coefficient \f$(C_{eW})_{ij}\f$ (Imaginary part).
    double CeB_11r_LNP = 0., CeB_12r_LNP = 0., CeB_13r_LNP = 0., CeB_21r_LNP = 0., CeB_22r_LNP = 0., CeB_23r_LNP = 0., CeB_31r_LNP = 0., CeB_32r_LNP = 0., CeB_33r_LNP = 0.; ///< The dimension-6 operator coefficient \f$(C_{eB})_{ij}\f$ (Real part and pure real operator).
    double CeB_11i_LNP = 0., CeB_12i_LNP = 0., CeB_13i_LNP = 0., CeB_21i_LNP = 0., CeB_22i_LNP = 0., CeB_23i_LNP = 0., CeB_31i_LNP = 0., CeB_32i_LNP = 0., CeB_33i_LNP = 0.; ///< The dimension-6 operator coefficient \f$(C_{eB})_{ij}\f$ (Imaginary part).
    double Cll_1111r_LNP = 0., Cll_1112r_LNP = 0., Cll_1113r_LNP = 0., Cll_1122r_LNP = 0., Cll_1123r_LNP = 0., Cll_1133r_LNP = 0., Cll_1212r_LNP = 0., Cll_1213r_LNP = 0., Cll_1221r_LNP = 0.,
            Cll_1222r_LNP = 0., Cll_1223r_LNP = 0., Cll_1231r_LNP = 0., Cll_1232r_LNP = 0., Cll_1233r_LNP = 0., Cll_1313r_LNP = 0., Cll_1322r_LNP = 0., Cll_1323r_LNP = 0., Cll_1331r_LNP = 0.,
            Cll_1332r_LNP = 0., Cll_1333r_LNP = 0., Cll_2222r_LNP = 0., Cll_2223r_LNP = 0., Cll_2233r_LNP = 0., Cll_2323r_LNP = 0., Cll_2332r_LNP = 0., Cll_2333r_LNP = 0., Cll_3333r_LNP = 0.; ///< The dimension-6 operator coefficient \f$(C_{ll})_{ijkm}\f$ (Real part and pure real operator).
    double Cll_1112i_LNP = 0., Cll_1113i_LNP = 0., Cll_1123i_LNP = 0., Cll_1212i_LNP = 0., Cll_1213i_LNP = 0., Cll_1222i_LNP = 0., Cll_1223i_LNP = 0., Cll_1231i_LNP = 0., Cll_1232i_LNP = 0.,
            Cll_1233i_LNP = 0., Cll_1313i_LNP = 0., Cll_1322i_LNP = 0., Cll_1323i_LNP = 0., Cll_1332i_LNP = 0., Cll_1333i_LNP = 0., Cll_2223i_LNP = 0., Cll_2323i_LNP = 0., Cll_2333i_LNP = 0.; ///< The dimension-6 operator coefficient \f$(C_{ll})_{ijkm}\f$ (Imaginary part).
    double Clq1_1111r_LNP = 0., Clq1_1112r_LNP = 0., Clq1_1113r_LNP = 0., Clq1_1122r_LNP = 0., Clq1_1123r_LNP = 0., Clq1_1133r_LNP = 0., Clq1_1211r_LNP = 0., Clq1_1212r_LNP = 0., Clq1_1213r_LNP = 0.,
            Clq1_1221r_LNP = 0., Clq1_1222r_LNP = 0., Clq1_1223r_LNP = 0., Clq1_1231r_LNP = 0., Clq1_1232r_LNP = 0., Clq1_1233r_LNP = 0., Clq1_1311r_LNP = 0., Clq1_1312r_LNP = 0., Clq1_1313r_LNP = 0.,
            Clq1_1321r_LNP = 0., Clq1_1322r_LNP = 0., Clq1_1323r_LNP = 0., Clq1_1331r_LNP = 0., Clq1_1332r_LNP = 0., Clq1_1333r_LNP = 0., Clq1_2211r_LNP = 0., Clq1_2212r_LNP = 0., Clq1_2213r_LNP = 0.,
            Clq1_2222r_LNP = 0., Clq1_2223r_LNP = 0., Clq1_2233r_LNP = 0., Clq1_2311r_LNP = 0., Clq1_2312r_LNP = 0., Clq1_2313r_LNP = 0., Clq1_2321r_LNP = 0., Clq1_2322r_LNP = 0., Clq1_2323r_LNP = 0.,
            Clq1_2331r_LNP = 0., Clq1_2332r_LNP = 0., Clq1_2333r_LNP = 0., Clq1_3311r_LNP = 0., Clq1_3312r_LNP = 0., Clq1_3313r_LNP = 0., Clq1_3322r_LNP = 0., Clq1_3323r_LNP = 0., Clq1_3333r_LNP = 0.; ///< The dimension-6 operator coefficient \f$(C_{lq}^{(1)})_{ijkm}\f$ (Real part and pure real operator).
    double Clq1_1112i_LNP = 0., Clq1_1113i_LNP = 0., Clq1_1123i_LNP = 0., Clq1_1211i_LNP = 0., Clq1_1212i_LNP = 0., Clq1_1213i_LNP = 0., Clq1_1221i_LNP = 0., Clq1_1222i_LNP = 0., Clq1_1223i_LNP = 0.,
            Clq1_1231i_LNP = 0., Clq1_1232i_LNP = 0., Clq1_1233i_LNP = 0., Clq1_1311i_LNP = 0., Clq1_1312i_LNP = 0., Clq1_1313i_LNP = 0., Clq1_1321i_LNP = 0., Clq1_1322i_LNP = 0., Clq1_1323i_LNP = 0.,
            Clq1_1331i_LNP = 0., Clq1_1332i_LNP = 0., Clq1_1333i_LNP = 0., Clq1_2212i_LNP = 0., Clq1_2213i_LNP = 0., Clq1_2223i_LNP = 0., Clq1_2311i_LNP = 0., Clq1_2312i_LNP = 0., Clq1_2313i_LNP = 0., Clq1_2321i_LNP = 0.,
            Clq1_2322i_LNP = 0., Clq1_2323i_LNP = 0., Clq1_2331i_LNP = 0., Clq1_2332i_LNP = 0., Clq1_2333i_LNP = 0., Clq1_3312i_LNP = 0., Clq1_3313i_LNP = 0., Clq1_3323i_LNP = 0.; ///< The dimension-6 operator coefficient \f$(C_{lq}^{(1)})_{ijkm}\f$ (Imaginary part).
    double Clq3_1111r_LNP = 0., Clq3_1112r_LNP = 0., Clq3_1113r_LNP = 0., Clq3_1122r_LNP = 0., Clq3_1123r_LNP = 0., Clq3_1133r_LNP = 0., Clq3_1211r_LNP = 0., Clq3_1212r_LNP = 0., Clq3_1213r_LNP = 0.,
            Clq3_1221r_LNP = 0., Clq3_1222r_LNP = 0., Clq3_1223r_LNP = 0., Clq3_1231r_LNP = 0., Clq3_1232r_LNP = 0., Clq3_1233r_LNP = 0., Clq3_1311r_LNP = 0., Clq3_1312r_LNP = 0., Clq3_1313r_LNP = 0.,
            Clq3_1321r_LNP = 0., Clq3_1322r_LNP = 0., Clq3_1323r_LNP = 0., Clq3_1331r_LNP = 0., Clq3_1332r_LNP = 0., Clq3_1333r_LNP = 0., Clq3_2211r_LNP = 0., Clq3_2212r_LNP = 0., Clq3_2213r_LNP = 0.,
            Clq3_2222r_LNP = 0., Clq3_2223r_LNP = 0., Clq3_2233r_LNP = 0., Clq3_2311r_LNP = 0., Clq3_2312r_LNP = 0., Clq3_2313r_LNP = 0., Clq3_2321r_LNP = 0., Clq3_2322r_LNP = 0., Clq3_2323r_LNP = 0.,
            Clq3_2331r_LNP = 0., Clq3_2332r_LNP = 0., Clq3_2333r_LNP = 0., Clq3_3311r_LNP = 0., Clq3_3312r_LNP = 0., Clq3_3313r_LNP = 0., Clq3_3322r_LNP = 0., Clq3_3323r_LNP = 0., Clq3_3333r_LNP = 0.; ///< The dimension-6 operator coefficient \f$(C_{lq}^{(3)})_{ijkm}\f$ (Real part and pure real operator).
    double Clq3_1112i_LNP = 0., Clq3_1113i_LNP = 0., Clq3_1123i_LNP = 0., Clq3_1211i_LNP = 0., Clq3_1212i_LNP = 0., Clq3_1213i_LNP = 0., Clq3_1221i_LNP = 0., Clq3_1222i_LNP = 0., Clq3_1223i_LNP = 0.,
            Clq3_1231i_LNP = 0., Clq3_1232i_LNP = 0., Clq3_1233i_LNP = 0., Clq3_1311i_LNP = 0., Clq3_1312i_LNP = 0., Clq3_1313i_LNP = 0., Clq3_1321i_LNP = 0., Clq3_1322i_LNP = 0., Clq3_1323i_LNP = 0.,
            Clq3_1331i_LNP = 0., Clq3_1332i_LNP = 0., Clq3_1333i_LNP = 0., Clq3_2212i_LNP = 0., Clq3_2213i_LNP = 0., Clq3_2223i_LNP = 0., Clq3_2311i_LNP = 0., Clq3_2312i_LNP = 0., Clq3_2313i_LNP = 0., Clq3_2321i_LNP = 0.,
            Clq3_2322i_LNP = 0., Clq3_2323i_LNP = 0., Clq3_2331i_LNP = 0., Clq3_2332i_LNP = 0., Clq3_2333i_LNP = 0., Clq3_3312i_LNP = 0., Clq3_3313i_LNP = 0., Clq3_3323i_LNP = 0.; ///< The dimension-6 operator coefficient \f$(C_{lq}^{(3)})_{ijkm}\f$ (Imaginary part).
    double Cee_1111r_LNP = 0., Cee_1112r_LNP = 0., Cee_1113r_LNP = 0., Cee_1122r_LNP = 0., Cee_1123r_LNP = 0., Cee_1133r_LNP = 0., Cee_1212r_LNP = 0., Cee_1213r_LNP = 0., Cee_1222r_LNP = 0.,
            Cee_1223r_LNP = 0., Cee_1232r_LNP = 0., Cee_1233r_LNP = 0., Cee_1313r_LNP = 0., Cee_1323r_LNP = 0., Cee_1333r_LNP = 0., Cee_2222r_LNP = 0., Cee_2223r_LNP = 0., Cee_2233r_LNP = 0.,
            Cee_2323r_LNP = 0., Cee_2333r_LNP = 0., Cee_3333r_LNP = 0.; ///< The dimension-6 operator coefficient \f$(C_{ee})_{ijkm}\f$ (Real part and pure real operator).
    double Cee_1112i_LNP = 0., Cee_1113i_LNP = 0., Cee_1123i_LNP = 0., Cee_1212i_LNP = 0., Cee_1213i_LNP = 0., Cee_1222i_LNP = 0., Cee_1223i_LNP = 0., Cee_1232i_LNP = 0., Cee_1233i_LNP = 0.,
            Cee_1313i_LNP = 0., Cee_1323i_LNP = 0., Cee_1333i_LNP = 0., Cee_2223i_LNP = 0., Cee_2323i_LNP = 0., Cee_2333i_LNP = 0.; ///< The dimension-6 operator coefficient \f$(C_{ee})_{ijkm}\f$ (Imaginary part).
    double Ceu_1111r_LNP = 0., Ceu_1112r_LNP = 0., Ceu_1113r_LNP = 0., Ceu_1122r_LNP = 0., Ceu_1123r_LNP = 0., Ceu_1133r_LNP = 0., Ceu_1211r_LNP = 0., Ceu_1212r_LNP = 0., Ceu_1213r_LNP = 0.,
            Ceu_1221r_LNP = 0., Ceu_1222r_LNP = 0., Ceu_1223r_LNP = 0., Ceu_1231r_LNP = 0., Ceu_1232r_LNP = 0., Ceu_1233r_LNP = 0., Ceu_1311r_LNP = 0., Ceu_1312r_LNP = 0., Ceu_1313r_LNP = 0.,
            Ceu_1321r_LNP = 0., Ceu_1322r_LNP = 0., Ceu_1323r_LNP = 0., Ceu_1331r_LNP = 0., Ceu_1332r_LNP = 0., Ceu_1333r_LNP = 0., Ceu_2211r_LNP = 0., Ceu_2212r_LNP = 0., Ceu_2213r_LNP = 0.,
            Ceu_2222r_LNP = 0., Ceu_2223r_LNP = 0., Ceu_2233r_LNP = 0., Ceu_2311r_LNP = 0., Ceu_2312r_LNP = 0., Ceu_2313r_LNP = 0., Ceu_2321r_LNP = 0., Ceu_2322r_LNP = 0., Ceu_2323r_LNP = 0.,
            Ceu_2331r_LNP = 0., Ceu_2332r_LNP = 0., Ceu_2333r_LNP = 0., Ceu_3311r_LNP = 0., Ceu_3312r_LNP = 0., Ceu_3313r_LNP = 0., Ceu_3322r_LNP = 0., Ceu_3323r_LNP = 0., Ceu_3333r_LNP = 0.; ///< The dimension-6 operator coefficient \f$(C_{eu})_{ijkm}\f$ (Real part and pure real operator).
    double Ceu_1112i_LNP = 0., Ceu_1113i_LNP = 0., Ceu_1123i_LNP = 0., Ceu_1211i_LNP = 0., Ceu_1212i_LNP = 0., Ceu_1213i_LNP = 0., Ceu_1221i_LNP = 0., Ceu_1222i_LNP = 0., Ceu_1223i_LNP = 0.,
            Ceu_1231i_LNP = 0., Ceu_1232i_LNP = 0., Ceu_1233i_LNP = 0., Ceu_1311i_LNP = 0., Ceu_1312i_LNP = 0., Ceu_1313i_LNP = 0., Ceu_1321i_LNP = 0., Ceu_1322i_LNP = 0., Ceu_1323i_LNP = 0.,
            Ceu_1331i_LNP = 0., Ceu_1332i_LNP = 0., Ceu_1333i_LNP = 0., Ceu_2212i_LNP = 0., Ceu_2213i_LNP = 0., Ceu_2223i_LNP = 0., Ceu_2312i_LNP = 0., Ceu_2313i_LNP = 0., Ceu_2321i_LNP = 0.,
            Ceu_2322i_LNP = 0., Ceu_2323i_LNP = 0., Ceu_2331i_LNP = 0., Ceu_2332i_LNP = 0., Ceu_2333i_LNP = 0., Ceu_2311i_LNP = 0., Ceu_3312i_LNP = 0., Ceu_3313i_LNP = 0., Ceu_3323i_LNP = 0.; ///< The dimension-6 operator coefficient \f$(C_{eu})_{ijkm}\f$ (Imaginary part).
    double Ced_1111r_LNP = 0., Ced_1112r_LNP = 0., Ced_1113r_LNP = 0., Ced_1122r_LNP = 0., Ced_1123r_LNP = 0., Ced_1133r_LNP = 0., Ced_1211r_LNP = 0., Ced_1212r_LNP = 0., Ced_1213r_LNP = 0.,
            Ced_1221r_LNP = 0., Ced_1222r_LNP = 0., Ced_1223r_LNP = 0., Ced_1231r_LNP = 0., Ced_1232r_LNP = 0., Ced_1233r_LNP = 0., Ced_1311r_LNP = 0., Ced_1312r_LNP = 0., Ced_1313r_LNP = 0.,
            Ced_1321r_LNP = 0., Ced_1322r_LNP = 0., Ced_1323r_LNP = 0., Ced_1331r_LNP = 0., Ced_1332r_LNP = 0., Ced_1333r_LNP = 0., Ced_2211r_LNP = 0., Ced_2212r_LNP = 0., Ced_2213r_LNP = 0.,
            Ced_2222r_LNP = 0., Ced_2223r_LNP = 0., Ced_2233r_LNP = 0., Ced_2311r_LNP = 0., Ced_2312r_LNP = 0., Ced_2313r_LNP = 0., Ced_2321r_LNP = 0., Ced_2322r_LNP = 0., Ced_2323r_LNP = 0.,
            Ced_2331r_LNP = 0., Ced_2332r_LNP = 0., Ced_2333r_LNP = 0., Ced_3311r_LNP = 0., Ced_3312r_LNP = 0., Ced_3313r_LNP = 0., Ced_3322r_LNP = 0., Ced_3323r_LNP = 0., Ced_3333r_LNP = 0.; ///< The dimension-6 operator coefficient \f$(C_{ed})_{ijkm}\f$ (Real part and pure real operator).
    double Ced_1112i_LNP = 0., Ced_1113i_LNP = 0., Ced_1123i_LNP = 0., Ced_1211i_LNP = 0., Ced_1212i_LNP = 0., Ced_1213i_LNP = 0., Ced_1221i_LNP = 0., Ced_1222i_LNP = 0., Ced_1223i_LNP = 0.,
            Ced_1231i_LNP = 0., Ced_1232i_LNP = 0., Ced_1233i_LNP = 0., Ced_1311i_LNP = 0., Ced_1312i_LNP = 0., Ced_1313i_LNP = 0., Ced_1321i_LNP = 0., Ced_1322i_LNP = 0., Ced_1323i_LNP = 0.,
            Ced_1331i_LNP = 0., Ced_1332i_LNP = 0., Ced_1333i_LNP = 0., Ced_2212i_LNP = 0., Ced_2213i_LNP = 0., Ced_2223i_LNP = 0., Ced_2312i_LNP = 0., Ced_2313i_LNP = 0., Ced_2321i_LNP = 0.,
            Ced_2322i_LNP = 0., Ced_2323i_LNP = 0., Ced_2331i_LNP = 0., Ced_2332i_LNP = 0., Ced_2333i_LNP = 0., Ced_2311i_LNP = 0., Ced_3312i_LNP = 0., Ced_3313i_LNP = 0., Ced_3323i_LNP = 0.; ///< The dimension-6 operator coefficient \f$(C_{ed})_{ijkm}\f$ (Imaginary part).
    double Cle_1111r_LNP = 0., Cle_1112r_LNP = 0., Cle_1113r_LNP = 0., Cle_1122r_LNP = 0., Cle_1123r_LNP = 0., Cle_1133r_LNP = 0., Cle_1211r_LNP = 0., Cle_1212r_LNP = 0., Cle_1213r_LNP = 0.,
            Cle_1221r_LNP = 0., Cle_1222r_LNP = 0., Cle_1223r_LNP = 0., Cle_1231r_LNP = 0., Cle_1232r_LNP = 0., Cle_1233r_LNP = 0., Cle_1311r_LNP = 0., Cle_1312r_LNP = 0., Cle_1313r_LNP = 0.,
            Cle_1321r_LNP = 0., Cle_1322r_LNP = 0., Cle_1323r_LNP = 0., Cle_1331r_LNP = 0., Cle_1332r_LNP = 0., Cle_1333r_LNP = 0., Cle_2211r_LNP = 0., Cle_2212r_LNP = 0., Cle_2213r_LNP = 0.,
            Cle_2222r_LNP = 0., Cle_2223r_LNP = 0., Cle_2233r_LNP = 0., Cle_2311r_LNP = 0., Cle_2312r_LNP = 0., Cle_2313r_LNP = 0., Cle_2321r_LNP = 0., Cle_2322r_LNP = 0., Cle_2323r_LNP = 0.,
            Cle_2331r_LNP = 0., Cle_2332r_LNP = 0., Cle_2333r_LNP = 0., Cle_3311r_LNP = 0., Cle_3312r_LNP = 0., Cle_3313r_LNP = 0., Cle_3322r_LNP = 0., Cle_3323r_LNP = 0., Cle_3333r_LNP = 0.; ///< The dimension-6 operator coefficient \f$(C_{le})_{ijkm}\f$ (Real part and pure real operator).
    double Cle_1112i_LNP = 0., Cle_1113i_LNP = 0., Cle_1123i_LNP = 0., Cle_1211i_LNP = 0., Cle_1212i_LNP = 0., Cle_1213i_LNP = 0., Cle_1221i_LNP = 0., Cle_1222i_LNP = 0., Cle_1223i_LNP = 0.,
            Cle_1231i_LNP = 0., Cle_1232i_LNP = 0., Cle_1233i_LNP = 0., Cle_1311i_LNP = 0., Cle_1312i_LNP = 0., Cle_1313i_LNP = 0., Cle_1321i_LNP = 0., Cle_1322i_LNP = 0., Cle_1323i_LNP = 0.,
            Cle_1331i_LNP = 0., Cle_1332i_LNP = 0., Cle_1333i_LNP = 0., Cle_2212i_LNP = 0., Cle_2213i_LNP = 0., Cle_2223i_LNP = 0., Cle_2312i_LNP = 0., Cle_2313i_LNP = 0., Cle_2321i_LNP = 0.,
            Cle_2322i_LNP = 0., Cle_2323i_LNP = 0., Cle_2331i_LNP = 0., Cle_2332i_LNP = 0., Cle_2333i_LNP = 0., Cle_2311i_LNP = 0., Cle_3312i_LNP = 0., Cle_3313i_LNP = 0., Cle_3323i_LNP = 0.; ///< The dimension-6 operator coefficient \f$(C_{le})_{ijkm}\f$ (Imaginary part).
    double Clu_1111r_LNP = 0., Clu_1112r_LNP = 0., Clu_1113r_LNP = 0., Clu_1122r_LNP = 0., Clu_1123r_LNP = 0., Clu_1133r_LNP = 0., Clu_1211r_LNP = 0., Clu_1212r_LNP = 0., Clu_1213r_LNP = 0.,
            Clu_1221r_LNP = 0., Clu_1222r_LNP = 0., Clu_1223r_LNP = 0., Clu_1231r_LNP = 0., Clu_1232r_LNP = 0., Clu_1233r_LNP = 0., Clu_1311r_LNP = 0., Clu_1312r_LNP = 0., Clu_1313r_LNP = 0.,
            Clu_1321r_LNP = 0., Clu_1322r_LNP = 0., Clu_1323r_LNP = 0., Clu_1331r_LNP = 0., Clu_1332r_LNP = 0., Clu_1333r_LNP = 0., Clu_2211r_LNP = 0., Clu_2212r_LNP = 0., Clu_2213r_LNP = 0.,
            Clu_2222r_LNP = 0., Clu_2223r_LNP = 0., Clu_2233r_LNP = 0., Clu_2311r_LNP = 0., Clu_2312r_LNP = 0., Clu_2313r_LNP = 0., Clu_2321r_LNP = 0., Clu_2322r_LNP = 0., Clu_2323r_LNP = 0.,
            Clu_2331r_LNP = 0., Clu_2332r_LNP = 0., Clu_2333r_LNP = 0., Clu_3311r_LNP = 0., Clu_3312r_LNP = 0., Clu_3313r_LNP = 0., Clu_3322r_LNP = 0., Clu_3323r_LNP = 0., Clu_3333r_LNP = 0.; ///< The dimension-6 operator coefficient \f$(C_{lu})_{ijkm}\f$ (Real part and pure real operator).
    double Clu_1112i_LNP = 0., Clu_1113i_LNP = 0., Clu_1123i_LNP = 0., Clu_1211i_LNP = 0., Clu_1212i_LNP = 0., Clu_1213i_LNP = 0., Clu_1221i_LNP = 0., Clu_1222i_LNP = 0., Clu_1223i_LNP = 0.,
            Clu_1231i_LNP = 0., Clu_1232i_LNP = 0., Clu_1233i_LNP = 0., Clu_1311i_LNP = 0., Clu_1312i_LNP = 0., Clu_1313i_LNP = 0., Clu_1321i_LNP = 0., Clu_1322i_LNP = 0., Clu_1323i_LNP = 0.,
            Clu_1331i_LNP = 0., Clu_1332i_LNP = 0., Clu_1333i_LNP = 0., Clu_2212i_LNP = 0., Clu_2213i_LNP = 0., Clu_2223i_LNP = 0., Clu_2312i_LNP = 0., Clu_2313i_LNP = 0., Clu_2321i_LNP = 0.,
            Clu_2322i_LNP = 0., Clu_2323i_LNP = 0., Clu_2331i_LNP = 0., Clu_2332i_LNP = 0., Clu_2333i_LNP = 0., Clu_2311i_LNP = 0., Clu_3312i_LNP = 0., Clu_3313i_LNP = 0., Clu_3323i_LNP = 0.; ///< The dimension-6 operator coefficient \f$(C_{lu})_{ijkm}\f$ (Imaginary part).
    double Cld_1111r_LNP = 0., Cld_1112r_LNP = 0., Cld_1113r_LNP = 0., Cld_1122r_LNP = 0., Cld_1123r_LNP = 0., Cld_1133r_LNP = 0., Cld_1211r_LNP = 0., Cld_1212r_LNP = 0., Cld_1213r_LNP = 0.,
            Cld_1221r_LNP = 0., Cld_1222r_LNP = 0., Cld_1223r_LNP = 0., Cld_1231r_LNP = 0., Cld_1232r_LNP = 0., Cld_1233r_LNP = 0., Cld_1311r_LNP = 0., Cld_1312r_LNP = 0., Cld_1313r_LNP = 0.,
            Cld_1321r_LNP = 0., Cld_1322r_LNP = 0., Cld_1323r_LNP = 0., Cld_1331r_LNP = 0., Cld_1332r_LNP = 0., Cld_1333r_LNP = 0., Cld_2211r_LNP = 0., Cld_2212r_LNP = 0., Cld_2213r_LNP = 0.,
            Cld_2222r_LNP = 0., Cld_2223r_LNP = 0., Cld_2233r_LNP = 0., Cld_2311r_LNP = 0., Cld_2312r_LNP = 0., Cld_2313r_LNP = 0., Cld_2321r_LNP = 0., Cld_2322r_LNP = 0., Cld_2323r_LNP = 0.,
            Cld_2331r_LNP = 0., Cld_2332r_LNP = 0., Cld_2333r_LNP = 0., Cld_3311r_LNP = 0., Cld_3312r_LNP = 0., Cld_3313r_LNP = 0., Cld_3322r_LNP = 0., Cld_3323r_LNP = 0., Cld_3333r_LNP = 0.; ///< The dimension-6 operator coefficient \f$(C_{ld})_{ijkm}\f$ (Real part and pure real operator).
    double Cld_1112i_LNP = 0., Cld_1113i_LNP = 0., Cld_1123i_LNP = 0., Cld_1211i_LNP = 0., Cld_1212i_LNP = 0., Cld_1213i_LNP = 0., Cld_1221i_LNP = 0., Cld_1222i_LNP = 0., Cld_1223i_LNP = 0.,
            Cld_1231i_LNP = 0., Cld_1232i_LNP = 0., Cld_1233i_LNP = 0., Cld_1311i_LNP = 0., Cld_1312i_LNP = 0., Cld_1313i_LNP = 0., Cld_1321i_LNP = 0., Cld_1322i_LNP = 0., Cld_1323i_LNP = 0.,
            Cld_1331i_LNP = 0., Cld_1332i_LNP = 0., Cld_1333i_LNP = 0., Cld_2212i_LNP = 0., Cld_2213i_LNP = 0., Cld_2223i_LNP = 0., Cld_2312i_LNP = 0., Cld_2313i_LNP = 0., Cld_2321i_LNP = 0.,
            Cld_2322i_LNP = 0., Cld_2323i_LNP = 0., Cld_2331i_LNP = 0., Cld_2332i_LNP = 0., Cld_2333i_LNP = 0., Cld_2311i_LNP = 0., Cld_3312i_LNP = 0., Cld_3313i_LNP = 0., Cld_3323i_LNP = 0.; ///< The dimension-6 operator coefficient \f$(C_{ld})_{ijkm}\f$ (Imaginary part).
    double Cqe_1111r_LNP = 0., Cqe_1112r_LNP = 0., Cqe_1113r_LNP = 0., Cqe_1122r_LNP = 0., Cqe_1123r_LNP = 0., Cqe_1133r_LNP = 0., Cqe_1211r_LNP = 0., Cqe_1212r_LNP = 0., Cqe_1213r_LNP = 0.,
            Cqe_1221r_LNP = 0., Cqe_1222r_LNP = 0., Cqe_1223r_LNP = 0., Cqe_1231r_LNP = 0., Cqe_1232r_LNP = 0., Cqe_1233r_LNP = 0., Cqe_1311r_LNP = 0., Cqe_1312r_LNP = 0., Cqe_1313r_LNP = 0.,
            Cqe_1321r_LNP = 0., Cqe_1322r_LNP = 0., Cqe_1323r_LNP = 0., Cqe_1331r_LNP = 0., Cqe_1332r_LNP = 0., Cqe_1333r_LNP = 0., Cqe_2211r_LNP = 0., Cqe_2212r_LNP = 0., Cqe_2213r_LNP = 0.,
            Cqe_2222r_LNP = 0., Cqe_2223r_LNP = 0., Cqe_2233r_LNP = 0., Cqe_2311r_LNP = 0., Cqe_2312r_LNP = 0., Cqe_2313r_LNP = 0., Cqe_2321r_LNP = 0., Cqe_2322r_LNP = 0., Cqe_2323r_LNP = 0.,
            Cqe_2331r_LNP = 0., Cqe_2332r_LNP = 0., Cqe_2333r_LNP = 0., Cqe_3311r_LNP = 0., Cqe_3312r_LNP = 0., Cqe_3313r_LNP = 0., Cqe_3322r_LNP = 0., Cqe_3323r_LNP = 0., Cqe_3333r_LNP = 0.; ///< The dimension-6 operator coefficient \f$(C_{qe})_{ijkm}\f$ (Real part and pure real operator).
    double Cqe_1112i_LNP = 0., Cqe_1113i_LNP = 0., Cqe_1123i_LNP = 0., Cqe_1211i_LNP = 0., Cqe_1212i_LNP = 0., Cqe_1213i_LNP = 0., Cqe_1221i_LNP = 0., Cqe_1222i_LNP = 0., Cqe_1223i_LNP = 0.,
            Cqe_1231i_LNP = 0., Cqe_1232i_LNP = 0., Cqe_1233i_LNP = 0., Cqe_1311i_LNP = 0., Cqe_1312i_LNP = 0., Cqe_1313i_LNP = 0., Cqe_1321i_LNP = 0., Cqe_1322i_LNP = 0., Cqe_1323i_LNP = 0.,
            Cqe_1331i_LNP = 0., Cqe_1332i_LNP = 0., Cqe_1333i_LNP = 0., Cqe_2212i_LNP = 0., Cqe_2213i_LNP = 0., Cqe_2223i_LNP = 0., Cqe_2312i_LNP = 0., Cqe_2313i_LNP = 0., Cqe_2321i_LNP = 0.,
            Cqe_2322i_LNP = 0., Cqe_2323i_LNP = 0., Cqe_2331i_LNP = 0., Cqe_2332i_LNP = 0., Cqe_2333i_LNP = 0., Cqe_2311i_LNP = 0., Cqe_3312i_LNP = 0., Cqe_3313i_LNP = 0., Cqe_3323i_LNP = 0.; ///< The dimension-6 operator coefficient \f$(C_{qe})_{ijkm}\f$ (Imaginary part).
    double Cledq_1111r_LNP = 0., Cledq_1112r_LNP = 0., Cledq_1113r_LNP = 0., Cledq_1121r_LNP = 0., Cledq_1122r_LNP = 0., Cledq_1123r_LNP = 0., Cledq_1131r_LNP = 0., Cledq_1132r_LNP = 0., Cledq_1133r_LNP = 0.,
            Cledq_1211r_LNP = 0., Cledq_1212r_LNP = 0., Cledq_1213r_LNP = 0., Cledq_1221r_LNP = 0., Cledq_1222r_LNP = 0., Cledq_1223r_LNP = 0., Cledq_1231r_LNP = 0., Cledq_1232r_LNP = 0., Cledq_1233r_LNP = 0.,
            Cledq_1311r_LNP = 0., Cledq_1312r_LNP = 0., Cledq_1313r_LNP = 0., Cledq_1321r_LNP = 0., Cledq_1322r_LNP = 0., Cledq_1323r_LNP = 0., Cledq_1331r_LNP = 0., Cledq_1332r_LNP = 0., Cledq_1333r_LNP = 0.,
            Cledq_2111r_LNP = 0., Cledq_2112r_LNP = 0., Cledq_2113r_LNP = 0., Cledq_2121r_LNP = 0., Cledq_2122r_LNP = 0., Cledq_2123r_LNP = 0., Cledq_2131r_LNP = 0., Cledq_2132r_LNP = 0., Cledq_2133r_LNP = 0.,
            Cledq_2211r_LNP = 0., Cledq_2212r_LNP = 0., Cledq_2213r_LNP = 0., Cledq_2221r_LNP = 0., Cledq_2222r_LNP = 0., Cledq_2223r_LNP = 0., Cledq_2231r_LNP = 0., Cledq_2232r_LNP = 0., Cledq_2233r_LNP = 0.,
            Cledq_2311r_LNP = 0., Cledq_2312r_LNP = 0., Cledq_2313r_LNP = 0., Cledq_2321r_LNP = 0., Cledq_2322r_LNP = 0., Cledq_2323r_LNP = 0., Cledq_2331r_LNP = 0., Cledq_2332r_LNP = 0., Cledq_2333r_LNP = 0.,
            Cledq_3111r_LNP = 0., Cledq_3112r_LNP = 0., Cledq_3113r_LNP = 0., Cledq_3121r_LNP = 0., Cledq_3122r_LNP = 0., Cledq_3123r_LNP = 0., Cledq_3131r_LNP = 0., Cledq_3132r_LNP = 0., Cledq_3133r_LNP = 0.,
            Cledq_3211r_LNP = 0., Cledq_3212r_LNP = 0., Cledq_3213r_LNP = 0., Cledq_3221r_LNP = 0., Cledq_3222r_LNP = 0., Cledq_3223r_LNP = 0., Cledq_3231r_LNP = 0., Cledq_3232r_LNP = 0., Cledq_3233r_LNP = 0.,
            Cledq_3311r_LNP = 0., Cledq_3312r_LNP = 0., Cledq_3313r_LNP = 0., Cledq_3321r_LNP = 0., Cledq_3322r_LNP = 0., Cledq_3323r_LNP = 0., Cledq_3331r_LNP = 0., Cledq_3332r_LNP = 0., Cledq_3333r_LNP = 0.; ///< The dimension-6 operator coefficient \f$(C_{ledq})_{ijkm}\f$ (Real part and pure real operator).
    double Cledq_1111i_LNP = 0., Cledq_1112i_LNP = 0., Cledq_1113i_LNP = 0., Cledq_1121i_LNP = 0., Cledq_1122i_LNP = 0., Cledq_1123i_LNP = 0., Cledq_1131i_LNP = 0., Cledq_1132i_LNP = 0., Cledq_1133i_LNP = 0.,
            Cledq_1211i_LNP = 0., Cledq_1212i_LNP = 0., Cledq_1213i_LNP = 0., Cledq_1221i_LNP = 0., Cledq_1222i_LNP = 0., Cledq_1223i_LNP = 0., Cledq_1231i_LNP = 0., Cledq_1232i_LNP = 0., Cledq_1233i_LNP = 0.,
            Cledq_1311i_LNP = 0., Cledq_1312i_LNP = 0., Cledq_1313i_LNP = 0., Cledq_1321i_LNP = 0., Cledq_1322i_LNP = 0., Cledq_1323i_LNP = 0., Cledq_1331i_LNP = 0., Cledq_1332i_LNP = 0., Cledq_1333i_LNP = 0.,
            Cledq_2111i_LNP = 0., Cledq_2112i_LNP = 0., Cledq_2113i_LNP = 0., Cledq_2121i_LNP = 0., Cledq_2122i_LNP = 0., Cledq_2123i_LNP = 0., Cledq_2131i_LNP = 0., Cledq_2132i_LNP = 0., Cledq_2133i_LNP = 0.,
            Cledq_2211i_LNP = 0., Cledq_2212i_LNP = 0., Cledq_2213i_LNP = 0., Cledq_2221i_LNP = 0., Cledq_2222i_LNP = 0., Cledq_2223i_LNP = 0., Cledq_2231i_LNP = 0., Cledq_2232i_LNP = 0., Cledq_2233i_LNP = 0.,
            Cledq_2311i_LNP = 0., Cledq_2312i_LNP = 0., Cledq_2313i_LNP = 0., Cledq_2321i_LNP = 0., Cledq_2322i_LNP = 0., Cledq_2323i_LNP = 0., Cledq_2331i_LNP = 0., Cledq_2332i_LNP = 0., Cledq_2333i_LNP = 0.,
            Cledq_3111i_LNP = 0., Cledq_3112i_LNP = 0., Cledq_3113i_LNP = 0., Cledq_3121i_LNP = 0., Cledq_3122i_LNP = 0., Cledq_3123i_LNP = 0., Cledq_3131i_LNP = 0., Cledq_3132i_LNP = 0., Cledq_3133i_LNP = 0.,
            Cledq_3211i_LNP = 0., Cledq_3212i_LNP = 0., Cledq_3213i_LNP = 0., Cledq_3221i_LNP = 0., Cledq_3222i_LNP = 0., Cledq_3223i_LNP = 0., Cledq_3231i_LNP = 0., Cledq_3232i_LNP = 0., Cledq_3233i_LNP = 0.,
            Cledq_3311i_LNP = 0., Cledq_3312i_LNP = 0., Cledq_3313i_LNP = 0., Cledq_3321i_LNP = 0., Cledq_3322i_LNP = 0., Cledq_3323i_LNP = 0., Cledq_3331i_LNP = 0., Cledq_3332i_LNP = 0., Cledq_3333i_LNP = 0.; ///< The dimension-6 operator coefficient \f$(C_{ledq})_{ijkm}\f$ (Imaginary part).
    double Cqq1_1111r_LNP = 0., Cqq1_1112r_LNP = 0., Cqq1_1113r_LNP = 0., Cqq1_1122r_LNP = 0., Cqq1_1123r_LNP = 0., Cqq1_1133r_LNP = 0., Cqq1_1212r_LNP = 0., Cqq1_1213r_LNP = 0., Cqq1_1221r_LNP = 0.,
            Cqq1_1222r_LNP = 0., Cqq1_1223r_LNP = 0., Cqq1_1231r_LNP = 0., Cqq1_1232r_LNP = 0., Cqq1_1233r_LNP = 0., Cqq1_1313r_LNP = 0., Cqq1_1322r_LNP = 0., Cqq1_1323r_LNP = 0., Cqq1_1331r_LNP = 0.,
            Cqq1_1332r_LNP = 0., Cqq1_1333r_LNP = 0., Cqq1_2222r_LNP = 0., Cqq1_2223r_LNP = 0., Cqq1_2233r_LNP = 0., Cqq1_2323r_LNP = 0., Cqq1_2332r_LNP = 0., Cqq1_2333r_LNP = 0., Cqq1_3333r_LNP = 0.; ///< The dimension-6 operator coefficient \f$(C_{qq}^{(1)})_{ijkm}\f$ (Real part and pure real operator).
    double Cqq1_1112i_LNP = 0., Cqq1_1113i_LNP = 0., Cqq1_1123i_LNP = 0., Cqq1_1212i_LNP = 0., Cqq1_1213i_LNP = 0., Cqq1_1222i_LNP = 0., Cqq1_1223i_LNP = 0., Cqq1_1231i_LNP = 0., Cqq1_1232i_LNP = 0.,
            Cqq1_1233i_LNP = 0., Cqq1_1313i_LNP = 0., Cqq1_1322i_LNP = 0., Cqq1_1323i_LNP = 0., Cqq1_1332i_LNP = 0., Cqq1_1333i_LNP = 0., Cqq1_2223i_LNP = 0., Cqq1_2323i_LNP = 0., Cqq1_2333i_LNP = 0.; ///< The dimension-6 operator coefficient \f$(C_{qq}^{(1)})_{ijkm}\f$ (Imaginary part).
    double Cqq3_1111r_LNP = 0., Cqq3_1112r_LNP = 0., Cqq3_1113r_LNP = 0., Cqq3_1122r_LNP = 0., Cqq3_1123r_LNP = 0., Cqq3_1133r_LNP = 0., Cqq3_1212r_LNP = 0., Cqq3_1213r_LNP = 0., Cqq3_1221r_LNP = 0.,
            Cqq3_1222r_LNP = 0., Cqq3_1223r_LNP = 0., Cqq3_1231r_LNP = 0., Cqq3_1232r_LNP = 0., Cqq3_1233r_LNP = 0., Cqq3_1313r_LNP = 0., Cqq3_1322r_LNP = 0., Cqq3_1323r_LNP = 0., Cqq3_1331r_LNP = 0.,
            Cqq3_1332r_LNP = 0., Cqq3_1333r_LNP = 0., Cqq3_2222r_LNP = 0., Cqq3_2223r_LNP = 0., Cqq3_2233r_LNP = 0., Cqq3_2323r_LNP = 0., Cqq3_2332r_LNP = 0., Cqq3_2333r_LNP = 0., Cqq3_3333r_LNP = 0.; ///< The dimension-6 operator coefficient \f$(C_{qq}^{(3)})_{ijkm}\f$ (Real part and pure real operator).
    double Cqq3_1112i_LNP = 0., Cqq3_1113i_LNP = 0., Cqq3_1123i_LNP = 0., Cqq3_1212i_LNP = 0., Cqq3_1213i_LNP = 0., Cqq3_1222i_LNP = 0., Cqq3_1223i_LNP = 0., Cqq3_1231i_LNP = 0., Cqq3_1232i_LNP = 0.,
            Cqq3_1233i_LNP = 0., Cqq3_1313i_LNP = 0., Cqq3_1322i_LNP = 0., Cqq3_1323i_LNP = 0., Cqq3_1332i_LNP = 0., Cqq3_1333i_LNP = 0., Cqq3_2223i_LNP = 0., Cqq3_2323i_LNP = 0., Cqq3_2333i_LNP = 0.; ///< The dimension-6 operator coefficient \f$(C_{qq}^{(3)})_{ijkm}\f$ (Imaginary part).
    double Cuu_1111r_LNP = 0., Cuu_1112r_LNP = 0., Cuu_1113r_LNP = 0., Cuu_1122r_LNP = 0., Cuu_1123r_LNP = 0., Cuu_1133r_LNP = 0., Cuu_1212r_LNP = 0., Cuu_1213r_LNP = 0., Cuu_1221r_LNP = 0.,
            Cuu_1222r_LNP = 0., Cuu_1223r_LNP = 0., Cuu_1231r_LNP = 0., Cuu_1232r_LNP = 0., Cuu_1233r_LNP = 0., Cuu_1313r_LNP = 0., Cuu_1322r_LNP = 0., Cuu_1323r_LNP = 0., Cuu_1331r_LNP = 0.,
            Cuu_1332r_LNP = 0., Cuu_1333r_LNP = 0., Cuu_2222r_LNP = 0., Cuu_2223r_LNP = 0., Cuu_2233r_LNP = 0., Cuu_2323r_LNP = 0., Cuu_2332r_LNP = 0., Cuu_2333r_LNP = 0., Cuu_3333r_LNP = 0.; ///< The dimension-6 operator coefficient \f$(C_{uu})_{ijkm}\f$ (Real part and pure real operator).
    double Cuu_1112i_LNP = 0., Cuu_1113i_LNP = 0., Cuu_1123i_LNP = 0., Cuu_1212i_LNP = 0., Cuu_1213i_LNP = 0., Cuu_1222i_LNP = 0., Cuu_1223i_LNP = 0., Cuu_1231i_LNP = 0., Cuu_1232i_LNP = 0.,
            Cuu_1233i_LNP = 0., Cuu_1313i_LNP = 0., Cuu_1322i_LNP = 0., Cuu_1323i_LNP = 0., Cuu_1332i_LNP = 0., Cuu_1333i_LNP = 0., Cuu_2223i_LNP = 0., Cuu_2323i_LNP = 0., Cuu_2333i_LNP = 0.; ///< The dimension-6 operator coefficient \f$(C_{uu})_{ijkm}\f$ (Imaginary part).
    double Cdd_1111r_LNP = 0., Cdd_1112r_LNP = 0., Cdd_1113r_LNP = 0., Cdd_1122r_LNP = 0., Cdd_1123r_LNP = 0., Cdd_1133r_LNP = 0., Cdd_1212r_LNP = 0., Cdd_1213r_LNP = 0., Cdd_1221r_LNP = 0.,
            Cdd_1222r_LNP = 0., Cdd_1223r_LNP = 0., Cdd_1231r_LNP = 0., Cdd_1232r_LNP = 0., Cdd_1233r_LNP = 0., Cdd_1313r_LNP = 0., Cdd_1322r_LNP = 0., Cdd_1323r_LNP = 0., Cdd_1331r_LNP = 0.,
            Cdd_1332r_LNP = 0., Cdd_1333r_LNP = 0., Cdd_2222r_LNP = 0., Cdd_2223r_LNP = 0., Cdd_2233r_LNP = 0., Cdd_2323r_LNP = 0., Cdd_2332r_LNP = 0., Cdd_2333r_LNP = 0., Cdd_3333r_LNP = 0.; ///< The dimension-6 operator coefficient \f$(C_{dd})_{ijkm}\f$ (Real part and pure real operator).
    double Cdd_1112i_LNP = 0., Cdd_1113i_LNP = 0., Cdd_1123i_LNP = 0., Cdd_1212i_LNP = 0., Cdd_1213i_LNP = 0., Cdd_1222i_LNP = 0., Cdd_1223i_LNP = 0., Cdd_1231i_LNP = 0., Cdd_1232i_LNP = 0.,
            Cdd_1233i_LNP = 0., Cdd_1313i_LNP = 0., Cdd_1322i_LNP = 0., Cdd_1323i_LNP = 0., Cdd_1332i_LNP = 0., Cdd_1333i_LNP = 0., Cdd_2223i_LNP = 0., Cdd_2323i_LNP = 0., Cdd_2333i_LNP = 0.; ///< The dimension-6 operator coefficient \f$(C_{dd})_{ijkm}\f$ (Imaginary part).
    double Cud1_1111r_LNP = 0., Cud1_1112r_LNP = 0., Cud1_1113r_LNP = 0., Cud1_1122r_LNP = 0., Cud1_1123r_LNP = 0., Cud1_1133r_LNP = 0., Cud1_1211r_LNP = 0., Cud1_1212r_LNP = 0., Cud1_1213r_LNP = 0.,
            Cud1_1221r_LNP = 0., Cud1_1222r_LNP = 0., Cud1_1223r_LNP = 0., Cud1_1231r_LNP = 0., Cud1_1232r_LNP = 0., Cud1_1233r_LNP = 0., Cud1_1311r_LNP = 0., Cud1_1312r_LNP = 0., Cud1_1313r_LNP = 0.,
            Cud1_1321r_LNP = 0., Cud1_1322r_LNP = 0., Cud1_1323r_LNP = 0., Cud1_1331r_LNP = 0., Cud1_1332r_LNP = 0., Cud1_1333r_LNP = 0., Cud1_2211r_LNP = 0., Cud1_2212r_LNP = 0., Cud1_2213r_LNP = 0.,
            Cud1_2222r_LNP = 0., Cud1_2223r_LNP = 0., Cud1_2233r_LNP = 0., Cud1_2311r_LNP = 0., Cud1_2312r_LNP = 0., Cud1_2313r_LNP = 0., Cud1_2321r_LNP = 0., Cud1_2322r_LNP = 0., Cud1_2323r_LNP = 0.,
            Cud1_2331r_LNP = 0., Cud1_2332r_LNP = 0., Cud1_2333r_LNP = 0., Cud1_3311r_LNP = 0., Cud1_3312r_LNP = 0., Cud1_3313r_LNP = 0., Cud1_3322r_LNP = 0., Cud1_3323r_LNP = 0., Cud1_3333r_LNP = 0.; ///< The dimension-6 operator coefficient \f$(C_{ud}^{(1)})_{ijkm}\f$ (Real part and pure real operator).
    double Cud1_1112i_LNP = 0., Cud1_1113i_LNP = 0., Cud1_1123i_LNP = 0., Cud1_1211i_LNP = 0., Cud1_1212i_LNP = 0., Cud1_1213i_LNP = 0., Cud1_1221i_LNP = 0., Cud1_1222i_LNP = 0., Cud1_1223i_LNP = 0.,
            Cud1_1231i_LNP = 0., Cud1_1232i_LNP = 0., Cud1_1233i_LNP = 0., Cud1_1311i_LNP = 0., Cud1_1312i_LNP = 0., Cud1_1313i_LNP = 0., Cud1_1321i_LNP = 0., Cud1_1322i_LNP = 0., Cud1_1323i_LNP = 0.,
            Cud1_1331i_LNP = 0., Cud1_1332i_LNP = 0., Cud1_1333i_LNP = 0., Cud1_2212i_LNP = 0., Cud1_2213i_LNP = 0., Cud1_2223i_LNP = 0., Cud1_2312i_LNP = 0., Cud1_2313i_LNP = 0., Cud1_2321i_LNP = 0.,
            Cud1_2322i_LNP = 0., Cud1_2323i_LNP = 0., Cud1_2331i_LNP = 0., Cud1_2332i_LNP = 0., Cud1_2333i_LNP = 0., Cud1_2311i_LNP = 0., Cud1_3312i_LNP = 0., Cud1_3313i_LNP = 0., Cud1_3323i_LNP = 0.; ///< The dimension-6 operator coefficient \f$(C_{ud}^{(1)})_{ijkm}\f$ (Imaginary part).
    double Cud8_1111r_LNP = 0., Cud8_1112r_LNP = 0., Cud8_1113r_LNP = 0., Cud8_1122r_LNP = 0., Cud8_1123r_LNP = 0., Cud8_1133r_LNP = 0., Cud8_1211r_LNP = 0., Cud8_1212r_LNP = 0., Cud8_1213r_LNP = 0.,
            Cud8_1221r_LNP = 0., Cud8_1222r_LNP = 0., Cud8_1223r_LNP = 0., Cud8_1231r_LNP = 0., Cud8_1232r_LNP = 0., Cud8_1233r_LNP = 0., Cud8_1311r_LNP = 0., Cud8_1312r_LNP = 0., Cud8_1313r_LNP = 0.,
            Cud8_1321r_LNP = 0., Cud8_1322r_LNP = 0., Cud8_1323r_LNP = 0., Cud8_1331r_LNP = 0., Cud8_1332r_LNP = 0., Cud8_1333r_LNP = 0., Cud8_2211r_LNP = 0., Cud8_2212r_LNP = 0., Cud8_2213r_LNP = 0.,
            Cud8_2222r_LNP = 0., Cud8_2223r_LNP = 0., Cud8_2233r_LNP = 0., Cud8_2311r_LNP = 0., Cud8_2312r_LNP = 0., Cud8_2313r_LNP = 0., Cud8_2321r_LNP = 0., Cud8_2322r_LNP = 0., Cud8_2323r_LNP = 0.,
            Cud8_2331r_LNP = 0., Cud8_2332r_LNP = 0., Cud8_2333r_LNP = 0., Cud8_3311r_LNP = 0., Cud8_3312r_LNP = 0., Cud8_3313r_LNP = 0., Cud8_3322r_LNP = 0., Cud8_3323r_LNP = 0., Cud8_3333r_LNP = 0.; ///< The dimension-6 operator coefficient \f$(C_{ud}^{(8)})_{ijkm}\f$ (Real part and pure real operator).
    double Cud8_1112i_LNP = 0., Cud8_1113i_LNP = 0., Cud8_1123i_LNP = 0., Cud8_1211i_LNP = 0., Cud8_1212i_LNP = 0., Cud8_1213i_LNP = 0., Cud8_1221i_LNP = 0., Cud8_1222i_LNP = 0., Cud8_1223i_LNP = 0.,
            Cud8_1231i_LNP = 0., Cud8_1232i_LNP = 0., Cud8_1233i_LNP = 0., Cud8_1311i_LNP = 0., Cud8_1312i_LNP = 0., Cud8_1313i_LNP = 0., Cud8_1321i_LNP = 0., Cud8_1322i_LNP = 0., Cud8_1323i_LNP = 0.,
            Cud8_1331i_LNP = 0., Cud8_1332i_LNP = 0., Cud8_1333i_LNP = 0., Cud8_2212i_LNP = 0., Cud8_2213i_LNP = 0., Cud8_2223i_LNP = 0., Cud8_2312i_LNP = 0., Cud8_2313i_LNP = 0., Cud8_2321i_LNP = 0.,
            Cud8_2322i_LNP = 0., Cud8_2323i_LNP = 0., Cud8_2331i_LNP = 0., Cud8_2332i_LNP = 0., Cud8_2333i_LNP = 0., Cud8_2311i_LNP = 0., Cud8_3312i_LNP = 0., Cud8_3313i_LNP = 0., Cud8_3323i_LNP = 0.; ///< The dimension-6 operator coefficient \f$(C_{ud}^{(8)})_{ijkm}\f$ (Imaginary part).
    double Cqu1_1111r_LNP = 0., Cqu1_1112r_LNP = 0., Cqu1_1113r_LNP = 0., Cqu1_1122r_LNP = 0., Cqu1_1123r_LNP = 0., Cqu1_1133r_LNP = 0., Cqu1_1211r_LNP = 0., Cqu1_1212r_LNP = 0., Cqu1_1213r_LNP = 0.,
            Cqu1_1221r_LNP = 0., Cqu1_1222r_LNP = 0., Cqu1_1223r_LNP = 0., Cqu1_1231r_LNP = 0., Cqu1_1232r_LNP = 0., Cqu1_1233r_LNP = 0., Cqu1_1311r_LNP = 0., Cqu1_1312r_LNP = 0., Cqu1_1313r_LNP = 0.,
            Cqu1_1321r_LNP = 0., Cqu1_1322r_LNP = 0., Cqu1_1323r_LNP = 0., Cqu1_1331r_LNP = 0., Cqu1_1332r_LNP = 0., Cqu1_1333r_LNP = 0., Cqu1_2211r_LNP = 0., Cqu1_2212r_LNP = 0., Cqu1_2213r_LNP = 0.,
            Cqu1_2222r_LNP = 0., Cqu1_2223r_LNP = 0., Cqu1_2233r_LNP = 0., Cqu1_2311r_LNP = 0., Cqu1_2312r_LNP = 0., Cqu1_2313r_LNP = 0., Cqu1_2321r_LNP = 0., Cqu1_2322r_LNP = 0., Cqu1_2323r_LNP = 0.,
            Cqu1_2331r_LNP = 0., Cqu1_2332r_LNP = 0., Cqu1_2333r_LNP = 0., Cqu1_3311r_LNP = 0., Cqu1_3312r_LNP = 0., Cqu1_3313r_LNP = 0., Cqu1_3322r_LNP = 0., Cqu1_3323r_LNP = 0., Cqu1_3333r_LNP = 0.; ///< The dimension-6 operator coefficient \f$(C_{qu}^{(1)})_{ijkm}\f$ (Real part and pure real operator).
    double Cqu1_1112i_LNP = 0., Cqu1_1113i_LNP = 0., Cqu1_1123i_LNP = 0., Cqu1_1211i_LNP = 0., Cqu1_1212i_LNP = 0., Cqu1_1213i_LNP = 0., Cqu1_1221i_LNP = 0., Cqu1_1222i_LNP = 0., Cqu1_1223i_LNP = 0.,
            Cqu1_1231i_LNP = 0., Cqu1_1232i_LNP = 0., Cqu1_1233i_LNP = 0., Cqu1_1311i_LNP = 0., Cqu1_1312i_LNP = 0., Cqu1_1313i_LNP = 0., Cqu1_1321i_LNP = 0., Cqu1_1322i_LNP = 0., Cqu1_1323i_LNP = 0.,
            Cqu1_1331i_LNP = 0., Cqu1_1332i_LNP = 0., Cqu1_1333i_LNP = 0., Cqu1_2212i_LNP = 0., Cqu1_2213i_LNP = 0., Cqu1_2223i_LNP = 0., Cqu1_2312i_LNP = 0., Cqu1_2313i_LNP = 0., Cqu1_2321i_LNP = 0.,
            Cqu1_2322i_LNP = 0., Cqu1_2323i_LNP = 0., Cqu1_2331i_LNP = 0., Cqu1_2332i_LNP = 0., Cqu1_2333i_LNP = 0., Cqu1_2311i_LNP = 0., Cqu1_3312i_LNP = 0., Cqu1_3313i_LNP = 0., Cqu1_3323i_LNP = 0.; ///< The dimension-6 operator coefficient \f$(C_{qu}^{(1)})_{ijkm}\f$ (Imaginary part).
    double Cqu8_1111r_LNP = 0., Cqu8_1112r_LNP = 0., Cqu8_1113r_LNP = 0., Cqu8_1122r_LNP = 0., Cqu8_1123r_LNP = 0., Cqu8_1133r_LNP = 0., Cqu8_1211r_LNP = 0., Cqu8_1212r_LNP = 0., Cqu8_1213r_LNP = 0.,
            Cqu8_1221r_LNP = 0., Cqu8_1222r_LNP = 0., Cqu8_1223r_LNP = 0., Cqu8_1231r_LNP = 0., Cqu8_1232r_LNP = 0., Cqu8_1233r_LNP = 0., Cqu8_1311r_LNP = 0., Cqu8_1312r_LNP = 0., Cqu8_1313r_LNP = 0.,
            Cqu8_1321r_LNP = 0., Cqu8_1322r_LNP = 0., Cqu8_1323r_LNP = 0., Cqu8_1331r_LNP = 0., Cqu8_1332r_LNP = 0., Cqu8_1333r_LNP = 0., Cqu8_2211r_LNP = 0., Cqu8_2212r_LNP = 0., Cqu8_2213r_LNP = 0.,
            Cqu8_2222r_LNP = 0., Cqu8_2223r_LNP = 0., Cqu8_2233r_LNP = 0., Cqu8_2311r_LNP = 0., Cqu8_2312r_LNP = 0., Cqu8_2313r_LNP = 0., Cqu8_2321r_LNP = 0., Cqu8_2322r_LNP = 0., Cqu8_2323r_LNP = 0.,
            Cqu8_2331r_LNP = 0., Cqu8_2332r_LNP = 0., Cqu8_2333r_LNP = 0., Cqu8_3311r_LNP = 0., Cqu8_3312r_LNP = 0., Cqu8_3313r_LNP = 0., Cqu8_3322r_LNP = 0., Cqu8_3323r_LNP = 0., Cqu8_3333r_LNP = 0.; ///< The dimension-6 operator coefficient \f$(C_{qu}^{(8)})_{ijkm}\f$ (Real part and pure real operator).
    double Cqu8_1112i_LNP = 0., Cqu8_1113i_LNP = 0., Cqu8_1123i_LNP = 0., Cqu8_1211i_LNP = 0., Cqu8_1212i_LNP = 0., Cqu8_1213i_LNP = 0., Cqu8_1221i_LNP = 0., Cqu8_1222i_LNP = 0., Cqu8_1223i_LNP = 0.,
            Cqu8_1231i_LNP = 0., Cqu8_1232i_LNP = 0., Cqu8_1233i_LNP = 0., Cqu8_1311i_LNP = 0., Cqu8_1312i_LNP = 0., Cqu8_1313i_LNP = 0., Cqu8_1321i_LNP = 0., Cqu8_1322i_LNP = 0., Cqu8_1323i_LNP = 0.,
            Cqu8_1331i_LNP = 0., Cqu8_1332i_LNP = 0., Cqu8_1333i_LNP = 0., Cqu8_2212i_LNP = 0., Cqu8_2213i_LNP = 0., Cqu8_2223i_LNP = 0., Cqu8_2312i_LNP = 0., Cqu8_2313i_LNP = 0., Cqu8_2321i_LNP = 0.,
            Cqu8_2322i_LNP = 0., Cqu8_2323i_LNP = 0., Cqu8_2331i_LNP = 0., Cqu8_2332i_LNP = 0., Cqu8_2333i_LNP = 0., Cqu8_2311i_LNP = 0., Cqu8_3312i_LNP = 0., Cqu8_3313i_LNP = 0., Cqu8_3323i_LNP = 0.; ///< The dimension-6 operator coefficient \f$(C_{qu}^{(8)})_{ijkm}\f$ (Imaginary part).
    double Cqd1_1111r_LNP = 0., Cqd1_1112r_LNP = 0., Cqd1_1113r_LNP = 0., Cqd1_1122r_LNP = 0., Cqd1_1123r_LNP = 0., Cqd1_1133r_LNP = 0., Cqd1_1211r_LNP = 0., Cqd1_1212r_LNP = 0., Cqd1_1213r_LNP = 0.,
            Cqd1_1221r_LNP = 0., Cqd1_1222r_LNP = 0., Cqd1_1223r_LNP = 0., Cqd1_1231r_LNP = 0., Cqd1_1232r_LNP = 0., Cqd1_1233r_LNP = 0., Cqd1_1311r_LNP = 0., Cqd1_1312r_LNP = 0., Cqd1_1313r_LNP = 0.,
            Cqd1_1321r_LNP = 0., Cqd1_1322r_LNP = 0., Cqd1_1323r_LNP = 0., Cqd1_1331r_LNP = 0., Cqd1_1332r_LNP = 0., Cqd1_1333r_LNP = 0., Cqd1_2211r_LNP = 0., Cqd1_2212r_LNP = 0., Cqd1_2213r_LNP = 0.,
            Cqd1_2222r_LNP = 0., Cqd1_2223r_LNP = 0., Cqd1_2233r_LNP = 0., Cqd1_2311r_LNP = 0., Cqd1_2312r_LNP = 0., Cqd1_2313r_LNP = 0., Cqd1_2321r_LNP = 0., Cqd1_2322r_LNP = 0., Cqd1_2323r_LNP = 0.,
            Cqd1_2331r_LNP = 0., Cqd1_2332r_LNP = 0., Cqd1_2333r_LNP = 0., Cqd1_3311r_LNP = 0., Cqd1_3312r_LNP = 0., Cqd1_3313r_LNP = 0., Cqd1_3322r_LNP = 0., Cqd1_3323r_LNP = 0., Cqd1_3333r_LNP = 0.; ///< The dimension-6 operator coefficient \f$(C_{qd}^{(1)})_{ijkm}\f$ (Real part and pure real operator).
    double Cqd1_1112i_LNP = 0., Cqd1_1113i_LNP = 0., Cqd1_1123i_LNP = 0., Cqd1_1211i_LNP = 0., Cqd1_1212i_LNP = 0., Cqd1_1213i_LNP = 0., Cqd1_1221i_LNP = 0., Cqd1_1222i_LNP = 0., Cqd1_1223i_LNP = 0.,
            Cqd1_1231i_LNP = 0., Cqd1_1232i_LNP = 0., Cqd1_1233i_LNP = 0., Cqd1_1311i_LNP = 0., Cqd1_1312i_LNP = 0., Cqd1_1313i_LNP = 0., Cqd1_1321i_LNP = 0., Cqd1_1322i_LNP = 0., Cqd1_1323i_LNP = 0.,
            Cqd1_1331i_LNP = 0., Cqd1_1332i_LNP = 0., Cqd1_1333i_LNP = 0., Cqd1_2212i_LNP = 0., Cqd1_2213i_LNP = 0., Cqd1_2223i_LNP = 0., Cqd1_2312i_LNP = 0., Cqd1_2313i_LNP = 0., Cqd1_2321i_LNP = 0.,
            Cqd1_2322i_LNP = 0., Cqd1_2323i_LNP = 0., Cqd1_2331i_LNP = 0., Cqd1_2332i_LNP = 0., Cqd1_2333i_LNP = 0., Cqd1_2311i_LNP = 0., Cqd1_3312i_LNP = 0., Cqd1_3313i_LNP = 0., Cqd1_3323i_LNP = 0.; ///< The dimension-6 operator coefficient \f$(C_{qd}^{(1)})_{ijkm}\f$ (Imaginary part).
    double Cqd8_1111r_LNP = 0., Cqd8_1112r_LNP = 0., Cqd8_1113r_LNP = 0., Cqd8_1122r_LNP = 0., Cqd8_1123r_LNP = 0., Cqd8_1133r_LNP = 0., Cqd8_1211r_LNP = 0., Cqd8_1212r_LNP = 0., Cqd8_1213r_LNP = 0.,
            Cqd8_1221r_LNP = 0., Cqd8_1222r_LNP = 0., Cqd8_1223r_LNP = 0., Cqd8_1231r_LNP = 0., Cqd8_1232r_LNP = 0., Cqd8_1233r_LNP = 0., Cqd8_1311r_LNP = 0., Cqd8_1312r_LNP = 0., Cqd8_1313r_LNP = 0.,
            Cqd8_1321r_LNP = 0., Cqd8_1322r_LNP = 0., Cqd8_1323r_LNP = 0., Cqd8_1331r_LNP = 0., Cqd8_1332r_LNP = 0., Cqd8_1333r_LNP = 0., Cqd8_2211r_LNP = 0., Cqd8_2212r_LNP = 0., Cqd8_2213r_LNP = 0.,
            Cqd8_2222r_LNP = 0., Cqd8_2223r_LNP = 0., Cqd8_2233r_LNP = 0., Cqd8_2311r_LNP = 0., Cqd8_2312r_LNP = 0., Cqd8_2313r_LNP = 0., Cqd8_2321r_LNP = 0., Cqd8_2322r_LNP = 0., Cqd8_2323r_LNP = 0.,
            Cqd8_2331r_LNP = 0., Cqd8_2332r_LNP = 0., Cqd8_2333r_LNP = 0., Cqd8_3311r_LNP = 0., Cqd8_3312r_LNP = 0., Cqd8_3313r_LNP = 0., Cqd8_3322r_LNP = 0., Cqd8_3323r_LNP = 0., Cqd8_3333r_LNP = 0.; ///< The dimension-6 operator coefficient \f$(C_{qd}^{(8)})_{ijkm}\f$ (Real part and pure real operator).
    double Cqd8_1112i_LNP = 0., Cqd8_1113i_LNP = 0., Cqd8_1123i_LNP = 0., Cqd8_1211i_LNP = 0., Cqd8_1212i_LNP = 0., Cqd8_1213i_LNP = 0., Cqd8_1221i_LNP = 0., Cqd8_1222i_LNP = 0., Cqd8_1223i_LNP = 0.,
            Cqd8_1231i_LNP = 0., Cqd8_1232i_LNP = 0., Cqd8_1233i_LNP = 0., Cqd8_1311i_LNP = 0., Cqd8_1312i_LNP = 0., Cqd8_1313i_LNP = 0., Cqd8_1321i_LNP = 0., Cqd8_1322i_LNP = 0., Cqd8_1323i_LNP = 0.,
            Cqd8_1331i_LNP = 0., Cqd8_1332i_LNP = 0., Cqd8_1333i_LNP = 0., Cqd8_2212i_LNP = 0., Cqd8_2213i_LNP = 0., Cqd8_2223i_LNP = 0., Cqd8_2312i_LNP = 0., Cqd8_2313i_LNP = 0., Cqd8_2321i_LNP = 0.,
            Cqd8_2322i_LNP = 0., Cqd8_2323i_LNP = 0., Cqd8_2331i_LNP = 0., Cqd8_2332i_LNP = 0., Cqd8_2333i_LNP = 0., Cqd8_2311i_LNP = 0., Cqd8_3312i_LNP = 0., Cqd8_3313i_LNP = 0., Cqd8_3323i_LNP = 0.; ///< The dimension-6 operator coefficient \f$(C_{qd}^{(8)})_{ijkm}\f$ (Imaginary part).
    double Cquqd1_1111r_LNP = 0., Cquqd1_1112r_LNP = 0., Cquqd1_1113r_LNP = 0., Cquqd1_1121r_LNP = 0., Cquqd1_1122r_LNP = 0., Cquqd1_1123r_LNP = 0., Cquqd1_1131r_LNP = 0., Cquqd1_1132r_LNP = 0., Cquqd1_1133r_LNP = 0.,
            Cquqd1_1211r_LNP = 0., Cquqd1_1212r_LNP = 0., Cquqd1_1213r_LNP = 0., Cquqd1_1221r_LNP = 0., Cquqd1_1222r_LNP = 0., Cquqd1_1223r_LNP = 0., Cquqd1_1231r_LNP = 0., Cquqd1_1232r_LNP = 0., Cquqd1_1233r_LNP = 0.,
            Cquqd1_1311r_LNP = 0., Cquqd1_1312r_LNP = 0., Cquqd1_1313r_LNP = 0., Cquqd1_1321r_LNP = 0., Cquqd1_1322r_LNP = 0., Cquqd1_1323r_LNP = 0., Cquqd1_1331r_LNP = 0., Cquqd1_1332r_LNP = 0., Cquqd1_1333r_LNP = 0.,
            Cquqd1_2111r_LNP = 0., Cquqd1_2112r_LNP = 0., Cquqd1_2113r_LNP = 0., Cquqd1_2121r_LNP = 0., Cquqd1_2122r_LNP = 0., Cquqd1_2123r_LNP = 0., Cquqd1_2131r_LNP = 0., Cquqd1_2132r_LNP = 0., Cquqd1_2133r_LNP = 0.,
            Cquqd1_2211r_LNP = 0., Cquqd1_2212r_LNP = 0., Cquqd1_2213r_LNP = 0., Cquqd1_2221r_LNP = 0., Cquqd1_2222r_LNP = 0., Cquqd1_2223r_LNP = 0., Cquqd1_2231r_LNP = 0., Cquqd1_2232r_LNP = 0., Cquqd1_2233r_LNP = 0.,
            Cquqd1_2311r_LNP = 0., Cquqd1_2312r_LNP = 0., Cquqd1_2313r_LNP = 0., Cquqd1_2321r_LNP = 0., Cquqd1_2322r_LNP = 0., Cquqd1_2323r_LNP = 0., Cquqd1_2331r_LNP = 0., Cquqd1_2332r_LNP = 0., Cquqd1_2333r_LNP = 0.,
            Cquqd1_3111r_LNP = 0., Cquqd1_3112r_LNP = 0., Cquqd1_3113r_LNP = 0., Cquqd1_3121r_LNP = 0., Cquqd1_3122r_LNP = 0., Cquqd1_3123r_LNP = 0., Cquqd1_3131r_LNP = 0., Cquqd1_3132r_LNP = 0., Cquqd1_3133r_LNP = 0.,
            Cquqd1_3211r_LNP = 0., Cquqd1_3212r_LNP = 0., Cquqd1_3213r_LNP = 0., Cquqd1_3221r_LNP = 0., Cquqd1_3222r_LNP = 0., Cquqd1_3223r_LNP = 0., Cquqd1_3231r_LNP = 0., Cquqd1_3232r_LNP = 0., Cquqd1_3233r_LNP = 0.,
            Cquqd1_3311r_LNP = 0., Cquqd1_3312r_LNP = 0., Cquqd1_3313r_LNP = 0., Cquqd1_3321r_LNP = 0., Cquqd1_3322r_LNP = 0., Cquqd1_3323r_LNP = 0., Cquqd1_3331r_LNP = 0., Cquqd1_3332r_LNP = 0., Cquqd1_3333r_LNP = 0.; ///< The dimension-6 operator coefficient \f$(C_{quqd}^{(1)})_{ijkm}\f$ (Real part and pure real operator).
    double Cquqd1_1111i_LNP = 0., Cquqd1_1112i_LNP = 0., Cquqd1_1113i_LNP = 0., Cquqd1_1121i_LNP = 0., Cquqd1_1122i_LNP = 0., Cquqd1_1123i_LNP = 0., Cquqd1_1131i_LNP = 0., Cquqd1_1132i_LNP = 0., Cquqd1_1133i_LNP = 0.,
            Cquqd1_1211i_LNP = 0., Cquqd1_1212i_LNP = 0., Cquqd1_1213i_LNP = 0., Cquqd1_1221i_LNP = 0., Cquqd1_1222i_LNP = 0., Cquqd1_1223i_LNP = 0., Cquqd1_1231i_LNP = 0., Cquqd1_1232i_LNP = 0., Cquqd1_1233i_LNP = 0.,
            Cquqd1_1311i_LNP = 0., Cquqd1_1312i_LNP = 0., Cquqd1_1313i_LNP = 0., Cquqd1_1321i_LNP = 0., Cquqd1_1322i_LNP = 0., Cquqd1_1323i_LNP = 0., Cquqd1_1331i_LNP = 0., Cquqd1_1332i_LNP = 0., Cquqd1_1333i_LNP = 0.,
            Cquqd1_2111i_LNP = 0., Cquqd1_2112i_LNP = 0., Cquqd1_2113i_LNP = 0., Cquqd1_2121i_LNP = 0., Cquqd1_2122i_LNP = 0., Cquqd1_2123i_LNP = 0., Cquqd1_2131i_LNP = 0., Cquqd1_2132i_LNP = 0., Cquqd1_2133i_LNP = 0.,
            Cquqd1_2211i_LNP = 0., Cquqd1_2212i_LNP = 0., Cquqd1_2213i_LNP = 0., Cquqd1_2221i_LNP = 0., Cquqd1_2222i_LNP = 0., Cquqd1_2223i_LNP = 0., Cquqd1_2231i_LNP = 0., Cquqd1_2232i_LNP = 0., Cquqd1_2233i_LNP = 0.,
            Cquqd1_2311i_LNP = 0., Cquqd1_2312i_LNP = 0., Cquqd1_2313i_LNP = 0., Cquqd1_2321i_LNP = 0., Cquqd1_2322i_LNP = 0., Cquqd1_2323i_LNP = 0., Cquqd1_2331i_LNP = 0., Cquqd1_2332i_LNP = 0., Cquqd1_2333i_LNP = 0.,
            Cquqd1_3111i_LNP = 0., Cquqd1_3112i_LNP = 0., Cquqd1_3113i_LNP = 0., Cquqd1_3121i_LNP = 0., Cquqd1_3122i_LNP = 0., Cquqd1_3123i_LNP = 0., Cquqd1_3131i_LNP = 0., Cquqd1_3132i_LNP = 0., Cquqd1_3133i_LNP = 0.,
            Cquqd1_3211i_LNP = 0., Cquqd1_3212i_LNP = 0., Cquqd1_3213i_LNP = 0., Cquqd1_3221i_LNP = 0., Cquqd1_3222i_LNP = 0., Cquqd1_3223i_LNP = 0., Cquqd1_3231i_LNP = 0., Cquqd1_3232i_LNP = 0., Cquqd1_3233i_LNP = 0.,
            Cquqd1_3311i_LNP = 0., Cquqd1_3312i_LNP = 0., Cquqd1_3313i_LNP = 0., Cquqd1_3321i_LNP = 0., Cquqd1_3322i_LNP = 0., Cquqd1_3323i_LNP = 0., Cquqd1_3331i_LNP = 0., Cquqd1_3332i_LNP = 0., Cquqd1_3333i_LNP = 0.; ///< The dimension-6 operator coefficient \f$(C_{quqd}^{(1)})_{ijkm}\f$ (Imaginary part).
    double Cquqd8_1111r_LNP = 0., Cquqd8_1112r_LNP = 0., Cquqd8_1113r_LNP = 0., Cquqd8_1121r_LNP = 0., Cquqd8_1122r_LNP = 0., Cquqd8_1123r_LNP = 0., Cquqd8_1131r_LNP = 0., Cquqd8_1132r_LNP = 0., Cquqd8_1133r_LNP = 0.,
            Cquqd8_1211r_LNP = 0., Cquqd8_1212r_LNP = 0., Cquqd8_1213r_LNP = 0., Cquqd8_1221r_LNP = 0., Cquqd8_1222r_LNP = 0., Cquqd8_1223r_LNP = 0., Cquqd8_1231r_LNP = 0., Cquqd8_1232r_LNP = 0., Cquqd8_1233r_LNP = 0.,
            Cquqd8_1311r_LNP = 0., Cquqd8_1312r_LNP = 0., Cquqd8_1313r_LNP = 0., Cquqd8_1321r_LNP = 0., Cquqd8_1322r_LNP = 0., Cquqd8_1323r_LNP = 0., Cquqd8_1331r_LNP = 0., Cquqd8_1332r_LNP = 0., Cquqd8_1333r_LNP = 0.,
            Cquqd8_2111r_LNP = 0., Cquqd8_2112r_LNP = 0., Cquqd8_2113r_LNP = 0., Cquqd8_2121r_LNP = 0., Cquqd8_2122r_LNP = 0., Cquqd8_2123r_LNP = 0., Cquqd8_2131r_LNP = 0., Cquqd8_2132r_LNP = 0., Cquqd8_2133r_LNP = 0.,
            Cquqd8_2211r_LNP = 0., Cquqd8_2212r_LNP = 0., Cquqd8_2213r_LNP = 0., Cquqd8_2221r_LNP = 0., Cquqd8_2222r_LNP = 0., Cquqd8_2223r_LNP = 0., Cquqd8_2231r_LNP = 0., Cquqd8_2232r_LNP = 0., Cquqd8_2233r_LNP = 0.,
            Cquqd8_2311r_LNP = 0., Cquqd8_2312r_LNP = 0., Cquqd8_2313r_LNP = 0., Cquqd8_2321r_LNP = 0., Cquqd8_2322r_LNP = 0., Cquqd8_2323r_LNP = 0., Cquqd8_2331r_LNP = 0., Cquqd8_2332r_LNP = 0., Cquqd8_2333r_LNP = 0.,
            Cquqd8_3111r_LNP = 0., Cquqd8_3112r_LNP = 0., Cquqd8_3113r_LNP = 0., Cquqd8_3121r_LNP = 0., Cquqd8_3122r_LNP = 0., Cquqd8_3123r_LNP = 0., Cquqd8_3131r_LNP = 0., Cquqd8_3132r_LNP = 0., Cquqd8_3133r_LNP = 0.,
            Cquqd8_3211r_LNP = 0., Cquqd8_3212r_LNP = 0., Cquqd8_3213r_LNP = 0., Cquqd8_3221r_LNP = 0., Cquqd8_3222r_LNP = 0., Cquqd8_3223r_LNP = 0., Cquqd8_3231r_LNP = 0., Cquqd8_3232r_LNP = 0., Cquqd8_3233r_LNP = 0.,
            Cquqd8_3311r_LNP = 0., Cquqd8_3312r_LNP = 0., Cquqd8_3313r_LNP = 0., Cquqd8_3321r_LNP = 0., Cquqd8_3322r_LNP = 0., Cquqd8_3323r_LNP = 0., Cquqd8_3331r_LNP = 0., Cquqd8_3332r_LNP = 0., Cquqd8_3333r_LNP = 0.; ///< The dimension-6 operator coefficient \f$(C_{quqd}^{(8)})_{ijkm}\f$ (Real part and pure real operator).
    double Cquqd8_1111i_LNP = 0., Cquqd8_1112i_LNP = 0., Cquqd8_1113i_LNP = 0., Cquqd8_1121i_LNP = 0., Cquqd8_1122i_LNP = 0., Cquqd8_1123i_LNP = 0., Cquqd8_1131i_LNP = 0., Cquqd8_1132i_LNP = 0., Cquqd8_1133i_LNP = 0.,
            Cquqd8_1211i_LNP = 0., Cquqd8_1212i_LNP = 0., Cquqd8_1213i_LNP = 0., Cquqd8_1221i_LNP = 0., Cquqd8_1222i_LNP = 0., Cquqd8_1223i_LNP = 0., Cquqd8_1231i_LNP = 0., Cquqd8_1232i_LNP = 0., Cquqd8_1233i_LNP = 0.,
            Cquqd8_1311i_LNP = 0., Cquqd8_1312i_LNP = 0., Cquqd8_1313i_LNP = 0., Cquqd8_1321i_LNP = 0., Cquqd8_1322i_LNP = 0., Cquqd8_1323i_LNP = 0., Cquqd8_1331i_LNP = 0., Cquqd8_1332i_LNP = 0., Cquqd8_1333i_LNP = 0.,
            Cquqd8_2111i_LNP = 0., Cquqd8_2112i_LNP = 0., Cquqd8_2113i_LNP = 0., Cquqd8_2121i_LNP = 0., Cquqd8_2122i_LNP = 0., Cquqd8_2123i_LNP = 0., Cquqd8_2131i_LNP = 0., Cquqd8_2132i_LNP = 0., Cquqd8_2133i_LNP = 0.,
            Cquqd8_2211i_LNP = 0., Cquqd8_2212i_LNP = 0., Cquqd8_2213i_LNP = 0., Cquqd8_2221i_LNP = 0., Cquqd8_2222i_LNP = 0., Cquqd8_2223i_LNP = 0., Cquqd8_2231i_LNP = 0., Cquqd8_2232i_LNP = 0., Cquqd8_2233i_LNP = 0.,
            Cquqd8_2311i_LNP = 0., Cquqd8_2312i_LNP = 0., Cquqd8_2313i_LNP = 0., Cquqd8_2321i_LNP = 0., Cquqd8_2322i_LNP = 0., Cquqd8_2323i_LNP = 0., Cquqd8_2331i_LNP = 0., Cquqd8_2332i_LNP = 0., Cquqd8_2333i_LNP = 0.,
            Cquqd8_3111i_LNP = 0., Cquqd8_3112i_LNP = 0., Cquqd8_3113i_LNP = 0., Cquqd8_3121i_LNP = 0., Cquqd8_3122i_LNP = 0., Cquqd8_3123i_LNP = 0., Cquqd8_3131i_LNP = 0., Cquqd8_3132i_LNP = 0., Cquqd8_3133i_LNP = 0.,
            Cquqd8_3211i_LNP = 0., Cquqd8_3212i_LNP = 0., Cquqd8_3213i_LNP = 0., Cquqd8_3221i_LNP = 0., Cquqd8_3222i_LNP = 0., Cquqd8_3223i_LNP = 0., Cquqd8_3231i_LNP = 0., Cquqd8_3232i_LNP = 0., Cquqd8_3233i_LNP = 0.,
            Cquqd8_3311i_LNP = 0., Cquqd8_3312i_LNP = 0., Cquqd8_3313i_LNP = 0., Cquqd8_3321i_LNP = 0., Cquqd8_3322i_LNP = 0., Cquqd8_3323i_LNP = 0., Cquqd8_3331i_LNP = 0., Cquqd8_3332i_LNP = 0., Cquqd8_3333i_LNP = 0.; ///< The dimension-6 operator coefficient \f$(C_{quqd}^{(8)})_{ijkm}\f$ (Imaginary part).
    double Clequ1_1111r_LNP = 0., Clequ1_1112r_LNP = 0., Clequ1_1113r_LNP = 0., Clequ1_1121r_LNP = 0., Clequ1_1122r_LNP = 0., Clequ1_1123r_LNP = 0., Clequ1_1131r_LNP = 0., Clequ1_1132r_LNP = 0., Clequ1_1133r_LNP = 0.,
            Clequ1_1211r_LNP = 0., Clequ1_1212r_LNP = 0., Clequ1_1213r_LNP = 0., Clequ1_1221r_LNP = 0., Clequ1_1222r_LNP = 0., Clequ1_1223r_LNP = 0., Clequ1_1231r_LNP = 0., Clequ1_1232r_LNP = 0., Clequ1_1233r_LNP = 0.,
            Clequ1_1311r_LNP = 0., Clequ1_1312r_LNP = 0., Clequ1_1313r_LNP = 0., Clequ1_1321r_LNP = 0., Clequ1_1322r_LNP = 0., Clequ1_1323r_LNP = 0., Clequ1_1331r_LNP = 0., Clequ1_1332r_LNP = 0., Clequ1_1333r_LNP = 0.,
            Clequ1_2111r_LNP = 0., Clequ1_2112r_LNP = 0., Clequ1_2113r_LNP = 0., Clequ1_2121r_LNP = 0., Clequ1_2122r_LNP = 0., Clequ1_2123r_LNP = 0., Clequ1_2131r_LNP = 0., Clequ1_2132r_LNP = 0., Clequ1_2133r_LNP = 0.,
            Clequ1_2211r_LNP = 0., Clequ1_2212r_LNP = 0., Clequ1_2213r_LNP = 0., Clequ1_2221r_LNP = 0., Clequ1_2222r_LNP = 0., Clequ1_2223r_LNP = 0., Clequ1_2231r_LNP = 0., Clequ1_2232r_LNP = 0., Clequ1_2233r_LNP = 0.,
            Clequ1_2311r_LNP = 0., Clequ1_2312r_LNP = 0., Clequ1_2313r_LNP = 0., Clequ1_2321r_LNP = 0., Clequ1_2322r_LNP = 0., Clequ1_2323r_LNP = 0., Clequ1_2331r_LNP = 0., Clequ1_2332r_LNP = 0., Clequ1_2333r_LNP = 0.,
            Clequ1_3111r_LNP = 0., Clequ1_3112r_LNP = 0., Clequ1_3113r_LNP = 0., Clequ1_3121r_LNP = 0., Clequ1_3122r_LNP = 0., Clequ1_3123r_LNP = 0., Clequ1_3131r_LNP = 0., Clequ1_3132r_LNP = 0., Clequ1_3133r_LNP = 0.,
            Clequ1_3211r_LNP = 0., Clequ1_3212r_LNP = 0., Clequ1_3213r_LNP = 0., Clequ1_3221r_LNP = 0., Clequ1_3222r_LNP = 0., Clequ1_3223r_LNP = 0., Clequ1_3231r_LNP = 0., Clequ1_3232r_LNP = 0., Clequ1_3233r_LNP = 0.,
            Clequ1_3311r_LNP = 0., Clequ1_3312r_LNP = 0., Clequ1_3313r_LNP = 0., Clequ1_3321r_LNP = 0., Clequ1_3322r_LNP = 0., Clequ1_3323r_LNP = 0., Clequ1_3331r_LNP = 0., Clequ1_3332r_LNP = 0., Clequ1_3333r_LNP = 0.; ///< The dimension-6 operator coefficient \f$(C_{lequ}^{(1)})_{ijkm}\f$ (Real part and pure real operator).
    double Clequ1_1111i_LNP = 0., Clequ1_1112i_LNP = 0., Clequ1_1113i_LNP = 0., Clequ1_1121i_LNP = 0., Clequ1_1122i_LNP = 0., Clequ1_1123i_LNP = 0., Clequ1_1131i_LNP = 0., Clequ1_1132i_LNP = 0., Clequ1_1133i_LNP = 0.,
            Clequ1_1211i_LNP = 0., Clequ1_1212i_LNP = 0., Clequ1_1213i_LNP = 0., Clequ1_1221i_LNP = 0., Clequ1_1222i_LNP = 0., Clequ1_1223i_LNP = 0., Clequ1_1231i_LNP = 0., Clequ1_1232i_LNP = 0., Clequ1_1233i_LNP = 0.,
            Clequ1_1311i_LNP = 0., Clequ1_1312i_LNP = 0., Clequ1_1313i_LNP = 0., Clequ1_1321i_LNP = 0., Clequ1_1322i_LNP = 0., Clequ1_1323i_LNP = 0., Clequ1_1331i_LNP = 0., Clequ1_1332i_LNP = 0., Clequ1_1333i_LNP = 0.,
            Clequ1_2111i_LNP = 0., Clequ1_2112i_LNP = 0., Clequ1_2113i_LNP = 0., Clequ1_2121i_LNP = 0., Clequ1_2122i_LNP = 0., Clequ1_2123i_LNP = 0., Clequ1_2131i_LNP = 0., Clequ1_2132i_LNP = 0., Clequ1_2133i_LNP = 0.,
            Clequ1_2211i_LNP = 0., Clequ1_2212i_LNP = 0., Clequ1_2213i_LNP = 0., Clequ1_2221i_LNP = 0., Clequ1_2222i_LNP = 0., Clequ1_2223i_LNP = 0., Clequ1_2231i_LNP = 0., Clequ1_2232i_LNP = 0., Clequ1_2233i_LNP = 0.,
            Clequ1_2311i_LNP = 0., Clequ1_2312i_LNP = 0., Clequ1_2313i_LNP = 0., Clequ1_2321i_LNP = 0., Clequ1_2322i_LNP = 0., Clequ1_2323i_LNP = 0., Clequ1_2331i_LNP = 0., Clequ1_2332i_LNP = 0., Clequ1_2333i_LNP = 0.,
            Clequ1_3111i_LNP = 0., Clequ1_3112i_LNP = 0., Clequ1_3113i_LNP = 0., Clequ1_3121i_LNP = 0., Clequ1_3122i_LNP = 0., Clequ1_3123i_LNP = 0., Clequ1_3131i_LNP = 0., Clequ1_3132i_LNP = 0., Clequ1_3133i_LNP = 0.,
            Clequ1_3211i_LNP = 0., Clequ1_3212i_LNP = 0., Clequ1_3213i_LNP = 0., Clequ1_3221i_LNP = 0., Clequ1_3222i_LNP = 0., Clequ1_3223i_LNP = 0., Clequ1_3231i_LNP = 0., Clequ1_3232i_LNP = 0., Clequ1_3233i_LNP = 0.,
            Clequ1_3311i_LNP = 0., Clequ1_3312i_LNP = 0., Clequ1_3313i_LNP = 0., Clequ1_3321i_LNP = 0., Clequ1_3322i_LNP = 0., Clequ1_3323i_LNP = 0., Clequ1_3331i_LNP = 0., Clequ1_3332i_LNP = 0., Clequ1_3333i_LNP = 0.; ///< The dimension-6 operator coefficient \f$(C_{lequ}^{(1)})_{ijkm}\f$ (Imaginary part).
    double Clequ3_1111r_LNP = 0., Clequ3_1112r_LNP = 0., Clequ3_1113r_LNP = 0., Clequ3_1121r_LNP = 0., Clequ3_1122r_LNP = 0., Clequ3_1123r_LNP = 0., Clequ3_1131r_LNP = 0., Clequ3_1132r_LNP = 0., Clequ3_1133r_LNP = 0.,
            Clequ3_1211r_LNP = 0., Clequ3_1212r_LNP = 0., Clequ3_1213r_LNP = 0., Clequ3_1221r_LNP = 0., Clequ3_1222r_LNP = 0., Clequ3_1223r_LNP = 0., Clequ3_1231r_LNP = 0., Clequ3_1232r_LNP = 0., Clequ3_1233r_LNP = 0.,
            Clequ3_1311r_LNP = 0., Clequ3_1312r_LNP = 0., Clequ3_1313r_LNP = 0., Clequ3_1321r_LNP = 0., Clequ3_1322r_LNP = 0., Clequ3_1323r_LNP = 0., Clequ3_1331r_LNP = 0., Clequ3_1332r_LNP = 0., Clequ3_1333r_LNP = 0.,
            Clequ3_2111r_LNP = 0., Clequ3_2112r_LNP = 0., Clequ3_2113r_LNP = 0., Clequ3_2121r_LNP = 0., Clequ3_2122r_LNP = 0., Clequ3_2123r_LNP = 0., Clequ3_2131r_LNP = 0., Clequ3_2132r_LNP = 0., Clequ3_2133r_LNP = 0.,
            Clequ3_2211r_LNP = 0., Clequ3_2212r_LNP = 0., Clequ3_2213r_LNP = 0., Clequ3_2221r_LNP = 0., Clequ3_2222r_LNP = 0., Clequ3_2223r_LNP = 0., Clequ3_2231r_LNP = 0., Clequ3_2232r_LNP = 0., Clequ3_2233r_LNP = 0.,
            Clequ3_2311r_LNP = 0., Clequ3_2312r_LNP = 0., Clequ3_2313r_LNP = 0., Clequ3_2321r_LNP = 0., Clequ3_2322r_LNP = 0., Clequ3_2323r_LNP = 0., Clequ3_2331r_LNP = 0., Clequ3_2332r_LNP = 0., Clequ3_2333r_LNP = 0.,
            Clequ3_3111r_LNP = 0., Clequ3_3112r_LNP = 0., Clequ3_3113r_LNP = 0., Clequ3_3121r_LNP = 0., Clequ3_3122r_LNP = 0., Clequ3_3123r_LNP = 0., Clequ3_3131r_LNP = 0., Clequ3_3132r_LNP = 0., Clequ3_3133r_LNP = 0.,
            Clequ3_3211r_LNP = 0., Clequ3_3212r_LNP = 0., Clequ3_3213r_LNP = 0., Clequ3_3221r_LNP = 0., Clequ3_3222r_LNP = 0., Clequ3_3223r_LNP = 0., Clequ3_3231r_LNP = 0., Clequ3_3232r_LNP = 0., Clequ3_3233r_LNP = 0.,
            Clequ3_3311r_LNP = 0., Clequ3_3312r_LNP = 0., Clequ3_3313r_LNP = 0., Clequ3_3321r_LNP = 0., Clequ3_3322r_LNP = 0., Clequ3_3323r_LNP = 0., Clequ3_3331r_LNP = 0., Clequ3_3332r_LNP = 0., Clequ3_3333r_LNP = 0.; ///< The dimension-6 operator coefficient \f$(C_{lequ}^{(3)})_{ijkm}\f$ (Real part and pure real operator).
    double Clequ3_1111i_LNP = 0., Clequ3_1112i_LNP = 0., Clequ3_1113i_LNP = 0., Clequ3_1121i_LNP = 0., Clequ3_1122i_LNP = 0., Clequ3_1123i_LNP = 0., Clequ3_1131i_LNP = 0., Clequ3_1132i_LNP = 0., Clequ3_1133i_LNP = 0.,
            Clequ3_1211i_LNP = 0., Clequ3_1212i_LNP = 0., Clequ3_1213i_LNP = 0., Clequ3_1221i_LNP = 0., Clequ3_1222i_LNP = 0., Clequ3_1223i_LNP = 0., Clequ3_1231i_LNP = 0., Clequ3_1232i_LNP = 0., Clequ3_1233i_LNP = 0.,
            Clequ3_1311i_LNP = 0., Clequ3_1312i_LNP = 0., Clequ3_1313i_LNP = 0., Clequ3_1321i_LNP = 0., Clequ3_1322i_LNP = 0., Clequ3_1323i_LNP = 0., Clequ3_1331i_LNP = 0., Clequ3_1332i_LNP = 0., Clequ3_1333i_LNP = 0.,
            Clequ3_2111i_LNP = 0., Clequ3_2112i_LNP = 0., Clequ3_2113i_LNP = 0., Clequ3_2121i_LNP = 0., Clequ3_2122i_LNP = 0., Clequ3_2123i_LNP = 0., Clequ3_2131i_LNP = 0., Clequ3_2132i_LNP = 0., Clequ3_2133i_LNP = 0.,
            Clequ3_2211i_LNP = 0., Clequ3_2212i_LNP = 0., Clequ3_2213i_LNP = 0., Clequ3_2221i_LNP = 0., Clequ3_2222i_LNP = 0., Clequ3_2223i_LNP = 0., Clequ3_2231i_LNP = 0., Clequ3_2232i_LNP = 0., Clequ3_2233i_LNP = 0.,
            Clequ3_2311i_LNP = 0., Clequ3_2312i_LNP = 0., Clequ3_2313i_LNP = 0., Clequ3_2321i_LNP = 0., Clequ3_2322i_LNP = 0., Clequ3_2323i_LNP = 0., Clequ3_2331i_LNP = 0., Clequ3_2332i_LNP = 0., Clequ3_2333i_LNP = 0.,
            Clequ3_3111i_LNP = 0., Clequ3_3112i_LNP = 0., Clequ3_3113i_LNP = 0., Clequ3_3121i_LNP = 0., Clequ3_3122i_LNP = 0., Clequ3_3123i_LNP = 0., Clequ3_3131i_LNP = 0., Clequ3_3132i_LNP = 0., Clequ3_3133i_LNP = 0.,
            Clequ3_3211i_LNP = 0., Clequ3_3212i_LNP = 0., Clequ3_3213i_LNP = 0., Clequ3_3221i_LNP = 0., Clequ3_3222i_LNP = 0., Clequ3_3223i_LNP = 0., Clequ3_3231i_LNP = 0., Clequ3_3232i_LNP = 0., Clequ3_3233i_LNP = 0.,
            Clequ3_3311i_LNP = 0., Clequ3_3312i_LNP = 0., Clequ3_3313i_LNP = 0., Clequ3_3321i_LNP = 0., Clequ3_3322i_LNP = 0., Clequ3_3323i_LNP = 0., Clequ3_3331i_LNP = 0., Clequ3_3332i_LNP = 0., Clequ3_3333i_LNP = 0.; ///< The dimension-6 operator coefficient \f$(C_{lequ}^{(3)})_{ijkm}\f$ (Imaginary part).
 
    double Lambda_NP; ///< The new physics scale [GeV].
    // The intrinsic and parametric theory relative errors in the Higgs observables.
    double eggFint = 0.; ///< Intrinsic relative theoretical error in ggF production. (Assumed to be constant in energy.)
    double eggFpar = 0.; ///< Parametric relative theoretical error in ggF production. (Assumed to be constant in energy.)
    double ettHint = 0.; ///< Intrinsic relative theoretical error in ttH production. (Assumed to be constant in energy.)
    double ettHpar = 0.; ///< Parametric relative theoretical error in ttH production. (Assumed to be constant in energy.)
    double eVBFint = 0.; ///< Intrinsic relative theoretical error in VBF production. (Assumed to be constant in energy.)
    double eVBFpar = 0.; ///< Parametric relative theoretical error in VBF production. (Assumed to be constant in energy.)
    double eWHint = 0.; ///< Intrinsic relative theoretical error in WH production. (Assumed to be constant in energy.)
    double eWHpar = 0.; ///< Parametric relative theoretical error in WH production. (Assumed to be constant in energy.)
    double eZHint = 0.; ///< Intrinsic relative theoretical error in ZH production. (Assumed to be constant in energy.)
    double eZHpar = 0.; ///< Parametric relative theoretical error in ZH production. (Assumed to be constant in energy.)
    double eeeWBFint = 0.; ///< Intrinsic relative theoretical error in \f$e^+ e^- \to H \bar{\nu} \nu\f$. (Assumed to be constant in energy.)
    double eeeWBFpar = 0.; ///< Parametric relative theoretical error in \f$e^+ e^- \to H \bar{\nu} \nu\f$. (Assumed to be constant in energy.)
    double eeeZHint = 0.; ///< Intrinsic relative theoretical error in \f$e^+ e^- \to Z H\f$. (Assumed to be constant in energy.)
    double eeeZHpar = 0.; ///< Parametric relative theoretical error in \f$e^+ e^- \to Z H\f$. (Assumed to be constant in energy.)
    double eeettHint = 0.; ///< Intrinsic relative theoretical error in \f$e^+ e^- \to t \bar{t} H\f$. (Assumed to be constant in energy.)
    double eeettHpar = 0.; ///< Parametric relative theoretical error in \f$e^+ e^- \to t \bar{t} H\f$. (Assumed to be constant in energy.)
    double eepWBFint = 0.; ///< Intrinsic relative theoretical error in \f$e^- p \to H e^- j\f$ via WBF. (Assumed to be constant in energy.)
    double eepWBFpar = 0.; ///< Parametric relative theoretical error in \f$e^- p \to H e^- j\f$ via WBF. (Assumed to be constant in energy.)    
    double eepZBFint = 0.; ///< Intrinsic relative theoretical error in \f$e^- p \to H e^- j\f$ via ZBF. (Assumed to be constant in energy.)
    double eepZBFpar = 0.; ///< Parametric relative theoretical error in \f$e^- p \to H e^- j\f$ via ZBF. (Assumed to be constant in energy.)
    double eHggint = 0.; ///< Intrinsic relative theoretical error in \f$H \to g g\f$.
    double eHggpar = 0.; ///< Parametric relative theoretical error in \f$H \to g g\f$.
    double eHWWint = 0.; ///< Intrinsic relative theoretical error in \f$H \to W W\f$.
    double eHWWpar = 0.; ///< Parametric relative theoretical error in \f$H \to W W\f$.
    double eHZZint = 0.; ///< Intrinsic relative theoretical error in \f$H \to Z Z\f$.
    double eHZZpar = 0.; ///< Parametric relative theoretical error in \f$H \to Z Z\f$.
    double eHZgaint = 0.; ///< Intrinsic relative theoretical error in \f$H \to Z \gamma\f$.
    double eHZgapar = 0.; ///< Parametric relative theoretical error in \f$H \to Z \gamma\f$.
    double eHgagaint = 0.; ///< Intrinsic relative theoretical error in \f$H \to \gamma\gamma\f$.
    double eHgagapar = 0.; ///< Parametric relative theoretical error in \f$H \to \gamma\gamma\f$.
    double eHmumuint = 0.; ///< Intrinsic relative theoretical error in \f$H \to \mu^+ \mu^-\f$.
    double eHmumupar = 0.; ///< Parametric relative theoretical error in \f$H \to \mu^+ \mu^-\f$.
    double eHtautauint = 0.; ///< Intrinsic relative theoretical error in \f$H \to \tau^+ \tau^-\f$.
    double eHtautaupar = 0.; ///< Parametric relative theoretical error in \f$H \to \tau^+ \tau^-\f$.
    double eHccint = 0.; ///< Intrinsic relative theoretical error in \f$H \to c\bar{c}\f$.
    double eHccpar = 0.; ///< Parametric relative theoretical error in \f$H \to c\bar{c}\f$.
    double eHssint = 0.; ///< Intrinsic relative theoretical error in \f$H \to s\bar{s}\f$.
    double eHsspar = 0.; ///< Parametric relative theoretical error in \f$H \to s\bar{s}\f$.
    double eHbbint = 0.; ///< Intrinsic relative theoretical error in \f$H \to b\bar{b}\f$.
    double eHbbpar = 0.; ///< Parametric relative theoretical error in \f$H \to b\bar{b}\f$.

    //  Intrinsic errors in \f$ee \to WW\f$ observables: error in total and differential cross section
    double eeeWWint = 0., edeeWWdcint = 0.; ///< Intrinsic relative theoretical error in \f$e^+ e^- \to W^+ W^-\f$: total cross section and \f$d/d\cos{\theta}\f$ distribution.

    double eggFHgaga = 0., eggFHZga = 0., eggFHZZ = 0., eggFHWW = 0., eggFHtautau = 0., eggFHbb = 0., eggFHmumu = 0.; ///< Total relative theoretical error in \f$gg \to H \to X\f$.   
    double eVBFHgaga = 0., eVBFHZga = 0., eVBFHZZ = 0., eVBFHWW = 0., eVBFHtautau = 0., eVBFHbb = 0., eVBFHmumu = 0.; ///< Total relative theoretical error in \f$pp \to Hjj (VBF) \to X jj\f$.   
    double eWHgaga = 0., eWHZga = 0., eWHZZ = 0., eWHWW = 0., eWHtautau = 0., eWHbb = 0., eWHmumu = 0.; ///< Total relative theoretical error in \f$pp \to WH \to W X\f$.    
    double eZHgaga = 0., eZHZga = 0., eZHZZ = 0., eZHWW = 0., eZHtautau = 0., eZHbb = 0., eZHmumu = 0.; ///< Total relative theoretical error in \f$pp \to ZH \to Z X\f$.
    double ettHgaga = 0., ettHZga = 0., ettHZZ = 0., ettHWW = 0., ettHtautau = 0., ettHbb = 0., ettHmumu = 0.; ///< Total relative theoretical error in \f$pp \to ttH \to tt X\f$.
    double eVBFHinv = 0., eVHinv = 0.; ///< Total relative theoretical error in \f$pp \to X H \to X + invisible\f$.

    double eHwidth = 0.; ///< Total relative theoretical error in the Higgs width

    double BrHinv = 0.; ///< The branching ratio of invisible Higgs decays.
    double BrHexo = 0.; ///< The branching ratio of exotic (not invisible) Higgs decays.

    //  ----------------------------------------------------------------------------
    //  Internal parameters    
    //  ----------------------------------------------------------------------------

    double LambdaNP2; ///< The square of the new physics scale [GeV\f$^2\f$].
    double LambdaNPm2; ///< The inverse square of the new physics scale [GeV\f$^{-2}\f$].

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

    double cNLOd6; ///< Parameter to control the inclusion of finite NLO corrections, when available.
    
    double cAsch, cWsch; ///< Parameters to control the SM EW input scheme: Alpha or MW.
    
//    double Yuke,Yukmu,Yuktau;///< SM lepton Yukawas
    double Yuku,Yukc,Yukt;///< SM u-quark Yukawas
//    double Yukd,Yuks,Yukb;///< SM d-quark Yukawas
    
    
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
    
    RGESolver SMEFTEvolEW;
    
    RGESolver SMEFTEvolMH, SMEFTEvol240, SMEFTEvol365, SMEFTEvol500;
    
    //These operators should be written in terms of the ones of the Warsaw basis -- TO BE FIXED!!!!!!!!
    double C2B = 0.; ///< The dimension-6 operator coefficient \f$C_{2W}(\Lambda_{\rm{EW}})\f$.
    double C2W = 0.; ///< The dimension-6 operator coefficient \f$C_{2B}(\Lambda_{\rm{EW}})\f$.
    double C2BS = 0.; ///< The dimension-6 operator coefficient \f$C_{2W}^{SILH}(\Lambda_{\rm{EW}})\f$.
    double C2WS = 0.; ///< The dimension-6 operator coefficient \f$C_{2B}^{SILH}(\Lambda_{\rm{EW}})\f$.
    double CDHB = 0.; ///< The dimension-6 operator coefficient \f$C_{DHB}(\Lambda_{\rm{EW}})\f$.
    double CDHW = 0.; ///< The dimension-6 operator coefficient \f$C_{DHW}(\Lambda_{\rm{EW}})\f$.
    double CDB = 0.; ///< The dimension-6 operator coefficient \f$C_{DB}(\Lambda_{\rm{EW}})\f$.
    double CDW = 0.; ///< The dimension-6 operator coefficient \f$C_{DW}(\Lambda_{\rm{EW}})\f$.
    double CT = 0.; ///< The dimension-6 operator coefficient \f$C_{T}(\Lambda_{\rm{EW}})\f$.
    
private:
    
    gslpp::matrix<gslpp::complex> VuL, VuLd, VuR, VuRd, VdL, VdLd, VdR, VdRd, VeL, VeLd, VeR, VeRd, MUQ, MDQ;

    bool FlagMWinput; ///< A boolean for the model flag %MWinput. 
    bool FlagQuadraticTerms; ///< A boolean flag that is true if the quadratic terms in cross sections and widths are switched on.
    bool FlagHiggsSM; ///< A boolean flag that is true if including dependence on small variations of the SM parameters (dependence is linearized). Available only in selected Higgs observables. 
    bool FlagLoopHd6; ///< A boolean flag that is true if including modifications in the SM loops in Higgs observables due to the dim 6 interactions.
    bool FlagLoopH3d6Quad; ///< A boolean flag that is true if including quadratic modifications in the SM loops in Higgs observables due to the dim 6 interactions that contribute to the trilinear Higgs coupling.
    bool FlagRGEci; ///< A boolean for the model flag %RGEci , to include RGE effects.
    bool FlagmultiScaleRGE; ///< A boolean for the model flag %multiScaleRGE , to include RGE effects at different scales.
    bool FlagfiniteNLO; ///< A boolean for the model flag %finiteNLO , to include finite NLO terms, when available. 
    bool FlagmatchLEFT; ///< A boolean for the model flag %matchLEFT , to include the matching between SMEFT and LEFT.

    gsl_integration_cquad_workspace * w_WW; /**< Gsl integral variable */
    double ZeroAle(double *dAle5h, double *params); //This was define as a const but we're changing the value of dAle5h so it shouldn't be a const function

};

#endif /* NPSMEFTD6GENERAL_H */

