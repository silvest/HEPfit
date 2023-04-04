/* 
 * Copyright (C) 2016 HEPfit Collaboration
 *
 *
 * For the licensing terms see doc/COPYING.
 */

#ifndef GENERALTHDM_H
#define	GENERALTHDM_H

#include "StandardModel.h"
#include "GeneralTHDMMatching.h"
#include "NPbase.h"



class GeneralTHDMcache; //forward reference to GeneralTHDMcache class

/**
 * @addtogroup GeneralTHDM
 * @brief A module for general Two-Higgs-Doublet models (without the assumption of a @f$Z_2@f$ symmetry).
 * @details This module contains code necessary for analyzing theoretical
 * and experimental constraints on general symmetric Two-Higgs-Doublet models.
 * While the theoretical bounds (positivity, unitarity and stability of the Higgs potential)
 * as well as the electroweak STU pseudo-observables
 * are available in the most general case,
 * the Higgs and flavour observables are currently only implemented in the flavour aligned limiting case.
 * @{
 * 
 * @class GeneralTHDM
 * @ingroup GeneralTHDM
 * @brief  A class for general Two-Higgs-Doublet models.
 * @author HEPfit Collaboration
 * @copyright GNU General Public License
 * @details 
 *
 * 
 * @anchor GeneralTHDMParameters
 * <h3>%Model parameters</h3>
 *
 * The model parameters of %GeneralTHDM are summarized below.
 * The current implementation allows for a Two-Higgs-Doublet model. The implementation allows for generañ
 * Yukawa matrices (even not allowed) and CP violation in the potenaital and in the Yukawa sector.
 * However some observables are only implemented in the aligned case (Aligned 2-Higgs doublet model).
 * The scalar 125 GeV resonance is assumed to be the one of the CP-even 
 * Higgs states of the model; the other is attributed to mH_2 in the GeneralTHDM configuration.
 * 
 * <table class="model">
 * <tr>
 *   <th>Label</th>
 *   <th>LaTeX symbol</th>
 *   <th>Description</th>
 * </tr>
 * <tr>
 *   <td class="mod_name">%logtb</td>
 *   <td class="mod_symb">@f$\log_{10} \tan \beta@f$</td>
 *   <td class="mod_desc">The decadic logarithm of the tangent of the mixing angle @f$\beta \equiv \arctan \frac{v_2}{v_1}@f$.</td>
 * </tr>
 * <tr>
 * <td class="mod_name">%mHp2</td>
 *   <td class="mod_symb">@f$m_{H^+}^2@f$</td>
 *   <td class="mod_desc">The squared mass of the charged Higgs.</td>
 * </tr>
 * <tr>
 *   <td class="mod_name">%mH2sq</td>
 *   <td class="mod_symb">@f$m_{H_2}^2@f$</td>
 *   <td class="mod_desc">The squared mass of the second Higgs @f$H_2@f$. Only used if flag use_sq_masses is set to true. </td>
 * </tr>
 * <tr>
 *   <td class="mod_name">%mH3sq</td>
 *   <td class="mod_symb">@f$m_{H_3}^2@f$</td>
 *   <td class="mod_desc">The squared mass of the third Higgs @f$H_3@f$. Only used if flag use_sq_masses is set to true. </td></td>
 * </tr>
 * <tr>
 * <td class="mod_name">%mH21</td>
 *   <td class="mod_symb">@f$m_{H_2}@f$</td>
 *   <td class="mod_desc">The masses of the second Higgs @f$H_2@f$. Only used if flag use_sq_masses is set to false. </td></td>
 * </tr>
 * <tr>
 *   <td class="mod_name">%mH31</td>
 *   <td class="mod_symb">@f$m_{H_3}@f$</td>
 *   <td class="mod_desc">The masses of the third Higgs @f$H_3@f$.Only used if flag use_sq_masses is set to false. </td></td>
 * </tr>
 * <tr>
 *  <td class="mod_name">%alpha1</td>
 *   <td class="mod_symb">@f$\alpha_1@f$</td>
 *   <td class="mod_desc">The  mixing angle @f$\alpha_1@f$.</td>
 * </tr>
 * <tr>
 * <td class="mod_name">%alpha2</td>
 *   <td class="mod_symb">@f$\alpha_2@f$</td>
 *   <td class="mod_desc">The  mixing angle @f$\alpha_2@f$.</td>
 * </tr>
 * <tr>
 *  <td class="mod_name">%alpha3</td>
 *   <td class="mod_symb">@f$\alpha_3@f$</td>
 *   <td class="mod_desc">The  mixing angle @f$\alpha_3@f$.</td>
 * </tr>
 * <tr>
 *   <td class="mod_name">%Relambda5</td>
 *   <td class="mod_symb">@f$Re(\lambda_{5})@f$</td>
 *   <td class="mod_desc">The real part of the potential paramter @f$\lambda_{5}@f$</td>
 * </tr>
 * <tr>
 * <td class="mod_name">%Imlambda5</td>
 *   <td class="mod_symb">@f$Im(\lambda_{5})@f$</td>
 *   <td class="mod_desc">The imaginary part of the potential paramter @f$\lambda_{5}@f$</td>
 * </tr>
 * <tr>
 * <td class="mod_name">%Relambda6</td>
 *   <td class="mod_symb">@f$Re(\lambda_{6})@f$</td>
 *   <td class="mod_desc">The real part of the potential paramter @f$\lambda_{6}@f$</td>
 * </tr>
 * <tr>
 *   <td class="mod_name">%Relambda7</td>
 *   <td class="mod_symb">@f$Re(\lambda_{7})@f$</td>
 *   <td class="mod_desc">The real part of the potential paramter @f$\lambda_{7}@f$</td>
 * </tr>
 * <tr>
 *   <td class="mod_name">%Q_GTHDM</td>
 *   <td class="mod_symb">@f$Q_{\text{GTHDM}}@f$</td>
 *   <td class="mod_desc">The GeneralTHDM scale.</td>
 * </tr>
 * <tr>
 *   <td class="mod_name">%RpepsGTHDM</td>
 *   <td class="mod_symb">@f$\epsilon_{\text{R'}}^{\text{GTHDM}}@f$</td>
 *   <td class="mod_desc">Minimal value for which the R' perturbativity criterion should be applied for the unitarity bounds.</td>
 * </tr>
 * <tr>
 *   <td class="mod_name">%NLOuniscaleGTHDM</td>
 *   <td class="mod_symb">@f$Q_{\text{min}}^{\text{NLOuni}}@f$</td>
 *   <td class="mod_desc">Minimal scale at which the NLO unitarity conditions are checked.</td>
 * </tr>
 * 
 * 
 * </tr>
 * <tr>
 *   <td class="mod_name">%Nu_11r</td>
 *   <td class="mod_symb">@f$Re(N_{u,11})@f$</td>
 *   <td class="mod_desc">Real part of the element (1,1) of the up-Yukawa matrix. Used only if use_sigma is set to true</td>
 * </tr>
 * <tr>
 *   <td class="mod_name">%Nu_11i</td>
 *   <td class="mod_symb">@f$Im(N_{u,11})@f$</td>
 *   <td class="mod_desc">Imaginary part of the element (1,1) of the up-Yukawa matrix. Used only if use_sigma is set to true</td>
 *  </tr>
 * <tr>
 *   <td class="mod_name">%Nu_12r</td>
 *   <td class="mod_symb">@f$Re(N_{u,12})@f$</td>
 *   <td class="mod_desc">Real part of the element (1,2) of the up-Yukawa matrix. Used only if use_sigma is set to true</td>
 * </tr>
 * <tr>
 *   <td class="mod_name">%Nu_12i</td>
 *   <td class="mod_symb">@f$Im(N_{u,12})@f$</td>
 *   <td class="mod_desc">Imaginary part of the element (1,2) of the up-Yukawa matrix. Used only if use_sigma is set to true</td>
 * </tr>
 * <tr>
 *   <td class="mod_name">%Nu_13r</td>
 *   <td class="mod_symb">@f$Re(N_{u,13})@f$</td>
 *   <td class="mod_desc">Real part of the element (1,3) of the up-Yukawa matrix. Used only if use_sigma is set to true</td>
 * </tr>
 * <tr>
 *   <td class="mod_name">%Nu_13i</td>
 *   <td class="mod_symb">@f$Im(N_{u,13})@f$</td>
 *   <td class="mod_desc">Imaginary part of the element (1,3) of the up-Yukawa matrix. Used only if use_sigma is set to true</td>
 * </tr>
 * <tr>
 *   <td class="mod_name">%Nu_21r</td>
 *   <td class="mod_symb">@f$Re(N_{u,21})@f$</td>
 *   <td class="mod_desc">Real part of the element (2,1) of the up-Yukawa matrix. Used only if use_sigma is set to true</td>
 * </tr>
 * <tr>
 *   <td class="mod_name">%Nu_21i</td>
 *   <td class="mod_symb">@f$Im(N_{u,21})@f$</td>
 *   <td class="mod_desc">Imaginary part of the element (2,1) of the up-Yukawa matrix. Used only if use_sigma is set to true</td>
 *  </tr>
 * <tr>
 *   <td class="mod_name">%Nu_22r</td>
 *   <td class="mod_symb">@f$Re(N_{u,22})@f$</td>
 *   <td class="mod_desc">Real part of the element (2,2) of the up-Yukawa matrix. Used only if use_sigma is set to true</td>
 * </tr>
 * <tr>
 *   <td class="mod_name">%Nu_22i</td>
 *   <td class="mod_symb">@f$Im(N_{u,22})@f$</td>
 *   <td class="mod_desc">Imaginary part of the element (2,2) of the up-Yukawa matrix. Used only if use_sigma is set to true</td>
 * </tr>
 * <tr>
 *   <td class="mod_name">%Nu_23r</td>
 *   <td class="mod_symb">@f$Re(N_{u,23})@f$</td>
 *   <td class="mod_desc">Real part of the element (2,3) of the up-Yukawa matrix. Used only if use_sigma is set to true</td>
 * </tr>
 * <tr>
 *   <td class="mod_name">%Nu_23i</td>
 *   <td class="mod_symb">@f$Im(N_{u,23})@f$</td>
 *   <td class="mod_desc">Imaginary part of the element (2,3) of the up-Yukawa matrix. Used only if use_sigma is set to true</td>
 *  </tr>
 * <tr>
 *   <td class="mod_name">%Nu_31r</td>
 *   <td class="mod_symb">@f$Re(N_{u,31})@f$</td>
 *   <td class="mod_desc">Real part of the element (3,1) of the up-Yukawa matrix. Used only if use_sigma is set to true</td>
 * </tr>
 * <tr>
 *   <td class="mod_name">%Nu_31i</td>
 *   <td class="mod_symb">@f$Im(N_{u,31})@f$</td>
 *   <td class="mod_desc">Imaginary part of the element (3,1) of the up-Yukawa matrix. Used only if use_sigma is set to true</td>
 *  </tr>
 * <tr>
 *   <td class="mod_name">%Nu_32r</td>
 *   <td class="mod_symb">@f$Re(N_{u,32})@f$</td>
 *   <td class="mod_desc">Real part of the element (3,2) of the up-Yukawa matrix. Used only if use_sigma is set to true</td>
 * </tr>
 * <tr>
 *   <td class="mod_name">%Nu_32i</td>
 *   <td class="mod_symb">@f$Im(N_{u,32})@f$</td>
 *   <td class="mod_desc">Imaginary part of the element (3,2) of the up-Yukawa matrix. Used only if use_sigma is set to true</td>
 * </tr>
 * <tr>
 *   <td class="mod_name">%Nu_33r</td>
 *   <td class="mod_symb">@f$Re(N_{u,33})@f$</td>
 *   <td class="mod_desc">Real part of the element (3,3) of the up-Yukawa matrix. Used only if use_sigma is set to true</td>
 * </tr>
 * <tr>
 *   <td class="mod_name">%Nu_33i</td>
 *   <td class="mod_symb">@f$Im(N_{u,33})@f$</td>
 *   <td class="mod_desc">Imaginary part of the element (3,3) of the up-Yukawa matrix. Used only if use_sigma is set to true</td>
 * </tr>
 * <tr>
 *   <td class="mod_name">%Nd_11r</td>
 *   <td class="mod_symb">@f$Re(N_{d,11})@f$</td>
 *   <td class="mod_desc">Real part of the element (1,1) of the down-Yukawa matrix. Used only if use_sigma is set to true</td>
 * </tr>
 * <tr>
 *   <td class="mod_name">%Nd_11i</td>
 *   <td class="mod_symb">@f$Im(N_{d,11})@f$</td>
 *   <td class="mod_desc">Imaginary part of the element (1,1) of the down-Yukawa matrix. Used only if use_sigma is set to true</td>
 *  </tr>
 * <tr>
 *   <td class="mod_name">%Nd_12r</td>
 *   <td class="mod_symb">@f$Re(N_{d,12})@f$</td>
 *   <td class="mod_desc">Real part of the element (1,2) of the down-Yukawa matrix. Used only if use_sigma is set to true</td>
 * </tr>
 * <tr>
 *   <td class="mod_name">%Nd_12i</td>
 *   <td class="mod_symb">@f$Im(N_{d,12})@f$</td>
 *   <td class="mod_desc">Imaginary part of the element (1,2) of the down-Yukawa matrix. Used only if use_sigma is set to true</td>
 * </tr>
 * <tr>
 *   <td class="mod_name">%Nd_13r</td>
 *   <td class="mod_symb">@f$Re(N_{d,13})@f$</td>
 *   <td class="mod_desc">Real part of the element (1,3) of the down-Yukawa matrix. Used only if use_sigma is set to true</td>
 * </tr>
 * <tr>
 *   <td class="mod_name">%Nd_13i</td>
 *   <td class="mod_symb">@f$Im(N_{d,13})@f$</td>
 *   <td class="mod_desc">Imaginary part of the element (1,3) of the down-Yukawa matrix. Used only if use_sigma is set to true</td>
 * </tr>
 * <tr>
 *   <td class="mod_name">%Nd_21r</td>
 *   <td class="mod_symb">@f$Re(N_{d,21})@f$</td>
 *   <td class="mod_desc">Real part of the element (2,1) of the down-Yukawa matrix. Used only if use_sigma is set to true</td>
 * </tr>
 * <tr>
 *   <td class="mod_name">%Nd_21i</td>
 *   <td class="mod_symb">@f$Im(N_{d,21})@f$</td>
 *   <td class="mod_desc">Imaginary part of the element (2,1) of the down-Yukawa matrix. Used only if use_sigma is set to true</td>
 *  </tr>
 * <tr>
 *   <td class="mod_name">%Nd_22r</td>
 *   <td class="mod_symb">@f$Re(N_{d,22})@f$</td>
 *   <td class="mod_desc">Real part of the element (2,2) of the down-Yukawa matrix. Used only if use_sigma is set to true</td>
 * </tr>
 * <tr>
 *   <td class="mod_name">%Nd_22i</td>
 *   <td class="mod_symb">@f$Im(N_{d,22})@f$</td>
 *   <td class="mod_desc">Imaginary part of the element (2,2) of the down-Yukawa matrix. Used only if use_sigma is set to true</td>
 * </tr>
 * <tr>
 *   <td class="mod_name">%Nd_23r</td>
 *   <td class="mod_symb">@f$Re(N_{d,23})@f$</td>
 *   <td class="mod_desc">Real part of the element (2,3) of the down-Yukawa matrix. Used only if use_sigma is set to true</td>
 * </tr>
 * <tr>
 *   <td class="mod_name">%Nd_23i</td>
 *   <td class="mod_symb">@f$Im(N_{d,23})@f$</td>
 *   <td class="mod_desc">Imaginary part of the element (2,3) of the down-Yukawa matrix. Used only if use_sigma is set to true</td>
 *  </tr>
 * <tr>
 *   <td class="mod_name">%Nd_31r</td>
 *   <td class="mod_symb">@f$Re(N_{d,31})@f$</td>
 *   <td class="mod_desc">Real part of the element (3,1) of the down-Yukawa matrix. Used only if use_sigma is set to true</td>
 * </tr>
 * <tr>
 *   <td class="mod_name">%Nd_31i</td>
 *   <td class="mod_symb">@f$Im(N_{d,31})@f$</td>
 *   <td class="mod_desc">Imaginary part of the element (3,1) of the down-Yukawa matrix. Used only if use_sigma is set to true</td>
 *  </tr>
 * <tr>
 *   <td class="mod_name">%Nd_32r</td>
 *   <td class="mod_symb">@f$Re(N_{d,32})@f$</td>
 *   <td class="mod_desc">Real part of the element (3,2) of the down-Yukawa matrix. Used only if use_sigma is set to true</td>
 * </tr>
 * <tr>
 *   <td class="mod_name">%Nd_32i</td>
 *   <td class="mod_symb">@f$Im(N_{d,32})@f$</td>
 *   <td class="mod_desc">Imaginary part of the element (3,2) of the down-Yukawa matrix. Used only if use_sigma is set to true</td>
 * </tr>
 * <tr>
 *   <td class="mod_name">%Nd_33r</td>
 *   <td class="mod_symb">@f$Re(N_{d,33})@f$</td>
 *   <td class="mod_desc">Real part of the element (3,3) of the down-Yukawa matrix. Used only if use_sigma is set to true</td>
 * </tr>
 * <tr>
 *   <td class="mod_name">%Nd_33i</td>
 *   <td class="mod_symb">@f$Im(N_{d,33})@f$</td>
 *   <td class="mod_desc">Imaginary part of the element (3,3) of the down-Yukawa matrix. Used only if use_sigma is set to true</td> 
 * </tr>
 * <tr>
 *   <td class="mod_name">%Nl_11r</td>
 *   <td class="mod_symb">@f$Re(N_{\ell,11})@f$</td>
 *   <td class="mod_desc">Real part of the element (1,1) of the leptonic Yukawa matrix. Used only if use_sigma is set to true</td>
 * </tr>
 * <tr>
 *   <td class="mod_name">%Nl_11i</td>
 *   <td class="mod_symb">@f$Im(N_{\ell,11})@f$</td>
 *   <td class="mod_desc">Imaginary part of the element (1,1) of the leptonic Yukawa matrix. Used only if use_sigma is set to true</td>
 *  </tr>
 * <tr>
 *   <td class="mod_name">%Nl_12r</td>
 *   <td class="mod_symb">@f$Re(N_{\ell,12})@f$</td>
 *   <td class="mod_desc">Real part of the element (1,2) of the leptonic Yukawa matrix. Used only if use_sigma is set to true</td>
 * </tr>
 * <tr>
 *   <td class="mod_name">%Nl_12i</td>
 *   <td class="mod_symb">@f$Im(N_{\ell,12})@f$</td>
 *   <td class="mod_desc">Imaginary part of the element (1,2) of the leptonic Yukawa matrix. Used only if use_sigma is set to true</td>
 * </tr>
 * <tr>
 *   <td class="mod_name">%Nl_13r</td>
 *   <td class="mod_symb">@f$Re(N_{\ell,13})@f$</td>
 *   <td class="mod_desc">Real part of the element (1,3) of the leptonic Yukawa matrix. Used only if use_sigma is set to true</td>
 * </tr>
 * <tr>
 *   <td class="mod_name">%Nl_13i</td>
 *   <td class="mod_symb">@f$Im(N_{\ell,13})@f$</td>
 *   <td class="mod_desc">Imaginary part of the element (1,3) of the leptonic Yukawa matrix. Used only if use_sigma is set to true</td>
 * </tr>
 * <tr>
 *   <td class="mod_name">%Nl_21r</td>
 *   <td class="mod_symb">@f$Re(N_{\ell,21})@f$</td>
 *   <td class="mod_desc">Real part of the element (2,1) of the leptonic Yukawa matrix. Used only if use_sigma is set to true</td>
 * </tr>
 * <tr>
 *   <td class="mod_name">%Nl_21i</td>
 *   <td class="mod_symb">@f$Im(N_{\ell,21})@f$</td>
 *   <td class="mod_desc">Imaginary part of the element (2,1) of the leptonic Yukawa matrix. Used only if use_sigma is set to true</td>
 *  </tr>
 * <tr>
 *   <td class="mod_name">%Nl_22r</td>
 *   <td class="mod_symb">@f$Re(N_{\ell,22})@f$</td>
 *   <td class="mod_desc">Real part of the element (2,2) of the leptonic Yukawa matrix. Used only if use_sigma is set to true</td>
 * </tr>
 * <tr>
 *   <td class="mod_name">%Nl_22i</td>
 *   <td class="mod_symb">@f$Im(N_{\ell,22})@f$</td>
 *   <td class="mod_desc">Imaginary part of the element (2,2) of the leptonic Yukawa matrix. Used only if use_sigma is set to true</td>
 * </tr>
 * <tr>
 *   <td class="mod_name">%Nl_23r</td>
 *   <td class="mod_symb">@f$Re(N_{\ell,23})@f$</td>
 *   <td class="mod_desc">Real part of the element (2,3) of the leptonic Yukawa matrix. Used only if use_sigma is set to true</td>
 * </tr>
 * <tr>
 *   <td class="mod_name">%Nl_23i</td>
 *   <td class="mod_symb">@f$Im(N_{\ell,23})@f$</td>
 *   <td class="mod_desc">Imaginary part of the element (2,3) of the leptonic Yukawa matrix. Used only if use_sigma is set to true</td>
 *  </tr>
 * <tr>
 *   <td class="mod_name">%Nl_31r</td>
 *   <td class="mod_symb">@f$Re(N_{\ell,31})@f$</td>
 *   <td class="mod_desc">Real part of the element (3,1) of the leptonic Yukawa matrix. Used only if use_sigma is set to true</td>
 * </tr>
 * <tr>
 *   <td class="mod_name">%Nl_31i</td>
 *   <td class="mod_symb">@f$Im(N_{\ell,31})@f$</td>
 *   <td class="mod_desc">Imaginary part of the element (3,1) of the leptonic Yukawa matrix. Used only if use_sigma is set to true</td>
 *  </tr>
 * <tr>
 *   <td class="mod_name">%Nl_32r</td>
 *   <td class="mod_symb">@f$Re(N_{\ell,32})@f$</td>
 *   <td class="mod_desc">Real part of the element (3,2) of the leptonic Yukawa matrix. Used only if use_sigma is set to true</td>
 * </tr>
 * <tr>
 *   <td class="mod_name">%Nl_32i</td>
 *   <td class="mod_symb">@f$Im(N_{\ell,32})@f$</td>
 *   <td class="mod_desc">Imaginary part of the element (3,2) of the leptonic Yukawa matrix. Used only if use_sigma is set to true</td>
 * </tr>
 * <tr>
 *   <td class="mod_name">%Nl_33r</td>
 *   <td class="mod_symb">@f$Re(N_{\ell,33})@f$</td>
 *   <td class="mod_desc">Real part of the element (3,3) of the leptonic Yukawa matrix. Used only if use_sigma is set to true</td>
 * </tr>
 * <tr>
 *   <td class="mod_name">%Nl_33i</td>
 *   <td class="mod_symb">@f$Im(N_{\ell,33})@f$</td>
 *   <td class="mod_desc">Imaginary part of the element (3,3) of the leptonic Yukawa matrix. Used only if use_sigma is set to true</td>
 * </tr>
 * <tr>
 *   <td class="mod_name">%yu1R_GTHDM</td>
 *   <td class="mod_symb">@f$y_{u}@f$</td>
 *   <td class="mod_desc">Yukawa coupling for up-quarks. Valid only in the aligned case and in the CP-conserved limit. @f$y_{u} = \cos \alpha_1 + \sin \alpha_1 Re(N_{u,11})@f$. Used only if use_sigma is set to false</td>
 * </tr>
 * <tr>
 *   <td class="mod_name">%yd1R_GTHDM</td>
 *   <td class="mod_symb">@f$y_{d}@f$</td>
 *   <td class="mod_desc">Yukawa coupling for down-quarks. Valid only in the aligned case and in the CP-conserved limit. @f$y_{d} = \cos \alpha_1 + \sin \alpha_1 Re(N{d,11})@f$. Used only if use_sigma is set to false</td>
 * </tr>
 * <tr>
 *   <td class="mod_name">%yl1R_GTHDM</td>
 *   <td class="mod_symb">@f$y_{\ell}@f$</td>
 *   <td class="mod_desc">Yukawa coupling for leptons. Valid only in the aligned case and in the CP-conserved limit.  @f$y_{\ell} = \cos \alpha_1 + \sin \alpha_1 Re(N_{\el,11}@f$=. Used only if use_sigma is set to false</td>
 * </table>
 * 
 * 
 * 
 *
 * @anchor GeneralTHDMFlags
 * <h3>%Model flags</h3>
 *
 * <table class="model">
 * <tr>
 *   <th>Label</th>
 *   <th>Value</th>
 *   <th>Description</th>
 * </tr>
 * <tr>
 *   <td class="mod_name">%ATHDMflag</td>
 *   <td class="mod_valu">true / false</td>
 *   <td class="mod_desc">This flag determines if the model is aligned (true) or not (false).</td>
 * </tr>
 * <tr>
 *  <td class="mod_name">%CPconservation</td>
 *   <td class="mod_valu">true / false</td>
 *   <td class="mod_desc">This flag determines if the model is CP-conserving (true) or not (false).</td>
 * </tr>
 * <tr>
 *   <td class="mod_name">%RGEorder</td>
 *   <td class="mod_valu">LO </td>
 *   <td class="mod_desc">This flag determines the order in perturbation theory of the renormalization group equations.</td>
 * </tr>
 * <tr>
 *   <td class="mod_name">%use_sq_masses</td>
 *   <td class="mod_valu">true / false</td>
 *   <td class="mod_desc">Whether to use the mass squared masses for the neutral scalars ot the mass.</td>
 * </tr>
 * <tr>
 *   <td class="mod_name">%use_sigma</td>
 *   <td class="mod_valu">true / false</td>
 *   <td class="mod_desc">Whether to use the elements  of the Yukawa couplings @f$N_{f,ij}@f$ or the couplings @f$y_f@f$. Valid only in the aligned case and in the CP-conserved limit</td>
 * </tr>
 *  * <tr>
 *   <td class="mod_name">%SM_Higgs</td>
 *   <td class="mod_valu">true / false</td>
 *   <td class="mod_desc">Whether the SM-Higgs of 125 GeV is chosen to be @f$\varphi_{1}@f$ (and the angle  @f$\tilde{\alpha}@f$) close to 0) or  the SM-Higgs of 125 GeV is chosen to be @f$\varphi_{2}@f$ (and the angle  @f$\tilde{\alpha}@f$ close to @f$\pi/2f$). So far only fully implemented at the aligend, CP-conserving limit </td>
 * </tr>
 * </table>
 */

class GeneralTHDM: public NPbase {
public:


    static const int NGeneralTHDMvars = 66;
    static std::string GeneralTHDMvars[NGeneralTHDMvars];

    /**
     * @brief GeneralTHDM constructor
     */
    GeneralTHDM();

    /**
     * @brief GeneralTHDM destructor
     */
    ~GeneralTHDM();
    
    /**
     * @brief A method to initialize the model.
     * @details This method, called via InputParser::ReadParameters(), allocates
     * memory to the pointers defined in the current class.
     * @return a boolean that is true if model initialization is successful
     */
    virtual bool InitializeModel();
    
    
    /**
     * @brief The post-update method for %GeneralTHDM.
     * @details This method runs all the procedures that are need to be executed
     * after the model is successfully updated.
     * @return a boolean that is true if the execution is successful
     */
    virtual bool PostUpdate();
    
    /**
     * @brief A method to check if all the mandatory parameters for %GeneralTHDM
     * have been provided in model initialization.
     * @param[in] DPars a map of the parameters that are being updated in the Monte Carlo run
     * (including parameters that are varied and those that are held constant)
     * @return a boolean that is true if the execution is successful
     */
    virtual bool CheckParameters(const std::map<std::string, double>& DPars);
    
     /**
     * @brief A get method to retrieve the matching scale @f$\mu_W@f$ around
     * the weak scale.
     * @return @f$\mu_W@f$ in GeV
     */
    double getMuw() const
    {
        return muw;
    }


    /**
     * @brief A get method to access the member reference of type GeneralTHDMMatching.
     * @return a reference to a GeneralTHDMMatching object
     */
    virtual GeneralTHDMMatching& getMatching() const
    {
        return GTHDMM.getObj();
    }


    ///////////////////////////////////////////////////////////////////////////
    // Flags

    /**
     * @brief A method to set a string flag of %GeneralTHDM.
     * @param[in] name name of a model flag
     * @param[in] value the string to be assigned to the flag specified by name
     * @return a boolean that is true if the execution is successful
     */
    virtual bool setFlagStr(const std::string name, const std::string value);
    
    /**
     * @brief A method to set a flag of %GeneralTHDM.
     * @param[in] name name of a model flag
     * @param[in] value the boolean to be assigned to the flag specified by name
     * @return a boolean that is true if the execution is successful
     */
    virtual bool setFlag(const std::string, const bool);

    GeneralTHDMcache* getMyGTHDMCache() const
    {
        return myGTHDMcache;
    }

    ///////////////////////////////////////////////////////////////////////////

    /**
     * @brief A getter for @f$\log_{10}(\tan \beta)@f$
     * @return @f$\log_{10}(\tan \beta)@f$
     */
    //double getlogtb() const {
    //    return logtb;
    //}

    /**
     * @brief A getter for @f$\tan \beta@f$
     * @return @f$\tan \beta@f$
     */
    //double gettanb() const {
    //    return tanb;
    //}

    /**
     * @brief A getter for @f$\sin \beta@f$
     * @return @f$\sin \beta@f$
     */
    //double getsinb() const {
    //    return sinb;
    //}

    /**
     * @brief A getter for @f$\cos \beta@f$
     * @return @f$\cos \beta@f$
     */
    //double getcosb() const {
    //    return cosb;
    //}

    /**
     *
     * @brief A getter for the squared charged Higgs mass
     * @return squared charged Higgs mass
     */
    double getmHp2() const {
        if(flag_use_sq_masses) {
            return mHp2;
        }
        else
        {
            return mHp1*mHp1;
        }
    }

     /**
     *
     * @brief A getter for the charged Higgs mass
     * @return charged Higgs mass
     */
    double getmHp() const {
        if(flag_use_sq_masses) {
            if(mHp2 < 0.) {
                throw std::runtime_error("error in GeneralTHDM: mHp2 < 0!");
            }
            else
            {
                return sqrt(mHp2);
            }
        }
        else
        {
                return mHp1;
        }
    }


    /**
     * @brief A getter for first neutral mass eigenvalue (The SM Higgs)
     * @return First neutral mass eigenvalue
     */
    double getmH1sq() const {
        if(flag_use_sq_masses) {
             if(mH1sq < 0.) {
                throw std::runtime_error("error in GeneralTHDM: mH1sq < 0!");
            }
                return mH1sq;
            }
         else
        {
            return mH1*mH1;
        }
    }
    
    
    
    /**
     * @brief A getter for second neutral mass eigenvalue
     * @return Second neutral mass eigenvalue
     */
    double getmH2sq() const {
        if(flag_use_sq_masses) {
             if(mH2sq < 0.) {
                throw std::runtime_error("error in GeneralTHDM: mH2sq < 0!");
            }
                return mH2sq;
            }
         else
        {
            return mH21*mH21;
        }
    }

      /**
     * @brief A getter for the third neutral mass eigenvalue
     * @return Third neutral mass eigenvalue
     */
      double getmH3sq() const {
        if(flag_use_sq_masses) {
             if(mH3sq < 0.) {
                throw std::runtime_error("error in GeneralTHDM: mH3sq < 0!");
            }
                return mH3sq;
            }
         else
        {
            return mH31*mH31;
        }
    }

    /**
     * @brief A getter for the rotation angle between the first and second neutral mass eigenstate
     * @return rotation angle between the first and second neutral mass eigenstate
     */
    double getalpha1() const {
        return alpha1;
    }

    /**
     * @brief A getter for cos of the rotation angle between the first and second neutral mass eigenstate
     * @return cos of the rotation angle between the first and second neutral mass eigenstate
     */
    double getcosalpha1() const {
        return cosalpha1;
    }

    /**
     * @brief A getter for sin of the rotation angle between the first and second neutral mass eigenstate
     * @return sin of the rotation angle between the first and second neutral mass eigenstate
     */
    double getsinalpha1() const {
        return sinalpha1;
    }
    
    /**
     * @brief A getter for tan of the rotation angle between the first and second neutral mass eigenstate
     * @return tan of the rotation angle between the first and second neutral mass eigenstate
     */
    double gettanalpha1() const {
        return tanalpha1;
    }
    
    
    /**
     * @brief A getter for the rotation angle between the first and third neutral mass eigenstate
     * @return rotation angle between the first and third neutral mass eigenstate
     */
    double getalpha2() const {
        if(flag_CPconservation==true) {
            return 0.0;
        }
        else {
            return alpha2;
        }
    }

    /**
     * @brief A getter for cos of the rotation angle between the first and third neutral mass eigenstate
     * @return cos of the rotation angle between the first and third neutral mass eigenstate
     */
    double getcosalpha2() const {
        if(flag_CPconservation==true) {
            return 1.0;
        }
        else {
            return cosalpha2;
        }
    }

    /**
     * @brief A getter for sin of the rotation angle between the first and third neutral mass eigenstate
     * @return sin of the rotation angle between the first and third neutral mass eigenstate
     */
    double getsinalpha2() const {
        if(flag_CPconservation==true) {
            return 0.0;
        }
        else {
            return sinalpha2;
        }
    }

    /**
     * @brief A getter for the rotation angle between the second and third neutral mass eigenstate
     * @return rotation angle between the second and third neutral mass eigenstate
     */
    double getalpha3() const {
        if(flag_CPconservation==true) {
            return 0.0;
        }
        else {
            return alpha3;
        }
    }

    /**
     * @brief A getter for cos of the rotation angle between the second and third neutral mass eigenstate
     * @return cos of the rotation angle between the second and third neutral mass eigenstate
     */
    double getcosalpha3() const {
        if(flag_CPconservation==true) {
            return 1.0;
        }
        else {
            return cosalpha3;
        }
    }

    /**
     * @brief A getter for sin of the rotation angle between the second and third neutral mass eigenstate
     * @return sin of the rotation angle between the second and third neutral mass eigenstate
     */
    double getsinalpha3() const {
        if(flag_CPconservation==true) {
            return 0.0;
        }
        else {
            return sinalpha3;
        }
    }
    
    
    /**
     * @brief A getter for the Higgs potential parameter @f$lambda_{2}@f$
     * @return Higgs potential parameter @f$lambda_{2}@f$
     */
    double getlambda1() const {
        return lambda1;
    }
    
    
    /**
     * @brief A getter for the Higgs potential parameter @f$lambda_{2}@f$
     * @return Higgs potential parameter @f$lambda_{2}@f$
     */
    double getlambda2() const {
        return lambda2;
    }
    
    
    /**
     * @brief A getter for the Higgs potential parameter @f$lambda_{3}@f$
     * @return Higgs potential parameter @f$lambda_{3}@f$
     */
    double getlambda3() const {
        return lambda3;
    }
    
    /**
     * @brief A getter for the Real part of the Higgs potential parameter @f$lambda_{5}@f$
     * @return Real part of the Higgs potential parameter @f$lambda_{5}@f$
     */
    double getlambda4() const {
        return lambda4;
    }
    
    /**
     * @brief A getter for the Real part of the Higgs potential parameter @f$lambda_{5}@f$
     * @return Real part of the Higgs potential parameter @f$lambda_{5}@f$
     */
    double getRelambda5() const {
        return Relambda5;
    }
    
    /**
     * @brief A getter for the Imaginary part of the Higgs potential parameter @f$lambda_{5}@f$
     * @return Imaginary part of the Higgs potential parameter @f$lambda_{5}@f$
     */
    double getImlambda5() const {
        if(flag_CPconservation==true) {
            return 0.0;
        }
        else {
            return Imlambda5;
        }
    }
    
    /**
     * @brief A getter for the Real part of the Higgs potential parameter @f$lambda_{6}@f$
     * @return Real part of the Higgs potential parameter @f$lambda_{6}@f$
     */
    double getRelambda6() const {
        return Relambda6;
    }

    
    /**
     * @brief A getter for the Imaginary part of the Higgs potential parameter @f$lambda_{5}@f$
     * @return Imaginary part of the Higgs potential parameter @f$lambda_{5}@f$
     */
    double getImlambda6() const {
        if(flag_CPconservation==true) {
            return 0.0;
        }
        else {
            return Imlambda6;
        }
    }
    
    /**
     * @brief A getter for the Real part of the Higgs potential parameter @f$lambda_{7}@f$
     * @return Real part of the Higgs potential parameter @f$lambda_{7}@f$
     */
    double getRelambda7() const {
        return Relambda7;
    }
    
    
    /**
     * @brief A getter for the Imaginary part of the Higgs potential parameter @f$lambda_{7}@f$
     * @return Imaginary part of the Higgs potential parameter @f$lambda_{7}@f$
     */
    double getImlambda7() const {
        if(flag_CPconservation==true) {
            return 0.0;
        }
        else {
            return Imlambda7;
        }
    }

    /**
     * 
     * @brief A getter for Nu matrix element 11
     * @return Nu matrix element 11
     */
    gslpp::complex getNu_11() const {
        gslpp::complex i = gslpp::complex::i();
         if(flag_sigma) {
        return (Nu_11r + i*Nu_11i);
        }
         else
         {
             if(alpha1==0)
               return 0 +i*0;
             else
                return (yu1R_GTHDM - cosalpha1)/(sinalpha1);
         }
    }
       
    /**
     * 
     * @brief A getter for Nu matrix element 12
     * @return Nu matrix element 12
     */
    gslpp::complex getNu_12() const {
        gslpp::complex i = gslpp::complex::i();
        return (Nu_12r + i*Nu_12i);
    }
    
    /**
     * 
     * @brief A getter for Nu matrix element 13
     * @return Nu matrix element 13
     */
    gslpp::complex getNu_13() const {
        gslpp::complex i = gslpp::complex::i();
        return (Nu_13r + i*Nu_13i);
    }
    
    /**
     * @brief A getter for Nu matrix element 21
     * @return Nu matrix element 21
     */
    gslpp::complex getNu_21() const {
        gslpp::complex i = gslpp::complex::i();
        return (Nu_21r + i*Nu_21i);
    }
    
    /**
     * @brief A getter for Nu matrix element 22
     * @return Nu matrix element 22
     */
    gslpp::complex getNu_22() const {
        gslpp::complex i = gslpp::complex::i();
        return (Nu_22r + i*Nu_22i);
    }
      
    /**
     * @brief A getter for Nu matrix element 23
     * @return Nu matrix element 23
     */
    gslpp::complex getNu_23() const {
        gslpp::complex i = gslpp::complex::i();
        return (Nu_23r + i*Nu_23i);
    }
    
    /**
     * @brief A getter for Nu matrix element 31
     * @return Nu matrix element 31
     */
    gslpp::complex getNu_31() const {
        gslpp::complex i = gslpp::complex::i();
        return (Nu_31r + i*Nu_31i);
    }
    
    /**
     * @brief A getter for Nu matrix element 32
     * @return Nu matrix element 32
     */
    gslpp::complex getNu_32() const {
        gslpp::complex i = gslpp::complex::i();
        return (Nu_32r + i*Nu_32i);
    }
    
    /**
     * @brief A getter for Nu matrix element 33
     * @return Nu matrix element 33
     */
    gslpp::complex getNu_33() const {
        gslpp::complex i = gslpp::complex::i();
        return (Nu_33r + i*Nu_33i);
    }
    
    /**
     * @brief A getter for Nd matrix element 11
     * @return Nd matrix element 11
     */
    gslpp::complex getNd_11() const {
        gslpp::complex i = gslpp::complex::i();
        if(flag_sigma) {
        return (Nd_11r + i*Nd_11i);
        }
     else
         {
            if(alpha1==0)
                return 0 +i*0;
            else
             return (yd1R_GTHDM - cosalpha1)/(sinalpha1);
         }
    }
    

    /**
     * @brief A getter for Nd matrix element 12
     * @return Nd matrix element 12
     */
    gslpp::complex getNd_12() const {
        gslpp::complex i = gslpp::complex::i();
        return (Nd_12r + i*Nd_12i);
    }
    
    /**
     * @brief A getter for Nd matrix element 13
     * @return Nd matrix element 13
     */
    gslpp::complex getNd_13() const {
        gslpp::complex i = gslpp::complex::i();
        return (Nd_13r + i*Nd_13i);
    }
    
    /**
     * @brief A getter for Nd matrix element 21
     * @return Nd matrix element 21
     */
    gslpp::complex getNd_21() const {
        gslpp::complex i = gslpp::complex::i();
        return (Nd_21r + i*Nd_21i);
    }
    
    /**
     * @brief A getter for Nd matrix element 22
     * @return Nd matrix element 22
     */
    gslpp::complex getNd_22() const {
        gslpp::complex i = gslpp::complex::i();
        return (Nd_22r + i*Nd_22i);
    }
    
    /**
     * @brief A getter for Nd matrix element 23
     * @return Nd matrix element 23
     */
    gslpp::complex getNd_23() const {
        gslpp::complex i = gslpp::complex::i();
        return (Nd_23r + i*Nd_23i);
    }
    
    /**
     * @brief A getter for Nd matrix element 31
     * @return Nd matrix element 31
     */
    gslpp::complex getNd_31() const {
        gslpp::complex i = gslpp::complex::i();
        return (Nd_31r + i*Nd_31i);
    }
    
    /**
     * @brief A getter for Nd matrix element 32
     * @return Nd matrix element 32
     */
    gslpp::complex getNd_32() const {
        gslpp::complex i = gslpp::complex::i();
        return (Nd_32r + i*Nd_32i);
    }
    
    /**
     * @brief A getter for Nd matrix element 33
     * @return Nd matrix element 33
     */
    gslpp::complex getNd_33() const {
        gslpp::complex i = gslpp::complex::i();
        return (Nd_33r + i*Nd_33i);
    }
    
    /**
     * @brief A getter for Nl matrix element 11
     * @return Nl matrix element 11
     */
    gslpp::complex getNl_11() const {
        gslpp::complex i = gslpp::complex::i();
        if(flag_sigma) {
        return (Nl_11r + i*Nl_11i);
        }
        else
         {
            if(alpha1==0)
                return 0 +i*0;
            else
             return (yl1R_GTHDM  -cosalpha1)/(sinalpha1);
         }
    }
    
    /**
     * @brief A getter for Nl matrix element 12
     * @return Nl matrix element 12
     */
    gslpp::complex getNl_12() const {
        gslpp::complex i = gslpp::complex::i();
        return (Nl_12r + i*Nl_12i);
    }
    
    /**
     * @brief A getter for Nl matrix element 13
     * @return Nl matrix element 13
     */
    gslpp::complex getNl_13() const {
        gslpp::complex i = gslpp::complex::i();
        return (Nl_13r + i*Nl_13i);
    }
    
    /**
     * @brief A getter for Nl matrix element 21
     * @return Nl matrix element 21
     */
    gslpp::complex getNl_21() const {
        gslpp::complex i = gslpp::complex::i();
        return (Nl_21r + i*Nl_21i);
    }
    
    /**
     * @brief A getter for Nl matrix element 22
     * @return Nl matrix element 22
     */
    gslpp::complex getNl_22() const {
        gslpp::complex i = gslpp::complex::i();
        return (Nl_22r + i*Nl_22i);
    }
    
    /**
     * @brief A getter for Nl matrix element 23
     * @return Nl matrix element 23
     */
    gslpp::complex getNl_23() const {
        gslpp::complex i = gslpp::complex::i();
        return (Nl_23r + i*Nl_23i);
    }
    
    /**
     * @brief A getter for Nl matrix element 31
     * @return Nl matrix element 31
     */
    gslpp::complex getNl_31() const {
        gslpp::complex i = gslpp::complex::i();
        return (Nl_31r + i*Nl_31i);
    }
    
    /**
     * @brief A getter for Nl matrix element 32
     * @return Nl matrix element 32
     */
    gslpp::complex getNl_32() const {
        gslpp::complex i = gslpp::complex::i();
        return (Nl_32r + i*Nl_32i);
    }
    
    /**
     * @brief A getter for Nl matrix element 33
     * @return Nl matrix element 33
     */
    gslpp::complex getNl_33() const {
        gslpp::complex i = gslpp::complex::i();
        return (Nl_33r + i*Nl_33i);
    }
    
    
    
    gslpp::complex getyu1() const {
        if(!flag_sigma) {
        return (yu1R_GTHDM);
        }
         else
         {
            gslpp::complex i = gslpp::complex::i();
            if(flag_SM_Higgs){
            double R11 = cosalpha1*cosalpha2;
            double R12 = sinalpha1*cosalpha2;
            double R13 = -sinalpha2;
            gslpp::complex yu1 = (R11 + (R12 - i*R13)*(Nu_11r + i*Nu_11i).conjugate());
           // if(yu1.abs2() > 4.) return std::numeric_limits<double>::quiet_NaN();
            return yu1;
            }
            else{
                double R21 = cosalpha1*sinalpha2*sinalpha3 - sinalpha1*cosalpha3;
                double R22 = cosalpha1*cosalpha3+sinalpha1*sinalpha2*sinalpha3;
                double R23 = cosalpha2*sinalpha3;
                gslpp::complex yu1 = (R21 + (R22 - i*R23)*(Nu_11r + i*Nu_11i).conjugate());
               // if(yu1.abs2() > 4.) return std::numeric_limits<double>::quiet_NaN();
                return yu1;
            }
         }
    }
        
     
    double getyu1R() const {
        if(!flag_sigma) {
        return (yu1R_GTHDM);
        }
         else
            {
             if(flag_SM_Higgs){
            double R11 = cosalpha1*cosalpha2;
            double R12 = sinalpha1*cosalpha2;
            double yu1R = (R11 + R12*Nu_11r);
          //  if(yu1R*yu1R > 4.) return std::numeric_limits<double>::quiet_NaN();
            return yu1R;
             }
            else{
                double R21 = - sinalpha1;
                double R22 = cosalpha1;
                double yu1R =(R21 + R22*Nu_11r);
            //    if(yu1R*yu1R > 4.) return std::numeric_limits<double>::quiet_NaN();
                return yu1R;
            }
       }
    }    
        
        
    gslpp::complex getyd1() const {
        if(!flag_sigma) {
        return (yd1R_GTHDM);
        }
         else
         {
            gslpp::complex i = gslpp::complex::i();
            if(flag_SM_Higgs){
            double R11 = cosalpha1*cosalpha2;
            double R12 = sinalpha1*cosalpha2;
            double R13 = -sinalpha2;
            gslpp::complex yd1 = (R11 + (R12 - i*R13)*(Nd_11r + i*Nd_11i));
          //  if(yd1.abs2() > 4.) return std::numeric_limits<double>::quiet_NaN();
            return yd1;
            }
            else{
                double R21 = cosalpha1*sinalpha2*sinalpha3 - sinalpha1*cosalpha3;
                double R22 = cosalpha1*cosalpha3+sinalpha1*sinalpha2*sinalpha3;
                double R23 = cosalpha2*sinalpha3;   
                gslpp::complex yd1 = (R21 + (R22 - i*R23)*(Nd_11r + i*Nd_11i));
           //     if(yd1.abs2() > 4.) return std::numeric_limits<double>::quiet_NaN();
                 return yd1;
            }
            }
         }
    
    
    double getyd1R() const {
        if(!flag_sigma) {
        return (yd1R_GTHDM);
        }
         else
         {
            if(flag_SM_Higgs){
            double R11 = cosalpha1*cosalpha2;
            double R12 = sinalpha1*cosalpha2;
            double yd1R = (R11 + R12*Nd_11r);
           // if(yd1R*yd1R > 4.) return std::numeric_limits<double>::quiet_NaN();
            return yd1R;
                 }
            else{
                double R21 = - sinalpha1;
                double R22 = cosalpha1;
                double yd1R = (R21 + R22*Nd_11r);
            //    if(yd1R*yd1R > 4.) return std::numeric_limits<double>::quiet_NaN();
                return yd1R;
                 }
         }
    }
    
        gslpp::complex getyl1() const {
         if(!flag_sigma) {
        return (yl1R_GTHDM);
        }
         else
         {
            gslpp::complex i = gslpp::complex::i();
            if(flag_SM_Higgs){
             double R11 = cosalpha1*cosalpha2;
             double R12 = sinalpha1*cosalpha2;
             double R13 = -sinalpha2;
             gslpp::complex yl1 = (R11 + (R12 - i*R13)*(Nl_11r + i*Nl_11i)); 
           // if(yl1.abs2() > 4.) return std::numeric_limits<double>::quiet_NaN();
            return yl1;
             
            }
            else{
                double R21 = cosalpha1*sinalpha2*sinalpha3 - sinalpha1*cosalpha3;
                double R22 = cosalpha1*cosalpha3+sinalpha1*sinalpha2*sinalpha3;
                double R23 = cosalpha2*sinalpha3;
                gslpp::complex yl1 = (R21 + (R22 - i*R23)*(Nl_11r + i*Nl_11i)); 
            //     if(yl1.abs2() > 4.) return std::numeric_limits<double>::quiet_NaN();
                 return yl1;
            }
         }
    }
    
    
    double getyl1R() const {
         if(!flag_sigma) {
        return (yl1R_GTHDM);
        }
         else
         {
           if(flag_SM_Higgs){
            double R11 = cosalpha1*cosalpha2;
            double R12 = sinalpha1*cosalpha2;            
            double yl1R = (R11 + R12*Nl_11r);
            //if(yl1R*yl1R > 4.) return std::numeric_limits<double>::quiet_NaN();
            return yl1R;
            
           }
           else{
             double R21 = - sinalpha1;
             double R22 = cosalpha1;   
             double yl1R = (R21 + R22*Nl_11r);
            //if(yl1R*yl1R > 4.) return std::numeric_limits<double>::quiet_NaN();
            return yl1R;  
           }
         }
    }
    
    /**
     * @brief A getter for the General THDM scale
     * @return General THDM scale
     */
    double getQ_GTHDM() const {
        return Q_GTHDM;
    }

    /**
     * @brief A getter for the minimal R' value
     * @return Minimal R' value
     */
    double getRpepsGTHDM() const {
        return RpepsGTHDM;
    }

    /**
     * @brief A getter for the minimal NLO unitarity check scale
     * @return Minimal NLO unitarity check scale
     */
    double getNLOuniscaleGTHDM() const {
        return NLOuniscaleGTHDM;
    }

    /**
     *
     * @brief A getter for the flag to switch on the Aligned Two Higgs Doublet Model
     * @return Flag to switch on the Aligned Two Higgs Doublet Model
     */
    bool getATHDMflag() const {
        return flag_ATHDM;
    }

    /**
     *
     * @brief A getter for the flag to switch on CP conservation
     * @return Flag to switch on CP conservation
     */
    bool getCPconservationflag() const {
        return flag_CPconservation;
    }

    /**
     *
     * @brief A switch for NLO RGE and approximate NLO RGE
     * @return Switch for NLO RGE and approximate NLO RGE
     */
    std::string getRGEorderflag() const {
        return flag_RGEorder;
    }
    
     /**
     *
      * @brief A getter for the choice if you want to use the THDM masses or 
      * rather their squares
     * @return Choose if you want to use the THDM masses or rather their squares
     */
    bool getsqmassesflag() const {
        return flag_use_sq_masses;
    }
/**
     *
     * @brief A getter for the flag to choose which scalar is the SM one.
 *   * Default value (true) is the SM-Higgs is the first one. Only fully impleneted in the aligned CP-conserving limit
     * @return Flag to select the SM Higgs of 125 GeV
     */
    bool getSMHiggs() const {
        return flag_SM_Higgs;
    }

    virtual double muggH(const double sqrt_s) const;
    virtual double muVBF(const double sqrt_s) const;
    virtual double mueeWBF(const double sqrt_s) const;
    virtual double muWH(const double sqrt_s) const;
    virtual double muZH(const double sqrt_s) const;
    virtual double mueeZH(const double sqrt_s) const;
    virtual double muVH(const double sqrt_s) const;
    virtual double muVBFpVH(const double sqrt_s) const;
    virtual double muttH(const double sqrt_s) const;
    virtual double GammaTotal() const;
    virtual double BrHggRatio() const;
    virtual double BrHWWRatio() const;
    virtual double BrHZZRatio() const;
    virtual double BrHZgaRatio() const;
    virtual double BrHgagaRatio() const;
    virtual double BrHmumuRatio() const;
    virtual double BrHtautauRatio() const;
    virtual double BrHccRatio() const;
    virtual double BrHbbRatio() const;
    virtual double muggHgaga(const double sqrt_s) const;
    virtual double muVBFHgaga(const double sqrt_s) const;
    virtual double muVHgaga(const double sqrt_s) const;
    virtual double muttHgaga(const double sqrt_s) const;
    virtual double muggHZZ(const double sqrt_s) const;
    virtual double muVBFHZZ(const double sqrt_s) const;
    virtual double muVHZZ(const double sqrt_s) const;
    virtual double muttHZZ(const double sqrt_s) const;
    virtual double muggHWW(const double sqrt_s) const;
    virtual double muVBFHWW(const double sqrt_s) const;
    virtual double muVHWW(const double sqrt_s) const;
    virtual double muttHWW(const double sqrt_s) const;
    virtual double muggHtautau(const double sqrt_s) const;
    virtual double muVBFHtautau(const double sqrt_s) const;
    virtual double muVHtautau(const double sqrt_s) const;
    virtual double muttHtautau(const double sqrt_s) const;
    virtual double muggHbb(const double sqrt_s) const;
    virtual double muVBFHbb(const double sqrt_s) const;
    virtual double muVHbb(const double sqrt_s) const;
    virtual double muttHbb(const double sqrt_s) const;
    virtual double muppHmumu(const double sqrt_s) const;
    virtual double muppHZga(const double sqrt_s) const;
    virtual double computeGammaTotalRatio() const;

    
    /**
     * @brief Function @f$F(m02,m12)@f$ used for %THDM. Remember that this function is
     * defined for %GeneralTHDM while for SUSY we have a multiplicative factor 2.
     * @param[in] m02 mass square @f$m_0^2@f$
     * @param[in] m12 mass square @f$m_1^2@f$
     * @return @f$F(m02,m12)@f$
     */
    double F(const double m02, const double m12) const;
    
    
    
    /**
     * @brief A method that returns the New Physics contribution to the oblique parameter S in the GeneralTHDM
     * @return a double with the New Physics contribution to the oblique parameter S.
     */
    double GTHDMDeltaS() const;
    
    /**
     * @brief A method that returns the New Physics contribution to the oblique parameter T in the GeneralTHDM
     * @return a double with the New Physics contribution to the oblique parameter T.
     */
    double GTHDMDeltaT() const;
    
    
    /**
     * @brief A method that returns the New Physics contribution to the oblique parameter U in the GeneralTHDM
     * @return a double with the New Physics contribution to the oblique parameter U.
     */
    double GTHDMDeltaU() const;
    
    
    //  Why was this done??? It doesn't seem to make any sense
    /* 
     * @brief Method overriding the NPbase Mw with the tree-level W mass
     */
    //virtual double Mw() const;
    
    //We substitute the tree-level value with the real SM value but the oblique corrections are missing
    /* 
     * @brief Method overriding the NPbase Mw with the SM W mass. We
     * did not include the oblique corrections for the moment but it 
     * should be done for consistency
     */
    virtual double Mw() const;

protected:

    /**
     * @brief A method to set the value of a parameter of %GeneralTHDM.
     * @param[in] name name of a model parameter
     * @param[in] value the value to be assigned to the parameter specified by name
     */
    virtual void setParameter(const std::string, const double&);

    /**
     * @brief A method to check if the RGE order name in string form is valid.
     * @param[in] GeneralTHDM RGE order
     * @return a boolean that is true if the RGE order string is valid
     */
    bool checkRGEorder(const std::string RGEorder) const
    {
        if (RGEorder.compare("LO") == 0)
//                || RGEorder.compare("approxNLO") == 0)
            return true;
        else
            return false;
    }

    mutable Matching<GeneralTHDMMatching,GeneralTHDM> GTHDMM; ///< An object of type Matching.

private:


    
    GeneralTHDMcache* myGTHDMcache;

    //Let's define here all the parameters of the model, including those who are linearly dependent
    //We'll also set all the parameters in the setParameter functions in such a way that we can use the formulae written in terms of any parameter
    //These parameters are private but used in the inherited classes, think a bit about if this is the best way of including them
    double /*logtb, tanb, sinb, cosb,*/mH1sq, mHp2, mH2sq, mH3sq, mH1, mHp1, mH21, mH31, alpha1, cosalpha1, sinalpha1, tanalpha1, alpha2, cosalpha2, sinalpha2,
            alpha3, cosalpha3, sinalpha3, mu2, lambda1, lambda2, lambda3, lambda4, Relambda5, Imlambda5, Relambda6, Imlambda6, Relambda7, Imlambda7,//we can trade one of the Im by the Im of lambda7, we include it any way since we've not decided yet which phase to keep
            yu1R_GTHDM, yd1R_GTHDM, yl1R_GTHDM, 
            Nu_11r, Nu_11i, Nu_12r, Nu_12i, Nu_13r, Nu_13i, 
            Nu_21r, Nu_21i, Nu_22r, Nu_22i, Nu_23r, Nu_23i, 
            Nu_31r, Nu_31i, Nu_32r, Nu_32i, Nu_33r, Nu_33i, 
            Nd_11r, Nd_11i, Nd_12r, Nd_12i, Nd_13r, Nd_13i, 
            Nd_21r, Nd_21i, Nd_22r, Nd_22i, Nd_23r, Nd_23i, 
            Nd_31r, Nd_31i, Nd_32r, Nd_32i, Nd_33r, Nd_33i, 
            Nl_11r, Nl_11i, Nl_12r, Nl_12i, Nl_13r, Nl_13i, 
            Nl_21r, Nl_21i, Nl_22r, Nl_22i, Nl_23r, Nl_23i, 
            Nl_31r, Nl_31i, Nl_32r, Nl_32i, Nl_33r, Nl_33i, 
            Q_GTHDM, RpepsGTHDM, NLOuniscaleGTHDM;
    bool flag_ATHDM, flag_CPconservation, flag_use_sq_masses, flag_sigma, flag_SM_Higgs;
    std::string flag_RGEorder;
};

/**
 * @}
 */

#endif	/* GENERALTHDM_H */
