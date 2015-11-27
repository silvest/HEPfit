/*
 * Copyright (C) 2012 HEPfit Collaboration
 * All rights reserved.
 *
 * For the licensing terms see doc/COPYING.
 */

#ifndef QCD_H
#define	QCD_H

#include "Model.h"
#include "Meson.h"
#include "OrderScheme.h"
#define MEPS 1.e-10 // mass precision

/**
 * @class QCD
 * @ingroup StandardModel
 * @brief A class for parameters related to %QCD, hadrons and quarks.
 * @author HEPfit Collaboration
 * @copyright GNU General Public License
 * @details This class is a Model class that assigns and updates parameters
 * related to and derived from %QCD. A complete list of parameters in the %QCD
 * class can be found below. This class includes, but is not limited to,
 * the running of the strong coupling constant (Full NNLO), running of the quark
 * masses and conversions between pole mass and \f$\overline{\mathrm{MS}}\f$ mass. All hadronization
 * parameters like the bag parameters for the mesons and their decay constants are
 * assigned and updated by this class.
 *
 *
 * @anchor QCDInitialization
 * <h3>Initialization</h3>
 *
 * The constructor QCD() sets the charge and isospin of the quarks. It also sets the
 * mass scale of the light quarks UP, DOWN and STRANGE to 2 GeV. The cache is initialized
 * too along with the computation of \f$\zeta(2)\f$ and \f$\zeta(3)\f$.
 *
 * The initializations and updates of the model parameters and flags are explained
 * below.
 *
 *
 * @anchor QCDParameters
 * <h3>%Model parameters</h3>
 *
 * The model parameters of %QCD are summarized below:
 * <table class="model">
 * <tr>
 *   <th>Label</th>
 *   <th>LaTeX symbol</th>
 *   <th>Description</th>
 * </tr>
 * <tr>
 *   <td class="mod_name">%AlsM</td>
 *   <td class="mod_symb">@f$\alpha_s(M_{\alpha_s})@f$</td>
 *   <td class="mod_desc">The strong coupling constant at the scale \f$M_{\alpha_s}\f$.</td>
 * </tr>
 * <tr>
 *   <td class="mod_name">%MAls</td>
 *   <td class="mod_symb">@f$M_{\alpha_s}@f$</td>
 *   <td class="mod_desc">The mass scale in GeV at which the strong coupling constant measurement is provided. Must be in the 5 flavour energy range.</td>
 * </tr>
 * <tr>
 *   <td class="mod_name">%mup</td>
 *   <td class="mod_symb">@f$m_{u}(2\,\mathrm{GeV})@f$</td>
 *   <td class="mod_desc">The \f$\overline{\mathrm{MS}}\f$ mass of the up quark at 2 GeV, \f$m_u(2\,\mathrm{GeV})\f$, in GeV.</td>
 * </tr>
 * <tr>
 *   <td class="mod_name">%mdown</td>
 *   <td class="mod_symb">@f$m_{d}(2\,\mathrm{GeV})@f$</td>
 *   <td class="mod_desc">The \f$\overline{\mathrm{MS}}\f$ mass of the down quark at 2 GeV, \f$m_d(2\,\mathrm{GeV})\f$, in GeV.</td>
 * </tr>
 * <tr>
 *   <td class="mod_name">%mcharm</td>
 *   <td class="mod_symb">@f$m_{c}(m_c)@f$</td>
 *   <td class="mod_desc">The \f$\overline{\mathrm{MS}}\f$ scale-invariant mass of the charm quark, \f$m_c(m_c)\f$, in GeV.</td>
 * </tr>
 * <tr>
 *   <td class="mod_name">%mstrange</td>
 *   <td class="mod_symb">@f$m_{s}(2\,\mathrm{GeV})@f$</td>
 *   <td class="mod_desc">The \f$\overline{\mathrm{MS}}\f$ mass of the strange quark at 2 GeV , \f$m_s(2\,\mathrm{GeV})\f$, in GeV.</td>
 * </tr>
 * <tr>
 *   <td class="mod_name">%mtop</td>
 *   <td class="mod_symb">@f$m_{t}@f$</td>
 *   <td class="mod_desc">The pole mass of the top quark in GeV.</td>
 * </tr>
 * <tr>
 *   <td class="mod_name">%mbottom</td>
 *   <td class="mod_symb">@f$m_{b}(m_b)@f$</td>
 *   <td class="mod_desc">The \f$\overline{\mathrm{MS}}\f$ scale-invariant mass of the bottom quark, \f$m_b(m_b)\f$, in GeV.</td>
 * </tr>
 * <tr>
 *   <td class="mod_name">%muc</td>
 *   <td class="mod_symb">@f$\mu_c@f$</td>
 *   <td class="mod_desc">The threshold between three- and four-flavour theory in GeV.</td>
 * </tr>
 * <tr>
 *   <td class="mod_name">%mub</td>
 *   <td class="mod_symb">@f$\mu_b@f$</td>
 *   <td class="mod_desc">The threshold between four- and five-flavour theory in GeV.</td>
 * </tr>
 * <tr>
 *   <td class="mod_name">%mut</td>
 *   <td class="mod_symb">@f$\mu_t@f$</td>
 *   <td class="mod_desc">The threshold between five- and six-flavour theory in GeV.</td>
 * </tr>
 * </table>
 *
 * @attention The parameters %AlsM and %MAls are not used in StandardModel
 * and the model classes inherited from it.
 *
 * 
 * The parameters below are associated with flavour observables
 * <table class="model">
 * <tr>
 *   <td class="mod_name">%MK0</td>
 *   <td class="mod_symb">@f$M_{K^0}@f$</td>
 *   <td class="mod_desc">The mass of the \f$ K^0 \f$ meson in GeV.</td>
 * </tr>
 * <tr>
 *   <td class="mod_name">%MKp</td>
 *   <td class="mod_symb">@f$M_{K^\pm}@f$</td>
 *   <td class="mod_desc">The mass of the \f$ K^\pm \f$ meson in GeV.</td>
 * </tr>
 * <tr>
 *   <td class="mod_name">%MKstar</td>
 *   <td class="mod_symb">@f$M_{K^*}@f$</td>
 *   <td class="mod_desc">The mass of the \f$ K^* \f$ meson in GeV.</td>
 * </tr>
 * <tr>
 *   <td class="mod_name">%Mphi</td>
 *   <td class="mod_symb">@f$M_{\phi}@f$</td>
 *   <td class="mod_desc">The mass of the \f$ \phi \f$ meson in GeV.</td>
 * </tr>
 * <tr>
 *   <td class="mod_name">%MD</td>
 *   <td class="mod_symb">@f$M_{D^0}@f$</td>
 *   <td class="mod_desc">The mass of the \f$ D^0 \f$ meson in GeV.</td>
 * </tr>
 * <tr>
 *   <td class="mod_name">%MBd</td>
 *   <td class="mod_symb">@f$M_{B_d}@f$</td>
 *   <td class="mod_desc">The mass of the \f$ B_d \f$ meson in GeV.</td>
 * </tr>
 * <tr>
 *   <td class="mod_name">%MBp</td>
 *   <td class="mod_symb">@f$M_{B^\pm}@f$</td>
 *   <td class="mod_desc">The mass of the \f$ B^\pm \f$ meson in GeV.</td>
 * </tr>
 * <tr>
 *   <td class="mod_name">%MBs</td>
 *   <td class="mod_symb">@f$M_{B_s}@f$</td>
 *   <td class="mod_desc">The mass of the \f$ B_s \f$ meson in GeV.</td>
 * </tr>
 * <tr>
 *   <td class="mod_name">%tKl</td>
 *   <td class="mod_symb">@f$\tau_{K_L}@f$</td>
 *   <td class="mod_desc">The lifetime of the \f$ K_L \f$ meson in \f$\mathrm{ps}^{-1} \f$.</td>
 * </tr>
 * <tr>
 *   <td class="mod_name">%tKp</td>
 *   <td class="mod_symb">@f$\tau_{K^\pm}@f$</td>
 *   <td class="mod_desc">The lifetime of the \f$ K^\pm \f$ meson in \f$\mathrm{ps}^{-1} \f$.</td>
 * </tr>
 * <tr>
 *   <td class="mod_name">%tKstar</td>
 *   <td class="mod_symb">@f$\tau_{K^*}@f$</td>
 *   <td class="mod_desc">The lifetime of the \f$ K^* \f$ meson in \f$\mathrm{ps}^{-1} \f$.</td>
 * </tr>
 * <tr>
 *   <td class="mod_name">%tphi</td>
 *   <td class="mod_symb">@f$\tau_{\phi}@f$</td>
 *   <td class="mod_desc">The lifetime of the \f$ \phi \f$ meson in \f$\mathrm{ps}^{-1} \f$.</td>
 * </tr>
 * <tr>
 *   <td class="mod_name">%tBd</td>
 *   <td class="mod_symb">@f$\tau_{B_d}@f$</td>
 *   <td class="mod_desc">The lifetime of the \f$ B_d \f$ meson in \f$\mathrm{ps}^{-1} \f$.</td>
 * </tr>
 * <tr>
 *   <td class="mod_name">%tBp</td>
 *   <td class="mod_symb">@f$\tau_{B^\pm}@f$</td>
 *   <td class="mod_desc">The lifetime of the \f$ B^\pm \f$ meson in \f$\mathrm{ps}^{-1} \f$.</td>
 * </tr>
 * <tr>
 *   <td class="mod_name">%tBs</td>
 *   <td class="mod_symb">@f$\tau_{B_s}@f$</td>
 *   <td class="mod_desc">The lifetime of the \f$ B_s \f$ meson in \f$\mathrm{ps}^{-1} \f$.</td>
 * </tr>
 * <tr>
 *   <td class="mod_name">%DGs_Gs</td>
 *   <td class="mod_symb"></td>
 *   <td class="mod_desc"></td>
 * </tr>
 * <tr>
 *   <td class="mod_name">%FK</td>
 *   <td class="mod_symb">@f$F_{K^-}@f$</td>
 *   <td class="mod_desc">The decay constant of the \f$ K^- \f$ meson in GeV.</td>
 * </tr>
 * <tr>
 *   <td class="mod_name">%FKstar</td>
 *   <td class="mod_symb">@f$F_{K^*}@f$</td>
 *   <td class="mod_desc">The decay constant of the \f$ K^* \f$ meson in GeV.</td>
 * </tr>
 * <tr>
 *   <td class="mod_name">%FKstarp</td>
 *   <td class="mod_symb">@f$F_{K^*}^{\perp}@f$</td>
 *   <td class="mod_desc">The decay constant of a transversely polarized \f$ K^* \f$ meson in GeV.</td>
 * </tr>
 * <tr>
 *   <td class="mod_name">%Fphi</td>
 *   <td class="mod_symb">@f$F_{\phi}@f$</td>
 *   <td class="mod_desc">The decay constant of the \f$ \phi \f$ meson in GeV.</td>
 * </tr>
 * <tr>
 *   <td class="mod_name">%FD</td>
 *   <td class="mod_symb">@f$F_{D^0}@f$</td>
 *   <td class="mod_desc">The decay constant of the \f$ D^0 \f$ meson in GeV.</td>
 * </tr>
 * <tr>
 *   <td class="mod_name">%FBs</td>
 *   <td class="mod_symb">@f$F_{B_s}@f$</td>
 *   <td class="mod_desc">The decay constant of the \f$ B_s \f$ meson in GeV.</td>
 * </tr>
 * <tr>
 *   <td class="mod_name">%FBsoFBd</td>
 *   <td class="mod_symb">@f$F_{B_d}/F_{B_d}@f$</td>
 *   <td class="mod_desc">The ratio \f$ F_{B_s}/F_{B_d} \f$ necessary to compute \f$ F_{B_s} \f$.</td>
 * </tr>
 * <tr>
 *   <td class="mod_name">%BK1 - %BK5</td>
 *   <td class="mod_symb">@f$B^1_{K} - B^5_{K}@f$</td>
 *   <td class="mod_desc">The bag parameter for \f$ O_1 - O_5\f$ in \f$ \Delta s = 2 \f$ processes in \f$ K^0 \f$.</td>
 * </tr>
 * <tr>
 *   <td class="mod_name">%BKscale</td>
 *   <td class="mod_symb">@f$\mu_K@f$</td>
 *   <td class="mod_desc">The scale at which the bag parameters are specified for the \f$ K^0 \f$ system.</td>
 * </tr>
 * <tr>
 *   <td class="mod_name">%BKscheme</td>
 *   <td class="mod_symb">@f$@f$</td>
 *   <td class="mod_desc">The scheme in which the bag parameters are specified for the \f$ K^0 \f$ system.</td>
 * </tr>
 * <tr>
 *   <td class="mod_name">%BD1 - %BD5</td>
 *   <td class="mod_symb">@f$B^1_{D} - B^5_{D}@f$</td>
 *   <td class="mod_desc">The bag parameter for \f$ O_1 - O_5\f$ in \f$ \Delta c = 2 \f$ processes in \f$ D^0 \f$.</td>
 * </tr>
 * <tr>
 *   <td class="mod_name">%BDscale</td>
 *   <td class="mod_symb">@f$\mu_D@f$</td>
 *   <td class="mod_desc">The scale at which the bag parameters are specified for the \f$ D_0 \f$ system.</td>
 * </tr>
 * <tr>
 *   <td class="mod_name">%BDscheme</td>
 *   <td class="mod_symb">@f$@f$</td>
 *   <td class="mod_desc">The scheme in which the bag parameters are specified for the \f$ D_0 \f$ system.</td>
 * </tr>
 * <tr>
 *   <td class="mod_name">%BBsoBBd</td>
 *   <td class="mod_symb">@f$B_{B_s}/B_{B_d}@f$</td>
 *   <td class="mod_desc">The ratio \f$ B_{B_s}/B_{B_d} \f$ necessary to compute \f$ B_{B_s} \f$.</td>
 * </tr>
 * <tr>
 *   <td class="mod_name">%BBs1 - %BBs5</td>
 *   <td class="mod_symb">@f$B^1_{B_s} - B^5_{B_s}@f$</td>
 *   <td class="mod_desc">The bag parameter for \f$ O_1 - O_5 \f$ in \f$ \Delta b = 2 \f$ processes in \f$ B_s \f$.</td>
 * </tr>
 * <tr>
 *   <td class="mod_name">%BBsscale</td>
 *   <td class="mod_symb">@f$\mu_{B_{s}}@f$</td>
 *   <td class="mod_desc">The scale at which the bag parameters are specified for the \f$ B_s \f$ system.</td>
 * </tr>
 * <tr>
 *   <td class="mod_name">%BBsscheme</td>
 *   <td class="mod_symb">@f$@f$</td>
 *   <td class="mod_desc">The scheme in which the bag parameters are specified for the \f$ B_s \f$ system.</td>
 * </tr>
 * <tr>
 *   <td class="mod_name">%BK(1/2)1 - %BK(1/2)10</td>
 *   <td class="mod_symb">@f$@f$</td>
 *   <td class="mod_desc"></td>
 * </tr>
 * <tr>
 *   <td class="mod_name">%BK(3/2)1 - %BK(3/2)10</td>
 *   <td class="mod_symb"></td>
 *   <td class="mod_desc"></td>
 * </tr>
 * <tr>
 *   <td class="mod_name">%BKd_scale</td>
 *   <td class="mod_symb">@f$@f$</td>
 *   <td class="mod_desc"></td>
 * </tr>
 * <tr>
 *   <td class="mod_name">%BKd_scheme</td>
 *   <td class="mod_symb">@f$@f$</td>
 *   <td class="mod_desc"></td>
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
 * <tr>
 *   <td class="mod_name">%Br_Kp_P0enu</td>
 *   <td class="mod_symb">@f$\mathrm{BR}(K^+\to\pi^0e^+\nu)@f$</td>
 *   <td class="mod_desc">The experimental value for the branching ratio of \f$K^+\to\pi^0e^+\nu\f$.</td>
 * </tr>
 * <tr>
 *   <td class="mod_name">%Br_Kp_munu</td>
 *   <td class="mod_symb">@f$\mathrm{BR}(K^+\to\mu^+\nu)@f$</td>
 *   <td class="mod_desc">The experimental value for the branching ratio of \f$K^+\to\mu^+\nu\f$.</td>
 * </tr>
 * <tr>
 *   <td class="mod_name">%Br_B_Xcenu</td>
 *   <td class="mod_symb">@f$\mathrm{BR}(B\to X_ce\nu)@f$</td>
 *   <td class="mod_desc">The experimental value for the branching ratio of \f$B\to X_c e\nu\f$.</td>
 * </tr>
 * <tr>
 *   <td class="mod_name">%DeltaP_cu</td>
 *   <td class="mod_symb">@f$@f$</td>
 *   <td class="mod_desc">The long-distance correction to the charm contribution of \f$K^+\to\pi^+\nu\bar{\nu}\f$.</td>
 * </tr>
 * <tr>
 *   <td class="mod_name">%IB_Kl</td>
 *   <td class="mod_symb">@f$@f$</td>
 *   <td class="mod_desc">the isospin breaking corrections between @f$K_L\to\pi^0\nu\bar{\nu}@f$ and \f$K^+\to\pi^0 e^+\nu\f$.</td>
 * </tr>
 * <tr>
 *   <td class="mod_name">%IB_Kp</td>
 *   <td class="mod_symb">@f$@f$</td>
 *   <td class="mod_desc">The isospin breaking corrections between @f$K^+\to\pi^+ \nu\bar{\nu}@f$ and \f$K^+\to\pi^0 e^+\nu\f$.</td>
 * </tr>
 * <tr>
 *   <td class="mod_name">%a_0V, %a_1V, %a_2V, %MRV</td>
 *   <td class="mod_symb">@f$a_0^V, a_1^V, a_2^V, \Delta m^V@f$</td>
 *   <td class="mod_desc">The fit parameters for the form factor @f$V@f$ of the @f$B\to K^*@f$.</td>
 * </tr>
 * <tr>
 *   <td class="mod_name">%a_0A0, %a_1A0, %a_2A0, %MRA0</td>
 *   <td class="mod_symb">@f$a_0^{A_0}, a_1^{A_0}, a_2^{A_0}, \Delta m^{A_0}@f$</td>
 *   <td class="mod_desc">The fit parameters for the form factor @f$A_0@f$ of the @f$B\to K^*@f$.</td>
 * </tr>
 * <tr>
 *   <td class="mod_name">%a_0A1, %a_1A1, %a_2A1, %MRA1</td>
 *   <td class="mod_symb">@f$a_0^{A_1}, a_1^{A_1}, a_2^{A_1}, \Delta m^{A_1}@f$</td>
 *   <td class="mod_desc">The fit parameters for the form factor @f$A_1@f$ of the @f$B\to K^*@f$.</td>
 * </tr>
 * <tr>
 *   <td class="mod_name">%a_0A12, %a_1A12, %a_2A12, %MRA12</td>
 *   <td class="mod_symb">@f$a_0^{A_{12}}, a_1^{A_{12}}, a_2^{A_{12}}, \Delta m^{A_{12}}@f$</td>
 *   <td class="mod_desc">The fit parameters for the form factor @f$A_{12}@f$ of the @f$B\to K^*@f$.</td>
 * </tr>
 * <tr>
 *   <td class="mod_name">%a_0T1, %a_1T1, %a_2T1, %MRA0</td>
 *   <td class="mod_symb">@f$a_0^{T_1}, a_1^{T_1}, a_2^{T_1}, \Delta m^{T_1}@f$</td>
 *   <td class="mod_desc">The fit parameters for the form factor @f$T_1@f$ of the @f$B\to K^*@f$.</td>
 * </tr>
 * <tr>
 *   <td class="mod_name">%a_0T2, %a_1T2, %a_2T2, %MRA1</td>
 *   <td class="mod_symb">@f$a_0^{T_2}, a_1^{T_2}, a_2^{T_2}, \Delta m^{T_2}@f$</td>
 *   <td class="mod_desc">The fit parameters for the form factor @f$T_2@f$ of the @f$B\to K^*@f$.</td>
 * </tr>
 * <tr>
 *   <td class="mod_name">%a_0T23, %a_1T23, %a_2T23, %MRA1</td>
 *   <td class="mod_symb">@f$a_0^{T_{23}}, a_1^{T_{23}}, a_2^{T_{23}}, \Delta m^{T_{23}}@f$</td>
 *   <td class="mod_desc">The fit parameters for the form factor @f$T_{23}@f$ of the @f$B\to K^*@f$.</td>
 * </tr>
 * <tr>
 *   <td class="mod_name">%a_0Vphi, %a_1Vphi, %a_2Vphi, %MRVphi</td>
 *   <td class="mod_symb">@f$a_0^V, a_1^V, a_2^V, \Delta m^V@f$</td>
 *   <td class="mod_desc">The fit parameters for the form factor @f$V@f$ of the @f$B\to\phi@f$.</td>
 * </tr>
 * <tr>
 *   <td class="mod_name">%a_0A0phi, %a_1A0phi, %a_2A0phi, %MRA0phi</td>
 *   <td class="mod_symb">@f$a_0^{A_0}, a_1^{A_0}, a_2^{A_0}, \Delta m^{A_0}@f$</td>
 *   <td class="mod_desc">The fit parameters for the form factor @f$A_0@f$ of the @f$B\to\phi@f$.</td>
 * </tr>
 * <tr>
 *   <td class="mod_name">%a_0A1phi, %a_1A1phi, %a_2A1phi, %MRA1phi</td>
 *   <td class="mod_symb">@f$a_0^{A_1}, a_1^{A_1}, a_2^{A_1}, \Delta m^{A_1}@f$</td>
 *   <td class="mod_desc">The fit parameters for the form factor @f$A_1@f$ of the @f$B\to\phi@f$.</td>
 * </tr>
 * <tr>
 *   <td class="mod_name">%a_0A1phi, %a_1A1phi, %a_2A1phi, %MRA1phi</td>
 *   <td class="mod_symb">@f$a_0^{A_{12}}, a_1^{A_{12}}, a_2^{A_{12}}, \Delta m^{A_{12}}@f$</td>
 *   <td class="mod_desc">The fit parameters for the form factor @f$A_{12}@f$ of the @f$B\to\phi@f$.</td>
 * </tr>
 * <tr>
 *   <td class="mod_name">%a_0T1phi, %a_1T1phi, %a_2T1phi, %MRA0phi</td>
 *   <td class="mod_symb">@f$a_0^{T_1}, a_1^{T_1}, a_2^{T_1}, \Delta m^{T_1}@f$</td>
 *   <td class="mod_desc">The fit parameters for the form factor @f$T_1@f$ of the @f$B\to\phi@f$.</td>
 * </tr>
 * <tr>
 *   <td class="mod_name">%a_0T2phi, %a_1T2phi, %a_2T2phi, %MRA1phi</td>
 *   <td class="mod_symb">@f$a_0^{T_2}, a_1^{T_2}, a_2^{T_2}, \Delta m^{T_2}@f$</td>
 *   <td class="mod_desc">The fit parameters for the form factor @f$T_2@f$ of the @f$B\to\phi@f$.</td>
 * </tr>
 * <tr>
 *   <td class="mod_name">%a_0T23phi, %a_1T23phi, %a_2T23phi, %MRA1phi</td>
 *   <td class="mod_symb">@f$a_0^{T_{23}}, a_1^{T_{23}}, a_2^{T_{23}}, \Delta m^{T_{23}}@f$</td>
 *   <td class="mod_desc">The fit parameters for the form factor @f$T_{23}@f$ of the @f$B\to\phi@f$.</td>
 * </tr>
 * <tr>
 *   <td class="mod_name">%reh_0, %reh_0_1, %reh_0_2</td>
 *   <td class="mod_symb">@f$\mathrm{Re}h_0^{(0)}, \mathrm{Re}h_0^{(1)}, \mathrm{Re}h_0^{(2)}@f$</td>
 *   <td class="mod_desc">The constant, linear and quadratic terms of the real part of the hadronic parameter @f$h_0@f$ of the @f$B\to K^*@f$.</td>
 * </tr>
 * <tr>
 *   <td class="mod_name">%imh_0, %imh_0_1, %imh_0_2</td>
 *   <td class="mod_symb">@f$\mathrm{Im}h_0^{(0)}, \mathrm{Im}h_0^{(1)}, \mathrm{Im}h_0^{(2)}@f$</td>
 *   <td class="mod_desc">The constant, linear and quadratic terms of the immaginary part of the hadronic parameter @f$h_0@f$ of the @f$B\to K^*@f$.</td>
 * </tr>
 * <tr>
 *   <td class="mod_name">%reh_1, %reh_1_1, %reh_1_2</td>
 *   <td class="mod_symb">@f$\mathrm{Re}h_+^{(0)}, \mathrm{Re}h_+^{(1)}, \mathrm{Re}h_+^{(2)}@f$</td>
 *   <td class="mod_desc">The constant, linear and quadratic terms of the real part of the hadronic parameter @f$h_+@f$ of the @f$B\to K^*@f$.</td>
 * </tr>
 * <tr>
 *   <td class="mod_name">%imh_1, %imh_1_1, %imh_1_2</td>
 *   <td class="mod_symb">@f$\mathrm{Im}h_+^{(0)}, \mathrm{Im}h_+^{(1)}, \mathrm{Im}h_+^{(2)}@f$</td>
 *   <td class="mod_desc">The constant, linear and quadratic terms of the immaginary part of the hadronic parameter @f$h_+@f$ of the @f$B\to K^*@f$.</td>
 * </tr>
 * <tr>
 *   <td class="mod_name">%reh_2, %reh_2_1, %reh_2_2</td>
 *   <td class="mod_symb">@f$\mathrm{Re}h_-^{(0)}, \mathrm{Re}h_-^{(1)}, \mathrm{Re}h_-^{(2)}@f$</td>
 *   <td class="mod_desc">The constant, linear and quadratic terms of the real part of the hadronic parameter @f$h_-@f$ of the @f$B\to K^*@f$.</td>
 * </tr>
 * <tr>
 *   <td class="mod_name">%imh_2, %imh_2_1, %imh_2_2</td>
 *   <td class="mod_symb">@f$\mathrm{Im}h_-^{(0)}, \mathrm{Im}h_-^{(1)}, \mathrm{Im}h_-^{(2)}@f$</td>
 *   <td class="mod_desc">The constant, linear and quadratic terms of the immaginary part of the hadronic parameter @f$h_-@f$ of the @f$B\to K^*@f$.</td>
 * </tr>
 * <tr>
 *   <td class="mod_name">%alpha1kst, %alpha2kst</td>
 *   <td class="mod_symb">@f$\alpha_1(\bar{K}^*), \alpha_2(\bar{K}^*)@f$</td>
 *   <td class="mod_desc">The Gegenbauer coefficients for the @f$\bar{K}^*@f$ meson.</td>
 * </tr>
 * <tr>
 *   <td class="mod_name">%lambdaB</td>
 *   <td class="mod_symb">@f$\Lambda_{B,+}@f$</td>
 *   <td class="mod_desc">The integrated leading twist light-cone distribution amplitudes of the B meson divided by the integral variable.</td>
 * </tr>
 * <tr>
 *   <td class="mod_name">%r_1_fplus, %r_2_fplus, %m_fit2_fplus</td>
 *   <td class="mod_symb">@f$r_1^{f_+}, r_2^{f_+}, m_{fit}^{2,f_+}@f$</td>
 *   <td class="mod_desc">The fit parameters for the LCSR form factor @f$f_+@f$ of the @f$B\to K@f$.</td>
 * </tr>
 * <tr>
 *   <td class="mod_name">%r_1_fT, %r_2_fT, %m_fit2_fT</td>
 *   <td class="mod_symb">@f$r_1^{f_T}, r_2^{f_T}, m_{fit}^{2,f_T}@f$</td>
 *   <td class="mod_desc">The fit parameters for the LCSR form factor @f$f_T@f$ of the @f$B\to K@f$.</td>
 * </tr>
 * <tr>
 *   <td class="mod_name">%r_2_f0, %m_fit2_f0</td>
 *   <td class="mod_symb">@f$r_2^{f_0}, m_{fit}^{2,f_0}@f$</td>
 *   <td class="mod_desc">The fit parameters for the LCSR form factor @f$f_0@f$ of the @f$B\to K@f$.</td>
 * </tr>
 * <tr>
 *   <td class="mod_name">%reh_0_MP, %reh_0_1_MP</td>
 *   <td class="mod_symb">@f$\mathrm{Re}h_0^{(0)}, \mathrm{Re}h_0^{(1)}@f$</td>
 *   <td class="mod_desc">The constant and linear terms of the real part of the hadronic parameter @f$h_0@f$ of the @f$B\to K@f$.</td>
 * </tr>
 * <tr>
 *   <td class="mod_name">%imh_0_MP, %imh_0_1_MP</td>
 *   <td class="mod_symb">@f$\mathrm{Im}h_0^{(0)}, \mathrm{Im}h_0^{(1)}@f$</td>
 *   <td class="mod_desc">The constant and linear terms of the immaginary part of the hadronic parameter @f$h_0@f$ of the @f$B\to K@f$.</td>
 * </tr>
 * <tr>
 *   <td class="mod_name">%bsgamma_E0</td>
 *   <td class="mod_symb">@f$E_0^{b\to s\gamma}@f$</td>
 *   <td class="mod_desc">The experimental energy cutoff @f$E_0@f$ of the @f$b\to s\gamma@f$.</td>
 * </tr>
 * <tr>
 *   <td class="mod_name">%BLNPcorr</td>
 *   <td class="mod_symb">@f$N_{b\to s\gamma}@f$</td>
 *   <td class="mod_desc">The non perturbative uncertainty associated to the @f$b\to s\gamma@f$ BR.</td>
 * </tr>
 * <tr>
 *   <td class="mod_name">%Gambino_mukin</td>
 *   <td class="mod_symb">@f$\mu^{\rm kin}@f$</td>
 *   <td class="mod_desc">The mass scale for the b quark in the kinetic scheme, employed in the @f$b\to s\gamma@f$ BR.</td>
 * </tr>
 * <tr>
 *   <td class="mod_name">%Gambino_BRsem</td>
 *   <td class="mod_symb">@f$\mathrm{BR}(B\to X_ce\nu)^{\rm Gambino}@f$</td>
 *   <td class="mod_desc">The branching ratio of \f$B\to X_c e\nu\f$ fitted by Paolo Gambino, employed in the @f$b\to s\gamma@f$ BR.</td>
 * </tr>
 * <tr>
 *   <td class="mod_name">%Gambino_Mbkin</td>
 *   <td class="mod_symb">@f$m_b^{\rm kin,\,Gambino}@f$</td>
 *   <td class="mod_desc">The b quark mass in the kinetic scheme fitted by Paolo Gambino, employed in the @f$b\to s\gamma@f$ BR.</td>
 * </tr>
 * <tr>
 *   <td class="mod_name">%Gambino_Mcatmuc</td>
 *   <td class="mod_symb">@f$m_c(\mu_c)^{\rm Gambino}@f$</td>
 *   <td class="mod_desc">The c quark at \f$\mu_c\f$ fitted by Paolo Gambino, employed in the @f$b\to s\gamma@f$ BR.</td>
 * </tr>
 * <tr>
 *   <td class="mod_name">%Gambino_mupi2, %Gambino_rhoD3, %Gambino_muG2, %Gambino_rhoLS3</td>
 *   <td class="mod_symb">@f$\mu_{\pi}^{2,\,\rm Gambino}, \rho_{D}^{3,\,\rm Gambino}, \mu_{G}^{2,\,\rm Gambino}, \rho_{LS}^{3,\,\rm Gambino}@f$</td>
 *   <td class="mod_desc">The B meson expecation values for the relevant dim. 5 and 6 operators fitted by Paolo Gambino, employed in the @f$b\to s\gamma@f$ BR.</td>
 * </tr>
 * </table>
 *
 * The set of the model parameters are initialized and updated with the methods
 * Init() and Update(), respectively, where the former calls the latter.
 * In Update(), the methods PreUpdate() and PostUpdate() are called to run all
 * the procedures that are need to be executed before and after the model parameters
 * are updated. The \f$\overline{\mathrm{MS}}\f$ mass for the top quark is computed and the scale set
 * in PostUpdate() with the updated parameters. Inside the Update() method, the
 * individual model parameter is assigned with the protected member function
 * setParameter().
 *
 * <H3>Computation of the strong coupling constant \f$\alpha_s\f$:  </H3>
 *
 * The strong coupling constant @f$\alpha_s@f$ at an arbitrary scale can be
 * computed with the member functions: 
 *
 * @li AlsWithInit(const double mu, const double alsi, const double mu_i, const orders order),
 * @li AlsWithLambda(const double mu, const orders order), 
 *
 * where another function
 * 
 * @li Als(const double mu, const orders order = FULLNLO)
 *
 * calls %AlsWithInit() for order=LO/FULLNLO(NLO),
 * and %AlsWithLambda() for order=FULLNNLO(NNLO).
 *
 * The function %AlsWithInit() computes @f$\alpha_s(\mu)@f$ with a given initial
 * value @f$\alpha_s(\mu_i)@f$: 
 * @f[
 *   \alpha_s(\mu)=\frac{\alpha_s(\mu_i)}{v(\mu)},\qquad
 *   \alpha_s(\mu)=\frac{\alpha_s(\mu_i)}{v(\mu)}
 *   \left[1-\frac{\beta_1}{\beta_0}\frac{\alpha_s(\mu_i)}{4\pi}
 *   \frac{\ln v(\mu)}{v(\mu)}\right],
 * @f]
 * at LO and FULLNLO (NLO) respectively, where
 * @f[
 *   v(\mu) = 1- \beta_0\frac{\alpha_s(\mu_i)}{2\pi}\ln\frac{\mu_i}{\mu},
 * @f]
 * and the one-loop and two-loop beta functions are given by 
 * @f[
 *  \beta_0 = \frac{11N_c-2N_f}{3},\qquad
 *  \beta_1 = \frac{34}{3}N_c^2-\frac{10}{3}N_cN_f-2C_FN_f.
 * @f]
 *
 * The function %AlsWithLambda() computes @f$\alpha_s(\mu)@f$ with the use of
 * \f$\Lambda_{\rm QCD}\f$, where the value of \f$\Lambda_{\rm QCD}\f$ for \f$N_f=5\f$
 * is derived from @f$\alpha_s(M_{\alpha_s})@f$ by solving
 * the equation (see e.g., @cite Chetyrkin:1997sg, which follows the convention
 * in @cite Buras:1977qg for \f$\ln(\mu^2/\Lambda^2)\f$): 
 * @f[
 *  \frac{\alpha_s(\mu)}{4\pi}=\frac{1}{\beta_0L}-\frac{\beta_1\ln L}{\beta_0^3L^2}+
 *  \frac{1}{\left(\beta_0L\right)^3}\left[\frac{\beta_1^2}{\beta_0^2}\left(\ln^2L-\ln L-1\right)+
 *  \frac{\beta_2}{\beta_0}\right]+O\left(\frac{\ln^3L}{L^4}\right),
 * @f]
 * where @f$L\equiv\ln(\mu^2/\Lambda^2)@f$, and the three-loop beta function
 * is given by
 * @f[
 *  \beta_2 = \frac{2857}{54}N_c^3+C_F^2N_f-\frac{205}{18}C_FN_cN_f
 *  -\frac{1415}{54}N_c^2N_f+\frac{11}{9}C_FN_f^2 + \frac{79}{54}N_cN_f^2.
 * @f]
 * For \f$N_f < 5\f$, \f$\Lambda_{QCD}\f$ can be obtained by solving the following 
 * matching condition at \f$\mu\f$ @cite Chetyrkin:1997sg @cite Chetyrkin:2000yt :
 * @f{eqnarray}
 *  \beta_0^{'}\ln\frac{\Lambda^{'2}}{\Lambda^2}
 *  &=& (\beta_0^{'}-\beta_0)L+\left(\frac{\beta_1^{'}}{\beta_0^{'}}-
 *  \frac{\beta_1}{\beta_0}\right)\ln L-\frac{\beta_1^{'}}{\beta_0^{'}}\ln\frac{\beta_0^{'}}{\beta_0} - C_1
 *  \\
 *  &&+
 *  \frac{1}{\beta_0L}\left[\frac{\beta_1}{\beta_0}\left(\frac{\beta_1^{'}}{\beta_0^{'}}-\frac{\beta_1}{\beta_0}
 *  \right)\ln L + \frac{\beta_1^{'2}}{\beta_0^{'2}}-\frac{\beta_1^2}{\beta_0^2}-\frac{\beta_2^{'2}}{\beta_0^{'2}}
 *  +\frac{\beta_2^2}{\beta_0^2}+\frac{\beta_1^{'}}{\beta_0^{'}}C_1-C_1^2-C_2\right]+
 *  O\left(\frac{\ln^2L}{L^2}\right),
 * @f}
 * where the primed (unprimed) quantities refer to those pertaining to \f$N_f-1\f$ \f$(N_f)\f$. 
 * The terms \f$C_1\f$ and \f$C_2\f$ are given by
 * @f[
 *  C_1 = \frac{2}{3}\ln\frac{\mu^2}{\mu_f^2},\qquad
 *  C_2 = -16\left(\frac{1}{36}\ln^2\frac{\mu^2}{\mu_f^2}-\frac{19}{24}\ln\frac{\mu^2}{\mu_f^2}+\frac{11}{72}\right).
 * @f]
 * where \f$\mu_f=m_f(\mu_f)\f$ denotes the \f$\overline{\mathrm{MS}}\f$ invariant
 * mass of the \f$N_f\f$-th flavour.
 * Moreover, the matching condition at a flavour threshold is given by
 * @f[
 *   \alpha_s^{(N_f-1)}(\mu) = (\zeta_g^f)^2\alpha_s^{(N_f)}(\mu),
 * @f]
 * where
 * @f[
 *   (\zeta_g^f)^2 = 1+\frac{\alpha_s^{(N_f)}(\mu)}{\pi}\left(-\frac{1}{6}\ln\frac{\mu^2}{\mu_f^2}\right)+
 *   \left(\frac{\alpha_s^{(N_f)}(\mu)}{\pi}\right)^2\left(\frac{1}{36}\ln^2\frac{\mu^2}{\mu_f^2} -
 *   \frac{19}{24}\ln\frac{\mu^2}{\mu_f^2} + \frac{11}{72}\right) + O\left(\left(\frac{\alpha_s}{\pi}
 *   \right)^3\right).
 * @f]
 * 
 * For the top quark mass, the pole mass @f$m_t@f$ is used as an input instead of the \f$\overline{\mathrm{MS}}\f$
 * invariant mass. Then the \f$\overline{\mathrm{MS}}\f$ invariant mass is computed from the pole mass with the
 * computed value of \f$\alpha_s^{(6)}(m_t)\f$, which in turn is computed from \f$\alpha_s^{(5)}(M_Z)\f$. Hence, the
 * matching condition from \f$N_f = 6\f$ to \f$N_f = 5\f$ has to be free from the \f$\overline{\mathrm{MS}}\f$
 * invariant mass of the top quark. In this case we use the following matching condition
 * @f[
 *   \alpha_s^{(6)}(\mu_t) = \frac{1}{(\zeta_g^{OS,t})^2}\alpha_s^{(5)}(\mu_t),
 * @f]
 * where
 * @f[
 *  \frac{1}{(\zeta_g^{OS,t})^2} = 1+\frac{\alpha_s^{(5)}(\mu)}{\pi}\left(\frac{1}{6}\ln\frac{\mu^2}{m_t^2}\right)+
 *  \left(\frac{\alpha_s^{(5)}(\mu)}{\pi}\right)^2\left(\frac{1}{36}\ln^2\frac{\mu^2}{m_t^2} +
 *  \frac{19}{24}\ln\frac{\mu^2}{m_t^2} + \frac{7}{24}\right) + O\left(\left(\frac{\alpha_s}{\pi}\right)^3\right).
 * @f]
 * Besides, \f$\Lambda_{QCD}\f$ for \f$N_f=6\f$ is derived from that for
 * \f$N_f=5\f$ with the relation:
 * @f{eqnarray}
 *  \beta_0\ln\frac{\Lambda^{2}}{\Lambda^{'2}} 
 *  &=& (\beta_0-\beta_0^{'})L^{'}+\left(\frac{\beta_1}{\beta_0}-
 *  \frac{\beta_1^{'}}{\beta_0^{'}}\right)\ln L^{'}-\frac{\beta_1}{\beta_0}\ln\frac{\beta_0}{\beta_0^{'}}-C_1^{'}
 *  \\
 *  && +
 *  \frac{1}{\beta_0^{'}L^{'}}\left[\frac{\beta_1^{'}}{\beta_0^{'}}\left(\frac{\beta_1}{\beta_0}-\frac{\beta_1^{'}}
 *  {\beta_0^{'}}\right)\ln L^{'}+\frac{\beta_1^{2}}{\beta_0^{2}}-\frac{\beta_1^{'2}}{\beta_0^{'2}}-
 *  \frac{\beta_2^{2}}{\beta_0^{2}}+\frac{\beta_2^{'2}}{\beta_0^{'2}}+\frac{\beta_1}{\beta_0}C_1^{'}-C_1^{'2}-
 *  C_2^{'}\right]+O\left(\frac{\ln^2L}{L^2}\right),
 * @f}
 * where
 * @f[
 *  C_1^{'} = -\frac{2}{3}\ln\frac{\mu^2}{m_t^2},\qquad
 *  C_2^{'} = -16\left(\frac{1}{36}\ln^2\frac{\mu^2}{m_t^2}+\frac{19}{24}\ln\frac{\mu^2}{m_t^2}+\frac{7}{24}\right).
 * @f]
 *
 * 
 * <H3>Computation of the running of the quark masses: </H3>
 *
 * In the \f$\overline{\mathrm{MS}}\f$ scheme the quark mass at a scale \f$\mu\f$ is given by 
 * (see e.g., @cite Chetyrkin:1997dh)
 * @f[
 *  m_q(\mu) = m_q(\mu_0)\left[\frac{\alpha_s(\mu)}{\alpha_s(\mu_0)}\right]^{\frac{\gamma^{(0)}_m}{2\beta_0}}
 *  \left\{1+A_1\frac{\alpha_s(\mu)-\alpha_s(\mu_0)}{4\pi}+\frac{A_1^2}{2}\left(\frac{\alpha_s(\mu)-\alpha_s(\mu_0)}{4\pi}
 *  \right)^2+\frac{A_2}{2}\left[\left(\frac{\alpha_s(\mu)}{4\pi}\right)^2-\left(\frac{\alpha_s(\mu_0)}{4\pi}
 *  \right)^2\right]\right\},
 * @f]
 * where
 * @f[
 *  A_1 = \frac{\gamma^{(1)}_m}{2\beta_0} - \frac{\beta_1\gamma^{(0)}_m}{2\beta_0^2},\qquad
 *  A_2 = \frac{\beta_1^2\gamma^{(0)}_m}{2\beta_0^3}-\frac{\beta_2\gamma^{(0)}_m}{2\beta_0^2}-
 *  \frac{\beta_1\gamma^{(1)}_m}{2\beta_0^2}+\frac{\gamma^{(2)}_m}{2\beta_0},
 * @f]
 * and 
 * @f[
 *  \gamma_m^{(0)} = 6C_F,\qquad
 *  \gamma_m^{(1)} = C_F\left(3C_F+\frac{97}{3}N_c-\frac{10}{3}N_f\right),\\
 *  \gamma_m^{(2)} = 129C_F^3-\frac{129}{2}C_F^2N_c+\frac{11413}{54}C_FN_c^2+C_F^2N_f(-46+48\zeta(3))+
 *  C_FN_CN_f\left(-\frac{556}{27}-48\zeta(3)\right)-\frac{70}{27}C_FN_f^2,
 * @f]
 * in the \f$\overline{\mathrm{MS}}\f$ scheme.
 * The threshold conditions are given by @cite Chetyrkin:2000yt
 * @f[
 *   m_q^{(N_f-1)}(\mu)=\zeta_m^f\, m_q^{(N_f)}(\mu),\qquad
 *   m_q^{(N_f)}(\mu)=\frac{1}{\zeta_m^f}\, m_q^{(N_f-1)}(\mu),
 * @f]
 * where
 * @f[
 *  \zeta^f_m=1+\left(\frac{\alpha_s^{(N_f)}(\mu)}{\pi}\right)^2\left(\frac{1}{12}\ln^2\frac{\mu^2}{\mu_f^2}
 *   -\frac{5}{36}\ln\frac{\mu^2}{\mu_f^2}+\frac{89}{432}\right),\\
 *   \frac{1}{\zeta^f_m}=1+\left(\frac{\alpha_s^{(N_f-1)}(\mu)}{\pi}\right)^2\left(-\frac{1}{12}\ln^2\frac{\mu^2}
 *   {\mu_f^2}+\frac{5}{36}\ln\frac{\mu^2}{\mu_f^2}-\frac{89}{432}\right), 
 * @f]
 * with \f$\mu_f=m_f(\mu_f)\f$.
 * 
 *
 * <H3>Pole mass vs. \f$ \overline{\mathrm{MS}} \f$ scale invariant mass:</H3>
 *
 * The pole mass \f$M_q\f$ of a heavy quark is related to the \f$\overline{\mathrm{MS}}\f$ scale invariant mass
 * \f$m_q(m_q)\f$ as @cite Chetyrkin:2000yt
 * @f[
 *  m_q(m_q)=M_q\left\{1-\frac{4}{3}\frac{\alpha_s(M_q)}{\pi}+\left[-\frac{2251}{268}-2\zeta(2)
 *  -\frac{2}{3}\zeta(2)\ln2+\frac{\zeta(3)}{6}+\frac{n_l}{3}\left(\zeta(2)+\frac{71}{48}\right)
 *  -\frac{4}{3}\sum_{1\le i \le n_l}\Delta\left(\frac{M_i}{M_q}\right)\right]\left(\frac{\alpha_s(M_q)}{\pi}
 *  \right)^2\right\},
 * @f]
 * where \f$n_l\f$ is the number of the light quarks, and \f$\Delta(x)\f$ is given by
 * @f[
 *  \Delta(x)=\frac{\pi^2}{8}x-0.597x^2+0.230x^3+O\left(x^4\right).
 * @f]
 * It should be noted that the above formula requires the light quark pole masses which are not well defined.
 * The pole mass in terms of the \f$\overline{\mathrm{MS}}\f$ mass is given by
 * @f[
 *  M_q=m_q(m_q)\left\{1+\frac{4}{3}\frac{\alpha_s(m_q)}{\pi}+\left[\frac{307}{32}+2\zeta(2)
 *  +\frac{2}{3}\zeta(2)\ln2-\frac{\zeta(3)}{6}-\frac{n_l}{3}\left(\zeta(2)+\frac{71}{48}\right)
 *  +\frac{4}{3}\sum_{1\le i \le n_l}\Delta\left(\frac{m_i}{m_q}\right)\right]\left(\frac{\alpha_s(m_q)}{\pi}
 *  \right)^2\right\}
 * @f]
 * To get the pole mass of the light quarks we solve the above equation numerically and derive it from the 
 * corresponding \f$\overline{\mathrm{MS}}\f$ mass. 
 * 
 */
class QCD : public Model {
public:

    /**
     * @brief An enum type for mesons.
     */
    enum meson {
        P_0, /**< @f$\pi^0@f$ meson */
        P_P, /**< @f$\pi^\pm@f$ meson */
        K_0, /**< @f$K^0@f$ meson */
        K_P, /**< @f$K^\pm@f$ meson */
        D_0, /**< @f$D^0@f$ meson */
        B_D, /**< @f$B_d@f$ meson */
        B_P, /**< @f$B^\pm@f$ meson */
        B_S, /**< @f$B_s@f$ meson */
        PHI, /**< @f$\phi@f$ meson */
        K_star, /**< @f$K^*@f$ meson */
        MESON_END /**< The size of this enum. */
    };

    /**
     * @brief An enum type for quarks.
     */
    enum quark {
        UP, /**< Up quark */
        DOWN, /**< Down quark */
        CHARM, /**< Charm quark */
        STRANGE, /**< Strange quark */
        TOP, /**< Top quark */
        BOTTOM /**< Bottom quark */
    };

    static const int NQCDvars = 186; ///< The number of model parameters in %QCD. 

    /**
     * @brief An array containing the labels under which all %QCD parameters are stored
     * in a vector of ModelParameter via InputParser::ReadParameters().
     */
    static const std::string QCDvars[NQCDvars];

    /**
     * @brief Constructor.
     */
    QCD();

    /**
     * @brief Converts an object of the enum type "orders" to the corresponding string.
     * @param[in] order an object of the enum type "orders"
     * @return the string of the given "order"
     */
    std::string orderToString(const orders order) const;

    ////////////////////////////////////////////////////////////////////////
    // Parameters

    /**
     * @brief Initializes the %QCD parameters found in the argument.
     * @param[in] DPars a map containing the parameters (all as double) to be used in Monte Carlo
     */
    virtual bool Init(const std::map<std::string, double>& DPars);

    /**
     * @brief The pre-update method for %QCD
     * @details This method resets the internal flags #requireYu, #requireYd,
     * #computeBd, #computeFBd and #computemt before updating the model parameters with the method Update().
     * @return a boolean that is true if the execution is successful
     */
    virtual bool PreUpdate();

    /**
     * @brief The update method for %QCD.
     * @details This method updates all the model parameters with given DPars.
     * @param[in] DPars a map of the parameters that are being updated in the Monte Carlo run
     * (including parameters that are varied and those that are held constant)
     * @return a boolean that is true if the execution is successful
     */
    virtual bool Update(const std::map<std::string, double>& DPars);

    /**
     * @brief The post-update method for %QCD.
     * @details This method runs all the procedures that are need to be executed
     * after the model is successfully updated. This includes 
     * \li computing the decay constant \f$F_{B_D}\f$ from \f$F_{B_s}\f$
     * \li computing the bag parameters \f$B_{B_d}\f$ from \f$B_{B_s}\f$
     * \li computing the \f$\overline{\rm MS}\f$ mass of the top quark at the \f$\overline{\rm MS}\f$ mass,
     * \f$m_t^{\overline{\rm MS}}(m_t^{\overline{\rm MS}})\f$ and setting the scale at the same value.
     * @return a boolean that is true if the execution is successful
     */
    virtual bool PostUpdate();

    /**
     * @brief A method to check if all the mandatory parameters for %StandardModel
     * have been provided in model initialization.
     * @param[in] DPars a map of the parameters that are being updated in the Monte Carlo run
     * (including parameters that are varied and those that are held constant)
     * @return a boolean that is true if the execution is successful
     */
    virtual bool CheckParameters(const std::map<std::string, double>& DPars);


    ////////////////////////////////////////////////////////////////////////
    // Flags

    /**
     * @brief A method to set a flag of %QCD.
     * @param[in] name name of a model flag
     * @param[in] value the boolean to be assigned to the flag specified by name
     * @return a boolean that is true if the execution is successful
     */
    virtual bool setFlag(const std::string name, const bool value);

    /**
     * @brief A method to set a flag of %QCD.
     * @param[in] name name of a model flag
     * @param[in] value the string to be assigned to the flag specified by name
     * @return a boolean that is true if the execution is successful
     */
    virtual bool setFlagStr(const std::string name, const std::string value);

    /**
     * @brief A method to check the sanity of the set of model flags.
     * @return a boolean that is true if the set of model flags is sane
     */
    virtual bool CheckFlags() const;


    ////////////////////////////////////////////////////////////////////////
    // get and set methods for class members

    /**
     * @brief A get method to access the number of colours \f$N_c\f$.
     * @return the number of colours
     */
    double getNc() const
    {
        return Nc;
    }

    /**
     * @brief A set method to change the number of colours \f$N_c\f$.
     * @param[in] Nc the number of colours
     */
    void setNc(double Nc)
    {
        this->Nc = Nc;
    }

    /**
     * @brief A get method to access a meson as an object of the type Meson.
     * @param[in] m the name of a meson
     * @return the object of the meson specified in the argument
     */
    Meson getMesons(const meson m) const
    {
        return mesons[m];
    }

    /**
     * @brief A get method to access a quark as an object of the type Particle.
     * @param[in] q the name of a quark
     * @return the object of the quark found in the argument
     */
    Particle getQuarks(const quark q) const
    {
        return quarks[q];
    }

    /**
     * @brief A get method to access the value of \f$\alpha_s(M_{\alpha_s})\f$.
     * @return the strong coupling constant at @f$M_{\alpha_s}@f$, @f$\alpha_s(M_{\alpha_s})@f$
     */
    double getAlsM() const
    {
        return AlsM;
    }

    /**
     * @brief A get method to access the mass scale \f$M_{\alpha_s}\f$ at which the strong coupling constant measurement is provided.
     * @return the mass scale in GeV @f$M_{\alpha_s}@f$ at which the strong coupling constant measurement is provided
     */
    double getMAls() const
    {
        return MAls;
    }

    /**
     * @brief A get method to access the threshold between six- and five-flavour theory in GeV.
     * @return the threshold \f$\mu_t\f$
     */
    double getMut() const
    {
        return mut;
    }

    /**
     * @brief A get method to access the threshold between five- and four-flavour theory in GeV.
     * @return the threshold \f$\mu_b\f$
     */
    double getMub() const
    {
        return mub;
    }

    /**
     * @brief A get method to access the threshold between four- and three-flavour theory in GeV.
     * @return the threshold \f$\mu_c\f$
     */
    double getMuc() const
    {
        return muc;
    }

    /**
     * @brief A get method to access the pole mass of the top quark.
     * @return the pole mass of the top quark \f$m_t^{pole}\f$
     */
    double getMtpole() const
    {
        return mtpole;
    }

    /**
     * @brief A get method to access the Casimir factor of %QCD.
     * @return the Casimir factor
     */
    double getCF() const
    {
        return CF;
    }

    /**
     * @brief For getting the bag parameters corresponding
     * to the operator basis \f$O_1 -O_5\f$ in \f$\Delta b = 2\f$
     * process in the \f$B_d\f$ meson system.
     * @return the vector of bag parameters
     */
    BParameter getBBd() const
    {
        return BBd;
    }

    /**
     * @brief For getting the bag parameters corresponding
     * to the operator basis \f$O_1 -O_5\f$ in \f$\Delta b = 2\f$
     * process in the \f$B_s\f$ meson system.
     * @return the vector of bag parameters
     */
    BParameter getBBs() const
    {
        return BBs;
    }

    /**
     * @brief For getting the bag parameters corresponding
     * to the operator basis \f$O_1 -O_5\f$ in \f$\Delta c = 2\f$
     * process in the \f$D^0\f$ meson system.
     * @return the vector of bag parameters
     */
    BParameter getBD() const
    {
        return BD;
    }

    /**
     * @brief For getting the bag parameters corresponding
     * to the operator basis \f$O_1 -O_5\f$ in \f$\Delta s = 2\f$
     * process in the \f$K^0\f$ meson system.
     * @return the vector of bag parameters
     */
    BParameter getBK() const
    {
        return BK;
    }

    /**
     * @return
     */
    BParameter getBKd1() const
    {
        return BKd1;
    }

    /**
     * @return
     */
    BParameter getBKd3() const
    {
        return BKd3;
    }

    /**
     * @return the experimental value of the real part of the amplitude for
     * @f$K^0\to\pi\pi@f$ with @f$\Delta I=0@f$
     */
    double getReA0_Kd() const
    {
        return ReA0_Kd;
    }

    /**
     * @return the experimental value of the real part of the amplitude for
     * @f$K^0\to\pi\pi@f$ with @f$\Delta I=2@f$
     */
    double getReA2_Kd() const
    {
        return ReA2_Kd;
    }

    /**
     * @return the isospin breaking contribution in @f$K^0\to\pi\pi@f$
     */
    double getOmega_eta_etap() const
    {
        return Omega_eta_etap;
    }

    /**
     * @return the experimental value for the branching ratio of @f$K^+\to\pi^0e^+\nu@f$
     */
    double getBr_Kp_P0enu() const
    {
        return Br_Kp_P0enu;
    }

    /**
     * @return the experimental value for the branching ratio of @f$K^+\to\mu^+\nu@f$
     */
    double getBr_Kp_munu() const
    {
        return Br_Kp_munu;
    }

    /**
     * @return the experimental value for the branching ratio of @f$B\to X_c e\nu@f$
     */
    double getBr_B_Xcenu() const
    {
        return Br_B_Xcenu;
    }

    /**
     * @return the long-distance correction to the charm contribution of @f$K^+\to\pi^+\nu\bar{\nu}@f$
     *
     * References:
     * [<A HREF="http://inspirehep.net/record/678222?ln=en" target="blank">Isidori et al.(2005)</A>],
     * [<A HREF="http://inspirehep.net/record/712083?ln=en" target="blank">Buras et al.(2006)</A>]
     */
    double getDeltaP_cu() const
    {
        return DeltaP_cu;
    }

    /**
     * @return the isospin breaking corrections between
     * @f$K_L\to\pi^0\nu\bar{\nu}@f$ and @f$K^+\to\pi^0 e^+\nu@f$
     */
    double getIB_Kl() const
    {
        return IB_Kl;
    }

    /**
     * @return the isospin breaking corrections between
     * @f$K^+\to\pi^+ \nu\bar{\nu}@f$ and @f$K^+\to\pi^0 e^+\nu@f$
     */
    double getIB_Kp() const
    {
        return IB_Kp;
    }
    
    /**
     * @return the real part of the constant term of the charm loop long distance contribution @f$h_0@f$ to @f$M\to V@f$
     */
    double getReh_0() const
    {
        return reh_0;
    }
    
    /**
     * @return the real part of the constant term of the charm loop long distance contribution @f$h_+@f$ to @f$M\to V@f$
     */
    double getReh_p() const
    {
        return reh_p;
    }
    
    /**
     * @return the real part of the constant term of the charm loop long distance contribution @f$h_-@f$ to @f$M\to V@f$
     */
    double getReh_m() const
    {
        return reh_m;
    }
    
    /**
     * @return the imaginary part of the constant term of the charm loop long distance contribution @f$h_0@f$ to @f$M\to V@f$
     */
    double getImh_0() const
    {
        return imh_0;
    }
    
    /**
     * @return the imaginary part of the constant term of the charm loop long distance contribution @f$h_+@f$ to @f$M\to V@f$
     */
    double getImh_p() const
    {
        return imh_p;
    }
    
    /**
     * @return the imaginary part of the constant term of the charm loop long distance contribution @f$h_-@f$ to @f$M\to V@f$
     */
    double getImh_m() const
    {
        return imh_m;
    }
    
    /**
     * @return the real part of the linear term of the charm loop long distance contribution @f$h_0@f$ to @f$M\to V@f$
     */
    double getReh_0_1() const
    {
        return reh_0_1;
    }
    
    /**
     * @return the real part of the linear term of the charm loop long distance contribution @f$h_+@f$ to @f$M\to V@f$
     */
    double getReh_p_1() const
    {
        return reh_p_1;
    }
    
    /**
     * @return the real part of the linear term of the charm loop long distance contribution @f$h_-@f$ to @f$M\to V@f$
     */
    double getReh_m_1() const
    {
        return reh_m_1;
    }
    
    /**
     * @return the imaginary part of the linear term of the charm loop long distance contribution @f$h_0@f$ to @f$M\to V@f$
     */
    double getImh_0_1() const
    {
        return imh_0_1;
    }
    
    /**
     * @return the imaginary part of the linear term of the charm loop long distance contribution @f$h_+@f$ to @f$M\to V@f$
     */
    double getImh_p_1() const
    {
        return imh_p_1;
    }
    
    /**
     * @return the imaginary part of the linear term of the charm loop long distance contribution @f$h_-@f$ to @f$M\to V@f$
     */
    double getImh_m_1() const
    {
        return imh_m_1;
    }
    
    /**
     * @return the real part of the constant term of the charm loop long distance contribution @f$h_0@f$ to @f$M\to P@f$
     */
    double getReh_0_MP() const
    {
        return reh_0_MP;
    }
    
    /**
     * @return the imaginary part of the constant term of the charm loop long distance contribution @f$h_0@f$ to @f$M\to P@f$
     */
    double getImh_0_MP() const
    {
        return imh_0_MP;
    }
    
    /**
     * @return the real part of the linear term of the charm loop long distance contribution @f$h_0@f$ to @f$M\to P@f$
     */
    double getReh_0_1_MP() const
    {
        return reh_0_1_MP;
    }
    
    /**
     * @return the imaginary part of the linear term of the charm loop long distance contribution @f$h_0@f$ to @f$M\to P@f$
     */
    double getImh_0_1_MP() const
    {
        return imh_0_1_MP;
    }
    
    /**
     * @return the constant term of the charm loop long distance contribution @f$h_0@f$ to @f$M\to V@f$
     */
    gslpp::complex geth_0() const
    {
        return myh_0;
    }

    /**
     * @return the constant term of the charm loop long distance contribution @f$h_+@f$ to @f$M\to V@f$
     */
    gslpp::complex geth_p() const
    {
        return myh_p;
    }
    
    /**
     * @return the constant term of the charm loop long distance contribution @f$h_-@f$ to @f$M\to V@f$
     */
    gslpp::complex geth_m() const
    {
        return myh_m;
    }
    
    /**
     * @return the linear term of the charm loop long distance contribution @f$h_0@f$ to @f$M\to V@f$
     */
    gslpp::complex geth_0_1() const
    {
        return myh_0_1;
    }

    /**
     * @return the linear term of the charm loop long distance contribution @f$h_+@f$ to @f$M\to V@f$
     */
    gslpp::complex geth_p_1() const
    {
        return myh_p_1;
    }
    
    /**
     * @return the linear term of the charm loop long distance contribution @f$h_-@f$ to @f$M\to V@f$
     */
    gslpp::complex geth_m_1() const
    {
        return myh_m_1;
    }
    
    /**
     * @return the quadratic term of the charm loop long distance contribution @f$h_0@f$ to @f$M\to V@f$
     */
    gslpp::complex geth_0_2() const
    {
        return myh_0_2;
    }

    /**
     * @return the quadratic term of the charm loop long distance contribution @f$h_+@f$ to @f$M\to V@f$
     */
    gslpp::complex geth_p_2() const
    {
        return myh_p_2;
    }
    
    /**
     * @return the quadratic term of the charm loop long distance contribution @f$h_-@f$ to @f$M\to V@f$
     */
    gslpp::complex geth_m_2() const
    {
        return myh_m_2;
    }
    
    /**
     * @return the constant term of the charm loop long distance contribution @f$h_0@f$ to @f$M\to P@f$
     */
    gslpp::complex geth_0_MP() const
    {
        return reh_0_MP + gslpp::complex::i() * imh_0_MP;
    }
    
    /**
     * @return the linear term of the charm loop long distance contribution @f$h_0@f$ to @f$M\to P@f$
     */
    gslpp::complex geth_0_1_MP() const
    {
        return reh_0_1_MP + gslpp::complex::i() * imh_0_1_MP;
    }
    
    /**
     * @return the LCSR fit parameter @f$a_0^V@f$ for @f$B\to K^*@f$ from arXiv:1503.05534v1
     */
    double geta_0V() const
    {
        return a_0V;
    }
    
    /**
     * @return the LCSR fit parameter @f$a_1^V@f$ for @f$B\to K^*@f$ from arXiv:1503.05534v1
     */
    double geta_1V() const
    {
        return a_1V;
    }
    
    /**
     * @return the LCSR fit parameter @f$a_2^V@f$ for @f$B\to K^*@f$ from arXiv:1503.05534v1
     */
    double geta_2V() const
    {
        return a_2V;
    }
    
    /**
     * @return the LCSR fit parameter @f$\Delta m^V@f$ for @f$B\to K^*@f$ from arXiv:1503.05534v1
     */
    double getMRV() const
    {
        return MRV;
    }
    
    /**
     * @return the LCSR fit parameter @f$a_0^{A0}@f$ for @f$B\to K^*@f$ from arXiv:1503.05534v1
     */
    double geta_0A0() const
    {
        return a_0A0;
    }
    
    /**
     * @return the LCSR fit parameter @f$a_1^{A0}@f$ for @f$B\to K^*@f$ from arXiv:1503.05534v1
     */
    double geta_1A0() const
    {
        return a_1A0;
    }
    
    /**
     * @return the LCSR fit parameter @f$a_2^{A0}@f$ for @f$B\to K^*@f$ from arXiv:1503.05534v1
     */
    double geta_2A0() const
    {
        return a_2A0;
    }
    
    /**
     * @return the LCSR fit parameter @f$\Delta m^{A0}@f$ for @f$B\to K^*@f$ from arXiv:1503.05534v1
     */
    double getMRA0() const
    {
        return MRA0;
    }
    
    /**
     * @return the LCSR fit parameter @f$a_0^{A1}@f$ for @f$B\to K^*@f$ from arXiv:1503.05534v1
     */
    double geta_0A1() const
    {
        return a_0A1;
    }
    
    /**
     * @return the LCSR fit parameter @f$a_1^{A1}@f$ for @f$B\to K^*@f$ from arXiv:1503.05534v1
     */
    double geta_1A1() const
    {
        return a_1A1;
    }
    
    /**
     * @return the LCSR fit parameter @f$a_2^{A1}@f$ for @f$B\to K^*@f$ from arXiv:1503.05534v1
     */
    double geta_2A1() const
    {
        return a_2A1;
    }
    
    /**
     * @return the LCSR fit parameter @f$\Delta m^{A1}@f$ for @f$B\to K^*@f$ from arXiv:1503.05534v1
     */
    double getMRA1() const
    {
        return MRA1;
    }
    
    /**
     * @return the LCSR fit parameter @f$a_0^{A12}@f$ for @f$B\to K^*@f$ from arXiv:1503.05534v1
     */
    double geta_0A12() const
    {
        return a_0A12;
    }
    
    /**
     * @return the LCSR fit parameter @f$a_1^{A12}@f$ for @f$B\to K^*@f$ from arXiv:1503.05534v1
     */
    double geta_1A12() const
    {
        return a_1A12;
    }
    
    /**
     * @return the LCSR fit parameter @f$a_2^{A12}@f$ for @f$B\to K^*@f$ from arXiv:1503.05534v1
     */
    double geta_2A12() const
    {
        return a_2A12;
    }
    
    /**
     * @return the LCSR fit parameter @f$\Delta m^{A12}@f$ for @f$B\to K^*@f$ from arXiv:1503.05534v1
     */
    double getMRA12() const
    {
        return MRA12;
    }
    
    /**
     * @return the LCSR fit parameter @f$a_0^{T1}@f$ for @f$B\to K^*@f$ from arXiv:1503.05534v1
     */
    double geta_0T1() const
    {
        return a_0T1;
    }
    
    /**
     * @return the LCSR fit parameter @f$a_1^{T1}@f$ for @f$B\to K^*@f$ from arXiv:1503.05534v1
     */
    double geta_1T1() const
    {
        return a_1T1;
    }
    
    /**
     * @return the LCSR fit parameter @f$a_2^{T1}@f$ for @f$B\to K^*@f$ from arXiv:1503.05534v1
     */
    double geta_2T1() const
    {
        return a_2T1;
    }
    
    /**
     * @return the LCSR fit parameter @f$\Delta m^{T1}@f$ for @f$B\to K^*@f$ from arXiv:1503.05534v1
     */
    double getMRT1() const
    {
        return MRT1;
    }
    
    /**
     * @return the LCSR fit parameter @f$a_0^{T2}@f$ for @f$B\to K^*@f$ from arXiv:1503.05534v1
     */
    double geta_0T2() const
    {
        return a_0T2;
    }
    
    /**
     * @return the LCSR fit parameter @f$a_1^{T2}@f$ for @f$B\to K^*@f$ from arXiv:1503.05534v1
     */
    double geta_1T2() const
    {
        return a_1T2;
    }
    
    /**
     * @return the LCSR fit parameter @f$a_2^{T2}@f$ for @f$B\to K^*@f$ from arXiv:1503.05534v1
     */
    double geta_2T2() const
    {
        return a_2T2;
    }
    
    /**
     * @return the LCSR fit parameter @f$\Delta m^{T2}@f$ for @f$B\to K^*@f$ from arXiv:1503.05534v1
     */
    double getMRT2() const
    {
        return MRT2;
    }
    
    /**
     * @return the LCSR fit parameter @f$a_0^{T23}@f$ for @f$B\to K^*@f$ from arXiv:1503.05534v1
     */
    double geta_0T23() const
    {
        return a_0T23;
    }
    
    /**
     * @return the LCSR fit parameter @f$a_1^{T23}@f$ for @f$B\to K^*@f$ from arXiv:1503.05534v1
     */
    double geta_1T23() const
    {
        return a_1T23;
    }
    
    /**
     * @return the LCSR fit parameter @f$a_2^{T23}@f$ for @f$B\to K^*@f$ from arXiv:1503.05534v1
     */
    double geta_2T23() const
    {
        return a_2T23;
    }
    
    /**
     * @return the LCSR fit parameter @f$\Delta m^{T23}@f$ for @f$B\to K^*@f$ from arXiv:1503.05534v1
     */
    double getMRT23() const
    {
        return MRT23;
    }

    /**
     * @return the LCSR fit parameter @f$a_0^V@f$ for @f$B\to\phi@f$ from arXiv:1503.05534v1
     */
    double geta_0Vphi() const
    {
        return a_0Vphi;
    }
    
    /**
     * @return the LCSR fit parameter @f$a_1^V@f$ for @f$B\to\phi@f$ from arXiv:1503.05534v1
     */
    double geta_1Vphi() const
    {
        return a_1Vphi;
    }
    
    /**
     * @return the LCSR fit parameter @f$a_2^V@f$ for @f$B\to\phi@f$ from arXiv:1503.05534v1
     */
    double geta_2Vphi() const
    {
        return a_2Vphi;
    }
    
    /**
     * @return the LCSR fit parameter @f$\Delta m^V@f$ for @f$B\to\phi@f$ from arXiv:1503.05534v1
     */
    double getMRVphi() const
    {
        return MRVphi;
    }
    
    /**
     * @return the LCSR fit parameter @f$a_0^{A0}@f$ for @f$B\to\phi@f$ from arXiv:1503.05534v1
     */
    double geta_0A0phi() const
    {
        return a_0A0phi;
    }
    
    /**
     * @return the LCSR fit parameter @f$a_1^{A0}@f$ for @f$B\to\phi@f$ from arXiv:1503.05534v1
     */
    double geta_1A0phi() const
    {
        return a_1A0phi;
    }
    
    /**
     * @return the LCSR fit parameter @f$a_2^{A0}@f$ for @f$B\to\phi@f$ from arXiv:1503.05534v1
     */
    double geta_2A0phi() const
    {
        return a_2A0phi;
    }
    
    /**
     * @return the LCSR fit parameter @f$\Delta m^{A0}@f$ for @f$B\to\phi@f$ from arXiv:1503.05534v1
     */
    double getMRA0phi() const
    {
        return MRA0phi;
    }
    
    /**
     * @return the LCSR fit parameter @f$a_0^{A1}@f$ for @f$B\to\phi@f$ from arXiv:1503.05534v1
     */
    double geta_0A1phi() const
    {
        return a_0A1phi;
    }
    
    /**
     * @return the LCSR fit parameter @f$a_1^{A1}@f$ for @f$B\to\phi@f$ from arXiv:1503.05534v1
     */
    double geta_1A1phi() const
    {
        return a_1A1phi;
    }
    
    /**
     * @return the LCSR fit parameter @f$a_2^{A1}@f$ for @f$B\to\phi@f$ from arXiv:1503.05534v1
     */
    double geta_2A1phi() const
    {
        return a_2A1phi;
    }
    
    /**
     * @return the LCSR fit parameter @f$\Delta m^{A1}@f$ for @f$B\to\phi@f$ from arXiv:1503.05534v1
     */
    double getMRA1phi() const
    {
        return MRA1phi;
    }
    
    /**
     * @return the LCSR fit parameter @f$a_0^{A12}@f$ for @f$B\to\phi@f$ from arXiv:1503.05534v1
     */
    double geta_0A12phi() const
    {
        return a_0A12phi;
    }
    
    /**
     * @return the LCSR fit parameter @f$a_1^{A12}@f$ for @f$B\to\phi@f$ from arXiv:1503.05534v1
     */
    double geta_1A12phi() const
    {
        return a_1A12phi;
    }
    
    /**
     * @return the LCSR fit parameter @f$a_2^{A12}@f$ for @f$B\to\phi@f$ from arXiv:1503.05534v1
     */
    double geta_2A12phi() const
    {
        return a_2A12phi;
    }
    
    /**
     * @return the LCSR fit parameter @f$\Delta m^{A12}@f$ for @f$B\to\phi@f$ from arXiv:1503.05534v1
     */
    double getMRA12phi() const
    {
        return MRA12phi;
    }
    
    /**
     * @return the LCSR fit parameter @f$a_0^{T1}@f$ for @f$B\to\phi@f$ from arXiv:1503.05534v1
     */
    double geta_0T1phi() const
    {
        return a_0T1phi;
    }
    
    /**
     * @return the LCSR fit parameter @f$a_1^{T1}@f$ for @f$B\to\phi@f$ from arXiv:1503.05534v1
     */
    double geta_1T1phi() const
    {
        return a_1T1phi;
    }
    
    /**
     * @return the LCSR fit parameter @f$a_2^{T1}@f$ for @f$B\to\phi@f$ from arXiv:1503.05534v1
     */
    double geta_2T1phi() const
    {
        return a_2T1phi;
    }
    
    /**
     * @return the LCSR fit parameter @f$\Delta m^{T1}@f$ for @f$B\to\phi@f$ from arXiv:1503.05534v1
     */
    double getMRT1phi() const
    {
        return MRT1phi;
    }
    
    /**
     * @return the LCSR fit parameter @f$a_0^{T2}@f$ for @f$B\to\phi@f$ from arXiv:1503.05534v1
     */
    double geta_0T2phi() const
    {
        return a_0T2phi;
    }
    
    /**
     * @return the LCSR fit parameter @f$a_1^{T2}@f$ for @f$B\to\phi@f$ from arXiv:1503.05534v1
     */
    double geta_1T2phi() const
    {
        return a_1T2phi;
    }
    
    /**
     * @return the LCSR fit parameter @f$a_2^{T2}@f$ for @f$B\to\phi@f$ from arXiv:1503.05534v1
     */
    double geta_2T2phi() const
    {
        return a_2T2phi;
    }
    
    /**
     * @return the LCSR fit parameter @f$\Delta m^{T2}@f$ for @f$B\to\phi@f$ from arXiv:1503.05534v1
     */
    double getMRT2phi() const
    {
        return MRT2phi;
    }
    
    /**
     * @return the LCSR fit parameter @f$a_0^{T23}@f$ for @f$B\to\phi@f$ from arXiv:1503.05534v1
     */
    double geta_0T23phi() const
    {
        return a_0T23phi;
    }
    
    /**
     * @return the LCSR fit parameter @f$a_1^{T23}@f$ for @f$B\to\phi@f$ from arXiv:1503.05534v1
     */
    double geta_1T23phi() const
    {
        return a_1T23phi;
    }
    
    /**
     * @return the LCSR fit parameter @f$a_2^{T23}@f$ for @f$B\to\phi@f$ from arXiv:1503.05534v1
     */
    double geta_2T23phi() const
    {
        return a_2T23phi;
    }
    
    /**
     * @return the LCSR fit parameter @f$\Delta m^{T23}@f$ for @f$B\to\phi@f$ from arXiv:1503.05534v1
     */
    double getMRT23phi() const
    {
        return MRT23phi;
    }
    
    /**
     * @return the LCSR fit parameter @f$r_1^{f_+}@f$ for @f$B\to K@f$ from arXiv:hep-ph/0406232v1
     */
    double getr_1_fplus() const
    {
        return r_1_fplus;
    }

    /**
     * @return the LCSR fit parameter @f$r_2^{f_+}@f$ for @f$B\to K@f$ from arXiv:hep-ph/0406232v1
     */
    double getr_2_fplus() const
    {
        return r_2_fplus;
    }
    
    /**
     * @return the LCSR fit parameter @f$m_{fit}^{2,f_+}@f$ for @f$B\to K@f$ from arXiv:hep-ph/0406232v1
     */
    double getm_fit2_fplus() const
    {
        return m_fit2_fplus;
    }
    
    /**
     * @return the LCSR fit parameter @f$r_1^{f_T}@f$ for @f$B\to K@f$ from arXiv:hep-ph/0406232v1
     */
    double getr_1_fT() const
    {
        return r_1_fT;
    }

    /**
     * @return the LCSR fit parameter @f$r_2^{f_T}@f$ for @f$B\to K@f$ from arXiv:hep-ph/0406232v1
     */
    double getr_2_fT() const
    {
        return r_2_fT;
    }
    
    /**
     * @return the LCSR fit parameter @f$m_{fit}^{2,f_T}@f$ for @f$B\to K@f$ from arXiv:hep-ph/0406232v1
     */
    double getm_fit2_fT() const
    {
        return m_fit2_fT;
    }
    
    /**
     * @return the LCSR fit parameter @f$r_2^{f_0}@f$ for @f$B\to K@f$ from arXiv:hep-ph/0406232v1
     */
    double getr_2_f0() const
    {
        return r_2_f0;
    }
    
    /**
     * @return the LCSR fit parameter @f$m_{fit}^{2,f_0}@f$ for @f$B\to K@f$ from arXiv:hep-ph/0406232v1
     */
    double getm_fit2_f0() const
    {
        return m_fit2_f0;
    }
    
    /**
     * @return the experimental energy cutoff for @f$b\to s \gamma@f$
     */
    double getbsgamma_E0() const
    {
        return bsgamma_E0;
    }
    
    /**
     * @return the Benzke, Lee, Neubert, Paz non perturbative correction for @f$b\to s \gamma@f$ from arXiv:1003.5012.
     */
    double getBLNPcorr() const
    {
        return BLNPcorr;
    }
    
    /**
     * @return the kinetic scale used in arXiv:1411.6560.
     */
    double getGambino_mukin() const
    {
        return Gambino_mukin;
    }
    
    /**
     * @return the fit value for the branching ratio of @f$B\to X_c e\nu@f$ computed as in arXiv:1411.6560, but with @f$\mu_c=2GeV@f$.
     */
    double getGambino_BRsem() const
    {
        return Gambino_BRsem;
    }
    
    /**
     * @return the fit value for the kinetic b mass @f$M_b^{\rm kin}(\mu^{\rm kin})@f$ computed as in arXiv:1411.6560, but with @f$\mu_c=2GeV@f$.
     */
    double getGambino_Mbkin() const
    {
        return Gambino_Mbkin;
    }
    
    /**
     * @return the fit value for the MSbar mass @f$M_c(\mu_c)@f$ computed as in arXiv:1411.6560, but with @f$\mu_c=2GeV@f$.
     */
    double getGambino_Mcatmuc() const
    {
        return Gambino_Mcatmuc;
    }
    
    /**
     * @return the fit value for @f$\mu_{\pi}^2@f$ computed as in arXiv:1411.6560, but with @f$\mu_c=2GeV@f$.
     */
    double getGambino_mupi2() const
    {
        return Gambino_mupi2;
    }
    
    /**
     * @return the fit value for @f$\rho_D^3@f$ computed as in arXiv:1411.6560, but with @f$\mu_c=2GeV@f$.
     */
    double getGambino_rhoD3() const
    {
        return Gambino_rhoD3;
    }
    
    /**
     * @return the fit value for @f$\mu_G^2@f$ computed as in arXiv:1411.6560, but with @f$\mu_c=2GeV@f$.
     */
    double getGambino_muG2() const
    {
        return Gambino_muG2;
    }
    
    /**
     * @return the fit value for@f$\rho_{LS}^3@f$ computed as in arXiv:1411.6560, but with @f$\mu_c=2GeV@f$.
     */
    double getGambino_rhoLS3() const
    {
        return Gambino_rhoLS3;
    }
    
    /**
     * @return the decay constant of a transversely polarized @f$K^*@f$ meson at 1 GeV
     */
    double getFKstarp() const
    {
        return FKstarp;
    }

    
    ////////////////////////////////////////////////////////////////////////

    /**
     * @brief For accessing the active flavour threshold scales.
     * @param[in] i the index referring to active flavour thresholds.
     * @return the threshold scale: 1.0E10 (i = 0), \f$\mu_t\f$ (i = 1),
     * \f$\mu_b\f$ (i = 2), \f$\mu_c\f$ (i = 3) and 0. (default)
     */
    double Thresholds(const int i) const;

    /**
     * @brief The active flavour threshold above the scale \f$\mu\f$
     * as defined in %QCD::Thresholds().
     * @param[in] mu a scale \f$\mu\f$ in GeV
     * @return the higher active flavour threshold
     */
    double AboveTh(const double mu) const;

    /**
     * @brief The active flavour threshold below the scale \f$\mu\f$
     * as defined in QCD::Thresholds().
     * @param[in] mu a scale \f$\mu\f$ in GeV
     * @return the lower active flavour threshold
     */
    double BelowTh(const double mu) const;

    /**
     * @brief The number of active flavour at scale @f$\mu@f$.
     * @param[in] mu a scale @f$\mu@f$ in GeV
     * @return active N_f
     */
    double Nf(const double mu) const;

    ////////////////////////////////////////////////////////////////////////

    /**
     * @brief The \f$\beta_0(n_f)\f$ coefficient for a certain number of flavours \f$n_f\f$.
     * @param[in] nf the number of active flavours \f$n_f\f$
     * @return @f$\beta_0(n_f)@f$
     */
    double Beta0(const double nf) const;

    /**
     * @brief The \f$\beta_1(n_f)\f$ coefficient for a certain number of flavours \f$n_f\f$.
     * @param[in] nf the number of active flavours \f$n_f\f$
     * @return @f$\beta_1(n_f)@f$
     */
    double Beta1(const double nf) const;

    /**
     * @brief The \f$\beta_2(n_f)\f$ coefficient for a certain number of flavours \f$n_f\f$.
     * @param[in] nf the number of active flavours \f$n_f\f$
     * @return @f$\beta_2(n_f)@f$
     */
    double Beta2(const double nf) const;

    /**
     * @brief Computes the running strong coupling @f$\alpha_s(\mu)@f$ from @f$\alpha_s(\mu_i)@f$
     * in the @f$\overline{\mathrm{MS}}@f$ scheme, where it is forbidden to across
     * a flavour threshold in the RG running from @f$\mu_i@f$ to @f$\mu@f$.
     * @param[in] mu a scale @f$\mu@f$ in GeV
     * @param[in] alsi the initial value for the coupling at the scale given below
     * @param[in] mu_i the initial scale @f$\mu_i@f$ in GeV
     * @param[in] order LO, NLO or FULLNLO in the @f$\alpha_s@f$ expansion defined in OrderScheme
     * @return the strong coupling constant @f$\alpha_s(\mu)@f$ in the
     * @f$\overline{\mathrm{MS}}@f$ scheme
     */
    double AlsWithInit(const double mu, const double alsi, const double mu_i,
            const orders order) const;

    /**
     * @brief Computes the running strong coupling @f$\alpha_s(\mu)@f$ in the
     * @f$\overline{\mathrm{MS}}@f$ scheme with the use of @f$\Lambda_{\rm QCD}@f$.
     * @param[in] mu A scale @f$\mu@f$ in GeV
     * @param[in] order LO, NLO, FULLNLO, NNLO or FULLNNLO in the @f$\alpha_s@f$ expansion defined in OrderScheme
     * @return the strong coupling constant @f$\alpha_s(\mu)@f$ in the
     * @f$\overline{\mathrm{MS}}@f$ scheme
     */
    double AlsWithLambda(const double mu, const orders order) const;

    /**
     * @brief Computes the running strong coupling @f$\alpha_s(\mu)@f$ in the
     * @f$\overline{\mathrm{MS}}@f$ scheme. In the cases of LO, NLO and FULLNNLO,
     * the coupling is computed with AlsWithInit(). On the other hand, in the
     * cases of NNLO and FULLNNLO, the coupling is computed with AlsWithLambda().
     * @param[in] mu a scale @f$\mu@f$ in GeV
     * @param[in] order LO, NLO, FULLNLO, NNLO or FULLNNLO in the @f$\alpha_s@f$ expansion defined in OrderScheme
     * @return the strong coupling constant @f$\alpha_s(\mu)@f$ in the
     * @f$\overline{\mathrm{MS}}@f$ scheme
     */
    double Als(const double mu, const orders order = FULLNLO) const;

    /**
     * @brief Computes @f$\ln\Lambda_\mathrm{QCD}@f$ with nf flavours in GeV.
     * @param[in] nf the number of active flavours \f$n_f\f$
     * @param[in] order LO, NLO, FULLNLO, NNLO or FULLNNLO in the @f$\alpha_s@f$ expansion defined in OrderScheme
     * @return @f$\ln\Lambda_\mathrm{QCD}@f$ with nf flavours in GeV
     */
    double logLambda(const double nf, orders order) const;

    /**
     * @brief The value of \f$\alpha_s^{\mathrm{FULLNLO}}\f$ at any scale \f$\mu\f$ with the number of flavours
     * \f$n_f = 4\f$.
     * @param[in] mu the scale at which \f$\alpha_s\f$ has to be computed
     * @return \f$\alpha_s^{\mathrm{FULLNLO}}(\mu)\f$ with \f$n_f = 4\f$
     *
     * @attention Temporary function waiting for the implementation of NNLO etact.
     */
    double Als4(const double mu) const;
    
    virtual double alphaMz() const = 0;

    /**
     * @brief The value of \f$\frac{\alpha_s^{\mathrm{FULLNLO}}}{4\pi}\f$ at any scale \f$\mu\f$ with the number of flavours
     * \f$n_f = 4\f$ and full EW corrections.
     * @param[in] mu the scale at which \f$\alpha_s\f$ has to be computed
     * @return \f$\alpha_s^{\mathrm{FULLNLO}}(\mu)\f$ with \f$n_f = 4\5$
     */
    double Alstilde5(const double mu) const;

    /**
     * @brief The running of a mass with the number of flavours \f$n_f = 4\f$.
     * @param[in] mu_f the final scale \f$\mu_f\f$ to which the mass is run
     * @param[in] mu_i the initial scale \f$\mu_i\f$ from which the mass is run
     * @param[in] m the mass at the scale \f$\mu_i\f$
     * @return the mass at the scale \f$\mu_f\f$ with \f$n_f = 4\f$
     *
     * @attention Temporary function waiting for the implementation of NNLO etact.
     */
    double Mrun4(const double mu_f, const double mu_i, const double m) const;

    ////////////////////////////////////////////////////////////////////////

    /**
     * @brief The @f$\gamma_0@f$ coefficient used to compute the running of a mass.
     * @param[in] nf the number of active flavours \f$n_f\f$
     * @return the @f$\gamma_0@f$ coefficient.
     */
    double Gamma0(const double nf) const;

    /**
     * @brief The @f$\gamma_1@f$ coefficient used to compute the running of a mass.
     * @param[in] nf the number of active flavours \f$n_f\f$
     * @return the @f$\gamma_1@f$ coefficient
     */
    double Gamma1(const double nf) const;

    /**
     * @brief The @f$\gamma_2@f$ coefficient used to compute the running of a mass.
     * @param[in] nf The number of active flavours \f$n_f\f$
     * @return the @f$\gamma_2@f$ coefficient
     */
    double Gamma2(const double nf) const;

    /**
     * @brief Computes a running quark mass @f$m(\mu)@f$ from @f$m(m)@f$.
     * @param[in] mu a scale @f$\mu@f$ in GeV
     * @param[in] m the @f$\overline{\mathrm{MS}}@f$ mass @f$m(m)@f$ in GeV
     * @param[in] order LO, NLO, FULLNLO, NNLO or FULLNNLO in the @f$\alpha_s@f$ expansion defined in OrderScheme
     * @return the running quark mass @f$m(\mu)@f$ in GeV
     */
    double Mrun(const double mu, const double m, const orders order = FULLNLO) const;

    /**
     * @brief Runs a quark mass from @f$\mu_i@f$ to @f$\mu_f@f$.
     * @param[in] mu_f a scale @f$\mu_f@f$ in GeV
     * @param[in] mu_i a scale @f$\mu_i@f$ in GeV
     * @param[in] m the @f$\overline{\mathrm{MS}}@f$ mass @f$m(\mu_i)@f$ in GeV
     * @param[in] order LO, NLO, FULLNLO, NNLO or FULLNNLO in the @f$\alpha_s@f$ expansion defined in OrderScheme
     * @return the running quark mass @f$m(\mu_f)@f$ in GeV
     */
    double Mrun(const double mu_f, const double mu_i, const double m,
            const orders order = FULLNLO) const;

    ////////////////////////////////////////////////////////////////////////

    /**
     * @brief Converts the @f$\overline{\mathrm{MS}}@f$ mass @f$m(m)@f$ to the pole mass.
     * @param[in] mbar the @f$\overline{\mathrm{MS}}@f$ mass @f$m(m)@f$ in GeV
     * @param[in] order LO, NLO, FULLNLO, NNLO or FULLNNLO in the @f$\alpha_s@f$ expansion defined in OrderScheme
     * @return the pole mass in GeV
     *
     * @attention Can only be used for conversion of mass of the top and bottom quarks.
     */
    double Mbar2Mp(const double mbar, const orders order = FULLNLO) const;

    /**
     * @brief Converts a quark pole mass to the corresponding @f$\overline{\mathrm{MS}}@f$ mass @f$m(m)@f$.
     * @param[in] mp the pole mass of the bottom or top quark in GeV
     * @param[in] order LO, NLO, FULLNLO, NNLO or FULLNNLO in the @f$\alpha_s@f$ expansion defined in OrderScheme
     * @return the @f$\overline{\mathrm{MS}}@f$ mass @f$m(m)@f$ in GeV
     */
    double Mp2Mbar(const double mp, const orders order = FULLNLO) const;

    /**
     * @brief Converts a quark mass from the @f$\overline{\mathrm{MS}}@f$ scheme to
     * the @f$\overline{\mathrm{DR}}@f$ scheme.
     * @param[in] MSscale the scale at which the @f$\overline{\mathrm{MS}}@f$ mass is defined
     * @param[in] MSbar the @f$\overline{\mathrm{MS}}@f$ mass @f$m(m)@f$ in GeV
     * @return the @f$\overline{\mathrm{DR}}@f$ mass @f$m(m)@f$ in GeV
     */
    double MS2DRqmass(const double MSscale, const double MSbar) const;

    /**
     * @brief Converts a quark mass from the @f$\overline{\mathrm{MS}}@f$ scheme to
     * the @f$\overline{\mathrm{DR}}@f$ scheme.
     * @param[in] MSbar The @f$\overline{\mathrm{MS}}@f$ mass @f$m(m)@f$ in GeV
     * @return The @f$\overline{\mathrm{DR}}@f$ mass @f$m(m)@f$ in GeV
     */
    double MS2DRqmass(const double MSbar) const;

    ////////////////////////////////////////////////////////////////////////

protected:

    /**
     * @brief A method to set the value of a parameter of %QCD.
     * @param[in] name name of a model parameter
     * @param[in] value the value to be assigned to the parameter specified by name
     */
    virtual void setParameter(const std::string name, const double& value);

    bool computemt; ///< Switch for computing the \f$\overline{\mathrm{MS}}\f$ mass of the top quark.
    bool requireYu; ///< Switch for generating the Yukawa couplings to the up-type quarks.
    bool requireYd; ///< Switch for generating the Yukawa couplings to the down-type quarks.

    // model parameters
    double AlsM; ///< The strong coupling constant at the mass scale MAls, \f$\alpha_s(M_{\alpha_s})\f$. 
    double MAls; ///< The mass scale in GeV at which the strong coupling measurement is provided. 
    double mtpole; ///< The pole mass of the top quark.
    double mut; ///< The threshold between six- and five-flavour theory in GeV. 
    double mub; ///< The threshold between five- and four-flavour theory in GeV. 
    double muc; ///< The threshold between four- and three-flavour theory in GeV. 
    double ReA0_Kd;
    double ReA2_Kd;
    double Omega_eta_etap;
    double Br_Kp_P0enu;
    double IB_Kl;
    double IB_Kp;
    double DeltaP_cu;
    double Br_Kp_munu;
    double Br_B_Xcenu;
    double BBsoBBd; ///< The ratio \f$ B_{B_s}/B_{B_d} \f$ necessary to compute \f$ B_{B_s} \f$. 
    double FBsoFBd; ///< The ratio \f$ F_{B_s}/F_{B_d} \f$ necessary to compute \f$ F_{B_s} \f$. 
    double reh_0, reh_p, reh_m, imh_0, imh_p, imh_m;
    double reh_0_1, reh_p_1, reh_m_1, imh_0_1, imh_p_1, imh_m_1;
    double reh_0_2, reh_p_2, reh_m_2, imh_0_2, imh_p_2, imh_m_2;
    double reh_0_MP, imh_0_MP, reh_0_1_MP, imh_0_1_MP;
    double a_0V, a_1V, a_2V, MRV;
    double a_0A0, a_1A0, a_2A0, MRA0;
    double a_0A1, a_1A1, a_2A1, MRA1;
    double a_0A12, a_1A12, a_2A12, MRA12;
    double a_0T1, a_1T1, a_2T1, MRT1;
    double a_0T2, a_1T2, a_2T2, MRT2;
    double a_0T23, a_1T23, a_2T23, MRT23;
    double a_0Vphi, a_1Vphi, a_2Vphi, MRVphi;
    double a_0A0phi, a_1A0phi, a_2A0phi, MRA0phi;
    double a_0A1phi, a_1A1phi, a_2A1phi, MRA1phi;
    double a_0A12phi, a_1A12phi, a_2A12phi, MRA12phi;
    double a_0T1phi, a_1T1phi, a_2T1phi, MRT1phi;
    double a_0T2phi, a_1T2phi, a_2T2phi, MRT2phi;
    double a_0T23phi, a_1T23phi, a_2T23phi, MRT23phi;
    double r_1_fplus, r_2_fplus, m_fit2_fplus;
    double r_1_fT, r_2_fT, m_R_fT, m_fit2_fT;
    double r_2_f0, m_fit2_f0;
    double bsgamma_E0, BLNPcorr;
    double Gambino_mukin, Gambino_BRsem, Gambino_Mbkin, Gambino_Mcatmuc, Gambino_mupi2, Gambino_rhoD3, Gambino_muG2, Gambino_rhoLS3;
    double FKstarp; //matrix element of tensor current for transverse polarization at 1 GeV
    
    //double r_2A0, r_2T1, r_2T2, r_2A0phi, r_2T1phi, r_2T2phi removed because they are fixed by form factors relations
    double Nc; ///< The number of colours.
    Particle quarks[6]; ///< The vector of all SM quarks.
    Meson mesons[MESON_END]; ///< The vector of defined mesons.
private:
    
    double CF; ///< The Casimir factor in the \f$SU(N_c)\f$ gauge theory.
    BParameter BBs; ///< The bag parameters for \f$\Delta b=2\f$ processes for the \f$B_s\f$ meson system.
    BParameter BBd; ///< The bag parameters for \f$\Delta b=2\f$ processes for the \f$B_d\f$ meson system.
    BParameter BD; ///< The bag parameters for \f$\Delta c=2\f$ processes for the \f$D^0\f$ meson system.
    BParameter BK; ///< The bag parameters for \f$\Delta s=2\f$ processes for the \f$K^0\f$ meson system.
    BParameter BKd1;
    BParameter BKd3;
    double BBsB0, BBsB1, BBsB2, BBsB3, BBsB4, BBsscale;
    //double BBdB0, BBdB1, BBdB2, BBdB3, BBdB4, BBdscale;
    double BDB0, BDB1, BDB2, BDB3, BDB4, BDscale;
    double BKB0, BKB1, BKB2, BKB3, BKB4, BKscale;
    double BKd1B0, BKd1B1, BKd1B2, BKd1B3, BKd1B4, BKd1B5, BKd1B6, BKd1B7, BKd1B8, BKd1B9, BKd_scale;
    double BKd3B0, BKd3B1, BKd3B2, BKd3B3, BKd3B4, BKd3B5, BKd3B6, BKd3B7, BKd3B8, BKd3B9;

    double zeta2; ///< \f$\zeta(2)\f$ computed with the <a href="http://www.gnu.org/software/gsl/" target=blank>GSL</a>.
    double zeta3; ///< \f$\zeta(3)\f$ computed with the <a href="http://www.gnu.org/software/gsl/" target=blank>GSL</a>.
    bool computeFBd; ///< Swith for computing \f$F_{B_d}\f$ from \f$F_{B_s}\f$.
    bool computeBd; ///< Swith for computing \f$B_{B_d}\f$ from \f$B_{B_s}\f$.
    static const int CacheSize = 5; ///< Defines the depth of the cache.
    mutable double als_cache[8][CacheSize]; ///< Cache for \f$\alpha_s\f$.
    mutable double logLambda5_cache[4][CacheSize];
    mutable double logLambdaNLO_cache[9][CacheSize];
    mutable double mrun_cache[10][CacheSize]; ///< Cache for running quark mass.
    mutable double mp2mbar_cache[5][CacheSize]; ///< Cache for pole mass to msbar mass conversion.
    bool unknownParameterWarning; ///< A flag to stop the unknown parameter warning after the first time.

    /**
     * @brief The strong coupling constant computed with using \f$\Lambda_{\rm QCD}\f$.
     * @param[in] mu the scale of the strong coupling constant
     * @param[in] logLambda \f$\log(\Lambda_{\rm QCD})\f$
     * @param[in] order the %QCD order at which \f$\alpha_s\f$ is required
     * @return \f$\alpha_s(\mu)\f$ for the specified order
     */
    double AlsWithLambda(const double mu, const double logLambda, const orders order) const;

    /**
     * @brief A member for calculating the difference in \f$\alpha_s^{\mathrm{FULLNLO}}\f$ across the six-five
     * flavour threshold using AlsWithLambda().
     * @param[in] logLambda6 \f$\log(\Lambda_{\rm QCD}^{(6)})\f$
     * @param[in] logLambda5_in \f$\log(\Lambda_{\rm QCD}^{(5)})\f$
     * @return the difference \f$\alpha_s^{\mathrm{FULLNLO}}(\mu_t+\epsilon)\f$ - \f$\alpha_s^{\mathrm{FULLNLO}}(\mu_t-\epsilon)\f$
     * with \f$\epsilon = 10^{-10}\f$
     */
    double ZeroNf6NLO(double *logLambda6, double *logLambda5_in) const;

    /**
     * @brief A member for calculating the difference in \f$\alpha_s\f$ using
     * AlsWithLambda() and the input value of \f$\alpha_s(M_{\alpha_s})\f$ given as a model parameter.
     * @param[in] logLambda5 \f$\log(\Lambda_{\rm QCD}^{(5)})\f$
     * @param[in] order the %QCD order of the calculation
     * @return AlsWithLambda(\f$M_{\alpha_s}\f$, *logLambda5, *order) - \f$\alpha_s(M_{\alpha_s})\f$
     */
    double ZeroNf5(double *logLambda5, double *order) const;

    /**
     * @brief A member for calculating the difference in \f$\alpha_s^{\mathrm{FULLNLO}}\f$ across the four-five
     * flavour threshold using AlsWithLambda().
     * @param[in] logLambda4 \f$\log(\Lambda_{\rm QCD}^{(4)})\f$
     * @param[in] logLambda5_in \f$\log(\Lambda_{\rm QCD}^{(5)})\f$
     * @return the difference \f$\alpha_s^{\mathrm{FULLNLO}}(\mu_b-\epsilon)\f$ - \f$\alpha_s^{\mathrm{FULLNLO}}(\mu_b+\epsilon)\f$
     * with \f$\epsilon = 10^{-10}\f$
     */
    double ZeroNf4NLO(double *logLambda4, double *logLambda5_in) const;

    /**
     * @brief A member for calculating the difference in \f$\alpha_s^{\mathrm{FULLNLO}}\f$ across the three-four
     * flavour threshold using AlsWithLambda().
     * @param[in] logLambda3 \f$\log(\Lambda_{\rm QCD}^{(3)})\f$
     * @param[in] logLambda4_in \f$\log(\Lambda_{\rm QCD}^{(4)})\f$
     * @return the difference \f$\alpha_s^{\mathrm{FULLNLO}}(\mu_c-\epsilon)\f$ - \f$\alpha_s^{\mathrm{FULLNLO}}(\mu_c+\epsilon)\f$
     * with \f$\epsilon = 10^{-10}\f$
     */
    double ZeroNf3NLO(double *logLambda3, double *logLambda4_in) const;

    /**
     * @brief \f$\log(\Lambda_{\rm QCD})\f$ for \f$n_f = 5\f$.
     * @param[in] order the %QCD order of the computation
     * @return \f$\log(\Lambda_{\rm QCD}^{(5)})\f$
     */
    double logLambda5(orders order) const;

    /**
     * @brief \f$\log(\Lambda_{\rm QCD})\f$ used for computation of \f$\alpha_s\f$ at FULLNLO.
     * @param[in] nfNEW the number of flavours after crossing the flavour threshold
     * @param[in] nfORG the number of flavours before crossing the flavour threshold
     * @param[in] logLambdaORG the value of \f$\log(\Lambda_{\rm QCD})\f$ with \f$n_f = nfORG\f$
     * @return \f$\log(\Lambda_{\rm QCD})\f$ for \f$n_f = nfNEW\f$
     */
    double logLambdaNLO(const double nfNEW, const double nfORG, const double logLambdaORG) const;

    /**
     * @brief \f$\log(\Lambda_{\rm QCD})\f$ used for computation of \f$\alpha_s\f$ at FULLNNLO.
     * @param[in] muMatching the scale at which the matching is done during crossing a flavour threshold
     * @param[in] mf the mass of the quark sitting at the flavour threshold being crossed
     * @param[in] nfNEW the number of flavours after crossing the flavour threshold
     * @param[in] nfORG the number of flavours before crossing the flavour threshold
     * @param[in] logLambdaORG the value of \f$\log(\Lambda_{\rm QCD})\f$ with \f$n_f = nfORG\f$
     * @param[in] order the %QCD order of the calculation
     * @return \f$\log(\Lambda_{\rm QCD})\f$ for \f$n_f = nfNEW\f$
     */
    double logLambda(const double muMatching, const double mf,
            const double nfNEW, const double nfORG,
            const double logLambdaORG, orders order) const;

    /**
     * @brief The threshold correction for running of a mass when crossing a flavour threshold.
     * @param[in] nf_f the number of flavours \f$n_f\f$ after crossing the threshold
     * @param[in] nf_i the number of flavours \f$n_i\f$ before crossing the threshold
     * @return the threshold correction factoe
     */
    double threCorrForMass(const double nf_f, const double nf_i) const;

    /**
     * @brief A function to calculate the running of the mass between flavour thresholds.
     * @param[in] mu_f the final scale \f$\mu_f\f$ to which the mass if run
     * @param[in] mu_i the initial scale \f$\mu_i\f$ from which the mass if run
     * @param[in] m the mass at the scale \f$\mu_i\f$
     * @param[in] order the %QCD order at which the running is being calculated
     * @return the mass run from \f$\mu_i\f$ to \f$\mu_f\f$
     */
    double MrunTMP(const double mu_f, const double mu_i, const double m, const orders order) const;

    /**
     * @brief The member used for finding the numerical solution to the pole mass from the \f$\overline{\rm MS}\f$
     * mass.
     * @param[in] mu a pointer to the \f$\overline{\rm MS}\f$ mass
     * @param[in] params a pointer to a vector containing the pole mass and the %QCD order of the computation
     * @return the difference in the pole mass and the pole mass as computed from the \f$\overline{\rm MS}\f$ mass
     */
    double Mp2MbarTMP(double *mu, double *params) const;

    /**
     * @brief A member used to manage the caching for this class.
     * @param[in] cache the cache to be moved
     * @param[in] n the dimension of the cache to be shifted
     */
    void CacheShift(double cache[][5], int n) const;

    gslpp::complex myh_0,myh_p,myh_m,myh_0_1,myh_p_1,myh_m_1,myh_0_2,myh_p_2,myh_m_2;
};

#endif	/* QCD_H */
