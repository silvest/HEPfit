/* 
 * Copyright (C) 2012-2014 SusyFit Collaboration
 * All rights reserved.
 *
 * For the licensing terms see doc/COPYING.
 */

#ifndef THFACTORY_H
#define	THFACTORY_H

#include <ThObservable.h>
#include <ModelObservable.h>
#include <StandardModel.h>
#include <EW.h>
#include <Flavour.h>

/**
 * @class ThFactory
 * @ingroup InputParser 
 * @brief A class for an interface to theoretical predictions for observables. 
 * @author SusyFit Collaboration
 * @copyright GNU General Public License
 * @details The list of available observables are as follows:
 *  
 * Electroweak Precision Observables
 * <table class="observables">
 * <tr>
 *   <th>Label</th>
 *   <th>LaTeX symbol</th>
 *   <th>Description</th>
 * </tr>
 * <tr>
 *   <td class="obs_name">%Mw</td>
 *   <td class="obs_symb">\f$ M_W \f$</td>
 *   <td class="obs_desc">The mass of the \f$ W \f$ boson.</td>
 * </tr>
 * <tr>
 *   <td class="obs_name">%GammaW</td>
 *   <td class="obs_symb">\f$ \Gamma_W \f$</td>
 *   <td class="obs_desc">The total decay width of the \f$W\f$ boson.</td>
 * </tr>
 * <tr>
 *   <td class="obs_name">%GammaZ</td>
 *   <td class="obs_symb">\f$ \Gamma_Z \f$</td>
 *   <td class="obs_desc">The total decay width of the \f$Z\f$ boson.</td>
 * </tr>
 * <tr>
 *   <td class="obs_name">%sigmaHadron</td>
 *   <td class="obs_symb"> @f$\sigma_h^0@f$</td>
 *   <td class="obs_desc">The \f$e^+ e^- \rightarrow Z \rightarrow \mathrm{hadrons}\f$ cross section at the \f$Z\f$ pole.</td>
 * </tr>
 * <tr>
 *   <td class="obs_name">%sin2thetaEff</td>
 *   <td class="obs_symb">\f$\sin^2{\theta_\mathrm{Eff}^\ell}\f$</td>
 *   <td class="obs_desc">The leptonic effective weak mixing angle.</td>
 * </tr>
 * <tr>
 *   <td class="obs_name">%PtauPol</td>
 *   <td class="obs_symb">\f$P_\tau^\mathrm{pol}\f$</td>
 *   <td class="obs_desc">The longitudinal polarization of the final state in \f$e^+ e^- \rightarrow Z \rightarrow \tau^+ \tau^-\f$</td>
 * </tr>
 * <tr>
 *   <td class="obs_name">%Alepton</td>
 *   <td class="obs_symb">\f$A_\ell\f$</td>
 *   <td class="obs_desc">The leptonic left-right asymmetry parameter at the \f$Z\f$ pole (asummes lepton universality).</td>
 * </tr>
 * <tr>
 *   <td class="obs_name">%Acharm</td>
 *   <td class="obs_symb">\f$A_c\f$</td>
 *   <td class="obs_desc">The charm quark left-right asymmetry parameter.</td>
 * </tr>
 * <tr>
 *   <td class="obs_name">%Abottom</td>
 *   <td class="obs_symb">\f$A_b\f$</td>
 *   <td class="obs_desc">The bottom quark left-right asymmetry parameter.</td>
 * </tr>
 * <tr>
 *   <td class="obs_name">%AFBlepton</td>
 *   <td class="obs_symb">\f$A_{FB}^{0,\ell}\f$</td>
 *   <td class="obs_desc">The leptonic forward-backward asymmetry at the \f$Z\f$ pole (asummes lepton universality).</td>
 * </tr>
 * <tr>
 *   <td class="obs_name">%AFBcharm</td>
 *   <td class="obs_symb">\f$A_{FB}^{0,c}\f$</td>
 *   <td class="obs_desc">The charm quark forward-backward asymmetry at the \f$Z\f$ pole.</td>
 * </tr>
 * <tr>
 *   <td class="obs_name">%AFBbottom</td>
 *   <td class="obs_symb">\f$A_{FB}^{0,b}\f$</td>
 *   <td class="obs_desc">The bottom quark forward-backward asymmetry at the \f$Z\f$ pole.</td>
 * </tr>
 * <tr>
 *   <td class="obs_name">%Rlepton</td>
 *   <td class="obs_symb">\f$R_\ell^0\f$</td>
 *   <td class="obs_desc">The ratio between the \f$Z\rightarrow \mathrm{hadrons}\f$ and \f$Z\rightarrow \ell^+ \ell^-\f$ decay widths (one lepton flavour, asummes lepton universality).</td>
 * </tr>
 * <tr>
 *   <td class="obs_name">%Rcharm</td>
 *   <td class="obs_symb">\f$R_c^0\f$</td>
 *   <td class="obs_desc">The ratio between the \f$Z\rightarrow c\bar{c}\f$ and \f$Z\rightarrow \mathrm{hadrons}\f$ decay widths.</td>
 * </tr>
 * <tr>
 *   <td class="obs_name">%Rbottom</td>
 *   <td class="obs_symb">\f$R_b^0\f$</td>
 *   <td class="obs_desc">he ratio between the \f$Z\rightarrow b\bar{b}\f$ and \f$Z\rightarrow \mathrm{hadrons}\f$ decay widths</td>
 * </tr>
 * </table>
 *
 * Epsilon Parameters
 * <table class="observables">
 * <tr>
 *   <th>Label</th>
 *   <th>LaTeX symbol</th>
 *   <th>Description</th>
 * </tr>
 * <tr>
 *   <td class="obs_name">%epsilon1</td>
 *   <td class="obs_symb"></td>
 *   <td class="obs_desc"></td>
 * </tr>
 * <tr>
 *   <td class="obs_name">%epsilon2</td>
 *   <td class="obs_symb"></td>
 *   <td class="obs_desc"></td>
 * </tr>
 * <tr>
 *   <td class="obs_name">%epsilon3</td>
 *   <td class="obs_symb"></td>
 *   <td class="obs_desc"></td>
 * </tr>
 * <tr>
 *   <td class="obs_name">%epsilonb</td>
 *   <td class="obs_symb"></td>
 *   <td class="obs_desc"></td>
 * </tr>
 * </table>
 *
 * LEP-II Two-Fermion Processes
 * <table class="observables">
 * <tr>
 *   <th>Label</th>
 *   <th>LaTeX symbol</th>
 *   <th>Description</th>
 * </tr>
 * <tr>
 *   <td class="obs_name">%sigmaqLEP2_</td>
 *   <td class="obs_symb"></td>
 *   <td class="obs_desc">The \f$e^+ e^-\rightarrow\mathrm{hadrons}\f$ cross section measured at LEP2.</td>
 * </tr>
 * <tr>
 *   <td class="obs_name">%sigmamuLEP2_</td>
 *   <td class="obs_symb"></td>
 *   <td class="obs_desc">The \f$e^+ e^-\rightarrow \mu^+ \mu^-\f$ cross section measured at LEP2.</td>
 * </tr>
 * <tr>
 *   <td class="obs_name">%sigmatauLEP2_</td>
 *   <td class="obs_symb"></td>
 *   <td class="obs_desc">The \f$e^+ e^-\rightarrow \tau^+ \tau^-\f$ cross section measured at LEP2.</td>
 * </tr>
 * <tr>
 *   <td class="obs_name">%AFBmuLEP2_</td>
 *   <td class="obs_symb"></td>
 *   <td class="obs_desc">The forward-backward asymmetry in \f$e^+ e^-\rightarrow \mu^+ \mu^-\f$ measured at LEP2.</td>
 * </tr>
 * <tr>
 *   <td class="obs_name">%AFBtauLEP2_</td>
 *   <td class="obs_symb"></td>
 *   <td class="obs_desc">The forward-backward asymmetry in \f$e^+ e^-\rightarrow \tau^+ \tau^-\f$ measured at LEP2.</td>
 * </tr>
 * <tr>
 *   <td class="obs_name">%AFBbottomLEP2_</td>
 *   <td class="obs_symb"></td>
 *   <td class="obs_desc">The forward-backward asymmetry in \f$e^+ e^-\rightarrow b \bar{b}\f$ measured at LEP2.</td>
 * </tr>
 * <tr>
 *   <td class="obs_name">%AFBcharmLEP2_</td>
 *   <td class="obs_symb"></td>
 *   <td class="obs_desc">The forward-backward asymmetry in \f$e^+ e^-\rightarrow c \bar{c}\f$ measured at LEP2.</td>
 * </tr>
 * <tr>
 *   <td class="obs_name">%RbottomLEP2_</td>
 *   <td class="obs_symb"></td>
 *   <td class="obs_desc">The ratio of the \f$e^+ e^-\rightarrow b \bar{b}\f$ and \f$e^+ e^-\rightarrow\mathrm{hadrons}\f$ cross sections measured at LEP2.</td>
 * </tr>
 * <tr>
 *   <td class="obs_name">%RcharmLEP2_</td>
 *   <td class="obs_symb"></td>
 *   <td class="obs_desc">The ratio of the \f$e^+ e^-\rightarrow c \bar{c}\f$ and \f$e^+ e^-\rightarrow\mathrm{hadrons}\f$ cross sections measured at LEP2.</td>
 * </tr>
 * </table>
 *
 * Flavour Observables
 * <table class="observables">
 * <tr>
 *   <th>Label</th>
 *   <th>LaTeX symbol</th>
 *   <th>Description</th>
 * </tr>
 * <tr>
 *   <td class="obs_name">%Dmd1</td>
 *   <td class="obs_symb"></td>
 *   <td class="obs_desc"></td>
 * </tr>
 * <tr>
 *   <td class="obs_name">%Dms1</td>
 *   <td class="obs_symb"></td>
 *   <td class="obs_desc"></td>
 * </tr>
 * <tr>
 *   <td class="obs_name">%M12D</td>
 *   <td class="obs_symb"></td>
 *   <td class="obs_desc"></td>
 * </tr>
 * <tr>
 *   <td class="obs_name">%ArgD</td>
 *   <td class="obs_symb"></td>
 *   <td class="obs_desc"></td>
 * </tr>
 * <tr>
 *   <td class="obs_name">%EpsilonK</td>
 *   <td class="obs_symb"></td>
 *   <td class="obs_desc"></td>
 * </tr>
 * <tr>
 *   <td class="obs_name">%EpsiloP_o_Epsilon</td>
 *   <td class="obs_symb"></td>
 *   <td class="obs_desc"></td>
 * </tr>
 * <tr>
 *   <td class="obs_name">%DmK</td>
 *   <td class="obs_symb"></td>
 *   <td class="obs_desc"></td>
 * </tr>
 * <tr>
 *   <td class="obs_name">%Vud</td>
 *   <td class="obs_symb"></td>
 *   <td class="obs_desc"></td>
 * </tr>
 * <tr>
 *   <td class="obs_name">%Vus</td>
 *   <td class="obs_symb"></td>
 *   <td class="obs_desc"></td>
 * </tr>
 * <tr>
 *   <td class="obs_name">%Vub</td>
 *   <td class="obs_symb"></td>
 *   <td class="obs_desc"></td>
 * </tr>
 * <tr>
 *   <td class="obs_name">%Vcb</td>
 *   <td class="obs_symb"></td>
 *   <td class="obs_desc"></td>
 * </tr>
 * <tr>
 *   <td class="obs_name">%alpha</td>
 *   <td class="obs_symb"></td>
 *   <td class="obs_desc"></td>
 * </tr>
 * <tr>
 *   <td class="obs_name">%alpha_2a</td>
 *   <td class="obs_symb"></td>
 *   <td class="obs_desc"></td>
 * </tr>
 * <tr>
 *   <td class="obs_name">%gamma</td>
 *   <td class="obs_symb"></td>
 *   <td class="obs_desc"></td>
 * </tr>
 * <tr>
 *   <td class="obs_name">%SJPsiK</td>
 *   <td class="obs_symb"></td>
 *   <td class="obs_desc"></td>
 * </tr>
 * <tr>
 *   <td class="obs_name">%SJPsiPhi</td>
 *   <td class="obs_symb"></td>
 *   <td class="obs_desc"></td>
 * </tr>
 * <tr>
 *   <td class="obs_name">%BR_Bdmumu</td>
 *   <td class="obs_symb"></td>
 *   <td class="obs_desc"></td>
 * </tr>
 * <tr>
 *   <td class="obs_name">%BRbar_Bdmumu</td>
 *   <td class="obs_symb"></td>
 *   <td class="obs_desc"></td>
 * </tr>
 * <tr>
 *   <td class="obs_name">%Amumu_Bd</td>
 *   <td class="obs_symb"></td>
 *   <td class="obs_desc"></td>
 * </tr>
 * <tr>
 *   <td class="obs_name">%Smumu_Bd</td>
 *   <td class="obs_symb"></td>
 *   <td class="obs_desc"></td>
 * </tr>
 * <tr>
 *   <td class="obs_name">%BR_Bsmumu</td>
 *   <td class="obs_symb"></td>
 *   <td class="obs_desc"></td>
 * </tr>
 * <tr>
 *   <td class="obs_name">%BRbar_Bsmumu</td>
 *   <td class="obs_symb"></td>
 *   <td class="obs_desc"></td>
 * </tr>
 * <tr>
 *   <td class="obs_name">%Amumu_Bs</td>
 *   <td class="obs_symb"></td>
 *   <td class="obs_desc"></td>
 * </tr>
 * <tr>
 *   <td class="obs_name">%Smumu_Bs</td>
 *   <td class="obs_symb"></td>
 *   <td class="obs_desc"></td>
 * </tr>
 * </table>
 *
 * Standard %Model Input Parameters
 * <table class="observables">
 * <tr>
 *   <th>Label</th>
 *   <th>LaTeX symbol</th>
 *   <th>Description</th>
 * </tr>
 * <tr>
 *   <td class="obs_name">%AlsMz</td>
 *   <td class="obs_symb">\f$\alpha_s(M_Z) \f$</td>
 *   <td class="obs_desc">The strong coupling constat at the \f$ Z \f$ pole.</td>
 * </tr>
 * <tr>
 *   <td class="obs_name">%dAle5Mz</td>
 *   <td class="obs_symb"></td>
 *   <td class="obs_desc"></td>
 * </tr>
 * <tr>
 *   <td class="obs_name">%Mz</td>
 *   <td class="obs_symb">\f$ M_Z \f$</td>
 *   <td class="obs_desc">The mass of the \f$ Z \f$ boson.</td>
 * </tr>
 * <tr>
 *   <td class="obs_name">%mtop</td>
 *   <td class="obs_symb">The mass of the top quark</td>
 *   <td class="obs_desc">\f$ M_{top} \f$</td>
 * </tr>
 * <tr>
 *   <td class="obs_name">%mHl</td>
 *   <td class="obs_symb"></td>
 *   <td class="obs_desc"></td>
 * </tr>
 * <tr>
 *   <td class="obs_name">%delMw</td>
 *   <td class="obs_symb"></td>
 *   <td class="obs_desc"></td>
 * </tr>
 * <tr>
 *   <td class="obs_name">%delSin2th_l</td>
 *   <td class="obs_symb"></td>
 *   <td class="obs_desc"></td>
 * </tr>
 * <tr>
 *   <td class="obs_name">%delGammaZ</td>
 *   <td class="obs_symb"></td>
 *   <td class="obs_desc"></td>
 * </tr>
 * </table>
 *
 * NP Input Parameters for NPEffective1 and NPeffective2
 * <table class="observables">
 * <tr>
 *   <th>Label</th>
 *   <th>LaTeX symbol</th>
 *   <th>Description</th>
 * </tr>
 * <tr>
 *   <td class="obs_name">%cHLp_NP</td>
 *   <td class="obs_symb"></td>
 *   <td class="obs_desc"></td>
 * </tr>
 * <tr>
 *   <td class="obs_name">%cHQp_NP</td>
 *   <td class="obs_symb"></td>
 *   <td class="obs_desc"></td>
 * </tr>
 * <tr>
 *   <td class="obs_name">%cHQ_NP</td>
 *   <td class="obs_symb"></td>
 *   <td class="obs_desc"></td>
 * </tr>
 * <tr>
 *   <td class="obs_name">%cHL_NP</td>
 *   <td class="obs_symb"></td>
 *   <td class="obs_desc"></td>
 * </tr>
 * <tr>
 *   <td class="obs_name">%cHE_NP</td>
 *   <td class="obs_symb"></td>
 *   <td class="obs_desc"></td>
 * </tr>
 * <tr>
 *   <td class="obs_name">%c_Ae_NP</td>
 *   <td class="obs_symb"></td>
 *   <td class="obs_desc"></td>
 * </tr>
 * <tr>
 *   <td class="obs_name">%c_GammaZ_uds_NP</td>
 *   <td class="obs_symb"></td>
 *   <td class="obs_desc"></td>
 * </tr>
 * <tr>
 *   <td class="obs_name">%cHU2_NP</td>
 *   <td class="obs_symb"></td>
 *   <td class="obs_desc"></td>
 * </tr>
 * <tr>
 *   <td class="obs_name">%cHD3_NP</td>
 *   <td class="obs_symb"></td>
 *   <td class="obs_desc"></td>
 * </tr>
 * <tr>
 *   <td class="obs_name">%cHQ1pPLUScHQ2p_NP</td>
 *   <td class="obs_symb"></td>
 *   <td class="obs_desc"></td>
 * </tr>
 * <tr>
 *   <td class="obs_name">%cHQ2pMINUScHQ2_NP</td>
 *   <td class="obs_symb"></td>
 *   <td class="obs_desc"></td>
 * </tr>
 * <tr>
 *   <td class="obs_name">%cHQ3pPLUScHQ3_NP</td>
 *   <td class="obs_symb"></td>
 *   <td class="obs_desc"></td>
 * </tr>
 * </table>
 *
 * NP Input Parameters for NPZbbbar
 * <table class="observables">
 * <tr>
 *   <th>Label</th>
 *   <th>LaTeX symbol</th>
 *   <th>Description</th>
 * </tr>
 * <tr>
 *   <td class="obs_name">%deltaGVb</td>
 *   <td class="obs_symb"></td>
 *   <td class="obs_desc"></td>
 * </tr>
 * <tr>
 *   <td class="obs_name">%deltaGAb</td>
 *   <td class="obs_symb"></td>
 *   <td class="obs_desc"></td>
 * </tr>
 * <tr>
 *   <td class="obs_name">%deltaGLb</td>
 *   <td class="obs_symb"></td>
 *   <td class="obs_desc"></td>
 * </tr>
 * <tr>
 *   <td class="obs_name">%deltaGRb</td>
 *   <td class="obs_symb"></td>
 *   <td class="obs_desc"></td>
 * </tr>
 * <tr>
 *   <td class="obs_name">%deltaRhoZb</td>
 *   <td class="obs_symb"></td>
 *   <td class="obs_desc"></td>
 * </tr>
 * <tr>
 *   <td class="obs_name">%deltaKappaZb</td>
 *   <td class="obs_symb"></td>
 *   <td class="obs_desc"></td>
 * </tr>
 * </table>
 *
 * SUSY Spectra and Observables
 * <table class="observables">
 * <tr>
 *   <th>Label</th>
 *   <th>LaTeX symbol</th>
 *   <th>Description</th>
 * </tr>
 * <tr>
 *   <td class="obs_name">%OutputSLHAfromFH</td>
 *   <td class="obs_symb"></td>
 *   <td class="obs_desc"></td>
 * </tr>
 * <tr>
 *   <td class="obs_name">%MHl</td>
 *   <td class="obs_symb"></td>
 *   <td class="obs_desc"></td>
 * </tr>
 * <tr>
 *   <td class="obs_name">%MHh</td>
 *   <td class="obs_symb"></td>
 *   <td class="obs_desc"></td>
 * </tr>
 * <tr>
 *   <td class="obs_name">%MHa</td>
 *   <td class="obs_symb"></td>
 *   <td class="obs_desc"></td>
 * </tr>
 * <tr>
 *   <td class="obs_name">%MHp</td>
 *   <td class="obs_symb"></td>
 *   <td class="obs_desc"></td>
 * </tr>
 * <tr>
 *   <td class="obs_name">%Msu1</td>
 *   <td class="obs_symb"></td>
 *   <td class="obs_desc"></td>
 * </tr>
 * <tr>
 *   <td class="obs_name">%Msu2</td>
 *   <td class="obs_symb"></td>
 *   <td class="obs_desc"></td>
 * </tr>
 * <tr>
 *   <td class="obs_name">%Msu3</td>
 *   <td class="obs_symb"></td>
 *   <td class="obs_desc"></td>
 * </tr>
 * <tr>
 *   <td class="obs_name">%Msu4</td>
 *   <td class="obs_symb"></td>
 *   <td class="obs_desc"></td>
 * </tr>
 * <tr>
 *   <td class="obs_name">%Msu5</td>
 *   <td class="obs_symb"></td>
 *   <td class="obs_desc"></td>
 * </tr>
 * <tr>
 *   <td class="obs_name">%Msu6</td>
 *   <td class="obs_symb"></td>
 *   <td class="obs_desc"></td>
 * </tr>
 * <tr>
 *   <td class="obs_name">%Msd1</td>
 *   <td class="obs_symb"></td>
 *   <td class="obs_desc"></td>
 * </tr>
 * <tr>
 *   <td class="obs_name">%Msd2</td>
 *   <td class="obs_symb"></td>
 *   <td class="obs_desc"></td>
 * </tr>
 * <tr>
 *   <td class="obs_name">%Msd3</td>
 *   <td class="obs_symb"></td>
 *   <td class="obs_desc"></td>
 * </tr>
 * <tr>
 *   <td class="obs_name">%Msd4</td>
 *   <td class="obs_symb"></td>
 *   <td class="obs_desc"></td>
 * </tr>
 * <tr>
 *   <td class="obs_name">%Msd5</td>
 *   <td class="obs_symb"></td>
 *   <td class="obs_desc"></td>
 * </tr>
 * <tr>
 *   <td class="obs_name">%Msd6</td>
 *   <td class="obs_symb"></td>
 *   <td class="obs_desc"></td>
 * </tr>
 * <tr>
 *   <td class="obs_name">%Mch1</td>
 *   <td class="obs_symb"></td>
 *   <td class="obs_desc"></td>
 * </tr>
 * <tr>
 *   <td class="obs_name">%Mch2</td>
 *   <td class="obs_symb"></td>
 *   <td class="obs_desc"></td>
 * </tr>
 * <tr>
 *   <td class="obs_name">%Mneu1</td>
 *   <td class="obs_symb"></td>
 *   <td class="obs_desc"></td>
 * </tr>
 * <tr>
 *   <td class="obs_name">%Mneu2</td>
 *   <td class="obs_symb"></td>
 *   <td class="obs_desc"></td>
 * </tr>
 * <tr>
 *   <td class="obs_name">%Mneu3</td>
 *   <td class="obs_symb"></td>
 *   <td class="obs_desc"></td>
 * </tr>
 * <tr>
 *   <td class="obs_name">%Mneu4</td>
 *   <td class="obs_symb"></td>
 *   <td class="obs_desc"></td>
 * </tr>
 * <tr>
 *   <td class="obs_name">%Mw_dRho</td>
 *   <td class="obs_symb"></td>
 *   <td class="obs_desc"></td>
 * </tr>
 * </table>
 * Masses and widhts are computed in GeV. Cross sections are given in nb.
 * 
 */
class ThFactory {
public:
    
    /**
     * @brief Constructor.
     * @details All the theory observables are built by the constructor. These
     * are implemented as a map between the names of the observables and a
     * pointer to an object of type ThObservable.
     * @param[in] myModel a reference to an object of the type StandardModel
     */
    ThFactory(const StandardModel& myModel);
    
    /**
     * @brief The default destructor. It destroys all the non-NULL pointers
     * in the thobs map.
     */
    virtual ~ThFactory();
    
    /**
     * @brief This method checks for the existence of an observable of a specific name in the map
     * thobs and returns a pointer to the particular ThObservable.
     * @param[in] name the name of the observable
     * @return a pointer to an object of the type ThObservable
     */
    ThObservable* getThMethod(const std::string& name);

private:
    std::map<std::string, ThObservable *> thobs; ///< A map between a pointer of type ThObservable and the name of the observable.
    EW myEW; ///< An object of type EW().
    Flavour myFlavour; ///< An object of type Flavour.
    ModelObservable myMO; ///< An object of type ModelObservable
};

#endif	/* THFACTORY_H */