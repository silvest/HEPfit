/*
 * Copyright (C) 2012 HEPfit Collaboration
 * All rights reserved.
 *
 * For the licensing terms see doc/COPYING.
 */

#ifndef STANDARDMODEL_H
#define	STANDARDMODEL_H

#include <gslpp.h>
#include "QCD.h"
#include "CKM.h"
#include "WilsonCoefficient.h"
#include "StandardModelMatching.h"

class EWSMcache;
class EWSMOneLoopEW;
class EWSMTwoLoopQCD;
class EWSMTwoLoopEW;
class EWSMThreeLoopQCD;
class EWSMThreeLoopEW2QCD;
class EWSMThreeLoopEW;
class EWSMApproximateFormulae;
class Flavour;
/** BEGIN: REMOVE FROM THE PACKAGE **/
class EWSMTwoFermionsLEP2;
class LeptonFlavour;
/** END: REMOVE FROM THE PACKAGE **/


/**
 * @class StandardModel
 * @ingroup StandardModel
 * @brief A model class for the Standard %Model.
 * @author HEPfit Collaboration
 * @copyright GNU General Public License
 * @details This is a Model class containing parameters and functions associated
 * with the Standard %Model. This class is inherited from the QCD class, which
 * defines parameters related to %QCD.
 *
 *
 * @anchor StandardModelInitialization
 * <h3>Initialization</h3>
 *
 * The constructor StandardModel() initializes some of the model flags to their
 * default values. After creating an instance of the current class,
 * it is required to call the initialization method InitializeModel(), which
 * allocates memory to the pointers defined in the current class.
 * These pointers are then used in computing %EW precision and flavour observables,
 * respectively. In the Monte Carlo run, the constructor as well as the initialization
 * method are called in InputParser::ReadParameters().
 *
 * The initializations and updates of the model parameters and flags are explained
 * below. 
 *
 * 
 * @anchor StandardModelParameters
 * <h3>%Model parameters</h3>
 *
 * The model parameters of %StandardModel are summarized below:
 * <table class="model">
 * <tr>
 *   <th>Label</th>
 *   <th>LaTeX symbol</th>
 *   <th>Description</th>
 * </tr>
 * <tr>
 *   <td class="mod_name">%Mz</td>
 *   <td class="mod_symb">@f$M_Z@f$</td>
 *   <td class="mod_desc">The mass of the \f$Z\f$ boson in GeV.</td>
 * </tr>
 * <tr>
 *   <td class="mod_name">%AlsMz</td>
 *   <td class="mod_symb">@f$\alpha_s(M_Z)@f$</td>
 *   <td class="mod_desc">The strong coupling constant at the Z-boson mass.</td>
 * </tr>
 * <tr>
 *   <td class="mod_name">%GF</td>
 *   <td class="mod_symb">@f$G_\mu@f$</td>
 *   <td class="mod_desc">The Fermi constant in @f${\rm GeV}^{-2}@f$, measured through muon decays.</td>
 * </tr>
 * <tr>
 *   <td class="mod_name">%ale</td>
 *   <td class="mod_symb">@f$\alpha@f$</td>
 *   <td class="mod_desc">The fine-structure constant.</td>
 * </tr>
 * <tr>
 *   <td class="mod_name">%dAle5Mz</td>
 *   <td class="mod_symb">@f$\Delta\alpha_{\mathrm{had}}^{(5)}(M_Z^2)@f$</td>
 *   <td class="mod_desc">The five-flavour hadronic contribution to the electromagnetic coupling.</td>
 * </tr>
 * <tr>
 *   <td class="mod_name">%mHl</td>
 *   <td class="mod_symb">@f$m_h@f$</td>
 *   <td class="mod_desc">The Higgs mass in GeV.</td>
 * </tr>
 * <tr>
 *   <td class="mod_name">%delMw</td>
 *   <td class="mod_symb">@f$\delta\,M_W@f$</td>
 *   <td class="mod_desc">The theoretical uncertainty in @f$M_W@f$ in GeV,
 *   which is applicable only when EWSMApproximateFormulae::Mw() is employed
 *   for @f$M_W@f$. See also the model flag @ref StandardModelFlags "Mw".</td>
 * </tr>
 * <tr>
 *   <td class="mod_name">%delSin2th_l</td>
 *   <td class="mod_symb">@f$\delta\sin^2\theta_{\rm eff}^{\rm lept}@f$</td>
 *   <td class="mod_desc">The theoretical uncertainty in @f$\sin^2\theta_{\rm eff}^{\rm lept}@f$,
 *   which is applicable only when EWSMApproximateFormulae::sin2thetaEff_l()
 *   is employed for @f$\sin^2\theta_{\rm eff}^{\rm lept}@f$.
 *   See also the model flag @ref StandardModelFlags "KappaZ".</td>
 * </tr>
 * <tr>
 *   <td class="mod_name">%delGammaZ</td>
 *   <td class="mod_symb">@f$\delta\,\Gamma_Z@f$</td>
 *   <td class="mod_desc">The theoretical uncertainty in @f$\Gamma_Z@f$ in GeV,
 *   which is applicable only when EWSMApproximateFormulae::X_extended() is 
 *   employed for @f$\Gamma_Z@f$. See also the model flag 
 *   @ref StandardModelFlags "NoApproximateGammaZ".</td>
 * </tr>
 * <tr>
 *   <td class="mod_name">%mneutrino_1</td>
 *   <td class="mod_symb">@f$m_{\nu_1}@f$</td>
 *   <td class="mod_desc">The mass of the first-generation neutrino in GeV.</td>
 * </tr>
 * <tr>
 *   <td class="mod_name">%mneutrino_2</td>
 *   <td class="mod_symb">@f$m_{\nu_2}@f$</td>
 *   <td class="mod_desc">The mass of the second-generation neutrino in GeV.</td>
 * </tr>
 * <tr>
 *   <td class="mod_name">%mneutrino_3</td>
 *   <td class="mod_symb">@f$m_{\nu_3}@f$</td>
 *   <td class="mod_desc">The mass of the third-generation neutrino in GeV.</td>
 * </tr>
 * <tr>
 *   <td class="mod_name">%melectron</td>
 *   <td class="mod_symb">@f$m_e@f$</td>
 *   <td class="mod_desc">The electron mass in GeV.</td>
 * </tr>
 * <tr>
 *   <td class="mod_name">%mmu</td>
 *   <td class="mod_symb">@f$m_\mu@f$</td>
 *   <td class="mod_desc">The muon mass in GeV.</td>
 * </tr>
 * <tr>
 *   <td class="mod_name">%mtau</td>
 *   <td class="mod_symb">@f$m_\tau@f$</td>
 *   <td class="mod_desc">The tau mass in GeV.</td>
 * </tr>
 * <tr>
 *   <td class="mod_name">%lambda</td>
 *   <td class="mod_symb">@f$\lambda@f$</td>
 *   <td class="mod_desc">The %CKM parameter @f$\lambda@f$ in the Wolfenstein parameterization.</td>
 * </tr>
 * <tr>
 *   <td class="mod_name">%A</td>
 *   <td class="mod_symb">@f$A@f$</td>
 *   <td class="mod_desc">The %CKM parameter @f$A@f$ in the Wolfenstein parameterization.</td>
 * </tr>
 * <tr>
 *   <td class="mod_name">%rhob</td>
 *   <td class="mod_symb">@f$\bar{\rho}@f$</td>
 *   <td class="mod_desc">The %CKM parameter @f$\bar{\rho}@f$ in the Wolfenstein parameterization.</td>
 * </tr>
 * <tr>
 *   <td class="mod_name">%etab</td>
 *   <td class="mod_symb">@f$\bar{\eta}@f$</td>
 *   <td class="mod_desc">The %CKM parameter @f$\bar{\eta}@f$ in the Wolfenstein parameterization.</td>
 * </tr>
 * <tr>
 *   <td class="mod_name">%muw</td>
 *   <td class="mod_symb">@f$\mu_W@f$</td>
 *   <td class="mod_desc">A matching scale around the weak scale in GeV.</td>
 * </tr>
 * </table>
 * 
 * The parameters below, associated with flavour observables, may be removed
 * from the current class later:
 * <table class="model">
 * <tr>
 *   <td class="mod_name">%EpsK</td>
 *   <td class="mod_symb">@f$\varepsilon_{K}@f$</td>
 *   <td class="mod_desc">The experimental value of @f$\varepsilon_{K}@f$.</td>
 * </tr>
 * <tr>
 *   <td class="mod_name">%phiEpsK</td>
 *   <td class="mod_symb"></td>
 *   <td class="mod_desc"></td>
 * </tr>
 * <tr>
 *   <td class="mod_name">%KbarEpsK</td>
 *   <td class="mod_symb"></td>
 *   <td class="mod_desc"></td>
 * </tr>
 * <tr>
 *   <td class="mod_name">%DeltaMK</td>
 *   <td class="mod_symb">@f$\Delta m_{K}@f$</td>
 *   <td class="mod_desc">The experimental value of @f$\Delta m_{K}@f$ in GeV.</td>
 * </tr>
 * <tr>
 *   <td class="mod_name">%Dmk</td>
 *   <td class="mod_symb">@f$\Delta m_{K}@f$</td>
 *   <td class="mod_desc">The SM contribution to @f$\Delta m_{K}@f$ in GeV.</td>
 * </tr>
 * <tr>
 *   <td class="mod_name">%SM_M12D</td>
 *   <td class="mod_symb">@f$@f$</td>
 *   <td class="mod_desc">The SM amplitude of the @f$D^{0}-\bar{D}^{0}@f$ mixing.</td>
 * </tr>
 * </table>
 *
 * The set of the model parameters are initialized and updated with the methods
 * Init() and Update(), respectively, where the former calls the latter actually.
 * In Update(), the methods PreUpdate() and PostUpdate() are called to run all
 * the procedures that are need to be executed before and after the model parameters
 * are updated. The %CKM and %PMNS matrices and the Yukawa matrices are recomputed
 * in PostUpdate() with the updated parameters. Inside the Update() method, the
 * individual model parameter is assigned with the protected member function
 * setParameter(). 
 *
 * The parameters delMw, delSin2th_l, delGammaZ represent theoretical uncertainties 
 * in the @f$W@f$-boson mass, the leptonic effective weak mixing angle at the 
 * @f$Z@f$-boson mass scale and the total decay width of the @f$Z@f$ boson, 
 * respectively, originating from missing higher-order corrections. The contributions
 * from these parameters are incorporated into their two-loop approximate formulae:
 * EWSMApproximateFormulae::Mw(), EWSMApproximateFormulae::sin2thetaEff_l() and
 * EWSMApproximateFormulae::X_extended("GammaZ"). Therefore, the parameters are
 * applicable only when the corresponding approximate formulae are employed.
 * See also the model flags below. 
 *
 *
 * @anchor StandardModelFlags
 * <h3>%Model flags</h3>
 *
 * The flags of %StandardModel are summarized below, where the values of the boolean 
 * flags (TRUE or FALSE) are case insensitive, while those of the other flags are
 * case sensitive. The default values of the flags are indicated in bold:
 * 
 * <table class="model">
 * <tr>
 *   <th>Label</th>
 *   <th>Value</th>
 *   <th>Description</th>
 * </tr>
 * <tr>
 *   <td class="mod_name">%CacheInStandardModel</td>
 *   <td class="mod_valu"><b>TRUE</b>&nbsp;/&nbsp;FALSE</td>
 *   <td class="mod_desc">This flag controls the use of the cashing method
 *   implemented in EWSM class. The default value is TRUE.</td>
 * </tr>
 * <tr>
 *   <td class="mod_name">%CacheInEWSMcache</td>
 *   <td class="mod_valu"><b>TRUE</b>&nbsp;/&nbsp;FALSE</td>
 *   <td class="mod_desc">This flag controls the use of the cashing method
 *   implemented in EWSMcache class. The default value is TRUE.</td>
 * </tr>
 * <tr>
 *   <td class="mod_name">%WithoutNonUniversalVC</td>
 *   <td class="mod_valu">TRUE&nbsp;/&nbsp;<b>FALSE</b></td>
 *   <td class="mod_desc">This flag controls if flavour non-universal vertex
 *    corrections are not added to the epsilon parameterization for the %EW
 *    precision observables. The default value is FALSE; the non-universal
 *    corrections are taken into account.</td>
 * </tr>
 * <tr>
 *   <td class="mod_name">%NoApproximateGammaZ</td>
 *   <td class="mod_valu">TRUE&nbsp;/&nbsp;<b>FALSE</b></td>
 *   <td class="mod_desc">This flag is set to true if the two-loop approximate
 *   formulae of the partial and total decay widths of the @f$Z@f$ boson 
 *   defined with the function EWSMApproximateFormulae::X_extended() are NOT
 *   employed. The default value is FALSE.</td>
 * </tr>
 * <tr>
 *   <td class="mod_name">%Mw</td>
 *   <td class="mod_valu">NORESUM / OMSI / INTERMEDIATE / OMSII / <b>APPROXIMATEFORMULA</b></td>
 *   <td class="mod_desc">This flag controls the formula used in computing the
 *   @f$W@f$-boson mass. The default flag is APPROXIMATEFORMULA.
 *   See EWSM::Mw_SM(), EWSM::resumMw() and EWSMApproximateFormulae::Mw() for detail.</td>
 * </tr>
 * <tr>
 *   <td class="mod_name">%RhoZ</td>
 *   <td class="mod_valu"><b>NORESUM</b> / OMSI / INTERMEDIATE / OMSII</td>
 *   <td class="mod_desc">This flag controls the formula used in computing the
 *   @f$Zf\bar{f}@f$ couplings @f$\rho_Z^f@f$. The default flag is NORESUM.
 *   See EWSM::rhoZ_l_SM(), EWSM::rhoZ_q_SM() and EWSM::resumRhoZ() for detail.</td>
 * </tr>
 * <tr>
 *   <td class="mod_name">%KappaZ</td>
 *   <td class="mod_valu">NORESUM / OMSI / INTERMEDIATE / OMSII / <b>APPROXIMATEFORMULA</b></td>
 *   <td class="mod_desc">This flag controls the formula used in computing the
 *   @f$Zf\bar{f}@f$ couplings @f$\kappa_Z^f@f$. The default flag is APPROXIMATEFORMULA.
 *   See EWSM::kappaZ_l_SM(), EWSM::kappaZ_q_SM() and EWSM::resumKappaZ() for detail.</td>
 * </tr>
 * </table>
 *
 * These flags can be set via the method setFlag() or setFlagStr(), where the
 * former is applicable for the boolean flags, while the latter is for the other
 * flags. The method CheckFlags() is responsible for checking whether the flags
 * are sane. The public member functions
 * IsFlagWithoutNonUniversalVC(), IsFlagNoApproximateGammaZ()
 * getFlagMw(), getFlagRhoZ() and getFlagKappaZ()
 * are used to retrieve the values of each flag.
 *
 * The first two flags CacheInStandardModel and CacheInEWSMcache for the cashing methods
 * in EWSM and EWSMcache classes are relevant to the computations of the electroweak 
 * precision observables. Those caches are effective when the @f$W@f$-boson mass,
 * the decay widths of the @f$Z@f$ boson and the @f$Zf\bar{f}@f$ effective couplings
 * @f$\kappa_Z^f@f$ are calculated without using their two-loop approximate formulae.
 *
 * <h3>Notation</h3>
 *
 * The on-mass-shell renormalization scheme @cite Sirlin:1980nh,
 * @cite Marciano:1980pb, @cite Bardin:1980fe, @cite Bardin:1981sv is adopted
 * for UV divergences, and the weak mixing angle is defined in terms of the
 * physical masses of the gauge bosons:
 * @f[
 * s_W^2 \equiv \sin^2\theta_W = 1 - \frac{M_W^2}{M_Z^2}\,,
 * @f]
 * and @f$c_W^2=1-s_W^2@f$.

 * The Fermi constant @f$G_\mu@f$ in @f$\mu@f$ decay is taken as an input
 * quantity instead of the @f$W@f$-boson mass, since the latter has not been
 * measured very precisely compared to the former. The relation between
 * @f$G_\mu@f$ and @f$M_W@f$ is written as
 * @f[
 * G_\mu = \frac{\pi\,\alpha}{\sqrt{2} s_W^2 M_W^2} (1+\Delta r)\,,
 * @f]
 * where @f$\Delta r@f$ represents radiative corrections. From this relation,
 * the @f$W@f$-boson mass is calculated as
 * @f[
 * M_W^2
 * = \frac{M_Z^2}{2}
 * \left( 1+\sqrt{1-\frac{4\pi\alpha}{\sqrt{2}G_\mu M_Z^2}\,(1+\Delta r)}\
 * \right).
 * @f]
 *
 * The interaction between the @f$Z@f$ boson and the neutral current can be
 * written in terms of the effective @f$Zf\bar{f}@f$ couplings @f$g_{V}^f@f$
 * and @f$g_{A}^f@f$, of @f$g_{R}^f@f$ and @f$g_{L}^f@f$, or of @f$\rho_Z^f@f$
 * and @f$\kappa_Z^f@f$:
 * @f{eqnarray}{
 * \mathcal{L}
 * &=&
 * \frac{e}{2 s_W c_W}\,
 * Z_\mu \sum_f \bar{f}
 * \left( g_{V}^f\gamma_\mu - g_{A}^f \gamma_\mu\gamma_5 \right)\, f\,,
 * \\
 * &=&
 * \frac{e}{2s_W c_W}\,
 * Z_\mu \sum_f \bar{f}
 * \left[ g_{R}^f \gamma_\mu (1 + \gamma_5)
 * + g_{L}^f \gamma_\mu (1 - \gamma_5) \right]\, f\,,
 * \\
 * &=&
 * \frac{e}{2 s_W c_W}\sqrt{\rho_Z^f}\,
 * Z_\mu \sum_f \bar{f}
 * \left[( I_3^f - 2Q_f\kappa_Z^f s_W^2)\gamma^\mu
 *   - I_3^f\gamma^\mu\gamma_5\right]\,f\,,
 * @f}
 * where @f$\rho_Z^f@f$ and @f$\kappa_Z^f@f$ are related to
 * @f$g_{V}^f@f$ and @f$g_{A}^f@f$ as the relations:
 * @f{eqnarray}{
 * g_V^f
 * &=&
 * \sqrt{\rho_Z^f} I_3^f (1 - 4|Q_f|\kappa_Z^fs_W^2)
 * = \sqrt{\rho_Z^f} (I_3^f - 2Q_f\kappa_Z^fs_W^2)\,,
 * \qquad
 * g_A^f
 * &=&
 * \sqrt{\rho_Z^f} I_3^f\,,
 * @f}
 * and
 * @f{eqnarray}{
 * \rho_Z^f &=& \left( \frac{g_A^f}{I_3^f} \right)^2,
 * \qquad
 * \kappa_Z^f &=& \frac{1}{4|Q_f|s_W^2}
 * \left( 1 - \frac{g_V^{f}}{g_A^{f}}\right).
 * @f}
 *
 *
 * <h3>Important member functions</h3>
 *
 * The current class handles the following quantities:
 *
 * @li @f$\Delta\alpha_{\mathrm{lept}}(s)@f$&nbsp;&nbsp; (with DeltaAlphaLepton()),
 * @li @f$\Delta\alpha^{\ell+5q}(M_Z^2)@f$&nbsp;&nbsp; (with DeltaAlphaL5q()),
 * @li @f$\Delta\alpha_{\mathrm{top}}(s)@f$&nbsp;&nbsp; (with DeltaAlphaTop()),
 * @li @f$\Delta\alpha(M_Z^2)@f$&nbsp;&nbsp; (with DeltaAlpha()),
 * @li @f$\alpha(M_Z^2)@f$&nbsp;&nbsp; (with alphaMz()),
 *
 * @li @f$M_W@f$&nbsp;&nbsp; (with Mw_SM()),
 * @li @f$\Delta r@f$&nbsp;&nbsp; (with DeltaR_SM()),
 * @li @f$c_W^2@f$ and @f$s_W^2@f$&nbsp;&nbsp; (with cW2_SM() and sW2_SM()),
 * @li @f$\Gamma_W@f$&nbsp;&nbsp; (with GammaW_SM()),
 *
 * @li @f$\rho_Z^f@f$&nbsp;&nbsp; (with rhoZ_l() and rhoZ_q()),
 * @li @f$\kappa_Z^f@f$&nbsp;&nbsp; (with kappaZ_l() and kappaZ_q()),
 * @li @f$g_V^f@f$&nbsp;&nbsp; (with gVl() and gVq()),
 * @li @f$g_A^f@f$&nbsp;&nbsp; (with gAl() and gAq()),
 *
 * @li @f$\varepsilon_{1,2,3,b}@f$&nbsp;&nbsp; (with epsilon1_SM(), epsilon2_SM(),
 * epsilon3_SM() and epsilonb_SM()).
 *
 * Moreover, the functions Mzbar(), MwbarFromMw(), MwFromMwbar() and DeltaRbar_SM()
 * can be used for the quantities in the complex-pole/fixed-width scheme.
 *
 *
 * <h3>Schemes</h3>
 *
 * The formulae used for the @f$W@f$-boson mass @f$M_W@f$ and the effective
 * couplings @f$\rho_Z^f@f$ and @f$\kappa_Z^f@f$ are controlled with the model
 * flags @ref StandardModelFlags "Mw", @ref StandardModelFlags "RhoZ" and
 * @ref StandardModelFlags "KappaZ" of StandardModel. For each flag, the
 * available schemes are as follows:
 *
 * @li NORESUM:&nbsp;&nbsp; No resummation is considered;
 * @li OMSI:&nbsp;&nbsp; the so-called OMS-I scheme is adopted;
 * @li INTERMEDIATE:&nbsp;&nbsp; an intermediate scheme between OMS-I and OMS-II is adopted;
 * @li OMSII:&nbsp;&nbsp; the so-called OMS-II scheme is adopted;
 * @li APPROXIMATEFORMULA:&nbsp;&nbsp; the approximate two-loop formula given
 * in EWSMApproximateFormulae class is employed.
 *
 * The scheme APPROXIMATEFORMULA provides the most accurate SM predictions for
 * @f$M_W@f$ and @f$\kappa_Z^f@f$, while the approximate two-loop formula is
 * not available for @f$\rho_Z^f@f$.
 *
 * See resumMw(), resumRhoZ() and resumKappaZ() for details on the other schemes.
 *
 *
 * <h3>Cashes</h3>
 *
 * This class contains caching methods for the following functions:
 * DeltaAlphaLepton(), DeltaAlpha(), Mw_SM(), GammaW_SM(),
 * rhoZ_l_SM(), rhoZ_q_SM(), kappaZ_l_SM() and kappaZ_q_SM(), to improve the
 * performance of the Monte Carlo run. The caching methods are implemented
 * with the function checkSMparams().
 *
 * The use of the caching methods can be controlled with the model flag
 * @ref StandardModelFlags "CacheInStandardModel" of StandardModel.
 *
 */
class StandardModel : public QCD {
public:

    /**
     * @brief An enum type for leptons.
     */
    enum lepton {
        NEUTRINO_1, /**< The 1st-generation neutrino */
        ELECTRON, /**< Electron */
        NEUTRINO_2, /**< The 2nd-generation neutrino */
        MU, /**< Muon */
        NEUTRINO_3, /**< The 3rd-generation neutrino */
        TAU /**< Tau */
    };

    /**
     * @brief An enumerated type representing perturbative orders of radiative
     * corrections to %EW precision observables.
     */
    enum orders_EW {
        EW1 = 0, ///< One-loop of @f$\mathcal{O}(\alpha)@f$.
        EW1QCD1, ///< Two-loop of @f$\mathcal{O}(\alpha\alpha_s)@f$.
        EW1QCD2, ///< Three-loop of @f$\mathcal{O}(\alpha\alpha_s^2)@f$.
        EW2, ///< Two-loop of @f$\mathcal{O}(\alpha^2)@f$.
        EW2QCD1, ///< Three-loop of @f$\mathcal{O}(\alpha^2\alpha_s)@f$.
        EW3, ///< Three-loop of @f$\mathcal{O}(\alpha^3)@f$.
        orders_EW_size ///< The size of this enum.
    };

    static const int NSMvars = 26; ///< The number of the model parameters in %StandardModel.
    /**
     * @brief  A string array containing the labels of the model parameters in %StandardModel.
     */
    static const std::string SMvars[NSMvars];

    static const double GeVminus2_to_nb;

    /**
     * @brief The target accuracy of the iterative calculation of the
     * @f$W@f$-boson mass in units of GeV.
     */
    static const double Mw_error;

    /**
     * @brief The default constructor.
     */
    StandardModel();

    /**
     * @brief The default destructor.
     */
    virtual ~StandardModel();


    ///////////////////////////////////////////////////////////////////////////
    // Initialization

    /**
     * @brief A method to initialize the model.
     * @details This method, called via InputParser::ReadParameters(), allocates
     * memory to the pointers defined in the current class.
     * @return a boolean that is true if model initialization is successful
     */
    virtual bool InitializeModel();


    ///////////////////////////////////////////////////////////////////////////
    // Model parameters

    /**
     * @brief A method to initialize the model parameters.
     * @param[in] DPars a map of the parameters that are being updated in the Monte Carlo run
     * (including parameters that are varied and those that are held constant)
     * @return a boolean that is true if the execution is successful
     */
    virtual bool Init(const std::map<std::string, double>& DPars);

    /**
     * @brief The pre-update method for %StandardModel.
     * @details This method initializes the internal flags #requireCKM, #requireYe
     * and #requireYn, and calls QCD::PreUpdate(), before updating the model
     * parameters with the method Update().
     * @return a boolean that is true if the execution is successful
     */
    virtual bool PreUpdate();

    /**
     * @brief The update method for %StandardModel.
     * @details This method updates all the model parameters with given DPars.
     * @param[in] DPars a map of the parameters that are being updated in the Monte Carlo run
     * (including parameters that are varied and those that are held constant)
     * @return a boolean that is true if the execution is successful
     */
    virtual bool Update(const std::map<std::string, double>& DPars);

    /**
     * @brief The post-update method for %StandardModel.
     * @details This method runs all the procedures that are need to be executed
     * after the model is successfully updated. This includes 
     * \li computing the updated CKM and PMNS matrices
     * \li computing the Yukawa matrices
     * \li updating the Standard %Model parameters in the StandardModelMatching class.
     * @return a boolean that is true if the execution is successful
     */
    virtual bool PostUpdate();
    
    int getIterationNo() const
    {
        return iterationNo;
    }

    /**
     * @brief A method to check if all the mandatory parameters for %StandardModel
     * have been provided in model initialization.
     * @param[in] DPars a map of the parameters that are being updated in the Monte Carlo run
     * (including parameters that are varied and those that are held constant)
     * @return a boolean that is true if the execution is successful
     */
    virtual bool CheckParameters(const std::map<std::string, double>& DPars);


    ///////////////////////////////////////////////////////////////////////////
    // Flags

    /**
     * @brief A method to set a flag of %StandardModel.
     * @param[in] name name of a model flag
     * @param[in] value the boolean to be assigned to the flag specified by name
     * @return a boolean that is true if the execution is successful
     */
    virtual bool setFlag(const std::string name, const bool value);

    /**
     * @brief A method to set a flag of %StandardModel.
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

    /**
     * @brief A method to retrieve the model flag WithoutNonUniversalVC.
     * @details See @ref StandardModelFlags for detail. 
     * @return a boolean that is true if flavour non-universal vertex corrections
     * are NOT added to the epsilon parameters describing new physics contribution
     *
     * @attention The flag FlagWithoutNonUniversalVC is applicable only for
     * the models StandardModel and NPEpsilons.
     */
    bool IsFlagWithoutNonUniversalVC() const
    {
        return FlagWithoutNonUniversalVC;
    }

    /**
     * @brief A method to retrieve the model flag NoApproximateGammaZ.
     * @details See @ref StandardModelFlags for detail.
     * @return a boolean that is true if the two-loop approximate
     * formulae of the partial and total decay widths of the @f$Z@f$ boson
     * defined with the function EWSMApproximateFormulae::X_extended() is NOT
     * employed
     */
    bool IsFlagNoApproximateGammaZ() const
    {
        return FlagNoApproximateGammaZ;
    }

    void setFlagNoApproximateGammaZ(bool FlagNoApproximateGammaZ)
    {
        this->FlagNoApproximateGammaZ = FlagNoApproximateGammaZ;
    }

    /**
     * @brief A method to retrieve the model flag %Mw.
     * @details See @ref StandardModelFlags for detail.
     * @return
     */
    std::string getFlagMw() const
    {
        return FlagMw;
    }

    /**
     * @brief A method to retrieve the model flag %RhoZ.
     * @details See @ref StandardModelFlags for detail.
     * @return
     */
    std::string getFlagRhoZ() const
    {
        return FlagRhoZ;
    }

    /**
     * @brief A method to retrieve the model flag %KappaZ.
     * @details See @ref StandardModelFlags for detail.
     * @return
     */
    std::string getFlagKappaZ() const
    {
        return FlagKappaZ;
    }

    /**
     * @brief A set method to change
     * the model flag @ref StandardModelFlags "CacheInStandardModel" of StandardModel.
     * @details Setting CacheInStandardModel to false, the caching methods defined in the
     * current class are not employed in numerical computations. The flag is set
     * to true in the constructor EWSM() by default.
     * @param[in] FlagCacheInStandardModel true (false) if the caching methods are turned
     * on (off);
     *
     * @sa @ref StandardModelFlags "the description of the StandardModel flags"
     */
    void setFlagCacheInStandardModel(bool FlagCacheInStandardModel)
    {
        this->FlagCacheInStandardModel = FlagCacheInStandardModel;
    }


    ///////////////////////////////////////////////////////////////////////////
    // get and set methods for class members

    /**
     * @brief A get method to retrieve the member object of a lepton.
     * @param[in] p name of a lepton
     * @return an object of the lepton specified by name
     */
    Particle getLeptons(const StandardModel::lepton p) const
    {
        return leptons[p];
    }

    /**
     * @brief A get method to access the mass of the \f$Z\f$ boson \f$M_Z\f$.
     * @return the @f$Z@f$-boson mass @f$M_Z@f$
     */
    double getMz() const
    {
        return Mz;
    }

    /**
     * @brief A get method to access the value of \f$\alpha_s(M_Z)\f$.
     * @return the strong coupling constant at @f$M_Z@f$, @f$\alpha_s(M_Z)@f$
     */
    double getAlsMz() const
    {
        return AlsMz;
    }

    /**
     * @brief A get method to retrieve the Fermi constant @f$G_\mu@f$.
     * @return @f$G_\mu@f$ in @f${\rm GeV}^{-2}@f$
     */
    double getGF() const
    {
        return GF;
    }

    /**
     * @brief A get method to retrieve the fine-structure constant @f$\alpha@f$.
     * @return @f$\alpha@f$
     */
    double getAle() const
    {
        return ale;
    }

    /**
     * @brief A get method to retrieve the five-flavour hadronic contribution
     * to the electromagnetic coupling,
     * @f$\Delta\alpha_{\mathrm{had}}^{(5)}(M_Z^2)@f$.
     * @return @f$\Delta\alpha_{\mathrm{had}}^{(5)}(M_Z^2)@f$
     */
    double getDAle5Mz() const
    {
        return dAle5Mz;
    }

    /**
     * @brief A get method to retrieve the Higgs mass @f$m_h@f$.
     * @return @f$m_h@f$ in GeV
     */
    double getMHl() const
    {
        return mHl;
    }

    /**
     * @brief A get method to retrieve the theoretical uncertainty in @f$M_W@f$,
     * denoted as @f$\delta\,M_W@f$.
     * @return @f$\delta\,M_W@f$ in GeV
     */
    double getDelMw() const
    {
        return delMw;
    }

    /**
     * @brief A get method to retrieve the theoretical uncertainty in
     * @f$\sin^2\theta_{\rm eff}^{\rm lept}@f$, denoted as
     * @f$\delta\sin^2\theta_{\rm eff}^{\rm lept}@f$.
     * @return @f$\delta\sin^2\theta_{\rm eff}^{\rm lept}@f$
     */
    double getDelSin2th_l() const
    {
        return delSin2th_l;
    }

    /**
     * @brief A get method to retrieve the theoretical uncertainty in
     * @f$\Gamma_Z@f$, denoted as @f$\delta\,\Gamma_Z@f$.
     * @return @f$\delta\,\Gamma_Z@f$ in GeV 
     */
    double getDelGammaZ() const
    {
        return delGammaZ;
    }

    /**
     * @brief A get method to retrieve the %CKM matrix. 
     * @return the %CKM matrix
     */
    gslpp::matrix<gslpp::complex> getVCKM() const
    {
        return VCKM;
    }

    /**
     * @brief A get method to retrieve the member object of type CKM.
     * @return the object of type CKM
     */
    CKM getCKM() const
    {
        return myCKM;
    }

    /**
     * @brief A get method to retrieve the %CKM element @f$\lambda@f$.
     * @return @f$\lambda@f$
     */
    double getLambda() const
    {
        return lambda;
    }

    /**
     * @brief A get method to retrieve the %CKM element @f$A@f$.
     * @return @f$A@f$
     */
    double getA() const
    {
        return A;
    }

    /**
     * @brief A get method to retrieve the %CKM element @f$\bar{\rho}@f$.
     * @return @f$\bar{\rho}@f$
     */
    double getRhob() const
    {
        return rhob;
    }

    /**
     * @brief A get method to retrieve the %CKM element @f$\bar{\eta}@f$.
     * @return @f$\bar{\eta}@f$
     */
    double getEtab() const
    {
        return etab;
    }

    /**
     * @brief A get method to retrieve the object of the %PMNS matrix.
     * @return the %PMNS matrix
     */
    gslpp::matrix<gslpp::complex> getUPMNS() const
    {
        return UPMNS;
    }

    /**
     * @brief A get method to retrieve the Yukawa matrix of the up-type quarks, 
     * @f$Y_u@f$.
     * @return @f$Y_u@f$
     */
    gslpp::matrix<gslpp::complex> getYu() const
    {
        return Yu;
    }

    /**
     * @brief A get method to retrieve the Yukawa matrix of the down-type quarks,
     * @f$Y_d@f$.
     * @return @f$Y_d@f$
     */
    gslpp::matrix<gslpp::complex> getYd() const
    {
        return Yd;
    }

    /**
     * @brief A get method to retrieve the Yukawa matrix of the neutrinos,
     * @f$Y_\nu@f$.
     * @return @f$Y_\nu@f$
     */
    gslpp::matrix<gslpp::complex> getYn() const
    {
        return Yn;
    }

    /**
     * @brief A get method to retrieve the Yukawa matrix of the charged leptons,
     * @f$Y_e@f$.
     * @return @f$Y_e@f$
     */
    gslpp::matrix<gslpp::complex> getYe() const
    {
        return Ye;
    }

    /**
     * @brief A get method to retrieve the matching scale @f$\mu_W@f$ around
     * the weak scale.
     * @return @f$\mu_W@f$ in GeV
     */
    double getMuw() const
    {
        return muw;
    }

    double getEpsK() const
    {
        return EpsK;
    }

    double getphiEpsK() const
    {
        return phiEpsK;
    }

    double getKbarEpsK() const
    {
        return KbarEpsK;
    }

    double getDeltaMK() const
    {
        return DeltaMK;
    }

    double getDmk() const
    {
        return Dmk;
    }

    double getSM_M12D() const
    {
        return SM_M12D;
    }

    virtual StandardModel getTrueSM() const
    {
        throw std::runtime_error("StandardModel::getTrueSM() must be overridden by the NP extension.");
    }

    /**
     * @brief A get method to access the member pointer of type StandardModelMatching.
     * @return the pointer #myStandardModelMatching
     */
    virtual StandardModelMatching* getMyMatching() const
    {
        return myStandardModelMatching;
    }

    /**
     * @brief A get method to retrieve the member pointer of type EWSMcache.
     * @return the pointer #myEWSMcache
     */
    EWSMcache* getMyEWSMcache() const
    {
        return myEWSMcache;
    }

    /**
     * @brief A get method to retrieve the member pointer of type EWSMOneLoopEW,
     * @return the pointer #myOneLoopEW
     */
    EWSMOneLoopEW* getMyOneLoopEW() const
    {
        return myOneLoopEW;
    }

    /**
     * @brief A get method to retrieve the member pointer of type EWSMApproximateFormulae.
     * @return the pointer #myApproximateFormulae
     */
    EWSMApproximateFormulae* getMyApproximateFormulae() const
    {
        return myApproximateFormulae;
    }

    /** BEGIN: REMOVE FROM THE PACKAGE **/
    /**
     * @brief A get method to retrieve the member pointer of type EWSMTwoFermionsLEP2.
     * @return the pointer #myTwoFermionsLEP2
     */
    EWSMTwoFermionsLEP2* getMyTwoFermionsLEP2() const
    {
        return myTwoFermionsLEP2;
    }
    /** END: REMOVE FROM THE PACKAGE **/

    EWSMThreeLoopEW* getMyThreeLoopEW() const
    {
        return myThreeLoopEW;
    }

    EWSMThreeLoopEW2QCD* getMyThreeLoopEW2QCD() const
    {
        return myThreeLoopEW2QCD;
    }

    EWSMThreeLoopQCD* getMyThreeLoopQCD() const
    {
        return myThreeLoopQCD;
    }

    EWSMTwoLoopEW* getMyTwoLoopEW() const
    {
        return myTwoLoopEW;
    }

    EWSMTwoLoopQCD* getMyTwoLoopQCD() const
    {
        return myTwoLoopQCD;
    }

    Flavour* getMyFlavour() const
    {
        return myFlavour;
    }

    /** BEGIN: REMOVE FROM THE PACKAGE **/
    LeptonFlavour* getMyLeptonFlavour() const
    {
        return myLeptonFlavour;
    }
    /** END: REMOVE FROM THE PACKAGE **/


    ////////////////////////////////////////////////////////////////////////
    // CKM parameters

    /**
     * @brief The %CKM angle @f$\beta@f$.
     * @return @f$\beta@f$ in radians
     */
    double computeBeta() const;

    /**
     * @brief The %CKM angle @f$\gamma@f$.
     * @return @f$\gamma@f$ in radians
     */
    double computeGamma() const;

    /**
     * @brief The %CKM angle @f$\alpha@f$.
     * @return @f$\alpha@f$ in radians
     */
    double computeAlpha() const;

    /**
     * @brief The %CKM angle @f$\beta_s@f$.
     * @return @f$\beta_s@f$ in radians
     */
    double computeBetas() const;

    /**
     * @brief The product of the %CKM elements @f$V_{td} V_{ts}^*@f$.
     * @return @f$V_{td} V_{ts}^*@f$
     */
    gslpp::complex computelamt() const;

    /**
     * @brief The product of the %CKM elements @f$V_{cd} V_{cs}^*@f$.
     * @return @f$V_{cd} V_{cs}^*@f$
     */
    gslpp::complex computelamc() const;

    /**
     * @brief The product of the %CKM elements @f$V_{ud} V_{us}^*@f$.
     * @return @f$V_{ud} V_{us}^*@f$
     */
    gslpp::complex computelamu() const;

    /**
     * @brief The product of the %CKM elements @f$V_{td} V_{tb}^*@f$.
     * @return @f$V_{td} V_{tb}^*@f$
     */
    gslpp::complex computelamt_d() const;

    /**
     * @brief The product of the %CKM elements @f$V_{cd} V_{cb}^*@f$.
     * @return @f$V_{cd} V_{cb}^*@f$
     */
    gslpp::complex computelamc_d() const;

    /**
     * @brief The product of the %CKM elements @f$V_{ud} V_{ub}^*@f$.
     * @return @f$V_{ud} V_{ub}^*@f$
     */
    gslpp::complex computelamu_d() const;

    /**
     * @brief The product of the %CKM elements @f$V_{ts} V_{tb}^*@f$.
     * @return @f$V_{ts} V_{tb}^*@f$
     */
    gslpp::complex computelamt_s() const;

    /**
     * @brief The product of the %CKM elements @f$V_{cs} V_{cb}^*@f$.
     * @return @f$V_{cs} V_{cb}^*@f$
     */
    gslpp::complex computelamc_s() const;

    /**
     * @brief The product of the %CKM elements @f$V_{us} V_{ub}^*@f$.
     * @return @f$V_{us} V_{ub}^*@f$
     */
    gslpp::complex computelamu_s() const;

    /**
     * @brief @f$R_t=|(V_{td} V_{tb}^*)/(V_{cd}V_{cb}^*)|@f$.
     * @return @f$R_t@f$
     */
    double computeRt() const;

    /**
     * @brief @f$R_{ts}=|(V_{ts}V_{tb}^*)/(V_{cs}V_{cb}^*)|@f$.
     * @return @f$R_{ts}@f$
     */
    double computeRts() const;

    /**
     * @brief @f$R_b=|(V_{ud}V_{ub}^*)/(V_{ud}V_{ub}^*)|@f$.
     * @return @f$R_b@f$
     */
    double computeRb() const;


    ////////////////////////////////////////////////////////////////////////
    // QED coupling

    /**
     * @brief The running electromagnetic coupling @f$\alpha(\mu)@f$ in the
     * on-shell schem.
     * @details See @cite Baikov:2012rr.
     * @param[in] mu renormalization scale @f$\mu@f$ in GeV.
     * @param[in] order LO/FULLNLO
     * @return @f$\alpha(\mu)@f$ in the on-shell scheme
     *
     * @attention This function is applicable to the scale where the three charged
     * leptons and the five quarks, not the top quark, run in the loops.
     */
    double ale_OS(const double mu, orders order = FULLNLO) const;

    /**
     * @brief Leptonic contribution to the electromagnetic coupling @f$\alpha@f$,
     * denoted as @f$\Delta\alpha_{\mathrm{lept}}(s)@f$.
     * @param[in] s invariant mass squared
     * @return @f$\Delta\alpha_{\mathrm{lept}}(s)@f$
     */
    double DeltaAlphaLepton(const double s) const;

    /**
     * @brief The sum of the leptonic and the five-flavour hadronic corrections
     * to the electromagnetic coupling @f$\alpha@f$ at the @f$Z@f$-mass scale,
     * denoted as @f$\Delta\alpha^{\ell+5q}(M_Z^2)@f$.
     * @details
     * @f[
     * \Delta\alpha^{\ell+5q}(M_Z^2)
     * = \Delta\alpha_{\rm lept}(M_Z^2)
     * + \Delta\alpha_{\rm had}^{(5)}(M_Z^2)\,.
     * @f]
     * @return @f$\Delta\alpha^{\ell+5q}(M_Z^2)@f$
     */
    double DeltaAlphaL5q() const;

    /**
     * @brief Top-quark contribution to the electromagnetic coupling @f$\alpha@f$,
     * denoted as @f$\Delta\alpha_{\mathrm{top}}(s)@f$.
     * @param[in] s invariant mass squared
     * @return @f$\Delta\alpha_{\mathrm{top}}(s)@f$
     */
    double DeltaAlphaTop(const double s) const;

    /**
     * @brief The total corrections to the electromagnetic coupling @f$\alpha@f$
     * at the @f$Z@f$-mass scale, denoted as @f$\Delta\alpha(M_Z^2)@f$.
     * @details
     * @f[
     * \Delta\alpha(M_Z^2) =
     * \Delta\alpha_{\rm lept}(M_Z^2)
     * + \Delta\alpha_{\rm had}^{(5)}(M_Z^2)
     * + \Delta\alpha_{\rm top}(M_Z^2)\,.
     * @f]
     * @return @f$\Delta\alpha(M_Z^2)@f$
     */
    double DeltaAlpha() const;

    /**
     * @brief The electromagnetic coupling at the @f$Z@f$-mass scale, 
     * @f$\alpha(M_Z^2)=\alpha/(1-\Delta\alpha(M_Z^2))@f$.
     * @details The radiative corrections are included with Dyson resummation:
     * @f[
     * \alpha(M_Z^2) = \frac{\alpha}{1 - \Delta\alpha(M_Z^2)}.
     * @f]
     * @return @f$\alpha(M_Z^2)@f$
     */
    double alphaMz() const;


    ///////////////////////////////////////////////////////////////////////////
    // Higgs VEV

    /**
     * @brief The Higgs vacuum expectation value.
     * @f[
     *    v = \left(\frac{1}{\sqrt{2} G_\mu}\right)^{1/2},
     * @f]
     * where @f$G_\mu@f$ is the Fermi constant, measured through muon decays.
     * @return @f$v@f$ in GeV
     */
    virtual double v() const;


    ///////////////////////////////////////////////////////////////////////////
    // The W-boson mass

    /**
     * @brief The tree-level mass of the @f$W@f$ boson, @f$M_W^{\mathrm{tree}}@f$.
     * @return @f$M_W^{\mathrm{tree}}@f$ in GeV.
     */
    virtual double Mw_tree() const;

    /**
     * @brief The square of the sine of the weak mixing angle @f$s_0^2@f$ defined
     * without weak radiative corrections.
     * @details The quantity @f$s_0^2@f$ is defined through
     * @f[
     * s_0^2 c_0^2 = \frac{\pi\,\alpha(M_Z^2)}{\sqrt{2}\,G_\mu M_Z^2}
     * \ \ \rightarrow\ \
     * s_0^2 = \frac{1}{2}
     * \left(1 - \sqrt{1
     *  - \frac{4\pi \alpha(M_Z^2)}{\sqrt{2}\,G_\mu M_Z^2}}\  \right)\,.
     * @f]
     *
     * See @cite Altarelli:1990zd and @cite Altarelli:1991fk.
     * @return @f$s_0^2@f$
     */
    double s02() const;

    /**
     * @brief The square of the cosine of the weak mixing angle @f$c_0^2@f$ defined
     * without weak radiative corrections.
     * @details The quantity @f$c_0^2@f$ is given by
     * @f[
     * c_0^2 = 1 - s_0^2\,,
     * @f]
     * where @f$s_0^2@f$ is defined in s02().
     *
     * See @cite Altarelli:1990zd and @cite Altarelli:1991fk.
     * @return @f$s_0^2@f$
     */
    double c02() const;

    /**
     * @brief The SM prediction for the @f$W@f$-boson mass in the on-shell scheme,
     * @f$M_{W,\mathrm{SM}}@f$.
     * @details
     * When the model flag @ref StandardModelFlags "Mw" of StandardModel is set
     * to APPROXIMATEFORMULA, the current function uses
     * the two-loop approximate formula in EWSMApproximateFormulae::Mw(),
     * which includes the full two-loop %EW contribution of
     * @f${\cal O}(\alpha^2)@f$ as well as the leading
     * @f${\cal O}(G_\mu^2\alpha_s m_t^4)@f$ and @f${\cal O}(G_\mu^3m_t^6)@f$
     * contributions.
     *
     * When the model flag @ref StandardModelFlags "Mw" is not set to APPROXIMATEFORMULA,
     * the @f$W@f$-boson mass is computed from @f$\Delta r(M_W)@f$ with an
     * iterative procedure. The target accuracy of the iterative calculation
     * is specified with the constant #Mw_error.
     * This function calls resumMw(), in which @f$M_W@f$ is computed with a given
     * @f$\Delta r@f$, equivalently with @f$\Delta\rho@f$ and
     * @f$\Delta r_{\mathrm{rem}}@f$
     * @return @f$M_{W,\mathrm{SM}}@f$ in GeV
     *
     * @sa resumMw()
     * @attention If the model flag @ref StandardModelFlags "CacheInStandardModel"
     * of StandardModel is set to true, the caching method implemented in the
     * current class is employed.
     */
    virtual double Mw() const;

    /**
     * @brief The square of the cosine of the weak mixing angle
     * in the on-shell scheme, denoted as @f$c_W^2@f$.
     * @details
     * @f[
     *   c_W^2=\cos^2{\theta_W}=\frac{M_W^2}{M_Z^2}.
     * @f]
     * @return @f$c_W^2@f$
     */
    virtual double cW2(const double Mw_i) const;
    virtual double cW2() const;

    /**
     * @brief The square of the sine of the weak mixing angle
     * in the on-shell scheme, denoted as @f$s_W^2@f$.
     * @details
     * @f[
     *   s_W^2=\sin^2{\theta_W}=1-\frac{M_W^2}{M_Z^2}.
     * @f]
     * @return @f$s_W^2@f$
     */
    virtual double sW2(const double Mw_i) const;
    double sW2() const;

    /**
     * @brief The SM prediction for @f$\Delta r@f$ derived from that for the
     * @f$W@f$ boson mass.
     * @details
     * If the model flag @ref StandardModelFlags "Mw" of StandardModel is set
     * to NORESUM or APPROXIMATEFORMULA, the quantity @f$\Delta r@f$ is computed
     * by using the following relation:
     * @f[
     * s_W^2 M_W^2 = \frac{\pi\,\alpha}{\sqrt{2}G_\mu}(1+\Delta r)\,.
     * @f]
     * Otherwise, the following relation is employed instead:
     * @f[
     * s_W^2 M_W^2 = \frac{\pi\,\alpha}{\sqrt{2}G_\mu(1-\Delta r)}\,,
     * @f]
     * where the resummation for @f$\Delta r@f$ is considered.
     * @return @f$\Delta r_{\mathrm{SM}}@f$
     *
     * @sa The corresponding quantity in the complex-pole/fixed-width scheme
     * (instead of the experimental/running-widthr scheme)
     * is defined in DeltaRbar_SM().
     */
    virtual double DeltaR() const;

    /**
     * @brief A method to collect @f$\Delta\rho@f$ computed via subclasses.
     * @details This function collects @f$\Delta\rho@f$
     * computed via EWSMOneLoopEW, EWSMTwoLoopQCD, EWSMTwoLoopEW, EWSMThreeLoopQCD,
     * EWSMThreeLoopEW2QCD and EWSMThreeLoopEW classes.
     * @param[in] Mw_i the @f$W@f$-boson mass
     * @param[out] DeltaRho Array of @f$\Delta\rho@f$
     */
    void ComputeDeltaRho(const double Mw_i, double DeltaRho[orders_EW_size]) const;

    /**
     * @brief A method to collect @f$\Delta r_{\mathrm{rem}}@f$ computed via subclasses.
     * @details This function collects @f$\Delta r_{\mathrm{rem}}@f$
     * computed via EWSMOneLoopEW, EWSMTwoLoopQCD, EWSMTwoLoopEW, EWSMThreeLoopQCD,
     * EWSMThreeLoopEW2QCD and EWSMThreeLoopEW classes.
     * @param[in] Mw_i the @f$W@f$-boson mass
     * @param[out] DeltaR_rem Array of @f$\Delta r_{\mathrm{rem}}@f$
     */
    void ComputeDeltaR_rem(const double Mw_i, double DeltaR_rem[orders_EW_size]) const;


    ////////////////////////////////////////////////////////////////////////
    // The W and Z masses in the complex-pole/fixed-width scheme

    /**
     * @brief The @f$Z@f$-boson mass @f$\overline{M}_Z@f$
     * in the complex-pole/fixed-width scheme.
     * @details The mass parameter @f$\overline{M}_Z@f$ in the
     * complex-pole/fixed-width scheme @cite Bardin:1988xt is given by
     * @f[
     *   \overline{M}_{Z} = M_{Z} - \frac{\Gamma_{Z}^2}{2M_{Z}}\,,
     * @f]
     * where @f$M_Z@f$ and @f$\Gamma_{Z}@f$ are the mass and width of the
     * @f$Z@f$ boson in the experimental/running-width scheme:
     * @f{align}{
     * \Gamma(Z\to f\bar{f})
     * = \frac{G_\mu M_Z^3}{24\sqrt{2}\pi}
     * \left[
     * \left( \frac{v_f}{a_f} \right)^2 + 1
     * \right]
     * \times
     * \left\{
     * \begin{array}{ll}
     * 1 & \mathrm{for}\quad f=\ell\,,
     * \\[2mm]
     * \displaystyle
     * N_c \left( 1 + \frac{\alpha_s(M_Z^2)}{\pi} \right)
     * & \mathrm{for}\quad f=q
     * \end{array}
     * \right.
     * @f}
     * with @f$v_f/a_f=1-4|Q_f|s_{W,\mathrm{tree}}^2@f$.
     * @return @f$\overline{M}_Z@f$ in GeV
     */
    double Mzbar() const;

    /**
     * @brief A method to convert the @f$W@f$-boson mass
     * in the experimental/running-width scheme
     * to that in the complex-pole/fixed-width scheme.
     * @details The mass parameter @f$\overline{M}_W@f$ in the
     * complex-pole/fixed-width scheme @cite Bardin:1988xt is given by
     * @f[
     *   \overline{M}_{W} = M_{W} - \frac{\Gamma_{W}^2}{2M_{W}}\,,
     * @f]
     * where @f$M_W@f$ and @f$\Gamma_{W}@f$ are the mass and width of the
     * @f$W@f$ boson in the experimental/running-width scheme:
     * @f[
     * \Gamma_W
     * =
     *   \frac{3G_\mu M_W^3}{2\sqrt{2}\pi}
     *   \left( 1 + \frac{2\alpha_s(M_W^2)}{3\pi} \right)\,.
     * @f]
     * @param[in] Mw the @f$W@f$-boson mass in the experimental/running-width scheme
     * @return @f$\overline{M}_W@f$ in GeV
     */
    double MwbarFromMw(const double Mw) const;

    /**
     * @brief A method to convert the @f$W@f$-boson mass
     * in the complex-pole/fixed-width scheme
     * to that in the experimental/running-width scheme.
     * @details The experimental mass @f$M_W@f$ is derived
     * @f[
     *   M_W = \overline{M}_W + \frac{\Gamma_{W}^2}{2\overline{M}_{W}}\,,
     * @f]
     * where @f$\overline{M}_W@f$ is the mass parameter in the
     * complex-pole/fixed-width scheme @cite Bardin:1988xt, and @f$\Gamma_{W}@f$
     * is the @f$W@f$-boson width in the experimental/running-width scheme:
     * @f[
     * \Gamma_W
     * =
     *   \frac{3G_\mu M_W^3}{2\sqrt{2}\pi}
     *   \left( 1 + \frac{2\alpha_s(M_W^2)}{3\pi} \right)
     * \approx
     *   \frac{3G_\mu \overline{M}_W^3}{2\sqrt{2}\pi}
     *   \left( 1 + \frac{2\alpha_s(\overline{M}_W^2)}{3\pi} \right)\,.
     * @f]
     * @param[in] Mwbar the @f$W@f$-boson mass in the complex-pole/fixed-width scheme
     * @return @f$M_W@f$ in GeV
     */
    double MwFromMwbar(const double Mwbar) const;

    /**
     * @brief The SM prediction for @f$\Delta \overline{r}@f$ derived from
     * that for the @f$W@f$-boson mass.
     * @details The quantity @f$\Delta \overline{r}@f$ is computed by using
     * the following relation:
     * @f[
     * \overline{s}_W^2 \overline{M}_W^2
     * = \frac{\pi\,\alpha}{\sqrt{2}G_\mu}(1+\Delta \overline{r})\,,
     * @f]
     * where @f$\overline{M}_W@f$ and @f$\overline{s}_W@f$ are the @f$W@f$-boson
     * mass and the sine of the weak mixing angle in the complex-pole/fixed-width
     * scheme @cite Bardin:1988xt.
     * @return @f$\Delta \overline{r}_{\mathrm{SM}}@f$
     *
     * @sa DeltaR_SM(), defining the corresponding quantity in the
     * experimental/running-width scheme.
     */
    virtual double DeltaRbar() const;


    ////////////////////////////////////////////////////////////////////////
    // The W-boson decay width

    /**
     * @brief @copybrief EWSMOneLoopEW::rho_GammaW()
     * @param[in] fi a lepton or quark
     * @param[in] fj a lepton or quark
     * @return @f$\rho^W_{ij}@f$
     *
     * @sa EWSMOneLoopEW::rho_GammaW()
     */
    virtual double rho_GammaW(const Particle fi, const Particle fj) const;

    /**
     * @brief A partial decay width of the @f$W@f$ boson decay into a SM fermion pair.
     * @details
     * @f[
     * \Gamma^W_{ij}
     * =
     * |U_{ij}|^2\,\frac{G_\mu M_W^3}{6\sqrt{2}\,\pi}\,\rho^W_{ij}
     * @f]
     * where @f$U@f$ denotes the %MNS matrix, and @f$\rho^W_{ij}@f$ represents
     * %EW radiative corrections.
     *
     * @f[
     * \Gamma^W_{ij}
     * =
     * 3 |V_{ij}|^2\,\frac{G_\mu M_W^3}{6\sqrt{2}\,\pi}\,\rho^W_{ij}
     * \left( 1 + \frac{\alpha_s(M_W^2)}{\pi} \right).
     * @f]
     * where @f$V@f$ denotes the %CKM matrix, and @f$\rho^W_{ij}@f$ represents
     * %EW radiative corrections.
     * @param[in] fi a lepton or quark
     * @param[in] fj a lepton or quark
     * @return @f$\Gamma^W_{ij}@f$
     *
     * @sa rho_GammaW_l_SM()
     * @attention Fermion masses are neglected.
     */
    virtual double GammaW(const Particle fi, const Particle fj) const;

    /**
     * @brief The total width of the @f$W@f$ boson, @f$\Gamma_W@f$.
     * @return @f$\Gamma_W@f$ in GeV
     */
    virtual double GammaW() const;


    ////////////////////////////////////////////////////////////////////////     
    // EWPO at Z-pole

    /**
     * @brief The left-right asymmetry in @f$e^+e^-\to Z\to \ell \bar{\ell}@f$ at the
     * @f$Z@f$-pole, @f$\mathcal{A}_\ell@f$.
     * @details The asymmetry  @f$\mathcal{A}_\ell@f$ is given by
     * @f[
     * \mathcal{A}_\ell =
     * \frac{2\, {\rm Re}\left(g_{V}^\ell/g_{A}^\ell\right)}
     * {1+\left[{\rm Re}\left(g_{V}^\ell/g_{A}^\ell\right)\right]^2}\,,
     * @f]
     * where the ratio of the effective couplings @f$g_{V}^\ell/g_{A}^\ell@f$ is
     * computed via the two-loop approximate formula of 
     * @f$\sin^2\theta_{\rm eff}^{\,\ell}@f$, EWSMApproximateFormulae::sin2thetaEff_l(),
     * when checkNPZff_linearized() returns true and
     * the model flag @ref StandardModelFlags "KappaZ" of StandardModel
     * is set to APPROXIMATEFORMULA.
     * @param[in] f a lepton or quark
     * @return @f$\mathcal{A}_\ell@f$
     */
    virtual double A_f(const Particle f) const;

    /**
     *
     * @param[in] f a lepton or quark
     * @return
     */
    virtual double AFB(const Particle f) const;

    /**
     * @brief The effective weak mixing angle @f$\sin^2\theta_{\rm eff}^{\,\ell}@f$
     * for @f$Z\ell\bar{\ell}@f$ at the the @f$Z@f$-mass scale.
     * @details
     * When checkNPZff_linearized() returns true and
     * the model flag @ref StandardModelFlags "KappaZ" of StandardModel 
     * is set to APPROXIMATEFORMULA, this function uses the two-loop approximate
     * formula of @f$\sin^2\theta_{\rm eff}^{\,\ell}@f$ via
     * EWSMApproximateFormulae::sin2thetaEff().
     * Otherwise, the effective weak mixing angle is calculated from the coupling
     * @f$\kappa_Z^\ell@f$:
     * @f[
     * \sin^2\theta_{\rm eff}^{\,\ell} = {\rm Re}(\kappa_Z^\ell)\,s_W^2\,.
     * @f]
     * @param[in] f a lepton or quark
     * @return @f$\sin^2\theta_{\rm eff}^{\,\ell}@f$
     *
     * @attention @f$\ell@f$ stands for both a neutrino and a charged lepton. 
     */
    virtual double sin2thetaEff(const Particle f) const;

    /**
     * @brief The @f$Z\to \ell\bar{\ell}@f$ partial decay width, @f$\Gamma_\ell@f$.
     * @details
     * When checkNPZff_linearized() returns true and the model flag
     * @ref StandardModelFlags "NoApproximateGammaZ" of StandardModel is set
     * to false, this function uses the two-loop approximate formula of
     * @f$\Gamma_\ell@f$ via EWSMApproximateFormulae::X_extended().
     * Otherwise, the partial width is calculated with
     * @f$\rho_Z^\ell@f$ and @f$g_{V}^\ell/g_{A}^\ell@f$ @cite Bardin:1999ak :
     * @f[
     * \Gamma_\ell =
     * \Gamma_0 \big|\rho_Z^f\big|
     * \sqrt{1-\frac{4m_\ell^2}{M_Z^2}}
     * \left[ \left(1+\frac{2m_\ell^2}{M_Z^2}\right)
     *   \left(\left|\frac{g_{V}^\ell}{g_{A}^\ell}\right|^2 + 1 \right)
     *   - \frac{6m_\ell^2}{M_Z^2}
     * \right]
     * \left( 1 + \frac{3}{4}\frac{\alpha(M_Z^2)}{\pi}\, Q_\ell^2 \right)
     * @f]
     * with @f$\Gamma_0=G_\mu M_Z^3/(24\sqrt{2}\pi)@f$. 
     * @param[in] f a lepton or quark
     * @return @f$\Gamma_\ell@f$ in GeV
     *
     * @attention @f$\ell@f$ stands for both a neutrino and a charged lepton. 
     */
    virtual double GammaZ(const Particle f) const;

    /**
     * @brief The invisible partial decay width of the @f$Z@f$ boson,
     * @f$\Gamma_{\mathrm{inv}}@f$.
     * @details
     * @f[
     * \Gamma_{\mathrm{inv}} = 3\,\Gamma_\nu\,,
     * @f]
     * where @f$\Gamma_{\nu}@f$ is the partial width for @f$Z\to\nu\bar{\nu}@f$.
     * @return @f$\Gamma_{\mathrm{inv}}@f$ in GeV
     */
    virtual double Gamma_inv() const;

    /**
     * @brief The hadronic decay width of the @f$Z@f$ boson, @f$\Gamma_{h}@f$.
     * @details
     * The hadronic width is given by the sum, 
     * @f[
     *  \Gamma_h = \Gamma_u + \Gamma_d + \Gamma_c + \Gamma_s + \Gamma_b\,.
     * @f]
     * Furthermore, the singlet vector corrections are added, following the
     * prescription in @cite Bardin:1997xq :
     * @f[
     * \Gamma_h = \sum_q \Gamma_q + 4N_c\Gamma_0 R_V^h\,.
     * @f]
     * @return @f$\Gamma_{h}@f$ in GeV
     */
    virtual double Gamma_had() const;

    /**
     * @brief The total decay width of the @f$Z@f$ boson, @f$\Gamma_Z@f$.
     * @details When checkNPZff_linearized() returns true and the model flag
     * @ref StandardModelFlags "NoApproximateGammaZ" of StandardModel is set
     * to false, this function uses the two-loop approximate formula of
     * @f$\Gamma_Z@f$ via EWSMApproximateFormulae::X_extended().
     * Otherwise, the total decay width is calculated with
     * @f[
     * \Gamma_Z = \Gamma_{e} + \Gamma_{\mu} + \Gamma_{\tau} 
     * + \Gamma_{\mathrm{inv}} + \Gamma_h\,.
     * @f]
     * @return @f$\Gamma_Z@f$ in GeV
     */
    virtual double Gamma_Z() const;

    /**
     * @brief The hadronic cross section for @f$e^+e^- \to Z \to \mathrm{hadrons}@f$
     * at the @f$Z@f$-pole, @f$\sigma_h^0@f$.
     * @details When checkNPZff_linearized() returns true and the model flag
     * @ref StandardModelFlags "NoApproximateGammaZ" of StandardModel is set
     * to false, this function uses the two-loop approximate formula of
     * @f$\sigma_h^0@f$ via EWSMApproximateFormulae::X_extended().
     * Otherwise, the hadronic cross section is calculated with
     * @f[
     * \sigma_h^0 = \frac{12\pi}{M_Z^2}\frac{\Gamma_e\Gamma_h}{\Gamma_Z^2}\,.
     * @f]
     * @return @f$\sigma_h^0@f$ in GeV@f$^{-2}@f$
     */
    virtual double sigma0_had() const;

    /**
     * @brief @copybrief Rlepton::computeThValue()
     * @details When checkNPZff_linearized() returns true and the model flag
     * @ref StandardModelFlags "NoApproximateGammaZ" of StandardModel is set
     * to false, this function uses the two-loop approximate formula of
     * @f$R_\ell^0@f$ via EWSMApproximateFormulae::X_extended().
     * Otherwise, @f$R_\ell^0@f$ is calculated with
     * @f[
     * R_\ell^0 = \frac{\Gamma_h}{\Gamma_\ell}\,.
     * @f],
     * where @f$\ell@f$ denotes a charged lepton.
     * @param[in] f a lepton or quark
     * @return @f$R_\ell^0 @f$
     */
    virtual double R0_f(const Particle f) const;


    ////////////////////////////////////////////////////////////////////////
    // Zff effective couplings

    /**
     * @brief The effective leptonic neutral-current vector coupling @f$g_V^l@f$ in the SM.
     * @details
     * @f[
     * g_V^l = g_A^l (1 - 4|Q_l|\kappa_Z^l s_W^2)\,.
     * @f]
     * @param[in] f a lepton or quark
     * @return @f$g_{V,\,\mathrm{SM}}^l@f$
     */
    virtual gslpp::complex gV_f(const Particle f) const;

    /**
     * @brief The effective leptonic neutral-current axial-vector coupling @f$g_A^l@f$ in the SM.
     * @details
     * @f[
     * g_A^l = \sqrt{\rho_Z^l}\, I_3^l\,.
     * @f]
     * @param[in] f a lepton or quark
     * @return @f$g_{A,\,\mathrm{SM}}^l@f$
     */
    virtual gslpp::complex gA_f(const Particle f) const;

    /**
     * @brief The effective leptonic neutral-current coupling @f$\rho_Z^l@f$ in the SM.
     * @details This function collects the radiative corrections to @f$\rho_Z^l@f$
     * computed via EWSMOneLoopEW, EWSMTwoLoopQCD, EWSMTwoLoopEW, EWSMThreeLoopQCD,
     * EWSMThreeLoopEW2QCD and EWSMThreeLoopEW classes. The real part is computed
     * with the function resumRhoZ(), while only the one-loop contribution is kept
     * in the imaginary part. 
     * @param[in] f a lepton or quark
     * @return @f$\rho_{Z,\,\mathrm{SM}}^l@f$
     * 
     * @sa resumRhoZ()
     * @attention If the model flag @ref StandardModelFlags "CacheInStandardModel"
     * of StandardModel is set to true, the caching method implemented in the
     * current class is employed.
     */
    virtual gslpp::complex rhoZ_f(const Particle f) const;

    /**
     * @brief The effective leptonic neutral-current coupling @f$\kappa_Z^l@f$ in the SM.
     * @details This function collects the radiative corrections to @f$\kappa_Z^l@f$
     * computed via EWSMOneLoopEW, EWSMTwoLoopQCD, EWSMTwoLoopEW, EWSMThreeLoopQCD,
     * EWSMThreeLoopEW2QCD and EWSMThreeLoopEW classes. The real part is computed
     * with the function resumKappaZ(), while only the one-loop contribution is kept
     * in the imaginary part.
     *
     * As a part of the two-loop %EW contribution, a correction associated with
     * the product of the imaginary part of @f$\Delta\alpha@f$ and that of
     * @f$\Pi_{Z\gamma}@f$ is included @cite Bardin:1999ak, @cite Bardin:1999yd :
     * @f{eqnarray}{
     * \Delta \kappa_Z^l
     * = - \frac{1}{s_W^2}\left( \frac{\alpha(M_Z^2)}{4\pi} \right)^2
     * {\rm Im}\,\overline{\Pi}_{\gamma\gamma}^{\rm fer}(M_Z^2)\,\,
     * {\rm Im}\,\overline{\Pi}_{Z\gamma}^{\rm fer}(M_Z^2)
     * = \frac{35\alpha^2(M_Z^2)}{18 s_W^2}\,
     * \left( 1 - \frac{8}{3}\, {\rm Re}(\kappa_Z^l) s_W^2 \right).
     * @f}
     * @param[in] f a lepton or quark
     * @return @f$\kappa_{Z,\,\mathrm{SM}}^l@f$
     *
     * @sa resumKappaZ()
     * @attention If the model flag @ref StandardModelFlags "CacheInStandardModel"
     * of StandardModel is set to true, the caching method implemented in the
     * current class is employed.
     */
    virtual gslpp::complex kappaZ_f(const Particle f) const;

    /**
     * @brief Flavour non-universal vertex corrections to @f$\rho_Z^l@f$,
     * denoted by @f$\Delta\rho_Z^l@f$.
     * @details The non-universal contribution @f$\Delta\rho_Z^l@f$ is given by
     * @f[
     * \Delta \rho_Z^l = \rho_Z^l - \rho_Z^e
     * = \frac{\alpha}{2\pi s_W^2}\left(u_l - u_e\right),
     * @f]
     * where @f$u_l@f$ is defined as
     * @f[
     * u_l = \frac{3v_l^2+a_l^2}{4c_W^2}\mathcal{F}_Z(M_Z^2) + \mathcal{F}_W^l(M_Z^2)
     * @f]
     * with the tree-level vector and axial-vector couplings
     * @f$v_l = I_3^l - 2Q_l s_W^2@f$ and @f$a_l = I_3^l@f$ and the form factors,
     * @f$\mathcal{F}_Z@f$ and @f$\mathcal{F}_W^l@f$.
     *
     * See @cite Ciuchini:2013pca and references therein. 
     * @param[in] f a lepton or quark
     * @return @f$\Delta\rho_Z^l@f$
     */
    virtual gslpp::complex deltaRhoZ_f(const Particle f) const;

    /**
     * @brief Flavour non-universal vertex corrections to @f$\kappa_Z^l@f$,
     * denoted by @f$\Delta\kappa_Z^l@f$.
     * @details The non-universal contribution @f$\Delta\kappa_Z^l@f$ is given by
     * @f[
     * \Delta \kappa_Z^l = \kappa_Z^l - \kappa_Z^e
     * = \frac{\alpha}{4\pi s_W^2}
     *  \left( \frac{\delta_l^2-\delta_e^2}{4c_W^2}\,\mathcal{F}_Z(M_Z^2)
     *  -u_l+u_e\right),
     * @f]
     * where @f$u_l@f$ and @f$\delta_l@f$ are defined as
     * @f[
     * u_l = \frac{3v_l^2+a_l^2}{4c_W^2}\mathcal{F}_Z(M_Z^2) + \mathcal{F}_W^l(M_Z^2)\,,
     * \qquad
     * \delta_l = v_l - a_l
     * @f]
     * with the tree-level vector and axial-vector couplings
     * @f$v_l = I_3^l - 2Q_l s_W^2@f$ and @f$a_l = I_3^l@f$, and the form factors
     * @f$\mathcal{F}_Z@f$ and @f$\mathcal{F}_W^l@f$.
     *
     * See @cite Ciuchini:2013pca and references therein.
     * @param[in] f a lepton or quark
     * @return @f$\Delta\kappa_Z^l@f$
     */
    virtual gslpp::complex deltaKappaZ_f(const Particle f) const;


    ////////////////////////////////////////////////////////////////////////
    // Epsilon parameters for EWPO

    /**
     * @brief The SM contribution to the epsilon parameter @f$\varepsilon_1@f$.
     * @details The parameters @f$\varepsilon_1@f$ is defined as
     * @f[
     * \varepsilon_1 = \Delta\rho'\,,
     * @f]
     * where @f$\Delta\rho'=2\left(\sqrt{{\rm Re}(\rho_Z^e)}-1\right)@f$.
     *
     * See @cite Altarelli:1990zd and @cite Altarelli:1991fk.
     * @return @f$\varepsilon_{1,\mathrm{SM}}@f$
     */
    virtual double epsilon1() const;

    /**
     * @brief The SM contribution to the epsilon parameter @f$\varepsilon_2@f$.
     * @details The parameters @f$\varepsilon_2@f$ is computed via the formula: 
     * @f[
     * \varepsilon_2 = c_0^2  \Delta\rho'
     * + \frac{s_0^2}{c_0^2 - s_0^2} \Delta r_W
     * - 2 s_0^2 \Delta\kappa'\,,
     * @f]
     * where @f$\Delta\rho'@f$, @f$\Delta r_W@f$ and @f$\Delta\kappa'@f$ are 
     * defined as 
     * @f{align}{
     * \Delta\rho'=2\left(\sqrt{{\rm Re}(\rho_Z^e)}-1\right),\qquad
     * \Delta r_W = 1 - \frac{\pi\,\alpha(M_Z^2)}{\sqrt{2}\,G_\mu M_Z^2 s_W^2 c_W^2},\qquad
     * \Delta\kappa' = \frac{\sin^2\theta_{\mathrm{eff}}^e}{s_0^2} - 1\,,
     * @f}
     * and @f$s_0^2@f$ and @f$c_0^2@f$ are given in s02() and c02(), respectively.
     *
     * See @cite Altarelli:1990zd and @cite Altarelli:1991fk.
     * @return @f$\varepsilon_{2,\mathrm{SM}}@f$
     */
    virtual double epsilon2() const;

    /**
     * @brief The SM contribution to the epsilon parameter @f$\varepsilon_3@f$.
     * @details The parameters @f$\varepsilon_3@f$ is computed via the formula:
     * @f[
     * \varepsilon_3 = c_0^2\Delta\rho' + (c_0^2-s_0^2)\Delta\kappa'\,,
     * @f]
     * where @f$\Delta\rho'@f$ and @f$\Delta\kappa'@f$ are
     * defined as
     * @f{align}{
     * \Delta\rho'=2\left(\sqrt{{\rm Re}(\rho_Z^e)}-1\right),\qquad
     * \Delta\kappa' = \frac{\sin^2\theta_{\mathrm{eff}}^e}{s_0^2} - 1\,,
     * @f}
     * and @f$s_0^2@f$ and @f$c_0^2@f$ are given in s02() and c02(), respectively.
     *
     * See @cite Altarelli:1990zd and @cite Altarelli:1991fk.
     * @return @f$\varepsilon_{3,\mathrm{SM}}@f$
     */
    virtual double epsilon3() const;

    /**
     * @brief The SM contribution to the epsilon parameter @f$\varepsilon_b@f$.
     * @details The parameters @f$\varepsilon_b@f$ is computed via the formula:
     * @f[
     * \epsilon_b = 
     * \frac{ {\rm Re}\left[ \kappa_Z^e + \Delta\kappa_Z^b \right]}
     * {{\rm Re}(\kappa_Z^b)} - 1\,,
     * @f]
     * where @f$\Delta\kappa_Z^b@f$, representing flavour non-universal vertex
     * corrections to the @f$Zb\bar{b}@f$ vertex, is neglected when the
     * model flag WithoutNonUniversalVC of StandardModel is set to true.
     * 
     * See @cite Altarelli:1990zd, @cite Altarelli:1991fk and @cite Altarelli:1993sz
     * for the @f$\varepsilon@f$ parameterization
     * and @cite Ciuchini:2013pca for the flavour non-universal vertex corrections. 
     * @return @f$\varepsilon_{b,\mathrm{SM}}@f$
     */
    virtual double epsilonb() const;


    ////////////////////////////////////////////////////////////////////////
    // For EWPO caches

    /**
     * @brief The number of the SM parameters that are relevant to the %EW
     * precision observables.
     * @details This constant is used for the cashing method.
     *
     * @sa checkSMparamsForEWPO()
     */
    static const int NumSMParamsForEWPO = 27;

    /**
     * @brief A method to check whether the parameters relevant to the EWPO
     * are updated. 
     * @details This function is used for the cashing methods implemented in
     * the current class:
     * DeltaAlphaLepton(), DeltaAlpha(), Mw_SM(), rhoZ_l_SM(), rhoZ_q_SM(),
     * kappaZ_l_SM(), kappaZ_q_SM() and GammaW_SM().
     * When the values of the StandardModel parameters are updated in the Monte
     * Carlo run and differ from those stored in the cache SMparamsForEWPO_cache,
     * this function updates the cache, and returns false.
     * @return a boolean that is true if the parameters are not updated. 
     *
     * @sa NumSMParamsForEWPO
     */
    bool checkSMparamsForEWPO();

    ////////////////////////////////////////////////////////////////////////        
    // Several Higgs-related quantities used in Higgs coupling analysis

    /**
     * @brief The ggH cross section in the Standard Model.
     * @details See Tables B.67 and B.74 in ref. @cite Heinemeyer:2013tqa
     * @param[in] sqrt_s the center-of-mass energy in TeV
     * @return ggH cross section in pb
     */
    double computeSigmaggH(const double sqrt_s) const
    {
        if (sqrt_s == 7.0) {
            return 15.11; // in pb for Mh=125.1 GeV
            //return 14.99; // in pb for Mh=125.6 GeV
        } else if (sqrt_s == 8.0) {
            return 19.24; // in pb for Mh=125.1 GeV
            //return 19.09; // in pb for Mh=125.6 GeV
        } else if (sqrt_s == 1.96) {
            return 0.9493; // in pb for Mh=125 GeV
        } else
            throw std::runtime_error("Bad argument in StandardModel::computeSigmaggH()");
    }

    /**
     * @brief The square of the top-quark contribution to the ggH cross section in the Standard Model.
     * @details The values have been obtained from  M. Spira. See also Table 38
     * in ref. @cite Heinemeyer:2013tqa, which are calculated with a scale choice of Mh.
     * @param[in] sqrt_s the center-of-mass energy in TeV
     * @return @f$\sigma_{ggH}^{tt}@f$ in pb
     */
    double computeSigmaggH_tt(const double sqrt_s) const
    {
        if (sqrt_s == 7.0) {
            //return 14.355; // in pb for Mh=125 GeV with a scale choice of Mh (corrected from the value in 1307.1347v2)
            //return 15.89; // in pb for Mh=125 GeV with a scale choice of Mh/2
            return 15.76; // in pb for Mh=125.5 GeV with a scale choice of Mh/2
        } else if (sqrt_s == 8.0) {
            //return 18.31; // in pb for Mh=125 GeV with a scale choice of Mh
            //return 20.18; // in pb for Mh=125 GeV with a scale choice of Mh/2
            return 20.02; // in pb for Mh=125.5 GeV with a scale choice of Mh/2
        } else
            throw std::runtime_error("Bad argument in StandardModel::computeSigmaggH_tt()");
    }

    /**
     * @brief The square of the bottom-quark contribution to the ggH cross section in the Standard Model.
     * @details The values have been obtained from  M. Spira. See also Table 38
     * in ref. @cite Heinemeyer:2013tqa, which are calculated with a scale choice of Mh.
     * @param[in] sqrt_s the center-of-mass energy in TeV
     * @return @f$\sigma_{ggH}^{bb}@f$ in pb
     */
    double computeSigmaggH_bb(const double sqrt_s) const
    {
        if (sqrt_s == 7.0) {
            //return 0.09528; // in pb for Mh=125 GeV with a scale choice of Mh
            //return 0.1079; // in pb for Mh=125 GeV with a scale choice of Mh/2
            return 0.1057; // in pb for Mh=125.5 GeV with a scale choice of Mh/2
        } else if (sqrt_s == 8.0) {
            //return 0.1206; // in pb for Mh=125 GeV with a scale choice of Mh
            //return 0.1357; // in pb for Mh=125 GeV with a scale choice of Mh/2
            return 0.1330; // in pb for Mh=125.5 GeV with a scale choice of Mh/2
        } else
            throw std::runtime_error("Bad argument in StandardModel::computeSigmaggH_bb()");
    }

    /**
     * @brief The top-bottom interference contribution to the ggH cross section in the Standard Model.
     * @details The values have been obtained from  M. Spira. See also Table 38
     * in ref. @cite Heinemeyer:2013tqa, which are calculated with a scale choice of Mh.
     * @param[in] sqrt_s the center-of-mass energy in TeV
     * @return @f$\sigma_{ggH}^{tb}@f$ in pb
     */
    double computeSigmaggH_tb(const double sqrt_s) const
    {
        if (sqrt_s == 7.0) {
            //return -0.8970; // in pb for Mh=125 GeV with a scale choice of Mh
            //return -0.9726; // in pb for Mh=125 GeV with a scale choice of Mh/2
            return -0.9609; // in pb for Mh=125.5 GeV with a scale choice of Mh/2
        } else if (sqrt_s == 8.0) {
            //return -1.125; // in pb for Mh=125 GeV with a scale choice of Mh
            //return -1.206; // in pb for Mh=125 GeV with a scale choice of Mh/2
            return -1.192; // in pb for Mh=125.5 GeV with a scale choice of Mh/2
        } else
            throw std::runtime_error("Bad argument in StandardModel::computeSigmaggH_tb()");
    }

    /**
     * @brief The VBF cross section in the Standard Model.
     * @details See Tables B.67 and B.74 in ref. @cite Heinemeyer:2013tqa
     * @param[in] sqrt_s the center-of-mass energy in TeV
     * @return VBF cross section in pb
     */
    double computeSigmaVBF(const double sqrt_s) const
    {
        if (sqrt_s == 7.0) {
            return 1.222; // in pb for Mh=125.1 GeV
            //return 1.214; // in pb for Mh=125.6 GeV
        } else if (sqrt_s == 8.0) {
            return 1.579; // in pb for Mh=125.1 GeV
            //return 1.572; // in pb for Mh=125.6 GeV
        } else if (sqrt_s == 1.96) {
            return 0.0653; // in pb for Mh=125 GeV
        } else
            throw std::runtime_error("Bad argument in StandardModel::computeSigmaVBF()");
    }

    /**
     * @brief The W fusion contribution @f$\sigma_{WF}@f$ to higgs-production
     * cross section in the Standard Model.
     * @details Currently it returns the value of tab 37 in ref. @cite Heinemeyer:2013tqa
     * @param[in] sqrt_s the center-of-mass energy in TeV
     * @return W fusion contribution @f$\sigma_{WF}@f$ to cross section in pb
     */
    double computeSigmaWF(const double sqrt_s) const
    {
        if (sqrt_s == 7.0) {
            return 0.938; // in pb for Mh=125 GeV
        } else if (sqrt_s == 8.0) {
            return 1.210; // in pb for Mh=125 GeV
        } else if (sqrt_s == 1.96) {
            return computeSigmaVBF(sqrt_s) / computeSigmaVBF(7.) * computeSigmaWF(7.); // in the absence of individual cross sections for TeVatron we rescale the LHC ones
        } else
            throw std::runtime_error("Bad argument in StandardModel::computeSigmaWF()");
    }

    /**
     * @brief The Z fusion contribution @f$\sigma_{ZF}@f$ to higgs-production
     * cross section in the Standard Model.
     * @details Currently it returns the value of tab 37 in ref. @cite Heinemeyer:2013tqa
     * @param[in] sqrt_s the center-of-mass energy in TeV
     * @return W fusion contribution @f$\sigma_{ZF}@f$ to cross section in pb
     */
    double computeSigmaZF(const double sqrt_s) const
    {
        if (sqrt_s == 7.0) {
            return 0.321; // in pb for Mh=125 GeV
        } else if (sqrt_s == 8.0) {
            return 0.417; // in pb for Mh=125 GeV
        } else if (sqrt_s == 1.96) {
            return computeSigmaVBF(sqrt_s) / computeSigmaVBF(7.) * computeSigmaZF(7.); // in the absence of individual cross sections for TeVatron we rescale the LHC ones
        } else
            throw std::runtime_error("Bad argument in StandardModel::computeSigmaZF()");
    }

    /**
     * @brief The Z W interference fusion contribution @f$\sigma_{ZWF}@f$ to
     * higgs-production cross section in the Standard Model.
     * @details Negligible (0.1%) in the Standard model.
     * @param[in] sqrt_s the center-of-mass energy in TeV
     * @return Z W interference fusion contribution @f$\sigma_{ZWF}@f$ to cross section in pb
     */
    double computeSigmaZWF(const double sqrt_s) const
    {
        return 0.;
    }

    /**
     * @brief The WH production cross section in the Standard Model.
     * @details See Tables B.67 and B.74 in ref. @cite Heinemeyer:2013tqa
     * @param[in] sqrt_s the center-of-mass energy in TeV
     * @return WH production cross section in pb
     */
    double computeSigmaWH(const double sqrt_s) const
    {
        if (sqrt_s == 7.0) {
            return 0.577; // in pb for Mh=125.1 GeV
            //return 0.5688; // in pb for Mh=125.6 GeV
        } else if (sqrt_s == 8.0) {
            return 0.7027; // in pb for Mh=125.1 GeV
            //return 0.6931; // in pb for Mh=125.6 GeV
        } else if (sqrt_s == 1.96) {
            return 0.1295; // in pb for Mh=125 GeV
        } else
            throw std::runtime_error("Bad argument in StandardModel::computeSigmaWH()");
    }

    /**
     * @brief The ZH production cross section in the Standard Model.
     * @details See Tables B.67 and B.74 in ref. @cite Heinemeyer:2013tqa
     * @param[in] sqrt_s the center-of-mass energy in TeV
     * @return ZH production cross section in pb
     */
    double computeSigmaZH(const double sqrt_s) const
    {
        if (sqrt_s == 7.0) {
            return 0.3341; // in pb for Mh=125.1 GeV
            //return 0.3299; // in pb for Mh=125.6 GeV
        } else if (sqrt_s == 8.0) {
            return 0.4142; // in pb for Mh=125.1 GeV
            //return 0.4091; // in pb for Mh=125.6 GeV
        } else if (sqrt_s == 1.96) {
            return 0.0785; // in pb for Mh=125 GeV
        } else
            throw std::runtime_error("Bad argument in StandardModel::computeSigmaZH()");
    }

    /**
     * @brief The ttH production cross section in the Standard Model.
     * @details See Tables B.67 and B.74 in ref. @cite Heinemeyer:2013tqa
     * @param[in] sqrt_s the center-of-mass energy in TeV
     * @return ttH production cross section in pb
     */
    double computeSigmattH(const double sqrt_s) const
    {
        if (sqrt_s == 7.0) {
            return 0.08611; // in pb for Mh=125.1 GeV
            //return 0.0851; // in pb for Mh=125.6 GeV
        } else if (sqrt_s == 8.0) {
            return 0.129; // in pb for Mh=125.1 GeV
            //return 0.1274; // in pb for Mh=125.6 GeV
        } else if (sqrt_s == 1.96) {
            return 0.0043; // in pb for Mh=125 GeV
        } else
            throw std::runtime_error("Bad argument in StandardModel::computeSigmattH()");
    }

    /**
     * @brief The Br@f$(H\to gg)@f$ in the Standard Model.
     * @details Currently it returns the value of tables in appendix A (Mh=125.1 GeV) in ref. @cite Heinemeyer:2013tqa
     * @return Br@f$(H\to gg)@f$
     */
    double computeBrHtogg() const
    {
        return 8.57e-2; // Mh=125.1 GeV
    }

    /**
     * @brief The Br@f$(H\to WW)@f$ in the Standard Model.
     * @details Currently it returns the value of tables in appendix A (Mh=125.1 GeV) in ref. @cite Heinemeyer:2013tqa
     * @return Br@f$(H\to WW)@f$ in the Standard Model
     */
    double computeBrHtoWW() const
    {
        //return 2.23e-1; // Mh=125.5 GeV
        return 2.15e-1; // Mh=125.1 GeV
    }

    /**
     * @brief The Br@f$(H\to ZZ)@f$ in the Standard Model.
     * @details Currently it returns the value of tables in appendix A (Mh=125.1 GeV) in ref. @cite Heinemeyer:2013tqa
     * @return Br@f$(H\to ZZ)@f$ in the Standard Model
     */
    double computeBrHtoZZ() const
    {
        return 2.64e-2; // Mh=125.5 GeV
        //return 2.79e-2; // Mh=125.6 GeV
    }

    /**
     * @brief The Br@f$(H\to Z\gamma)@f$ in the Standard Model.
     * @details Currently it returns the value of tables in appendix A (Mh=125.1 GeV) in ref. @cite Heinemeyer:2013tqa
     * @return Br@f$(H\to Z\gamma)@f$ in the Standard Model
     */
    double computeBrHtoZga() const
    {
        return 1.54e-3; // Mh=125.1 GeV
        //return 1.59e-3; // Mh=125.6 GeV
    }

    /**
     * @brief The Br@f$(H\to\gamma\gamma)@f$ in the Standard Model.
     * @details Currently it returns the value of tables in appendix A (Mh=125.1 GeV) in ref. @cite Heinemeyer:2013tqa
     * @return Br@f$(H\to\gamma\gamma)@f$ in the Standard Model
     */
    double computeBrHtogaga() const
    {
        return 2.28e-3; // Mh=125.1 GeV
    }

    /**
     * @brief The Br@f$(H\to \tau\tau)@f$ in the Standard Model.
     * @details Currently it returns the value of tables in appendix A (Mh=125.1 GeV) in ref. @cite Heinemeyer:2013tqa
     * @return Br@f$(H\to \tau\tau)@f$ in the Standard Model
     */
    double computeBrHtotautau() const
    {
        return 6.32e-2; // Mh=125.1 GeV
        //return 6.22e-2; // Mh=125.6 GeV
    }

    /**
     * @brief The Br@f$(H\to cc)@f$ in the Standard Model.
     * @details Currently it returns the value of tables in appendix A (Mh=125.1 GeV) in ref. @cite Heinemeyer:2013tqa
     * @return Br@f$(H\to cc)@f$ in the Standard Model
     */
    double computeBrHtocc() const
    {
        return 2.91e-2; // Mh=125.1 GeV
        //return 2.86e-2; // Mh=125.6 GeV
    }

    /**
     * @brief The Br@f$(H\to bb)@f$ in the Standard Model.
     * @details Currently it returns the value of tables in appendix A (Mh=125.1 GeV) in ref. @cite Heinemeyer:2013tqa
     * @return Br@f$(H\to bb)@f$ in the Standard Model
     */
    double computeBrHtobb() const
    {
        return 5.77e-1; // Mh=125.1 GeV
        //return 5.67e-1; // Mh=125.6 GeV
    }

    /**
     * @brief The Higgs total width in the Standard Model.
     * @details Currently it returns the value of tables in appendix A (Mh=125.1 GeV) in ref. @cite Heinemeyer:2013tqa
     * @return @f$\Gamma_h@f$ in GeV in the Standard Model
     */
    double computeGammaHTotal() const
    {
        return 4.07e-3; // Mh=125.1 GeV
        //return 4.15e-3; // Mh=125.6 GeV
    }

    /**
     * @brief The top loop contribution to @f$H\to gg@f$ in the Standard Model.
     * @details Currently it returns the value of tab 39 in ref. @cite Heinemeyer:2013tqa
     * @return Width of @f$H\to gg@f$ (top loop contribution squared) in keV
     */
    double computeGammaHgg_tt() const
    {
        return 380.8; // in keV for Mh=125 GeV
        //return 389.6; // in keV for Mh=126 GeV
    }

    /**
     * @brief The bottom loop contribution to @f$H\to gg@f$ in the Standard Model.
     * @details Currently it returns the value of tab 39 in ref. @cite Heinemeyer:2013tqa
     * @return Width of @f$H\to gg@f$ (bottom loop contribution squared) in keV
     */
    double computeGammaHgg_bb() const
    {
        return 3.96; // in keV for Mh=125 GeV
        //return 3.95; // in keV for Mh=126 GeV
    }

    /**
     * @brief The top-bottom interference contribution to @f$H\to gg@f$ in the Standard Model.
     * @details Currently it returns the value of tab 39 in ref. @cite Heinemeyer:2013tqa
     * @return Width of @f$H\to gg@f$ (top-bottom interference contribution) in keV
     */
    double computeGammaHgg_tb() const
    {
        return -42.1; // in keV for Mh=125 GeV
        //return -42.7; // in keV for Mh=126 GeV
    }

    /**
     * @brief The top loop contribution to @f$H\to Z\gamma@f$ in the Standard Model.
     * @details Currently it returns the value of tab 41 in ref. @cite Heinemeyer:2013tqa
     * @return Width of @f$H\to Z\gamma@f$ (top loop contribution squared) in eV
     */
    double computeGammaHZga_tt() const
    {
        return 21.74; // in eV for Mh=125 GeV
        //return 23.51; // in eV for Mh=126 GeV
    }

    /**
     * @brief The @f$W@f$ loop contribution to @f$H\to Z\gamma@f$ in the Standard Model.
     * Currently it returns the value of tab 41 in ref. @cite Heinemeyer:2013tqa
     * @return Width of @f$H\to Z\gamma@f$ (W loop contribution squared) in eV
     */
    double computeGammaHZga_WW() const
    {
        return 7005.6; // in eV for Mh=125 GeV
        //return 7648.4; // in eV for Mh=126 GeV
    }

    /**
     * @brief The mixed @f$t-W@f$ loop contribution to @f$H\to Z\gamma@f$ in the Standard Model.
     * @details Currently it returns the value of tab 41 in ref. @cite Heinemeyer:2013tqa
     * @return Width of @f$H\to Z\gamma@f$ (top W loop interference) in eV
     */
    double computeGammaHZga_tW() const
    {
        return -780.4; // in eV for Mh=125 GeV
        //return -848.1; // in eV for Mh=126 GeV
    }

    /**
     * @brief The top loop contribution to @f$H\to\gamma\gamma@f$ in the Standard Model.
     * @details Currently it returns the value of tab 40 in ref. @cite Heinemeyer:2013tqa
     * @return Width of @f$H\to\gamma\gamma@f$ (top loop contribution squared) in eV
     */
    double computeGammaHgaga_tt() const
    {
        return 662.84; // in eV for Mh=125 GeV
        //return 680.39; // in eV for Mh=126 GeV
    }

    /**
     * @brief The @f$W@f$ loop contribution to @f$H\to\gamma\gamma@f$ in the Standard Model.
     * @details Currently it returns the value of tab 40 in ref. @cite Heinemeyer:2013tqa
     * @return Width of @f$H\to\gamma\gamma@f$ (W loop contribution squared) in eV
     */
    double computeGammaHgaga_WW() const
    {
        return 14731.86; // in eV for Mh=125 GeV
        //return 15221.98; // in eV for Mh=126 GeV
    }

    /**
     * @brief The mixed @f$t-W@f$ loop contribution to @f$H\to\gamma\gamma@f$ in the Standard Model.
     * @details Currently it returns the value of tab 40 in ref. @cite Heinemeyer:2013tqa
     * @return Width of @f$H\to\gamma\gamma@f$ (top W loop interference) in eV
     */
    double computeGammaHgaga_tW() const
    {
        return -6249.93; // in eV for Mh=125 GeV
        //return -6436.35; // in eV for Mh=126 GeV
    }

    ////////////////////////////////////////////////////////////////////////
protected:

    /**
     * @brief A method to set the value of a parameter of %StandardModel.
     * @param[in] name name of a model parameter
     * @param[in] value the value to be assigned to the parameter specified by name
     */
    virtual void setParameter(const std::string name, const double& value);

    /**
     * @brief The method to compute the %CKM matrix.
     */
    virtual void computeCKM();

    /**
     * @brief The method to compute the Yukawa matrices.
     * 
     * @attention This function has not been correctly implemented yet. 
     */
    virtual void computeYukawas();

    Particle leptons[6]; ///< An array of Particle objects for the leptons. 
    CKM myCKM; ///< An object of type CKM. 
    gslpp::matrix<gslpp::complex> VCKM; ///< The %CKM matrix.
    gslpp::matrix<gslpp::complex> UPMNS; ///<  The %PMNS matrix.
    gslpp::matrix<gslpp::complex> Yu; ///< The Yukawa matrix of the up-type quarks.
    gslpp::matrix<gslpp::complex> Yd; ///< The Yukawa matrix of the down-type quarks.
    gslpp::matrix<gslpp::complex> Yn; ///< The Yukawa matrix of the neutrinos.
    gslpp::matrix<gslpp::complex> Ye; ///< The Yukawa matrix of the charged leptons.

    // model parameters
    double AlsMz; ///< The strong coupling constant at the Z-boson mass, \f$\alpha_s(M_Z)\f$.
    double Mz; ///< The mass of the \f$Z\f$ boson in GeV. 
    double GF; ///< The Fermi constant @f$G_\mu@f$ in @f${\rm GeV}^{-2}@f$.
    double ale; ///< The fine-structure constant @f$\alpha@f$.
    double dAle5Mz; ///< The five-flavour hadronic contribution to the electromagnetic coupling, @f$\Delta\alpha_{\mathrm{had}}^{(5)}(M_Z^2)@f$.
    double mHl; ///< The Higgs mass @f$m_h@f$ in GeV. 
    double delMw; ///< The theoretical uncertainty in @f$M_W@f$, denoted as @f$\delta\,M_W@f$, in GeV.
    double delSin2th_l; ///< The theoretical uncertainty in @f$\sin^2\theta_{\rm eff}^{\rm lept}@f$, denoted as @f$\delta\sin^2\theta_{\rm eff}^{\rm lept}@f$.
    double delGammaZ; ///< The theoretical uncertainty in @f$\Gamma_Z@f$, denoted as @f$\delta\,\Gamma_Z@f$, in GeV.
    double lambda; ///< The %CKM parameter @f$\lambda@f$ in the Wolfenstein parameterization.
    double A; ///< The %CKM parameter @f$A@f$ in the Wolfenstein parameterization.
    double rhob; ///< The %CKM parameter @f$\bar{\rho}@f$ in the Wolfenstein parameterization.
    double etab; ///< The %CKM parameter @f$\bar{\eta}@f$ in the Wolfenstein parameterization.
    double muw; ///< A matching scale @f$\mu_W@f$ around the weak scale in GeV.

    double EpsK;
    double phiEpsK;
    double DeltaMK;
    double KbarEpsK;
    double Dmk;
    double SM_M12D;


    ////////////////////////////////////////////////////////////////////////
    // For EWPO

    /**
     * @brief An array of internal flags controlling the inclusions of higher-order
     * corrections.
     * @details These flags are prepared for debugging.
     * The flags are initialized in the constructor EWSM().
     */
    bool flag_order[orders_EW_size];

    /**
     * @brief A method to convert a given scheme name in string form into a
     * floating-point number with double precision.
     * @details This method is used in EWSM::checkSMparams() for caching the
     * schemes used in computing @f$M_W@f$, @f$\rho_Z^f@f$ and @f$\kappa_Z^f@f$.
     * @param[in] scheme scheme name that is used in computing @f$M_W@f$,
     * @f$\rho_Z^f@f$ or @f$\kappa_Z^f@f$
     * @return a floating-point number with double precision corresponding to
     * the given scheme name
     */
    double SchemeToDouble(const std::string scheme) const
    {
        if (scheme.compare("NORESUM") == 0)
            return 0.0;
        else if (scheme.compare("OMSI") == 0)
            return 1.0;
        else if (scheme.compare("INTERMEDIATE") == 0)
            return 2.0;
        else if (scheme.compare("OMSII") == 0)
            return 3.0;
        else if (scheme.compare("APPROXIMATEFORMULA") == 0)
            return 4.0;
        else
            throw std::runtime_error("EWSM::SchemeToDouble: bad scheme");
    }

    /**
     * @brief A method to check if a given scheme name in string form is valid.
     * @param[in] scheme scheme name for @f$M_W@f$, @f$\rho_Z^f@f$ or @f$\kappa_Z^f@f$
     * @return a boolean that is true if the scheme name is valid
     */
    bool checkEWPOscheme(const std::string scheme) const
    {
        if (scheme.compare("NORESUM") == 0
                || scheme.compare("OMSI") == 0
                || scheme.compare("INTERMEDIATE") == 0
                || scheme.compare("OMSII") == 0
                || scheme.compare("APPROXIMATEFORMULA") == 0)
            return true;
        else
            return false;
    }

    /**
     * @brief A method to compute the @f$W@f$-boson mass from
     * @f$\Delta\rho@f$ and @f$\Delta r_{\mathrm{rem}}@f$.
     * @details This function computes the @f$W@f$-boson mass without or with
     * resummation of @f$\Delta r@f$,
     * depending on the model flag @ref StandardModelFlags "Mw" of StandardModel:
     *
     * @li NORESUM (recommended):&nbsp;&nbsp; no resummation is considered;
     * @li OMSI:&nbsp;&nbsp; the so-called OMS-I scheme is adopted;
     * @li INTERMEDIATE:&nbsp;&nbsp; an intermediate scheme between OMS-I and OMS-II is adopted;
     * @li OMSII:&nbsp;&nbsp; the so-called OMS-II scheme is adopted;
     * @li APPROXIMATEFORMULA:&nbsp;&nbsp; this is not applicable to the current function.
     *
     * where the OMS-I, INTERMEDIATE and OMS-II schemes are adopted in ZFITTER
     * @cite Bardin:1999yd (see also @cite Degrassi:1996mg, @cite Degrassi:1996ps,
     * @cite Degrassi:1999jd, @cite Bardin:1999ak),
     * and used for making comparisons to the outputs of ZFITTER.
     * The full two-loop %EW contribution is included in the case of "NORESUM",
     * while the large-@f$m_t@f$ expansion for the two-loop contribution is
     * adopted in the other cases.
     *
     * In the case of "NORESUM", the two-loop %EW contribution to @f$\Delta r@f$
     * is calculated via the function EWSMApproximateFormulae::DeltaR_TwoLoopEW_rem(),
     * given in the complex-pole/fixed-width scheme. The @f$W@f$-boson mass in
     * the complex-pole/fixed-width scheme, obtained from @f$\Delta r@f$, is
     * converted into the one in the experimental/running-width scheme with the
     * function MwFromMwbar().
     * @param[in] Mw_i the @f$W@f$-boson mass
     * @param[in] DeltaRho Array of @f$\Delta\rho@f$
     * @param[in] DeltaR_rem Array of @f$\Delta r_{\mathrm{rem}}@f$
     * @return @f$M_W@f$
     */
    double resumMw(const double Mw_i, const double DeltaRho[orders_EW_size],
            const double DeltaR_rem[orders_EW_size]) const;

    /**
     * @brief A method to compute the real part of the effective coupling
     * @f$\rho_Z^f@f$ from @f$\Delta\rho@f$, @f$\delta\rho_{\rm rem}^{f}@f$
     * and @f$\Delta r_{\mathrm{rem}}@f$.
     * @details This function computes @f$\rho_Z^f@f$ without or with
     * resummation of @f$\Delta\rho@f$, depending on
     * the model flag @ref StandardModelFlags "RhoZ" of StandardModel:
     *
     * @li NORESUM (recommended):&nbsp;&nbsp; no resummation is considered;
     * @li OMSI:&nbsp;&nbsp; the so-called OMS-I scheme is adopted;
     * @li INTERMEDIATE:&nbsp;&nbsp; an intermediate scheme between OMS-I and OMS-II is adopted;
     * @li OMSII:&nbsp;&nbsp; the so-called OMS-II scheme is adopted;
     * @li APPROXIMATEFORMULA:&nbsp;&nbsp; this is not applicable to the current function.
     *
     * where the OMS-I, INTERMEDIATE and OMS-II schemes are adopted in ZFITTER
     * @cite Bardin:1999yd (see also @cite Degrassi:1996mg, @cite Degrassi:1996ps,
     * @cite Degrassi:1999jd, @cite Bardin:1999ak),
     * and used for making comparisons to the outputs of ZFITTER.
     * In all the cases, the two-loop %EW corrections are calculated in the
     * large-@f$m_t@f$ expansion.
     * @param[in] DeltaRho Array of @f$\Delta\rho@f$
     * @param[in] deltaRho_rem Array of @f$\delta\rho_{\rm rem}^{f}@f$
     * @param[in] DeltaRbar_rem Array of @f$\Delta \bar{r}_{\rm rem}@f$
     * @param[in] bool_Zbb true for @f$Zb\bar{b}@f$
     * @return @f$\mathrm{Re}(\rho_Z^f)@f$
     */
    double resumRhoZ(const double DeltaRho[orders_EW_size],
            const double deltaRho_rem[orders_EW_size],
            const double DeltaRbar_rem, const bool bool_Zbb) const;

    /**
     * @brief A method to compute the real part of the effetvive coupling
     * @f$\kappa_Z^f@f$ from @f$\Delta\rho@f$, @f$\delta\rho_{\rm rem}^{f}@f$
     * and @f$\Delta r_{\mathrm{rem}}@f$.
     * @details This function computes @f$\kappa_Z^f@f$ without or with
     * resummation of @f$\Delta\rho@f$, depending on
     * the model flag @ref StandardModelFlags "KappaZ" of StandardModel:
     *
     * @li NORESUM (recommended):&nbsp;&nbsp; no resummation is considered;
     * @li OMSI:&nbsp;&nbsp; the so-called OMS-I scheme is adopted;
     * @li INTERMEDIATE:&nbsp;&nbsp; an intermediate scheme between OMS-I and OMS-II is adopted;
     * @li OMSII:&nbsp;&nbsp; the so-called OMS-II scheme is adopted;
     * @li APPROXIMATEFORMULA:&nbsp;&nbsp; this is not applicable to the current function.
     *
     * where the OMS-I, INTERMEDIATE and OMS-II schemes are adopted in ZFITTER
     * @cite Bardin:1999yd (see also @cite Degrassi:1996mg, @cite Degrassi:1996ps,
     * @cite Degrassi:1999jd, @cite Bardin:1999ak),
     * and used for making comparisons to the outputs of ZFITTER.
     * In all the cases, the two-loop %EW corrections are calculated in the
     * large-@f$m_t@f$ expansion.
     * @param[in] DeltaRho Array of @f$\Delta\rho@f$
     * @param[in] deltaKappa_rem Array of @f$\delta\kappa_{\rm rem}^{f}@f$
     * @param[in] DeltaRbar_rem Array of @f$\Delta \bar{r}_{\rm rem}@f$
     * @param[in] bool_Zbb true for @f$Zb\bar{b}@f$
     * @return @f$\mathrm{Re}(\kappa_Z^f)@f$
     */
    double resumKappaZ(const double DeltaRho[orders_EW_size],
            const double deltaKappa_rem[orders_EW_size],
            const double DeltaRbar_rem, const bool bool_Zbb) const;

    /**
     * @brief Top-mass corrections to the @f$Zb\bar{b}@f$ vertex, denoted by
     * @f$\tau_b@f$.
     * @details The large top-quark mass gives important corrections to the
     * %EW observables through the gauge-boson self-energies, i.e.,
     * @f$\Delta\rho@f$, and through the @f$Zb\bar{b}@f$ vertex. The latter
     * contribution is parameterised by the quantity @f$\tau_b@f$:
     * @f[
     * \tau_{b} =
     * -2\, X_t^{G_\mu}
     * \left[ 1 - \frac{\pi}{3}\alpha_s(M^2_t)
     *   + X_t^{G_\mu} \tau^{(2)}
     *     \left( \frac{M_t^2}{m_h^2} \right)
     * \right],
     * @f]
     * where the @f$O(G_\mu\alpha_s m_t^2)@f$ term was calculated in
     * @cite Fleischer:1992fq, @cite Buchalla:1992zm, @cite Degrassi:1993ij,
     * @cite Chetyrkin:1993jp, and the @f$O(G_\mu^2 m_t^4)@f$ term can be found
     * in @cite Barbieri:1992nz, @cite Barbieri:1992dq, @cite Fleischer:1993ub,
     * @cite Fleischer:1994cb.
     * @return @f$\tau_b@f$
     */
    double taub() const;

    /**
     * @brief The non-factorizable EW-%QCD corrections to the partial widths
     * for @f$Z\to q\bar{q}@f$, denoted as @f$\Delta_{\mathrm{EW/QCD}}@f$.
     * @details
     * See @cite Czarnecki:1996ei and @cite Harlander:1997zb.
     * @param[in] q name of a quark (see QCD::quark)
     * @return @f$\Delta_{\mathrm{EW/QCD}}@f$ in GeV
     */
    double Delta_EWQCD(const QCD::quark q) const;

    /**
     * @brief The radiator factor associated with the final-state QED and %QCD
     * corrections to the the vector-current interactions, @f$R_V^q(M_Z^2)@f$.
     * @details
     * See @cite Chetyrkin:1994js, @cite Bardin:1999ak, @cite Bardin:1999yd
     * and references therein.
     * @param[in] q name of a quark (see QCD::quark)
     * @return @f$R_V^q(M_Z^2)@f$
     */
    double RVq(const QCD::quark q) const;

    /**
     * @brief The radiator factor associated with the final-state QED and %QCD
     * corrections to the the axial-vector-current interactions, @f$R_A^q(M_Z^2)@f$.
     * @details
     * See @cite Chetyrkin:1994js, @cite Bardin:1999ak, @cite Bardin:1999yd,
     * @cite Baikov:2012er and references therein.
     * @param[in] q name of a quark (see QCD::quark)
     * @return @f$R_A^q(M_Z^2)@f$
     */
    double RAq(const QCD::quark q) const;

    /**
     * @brief The singlet vector corrections to the hadronic @f$Z@f$-boson width,
     * denoted as @f$R_V^h@f$.
     * @details In addition to the final-state corrections represented by
     * the radiator factors @f$R_V^q(M_Z^2)@f$ and @f$R_A^q(M_Z^2)@f$,
     * there exist singlet vector corrections to the total hadronic width
     * @cite Chetyrkin:1994js, @cite Baikov:2012er, which is much smaller than
     * the other corrections.
     *
     * The assignment of the singlet vector corrections to the partial widths
     * is ambiguous @cite Bardin:1997xq. See Gamma_had() for our prescription.
     * @return @f$R_V^h@f$
     */
    double RVh() const;

    bool requireCKM; ///< An internal flag to control whether the %CKM matrix has to be recomputed.
    bool requireYe; ///< An internal flag to control whether the charged-lepton Yukawa matrix has to be recomputed.
    bool requireYn; ///<  An internal flag to control whether the neutrino Yukawa matrix has to be recomputed.


    ////////////////////////////////////////////////////////////////////////    
private:

    EWSMcache* myEWSMcache; ///< A pointer to an object of type EWSMcache.
    EWSMOneLoopEW* myOneLoopEW; ///< A pointer to an object of type EWSMOneLoopEW.
    EWSMTwoLoopQCD* myTwoLoopQCD; ///< A pointer to an object of type EWSMTwoLoopQCD.
    EWSMThreeLoopQCD* myThreeLoopQCD; ///< A pointer to an object of type EWSMThreeLoopQCD.
    EWSMTwoLoopEW* myTwoLoopEW; ///< A pointer to an object of type EWSMTwoLoopEW.
    EWSMThreeLoopEW2QCD* myThreeLoopEW2QCD; ///< A pointer to an object of type EWSMThreeLoopEW2QCD.
    EWSMThreeLoopEW* myThreeLoopEW; ///< A pointer to an object of type EWSMThreeLoopEW.
    EWSMApproximateFormulae* myApproximateFormulae; ///< A pointer to an object of type EWSMApproximateFormulae.
    Flavour* myFlavour; ///< A pointer to an object of the type Flavour.
    /** BEGIN: REMOVE FROM THE PACKAGE **/
    EWSMTwoFermionsLEP2* myTwoFermionsLEP2; ///< A pointer to an object of type EWSMTwoFermionsLEP2.
    LeptonFlavour* myLeptonFlavour; ///< A pointer to an object of the type LeptonFlavour.
    /** END: REMOVE FROM THE PACKAGE **/

    StandardModelMatching* myStandardModelMatching; ///< A pointer to an object of type StandardModelMatching.

    bool FlagWithoutNonUniversalVC; ///< A boolean for the model flag %WithoutNonUniversalVC.
    bool FlagNoApproximateGammaZ; ///< A boolean for the model flag %NoApproximateGammaZ.
    std::string FlagMw; ///< A string for the model flag %Mw.
    std::string FlagRhoZ; ///< A string for the model flag %RhoZ.
    std::string FlagKappaZ; ///< A string for the model flag %KappaZ.



    ////////////////////////////////////////////////////////////////////////     
    // Caches for EWPO

    bool FlagCacheInStandardModel; ///< A flag for caching (true by default).
    mutable double SMparamsForEWPO_cache[NumSMParamsForEWPO];
    mutable double DeltaAlphaLepton_cache; ///< A cache of the value of @f$\Delta\alpha_{\mathrm{lept}}(M_Z^2)@f$.
    mutable double DeltaAlpha_cache; ///< A cache of the value of @f$\Delta\alpha(M_Z^2)@f$.
    mutable double Mw_cache; ///< A cache of the value of @f$M_W@f$.
    mutable double GammaW_cache; ///< A cache of the value of @f$\Gamma_W@f$.
    mutable gslpp::complex rhoZ_f_cache[12]; ///< A cache of the value of @f$\rho_Z^l@f$.
    mutable gslpp::complex kappaZ_f_cache[12]; ///< A cache of the value of @f$\kappa_Z^l@f$.
    mutable bool useDeltaAlphaLepton_cache;
    mutable bool useDeltaAlpha_cache;
    mutable bool useMw_cache;
    mutable bool useGammaW_cache;
    mutable bool useRhoZ_f_cache[12];
    mutable bool useKappaZ_f_cache[12];
    
    int iterationNo;

};

#endif	/* STANDARDMODEL_H */
