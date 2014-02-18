/*
 * Copyright (C) 2012-2014 SusyFit Collaboration
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

using namespace gslpp;
class EWSM; // forward reference to EWSM class

/**
 * @class StandardModel
 * @ingroup StandardModel
 * @brief A model class for the Standard %Model.
 * @author SusyFit Collaboration
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
 * allocates memory to the pointers #myEWSM and #myStandardModelMatching. 
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
 * The flags of %StandardModel, which can be specified in the model configuration
 * file, are summarized below, where the values of the boolean flags (TRUE or
 * FALSE) are case insensitive, while those of the other flags are case sensitive.
 * The default values of the flags are indicated in bold: 
 * 
 * <table class="model">
 * <tr>
 *   <th>Label</th>
 *   <th>Value</th>
 *   <th>Description</th>
 * </tr>
 * <tr>
 *   <td class="mod_name">%CacheInEWSM</td>
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
 * The first two flags CacheInEWSM and CacheInEWSMcache for the cashing methods
 * in EWSM and EWSMcache classes are relevant to the computations of the electroweak 
 * precision observables. Those caches are effective when the @f$W@f$-boson mass,
 * the decay widths of the @f$Z@f$ boson and the @f$Zf\bar{f}@f$ effective couplings
 * @f$\kappa_Z^f@f$ are calculated without using their two-loop approximate formulae.
 *
 *
 */
class StandardModel: public QCD {
public:

    /**
     * @brief An enum type for leptons.
     */
    enum lepton
    {
        NEUTRINO_1, /**< The 1st-generation neutrino */
        ELECTRON, /**< Electron */
        NEUTRINO_2, /**< The 2nd-generation neutrino */
        MU, /**< Muon */
        NEUTRINO_3, /**< The 3rd-generation neutrino */
        TAU /**< Tau */
    };

    static const int NSMvars = 24;///< The number of the model parameters in %StandardModel.

    /**
     * @brief  A string array containing the labels of the model parameters in %StandardModel.
     */
    static const std::string SMvars[NSMvars];

    /**
     * @brief The default constructor.
     */
    StandardModel();

    /**
     * @brief A method to fetch the name of %StandardModel.
     * @return the name of the model as a string 
     */
    virtual std::string ModelName() const
    {
        return "StandardModel";
    }


    ///////////////////////////////////////////////////////////////////////////
    // Initialization

    /**
     * @brief A method to initialize the model.
     * @details This method, called via InputParser::ReadParameters(), allocates
     * memory to the pointers #myEWSM and #myStandardModelMatching, which are used
     * for %EW precision and flavour observables, respectively. 
     * @return a boolean that is true if model initialization is successful
     */
    virtual bool InitializeModel();

    /**
     * @brief A get method to access the member pointer of type EWSM.
     * @return the pointer #myEWSM
     */
    EWSM* getEWSM() const
    {
        return myEWSM;
    }

    /**
     * @brief A get method to access the member pointer of type StandardModelMatching.
     * @return the pointer #myStandardModelMatching
     */
    virtual StandardModelMatching* getMyMatching() const
    {
        return myStandardModelMatching;
    }


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

    /**
     * @brief A method to check if all the mandatory parameters for %StandardModel
     * have been provided in the model configuration file.
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
    matrix<complex> getVCKM() const
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
    matrix<complex> getUPMNS() const
    {
        return UPMNS;
    }

    /**
     * @brief A get method to retrieve the Yukawa matrix of the up-type quarks, 
     * @f$Y_u@f$.
     * @return @f$Y_u@f$
     */
    matrix<complex> getYu() const
    {
        return Yu;
    }


    /**
     * @brief A get method to retrieve the Yukawa matrix of the down-type quarks,
     * @f$Y_d@f$.
     * @return @f$Y_d@f$
     */
    matrix<complex> getYd() const
    {
        return Yd;
    }

    /**
     * @brief A get method to retrieve the Yukawa matrix of the neutrinos,
     * @f$Y_\nu@f$.
     * @return @f$Y_\nu@f$
     */
    matrix<complex> getYn() const
    {
        return Yn;
    }

    /**
     * @brief A get method to retrieve the Yukawa matrix of the charged leptons,
     * @f$Y_e@f$.
     * @return @f$Y_e@f$
     */
    matrix<complex> getYe() const
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

    double getEpsK() const {
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


    ///////////////////////////////////////////////////////////////////////////

    /**
     * @brief The Higgs vacuum expectation value.
     * @f[
     *    v = \left(\frac{1}{\sqrt{2} G_\mu}\right)^{1/2},
     * @f]
     * where @f$G_\mu@f$ is the Fermi constant, measured through muon decays. 
     * @return @f$v@f$ in GeV
     */
    virtual double v() const;

    /**
     * @brief The tree-level mass of the @f$W@f$ boson, @f$M_W^{\mathrm{tree}}@f$.
     * @return @f$M_W^{\mathrm{tree}}@f$ in GeV.
     */
    virtual double Mw_tree() const;

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
    double ale_OS(const double mu, orders order=FULLNLO) const;


    ////////////////////////////////////////////////////////////////////////

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

    /**
     * @brief The @f$W@f$-boson mass, @f$M_W@f$.
     * @return @f$M_W@f$ in GeV
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
    virtual double sW2() const;

    /**
     * @brief The total width of the @f$W@f$ boson, @f$\Gamma_W@f$.
     * @return @f$\Gamma_W@f$ in GeV
     */
    virtual double GammaW() const;


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
    complex computelamt() const;

    /**
     * @brief The product of the %CKM elements @f$V_{cd} V_{cs}^*@f$.
     * @return @f$V_{cd} V_{cs}^*@f$
     */
    complex computelamc() const;

    /**
     * @brief The product of the %CKM elements @f$V_{ud} V_{us}^*@f$.
     * @return @f$V_{ud} V_{us}^*@f$
     */
    complex computelamu() const;

    /**
     * @brief The product of the %CKM elements @f$V_{td} V_{tb}^*@f$.
     * @return @f$V_{td} V_{tb}^*@f$
     */
    complex computelamt_d() const;

    /**
     * @brief The product of the %CKM elements @f$V_{cd} V_{cb}^*@f$.
     * @return @f$V_{cd} V_{cb}^*@f$
     */
    complex computelamc_d() const;

    /**
     * @brief The product of the %CKM elements @f$V_{ud} V_{ub}^*@f$.
     * @return @f$V_{ud} V_{ub}^*@f$
     */
    complex computelamu_d() const;

    /**
     * @brief The product of the %CKM elements @f$V_{ts} V_{tb}^*@f$.
     * @return @f$V_{ts} V_{tb}^*@f$
     */
    complex computelamt_s() const;

    /**
     * @brief The product of the %CKM elements @f$V_{cs} V_{cb}^*@f$.
     * @return @f$V_{cs} V_{cb}^*@f$
     */
    complex computelamc_s() const;

    /**
     * @brief The product of the %CKM elements @f$V_{us} V_{ub}^*@f$.
     * @return @f$V_{us} V_{ub}^*@f$
     */
    complex computelamu_s() const;

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

    EWSM* myEWSM;///< A pointer to an object of type EWSM.

    Particle leptons[6];///< An array of Particle objects for the leptons. 
    CKM myCKM;///< An object of type CKM. 
    matrix<complex> VCKM;///< The %CKM matrix.
    matrix<complex> UPMNS;///<  The %PMNS matrix.
    matrix<complex> Yu;///< The Yukawa matrix of the up-type quarks.
    matrix<complex> Yd;///< The Yukawa matrix of the down-type quarks.
    matrix<complex> Yn;///< The Yukawa matrix of the neutrinos.
    matrix<complex> Ye;///< The Yukawa matrix of the charged leptons.

    // model parameters
    double GF;///< The Fermi constant @f$G_\mu@f$ in @f${\rm GeV}^{-2}@f$.
    double ale;///< The fine-structure constant @f$\alpha@f$.
    double dAle5Mz;///< The five-flavour hadronic contribution to the electromagnetic coupling, @f$\Delta\alpha_{\mathrm{had}}^{(5)}(M_Z^2)@f$.
    double mHl;///< The Higgs mass @f$m_h@f$ in GeV. 
    double delMw;///< The theoretical uncertainty in @f$M_W@f$, denoted as @f$\delta\,M_W@f$, in GeV.
    double delSin2th_l;///< The theoretical uncertainty in @f$\sin^2\theta_{\rm eff}^{\rm lept}@f$, denoted as @f$\delta\sin^2\theta_{\rm eff}^{\rm lept}@f$.
    double delGammaZ;///< The theoretical uncertainty in @f$\Gamma_Z@f$, denoted as @f$\delta\,\Gamma_Z@f$, in GeV.
    double lambda;///< The %CKM parameter @f$\lambda@f$ in the Wolfenstein parameterization.
    double A;///< The %CKM parameter @f$A@f$ in the Wolfenstein parameterization.
    double rhob;///< The %CKM parameter @f$\bar{\rho}@f$ in the Wolfenstein parameterization.
    double etab;///< The %CKM parameter @f$\bar{\eta}@f$ in the Wolfenstein parameterization.
    double muw;///< A matching scale @f$\mu_W@f$ around the weak scale in GeV.

    double EpsK;
    double phiEpsK;
    double DeltaMK;
    double KbarEpsK;
    double Dmk;
    double SM_M12D;
    
    
    ////////////////////////////////////////////////////////////////////////    
private:

    StandardModelMatching* myStandardModelMatching;///< A pointer to an object of type StandardModelMatching.

    bool FlagWithoutNonUniversalVC;///< A boolean for the model flag %WithoutNonUniversalVC.
    bool FlagNoApproximateGammaZ;///< A boolean for the model flag %NoApproximateGammaZ.
    std::string FlagMw;///< A string for the model flag %Mw.
    std::string FlagRhoZ;///< A string for the model flag %RhoZ.
    std::string FlagKappaZ;///< A string for the model flag %KappaZ. 

    bool requireCKM;///< An internal flag to control whether the %CKM matrix has to be recomputed.
    bool requireYe;///< An internal flag to control whether the charged-lepton Yukawa matrix has to be recomputed.
    bool requireYn;///<  An internal flag to control whether the neutrino Yukawa matrix has to be recomputed.

};

#endif	/* STANDARDMODEL_H */
