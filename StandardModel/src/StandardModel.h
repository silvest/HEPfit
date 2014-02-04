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
 * @brief A class for the Standard %Model.
 * @author SusyFit Collaboration
 * @copyright GNU General Public License
 * @details This class is ....
 *
 *
 * The model parameters of StandardModel are summarized below: 
 * @anchor StandardModelParameters
 * <table class="model">
 * <tr>
 *   <th>Parameter</th>
 *   <th>LaTeX symbol</th>
 *   <th>Description</th>
 * </tr>
 * <tr>
 *   <td class="mod_name">#GF</td>
 *   <td class="mod_symb">@f$G_\mu@f$</td>
 *   <td class="mod_desc">The Fermi constant in @f${\rm GeV}^{-2}@f$, measured through muon decays.</td>
 * </tr>
 * <tr>
 *   <td class="mod_name">#ale</td>
 *   <td class="mod_symb">@f$\alpha@f$</td>
 *   <td class="mod_desc">The fine-structure constant.</td>
 * </tr>
 * <tr>
 *   <td class="mod_name">#dAle5Mz</td>
 *   <td class="mod_symb">@f$\Delta\alpha_{\mathrm{had}}^{(5)}(M_Z^2)@f$</td>
 *   <td class="mod_desc">The five-flavour hadronic contribution to the electromagnetic coupling.</td>
 * </tr>
 * <tr>
 *   <td class="mod_name">#mHl</td>
 *   <td class="mod_symb">@f$m_h@f$</td>
 *   <td class="mod_desc">The Higgs mass in GeV.</td>
 * </tr>
 * <tr>
 *   <td class="mod_name">#delMw</td>
 *   <td class="mod_symb">@f$\delta\,M_W@f$</td>
 *   <td class="mod_desc">The theoretical uncertainty in @f$M_W@f$ in GeV.</td>
 * </tr>
 * <tr>
 *   <td class="mod_name">#delSin2th_l</td>
 *   <td class="mod_symb">@f$\delta\sin^2\theta_{\rm eff}^{\rm lept}@f$</td>
 *   <td class="mod_desc">The theoretical uncertainty in @f$\sin^2\theta_{\rm eff}^{\rm lept}@f$.</td>
 * </tr>
 * <tr>
 *   <td class="mod_name">#delGammaZ</td>
 *   <td class="mod_symb">@f$\delta\,\Gamma_Z@f$</td>
 *   <td class="mod_desc">The theoretical uncertainty in @f$\Gamma_Z@f$ in GeV.</td>
 * </tr>
 * <tr>
 *   <td class="mod_name">mneutrino_1</td>
 *   <td class="mod_symb">@f$m_{\nu_1}@f$</td>
 *   <td class="mod_desc">The mass of the first-generation neutrino in GeV.</td>
 * </tr>
 * <tr>
 *   <td class="mod_name">mneutrino_2</td>
 *   <td class="mod_symb">@f$m_{\nu_2}@f$</td>
 *   <td class="mod_desc">The mass of the second-generation neutrino in GeV.</td>
 * </tr>
 * <tr>
 *   <td class="mod_name">mneutrino_3</td>
 *   <td class="mod_symb">@f$m_{\nu_3}@f$</td>
 *   <td class="mod_desc">The mass of the third-generation neutrino in GeV.</td>
 * </tr>
 * <tr>
 *   <td class="mod_name">melectron</td>
 *   <td class="mod_symb">@f$m_e@f$</td>
 *   <td class="mod_desc">The electron mass in GeV.</td>
 * </tr>
 * <tr>
 *   <td class="mod_name">mmu</td>
 *   <td class="mod_symb">@f$m_\mu@f$</td>
 *   <td class="mod_desc">The muon mass in GeV.</td>
 * </tr>
 * <tr>
 *   <td class="mod_name">mtau</td>
 *   <td class="mod_symb">@f$m_\tau@f$</td>
 *   <td class="mod_desc">The tau mass in GeV.</td>
 * </tr>
 * <tr>
 *   <td class="mod_name">#lambda</td>
 *   <td class="mod_symb">@f$\lambda@f$</td>
 *   <td class="mod_desc">The %CKM parameter @f$\lambda@f$ in the Wolfenstein parameterization.</td>
 * </tr>
 * <tr>
 *   <td class="mod_name">#A</td>
 *   <td class="mod_symb">@f$A@f$</td>
 *   <td class="mod_desc">The %CKM parameter @f$A@f$ in the Wolfenstein parameterization.</td>
 * </tr>
 * <tr>
 *   <td class="mod_name">#rhob</td>
 *   <td class="mod_symb">@f$\bar{\rho}@f$</td>
 *   <td class="mod_desc">The %CKM parameter @f$\bar{\rho}@f$ in the Wolfenstein parameterization.</td>
 * </tr>
 * <tr>
 *   <td class="mod_name">#etab</td>
 *   <td class="mod_symb">@f$\bar{\eta}@f$</td>
 *   <td class="mod_desc">The %CKM parameter @f$\bar{\eta}@f$ in the Wolfenstein parameterization.</td>
 * </tr>
 * <tr>
 *   <td class="mod_name">#muw</td>
 *   <td class="mod_symb">@f$\mu_W@f$</td>
 *   <td class="mod_desc">A matching scale around the @f$W@f$-boson mass scale.</td>
 * </tr>
 * </table>
 *
 * 
 * Other model parameters (to be removed?):
 * \li \b EpsK:&nbsp; the experimental value of @f$\varepsilon_{K}@f$,
 * \li \b phiEpsK:&nbsp; the experimental value of @f$\Delta M_{K}/(\Delta\Gamma_{K}/2)@f$,
 * \li \b DeltaMK:&nbsp; the experimental value of @f$\Delta M_{K}@f$ in GeV,
 * \li \b KbarEpsK:&nbsp;
 * \li \b Dmk:&nbsp; the SM contribution to @f$\Delta m_{K}@f$ in GeV,
 * \li \b SM_M12D:&nbsp; the SM amplitude of the @f$D^{0}-\bar{D}^{0}@f$ mixing,
 *
 *
 *
 * The flags of StandardModel: 
 * @anchor StandardModelFlags
 * <table class="model">
 * <tr>
 *   <th>Flag</th>
 *   <th>Value</th>
 *   <th>Description</th>
 * </tr>
 * <tr>
 *   <td class="mod_name">CacheInEWSM</td>
 *   <td class="mod_valu">true / false</td>
 *   <td class="mod_desc">This flag controls the use of the cashing method
 *   implemented in EWSM class. </td>
 * </tr>
 * <tr>
 *   <td class="mod_name">CacheInEWSMcache</td>
 *   <td class="mod_valu">true / false</td>
 *   <td class="mod_desc"></td>
 * </tr>
 * <tr>
 *   <td class="mod_name">WithoutNonUniversalVC</td>
 *   <td class="mod_valu">true / false</td>
 *   <td class="mod_desc"></td>
 * </tr>
 * <tr>
 *   <td class="mod_name">NoApproximateGammaZ</td>
 *   <td class="mod_valu">true / false</td>
 *   <td class="mod_desc"></td>
 * </tr>
 * <tr>
 *   <td class="mod_name">NoApproximateSigmaH</td>
 *   <td class="mod_valu">true / false</td>
 *   <td class="mod_desc"></td>
 * </tr>
 * <tr>
 *   <td class="mod_name">NoApproximateRl</td>
 *   <td class="mod_valu">true / false</td>
 *   <td class="mod_desc"></td>
 * </tr>
 * <tr>
 *   <td class="mod_name">NoApproximateRc</td>
 *   <td class="mod_valu">true / false</td>
 *   <td class="mod_desc"></td>
 * </tr>
 * <tr>
 *   <td class="mod_name">NoApproximateRb</td>
 *   <td class="mod_valu">true / false</td>
 *   <td class="mod_desc"></td>
 * </tr>
 * <tr>
 *   <td class="mod_name">Mw</td>
 *   <td class="mod_valu">NORESUM / OMSI / INTERMEDIATE / OMSII / APPROXIMATEFORMULA</td>
 *   <td class="mod_desc"></td>
 * </tr>
 * <tr>
 *   <td class="mod_name">RhoZ</td>
 *   <td class="mod_valu">NORESUM / OMSI / INTERMEDIATE / OMSII</td>
 *   <td class="mod_desc"></td>
 * </tr>
 * <tr>
 *   <td class="mod_name">KappaZ</td>
 *   <td class="mod_valu">NORESUM / OMSI / INTERMEDIATE / OMSII / APPROXIMATEFORMULA</td>
 *   <td class="mod_desc"></td>
 * </tr>
 * </table>
 *
 *
 */
class StandardModel: public QCD {
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

    static const int NSMvars = 24;///< The number of model parameters in the current model.

    /**
     * @brief An array containing the labels under which all StandardModel
     * parameters are stored in a vector of ModelParameter via
     * InputParser::ReadParameters().
     */
    static const std::string SMvars[NSMvars];

    /**
     * @brief The default constructor.
     */
    StandardModel();

    /**
     * @brief A method to fetch the name of the current model.
     * @return the name of the model as a string 
     */
    virtual std::string ModelName() const
    {
        return "StandardModel";
    }


    ///////////////////////////////////////////////////////////////////////////
    // Initialization

    /**
     * @brief A method to initialize the current model.
     * @details This method, called via InputParser::ReadParameters(), allocates
     * memory to the pointers #myEWSM and #myStandardModelMatching, which are used
     * for %EW precision and flavour observables, respectively. 
     * @return a boolean that is true if model initialization is successful
     */
    virtual bool InitializeModel();

    /**
     * @brief A get method to access the pointer #myEWSM.
     * @return #myEWSM
     */
    EWSM* getEWSM() const
    {
        return myEWSM;
    }

    /**
     * @brief A get method to access the pointer #myStandardModelMatching.
     * @return #myStandardModelMatching
     */
    virtual StandardModelMatching* getMyMatching() const
    {
        return myStandardModelMatching;
    }


    ///////////////////////////////////////////////////////////////////////////
    // Model parameters

    /**
     * @brief A method to initialize the model parameters.
     * @param[in] Dpars a map of parameters that are being updated in the Monte Carlo run
     * @return a boolean that is true if the execution is successful
     */
    virtual bool Init(const std::map<std::string, double>& DPars);

    /**
     * @brief The pre-update method for the current model. 
     * @details This method initializes the internal flags #requireCKM, #requireYe
     * and #requireYn, and calls QCD::PreUpdate(), before updating the model
     * parameters with the method Update().
     * @return a boolean that is true if the execution is successful
     */
    virtual bool PreUpdate();

    /**
     * @brief The update method for the current model.
     * @details This method updates all the model parameters with giving Dpars.
     * @param[in] Dpar a map of parameters that are being updated in the Monte Carlo run
     * (including parameters that are varied and those that are held constant)
     * @return a boolean that is true if the execution is successful
     */
    virtual bool Update(const std::map<std::string, double>& DPars);

    /**
     * @brief The post-update method for the current model.
     * @details This method runs all the procedures that are need to be executed
     * after the model is successfully updated. This includes updating
     * any other variable that needs to be updated at this time due to the update
     * of the model parameters
     * @return a boolean that is true if the execution is successful
     */
    virtual bool PostUpdate();

    /**
     * @brief A method to check if all the mandatory parameters for the model 
     * have been provided in the model configuration file.
     * @param[in] Dpar a map of parameters that are being updated in the Monte Carlo run
     * (including parameters that are varied and those that are held constant)
     * @return a boolean that is true if the execution is successful
     */
    virtual bool CheckParameters(const std::map<std::string, double>& DPars);


    ///////////////////////////////////////////////////////////////////////////
    // Flags

    /**
     * @brief A method to set a flag of the current model.
     * @param[in] name name of a flag
     * @param[in] value the boolean to be assigned to the flag specified by name
     * @return a boolean that is true if the execution is successful
     */
    virtual bool setFlag(const std::string name, const bool value);

    /**
     * @brief A method to set a flag of the current model.
     * @param[in] name name of a flag
     * @param[in] value the string to be assigned to the flag specified by name
     * @return a boolean that is true if the execution is successful
     */
    virtual bool setFlagStr(const std::string name, const std::string value);

    /**
     * @brief A method to check the sanity of the set of flags.
     * @return true if the set of flags is sane
     */
    virtual bool CheckFlags() const;

    /**
     * @brief A method to control
     * @attention The flag FlagWithoutNonUniversalVC is applicable for the model
     * NPEpsilons. 
     * @return a boolean that is true if flavour non-universal vertex corrections
     * are NOT added to the epsilon parameters describing new physics contribution.
     */
    bool IsFlagWithoutNonUniversalVC() const
    {
        return FlagWithoutNonUniversalVC;
    }

    /**
     * @brief A method to access the boolean value of the flag #FlagNoApproximateGammaZ.
     * @return a boolean that is true if
     */
    bool IsFlagNoApproximateGammaZ() const
    {
        return FlagNoApproximateGammaZ;
    }

    /**
     * @brief A method to access the boolean value of the flag #FlagNoApproximateSigmaH.
     * @return #FlagNoApproximateSigmaH
     */
    bool IsFlagNoApproximateSigmaH() const
    {
        return FlagNoApproximateSigmaH;
    }

    /**
     * @brief A method to access the boolean value of the flag #FlagNoApproximateRl.
     * @return #FlagNoApproximateRl
     */
    bool IsFlagNoApproximateRl() const
    {
        return FlagNoApproximateRl;
    }

    /**
     * @brief A method to access the boolean value of the flag #FlagNoApproximateRc.
     * @return #FlagNoApproximateRc
     */
    bool IsFlagNoApproximateRc() const
    {
        return FlagNoApproximateRc;
    }

    /**
     * @brief A method to access the boolean value of the flag #FlagNoApproximateRb.
     * @return #FlagNoApproximateRb
     */
    bool IsFlagNoApproximateRb() const
    {
        return FlagNoApproximateRb;
    }

    std::string getFlagMw() const
    {
        return FlagMw;
    }

    std::string getFlagRhoZ() const
    {
        return FlagRhoZ;
    }

    std::string getFlagKappaZ() const
    {
        return FlagKappaZ;
    }

    
    ///////////////////////////////////////////////////////////////////////////
    // get and set methods for class members

    Particle getLeptons(const StandardModel::lepton p) const
    {
        return leptons[p];
    }

    /**
     * @return the electromagnetic coupling
     */
    double getAle() const
    {
        return ale;
    }

    /**
     * @return @f$\Delta\alpha_\mathrm{had}^5(M_Z)@f$.
     */
    double getDAle5Mz() const
    {
        return dAle5Mz;
    }

    /**
     * @return The Fermi constant.
     */
    double getGF() const
    {
        return GF;
    }

    /**
     * @return The Higgs mass.
     */
    double getMHl() const
    {
        return mHl;
    }

    /**
     * @return Theoretical uncertainty in the approximate formula for @f$M_W@f$.
     */
    double getDelMw() const
    {
        return delMw;
    }

    /**
     * @return Theoretical uncertainty in the approximate formula for the leptonic weak mixing angle.
     */
    double getDelSin2th_l() const
    {
        return delSin2th_l;
    }

    /**
     * @return Theoretical uncertainty in the total width of Z boson.
     */
    double getDelGammaZ() const
    {
        return delGammaZ;
    }

    /**
     * @return The CKM matrix.
     */
    matrix<complex> getVCKM() const
    {
        return VCKM;
    }

    CKM getCKM() const
    {
        return myCKM;
    }

    double getLambda() const 
    {
        return lambda;
    }

    double getA() const
    {
        return A;
    }

    double getEtab() const
    {
        return etab;
    }
    double getRhob() const
    {
        return rhob;
    }

    /**
     * @return The PMNS matrix.
     */
    matrix<complex> getUPMNS() const
    {
        return UPMNS;
    }

    /**
     * @return The up Yukawa matrix.
     */
    matrix<complex> getYu() const
    {
        return Yu;
    }

    /**
     * @return The down Yukawa matrix.
     */
    matrix<complex> getYd() const
    {
        return Yd;
    }

    /**
     * @return The neutrino Yukawa matrix.
     */
    matrix<complex> getYn() const
    {
        return Yn;
    }

    /**
     * @return The charged lepton Yukawa matrix.
     */
    matrix<complex> getYe() const
    {
        return Ye;
    }

    /**
     * @return The Standard Model contribution to @f$ \Delta m_{K} @f$.
     */
    double getDmk() const
    {
        return Dmk;
    }

    /**
     * @return The Standard Model amplitude of the @f$ D^{0} - \bar{D}^{0} @f$ mixing.
     */
    double getSM_M12D() const
    {
        return SM_M12D;
    }

    double getMuw() const
    {
        return muw;
    }

    double getKbarEpsK() const
    {
        return KbarEpsK;
    }

    /**
     * @return The experimental value of @f$ \Delta M_{K}/(\Delta\Gamma_{K}/2) @f$.
     */
    double getphiEpsK() const
    {
        return phiEpsK;
    }

    /**
     * @return The experimental value of @f$ \Delta M_{K} @f$.
     */
    double getDeltaMK() const
    {
        return DeltaMK;
    }

    /**
     * @return The experimental value of @f$ \varepsilon_{K} @f$.
     */
    double getEpsK() const {
        return EpsK;
    }


    ///////////////////////////////////////////////////////////////////////////

    /**
     * @return The VEV.
     */
    virtual double v() const;

    /**
     * @return The W boson mass at tree level.
     */
    virtual double Mw_tree() const;

    /**
     * Computes the running electromagnetic coupling alpha(mu) in the on-shell
     * scheme, where the top-quark contribution is not included.
     * @param[in] mu A scale @f$\mu@f$ in GeV.
     * @param[in] order (=LO, FULLNLO)
     * @return @f$\alpha(\mu)@f$ in the on-shell scheme.
     */
    double ale_OS(const double mu, orders order=FULLNLO) const;


    ////////////////////////////////////////////////////////////////////////

    /**
     * @param[in] s An invariant mass squared.
     * @return The leptonic corrections to alpha at @f$M_Z@f$.
     */
    double DeltaAlphaLepton(const double s) const;

    /**
     * @return The sum of the leptonic and hadronic corrections to alpha at @f$M_Z@f$.
     */
    double DeltaAlphaL5q() const;

    /**
     * @return The total (leptonic+hadronic+top) corrections to alpha at @f$M_Z@f$.
     */
    double DeltaAlpha() const;

    /**
     * @return The electromagnetic coupling at @f$M_Z@f$, @f$\alpha(M_Z)@f$ in the on-shell scheme.
     */
    double alphaMz() const;

    /**
     * @return The W boson mass in the on-shell scheme.
     */
    virtual double Mw() const;

    /**
     * @return @f$M_W^2/M_Z^2@f$ in the on-shell scheme.
     */
    virtual double cW2() const;

    /**
     * @return @f$1-M_W^2/M_Z^2@f$ in the on-shell scheme.
     */
    virtual double sW2() const;

    /**
     * @return The total width of the W boson.
     */
    virtual double GammaW() const;


    ////////////////////////////////////////////////////////////////////////
    // CKM parameters

    // Angles
    double computeBeta() const;
    double computeGamma() const;
    double computeAlpha() const;
    double computeBetas() const;

    // Lambda_q
    complex computelamt() const;
    complex computelamc() const;
    complex computelamu() const;

    complex computelamt_d() const;
    complex computelamc_d() const;
    complex computelamu_d() const;

    complex computelamt_s() const;
    complex computelamc_s() const;
    complex computelamu_s() const;

    // Sides
    double getRt() const;
    double getRts() const;
    double getRb() const;


    ////////////////////////////////////////////////////////////////////////
protected:

    /**
     * @brief A set method to change the value of a parameter in the current model.
     * @param[in] name name of a parameter
     * @param[in] value the value to be assigned to the parameter specified by name
     */
    virtual void setParameter(const std::string name, const double& value);
    
    virtual void computeCKM();
    
    virtual void computeYukawas();

    EWSM* myEWSM;///< A pointer to an object of type EWSM.

    Particle leptons[6];///<
    CKM myCKM;///<
    matrix<complex> VCKM;///<
    matrix<complex> UPMNS;///< 
    matrix<complex> Yu;///<
    matrix<complex> Yd;///<
    matrix<complex> Yn;///<
    matrix<complex> Ye;///<

    // model parameters
    double ale;///<
    double dAle5Mz;///<
    double GF;///<
    double mHl;///<
    double delMw;///<
    double delSin2th_l;///<
    double delGammaZ;///<
    double muw;///<
    double lambda;///<
    double A;///<
    double rhob;///<
    double etab;///<
    double EpsK;///<
    double phiEpsK;///<
    double DeltaMK;///<
    double KbarEpsK;///<
    double Dmk;///<
    double SM_M12D;///< 
    
    
    ////////////////////////////////////////////////////////////////////////    
private:
    StandardModelMatching* myStandardModelMatching;///< A pointer to an object of type StandardModelMatching.

    bool FlagWithoutNonUniversalVC;///<
    bool FlagNoApproximateGammaZ;///<
    bool FlagNoApproximateSigmaH;///<
    bool FlagNoApproximateRl;///<
    bool FlagNoApproximateRc;///<
    bool FlagNoApproximateRb;///<
    std::string FlagMw;///<
    std::string FlagRhoZ;///<
    std::string FlagKappaZ;///<

    bool requireCKM;///<
    bool requireYe;///<
    bool requireYn;///< 

};

#endif	/* STANDARDMODEL_H */
