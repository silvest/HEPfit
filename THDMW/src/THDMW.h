/* 
 * Copyright (C) 2017 HEPfit Collaboration
 *
 *
 * For the licensing terms see doc/COPYING.
 */

#ifndef THDMW_H
#define	THDMW_H

#include "StandardModel.h"
#include "THDMWMatching.h"
#include "NPbase.h"

class THDMWcache; //forward reference to THDMWcache class

/**
 * @class THDMW
 * @ingroup THDMW
 * @brief A base class for @f$Z_2@f$ symmetric Two-Higgs-Doublet-Manohar-Wise models.
 * @author HEPfit Collaboration
 * @copyright GNU General Public License
 * @details We use the parametrization of Li/Valencia (see eq. (1) in 1703.03445 @cite Cheng:2017tbn ).
 *
 * 
 * @anchor THDMWParameters
 * <h3>%Model parameters</h3>
 *
 * The model parameters of %THDMW are summarized below.
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
 *   <td class="mod_name">%bma</td>
 *   <td class="mod_symb">@f$\beta - \alpha@f$</td>
 *   <td class="mod_desc">The difference of the mixing angles @f$\beta@f$ and @f$\alpha@f$.</td>
 * </tr>
 * <tr>
 *   <td class="mod_name">%THDMW_lambda1</td>
 *   <td class="mod_symb">@f$\lambda_{1}^{\text{THDMW}}@f$</td>
 *   <td class="mod_desc"> The potential parameter @f$\lambda_{1}@f$ using the notation of @cite Cheng:2016tlc (arxiv:1606.01298)  </td>
 * </tr>
 * <tr>
 *   <td class="mod_name">%THDMW_lambda2</td>
 *   <td class="mod_symb">@f$\lambda_{2}^{\text{THDMW}}@f$</td>
 *   <td class="mod_desc"> The potential parameter @f$\lambda_{2}@f$ using the notation of @cite Cheng:2016tlc (arxiv:1606.01298)  </td>
 * </tr>
 * <tr>
 *   <td class="mod_name">%THDMW_lambda3</td>
 *   <td class="mod_symb">@f$\lambda_{3}^{\text{THDMW}}@f$</td>
 *   <td class="mod_desc"> The potential parameter @f$\lambda_{3}@f$ using the notation of @cite Cheng:2016tlc (arxiv:1606.01298)  </td>
 * </tr>
 * <tr>
 *   <td class="mod_name">%THDMW_lambda4</td>
 *   <td class="mod_symb">@f$\lambda_{4}^{\text{THDMW}}@f$</td>
 *   <td class="mod_desc"> The potential parameter @f$\lambda_{4}@f$ using the notation of @cite Cheng:2016tlc (arxiv:1606.01298)  </td>
 * </tr>
 * <tr>
 *   <td class="mod_name">%THDMW_lambda5</td>
 *   <td class="mod_symb">@f$\lambda_{5}^{\text{THDMW}}@f$</td>
 *   <td class="mod_desc"> The potential parameter @f$\lambda_{5}@f$ using the notation of @cite Cheng:2016tlc (arxiv:1606.01298)  </td>
 * </tr>
 * <tr>
 *   <td class="mod_name">%THDMW_mS2</td>
 *   <td class="mod_symb">@f$m_{S}^{2}@f$</td>
 *   <td class="mod_desc"> The mass parameter for the colour scalars. Note that this is not, in general, the physical mass because it gets contributions from the vev </td>
 * </tr>
 * <tr>
 *   <td class="mod_name">%THDMW_mu1</td>
 *   <td class="mod_symb">@f$\mu_{1}^{\text{THDMW}}@f$</td>
 *   <td class="mod_desc"> The potential parameter @f$\mu_{1}@f$ using the notation of @cite Cheng:2016tlc (arxiv:1606.01298)  </td>
 * </tr>
 * <tr>
 *   <td class="mod_name">%THDMW_mu2</td>
 *   <td class="mod_symb">@f$\mu_{2}^{\text{THDMW}}@f$</td>
 *   <td class="mod_desc"> The potential parameter @f$\mu_{2}@f$ using the notation of @cite Cheng:2016tlc (arxiv:1606.01298)  </td>
 * </tr>
 * <tr>
 *   <td class="mod_name">%THDMW_mu3</td>
 *   <td class="mod_symb">@f$\mu_{3}^{\text{THDMW}}@f$</td>
 *   <td class="mod_desc"> The potential parameter @f$\mu_{3}@f$ using the notation of @cite Cheng:2016tlc (arxiv:1606.01298)  </td>
 * </tr>
 * <tr>
 *   <td class="mod_name">%THDMW_mu4</td>
 *   <td class="mod_symb">@f$\mu_{4}^{\text{THDMW}}@f$</td>
 *   <td class="mod_desc"> The potential parameter @f$\mu_{4}@f$ using the notation of @cite Cheng:2016tlc (arxiv:1606.01298)  </td>
 * </tr>
 * <tr>
 *   <td class="mod_name">%THDMW_mu5</td>
 *   <td class="mod_symb">@f$\mu_{5}^{\text{THDMW}}@f$</td>
 *   <td class="mod_desc"> The potential parameter @f$\mu_{5}@f$ using the notation of @cite Cheng:2016tlc (arxiv:1606.01298)  </td>
 * </tr>
 * <tr>
 *   <td class="mod_name">%THDMW_mu6</td>
 *   <td class="mod_symb">@f$\mu_{6}^{\text{THDMW}}@f$</td>
 *   <td class="mod_desc"> The potential parameter @f$\mu_{6}@f$ using the notation of @cite Cheng:2016tlc (arxiv:1606.01298)  </td>
 * </tr>
 * <tr>
 *   <td class="mod_name">%THDMW_nu1</td>
 *   <td class="mod_symb">@f$nu_{1}^{\text{THDMW}}@f$</td>
 *   <td class="mod_desc"> The potential parameter @f$nu_{1}@f$ using the notation of @cite Cheng:2016tlc (arxiv:1606.01298)  </td>
 * </tr>
 * <tr>
 *   <td class="mod_name">%THDMW_nu2</td>
 *   <td class="mod_symb">@f$nu_{2}^{\text{THDMW}}@f$</td>
 *   <td class="mod_desc"> The potential parameter @f$nu_{2}@f$ using the notation of @cite Cheng:2016tlc (arxiv:1606.01298)  </td>
 * </tr>
 * <tr>
 *   <td class="mod_name">%THDMW_nu3</td>
 *   <td class="mod_symb">@f$nu_{3}^{\text{THDMW}}@f$</td>
 *   <td class="mod_desc"> The potential parameter @f$nu_{3}@f$ using the notation of @cite Cheng:2016tlc (arxiv:1606.01298)  </td>
 * </tr>
 * <tr>
 *   <td class="mod_name">%THDMW_nu4</td>
 *   <td class="mod_symb">@f$nu_{4}^{\text{THDMW}}@f$</td>
 *   <td class="mod_desc"> The potential parameter @f$nu_{4}@f$ using the notation of @cite Cheng:2016tlc (arxiv:1606.01298)  </td>
 * </tr>
 * <tr>
 *   <td class="mod_name">%THDMW_nu5</td>
 *   <td class="mod_symb">@f$nu_{5}^{\text{THDMW}}@f$</td>
 *   <td class="mod_desc"> The potential parameter @f$nu_{5}@f$ using the notation of @cite Cheng:2016tlc (arxiv:1606.01298)  </td>
 * </tr>
 * <tr>
 *   <td class="mod_name">%THDMW_omega1</td>
 *   <td class="mod_symb">@f$\omega_{1}^{\text{THDMW}}@f$</td>
 *   <td class="mod_desc"> The potential parameter @f$\omega_{1}@f$ using the notation of @cite Cheng:2016tlc (arxiv:1606.01298)  </td>
 * </tr>
 * <tr>
 *   <td class="mod_name">%THDMW_omega2</td>
 *   <td class="mod_symb">@f$\omega_{2}^{\text{THDMW}}@f$</td>
 *   <td class="mod_desc"> The potential parameter @f$\omega_{2}@f$ using the notation of @cite Cheng:2016tlc (arxiv:1606.01298)  </td>
 * </tr>
 * <tr>
 *   <td class="mod_name">%THDMW_omega3</td>
 *   <td class="mod_symb">@f$\omega_{3}^{\text{THDMW}}@f$</td>
 *   <td class="mod_desc"> The potential parameter @f$\omega_{3}@f$ using the notation of @cite Cheng:2016tlc (arxiv:1606.01298)  </td>
 * </tr>
 * <tr>
 *   <td class="mod_name">%THDMW_omega4</td>
 *   <td class="mod_symb">@f$\omega_{4}^{\text{THDMW}}@f$</td>
 *   <td class="mod_desc"> The potential parameter @f$\omega_{4}@f$ using the notation of @cite Cheng:2016tlc (arxiv:1606.01298)  </td>
 * </tr>
 * <tr>
 *   <td class="mod_name">%THDMW_omega5</td>
 *   <td class="mod_symb">@f$\omega_{5}^{\text{THDMW}}@f$</td>
 *   <td class="mod_desc"> The potential parameter @f$\omega_{5}@f$ using the notation of @cite Cheng:2016tlc (arxiv:1606.01298)  </td>
 * </tr>
 * <tr>
 *   <td class="mod_name">%THDMW_kappa1</td>
 *   <td class="mod_symb">@f$\kappa_{1}^{\text{THDMW}}@f$</td>
 *   <td class="mod_desc"> The potential parameter @f$\kappa_{1}@f$ using the notation of @cite Cheng:2016tlc (arxiv:1606.01298)  </td>
 * </tr>
 * <tr>
 *   <td class="mod_name">%THDMW_kappa2</td>
 *   <td class="mod_symb">@f$\kappa_{2}^{\text{THDMW}}@f$</td>
 *   <td class="mod_desc"> The potential parameter @f$\kappa_{2}@f$ using the notation of @cite Cheng:2016tlc (arxiv:1606.01298)  </td>
 * </tr>
 * <tr>
 *   <td class="mod_name">%THDMW_kappa3</td>
 *   <td class="mod_symb">@f$\kappa_{3}^{\text{THDMW}}@f$</td>
 *   <td class="mod_desc"> The potential parameter @f$\kappa_{3}@f$ using the notation of @cite Cheng:2016tlc (arxiv:1606.01298)  </td>
 * </tr>
 * <tr>
 *   <td class="mod_name">%THDMW_etaU</td>
 *   <td class="mod_symb">@f$\eta_{U}^{\text{THDMW}}@f$</td>
 *   <td class="mod_desc"> The up-type Yukawa proportionality factor for the colour scalars  </td>
 * </tr>
 * <tr>
 *   <td class="mod_name">%THDMW_etaD</td>
 *   <td class="mod_symb">@f$\eta_{D}^{\text{THDMW}}@f$</td>
 *   <td class="mod_desc"> The down-type Yukawa proportionality factor for the colour scalars  </td>
 * </tr>
 * <tr>
 *   <td class="mod_name">%THDMW_rho_b</td>
 *   <td class="mod_symb">@f$\rho_{b}^{SM}@f$</td>
 *   <td class="mod_desc"> The radiative parameter @f$\rho_{b}@f$ which corrects the @f$Zq\bar{q}@f$ couplings. This is, the parameter @f$\rho_{b}@f$ of arxiv:0509008  </td>
 * </tr>
 * <tr>
 *   <td class="mod_name">%THDMW_S_b</td>
 *   <td class="mod_symb">@f$S_{b}@f$</td>
 *   <td class="mod_desc"> The sum of the @f$Zq\bar{q}@f$ couplings for the up, down and charm quarks including the QED corrections. Look at arxiv:1002.1071 for more details  </td>
 * </tr>
 * <tr>
 *   <td class="mod_name">%Q_THDMW</td>
 *   <td class="mod_symb">@f$Q_{\text{THDMW}}@f$</td>
 *   <td class="mod_desc">The THDM scale.</td>
 * </tr>
 * <tr>
 *   <td class="mod_name">%RpepsTHDMW</td>
 *   <td class="mod_symb">@f$\varepsilon_{\text{R'}}^{\text{THDMW}}@f$</td>
 *   <td class="mod_desc">Minimal value for which the R' perturbativity criterion should be applied for the unitarity bounds.</td>
 * </tr>
 * <tr>
 *   <td class="mod_name">%NLOuniscaleTHDMW</td>
 *   <td class="mod_symb">@f$Q_{\text{min}}^{\text{NLOuni}}@f$</td>
 *   <td class="mod_desc">Minimal scale at which the NLO unitarity conditions are checked.</td>
 * </tr>
 * </table>
 * 
 * 
 * 
 * @anchor THDMWFlags
 * <h3>%Model flags</h3>
 *
 * <table class="model">
 * <tr>
 *   <th>Label</th>
 *   <th>Value</th>
 *   <th>Description</th>
 * </tr>
 * <tr>
 *   <td class="mod_name">%modelTypeflag</td>
 *   <td class="mod_valu">custodial1 / custodial2 / ManoharWise / custodialMW </td>
 *   <td class="mod_desc">This flag determines the type of model we are using. </td>
 * </tr>
 * <tr>
 *   <td class="mod_name">%RGEorder</td>
 *   <td class="mod_valu">LO / approxNLO </td>
 *   <td class="mod_desc">This flag determines the order in perturbation theory of the renormalization group equations.</td>
 * </tr>
 * </table>
 *
 */
class THDMW: public NPbase {
public:

    static const int NTHDMWvars = 34;
    static const std::string THDMWvars[NTHDMWvars];
    
    /**
     * @brief %THDMW constructor
     */
    THDMW();
    
    /**
     * @brief %THDMW destructor
     */
    ~THDMW();
    
    /**
     * @brief A method to initialize the model.
     * @details This method, called via InputParser::ReadParameters(), allocates
     * memory to the pointers defined in the current class.
     * @return a boolean that is true if model initialization is successful
     */
    virtual bool InitializeModel();
    
    /**
     * @brief Initializes the %THDMW parameters found in the argument.
     * @param[in] DPars a map containing the parameters (all as double) to be used in Monte Carlo
     */
    virtual bool Init(const std::map<std::string, double>& DPars);
    
    /**
     * @brief The pre-update method for %THDMW
     * @return a boolean that is true if the execution is successful
     */
    virtual bool PreUpdate();
    
    /**
     * @brief The update method for %THDMW.
     * @details This method updates all the model parameters with given DPars.
     * @param[in] DPars a map of the parameters that are being updated in the Monte Carlo run
     * @return a boolean that is true if the execution is successful
     */
    virtual bool Update(const std::map<std::string, double>& DPars);
    
    /**
     * @brief The post-update method for %THDMW.
     * @details This method runs all the procedures that are need to be executed
     * after the model is successfully updated.
     * @return a boolean that is true if the execution is successful
     */
    virtual bool PostUpdate();
    
    /**
     * @brief A method to check if all the mandatory parameters for %THDMW
     * have been provided in model initialization.
     * @param[in] DPars a map of the parameters that are being updated in the Monte Carlo run
     * (including parameters that are varied and those that are held constant)
     * @return a boolean that is true if the execution is successful
     */
    virtual bool CheckParameters(const std::map<std::string, double>& DPars);

    /**
     * @brief A get method to access the member reference of type THDMWMatching.
     * @return a reference to a THDMWMatching object
     */
    virtual THDMWMatching& getMatching() const
    {
        return THDMWM.getObj();
    }

    ///////////////////////////////////////////////////////////////////////////
    // Flags

    /**
     * @brief A method to set a string flag of %THDMW.
     * @param[in] name name of a model flag
     * @param[in] value the string to be assigned to the flag specified by name
     * @return a boolean that is true if the execution is successful
     */
    virtual bool setFlagStr(const std::string name, const std::string value);

    /**
     * @brief A method get the THDMWCache
     * @return a object of the type %THDMWCache
     */
    THDMWcache* getMyTHDMWCache() const
    {
        return myTHDMWcache;
    }

    ///////////////////////////////////////////////////////////////////////////

    /**
     * @brief A getter for the %THDMW model type
     * @return %THDMW model type
     */
    std::string getModelTypeTHDMWflag() const {
        return flag_modelTHDMW;
    }

    /**
     * @brief A getter for the switch for NLO RGE and approximate NLO RGE
     * @return Switch for NLO RGE and approximate NLO RGE
     */
    std::string getRGEorderflag() const {
        return flag_RGEorder;
    }

    /**
     * @brief A getter for @f$\log_{10}(\tan \beta)@f$
     * @return @f$\log_{10}(\tan \beta)@f$
     */
    const double getTHDMW_logtb() const {
        return THDMW_logtb;
    }

    /**
     * @brief A getter for @f$\tan \beta@f$
     * @return @f$\tan \beta@f$
     */
    const double getTHDMW_tanb() const {
        return THDMW_tanb;
    }

    /**
     * @brief A getter for @f$\sin \beta@f$
     * @return @f$\sin \beta@f$
     */
    const double getTHDMW_sinb() const {
        return THDMW_sinb;
    }

    /**
     * @brief A getter for @f$\cos \beta@f$
     * @return @f$\cos \beta@f$
     */
    const double getTHDMW_cosb() const {
        return THDMW_cosb;
    }

    /**
     * @brief A getter for @f$\beta-\alpha@f$
     * @return @f$\beta-\alpha@f$
     */
    const double getTHDMW_bma() const {
        return THDMW_bma;
    }

    /**
     * @brief A getter for @f$\sin(\beta-\alpha)@f$
     * @return @f$\sin(\beta-\alpha)@f$
     */
    const double getTHDMW_sin_ba() const {
        return THDMW_sin_ba;
    }

    /**
     * @brief A getter for @f$\cos \alpha@f$
     * @return @f$\cos \alpha@f$
     */
    const double getTHDMW_cosa() const{
        return cos(atan(pow(10.,THDMW_logtb))-THDMW_bma);
    }

    /**
     * @brief A getter for @f$\sin \alpha@f$
     * @return @f$\sin \alpha@f$
     */
    const double getTHDMW_sina() const{
        return sin(atan(pow(10.,THDMW_logtb))-THDMW_bma);
    }

    /**
     * @brief A getter for @f$\lambda_1@f$
     * @return @f$\lambda_1@f$
     */
    const double getTHDMW_lambda1() const {
        return THDMW_lambda1;
    }

    /**
     * @brief A getter for @f$\lambda_2@f$
     * @return @f$\lambda_2@f$
     */
    const double getTHDMW_lambda2() const {
        if (flag_modelTHDMW == "ManoharWise" || flag_modelTHDMW == "custodialMW") {
            return 0.0;
        }
        else {
            return THDMW_lambda2;
        }
    }

    /**
     * @brief A getter for @f$\lambda_3@f$
     * @return @f$\lambda_3@f$
     */
    const double getTHDMW_lambda3() const {
        if (flag_modelTHDMW == "ManoharWise" || flag_modelTHDMW == "custodialMW") {
            return 0.0;
        }
        else {
            return THDMW_lambda3;
        }
    }

    /**
     * @brief A getter for @f$\lambda_4@f$
     * @return @f$\lambda_4@f$
     */
    const double getTHDMW_lambda4() const {
        if (flag_modelTHDMW == "ManoharWise" || flag_modelTHDMW == "custodialMW") {
            return 0.0;
        }
        else {
            return THDMW_lambda4;
        }
    }

    /**
     * @brief A getter for @f$\lambda_5@f$
     * @return @f$\lambda_5@f$
     */
    const double getTHDMW_lambda5() const {
        if (flag_modelTHDMW == "custodial1") {
            return THDMW_lambda4;
        }
        else if (flag_modelTHDMW == "ManoharWise" || flag_modelTHDMW == "custodialMW") {
            return 0.0;
        }
        else {
            return THDMW_lambda5;
        }
    }

    /**
     * @brief A getter for @f$m_S^2@f$
     * @return @f$m_S^2@f$
     */
    const double getTHDMW_mS2() const {
        return THDMW_mS2;
    }

    /**
     * @brief A getter for @f$\mu_1@f$
     * @return @f$\mu_1@f$
     */
    const double getTHDMW_mu1() const {
        return THDMW_mu1;
    }

    /**
     * @brief A getter for @f$\mu_2@f$
     * @return @f$\mu_2@f$
     */
    const double getTHDMW_mu2() const {
        if (flag_modelTHDMW == "custodial1" || flag_modelTHDMW == "custodialMW") {
            return THDMW_mu1;
        }
        else {
            return THDMW_mu2;
        }
    }

    /**
     * @brief A getter for @f$\mu_3@f$
     * @return @f$\mu_3@f$
     */
    const double getTHDMW_mu3() const {
        return THDMW_mu3;
    }
    /**
     * @brief A getter for @f$\mu_4@f$
     * @return @f$\mu_4@f$
     */
    const double getTHDMW_mu4() const {
        return THDMW_mu4;
    }

    /**
     * @brief A getter for @f$\mu_5@f$
     * @return @f$\mu_5@f$
     */
    const double getTHDMW_mu5() const {
        if (flag_modelTHDMW == "custodial1" || flag_modelTHDMW == "custodialMW") {
            return THDMW_mu4;
        }
        else {
            return THDMW_mu5;
        }
    }

    /**
     * @brief A getter for @f$\mu_6@f$
     * @return @f$\mu_6@f$
     */
    const double getTHDMW_mu6() const {
        if (flag_modelTHDMW == "custodial1" || flag_modelTHDMW == "custodialMW") {
            return 2.0*THDMW_mu1;
        }
        else {
            return THDMW_mu6;
        }
    }

    /**
     * @brief A getter for @f$\nu_1@f$
     * @return @f$\nu_1@f$
     */
    const double getTHDMW_nu1() const {
        return THDMW_nu1;
    }

    /**
     * @brief A getter for @f$\nu_2@f$
     * @return @f$\nu_2@f$
     */
    const double getTHDMW_nu2() const {
        return THDMW_nu2;
    }

    /**
     * @brief A getter for @f$\nu_3@f$
     * @return @f$\nu_3@f$
     */
    const double getTHDMW_nu3() const {
        if (flag_modelTHDMW == "custodial1" || flag_modelTHDMW == "custodialMW") {
            return 0.5*THDMW_nu2;
        }
        else {
            return THDMW_nu3;
        }
    }
    /**
     * @brief A getter for @f$\nu_4@f$
     * @return @f$\nu_4@f$
     */
    const double getTHDMW_nu4() const {
        return THDMW_nu4;
    }

    /**
     * @brief A getter for @f$\nu_5@f$
     * @return @f$\nu_5@f$
     */
    const double getTHDMW_nu5() const {
        if (flag_modelTHDMW == "custodial1" || flag_modelTHDMW == "custodialMW") {
            return THDMW_nu4;
        }
        else {
            return THDMW_nu5;
        }
    }

    /**
     * @brief A getter for @f$\omega_1@f$
     * @return @f$\omega_1@f$
     */
    const double getTHDMW_omega1() const {
        if (flag_modelTHDMW == "ManoharWise" || flag_modelTHDMW == "custodialMW") {
            return 0.0;
        }
        else {
            return THDMW_omega1;
        }
    }

    /**
     * @brief A getter for @f$\omega_2@f$
     * @return @f$\omega_2@f$
     */
    const double getTHDMW_omega2() const {
        if (flag_modelTHDMW == "ManoharWise" || flag_modelTHDMW == "custodialMW") {
            return 0.0;
        }
        else {
            return THDMW_omega2;
        }
    }

    /**
     * @brief A getter for @f$\omega_3@f$
     * @return @f$\omega_3@f$
     */
    const double getTHDMW_omega3() const {
        if (flag_modelTHDMW == "custodial1") {
            return 0.5*THDMW_omega2;
        }
        else if (flag_modelTHDMW == "ManoharWise" || flag_modelTHDMW == "custodialMW") {
            return 0.0;
        }
        else {
            return THDMW_omega3;
        }
    }

    /**
     * @brief A getter for @f$\omega_4@f$
     * @return @f$\omega_4@f$
     */
    const double getTHDMW_omega4() const {
        if (flag_modelTHDMW == "ManoharWise" || flag_modelTHDMW == "custodialMW") {
            return 0.0;
        }
        else {
            return THDMW_omega4;
        }
    }

    /**
     * @brief A getter for  @f$\omega_5@f$
     * @return @f$\omega_5@f$
     */
    const double getTHDMW_omega5() const {
        if (flag_modelTHDMW == "custodial1") {
            return THDMW_omega4;
        }
        else if (flag_modelTHDMW == "ManoharWise" || flag_modelTHDMW == "custodialMW") {
            return 0.0;
        }
        else {
            return THDMW_omega5;
        }
    }

    /**
     * @brief A getter for  @f$\kappa_1@f$
     * @return @f$\kappa_1@f$
     */
    const double getTHDMW_kappa1() const {
        if (flag_modelTHDMW == "ManoharWise" || flag_modelTHDMW == "custodialMW") {
            return 0.0;
        }
        else {
            return THDMW_kappa1;
        }
    }

    /**
     * @brief A getter for @f$\kappa_2@f$
     * @return @f$\kappa_2@f$
     */
    const double getTHDMW_kappa2() const {
        if (flag_modelTHDMW == "ManoharWise" || flag_modelTHDMW == "custodialMW") {
            return 0.0;
        }
        else {
            return THDMW_kappa2;
        }
    }

    /**
     * @brief A getter for @f$\kappa_3@f$ 
     * @return @f$\kappa_3@f$
     */
    const double getTHDMW_kappa3() const {
        if (flag_modelTHDMW == "custodial1") {
            return THDMW_kappa2;
        }
        else if (flag_modelTHDMW == "ManoharWise" || flag_modelTHDMW == "custodialMW") {
            return 0.0;
        }
        else {
            return THDMW_kappa3;
        }
    }

    /**
     * @brief A getter for @f$\eta_U@f$
     * @return @f$\eta_U@f$
     */
    const double getTHDMW_etaU() const {
        return THDMW_etaU;
    }

    /**
     * @brief A getter for @f$\eta_D@f$
     * @return @f$\eta_D@f$
     */
    const double getTHDMW_etaD() const {
        return THDMW_etaD;
    }
    
    
    /**
     * @brief A getter for @f$\eta_U@f$
     * @return @f$\eta_U@f$
     */
    const double getTHDMW_rho_b() const {
        return THDMW_rho_b;
    }

    /**
     * @brief A getter for @f$\eta_D@f$
     * @return @f$\eta_D@f$
     */
    const double getTHDMW_S_b() const {
        return THDMW_S_b;
    }
    
    

    /**
     * @brief A getter for the %THDMW scale
     * @return %THDMW scale
     */
    const double getQ_THDMW() const {
        return Q_THDMW;
    }

    /**
     * @brief A getter for the minimal R' value
     * @return Minimal R' value
     */
    const double getRpepsTHDMW() const {
        return RpepsTHDMW;
    }

    /**
     * @brief A getter for the minimal NLO unitarity check scale
     * @return Minimal NLO unitarity check scale
     */
    const double getNLOuniscaleTHDMW() const {
        return NLOuniscaleTHDMW;
    }

    virtual const double muggH(const double sqrt_s) const;
    virtual const double muVBF(const double sqrt_s) const;
    virtual const double mueeWBF(const double sqrt_s) const;
    virtual const double muWH(const double sqrt_s) const;
    virtual const double muZH(const double sqrt_s) const;
    virtual const double mueeZH(const double sqrt_s) const;
    virtual const double muVH(const double sqrt_s) const;
    virtual const double muVBFpVH(const double sqrt_s) const;
    virtual const double muttH(const double sqrt_s) const;
    virtual const double GammaTotal() const;
    virtual const double BrHggRatio() const;
    virtual const double BrHWWRatio() const;
    virtual const double BrHZZRatio() const;
    virtual const double BrHZgaRatio() const;
    virtual const double BrHgagaRatio() const;
    virtual const double BrHmumuRatio() const;
    virtual const double BrHtautauRatio() const;
    virtual const double BrHccRatio() const;
    virtual const double BrHbbRatio() const;
    virtual const double muggHgaga(const double sqrt_s) const;
    virtual const double muVBFHgaga(const double sqrt_s) const;
    virtual const double muVHgaga(const double sqrt_s) const;
    virtual const double muttHgaga(const double sqrt_s) const;
    virtual const double muggHZZ(const double sqrt_s) const;
    virtual const double muVBFHZZ(const double sqrt_s) const;
    virtual const double muVHZZ(const double sqrt_s) const;
    virtual const double muttHZZ(const double sqrt_s) const;
    virtual const double muggHWW(const double sqrt_s) const;
    virtual const double muVBFHWW(const double sqrt_s) const;
    virtual const double muVHWW(const double sqrt_s) const;
    virtual const double muttHWW(const double sqrt_s) const;
    virtual const double muggHtautau(const double sqrt_s) const;
    virtual const double muVBFHtautau(const double sqrt_s) const;
    virtual const double muVHtautau(const double sqrt_s) const;
    virtual const double muttHtautau(const double sqrt_s) const;
    virtual const double muggHbb(const double sqrt_s) const;
    virtual const double muVBFHbb(const double sqrt_s) const;
    virtual const double muVHbb(const double sqrt_s) const;
    virtual const double muttHbb(const double sqrt_s) const;
    virtual const double muppHmumu(const double sqrt_s) const;
    virtual const double muppHZga(const double sqrt_s) const;
    virtual const double computeGammaTotalRatio() const;

    /*
     * @brief Method overriding the NPbase Mw with the tree-level W mass
     */
    virtual const double Mw() const;

protected:

    virtual void setParameter(const std::string, const double&);

    /**
     * @brief A method to check if the model type name in string form is valid.
     * @param[in] modeltype THDMW model type name
     * @return a boolean that is true if the model type name is valid
     */
    bool checkmodelTypeTHDMW(const std::string modeltype) const
    {
        if (modeltype.compare("custodial1") == 0)
            return true;
        else if (modeltype.compare("ManoharWise") == 0)
            return true;
        else if (modeltype.compare("custodialMW") == 0)
            return true;
        else if (modeltype.compare("custodial2") == 0)
            return true;
        else
            return false;
    }

    mutable Matching<THDMWMatching,THDMW> THDMWM; ///< An object of type Matching.

    /**
     * @brief A method to check if the RGE order name in string form is valid.
     * @param[in] RGEorder THDMW RGE order
     * @return a boolean that is true if the RGE order string is valid
     */
    bool checkRGEorder(const std::string RGEorder) const
    {
        if (RGEorder.compare("LO") == 0
                || RGEorder.compare("approxNLO") == 0)
            return true;
        else
            return false;
    }

private:

    THDMWcache* myTHDMWcache;

    double THDMW_logtb, THDMW_tanb, THDMW_sinb, THDMW_cosb, THDMW_bma, THDMW_sin_ba;
    double THDMW_lambda1, THDMW_lambda2, THDMW_lambda3, THDMW_lambda4, THDMW_lambda5;
    double THDMW_mS2, THDMW_mu1, THDMW_mu2, THDMW_mu3, THDMW_mu4, THDMW_mu5, THDMW_mu6;
    double THDMW_nu1, THDMW_nu2, THDMW_nu3, THDMW_nu4, THDMW_nu5;
    double THDMW_omega1, THDMW_omega2, THDMW_omega3, THDMW_omega4, THDMW_omega5;
    double THDMW_kappa1, THDMW_kappa2, THDMW_kappa3;
    double THDMW_etaU, THDMW_etaD;
    double THDMW_rho_b, THDMW_S_b;
    double Q_THDMW, RpepsTHDMW, NLOuniscaleTHDMW;
    std::string flag_modelTHDMW, flag_RGEorder;
};

#endif	/* THDMW_H */
