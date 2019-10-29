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
 * @details We use the parametrization of Li/Valencia (see eq. (1) in 1703.03445).
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
 *   <td class="mod_symb">@f$lambda_{1}^{\text{THDMW}}@f$</td>
 *   <td class="mod_desc"> The potential parameter @f$lambda_{1}@f$ using the notation of @cite Cheng:2016tlc (arxiv:1606.01298)  </td>
 * </tr>
 * <tr>
 *   <td class="mod_name">%THDMW_lambda2</td>
 *   <td class="mod_symb">@f$lambda_{2}^{\text{THDMW}}@f$</td>
 *   <td class="mod_desc"> The potential parameter @f$lambda_{2}@f$ using the notation of @cite Cheng:2016tlc (arxiv:1606.01298)  </td>
 * </tr>
 * <tr>
 *   <td class="mod_name">%THDMW_lambda3</td>
 *   <td class="mod_symb">@f$lambda_{3}^{\text{THDMW}}@f$</td>
 *   <td class="mod_desc"> The potential parameter @f$lambda_{3}@f$ using the notation of @cite Cheng:2016tlc (arxiv:1606.01298)  </td>
 * </tr>
 * <tr>
 *   <td class="mod_name">%THDMW_lambda4</td>
 *   <td class="mod_symb">@f$lambda_{4}^{\text{THDMW}}@f$</td>
 *   <td class="mod_desc"> The potential parameter @f$lambda_{4}@f$ using the notation of @cite Cheng:2016tlc (arxiv:1606.01298)  </td>
 * </tr>
 * <tr>
 *   <td class="mod_name">%THDMW_lambda5</td>
 *   <td class="mod_symb">@f$lambda_{5}^{\text{THDMW}}@f$</td>
 *   <td class="mod_desc"> The potential parameter @f$lambda_{5}@f$ using the notation of @cite Cheng:2016tlc (arxiv:1606.01298)  </td>
 * </tr>
 * <tr>
 *   <td class="mod_name">%THDMW_mS2</td>
 *   <td class="mod_symb">@f$m_{S}^{2}@f$</td>
 *   <td class="mod_desc"> The mass parameter for the colour scalars. Note that this is not, in general, the physical mass because it gets contributions from the vev </td>
 * </tr>
 * <tr>
 *   <td class="mod_name">%THDMW_mu1</td>
 *   <td class="mod_symb">@f$mu_{1}^{\text{THDMW}}@f$</td>
 *   <td class="mod_desc"> The potential parameter @f$mu_{1}@f$ using the notation of @cite Cheng:2016tlc (arxiv:1606.01298)  </td>
 * </tr>
 * <tr>
 *   <td class="mod_name">%THDMW_mu2</td>
 *   <td class="mod_symb">@f$mu_{2}^{\text{THDMW}}@f$</td>
 *   <td class="mod_desc"> The potential parameter @f$mu_{2}@f$ using the notation of @cite Cheng:2016tlc (arxiv:1606.01298)  </td>
 * </tr>
 * <tr>
 *   <td class="mod_name">%THDMW_mu3</td>
 *   <td class="mod_symb">@f$mu_{3}^{\text{THDMW}}@f$</td>
 *   <td class="mod_desc"> The potential parameter @f$mu_{3}@f$ using the notation of @cite Cheng:2016tlc (arxiv:1606.01298)  </td>
 * </tr>
 * <tr>
 *   <td class="mod_name">%THDMW_mu4</td>
 *   <td class="mod_symb">@f$mu_{4}^{\text{THDMW}}@f$</td>
 *   <td class="mod_desc"> The potential parameter @f$mu_{4}@f$ using the notation of @cite Cheng:2016tlc (arxiv:1606.01298)  </td>
 * </tr>
 * <tr>
 *   <td class="mod_name">%THDMW_mu5</td>
 *   <td class="mod_symb">@f$mu_{5}^{\text{THDMW}}@f$</td>
 *   <td class="mod_desc"> The potential parameter @f$mu_{5}@f$ using the notation of @cite Cheng:2016tlc (arxiv:1606.01298)  </td>
 * </tr>
 * <tr>
 *   <td class="mod_name">%THDMW_mu6</td>
 *   <td class="mod_symb">@f$mu_{6}^{\text{THDMW}}@f$</td>
 *   <td class="mod_desc"> The potential parameter @f$mu_{6}@f$ using the notation of @cite Cheng:2016tlc (arxiv:1606.01298)  </td>
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
 *   <td class="mod_symb">@f$omega_{1}^{\text{THDMW}}@f$</td>
 *   <td class="mod_desc"> The potential parameter @f$omega_{1}@f$ using the notation of @cite Cheng:2016tlc (arxiv:1606.01298)  </td>
 * </tr>
 * <tr>
 *   <td class="mod_name">%THDMW_omega2</td>
 *   <td class="mod_symb">@f$omega_{2}^{\text{THDMW}}@f$</td>
 *   <td class="mod_desc"> The potential parameter @f$omega_{2}@f$ using the notation of @cite Cheng:2016tlc (arxiv:1606.01298)  </td>
 * </tr>
 * <tr>
 *   <td class="mod_name">%THDMW_omega3</td>
 *   <td class="mod_symb">@f$omega_{3}^{\text{THDMW}}@f$</td>
 *   <td class="mod_desc"> The potential parameter @f$omega_{3}@f$ using the notation of @cite Cheng:2016tlc (arxiv:1606.01298)  </td>
 * </tr>
 * <tr>
 *   <td class="mod_name">%THDMW_omega4</td>
 *   <td class="mod_symb">@f$omega_{4}^{\text{THDMW}}@f$</td>
 *   <td class="mod_desc"> The potential parameter @f$omega_{4}@f$ using the notation of @cite Cheng:2016tlc (arxiv:1606.01298)  </td>
 * </tr>
 * <tr>
 *   <td class="mod_name">%THDMW_omega5</td>
 *   <td class="mod_symb">@f$omega_{5}^{\text{THDMW}}@f$</td>
 *   <td class="mod_desc"> The potential parameter @f$omega_{5}@f$ using the notation of @cite Cheng:2016tlc (arxiv:1606.01298)  </td>
 * </tr>
 * <tr>
 *   <td class="mod_name">%THDMW_kappa1</td>
 *   <td class="mod_symb">@f$kappa_{1}^{\text{THDMW}}@f$</td>
 *   <td class="mod_desc"> The potential parameter @f$kappa_{1}@f$ using the notation of @cite Cheng:2016tlc (arxiv:1606.01298)  </td>
 * </tr>
 * <tr>
 *   <td class="mod_name">%THDMW_kappa2</td>
 *   <td class="mod_symb">@f$kappa_{2}^{\text{THDMW}}@f$</td>
 *   <td class="mod_desc"> The potential parameter @f$kappa_{2}@f$ using the notation of @cite Cheng:2016tlc (arxiv:1606.01298)  </td>
 * </tr>
 * <tr>
 *   <td class="mod_name">%THDMW_kappa3</td>
 *   <td class="mod_symb">@f$kappa_{3}^{\text{THDMW}}@f$</td>
 *   <td class="mod_desc"> The potential parameter @f$kappa_{3}@f$ using the notation of @cite Cheng:2016tlc (arxiv:1606.01298)  </td>
 * </tr>
 * <tr>
 *   <td class="mod_name">%THDMW_etaU</td>
 *   <td class="mod_symb">@f$eta_{U}^{\text{THDMW}}@f$</td>
 *   <td class="mod_desc"> The up-type Yukawa proportionality factor for the colour scalars  </td>
 * </tr>
 * <tr>
 *   <td class="mod_name">%THDMW_etaD</td>
 *   <td class="mod_symb">@f$eta_{D}^{\text{THDMW}}@f$</td>
 *   <td class="mod_desc"> The down-type Yukawa proportionality factor for the colour scalars  </td>
 * </tr>
 * <tr>
 *   <td class="mod_name">%THDMW_rho_b</td>
 *   <td class="mod_symb">@f$rho_{b}^{SM}@f$</td>
 *   <td class="mod_desc"> The radiative parameter @f$rho_{b}@f$ which corrects the @f$Zqq\bar@f$ couplings. This is, the parameter @f$rho_{b}@f$ of arxiv:0509008  </td>
 * </tr>
 * <tr>
 *   <td class="mod_name">%THDMW_S_b</td>
 *   <td class="mod_symb">@f$S_{b}@f$</td>
 *   <td class="mod_desc"> The sum of the @f$Zqq\bar@f$ couplings for the up, down and charm quarks including the QED corrections. Look at arxiv:1002.1071 for more details  </td>
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
 * @anchor THDMFlags
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
     * @brief THDMW constructor
     */
    THDMW();
    
    /**
     * @brief THDMW destructor
     */
    ~THDMW();
    
    virtual bool InitializeModel();
    
    virtual bool Init(const std::map<std::string, double>& DPars);
    
    virtual bool PreUpdate();
    
    virtual bool Update(const std::map<std::string, double>& DPars);
    
    virtual bool PostUpdate();
    
    virtual bool CheckParameters(const std::map<std::string, double>& DPars);

    /**
     * @brief A get method to access the member reference of type StandardModelMatching.
     * @return a reference to a StandardModelMatching object
     */
    virtual THDMWMatching& getMatching() const
    {
        return THDMWM.getObj();
    }

    ///////////////////////////////////////////////////////////////////////////
    // Flags

    virtual bool setFlagStr(const std::string name, const std::string value);
//    virtual bool setFlag(const std::string, const bool);

    THDMWcache* getMyTHDMWCache() const
    {
        return myTHDMWcache;
    }

    ///////////////////////////////////////////////////////////////////////////

    /**
     *
     * @return THDMW model type
     */
    std::string getModelTypeTHDMWflag() const {
        return flag_modelTHDMW;
    }

    /**
     *
     * @return Switch for NLO RGE and approximate NLO RGE
     */
    std::string getRGEorderflag() const {
        return flag_RGEorder;
    }

    /**
     *
     * @return @f$\log_{10}(\tan \beta)@f$
     */
    double getTHDMW_logtb() const {
        return THDMW_logtb;
    }

    /**
     *
     * @return @f$\tan \beta@f$
     */
    double getTHDMW_tanb() const {
        return THDMW_tanb;
    }

    /**
     *
     * @return @f$\sin \beta@f$
     */
    double getTHDMW_sinb() const {
        return THDMW_sinb;
    }

    /**
     *
     * @return @f$\cos \beta@f$
     */
    double getTHDMW_cosb() const {
        return THDMW_cosb;
    }

    /**
     *
     * @return @f$\beta-\alpha@f$
     */
    double getTHDMW_bma() const {
        return THDMW_bma;
    }

    /**
     *
     * @return @f$\sin(\beta-\alpha)@f$
     */
    double getTHDMW_sin_ba() const {
        return THDMW_sin_ba;
    }

    /**
     *
     * @return @f$\cos \alpha@f$
     */
    double getTHDMW_cosa() const{
        return cos(atan(pow(10.,THDMW_logtb))-THDMW_bma);
    }

    /**
     *
     * @return @f$\sin \alpha@f$
     */
    double getTHDMW_sina() const{
        return sin(atan(pow(10.,THDMW_logtb))-THDMW_bma);
    }

    /**
     *
     * @return @f$\lambda_1@f$
     */
    double getTHDMW_lambda1() const {
        return THDMW_lambda1;
    }

    /**
     *
     * @return @f$\lambda_2@f$
     */
    double getTHDMW_lambda2() const {
        if (flag_modelTHDMW == "ManoharWise" || flag_modelTHDMW == "custodialMW") {
            return 0.0;
        }
        else {
            return THDMW_lambda2;
        }
    }

    /**
     *
     * @return @f$\lambda_3@f$
     */
    double getTHDMW_lambda3() const {
        if (flag_modelTHDMW == "ManoharWise" || flag_modelTHDMW == "custodialMW") {
            return 0.0;
        }
        else {
            return THDMW_lambda3;
        }
    }

    /**
     *
     * @return @f$\lambda_4@f$
     */
    double getTHDMW_lambda4() const {
        if (flag_modelTHDMW == "ManoharWise" || flag_modelTHDMW == "custodialMW") {
            return 0.0;
        }
        else {
            return THDMW_lambda4;
        }
    }

    /**
     *
     * @return @f$\lambda_5@f$
     */
    double getTHDMW_lambda5() const {
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
     *
     * @return @f$m_S^2@f$
     */
    double getTHDMW_mS2() const {
        return THDMW_mS2;
    }

    /**
     *
     * @return @f$\mu_1@f$
     */
    double getTHDMW_mu1() const {
        return THDMW_mu1;
    }

    /**
     *
     * @return @f$\mu_2@f$
     */
    double getTHDMW_mu2() const {
        if (flag_modelTHDMW == "custodial1" || flag_modelTHDMW == "custodialMW") {
            return THDMW_mu1;
        }
        else {
            return THDMW_mu2;
        }
    }

    /**
     *
     * @return @f$\mu_3@f$
     */
    double getTHDMW_mu3() const {
        return THDMW_mu3;
    }
    /**
     *
     * @return @f$\mu_4@f$
     */
    double getTHDMW_mu4() const {
        return THDMW_mu4;
    }

    /**
     *
     * @return @f$\mu_5@f$
     */
    double getTHDMW_mu5() const {
        if (flag_modelTHDMW == "custodial1" || flag_modelTHDMW == "custodialMW") {
            return THDMW_mu4;
        }
        else {
            return THDMW_mu5;
        }
    }

    /**
     *
     * @return @f$\mu_6@f$
     */
    double getTHDMW_mu6() const {
        if (flag_modelTHDMW == "custodial1" || flag_modelTHDMW == "custodialMW") {
            return 2.0*THDMW_mu1;
        }
        else {
            return THDMW_mu6;
        }
    }

    /**
     *
     * @return @f$\nu_1@f$
     */
    double getTHDMW_nu1() const {
        return THDMW_nu1;
    }

    /**
     *
     * @return @f$\nu_2@f$
     */
    double getTHDMW_nu2() const {
        return THDMW_nu2;
    }

    /**
     *
     * @return @f$\nu_3@f$
     */
    double getTHDMW_nu3() const {
        if (flag_modelTHDMW == "custodial1" || flag_modelTHDMW == "custodialMW") {
            return 0.5*THDMW_nu2;
        }
        else {
            return THDMW_nu3;
        }
    }
    /**
     *
     * @return @f$\nu_4@f$
     */
    double getTHDMW_nu4() const {
        return THDMW_nu4;
    }

    /**
     *
     * @return @f$\nu_5@f$
     */
    double getTHDMW_nu5() const {
        if (flag_modelTHDMW == "custodial1" || flag_modelTHDMW == "custodialMW") {
            return THDMW_nu4;
        }
        else {
            return THDMW_nu5;
        }
    }

    /**
     *
     * @return @f$\omega_1@f$
     */
    double getTHDMW_omega1() const {
        if (flag_modelTHDMW == "ManoharWise" || flag_modelTHDMW == "custodialMW") {
            return 0.0;
        }
        else {
            return THDMW_omega1;
        }
    }

    /**
     *
     * @return @f$\omega_2@f$
     */
    double getTHDMW_omega2() const {
        if (flag_modelTHDMW == "ManoharWise" || flag_modelTHDMW == "custodialMW") {
            return 0.0;
        }
        else {
            return THDMW_omega2;
        }
    }

    /**
     *
     * @return @f$\omega_3@f$
     */
    double getTHDMW_omega3() const {
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
     *
     * @return @f$\omega_4@f$
     */
    double getTHDMW_omega4() const {
        if (flag_modelTHDMW == "ManoharWise" || flag_modelTHDMW == "custodialMW") {
            return 0.0;
        }
        else {
            return THDMW_omega4;
        }
    }

    /**
     *
     * @return @f$\omega_5@f$
     */
    double getTHDMW_omega5() const {
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
     *
     * @return @f$\kappa_1@f$
     */
    double getTHDMW_kappa1() const {
        if (flag_modelTHDMW == "ManoharWise" || flag_modelTHDMW == "custodialMW") {
            return 0.0;
        }
        else {
            return THDMW_kappa1;
        }
    }

    /**
     *
     * @return @f$\kappa_2@f$
     */
    double getTHDMW_kappa2() const {
        if (flag_modelTHDMW == "ManoharWise" || flag_modelTHDMW == "custodialMW") {
            return 0.0;
        }
        else {
            return THDMW_kappa2;
        }
    }

    /**
     *
     * @return @f$\kappa_3@f$
     */
    double getTHDMW_kappa3() const {
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
     *
     * @return @f$\eta_U@f$
     */
    double getTHDMW_etaU() const {
        return THDMW_etaU;
    }

    /**
     *
     * @return @f$\eta_D@f$
     */
    double getTHDMW_etaD() const {
        return THDMW_etaD;
    }
    
    
    /**
     *
     * @return @f$\eta_U@f$
     */
    double getTHDMW_rho_b() const {
        return THDMW_rho_b;
    }

    /**
     *
     * @return @f$\eta_D@f$
     */
    double getTHDMW_S_b() const {
        return THDMW_S_b;
    }
    
    

    /**
     *
     * @return THDMW scale
     */
    double getQ_THDMW() const {
        return Q_THDMW;
    }

    /**
     *
     * @return Minimal R' value
     */
    double getRpepsTHDMW() const {
        return RpepsTHDMW;
    }

    /**
     *
     * @return Minimal NLO unitarity check scale
     */
    double getNLOuniscaleTHDMW() const {
        return NLOuniscaleTHDMW;
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

    /*
     * @brief Method overriding the NPbase Mw with the tree-level W mass
     */
    virtual double Mw() const;

protected:

    virtual void setParameter(const std::string, const double&);

    /**
     * @brief A method to check if the model type name in string form is valid.
     * @param[in] THDMW model type name
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
     * @param[in] THDMW RGE order
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
