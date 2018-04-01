/* 
 * Copyright (C) 2017 HEPfit Collaboration
 *
 *
 * For the licensing terms see doc/COPYING.
 */

#ifndef THDMW_H
#define	THDMW_H

#include "StandardModel.h"
//#include "THDMWMatching.h"

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
 *   <td class="mod_name">%m12_2</td>
 *   <td class="mod_symb">@f$m_{12}^2@f$</td>
 *   <td class="mod_desc"></td>
 * </tr>
 * </table>
 */
class THDMW: public StandardModel {
public:

    static const int NTHDMWvars = 30;
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
//    virtual THDMWMatching& getMatching() const
//    {
//        return THDMWM.getObj();
//    }

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
        return THDMW_lambda2;
    }

    /**
     *
     * @return @f$\lambda_3@f$
     */
    double getTHDMW_lambda3() const {
        return THDMW_lambda3;
    }

    /**
     *
     * @return @f$\lambda_4@f$
     */
    double getTHDMW_lambda4() const {
        return THDMW_lambda4;
    }

    /**
     *
     * @return @f$\lambda_5@f$
     */
    double getTHDMW_lambda5() const {
        if (flag_modelTHDMW == "custodial1") {
            return THDMW_lambda4;
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
        if (flag_modelTHDMW == "custodial1") {
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
        if (flag_modelTHDMW == "custodial1") {
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
        if (flag_modelTHDMW == "custodial1") {
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
        if (flag_modelTHDMW == "custodial1") {
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
        if (flag_modelTHDMW == "custodial1") {
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
        return THDMW_omega1;
    }

    /**
     *
     * @return @f$\omega_2@f$
     */
    double getTHDMW_omega2() const {
        return THDMW_omega2;
    }

    /**
     *
     * @return @f$\omega_3@f$
     */
    double getTHDMW_omega3() const {
        if (flag_modelTHDMW == "custodial1") {
            return 0.5*THDMW_omega2;
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
        return THDMW_omega4;
    }

    /**
     *
     * @return @f$\omega_5@f$
     */
    double getTHDMW_omega5() const {
        if (flag_modelTHDMW == "custodial1") {
            return THDMW_omega4;
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
        return THDMW_kappa1;
    }

    /**
     *
     * @return @f$\kappa_2@f$
     */
    double getTHDMW_kappa2() const {
        return THDMW_kappa2;
    }

    /**
     *
     * @return @f$\kappa_3@f$
     */
    double getTHDMW_kappa3() const {
        if (flag_modelTHDMW == "custodial1") {
            return THDMW_kappa2;
        }
        else {
            return THDMW_kappa3;
        }
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
        else if (modeltype.compare("custodial2") == 0)
            return true;
        else
            return false;
    }

//    mutable Matching<THDMWMatching,THDMW> THDMWM; ///< An object of type Matching.

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
    double Q_THDMW, RpepsTHDMW, NLOuniscaleTHDMW;
    std::string flag_modelTHDMW, flag_RGEorder;
};

#endif	/* THDMW_H */
