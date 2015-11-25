/* 
 * Copyright (C) 2012 HEPfit Collaboration
 * All rights reserved.
 *
 * For the licensing terms see doc/COPYING.
 */

#ifndef THDM_H
#define	THDM_H

#include "StandardModel.h"
#include "THDMMatching.h"
#include "THDMcache.h"

//class THDMcache; //forward reference to THDMcache class

/**
 * @class THDM
 * @ingroup THDM
 * @brief A base class @f$Z_2@f$ symmetric Two-Higgs-Doublet models. 
 * @author HEPfit Collaboration
 * @copyright GNU General Public License
 * @details 
 *
 * 
 * @anchor THDMParameters
 * <h3>%Model parameters</h3>
 *
 * The model parameters of %THDM are summarized below.
 * The current implementation only allows for a @f$Z_2@f$ symmetric Two-Higgs-Doublet model without CP violation in the Higgs potential.
 * The scalar 125 GeV resonance is assumed to be the lighter CP-even Higgs state of the model.
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
 *   <td class="mod_name">%mHh2</td>
 *   <td class="mod_symb">@f$m_H^2@f$</td>
 *   <td class="mod_desc">The mass square of the heavier CP-even Higgs state.</td>
 * </tr>
 * <tr>
 *   <td class="mod_name">%mHh2mmA2</td>
 *   <td class="mod_symb">@f$m_H^2-m_A^2@f$</td>
 *   <td class="mod_desc">The difference of the squared masses of the heavier CP-even Higgs and the CP-odd Higgs.</td>
 * </tr>
 * <tr>
 *   <td class="mod_name">%mHh2mmHp2</td>
 *   <td class="mod_symb">@f$m_H^2-m_{H^+}^2@f$</td>
 *   <td class="mod_desc">The difference of the squared masses of the heavier CP-even Higgs and the charged Higgs.</td>
 * </tr>
 * <tr>
 *   <td class="mod_name">%m12_2</td>
 *   <td class="mod_symb">@f$m_{12}^2@f$</td>
 *   <td class="mod_desc">The soft @f$Z_2@f$ symmetry breaking parameter.</td>
 * </tr>
 * <tr>
 *   <td class="mod_name">%bsgamma_theoryerror</td>
 *   <td class="mod_symb">@f$\theta_{b\to s\gamma}^{\text{theo}}@f$</td>
 *   <td class="mod_desc">A nuisance parameter between -1 and +1 for the theoretical error in the determination of BR(B\to X_s \gamma).</td>
 * </tr>
 * </table>
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
 *   <td class="mod_valu">typeI / <b>typeII</b> / typeX / typeY</td>
 *   <td class="mod_desc">This flag determines the type of @f$Z_2@f$ symmetry.</td>
 * </tr>
 * </table>
 *
 */
class THDM: public StandardModel {
public:

    static const int NTHDMvars = 7;
    static const std::string THDMvars[NTHDMvars];
    
    /**
     * @brief THDM constructor
     */
    THDM();
    
    /**
     * @brief THDM destructor
     */
    ~THDM();
    
    virtual bool InitializeModel();
    
    virtual bool Init(const std::map<std::string, double>& DPars);
    
    virtual bool PreUpdate();
    
    virtual bool Update(const std::map<std::string, double>& DPars);
    
    virtual bool PostUpdate();
    
    virtual bool CheckParameters(const std::map<std::string, double>& DPars);

    
    ///////////////////////////////////////////////////////////////////////////
    // Flags

    virtual bool setFlagStr(const std::string name, const std::string value);

    virtual THDMMatching* getMyMatching() const
    {
        return myTHDMMatching;
    }

    
    ///////////////////////////////////////////////////////////////////////////
    
    /**
     * 
     * @return the VEV @f$v_1@f$
     */
    double getv1() const {
        return v() * cosb;
    }

    /**
     *
     * @return the VEV @f$v_2@f$
     */
    double getv2() const {
        return v() * sinb;
    }
    
    ///////////////////////////////////////////////////////////////////////////

    /**
     *
     * @return @f$\log(\tan \beta)@f$
     */
    std::string getModelTypeflag() const {
        return flag_model;
    }

    /**
     *
     * @return @f$\log_{10}(\tan \beta)@f$
     */
    double getlogtb() const {
        return logtb;
    }

    /**
     *
     * @return @f$\tan \beta@f$
     */
    double gettanb() const {
        return tanb;
    }

    /**
     *
     * @return @f$\sin \beta@f$
     */
    double getsinb() const {
        return sinb;
    }

    /**
     *
     * @return @f$\cos \beta@f$
     */
    double getcosb() const {
        return cosb;
    }

    /**
     *
     * @return @f$\beta-\alpha@f$
     */
    double getbma() const {
        return bma;
    }

    /**
     *
     * @return @f$\sin(\beta-\alpha)@f$
     */
    double getsin_ba() const {
        return sin_ba;
    }

    /**
     *
     * @return @f$\cos \alpha@f$
     */
    double getcosa() const{
        return cos(atan(pow(10.,logtb))-bma);
    }

    /**
     *
     * @return @f$\sin \alpha@f$
     */
    double getsina() const{
        return sin(atan(pow(10.,logtb))-bma);
    }

    /**
     *
     * @return mass squared heavy neutral scalar Higgs
     */
    double getmHh2() const {
        return mHh2;
    }

    /**
     *
     * @return mass heavy neutral scalar Higgs
     */
    double getmHh() const {
        return sqrt(mHh2);
    }

    /**
     *
     * @return mass squared pseudoscalar Higgs A
     */
    double getmA2() const {
        return mHh2-mHh2mmA2;
    }

    /**
     *
     * @return mass pseudoscalar Higgs A
     */
    double getmA() const {
        return sqrt(mHh2-mHh2mmA2);
    }

    /**
     *
     * @return charged Higgs mass squared
     */
    double getmHp2() const {
        return mHh2-mHh2mmHp2;
    }

    /**
     *
     * @return charged Higgs mass
     */
    double getmHp() const {
        return sqrt(mHh2-mHh2mmHp2);
    }

    /**
     *
     * @return parameter of the Higgs potential @f$m_{12}^2@f$ 
     */
    double getm12_2() const {
        return m12_2;
    }

    /**
     *
     * @return nuisance parameter for the theoretical error on bsgamma
     */
    double getbsgamma_theoryerror() const {
        return bsgamma_theoryerror;
    }

protected: 

    virtual void setParameter(const std::string, const double&);
    THDMcache * mycache;

    /**
     * @brief A method to check if the model type name in string form is valid.
     * @param[in] THDM model type name
     * @return a boolean that is true if the model type name is valid
     */
    bool checkmodelType(const std::string modeltype) const
    {
        if (modeltype.compare("type1") == 0
                || modeltype.compare("type2") == 0
                || modeltype.compare("typeX") == 0
                || modeltype.compare("typeY") == 0)
            return true;
        else
            return false;
    }

private:

    THDMMatching* myTHDMMatching;

    double logtb, tanb, sinb, cosb, bma, sin_ba, mHh2, mHh2mmA2, mHh2mmHp2, m12_2, bsgamma_theoryerror;
    std::string flag_model;
};

#endif	/* THDM_H */
