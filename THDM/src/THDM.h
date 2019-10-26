/* 
 * Copyright (C) 2012 HEPfit Collaboration
 *
 *
 * For the licensing terms see doc/COPYING.
 */

#ifndef THDM_H
#define	THDM_H

#include "StandardModel.h"
#include "THDMMatching.h"

class THDMcache; //forward reference to THDMcache class

/**
 * @class THDM
 * @ingroup THDM
 * @brief A base class for @f$Z_2@f$ symmetric Two-Higgs-Doublet models. 
 * @author HEPfit Collaboration
 * @copyright GNU General Public License
 * @details 
 *
 * 
 * @anchor THDMParameters
 * <h3>%Model parameters</h3>
 *
 * The model parameters of %THDM are summarized below.
 * The current implementation only allows for a Two-Higgs-Doublet model with a softly broken @f$Z_2@f$ symmetry and without CP violation in the Higgs potential.
 * The scalar 125 GeV resonance is assumed to be the one of the CP-even Higgs states of the model; the other is attributed to mHh in the THDM configuration.
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
 *   <td class="mod_desc">The mass square of the "non-125 GeV" CP-even Higgs state.</td>
 * </tr>
 * <tr>
 *   <td class="mod_name">%mA2</td>
 *   <td class="mod_symb">@f$m_A^2@f$</td>
 *   <td class="mod_desc">The squared masses of the CP-odd Higgs @f$A@f$.</td>
 * </tr>
 * <tr>
 *   <td class="mod_name">%mHp2</td>
 *   <td class="mod_symb">@f$m_{H^+}^2@f$</td>
 *   <td class="mod_desc">The squared masses of the charged Higgs.</td>
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
 * <tr>
 *   <td class="mod_name">%Q_THDM</td>
 *   <td class="mod_symb">@f$Q_{\text{THDM}}@f$</td>
 *   <td class="mod_desc">The THDM scale.</td>
 * </tr>
 * <tr>
 *   <td class="mod_name">%Rpeps</td>
 *   <td class="mod_symb">@f$\varepsilon_{\text{R'}}@f$</td>
 *   <td class="mod_desc">Minimal value for which the R' perturbativity criterion should be applied for the unitarity bounds.</td>
 * </tr>
 * <tr>
 *   <td class="mod_name">%NLOuniscale</td>
 *   <td class="mod_symb">@f$Q_{\text{min}}^{\text{NLOuni}}@f$</td>
 *   <td class="mod_desc">Minimal scale at which the NLO unitarity conditions are checked.</td>
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
 *   <td class="mod_valu">typeI / typeII / typeX / typeY</td>
 *   <td class="mod_desc">This flag determines the type of @f$Z_2@f$ symmetry.</td>
 * </tr>
 * <tr>
 *   <td class="mod_name">%RGEorder</td>
 *   <td class="mod_valu">LO / approxNLO / NLO</td>
 *   <td class="mod_desc">This flag determines the order in perturbation theory of the renormalization group equations.</td>
 * </tr>
 * <tr>
 *   <td class="mod_name">%wavefunctionrenormalization</td>
 *   <td class="mod_valu">true / false</td>
 *   <td class="mod_desc">Whether to use wavefunction renormalization for NLO unitarity constraints.</td>
 * </tr>
 * </table>
 *
 */
class THDM: public StandardModel {
public:

    static const int NTHDMvars = 16;
    static std::string THDMvars[NTHDMvars];
    
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

    /**
     * @brief A get method to access the member reference of type StandardModelMatching.
     * @return a reference to a StandardModelMatching object
     */
    virtual THDMMatching& getMatching() const
    {
        return THDMM.getObj();
    }

    ///////////////////////////////////////////////////////////////////////////
    // Flags

    virtual bool setFlagStr(const std::string name, const std::string value);
    virtual bool setFlag(const std::string, const bool);

    THDMcache* getMyTHDMCache() const
    {
        return myTHDMcache;
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
     * @return THDM model type
     */
    std::string getModelTypeflag() const {
        return flag_model;
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
     * @return Choose if you want to use the THDM masses or rather their squares
     */
    bool getsqmassesflag() const {
        return flag_use_sq_masses;
    }

    /**
     *
     * @return Flag to switch on wavefunction renormalization for the NLO unitarity conditions
     */
    bool getWFRflag() const {
        return flag_wfr;
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
     * @return squared mass of the lighter neutral scalar Higgs
     */
    double getmHl2() const {
        if(flag_use_sq_masses) {
            if(mHh2 < mHl2) {
                return mHh2;
            }
            else
            {
                return mHl2;
            }
        }
        else
        {
            if(mHh1*mHh1 < mHl2) {
                return mHh1*mHh1;
            }
            else
            {
                return mHl2;
            }
        }
    }

    /**
     *
     * @return squared mass of the "non-125 GeV" neutral scalar Higgs
     */
    double getmHh2() const {
        if(flag_use_sq_masses) {
            if(mHh2 < 0.) {
                return 0.;
            }
            else if(mHh2 < mHl2) {
                return mHl2;
            }
            else
            {
                return mHh2;
            }
        }
        else
        {
            if(mHh1*mHh1 < mHl2) {
                return mHl2;
            }
            else
            {
                return mHh1*mHh1;
            }
        }
    }

    /**
     *
     * @return mass of the "non-125 GeV" neutral scalar Higgs
     */
    double getmHh() const {
        if(flag_use_sq_masses) {
            if(mHh2 < 0.) {
                return 0.;
            }
            else if(mHh2 < mHl2) {
                return sqrt(mHl2);
            }
            else
            {
                return sqrt(mHh2);
            }
        }
        else
        {
            if(mHh1*mHh1 < mHl2) {
                return sqrt(mHl2);
            }
            else
            {
                return mHh1;
            }
        }
    }

    /**
     *
     * @return squared mass of the pseudoscalar Higgs A
     */
    double getmA2() const {
        if(flag_use_sq_masses) {
            return mA2;
        }
        else
        {
            return mA1*mA1;
        }
    }

    /**
     *
     * @return mass of the pseudoscalar Higgs A
     */
    double getmA() const {
        if(flag_use_sq_masses) {
            if(mA2 < 0.) {
                return 0.;
            }
            else
            {
                return sqrt(mA2);
            }
        }
        else
        {
                return mA1;
        }
    }

    /**
     *
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
     * @return charged Higgs mass
     */
    double getmHp() const {
        if(flag_use_sq_masses) {
            if(mHp2 < 0.) {
                return 0.;
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
     *
     * @return parameter of the Higgs potential @f$m_{12}^2@f$ 
     */
    double getm12_2() const {
        return m12_2;
    }

    /**
     *
     * @return BDtaunu SM expectation
     */
    double getBDtaunu_SM() const {
        return BDtaunu_SM;
    }

    /**
     *
     * @return BDtaunu coefficient A
     */
    double getBDtaunu_A() const {
        return BDtaunu_A;
    }

    /**
     *
     * @return BDtaunu coefficient B
     */
    double getBDtaunu_B() const {
        return BDtaunu_B;
    }

    /**
     *
     * @return BDstartaunu SM expectation
     */
    double getBDstartaunu_SM() const {
        return BDstartaunu_SM;
    }

    /**
     *
     * @return BDstartaunu coefficient A
     */
    double getBDstartaunu_A() const {
        return BDstartaunu_A;
    }

    /**
     *
     * @return BDtaunu coefficient B
     */
    double getBDstartaunu_B() const {
        return BDstartaunu_B;
    }

    /**
     *
     * @return nuisance parameter for the theoretical error on bsgamma
     */
    double getbsgamma_theoryerror() const {
        return bsgamma_theoryerror;
    }

    /**
     *
     * @return THDM scale
     */
    double getQ_THDM() const {
        return Q_THDM;
    }

    /**
     *
     * @return Minimal R' value
     */
    double getRpeps() const {
        return Rpeps;
    }

    /**
     *
     * @return Minimal NLO unitarity check scale
     */
    double getNLOuniscale() const {
        return NLOuniscale;
    }

protected: 

    virtual void setParameter(const std::string, const double&);
//    THDMcache * mycache;

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

    /**
     * @brief A method to check if the RGE order name in string form is valid.
     * @param[in] THDM RGE order
     * @return a boolean that is true if the RGE order string is valid
     */
    bool checkRGEorder(const std::string RGEorder) const
    {
        if (RGEorder.compare("LO") == 0
                || RGEorder.compare("approxNLO") == 0
                || RGEorder.compare("NLO") == 0)
            return true;
        else
            return false;
    }

    mutable Matching<THDMMatching,THDM> THDMM; ///< An object of type Matching.

private:

    THDMcache* myTHDMcache;

    double logtb, tanb, sinb, cosb, bma, sin_ba, mHh1, mA1, mHp1, mHh2, mA2, mHp2, m12_2, bsgamma_theoryerror, Q_THDM, Rpeps, NLOuniscale;
    double mHl2;
    double BDtaunu_SM, BDtaunu_A, BDtaunu_B, BDstartaunu_SM, BDstartaunu_A, BDstartaunu_B;
    std::string flag_model, flag_RGEorder;
    bool flag_use_sq_masses, flag_wfr;
};

#endif	/* THDM_H */
