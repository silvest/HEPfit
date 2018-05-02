/*
 * Copyright (C) 2014 HEPfit Collaboration
 *
 *
 * For the licensing terms see doc/COPYING.
 */

#ifndef HIGGSKIGEN_H
#define	HIGGSKIGEN_H
#include "NPbase.h"

/**
 * @class HiggsKigen
 * @ingroup HiggsExtensions
 * @brief A model class extending the %StandardModel Higgs sector with seven flavour-universal couplings.
 * @author HEPfit Collaboration
 * @copyright GNU General Public License
 * @details This is a Model class containing parameters and functions associated
 * to rescaling the Higgs decay into vector bosons (@f$\kappa_v@f$), 
 * gluons (@f$\kappa_g@f$), photons (@f$\kappa_ga@f$) and fermions (@f$\kappa_u, \kappa_d, \kappa_l@f$)
 * (as well as the corresponding production mechanisms) with respect to the %StandardModel.
 * The invisible decay width is also parametrized independently by Br@f$(H\to invisible)@f$.
 * This class inherits from the %NPbase class, which defines parameters related to generic
 * extensions of the %StandardModel Higgs sector.
 *
 *
 * @anchor HiggsKigenInitialization
 * <h3>Initialization</h3>
 *
 * After creating an instance of the current class,
 * it is required to call the initialization method InitializeModel(), which
 * is needed by the base class. 
 *
 * The initializations and updates of the model parameters are explained
 * below. 
 *
 * 
 * @anchor HiggsKigenParameters
 * <h3>%Model parameters</h3>
 *
 * The model parameters of %HiggsKigen (in addition to the %StandardModel ones) are summarized below:
 * <table class="model">
 * <tr>
 *   <th>Label</th>
 *   <th>LaTeX symbol</th>
 *   <th>Description</th>
 * </tr>
 * <tr>
 *   <td class="mod_name">%Kw</td>
 *   <td class="mod_symb">@f$\kappa_W@f$</td>
 *   <td class="mod_desc">The factor rescaling all Higgs couplings to the W vector bosons with respect to the SM.</td>
 * </tr>
 * <tr>
 *   <td class="mod_name">%Kz</td>
 *   <td class="mod_symb">@f$\kappa_Z@f$</td>
 *   <td class="mod_desc">The factor rescaling all Higgs couplings to the Z vector bosons with respect to the SM.</td>
 * </tr>
 * <tr>
 *   <td class="mod_name">%Kg</td>
 *   <td class="mod_symb">@f$\kappa_g@f$</td>
 *   <td class="mod_desc">The factor rescaling all Higgs couplings to gluons with respect to the SM.</td>
 * </tr>
 * <tr>
 *   <td class="mod_name">%Kga</td>
 *   <td class="mod_symb">@f$\kappa_ga@f$</td>
 *   <td class="mod_desc">The factor rescaling all Higgs couplings to photons with respect to the SM.</td>
 * </tr>
 * <tr>
 *   <td class="mod_name">%Kzga</td>
 *   <td class="mod_symb">@f$\kappa_{zga}@f$</td>
 *   <td class="mod_desc">The factor rescaling all Higgs couplings to the Z and a photon with respect to the SM.</td>
 * </tr>
 * <tr>
 *   <td class="mod_name">%Ku</td>
 *   <td class="mod_symb">@f$\kappa_u@f$</td>
 *   <td class="mod_desc">The factor rescaling all Higgs couplings to the up quark with respect to the SM.</td>
 * </tr>
 * <tr>
 *   <td class="mod_name">%Kc</td>
 *   <td class="mod_symb">@f$\kappa_c@f$</td>
 *   <td class="mod_desc">The factor rescaling all Higgs couplings to the charm quark with respect to the SM.</td>
 * </tr>
 * <tr>
 *   <td class="mod_name">%Kt</td>
 *   <td class="mod_symb">@f$\kappa_t@f$</td>
 *   <td class="mod_desc">The factor rescaling all Higgs couplings to the top quark with respect to the SM.</td>
 * </tr>
 * <tr>
 *   <td class="mod_name">%Kd</td>
 *   <td class="mod_symb">@f$\kappa_d@f$</td>
 *   <td class="mod_desc">The factor rescaling all Higgs couplings to the down quark with respect to the SM.</td>
 * </tr>
 * <tr>
 *   <td class="mod_name">%Ks</td>
 *   <td class="mod_symb">@f$\kappa_s@f$</td>
 *   <td class="mod_desc">The factor rescaling all Higgs couplings to the strange quark with respect to the SM.</td>
 * </tr>
 * <tr>
 *   <td class="mod_name">%Kb</td>
 *   <td class="mod_symb">@f$\kappa_b@f$</td>
 *   <td class="mod_desc">The factor rescaling all Higgs couplings to the bottom quark with respect to the SM.</td>
 * </tr>
 * <tr>
 *   <td class="mod_name">%Ke</td>
 *   <td class="mod_symb">@f$\kappa_e@f$</td>
 *   <td class="mod_desc">The factor rescaling all Higgs couplings to electrons with respect to the SM.</td>
 * </tr>
 * <tr>
 *   <td class="mod_name">%Kmu</td>
 *   <td class="mod_symb">@f$\kappa_\mu@f$</td>
 *   <td class="mod_desc">The factor rescaling all Higgs couplings to muons with respect to the SM.</td>
 * </tr>
 * <tr>
 *   <td class="mod_name">%Ktau</td>
 *   <td class="mod_symb">@f$\kappa_\tau@f$</td>
 *   <td class="mod_desc">The factor rescaling all Higgs couplings to taus with respect to the SM.</td>
 * </tr>
 * <tr>
 *   <td class="mod_name">%KH</td>
 *   <td class="mod_symb">@f$\kappa_H@f$</td>
 *   <td class="mod_desc">The factor rescaling the total Higgs width respect to the SM.</td>
 * </tr>
 * </table>
 * 
 *
 * @anchor HiggsKigenFlags
 * <h3>%Model flags</h3>
 *
 * The Flags of HiggsKigen are summarized below:
 * <table class="model">
 * <tr>
 *   <th>Label</th>
 *   <th>Value</th>
 *   <th>Description</th>
 * </tr>
 * <tr>
 *   <td class="mod_name">%ExoDec</td>
 *   <td class="mod_valu">TRUE&nbsp;/&nbsp;<b>FALSE</b></td>
 *   <td class="mod_desc">This flag is set to TRUE if exotic Higgs decays are allowed.
 *   The default value is FALSE.</td>
 * </tr>
 * <tr>
 *   <td class="mod_name">%Custodial</td>
 *   <td class="mod_valu">TRUE&nbsp;/&nbsp;<b>FALSE</b></td>
 *   <td class="mod_desc">This flag is set to TRUE if KZ=KW is assumed.
 *   The default value is FALSE.</td>
 * </tr>
 * </table>
 *
 * Please read information about parameter initialization and update in the documentation of the %StandardModel class.
 */
class HiggsKigen : public NPbase {
public:

    static const int NHKvKfgenvars = 15; ///< The number of the model parameters in %HiggsKigen.

    /**
     * @brief A string array containing the labels of the model parameters in %HiggsKigen.
     */
    static const std::string HKvKfgenvars[NHKvKfgenvars];

    /**
     * @brief The default constructor.
     */
    HiggsKigen();

    /**
     * @brief The default destructor.
     */
    virtual ~HiggsKigen()
    {
    };
    
    /**
     * @brief The post-update method for %HiggsKigen.
     * @details This method runs all the procedures that are need to be executed
     * after the model is successfully updated.
     * @return a boolean that is true if the execution is successful
     */
    virtual bool PostUpdate();
    
    /**
     * @brief A method to set a flag of %HiggsKigen.
     * @param[in] name name of a model flag
     * @param[in] value the boolean to be assigned to the flag specified by name
     * @return a boolean that is true if the execution is successful
     */
    virtual bool setFlag(const std::string name, const bool value);

    /**
     * @brief A get method to retrieve the factor rescaling the Higgs coupling
     * to the W vector bosons with respect to the SM @f$K_W@f$.
     * @return @f$K_W@f$
     */
    double getKw() const
    {
        return Kw;
    }

    /**
     * @brief A set method to change the factor rescaling the Higgs coupling
     * to the W vector bosons with respect to the SM @f$K_W@f$.
     * @param[in] @f$K_W@f$ the factor rescaling the Higgs coupling to the W vector bosons.
     */
    void setKw(double Kw)
    {
        this->Kw = Kw;
    }
       
    /**
     * @brief A get method to retrieve the factor rescaling the Higgs coupling
     * to the Z vector bosons with respect to the SM @f$K_Z@f$.
     * @return @f$K_Z@f$
     */
    double getKz() const
    {
        return Kz;
    }

    /**
     * @brief A set method to change the factor rescaling the Higgs coupling
     * to the Z vector bosons with respect to the SM @f$K_Z@f$.
     * @param[in] @f$K_Z@f$ the factor rescaling the Higgs coupling to the Z vector bosons.
     */
    void setKz(double Kz)
    {
        this->Kz = Kz;
    }
      
    /**
     * @brief A get method to retrieve the factor rescaling the Higgs coupling
     * to gluons with respect to the SM @f$K_g@f$.
     * @return @f$K_g@f$
     */
    double getKg() const
    {
        return Kg;
    }

    /**
     * @brief A set method to change the factor rescaling the Higgs coupling
     * to gluons with respect to the SM @f$K_g@f$.
     * @param[in] @f$K_g@f$ the factor rescaling the Higgs coupling to gluons.
     */
    void setKg(double Kg)
    {
        this->Kg = Kg;
    }
      
    /**
     * @brief A get method to retrieve the factor rescaling the Higgs coupling
     * to photons with respect to the SM @f$K_ga@f$.
     * @return @f$K_ga@f$
     */
    double getKga() const
    {
        return Kga;
    }

    /**
     * @brief A set method to change the factor rescaling the Higgs coupling
     * to photons with respect to the SM @f$K_ga@f$.
     * @param[in] @f$K_ga@f$ the factor rescaling the Higgs coupling to photons.
     */
    void setKga(double Kga)
    {
        this->Kga = Kga;
    }

    /**
     * @brief A get method to retrieve the factor rescaling the Higgs coupling
     * to Z and photon with respect to the SM @f$K_{zga}@f$.
     * @return @f$K_{zga}@f$
     */
    double getKzga() const
    {
        return Kzga;
    }

    /**
     * @brief A set method to change the factor rescaling the Higgs coupling
     * to Z and photon with respect to the SM @f$K_{zga}@f$.
     * @param[in] @f$K_{zga}@f$ the factor rescaling the Higgs coupling to photons.
     */
    void setKzga(double Kzga)
    {
        this->Kzga = Kzga;
    }

    /**
     * @brief A get method to retrieve the factor rescaling the Higgs coupling
     * to d quarks with respect to the SM @f$K_d@f$.
     * @return @f$K_d@f$
     */
    double getKd() const
    {
        return Kd;
    }

    /**
     * @brief A set method to change the factor rescaling the Higgs coupling
     * to d quarks with respect to the SM @f$K_d@f$.
     * @param[in] @f$K_d@f$ the factor rescaling the Higgs coupling to d quarks.
     */
    void setKd(double Kd)
    {
        this->Kd = Kd;
    }
        
    /**
     * @brief A get method to retrieve the factor rescaling the Higgs coupling
     * to s quarks with respect to the SM @f$K_s@f$.
     * @return @f$K_s@f$
     */
    double getKs() const
    {
        return Ks;
    }

    /**
     * @brief A set method to change the factor rescaling the Higgs coupling
     * to s quarks with respect to the SM @f$K_s@f$.
     * @param[in] @f$K_s@f$ the factor rescaling the Higgs coupling to s quarks.
     */
    void setKs(double Ks)
    {
        this->Ks = Ks;
    }    
    
    /**
     * @brief A get method to retrieve the factor rescaling the Higgs coupling
     * to b quarks with respect to the SM @f$K_b@f$.
     * @return @f$K_b@f$
     */
    double getKb() const
    {
        return Kb;
    }

    /**
     * @brief A set method to change the factor rescaling the Higgs coupling
     * to b quarks with respect to the SM @f$K_b@f$.
     * @param[in] @f$K_b@f$ the factor rescaling the Higgs coupling to b quarks.
     */
    void setKb(double Kb)
    {
        this->Kb = Kb;
    }

    /**
     * @brief A get method to retrieve the factor rescaling the Higgs coupling
     * to electrons with respect to the SM @f$K_e@f$.
     * @return @f$K_e@f$
     */
    double getKe() const
    {
        return Ke;
    }

    /**
     * @brief A set method to change the factor rescaling the Higgs coupling
     * to electrons with respect to the SM @f$K_e@f$.
     * @param[in] @f$K_e@f$ the factor rescaling the Higgs coupling to electrons.
     */
    void setKe(double Ke)
    {
        this->Ke = Ke;
    }
    
    /**
     * @brief A get method to retrieve the factor rescaling the Higgs coupling
     * to muons with respect to the SM @f$K_\mu@f$.
     * @return @f$K_\mu@f$
     */
    double getKmu() const
    {
        return Kmu;
    }

    /**
     * @brief A set method to change the factor rescaling the Higgs coupling
     * to muons with respect to the SM @f$K_\mu@f$.
     * @param[in] @f$K_\mu@f$ the factor rescaling the Higgs coupling to muons.
     */
    void setKmu(double Kmu)
    {
        this->Kmu = Kmu;
    }
    
    /**
     * @brief A get method to retrieve the factor rescaling the Higgs coupling
     * to taus with respect to the SM @f$K_\tau@f$.
     * @return @f$K_\tau@f$
     */
    double getKtau() const
    {
        return Ktau;
    }

    /**
     * @brief A set method to change the factor rescaling the Higgs coupling
     * to taus with respect to the SM @f$K_\tau@f$.
     * @param[in] @f$K_\tau@f$ the factor rescaling the Higgs coupling to taus.
     */
    void setKtau(double Ktau)
    {
        this->Ktau = Ktau;
    }

    /**
     * @brief A get method to retrieve the factor rescaling the Higgs coupling
     * to u quarks with respect to the SM @f$K_u@f$.
     * @return @f$K_u@f$
     */
    double getKu() const
    {
        return Ku;
    }

    /**
     * @brief A set method to change the factor rescaling the Higgs coupling
     * to u quarks with respect to the SM @f$K_u@f$.
     * @param[in] @f$K_u@f$ the factor rescaling the Higgs coupling to u quarks.
     */
    void setKu(double Ku)
    {
        this->Ku = Ku;
    }
    
    /**
     * @brief A get method to retrieve the factor rescaling the Higgs coupling
     * to c quarks with respect to the SM @f$K_c@f$.
     * @return @f$K_c@f$
     */
    double getKc() const
    {
        return Kc;
    }

    /**
     * @brief A set method to change the factor rescaling the Higgs coupling
     * to c quarks with respect to the SM @f$K_c@f$.
     * @param[in] @f$K_c@f$ the factor rescaling the Higgs coupling to c quarks.
     */
    void setKc(double Kc)
    {
        this->Kc = Kc;
    }    
    
    /**
     * @brief A get method to retrieve the factor rescaling the Higgs coupling
     * to t quarks with respect to the SM @f$K_t@f$.
     * @return @f$K_t@f$
     */
    double getKt() const
    {
        return Kt;
    }

    /**
     * @brief A set method to change the factor rescaling the Higgs coupling
     * to t quarks with respect to the SM @f$K_t@f$.
     * @param[in] @f$K_t@f$ the factor rescaling the Higgs coupling to t quarks.
     */
    void setKt(double Kt)
    {
        this->Kt = Kt;
    }   
    
    
    /**
     * @brief A get method to retrieve the factor rescaling the Higgs width
     * respect to the SM @f$K_H@f$.
     * @return @f$K_H@f$
     */
    double getKH() const
    {
        return KH;
    }

    /**
     * @brief A set method to change the factor rescaling the Higgs width
     * with respect to the SM @f$K_H@f$.
     * @param[in] @f$K_H@f$ the factor rescaling the Higgs width.
     */
    void setKH(double KH)
    {
        this->KH = KH;
    } 

    /**
     * @brief A method to check if all the mandatory parameters for %HiggsKigen
     * have been provided in model initialization.
     * @param[in] DPars a map of the parameters that are being updated in the Monte Carlo run
     * (including parameters that are varied and those that are held constant)
     * @return a boolean that is true if the execution is successful
     */
    virtual bool CheckParameters(const std::map<std::string, double>& DPars);

    ////////////////////////////////////////////////////////////////////////

    /**
     * @brief The oblique parameter @f$S@f$.
     * @return @f$S@f$
     */
    virtual double obliqueS() const;

    /**
     * @brief The oblique parameter @f$T@f$.
     * @return @f$T@f$
     */
    virtual double obliqueT() const;

    /**
     * @brief The oblique parameter @f$U@f$.
     * @return @f$U=0@f$
     */
    virtual double obliqueU() const;

    /**
     * @brief The ratio @f$\mu_{ggH}@f$ between the gluon-gluon fusion Higgs
     * production cross-section in the current model and in the Standard Model.
     * @param[in] sqrt_s the center-of-mass energy in TeV
     * @return @f$\mu_{ggH}@f$
     */
    virtual double muggH(const double sqrt_s) const;
    /**
     * @brief The ratio @f$\mu_{VBF}@f$ between the vector-boson fusion Higgs
     * production cross-section in the current model and in the Standard Model. 
     * @param[in] sqrt_s the center-of-mass energy in TeV
     * @return @f$\mu_{VBF}@f$
     */
    virtual double muVBF(const double sqrt_s) const;
    /**
     * @brief The ratio @f$\mu_{eeWBF}@f$ between the 
     * @f$ e^{+}e^{-}\to \nu\bar{\nu} H @f$ production
     * cross-section in the current model and in the Standard Model.
     * @param[in] sqrt_s the center-of-mass energy in TeV
     * @return @f$\mu_{eeWBF}@f$
     */
    virtual double mueeWBF(const double sqrt_s) const;
    /**
     * @brief The ratio @f$\mu_{epWBF}@f$ between the 
     * @f$ e^{-} p\to \nu j H @f$ production
     * cross-section in the current model and in the Standard Model.
     * @param[in] sqrt_s the center-of-mass energy in TeV
     * @return @f$\mu_{epWBF}@f$
     */
    virtual double muepWBF(const double sqrt_s) const;
    /**
     * @brief The ratio @f$\mu_{epZBF}@f$ between the 
     * @f$ e^{-} p\to e^{-} j H @f$ production
     * cross-section in the current model and in the Standard Model.
     * @param[in] sqrt_s the center-of-mass energy in TeV
     * @return @f$\mu_{epZBF}@f$
     */
    virtual double muepZBF(const double sqrt_s) const;
    /**
     * @brief The ratio @f$\mu_{WH}@f$ between the W-Higgs associated production
     * cross-section in the current model and in the Standard Model. 
     * @param[in] sqrt_s the center-of-mass energy in TeV
     * @return @f$\mu_{WH}@f$
     */
    virtual double muWH(const double sqrt_s) const;
    /**
     * @brief The ratio @f$\mu_{ZH}@f$ between the Z-Higgs associated production
     * cross-section in the current model and in the Standard Model.
     * @param[in] sqrt_s the center-of-mass energy in TeV
     * @return @f$\mu_{ZH}@f$
     */
    virtual double muZH(const double sqrt_s) const;
    /**
     * @brief The ratio @f$\mu_{eeZH}@f$ between the 
     * @f$e^{+}e^{-}\to ZH@f$ associated production
     * cross-section in the current model and in the Standard Model.
     * @param[in] sqrt_s the center-of-mass energy in TeV
     * @return @f$\mu_{eeZH}@f$
     */
    virtual double mueeZH(const double sqrt_s) const;
    /**
     * @brief The ratio @f$\mu_{VH}@f$ between the WH+ZH associated production
     * cross-section in the current model and in the Standard Model.
     * @param[in] sqrt_s the center-of-mass energy in TeV
     * @return @f$\mu_{VH}@f$
     */
    virtual double muVH(const double sqrt_s) const;
    /**
     * @brief The ratio @f$\mu_{VBF+VH}@f$ between the sum of VBF and WH+ZH associated production
     * cross-section in the current model and in the Standard Model.
     * @param[in] sqrt_s the center-of-mass energy in TeV
     * @return @f$\mu_{VBF+VH}@f$
     */
    virtual double muVBFpVH(const double sqrt_s) const;
    /**
     * @brief The ratio @f$\mu_{ttH}@f$ between the t-tbar-Higgs associated 
     * production cross-section in the current model and in the Standard Model.
     * @param[in] sqrt_s the center-of-mass energy in TeV
     * @return @f$\mu_{ttH}@f$
     */
    virtual double muttH(const double sqrt_s) const;
    /**
     * @brief The ratio @f$\mu_{ggH+ttH}@f$ between the sum of gluon-gluon fusion
     * and t-tbar-Higgs associated 
     * production cross-section in the current model and in the Standard Model.
     * @param[in] sqrt_s the center-of-mass energy in TeV
     * @return @f$\mu_{ggH+ttH}@f$
     */
    virtual double muggHpttH(const double sqrt_s) const;
    /**
     * @brief The ratio @f$\mu_{eettH}@f$ between the 
     * @f$ e^{+}e^{-}\to t\bar{t} H @f$ production
     * cross-section in the current model and in the Standard Model.
     * @param[in] sqrt_s the center-of-mass energy in TeV
     * @return @f$\mu_{eettH}@f$
     */
    virtual double mueettH(const double sqrt_s) const;
    /**
     * @brief The ratio of the Br@f$(H\to gg)@f$ in the current model
     * and in the Standard Model.
     * @return Br@f$(H\to gg)@f$/Br@f$(H\to gg)_{\mathrm{SM}}@f$
     */
    virtual double BrHggRatio() const;
    /**
     * @brief The ratio of the Br@f$(H\to WW)@f$ in the current model
     * and in the Standard Model.
     * @return Br@f$(H\to WW)@f$/Br@f$(H\to WW)_{\mathrm{SM}}@f$
     */
    virtual double BrHWWRatio() const;
    /**
     * @brief The ratio of the Br@f$(H\to ZZ)@f$ in the current model
     * and in the Standard Model.
     * @return Br@f$(H\to ZZ)@f$/Br@f$(H\to ZZ)_{\mathrm{SM}}@f$
     */
    virtual double BrHZZRatio() const;
    /**
     * @brief The ratio of the Br@f$(H\to Z\gamma)@f$ in the current model
     * and in the Standard Model.
     * @return Br@f$(H\to Z\gamma)@f$/Br@f$(H\to Z\gamma)_{\mathrm{SM}}@f$
     */
    virtual double BrHZgaRatio() const;
    /**
     * @brief The ratio of the Br@f$(H\to \gamma\gamma)@f$ in the current model
     * and in the Standard Model.
     * @return Br@f$(H\to \gamma\gamma)@f$/Br@f$(H\to \gamma\gamma)_{\mathrm{SM}}@f$
     */
    virtual double BrHgagaRatio() const;
    /**
     * @brief The ratio of the Br@f$(H\to \mu^+\mu^-)@f$ in the current model
     * and in the Standard Model.
     * @return Br@f$(H\to \mu^+\mu^-)@f$/Br@f$(H\to \mu^+\mu^-)_{\mathrm{SM}}@f$
     */
    virtual double BrHmumuRatio() const;
    /**
     * @brief The ratio of the Br@f$(H\to \tau^+\tau^-)@f$ in the current model
     * and in the Standard Model.
     * @return Br@f$(H\to \tau^+\tau^-)@f$/Br@f$(H\to \tau^+\tau^-)_{\mathrm{SM}}@f$
     */
    virtual double BrHtautauRatio() const;
    /**
     * @brief The ratio of the Br@f$(H\to c\bar{c})@f$ in the current model
     * and in the Standard Model.
     * @return Br@f$(H\to c\bar{c})@f$/Br@f$(H\to c\bar{c})_{\mathrm{SM}}@f$
     */
    virtual double BrHccRatio() const;
    /**
     * @brief The ratio of the Br@f$(H\to b\bar{b})@f$ in the current model
     * and in the Standard Model.
     * @return Br@f$(H\to b\bar{b})@f$/Br@f$(H\to b\bar{b})_{\mathrm{SM}}@f$
     */
    virtual double BrHbbRatio() const;
    /**
     * @brief The ratio @f$\sigma(ttH)/\sigma(ttZ)@f$ 
     * in the @f$H,Z\to b\bar{b}@f$ channel in the current model and in the Standard Model.
     * @param[in] sqrt_s the center-of-mass energy in TeV
     * @return @f$\sigma(ttH)/\sigma(ttZ)@f$ normalized to the SM
     */
    virtual double muttHZbbboost(const double sqrt_s) const;

    virtual double muttHbb(const double sqrt_s) const;
    virtual double muttHWW(const double sqrt_s) const;
    virtual double muttHtautau(const double sqrt_s) const;
    /**
     * @brief The ratio of the @f$\Gamma(H)@f$ in the current model
     * and in the Standard Model.
     * @return @f$\Gamma(H)@f$/@f$\Gamma(H)_{\mathrm{SM}}@f$
     */
    virtual double computeGammaTotalRatio() const;
    
    /**
     * @brief The branching ratio of the of the Higgs into exotic particles (invisible or not).
     * @return Br@f$(H\to exotic)@f$
     */
    virtual double BrHexo() const;

    ////////////////////////////////////////////////////////////////////////
protected:

    /**
     * @brief A method to set the value of a parameter of %HiggsKigen.
     * @param[in] name name of a model parameter
     * @param[in] value the value to be assigned to the parameter specified by name
     */
    virtual void setParameter(const std::string name, const double& value);

    /**
     * @brief A method to compute the ratio of the @f$Hgg@f$ coupling in the current model and in the SM.
     * @return the ratio of the @f$Hgg@f$ coupling in the current model and in the SM
     */
    virtual double computeKg() const;

    /**
     * @brief A method to compute the ratio of the @f$HWW@f$ coupling in the current model and in the SM.
     * @return the ratio of the @f$HWW@f$ coupling in the current model and in the SM
     */
    virtual double computeKW() const;

    /**
     * @brief A method to compute the ratio of the @f$HZZ@f$ coupling in the current model and in the SM.
     * @return the ratio of the @f$HZZ@f$ coupling in the current model and in the SM
     */
    virtual double computeKZ() const;

    /**
     * @brief A method to compute the ratio of the @f$HZ\gamma@f$ coupling in the current model and in the SM.
     * @return the ratio of the @f$HZ\gamma@f$ coupling in the current model and in the SM
     */
    virtual double computeKZga() const;

    /**
     * @brief A method to compute the ratio of the @f$H\gamma\gamma@f$ coupling in the current model and in the SM.
     * @return the ratio of the @f$H\gamma\gamma@f$ coupling in the current model and in the SM
     */
    virtual double computeKgaga() const;
    
    /**
     * @brief A method to compute the ratio of the @f$H\mu\mu@f$ coupling in the current model and in the SM.
     * @return the ratio of the @f$H\mu\mu@f$ coupling in the current model and in the SM
     */
    virtual double computeKmu() const;

    /**
     * @brief A method to compute the ratio of the @f$H\tau\tau@f$ coupling in the current model and in the SM.
     * @return the ratio of the @f$H\tau\tau@f$ coupling in the current model and in the SM
     */
    virtual double computeKtau() const;

    /**
     * @brief A method to compute the ratio of the @f$Hcc@f$ coupling in the current model and in the SM.
     * @return the ratio of the @f$Hcc@f$ coupling in the current model and in the SM
     */
    virtual double computeKc() const;

    /**
     * @brief A method to compute the ratio of the @f$Htt@f$ coupling in the current model and in the SM.
     * @return the ratio of the @f$Htt@f$ coupling in the current model and in the SM
     */
    virtual double computeKt() const;

    /**
     * @brief A method to compute the ratio of the @f$Hbb@f$ coupling in the current model and in the SM.
     * @return the ratio of the @f$Hbb@f$ coupling in the current model and in the SM
     */
    virtual double computeKb() const;

    ////////////////////////////////////////////////////////////////////////
private:
    double Kw; ///< The factor rescaling all Higgs couplings to the W vector bosons with respect to the SM.
    double Kz; ///< The factor rescaling all Higgs couplings to the Z vector bosons with respect to the SM.
    double Kg; ///< The factor rescaling all Higgs couplings to gluons with respect to the SM.
    double Kga; ///< The factor rescaling all Higgs couplings to photons with respect to the SM.
    double Kzga; ///< The factor rescaling all Higgs couplings to the Z and a photon with respect to the SM.
    double Ku; ///< The factor rescaling all Higgs couplings to the up quark with respect to the SM.
    double Kc; ///< The factor rescaling all Higgs couplings to the charm quark with respect to the SM.
    double Kt; ///< The factor rescaling all Higgs couplings to the top quark with respect to the SM.
    double Kd; ///< The factor rescaling all Higgs couplings to the down quark with respect to the SM.
    double Ks; ///< The factor rescaling all Higgs couplings to the strange quark with respect to the SM.
    double Kb; ///< The factor rescaling all Higgs couplings to the bottom quark with respect to the SM.
    double Ke; ///< The factor rescaling all Higgs couplings to electrons with respect to the SM.
    double Kmu; ///< The factor rescaling all Higgs couplings to muons with respect to the SM.
    double Ktau; ///< The factor rescaling all Higgs couplings to taus with respect to the SM.
    double KH; ///< The factor rescaling the total Higgs width respect to the SM.
    
    bool FlagExoDec; ///< A boolean flag that is true if one allows exotic Higgs decays (invisible or not).
    bool FlagCustodial; ///< A boolean flag that is true if KZ=KW.

};

#endif	/* HIGGSKIGEN_H */

