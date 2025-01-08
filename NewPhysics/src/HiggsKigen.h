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
 * @ingroup NewPhysics
 * @brief A model class extending the %StandardModel Higgs sector with seven flavour-universal couplings.
 * @author HEPfit Collaboration
 * @copyright GNU General Public License
 * @details This is a Model class containing parameters and functions associated
 * to rescaling the Higgs decay into vector bosons (@f$\kappa_W,\kappa_Z@f$), 
 * gluons (@f$\kappa_g@f$), photons (@f$\kappa_\gamma@f$), Z and photons (@f$\kappa_{Z\gamma}@f$) 
 * and fermions (@f$\kappa_{u,c,t}, \kappa_{d,s,b}, \kappa_{e\mu\tau}@f$)
 * (as well as the corresponding production mechanisms) with respect to the %StandardModel.
 * The possibility of extra decay width is also parametrized independently by Br@f$(H\to invisible)@f$ and Br@f$(H\to exotic)@f$.
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
 *   <td class="mod_desc">The factor rescaling the Higgs total width with respect to the SM. (Only applies when flag UseKH is true. )</td>
 * </tr>
 * <tr>
 *   <td class="mod_name">%BrHinv</td>
 *   <td class="mod_symb">Br@f$(H\to invisible)@f$</td>
 *   <td class="mod_desc">The branching ratio of invisible Higgs decays. Only the absolute value of this parameter is considered. </td>
 * </tr>
 * <tr>
 *   <td class="mod_name">%BrHexo</td>
 *   <td class="mod_symb">Br@f$(H\to exotic (visible))@f$</td>
 *   <td class="mod_desc">The branching ratio of exotic, not invisible, Higgs decays. Only the absolute value of this parameter is considered. 
 *   (If flag UseKH is true this parameter is ignored and it is recomputed as a function of KH. )</td>
 * </tr>
 * <tr>
 *   <td class="mod_name">%eXint </td>
 *   <td class="mod_symb">\f$\varepsilon_{X}^{int}\f$</td>
 *   <td class="mod_desc">The relative intrinsic theoretical uncertainty for the process X.
 *    (Only for Higgs observables and assumed to be constant in the energy.)</td>
 * </tr>
 * <tr>
 *   <td class="mod_name">%eXpar </td>
 *   <td class="mod_symb">\f$\varepsilon_{X}^{par}\f$</td>
 *   <td class="mod_desc">The relative parametric theoretical uncertainty for the process X.
 *    (Only for Higgs observables and assumed to be constant in the energy.)</td>
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
 *   <td class="mod_name">%KiLoop</td>
 *   <td class="mod_valu">TRUE&nbsp;/&nbsp;<b>FALSE</b></td>
 *   <td class="mod_desc">This flag is set to TRUE if independent kappa's are used
 *   for the loop induced processes (g,ga,Zga).
 *   The default value is TRUE.</td>
 * </tr>
 * <tr>
 *   <td class="mod_name">%Custodial</td>
 *   <td class="mod_valu">TRUE&nbsp;/&nbsp;<b>FALSE</b></td>
 *   <td class="mod_desc">This flag is set to TRUE if KZ=KW is assumed.
 *   The default value is FALSE.</td>
 * </tr>
 * <tr>
 *   <td class="mod_name">%UniversalKf</td>
 *   <td class="mod_valu">TRUE&nbsp;/&nbsp;<b>FALSE</b></td>
 *   <td class="mod_desc">This flag is set to TRUE if all Kf are assumed to take the same 
 *   universal value. The value is controlled by the parameter Kt.
 *   The default value is FALSE.</td>
 * </tr>
 * <tr>
 *   <td class="mod_name">%UniversalK</td>
 *   <td class="mod_valu">TRUE&nbsp;/&nbsp;<b>FALSE</b></td>
 *   <td class="mod_desc">This flag is set to TRUE if all K are assumed to take the same 
 *   universal value. The value is controlled by the parameter Kt. (It does not apply to the extra BR.)
 *   The default value is FALSE.</td>
 * </tr>
 * <tr>
 *   <td class="mod_name">%UseKH</td>
 *   <td class="mod_valu">TRUE&nbsp;/&nbsp;<b>FALSE</b></td>
 *   <td class="mod_desc">This flag is set to TRUE if the total Higgs width is rescaled by KH. Otherwise, it is computed
 *   from the other SM decays, BrHinv and BrHexo. (If true BrHexo is ignored and computed from KH.)
 *   The default value is FALSE.</td>
 * </tr>
 * </table>
 *
 * Please read information about parameter initialization and update in the documentation of the %StandardModel class.
 */
class HiggsKigen : public NPbase {
public:

    static const int NHKvKfgenvars = 94; ///< The number of the model parameters in %HiggsKigen: 17 modelpars + 39 + 38 th. err. pars.

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
     * @param[in] Kw @f$K_W@f$ the factor rescaling the Higgs coupling to the W vector bosons.
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
     * @param[in] Kz @f$K_Z@f$ the factor rescaling the Higgs coupling to the Z vector bosons.
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
     * @param[in] Kg @f$K_g@f$ the factor rescaling the Higgs coupling to gluons.
     */
    void setKg(double Kg)
    {
        this->Kg = Kg;
    }
      
    /**
     * @brief A get method to retrieve the factor rescaling the Higgs coupling
     * to photons with respect to the SM @f$K_ga@f$.
     * @return Kga @f$K_ga@f$ the factor rescaling the Higgs coupling to photons.
     */
    double getKga() const
    {
        return Kga;
    }

    /**
     * @brief A set method to change the factor rescaling the Higgs coupling
     * to photons with respect to the SM @f$K_ga@f$.
     * @param[in] Kga @f$K_ga@f$ the factor rescaling the Higgs coupling to photons.
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
     * @param[in] Kzga @f$K_{zga}@f$ the factor rescaling the Higgs coupling to photons.
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
     * @param[in] Kd @f$K_d@f$ the factor rescaling the Higgs coupling to d quarks.
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
     * @param[in] Ks @f$K_s@f$ the factor rescaling the Higgs coupling to s quarks.
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
     * @param[in] Kb @f$K_b@f$ the factor rescaling the Higgs coupling to b quarks.
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
     * @param[in] Ke @f$K_e@f$ the factor rescaling the Higgs coupling to electrons.
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
     * @param[in] Kmu @f$K_\mu@f$ the factor rescaling the Higgs coupling to muons.
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
     * @param[in] Ktau @f$K_\tau@f$ the factor rescaling the Higgs coupling to taus.
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
     * @param[in] Ku @f$K_u@f$ the factor rescaling the Higgs coupling to u quarks.
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
     * @param[in] Kc @f$K_c@f$ the factor rescaling the Higgs coupling to c quarks.
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
     * @param[in] Kt @f$K_t@f$ the factor rescaling the Higgs coupling to t quarks.
     */
    void setKt(double Kt)
    {
        this->Kt = Kt;
    }   
    
    /**
     * @brief A get method to retrieve the factor rescaling the total Higgs width
     * with respect to the SM @f$K_H@f$.
     * @return @f$K_H@f$
     */
    double getKH() const
    {
        return KH;
    }

    /**
     * @brief A set method to change the factor rescaling the total Higgs width
     * with respect to the SM @f$K_H@f$.
     * @param[in] KH @f$K_H@f$ the factor rescaling the total Higgs width.
     */
    void setKH(double KH)
    {
        this->KH = KH;
    }
    
    
    /**
     * @brief A get method to retrieve the invisible Higgs branching ratio.
     * @return @f$BR_{inv}@f$ the invisible Higgs branching ratio.
     */
    double getBrHinv() const
    {
        return BrHinv;
    }

    /**
     * @brief A set method to change the invisible Higgs branching ratio.
     * @param[in] BrHinv @f$BR_{inv}@f$ the invisible Higgs branching ratio.
     */
    void setBrHinv(double BrHinv)
    {
        this->BrHinv = BrHinv;
    } 
    
    /**
     * @brief A get method to retrieve the exotic (not invisible) Higgs branching ratio.
     * @return @f$BR_{exo}@f$ the exotic (not invisible) Higgs branching ratio.
     */
    double getBrHexo() const
    {
        return BrHexo;
    }
    
    /**
     * @brief A set method to change the exotic (not invisible) Higgs branching ratio.
     * @param[in] BrHexo @f$BR_{exo}@f$ the exotic (not invisible) Higgs branching ratio.
     */
    void setBrHexo(double BrHexo)
    {
        this->BrHexo = BrHexo;
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
    virtual const double obliqueS() const;

    /**
     * @brief The oblique parameter @f$T@f$.
     * @return @f$T@f$
     */
    virtual const double obliqueT() const;

    /**
     * @brief The oblique parameter @f$U@f$.
     * @return @f$U=0@f$
     */
    virtual const double obliqueU() const;

    /**
     * @brief The ratio @f$\mu_{ggH}@f$ between the gluon-gluon fusion Higgs
     * production cross-section in the current model and in the Standard Model.
     * @param[in] sqrt_s the center-of-mass energy in TeV
     * @return @f$\mu_{ggH}@f$
     */
    virtual const double muggH(const double sqrt_s) const;
    /**
     * @brief The ratio @f$\mu_{VBF}@f$ between the vector-boson fusion Higgs
     * production cross-section in the current model and in the Standard Model. 
     * @param[in] sqrt_s the center-of-mass energy in TeV
     * @return @f$\mu_{VBF}@f$
     */
    virtual const double muVBF(const double sqrt_s) const;
    /**
     * @brief The ratio @f$\mu_{VBF+\gamma}@f$ between the vector-boson fusion Higgs
     * production cross-section in association with a hard photon in the current model
     * and in the Standard Model. 
     * @param[in] sqrt_s the center-of-mass energy in TeV
     * @return @f$\mu_{VBF+\gamma}@f$
     */
    virtual const double muVBFgamma(const double sqrt_s) const;
    /**
     * @brief The ratio @f$\mu_{eeWBF}@f$ between the 
     * @f$ e^{+}e^{-}\to \nu\bar{\nu} H @f$ production
     * cross-section in the current model and in the Standard Model.
     * @param[in] sqrt_s the center-of-mass energy in TeV
     * @return @f$\mu_{eeWBF}@f$
     */
    virtual const double mueeWBF(const double sqrt_s) const;
    /**
     * @brief The ratio @f$\mu_{eeWBF}@f$ between the 
     * @f$ e^{+}e^{-}\to \nu\bar{\nu} H @f$ production
     * cross-section in the current model and in the Standard Model.
     * @param[in] sqrt_s the center-of-mass energy in TeV,  
     * @param[in] Pol_em polarization of electrons
     * @param[in ]Pol_ep polarization of positrons
     * @return @f$\mu_{eeWBF}@f$
     */
    virtual const double mueeWBFPol(const double sqrt_s, const double Pol_em, const double Pol_ep) const;
    /**
     * @brief The ratio @f$\mu_{e^+e^- \to H\nu\bar{\nu}}@f$ between the 
     * @f$ e^+e^- \to H\nu\bar{\nu} @f$ associated production
     * cross-section in the current model and in the Standard Model.
     * @param[in] sqrt_s the center-of-mass energy in TeV
     * @return @f$\mu_{e^+e^- \to H\nu\bar{\nu}}@f$
     */
    virtual const double mueeHvv(const double sqrt_s) const;
    /**
     * @brief The ratio @f$\mu_{e^+e^- \to H\nu\bar{\nu}}@f$ between the 
     * @f$ e^+e^- \to H\nu\bar{\nu} @f$ associated production
     * cross-section in the current model and in the Standard Model.
     * @param[in] sqrt_s the center-of-mass energy in TeV, 
     * @param[in] Pol_em polarization of electrons
     * @param[in ]Pol_ep polarization of positrons
     * @return @f$\mu_{e^+e^- \to H\nu\bar{\nu}}@f$
     */
    virtual const double mueeHvvPol(const double sqrt_s, const double Pol_em, const double Pol_ep) const;
    /**
     * @brief The ratio @f$\mu_{eeZBF}@f$ between the 
     * @f$ e^{+}e^{-}\to e^{+}e^{-} H @f$ production
     * cross-section in the current model and in the Standard Model.
     * @param[in] sqrt_s the center-of-mass energy in TeV
     * @return @f$\mu_{eeZBF}@f$
     */
    virtual const double mueeZBF(const double sqrt_s) const;
    /**
     * @brief The ratio @f$\mu_{eeZBF}@f$ between the 
     * @f$ e^{+}e^{-}\to e^{+}e^{-} H @f$ production
     * cross-section in the current model and in the Standard Model.
     * @param[in] sqrt_s the center-of-mass energy in TeV,  
     * @param[in] Pol_em polarization of electrons
     * @param[in ]Pol_ep polarization of positrons
     * @return @f$\mu_{eeZBF}@f$
     */
    virtual const double mueeZBFPol(const double sqrt_s, const double Pol_em, const double Pol_ep) const;
    /**
     * @brief The ratio @f$\mu_{epWBF}@f$ between the 
     * @f$ e^{-} p\to \nu j H @f$ production
     * cross-section in the current model and in the Standard Model.
     * @param[in] sqrt_s the center-of-mass energy in TeV
     * @return @f$\mu_{epWBF}@f$
     */
    virtual const double muepWBF(const double sqrt_s) const;
    /**
     * @brief The ratio @f$\mu_{epZBF}@f$ between the 
     * @f$ e^{-} p\to e^{-} j H @f$ production
     * cross-section in the current model and in the Standard Model.
     * @param[in] sqrt_s the center-of-mass energy in TeV
     * @return @f$\mu_{epZBF}@f$
     */
    virtual const double muepZBF(const double sqrt_s) const;
    /**
     * @brief The ratio @f$\mu_{WH}@f$ between the W-Higgs associated production
     * cross-section in the current model and in the Standard Model. 
     * @param[in] sqrt_s the center-of-mass energy in TeV
     * @return @f$\mu_{WH}@f$
     */
    virtual const double muWH(const double sqrt_s) const;
    /**
     * @brief The ratio @f$\mu_{WH}@f$ between the W-Higgs associated production
     * cross-section in the current model and in the Standard Model, with @f$p_{T,H}>250@f$ GeV.
     * @param[in] sqrt_s the center-of-mass energy in TeV
     * @return @f$\mu_{WH}@f$
     */
    virtual const double muWHpT250(const double sqrt_s) const;
    /**
     * @brief The ratio @f$\mu_{ZH}@f$ between the Z-Higgs associated production
     * cross-section in the current model and in the Standard Model.
     * @param[in] sqrt_s the center-of-mass energy in TeV
     * @return @f$\mu_{ZH}@f$
     */
    virtual const double muZH(const double sqrt_s) const;
    /**
     * @brief The ratio @f$\mu_{ZH}@f$ between the Z-Higgs associated production
     * cross-section in the current model and in the Standard Model, with @f$p_{T,H}>250@f$ GeV.
     * @param[in] sqrt_s the center-of-mass energy in TeV
     * @return @f$\mu_{ZH}@f$
     */
    virtual const double muZHpT250(const double sqrt_s) const;
    /**
     * @brief The ratio @f$\mu_{eeZH}@f$ between the 
     * @f$e^{+}e^{-}\to ZH@f$ associated production
     * cross-section in the current model and in the Standard Model.
     * @param[in] sqrt_s the center-of-mass energy in TeV, Pol_em and Pol_ep
     * are the polarization of electrons and positrons, respectively
     * @return @f$\mu_{eeZH}@f$
     */
    virtual const double mueeZH(const double sqrt_s, const double Pol_em, const double Pol_ep) const;
    /**
     * @brief The ratio @f$\mu_{eeZH, Z \to e^+ e^-, \mu^+ \mu^-}@f$ between the 
     * @f$ e^{+}e^{-}\to ZH, Z \to e^+ e^-, \mu^+ \mu^- @f$ associated production
     * cross-section in the current model and in the Standard Model.
     * @param[in] sqrt_s the center-of-mass energy in TeV
     * @return @f$\mu_{eeZH, Z \to e^+ e^-, \mu^+ \mu^-}@f$
     */
    virtual const double mueeZllH(const double sqrt_s) const;
    /**
     * @brief The ratio @f$\mu_{eeZH, Z \to q \bar{q}}@f$ between the 
     * @f$ e^{+}e^{-}\to ZH, Z \to q \bar{q} @f$ associated production
     * cross-section in the current model and in the Standard Model.
     * @param[in] sqrt_s the center-of-mass energy in TeV
     * @return @f$\mu_{eeZH, Z \to q \bar{q}}@f$
     */
    virtual const double mueeZqqH(const double sqrt_s) const; 
    /**
     * @brief The ratio @f$\mu_{eeZH, Z \to e^+ e^-, \mu^+ \mu^-}@f$ between the 
     * @f$ e^{+}e^{-}\to ZH, Z \to e^+ e^-, \mu^+ \mu^- @f$ associated production
     * cross-section in the current model and in the Standard Model.
     * @param[in] sqrt_s the center-of-mass energy in TeV,  
     * @param[in] Pol_em polarization of electrons
     * @param[in ]Pol_ep polarization of positrons
     * @return @f$\mu_{eeZH, Z \to e^+ e^-, \mu^+ \mu^-}@f$
     */
    virtual const double mueeZllHPol(const double sqrt_s, const double Pol_em, const double Pol_ep) const;    
    /**
     * @brief The ratio @f$\mu_{eeZH, Z \to q \bar{q}}@f$ between the 
     * @f$ e^{+}e^{-}\to ZH, Z \to q \bar{q} @f$ associated production
     * cross-section in the current model and in the Standard Model.
     * @param[in] sqrt_s the center-of-mass energy in TeV,  
     * @param[in] Pol_em polarization of electrons
     * @param[in ]Pol_ep polarization of positrons 
     * @return @f$\mu_{eeZH, Z \to q \bar{q}}@f$
     */
    virtual const double mueeZqqHPol(const double sqrt_s, const double Pol_em, const double Pol_ep) const;    
    /**
     * @brief The ratio @f$\mu_{VH}@f$ between the WH+ZH associated production
     * cross-section in the current model and in the Standard Model.
     * @param[in] sqrt_s the center-of-mass energy in TeV
     * @return @f$\mu_{VH}@f$
     */
    virtual const double muVH(const double sqrt_s) const;
    /**
     * @brief The ratio @f$\mu_{VH}@f$ between the WH+ZH associated production
     * cross-section in the current model and in the Standard Model, with @f$p_{T,H}>250@f$ GeV.
     * @param[in] sqrt_s the center-of-mass energy in TeV
     * @return @f$\mu_{VH}@f$
     */
    virtual const double muVHpT250(const double sqrt_s) const;
    /**
     * @brief The ratio @f$\mu_{VBF+VH}@f$ between the sum of VBF and WH+ZH associated production
     * cross-section in the current model and in the Standard Model.
     * @param[in] sqrt_s the center-of-mass energy in TeV
     * @return @f$\mu_{VBF+VH}@f$
     */
    virtual const double muVBFpVH(const double sqrt_s) const;
    /**
     * @brief The ratio @f$\mu_{ttH}@f$ between the t-tbar-Higgs associated 
     * production cross-section in the current model and in the Standard Model.
     * @param[in] sqrt_s the center-of-mass energy in TeV
     * @return @f$\mu_{ttH}@f$
     */
    virtual const double muttH(const double sqrt_s) const;
    /**
     * @brief The ratio @f$\mu_{tHq}@f$ between the t-q-Higgs associated 
     * production cross-section in the current model and in the Standard Model.
     * @param[in] sqrt_s the center-of-mass energy in TeV
     * @return @f$\mu_{tHq}@f$
     */
    virtual const double mutHq(const double sqrt_s) const;
     /**
     * @brief The ratio @f$\mu_{\mu\mu H}@f$ between the @f$\sigma(\mu \mu \to H)}@f$
     * production cross-section in the current model and in the Standard Model.
     * @param[in] sqrt_s the center-of-mass energy in TeV
     * @return @f$\mu_{\mu\mu H}@f$
     */
    virtual const double mummH(const double sqrt_s) const;
    /**
     * @brief The ratio @f$\mu_{\mu\mu H}@f$ between the @f$\sigma(\mu \mu \to H)}@f$
     * production cross-section in the current model and in the Standard Model, 
     * in the narrow width approximation.
     * @param[in] sqrt_s the center-of-mass energy in TeV
     * @return @f$\mu_{\mu\mu H}@f$
     */
    virtual const double mummHNWA(const double sqrt_s) const; 
    /**
     * @brief The ratio @f$\mu_{\mu\mu ZH}@f$ between the @f$\sigma(\mu \mu \to Z H)}@f$
     * production cross-section in the current model and in the Standard Model.
     * @param[in] sqrt_s the center-of-mass energy in TeV
     * @return @f$\mu_{\mu\mu ZH}@f$
     */
    virtual const double mummZH(const double sqrt_s) const;
    /**
     * @brief The ratio @f$\mu_{\mu\mu H\nu\nu}@f$ between the @f$\sigma(\mu \mu \to H \nu \nu)}@f$
     * production cross-section in the current model and in the Standard Model.
     * @param[in] sqrt_s the center-of-mass energy in TeV
     * @return @f$\mu_{\mu\mu H\nu\nu}@f$
     */
    virtual const double mummHvv(const double sqrt_s) const;
    /**
     * @brief The ratio @f$\mu_{\mu\mu H\mu\mu}@f$ between the @f$\sigma(\mu \mu \to H \mu \mu)}@f$
     * production cross-section in the current model and in the Standard Model.
     * @param[in] sqrt_s the center-of-mass energy in TeV
     * @return @f$\mu_{\mu\mu H\mu\mu}@f$
     */
    virtual const double mummHmm(const double sqrt_s) const;
    /**
     * @brief The ratio @f$\mu_{\mu\mu ttH}@f$ between the @f$\sigma(\mu \mu \to t\bar{t} H )}@f$
     * production cross-section in the current model and in the Standard Model.
     * @param[in] sqrt_s the center-of-mass energy in TeV
     * @return @f$\mu_{\mu\mu ttH}@f$
     */
    virtual const double mummttH(const double sqrt_s) const;
    /**
     * @brief The ratio @f$\mu_{ggH+ttH}@f$ between the sum of gluon-gluon fusion
     * and t-tbar-Higgs associated 
     * production cross-section in the current model and in the Standard Model.
     * @param[in] sqrt_s the center-of-mass energy in TeV
     * @return @f$\mu_{ggH+ttH}@f$
     */
    virtual const double muggHpttH(const double sqrt_s) const;
    /**
     * @brief The ratio @f$\mu_{eettH}@f$ between the 
     * @f$ e^{+}e^{-}\to t\bar{t} H @f$ production
     * cross-section in the current model and in the Standard Model.
     * @param[in] sqrt_s the center-of-mass energy in TeV
     * @return @f$\mu_{eettH}@f$
     */
    virtual const double mueettH(const double sqrt_s) const;
    /**
     * @brief The ratio @f$\mu_{eettH}@f$ between the 
     * @f$ e^{+}e^{-}\to t\bar{t} H @f$ production
     * cross-section in the current model and in the Standard Model.
     * @param[in] sqrt_s the center-of-mass energy in TeV,  
     * @param[in] Pol_em polarization of electrons
     * @param[in ]Pol_ep polarization of positrons
     * @return @f$\mu_{eettH}@f$
     */
    virtual const double mueettHPol(const double sqrt_s, const double Pol_em, const double Pol_ep) const;
    /**
     * @brief The ratio @f$\sigma(ttH)/\sigma(ttZ)@f$ 
     * in the @f$H,Z\to b\bar{b}@f$ channel in the current model and in the Standard Model.
     * @param[in] sqrt_s the center-of-mass energy in TeV
     * @return @f$\sigma(ttH)/\sigma(ttZ)@f$ normalized to the SM
     */
    virtual const double muttHZbbboost(const double sqrt_s) const;

    
    ///////////////////////FULL SIGNAL STRENGTHS/////////////////////////
    
    /**
     * @brief The ratio @f$\mu_{ggH,\gamma\gamma}@f$ between the gluon-gluon fusion Higgs
     * production cross-section with subsequent decay into 2 photons in the
     * current model and in the Standard Model.
     * @param[in] sqrt_s the center-of-mass energy in TeV
     * @return @f$\mu_{ggH,\gamma\gamma}@f$
     */    
    virtual const double muggHgaga(const double sqrt_s) const;    
    /**
     * @brief The ratio @f$\mu_{ggH,\gamma\gamma}@f$ between the gluon-gluon fusion Higgs
     * production cross-section with subsequent decay into 2 photons in the
     * current model and in the Standard Model. Includes interference effects
     * with the background, following arXiv:1704.08259
     * @param[in] sqrt_s the center-of-mass energy in TeV
     * @return @f$\mu_{ggH,\gamma\gamma}@f$
     */
    virtual const double muggHgagaInt(const double sqrt_s) const;
    /**
     * @brief The ratio @f$\mu_{VBF,\gamma\gamma}@f$ between the VBF Higgs
     * production cross-section with subsequent decay into 2 photons in the
     * current model and in the Standard Model.
     * @param[in] sqrt_s the center-of-mass energy in TeV
     * @return @f$\mu_{VBF,\gamma\gamma}@f$
     */    
    virtual const double muVBFHgaga(const double sqrt_s) const;
    /**
     * @brief The ratio @f$\mu_{ZH,\gamma\gamma}@f$ between the ZH
     * production cross-section with subsequent decay into 2 photons in the
     * current model and in the Standard Model.
     * @param[in] sqrt_s the center-of-mass energy in TeV
     * @return @f$\mu_{ZH,\gamma\gamma}@f$
     */
    virtual const double muZHgaga(const double sqrt_s) const;
    /**
     * @brief The ratio @f$\mu_{WH,\gamma\gamma}@f$ between the WH
     * production cross-section with subsequent decay into 2 photons in the
     * current model and in the Standard Model.
     * @param[in] sqrt_s the center-of-mass energy in TeV
     * @return @f$\mu_{WH,\gamma\gamma}@f$
     */
    virtual const double muWHgaga(const double sqrt_s) const;
    /**
     * @brief The ratio @f$\mu_{VH,\gamma\gamma}@f$ between the VH
     * production cross-section with subsequent decay into 2 photons in the
     * current model and in the Standard Model.
     * @param[in] sqrt_s the center-of-mass energy in TeV
     * @return @f$\mu_{VH,\gamma\gamma}@f$
     */
    virtual const double muVHgaga(const double sqrt_s) const;
    /**
     * @brief The ratio @f$\mu_{ttH,\gamma\gamma}@f$ between the ttH
     * production cross-section with subsequent decay into 2 photons in the
     * current model and in the Standard Model.
     * @param[in] sqrt_s the center-of-mass energy in TeV
     * @return @f$\mu_{ttH,\gamma\gamma}@f$
     */
    virtual const double muttHgaga(const double sqrt_s) const;
    /**
     * @brief The ratio @f$\mu_{ggH,Z\gamma}@f$ between the gluon-gluon fusion Higgs
     * production cross-section with subsequent decay into @f$Z \gamma@f$ in the
     * current model and in the Standard Model.
     * @param[in] sqrt_s the center-of-mass energy in TeV
     * @return @f$\mu_{ggH,Z\gamma}@f$
     */
    virtual const double muggHZga(const double sqrt_s) const;
    /**
     * @brief The ratio @f$\mu_{VBF,Z\gamma}@f$ between the VBF Higgs
     * production cross-section with subsequent decay into @f$Z \gamma@f$ in the
     * current model and in the Standard Model.
     * @param[in] sqrt_s the center-of-mass energy in TeV
     * @return @f$\mu_{VBF,Z\gamma}@f$
     */
    virtual const double muVBFHZga(const double sqrt_s) const;    
    /**
     * @brief The ratio @f$\mu_{ZH,Z\gamma}@f$ between the ZH
     * production cross-section with subsequent decay into @f$Z \gamma@f$ in the
     * current model and in the Standard Model.
     * @param[in] sqrt_s the center-of-mass energy in TeV
     * @return @f$\mu_{ZH,Z\gamma}@f$
     */
    virtual const double muZHZga(const double sqrt_s) const;
    /**
     * @brief The ratio @f$\mu_{WH,Z\gamma}@f$ between the WH
     * production cross-section with subsequent decay into @f$Z \gamma@f$ in the
     * current model and in the Standard Model.
     * @param[in] sqrt_s the center-of-mass energy in TeV
     * @return @f$\mu_{WH,Z\gamma}@f$
     */
    virtual const double muWHZga(const double sqrt_s) const;    
    /**
     * @brief The ratio @f$\mu_{VH,Z\gamma}@f$ between the VH
     * production cross-section with subsequent decay into @f$Z \gamma@f$ in the
     * current model and in the Standard Model.
     * @param[in] sqrt_s the center-of-mass energy in TeV
     * @return @f$\mu_{VH,Z\gamma}@f$
     */
    virtual const double muVHZga(const double sqrt_s) const;
    /**
     * @brief The ratio @f$\mu_{ttH,Z\gamma}@f$ between the ttH
     * production cross-section with subsequent decay into @f$Z \gamma@f$ in the
     * current model and in the Standard Model.
     * @param[in] sqrt_s the center-of-mass energy in TeV
     * @return @f$\mu_{ttH,Z\gamma}@f$
     */
    virtual const double muttHZga(const double sqrt_s) const; 
    /**
     * @brief The ratio @f$\mu_{ggH,ZZ}@f$ between the gluon-gluon fusion Higgs
     * production cross-section with subsequent decay into @f$Z Z^*@f$ in the
     * current model and in the Standard Model.
     * @param[in] sqrt_s the center-of-mass energy in TeV
     * @return @f$\mu_{ggH,ZZ}@f$
     */
    virtual const double muggHZZ(const double sqrt_s) const;
    /**
     * @brief The ratio @f$\mu_{VBF,ZZ}@f$ between the VBF Higgs
     * production cross-section with subsequent decay into @f$Z Z^*@f$ in the
     * current model and in the Standard Model.
     * @param[in] sqrt_s the center-of-mass energy in TeV
     * @return @f$\mu_{VBF,ZZ}@f$
     */
    virtual const double muVBFHZZ(const double sqrt_s) const;
    /**
     * @brief The ratio @f$\mu_{ZH,ZZ}@f$ between the ZH
     * production cross-section with subsequent decay into @f$Z Z^*@f$ in the
     * current model and in the Standard Model.
     * @param[in] sqrt_s the center-of-mass energy in TeV
     * @return @f$\mu_{ZH,ZZ}@f$
     */
    virtual const double muZHZZ(const double sqrt_s) const;
    /**
     * @brief The ratio @f$\mu_{WH,ZZ}@f$ between the WH
     * production cross-section with subsequent decay into @f$Z Z^*@f$ in the
     * current model and in the Standard Model.
     * @param[in] sqrt_s the center-of-mass energy in TeV
     * @return @f$\mu_{WH,ZZ}@f$
     */
    virtual const double muWHZZ(const double sqrt_s) const;
    /**
     * @brief The ratio @f$\mu_{VH,ZZ}@f$ between the VH
     * production cross-section with subsequent decay into @f$Z Z^*@f$ in the
     * current model and in the Standard Model.
     * @param[in] sqrt_s the center-of-mass energy in TeV
     * @return @f$\mu_{VH,ZZ}@f$
     */
    virtual const double muVHZZ(const double sqrt_s) const;
    /**
     * @brief The ratio @f$\mu_{ttH,ZZ}@f$ between the ttH
     * production cross-section with subsequent decay into @f$Z Z^*@f$ in the
     * current model and in the Standard Model.
     * @param[in] sqrt_s the center-of-mass energy in TeV
     * @return @f$\mu_{ttH,ZZ}@f$
     */
    virtual const double muttHZZ(const double sqrt_s) const;
    
    /**
     * @brief The ratio @f$\mu_{ggH,ZZ\to 4l}@f$ between the gluon-gluon fusion Higgs
     * production cross-section with subsequent decay into @f$Z Z^*\to 4l@f$ in the
     * current model and in the Standard Model.
     * @param[in] sqrt_s the center-of-mass energy in TeV
     * @return @f$\mu_{ggH,ZZ\to 4l}@f$
     */
    virtual const double muggHZZ4l(const double sqrt_s) const;
    /**
     * @brief The ratio @f$\mu_{VBF,ZZ\to 4l}@f$ between the VBF Higgs
     * production cross-section with subsequent decay into @f$Z Z^*\to 4l@f$ in the
     * current model and in the Standard Model.
     * @param[in] sqrt_s the center-of-mass energy in TeV
     * @return @f$\mu_{VBF,ZZ\to 4l}@f$
     */
    virtual const double muVBFHZZ4l(const double sqrt_s) const;    
    /**
     * @brief The ratio @f$\mu_{ZH,ZZ\to 4l}@f$ between the ZH
     * production cross-section with subsequent decay into @f$Z Z^*\to 4l@f$ in the
     * current model and in the Standard Model.
     * @param[in] sqrt_s the center-of-mass energy in TeV
     * @return @f$\mu_{ZH,ZZ\to 4l}@f$
     */
    virtual const double muZHZZ4l(const double sqrt_s) const;
    /**
     * @brief The ratio @f$\mu_{WH,ZZ\to 4l}@f$ between the WH
     * production cross-section with subsequent decay into @f$Z Z^*\to 4l@f$ in the
     * current model and in the Standard Model.
     * @param[in] sqrt_s the center-of-mass energy in TeV
     * @return @f$\mu_{WH,ZZ\to 4l}@f$
     */
    virtual const double muWHZZ4l(const double sqrt_s) const;    
    /**
     * @brief The ratio @f$\mu_{VH,ZZ\to 4l}@f$ between the VH
     * production cross-section with subsequent decay into @f$Z Z^*\to 4l@f$ in the
     * current model and in the Standard Model.
     * @param[in] sqrt_s the center-of-mass energy in TeV
     * @return @f$\mu_{VH,ZZ\to 4l}@f$
     */
    virtual const double muVHZZ4l(const double sqrt_s) const;
    /**
     * @brief The ratio @f$\mu_{ttH,ZZ\to 4l}@f$ between the ttH
     * production cross-section with subsequent decay into @f$Z Z^*\to 4l@f$ in the
     * current model and in the Standard Model.
     * @param[in] sqrt_s the center-of-mass energy in TeV
     * @return @f$\mu_{ttH,ZZ\to 4l}@f$
     */
    virtual const double muttHZZ4l(const double sqrt_s) const;
    
    /**
     * @brief The ratio @f$\mu_{ggH,WW}@f$ between the gluon-gluon fusion Higgs
     * production cross-section with subsequent decay into @f$W W^*@f$ in the
     * current model and in the Standard Model.
     * @param[in] sqrt_s the center-of-mass energy in TeV
     * @return @f$\mu_{ggH,WW}@f$
     */
    virtual const double muggHWW(const double sqrt_s) const;
    /**
     * @brief The ratio @f$\mu_{VBF,WW}@f$ between the VBF Higgs
     * production cross-section with subsequent decay into @f$W W^*@f$ in the
     * current model and in the Standard Model.
     * @param[in] sqrt_s the center-of-mass energy in TeV
     * @return @f$\mu_{VBF,WW}@f$
     */
    virtual const double muVBFHWW(const double sqrt_s) const;
    /**
     * @brief The ratio @f$\mu_{ZH,WW}@f$ between the ZH
     * production cross-section with subsequent decay into @f$W W^*@f$ in the
     * current model and in the Standard Model.
     * @param[in] sqrt_s the center-of-mass energy in TeV
     * @return @f$\mu_{ZH,WW}@f$
     */
    virtual const double muZHWW(const double sqrt_s) const;
    /**
     * @brief The ratio @f$\mu_{WH,WW}@f$ between the WH
     * production cross-section with subsequent decay into @f$W W^*@f$ in the
     * current model and in the Standard Model.
     * @param[in] sqrt_s the center-of-mass energy in TeV
     * @return @f$\mu_{WH,WW}@f$
     */
    virtual const double muWHWW(const double sqrt_s) const; 
    /**
     * @brief The ratio @f$\mu_{VH,WW}@f$ between the VH
     * production cross-section with subsequent decay into @f$W W^*@f$ in the
     * current model and in the Standard Model.
     * @param[in] sqrt_s the center-of-mass energy in TeV
     * @return @f$\mu_{VH,WW}@f$
     */
    virtual const double muVHWW(const double sqrt_s) const;
    /**
     * @brief The ratio @f$\mu_{ttH,WW}@f$ between the ttH
     * production cross-section with subsequent decay into @f$W W^*@f$ in the
     * current model and in the Standard Model.
     * @param[in] sqrt_s the center-of-mass energy in TeV
     * @return @f$\mu_{ttH,WW}@f$
     */
    virtual const double muttHWW(const double sqrt_s) const;
    /**
     * @brief The ratio @f$\mu_{ggH,WW\to 2l2\nu}@f$ between the gluon-gluon fusion Higgs
     * production cross-section with subsequent decay into @f$W W^*\to 2l2\nu@f$ in the
     * current model and in the Standard Model.
     * @param[in] sqrt_s the center-of-mass energy in TeV
     * @return @f$\mu_{ggH,WW\to 2l2\nu}@f$
     */
    virtual const double muggHWW2l2v(const double sqrt_s) const;
    /**
     * @brief The ratio @f$\mu_{VBF,WW\to 2l2\nu}@f$ between the VBF Higgs
     * production cross-section with subsequent decay into @f$W W^*\to 2l2\nu@f$ in the
     * current model and in the Standard Model.
     * @param[in] sqrt_s the center-of-mass energy in TeV
     * @return @f$\mu_{VBF,WW\to 2l2\nu}@f$
     */
    virtual const double muVBFHWW2l2v(const double sqrt_s) const;   
    /**
     * @brief The ratio @f$\mu_{ZH,WW\to 2l2\nu}@f$ between the ZH
     * production cross-section with subsequent decay into @f$W W^*\to 2l2\nu@f$ in the
     * current model and in the Standard Model.
     * @param[in] sqrt_s the center-of-mass energy in TeV
     * @return @f$\mu_{ZH,WW\to 2l2\nu}@f$
     */
    virtual const double muZHWW2l2v(const double sqrt_s) const;
    /**
     * @brief The ratio @f$\mu_{WH,WW\to 2l2\nu}@f$ between the WH
     * production cross-section with subsequent decay into @f$W W^*\to 2l2\nu@f$ in the
     * current model and in the Standard Model.
     * @param[in] sqrt_s the center-of-mass energy in TeV
     * @return @f$\mu_{WH,WW\to 2l2\nu}@f$
     */
    virtual const double muWHWW2l2v(const double sqrt_s) const;    
    /**
     * @brief The ratio @f$\mu_{VH,WW\to 2l2\nu}@f$ between the VH
     * production cross-section with subsequent decay into @f$W W^*\to 2l2\nu@f$ in the
     * current model and in the Standard Model.
     * @param[in] sqrt_s the center-of-mass energy in TeV
     * @return @f$\mu_{VH,WW\to 2l2\nu}@f$
     */
    virtual const double muVHWW2l2v(const double sqrt_s) const;
    /**
     * @brief The ratio @f$\mu_{ttH,WW\to 2l2\nu}@f$ between the ttH
     * production cross-section with subsequent decay into @f$W W^*\to 2l2\nu@f$ in the
     * current model and in the Standard Model.
     * @param[in] sqrt_s the center-of-mass energy in TeV
     * @return @f$\mu_{ttH,WW\to 2l2\nu}@f$
     */
    virtual const double muttHWW2l2v(const double sqrt_s) const;
    /**
     * @brief The ratio @f$\mu_{ggH,\mu\mu}@f$ between the gluon-gluon fusion Higgs
     * production cross-section with subsequent decay into @f$\mu\mu@f$ in the
     * current model and in the Standard Model.
     * @param[in] sqrt_s the center-of-mass energy in TeV
     * @return @f$\mu_{ggH,\mu\mu}@f$
     */
    virtual const double muggHmumu(const double sqrt_s) const;
    /**
     * @brief The ratio @f$\mu_{VBF,\mu\mu}@f$ between the VBF Higgs
     * production cross-section with subsequent decay into @f$\mu\mu@f$ in the
     * current model and in the Standard Model.
     * @param[in] sqrt_s the center-of-mass energy in TeV
     * @return @f$\mu_{VBF,\mu\mu}@f$
     */
    virtual const double muVBFHmumu(const double sqrt_s) const;    
    /**
     * @brief The ratio @f$\mu_{ZH,\mu\mu}@f$ between the ZH
     * production cross-section with subsequent decay into @f$\mu\mu@f$ in the
     * current model and in the Standard Model.
     * @param[in] sqrt_s the center-of-mass energy in TeV
     * @return @f$\mu_{ZH,\mu\mu}@f$
     */
    virtual const double muZHmumu(const double sqrt_s) const;
    /**
     * @brief The ratio @f$\mu_{WH,\mu\mu}@f$ between the WH
     * production cross-section with subsequent decay into @f$\mu\mu@f$ in the
     * current model and in the Standard Model.
     * @param[in] sqrt_s the center-of-mass energy in TeV
     * @return @f$\mu_{WH,\mu\mu}@f$
     */
    virtual const double muWHmumu(const double sqrt_s) const;    
    /**
     * @brief The ratio @f$\mu_{VH,\mu\mu}@f$ between the VH
     * production cross-section with subsequent decay into @f$\mu\mu@f$ in the
     * current model and in the Standard Model.
     * @param[in] sqrt_s the center-of-mass energy in TeV
     * @return @f$\mu_{VH,\mu\mu}@f$
     */
    virtual const double muVHmumu(const double sqrt_s) const;
    /**
     * @brief The ratio @f$\mu_{ttH,\mu\mu}@f$ between the ttH
     * production cross-section with subsequent decay into @f$\mu\mu@f$ in the
     * current model and in the Standard Model.
     * @param[in] sqrt_s the center-of-mass energy in TeV
     * @return @f$\mu_{ttH,\mu\mu}@f$
     */
    virtual const double muttHmumu(const double sqrt_s) const;
    /**
     * @brief The ratio @f$\mu_{ggH,\tau\tau}@f$ between the gluon-gluon fusion Higgs
     * production cross-section with subsequent decay into @f$\tau\tau@f$ in the
     * current model and in the Standard Model.
     * @param[in] sqrt_s the center-of-mass energy in TeV
     * @return @f$\mu_{ggH,\tau\tau}@f$
     */
    virtual const double muggHtautau(const double sqrt_s) const;
    /**
     * @brief The ratio @f$\mu_{VBF,\tau\tau}@f$ between the VBF Higgs
     * production cross-section with subsequent decay into @f$\tau\tau@f$ in the
     * current model and in the Standard Model.
     * @param[in] sqrt_s the center-of-mass energy in TeV
     * @return @f$\mu_{VBF,\tau\tau}@f$
     */
    virtual const double muVBFHtautau(const double sqrt_s) const;
    /**
     * @brief The ratio @f$\mu_{ZH,\tau\tau}@f$ between the ZH
     * production cross-section with subsequent decay into @f$\tau\tau@f$ in the
     * current model and in the Standard Model.
     * @param[in] sqrt_s the center-of-mass energy in TeV
     * @return @f$\mu_{ZH,\tau\tau}@f$
     */
    virtual const double muZHtautau(const double sqrt_s) const;
    /**
     * @brief The ratio @f$\mu_{WH,\tau\tau}@f$ between the WH
     * production cross-section with subsequent decay into @f$\tau\tau@f$ in the
     * current model and in the Standard Model.
     * @param[in] sqrt_s the center-of-mass energy in TeV
     * @return @f$\mu_{WH,\tau\tau}@f$
     */
    virtual const double muWHtautau(const double sqrt_s) const;
    /**
     * @brief The ratio @f$\mu_{VH,\tau\tau}@f$ between the VH
     * production cross-section with subsequent decay into @f$\tau\tau@f$ in the
     * current model and in the Standard Model.
     * @param[in] sqrt_s the center-of-mass energy in TeV
     * @return @f$\mu_{VH,\tau\tau}@f$
     */
    virtual const double muVHtautau(const double sqrt_s) const;
    /**
     * @brief The ratio @f$\mu_{ttH,\tau\tau}@f$ between the ttH
     * production cross-section with subsequent decay into @f$\tau\tau@f$ in the
     * current model and in the Standard Model.
     * @param[in] sqrt_s the center-of-mass energy in TeV
     * @return @f$\mu_{ttH,\tau\tau}@f$
     */
    virtual const double muttHtautau(const double sqrt_s) const;
    /**
     * @brief The ratio @f$\mu_{ggH,bb}@f$ between the gluon-gluon fusion Higgs
     * production cross-section with subsequent decay into @f$bb@f$ in the
     * current model and in the Standard Model.
     * @param[in] sqrt_s the center-of-mass energy in TeV
     * @return @f$\mu_{ggH,bb}@f$
     */
    virtual const double muggHbb(const double sqrt_s) const;
    /**
     * @brief The ratio @f$\mu_{VBF,bb}@f$ between the VBF Higgs
     * production cross-section with subsequent decay into @f$bb@f$ in the
     * current model and in the Standard Model.
     * @param[in] sqrt_s the center-of-mass energy in TeV
     * @return @f$\mu_{VBF,bb}@f$
     */
    virtual const double muVBFHbb(const double sqrt_s) const;
    /**
     * @brief The ratio @f$\mu_{ZH,bb}@f$ between the ZH
     * production cross-section with subsequent decay into @f$bb@f$ in the
     * current model and in the Standard Model.
     * @param[in] sqrt_s the center-of-mass energy in TeV
     * @return @f$\mu_{ZH,bb}@f$
     */
    virtual const double muZHbb(const double sqrt_s) const;
    /**
     * @brief The ratio @f$\mu_{WH,bb}@f$ between the WH
     * production cross-section with subsequent decay into @f$bb@f$ in the
     * current model and in the Standard Model.
     * @param[in] sqrt_s the center-of-mass energy in TeV
     * @return @f$\mu_{WH,bb}@f$
     */
    virtual const double muWHbb(const double sqrt_s) const;
    /**
     * @brief The ratio @f$\mu_{VH,bb}@f$ between the VH
     * production cross-section with subsequent decay into @f$bb@f$ in the
     * current model and in the Standard Model.
     * @param[in] sqrt_s the center-of-mass energy in TeV
     * @return @f$\mu_{VH,bb}@f$
     */
    virtual const double muVHbb(const double sqrt_s) const;
    /**
     * @brief The ratio @f$\mu_{ttH,bb}@f$ between the ttH
     * production cross-section with subsequent decay into @f$bb@f$ in the
     * current model and in the Standard Model.
     * @param[in] sqrt_s the center-of-mass energy in TeV
     * @return @f$\mu_{ttH,bb}@f$
     */
    virtual const double muttHbb(const double sqrt_s) const;

    ///////////////////////HIGGS DECAY WIDTHS/////////////////////////
    
    /**
     * @brief The ratio of the @f$\Gamma(H)@f$ in the current model
     * and in the Standard Model.
     * @return @f$\Gamma(H)@f$/@f$\Gamma(H)_{\mathrm{SM}}@f$
     */
    virtual const double computeGammaTotalRatio() const;
    
    /**
     * @brief The ratio of the @f$\Gamma(H)@f$ in the current model
     * and in the Standard Model and neglecting theory errors.
     * @return @f$\Gamma(H)@f$/@f$\Gamma(H)_{\mathrm{SM}}@f$
     */
    virtual const double computeGammaTotalnoErrorRatio() const;
    
    /**
     * @brief The ratio of the @f$\Gamma(H\to SM)@f$ in the current model
     * and in the Standard Model, including only SM decays.
     * @return @f$\Gamma(H\to SM)@f$/@f$\Gamma(H\to SM)_{\mathrm{SM}}@f$
     */
    virtual const double computeGammaTotalSMRatio() const;
    
    /**
     * @brief The ratio of the @f$\Gamma(H\to SM)@f$ in the current model
     * and in the Standard Model, including only SM decays and neglecting theory errors.
     * @return @f$\Gamma(H\to SM)@f$/@f$\Gamma(H\to SM)_{\mathrm{SM}}@f$
     */
    virtual const double computeGammaTotalSMnoErrorRatio() const;
    
    ///////////////////////HIGGS BRANCHING RATIOS/////////////////////////
    
    /**
     * @brief The ratio of the Br@f$(H\to gg)@f$ in the current model
     * and in the Standard Model.
     * @return Br@f$(H\to gg)@f$/Br@f$(H\to gg)_{\mathrm{SM}}@f$
     */
    virtual const double BrHggRatio() const;
    /**
     * @brief The ratio of the Br@f$(H\to WW)@f$ in the current model
     * and in the Standard Model.
     * @return Br@f$(H\to WW)@f$/Br@f$(H\to WW)_{\mathrm{SM}}@f$
     */
    virtual const double BrHWWRatio() const;
    /**
     * @brief The ratio of the Br@f$(H\to ZZ)@f$ in the current model
     * and in the Standard Model.
     * @return Br@f$(H\to ZZ)@f$/Br@f$(H\to ZZ)_{\mathrm{SM}}@f$
     */
    virtual const double BrHZZRatio() const;
    /**
     * @brief The ratio of the Br@f$(H\to VV)@f$ in the current model
     * and in the Standard Model.
     * @return Br@f$(H\to VV)@f$/Br@f$(H\to VV)_{\mathrm{SM}}@f$
     */
    virtual const double BrHVVRatio() const;
    /**
     * @brief The ratio of the Br@f$(H\to Z\gamma)@f$ in the current model
     * and in the Standard Model.
     * @return Br@f$(H\to Z\gamma)@f$/Br@f$(H\to Z\gamma)_{\mathrm{SM}}@f$
     */
    virtual const double BrHZgaRatio() const;
    /**
     * @brief The ratio of the Br@f$(H\to Z\gamma\to ll\gamma)@f$ (@f$l=e,\mu @f$) in the current model
     * and in the Standard Model.
     * @return Br@f$(H\to Z\gamma\to ll\gamma)@f$/Br@f$(H\to Z\gamma\to ll\gamma)_{\mathrm{SM}}@f$
     */
    virtual const double BrHZgallRatio() const;
    /**
     * @brief The ratio of the Br@f$(H\to Z\gamma\to ee\gamma)@f$ in the current model
     * and in the Standard Model.
     * @return Br@f$(H\to Z\gamma\to ee\gamma)@f$/Br@f$(H\to Z\gamma\to ee\gamma)_{\mathrm{SM}}@f$
     */
    virtual const double BrHZgaeeRatio() const;
    /**
     * @brief The ratio of the Br@f$(H\to Z\gamma\to \mu\mu\gamma)@f$ in the current model
     * and in the Standard Model.
     * @return Br@f$(H\to Z\gamma\to \mu\mu\gamma)@f$/Br@f$(H\to Z\gamma\to \mu\mu\gamma)_{\mathrm{SM}}@f$
     */
    virtual const double BrHZgamumuRatio() const;
    /**
     * @brief The ratio of the Br@f$(H\to \gamma\gamma)@f$ in the current model
     * and in the Standard Model.
     * @return Br@f$(H\to \gamma\gamma)@f$/Br@f$(H\to \gamma\gamma)_{\mathrm{SM}}@f$
     */
    virtual const double BrHgagaRatio() const;
    /**
     * @brief The ratio of the Br@f$(H\to \mu^+\mu^-)@f$ in the current model
     * and in the Standard Model.
     * @return Br@f$(H\to \mu^+\mu^-)@f$/Br@f$(H\to \mu^+\mu^-)_{\mathrm{SM}}@f$
     */
    virtual const double BrHmumuRatio() const;
    /**
     * @brief The ratio of the Br@f$(H\to \tau^+\tau^-)@f$ in the current model
     * and in the Standard Model.
     * @return Br@f$(H\to \tau^+\tau^-)@f$/Br@f$(H\to \tau^+\tau^-)_{\mathrm{SM}}@f$
     */
    virtual const double BrHtautauRatio() const;
    /**
     * @brief The ratio of the Br@f$(H\to c\bar{c})@f$ in the current model
     * and in the Standard Model.
     * @return Br@f$(H\to c\bar{c})@f$/Br@f$(H\to c\bar{c})_{\mathrm{SM}}@f$
     */
    virtual const double BrHccRatio() const;
    /**
     * @brief The ratio of the Br@f$(H\to s\bar{s})@f$ in the current model
     * and in the Standard Model.
     * @return Br@f$(H\to s\bar{s})@f$/Br@f$(H\to s\bar{s})_{\mathrm{SM}}@f$
     */
    virtual const double BrHssRatio() const;
    /**
     * @brief The ratio of the Br@f$(H\to b\bar{b})@f$ in the current model
     * and in the Standard Model.
     * @return Br@f$(H\to b\bar{b})@f$/Br@f$(H\to b\bar{b})_{\mathrm{SM}}@f$
     */
    virtual const double BrHbbRatio() const;
    
    /**
     * @brief The branching ratio of the of the Higgs into exotic particles (invisible or not).
     * @return Br@f$(H\to exotic)@f$
     */
    virtual const double Br_H_exo() const;
    
    /**
     * @brief The branching ratio of the of the Higgs into invisible particles.
     * @return Br@f$(H\to invisible)@f$
     */
    virtual const double Br_H_inv() const;
    
    /**
     * @brief The branching ratio of the of the Higgs into invisible particles 
     * (only invisible new particles).
     * @return Br@f$(H\to invisible,NP)@f$
     */
    virtual const double Br_H_inv_NP() const;
    
    /**
     * @brief The ratio of the Br@f$(H\to visible)@f$ in the current model
     * and in the Standard Model.
     * @return Br@f$(H\to visible)@f$/Br@f$(H\to visible)_{\mathrm{SM}}@f$
     */
    virtual const double BrHvisRatio() const;
    
    /**
     * @brief The ratio of the Br@f$(H\to invisible)@f$ in the current model
     * and in the Standard Model.
     * @return Br@f$(H\to invisible)@f$/Br@f$(H\to ZZ \to invisible)_{\mathrm{SM}}@f$
     */
    virtual const double BrHtoinvRatio() const;
    
    /////////////////////// HIGGS TO 4 FERMION DECAYS /////////////////////////   
    
    /**
     * @brief The ratio of the Br@f$(H\to 2L2L')@f$ (@f$L,L'=e,\mu,\tau@f$) in the current model
     * and in the Standard Model.
     * @return Br@f$(H\to 2L2L')@f$/Br@f$(H\to 2L2L')_{\mathrm{SM}}@f$
     */
    virtual const double BrH2L2LRatio() const;
    
    /**
     * @brief The ratio of the Br@f$(H\to 2e 2\mu)@f$ in the current model
     * and in the Standard Model.
     * @return Br@f$(H\to 2L2L)@f$/Br@f$(H\to 2e 2\mu)_{\mathrm{SM}}@f$
     */
    virtual const double BrH2e2muRatio() const;    

    /**
     * @brief The ratio of the Br@f$(H\to 2v2v)@f$ in the current model
     * and in the Standard Model.
     * @return Br@f$(H\to 2v2v)@f$/Br@f$(H\to 2v2v)_{\mathrm{SM}}@f$
     */
    virtual const double BrH2v2vRatio() const;

    /**
     * @brief The ratio of the Br@f$(H\to 2L2v)@f$ (@f$L=e,\mu,\tau@f$) in the current model
     * and in the Standard Model.
     * @return Br@f$(H\to 2L2v)@f$/Br@f$(H\to 2L2v)_{\mathrm{SM}}@f$
     */
    virtual const double BrH2L2vRatio() const;
    
    /**
     * @brief The ratio of the Br@f$(H\to 2L2v)@f$ (@f$L=e,\mu@f$) in the current model
     * and in the Standard Model.
     * @return Br@f$(H\to 2L2v)@f$/Br@f$(H\to 2L2v)_{\mathrm{SM}}@f$
     */
    virtual const double BrH2L2v2Ratio() const;
    
    /**
     * @brief The ratio of the Br@f$(H\to 2e2v)@f$ in the current model
     * and in the Standard Model.
     * @return Br@f$(H\to 2e2v)@f$/Br@f$(H\to 2e2v)_{\mathrm{SM}}@f$
     */
    virtual const double BrH2e2vRatio() const;
    
    /**
     * @brief The ratio of the Br@f$(H\to 2\mu 2v)@f$ in the current model
     * and in the Standard Model.
     * @return Br@f$(H\to 2\mu 2v)@f$/Br@f$(H\to 2\mu 2v)_{\mathrm{SM}}@f$
     */
    virtual const double BrH2mu2vRatio() const;
    
    /**
     * @brief The ratio of the Br@f$(H\to 2u2u)@f$ in the current model
     * and in the Standard Model.
     * @return Br@f$(H\to 2u2u)@f$/Br@f$(H\to 2u2u)_{\mathrm{SM}}@f$
     */
    virtual const double BrH2u2uRatio() const;
    
    /**
     * @brief The ratio of the Br@f$(H\to 2d2d)@f$ in the current model
     * and in the Standard Model.
     * @return Br@f$(H\to 2d2d)@f$/Br@f$(H\to 2d2d)_{\mathrm{SM}}@f$
     */
    virtual const double BrH2d2dRatio() const;
    
    /**
     * @brief The ratio of the Br@f$(H\to 2u2d)@f$ in the current model
     * and in the Standard Model.
     * @return Br@f$(H\to 2u2d)@f$/Br@f$(H\to 2u2d)_{\mathrm{SM}}@f$
     */
    virtual const double BrH2u2dRatio() const;
    
    /**
     * @brief The ratio of the Br@f$(H\to 2L2u)@f$ (@f$L=e,\mu,\tau@f$) in the current model
     * and in the Standard Model.
     * @return Br@f$(H\to 2L2u)@f$/Br@f$(H\to 2L2u)_{\mathrm{SM}}@f$
     */
    virtual const double BrH2L2uRatio() const;
    
    /**
     * @brief The ratio of the Br@f$(H\to 2L2d)@f$ (@f$L=e,\mu,\tau@f$) in the current model
     * and in the Standard Model.
     * @return Br@f$(H\to 2L2d)@f$/Br@f$(H\to 2L2d)_{\mathrm{SM}}@f$
     */
    virtual const double BrH2L2dRatio() const;
    
    /**
     * @brief The ratio of the Br@f$(H\to 2v2u)@f$ in the current model
     * and in the Standard Model.
     * @return Br@f$(H\to 2v2u)@f$/Br@f$(H\to 2v2u)_{\mathrm{SM}}@f$
     */
    virtual const double BrH2v2uRatio() const;
    
    /**
     * @brief The ratio of the Br@f$(H\to 2v2d)@f$ in the current model
     * and in the Standard Model.
     * @return Br@f$(H\to 2v2d)@f$/Br@f$(H\to 2v2d)_{\mathrm{SM}}@f$
     */
    virtual const double BrH2v2dRatio() const;
    
    /**
     * @brief The ratio of the Br@f$(H\to 4L)@f$ (@f$L=e,\mu,\tau@f$) in the current model
     * and in the Standard Model.
     * @return Br@f$(H\to 4L)@f$/Br@f$(H\to 4L)_{\mathrm{SM}}@f$
     */
    virtual const double BrH4LRatio() const;
    
    /**
     * @brief The ratio of the Br@f$(H\to 4L)@f$ (@f$L=e,\mu@f$) in the current model
     * and in the Standard Model.
     * @return Br@f$(H\to 4L)@f$/Br@f$(H\to 4L)_{\mathrm{SM}}@f$
     */
    virtual const double BrH4L2Ratio() const;
    
    /**
     * @brief The ratio of the Br@f$(H\to 4e)@f$ in the current model
     * and in the Standard Model.
     * @return Br@f$(H\to 4e)@f$/Br@f$(H\to 4e)_{\mathrm{SM}}@f$
     */
    virtual const double BrH4eRatio() const;
    
    /**
     * @brief The ratio of the Br@f$(H\to 4\mu)@f$ in the current model
     * and in the Standard Model.
     * @return Br@f$(H\to 4\mu)@f$/Br@f$(H\to 4\mu)_{\mathrm{SM}}@f$
     */
    virtual const double BrH4muRatio() const;
    
    /**
     * @brief The ratio of the Br@f$(H\to 4v)@f$ in the current model
     * and in the Standard Model.
     * @return Br@f$(H\to 4v)@f$/Br@f$(H\to 4v)_{\mathrm{SM}}@f$
     */
    virtual const double BrH4vRatio() const;
    
    /**
     * @brief The ratio of the Br@f$(H\to 4u)@f$ in the current model
     * and in the Standard Model.
     * @return Br@f$(H\to 4u)@f$/Br@f$(H\to 4u)_{\mathrm{SM}}@f$
     */
    virtual const double BrH4uRatio() const;
    
    /**
     * @brief The ratio of the Br@f$(H\to 4d)@f$ in the current model
     * and in the Standard Model.
     * @return Br@f$(H\to 4d)@f$/Br@f$(H\to 4d)_{\mathrm{SM}}@f$
     */
    virtual const double BrH4dRatio() const;
    
    /**
     * @brief The ratio of the Br@f$(H\to LvvL)@f$ (@f$L=e,\mu,\tau@f$) in the current model
     * and in the Standard Model.
     * @return Br@f$(H\to LvvL)@f$/Br@f$(H\to LvvL)_{\mathrm{SM}}@f$
     */
    virtual const double BrHLvvLRatio() const;
    
    /**
     * @brief The ratio of the Br@f$(H\to e\nu \mu\nu)@f$ in the current model
     * and in the Standard Model.
     * @return Br@f$(H\to e\nu \mu\nu)@f$/Br@f$(H\to e\nu \mu\nu)_{\mathrm{SM}}@f$
     */
    virtual const double BrHevmuvRatio() const;
    
    /**
     * @brief The ratio of the Br@f$(H\to uddu)@f$ in the current model
     * and in the Standard Model.
     * @return Br@f$(H\to uddu)@f$/Br@f$(H\to uddu)_{\mathrm{SM}}@f$
     */
    virtual const double BrHudduRatio() const;
    
    /**
     * @brief The ratio of the Br@f$(H\to Lvud)@f$ (@f$L=e,\mu,\tau@f$) in the current model
     * and in the Standard Model.
     * @return Br@f$(H\to Lvud)@f$/Br@f$(H\to Lvud)_{\mathrm{SM}}@f$
     */
    virtual const double BrHLvudRatio() const;
    
    /**
     * @brief The ratio of the Br@f$(H\to 2ud)@f$ in the current model
     * and in the Standard Model.
     * @return Br@f$(H\to 2ud)@f$/Br@f$(H\to 2ud)_{\mathrm{SM}}@f$
     */
    virtual const double BrH2udRatio() const;
    
    /**
     * @brief The ratio of the Br@f$(H\to 2Lv)@f$ (@f$L=e,\mu,\tau@f$) in the current model
     * and in the Standard Model.
     * @return Br@f$(H\to 2Lv)@f$/Br@f$(H\to 2Lv)_{\mathrm{SM}}@f$
     */
    virtual const double BrH2LvRatio() const;   
    
    /**
     * @brief The ratio of the Br@f$(H\to 2Lv)@f$ (@f$L=e,\mu@f$) in the current model
     * and in the Standard Model.
     * @return Br@f$(H\to 2Lv)@f$/Br@f$(H\to 2Lv)_{\mathrm{SM}}@f$
     */
    virtual const double BrH2Lv2Ratio() const;
    
    /**
     * @brief The ratio of the Br@f$(H\to 2ev)@f$ in the current model
     * and in the Standard Model.
     * @return Br@f$(H\to 2ev)@f$/Br@f$(H\to 2ev)_{\mathrm{SM}}@f$
     */
    virtual const double BrH2evRatio() const;    
    
    /**
     * @brief The ratio of the Br@f$(H\to 2ev)@f$ in the current model
     * and in the Standard Model.
     * @return Br@f$(H\to 2\mu v)@f$/Br@f$(H\to 2\mu v)_{\mathrm{SM}}@f$
     */
    virtual const double BrH2muvRatio() const;    
    
    /**
     * @brief The ratio of the Br@f$(H\to 4f)@f$ in the current model
     * and in the Standard Model.
     * @return Br@f$(H\to 4f)@f$/Br@f$(H\to 4f)_{\mathrm{SM}}@f$
     */
    virtual const double BrH4fRatio() const;
    
    // DECAYS INVOLVING ONLY ELECTRONS, MUONS OR NEUTRINOS IN THE FINAL STATES 

    /**
     * @brief The ratio of the Br@f$(H\to 4l)@f$ (@f$l=e,\mu@f$) in the current model
     * and in the Standard Model.
     * @return Br@f$(H\to 4l)@f$/Br@f$(H\to 4l)_{\mathrm{SM}}@f$
     */
    virtual const double BrH4lRatio() const;
    
    /**
     * @brief The ratio of the Br@f$(H\to 2l2v)@f$ (@f$l=e,\mu@f$) in the current model
     * and in the Standard Model.
     * @return Br@f$(H\to 2l2v)@f$/Br@f$(H\to 2l2v)_{\mathrm{SM}}@f$
     */
    virtual const double BrH2l2vRatio() const;
    
    ///////////////////////OTHER HIGGS 4 FERMION BRANCHING RATIOS/////////////////////////
    
    /**
     * @brief The ratio of the Br@f$(H\to l l j j)@f$, @f$l=e,\mu,~~j\not=b)@f$ in the current model
     * and in the Standard Model.
     * @return Br@f$(H\to l l j j)@f$/Br@f$(H\to l l j j)_{\mathrm{SM}}@f$
     */
    virtual const double BrHlljjRatio() const
    {
        return BrHZZRatio(); 
    }
    
    /**
     * @brief The ratio of the Br@f$(H\to l \nu j j)@f$, @f$l=e,\mu,~~j\not=b)@f$ in the current model
     * and in the Standard Model.
     * @return Br@f$(H\to l \nu j j)@f$/Br@f$(H\to l \nu j j)_{\mathrm{SM}}@f$
     */
    virtual const double BrHlvjjRatio() const
    {
        return BrHWWRatio(); 
    }
    
    /**
     * @brief The ratio of the Br@f$(H\to l \nu l \nu, l \nu j j)@f$, @f$l=e,\mu,~~j\not=b)@f$ in the current model
     * and in the Standard Model.
     * @return Br@f$(H\to l \nu l \nu, l \nu j j)@f$/Br@f$(H\to l \nu l \nu, l \nu j j)_{\mathrm{SM}}@f$
     */
    virtual const double BrHlv_lvorjjRatio() const
    {
        return BrHWWRatio();
    }   
    
    /**
     * @brief The ratio of the Br@f$(H\to l l \nu\nu, l l j j)@f$, @f$l=e,\mu,~~j\not=b)@f$ in the current model
     * and in the Standard Model.
     * @return Br@f$(H\to l l \nu\nu, l l j j)@f$/Br@f$(H\to l l \nu\nu, l l j j)_{\mathrm{SM}}@f$
     */
    virtual const double BrHll_vvorjjRatio() const
    {
        return BrHZZRatio();
    }
    
////////////////////////////////////////////////////////////////////////////////////////////
//-----------------------------------------------------------------------------------------
//-- Special Hadron collider signal strengths with separate full TH unc U(prod x decay) ---
//-----------------------------------------------------------------------------------------
//////////////////////////////////////////////////////////////////////////////////////////// 
    
    /**
     * @brief The ratio @f$\mu_{ggH,\gamma\gamma}@f$ between the gluon-gluon fusion Higgs
     * production cross-section with subsequent decay into 2 photons in the
     * current model and in the Standard Model.
     * @param[in] sqrt_s the center-of-mass energy in TeV
     * @return @f$\mu_{ggH,\gamma\gamma}@f$
     */    
    virtual const double muTHUggHgaga(const double sqrt_s) const;    
    /**
     * @brief The ratio @f$\mu_{VBF,\gamma\gamma}@f$ between the VBF Higgs
     * production cross-section with subsequent decay into 2 photons in the
     * current model and in the Standard Model.
     * @param[in] sqrt_s the center-of-mass energy in TeV
     * @return @f$\mu_{VBF,\gamma\gamma}@f$
     */    
    virtual const double muTHUVBFHgaga(const double sqrt_s) const;
    /**
     * @brief The ratio @f$\mu_{ZH,\gamma\gamma}@f$ between the ZH
     * production cross-section with subsequent decay into 2 photons in the
     * current model and in the Standard Model.
     * @param[in] sqrt_s the center-of-mass energy in TeV
     * @return @f$\mu_{ZH,\gamma\gamma}@f$
     */
    virtual const double muTHUZHgaga(const double sqrt_s) const;
    /**
     * @brief The ratio @f$\mu_{WH,\gamma\gamma}@f$ between the WH
     * production cross-section with subsequent decay into 2 photons in the
     * current model and in the Standard Model.
     * @param[in] sqrt_s the center-of-mass energy in TeV
     * @return @f$\mu_{WH,\gamma\gamma}@f$
     */
    virtual const double muTHUWHgaga(const double sqrt_s) const;
    /**
     * @brief The ratio @f$\mu_{VH,\gamma\gamma}@f$ between the VH
     * production cross-section with subsequent decay into 2 photons in the
     * current model and in the Standard Model.
     * @param[in] sqrt_s the center-of-mass energy in TeV
     * @return @f$\mu_{VH,\gamma\gamma}@f$
     */
    virtual const double muTHUVHgaga(const double sqrt_s) const;
    /**
     * @brief The ratio @f$\mu_{ttH,\gamma\gamma}@f$ between the ttH
     * production cross-section with subsequent decay into 2 photons in the
     * current model and in the Standard Model.
     * @param[in] sqrt_s the center-of-mass energy in TeV
     * @return @f$\mu_{ttH,\gamma\gamma}@f$
     */
    virtual const double muTHUttHgaga(const double sqrt_s) const;
    /**
     * @brief The ratio @f$\mu_{ggH,Z\gamma}@f$ between the gluon-gluon fusion Higgs
     * production cross-section with subsequent decay into @f$Z \gamma@f$ in the
     * current model and in the Standard Model.
     * @param[in] sqrt_s the center-of-mass energy in TeV
     * @return @f$\mu_{ggH,Z\gamma}@f$
     */
    virtual const double muTHUggHZga(const double sqrt_s) const;
    /**
     * @brief The ratio @f$\mu_{VBF,Z\gamma}@f$ between the VBF Higgs
     * production cross-section with subsequent decay into @f$Z \gamma@f$ in the
     * current model and in the Standard Model.
     * @param[in] sqrt_s the center-of-mass energy in TeV
     * @return @f$\mu_{VBF,Z\gamma}@f$
     */
    virtual const double muTHUVBFHZga(const double sqrt_s) const;    
    /**
     * @brief The ratio @f$\mu_{ZH,Z\gamma}@f$ between the ZH
     * production cross-section with subsequent decay into @f$Z \gamma@f$ in the
     * current model and in the Standard Model.
     * @param[in] sqrt_s the center-of-mass energy in TeV
     * @return @f$\mu_{ZH,Z\gamma}@f$
     */
    virtual const double muTHUZHZga(const double sqrt_s) const;
    /**
     * @brief The ratio @f$\mu_{WH,Z\gamma}@f$ between the WH
     * production cross-section with subsequent decay into @f$Z \gamma@f$ in the
     * current model and in the Standard Model.
     * @param[in] sqrt_s the center-of-mass energy in TeV
     * @return @f$\mu_{WH,Z\gamma}@f$
     */
    virtual const double muTHUWHZga(const double sqrt_s) const;    
    /**
     * @brief The ratio @f$\mu_{VH,Z\gamma}@f$ between the VH
     * production cross-section with subsequent decay into @f$Z \gamma@f$ in the
     * current model and in the Standard Model.
     * @param[in] sqrt_s the center-of-mass energy in TeV
     * @return @f$\mu_{VH,Z\gamma}@f$
     */
    virtual const double muTHUVHZga(const double sqrt_s) const;
    /**
     * @brief The ratio @f$\mu_{ttH,Z\gamma}@f$ between the ttH
     * production cross-section with subsequent decay into @f$Z \gamma@f$ in the
     * current model and in the Standard Model.
     * @param[in] sqrt_s the center-of-mass energy in TeV
     * @return @f$\mu_{ttH,Z\gamma}@f$
     */
    virtual const double muTHUttHZga(const double sqrt_s) const; 
    /**
     * @brief The ratio @f$\mu_{ggH,ZZ}@f$ between the gluon-gluon fusion Higgs
     * production cross-section with subsequent decay into @f$Z Z^*@f$ in the
     * current model and in the Standard Model.
     * @param[in] sqrt_s the center-of-mass energy in TeV
     * @return @f$\mu_{ggH,ZZ}@f$
     */
    virtual const double muTHUggHZZ(const double sqrt_s) const;
    /**
     * @brief The ratio @f$\mu_{VBF,ZZ}@f$ between the VBF Higgs
     * production cross-section with subsequent decay into @f$Z Z^*@f$ in the
     * current model and in the Standard Model.
     * @param[in] sqrt_s the center-of-mass energy in TeV
     * @return @f$\mu_{VBF,ZZ}@f$
     */
    virtual const double muTHUVBFHZZ(const double sqrt_s) const;
    /**
     * @brief The ratio @f$\mu_{ZH,ZZ}@f$ between the ZH
     * production cross-section with subsequent decay into @f$Z Z^*@f$ in the
     * current model and in the Standard Model.
     * @param[in] sqrt_s the center-of-mass energy in TeV
     * @return @f$\mu_{ZH,ZZ}@f$
     */
    virtual const double muTHUZHZZ(const double sqrt_s) const;
    /**
     * @brief The ratio @f$\mu_{WH,ZZ}@f$ between the WH
     * production cross-section with subsequent decay into @f$Z Z^*@f$ in the
     * current model and in the Standard Model.
     * @param[in] sqrt_s the center-of-mass energy in TeV
     * @return @f$\mu_{WH,ZZ}@f$
     */
    virtual const double muTHUWHZZ(const double sqrt_s) const;
    /**
     * @brief The ratio @f$\mu_{VH,ZZ}@f$ between the VH
     * production cross-section with subsequent decay into @f$Z Z^*@f$ in the
     * current model and in the Standard Model.
     * @param[in] sqrt_s the center-of-mass energy in TeV
     * @return @f$\mu_{VH,ZZ}@f$
     */
    virtual const double muTHUVHZZ(const double sqrt_s) const;
    /**
     * @brief The ratio @f$\mu_{ttH,ZZ}@f$ between the ttH
     * production cross-section with subsequent decay into @f$Z Z^*@f$ in the
     * current model and in the Standard Model.
     * @param[in] sqrt_s the center-of-mass energy in TeV
     * @return @f$\mu_{ttH,ZZ}@f$
     */
    virtual const double muTHUttHZZ(const double sqrt_s) const;
    
    /**
     * @brief The ratio @f$\mu_{ggH,ZZ\to 4l}@f$ between the gluon-gluon fusion Higgs
     * production cross-section with subsequent decay into @f$Z Z^*\to 4l@f$ in the
     * current model and in the Standard Model.
     * @param[in] sqrt_s the center-of-mass energy in TeV
     * @return @f$\mu_{ggH,ZZ\to 4l}@f$
     */
    virtual const double muTHUggHZZ4l(const double sqrt_s) const;
    /**
     * @brief The ratio @f$\mu_{VBF,ZZ\to 4l}@f$ between the VBF Higgs
     * production cross-section with subsequent decay into @f$Z Z^*\to 4l@f$ in the
     * current model and in the Standard Model.
     * @param[in] sqrt_s the center-of-mass energy in TeV
     * @return @f$\mu_{VBF,ZZ\to 4l}@f$
     */
    virtual const double muTHUVBFHZZ4l(const double sqrt_s) const;    
    /**
     * @brief The ratio @f$\mu_{ZH,ZZ\to 4l}@f$ between the ZH
     * production cross-section with subsequent decay into @f$Z Z^*\to 4l@f$ in the
     * current model and in the Standard Model.
     * @param[in] sqrt_s the center-of-mass energy in TeV
     * @return @f$\mu_{ZH,ZZ\to 4l}@f$
     */
    virtual const double muTHUZHZZ4l(const double sqrt_s) const;
    /**
     * @brief The ratio @f$\mu_{WH,ZZ\to 4l}@f$ between the WH
     * production cross-section with subsequent decay into @f$Z Z^*\to 4l@f$ in the
     * current model and in the Standard Model.
     * @param[in] sqrt_s the center-of-mass energy in TeV
     * @return @f$\mu_{WH,ZZ\to 4l}@f$
     */
    virtual const double muTHUWHZZ4l(const double sqrt_s) const;    
    /**
     * @brief The ratio @f$\mu_{VH,ZZ\to 4l}@f$ between the VH
     * production cross-section with subsequent decay into @f$Z Z^*\to 4l@f$ in the
     * current model and in the Standard Model.
     * @param[in] sqrt_s the center-of-mass energy in TeV
     * @return @f$\mu_{VH,ZZ\to 4l}@f$
     */
    virtual const double muTHUVHZZ4l(const double sqrt_s) const;
    /**
     * @brief The ratio @f$\mu_{ttH,ZZ\to 4l}@f$ between the ttH
     * production cross-section with subsequent decay into @f$Z Z^*\to 4l@f$ in the
     * current model and in the Standard Model.
     * @param[in] sqrt_s the center-of-mass energy in TeV
     * @return @f$\mu_{ttH,ZZ\to 4l}@f$
     */
    virtual const double muTHUttHZZ4l(const double sqrt_s) const;
    
    /**
     * @brief The ratio @f$\mu_{ggH,WW}@f$ between the gluon-gluon fusion Higgs
     * production cross-section with subsequent decay into @f$W W^*@f$ in the
     * current model and in the Standard Model.
     * @param[in] sqrt_s the center-of-mass energy in TeV
     * @return @f$\mu_{ggH,WW}@f$
     */
    virtual const double muTHUggHWW(const double sqrt_s) const;
    /**
     * @brief The ratio @f$\mu_{VBF,WW}@f$ between the VBF Higgs
     * production cross-section with subsequent decay into @f$W W^*@f$ in the
     * current model and in the Standard Model.
     * @param[in] sqrt_s the center-of-mass energy in TeV
     * @return @f$\mu_{VBF,WW}@f$
     */
    virtual const double muTHUVBFHWW(const double sqrt_s) const;
    /**
     * @brief The ratio @f$\mu_{ZH,WW}@f$ between the ZH
     * production cross-section with subsequent decay into @f$W W^*@f$ in the
     * current model and in the Standard Model.
     * @param[in] sqrt_s the center-of-mass energy in TeV
     * @return @f$\mu_{ZH,WW}@f$
     */
    virtual const double muTHUZHWW(const double sqrt_s) const;
    /**
     * @brief The ratio @f$\mu_{WH,WW}@f$ between the WH
     * production cross-section with subsequent decay into @f$W W^*@f$ in the
     * current model and in the Standard Model.
     * @param[in] sqrt_s the center-of-mass energy in TeV
     * @return @f$\mu_{WH,WW}@f$
     */
    virtual const double muTHUWHWW(const double sqrt_s) const; 
    /**
     * @brief The ratio @f$\mu_{VH,WW}@f$ between the VH
     * production cross-section with subsequent decay into @f$W W^*@f$ in the
     * current model and in the Standard Model.
     * @param[in] sqrt_s the center-of-mass energy in TeV
     * @return @f$\mu_{VH,WW}@f$
     */
    virtual const double muTHUVHWW(const double sqrt_s) const;
    /**
     * @brief The ratio @f$\mu_{ttH,WW}@f$ between the ttH
     * production cross-section with subsequent decay into @f$W W^*@f$ in the
     * current model and in the Standard Model.
     * @param[in] sqrt_s the center-of-mass energy in TeV
     * @return @f$\mu_{ttH,WW}@f$
     */
    virtual const double muTHUttHWW(const double sqrt_s) const;
    /**
     * @brief The ratio @f$\mu_{ggH,WW\to 2l2\nu}@f$ between the gluon-gluon fusion Higgs
     * production cross-section with subsequent decay into @f$W W^*\to 2l2\nu@f$ in the
     * current model and in the Standard Model.
     * @param[in] sqrt_s the center-of-mass energy in TeV
     * @return @f$\mu_{ggH,WW\to 2l2\nu}@f$
     */
    virtual const double muTHUggHWW2l2v(const double sqrt_s) const;
    /**
     * @brief The ratio @f$\mu_{VBF,WW\to 2l2\nu}@f$ between the VBF Higgs
     * production cross-section with subsequent decay into @f$W W^*\to 2l2\nu@f$ in the
     * current model and in the Standard Model.
     * @param[in] sqrt_s the center-of-mass energy in TeV
     * @return @f$\mu_{VBF,WW\to 2l2\nu}@f$
     */
    virtual const double muTHUVBFHWW2l2v(const double sqrt_s) const;   
    /**
     * @brief The ratio @f$\mu_{ZH,WW\to 2l2\nu}@f$ between the ZH
     * production cross-section with subsequent decay into @f$W W^*\to 2l2\nu@f$ in the
     * current model and in the Standard Model.
     * @param[in] sqrt_s the center-of-mass energy in TeV
     * @return @f$\mu_{ZH,WW\to 2l2\nu}@f$
     */
    virtual const double muTHUZHWW2l2v(const double sqrt_s) const;
    /**
     * @brief The ratio @f$\mu_{WH,WW\to 2l2\nu}@f$ between the WH
     * production cross-section with subsequent decay into @f$W W^*\to 2l2\nu@f$ in the
     * current model and in the Standard Model.
     * @param[in] sqrt_s the center-of-mass energy in TeV
     * @return @f$\mu_{WH,WW\to 2l2\nu}@f$
     */
    virtual const double muTHUWHWW2l2v(const double sqrt_s) const;    
    /**
     * @brief The ratio @f$\mu_{VH,WW\to 2l2\nu}@f$ between the VH
     * production cross-section with subsequent decay into @f$W W^*\to 2l2\nu@f$ in the
     * current model and in the Standard Model.
     * @param[in] sqrt_s the center-of-mass energy in TeV
     * @return @f$\mu_{VH,WW\to 2l2\nu}@f$
     */
    virtual const double muTHUVHWW2l2v(const double sqrt_s) const;
    /**
     * @brief The ratio @f$\mu_{ttH,WW\to 2l2\nu}@f$ between the ttH
     * production cross-section with subsequent decay into @f$W W^*\to 2l2\nu@f$ in the
     * current model and in the Standard Model.
     * @param[in] sqrt_s the center-of-mass energy in TeV
     * @return @f$\mu_{ttH,WW\to 2l2\nu}@f$
     */
    virtual const double muTHUttHWW2l2v(const double sqrt_s) const;
    /**
     * @brief The ratio @f$\mu_{ggH,\mu\mu}@f$ between the gluon-gluon fusion Higgs
     * production cross-section with subsequent decay into @f$\mu\mu@f$ in the
     * current model and in the Standard Model.
     * @param[in] sqrt_s the center-of-mass energy in TeV
     * @return @f$\mu_{ggH,\mu\mu}@f$
     */
    virtual const double muTHUggHmumu(const double sqrt_s) const;
    /**
     * @brief The ratio @f$\mu_{VBF,\mu\mu}@f$ between the VBF Higgs
     * production cross-section with subsequent decay into @f$\mu\mu@f$ in the
     * current model and in the Standard Model.
     * @param[in] sqrt_s the center-of-mass energy in TeV
     * @return @f$\mu_{VBF,\mu\mu}@f$
     */
    virtual const double muTHUVBFHmumu(const double sqrt_s) const;    
    /**
     * @brief The ratio @f$\mu_{ZH,\mu\mu}@f$ between the ZH
     * production cross-section with subsequent decay into @f$\mu\mu@f$ in the
     * current model and in the Standard Model.
     * @param[in] sqrt_s the center-of-mass energy in TeV
     * @return @f$\mu_{ZH,\mu\mu}@f$
     */
    virtual const double muTHUZHmumu(const double sqrt_s) const;
    /**
     * @brief The ratio @f$\mu_{WH,\mu\mu}@f$ between the WH
     * production cross-section with subsequent decay into @f$\mu\mu@f$ in the
     * current model and in the Standard Model.
     * @param[in] sqrt_s the center-of-mass energy in TeV
     * @return @f$\mu_{WH,\mu\mu}@f$
     */
    virtual const double muTHUWHmumu(const double sqrt_s) const;    
    /**
     * @brief The ratio @f$\mu_{VH,\mu\mu}@f$ between the VH
     * production cross-section with subsequent decay into @f$\mu\mu@f$ in the
     * current model and in the Standard Model.
     * @param[in] sqrt_s the center-of-mass energy in TeV
     * @return @f$\mu_{VH,\mu\mu}@f$
     */
    virtual const double muTHUVHmumu(const double sqrt_s) const;
    /**
     * @brief The ratio @f$\mu_{ttH,\mu\mu}@f$ between the ttH
     * production cross-section with subsequent decay into @f$\mu\mu@f$ in the
     * current model and in the Standard Model.
     * @param[in] sqrt_s the center-of-mass energy in TeV
     * @return @f$\mu_{ttH,\mu\mu}@f$
     */
    virtual const double muTHUttHmumu(const double sqrt_s) const;
    /**
     * @brief The ratio @f$\mu_{ggH,\tau\tau}@f$ between the gluon-gluon fusion Higgs
     * production cross-section with subsequent decay into @f$\tau\tau@f$ in the
     * current model and in the Standard Model.
     * @param[in] sqrt_s the center-of-mass energy in TeV
     * @return @f$\mu_{ggH,\tau\tau}@f$
     */
    virtual const double muTHUggHtautau(const double sqrt_s) const;
    /**
     * @brief The ratio @f$\mu_{VBF,\tau\tau}@f$ between the VBF Higgs
     * production cross-section with subsequent decay into @f$\tau\tau@f$ in the
     * current model and in the Standard Model.
     * @param[in] sqrt_s the center-of-mass energy in TeV
     * @return @f$\mu_{VBF,\tau\tau}@f$
     */
    virtual const double muTHUVBFHtautau(const double sqrt_s) const;
    /**
     * @brief The ratio @f$\mu_{ZH,\tau\tau}@f$ between the ZH
     * production cross-section with subsequent decay into @f$\tau\tau@f$ in the
     * current model and in the Standard Model.
     * @param[in] sqrt_s the center-of-mass energy in TeV
     * @return @f$\mu_{ZH,\tau\tau}@f$
     */
    virtual const double muTHUZHtautau(const double sqrt_s) const;
    /**
     * @brief The ratio @f$\mu_{WH,\tau\tau}@f$ between the WH
     * production cross-section with subsequent decay into @f$\tau\tau@f$ in the
     * current model and in the Standard Model.
     * @param[in] sqrt_s the center-of-mass energy in TeV
     * @return @f$\mu_{WH,\tau\tau}@f$
     */
    virtual const double muTHUWHtautau(const double sqrt_s) const;
    /**
     * @brief The ratio @f$\mu_{VH,\tau\tau}@f$ between the VH
     * production cross-section with subsequent decay into @f$\tau\tau@f$ in the
     * current model and in the Standard Model.
     * @param[in] sqrt_s the center-of-mass energy in TeV
     * @return @f$\mu_{VH,\tau\tau}@f$
     */
    virtual const double muTHUVHtautau(const double sqrt_s) const;
    /**
     * @brief The ratio @f$\mu_{ttH,\tau\tau}@f$ between the ttH
     * production cross-section with subsequent decay into @f$\tau\tau@f$ in the
     * current model and in the Standard Model.
     * @param[in] sqrt_s the center-of-mass energy in TeV
     * @return @f$\mu_{ttH,\tau\tau}@f$
     */
    virtual const double muTHUttHtautau(const double sqrt_s) const;
    /**
     * @brief The ratio @f$\mu_{ggH,bb}@f$ between the gluon-gluon fusion Higgs
     * production cross-section with subsequent decay into @f$bb@f$ in the
     * current model and in the Standard Model.
     * @param[in] sqrt_s the center-of-mass energy in TeV
     * @return @f$\mu_{ggH,bb}@f$
     */
    virtual const double muTHUggHbb(const double sqrt_s) const;
    /**
     * @brief The ratio @f$\mu_{VBF,bb}@f$ between the VBF Higgs
     * production cross-section with subsequent decay into @f$bb@f$ in the
     * current model and in the Standard Model.
     * @param[in] sqrt_s the center-of-mass energy in TeV
     * @return @f$\mu_{VBF,bb}@f$
     */
    virtual const double muTHUVBFHbb(const double sqrt_s) const;
    /**
     * @brief The ratio @f$\mu_{ZH,bb}@f$ between the ZH
     * production cross-section with subsequent decay into @f$bb@f$ in the
     * current model and in the Standard Model.
     * @param[in] sqrt_s the center-of-mass energy in TeV
     * @return @f$\mu_{ZH,bb}@f$
     */
    virtual const double muTHUZHbb(const double sqrt_s) const;
    /**
     * @brief The ratio @f$\mu_{WH,bb}@f$ between the WH
     * production cross-section with subsequent decay into @f$bb@f$ in the
     * current model and in the Standard Model.
     * @param[in] sqrt_s the center-of-mass energy in TeV
     * @return @f$\mu_{WH,bb}@f$
     */
    virtual const double muTHUWHbb(const double sqrt_s) const;
    /**
     * @brief The ratio @f$\mu_{VH,bb}@f$ between the VH
     * production cross-section with subsequent decay into @f$bb@f$ in the
     * current model and in the Standard Model.
     * @param[in] sqrt_s the center-of-mass energy in TeV
     * @return @f$\mu_{VH,bb}@f$
     */
    virtual const double muTHUVHbb(const double sqrt_s) const;
    /**
     * @brief The ratio @f$\mu_{ttH,bb}@f$ between the ttH
     * production cross-section with subsequent decay into @f$bb@f$ in the
     * current model and in the Standard Model.
     * @param[in] sqrt_s the center-of-mass energy in TeV
     * @return @f$\mu_{ttH,bb}@f$
     */
    virtual const double muTHUttHbb(const double sqrt_s) const;
    
    /**
     * @brief The ratio @f$\mu_{VBF}@f$ between the VBF
     * production cross-section in the
     * current model and in the Standard Model, multiplied by the 
     * total (SM+new physics) invisible decay branching ratio.
     * @param[in] sqrt_s the center-of-mass energy in TeV
     * @return @f$\mu_{VBF}BR_{inv}@f$
     */    
    virtual const double muTHUVBFBRinv(const double sqrt_s) const;
    /**
     * @brief The ratio @f$\mu_{VBF,inv}@f$ between the VBF
     * production cross-section with subsequent decay into invisible states in the
     * current model and in the Standard Model.
     * @param[in] sqrt_s the center-of-mass energy in TeV
     * @return @f$\mu_{VBF,inv}@f$
     */     
    virtual const double muTHUVBFHinv(const double sqrt_s) const;
    /**
     * @brief The ratio @f$\mu_{VH}@f$ between the VH
     * production cross-section in the
     * current model and in the Standard Model, multiplied by the 
     * total (SM+new physics) invisible decay branching ratio.
     * @param[in] sqrt_s the center-of-mass energy in TeV
     * @return @f$\mu_{VH}BR_{inv}@f$
     */     
    virtual const double muTHUVHBRinv(const double sqrt_s) const;
    /**
     * @brief The ratio @f$\mu_{VH,inv}@f$ between the VH
     * production cross-section with subsequent decay into invisible states in the
     * current model and in the Standard Model.
     * @param[in] sqrt_s the center-of-mass energy in TeV
     * @return @f$\mu_{VH,inv}@f$
     */    
    virtual const double muTHUVHinv(const double sqrt_s) const;
    
    /**
     * @brief The ratio @f$\mu_{ggH,ZZ\to 4\mu}@f$ between the gluon-gluon fusion Higgs
     * production cross-section with subsequent decay into @f$Z Z^*\to 4\mu@f$ in the
     * current model and in the Standard Model.
     * @param[in] sqrt_s the center-of-mass energy in TeV
     * @return @f$\mu_{ggH,ZZ\to 4\mu}@f$
     */    
    virtual const double muTHUggHZZ4mu(const double sqrt_s) const;    
    /**
     * @brief The ratio @f$\mu_{ggH,Z\gamma\to \gamma 2\mu}@f$ between the gluon-gluon fusion Higgs
     * production cross-section with subsequent decay into @f$Z \gamma\to \gamma 2\mu@f$ in the
     * current model and in the Standard Model.
     * @param[in] sqrt_s the center-of-mass energy in TeV
     * @return @f$\mu_{ggH,Z\gamma\to \gamma 2\mu}@f$
     */    
    virtual const double muTHUggHZgamumu(const double sqrt_s) const;    

    ////////////////////////////////////////////////////////////////////////

    /**
     * @brief The new physics contribution to the coupling of the effective interaction @f$H G_{\mu\nu}^AG^{A \mu\nu}@f$.
     * @return @f$\delta g_{HGG}@f$
     */
    virtual const double deltaG_hgg() const;
    /**
     * @brief The full new physics contribution to the coupling of the effective interaction @f$H G_{\mu\nu}^AG^{A \mu\nu}@f$,
     * including new local terms and modifications on the SM-loops. Normalized to the SM value.
     * @return @f$\delta g_{HGG}/g_{HGG}^SM}@f$
     */
    virtual const double deltaG_hggRatio() const;
    /**
     * @brief The new physics contribution to the coupling of the effective interaction @f$H W_{\mu\nu}^\dagger W^{\mu\nu}@f$.
     * @return @f$\delta g_{HWW}^{(1)}@f$
     */
    virtual const double deltaG1_hWW() const;
    /**
     * @brief The new physics contribution to the coupling of the effective interaction @f$H W_{\nu}^\dagger \partial^\mu W^{\mu\nu}@f$.
     * @return @f$\delta g_{HWW}^{(2)}@f$
     */
    virtual const double deltaG2_hWW() const;
    /**
     * @brief The new physics contribution to the coupling of the effective interaction @f$H W_{\mu}^\dagger W^{\mu}@f$.
     * @return @f$\delta g_{HWW}^{(3)}@f$
     */
    virtual const double deltaG3_hWW() const;
    /**
     * @brief The new physics contribution to the coupling of the effective interaction @f$H Z_{\mu\nu} Z^{\mu\nu}@f$.
     * @return @f$\delta g_{HZZ}^{(1)}@f$
     */
    virtual const double deltaG1_hZZ() const;
    /**
     * @brief The new physics contribution to the coupling of the effective interaction @f$H Z_{\nu} \partial^\mu Z^{\mu\nu}@f$.
     * @return @f$\delta g_{HZZ}^{(2)}@f$
     */
    virtual const double deltaG2_hZZ() const;
    /**
     * @brief The new physics contribution to the coupling of the effective interaction @f$H Z_{\mu} Z^{\mu}@f$.
     * @return @f$\delta g_{HZZ}^{(3)}@f$
     */
    virtual const double deltaG3_hZZ() const;
    /**
     * @brief The new physics contribution to the coupling of the effective interaction @f$H Z_{\mu\nu} F^{\mu\nu}@f$.
     * @return @f$\delta g_{HZA}^{(1)}@f$
     */
    virtual const double deltaG1_hZA() const;
    /**
     * @brief The full new physics contribution to the coupling of the effective interaction @f$H Z_{\mu\nu} F^{A \mu\nu}@f$,
     * including new local terms and modifications on the SM-loops. Normalized to the SM value.
     * @return @f$\delta g_{HZA}^{(1)}/g_{HZA}^{(1),SM}@f$
     */
    virtual const double deltaG1_hZARatio() const;
    /**
     * @brief The new physics contribution to the coupling of the effective interaction @f$H Z_{\nu} \partial^\mu F^{\mu\nu}@f$.
     * @return @f$\delta g_{HZA}^{(2)}@f$
     */
    virtual const double deltaG2_hZA() const;
    /**
     * @brief The new physics contribution to the coupling of the effective interaction @f$H F_{\mu\nu} F^{\mu\nu}@f$.
     * @return @f$\delta g_{HAA}@f$
     */
    virtual const double deltaG_hAA() const;
    /**
     * @brief The full new physics contribution to the coupling of the effective interaction @f$H F_{\mu\nu} F^{\mu\nu}@f$,
     * including new local terms and modifications on the SM-loops. Normalized to the SM value.
     * @return @f$\delta g_{HAA}/g_{HAA}^SM}@f$
     */
    virtual const double deltaG_hAARatio() const;
    /**
     * @brief The new physics contribution to the coupling of the effective interaction @f$H f\bar{f}@f$.
     * @param[in] p a lepton or quark
     * @return @f$\delta g_{Hff}@f$
     */
    // no generation mixing
    virtual gslpp::complex deltaG_hff(const Particle p) const;   
   
    ////////////////////////////////////////////////////////////////////////
    
    /**
     * @brief The effective coupling @f$\kappa_{\mu,eff}=\sqrt{\Gamma_{H\mu\mu}/\Gamma_{H\mu\mu}^{SM}}@f$.
     * @return @f$\kappa_{\mu,eff}@f$
     */
    virtual const double kappamueff() const;
    
    /**
     * @brief The effective coupling @f$\kappa_{\tau,eff}=\sqrt{\Gamma_{H\tau\tau}/\Gamma_{H\tau\tau}^{SM}}@f$.
     * @return @f$\kappa_{\tau,eff}@f$
     */
    virtual const double kappataueff() const;
    
    /**
     * @brief The effective coupling @f$\kappa_{c,eff}=\sqrt{\Gamma_{Hcc}/\Gamma_{Hcc}^{SM}}@f$.
     * @return @f$\kappa_{c,eff}@f$
     */
    virtual const double kappaceff() const;
    
    /**
     * @brief The effective coupling @f$\kappa_{b,eff}=\sqrt{\Gamma_{Hbb}/\Gamma_{Hbb}^{SM}}@f$.
     * @return @f$\kappa_{b,eff}@f$
     */
    virtual const double kappabeff() const;
    
    /**
     * @brief The effective coupling @f$\kappa_{G,eff}=\sqrt{\Gamma_{HGG}/\Gamma_{HGG}^{SM}}@f$.
     * @return @f$\kappa_{G,eff}@f$
     */
    virtual const double kappaGeff() const;
    
    /**
     * @brief The effective coupling @f$\kappa_{Z,eff}=\sqrt{\Gamma_{HZZ}/\Gamma_{HZZ}^{SM}}@f$.
     * @return @f$\kappa_{Z,eff}@f$
     */
    virtual const double kappaZeff() const;
    
    /**
     * @brief The effective coupling @f$\kappa_{W,eff}=\sqrt{\Gamma_{HWW}/\Gamma_{HWW}^{SM}}@f$.
     * @return @f$\kappa_{W,eff}@f$
     */
    virtual const double kappaWeff() const;
    
    /**
     * @brief The effective coupling @f$\kappa_{A,eff}=\sqrt{\Gamma_{HAA}/\Gamma_{HAA}^{SM}}@f$.
     * @return @f$\kappa_{A,eff}@f$
     */
    virtual const double kappaAeff() const;
    
    /**
     * @brief The effective coupling @f$\kappa_{ZA,eff}=\sqrt{\Gamma_{HZA}/\Gamma_{HZA}^{SM}}@f$.
     * @return @f$\kappa_{ZA,eff}@f$
     */
    virtual const double kappaZAeff() const;

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
    virtual const double computeKg() const;

    /**
     * @brief A method to compute the ratio of the @f$HWW@f$ coupling in the current model and in the SM.
     * @return the ratio of the @f$HWW@f$ coupling in the current model and in the SM
     */
    virtual const double computeKW() const;

    /**
     * @brief A method to compute the ratio of the @f$HZZ@f$ coupling in the current model and in the SM.
     * @return the ratio of the @f$HZZ@f$ coupling in the current model and in the SM
     */
    virtual const double computeKZ() const;

    /**
     * @brief A method to compute the ratio of the @f$HZ\gamma@f$ coupling in the current model and in the SM.
     * @return the ratio of the @f$HZ\gamma@f$ coupling in the current model and in the SM
     */
    virtual const double computeKZga() const;

    /**
     * @brief A method to compute the ratio of the @f$H\gamma\gamma@f$ coupling in the current model and in the SM.
     * @return the ratio of the @f$H\gamma\gamma@f$ coupling in the current model and in the SM
     */
    virtual const double computeKgaga() const;
    
    /**
     * @brief A method to compute the ratio of the @f$Hee@f$ coupling in the current model and in the SM.
     * @return the ratio of the @f$Hee@f$ coupling in the current model and in the SM
     */
    virtual const double computeKe() const;
    
    /**
     * @brief A method to compute the ratio of the @f$H\mu\mu@f$ coupling in the current model and in the SM.
     * @return the ratio of the @f$H\mu\mu@f$ coupling in the current model and in the SM
     */
    virtual const double computeKmu() const;

    /**
     * @brief A method to compute the ratio of the @f$H\tau\tau@f$ coupling in the current model and in the SM.
     * @return the ratio of the @f$H\tau\tau@f$ coupling in the current model and in the SM
     */
    virtual const double computeKtau() const;
    
    /**
     * @brief A method to compute the ratio of the @f$Huu@f$ coupling in the current model and in the SM.
     * @return the ratio of the @f$Huu@f$ coupling in the current model and in the SM
     */
    virtual const double computeKu() const;

    /**
     * @brief A method to compute the ratio of the @f$Hcc@f$ coupling in the current model and in the SM.
     * @return the ratio of the @f$Hcc@f$ coupling in the current model and in the SM
     */
    virtual const double computeKc() const;

    /**
     * @brief A method to compute the ratio of the @f$Htt@f$ coupling in the current model and in the SM.
     * @return the ratio of the @f$Htt@f$ coupling in the current model and in the SM
     */
    virtual const double computeKt() const;
    
    /**
     * @brief A method to compute the ratio of the @f$Hdd@f$ coupling in the current model and in the SM.
     * @return the ratio of the @f$Hdd@f$ coupling in the current model and in the SM
     */
    virtual const double computeKd() const;
    
    /**
     * @brief A method to compute the ratio of the @f$Hss@f$ coupling in the current model and in the SM.
     * @return the ratio of the @f$Hss@f$ coupling in the current model and in the SM
     */
    virtual const double computeKs() const;

    /**
     * @brief A method to compute the ratio of the @f$Hbb@f$ coupling in the current model and in the SM.
     * @return the ratio of the @f$Hbb@f$ coupling in the current model and in the SM
     */
    virtual const double computeKb() const;
    
    /**
     * @brief A method to compute the ratio of the @f$Hgg@f$ coupling in the current model and in the SM.
     * @return the ratio of the @f$Hgg@f$ coupling in the current model and in the SM
     */
    virtual const double computeKgLoop() const;
    
    /**
     * @brief A method to compute the ratio of the @f$HZ\gamma@f$ coupling in the current model and in the SM.
     * @return the ratio of the @f$HZ\gamma@f$ coupling in the current model and in the SM
     */
    virtual const double computeKZgaLoop() const;

    /**
     * @brief A method to compute the ratio of the @f$H\gamma\gamma@f$ coupling in the current model and in the SM.
     * @return the ratio of the @f$H\gamma\gamma@f$ coupling in the current model and in the SM
     */
    virtual const double computeKgagaLoop() const;

    ////////////////////////////////////////////////////////////////////////
private:
    gslpp::complex f_func(const double x) const;
    gslpp::complex g_func(const double x) const;
    gslpp::complex Int1(const double tau, const double lambda) const;
    gslpp::complex Int2(const double tau, const double lambda) const;

    // Model parameters    
    double Kw; ///< The factor rescaling the SM-like Higgs couplings to the W vector bosons with respect to the SM.
    double Kz; ///< The factor rescaling the SM-like Higgs couplings to the Z vector bosons with respect to the SM.
    double Kg; ///< The factor rescaling the effective Higgs couplings to gluons with respect to the SM.
    double Kga; ///< The factor rescaling the effective Higgs couplings to photons with respect to the SM.
    double Kzga; ///< The factor rescaling the effective Higgs couplings to the Z and a photon with respect to the SM.
    double Ku; ///< The factor rescaling the Higgs couplings to the up quark with respect to the SM.
    double Kc; ///< The factor rescaling the Higgs couplings to the charm quark with respect to the SM.
    double Kt; ///< The factor rescaling the Higgs couplings to the top quark with respect to the SM.
    double Kd; ///< The factor rescaling the Higgs couplings to the down quark with respect to the SM.
    double Ks; ///< The factor rescaling the Higgs couplings to the strange quark with respect to the SM.
    double Kb; ///< The factor rescaling the Higgs couplings to the bottom quark with respect to the SM.
    double Ke; ///< The factor rescaling the Higgs couplings to electrons with respect to the SM.
    double Kmu; ///< The factor rescaling the Higgs couplings to muons with respect to the SM.
    double Ktau; ///< The factor rescaling the Higgs couplings to taus with respect to the SM.
    double KH; ///< The factor rescaling the total Higgs width with respect to the SM.
    double BrHinv; ///< The branching ratio of invisible Higgs decays.
    double BrHexo; ///< The branching ratio of exotic (not invisible) Higgs decays.
// The intrinsic and parametric theory relative errors in the Higgs observables.
    double eggFint; ///< Intrinsic relative theoretical error in ggF production. (Assumed to be constant in energy.)
    double eggFpar; ///< Parametric relative theoretical error in ggF production. (Assumed to be constant in energy.)
    double ettHint; ///< Intrinsic relative theoretical error in ttH production. (Assumed to be constant in energy.)
    double ettHpar; ///< Parametric relative theoretical error in ttH production. (Assumed to be constant in energy.)
    double eVBFint; ///< Intrinsic relative theoretical error in VBF production. (Assumed to be constant in energy.)
    double eVBFpar; ///< Parametric relative theoretical error in VBF production. (Assumed to be constant in energy.)
    double eWHint; ///< Intrinsic relative theoretical error in WH production. (Assumed to be constant in energy.)
    double eWHpar; ///< Parametric relative theoretical error in WH production. (Assumed to be constant in energy.)
    double eZHint; ///< Intrinsic relative theoretical error in ZH production. (Assumed to be constant in energy.)
    double eZHpar; ///< Parametric relative theoretical error in ZH production. (Assumed to be constant in energy.)
    double eeeWBFint; ///< Intrinsic relative theoretical error in \f$e^+ e^- \to H \bar{\nu} \nu\f$. (Assumed to be constant in energy.)
    double eeeWBFpar; ///< Parametric relative theoretical error in \f$e^+ e^- \to H \bar{\nu} \nu\f$. (Assumed to be constant in energy.)
    double eeeZHint; ///< Intrinsic relative theoretical error in \f$e^+ e^- \to Z H\f$. (Assumed to be constant in energy.)
    double eeeZHpar; ///< Parametric relative theoretical error in \f$e^+ e^- \to Z H\f$. (Assumed to be constant in energy.)
    double eeettHint; ///< Intrinsic relative theoretical error in \f$e^+ e^- \to t \bar{t} H\f$. (Assumed to be constant in energy.)
    double eeettHpar; ///< Parametric relative theoretical error in \f$e^+ e^- \to t \bar{t} H\f$. (Assumed to be constant in energy.)
    double eepWBFint; ///< Intrinsic relative theoretical error in \f$e^- p \to H e^- j\f$ via WBF. (Assumed to be constant in energy.)
    double eepWBFpar; ///< Parametric relative theoretical error in \f$e^- p \to H e^- j\f$ via WBF. (Assumed to be constant in energy.)    
    double eepZBFint; ///< Intrinsic relative theoretical error in \f$e^- p \to H e^- j\f$ via ZBF. (Assumed to be constant in energy.)
    double eepZBFpar; ///< Parametric relative theoretical error in \f$e^- p \to H e^- j\f$ via ZBF. (Assumed to be constant in energy.)
    double eHggint; ///< Intrinsic relative theoretical error in \f$H \to g g\f$.
    double eHggpar; ///< Parametric relative theoretical error in \f$H \to g g\f$.
    double eHWWint; ///< Intrinsic relative theoretical error in \f$H \to W W\f$.
    double eHWWpar; ///< Parametric relative theoretical error in \f$H \to W W\f$.
    double eHZZint; ///< Intrinsic relative theoretical error in \f$H \to Z Z\f$.
    double eHZZpar; ///< Parametric relative theoretical error in \f$H \to Z Z\f$.
    double eHZgaint; ///< Intrinsic relative theoretical error in \f$H \to Z \gamma\f$.
    double eHZgapar; ///< Parametric relative theoretical error in \f$H \to Z \gamma\f$.
    double eHgagaint; ///< Intrinsic relative theoretical error in \f$H \to \gamma\gamma\f$.
    double eHgagapar; ///< Parametric relative theoretical error in \f$H \to \gamma\gamma\f$.
    double eHmumuint; ///< Intrinsic relative theoretical error in \f$H \to \mu^+ \mu^-\f$.
    double eHmumupar; ///< Parametric relative theoretical error in \f$H \to \mu^+ \mu^-\f$.
    double eHtautauint; ///< Intrinsic relative theoretical error in \f$H \to \tau^+ \tau^-\f$.
    double eHtautaupar; ///< Parametric relative theoretical error in \f$H \to \tau^+ \tau^-\f$.
    double eHccint; ///< Intrinsic relative theoretical error in \f$H \to c\bar{c}\f$.
    double eHccpar; ///< Parametric relative theoretical error in \f$H \to c\bar{c}\f$.
    double eHssint; ///< Intrinsic relative theoretical error in \f$H \to s\bar{s}\f$.
    double eHsspar; ///< Parametric relative theoretical error in \f$H \to s\bar{s}\f$.
    double eHbbint; ///< Intrinsic relative theoretical error in \f$H \to b\bar{b}\f$.
    double eHbbpar; ///< Parametric relative theoretical error in \f$H \to b\bar{b}\f$.
    
    double eggFHgaga,eggFHZga,eggFHZZ,eggFHWW,eggFHtautau,eggFHbb,eggFHmumu; ///< Total relative theoretical error in \f$gg \to H \to X\f$.   
    double eVBFHgaga,eVBFHZga,eVBFHZZ,eVBFHWW,eVBFHtautau,eVBFHbb,eVBFHmumu; ///< Total relative theoretical error in \f$pp \to Hjj (VBF) \to X jj\f$.   
    double eWHgaga,eWHZga,eWHZZ,eWHWW,eWHtautau,eWHbb,eWHmumu; ///< Total relative theoretical error in \f$pp \to WH \to W X\f$.    
    double eZHgaga,eZHZga,eZHZZ,eZHWW,eZHtautau,eZHbb,eZHmumu; ///< Total relative theoretical error in \f$pp \to ZH \to Z X\f$.
    double ettHgaga,ettHZga,ettHZZ,ettHWW,ettHtautau,ettHbb,ettHmumu; ///< Total relative theoretical error in \f$pp \to ttH \to tt X\f$.
    double eVBFHinv,eVHinv; ///< Total relative theoretical error in \f$pp \to X H \to X + invisible\f$.
    
    double eHwidth; ///< Total relative theoretical error in the Higgs width
    
    // Cache variables
    double GammaHTotR;///< Total Higgs width ratio with respect to SM
    double GammaHTotSMR;///< Total Higgs width ratio with respect to SM (only SM channels)
    
    double GammaHTotnoER;///< Total Higgs width ratio with respect to SM, neglecting theory errors
    double GammaHTotSMnoER;///< Total Higgs width ratio with respect to SM, neglecting theory errors (only SM channels)
    
    // Flags    
    bool FlagKiLoop; ///< A boolean flag that is true if one allows independent kappa's for the loop induced processes (g,ga,Zga)
    bool FlagCustodial; ///< A boolean flag that is true if KZ=KW.
    bool FlagUniversalKf; ///< A boolean flag that is true if all Kf take the same universal value.
    bool FlagUniversalK; ///< A boolean flag that is true if all K take the same universal value.
    bool FlagUseKH; ///< A boolean flag that is true if the Higgs width is rescaled by KH. (See above for the description of flag.)
    
};

#endif	/* HIGGSKIGEN_H */

