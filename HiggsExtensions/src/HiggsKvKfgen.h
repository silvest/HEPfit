/*
 * Copyright (C) 2014 HEPfit Collaboration
 * All rights reserved.
 *
 * For the licensing terms see doc/COPYING.
 */

#ifndef HIGGSKVKFGEN_H
#define	HIGGSKVKFGEN_H
#include "NPbase.h"

/**
 * @class HiggsKvKfgen
 * @ingroup HiggsExtensions
 * @brief A model class extending the %StandardModel Higgs sector with five flavour-universal couplings.
 * @author HEPfit Collaboration
 * @copyright GNU General Public License
 * @details This is a Model class containing parameters and functions associated
 * with an extension of the %StandardModel where Higgs couplings to all vector bosons
 * are rescaled by @f$\kappa_v@f$ and Higgs couplings to all up, down and lepton fermions
 * are rescaled by @f$\kappa_u, \kappa_d, \kappa_l@f$ respectively.
 * The invisible decay width is also parametrized independently by Br@f$(H\to invisible)@f$.
 * This class inherits from the %NPbase class, which defines parameters related to generic
 * extensions of the %StandardModel Higgs sector.
 *
 *
 * @anchor HiggsKvKfgenInitialization
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
 * @anchor HiggsKvKfgenParameters
 * <h3>%Model parameters</h3>
 *
 * The model parameters of %HiggsKvKfgen (in addition to the %StandardModel ones) are summarized below:
 * <table class="model">
 * <tr>
 *   <th>Label</th>
 *   <th>LaTeX symbol</th>
 *   <th>Description</th>
 * </tr>
 * <tr>
 *   <td class="mod_name">%Kv</td>
 *   <td class="mod_symb">@f$\kappa_V@f$</td>
 *   <td class="mod_desc">The factor rescaling all Higgs couplings to vector bosons with respect to the SM.</td>
 * </tr>
 * <tr>
 *   <td class="mod_name">%Ku</td>
 *   <td class="mod_symb">@f$\kappa_u@f$</td>
 *   <td class="mod_desc">The factor rescaling all Higgs couplings to up-type fermions with respect to the SM.</td>
 * </tr>
 * <tr>
 *   <td class="mod_name">%Kd</td>
 *   <td class="mod_symb">@f$\kappa_d@f$</td>
 *   <td class="mod_desc">The factor rescaling all Higgs couplings to down-type fermions with respect to the SM.</td>
 * </tr>
 * <tr>
 *   <td class="mod_name">%Kl</td>
 *   <td class="mod_symb">@f$\kappa_l@f$</td>
 *   <td class="mod_desc">The factor rescaling all Higgs couplings to leptons with respect to the SM.</td>
 * </tr>
 * <tr>
 *   <td class="mod_name">%BrHinv</td>
 *   <td class="mod_symb">Br@f$(H\to invisible)@f$</td>
 *   <td class="mod_desc">The branching ratio of invisible Higgs decays.</td>
 * </tr>
 * </table>
 * 
 * Please read information about parameter initialization and update in the documentation of the %StandardModel class.
 */
class HiggsKvKfgen : public NPbase {
public:

    static const int NHKvKfgenvars = 5; ///< The number of the model parameters in %HiggsKvKf.

    /**
     * @brief A string array containing the labels of the model parameters in %HiggsKvKf.
     */
    static const std::string HKvKfgenvars[NHKvKfgenvars];

    /**
     * @brief The default constructor.
     */
    HiggsKvKfgen();

    /**
     * @brief The default destructor.
     */
    virtual ~HiggsKvKfgen()
    {
    };

    /**
     * @brief A get method to retrieve the factor rescaling the Higgs coupling
     * to EW vector bosons with respect to the SM @f$K_V@f$.
     * @return @f$K_V@f$
     */
    double getKv() const
    {
        return Kv;
    }

    /**
     * @brief A set method to change the factor rescaling the Higgs coupling
     * to EW vector bosons with respect to the SM @f$K_V@f$.
     * @param[in] @f$K_V@f$ the factor rescaling the Higgs coupling to EW vector bosons.
     */
    void setKv(double Kv)
    {
        this->Kv = Kv;
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
     * to charged leptons with respect to the SM @f$K_l@f$.
     * @return @f$K_l@f$
     */
    double getKl() const
    {
        return Kl;
    }

    /**
     * @brief A set method to change the factor rescaling the Higgs coupling
     * to charged leptons with respect to the SM @f$K_l@f$.
     * @param[in] @f$K_f@f$ the factor rescaling the Higgs coupling to charged leptons.
     */
    void setKl(double Kl)
    {
        this->Kl = Kl;
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
     * @brief A method to check if all the mandatory parameters for %HiggsKvKfgen
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
     * @brief The ratio of the @f$\Gamma(H)@f$ in the current model
     * and in the Standard Model.
     * @return @f$\Gamma(H)@f$/@f$\Gamma(H)_{\mathrm{SM}}@f$
     */
    virtual double computeGammaTotalRatio() const;

    ////////////////////////////////////////////////////////////////////////
protected:

    /**
     * @brief A method to set the value of a parameter of %HiggsKvKf.
     * @param[in] name name of a model parameter
     * @param[in] value the value to be assigned to the parameter specified by name
     */
    virtual void setParameter(const std::string name, const double& value);

    /**
     * @brief A method to compute the ratio of the @f$Hgg@f$ coupling in the current model and in the SM.
     * @return the ratio of the @f$Hgg@f$ coupling in the current model and in the SM
     */
    virtual double computeKg(const double sqrt_s) const;

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
    double Kv; ///< The factor rescaling all Higgs couplings to vector bosons with respect to the SM. 
    double Ku; ///< The factor rescaling all Higgs couplings to up-type fermions with respect to the SM.
    double Kd; ///< The factor rescaling all Higgs couplings to down-type fermions with respect to the SM.
    double Kl; ///< The factor rescaling all Higgs couplings to leptons with respect to the SM.
    double BrHinv; ///< The branching ratio of invisible Higgs decays.
};

#endif	/* HIGGSKVKFGEN_H */

