/*
 * Copyright (C) 2014 HEPfit Collaboration
 * All rights reserved.
 *
 * For the licensing terms see doc/COPYING.
 */

#ifndef HIGGSKVGENKF_H
#define	HIGGSKVGENKF_H
#include <NPbase.h>

/**
 * @class HiggsKvgenKf
 * @ingroup HiggsExtensions
 * @brief A model class extending the %StandardModel Higgs sector with three couplings.
 * @author HEPfit Collaboration
 * @copyright GNU General Public License
 * @details This is a Model class containing parameters and functions associated
 * with an extension of the %StandardModel where Higgs couplings to the @f$Z@f$ boson
 * are rescaled by @f$K_Z@f$, Higgs couplings to the @f$W@f$ boson
 * are rescaled by @f$K_W@f$ and Higgs couplings to all fermions are rescaled by @f$K_f@f$.
 * This class inherits from the %HiggsExtensionModel class, which defines parameters related to generic
 * extensions of the %StandardModel Higgs sector.
 *
 *
 * @anchor HiggsKvgenKfInitialization
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
 * @anchor HiggsKvgenKfParameters
 * <h3>%Model parameters</h3>
 *
 * The model parameters of %HiggsKvgenKf (in addition to the %StandardModel ones) are summarized below:
 * <table class="model">
 * <tr>
 *   <th>Label</th>
 *   <th>LaTeX symbol</th>
 *   <th>Description</th>
 * </tr>
 * <tr>
 *   <td class="mod_name">%KZ</td>
 *   <td class="mod_symb">@f$\kappa_Z@f$</td>
 *   <td class="mod_desc">The factor rescaling Higgs couplings to @f$Z@f$ bosons with respect to the SM.</td>
 * </tr>
 * <tr>
 *   <td class="mod_name">%KW</td>
 *   <td class="mod_symb">@f$\kappa_W@f$</td>
 *   <td class="mod_desc">The factor rescaling Higgs couplings to @f$W@f$ bosons with respect to the SM.</td>
 * </tr>
 * <tr>
 *   <td class="mod_name">%Kf</td>
 *   <td class="mod_symb">@f$\kappa_f@f$</td>
 *   <td class="mod_desc">The factor rescaling all Higgs couplings to fermions with respect to the SM.</td>
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
class HiggsKvgenKf : public NPbase {
public:

    static const int NHKvgenKfvars = 4; ///< The number of the model parameters in %HiggsKvgenKf.

    /**
     * @brief A string array containing the labels of the model parameters in %HiggsKvgenKf.
     */
    static const std::string HKvgenKfvars[NHKvgenKfvars];

    HiggsKvgenKf();

    virtual ~HiggsKvgenKf()
    {
    };

    double getKf() const
    {
        return Kf;
    }

    void setKf(double Kf)
    {
        this->Kf = Kf;
    }

    double getKW() const
    {
        return KW;
    }

    void setKW(double KW)
    {
        this->KW = KW;
    }

    double getKZ() const
    {
        return KZ;
    }

    void setKZ(double KZ)
    {
        this->KZ = KZ;
    }

    /**
     * @brief A method to check if all the mandatory parameters for %HiggsKvgenKf
     * have been provided in model initialization.
     * @param[in] DPars a map of the parameters that are being updated in the Monte Carlo run
     * (including parameters that are varied and those that are held constant)
     * @return a boolean that is true if the execution is successful
     */
    virtual bool CheckParameters(const std::map<std::string, double>& DPars);

    ////////////////////////////////////////////////////////////////////////

    virtual double muggH(const double sqrt_s) const;
    virtual double muVBF(const double sqrt_s) const;
    virtual double muWH(const double sqrt_s) const;
    virtual double muZH(const double sqrt_s) const;
    virtual double muVH(const double sqrt_s) const;
    virtual double muVBFpVH(const double sqrt_s) const;
    virtual double muttH(const double sqrt_s) const;
    virtual double muggHpttH(const double sqrt_s) const;
    virtual double BrHggRatio() const;
    virtual double BrHWWRatio() const;
    virtual double BrHZZRatio() const;
    virtual double BrHZgaRatio() const;
    virtual double BrHgagaRatio() const;
    virtual double BrHtautauRatio() const;
    virtual double BrHccRatio() const;
    virtual double BrHbbRatio() const;
    virtual double computeGammaTotalRatio() const;

    ////////////////////////////////////////////////////////////////////////
protected:

    /**
     * @brief A method to set the value of a parameter of %HiggsKvKf.
     * @param[in] name name of a model parameter
     * @param[in] value the value to be assigned to the parameter specified by name
     */
    virtual void setParameter(const std::string name, const double& value);

    virtual double computeKg() const;

    virtual double computeKW() const;

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

    virtual double computeKtau() const;

    virtual double computeKc() const;

    virtual double computeKt() const;

    virtual double computeKb() const;

    ////////////////////////////////////////////////////////////////////////
private:
    double KW, KZ, Kf, BrHinv;
};

#endif	/* HIGGSKVGENKF_H */

