/*
 * Copyright (C) 2014 SusyFit Collaboration
 * All rights reserved.
 *
 * For the licensing terms see doc/COPYING.
 */

#ifndef HIGGSKVKF_H
#define	HIGGSKVKF_H
#include "HiggsExtensionModel.h"
/**
 * @class HiggsKvKf
 * @ingroup HiggsExtensions
 * @brief A model class extending the %StandardModel Higgs sector with two universal couplings.
 * @author SusyFit Collaboration
 * @copyright GNU General Public License
 * @details This is a Model class containing parameters and functions associated
 * with an extension of the %StandardModel where Higgs couplings to all vector bosons
 * are rescaled by @f[K_v@f] and Higgs couplings to all fermions are rescaled by @f[K_f@f]. 
 * This class inherits from the %HiggsExtensionModel class, which defines parameters related to generic
 * extensions of the %StandardModel Higgs sector.
 *
 *
 * @anchor HiggsKvKfInitialization
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
 * @anchor HiggsKvKfParameters
 * <h3>%Model parameters</h3>
 *
 * The model parameters of %HiggsKvKf (in addition to the %StandardModel ones) are summarized below:
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
 *   <td class="mod_name">%Kf</td>
 *   <td class="mod_symb">@f$\kappa_f@f$</td>
 *   <td class="mod_desc">The factor rescaling all Higgs couplings to fermions with respect to the SM.</td>
 * </tr>
 * </table>
 * 
 * Please read information about parameter initialization and update in the documentation of the %StandardModel class.
 */ 
class HiggsKvKf : public HiggsExtensionModel {
public:

    static const int NHKvKfvars = 2;///< The number of the model parameters in %HiggsKvKf.

    /**
     * @brief  A string array containing the labels of the model parameters in %HiggsKvKf.
     */
    static const std::string HKvKfvars[NHKvKfvars];

    HiggsKvKf() : HiggsExtensionModel() 
    {
    };
    HiggsKvKf(const HiggsKvKf& orig);
    virtual ~HiggsKvKf() {};
    
    virtual bool InitializeModel();
    
     ///////////////////////////////////////////////////////////////////////////
    // Model parameters

    /**
     * @brief A method to check if all the mandatory parameters for %HiggsKvKf
     * have been provided in model initialization.
     * @param[in] DPars a map of the parameters that are being updated in the Monte Carlo run
     * (including parameters that are varied and those that are held constant)
     * @return a boolean that is true if the execution is successful
     */
    virtual bool CheckParameters(const std::map<std::string, double>& DPars);


   virtual double computeKW() const
    {
        return Kv;
    }
    virtual double computeKZ() const
    {
        return Kv;
    }

    /**
     * @brief A method to compute the ratio of the @f[HZ\gamma@f] coupling in the current model and in the SM.
     * @return the ratio of the @f[HZ\gamma@f] coupling in the current model and in the SM
     */
    virtual double computeKZga() const
    {
        double gtt = computeGammaZgatt();
        double gWW = computeGammaZgaWW();
        double gtW = computeGammaZgatW();
        return sqrt((gtt*Kf*Kf + gWW*Kv*Kv + gtW*Kf*Kv)/(gtt + gWW + gtW));
    }

     /**
     * @brief A method to compute the ratio of the @f[H\gamma\gamma@f] coupling in the current model and in the SM.
     * @return the ratio of the @f[H\gamma\gamma@f] coupling in the current model and in the SM
     */
    virtual double computeKgaga() const
    {
        double gtt = computeGammagagatt();
        double gWW = computeGammagagaWW();
        double gtW = computeGammagagatW();
        return sqrt((gtt*Kf*Kf + gWW*Kv*Kv + gtW*Kf*Kv)/(gtt + gWW + gtW));
    }

   virtual double computeKb() const
    {
        return Kf;
    }

   virtual double computeKc() const
    {
        return Kf;
    }

    virtual double computeKglgl() const
    {
        return Kf;
    }

    virtual double computeKt() const
    {
        return Kf;
    }

    virtual double computeKtau() const
    {
        return Kf;
    }
    /**
     * @brief This method computes the ratio of the total Higgs width w.r.t SM.
     * @return The he ratio of the total Higgs width w.r.t SM
     */
    double computeGTotalRatio() const
    {
        return computeKW()*computeKW()*computeBRWW()+
               computeKZ()*computeKZ()*computeBRZZ()+
               computeKgaga()*computeKgaga()*computeBRgaga()+
               computeKglgl()*computeKglgl()*computeBRglgl()+
               computeKb()*computeKb()*computeBRbb()+
               computeKc()*computeKc()*computeBRcc()+
               computeKtau()*computeKtau()*computeBRtautau();
    }
    
    double getKf() const {
        return Kf;
    }

    void setKf(double Kf) {
        this->Kf = Kf;
    }

    double getKv() const {
        return Kv;
    }

    void setKv(double Kv) {
        this->Kv = Kv;
    }

    protected:
        
     /**
     * @brief A method to set the value of a parameter of %HiggsKvKf.
     * @param[in] name name of a model parameter
     * @param[in] value the value to be assigned to the parameter specified by name
     */
    virtual void setParameter(const std::string name, const double& value);

    private:
        double Kv, Kf;
};

#endif	/* HIGGSKVKF_H */

