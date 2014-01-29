/*
 * Copyright (C) 2013-2014 SusyFit Collaboration
 * All rights reserved.
 *
 * For the licensing terms see doc/COPYING.
 */

#ifndef NPBASE_H
#define	NPBASE_H

#include <StandardModel.h>

/**
 * @class NPbase
 * @brief A base class for general new physics corrections to electroweak precision
 * observables
 * @ingroup NewPhysics
 * @author SusyFit Collaboration
 * @copyright GNU General Public License
 * @details The methods in this class contain the basic structure to compute
 * new physics corrections to the neutral-current couplings to fermions in the other
 * \b NewPhysics classes. In several of these classes some of these methods are reimplemented
 * to account for the details of more specific scenarios. 
 * 
 * Flags:
 * \li \b FlagFixSMcontribution:&nbsp; (True) (False).
 */
class NPbase : public StandardModel {
public:
    /**
     * @brief The number of flags in this class.
     */
    static const int NNPbaseflags = 1;
    /**
     * @brief A string array with the names of the flags in the class.
     */
    static const std::string NPbaseflags[NNPbaseflags];

    /**
     * @brief Constructor. 
     */
    NPbase();

    /**
     * @brief The name of the model.
     * @return the name of the model as a string
     */
    virtual std::string ModelName() const
    {
        return "NPbase";
    }

    /**
     * @brief A method to initialize the model.
     * @return true is model initialization is successful
     */
    virtual bool InitializeModel();
    
    /**
     * @brief A set method to fix the flags EWSM for the model.
     * @param[in] myEWSM flag value for the EWSM flag
     * @return a boolean to designate the success or failure of this procedure
     */
    virtual void setEWSMflags(EWSM& myEWSM);

    /**
     * @brief A method to initialize the model.
     * @param[in] Dpars a map of parameters that are being updated in the Monte Carlo run
     * @return true is model initialization is successful
     */
    virtual bool Init(const std::map<std::string, double>& DPars);
    
    /**
     * @brief The update method for the model class.
     * @details This method updates all the parameters of the model every time a
     * new set of parameters is generated.
     * @param[in] Dpars a map of parameters that are being updated in the Monte Carlo run
     * (including parameters that are varied and those that are held constant)
     * @return a boolean that is true if the execution is successful.
     */
    virtual bool Update(const std::map<std::string, double>& DPars);
    
    /**
     * @brief A method to check if all the mandatory parameters for the model have been
     * provided in the model configuration file.
     * @param[in] Dpars a map of parameters that are being updated in the Monte Carlo run
     * (including parameters that are varied and those that are held constant)
     */
    virtual bool CheckParameters(const std::map<std::string, double>& DPars);

    /**
     * @brief A set method to fix the flags for the model.
     * @param[in] flagname the name of the flag
     * @param[in] flagvalue the value of the flag that can be true or false
     * @return a boolean to designate the success or failure of this procedure
     */
    virtual bool setFlag(const std::string, const bool&);
    
    /**
     * @brief A method to check the sanity of the set of flags.
     * @return true if the set of flags is sane.
     */
    virtual bool CheckFlags() const;

    
    ////////////////////////////////////////////////////////////////////////

    /**
     * @brief A method to check the value of the flag \b FlagFixSMcontribution
     * @return the value of \b FlagFixSMcontribution
     */
    bool IsFlagFixSMcontribution() const
    {
        return FlagFixSMcontribution;
    }


    ////////////////////////////////////////////////////////////////////////

    /**
     * @return the oblique parameter @f$S@f$
     */
    virtual double obliqueS() const;

    /**
     * @return the oblique parameter @f$T@f$
     */
    virtual double obliqueT() const;

    /**
     * @return the oblique parameter @f$U@f$
     */
    virtual double obliqueU() const;

    /**
     * @return the oblique parameter \f$\hat{S}=\frac{\alpha}{4\sin^2{\theta_W}}S\f$
     */
    virtual double obliqueShat() const;

    /**
     * @return the oblique parameter \f$\hat{T}=\alpha T\f$
     */
    virtual double obliqueThat() const;

    /**
     * @return the oblique parameter \f$\hat{U}=-\frac{\alpha}{4\sin^2{\theta_W}}U\f$
     */
    virtual double obliqueUhat() const;

    /**
     * @return the oblique parameter \f$V\f$
     */
    virtual double obliqueV() const;

    /**
     * @return the oblique parameter \f$W\f$
     */
    virtual double obliqueW() const;

    /**
     * @return the oblique parameter \f$X\f$
     */
    virtual double obliqueX() const;

    /**
     * @return the oblique parameter \f$Y\f$
     */
    virtual double obliqueY() const;
    
    /**
     * @return the new physics contribution to the Fermi constant, @f$G_{\mu}=G_{\mu}^{\rm SM}(1+\Delta GF)@f$
     */
    virtual double DeltaGF() const;

    /**
     * @param[in] l name of a lepton
     * @return the new physics contribution to neutral-current vector coupling @f$\delta g_V^l@f$
     */
    virtual double deltaGVl(StandardModel::lepton l) const;
    
    /**
     * @param[in] q name of a quark
     * @return the new physics contribution to neutral-current vector coupling @f$\delta g_V^q@f$
     */
    virtual double deltaGVq(StandardModel::quark q) const;

    /**
     * @param[in] l name of a lepton
     * @return the new physics contribution to neutral-current axial-vector coupling @f$\delta g_A^l@f$
     */   
    virtual double deltaGAl(StandardModel::lepton l) const;

    /**
     * @param[in] q name of a quark
     * @return the new physics contribution to neutral-current axial-vector coupling @f$\delta g_A^q@f$
     */ 
    virtual double deltaGAq(StandardModel::quark q) const;

    
    ////////////////////////////////////////////////////////////////////////
protected:
     /**
     * @brief A set method to fix the parameters of the model.
     * @param[in] name a string with the parameter name
     * @param[in] value the value to be asigned to the parameter specified by name
     */
    virtual void setParameter(const std::string name, const double& value);


    ////////////////////////////////////////////////////////////////////////
private:
    bool FlagFixSMcontribution;
    
};

#endif	/* NPBASE_H */

