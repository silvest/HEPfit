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
     * @brief The oblique paramenter \f$S\f$.
     * @return the value of the oblique parameter @f$S@f$
     */
    virtual double obliqueS() const;

    /**
     * @brief The oblique paramenter \f$T\f$.
     * @return the value of the oblique parameter @f$T@f$
     */
    virtual double obliqueT() const;

    /**
     * @brief The oblique paramenter \f$U\f$.
     * @return the value of the oblique parameter @f$U@f$
     */
    virtual double obliqueU() const;

    /**
     * @brief The oblique paramenter \f$\hat{S}\f$.
     * @return the value of the oblique parameter \f$\hat{S}=\frac{\alpha}{4\sin^2{\theta_W}}S\f$
     */
    virtual double obliqueShat() const;

    /**
     * @brief The oblique paramenter \f$\hat{T}\f$.
     * @return the value of the oblique parameter \f$\hat{T}=\alpha T\f$
     */
    virtual double obliqueThat() const;

    /**
     * @brief The oblique paramenter \f$\hat{U}\f$.
     * @return the value of the oblique parameter \f$\hat{U}=-\frac{\alpha}{4\sin^2{\theta_W}}U\f$
     */
    virtual double obliqueUhat() const;

    /**
     * @brief The oblique paramenter \f$V\f$.
     * @return the value of the oblique parameter \f$V\f$
     */
    virtual double obliqueV() const;

    /**
     * @brief The oblique paramenter \f$W\f$.
     * @return the value of the oblique parameter \f$W\f$
     */
    virtual double obliqueW() const;

    /**
     * @brief The oblique paramenter \f$X\f$.
     * @return the value of the oblique parameter \f$X\f$
     */
    virtual double obliqueX() const;

    /**
     * @brief The oblique paramenter \f$Y\f$.
     * @return the value of the oblique parameter \f$Y\f$
     */
    virtual double obliqueY() const;
    
    /**
     * @brief The new physics corrections to the Fermi constant.
     * the new physics correction to the Fermi constant, 
     * \f$ G_F\equiv G_F^\mathrm{SM}(1+\Delta G_F)\f$
     */
    virtual double DeltaGF() const;

    /**
     * @brief The new physics correction to @f$g_V^l@f$.
     * @param[in] l name of a lepton
     * @return the new physics contribution to neutral-current vector coupling @f$\delta g_V^l@f$
     */
    virtual double deltaGVl(StandardModel::lepton l) const;
    
    /**
     * @brief The new physics correction to @f$g_V^q@f$.
     * @param[in] q name of a quark
     * @return the new physics contribution to neutral-current vector coupling @f$\delta g_V^q@f$
     */
    virtual double deltaGVq(StandardModel::quark q) const;

    /**
     * @brief The new physics correction to @f$g_A^l@f$.
     * @param[in] l name of a lepton
     * @return the new physics contribution to neutral-current axial-vector coupling @f$\delta g_A^l@f$
     */   
    virtual double deltaGAl(StandardModel::lepton l) const;

    /**
     * @brief The new physics correction to @f$g_A^q@f$.
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

