/* 
 * Copyright (C) 2012-2014 SusyFit Collaboration
 * All rights reserved.
 *
 * For the licensing terms see doc/COPYING.
 */

#ifndef NPSTU_H
#define	NPSTU_H

#include <EW_BURGESS.h>
#include "NPbase.h"

/**
 * @class NPSTU
 * @brief A class for new physics in the form of contributions to the oblique 
 * parameters \f$S,~T\f$ and \f$U\f$. 
 * @ingroup NewPhysics
 * @author SusyFit Collaboration
 * @copyright GNU General Public License
 * @details This class contains the necessary functions to compute new physics 
 * tree-level corrections to electroweak precision observables, in the form of 
 * contributions to the Peskin-Takeuchi oblique parameters \cite Peskin:1991sw. These corrections
 * are parameterized in terms of the \f$S,~T\f$ and \f$U\f$ contributions to \f$M_W\f$,
 * and to \f$Z\f$-pole observables through the corrections to the different neutral-current
 * effective couplings to leptons and quarks. The contributions to the later are 
 * implemented in the \b NPbase class.
 */
class NPSTU : public NPbase {
public:
    /**
     * @brief The number of new physics parameters in the model.
     */
    static const int NSTUvars = 3;
    /**
     * @brief A string array with the names of the new physics parameters in the model.
     */
    static const std::string STUvars[NSTUvars];
    
    /**
     * @brief Constructor.
     */
    NPSTU();

    /**
     * @brief The name of the model.
     * @return the name of the model as a string
     */
    virtual std::string ModelName() const 
    {
        return "NPSTU";
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
     * @brief The oblique parameter \f$S\f$.
     * @return the value of the oblique parameter \f$S\f$
     */
    virtual double obliqueS() const 
    {
        return myObliqueS;
    }

    /**
     * @brief The oblique parameter \f$T\f$.
     * @return the value of the oblique parameter \f$T\f$
     */
    virtual double obliqueT() const 
    {
        return myObliqueT;
    }

    /**
     * @brief The oblique parameter \f$U\f$.
     * @return the value of the oblique parameter \f$U\f$
     */
    virtual double obliqueU() const 
    {
        return myObliqueU;
    }


    ////////////////////////////////////////////////////////////////////////

    /**
     * @brief The parameter \f$\varepsilon_1\f$.
     * @return the value of the @f$\epsilon_1@f$ parameter (SM plus new physics
     * corrections)
     */
    double epsilon1() const;

    /**
     * @brief The parameter \f$\varepsilon_2\f$.
     * @return the value of the @f$\epsilon_2@f$ parameter (SM plus new physics
     * corrections)
     */
    double epsilon2() const;

    /**
     * @brief The parameter \f$\varepsilon_3\f$.
     * @return the value of the @f$\epsilon_3@f$ parameter (SM plus new physics
     * corrections)
     */
    double epsilon3() const;

    /**
     * @brief The parameter \f$\varepsilon_b\f$.
     * @return the SM value of the @f$\epsilon_b@f$ parameter
     */
    double epsilonb() const;

    
    ////////////////////////////////////////////////////////////////////////

    /**
     * @brief The \f$W\f$ boson mass.
     * @return the \f$W\f$-boson mass in GeV
     */
    virtual double Mw() const;

    /**
     * @brief The (square of the) cosine of the weak angle \f$\cos^2{\theta_W}\f$.
     * @return the value of \f$\cos^2{\theta_W}\f$ in the On-mass-shell renormalization scheme,
     *  \f$\cos^2{\theta_W}=\frac{M_W^2}{M_Z^2}\f$
     */
    virtual double cW2() const;

    /**
     * @brief The (square of the) sine of the weak angle \f$\sin^2{\theta_W}\f$.
     * @return the value of \f$\sin^2{\theta_W}\f$ in the On-mass-shell renormalization scheme,
     *  \f$\sin^2{\theta_W}=1-\frac{M_W^2}{M_Z^2}\f$
     */
    virtual double sW2() const;

    /**
     * @brief The \f$W\f$ decay width \f$\Gamma_W\f$.
     * @return the total width of the \f$W\f$ boson in GeV
     */
    virtual double GammaW() const;

    
    ////////////////////////////////////////////////////////////////////////
protected:    
    double myObliqueS, myObliqueT, myObliqueU;
     /**
     * @brief A set method to fix the parameters of the model.
     * @param[in] name a string with the parameter name
     * @param[in] value the value to be asigned to the parameter specified by name
     */
    virtual void setParameter(const std::string name, const double& value);

};

#endif	/* NPSTU_H */

