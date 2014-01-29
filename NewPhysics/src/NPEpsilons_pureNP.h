/*
 * Copyright (C) 2013-2014 SusyFit Collaboration
 * All rights reserved.
 *
 * For the licensing terms see doc/COPYING.
 */

#ifndef NPEPSILONS_PURENP_H
#define	NPEPSILONS_PURENP_H

#include "NPbase.h"

/**
 * @class NPEpsilons_pureNP
 * @brief A class for new physics in the form of contributions to the \f$\varepsilon_{1,2,3,b}\f$ parameters.
 * 
 * @ingroup NewPhysics
 * @author SusyFit Collaboration
 * @copyright GNU General Public License
 * @details This class contains the necessary functions to work with new physics 
 * corrections to electroweak precision observables, in the form of contributions
 * to the \f$\varepsilon_{1,2,3,b}\f$ parameters \cite. 
 * Only new physics contributions to \f$\varepsilon_i\f$ are parameterized.
 */
class NPEpsilons_pureNP : public NPbase {
public:
    /**
     * @brief The number of new physics parameters in the model.
     */
    static const int NEPSILONpureNPvars = 4;
    /**
     * @brief A string array with the names of the new physics parameters in the model.
     */
    static const std::string EPSILONpureNPvars[NEPSILONpureNPvars];

    /**
     * @brief The default constructor. 
     */
    NPEpsilons_pureNP();

    /**
     * @brief The name of the model.
     * @return the name of the model as a string
     */
    virtual std::string ModelName() const
    {
        return "NPEpsilons_pureNP";
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
     * @return the pure new physics contribution to \f$\varepsilon_1\f$
     */
    virtual double epsilon1() const;

    /**
     * @return the pure new physics contribution to \f$\varepsilon_2\f$
     */
    virtual double epsilon2() const;
    
    /**
     * @return the pure new physics contribution to \f$\varepsilon_3\f$
     */
    virtual double epsilon3() const;
    
    /**
     * @return the pure new physics contribution to \f$\varepsilon_b\f$
     */
    virtual double epsilonb() const;


    ////////////////////////////////////////////////////////////////////////

    /**
     * @return the \f$W\f$-boson mass in GeV
     */
    virtual double Mw() const;

    /**
     * @return the (square of the) cosine of the weak angle in the On-mass-shell renormalization scheme,
     *  \f$\cos^2{\theta_W}=\frac{M_W^2}{M_Z^2}\f$
     */
    virtual double cW2() const;

    /**
     * @return the (square of the) sine of the weak angle in the On-mass-shell renormalization scheme,
     *  \f$\sin^2{\theta_W}=1-\frac{M_W^2}{M_Z^2}\f$
     */
    virtual double sW2() const;

    /**
     * @return the total width of the \f$W\f$ boson in GeV [NOT IMPLEMENTED YET]
     */
    virtual double GammaW() const;


    ////////////////////////////////////////////////////////////////////////

  /**
     * @param[in] l name of a lepton
     * @return the new physics correction to the neutral-current vector coupling @f$g_A^l@f$
     */
    virtual double deltaGVl(StandardModel::lepton l) const;

   /**
     * @param[in] q name of a quark
     * @return the new physics correction to the neutral-current vector coupling @f$g_V^q@f$
     */  
    virtual double deltaGVq(StandardModel::quark q) const;

   /**
     * @param[in] l name of a lepton
     * @return the new physics correction to the neutral-current axial-vector coupling @f$g_A^l@f$
     */  
    virtual double deltaGAl(StandardModel::lepton l) const;

   /**
     * @param[in] q name of a quark
     * @return the new physics correction to the neutral-current axial-vector coupling @f$g_A^q@f$
     */
    virtual double deltaGAq(StandardModel::quark q) const;


    ////////////////////////////////////////////////////////////////////////
protected:
    double deltaEps_1, deltaEps_2, deltaEps_3, deltaEps_b;
    /**
     * @brief A set method to fix the parameters of the model.
     * @param[in] name a string with the parameter name
     * @param[in] value the value to be asigned to the parameter specified by name
     */
    virtual void setParameter(const std::string name, const double& value);

};

#endif	/* NPEPSILONS_PURENP_H */

