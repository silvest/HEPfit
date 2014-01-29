/* 
 * Copyright (C) 2013-2014 SusyFit Collaboration
 * All rights reserved.
 *
 * For the licensing terms see doc/COPYING.
 */

#ifndef NPHIGGSST_H
#define	NPHIGGSST_H

#include "NPbase.h"

/**
 * @class NPHiggsST
 * @brief A class for new physics with a non-standard @f$HVV@f$ couplings. 
 * @ingroup NewPhysics
 * @author SusyFit Collaboration
 * @copyright GNU General Public License
 * @details This class is based in new physics models described by the EW chiral
 * Lagrangian plus all interactions involving a light neutral scalar \f$H\f$. The
 * effective Lagrangian expansion is truncated at the level of two-derivatives \cite.
 * Currently, only the dominant effects in EWPD, coming from non-standard @f$HVV@f$
 * interactions, are considered.
 */
class NPHiggsST : public NPbase {
public:
    /**
     * @brief The number of new physics parameters in the model.
     */
    static const int NNPHIGGSSTvars = 8;
    /**
     * @brief A string array with the names of the new physics parameters in the model.
     */
    static const std::string NPHIGGSSTvars[NNPHIGGSSTvars];

    /**
     * @brief The default constructor.
     */
    NPHiggsST();

    /**
     * @brief The name of the model.
     * @return the name of the model as a string
     */
    virtual std::string ModelName() const 
    {
        return "NPHiggsST";
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
     * @return the contribution to the oblique parameter \f$S\f$ in this model
     */
    virtual double obliqueS() const;
        
    /**
     * @return the contribution to the oblique parameter \f$T\f$ in this model
     */
    virtual double obliqueT() const;
    
    /**
     * @return the contribution to the oblique parameter \f$U\f$ in this model 
     */
    virtual double obliqueU() const;


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
     * @return the total width of the \f$W\f$ boson in GeV
     */
    virtual double GammaW() const;
    

    ////////////////////////////////////////////////////////////////////////    
protected:    
    double a, b, c_u, c_d, c_e, d_3, d_4, LambdaNP_in;
     /**
     * @brief A set method to fix the parameters of the model.
     * @param[in] name a string with the parameter name
     * @param[in] value the value to be asigned to the parameter specified by name
     */
    virtual void setParameter(const std::string name, const double& value);

};

#endif	/* NPHIGGSST_H */

