/* 
 * Copyright (C) 2013 SusyFit Collaboration
 * All rights reserved.
 *
 * For the licensing terms see doc/COPYING.
 */

#ifndef NPHIGGSST_H
#define	NPHIGGSST_H

#include "NPbase.h"

/**
 * @class NPHiggsST
 * @brief A class for new physics with a non-standard @f$HVV@f$ coupling. 
 * @ingroup NewPhysics
 * @author SusyFit Collaboration
 * @copyright GNU General Public License
 * @details  
 */
class NPHiggsST : public NPbase {
public:
    static const int NNPHIGGSSTvars = 8;
    static const std::string NPHIGGSSTvars[NNPHIGGSSTvars];
    
    NPHiggsST();

    virtual std::string ModelName() const 
    {
        return "NPHiggsST";
    }

    virtual bool InitializeModel();
    virtual void setEWSMflags(EWSM& myEWSM);

    virtual bool Init(const std::map<std::string, double>& DPars);    
    virtual bool Update(const std::map<std::string, double>& DPars);
    virtual bool CheckParameters(const std::map<std::string, double>& DPars);

    virtual bool setFlag(const std::string, const bool&); 
    virtual bool CheckFlags() const;

    
    ////////////////////////////////////////////////////////////////////////     
    
    /**
     * @return the value of the oblique parameter \f$S\f$.
     */
    virtual double obliqueS() const;
        
    /**
     * @return the value of the oblique parameter \f$T\f$
     */
    virtual double obliqueT() const;
    
    /**
     * @return the value of the oblique parameter \f$U\f$ (0 in this class of new
     * physics models)
     */
    virtual double obliqueU() const;


    ////////////////////////////////////////////////////////////////////////

   /**
     *@brief Auxiliary \f$\varepsilon_1\f$ function. Needed for consistency within
     * the structure of the code.
     */        
    virtual double epsilon1() const;

   /**
     *@brief Auxiliary \f$\varepsilon_2\f$ function. Needed for consistency within
     * the structure of the code.
     */    
    virtual double epsilon2() const;

   /**
     *@brief Auxiliary \f$\varepsilon_3\f$ function. Needed for consistency within
     * the structure of the code.
     */     
    virtual double epsilon3() const;

   /**
     *@brief Auxiliary \f$\varepsilon_b\f$ function. Needed for consistency within
     * the structure of the code.
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
     * @return the total width of the \f$W\f$ boson in GeV
     */
    virtual double GammaW() const;
    

    ////////////////////////////////////////////////////////////////////////    
protected:    
    double a, b, c_u, c_d, c_e, d_3, d_4, LambdaNP_in;
    virtual void setParameter(const std::string name, const double& value);

};

#endif	/* NPHIGGSST_H */

