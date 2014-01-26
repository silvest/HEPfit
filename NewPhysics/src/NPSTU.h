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
 * parameters \f$S,~T\f$ and \f$U\f$ 
 * @ingroup NewPhysics
 * @author SusyFit Collaboration
 * @copyright GNU General Public License
 * @details This class contains the necessary functions to compute new physics tree-level corrections to electroweak precision
 * observables, in the form of contributions to the Peskin-Takeuchi oblique parameters \cite . These corrections are
 * parameterized in terms of the \f$S,~T\f$ and \f$U\f$ contributions to \f$M_W\f$, and to \f$Z\f$-pole observables
 * through the corrections to the different neutral-current effective couplings to leptons and quarks. The contributions
 * to the later are implemented in the \b NPbase class.
 */
class NPSTU : public NPbase {
public:
    static const int NSTUvars = 3;
    static const std::string STUvars[NSTUvars];
    
    /**
     * @brief Constructor.
     */
    NPSTU();

    virtual std::string ModelName() const 
    {
        return "NPSTU";
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
     * @return the oblique parameter S
     */
    virtual double obliqueS() const 
    {
        return myObliqueS;
    }

    /**
     * @return the oblique parameter T
     */
    virtual double obliqueT() const 
    {
        return myObliqueT;
    }

    /**
     * @return the oblique parameter U
     */
    virtual double obliqueU() const 
    {
        return myObliqueU;
    }


    ////////////////////////////////////////////////////////////////////////

    /**
     * @return the value of the @f$\epsilon_1@f$ parameter
     */
    double epsilon1() const;

    /**
     * @return the value of the @f$\epsilon_2@f$ parameter
     */
    double epsilon2() const;

    /**
     * @return the value of the @f$\epsilon_3@f$ parameter
     */
    double epsilon3() const;

    /**
     * @return the value of the @f$\epsilon_b@f$ parameter
     */
    double epsilonb() const;

    
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
    double myObliqueS, myObliqueT, myObliqueU;
    virtual void setParameter(const std::string name, const double& value);

};

#endif	/* NPSTU_H */

