/* 
 * Copyright (C) 2012-2013 SusyFit Collaboration
 * All rights reserved.
 *
 * For the licensing terms see doc/COPYING.
 */

#ifndef NPEPSILONS_H
#define	NPEPSILONS_H

#include <EWSM.h>
#include "NPbase.h"

/**
 * @class NPEpsilons
 * @brief A class for new physics in the form of contributions to the \f$\varepsilon_{1,2,3,b}\f$ parameters.
 * Both SM and new physics contributions to \f$\varepsilon_i\f$ are parameterized. 
 * @ingroup NewPhysics
 * @author SusyFit Collaboration
 * @copyright GNU General Public License
 * @details  
 *
 * Flags:
 * \li \b FlagEpsilon1SM:&nbsp; Return the SM (True) or SM plus new physics (False)
 * contribution to \f$\varepsilon_1\f$.
 * \li \b FlagEpsilon2SM:&nbsp; Return the SM (True) or SM plus new physics (False)
 * contribution to \f$\varepsilon_2\f$.
 * \li \b FlagEpsilon3SM:&nbsp; Return the SM (True) or SM plus new physics (False)
 * contribution to \f$\varepsilon_3\f$.
 * \li \b FlagEpsilonbSM:&nbsp; Return the SM (True) or SM plus new physics (False)
 * contribution to \f$\varepsilon_b\f$.
 */
class NPEpsilons : public NPbase  {
public:
    static const int NEPSILONvars = 4;
    static const std::string EPSILONvars[NEPSILONvars];
    static const int NEPSILONflags = 4;
    static const std::string EPSILONflags[NEPSILONflags];
    
    /**
     * @brief NPEpsilons constructor
     */
    NPEpsilons();

    virtual std::string ModelName() const 
    {
        return "NPEpsilons";
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
     * @return the SM value (FlagEpsilon1SM=True)
     * or the SM plus new physics value (FlagEpsilon1SM=False) of \f$\varepsilon_1\f$
     */
    virtual double epsilon1() const 
    {
        if (FlagEpsilon1SM) 
            return myEWSM->epsilon1_SM();
        else
            return myEpsilon_1;
    }

    /**
     * @return the SM value (FlagEpsilon2SM=True)
     * or the SM plus new physics value (FlagEpsilon2SM=False) of \f$\varepsilon_2\f$
     */
    virtual double epsilon2() const 
    {
        if (FlagEpsilon2SM) 
            return myEWSM->epsilon2_SM();
        else
            return myEpsilon_2;
    }

    /**
     * @return the SM value (FlagEpsilon3SM=True)
     * or the SM plus new physics value (FlagEpsilon3SM=False) of \f$\varepsilon_3\f$
     */
    virtual double epsilon3() const 
    {
        if (FlagEpsilon3SM) 
            return myEWSM->epsilon3_SM();
        else
            return myEpsilon_3;
    }
 
    /**
     * @return the SM value (FlagEpsilonbSM=True)
     * or the SM plus new physics value (FlagEpsilonbSM=False) of \f$\varepsilon_b\f$
     */
    virtual double epsilonb() const 
    {
        if (FlagEpsilonbSM) 
            return myEWSM->epsilonb_SM();
        else
            return myEpsilon_b;
    }    

    
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
protected:    
    double myEpsilon_1, myEpsilon_2, myEpsilon_3, myEpsilon_b;
    virtual void setParameter(const std::string name, const double& value);
    
    
    ////////////////////////////////////////////////////////////////////////         
private:
    bool FlagEpsilon1SM, FlagEpsilon2SM, FlagEpsilon3SM, FlagEpsilonbSM;
    
};

#endif	/* NPEPSILONS_H */

