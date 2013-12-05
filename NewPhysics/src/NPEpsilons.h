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
 * @brief A class for new physics with the @f$\varepsilon@f$ parameters. 
 * @ingroup NewPhysics
 * @author SusyFit Collaboration
 * @copyright GNU General Public License
 * @details  
 *
 * Flags:
 * \li \b FlagEpsilon1SM:&nbsp;
 * \li \b FlagEpsilon2SM:&nbsp;
 * \li \b FlagEpsilon3SM:&nbsp;
 * \li \b FlagEpsilonbSM:&nbsp;
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

    virtual double epsilon1() const 
    {
        if (FlagEpsilon1SM) 
            return myEWSM->epsilon1_SM();
        else
            return myEpsilon_1;
    }

    virtual double epsilon2() const 
    {
        if (FlagEpsilon2SM) 
            return myEWSM->epsilon2_SM();
        else
            return myEpsilon_2;
    }

    virtual double epsilon3() const 
    {
        if (FlagEpsilon3SM) 
            return myEWSM->epsilon3_SM();
        else
            return myEpsilon_3;
    }

    virtual double epsilonb() const 
    {
        if (FlagEpsilonbSM) 
            return myEWSM->epsilonb_SM();
        else
            return myEpsilon_b;
    }    

    
    ////////////////////////////////////////////////////////////////////////     
    
    /**
     * @return the W boson mass
     */
    virtual double Mw() const;

    /**
     * @return Mw^2/Mz^2
     */
    virtual double cW2() const;
    
    /**
     * @return 1-Mw^2/Mz^2
     */
    virtual double sW2() const;

    /**
     * @return the total width of the W boson
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

