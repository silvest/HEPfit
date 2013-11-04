/* 
 * Copyright (C) 2012-2013 SusyFit Collaboration
 * All rights reserved.
 *
 * For the licensing terms see doc/COPYING.
 */

#ifndef NPSTU_H
#define	NPSTU_H

#include <EW_BURGESS.h>
#include "StandardModel.h"

/**
 * @class NPSTU
 * @brief A class for new physics with the oblique parameters. 
 * @ingroup NewPhysics
 * @author SusyFit Collaboration
 * @copyright GNU General Public License
 * @details  
 */
class NPSTU : public StandardModel {
public:
    static const int NSTUvars = 3;
    static const std::string STUvars[NSTUvars];
    static const int NSTUflags = 2;
    static const std::string STUflags[NSTUflags];
    
    /**
     * @brief NPSTU constructor.
     */
    NPSTU();

    virtual std::string ModelName() const 
    {
        return "NPSTU";
    }

    virtual bool Update(const std::map<std::string, double>& DPars);
    virtual bool Init(const std::map<std::string, double>& DPars);    
    virtual bool CheckParameters(const std::map<std::string, double>& DPars);

    virtual bool InitializeModel();  
    virtual void SetEWSMflags(EWSM& myEWSM);    

    virtual bool SetFlag(const std::string, const bool&); 

    ////////////////////////////////////////////////////////////////////////

    bool IsFlagEWCHMN() const
    {
        return FlagEWCHMN;
    }
    
    bool IsFlagEWBURGESS() const
    {
        return FlagEWBURGESS;
    }

    ////////////////////////////////////////////////////////////////////////

    /**
     * @return Oblique parameter S.
     */
    virtual double obliqueS() const 
    {
        return myObliqueS;
    }

    /**
     * @return Oblique parameter T.
     */
    virtual double obliqueT() const 
    {
        return myObliqueT;
    }

    /**
     * @return Oblique parameter U.
     */
    virtual double obliqueU() const 
    {
        return myObliqueU;
    }


    ////////////////////////////////////////////////////////////////////////

    virtual double epsilon1() const;

    virtual double epsilon2() const;

    virtual double epsilon3() const;

    virtual double epsilonb() const;


    ////////////////////////////////////////////////////////////////////////     

    /**
     * @return The W boson mass.
     */
    virtual double Mw() const;

    /**
     * @return @f$M_W^2/M_Z^2@f$.
     */
    virtual double cW2() const;

    /**
     * @return @f$1-M_W^2/M_Z^2@f$.
     */
    virtual double sW2() const;

    /**
     * @return The total width of the W boson.
     */
    virtual double GammaW() const;

    
    ////////////////////////////////////////////////////////////////////////
    
protected:    
    double myObliqueS, myObliqueT, myObliqueU;
    virtual void setParameters(const std::string name, const double& value);

    
    ////////////////////////////////////////////////////////////////////////     
    
private:
    bool FlagEWCHMN;
    bool FlagEWBURGESS;
    const EW_BURGESS myEW_BURGESS;
    
};

#endif	/* NPSTU_H */

