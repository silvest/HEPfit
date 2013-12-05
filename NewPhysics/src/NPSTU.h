/* 
 * Copyright (C) 2012-2013 SusyFit Collaboration
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
 * @brief A class for new physics with the oblique parameters. 
 * @ingroup NewPhysics
 * @author SusyFit Collaboration
 * @copyright GNU General Public License
 * @details  
 */
class NPSTU : public NPbase {
public:
    static const int NSTUvars = 3;
    static const std::string STUvars[NSTUvars];
    
    /**
     * @brief NPSTU constructor.
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
    virtual void setParameter(const std::string name, const double& value);

};

#endif	/* NPSTU_H */

