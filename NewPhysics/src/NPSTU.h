/* 
 * Copyright (C) 2012-2013 SusyFit Collaboration
 * All rights reserved.
 *
 * For the licensing terms see doc/COPYING.
 */

#ifndef NPSTU_H
#define	NPSTU_H

#include "NPZbbbar.h"

/**
 * @class NPSTU
 * @brief A class for new physics with the oblique parameters. 
 * @ingroup NewPhysics
 * @author SusyFit Collaboration
 * @copyright GNU General Public License
 * @details  
 */
class NPSTU : public NPZbbbar {
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

    virtual bool Update(const std::map<std::string, double>& DPars);
    virtual bool Init(const std::map<std::string, double>& DPars);    
    virtual bool CheckParameters(const std::map<std::string, double>& DPars);

    virtual bool InitializeModel();  
    virtual void SetEWSMflags(EWSM& myEWSM);    

    virtual bool SetFlag(const std::string, const bool&); 

    
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
    
protected:    
    double myObliqueS, myObliqueT, myObliqueU;
    virtual void SetParameter(const std::string name, const double& value);

    
    ////////////////////////////////////////////////////////////////////////     
    
private:

};

#endif	/* NPSTU_H */

