/* 
 * Copyright (C) 2013 SusyFit Collaboration
 * All rights reserved.
 *
 * For the licensing terms see doc/COPYING.
 */

#ifndef NPHIGGSST_H
#define	NPHIGGSST_H

#include "NPZbbbar.h"

/**
 * @class NPHiggsST
 * @brief A class for new physics with a non-standard @f$HVV@f$ coupling. 
 * @ingroup NewPhysics
 * @author SusyFit Collaboration
 * @copyright GNU General Public License
 * @details  
 */
class NPHiggsST : public NPZbbbar {
public:
    static const int NNPHIGGSSTvars = 7;
    static const std::string NPHIGGSSTvars[NNPHIGGSSTvars];
    
    NPHiggsST();

    virtual std::string ModelName() const 
    {
        return "NPHiggsST";
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
    virtual double obliqueS() const;
        
    /**
     * @return Oblique parameter T.
     */
    virtual double obliqueT() const;
    
    /**
     * @return Oblique parameter U.
     */
    virtual double obliqueU() const;

    
    ////////////////////////////////////////////////////////////////////////     

    virtual double epsilon1() const;

    virtual double epsilon2() const;

    virtual double epsilon3() const;

    virtual double epsilonb() const;


    ////////////////////////////////////////////////////////////////////////
    
protected:    
    double a, b, c_u, c_d, c_e, d_3, d_4;
    virtual void SetParameter(const std::string name, const double& value);
    
    
    ////////////////////////////////////////////////////////////////////////   

private:

};

#endif	/* NPHIGGSST_H */

