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
    double a, b, c_u, c_d, c_e, d_3, d_4, LambdaNP_in;
    virtual void setParameter(const std::string name, const double& value);

};

#endif	/* NPHIGGSST_H */

