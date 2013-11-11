/* 
 * Copyright (C) 2012-2013 SusyFit Collaboration
 * All rights reserved.
 *
 * For the licensing terms see doc/COPYING.
 */

#ifndef NPSTUVWXY_H
#define	NPSTUVWXY_H

#include <cmath>
#include "NPbase.h"

/**
 * @class NPSTUVWXY
 * @brief A class for new physics with the extended oblique parameters. 
 * @ingroup NewPhysics
 * @author SusyFit Collaboration
 * @copyright GNU General Public License
 * @details  
 */
class NPSTUVWXY : public NPbase {
public:
    static const int NSTUVWXYvars = 7;
    static const std::string STUVWXYvars[NSTUVWXYvars];
    
    /**
     * @brief NPSTUVWXY constructor
     */
    NPSTUVWXY();

    virtual std::string ModelName() const 
    {
        return "NPSTUVWXY";
    }

    virtual bool Update(const std::map<std::string, double>& DPars);
    virtual bool Init(const std::map<std::string, double>& DPars);    
    virtual bool CheckParameters(const std::map<std::string, double>& DPars);

    virtual bool InitializeModel();  
    virtual void setEWSMflags(EWSM& myEWSM);    

    virtual bool setFlag(const std::string, const bool&); 
    virtual bool CheckFlags() const;
    
    
    ////////////////////////////////////////////////////////////////////////     
    
    /**
     * @return Oblique parameter \hat{S}
     */
    virtual double obliqueShat() const
    {
        return myObliqueShat;
    }

    /**
     * @return Oblique parameter \hat{T}
     */
    virtual double obliqueThat() const 
    {
        return myObliqueThat;
    }

    /**
     * @return Oblique parameter \hat{U}
     */
    virtual double obliqueUhat() const 
    {
        return myObliqueUhat;
    }

    /**
     * @return Oblique parameter V
     */
    virtual double obliqueV() const 
    {
        return myObliqueV;
    }

    /**
     * @return Oblique parameter W
     */
    virtual double obliqueW() const 
    {
        return myObliqueW;
    }

    /**
     * @return Oblique parameter X
     */
    virtual double obliqueX() const 
    {
        return myObliqueX;
    }

    /**
     * @return Oblique parameter Y
     */
    virtual double obliqueY() const 
    {
        return myObliqueY;
    }


    ////////////////////////////////////////////////////////////////////////

    virtual double epsilon1() const;

    virtual double epsilon2() const;

    virtual double epsilon3() const;

    virtual double epsilonb() const;

    
    ////////////////////////////////////////////////////////////////////////     
    // Combinations of the extended oblique parameters
    
    /**
     * @return Oblique parameter S'
     */
    virtual double obliqueS() const;

    /**
     * @return Oblique parameter T'
     */
    virtual double obliqueT() const;

    /**
     * @return Oblique parameter U'
     */
    virtual double obliqueU() const;


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
    virtual void setParameter(const std::string name, const double& value);
    double myObliqueShat, myObliqueThat, myObliqueUhat;
    double myObliqueV, myObliqueW, myObliqueX, myObliqueY;

    
    ////////////////////////////////////////////////////////////////////////     
    
private:

};

#endif	/* NPSTUVWXY_H */

