/* 
 * Copyright (C) 2012-2013 SusyFit Collaboration
 * All rights reserved.
 *
 * For the licensing terms see doc/COPYING.
 */

#ifndef NPSTUVWXY_H
#define	NPSTUVWXY_H

#include <cmath>
#include "NPZbbbar.h"

/**
 * @class NPSTUVWXY
 * @brief A class for new physics with the extended oblique parameters. 
 * @ingroup NewPhysics
 * @author SusyFit Collaboration
 * @copyright GNU General Public License
 * @details  
 */
class NPSTUVWXY : public NPZbbbar {
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
    
    
    ////////////////////////////////////////////////////////////////////////     

    bool SetFlag(const std::string, const bool&); 
    
    
    ////////////////////////////////////////////////////////////////////////     
    
    /**
     * @return Oblique parameter \hat{S}
     */
    double obliqueShat() const
    {
        return myObliqueShat;
    }

    /**
     * @return Oblique parameter \hat{T}
     */
    double obliqueThat() const 
    {
        return myObliqueThat;
    }

    /**
     * @return Oblique parameter \hat{U}
     */
    double obliqueUhat() const 
    {
        return myObliqueUhat;
    }

    /**
     * @return Oblique parameter V
     */
    double obliqueV() const 
    {
        return myObliqueV;
    }

    /**
     * @return Oblique parameter W
     */
    double obliqueW() const 
    {
        return myObliqueW;
    }

    /**
     * @return Oblique parameter X
     */
    double obliqueX() const 
    {
        return myObliqueX;
    }

    /**
     * @return Oblique parameter Y
     */
    double obliqueY() const 
    {
        return myObliqueY;
    }
    
    
    ////////////////////////////////////////////////////////////////////////     
    // Combinations of the extended oblique parameters
    
    /**
     * @return Oblique parameter S'
     */
    double obliqueS() const 
    {
        double s0 = sqrt(s02());
        double c0 = sqrt(c02());
        return ( ( myObliqueShat - myObliqueW + myObliqueX/(s0*c0) - myObliqueY )
                 * 4.0*s02()/alphaMz() );
    }

    /**
     * @return Oblique parameter T'
     */
    double obliqueT() const 
    {
        double s0 = sqrt(s02());
        double c0 = sqrt(c02());
        return ( ( myObliqueThat - myObliqueW + 2.0*s0/c0*myObliqueX 
                   - s02()/c02()*myObliqueY )/alphaMz() );
    }

    /**
     * @return Oblique parameter U'
     */
    double obliqueU() const 
    {
        double s0 = sqrt(s02());
        double c0 = sqrt(c02());
        return ( ( - myObliqueUhat + myObliqueV + myObliqueW 
                   - 2.0*s0/c0*myObliqueX )*4.0*s02()/alphaMz() );
    }
    
    
    ////////////////////////////////////////////////////////////////////////     
    
protected:    
    virtual void SetParameter(const std::string name, const double& value);
    double myObliqueShat, myObliqueThat, myObliqueUhat;
    double myObliqueV, myObliqueW, myObliqueX, myObliqueY;

    ////////////////////////////////////////////////////////////////////////     
    
private:

};

#endif	/* NPSTUVWXY_H */

