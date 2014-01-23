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

    virtual bool InitializeModel();
    virtual void setEWSMflags(EWSM& myEWSM);

    virtual bool Init(const std::map<std::string, double>& DPars);    
    virtual bool Update(const std::map<std::string, double>& DPars);
    virtual bool CheckParameters(const std::map<std::string, double>& DPars);

    virtual bool setFlag(const std::string, const bool&); 
    virtual bool CheckFlags() const;
    
    
    ////////////////////////////////////////////////////////////////////////     
    
    /**
     * @return the value of the oblique parameter \f$\hat{S}\f$
     */
    virtual double obliqueShat() const
    {
        return myObliqueShat;
    }

    /**
     * @return the value of the oblique parameter \f$\hat{T}\f$
     */
    virtual double obliqueThat() const 
    {
        return myObliqueThat;
    }

    /**
     * @return the value of the oblique parameter \f$\hat{U}\f$
     */
    virtual double obliqueUhat() const 
    {
        return myObliqueUhat;
    }

    /**
     * @return the value of the oblique parameter \f$V\f$
     */
    virtual double obliqueV() const 
    {
        return myObliqueV;
    }

    /**
     * @return the value of the oblique parameter \f$W\f$
     */
    virtual double obliqueW() const 
    {
        return myObliqueW;
    }

    /**
     * @return the value of the oblique parameter \f$X\f$
     */
    virtual double obliqueX() const 
    {
        return myObliqueX;
    }

    /**
     * @return the value of the oblique parameter \f$Y\f$
     */
    virtual double obliqueY() const 
    {
        return myObliqueY;
    }


    ////////////////////////////////////////////////////////////////////////

    /**
     * @return the value of the \f$\varepsilon_1\f$ parameter including new physics
     * corrections
     */
    virtual double epsilon1() const;

    /**
     * @return the value of the \f$\varepsilon_2\f$ parameter including new physics
     * corrections
     */
    virtual double epsilon2() const;
    
    /**
      * @return the value of the \f$\varepsilon_3\f$ parameter including new physics
     * corrections
     */
    virtual double epsilon3() const;
    
    /**
     * @return the SM value of the \f$\varepsilon_b\f$ parameter
     */
    virtual double epsilonb() const;

    
    ////////////////////////////////////////////////////////////////////////     
    // Combinations of the extended oblique parameters
    
    /**
     * @return the value of the oblique parameter \f$S'\f$
     */
    virtual double obliqueS() const;

    /**
     * @return the value of the oblique parameter \f$T'\f$
     */
    virtual double obliqueT() const;

    /**
     * @return the value of the oblique parameter \f$U'\f$
     */
    virtual double obliqueU() const;


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
    virtual void setParameter(const std::string name, const double& value);
    double myObliqueShat, myObliqueThat, myObliqueUhat;
    double myObliqueV, myObliqueW, myObliqueX, myObliqueY;

};

#endif	/* NPSTUVWXY_H */

