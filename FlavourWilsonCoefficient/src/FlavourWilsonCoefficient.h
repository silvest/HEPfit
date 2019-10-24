/* 
 * Copyright (C) 2015 HEPfit Collaboration
 *
 *
 * For the licensing terms see doc/COPYING.
 */

#ifndef FLAVOURWILSONCOEFFICIENT_H
#define	FLAVOURWILSONCOEFFICIENT_H

#include "StandardModel.h"
#include "gslpp.h"
#include "FlavourWilsonCoefficientMatching.h"

/**
 * @addtogroup FlavourWilsonCoefficent
 * @brief A module for EFT models for flavour observaqbles.
 * @details This module has implementations of several EFTs for flavour
 * analyses like the Weak Effective Theory at low energy and generalized
 * EFT models for @f$\Delta F = 1@f$ and @f$\Delta F = 2@f$.
 * @{
 */

/**
 * @class FlavourWilsonCoefficient
 * @brief Model for NP contributions to flavour.
 */
class FlavourWilsonCoefficient: public StandardModel {
public:

    static const int NFlavourWilsonCoefficientvars = 37;

    static const std::string FlavourWilsonCoefficientvars[NFlavourWilsonCoefficientvars];
    
    /**
     * @brief FlavourWilsonCoefficient constructor
     */
    FlavourWilsonCoefficient();
    
    /**
     * @brief FlavourWilsonCoefficient destructor
     */
    ~FlavourWilsonCoefficient();
    
    virtual bool InitializeModel();
    
    virtual bool Init(const std::map<std::string, double>& DPars);
    
    virtual bool PreUpdate();
    
    virtual bool Update(const std::map<std::string, double>& DPars);
    
    virtual bool PostUpdate();
    
    virtual bool CheckParameters(const std::map<std::string, double>& DPars);
    
    virtual bool setFlag(const std::string name, const bool value);
    
    virtual FlavourWilsonCoefficientMatching& getMatching() const
    {
        return FWCM.getObj();
    }
    
    /**
     *
     * @return \f$\Delta C_1$\f
     */
    gslpp::complex getDC1() const {
        return DC1;
    }
    
    /**
     *
     * @return \f$\Delta C_2$\f
     */
    gslpp::complex getDC2() const {
        return DC2;
    }
    
    /**
     *
     * @return \f$\Delta C_3$\f
     */
    gslpp::complex getDC3() const {
        return DC3;
    }
    
    /**
     *
     * @return \f$\Delta C_4$\f
     */
    gslpp::complex getDC4() const {
        return DC4;
    }
    
    /**
     *
     * @return \f$\Delta C_5$\f
     */
    gslpp::complex getDC5() const {
        return DC5;
    }
    
    /**
     *
     * @return \f$\Delta C_6$\f
     */
    gslpp::complex getDC6() const {
        return DC6;
    }
    
    /**
     *
     * @return \f$\Delta C_7$\f
     */
    gslpp::complex getDC7() const {
        return DC7;
    }
    
    /**
     *
     * @return \f$\Delta C_8$\f
     */
    gslpp::complex getDC8() const {
        return DC8;
    }
    
    /**
     *
     * @return \f$\Delta C_9$\f
     */
    gslpp::complex getDC9() const {
        return DC9;
    }
    
    /**
     *
     * @return \f$\Delta C_10$\f
     */
    gslpp::complex getDC10() const {
        return DC10;
    }
    
    /**
     *
     * @return \f$\Delta C_{7\gamma}$\f
     */
    gslpp::complex getDC7g() const {
        return DC7g;
    }
    
    /**
     *
     * @return \f$\Delta C_{8g}$\f
     */
    gslpp::complex getDC8g() const {
        return DC8g;
    }
    
    /**
     *
     * @return \f$\Delta C'_7$\f
     */
    gslpp::complex getDC7p() const {
        return DC7p;
    }
    
    /**
     *
     * @return \f$\Delta C'_8$\f
     */
    gslpp::complex getDC8p() const {
        return DC8p;
    }
    
    /**
     *
     * @return \f$\Delta C'_9$\f
     */
    gslpp::complex getDC9p() const {
        return DC9;
    }
    
    /**
     *
     * @return \f$\Delta C'_10$\f
     */
    gslpp::complex getDC10p() const {
        return DC10p;
    }
    
    /**
     *
     * @return \f$\Delta C'_{7\gamma}$\f
     */
    gslpp::complex getDC7gp() const {
        return DC7gp;
    }
    
    /**
     *
     * @return \f$\Delta C'_{8g}$\f
     */
    gslpp::complex getDC8gp() const {
        return DC8gp;
    }
    
    /**
     *
     * @return the matching scale of the Wilson coefficients
     */
    double getWCscale() const {
        return WCscale;
    }
    
protected: 
    
    virtual void setParameter(const std::string, const double&);
    mutable Matching<FlavourWilsonCoefficientMatching,FlavourWilsonCoefficient> FWCM;

private:

    
    gslpp::complex DC1, DC2, DC3, DC4, DC5, DC6, DC7, DC8, DC9, DC10, DC7g, DC8g;
    gslpp::complex DC7p, DC8p, DC9p, DC10p, DC7gp, DC8gp;
    double reDC1, reDC2, reDC3, reDC4, reDC5, reDC6, reDC7, reDC8, reDC9, reDC10, reDC7g, reDC8g;
    double imDC1, imDC2, imDC3, imDC4, imDC5, imDC6, imDC7, imDC8, imDC9, imDC10, imDC7g, imDC8g;
    double reDC7p, reDC8p, reDC9p, reDC10p, reDC7gp, reDC8gp;
    double imDC7p, imDC8p, imDC9p, imDC10p, imDC7gp, imDC8gp;
    double WCscale;
    
      
};

/**
 * @}
 */

#endif	/* FLAVOURWILSONCOEFFICIENT_H */

