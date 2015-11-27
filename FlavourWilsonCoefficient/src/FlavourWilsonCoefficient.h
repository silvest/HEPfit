/* 
 * Copyright (C) 2015 HEPfit Collaboration
 * All rights reserved.
 *
 * For the licensing terms see doc/COPYING.
 */

#ifndef FLAVOURWILSONCOEFFICIENT_H
#define	FLAVOURWILSONCOEFFICIENT_H

#include "StandardModel.h"
#include "gslpp.h"
#include "FlavourWilsonCoefficientMatching.h"

/**
 * @class FlavourWilsonCoefficient
 * @brief Model for NP contributions to flavour.
 */
class FlavourWilsonCoefficient: public StandardModel {
public:

    static const int NFlavourWilsonCoefficientvars = 36;

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
    
    virtual FlavourWilsonCoefficientMatching* getMyMatching() const
    {
        return myFlavourWilsonCoefficientMatching;
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
    
protected: 
    
    virtual void setParameter(const std::string, const double&);

private:

    FlavourWilsonCoefficientMatching* myFlavourWilsonCoefficientMatching;
    
    gslpp::complex DC1, DC2, DC3, DC4, DC5, DC6, DC7, DC8, DC9, DC10, DC7g, DC8g;
    gslpp::complex DC7p, DC8p, DC9p, DC10p, DC7gp, DC8gp;
    double absDC1, absDC2, absDC3, absDC4, absDC5, absDC6, absDC7, absDC8, absDC9, absDC10, absDC7g, absDC8g;
    double argDC1, argDC2, argDC3, argDC4, argDC5, argDC6, argDC7, argDC8, argDC9, argDC10, argDC7g, argDC8g;
    double absDC7p, absDC8p, absDC9p, absDC10p, absDC7gp, absDC8gp;
    double argDC7p, argDC8p, argDC9p, argDC10p, argDC7gp, argDC8gp;
    
      
};

#endif	/* FLAVOURWILSONCOEFFICIENT_H */

