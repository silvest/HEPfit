/* 
 * Copyright (C) 2015 HEPfit Collaboration
 * All rights reserved.
 *
 * For the licensing terms see doc/COPYING.
 */

#ifndef FLAVOURWILSONCOEFFICIENT_H
#define	FLAVOURWILSONCOEFFICIENT_H

#include <StandardModel.h>
#include "FlavourWilsonCoefficientMatching.h"

/**
 * @class FlavourWilsonCoefficient
 * @brief Model for NP contributions to flavour.
 */
class FlavourWilsonCoefficient: public StandardModel {
public:

    static const int NFlavourWilsonCoefficientvars = 12;

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

    /**
     *
     * @return \f$\Delta C_1$\f
     */
    double getDC1() const {
        return DC1;
    }
    
    /**
     *
     * @return \f$\Delta C_2$\f
     */
    double getDC2() const {
        return DC2;
    }
    
    /**
     *
     * @return \f$\Delta C_3$\f
     */
    double getDC3() const {
        return DC3;
    }
    
    /**
     *
     * @return \f$\Delta C_4$\f
     */
    double getDC4() const {
        return DC4;
    }
    
    /**
     *
     * @return \f$\Delta C_5$\f
     */
    double getDC5() const {
        return DC5;
    }
    
    /**
     *
     * @return \f$\Delta C_6$\f
     */
    double getDC6() const {
        return DC6;
    }
    
    /**
     *
     * @return \f$\Delta C_7$\f
     */
    double getDC7() const {
        return DC7;
    }
    
    /**
     *
     * @return \f$\Delta C_8$\f
     */
    double getDC8() const {
        return DC8;
    }
    
    /**
     *
     * @return \f$\Delta C_9$\f
     */
    double getDC9() const {
        return DC9;
    }
    
    /**
     *
     * @return \f$\Delta C_10$\f
     */
    double getDC10() const {
        return DC10;
    }
    
    /**
     *
     * @return \f$\Delta C_{7\gamma}$\f
     */
    double getDC7g() const {
        return DC7g;
    }
    
    /**
     *
     * @return \f$\Delta C_{8g}$\f
     */
    double getDC8g() const {
        return DC8g;
    }
    
protected: 
    
    virtual void setParameter(const std::string, const double&);

private:

    FlavourWilsonCoefficientMatching* myFlavourWilsonCoefficientMatching;
    
    double DC1, DC2, DC3, DC4, DC5, DC6, DC7, DC8, DC9, DC10, DC7g, DC8g;
    
      
};

#endif	/* FLAVOURWILSONCOEFFICIENT_H */

