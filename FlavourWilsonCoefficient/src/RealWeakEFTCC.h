/* 
 * Copyright (C) 2018 HEPfit Collaboration
 *
 *
 * For the licensing terms see doc/COPYING.
 */

#ifndef REALWEAKEFTCC_H
#define REALWEAKEFTCC_H

#include "StandardModel.h"
#include "gslpp.h"

/**
 * @class RealWeakEFTCC
 * @brief Model for NP contributions to flavour.
 */
class RealWeakEFTCC: public StandardModel {
public:

    static const int NRealWeakEFTCCvars = 5;

    static const std::string RealWeakEFTCCvars[NRealWeakEFTCCvars];
    static const std::string RealWeakEFTCCPMvars[NRealWeakEFTCCvars];
    
    /**
     * @brief FlavourWilsonCoefficient constructor
     */
    RealWeakEFTCC(const unsigned int basis);
    
    /**
     * @brief FlavourWilsonCoefficient destructor
     */
    ~RealWeakEFTCC();
    
    virtual bool InitializeModel();
    
    virtual bool Init(const std::map<std::string, double>& DPars);
    
    virtual bool PreUpdate();
    
    virtual bool Update(const std::map<std::string, double>& DPars);
    
    virtual bool PostUpdate();
    
    virtual bool CheckParameters(const std::map<std::string, double>& DPars);
    
    virtual bool setFlag(const std::string name, const bool value);
    
    /**
     *
     * @return \f$ C_S or C_{S}^{L} $\f
     */
    double getCCC1() const {
        return C1;
    }
    
    /**
     *
     * @return \f$ C_P or C_{S}^{R} $\f
     */
    double getCCC2() const {
        return C2;
    }
    
    /**
     *
     * @return \f$ C_V or C_{V}^{L} $\f
     */
    double getCCC3() const {
        return C3;
    }
    
        /**
     *
     * @return \f$ C_A or or C_{V}^{R}$\f
     */
    double getCCC4() const {
        return C4;
    }
    
    /**
     *
     * @return \f$ C_T $\f
     */
    double getC5() const {
        return C5;
    }
    
protected: 
    
    virtual void setParameter(const std::string, const double&);

private:

    double C1, C2, C3, C4, C5;   
    const unsigned int basis;
};

#endif /* REALWEAKEFTCC_H */

