/* 
 * Copyright (C) 2018 HEPfit Collaboration
 *
 *
 * For the licensing terms see doc/COPYING.
 */

#ifndef REALWEAKEFTLFV_H
#define REALWEAKEFTLFV_H

#include "StandardModel.h"
#include "gslpp.h"
#include "RealWeakEFTLFVMatching.h"

/**
 * @class RealWeakEFTLFV
 * @brief Model for NP contributions to flavour.
 */
class RealWeakEFTLFV: public StandardModel {
public:

    static const int NRealWeakEFTLFVvars = 21;

    static const std::string RealWeakEFTLFVvars[NRealWeakEFTLFVvars];
    
    /**
     * @brief FlavourWilsonCoefficient constructor
     */
    RealWeakEFTLFV();
    
    /**
     * @brief FlavourWilsonCoefficient destructor
     */
    ~RealWeakEFTLFV();
    
    virtual bool InitializeModel();
    
    virtual bool Init(const std::map<std::string, double>& DPars);
    
    virtual bool PreUpdate();
    
    virtual bool Update(const std::map<std::string, double>& DPars);
    
    virtual bool PostUpdate();
    
    virtual bool CheckParameters(const std::map<std::string, double>& DPars);
    
    virtual bool setFlag(const std::string name, const bool value);
    
    virtual RealWeakEFTLFVMatching& getMatching() const
    {
        return ReWEFTM.getObj();
    }
    
    /**
     *
     * @return \f$ C_7 $\f
     */
    double getC7() const {
        return C7;
    }
    
    /**
     *
     * @return \f$ C_7' $\f
     */
    double getC7p() const {
        return C7p;
    }
    
    /**
     *
     * @return \f$ C_8 $\f
     */
    double getC8() const {
        return C8;
    }
    
        /**
     *
     * @return \f$ C_8'$\f
     */
    double getC8p() const {
        return C8p;
    }
    
    /**
     *
     * @return \f$ C_9^{11} $\f
     */
    double getC9_11() const {
        return C9_11;
    }
    
        /**
     *
     * @return \f$ C_9'^{11} $\f
     */
    double getC9p_11() const {
        return C9p_11;
    }
    
    /**
     *
     * @return \f$ C_{10}^{11} $\f
     */
    double getC10_11() const {
        return C10_11;
    }
    
        /**
     *
     * @return \f$ C_{10}'^{11} $\f
     */
    double getC10p_11() const {
        return C10p_11;
    }
    
    /**
     *
     * @return \f$ C_S^{11} $\f
     */
    double getCS_11() const {
        return CS_11;
    }
    
        /**
     *
     * @return \f$ C_S'^{11} $\f
     */
    double getCSp_11() const {
        return CSp_11;
    }
    
    /**
     *
     * @return \f$ C_P^{11} $\f
     */
    double getCP_11() const {
        return CP_11;
    }
    
        /**
     *
     * @return \f$ C_P'^{11} $\f
     */
    double getCPp_11() const {
        return CPp_11;
    }
    
    /**
     *
     * @return \f$ C_9^{22} $\f
     */
    double getC9_22() const {
        return C9_22;
    }
    
        /**
     *
     * @return \f$ C_9'^{22} $\f
     */
    double getC9p_22() const {
        return C9p_22;
    }
    
    /**
     *
     * @return \f$ C_{10}^{22} $\f
     */
    double getC10_22() const {
        return C10_22;
    }
    
        /**
     *
     * @return \f$ C_{10}'^{22} $\f
     */
    double getC10p_22() const {
        return C10p_22;
    }
    
    /**
     *
     * @return \f$ C_S^{22} $\f
     */
    double getCS_22() const {
        return CS_22;
    }
    
        /**
     *
     * @return \f$ C_S'^{22} $\f
     */
    double getCSp_22() const {
        return CSp_22;
    }
    
    /**
     *
     * @return \f$ C_P^{22} $\f
     */
    double getCP_22() const {
        return CP_22;
    }
    
        /**
     *
     * @return \f$ C_P'^{22} $\f
     */
    double getCPp_22() const {
        return CPp_22;
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
    mutable Matching<RealWeakEFTLFVMatching,RealWeakEFTLFV> ReWEFTM;

private:

    double C7, C7p, C8, C8p;
    double C9_11, C9p_11, C10_11, C10p_11, CS_11, CSp_11, CP_11, CPp_11;
    double C9_22, C9p_22, C10_22, C10p_22, CS_22, CSp_22, CP_22, CPp_22;
    
    double WCscale;
    
      
};

#endif /* REALWEAKEFTLFV_H */

