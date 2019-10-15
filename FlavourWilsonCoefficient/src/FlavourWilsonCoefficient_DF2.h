/* 
 * Copyright (C) 2019 HEPfit Collaboration
 *
 *
 * For the licensing terms see doc/COPYING.
 */


#ifndef FLAVOURWILSONCOEFFICIENT_DF2_H
#define FLAVOURWILSONCOEFFICIENT_DF2_H

#include "StandardModel.h"
#include "gslpp.h"
#include "FlavourWilsonCoefficient_DF2Matching.h"


class FlavourWilsonCoefficient_DF2: public StandardModel {
public:
    
    static const int NFlavourWilsonCoefficient_DF2vars = 44;
    
    static const std::string FlavourWilsonCoefficient_DF2vars[NFlavourWilsonCoefficient_DF2vars];
    
    FlavourWilsonCoefficient_DF2();
    
    virtual bool InitializeModel();
        
    virtual bool PostUpdate();
    
    virtual bool CheckParameters(const std::map<std::string, double>& DPars);
    
    virtual FlavourWilsonCoefficient_DF2Matching& getMatching() const
    {
        return FWCM.getObj();
    }
        
    gslpp::vector<gslpp::complex> getC_bd() const
    {
        return C_bd;
    }

    gslpp::vector<gslpp::complex> getC_bs() const
    {
        return C_bs;
    }

    gslpp::vector<gslpp::complex> getC_c() const
    {
        return C_c;
    }

    gslpp::vector<gslpp::complex> getC_s() const
    {
        return C_s;
    }

    /**
     * 
     * @return the scale at which the NP Wilson coefficients in the \f$bd$\f sector are defined
     */
    double getWCscale_bd() const
    {
        return WCscale_bd;
    }

    /**
     * 
     * @return the scale at which the NP Wilson coefficients in the \f$bs$\f sector are defined
     */
    double getWCscale_bs() const
    {
        return WCscale_bs;
    }

    /**
     * 
     * @return the scale at which the NP Wilson coefficients in the \f$cu$\f sector are defined
     */
    double getWCscale_c() const
    {
        return WCscale_c;
    }

    /**
     * 
     * @return the scale at which the NP Wilson coefficients in the \f$sd$\f sector are defined
     */
    double getWCscale_s() const
    {
        return WCscale_s;
    }

protected: 
    
    virtual void setParameter(const std::string, const double&);
    mutable Matching<FlavourWilsonCoefficient_DF2Matching,FlavourWilsonCoefficient_DF2> FWCM;

private:
    double reC1_s, reC2_s, reC3_s, reC4_s, reC5_s;
    double reC1_c, reC2_c, reC3_c, reC4_c, reC5_c;
    double reC1_bd, reC2_bd, reC3_bd, reC4_bd, reC5_bd;
    double reC1_bs, reC2_bs, reC3_bs, reC4_bs, reC5_bs;
    double imC1_s, imC2_s, imC3_s, imC4_s, imC5_s;
    double imC1_c, imC2_c, imC3_c, imC4_c, imC5_c;
    double imC1_bd, imC2_bd, imC3_bd, imC4_bd, imC5_bd;
    double imC1_bs, imC2_bs, imC3_bs, imC4_bs, imC5_bs;
    gslpp::vector<gslpp::complex> C_s, C_c, C_bd, C_bs;
    double WCscale_s, WCscale_c, WCscale_bd, WCscale_bs;

};

#endif /* FLAVOURWILSONCOEFFICIENT_DF2_H */

