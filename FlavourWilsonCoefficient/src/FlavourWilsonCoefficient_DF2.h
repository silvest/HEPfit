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
    
    /**
     * 
     * @return \f $C_{1}^{bd}$ \f
     */
    gslpp::complex GetC1_bd() const
    {
        return C1_bd;
    }

    /**
     * 
     * @return  \f $C_{1}^{bs}$ \f
     */
    gslpp::complex GetC1_bs() const
    {
        return C1_bs;
    }

    /**
     * 
     * @return  \f $C_{1}^{cu}$ \f
     */
    gslpp::complex GetC1_c() const
    {
        return C1_c;
    }

    /**
     * 
     * @return  \f $C_{1}^{sd}$ \f
     */
    gslpp::complex GetC1_s() const
    {
        return C1_s;
    }

    /**
     * 
     * @return  \f $C_{2}^{bd}$ \f
     */
    gslpp::complex GetC2_bd() const
    {
        return C2_bd;
    }

    /**
     * 
     * @return  \f $C_{2}^{bs}$ \f
     */
    gslpp::complex GetC2_bs() const
    {
        return C2_bs;
    }

    /**
     * 
     * @return  \f $C_{2}^{cu}$ \f
     */
    gslpp::complex GetC2_c() const
    {
        return C2_c;
    }

    /**
     * 
     * @return  \f $C_{2}^{sd}$ \f
     */
    gslpp::complex GetC2_s() const
    {
        return C2_s;
    }

    /**
     * 
     * @return  \f $C_{3}^{bd}$ \f
     */
    gslpp::complex GetC3_bd() const
    {
        return C3_bd;
    }

    /**
     * 
     * @return  \f $C_{3}^{bs}$ \f
     */
    gslpp::complex GetC3_bs() const
    {
        return C3_bs;
    }

    /**
     * 
     * @return  \f $C_{3}^{cu}$ \f
     */
    gslpp::complex GetC3_c() const
    {
        return C3_c;
    }

    /**
     * 
     * @return  \f $C_{3}^{sd}$ \f
     */
    gslpp::complex GetC3_s() const
    {
        return C3_s;
    }

    /**
     * 
     * @return  \f $C_{4}^{bd}$ \f
     */
    gslpp::complex GetC4_bd() const
    {
        return C4_bd;
    }

    /**
     * 
     * @return  \f $C_{4}^{bs}$ \f
     */
    gslpp::complex GetC4_bs() const
    {
        return C4_bs;
    }

    /**
     * 
     * @return  \f $C_{4}^{cu}$ \f
     */
    gslpp::complex GetC4_c() const
    {
        return C4_c;
    }

    /**
     * 
     * @return  \f $C_{4}^{sd}$ \f
     */
    gslpp::complex GetC4_s() const
    {
        return C4_s;
    }

    /**
     * 
     * @return  \f $C_{5}^{bd}$ \f
     */
    gslpp::complex GetC5_bd() const
    {
        return C5_bd;
    }

    /**
     * 
     * @return  \f $C_{5}^{bs}$ \f
     */
    gslpp::complex GetC5_bs() const
    {
        return C5_bs;
    }

    /**
     * 
     * @return  \f $C_{5}^{cu}$ \f
     */
    gslpp::complex GetC5_c() const
    {
        return C5_c;
    }

    /**
     * 
     * @return  \f $C_{5}^{sd}$ \f
     */
    gslpp::complex GetC5_s() const
    {
        return C5_s;
    }

    /**
     * 
     * @return the scale at which the NP Wilson coefficients in the \f$bd$\f sector are defined
     */
    double GetWCscale_bd() const
    {
        return WCscale_bd;
    }

    /**
     * 
     * @return the scale at which the NP Wilson coefficients in the \f$bs$\f sector are defined
     */
    double GetWCscale_bs() const
    {
        return WCscale_bs;
    }

    /**
     * 
     * @return the scale at which the NP Wilson coefficients in the \f$cu$\f sector are defined
     */
    double GetWCscale_c() const
    {
        return WCscale_c;
    }

    /**
     * 
     * @return the scale at which the NP Wilson coefficients in the \f$sd$\f sector are defined
     */
    double GetWCscale_s() const
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
    gslpp::complex C1_s, C2_s, C3_s, C4_s, C5_s;
    gslpp::complex C1_c, C2_c, C3_c, C4_c, C5_c;
    gslpp::complex C1_bd, C2_bd, C3_bd, C4_bd, C5_bd;
    gslpp::complex C1_bs, C2_bs, C3_bs, C4_bs, C5_bs;
    double WCscale_s, WCscale_c, WCscale_bd, WCscale_bs;

};

#endif /* FLAVOURWILSONCOEFFICIENT_DF2_H */

