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
    
    /**
     * @brief FlavourWilsonCoefficient constructor
     */
    FlavourWilsonCoefficient_DF2();
    
    /**
     * @brief A method to initialize the model.
     * @details This method, called via InputParser::ReadParameters(), allocates
     * memory to the pointers defined in the current class.
     * @return a boolean that is true if model initialization is successful
     */
    virtual bool InitializeModel();
        
    /**
     * @brief The post-update method for %FlavourWilsonCoefficient_DF2.
     * @details This method runs all the procedures that are need to be executed
     * after the model is successfully updated.
     * @return a boolean that is true if the execution is successful
     */
    virtual bool PostUpdate();
    
    /**
     * @brief A method to check if all the mandatory parameters for %FlavourWilsonCoefficient_DF2
     * have been provided in model initialization.
     * @param[in] DPars a map of the parameters that are being updated in the Monte Carlo run
     * (including parameters that are varied and those that are held constant)
     * @return a boolean that is true if the execution is successful
     */
    virtual bool CheckParameters(const std::map<std::string, double>& DPars);
    
    /**
     * @brief A get method to access the member reference of type %FlavourWilsonCoefficient_DF2Matching.
     * @return a reference to a %FlavourWilsonCoefficient_DF2Matching object
     */
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
     * @return the scale at which the NP Wilson coefficients in the \f$bd\f$ sector are defined
     */
    double getWCscale_bd() const
    {
        return WCscale_bd;
    }

    /**
     * 
     * @return the scale at which the NP Wilson coefficients in the \f$bs\f$ sector are defined
     */
    double getWCscale_bs() const
    {
        return WCscale_bs;
    }

    /**
     * 
     * @return the scale at which the NP Wilson coefficients in the \f$cu\f$ sector are defined
     */
    double getWCscale_c() const
    {
        return WCscale_c;
    }

    /**
     * 
     * @return the scale at which the NP Wilson coefficients in the \f$sd\f$ sector are defined
     */
    double getWCscale_s() const
    {
        return WCscale_s;
    }

protected: 
    
    /**
     * @brief A method to set the value of a parameter of %FlavourWilsonCoefficient_DF2.
     * @param[in] name name of a model parameter
     * @param[in] value the value to be assigned to the parameter specified by name
     */
    virtual void setParameter(const std::string, const double&);
    mutable Matching<FlavourWilsonCoefficient_DF2Matching,FlavourWilsonCoefficient_DF2> FWCM;///< The FlavourWilsonCoefficientMatching_DF2 object.

private:
    double reC1_s, reC2_s, reC3_s, reC4_s, reC5_s;///< The real parts of the  Wilson Coefficients.
    double reC1_c, reC2_c, reC3_c, reC4_c, reC5_c;///< The real parts of the  Wilson Coefficients.
    double reC1_bd, reC2_bd, reC3_bd, reC4_bd, reC5_bd;///< The real parts of the  Wilson Coefficients.
    double reC1_bs, reC2_bs, reC3_bs, reC4_bs, reC5_bs;///< The real parts of the  Wilson Coefficients.
    double imC1_s, imC2_s, imC3_s, imC4_s, imC5_s;///< The imaginary parts of the  Wilson Coefficients.
    double imC1_c, imC2_c, imC3_c, imC4_c, imC5_c;///< The imaginary parts of the  Wilson Coefficients.
    double imC1_bd, imC2_bd, imC3_bd, imC4_bd, imC5_bd;///< The imaginary parts of the  Wilson Coefficients.
    double imC1_bs, imC2_bs, imC3_bs, imC4_bs, imC5_bs;///< The imaginary parts of the  Wilson Coefficients.
    gslpp::vector<gslpp::complex> C_s, C_c, C_bd, C_bs;///< The complex Wilson Coefficients.
    double WCscale_s, WCscale_c, WCscale_bd, WCscale_bs;///< The scale of the  Wilson Coefficients.

};

#endif /* FLAVOURWILSONCOEFFICIENT_DF2_H */

