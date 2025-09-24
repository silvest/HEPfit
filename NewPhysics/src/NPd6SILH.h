/* Copyright (C) 2025 HEPfit Collaboration
* 
* For the licensing terms see doc/COPYING.
*
* Created on Tue 29 Jul 2025 16:05:01
*
*/

#ifndef NPD6SILH_H
#define NPD6SILH_H

#include "NPSMEFTd6General.h"

class NPd6SILH:public NPSMEFTd6General {
public:

    static const int NNPd6SILHVars = 31+1;

    static std::string NPd6SILHVars[NNPd6SILHVars];

    NPd6SILH();


    /*
     * @brief The post-update method for %NPd6SILH.
     * @details This method runs all the procedures that are need to be executed
     * after the model is successfully updated.
     * @return a boolean that is true if the execution is successful
     */
    virtual bool PostUpdate();


    /*
     * @brief A method to set a flag of %NPd6SILH.
     * @param[in] name name of a model flag
     * @param[in] value the boolean to be assigned to the flag specified by name
     * @return a boolean that is true if the execution is successful
     */
    virtual bool setFlag(const std::string name, const bool value);
    
      
    /**
     * @brief The oblique parameter \f$S\f$.
     * (Simplified implementation. Contribution only from @f$O_{HWB}@f$.)
     * @return the value of @f$S@f$
     */
    virtual const double obliqueS() const;

    /**
     * @brief The oblique parameter \f$T\f$.
     * (Simplified implementation. Contribution only from @f$O_{HD}@f$.)
     * @return the value of @f$T@f$
     */
    virtual const double obliqueT() const;

    /**
     * @brief The oblique parameter \f$W\f$.
     * (Simplified implementation. Contribution only from @f$O_{2W}@f$.)
     * @return the value of @f$W@f$
     */
    virtual const double obliqueW() const;

    /**
     * @brief The oblique parameter \f$Y\f$.
     * (Simplified implementation. Contribution only from @f$O_{2B}@f$.)
     * @return the value of @f$Y@f$
     */
    virtual const double obliqueY() const;
    
    /**
     * @brief Auxiliary observable AuxObs_NP30
     * @return AuxObs_NP30
     */
    virtual const double AuxObs_NP30() const;


protected:

    double cH_LNP = 0.0;
    double cT_LNP = 0.0;
    double c6_LNP = 0.0;
    double cB_LNP = 0.0;
    double cW_LNP = 0.0;
    double c2B_LNP = 0.0;
    double c2W_LNP = 0.0;
    double c2G_LNP = 0.0;
    double c3W_LNP = 0.0;
    double c3G_LNP = 0.0;
    double cHW_LNP = 0.0;
    double cHB_LNP = 0.0;
    double cgam_LNP = 0.0;
    double cg_LNP = 0.0;
    double cHq1_LNP = 0.0;
    double cHq3_LNP = 0.0;
    double cHt_LNP = 0.0;
    double ctD_LNP = 0.0;
    double cqD1_LNP = 0.0;
    double cqD3_LNP = 0.0;
    double cqq1_LNP = 0.0;
    double cqq3_LNP = 0.0;
    double cqt1_LNP = 0.0;
    double cqt8_LNP = 0.0;
    double ctt_LNP = 0.0;
    double ctG_LNP = 0.0;
    double ctB_LNP = 0.0;
    double ctW_LNP = 0.0;
    double cu_LNP = 0.0;
    double cd_LNP = 0.0;
    double ce_LNP = 0.0;
    
    double v2LambdaNP2; ///< The ratio between the EW vev and the new physics scale, squared \f$v^2/\Lambda^2\f$.

    virtual void setParameter(const std::string name, const double& value);


/*
* @brief An auxiliary method to set the WC of the general class
*/
    void setNPSMEFTd6GeneralParameters();


private:

    double g1UV = 0., g2UV = 0., g3UV = 0., lambdaHUV = 0.; ///< SM couplings at the UV scale. 
    double g1UV2 = 0., g2UV2 = 0., g3UV2 = 0.; ///< SM gauge couplings at the UV scale (squared). 
    gslpp::matrix<gslpp::complex> YuUV, YuUVhc, YdUV, YdUVhc, YeUV, YeUVhc; ///< Yukawa matrices and h.c. at the UV scale. 

    bool FlagRGEci; ///< A boolean for the model flag %RGEci , to include RGE effects. (Overwrittes NPSMEFTGeneral.)
};


#endif /* NPD6SILH_H */
