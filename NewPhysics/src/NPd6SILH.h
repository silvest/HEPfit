/* Copyright (C) 2025 HEPfit Collaboration
* 
* For the licensing terms see doc/COPYING.
*
* Created on Sun 27 Jul 2025 17:40:34
*
*/

#ifndef NPD6SILH_H
#define NPD6SILH_H

#include "NPSMEFTd6General.h"

class NPd6SILH:public NPSMEFTd6General {
public:

    static const int NNPd6SILHVars = 28+1;

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


protected:

    double cH = 0.0;
    double cT = 0.0;
    double c6 = 0.0;
    double cB = 0.0;
    double cW = 0.0;
    double c2B = 0.0;
    double c2W = 0.0;
    double c2G = 0.0;
    double c3W = 0.0;
    double c3G = 0.0;
    double cHW = 0.0;
    double cHB = 0.0;
    double cgam = 0.0;
    double cg = 0.0;
    double ctD = 0.0;
    double cqD1 = 0.0;
    double cqD3 = 0.0;
    double cqq1 = 0.0;
    double cqq3 = 0.0;
    double cqt1 = 0.0;
    double cqt8 = 0.0;
    double ctt = 0.0;
    double ctG = 0.0;
    double ctB = 0.0;
    double ctW = 0.0;
    double cu = 0.0;
    double cd = 0.0;
    double ce = 0.0;


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
