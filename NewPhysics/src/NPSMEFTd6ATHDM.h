/*
 * Click nbfs://nbhost/SystemFileSystem/Templates/Licenses/license-default.txt to change this license
 * Click nbfs://nbhost/SystemFileSystem/Templates/cppFiles/file.h to edit this template
 */

/* 
 * File:   NPSMEFTd6ATHDM.h
 * Author: miralles
 *
 * Created on April 17, 2026, 10:17 AM
 */

#ifndef NPSMEFTD6ATHDM_H
#define NPSMEFTD6ATHDM_H



#include "NPSMEFTd6MFV.h"


class NPSMEFTd6ATHDM: public NPSMEFTd6MFV {
public:
    
    static const int NNPSMEFTd6ATHDMVars = 9+1;
    
    static const std::string NPSMEFTd6ATHDMVars[NNPSMEFTd6ATHDMVars];
    
    NPSMEFTd6ATHDM();

    /**
     * @brief The post-update method for %NPSMEFTd6General.
     * @details This method runs all the procedures that are need to be executed
     * after the model is successfully updated.
     * @return a boolean that is true if the execution is successful
     */
    virtual bool PostUpdate();
    
    
protected:
    
    // These Higgs parameters are not exactly those of Eq. 4 of [2401.12279] but a small rotation of them
    // These would be exactly the parameters of Eq. 4 of [2401.12279] in the basis in which Y3=0, which
    // is the relevant basis to perform the matchings. Therefore, these wouldn't be the parameters or the
    // Higgs basis but a small rotation of the Higgs basis in which Y3=0. The relation of these parameters 
    // and thos of the Higgs basis, expanding at order 1/Y2 would be:
    // Y1 = Y1Higgs - Y3Higgs^2/(Y2Higgs−Y1Higgs)
    // Y2 = Y2Higgs - Y3Higgs^2/(Y2Higgs−Y1Higgs)
    // Y3 = 0
    // Z1 = Z1Higgs − 4*(Y3Higgs/(Y2Higgs−Y1Higgs))*Z6Higgs
    // Z2 = Z2Higgs + 4*(Y3Higgs/(Y2Higgs−Y1Higgs))*Z7Higgs
    // Z3 = Z3Higgs + 2*(Y3Higgs/(Y2Higgs−Y1Higgs))*(Z6Higgs-Z7Higgs)
    // Z4 = Z4Higgs + 2*(Y3Higgs/(Y2Higgs−Y1Higgs))*(Z6Higgs-Z7Higgs)
    // Z5 = Z5Higgs + 2*(Y3Higgs/(Y2Higgs−Y1Higgs))*(Z6Higgs-Z7Higgs)
    // Z6 = Z6Higgs + 4*(Y3Higgs/(Y2Higgs−Y1Higgs))*(Z1Higgs-Z3Higgs-Z4Higgs-Z5Higgs)
    // Z7 = Z7Higgs + 4*(Y3Higgs/(Y2Higgs−Y1Higgs))*(Z3Higgs+Z4Higgs+Z5Higgs-Z1Higgs)
    // Also, we have not included Y1 and Z1 as parameters since they are related to the SM parameters:
    // Z1 = 2\lambda
    // Y1 = -\mu^2
    // being \lambda and \mu the parameters of the SM Higgs potential
    
    double sqY2 = 0.; ///< The square root of parameter of the Higgs potential Y2.
    
    double Z2 = 0.; ///< The parameter of the Higgs potential.
    double Z3 = 0.; ///< The parameter of the Higgs potential.
    double Z4 = 0.; ///< The parameter of the Higgs potential.
    double Z5 = 0.; ///< The parameter of the Higgs potential.
    double Z6 = 0.; ///< The parameter of the Higgs potential.
    double Z7 = 0.; ///< The parameter of the Higgs potential.
    
    // The yukawa proportionality parameters are also defined in the basis in which Y3=0, so the 
    // matching basis. Furthermore, they're not defined exactly as in Eq. 3 of [2401.12279] but 
    // we reabsorbe 1/tan\beta inside them such that the new Yukawa have the form
    // \eta_u Y_u *\bar{q}_L \tilde{H}_2 u_R and so on.
    
    double etau = 0.; ///< The Yukawa proprtionality parameter of the new doublet
    double etad = 0.; ///< The Yukawa proprtionality parameter of the new doublet
    double etae = 0.; ///< The Yukawa proprtionality parameter of the new doublet
    
    
    virtual void setParameter(const std::string name, const double& value);
    
    
    /**
     * @brief An auxiliary method to set the WC of the MFV class
     */
    void setNPSMEFTd6MFVParameters();
    
private:

};





#endif /* ATHDMMATCHINGNPSMEFTD6MFV_H */

