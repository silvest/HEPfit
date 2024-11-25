
/* 
 * File:   NPSMEFTd6U2qU1le.h
 * Author: de Blas
 *
 * Created on November 25, 2024
 */

#ifndef NPSMEFTD6U2QU1LE_H
#define NPSMEFTD6U2QU1LE_H

#include "NPSMEFTd6General.h"



class NPSMEFTd6U2qU1le: public NPSMEFTd6General {
public:
    
    static const int NNPSMEFTd6U2qU1leVars = 168+1;
    
    static std::string NPSMEFTd6U2qU1leVars[NNPSMEFTd6U2qU1leVars];
    
    NPSMEFTd6U2qU1le();

    
    
    /**
     * @brief The post-update method for %NPSMEFTd6General.
     * @details This method runs all the procedures that are need to be executed
     * after the model is successfully updated.
     * @return a boolean that is true if the execution is successful
     */
    virtual bool PostUpdate();
    
    
protected:

//If we define here the WC which have the same name as those from NPSMEFTd6General the code fails
//since those variables in NPSMEFTd6General (which are the ones used in the observables) will not
//be properly assigned with the right value  
double CHq1_aar_LNP = 0.;//, CHq1_33r_LNP = 0.; ///< The dimension-6 operator coefficient \f$(C_{Hq}^{(1)})_{ij}(\Lambda_{\rm{NP}})\f$ (Real part and pure real operator).
double CHq3_aar_LNP = 0.;//, CHq3_33r_LNP = 0.; ///< The dimension-6 operator coefficient \f$(C_{Hq}^{(3)})_{ij}(\Lambda_{\rm{NP}})\f$ (Real part and pure real operator).
double CHu_aar_LNP = 0.;//, CHu_33r_LNP = 0.; ///< The dimension-6 operator coefficient \f$(C_{Hu})_{ij}(\Lambda_{\rm{NP}})\f$ (Real part and pure real operator).
double CHd_aar_LNP = 0.;//, CHd_33r_LNP = 0.; ///< The dimension-6 operator coefficient \f$(C_{Hd})_{ij}(\Lambda_{\rm{NP}})\f$ (Real part and pure real operator).

double Clq1_11aar_LNP = 0., Clq1_22aar_LNP = 0., Clq1_33aar_LNP = 0.; ///< The dimension-6 operator coefficient \f$(C_{lq}^{(1)})_{ijkm}(\Lambda_{\rm{NP}})\f$ (Real part and pure real operator).
double Clq3_11aar_LNP = 0., Clq3_22aar_LNP = 0., Clq3_33aar_LNP = 0.; ///< The dimension-6 operator coefficient \f$(C_{lq}^{(3)})_{ijkm}(\Lambda_{\rm{NP}})\f$ (Real part and pure real operator).
double Ceu_11aar_LNP = 0., Ceu_22aar_LNP = 0., Ceu_33aar_LNP = 0.; ///< The dimension-6 operator coefficient \f$(C_{eu})_{ijkm}(\Lambda_{\rm{NP}})\f$ (Real part and pure real operator).
double Ced_11aar_LNP = 0., Ced_22aar_LNP = 0., Ced_33aar_LNP = 0.; ///< The dimension-6 operator coefficient \f$(C_{ed})_{ijkm}(\Lambda_{\rm{NP}})\f$ (Real part and pure real operator).
double Clu_11aar_LNP = 0., Clu_22aar_LNP = 0., Clu_33aar_LNP = 0.; ///< The dimension-6 operator coefficient \f$(C_{lu})_{ijkm}(\Lambda_{\rm{NP}})\f$ (Real part and pure real operator).
double Cld_11aar_LNP = 0., Cld_22aar_LNP = 0., Cld_33aar_LNP = 0.; ///< The dimension-6 operator coefficient \f$(C_{ld})_{ijkm}(\Lambda_{\rm{NP}})\f$ (Real part and pure real operator).
double Cqe_aa11r_LNP = 0., Cqe_aa22r_LNP = 0., Cqe_aa33r_LNP = 0.; ///< The dimension-6 operator coefficient \f$(C_{qe})_{ijkm}(\Lambda_{\rm{NP}})\f$ (Real part and pure real operator).
double Cqq1_aabbr_LNP = 0., Cqq1_abbar_LNP = 0., Cqq1_aa33r_LNP = 0., Cqq1_a33ar_LNP = 0.; ///< The dimension-6 operator coefficient \f$(C_{ll})_{ijkm}(\Lambda_{\rm{NP}})\f$ (Real part and pure real operator).
double Cqq3_aabbr_LNP = 0., Cqq3_abbar_LNP = 0., Cqq3_aa33r_LNP = 0., Cqq3_a33ar_LNP = 0.; ///< The dimension-6 operator coefficient \f$(C_{ll})_{ijkm}(\Lambda_{\rm{NP}})\f$ (Real part and pure real operator).
double Cuu_aabbr_LNP = 0., Cuu_abbar_LNP = 0., Cuu_aa33r_LNP = 0., Cuu_a33ar_LNP = 0.; ///< The dimension-6 operator coefficient \f$(C_{ll})_{ijkm}(\Lambda_{\rm{NP}})\f$ (Real part and pure real operator).
double Cdd_aabbr_LNP = 0., Cdd_abbar_LNP = 0., Cdd_aa33r_LNP = 0., Cdd_a33ar_LNP = 0.; ///< The dimension-6 operator coefficient \f$(C_{ll})_{ijkm}(\Lambda_{\rm{NP}})\f$ (Real part and pure real operator).
double Cud1_aabbr_LNP = 0., Cud1_aa33r_LNP = 0., Cud1_33aar_LNP = 0.; ///< The dimension-6 operator coefficient \f$(C_{ud}^{(1)})_{ijkm}(\Lambda_{\rm{NP}})\f$ (Real part and pure real operator).
double Cud8_aabbr_LNP = 0., Cud8_aa33r_LNP = 0., Cud8_33aar_LNP = 0.; ///< The dimension-6 operator coefficient \f$(C_{ud}^{(8)})_{ijkm}(\Lambda_{\rm{NP}})\f$ (Real part and pure real operator).
double Cqu1_aabbr_LNP = 0., Cqu1_aa33r_LNP = 0., Cqu1_33aar_LNP = 0.; ///< The dimension-6 operator coefficient \f$(C_{qu}^{(1)})_{ijkm}(\Lambda_{\rm{NP}})\f$ (Real part and pure real operator).
double Cqu8_aabbr_LNP = 0., Cqu8_aa33r_LNP = 0., Cqu8_33aar_LNP = 0.; ///< The dimension-6 operator coefficient \f$(C_{qu}^{(8)})_{ijkm}(\Lambda_{\rm{NP}})\f$ (Real part and pure real operator).
double Cqd1_aabbr_LNP = 0., Cqd1_aa33r_LNP = 0., Cqd1_33aar_LNP = 0.; ///< The dimension-6 operator coefficient \f$(C_{qd}^{(1)})_{ijkm}(\Lambda_{\rm{NP}})\f$ (Real part and pure real operator).
double Cqd8_aabbr_LNP = 0., Cqd8_aa33r_LNP = 0., Cqd8_33aar_LNP = 0.; ///< The dimension-6 operator coefficient \f$(C_{qd}^{(8)})_{ijkm}(\Lambda_{\rm{NP}})\f$ (Real part and pure real operator).


virtual void setParameter(const std::string name, const double& value);


/**
* @brief An auxiliary method to set the WC of the general class
*/
void setNPSMEFTd6GeneralParameters();

private:

};


#endif /* NPSMEFTD6U2QU1LE_H */


