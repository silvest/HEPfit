/*
 * Click nbfs://nbhost/SystemFileSystem/Templates/Licenses/license-default.txt to change this license
 * Click nbfs://nbhost/SystemFileSystem/Templates/cppFiles/file_header.h to edit this template
 */

/* 
 * File:   NPSMEFTd6U2.h
 * Author: miralles
 *
 * Created on 18 September 2023, 16:23
 */

#ifndef NPSMEFTD6U2_H
#define NPSMEFTD6U2_H

#include "NPSMEFTd6General.h"



class NPSMEFTd6U2: public NPSMEFTd6General {
public:
    
    static const int NNPSMEFTd6U2Vars = 124+1;
    
    static std::string NPSMEFTd6U2Vars[NNPSMEFTd6U2Vars];
    
    NPSMEFTd6U2();

    
    
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
//double CG_LNP = 0.; ///< The dimension-6 operator coefficient \f$C_{G}(\Lambda_{\rm{NP}})\f$.
//double CW_LNP = 0.; ///< The dimension-6 operator coefficient \f$C_{W}(\Lambda_{\rm{NP}})\f$.
//double CHG_LNP = 0.; ///< The dimension-6 operator coefficient \f$C_{HG}(\Lambda_{\rm{NP}})\f$.
//double CHW_LNP = 0.; ///< The dimension-6 operator coefficient \f$C_{HW}(\Lambda_{\rm{NP}})\f$.
//double CHB_LNP = 0.; ///< The dimension-6 operator coefficient \f$C_{HB}(\Lambda_{\rm{NP}})\f$.
//double CHWB_LNP = 0.; ///< The dimension-6 operator coefficient \f$C_{HWB}(\Lambda_{\rm{NP}})\f$.
//double CHD_LNP = 0.; ///< The dimension-6 operator coefficient \f$C_{HD}(\Lambda_{\rm{NP}})\f$.
//double CHbox_LNP = 0.; ///< The dimension-6 operator coefficient \f$C_{H\Box}(\Lambda_{\rm{NP}})\f$.
//double CH_LNP = 0.; ///< The dimension-6 operator coefficient \f$C_{H}(\Lambda_{\rm{NP}})\f$.
////double CGtilde_LNP = 0.; ///< The dimension-6 operator coefficient \f$C_{\tilde{G}}(\Lambda_{\rm{NP}})\f$.
////double CWtilde_LNP = 0.; ///< The dimension-6 operator coefficient \f$C_{\tilde{W}}(\Lambda_{\rm{NP}})\f$.
////double CHGtilde_LNP = 0.; ///< The dimension-6 operator coefficient \f$C_{H\tilde{G}}(\Lambda_{\rm{NP}})\f$.
////double CHWtilde_LNP = 0.; ///< The dimension-6 operator coefficient \f$C_{H\tilde{W}}(\Lambda_{\rm{NP}})\f$.
////double CHBtilde_LNP = 0.; ///< The dimension-6 operator coefficient \f$C_{H\tilde{B}}(\Lambda_{\rm{NP}})\f$.
////double CHWtildeB_LNP = 0.; ///< The dimension-6 operator coefficient \f$C_{H\tilde{W}B}(\Lambda_{\rm{NP}})\f$.
double CHl1_aar_LNP = 0.;//, CHl1_33r_LNP = 0.; ///< The dimension-6 operator coefficient \f$(C_{Hl}^{(1)})_{ij}(\Lambda_{\rm{NP}})\f$ (Real part and pure real operator).
double CHl3_aar_LNP = 0.;//, CHl3_33r_LNP = 0.; ///< The dimension-6 operator coefficient \f$(C_{Hl}^{(3)})_{ij}(\Lambda_{\rm{NP}})\f$ (Real part and pure real operator).
double CHe_aar_LNP = 0.;//, CHe_33r_LNP = 0.; ///< The dimension-6 operator coefficient \f$(C_{He})_{ij}(\Lambda_{\rm{NP}})\f$ (Real part and pure real operator).
double CHq1_aar_LNP = 0.;//, CHq1_33r_LNP = 0.; ///< The dimension-6 operator coefficient \f$(C_{Hq}^{(1)})_{ij}(\Lambda_{\rm{NP}})\f$ (Real part and pure real operator).
double CHq3_aar_LNP = 0.;//, CHq3_33r_LNP = 0.; ///< The dimension-6 operator coefficient \f$(C_{Hq}^{(3)})_{ij}(\Lambda_{\rm{NP}})\f$ (Real part and pure real operator).
double CHu_aar_LNP = 0.;//, CHu_33r_LNP = 0.; ///< The dimension-6 operator coefficient \f$(C_{Hu})_{ij}(\Lambda_{\rm{NP}})\f$ (Real part and pure real operator).
double CHd_aar_LNP = 0.;//, CHd_33r_LNP = 0.; ///< The dimension-6 operator coefficient \f$(C_{Hd})_{ij}(\Lambda_{\rm{NP}})\f$ (Real part and pure real operator).
//double CHud_33r_LNP = 0.; ///< The dimension-6 operator coefficient \f$(C_{Hud})_{ij}(\Lambda_{\rm{NP}})\f$ (Real part and pure real operator).
//double CeH_33r_LNP = 0.; ///< The dimension-6 operator coefficient \f$(C_{eH})_{ij}(\Lambda_{\rm{NP}})\f$ (Real part and pure real operator).
//double CuH_33r_LNP = 0.; ///< The dimension-6 operator coefficient \f$(C_{uH})_{ij}(\Lambda_{\rm{NP}})\f$ (Real part and pure real operator).
//double CdH_33r_LNP = 0.; ///< The dimension-6 operator coefficient \f$(C_{dH})_{ij}(\Lambda_{\rm{NP}})\f$ (Real part and pure real operator).
//double CuG_33r_LNP = 0.; ///< The dimension-6 operator coefficient \f$(C_{uG})_{ij}(\Lambda_{\rm{NP}})\f$ (Real part and pure real operator).
//double CuW_33r_LNP = 0.; ///< The dimension-6 operator coefficient \f$(C_{uW})_{ij}(\Lambda_{\rm{NP}})\f$ (Real part and pure real operator).
//double CuB_33r_LNP = 0.; ///< The dimension-6 operator coefficient \f$(C_{uB})_{ij}(\Lambda_{\rm{NP}})\f$ (Real part and pure real operator).
//double CdG_33r_LNP = 0.; ///< The dimension-6 operator coefficient \f$(C_{dG})_{ij}(\Lambda_{\rm{NP}})\f$ (Real part and pure real operator).
//double CdW_33r_LNP = 0.; ///< The dimension-6 operator coefficient \f$(C_{dW})_{ij}(\Lambda_{\rm{NP}})\f$ (Real part and pure real operator).
//double CdB_33r_LNP = 0.; ///< The dimension-6 operator coefficient \f$(C_{dB})_{ij}(\Lambda_{\rm{NP}})\f$ (Real part and pure real operator).
//double CeW_33r_LNP = 0.; ///< The dimension-6 operator coefficient \f$(C_{eW})_{ij}(\Lambda_{\rm{NP}})\f$ (Real part and pure real operator).
//double CeB_33r_LNP = 0.; ///< The dimension-6 operator coefficient \f$(C_{eB})_{ij}(\Lambda_{\rm{NP}})\f$ (Real part and pure real operator).


double Cll_aabbr_LNP = 0., Cll_abbar_LNP = 0., Cll_aa33r_LNP = 0., Cll_a33ar_LNP = 0.;//, Cll_3333r_LNP = 0.; ///< The dimension-6 operator coefficient \f$(C_{ll})_{ijkm}(\Lambda_{\rm{NP}})\f$ (Real part and pure real operator).
double Clq1_aabbr_LNP = 0., Clq1_aa33r_LNP = 0., Clq1_33aar_LNP = 0.;//, Clq1_3333r_LNP = 0.; ///< The dimension-6 operator coefficient \f$(C_{lq}^{(1)})_{ijkm}(\Lambda_{\rm{NP}})\f$ (Real part and pure real operator).
double Clq3_aabbr_LNP = 0., Clq3_aa33r_LNP = 0., Clq3_33aar_LNP = 0.;//, Clq3_3333r_LNP = 0.; ///< The dimension-6 operator coefficient \f$(C_{lq}^{(3)})_{ijkm}(\Lambda_{\rm{NP}})\f$ (Real part and pure real operator).
double Cee_aabbr_LNP = 0., Cee_aa33r_LNP = 0.;//, Cee_3333r_LNP = 0.; ///< The dimension-6 operator coefficient \f$(C_{ee})_{ijkm}(\Lambda_{\rm{NP}})\f$ (Real part and pure real operator).
double Ceu_aabbr_LNP = 0., Ceu_aa33r_LNP = 0., Ceu_33aar_LNP = 0.;//, Ceu_3333r_LNP = 0.; ///< The dimension-6 operator coefficient \f$(C_{eu})_{ijkm}(\Lambda_{\rm{NP}})\f$ (Real part and pure real operator).
double Ced_aabbr_LNP = 0., Ced_aa33r_LNP = 0., Ced_33aar_LNP = 0.;//, Ced_3333r_LNP = 0.; ///< The dimension-6 operator coefficient \f$(C_{ed})_{ijkm}(\Lambda_{\rm{NP}})\f$ (Real part and pure real operator).
double Cle_aabbr_LNP = 0., Cle_aa33r_LNP = 0., Cle_33aar_LNP = 0.;//, Cle_3333r_LNP = 0.; ///< The dimension-6 operator coefficient \f$(C_{le})_{ijkm}(\Lambda_{\rm{NP}})\f$ (Real part and pure real operator).
double Clu_aabbr_LNP = 0., Clu_aa33r_LNP = 0., Clu_33aar_LNP = 0.;//, Clu_3333r_LNP = 0.; ///< The dimension-6 operator coefficient \f$(C_{lu})_{ijkm}(\Lambda_{\rm{NP}})\f$ (Real part and pure real operator).
double Cld_aabbr_LNP = 0., Cld_aa33r_LNP = 0., Cld_33aar_LNP = 0.;//, Cld_3333r_LNP = 0.; ///< The dimension-6 operator coefficient \f$(C_{ld})_{ijkm}(\Lambda_{\rm{NP}})\f$ (Real part and pure real operator).
double Cqe_aabbr_LNP = 0., Cqe_aa33r_LNP = 0., Cqe_33aar_LNP = 0.;//, Cqe_3333r_LNP = 0.; ///< The dimension-6 operator coefficient \f$(C_{qe})_{ijkm}(\Lambda_{\rm{NP}})\f$ (Real part and pure real operator).
//double Cledq_3333r_LNP = 0.; ///< The dimension-6 operator coefficient \f$(C_{ledq})_{ijkm}(\Lambda_{\rm{NP}})\f$ (Real part and pure real operator).
double Cqq1_aabbr_LNP = 0., Cqq1_abbar_LNP = 0., Cqq1_aa33r_LNP = 0., Cqq1_a33ar_LNP = 0.;//, Cqq1_3333r_LNP = 0.; ///< The dimension-6 operator coefficient \f$(C_{ll})_{ijkm}(\Lambda_{\rm{NP}})\f$ (Real part and pure real operator).
double Cqq3_aabbr_LNP = 0., Cqq3_abbar_LNP = 0., Cqq3_aa33r_LNP = 0., Cqq3_a33ar_LNP = 0.;//, Cqq3_3333r_LNP = 0.; ///< The dimension-6 operator coefficient \f$(C_{ll})_{ijkm}(\Lambda_{\rm{NP}})\f$ (Real part and pure real operator).
double Cuu_aabbr_LNP = 0., Cuu_abbar_LNP = 0., Cuu_aa33r_LNP = 0., Cuu_a33ar_LNP = 0.;//, Cuu_3333r_LNP = 0.; ///< The dimension-6 operator coefficient \f$(C_{ll})_{ijkm}(\Lambda_{\rm{NP}})\f$ (Real part and pure real operator).
double Cdd_aabbr_LNP = 0., Cdd_abbar_LNP = 0., Cdd_aa33r_LNP = 0., Cdd_a33ar_LNP = 0.;//, Cdd_3333r_LNP = 0.; ///< The dimension-6 operator coefficient \f$(C_{ll})_{ijkm}(\Lambda_{\rm{NP}})\f$ (Real part and pure real operator).
double Cud1_aabbr_LNP = 0., Cud1_aa33r_LNP = 0., Cud1_33aar_LNP = 0.;//, Cud1_3333r_LNP = 0.; ///< The dimension-6 operator coefficient \f$(C_{ud}^{(1)})_{ijkm}(\Lambda_{\rm{NP}})\f$ (Real part and pure real operator).
double Cud8_aabbr_LNP = 0., Cud8_aa33r_LNP = 0., Cud8_33aar_LNP = 0.;//, Cud8_3333r_LNP = 0.; ///< The dimension-6 operator coefficient \f$(C_{ud}^{(8)})_{ijkm}(\Lambda_{\rm{NP}})\f$ (Real part and pure real operator).
double Cqu1_aabbr_LNP = 0., Cqu1_aa33r_LNP = 0., Cqu1_33aar_LNP = 0.;//, Cqu1_3333r_LNP = 0.; ///< The dimension-6 operator coefficient \f$(C_{qu}^{(1)})_{ijkm}(\Lambda_{\rm{NP}})\f$ (Real part and pure real operator).
double Cqu8_aabbr_LNP = 0., Cqu8_aa33r_LNP = 0., Cqu8_33aar_LNP = 0.;//, Cqu8_3333r_LNP = 0.; ///< The dimension-6 operator coefficient \f$(C_{qu}^{(8)})_{ijkm}(\Lambda_{\rm{NP}})\f$ (Real part and pure real operator).
double Cqd1_aabbr_LNP = 0., Cqd1_aa33r_LNP = 0., Cqd1_33aar_LNP = 0.;//, Cqd1_3333r_LNP = 0.; ///< The dimension-6 operator coefficient \f$(C_{qd}^{(1)})_{ijkm}(\Lambda_{\rm{NP}})\f$ (Real part and pure real operator).
double Cqd8_aabbr_LNP = 0., Cqd8_aa33r_LNP = 0., Cqd8_33aar_LNP = 0.;//, Cqd8_3333r_LNP = 0.; ///< The dimension-6 operator coefficient \f$(C_{qd}^{(8)})_{ijkm}(\Lambda_{\rm{NP}})\f$ (Real part and pure real operator).
//double Cquqd1_3333r_LNP = 0.; ///< The dimension-6 operator coefficient \f$(C_{quqd}^{(1)})_{ijkm}(\Lambda_{\rm{NP}})\f$ (Real part and pure real operator).
//double Cquqd8_3333r_LNP = 0.; ///< The dimension-6 operator coefficient \f$(C_{quqd}^{(8)})_{ijkm}(\Lambda_{\rm{NP}})\f$ (Real part and pure real operator).
//double Clequ1_3333r_LNP = 0.; ///< The dimension-6 operator coefficient \f$(C_{lequ}^{(1)})_{ijkm}(\Lambda_{\rm{NP}})\f$ (Real part and pure real operator).
//double Clequ3_3333r_LNP = 0.; ///< The dimension-6 operator coefficient \f$(C_{lequ}^{(3)})_{ijkm}(\Lambda_{\rm{NP}})\f$ (Real part and pure real operator).



virtual void setParameter(const std::string name, const double& value);


/**
* @brief An auxiliary method to set the WC of the general class
*/
void setNPSMEFTd6GeneralParameters();

private:

};


#endif /* NPSMEFTD6U2_H */


