/* 
 * Copyright (C) 2026 HEPfit Collaboration
 *
 *
 * For the licensing terms see doc/COPYING.
 * 
 * Created on 18 September 2023, 16:23
 * 
 */

/* 
 * File:   NPSMEFTd6MFP.h
 * Author: silvest
 *
 * Created on 25 May 2026, 12:16
 */

#ifndef NPSMEFTD6MFP_H
#define NPSMEFTD6MFP_H

#include "NPSMEFTd6General.h"



class NPSMEFTd6MFP: public NPSMEFTd6General {
public:
    
    static const int NNPSMEFTd6MFPVars = 185+1;
    
    static std::string NPSMEFTd6MFPVars[NNPSMEFTd6MFPVars];
    
    NPSMEFTd6MFP();

    /**
     * @brief A method to initialize the model parameters.
     * @param[in] DPars a map of the parameters that are being updated in the Monte Carlo run
     * (including parameters that are varied and those that are held constant)
     * @return a boolean that is true if the execution is successful
     */
    virtual bool Init(const std::map<std::string, double>& DPars);
    
    
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
double CHq1_aar_LNP = 0.; ///< The dimension-6 operator coefficient \f$(C_{Hq}^{(1)})_{aa}(\Lambda_{\rm{NP}})\f$ (light-quark U(2) aggregate).
double CHq3_aar_LNP = 0.; ///< The dimension-6 operator coefficient \f$(C_{Hq}^{(3)})_{aa}(\Lambda_{\rm{NP}})\f$ (light-quark U(2) aggregate).
double Clq1_11aar_LNP = 0., Clq1_22aar_LNP = 0., Clq1_33aar_LNP = 0.; ///< \f$(C_{lq}^{(1)})_{llaa}\f$: per-lepton-gen, light-quark U(2) aggregate.
double Clq3_11aar_LNP = 0., Clq3_22aar_LNP = 0., Clq3_33aar_LNP = 0.; ///< \f$(C_{lq}^{(3)})_{llaa}\f$: per-lepton-gen, light-quark U(2) aggregate.
double Cqe_aa11r_LNP = 0., Cqe_aa22r_LNP = 0., Cqe_aa33r_LNP = 0.; ///< \f$(C_{qe})_{aall}\f$: light-quark U(2) aggregate, per-lepton-gen.
double Cqq1_aabbr_LNP = 0., Cqq1_abbar_LNP = 0., Cqq1_a33ar_LNP = 0., Cqq1_aa33r_LNP = 0.; ///< \f$(C_{qq}^{(1)})\f$ light-quark aggregates.
double Cqq3_aabbr_LNP = 0., Cqq3_abbar_LNP = 0., Cqq3_a33ar_LNP = 0., Cqq3_aa33r_LNP = 0.; ///< \f$(C_{qq}^{(3)})\f$ light-quark aggregates.
double Cqu1_aa11r_LNP = 0., Cqu1_aa22r_LNP = 0., Cqu1_aa33r_LNP = 0.; ///< \f$(C_{qu}^{(1)})_{aaii}\f$: light-quark U(2) aggregate, per-u_R-gen.
double Cqu8_aa11r_LNP = 0., Cqu8_aa22r_LNP = 0., Cqu8_aa33r_LNP = 0.; ///< \f$(C_{qu}^{(8)})_{aaii}\f$: light-quark U(2) aggregate, per-u_R-gen.
double Cqd1_aa11r_LNP = 0., Cqd1_aa22r_LNP = 0., Cqd1_aa33r_LNP = 0.; ///< \f$(C_{qd}^{(1)})_{aaii}\f$: light-quark U(2) aggregate, per-d_R-gen.
double Cqd8_aa11r_LNP = 0., Cqd8_aa22r_LNP = 0., Cqd8_aa33r_LNP = 0.; ///< \f$(C_{qd}^{(8)})_{aaii}\f$: light-quark U(2) aggregate, per-d_R-gen.
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


#endif /* NPSMEFTD6MFP_H */
