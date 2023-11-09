/*
 * Click nbfs://nbhost/SystemFileSystem/Templates/Licenses/license-default.txt to change this license
 * Click nbfs://nbhost/SystemFileSystem/Templates/cppFiles/class.h to edit this template
 */

/* 
 * File:   NPSMEFTd6U3.h
 * Author: silvest
 *
 * Created on 5 settembre 2023, 15.54
 */

#ifndef NPSMEFTD6MFVLONOCKM_H
#define NPSMEFTD6MFVLONOCKM_H

#include "NPSMEFTd6General.h"

class NPSMEFTd6U3: public NPSMEFTd6General {
public:
    
    static const int NNPSMEFTd6U3Vars = 41+1;
    
    static const std::string NPSMEFTd6U3Vars[NNPSMEFTd6U3Vars];
    
    NPSMEFTd6U3();

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
    //double CG_LNP = 0.; ///< The dimension-6 operator coefficient \f$C_{G}\f$.
    //double CW_LNP = 0.; ///< The dimension-6 operator coefficient \f$C_{W}\f$.
    //double CHG_LNP = 0.; ///< The dimension-6 operator coefficient \f$C_{HG}\f$.
    //double CHW_LNP = 0.; ///< The dimension-6 operator coefficient \f$C_{HW}\f$.
    //double CHB_LNP = 0.; ///< The dimension-6 operator coefficient \f$C_{HB}\f$.
    //double CHWB_LNP = 0.; ///< The dimension-6 operator coefficient \f$C_{HWB}\f$.
    //double CHD_LNP = 0.; ///< The dimension-6 operator coefficient \f$C_{HD}\f$.
    //double CHbox_LNP = 0.; ///< The dimension-6 operator coefficient \f$C_{H\Box}\f$.
    //double CH_LNP = 0.; ///< The dimension-6 operator coefficient \f$C_{H}\f$.
    
    
    double CHl1_LNP = 0.; ///< The dimension-6 operator coefficient \f$(C_{HL}^{(1)})_{ij}\f$.
    double CHl3_LNP = 0.; ///< The dimension-6 operator coefficient \f$(C_{HL}^{(3)})_{ij}\f$.
    double CHe_LNP = 0.; ///< The dimension-6 operator coefficient \f$(C_{He})_{ij}\f$.
    double CHq1_LNP = 0.; ///< The dimension-6 operator coefficient \f$(C_{HQ}^{(1)})_{ij}\f$.
    double CHq3_LNP = 0.; ///< The dimension-6 operator coefficient \f$(C_{HQ}^{(3)})_{ij}\f$.
    double CHu_LNP = 0.; ///< The dimension-6 operator coefficient \f$(C_{Hu})_{ij}\f$.
    double CHd_LNP = 0.; ///< The dimension-6 operator coefficient \f$(C_{Hd})_{ij}\f$.
    
    double Cll_aabb_LNP = 0.; ///< The dimension-6 operator coefficient \f$(C_{ll})_{ijkm}\f$ (Real part and pure real operator).
    double Cll_abba_LNP = 0.; ///< The dimension-6 operator coefficient \f$(C_{ll})_{ijkm}\f$ (Real part and pure real operator).

    double Clq1_LNP = 0.; ///< The dimension-6 operator coefficient \f$(C_{lq}^{(1)})_{ijkm}\f$ (Imaginary part).
    double Clq3_LNP = 0.; ///< The dimension-6 operator coefficient \f$(C_{lq}^{(3)})_{ijkm}\f$ (Real part and pure real operator).
    double Cee_LNP = 0.; ///< The dimension-6 operator coefficient \f$(C_{ee})_{ijkm}\f$ (Imaginary part).
    double Ceu_LNP = 0.; ///< The dimension-6 operator coefficient \f$(C_{eu})_{ijkm}\f$ (Real part and pure real operator).
    double Ced_LNP = 0.; ///< The dimension-6 operator coefficient \f$(C_{ed})_{ijkm}\f$ (Imaginary part).
    double Cle_LNP = 0.; ///< The dimension-6 operator coefficient \f$(C_{le})_{ijkm}\f$ (Real part and pure real operator).
    double Clu_LNP = 0.; ///< The dimension-6 operator coefficient \f$(C_{lu})_{ijkm}\f$ (Imaginary part).
    double Cld_LNP = 0.; ///< The dimension-6 operator coefficient \f$(C_{ld})_{ijkm}\f$ (Real part and pure real operator).
    double Cqe_LNP = 0.; ///< The dimension-6 operator coefficient \f$(C_{qe})_{ijkm}\f$ (Imaginary part).
    
    double Cqq1_aabb_LNP = 0.; ///< The dimension-6 operator coefficient \f$(C_{qq}^{(1)})_{ijkm}\f$ (Imaginary part).
    double Cqq1_abba_LNP = 0.; ///< The dimension-6 operator coefficient \f$(C_{qq}^{(1)})_{ijkm}\f$ (Imaginary part).
    double Cqq3_aabb_LNP = 0.; ///< The dimension-6 operator coefficient \f$(C_{qq}^{(3)})_{ijkm}\f$ (Real part and pure real operator).
    double Cqq3_abba_LNP = 0.; ///< The dimension-6 operator coefficient \f$(C_{qq}^{(3)})_{ijkm}\f$ (Real part and pure real operator).
    double Cuu_aabb_LNP = 0.; ///< The dimension-6 operator coefficient \f$(C_{uu})_{ijkm}\f$ (Imaginary part).
    double Cuu_abba_LNP = 0.; ///< The dimension-6 operator coefficient \f$(C_{uu})_{ijkm}\f$ (Imaginary part).
    double Cdd_aabb_LNP = 0.; ///< The dimension-6 operator coefficient \f$(C_{dd})_{ijkm}\f$ (Real part and pure real operator).
    double Cdd_abba_LNP = 0.; ///< The dimension-6 operator coefficient \f$(C_{dd})_{ijkm}\f$ (Real part and pure real operator).

    double Cud1_LNP = 0.; ///< The dimension-6 operator coefficient \f$(C_{ud}^{(1)})_{ijkm}\f$ (Imaginary part).
    double Cud8_LNP = 0.; ///< The dimension-6 operator coefficient \f$(C_{ud}^{(8)})_{ijkm}\f$ (Real part and pure real operator).
    double Cqu1_LNP = 0.; ///< The dimension-6 operator coefficient \f$(C_{qu}^{(1)})_{ijkm}\f$ (Imaginary part).
    double Cqu8_LNP = 0.; ///< The dimension-6 operator coefficient \f$(C_{qu}^{(8)})_{ijkm}\f$ (Real part and pure real operator).
    double Cqd1_LNP = 0.; ///< The dimension-6 operator coefficient \f$(C_{qd}^{(1)})_{ijkm}\f$ (Imaginary part).
    double Cqd8_LNP = 0.; ///< The dimension-6 operator coefficient \f$(C_{qd}^{(8)})_{ijkm}\f$ (Real part and pure real operator).
    
    

    
    virtual void setParameter(const std::string name, const double& value);
    
    /**
     * @brief An auxiliary method to set the WC of the general class
     */
    void setNPSMEFTd6GeneralParameters();
    
private:

};

#endif /* NPSMEFTD6MFVLONOCKM_H */

