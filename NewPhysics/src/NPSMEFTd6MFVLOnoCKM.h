/*
 * Click nbfs://nbhost/SystemFileSystem/Templates/Licenses/license-default.txt to change this license
 * Click nbfs://nbhost/SystemFileSystem/Templates/cppFiles/class.h to edit this template
 */

/* 
 * File:   NPSMEFTd6MFVLOnoCKM.h
 * Author: silvest
 *
 * Created on 5 settembre 2023, 15.54
 */

#ifndef NPSMEFTD6MFVLONOCKM_H
#define NPSMEFTD6MFVLONOCKM_H

#include "NPSMEFTd6General.h"

class NPSMEFTd6MFVLOnoCKM: public NPSMEFTd6General {
public:
    
    static const int NNPSMEFTd6MFVLOnoCKMVars = 53;
    
    static std::string NPSMEFTd6MFVLOnoCKMVars[NNPSMEFTd6MFVLOnoCKMVars];
    
    NPSMEFTd6MFVLOnoCKM();
    NPSMEFTd6MFVLOnoCKM(const NPSMEFTd6MFVLOnoCKM& orig);
    virtual ~NPSMEFTd6MFVLOnoCKM();
    
    
protected:
    
    
    
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
    double CHud_LNP = 0.; ///< The dimension-6 operator coefficient \f$(C_{Hud})_{ij}\f$.
    double CeH_LNP = 0.; ///< The dimension-6 operator coefficient \f$(C_{eH})_{ij}\f$.
    double CuH_LNP = 0.; ///< The dimension-6 operator coefficient \f$(C_{uH})_{ij}\f$.
    double CdH_LNP = 0.; ///< The dimension-6 operator coefficient \f$(C_{dH})_{ij}\f$.
    double CuG_LNP = 0.; ///< The dimension-6 operator coefficient \f$(C_{uG})_{ij}\f$.
    double CuW_LNP = 0.; ///< The dimension-6 operator coefficient \f$(C_{uW})_{ij}\f$.
    double CuB_LNP = 0.; ///< The dimension-6 operator coefficient \f$(C_{uB})_{ij}\f$.
    double CdG_LNP = 0.; ///< The dimension-6 operator coefficient \f$(C_{dG})_{ij}\f$.
    double CdW_LNP = 0.; ///< The dimension-6 operator coefficient \f$(C_{dW})_{ij}\f$.
    double CdB_LNP = 0.; ///< The dimension-6 operator coefficient \f$(C_{dB})_{ij}\f$.
    double CeW_LNP = 0.; ///< The dimension-6 operator coefficient \f$(C_{eW})_{ij}\f$.
    double CeB_LNP = 0.; ///< The dimension-6 operator coefficient \f$(C_{eB})_{ij}\f$.

    
    double Cll_LNP = 0.; ///< The dimension-6 operator coefficient \f$(C_{ll})_{ijkm}\f$ (Real part and pure real operator).
    double Clq1_LNP = 0.; ///< The dimension-6 operator coefficient \f$(C_{lq}^{(1)})_{ijkm}\f$ (Imaginary part).
    double Clq3_LNP = 0.; ///< The dimension-6 operator coefficient \f$(C_{lq}^{(3)})_{ijkm}\f$ (Real part and pure real operator).
    double Cee_LNP = 0.; ///< The dimension-6 operator coefficient \f$(C_{ee})_{ijkm}\f$ (Imaginary part).
    double Ceu_LNP = 0.; ///< The dimension-6 operator coefficient \f$(C_{eu})_{ijkm}\f$ (Real part and pure real operator).
    double Ced_LNP = 0.; ///< The dimension-6 operator coefficient \f$(C_{ed})_{ijkm}\f$ (Imaginary part).
    double Cle_LNP = 0.; ///< The dimension-6 operator coefficient \f$(C_{le})_{ijkm}\f$ (Real part and pure real operator).
    double Clu_LNP = 0.; ///< The dimension-6 operator coefficient \f$(C_{lu})_{ijkm}\f$ (Imaginary part).
    double Cld_LNP = 0.; ///< The dimension-6 operator coefficient \f$(C_{ld})_{ijkm}\f$ (Real part and pure real operator).
    double Cqe_LNP = 0.; ///< The dimension-6 operator coefficient \f$(C_{qe})_{ijkm}\f$ (Imaginary part).
    double Cledq_LNP = 0.; ///< The dimension-6 operator coefficient \f$(C_{ledq})_{ijkm}\f$ (Real part and pure real operator).
    double Cqq1_LNP = 0.; ///< The dimension-6 operator coefficient \f$(C_{qq}^{(1)})_{ijkm}\f$ (Imaginary part).
    double Cqq3_LNP = 0.; ///< The dimension-6 operator coefficient \f$(C_{qq}^{(3)})_{ijkm}\f$ (Real part and pure real operator).
    double Cuu_LNP = 0.; ///< The dimension-6 operator coefficient \f$(C_{uu})_{ijkm}\f$ (Imaginary part).
    double Cdd_LNP = 0.; ///< The dimension-6 operator coefficient \f$(C_{dd})_{ijkm}\f$ (Real part and pure real operator).
    double Cud1_LNP = 0.; ///< The dimension-6 operator coefficient \f$(C_{ud}^{(1)})_{ijkm}\f$ (Imaginary part).
    double Cud8_LNP = 0.; ///< The dimension-6 operator coefficient \f$(C_{ud}^{(8)})_{ijkm}\f$ (Real part and pure real operator).
    double Cqu1_LNP = 0.; ///< The dimension-6 operator coefficient \f$(C_{qu}^{(1)})_{ijkm}\f$ (Imaginary part).
    double Cqu8_LNP = 0.; ///< The dimension-6 operator coefficient \f$(C_{qu}^{(8)})_{ijkm}\f$ (Real part and pure real operator).
    double Cqd1_LNP = 0.; ///< The dimension-6 operator coefficient \f$(C_{qd}^{(1)})_{ijkm}\f$ (Imaginary part).
    double Cqd8_LNP = 0.; ///< The dimension-6 operator coefficient \f$(C_{qd}^{(8)})_{ijkm}\f$ (Real part and pure real operator).
    double Cquqd1_LNP = 0.; ///< The dimension-6 operator coefficient \f$(C_{quqd}^{(1)})_{ijkm}\f$ (Imaginary part).
    double Cquqd8_LNP = 0.; ///< The dimension-6 operator coefficient \f$(C_{quqd}^{(8)})_{ijkm}\f$ (Real part and pure real operator).
    double Clequ1_LNP = 0.; ///< The dimension-6 operator coefficient \f$(C_{lequ}^{(1)})_{ijkm}\f$ (Real part and pure real operator).
    double Clequ3_LNP = 0.; ///< The dimension-6 operator coefficient \f$(C_{lequ}^{(3)})_{ijkm}\f$ (Imaginary part).

    
private:

};

#endif /* NPSMEFTD6MFVLONOCKM_H */

