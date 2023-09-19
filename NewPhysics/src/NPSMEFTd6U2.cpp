/* 
 * Copyright (C) 2012 HEPfit Collaboration
 *
 *
 * For the licensing terms see doc/COPYING.
 * 
 * Created on 18 September 2023, 16:23
 * 
 */


#include "NPSMEFTd6U2.h"




std::string NPSMEFTd6U2::NPSMEFTd6U2Vars[NNPSMEFTd6U2Vars] = {
    "CG", "CW", "CHG", "CHW", "CHB", "CHWB", "CHD", "CT", "CHbox", "CH", 
//"CGtilde", "CWtilde", "CHGtilde", "CHWtilde", "CHBtilde", "CHWtildeB", 
"CHl1_aar_LNP", "CHl1_33r_LNP", "CHl3_aar_LNP", "CHl3_33r_LNP", 
"CHe_aar_LNP", "CHe_33r_LNP", "CHq1_aar_LNP", "CHq1_33r_LNP", "CHq3_aar_LNP", "CHq3_33r_LNP", 
"CHu_aar_LNP", "CHu_33r_LNP", "CHd_aar_LNP", "CHd_33r_LNP", "CHud_aar_LNP", "CHud_33r_LNP", 
"CeH_aar_LNP", "CeH_33r_LNP", "CuH_aar_LNP", "CuH_33r_LNP", "CdH_aar_LNP", "CdH_33r_LNP", 
"CuG_aar_LNP", "CuG_33r_LNP", "CuW_aar_LNP", "CuW_33r_LNP", "CuB_aar_LNP", "CuB_33r_LNP", 
"CdG_aar_LNP", "CdG_33r_LNP", "CdW_aar_LNP", "CdW_33r_LNP", "CdB_aar_LNP", "CdB_33r_LNP", 
"CeW_aar_LNP", "CeW_33r_LNP", "CeB_aar_LNP", "CeB_33r_LNP", 
"Cll_aabbr_LNP", "Cll_aa33r_LNP", "Cll_33aar_LNP", "Cll_3aa3r_LNP", "Cll_3333r_LNP", 
"Clq1_aabbr_LNP", "Clq1_aa33r_LNP", "Clq1_33aar_LNP", "Clq1_3333r_LNP", 
"Clq3_aabbr_LNP", "Clq3_aa33r_LNP", "Clq3_33aar_LNP", "Clq3_3333r_LNP", 
"Cee_aabbr_LNP", "Cee_aa33r_LNP", "Cee_33aar_LNP", "Cee_3aa3r_LNP", "Cee_3333r_LNP", 
"Ceu_aabbr_LNP", "Ceu_aa33r_LNP", "Ceu_33aar_LNP", "Ceu_3333r_LNP", 
"Ced_aabbr_LNP", "Ced_aa33r_LNP", "Ced_33aar_LNP", "Ced_3333r_LNP", 
"Cle_aabbr_LNP", "Cle_aa33r_LNP", "Cle_33aar_LNP", "Cle_3333r_LNP", 
"Clu_aabbr_LNP", "Clu_aa33r_LNP", "Clu_33aar_LNP", "Clu_3333r_LNP", 
"Cld_aabbr_LNP", "Cld_aa33r_LNP", "Cld_33aar_LNP", "Cld_3333r_LNP", 
"Cqe_aabbr_LNP", "Cqe_aa33r_LNP", "Cqe_33aar_LNP", "Cqe_3333r_LNP", 
"Cledq_aabbr_LNP", "Cledq_aa33r_LNP", "Cledq_33aar_LNP", "Cledq_3333r_LNP", 
"Cqq1_aabbr_LNP", "Cqq1_aa33r_LNP", "Cqq1_33aar_LNP", "Cqq1_3aa3r_LNP", "Cqq1_3333r_LNP", 
"Cqq3_aabbr_LNP", "Cqq3_aa33r_LNP", "Cqq3_33aar_LNP", "Cqq3_3aa3r_LNP", "Cqq3_3333r_LNP", 
"Cuu_aabbr_LNP", "Cuu_aa33r_LNP", "Cuu_33aar_LNP", "Cuu_3aa3r_LNP", "Cuu_3333r_LNP", 
"Cdd_aabbr_LNP", "Cdd_aa33r_LNP", "Cdd_33aar_LNP", "Cdd_3aa3r_LNP", "Cdd_3333r_LNP", 
"Cud1_aabbr_LNP", "Cud1_aa33r_LNP", "Cud1_33aar_LNP", "Cud1_3333r_LNP", 
"Cud8_aabbr_LNP", "Cud8_aa33r_LNP", "Cud8_33aar_LNP", "Cud8_3333r_LNP", 
"Cqu1_aabbr_LNP", "Cqu1_aa33r_LNP", "Cqu1_33aar_LNP", "Cqu1_3333r_LNP", 
"Cqu8_aabbr_LNP", "Cqu8_aa33r_LNP", "Cqu8_33aar_LNP", "Cqu8_3333r_LNP", 
"Cqd1_aabbr_LNP", "Cqd1_aa33r_LNP", "Cqd1_33aar_LNP", "Cqd1_3333r_LNP", 
"Cqd8_aabbr_LNP", "Cqd8_aa33r_LNP", "Cqd8_33aar_LNP", "Cqd8_3333r_LNP", 
"Cquqd1_aabbr_LNP", "Cquqd1_aa33r_LNP", "Cquqd1_33aar_LNP", "Cquqd1_3333r_LNP", 
"Cquqd8_aabbr_LNP", "Cquqd8_aa33r_LNP", "Cquqd8_33aar_LNP", "Cquqd8_3333r_LNP", 
"Clequ1_aabbr_LNP", "Clequ1_aa33r_LNP", "Clequ1_33aar_LNP", "Clequ1_3333r_LNP", 
"Clequ3_aabbr_LNP", "Clequ3_aa33r_LNP", "Clequ3_33aar_LNP", "Clequ3_3333r_LNP", 
"Lambda_NP"
};


NPSMEFTd6U2::NPSMEFTd6U2()
: NPSMEFTd6General()
{
    
    ModelParamMap.insert(std::make_pair("CG_LNP", std::cref(CG_LNP))); 
ModelParamMap.insert(std::make_pair("CW_LNP", std::cref(CW_LNP))); 
ModelParamMap.insert(std::make_pair("C2B_LNP", std::cref(C2B_LNP))); 
ModelParamMap.insert(std::make_pair("C2W_LNP", std::cref(C2W_LNP))); 
ModelParamMap.insert(std::make_pair("C2BS_LNP", std::cref(C2BS_LNP))); 
ModelParamMap.insert(std::make_pair("C2WS_LNP", std::cref(C2WS_LNP))); 
ModelParamMap.insert(std::make_pair("CHG_LNP", std::cref(CHG_LNP))); 
ModelParamMap.insert(std::make_pair("CHW_LNP", std::cref(CHW_LNP))); 
ModelParamMap.insert(std::make_pair("CHB_LNP", std::cref(CHB_LNP))); 
ModelParamMap.insert(std::make_pair("CDHB_LNP", std::cref(CDHB_LNP))); 
ModelParamMap.insert(std::make_pair("CDHW_LNP", std::cref(CDHW_LNP))); 
ModelParamMap.insert(std::make_pair("CDB_LNP", std::cref(CDB_LNP))); 
ModelParamMap.insert(std::make_pair("CDW_LNP", std::cref(CDW_LNP))); 
ModelParamMap.insert(std::make_pair("CHWB_LNP", std::cref(CHWB_LNP))); 
ModelParamMap.insert(std::make_pair("CHD_LNP", std::cref(CHD_LNP))); 
ModelParamMap.insert(std::make_pair("CT_LNP", std::cref(CT_LNP))); 
ModelParamMap.insert(std::make_pair("CHbox_LNP", std::cref(CHbox_LNP))); 
ModelParamMap.insert(std::make_pair("CH_LNP", std::cref(CH_LNP))); 
ModelParamMap.insert(std::make_pair("CGtilde_LNP", std::cref(CGtilde_LNP))); 
ModelParamMap.insert(std::make_pair("CWtilde_LNP", std::cref(CWtilde_LNP))); 
ModelParamMap.insert(std::make_pair("CHGtilde_LNP", std::cref(CHGtilde_LNP))); 
ModelParamMap.insert(std::make_pair("CHWtilde_LNP", std::cref(CHWtilde_LNP))); 
ModelParamMap.insert(std::make_pair("CHBtilde_LNP", std::cref(CHBtilde_LNP))); 
ModelParamMap.insert(std::make_pair("CHWtildeB_LNP", std::cref(CHWtildeB_LNP))); 


ModelParamMap.insert(std::make_pair("CHl1_aar_LNP", std::cref(CHl1_aar_LNP))); 
ModelParamMap.insert(std::make_pair("CHl1_33r_LNP", std::cref(CHl1_33r_LNP))); 
ModelParamMap.insert(std::make_pair("CHl3_aar_LNP", std::cref(CHl3_aar_LNP))); 
ModelParamMap.insert(std::make_pair("CHl3_33r_LNP", std::cref(CHl3_33r_LNP))); 
ModelParamMap.insert(std::make_pair("CHe_aar_LNP", std::cref(CHe_aar_LNP))); 
ModelParamMap.insert(std::make_pair("CHe_33r_LNP", std::cref(CHe_33r_LNP))); 
ModelParamMap.insert(std::make_pair("CeH_aar_LNP", std::cref(CeH_aar_LNP))); 
ModelParamMap.insert(std::make_pair("CeH_33r_LNP", std::cref(CeH_33r_LNP))); 
ModelParamMap.insert(std::make_pair("Cll_aabbr_LNP", std::cref(Cll_aabbr_LNP))); 
ModelParamMap.insert(std::make_pair("Cll_aa33r_LNP", std::cref(Cll_aa33r_LNP))); 
ModelParamMap.insert(std::make_pair("Cll_33aar_LNP", std::cref(Cll_33aar_LNP))); 
ModelParamMap.insert(std::make_pair("Cll_3aa3r_LNP", std::cref(Cll_3aa3r_LNP))); 
ModelParamMap.insert(std::make_pair("Cll_3333r_LNP", std::cref(Cll_3333r_LNP))); 
ModelParamMap.insert(std::make_pair("Cee_aabbr_LNP", std::cref(Cee_aabbr_LNP))); 
ModelParamMap.insert(std::make_pair("Cee_aa33r_LNP", std::cref(Cee_aa33r_LNP))); 
ModelParamMap.insert(std::make_pair("Cee_33aar_LNP", std::cref(Cee_33aar_LNP))); 
ModelParamMap.insert(std::make_pair("Cee_3aa3r_LNP", std::cref(Cee_3aa3r_LNP))); 
ModelParamMap.insert(std::make_pair("Cee_3333r_LNP", std::cref(Cee_3333r_LNP))); 
ModelParamMap.insert(std::make_pair("Cle_aabbr_LNP", std::cref(Cle_aabbr_LNP))); 
ModelParamMap.insert(std::make_pair("Cle_aa33r_LNP", std::cref(Cle_aa33r_LNP))); 
ModelParamMap.insert(std::make_pair("Cle_33aar_LNP", std::cref(Cle_33aar_LNP))); 
ModelParamMap.insert(std::make_pair("Cle_3333r_LNP", std::cref(Cle_3333r_LNP))); 



ModelParamMap.insert(std::make_pair("CHq1_aar_LNP", std::cref(CHq1_aar_LNP))); 
ModelParamMap.insert(std::make_pair("CHq1_33r_LNP", std::cref(CHq1_33r_LNP))); 
ModelParamMap.insert(std::make_pair("CHq3_aar_LNP", std::cref(CHq3_aar_LNP))); 
ModelParamMap.insert(std::make_pair("CHq3_33r_LNP", std::cref(CHq3_33r_LNP))); 
ModelParamMap.insert(std::make_pair("CHu_aar_LNP", std::cref(CHu_aar_LNP))); 
ModelParamMap.insert(std::make_pair("CHu_33r_LNP", std::cref(CHu_33r_LNP))); 
ModelParamMap.insert(std::make_pair("CHd_aar_LNP", std::cref(CHd_aar_LNP))); 
ModelParamMap.insert(std::make_pair("CHd_33r_LNP", std::cref(CHd_33r_LNP))); 
ModelParamMap.insert(std::make_pair("CHud_aar_LNP", std::cref(CHud_aar_LNP))); 
ModelParamMap.insert(std::make_pair("CHud_33r_LNP", std::cref(CHud_33r_LNP))); 
ModelParamMap.insert(std::make_pair("CuH_aar_LNP", std::cref(CuH_aar_LNP))); 
ModelParamMap.insert(std::make_pair("CuH_33r_LNP", std::cref(CuH_33r_LNP))); 
ModelParamMap.insert(std::make_pair("CdH_aar_LNP", std::cref(CdH_aar_LNP))); 
ModelParamMap.insert(std::make_pair("CdH_33r_LNP", std::cref(CdH_33r_LNP))); 
ModelParamMap.insert(std::make_pair("CuG_aar_LNP", std::cref(CuG_aar_LNP))); 
ModelParamMap.insert(std::make_pair("CuG_33r_LNP", std::cref(CuG_33r_LNP))); 
ModelParamMap.insert(std::make_pair("CuW_aar_LNP", std::cref(CuW_aar_LNP))); 
ModelParamMap.insert(std::make_pair("CuW_33r_LNP", std::cref(CuW_33r_LNP))); 
ModelParamMap.insert(std::make_pair("CuB_aar_LNP", std::cref(CuB_aar_LNP))); 
ModelParamMap.insert(std::make_pair("CuB_33r_LNP", std::cref(CuB_33r_LNP))); 
ModelParamMap.insert(std::make_pair("CdG_aar_LNP", std::cref(CdG_aar_LNP))); 
ModelParamMap.insert(std::make_pair("CdG_33r_LNP", std::cref(CdG_33r_LNP))); 
ModelParamMap.insert(std::make_pair("CdW_aar_LNP", std::cref(CdW_aar_LNP))); 
ModelParamMap.insert(std::make_pair("CdW_33r_LNP", std::cref(CdW_33r_LNP))); 
ModelParamMap.insert(std::make_pair("CdB_aar_LNP", std::cref(CdB_aar_LNP))); 
ModelParamMap.insert(std::make_pair("CdB_33r_LNP", std::cref(CdB_33r_LNP))); 
ModelParamMap.insert(std::make_pair("CeW_aar_LNP", std::cref(CeW_aar_LNP))); 
ModelParamMap.insert(std::make_pair("CeW_33r_LNP", std::cref(CeW_33r_LNP))); 
ModelParamMap.insert(std::make_pair("CeB_aar_LNP", std::cref(CeB_aar_LNP))); 
ModelParamMap.insert(std::make_pair("CeB_33r_LNP", std::cref(CeB_33r_LNP))); 
ModelParamMap.insert(std::make_pair("Cqq1_aabbr_LNP", std::cref(Cqq1_aabbr_LNP))); 
ModelParamMap.insert(std::make_pair("Cqq1_aa33r_LNP", std::cref(Cqq1_aa33r_LNP))); 
ModelParamMap.insert(std::make_pair("Cqq1_33aar_LNP", std::cref(Cqq1_33aar_LNP))); 
ModelParamMap.insert(std::make_pair("Cqq1_3aa3r_LNP", std::cref(Cqq1_3aa3r_LNP))); 
ModelParamMap.insert(std::make_pair("Cqq1_3333r_LNP", std::cref(Cqq1_3333r_LNP))); 
ModelParamMap.insert(std::make_pair("Cqq3_aabbr_LNP", std::cref(Cqq3_aabbr_LNP))); 
ModelParamMap.insert(std::make_pair("Cqq3_aa33r_LNP", std::cref(Cqq3_aa33r_LNP))); 
ModelParamMap.insert(std::make_pair("Cqq3_33aar_LNP", std::cref(Cqq3_33aar_LNP))); 
ModelParamMap.insert(std::make_pair("Cqq3_3aa3r_LNP", std::cref(Cqq3_3aa3r_LNP))); 
ModelParamMap.insert(std::make_pair("Cqq3_3333r_LNP", std::cref(Cqq3_3333r_LNP))); 
ModelParamMap.insert(std::make_pair("Cuu_aabbr_LNP", std::cref(Cuu_aabbr_LNP))); 
ModelParamMap.insert(std::make_pair("Cuu_aa33r_LNP", std::cref(Cuu_aa33r_LNP))); 
ModelParamMap.insert(std::make_pair("Cuu_33aar_LNP", std::cref(Cuu_33aar_LNP))); 
ModelParamMap.insert(std::make_pair("Cuu_3aa3r_LNP", std::cref(Cuu_3aa3r_LNP))); 
ModelParamMap.insert(std::make_pair("Cuu_3333r_LNP", std::cref(Cuu_3333r_LNP))); 
ModelParamMap.insert(std::make_pair("Cdd_aabbr_LNP", std::cref(Cdd_aabbr_LNP))); 
ModelParamMap.insert(std::make_pair("Cdd_aa33r_LNP", std::cref(Cdd_aa33r_LNP))); 
ModelParamMap.insert(std::make_pair("Cdd_33aar_LNP", std::cref(Cdd_33aar_LNP))); 
ModelParamMap.insert(std::make_pair("Cdd_3aa3r_LNP", std::cref(Cdd_3aa3r_LNP))); 
ModelParamMap.insert(std::make_pair("Cdd_3333r_LNP", std::cref(Cdd_3333r_LNP))); 
ModelParamMap.insert(std::make_pair("Cud1_aabbr_LNP", std::cref(Cud1_aabbr_LNP))); 
ModelParamMap.insert(std::make_pair("Cud1_aa33r_LNP", std::cref(Cud1_aa33r_LNP))); 
ModelParamMap.insert(std::make_pair("Cud1_33aar_LNP", std::cref(Cud1_33aar_LNP))); 
ModelParamMap.insert(std::make_pair("Cud1_3333r_LNP", std::cref(Cud1_3333r_LNP))); 
ModelParamMap.insert(std::make_pair("Cud8_aabbr_LNP", std::cref(Cud8_aabbr_LNP))); 
ModelParamMap.insert(std::make_pair("Cud8_aa33r_LNP", std::cref(Cud8_aa33r_LNP))); 
ModelParamMap.insert(std::make_pair("Cud8_33aar_LNP", std::cref(Cud8_33aar_LNP))); 
ModelParamMap.insert(std::make_pair("Cud8_3333r_LNP", std::cref(Cud8_3333r_LNP))); 
ModelParamMap.insert(std::make_pair("Cqu1_aabbr_LNP", std::cref(Cqu1_aabbr_LNP))); 
ModelParamMap.insert(std::make_pair("Cqu1_aa33r_LNP", std::cref(Cqu1_aa33r_LNP))); 
ModelParamMap.insert(std::make_pair("Cqu1_33aar_LNP", std::cref(Cqu1_33aar_LNP))); 
ModelParamMap.insert(std::make_pair("Cqu1_3333r_LNP", std::cref(Cqu1_3333r_LNP))); 
ModelParamMap.insert(std::make_pair("Cqu8_aabbr_LNP", std::cref(Cqu8_aabbr_LNP))); 
ModelParamMap.insert(std::make_pair("Cqu8_aa33r_LNP", std::cref(Cqu8_aa33r_LNP))); 
ModelParamMap.insert(std::make_pair("Cqu8_33aar_LNP", std::cref(Cqu8_33aar_LNP))); 
ModelParamMap.insert(std::make_pair("Cqu8_3333r_LNP", std::cref(Cqu8_3333r_LNP))); 
ModelParamMap.insert(std::make_pair("Cqd1_aabbr_LNP", std::cref(Cqd1_aabbr_LNP))); 
ModelParamMap.insert(std::make_pair("Cqd1_aa33r_LNP", std::cref(Cqd1_aa33r_LNP))); 
ModelParamMap.insert(std::make_pair("Cqd1_33aar_LNP", std::cref(Cqd1_33aar_LNP))); 
ModelParamMap.insert(std::make_pair("Cqd1_3333r_LNP", std::cref(Cqd1_3333r_LNP))); 
ModelParamMap.insert(std::make_pair("Cqd8_aabbr_LNP", std::cref(Cqd8_aabbr_LNP))); 
ModelParamMap.insert(std::make_pair("Cqd8_aa33r_LNP", std::cref(Cqd8_aa33r_LNP))); 
ModelParamMap.insert(std::make_pair("Cqd8_33aar_LNP", std::cref(Cqd8_33aar_LNP))); 
ModelParamMap.insert(std::make_pair("Cqd8_3333r_LNP", std::cref(Cqd8_3333r_LNP))); 
ModelParamMap.insert(std::make_pair("Cquqd1_aabbr_LNP", std::cref(Cquqd1_aabbr_LNP))); 
ModelParamMap.insert(std::make_pair("Cquqd1_aa33r_LNP", std::cref(Cquqd1_aa33r_LNP))); 
ModelParamMap.insert(std::make_pair("Cquqd1_33aar_LNP", std::cref(Cquqd1_33aar_LNP))); 
ModelParamMap.insert(std::make_pair("Cquqd1_3333r_LNP", std::cref(Cquqd1_3333r_LNP))); 
ModelParamMap.insert(std::make_pair("Cquqd8_aabbr_LNP", std::cref(Cquqd8_aabbr_LNP))); 
ModelParamMap.insert(std::make_pair("Cquqd8_aa33r_LNP", std::cref(Cquqd8_aa33r_LNP))); 
ModelParamMap.insert(std::make_pair("Cquqd8_33aar_LNP", std::cref(Cquqd8_33aar_LNP))); 
ModelParamMap.insert(std::make_pair("Cquqd8_3333r_LNP", std::cref(Cquqd8_3333r_LNP))); 


ModelParamMap.insert(std::make_pair("Clq1_aabbr_LNP", std::cref(Clq1_aabbr_LNP))); 
ModelParamMap.insert(std::make_pair("Clq1_aa33r_LNP", std::cref(Clq1_aa33r_LNP))); 
ModelParamMap.insert(std::make_pair("Clq1_33aar_LNP", std::cref(Clq1_33aar_LNP))); 
ModelParamMap.insert(std::make_pair("Clq1_3333r_LNP", std::cref(Clq1_3333r_LNP))); 
ModelParamMap.insert(std::make_pair("Clq3_aabbr_LNP", std::cref(Clq3_aabbr_LNP))); 
ModelParamMap.insert(std::make_pair("Clq3_aa33r_LNP", std::cref(Clq3_aa33r_LNP))); 
ModelParamMap.insert(std::make_pair("Clq3_33aar_LNP", std::cref(Clq3_33aar_LNP))); 
ModelParamMap.insert(std::make_pair("Clq3_3333r_LNP", std::cref(Clq3_3333r_LNP))); 
ModelParamMap.insert(std::make_pair("Ceu_aabbr_LNP", std::cref(Ceu_aabbr_LNP))); 
ModelParamMap.insert(std::make_pair("Ceu_aa33r_LNP", std::cref(Ceu_aa33r_LNP))); 
ModelParamMap.insert(std::make_pair("Ceu_33aar_LNP", std::cref(Ceu_33aar_LNP))); 
ModelParamMap.insert(std::make_pair("Ceu_3333r_LNP", std::cref(Ceu_3333r_LNP))); 
ModelParamMap.insert(std::make_pair("Ced_aabbr_LNP", std::cref(Ced_aabbr_LNP))); 
ModelParamMap.insert(std::make_pair("Ced_aa33r_LNP", std::cref(Ced_aa33r_LNP))); 
ModelParamMap.insert(std::make_pair("Ced_33aar_LNP", std::cref(Ced_33aar_LNP))); 
ModelParamMap.insert(std::make_pair("Ced_3333r_LNP", std::cref(Ced_3333r_LNP))); 
ModelParamMap.insert(std::make_pair("Clu_aabbr_LNP", std::cref(Clu_aabbr_LNP))); 
ModelParamMap.insert(std::make_pair("Clu_aa33r_LNP", std::cref(Clu_aa33r_LNP))); 
ModelParamMap.insert(std::make_pair("Clu_33aar_LNP", std::cref(Clu_33aar_LNP))); 
ModelParamMap.insert(std::make_pair("Clu_3333r_LNP", std::cref(Clu_3333r_LNP))); 
ModelParamMap.insert(std::make_pair("Cld_aabbr_LNP", std::cref(Cld_aabbr_LNP))); 
ModelParamMap.insert(std::make_pair("Cld_aa33r_LNP", std::cref(Cld_aa33r_LNP))); 
ModelParamMap.insert(std::make_pair("Cld_33aar_LNP", std::cref(Cld_33aar_LNP))); 
ModelParamMap.insert(std::make_pair("Cld_3333r_LNP", std::cref(Cld_3333r_LNP))); 
ModelParamMap.insert(std::make_pair("Cqe_aabbr_LNP", std::cref(Cqe_aabbr_LNP))); 
ModelParamMap.insert(std::make_pair("Cqe_aa33r_LNP", std::cref(Cqe_aa33r_LNP))); 
ModelParamMap.insert(std::make_pair("Cqe_33aar_LNP", std::cref(Cqe_33aar_LNP))); 
ModelParamMap.insert(std::make_pair("Cqe_3333r_LNP", std::cref(Cqe_3333r_LNP))); 
ModelParamMap.insert(std::make_pair("Cledq_aabbr_LNP", std::cref(Cledq_aabbr_LNP))); 
ModelParamMap.insert(std::make_pair("Cledq_aa33r_LNP", std::cref(Cledq_aa33r_LNP))); 
ModelParamMap.insert(std::make_pair("Cledq_33aar_LNP", std::cref(Cledq_33aar_LNP))); 
ModelParamMap.insert(std::make_pair("Cledq_3333r_LNP", std::cref(Cledq_3333r_LNP))); 
ModelParamMap.insert(std::make_pair("Clequ1_aabbr_LNP", std::cref(Clequ1_aabbr_LNP))); 
ModelParamMap.insert(std::make_pair("Clequ1_aa33r_LNP", std::cref(Clequ1_aa33r_LNP))); 
ModelParamMap.insert(std::make_pair("Clequ1_33aar_LNP", std::cref(Clequ1_33aar_LNP))); 
ModelParamMap.insert(std::make_pair("Clequ1_3333r_LNP", std::cref(Clequ1_3333r_LNP))); 
ModelParamMap.insert(std::make_pair("Clequ3_aabbr_LNP", std::cref(Clequ3_aabbr_LNP))); 
ModelParamMap.insert(std::make_pair("Clequ3_aa33r_LNP", std::cref(Clequ3_aa33r_LNP))); 
ModelParamMap.insert(std::make_pair("Clequ3_33aar_LNP", std::cref(Clequ3_33aar_LNP))); 
ModelParamMap.insert(std::make_pair("Clequ3_3333r_LNP", std::cref(Clequ3_3333r_LNP))); 

    

}



void NPSMEFTd6U2::setParameter(const std::string name, const double& value)
{
    
 if (name.compare("CG_LNP") == 0) { 
	CG_LNP =  value;
	SMEFTEvol.SetCoefficient("CG", value * LambdaNPm2);
} else if (name.compare("CW_LNP") == 0) { 
	CW_LNP =  value;
	SMEFTEvol.SetCoefficient("CW", value * LambdaNPm2);
} else if (name.compare("C2B_LNP") == 0) { 
	C2B_LNP =  value;
	SMEFTEvol.SetCoefficient("C2B", value * LambdaNPm2);
} else if (name.compare("C2W_LNP") == 0) { 
	C2W_LNP =  value;
	SMEFTEvol.SetCoefficient("C2W", value * LambdaNPm2);
} else if (name.compare("C2BS_LNP") == 0) { 
	C2BS_LNP =  value;
	SMEFTEvol.SetCoefficient("C2BS", value * LambdaNPm2);
} else if (name.compare("C2WS_LNP") == 0) { 
	C2WS_LNP =  value;
	SMEFTEvol.SetCoefficient("C2WS", value * LambdaNPm2);
} else if (name.compare("CHG_LNP") == 0) { 
	CHG_LNP =  value;
	SMEFTEvol.SetCoefficient("CHG", value * LambdaNPm2);
} else if (name.compare("CHW_LNP") == 0) { 
	CHW_LNP =  value;
	SMEFTEvol.SetCoefficient("CHW", value * LambdaNPm2);
} else if (name.compare("CHB_LNP") == 0) { 
	CHB_LNP =  value;
	SMEFTEvol.SetCoefficient("CHB", value * LambdaNPm2);
} else if (name.compare("CDHB_LNP") == 0) { 
	CDHB_LNP =  value;
	SMEFTEvol.SetCoefficient("CDHB", value * LambdaNPm2);
} else if (name.compare("CDHW_LNP") == 0) { 
	CDHW_LNP =  value;
	SMEFTEvol.SetCoefficient("CDHW", value * LambdaNPm2);
} else if (name.compare("CDB_LNP") == 0) { 
	CDB_LNP =  value;
	SMEFTEvol.SetCoefficient("CDB", value * LambdaNPm2);
} else if (name.compare("CDW_LNP") == 0) { 
	CDW_LNP =  value;
	SMEFTEvol.SetCoefficient("CDW", value * LambdaNPm2);
} else if (name.compare("CHWB_LNP") == 0) { 
	CHWB_LNP =  value;
	SMEFTEvol.SetCoefficient("CHWB", value * LambdaNPm2);
} else if (name.compare("CHD_LNP") == 0) { 
	CHD_LNP =  value;
	SMEFTEvol.SetCoefficient("CHD", value * LambdaNPm2);
} else if (name.compare("CT_LNP") == 0) { 
	CT_LNP =  value;
	SMEFTEvol.SetCoefficient("CT", value * LambdaNPm2);
} else if (name.compare("CHbox_LNP") == 0) { 
	CHbox_LNP =  value;
	SMEFTEvol.SetCoefficient("CHbox", value * LambdaNPm2);
} else if (name.compare("CH_LNP") == 0) { 
	CH_LNP =  value;
	SMEFTEvol.SetCoefficient("CH", value * LambdaNPm2);
} else if (name.compare("CGtilde_LNP") == 0) { 
	CGtilde_LNP =  value;
	SMEFTEvol.SetCoefficient("CGtilde", value * LambdaNPm2);
} else if (name.compare("CWtilde_LNP") == 0) { 
	CWtilde_LNP =  value;
	SMEFTEvol.SetCoefficient("CWtilde", value * LambdaNPm2);
} else if (name.compare("CHGtilde_LNP") == 0) { 
	CHGtilde_LNP =  value;
	SMEFTEvol.SetCoefficient("CHGtilde", value * LambdaNPm2);
} else if (name.compare("CHWtilde_LNP") == 0) { 
	CHWtilde_LNP =  value;
	SMEFTEvol.SetCoefficient("CHWtilde", value * LambdaNPm2);
} else if (name.compare("CHBtilde_LNP") == 0) { 
	CHBtilde_LNP =  value;
	SMEFTEvol.SetCoefficient("CHBtilde", value * LambdaNPm2);
} else if (name.compare("CHWtildeB_LNP") == 0) { 
	CHWtildeB_LNP =  value;
	SMEFTEvol.SetCoefficient("CHWtildeB", value * LambdaNPm2);
} else if (name.compare("CHl1_aar_LNP") == 0) {
	CHl1_aar_LNP =  value;
	SMEFTEvol.SetCoefficient("CHl1R", value * LambdaNPm2, 0, 0);
	SMEFTEvol.SetCoefficient("CHl1R", value * LambdaNPm2, 1, 1);
} else if (name.compare("CHl1_33r_LNP") == 0) {
	CHl1_33r_LNP =  value;
	SMEFTEvol.SetCoefficient("CHl1R", value * LambdaNPm2, 2, 2);
} else if (name.compare("CHl3_aar_LNP") == 0) {
	CHl3_aar_LNP =  value;
	SMEFTEvol.SetCoefficient("CHl3R", value * LambdaNPm2, 0, 0);
	SMEFTEvol.SetCoefficient("CHl3R", value * LambdaNPm2, 1, 1);
} else if (name.compare("CHl3_33r_LNP") == 0) {
	CHl3_33r_LNP =  value;
	SMEFTEvol.SetCoefficient("CHl3R", value * LambdaNPm2, 2, 2);
} else if (name.compare("CHe_aar_LNP") == 0) {
	CHe_aar_LNP =  value;
	SMEFTEvol.SetCoefficient("CHeR", value * LambdaNPm2, 0, 0);
	SMEFTEvol.SetCoefficient("CHeR", value * LambdaNPm2, 1, 1);
} else if (name.compare("CHe_33r_LNP") == 0) {
	CHe_33r_LNP =  value;
	SMEFTEvol.SetCoefficient("CHeR", value * LambdaNPm2, 2, 2);
} else if (name.compare("CHq1_aar_LNP") == 0) {
	CHq1_aar_LNP =  value;
	SMEFTEvol.SetCoefficient("CHq1R", value * LambdaNPm2, 0, 0);
	SMEFTEvol.SetCoefficient("CHq1R", value * LambdaNPm2, 1, 1);
} else if (name.compare("CHq1_33r_LNP") == 0) {
	CHq1_33r_LNP =  value;
	SMEFTEvol.SetCoefficient("CHq1R", value * LambdaNPm2, 2, 2);
} else if (name.compare("CHq3_aar_LNP") == 0) {
	CHq3_aar_LNP =  value;
	SMEFTEvol.SetCoefficient("CHq3R", value * LambdaNPm2, 0, 0);
	SMEFTEvol.SetCoefficient("CHq3R", value * LambdaNPm2, 1, 1);
} else if (name.compare("CHq3_33r_LNP") == 0) {
	CHq3_33r_LNP =  value;
	SMEFTEvol.SetCoefficient("CHq3R", value * LambdaNPm2, 2, 2);
} else if (name.compare("CHu_aar_LNP") == 0) {
	CHu_aar_LNP =  value;
	SMEFTEvol.SetCoefficient("CHuR", value * LambdaNPm2, 0, 0);
	SMEFTEvol.SetCoefficient("CHuR", value * LambdaNPm2, 1, 1);
} else if (name.compare("CHu_33r_LNP") == 0) {
	CHu_33r_LNP =  value;
	SMEFTEvol.SetCoefficient("CHuR", value * LambdaNPm2, 2, 2);
} else if (name.compare("CHd_aar_LNP") == 0) {
	CHd_aar_LNP =  value;
	SMEFTEvol.SetCoefficient("CHdR", value * LambdaNPm2, 0, 0);
	SMEFTEvol.SetCoefficient("CHdR", value * LambdaNPm2, 1, 1);
} else if (name.compare("CHd_33r_LNP") == 0) {
	CHd_33r_LNP =  value;
	SMEFTEvol.SetCoefficient("CHdR", value * LambdaNPm2, 2, 2);
} else if (name.compare("CHud_aar_LNP") == 0) {
	CHud_aar_LNP =  value;
	SMEFTEvol.SetCoefficient("CHudR", value * LambdaNPm2, 0, 0);
	SMEFTEvol.SetCoefficient("CHudR", value * LambdaNPm2, 1, 1);
} else if (name.compare("CHud_33r_LNP") == 0) {
	CHud_33r_LNP =  value;
	SMEFTEvol.SetCoefficient("CHudR", value * LambdaNPm2, 2, 2);
} else if (name.compare("CeH_aar_LNP") == 0) {
	CeH_aar_LNP =  value;
	SMEFTEvol.SetCoefficient("CeHR", value * LambdaNPm2, 0, 0);
	SMEFTEvol.SetCoefficient("CeHR", value * LambdaNPm2, 1, 1);
} else if (name.compare("CeH_33r_LNP") == 0) {
	CeH_33r_LNP =  value;
	SMEFTEvol.SetCoefficient("CeHR", value * LambdaNPm2, 2, 2);
} else if (name.compare("CuH_aar_LNP") == 0) {
	CuH_aar_LNP =  value;
	SMEFTEvol.SetCoefficient("CuHR", value * LambdaNPm2, 0, 0);
	SMEFTEvol.SetCoefficient("CuHR", value * LambdaNPm2, 1, 1);
} else if (name.compare("CuH_33r_LNP") == 0) {
	CuH_33r_LNP =  value;
	SMEFTEvol.SetCoefficient("CuHR", value * LambdaNPm2, 2, 2);
} else if (name.compare("CdH_aar_LNP") == 0) {
	CdH_aar_LNP =  value;
	SMEFTEvol.SetCoefficient("CdHR", value * LambdaNPm2, 0, 0);
	SMEFTEvol.SetCoefficient("CdHR", value * LambdaNPm2, 1, 1);
} else if (name.compare("CdH_33r_LNP") == 0) {
	CdH_33r_LNP =  value;
	SMEFTEvol.SetCoefficient("CdHR", value * LambdaNPm2, 2, 2);
} else if (name.compare("CuG_aar_LNP") == 0) {
	CuG_aar_LNP =  value;
	SMEFTEvol.SetCoefficient("CuGR", value * LambdaNPm2, 0, 0);
	SMEFTEvol.SetCoefficient("CuGR", value * LambdaNPm2, 1, 1);
} else if (name.compare("CuG_33r_LNP") == 0) {
	CuG_33r_LNP =  value;
	SMEFTEvol.SetCoefficient("CuGR", value * LambdaNPm2, 2, 2);
} else if (name.compare("CuW_aar_LNP") == 0) {
	CuW_aar_LNP =  value;
	SMEFTEvol.SetCoefficient("CuWR", value * LambdaNPm2, 0, 0);
	SMEFTEvol.SetCoefficient("CuWR", value * LambdaNPm2, 1, 1);
} else if (name.compare("CuW_33r_LNP") == 0) {
	CuW_33r_LNP =  value;
	SMEFTEvol.SetCoefficient("CuWR", value * LambdaNPm2, 2, 2);
} else if (name.compare("CuB_aar_LNP") == 0) {
	CuB_aar_LNP =  value;
	SMEFTEvol.SetCoefficient("CuBR", value * LambdaNPm2, 0, 0);
	SMEFTEvol.SetCoefficient("CuBR", value * LambdaNPm2, 1, 1);
} else if (name.compare("CuB_33r_LNP") == 0) {
	CuB_33r_LNP =  value;
	SMEFTEvol.SetCoefficient("CuBR", value * LambdaNPm2, 2, 2);
} else if (name.compare("CdG_aar_LNP") == 0) {
	CdG_aar_LNP =  value;
	SMEFTEvol.SetCoefficient("CdGR", value * LambdaNPm2, 0, 0);
	SMEFTEvol.SetCoefficient("CdGR", value * LambdaNPm2, 1, 1);
} else if (name.compare("CdG_33r_LNP") == 0) {
	CdG_33r_LNP =  value;
	SMEFTEvol.SetCoefficient("CdGR", value * LambdaNPm2, 2, 2);
} else if (name.compare("CdW_aar_LNP") == 0) {
	CdW_aar_LNP =  value;
	SMEFTEvol.SetCoefficient("CdWR", value * LambdaNPm2, 0, 0);
	SMEFTEvol.SetCoefficient("CdWR", value * LambdaNPm2, 1, 1);
} else if (name.compare("CdW_33r_LNP") == 0) {
	CdW_33r_LNP =  value;
	SMEFTEvol.SetCoefficient("CdWR", value * LambdaNPm2, 2, 2);
} else if (name.compare("CdB_aar_LNP") == 0) {
	CdB_aar_LNP =  value;
	SMEFTEvol.SetCoefficient("CdBR", value * LambdaNPm2, 0, 0);
	SMEFTEvol.SetCoefficient("CdBR", value * LambdaNPm2, 1, 1);
} else if (name.compare("CdB_33r_LNP") == 0) {
	CdB_33r_LNP =  value;
	SMEFTEvol.SetCoefficient("CdBR", value * LambdaNPm2, 2, 2);
} else if (name.compare("CeW_aar_LNP") == 0) {
	CeW_aar_LNP =  value;
	SMEFTEvol.SetCoefficient("CeWR", value * LambdaNPm2, 0, 0);
	SMEFTEvol.SetCoefficient("CeWR", value * LambdaNPm2, 1, 1);
} else if (name.compare("CeW_33r_LNP") == 0) {
	CeW_33r_LNP =  value;
	SMEFTEvol.SetCoefficient("CeWR", value * LambdaNPm2, 2, 2);
} else if (name.compare("CeB_aar_LNP") == 0) {
	CeB_aar_LNP =  value;
	SMEFTEvol.SetCoefficient("CeBR", value * LambdaNPm2, 0, 0);
	SMEFTEvol.SetCoefficient("CeBR", value * LambdaNPm2, 1, 1);
} else if (name.compare("CeB_33r_LNP") == 0) {
	CeB_33r_LNP =  value;
	SMEFTEvol.SetCoefficient("CeBR", value * LambdaNPm2, 2, 2);
} else if (name.compare("Cll_aabbr_LNP") == 0) {
	Cll_aabbr_LNP =  value;
	SMEFTEvol.SetCoefficient("CllR", value * LambdaNPm2, 0, 0, 0, 0);
	SMEFTEvol.SetCoefficient("CllR", value * LambdaNPm2, 0, 0, 1, 1);
	SMEFTEvol.SetCoefficient("CllR", value * LambdaNPm2, 0, 1, 1, 0);
	SMEFTEvol.SetCoefficient("CllR", value * LambdaNPm2, 1, 1, 1, 1);
} else if (name.compare("Cll_aa33r_LNP") == 0) {
	Cll_aa33r_LNP =  value;
	SMEFTEvol.SetCoefficient("CllR", value * LambdaNPm2, 0, 0, 2, 2);
	SMEFTEvol.SetCoefficient("CllR", value * LambdaNPm2, 1, 1, 2, 2);
} else if (name.compare("Cll_33aar_LNP") == 0) {
	Cll_33aar_LNP =  value;
} else if (name.compare("Cll_3aa3r_LNP") == 0) {
	Cll_3aa3r_LNP =  value;
} else if (name.compare("Cll_3333r_LNP") == 0) {
	Cll_3333r_LNP =  value;
	SMEFTEvol.SetCoefficient("CllR", value * LambdaNPm2, 2, 2, 2, 2);
} else if (name.compare("Clq1_aabbr_LNP") == 0) {
	Clq1_aabbr_LNP =  value;
	SMEFTEvol.SetCoefficient("Clq1R", value * LambdaNPm2, 0, 0, 0, 0);
	SMEFTEvol.SetCoefficient("Clq1R", value * LambdaNPm2, 0, 0, 1, 1);
	SMEFTEvol.SetCoefficient("Clq1R", value * LambdaNPm2, 1, 1, 0, 0);
} else if (name.compare("Clq1_aa33r_LNP") == 0) {
	Clq1_aa33r_LNP =  value;
	SMEFTEvol.SetCoefficient("Clq1R", value * LambdaNPm2, 0, 0, 2, 2);
	SMEFTEvol.SetCoefficient("Clq1R", value * LambdaNPm2, 1, 1, 2, 2);
} else if (name.compare("Clq1_33aar_LNP") == 0) {
	Clq1_33aar_LNP =  value;
	SMEFTEvol.SetCoefficient("Clq1R", value * LambdaNPm2, 2, 2, 0, 0);
} else if (name.compare("Clq1_3333r_LNP") == 0) {
	Clq1_3333r_LNP =  value;
	SMEFTEvol.SetCoefficient("Clq1R", value * LambdaNPm2, 2, 2, 2, 2);
} else if (name.compare("Clq3_aabbr_LNP") == 0) {
	Clq3_aabbr_LNP =  value;
	SMEFTEvol.SetCoefficient("Clq3R", value * LambdaNPm2, 0, 0, 0, 0);
	SMEFTEvol.SetCoefficient("Clq3R", value * LambdaNPm2, 0, 0, 1, 1);
	SMEFTEvol.SetCoefficient("Clq3R", value * LambdaNPm2, 1, 1, 0, 0);
} else if (name.compare("Clq3_aa33r_LNP") == 0) {
	Clq3_aa33r_LNP =  value;
	SMEFTEvol.SetCoefficient("Clq3R", value * LambdaNPm2, 0, 0, 2, 2);
	SMEFTEvol.SetCoefficient("Clq3R", value * LambdaNPm2, 1, 1, 2, 2);
} else if (name.compare("Clq3_33aar_LNP") == 0) {
	Clq3_33aar_LNP =  value;
	SMEFTEvol.SetCoefficient("Clq3R", value * LambdaNPm2, 2, 2, 0, 0);
} else if (name.compare("Clq3_3333r_LNP") == 0) {
	Clq3_3333r_LNP =  value;
	SMEFTEvol.SetCoefficient("Clq3R", value * LambdaNPm2, 2, 2, 2, 2);
} else if (name.compare("Cee_aabbr_LNP") == 0) {
	Cee_aabbr_LNP =  value;
	SMEFTEvol.SetCoefficient("CeeR", value * LambdaNPm2, 0, 0, 0, 0);
	SMEFTEvol.SetCoefficient("CeeR", value * LambdaNPm2, 0, 0, 1, 1);
	SMEFTEvol.SetCoefficient("CeeR", value * LambdaNPm2, 1, 1, 1, 1);
} else if (name.compare("Cee_aa33r_LNP") == 0) {
	Cee_aa33r_LNP =  value;
	SMEFTEvol.SetCoefficient("CeeR", value * LambdaNPm2, 0, 0, 2, 2);
	SMEFTEvol.SetCoefficient("CeeR", value * LambdaNPm2, 1, 1, 2, 2);
} else if (name.compare("Cee_33aar_LNP") == 0) {
	Cee_33aar_LNP =  value;
} else if (name.compare("Cee_3aa3r_LNP") == 0) {
	Cee_3aa3r_LNP =  value;
} else if (name.compare("Cee_3333r_LNP") == 0) {
	Cee_3333r_LNP =  value;
	SMEFTEvol.SetCoefficient("CeeR", value * LambdaNPm2, 2, 2, 2, 2);
} else if (name.compare("Ceu_aabbr_LNP") == 0) {
	Ceu_aabbr_LNP =  value;
	SMEFTEvol.SetCoefficient("CeuR", value * LambdaNPm2, 0, 0, 0, 0);
	SMEFTEvol.SetCoefficient("CeuR", value * LambdaNPm2, 0, 0, 1, 1);
	SMEFTEvol.SetCoefficient("CeuR", value * LambdaNPm2, 1, 1, 0, 0);
} else if (name.compare("Ceu_aa33r_LNP") == 0) {
	Ceu_aa33r_LNP =  value;
	SMEFTEvol.SetCoefficient("CeuR", value * LambdaNPm2, 0, 0, 2, 2);
	SMEFTEvol.SetCoefficient("CeuR", value * LambdaNPm2, 1, 1, 2, 2);
} else if (name.compare("Ceu_33aar_LNP") == 0) {
	Ceu_33aar_LNP =  value;
	SMEFTEvol.SetCoefficient("CeuR", value * LambdaNPm2, 2, 2, 0, 0);
} else if (name.compare("Ceu_3333r_LNP") == 0) {
	Ceu_3333r_LNP =  value;
	SMEFTEvol.SetCoefficient("CeuR", value * LambdaNPm2, 2, 2, 2, 2);
} else if (name.compare("Ced_aabbr_LNP") == 0) {
	Ced_aabbr_LNP =  value;
	SMEFTEvol.SetCoefficient("CedR", value * LambdaNPm2, 0, 0, 0, 0);
	SMEFTEvol.SetCoefficient("CedR", value * LambdaNPm2, 0, 0, 1, 1);
	SMEFTEvol.SetCoefficient("CedR", value * LambdaNPm2, 1, 1, 0, 0);
} else if (name.compare("Ced_aa33r_LNP") == 0) {
	Ced_aa33r_LNP =  value;
	SMEFTEvol.SetCoefficient("CedR", value * LambdaNPm2, 0, 0, 2, 2);
	SMEFTEvol.SetCoefficient("CedR", value * LambdaNPm2, 1, 1, 2, 2);
} else if (name.compare("Ced_33aar_LNP") == 0) {
	Ced_33aar_LNP =  value;
	SMEFTEvol.SetCoefficient("CedR", value * LambdaNPm2, 2, 2, 0, 0);
} else if (name.compare("Ced_3333r_LNP") == 0) {
	Ced_3333r_LNP =  value;
	SMEFTEvol.SetCoefficient("CedR", value * LambdaNPm2, 2, 2, 2, 2);
} else if (name.compare("Cle_aabbr_LNP") == 0) {
	Cle_aabbr_LNP =  value;
	SMEFTEvol.SetCoefficient("CleR", value * LambdaNPm2, 0, 0, 0, 0);
	SMEFTEvol.SetCoefficient("CleR", value * LambdaNPm2, 0, 0, 1, 1);
	SMEFTEvol.SetCoefficient("CleR", value * LambdaNPm2, 1, 1, 0, 0);
} else if (name.compare("Cle_aa33r_LNP") == 0) {
	Cle_aa33r_LNP =  value;
	SMEFTEvol.SetCoefficient("CleR", value * LambdaNPm2, 0, 0, 2, 2);
	SMEFTEvol.SetCoefficient("CleR", value * LambdaNPm2, 1, 1, 2, 2);
} else if (name.compare("Cle_33aar_LNP") == 0) {
	Cle_33aar_LNP =  value;
	SMEFTEvol.SetCoefficient("CleR", value * LambdaNPm2, 2, 2, 0, 0);
} else if (name.compare("Cle_3333r_LNP") == 0) {
	Cle_3333r_LNP =  value;
	SMEFTEvol.SetCoefficient("CleR", value * LambdaNPm2, 2, 2, 2, 2);
} else if (name.compare("Clu_aabbr_LNP") == 0) {
	Clu_aabbr_LNP =  value;
	SMEFTEvol.SetCoefficient("CluR", value * LambdaNPm2, 0, 0, 0, 0);
	SMEFTEvol.SetCoefficient("CluR", value * LambdaNPm2, 0, 0, 1, 1);
	SMEFTEvol.SetCoefficient("CluR", value * LambdaNPm2, 1, 1, 0, 0);
} else if (name.compare("Clu_aa33r_LNP") == 0) {
	Clu_aa33r_LNP =  value;
	SMEFTEvol.SetCoefficient("CluR", value * LambdaNPm2, 0, 0, 2, 2);
	SMEFTEvol.SetCoefficient("CluR", value * LambdaNPm2, 1, 1, 2, 2);
} else if (name.compare("Clu_33aar_LNP") == 0) {
	Clu_33aar_LNP =  value;
	SMEFTEvol.SetCoefficient("CluR", value * LambdaNPm2, 2, 2, 0, 0);
} else if (name.compare("Clu_3333r_LNP") == 0) {
	Clu_3333r_LNP =  value;
	SMEFTEvol.SetCoefficient("CluR", value * LambdaNPm2, 2, 2, 2, 2);
} else if (name.compare("Cld_aabbr_LNP") == 0) {
	Cld_aabbr_LNP =  value;
	SMEFTEvol.SetCoefficient("CldR", value * LambdaNPm2, 0, 0, 0, 0);
	SMEFTEvol.SetCoefficient("CldR", value * LambdaNPm2, 0, 0, 1, 1);
	SMEFTEvol.SetCoefficient("CldR", value * LambdaNPm2, 1, 1, 0, 0);
} else if (name.compare("Cld_aa33r_LNP") == 0) {
	Cld_aa33r_LNP =  value;
	SMEFTEvol.SetCoefficient("CldR", value * LambdaNPm2, 0, 0, 2, 2);
	SMEFTEvol.SetCoefficient("CldR", value * LambdaNPm2, 1, 1, 2, 2);
} else if (name.compare("Cld_33aar_LNP") == 0) {
	Cld_33aar_LNP =  value;
	SMEFTEvol.SetCoefficient("CldR", value * LambdaNPm2, 2, 2, 0, 0);
} else if (name.compare("Cld_3333r_LNP") == 0) {
	Cld_3333r_LNP =  value;
	SMEFTEvol.SetCoefficient("CldR", value * LambdaNPm2, 2, 2, 2, 2);
} else if (name.compare("Cqe_aabbr_LNP") == 0) {
	Cqe_aabbr_LNP =  value;
	SMEFTEvol.SetCoefficient("CqeR", value * LambdaNPm2, 0, 0, 0, 0);
	SMEFTEvol.SetCoefficient("CqeR", value * LambdaNPm2, 0, 0, 1, 1);
	SMEFTEvol.SetCoefficient("CqeR", value * LambdaNPm2, 1, 1, 0, 0);
} else if (name.compare("Cqe_aa33r_LNP") == 0) {
	Cqe_aa33r_LNP =  value;
	SMEFTEvol.SetCoefficient("CqeR", value * LambdaNPm2, 0, 0, 2, 2);
	SMEFTEvol.SetCoefficient("CqeR", value * LambdaNPm2, 1, 1, 2, 2);
} else if (name.compare("Cqe_33aar_LNP") == 0) {
	Cqe_33aar_LNP =  value;
	SMEFTEvol.SetCoefficient("CqeR", value * LambdaNPm2, 2, 2, 0, 0);
} else if (name.compare("Cqe_3333r_LNP") == 0) {
	Cqe_3333r_LNP =  value;
	SMEFTEvol.SetCoefficient("CqeR", value * LambdaNPm2, 2, 2, 2, 2);
} else if (name.compare("Cledq_aabbr_LNP") == 0) {
	Cledq_aabbr_LNP =  value;
	SMEFTEvol.SetCoefficient("CledqR", value * LambdaNPm2, 0, 0, 0, 0);
	SMEFTEvol.SetCoefficient("CledqR", value * LambdaNPm2, 0, 0, 1, 1);
	SMEFTEvol.SetCoefficient("CledqR", value * LambdaNPm2, 1, 1, 0, 0);
	SMEFTEvol.SetCoefficient("CledqR", value * LambdaNPm2, 1, 1, 1, 1);
} else if (name.compare("Cledq_aa33r_LNP") == 0) {
	Cledq_aa33r_LNP =  value;
	SMEFTEvol.SetCoefficient("CledqR", value * LambdaNPm2, 0, 0, 2, 2);
	SMEFTEvol.SetCoefficient("CledqR", value * LambdaNPm2, 1, 1, 2, 2);
} else if (name.compare("Cledq_33aar_LNP") == 0) {
	Cledq_33aar_LNP =  value;
} else if (name.compare("Cledq_3333r_LNP") == 0) {
	Cledq_3333r_LNP =  value;
	SMEFTEvol.SetCoefficient("CledqR", value * LambdaNPm2, 2, 2, 2, 2);
} else if (name.compare("Cqq1_aabbr_LNP") == 0) {
	Cqq1_aabbr_LNP =  value;
	SMEFTEvol.SetCoefficient("Cqq1R", value * LambdaNPm2, 0, 0, 0, 0);
	SMEFTEvol.SetCoefficient("Cqq1R", value * LambdaNPm2, 0, 0, 1, 1);
	SMEFTEvol.SetCoefficient("Cqq1R", value * LambdaNPm2, 0, 1, 1, 0);
	SMEFTEvol.SetCoefficient("Cqq1R", value * LambdaNPm2, 1, 1, 1, 1);
} else if (name.compare("Cqq1_aa33r_LNP") == 0) {
	Cqq1_aa33r_LNP =  value;
	SMEFTEvol.SetCoefficient("Cqq1R", value * LambdaNPm2, 0, 0, 2, 2);
	SMEFTEvol.SetCoefficient("Cqq1R", value * LambdaNPm2, 1, 1, 2, 2);
} else if (name.compare("Cqq1_33aar_LNP") == 0) {
	Cqq1_33aar_LNP =  value;
} else if (name.compare("Cqq1_3aa3r_LNP") == 0) {
	Cqq1_3aa3r_LNP =  value;
} else if (name.compare("Cqq1_3333r_LNP") == 0) {
	Cqq1_3333r_LNP =  value;
	SMEFTEvol.SetCoefficient("Cqq1R", value * LambdaNPm2, 2, 2, 2, 2);
} else if (name.compare("Cqq3_aabbr_LNP") == 0) {
	Cqq3_aabbr_LNP =  value;
	SMEFTEvol.SetCoefficient("Cqq3R", value * LambdaNPm2, 0, 0, 0, 0);
	SMEFTEvol.SetCoefficient("Cqq3R", value * LambdaNPm2, 0, 0, 1, 1);
	SMEFTEvol.SetCoefficient("Cqq3R", value * LambdaNPm2, 0, 1, 1, 0);
	SMEFTEvol.SetCoefficient("Cqq3R", value * LambdaNPm2, 1, 1, 1, 1);
} else if (name.compare("Cqq3_aa33r_LNP") == 0) {
	Cqq3_aa33r_LNP =  value;
	SMEFTEvol.SetCoefficient("Cqq3R", value * LambdaNPm2, 0, 0, 2, 2);
	SMEFTEvol.SetCoefficient("Cqq3R", value * LambdaNPm2, 1, 1, 2, 2);
} else if (name.compare("Cqq3_33aar_LNP") == 0) {
	Cqq3_33aar_LNP =  value;
} else if (name.compare("Cqq3_3aa3r_LNP") == 0) {
	Cqq3_3aa3r_LNP =  value;
} else if (name.compare("Cqq3_3333r_LNP") == 0) {
	Cqq3_3333r_LNP =  value;
	SMEFTEvol.SetCoefficient("Cqq3R", value * LambdaNPm2, 2, 2, 2, 2);
} else if (name.compare("Cuu_aabbr_LNP") == 0) {
	Cuu_aabbr_LNP =  value;
	SMEFTEvol.SetCoefficient("CuuR", value * LambdaNPm2, 0, 0, 0, 0);
	SMEFTEvol.SetCoefficient("CuuR", value * LambdaNPm2, 0, 0, 1, 1);
	SMEFTEvol.SetCoefficient("CuuR", value * LambdaNPm2, 0, 1, 1, 0);
	SMEFTEvol.SetCoefficient("CuuR", value * LambdaNPm2, 1, 1, 1, 1);
} else if (name.compare("Cuu_aa33r_LNP") == 0) {
	Cuu_aa33r_LNP =  value;
	SMEFTEvol.SetCoefficient("CuuR", value * LambdaNPm2, 0, 0, 2, 2);
	SMEFTEvol.SetCoefficient("CuuR", value * LambdaNPm2, 1, 1, 2, 2);
} else if (name.compare("Cuu_33aar_LNP") == 0) {
	Cuu_33aar_LNP =  value;
} else if (name.compare("Cuu_3aa3r_LNP") == 0) {
	Cuu_3aa3r_LNP =  value;
} else if (name.compare("Cuu_3333r_LNP") == 0) {
	Cuu_3333r_LNP =  value;
	SMEFTEvol.SetCoefficient("CuuR", value * LambdaNPm2, 2, 2, 2, 2);
} else if (name.compare("Cdd_aabbr_LNP") == 0) {
	Cdd_aabbr_LNP =  value;
	SMEFTEvol.SetCoefficient("CddR", value * LambdaNPm2, 0, 0, 0, 0);
	SMEFTEvol.SetCoefficient("CddR", value * LambdaNPm2, 0, 0, 1, 1);
	SMEFTEvol.SetCoefficient("CddR", value * LambdaNPm2, 0, 1, 1, 0);
	SMEFTEvol.SetCoefficient("CddR", value * LambdaNPm2, 1, 1, 1, 1);
} else if (name.compare("Cdd_aa33r_LNP") == 0) {
	Cdd_aa33r_LNP =  value;
	SMEFTEvol.SetCoefficient("CddR", value * LambdaNPm2, 0, 0, 2, 2);
	SMEFTEvol.SetCoefficient("CddR", value * LambdaNPm2, 1, 1, 2, 2);
} else if (name.compare("Cdd_33aar_LNP") == 0) {
	Cdd_33aar_LNP =  value;
} else if (name.compare("Cdd_3aa3r_LNP") == 0) {
	Cdd_3aa3r_LNP =  value;
} else if (name.compare("Cdd_3333r_LNP") == 0) {
	Cdd_3333r_LNP =  value;
	SMEFTEvol.SetCoefficient("CddR", value * LambdaNPm2, 2, 2, 2, 2);
} else if (name.compare("Cud1_aabbr_LNP") == 0) {
	Cud1_aabbr_LNP =  value;
	SMEFTEvol.SetCoefficient("Cud1R", value * LambdaNPm2, 0, 0, 0, 0);
	SMEFTEvol.SetCoefficient("Cud1R", value * LambdaNPm2, 0, 0, 1, 1);
	SMEFTEvol.SetCoefficient("Cud1R", value * LambdaNPm2, 1, 1, 0, 0);
} else if (name.compare("Cud1_aa33r_LNP") == 0) {
	Cud1_aa33r_LNP =  value;
	SMEFTEvol.SetCoefficient("Cud1R", value * LambdaNPm2, 0, 0, 2, 2);
	SMEFTEvol.SetCoefficient("Cud1R", value * LambdaNPm2, 1, 1, 2, 2);
} else if (name.compare("Cud1_33aar_LNP") == 0) {
	Cud1_33aar_LNP =  value;
	SMEFTEvol.SetCoefficient("Cud1R", value * LambdaNPm2, 2, 2, 0, 0);
} else if (name.compare("Cud1_3333r_LNP") == 0) {
	Cud1_3333r_LNP =  value;
	SMEFTEvol.SetCoefficient("Cud1R", value * LambdaNPm2, 2, 2, 2, 2);
} else if (name.compare("Cud8_aabbr_LNP") == 0) {
	Cud8_aabbr_LNP =  value;
	SMEFTEvol.SetCoefficient("Cud8R", value * LambdaNPm2, 0, 0, 0, 0);
	SMEFTEvol.SetCoefficient("Cud8R", value * LambdaNPm2, 0, 0, 1, 1);
	SMEFTEvol.SetCoefficient("Cud8R", value * LambdaNPm2, 1, 1, 0, 0);
} else if (name.compare("Cud8_aa33r_LNP") == 0) {
	Cud8_aa33r_LNP =  value;
	SMEFTEvol.SetCoefficient("Cud8R", value * LambdaNPm2, 0, 0, 2, 2);
	SMEFTEvol.SetCoefficient("Cud8R", value * LambdaNPm2, 1, 1, 2, 2);
} else if (name.compare("Cud8_33aar_LNP") == 0) {
	Cud8_33aar_LNP =  value;
	SMEFTEvol.SetCoefficient("Cud8R", value * LambdaNPm2, 2, 2, 0, 0);
} else if (name.compare("Cud8_3333r_LNP") == 0) {
	Cud8_3333r_LNP =  value;
	SMEFTEvol.SetCoefficient("Cud8R", value * LambdaNPm2, 2, 2, 2, 2);
} else if (name.compare("Cqu1_aabbr_LNP") == 0) {
	Cqu1_aabbr_LNP =  value;
	SMEFTEvol.SetCoefficient("Cqu1R", value * LambdaNPm2, 0, 0, 0, 0);
	SMEFTEvol.SetCoefficient("Cqu1R", value * LambdaNPm2, 0, 0, 1, 1);
	SMEFTEvol.SetCoefficient("Cqu1R", value * LambdaNPm2, 1, 1, 0, 0);
} else if (name.compare("Cqu1_aa33r_LNP") == 0) {
	Cqu1_aa33r_LNP =  value;
	SMEFTEvol.SetCoefficient("Cqu1R", value * LambdaNPm2, 0, 0, 2, 2);
	SMEFTEvol.SetCoefficient("Cqu1R", value * LambdaNPm2, 1, 1, 2, 2);
} else if (name.compare("Cqu1_33aar_LNP") == 0) {
	Cqu1_33aar_LNP =  value;
	SMEFTEvol.SetCoefficient("Cqu1R", value * LambdaNPm2, 2, 2, 0, 0);
} else if (name.compare("Cqu1_3333r_LNP") == 0) {
	Cqu1_3333r_LNP =  value;
	SMEFTEvol.SetCoefficient("Cqu1R", value * LambdaNPm2, 2, 2, 2, 2);
} else if (name.compare("Cqu8_aabbr_LNP") == 0) {
	Cqu8_aabbr_LNP =  value;
	SMEFTEvol.SetCoefficient("Cqu8R", value * LambdaNPm2, 0, 0, 0, 0);
	SMEFTEvol.SetCoefficient("Cqu8R", value * LambdaNPm2, 0, 0, 1, 1);
	SMEFTEvol.SetCoefficient("Cqu8R", value * LambdaNPm2, 1, 1, 0, 0);
} else if (name.compare("Cqu8_aa33r_LNP") == 0) {
	Cqu8_aa33r_LNP =  value;
	SMEFTEvol.SetCoefficient("Cqu8R", value * LambdaNPm2, 0, 0, 2, 2);
	SMEFTEvol.SetCoefficient("Cqu8R", value * LambdaNPm2, 1, 1, 2, 2);
} else if (name.compare("Cqu8_33aar_LNP") == 0) {
	Cqu8_33aar_LNP =  value;
	SMEFTEvol.SetCoefficient("Cqu8R", value * LambdaNPm2, 2, 2, 0, 0);
} else if (name.compare("Cqu8_3333r_LNP") == 0) {
	Cqu8_3333r_LNP =  value;
	SMEFTEvol.SetCoefficient("Cqu8R", value * LambdaNPm2, 2, 2, 2, 2);
} else if (name.compare("Cqd1_aabbr_LNP") == 0) {
	Cqd1_aabbr_LNP =  value;
	SMEFTEvol.SetCoefficient("Cqd1R", value * LambdaNPm2, 0, 0, 0, 0);
	SMEFTEvol.SetCoefficient("Cqd1R", value * LambdaNPm2, 0, 0, 1, 1);
	SMEFTEvol.SetCoefficient("Cqd1R", value * LambdaNPm2, 1, 1, 0, 0);
} else if (name.compare("Cqd1_aa33r_LNP") == 0) {
	Cqd1_aa33r_LNP =  value;
	SMEFTEvol.SetCoefficient("Cqd1R", value * LambdaNPm2, 0, 0, 2, 2);
	SMEFTEvol.SetCoefficient("Cqd1R", value * LambdaNPm2, 1, 1, 2, 2);
} else if (name.compare("Cqd1_33aar_LNP") == 0) {
	Cqd1_33aar_LNP =  value;
	SMEFTEvol.SetCoefficient("Cqd1R", value * LambdaNPm2, 2, 2, 0, 0);
} else if (name.compare("Cqd1_3333r_LNP") == 0) {
	Cqd1_3333r_LNP =  value;
	SMEFTEvol.SetCoefficient("Cqd1R", value * LambdaNPm2, 2, 2, 2, 2);
} else if (name.compare("Cqd8_aabbr_LNP") == 0) {
	Cqd8_aabbr_LNP =  value;
	SMEFTEvol.SetCoefficient("Cqd8R", value * LambdaNPm2, 0, 0, 0, 0);
	SMEFTEvol.SetCoefficient("Cqd8R", value * LambdaNPm2, 0, 0, 1, 1);
	SMEFTEvol.SetCoefficient("Cqd8R", value * LambdaNPm2, 1, 1, 0, 0);
} else if (name.compare("Cqd8_aa33r_LNP") == 0) {
	Cqd8_aa33r_LNP =  value;
	SMEFTEvol.SetCoefficient("Cqd8R", value * LambdaNPm2, 0, 0, 2, 2);
	SMEFTEvol.SetCoefficient("Cqd8R", value * LambdaNPm2, 1, 1, 2, 2);
} else if (name.compare("Cqd8_33aar_LNP") == 0) {
	Cqd8_33aar_LNP =  value;
	SMEFTEvol.SetCoefficient("Cqd8R", value * LambdaNPm2, 2, 2, 0, 0);
} else if (name.compare("Cqd8_3333r_LNP") == 0) {
	Cqd8_3333r_LNP =  value;
	SMEFTEvol.SetCoefficient("Cqd8R", value * LambdaNPm2, 2, 2, 2, 2);
} else if (name.compare("Cquqd1_aabbr_LNP") == 0) {
	Cquqd1_aabbr_LNP =  value;
	SMEFTEvol.SetCoefficient("Cquqd1R", value * LambdaNPm2, 0, 0, 0, 0);
	SMEFTEvol.SetCoefficient("Cquqd1R", value * LambdaNPm2, 0, 0, 1, 1);
	SMEFTEvol.SetCoefficient("Cquqd1R", value * LambdaNPm2, 1, 1, 0, 0);
	SMEFTEvol.SetCoefficient("Cquqd1R", value * LambdaNPm2, 1, 1, 1, 1);
} else if (name.compare("Cquqd1_aa33r_LNP") == 0) {
	Cquqd1_aa33r_LNP =  value;
	SMEFTEvol.SetCoefficient("Cquqd1R", value * LambdaNPm2, 0, 0, 2, 2);
	SMEFTEvol.SetCoefficient("Cquqd1R", value * LambdaNPm2, 1, 1, 2, 2);
} else if (name.compare("Cquqd1_33aar_LNP") == 0) {
	Cquqd1_33aar_LNP =  value;
} else if (name.compare("Cquqd1_3333r_LNP") == 0) {
	Cquqd1_3333r_LNP =  value;
	SMEFTEvol.SetCoefficient("Cquqd1R", value * LambdaNPm2, 2, 2, 2, 2);
} else if (name.compare("Cquqd8_aabbr_LNP") == 0) {
	Cquqd8_aabbr_LNP =  value;
	SMEFTEvol.SetCoefficient("Cquqd8R", value * LambdaNPm2, 0, 0, 0, 0);
	SMEFTEvol.SetCoefficient("Cquqd8R", value * LambdaNPm2, 0, 0, 1, 1);
	SMEFTEvol.SetCoefficient("Cquqd8R", value * LambdaNPm2, 1, 1, 0, 0);
	SMEFTEvol.SetCoefficient("Cquqd8R", value * LambdaNPm2, 1, 1, 1, 1);
} else if (name.compare("Cquqd8_aa33r_LNP") == 0) {
	Cquqd8_aa33r_LNP =  value;
	SMEFTEvol.SetCoefficient("Cquqd8R", value * LambdaNPm2, 0, 0, 2, 2);
	SMEFTEvol.SetCoefficient("Cquqd8R", value * LambdaNPm2, 1, 1, 2, 2);
} else if (name.compare("Cquqd8_33aar_LNP") == 0) {
	Cquqd8_33aar_LNP =  value;
} else if (name.compare("Cquqd8_3333r_LNP") == 0) {
	Cquqd8_3333r_LNP =  value;
	SMEFTEvol.SetCoefficient("Cquqd8R", value * LambdaNPm2, 2, 2, 2, 2);
} else if (name.compare("Clequ1_aabbr_LNP") == 0) {
	Clequ1_aabbr_LNP =  value;
	SMEFTEvol.SetCoefficient("Clequ1R", value * LambdaNPm2, 0, 0, 0, 0);
	SMEFTEvol.SetCoefficient("Clequ1R", value * LambdaNPm2, 0, 0, 1, 1);
	SMEFTEvol.SetCoefficient("Clequ1R", value * LambdaNPm2, 1, 1, 0, 0);
	SMEFTEvol.SetCoefficient("Clequ1R", value * LambdaNPm2, 1, 1, 1, 1);
} else if (name.compare("Clequ1_aa33r_LNP") == 0) {
	Clequ1_aa33r_LNP =  value;
	SMEFTEvol.SetCoefficient("Clequ1R", value * LambdaNPm2, 0, 0, 2, 2);
	SMEFTEvol.SetCoefficient("Clequ1R", value * LambdaNPm2, 1, 1, 2, 2);
} else if (name.compare("Clequ1_33aar_LNP") == 0) {
	Clequ1_33aar_LNP =  value;
} else if (name.compare("Clequ1_3333r_LNP") == 0) {
	Clequ1_3333r_LNP =  value;
	SMEFTEvol.SetCoefficient("Clequ1R", value * LambdaNPm2, 2, 2, 2, 2);
} else if (name.compare("Clequ3_aabbr_LNP") == 0) {
	Clequ3_aabbr_LNP =  value;
	SMEFTEvol.SetCoefficient("Clequ3R", value * LambdaNPm2, 0, 0, 0, 0);
	SMEFTEvol.SetCoefficient("Clequ3R", value * LambdaNPm2, 0, 0, 1, 1);
	SMEFTEvol.SetCoefficient("Clequ3R", value * LambdaNPm2, 1, 1, 0, 0);
	SMEFTEvol.SetCoefficient("Clequ3R", value * LambdaNPm2, 1, 1, 1, 1);
} else if (name.compare("Clequ3_aa33r_LNP") == 0) {
	Clequ3_aa33r_LNP =  value;
	SMEFTEvol.SetCoefficient("Clequ3R", value * LambdaNPm2, 0, 0, 2, 2);
	SMEFTEvol.SetCoefficient("Clequ3R", value * LambdaNPm2, 1, 1, 2, 2);
} else if (name.compare("Clequ3_33aar_LNP") == 0) {
	Clequ3_33aar_LNP =  value;
} else if (name.compare("Clequ3_3333r_LNP") == 0) {
	Clequ3_3333r_LNP =  value;
	SMEFTEvol.SetCoefficient("Clequ3R", value * LambdaNPm2, 2, 2, 2, 2);
} else if (name.compare("Lambda_NP") == 0)
        Lambda_NP = value;
    //Let's go directly for the NPbase setParameter, we should change it to the NPSMEFTd6General removing reading the unnecessary parameters
    else
        NPbase::setParameter(name, value);

}