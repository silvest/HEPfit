/* 
 * Copyright (C) 2012 HEPfit Collaboration
 *
 *
 * For the licensing terms see doc/COPYING.
 * 
 * Created on 18 September 2023, 16:23
 * 
 */


#include "NPSMEFTd6MFP.h"



std::string NPSMEFTd6MFP::NPSMEFTd6MFPVars[NNPSMEFTd6MFPVars] = {
    "CG_LNP", "CW_LNP", "CHG_LNP", "CHW_LNP", "CHB_LNP", "CHWB_LNP", "CHD_LNP", "CHbox_LNP", "CH_LNP",
    "CHl1_11r_LNP", "CHl1_22r_LNP", "CHl1_33r_LNP", "CHl3_11r_LNP", "CHl3_22r_LNP", "CHl3_33r_LNP", 
    "CHe_11r_LNP", "CHe_22r_LNP", "CHe_33r_LNP",
    "CHq1_aar_LNP", "CHq1_33r_LNP", "CHq3_aar_LNP", "CHq3_33r_LNP", 
    "CHu_11r_LNP", "CHu_22r_LNP", "CHu_33r_LNP", "CHd_11r_LNP", "CHd_22r_LNP", "CHd_33r_LNP", "CHud_11r_LNP", 
    "CuH_33r_LNP",
    "CuG_33r_LNP", "CuW_33r_LNP", "CuB_33r_LNP", 
    "Cll_1111r_LNP", "Cll_2222r_LNP", "Cll_3333r_LNP", "Cll_1122r_LNP", "Cll_1133r_LNP", "Cll_2233r_LNP", 
    "Cll_1221r_LNP", "Cll_1331r_LNP", "Cll_2332r_LNP", "Cll_1213r_LNP",
    "Clq1_11aar_LNP", "Clq1_22aar_LNP", "Clq1_33aar_LNP", "Clq1_1133r_LNP", "Clq1_2233r_LNP", "Clq1_3333r_LNP",
    "Clq3_11aar_LNP", "Clq3_22aar_LNP", "Clq3_33aar_LNP", "Clq3_1133r_LNP", "Clq3_2233r_LNP", "Clq3_3333r_LNP",
    "Cee_1111r_LNP", "Cee_2222r_LNP", "Cee_3333r_LNP", "Cee_1122r_LNP", "Cee_1133r_LNP", "Cee_2233r_LNP",
    "Ceu_1111r_LNP", "Ceu_2222r_LNP", "Ceu_3333r_LNP", "Ceu_1122r_LNP", "Ceu_1133r_LNP", "Ceu_2233r_LNP", 
    "Ceu_2211r_LNP", "Ceu_3311r_LNP", "Ceu_3322r_LNP", "Ceu_1332r_LNP",
    "Ced_1111r_LNP", "Ced_2222r_LNP", "Ced_3333r_LNP", "Ced_1122r_LNP", "Ced_1133r_LNP", "Ced_2233r_LNP", 
    "Ced_2211r_LNP", "Ced_3311r_LNP", "Ced_3322r_LNP", 
    "Cle_1111r_LNP", "Cle_2222r_LNP", "Cle_3333r_LNP", "Cle_1122r_LNP", "Cle_1133r_LNP", "Cle_2233r_LNP", 
    "Cle_2211r_LNP", "Cle_3311r_LNP", "Cle_3322r_LNP", 
    "Clu_1111r_LNP", "Clu_2222r_LNP", "Clu_3333r_LNP", "Clu_1122r_LNP", "Clu_1133r_LNP", "Clu_2233r_LNP", 
    "Clu_2211r_LNP", "Clu_3311r_LNP", "Clu_3322r_LNP", 
    "Cld_1111r_LNP", "Cld_2222r_LNP", "Cld_3333r_LNP", "Cld_1122r_LNP", "Cld_1133r_LNP", "Cld_2233r_LNP", 
    "Cld_2211r_LNP", "Cld_3311r_LNP", "Cld_3322r_LNP", "Cld_2332r_LNP",
    "Cqe_aa11r_LNP", "Cqe_aa22r_LNP", "Cqe_aa33r_LNP", "Cqe_3311r_LNP", "Cqe_3322r_LNP", "Cqe_3333r_LNP",
    "Cledq_3333r_LNP",
    "Cqq1_aabbr_LNP", "Cqq1_abbar_LNP", "Cqq1_a33ar_LNP","Cqq1_aa33r_LNP", "Cqq1_3333r_LNP",
    "Cqq3_aabbr_LNP", "Cqq3_abbar_LNP", "Cqq3_a33ar_LNP","Cqq3_aa33r_LNP", "Cqq3_3333r_LNP",
    "Cuu_1111r_LNP", "Cuu_2222r_LNP", "Cuu_3333r_LNP", "Cuu_1122r_LNP", "Cuu_1133r_LNP", "Cuu_2233r_LNP", "Cuu_1221r_LNP", "Cuu_1331r_LNP", "Cuu_2332r_LNP",
    "Cdd_1111r_LNP", "Cdd_2222r_LNP", "Cdd_3333r_LNP", "Cdd_1122r_LNP", "Cdd_1133r_LNP", "Cdd_2233r_LNP", "Cdd_1221r_LNP", "Cdd_1331r_LNP", "Cdd_2332r_LNP",
    "Cud1_1111r_LNP", "Cud1_2222r_LNP", "Cud1_3333r_LNP", "Cud1_1122r_LNP", "Cud1_2211r_LNP", "Cud1_1133r_LNP", "Cud1_3311r_LNP", "Cud1_2233r_LNP", "Cud1_3322r_LNP",
    "Cud8_1111r_LNP", "Cud8_2222r_LNP", "Cud8_3333r_LNP", "Cud8_1122r_LNP", "Cud8_2211r_LNP", "Cud8_1133r_LNP", "Cud8_3311r_LNP", "Cud8_2233r_LNP", "Cud8_3322r_LNP",
    "Cqu1_aa11r_LNP", "Cqu1_aa22r_LNP", "Cqu1_aa33r_LNP", "Cqu1_3311r_LNP", "Cqu1_3322r_LNP", "Cqu1_3333r_LNP",
    "Cqu8_aa11r_LNP", "Cqu8_aa22r_LNP", "Cqu8_aa33r_LNP", "Cqu8_3311r_LNP", "Cqu8_3322r_LNP", "Cqu8_3333r_LNP",
    "Cqd1_aa11r_LNP", "Cqd1_aa22r_LNP", "Cqd1_aa33r_LNP", "Cqd1_3311r_LNP", "Cqd1_3322r_LNP", "Cqd1_3333r_LNP",
    "Cqd8_aa11r_LNP", "Cqd8_aa22r_LNP", "Cqd8_aa33r_LNP", "Cqd8_3311r_LNP", "Cqd8_3322r_LNP", "Cqd8_3333r_LNP",
    "Lambda_NP"
};

NPSMEFTd6MFP::NPSMEFTd6MFP()
: NPSMEFTd6General() {
    setModelName("NPSMEFTd6MFP");
    SMEFTBasisFlag = "UP";
    ModelParamMap.insert(std::make_pair("CG_LNP", std::cref(CG_LNP)));
    ModelParamMap.insert(std::make_pair("CW_LNP", std::cref(CW_LNP)));
    ModelParamMap.insert(std::make_pair("CHG_LNP", std::cref(CHG_LNP)));
    ModelParamMap.insert(std::make_pair("CHW_LNP", std::cref(CHW_LNP)));
    ModelParamMap.insert(std::make_pair("CHB_LNP", std::cref(CHB_LNP)));
    ModelParamMap.insert(std::make_pair("CHWB_LNP", std::cref(CHWB_LNP)));
    ModelParamMap.insert(std::make_pair("CHD_LNP", std::cref(CHD_LNP)));
    ModelParamMap.insert(std::make_pair("CHbox_LNP", std::cref(CHbox_LNP)));
    ModelParamMap.insert(std::make_pair("CH_LNP", std::cref(CH_LNP)));
    ModelParamMap.insert(std::make_pair("CHl1_11r_LNP", std::cref(CHl1_11r_LNP)));
    ModelParamMap.insert(std::make_pair("CHl1_22r_LNP", std::cref(CHl1_22r_LNP)));
    ModelParamMap.insert(std::make_pair("CHl1_33r_LNP", std::cref(CHl1_33r_LNP)));
    ModelParamMap.insert(std::make_pair("CHl3_11r_LNP", std::cref(CHl3_11r_LNP)));
    ModelParamMap.insert(std::make_pair("CHl3_22r_LNP", std::cref(CHl3_22r_LNP)));
    ModelParamMap.insert(std::make_pair("CHl3_33r_LNP", std::cref(CHl3_33r_LNP)));
    ModelParamMap.insert(std::make_pair("CHe_11r_LNP", std::cref(CHe_11r_LNP)));
    ModelParamMap.insert(std::make_pair("CHe_22r_LNP", std::cref(CHe_22r_LNP)));
    ModelParamMap.insert(std::make_pair("CHe_33r_LNP", std::cref(CHe_33r_LNP)));
    ModelParamMap.insert(std::make_pair("CHq1_aar_LNP", std::cref(CHq1_aar_LNP)));
    ModelParamMap.insert(std::make_pair("CHq1_33r_LNP", std::cref(CHq1_33r_LNP)));
    ModelParamMap.insert(std::make_pair("CHq3_aar_LNP", std::cref(CHq3_aar_LNP)));
    ModelParamMap.insert(std::make_pair("CHq3_33r_LNP", std::cref(CHq3_33r_LNP)));
    ModelParamMap.insert(std::make_pair("CHu_11r_LNP", std::cref(CHu_11r_LNP)));
    ModelParamMap.insert(std::make_pair("CHu_22r_LNP", std::cref(CHu_22r_LNP)));
    ModelParamMap.insert(std::make_pair("CHu_33r_LNP", std::cref(CHu_33r_LNP)));
    ModelParamMap.insert(std::make_pair("CHd_11r_LNP", std::cref(CHd_11r_LNP)));
    ModelParamMap.insert(std::make_pair("CHd_22r_LNP", std::cref(CHd_22r_LNP)));
    ModelParamMap.insert(std::make_pair("CHd_33r_LNP", std::cref(CHd_33r_LNP)));
    ModelParamMap.insert(std::make_pair("CHud_11r_LNP", std::cref(CHud_11r_LNP)));
    ModelParamMap.insert(std::make_pair("CuH_33r_LNP", std::cref(CuH_33r_LNP)));
    ModelParamMap.insert(std::make_pair("CuG_33r_LNP", std::cref(CuG_33r_LNP)));
    ModelParamMap.insert(std::make_pair("CuW_33r_LNP", std::cref(CuW_33r_LNP)));
    ModelParamMap.insert(std::make_pair("CuB_33r_LNP", std::cref(CuB_33r_LNP)));
    ModelParamMap.insert(std::make_pair("Cll_1111r_LNP", std::cref(Cll_1111r_LNP)));
    ModelParamMap.insert(std::make_pair("Cll_2222r_LNP", std::cref(Cll_2222r_LNP)));
    ModelParamMap.insert(std::make_pair("Cll_3333r_LNP", std::cref(Cll_3333r_LNP)));
    ModelParamMap.insert(std::make_pair("Cll_1122r_LNP", std::cref(Cll_1122r_LNP)));
    ModelParamMap.insert(std::make_pair("Cll_1133r_LNP", std::cref(Cll_1133r_LNP)));
    ModelParamMap.insert(std::make_pair("Cll_2233r_LNP", std::cref(Cll_2233r_LNP)));
    ModelParamMap.insert(std::make_pair("Cll_1221r_LNP", std::cref(Cll_1221r_LNP)));
    ModelParamMap.insert(std::make_pair("Cll_1331r_LNP", std::cref(Cll_1331r_LNP)));
    ModelParamMap.insert(std::make_pair("Cll_2332r_LNP", std::cref(Cll_2332r_LNP)));
    ModelParamMap.insert(std::make_pair("Cll_1213r_LNP", std::cref(Cll_1213r_LNP)));
    ModelParamMap.insert(std::make_pair("Clq1_11aar_LNP", std::cref(Clq1_11aar_LNP)));
    ModelParamMap.insert(std::make_pair("Clq1_22aar_LNP", std::cref(Clq1_22aar_LNP)));
    ModelParamMap.insert(std::make_pair("Clq1_33aar_LNP", std::cref(Clq1_33aar_LNP)));
    ModelParamMap.insert(std::make_pair("Clq1_1133r_LNP", std::cref(Clq1_1133r_LNP)));
    ModelParamMap.insert(std::make_pair("Clq1_2233r_LNP", std::cref(Clq1_2233r_LNP)));
    ModelParamMap.insert(std::make_pair("Clq1_3333r_LNP", std::cref(Clq1_3333r_LNP)));
    ModelParamMap.insert(std::make_pair("Clq3_11aar_LNP", std::cref(Clq3_11aar_LNP)));
    ModelParamMap.insert(std::make_pair("Clq3_22aar_LNP", std::cref(Clq3_22aar_LNP)));
    ModelParamMap.insert(std::make_pair("Clq3_33aar_LNP", std::cref(Clq3_33aar_LNP)));
    ModelParamMap.insert(std::make_pair("Clq3_1133r_LNP", std::cref(Clq3_1133r_LNP)));
    ModelParamMap.insert(std::make_pair("Clq3_2233r_LNP", std::cref(Clq3_2233r_LNP)));
    ModelParamMap.insert(std::make_pair("Clq3_3333r_LNP", std::cref(Clq3_3333r_LNP))); 
    ModelParamMap.insert(std::make_pair("Cee_1111r_LNP", std::cref(Cee_1111r_LNP)));
    ModelParamMap.insert(std::make_pair("Cee_2222r_LNP", std::cref(Cee_2222r_LNP)));
    ModelParamMap.insert(std::make_pair("Cee_3333r_LNP", std::cref(Cee_3333r_LNP)));
    ModelParamMap.insert(std::make_pair("Cee_1122r_LNP", std::cref(Cee_1122r_LNP)));
    ModelParamMap.insert(std::make_pair("Cee_1133r_LNP", std::cref(Cee_1133r_LNP)));
    ModelParamMap.insert(std::make_pair("Cee_2233r_LNP", std::cref(Cee_2233r_LNP)));
    ModelParamMap.insert(std::make_pair("Ceu_1111r_LNP", std::cref(Ceu_1111r_LNP)));
    ModelParamMap.insert(std::make_pair("Ceu_1122r_LNP", std::cref(Ceu_1122r_LNP)));
    ModelParamMap.insert(std::make_pair("Ceu_1133r_LNP", std::cref(Ceu_1133r_LNP)));
    ModelParamMap.insert(std::make_pair("Ceu_2222r_LNP", std::cref(Ceu_2222r_LNP)));
    ModelParamMap.insert(std::make_pair("Ceu_2233r_LNP", std::cref(Ceu_2233r_LNP)));
    ModelParamMap.insert(std::make_pair("Ceu_3333r_LNP", std::cref(Ceu_3333r_LNP)));
    ModelParamMap.insert(std::make_pair("Ceu_2211r_LNP", std::cref(Ceu_2211r_LNP)));
    ModelParamMap.insert(std::make_pair("Ceu_3311r_LNP", std::cref(Ceu_3311r_LNP)));
    ModelParamMap.insert(std::make_pair("Ceu_3322r_LNP", std::cref(Ceu_3322r_LNP)));
    ModelParamMap.insert(std::make_pair("Ceu_1332r_LNP", std::cref(Ceu_1332r_LNP)));
    ModelParamMap.insert(std::make_pair("Ced_1111r_LNP", std::cref(Ced_1111r_LNP)));
    ModelParamMap.insert(std::make_pair("Ced_1122r_LNP", std::cref(Ced_1122r_LNP)));
    ModelParamMap.insert(std::make_pair("Ced_1133r_LNP", std::cref(Ced_1133r_LNP)));
    ModelParamMap.insert(std::make_pair("Ced_2222r_LNP", std::cref(Ced_2222r_LNP)));
    ModelParamMap.insert(std::make_pair("Ced_2233r_LNP", std::cref(Ced_2233r_LNP)));
    ModelParamMap.insert(std::make_pair("Ced_3333r_LNP", std::cref(Ced_3333r_LNP)));
    ModelParamMap.insert(std::make_pair("Cle_1111r_LNP", std::cref(Cle_1111r_LNP)));
    ModelParamMap.insert(std::make_pair("Cle_1122r_LNP", std::cref(Cle_1122r_LNP)));
    ModelParamMap.insert(std::make_pair("Cle_1133r_LNP", std::cref(Cle_1133r_LNP)));
    ModelParamMap.insert(std::make_pair("Cle_2222r_LNP", std::cref(Cle_2222r_LNP)));
    ModelParamMap.insert(std::make_pair("Cle_2233r_LNP", std::cref(Cle_2233r_LNP)));
    ModelParamMap.insert(std::make_pair("Cle_3333r_LNP", std::cref(Cle_3333r_LNP)));
    ModelParamMap.insert(std::make_pair("Cle_2211r_LNP", std::cref(Cle_2211r_LNP)));
    ModelParamMap.insert(std::make_pair("Cle_3311r_LNP", std::cref(Cle_3311r_LNP)));
    ModelParamMap.insert(std::make_pair("Cle_3322r_LNP", std::cref(Cle_3322r_LNP)));
    ModelParamMap.insert(std::make_pair("Clu_1111r_LNP", std::cref(Clu_1111r_LNP)));
    ModelParamMap.insert(std::make_pair("Clu_1122r_LNP", std::cref(Clu_1122r_LNP)));
    ModelParamMap.insert(std::make_pair("Clu_1133r_LNP", std::cref(Clu_1133r_LNP)));
    ModelParamMap.insert(std::make_pair("Clu_2222r_LNP", std::cref(Clu_2222r_LNP)));
    ModelParamMap.insert(std::make_pair("Clu_2233r_LNP", std::cref(Clu_2233r_LNP)));
    ModelParamMap.insert(std::make_pair("Clu_3333r_LNP", std::cref(Clu_3333r_LNP)));
    ModelParamMap.insert(std::make_pair("Clu_2211r_LNP", std::cref(Clu_2211r_LNP)));
    ModelParamMap.insert(std::make_pair("Clu_3311r_LNP", std::cref(Clu_3311r_LNP)));
    ModelParamMap.insert(std::make_pair("Clu_3322r_LNP", std::cref(Clu_3322r_LNP)));
    ModelParamMap.insert(std::make_pair("Cld_1111r_LNP", std::cref(Cld_1111r_LNP)));
    ModelParamMap.insert(std::make_pair("Cld_1122r_LNP", std::cref(Cld_1122r_LNP)));
    ModelParamMap.insert(std::make_pair("Cld_1133r_LNP", std::cref(Cld_1133r_LNP)));
    ModelParamMap.insert(std::make_pair("Cld_2222r_LNP", std::cref(Cld_2222r_LNP)));
    ModelParamMap.insert(std::make_pair("Cld_2233r_LNP", std::cref(Cld_2233r_LNP)));
    ModelParamMap.insert(std::make_pair("Cld_3333r_LNP", std::cref(Cld_3333r_LNP)));
    ModelParamMap.insert(std::make_pair("Cld_2211r_LNP", std::cref(Cld_2211r_LNP)));
    ModelParamMap.insert(std::make_pair("Cld_3311r_LNP", std::cref(Cld_3311r_LNP)));
    ModelParamMap.insert(std::make_pair("Cld_3322r_LNP", std::cref(Cld_3322r_LNP)));
    ModelParamMap.insert(std::make_pair("Cld_2332r_LNP", std::cref(Cld_2332r_LNP)));
    ModelParamMap.insert(std::make_pair("Cqe_aa11r_LNP", std::cref(Cqe_aa11r_LNP)));
    ModelParamMap.insert(std::make_pair("Cqe_aa22r_LNP", std::cref(Cqe_aa22r_LNP)));
    ModelParamMap.insert(std::make_pair("Cqe_aa33r_LNP", std::cref(Cqe_aa33r_LNP)));
    ModelParamMap.insert(std::make_pair("Cqe_3311r_LNP", std::cref(Cqe_3311r_LNP)));
    ModelParamMap.insert(std::make_pair("Cqe_3322r_LNP", std::cref(Cqe_3322r_LNP)));
    ModelParamMap.insert(std::make_pair("Cqe_3333r_LNP", std::cref(Cqe_3333r_LNP)));
    ModelParamMap.insert(std::make_pair("Cledq_3333r_LNP", std::cref(Cledq_3333r_LNP)));
    ModelParamMap.insert(std::make_pair("Cqq1_aabbr_LNP", std::cref(Cqq1_aabbr_LNP)));
    ModelParamMap.insert(std::make_pair("Cqq1_abbar_LNP", std::cref(Cqq1_abbar_LNP)));  
    ModelParamMap.insert(std::make_pair("Cqq1_a33ar_LNP", std::cref(Cqq1_a33ar_LNP)));
    ModelParamMap.insert(std::make_pair("Cqq1_aa33r_LNP", std::cref(Cqq1_aa33r_LNP)));
    ModelParamMap.insert(std::make_pair("Cqq1_3333r_LNP", std::cref(Cqq1_3333r_LNP)));
    ModelParamMap.insert(std::make_pair("Cqq3_aabbr_LNP", std::cref(Cqq3_aabbr_LNP)));
    ModelParamMap.insert(std::make_pair("Cqq3_abbar_LNP", std::cref(Cqq3_abbar_LNP)));
    ModelParamMap.insert(std::make_pair("Cqq3_a33ar_LNP", std::cref(Cqq3_a33ar_LNP)));
    ModelParamMap.insert(std::make_pair("Cqq3_aa33r_LNP", std::cref(Cqq3_aa33r_LNP)));
    ModelParamMap.insert(std::make_pair("Cqq3_3333r_LNP", std::cref(Cqq3_3333r_LNP)));
    ModelParamMap.insert(std::make_pair("Cuu_1111r_LNP", std::cref(Cuu_1111r_LNP)));
    ModelParamMap.insert(std::make_pair("Cuu_2222r_LNP", std::cref(Cuu_2222r_LNP)));
    ModelParamMap.insert(std::make_pair("Cuu_3333r_LNP", std::cref(Cuu_3333r_LNP)));
    ModelParamMap.insert(std::make_pair("Cuu_1122r_LNP", std::cref(Cuu_1122r_LNP)));
    ModelParamMap.insert(std::make_pair("Cuu_1133r_LNP", std::cref(Cuu_1133r_LNP)));
    ModelParamMap.insert(std::make_pair("Cuu_2233r_LNP", std::cref(Cuu_2233r_LNP)));
    ModelParamMap.insert(std::make_pair("Cuu_1221r_LNP", std::cref(Cuu_1221r_LNP)));
    ModelParamMap.insert(std::make_pair("Cuu_1331r_LNP", std::cref(Cuu_1331r_LNP)));
    ModelParamMap.insert(std::make_pair("Cuu_2332r_LNP", std::cref(Cuu_2332r_LNP)));
    ModelParamMap.insert(std::make_pair("Cdd_1111r_LNP", std::cref(Cdd_1111r_LNP)));
    ModelParamMap.insert(std::make_pair("Cdd_2222r_LNP", std::cref(Cdd_2222r_LNP)));
    ModelParamMap.insert(std::make_pair("Cdd_3333r_LNP", std::cref(Cdd_3333r_LNP)));
    ModelParamMap.insert(std::make_pair("Cdd_1122r_LNP", std::cref(Cdd_1122r_LNP)));
    ModelParamMap.insert(std::make_pair("Cdd_1133r_LNP", std::cref(Cdd_1133r_LNP)));
    ModelParamMap.insert(std::make_pair("Cdd_2233r_LNP", std::cref(Cdd_2233r_LNP)));
    ModelParamMap.insert(std::make_pair("Cdd_1221r_LNP", std::cref(Cdd_1221r_LNP)));
    ModelParamMap.insert(std::make_pair("Cdd_1331r_LNP", std::cref(Cdd_1331r_LNP)));
    ModelParamMap.insert(std::make_pair("Cdd_2332r_LNP", std::cref(Cdd_2332r_LNP)));
    ModelParamMap.insert(std::make_pair("Cud1_1111r_LNP", std::cref(Cud1_1111r_LNP)));
    ModelParamMap.insert(std::make_pair("Cud1_2222r_LNP", std::cref(Cud1_2222r_LNP)));
    ModelParamMap.insert(std::make_pair("Cud1_3333r_LNP", std::cref(Cud1_3333r_LNP)));
    ModelParamMap.insert(std::make_pair("Cud1_1122r_LNP", std::cref(Cud1_1122r_LNP)));
    ModelParamMap.insert(std::make_pair("Cud1_1133r_LNP", std::cref(Cud1_1133r_LNP)));
    ModelParamMap.insert(std::make_pair("Cud1_2233r_LNP", std::cref(Cud1_2233r_LNP)));
    ModelParamMap.insert(std::make_pair("Cud1_2211r_LNP", std::cref(Cud1_2211r_LNP)));
    ModelParamMap.insert(std::make_pair("Cud1_3311r_LNP", std::cref(Cud1_3311r_LNP)));
    ModelParamMap.insert(std::make_pair("Cud1_3322r_LNP", std::cref(Cud1_3322r_LNP)));
    ModelParamMap.insert(std::make_pair("Cud8_1111r_LNP", std::cref(Cud8_1111r_LNP)));
    ModelParamMap.insert(std::make_pair("Cud8_2222r_LNP", std::cref(Cud8_2222r_LNP)));
    ModelParamMap.insert(std::make_pair("Cud8_3333r_LNP", std::cref(Cud8_3333r_LNP)));
    ModelParamMap.insert(std::make_pair("Cud8_1122r_LNP", std::cref(Cud8_1122r_LNP)));
    ModelParamMap.insert(std::make_pair("Cud8_1133r_LNP", std::cref(Cud8_1133r_LNP)));
    ModelParamMap.insert(std::make_pair("Cud8_2233r_LNP", std::cref(Cud8_2233r_LNP)));
    ModelParamMap.insert(std::make_pair("Cud8_2211r_LNP", std::cref(Cud8_2211r_LNP)));
    ModelParamMap.insert(std::make_pair("Cud8_3311r_LNP", std::cref(Cud8_3311r_LNP)));
    ModelParamMap.insert(std::make_pair("Cud8_3322r_LNP", std::cref(Cud8_3322r_LNP)));
    ModelParamMap.insert(std::make_pair("Cqu1_aa11r_LNP", std::cref(Cqu1_aa11r_LNP)));
    ModelParamMap.insert(std::make_pair("Cqu1_aa22r_LNP", std::cref(Cqu1_aa22r_LNP)));
    ModelParamMap.insert(std::make_pair("Cqu1_aa33r_LNP", std::cref(Cqu1_aa33r_LNP)));
    ModelParamMap.insert(std::make_pair("Cqu1_3311r_LNP", std::cref(Cqu1_3311r_LNP)));
    ModelParamMap.insert(std::make_pair("Cqu1_3322r_LNP", std::cref(Cqu1_3322r_LNP)));
    ModelParamMap.insert(std::make_pair("Cqu1_3333r_LNP", std::cref(Cqu1_3333r_LNP)));
    ModelParamMap.insert(std::make_pair("Cqu8_aa11r_LNP", std::cref(Cqu8_aa11r_LNP)));
    ModelParamMap.insert(std::make_pair("Cqu8_aa22r_LNP", std::cref(Cqu8_aa22r_LNP)));
    ModelParamMap.insert(std::make_pair("Cqu8_aa33r_LNP", std::cref(Cqu8_aa33r_LNP)));
    ModelParamMap.insert(std::make_pair("Cqu8_3311r_LNP", std::cref(Cqu8_3311r_LNP)));
    ModelParamMap.insert(std::make_pair("Cqu8_3322r_LNP", std::cref(Cqu8_3322r_LNP)));
    ModelParamMap.insert(std::make_pair("Cqu8_3333r_LNP", std::cref(Cqu8_3333r_LNP)));
    ModelParamMap.insert(std::make_pair("Cqd1_aa11r_LNP", std::cref(Cqd1_aa11r_LNP)));
    ModelParamMap.insert(std::make_pair("Cqd1_aa22r_LNP", std::cref(Cqd1_aa22r_LNP))); 
    ModelParamMap.insert(std::make_pair("Cqd1_aa33r_LNP", std::cref(Cqd1_aa33r_LNP)));
    ModelParamMap.insert(std::make_pair("Cqd1_3311r_LNP", std::cref(Cqd1_3311r_LNP)));
    ModelParamMap.insert(std::make_pair("Cqd1_3322r_LNP", std::cref(Cqd1_3322r_LNP)));
    ModelParamMap.insert(std::make_pair("Cqd1_3333r_LNP", std::cref(Cqd1_3333r_LNP)));
    ModelParamMap.insert(std::make_pair("Cqd8_aa11r_LNP", std::cref(Cqd8_aa11r_LNP)));
    ModelParamMap.insert(std::make_pair("Cqd8_aa22r_LNP", std::cref(Cqd8_aa22r_LNP)));
    ModelParamMap.insert(std::make_pair("Cqd8_aa33r_LNP", std::cref(Cqd8_aa33r_LNP)));
    ModelParamMap.insert(std::make_pair("Cqd8_3311r_LNP", std::cref(Cqd8_3311r_LNP)));
    ModelParamMap.insert(std::make_pair("Cqd8_3322r_LNP", std::cref(Cqd8_3322r_LNP)));
    ModelParamMap.insert(std::make_pair("Cqd8_3333r_LNP", std::cref(Cqd8_3333r_LNP)));
    ModelParamMap.insert(std::make_pair("Lambda_NP", std::cref(Lambda_NP)));
}

bool NPSMEFTd6MFP::Init(const std::map<std::string, double>& DPars) {

    if (SMEFTBasisFlag != "UP") {
        std::cerr << "Error in NPSMEFTd6MFP::Init() : SMEFT basis flag must be 'UP' for this model." << std::endl;
        return false;
    }

    return (NPSMEFTd6General::Init(DPars));
}

void NPSMEFTd6MFP::setParameter(const std::string name, const double& value)
{
    
 if(name.compare("CG_LNP") == 0 ) { 
	
	CG_LNP = value;
	
} else if(name.compare("CW_LNP") == 0 ) { 
	
	CW_LNP = value;
	
} else if(name.compare("CHG_LNP") == 0 ) { 
	
	CHG_LNP = value;
	
} else if(name.compare("CHW_LNP") == 0 ) { 
	
	CHW_LNP = value;
	
} else if(name.compare("CHB_LNP") == 0 ) { 
	
	CHB_LNP = value;
	
} else if(name.compare("CHWB_LNP") == 0 ) { 
	
	CHWB_LNP = value;
	
} else if(name.compare("CHD_LNP") == 0 ) { 
	
	CHD_LNP = value;
	
} else if(name.compare("CHbox_LNP") == 0 ) { 
	
	CHbox_LNP = value;
	
} else if(name.compare("CH_LNP") == 0 ) { 
	
	CH_LNP = value;
	
} else if(name.compare("CHl1_11r_LNP") == 0 ) { 
	
	CHl1_11r_LNP = value;
	
} else if(name.compare("CHl1_22r_LNP") == 0 ) { 
	
	CHl1_22r_LNP = value;
	
} else if(name.compare("CHl1_33r_LNP") == 0 ) { 
	
	CHl1_33r_LNP = value;
	
} else if(name.compare("CHl3_11r_LNP") == 0 ) { 
	
	CHl3_11r_LNP = value;
	
} else if(name.compare("CHl3_22r_LNP") == 0 ) { 
	
	CHl3_22r_LNP = value;
	
} else if(name.compare("CHl3_33r_LNP") == 0 ) { 
	
	CHl3_33r_LNP = value;
	
} else if(name.compare("CHe_11r_LNP") == 0 ) { 
	
	CHe_11r_LNP = value;
	
} else if(name.compare("CHe_22r_LNP") == 0 ) { 
	
	CHe_22r_LNP = value;
	
} else if(name.compare("CHe_33r_LNP") == 0 ) { 
	
	CHe_33r_LNP = value;
	
} else if(name.compare("CHq1_aar_LNP") == 0 ) { 
	
	CHq1_aar_LNP = value;
	
} else if(name.compare("CHq1_33r_LNP") == 0 ) { 
	
	CHq1_33r_LNP = value;
	
} else if(name.compare("CHq3_aar_LNP") == 0 ) { 
	
	CHq3_aar_LNP = value;
	
} else if(name.compare("CHq3_33r_LNP") == 0 ) { 
	
	CHq3_33r_LNP = value;
	
} else if(name.compare("CHu_11r_LNP") == 0 ) { 
	
	CHu_11r_LNP = value;
	
} else if(name.compare("CHu_22r_LNP") == 0 ) { 
	
	CHu_22r_LNP = value;
	
} else if(name.compare("CHu_33r_LNP") == 0 ) { 
	
	CHu_33r_LNP = value;
	
} else if(name.compare("CHd_11r_LNP") == 0 ) { 
	
	CHd_11r_LNP = value;
	
} else if(name.compare("CHd_22r_LNP") == 0 ) { 
	
	CHd_22r_LNP = value;
	
} else if(name.compare("CHd_33r_LNP") == 0 ) { 
	
	CHd_33r_LNP = value;
	
} else if(name.compare("CHud_11r_LNP") == 0 ) { 
	
	CHud_11r_LNP = value;
	
} else if(name.compare("CuH_33r_LNP") == 0 ) { 
	
	CuH_33r_LNP = value;
	
} else if(name.compare("CuG_33r_LNP") == 0 ) { 
	
	CuG_33r_LNP = value;
	
} else if(name.compare("CuW_33r_LNP") == 0 ) { 
	
	CuW_33r_LNP = value;
	
} else if(name.compare("CuB_33r_LNP") == 0 ) { 
	
	CuB_33r_LNP = value;
	
} else if(name.compare("Cll_1111r_LNP") == 0 ) { 
	
	Cll_1111r_LNP = value;
	
} else if(name.compare("Cll_2222r_LNP") == 0 ) { 
	
	Cll_2222r_LNP = value;
	
} else if(name.compare("Cll_3333r_LNP") == 0 ) { 
	
	Cll_3333r_LNP = value;
	
} else if(name.compare("Cll_1122r_LNP") == 0 ) { 
	
	Cll_1122r_LNP = value;
	
} else if(name.compare("Cll_1133r_LNP") == 0 ) { 
	
	Cll_1133r_LNP = value;
	
} else if(name.compare("Cll_2233r_LNP") == 0 ) { 
	
	Cll_2233r_LNP = value;
	
} else if(name.compare("Cll_1221r_LNP") == 0 ) { 
	
	Cll_1221r_LNP = value;
	
} else if(name.compare("Cll_1331r_LNP") == 0 ) { 
	
	Cll_1331r_LNP = value;
	
} else if(name.compare("Cll_2332r_LNP") == 0 ) { 
	
	Cll_2332r_LNP = value;
	
} else if(name.compare("Cll_1213r_LNP") == 0 ) { 
	
	Cll_1213r_LNP = value;
	
} else if(name.compare("Clq1_11aar_LNP") == 0 ) { 
	
	Clq1_11aar_LNP = value;
	
} else if(name.compare("Clq1_22aar_LNP") == 0 ) { 
	
	Clq1_22aar_LNP = value;
	
} else if(name.compare("Clq1_33aar_LNP") == 0 ) { 
	
	Clq1_33aar_LNP = value;
	
} else if(name.compare("Clq1_1133r_LNP") == 0 ) { 
	
	Clq1_1133r_LNP = value;
	
} else if(name.compare("Clq1_2233r_LNP") == 0 ) { 
	
	Clq1_2233r_LNP = value;
	
} else if(name.compare("Clq1_3333r_LNP") == 0 ) { 
	
	Clq1_3333r_LNP = value;
	
} else if(name.compare("Clq3_11aar_LNP") == 0 ) { 
	
	Clq3_11aar_LNP = value;
	
} else if(name.compare("Clq3_22aar_LNP") == 0 ) { 
	
	Clq3_22aar_LNP = value;
	
} else if(name.compare("Clq3_33aar_LNP") == 0 ) { 
	
	Clq3_33aar_LNP = value;
	
} else if(name.compare("Clq3_1133r_LNP") == 0 ) { 
	
	Clq3_1133r_LNP = value;
	
} else if(name.compare("Clq3_2233r_LNP") == 0 ) { 
	
	Clq3_2233r_LNP = value;
	
} else if(name.compare("Clq3_3333r_LNP") == 0 ) { 
	
	Clq3_3333r_LNP = value;
	
} else if(name.compare("Cee_1111r_LNP") == 0 ) { 
	
	Cee_1111r_LNP = value;
	
} else if(name.compare("Cee_2222r_LNP") == 0 ) { 
	
	Cee_2222r_LNP = value;
	
} else if(name.compare("Cee_3333r_LNP") == 0 ) { 
	
	Cee_3333r_LNP = value;
	
} else if(name.compare("Cee_1122r_LNP") == 0 ) { 
	
	Cee_1122r_LNP = value;
	
} else if(name.compare("Cee_1133r_LNP") == 0 ) { 
	
	Cee_1133r_LNP = value;
	
} else if(name.compare("Cee_2233r_LNP") == 0 ) { 
	
	Cee_2233r_LNP = value;
	
} else if(name.compare("Ceu_1111r_LNP") == 0 ) { 
	
	Ceu_1111r_LNP = value;
	
} else if(name.compare("Ceu_2222r_LNP") == 0 ) { 
	
	Ceu_2222r_LNP = value;
	
} else if(name.compare("Ceu_3333r_LNP") == 0 ) { 
	
	Ceu_3333r_LNP = value;
	
} else if(name.compare("Ceu_1122r_LNP") == 0 ) { 
	
	Ceu_1122r_LNP = value;
	
} else if(name.compare("Ceu_1133r_LNP") == 0 ) { 
	
	Ceu_1133r_LNP = value;
	
} else if(name.compare("Ceu_2233r_LNP") == 0 ) { 
	
	Ceu_2233r_LNP = value;
	
} else if(name.compare("Ceu_2211r_LNP") == 0 ) { 
	
	Ceu_2211r_LNP = value;
	
} else if(name.compare("Ceu_3311r_LNP") == 0 ) { 
	
	Ceu_3311r_LNP = value;
	
} else if(name.compare("Ceu_3322r_LNP") == 0 ) { 
	
	Ceu_3322r_LNP = value;
	
} else if(name.compare("Ceu_1332r_LNP") == 0 ) { 
	
	Ceu_1332r_LNP = value;
	
} else if(name.compare("Ced_1111r_LNP") == 0 ) { 
	
	Ced_1111r_LNP = value;
	
} else if(name.compare("Ced_2222r_LNP") == 0 ) { 
	
	Ced_2222r_LNP = value;
	
} else if(name.compare("Ced_3333r_LNP") == 0 ) { 
	
	Ced_3333r_LNP = value;
	
} else if(name.compare("Ced_1122r_LNP") == 0 ) { 
	
	Ced_1122r_LNP = value;
	
} else if(name.compare("Ced_1133r_LNP") == 0 ) { 
	
	Ced_1133r_LNP = value;
	
} else if(name.compare("Ced_2233r_LNP") == 0 ) { 
	
	Ced_2233r_LNP = value;
	
} else if(name.compare("Ced_2211r_LNP") == 0 ) { 
	
	Ced_2211r_LNP = value;
	
} else if(name.compare("Ced_3311r_LNP") == 0 ) { 
	
	Ced_3311r_LNP = value;
	
} else if(name.compare("Ced_3322r_LNP") == 0 ) { 
	
	Ced_3322r_LNP = value;
	
} else if(name.compare("Cle_1111r_LNP") == 0 ) { 
	
	Cle_1111r_LNP = value;
	
} else if(name.compare("Cle_2222r_LNP") == 0 ) { 
	
	Cle_2222r_LNP = value;
	
} else if(name.compare("Cle_3333r_LNP") == 0 ) { 
	
	Cle_3333r_LNP = value;
	
} else if(name.compare("Cle_1122r_LNP") == 0 ) { 
	
	Cle_1122r_LNP = value;
	
} else if(name.compare("Cle_1133r_LNP") == 0 ) { 
	
	Cle_1133r_LNP = value;
	
} else if(name.compare("Cle_2233r_LNP") == 0 ) { 
	
	Cle_2233r_LNP = value;
	
} else if(name.compare("Cle_2211r_LNP") == 0 ) { 
	
	Cle_2211r_LNP = value;
	
} else if(name.compare("Cle_3311r_LNP") == 0 ) { 
	
	Cle_3311r_LNP = value;
	
} else if(name.compare("Cle_3322r_LNP") == 0 ) { 
	
	Cle_3322r_LNP = value;
	
} else if(name.compare("Clu_1111r_LNP") == 0 ) { 
	
	Clu_1111r_LNP = value;
	
} else if(name.compare("Clu_2222r_LNP") == 0 ) { 
	
	Clu_2222r_LNP = value;
	
} else if(name.compare("Clu_3333r_LNP") == 0 ) { 
	
	Clu_3333r_LNP = value;
	
} else if(name.compare("Clu_1122r_LNP") == 0 ) { 
	
	Clu_1122r_LNP = value;
	
} else if(name.compare("Clu_1133r_LNP") == 0 ) { 
	
	Clu_1133r_LNP = value;
	
} else if(name.compare("Clu_2233r_LNP") == 0 ) { 
	
	Clu_2233r_LNP = value;
	
} else if(name.compare("Clu_2211r_LNP") == 0 ) { 
	
	Clu_2211r_LNP = value;
	
} else if(name.compare("Clu_3311r_LNP") == 0 ) { 
	
	Clu_3311r_LNP = value;
	
} else if(name.compare("Clu_3322r_LNP") == 0 ) { 
	
	Clu_3322r_LNP = value;
	
} else if(name.compare("Cld_1111r_LNP") == 0 ) { 
	
	Cld_1111r_LNP = value;
	
} else if(name.compare("Cld_2222r_LNP") == 0 ) { 
	
	Cld_2222r_LNP = value;
	
} else if(name.compare("Cld_3333r_LNP") == 0 ) { 
	
	Cld_3333r_LNP = value;
	
} else if(name.compare("Cld_1122r_LNP") == 0 ) { 
	
	Cld_1122r_LNP = value;
	
} else if(name.compare("Cld_1133r_LNP") == 0 ) { 
	
	Cld_1133r_LNP = value;
	
} else if(name.compare("Cld_2233r_LNP") == 0 ) { 
	
	Cld_2233r_LNP = value;
	
} else if(name.compare("Cld_2211r_LNP") == 0 ) { 
	
	Cld_2211r_LNP = value;
	
} else if(name.compare("Cld_3311r_LNP") == 0 ) { 
	
	Cld_3311r_LNP = value;
	
} else if(name.compare("Cld_3322r_LNP") == 0 ) { 
	
	Cld_3322r_LNP = value;
	
} else if(name.compare("Cld_2332r_LNP") == 0 ) { 
	
	Cld_2332r_LNP = value;
	
} else if(name.compare("Cqe_aa11r_LNP") == 0 ) { 
	
	Cqe_aa11r_LNP = value;
	
} else if(name.compare("Cqe_aa22r_LNP") == 0 ) { 
	
	Cqe_aa22r_LNP = value;
	
} else if(name.compare("Cqe_aa33r_LNP") == 0 ) { 
	
	Cqe_aa33r_LNP = value;
	
} else if(name.compare("Cqe_3311r_LNP") == 0 ) { 
	
	Cqe_3311r_LNP = value;
	
} else if(name.compare("Cqe_3322r_LNP") == 0 ) { 
	
	Cqe_3322r_LNP = value;
	
} else if(name.compare("Cqe_3333r_LNP") == 0 ) { 
	
	Cqe_3333r_LNP = value;
	
} else if(name.compare("Cledq_3333r_LNP") == 0 ) { 
	
	Cledq_3333r_LNP = value;
	
} else if(name.compare("Cqq1_aabbr_LNP") == 0 ) { 
	
	Cqq1_aabbr_LNP = value;
	
} else if(name.compare("Cqq1_abbar_LNP") == 0 ) { 
	
	Cqq1_abbar_LNP = value;
	
} else if(name.compare("Cqq1_a33ar_LNP") == 0 ) { 
	
	Cqq1_a33ar_LNP = value;
	
} else if(name.compare("Cqq1_aa33r_LNP") == 0 ) { 
	
	Cqq1_aa33r_LNP = value;
	
} else if(name.compare("Cqq1_3333r_LNP") == 0 ) { 
	
	Cqq1_3333r_LNP = value;
	
} else if(name.compare("Cqq3_aabbr_LNP") == 0 ) { 
	
	Cqq3_aabbr_LNP = value;
	
} else if(name.compare("Cqq3_abbar_LNP") == 0 ) { 
	
	Cqq3_abbar_LNP = value;
	
} else if(name.compare("Cqq3_a33ar_LNP") == 0 ) { 
	
	Cqq3_a33ar_LNP = value;
	
} else if(name.compare("Cqq3_aa33r_LNP") == 0 ) { 
	
	Cqq3_aa33r_LNP = value;
	
} else if(name.compare("Cqq3_3333r_LNP") == 0 ) { 
	
	Cqq3_3333r_LNP = value;
	
} else if(name.compare("Cuu_1111r_LNP") == 0 ) { 
	
	Cuu_1111r_LNP = value;
	
} else if(name.compare("Cuu_2222r_LNP") == 0 ) { 
	
	Cuu_2222r_LNP = value;
	
} else if(name.compare("Cuu_3333r_LNP") == 0 ) { 
	
	Cuu_3333r_LNP = value;
	
} else if(name.compare("Cuu_1122r_LNP") == 0 ) { 
	
	Cuu_1122r_LNP = value;
	
} else if(name.compare("Cuu_1133r_LNP") == 0 ) { 
	
	Cuu_1133r_LNP = value;
	
} else if(name.compare("Cuu_2233r_LNP") == 0 ) { 
	
	Cuu_2233r_LNP = value;
	
} else if(name.compare("Cuu_1221r_LNP") == 0 ) { 
	
	Cuu_1221r_LNP = value;
	
} else if(name.compare("Cuu_1331r_LNP") == 0 ) { 
	
	Cuu_1331r_LNP = value;
	
} else if(name.compare("Cuu_2332r_LNP") == 0 ) { 
	
	Cuu_2332r_LNP = value;
	
} else if(name.compare("Cdd_1111r_LNP") == 0 ) { 
	
	Cdd_1111r_LNP = value;
	
} else if(name.compare("Cdd_2222r_LNP") == 0 ) { 
	
	Cdd_2222r_LNP = value;
	
} else if(name.compare("Cdd_3333r_LNP") == 0 ) { 
	
	Cdd_3333r_LNP = value;
	
} else if(name.compare("Cdd_1122r_LNP") == 0 ) { 
	
	Cdd_1122r_LNP = value;
	
} else if(name.compare("Cdd_1133r_LNP") == 0 ) { 
	
	Cdd_1133r_LNP = value;
	
} else if(name.compare("Cdd_2233r_LNP") == 0 ) { 
	
	Cdd_2233r_LNP = value;
	
} else if(name.compare("Cdd_1221r_LNP") == 0 ) { 
	
	Cdd_1221r_LNP = value;
	
} else if(name.compare("Cdd_1331r_LNP") == 0 ) { 
	
	Cdd_1331r_LNP = value;
	
} else if(name.compare("Cdd_2332r_LNP") == 0 ) { 
	
	Cdd_2332r_LNP = value;
	
} else if(name.compare("Cud1_1111r_LNP") == 0 ) { 
	
	Cud1_1111r_LNP = value;
	
} else if(name.compare("Cud1_2222r_LNP") == 0 ) { 
	
	Cud1_2222r_LNP = value;
	
} else if(name.compare("Cud1_3333r_LNP") == 0 ) { 
	
	Cud1_3333r_LNP = value;
	
} else if(name.compare("Cud1_1122r_LNP") == 0 ) { 
	
	Cud1_1122r_LNP = value;
	
} else if(name.compare("Cud1_2211r_LNP") == 0 ) { 
	
	Cud1_2211r_LNP = value;
	
} else if(name.compare("Cud1_1133r_LNP") == 0 ) { 
	
	Cud1_1133r_LNP = value;
	
} else if(name.compare("Cud1_3311r_LNP") == 0 ) { 
	
	Cud1_3311r_LNP = value;
	
} else if(name.compare("Cud1_2233r_LNP") == 0 ) { 
	
	Cud1_2233r_LNP = value;
	
} else if(name.compare("Cud1_3322r_LNP") == 0 ) { 
	
	Cud1_3322r_LNP = value;
	
} else if(name.compare("Cud8_1111r_LNP") == 0 ) { 
	
	Cud8_1111r_LNP = value;
	
} else if(name.compare("Cud8_2222r_LNP") == 0 ) { 
	
	Cud8_2222r_LNP = value;
	
} else if(name.compare("Cud8_3333r_LNP") == 0 ) { 
	
	Cud8_3333r_LNP = value;
	
} else if(name.compare("Cud8_1122r_LNP") == 0 ) { 
	
	Cud8_1122r_LNP = value;
	
} else if(name.compare("Cud8_2211r_LNP") == 0 ) { 
	
	Cud8_2211r_LNP = value;
	
} else if(name.compare("Cud8_1133r_LNP") == 0 ) { 
	
	Cud8_1133r_LNP = value;
	
} else if(name.compare("Cud8_3311r_LNP") == 0 ) { 
	
	Cud8_3311r_LNP = value;
	
} else if(name.compare("Cud8_2233r_LNP") == 0 ) { 
	
	Cud8_2233r_LNP = value;
	
} else if(name.compare("Cud8_3322r_LNP") == 0 ) { 
	
	Cud8_3322r_LNP = value;
	
} else if(name.compare("Cqu1_aa11r_LNP") == 0 ) { 
	
	Cqu1_aa11r_LNP = value;
	
} else if(name.compare("Cqu1_aa22r_LNP") == 0 ) { 
	
	Cqu1_aa22r_LNP = value;
	
} else if(name.compare("Cqu1_aa33r_LNP") == 0 ) { 
	
	Cqu1_aa33r_LNP = value;
	
} else if(name.compare("Cqu1_3311r_LNP") == 0 ) { 
	
	Cqu1_3311r_LNP = value;
	
} else if(name.compare("Cqu1_3322r_LNP") == 0 ) { 
	
	Cqu1_3322r_LNP = value;
	
} else if(name.compare("Cqu1_3333r_LNP") == 0 ) { 
	
	Cqu1_3333r_LNP = value;
	
} else if(name.compare("Cqu8_aa11r_LNP") == 0 ) { 
	
	Cqu8_aa11r_LNP = value;
	
} else if(name.compare("Cqu8_aa22r_LNP") == 0 ) { 
	
	Cqu8_aa22r_LNP = value;
	
} else if(name.compare("Cqu8_aa33r_LNP") == 0 ) { 
	
	Cqu8_aa33r_LNP = value;
	
} else if(name.compare("Cqu8_3311r_LNP") == 0 ) { 
	
	Cqu8_3311r_LNP = value;
	
} else if(name.compare("Cqu8_3322r_LNP") == 0 ) { 
	
	Cqu8_3322r_LNP = value;
	
} else if(name.compare("Cqu8_3333r_LNP") == 0 ) { 
	
	Cqu8_3333r_LNP = value;
	
} else if(name.compare("Cqd1_aa11r_LNP") == 0 ) { 
	
	Cqd1_aa11r_LNP = value;
	
} else if(name.compare("Cqd1_aa22r_LNP") == 0 ) { 
	
	Cqd1_aa22r_LNP = value;
	
} else if(name.compare("Cqd1_aa33r_LNP") == 0 ) { 
	
	Cqd1_aa33r_LNP = value;
	
} else if(name.compare("Cqd1_3311r_LNP") == 0 ) { 
	
	Cqd1_3311r_LNP = value;
	
} else if(name.compare("Cqd1_3322r_LNP") == 0 ) { 
	
	Cqd1_3322r_LNP = value;
	
} else if(name.compare("Cqd1_3333r_LNP") == 0 ) { 
	
	Cqd1_3333r_LNP = value;
	
} else if(name.compare("Cqd8_aa11r_LNP") == 0 ) { 
	
	Cqd8_aa11r_LNP = value;
	
} else if(name.compare("Cqd8_aa22r_LNP") == 0 ) { 
	
	Cqd8_aa22r_LNP = value;
	
} else if(name.compare("Cqd8_aa33r_LNP") == 0 ) { 
	
	Cqd8_aa33r_LNP = value;
	
} else if(name.compare("Cqd8_3311r_LNP") == 0 ) { 
	
	Cqd8_3311r_LNP = value;
	
} else if(name.compare("Cqd8_3322r_LNP") == 0 ) { 
	
	Cqd8_3322r_LNP = value;
	
} else if(name.compare("Cqd8_3333r_LNP") == 0 ) { 
	
	Cqd8_3333r_LNP = value;
	
} else if (name.compare("Lambda_NP") == 0 )
        Lambda_NP = value;
    else
        NPSMEFTd6General::setParameter(name, value);

}


void NPSMEFTd6MFP::setNPSMEFTd6GeneralParameters()
{
    // Map MFP aggregate parameters to individual-index WC in NPSMEFTd6General.
    // Only parameters that are computed from aggregates are set here;
    // direct inputs (fully-indexed WC in NPSMEFTd6MFPVars) are set by setParameter.

    // CHq1, CHq3: light-quark SU(2) aggregate → gen 1 and gen 2
    CHq1_11r_LNP = CHq1_aar_LNP;
    CHq1_22r_LNP = CHq1_aar_LNP;

    CHq3_11r_LNP = CHq3_aar_LNP;
    CHq3_22r_LNP = CHq3_aar_LNP;

    // Clq1: l independent per generation; q has SU(2) → gen 1 = gen 2
    Clq1_1111r_LNP = Clq1_11aar_LNP;
    Clq1_1122r_LNP = Clq1_11aar_LNP;

    Clq1_2211r_LNP = Clq1_22aar_LNP;
    Clq1_2222r_LNP = Clq1_22aar_LNP;

    Clq1_3311r_LNP = Clq1_33aar_LNP;
    Clq1_3322r_LNP = Clq1_33aar_LNP;

    // Clq3: same structure as Clq1
    Clq3_1111r_LNP = Clq3_11aar_LNP;
    Clq3_1122r_LNP = Clq3_11aar_LNP;

    Clq3_2211r_LNP = Clq3_22aar_LNP;
    Clq3_2222r_LNP = Clq3_22aar_LNP;

    Clq3_3311r_LNP = Clq3_33aar_LNP;
    Clq3_3322r_LNP = Clq3_33aar_LNP;

    // Cqe: q has SU(2) aggregate; e independent per generation
    Cqe_1111r_LNP = Cqe_aa11r_LNP;
    Cqe_2211r_LNP = Cqe_aa11r_LNP;

    Cqe_1122r_LNP = Cqe_aa22r_LNP;
    Cqe_2222r_LNP = Cqe_aa22r_LNP;

    Cqe_1133r_LNP = Cqe_aa33r_LNP;
    Cqe_2233r_LNP = Cqe_aa33r_LNP;

    // Cqq1: light-light and light-top sectors
    Cqq1_1111r_LNP = Cqq1_aabbr_LNP + Cqq1_abbar_LNP;
    Cqq1_2222r_LNP = Cqq1_aabbr_LNP + Cqq1_abbar_LNP;

    Cqq1_1122r_LNP = Cqq1_aabbr_LNP;

    Cqq1_1221r_LNP = Cqq1_abbar_LNP;

    Cqq1_1331r_LNP = Cqq1_a33ar_LNP;
    Cqq1_2332r_LNP = Cqq1_a33ar_LNP;

    Cqq1_1133r_LNP = Cqq1_aa33r_LNP;
    Cqq1_2233r_LNP = Cqq1_aa33r_LNP;

    // Cqq3: same structure as Cqq1
    Cqq3_1111r_LNP = Cqq3_aabbr_LNP + Cqq3_abbar_LNP;
    Cqq3_2222r_LNP = Cqq3_aabbr_LNP + Cqq3_abbar_LNP;

    Cqq3_1122r_LNP = Cqq3_aabbr_LNP;

    Cqq3_1221r_LNP = Cqq3_abbar_LNP;

    Cqq3_1331r_LNP = Cqq3_a33ar_LNP;
    Cqq3_2332r_LNP = Cqq3_a33ar_LNP;

    Cqq3_1133r_LNP = Cqq3_aa33r_LNP;
    Cqq3_2233r_LNP = Cqq3_aa33r_LNP;

    // Cqu1: q=light U(2) aggregate; u independent per generation
    Cqu1_1111r_LNP = Cqu1_aa11r_LNP;
    Cqu1_2211r_LNP = Cqu1_aa11r_LNP;

    Cqu1_1122r_LNP = Cqu1_aa22r_LNP;
    Cqu1_2222r_LNP = Cqu1_aa22r_LNP;

    Cqu1_1133r_LNP = Cqu1_aa33r_LNP;
    Cqu1_2233r_LNP = Cqu1_aa33r_LNP;

    // Cqu8: same structure as Cqu1
    Cqu8_1111r_LNP = Cqu8_aa11r_LNP;
    Cqu8_2211r_LNP = Cqu8_aa11r_LNP;

    Cqu8_1122r_LNP = Cqu8_aa22r_LNP;
    Cqu8_2222r_LNP = Cqu8_aa22r_LNP;

    Cqu8_1133r_LNP = Cqu8_aa33r_LNP;
    Cqu8_2233r_LNP = Cqu8_aa33r_LNP;

    // Cqd1: q=light U(2) aggregate; d independent per generation
    Cqd1_1111r_LNP = Cqd1_aa11r_LNP;
    Cqd1_2211r_LNP = Cqd1_aa11r_LNP;

    Cqd1_1122r_LNP = Cqd1_aa22r_LNP;
    Cqd1_2222r_LNP = Cqd1_aa22r_LNP;

    Cqd1_1133r_LNP = Cqd1_aa33r_LNP;
    Cqd1_2233r_LNP = Cqd1_aa33r_LNP;

    // Cqd8: same structure as Cqd1
    Cqd8_1111r_LNP = Cqd8_aa11r_LNP;
    Cqd8_2211r_LNP = Cqd8_aa11r_LNP;

    Cqd8_1122r_LNP = Cqd8_aa22r_LNP;
    Cqd8_2222r_LNP = Cqd8_aa22r_LNP;

    Cqd8_1133r_LNP = Cqd8_aa33r_LNP;
    Cqd8_2233r_LNP = Cqd8_aa33r_LNP;
}



bool NPSMEFTd6MFP::PostUpdate()
{
    GenerateSMInitialConditions();
    
    setNPSMEFTd6GeneralParameters();
    
    if (!NPSMEFTd6General::PostUpdate()) return (false);
    
    return (true);
        
}
