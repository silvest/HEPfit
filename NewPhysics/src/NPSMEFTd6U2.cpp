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
    "CG_LNP", "CW_LNP", "CHG_LNP", "CHW_LNP", "CHB_LNP", "CHWB_LNP", "CHD_LNP", "CHbox_LNP", "CH_LNP",
    "CHl1_aar_LNP", "CHl1_33r_LNP", "CHl3_aar_LNP", "CHl3_33r_LNP", "CHe_aar_LNP", "CHe_33r_LNP",
    "CHq1_aar_LNP", "CHq1_33r_LNP", "CHq3_aar_LNP", "CHq3_33r_LNP", "CHu_aar_LNP", "CHu_33r_LNP",
    "CHd_aar_LNP", "CHd_33r_LNP", "CHud_33r_LNP", "CeH_33r_LNP", "CuH_33r_LNP", "CdH_33r_LNP",
    "CuG_33r_LNP", "CuW_33r_LNP", "CuB_33r_LNP", "CdG_33r_LNP", "CdW_33r_LNP", "CdB_33r_LNP",
    "CeW_33r_LNP", "CeB_33r_LNP", "Cll_aabbr_LNP", "Cll_abbar_LNP", "Cll_aa33r_LNP", "Cll_a33ar_LNP", "Cll_3333r_LNP",
    "Clq1_aabbr_LNP", "Clq1_aa33r_LNP", "Clq1_33aar_LNP", "Clq1_3333r_LNP", "Clq3_aabbr_LNP", "Clq3_aa33r_LNP", "Clq3_33aar_LNP", "Clq3_3333r_LNP",
    "Cee_aabbr_LNP", "Cee_aa33r_LNP", "Cee_3333r_LNP", "Ceu_aabbr_LNP", "Ceu_aa33r_LNP", "Ceu_33aar_LNP", "Ceu_3333r_LNP",
    "Ced_aabbr_LNP", "Ced_aa33r_LNP", "Ced_33aar_LNP", "Ced_3333r_LNP", "Cle_aabbr_LNP", "Cle_aa33r_LNP", "Cle_33aar_LNP", "Cle_3333r_LNP",
    "Clu_aabbr_LNP", "Clu_aa33r_LNP", "Clu_33aar_LNP", "Clu_3333r_LNP", "Cld_aabbr_LNP", "Cld_aa33r_LNP", "Cld_33aar_LNP", "Cld_3333r_LNP",
    "Cqe_aabbr_LNP", "Cqe_aa33r_LNP", "Cqe_33aar_LNP", "Cqe_3333r_LNP", "Cledq_3333r_LNP",
    "Cqq1_aabbr_LNP", "Cqq1_abbar_LNP", "Cqq1_aa33r_LNP", "Cqq1_a33ar_LNP", "Cqq1_3333r_LNP",
    "Cqq3_aabbr_LNP", "Cqq3_abbar_LNP", "Cqq3_aa33r_LNP", "Cqq3_a33ar_LNP", "Cqq3_3333r_LNP",
    "Cuu_aabbr_LNP", "Cuu_abbar_LNP", "Cuu_aa33r_LNP", "Cuu_a33ar_LNP", "Cuu_3333r_LNP",
    "Cdd_aabbr_LNP", "Cdd_abbar_LNP", "Cdd_aa33r_LNP", "Cdd_a33ar_LNP", "Cdd_3333r_LNP",
    "Cud1_aabbr_LNP", "Cud1_aa33r_LNP", "Cud1_33aar_LNP", "Cud1_3333r_LNP", "Cud8_aabbr_LNP", "Cud8_aa33r_LNP", "Cud8_33aar_LNP", "Cud8_3333r_LNP",
    "Cqu1_aabbr_LNP", "Cqu1_aa33r_LNP", "Cqu1_33aar_LNP", "Cqu1_3333r_LNP", "Cqu8_aabbr_LNP", "Cqu8_aa33r_LNP", "Cqu8_33aar_LNP", "Cqu8_3333r_LNP",
    "Cqd1_aabbr_LNP", "Cqd1_aa33r_LNP", "Cqd1_33aar_LNP", "Cqd1_3333r_LNP", "Cqd8_aabbr_LNP", "Cqd8_aa33r_LNP", "Cqd8_33aar_LNP", "Cqd8_3333r_LNP",
    "Cquqd1_3333r_LNP", "Cquqd8_3333r_LNP", "Clequ1_3333r_LNP", "Clequ3_3333r_LNP",
    "Lambda_NP"
};

NPSMEFTd6U2::NPSMEFTd6U2()
: NPSMEFTd6General() {
    setModelName("NPSMEFTd6U2");
    ModelParamMap.insert(std::make_pair("CG_LNP", std::cref(CG_LNP)));
    ModelParamMap.insert(std::make_pair("CW_LNP", std::cref(CW_LNP)));
    ModelParamMap.insert(std::make_pair("CHG_LNP", std::cref(CHG_LNP)));
    ModelParamMap.insert(std::make_pair("CHW_LNP", std::cref(CHW_LNP)));
    ModelParamMap.insert(std::make_pair("CHB_LNP", std::cref(CHB_LNP)));
    ModelParamMap.insert(std::make_pair("CHWB_LNP", std::cref(CHWB_LNP)));
    ModelParamMap.insert(std::make_pair("CHD_LNP", std::cref(CHD_LNP)));
    ModelParamMap.insert(std::make_pair("CHbox_LNP", std::cref(CHbox_LNP)));
    ModelParamMap.insert(std::make_pair("CH_LNP", std::cref(CH_LNP)));
    ModelParamMap.insert(std::make_pair("CHl1_aar_LNP", std::cref(CHl1_aar_LNP)));
    ModelParamMap.insert(std::make_pair("CHl1_33r_LNP", std::cref(CHl1_33r_LNP)));
    ModelParamMap.insert(std::make_pair("CHl3_aar_LNP", std::cref(CHl3_aar_LNP)));
    ModelParamMap.insert(std::make_pair("CHl3_33r_LNP", std::cref(CHl3_33r_LNP)));
    ModelParamMap.insert(std::make_pair("CHe_aar_LNP", std::cref(CHe_aar_LNP)));
    ModelParamMap.insert(std::make_pair("CHe_33r_LNP", std::cref(CHe_33r_LNP)));
    ModelParamMap.insert(std::make_pair("CHq1_aar_LNP", std::cref(CHq1_aar_LNP)));
    ModelParamMap.insert(std::make_pair("CHq1_33r_LNP", std::cref(CHq1_33r_LNP)));
    ModelParamMap.insert(std::make_pair("CHq3_aar_LNP", std::cref(CHq3_aar_LNP)));
    ModelParamMap.insert(std::make_pair("CHq3_33r_LNP", std::cref(CHq3_33r_LNP)));
    ModelParamMap.insert(std::make_pair("CHu_aar_LNP", std::cref(CHu_aar_LNP)));
    ModelParamMap.insert(std::make_pair("CHu_33r_LNP", std::cref(CHu_33r_LNP)));
    ModelParamMap.insert(std::make_pair("CHd_aar_LNP", std::cref(CHd_aar_LNP)));
    ModelParamMap.insert(std::make_pair("CHd_33r_LNP", std::cref(CHd_33r_LNP)));
    ModelParamMap.insert(std::make_pair("CHud_33r_LNP", std::cref(CHud_33r_LNP)));
    ModelParamMap.insert(std::make_pair("CeH_33r_LNP", std::cref(CeH_33r_LNP)));
    ModelParamMap.insert(std::make_pair("CuH_33r_LNP", std::cref(CuH_33r_LNP)));
    ModelParamMap.insert(std::make_pair("CdH_33r_LNP", std::cref(CdH_33r_LNP)));
    ModelParamMap.insert(std::make_pair("CuG_33r_LNP", std::cref(CuG_33r_LNP)));
    ModelParamMap.insert(std::make_pair("CuW_33r_LNP", std::cref(CuW_33r_LNP)));
    ModelParamMap.insert(std::make_pair("CuB_33r_LNP", std::cref(CuB_33r_LNP)));
    ModelParamMap.insert(std::make_pair("CdG_33r_LNP", std::cref(CdG_33r_LNP)));
    ModelParamMap.insert(std::make_pair("CdW_33r_LNP", std::cref(CdW_33r_LNP)));
    ModelParamMap.insert(std::make_pair("CdB_33r_LNP", std::cref(CdB_33r_LNP)));
    ModelParamMap.insert(std::make_pair("CeW_33r_LNP", std::cref(CeW_33r_LNP)));
    ModelParamMap.insert(std::make_pair("CeB_33r_LNP", std::cref(CeB_33r_LNP)));
    ModelParamMap.insert(std::make_pair("Cll_aabbr_LNP", std::cref(Cll_aabbr_LNP)));
    ModelParamMap.insert(std::make_pair("Cll_abbar_LNP", std::cref(Cll_abbar_LNP)));
    ModelParamMap.insert(std::make_pair("Cll_aa33r_LNP", std::cref(Cll_aa33r_LNP)));
    ModelParamMap.insert(std::make_pair("Cll_a33ar_LNP", std::cref(Cll_a33ar_LNP)));
    ModelParamMap.insert(std::make_pair("Cll_3333r_LNP", std::cref(Cll_3333r_LNP)));
    ModelParamMap.insert(std::make_pair("Clq1_aabbr_LNP", std::cref(Clq1_aabbr_LNP)));
    ModelParamMap.insert(std::make_pair("Clq1_aa33r_LNP", std::cref(Clq1_aa33r_LNP)));
    ModelParamMap.insert(std::make_pair("Clq1_33aar_LNP", std::cref(Clq1_33aar_LNP)));
    ModelParamMap.insert(std::make_pair("Clq1_3333r_LNP", std::cref(Clq1_3333r_LNP)));
    ModelParamMap.insert(std::make_pair("Clq3_aabbr_LNP", std::cref(Clq3_aabbr_LNP)));
    ModelParamMap.insert(std::make_pair("Clq3_aa33r_LNP", std::cref(Clq3_aa33r_LNP)));
    ModelParamMap.insert(std::make_pair("Clq3_33aar_LNP", std::cref(Clq3_33aar_LNP)));
    ModelParamMap.insert(std::make_pair("Clq3_3333r_LNP", std::cref(Clq3_3333r_LNP)));
    ModelParamMap.insert(std::make_pair("Cee_aabbr_LNP", std::cref(Cee_aabbr_LNP)));
    ModelParamMap.insert(std::make_pair("Cee_aa33r_LNP", std::cref(Cee_aa33r_LNP)));
    ModelParamMap.insert(std::make_pair("Cee_3333r_LNP", std::cref(Cee_3333r_LNP)));
    ModelParamMap.insert(std::make_pair("Ceu_aabbr_LNP", std::cref(Ceu_aabbr_LNP)));
    ModelParamMap.insert(std::make_pair("Ceu_aa33r_LNP", std::cref(Ceu_aa33r_LNP)));
    ModelParamMap.insert(std::make_pair("Ceu_33aar_LNP", std::cref(Ceu_33aar_LNP)));
    ModelParamMap.insert(std::make_pair("Ceu_3333r_LNP", std::cref(Ceu_3333r_LNP)));
    ModelParamMap.insert(std::make_pair("Ced_aabbr_LNP", std::cref(Ced_aabbr_LNP)));
    ModelParamMap.insert(std::make_pair("Ced_aa33r_LNP", std::cref(Ced_aa33r_LNP)));
    ModelParamMap.insert(std::make_pair("Ced_33aar_LNP", std::cref(Ced_33aar_LNP)));
    ModelParamMap.insert(std::make_pair("Ced_3333r_LNP", std::cref(Ced_3333r_LNP)));
    ModelParamMap.insert(std::make_pair("Cle_aabbr_LNP", std::cref(Cle_aabbr_LNP)));
    ModelParamMap.insert(std::make_pair("Cle_aa33r_LNP", std::cref(Cle_aa33r_LNP)));
    ModelParamMap.insert(std::make_pair("Cle_33aar_LNP", std::cref(Cle_33aar_LNP)));
    ModelParamMap.insert(std::make_pair("Cle_3333r_LNP", std::cref(Cle_3333r_LNP)));
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
    ModelParamMap.insert(std::make_pair("Cledq_3333r_LNP", std::cref(Cledq_3333r_LNP)));
    ModelParamMap.insert(std::make_pair("Cqq1_aabbr_LNP", std::cref(Cqq1_aabbr_LNP)));
    ModelParamMap.insert(std::make_pair("Cqq1_abbar_LNP", std::cref(Cqq1_abbar_LNP)));
    ModelParamMap.insert(std::make_pair("Cqq1_aa33r_LNP", std::cref(Cqq1_aa33r_LNP)));
    ModelParamMap.insert(std::make_pair("Cqq1_a33ar_LNP", std::cref(Cqq1_a33ar_LNP)));
    ModelParamMap.insert(std::make_pair("Cqq1_3333r_LNP", std::cref(Cqq1_3333r_LNP)));
    ModelParamMap.insert(std::make_pair("Cqq3_aabbr_LNP", std::cref(Cqq3_aabbr_LNP)));
    ModelParamMap.insert(std::make_pair("Cqq3_abbar_LNP", std::cref(Cqq3_abbar_LNP)));
    ModelParamMap.insert(std::make_pair("Cqq3_aa33r_LNP", std::cref(Cqq3_aa33r_LNP)));
    ModelParamMap.insert(std::make_pair("Cqq3_a33ar_LNP", std::cref(Cqq3_a33ar_LNP)));
    ModelParamMap.insert(std::make_pair("Cqq3_3333r_LNP", std::cref(Cqq3_3333r_LNP)));
    ModelParamMap.insert(std::make_pair("Cuu_aabbr_LNP", std::cref(Cuu_aabbr_LNP)));
    ModelParamMap.insert(std::make_pair("Cuu_abbar_LNP", std::cref(Cuu_abbar_LNP)));
    ModelParamMap.insert(std::make_pair("Cuu_aa33r_LNP", std::cref(Cuu_aa33r_LNP)));
    ModelParamMap.insert(std::make_pair("Cuu_a33ar_LNP", std::cref(Cuu_a33ar_LNP)));
    ModelParamMap.insert(std::make_pair("Cuu_3333r_LNP", std::cref(Cuu_3333r_LNP)));
    ModelParamMap.insert(std::make_pair("Cdd_aabbr_LNP", std::cref(Cdd_aabbr_LNP)));
    ModelParamMap.insert(std::make_pair("Cdd_abbar_LNP", std::cref(Cdd_abbar_LNP)));
    ModelParamMap.insert(std::make_pair("Cdd_aa33r_LNP", std::cref(Cdd_aa33r_LNP)));
    ModelParamMap.insert(std::make_pair("Cdd_a33ar_LNP", std::cref(Cdd_a33ar_LNP)));
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
    ModelParamMap.insert(std::make_pair("Cquqd1_3333r_LNP", std::cref(Cquqd1_3333r_LNP)));
    ModelParamMap.insert(std::make_pair("Cquqd8_3333r_LNP", std::cref(Cquqd8_3333r_LNP)));
    ModelParamMap.insert(std::make_pair("Clequ1_3333r_LNP", std::cref(Clequ1_3333r_LNP)));
    ModelParamMap.insert(std::make_pair("Clequ3_3333r_LNP", std::cref(Clequ3_3333r_LNP)));

}



void NPSMEFTd6U2::setParameter(const std::string name, const double& value)
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
	
} else if(name.compare("CHl1_aar_LNP") == 0 ) { 
	
	CHl1_aar_LNP = value;
	
} else if(name.compare("CHl1_33r_LNP") == 0 ) { 
	
	CHl1_33r_LNP = value;
	
} else if(name.compare("CHl3_aar_LNP") == 0 ) { 
	
	CHl3_aar_LNP = value;
	
} else if(name.compare("CHl3_33r_LNP") == 0 ) { 
	
	CHl3_33r_LNP = value;
	
} else if(name.compare("CHe_aar_LNP") == 0 ) { 
	
	CHe_aar_LNP = value;
	
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
	
} else if(name.compare("CHu_aar_LNP") == 0 ) { 
	
	CHu_aar_LNP = value;
	
} else if(name.compare("CHu_33r_LNP") == 0 ) { 
	
	CHu_33r_LNP = value;
	
} else if(name.compare("CHd_aar_LNP") == 0 ) { 
	
	CHd_aar_LNP = value;
	
} else if(name.compare("CHd_33r_LNP") == 0 ) { 
	
	CHd_33r_LNP = value;
	
} else if(name.compare("CHud_33r_LNP") == 0 ) { 
	
	CHud_33r_LNP = value;
	
} else if(name.compare("CeH_33r_LNP") == 0 ) { 
	
	CeH_33r_LNP = value;
	
} else if(name.compare("CuH_33r_LNP") == 0 ) { 
	
	CuH_33r_LNP = value;
	
} else if(name.compare("CdH_33r_LNP") == 0 ) { 
	
	CdH_33r_LNP = value;
	
} else if(name.compare("CuG_33r_LNP") == 0 ) { 
	
	CuG_33r_LNP = value;
	
} else if(name.compare("CuW_33r_LNP") == 0 ) { 
	
	CuW_33r_LNP = value;
	
} else if(name.compare("CuB_33r_LNP") == 0 ) { 
	
	CuB_33r_LNP = value;
	
} else if(name.compare("CdG_33r_LNP") == 0 ) { 
	
	CdG_33r_LNP = value;
	
} else if(name.compare("CdW_33r_LNP") == 0 ) { 
	
	CdW_33r_LNP = value;
	
} else if(name.compare("CdB_33r_LNP") == 0 ) { 
	
	CdB_33r_LNP = value;
	
} else if(name.compare("CeW_33r_LNP") == 0 ) { 
	
	CeW_33r_LNP = value;
	
} else if(name.compare("CeB_33r_LNP") == 0 ) { 
	
	CeB_33r_LNP = value;
	
} else if(name.compare("Cll_aabbr_LNP") == 0 ) { 
	
	Cll_aabbr_LNP = value;
	
} else if(name.compare("Cll_abbar_LNP") == 0 ) { 
	
	Cll_abbar_LNP = value;
	
} else if(name.compare("Cll_aa33r_LNP") == 0 ) { 
	
	Cll_aa33r_LNP = value;
	
} else if(name.compare("Cll_a33ar_LNP") == 0 ) { 
	
	Cll_a33ar_LNP = value;
	
} else if(name.compare("Cll_3333r_LNP") == 0 ) { 
	
	Cll_3333r_LNP = value;
	
} else if(name.compare("Clq1_aabbr_LNP") == 0 ) { 
	
	Clq1_aabbr_LNP = value;
	
} else if(name.compare("Clq1_aa33r_LNP") == 0 ) { 
	
	Clq1_aa33r_LNP = value;
	
} else if(name.compare("Clq1_33aar_LNP") == 0 ) { 
	
	Clq1_33aar_LNP = value;
	
} else if(name.compare("Clq1_3333r_LNP") == 0 ) { 
	
	Clq1_3333r_LNP = value;
	
} else if(name.compare("Clq3_aabbr_LNP") == 0 ) { 
	
	Clq3_aabbr_LNP = value;
	
} else if(name.compare("Clq3_aa33r_LNP") == 0 ) { 
	
	Clq3_aa33r_LNP = value;
	
} else if(name.compare("Clq3_33aar_LNP") == 0 ) { 
	
	Clq3_33aar_LNP = value;
	
} else if(name.compare("Clq3_3333r_LNP") == 0 ) { 
	
	Clq3_3333r_LNP = value;
	
} else if(name.compare("Cee_aabbr_LNP") == 0 ) { 
	
	Cee_aabbr_LNP = value;
	
} else if(name.compare("Cee_aa33r_LNP") == 0 ) { 
	
	Cee_aa33r_LNP = value;
	
} else if(name.compare("Cee_3333r_LNP") == 0 ) { 
	
	Cee_3333r_LNP = value;
	
} else if(name.compare("Ceu_aabbr_LNP") == 0 ) { 
	
	Ceu_aabbr_LNP = value;
	
} else if(name.compare("Ceu_aa33r_LNP") == 0 ) { 
	
	Ceu_aa33r_LNP = value;
	
} else if(name.compare("Ceu_33aar_LNP") == 0 ) { 
	
	Ceu_33aar_LNP = value;
	
} else if(name.compare("Ceu_3333r_LNP") == 0 ) { 
	
	Ceu_3333r_LNP = value;
	
} else if(name.compare("Ced_aabbr_LNP") == 0 ) { 
	
	Ced_aabbr_LNP = value;
	
} else if(name.compare("Ced_aa33r_LNP") == 0 ) { 
	
	Ced_aa33r_LNP = value;
	
} else if(name.compare("Ced_33aar_LNP") == 0 ) { 
	
	Ced_33aar_LNP = value;
	
} else if(name.compare("Ced_3333r_LNP") == 0 ) { 
	
	Ced_3333r_LNP = value;
	
} else if(name.compare("Cle_aabbr_LNP") == 0 ) { 
	
	Cle_aabbr_LNP = value;
	
} else if(name.compare("Cle_aa33r_LNP") == 0 ) { 
	
	Cle_aa33r_LNP = value;
	
} else if(name.compare("Cle_33aar_LNP") == 0 ) { 
	
	Cle_33aar_LNP = value;
	
} else if(name.compare("Cle_3333r_LNP") == 0 ) { 
	
	Cle_3333r_LNP = value;
	
} else if(name.compare("Clu_aabbr_LNP") == 0 ) { 
	
	Clu_aabbr_LNP = value;
	
} else if(name.compare("Clu_aa33r_LNP") == 0 ) { 
	
	Clu_aa33r_LNP = value;
	
} else if(name.compare("Clu_33aar_LNP") == 0 ) { 
	
	Clu_33aar_LNP = value;
	
} else if(name.compare("Clu_3333r_LNP") == 0 ) { 
	
	Clu_3333r_LNP = value;
	
} else if(name.compare("Cld_aabbr_LNP") == 0 ) { 
	
	Cld_aabbr_LNP = value;
	
} else if(name.compare("Cld_aa33r_LNP") == 0 ) { 
	
	Cld_aa33r_LNP = value;
	
} else if(name.compare("Cld_33aar_LNP") == 0 ) { 
	
	Cld_33aar_LNP = value;
	
} else if(name.compare("Cld_3333r_LNP") == 0 ) { 
	
	Cld_3333r_LNP = value;
	
} else if(name.compare("Cqe_aabbr_LNP") == 0 ) { 
	
	Cqe_aabbr_LNP = value;
	
} else if(name.compare("Cqe_aa33r_LNP") == 0 ) { 
	
	Cqe_aa33r_LNP = value;
	
} else if(name.compare("Cqe_33aar_LNP") == 0 ) { 
	
	Cqe_33aar_LNP = value;
	
} else if(name.compare("Cqe_3333r_LNP") == 0 ) { 
	
	Cqe_3333r_LNP = value;
	
} else if(name.compare("Cledq_3333r_LNP") == 0 ) { 
	
	Cledq_3333r_LNP = value;
	
} else if(name.compare("Cqq1_aabbr_LNP") == 0 ) { 
	
	Cqq1_aabbr_LNP = value;
	
} else if(name.compare("Cqq1_abbar_LNP") == 0 ) { 
	
	Cqq1_abbar_LNP = value;
	
} else if(name.compare("Cqq1_aa33r_LNP") == 0 ) { 
	
	Cqq1_aa33r_LNP = value;
	
} else if(name.compare("Cqq1_a33ar_LNP") == 0 ) { 
	
	Cqq1_a33ar_LNP = value;
	
} else if(name.compare("Cqq1_3333r_LNP") == 0 ) { 
	
	Cqq1_3333r_LNP = value;
	
} else if(name.compare("Cqq3_aabbr_LNP") == 0 ) { 
	
	Cqq3_aabbr_LNP = value;
	
} else if(name.compare("Cqq3_abbar_LNP") == 0 ) { 
	
	Cqq3_abbar_LNP = value;
	
} else if(name.compare("Cqq3_aa33r_LNP") == 0 ) { 
	
	Cqq3_aa33r_LNP = value;
	
} else if(name.compare("Cqq3_a33ar_LNP") == 0 ) { 
	
	Cqq3_a33ar_LNP = value;
	
} else if(name.compare("Cqq3_3333r_LNP") == 0 ) { 
	
	Cqq3_3333r_LNP = value;
	
} else if(name.compare("Cuu_aabbr_LNP") == 0 ) { 
	
	Cuu_aabbr_LNP = value;
	
} else if(name.compare("Cuu_abbar_LNP") == 0 ) { 
	
	Cuu_abbar_LNP = value;
	
} else if(name.compare("Cuu_aa33r_LNP") == 0 ) { 
	
	Cuu_aa33r_LNP = value;
	
} else if(name.compare("Cuu_a33ar_LNP") == 0 ) { 
	
	Cuu_a33ar_LNP = value;
	
} else if(name.compare("Cuu_3333r_LNP") == 0 ) { 
	
	Cuu_3333r_LNP = value;
	
} else if(name.compare("Cdd_aabbr_LNP") == 0 ) { 
	
	Cdd_aabbr_LNP = value;
	
} else if(name.compare("Cdd_abbar_LNP") == 0 ) { 
	
	Cdd_abbar_LNP = value;
	
} else if(name.compare("Cdd_aa33r_LNP") == 0 ) { 
	
	Cdd_aa33r_LNP = value;
	
} else if(name.compare("Cdd_a33ar_LNP") == 0 ) { 
	
	Cdd_a33ar_LNP = value;
	
} else if(name.compare("Cdd_3333r_LNP") == 0 ) { 
	
	Cdd_3333r_LNP = value;
	
} else if(name.compare("Cud1_aabbr_LNP") == 0 ) { 
	
	Cud1_aabbr_LNP = value;
	
} else if(name.compare("Cud1_aa33r_LNP") == 0 ) { 
	
	Cud1_aa33r_LNP = value;
	
} else if(name.compare("Cud1_33aar_LNP") == 0 ) { 
	
	Cud1_33aar_LNP = value;
	
} else if(name.compare("Cud1_3333r_LNP") == 0 ) { 
	
	Cud1_3333r_LNP = value;
	
} else if(name.compare("Cud8_aabbr_LNP") == 0 ) { 
	
	Cud8_aabbr_LNP = value;
	
} else if(name.compare("Cud8_aa33r_LNP") == 0 ) { 
	
	Cud8_aa33r_LNP = value;
	
} else if(name.compare("Cud8_33aar_LNP") == 0 ) { 
	
	Cud8_33aar_LNP = value;
	
} else if(name.compare("Cud8_3333r_LNP") == 0 ) { 
	
	Cud8_3333r_LNP = value;
	
} else if(name.compare("Cqu1_aabbr_LNP") == 0 ) { 
	
	Cqu1_aabbr_LNP = value;
	
} else if(name.compare("Cqu1_aa33r_LNP") == 0 ) { 
	
	Cqu1_aa33r_LNP = value;
	
} else if(name.compare("Cqu1_33aar_LNP") == 0 ) { 
	
	Cqu1_33aar_LNP = value;
	
} else if(name.compare("Cqu1_3333r_LNP") == 0 ) { 
	
	Cqu1_3333r_LNP = value;
	
} else if(name.compare("Cqu8_aabbr_LNP") == 0 ) { 
	
	Cqu8_aabbr_LNP = value;
	
} else if(name.compare("Cqu8_aa33r_LNP") == 0 ) { 
	
	Cqu8_aa33r_LNP = value;
	
} else if(name.compare("Cqu8_33aar_LNP") == 0 ) { 
	
	Cqu8_33aar_LNP = value;
	
} else if(name.compare("Cqu8_3333r_LNP") == 0 ) { 
	
	Cqu8_3333r_LNP = value;
	
} else if(name.compare("Cqd1_aabbr_LNP") == 0 ) { 
	
	Cqd1_aabbr_LNP = value;
	
} else if(name.compare("Cqd1_aa33r_LNP") == 0 ) { 
	
	Cqd1_aa33r_LNP = value;
	
} else if(name.compare("Cqd1_33aar_LNP") == 0 ) { 
	
	Cqd1_33aar_LNP = value;
	
} else if(name.compare("Cqd1_3333r_LNP") == 0 ) { 
	
	Cqd1_3333r_LNP = value;
	
} else if(name.compare("Cqd8_aabbr_LNP") == 0 ) { 
	
	Cqd8_aabbr_LNP = value;
	
} else if(name.compare("Cqd8_aa33r_LNP") == 0 ) { 
	
	Cqd8_aa33r_LNP = value;
	
} else if(name.compare("Cqd8_33aar_LNP") == 0 ) { 
	
	Cqd8_33aar_LNP = value;
	
} else if(name.compare("Cqd8_3333r_LNP") == 0 ) { 
	
	Cqd8_3333r_LNP = value;
	
} else if(name.compare("Cquqd1_3333r_LNP") == 0 ) { 
	
	Cquqd1_3333r_LNP = value;
	
} else if(name.compare("Cquqd8_3333r_LNP") == 0 ) { 
	
	Cquqd8_3333r_LNP = value;
	
} else if(name.compare("Clequ1_3333r_LNP") == 0 ) { 
	
	Clequ1_3333r_LNP = value;
	
} else if(name.compare("Clequ3_3333r_LNP") == 0 ) { 
	
	Clequ3_3333r_LNP = value;
        
} else if (name.compare("Lambda_NP") == 0 )
        Lambda_NP = value;
    else
        NPSMEFTd6General::setParameter(name, value);

}


void NPSMEFTd6U2::setNPSMEFTd6GeneralParameters()
{
    //The names of some WC are the same in both classes, we could try to define again only those whose name do not coincide
    //with the one of the General class. However, I think it's clearer to define again all the WC which are relevant for the
    //model with the new symmetry. Since we redefine them we use the scope resolution operator to set the values of the WC in
    //the General class
    
    CHl1_11r_LNP = CHl1_aar_LNP;
    CHl1_22r_LNP = CHl1_aar_LNP;
    
    CHl3_11r_LNP = CHl3_aar_LNP;
    CHl3_22r_LNP = CHl3_aar_LNP;
    
    CHe_11r_LNP = CHe_aar_LNP;
    CHe_22r_LNP = CHe_aar_LNP;
    
    CHq1_11r_LNP = CHq1_aar_LNP;
    CHq1_22r_LNP = CHq1_aar_LNP;
    
    CHq3_11r_LNP = CHq3_aar_LNP;
    CHq3_22r_LNP = CHq3_aar_LNP;
    
    CHu_11r_LNP = CHu_aar_LNP;
    CHu_22r_LNP = CHu_aar_LNP;
    
    CHd_11r_LNP = CHd_aar_LNP;
    CHd_22r_LNP = CHd_aar_LNP;
    
    Cll_1111r_LNP = Cll_aabbr_LNP+Cll_abbar_LNP;
    Cll_2222r_LNP = Cll_aabbr_LNP+Cll_abbar_LNP;
    
    Cll_1122r_LNP = Cll_aabbr_LNP;
    
    Cll_1133r_LNP = Cll_aa33r_LNP;
    Cll_2233r_LNP = Cll_aa33r_LNP;
    
    Cll_1221r_LNP = Cll_abbar_LNP;
    
    Cll_1331r_LNP = Cll_a33ar_LNP;
    Cll_2332r_LNP = Cll_a33ar_LNP;
    
    Clq1_1111r_LNP = Clq1_aabbr_LNP;
    Clq1_1122r_LNP = Clq1_aabbr_LNP;
    Clq1_2211r_LNP = Clq1_aabbr_LNP;
    Clq1_2222r_LNP = Clq1_aabbr_LNP;
    
    Clq1_1133r_LNP = Clq1_aa33r_LNP;
    Clq1_2233r_LNP = Clq1_aa33r_LNP;
    
    Clq1_3311r_LNP = Clq1_33aar_LNP;
    Clq1_3322r_LNP = Clq1_33aar_LNP;
    
    Clq3_1111r_LNP = Clq3_aabbr_LNP;
    Clq3_1122r_LNP = Clq3_aabbr_LNP;
    Clq3_2211r_LNP = Clq3_aabbr_LNP;
    Clq3_2222r_LNP = Clq3_aabbr_LNP;
    
    Clq3_1133r_LNP = Clq3_aa33r_LNP;
    Clq3_2233r_LNP = Clq3_aa33r_LNP;
    
    Clq3_3311r_LNP = Clq3_33aar_LNP;
    Clq3_3322r_LNP = Clq3_33aar_LNP;
    
    Cee_1111r_LNP = Cee_aabbr_LNP;
    Cee_1122r_LNP = Cee_aabbr_LNP;
    Cee_2222r_LNP = Cee_aabbr_LNP;
    
    Cee_1133r_LNP = Cee_aa33r_LNP;
    Cee_2233r_LNP = Cee_aa33r_LNP;
    
    Ceu_1111r_LNP = Ceu_aabbr_LNP;
    Ceu_1122r_LNP = Ceu_aabbr_LNP;
    
    Ceu_2211r_LNP = Ceu_aabbr_LNP;
    Ceu_2222r_LNP = Ceu_aabbr_LNP;
    
    Ceu_1133r_LNP = Ceu_aa33r_LNP;
    Ceu_2233r_LNP = Ceu_aa33r_LNP;
    
    Ceu_3311r_LNP = Ceu_33aar_LNP;
    Ceu_3322r_LNP = Ceu_33aar_LNP;
    
    Ced_1111r_LNP = Ced_aabbr_LNP;
    Ced_1122r_LNP = Ced_aabbr_LNP;
    
    Ced_2211r_LNP = Ced_aabbr_LNP;
    Ced_2222r_LNP = Ced_aabbr_LNP;
    
    Ced_1133r_LNP = Ced_aa33r_LNP;
    Ced_2233r_LNP = Ced_aa33r_LNP;
    
    Ced_3311r_LNP = Ced_33aar_LNP;
    Ced_3322r_LNP = Ced_33aar_LNP;
    
    Cle_1111r_LNP = Cle_aabbr_LNP;
    Cle_1122r_LNP = Cle_aabbr_LNP;
    
    Cle_2211r_LNP = Cle_aabbr_LNP;
    Cle_2222r_LNP = Cle_aabbr_LNP;
    
    Cle_1133r_LNP = Cle_aa33r_LNP;
    Cle_2233r_LNP = Cle_aa33r_LNP;
    
    Cle_3311r_LNP = Cle_33aar_LNP;
    Cle_3322r_LNP = Cle_33aar_LNP;
    
    Clu_1111r_LNP = Clu_aabbr_LNP;
    Clu_1122r_LNP = Clu_aabbr_LNP;
    
    Clu_2211r_LNP = Clu_aabbr_LNP;
    Clu_2222r_LNP = Clu_aabbr_LNP;
    
    Clu_1133r_LNP = Clu_aa33r_LNP;
    Clu_2233r_LNP = Clu_aa33r_LNP;
    
    Clu_3311r_LNP = Clu_33aar_LNP;
    Clu_3322r_LNP = Clu_33aar_LNP;
    
    Cld_1111r_LNP = Cld_aabbr_LNP;
    Cld_1122r_LNP = Cld_aabbr_LNP;
    
    Cld_2211r_LNP = Cld_aabbr_LNP;
    Cld_2222r_LNP = Cld_aabbr_LNP;
    
    Cld_1133r_LNP = Cld_aa33r_LNP;
    Cld_2233r_LNP = Cld_aa33r_LNP;
    
    Cld_3311r_LNP = Cld_33aar_LNP;
    Cld_3322r_LNP = Cld_33aar_LNP;
    
    Cqe_1111r_LNP = Cqe_aabbr_LNP;
    Cqe_1122r_LNP = Cqe_aabbr_LNP;
    
    Cqe_2211r_LNP = Cqe_aabbr_LNP;
    Cqe_2222r_LNP = Cqe_aabbr_LNP;
    
    Cqe_1133r_LNP = Cqe_aa33r_LNP;
    Cqe_2233r_LNP = Cqe_aa33r_LNP;
    
    Cqe_3311r_LNP = Cqe_33aar_LNP;
    Cqe_3322r_LNP = Cqe_33aar_LNP;
    
    Cqq1_1111r_LNP = Cqq1_aabbr_LNP+Cqq1_abbar_LNP;
    Cqq1_2222r_LNP = Cqq1_aabbr_LNP+Cqq1_abbar_LNP;
    
    Cqq1_1122r_LNP = Cqq1_aabbr_LNP;
    
    Cqq1_1133r_LNP = Cqq1_aa33r_LNP;
    Cqq1_2233r_LNP = Cqq1_aa33r_LNP;
    Cqq1_1221r_LNP = Cqq1_abbar_LNP;
    
    Cqq1_1331r_LNP = Cqq1_a33ar_LNP;
    Cqq1_2332r_LNP = Cqq1_a33ar_LNP;
     
    Cqq3_1111r_LNP = Cqq3_aabbr_LNP+Cqq3_abbar_LNP;
    Cqq3_2222r_LNP = Cqq3_aabbr_LNP+Cqq3_abbar_LNP;
    
    Cqq3_1122r_LNP = Cqq3_aabbr_LNP;
    
    Cqq3_1133r_LNP = Cqq3_aa33r_LNP;
    Cqq3_2233r_LNP = Cqq3_aa33r_LNP;

    Cqq3_1221r_LNP = Cqq3_abbar_LNP;
    
    Cqq3_1331r_LNP = Cqq3_a33ar_LNP;
    Cqq3_2332r_LNP = Cqq3_a33ar_LNP;
   
    Cuu_1111r_LNP = Cuu_aabbr_LNP+Cuu_abbar_LNP;
    Cuu_2222r_LNP = Cuu_aabbr_LNP+Cuu_abbar_LNP;
    
    Cuu_1122r_LNP = Cuu_aabbr_LNP;
    
    Cuu_1133r_LNP = Cuu_aa33r_LNP;
    Cuu_2233r_LNP = Cuu_aa33r_LNP;

    Cuu_1221r_LNP = Cuu_abbar_LNP;
    
    Cuu_1331r_LNP = Cuu_a33ar_LNP;
    Cuu_2332r_LNP = Cuu_a33ar_LNP;   
    
    Cdd_1111r_LNP = Cdd_aabbr_LNP+Cdd_abbar_LNP;
    Cdd_2222r_LNP = Cdd_aabbr_LNP+Cdd_abbar_LNP;
    
    Cdd_1122r_LNP = Cdd_aabbr_LNP;
    
    Cdd_1133r_LNP = Cdd_aa33r_LNP;
    Cdd_2233r_LNP = Cdd_aa33r_LNP;

    Cdd_1221r_LNP = Cdd_abbar_LNP;
    
    Cdd_1331r_LNP = Cdd_a33ar_LNP;
    Cdd_2332r_LNP = Cdd_a33ar_LNP;    
    
    Cud1_1111r_LNP = Cud1_aabbr_LNP;
    Cud1_1122r_LNP = Cud1_aabbr_LNP;
    
    Cud1_2211r_LNP = Cud1_aabbr_LNP;
    Cud1_2222r_LNP = Cud1_aabbr_LNP;
    
    Cud1_1133r_LNP = Cud1_aa33r_LNP;
    Cud1_2233r_LNP = Cud1_aa33r_LNP;
    
    Cud1_3311r_LNP = Cud1_33aar_LNP;
    Cud1_3322r_LNP = Cud1_33aar_LNP;
    
    Cud8_1111r_LNP = Cud8_aabbr_LNP;
    Cud8_1122r_LNP = Cud8_aabbr_LNP;
    
    Cud8_2211r_LNP = Cud8_aabbr_LNP;
    Cud8_2222r_LNP = Cud8_aabbr_LNP;
    
    Cud8_1133r_LNP = Cud8_aa33r_LNP;
    Cud8_2233r_LNP = Cud8_aa33r_LNP;
    
    Cud8_3311r_LNP = Cud8_33aar_LNP;
    Cud8_3322r_LNP = Cud8_33aar_LNP;
    
    Cqu1_1111r_LNP = Cqu1_aabbr_LNP;
    Cqu1_1122r_LNP = Cqu1_aabbr_LNP;
    
    Cqu1_2211r_LNP = Cqu1_aabbr_LNP;
    Cqu1_2222r_LNP = Cqu1_aabbr_LNP;
    
    Cqu1_1133r_LNP = Cqu1_aa33r_LNP;
    Cqu1_2233r_LNP = Cqu1_aa33r_LNP;
    
    Cqu1_3311r_LNP = Cqu1_33aar_LNP;
    Cqu1_3322r_LNP = Cqu1_33aar_LNP;
    
    Cqu8_1111r_LNP = Cqu8_aabbr_LNP;
    Cqu8_1122r_LNP = Cqu8_aabbr_LNP;
    
    Cqu8_2211r_LNP = Cqu8_aabbr_LNP;
    Cqu8_2222r_LNP = Cqu8_aabbr_LNP;
    
    Cqu8_1133r_LNP = Cqu8_aa33r_LNP;
    Cqu8_2233r_LNP = Cqu8_aa33r_LNP;
    
    Cqu8_3311r_LNP = Cqu8_33aar_LNP;
    Cqu8_3322r_LNP = Cqu8_33aar_LNP;
    
    Cqd1_1111r_LNP = Cqd1_aabbr_LNP;
    Cqd1_1122r_LNP = Cqd1_aabbr_LNP;
    
    Cqd1_2211r_LNP = Cqd1_aabbr_LNP;
    Cqd1_2222r_LNP = Cqd1_aabbr_LNP;
    
    Cqd1_1133r_LNP = Cqd1_aa33r_LNP;
    Cqd1_2233r_LNP = Cqd1_aa33r_LNP;
    
    Cqd1_3311r_LNP = Cqd1_33aar_LNP;
    Cqd1_3322r_LNP = Cqd1_33aar_LNP;
    
    Cqd8_1111r_LNP = Cqd8_aabbr_LNP;
    Cqd8_1122r_LNP = Cqd8_aabbr_LNP;
    
    Cqd8_2211r_LNP = Cqd8_aabbr_LNP;
    Cqd8_2222r_LNP = Cqd8_aabbr_LNP;
    
    Cqd8_1133r_LNP = Cqd8_aa33r_LNP;
    Cqd8_2233r_LNP = Cqd8_aa33r_LNP;
    
    Cqd8_3311r_LNP = Cqd8_33aar_LNP;
    Cqd8_3322r_LNP = Cqd8_33aar_LNP;    
    
    
}



bool NPSMEFTd6U2::PostUpdate()
{
    NPSMEFTd6General::GenerateSMInitialConditions();
    
    setNPSMEFTd6GeneralParameters();
    
    if (!NPSMEFTd6General::PostUpdate()) return (false);
    
    return (true);
        
}
