/* 
 * Copyright (C) 2012 HEPfit Collaboration
 *
 *
 * For the licensing terms see doc/COPYING.
 * 
 * Created on 18 September 2023, 16:23
 * 
 */


#include "NPSMEFTd6MFV.h"


std::string NPSMEFTd6MFV::NPSMEFTd6MFVVars[NNPSMEFTd6MFVVars] = {
	"CG_LNP","CW_LNP","CHG_LNP","CHW_LNP","CHB_LNP","CHWB_LNP","CHD_LNP","CHbox_LNP",
	"CH_LNP","CHl1_LNP","CHl3_LNP","CHe_LNP","Cll_aabb_LNP","Cll_abba_LNP","Cee_LNP","Cle_LNP",
	"CuH_0_LNP","CuH_u_LNP","CuH_d_LNP","CuG_0_LNP","CuG_u_LNP","CuG_d_LNP","CuW_0_LNP","CuW_u_LNP",
	"CuW_d_LNP","CuB_0_LNP","CuB_u_LNP","CuB_d_LNP","CdH_0_LNP","CdH_u_LNP","CdH_d_LNP","CdG_0_LNP",
	"CdG_u_LNP","CdG_d_LNP","CdW_0_LNP","CdW_u_LNP","CdW_d_LNP","CdB_0_LNP","CdB_u_LNP","CdB_d_LNP",
	"CHq1_0_LNP","CHq1_u_LNP","CHq1_d_LNP","CHq3_0_LNP","CHq3_u_LNP","CHq3_d_LNP","CHu_0_LNP","CHu_u_LNP",
	"CHd_0_LNP","CHd_d_LNP","CHud_ud_LNP","Clq1_0_LNP","Clq1_u_LNP","Clq1_d_LNP","Clq3_0_LNP","Clq3_u_LNP",
	"Clq3_d_LNP","Cqe_0_LNP","Cqe_u_LNP","Cqe_d_LNP","Clu_0_LNP","Clu_u_LNP","Ceu_0_LNP","Ceu_u_LNP",
	"Cld_0_LNP","Cld_d_LNP","Ced_0_LNP","Ced_d_LNP","Cqq1_00_LNP","Cqq1_0u_LNP","Cqq1_0d_LNP","Cqq1_u0_LNP",
	"Cqq1_uu_LNP","Cqq1_ud_LNP","Cqq1_d0_LNP","Cqq1_du_LNP","Cqq1_dd_LNP","Cqq3_00_LNP","Cqq3_0u_LNP","Cqq3_0d_LNP",
	"Cqq3_u0_LNP","Cqq3_uu_LNP","Cqq3_ud_LNP","Cqq3_d0_LNP","Cqq3_du_LNP","Cqq3_dd_LNP","Cuu_00_LNP","Cuu_0u_LNP",
	"Cuu_u0_LNP","Cuu_uu_LNP","Cdd_00_LNP","Cdd_0d_LNP","Cdd_d0_LNP","Cdd_dd_LNP","Cud1_00_LNP","Cud1_u0_LNP",
	"Cud1_0d_LNP","Cud1_ud_LNP","Cud8_00_LNP","Cud8_u0_LNP","Cud8_0d_LNP","Cud8_ud_LNP","Cqu1_00_LNP","Cqu1_u0_LNP",
	"Cqu1_d0_LNP","Cqu1_0u_LNP","Cqu1_uu_LNP","Cqu1_du_LNP","Cqu8_00_LNP","Cqu8_u0_LNP","Cqu8_d0_LNP","Cqu8_0u_LNP",
	"Cqu8_uu_LNP","Cqu8_du_LNP","Cqd1_00_LNP","Cqd1_u0_LNP","Cqd1_d0_LNP","Cqd1_0d_LNP","Cqd1_ud_LNP","Cqd1_dd_LNP",
	"Cqd8_00_LNP","Cqd8_u0_LNP","Cqd8_d0_LNP","Cqd8_0d_LNP","Cqd8_ud_LNP","Cqd8_dd_LNP","Cquqd1_00_LNP","Cquqd8_00_LNP",
    "Lambda_NP"
};

NPSMEFTd6MFV::NPSMEFTd6MFV(): NPSMEFTd6General() {
    setModelName("NPSMEFTd6MFV");
	ModelParamMap.insert(std::make_pair("CG_LNP",std::cref(CG_LNP)));
	ModelParamMap.insert(std::make_pair("CW_LNP",std::cref(CW_LNP)));
	ModelParamMap.insert(std::make_pair("CHG_LNP",std::cref(CHG_LNP)));
	ModelParamMap.insert(std::make_pair("CHW_LNP",std::cref(CHW_LNP)));
	ModelParamMap.insert(std::make_pair("CHB_LNP",std::cref(CHB_LNP)));
	ModelParamMap.insert(std::make_pair("CHWB_LNP",std::cref(CHWB_LNP)));
	ModelParamMap.insert(std::make_pair("CHD_LNP",std::cref(CHD_LNP)));
	ModelParamMap.insert(std::make_pair("CHbox_LNP",std::cref(CHbox_LNP)));
	ModelParamMap.insert(std::make_pair("CH_LNP",std::cref(CH_LNP)));
	ModelParamMap.insert(std::make_pair("CHl1_LNP",std::cref(CHl1_LNP)));
	ModelParamMap.insert(std::make_pair("CHl3_LNP",std::cref(CHl3_LNP)));
	ModelParamMap.insert(std::make_pair("CHe_LNP",std::cref(CHe_LNP)));
	ModelParamMap.insert(std::make_pair("Cll_aabb_LNP",std::cref(Cll_aabb_LNP)));
	ModelParamMap.insert(std::make_pair("Cll_abba_LNP",std::cref(Cll_abba_LNP)));
	ModelParamMap.insert(std::make_pair("Cee_LNP",std::cref(Cee_LNP)));
	ModelParamMap.insert(std::make_pair("Cle_LNP",std::cref(Cle_LNP)));
	ModelParamMap.insert(std::make_pair("CuH_0_LNP",std::cref(CuH_0_LNP)));
	ModelParamMap.insert(std::make_pair("CuH_u_LNP",std::cref(CuH_u_LNP)));
	ModelParamMap.insert(std::make_pair("CuH_d_LNP",std::cref(CuH_d_LNP)));
	ModelParamMap.insert(std::make_pair("CuG_0_LNP",std::cref(CuG_0_LNP)));
	ModelParamMap.insert(std::make_pair("CuG_u_LNP",std::cref(CuG_u_LNP)));
	ModelParamMap.insert(std::make_pair("CuG_d_LNP",std::cref(CuG_d_LNP)));
	ModelParamMap.insert(std::make_pair("CuW_0_LNP",std::cref(CuW_0_LNP)));
	ModelParamMap.insert(std::make_pair("CuW_u_LNP",std::cref(CuW_u_LNP)));
	ModelParamMap.insert(std::make_pair("CuW_d_LNP",std::cref(CuW_d_LNP)));
	ModelParamMap.insert(std::make_pair("CuB_0_LNP",std::cref(CuB_0_LNP)));
	ModelParamMap.insert(std::make_pair("CuB_u_LNP",std::cref(CuB_u_LNP)));
	ModelParamMap.insert(std::make_pair("CuB_d_LNP",std::cref(CuB_d_LNP)));
	ModelParamMap.insert(std::make_pair("CdH_0_LNP",std::cref(CdH_0_LNP)));
	ModelParamMap.insert(std::make_pair("CdH_u_LNP",std::cref(CdH_u_LNP)));
	ModelParamMap.insert(std::make_pair("CdH_d_LNP",std::cref(CdH_d_LNP)));
	ModelParamMap.insert(std::make_pair("CdG_0_LNP",std::cref(CdG_0_LNP)));
	ModelParamMap.insert(std::make_pair("CdG_u_LNP",std::cref(CdG_u_LNP)));
	ModelParamMap.insert(std::make_pair("CdG_d_LNP",std::cref(CdG_d_LNP)));
	ModelParamMap.insert(std::make_pair("CdW_0_LNP",std::cref(CdW_0_LNP)));
	ModelParamMap.insert(std::make_pair("CdW_u_LNP",std::cref(CdW_u_LNP)));
	ModelParamMap.insert(std::make_pair("CdW_d_LNP",std::cref(CdW_d_LNP)));
	ModelParamMap.insert(std::make_pair("CdB_0_LNP",std::cref(CdB_0_LNP)));
	ModelParamMap.insert(std::make_pair("CdB_u_LNP",std::cref(CdB_u_LNP)));
	ModelParamMap.insert(std::make_pair("CdB_d_LNP",std::cref(CdB_d_LNP)));
	ModelParamMap.insert(std::make_pair("CHq1_0_LNP",std::cref(CHq1_0_LNP)));
	ModelParamMap.insert(std::make_pair("CHq1_u_LNP",std::cref(CHq1_u_LNP)));
	ModelParamMap.insert(std::make_pair("CHq1_d_LNP",std::cref(CHq1_d_LNP)));
	ModelParamMap.insert(std::make_pair("CHq3_0_LNP",std::cref(CHq3_0_LNP)));
	ModelParamMap.insert(std::make_pair("CHq3_u_LNP",std::cref(CHq3_u_LNP)));
	ModelParamMap.insert(std::make_pair("CHq3_d_LNP",std::cref(CHq3_d_LNP)));
	ModelParamMap.insert(std::make_pair("CHu_0_LNP",std::cref(CHu_0_LNP)));
	ModelParamMap.insert(std::make_pair("CHu_u_LNP",std::cref(CHu_u_LNP)));
	ModelParamMap.insert(std::make_pair("CHd_0_LNP",std::cref(CHd_0_LNP)));
	ModelParamMap.insert(std::make_pair("CHd_d_LNP",std::cref(CHd_d_LNP)));
	ModelParamMap.insert(std::make_pair("CHud_ud_LNP",std::cref(CHud_ud_LNP)));
	ModelParamMap.insert(std::make_pair("Clq1_0_LNP",std::cref(Clq1_0_LNP)));
	ModelParamMap.insert(std::make_pair("Clq1_u_LNP",std::cref(Clq1_u_LNP)));
	ModelParamMap.insert(std::make_pair("Clq1_d_LNP",std::cref(Clq1_d_LNP)));
	ModelParamMap.insert(std::make_pair("Clq3_0_LNP",std::cref(Clq3_0_LNP)));
	ModelParamMap.insert(std::make_pair("Clq3_u_LNP",std::cref(Clq3_u_LNP)));
	ModelParamMap.insert(std::make_pair("Clq3_d_LNP",std::cref(Clq3_d_LNP)));
	ModelParamMap.insert(std::make_pair("Cqe_0_LNP",std::cref(Cqe_0_LNP)));
	ModelParamMap.insert(std::make_pair("Cqe_u_LNP",std::cref(Cqe_u_LNP)));
	ModelParamMap.insert(std::make_pair("Cqe_d_LNP",std::cref(Cqe_d_LNP)));
	ModelParamMap.insert(std::make_pair("Clu_0_LNP",std::cref(Clu_0_LNP)));
	ModelParamMap.insert(std::make_pair("Clu_u_LNP",std::cref(Clu_u_LNP)));
	ModelParamMap.insert(std::make_pair("Ceu_0_LNP",std::cref(Ceu_0_LNP)));
	ModelParamMap.insert(std::make_pair("Ceu_u_LNP",std::cref(Ceu_u_LNP)));
	ModelParamMap.insert(std::make_pair("Cld_0_LNP",std::cref(Cld_0_LNP)));
	ModelParamMap.insert(std::make_pair("Cld_d_LNP",std::cref(Cld_d_LNP)));
	ModelParamMap.insert(std::make_pair("Ced_0_LNP",std::cref(Ced_0_LNP)));
	ModelParamMap.insert(std::make_pair("Ced_d_LNP",std::cref(Ced_d_LNP)));
	ModelParamMap.insert(std::make_pair("Cqq1_00_LNP",std::cref(Cqq1_00_LNP)));
	ModelParamMap.insert(std::make_pair("Cqq1_0u_LNP",std::cref(Cqq1_0u_LNP)));
	ModelParamMap.insert(std::make_pair("Cqq1_0d_LNP",std::cref(Cqq1_0d_LNP)));
	ModelParamMap.insert(std::make_pair("Cqq1_u0_LNP",std::cref(Cqq1_u0_LNP)));
	ModelParamMap.insert(std::make_pair("Cqq1_uu_LNP",std::cref(Cqq1_uu_LNP)));
	ModelParamMap.insert(std::make_pair("Cqq1_ud_LNP",std::cref(Cqq1_ud_LNP)));
	ModelParamMap.insert(std::make_pair("Cqq1_d0_LNP",std::cref(Cqq1_d0_LNP)));
	ModelParamMap.insert(std::make_pair("Cqq1_du_LNP",std::cref(Cqq1_du_LNP)));
	ModelParamMap.insert(std::make_pair("Cqq1_dd_LNP",std::cref(Cqq1_dd_LNP)));
	ModelParamMap.insert(std::make_pair("Cqq3_00_LNP",std::cref(Cqq3_00_LNP)));
	ModelParamMap.insert(std::make_pair("Cqq3_0u_LNP",std::cref(Cqq3_0u_LNP)));
	ModelParamMap.insert(std::make_pair("Cqq3_0d_LNP",std::cref(Cqq3_0d_LNP)));
	ModelParamMap.insert(std::make_pair("Cqq3_u0_LNP",std::cref(Cqq3_u0_LNP)));
	ModelParamMap.insert(std::make_pair("Cqq3_uu_LNP",std::cref(Cqq3_uu_LNP)));
	ModelParamMap.insert(std::make_pair("Cqq3_ud_LNP",std::cref(Cqq3_ud_LNP)));
	ModelParamMap.insert(std::make_pair("Cqq3_d0_LNP",std::cref(Cqq3_d0_LNP)));
	ModelParamMap.insert(std::make_pair("Cqq3_du_LNP",std::cref(Cqq3_du_LNP)));
	ModelParamMap.insert(std::make_pair("Cqq3_dd_LNP",std::cref(Cqq3_dd_LNP)));
	ModelParamMap.insert(std::make_pair("Cuu_00_LNP",std::cref(Cuu_00_LNP)));
	ModelParamMap.insert(std::make_pair("Cuu_0u_LNP",std::cref(Cuu_0u_LNP)));
	ModelParamMap.insert(std::make_pair("Cuu_u0_LNP",std::cref(Cuu_u0_LNP)));
	ModelParamMap.insert(std::make_pair("Cuu_uu_LNP",std::cref(Cuu_uu_LNP)));
	ModelParamMap.insert(std::make_pair("Cdd_00_LNP",std::cref(Cdd_00_LNP)));
	ModelParamMap.insert(std::make_pair("Cdd_0d_LNP",std::cref(Cdd_0d_LNP)));
	ModelParamMap.insert(std::make_pair("Cdd_d0_LNP",std::cref(Cdd_d0_LNP)));
	ModelParamMap.insert(std::make_pair("Cdd_dd_LNP",std::cref(Cdd_dd_LNP)));
	ModelParamMap.insert(std::make_pair("Cud1_00_LNP",std::cref(Cud1_00_LNP)));
	ModelParamMap.insert(std::make_pair("Cud1_u0_LNP",std::cref(Cud1_u0_LNP)));
	ModelParamMap.insert(std::make_pair("Cud1_0d_LNP",std::cref(Cud1_0d_LNP)));
	ModelParamMap.insert(std::make_pair("Cud1_ud_LNP",std::cref(Cud1_ud_LNP)));
	ModelParamMap.insert(std::make_pair("Cud8_00_LNP",std::cref(Cud8_00_LNP)));
	ModelParamMap.insert(std::make_pair("Cud8_u0_LNP",std::cref(Cud8_u0_LNP)));
	ModelParamMap.insert(std::make_pair("Cud8_0d_LNP",std::cref(Cud8_0d_LNP)));
	ModelParamMap.insert(std::make_pair("Cud8_ud_LNP",std::cref(Cud8_ud_LNP)));
	ModelParamMap.insert(std::make_pair("Cqu1_00_LNP",std::cref(Cqu1_00_LNP)));
	ModelParamMap.insert(std::make_pair("Cqu1_u0_LNP",std::cref(Cqu1_u0_LNP)));
	ModelParamMap.insert(std::make_pair("Cqu1_d0_LNP",std::cref(Cqu1_d0_LNP)));
	ModelParamMap.insert(std::make_pair("Cqu1_0u_LNP",std::cref(Cqu1_0u_LNP)));
	ModelParamMap.insert(std::make_pair("Cqu1_uu_LNP",std::cref(Cqu1_uu_LNP)));
	ModelParamMap.insert(std::make_pair("Cqu1_du_LNP",std::cref(Cqu1_du_LNP)));
	ModelParamMap.insert(std::make_pair("Cqu8_00_LNP",std::cref(Cqu8_00_LNP)));
	ModelParamMap.insert(std::make_pair("Cqu8_u0_LNP",std::cref(Cqu8_u0_LNP)));
	ModelParamMap.insert(std::make_pair("Cqu8_d0_LNP",std::cref(Cqu8_d0_LNP)));
	ModelParamMap.insert(std::make_pair("Cqu8_0u_LNP",std::cref(Cqu8_0u_LNP)));
	ModelParamMap.insert(std::make_pair("Cqu8_uu_LNP",std::cref(Cqu8_uu_LNP)));
	ModelParamMap.insert(std::make_pair("Cqu8_du_LNP",std::cref(Cqu8_du_LNP)));
	ModelParamMap.insert(std::make_pair("Cqd1_00_LNP",std::cref(Cqd1_00_LNP)));
	ModelParamMap.insert(std::make_pair("Cqd1_u0_LNP",std::cref(Cqd1_u0_LNP)));
	ModelParamMap.insert(std::make_pair("Cqd1_d0_LNP",std::cref(Cqd1_d0_LNP)));
	ModelParamMap.insert(std::make_pair("Cqd1_0d_LNP",std::cref(Cqd1_0d_LNP)));
	ModelParamMap.insert(std::make_pair("Cqd1_ud_LNP",std::cref(Cqd1_ud_LNP)));
	ModelParamMap.insert(std::make_pair("Cqd1_dd_LNP",std::cref(Cqd1_dd_LNP)));
	ModelParamMap.insert(std::make_pair("Cqd8_00_LNP",std::cref(Cqd8_00_LNP)));
	ModelParamMap.insert(std::make_pair("Cqd8_u0_LNP",std::cref(Cqd8_u0_LNP)));
	ModelParamMap.insert(std::make_pair("Cqd8_d0_LNP",std::cref(Cqd8_d0_LNP)));
	ModelParamMap.insert(std::make_pair("Cqd8_0d_LNP",std::cref(Cqd8_0d_LNP)));
	ModelParamMap.insert(std::make_pair("Cqd8_ud_LNP",std::cref(Cqd8_ud_LNP)));
	ModelParamMap.insert(std::make_pair("Cqd8_dd_LNP",std::cref(Cqd8_dd_LNP)));
	ModelParamMap.insert(std::make_pair("Cquqd1_00_LNP",std::cref(Cquqd1_00_LNP)));
	ModelParamMap.insert(std::make_pair("Cquqd8_00_LNP",std::cref(Cquqd8_00_LNP)));

}

bool NPSMEFTd6MFV::Init(const std::map<std::string, double>& DPars) {
	if (SMEFTBasisFlag == "UP")
		throw std::runtime_error("Bad argument in SMEFTBasisFlag. (Only DOWN is allowed for this Model)");

    return (NPSMEFTd6General::Init(DPars));
}

void NPSMEFTd6MFV::setParameter(const std::string name, const double& value) {
	if (name.compare("CG_LNP") == 0) {
		CG_LNP = value;
	} else if (name.compare("CW_LNP") == 0) {
		CW_LNP = value;
	} else if (name.compare("CHG_LNP") == 0) {
		CHG_LNP = value;
	} else if (name.compare("CHW_LNP") == 0) {
		CHW_LNP = value;
	} else if (name.compare("CHB_LNP") == 0) {
		CHB_LNP = value;
	} else if (name.compare("CHWB_LNP") == 0) {
		CHWB_LNP = value;
	} else if (name.compare("CHD_LNP") == 0) {
		CHD_LNP = value;
	} else if (name.compare("CHbox_LNP") == 0) {
		CHbox_LNP = value;
	} else if (name.compare("CH_LNP") == 0) {
		CH_LNP = value;
	} else if (name.compare("CHl1_LNP") == 0) {
		CHl1_LNP = value;
	} else if (name.compare("CHl3_LNP") == 0) {
		CHl3_LNP = value;
	} else if (name.compare("CHe_LNP") == 0) {
		CHe_LNP = value;
	} else if (name.compare("Cll_aabb_LNP") == 0) {
		Cll_aabb_LNP = value;
	} else if (name.compare("Cll_abba_LNP") == 0) {
		Cll_abba_LNP = value;
	} else if (name.compare("Cee_LNP") == 0) {
		Cee_LNP = value;
	} else if (name.compare("Cle_LNP") == 0) {
		Cle_LNP = value;
	} else if (name.compare("CuH_0_LNP") == 0) {
		CuH_0_LNP = value;
	} else if (name.compare("CuH_u_LNP") == 0) {
		CuH_u_LNP = value;
	} else if (name.compare("CuH_d_LNP") == 0) {
		CuH_d_LNP = value;
	} else if (name.compare("CuG_0_LNP") == 0) {
		CuG_0_LNP = value;
	} else if (name.compare("CuG_u_LNP") == 0) {
		CuG_u_LNP = value;
	} else if (name.compare("CuG_d_LNP") == 0) {
		CuG_d_LNP = value;
	} else if (name.compare("CuW_0_LNP") == 0) {
		CuW_0_LNP = value;
	} else if (name.compare("CuW_u_LNP") == 0) {
		CuW_u_LNP = value;
	} else if (name.compare("CuW_d_LNP") == 0) {
		CuW_d_LNP = value;
	} else if (name.compare("CuB_0_LNP") == 0) {
		CuB_0_LNP = value;
	} else if (name.compare("CuB_u_LNP") == 0) {
		CuB_u_LNP = value;
	} else if (name.compare("CuB_d_LNP") == 0) {
		CuB_d_LNP = value;
	} else if (name.compare("CdH_0_LNP") == 0) {
		CdH_0_LNP = value;
	} else if (name.compare("CdH_u_LNP") == 0) {
		CdH_u_LNP = value;
	} else if (name.compare("CdH_d_LNP") == 0) {
		CdH_d_LNP = value;
	} else if (name.compare("CdG_0_LNP") == 0) {
		CdG_0_LNP = value;
	} else if (name.compare("CdG_u_LNP") == 0) {
		CdG_u_LNP = value;
	} else if (name.compare("CdG_d_LNP") == 0) {
		CdG_d_LNP = value;
	} else if (name.compare("CdW_0_LNP") == 0) {
		CdW_0_LNP = value;
	} else if (name.compare("CdW_u_LNP") == 0) {
		CdW_u_LNP = value;
	} else if (name.compare("CdW_d_LNP") == 0) {
		CdW_d_LNP = value;
	} else if (name.compare("CdB_0_LNP") == 0) {
		CdB_0_LNP = value;
	} else if (name.compare("CdB_u_LNP") == 0) {
		CdB_u_LNP = value;
	} else if (name.compare("CdB_d_LNP") == 0) {
		CdB_d_LNP = value;
	} else if (name.compare("CHq1_0_LNP") == 0) {
		CHq1_0_LNP = value;
	} else if (name.compare("CHq1_u_LNP") == 0) {
		CHq1_u_LNP = value;
	} else if (name.compare("CHq1_d_LNP") == 0) {
		CHq1_d_LNP = value;
	} else if (name.compare("CHq3_0_LNP") == 0) {
		CHq3_0_LNP = value;
	} else if (name.compare("CHq3_u_LNP") == 0) {
		CHq3_u_LNP = value;
	} else if (name.compare("CHq3_d_LNP") == 0) {
		CHq3_d_LNP = value;
	} else if (name.compare("CHu_0_LNP") == 0) {
		CHu_0_LNP = value;
	} else if (name.compare("CHu_u_LNP") == 0) {
		CHu_u_LNP = value;
	} else if (name.compare("CHd_0_LNP") == 0) {
		CHd_0_LNP = value;
	} else if (name.compare("CHd_d_LNP") == 0) {
		CHd_d_LNP = value;
	} else if (name.compare("CHud_ud_LNP") == 0) {
		CHud_ud_LNP = value;
	} else if (name.compare("Clq1_0_LNP") == 0) {
		Clq1_0_LNP = value;
	} else if (name.compare("Clq1_u_LNP") == 0) {
		Clq1_u_LNP = value;
	} else if (name.compare("Clq1_d_LNP") == 0) {
		Clq1_d_LNP = value;
	} else if (name.compare("Clq3_0_LNP") == 0) {
		Clq3_0_LNP = value;
	} else if (name.compare("Clq3_u_LNP") == 0) {
		Clq3_u_LNP = value;
	} else if (name.compare("Clq3_d_LNP") == 0) {
		Clq3_d_LNP = value;
	} else if (name.compare("Cqe_0_LNP") == 0) {
		Cqe_0_LNP = value;
	} else if (name.compare("Cqe_u_LNP") == 0) {
		Cqe_u_LNP = value;
	} else if (name.compare("Cqe_d_LNP") == 0) {
		Cqe_d_LNP = value;
	} else if (name.compare("Clu_0_LNP") == 0) {
		Clu_0_LNP = value;
	} else if (name.compare("Clu_u_LNP") == 0) {
		Clu_u_LNP = value;
	} else if (name.compare("Ceu_0_LNP") == 0) {
		Ceu_0_LNP = value;
	} else if (name.compare("Ceu_u_LNP") == 0) {
		Ceu_u_LNP = value;
	} else if (name.compare("Cld_0_LNP") == 0) {
		Cld_0_LNP = value;
	} else if (name.compare("Cld_d_LNP") == 0) {
		Cld_d_LNP = value;
	} else if (name.compare("Ced_0_LNP") == 0) {
		Ced_0_LNP = value;
	} else if (name.compare("Ced_d_LNP") == 0) {
		Ced_d_LNP = value;
	} else if (name.compare("Cqq1_00_LNP") == 0) {
		Cqq1_00_LNP = value;
	} else if (name.compare("Cqq1_0u_LNP") == 0) {
		Cqq1_0u_LNP = value;
	} else if (name.compare("Cqq1_0d_LNP") == 0) {
		Cqq1_0d_LNP = value;
	} else if (name.compare("Cqq1_u0_LNP") == 0) {
		Cqq1_u0_LNP = value;
	} else if (name.compare("Cqq1_uu_LNP") == 0) {
		Cqq1_uu_LNP = value;
	} else if (name.compare("Cqq1_ud_LNP") == 0) {
		Cqq1_ud_LNP = value;
	} else if (name.compare("Cqq1_d0_LNP") == 0) {
		Cqq1_d0_LNP = value;
	} else if (name.compare("Cqq1_du_LNP") == 0) {
		Cqq1_du_LNP = value;
	} else if (name.compare("Cqq1_dd_LNP") == 0) {
		Cqq1_dd_LNP = value;
	} else if (name.compare("Cqq3_00_LNP") == 0) {
		Cqq3_00_LNP = value;
	} else if (name.compare("Cqq3_0u_LNP") == 0) {
		Cqq3_0u_LNP = value;
	} else if (name.compare("Cqq3_0d_LNP") == 0) {
		Cqq3_0d_LNP = value;
	} else if (name.compare("Cqq3_u0_LNP") == 0) {
		Cqq3_u0_LNP = value;
	} else if (name.compare("Cqq3_uu_LNP") == 0) {
		Cqq3_uu_LNP = value;
	} else if (name.compare("Cqq3_ud_LNP") == 0) {
		Cqq3_ud_LNP = value;
	} else if (name.compare("Cqq3_d0_LNP") == 0) {
		Cqq3_d0_LNP = value;
	} else if (name.compare("Cqq3_du_LNP") == 0) {
		Cqq3_du_LNP = value;
	} else if (name.compare("Cqq3_dd_LNP") == 0) {
		Cqq3_dd_LNP = value;
	} else if (name.compare("Cuu_00_LNP") == 0) {
		Cuu_00_LNP = value;
	} else if (name.compare("Cuu_0u_LNP") == 0) {
		Cuu_0u_LNP = value;
	} else if (name.compare("Cuu_u0_LNP") == 0) {
		Cuu_u0_LNP = value;
	} else if (name.compare("Cuu_uu_LNP") == 0) {
		Cuu_uu_LNP = value;
	} else if (name.compare("Cdd_00_LNP") == 0) {
		Cdd_00_LNP = value;
	} else if (name.compare("Cdd_0d_LNP") == 0) {
		Cdd_0d_LNP = value;
	} else if (name.compare("Cdd_d0_LNP") == 0) {
		Cdd_d0_LNP = value;
	} else if (name.compare("Cdd_dd_LNP") == 0) {
		Cdd_dd_LNP = value;
	} else if (name.compare("Cud1_00_LNP") == 0) {
		Cud1_00_LNP = value;
	} else if (name.compare("Cud1_u0_LNP") == 0) {
		Cud1_u0_LNP = value;
	} else if (name.compare("Cud1_0d_LNP") == 0) {
		Cud1_0d_LNP = value;
	} else if (name.compare("Cud1_ud_LNP") == 0) {
		Cud1_ud_LNP = value;
	} else if (name.compare("Cud8_00_LNP") == 0) {
		Cud8_00_LNP = value;
	} else if (name.compare("Cud8_u0_LNP") == 0) {
		Cud8_u0_LNP = value;
	} else if (name.compare("Cud8_0d_LNP") == 0) {
		Cud8_0d_LNP = value;
	} else if (name.compare("Cud8_ud_LNP") == 0) {
		Cud8_ud_LNP = value;
	} else if (name.compare("Cqu1_00_LNP") == 0) {
		Cqu1_00_LNP = value;
	} else if (name.compare("Cqu1_u0_LNP") == 0) {
		Cqu1_u0_LNP = value;
	} else if (name.compare("Cqu1_d0_LNP") == 0) {
		Cqu1_d0_LNP = value;
	} else if (name.compare("Cqu1_0u_LNP") == 0) {
		Cqu1_0u_LNP = value;
	} else if (name.compare("Cqu1_uu_LNP") == 0) {
		Cqu1_uu_LNP = value;
	} else if (name.compare("Cqu1_du_LNP") == 0) {
		Cqu1_du_LNP = value;
	} else if (name.compare("Cqu8_00_LNP") == 0) {
		Cqu8_00_LNP = value;
	} else if (name.compare("Cqu8_u0_LNP") == 0) {
		Cqu8_u0_LNP = value;
	} else if (name.compare("Cqu8_d0_LNP") == 0) {
		Cqu8_d0_LNP = value;
	} else if (name.compare("Cqu8_0u_LNP") == 0) {
		Cqu8_0u_LNP = value;
	} else if (name.compare("Cqu8_uu_LNP") == 0) {
		Cqu8_uu_LNP = value;
	} else if (name.compare("Cqu8_du_LNP") == 0) {
		Cqu8_du_LNP = value;
	} else if (name.compare("Cqd1_00_LNP") == 0) {
		Cqd1_00_LNP = value;
	} else if (name.compare("Cqd1_u0_LNP") == 0) {
		Cqd1_u0_LNP = value;
	} else if (name.compare("Cqd1_d0_LNP") == 0) {
		Cqd1_d0_LNP = value;
	} else if (name.compare("Cqd1_0d_LNP") == 0) {
		Cqd1_0d_LNP = value;
	} else if (name.compare("Cqd1_ud_LNP") == 0) {
		Cqd1_ud_LNP = value;
	} else if (name.compare("Cqd1_dd_LNP") == 0) {
		Cqd1_dd_LNP = value;
	} else if (name.compare("Cqd8_00_LNP") == 0) {
		Cqd8_00_LNP = value;
	} else if (name.compare("Cqd8_u0_LNP") == 0) {
		Cqd8_u0_LNP = value;
	} else if (name.compare("Cqd8_d0_LNP") == 0) {
		Cqd8_d0_LNP = value;
	} else if (name.compare("Cqd8_0d_LNP") == 0) {
		Cqd8_0d_LNP = value;
	} else if (name.compare("Cqd8_ud_LNP") == 0) {
		Cqd8_ud_LNP = value;
	} else if (name.compare("Cqd8_dd_LNP") == 0) {
		Cqd8_dd_LNP = value;
	} else if (name.compare("Cquqd1_00_LNP") == 0) {
		Cquqd1_00_LNP = value;
	} else if (name.compare("Cquqd8_00_LNP") == 0) {
		Cquqd8_00_LNP = value;
	}  else if (name.compare("Lambda_NP") == 0) {
        Lambda_NP = value;
    } else {
        NPSMEFTd6General::setParameter(name, value);
	}
}

void NPSMEFTd6MFV::setNPSMEFTd6GeneralParameters() {
	gslpp::vector<gslpp::complex> YdL(3, 0.), YddL(3, 0.), YdddL(3, 0.);
	gslpp::matrix<gslpp::complex> YuL(3, 3, 0.), YdML(3, 3, 0.);
	gslpp::matrix<gslpp::complex> YucuL(3, 3, 0.), YuucL(3, 3, 0.), YudL(3, 3, 0.);
	gslpp::matrix<gslpp::complex> YucuuL(3, 3, 0.), YdduL(3, 3, 0.), YucudL(3, 3, 0.);
	
	// Create the up-Yukawa matrix
    for (int i = 0; i < 3; i++) {
        for (int j = 0; j < 3; j++) {
            YuL.assignre(i, j, getSMEFTCoeffEW("YuR", i, j));
            YuL.assignim(i, j, getSMEFTCoeffEW("YuI", i, j));
        }
	}
	
	// Create the down-Yukawa diagonal matrix and its powers
	for (int i = 0; i < 3; i++) {
		double ydi = getSMEFTCoeffEW("YdR", i, i);
		YdML.assign(i, i, ydi);
		YdL.assign(i, ydi);
		YddL.assign(i, ydi*ydi);
		YdddL.assign(i, ydi*ydi*ydi);
	}
	
	// Products of two Yukawas
	YucuL = YuL.hconjugate() * YuL;
	YuucL = YuL * YuL.hconjugate();
	YudL = YuL * YdML;
	
	// Products of three Yukawas
	YucuuL = YucuL * YuL;
	YdduL = YdML * YdML * YuL;
	YucudL = YucuL * YdML;

	// Operator assignments

	CHl1_11r_LNP = CHl1_LNP;
	CHl1_22r_LNP = CHl1_LNP;
	CHl1_33r_LNP = CHl1_LNP;

	CHl3_11r_LNP = CHl3_LNP;
	CHl3_22r_LNP = CHl3_LNP;
	CHl3_33r_LNP = CHl3_LNP;

	CHe_11r_LNP = CHe_LNP;
	CHe_22r_LNP = CHe_LNP;
	CHe_33r_LNP = CHe_LNP;

	Cll_1111r_LNP = Cll_aabb_LNP + Cll_abba_LNP;
	Cll_2222r_LNP = Cll_aabb_LNP + Cll_abba_LNP;
	Cll_3333r_LNP = Cll_aabb_LNP + Cll_abba_LNP;
	Cll_1122r_LNP = Cll_aabb_LNP;
	Cll_1133r_LNP = Cll_aabb_LNP;
	Cll_2233r_LNP = Cll_aabb_LNP;
	Cll_1221r_LNP = Cll_abba_LNP;
	Cll_1331r_LNP = Cll_abba_LNP;
	Cll_2332r_LNP = Cll_abba_LNP;

	CuH_11r_LNP = (CuH_d_LNP*YdduL(0,0) + CuH_u_LNP*YucuuL(0,0) + CuH_0_LNP*YuL(0,0)).real();
	CuH_11i_LNP = (CuH_d_LNP*YdduL(0,0) + CuH_u_LNP*YucuuL(0,0) + CuH_0_LNP*YuL(0,0)).imag();
	CuH_12r_LNP = (CuH_d_LNP*YdduL(0,1) + CuH_u_LNP*YucuuL(0,1) + CuH_0_LNP*YuL(0,1)).real();
	CuH_12i_LNP = (CuH_d_LNP*YdduL(0,1) + CuH_u_LNP*YucuuL(0,1) + CuH_0_LNP*YuL(0,1)).imag();
	CuH_13r_LNP = (CuH_d_LNP*YdduL(0,2) + CuH_u_LNP*YucuuL(0,2) + CuH_0_LNP*YuL(0,2)).real();
	CuH_13i_LNP = (CuH_d_LNP*YdduL(0,2) + CuH_u_LNP*YucuuL(0,2) + CuH_0_LNP*YuL(0,2)).imag();
	CuH_21r_LNP = (CuH_d_LNP*YdduL(1,0) + CuH_u_LNP*YucuuL(1,0) + CuH_0_LNP*YuL(1,0)).real();
	CuH_21i_LNP = (CuH_d_LNP*YdduL(1,0) + CuH_u_LNP*YucuuL(1,0) + CuH_0_LNP*YuL(1,0)).imag();
	CuH_22r_LNP = (CuH_d_LNP*YdduL(1,1) + CuH_u_LNP*YucuuL(1,1) + CuH_0_LNP*YuL(1,1)).real();
	CuH_22i_LNP = (CuH_d_LNP*YdduL(1,1) + CuH_u_LNP*YucuuL(1,1) + CuH_0_LNP*YuL(1,1)).imag();
	CuH_23r_LNP = (CuH_d_LNP*YdduL(1,2) + CuH_u_LNP*YucuuL(1,2) + CuH_0_LNP*YuL(1,2)).real();
	CuH_23i_LNP = (CuH_d_LNP*YdduL(1,2) + CuH_u_LNP*YucuuL(1,2) + CuH_0_LNP*YuL(1,2)).imag();
	CuH_31r_LNP = (CuH_d_LNP*YdduL(2,0) + CuH_u_LNP*YucuuL(2,0) + CuH_0_LNP*YuL(2,0)).real();
	CuH_31i_LNP = (CuH_d_LNP*YdduL(2,0) + CuH_u_LNP*YucuuL(2,0) + CuH_0_LNP*YuL(2,0)).imag();
	CuH_32r_LNP = (CuH_d_LNP*YdduL(2,1) + CuH_u_LNP*YucuuL(2,1) + CuH_0_LNP*YuL(2,1)).real();
	CuH_32i_LNP = (CuH_d_LNP*YdduL(2,1) + CuH_u_LNP*YucuuL(2,1) + CuH_0_LNP*YuL(2,1)).imag();
	CuH_33r_LNP = (CuH_d_LNP*YdduL(2,2) + CuH_u_LNP*YucuuL(2,2) + CuH_0_LNP*YuL(2,2)).real();
	CuH_33i_LNP = (CuH_d_LNP*YdduL(2,2) + CuH_u_LNP*YucuuL(2,2) + CuH_0_LNP*YuL(2,2)).imag();

	CuG_11r_LNP = (CuG_d_LNP*YdduL(0,0) + CuG_u_LNP*YucuuL(0,0) + CuG_0_LNP*YuL(0,0)).real();
	CuG_11i_LNP = (CuG_d_LNP*YdduL(0,0) + CuG_u_LNP*YucuuL(0,0) + CuG_0_LNP*YuL(0,0)).imag();
	CuG_12r_LNP = (CuG_d_LNP*YdduL(0,1) + CuG_u_LNP*YucuuL(0,1) + CuG_0_LNP*YuL(0,1)).real();
	CuG_12i_LNP = (CuG_d_LNP*YdduL(0,1) + CuG_u_LNP*YucuuL(0,1) + CuG_0_LNP*YuL(0,1)).imag();
	CuG_13r_LNP = (CuG_d_LNP*YdduL(0,2) + CuG_u_LNP*YucuuL(0,2) + CuG_0_LNP*YuL(0,2)).real();
	CuG_13i_LNP = (CuG_d_LNP*YdduL(0,2) + CuG_u_LNP*YucuuL(0,2) + CuG_0_LNP*YuL(0,2)).imag();
	CuG_21r_LNP = (CuG_d_LNP*YdduL(1,0) + CuG_u_LNP*YucuuL(1,0) + CuG_0_LNP*YuL(1,0)).real();
	CuG_21i_LNP = (CuG_d_LNP*YdduL(1,0) + CuG_u_LNP*YucuuL(1,0) + CuG_0_LNP*YuL(1,0)).imag();
	CuG_22r_LNP = (CuG_d_LNP*YdduL(1,1) + CuG_u_LNP*YucuuL(1,1) + CuG_0_LNP*YuL(1,1)).real();
	CuG_22i_LNP = (CuG_d_LNP*YdduL(1,1) + CuG_u_LNP*YucuuL(1,1) + CuG_0_LNP*YuL(1,1)).imag();
	CuG_23r_LNP = (CuG_d_LNP*YdduL(1,2) + CuG_u_LNP*YucuuL(1,2) + CuG_0_LNP*YuL(1,2)).real();
	CuG_23i_LNP = (CuG_d_LNP*YdduL(1,2) + CuG_u_LNP*YucuuL(1,2) + CuG_0_LNP*YuL(1,2)).imag();
	CuG_31r_LNP = (CuG_d_LNP*YdduL(2,0) + CuG_u_LNP*YucuuL(2,0) + CuG_0_LNP*YuL(2,0)).real();
	CuG_31i_LNP = (CuG_d_LNP*YdduL(2,0) + CuG_u_LNP*YucuuL(2,0) + CuG_0_LNP*YuL(2,0)).imag();
	CuG_32r_LNP = (CuG_d_LNP*YdduL(2,1) + CuG_u_LNP*YucuuL(2,1) + CuG_0_LNP*YuL(2,1)).real();
	CuG_32i_LNP = (CuG_d_LNP*YdduL(2,1) + CuG_u_LNP*YucuuL(2,1) + CuG_0_LNP*YuL(2,1)).imag();
	CuG_33r_LNP = (CuG_d_LNP*YdduL(2,2) + CuG_u_LNP*YucuuL(2,2) + CuG_0_LNP*YuL(2,2)).real();
	CuG_33i_LNP = (CuG_d_LNP*YdduL(2,2) + CuG_u_LNP*YucuuL(2,2) + CuG_0_LNP*YuL(2,2)).imag();

	CuW_11r_LNP = (CuW_d_LNP*YdduL(0,0) + CuW_u_LNP*YucuuL(0,0) + CuW_0_LNP*YuL(0,0)).real();
	CuW_11i_LNP = (CuW_d_LNP*YdduL(0,0) + CuW_u_LNP*YucuuL(0,0) + CuW_0_LNP*YuL(0,0)).imag();
	CuW_12r_LNP = (CuW_d_LNP*YdduL(0,1) + CuW_u_LNP*YucuuL(0,1) + CuW_0_LNP*YuL(0,1)).real();
	CuW_12i_LNP = (CuW_d_LNP*YdduL(0,1) + CuW_u_LNP*YucuuL(0,1) + CuW_0_LNP*YuL(0,1)).imag();
	CuW_13r_LNP = (CuW_d_LNP*YdduL(0,2) + CuW_u_LNP*YucuuL(0,2) + CuW_0_LNP*YuL(0,2)).real();
	CuW_13i_LNP = (CuW_d_LNP*YdduL(0,2) + CuW_u_LNP*YucuuL(0,2) + CuW_0_LNP*YuL(0,2)).imag();
	CuW_21r_LNP = (CuW_d_LNP*YdduL(1,0) + CuW_u_LNP*YucuuL(1,0) + CuW_0_LNP*YuL(1,0)).real();
	CuW_21i_LNP = (CuW_d_LNP*YdduL(1,0) + CuW_u_LNP*YucuuL(1,0) + CuW_0_LNP*YuL(1,0)).imag();
	CuW_22r_LNP = (CuW_d_LNP*YdduL(1,1) + CuW_u_LNP*YucuuL(1,1) + CuW_0_LNP*YuL(1,1)).real();
	CuW_22i_LNP = (CuW_d_LNP*YdduL(1,1) + CuW_u_LNP*YucuuL(1,1) + CuW_0_LNP*YuL(1,1)).imag();
	CuW_23r_LNP = (CuW_d_LNP*YdduL(1,2) + CuW_u_LNP*YucuuL(1,2) + CuW_0_LNP*YuL(1,2)).real();
	CuW_23i_LNP = (CuW_d_LNP*YdduL(1,2) + CuW_u_LNP*YucuuL(1,2) + CuW_0_LNP*YuL(1,2)).imag();
	CuW_31r_LNP = (CuW_d_LNP*YdduL(2,0) + CuW_u_LNP*YucuuL(2,0) + CuW_0_LNP*YuL(2,0)).real();
	CuW_31i_LNP = (CuW_d_LNP*YdduL(2,0) + CuW_u_LNP*YucuuL(2,0) + CuW_0_LNP*YuL(2,0)).imag();
	CuW_32r_LNP = (CuW_d_LNP*YdduL(2,1) + CuW_u_LNP*YucuuL(2,1) + CuW_0_LNP*YuL(2,1)).real();
	CuW_32i_LNP = (CuW_d_LNP*YdduL(2,1) + CuW_u_LNP*YucuuL(2,1) + CuW_0_LNP*YuL(2,1)).imag();
	CuW_33r_LNP = (CuW_d_LNP*YdduL(2,2) + CuW_u_LNP*YucuuL(2,2) + CuW_0_LNP*YuL(2,2)).real();
	CuW_33i_LNP = (CuW_d_LNP*YdduL(2,2) + CuW_u_LNP*YucuuL(2,2) + CuW_0_LNP*YuL(2,2)).imag();

	CuB_11r_LNP = (CuB_d_LNP*YdduL(0,0) + CuB_u_LNP*YucuuL(0,0) + CuB_0_LNP*YuL(0,0)).real();
	CuB_11i_LNP = (CuB_d_LNP*YdduL(0,0) + CuB_u_LNP*YucuuL(0,0) + CuB_0_LNP*YuL(0,0)).imag();
	CuB_12r_LNP = (CuB_d_LNP*YdduL(0,1) + CuB_u_LNP*YucuuL(0,1) + CuB_0_LNP*YuL(0,1)).real();
	CuB_12i_LNP = (CuB_d_LNP*YdduL(0,1) + CuB_u_LNP*YucuuL(0,1) + CuB_0_LNP*YuL(0,1)).imag();
	CuB_13r_LNP = (CuB_d_LNP*YdduL(0,2) + CuB_u_LNP*YucuuL(0,2) + CuB_0_LNP*YuL(0,2)).real();
	CuB_13i_LNP = (CuB_d_LNP*YdduL(0,2) + CuB_u_LNP*YucuuL(0,2) + CuB_0_LNP*YuL(0,2)).imag();
	CuB_21r_LNP = (CuB_d_LNP*YdduL(1,0) + CuB_u_LNP*YucuuL(1,0) + CuB_0_LNP*YuL(1,0)).real();
	CuB_21i_LNP = (CuB_d_LNP*YdduL(1,0) + CuB_u_LNP*YucuuL(1,0) + CuB_0_LNP*YuL(1,0)).imag();
	CuB_22r_LNP = (CuB_d_LNP*YdduL(1,1) + CuB_u_LNP*YucuuL(1,1) + CuB_0_LNP*YuL(1,1)).real();
	CuB_22i_LNP = (CuB_d_LNP*YdduL(1,1) + CuB_u_LNP*YucuuL(1,1) + CuB_0_LNP*YuL(1,1)).imag();
	CuB_23r_LNP = (CuB_d_LNP*YdduL(1,2) + CuB_u_LNP*YucuuL(1,2) + CuB_0_LNP*YuL(1,2)).real();
	CuB_23i_LNP = (CuB_d_LNP*YdduL(1,2) + CuB_u_LNP*YucuuL(1,2) + CuB_0_LNP*YuL(1,2)).imag();
	CuB_31r_LNP = (CuB_d_LNP*YdduL(2,0) + CuB_u_LNP*YucuuL(2,0) + CuB_0_LNP*YuL(2,0)).real();
	CuB_31i_LNP = (CuB_d_LNP*YdduL(2,0) + CuB_u_LNP*YucuuL(2,0) + CuB_0_LNP*YuL(2,0)).imag();
	CuB_32r_LNP = (CuB_d_LNP*YdduL(2,1) + CuB_u_LNP*YucuuL(2,1) + CuB_0_LNP*YuL(2,1)).real();
	CuB_32i_LNP = (CuB_d_LNP*YdduL(2,1) + CuB_u_LNP*YucuuL(2,1) + CuB_0_LNP*YuL(2,1)).imag();
	CuB_33r_LNP = (CuB_d_LNP*YdduL(2,2) + CuB_u_LNP*YucuuL(2,2) + CuB_0_LNP*YuL(2,2)).real();
	CuB_33i_LNP = (CuB_d_LNP*YdduL(2,2) + CuB_u_LNP*YucuuL(2,2) + CuB_0_LNP*YuL(2,2)).imag();

	CdH_11r_LNP = (CdH_d_LNP*YdddL(0) + CdH_0_LNP*YdL(0) + CdH_u_LNP*YucudL(0,0)).real();
	CdH_11i_LNP = (CdH_d_LNP*YdddL(0) + CdH_0_LNP*YdL(0) + CdH_u_LNP*YucudL(0,0)).imag();
	CdH_12r_LNP = (CdH_u_LNP*YucudL(0,1)).real();
	CdH_12i_LNP = (CdH_u_LNP*YucudL(0,1)).imag();
	CdH_13r_LNP = (CdH_u_LNP*YucudL(0,2)).real();
	CdH_13i_LNP = (CdH_u_LNP*YucudL(0,2)).imag();
	CdH_21r_LNP = (CdH_u_LNP*YucudL(1,0)).real();
	CdH_21i_LNP = (CdH_u_LNP*YucudL(1,0)).imag();
	CdH_22r_LNP = (CdH_d_LNP*YdddL(1) + CdH_0_LNP*YdL(1) + CdH_u_LNP*YucudL(1,1)).real();
	CdH_22i_LNP = (CdH_d_LNP*YdddL(1) + CdH_0_LNP*YdL(1) + CdH_u_LNP*YucudL(1,1)).imag();
	CdH_23r_LNP = (CdH_u_LNP*YucudL(1,2)).real();
	CdH_23i_LNP = (CdH_u_LNP*YucudL(1,2)).imag();
	CdH_31r_LNP = (CdH_u_LNP*YucudL(2,0)).real();
	CdH_31i_LNP = (CdH_u_LNP*YucudL(2,0)).imag();
	CdH_32r_LNP = (CdH_u_LNP*YucudL(2,1)).real();
	CdH_32i_LNP = (CdH_u_LNP*YucudL(2,1)).imag();
	CdH_33r_LNP = (CdH_d_LNP*YdddL(2) + CdH_0_LNP*YdL(2) + CdH_u_LNP*YucudL(2,2)).real();
	CdH_33i_LNP = (CdH_d_LNP*YdddL(2) + CdH_0_LNP*YdL(2) + CdH_u_LNP*YucudL(2,2)).imag();

	CdG_11r_LNP = (CdG_d_LNP*YdddL(0) + CdG_0_LNP*YdL(0) + CdG_u_LNP*YucudL(0,0)).real();
	CdG_11i_LNP = (CdG_d_LNP*YdddL(0) + CdG_0_LNP*YdL(0) + CdG_u_LNP*YucudL(0,0)).imag();
	CdG_12r_LNP = (CdG_u_LNP*YucudL(0,1)).real();
	CdG_12i_LNP = (CdG_u_LNP*YucudL(0,1)).imag();
	CdG_13r_LNP = (CdG_u_LNP*YucudL(0,2)).real();
	CdG_13i_LNP = (CdG_u_LNP*YucudL(0,2)).imag();
	CdG_21r_LNP = (CdG_u_LNP*YucudL(1,0)).real();
	CdG_21i_LNP = (CdG_u_LNP*YucudL(1,0)).imag();
	CdG_22r_LNP = (CdG_d_LNP*YdddL(1) + CdG_0_LNP*YdL(1) + CdG_u_LNP*YucudL(1,1)).real();
	CdG_22i_LNP = (CdG_d_LNP*YdddL(1) + CdG_0_LNP*YdL(1) + CdG_u_LNP*YucudL(1,1)).imag();
	CdG_23r_LNP = (CdG_u_LNP*YucudL(1,2)).real();
	CdG_23i_LNP = (CdG_u_LNP*YucudL(1,2)).imag();
	CdG_31r_LNP = (CdG_u_LNP*YucudL(2,0)).real();
	CdG_31i_LNP = (CdG_u_LNP*YucudL(2,0)).imag();
	CdG_32r_LNP = (CdG_u_LNP*YucudL(2,1)).real();
	CdG_32i_LNP = (CdG_u_LNP*YucudL(2,1)).imag();
	CdG_33r_LNP = (CdG_d_LNP*YdddL(2) + CdG_0_LNP*YdL(2) + CdG_u_LNP*YucudL(2,2)).real();
	CdG_33i_LNP = (CdG_d_LNP*YdddL(2) + CdG_0_LNP*YdL(2) + CdG_u_LNP*YucudL(2,2)).imag();

	CdW_11r_LNP = (CdW_d_LNP*YdddL(0) + CdW_0_LNP*YdL(0) + CdW_u_LNP*YucudL(0,0)).real();
	CdW_11i_LNP = (CdW_d_LNP*YdddL(0) + CdW_0_LNP*YdL(0) + CdW_u_LNP*YucudL(0,0)).imag();
	CdW_12r_LNP = (CdW_u_LNP*YucudL(0,1)).real();
	CdW_12i_LNP = (CdW_u_LNP*YucudL(0,1)).imag();
	CdW_13r_LNP = (CdW_u_LNP*YucudL(0,2)).real();
	CdW_13i_LNP = (CdW_u_LNP*YucudL(0,2)).imag();
	CdW_21r_LNP = (CdW_u_LNP*YucudL(1,0)).real();
	CdW_21i_LNP = (CdW_u_LNP*YucudL(1,0)).imag();
	CdW_22r_LNP = (CdW_d_LNP*YdddL(1) + CdW_0_LNP*YdL(1) + CdW_u_LNP*YucudL(1,1)).real();
	CdW_22i_LNP = (CdW_d_LNP*YdddL(1) + CdW_0_LNP*YdL(1) + CdW_u_LNP*YucudL(1,1)).imag();
	CdW_23r_LNP = (CdW_u_LNP*YucudL(1,2)).real();
	CdW_23i_LNP = (CdW_u_LNP*YucudL(1,2)).imag();
	CdW_31r_LNP = (CdW_u_LNP*YucudL(2,0)).real();
	CdW_31i_LNP = (CdW_u_LNP*YucudL(2,0)).imag();
	CdW_32r_LNP = (CdW_u_LNP*YucudL(2,1)).real();
	CdW_32i_LNP = (CdW_u_LNP*YucudL(2,1)).imag();
	CdW_33r_LNP = (CdW_d_LNP*YdddL(2) + CdW_0_LNP*YdL(2) + CdW_u_LNP*YucudL(2,2)).real();
	CdW_33i_LNP = (CdW_d_LNP*YdddL(2) + CdW_0_LNP*YdL(2) + CdW_u_LNP*YucudL(2,2)).imag();

	CdB_11r_LNP = (CdB_d_LNP*YdddL(0) + CdB_0_LNP*YdL(0) + CdB_u_LNP*YucudL(0,0)).real();
	CdB_11i_LNP = (CdB_d_LNP*YdddL(0) + CdB_0_LNP*YdL(0) + CdB_u_LNP*YucudL(0,0)).imag();
	CdB_12r_LNP = (CdB_u_LNP*YucudL(0,1)).real();
	CdB_12i_LNP = (CdB_u_LNP*YucudL(0,1)).imag();
	CdB_13r_LNP = (CdB_u_LNP*YucudL(0,2)).real();
	CdB_13i_LNP = (CdB_u_LNP*YucudL(0,2)).imag();
	CdB_21r_LNP = (CdB_u_LNP*YucudL(1,0)).real();
	CdB_21i_LNP = (CdB_u_LNP*YucudL(1,0)).imag();
	CdB_22r_LNP = (CdB_d_LNP*YdddL(1) + CdB_0_LNP*YdL(1) + CdB_u_LNP*YucudL(1,1)).real();
	CdB_22i_LNP = (CdB_d_LNP*YdddL(1) + CdB_0_LNP*YdL(1) + CdB_u_LNP*YucudL(1,1)).imag();
	CdB_23r_LNP = (CdB_u_LNP*YucudL(1,2)).real();
	CdB_23i_LNP = (CdB_u_LNP*YucudL(1,2)).imag();
	CdB_31r_LNP = (CdB_u_LNP*YucudL(2,0)).real();
	CdB_31i_LNP = (CdB_u_LNP*YucudL(2,0)).imag();
	CdB_32r_LNP = (CdB_u_LNP*YucudL(2,1)).real();
	CdB_32i_LNP = (CdB_u_LNP*YucudL(2,1)).imag();
	CdB_33r_LNP = (CdB_d_LNP*YdddL(2) + CdB_0_LNP*YdL(2) + CdB_u_LNP*YucudL(2,2)).real();
	CdB_33i_LNP = (CdB_d_LNP*YdddL(2) + CdB_0_LNP*YdL(2) + CdB_u_LNP*YucudL(2,2)).imag();

	CHq1_11r_LNP = (CHq1_0_LNP + CHq1_d_LNP*YddL(0) + CHq1_u_LNP*YucuL(0,0)).real();
	CHq1_12r_LNP = (CHq1_u_LNP*YucuL(0,1)).real();
	CHq1_12i_LNP = (CHq1_u_LNP*YucuL(0,1)).imag();
	CHq1_13r_LNP = (CHq1_u_LNP*YucuL(0,2)).real();
	CHq1_13i_LNP = (CHq1_u_LNP*YucuL(0,2)).imag();
	CHq1_22r_LNP = (CHq1_0_LNP + CHq1_d_LNP*YddL(1) + CHq1_u_LNP*YucuL(1,1)).real();
	CHq1_23r_LNP = (CHq1_u_LNP*YucuL(1,2)).real();
	CHq1_23i_LNP = (CHq1_u_LNP*YucuL(1,2)).imag();
	CHq1_33r_LNP = (CHq1_0_LNP + CHq1_d_LNP*YddL(2) + CHq1_u_LNP*YucuL(2,2)).real();

	CHq3_11r_LNP = (CHq3_0_LNP + CHq3_d_LNP*YddL(0) + CHq3_u_LNP*YucuL(0,0)).real();
	CHq3_12r_LNP = (CHq3_u_LNP*YucuL(0,1)).real();
	CHq3_12i_LNP = (CHq3_u_LNP*YucuL(0,1)).imag();
	CHq3_13r_LNP = (CHq3_u_LNP*YucuL(0,2)).real();
	CHq3_13i_LNP = (CHq3_u_LNP*YucuL(0,2)).imag();
	CHq3_22r_LNP = (CHq3_0_LNP + CHq3_d_LNP*YddL(1) + CHq3_u_LNP*YucuL(1,1)).real();
	CHq3_23r_LNP = (CHq3_u_LNP*YucuL(1,2)).real();
	CHq3_23i_LNP = (CHq3_u_LNP*YucuL(1,2)).imag();
	CHq3_33r_LNP = (CHq3_0_LNP + CHq3_d_LNP*YddL(2) + CHq3_u_LNP*YucuL(2,2)).real();

	CHu_11r_LNP = (CHu_0_LNP + CHu_u_LNP*YuucL(0,0)).real();
	CHu_12r_LNP = (CHu_u_LNP*YuucL(0,1)).real();
	CHu_12i_LNP = (CHu_u_LNP*YuucL(0,1)).imag();
	CHu_13r_LNP = (CHu_u_LNP*YuucL(0,2)).real();
	CHu_13i_LNP = (CHu_u_LNP*YuucL(0,2)).imag();
	CHu_22r_LNP = (CHu_0_LNP + CHu_u_LNP*YuucL(1,1)).real();
	CHu_23r_LNP = (CHu_u_LNP*YuucL(1,2)).real();
	CHu_23i_LNP = (CHu_u_LNP*YuucL(1,2)).imag();
	CHu_33r_LNP = (CHu_0_LNP + CHu_u_LNP*YuucL(2,2)).real();

	CHd_11r_LNP = (CHd_0_LNP + CHd_d_LNP*YddL(0)).real();
	CHd_22r_LNP = (CHd_0_LNP + CHd_d_LNP*YddL(1)).real();
	CHd_33r_LNP = (CHd_0_LNP + CHd_d_LNP*YddL(2)).real();

	CHud_11r_LNP = (CHud_ud_LNP*YudL(0,0)).real();
	CHud_11i_LNP = (CHud_ud_LNP*YudL(0,0)).imag();
	CHud_12r_LNP = (CHud_ud_LNP*YudL(0,1)).real();
	CHud_12i_LNP = (CHud_ud_LNP*YudL(0,1)).imag();
	CHud_13r_LNP = (CHud_ud_LNP*YudL(0,2)).real();
	CHud_13i_LNP = (CHud_ud_LNP*YudL(0,2)).imag();
	CHud_21r_LNP = (CHud_ud_LNP*YudL(1,0)).real();
	CHud_21i_LNP = (CHud_ud_LNP*YudL(1,0)).imag();
	CHud_22r_LNP = (CHud_ud_LNP*YudL(1,1)).real();
	CHud_22i_LNP = (CHud_ud_LNP*YudL(1,1)).imag();
	CHud_23r_LNP = (CHud_ud_LNP*YudL(1,2)).real();
	CHud_23i_LNP = (CHud_ud_LNP*YudL(1,2)).imag();
	CHud_31r_LNP = (CHud_ud_LNP*YudL(2,0)).real();
	CHud_31i_LNP = (CHud_ud_LNP*YudL(2,0)).imag();
	CHud_32r_LNP = (CHud_ud_LNP*YudL(2,1)).real();
	CHud_32i_LNP = (CHud_ud_LNP*YudL(2,1)).imag();
	CHud_33r_LNP = (CHud_ud_LNP*YudL(2,2)).real();
	CHud_33i_LNP = (CHud_ud_LNP*YudL(2,2)).imag();

	Clq1_1111r_LNP = (Clq1_0_LNP + Clq1_d_LNP*YddL(0) + Clq1_u_LNP*YucuL(0,0)).real();
	Clq1_1112r_LNP = (Clq1_u_LNP*YucuL(0,1)).real();
	Clq1_1112i_LNP = (Clq1_u_LNP*YucuL(0,1)).imag();
	Clq1_1113r_LNP = (Clq1_u_LNP*YucuL(0,2)).real();
	Clq1_1113i_LNP = (Clq1_u_LNP*YucuL(0,2)).imag();
	Clq1_1122r_LNP = (Clq1_0_LNP + Clq1_d_LNP*YddL(1) + Clq1_u_LNP*YucuL(1,1)).real();
	Clq1_1123r_LNP = (Clq1_u_LNP*YucuL(1,2)).real();
	Clq1_1123i_LNP = (Clq1_u_LNP*YucuL(1,2)).imag();
	Clq1_1133r_LNP = (Clq1_0_LNP + Clq1_d_LNP*YddL(2) + Clq1_u_LNP*YucuL(2,2)).real();
	Clq1_2211r_LNP = (Clq1_0_LNP + Clq1_d_LNP*YddL(0) + Clq1_u_LNP*YucuL(0,0)).real();
	Clq1_2212r_LNP = (Clq1_u_LNP*YucuL(0,1)).real();
	Clq1_2212i_LNP = (Clq1_u_LNP*YucuL(0,1)).imag();
	Clq1_2213r_LNP = (Clq1_u_LNP*YucuL(0,2)).real();
	Clq1_2213i_LNP = (Clq1_u_LNP*YucuL(0,2)).imag();
	Clq1_2222r_LNP = (Clq1_0_LNP + Clq1_d_LNP*YddL(1) + Clq1_u_LNP*YucuL(1,1)).real();
	Clq1_2223r_LNP = (Clq1_u_LNP*YucuL(1,2)).real();
	Clq1_2223i_LNP = (Clq1_u_LNP*YucuL(1,2)).imag();
	Clq1_2233r_LNP = (Clq1_0_LNP + Clq1_d_LNP*YddL(2) + Clq1_u_LNP*YucuL(2,2)).real();
	Clq1_3311r_LNP = (Clq1_0_LNP + Clq1_d_LNP*YddL(0) + Clq1_u_LNP*YucuL(0,0)).real();
	Clq1_3312r_LNP = (Clq1_u_LNP*YucuL(0,1)).real();
	Clq1_3312i_LNP = (Clq1_u_LNP*YucuL(0,1)).imag();
	Clq1_3313r_LNP = (Clq1_u_LNP*YucuL(0,2)).real();
	Clq1_3313i_LNP = (Clq1_u_LNP*YucuL(0,2)).imag();
	Clq1_3322r_LNP = (Clq1_0_LNP + Clq1_d_LNP*YddL(1) + Clq1_u_LNP*YucuL(1,1)).real();
	Clq1_3323r_LNP = (Clq1_u_LNP*YucuL(1,2)).real();
	Clq1_3323i_LNP = (Clq1_u_LNP*YucuL(1,2)).imag();
	Clq1_3333r_LNP = (Clq1_0_LNP + Clq1_d_LNP*YddL(2) + Clq1_u_LNP*YucuL(2,2)).real();

	Clq3_1111r_LNP = (Clq3_0_LNP + Clq3_d_LNP*YddL(0) + Clq3_u_LNP*YucuL(0,0)).real();
	Clq3_1112r_LNP = (Clq3_u_LNP*YucuL(0,1)).real();
	Clq3_1112i_LNP = (Clq3_u_LNP*YucuL(0,1)).imag();
	Clq3_1113r_LNP = (Clq3_u_LNP*YucuL(0,2)).real();
	Clq3_1113i_LNP = (Clq3_u_LNP*YucuL(0,2)).imag();
	Clq3_1122r_LNP = (Clq3_0_LNP + Clq3_d_LNP*YddL(1) + Clq3_u_LNP*YucuL(1,1)).real();
	Clq3_1123r_LNP = (Clq3_u_LNP*YucuL(1,2)).real();
	Clq3_1123i_LNP = (Clq3_u_LNP*YucuL(1,2)).imag();
	Clq3_1133r_LNP = (Clq3_0_LNP + Clq3_d_LNP*YddL(2) + Clq3_u_LNP*YucuL(2,2)).real();
	Clq3_2211r_LNP = (Clq3_0_LNP + Clq3_d_LNP*YddL(0) + Clq3_u_LNP*YucuL(0,0)).real();
	Clq3_2212r_LNP = (Clq3_u_LNP*YucuL(0,1)).real();
	Clq3_2212i_LNP = (Clq3_u_LNP*YucuL(0,1)).imag();
	Clq3_2213r_LNP = (Clq3_u_LNP*YucuL(0,2)).real();
	Clq3_2213i_LNP = (Clq3_u_LNP*YucuL(0,2)).imag();
	Clq3_2222r_LNP = (Clq3_0_LNP + Clq3_d_LNP*YddL(1) + Clq3_u_LNP*YucuL(1,1)).real();
	Clq3_2223r_LNP = (Clq3_u_LNP*YucuL(1,2)).real();
	Clq3_2223i_LNP = (Clq3_u_LNP*YucuL(1,2)).imag();
	Clq3_2233r_LNP = (Clq3_0_LNP + Clq3_d_LNP*YddL(2) + Clq3_u_LNP*YucuL(2,2)).real();
	Clq3_3311r_LNP = (Clq3_0_LNP + Clq3_d_LNP*YddL(0) + Clq3_u_LNP*YucuL(0,0)).real();
	Clq3_3312r_LNP = (Clq3_u_LNP*YucuL(0,1)).real();
	Clq3_3312i_LNP = (Clq3_u_LNP*YucuL(0,1)).imag();
	Clq3_3313r_LNP = (Clq3_u_LNP*YucuL(0,2)).real();
	Clq3_3313i_LNP = (Clq3_u_LNP*YucuL(0,2)).imag();
	Clq3_3322r_LNP = (Clq3_0_LNP + Clq3_d_LNP*YddL(1) + Clq3_u_LNP*YucuL(1,1)).real();
	Clq3_3323r_LNP = (Clq3_u_LNP*YucuL(1,2)).real();
	Clq3_3323i_LNP = (Clq3_u_LNP*YucuL(1,2)).imag();
	Clq3_3333r_LNP = (Clq3_0_LNP + Clq3_d_LNP*YddL(2) + Clq3_u_LNP*YucuL(2,2)).real();

	Cqe_1111r_LNP = (Cqe_0_LNP + Cqe_d_LNP*YddL(0) + Cqe_u_LNP*YucuL(0,0)).real();
	Cqe_1122r_LNP = (Cqe_0_LNP + Cqe_d_LNP*YddL(0) + Cqe_u_LNP*YucuL(0,0)).real();
	Cqe_1133r_LNP = (Cqe_0_LNP + Cqe_d_LNP*YddL(0) + Cqe_u_LNP*YucuL(0,0)).real();
	Cqe_1211r_LNP = (Cqe_u_LNP*YucuL(0,1)).real();
	Cqe_1211i_LNP = (Cqe_u_LNP*YucuL(0,1)).imag();
	Cqe_1222r_LNP = (Cqe_u_LNP*YucuL(0,1)).real();
	Cqe_1222i_LNP = (Cqe_u_LNP*YucuL(0,1)).imag();
	Cqe_1233r_LNP = (Cqe_u_LNP*YucuL(0,1)).real();
	Cqe_1233i_LNP = (Cqe_u_LNP*YucuL(0,1)).imag();
	Cqe_1311r_LNP = (Cqe_u_LNP*YucuL(0,2)).real();
	Cqe_1311i_LNP = (Cqe_u_LNP*YucuL(0,2)).imag();
	Cqe_1322r_LNP = (Cqe_u_LNP*YucuL(0,2)).real();
	Cqe_1322i_LNP = (Cqe_u_LNP*YucuL(0,2)).imag();
	Cqe_1333r_LNP = (Cqe_u_LNP*YucuL(0,2)).real();
	Cqe_1333i_LNP = (Cqe_u_LNP*YucuL(0,2)).imag();
	Cqe_2211r_LNP = (Cqe_0_LNP + Cqe_d_LNP*YddL(1) + Cqe_u_LNP*YucuL(1,1)).real();
	Cqe_2222r_LNP = (Cqe_0_LNP + Cqe_d_LNP*YddL(1) + Cqe_u_LNP*YucuL(1,1)).real();
	Cqe_2233r_LNP = (Cqe_0_LNP + Cqe_d_LNP*YddL(1) + Cqe_u_LNP*YucuL(1,1)).real();
	Cqe_2311r_LNP = (Cqe_u_LNP*YucuL(1,2)).real();
	Cqe_2311i_LNP = (Cqe_u_LNP*YucuL(1,2)).imag();
	Cqe_2322r_LNP = (Cqe_u_LNP*YucuL(1,2)).real();
	Cqe_2322i_LNP = (Cqe_u_LNP*YucuL(1,2)).imag();
	Cqe_2333r_LNP = (Cqe_u_LNP*YucuL(1,2)).real();
	Cqe_2333i_LNP = (Cqe_u_LNP*YucuL(1,2)).imag();
	Cqe_3311r_LNP = (Cqe_0_LNP + Cqe_d_LNP*YddL(2) + Cqe_u_LNP*YucuL(2,2)).real();
	Cqe_3322r_LNP = (Cqe_0_LNP + Cqe_d_LNP*YddL(2) + Cqe_u_LNP*YucuL(2,2)).real();
	Cqe_3333r_LNP = (Cqe_0_LNP + Cqe_d_LNP*YddL(2) + Cqe_u_LNP*YucuL(2,2)).real();

	Clu_1111r_LNP = (Clu_0_LNP + Clu_u_LNP*YuucL(0,0)).real();
	Clu_1112r_LNP = (Clu_u_LNP*YuucL(0,1)).real();
	Clu_1112i_LNP = (Clu_u_LNP*YuucL(0,1)).imag();
	Clu_1113r_LNP = (Clu_u_LNP*YuucL(0,2)).real();
	Clu_1113i_LNP = (Clu_u_LNP*YuucL(0,2)).imag();
	Clu_1122r_LNP = (Clu_0_LNP + Clu_u_LNP*YuucL(1,1)).real();
	Clu_1123r_LNP = (Clu_u_LNP*YuucL(1,2)).real();
	Clu_1123i_LNP = (Clu_u_LNP*YuucL(1,2)).imag();
	Clu_1133r_LNP = (Clu_0_LNP + Clu_u_LNP*YuucL(2,2)).real();
	Clu_2211r_LNP = (Clu_0_LNP + Clu_u_LNP*YuucL(0,0)).real();
	Clu_2212r_LNP = (Clu_u_LNP*YuucL(0,1)).real();
	Clu_2212i_LNP = (Clu_u_LNP*YuucL(0,1)).imag();
	Clu_2213r_LNP = (Clu_u_LNP*YuucL(0,2)).real();
	Clu_2213i_LNP = (Clu_u_LNP*YuucL(0,2)).imag();
	Clu_2222r_LNP = (Clu_0_LNP + Clu_u_LNP*YuucL(1,1)).real();
	Clu_2223r_LNP = (Clu_u_LNP*YuucL(1,2)).real();
	Clu_2223i_LNP = (Clu_u_LNP*YuucL(1,2)).imag();
	Clu_2233r_LNP = (Clu_0_LNP + Clu_u_LNP*YuucL(2,2)).real();
	Clu_3311r_LNP = (Clu_0_LNP + Clu_u_LNP*YuucL(0,0)).real();
	Clu_3312r_LNP = (Clu_u_LNP*YuucL(0,1)).real();
	Clu_3312i_LNP = (Clu_u_LNP*YuucL(0,1)).imag();
	Clu_3313r_LNP = (Clu_u_LNP*YuucL(0,2)).real();
	Clu_3313i_LNP = (Clu_u_LNP*YuucL(0,2)).imag();
	Clu_3322r_LNP = (Clu_0_LNP + Clu_u_LNP*YuucL(1,1)).real();
	Clu_3323r_LNP = (Clu_u_LNP*YuucL(1,2)).real();
	Clu_3323i_LNP = (Clu_u_LNP*YuucL(1,2)).imag();
	Clu_3333r_LNP = (Clu_0_LNP + Clu_u_LNP*YuucL(2,2)).real();

	Ceu_1111r_LNP = (Ceu_0_LNP + Ceu_u_LNP*YuucL(0,0)).real();
	Ceu_1112r_LNP = (Ceu_u_LNP*YuucL(0,1)).real();
	Ceu_1112i_LNP = (Ceu_u_LNP*YuucL(0,1)).imag();
	Ceu_1113r_LNP = (Ceu_u_LNP*YuucL(0,2)).real();
	Ceu_1113i_LNP = (Ceu_u_LNP*YuucL(0,2)).imag();
	Ceu_1122r_LNP = (Ceu_0_LNP + Ceu_u_LNP*YuucL(1,1)).real();
	Ceu_1123r_LNP = (Ceu_u_LNP*YuucL(1,2)).real();
	Ceu_1123i_LNP = (Ceu_u_LNP*YuucL(1,2)).imag();
	Ceu_1133r_LNP = (Ceu_0_LNP + Ceu_u_LNP*YuucL(2,2)).real();
	Ceu_2211r_LNP = (Ceu_0_LNP + Ceu_u_LNP*YuucL(0,0)).real();
	Ceu_2212r_LNP = (Ceu_u_LNP*YuucL(0,1)).real();
	Ceu_2212i_LNP = (Ceu_u_LNP*YuucL(0,1)).imag();
	Ceu_2213r_LNP = (Ceu_u_LNP*YuucL(0,2)).real();
	Ceu_2213i_LNP = (Ceu_u_LNP*YuucL(0,2)).imag();
	Ceu_2222r_LNP = (Ceu_0_LNP + Ceu_u_LNP*YuucL(1,1)).real();
	Ceu_2223r_LNP = (Ceu_u_LNP*YuucL(1,2)).real();
	Ceu_2223i_LNP = (Ceu_u_LNP*YuucL(1,2)).imag();
	Ceu_2233r_LNP = (Ceu_0_LNP + Ceu_u_LNP*YuucL(2,2)).real();
	Ceu_3311r_LNP = (Ceu_0_LNP + Ceu_u_LNP*YuucL(0,0)).real();
	Ceu_3312r_LNP = (Ceu_u_LNP*YuucL(0,1)).real();
	Ceu_3312i_LNP = (Ceu_u_LNP*YuucL(0,1)).imag();
	Ceu_3313r_LNP = (Ceu_u_LNP*YuucL(0,2)).real();
	Ceu_3313i_LNP = (Ceu_u_LNP*YuucL(0,2)).imag();
	Ceu_3322r_LNP = (Ceu_0_LNP + Ceu_u_LNP*YuucL(1,1)).real();
	Ceu_3323r_LNP = (Ceu_u_LNP*YuucL(1,2)).real();
	Ceu_3323i_LNP = (Ceu_u_LNP*YuucL(1,2)).imag();
	Ceu_3333r_LNP = (Ceu_0_LNP + Ceu_u_LNP*YuucL(2,2)).real();

	Cld_1111r_LNP = (Cld_0_LNP + Cld_d_LNP*YddL(0)).real();
	Cld_1122r_LNP = (Cld_0_LNP + Cld_d_LNP*YddL(1)).real();
	Cld_1133r_LNP = (Cld_0_LNP + Cld_d_LNP*YddL(2)).real();
	Cld_2211r_LNP = (Cld_0_LNP + Cld_d_LNP*YddL(0)).real();
	Cld_2222r_LNP = (Cld_0_LNP + Cld_d_LNP*YddL(1)).real();
	Cld_2233r_LNP = (Cld_0_LNP + Cld_d_LNP*YddL(2)).real();
	Cld_3311r_LNP = (Cld_0_LNP + Cld_d_LNP*YddL(0)).real();
	Cld_3322r_LNP = (Cld_0_LNP + Cld_d_LNP*YddL(1)).real();
	Cld_3333r_LNP = (Cld_0_LNP + Cld_d_LNP*YddL(2)).real();

	Ced_1111r_LNP = (Ced_0_LNP + Ced_d_LNP*YddL(0)).real();
	Ced_1122r_LNP = (Ced_0_LNP + Ced_d_LNP*YddL(1)).real();
	Ced_1133r_LNP = (Ced_0_LNP + Ced_d_LNP*YddL(2)).real();
	Ced_2211r_LNP = (Ced_0_LNP + Ced_d_LNP*YddL(0)).real();
	Ced_2222r_LNP = (Ced_0_LNP + Ced_d_LNP*YddL(1)).real();
	Ced_2233r_LNP = (Ced_0_LNP + Ced_d_LNP*YddL(2)).real();
	Ced_3311r_LNP = (Ced_0_LNP + Ced_d_LNP*YddL(0)).real();
	Ced_3322r_LNP = (Ced_0_LNP + Ced_d_LNP*YddL(1)).real();
	Ced_3333r_LNP = (Ced_0_LNP + Ced_d_LNP*YddL(2)).real();

	Cqq1_1111r_LNP = (2*Cqq1_00_LNP + 2*Cqq1_0d_LNP*YddL(0) + 2*Cqq1_d0_LNP*YddL(0) + 2*Cqq1_dd_LNP*YddL(0)*YddL(0) + 2*Cqq1_0u_LNP*YucuL(0,0) + 2*Cqq1_u0_LNP*YucuL(0,0) + 2*Cqq1_du_LNP*YddL(0)*YucuL(0,0) + 2*Cqq1_ud_LNP*YddL(0)*YucuL(0,0) + 2*Cqq1_uu_LNP*YucuL(0,0)*YucuL(0,0)).real();
	Cqq1_1112r_LNP = (Cqq1_0u_LNP*YucuL(0,1) + Cqq1_u0_LNP*YucuL(0,1) + Cqq1_du_LNP*YddL(0)*YucuL(0,1) + Cqq1_ud_LNP*YddL(0)*YucuL(0,1) + 2*Cqq1_uu_LNP*YucuL(0,0)*YucuL(0,1)).real();
	Cqq1_1112i_LNP = (Cqq1_0u_LNP*YucuL(0,1) + Cqq1_u0_LNP*YucuL(0,1) + Cqq1_du_LNP*YddL(0)*YucuL(0,1) + Cqq1_ud_LNP*YddL(0)*YucuL(0,1) + 2*Cqq1_uu_LNP*YucuL(0,0)*YucuL(0,1)).imag();
	Cqq1_1113r_LNP = (Cqq1_0u_LNP*YucuL(0,2) + Cqq1_u0_LNP*YucuL(0,2) + Cqq1_du_LNP*YddL(0)*YucuL(0,2) + Cqq1_ud_LNP*YddL(0)*YucuL(0,2) + 2*Cqq1_uu_LNP*YucuL(0,0)*YucuL(0,2)).real();
	Cqq1_1113i_LNP = (Cqq1_0u_LNP*YucuL(0,2) + Cqq1_u0_LNP*YucuL(0,2) + Cqq1_du_LNP*YddL(0)*YucuL(0,2) + Cqq1_ud_LNP*YddL(0)*YucuL(0,2) + 2*Cqq1_uu_LNP*YucuL(0,0)*YucuL(0,2)).imag();
	Cqq1_1122r_LNP = (Cqq1_00_LNP + Cqq1_d0_LNP*YddL(0) + Cqq1_0d_LNP*YddL(1) + Cqq1_dd_LNP*YddL(0)*YddL(1) + Cqq1_u0_LNP*YucuL(0,0) + Cqq1_ud_LNP*YddL(1)*YucuL(0,0) + Cqq1_uu_LNP*YucuL(0,1)*YucuL(1,0) + Cqq1_0u_LNP*YucuL(1,1) + Cqq1_du_LNP*YddL(0)*YucuL(1,1) + Cqq1_uu_LNP*YucuL(0,0)*YucuL(1,1)).real();
	Cqq1_1123r_LNP = (Cqq1_uu_LNP*YucuL(0,2)*YucuL(1,0) + Cqq1_0u_LNP*YucuL(1,2) + Cqq1_du_LNP*YddL(0)*YucuL(1,2) + Cqq1_uu_LNP*YucuL(0,0)*YucuL(1,2)).real();
	Cqq1_1123i_LNP = (Cqq1_uu_LNP*YucuL(0,2)*YucuL(1,0) + Cqq1_0u_LNP*YucuL(1,2) + Cqq1_du_LNP*YddL(0)*YucuL(1,2) + Cqq1_uu_LNP*YucuL(0,0)*YucuL(1,2)).imag();
	Cqq1_1133r_LNP = (Cqq1_00_LNP + Cqq1_d0_LNP*YddL(0) + Cqq1_0d_LNP*YddL(2) + Cqq1_dd_LNP*YddL(0)*YddL(2) + Cqq1_u0_LNP*YucuL(0,0) + Cqq1_ud_LNP*YddL(2)*YucuL(0,0) + Cqq1_uu_LNP*YucuL(0,2)*YucuL(2,0) + Cqq1_0u_LNP*YucuL(2,2) + Cqq1_du_LNP*YddL(0)*YucuL(2,2) + Cqq1_uu_LNP*YucuL(0,0)*YucuL(2,2)).real();
	Cqq1_1212r_LNP = (2*Cqq1_uu_LNP*YucuL(0,1)*YucuL(0,1)).real();
	Cqq1_1212i_LNP = (2*Cqq1_uu_LNP*YucuL(0,1)*YucuL(0,1)).imag();
	Cqq1_1213r_LNP = (2*Cqq1_uu_LNP*YucuL(0,1)*YucuL(0,2)).real();
	Cqq1_1213i_LNP = (2*Cqq1_uu_LNP*YucuL(0,1)*YucuL(0,2)).imag();
	Cqq1_1221r_LNP = (Cqq1_00_LNP + Cqq1_d0_LNP*YddL(0) + Cqq1_0d_LNP*YddL(1) + Cqq1_dd_LNP*YddL(0)*YddL(1) + Cqq1_u0_LNP*YucuL(0,0) + Cqq1_ud_LNP*YddL(1)*YucuL(0,0) + Cqq1_uu_LNP*YucuL(0,1)*YucuL(1,0) + Cqq1_0u_LNP*YucuL(1,1) + Cqq1_du_LNP*YddL(0)*YucuL(1,1) + Cqq1_uu_LNP*YucuL(0,0)*YucuL(1,1)).real();
	Cqq1_1222r_LNP = (2*Cqq1_u0_LNP*YucuL(0,1) + 2*Cqq1_ud_LNP*YddL(1)*YucuL(0,1) + 2*Cqq1_uu_LNP*YucuL(0,1)*YucuL(1,1)).real();
	Cqq1_1222i_LNP = (2*Cqq1_u0_LNP*YucuL(0,1) + 2*Cqq1_ud_LNP*YddL(1)*YucuL(0,1) + 2*Cqq1_uu_LNP*YucuL(0,1)*YucuL(1,1)).imag();
	Cqq1_1223r_LNP = (Cqq1_u0_LNP*YucuL(0,2) + Cqq1_ud_LNP*YddL(1)*YucuL(0,2) + Cqq1_uu_LNP*YucuL(0,2)*YucuL(1,1) + Cqq1_uu_LNP*YucuL(0,1)*YucuL(1,2)).real();
	Cqq1_1223i_LNP = (Cqq1_u0_LNP*YucuL(0,2) + Cqq1_ud_LNP*YddL(1)*YucuL(0,2) + Cqq1_uu_LNP*YucuL(0,2)*YucuL(1,1) + Cqq1_uu_LNP*YucuL(0,1)*YucuL(1,2)).imag();
	Cqq1_1231r_LNP = (Cqq1_uu_LNP*YucuL(0,1)*YucuL(2,0) + Cqq1_0u_LNP*YucuL(2,1) + Cqq1_du_LNP*YddL(0)*YucuL(2,1) + Cqq1_uu_LNP*YucuL(0,0)*YucuL(2,1)).real();
	Cqq1_1231i_LNP = (Cqq1_uu_LNP*YucuL(0,1)*YucuL(2,0) + Cqq1_0u_LNP*YucuL(2,1) + Cqq1_du_LNP*YddL(0)*YucuL(2,1) + Cqq1_uu_LNP*YucuL(0,0)*YucuL(2,1)).imag();
	Cqq1_1232r_LNP = (2*Cqq1_uu_LNP*YucuL(0,1)*YucuL(2,1)).real();
	Cqq1_1232i_LNP = (2*Cqq1_uu_LNP*YucuL(0,1)*YucuL(2,1)).imag();
	Cqq1_1233r_LNP = (Cqq1_u0_LNP*YucuL(0,1) + Cqq1_ud_LNP*YddL(2)*YucuL(0,1) + Cqq1_uu_LNP*YucuL(0,2)*YucuL(2,1) + Cqq1_uu_LNP*YucuL(0,1)*YucuL(2,2)).real();
	Cqq1_1233i_LNP = (Cqq1_u0_LNP*YucuL(0,1) + Cqq1_ud_LNP*YddL(2)*YucuL(0,1) + Cqq1_uu_LNP*YucuL(0,2)*YucuL(2,1) + Cqq1_uu_LNP*YucuL(0,1)*YucuL(2,2)).imag();
	Cqq1_1313r_LNP = (2*Cqq1_uu_LNP*YucuL(0,2)*YucuL(0,2)).real();
	Cqq1_1313i_LNP = (2*Cqq1_uu_LNP*YucuL(0,2)*YucuL(0,2)).imag();
	Cqq1_1322r_LNP = (Cqq1_u0_LNP*YucuL(0,2) + Cqq1_ud_LNP*YddL(1)*YucuL(0,2) + Cqq1_uu_LNP*YucuL(0,2)*YucuL(1,1) + Cqq1_uu_LNP*YucuL(0,1)*YucuL(1,2)).real();
	Cqq1_1322i_LNP = (Cqq1_u0_LNP*YucuL(0,2) + Cqq1_ud_LNP*YddL(1)*YucuL(0,2) + Cqq1_uu_LNP*YucuL(0,2)*YucuL(1,1) + Cqq1_uu_LNP*YucuL(0,1)*YucuL(1,2)).imag();
	Cqq1_1323r_LNP = (2*Cqq1_uu_LNP*YucuL(0,2)*YucuL(1,2)).real();
	Cqq1_1323i_LNP = (2*Cqq1_uu_LNP*YucuL(0,2)*YucuL(1,2)).imag();
	Cqq1_1331r_LNP = (Cqq1_00_LNP + Cqq1_d0_LNP*YddL(0) + Cqq1_0d_LNP*YddL(2) + Cqq1_dd_LNP*YddL(0)*YddL(2) + Cqq1_u0_LNP*YucuL(0,0) + Cqq1_ud_LNP*YddL(2)*YucuL(0,0) + Cqq1_uu_LNP*YucuL(0,2)*YucuL(2,0) + Cqq1_0u_LNP*YucuL(2,2) + Cqq1_du_LNP*YddL(0)*YucuL(2,2) + Cqq1_uu_LNP*YucuL(0,0)*YucuL(2,2)).real();
	Cqq1_1332r_LNP = (Cqq1_u0_LNP*YucuL(0,1) + Cqq1_ud_LNP*YddL(2)*YucuL(0,1) + Cqq1_uu_LNP*YucuL(0,2)*YucuL(2,1) + Cqq1_uu_LNP*YucuL(0,1)*YucuL(2,2)).real();
	Cqq1_1332i_LNP = (Cqq1_u0_LNP*YucuL(0,1) + Cqq1_ud_LNP*YddL(2)*YucuL(0,1) + Cqq1_uu_LNP*YucuL(0,2)*YucuL(2,1) + Cqq1_uu_LNP*YucuL(0,1)*YucuL(2,2)).imag();
	Cqq1_1333r_LNP = (2*Cqq1_u0_LNP*YucuL(0,2) + 2*Cqq1_ud_LNP*YddL(2)*YucuL(0,2) + 2*Cqq1_uu_LNP*YucuL(0,2)*YucuL(2,2)).real();
	Cqq1_1333i_LNP = (2*Cqq1_u0_LNP*YucuL(0,2) + 2*Cqq1_ud_LNP*YddL(2)*YucuL(0,2) + 2*Cqq1_uu_LNP*YucuL(0,2)*YucuL(2,2)).imag();
	Cqq1_2222r_LNP = (2*Cqq1_00_LNP + 2*Cqq1_0d_LNP*YddL(1) + 2*Cqq1_d0_LNP*YddL(1) + 2*Cqq1_dd_LNP*YddL(1)*YddL(1) + 2*Cqq1_0u_LNP*YucuL(1,1) + 2*Cqq1_u0_LNP*YucuL(1,1) + 2*Cqq1_du_LNP*YddL(1)*YucuL(1,1) + 2*Cqq1_ud_LNP*YddL(1)*YucuL(1,1) + 2*Cqq1_uu_LNP*YucuL(1,1)*YucuL(1,1)).real();
	Cqq1_2223r_LNP = (Cqq1_0u_LNP*YucuL(1,2) + Cqq1_u0_LNP*YucuL(1,2) + Cqq1_du_LNP*YddL(1)*YucuL(1,2) + Cqq1_ud_LNP*YddL(1)*YucuL(1,2) + 2*Cqq1_uu_LNP*YucuL(1,1)*YucuL(1,2)).real();
	Cqq1_2223i_LNP = (Cqq1_0u_LNP*YucuL(1,2) + Cqq1_u0_LNP*YucuL(1,2) + Cqq1_du_LNP*YddL(1)*YucuL(1,2) + Cqq1_ud_LNP*YddL(1)*YucuL(1,2) + 2*Cqq1_uu_LNP*YucuL(1,1)*YucuL(1,2)).imag();
	Cqq1_2233r_LNP = (Cqq1_00_LNP + Cqq1_d0_LNP*YddL(1) + Cqq1_0d_LNP*YddL(2) + Cqq1_dd_LNP*YddL(1)*YddL(2) + Cqq1_u0_LNP*YucuL(1,1) + Cqq1_ud_LNP*YddL(2)*YucuL(1,1) + Cqq1_uu_LNP*YucuL(1,2)*YucuL(2,1) + Cqq1_0u_LNP*YucuL(2,2) + Cqq1_du_LNP*YddL(1)*YucuL(2,2) + Cqq1_uu_LNP*YucuL(1,1)*YucuL(2,2)).real();
	Cqq1_2323r_LNP = (2*Cqq1_uu_LNP*YucuL(1,2)*YucuL(1,2)).real();
	Cqq1_2323i_LNP = (2*Cqq1_uu_LNP*YucuL(1,2)*YucuL(1,2)).imag();
	Cqq1_2332r_LNP = (Cqq1_00_LNP + Cqq1_d0_LNP*YddL(1) + Cqq1_0d_LNP*YddL(2) + Cqq1_dd_LNP*YddL(1)*YddL(2) + Cqq1_u0_LNP*YucuL(1,1) + Cqq1_ud_LNP*YddL(2)*YucuL(1,1) + Cqq1_uu_LNP*YucuL(1,2)*YucuL(2,1) + Cqq1_0u_LNP*YucuL(2,2) + Cqq1_du_LNP*YddL(1)*YucuL(2,2) + Cqq1_uu_LNP*YucuL(1,1)*YucuL(2,2)).real();
	Cqq1_2333r_LNP = (2*Cqq1_u0_LNP*YucuL(1,2) + 2*Cqq1_ud_LNP*YddL(2)*YucuL(1,2) + 2*Cqq1_uu_LNP*YucuL(1,2)*YucuL(2,2)).real();
	Cqq1_2333i_LNP = (2*Cqq1_u0_LNP*YucuL(1,2) + 2*Cqq1_ud_LNP*YddL(2)*YucuL(1,2) + 2*Cqq1_uu_LNP*YucuL(1,2)*YucuL(2,2)).imag();
	Cqq1_3333r_LNP = (2*Cqq1_00_LNP + 2*Cqq1_0d_LNP*YddL(2) + 2*Cqq1_d0_LNP*YddL(2) + 2*Cqq1_dd_LNP*YddL(2)*YddL(2) + 2*Cqq1_0u_LNP*YucuL(2,2) + 2*Cqq1_u0_LNP*YucuL(2,2) + 2*Cqq1_du_LNP*YddL(2)*YucuL(2,2) + 2*Cqq1_ud_LNP*YddL(2)*YucuL(2,2) + 2*Cqq1_uu_LNP*YucuL(2,2)*YucuL(2,2)).real();

	Cqq3_1111r_LNP = (2*Cqq3_00_LNP + 2*Cqq3_0d_LNP*YddL(0) + 2*Cqq3_d0_LNP*YddL(0) + 2*Cqq3_dd_LNP*YddL(0)*YddL(0) + 2*Cqq3_0u_LNP*YucuL(0,0) + 2*Cqq3_u0_LNP*YucuL(0,0) + 2*Cqq3_du_LNP*YddL(0)*YucuL(0,0) + 2*Cqq3_ud_LNP*YddL(0)*YucuL(0,0) + 2*Cqq3_uu_LNP*YucuL(0,0)*YucuL(0,0)).real();
	Cqq3_1112r_LNP = (Cqq3_0u_LNP*YucuL(0,1) + Cqq3_u0_LNP*YucuL(0,1) + Cqq3_du_LNP*YddL(0)*YucuL(0,1) + Cqq3_ud_LNP*YddL(0)*YucuL(0,1) + 2*Cqq3_uu_LNP*YucuL(0,0)*YucuL(0,1)).real();
	Cqq3_1112i_LNP = (Cqq3_0u_LNP*YucuL(0,1) + Cqq3_u0_LNP*YucuL(0,1) + Cqq3_du_LNP*YddL(0)*YucuL(0,1) + Cqq3_ud_LNP*YddL(0)*YucuL(0,1) + 2*Cqq3_uu_LNP*YucuL(0,0)*YucuL(0,1)).imag();
	Cqq3_1113r_LNP = (Cqq3_0u_LNP*YucuL(0,2) + Cqq3_u0_LNP*YucuL(0,2) + Cqq3_du_LNP*YddL(0)*YucuL(0,2) + Cqq3_ud_LNP*YddL(0)*YucuL(0,2) + 2*Cqq3_uu_LNP*YucuL(0,0)*YucuL(0,2)).real();
	Cqq3_1113i_LNP = (Cqq3_0u_LNP*YucuL(0,2) + Cqq3_u0_LNP*YucuL(0,2) + Cqq3_du_LNP*YddL(0)*YucuL(0,2) + Cqq3_ud_LNP*YddL(0)*YucuL(0,2) + 2*Cqq3_uu_LNP*YucuL(0,0)*YucuL(0,2)).imag();
	Cqq3_1122r_LNP = (Cqq3_00_LNP + Cqq3_d0_LNP*YddL(0) + Cqq3_0d_LNP*YddL(1) + Cqq3_dd_LNP*YddL(0)*YddL(1) + Cqq3_u0_LNP*YucuL(0,0) + Cqq3_ud_LNP*YddL(1)*YucuL(0,0) + Cqq3_uu_LNP*YucuL(0,1)*YucuL(1,0) + Cqq3_0u_LNP*YucuL(1,1) + Cqq3_du_LNP*YddL(0)*YucuL(1,1) + Cqq3_uu_LNP*YucuL(0,0)*YucuL(1,1)).real();
	Cqq3_1123r_LNP = (Cqq3_uu_LNP*YucuL(0,2)*YucuL(1,0) + Cqq3_0u_LNP*YucuL(1,2) + Cqq3_du_LNP*YddL(0)*YucuL(1,2) + Cqq3_uu_LNP*YucuL(0,0)*YucuL(1,2)).real();
	Cqq3_1123i_LNP = (Cqq3_uu_LNP*YucuL(0,2)*YucuL(1,0) + Cqq3_0u_LNP*YucuL(1,2) + Cqq3_du_LNP*YddL(0)*YucuL(1,2) + Cqq3_uu_LNP*YucuL(0,0)*YucuL(1,2)).imag();
	Cqq3_1133r_LNP = (Cqq3_00_LNP + Cqq3_d0_LNP*YddL(0) + Cqq3_0d_LNP*YddL(2) + Cqq3_dd_LNP*YddL(0)*YddL(2) + Cqq3_u0_LNP*YucuL(0,0) + Cqq3_ud_LNP*YddL(2)*YucuL(0,0) + Cqq3_uu_LNP*YucuL(0,2)*YucuL(2,0) + Cqq3_0u_LNP*YucuL(2,2) + Cqq3_du_LNP*YddL(0)*YucuL(2,2) + Cqq3_uu_LNP*YucuL(0,0)*YucuL(2,2)).real();
	Cqq3_1212r_LNP = (2*Cqq3_uu_LNP*YucuL(0,1)*YucuL(0,1)).real();
	Cqq3_1212i_LNP = (2*Cqq3_uu_LNP*YucuL(0,1)*YucuL(0,1)).imag();
	Cqq3_1213r_LNP = (2*Cqq3_uu_LNP*YucuL(0,1)*YucuL(0,2)).real();
	Cqq3_1213i_LNP = (2*Cqq3_uu_LNP*YucuL(0,1)*YucuL(0,2)).imag();
	Cqq3_1221r_LNP = (Cqq3_00_LNP + Cqq3_d0_LNP*YddL(0) + Cqq3_0d_LNP*YddL(1) + Cqq3_dd_LNP*YddL(0)*YddL(1) + Cqq3_u0_LNP*YucuL(0,0) + Cqq3_ud_LNP*YddL(1)*YucuL(0,0) + Cqq3_uu_LNP*YucuL(0,1)*YucuL(1,0) + Cqq3_0u_LNP*YucuL(1,1) + Cqq3_du_LNP*YddL(0)*YucuL(1,1) + Cqq3_uu_LNP*YucuL(0,0)*YucuL(1,1)).real();
	Cqq3_1222r_LNP = (2*Cqq3_u0_LNP*YucuL(0,1) + 2*Cqq3_ud_LNP*YddL(1)*YucuL(0,1) + 2*Cqq3_uu_LNP*YucuL(0,1)*YucuL(1,1)).real();
	Cqq3_1222i_LNP = (2*Cqq3_u0_LNP*YucuL(0,1) + 2*Cqq3_ud_LNP*YddL(1)*YucuL(0,1) + 2*Cqq3_uu_LNP*YucuL(0,1)*YucuL(1,1)).imag();
	Cqq3_1223r_LNP = (Cqq3_u0_LNP*YucuL(0,2) + Cqq3_ud_LNP*YddL(1)*YucuL(0,2) + Cqq3_uu_LNP*YucuL(0,2)*YucuL(1,1) + Cqq3_uu_LNP*YucuL(0,1)*YucuL(1,2)).real();
	Cqq3_1223i_LNP = (Cqq3_u0_LNP*YucuL(0,2) + Cqq3_ud_LNP*YddL(1)*YucuL(0,2) + Cqq3_uu_LNP*YucuL(0,2)*YucuL(1,1) + Cqq3_uu_LNP*YucuL(0,1)*YucuL(1,2)).imag();
	Cqq3_1231r_LNP = (Cqq3_uu_LNP*YucuL(0,1)*YucuL(2,0) + Cqq3_0u_LNP*YucuL(2,1) + Cqq3_du_LNP*YddL(0)*YucuL(2,1) + Cqq3_uu_LNP*YucuL(0,0)*YucuL(2,1)).real();
	Cqq3_1231i_LNP = (Cqq3_uu_LNP*YucuL(0,1)*YucuL(2,0) + Cqq3_0u_LNP*YucuL(2,1) + Cqq3_du_LNP*YddL(0)*YucuL(2,1) + Cqq3_uu_LNP*YucuL(0,0)*YucuL(2,1)).imag();
	Cqq3_1232r_LNP = (2*Cqq3_uu_LNP*YucuL(0,1)*YucuL(2,1)).real();
	Cqq3_1232i_LNP = (2*Cqq3_uu_LNP*YucuL(0,1)*YucuL(2,1)).imag();
	Cqq3_1233r_LNP = (Cqq3_u0_LNP*YucuL(0,1) + Cqq3_ud_LNP*YddL(2)*YucuL(0,1) + Cqq3_uu_LNP*YucuL(0,2)*YucuL(2,1) + Cqq3_uu_LNP*YucuL(0,1)*YucuL(2,2)).real();
	Cqq3_1233i_LNP = (Cqq3_u0_LNP*YucuL(0,1) + Cqq3_ud_LNP*YddL(2)*YucuL(0,1) + Cqq3_uu_LNP*YucuL(0,2)*YucuL(2,1) + Cqq3_uu_LNP*YucuL(0,1)*YucuL(2,2)).imag();
	Cqq3_1313r_LNP = (2*Cqq3_uu_LNP*YucuL(0,2)*YucuL(0,2)).real();
	Cqq3_1313i_LNP = (2*Cqq3_uu_LNP*YucuL(0,2)*YucuL(0,2)).imag();
	Cqq3_1322r_LNP = (Cqq3_u0_LNP*YucuL(0,2) + Cqq3_ud_LNP*YddL(1)*YucuL(0,2) + Cqq3_uu_LNP*YucuL(0,2)*YucuL(1,1) + Cqq3_uu_LNP*YucuL(0,1)*YucuL(1,2)).real();
	Cqq3_1322i_LNP = (Cqq3_u0_LNP*YucuL(0,2) + Cqq3_ud_LNP*YddL(1)*YucuL(0,2) + Cqq3_uu_LNP*YucuL(0,2)*YucuL(1,1) + Cqq3_uu_LNP*YucuL(0,1)*YucuL(1,2)).imag();
	Cqq3_1323r_LNP = (2*Cqq3_uu_LNP*YucuL(0,2)*YucuL(1,2)).real();
	Cqq3_1323i_LNP = (2*Cqq3_uu_LNP*YucuL(0,2)*YucuL(1,2)).imag();
	Cqq3_1331r_LNP = (Cqq3_00_LNP + Cqq3_d0_LNP*YddL(0) + Cqq3_0d_LNP*YddL(2) + Cqq3_dd_LNP*YddL(0)*YddL(2) + Cqq3_u0_LNP*YucuL(0,0) + Cqq3_ud_LNP*YddL(2)*YucuL(0,0) + Cqq3_uu_LNP*YucuL(0,2)*YucuL(2,0) + Cqq3_0u_LNP*YucuL(2,2) + Cqq3_du_LNP*YddL(0)*YucuL(2,2) + Cqq3_uu_LNP*YucuL(0,0)*YucuL(2,2)).real();
	Cqq3_1332r_LNP = (Cqq3_u0_LNP*YucuL(0,1) + Cqq3_ud_LNP*YddL(2)*YucuL(0,1) + Cqq3_uu_LNP*YucuL(0,2)*YucuL(2,1) + Cqq3_uu_LNP*YucuL(0,1)*YucuL(2,2)).real();
	Cqq3_1332i_LNP = (Cqq3_u0_LNP*YucuL(0,1) + Cqq3_ud_LNP*YddL(2)*YucuL(0,1) + Cqq3_uu_LNP*YucuL(0,2)*YucuL(2,1) + Cqq3_uu_LNP*YucuL(0,1)*YucuL(2,2)).imag();
	Cqq3_1333r_LNP = (2*Cqq3_u0_LNP*YucuL(0,2) + 2*Cqq3_ud_LNP*YddL(2)*YucuL(0,2) + 2*Cqq3_uu_LNP*YucuL(0,2)*YucuL(2,2)).real();
	Cqq3_1333i_LNP = (2*Cqq3_u0_LNP*YucuL(0,2) + 2*Cqq3_ud_LNP*YddL(2)*YucuL(0,2) + 2*Cqq3_uu_LNP*YucuL(0,2)*YucuL(2,2)).imag();
	Cqq3_2222r_LNP = (2*Cqq3_00_LNP + 2*Cqq3_0d_LNP*YddL(1) + 2*Cqq3_d0_LNP*YddL(1) + 2*Cqq3_dd_LNP*YddL(1)*YddL(1) + 2*Cqq3_0u_LNP*YucuL(1,1) + 2*Cqq3_u0_LNP*YucuL(1,1) + 2*Cqq3_du_LNP*YddL(1)*YucuL(1,1) + 2*Cqq3_ud_LNP*YddL(1)*YucuL(1,1) + 2*Cqq3_uu_LNP*YucuL(1,1)*YucuL(1,1)).real();
	Cqq3_2223r_LNP = (Cqq3_0u_LNP*YucuL(1,2) + Cqq3_u0_LNP*YucuL(1,2) + Cqq3_du_LNP*YddL(1)*YucuL(1,2) + Cqq3_ud_LNP*YddL(1)*YucuL(1,2) + 2*Cqq3_uu_LNP*YucuL(1,1)*YucuL(1,2)).real();
	Cqq3_2223i_LNP = (Cqq3_0u_LNP*YucuL(1,2) + Cqq3_u0_LNP*YucuL(1,2) + Cqq3_du_LNP*YddL(1)*YucuL(1,2) + Cqq3_ud_LNP*YddL(1)*YucuL(1,2) + 2*Cqq3_uu_LNP*YucuL(1,1)*YucuL(1,2)).imag();
	Cqq3_2233r_LNP = (Cqq3_00_LNP + Cqq3_d0_LNP*YddL(1) + Cqq3_0d_LNP*YddL(2) + Cqq3_dd_LNP*YddL(1)*YddL(2) + Cqq3_u0_LNP*YucuL(1,1) + Cqq3_ud_LNP*YddL(2)*YucuL(1,1) + Cqq3_uu_LNP*YucuL(1,2)*YucuL(2,1) + Cqq3_0u_LNP*YucuL(2,2) + Cqq3_du_LNP*YddL(1)*YucuL(2,2) + Cqq3_uu_LNP*YucuL(1,1)*YucuL(2,2)).real();
	Cqq3_2323r_LNP = (2*Cqq3_uu_LNP*YucuL(1,2)*YucuL(1,2)).real();
	Cqq3_2323i_LNP = (2*Cqq3_uu_LNP*YucuL(1,2)*YucuL(1,2)).imag();
	Cqq3_2332r_LNP = (Cqq3_00_LNP + Cqq3_d0_LNP*YddL(1) + Cqq3_0d_LNP*YddL(2) + Cqq3_dd_LNP*YddL(1)*YddL(2) + Cqq3_u0_LNP*YucuL(1,1) + Cqq3_ud_LNP*YddL(2)*YucuL(1,1) + Cqq3_uu_LNP*YucuL(1,2)*YucuL(2,1) + Cqq3_0u_LNP*YucuL(2,2) + Cqq3_du_LNP*YddL(1)*YucuL(2,2) + Cqq3_uu_LNP*YucuL(1,1)*YucuL(2,2)).real();
	Cqq3_2333r_LNP = (2*Cqq3_u0_LNP*YucuL(1,2) + 2*Cqq3_ud_LNP*YddL(2)*YucuL(1,2) + 2*Cqq3_uu_LNP*YucuL(1,2)*YucuL(2,2)).real();
	Cqq3_2333i_LNP = (2*Cqq3_u0_LNP*YucuL(1,2) + 2*Cqq3_ud_LNP*YddL(2)*YucuL(1,2) + 2*Cqq3_uu_LNP*YucuL(1,2)*YucuL(2,2)).imag();
	Cqq3_3333r_LNP = (2*Cqq3_00_LNP + 2*Cqq3_0d_LNP*YddL(2) + 2*Cqq3_d0_LNP*YddL(2) + 2*Cqq3_dd_LNP*YddL(2)*YddL(2) + 2*Cqq3_0u_LNP*YucuL(2,2) + 2*Cqq3_u0_LNP*YucuL(2,2) + 2*Cqq3_du_LNP*YddL(2)*YucuL(2,2) + 2*Cqq3_ud_LNP*YddL(2)*YucuL(2,2) + 2*Cqq3_uu_LNP*YucuL(2,2)*YucuL(2,2)).real();

	Cuu_1111r_LNP = (2*Cuu_00_LNP + 2*Cuu_0u_LNP*YuucL(0,0) + 2*Cuu_u0_LNP*YuucL(0,0) + 2*Cuu_uu_LNP*YuucL(0,0)*YuucL(0,0)).real();
	Cuu_1112r_LNP = (Cuu_0u_LNP*YuucL(0,1) + Cuu_u0_LNP*YuucL(0,1) + 2*Cuu_uu_LNP*YuucL(0,0)*YuucL(0,1)).real();
	Cuu_1112i_LNP = (Cuu_0u_LNP*YuucL(0,1) + Cuu_u0_LNP*YuucL(0,1) + 2*Cuu_uu_LNP*YuucL(0,0)*YuucL(0,1)).imag();
	Cuu_1113r_LNP = (Cuu_0u_LNP*YuucL(0,2) + Cuu_u0_LNP*YuucL(0,2) + 2*Cuu_uu_LNP*YuucL(0,0)*YuucL(0,2)).real();
	Cuu_1113i_LNP = (Cuu_0u_LNP*YuucL(0,2) + Cuu_u0_LNP*YuucL(0,2) + 2*Cuu_uu_LNP*YuucL(0,0)*YuucL(0,2)).imag();
	Cuu_1122r_LNP = (Cuu_00_LNP + Cuu_u0_LNP*YuucL(0,0) + Cuu_uu_LNP*YuucL(0,1)*YuucL(1,0) + Cuu_0u_LNP*YuucL(1,1) + Cuu_uu_LNP*YuucL(0,0)*YuucL(1,1)).real();
	Cuu_1123r_LNP = (Cuu_uu_LNP*YuucL(0,2)*YuucL(1,0) + Cuu_0u_LNP*YuucL(1,2) + Cuu_uu_LNP*YuucL(0,0)*YuucL(1,2)).real();
	Cuu_1123i_LNP = (Cuu_uu_LNP*YuucL(0,2)*YuucL(1,0) + Cuu_0u_LNP*YuucL(1,2) + Cuu_uu_LNP*YuucL(0,0)*YuucL(1,2)).imag();
	Cuu_1133r_LNP = (Cuu_00_LNP + Cuu_u0_LNP*YuucL(0,0) + Cuu_uu_LNP*YuucL(0,2)*YuucL(2,0) + Cuu_0u_LNP*YuucL(2,2) + Cuu_uu_LNP*YuucL(0,0)*YuucL(2,2)).real();
	Cuu_1212r_LNP = (2*Cuu_uu_LNP*YuucL(0,1)*YuucL(0,1)).real();
	Cuu_1212i_LNP = (2*Cuu_uu_LNP*YuucL(0,1)*YuucL(0,1)).imag();
	Cuu_1213r_LNP = (2*Cuu_uu_LNP*YuucL(0,1)*YuucL(0,2)).real();
	Cuu_1213i_LNP = (2*Cuu_uu_LNP*YuucL(0,1)*YuucL(0,2)).imag();
	Cuu_1221r_LNP = (Cuu_00_LNP + Cuu_u0_LNP*YuucL(0,0) + Cuu_uu_LNP*YuucL(0,1)*YuucL(1,0) + Cuu_0u_LNP*YuucL(1,1) + Cuu_uu_LNP*YuucL(0,0)*YuucL(1,1)).real();
	Cuu_1222r_LNP = (2*Cuu_u0_LNP*YuucL(0,1) + 2*Cuu_uu_LNP*YuucL(0,1)*YuucL(1,1)).real();
	Cuu_1222i_LNP = (2*Cuu_u0_LNP*YuucL(0,1) + 2*Cuu_uu_LNP*YuucL(0,1)*YuucL(1,1)).imag();
	Cuu_1223r_LNP = (Cuu_u0_LNP*YuucL(0,2) + Cuu_uu_LNP*YuucL(0,2)*YuucL(1,1) + Cuu_uu_LNP*YuucL(0,1)*YuucL(1,2)).real();
	Cuu_1223i_LNP = (Cuu_u0_LNP*YuucL(0,2) + Cuu_uu_LNP*YuucL(0,2)*YuucL(1,1) + Cuu_uu_LNP*YuucL(0,1)*YuucL(1,2)).imag();
	Cuu_1231r_LNP = (Cuu_uu_LNP*YuucL(0,1)*YuucL(2,0) + Cuu_0u_LNP*YuucL(2,1) + Cuu_uu_LNP*YuucL(0,0)*YuucL(2,1)).real();
	Cuu_1231i_LNP = (Cuu_uu_LNP*YuucL(0,1)*YuucL(2,0) + Cuu_0u_LNP*YuucL(2,1) + Cuu_uu_LNP*YuucL(0,0)*YuucL(2,1)).imag();
	Cuu_1232r_LNP = (2*Cuu_uu_LNP*YuucL(0,1)*YuucL(2,1)).real();
	Cuu_1232i_LNP = (2*Cuu_uu_LNP*YuucL(0,1)*YuucL(2,1)).imag();
	Cuu_1233r_LNP = (Cuu_u0_LNP*YuucL(0,1) + Cuu_uu_LNP*YuucL(0,2)*YuucL(2,1) + Cuu_uu_LNP*YuucL(0,1)*YuucL(2,2)).real();
	Cuu_1233i_LNP = (Cuu_u0_LNP*YuucL(0,1) + Cuu_uu_LNP*YuucL(0,2)*YuucL(2,1) + Cuu_uu_LNP*YuucL(0,1)*YuucL(2,2)).imag();
	Cuu_1313r_LNP = (2*Cuu_uu_LNP*YuucL(0,2)*YuucL(0,2)).real();
	Cuu_1313i_LNP = (2*Cuu_uu_LNP*YuucL(0,2)*YuucL(0,2)).imag();
	Cuu_1322r_LNP = (Cuu_u0_LNP*YuucL(0,2) + Cuu_uu_LNP*YuucL(0,2)*YuucL(1,1) + Cuu_uu_LNP*YuucL(0,1)*YuucL(1,2)).real();
	Cuu_1322i_LNP = (Cuu_u0_LNP*YuucL(0,2) + Cuu_uu_LNP*YuucL(0,2)*YuucL(1,1) + Cuu_uu_LNP*YuucL(0,1)*YuucL(1,2)).imag();
	Cuu_1323r_LNP = (2*Cuu_uu_LNP*YuucL(0,2)*YuucL(1,2)).real();
	Cuu_1323i_LNP = (2*Cuu_uu_LNP*YuucL(0,2)*YuucL(1,2)).imag();
	Cuu_1331r_LNP = (Cuu_00_LNP + Cuu_u0_LNP*YuucL(0,0) + Cuu_uu_LNP*YuucL(0,2)*YuucL(2,0) + Cuu_0u_LNP*YuucL(2,2) + Cuu_uu_LNP*YuucL(0,0)*YuucL(2,2)).real();
	Cuu_1332r_LNP = (Cuu_u0_LNP*YuucL(0,1) + Cuu_uu_LNP*YuucL(0,2)*YuucL(2,1) + Cuu_uu_LNP*YuucL(0,1)*YuucL(2,2)).real();
	Cuu_1332i_LNP = (Cuu_u0_LNP*YuucL(0,1) + Cuu_uu_LNP*YuucL(0,2)*YuucL(2,1) + Cuu_uu_LNP*YuucL(0,1)*YuucL(2,2)).imag();
	Cuu_1333r_LNP = (2*Cuu_u0_LNP*YuucL(0,2) + 2*Cuu_uu_LNP*YuucL(0,2)*YuucL(2,2)).real();
	Cuu_1333i_LNP = (2*Cuu_u0_LNP*YuucL(0,2) + 2*Cuu_uu_LNP*YuucL(0,2)*YuucL(2,2)).imag();
	Cuu_2222r_LNP = (2*Cuu_00_LNP + 2*Cuu_0u_LNP*YuucL(1,1) + 2*Cuu_u0_LNP*YuucL(1,1) + 2*Cuu_uu_LNP*YuucL(1,1)*YuucL(1,1)).real();
	Cuu_2223r_LNP = (Cuu_0u_LNP*YuucL(1,2) + Cuu_u0_LNP*YuucL(1,2) + 2*Cuu_uu_LNP*YuucL(1,1)*YuucL(1,2)).real();
	Cuu_2223i_LNP = (Cuu_0u_LNP*YuucL(1,2) + Cuu_u0_LNP*YuucL(1,2) + 2*Cuu_uu_LNP*YuucL(1,1)*YuucL(1,2)).imag();
	Cuu_2233r_LNP = (Cuu_00_LNP + Cuu_u0_LNP*YuucL(1,1) + Cuu_uu_LNP*YuucL(1,2)*YuucL(2,1) + Cuu_0u_LNP*YuucL(2,2) + Cuu_uu_LNP*YuucL(1,1)*YuucL(2,2)).real();
	Cuu_2323r_LNP = (2*Cuu_uu_LNP*YuucL(1,2)*YuucL(1,2)).real();
	Cuu_2323i_LNP = (2*Cuu_uu_LNP*YuucL(1,2)*YuucL(1,2)).imag();
	Cuu_2332r_LNP = (Cuu_00_LNP + Cuu_u0_LNP*YuucL(1,1) + Cuu_uu_LNP*YuucL(1,2)*YuucL(2,1) + Cuu_0u_LNP*YuucL(2,2) + Cuu_uu_LNP*YuucL(1,1)*YuucL(2,2)).real();
	Cuu_2333r_LNP = (2*Cuu_u0_LNP*YuucL(1,2) + 2*Cuu_uu_LNP*YuucL(1,2)*YuucL(2,2)).real();
	Cuu_2333i_LNP = (2*Cuu_u0_LNP*YuucL(1,2) + 2*Cuu_uu_LNP*YuucL(1,2)*YuucL(2,2)).imag();
	Cuu_3333r_LNP = (2*Cuu_00_LNP + 2*Cuu_0u_LNP*YuucL(2,2) + 2*Cuu_u0_LNP*YuucL(2,2) + 2*Cuu_uu_LNP*YuucL(2,2)*YuucL(2,2)).real();

	Cdd_1111r_LNP = (2*Cdd_00_LNP + 2*Cdd_0d_LNP*YddL(0) + 2*Cdd_d0_LNP*YddL(0) + 2*Cdd_dd_LNP*YddL(0)*YddL(0)).real();
	Cdd_1122r_LNP = (Cdd_00_LNP + Cdd_d0_LNP*YddL(0) + Cdd_0d_LNP*YddL(1) + Cdd_dd_LNP*YddL(0)*YddL(1)).real();
	Cdd_1133r_LNP = (Cdd_00_LNP + Cdd_d0_LNP*YddL(0) + Cdd_0d_LNP*YddL(2) + Cdd_dd_LNP*YddL(0)*YddL(2)).real();
	Cdd_1221r_LNP = (Cdd_00_LNP + Cdd_d0_LNP*YddL(0) + Cdd_0d_LNP*YddL(1) + Cdd_dd_LNP*YddL(0)*YddL(1)).real();
	Cdd_1331r_LNP = (Cdd_00_LNP + Cdd_d0_LNP*YddL(0) + Cdd_0d_LNP*YddL(2) + Cdd_dd_LNP*YddL(0)*YddL(2)).real();
	Cdd_2222r_LNP = (2*Cdd_00_LNP + 2*Cdd_0d_LNP*YddL(1) + 2*Cdd_d0_LNP*YddL(1) + 2*Cdd_dd_LNP*YddL(1)*YddL(1)).real();
	Cdd_2233r_LNP = (Cdd_00_LNP + Cdd_d0_LNP*YddL(1) + Cdd_0d_LNP*YddL(2) + Cdd_dd_LNP*YddL(1)*YddL(2)).real();
	Cdd_2332r_LNP = (Cdd_00_LNP + Cdd_d0_LNP*YddL(1) + Cdd_0d_LNP*YddL(2) + Cdd_dd_LNP*YddL(1)*YddL(2)).real();
	Cdd_3333r_LNP = (2*Cdd_00_LNP + 2*Cdd_0d_LNP*YddL(2) + 2*Cdd_d0_LNP*YddL(2) + 2*Cdd_dd_LNP*YddL(2)*YddL(2)).real();

	Cud1_1111r_LNP = (Cud1_00_LNP + Cud1_0d_LNP*YddL(0) + Cud1_u0_LNP*YuucL(0,0) + Cud1_ud_LNP*YddL(0)*YuucL(0,0)).real();
	Cud1_1122r_LNP = (Cud1_00_LNP + Cud1_0d_LNP*YddL(1) + Cud1_u0_LNP*YuucL(0,0) + Cud1_ud_LNP*YddL(1)*YuucL(0,0)).real();
	Cud1_1133r_LNP = (Cud1_00_LNP + Cud1_0d_LNP*YddL(2) + Cud1_u0_LNP*YuucL(0,0) + Cud1_ud_LNP*YddL(2)*YuucL(0,0)).real();
	Cud1_1211r_LNP = (Cud1_u0_LNP*YuucL(0,1) + Cud1_ud_LNP*YddL(0)*YuucL(0,1)).real();
	Cud1_1211i_LNP = (Cud1_u0_LNP*YuucL(0,1) + Cud1_ud_LNP*YddL(0)*YuucL(0,1)).imag();
	Cud1_1222r_LNP = (Cud1_u0_LNP*YuucL(0,1) + Cud1_ud_LNP*YddL(1)*YuucL(0,1)).real();
	Cud1_1222i_LNP = (Cud1_u0_LNP*YuucL(0,1) + Cud1_ud_LNP*YddL(1)*YuucL(0,1)).imag();
	Cud1_1233r_LNP = (Cud1_u0_LNP*YuucL(0,1) + Cud1_ud_LNP*YddL(2)*YuucL(0,1)).real();
	Cud1_1233i_LNP = (Cud1_u0_LNP*YuucL(0,1) + Cud1_ud_LNP*YddL(2)*YuucL(0,1)).imag();
	Cud1_1311r_LNP = (Cud1_u0_LNP*YuucL(0,2) + Cud1_ud_LNP*YddL(0)*YuucL(0,2)).real();
	Cud1_1311i_LNP = (Cud1_u0_LNP*YuucL(0,2) + Cud1_ud_LNP*YddL(0)*YuucL(0,2)).imag();
	Cud1_1322r_LNP = (Cud1_u0_LNP*YuucL(0,2) + Cud1_ud_LNP*YddL(1)*YuucL(0,2)).real();
	Cud1_1322i_LNP = (Cud1_u0_LNP*YuucL(0,2) + Cud1_ud_LNP*YddL(1)*YuucL(0,2)).imag();
	Cud1_1333r_LNP = (Cud1_u0_LNP*YuucL(0,2) + Cud1_ud_LNP*YddL(2)*YuucL(0,2)).real();
	Cud1_1333i_LNP = (Cud1_u0_LNP*YuucL(0,2) + Cud1_ud_LNP*YddL(2)*YuucL(0,2)).imag();
	Cud1_2211r_LNP = (Cud1_00_LNP + Cud1_0d_LNP*YddL(0) + Cud1_u0_LNP*YuucL(1,1) + Cud1_ud_LNP*YddL(0)*YuucL(1,1)).real();
	Cud1_2222r_LNP = (Cud1_00_LNP + Cud1_0d_LNP*YddL(1) + Cud1_u0_LNP*YuucL(1,1) + Cud1_ud_LNP*YddL(1)*YuucL(1,1)).real();
	Cud1_2233r_LNP = (Cud1_00_LNP + Cud1_0d_LNP*YddL(2) + Cud1_u0_LNP*YuucL(1,1) + Cud1_ud_LNP*YddL(2)*YuucL(1,1)).real();
	Cud1_2311r_LNP = (Cud1_u0_LNP*YuucL(1,2) + Cud1_ud_LNP*YddL(0)*YuucL(1,2)).real();
	Cud1_2311i_LNP = (Cud1_u0_LNP*YuucL(1,2) + Cud1_ud_LNP*YddL(0)*YuucL(1,2)).imag();
	Cud1_2322r_LNP = (Cud1_u0_LNP*YuucL(1,2) + Cud1_ud_LNP*YddL(1)*YuucL(1,2)).real();
	Cud1_2322i_LNP = (Cud1_u0_LNP*YuucL(1,2) + Cud1_ud_LNP*YddL(1)*YuucL(1,2)).imag();
	Cud1_2333r_LNP = (Cud1_u0_LNP*YuucL(1,2) + Cud1_ud_LNP*YddL(2)*YuucL(1,2)).real();
	Cud1_2333i_LNP = (Cud1_u0_LNP*YuucL(1,2) + Cud1_ud_LNP*YddL(2)*YuucL(1,2)).imag();
	Cud1_3311r_LNP = (Cud1_00_LNP + Cud1_0d_LNP*YddL(0) + Cud1_u0_LNP*YuucL(2,2) + Cud1_ud_LNP*YddL(0)*YuucL(2,2)).real();
	Cud1_3322r_LNP = (Cud1_00_LNP + Cud1_0d_LNP*YddL(1) + Cud1_u0_LNP*YuucL(2,2) + Cud1_ud_LNP*YddL(1)*YuucL(2,2)).real();
	Cud1_3333r_LNP = (Cud1_00_LNP + Cud1_0d_LNP*YddL(2) + Cud1_u0_LNP*YuucL(2,2) + Cud1_ud_LNP*YddL(2)*YuucL(2,2)).real();

	Cud8_1111r_LNP = (Cud8_00_LNP + Cud8_0d_LNP*YddL(0) + Cud8_u0_LNP*YuucL(0,0) + Cud8_ud_LNP*YddL(0)*YuucL(0,0)).real();
	Cud8_1122r_LNP = (Cud8_00_LNP + Cud8_0d_LNP*YddL(1) + Cud8_u0_LNP*YuucL(0,0) + Cud8_ud_LNP*YddL(1)*YuucL(0,0)).real();
	Cud8_1133r_LNP = (Cud8_00_LNP + Cud8_0d_LNP*YddL(2) + Cud8_u0_LNP*YuucL(0,0) + Cud8_ud_LNP*YddL(2)*YuucL(0,0)).real();
	Cud8_1211r_LNP = (Cud8_u0_LNP*YuucL(0,1) + Cud8_ud_LNP*YddL(0)*YuucL(0,1)).real();
	Cud8_1211i_LNP = (Cud8_u0_LNP*YuucL(0,1) + Cud8_ud_LNP*YddL(0)*YuucL(0,1)).imag();
	Cud8_1222r_LNP = (Cud8_u0_LNP*YuucL(0,1) + Cud8_ud_LNP*YddL(1)*YuucL(0,1)).real();
	Cud8_1222i_LNP = (Cud8_u0_LNP*YuucL(0,1) + Cud8_ud_LNP*YddL(1)*YuucL(0,1)).imag();
	Cud8_1233r_LNP = (Cud8_u0_LNP*YuucL(0,1) + Cud8_ud_LNP*YddL(2)*YuucL(0,1)).real();
	Cud8_1233i_LNP = (Cud8_u0_LNP*YuucL(0,1) + Cud8_ud_LNP*YddL(2)*YuucL(0,1)).imag();
	Cud8_1311r_LNP = (Cud8_u0_LNP*YuucL(0,2) + Cud8_ud_LNP*YddL(0)*YuucL(0,2)).real();
	Cud8_1311i_LNP = (Cud8_u0_LNP*YuucL(0,2) + Cud8_ud_LNP*YddL(0)*YuucL(0,2)).imag();
	Cud8_1322r_LNP = (Cud8_u0_LNP*YuucL(0,2) + Cud8_ud_LNP*YddL(1)*YuucL(0,2)).real();
	Cud8_1322i_LNP = (Cud8_u0_LNP*YuucL(0,2) + Cud8_ud_LNP*YddL(1)*YuucL(0,2)).imag();
	Cud8_1333r_LNP = (Cud8_u0_LNP*YuucL(0,2) + Cud8_ud_LNP*YddL(2)*YuucL(0,2)).real();
	Cud8_1333i_LNP = (Cud8_u0_LNP*YuucL(0,2) + Cud8_ud_LNP*YddL(2)*YuucL(0,2)).imag();
	Cud8_2211r_LNP = (Cud8_00_LNP + Cud8_0d_LNP*YddL(0) + Cud8_u0_LNP*YuucL(1,1) + Cud8_ud_LNP*YddL(0)*YuucL(1,1)).real();
	Cud8_2222r_LNP = (Cud8_00_LNP + Cud8_0d_LNP*YddL(1) + Cud8_u0_LNP*YuucL(1,1) + Cud8_ud_LNP*YddL(1)*YuucL(1,1)).real();
	Cud8_2233r_LNP = (Cud8_00_LNP + Cud8_0d_LNP*YddL(2) + Cud8_u0_LNP*YuucL(1,1) + Cud8_ud_LNP*YddL(2)*YuucL(1,1)).real();
	Cud8_2311r_LNP = (Cud8_u0_LNP*YuucL(1,2) + Cud8_ud_LNP*YddL(0)*YuucL(1,2)).real();
	Cud8_2311i_LNP = (Cud8_u0_LNP*YuucL(1,2) + Cud8_ud_LNP*YddL(0)*YuucL(1,2)).imag();
	Cud8_2322r_LNP = (Cud8_u0_LNP*YuucL(1,2) + Cud8_ud_LNP*YddL(1)*YuucL(1,2)).real();
	Cud8_2322i_LNP = (Cud8_u0_LNP*YuucL(1,2) + Cud8_ud_LNP*YddL(1)*YuucL(1,2)).imag();
	Cud8_2333r_LNP = (Cud8_u0_LNP*YuucL(1,2) + Cud8_ud_LNP*YddL(2)*YuucL(1,2)).real();
	Cud8_2333i_LNP = (Cud8_u0_LNP*YuucL(1,2) + Cud8_ud_LNP*YddL(2)*YuucL(1,2)).imag();
	Cud8_3311r_LNP = (Cud8_00_LNP + Cud8_0d_LNP*YddL(0) + Cud8_u0_LNP*YuucL(2,2) + Cud8_ud_LNP*YddL(0)*YuucL(2,2)).real();
	Cud8_3322r_LNP = (Cud8_00_LNP + Cud8_0d_LNP*YddL(1) + Cud8_u0_LNP*YuucL(2,2) + Cud8_ud_LNP*YddL(1)*YuucL(2,2)).real();
	Cud8_3333r_LNP = (Cud8_00_LNP + Cud8_0d_LNP*YddL(2) + Cud8_u0_LNP*YuucL(2,2) + Cud8_ud_LNP*YddL(2)*YuucL(2,2)).real();

	Cqu1_1111r_LNP = (Cqu1_00_LNP + Cqu1_d0_LNP*YddL(0) + Cqu1_u0_LNP*YucuL(0,0) + Cqu1_0u_LNP*YuucL(0,0) + Cqu1_du_LNP*YddL(0)*YuucL(0,0) + Cqu1_uu_LNP*YucuL(0,0)*YuucL(0,0)).real();
	Cqu1_1112r_LNP = (Cqu1_0u_LNP*YuucL(0,1) + Cqu1_du_LNP*YddL(0)*YuucL(0,1) + Cqu1_uu_LNP*YucuL(0,0)*YuucL(0,1)).real();
	Cqu1_1112i_LNP = (Cqu1_0u_LNP*YuucL(0,1) + Cqu1_du_LNP*YddL(0)*YuucL(0,1) + Cqu1_uu_LNP*YucuL(0,0)*YuucL(0,1)).imag();
	Cqu1_1113r_LNP = (Cqu1_0u_LNP*YuucL(0,2) + Cqu1_du_LNP*YddL(0)*YuucL(0,2) + Cqu1_uu_LNP*YucuL(0,0)*YuucL(0,2)).real();
	Cqu1_1113i_LNP = (Cqu1_0u_LNP*YuucL(0,2) + Cqu1_du_LNP*YddL(0)*YuucL(0,2) + Cqu1_uu_LNP*YucuL(0,0)*YuucL(0,2)).imag();
	Cqu1_1122r_LNP = (Cqu1_00_LNP + Cqu1_d0_LNP*YddL(0) + Cqu1_u0_LNP*YucuL(0,0) + Cqu1_0u_LNP*YuucL(1,1) + Cqu1_du_LNP*YddL(0)*YuucL(1,1) + Cqu1_uu_LNP*YucuL(0,0)*YuucL(1,1)).real();
	Cqu1_1123r_LNP = (Cqu1_0u_LNP*YuucL(1,2) + Cqu1_du_LNP*YddL(0)*YuucL(1,2) + Cqu1_uu_LNP*YucuL(0,0)*YuucL(1,2)).real();
	Cqu1_1123i_LNP = (Cqu1_0u_LNP*YuucL(1,2) + Cqu1_du_LNP*YddL(0)*YuucL(1,2) + Cqu1_uu_LNP*YucuL(0,0)*YuucL(1,2)).imag();
	Cqu1_1133r_LNP = (Cqu1_00_LNP + Cqu1_d0_LNP*YddL(0) + Cqu1_u0_LNP*YucuL(0,0) + Cqu1_0u_LNP*YuucL(2,2) + Cqu1_du_LNP*YddL(0)*YuucL(2,2) + Cqu1_uu_LNP*YucuL(0,0)*YuucL(2,2)).real();
	Cqu1_1211r_LNP = (Cqu1_u0_LNP*YucuL(0,1) + Cqu1_uu_LNP*YucuL(0,1)*YuucL(0,0)).real();
	Cqu1_1211i_LNP = (Cqu1_u0_LNP*YucuL(0,1) + Cqu1_uu_LNP*YucuL(0,1)*YuucL(0,0)).imag();
	Cqu1_1212r_LNP = (Cqu1_uu_LNP*YucuL(0,1)*YuucL(0,1)).real();
	Cqu1_1212i_LNP = (Cqu1_uu_LNP*YucuL(0,1)*YuucL(0,1)).imag();
	Cqu1_1213r_LNP = (Cqu1_uu_LNP*YucuL(0,1)*YuucL(0,2)).real();
	Cqu1_1213i_LNP = (Cqu1_uu_LNP*YucuL(0,1)*YuucL(0,2)).imag();
	Cqu1_1221r_LNP = (Cqu1_uu_LNP*YucuL(0,1)*YuucL(1,0)).real();
	Cqu1_1221i_LNP = (Cqu1_uu_LNP*YucuL(0,1)*YuucL(1,0)).imag();
	Cqu1_1222r_LNP = (Cqu1_u0_LNP*YucuL(0,1) + Cqu1_uu_LNP*YucuL(0,1)*YuucL(1,1)).real();
	Cqu1_1222i_LNP = (Cqu1_u0_LNP*YucuL(0,1) + Cqu1_uu_LNP*YucuL(0,1)*YuucL(1,1)).imag();
	Cqu1_1223r_LNP = (Cqu1_uu_LNP*YucuL(0,1)*YuucL(1,2)).real();
	Cqu1_1223i_LNP = (Cqu1_uu_LNP*YucuL(0,1)*YuucL(1,2)).imag();
	Cqu1_1231r_LNP = (Cqu1_uu_LNP*YucuL(0,1)*YuucL(2,0)).real();
	Cqu1_1231i_LNP = (Cqu1_uu_LNP*YucuL(0,1)*YuucL(2,0)).imag();
	Cqu1_1232r_LNP = (Cqu1_uu_LNP*YucuL(0,1)*YuucL(2,1)).real();
	Cqu1_1232i_LNP = (Cqu1_uu_LNP*YucuL(0,1)*YuucL(2,1)).imag();
	Cqu1_1233r_LNP = (Cqu1_u0_LNP*YucuL(0,1) + Cqu1_uu_LNP*YucuL(0,1)*YuucL(2,2)).real();
	Cqu1_1233i_LNP = (Cqu1_u0_LNP*YucuL(0,1) + Cqu1_uu_LNP*YucuL(0,1)*YuucL(2,2)).imag();
	Cqu1_1311r_LNP = (Cqu1_u0_LNP*YucuL(0,2) + Cqu1_uu_LNP*YucuL(0,2)*YuucL(0,0)).real();
	Cqu1_1311i_LNP = (Cqu1_u0_LNP*YucuL(0,2) + Cqu1_uu_LNP*YucuL(0,2)*YuucL(0,0)).imag();
	Cqu1_1312r_LNP = (Cqu1_uu_LNP*YucuL(0,2)*YuucL(0,1)).real();
	Cqu1_1312i_LNP = (Cqu1_uu_LNP*YucuL(0,2)*YuucL(0,1)).imag();
	Cqu1_1313r_LNP = (Cqu1_uu_LNP*YucuL(0,2)*YuucL(0,2)).real();
	Cqu1_1313i_LNP = (Cqu1_uu_LNP*YucuL(0,2)*YuucL(0,2)).imag();
	Cqu1_1321r_LNP = (Cqu1_uu_LNP*YucuL(0,2)*YuucL(1,0)).real();
	Cqu1_1321i_LNP = (Cqu1_uu_LNP*YucuL(0,2)*YuucL(1,0)).imag();
	Cqu1_1322r_LNP = (Cqu1_u0_LNP*YucuL(0,2) + Cqu1_uu_LNP*YucuL(0,2)*YuucL(1,1)).real();
	Cqu1_1322i_LNP = (Cqu1_u0_LNP*YucuL(0,2) + Cqu1_uu_LNP*YucuL(0,2)*YuucL(1,1)).imag();
	Cqu1_1323r_LNP = (Cqu1_uu_LNP*YucuL(0,2)*YuucL(1,2)).real();
	Cqu1_1323i_LNP = (Cqu1_uu_LNP*YucuL(0,2)*YuucL(1,2)).imag();
	Cqu1_1331r_LNP = (Cqu1_uu_LNP*YucuL(0,2)*YuucL(2,0)).real();
	Cqu1_1331i_LNP = (Cqu1_uu_LNP*YucuL(0,2)*YuucL(2,0)).imag();
	Cqu1_1332r_LNP = (Cqu1_uu_LNP*YucuL(0,2)*YuucL(2,1)).real();
	Cqu1_1332i_LNP = (Cqu1_uu_LNP*YucuL(0,2)*YuucL(2,1)).imag();
	Cqu1_1333r_LNP = (Cqu1_u0_LNP*YucuL(0,2) + Cqu1_uu_LNP*YucuL(0,2)*YuucL(2,2)).real();
	Cqu1_1333i_LNP = (Cqu1_u0_LNP*YucuL(0,2) + Cqu1_uu_LNP*YucuL(0,2)*YuucL(2,2)).imag();
	Cqu1_2211r_LNP = (Cqu1_00_LNP + Cqu1_d0_LNP*YddL(1) + Cqu1_u0_LNP*YucuL(1,1) + Cqu1_0u_LNP*YuucL(0,0) + Cqu1_du_LNP*YddL(1)*YuucL(0,0) + Cqu1_uu_LNP*YucuL(1,1)*YuucL(0,0)).real();
	Cqu1_2212r_LNP = (Cqu1_0u_LNP*YuucL(0,1) + Cqu1_du_LNP*YddL(1)*YuucL(0,1) + Cqu1_uu_LNP*YucuL(1,1)*YuucL(0,1)).real();
	Cqu1_2212i_LNP = (Cqu1_0u_LNP*YuucL(0,1) + Cqu1_du_LNP*YddL(1)*YuucL(0,1) + Cqu1_uu_LNP*YucuL(1,1)*YuucL(0,1)).imag();
	Cqu1_2213r_LNP = (Cqu1_0u_LNP*YuucL(0,2) + Cqu1_du_LNP*YddL(1)*YuucL(0,2) + Cqu1_uu_LNP*YucuL(1,1)*YuucL(0,2)).real();
	Cqu1_2213i_LNP = (Cqu1_0u_LNP*YuucL(0,2) + Cqu1_du_LNP*YddL(1)*YuucL(0,2) + Cqu1_uu_LNP*YucuL(1,1)*YuucL(0,2)).imag();
	Cqu1_2222r_LNP = (Cqu1_00_LNP + Cqu1_d0_LNP*YddL(1) + Cqu1_u0_LNP*YucuL(1,1) + Cqu1_0u_LNP*YuucL(1,1) + Cqu1_du_LNP*YddL(1)*YuucL(1,1) + Cqu1_uu_LNP*YucuL(1,1)*YuucL(1,1)).real();
	Cqu1_2223r_LNP = (Cqu1_0u_LNP*YuucL(1,2) + Cqu1_du_LNP*YddL(1)*YuucL(1,2) + Cqu1_uu_LNP*YucuL(1,1)*YuucL(1,2)).real();
	Cqu1_2223i_LNP = (Cqu1_0u_LNP*YuucL(1,2) + Cqu1_du_LNP*YddL(1)*YuucL(1,2) + Cqu1_uu_LNP*YucuL(1,1)*YuucL(1,2)).imag();
	Cqu1_2233r_LNP = (Cqu1_00_LNP + Cqu1_d0_LNP*YddL(1) + Cqu1_u0_LNP*YucuL(1,1) + Cqu1_0u_LNP*YuucL(2,2) + Cqu1_du_LNP*YddL(1)*YuucL(2,2) + Cqu1_uu_LNP*YucuL(1,1)*YuucL(2,2)).real();
	Cqu1_2311r_LNP = (Cqu1_u0_LNP*YucuL(1,2) + Cqu1_uu_LNP*YucuL(1,2)*YuucL(0,0)).real();
	Cqu1_2311i_LNP = (Cqu1_u0_LNP*YucuL(1,2) + Cqu1_uu_LNP*YucuL(1,2)*YuucL(0,0)).imag();
	Cqu1_2312r_LNP = (Cqu1_uu_LNP*YucuL(1,2)*YuucL(0,1)).real();
	Cqu1_2312i_LNP = (Cqu1_uu_LNP*YucuL(1,2)*YuucL(0,1)).imag();
	Cqu1_2313r_LNP = (Cqu1_uu_LNP*YucuL(1,2)*YuucL(0,2)).real();
	Cqu1_2313i_LNP = (Cqu1_uu_LNP*YucuL(1,2)*YuucL(0,2)).imag();
	Cqu1_2321r_LNP = (Cqu1_uu_LNP*YucuL(1,2)*YuucL(1,0)).real();
	Cqu1_2321i_LNP = (Cqu1_uu_LNP*YucuL(1,2)*YuucL(1,0)).imag();
	Cqu1_2322r_LNP = (Cqu1_u0_LNP*YucuL(1,2) + Cqu1_uu_LNP*YucuL(1,2)*YuucL(1,1)).real();
	Cqu1_2322i_LNP = (Cqu1_u0_LNP*YucuL(1,2) + Cqu1_uu_LNP*YucuL(1,2)*YuucL(1,1)).imag();
	Cqu1_2323r_LNP = (Cqu1_uu_LNP*YucuL(1,2)*YuucL(1,2)).real();
	Cqu1_2323i_LNP = (Cqu1_uu_LNP*YucuL(1,2)*YuucL(1,2)).imag();
	Cqu1_2331r_LNP = (Cqu1_uu_LNP*YucuL(1,2)*YuucL(2,0)).real();
	Cqu1_2331i_LNP = (Cqu1_uu_LNP*YucuL(1,2)*YuucL(2,0)).imag();
	Cqu1_2332r_LNP = (Cqu1_uu_LNP*YucuL(1,2)*YuucL(2,1)).real();
	Cqu1_2332i_LNP = (Cqu1_uu_LNP*YucuL(1,2)*YuucL(2,1)).imag();
	Cqu1_2333r_LNP = (Cqu1_u0_LNP*YucuL(1,2) + Cqu1_uu_LNP*YucuL(1,2)*YuucL(2,2)).real();
	Cqu1_2333i_LNP = (Cqu1_u0_LNP*YucuL(1,2) + Cqu1_uu_LNP*YucuL(1,2)*YuucL(2,2)).imag();
	Cqu1_3311r_LNP = (Cqu1_00_LNP + Cqu1_d0_LNP*YddL(2) + Cqu1_u0_LNP*YucuL(2,2) + Cqu1_0u_LNP*YuucL(0,0) + Cqu1_du_LNP*YddL(2)*YuucL(0,0) + Cqu1_uu_LNP*YucuL(2,2)*YuucL(0,0)).real();
	Cqu1_3312r_LNP = (Cqu1_0u_LNP*YuucL(0,1) + Cqu1_du_LNP*YddL(2)*YuucL(0,1) + Cqu1_uu_LNP*YucuL(2,2)*YuucL(0,1)).real();
	Cqu1_3312i_LNP = (Cqu1_0u_LNP*YuucL(0,1) + Cqu1_du_LNP*YddL(2)*YuucL(0,1) + Cqu1_uu_LNP*YucuL(2,2)*YuucL(0,1)).imag();
	Cqu1_3313r_LNP = (Cqu1_0u_LNP*YuucL(0,2) + Cqu1_du_LNP*YddL(2)*YuucL(0,2) + Cqu1_uu_LNP*YucuL(2,2)*YuucL(0,2)).real();
	Cqu1_3313i_LNP = (Cqu1_0u_LNP*YuucL(0,2) + Cqu1_du_LNP*YddL(2)*YuucL(0,2) + Cqu1_uu_LNP*YucuL(2,2)*YuucL(0,2)).imag();
	Cqu1_3322r_LNP = (Cqu1_00_LNP + Cqu1_d0_LNP*YddL(2) + Cqu1_u0_LNP*YucuL(2,2) + Cqu1_0u_LNP*YuucL(1,1) + Cqu1_du_LNP*YddL(2)*YuucL(1,1) + Cqu1_uu_LNP*YucuL(2,2)*YuucL(1,1)).real();
	Cqu1_3323r_LNP = (Cqu1_0u_LNP*YuucL(1,2) + Cqu1_du_LNP*YddL(2)*YuucL(1,2) + Cqu1_uu_LNP*YucuL(2,2)*YuucL(1,2)).real();
	Cqu1_3323i_LNP = (Cqu1_0u_LNP*YuucL(1,2) + Cqu1_du_LNP*YddL(2)*YuucL(1,2) + Cqu1_uu_LNP*YucuL(2,2)*YuucL(1,2)).imag();
	Cqu1_3333r_LNP = (Cqu1_00_LNP + Cqu1_d0_LNP*YddL(2) + Cqu1_u0_LNP*YucuL(2,2) + Cqu1_0u_LNP*YuucL(2,2) + Cqu1_du_LNP*YddL(2)*YuucL(2,2) + Cqu1_uu_LNP*YucuL(2,2)*YuucL(2,2)).real();

	Cqu8_1111r_LNP = (Cqu8_00_LNP + Cqu8_d0_LNP*YddL(0) + Cqu8_u0_LNP*YucuL(0,0) + Cqu8_0u_LNP*YuucL(0,0) + Cqu8_du_LNP*YddL(0)*YuucL(0,0) + Cqu8_uu_LNP*YucuL(0,0)*YuucL(0,0)).real();
	Cqu8_1112r_LNP = (Cqu8_0u_LNP*YuucL(0,1) + Cqu8_du_LNP*YddL(0)*YuucL(0,1) + Cqu8_uu_LNP*YucuL(0,0)*YuucL(0,1)).real();
	Cqu8_1112i_LNP = (Cqu8_0u_LNP*YuucL(0,1) + Cqu8_du_LNP*YddL(0)*YuucL(0,1) + Cqu8_uu_LNP*YucuL(0,0)*YuucL(0,1)).imag();
	Cqu8_1113r_LNP = (Cqu8_0u_LNP*YuucL(0,2) + Cqu8_du_LNP*YddL(0)*YuucL(0,2) + Cqu8_uu_LNP*YucuL(0,0)*YuucL(0,2)).real();
	Cqu8_1113i_LNP = (Cqu8_0u_LNP*YuucL(0,2) + Cqu8_du_LNP*YddL(0)*YuucL(0,2) + Cqu8_uu_LNP*YucuL(0,0)*YuucL(0,2)).imag();
	Cqu8_1122r_LNP = (Cqu8_00_LNP + Cqu8_d0_LNP*YddL(0) + Cqu8_u0_LNP*YucuL(0,0) + Cqu8_0u_LNP*YuucL(1,1) + Cqu8_du_LNP*YddL(0)*YuucL(1,1) + Cqu8_uu_LNP*YucuL(0,0)*YuucL(1,1)).real();
	Cqu8_1123r_LNP = (Cqu8_0u_LNP*YuucL(1,2) + Cqu8_du_LNP*YddL(0)*YuucL(1,2) + Cqu8_uu_LNP*YucuL(0,0)*YuucL(1,2)).real();
	Cqu8_1123i_LNP = (Cqu8_0u_LNP*YuucL(1,2) + Cqu8_du_LNP*YddL(0)*YuucL(1,2) + Cqu8_uu_LNP*YucuL(0,0)*YuucL(1,2)).imag();
	Cqu8_1133r_LNP = (Cqu8_00_LNP + Cqu8_d0_LNP*YddL(0) + Cqu8_u0_LNP*YucuL(0,0) + Cqu8_0u_LNP*YuucL(2,2) + Cqu8_du_LNP*YddL(0)*YuucL(2,2) + Cqu8_uu_LNP*YucuL(0,0)*YuucL(2,2)).real();
	Cqu8_1211r_LNP = (Cqu8_u0_LNP*YucuL(0,1) + Cqu8_uu_LNP*YucuL(0,1)*YuucL(0,0)).real();
	Cqu8_1211i_LNP = (Cqu8_u0_LNP*YucuL(0,1) + Cqu8_uu_LNP*YucuL(0,1)*YuucL(0,0)).imag();
	Cqu8_1212r_LNP = (Cqu8_uu_LNP*YucuL(0,1)*YuucL(0,1)).real();
	Cqu8_1212i_LNP = (Cqu8_uu_LNP*YucuL(0,1)*YuucL(0,1)).imag();
	Cqu8_1213r_LNP = (Cqu8_uu_LNP*YucuL(0,1)*YuucL(0,2)).real();
	Cqu8_1213i_LNP = (Cqu8_uu_LNP*YucuL(0,1)*YuucL(0,2)).imag();
	Cqu8_1221r_LNP = (Cqu8_uu_LNP*YucuL(0,1)*YuucL(1,0)).real();
	Cqu8_1221i_LNP = (Cqu8_uu_LNP*YucuL(0,1)*YuucL(1,0)).imag();
	Cqu8_1222r_LNP = (Cqu8_u0_LNP*YucuL(0,1) + Cqu8_uu_LNP*YucuL(0,1)*YuucL(1,1)).real();
	Cqu8_1222i_LNP = (Cqu8_u0_LNP*YucuL(0,1) + Cqu8_uu_LNP*YucuL(0,1)*YuucL(1,1)).imag();
	Cqu8_1223r_LNP = (Cqu8_uu_LNP*YucuL(0,1)*YuucL(1,2)).real();
	Cqu8_1223i_LNP = (Cqu8_uu_LNP*YucuL(0,1)*YuucL(1,2)).imag();
	Cqu8_1231r_LNP = (Cqu8_uu_LNP*YucuL(0,1)*YuucL(2,0)).real();
	Cqu8_1231i_LNP = (Cqu8_uu_LNP*YucuL(0,1)*YuucL(2,0)).imag();
	Cqu8_1232r_LNP = (Cqu8_uu_LNP*YucuL(0,1)*YuucL(2,1)).real();
	Cqu8_1232i_LNP = (Cqu8_uu_LNP*YucuL(0,1)*YuucL(2,1)).imag();
	Cqu8_1233r_LNP = (Cqu8_u0_LNP*YucuL(0,1) + Cqu8_uu_LNP*YucuL(0,1)*YuucL(2,2)).real();
	Cqu8_1233i_LNP = (Cqu8_u0_LNP*YucuL(0,1) + Cqu8_uu_LNP*YucuL(0,1)*YuucL(2,2)).imag();
	Cqu8_1311r_LNP = (Cqu8_u0_LNP*YucuL(0,2) + Cqu8_uu_LNP*YucuL(0,2)*YuucL(0,0)).real();
	Cqu8_1311i_LNP = (Cqu8_u0_LNP*YucuL(0,2) + Cqu8_uu_LNP*YucuL(0,2)*YuucL(0,0)).imag();
	Cqu8_1312r_LNP = (Cqu8_uu_LNP*YucuL(0,2)*YuucL(0,1)).real();
	Cqu8_1312i_LNP = (Cqu8_uu_LNP*YucuL(0,2)*YuucL(0,1)).imag();
	Cqu8_1313r_LNP = (Cqu8_uu_LNP*YucuL(0,2)*YuucL(0,2)).real();
	Cqu8_1313i_LNP = (Cqu8_uu_LNP*YucuL(0,2)*YuucL(0,2)).imag();
	Cqu8_1321r_LNP = (Cqu8_uu_LNP*YucuL(0,2)*YuucL(1,0)).real();
	Cqu8_1321i_LNP = (Cqu8_uu_LNP*YucuL(0,2)*YuucL(1,0)).imag();
	Cqu8_1322r_LNP = (Cqu8_u0_LNP*YucuL(0,2) + Cqu8_uu_LNP*YucuL(0,2)*YuucL(1,1)).real();
	Cqu8_1322i_LNP = (Cqu8_u0_LNP*YucuL(0,2) + Cqu8_uu_LNP*YucuL(0,2)*YuucL(1,1)).imag();
	Cqu8_1323r_LNP = (Cqu8_uu_LNP*YucuL(0,2)*YuucL(1,2)).real();
	Cqu8_1323i_LNP = (Cqu8_uu_LNP*YucuL(0,2)*YuucL(1,2)).imag();
	Cqu8_1331r_LNP = (Cqu8_uu_LNP*YucuL(0,2)*YuucL(2,0)).real();
	Cqu8_1331i_LNP = (Cqu8_uu_LNP*YucuL(0,2)*YuucL(2,0)).imag();
	Cqu8_1332r_LNP = (Cqu8_uu_LNP*YucuL(0,2)*YuucL(2,1)).real();
	Cqu8_1332i_LNP = (Cqu8_uu_LNP*YucuL(0,2)*YuucL(2,1)).imag();
	Cqu8_1333r_LNP = (Cqu8_u0_LNP*YucuL(0,2) + Cqu8_uu_LNP*YucuL(0,2)*YuucL(2,2)).real();
	Cqu8_1333i_LNP = (Cqu8_u0_LNP*YucuL(0,2) + Cqu8_uu_LNP*YucuL(0,2)*YuucL(2,2)).imag();
	Cqu8_2211r_LNP = (Cqu8_00_LNP + Cqu8_d0_LNP*YddL(1) + Cqu8_u0_LNP*YucuL(1,1) + Cqu8_0u_LNP*YuucL(0,0) + Cqu8_du_LNP*YddL(1)*YuucL(0,0) + Cqu8_uu_LNP*YucuL(1,1)*YuucL(0,0)).real();
	Cqu8_2212r_LNP = (Cqu8_0u_LNP*YuucL(0,1) + Cqu8_du_LNP*YddL(1)*YuucL(0,1) + Cqu8_uu_LNP*YucuL(1,1)*YuucL(0,1)).real();
	Cqu8_2212i_LNP = (Cqu8_0u_LNP*YuucL(0,1) + Cqu8_du_LNP*YddL(1)*YuucL(0,1) + Cqu8_uu_LNP*YucuL(1,1)*YuucL(0,1)).imag();
	Cqu8_2213r_LNP = (Cqu8_0u_LNP*YuucL(0,2) + Cqu8_du_LNP*YddL(1)*YuucL(0,2) + Cqu8_uu_LNP*YucuL(1,1)*YuucL(0,2)).real();
	Cqu8_2213i_LNP = (Cqu8_0u_LNP*YuucL(0,2) + Cqu8_du_LNP*YddL(1)*YuucL(0,2) + Cqu8_uu_LNP*YucuL(1,1)*YuucL(0,2)).imag();
	Cqu8_2222r_LNP = (Cqu8_00_LNP + Cqu8_d0_LNP*YddL(1) + Cqu8_u0_LNP*YucuL(1,1) + Cqu8_0u_LNP*YuucL(1,1) + Cqu8_du_LNP*YddL(1)*YuucL(1,1) + Cqu8_uu_LNP*YucuL(1,1)*YuucL(1,1)).real();
	Cqu8_2223r_LNP = (Cqu8_0u_LNP*YuucL(1,2) + Cqu8_du_LNP*YddL(1)*YuucL(1,2) + Cqu8_uu_LNP*YucuL(1,1)*YuucL(1,2)).real();
	Cqu8_2223i_LNP = (Cqu8_0u_LNP*YuucL(1,2) + Cqu8_du_LNP*YddL(1)*YuucL(1,2) + Cqu8_uu_LNP*YucuL(1,1)*YuucL(1,2)).imag();
	Cqu8_2233r_LNP = (Cqu8_00_LNP + Cqu8_d0_LNP*YddL(1) + Cqu8_u0_LNP*YucuL(1,1) + Cqu8_0u_LNP*YuucL(2,2) + Cqu8_du_LNP*YddL(1)*YuucL(2,2) + Cqu8_uu_LNP*YucuL(1,1)*YuucL(2,2)).real();
	Cqu8_2311r_LNP = (Cqu8_u0_LNP*YucuL(1,2) + Cqu8_uu_LNP*YucuL(1,2)*YuucL(0,0)).real();
	Cqu8_2311i_LNP = (Cqu8_u0_LNP*YucuL(1,2) + Cqu8_uu_LNP*YucuL(1,2)*YuucL(0,0)).imag();
	Cqu8_2312r_LNP = (Cqu8_uu_LNP*YucuL(1,2)*YuucL(0,1)).real();
	Cqu8_2312i_LNP = (Cqu8_uu_LNP*YucuL(1,2)*YuucL(0,1)).imag();
	Cqu8_2313r_LNP = (Cqu8_uu_LNP*YucuL(1,2)*YuucL(0,2)).real();
	Cqu8_2313i_LNP = (Cqu8_uu_LNP*YucuL(1,2)*YuucL(0,2)).imag();
	Cqu8_2321r_LNP = (Cqu8_uu_LNP*YucuL(1,2)*YuucL(1,0)).real();
	Cqu8_2321i_LNP = (Cqu8_uu_LNP*YucuL(1,2)*YuucL(1,0)).imag();
	Cqu8_2322r_LNP = (Cqu8_u0_LNP*YucuL(1,2) + Cqu8_uu_LNP*YucuL(1,2)*YuucL(1,1)).real();
	Cqu8_2322i_LNP = (Cqu8_u0_LNP*YucuL(1,2) + Cqu8_uu_LNP*YucuL(1,2)*YuucL(1,1)).imag();
	Cqu8_2323r_LNP = (Cqu8_uu_LNP*YucuL(1,2)*YuucL(1,2)).real();
	Cqu8_2323i_LNP = (Cqu8_uu_LNP*YucuL(1,2)*YuucL(1,2)).imag();
	Cqu8_2331r_LNP = (Cqu8_uu_LNP*YucuL(1,2)*YuucL(2,0)).real();
	Cqu8_2331i_LNP = (Cqu8_uu_LNP*YucuL(1,2)*YuucL(2,0)).imag();
	Cqu8_2332r_LNP = (Cqu8_uu_LNP*YucuL(1,2)*YuucL(2,1)).real();
	Cqu8_2332i_LNP = (Cqu8_uu_LNP*YucuL(1,2)*YuucL(2,1)).imag();
	Cqu8_2333r_LNP = (Cqu8_u0_LNP*YucuL(1,2) + Cqu8_uu_LNP*YucuL(1,2)*YuucL(2,2)).real();
	Cqu8_2333i_LNP = (Cqu8_u0_LNP*YucuL(1,2) + Cqu8_uu_LNP*YucuL(1,2)*YuucL(2,2)).imag();
	Cqu8_3311r_LNP = (Cqu8_00_LNP + Cqu8_d0_LNP*YddL(2) + Cqu8_u0_LNP*YucuL(2,2) + Cqu8_0u_LNP*YuucL(0,0) + Cqu8_du_LNP*YddL(2)*YuucL(0,0) + Cqu8_uu_LNP*YucuL(2,2)*YuucL(0,0)).real();
	Cqu8_3312r_LNP = (Cqu8_0u_LNP*YuucL(0,1) + Cqu8_du_LNP*YddL(2)*YuucL(0,1) + Cqu8_uu_LNP*YucuL(2,2)*YuucL(0,1)).real();
	Cqu8_3312i_LNP = (Cqu8_0u_LNP*YuucL(0,1) + Cqu8_du_LNP*YddL(2)*YuucL(0,1) + Cqu8_uu_LNP*YucuL(2,2)*YuucL(0,1)).imag();
	Cqu8_3313r_LNP = (Cqu8_0u_LNP*YuucL(0,2) + Cqu8_du_LNP*YddL(2)*YuucL(0,2) + Cqu8_uu_LNP*YucuL(2,2)*YuucL(0,2)).real();
	Cqu8_3313i_LNP = (Cqu8_0u_LNP*YuucL(0,2) + Cqu8_du_LNP*YddL(2)*YuucL(0,2) + Cqu8_uu_LNP*YucuL(2,2)*YuucL(0,2)).imag();
	Cqu8_3322r_LNP = (Cqu8_00_LNP + Cqu8_d0_LNP*YddL(2) + Cqu8_u0_LNP*YucuL(2,2) + Cqu8_0u_LNP*YuucL(1,1) + Cqu8_du_LNP*YddL(2)*YuucL(1,1) + Cqu8_uu_LNP*YucuL(2,2)*YuucL(1,1)).real();
	Cqu8_3323r_LNP = (Cqu8_0u_LNP*YuucL(1,2) + Cqu8_du_LNP*YddL(2)*YuucL(1,2) + Cqu8_uu_LNP*YucuL(2,2)*YuucL(1,2)).real();
	Cqu8_3323i_LNP = (Cqu8_0u_LNP*YuucL(1,2) + Cqu8_du_LNP*YddL(2)*YuucL(1,2) + Cqu8_uu_LNP*YucuL(2,2)*YuucL(1,2)).imag();
	Cqu8_3333r_LNP = (Cqu8_00_LNP + Cqu8_d0_LNP*YddL(2) + Cqu8_u0_LNP*YucuL(2,2) + Cqu8_0u_LNP*YuucL(2,2) + Cqu8_du_LNP*YddL(2)*YuucL(2,2) + Cqu8_uu_LNP*YucuL(2,2)*YuucL(2,2)).real();

	Cqd1_1111r_LNP = (Cqd1_00_LNP + Cqd1_0d_LNP*YddL(0) + Cqd1_d0_LNP*YddL(0) + Cqd1_dd_LNP*YddL(0)*YddL(0) + Cqd1_u0_LNP*YucuL(0,0) + Cqd1_ud_LNP*YddL(0)*YucuL(0,0)).real();
	Cqd1_1122r_LNP = (Cqd1_00_LNP + Cqd1_d0_LNP*YddL(0) + Cqd1_0d_LNP*YddL(1) + Cqd1_dd_LNP*YddL(0)*YddL(1) + Cqd1_u0_LNP*YucuL(0,0) + Cqd1_ud_LNP*YddL(1)*YucuL(0,0)).real();
	Cqd1_1133r_LNP = (Cqd1_00_LNP + Cqd1_d0_LNP*YddL(0) + Cqd1_0d_LNP*YddL(2) + Cqd1_dd_LNP*YddL(0)*YddL(2) + Cqd1_u0_LNP*YucuL(0,0) + Cqd1_ud_LNP*YddL(2)*YucuL(0,0)).real();
	Cqd1_1211r_LNP = (Cqd1_u0_LNP*YucuL(0,1) + Cqd1_ud_LNP*YddL(0)*YucuL(0,1)).real();
	Cqd1_1211i_LNP = (Cqd1_u0_LNP*YucuL(0,1) + Cqd1_ud_LNP*YddL(0)*YucuL(0,1)).imag();
	Cqd1_1222r_LNP = (Cqd1_u0_LNP*YucuL(0,1) + Cqd1_ud_LNP*YddL(1)*YucuL(0,1)).real();
	Cqd1_1222i_LNP = (Cqd1_u0_LNP*YucuL(0,1) + Cqd1_ud_LNP*YddL(1)*YucuL(0,1)).imag();
	Cqd1_1233r_LNP = (Cqd1_u0_LNP*YucuL(0,1) + Cqd1_ud_LNP*YddL(2)*YucuL(0,1)).real();
	Cqd1_1233i_LNP = (Cqd1_u0_LNP*YucuL(0,1) + Cqd1_ud_LNP*YddL(2)*YucuL(0,1)).imag();
	Cqd1_1311r_LNP = (Cqd1_u0_LNP*YucuL(0,2) + Cqd1_ud_LNP*YddL(0)*YucuL(0,2)).real();
	Cqd1_1311i_LNP = (Cqd1_u0_LNP*YucuL(0,2) + Cqd1_ud_LNP*YddL(0)*YucuL(0,2)).imag();
	Cqd1_1322r_LNP = (Cqd1_u0_LNP*YucuL(0,2) + Cqd1_ud_LNP*YddL(1)*YucuL(0,2)).real();
	Cqd1_1322i_LNP = (Cqd1_u0_LNP*YucuL(0,2) + Cqd1_ud_LNP*YddL(1)*YucuL(0,2)).imag();
	Cqd1_1333r_LNP = (Cqd1_u0_LNP*YucuL(0,2) + Cqd1_ud_LNP*YddL(2)*YucuL(0,2)).real();
	Cqd1_1333i_LNP = (Cqd1_u0_LNP*YucuL(0,2) + Cqd1_ud_LNP*YddL(2)*YucuL(0,2)).imag();
	Cqd1_2211r_LNP = (Cqd1_00_LNP + Cqd1_0d_LNP*YddL(0) + Cqd1_d0_LNP*YddL(1) + Cqd1_dd_LNP*YddL(0)*YddL(1) + Cqd1_u0_LNP*YucuL(1,1) + Cqd1_ud_LNP*YddL(0)*YucuL(1,1)).real();
	Cqd1_2222r_LNP = (Cqd1_00_LNP + Cqd1_0d_LNP*YddL(1) + Cqd1_d0_LNP*YddL(1) + Cqd1_dd_LNP*YddL(1)*YddL(1) + Cqd1_u0_LNP*YucuL(1,1) + Cqd1_ud_LNP*YddL(1)*YucuL(1,1)).real();
	Cqd1_2233r_LNP = (Cqd1_00_LNP + Cqd1_d0_LNP*YddL(1) + Cqd1_0d_LNP*YddL(2) + Cqd1_dd_LNP*YddL(1)*YddL(2) + Cqd1_u0_LNP*YucuL(1,1) + Cqd1_ud_LNP*YddL(2)*YucuL(1,1)).real();
	Cqd1_2311r_LNP = (Cqd1_u0_LNP*YucuL(1,2) + Cqd1_ud_LNP*YddL(0)*YucuL(1,2)).real();
	Cqd1_2311i_LNP = (Cqd1_u0_LNP*YucuL(1,2) + Cqd1_ud_LNP*YddL(0)*YucuL(1,2)).imag();
	Cqd1_2322r_LNP = (Cqd1_u0_LNP*YucuL(1,2) + Cqd1_ud_LNP*YddL(1)*YucuL(1,2)).real();
	Cqd1_2322i_LNP = (Cqd1_u0_LNP*YucuL(1,2) + Cqd1_ud_LNP*YddL(1)*YucuL(1,2)).imag();
	Cqd1_2333r_LNP = (Cqd1_u0_LNP*YucuL(1,2) + Cqd1_ud_LNP*YddL(2)*YucuL(1,2)).real();
	Cqd1_2333i_LNP = (Cqd1_u0_LNP*YucuL(1,2) + Cqd1_ud_LNP*YddL(2)*YucuL(1,2)).imag();
	Cqd1_3311r_LNP = (Cqd1_00_LNP + Cqd1_0d_LNP*YddL(0) + Cqd1_d0_LNP*YddL(2) + Cqd1_dd_LNP*YddL(0)*YddL(2) + Cqd1_u0_LNP*YucuL(2,2) + Cqd1_ud_LNP*YddL(0)*YucuL(2,2)).real();
	Cqd1_3322r_LNP = (Cqd1_00_LNP + Cqd1_0d_LNP*YddL(1) + Cqd1_d0_LNP*YddL(2) + Cqd1_dd_LNP*YddL(1)*YddL(2) + Cqd1_u0_LNP*YucuL(2,2) + Cqd1_ud_LNP*YddL(1)*YucuL(2,2)).real();
	Cqd1_3333r_LNP = (Cqd1_00_LNP + Cqd1_0d_LNP*YddL(2) + Cqd1_d0_LNP*YddL(2) + Cqd1_dd_LNP*YddL(2)*YddL(2) + Cqd1_u0_LNP*YucuL(2,2) + Cqd1_ud_LNP*YddL(2)*YucuL(2,2)).real();

	Cqd8_1111r_LNP = (Cqd8_00_LNP + Cqd8_0d_LNP*YddL(0) + Cqd8_d0_LNP*YddL(0) + Cqd8_dd_LNP*YddL(0)*YddL(0) + Cqd8_u0_LNP*YucuL(0,0) + Cqd8_ud_LNP*YddL(0)*YucuL(0,0)).real();
	Cqd8_1122r_LNP = (Cqd8_00_LNP + Cqd8_d0_LNP*YddL(0) + Cqd8_0d_LNP*YddL(1) + Cqd8_dd_LNP*YddL(0)*YddL(1) + Cqd8_u0_LNP*YucuL(0,0) + Cqd8_ud_LNP*YddL(1)*YucuL(0,0)).real();
	Cqd8_1133r_LNP = (Cqd8_00_LNP + Cqd8_d0_LNP*YddL(0) + Cqd8_0d_LNP*YddL(2) + Cqd8_dd_LNP*YddL(0)*YddL(2) + Cqd8_u0_LNP*YucuL(0,0) + Cqd8_ud_LNP*YddL(2)*YucuL(0,0)).real();
	Cqd8_1211r_LNP = (Cqd8_u0_LNP*YucuL(0,1) + Cqd8_ud_LNP*YddL(0)*YucuL(0,1)).real();
	Cqd8_1211i_LNP = (Cqd8_u0_LNP*YucuL(0,1) + Cqd8_ud_LNP*YddL(0)*YucuL(0,1)).imag();
	Cqd8_1222r_LNP = (Cqd8_u0_LNP*YucuL(0,1) + Cqd8_ud_LNP*YddL(1)*YucuL(0,1)).real();
	Cqd8_1222i_LNP = (Cqd8_u0_LNP*YucuL(0,1) + Cqd8_ud_LNP*YddL(1)*YucuL(0,1)).imag();
	Cqd8_1233r_LNP = (Cqd8_u0_LNP*YucuL(0,1) + Cqd8_ud_LNP*YddL(2)*YucuL(0,1)).real();
	Cqd8_1233i_LNP = (Cqd8_u0_LNP*YucuL(0,1) + Cqd8_ud_LNP*YddL(2)*YucuL(0,1)).imag();
	Cqd8_1311r_LNP = (Cqd8_u0_LNP*YucuL(0,2) + Cqd8_ud_LNP*YddL(0)*YucuL(0,2)).real();
	Cqd8_1311i_LNP = (Cqd8_u0_LNP*YucuL(0,2) + Cqd8_ud_LNP*YddL(0)*YucuL(0,2)).imag();
	Cqd8_1322r_LNP = (Cqd8_u0_LNP*YucuL(0,2) + Cqd8_ud_LNP*YddL(1)*YucuL(0,2)).real();
	Cqd8_1322i_LNP = (Cqd8_u0_LNP*YucuL(0,2) + Cqd8_ud_LNP*YddL(1)*YucuL(0,2)).imag();
	Cqd8_1333r_LNP = (Cqd8_u0_LNP*YucuL(0,2) + Cqd8_ud_LNP*YddL(2)*YucuL(0,2)).real();
	Cqd8_1333i_LNP = (Cqd8_u0_LNP*YucuL(0,2) + Cqd8_ud_LNP*YddL(2)*YucuL(0,2)).imag();
	Cqd8_2211r_LNP = (Cqd8_00_LNP + Cqd8_0d_LNP*YddL(0) + Cqd8_d0_LNP*YddL(1) + Cqd8_dd_LNP*YddL(0)*YddL(1) + Cqd8_u0_LNP*YucuL(1,1) + Cqd8_ud_LNP*YddL(0)*YucuL(1,1)).real();
	Cqd8_2222r_LNP = (Cqd8_00_LNP + Cqd8_0d_LNP*YddL(1) + Cqd8_d0_LNP*YddL(1) + Cqd8_dd_LNP*YddL(1)*YddL(1) + Cqd8_u0_LNP*YucuL(1,1) + Cqd8_ud_LNP*YddL(1)*YucuL(1,1)).real();
	Cqd8_2233r_LNP = (Cqd8_00_LNP + Cqd8_d0_LNP*YddL(1) + Cqd8_0d_LNP*YddL(2) + Cqd8_dd_LNP*YddL(1)*YddL(2) + Cqd8_u0_LNP*YucuL(1,1) + Cqd8_ud_LNP*YddL(2)*YucuL(1,1)).real();
	Cqd8_2311r_LNP = (Cqd8_u0_LNP*YucuL(1,2) + Cqd8_ud_LNP*YddL(0)*YucuL(1,2)).real();
	Cqd8_2311i_LNP = (Cqd8_u0_LNP*YucuL(1,2) + Cqd8_ud_LNP*YddL(0)*YucuL(1,2)).imag();
	Cqd8_2322r_LNP = (Cqd8_u0_LNP*YucuL(1,2) + Cqd8_ud_LNP*YddL(1)*YucuL(1,2)).real();
	Cqd8_2322i_LNP = (Cqd8_u0_LNP*YucuL(1,2) + Cqd8_ud_LNP*YddL(1)*YucuL(1,2)).imag();
	Cqd8_2333r_LNP = (Cqd8_u0_LNP*YucuL(1,2) + Cqd8_ud_LNP*YddL(2)*YucuL(1,2)).real();
	Cqd8_2333i_LNP = (Cqd8_u0_LNP*YucuL(1,2) + Cqd8_ud_LNP*YddL(2)*YucuL(1,2)).imag();
	Cqd8_3311r_LNP = (Cqd8_00_LNP + Cqd8_0d_LNP*YddL(0) + Cqd8_d0_LNP*YddL(2) + Cqd8_dd_LNP*YddL(0)*YddL(2) + Cqd8_u0_LNP*YucuL(2,2) + Cqd8_ud_LNP*YddL(0)*YucuL(2,2)).real();
	Cqd8_3322r_LNP = (Cqd8_00_LNP + Cqd8_0d_LNP*YddL(1) + Cqd8_d0_LNP*YddL(2) + Cqd8_dd_LNP*YddL(1)*YddL(2) + Cqd8_u0_LNP*YucuL(2,2) + Cqd8_ud_LNP*YddL(1)*YucuL(2,2)).real();
	Cqd8_3333r_LNP = (Cqd8_00_LNP + Cqd8_0d_LNP*YddL(2) + Cqd8_d0_LNP*YddL(2) + Cqd8_dd_LNP*YddL(2)*YddL(2) + Cqd8_u0_LNP*YucuL(2,2) + Cqd8_ud_LNP*YddL(2)*YucuL(2,2)).real();

	Cquqd1_1111r_LNP = (2*Cquqd1_00_LNP*YdL(0)*YuL(0,0)).real();
	Cquqd1_1111i_LNP = (2*Cquqd1_00_LNP*YdL(0)*YuL(0,0)).imag();
	Cquqd1_1121r_LNP = (Cquqd1_00_LNP*YdL(0)*YuL(1,0)).real();
	Cquqd1_1121i_LNP = (Cquqd1_00_LNP*YdL(0)*YuL(1,0)).imag();
	Cquqd1_1122r_LNP = (Cquqd1_00_LNP*YdL(1)*YuL(0,0)).real();
	Cquqd1_1122i_LNP = (Cquqd1_00_LNP*YdL(1)*YuL(0,0)).imag();
	Cquqd1_1131r_LNP = (Cquqd1_00_LNP*YdL(0)*YuL(2,0)).real();
	Cquqd1_1131i_LNP = (Cquqd1_00_LNP*YdL(0)*YuL(2,0)).imag();
	Cquqd1_1133r_LNP = (Cquqd1_00_LNP*YdL(2)*YuL(0,0)).real();
	Cquqd1_1133i_LNP = (Cquqd1_00_LNP*YdL(2)*YuL(0,0)).imag();
	Cquqd1_1211r_LNP = (2*Cquqd1_00_LNP*YdL(0)*YuL(0,1)).real();
	Cquqd1_1211i_LNP = (2*Cquqd1_00_LNP*YdL(0)*YuL(0,1)).imag();
	Cquqd1_1221r_LNP = (Cquqd1_00_LNP*YdL(0)*YuL(1,1)).real();
	Cquqd1_1221i_LNP = (Cquqd1_00_LNP*YdL(0)*YuL(1,1)).imag();
	Cquqd1_1222r_LNP = (Cquqd1_00_LNP*YdL(1)*YuL(0,1)).real();
	Cquqd1_1222i_LNP = (Cquqd1_00_LNP*YdL(1)*YuL(0,1)).imag();
	Cquqd1_1231r_LNP = (Cquqd1_00_LNP*YdL(0)*YuL(2,1)).real();
	Cquqd1_1231i_LNP = (Cquqd1_00_LNP*YdL(0)*YuL(2,1)).imag();
	Cquqd1_1233r_LNP = (Cquqd1_00_LNP*YdL(2)*YuL(0,1)).real();
	Cquqd1_1233i_LNP = (Cquqd1_00_LNP*YdL(2)*YuL(0,1)).imag();
	Cquqd1_1311r_LNP = (2*Cquqd1_00_LNP*YdL(0)*YuL(0,2)).real();
	Cquqd1_1311i_LNP = (2*Cquqd1_00_LNP*YdL(0)*YuL(0,2)).imag();
	Cquqd1_1321r_LNP = (Cquqd1_00_LNP*YdL(0)*YuL(1,2)).real();
	Cquqd1_1321i_LNP = (Cquqd1_00_LNP*YdL(0)*YuL(1,2)).imag();
	Cquqd1_1322r_LNP = (Cquqd1_00_LNP*YdL(1)*YuL(0,2)).real();
	Cquqd1_1322i_LNP = (Cquqd1_00_LNP*YdL(1)*YuL(0,2)).imag();
	Cquqd1_1331r_LNP = (Cquqd1_00_LNP*YdL(0)*YuL(2,2)).real();
	Cquqd1_1331i_LNP = (Cquqd1_00_LNP*YdL(0)*YuL(2,2)).imag();
	Cquqd1_1333r_LNP = (Cquqd1_00_LNP*YdL(2)*YuL(0,2)).real();
	Cquqd1_1333i_LNP = (Cquqd1_00_LNP*YdL(2)*YuL(0,2)).imag();
	Cquqd1_2111r_LNP = (Cquqd1_00_LNP*YdL(0)*YuL(1,0)).real();
	Cquqd1_2111i_LNP = (Cquqd1_00_LNP*YdL(0)*YuL(1,0)).imag();
	Cquqd1_2112r_LNP = (Cquqd1_00_LNP*YdL(1)*YuL(0,0)).real();
	Cquqd1_2112i_LNP = (Cquqd1_00_LNP*YdL(1)*YuL(0,0)).imag();
	Cquqd1_2122r_LNP = (2*Cquqd1_00_LNP*YdL(1)*YuL(1,0)).real();
	Cquqd1_2122i_LNP = (2*Cquqd1_00_LNP*YdL(1)*YuL(1,0)).imag();
	Cquqd1_2132r_LNP = (Cquqd1_00_LNP*YdL(1)*YuL(2,0)).real();
	Cquqd1_2132i_LNP = (Cquqd1_00_LNP*YdL(1)*YuL(2,0)).imag();
	Cquqd1_2133r_LNP = (Cquqd1_00_LNP*YdL(2)*YuL(1,0)).real();
	Cquqd1_2133i_LNP = (Cquqd1_00_LNP*YdL(2)*YuL(1,0)).imag();
	Cquqd1_2211r_LNP = (Cquqd1_00_LNP*YdL(0)*YuL(1,1)).real();
	Cquqd1_2211i_LNP = (Cquqd1_00_LNP*YdL(0)*YuL(1,1)).imag();
	Cquqd1_2212r_LNP = (Cquqd1_00_LNP*YdL(1)*YuL(0,1)).real();
	Cquqd1_2212i_LNP = (Cquqd1_00_LNP*YdL(1)*YuL(0,1)).imag();
	Cquqd1_2222r_LNP = (2*Cquqd1_00_LNP*YdL(1)*YuL(1,1)).real();
	Cquqd1_2222i_LNP = (2*Cquqd1_00_LNP*YdL(1)*YuL(1,1)).imag();
	Cquqd1_2232r_LNP = (Cquqd1_00_LNP*YdL(1)*YuL(2,1)).real();
	Cquqd1_2232i_LNP = (Cquqd1_00_LNP*YdL(1)*YuL(2,1)).imag();
	Cquqd1_2233r_LNP = (Cquqd1_00_LNP*YdL(2)*YuL(1,1)).real();
	Cquqd1_2233i_LNP = (Cquqd1_00_LNP*YdL(2)*YuL(1,1)).imag();
	Cquqd1_2311r_LNP = (Cquqd1_00_LNP*YdL(0)*YuL(1,2)).real();
	Cquqd1_2311i_LNP = (Cquqd1_00_LNP*YdL(0)*YuL(1,2)).imag();
	Cquqd1_2312r_LNP = (Cquqd1_00_LNP*YdL(1)*YuL(0,2)).real();
	Cquqd1_2312i_LNP = (Cquqd1_00_LNP*YdL(1)*YuL(0,2)).imag();
	Cquqd1_2322r_LNP = (2*Cquqd1_00_LNP*YdL(1)*YuL(1,2)).real();
	Cquqd1_2322i_LNP = (2*Cquqd1_00_LNP*YdL(1)*YuL(1,2)).imag();
	Cquqd1_2332r_LNP = (Cquqd1_00_LNP*YdL(1)*YuL(2,2)).real();
	Cquqd1_2332i_LNP = (Cquqd1_00_LNP*YdL(1)*YuL(2,2)).imag();
	Cquqd1_2333r_LNP = (Cquqd1_00_LNP*YdL(2)*YuL(1,2)).real();
	Cquqd1_2333i_LNP = (Cquqd1_00_LNP*YdL(2)*YuL(1,2)).imag();
	Cquqd1_3111r_LNP = (Cquqd1_00_LNP*YdL(0)*YuL(2,0)).real();
	Cquqd1_3111i_LNP = (Cquqd1_00_LNP*YdL(0)*YuL(2,0)).imag();
	Cquqd1_3113r_LNP = (Cquqd1_00_LNP*YdL(2)*YuL(0,0)).real();
	Cquqd1_3113i_LNP = (Cquqd1_00_LNP*YdL(2)*YuL(0,0)).imag();
	Cquqd1_3122r_LNP = (Cquqd1_00_LNP*YdL(1)*YuL(2,0)).real();
	Cquqd1_3122i_LNP = (Cquqd1_00_LNP*YdL(1)*YuL(2,0)).imag();
	Cquqd1_3123r_LNP = (Cquqd1_00_LNP*YdL(2)*YuL(1,0)).real();
	Cquqd1_3123i_LNP = (Cquqd1_00_LNP*YdL(2)*YuL(1,0)).imag();
	Cquqd1_3133r_LNP = (2*Cquqd1_00_LNP*YdL(2)*YuL(2,0)).real();
	Cquqd1_3133i_LNP = (2*Cquqd1_00_LNP*YdL(2)*YuL(2,0)).imag();
	Cquqd1_3211r_LNP = (Cquqd1_00_LNP*YdL(0)*YuL(2,1)).real();
	Cquqd1_3211i_LNP = (Cquqd1_00_LNP*YdL(0)*YuL(2,1)).imag();
	Cquqd1_3213r_LNP = (Cquqd1_00_LNP*YdL(2)*YuL(0,1)).real();
	Cquqd1_3213i_LNP = (Cquqd1_00_LNP*YdL(2)*YuL(0,1)).imag();
	Cquqd1_3222r_LNP = (Cquqd1_00_LNP*YdL(1)*YuL(2,1)).real();
	Cquqd1_3222i_LNP = (Cquqd1_00_LNP*YdL(1)*YuL(2,1)).imag();
	Cquqd1_3223r_LNP = (Cquqd1_00_LNP*YdL(2)*YuL(1,1)).real();
	Cquqd1_3223i_LNP = (Cquqd1_00_LNP*YdL(2)*YuL(1,1)).imag();
	Cquqd1_3233r_LNP = (2*Cquqd1_00_LNP*YdL(2)*YuL(2,1)).real();
	Cquqd1_3233i_LNP = (2*Cquqd1_00_LNP*YdL(2)*YuL(2,1)).imag();
	Cquqd1_3311r_LNP = (Cquqd1_00_LNP*YdL(0)*YuL(2,2)).real();
	Cquqd1_3311i_LNP = (Cquqd1_00_LNP*YdL(0)*YuL(2,2)).imag();
	Cquqd1_3313r_LNP = (Cquqd1_00_LNP*YdL(2)*YuL(0,2)).real();
	Cquqd1_3313i_LNP = (Cquqd1_00_LNP*YdL(2)*YuL(0,2)).imag();
	Cquqd1_3322r_LNP = (Cquqd1_00_LNP*YdL(1)*YuL(2,2)).real();
	Cquqd1_3322i_LNP = (Cquqd1_00_LNP*YdL(1)*YuL(2,2)).imag();
	Cquqd1_3323r_LNP = (Cquqd1_00_LNP*YdL(2)*YuL(1,2)).real();
	Cquqd1_3323i_LNP = (Cquqd1_00_LNP*YdL(2)*YuL(1,2)).imag();
	Cquqd1_3333r_LNP = (2*Cquqd1_00_LNP*YdL(2)*YuL(2,2)).real();
	Cquqd1_3333i_LNP = (2*Cquqd1_00_LNP*YdL(2)*YuL(2,2)).imag();

	Cquqd8_1111r_LNP = (2*Cquqd8_00_LNP*YdL(0)*YuL(0,0)).real();
	Cquqd8_1111i_LNP = (2*Cquqd8_00_LNP*YdL(0)*YuL(0,0)).imag();
	Cquqd8_1121r_LNP = (Cquqd8_00_LNP*YdL(0)*YuL(1,0)).real();
	Cquqd8_1121i_LNP = (Cquqd8_00_LNP*YdL(0)*YuL(1,0)).imag();
	Cquqd8_1122r_LNP = (Cquqd8_00_LNP*YdL(1)*YuL(0,0)).real();
	Cquqd8_1122i_LNP = (Cquqd8_00_LNP*YdL(1)*YuL(0,0)).imag();
	Cquqd8_1131r_LNP = (Cquqd8_00_LNP*YdL(0)*YuL(2,0)).real();
	Cquqd8_1131i_LNP = (Cquqd8_00_LNP*YdL(0)*YuL(2,0)).imag();
	Cquqd8_1133r_LNP = (Cquqd8_00_LNP*YdL(2)*YuL(0,0)).real();
	Cquqd8_1133i_LNP = (Cquqd8_00_LNP*YdL(2)*YuL(0,0)).imag();
	Cquqd8_1211r_LNP = (2*Cquqd8_00_LNP*YdL(0)*YuL(0,1)).real();
	Cquqd8_1211i_LNP = (2*Cquqd8_00_LNP*YdL(0)*YuL(0,1)).imag();
	Cquqd8_1221r_LNP = (Cquqd8_00_LNP*YdL(0)*YuL(1,1)).real();
	Cquqd8_1221i_LNP = (Cquqd8_00_LNP*YdL(0)*YuL(1,1)).imag();
	Cquqd8_1222r_LNP = (Cquqd8_00_LNP*YdL(1)*YuL(0,1)).real();
	Cquqd8_1222i_LNP = (Cquqd8_00_LNP*YdL(1)*YuL(0,1)).imag();
	Cquqd8_1231r_LNP = (Cquqd8_00_LNP*YdL(0)*YuL(2,1)).real();
	Cquqd8_1231i_LNP = (Cquqd8_00_LNP*YdL(0)*YuL(2,1)).imag();
	Cquqd8_1233r_LNP = (Cquqd8_00_LNP*YdL(2)*YuL(0,1)).real();
	Cquqd8_1233i_LNP = (Cquqd8_00_LNP*YdL(2)*YuL(0,1)).imag();
	Cquqd8_1311r_LNP = (2*Cquqd8_00_LNP*YdL(0)*YuL(0,2)).real();
	Cquqd8_1311i_LNP = (2*Cquqd8_00_LNP*YdL(0)*YuL(0,2)).imag();
	Cquqd8_1321r_LNP = (Cquqd8_00_LNP*YdL(0)*YuL(1,2)).real();
	Cquqd8_1321i_LNP = (Cquqd8_00_LNP*YdL(0)*YuL(1,2)).imag();
	Cquqd8_1322r_LNP = (Cquqd8_00_LNP*YdL(1)*YuL(0,2)).real();
	Cquqd8_1322i_LNP = (Cquqd8_00_LNP*YdL(1)*YuL(0,2)).imag();
	Cquqd8_1331r_LNP = (Cquqd8_00_LNP*YdL(0)*YuL(2,2)).real();
	Cquqd8_1331i_LNP = (Cquqd8_00_LNP*YdL(0)*YuL(2,2)).imag();
	Cquqd8_1333r_LNP = (Cquqd8_00_LNP*YdL(2)*YuL(0,2)).real();
	Cquqd8_1333i_LNP = (Cquqd8_00_LNP*YdL(2)*YuL(0,2)).imag();
	Cquqd8_2111r_LNP = (Cquqd8_00_LNP*YdL(0)*YuL(1,0)).real();
	Cquqd8_2111i_LNP = (Cquqd8_00_LNP*YdL(0)*YuL(1,0)).imag();
	Cquqd8_2112r_LNP = (Cquqd8_00_LNP*YdL(1)*YuL(0,0)).real();
	Cquqd8_2112i_LNP = (Cquqd8_00_LNP*YdL(1)*YuL(0,0)).imag();
	Cquqd8_2122r_LNP = (2*Cquqd8_00_LNP*YdL(1)*YuL(1,0)).real();
	Cquqd8_2122i_LNP = (2*Cquqd8_00_LNP*YdL(1)*YuL(1,0)).imag();
	Cquqd8_2132r_LNP = (Cquqd8_00_LNP*YdL(1)*YuL(2,0)).real();
	Cquqd8_2132i_LNP = (Cquqd8_00_LNP*YdL(1)*YuL(2,0)).imag();
	Cquqd8_2133r_LNP = (Cquqd8_00_LNP*YdL(2)*YuL(1,0)).real();
	Cquqd8_2133i_LNP = (Cquqd8_00_LNP*YdL(2)*YuL(1,0)).imag();
	Cquqd8_2211r_LNP = (Cquqd8_00_LNP*YdL(0)*YuL(1,1)).real();
	Cquqd8_2211i_LNP = (Cquqd8_00_LNP*YdL(0)*YuL(1,1)).imag();
	Cquqd8_2212r_LNP = (Cquqd8_00_LNP*YdL(1)*YuL(0,1)).real();
	Cquqd8_2212i_LNP = (Cquqd8_00_LNP*YdL(1)*YuL(0,1)).imag();
	Cquqd8_2222r_LNP = (2*Cquqd8_00_LNP*YdL(1)*YuL(1,1)).real();
	Cquqd8_2222i_LNP = (2*Cquqd8_00_LNP*YdL(1)*YuL(1,1)).imag();
	Cquqd8_2232r_LNP = (Cquqd8_00_LNP*YdL(1)*YuL(2,1)).real();
	Cquqd8_2232i_LNP = (Cquqd8_00_LNP*YdL(1)*YuL(2,1)).imag();
	Cquqd8_2233r_LNP = (Cquqd8_00_LNP*YdL(2)*YuL(1,1)).real();
	Cquqd8_2233i_LNP = (Cquqd8_00_LNP*YdL(2)*YuL(1,1)).imag();
	Cquqd8_2311r_LNP = (Cquqd8_00_LNP*YdL(0)*YuL(1,2)).real();
	Cquqd8_2311i_LNP = (Cquqd8_00_LNP*YdL(0)*YuL(1,2)).imag();
	Cquqd8_2312r_LNP = (Cquqd8_00_LNP*YdL(1)*YuL(0,2)).real();
	Cquqd8_2312i_LNP = (Cquqd8_00_LNP*YdL(1)*YuL(0,2)).imag();
	Cquqd8_2322r_LNP = (2*Cquqd8_00_LNP*YdL(1)*YuL(1,2)).real();
	Cquqd8_2322i_LNP = (2*Cquqd8_00_LNP*YdL(1)*YuL(1,2)).imag();
	Cquqd8_2332r_LNP = (Cquqd8_00_LNP*YdL(1)*YuL(2,2)).real();
	Cquqd8_2332i_LNP = (Cquqd8_00_LNP*YdL(1)*YuL(2,2)).imag();
	Cquqd8_2333r_LNP = (Cquqd8_00_LNP*YdL(2)*YuL(1,2)).real();
	Cquqd8_2333i_LNP = (Cquqd8_00_LNP*YdL(2)*YuL(1,2)).imag();
	Cquqd8_3111r_LNP = (Cquqd8_00_LNP*YdL(0)*YuL(2,0)).real();
	Cquqd8_3111i_LNP = (Cquqd8_00_LNP*YdL(0)*YuL(2,0)).imag();
	Cquqd8_3113r_LNP = (Cquqd8_00_LNP*YdL(2)*YuL(0,0)).real();
	Cquqd8_3113i_LNP = (Cquqd8_00_LNP*YdL(2)*YuL(0,0)).imag();
	Cquqd8_3122r_LNP = (Cquqd8_00_LNP*YdL(1)*YuL(2,0)).real();
	Cquqd8_3122i_LNP = (Cquqd8_00_LNP*YdL(1)*YuL(2,0)).imag();
	Cquqd8_3123r_LNP = (Cquqd8_00_LNP*YdL(2)*YuL(1,0)).real();
	Cquqd8_3123i_LNP = (Cquqd8_00_LNP*YdL(2)*YuL(1,0)).imag();
	Cquqd8_3133r_LNP = (2*Cquqd8_00_LNP*YdL(2)*YuL(2,0)).real();
	Cquqd8_3133i_LNP = (2*Cquqd8_00_LNP*YdL(2)*YuL(2,0)).imag();
	Cquqd8_3211r_LNP = (Cquqd8_00_LNP*YdL(0)*YuL(2,1)).real();
	Cquqd8_3211i_LNP = (Cquqd8_00_LNP*YdL(0)*YuL(2,1)).imag();
	Cquqd8_3213r_LNP = (Cquqd8_00_LNP*YdL(2)*YuL(0,1)).real();
	Cquqd8_3213i_LNP = (Cquqd8_00_LNP*YdL(2)*YuL(0,1)).imag();
	Cquqd8_3222r_LNP = (Cquqd8_00_LNP*YdL(1)*YuL(2,1)).real();
	Cquqd8_3222i_LNP = (Cquqd8_00_LNP*YdL(1)*YuL(2,1)).imag();
	Cquqd8_3223r_LNP = (Cquqd8_00_LNP*YdL(2)*YuL(1,1)).real();
	Cquqd8_3223i_LNP = (Cquqd8_00_LNP*YdL(2)*YuL(1,1)).imag();
	Cquqd8_3233r_LNP = (2*Cquqd8_00_LNP*YdL(2)*YuL(2,1)).real();
	Cquqd8_3233i_LNP = (2*Cquqd8_00_LNP*YdL(2)*YuL(2,1)).imag();
	Cquqd8_3311r_LNP = (Cquqd8_00_LNP*YdL(0)*YuL(2,2)).real();
	Cquqd8_3311i_LNP = (Cquqd8_00_LNP*YdL(0)*YuL(2,2)).imag();
	Cquqd8_3313r_LNP = (Cquqd8_00_LNP*YdL(2)*YuL(0,2)).real();
	Cquqd8_3313i_LNP = (Cquqd8_00_LNP*YdL(2)*YuL(0,2)).imag();
	Cquqd8_3322r_LNP = (Cquqd8_00_LNP*YdL(1)*YuL(2,2)).real();
	Cquqd8_3322i_LNP = (Cquqd8_00_LNP*YdL(1)*YuL(2,2)).imag();
	Cquqd8_3323r_LNP = (Cquqd8_00_LNP*YdL(2)*YuL(1,2)).real();
	Cquqd8_3323i_LNP = (Cquqd8_00_LNP*YdL(2)*YuL(1,2)).imag();
	Cquqd8_3333r_LNP = (2*Cquqd8_00_LNP*YdL(2)*YuL(2,2)).real();
	Cquqd8_3333i_LNP = (2*Cquqd8_00_LNP*YdL(2)*YuL(2,2)).imag();

}

bool NPSMEFTd6MFV::PostUpdate() {

	GenerateSMInitialConditions();
	
    setNPSMEFTd6GeneralParameters();
    
    if (!NPSMEFTd6General::PostUpdate()) return (false);
    
    return (true);    
}