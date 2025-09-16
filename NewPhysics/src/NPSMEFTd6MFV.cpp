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
	"CG","CW","CHG","CHW","CHB","CHWB","CHD","CHbox","CH","CHl1","CHl3","CHe","Cll","Cll1","Cee","Cle","CuH_0","CuH_u","CuH_d","CuG_0","CuG_u","CuG_d","CuW_0","CuW_u","CuW_d","CuB_0","CuB_u","CuB_d","CdH_0","CdH_u","CdH_d","CdG_0","CdG_u","CdG_d","CdW_0","CdW_u","CdW_d","CdB_0","CdB_u","CdB_d","CHq1_0","CHq1_u","CHq1_d","CHq3_0","CHq3_u","CHq3_d","CHu_0","CHu_u","CHd_0","CHd_d","CHud_ud","Clq1_0","Clq1_u","Clq1_d","Clq3_0","Clq3_u","Clq3_d","Cqe_0","Cqe_u","Cqe_d","Clu_0","Clu_u","Ceu_0","Ceu_u","Cld_0","Cld_d","Ced_0","Ced_d","Cqq1_00","Cqq1_0u","Cqq1_0d","Cqq1_u0","Cqq1_uu","Cqq1_ud","Cqq1_d0","Cqq1_du","Cqq1_dd","Cqq3_00","Cqq3_0u","Cqq3_0d","Cqq3_u0","Cqq3_uu","Cqq3_ud","Cqq3_d0","Cqq3_du","Cqq3_dd","Cuu_00","Cuu_0u","Cuu_u0","Cuu_uu","Cdd_00","Cdd_0d","Cdd_d0","Cdd_dd","Cud1_00","Cud1_u0","Cud1_0d","Cud1_ud","Cud8_00","Cud8_u0","Cud8_0d","Cud8_ud","Cqu1_00","Cqu1_u0","Cqu1_d0","Cqu1_0u","Cqu1_uu","Cqu1_du","Cqu8_00","Cqu8_u0","Cqu8_d0","Cqu8_0u","Cqu8_uu","Cqu8_du","Cqd1_00","Cqd1_u0","Cqd1_d0","Cqd1_0d","Cqd1_ud","Cqd1_dd","Cqd8_00","Cqd8_u0","Cqd8_d0","Cqd8_0d","Cqd8_ud","Cqd8_dd","Cquqd1","Cquqd8",
    "Lambda_NP"
};

NPSMEFTd6MFV::NPSMEFTd6MFV(): NPSMEFTd6General() {
    setModelName("NPSMEFTd6MFV");
	ModelParamMap.insert(std::make_pair("CG",std::cref(CG)));
	ModelParamMap.insert(std::make_pair("CW",std::cref(CW)));
	ModelParamMap.insert(std::make_pair("CHG",std::cref(CHG)));
	ModelParamMap.insert(std::make_pair("CHW",std::cref(CHW)));
	ModelParamMap.insert(std::make_pair("CHB",std::cref(CHB)));
	ModelParamMap.insert(std::make_pair("CHWB",std::cref(CHWB)));
	ModelParamMap.insert(std::make_pair("CHD",std::cref(CHD)));
	ModelParamMap.insert(std::make_pair("CHbox",std::cref(CHbox)));
	ModelParamMap.insert(std::make_pair("CH",std::cref(CH)));
	ModelParamMap.insert(std::make_pair("CHl1",std::cref(CHl1)));
	ModelParamMap.insert(std::make_pair("CHl3",std::cref(CHl3)));
	ModelParamMap.insert(std::make_pair("CHe",std::cref(CHe)));
	ModelParamMap.insert(std::make_pair("Cll",std::cref(Cll)));
	ModelParamMap.insert(std::make_pair("Cll1",std::cref(Cll1)));
	ModelParamMap.insert(std::make_pair("Cee",std::cref(Cee)));
	ModelParamMap.insert(std::make_pair("Cle",std::cref(Cle)));
	ModelParamMap.insert(std::make_pair("CuH_0",std::cref(CuH_0)));
	ModelParamMap.insert(std::make_pair("CuH_u",std::cref(CuH_u)));
	ModelParamMap.insert(std::make_pair("CuH_d",std::cref(CuH_d)));
	ModelParamMap.insert(std::make_pair("CuG_0",std::cref(CuG_0)));
	ModelParamMap.insert(std::make_pair("CuG_u",std::cref(CuG_u)));
	ModelParamMap.insert(std::make_pair("CuG_d",std::cref(CuG_d)));
	ModelParamMap.insert(std::make_pair("CuW_0",std::cref(CuW_0)));
	ModelParamMap.insert(std::make_pair("CuW_u",std::cref(CuW_u)));
	ModelParamMap.insert(std::make_pair("CuW_d",std::cref(CuW_d)));
	ModelParamMap.insert(std::make_pair("CuB_0",std::cref(CuB_0)));
	ModelParamMap.insert(std::make_pair("CuB_u",std::cref(CuB_u)));
	ModelParamMap.insert(std::make_pair("CuB_d",std::cref(CuB_d)));
	ModelParamMap.insert(std::make_pair("CdH_0",std::cref(CdH_0)));
	ModelParamMap.insert(std::make_pair("CdH_u",std::cref(CdH_u)));
	ModelParamMap.insert(std::make_pair("CdH_d",std::cref(CdH_d)));
	ModelParamMap.insert(std::make_pair("CdG_0",std::cref(CdG_0)));
	ModelParamMap.insert(std::make_pair("CdG_u",std::cref(CdG_u)));
	ModelParamMap.insert(std::make_pair("CdG_d",std::cref(CdG_d)));
	ModelParamMap.insert(std::make_pair("CdW_0",std::cref(CdW_0)));
	ModelParamMap.insert(std::make_pair("CdW_u",std::cref(CdW_u)));
	ModelParamMap.insert(std::make_pair("CdW_d",std::cref(CdW_d)));
	ModelParamMap.insert(std::make_pair("CdB_0",std::cref(CdB_0)));
	ModelParamMap.insert(std::make_pair("CdB_u",std::cref(CdB_u)));
	ModelParamMap.insert(std::make_pair("CdB_d",std::cref(CdB_d)));
	ModelParamMap.insert(std::make_pair("CHq1_0",std::cref(CHq1_0)));
	ModelParamMap.insert(std::make_pair("CHq1_u",std::cref(CHq1_u)));
	ModelParamMap.insert(std::make_pair("CHq1_d",std::cref(CHq1_d)));
	ModelParamMap.insert(std::make_pair("CHq3_0",std::cref(CHq3_0)));
	ModelParamMap.insert(std::make_pair("CHq3_u",std::cref(CHq3_u)));
	ModelParamMap.insert(std::make_pair("CHq3_d",std::cref(CHq3_d)));
	ModelParamMap.insert(std::make_pair("CHu_0",std::cref(CHu_0)));
	ModelParamMap.insert(std::make_pair("CHu_u",std::cref(CHu_u)));
	ModelParamMap.insert(std::make_pair("CHd_0",std::cref(CHd_0)));
	ModelParamMap.insert(std::make_pair("CHd_d",std::cref(CHd_d)));
	ModelParamMap.insert(std::make_pair("CHud_ud",std::cref(CHud_ud)));
	ModelParamMap.insert(std::make_pair("Clq1_0",std::cref(Clq1_0)));
	ModelParamMap.insert(std::make_pair("Clq1_u",std::cref(Clq1_u)));
	ModelParamMap.insert(std::make_pair("Clq1_d",std::cref(Clq1_d)));
	ModelParamMap.insert(std::make_pair("Clq3_0",std::cref(Clq3_0)));
	ModelParamMap.insert(std::make_pair("Clq3_u",std::cref(Clq3_u)));
	ModelParamMap.insert(std::make_pair("Clq3_d",std::cref(Clq3_d)));
	ModelParamMap.insert(std::make_pair("Cqe_0",std::cref(Cqe_0)));
	ModelParamMap.insert(std::make_pair("Cqe_u",std::cref(Cqe_u)));
	ModelParamMap.insert(std::make_pair("Cqe_d",std::cref(Cqe_d)));
	ModelParamMap.insert(std::make_pair("Clu_0",std::cref(Clu_0)));
	ModelParamMap.insert(std::make_pair("Clu_u",std::cref(Clu_u)));
	ModelParamMap.insert(std::make_pair("Ceu_0",std::cref(Ceu_0)));
	ModelParamMap.insert(std::make_pair("Ceu_u",std::cref(Ceu_u)));
	ModelParamMap.insert(std::make_pair("Cld_0",std::cref(Cld_0)));
	ModelParamMap.insert(std::make_pair("Cld_d",std::cref(Cld_d)));
	ModelParamMap.insert(std::make_pair("Ced_0",std::cref(Ced_0)));
	ModelParamMap.insert(std::make_pair("Ced_d",std::cref(Ced_d)));
	ModelParamMap.insert(std::make_pair("Cqq1_00",std::cref(Cqq1_00)));
	ModelParamMap.insert(std::make_pair("Cqq1_0u",std::cref(Cqq1_0u)));
	ModelParamMap.insert(std::make_pair("Cqq1_0d",std::cref(Cqq1_0d)));
	ModelParamMap.insert(std::make_pair("Cqq1_u0",std::cref(Cqq1_u0)));
	ModelParamMap.insert(std::make_pair("Cqq1_uu",std::cref(Cqq1_uu)));
	ModelParamMap.insert(std::make_pair("Cqq1_ud",std::cref(Cqq1_ud)));
	ModelParamMap.insert(std::make_pair("Cqq1_d0",std::cref(Cqq1_d0)));
	ModelParamMap.insert(std::make_pair("Cqq1_du",std::cref(Cqq1_du)));
	ModelParamMap.insert(std::make_pair("Cqq1_dd",std::cref(Cqq1_dd)));
	ModelParamMap.insert(std::make_pair("Cqq3_00",std::cref(Cqq3_00)));
	ModelParamMap.insert(std::make_pair("Cqq3_0u",std::cref(Cqq3_0u)));
	ModelParamMap.insert(std::make_pair("Cqq3_0d",std::cref(Cqq3_0d)));
	ModelParamMap.insert(std::make_pair("Cqq3_u0",std::cref(Cqq3_u0)));
	ModelParamMap.insert(std::make_pair("Cqq3_uu",std::cref(Cqq3_uu)));
	ModelParamMap.insert(std::make_pair("Cqq3_ud",std::cref(Cqq3_ud)));
	ModelParamMap.insert(std::make_pair("Cqq3_d0",std::cref(Cqq3_d0)));
	ModelParamMap.insert(std::make_pair("Cqq3_du",std::cref(Cqq3_du)));
	ModelParamMap.insert(std::make_pair("Cqq3_dd",std::cref(Cqq3_dd)));
	ModelParamMap.insert(std::make_pair("Cuu_00",std::cref(Cuu_00)));
	ModelParamMap.insert(std::make_pair("Cuu_0u",std::cref(Cuu_0u)));
	ModelParamMap.insert(std::make_pair("Cuu_u0",std::cref(Cuu_u0)));
	ModelParamMap.insert(std::make_pair("Cuu_uu",std::cref(Cuu_uu)));
	ModelParamMap.insert(std::make_pair("Cdd_00",std::cref(Cdd_00)));
	ModelParamMap.insert(std::make_pair("Cdd_0d",std::cref(Cdd_0d)));
	ModelParamMap.insert(std::make_pair("Cdd_d0",std::cref(Cdd_d0)));
	ModelParamMap.insert(std::make_pair("Cdd_dd",std::cref(Cdd_dd)));
	ModelParamMap.insert(std::make_pair("Cud1_00",std::cref(Cud1_00)));
	ModelParamMap.insert(std::make_pair("Cud1_u0",std::cref(Cud1_u0)));
	ModelParamMap.insert(std::make_pair("Cud1_0d",std::cref(Cud1_0d)));
	ModelParamMap.insert(std::make_pair("Cud1_ud",std::cref(Cud1_ud)));
	ModelParamMap.insert(std::make_pair("Cud8_00",std::cref(Cud8_00)));
	ModelParamMap.insert(std::make_pair("Cud8_u0",std::cref(Cud8_u0)));
	ModelParamMap.insert(std::make_pair("Cud8_0d",std::cref(Cud8_0d)));
	ModelParamMap.insert(std::make_pair("Cud8_ud",std::cref(Cud8_ud)));
	ModelParamMap.insert(std::make_pair("Cqu1_00",std::cref(Cqu1_00)));
	ModelParamMap.insert(std::make_pair("Cqu1_u0",std::cref(Cqu1_u0)));
	ModelParamMap.insert(std::make_pair("Cqu1_d0",std::cref(Cqu1_d0)));
	ModelParamMap.insert(std::make_pair("Cqu1_0u",std::cref(Cqu1_0u)));
	ModelParamMap.insert(std::make_pair("Cqu1_uu",std::cref(Cqu1_uu)));
	ModelParamMap.insert(std::make_pair("Cqu1_du",std::cref(Cqu1_du)));
	ModelParamMap.insert(std::make_pair("Cqu8_00",std::cref(Cqu8_00)));
	ModelParamMap.insert(std::make_pair("Cqu8_u0",std::cref(Cqu8_u0)));
	ModelParamMap.insert(std::make_pair("Cqu8_d0",std::cref(Cqu8_d0)));
	ModelParamMap.insert(std::make_pair("Cqu8_0u",std::cref(Cqu8_0u)));
	ModelParamMap.insert(std::make_pair("Cqu8_uu",std::cref(Cqu8_uu)));
	ModelParamMap.insert(std::make_pair("Cqu8_du",std::cref(Cqu8_du)));
	ModelParamMap.insert(std::make_pair("Cqd1_00",std::cref(Cqd1_00)));
	ModelParamMap.insert(std::make_pair("Cqd1_u0",std::cref(Cqd1_u0)));
	ModelParamMap.insert(std::make_pair("Cqd1_d0",std::cref(Cqd1_d0)));
	ModelParamMap.insert(std::make_pair("Cqd1_0d",std::cref(Cqd1_0d)));
	ModelParamMap.insert(std::make_pair("Cqd1_ud",std::cref(Cqd1_ud)));
	ModelParamMap.insert(std::make_pair("Cqd1_dd",std::cref(Cqd1_dd)));
	ModelParamMap.insert(std::make_pair("Cqd8_00",std::cref(Cqd8_00)));
	ModelParamMap.insert(std::make_pair("Cqd8_u0",std::cref(Cqd8_u0)));
	ModelParamMap.insert(std::make_pair("Cqd8_d0",std::cref(Cqd8_d0)));
	ModelParamMap.insert(std::make_pair("Cqd8_0d",std::cref(Cqd8_0d)));
	ModelParamMap.insert(std::make_pair("Cqd8_ud",std::cref(Cqd8_ud)));
	ModelParamMap.insert(std::make_pair("Cqd8_dd",std::cref(Cqd8_dd)));
	ModelParamMap.insert(std::make_pair("Cquqd1",std::cref(Cquqd1)));
	ModelParamMap.insert(std::make_pair("Cquqd8",std::cref(Cquqd8)));

}

bool NPSMEFTd6MFV::Init(const std::map<std::string, double>& DPars) {
	if (SMEFTBasisFlag == "UP") return false;
    return (NPSMEFTd6General::Init(DPars));
}

void NPSMEFTd6MFV::setParameter(const std::string name, const double& value) {
	if (name.compare("CG") == 0) {
		CG = value;
	} else if (name.compare("CW") == 0) {
		CW = value;
	} else if (name.compare("CHG") == 0) {
		CHG = value;
	} else if (name.compare("CHW") == 0) {
		CHW = value;
	} else if (name.compare("CHB") == 0) {
		CHB = value;
	} else if (name.compare("CHWB") == 0) {
		CHWB = value;
	} else if (name.compare("CHD") == 0) {
		CHD = value;
	} else if (name.compare("CHbox") == 0) {
		CHbox = value;
	} else if (name.compare("CH") == 0) {
		CH = value;
	} else if (name.compare("CHl1") == 0) {
		CHl1 = value;
	} else if (name.compare("CHl3") == 0) {
		CHl3 = value;
	} else if (name.compare("CHe") == 0) {
		CHe = value;
	} else if (name.compare("Cll") == 0) {
		Cll = value;
	} else if (name.compare("Cll1") == 0) {
		Cll1 = value;
	} else if (name.compare("Cee") == 0) {
		Cee = value;
	} else if (name.compare("Cle") == 0) {
		Cle = value;
	} else if (name.compare("CuH_0") == 0) {
		CuH_0 = value;
	} else if (name.compare("CuH_u") == 0) {
		CuH_u = value;
	} else if (name.compare("CuH_d") == 0) {
		CuH_d = value;
	} else if (name.compare("CuG_0") == 0) {
		CuG_0 = value;
	} else if (name.compare("CuG_u") == 0) {
		CuG_u = value;
	} else if (name.compare("CuG_d") == 0) {
		CuG_d = value;
	} else if (name.compare("CuW_0") == 0) {
		CuW_0 = value;
	} else if (name.compare("CuW_u") == 0) {
		CuW_u = value;
	} else if (name.compare("CuW_d") == 0) {
		CuW_d = value;
	} else if (name.compare("CuB_0") == 0) {
		CuB_0 = value;
	} else if (name.compare("CuB_u") == 0) {
		CuB_u = value;
	} else if (name.compare("CuB_d") == 0) {
		CuB_d = value;
	} else if (name.compare("CdH_0") == 0) {
		CdH_0 = value;
	} else if (name.compare("CdH_u") == 0) {
		CdH_u = value;
	} else if (name.compare("CdH_d") == 0) {
		CdH_d = value;
	} else if (name.compare("CdG_0") == 0) {
		CdG_0 = value;
	} else if (name.compare("CdG_u") == 0) {
		CdG_u = value;
	} else if (name.compare("CdG_d") == 0) {
		CdG_d = value;
	} else if (name.compare("CdW_0") == 0) {
		CdW_0 = value;
	} else if (name.compare("CdW_u") == 0) {
		CdW_u = value;
	} else if (name.compare("CdW_d") == 0) {
		CdW_d = value;
	} else if (name.compare("CdB_0") == 0) {
		CdB_0 = value;
	} else if (name.compare("CdB_u") == 0) {
		CdB_u = value;
	} else if (name.compare("CdB_d") == 0) {
		CdB_d = value;
	} else if (name.compare("CHq1_0") == 0) {
		CHq1_0 = value;
	} else if (name.compare("CHq1_u") == 0) {
		CHq1_u = value;
	} else if (name.compare("CHq1_d") == 0) {
		CHq1_d = value;
	} else if (name.compare("CHq3_0") == 0) {
		CHq3_0 = value;
	} else if (name.compare("CHq3_u") == 0) {
		CHq3_u = value;
	} else if (name.compare("CHq3_d") == 0) {
		CHq3_d = value;
	} else if (name.compare("CHu_0") == 0) {
		CHu_0 = value;
	} else if (name.compare("CHu_u") == 0) {
		CHu_u = value;
	} else if (name.compare("CHd_0") == 0) {
		CHd_0 = value;
	} else if (name.compare("CHd_d") == 0) {
		CHd_d = value;
	} else if (name.compare("CHud_ud") == 0) {
		CHud_ud = value;
	} else if (name.compare("Clq1_0") == 0) {
		Clq1_0 = value;
	} else if (name.compare("Clq1_u") == 0) {
		Clq1_u = value;
	} else if (name.compare("Clq1_d") == 0) {
		Clq1_d = value;
	} else if (name.compare("Clq3_0") == 0) {
		Clq3_0 = value;
	} else if (name.compare("Clq3_u") == 0) {
		Clq3_u = value;
	} else if (name.compare("Clq3_d") == 0) {
		Clq3_d = value;
	} else if (name.compare("Cqe_0") == 0) {
		Cqe_0 = value;
	} else if (name.compare("Cqe_u") == 0) {
		Cqe_u = value;
	} else if (name.compare("Cqe_d") == 0) {
		Cqe_d = value;
	} else if (name.compare("Clu_0") == 0) {
		Clu_0 = value;
	} else if (name.compare("Clu_u") == 0) {
		Clu_u = value;
	} else if (name.compare("Ceu_0") == 0) {
		Ceu_0 = value;
	} else if (name.compare("Ceu_u") == 0) {
		Ceu_u = value;
	} else if (name.compare("Cld_0") == 0) {
		Cld_0 = value;
	} else if (name.compare("Cld_d") == 0) {
		Cld_d = value;
	} else if (name.compare("Ced_0") == 0) {
		Ced_0 = value;
	} else if (name.compare("Ced_d") == 0) {
		Ced_d = value;
	} else if (name.compare("Cqq1_00") == 0) {
		Cqq1_00 = value;
	} else if (name.compare("Cqq1_0u") == 0) {
		Cqq1_0u = value;
	} else if (name.compare("Cqq1_0d") == 0) {
		Cqq1_0d = value;
	} else if (name.compare("Cqq1_u0") == 0) {
		Cqq1_u0 = value;
	} else if (name.compare("Cqq1_uu") == 0) {
		Cqq1_uu = value;
	} else if (name.compare("Cqq1_ud") == 0) {
		Cqq1_ud = value;
	} else if (name.compare("Cqq1_d0") == 0) {
		Cqq1_d0 = value;
	} else if (name.compare("Cqq1_du") == 0) {
		Cqq1_du = value;
	} else if (name.compare("Cqq1_dd") == 0) {
		Cqq1_dd = value;
	} else if (name.compare("Cqq3_00") == 0) {
		Cqq3_00 = value;
	} else if (name.compare("Cqq3_0u") == 0) {
		Cqq3_0u = value;
	} else if (name.compare("Cqq3_0d") == 0) {
		Cqq3_0d = value;
	} else if (name.compare("Cqq3_u0") == 0) {
		Cqq3_u0 = value;
	} else if (name.compare("Cqq3_uu") == 0) {
		Cqq3_uu = value;
	} else if (name.compare("Cqq3_ud") == 0) {
		Cqq3_ud = value;
	} else if (name.compare("Cqq3_d0") == 0) {
		Cqq3_d0 = value;
	} else if (name.compare("Cqq3_du") == 0) {
		Cqq3_du = value;
	} else if (name.compare("Cqq3_dd") == 0) {
		Cqq3_dd = value;
	} else if (name.compare("Cuu_00") == 0) {
		Cuu_00 = value;
	} else if (name.compare("Cuu_0u") == 0) {
		Cuu_0u = value;
	} else if (name.compare("Cuu_u0") == 0) {
		Cuu_u0 = value;
	} else if (name.compare("Cuu_uu") == 0) {
		Cuu_uu = value;
	} else if (name.compare("Cdd_00") == 0) {
		Cdd_00 = value;
	} else if (name.compare("Cdd_0d") == 0) {
		Cdd_0d = value;
	} else if (name.compare("Cdd_d0") == 0) {
		Cdd_d0 = value;
	} else if (name.compare("Cdd_dd") == 0) {
		Cdd_dd = value;
	} else if (name.compare("Cud1_00") == 0) {
		Cud1_00 = value;
	} else if (name.compare("Cud1_u0") == 0) {
		Cud1_u0 = value;
	} else if (name.compare("Cud1_0d") == 0) {
		Cud1_0d = value;
	} else if (name.compare("Cud1_ud") == 0) {
		Cud1_ud = value;
	} else if (name.compare("Cud8_00") == 0) {
		Cud8_00 = value;
	} else if (name.compare("Cud8_u0") == 0) {
		Cud8_u0 = value;
	} else if (name.compare("Cud8_0d") == 0) {
		Cud8_0d = value;
	} else if (name.compare("Cud8_ud") == 0) {
		Cud8_ud = value;
	} else if (name.compare("Cqu1_00") == 0) {
		Cqu1_00 = value;
	} else if (name.compare("Cqu1_u0") == 0) {
		Cqu1_u0 = value;
	} else if (name.compare("Cqu1_d0") == 0) {
		Cqu1_d0 = value;
	} else if (name.compare("Cqu1_0u") == 0) {
		Cqu1_0u = value;
	} else if (name.compare("Cqu1_uu") == 0) {
		Cqu1_uu = value;
	} else if (name.compare("Cqu1_du") == 0) {
		Cqu1_du = value;
	} else if (name.compare("Cqu8_00") == 0) {
		Cqu8_00 = value;
	} else if (name.compare("Cqu8_u0") == 0) {
		Cqu8_u0 = value;
	} else if (name.compare("Cqu8_d0") == 0) {
		Cqu8_d0 = value;
	} else if (name.compare("Cqu8_0u") == 0) {
		Cqu8_0u = value;
	} else if (name.compare("Cqu8_uu") == 0) {
		Cqu8_uu = value;
	} else if (name.compare("Cqu8_du") == 0) {
		Cqu8_du = value;
	} else if (name.compare("Cqd1_00") == 0) {
		Cqd1_00 = value;
	} else if (name.compare("Cqd1_u0") == 0) {
		Cqd1_u0 = value;
	} else if (name.compare("Cqd1_d0") == 0) {
		Cqd1_d0 = value;
	} else if (name.compare("Cqd1_0d") == 0) {
		Cqd1_0d = value;
	} else if (name.compare("Cqd1_ud") == 0) {
		Cqd1_ud = value;
	} else if (name.compare("Cqd1_dd") == 0) {
		Cqd1_dd = value;
	} else if (name.compare("Cqd8_00") == 0) {
		Cqd8_00 = value;
	} else if (name.compare("Cqd8_u0") == 0) {
		Cqd8_u0 = value;
	} else if (name.compare("Cqd8_d0") == 0) {
		Cqd8_d0 = value;
	} else if (name.compare("Cqd8_0d") == 0) {
		Cqd8_0d = value;
	} else if (name.compare("Cqd8_ud") == 0) {
		Cqd8_ud = value;
	} else if (name.compare("Cqd8_dd") == 0) {
		Cqd8_dd = value;
	} else if (name.compare("Cquqd1") == 0) {
		Cquqd1 = value;
	} else if (name.compare("Cquqd8") == 0) {
		Cquqd8 = value;
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

	CG_LNP = CG;

	CW_LNP = CW;

	CHG_LNP = CHG;

	CHW_LNP = CHW;

	CHB_LNP = CHB;

	CHWB_LNP = CHWB;

	CHD_LNP = CHD;

	CHbox_LNP = CHbox;

	CH_LNP = CH;

	CHl1_11r_LNP = CHl1;
	CHl1_22r_LNP = CHl1;
	CHl1_33r_LNP = CHl1;

	CHl3_11r_LNP = CHl3;
	CHl3_22r_LNP = CHl3;
	CHl3_33r_LNP = CHl3;

	CHe_11r_LNP = CHe;
	CHe_22r_LNP = CHe;
	CHe_33r_LNP = CHe;

	Cll_1111r_LNP = Cll+Cll1;
	Cll_2222r_LNP = Cll+Cll1;
	Cll_3333r_LNP = Cll+Cll1;
	Cll_1122r_LNP = Cll;
	Cll_1133r_LNP = Cll;
	Cll_2233r_LNP = Cll;
	Cll_1221r_LNP = Cll1;
	Cll_1331r_LNP = Cll1;
	Cll_2332r_LNP = Cll1;

	CuH_11r_LNP = (CuH_d*YdduL(0,0) + CuH_u*YucuuL(0,0) + CuH_0*YuL(0,0)).real();
	CuH_11i_LNP = (CuH_d*YdduL(0,0) + CuH_u*YucuuL(0,0) + CuH_0*YuL(0,0)).imag();
	CuH_12r_LNP = (CuH_d*YdduL(0,1) + CuH_u*YucuuL(0,1) + CuH_0*YuL(0,1)).real();
	CuH_12i_LNP = (CuH_d*YdduL(0,1) + CuH_u*YucuuL(0,1) + CuH_0*YuL(0,1)).imag();
	CuH_13r_LNP = (CuH_d*YdduL(0,2) + CuH_u*YucuuL(0,2) + CuH_0*YuL(0,2)).real();
	CuH_13i_LNP = (CuH_d*YdduL(0,2) + CuH_u*YucuuL(0,2) + CuH_0*YuL(0,2)).imag();
	CuH_21r_LNP = (CuH_d*YdduL(1,0) + CuH_u*YucuuL(1,0) + CuH_0*YuL(1,0)).real();
	CuH_21i_LNP = (CuH_d*YdduL(1,0) + CuH_u*YucuuL(1,0) + CuH_0*YuL(1,0)).imag();
	CuH_22r_LNP = (CuH_d*YdduL(1,1) + CuH_u*YucuuL(1,1) + CuH_0*YuL(1,1)).real();
	CuH_22i_LNP = (CuH_d*YdduL(1,1) + CuH_u*YucuuL(1,1) + CuH_0*YuL(1,1)).imag();
	CuH_23r_LNP = (CuH_d*YdduL(1,2) + CuH_u*YucuuL(1,2) + CuH_0*YuL(1,2)).real();
	CuH_23i_LNP = (CuH_d*YdduL(1,2) + CuH_u*YucuuL(1,2) + CuH_0*YuL(1,2)).imag();
	CuH_31r_LNP = (CuH_d*YdduL(2,0) + CuH_u*YucuuL(2,0) + CuH_0*YuL(2,0)).real();
	CuH_31i_LNP = (CuH_d*YdduL(2,0) + CuH_u*YucuuL(2,0) + CuH_0*YuL(2,0)).imag();
	CuH_32r_LNP = (CuH_d*YdduL(2,1) + CuH_u*YucuuL(2,1) + CuH_0*YuL(2,1)).real();
	CuH_32i_LNP = (CuH_d*YdduL(2,1) + CuH_u*YucuuL(2,1) + CuH_0*YuL(2,1)).imag();
	CuH_33r_LNP = (CuH_d*YdduL(2,2) + CuH_u*YucuuL(2,2) + CuH_0*YuL(2,2)).real();
	CuH_33i_LNP = (CuH_d*YdduL(2,2) + CuH_u*YucuuL(2,2) + CuH_0*YuL(2,2)).imag();

	CuG_11r_LNP = (CuG_d*YdduL(0,0) + CuG_u*YucuuL(0,0) + CuG_0*YuL(0,0)).real();
	CuG_11i_LNP = (CuG_d*YdduL(0,0) + CuG_u*YucuuL(0,0) + CuG_0*YuL(0,0)).imag();
	CuG_12r_LNP = (CuG_d*YdduL(0,1) + CuG_u*YucuuL(0,1) + CuG_0*YuL(0,1)).real();
	CuG_12i_LNP = (CuG_d*YdduL(0,1) + CuG_u*YucuuL(0,1) + CuG_0*YuL(0,1)).imag();
	CuG_13r_LNP = (CuG_d*YdduL(0,2) + CuG_u*YucuuL(0,2) + CuG_0*YuL(0,2)).real();
	CuG_13i_LNP = (CuG_d*YdduL(0,2) + CuG_u*YucuuL(0,2) + CuG_0*YuL(0,2)).imag();
	CuG_21r_LNP = (CuG_d*YdduL(1,0) + CuG_u*YucuuL(1,0) + CuG_0*YuL(1,0)).real();
	CuG_21i_LNP = (CuG_d*YdduL(1,0) + CuG_u*YucuuL(1,0) + CuG_0*YuL(1,0)).imag();
	CuG_22r_LNP = (CuG_d*YdduL(1,1) + CuG_u*YucuuL(1,1) + CuG_0*YuL(1,1)).real();
	CuG_22i_LNP = (CuG_d*YdduL(1,1) + CuG_u*YucuuL(1,1) + CuG_0*YuL(1,1)).imag();
	CuG_23r_LNP = (CuG_d*YdduL(1,2) + CuG_u*YucuuL(1,2) + CuG_0*YuL(1,2)).real();
	CuG_23i_LNP = (CuG_d*YdduL(1,2) + CuG_u*YucuuL(1,2) + CuG_0*YuL(1,2)).imag();
	CuG_31r_LNP = (CuG_d*YdduL(2,0) + CuG_u*YucuuL(2,0) + CuG_0*YuL(2,0)).real();
	CuG_31i_LNP = (CuG_d*YdduL(2,0) + CuG_u*YucuuL(2,0) + CuG_0*YuL(2,0)).imag();
	CuG_32r_LNP = (CuG_d*YdduL(2,1) + CuG_u*YucuuL(2,1) + CuG_0*YuL(2,1)).real();
	CuG_32i_LNP = (CuG_d*YdduL(2,1) + CuG_u*YucuuL(2,1) + CuG_0*YuL(2,1)).imag();
	CuG_33r_LNP = (CuG_d*YdduL(2,2) + CuG_u*YucuuL(2,2) + CuG_0*YuL(2,2)).real();
	CuG_33i_LNP = (CuG_d*YdduL(2,2) + CuG_u*YucuuL(2,2) + CuG_0*YuL(2,2)).imag();

	CuW_11r_LNP = (CuW_d*YdduL(0,0) + CuW_u*YucuuL(0,0) + CuW_0*YuL(0,0)).real();
	CuW_11i_LNP = (CuW_d*YdduL(0,0) + CuW_u*YucuuL(0,0) + CuW_0*YuL(0,0)).imag();
	CuW_12r_LNP = (CuW_d*YdduL(0,1) + CuW_u*YucuuL(0,1) + CuW_0*YuL(0,1)).real();
	CuW_12i_LNP = (CuW_d*YdduL(0,1) + CuW_u*YucuuL(0,1) + CuW_0*YuL(0,1)).imag();
	CuW_13r_LNP = (CuW_d*YdduL(0,2) + CuW_u*YucuuL(0,2) + CuW_0*YuL(0,2)).real();
	CuW_13i_LNP = (CuW_d*YdduL(0,2) + CuW_u*YucuuL(0,2) + CuW_0*YuL(0,2)).imag();
	CuW_21r_LNP = (CuW_d*YdduL(1,0) + CuW_u*YucuuL(1,0) + CuW_0*YuL(1,0)).real();
	CuW_21i_LNP = (CuW_d*YdduL(1,0) + CuW_u*YucuuL(1,0) + CuW_0*YuL(1,0)).imag();
	CuW_22r_LNP = (CuW_d*YdduL(1,1) + CuW_u*YucuuL(1,1) + CuW_0*YuL(1,1)).real();
	CuW_22i_LNP = (CuW_d*YdduL(1,1) + CuW_u*YucuuL(1,1) + CuW_0*YuL(1,1)).imag();
	CuW_23r_LNP = (CuW_d*YdduL(1,2) + CuW_u*YucuuL(1,2) + CuW_0*YuL(1,2)).real();
	CuW_23i_LNP = (CuW_d*YdduL(1,2) + CuW_u*YucuuL(1,2) + CuW_0*YuL(1,2)).imag();
	CuW_31r_LNP = (CuW_d*YdduL(2,0) + CuW_u*YucuuL(2,0) + CuW_0*YuL(2,0)).real();
	CuW_31i_LNP = (CuW_d*YdduL(2,0) + CuW_u*YucuuL(2,0) + CuW_0*YuL(2,0)).imag();
	CuW_32r_LNP = (CuW_d*YdduL(2,1) + CuW_u*YucuuL(2,1) + CuW_0*YuL(2,1)).real();
	CuW_32i_LNP = (CuW_d*YdduL(2,1) + CuW_u*YucuuL(2,1) + CuW_0*YuL(2,1)).imag();
	CuW_33r_LNP = (CuW_d*YdduL(2,2) + CuW_u*YucuuL(2,2) + CuW_0*YuL(2,2)).real();
	CuW_33i_LNP = (CuW_d*YdduL(2,2) + CuW_u*YucuuL(2,2) + CuW_0*YuL(2,2)).imag();

	CuB_11r_LNP = (CuB_d*YdduL(0,0) + CuB_u*YucuuL(0,0) + CuB_0*YuL(0,0)).real();
	CuB_11i_LNP = (CuB_d*YdduL(0,0) + CuB_u*YucuuL(0,0) + CuB_0*YuL(0,0)).imag();
	CuB_12r_LNP = (CuB_d*YdduL(0,1) + CuB_u*YucuuL(0,1) + CuB_0*YuL(0,1)).real();
	CuB_12i_LNP = (CuB_d*YdduL(0,1) + CuB_u*YucuuL(0,1) + CuB_0*YuL(0,1)).imag();
	CuB_13r_LNP = (CuB_d*YdduL(0,2) + CuB_u*YucuuL(0,2) + CuB_0*YuL(0,2)).real();
	CuB_13i_LNP = (CuB_d*YdduL(0,2) + CuB_u*YucuuL(0,2) + CuB_0*YuL(0,2)).imag();
	CuB_21r_LNP = (CuB_d*YdduL(1,0) + CuB_u*YucuuL(1,0) + CuB_0*YuL(1,0)).real();
	CuB_21i_LNP = (CuB_d*YdduL(1,0) + CuB_u*YucuuL(1,0) + CuB_0*YuL(1,0)).imag();
	CuB_22r_LNP = (CuB_d*YdduL(1,1) + CuB_u*YucuuL(1,1) + CuB_0*YuL(1,1)).real();
	CuB_22i_LNP = (CuB_d*YdduL(1,1) + CuB_u*YucuuL(1,1) + CuB_0*YuL(1,1)).imag();
	CuB_23r_LNP = (CuB_d*YdduL(1,2) + CuB_u*YucuuL(1,2) + CuB_0*YuL(1,2)).real();
	CuB_23i_LNP = (CuB_d*YdduL(1,2) + CuB_u*YucuuL(1,2) + CuB_0*YuL(1,2)).imag();
	CuB_31r_LNP = (CuB_d*YdduL(2,0) + CuB_u*YucuuL(2,0) + CuB_0*YuL(2,0)).real();
	CuB_31i_LNP = (CuB_d*YdduL(2,0) + CuB_u*YucuuL(2,0) + CuB_0*YuL(2,0)).imag();
	CuB_32r_LNP = (CuB_d*YdduL(2,1) + CuB_u*YucuuL(2,1) + CuB_0*YuL(2,1)).real();
	CuB_32i_LNP = (CuB_d*YdduL(2,1) + CuB_u*YucuuL(2,1) + CuB_0*YuL(2,1)).imag();
	CuB_33r_LNP = (CuB_d*YdduL(2,2) + CuB_u*YucuuL(2,2) + CuB_0*YuL(2,2)).real();
	CuB_33i_LNP = (CuB_d*YdduL(2,2) + CuB_u*YucuuL(2,2) + CuB_0*YuL(2,2)).imag();

	CdH_11r_LNP = (CdH_d*YdddL(0) + CdH_0*YdL(0) + CdH_u*YucudL(0,0)).real();
	CdH_11i_LNP = (CdH_d*YdddL(0) + CdH_0*YdL(0) + CdH_u*YucudL(0,0)).imag();
	CdH_12r_LNP = (CdH_u*YucudL(0,1)).real();
	CdH_12i_LNP = (CdH_u*YucudL(0,1)).imag();
	CdH_13r_LNP = (CdH_u*YucudL(0,2)).real();
	CdH_13i_LNP = (CdH_u*YucudL(0,2)).imag();
	CdH_21r_LNP = (CdH_u*YucudL(1,0)).real();
	CdH_21i_LNP = (CdH_u*YucudL(1,0)).imag();
	CdH_22r_LNP = (CdH_d*YdddL(1) + CdH_0*YdL(1) + CdH_u*YucudL(1,1)).real();
	CdH_22i_LNP = (CdH_d*YdddL(1) + CdH_0*YdL(1) + CdH_u*YucudL(1,1)).imag();
	CdH_23r_LNP = (CdH_u*YucudL(1,2)).real();
	CdH_23i_LNP = (CdH_u*YucudL(1,2)).imag();
	CdH_31r_LNP = (CdH_u*YucudL(2,0)).real();
	CdH_31i_LNP = (CdH_u*YucudL(2,0)).imag();
	CdH_32r_LNP = (CdH_u*YucudL(2,1)).real();
	CdH_32i_LNP = (CdH_u*YucudL(2,1)).imag();
	CdH_33r_LNP = (CdH_d*YdddL(2) + CdH_0*YdL(2) + CdH_u*YucudL(2,2)).real();
	CdH_33i_LNP = (CdH_d*YdddL(2) + CdH_0*YdL(2) + CdH_u*YucudL(2,2)).imag();

	CdG_11r_LNP = (CdG_d*YdddL(0) + CdG_0*YdL(0) + CdG_u*YucudL(0,0)).real();
	CdG_11i_LNP = (CdG_d*YdddL(0) + CdG_0*YdL(0) + CdG_u*YucudL(0,0)).imag();
	CdG_12r_LNP = (CdG_u*YucudL(0,1)).real();
	CdG_12i_LNP = (CdG_u*YucudL(0,1)).imag();
	CdG_13r_LNP = (CdG_u*YucudL(0,2)).real();
	CdG_13i_LNP = (CdG_u*YucudL(0,2)).imag();
	CdG_21r_LNP = (CdG_u*YucudL(1,0)).real();
	CdG_21i_LNP = (CdG_u*YucudL(1,0)).imag();
	CdG_22r_LNP = (CdG_d*YdddL(1) + CdG_0*YdL(1) + CdG_u*YucudL(1,1)).real();
	CdG_22i_LNP = (CdG_d*YdddL(1) + CdG_0*YdL(1) + CdG_u*YucudL(1,1)).imag();
	CdG_23r_LNP = (CdG_u*YucudL(1,2)).real();
	CdG_23i_LNP = (CdG_u*YucudL(1,2)).imag();
	CdG_31r_LNP = (CdG_u*YucudL(2,0)).real();
	CdG_31i_LNP = (CdG_u*YucudL(2,0)).imag();
	CdG_32r_LNP = (CdG_u*YucudL(2,1)).real();
	CdG_32i_LNP = (CdG_u*YucudL(2,1)).imag();
	CdG_33r_LNP = (CdG_d*YdddL(2) + CdG_0*YdL(2) + CdG_u*YucudL(2,2)).real();
	CdG_33i_LNP = (CdG_d*YdddL(2) + CdG_0*YdL(2) + CdG_u*YucudL(2,2)).imag();

	CdW_11r_LNP = (CdW_d*YdddL(0) + CdW_0*YdL(0) + CdW_u*YucudL(0,0)).real();
	CdW_11i_LNP = (CdW_d*YdddL(0) + CdW_0*YdL(0) + CdW_u*YucudL(0,0)).imag();
	CdW_12r_LNP = (CdW_u*YucudL(0,1)).real();
	CdW_12i_LNP = (CdW_u*YucudL(0,1)).imag();
	CdW_13r_LNP = (CdW_u*YucudL(0,2)).real();
	CdW_13i_LNP = (CdW_u*YucudL(0,2)).imag();
	CdW_21r_LNP = (CdW_u*YucudL(1,0)).real();
	CdW_21i_LNP = (CdW_u*YucudL(1,0)).imag();
	CdW_22r_LNP = (CdW_d*YdddL(1) + CdW_0*YdL(1) + CdW_u*YucudL(1,1)).real();
	CdW_22i_LNP = (CdW_d*YdddL(1) + CdW_0*YdL(1) + CdW_u*YucudL(1,1)).imag();
	CdW_23r_LNP = (CdW_u*YucudL(1,2)).real();
	CdW_23i_LNP = (CdW_u*YucudL(1,2)).imag();
	CdW_31r_LNP = (CdW_u*YucudL(2,0)).real();
	CdW_31i_LNP = (CdW_u*YucudL(2,0)).imag();
	CdW_32r_LNP = (CdW_u*YucudL(2,1)).real();
	CdW_32i_LNP = (CdW_u*YucudL(2,1)).imag();
	CdW_33r_LNP = (CdW_d*YdddL(2) + CdW_0*YdL(2) + CdW_u*YucudL(2,2)).real();
	CdW_33i_LNP = (CdW_d*YdddL(2) + CdW_0*YdL(2) + CdW_u*YucudL(2,2)).imag();

	CdB_11r_LNP = (CdB_d*YdddL(0) + CdB_0*YdL(0) + CdB_u*YucudL(0,0)).real();
	CdB_11i_LNP = (CdB_d*YdddL(0) + CdB_0*YdL(0) + CdB_u*YucudL(0,0)).imag();
	CdB_12r_LNP = (CdB_u*YucudL(0,1)).real();
	CdB_12i_LNP = (CdB_u*YucudL(0,1)).imag();
	CdB_13r_LNP = (CdB_u*YucudL(0,2)).real();
	CdB_13i_LNP = (CdB_u*YucudL(0,2)).imag();
	CdB_21r_LNP = (CdB_u*YucudL(1,0)).real();
	CdB_21i_LNP = (CdB_u*YucudL(1,0)).imag();
	CdB_22r_LNP = (CdB_d*YdddL(1) + CdB_0*YdL(1) + CdB_u*YucudL(1,1)).real();
	CdB_22i_LNP = (CdB_d*YdddL(1) + CdB_0*YdL(1) + CdB_u*YucudL(1,1)).imag();
	CdB_23r_LNP = (CdB_u*YucudL(1,2)).real();
	CdB_23i_LNP = (CdB_u*YucudL(1,2)).imag();
	CdB_31r_LNP = (CdB_u*YucudL(2,0)).real();
	CdB_31i_LNP = (CdB_u*YucudL(2,0)).imag();
	CdB_32r_LNP = (CdB_u*YucudL(2,1)).real();
	CdB_32i_LNP = (CdB_u*YucudL(2,1)).imag();
	CdB_33r_LNP = (CdB_d*YdddL(2) + CdB_0*YdL(2) + CdB_u*YucudL(2,2)).real();
	CdB_33i_LNP = (CdB_d*YdddL(2) + CdB_0*YdL(2) + CdB_u*YucudL(2,2)).imag();

	CHq1_11r_LNP = (CHq1_0 + CHq1_d*YddL(0) + CHq1_u*YucuL(0,0)).real();
	CHq1_12r_LNP = (CHq1_u*YucuL(0,1)).real();
	CHq1_12i_LNP = (CHq1_u*YucuL(0,1)).imag();
	CHq1_13r_LNP = (CHq1_u*YucuL(0,2)).real();
	CHq1_13i_LNP = (CHq1_u*YucuL(0,2)).imag();
	CHq1_22r_LNP = (CHq1_0 + CHq1_d*YddL(1) + CHq1_u*YucuL(1,1)).real();
	CHq1_23r_LNP = (CHq1_u*YucuL(1,2)).real();
	CHq1_23i_LNP = (CHq1_u*YucuL(1,2)).imag();
	CHq1_33r_LNP = (CHq1_0 + CHq1_d*YddL(2) + CHq1_u*YucuL(2,2)).real();

	CHq3_11r_LNP = (CHq3_0 + CHq3_d*YddL(0) + CHq3_u*YucuL(0,0)).real();
	CHq3_12r_LNP = (CHq3_u*YucuL(0,1)).real();
	CHq3_12i_LNP = (CHq3_u*YucuL(0,1)).imag();
	CHq3_13r_LNP = (CHq3_u*YucuL(0,2)).real();
	CHq3_13i_LNP = (CHq3_u*YucuL(0,2)).imag();
	CHq3_22r_LNP = (CHq3_0 + CHq3_d*YddL(1) + CHq3_u*YucuL(1,1)).real();
	CHq3_23r_LNP = (CHq3_u*YucuL(1,2)).real();
	CHq3_23i_LNP = (CHq3_u*YucuL(1,2)).imag();
	CHq3_33r_LNP = (CHq3_0 + CHq3_d*YddL(2) + CHq3_u*YucuL(2,2)).real();

	CHu_11r_LNP = (CHu_0 + CHu_u*YuucL(0,0)).real();
	CHu_12r_LNP = (CHu_u*YuucL(0,1)).real();
	CHu_12i_LNP = (CHu_u*YuucL(0,1)).imag();
	CHu_13r_LNP = (CHu_u*YuucL(0,2)).real();
	CHu_13i_LNP = (CHu_u*YuucL(0,2)).imag();
	CHu_22r_LNP = (CHu_0 + CHu_u*YuucL(1,1)).real();
	CHu_23r_LNP = (CHu_u*YuucL(1,2)).real();
	CHu_23i_LNP = (CHu_u*YuucL(1,2)).imag();
	CHu_33r_LNP = (CHu_0 + CHu_u*YuucL(2,2)).real();

	CHd_11r_LNP = (CHd_0 + CHd_d*YddL(0)).real();
	CHd_12r_LNP = 0;
	CHd_12i_LNP = 0;
	CHd_13r_LNP = 0;
	CHd_13i_LNP = 0;
	CHd_22r_LNP = (CHd_0 + CHd_d*YddL(1)).real();
	CHd_23r_LNP = 0;
	CHd_23i_LNP = 0;
	CHd_33r_LNP = (CHd_0 + CHd_d*YddL(2)).real();

	CHud_11r_LNP = (CHud_ud*YudL(0,0)).real();
	CHud_11i_LNP = (CHud_ud*YudL(0,0)).imag();
	CHud_12r_LNP = (CHud_ud*YudL(0,1)).real();
	CHud_12i_LNP = (CHud_ud*YudL(0,1)).imag();
	CHud_13r_LNP = (CHud_ud*YudL(0,2)).real();
	CHud_13i_LNP = (CHud_ud*YudL(0,2)).imag();
	CHud_21r_LNP = (CHud_ud*YudL(1,0)).real();
	CHud_21i_LNP = (CHud_ud*YudL(1,0)).imag();
	CHud_22r_LNP = (CHud_ud*YudL(1,1)).real();
	CHud_22i_LNP = (CHud_ud*YudL(1,1)).imag();
	CHud_23r_LNP = (CHud_ud*YudL(1,2)).real();
	CHud_23i_LNP = (CHud_ud*YudL(1,2)).imag();
	CHud_31r_LNP = (CHud_ud*YudL(2,0)).real();
	CHud_31i_LNP = (CHud_ud*YudL(2,0)).imag();
	CHud_32r_LNP = (CHud_ud*YudL(2,1)).real();
	CHud_32i_LNP = (CHud_ud*YudL(2,1)).imag();
	CHud_33r_LNP = (CHud_ud*YudL(2,2)).real();
	CHud_33i_LNP = (CHud_ud*YudL(2,2)).imag();

	Clq1_1111r_LNP = (Clq1_0 + Clq1_d*YddL(0) + Clq1_u*YucuL(0,0)).real();
	Clq1_1112r_LNP = (Clq1_u*YucuL(0,1)).real();
	Clq1_1112i_LNP = (Clq1_u*YucuL(0,1)).imag();
	Clq1_1113r_LNP = (Clq1_u*YucuL(0,2)).real();
	Clq1_1113i_LNP = (Clq1_u*YucuL(0,2)).imag();
	Clq1_1122r_LNP = (Clq1_0 + Clq1_d*YddL(1) + Clq1_u*YucuL(1,1)).real();
	Clq1_1123r_LNP = (Clq1_u*YucuL(1,2)).real();
	Clq1_1123i_LNP = (Clq1_u*YucuL(1,2)).imag();
	Clq1_1133r_LNP = (Clq1_0 + Clq1_d*YddL(2) + Clq1_u*YucuL(2,2)).real();
	Clq1_2211r_LNP = (Clq1_0 + Clq1_d*YddL(0) + Clq1_u*YucuL(0,0)).real();
	Clq1_2212r_LNP = (Clq1_u*YucuL(0,1)).real();
	Clq1_2212i_LNP = (Clq1_u*YucuL(0,1)).imag();
	Clq1_2213r_LNP = (Clq1_u*YucuL(0,2)).real();
	Clq1_2213i_LNP = (Clq1_u*YucuL(0,2)).imag();
	Clq1_2222r_LNP = (Clq1_0 + Clq1_d*YddL(1) + Clq1_u*YucuL(1,1)).real();
	Clq1_2223r_LNP = (Clq1_u*YucuL(1,2)).real();
	Clq1_2223i_LNP = (Clq1_u*YucuL(1,2)).imag();
	Clq1_2233r_LNP = (Clq1_0 + Clq1_d*YddL(2) + Clq1_u*YucuL(2,2)).real();
	Clq1_3311r_LNP = (Clq1_0 + Clq1_d*YddL(0) + Clq1_u*YucuL(0,0)).real();
	Clq1_3312r_LNP = (Clq1_u*YucuL(0,1)).real();
	Clq1_3312i_LNP = (Clq1_u*YucuL(0,1)).imag();
	Clq1_3313r_LNP = (Clq1_u*YucuL(0,2)).real();
	Clq1_3313i_LNP = (Clq1_u*YucuL(0,2)).imag();
	Clq1_3322r_LNP = (Clq1_0 + Clq1_d*YddL(1) + Clq1_u*YucuL(1,1)).real();
	Clq1_3323r_LNP = (Clq1_u*YucuL(1,2)).real();
	Clq1_3323i_LNP = (Clq1_u*YucuL(1,2)).imag();
	Clq1_3333r_LNP = (Clq1_0 + Clq1_d*YddL(2) + Clq1_u*YucuL(2,2)).real();

	Clq3_1111r_LNP = (Clq3_0 + Clq3_d*YddL(0) + Clq3_u*YucuL(0,0)).real();
	Clq3_1112r_LNP = (Clq3_u*YucuL(0,1)).real();
	Clq3_1112i_LNP = (Clq3_u*YucuL(0,1)).imag();
	Clq3_1113r_LNP = (Clq3_u*YucuL(0,2)).real();
	Clq3_1113i_LNP = (Clq3_u*YucuL(0,2)).imag();
	Clq3_1122r_LNP = (Clq3_0 + Clq3_d*YddL(1) + Clq3_u*YucuL(1,1)).real();
	Clq3_1123r_LNP = (Clq3_u*YucuL(1,2)).real();
	Clq3_1123i_LNP = (Clq3_u*YucuL(1,2)).imag();
	Clq3_1133r_LNP = (Clq3_0 + Clq3_d*YddL(2) + Clq3_u*YucuL(2,2)).real();
	Clq3_2211r_LNP = (Clq3_0 + Clq3_d*YddL(0) + Clq3_u*YucuL(0,0)).real();
	Clq3_2212r_LNP = (Clq3_u*YucuL(0,1)).real();
	Clq3_2212i_LNP = (Clq3_u*YucuL(0,1)).imag();
	Clq3_2213r_LNP = (Clq3_u*YucuL(0,2)).real();
	Clq3_2213i_LNP = (Clq3_u*YucuL(0,2)).imag();
	Clq3_2222r_LNP = (Clq3_0 + Clq3_d*YddL(1) + Clq3_u*YucuL(1,1)).real();
	Clq3_2223r_LNP = (Clq3_u*YucuL(1,2)).real();
	Clq3_2223i_LNP = (Clq3_u*YucuL(1,2)).imag();
	Clq3_2233r_LNP = (Clq3_0 + Clq3_d*YddL(2) + Clq3_u*YucuL(2,2)).real();
	Clq3_3311r_LNP = (Clq3_0 + Clq3_d*YddL(0) + Clq3_u*YucuL(0,0)).real();
	Clq3_3312r_LNP = (Clq3_u*YucuL(0,1)).real();
	Clq3_3312i_LNP = (Clq3_u*YucuL(0,1)).imag();
	Clq3_3313r_LNP = (Clq3_u*YucuL(0,2)).real();
	Clq3_3313i_LNP = (Clq3_u*YucuL(0,2)).imag();
	Clq3_3322r_LNP = (Clq3_0 + Clq3_d*YddL(1) + Clq3_u*YucuL(1,1)).real();
	Clq3_3323r_LNP = (Clq3_u*YucuL(1,2)).real();
	Clq3_3323i_LNP = (Clq3_u*YucuL(1,2)).imag();
	Clq3_3333r_LNP = (Clq3_0 + Clq3_d*YddL(2) + Clq3_u*YucuL(2,2)).real();

	Cqe_1111r_LNP = (Cqe_0 + Cqe_d*YddL(0) + Cqe_u*YucuL(0,0)).real();
	Cqe_1122r_LNP = (Cqe_0 + Cqe_d*YddL(0) + Cqe_u*YucuL(0,0)).real();
	Cqe_1133r_LNP = (Cqe_0 + Cqe_d*YddL(0) + Cqe_u*YucuL(0,0)).real();
	Cqe_1211r_LNP = (Cqe_u*YucuL(0,1)).real();
	Cqe_1211i_LNP = (Cqe_u*YucuL(0,1)).imag();
	Cqe_1222r_LNP = (Cqe_u*YucuL(0,1)).real();
	Cqe_1222i_LNP = (Cqe_u*YucuL(0,1)).imag();
	Cqe_1233r_LNP = (Cqe_u*YucuL(0,1)).real();
	Cqe_1233i_LNP = (Cqe_u*YucuL(0,1)).imag();
	Cqe_1311r_LNP = (Cqe_u*YucuL(0,2)).real();
	Cqe_1311i_LNP = (Cqe_u*YucuL(0,2)).imag();
	Cqe_1322r_LNP = (Cqe_u*YucuL(0,2)).real();
	Cqe_1322i_LNP = (Cqe_u*YucuL(0,2)).imag();
	Cqe_1333r_LNP = (Cqe_u*YucuL(0,2)).real();
	Cqe_1333i_LNP = (Cqe_u*YucuL(0,2)).imag();
	Cqe_2211r_LNP = (Cqe_0 + Cqe_d*YddL(1) + Cqe_u*YucuL(1,1)).real();
	Cqe_2222r_LNP = (Cqe_0 + Cqe_d*YddL(1) + Cqe_u*YucuL(1,1)).real();
	Cqe_2233r_LNP = (Cqe_0 + Cqe_d*YddL(1) + Cqe_u*YucuL(1,1)).real();
	Cqe_2311r_LNP = (Cqe_u*YucuL(1,2)).real();
	Cqe_2311i_LNP = (Cqe_u*YucuL(1,2)).imag();
	Cqe_2322r_LNP = (Cqe_u*YucuL(1,2)).real();
	Cqe_2322i_LNP = (Cqe_u*YucuL(1,2)).imag();
	Cqe_2333r_LNP = (Cqe_u*YucuL(1,2)).real();
	Cqe_2333i_LNP = (Cqe_u*YucuL(1,2)).imag();
	Cqe_3311r_LNP = (Cqe_0 + Cqe_d*YddL(2) + Cqe_u*YucuL(2,2)).real();
	Cqe_3322r_LNP = (Cqe_0 + Cqe_d*YddL(2) + Cqe_u*YucuL(2,2)).real();
	Cqe_3333r_LNP = (Cqe_0 + Cqe_d*YddL(2) + Cqe_u*YucuL(2,2)).real();

	Clu_1111r_LNP = (Clu_0 + Clu_u*YuucL(0,0)).real();
	Clu_1112r_LNP = (Clu_u*YuucL(0,1)).real();
	Clu_1112i_LNP = (Clu_u*YuucL(0,1)).imag();
	Clu_1113r_LNP = (Clu_u*YuucL(0,2)).real();
	Clu_1113i_LNP = (Clu_u*YuucL(0,2)).imag();
	Clu_1122r_LNP = (Clu_0 + Clu_u*YuucL(1,1)).real();
	Clu_1123r_LNP = (Clu_u*YuucL(1,2)).real();
	Clu_1123i_LNP = (Clu_u*YuucL(1,2)).imag();
	Clu_1133r_LNP = (Clu_0 + Clu_u*YuucL(2,2)).real();
	Clu_2211r_LNP = (Clu_0 + Clu_u*YuucL(0,0)).real();
	Clu_2212r_LNP = (Clu_u*YuucL(0,1)).real();
	Clu_2212i_LNP = (Clu_u*YuucL(0,1)).imag();
	Clu_2213r_LNP = (Clu_u*YuucL(0,2)).real();
	Clu_2213i_LNP = (Clu_u*YuucL(0,2)).imag();
	Clu_2222r_LNP = (Clu_0 + Clu_u*YuucL(1,1)).real();
	Clu_2223r_LNP = (Clu_u*YuucL(1,2)).real();
	Clu_2223i_LNP = (Clu_u*YuucL(1,2)).imag();
	Clu_2233r_LNP = (Clu_0 + Clu_u*YuucL(2,2)).real();
	Clu_3311r_LNP = (Clu_0 + Clu_u*YuucL(0,0)).real();
	Clu_3312r_LNP = (Clu_u*YuucL(0,1)).real();
	Clu_3312i_LNP = (Clu_u*YuucL(0,1)).imag();
	Clu_3313r_LNP = (Clu_u*YuucL(0,2)).real();
	Clu_3313i_LNP = (Clu_u*YuucL(0,2)).imag();
	Clu_3322r_LNP = (Clu_0 + Clu_u*YuucL(1,1)).real();
	Clu_3323r_LNP = (Clu_u*YuucL(1,2)).real();
	Clu_3323i_LNP = (Clu_u*YuucL(1,2)).imag();
	Clu_3333r_LNP = (Clu_0 + Clu_u*YuucL(2,2)).real();

	Ceu_1111r_LNP = (Ceu_0 + Ceu_u*YuucL(0,0)).real();
	Ceu_1112r_LNP = (Ceu_u*YuucL(0,1)).real();
	Ceu_1112i_LNP = (Ceu_u*YuucL(0,1)).imag();
	Ceu_1113r_LNP = (Ceu_u*YuucL(0,2)).real();
	Ceu_1113i_LNP = (Ceu_u*YuucL(0,2)).imag();
	Ceu_1122r_LNP = (Ceu_0 + Ceu_u*YuucL(1,1)).real();
	Ceu_1123r_LNP = (Ceu_u*YuucL(1,2)).real();
	Ceu_1123i_LNP = (Ceu_u*YuucL(1,2)).imag();
	Ceu_1133r_LNP = (Ceu_0 + Ceu_u*YuucL(2,2)).real();
	Ceu_2211r_LNP = (Ceu_0 + Ceu_u*YuucL(0,0)).real();
	Ceu_2212r_LNP = (Ceu_u*YuucL(0,1)).real();
	Ceu_2212i_LNP = (Ceu_u*YuucL(0,1)).imag();
	Ceu_2213r_LNP = (Ceu_u*YuucL(0,2)).real();
	Ceu_2213i_LNP = (Ceu_u*YuucL(0,2)).imag();
	Ceu_2222r_LNP = (Ceu_0 + Ceu_u*YuucL(1,1)).real();
	Ceu_2223r_LNP = (Ceu_u*YuucL(1,2)).real();
	Ceu_2223i_LNP = (Ceu_u*YuucL(1,2)).imag();
	Ceu_2233r_LNP = (Ceu_0 + Ceu_u*YuucL(2,2)).real();
	Ceu_3311r_LNP = (Ceu_0 + Ceu_u*YuucL(0,0)).real();
	Ceu_3312r_LNP = (Ceu_u*YuucL(0,1)).real();
	Ceu_3312i_LNP = (Ceu_u*YuucL(0,1)).imag();
	Ceu_3313r_LNP = (Ceu_u*YuucL(0,2)).real();
	Ceu_3313i_LNP = (Ceu_u*YuucL(0,2)).imag();
	Ceu_3322r_LNP = (Ceu_0 + Ceu_u*YuucL(1,1)).real();
	Ceu_3323r_LNP = (Ceu_u*YuucL(1,2)).real();
	Ceu_3323i_LNP = (Ceu_u*YuucL(1,2)).imag();
	Ceu_3333r_LNP = (Ceu_0 + Ceu_u*YuucL(2,2)).real();

	Cld_1111r_LNP = (Cld_0 + Cld_d*YddL(0)).real();
	Cld_1112r_LNP = 0;
	Cld_1112i_LNP = 0;
	Cld_1113r_LNP = 0;
	Cld_1113i_LNP = 0;
	Cld_1122r_LNP = (Cld_0 + Cld_d*YddL(1)).real();
	Cld_1123r_LNP = 0;
	Cld_1123i_LNP = 0;
	Cld_1133r_LNP = (Cld_0 + Cld_d*YddL(2)).real();
	Cld_2211r_LNP = (Cld_0 + Cld_d*YddL(0)).real();
	Cld_2212r_LNP = 0;
	Cld_2212i_LNP = 0;
	Cld_2213r_LNP = 0;
	Cld_2213i_LNP = 0;
	Cld_2222r_LNP = (Cld_0 + Cld_d*YddL(1)).real();
	Cld_2223r_LNP = 0;
	Cld_2223i_LNP = 0;
	Cld_2233r_LNP = (Cld_0 + Cld_d*YddL(2)).real();
	Cld_3311r_LNP = (Cld_0 + Cld_d*YddL(0)).real();
	Cld_3312r_LNP = 0;
	Cld_3312i_LNP = 0;
	Cld_3313r_LNP = 0;
	Cld_3313i_LNP = 0;
	Cld_3322r_LNP = (Cld_0 + Cld_d*YddL(1)).real();
	Cld_3323r_LNP = 0;
	Cld_3323i_LNP = 0;
	Cld_3333r_LNP = (Cld_0 + Cld_d*YddL(2)).real();

	Ced_1111r_LNP = (Ced_0 + Ced_d*YddL(0)).real();
	Ced_1112r_LNP = 0;
	Ced_1112i_LNP = 0;
	Ced_1113r_LNP = 0;
	Ced_1113i_LNP = 0;
	Ced_1122r_LNP = (Ced_0 + Ced_d*YddL(1)).real();
	Ced_1123r_LNP = 0;
	Ced_1123i_LNP = 0;
	Ced_1133r_LNP = (Ced_0 + Ced_d*YddL(2)).real();
	Ced_2211r_LNP = (Ced_0 + Ced_d*YddL(0)).real();
	Ced_2212r_LNP = 0;
	Ced_2212i_LNP = 0;
	Ced_2213r_LNP = 0;
	Ced_2213i_LNP = 0;
	Ced_2222r_LNP = (Ced_0 + Ced_d*YddL(1)).real();
	Ced_2223r_LNP = 0;
	Ced_2223i_LNP = 0;
	Ced_2233r_LNP = (Ced_0 + Ced_d*YddL(2)).real();
	Ced_3311r_LNP = (Ced_0 + Ced_d*YddL(0)).real();
	Ced_3312r_LNP = 0;
	Ced_3312i_LNP = 0;
	Ced_3313r_LNP = 0;
	Ced_3313i_LNP = 0;
	Ced_3322r_LNP = (Ced_0 + Ced_d*YddL(1)).real();
	Ced_3323r_LNP = 0;
	Ced_3323i_LNP = 0;
	Ced_3333r_LNP = (Ced_0 + Ced_d*YddL(2)).real();

	Cqq1_1111r_LNP = (2*Cqq1_00 + 2*Cqq1_0d*YddL(0) + 2*Cqq1_d0*YddL(0) + 2*Cqq1_dd*YddL(0)*YddL(0) + 2*Cqq1_0u*YucuL(0,0) + 2*Cqq1_u0*YucuL(0,0) + 2*Cqq1_du*YddL(0)*YucuL(0,0) + 2*Cqq1_ud*YddL(0)*YucuL(0,0) + 2*Cqq1_uu*YucuL(0,0)*YucuL(0,0)).real();
	Cqq1_1112r_LNP = (Cqq1_0u*YucuL(0,1) + Cqq1_u0*YucuL(0,1) + Cqq1_du*YddL(0)*YucuL(0,1) + Cqq1_ud*YddL(0)*YucuL(0,1) + 2*Cqq1_uu*YucuL(0,0)*YucuL(0,1)).real();
	Cqq1_1112i_LNP = (Cqq1_0u*YucuL(0,1) + Cqq1_u0*YucuL(0,1) + Cqq1_du*YddL(0)*YucuL(0,1) + Cqq1_ud*YddL(0)*YucuL(0,1) + 2*Cqq1_uu*YucuL(0,0)*YucuL(0,1)).imag();
	Cqq1_1113r_LNP = (Cqq1_0u*YucuL(0,2) + Cqq1_u0*YucuL(0,2) + Cqq1_du*YddL(0)*YucuL(0,2) + Cqq1_ud*YddL(0)*YucuL(0,2) + 2*Cqq1_uu*YucuL(0,0)*YucuL(0,2)).real();
	Cqq1_1113i_LNP = (Cqq1_0u*YucuL(0,2) + Cqq1_u0*YucuL(0,2) + Cqq1_du*YddL(0)*YucuL(0,2) + Cqq1_ud*YddL(0)*YucuL(0,2) + 2*Cqq1_uu*YucuL(0,0)*YucuL(0,2)).imag();
	Cqq1_1122r_LNP = (Cqq1_00 + Cqq1_d0*YddL(0) + Cqq1_0d*YddL(1) + Cqq1_dd*YddL(0)*YddL(1) + Cqq1_u0*YucuL(0,0) + Cqq1_ud*YddL(1)*YucuL(0,0) + Cqq1_uu*YucuL(0,1)*YucuL(1,0) + Cqq1_0u*YucuL(1,1) + Cqq1_du*YddL(0)*YucuL(1,1) + Cqq1_uu*YucuL(0,0)*YucuL(1,1)).real();
	Cqq1_1123r_LNP = (Cqq1_uu*YucuL(0,2)*YucuL(1,0) + Cqq1_0u*YucuL(1,2) + Cqq1_du*YddL(0)*YucuL(1,2) + Cqq1_uu*YucuL(0,0)*YucuL(1,2)).real();
	Cqq1_1123i_LNP = (Cqq1_uu*YucuL(0,2)*YucuL(1,0) + Cqq1_0u*YucuL(1,2) + Cqq1_du*YddL(0)*YucuL(1,2) + Cqq1_uu*YucuL(0,0)*YucuL(1,2)).imag();
	Cqq1_1133r_LNP = (Cqq1_00 + Cqq1_d0*YddL(0) + Cqq1_0d*YddL(2) + Cqq1_dd*YddL(0)*YddL(2) + Cqq1_u0*YucuL(0,0) + Cqq1_ud*YddL(2)*YucuL(0,0) + Cqq1_uu*YucuL(0,2)*YucuL(2,0) + Cqq1_0u*YucuL(2,2) + Cqq1_du*YddL(0)*YucuL(2,2) + Cqq1_uu*YucuL(0,0)*YucuL(2,2)).real();
	Cqq1_1212r_LNP = (2*Cqq1_uu*YucuL(0,1)*YucuL(0,1)).real();
	Cqq1_1212i_LNP = (2*Cqq1_uu*YucuL(0,1)*YucuL(0,1)).imag();
	Cqq1_1213r_LNP = (2*Cqq1_uu*YucuL(0,1)*YucuL(0,2)).real();
	Cqq1_1213i_LNP = (2*Cqq1_uu*YucuL(0,1)*YucuL(0,2)).imag();
	Cqq1_1221r_LNP = (Cqq1_00 + Cqq1_d0*YddL(0) + Cqq1_0d*YddL(1) + Cqq1_dd*YddL(0)*YddL(1) + Cqq1_u0*YucuL(0,0) + Cqq1_ud*YddL(1)*YucuL(0,0) + Cqq1_uu*YucuL(0,1)*YucuL(1,0) + Cqq1_0u*YucuL(1,1) + Cqq1_du*YddL(0)*YucuL(1,1) + Cqq1_uu*YucuL(0,0)*YucuL(1,1)).real();
	Cqq1_1222r_LNP = (2*Cqq1_u0*YucuL(0,1) + 2*Cqq1_ud*YddL(1)*YucuL(0,1) + 2*Cqq1_uu*YucuL(0,1)*YucuL(1,1)).real();
	Cqq1_1222i_LNP = (2*Cqq1_u0*YucuL(0,1) + 2*Cqq1_ud*YddL(1)*YucuL(0,1) + 2*Cqq1_uu*YucuL(0,1)*YucuL(1,1)).imag();
	Cqq1_1223r_LNP = (Cqq1_u0*YucuL(0,2) + Cqq1_ud*YddL(1)*YucuL(0,2) + Cqq1_uu*YucuL(0,2)*YucuL(1,1) + Cqq1_uu*YucuL(0,1)*YucuL(1,2)).real();
	Cqq1_1223i_LNP = (Cqq1_u0*YucuL(0,2) + Cqq1_ud*YddL(1)*YucuL(0,2) + Cqq1_uu*YucuL(0,2)*YucuL(1,1) + Cqq1_uu*YucuL(0,1)*YucuL(1,2)).imag();
	Cqq1_1231r_LNP = (Cqq1_uu*YucuL(0,1)*YucuL(2,0) + Cqq1_0u*YucuL(2,1) + Cqq1_du*YddL(0)*YucuL(2,1) + Cqq1_uu*YucuL(0,0)*YucuL(2,1)).real();
	Cqq1_1231i_LNP = (Cqq1_uu*YucuL(0,1)*YucuL(2,0) + Cqq1_0u*YucuL(2,1) + Cqq1_du*YddL(0)*YucuL(2,1) + Cqq1_uu*YucuL(0,0)*YucuL(2,1)).imag();
	Cqq1_1232r_LNP = (2*Cqq1_uu*YucuL(0,1)*YucuL(2,1)).real();
	Cqq1_1232i_LNP = (2*Cqq1_uu*YucuL(0,1)*YucuL(2,1)).imag();
	Cqq1_1233r_LNP = (Cqq1_u0*YucuL(0,1) + Cqq1_ud*YddL(2)*YucuL(0,1) + Cqq1_uu*YucuL(0,2)*YucuL(2,1) + Cqq1_uu*YucuL(0,1)*YucuL(2,2)).real();
	Cqq1_1233i_LNP = (Cqq1_u0*YucuL(0,1) + Cqq1_ud*YddL(2)*YucuL(0,1) + Cqq1_uu*YucuL(0,2)*YucuL(2,1) + Cqq1_uu*YucuL(0,1)*YucuL(2,2)).imag();
	Cqq1_1313r_LNP = (2*Cqq1_uu*YucuL(0,2)*YucuL(0,2)).real();
	Cqq1_1313i_LNP = (2*Cqq1_uu*YucuL(0,2)*YucuL(0,2)).imag();
	Cqq1_1322r_LNP = (Cqq1_u0*YucuL(0,2) + Cqq1_ud*YddL(1)*YucuL(0,2) + Cqq1_uu*YucuL(0,2)*YucuL(1,1) + Cqq1_uu*YucuL(0,1)*YucuL(1,2)).real();
	Cqq1_1322i_LNP = (Cqq1_u0*YucuL(0,2) + Cqq1_ud*YddL(1)*YucuL(0,2) + Cqq1_uu*YucuL(0,2)*YucuL(1,1) + Cqq1_uu*YucuL(0,1)*YucuL(1,2)).imag();
	Cqq1_1323r_LNP = (2*Cqq1_uu*YucuL(0,2)*YucuL(1,2)).real();
	Cqq1_1323i_LNP = (2*Cqq1_uu*YucuL(0,2)*YucuL(1,2)).imag();
	Cqq1_1331r_LNP = (Cqq1_00 + Cqq1_d0*YddL(0) + Cqq1_0d*YddL(2) + Cqq1_dd*YddL(0)*YddL(2) + Cqq1_u0*YucuL(0,0) + Cqq1_ud*YddL(2)*YucuL(0,0) + Cqq1_uu*YucuL(0,2)*YucuL(2,0) + Cqq1_0u*YucuL(2,2) + Cqq1_du*YddL(0)*YucuL(2,2) + Cqq1_uu*YucuL(0,0)*YucuL(2,2)).real();
	Cqq1_1332r_LNP = (Cqq1_u0*YucuL(0,1) + Cqq1_ud*YddL(2)*YucuL(0,1) + Cqq1_uu*YucuL(0,2)*YucuL(2,1) + Cqq1_uu*YucuL(0,1)*YucuL(2,2)).real();
	Cqq1_1332i_LNP = (Cqq1_u0*YucuL(0,1) + Cqq1_ud*YddL(2)*YucuL(0,1) + Cqq1_uu*YucuL(0,2)*YucuL(2,1) + Cqq1_uu*YucuL(0,1)*YucuL(2,2)).imag();
	Cqq1_1333r_LNP = (2*Cqq1_u0*YucuL(0,2) + 2*Cqq1_ud*YddL(2)*YucuL(0,2) + 2*Cqq1_uu*YucuL(0,2)*YucuL(2,2)).real();
	Cqq1_1333i_LNP = (2*Cqq1_u0*YucuL(0,2) + 2*Cqq1_ud*YddL(2)*YucuL(0,2) + 2*Cqq1_uu*YucuL(0,2)*YucuL(2,2)).imag();
	Cqq1_2222r_LNP = (2*Cqq1_00 + 2*Cqq1_0d*YddL(1) + 2*Cqq1_d0*YddL(1) + 2*Cqq1_dd*YddL(1)*YddL(1) + 2*Cqq1_0u*YucuL(1,1) + 2*Cqq1_u0*YucuL(1,1) + 2*Cqq1_du*YddL(1)*YucuL(1,1) + 2*Cqq1_ud*YddL(1)*YucuL(1,1) + 2*Cqq1_uu*YucuL(1,1)*YucuL(1,1)).real();
	Cqq1_2223r_LNP = (Cqq1_0u*YucuL(1,2) + Cqq1_u0*YucuL(1,2) + Cqq1_du*YddL(1)*YucuL(1,2) + Cqq1_ud*YddL(1)*YucuL(1,2) + 2*Cqq1_uu*YucuL(1,1)*YucuL(1,2)).real();
	Cqq1_2223i_LNP = (Cqq1_0u*YucuL(1,2) + Cqq1_u0*YucuL(1,2) + Cqq1_du*YddL(1)*YucuL(1,2) + Cqq1_ud*YddL(1)*YucuL(1,2) + 2*Cqq1_uu*YucuL(1,1)*YucuL(1,2)).imag();
	Cqq1_2233r_LNP = (Cqq1_00 + Cqq1_d0*YddL(1) + Cqq1_0d*YddL(2) + Cqq1_dd*YddL(1)*YddL(2) + Cqq1_u0*YucuL(1,1) + Cqq1_ud*YddL(2)*YucuL(1,1) + Cqq1_uu*YucuL(1,2)*YucuL(2,1) + Cqq1_0u*YucuL(2,2) + Cqq1_du*YddL(1)*YucuL(2,2) + Cqq1_uu*YucuL(1,1)*YucuL(2,2)).real();
	Cqq1_2323r_LNP = (2*Cqq1_uu*YucuL(1,2)*YucuL(1,2)).real();
	Cqq1_2323i_LNP = (2*Cqq1_uu*YucuL(1,2)*YucuL(1,2)).imag();
	Cqq1_2332r_LNP = (Cqq1_00 + Cqq1_d0*YddL(1) + Cqq1_0d*YddL(2) + Cqq1_dd*YddL(1)*YddL(2) + Cqq1_u0*YucuL(1,1) + Cqq1_ud*YddL(2)*YucuL(1,1) + Cqq1_uu*YucuL(1,2)*YucuL(2,1) + Cqq1_0u*YucuL(2,2) + Cqq1_du*YddL(1)*YucuL(2,2) + Cqq1_uu*YucuL(1,1)*YucuL(2,2)).real();
	Cqq1_2333r_LNP = (2*Cqq1_u0*YucuL(1,2) + 2*Cqq1_ud*YddL(2)*YucuL(1,2) + 2*Cqq1_uu*YucuL(1,2)*YucuL(2,2)).real();
	Cqq1_2333i_LNP = (2*Cqq1_u0*YucuL(1,2) + 2*Cqq1_ud*YddL(2)*YucuL(1,2) + 2*Cqq1_uu*YucuL(1,2)*YucuL(2,2)).imag();
	Cqq1_3333r_LNP = (2*Cqq1_00 + 2*Cqq1_0d*YddL(2) + 2*Cqq1_d0*YddL(2) + 2*Cqq1_dd*YddL(2)*YddL(2) + 2*Cqq1_0u*YucuL(2,2) + 2*Cqq1_u0*YucuL(2,2) + 2*Cqq1_du*YddL(2)*YucuL(2,2) + 2*Cqq1_ud*YddL(2)*YucuL(2,2) + 2*Cqq1_uu*YucuL(2,2)*YucuL(2,2)).real();

	Cqq3_1111r_LNP = (2*Cqq3_00 + 2*Cqq3_0d*YddL(0) + 2*Cqq3_d0*YddL(0) + 2*Cqq3_dd*YddL(0)*YddL(0) + 2*Cqq3_0u*YucuL(0,0) + 2*Cqq3_u0*YucuL(0,0) + 2*Cqq3_du*YddL(0)*YucuL(0,0) + 2*Cqq3_ud*YddL(0)*YucuL(0,0) + 2*Cqq3_uu*YucuL(0,0)*YucuL(0,0)).real();
	Cqq3_1112r_LNP = (Cqq3_0u*YucuL(0,1) + Cqq3_u0*YucuL(0,1) + Cqq3_du*YddL(0)*YucuL(0,1) + Cqq3_ud*YddL(0)*YucuL(0,1) + 2*Cqq3_uu*YucuL(0,0)*YucuL(0,1)).real();
	Cqq3_1112i_LNP = (Cqq3_0u*YucuL(0,1) + Cqq3_u0*YucuL(0,1) + Cqq3_du*YddL(0)*YucuL(0,1) + Cqq3_ud*YddL(0)*YucuL(0,1) + 2*Cqq3_uu*YucuL(0,0)*YucuL(0,1)).imag();
	Cqq3_1113r_LNP = (Cqq3_0u*YucuL(0,2) + Cqq3_u0*YucuL(0,2) + Cqq3_du*YddL(0)*YucuL(0,2) + Cqq3_ud*YddL(0)*YucuL(0,2) + 2*Cqq3_uu*YucuL(0,0)*YucuL(0,2)).real();
	Cqq3_1113i_LNP = (Cqq3_0u*YucuL(0,2) + Cqq3_u0*YucuL(0,2) + Cqq3_du*YddL(0)*YucuL(0,2) + Cqq3_ud*YddL(0)*YucuL(0,2) + 2*Cqq3_uu*YucuL(0,0)*YucuL(0,2)).imag();
	Cqq3_1122r_LNP = (Cqq3_00 + Cqq3_d0*YddL(0) + Cqq3_0d*YddL(1) + Cqq3_dd*YddL(0)*YddL(1) + Cqq3_u0*YucuL(0,0) + Cqq3_ud*YddL(1)*YucuL(0,0) + Cqq3_uu*YucuL(0,1)*YucuL(1,0) + Cqq3_0u*YucuL(1,1) + Cqq3_du*YddL(0)*YucuL(1,1) + Cqq3_uu*YucuL(0,0)*YucuL(1,1)).real();
	Cqq3_1123r_LNP = (Cqq3_uu*YucuL(0,2)*YucuL(1,0) + Cqq3_0u*YucuL(1,2) + Cqq3_du*YddL(0)*YucuL(1,2) + Cqq3_uu*YucuL(0,0)*YucuL(1,2)).real();
	Cqq3_1123i_LNP = (Cqq3_uu*YucuL(0,2)*YucuL(1,0) + Cqq3_0u*YucuL(1,2) + Cqq3_du*YddL(0)*YucuL(1,2) + Cqq3_uu*YucuL(0,0)*YucuL(1,2)).imag();
	Cqq3_1133r_LNP = (Cqq3_00 + Cqq3_d0*YddL(0) + Cqq3_0d*YddL(2) + Cqq3_dd*YddL(0)*YddL(2) + Cqq3_u0*YucuL(0,0) + Cqq3_ud*YddL(2)*YucuL(0,0) + Cqq3_uu*YucuL(0,2)*YucuL(2,0) + Cqq3_0u*YucuL(2,2) + Cqq3_du*YddL(0)*YucuL(2,2) + Cqq3_uu*YucuL(0,0)*YucuL(2,2)).real();
	Cqq3_1212r_LNP = (2*Cqq3_uu*YucuL(0,1)*YucuL(0,1)).real();
	Cqq3_1212i_LNP = (2*Cqq3_uu*YucuL(0,1)*YucuL(0,1)).imag();
	Cqq3_1213r_LNP = (2*Cqq3_uu*YucuL(0,1)*YucuL(0,2)).real();
	Cqq3_1213i_LNP = (2*Cqq3_uu*YucuL(0,1)*YucuL(0,2)).imag();
	Cqq3_1221r_LNP = (Cqq3_00 + Cqq3_d0*YddL(0) + Cqq3_0d*YddL(1) + Cqq3_dd*YddL(0)*YddL(1) + Cqq3_u0*YucuL(0,0) + Cqq3_ud*YddL(1)*YucuL(0,0) + Cqq3_uu*YucuL(0,1)*YucuL(1,0) + Cqq3_0u*YucuL(1,1) + Cqq3_du*YddL(0)*YucuL(1,1) + Cqq3_uu*YucuL(0,0)*YucuL(1,1)).real();
	Cqq3_1222r_LNP = (2*Cqq3_u0*YucuL(0,1) + 2*Cqq3_ud*YddL(1)*YucuL(0,1) + 2*Cqq3_uu*YucuL(0,1)*YucuL(1,1)).real();
	Cqq3_1222i_LNP = (2*Cqq3_u0*YucuL(0,1) + 2*Cqq3_ud*YddL(1)*YucuL(0,1) + 2*Cqq3_uu*YucuL(0,1)*YucuL(1,1)).imag();
	Cqq3_1223r_LNP = (Cqq3_u0*YucuL(0,2) + Cqq3_ud*YddL(1)*YucuL(0,2) + Cqq3_uu*YucuL(0,2)*YucuL(1,1) + Cqq3_uu*YucuL(0,1)*YucuL(1,2)).real();
	Cqq3_1223i_LNP = (Cqq3_u0*YucuL(0,2) + Cqq3_ud*YddL(1)*YucuL(0,2) + Cqq3_uu*YucuL(0,2)*YucuL(1,1) + Cqq3_uu*YucuL(0,1)*YucuL(1,2)).imag();
	Cqq3_1231r_LNP = (Cqq3_uu*YucuL(0,1)*YucuL(2,0) + Cqq3_0u*YucuL(2,1) + Cqq3_du*YddL(0)*YucuL(2,1) + Cqq3_uu*YucuL(0,0)*YucuL(2,1)).real();
	Cqq3_1231i_LNP = (Cqq3_uu*YucuL(0,1)*YucuL(2,0) + Cqq3_0u*YucuL(2,1) + Cqq3_du*YddL(0)*YucuL(2,1) + Cqq3_uu*YucuL(0,0)*YucuL(2,1)).imag();
	Cqq3_1232r_LNP = (2*Cqq3_uu*YucuL(0,1)*YucuL(2,1)).real();
	Cqq3_1232i_LNP = (2*Cqq3_uu*YucuL(0,1)*YucuL(2,1)).imag();
	Cqq3_1233r_LNP = (Cqq3_u0*YucuL(0,1) + Cqq3_ud*YddL(2)*YucuL(0,1) + Cqq3_uu*YucuL(0,2)*YucuL(2,1) + Cqq3_uu*YucuL(0,1)*YucuL(2,2)).real();
	Cqq3_1233i_LNP = (Cqq3_u0*YucuL(0,1) + Cqq3_ud*YddL(2)*YucuL(0,1) + Cqq3_uu*YucuL(0,2)*YucuL(2,1) + Cqq3_uu*YucuL(0,1)*YucuL(2,2)).imag();
	Cqq3_1313r_LNP = (2*Cqq3_uu*YucuL(0,2)*YucuL(0,2)).real();
	Cqq3_1313i_LNP = (2*Cqq3_uu*YucuL(0,2)*YucuL(0,2)).imag();
	Cqq3_1322r_LNP = (Cqq3_u0*YucuL(0,2) + Cqq3_ud*YddL(1)*YucuL(0,2) + Cqq3_uu*YucuL(0,2)*YucuL(1,1) + Cqq3_uu*YucuL(0,1)*YucuL(1,2)).real();
	Cqq3_1322i_LNP = (Cqq3_u0*YucuL(0,2) + Cqq3_ud*YddL(1)*YucuL(0,2) + Cqq3_uu*YucuL(0,2)*YucuL(1,1) + Cqq3_uu*YucuL(0,1)*YucuL(1,2)).imag();
	Cqq3_1323r_LNP = (2*Cqq3_uu*YucuL(0,2)*YucuL(1,2)).real();
	Cqq3_1323i_LNP = (2*Cqq3_uu*YucuL(0,2)*YucuL(1,2)).imag();
	Cqq3_1331r_LNP = (Cqq3_00 + Cqq3_d0*YddL(0) + Cqq3_0d*YddL(2) + Cqq3_dd*YddL(0)*YddL(2) + Cqq3_u0*YucuL(0,0) + Cqq3_ud*YddL(2)*YucuL(0,0) + Cqq3_uu*YucuL(0,2)*YucuL(2,0) + Cqq3_0u*YucuL(2,2) + Cqq3_du*YddL(0)*YucuL(2,2) + Cqq3_uu*YucuL(0,0)*YucuL(2,2)).real();
	Cqq3_1332r_LNP = (Cqq3_u0*YucuL(0,1) + Cqq3_ud*YddL(2)*YucuL(0,1) + Cqq3_uu*YucuL(0,2)*YucuL(2,1) + Cqq3_uu*YucuL(0,1)*YucuL(2,2)).real();
	Cqq3_1332i_LNP = (Cqq3_u0*YucuL(0,1) + Cqq3_ud*YddL(2)*YucuL(0,1) + Cqq3_uu*YucuL(0,2)*YucuL(2,1) + Cqq3_uu*YucuL(0,1)*YucuL(2,2)).imag();
	Cqq3_1333r_LNP = (2*Cqq3_u0*YucuL(0,2) + 2*Cqq3_ud*YddL(2)*YucuL(0,2) + 2*Cqq3_uu*YucuL(0,2)*YucuL(2,2)).real();
	Cqq3_1333i_LNP = (2*Cqq3_u0*YucuL(0,2) + 2*Cqq3_ud*YddL(2)*YucuL(0,2) + 2*Cqq3_uu*YucuL(0,2)*YucuL(2,2)).imag();
	Cqq3_2222r_LNP = (2*Cqq3_00 + 2*Cqq3_0d*YddL(1) + 2*Cqq3_d0*YddL(1) + 2*Cqq3_dd*YddL(1)*YddL(1) + 2*Cqq3_0u*YucuL(1,1) + 2*Cqq3_u0*YucuL(1,1) + 2*Cqq3_du*YddL(1)*YucuL(1,1) + 2*Cqq3_ud*YddL(1)*YucuL(1,1) + 2*Cqq3_uu*YucuL(1,1)*YucuL(1,1)).real();
	Cqq3_2223r_LNP = (Cqq3_0u*YucuL(1,2) + Cqq3_u0*YucuL(1,2) + Cqq3_du*YddL(1)*YucuL(1,2) + Cqq3_ud*YddL(1)*YucuL(1,2) + 2*Cqq3_uu*YucuL(1,1)*YucuL(1,2)).real();
	Cqq3_2223i_LNP = (Cqq3_0u*YucuL(1,2) + Cqq3_u0*YucuL(1,2) + Cqq3_du*YddL(1)*YucuL(1,2) + Cqq3_ud*YddL(1)*YucuL(1,2) + 2*Cqq3_uu*YucuL(1,1)*YucuL(1,2)).imag();
	Cqq3_2233r_LNP = (Cqq3_00 + Cqq3_d0*YddL(1) + Cqq3_0d*YddL(2) + Cqq3_dd*YddL(1)*YddL(2) + Cqq3_u0*YucuL(1,1) + Cqq3_ud*YddL(2)*YucuL(1,1) + Cqq3_uu*YucuL(1,2)*YucuL(2,1) + Cqq3_0u*YucuL(2,2) + Cqq3_du*YddL(1)*YucuL(2,2) + Cqq3_uu*YucuL(1,1)*YucuL(2,2)).real();
	Cqq3_2323r_LNP = (2*Cqq3_uu*YucuL(1,2)*YucuL(1,2)).real();
	Cqq3_2323i_LNP = (2*Cqq3_uu*YucuL(1,2)*YucuL(1,2)).imag();
	Cqq3_2332r_LNP = (Cqq3_00 + Cqq3_d0*YddL(1) + Cqq3_0d*YddL(2) + Cqq3_dd*YddL(1)*YddL(2) + Cqq3_u0*YucuL(1,1) + Cqq3_ud*YddL(2)*YucuL(1,1) + Cqq3_uu*YucuL(1,2)*YucuL(2,1) + Cqq3_0u*YucuL(2,2) + Cqq3_du*YddL(1)*YucuL(2,2) + Cqq3_uu*YucuL(1,1)*YucuL(2,2)).real();
	Cqq3_2333r_LNP = (2*Cqq3_u0*YucuL(1,2) + 2*Cqq3_ud*YddL(2)*YucuL(1,2) + 2*Cqq3_uu*YucuL(1,2)*YucuL(2,2)).real();
	Cqq3_2333i_LNP = (2*Cqq3_u0*YucuL(1,2) + 2*Cqq3_ud*YddL(2)*YucuL(1,2) + 2*Cqq3_uu*YucuL(1,2)*YucuL(2,2)).imag();
	Cqq3_3333r_LNP = (2*Cqq3_00 + 2*Cqq3_0d*YddL(2) + 2*Cqq3_d0*YddL(2) + 2*Cqq3_dd*YddL(2)*YddL(2) + 2*Cqq3_0u*YucuL(2,2) + 2*Cqq3_u0*YucuL(2,2) + 2*Cqq3_du*YddL(2)*YucuL(2,2) + 2*Cqq3_ud*YddL(2)*YucuL(2,2) + 2*Cqq3_uu*YucuL(2,2)*YucuL(2,2)).real();

	Cuu_1111r_LNP = (2*Cuu_00 + 2*Cuu_0u*YuucL(0,0) + 2*Cuu_u0*YuucL(0,0) + 2*Cuu_uu*YuucL(0,0)*YuucL(0,0)).real();
	Cuu_1112r_LNP = (Cuu_0u*YuucL(0,1) + Cuu_u0*YuucL(0,1) + 2*Cuu_uu*YuucL(0,0)*YuucL(0,1)).real();
	Cuu_1112i_LNP = (Cuu_0u*YuucL(0,1) + Cuu_u0*YuucL(0,1) + 2*Cuu_uu*YuucL(0,0)*YuucL(0,1)).imag();
	Cuu_1113r_LNP = (Cuu_0u*YuucL(0,2) + Cuu_u0*YuucL(0,2) + 2*Cuu_uu*YuucL(0,0)*YuucL(0,2)).real();
	Cuu_1113i_LNP = (Cuu_0u*YuucL(0,2) + Cuu_u0*YuucL(0,2) + 2*Cuu_uu*YuucL(0,0)*YuucL(0,2)).imag();
	Cuu_1122r_LNP = (Cuu_00 + Cuu_u0*YuucL(0,0) + Cuu_uu*YuucL(0,1)*YuucL(1,0) + Cuu_0u*YuucL(1,1) + Cuu_uu*YuucL(0,0)*YuucL(1,1)).real();
	Cuu_1123r_LNP = (Cuu_uu*YuucL(0,2)*YuucL(1,0) + Cuu_0u*YuucL(1,2) + Cuu_uu*YuucL(0,0)*YuucL(1,2)).real();
	Cuu_1123i_LNP = (Cuu_uu*YuucL(0,2)*YuucL(1,0) + Cuu_0u*YuucL(1,2) + Cuu_uu*YuucL(0,0)*YuucL(1,2)).imag();
	Cuu_1133r_LNP = (Cuu_00 + Cuu_u0*YuucL(0,0) + Cuu_uu*YuucL(0,2)*YuucL(2,0) + Cuu_0u*YuucL(2,2) + Cuu_uu*YuucL(0,0)*YuucL(2,2)).real();
	Cuu_1212r_LNP = (2*Cuu_uu*YuucL(0,1)*YuucL(0,1)).real();
	Cuu_1212i_LNP = (2*Cuu_uu*YuucL(0,1)*YuucL(0,1)).imag();
	Cuu_1213r_LNP = (2*Cuu_uu*YuucL(0,1)*YuucL(0,2)).real();
	Cuu_1213i_LNP = (2*Cuu_uu*YuucL(0,1)*YuucL(0,2)).imag();
	Cuu_1221r_LNP = (Cuu_00 + Cuu_u0*YuucL(0,0) + Cuu_uu*YuucL(0,1)*YuucL(1,0) + Cuu_0u*YuucL(1,1) + Cuu_uu*YuucL(0,0)*YuucL(1,1)).real();
	Cuu_1222r_LNP = (2*Cuu_u0*YuucL(0,1) + 2*Cuu_uu*YuucL(0,1)*YuucL(1,1)).real();
	Cuu_1222i_LNP = (2*Cuu_u0*YuucL(0,1) + 2*Cuu_uu*YuucL(0,1)*YuucL(1,1)).imag();
	Cuu_1223r_LNP = (Cuu_u0*YuucL(0,2) + Cuu_uu*YuucL(0,2)*YuucL(1,1) + Cuu_uu*YuucL(0,1)*YuucL(1,2)).real();
	Cuu_1223i_LNP = (Cuu_u0*YuucL(0,2) + Cuu_uu*YuucL(0,2)*YuucL(1,1) + Cuu_uu*YuucL(0,1)*YuucL(1,2)).imag();
	Cuu_1231r_LNP = (Cuu_uu*YuucL(0,1)*YuucL(2,0) + Cuu_0u*YuucL(2,1) + Cuu_uu*YuucL(0,0)*YuucL(2,1)).real();
	Cuu_1231i_LNP = (Cuu_uu*YuucL(0,1)*YuucL(2,0) + Cuu_0u*YuucL(2,1) + Cuu_uu*YuucL(0,0)*YuucL(2,1)).imag();
	Cuu_1232r_LNP = (2*Cuu_uu*YuucL(0,1)*YuucL(2,1)).real();
	Cuu_1232i_LNP = (2*Cuu_uu*YuucL(0,1)*YuucL(2,1)).imag();
	Cuu_1233r_LNP = (Cuu_u0*YuucL(0,1) + Cuu_uu*YuucL(0,2)*YuucL(2,1) + Cuu_uu*YuucL(0,1)*YuucL(2,2)).real();
	Cuu_1233i_LNP = (Cuu_u0*YuucL(0,1) + Cuu_uu*YuucL(0,2)*YuucL(2,1) + Cuu_uu*YuucL(0,1)*YuucL(2,2)).imag();
	Cuu_1313r_LNP = (2*Cuu_uu*YuucL(0,2)*YuucL(0,2)).real();
	Cuu_1313i_LNP = (2*Cuu_uu*YuucL(0,2)*YuucL(0,2)).imag();
	Cuu_1322r_LNP = (Cuu_u0*YuucL(0,2) + Cuu_uu*YuucL(0,2)*YuucL(1,1) + Cuu_uu*YuucL(0,1)*YuucL(1,2)).real();
	Cuu_1322i_LNP = (Cuu_u0*YuucL(0,2) + Cuu_uu*YuucL(0,2)*YuucL(1,1) + Cuu_uu*YuucL(0,1)*YuucL(1,2)).imag();
	Cuu_1323r_LNP = (2*Cuu_uu*YuucL(0,2)*YuucL(1,2)).real();
	Cuu_1323i_LNP = (2*Cuu_uu*YuucL(0,2)*YuucL(1,2)).imag();
	Cuu_1331r_LNP = (Cuu_00 + Cuu_u0*YuucL(0,0) + Cuu_uu*YuucL(0,2)*YuucL(2,0) + Cuu_0u*YuucL(2,2) + Cuu_uu*YuucL(0,0)*YuucL(2,2)).real();
	Cuu_1332r_LNP = (Cuu_u0*YuucL(0,1) + Cuu_uu*YuucL(0,2)*YuucL(2,1) + Cuu_uu*YuucL(0,1)*YuucL(2,2)).real();
	Cuu_1332i_LNP = (Cuu_u0*YuucL(0,1) + Cuu_uu*YuucL(0,2)*YuucL(2,1) + Cuu_uu*YuucL(0,1)*YuucL(2,2)).imag();
	Cuu_1333r_LNP = (2*Cuu_u0*YuucL(0,2) + 2*Cuu_uu*YuucL(0,2)*YuucL(2,2)).real();
	Cuu_1333i_LNP = (2*Cuu_u0*YuucL(0,2) + 2*Cuu_uu*YuucL(0,2)*YuucL(2,2)).imag();
	Cuu_2222r_LNP = (2*Cuu_00 + 2*Cuu_0u*YuucL(1,1) + 2*Cuu_u0*YuucL(1,1) + 2*Cuu_uu*YuucL(1,1)*YuucL(1,1)).real();
	Cuu_2223r_LNP = (Cuu_0u*YuucL(1,2) + Cuu_u0*YuucL(1,2) + 2*Cuu_uu*YuucL(1,1)*YuucL(1,2)).real();
	Cuu_2223i_LNP = (Cuu_0u*YuucL(1,2) + Cuu_u0*YuucL(1,2) + 2*Cuu_uu*YuucL(1,1)*YuucL(1,2)).imag();
	Cuu_2233r_LNP = (Cuu_00 + Cuu_u0*YuucL(1,1) + Cuu_uu*YuucL(1,2)*YuucL(2,1) + Cuu_0u*YuucL(2,2) + Cuu_uu*YuucL(1,1)*YuucL(2,2)).real();
	Cuu_2323r_LNP = (2*Cuu_uu*YuucL(1,2)*YuucL(1,2)).real();
	Cuu_2323i_LNP = (2*Cuu_uu*YuucL(1,2)*YuucL(1,2)).imag();
	Cuu_2332r_LNP = (Cuu_00 + Cuu_u0*YuucL(1,1) + Cuu_uu*YuucL(1,2)*YuucL(2,1) + Cuu_0u*YuucL(2,2) + Cuu_uu*YuucL(1,1)*YuucL(2,2)).real();
	Cuu_2333r_LNP = (2*Cuu_u0*YuucL(1,2) + 2*Cuu_uu*YuucL(1,2)*YuucL(2,2)).real();
	Cuu_2333i_LNP = (2*Cuu_u0*YuucL(1,2) + 2*Cuu_uu*YuucL(1,2)*YuucL(2,2)).imag();
	Cuu_3333r_LNP = (2*Cuu_00 + 2*Cuu_0u*YuucL(2,2) + 2*Cuu_u0*YuucL(2,2) + 2*Cuu_uu*YuucL(2,2)*YuucL(2,2)).real();

	Cdd_1111r_LNP = (2*Cdd_00 + 2*Cdd_0d*YddL(0) + 2*Cdd_d0*YddL(0) + 2*Cdd_dd*YddL(0)*YddL(0)).real();
	Cdd_1112r_LNP = 0;
	Cdd_1112i_LNP = 0;
	Cdd_1113r_LNP = 0;
	Cdd_1113i_LNP = 0;
	Cdd_1122r_LNP = (Cdd_00 + Cdd_d0*YddL(0) + Cdd_0d*YddL(1) + Cdd_dd*YddL(0)*YddL(1)).real();
	Cdd_1123r_LNP = 0;
	Cdd_1123i_LNP = 0;
	Cdd_1133r_LNP = (Cdd_00 + Cdd_d0*YddL(0) + Cdd_0d*YddL(2) + Cdd_dd*YddL(0)*YddL(2)).real();
	Cdd_1212r_LNP = 0;
	Cdd_1212i_LNP = 0;
	Cdd_1213r_LNP = 0;
	Cdd_1213i_LNP = 0;
	Cdd_1221r_LNP = (Cdd_00 + Cdd_d0*YddL(0) + Cdd_0d*YddL(1) + Cdd_dd*YddL(0)*YddL(1)).real();
	Cdd_1222r_LNP = 0;
	Cdd_1222i_LNP = 0;
	Cdd_1223r_LNP = 0;
	Cdd_1223i_LNP = 0;
	Cdd_1231r_LNP = 0;
	Cdd_1231i_LNP = 0;
	Cdd_1232r_LNP = 0;
	Cdd_1232i_LNP = 0;
	Cdd_1233r_LNP = 0;
	Cdd_1233i_LNP = 0;
	Cdd_1313r_LNP = 0;
	Cdd_1313i_LNP = 0;
	Cdd_1322r_LNP = 0;
	Cdd_1322i_LNP = 0;
	Cdd_1323r_LNP = 0;
	Cdd_1323i_LNP = 0;
	Cdd_1331r_LNP = (Cdd_00 + Cdd_d0*YddL(0) + Cdd_0d*YddL(2) + Cdd_dd*YddL(0)*YddL(2)).real();
	Cdd_1332r_LNP = 0;
	Cdd_1332i_LNP = 0;
	Cdd_1333r_LNP = 0;
	Cdd_1333i_LNP = 0;
	Cdd_2222r_LNP = (2*Cdd_00 + 2*Cdd_0d*YddL(1) + 2*Cdd_d0*YddL(1) + 2*Cdd_dd*YddL(1)*YddL(1)).real();
	Cdd_2223r_LNP = 0;
	Cdd_2223i_LNP = 0;
	Cdd_2233r_LNP = (Cdd_00 + Cdd_d0*YddL(1) + Cdd_0d*YddL(2) + Cdd_dd*YddL(1)*YddL(2)).real();
	Cdd_2323r_LNP = 0;
	Cdd_2323i_LNP = 0;
	Cdd_2332r_LNP = (Cdd_00 + Cdd_d0*YddL(1) + Cdd_0d*YddL(2) + Cdd_dd*YddL(1)*YddL(2)).real();
	Cdd_2333r_LNP = 0;
	Cdd_2333i_LNP = 0;
	Cdd_3333r_LNP = (2*Cdd_00 + 2*Cdd_0d*YddL(2) + 2*Cdd_d0*YddL(2) + 2*Cdd_dd*YddL(2)*YddL(2)).real();

	Cud1_1111r_LNP = (Cud1_00 + Cud1_0d*YddL(0) + Cud1_u0*YuucL(0,0) + Cud1_ud*YddL(0)*YuucL(0,0)).real();
	Cud1_1112r_LNP = 0;
	Cud1_1112i_LNP = 0;
	Cud1_1113r_LNP = 0;
	Cud1_1113i_LNP = 0;
	Cud1_1122r_LNP = (Cud1_00 + Cud1_0d*YddL(1) + Cud1_u0*YuucL(0,0) + Cud1_ud*YddL(1)*YuucL(0,0)).real();
	Cud1_1123r_LNP = 0;
	Cud1_1123i_LNP = 0;
	Cud1_1133r_LNP = (Cud1_00 + Cud1_0d*YddL(2) + Cud1_u0*YuucL(0,0) + Cud1_ud*YddL(2)*YuucL(0,0)).real();
	Cud1_1211r_LNP = (Cud1_u0*YuucL(0,1) + Cud1_ud*YddL(0)*YuucL(0,1)).real();
	Cud1_1211i_LNP = (Cud1_u0*YuucL(0,1) + Cud1_ud*YddL(0)*YuucL(0,1)).imag();
	Cud1_1212r_LNP = 0;
	Cud1_1212i_LNP = 0;
	Cud1_1213r_LNP = 0;
	Cud1_1213i_LNP = 0;
	Cud1_1221r_LNP = 0;
	Cud1_1221i_LNP = 0;
	Cud1_1222r_LNP = (Cud1_u0*YuucL(0,1) + Cud1_ud*YddL(1)*YuucL(0,1)).real();
	Cud1_1222i_LNP = (Cud1_u0*YuucL(0,1) + Cud1_ud*YddL(1)*YuucL(0,1)).imag();
	Cud1_1223r_LNP = 0;
	Cud1_1223i_LNP = 0;
	Cud1_1231r_LNP = 0;
	Cud1_1231i_LNP = 0;
	Cud1_1232r_LNP = 0;
	Cud1_1232i_LNP = 0;
	Cud1_1233r_LNP = (Cud1_u0*YuucL(0,1) + Cud1_ud*YddL(2)*YuucL(0,1)).real();
	Cud1_1233i_LNP = (Cud1_u0*YuucL(0,1) + Cud1_ud*YddL(2)*YuucL(0,1)).imag();
	Cud1_1311r_LNP = (Cud1_u0*YuucL(0,2) + Cud1_ud*YddL(0)*YuucL(0,2)).real();
	Cud1_1311i_LNP = (Cud1_u0*YuucL(0,2) + Cud1_ud*YddL(0)*YuucL(0,2)).imag();
	Cud1_1312r_LNP = 0;
	Cud1_1312i_LNP = 0;
	Cud1_1313r_LNP = 0;
	Cud1_1313i_LNP = 0;
	Cud1_1321r_LNP = 0;
	Cud1_1321i_LNP = 0;
	Cud1_1322r_LNP = (Cud1_u0*YuucL(0,2) + Cud1_ud*YddL(1)*YuucL(0,2)).real();
	Cud1_1322i_LNP = (Cud1_u0*YuucL(0,2) + Cud1_ud*YddL(1)*YuucL(0,2)).imag();
	Cud1_1323r_LNP = 0;
	Cud1_1323i_LNP = 0;
	Cud1_1331r_LNP = 0;
	Cud1_1331i_LNP = 0;
	Cud1_1332r_LNP = 0;
	Cud1_1332i_LNP = 0;
	Cud1_1333r_LNP = (Cud1_u0*YuucL(0,2) + Cud1_ud*YddL(2)*YuucL(0,2)).real();
	Cud1_1333i_LNP = (Cud1_u0*YuucL(0,2) + Cud1_ud*YddL(2)*YuucL(0,2)).imag();
	Cud1_2211r_LNP = (Cud1_00 + Cud1_0d*YddL(0) + Cud1_u0*YuucL(1,1) + Cud1_ud*YddL(0)*YuucL(1,1)).real();
	Cud1_2212r_LNP = 0;
	Cud1_2212i_LNP = 0;
	Cud1_2213r_LNP = 0;
	Cud1_2213i_LNP = 0;
	Cud1_2222r_LNP = (Cud1_00 + Cud1_0d*YddL(1) + Cud1_u0*YuucL(1,1) + Cud1_ud*YddL(1)*YuucL(1,1)).real();
	Cud1_2223r_LNP = 0;
	Cud1_2223i_LNP = 0;
	Cud1_2233r_LNP = (Cud1_00 + Cud1_0d*YddL(2) + Cud1_u0*YuucL(1,1) + Cud1_ud*YddL(2)*YuucL(1,1)).real();
	Cud1_2311r_LNP = (Cud1_u0*YuucL(1,2) + Cud1_ud*YddL(0)*YuucL(1,2)).real();
	Cud1_2311i_LNP = (Cud1_u0*YuucL(1,2) + Cud1_ud*YddL(0)*YuucL(1,2)).imag();
	Cud1_2312r_LNP = 0;
	Cud1_2312i_LNP = 0;
	Cud1_2313r_LNP = 0;
	Cud1_2313i_LNP = 0;
	Cud1_2321r_LNP = 0;
	Cud1_2321i_LNP = 0;
	Cud1_2322r_LNP = (Cud1_u0*YuucL(1,2) + Cud1_ud*YddL(1)*YuucL(1,2)).real();
	Cud1_2322i_LNP = (Cud1_u0*YuucL(1,2) + Cud1_ud*YddL(1)*YuucL(1,2)).imag();
	Cud1_2323r_LNP = 0;
	Cud1_2323i_LNP = 0;
	Cud1_2331r_LNP = 0;
	Cud1_2331i_LNP = 0;
	Cud1_2332r_LNP = 0;
	Cud1_2332i_LNP = 0;
	Cud1_2333r_LNP = (Cud1_u0*YuucL(1,2) + Cud1_ud*YddL(2)*YuucL(1,2)).real();
	Cud1_2333i_LNP = (Cud1_u0*YuucL(1,2) + Cud1_ud*YddL(2)*YuucL(1,2)).imag();
	Cud1_3311r_LNP = (Cud1_00 + Cud1_0d*YddL(0) + Cud1_u0*YuucL(2,2) + Cud1_ud*YddL(0)*YuucL(2,2)).real();
	Cud1_3312r_LNP = 0;
	Cud1_3312i_LNP = 0;
	Cud1_3313r_LNP = 0;
	Cud1_3313i_LNP = 0;
	Cud1_3322r_LNP = (Cud1_00 + Cud1_0d*YddL(1) + Cud1_u0*YuucL(2,2) + Cud1_ud*YddL(1)*YuucL(2,2)).real();
	Cud1_3323r_LNP = 0;
	Cud1_3323i_LNP = 0;
	Cud1_3333r_LNP = (Cud1_00 + Cud1_0d*YddL(2) + Cud1_u0*YuucL(2,2) + Cud1_ud*YddL(2)*YuucL(2,2)).real();

	Cud8_1111r_LNP = (Cud8_00 + Cud8_0d*YddL(0) + Cud8_u0*YuucL(0,0) + Cud8_ud*YddL(0)*YuucL(0,0)).real();
	Cud8_1112r_LNP = 0;
	Cud8_1112i_LNP = 0;
	Cud8_1113r_LNP = 0;
	Cud8_1113i_LNP = 0;
	Cud8_1122r_LNP = (Cud8_00 + Cud8_0d*YddL(1) + Cud8_u0*YuucL(0,0) + Cud8_ud*YddL(1)*YuucL(0,0)).real();
	Cud8_1123r_LNP = 0;
	Cud8_1123i_LNP = 0;
	Cud8_1133r_LNP = (Cud8_00 + Cud8_0d*YddL(2) + Cud8_u0*YuucL(0,0) + Cud8_ud*YddL(2)*YuucL(0,0)).real();
	Cud8_1211r_LNP = (Cud8_u0*YuucL(0,1) + Cud8_ud*YddL(0)*YuucL(0,1)).real();
	Cud8_1211i_LNP = (Cud8_u0*YuucL(0,1) + Cud8_ud*YddL(0)*YuucL(0,1)).imag();
	Cud8_1212r_LNP = 0;
	Cud8_1212i_LNP = 0;
	Cud8_1213r_LNP = 0;
	Cud8_1213i_LNP = 0;
	Cud8_1221r_LNP = 0;
	Cud8_1221i_LNP = 0;
	Cud8_1222r_LNP = (Cud8_u0*YuucL(0,1) + Cud8_ud*YddL(1)*YuucL(0,1)).real();
	Cud8_1222i_LNP = (Cud8_u0*YuucL(0,1) + Cud8_ud*YddL(1)*YuucL(0,1)).imag();
	Cud8_1223r_LNP = 0;
	Cud8_1223i_LNP = 0;
	Cud8_1231r_LNP = 0;
	Cud8_1231i_LNP = 0;
	Cud8_1232r_LNP = 0;
	Cud8_1232i_LNP = 0;
	Cud8_1233r_LNP = (Cud8_u0*YuucL(0,1) + Cud8_ud*YddL(2)*YuucL(0,1)).real();
	Cud8_1233i_LNP = (Cud8_u0*YuucL(0,1) + Cud8_ud*YddL(2)*YuucL(0,1)).imag();
	Cud8_1311r_LNP = (Cud8_u0*YuucL(0,2) + Cud8_ud*YddL(0)*YuucL(0,2)).real();
	Cud8_1311i_LNP = (Cud8_u0*YuucL(0,2) + Cud8_ud*YddL(0)*YuucL(0,2)).imag();
	Cud8_1312r_LNP = 0;
	Cud8_1312i_LNP = 0;
	Cud8_1313r_LNP = 0;
	Cud8_1313i_LNP = 0;
	Cud8_1321r_LNP = 0;
	Cud8_1321i_LNP = 0;
	Cud8_1322r_LNP = (Cud8_u0*YuucL(0,2) + Cud8_ud*YddL(1)*YuucL(0,2)).real();
	Cud8_1322i_LNP = (Cud8_u0*YuucL(0,2) + Cud8_ud*YddL(1)*YuucL(0,2)).imag();
	Cud8_1323r_LNP = 0;
	Cud8_1323i_LNP = 0;
	Cud8_1331r_LNP = 0;
	Cud8_1331i_LNP = 0;
	Cud8_1332r_LNP = 0;
	Cud8_1332i_LNP = 0;
	Cud8_1333r_LNP = (Cud8_u0*YuucL(0,2) + Cud8_ud*YddL(2)*YuucL(0,2)).real();
	Cud8_1333i_LNP = (Cud8_u0*YuucL(0,2) + Cud8_ud*YddL(2)*YuucL(0,2)).imag();
	Cud8_2211r_LNP = (Cud8_00 + Cud8_0d*YddL(0) + Cud8_u0*YuucL(1,1) + Cud8_ud*YddL(0)*YuucL(1,1)).real();
	Cud8_2212r_LNP = 0;
	Cud8_2212i_LNP = 0;
	Cud8_2213r_LNP = 0;
	Cud8_2213i_LNP = 0;
	Cud8_2222r_LNP = (Cud8_00 + Cud8_0d*YddL(1) + Cud8_u0*YuucL(1,1) + Cud8_ud*YddL(1)*YuucL(1,1)).real();
	Cud8_2223r_LNP = 0;
	Cud8_2223i_LNP = 0;
	Cud8_2233r_LNP = (Cud8_00 + Cud8_0d*YddL(2) + Cud8_u0*YuucL(1,1) + Cud8_ud*YddL(2)*YuucL(1,1)).real();
	Cud8_2311r_LNP = (Cud8_u0*YuucL(1,2) + Cud8_ud*YddL(0)*YuucL(1,2)).real();
	Cud8_2311i_LNP = (Cud8_u0*YuucL(1,2) + Cud8_ud*YddL(0)*YuucL(1,2)).imag();
	Cud8_2312r_LNP = 0;
	Cud8_2312i_LNP = 0;
	Cud8_2313r_LNP = 0;
	Cud8_2313i_LNP = 0;
	Cud8_2321r_LNP = 0;
	Cud8_2321i_LNP = 0;
	Cud8_2322r_LNP = (Cud8_u0*YuucL(1,2) + Cud8_ud*YddL(1)*YuucL(1,2)).real();
	Cud8_2322i_LNP = (Cud8_u0*YuucL(1,2) + Cud8_ud*YddL(1)*YuucL(1,2)).imag();
	Cud8_2323r_LNP = 0;
	Cud8_2323i_LNP = 0;
	Cud8_2331r_LNP = 0;
	Cud8_2331i_LNP = 0;
	Cud8_2332r_LNP = 0;
	Cud8_2332i_LNP = 0;
	Cud8_2333r_LNP = (Cud8_u0*YuucL(1,2) + Cud8_ud*YddL(2)*YuucL(1,2)).real();
	Cud8_2333i_LNP = (Cud8_u0*YuucL(1,2) + Cud8_ud*YddL(2)*YuucL(1,2)).imag();
	Cud8_3311r_LNP = (Cud8_00 + Cud8_0d*YddL(0) + Cud8_u0*YuucL(2,2) + Cud8_ud*YddL(0)*YuucL(2,2)).real();
	Cud8_3312r_LNP = 0;
	Cud8_3312i_LNP = 0;
	Cud8_3313r_LNP = 0;
	Cud8_3313i_LNP = 0;
	Cud8_3322r_LNP = (Cud8_00 + Cud8_0d*YddL(1) + Cud8_u0*YuucL(2,2) + Cud8_ud*YddL(1)*YuucL(2,2)).real();
	Cud8_3323r_LNP = 0;
	Cud8_3323i_LNP = 0;
	Cud8_3333r_LNP = (Cud8_00 + Cud8_0d*YddL(2) + Cud8_u0*YuucL(2,2) + Cud8_ud*YddL(2)*YuucL(2,2)).real();

	Cqu1_1111r_LNP = (Cqu1_00 + Cqu1_d0*YddL(0) + Cqu1_u0*YucuL(0,0) + Cqu1_0u*YuucL(0,0) + Cqu1_du*YddL(0)*YuucL(0,0) + Cqu1_uu*YucuL(0,0)*YuucL(0,0)).real();
	Cqu1_1112r_LNP = (Cqu1_0u*YuucL(0,1) + Cqu1_du*YddL(0)*YuucL(0,1) + Cqu1_uu*YucuL(0,0)*YuucL(0,1)).real();
	Cqu1_1112i_LNP = (Cqu1_0u*YuucL(0,1) + Cqu1_du*YddL(0)*YuucL(0,1) + Cqu1_uu*YucuL(0,0)*YuucL(0,1)).imag();
	Cqu1_1113r_LNP = (Cqu1_0u*YuucL(0,2) + Cqu1_du*YddL(0)*YuucL(0,2) + Cqu1_uu*YucuL(0,0)*YuucL(0,2)).real();
	Cqu1_1113i_LNP = (Cqu1_0u*YuucL(0,2) + Cqu1_du*YddL(0)*YuucL(0,2) + Cqu1_uu*YucuL(0,0)*YuucL(0,2)).imag();
	Cqu1_1122r_LNP = (Cqu1_00 + Cqu1_d0*YddL(0) + Cqu1_u0*YucuL(0,0) + Cqu1_0u*YuucL(1,1) + Cqu1_du*YddL(0)*YuucL(1,1) + Cqu1_uu*YucuL(0,0)*YuucL(1,1)).real();
	Cqu1_1123r_LNP = (Cqu1_0u*YuucL(1,2) + Cqu1_du*YddL(0)*YuucL(1,2) + Cqu1_uu*YucuL(0,0)*YuucL(1,2)).real();
	Cqu1_1123i_LNP = (Cqu1_0u*YuucL(1,2) + Cqu1_du*YddL(0)*YuucL(1,2) + Cqu1_uu*YucuL(0,0)*YuucL(1,2)).imag();
	Cqu1_1133r_LNP = (Cqu1_00 + Cqu1_d0*YddL(0) + Cqu1_u0*YucuL(0,0) + Cqu1_0u*YuucL(2,2) + Cqu1_du*YddL(0)*YuucL(2,2) + Cqu1_uu*YucuL(0,0)*YuucL(2,2)).real();
	Cqu1_1211r_LNP = (Cqu1_u0*YucuL(0,1) + Cqu1_uu*YucuL(0,1)*YuucL(0,0)).real();
	Cqu1_1211i_LNP = (Cqu1_u0*YucuL(0,1) + Cqu1_uu*YucuL(0,1)*YuucL(0,0)).imag();
	Cqu1_1212r_LNP = (Cqu1_uu*YucuL(0,1)*YuucL(0,1)).real();
	Cqu1_1212i_LNP = (Cqu1_uu*YucuL(0,1)*YuucL(0,1)).imag();
	Cqu1_1213r_LNP = (Cqu1_uu*YucuL(0,1)*YuucL(0,2)).real();
	Cqu1_1213i_LNP = (Cqu1_uu*YucuL(0,1)*YuucL(0,2)).imag();
	Cqu1_1221r_LNP = (Cqu1_uu*YucuL(0,1)*YuucL(1,0)).real();
	Cqu1_1221i_LNP = (Cqu1_uu*YucuL(0,1)*YuucL(1,0)).imag();
	Cqu1_1222r_LNP = (Cqu1_u0*YucuL(0,1) + Cqu1_uu*YucuL(0,1)*YuucL(1,1)).real();
	Cqu1_1222i_LNP = (Cqu1_u0*YucuL(0,1) + Cqu1_uu*YucuL(0,1)*YuucL(1,1)).imag();
	Cqu1_1223r_LNP = (Cqu1_uu*YucuL(0,1)*YuucL(1,2)).real();
	Cqu1_1223i_LNP = (Cqu1_uu*YucuL(0,1)*YuucL(1,2)).imag();
	Cqu1_1231r_LNP = (Cqu1_uu*YucuL(0,1)*YuucL(2,0)).real();
	Cqu1_1231i_LNP = (Cqu1_uu*YucuL(0,1)*YuucL(2,0)).imag();
	Cqu1_1232r_LNP = (Cqu1_uu*YucuL(0,1)*YuucL(2,1)).real();
	Cqu1_1232i_LNP = (Cqu1_uu*YucuL(0,1)*YuucL(2,1)).imag();
	Cqu1_1233r_LNP = (Cqu1_u0*YucuL(0,1) + Cqu1_uu*YucuL(0,1)*YuucL(2,2)).real();
	Cqu1_1233i_LNP = (Cqu1_u0*YucuL(0,1) + Cqu1_uu*YucuL(0,1)*YuucL(2,2)).imag();
	Cqu1_1311r_LNP = (Cqu1_u0*YucuL(0,2) + Cqu1_uu*YucuL(0,2)*YuucL(0,0)).real();
	Cqu1_1311i_LNP = (Cqu1_u0*YucuL(0,2) + Cqu1_uu*YucuL(0,2)*YuucL(0,0)).imag();
	Cqu1_1312r_LNP = (Cqu1_uu*YucuL(0,2)*YuucL(0,1)).real();
	Cqu1_1312i_LNP = (Cqu1_uu*YucuL(0,2)*YuucL(0,1)).imag();
	Cqu1_1313r_LNP = (Cqu1_uu*YucuL(0,2)*YuucL(0,2)).real();
	Cqu1_1313i_LNP = (Cqu1_uu*YucuL(0,2)*YuucL(0,2)).imag();
	Cqu1_1321r_LNP = (Cqu1_uu*YucuL(0,2)*YuucL(1,0)).real();
	Cqu1_1321i_LNP = (Cqu1_uu*YucuL(0,2)*YuucL(1,0)).imag();
	Cqu1_1322r_LNP = (Cqu1_u0*YucuL(0,2) + Cqu1_uu*YucuL(0,2)*YuucL(1,1)).real();
	Cqu1_1322i_LNP = (Cqu1_u0*YucuL(0,2) + Cqu1_uu*YucuL(0,2)*YuucL(1,1)).imag();
	Cqu1_1323r_LNP = (Cqu1_uu*YucuL(0,2)*YuucL(1,2)).real();
	Cqu1_1323i_LNP = (Cqu1_uu*YucuL(0,2)*YuucL(1,2)).imag();
	Cqu1_1331r_LNP = (Cqu1_uu*YucuL(0,2)*YuucL(2,0)).real();
	Cqu1_1331i_LNP = (Cqu1_uu*YucuL(0,2)*YuucL(2,0)).imag();
	Cqu1_1332r_LNP = (Cqu1_uu*YucuL(0,2)*YuucL(2,1)).real();
	Cqu1_1332i_LNP = (Cqu1_uu*YucuL(0,2)*YuucL(2,1)).imag();
	Cqu1_1333r_LNP = (Cqu1_u0*YucuL(0,2) + Cqu1_uu*YucuL(0,2)*YuucL(2,2)).real();
	Cqu1_1333i_LNP = (Cqu1_u0*YucuL(0,2) + Cqu1_uu*YucuL(0,2)*YuucL(2,2)).imag();
	Cqu1_2211r_LNP = (Cqu1_00 + Cqu1_d0*YddL(1) + Cqu1_u0*YucuL(1,1) + Cqu1_0u*YuucL(0,0) + Cqu1_du*YddL(1)*YuucL(0,0) + Cqu1_uu*YucuL(1,1)*YuucL(0,0)).real();
	Cqu1_2212r_LNP = (Cqu1_0u*YuucL(0,1) + Cqu1_du*YddL(1)*YuucL(0,1) + Cqu1_uu*YucuL(1,1)*YuucL(0,1)).real();
	Cqu1_2212i_LNP = (Cqu1_0u*YuucL(0,1) + Cqu1_du*YddL(1)*YuucL(0,1) + Cqu1_uu*YucuL(1,1)*YuucL(0,1)).imag();
	Cqu1_2213r_LNP = (Cqu1_0u*YuucL(0,2) + Cqu1_du*YddL(1)*YuucL(0,2) + Cqu1_uu*YucuL(1,1)*YuucL(0,2)).real();
	Cqu1_2213i_LNP = (Cqu1_0u*YuucL(0,2) + Cqu1_du*YddL(1)*YuucL(0,2) + Cqu1_uu*YucuL(1,1)*YuucL(0,2)).imag();
	Cqu1_2222r_LNP = (Cqu1_00 + Cqu1_d0*YddL(1) + Cqu1_u0*YucuL(1,1) + Cqu1_0u*YuucL(1,1) + Cqu1_du*YddL(1)*YuucL(1,1) + Cqu1_uu*YucuL(1,1)*YuucL(1,1)).real();
	Cqu1_2223r_LNP = (Cqu1_0u*YuucL(1,2) + Cqu1_du*YddL(1)*YuucL(1,2) + Cqu1_uu*YucuL(1,1)*YuucL(1,2)).real();
	Cqu1_2223i_LNP = (Cqu1_0u*YuucL(1,2) + Cqu1_du*YddL(1)*YuucL(1,2) + Cqu1_uu*YucuL(1,1)*YuucL(1,2)).imag();
	Cqu1_2233r_LNP = (Cqu1_00 + Cqu1_d0*YddL(1) + Cqu1_u0*YucuL(1,1) + Cqu1_0u*YuucL(2,2) + Cqu1_du*YddL(1)*YuucL(2,2) + Cqu1_uu*YucuL(1,1)*YuucL(2,2)).real();
	Cqu1_2311r_LNP = (Cqu1_u0*YucuL(1,2) + Cqu1_uu*YucuL(1,2)*YuucL(0,0)).real();
	Cqu1_2311i_LNP = (Cqu1_u0*YucuL(1,2) + Cqu1_uu*YucuL(1,2)*YuucL(0,0)).imag();
	Cqu1_2312r_LNP = (Cqu1_uu*YucuL(1,2)*YuucL(0,1)).real();
	Cqu1_2312i_LNP = (Cqu1_uu*YucuL(1,2)*YuucL(0,1)).imag();
	Cqu1_2313r_LNP = (Cqu1_uu*YucuL(1,2)*YuucL(0,2)).real();
	Cqu1_2313i_LNP = (Cqu1_uu*YucuL(1,2)*YuucL(0,2)).imag();
	Cqu1_2321r_LNP = (Cqu1_uu*YucuL(1,2)*YuucL(1,0)).real();
	Cqu1_2321i_LNP = (Cqu1_uu*YucuL(1,2)*YuucL(1,0)).imag();
	Cqu1_2322r_LNP = (Cqu1_u0*YucuL(1,2) + Cqu1_uu*YucuL(1,2)*YuucL(1,1)).real();
	Cqu1_2322i_LNP = (Cqu1_u0*YucuL(1,2) + Cqu1_uu*YucuL(1,2)*YuucL(1,1)).imag();
	Cqu1_2323r_LNP = (Cqu1_uu*YucuL(1,2)*YuucL(1,2)).real();
	Cqu1_2323i_LNP = (Cqu1_uu*YucuL(1,2)*YuucL(1,2)).imag();
	Cqu1_2331r_LNP = (Cqu1_uu*YucuL(1,2)*YuucL(2,0)).real();
	Cqu1_2331i_LNP = (Cqu1_uu*YucuL(1,2)*YuucL(2,0)).imag();
	Cqu1_2332r_LNP = (Cqu1_uu*YucuL(1,2)*YuucL(2,1)).real();
	Cqu1_2332i_LNP = (Cqu1_uu*YucuL(1,2)*YuucL(2,1)).imag();
	Cqu1_2333r_LNP = (Cqu1_u0*YucuL(1,2) + Cqu1_uu*YucuL(1,2)*YuucL(2,2)).real();
	Cqu1_2333i_LNP = (Cqu1_u0*YucuL(1,2) + Cqu1_uu*YucuL(1,2)*YuucL(2,2)).imag();
	Cqu1_3311r_LNP = (Cqu1_00 + Cqu1_d0*YddL(2) + Cqu1_u0*YucuL(2,2) + Cqu1_0u*YuucL(0,0) + Cqu1_du*YddL(2)*YuucL(0,0) + Cqu1_uu*YucuL(2,2)*YuucL(0,0)).real();
	Cqu1_3312r_LNP = (Cqu1_0u*YuucL(0,1) + Cqu1_du*YddL(2)*YuucL(0,1) + Cqu1_uu*YucuL(2,2)*YuucL(0,1)).real();
	Cqu1_3312i_LNP = (Cqu1_0u*YuucL(0,1) + Cqu1_du*YddL(2)*YuucL(0,1) + Cqu1_uu*YucuL(2,2)*YuucL(0,1)).imag();
	Cqu1_3313r_LNP = (Cqu1_0u*YuucL(0,2) + Cqu1_du*YddL(2)*YuucL(0,2) + Cqu1_uu*YucuL(2,2)*YuucL(0,2)).real();
	Cqu1_3313i_LNP = (Cqu1_0u*YuucL(0,2) + Cqu1_du*YddL(2)*YuucL(0,2) + Cqu1_uu*YucuL(2,2)*YuucL(0,2)).imag();
	Cqu1_3322r_LNP = (Cqu1_00 + Cqu1_d0*YddL(2) + Cqu1_u0*YucuL(2,2) + Cqu1_0u*YuucL(1,1) + Cqu1_du*YddL(2)*YuucL(1,1) + Cqu1_uu*YucuL(2,2)*YuucL(1,1)).real();
	Cqu1_3323r_LNP = (Cqu1_0u*YuucL(1,2) + Cqu1_du*YddL(2)*YuucL(1,2) + Cqu1_uu*YucuL(2,2)*YuucL(1,2)).real();
	Cqu1_3323i_LNP = (Cqu1_0u*YuucL(1,2) + Cqu1_du*YddL(2)*YuucL(1,2) + Cqu1_uu*YucuL(2,2)*YuucL(1,2)).imag();
	Cqu1_3333r_LNP = (Cqu1_00 + Cqu1_d0*YddL(2) + Cqu1_u0*YucuL(2,2) + Cqu1_0u*YuucL(2,2) + Cqu1_du*YddL(2)*YuucL(2,2) + Cqu1_uu*YucuL(2,2)*YuucL(2,2)).real();

	Cqu8_1111r_LNP = (Cqu8_00 + Cqu8_d0*YddL(0) + Cqu8_u0*YucuL(0,0) + Cqu8_0u*YuucL(0,0) + Cqu8_du*YddL(0)*YuucL(0,0) + Cqu8_uu*YucuL(0,0)*YuucL(0,0)).real();
	Cqu8_1112r_LNP = (Cqu8_0u*YuucL(0,1) + Cqu8_du*YddL(0)*YuucL(0,1) + Cqu8_uu*YucuL(0,0)*YuucL(0,1)).real();
	Cqu8_1112i_LNP = (Cqu8_0u*YuucL(0,1) + Cqu8_du*YddL(0)*YuucL(0,1) + Cqu8_uu*YucuL(0,0)*YuucL(0,1)).imag();
	Cqu8_1113r_LNP = (Cqu8_0u*YuucL(0,2) + Cqu8_du*YddL(0)*YuucL(0,2) + Cqu8_uu*YucuL(0,0)*YuucL(0,2)).real();
	Cqu8_1113i_LNP = (Cqu8_0u*YuucL(0,2) + Cqu8_du*YddL(0)*YuucL(0,2) + Cqu8_uu*YucuL(0,0)*YuucL(0,2)).imag();
	Cqu8_1122r_LNP = (Cqu8_00 + Cqu8_d0*YddL(0) + Cqu8_u0*YucuL(0,0) + Cqu8_0u*YuucL(1,1) + Cqu8_du*YddL(0)*YuucL(1,1) + Cqu8_uu*YucuL(0,0)*YuucL(1,1)).real();
	Cqu8_1123r_LNP = (Cqu8_0u*YuucL(1,2) + Cqu8_du*YddL(0)*YuucL(1,2) + Cqu8_uu*YucuL(0,0)*YuucL(1,2)).real();
	Cqu8_1123i_LNP = (Cqu8_0u*YuucL(1,2) + Cqu8_du*YddL(0)*YuucL(1,2) + Cqu8_uu*YucuL(0,0)*YuucL(1,2)).imag();
	Cqu8_1133r_LNP = (Cqu8_00 + Cqu8_d0*YddL(0) + Cqu8_u0*YucuL(0,0) + Cqu8_0u*YuucL(2,2) + Cqu8_du*YddL(0)*YuucL(2,2) + Cqu8_uu*YucuL(0,0)*YuucL(2,2)).real();
	Cqu8_1211r_LNP = (Cqu8_u0*YucuL(0,1) + Cqu8_uu*YucuL(0,1)*YuucL(0,0)).real();
	Cqu8_1211i_LNP = (Cqu8_u0*YucuL(0,1) + Cqu8_uu*YucuL(0,1)*YuucL(0,0)).imag();
	Cqu8_1212r_LNP = (Cqu8_uu*YucuL(0,1)*YuucL(0,1)).real();
	Cqu8_1212i_LNP = (Cqu8_uu*YucuL(0,1)*YuucL(0,1)).imag();
	Cqu8_1213r_LNP = (Cqu8_uu*YucuL(0,1)*YuucL(0,2)).real();
	Cqu8_1213i_LNP = (Cqu8_uu*YucuL(0,1)*YuucL(0,2)).imag();
	Cqu8_1221r_LNP = (Cqu8_uu*YucuL(0,1)*YuucL(1,0)).real();
	Cqu8_1221i_LNP = (Cqu8_uu*YucuL(0,1)*YuucL(1,0)).imag();
	Cqu8_1222r_LNP = (Cqu8_u0*YucuL(0,1) + Cqu8_uu*YucuL(0,1)*YuucL(1,1)).real();
	Cqu8_1222i_LNP = (Cqu8_u0*YucuL(0,1) + Cqu8_uu*YucuL(0,1)*YuucL(1,1)).imag();
	Cqu8_1223r_LNP = (Cqu8_uu*YucuL(0,1)*YuucL(1,2)).real();
	Cqu8_1223i_LNP = (Cqu8_uu*YucuL(0,1)*YuucL(1,2)).imag();
	Cqu8_1231r_LNP = (Cqu8_uu*YucuL(0,1)*YuucL(2,0)).real();
	Cqu8_1231i_LNP = (Cqu8_uu*YucuL(0,1)*YuucL(2,0)).imag();
	Cqu8_1232r_LNP = (Cqu8_uu*YucuL(0,1)*YuucL(2,1)).real();
	Cqu8_1232i_LNP = (Cqu8_uu*YucuL(0,1)*YuucL(2,1)).imag();
	Cqu8_1233r_LNP = (Cqu8_u0*YucuL(0,1) + Cqu8_uu*YucuL(0,1)*YuucL(2,2)).real();
	Cqu8_1233i_LNP = (Cqu8_u0*YucuL(0,1) + Cqu8_uu*YucuL(0,1)*YuucL(2,2)).imag();
	Cqu8_1311r_LNP = (Cqu8_u0*YucuL(0,2) + Cqu8_uu*YucuL(0,2)*YuucL(0,0)).real();
	Cqu8_1311i_LNP = (Cqu8_u0*YucuL(0,2) + Cqu8_uu*YucuL(0,2)*YuucL(0,0)).imag();
	Cqu8_1312r_LNP = (Cqu8_uu*YucuL(0,2)*YuucL(0,1)).real();
	Cqu8_1312i_LNP = (Cqu8_uu*YucuL(0,2)*YuucL(0,1)).imag();
	Cqu8_1313r_LNP = (Cqu8_uu*YucuL(0,2)*YuucL(0,2)).real();
	Cqu8_1313i_LNP = (Cqu8_uu*YucuL(0,2)*YuucL(0,2)).imag();
	Cqu8_1321r_LNP = (Cqu8_uu*YucuL(0,2)*YuucL(1,0)).real();
	Cqu8_1321i_LNP = (Cqu8_uu*YucuL(0,2)*YuucL(1,0)).imag();
	Cqu8_1322r_LNP = (Cqu8_u0*YucuL(0,2) + Cqu8_uu*YucuL(0,2)*YuucL(1,1)).real();
	Cqu8_1322i_LNP = (Cqu8_u0*YucuL(0,2) + Cqu8_uu*YucuL(0,2)*YuucL(1,1)).imag();
	Cqu8_1323r_LNP = (Cqu8_uu*YucuL(0,2)*YuucL(1,2)).real();
	Cqu8_1323i_LNP = (Cqu8_uu*YucuL(0,2)*YuucL(1,2)).imag();
	Cqu8_1331r_LNP = (Cqu8_uu*YucuL(0,2)*YuucL(2,0)).real();
	Cqu8_1331i_LNP = (Cqu8_uu*YucuL(0,2)*YuucL(2,0)).imag();
	Cqu8_1332r_LNP = (Cqu8_uu*YucuL(0,2)*YuucL(2,1)).real();
	Cqu8_1332i_LNP = (Cqu8_uu*YucuL(0,2)*YuucL(2,1)).imag();
	Cqu8_1333r_LNP = (Cqu8_u0*YucuL(0,2) + Cqu8_uu*YucuL(0,2)*YuucL(2,2)).real();
	Cqu8_1333i_LNP = (Cqu8_u0*YucuL(0,2) + Cqu8_uu*YucuL(0,2)*YuucL(2,2)).imag();
	Cqu8_2211r_LNP = (Cqu8_00 + Cqu8_d0*YddL(1) + Cqu8_u0*YucuL(1,1) + Cqu8_0u*YuucL(0,0) + Cqu8_du*YddL(1)*YuucL(0,0) + Cqu8_uu*YucuL(1,1)*YuucL(0,0)).real();
	Cqu8_2212r_LNP = (Cqu8_0u*YuucL(0,1) + Cqu8_du*YddL(1)*YuucL(0,1) + Cqu8_uu*YucuL(1,1)*YuucL(0,1)).real();
	Cqu8_2212i_LNP = (Cqu8_0u*YuucL(0,1) + Cqu8_du*YddL(1)*YuucL(0,1) + Cqu8_uu*YucuL(1,1)*YuucL(0,1)).imag();
	Cqu8_2213r_LNP = (Cqu8_0u*YuucL(0,2) + Cqu8_du*YddL(1)*YuucL(0,2) + Cqu8_uu*YucuL(1,1)*YuucL(0,2)).real();
	Cqu8_2213i_LNP = (Cqu8_0u*YuucL(0,2) + Cqu8_du*YddL(1)*YuucL(0,2) + Cqu8_uu*YucuL(1,1)*YuucL(0,2)).imag();
	Cqu8_2222r_LNP = (Cqu8_00 + Cqu8_d0*YddL(1) + Cqu8_u0*YucuL(1,1) + Cqu8_0u*YuucL(1,1) + Cqu8_du*YddL(1)*YuucL(1,1) + Cqu8_uu*YucuL(1,1)*YuucL(1,1)).real();
	Cqu8_2223r_LNP = (Cqu8_0u*YuucL(1,2) + Cqu8_du*YddL(1)*YuucL(1,2) + Cqu8_uu*YucuL(1,1)*YuucL(1,2)).real();
	Cqu8_2223i_LNP = (Cqu8_0u*YuucL(1,2) + Cqu8_du*YddL(1)*YuucL(1,2) + Cqu8_uu*YucuL(1,1)*YuucL(1,2)).imag();
	Cqu8_2233r_LNP = (Cqu8_00 + Cqu8_d0*YddL(1) + Cqu8_u0*YucuL(1,1) + Cqu8_0u*YuucL(2,2) + Cqu8_du*YddL(1)*YuucL(2,2) + Cqu8_uu*YucuL(1,1)*YuucL(2,2)).real();
	Cqu8_2311r_LNP = (Cqu8_u0*YucuL(1,2) + Cqu8_uu*YucuL(1,2)*YuucL(0,0)).real();
	Cqu8_2311i_LNP = (Cqu8_u0*YucuL(1,2) + Cqu8_uu*YucuL(1,2)*YuucL(0,0)).imag();
	Cqu8_2312r_LNP = (Cqu8_uu*YucuL(1,2)*YuucL(0,1)).real();
	Cqu8_2312i_LNP = (Cqu8_uu*YucuL(1,2)*YuucL(0,1)).imag();
	Cqu8_2313r_LNP = (Cqu8_uu*YucuL(1,2)*YuucL(0,2)).real();
	Cqu8_2313i_LNP = (Cqu8_uu*YucuL(1,2)*YuucL(0,2)).imag();
	Cqu8_2321r_LNP = (Cqu8_uu*YucuL(1,2)*YuucL(1,0)).real();
	Cqu8_2321i_LNP = (Cqu8_uu*YucuL(1,2)*YuucL(1,0)).imag();
	Cqu8_2322r_LNP = (Cqu8_u0*YucuL(1,2) + Cqu8_uu*YucuL(1,2)*YuucL(1,1)).real();
	Cqu8_2322i_LNP = (Cqu8_u0*YucuL(1,2) + Cqu8_uu*YucuL(1,2)*YuucL(1,1)).imag();
	Cqu8_2323r_LNP = (Cqu8_uu*YucuL(1,2)*YuucL(1,2)).real();
	Cqu8_2323i_LNP = (Cqu8_uu*YucuL(1,2)*YuucL(1,2)).imag();
	Cqu8_2331r_LNP = (Cqu8_uu*YucuL(1,2)*YuucL(2,0)).real();
	Cqu8_2331i_LNP = (Cqu8_uu*YucuL(1,2)*YuucL(2,0)).imag();
	Cqu8_2332r_LNP = (Cqu8_uu*YucuL(1,2)*YuucL(2,1)).real();
	Cqu8_2332i_LNP = (Cqu8_uu*YucuL(1,2)*YuucL(2,1)).imag();
	Cqu8_2333r_LNP = (Cqu8_u0*YucuL(1,2) + Cqu8_uu*YucuL(1,2)*YuucL(2,2)).real();
	Cqu8_2333i_LNP = (Cqu8_u0*YucuL(1,2) + Cqu8_uu*YucuL(1,2)*YuucL(2,2)).imag();
	Cqu8_3311r_LNP = (Cqu8_00 + Cqu8_d0*YddL(2) + Cqu8_u0*YucuL(2,2) + Cqu8_0u*YuucL(0,0) + Cqu8_du*YddL(2)*YuucL(0,0) + Cqu8_uu*YucuL(2,2)*YuucL(0,0)).real();
	Cqu8_3312r_LNP = (Cqu8_0u*YuucL(0,1) + Cqu8_du*YddL(2)*YuucL(0,1) + Cqu8_uu*YucuL(2,2)*YuucL(0,1)).real();
	Cqu8_3312i_LNP = (Cqu8_0u*YuucL(0,1) + Cqu8_du*YddL(2)*YuucL(0,1) + Cqu8_uu*YucuL(2,2)*YuucL(0,1)).imag();
	Cqu8_3313r_LNP = (Cqu8_0u*YuucL(0,2) + Cqu8_du*YddL(2)*YuucL(0,2) + Cqu8_uu*YucuL(2,2)*YuucL(0,2)).real();
	Cqu8_3313i_LNP = (Cqu8_0u*YuucL(0,2) + Cqu8_du*YddL(2)*YuucL(0,2) + Cqu8_uu*YucuL(2,2)*YuucL(0,2)).imag();
	Cqu8_3322r_LNP = (Cqu8_00 + Cqu8_d0*YddL(2) + Cqu8_u0*YucuL(2,2) + Cqu8_0u*YuucL(1,1) + Cqu8_du*YddL(2)*YuucL(1,1) + Cqu8_uu*YucuL(2,2)*YuucL(1,1)).real();
	Cqu8_3323r_LNP = (Cqu8_0u*YuucL(1,2) + Cqu8_du*YddL(2)*YuucL(1,2) + Cqu8_uu*YucuL(2,2)*YuucL(1,2)).real();
	Cqu8_3323i_LNP = (Cqu8_0u*YuucL(1,2) + Cqu8_du*YddL(2)*YuucL(1,2) + Cqu8_uu*YucuL(2,2)*YuucL(1,2)).imag();
	Cqu8_3333r_LNP = (Cqu8_00 + Cqu8_d0*YddL(2) + Cqu8_u0*YucuL(2,2) + Cqu8_0u*YuucL(2,2) + Cqu8_du*YddL(2)*YuucL(2,2) + Cqu8_uu*YucuL(2,2)*YuucL(2,2)).real();

	Cqd1_1111r_LNP = (Cqd1_00 + Cqd1_0d*YddL(0) + Cqd1_d0*YddL(0) + Cqd1_dd*YddL(0)*YddL(0) + Cqd1_u0*YucuL(0,0) + Cqd1_ud*YddL(0)*YucuL(0,0)).real();
	Cqd1_1112r_LNP = 0;
	Cqd1_1112i_LNP = 0;
	Cqd1_1113r_LNP = 0;
	Cqd1_1113i_LNP = 0;
	Cqd1_1122r_LNP = (Cqd1_00 + Cqd1_d0*YddL(0) + Cqd1_0d*YddL(1) + Cqd1_dd*YddL(0)*YddL(1) + Cqd1_u0*YucuL(0,0) + Cqd1_ud*YddL(1)*YucuL(0,0)).real();
	Cqd1_1123r_LNP = 0;
	Cqd1_1123i_LNP = 0;
	Cqd1_1133r_LNP = (Cqd1_00 + Cqd1_d0*YddL(0) + Cqd1_0d*YddL(2) + Cqd1_dd*YddL(0)*YddL(2) + Cqd1_u0*YucuL(0,0) + Cqd1_ud*YddL(2)*YucuL(0,0)).real();
	Cqd1_1211r_LNP = (Cqd1_u0*YucuL(0,1) + Cqd1_ud*YddL(0)*YucuL(0,1)).real();
	Cqd1_1211i_LNP = (Cqd1_u0*YucuL(0,1) + Cqd1_ud*YddL(0)*YucuL(0,1)).imag();
	Cqd1_1212r_LNP = 0;
	Cqd1_1212i_LNP = 0;
	Cqd1_1213r_LNP = 0;
	Cqd1_1213i_LNP = 0;
	Cqd1_1221r_LNP = 0;
	Cqd1_1221i_LNP = 0;
	Cqd1_1222r_LNP = (Cqd1_u0*YucuL(0,1) + Cqd1_ud*YddL(1)*YucuL(0,1)).real();
	Cqd1_1222i_LNP = (Cqd1_u0*YucuL(0,1) + Cqd1_ud*YddL(1)*YucuL(0,1)).imag();
	Cqd1_1223r_LNP = 0;
	Cqd1_1223i_LNP = 0;
	Cqd1_1231r_LNP = 0;
	Cqd1_1231i_LNP = 0;
	Cqd1_1232r_LNP = 0;
	Cqd1_1232i_LNP = 0;
	Cqd1_1233r_LNP = (Cqd1_u0*YucuL(0,1) + Cqd1_ud*YddL(2)*YucuL(0,1)).real();
	Cqd1_1233i_LNP = (Cqd1_u0*YucuL(0,1) + Cqd1_ud*YddL(2)*YucuL(0,1)).imag();
	Cqd1_1311r_LNP = (Cqd1_u0*YucuL(0,2) + Cqd1_ud*YddL(0)*YucuL(0,2)).real();
	Cqd1_1311i_LNP = (Cqd1_u0*YucuL(0,2) + Cqd1_ud*YddL(0)*YucuL(0,2)).imag();
	Cqd1_1312r_LNP = 0;
	Cqd1_1312i_LNP = 0;
	Cqd1_1313r_LNP = 0;
	Cqd1_1313i_LNP = 0;
	Cqd1_1321r_LNP = 0;
	Cqd1_1321i_LNP = 0;
	Cqd1_1322r_LNP = (Cqd1_u0*YucuL(0,2) + Cqd1_ud*YddL(1)*YucuL(0,2)).real();
	Cqd1_1322i_LNP = (Cqd1_u0*YucuL(0,2) + Cqd1_ud*YddL(1)*YucuL(0,2)).imag();
	Cqd1_1323r_LNP = 0;
	Cqd1_1323i_LNP = 0;
	Cqd1_1331r_LNP = 0;
	Cqd1_1331i_LNP = 0;
	Cqd1_1332r_LNP = 0;
	Cqd1_1332i_LNP = 0;
	Cqd1_1333r_LNP = (Cqd1_u0*YucuL(0,2) + Cqd1_ud*YddL(2)*YucuL(0,2)).real();
	Cqd1_1333i_LNP = (Cqd1_u0*YucuL(0,2) + Cqd1_ud*YddL(2)*YucuL(0,2)).imag();
	Cqd1_2211r_LNP = (Cqd1_00 + Cqd1_0d*YddL(0) + Cqd1_d0*YddL(1) + Cqd1_dd*YddL(0)*YddL(1) + Cqd1_u0*YucuL(1,1) + Cqd1_ud*YddL(0)*YucuL(1,1)).real();
	Cqd1_2212r_LNP = 0;
	Cqd1_2212i_LNP = 0;
	Cqd1_2213r_LNP = 0;
	Cqd1_2213i_LNP = 0;
	Cqd1_2222r_LNP = (Cqd1_00 + Cqd1_0d*YddL(1) + Cqd1_d0*YddL(1) + Cqd1_dd*YddL(1)*YddL(1) + Cqd1_u0*YucuL(1,1) + Cqd1_ud*YddL(1)*YucuL(1,1)).real();
	Cqd1_2223r_LNP = 0;
	Cqd1_2223i_LNP = 0;
	Cqd1_2233r_LNP = (Cqd1_00 + Cqd1_d0*YddL(1) + Cqd1_0d*YddL(2) + Cqd1_dd*YddL(1)*YddL(2) + Cqd1_u0*YucuL(1,1) + Cqd1_ud*YddL(2)*YucuL(1,1)).real();
	Cqd1_2311r_LNP = (Cqd1_u0*YucuL(1,2) + Cqd1_ud*YddL(0)*YucuL(1,2)).real();
	Cqd1_2311i_LNP = (Cqd1_u0*YucuL(1,2) + Cqd1_ud*YddL(0)*YucuL(1,2)).imag();
	Cqd1_2312r_LNP = 0;
	Cqd1_2312i_LNP = 0;
	Cqd1_2313r_LNP = 0;
	Cqd1_2313i_LNP = 0;
	Cqd1_2321r_LNP = 0;
	Cqd1_2321i_LNP = 0;
	Cqd1_2322r_LNP = (Cqd1_u0*YucuL(1,2) + Cqd1_ud*YddL(1)*YucuL(1,2)).real();
	Cqd1_2322i_LNP = (Cqd1_u0*YucuL(1,2) + Cqd1_ud*YddL(1)*YucuL(1,2)).imag();
	Cqd1_2323r_LNP = 0;
	Cqd1_2323i_LNP = 0;
	Cqd1_2331r_LNP = 0;
	Cqd1_2331i_LNP = 0;
	Cqd1_2332r_LNP = 0;
	Cqd1_2332i_LNP = 0;
	Cqd1_2333r_LNP = (Cqd1_u0*YucuL(1,2) + Cqd1_ud*YddL(2)*YucuL(1,2)).real();
	Cqd1_2333i_LNP = (Cqd1_u0*YucuL(1,2) + Cqd1_ud*YddL(2)*YucuL(1,2)).imag();
	Cqd1_3311r_LNP = (Cqd1_00 + Cqd1_0d*YddL(0) + Cqd1_d0*YddL(2) + Cqd1_dd*YddL(0)*YddL(2) + Cqd1_u0*YucuL(2,2) + Cqd1_ud*YddL(0)*YucuL(2,2)).real();
	Cqd1_3312r_LNP = 0;
	Cqd1_3312i_LNP = 0;
	Cqd1_3313r_LNP = 0;
	Cqd1_3313i_LNP = 0;
	Cqd1_3322r_LNP = (Cqd1_00 + Cqd1_0d*YddL(1) + Cqd1_d0*YddL(2) + Cqd1_dd*YddL(1)*YddL(2) + Cqd1_u0*YucuL(2,2) + Cqd1_ud*YddL(1)*YucuL(2,2)).real();
	Cqd1_3323r_LNP = 0;
	Cqd1_3323i_LNP = 0;
	Cqd1_3333r_LNP = (Cqd1_00 + Cqd1_0d*YddL(2) + Cqd1_d0*YddL(2) + Cqd1_dd*YddL(2)*YddL(2) + Cqd1_u0*YucuL(2,2) + Cqd1_ud*YddL(2)*YucuL(2,2)).real();

	Cqd8_1111r_LNP = (Cqd8_00 + Cqd8_0d*YddL(0) + Cqd8_d0*YddL(0) + Cqd8_dd*YddL(0)*YddL(0) + Cqd8_u0*YucuL(0,0) + Cqd8_ud*YddL(0)*YucuL(0,0)).real();
	Cqd8_1112r_LNP = 0;
	Cqd8_1112i_LNP = 0;
	Cqd8_1113r_LNP = 0;
	Cqd8_1113i_LNP = 0;
	Cqd8_1122r_LNP = (Cqd8_00 + Cqd8_d0*YddL(0) + Cqd8_0d*YddL(1) + Cqd8_dd*YddL(0)*YddL(1) + Cqd8_u0*YucuL(0,0) + Cqd8_ud*YddL(1)*YucuL(0,0)).real();
	Cqd8_1123r_LNP = 0;
	Cqd8_1123i_LNP = 0;
	Cqd8_1133r_LNP = (Cqd8_00 + Cqd8_d0*YddL(0) + Cqd8_0d*YddL(2) + Cqd8_dd*YddL(0)*YddL(2) + Cqd8_u0*YucuL(0,0) + Cqd8_ud*YddL(2)*YucuL(0,0)).real();
	Cqd8_1211r_LNP = (Cqd8_u0*YucuL(0,1) + Cqd8_ud*YddL(0)*YucuL(0,1)).real();
	Cqd8_1211i_LNP = (Cqd8_u0*YucuL(0,1) + Cqd8_ud*YddL(0)*YucuL(0,1)).imag();
	Cqd8_1212r_LNP = 0;
	Cqd8_1212i_LNP = 0;
	Cqd8_1213r_LNP = 0;
	Cqd8_1213i_LNP = 0;
	Cqd8_1221r_LNP = 0;
	Cqd8_1221i_LNP = 0;
	Cqd8_1222r_LNP = (Cqd8_u0*YucuL(0,1) + Cqd8_ud*YddL(1)*YucuL(0,1)).real();
	Cqd8_1222i_LNP = (Cqd8_u0*YucuL(0,1) + Cqd8_ud*YddL(1)*YucuL(0,1)).imag();
	Cqd8_1223r_LNP = 0;
	Cqd8_1223i_LNP = 0;
	Cqd8_1231r_LNP = 0;
	Cqd8_1231i_LNP = 0;
	Cqd8_1232r_LNP = 0;
	Cqd8_1232i_LNP = 0;
	Cqd8_1233r_LNP = (Cqd8_u0*YucuL(0,1) + Cqd8_ud*YddL(2)*YucuL(0,1)).real();
	Cqd8_1233i_LNP = (Cqd8_u0*YucuL(0,1) + Cqd8_ud*YddL(2)*YucuL(0,1)).imag();
	Cqd8_1311r_LNP = (Cqd8_u0*YucuL(0,2) + Cqd8_ud*YddL(0)*YucuL(0,2)).real();
	Cqd8_1311i_LNP = (Cqd8_u0*YucuL(0,2) + Cqd8_ud*YddL(0)*YucuL(0,2)).imag();
	Cqd8_1312r_LNP = 0;
	Cqd8_1312i_LNP = 0;
	Cqd8_1313r_LNP = 0;
	Cqd8_1313i_LNP = 0;
	Cqd8_1321r_LNP = 0;
	Cqd8_1321i_LNP = 0;
	Cqd8_1322r_LNP = (Cqd8_u0*YucuL(0,2) + Cqd8_ud*YddL(1)*YucuL(0,2)).real();
	Cqd8_1322i_LNP = (Cqd8_u0*YucuL(0,2) + Cqd8_ud*YddL(1)*YucuL(0,2)).imag();
	Cqd8_1323r_LNP = 0;
	Cqd8_1323i_LNP = 0;
	Cqd8_1331r_LNP = 0;
	Cqd8_1331i_LNP = 0;
	Cqd8_1332r_LNP = 0;
	Cqd8_1332i_LNP = 0;
	Cqd8_1333r_LNP = (Cqd8_u0*YucuL(0,2) + Cqd8_ud*YddL(2)*YucuL(0,2)).real();
	Cqd8_1333i_LNP = (Cqd8_u0*YucuL(0,2) + Cqd8_ud*YddL(2)*YucuL(0,2)).imag();
	Cqd8_2211r_LNP = (Cqd8_00 + Cqd8_0d*YddL(0) + Cqd8_d0*YddL(1) + Cqd8_dd*YddL(0)*YddL(1) + Cqd8_u0*YucuL(1,1) + Cqd8_ud*YddL(0)*YucuL(1,1)).real();
	Cqd8_2212r_LNP = 0;
	Cqd8_2212i_LNP = 0;
	Cqd8_2213r_LNP = 0;
	Cqd8_2213i_LNP = 0;
	Cqd8_2222r_LNP = (Cqd8_00 + Cqd8_0d*YddL(1) + Cqd8_d0*YddL(1) + Cqd8_dd*YddL(1)*YddL(1) + Cqd8_u0*YucuL(1,1) + Cqd8_ud*YddL(1)*YucuL(1,1)).real();
	Cqd8_2223r_LNP = 0;
	Cqd8_2223i_LNP = 0;
	Cqd8_2233r_LNP = (Cqd8_00 + Cqd8_d0*YddL(1) + Cqd8_0d*YddL(2) + Cqd8_dd*YddL(1)*YddL(2) + Cqd8_u0*YucuL(1,1) + Cqd8_ud*YddL(2)*YucuL(1,1)).real();
	Cqd8_2311r_LNP = (Cqd8_u0*YucuL(1,2) + Cqd8_ud*YddL(0)*YucuL(1,2)).real();
	Cqd8_2311i_LNP = (Cqd8_u0*YucuL(1,2) + Cqd8_ud*YddL(0)*YucuL(1,2)).imag();
	Cqd8_2312r_LNP = 0;
	Cqd8_2312i_LNP = 0;
	Cqd8_2313r_LNP = 0;
	Cqd8_2313i_LNP = 0;
	Cqd8_2321r_LNP = 0;
	Cqd8_2321i_LNP = 0;
	Cqd8_2322r_LNP = (Cqd8_u0*YucuL(1,2) + Cqd8_ud*YddL(1)*YucuL(1,2)).real();
	Cqd8_2322i_LNP = (Cqd8_u0*YucuL(1,2) + Cqd8_ud*YddL(1)*YucuL(1,2)).imag();
	Cqd8_2323r_LNP = 0;
	Cqd8_2323i_LNP = 0;
	Cqd8_2331r_LNP = 0;
	Cqd8_2331i_LNP = 0;
	Cqd8_2332r_LNP = 0;
	Cqd8_2332i_LNP = 0;
	Cqd8_2333r_LNP = (Cqd8_u0*YucuL(1,2) + Cqd8_ud*YddL(2)*YucuL(1,2)).real();
	Cqd8_2333i_LNP = (Cqd8_u0*YucuL(1,2) + Cqd8_ud*YddL(2)*YucuL(1,2)).imag();
	Cqd8_3311r_LNP = (Cqd8_00 + Cqd8_0d*YddL(0) + Cqd8_d0*YddL(2) + Cqd8_dd*YddL(0)*YddL(2) + Cqd8_u0*YucuL(2,2) + Cqd8_ud*YddL(0)*YucuL(2,2)).real();
	Cqd8_3312r_LNP = 0;
	Cqd8_3312i_LNP = 0;
	Cqd8_3313r_LNP = 0;
	Cqd8_3313i_LNP = 0;
	Cqd8_3322r_LNP = (Cqd8_00 + Cqd8_0d*YddL(1) + Cqd8_d0*YddL(2) + Cqd8_dd*YddL(1)*YddL(2) + Cqd8_u0*YucuL(2,2) + Cqd8_ud*YddL(1)*YucuL(2,2)).real();
	Cqd8_3323r_LNP = 0;
	Cqd8_3323i_LNP = 0;
	Cqd8_3333r_LNP = (Cqd8_00 + Cqd8_0d*YddL(2) + Cqd8_d0*YddL(2) + Cqd8_dd*YddL(2)*YddL(2) + Cqd8_u0*YucuL(2,2) + Cqd8_ud*YddL(2)*YucuL(2,2)).real();

	Cquqd1_1111r_LNP = (2*YdL(0)*YuL(0,0)).real();
	Cquqd1_1111i_LNP = (2*YdL(0)*YuL(0,0)).imag();
	Cquqd1_1112r_LNP = 0;
	Cquqd1_1112i_LNP = 0;
	Cquqd1_1113r_LNP = 0;
	Cquqd1_1113i_LNP = 0;
	Cquqd1_1121r_LNP = (YdL(0)*YuL(1,0)).real();
	Cquqd1_1121i_LNP = (YdL(0)*YuL(1,0)).imag();
	Cquqd1_1122r_LNP = (YdL(1)*YuL(0,0)).real();
	Cquqd1_1122i_LNP = (YdL(1)*YuL(0,0)).imag();
	Cquqd1_1123r_LNP = 0;
	Cquqd1_1123i_LNP = 0;
	Cquqd1_1131r_LNP = (YdL(0)*YuL(2,0)).real();
	Cquqd1_1131i_LNP = (YdL(0)*YuL(2,0)).imag();
	Cquqd1_1132r_LNP = 0;
	Cquqd1_1132i_LNP = 0;
	Cquqd1_1133r_LNP = (YdL(2)*YuL(0,0)).real();
	Cquqd1_1133i_LNP = (YdL(2)*YuL(0,0)).imag();
	Cquqd1_1211r_LNP = (2*YdL(0)*YuL(0,1)).real();
	Cquqd1_1211i_LNP = (2*YdL(0)*YuL(0,1)).imag();
	Cquqd1_1212r_LNP = 0;
	Cquqd1_1212i_LNP = 0;
	Cquqd1_1213r_LNP = 0;
	Cquqd1_1213i_LNP = 0;
	Cquqd1_1221r_LNP = (YdL(0)*YuL(1,1)).real();
	Cquqd1_1221i_LNP = (YdL(0)*YuL(1,1)).imag();
	Cquqd1_1222r_LNP = (YdL(1)*YuL(0,1)).real();
	Cquqd1_1222i_LNP = (YdL(1)*YuL(0,1)).imag();
	Cquqd1_1223r_LNP = 0;
	Cquqd1_1223i_LNP = 0;
	Cquqd1_1231r_LNP = (YdL(0)*YuL(2,1)).real();
	Cquqd1_1231i_LNP = (YdL(0)*YuL(2,1)).imag();
	Cquqd1_1232r_LNP = 0;
	Cquqd1_1232i_LNP = 0;
	Cquqd1_1233r_LNP = (YdL(2)*YuL(0,1)).real();
	Cquqd1_1233i_LNP = (YdL(2)*YuL(0,1)).imag();
	Cquqd1_1311r_LNP = (2*YdL(0)*YuL(0,2)).real();
	Cquqd1_1311i_LNP = (2*YdL(0)*YuL(0,2)).imag();
	Cquqd1_1312r_LNP = 0;
	Cquqd1_1312i_LNP = 0;
	Cquqd1_1313r_LNP = 0;
	Cquqd1_1313i_LNP = 0;
	Cquqd1_1321r_LNP = (YdL(0)*YuL(1,2)).real();
	Cquqd1_1321i_LNP = (YdL(0)*YuL(1,2)).imag();
	Cquqd1_1322r_LNP = (YdL(1)*YuL(0,2)).real();
	Cquqd1_1322i_LNP = (YdL(1)*YuL(0,2)).imag();
	Cquqd1_1323r_LNP = 0;
	Cquqd1_1323i_LNP = 0;
	Cquqd1_1331r_LNP = (YdL(0)*YuL(2,2)).real();
	Cquqd1_1331i_LNP = (YdL(0)*YuL(2,2)).imag();
	Cquqd1_1332r_LNP = 0;
	Cquqd1_1332i_LNP = 0;
	Cquqd1_1333r_LNP = (YdL(2)*YuL(0,2)).real();
	Cquqd1_1333i_LNP = (YdL(2)*YuL(0,2)).imag();
	Cquqd1_2111r_LNP = (YdL(0)*YuL(1,0)).real();
	Cquqd1_2111i_LNP = (YdL(0)*YuL(1,0)).imag();
	Cquqd1_2112r_LNP = (YdL(1)*YuL(0,0)).real();
	Cquqd1_2112i_LNP = (YdL(1)*YuL(0,0)).imag();
	Cquqd1_2113r_LNP = 0;
	Cquqd1_2113i_LNP = 0;
	Cquqd1_2121r_LNP = 0;
	Cquqd1_2121i_LNP = 0;
	Cquqd1_2122r_LNP = (2*YdL(1)*YuL(1,0)).real();
	Cquqd1_2122i_LNP = (2*YdL(1)*YuL(1,0)).imag();
	Cquqd1_2123r_LNP = 0;
	Cquqd1_2123i_LNP = 0;
	Cquqd1_2131r_LNP = 0;
	Cquqd1_2131i_LNP = 0;
	Cquqd1_2132r_LNP = (YdL(1)*YuL(2,0)).real();
	Cquqd1_2132i_LNP = (YdL(1)*YuL(2,0)).imag();
	Cquqd1_2133r_LNP = (YdL(2)*YuL(1,0)).real();
	Cquqd1_2133i_LNP = (YdL(2)*YuL(1,0)).imag();
	Cquqd1_2211r_LNP = (YdL(0)*YuL(1,1)).real();
	Cquqd1_2211i_LNP = (YdL(0)*YuL(1,1)).imag();
	Cquqd1_2212r_LNP = (YdL(1)*YuL(0,1)).real();
	Cquqd1_2212i_LNP = (YdL(1)*YuL(0,1)).imag();
	Cquqd1_2213r_LNP = 0;
	Cquqd1_2213i_LNP = 0;
	Cquqd1_2221r_LNP = 0;
	Cquqd1_2221i_LNP = 0;
	Cquqd1_2222r_LNP = (2*YdL(1)*YuL(1,1)).real();
	Cquqd1_2222i_LNP = (2*YdL(1)*YuL(1,1)).imag();
	Cquqd1_2223r_LNP = 0;
	Cquqd1_2223i_LNP = 0;
	Cquqd1_2231r_LNP = 0;
	Cquqd1_2231i_LNP = 0;
	Cquqd1_2232r_LNP = (YdL(1)*YuL(2,1)).real();
	Cquqd1_2232i_LNP = (YdL(1)*YuL(2,1)).imag();
	Cquqd1_2233r_LNP = (YdL(2)*YuL(1,1)).real();
	Cquqd1_2233i_LNP = (YdL(2)*YuL(1,1)).imag();
	Cquqd1_2311r_LNP = (YdL(0)*YuL(1,2)).real();
	Cquqd1_2311i_LNP = (YdL(0)*YuL(1,2)).imag();
	Cquqd1_2312r_LNP = (YdL(1)*YuL(0,2)).real();
	Cquqd1_2312i_LNP = (YdL(1)*YuL(0,2)).imag();
	Cquqd1_2313r_LNP = 0;
	Cquqd1_2313i_LNP = 0;
	Cquqd1_2321r_LNP = 0;
	Cquqd1_2321i_LNP = 0;
	Cquqd1_2322r_LNP = (2*YdL(1)*YuL(1,2)).real();
	Cquqd1_2322i_LNP = (2*YdL(1)*YuL(1,2)).imag();
	Cquqd1_2323r_LNP = 0;
	Cquqd1_2323i_LNP = 0;
	Cquqd1_2331r_LNP = 0;
	Cquqd1_2331i_LNP = 0;
	Cquqd1_2332r_LNP = (YdL(1)*YuL(2,2)).real();
	Cquqd1_2332i_LNP = (YdL(1)*YuL(2,2)).imag();
	Cquqd1_2333r_LNP = (YdL(2)*YuL(1,2)).real();
	Cquqd1_2333i_LNP = (YdL(2)*YuL(1,2)).imag();
	Cquqd1_3111r_LNP = (YdL(0)*YuL(2,0)).real();
	Cquqd1_3111i_LNP = (YdL(0)*YuL(2,0)).imag();
	Cquqd1_3112r_LNP = 0;
	Cquqd1_3112i_LNP = 0;
	Cquqd1_3113r_LNP = (YdL(2)*YuL(0,0)).real();
	Cquqd1_3113i_LNP = (YdL(2)*YuL(0,0)).imag();
	Cquqd1_3121r_LNP = 0;
	Cquqd1_3121i_LNP = 0;
	Cquqd1_3122r_LNP = (YdL(1)*YuL(2,0)).real();
	Cquqd1_3122i_LNP = (YdL(1)*YuL(2,0)).imag();
	Cquqd1_3123r_LNP = (YdL(2)*YuL(1,0)).real();
	Cquqd1_3123i_LNP = (YdL(2)*YuL(1,0)).imag();
	Cquqd1_3131r_LNP = 0;
	Cquqd1_3131i_LNP = 0;
	Cquqd1_3132r_LNP = 0;
	Cquqd1_3132i_LNP = 0;
	Cquqd1_3133r_LNP = (2*YdL(2)*YuL(2,0)).real();
	Cquqd1_3133i_LNP = (2*YdL(2)*YuL(2,0)).imag();
	Cquqd1_3211r_LNP = (YdL(0)*YuL(2,1)).real();
	Cquqd1_3211i_LNP = (YdL(0)*YuL(2,1)).imag();
	Cquqd1_3212r_LNP = 0;
	Cquqd1_3212i_LNP = 0;
	Cquqd1_3213r_LNP = (YdL(2)*YuL(0,1)).real();
	Cquqd1_3213i_LNP = (YdL(2)*YuL(0,1)).imag();
	Cquqd1_3221r_LNP = 0;
	Cquqd1_3221i_LNP = 0;
	Cquqd1_3222r_LNP = (YdL(1)*YuL(2,1)).real();
	Cquqd1_3222i_LNP = (YdL(1)*YuL(2,1)).imag();
	Cquqd1_3223r_LNP = (YdL(2)*YuL(1,1)).real();
	Cquqd1_3223i_LNP = (YdL(2)*YuL(1,1)).imag();
	Cquqd1_3231r_LNP = 0;
	Cquqd1_3231i_LNP = 0;
	Cquqd1_3232r_LNP = 0;
	Cquqd1_3232i_LNP = 0;
	Cquqd1_3233r_LNP = (2*YdL(2)*YuL(2,1)).real();
	Cquqd1_3233i_LNP = (2*YdL(2)*YuL(2,1)).imag();
	Cquqd1_3311r_LNP = (YdL(0)*YuL(2,2)).real();
	Cquqd1_3311i_LNP = (YdL(0)*YuL(2,2)).imag();
	Cquqd1_3312r_LNP = 0;
	Cquqd1_3312i_LNP = 0;
	Cquqd1_3313r_LNP = (YdL(2)*YuL(0,2)).real();
	Cquqd1_3313i_LNP = (YdL(2)*YuL(0,2)).imag();
	Cquqd1_3321r_LNP = 0;
	Cquqd1_3321i_LNP = 0;
	Cquqd1_3322r_LNP = (YdL(1)*YuL(2,2)).real();
	Cquqd1_3322i_LNP = (YdL(1)*YuL(2,2)).imag();
	Cquqd1_3323r_LNP = (YdL(2)*YuL(1,2)).real();
	Cquqd1_3323i_LNP = (YdL(2)*YuL(1,2)).imag();
	Cquqd1_3331r_LNP = 0;
	Cquqd1_3331i_LNP = 0;
	Cquqd1_3332r_LNP = 0;
	Cquqd1_3332i_LNP = 0;
	Cquqd1_3333r_LNP = (2*YdL(2)*YuL(2,2)).real();
	Cquqd1_3333i_LNP = (2*YdL(2)*YuL(2,2)).imag();

	Cquqd8_1111r_LNP = (2*YdL(0)*YuL(0,0)).real();
	Cquqd8_1111i_LNP = (2*YdL(0)*YuL(0,0)).imag();
	Cquqd8_1112r_LNP = 0;
	Cquqd8_1112i_LNP = 0;
	Cquqd8_1113r_LNP = 0;
	Cquqd8_1113i_LNP = 0;
	Cquqd8_1121r_LNP = (YdL(0)*YuL(1,0)).real();
	Cquqd8_1121i_LNP = (YdL(0)*YuL(1,0)).imag();
	Cquqd8_1122r_LNP = (YdL(1)*YuL(0,0)).real();
	Cquqd8_1122i_LNP = (YdL(1)*YuL(0,0)).imag();
	Cquqd8_1123r_LNP = 0;
	Cquqd8_1123i_LNP = 0;
	Cquqd8_1131r_LNP = (YdL(0)*YuL(2,0)).real();
	Cquqd8_1131i_LNP = (YdL(0)*YuL(2,0)).imag();
	Cquqd8_1132r_LNP = 0;
	Cquqd8_1132i_LNP = 0;
	Cquqd8_1133r_LNP = (YdL(2)*YuL(0,0)).real();
	Cquqd8_1133i_LNP = (YdL(2)*YuL(0,0)).imag();
	Cquqd8_1211r_LNP = (2*YdL(0)*YuL(0,1)).real();
	Cquqd8_1211i_LNP = (2*YdL(0)*YuL(0,1)).imag();
	Cquqd8_1212r_LNP = 0;
	Cquqd8_1212i_LNP = 0;
	Cquqd8_1213r_LNP = 0;
	Cquqd8_1213i_LNP = 0;
	Cquqd8_1221r_LNP = (YdL(0)*YuL(1,1)).real();
	Cquqd8_1221i_LNP = (YdL(0)*YuL(1,1)).imag();
	Cquqd8_1222r_LNP = (YdL(1)*YuL(0,1)).real();
	Cquqd8_1222i_LNP = (YdL(1)*YuL(0,1)).imag();
	Cquqd8_1223r_LNP = 0;
	Cquqd8_1223i_LNP = 0;
	Cquqd8_1231r_LNP = (YdL(0)*YuL(2,1)).real();
	Cquqd8_1231i_LNP = (YdL(0)*YuL(2,1)).imag();
	Cquqd8_1232r_LNP = 0;
	Cquqd8_1232i_LNP = 0;
	Cquqd8_1233r_LNP = (YdL(2)*YuL(0,1)).real();
	Cquqd8_1233i_LNP = (YdL(2)*YuL(0,1)).imag();
	Cquqd8_1311r_LNP = (2*YdL(0)*YuL(0,2)).real();
	Cquqd8_1311i_LNP = (2*YdL(0)*YuL(0,2)).imag();
	Cquqd8_1312r_LNP = 0;
	Cquqd8_1312i_LNP = 0;
	Cquqd8_1313r_LNP = 0;
	Cquqd8_1313i_LNP = 0;
	Cquqd8_1321r_LNP = (YdL(0)*YuL(1,2)).real();
	Cquqd8_1321i_LNP = (YdL(0)*YuL(1,2)).imag();
	Cquqd8_1322r_LNP = (YdL(1)*YuL(0,2)).real();
	Cquqd8_1322i_LNP = (YdL(1)*YuL(0,2)).imag();
	Cquqd8_1323r_LNP = 0;
	Cquqd8_1323i_LNP = 0;
	Cquqd8_1331r_LNP = (YdL(0)*YuL(2,2)).real();
	Cquqd8_1331i_LNP = (YdL(0)*YuL(2,2)).imag();
	Cquqd8_1332r_LNP = 0;
	Cquqd8_1332i_LNP = 0;
	Cquqd8_1333r_LNP = (YdL(2)*YuL(0,2)).real();
	Cquqd8_1333i_LNP = (YdL(2)*YuL(0,2)).imag();
	Cquqd8_2111r_LNP = (YdL(0)*YuL(1,0)).real();
	Cquqd8_2111i_LNP = (YdL(0)*YuL(1,0)).imag();
	Cquqd8_2112r_LNP = (YdL(1)*YuL(0,0)).real();
	Cquqd8_2112i_LNP = (YdL(1)*YuL(0,0)).imag();
	Cquqd8_2113r_LNP = 0;
	Cquqd8_2113i_LNP = 0;
	Cquqd8_2121r_LNP = 0;
	Cquqd8_2121i_LNP = 0;
	Cquqd8_2122r_LNP = (2*YdL(1)*YuL(1,0)).real();
	Cquqd8_2122i_LNP = (2*YdL(1)*YuL(1,0)).imag();
	Cquqd8_2123r_LNP = 0;
	Cquqd8_2123i_LNP = 0;
	Cquqd8_2131r_LNP = 0;
	Cquqd8_2131i_LNP = 0;
	Cquqd8_2132r_LNP = (YdL(1)*YuL(2,0)).real();
	Cquqd8_2132i_LNP = (YdL(1)*YuL(2,0)).imag();
	Cquqd8_2133r_LNP = (YdL(2)*YuL(1,0)).real();
	Cquqd8_2133i_LNP = (YdL(2)*YuL(1,0)).imag();
	Cquqd8_2211r_LNP = (YdL(0)*YuL(1,1)).real();
	Cquqd8_2211i_LNP = (YdL(0)*YuL(1,1)).imag();
	Cquqd8_2212r_LNP = (YdL(1)*YuL(0,1)).real();
	Cquqd8_2212i_LNP = (YdL(1)*YuL(0,1)).imag();
	Cquqd8_2213r_LNP = 0;
	Cquqd8_2213i_LNP = 0;
	Cquqd8_2221r_LNP = 0;
	Cquqd8_2221i_LNP = 0;
	Cquqd8_2222r_LNP = (2*YdL(1)*YuL(1,1)).real();
	Cquqd8_2222i_LNP = (2*YdL(1)*YuL(1,1)).imag();
	Cquqd8_2223r_LNP = 0;
	Cquqd8_2223i_LNP = 0;
	Cquqd8_2231r_LNP = 0;
	Cquqd8_2231i_LNP = 0;
	Cquqd8_2232r_LNP = (YdL(1)*YuL(2,1)).real();
	Cquqd8_2232i_LNP = (YdL(1)*YuL(2,1)).imag();
	Cquqd8_2233r_LNP = (YdL(2)*YuL(1,1)).real();
	Cquqd8_2233i_LNP = (YdL(2)*YuL(1,1)).imag();
	Cquqd8_2311r_LNP = (YdL(0)*YuL(1,2)).real();
	Cquqd8_2311i_LNP = (YdL(0)*YuL(1,2)).imag();
	Cquqd8_2312r_LNP = (YdL(1)*YuL(0,2)).real();
	Cquqd8_2312i_LNP = (YdL(1)*YuL(0,2)).imag();
	Cquqd8_2313r_LNP = 0;
	Cquqd8_2313i_LNP = 0;
	Cquqd8_2321r_LNP = 0;
	Cquqd8_2321i_LNP = 0;
	Cquqd8_2322r_LNP = (2*YdL(1)*YuL(1,2)).real();
	Cquqd8_2322i_LNP = (2*YdL(1)*YuL(1,2)).imag();
	Cquqd8_2323r_LNP = 0;
	Cquqd8_2323i_LNP = 0;
	Cquqd8_2331r_LNP = 0;
	Cquqd8_2331i_LNP = 0;
	Cquqd8_2332r_LNP = (YdL(1)*YuL(2,2)).real();
	Cquqd8_2332i_LNP = (YdL(1)*YuL(2,2)).imag();
	Cquqd8_2333r_LNP = (YdL(2)*YuL(1,2)).real();
	Cquqd8_2333i_LNP = (YdL(2)*YuL(1,2)).imag();
	Cquqd8_3111r_LNP = (YdL(0)*YuL(2,0)).real();
	Cquqd8_3111i_LNP = (YdL(0)*YuL(2,0)).imag();
	Cquqd8_3112r_LNP = 0;
	Cquqd8_3112i_LNP = 0;
	Cquqd8_3113r_LNP = (YdL(2)*YuL(0,0)).real();
	Cquqd8_3113i_LNP = (YdL(2)*YuL(0,0)).imag();
	Cquqd8_3121r_LNP = 0;
	Cquqd8_3121i_LNP = 0;
	Cquqd8_3122r_LNP = (YdL(1)*YuL(2,0)).real();
	Cquqd8_3122i_LNP = (YdL(1)*YuL(2,0)).imag();
	Cquqd8_3123r_LNP = (YdL(2)*YuL(1,0)).real();
	Cquqd8_3123i_LNP = (YdL(2)*YuL(1,0)).imag();
	Cquqd8_3131r_LNP = 0;
	Cquqd8_3131i_LNP = 0;
	Cquqd8_3132r_LNP = 0;
	Cquqd8_3132i_LNP = 0;
	Cquqd8_3133r_LNP = (2*YdL(2)*YuL(2,0)).real();
	Cquqd8_3133i_LNP = (2*YdL(2)*YuL(2,0)).imag();
	Cquqd8_3211r_LNP = (YdL(0)*YuL(2,1)).real();
	Cquqd8_3211i_LNP = (YdL(0)*YuL(2,1)).imag();
	Cquqd8_3212r_LNP = 0;
	Cquqd8_3212i_LNP = 0;
	Cquqd8_3213r_LNP = (YdL(2)*YuL(0,1)).real();
	Cquqd8_3213i_LNP = (YdL(2)*YuL(0,1)).imag();
	Cquqd8_3221r_LNP = 0;
	Cquqd8_3221i_LNP = 0;
	Cquqd8_3222r_LNP = (YdL(1)*YuL(2,1)).real();
	Cquqd8_3222i_LNP = (YdL(1)*YuL(2,1)).imag();
	Cquqd8_3223r_LNP = (YdL(2)*YuL(1,1)).real();
	Cquqd8_3223i_LNP = (YdL(2)*YuL(1,1)).imag();
	Cquqd8_3231r_LNP = 0;
	Cquqd8_3231i_LNP = 0;
	Cquqd8_3232r_LNP = 0;
	Cquqd8_3232i_LNP = 0;
	Cquqd8_3233r_LNP = (2*YdL(2)*YuL(2,1)).real();
	Cquqd8_3233i_LNP = (2*YdL(2)*YuL(2,1)).imag();
	Cquqd8_3311r_LNP = (YdL(0)*YuL(2,2)).real();
	Cquqd8_3311i_LNP = (YdL(0)*YuL(2,2)).imag();
	Cquqd8_3312r_LNP = 0;
	Cquqd8_3312i_LNP = 0;
	Cquqd8_3313r_LNP = (YdL(2)*YuL(0,2)).real();
	Cquqd8_3313i_LNP = (YdL(2)*YuL(0,2)).imag();
	Cquqd8_3321r_LNP = 0;
	Cquqd8_3321i_LNP = 0;
	Cquqd8_3322r_LNP = (YdL(1)*YuL(2,2)).real();
	Cquqd8_3322i_LNP = (YdL(1)*YuL(2,2)).imag();
	Cquqd8_3323r_LNP = (YdL(2)*YuL(1,2)).real();
	Cquqd8_3323i_LNP = (YdL(2)*YuL(1,2)).imag();
	Cquqd8_3331r_LNP = 0;
	Cquqd8_3331i_LNP = 0;
	Cquqd8_3332r_LNP = 0;
	Cquqd8_3332i_LNP = 0;
	Cquqd8_3333r_LNP = (2*YdL(2)*YuL(2,2)).real();
	Cquqd8_3333i_LNP = (2*YdL(2)*YuL(2,2)).imag();

}

bool NPSMEFTd6MFV::PostUpdate() {

	NPSMEFTd6General::GenerateSMInitialConditions();
	
    setNPSMEFTd6GeneralParameters();
    
    if (!NPSMEFTd6General::PostUpdate()) return (false);
    
    return (true);    
}