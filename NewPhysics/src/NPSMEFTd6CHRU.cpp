/* 
 * Copyright (C) 2012 HEPfit Collaboration
 *
 *
 * For the licensing terms see doc/COPYING.
 * 
 */


#include "NPSMEFTd6CHRU.h"


/*
TODO:

 - testare se flagPLR funziona

 - CT

 - definire osservabile eps_L, nel file di conf AsyGausObservable

*/



std::string NPSMEFTd6CHRU::NPSMEFTd6CHRUVars[NNPSMEFTd6CHRUVars] = {
    "mstar", "gstar", "eps_u", "eps_d",

    "cH", "cB", "cW", "c2B", "c2W", "c2G", "cHB", "cHW", "cga", "cg", "c3W", "c3G", 
    "cuG", "cuW", "cuB", "cdG", "cdW", "cdB", 
    "cuG_u", "cuW_u", "cuB_u", "cdG_u", "cdW_u", "cdB_u", 
    "cuG_d", "cuW_d", "cuB_d", "cdG_d", "cdW_d", "cdB_d",
    "cqD1_u", "cqD1_d", "cqD3_u", "cqD3_d", "cuD", "cdD", 
    "cHu", "cHd", "cHq1_u", "cHq1_d", "cHq3_u", "cHq3_d", 
    "cqq1_uu", "cqq1_ud", "cqq1_dd", "cqq1p_uu", "cqq1p_ud", "cqq1p_dd", 
    "cqq3_uu", "cqq3_ud", "cqq3_dd", "cqq3p_uu", "cqq3p_ud", "cqq3p_dd", 
    "cuu", "cuup", "cdd", "cddp", "cud1", "cud8", 
    "cqu1_u", "cqu1_d", "cqu1_y", "cqu8_u", "cqu8_d", "cqu8_y", 
    "cqd1_u", "cqd1_d", "cqd1_y", "cqd8_u", "cqd8_d", "cqd8_y", 
    "cquqd1", "cquqd1p", "cquqd8", "cquqd8p",

    "scH", "scB", "scW", "sc2B", "sc2W", "sc2G", "scHB", "scHW", "scga", "scg", "sc3W", "sc3G",
    "scuG", "scuW", "scuB", "scdG", "scdW", "scdB",
    "scuG_u", "scuW_u", "scuB_u", "scdG_u", "scdW_u", "scdB_u",
    "scuG_d", "scuW_d", "scuB_d", "scdG_d", "scdW_d", "scdB_d",
    "scqD1_u", "scqD1_d", "scqD3_u", "scqD3_d", "scuD", "scdD",
    "scHu", "scHd", "scHq1_u", "scHq1_d", "scHq3_u", "scHq3_d",
    "scqq1_uu", "scqq1_ud", "scqq1_dd", "scqq1p_uu", "scqq1p_ud", "scqq1p_dd",
    "scqq3_uu", "scqq3_ud", "scqq3_dd", "scqq3p_uu", "scqq3p_ud", "scqq3p_dd",
    "scuu", "scuup", "scdd", "scddp", "scud1", "scud8",
    "scqu1_u", "scqu1_d", "scqu1_y", "scqu8_u", "scqu8_d", "scqu8_y",
    "scqd1_u", "scqd1_d", "scqd1_y", "scqd8_u", "scqd8_d", "scqd8_y",
    "scquqd1", "scquqd1p","scquqd8", "scquqd8p",

    "Lambda_NP"
};

NPSMEFTd6CHRU::NPSMEFTd6CHRU(): NPSMEFTd6MFV() {
    setModelName("NPSMEFTd6CHRU");

    ModelParamMap.insert(std::make_pair("mstar",std::cref(mstar)));
    ModelParamMap.insert(std::make_pair("gstar",std::cref(gstar)));
    ModelParamMap.insert(std::make_pair("eps_u",std::cref(eps_u)));
    ModelParamMap.insert(std::make_pair("eps_d",std::cref(eps_d)));

    ModelParamMap.insert(std::make_pair("cH",std::cref(cH)));
    ModelParamMap.insert(std::make_pair("cB",std::cref(cB)));
    ModelParamMap.insert(std::make_pair("cW",std::cref(cW)));
    ModelParamMap.insert(std::make_pair("c2B",std::cref(c2B)));
    ModelParamMap.insert(std::make_pair("c2W",std::cref(c2W)));
    ModelParamMap.insert(std::make_pair("c2G",std::cref(c2G)));
    ModelParamMap.insert(std::make_pair("cHB",std::cref(cHB)));
    ModelParamMap.insert(std::make_pair("cHW",std::cref(cHW)));
    ModelParamMap.insert(std::make_pair("cga",std::cref(cga)));
    ModelParamMap.insert(std::make_pair("cg",std::cref(cg)));
    ModelParamMap.insert(std::make_pair("c3W",std::cref(c3W)));
    ModelParamMap.insert(std::make_pair("c3G",std::cref(c3G)));
    ModelParamMap.insert(std::make_pair("cuG",std::cref(cuG)));
    ModelParamMap.insert(std::make_pair("cuW",std::cref(cuW)));
    ModelParamMap.insert(std::make_pair("cuB",std::cref(cuB)));
    ModelParamMap.insert(std::make_pair("cdG",std::cref(cdG)));
    ModelParamMap.insert(std::make_pair("cdW",std::cref(cdW)));
    ModelParamMap.insert(std::make_pair("cdB",std::cref(cdB)));
    ModelParamMap.insert(std::make_pair("cuG_u",std::cref(cuG_u)));
    ModelParamMap.insert(std::make_pair("cuW_u",std::cref(cuW_u)));
    ModelParamMap.insert(std::make_pair("cuB_u",std::cref(cuB_u)));
    ModelParamMap.insert(std::make_pair("cdG_u",std::cref(cdG_u)));
    ModelParamMap.insert(std::make_pair("cdW_u",std::cref(cdW_u)));
    ModelParamMap.insert(std::make_pair("cdB_u",std::cref(cdB_u)));
    ModelParamMap.insert(std::make_pair("cuG_d",std::cref(cuG_d)));
    ModelParamMap.insert(std::make_pair("cuW_d",std::cref(cuW_d)));
    ModelParamMap.insert(std::make_pair("cuB_d",std::cref(cuB_d)));
    ModelParamMap.insert(std::make_pair("cdG_d",std::cref(cdG_d)));
    ModelParamMap.insert(std::make_pair("cdW_d",std::cref(cdW_d)));
    ModelParamMap.insert(std::make_pair("cdB_d",std::cref(cdB_d)));
    ModelParamMap.insert(std::make_pair("cqD1_u",std::cref(cqD1_u)));
    ModelParamMap.insert(std::make_pair("cqD1_d",std::cref(cqD1_d)));
    ModelParamMap.insert(std::make_pair("cqD3_u",std::cref(cqD3_u)));
    ModelParamMap.insert(std::make_pair("cqD3_d",std::cref(cqD3_d)));
    ModelParamMap.insert(std::make_pair("cuD",std::cref(cuD)));
    ModelParamMap.insert(std::make_pair("cdD",std::cref(cdD)));
    ModelParamMap.insert(std::make_pair("cHu",std::cref(cHu)));
    ModelParamMap.insert(std::make_pair("cHd",std::cref(cHd)));
    ModelParamMap.insert(std::make_pair("cHq1_u",std::cref(cHq1_u)));
    ModelParamMap.insert(std::make_pair("cHq1_d",std::cref(cHq1_d)));
    ModelParamMap.insert(std::make_pair("cHq3_u",std::cref(cHq3_u)));
    ModelParamMap.insert(std::make_pair("cHq3_d",std::cref(cHq3_d)));
    ModelParamMap.insert(std::make_pair("cqq1_uu",std::cref(cqq1_uu)));
    ModelParamMap.insert(std::make_pair("cqq1_ud",std::cref(cqq1_ud)));
    ModelParamMap.insert(std::make_pair("cqq1_dd",std::cref(cqq1_dd)));
    ModelParamMap.insert(std::make_pair("cqq1p_uu",std::cref(cqq1p_uu)));
    ModelParamMap.insert(std::make_pair("cqq1p_ud",std::cref(cqq1p_ud)));
    ModelParamMap.insert(std::make_pair("cqq1p_dd",std::cref(cqq1p_dd)));
    ModelParamMap.insert(std::make_pair("cqq3_uu",std::cref(cqq3_uu)));
    ModelParamMap.insert(std::make_pair("cqq3_ud",std::cref(cqq3_ud)));
    ModelParamMap.insert(std::make_pair("cqq3_dd",std::cref(cqq3_dd)));
    ModelParamMap.insert(std::make_pair("cqq3p_uu",std::cref(cqq3p_uu)));
    ModelParamMap.insert(std::make_pair("cqq3p_ud",std::cref(cqq3p_ud)));
    ModelParamMap.insert(std::make_pair("cqq3p_dd",std::cref(cqq3p_dd)));
    ModelParamMap.insert(std::make_pair("cuu",std::cref(cuu)));
    ModelParamMap.insert(std::make_pair("cuup",std::cref(cuup)));
    ModelParamMap.insert(std::make_pair("cdd",std::cref(cdd)));
    ModelParamMap.insert(std::make_pair("cddp",std::cref(cddp)));
    ModelParamMap.insert(std::make_pair("cud1",std::cref(cud1)));
    ModelParamMap.insert(std::make_pair("cud8",std::cref(cud8)));
    ModelParamMap.insert(std::make_pair("cqu1_u",std::cref(cqu1_u)));
    ModelParamMap.insert(std::make_pair("cqu1_d",std::cref(cqu1_d)));
    ModelParamMap.insert(std::make_pair("cqu1_y",std::cref(cqu1_y)));
    ModelParamMap.insert(std::make_pair("cqu8_u",std::cref(cqu8_u)));
    ModelParamMap.insert(std::make_pair("cqu8_d",std::cref(cqu8_d)));
    ModelParamMap.insert(std::make_pair("cqu8_y",std::cref(cqu8_y)));
    ModelParamMap.insert(std::make_pair("cqd1_u",std::cref(cqd1_u)));
    ModelParamMap.insert(std::make_pair("cqd1_d",std::cref(cqd1_d)));
    ModelParamMap.insert(std::make_pair("cqd1_y",std::cref(cqd1_y)));
    ModelParamMap.insert(std::make_pair("cqd8_u",std::cref(cqd8_u)));
    ModelParamMap.insert(std::make_pair("cqd8_d",std::cref(cqd8_d)));
    ModelParamMap.insert(std::make_pair("cqd8_y",std::cref(cqd8_y)));
    ModelParamMap.insert(std::make_pair("cquqd1",std::cref(cquqd1)));
    ModelParamMap.insert(std::make_pair("cquqd1p",std::cref(cquqd1p)));
    ModelParamMap.insert(std::make_pair("cquqd8",std::cref(cquqd8)));
    ModelParamMap.insert(std::make_pair("cquqd8p",std::cref(cquqd8p)));

    ModelParamMap.insert(std::make_pair("scH",std::cref(scH)));
    ModelParamMap.insert(std::make_pair("scB",std::cref(scB)));
    ModelParamMap.insert(std::make_pair("scW",std::cref(scW)));
    ModelParamMap.insert(std::make_pair("sc2B",std::cref(sc2B)));
    ModelParamMap.insert(std::make_pair("sc2W",std::cref(sc2W)));
    ModelParamMap.insert(std::make_pair("sc2G",std::cref(sc2G)));
    ModelParamMap.insert(std::make_pair("scHB",std::cref(scHB)));
    ModelParamMap.insert(std::make_pair("scHW",std::cref(scHW)));
    ModelParamMap.insert(std::make_pair("scga",std::cref(scga)));
    ModelParamMap.insert(std::make_pair("scg",std::cref(scg)));
    ModelParamMap.insert(std::make_pair("sc3W",std::cref(sc3W)));
    ModelParamMap.insert(std::make_pair("sc3G",std::cref(sc3G)));
    ModelParamMap.insert(std::make_pair("scuG",std::cref(scuG)));
    ModelParamMap.insert(std::make_pair("scuW",std::cref(scuW)));
    ModelParamMap.insert(std::make_pair("scuB",std::cref(scuB)));
    ModelParamMap.insert(std::make_pair("scdG",std::cref(scdG)));
    ModelParamMap.insert(std::make_pair("scdW",std::cref(scdW)));
    ModelParamMap.insert(std::make_pair("scdB",std::cref(scdB)));
    ModelParamMap.insert(std::make_pair("scuG_u",std::cref(scuG_u)));
    ModelParamMap.insert(std::make_pair("scuW_u",std::cref(scuW_u)));
    ModelParamMap.insert(std::make_pair("scuB_u",std::cref(scuB_u)));
    ModelParamMap.insert(std::make_pair("scdG_u",std::cref(scdG_u)));
    ModelParamMap.insert(std::make_pair("scdW_u",std::cref(scdW_u)));
    ModelParamMap.insert(std::make_pair("scdB_u",std::cref(scdB_u)));
    ModelParamMap.insert(std::make_pair("scuG_d",std::cref(scuG_d)));
    ModelParamMap.insert(std::make_pair("scuW_d",std::cref(scuW_d)));
    ModelParamMap.insert(std::make_pair("scuB_d",std::cref(scuB_d)));
    ModelParamMap.insert(std::make_pair("scdG_d",std::cref(scdG_d)));
    ModelParamMap.insert(std::make_pair("scdW_d",std::cref(scdW_d)));
    ModelParamMap.insert(std::make_pair("scdB_d",std::cref(scdB_d)));
    ModelParamMap.insert(std::make_pair("scqD1_u",std::cref(scqD1_u)));
    ModelParamMap.insert(std::make_pair("scqD1_d",std::cref(scqD1_d)));
    ModelParamMap.insert(std::make_pair("scqD3_u",std::cref(scqD3_u)));
    ModelParamMap.insert(std::make_pair("scqD3_d",std::cref(scqD3_d)));
    ModelParamMap.insert(std::make_pair("scuD",std::cref(scuD)));
    ModelParamMap.insert(std::make_pair("scdD",std::cref(scdD)));
    ModelParamMap.insert(std::make_pair("scHu",std::cref(scHu)));
    ModelParamMap.insert(std::make_pair("scHd",std::cref(scHd)));
    ModelParamMap.insert(std::make_pair("scHq1_u",std::cref(scHq1_u)));
    ModelParamMap.insert(std::make_pair("scHq1_d",std::cref(scHq1_d)));
    ModelParamMap.insert(std::make_pair("scHq3_u",std::cref(scHq3_u)));
    ModelParamMap.insert(std::make_pair("scHq3_d",std::cref(scHq3_d)));
    ModelParamMap.insert(std::make_pair("scqq1_uu",std::cref(scqq1_uu)));
    ModelParamMap.insert(std::make_pair("scqq1_ud",std::cref(scqq1_ud)));
    ModelParamMap.insert(std::make_pair("scqq1_dd",std::cref(scqq1_dd)));
    ModelParamMap.insert(std::make_pair("scqq1p_uu",std::cref(scqq1p_uu)));
    ModelParamMap.insert(std::make_pair("scqq1p_ud",std::cref(scqq1p_ud)));
    ModelParamMap.insert(std::make_pair("scqq1p_dd",std::cref(scqq1p_dd)));
    ModelParamMap.insert(std::make_pair("scqq3_uu",std::cref(scqq3_uu)));
    ModelParamMap.insert(std::make_pair("scqq3_ud",std::cref(scqq3_ud)));
    ModelParamMap.insert(std::make_pair("scqq3_dd",std::cref(scqq3_dd)));
    ModelParamMap.insert(std::make_pair("scqq3p_uu",std::cref(scqq3p_uu)));
    ModelParamMap.insert(std::make_pair("scqq3p_ud",std::cref(scqq3p_ud)));
    ModelParamMap.insert(std::make_pair("scqq3p_dd",std::cref(scqq3p_dd)));
    ModelParamMap.insert(std::make_pair("scuu",std::cref(scuu)));
    ModelParamMap.insert(std::make_pair("scuup",std::cref(scuup)));
    ModelParamMap.insert(std::make_pair("scdd",std::cref(scdd)));
    ModelParamMap.insert(std::make_pair("scddp",std::cref(scddp)));
    ModelParamMap.insert(std::make_pair("scud1",std::cref(scud1)));
    ModelParamMap.insert(std::make_pair("scud8",std::cref(scud8)));
    ModelParamMap.insert(std::make_pair("scqu1_u",std::cref(scqu1_u)));
    ModelParamMap.insert(std::make_pair("scqu1_d",std::cref(scqu1_d)));
    ModelParamMap.insert(std::make_pair("scqu1_y",std::cref(scqu1_y)));
    ModelParamMap.insert(std::make_pair("scqu8_u",std::cref(scqu8_u)));
    ModelParamMap.insert(std::make_pair("scqu8_d",std::cref(scqu8_d)));
    ModelParamMap.insert(std::make_pair("scqu8_y",std::cref(scqu8_y)));
    ModelParamMap.insert(std::make_pair("scqd1_u",std::cref(scqd1_u)));
    ModelParamMap.insert(std::make_pair("scqd1_d",std::cref(scqd1_d)));
    ModelParamMap.insert(std::make_pair("scqd1_y",std::cref(scqd1_y)));
    ModelParamMap.insert(std::make_pair("scqd8_u",std::cref(scqd8_u)));
    ModelParamMap.insert(std::make_pair("scqd8_d",std::cref(scqd8_d)));
    ModelParamMap.insert(std::make_pair("scqd8_y",std::cref(scqd8_y)));
    ModelParamMap.insert(std::make_pair("scquqd1",std::cref(scquqd1)));
    ModelParamMap.insert(std::make_pair("scquqd1p",std::cref(scquqd1p)));
    ModelParamMap.insert(std::make_pair("scquqd8",std::cref(scquqd8)));
    ModelParamMap.insert(std::make_pair("scquqd8p",std::cref(scquqd8p)));

}

bool NPSMEFTd6CHRU::Init(const std::map<std::string, double>& DPars) {

    if (FlagPLR) {
        std::cout << " -- Warning! PLR flag is True. cHd and cHq1 parameters will be ignored! -- " << std::endl;
    }

    return (NPSMEFTd6MFV::Init(DPars));
}

bool NPSMEFTd6CHRU::setFlag(const std::string name, const bool value) {
    bool res = false;
    if(name.compare("PLRFlag") == 0) {
        std::cout<<"PLRFlag = "<< value <<std::endl;
        FlagPLR = value;
        res = true;
    } else {
        res = NPSMEFTd6MFV::setFlag(name,value);
    }

    return res;
}

void NPSMEFTd6CHRU::setParameter(const std::string name, const double& value) {
    if (name.compare("mstar") == 0) {
		mstar = value;
	} else if (name.compare("gstar") == 0) {
		gstar = value;
	} else if (name.compare("eps_u") == 0) {
		eps_u = value;
	} else if (name.compare("eps_d") == 0) {
        eps_d = value;
    } else if (name.compare("cH") == 0) {
		cH = value;
	} else if (name.compare("cB") == 0) {
		cB = value;
	} else if (name.compare("cW") == 0) {
		cW = value;
	} else if (name.compare("c2B") == 0) {
		c2B = value;
	} else if (name.compare("c2W") == 0) {
		c2W = value;
	} else if (name.compare("c2G") == 0) {
		c2G = value;
	} else if (name.compare("cHB") == 0) {
		cHB = value;
	} else if (name.compare("cHW") == 0) {
		cHW = value;
	} else if (name.compare("cga") == 0) {
		cga = value;
	} else if (name.compare("cg") == 0) {
		cg = value;
	} else if (name.compare("c3W") == 0) {
		c3W = value;
	} else if (name.compare("c3G") == 0) {
		c3G = value;
	} else if (name.compare("cuG") == 0) {
		cuG = value;
	} else if (name.compare("cuW") == 0) {
		cuW = value;
	} else if (name.compare("cuB") == 0) {
		cuB = value;
	} else if (name.compare("cdG") == 0) {
		cdG = value;
	} else if (name.compare("cdW") == 0) {
		cdW = value;
	} else if (name.compare("cdB") == 0) {
		cdB = value;
	} else if (name.compare("cuG_u") == 0) {
		cuG_u = value;
	} else if (name.compare("cuW_u") == 0) {
		cuW_u = value;
	} else if (name.compare("cuB_u") == 0) {
		cuB_u = value;
	} else if (name.compare("cdG_u") == 0) {
		cdG_u = value;
	} else if (name.compare("cdW_u") == 0) {
		cdW_u = value;
	} else if (name.compare("cdB_u") == 0) {
		cdB_u = value;
	} else if (name.compare("cuG_d") == 0) {
		cuG_d = value;
	} else if (name.compare("cuW_d") == 0) {
		cuW_d = value;
	} else if (name.compare("cuB_d") == 0) {
		cuB_d = value;
	} else if (name.compare("cdG_d") == 0) {
		cdG_d = value;
	} else if (name.compare("cdW_d") == 0) {
		cdW_d = value;
	} else if (name.compare("cdB_d") == 0) {
		cdB_d = value;
	} else if (name.compare("cqD1_u") == 0) {
		cqD1_u = value;
	} else if (name.compare("cqD1_d") == 0) {
		cqD1_d = value;
	} else if (name.compare("cqD3_u") == 0) {
		cqD3_u = value;
	} else if (name.compare("cqD3_d") == 0) {
		cqD3_d = value;
	} else if (name.compare("cuD") == 0) {
		cuD = value;
	} else if (name.compare("cdD") == 0) {
		cdD = value;
	} else if (name.compare("cHu") == 0) {
		cHu = value;
	} else if (name.compare("cHd") == 0) {
		cHd = value;
	} else if (name.compare("cHq1_u") == 0) {
		cHq1_u = value;
	} else if (name.compare("cHq1_d") == 0) {
		cHq1_d = value;
	} else if (name.compare("cHq3_u") == 0) {
		cHq3_u = value;
	} else if (name.compare("cHq3_d") == 0) {
		cHq3_d = value;
	} else if (name.compare("cqq1_uu") == 0) {
		cqq1_uu = value;
	} else if (name.compare("cqq1_ud") == 0) {
		cqq1_ud = value;
	} else if (name.compare("cqq1_dd") == 0) {
		cqq1_dd = value;
	} else if (name.compare("cqq1p_uu") == 0) {
		cqq1p_uu = value;
	} else if (name.compare("cqq1p_ud") == 0) {
		cqq1p_ud = value;
	} else if (name.compare("cqq1p_dd") == 0) {
		cqq1p_dd = value;
	} else if (name.compare("cqq3_uu") == 0) {
		cqq3_uu = value;
	} else if (name.compare("cqq3_ud") == 0) {
		cqq3_ud = value;
	} else if (name.compare("cqq3_dd") == 0) {
		cqq3_dd = value;
	} else if (name.compare("cqq3p_uu") == 0) {
		cqq3p_uu = value;
	} else if (name.compare("cqq3p_ud") == 0) {
		cqq3p_ud = value;
	} else if (name.compare("cqq3p_dd") == 0) {
		cqq3p_dd = value;
	} else if (name.compare("cuu") == 0) {
		cuu = value;
	} else if (name.compare("cuup") == 0) {
		cuup = value;
	} else if (name.compare("cdd") == 0) {
		cdd = value;
	} else if (name.compare("cddp") == 0) {
		cddp = value;
	} else if (name.compare("cud1") == 0) {
		cud1 = value;
	} else if (name.compare("cud8") == 0) {
		cud8 = value;
	} else if (name.compare("cqu1_u") == 0) {
		cqu1_u = value;
	} else if (name.compare("cqu1_d") == 0) {
		cqu1_d = value;
	} else if (name.compare("cqu1_y") == 0) {
		cqu1_y = value;
	} else if (name.compare("cqu8_u") == 0) {
		cqu8_u = value;
	} else if (name.compare("cqu8_d") == 0) {
		cqu8_d = value;
	} else if (name.compare("cqu8_y") == 0) {
		cqu8_y = value;
	} else if (name.compare("cqd1_u") == 0) {
		cqd1_u = value;
	} else if (name.compare("cqd1_d") == 0) {
		cqd1_d = value;
	} else if (name.compare("cqd1_y") == 0) {
		cqd1_y = value;
	} else if (name.compare("cqd8_u") == 0) {
		cqd8_u = value;
	} else if (name.compare("cqd8_d") == 0) {
		cqd8_d = value;
	} else if (name.compare("cqd8_y") == 0) {
		cqd8_y = value;
	} else if (name.compare("cquqd1") == 0) {
		cquqd1 = value;
	} else if (name.compare("cquqd1p") == 0) {
		cquqd1p = value;
	} else if (name.compare("cquqd8") == 0) {
		cquqd8 = value;
	} else if (name.compare("cquqd8p") == 0) {
		cquqd8p = value;
	} else if (name.compare("scH") == 0) {
		scH = value;
	} else if (name.compare("scB") == 0) {
		scB = value;
	} else if (name.compare("scW") == 0) {
		scW = value;
	} else if (name.compare("sc2B") == 0) {
		sc2B = value;
	} else if (name.compare("sc2W") == 0) {
		sc2W = value;
	} else if (name.compare("sc2G") == 0) {
		sc2G = value;
	} else if (name.compare("scHB") == 0) {
		scHB = value;
	} else if (name.compare("scHW") == 0) {
		scHW = value;
	} else if (name.compare("scga") == 0) {
		scga = value;
	} else if (name.compare("scg") == 0) {
		scg = value;
	} else if (name.compare("sc3W") == 0) {
		sc3W = value;
	} else if (name.compare("sc3G") == 0) {
		sc3G = value;
	} else if (name.compare("scuG") == 0) {
		scuG = value;
	} else if (name.compare("scuW") == 0) {
		scuW = value;
	} else if (name.compare("scuB") == 0) {
		scuB = value;
	} else if (name.compare("scdG") == 0) {
		scdG = value;
	} else if (name.compare("scdW") == 0) {
		scdW = value;
	} else if (name.compare("scdB") == 0) {
		scdB = value;
	} else if (name.compare("scuG_u") == 0) {
		scuG_u = value;
	} else if (name.compare("scuW_u") == 0) {
		scuW_u = value;
	} else if (name.compare("scuB_u") == 0) {
		scuB_u = value;
	} else if (name.compare("scdG_u") == 0) {
		scdG_u = value;
	} else if (name.compare("scdW_u") == 0) {
		scdW_u = value;
	} else if (name.compare("scdB_u") == 0) {
		scdB_u = value;
	} else if (name.compare("scuG_d") == 0) {
		scuG_d = value;
	} else if (name.compare("scuW_d") == 0) {
		scuW_d = value;
	} else if (name.compare("scuB_d") == 0) {
		scuB_d = value;
	} else if (name.compare("scdG_d") == 0) {
		scdG_d = value;
	} else if (name.compare("scdW_d") == 0) {
		scdW_d = value;
	} else if (name.compare("scdB_d") == 0) {
		scdB_d = value;
	} else if (name.compare("scqD1_u") == 0) {
		scqD1_u = value;
	} else if (name.compare("scqD1_d") == 0) {
		scqD1_d = value;
	} else if (name.compare("scqD3_u") == 0) {
		scqD3_u = value;
	} else if (name.compare("scqD3_d") == 0) {
		scqD3_d = value;
	} else if (name.compare("scuD") == 0) {
		scuD = value;
	} else if (name.compare("scdD") == 0) {
		scdD = value;
	} else if (name.compare("scHu") == 0) {
		scHu = value;
	} else if (name.compare("scHd") == 0) {
		scHd = value;
	} else if (name.compare("scHq1_u") == 0) {
		scHq1_u = value;
	} else if (name.compare("scHq1_d") == 0) {
		scHq1_d = value;
	} else if (name.compare("scHq3_u") == 0) {
		scHq3_u = value;
	} else if (name.compare("scHq3_d") == 0) {
		scHq3_d = value;
	} else if (name.compare("scqq1_uu") == 0) {
		scqq1_uu = value;
	} else if (name.compare("scqq1_ud") == 0) {
		scqq1_ud = value;
	} else if (name.compare("scqq1_dd") == 0) {
		scqq1_dd = value;
	} else if (name.compare("scqq1p_uu") == 0) {
		scqq1p_uu = value;
	} else if (name.compare("scqq1p_ud") == 0) {
		scqq1p_ud = value;
	} else if (name.compare("scqq1p_dd") == 0) {
		scqq1p_dd = value;
	} else if (name.compare("scqq3_uu") == 0) {
		scqq3_uu = value;
	} else if (name.compare("scqq3_ud") == 0) {
		scqq3_ud = value;
	} else if (name.compare("scqq3_dd") == 0) {
		scqq3_dd = value;
	} else if (name.compare("scqq3p_uu") == 0) {
		scqq3p_uu = value;
	} else if (name.compare("scqq3p_ud") == 0) {
		scqq3p_ud = value;
	} else if (name.compare("scqq3p_dd") == 0) {
		scqq3p_dd = value;
	} else if (name.compare("scuu") == 0) {
		scuu = value;
	} else if (name.compare("scuup") == 0) {
		scuup = value;
	} else if (name.compare("scdd") == 0) {
		scdd = value;
	} else if (name.compare("scddp") == 0) {
		scddp = value;
	} else if (name.compare("scud1") == 0) {
		scud1 = value;
	} else if (name.compare("scud8") == 0) {
		scud8 = value;
	} else if (name.compare("scqu1_u") == 0) {
		scqu1_u = value;
	} else if (name.compare("scqu1_d") == 0) {
		scqu1_d = value;
	} else if (name.compare("scqu1_y") == 0) {
		scqu1_y = value;
	} else if (name.compare("scqu8_u") == 0) {
		scqu8_u = value;
	} else if (name.compare("scqu8_d") == 0) {
		scqu8_d = value;
	} else if (name.compare("scqu8_y") == 0) {
		scqu8_y = value;
	} else if (name.compare("scqd1_u") == 0) {
		scqd1_u = value;
	} else if (name.compare("scqd1_d") == 0) {
		scqd1_d = value;
	} else if (name.compare("scqd1_y") == 0) {
		scqd1_y = value;
	} else if (name.compare("scqd8_u") == 0) {
		scqd8_u = value;
	} else if (name.compare("scqd8_d") == 0) {
		scqd8_d = value;
	} else if (name.compare("scqd8_y") == 0) {
		scqd8_y = value;
	} else if (name.compare("scquqd1") == 0) {
		scquqd1 = value;
	} else if (name.compare("scquqd1p") == 0) {
		scquqd1p = value;
	} else if (name.compare("scquqd8") == 0) {
		scquqd8 = value;
	} else if (name.compare("scquqd8p") == 0) {
		scquqd8p = value;
	} else {
        NPSMEFTd6General::setParameter(name, value);
	}
}

void NPSMEFTd6CHRU::setNPSMEFTd6GeneralParameters() {
	double g1, g2, g3;
    double g1_2, g1_3, g1_4;
    double g2_2, g2_3, g2_4;
    double g3_2, g3_3, g3_4;

    double mstar_2;
    double gstar_2;
    double eps_u2, eps_u4;
    double eps_d2, eps_d4;

    double loop;

    double CH, CT, CB, CW, C2B, C2W, C2G, Cg, CHB, CHW, Cga, C3W, C3G;
    double CqD1_u, CqD1_d, CqD3_u, CqD3_d, CuD, CdD;

    loop   = 1/(16*M_PI*M_PI);

    mstar_2 = mstar*mstar;
    gstar_2 = gstar*gstar;

    eps_u2 = eps_u*eps_u;
    eps_u4 = eps_u*eps_u*eps_u*eps_u;
    eps_d2 = eps_d*eps_d;
    eps_d4 = eps_d*eps_d*eps_d*eps_d;

    g1 = getSMEFTCoeffEW("g1");
    g2 = getSMEFTCoeffEW("g2");
    g3 = getSMEFTCoeffEW("g3");

    g1_2 = g1*g1;
    g1_3 = g1*g1*g1;
    g1_4 = g1*g1*g1*g1;

    g2_2 = g2*g2;
    g2_3 = g2*g2*g2;
    g2_4 = g2*g2*g2*g2;

    g3_2 = g3*g3;
    g3_3 = g3*g3*g3;
    g3_4 = g3*g3*g3*g3;





    // SILH universal operators

    CH  = cH*gstar_2/(mstar_2);
    CT  = 0; // TODO: CUSTODIAL BREAKING TERM, DEPENDS ON THE REPRESENTATIONS AND FLAVOR

    CB  = sign(scB)*cB/(mstar_2);
    CW  = sign(scW)*cW/(mstar_2);

    C2B = sign(sc2B)*c2B/(gstar_2*mstar_2);
    C2W = sign(sc2W)*c2W/(gstar_2*mstar_2);
    C2G = sign(sc2G)*c2G/(gstar_2*mstar_2);

    CHB = sign(scHB)*cHB*loop/(mstar_2);
    CHW = sign(scHW)*cHW*loop/(mstar_2);

    Cga = sign(scga)*cga*loop/(mstar_2); // TODO: add yt^2
    Cg  = sign(scg)*cg*loop/(mstar_2);

    C3W = sign(sc3W)*c3W*loop/(mstar_2);
    C3G = sign(sc3G)*c3G*loop/(mstar_2);

    // Subleading EW vertices operators

    CqD1_u = sign(scqD1_u)*cqD1_u/(eps_u2*gstar_2*mstar_2);
    CqD1_d = sign(scqD1_d)*cqD1_d/(eps_d2*gstar_2*mstar_2);

    CqD3_u = sign(scqD3_u)*cqD3_u/(eps_u2*gstar_2*mstar_2);
    CqD3_d = sign(scqD3_d)*cqD3_d/(eps_d2*gstar_2*mstar_2);

    CuD    = sign(scuD)*cuD*eps_u2/mstar_2;
    CdD    = sign(scdD)*cdD*eps_d2/mstar_2;




    // SMEFT operators

    // Bosonic operators (9)

    CG_LNP		= C3G*g3_3;
    CW_LNP		= C3W*g2_3;
    CHG_LNP		= Cg*g3_2;
    CHW_LNP		= -CHW*g2_2/4;
    CHB_LNP		= -(CHB-4*Cga)*g1_2/4;
    CHWB_LNP	= -(CHB+CHW)*g1*g2/4;
    CHD_LNP		= -2*CT+CB*g1_2+CHB*g1_2-C2B*g1_4/2;
    CHbox_LNP	= (-4*CH-4*CT+2*CB*g1_2+2*CHB*g1_2-C2B*g1_4+6*CHW*g2_2+6*CW*g2_2-3*C2W*g2_4)/8;
    CH_LNP		= 0; // not 0 in principle


    // Yukawa corrections (7 = 1 leptons + 6 quarks)

    CuH_0_LNP	= (2*CHW+2*CW-C2W*g2_2)*g2_2/4;
    CdH_0_LNP	= (2*CHW+2*CW-C2W*g2_2)*g2_2/4;
    CeH_0_LNP	= (2*CHW+2*CW-C2W*g2_2)*g2_2/4;

    /*
    cuH_u_LNP	= 0;
    cuH_d_LNP	= 0;

    cdH_d_LNP	= 0;
    cuH_u_LNP	= 0;
    */

    // Dipoles (20 = 2 leptons + 18 quarks)
    // Generated at 1-loop, skip for now


    /*
    CeW_0_LNP = 0;
    CeB_0_LNP = 0;
    */

    CuG_0_LNP = sign(scuG)*cuG*loop*g3*gstar_2/(mstar_2);		
    CuW_0_LNP = sign(scuW)*cuW*loop*g2*gstar_2/(mstar_2);
    CuB_0_LNP = sign(scuB)*cuB*loop*g1*gstar_2/(mstar_2);

    CdG_0_LNP = sign(scdG)*cdG*loop*g3*gstar_2/(mstar_2);		
    CdW_0_LNP = sign(scdW)*cdW*loop*g2*gstar_2/(mstar_2);		
    CdB_0_LNP = sign(scdB)*cdB*loop*g1*gstar_2/(mstar_2);		

    CuG_u_LNP = sign(scuG_u)*cuG_u*loop*g3/(mstar_2*eps_u2);
    CuW_u_LNP = sign(scuW_u)*cuW_u*loop*g2/(mstar_2*eps_u2);
    CuB_u_LNP = sign(scuB_u)*cuB_u*loop*g1/(mstar_2*eps_u2);

    CdG_u_LNP = sign(scdG_u)*cdG_u*loop*g3/(mstar_2*eps_u2);
    CdW_u_LNP = sign(scdW_u)*cdW_u*loop*g2/(mstar_2*eps_u2);
    CdB_u_LNP = sign(scdB_u)*cdB_u*loop*g1/(mstar_2*eps_u2);

    CuG_d_LNP = sign(scuG_d)*cuG_d*loop*g3/(mstar_2*eps_d2);	
    CuW_d_LNP = sign(scuW_d)*cuW_d*loop*g2/(mstar_2*eps_d2);
    CuB_d_LNP = sign(scuB_d)*cuB_d*loop*g1/(mstar_2*eps_d2);

    CdG_d_LNP = sign(scdG_d)*cdG_d*loop*g3/(mstar_2*eps_d2);
    CdW_d_LNP = sign(scdW_d)*cdW_d*loop*g2/(mstar_2*eps_d2);
    CdB_d_LNP = sign(scdB_d)*cdB_d*loop*g1/(mstar_2*eps_d2);

    // Higgs Currents (17 = 6 leptons + 11 quarks)

    // P_LR: either cHu or cHd can be set to 0 with P_LR.
    // P_LR: also cHq3 = -cHq1  

    CHl1_0_LNP		= (-CB-CHB+C2B*g1_2)*g1_2/4;
    CHl3_0_LNP		= (CHW+CW-C2W*g2_2)*g2_2/4;
    CHe_0_LNP		= (-CB-CHB+C2B*g1_2)*g1_2/2;

    CHq1_0_LNP		= (CB+CHB-C2B*g1_2)*g1_2/12;
    CHq3_0_LNP		= (CHW+CW-C2W*g2_2)*g2_2/4;
    CHu_0_LNP		= (CB+CHB-C2B*g1_2)*g1_2/3 + sign(scHu)*cHu*eps_u2*gstar_2/(mstar_2) + CuD*g1_2/2;

    if (FlagPLR) {
        CHd_0_LNP	= (-CB-CHB+C2B*g1_2)*g1_2/6 + CdD*g1_2/2;
    } else {
        CHd_0_LNP	= (-CB-CHB+C2B*g1_2)*g1_2/6 + sign(scHd)*cHd*eps_d2*gstar_2/(mstar_2) + CdD*g1_2/2;
    }

    if (FlagPLR) { 

        CHq1_u_LNP		= -sign(scHq3_u)*cHq3_u/(mstar_2*eps_u2) + CqD1_u*g1_2/2;
        CHq3_u_LNP		= sign(scHq3_u)*cHq3_u/(mstar_2*eps_u2) + CqD3_u*g2_2/2;

        CHq1_d_LNP		= -sign(scHq3_d)*cHq3_d/(mstar_2*eps_d2) + CqD1_d*g1_2/2; 
        CHq3_d_LNP		= sign(scHq3_d)*cHq3_d/(mstar_2*eps_d2) + CqD3_d*g2_2/2;

    } else {

        CHq1_u_LNP		= sign(scHq1_u)*cHq1_u/(mstar_2*eps_u2) + CqD1_u*g1_2/2;
        CHq3_u_LNP		= sign(scHq3_u)*cHq3_u/(mstar_2*eps_u2) + CqD3_u*g2_2/2;

        CHq1_d_LNP		= sign(scHq1_d)*cHq1_d/(mstar_2*eps_d2) + CqD1_d*g1_2/2; 
        CHq3_d_LNP		= sign(scHq3_d)*cHq3_d/(mstar_2*eps_d2) + CqD3_d*g2_2/2;
    }


    // Four Fermions LL LL (4 leptons + 8 mixed + 24 quarks = 36)

    Cll_00_LNP		= (-C2B*g1_4+C2W*g2_4)/8;
    Cllp_00_LNP		= -C2W*g2_4/4;

    /*
    Cll_l0_LNP		= 0;
    Cllp_l0_LNP		= 0;
    */



    Clq1_00_LNP		= C2B*g1_4/12;
    Clq3_00_LNP		= -C2W*g2_4/4;

    /*
    Clq1_l0_LNP		= 0;
    Clq3_l0_LNP		= 0;
    */

    Clq1_0u_LNP		= -CqD1_u*g1_2/2;
    Clq3_0u_LNP		= CqD3_u*g2_2/2;

    Clq1_0d_LNP		= -CqD1_d*g1_2/2;
    Clq3_0d_LNP		= CqD3_d*g2_2/2;



    Cqq1_00_LNP		= (-C2B*g1_4+6*C2G*g3_4)/72;
    Cqq1p_00_LNP	= -C2G*g3_4/8;
    Cqq3_00_LNP		= -C2W*g2_4/8;
    Cqq3p_00_LNP	= -C2G*g3_4/8;

    Cqq1_u0_LNP		= CqD1_u*g1_2/6/2; // Extra factors of 2 for symmetrization
    Cqq1p_u0_LNP	= 0;
    Cqq3_u0_LNP		= CqD3_u*g2_2/2/2;
    Cqq3p_u0_LNP	= 0;

    Cqq1_d0_LNP		= CqD1_d*g1_2/6/2;
    Cqq1p_d0_LNP	= 0;
    Cqq3_d0_LNP		= CqD3_d*g2_2/2/2;
    Cqq3p_d0_LNP	= 0;

    Cqq1_uu_LNP		= sign(scqq1_uu)*cqq1_uu/(gstar_2*mstar_2*eps_u4);
    Cqq1_ud_LNP		= sign(scqq1_ud)*cqq1_ud/(gstar_2*mstar_2*eps_u2*eps_d2)/2;
    Cqq1_dd_LNP		= sign(scqq1_dd)*cqq1_dd/(gstar_2*mstar_2*eps_d4);
    Cqq1p_uu_LNP	= sign(scqq1p_uu)*cqq1p_uu/(gstar_2*mstar_2*eps_u4);
    Cqq1p_ud_LNP	= sign(scqq1p_ud)*cqq1p_ud/(gstar_2*mstar_2*eps_u2*eps_d2)/2;
    Cqq1p_dd_LNP	= sign(scqq1p_dd)*cqq1p_dd/(gstar_2*mstar_2*eps_d4);
    Cqq3_uu_LNP		= sign(scqq3_uu)*cqq3_uu/(gstar_2*mstar_2*eps_u4);
    Cqq3_ud_LNP		= sign(scqq3_ud)*cqq3_ud/(gstar_2*mstar_2*eps_u2*eps_d2)/2;
    Cqq3_dd_LNP		= sign(scqq3_dd)*cqq3_dd/(gstar_2*mstar_2*eps_d4);
    Cqq3p_uu_LNP	= sign(scqq3p_uu)*cqq3p_uu/(gstar_2*mstar_2*eps_u4);
    Cqq3p_ud_LNP	= sign(scqq3p_ud)*cqq3p_ud/(gstar_2*mstar_2*eps_u2*eps_d2)/2;
    Cqq3p_dd_LNP	= sign(scqq3p_dd)*cqq3p_dd/(gstar_2*mstar_2*eps_d4);


    // Four Fermions RR RR (2 leptons + 6 mixed + 22 quarks = 30)


    Cee_00_LNP		= -C2B*g1_4/2;
    /*
    Cee_e0_LNP		= 0;
    */



    Ceu_00_LNP		= C2B*2*g1_4/3 - CdD*g1_2;
    Ced_00_LNP		= -C2B*g1_4/3 - CdD*g1_2;
    /*
    Ceu_e0_LNP		= 0;
    Ceu_0u_LNP		= 0;

    Ced_e0_LNP		= 0;
    Ced_0d_LNP		= 0;
    */



    Cuu_00_LNP		= (-8*C2B*g1_4+3*C2G*g3_4)/36 + 2*CuD*g1_2/3 + sign(cuu)*cuu*eps_u4*gstar_2/(mstar_2);
    Cuup_00_LNP		= -C2G*g3_4/4 + sign(cuup)*cuup*eps_u4*gstar_2/(mstar_2);

    Cdd_00_LNP		= (-2*C2B*g1_4+3*C2G*g3_4)/36 - CdD*g1_2/3 + sign(cdd)*cdd*eps_d4*gstar_2/(mstar_2);
    Cddp_00_LNP		= -C2G*g3_4/4 + sign(cddp)*cddp*eps_d4*gstar_2/(mstar_2);

    Cud1_00_LNP		= 2*C2B*g1_4/9 + (2*CdD-CuD)*g1_2/3 + sign(cud1)*cud1*eps_u2*eps_d2*gstar_2/(mstar_2);
    Cud8_00_LNP		= -C2G*g3_4 + sign(cud8)*cud8*eps_u2*eps_d2*gstar_2/(mstar_2);

    /*
    Cuu_u0_LNP		= 0;
    Cuup_u0_LNP		= 0;

    Cdd_d0_LNP		= 0;
    Cddp_d0_LNP		= 0;

    Cud1_u0_LNP		= 0;
    Cud1_0d_LNP		= 0;
    Cud8_u0_LNP		= 0;
    Cud8_0d_LNP		= 0;

    Cuu_uu_LNP		= 0;
    Cuup_uu_LNP		= 0;

    Cdd_dd_LNP		= 0;
    Cddp_dd_LNP		= 0;

    Cud1_ud_LNP		= 0;
    Cud1p_ud_LNP	= 0;
    Cud8_ud_LNP		= 0;
    Cud8p_ud_LNP	= 0;
    */



    // Four Fermions LL RR (4 leptons + 10 mixed + 44 quarks = 58)

    Cle_00_LNP		= -C2B*g1_4/2;

    /*
    Cle_l0_LNP		= 0;
    Cle_0e_LNP		= 0;
    Cle_y_LNP		= 0;
    */



    Cqe_00_LNP		= C2B*g1_4/6;
    Clu_00_LNP		= C2B*g1_4/3 - CuD*g1_2/2;
    Cld_00_LNP		= -C2B*g1_4/6 - CdD*g1_2/2;

    Cqe_u0_LNP		= -CqD1_u*g1_2;
    Cqe_d0_LNP		= -CqD1_d*g1_2;
    /*
    Cqe_0e_LNP		= 0;

    Clu_l0_LNP		= 0;
    Clu_0u_LNP		= 0;

    Cld_l0_LNP		= 0;
    Cld_0d_LNP		= 0;
    */




    Cqu1_00_LNP		= -C2B*g1_4/9 + CuD*g1_2/6;
    Cqu8_00_LNP		= -C2G*g3_4;
    Cqd1_00_LNP		= C2B*g1_4/18 + CdD*g1_2/6;
    Cqd8_00_LNP		= -C2G*g3_4;

    Cqu1_u0_LNP		= 2*CqD1_u/3 + sign(cqu1_u)*cqu1_u/(mstar_2);
    Cqu1_d0_LNP		= 2*CqD1_d/3 + sign(cqu1_d)*cqu1_d*eps_u2/(mstar_2*eps_d2);
    Cqu1_0u_LNP		= 0;
    Cqu1_y_LNP		= sign(cqu1_y)*cqu1_y/(mstar_2);

    Cqu8_u0_LNP		= sign(cqu8_u)*cqu8_u/(mstar_2);
    Cqu8_d0_LNP		= sign(cqu8_d)*cqu8_d*eps_u2/(mstar_2*eps_d2);
    Cqu8_0u_LNP		= 0;
    Cqu8_y_LNP		= sign(cqu8_y)*cqu8_y/(mstar_2);

    Cqd1_u0_LNP		= -CqD1_u/3 + sign(cqd1_u)*cqd1_u*eps_d2/(mstar_2*eps_u2);
    Cqd1_d0_LNP		= -CqD1_d/3 + sign(cqd1_d)*cqd1_d/(mstar_2);
    Cqd1_0d_LNP		= 0;
    Cqd1_y_LNP		= sign(cqd1_y)*cqd1_y/(mstar_2);

    Cqd8_u0_LNP		= sign(cqd8_u)*cqd8_u*eps_d2/(mstar_2*eps_u2);
    Cqd8_d0_LNP		= sign(cqd8_d)*cqd8_d/(mstar_2);
    Cqd8_0d_LNP		= 0;
    Cqd8_y_LNP		= sign(cqd8_y)*cqd8_y/(mstar_2);

    /*
    Cqu1_uu_LNP		= 0;
    Cqu1_du_LNP		= 0;
    Cqu1_uy_LNP		= 0;
    Cqu1_dy_LNP		= 0;
    Cqu1_yu_LNP		= 0;
    Cqu1_yd_LNP		= 0;

    Cqu8_uu_LNP		= 0;
    Cqu8_du_LNP		= 0;
    Cqu8_uy_LNP		= 0;
    Cqu8_dy_LNP		= 0;
    Cqu8_yu_LNP		= 0;
    Cqu8_yd_LNP		= 0;

    Cqd1_ud_LNP		= 0;
    Cqd1_dd_LNP		= 0;
    Cqd1_uy_LNP		= 0;
    Cqd1_dy_LNP		= 0;
    Cqd1_yu_LNP		= 0;
    Cqd1_yd_LNP		= 0;

    Cqd8_ud_LNP		= 0;
    Cqd8_dd_LNP		= 0;
    Cqd8_uy_LNP		= 0;
    Cqd8_dy_LNP		= 0;
    Cqd8_yu_LNP		= 0;
    Cqd8_yd_LNP		= 0;
    */


    // Four Fermions LR LR (3 mixed + 20 quarks = 23)


    /*
    Cledq_00_LNP 		= 0;
    Clequ1_00_LNP 		= 0;
    Clequ3_00_LNP 		= 0;
    */



    Cquqd1_00_LNP		= sign(cquqd1)*cquqd1/(mstar_2); 
    Cquqd1p_00_LNP		= sign(cquqd1p)*cquqd1p/(mstar_2);
    Cquqd8_00_LNP		= sign(cquqd8)*cquqd8/(mstar_2);
    Cquqd8p_00_LNP		= sign(cquqd8p)*cquqd8p/(mstar_2);

    /*
    Cquqd1_u0_LNP 		= 0;
    Cquqd1_d0_LNP 		= 0;
    Cquqd1_0u_LNP 		= 0;
    Cquqd1_0d_LNP 		= 0;

    Cquqd1p_u0_LNP 		= 0;
    Cquqd1p_d0_LNP 		= 0;
    Cquqd1p_0u_LNP 		= 0;
    Cquqd1p_0d_LNP 		= 0;

    Cquqd8_u0_LNP 		= 0;
    Cquqd8_d0_LNP 		= 0;
    Cquqd8_0u_LNP 		= 0;
    Cquqd8_0d_LNP 		= 0;

    Cquqd8p_u0_LNP 		= 0;
    Cquqd8p_d0_LNP 		= 0;
    Cquqd8p_0u_LNP 		= 0;
    Cquqd8p_0d_LNP 		= 0;
    */

}

bool NPSMEFTd6CHRU::PostUpdate() {
    
    if (!NPSMEFTd6MFV::PostUpdate()) return (false);
    
    return (true);    
}







CHRU_eps_uL::CHRU_eps_uL(const StandardModel& SM_i) : ThObservable(SM_i), myNPSMEFTd6CHRU(static_cast<const NPSMEFTd6CHRU&> (SM_i)) {}

double CHRU_eps_uL::computeThValue() {
    return myNPSMEFTd6CHRU.get_eps_uL();
}

CHRU_eps_dL::CHRU_eps_dL(const StandardModel& SM_i) : ThObservable(SM_i), myNPSMEFTd6CHRU(static_cast<const NPSMEFTd6CHRU&> (SM_i)) {}

double CHRU_eps_dL::computeThValue() {
    return myNPSMEFTd6CHRU.get_eps_dL();
}