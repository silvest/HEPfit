/* 
 * Copyright (C) 2012 HEPfit Collaboration
 *
 *
 * For the licensing terms see doc/COPYING.
 */

#include "NPSMEFTd6MFVLOnoCKM.h"

std::string NPSMEFTd6MFVLOnoCKM::NPSMEFTd6MFVLOnoCKMVars[NNPSMEFTd6MFVLOnoCKMVars] = {
    "CG_LNP", "CW_LNP", "CHG_LNP", "CHW_LNP", "CHB_LNP",
    "CHWB_LNP", "CHD_LNP", "CHbox_LNP", "CH_LNP",
    "CHl1_LNP","CHl3_LNP","CHe_LNP"."CHq1_LNP","CHq3_LNP","CHu_LNP","CHd_LNP","CHud_LNP",
    "CeH_LNP","CuH_LNP","CdH_LNP","CuG_LNP","CuW_LNP","CuB_LNP","CdG_LNP","CdW_LNP","CdB_LNP","CeW_LNP","CeB_LNP",
    "Cll_LNP","Clq1_LNP","Clq3_LNP","Cee_LNP","Ceu_LNP","Ced_LNP","Cle_LNP","Clu_LNP","Cld_LNP","Cqe_LNP",
    "Cledq_LNP", "Cqq1_LNP", "Cqq3_LNP", "Cuu_LNP", "Cdd_LNP", "Cud1_LNP", "Cud8_LNP", "Cqu1_LNP", "Cqu8_LNP", "Cqd1_LNP", "Cqd8_LNP",
    "Cquqd1_LNP", "Cquqd8_LNP", "Clequ1_LNP", "Clequ3_LNP"
};

NPSMEFTd6MFVLOnoCKM::NPSMEFTd6MFVLOnoCKM()
: NPSMEFTd6General()
{
   
    ModelParamMap.insert(std::make_pair("Mz", std::cref(Mz)));
    ModelParamMap.insert(std::make_pair("AlsMz", std::cref(AlsMz)));
    ModelParamMap.insert(std::make_pair("GF", std::cref(GF)));
    ModelParamMap.insert(std::make_pair("ale", std::cref(ale)));
    ModelParamMap.insert(std::make_pair("dAle5Mz", std::cref(dAle5Mz)));
//    ModelParamMap.insert(std::make_pair("Mw_inp", std::cref(Mw_inp)));
    ModelParamMap.insert(std::make_pair("mHl", std::cref(mHl)));
    ModelParamMap.insert(std::make_pair("delMw", std::cref(delMw)));
    ModelParamMap.insert(std::make_pair("delSin2th_l", std::cref(delSin2th_l)));
    ModelParamMap.insert(std::make_pair("delSin2th_q", std::cref(delSin2th_q)));
    ModelParamMap.insert(std::make_pair("delSin2th_b", std::cref(delSin2th_b)));
    ModelParamMap.insert(std::make_pair("delGammaZ", std::cref(delGammaZ)));
    ModelParamMap.insert(std::make_pair("delsigma0H", std::cref(delsigma0H)));
    ModelParamMap.insert(std::make_pair("delR0l", std::cref(delR0l)));
    ModelParamMap.insert(std::make_pair("delR0c", std::cref(delR0c)));
    ModelParamMap.insert(std::make_pair("delR0b", std::cref(delR0b)));
    ModelParamMap.insert(std::make_pair("mneutrino_1", std::cref(leptons[NEUTRINO_1].getMass())));
    ModelParamMap.insert(std::make_pair("mneutrino_2", std::cref(leptons[NEUTRINO_2].getMass())));
    ModelParamMap.insert(std::make_pair("mneutrino_3", std::cref(leptons[NEUTRINO_3].getMass())));
    ModelParamMap.insert(std::make_pair("melectron", std::cref(leptons[ELECTRON].getMass())));
    ModelParamMap.insert(std::make_pair("mmu", std::cref(leptons[MU].getMass())));
    ModelParamMap.insert(std::make_pair("mtau", std::cref(leptons[TAU].getMass())));
    ModelParamMap.insert(std::make_pair("lambda", std::cref(lambda)));
    ModelParamMap.insert(std::make_pair("A", std::cref(A)));
    ModelParamMap.insert(std::make_pair("rhob", std::cref(rhob)));
    ModelParamMap.insert(std::make_pair("etab", std::cref(etab)));
    ModelParamMap.insert(std::make_pair("muw", std::cref(muw)));
    
}

