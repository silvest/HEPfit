/*
 * Copyright (C) 2014 HEPfit Collaboration
 *
 *
 * For the licensing terms see doc/COPYING.
 */

#include "NPEffectiveGIMR.h"
#include <limits>

const std::string NPEffectiveGIMR::NPEffectiveGIMRVars[NNPEffectiveGIMRVars]
        = {"CG", "CW", "CHG", "CHW", "CHB", "CHWB", "CHD", "CHbox", "CH",
    "CHL1_11", "CHL1_12r", "CHL1_13r", "CHL1_22", "CHL1_23r", "CHL1_33",
    "CHL1_12i", "CHL1_13i", "CHL1_23i",
    "CHL3_11", "CHL3_12r", "CHL3_13r", "CHL3_22", "CHL3_23r", "CHL3_33",
    "CHL3_12i", "CHL3_13i", "CHL3_23i",
    "CHe_11", "CHe_12r", "CHe_13r", "CHe_22", "CHe_23r", "CHe_33",
    "CHe_12i", "CHe_13i", "CHe_23i",
    "CHQ1_11", "CHQ1_12r", "CHQ1_13r", "CHQ1_22", "CHQ1_23r", "CHQ1_33",
    "CHQ1_12i", "CHQ1_13i", "CHQ1_23i",
    "CHQ3_11", "CHQ3_12r", "CHQ3_13r", "CHQ3_22", "CHQ3_23r", "CHQ3_33",
    "CHQ3_12i", "CHQ3_13i", "CHQ3_23i",
    "CHu_11", "CHu_12r", "CHu_13r", "CHu_22", "CHu_23r", "CHu_33",
    "CHu_12i", "CHu_13i", "CHu_23i",
    "CHd_11", "CHd_12r", "CHd_13r", "CHd_22", "CHd_23r", "CHd_33",
    "CHd_12i", "CHd_13i", "CHd_23i",
    "CHud_11r", "CHud_12r", "CHud_13r", "CHud_22r", "CHud_23r", "CHud_33r",
    "CHud_11i", "CHud_12i", "CHud_13i", "CHud_22i", "CHud_23i", "CHud_33i",
    "CeH_11r", "CeH_12r", "CeH_13r", "CeH_22r", "CeH_23r", "CeH_33r",
    "CeH_11i", "CeH_12i", "CeH_13i", "CeH_22i", "CeH_23i", "CeH_33i",
    "CuH_11r", "CuH_12r", "CuH_13r", "CuH_22r", "CuH_23r", "CuH_33r",
    "CuH_11i", "CuH_12i", "CuH_13i", "CuH_22i", "CuH_23i", "CuH_33i",
    "CdH_11r", "CdH_12r", "CdH_13r", "CdH_22r", "CdH_23r", "CdH_33r",
    "CdH_11i", "CdH_12i", "CdH_13i", "CdH_22i", "CdH_23i", "CdH_33i",
    "CuG_11r", "CuG_12r", "CuG_13r", "CuG_22r", "CuG_23r", "CuG_33r",
    "CuG_11i", "CuG_12i", "CuG_13i", "CuG_22i", "CuG_23i", "CuG_33i",
    "CuW_11r", "CuW_12r", "CuW_13r", "CuW_22r", "CuW_23r", "CuW_33r",
    "CuW_11i", "CuW_12i", "CuW_13i", "CuW_22i", "CuW_23i", "CuW_33i",
    "CuB_11r", "CuB_12r", "CuB_13r", "CuB_22r", "CuB_23r", "CuB_33r",
    "CuB_11i", "CuB_12i", "CuB_13i", "CuB_22i", "CuB_23i", "CuB_33i",
    "CLL_1111","CLL_1221","CLL_2112","CLL_2211","CLL_1122",
    "CLL_1133","CLL_3311","CLL_1331","CLL_3113",
    "CLQ1_1111","CLQ1_1122","CLQ1_2211","CLQ1_1221","CLQ1_2112",
    "CLQ1_1133","CLQ1_3311","CLQ1_1331","CLQ1_3113",
    "CLQ1_1123","CLQ1_2223","CLQ1_3323",
    "CLQ1_1132","CLQ1_2232","CLQ1_3332",
    "CLQ3_1111","CLQ3_1122","CLQ3_2211","CLQ3_1221","CLQ3_2112",
    "CLQ3_1133","CLQ3_3311","CLQ3_1331","CLQ3_3113",
    "CLQ3_1123","CLQ3_2223","CLQ3_3323",
    "CLQ3_1132","CLQ3_2232","CLQ3_3332",
    "Cee_1111","Cee_1122","Cee_2211","Cee_1133","Cee_3311",
    "Ceu_1111","Ceu_1122","Ceu_2211","Ceu_1133","Ceu_3311",
    "Ced_1111","Ced_1122","Ced_2211","Ced_1133","Ced_3311",
    "Ced_1123","Ced_2223","Ced_3323",
    "Ced_1132","Ced_2232","Ced_3332",
    "CLe_1111","CLe_1122","CLe_2211","CLe_1133","CLe_3311",
    "CLu_1111","CLu_1122","CLu_2211","CLu_1133","CLu_3311",
    "CLd_1111","CLd_1122","CLd_2211","CLd_1133","CLd_3311",
    "CLd_1123","CLd_2223","CLd_3323",
    "CLd_1132","CLd_2232","CLd_3332",
    "CQe_1111","CQe_1122","CQe_2211","CQe_1133","CQe_3311",
    "CQe_2311","CQe_2322","CQe_2333",
    "CQe_3211","CQe_3222","CQe_3233",
    "Lambda_NP",
    "eVBF2_HZZ1", "eVBF2_HZZ2", "eVBF2_HZZ3", "eVBF2_HZA1", "eVBF2_HZA2", "eVBF2_HAA",
    "eVBF2_HWW1", "eVBF2_HWW2", "eVBF2_HWW3", "eVBF2_Hgg", "eVBF2_HZuL", "eVBF2_HZuR",
    "eVBF2_HZdL", "eVBF2_HZdR", "eVBF2_HWud", "eVBF2_ZuL", "eVBF2_ZuR", "eVBF2_ZdL",
    "eVBF2_ZdR", "eVBF2_Wud",
    "eVBF78_HZZ1", "eVBF78_HZZ2", "eVBF78_HZZ3", "eVBF78_HZA1", "eVBF78_HZA2", "eVBF78_HAA",
    "eVBF78_HWW1", "eVBF78_HWW2", "eVBF78_HWW3", "eVBF78_Hgg", "eVBF78_HZuL", "eVBF78_HZuR",
    "eVBF78_HZdL", "eVBF78_HZdR", "eVBF78_HWud", "eVBF78_ZuL", "eVBF78_ZuR", "eVBF78_ZdL",
    "eVBF78_ZdR", "eVBF78_Wud",
    "eWH2_HWW1", "eWH2_HWW2", "eWH2_HWW3", "eWH2_HWud", "eWH2_Wud",
    "eWH78_HWW1", "eWH78_HWW2", "eWH78_HWW3", "eWH78_HWud", "eWH78_Wud",
    "eZH2_HZZ1", "eZH2_HZZ2", "eZH2_HZZ3", "eZH2_HZA1", "eZH2_HZA2", "eZH2_HZuL", "eZH2_HZuR",
    "eZH2_HZdL", "eZH2_HZdR", "eZH2_ZuL", "eZH2_ZuR", "eZH2_ZdL", "eZH2_ZdR",
    "eZH78_HZZ1", "eZH78_HZZ2", "eZH78_HZZ3", "eZH78_HZA1", "eZH78_HZA2", "eZH78_HZuL", "eZH78_HZuR",
    "eZH78_HZdL", "eZH78_HZdR", "eZH78_ZuL", "eZH78_ZuR", "eZH78_ZdL", "eZH78_ZdR",
    "ettH2_Htt", "ettH2_Hgg",
    "ettH78_Htt", "ettH78_Hgg"};

const std::string NPEffectiveGIMR::NPEffectiveGIMRVars_LFU_QFU[NNPEffectiveGIMRVars_LFU_QFU]
        = {"CG", "CW", "CHG", "CHW", "CHB", "CHWB", "CHD", "CHbox", "CH",
    "CHL1", "CHL3", "CHe", "CHQ1", "CHQ3", "CHu", "CHd", "CHud_r", "CHud_i",
    "CeH_r", "CeH_i", "CuH_r", "CuH_i", "CdH_r", "CdH_i",
    "CuG_r", "CuG_i", "CuW_r", "CuW_i", "CuB_r", "CuB_i",
    "CLL", "CLQ1", "CLQ3",
    "Cee", "Ceu", "Ced", "CLe", "CLu", "CLd", "CQe","Lambda_NP",
    "eVBF2_HZZ1", "eVBF2_HZZ2", "eVBF2_HZZ3", "eVBF2_HZA1", "eVBF2_HZA2", "eVBF2_HAA",
    "eVBF2_HWW1", "eVBF2_HWW2", "eVBF2_HWW3", "eVBF2_Hgg", "eVBF2_HZuL", "eVBF2_HZuR",
    "eVBF2_HZdL", "eVBF2_HZdR", "eVBF2_HWud", "eVBF2_ZuL", "eVBF2_ZuR", "eVBF2_ZdL",
    "eVBF2_ZdR", "eVBF2_Wud",
    "eVBF78_HZZ1", "eVBF78_HZZ2", "eVBF78_HZZ3", "eVBF78_HZA1", "eVBF78_HZA2", "eVBF78_HAA",
    "eVBF78_HWW1", "eVBF78_HWW2", "eVBF78_HWW3", "eVBF78_Hgg", "eVBF78_HZuL", "eVBF78_HZuR",
    "eVBF78_HZdL", "eVBF78_HZdR", "eVBF78_HWud", "eVBF78_ZuL", "eVBF78_ZuR", "eVBF78_ZdL",
    "eVBF78_ZdR", "eVBF78_Wud",
    "eWH2_HWW1", "eWH2_HWW2", "eWH2_HWW3", "eWH2_HWud", "eWH2_Wud",
    "eWH78_HWW1", "eWH78_HWW2", "eWH78_HWW3", "eWH78_HWud", "eWH78_Wud",
    "eZH2_HZZ1", "eZH2_HZZ2", "eZH2_HZZ3", "eZH2_HZA1", "eZH2_HZA2", "eZH2_HZuL", "eZH2_HZuR",
    "eZH2_HZdL", "eZH2_HZdR", "eZH2_ZuL", "eZH2_ZuR", "eZH2_ZdL", "eZH2_ZdR",
    "eZH78_HZZ1", "eZH78_HZZ2", "eZH78_HZZ3", "eZH78_HZA1", "eZH78_HZA2", "eZH78_HZuL", "eZH78_HZuR",
    "eZH78_HZdL", "eZH78_HZdR", "eZH78_ZuL", "eZH78_ZuR", "eZH78_ZdL", "eZH78_ZdR",
    "ettH2_Htt", "ettH2_Hgg",
    "ettH78_Htt", "ettH78_Hgg"};

NPEffectiveGIMR::NPEffectiveGIMR(const bool FlagLeptonUniversal_in, const bool FlagQuarkUniversal_in)
: NPbase(), myLEP2GIMR(trueSM), FlagLeptonUniversal(FlagLeptonUniversal_in), FlagQuarkUniversal(FlagQuarkUniversal_in)
{
    if ((!FlagLeptonUniversal && FlagQuarkUniversal)
            || (FlagLeptonUniversal && !FlagQuarkUniversal))
        throw std::runtime_error("Invalid arguments for NPEffectiveGIMR::NPEffectiveGIMR()");

    FlagMwInput = false;
    FlagQuadraticTerms = false;
    setModelLinearized();

    ModelParamMap.insert(std::pair<std::string, boost::reference_wrapper<const double> >("CG", boost::cref(CG)));
    ModelParamMap.insert(std::pair<std::string, boost::reference_wrapper<const double> >("CW", boost::cref(CW)));
    ModelParamMap.insert(std::pair<std::string, boost::reference_wrapper<const double> >("CHG", boost::cref(CHG)));
    ModelParamMap.insert(std::pair<std::string, boost::reference_wrapper<const double> >("CHW", boost::cref(CHW)));
    ModelParamMap.insert(std::pair<std::string, boost::reference_wrapper<const double> >("CHB", boost::cref(CHB)));
    ModelParamMap.insert(std::pair<std::string, boost::reference_wrapper<const double> >("CHWB", boost::cref(CHWB)));
    ModelParamMap.insert(std::pair<std::string, boost::reference_wrapper<const double> >("CHD", boost::cref(CHD)));
    ModelParamMap.insert(std::pair<std::string, boost::reference_wrapper<const double> >("CHbox", boost::cref(CHbox)));
    ModelParamMap.insert(std::pair<std::string, boost::reference_wrapper<const double> >("CH", boost::cref(CH)));
    if (FlagLeptonUniversal) {
        ModelParamMap.insert(std::pair<std::string, boost::reference_wrapper<const double> >("CHL1", boost::cref(CHL1_11)));
        ModelParamMap.insert(std::pair<std::string, boost::reference_wrapper<const double> >("CHL3", boost::cref(CHL3_11)));
        ModelParamMap.insert(std::pair<std::string, boost::reference_wrapper<const double> >("CHe", boost::cref(CHe_11)));
        ModelParamMap.insert(std::pair<std::string, boost::reference_wrapper<const double> >("CeH_r", boost::cref(CeH_11r)));
        ModelParamMap.insert(std::pair<std::string, boost::reference_wrapper<const double> >("CeH_i", boost::cref(CeH_11i)));
        ModelParamMap.insert(std::pair<std::string, boost::reference_wrapper<const double> >("CLL", boost::cref(CLL_1221)));
        ModelParamMap.insert(std::pair<std::string, boost::reference_wrapper<const double> >("Cee", boost::cref(Cee_1111)));
        ModelParamMap.insert(std::pair<std::string, boost::reference_wrapper<const double> >("CLe", boost::cref(CLe_1111)));
    } else {
        ModelParamMap.insert(std::pair<std::string, boost::reference_wrapper<const double> >("CHL1_11", boost::cref(CHL1_11)));
        ModelParamMap.insert(std::pair<std::string, boost::reference_wrapper<const double> >("CHL1_12r", boost::cref(CHL1_12r)));
        ModelParamMap.insert(std::pair<std::string, boost::reference_wrapper<const double> >("CHL1_13r", boost::cref(CHL1_13r)));
        ModelParamMap.insert(std::pair<std::string, boost::reference_wrapper<const double> >("CHL1_22", boost::cref(CHL1_22)));
        ModelParamMap.insert(std::pair<std::string, boost::reference_wrapper<const double> >("CHL1_23r", boost::cref(CHL1_23r)));
        ModelParamMap.insert(std::pair<std::string, boost::reference_wrapper<const double> >("CHL1_33", boost::cref(CHL1_33)));
        ModelParamMap.insert(std::pair<std::string, boost::reference_wrapper<const double> >("CHL1_12i", boost::cref(CHL1_12i)));
        ModelParamMap.insert(std::pair<std::string, boost::reference_wrapper<const double> >("CHL1_13i", boost::cref(CHL1_13i)));
        ModelParamMap.insert(std::pair<std::string, boost::reference_wrapper<const double> >("CHL1_23i", boost::cref(CHL1_23i)));
        ModelParamMap.insert(std::pair<std::string, boost::reference_wrapper<const double> >("CHL3_11", boost::cref(CHL3_11)));
        ModelParamMap.insert(std::pair<std::string, boost::reference_wrapper<const double> >("CHL3_12r", boost::cref(CHL3_12r)));
        ModelParamMap.insert(std::pair<std::string, boost::reference_wrapper<const double> >("CHL3_13r", boost::cref(CHL3_13r)));
        ModelParamMap.insert(std::pair<std::string, boost::reference_wrapper<const double> >("CHL3_22", boost::cref(CHL3_22)));
        ModelParamMap.insert(std::pair<std::string, boost::reference_wrapper<const double> >("CHL3_23r", boost::cref(CHL3_23r)));
        ModelParamMap.insert(std::pair<std::string, boost::reference_wrapper<const double> >("CHL3_33", boost::cref(CHL3_33)));
        ModelParamMap.insert(std::pair<std::string, boost::reference_wrapper<const double> >("CHL3_12i", boost::cref(CHL3_12i)));
        ModelParamMap.insert(std::pair<std::string, boost::reference_wrapper<const double> >("CHL3_13i", boost::cref(CHL3_13i)));
        ModelParamMap.insert(std::pair<std::string, boost::reference_wrapper<const double> >("CHL3_23i", boost::cref(CHL3_23i)));
        ModelParamMap.insert(std::pair<std::string, boost::reference_wrapper<const double> >("CHe_11", boost::cref(CHe_11)));
        ModelParamMap.insert(std::pair<std::string, boost::reference_wrapper<const double> >("CHe_12r", boost::cref(CHe_12r)));
        ModelParamMap.insert(std::pair<std::string, boost::reference_wrapper<const double> >("CHe_13r", boost::cref(CHe_13r)));
        ModelParamMap.insert(std::pair<std::string, boost::reference_wrapper<const double> >("CHe_22", boost::cref(CHe_22)));
        ModelParamMap.insert(std::pair<std::string, boost::reference_wrapper<const double> >("CHe_23r", boost::cref(CHe_23r)));
        ModelParamMap.insert(std::pair<std::string, boost::reference_wrapper<const double> >("CHe_33", boost::cref(CHe_33)));
        ModelParamMap.insert(std::pair<std::string, boost::reference_wrapper<const double> >("CHe_12i", boost::cref(CHe_12i)));
        ModelParamMap.insert(std::pair<std::string, boost::reference_wrapper<const double> >("CHe_13i", boost::cref(CHe_13i)));
        ModelParamMap.insert(std::pair<std::string, boost::reference_wrapper<const double> >("CHe_23i", boost::cref(CHe_23i)));
        ModelParamMap.insert(std::pair<std::string, boost::reference_wrapper<const double> >("CeH_11r", boost::cref(CeH_11r)));
        ModelParamMap.insert(std::pair<std::string, boost::reference_wrapper<const double> >("CeH_12r", boost::cref(CeH_12r)));
        ModelParamMap.insert(std::pair<std::string, boost::reference_wrapper<const double> >("CeH_13r", boost::cref(CeH_13r)));
        ModelParamMap.insert(std::pair<std::string, boost::reference_wrapper<const double> >("CeH_22r", boost::cref(CeH_22r)));
        ModelParamMap.insert(std::pair<std::string, boost::reference_wrapper<const double> >("CeH_23r", boost::cref(CeH_23r)));
        ModelParamMap.insert(std::pair<std::string, boost::reference_wrapper<const double> >("CeH_33r", boost::cref(CeH_33r)));
        ModelParamMap.insert(std::pair<std::string, boost::reference_wrapper<const double> >("CeH_11i", boost::cref(CeH_11i)));
        ModelParamMap.insert(std::pair<std::string, boost::reference_wrapper<const double> >("CeH_12i", boost::cref(CeH_12i)));
        ModelParamMap.insert(std::pair<std::string, boost::reference_wrapper<const double> >("CeH_13i", boost::cref(CeH_13i)));
        ModelParamMap.insert(std::pair<std::string, boost::reference_wrapper<const double> >("CeH_22i", boost::cref(CeH_22i)));
        ModelParamMap.insert(std::pair<std::string, boost::reference_wrapper<const double> >("CeH_23i", boost::cref(CeH_23i)));
        ModelParamMap.insert(std::pair<std::string, boost::reference_wrapper<const double> >("CeH_33i", boost::cref(CeH_33i)));
        ModelParamMap.insert(std::pair<std::string, boost::reference_wrapper<const double> >("CLL_1111", boost::cref(CLL_1111)));
        ModelParamMap.insert(std::pair<std::string, boost::reference_wrapper<const double> >("CLL_1221", boost::cref(CLL_1221)));
        ModelParamMap.insert(std::pair<std::string, boost::reference_wrapper<const double> >("CLL_2112", boost::cref(CLL_2112)));
        ModelParamMap.insert(std::pair<std::string, boost::reference_wrapper<const double> >("CLL_2211", boost::cref(CLL_2211)));
        ModelParamMap.insert(std::pair<std::string, boost::reference_wrapper<const double> >("CLL_1122", boost::cref(CLL_1122)));
        ModelParamMap.insert(std::pair<std::string, boost::reference_wrapper<const double> >("CLL_1331", boost::cref(CLL_1331)));
        ModelParamMap.insert(std::pair<std::string, boost::reference_wrapper<const double> >("CLL_3113", boost::cref(CLL_3113)));
        ModelParamMap.insert(std::pair<std::string, boost::reference_wrapper<const double> >("CLL_3311", boost::cref(CLL_3311)));
        ModelParamMap.insert(std::pair<std::string, boost::reference_wrapper<const double> >("CLL_1133", boost::cref(CLL_1133)));
        ModelParamMap.insert(std::pair<std::string, boost::reference_wrapper<const double> >("Cee_1111", boost::cref(Cee_1111)));
        ModelParamMap.insert(std::pair<std::string, boost::reference_wrapper<const double> >("Cee_1122", boost::cref(Cee_1122)));
        ModelParamMap.insert(std::pair<std::string, boost::reference_wrapper<const double> >("Cee_2211", boost::cref(Cee_2211)));
        ModelParamMap.insert(std::pair<std::string, boost::reference_wrapper<const double> >("Cee_1133", boost::cref(Cee_1133)));
        ModelParamMap.insert(std::pair<std::string, boost::reference_wrapper<const double> >("Cee_3311", boost::cref(Cee_3311)));
        ModelParamMap.insert(std::pair<std::string, boost::reference_wrapper<const double> >("CLe_1111", boost::cref(CLe_1111)));
        ModelParamMap.insert(std::pair<std::string, boost::reference_wrapper<const double> >("CLe_1122", boost::cref(CLe_1122)));
        ModelParamMap.insert(std::pair<std::string, boost::reference_wrapper<const double> >("CLe_2211", boost::cref(CLe_2211)));
        ModelParamMap.insert(std::pair<std::string, boost::reference_wrapper<const double> >("CLe_1133", boost::cref(CLe_1133)));
        ModelParamMap.insert(std::pair<std::string, boost::reference_wrapper<const double> >("CLe_3311", boost::cref(CLe_3311)));
        ModelParamMap.insert(std::pair<std::string, boost::reference_wrapper<const double> >("CLQ1_1123", boost::cref(CLQ1_1123)));
        ModelParamMap.insert(std::pair<std::string, boost::reference_wrapper<const double> >("CLQ1_2223", boost::cref(CLQ1_2223)));
        ModelParamMap.insert(std::pair<std::string, boost::reference_wrapper<const double> >("CLQ1_3323", boost::cref(CLQ1_3323)));
        ModelParamMap.insert(std::pair<std::string, boost::reference_wrapper<const double> >("CLQ1_1132", boost::cref(CLQ1_1132)));
        ModelParamMap.insert(std::pair<std::string, boost::reference_wrapper<const double> >("CLQ1_2232", boost::cref(CLQ1_2232)));       
        ModelParamMap.insert(std::pair<std::string, boost::reference_wrapper<const double> >("CLQ1_3332", boost::cref(CLQ1_3332)));
        ModelParamMap.insert(std::pair<std::string, boost::reference_wrapper<const double> >("CLQ3_1123", boost::cref(CLQ3_1123)));
        ModelParamMap.insert(std::pair<std::string, boost::reference_wrapper<const double> >("CLQ3_2223", boost::cref(CLQ3_2223)));
        ModelParamMap.insert(std::pair<std::string, boost::reference_wrapper<const double> >("CLQ3_3323", boost::cref(CLQ3_3323)));
        ModelParamMap.insert(std::pair<std::string, boost::reference_wrapper<const double> >("CLQ3_1132", boost::cref(CLQ3_1132)));
        ModelParamMap.insert(std::pair<std::string, boost::reference_wrapper<const double> >("CLQ3_2232", boost::cref(CLQ3_2232)));       
        ModelParamMap.insert(std::pair<std::string, boost::reference_wrapper<const double> >("CLQ3_3332", boost::cref(CLQ3_3332)));
        ModelParamMap.insert(std::pair<std::string, boost::reference_wrapper<const double> >("Ced_1123", boost::cref(Ced_1123)));
        ModelParamMap.insert(std::pair<std::string, boost::reference_wrapper<const double> >("Ced_2223", boost::cref(Ced_2223)));
        ModelParamMap.insert(std::pair<std::string, boost::reference_wrapper<const double> >("Ced_3323", boost::cref(Ced_3323)));
        ModelParamMap.insert(std::pair<std::string, boost::reference_wrapper<const double> >("Ced_1132", boost::cref(Ced_1132)));
        ModelParamMap.insert(std::pair<std::string, boost::reference_wrapper<const double> >("Ced_2232", boost::cref(Ced_2232)));       
        ModelParamMap.insert(std::pair<std::string, boost::reference_wrapper<const double> >("Ced_3332", boost::cref(Ced_3332)));
        ModelParamMap.insert(std::pair<std::string, boost::reference_wrapper<const double> >("CLd_1123", boost::cref(CLd_1123)));
        ModelParamMap.insert(std::pair<std::string, boost::reference_wrapper<const double> >("CLd_2223", boost::cref(CLd_2223)));
        ModelParamMap.insert(std::pair<std::string, boost::reference_wrapper<const double> >("CLd_3323", boost::cref(CLd_3323)));
        ModelParamMap.insert(std::pair<std::string, boost::reference_wrapper<const double> >("CLd_1132", boost::cref(CLd_1132)));
        ModelParamMap.insert(std::pair<std::string, boost::reference_wrapper<const double> >("CLd_2232", boost::cref(CLd_2232)));       
        ModelParamMap.insert(std::pair<std::string, boost::reference_wrapper<const double> >("CLd_3332", boost::cref(CLd_3332)));
        ModelParamMap.insert(std::pair<std::string, boost::reference_wrapper<const double> >("CQe_2311", boost::cref(CQe_2311)));
        ModelParamMap.insert(std::pair<std::string, boost::reference_wrapper<const double> >("CQe_2322", boost::cref(CQe_2322)));
        ModelParamMap.insert(std::pair<std::string, boost::reference_wrapper<const double> >("CQe_2333", boost::cref(CQe_2333)));
        ModelParamMap.insert(std::pair<std::string, boost::reference_wrapper<const double> >("CQe_3211", boost::cref(CQe_3211)));
        ModelParamMap.insert(std::pair<std::string, boost::reference_wrapper<const double> >("CQe_3222", boost::cref(CQe_3222)));
        ModelParamMap.insert(std::pair<std::string, boost::reference_wrapper<const double> >("CQe_3233", boost::cref(CQe_3233)));
    }
    if (FlagQuarkUniversal) {
        ModelParamMap.insert(std::pair<std::string, boost::reference_wrapper<const double> >("CHQ1", boost::cref(CHQ1_11)));
        ModelParamMap.insert(std::pair<std::string, boost::reference_wrapper<const double> >("CHQ3", boost::cref(CHQ3_11)));
        ModelParamMap.insert(std::pair<std::string, boost::reference_wrapper<const double> >("CHu", boost::cref(CHu_11)));
        ModelParamMap.insert(std::pair<std::string, boost::reference_wrapper<const double> >("CHd", boost::cref(CHd_11)));
        ModelParamMap.insert(std::pair<std::string, boost::reference_wrapper<const double> >("CHud_r", boost::cref(CHud_11r)));
        ModelParamMap.insert(std::pair<std::string, boost::reference_wrapper<const double> >("CHud_i", boost::cref(CHud_11i)));
        ModelParamMap.insert(std::pair<std::string, boost::reference_wrapper<const double> >("CuH_r", boost::cref(CuH_11r)));
        ModelParamMap.insert(std::pair<std::string, boost::reference_wrapper<const double> >("CuH_i", boost::cref(CuH_11i)));
        ModelParamMap.insert(std::pair<std::string, boost::reference_wrapper<const double> >("CdH_r", boost::cref(CdH_11r)));
        ModelParamMap.insert(std::pair<std::string, boost::reference_wrapper<const double> >("CdH_i", boost::cref(CdH_11i)));
        ModelParamMap.insert(std::pair<std::string, boost::reference_wrapper<const double> >("CuG_r", boost::cref(CuG_11r)));
        ModelParamMap.insert(std::pair<std::string, boost::reference_wrapper<const double> >("CuG_i", boost::cref(CuG_11i)));
        ModelParamMap.insert(std::pair<std::string, boost::reference_wrapper<const double> >("CuW_r", boost::cref(CuW_11r)));
        ModelParamMap.insert(std::pair<std::string, boost::reference_wrapper<const double> >("CuW_i", boost::cref(CuW_11i)));
        ModelParamMap.insert(std::pair<std::string, boost::reference_wrapper<const double> >("CuB_r", boost::cref(CuB_11r)));
        ModelParamMap.insert(std::pair<std::string, boost::reference_wrapper<const double> >("CuB_i", boost::cref(CuB_11i)));
    } else {
        ModelParamMap.insert(std::pair<std::string, boost::reference_wrapper<const double> >("CHQ1_11", boost::cref(CHQ1_11)));
        ModelParamMap.insert(std::pair<std::string, boost::reference_wrapper<const double> >("CHQ1_12r", boost::cref(CHQ1_12r)));
        ModelParamMap.insert(std::pair<std::string, boost::reference_wrapper<const double> >("CHQ1_13r", boost::cref(CHQ1_13r)));
        ModelParamMap.insert(std::pair<std::string, boost::reference_wrapper<const double> >("CHQ1_22", boost::cref(CHQ1_22)));
        ModelParamMap.insert(std::pair<std::string, boost::reference_wrapper<const double> >("CHQ1_23r", boost::cref(CHQ1_23r)));
        ModelParamMap.insert(std::pair<std::string, boost::reference_wrapper<const double> >("CHQ1_33", boost::cref(CHQ1_33)));
        ModelParamMap.insert(std::pair<std::string, boost::reference_wrapper<const double> >("CHQ1_12i", boost::cref(CHQ1_12i)));
        ModelParamMap.insert(std::pair<std::string, boost::reference_wrapper<const double> >("CHQ1_13i", boost::cref(CHQ1_13i)));
        ModelParamMap.insert(std::pair<std::string, boost::reference_wrapper<const double> >("CHQ1_23i", boost::cref(CHQ1_23i)));
        ModelParamMap.insert(std::pair<std::string, boost::reference_wrapper<const double> >("CHQ3_11", boost::cref(CHQ3_11)));
        ModelParamMap.insert(std::pair<std::string, boost::reference_wrapper<const double> >("CHQ3_12r", boost::cref(CHQ3_12r)));
        ModelParamMap.insert(std::pair<std::string, boost::reference_wrapper<const double> >("CHQ3_13r", boost::cref(CHQ3_13r)));
        ModelParamMap.insert(std::pair<std::string, boost::reference_wrapper<const double> >("CHQ3_22", boost::cref(CHQ3_22)));
        ModelParamMap.insert(std::pair<std::string, boost::reference_wrapper<const double> >("CHQ3_23r", boost::cref(CHQ3_23r)));
        ModelParamMap.insert(std::pair<std::string, boost::reference_wrapper<const double> >("CHQ3_33", boost::cref(CHQ3_33)));
        ModelParamMap.insert(std::pair<std::string, boost::reference_wrapper<const double> >("CHQ3_12i", boost::cref(CHQ3_12i)));
        ModelParamMap.insert(std::pair<std::string, boost::reference_wrapper<const double> >("CHQ3_13i", boost::cref(CHQ3_13i)));
        ModelParamMap.insert(std::pair<std::string, boost::reference_wrapper<const double> >("CHQ3_23i", boost::cref(CHQ3_23i)));
        ModelParamMap.insert(std::pair<std::string, boost::reference_wrapper<const double> >("CHu_11", boost::cref(CHu_11)));
        ModelParamMap.insert(std::pair<std::string, boost::reference_wrapper<const double> >("CHu_12r", boost::cref(CHu_12r)));
        ModelParamMap.insert(std::pair<std::string, boost::reference_wrapper<const double> >("CHu_13r", boost::cref(CHu_13r)));
        ModelParamMap.insert(std::pair<std::string, boost::reference_wrapper<const double> >("CHu_22", boost::cref(CHu_22)));
        ModelParamMap.insert(std::pair<std::string, boost::reference_wrapper<const double> >("CHu_23r", boost::cref(CHu_23r)));
        ModelParamMap.insert(std::pair<std::string, boost::reference_wrapper<const double> >("CHu_33", boost::cref(CHu_33)));
        ModelParamMap.insert(std::pair<std::string, boost::reference_wrapper<const double> >("CHu_12i", boost::cref(CHu_12i)));
        ModelParamMap.insert(std::pair<std::string, boost::reference_wrapper<const double> >("CHu_13i", boost::cref(CHu_13i)));
        ModelParamMap.insert(std::pair<std::string, boost::reference_wrapper<const double> >("CHu_23i", boost::cref(CHu_23i)));
        ModelParamMap.insert(std::pair<std::string, boost::reference_wrapper<const double> >("CHd_11", boost::cref(CHd_11)));
        ModelParamMap.insert(std::pair<std::string, boost::reference_wrapper<const double> >("CHd_12r", boost::cref(CHd_12r)));
        ModelParamMap.insert(std::pair<std::string, boost::reference_wrapper<const double> >("CHd_13r", boost::cref(CHd_13r)));
        ModelParamMap.insert(std::pair<std::string, boost::reference_wrapper<const double> >("CHd_22", boost::cref(CHd_22)));
        ModelParamMap.insert(std::pair<std::string, boost::reference_wrapper<const double> >("CHd_23r", boost::cref(CHd_23r)));
        ModelParamMap.insert(std::pair<std::string, boost::reference_wrapper<const double> >("CHd_33", boost::cref(CHd_33)));
        ModelParamMap.insert(std::pair<std::string, boost::reference_wrapper<const double> >("CHd_12i", boost::cref(CHd_12i)));
        ModelParamMap.insert(std::pair<std::string, boost::reference_wrapper<const double> >("CHd_13i", boost::cref(CHd_13i)));
        ModelParamMap.insert(std::pair<std::string, boost::reference_wrapper<const double> >("CHd_23i", boost::cref(CHd_23i)));
        ModelParamMap.insert(std::pair<std::string, boost::reference_wrapper<const double> >("CHud_11r", boost::cref(CHud_11r)));
        ModelParamMap.insert(std::pair<std::string, boost::reference_wrapper<const double> >("CHud_12r", boost::cref(CHud_12r)));
        ModelParamMap.insert(std::pair<std::string, boost::reference_wrapper<const double> >("CHud_13r", boost::cref(CHud_13r)));
        ModelParamMap.insert(std::pair<std::string, boost::reference_wrapper<const double> >("CHud_22r", boost::cref(CHud_22r)));
        ModelParamMap.insert(std::pair<std::string, boost::reference_wrapper<const double> >("CHud_23r", boost::cref(CHud_23r)));
        ModelParamMap.insert(std::pair<std::string, boost::reference_wrapper<const double> >("CHud_33r", boost::cref(CHud_33r)));
        ModelParamMap.insert(std::pair<std::string, boost::reference_wrapper<const double> >("CHud_11i", boost::cref(CHud_11i)));
        ModelParamMap.insert(std::pair<std::string, boost::reference_wrapper<const double> >("CHud_12i", boost::cref(CHud_12i)));
        ModelParamMap.insert(std::pair<std::string, boost::reference_wrapper<const double> >("CHud_13i", boost::cref(CHud_13i)));
        ModelParamMap.insert(std::pair<std::string, boost::reference_wrapper<const double> >("CHud_22i", boost::cref(CHud_22i)));
        ModelParamMap.insert(std::pair<std::string, boost::reference_wrapper<const double> >("CHud_23i", boost::cref(CHud_23i)));
        ModelParamMap.insert(std::pair<std::string, boost::reference_wrapper<const double> >("CHud_33i", boost::cref(CHud_33i)));
        ModelParamMap.insert(std::pair<std::string, boost::reference_wrapper<const double> >("CuH_11r", boost::cref(CuH_11r)));
        ModelParamMap.insert(std::pair<std::string, boost::reference_wrapper<const double> >("CuH_12r", boost::cref(CuH_12r)));
        ModelParamMap.insert(std::pair<std::string, boost::reference_wrapper<const double> >("CuH_13r", boost::cref(CuH_13r)));
        ModelParamMap.insert(std::pair<std::string, boost::reference_wrapper<const double> >("CuH_22r", boost::cref(CuH_22r)));
        ModelParamMap.insert(std::pair<std::string, boost::reference_wrapper<const double> >("CuH_23r", boost::cref(CuH_23r)));
        ModelParamMap.insert(std::pair<std::string, boost::reference_wrapper<const double> >("CuH_33r", boost::cref(CuH_33r)));
        ModelParamMap.insert(std::pair<std::string, boost::reference_wrapper<const double> >("CuH_11i", boost::cref(CuH_11i)));
        ModelParamMap.insert(std::pair<std::string, boost::reference_wrapper<const double> >("CuH_12i", boost::cref(CuH_12i)));
        ModelParamMap.insert(std::pair<std::string, boost::reference_wrapper<const double> >("CuH_13i", boost::cref(CuH_13i)));
        ModelParamMap.insert(std::pair<std::string, boost::reference_wrapper<const double> >("CuH_22i", boost::cref(CuH_22i)));
        ModelParamMap.insert(std::pair<std::string, boost::reference_wrapper<const double> >("CuH_23i", boost::cref(CuH_23i)));
        ModelParamMap.insert(std::pair<std::string, boost::reference_wrapper<const double> >("CuH_33i", boost::cref(CuH_33i)));
        ModelParamMap.insert(std::pair<std::string, boost::reference_wrapper<const double> >("CdH_11r", boost::cref(CdH_11r)));
        ModelParamMap.insert(std::pair<std::string, boost::reference_wrapper<const double> >("CdH_12r", boost::cref(CdH_12r)));
        ModelParamMap.insert(std::pair<std::string, boost::reference_wrapper<const double> >("CdH_13r", boost::cref(CdH_13r)));
        ModelParamMap.insert(std::pair<std::string, boost::reference_wrapper<const double> >("CdH_22r", boost::cref(CdH_22r)));
        ModelParamMap.insert(std::pair<std::string, boost::reference_wrapper<const double> >("CdH_23r", boost::cref(CdH_23r)));
        ModelParamMap.insert(std::pair<std::string, boost::reference_wrapper<const double> >("CdH_33r", boost::cref(CdH_33r)));
        ModelParamMap.insert(std::pair<std::string, boost::reference_wrapper<const double> >("CdH_11i", boost::cref(CdH_11i)));
        ModelParamMap.insert(std::pair<std::string, boost::reference_wrapper<const double> >("CdH_12i", boost::cref(CdH_12i)));
        ModelParamMap.insert(std::pair<std::string, boost::reference_wrapper<const double> >("CdH_13i", boost::cref(CdH_13i)));
        ModelParamMap.insert(std::pair<std::string, boost::reference_wrapper<const double> >("CdH_22i", boost::cref(CdH_22i)));
        ModelParamMap.insert(std::pair<std::string, boost::reference_wrapper<const double> >("CdH_23i", boost::cref(CdH_23i)));
        ModelParamMap.insert(std::pair<std::string, boost::reference_wrapper<const double> >("CdH_33i", boost::cref(CdH_33i)));
        ModelParamMap.insert(std::pair<std::string, boost::reference_wrapper<const double> >("CuG_11r", boost::cref(CuG_11r)));
        ModelParamMap.insert(std::pair<std::string, boost::reference_wrapper<const double> >("CuG_12r", boost::cref(CuG_12r)));
        ModelParamMap.insert(std::pair<std::string, boost::reference_wrapper<const double> >("CuG_13r", boost::cref(CuG_13r)));
        ModelParamMap.insert(std::pair<std::string, boost::reference_wrapper<const double> >("CuG_22r", boost::cref(CuG_22r)));
        ModelParamMap.insert(std::pair<std::string, boost::reference_wrapper<const double> >("CuG_23r", boost::cref(CuG_23r)));
        ModelParamMap.insert(std::pair<std::string, boost::reference_wrapper<const double> >("CuG_33r", boost::cref(CuG_33r)));
        ModelParamMap.insert(std::pair<std::string, boost::reference_wrapper<const double> >("CuG_11i", boost::cref(CuG_11i)));
        ModelParamMap.insert(std::pair<std::string, boost::reference_wrapper<const double> >("CuG_12i", boost::cref(CuG_12i)));
        ModelParamMap.insert(std::pair<std::string, boost::reference_wrapper<const double> >("CuG_13i", boost::cref(CuG_13i)));
        ModelParamMap.insert(std::pair<std::string, boost::reference_wrapper<const double> >("CuG_22i", boost::cref(CuG_22i)));
        ModelParamMap.insert(std::pair<std::string, boost::reference_wrapper<const double> >("CuG_23i", boost::cref(CuG_23i)));
        ModelParamMap.insert(std::pair<std::string, boost::reference_wrapper<const double> >("CuG_33i", boost::cref(CuG_33i)));
        ModelParamMap.insert(std::pair<std::string, boost::reference_wrapper<const double> >("CuW_11r", boost::cref(CuW_11r)));
        ModelParamMap.insert(std::pair<std::string, boost::reference_wrapper<const double> >("CuW_12r", boost::cref(CuW_12r)));
        ModelParamMap.insert(std::pair<std::string, boost::reference_wrapper<const double> >("CuW_13r", boost::cref(CuW_13r)));
        ModelParamMap.insert(std::pair<std::string, boost::reference_wrapper<const double> >("CuW_22r", boost::cref(CuW_22r)));
        ModelParamMap.insert(std::pair<std::string, boost::reference_wrapper<const double> >("CuW_23r", boost::cref(CuW_23r)));
        ModelParamMap.insert(std::pair<std::string, boost::reference_wrapper<const double> >("CuW_33r", boost::cref(CuW_33r)));
        ModelParamMap.insert(std::pair<std::string, boost::reference_wrapper<const double> >("CuW_11i", boost::cref(CuW_11i)));
        ModelParamMap.insert(std::pair<std::string, boost::reference_wrapper<const double> >("CuW_12i", boost::cref(CuW_12i)));
        ModelParamMap.insert(std::pair<std::string, boost::reference_wrapper<const double> >("CuW_13i", boost::cref(CuW_13i)));
        ModelParamMap.insert(std::pair<std::string, boost::reference_wrapper<const double> >("CuW_22i", boost::cref(CuW_22i)));
        ModelParamMap.insert(std::pair<std::string, boost::reference_wrapper<const double> >("CuW_23i", boost::cref(CuW_23i)));
        ModelParamMap.insert(std::pair<std::string, boost::reference_wrapper<const double> >("CuW_33i", boost::cref(CuW_33i)));
        ModelParamMap.insert(std::pair<std::string, boost::reference_wrapper<const double> >("CuB_11r", boost::cref(CuB_11r)));
        ModelParamMap.insert(std::pair<std::string, boost::reference_wrapper<const double> >("CuB_12r", boost::cref(CuB_12r)));
        ModelParamMap.insert(std::pair<std::string, boost::reference_wrapper<const double> >("CuB_13r", boost::cref(CuB_13r)));
        ModelParamMap.insert(std::pair<std::string, boost::reference_wrapper<const double> >("CuB_22r", boost::cref(CuB_22r)));
        ModelParamMap.insert(std::pair<std::string, boost::reference_wrapper<const double> >("CuB_23r", boost::cref(CuB_23r)));
        ModelParamMap.insert(std::pair<std::string, boost::reference_wrapper<const double> >("CuB_33r", boost::cref(CuB_33r)));
        ModelParamMap.insert(std::pair<std::string, boost::reference_wrapper<const double> >("CuB_11i", boost::cref(CuB_11i)));
        ModelParamMap.insert(std::pair<std::string, boost::reference_wrapper<const double> >("CuB_12i", boost::cref(CuB_12i)));
        ModelParamMap.insert(std::pair<std::string, boost::reference_wrapper<const double> >("CuB_13i", boost::cref(CuB_13i)));
        ModelParamMap.insert(std::pair<std::string, boost::reference_wrapper<const double> >("CuB_22i", boost::cref(CuB_22i)));
        ModelParamMap.insert(std::pair<std::string, boost::reference_wrapper<const double> >("CuB_23i", boost::cref(CuB_23i)));
        ModelParamMap.insert(std::pair<std::string, boost::reference_wrapper<const double> >("CuB_33i", boost::cref(CuB_33i)));
    }
    if(FlagLeptonUniversal && FlagQuarkUniversal){
        ModelParamMap.insert(std::pair<std::string, boost::reference_wrapper<const double> >("CLQ1", boost::cref(CLQ1_1111)));  
        ModelParamMap.insert(std::pair<std::string, boost::reference_wrapper<const double> >("CLQ3", boost::cref(CLQ3_1111))); 
        ModelParamMap.insert(std::pair<std::string, boost::reference_wrapper<const double> >("Ceu", boost::cref(Ceu_1111)));
        ModelParamMap.insert(std::pair<std::string, boost::reference_wrapper<const double> >("Ced", boost::cref(Ced_1111)));
        ModelParamMap.insert(std::pair<std::string, boost::reference_wrapper<const double> >("CLu", boost::cref(CLu_1111)));
        ModelParamMap.insert(std::pair<std::string, boost::reference_wrapper<const double> >("CLd", boost::cref(CLd_1111)));
        ModelParamMap.insert(std::pair<std::string, boost::reference_wrapper<const double> >("CQe", boost::cref(CQe_1111)));
    } else {
        ModelParamMap.insert(std::pair<std::string, boost::reference_wrapper<const double> >("CLQ1_1111", boost::cref(CLQ1_1111)));
        ModelParamMap.insert(std::pair<std::string, boost::reference_wrapper<const double> >("CLQ1_1122", boost::cref(CLQ1_1122)));
        ModelParamMap.insert(std::pair<std::string, boost::reference_wrapper<const double> >("CLQ1_2211", boost::cref(CLQ1_2211)));
        ModelParamMap.insert(std::pair<std::string, boost::reference_wrapper<const double> >("CLQ1_1221", boost::cref(CLQ1_1221)));
        ModelParamMap.insert(std::pair<std::string, boost::reference_wrapper<const double> >("CLQ1_2112", boost::cref(CLQ1_2112)));
        ModelParamMap.insert(std::pair<std::string, boost::reference_wrapper<const double> >("CLQ1_1133", boost::cref(CLQ1_1133)));
        ModelParamMap.insert(std::pair<std::string, boost::reference_wrapper<const double> >("CLQ1_3311", boost::cref(CLQ1_3311)));
        ModelParamMap.insert(std::pair<std::string, boost::reference_wrapper<const double> >("CLQ1_1331", boost::cref(CLQ1_1331)));
        ModelParamMap.insert(std::pair<std::string, boost::reference_wrapper<const double> >("CLQ1_3113", boost::cref(CLQ1_3113)));
        ModelParamMap.insert(std::pair<std::string, boost::reference_wrapper<const double> >("CLQ3_1111", boost::cref(CLQ3_1111)));
        ModelParamMap.insert(std::pair<std::string, boost::reference_wrapper<const double> >("CLQ3_1122", boost::cref(CLQ3_1122)));
        ModelParamMap.insert(std::pair<std::string, boost::reference_wrapper<const double> >("CLQ3_2211", boost::cref(CLQ3_2211)));
        ModelParamMap.insert(std::pair<std::string, boost::reference_wrapper<const double> >("CLQ3_1221", boost::cref(CLQ3_1221)));
        ModelParamMap.insert(std::pair<std::string, boost::reference_wrapper<const double> >("CLQ3_2112", boost::cref(CLQ3_2112)));
        ModelParamMap.insert(std::pair<std::string, boost::reference_wrapper<const double> >("CLQ3_1133", boost::cref(CLQ3_1133)));
        ModelParamMap.insert(std::pair<std::string, boost::reference_wrapper<const double> >("CLQ3_3311", boost::cref(CLQ3_3311)));
        ModelParamMap.insert(std::pair<std::string, boost::reference_wrapper<const double> >("CLQ3_1331", boost::cref(CLQ3_1331)));
        ModelParamMap.insert(std::pair<std::string, boost::reference_wrapper<const double> >("CLQ3_3113", boost::cref(CLQ3_3113)));
        ModelParamMap.insert(std::pair<std::string, boost::reference_wrapper<const double> >("Ceu_1111", boost::cref(Ceu_1111)));
        ModelParamMap.insert(std::pair<std::string, boost::reference_wrapper<const double> >("Ceu_1122", boost::cref(Ceu_1122)));
        ModelParamMap.insert(std::pair<std::string, boost::reference_wrapper<const double> >("Ceu_2211", boost::cref(Ceu_2211)));
        ModelParamMap.insert(std::pair<std::string, boost::reference_wrapper<const double> >("Ceu_1133", boost::cref(Ceu_1133)));
        ModelParamMap.insert(std::pair<std::string, boost::reference_wrapper<const double> >("Ceu_3311", boost::cref(Ceu_3311)));
        ModelParamMap.insert(std::pair<std::string, boost::reference_wrapper<const double> >("Ced_1111", boost::cref(Ced_1111)));
        ModelParamMap.insert(std::pair<std::string, boost::reference_wrapper<const double> >("Ced_1122", boost::cref(Ced_1122)));
        ModelParamMap.insert(std::pair<std::string, boost::reference_wrapper<const double> >("Ced_2211", boost::cref(Ced_2211)));
        ModelParamMap.insert(std::pair<std::string, boost::reference_wrapper<const double> >("Ced_1133", boost::cref(Ced_1133)));
        ModelParamMap.insert(std::pair<std::string, boost::reference_wrapper<const double> >("Ced_3311", boost::cref(Ced_3311)));
        ModelParamMap.insert(std::pair<std::string, boost::reference_wrapper<const double> >("CLe_1111", boost::cref(CLe_1111)));
        ModelParamMap.insert(std::pair<std::string, boost::reference_wrapper<const double> >("CLe_1122", boost::cref(CLe_1122)));
        ModelParamMap.insert(std::pair<std::string, boost::reference_wrapper<const double> >("CLe_2211", boost::cref(CLe_2211)));
        ModelParamMap.insert(std::pair<std::string, boost::reference_wrapper<const double> >("CLe_1133", boost::cref(CLe_1133)));
        ModelParamMap.insert(std::pair<std::string, boost::reference_wrapper<const double> >("CLe_3311", boost::cref(CLe_3311)));
        ModelParamMap.insert(std::pair<std::string, boost::reference_wrapper<const double> >("CLu_1111", boost::cref(CLu_1111)));
        ModelParamMap.insert(std::pair<std::string, boost::reference_wrapper<const double> >("CLu_1122", boost::cref(CLu_1122)));
        ModelParamMap.insert(std::pair<std::string, boost::reference_wrapper<const double> >("CLu_2211", boost::cref(CLu_2211)));
        ModelParamMap.insert(std::pair<std::string, boost::reference_wrapper<const double> >("CLu_1133", boost::cref(CLu_1133)));
        ModelParamMap.insert(std::pair<std::string, boost::reference_wrapper<const double> >("CLu_3311", boost::cref(CLu_3311)));
        ModelParamMap.insert(std::pair<std::string, boost::reference_wrapper<const double> >("CLd_1111", boost::cref(CLd_1111)));
        ModelParamMap.insert(std::pair<std::string, boost::reference_wrapper<const double> >("CLd_1122", boost::cref(CLd_1122)));
        ModelParamMap.insert(std::pair<std::string, boost::reference_wrapper<const double> >("CLd_2211", boost::cref(CLd_2211)));
        ModelParamMap.insert(std::pair<std::string, boost::reference_wrapper<const double> >("CLd_1133", boost::cref(CLd_1133)));
        ModelParamMap.insert(std::pair<std::string, boost::reference_wrapper<const double> >("CLd_3311", boost::cref(CLd_3311)));
        ModelParamMap.insert(std::pair<std::string, boost::reference_wrapper<const double> >("CQe_1111", boost::cref(CQe_1111)));
        ModelParamMap.insert(std::pair<std::string, boost::reference_wrapper<const double> >("CQe_1122", boost::cref(CQe_1122)));
        ModelParamMap.insert(std::pair<std::string, boost::reference_wrapper<const double> >("CQe_2211", boost::cref(CQe_2211)));
        ModelParamMap.insert(std::pair<std::string, boost::reference_wrapper<const double> >("CQe_1133", boost::cref(CQe_1133)));
        ModelParamMap.insert(std::pair<std::string, boost::reference_wrapper<const double> >("CQe_3311", boost::cref(CQe_3311)));
    }
    ModelParamMap.insert(std::pair<std::string, boost::reference_wrapper<const double> >("Lambda_NP", boost::cref(Lambda_NP)));
    ModelParamMap.insert(std::pair<std::string, boost::reference_wrapper<const double> >("eVBF2_HZZ1", boost::cref(eVBF2_HZZ1)));
    ModelParamMap.insert(std::pair<std::string, boost::reference_wrapper<const double> >("eVBF2_HZZ2", boost::cref(eVBF2_HZZ2)));
    ModelParamMap.insert(std::pair<std::string, boost::reference_wrapper<const double> >("eVBF2_HZZ3", boost::cref(eVBF2_HZZ3)));
    ModelParamMap.insert(std::pair<std::string, boost::reference_wrapper<const double> >("eVBF2_HZA1", boost::cref(eVBF2_HZA1)));
    ModelParamMap.insert(std::pair<std::string, boost::reference_wrapper<const double> >("eVBF2_HZA2", boost::cref(eVBF2_HZA2)));
    ModelParamMap.insert(std::pair<std::string, boost::reference_wrapper<const double> >("eVBF2_HAA", boost::cref(eVBF2_HAA)));
    ModelParamMap.insert(std::pair<std::string, boost::reference_wrapper<const double> >("eVBF2_HWW1", boost::cref(eVBF2_HWW1)));
    ModelParamMap.insert(std::pair<std::string, boost::reference_wrapper<const double> >("eVBF2_HWW2", boost::cref(eVBF2_HWW2)));
    ModelParamMap.insert(std::pair<std::string, boost::reference_wrapper<const double> >("eVBF2_HWW3", boost::cref(eVBF2_HWW3)));
    ModelParamMap.insert(std::pair<std::string, boost::reference_wrapper<const double> >("eVBF2_Hgg", boost::cref(eVBF2_Hgg)));
    ModelParamMap.insert(std::pair<std::string, boost::reference_wrapper<const double> >("eVBF2_HZuL", boost::cref(eVBF2_HZuL)));
    ModelParamMap.insert(std::pair<std::string, boost::reference_wrapper<const double> >("eVBF2_HZuR", boost::cref(eVBF2_HZuR)));
    ModelParamMap.insert(std::pair<std::string, boost::reference_wrapper<const double> >("eVBF2_HZdL", boost::cref(eVBF2_HZdL)));
    ModelParamMap.insert(std::pair<std::string, boost::reference_wrapper<const double> >("eVBF2_HZdR", boost::cref(eVBF2_HZdR)));
    ModelParamMap.insert(std::pair<std::string, boost::reference_wrapper<const double> >("eVBF2_HWud", boost::cref(eVBF2_HWud)));
    ModelParamMap.insert(std::pair<std::string, boost::reference_wrapper<const double> >("eVBF2_ZuL", boost::cref(eVBF2_ZuL)));
    ModelParamMap.insert(std::pair<std::string, boost::reference_wrapper<const double> >("eVBF2_ZuR", boost::cref(eVBF2_ZuR)));
    ModelParamMap.insert(std::pair<std::string, boost::reference_wrapper<const double> >("eVBF2_ZdL", boost::cref(eVBF2_ZdL)));
    ModelParamMap.insert(std::pair<std::string, boost::reference_wrapper<const double> >("eVBF2_ZdR", boost::cref(eVBF2_ZdR)));
    ModelParamMap.insert(std::pair<std::string, boost::reference_wrapper<const double> >("eVBF2_Wud", boost::cref(eVBF2_Wud)));
    ModelParamMap.insert(std::pair<std::string, boost::reference_wrapper<const double> >("eVBF78_HZZ1", boost::cref(eVBF78_HZZ1)));
    ModelParamMap.insert(std::pair<std::string, boost::reference_wrapper<const double> >("eVBF78_HZZ2", boost::cref(eVBF78_HZZ2)));
    ModelParamMap.insert(std::pair<std::string, boost::reference_wrapper<const double> >("eVBF78_HZZ3", boost::cref(eVBF78_HZZ3)));
    ModelParamMap.insert(std::pair<std::string, boost::reference_wrapper<const double> >("eVBF78_HZA1", boost::cref(eVBF78_HZA1)));
    ModelParamMap.insert(std::pair<std::string, boost::reference_wrapper<const double> >("eVBF78_HZA2", boost::cref(eVBF78_HZA2)));
    ModelParamMap.insert(std::pair<std::string, boost::reference_wrapper<const double> >("eVBF78_HAA", boost::cref(eVBF78_HAA)));
    ModelParamMap.insert(std::pair<std::string, boost::reference_wrapper<const double> >("eVBF78_HWW1", boost::cref(eVBF78_HWW1)));
    ModelParamMap.insert(std::pair<std::string, boost::reference_wrapper<const double> >("eVBF78_HWW2", boost::cref(eVBF78_HWW2)));
    ModelParamMap.insert(std::pair<std::string, boost::reference_wrapper<const double> >("eVBF78_HWW3", boost::cref(eVBF78_HWW3)));
    ModelParamMap.insert(std::pair<std::string, boost::reference_wrapper<const double> >("eVBF78_Hgg", boost::cref(eVBF78_Hgg)));
    ModelParamMap.insert(std::pair<std::string, boost::reference_wrapper<const double> >("eVBF78_HZuL", boost::cref(eVBF78_HZuL)));
    ModelParamMap.insert(std::pair<std::string, boost::reference_wrapper<const double> >("eVBF78_HZuR", boost::cref(eVBF78_HZuR)));
    ModelParamMap.insert(std::pair<std::string, boost::reference_wrapper<const double> >("eVBF78_HZdL", boost::cref(eVBF78_HZdL)));
    ModelParamMap.insert(std::pair<std::string, boost::reference_wrapper<const double> >("eVBF78_HZdR", boost::cref(eVBF78_HZdR)));
    ModelParamMap.insert(std::pair<std::string, boost::reference_wrapper<const double> >("eVBF78_HWud", boost::cref(eVBF78_HWud)));
    ModelParamMap.insert(std::pair<std::string, boost::reference_wrapper<const double> >("eVBF78_ZuL", boost::cref(eVBF78_ZuL)));
    ModelParamMap.insert(std::pair<std::string, boost::reference_wrapper<const double> >("eVBF78_ZuR", boost::cref(eVBF78_ZuR)));
    ModelParamMap.insert(std::pair<std::string, boost::reference_wrapper<const double> >("eVBF78_ZdL", boost::cref(eVBF78_ZdL)));
    ModelParamMap.insert(std::pair<std::string, boost::reference_wrapper<const double> >("eVBF78_ZdR", boost::cref(eVBF78_ZdR)));
    ModelParamMap.insert(std::pair<std::string, boost::reference_wrapper<const double> >("eVBF78_Wud", boost::cref(eVBF78_Wud)));
    ModelParamMap.insert(std::pair<std::string, boost::reference_wrapper<const double> >("eWH2_HWW1", boost::cref(eWH2_HWW1)));
    ModelParamMap.insert(std::pair<std::string, boost::reference_wrapper<const double> >("eWH2_HWW2", boost::cref(eWH2_HWW2)));
    ModelParamMap.insert(std::pair<std::string, boost::reference_wrapper<const double> >("eWH2_HWW3", boost::cref(eWH2_HWW3)));
    ModelParamMap.insert(std::pair<std::string, boost::reference_wrapper<const double> >("eWH2_HWud", boost::cref(eWH2_HWud)));
    ModelParamMap.insert(std::pair<std::string, boost::reference_wrapper<const double> >("eWH2_Wud", boost::cref(eWH2_Wud)));
    ModelParamMap.insert(std::pair<std::string, boost::reference_wrapper<const double> >("eWH78_HWW1", boost::cref(eWH78_HWW1)));
    ModelParamMap.insert(std::pair<std::string, boost::reference_wrapper<const double> >("eWH78_HWW2", boost::cref(eWH78_HWW2)));
    ModelParamMap.insert(std::pair<std::string, boost::reference_wrapper<const double> >("eWH78_HWW3", boost::cref(eWH78_HWW3)));
    ModelParamMap.insert(std::pair<std::string, boost::reference_wrapper<const double> >("eWH78_HWud", boost::cref(eWH78_HWud)));
    ModelParamMap.insert(std::pair<std::string, boost::reference_wrapper<const double> >("eWH78_Wud", boost::cref(eWH78_Wud)));
    ModelParamMap.insert(std::pair<std::string, boost::reference_wrapper<const double> >("eZH2_HZZ1", boost::cref(eZH2_HZZ1)));
    ModelParamMap.insert(std::pair<std::string, boost::reference_wrapper<const double> >("eZH2_HZZ2", boost::cref(eZH2_HZZ2)));
    ModelParamMap.insert(std::pair<std::string, boost::reference_wrapper<const double> >("eZH2_HZZ3", boost::cref(eZH2_HZZ3)));
    ModelParamMap.insert(std::pair<std::string, boost::reference_wrapper<const double> >("eZH2_HZA1", boost::cref(eZH2_HZA1)));
    ModelParamMap.insert(std::pair<std::string, boost::reference_wrapper<const double> >("eZH2_HZA2", boost::cref(eZH2_HZA2)));
    ModelParamMap.insert(std::pair<std::string, boost::reference_wrapper<const double> >("eZH2_HZuL", boost::cref(eZH2_HZuL)));
    ModelParamMap.insert(std::pair<std::string, boost::reference_wrapper<const double> >("eZH2_HZuR", boost::cref(eZH2_HZuR)));
    ModelParamMap.insert(std::pair<std::string, boost::reference_wrapper<const double> >("eZH2_HZdL", boost::cref(eZH2_HZdL)));
    ModelParamMap.insert(std::pair<std::string, boost::reference_wrapper<const double> >("eZH2_HZdR", boost::cref(eZH2_HZdR)));
    ModelParamMap.insert(std::pair<std::string, boost::reference_wrapper<const double> >("eZH2_ZuL", boost::cref(eZH2_ZuL)));
    ModelParamMap.insert(std::pair<std::string, boost::reference_wrapper<const double> >("eZH2_ZuR", boost::cref(eZH2_ZuR)));
    ModelParamMap.insert(std::pair<std::string, boost::reference_wrapper<const double> >("eZH2_ZdL", boost::cref(eZH2_ZdL)));
    ModelParamMap.insert(std::pair<std::string, boost::reference_wrapper<const double> >("eZH2_ZdR", boost::cref(eZH2_ZdR)));
    ModelParamMap.insert(std::pair<std::string, boost::reference_wrapper<const double> >("eZH78_HZZ1", boost::cref(eZH78_HZZ1)));
    ModelParamMap.insert(std::pair<std::string, boost::reference_wrapper<const double> >("eZH78_HZZ2", boost::cref(eZH78_HZZ2)));
    ModelParamMap.insert(std::pair<std::string, boost::reference_wrapper<const double> >("eZH78_HZZ3", boost::cref(eZH78_HZZ3)));
    ModelParamMap.insert(std::pair<std::string, boost::reference_wrapper<const double> >("eZH78_HZA1", boost::cref(eZH78_HZA1)));
    ModelParamMap.insert(std::pair<std::string, boost::reference_wrapper<const double> >("eZH78_HZA2", boost::cref(eZH78_HZA2)));
    ModelParamMap.insert(std::pair<std::string, boost::reference_wrapper<const double> >("eZH78_HZuL", boost::cref(eZH78_HZuL)));
    ModelParamMap.insert(std::pair<std::string, boost::reference_wrapper<const double> >("eZH78_HZuR", boost::cref(eZH78_HZuR)));
    ModelParamMap.insert(std::pair<std::string, boost::reference_wrapper<const double> >("eZH78_HZdL", boost::cref(eZH78_HZdL)));
    ModelParamMap.insert(std::pair<std::string, boost::reference_wrapper<const double> >("eZH78_HZdR", boost::cref(eZH78_HZdR)));
    ModelParamMap.insert(std::pair<std::string, boost::reference_wrapper<const double> >("eZH78_ZuL", boost::cref(eZH78_ZuL)));
    ModelParamMap.insert(std::pair<std::string, boost::reference_wrapper<const double> >("eZH78_ZuR", boost::cref(eZH78_ZuR)));
    ModelParamMap.insert(std::pair<std::string, boost::reference_wrapper<const double> >("eZH78_ZdL", boost::cref(eZH78_ZdL)));
    ModelParamMap.insert(std::pair<std::string, boost::reference_wrapper<const double> >("eZH78_ZdR", boost::cref(eZH78_ZdR)));
    ModelParamMap.insert(std::pair<std::string, boost::reference_wrapper<const double> >("ettH2_Htt", boost::cref(ettH2_Htt)));
    ModelParamMap.insert(std::pair<std::string, boost::reference_wrapper<const double> >("ettH2_Hgg", boost::cref(ettH2_Hgg)));
    ModelParamMap.insert(std::pair<std::string, boost::reference_wrapper<const double> >("ettH78_Htt", boost::cref(ettH78_Htt)));
    ModelParamMap.insert(std::pair<std::string, boost::reference_wrapper<const double> >("ettH78_Hgg", boost::cref(ettH78_Hgg)));
    if (FlagMwInput) {
        ModelParamMap.insert(std::pair<std::string, boost::reference_wrapper<const double> >("MwInput", boost::cref(MwInput)));
    }

    if (FlagLeptonUniversal) {        
//  bsll/sbll entries only interesting (for the moment) if non-lepton universal. Set to 0 otherwise
        CLQ1_1123 = 0.0;
        CLQ1_2223 = 0.0;
        CLQ1_3323 = 0.0;
        CLQ1_1132 = 0.0;
        CLQ1_2232 = 0.0;
        CLQ1_3332 = 0.0;
        
        CLQ3_1123 = 0.0;
        CLQ3_2223 = 0.0;
        CLQ3_3323 = 0.0;
        CLQ3_1132 = 0.0;
        CLQ3_2232 = 0.0;
        CLQ3_3332 = 0.0;
        
        Ced_1123 = 0.0;
        Ced_2223 = 0.0;
        Ced_3323 = 0.0;
        Ced_1132 = 0.0;
        Ced_2232 = 0.0;
        Ced_3332 = 0.0;  
        
        CLd_1123 = 0.0;
        CLd_2223 = 0.0;
        CLd_3323 = 0.0;
        CLd_1132 = 0.0;
        CLd_2232 = 0.0;
        CLd_3332 = 0.0;
        
        CQe_2311 = 0.0;
        CQe_2322 = 0.0;
        CQe_2333 = 0.0;
        CQe_3211 = 0.0;
        CQe_3222 = 0.0;
        CQe_3233 = 0.0;
    }
    
}

bool NPEffectiveGIMR::PostUpdate()
{
    if (!NPbase::PostUpdate()) return (false);

    LambdaNP2 = Lambda_NP * Lambda_NP;
    v2_over_LambdaNP2 = v() * v() / LambdaNP2;
    if (FlagMwInput)
        cW_tree = MwInput / Mz;
    else
        cW_tree = Mw_tree() / Mz;
    cW2_tree = cW_tree * cW_tree;
    sW2_tree = 1.0 - cW2_tree;
    sW_tree = sqrt(sW2_tree);

    delta_ZZ = (cW2_tree * CHW + sW2_tree * CHB + sW_tree * cW_tree * CHWB) * v2_over_LambdaNP2;
    delta_AA = (sW2_tree * CHW + cW2_tree * CHB - sW_tree * cW_tree * CHWB) * v2_over_LambdaNP2;
    delta_AZ = 2.0 * sW_tree * cW_tree * (CHW - CHB) * v2_over_LambdaNP2
            - (cW2_tree - sW2_tree) * CHWB * v2_over_LambdaNP2;
    delta_h = (-CHD / 4.0 + CHbox) * v2_over_LambdaNP2;

    return (true);
}

void NPEffectiveGIMR::setParameter(const std::string name, const double& value)
{
    if (name.compare("CG") == 0)
        CG = value;
    else if (name.compare("CW") == 0)
        CW = value;
    else if (name.compare("CHG") == 0)
        CHG = value;
    else if (name.compare("CHW") == 0)
        CHW = value;
    else if (name.compare("CHB") == 0)
        CHB = value;
    else if (name.compare("CHWB") == 0)
        CHWB = value;
    else if (name.compare("CHD") == 0)
        CHD = value;
    else if (name.compare("CHbox") == 0)
        CHbox = value;
    else if (name.compare("CH") == 0)
        CH = value;
    else if (name.compare("CHL1_11") == 0)
        CHL1_11 = value;
    else if (name.compare("CHL1_12r") == 0)
        CHL1_12r = value;
    else if (name.compare("CHL1_13r") == 0)
        CHL1_13r = value;
    else if (name.compare("CHL1_22") == 0)
        CHL1_22 = value;
    else if (name.compare("CHL1_23r") == 0)
        CHL1_23r = value;
    else if (name.compare("CHL1_33") == 0)
        CHL1_33 = value;
    else if (name.compare("CHL1_12i") == 0)
        CHL1_12i = value;
    else if (name.compare("CHL1_13i") == 0)
        CHL1_13i = value;
    else if (name.compare("CHL1_23i") == 0)
        CHL1_23i = value;
    else if (name.compare("CHL1") == 0) {
        CHL1_11 = value;
        CHL1_12r = 0.0;
        CHL1_13r = 0.0;
        CHL1_22 = value;
        CHL1_23r = 0.0;
        CHL1_33 = value;
        CHL1_12i = 0.0;
        CHL1_13i = 0.0;
        CHL1_23i = 0.0;
    } else if (name.compare("CHL3_11") == 0)
        CHL3_11 = value;
    else if (name.compare("CHL3_12r") == 0)
        CHL3_12r = value;
    else if (name.compare("CHL3_13r") == 0)
        CHL3_13r = value;
    else if (name.compare("CHL3_22") == 0)
        CHL3_22 = value;
    else if (name.compare("CHL3_23r") == 0)
        CHL3_23r = value;
    else if (name.compare("CHL3_33") == 0)
        CHL3_33 = value;
    else if (name.compare("CHL3_12i") == 0)
        CHL3_12i = value;
    else if (name.compare("CHL3_13i") == 0)
        CHL3_13i = value;
    else if (name.compare("CHL3_23i") == 0)
        CHL3_23i = value;
    else if (name.compare("CHL3") == 0) {
        CHL3_11 = value;
        CHL3_12r = 0.0;
        CHL3_13r = 0.0;
        CHL3_22 = value;
        CHL3_23r = 0.0;
        CHL3_33 = value;
        CHL3_12i = 0.0;
        CHL3_13i = 0.0;
        CHL3_23i = 0.0;
    } else if (name.compare("CHe_11") == 0)
        CHe_11 = value;
    else if (name.compare("CHe_12r") == 0)
        CHe_12r = value;
    else if (name.compare("CHe_13r") == 0)
        CHe_13r = value;
    else if (name.compare("CHe_22") == 0)
        CHe_22 = value;
    else if (name.compare("CHe_23r") == 0)
        CHe_23r = value;
    else if (name.compare("CHe_33") == 0)
        CHe_33 = value;
    else if (name.compare("CHe_12i") == 0)
        CHe_12i = value;
    else if (name.compare("CHe_13i") == 0)
        CHe_13i = value;
    else if (name.compare("CHe_23i") == 0)
        CHe_23i = value;
    else if (name.compare("CHe") == 0) {
        CHe_11 = value;
        CHe_12r = 0.0;
        CHe_13r = 0.0;
        CHe_22 = value;
        CHe_23r = 0.0;
        CHe_33 = value;
        CHe_12i = 0.0;
        CHe_13i = 0.0;
        CHe_23i = 0.0;
    } else if (name.compare("CHQ1_11") == 0)
        CHQ1_11 = value;
    else if (name.compare("CHQ1_12r") == 0)
        CHQ1_12r = value;
    else if (name.compare("CHQ1_13r") == 0)
        CHQ1_13r = value;
    else if (name.compare("CHQ1_22") == 0)
        CHQ1_22 = value;
    else if (name.compare("CHQ1_23r") == 0)
        CHQ1_23r = value;
    else if (name.compare("CHQ1_33") == 0)
        CHQ1_33 = value;
    else if (name.compare("CHQ1_12i") == 0)
        CHQ1_12i = value;
    else if (name.compare("CHQ1_13i") == 0)
        CHQ1_13i = value;
    else if (name.compare("CHQ1_23i") == 0)
        CHQ1_23i = value;
    else if (name.compare("CHQ1") == 0) {
        CHQ1_11 = value;
        CHQ1_12r = 0.0;
        CHQ1_13r = 0.0;
        CHQ1_22 = value;
        CHQ1_23r = 0.0;
        CHQ1_33 = value;
        CHQ1_12i = 0.0;
        CHQ1_13i = 0.0;
        CHQ1_23i = 0.0;
    } else if (name.compare("CHQ3_11") == 0)
        CHQ3_11 = value;
    else if (name.compare("CHQ3_12r") == 0)
        CHQ3_12r = value;
    else if (name.compare("CHQ3_13r") == 0)
        CHQ3_13r = value;
    else if (name.compare("CHQ3_22") == 0)
        CHQ3_22 = value;
    else if (name.compare("CHQ3_23r") == 0)
        CHQ3_23r = value;
    else if (name.compare("CHQ3_33") == 0)
        CHQ3_33 = value;
    else if (name.compare("CHQ3_12i") == 0)
        CHQ3_12i = value;
    else if (name.compare("CHQ3_13i") == 0)
        CHQ3_13i = value;
    else if (name.compare("CHQ3_23i") == 0)
        CHQ3_23i = value;
    else if (name.compare("CHQ3") == 0) {
        CHQ3_11 = value;
        CHQ3_12r = 0.0;
        CHQ3_13r = 0.0;
        CHQ3_22 = value;
        CHQ3_23r = 0.0;
        CHQ3_33 = value;
        CHQ3_12i = 0.0;
        CHQ3_13i = 0.0;
        CHQ3_23i = 0.0;
    } else if (name.compare("CHu_11") == 0)
        CHu_11 = value;
    else if (name.compare("CHu_12r") == 0)
        CHu_12r = value;
    else if (name.compare("CHu_13r") == 0)
        CHu_13r = value;
    else if (name.compare("CHu_22") == 0)
        CHu_22 = value;
    else if (name.compare("CHu_23r") == 0)
        CHu_23r = value;
    else if (name.compare("CHu_33") == 0)
        CHu_33 = value;
    else if (name.compare("CHu_12i") == 0)
        CHu_12i = value;
    else if (name.compare("CHu_13i") == 0)
        CHu_13i = value;
    else if (name.compare("CHu_23i") == 0)
        CHu_23i = value;
    else if (name.compare("CHu") == 0) {
        CHu_11 = value;
        CHu_12r = 0.0;
        CHu_13r = 0.0;
        CHu_22 = value;
        CHu_23r = 0.0;
        CHu_33 = value;
        CHu_12i = 0.0;
        CHu_13i = 0.0;
        CHu_23i = 0.0;
    } else if (name.compare("CHd_11") == 0)
        CHd_11 = value;
    else if (name.compare("CHd_12r") == 0)
        CHd_12r = value;
    else if (name.compare("CHd_13r") == 0)
        CHd_13r = value;
    else if (name.compare("CHd_22") == 0)
        CHd_22 = value;
    else if (name.compare("CHd_23r") == 0)
        CHd_23r = value;
    else if (name.compare("CHd_33") == 0)
        CHd_33 = value;
    else if (name.compare("CHd_12i") == 0)
        CHd_12i = value;
    else if (name.compare("CHd_13i") == 0)
        CHd_13i = value;
    else if (name.compare("CHd_23i") == 0)
        CHd_23i = value;
    else if (name.compare("CHd") == 0) {
        CHd_11 = value;
        CHd_12r = 0.0;
        CHd_13r = 0.0;
        CHd_22 = value;
        CHd_23r = 0.0;
        CHd_33 = value;
        CHd_12i = 0.0;
        CHd_13i = 0.0;
        CHd_23i = 0.0;
    } else if (name.compare("CHud_11r") == 0)
        CHud_11r = value;
    else if (name.compare("CHud_12r") == 0)
        CHud_12r = value;
    else if (name.compare("CHud_13r") == 0)
        CHud_13r = value;
    else if (name.compare("CHud_22r") == 0)
        CHud_22r = value;
    else if (name.compare("CHud_23r") == 0)
        CHud_23r = value;
    else if (name.compare("CHud_33r") == 0)
        CHud_33r = value;
    else if (name.compare("CHud_r") == 0) {
        CHud_11r = value;
        CHud_12r = 0.0;
        CHud_13r = 0.0;
        CHud_22r = value;
        CHud_23r = 0.0;
        CHud_33r = value;
    } else if (name.compare("CHud_11i") == 0)
        CHud_11i = value;
    else if (name.compare("CHud_12i") == 0)
        CHud_12i = value;
    else if (name.compare("CHud_13i") == 0)
        CHud_13i = value;
    else if (name.compare("CHud_22i") == 0)
        CHud_22i = value;
    else if (name.compare("CHud_23i") == 0)
        CHud_23i = value;
    else if (name.compare("CHud_33i") == 0)
        CHud_33i = value;
    else if (name.compare("CHud_i") == 0) {
        CHud_11i = value;
        CHud_12i = 0.0;
        CHud_13i = 0.0;
        CHud_22i = value;
        CHud_23i = 0.0;
        CHud_33i = value;
    } else if (name.compare("CeH_11r") == 0)
        CeH_11r = value;
    else if (name.compare("CeH_12r") == 0)
        CeH_12r = value;
    else if (name.compare("CeH_13r") == 0)
        CeH_13r = value;
    else if (name.compare("CeH_22r") == 0)
        CeH_22r = value;
    else if (name.compare("CeH_23r") == 0)
        CeH_23r = value;
    else if (name.compare("CeH_33r") == 0)
        CeH_33r = value;
    else if (name.compare("CeH_r") == 0) {
        CeH_11r = value;
        CeH_12r = 0.0;
        CeH_13r = 0.0;
        CeH_22r = value;
        CeH_23r = 0.0;
        CeH_33r = value;
    } else if (name.compare("CeH_11i") == 0)
        CeH_11i = value;
    else if (name.compare("CeH_12i") == 0)
        CeH_12i = value;
    else if (name.compare("CeH_13i") == 0)
        CeH_13i = value;
    else if (name.compare("CeH_22i") == 0)
        CeH_22i = value;
    else if (name.compare("CeH_23i") == 0)
        CeH_23i = value;
    else if (name.compare("CeH_33i") == 0)
        CeH_33i = value;
    else if (name.compare("CeH_i") == 0) {
        CeH_11i = value;
        CeH_12i = 0.0;
        CeH_13i = 0.0;
        CeH_22i = value;
        CeH_23i = 0.0;
        CeH_33i = value;
    } else if (name.compare("CuH_11r") == 0)
        CuH_11r = value;
    else if (name.compare("CuH_12r") == 0)
        CuH_12r = value;
    else if (name.compare("CuH_13r") == 0)
        CuH_13r = value;
    else if (name.compare("CuH_22r") == 0)
        CuH_22r = value;
    else if (name.compare("CuH_23r") == 0)
        CuH_23r = value;
    else if (name.compare("CuH_33r") == 0)
        CuH_33r = value;
    else if (name.compare("CuH_r") == 0) {
        CuH_11r = value;
        CuH_12r = 0.0;
        CuH_13r = 0.0;
        CuH_22r = value;
        CuH_23r = 0.0;
        CuH_33r = value;
    } else if (name.compare("CuH_11i") == 0)
        CuH_11i = value;
    else if (name.compare("CuH_12i") == 0)
        CuH_12i = value;
    else if (name.compare("CuH_13i") == 0)
        CuH_13i = value;
    else if (name.compare("CuH_22i") == 0)
        CuH_22i = value;
    else if (name.compare("CuH_23i") == 0)
        CuH_23i = value;
    else if (name.compare("CuH_33i") == 0)
        CuH_33i = value;
    else if (name.compare("CuH_i") == 0) {
        CuH_11i = value;
        CuH_12i = 0.0;
        CuH_13i = 0.0;
        CuH_22i = value;
        CuH_23i = 0.0;
        CuH_33i = value;
    } else if (name.compare("CdH_11r") == 0)
        CdH_11r = value;
    else if (name.compare("CdH_12r") == 0)
        CdH_12r = value;
    else if (name.compare("CdH_13r") == 0)
        CdH_13r = value;
    else if (name.compare("CdH_22r") == 0)
        CdH_22r = value;
    else if (name.compare("CdH_23r") == 0)
        CdH_23r = value;
    else if (name.compare("CdH_33r") == 0)
        CdH_33r = value;
    else if (name.compare("CdH_r") == 0) {
        CdH_11r = value;
        CdH_12r = 0.0;
        CdH_13r = 0.0;
        CdH_22r = value;
        CdH_23r = 0.0;
        CdH_33r = value;
    } else if (name.compare("CdH_11i") == 0)
        CdH_11i = value;
    else if (name.compare("CdH_12i") == 0)
        CdH_12i = value;
    else if (name.compare("CdH_13i") == 0)
        CdH_13i = value;
    else if (name.compare("CdH_22i") == 0)
        CdH_22i = value;
    else if (name.compare("CdH_23i") == 0)
        CdH_23i = value;
    else if (name.compare("CdH_33i") == 0)
        CdH_33i = value;
    else if (name.compare("CdH_i") == 0) {
        CdH_11i = value;
        CdH_12i = 0.0;
        CdH_13i = 0.0;
        CdH_22i = value;
        CdH_23i = 0.0;
        CdH_33i = value;
    } else if (name.compare("CuG_11r") == 0)
        CuG_11r = value;
    else if (name.compare("CuG_12r") == 0)
        CuG_12r = value;
    else if (name.compare("CuG_13r") == 0)
        CuG_13r = value;
    else if (name.compare("CuG_22r") == 0)
        CuG_22r = value;
    else if (name.compare("CuG_23r") == 0)
        CuG_23r = value;
    else if (name.compare("CuG_33r") == 0)
        CuG_33r = value;
    else if (name.compare("CuG_r") == 0) {
        CuG_11r = value;
        CuG_12r = 0.0;
        CuG_13r = 0.0;
        CuG_22r = value;
        CuG_23r = 0.0;
        CuG_33r = value;
    } else if (name.compare("CuG_11i") == 0)
        CuG_11i = value;
    else if (name.compare("CuG_12i") == 0)
        CuG_12i = value;
    else if (name.compare("CuG_13i") == 0)
        CuG_13i = value;
    else if (name.compare("CuG_22i") == 0)
        CuG_22i = value;
    else if (name.compare("CuG_23i") == 0)
        CuG_23i = value;
    else if (name.compare("CuG_33i") == 0)
        CuG_33i = value;
    else if (name.compare("CuG_i") == 0) {
        CuG_11i = value;
        CuG_12i = 0.0;
        CuG_13i = 0.0;
        CuG_22i = value;
        CuG_23i = 0.0;
        CuG_33i = value;
    } else if (name.compare("CuW_11r") == 0)
        CuW_11r = value;
    else if (name.compare("CuW_12r") == 0)
        CuW_12r = value;
    else if (name.compare("CuW_13r") == 0)
        CuW_13r = value;
    else if (name.compare("CuW_22r") == 0)
        CuW_22r = value;
    else if (name.compare("CuW_23r") == 0)
        CuW_23r = value;
    else if (name.compare("CuW_33r") == 0)
        CuW_33r = value;
    else if (name.compare("CuW_r") == 0) {
        CuW_11r = value;
        CuW_12r = 0.0;
        CuW_13r = 0.0;
        CuW_22r = value;
        CuW_23r = 0.0;
        CuW_33r = value;
    } else if (name.compare("CuW_11i") == 0)
        CuW_11i = value;
    else if (name.compare("CuW_12i") == 0)
        CuW_12i = value;
    else if (name.compare("CuW_13i") == 0)
        CuW_13i = value;
    else if (name.compare("CuW_22i") == 0)
        CuW_22i = value;
    else if (name.compare("CuW_23i") == 0)
        CuW_23i = value;
    else if (name.compare("CuW_33i") == 0)
        CuW_33i = value;
    else if (name.compare("CuW_i") == 0) {
        CuW_11i = value;
        CuW_12i = 0.0;
        CuW_13i = 0.0;
        CuW_22i = value;
        CuW_23i = 0.0;
        CuW_33i = value;
    } else if (name.compare("CuB_11r") == 0)
        CuB_11r = value;
    else if (name.compare("CuB_12r") == 0)
        CuB_12r = value;
    else if (name.compare("CuB_13r") == 0)
        CuB_13r = value;
    else if (name.compare("CuB_22r") == 0)
        CuB_22r = value;
    else if (name.compare("CuB_23r") == 0)
        CuB_23r = value;
    else if (name.compare("CuB_33r") == 0)
        CuB_33r = value;
    else if (name.compare("CuB_r") == 0) {
        CuB_11r = value;
        CuB_12r = 0.0;
        CuB_13r = 0.0;
        CuB_22r = value;
        CuB_23r = 0.0;
        CuB_33r = value;
    } else if (name.compare("CuB_11i") == 0)
        CuB_11i = value;
    else if (name.compare("CuB_12i") == 0)
        CuB_12i = value;
    else if (name.compare("CuB_13i") == 0)
        CuB_13i = value;
    else if (name.compare("CuB_22i") == 0)
        CuB_22i = value;
    else if (name.compare("CuB_23i") == 0)
        CuB_23i = value;
    else if (name.compare("CuB_33i") == 0)
        CuB_33i = value;
    else if (name.compare("CuB_i") == 0) {
        CuB_11i = value;
        CuB_12i = 0.0;
        CuB_13i = 0.0;
        CuB_22i = value;
        CuB_23i = 0.0;
        CuB_33i = value;
    } else if (name.compare("CLL_1111") == 0) {
        CLL_1111 = value;
    } else if (name.compare("CLL_1221") == 0) {
        CLL_1221 = value;
        CLL_2112 = value;
        CLL_2211 = value;
        CLL_1122 = value;
    } else if (name.compare("CLL_1133") == 0) {
        CLL_1331 = value;
        CLL_3113 = value;
        CLL_3311 = value;
        CLL_1133 = value;
    } else if (name.compare("CLL") == 0) {
        CLL_1111 = value;
        CLL_1221 = value;
        CLL_2112 = value;
        CLL_2211 = value;
        CLL_1122 = value;
        CLL_3311 = value;
        CLL_1133 = value;
        CLL_1331 = value;
        CLL_3113 = value;
    } else if (name.compare("CLQ1_1111") == 0) {
        CLQ1_1111 = value;
    } else if (name.compare("CLQ1_1122") == 0) {
        CLQ1_1122 = value;
        CLQ1_2211 = value;
        CLQ1_1221 = value;
        CLQ1_2112 = value;
    } else if (name.compare("CLQ1_1133") == 0) {
        CLQ1_1133 = value;
        CLQ1_3311 = value;
        CLQ1_1331 = value;
        CLQ1_3113 = value;
    } else if (name.compare("CLQ1") == 0) {
        CLQ1_1111 = value;
        CLQ1_1122 = value;
        CLQ1_2211 = value;
        CLQ1_1221 = value;
        CLQ1_2112 = value;
        CLQ1_1133 = value;
        CLQ1_3311 = value;
        CLQ1_1331 = value;
        CLQ1_3113 = value;
    } else if (name.compare("CLQ1_1111") == 0) {
        CLQ1_1111 = value;
    } else if (name.compare("CLQ1_1122") == 0) {
        CLQ1_1122 = value;
        CLQ1_2211 = value;
        CLQ1_2112 = value;
        CLQ1_1221 = value;
    } else if (name.compare("CLQ1_1133") == 0) {
        CLQ1_1133 = value;
        CLQ1_3311 = value;
        CLQ1_3113 = value;
        CLQ1_1331 = value;
    } else if (name.compare("CLQ1_1123") == 0) {
        CLQ1_1123 = value;
    } else if (name.compare("CLQ1_2223") == 0) {
        CLQ1_2223 = value;
    } else if (name.compare("CLQ1_3323") == 0) {
        CLQ1_3323 = value;
    } else if (name.compare("CLQ1_1132") == 0) {
        CLQ1_1132 = value;
    } else if (name.compare("CLQ1_2232") == 0) {
        CLQ1_2232 = value;
    } else if (name.compare("CLQ1_3332") == 0) {
        CLQ1_3332 = value;
    } else if (name.compare("CLQ3") == 0) {
        CLQ3_1111 = value;
        CLQ3_1122 = value;
        CLQ3_2211 = value;
        CLQ3_1221 = value;
        CLQ3_2112 = value;
        CLQ3_1133 = value;
        CLQ3_3311 = value;
        CLQ3_1331 = value;
        CLQ3_3113 = value;
    } else if (name.compare("CLQ3_1111") == 0) {
        CLQ3_1111 = value;
    } else if (name.compare("CLQ3_1122") == 0) {
        CLQ3_1122 = value;
        CLQ3_2211 = value;
        CLQ3_2112 = value;
        CLQ3_1221 = value;
    } else if (name.compare("CLQ3_1133") == 0) {
        CLQ3_1133 = value;
        CLQ3_3311 = value;
        CLQ3_3113 = value;
        CLQ3_1331 = value;
    } else if (name.compare("CLQ3_1123") == 0) {
        CLQ3_1123 = value;
    } else if (name.compare("CLQ3_2223") == 0) {
        CLQ3_2223 = value;
    } else if (name.compare("CLQ3_3323") == 0) {
        CLQ3_3323 = value;
    } else if (name.compare("CLQ3_1132") == 0) {
        CLQ3_1132 = value;
    } else if (name.compare("CLQ3_2232") == 0) {
        CLQ3_2232 = value;
    } else if (name.compare("CLQ3_3332") == 0) {
        CLQ3_3332 = value;
    } else if (name.compare("Cee") == 0) {
        Cee_1111 = value;
        Cee_1122 = value;
        Cee_2211 = value;
        Cee_1133 = value;
        Cee_3311 = value;
    } else if (name.compare("Cee_1111") == 0) {
        Cee_1111 = value;
    } else if (name.compare("Cee_1122") == 0) {
        Cee_1122 = value;
        Cee_2211 = value;
    } else if (name.compare("Cee_1133") == 0) {
        Cee_1133 = value;
        Cee_3311 = value;
    } else if (name.compare("Ceu") == 0) {
        Ceu_1111 = value;
        Ceu_1122 = value;
        Ceu_2211 = value;
        Ceu_1133 = value;
        Ceu_3311 = value;
    } else if (name.compare("Ceu_1111") == 0) {
        Ceu_1111 = value;
    } else if (name.compare("Ceu_1122") == 0) {
        Ceu_1122 = value;
        Ceu_2211 = value;
    } else if (name.compare("Ceu_1133") == 0) {
        Ceu_1133 = value;
        Ceu_3311 = value;
    } else if (name.compare("Ced") == 0) {
        Ced_1111 = value;
        Ced_1122 = value;
        Ced_2211 = value;
        Ced_1133 = value;
        Ced_3311 = value;
    } else if (name.compare("Ced_1111") == 0) {
        Ced_1111 = value;
    } else if (name.compare("Ced_1122") == 0) {
        Ced_1122 = value;
        Ced_2211 = value;
    } else if (name.compare("Ced_1133") == 0) {
        Ced_1133 = value;
        Ced_3311 = value;
    } else if (name.compare("Ced_1123") == 0) {
        Ced_1123 = value;
    } else if (name.compare("Ced_2223") == 0) {
        Ced_2223 = value;
    } else if (name.compare("Ced_3323") == 0) {
        Ced_3323 = value;
    } else if (name.compare("Ced_1132") == 0) {
        Ced_1132 = value;
    } else if (name.compare("Ced_2232") == 0) {
        Ced_2232 = value;
    } else if (name.compare("Ced_3332") == 0) {
        Ced_3332 = value;
    } else if (name.compare("CLe") == 0) {
        CLe_1111 = value;
        CLe_1122 = value;
        CLe_2211 = value;
        CLe_1133 = value;
        CLe_3311 = value;
    } else if (name.compare("CLe_1111") == 0) {
        CLe_1111 = value;
    } else if (name.compare("CLe_1122") == 0) {
        CLe_1122 = value;
        CLe_2211 = value;
    } else if (name.compare("CLe_1133") == 0) {
        CLe_1133 = value;
        CLe_3311 = value;
    } else if (name.compare("CLu") == 0) {
        CLu_1111 = value;
        CLu_1122 = value;
        CLu_2211 = value;
        CLu_1133 = value;
        CLu_3311 = value;
    } else if (name.compare("CLu_1111") == 0) {
        CLu_1111 = value;
    } else if (name.compare("CLu_1122") == 0) {
        CLu_1122 = value;
        CLu_2211 = value;
    } else if (name.compare("CLu_1133") == 0) {
        CLu_1133 = value;
        CLu_3311 = value;
    } else if (name.compare("CLd") == 0) {
        CLd_1111 = value;
        CLd_1122 = value;
        CLd_2211 = value;
        CLd_1133 = value;
        CLd_3311 = value;
    } else if (name.compare("CLd_1111") == 0) {
        CLd_1111 = value;
    } else if (name.compare("CLd_1122") == 0) {
        CLd_1122 = value;
        CLd_2211 = value;
    } else if (name.compare("CLd_1133") == 0) {
        CLd_1133 = value;
        CLd_3311 = value;
    } else if (name.compare("CLd_1123") == 0) {
        CLd_1123 = value;
    } else if (name.compare("CLd_2223") == 0) {
        CLd_2223 = value;
    } else if (name.compare("CLd_3323") == 0) {
        CLd_3323 = value;
    } else if (name.compare("CLd_1132") == 0) {
        CLd_1132 = value;
    } else if (name.compare("CLd_2232") == 0) {
        CLd_2232 = value;
    } else if (name.compare("CLd_3332") == 0) {
        CLd_3332 = value;
    } else if (name.compare("CQe") == 0) {
        CQe_1111 = value;
        CQe_1122 = value;
        CQe_2211 = value;
        CQe_1133 = value;
        CQe_3311 = value;
    } else if (name.compare("CQe_1111") == 0) {
        CQe_1111 = value;
    } else if (name.compare("CQe_1122") == 0) {
        CQe_1122 = value;
        CQe_2211 = value;
    } else if (name.compare("CQe_1133") == 0) {
        CQe_1133 = value;
        CQe_3311 = value;
    } else if (name.compare("CQe_2311") == 0) {
        CQe_2311 = value;
    } else if (name.compare("CQe_2322") == 0) {
        CQe_2322 = value;
    } else if (name.compare("CQe_2333") == 0) {
        CQe_2333 = value;
    } else if (name.compare("CQe_3211") == 0) {
        CQe_3211 = value;
    } else if (name.compare("CQe_3222") == 0) {
        CQe_3222 = value;
    } else if (name.compare("CQe_3233") == 0) {
        CQe_3233 = value;
    } else if (name.compare("Lambda_NP") == 0) {
        Lambda_NP = value;
    } else if (name.compare("eVBF2_HZZ1") == 0) {
         eVBF2_HZZ1 = value;
    } else if (name.compare("eVBF2_HZZ2") == 0) {
         eVBF2_HZZ2 = value;
    } else if (name.compare("eVBF2_HZZ3") == 0) {
         eVBF2_HZZ3 = value;
    } else if (name.compare("eVBF2_HZA1") == 0) {
         eVBF2_HZA1 = value;
    } else if (name.compare("eVBF2_HZA2") == 0) {
         eVBF2_HZA2 = value;
    } else if (name.compare("eVBF2_HAA") == 0) {
         eVBF2_HAA = value;
    } else if (name.compare("eVBF2_HWW1") == 0) {
         eVBF2_HWW1 = value;
    } else if (name.compare("eVBF2_HWW2") == 0) {
         eVBF2_HWW2 = value;
    } else if (name.compare("eVBF2_HWW3") == 0) {
         eVBF2_HWW3 = value;
    } else if (name.compare("eVBF2_Hgg") == 0) {
         eVBF2_Hgg = value;
    } else if (name.compare("eVBF2_HZuL") == 0) {
         eVBF2_HZuL = value;
    } else if (name.compare("eVBF2_HZuR") == 0) {
         eVBF2_HZuR = value;
    } else if (name.compare("eVBF2_HZdL") == 0) {
         eVBF2_HZdL = value;
    } else if (name.compare("eVBF2_HZdR") == 0) {
         eVBF2_HZdR = value;
    } else if (name.compare("eVBF2_HWud") == 0) {
         eVBF2_HWud = value;
    } else if (name.compare("eVBF2_ZuL") == 0) {
         eVBF2_ZuL = value;
    } else if (name.compare("eVBF2_ZuR") == 0) {
         eVBF2_ZuR = value;
    } else if (name.compare("eVBF2_ZdL") == 0) {
         eVBF2_ZdL = value;
    } else if (name.compare("eVBF2_ZdR") == 0) {
         eVBF2_ZdR = value;
    } else if (name.compare("eVBF2_Wud") == 0) {
         eVBF2_Wud = value;
    } else if (name.compare("eVBF78_HZZ1") == 0) {
         eVBF78_HZZ1 = value;
    } else if (name.compare("eVBF78_HZZ2") == 0) {
         eVBF78_HZZ2 = value;
    } else if (name.compare("eVBF78_HZZ3") == 0) {
         eVBF78_HZZ3 = value;
    } else if (name.compare("eVBF78_HZA1") == 0) {
         eVBF78_HZA1 = value;
    } else if (name.compare("eVBF78_HZA2") == 0) {
         eVBF78_HZA2 = value;
    } else if (name.compare("eVBF78_HAA") == 0) {
         eVBF78_HAA = value;
    } else if (name.compare("eVBF78_HWW1") == 0) {
         eVBF78_HWW1 = value;
    } else if (name.compare("eVBF78_HWW2") == 0) {
         eVBF78_HWW2 = value;
    } else if (name.compare("eVBF78_HWW3") == 0) {
         eVBF78_HWW3 = value;
    } else if (name.compare("eVBF78_Hgg") == 0) {
         eVBF78_Hgg = value;
    } else if (name.compare("eVBF78_HZuL") == 0) {
         eVBF78_HZuL = value;
    } else if (name.compare("eVBF78_HZuR") == 0) {
         eVBF78_HZuR = value;
    } else if (name.compare("eVBF78_HZdL") == 0) {
         eVBF78_HZdL = value;
    } else if (name.compare("eVBF78_HZdR") == 0) {
         eVBF78_HZdR = value;
    } else if (name.compare("eVBF78_HWud") == 0) {
         eVBF78_HWud = value;
    } else if (name.compare("eVBF78_ZuL") == 0) {
         eVBF78_ZuL = value;
    } else if (name.compare("eVBF78_ZuR") == 0) {
         eVBF78_ZuR = value;
    } else if (name.compare("eVBF78_ZdL") == 0) {
         eVBF78_ZdL = value;
    } else if (name.compare("eVBF78_ZdR") == 0) {
         eVBF78_ZdR = value;
    } else if (name.compare("eVBF78_Wud") == 0) {
         eVBF78_Wud = value;
    } else if (name.compare("eWH2_HWW1") == 0) {
         eWH2_HWW1 = value;
    } else if (name.compare("eWH2_HWW2") == 0) {
         eWH2_HWW2 = value;
    } else if (name.compare("eWH2_HWW3") == 0) {
         eWH2_HWW3 = value;
    } else if (name.compare("eWH2_HWud") == 0) {
         eWH2_HWud = value;
    } else if (name.compare("eWH2_Wud") == 0) {
         eWH2_Wud = value;
    } else if (name.compare("eWH78_HWW1") == 0) {
         eWH78_HWW1 = value;
    } else if (name.compare("eWH78_HWW2") == 0) {
         eWH78_HWW2 = value;
    } else if (name.compare("eWH78_HWW3") == 0) {
         eWH78_HWW3 = value;
    } else if (name.compare("eWH78_HWud") == 0) {
         eWH78_HWud = value;
    } else if (name.compare("eWH78_Wud") == 0) {
         eWH78_Wud = value;
    } else if (name.compare("eZH2_HZZ1") == 0) {
         eZH2_HZZ1 = value;
    } else if (name.compare("eZH2_HZZ2") == 0) {
         eZH2_HZZ2 = value;
    } else if (name.compare("eZH2_HZZ3") == 0) {
         eZH2_HZZ3 = value;
    } else if (name.compare("eZH2_HZA1") == 0) {
         eZH2_HZA1 = value;
    } else if (name.compare("eZH2_HZA2") == 0) {
         eZH2_HZA2 = value;
    } else if (name.compare("eZH2_HZuL") == 0) {
         eZH2_HZuL = value;
    } else if (name.compare("eZH2_HZuR") == 0) {
         eZH2_HZuR = value;
    } else if (name.compare("eZH2_HZdL") == 0) {
         eZH2_HZdL = value;
    } else if (name.compare("eZH2_HZdR") == 0) {
         eZH2_HZdR = value;
    } else if (name.compare("eZH2_ZuL") == 0) {
         eZH2_ZuL = value;
    } else if (name.compare("eZH2_ZuR") == 0) {
         eZH2_ZuR = value;
    } else if (name.compare("eZH2_ZdL") == 0) {
         eZH2_ZdL = value;
    } else if (name.compare("eZH2_ZdR") == 0) {
         eZH2_ZdR = value;   
    } else if (name.compare("eZH78_HZZ1") == 0) {
         eZH78_HZZ1 = value;
    } else if (name.compare("eZH78_HZZ2") == 0) {
         eZH78_HZZ2 = value;
    } else if (name.compare("eZH78_HZZ3") == 0) {
         eZH78_HZZ3 = value;
    } else if (name.compare("eZH78_HZA1") == 0) {
         eZH78_HZA1 = value;
    } else if (name.compare("eZH78_HZA2") == 0) {
         eZH78_HZA2 = value;
    } else if (name.compare("eZH78_HZuL") == 0) {
         eZH78_HZuL = value;
    } else if (name.compare("eZH78_HZuR") == 0) {
         eZH78_HZuR = value;
    } else if (name.compare("eZH78_HZdL") == 0) {
         eZH78_HZdL = value;
    } else if (name.compare("eZH78_HZdR") == 0) {
         eZH78_HZdR = value;
    } else if (name.compare("eZH78_ZuL") == 0) {
         eZH78_ZuL = value;
    } else if (name.compare("eZH78_ZuR") == 0) {
         eZH78_ZuR = value;
    } else if (name.compare("eZH78_ZdL") == 0) {
         eZH78_ZdL = value;
    } else if (name.compare("eZH78_ZdR") == 0) {
         eZH78_ZdR = value;     
    } else if (name.compare("ettH2_Htt") == 0) {
         ettH2_Htt = value;
    } else if (name.compare("ettH2_Hgg") == 0) {
         ettH2_Hgg = value;
    } else if (name.compare("ettH78_Htt") == 0) {
         ettH78_Htt = value;
    } else if (name.compare("ettH78_Hgg") == 0) {
         ettH78_Hgg = value;
    } else if (name.compare("MwInput") == 0)
        MwInput = value;
    else
        NPbase::setParameter(name, value);
}

bool NPEffectiveGIMR::CheckParameters(const std::map<std::string, double>& DPars)
{
    if (FlagLeptonUniversal && FlagQuarkUniversal) {
        if (FlagMwInput) {
            if (DPars.find("MwInput") == DPars.end()) {
                std::cout << "ERROR: Missing mandatory NPEffectiveGIMR_LFU_QFU parameter MwInput" << std::endl;
                return false;
            }
        }
        for (int i = 0; i < NNPEffectiveGIMRVars_LFU_QFU; i++) {
            if (DPars.find(NPEffectiveGIMRVars_LFU_QFU[i]) == DPars.end()) {
                std::cout << "ERROR: Missing mandatory NPEffectiveGIMR_LFU_QFU parameter "
                        << NPEffectiveGIMRVars_LFU_QFU[i] << std::endl;
                return false;
            }
        }
        //} else if (FlagLeptonUniversal && !FlagQuarkUniversal) {
        //} else if (!FlagLeptonUniversal && FlagQuarkUniversal) {
    } else if (!FlagLeptonUniversal && !FlagQuarkUniversal) {
        if (FlagMwInput) {
            if (DPars.find("MwInput") == DPars.end()) {
                std::cout << "ERROR: Missing mandatory NPEffectiveGIMR parameter MwInput" << std::endl;
                return false;
            }
        }
        for (int i = 0; i < NNPEffectiveGIMRVars; i++) {
            if (DPars.find(NPEffectiveGIMRVars[i]) == DPars.end()) {
                std::cout << "ERROR: Missing mandatory NPEffectiveGIMR parameter"
                        << NPEffectiveGIMRVars[i] << std::endl;
                return false;
            }
        }
    } else
        throw std::runtime_error("Error in NPEffectiveGIMR::CheckParameters()");

    return (NPbase::CheckParameters(DPars));
}

bool NPEffectiveGIMR::setFlag(const std::string name, const bool value)
{
    bool res = false;
    if (name.compare("MwInput") == 0) {
        FlagMwInput = value;
        res = true;
    } else if (name.compare("QuadraticTerms") == 0) {
        FlagQuadraticTerms = value;
        if(value) setModelLinearized(false);
        res = true;
    } else
        res = NPbase::setFlag(name, value);

    return (res);
}


////////////////////////////////////////////////////////////////////////

double NPEffectiveGIMR::CHF1_diag(const Particle F) const
{
    if (F.is("NEUTRINO_1") || F.is("ELECTRON"))
        return CHL1_11;
    else if (F.is("NEUTRINO_2") || F.is("MU"))
        return CHL1_22;
    else if (F.is("NEUTRINO_3") || F.is("TAU"))
        return CHL1_33;
    else if (F.is("UP") || F.is("DOWN"))
        return CHQ1_11;
    else if (F.is("CHARM") || F.is("STRANGE"))
        return CHQ1_22;
    else if (F.is("TOP") || F.is("BOTTOM"))
        return CHQ1_33;
    else
        throw std::runtime_error("NPEffectiveGIMR::CHF1_diag(): wrong argument");
}

double NPEffectiveGIMR::CHF3_diag(const Particle F) const
{
    if (F.is("NEUTRINO_1") || F.is("ELECTRON"))
        return CHL3_11;
    else if (F.is("NEUTRINO_2") || F.is("MU"))
        return CHL3_22;
    else if (F.is("NEUTRINO_3") || F.is("TAU"))
        return CHL3_33;
    else if (F.is("UP") || F.is("DOWN"))
        return CHQ3_11;
    else if (F.is("CHARM") || F.is("STRANGE"))
        return CHQ3_22;
    else if (F.is("TOP") || F.is("BOTTOM"))
        return CHQ3_33;
    else
        throw std::runtime_error("NPEffectiveGIMR::CHF3_diag(): wrong argument");
}

double NPEffectiveGIMR::CHf_diag(const Particle f) const
{
    if (f.is("NEUTRINO_1") || f.is("NEUTRINO_2") || f.is("NEUTRINO_3"))
        return 0.0;
    else if (f.is("ELECTRON"))
        return CHe_11;
    else if (f.is("MU"))
        return CHe_22;
    else if (f.is("TAU"))
        return CHe_33;
    else if (f.is("UP"))
        return CHu_11;
    else if (f.is("CHARM"))
        return CHu_22;
    else if (f.is("TOP"))
        return CHu_33;
    else if (f.is("DOWN"))
        return CHd_11;
    else if (f.is("STRANGE"))
        return CHd_22;
    else if (f.is("BOTTOM"))
        return CHd_33;
    else
        throw std::runtime_error("NPEffectiveGIMR::CHf_diag(): wrong argument");
}

gslpp::complex NPEffectiveGIMR::CHud_diag(const Particle u) const
{
    if (!u.is("QUARK") || u.getIndex() % 2 != 0)
        throw std::runtime_error("NPEffectiveGIMR::CHud_diag(): wrong argument");

    if (u.is("UP"))
        return gslpp::complex(CHud_11r, CHud_11i, false);
    else if (u.is("CHARM"))
        return gslpp::complex(CHud_22r, CHud_22i, false);
    else if (u.is("TOP"))
        return gslpp::complex(CHud_22r, CHud_33i, false);
    else
        throw std::runtime_error("NPEffectiveGIMR::CHud_diag(): wrong argument");
}

gslpp::complex NPEffectiveGIMR::CfH_diag(const Particle f) const
{
    if (f.is("NEUTRINO_1") || f.is("NEUTRINO_2") || f.is("NEUTRINO_3"))
        return 0.0;
    else if (f.is("ELECTRON"))
        return gslpp::complex(CeH_11r, CeH_11i, false);
    else if (f.is("MU"))
        return gslpp::complex(CeH_22r, CeH_22i, false);
    else if (f.is("TAU"))
        return gslpp::complex(CeH_33r, CeH_33i, false);
    else if (f.is("UP"))
        return gslpp::complex(CuH_11r, CuH_11i, false);
    else if (f.is("CHARM"))
        return gslpp::complex(CuH_22r, CuH_22i, false);
    else if (f.is("TOP"))
        return gslpp::complex(CuH_33r, CuH_33i, false);
    else if (f.is("DOWN"))
        return gslpp::complex(CdH_11r, CdH_11i, false);
    else if (f.is("STRANGE"))
        return gslpp::complex(CdH_22r, CdH_22i, false);
    else if (f.is("BOTTOM"))
        return gslpp::complex(CdH_33r, CdH_33i, false);
    else
        throw std::runtime_error("NPEffectiveGIMR::CfH_diag(): wrong argument");
}

gslpp::complex NPEffectiveGIMR::CfG_diag(const Particle f) const
{
    if (f.is("NEUTRINO_1") || f.is("NEUTRINO_2") || f.is("NEUTRINO_3"))
        return 0.0;
    else if (f.is("ELECTRON"))
        return 0.0;
    else if (f.is("MU"))
        return 0.0;
    else if (f.is("TAU"))
        return 0.0;
    else if (f.is("UP"))
        return gslpp::complex(CuG_11r, CuG_11i, false);
    else if (f.is("CHARM"))
        return gslpp::complex(CuG_22r, CuG_22i, false);
    else if (f.is("TOP"))
        return gslpp::complex(CuG_33r, CuG_33i, false);
    else if (f.is("DOWN"))
        return 0.0;
    else if (f.is("STRANGE"))
        return 0.0;
    else if (f.is("BOTTOM"))
        return 0.0;
    else
        throw std::runtime_error("NPEffectiveGIMR::CfG_diag(): wrong argument");
}
      
gslpp::complex NPEffectiveGIMR::CfW_diag(const Particle f) const
{
    if (f.is("NEUTRINO_1") || f.is("NEUTRINO_2") || f.is("NEUTRINO_3"))
        return 0.0;
    else if (f.is("ELECTRON"))
        return 0.0;
    else if (f.is("MU"))
        return 0.0;
    else if (f.is("TAU"))
        return 0.0;
    else if (f.is("UP"))
        return gslpp::complex(CuW_11r, CuW_11i, false);
    else if (f.is("CHARM"))
        return gslpp::complex(CuW_22r, CuW_22i, false);
    else if (f.is("TOP"))
        return gslpp::complex(CuW_33r, CuW_33i, false);
    else if (f.is("DOWN"))
        return 0.0;
    else if (f.is("STRANGE"))
        return 0.0;
    else if (f.is("BOTTOM"))
        return 0.0;
    else
        throw std::runtime_error("NPEffectiveGIMR::CfW_diag(): wrong argument");
}
      
gslpp::complex NPEffectiveGIMR::CfB_diag(const Particle f) const
{
    if (f.is("NEUTRINO_1") || f.is("NEUTRINO_2") || f.is("NEUTRINO_3"))
        return 0.0;
    else if (f.is("ELECTRON"))
        return 0.0;
    else if (f.is("MU"))
        return 0.0;
    else if (f.is("TAU"))
        return 0.0;
    else if (f.is("UP"))
        return gslpp::complex(CuB_11r, CuB_11i, false);
    else if (f.is("CHARM"))
        return gslpp::complex(CuB_22r, CuB_22i, false);
    else if (f.is("TOP"))
        return gslpp::complex(CuB_33r, CuB_33i, false);
    else if (f.is("DOWN"))
        return 0.0;
    else if (f.is("STRANGE"))
        return 0.0;
    else if (f.is("BOTTOM"))
        return 0.0;
    else
        throw std::runtime_error("NPEffectiveGIMR::CfB_diag(): wrong argument");
}


////////////////////////////////////////////////////////////////////////

double NPEffectiveGIMR::DeltaGF() const
{
    return ((CHL3_11 + CHL3_22 - 0.5 * (CLL_1221 + CLL_2112)) * v2_over_LambdaNP2);
}

double NPEffectiveGIMR::obliqueS() const
{
    return (4.0 * sW_tree * cW_tree * CHWB / alphaMz() * v2_over_LambdaNP2);
}

double NPEffectiveGIMR::obliqueT() const
{
    return (-CHD / 2.0 / alphaMz() * v2_over_LambdaNP2);
}

double NPEffectiveGIMR::obliqueU() const
{
    return 0.0;
}

double NPEffectiveGIMR::deltaMz2() const
{
    return (2.0*(cW2_tree*CHW+sW2_tree*CHB+cW_tree*sW_tree*CHWB)*v2_over_LambdaNP2);
}

double NPEffectiveGIMR::Mw() const
{
    if (FlagMwInput)
        return MwInput;
    else
        return (trueSM.Mw() - Mw_tree() / 4.0 / (cW2_tree - sW2_tree)
            *(4.0 * sW_tree * cW_tree * CHWB * v2_over_LambdaNP2
            + cW2_tree * CHD * v2_over_LambdaNP2
            + 2.0 * sW2_tree * DeltaGF()));
}

double NPEffectiveGIMR::GammaW() const
{
    double G0 = GF * pow(Mw(), 3.0) / 6.0 / sqrt(2.0) / M_PI;
    double GammaW_tree = (3.0 + 2.0 * Nc) * G0;

    if (FlagMwInput)
        throw std::runtime_error("Write codes in NPEffectiveGIMR::GammaW()!");
    else
        return (trueSM.GammaW()
            - 3.0 * GammaW_tree / 4.0 / (cW2_tree - sW2_tree)
            *(4.0 * sW_tree * cW_tree * CHWB * v2_over_LambdaNP2
            + cW2_tree * CHD * v2_over_LambdaNP2
            + 2.0 * (1.0 + cW2_tree) / 3.0 * DeltaGF())
            + 2.0 * GammaW_tree / 3.0 * (CHL3_11 + CHQ3_11 + CHQ3_22) * v2_over_LambdaNP2);
}

double NPEffectiveGIMR::deltaGV_f(const Particle p) const
{
    return (deltaGL_f(p) + deltaGR_f(p));
}

double NPEffectiveGIMR::deltaGA_f(const Particle p) const
{
    return (deltaGL_f(p) - deltaGR_f(p));
}

double NPEffectiveGIMR::deltaGL_f(const Particle p) const
{
    double I3p = p.getIsospin(), Qp = p.getCharge();
    double CHF1 = CHF1_diag(p);
    double CHF3 = CHF3_diag(p);
    double NPindirect;
    if (FlagMwInput) {
        NPindirect = -I3p / 4.0 * (CHD * v2_over_LambdaNP2 + 2.0 * DeltaGF())
                + Qp * sW2_tree
                * ((cW_tree / sW_tree * CHWB + (1.0 + cW2_tree) / 4.0 / sW2_tree * CHD) * v2_over_LambdaNP2 + 0.5 * DeltaGF());
    } else {
        NPindirect = -I3p / 4.0 * (CHD * v2_over_LambdaNP2 + 2.0 * DeltaGF())
                - Qp * sW2_tree / 4.0 / (cW2_tree - sW2_tree)
                *((4.0 * cW_tree / sW_tree * CHWB + CHD) * v2_over_LambdaNP2 + 2.0 * DeltaGF());
    }
    double NPdirect = -0.5 * (CHF1 - 2.0 * I3p * CHF3) * v2_over_LambdaNP2;
    return (NPindirect + NPdirect);
}

double NPEffectiveGIMR::deltaGR_f(const Particle p) const
{
    double Qp = p.getCharge();
    double CHf = CHf_diag(p);
    double NPindirect;
    if (FlagMwInput) {
        NPindirect = Qp * sW2_tree
                * ((cW_tree / sW_tree * CHWB + (1.0 + cW2_tree) / 4.0 / sW2_tree * CHD) * v2_over_LambdaNP2 + 0.5 * DeltaGF());
    } else {
        NPindirect = -Qp * sW2_tree / 4.0 / (cW2_tree - sW2_tree)
                *((4.0 * cW_tree / sW_tree * CHWB + CHD) * v2_over_LambdaNP2 + 2.0 * DeltaGF());
    }
    double NPdirect = -0.5 * CHf*v2_over_LambdaNP2;
    return (NPindirect + NPdirect);
}


////////////////////////////////////////////////////////////////////////

gslpp::complex NPEffectiveGIMR::deltaGL_Wff(const Particle pbar, const Particle p) const
{
    if (pbar.getIndex() + 1 != p.getIndex() || pbar.getIndex() % 2 != 0)
        throw std::runtime_error("NPEffectiveGIMR::deltaGL_Wff(): Not implemented");

    double CHF3 = CHF3_diag(pbar);
    double NPindirect;
    if (FlagMwInput) {
        NPindirect = -0.5 * DeltaGF();
    } else {
        NPindirect = -cW2_tree / 4.0 / (cW2_tree - sW2_tree)
                * ((4.0 * sW_tree / cW_tree * CHWB + CHD) * v2_over_LambdaNP2 + 2.0 * DeltaGF());
    }
    double NPdirect = CHF3 * v2_over_LambdaNP2;
    return (NPindirect + NPdirect);
}

gslpp::complex NPEffectiveGIMR::deltaGR_Wff(const Particle pbar, const Particle p) const
{
    if (pbar.getIndex() + 1 != p.getIndex() || pbar.getIndex() % 2 != 0)
        throw std::runtime_error("NPEffectiveGIMR::deltaGR_Wff(): Not implemented");

    gslpp::complex CHud = CHud_diag(pbar);
    return (0.5 * CHud * v2_over_LambdaNP2);
}

double NPEffectiveGIMR::deltaG_hgg() const
{
    return (CHG * v2_over_LambdaNP2 / v());
}

double NPEffectiveGIMR::deltaG1_hWW() const
{
    return (2.0 * CHW * v2_over_LambdaNP2 / v());
}

double NPEffectiveGIMR::deltaG2_hWW() const
{
    return 0.0;
}

double NPEffectiveGIMR::deltaG3_hWW() const
{
    double NPindirect;
    if (FlagMwInput) {
        NPindirect = 2.0 * MwInput * MwInput / v() * (delta_h - 0.5 * DeltaGF());
    } else {
        NPindirect = 2.0 * cW2_tree * Mz * Mz / v()
                * (delta_h - 1.0 / 2.0 / (cW2_tree - sW2_tree)
                * ((4.0 * sW_tree * cW_tree * CHWB + cW2_tree * CHD) * v2_over_LambdaNP2 + DeltaGF()));
    }
    return NPindirect;
}

double NPEffectiveGIMR::deltaG1_hZZ() const
{
    return (delta_ZZ / v());
}

double NPEffectiveGIMR::deltaG2_hZZ() const
{
    return 0.0;
}

double NPEffectiveGIMR::deltaG3_hZZ() const
{
    double NPindirect = Mz * Mz / v() * (-0.5 * CHD * v2_over_LambdaNP2 + delta_h - 0.5 * DeltaGF());
    double NPdirect = Mz * Mz / v() * CHD * v2_over_LambdaNP2;
    return (NPindirect + NPdirect);
}

double NPEffectiveGIMR::deltaG1_hZA() const
{
    return (delta_AZ / v());
}

double NPEffectiveGIMR::deltaG2_hZA() const
{
    return 0.0;
}

double NPEffectiveGIMR::deltaG_hAA() const
{
    return (delta_AA / v());
}

gslpp::complex NPEffectiveGIMR::deltaG_hff(const Particle p) const
{
    /* The effects of the RG running are neglected. */
    double mf;
    if (p.is("TOP"))
        //mf = p.getMass(); // m_t(m_t)
        mf = mtpole; // pole mass
    else
        mf = p.getMass();
    gslpp::complex CfH = CfH_diag(p);
    return (-mf / v() * (delta_h - 0.5 * DeltaGF())
            + CfH * v2_over_LambdaNP2 / sqrt(2.0));
}

gslpp::complex NPEffectiveGIMR::deltaGL_Wffh(const Particle pbar, const Particle p) const
{
    if (pbar.getIndex() + 1 != p.getIndex() || pbar.getIndex() % 2 != 0)
        throw std::runtime_error("NPEffectiveGIMR::deltaGL_Wffh(): Not implemented");

    double CHF3 = CHF3_diag(pbar);
    return (2.0 * sqrt(2.0) * Mz * cW_tree / v() / v() * CHF3 * v2_over_LambdaNP2);
}

gslpp::complex NPEffectiveGIMR::deltaGR_Wffh(const Particle pbar, const Particle p) const
{
    if (pbar.getIndex() + 1 != p.getIndex() || pbar.getIndex() % 2 != 0)
        throw std::runtime_error("NPEffectiveGIMR::deltaGR_Wffh(): Not implemented");

    gslpp::complex CHud = CHud_diag(pbar);
    return (sqrt(2.0) * Mz * cW_tree / v() / v() * CHud * v2_over_LambdaNP2);
}

double NPEffectiveGIMR::deltaGL_Zffh(const Particle p) const
{
    double I3p = p.getIsospin();
    double CHF1 = CHF1_diag(p);
    double CHF3 = CHF3_diag(p);
    return (-2.0 * Mz / v() / v() * (CHF1 - 2.0 * I3p * CHF3) * v2_over_LambdaNP2);
}

double NPEffectiveGIMR::deltaGR_Zffh(const Particle p) const
{
    double CHf = CHf_diag(p);
    return (-2.0 * Mz / v() / v() * CHf * v2_over_LambdaNP2);
}

gslpp::complex NPEffectiveGIMR::deltaG_hGff(const Particle p) const
{
    /* Set to 0. for the moment */

    return 0.;
}
      
gslpp::complex NPEffectiveGIMR::deltaG_hZff(const Particle p) const
{
    /* Set to 0. for the moment */

    return 0.;
}
      
gslpp::complex NPEffectiveGIMR::deltaG_hAff(const Particle p) const
{
    /* Set to 0. for the moment */

    return 0.;
}
      
gslpp::complex NPEffectiveGIMR::deltaG_Gff(const Particle p) const
{
    /* Set to 0. for the moment */

    return 0.;
}
      
gslpp::complex NPEffectiveGIMR::deltaG_Zff(const Particle p) const
{
    /* Set to 0. for the moment */

    return 0.;
}
      
gslpp::complex NPEffectiveGIMR::deltaG_Aff(const Particle p) const
{
    /* Set to 0. for the moment */

    return 0.;
}
      
double NPEffectiveGIMR::deltag3G() const
{
    /* Set to 0. for the moment */

    return 0.;
}


////////////////////////////////////////////////////////////////////////

gslpp::complex NPEffectiveGIMR::f_triangle(const double tau) const
{
    gslpp::complex tmp;
    if (tau >= 1.0) {
        tmp = asin(1.0 / sqrt(tau));
        return (tmp * tmp);
    } else {
        tmp = log((1.0 + sqrt(1.0 - tau)) / (1.0 - sqrt(1.0 - tau))) - M_PI * gslpp::complex::i();
        return (-0.25 * tmp * tmp);
    }
}

gslpp::complex NPEffectiveGIMR::AH_f(const double tau) const
{
    return (2.0 * tau * (1.0 + (1.0 - tau) * f_triangle(tau)));
}

double NPEffectiveGIMR::muggH(const double sqrt_s) const
{
    double m_t = mtpole;
    //doulbe m_t = quarks[TOP].getMass();
    //double m_b = quarks[BOTTOM].getMass();

    gslpp::complex dKappa_t = deltaG_hff(quarks[TOP]) / (-m_t / v());
    //gslpp::complex dKappa_b = deltaG_hff(quarks[BOTTOM]) / (-m_b / v());

    /* L_eff = G_eff_t_SM*hGG */
    gslpp::complex G_eff_t_SM = AlsMz / 16.0 / M_PI / v() * AH_f(4.0 * m_t * m_t / mHl / mHl);

    //double sigma_tt_SM = trueSM.computeSigmaggH_tt(sqrt_s);
    //double sigma_bb_SM = trueSM.computeSigmaggH_bb(sqrt_s);
    //double sigma_tb_SM = trueSM.computeSigmaggH_tb(sqrt_s);
    //gslpp::complex tmp = (2.0 * dKappa_t * sigma_tt_SM
    //        + 2.0 * dKappa_b * sigma_bb_SM
    //        + (dKappa_t + dKappa_b) * sigma_tb_SM)
    //        / (sigma_tt_SM + sigma_bb_SM + sigma_tb_SM);
    gslpp::complex tmp = CHG / v() * v2_over_LambdaNP2 / G_eff_t_SM;
    
    double mu = (1.0 + 2.0 * ( dKappa_t.real() + tmp.real() ) );
    
    if (FlagQuadraticTerms) {
            //Add contributions that are quadratic in the effective coefficients
        gslpp::complex tmp2 = dKappa_t + tmp;
        
        mu += tmp2.abs2();
        
    }
    
    if (mu < 0) return std::numeric_limits<double>::quiet_NaN();
    
    return mu;
}

double NPEffectiveGIMR::muVBF(const double sqrt_s) const
{
    double mu = 1.0;
    if (sqrt_s == 1.96) {
        mu +=  +1.123 * (1. + eVBF2_ZuL ) * deltaGL_f(quarks[UP])
                -0.531 * (1. + eVBF2_ZuR ) * deltaGR_f(quarks[UP])
                -0.705 * (1. + eVBF2_ZdL ) * deltaGL_f(quarks[DOWN])
                +0.136 * (1. + eVBF2_ZdR ) * deltaGR_f(quarks[DOWN])
	        +2.662 * (1. + eVBF2_Wud ) * deltaGL_Wff(quarks[UP],quarks[DOWN]).real()
                -1407.72 * (1. + eVBF2_HWud ) * deltaGL_Wffh(quarks[UP], quarks[DOWN]).real()
                +14928.1 * (1. + eVBF2_Hgg ) * deltaG_hgg()
                -12.451 * (1. + eVBF2_HAA ) * deltaG_hAA()
                -21.274 * (1. + eVBF2_HZA1 ) * deltaG1_hZA()
                +45.617 * (1. + eVBF2_HZA2 ) * deltaG2_hZA()
                -84.016 * (1. + eVBF2_HWW1 ) * deltaG1_hWW()
                +390.524 * (1. + eVBF2_HWW2 ) * deltaG2_hWW()
                +0.026 * (1. + eVBF2_HWW3 ) * deltaG3_hWW()
                -45.832 * (1. + eVBF2_HZZ1 ) * deltaG1_hZZ()
                +88.358 * (1. + eVBF2_HZZ2 ) * deltaG2_hZZ()
                +0.012 * (1. + eVBF2_HZZ3 ) * deltaG3_hZZ()
                -129.338 * (1. + eVBF2_HZuL ) * deltaGL_Zffh(quarks[UP])
                +84.325 * (1. + eVBF2_HZuR ) * deltaGR_Zffh(quarks[UP])
                +164.195 * (1. + eVBF2_HZdL ) * deltaGL_Zffh(quarks[DOWN])
                -32.751 * (1. + eVBF2_HZdR ) * deltaGR_Zffh(quarks[DOWN]);
        
        if (FlagQuadraticTerms) {
            //Add contributions that are quadratic in the effective coefficients
            //(Only valid under the assumptions of one dim 6 operator at a time)
            mu +=  +2.478 * pow(deltaGL_f(quarks[UP]),2.0)
                +1.878 * pow(deltaGR_f(quarks[UP]),2.0)
                +1.214 * pow(deltaGL_f(quarks[DOWN]),2.0)
                +0.898 * pow(deltaGR_f(quarks[DOWN]),2.0)
	        +2.659 * pow(deltaGL_Wff(quarks[UP],quarks[DOWN]).real(),2.0)
                +1917816. * pow(deltaGL_Wffh(quarks[UP], quarks[DOWN]).real(),2.0)
                +524312994. * pow(deltaG_hgg(),2.0)
                +831253. * pow(deltaG_hAA(),2.0)
                +151140. * pow(deltaG1_hZA(),2.0)
                +58067.7 * pow(deltaG2_hZA(),2.0)
                +106835. * pow(deltaG1_hWW(),2.0)
                +219369. * pow(deltaG2_hWW(),2.0)
                +145840. * pow(deltaG1_hZZ(),2.0)
                +66461.2 * pow(deltaG2_hZZ(),2.0)
                +1608277. * pow(deltaGL_Zffh(quarks[UP]),2.0)
                +1449825. * pow(deltaGR_Zffh(quarks[UP]),2.0)
                +409700. * pow(deltaGL_Zffh(quarks[DOWN]),2.0)
                +385965. * pow(deltaGR_Zffh(quarks[DOWN]),2.0);            
        }
        
    } else if (sqrt_s == 7.0) {
        mu +=  +1.188 * (1. + eVBF78_ZuL ) * deltaGL_f(quarks[UP])
                -0.536 * (1. + eVBF78_ZuR ) * deltaGR_f(quarks[UP])
                -0.976 * (1. + eVBF78_ZdL ) * deltaGL_f(quarks[DOWN])
                +0.179 * (1. + eVBF78_ZdR ) * deltaGR_f(quarks[DOWN])
	        +2.592 * (1. + eVBF78_Wud ) * deltaGL_Wff(quarks[UP],quarks[DOWN]).real()
                -1826.63 * (1. + eVBF78_HWud ) * deltaGL_Wffh(quarks[UP], quarks[DOWN]).real()
                +14265.8 * (1. + eVBF78_Hgg ) * deltaG_hgg()
                -40.051 * (1. + eVBF78_HAA ) * deltaG_hAA()
                -42.43 * (1. + eVBF78_HZA1 ) * deltaG1_hZA()
                +88.972 * (1. + eVBF78_HZA2 ) * deltaG2_hZA()
                -108.107 * (1. + eVBF78_HWW1 ) * deltaG1_hWW()
                +547.508 * (1. + eVBF78_HWW2 ) * deltaG2_hWW()
                +0.026 * (1. + eVBF78_HWW3 ) * deltaG3_hWW()
                -67.672 * (1. + eVBF78_HZZ1 ) * deltaG1_hZZ()
                +168.86 * (1. + eVBF78_HZZ2 ) * deltaG2_hZZ()
                +0.014 * (1. + eVBF78_HZZ3 ) * deltaG3_hZZ()
                -466.198 * (1. + eVBF78_HZuL ) * deltaGL_Zffh(quarks[UP])
                +211.308 * (1. + eVBF78_HZuR ) * deltaGR_Zffh(quarks[UP])
                +374.597 * (1. + eVBF78_HZdL ) * deltaGL_Zffh(quarks[DOWN])
                -69.916 * (1. + eVBF78_HZdR ) * deltaGR_Zffh(quarks[DOWN]);
        
        if (FlagQuadraticTerms) {
            //Add contributions that are quadratic in the effective coefficients
            //(Only valid under the assumptions of one dim 6 operator at a time)
            mu +=  +2.534 * pow(deltaGL_f(quarks[UP]),2.0)
                +1.9 * pow(deltaGR_f(quarks[UP]),2.0)
                +1.695 * pow(deltaGL_f(quarks[DOWN]),2.0)
                +1.177 * pow(deltaGR_f(quarks[DOWN]),2.0)
	        +2.608 * pow(deltaGL_Wff(quarks[UP],quarks[DOWN]).real(),2.0)
                +2862580. * pow(deltaGL_Wffh(quarks[UP], quarks[DOWN]).real(),2.0)
                +519301209. * pow(deltaG_hgg(),2.0)
                +777159. * pow(deltaG_hAA(),2.0)
                +206157. * pow(deltaG1_hZA(),2.0)
                +94511.2 * pow(deltaG2_hZA(),2.0)
                +174828. * pow(deltaG1_hWW(),2.0)
                +414624. * pow(deltaG2_hWW(),2.0)
                +209132. * pow(deltaG1_hZZ(),2.0)
                +120250. * pow(deltaG2_hZZ(),2.0)
                +1311032. * pow(deltaGL_Zffh(quarks[UP]),2.0)
                +1130789. * pow(deltaGR_Zffh(quarks[UP]),2.0)
                +757088. * pow(deltaGL_Zffh(quarks[DOWN]),2.0)
                +651756. * pow(deltaGR_Zffh(quarks[DOWN]),2.0);            
        }
        
    } else if (sqrt_s == 8.0) {
        mu +=  +1.179 * (1. + eVBF78_ZuL ) * deltaGL_f(quarks[UP])
                -0.532 * (1. + eVBF78_ZuR ) * deltaGR_f(quarks[UP])
                -0.984 * (1. + eVBF78_ZdL ) * deltaGL_f(quarks[DOWN])
                +0.181 * (1. + eVBF78_ZdR ) * deltaGR_f(quarks[DOWN])
	            +2.591 * (1. + eVBF78_Wud ) * deltaGL_Wff(quarks[UP],quarks[DOWN]).real()
                -1858.03 * (1. + eVBF78_HWud ) * deltaGL_Wffh(quarks[UP], quarks[DOWN]).real()
                +14247.4 * (1. + eVBF78_Hgg ) * deltaG_hgg()
                -40.46 * (1. + eVBF78_HAA ) * deltaG_hAA()
                -41.713 * (1. + eVBF78_HZA1 ) * deltaG1_hZA()
                +90.462 * (1. + eVBF78_HZA2 ) * deltaG2_hZA()
                -106.576 * (1. + eVBF78_HWW1 ) * deltaG1_hWW()
                +562.98 * (1. + eVBF78_HWW2 ) * deltaG2_hWW()
                +0.026 * (1. + eVBF78_HWW3 ) * deltaG3_hWW()
                -67.57 * (1. + eVBF78_HZZ1 ) * deltaG1_hZZ()
                +174.474 * (1. + eVBF78_HZZ2 ) * deltaG2_hZZ()
                +0.014 * (1. + eVBF78_HZZ3 ) * deltaG3_hZZ()
                -472.887 * (1. + eVBF78_HZuL ) * deltaGL_Zffh(quarks[UP])
                +214.739 * (1. + eVBF78_HZuR ) * deltaGR_Zffh(quarks[UP])
                +386.582 * (1. + eVBF78_HZdL ) * deltaGL_Zffh(quarks[DOWN])
                -72.228 * (1. + eVBF78_HZdR ) * deltaGR_Zffh(quarks[DOWN]);
        
        if (FlagQuadraticTerms) {
            //Add contributions that are quadratic in the effective coefficients
            //(Only valid under the assumptions of one dim 6 operator at a time)
            mu +=  +2.503 * pow(deltaGL_f(quarks[UP]),2.0)
                +1.877 * pow(deltaGR_f(quarks[UP]),2.0)
                +1.712 * pow(deltaGL_f(quarks[DOWN]),2.0)
                +1.191 * pow(deltaGR_f(quarks[DOWN]),2.0)
	        +2.606 * pow(deltaGL_Wff(quarks[UP],quarks[DOWN]).real(),2.0)
                +3057041. * pow(deltaGL_Wffh(quarks[UP], quarks[DOWN]).real(),2.0)
                +517064803. * pow(deltaG_hgg(),2.0)
                +766750. * pow(deltaG_hAA(),2.0)
                +207500. * pow(deltaG1_hZA(),2.0)
                +101779. * pow(deltaG2_hZA(),2.0)
                +177714. * pow(deltaG1_hWW(),2.0)
                +454117. * pow(deltaG2_hWW(),2.0)
                +210212. * pow(deltaG1_hZZ(),2.0)
                +131594. * pow(deltaG2_hZZ(),2.0)
                +1399281. * pow(deltaGL_Zffh(quarks[UP]),2.0)
                +1231240. * pow(deltaGR_Zffh(quarks[UP]),2.0)
                +820259. * pow(deltaGL_Zffh(quarks[DOWN]),2.0)
                +713820. * pow(deltaGR_Zffh(quarks[DOWN]),2.0);            
        }
        
    } else if (sqrt_s == 14.0) {
        mu +=  +1.129 * deltaGL_f(quarks[UP])
                -0.505 * deltaGR_f(quarks[UP])
                -1.05 * deltaGL_f(quarks[DOWN])
                +0.191 * deltaGR_f(quarks[DOWN])
	        +2.586 * deltaGL_Wff(quarks[UP],quarks[DOWN]).real()
                -1989.34 * deltaGL_Wffh(quarks[UP], quarks[DOWN]).real()
                +14228.8 * deltaG_hgg()
                -35.554 * deltaG_hAA()
                -39.847 * deltaG1_hZA()
                +98.522 * deltaG2_hZA()
                -99.287 * deltaG1_hWW()
                +622.352 * deltaG2_hWW()
                +0.026 * deltaG3_hWW()
                -66.196 * deltaG1_hZZ()
                +196.676 * deltaG2_hZZ()
                +0.014 * deltaG3_hZZ()
                -493.198 * deltaGL_Zffh(quarks[UP])
                +217.017 * deltaGR_Zffh(quarks[UP])
                +447.396 * deltaGL_Zffh(quarks[DOWN])
                -82.396 * deltaGR_Zffh(quarks[DOWN]); 
        
        if (FlagQuadraticTerms) {
            //Add contributions that are quadratic in the effective coefficients
            //(Only valid under the assumptions of one dim 6 operator at a time)
            mu +=  +2.319 * pow(deltaGL_f(quarks[UP]),2.0)
                +1.783 * pow(deltaGR_f(quarks[UP]),2.0)
                +1.849 * pow(deltaGL_f(quarks[DOWN]),2.0)
                +1.263 * pow(deltaGR_f(quarks[DOWN]),2.0)
	        +2.592 * pow(deltaGL_Wff(quarks[UP],quarks[DOWN]).real(),2.0)
                +4077238. * pow(deltaGL_Wffh(quarks[UP], quarks[DOWN]).real(),2.0)
                +507787376. * pow(deltaG_hgg(),2.0)
                +702353. * pow(deltaG_hAA(),2.0)
                +212082. * pow(deltaG1_hZA(),2.0)
                +141422. * pow(deltaG2_hZA(),2.0)
                +195770. * pow(deltaG1_hWW(),2.0)
                +655804. * pow(deltaG2_hWW(),2.0)
                +0. * pow(deltaG3_hWW(),2.0)
                +240333. * pow(deltaG1_hZZ(),2.0)
                +192371. * pow(deltaG2_hZZ(),2.0)
                +0. * pow(deltaG3_hZZ(),2.0)
                +1904757. * pow(deltaGL_Zffh(quarks[UP]),2.0)
                +1743849. * pow(deltaGR_Zffh(quarks[UP]),2.0)
                +1185212. * pow(deltaGL_Zffh(quarks[DOWN]),2.0)
                +1061826. * pow(deltaGR_Zffh(quarks[DOWN]),2.0); 
        }
        
    } else
        throw std::runtime_error("Bad argument in NPEffectiveGIMR::muVBF()");

    if (mu < 0) return std::numeric_limits<double>::quiet_NaN();
    
    return mu;
}

double NPEffectiveGIMR::mueeWBF(const double sqrt_s) const
{
    double mu = 1.0;
    if (sqrt_s == 0.24) {
        mu +=  985.974 * deltaG1_hZZ()
                +23.622 * deltaG2_hZZ()
                +0. * deltaG3_hZZ()
                +377.441 * deltaG1_hZA()
                -394.144 * deltaG2_hZA()
                -30.997 * deltaG1_hWW()
                +273.526 * deltaG2_hWW()
                +0.038 * deltaG3_hWW()
                -62.505 * deltaGL_Zffh(leptons[NEUTRINO_1])
                -1781.02 * deltaGL_Zffh(leptons[ELECTRON])
                +183.495 * deltaGR_Zffh(leptons[ELECTRON])
                -4.991 * deltaGL_Wffh(leptons[NEUTRINO_1],leptons[ELECTRON]).real()
                -0.062 * deltaGL_f(leptons[NEUTRINO_1])
                -1.436 * deltaGL_f(leptons[ELECTRON])
                +0.004 * deltaGR_f(leptons[ELECTRON])
                +4.07 * deltaGL_Wff(leptons[NEUTRINO_1],leptons[ELECTRON]).real();
        
        if (FlagQuadraticTerms) {
            //Add contributions that are quadratic in the effective coefficients
            //(Only valid under the assumptions of one dim 6 operator at a time)
            mu +=  +0.0;            
        }
          
    } else if (sqrt_s == 0.25) {
        mu +=  903.947 * deltaG1_hZZ()
                +93.416 * deltaG2_hZZ()
                -0.002 * deltaG3_hZZ()
                +277.754 * deltaG1_hZA()
                -271.135 * deltaG2_hZA()
                -29.647 * deltaG1_hWW()
                +298.034 * deltaG2_hWW()
                +0.039 * deltaG3_hWW()
                -61.73 * deltaGL_Zffh(leptons[NEUTRINO_1])
                -1409.1 * deltaGL_Zffh(leptons[ELECTRON])
                +182.068 * deltaGR_Zffh(leptons[ELECTRON])
                -5.518 * deltaGL_Wffh(leptons[NEUTRINO_1],leptons[ELECTRON]).real()
                -0.174 * deltaGL_f(leptons[NEUTRINO_1])
                -1.003 * deltaGL_f(leptons[ELECTRON])
                +0. * deltaGR_f(leptons[ELECTRON])
                +4.203 * deltaGL_Wff(leptons[NEUTRINO_1],leptons[ELECTRON]).real();
        
        if (FlagQuadraticTerms) {
            //Add contributions that are quadratic in the effective coefficients
            //(Only valid under the assumptions of one dim 6 operator at a time)
            mu +=  +0.0;
        }
        
    } else if (sqrt_s == 0.35) {
        mu +=  -63.056 * deltaG1_hZZ()
                +99.834 * deltaG2_hZZ()
                -0.001 * deltaG3_hZZ()
                -100.006 * deltaG1_hZA()
                +120.858 * deltaG2_hZA()
                -77.085 * deltaG1_hWW()
                +424.21 * deltaG2_hWW()
                +0.039 * deltaG3_hWW()
                -36.176 * deltaGL_Zffh(leptons[NEUTRINO_1])
                +108.441 * deltaGL_Zffh(leptons[ELECTRON])
                +137.702 * deltaGR_Zffh(leptons[ELECTRON])
                -8.003 * deltaGL_Wffh(leptons[NEUTRINO_1],leptons[ELECTRON]).real()
                -0.091 * deltaGL_f(leptons[NEUTRINO_1])
                +0.074 * deltaGL_f(leptons[ELECTRON])
                +0.033 * deltaGR_f(leptons[ELECTRON])
                +4.102 * deltaGL_Wff(leptons[NEUTRINO_1],leptons[ELECTRON]).real();
        
        if (FlagQuadraticTerms) {
            //Add contributions that are quadratic in the effective coefficients
            //(Only valid under the assumptions of one dim 6 operator at a time)
            mu +=  +0.0;
        }
        
    } else if (sqrt_s == 0.5) {
        mu +=  -82.771 * deltaG1_hZZ()
                +48.73 * deltaG2_hZZ()
                +0. * deltaG3_hZZ()
                -78.056 * deltaG1_hZA()
                +78.264 * deltaG2_hZA()
                -98.794 * deltaG1_hWW()
                +579.5 * deltaG2_hWW()
                +0.039 * deltaG3_hWW()
                -26.448 * deltaGL_Zffh(leptons[NEUTRINO_1])
                +163.236 * deltaGL_Zffh(leptons[ELECTRON])
                +56.583 * deltaGR_Zffh(leptons[ELECTRON])
                -9.491 * deltaGL_Wffh(leptons[NEUTRINO_1],leptons[ELECTRON]).real()
                -0.02 * deltaGL_f(leptons[NEUTRINO_1])
                +0.037 * deltaGL_f(leptons[ELECTRON])
                +0.009 * deltaGR_f(leptons[ELECTRON])
                +4.025 * deltaGL_Wff(leptons[NEUTRINO_1],leptons[ELECTRON]).real();
        
        if (FlagQuadraticTerms) {
            //Add contributions that are quadratic in the effective coefficients
            //(Only valid under the assumptions of one dim 6 operator at a time)
            mu +=  +0.0;
        }
        
    } else if (sqrt_s == 1.0) {
        mu +=  -32.198 * deltaG1_hZZ()
                +13.389 * deltaG2_hZZ()
                +0. * deltaG3_hZZ()
                -27.018 * deltaG1_hZA()
                +18.957 * deltaG2_hZA()
                -100.42 * deltaG1_hWW()
                +884.402 * deltaG2_hWW()
                +0.039 * deltaG3_hWW()
                -11.556 * deltaGL_Zffh(leptons[NEUTRINO_1])
                +56.187 * deltaGL_Zffh(leptons[ELECTRON])
                +16.258 * deltaGR_Zffh(leptons[ELECTRON])
                -12.168 * deltaGL_Wffh(leptons[NEUTRINO_1],leptons[ELECTRON]).real()
                -0.001 * deltaGL_f(leptons[NEUTRINO_1])
                +0.004 * deltaGL_f(leptons[ELECTRON])
                +0.002 * deltaGR_f(leptons[ELECTRON])
                +4.003 * deltaGL_Wff(leptons[NEUTRINO_1],leptons[ELECTRON]).real();
        
        if (FlagQuadraticTerms) {
            //Add contributions that are quadratic in the effective coefficients
            //(Only valid under the assumptions of one dim 6 operator at a time)
            mu +=  +0.0;
        }
        
    } else
        throw std::runtime_error("Bad argument in NPEffectiveGIMR::mueeWBF()");

    if (mu < 0) return std::numeric_limits<double>::quiet_NaN();
    
    return mu;
}

double NPEffectiveGIMR::muWH(const double sqrt_s) const
{
    double mu = 1.0;
    if (sqrt_s == 1.96) {
        mu += +2.032 * (1. + eWH2_Wud ) * deltaGL_Wff(quarks[UP], quarks[DOWN]).real()
                +1738.87 * (1. + eWH2_HWW1 ) * deltaG1_hWW()
                -3432.64 * (1. + eWH2_HWW2 ) * deltaG2_hWW()
                +0.039 * (1. + eWH2_HWW3 ) * deltaG3_hWW()
                +6523.35 * (1. + eWH2_HWud ) * deltaGL_Wffh(quarks[UP], quarks[DOWN]).real();
        
        if (FlagQuadraticTerms) {
            //Add contributions that are quadratic in the effective coefficients
            //(Only valid under the assumptions of one dim 6 operator at a time)
            mu += +1.042 * pow(deltaGL_Wff(quarks[UP], quarks[DOWN]).real(),2.0)
                +1075949.  * pow(deltaG1_hWW(),2.0)
                +3978950.  * pow(deltaG2_hWW(),2.0)
                +15887131. * pow(deltaGL_Wffh(quarks[UP], quarks[DOWN]).real(),2.0);
        }
        
    } else if (sqrt_s == 7.0) {
        mu += +1.979 * (1. + eWH78_Wud ) * deltaGL_Wff(quarks[UP], quarks[DOWN]).real()
                +1777.77 * (1. + eWH78_HWW1 ) * deltaG1_hWW()
                -3890.65 * (1. + eWH78_HWW2 ) * deltaG2_hWW()
                +0.039 * (1. + eWH78_HWW3 ) * deltaG3_hWW()
                +7344.73 * (1. + eWH78_HWud ) * deltaGL_Wffh(quarks[UP], quarks[DOWN]).real();
        
        if (FlagQuadraticTerms) {
            //Add contributions that are quadratic in the effective coefficients
            //(Only valid under the assumptions of one dim 6 operator at a time)
        mu += +1.015 * pow(deltaGL_Wff(quarks[UP], quarks[DOWN]).real(),2.0)
                +1294405.  * pow(deltaG1_hWW(),2.0)
                +7356224.  * pow(deltaG2_hWW(),2.0)
                +31355627. * pow(deltaGL_Wffh(quarks[UP], quarks[DOWN]).real(),2.0);
        }
        
    } else if (sqrt_s == 8.0) {
        mu += +1.978 * (1. + eWH78_Wud ) * deltaGL_Wff(quarks[UP], quarks[DOWN]).real()
                +1784.47 * (1. + eWH78_HWW1 ) * deltaG1_hWW()
                -3967.38 * (1. + eWH78_HWW2 ) * deltaG2_hWW()
                +0.039 * (1. + eWH78_HWW3 ) * deltaG3_hWW()
                +7507.02 * (1. + eWH78_HWud ) * deltaGL_Wffh(quarks[UP], quarks[DOWN]).real();
        
        if (FlagQuadraticTerms) {
            //Add contributions that are quadratic in the effective coefficients
            //(Only valid under the assumptions of one dim 6 operator at a time)
            mu += +1.016 * pow(deltaGL_Wff(quarks[UP], quarks[DOWN]).real(),2.0)
                +1331512.  * pow(deltaG1_hWW(),2.0)
                +8168916.  * pow(deltaG2_hWW(),2.0)
                +35201222. * pow(deltaGL_Wffh(quarks[UP], quarks[DOWN]).real(),2.0);
        }
        
    } else if (sqrt_s == 14.0) {
        mu += +1.963 * deltaGL_Wff(quarks[UP], quarks[DOWN]).real()
                +1799.45 * deltaG1_hWW()
                -4252.03 * deltaG2_hWW()
                +0.039 * deltaG3_hWW()
                +8047.59 * deltaGL_Wffh(quarks[UP], quarks[DOWN]).real(); 
        
        if (FlagQuadraticTerms) {
            //Add contributions that are quadratic in the effective coefficients
            //(Only valid under the assumptions of one dim 6 operator at a time)
            mu +=  +1.007 * pow(deltaGL_Wff(quarks[UP], quarks[DOWN]).real(),2.0)
                +1467903.  * pow(deltaG1_hWW(),2.0)
                +13173439.  * pow(deltaG2_hWW(),2.0)
                +58780336.  * pow(deltaGL_Wffh(quarks[UP], quarks[DOWN]).real(),2.0); 
        }
        
    } else
        throw std::runtime_error("Bad argument in NPEffectiveGIMR::muWH()");

    if (mu < 0) return std::numeric_limits<double>::quiet_NaN();
    
    return mu;
}

double NPEffectiveGIMR::muZH(const double sqrt_s) const
{
    double mu = 1.0;
    if (sqrt_s == 1.96) {
        mu += +3.529 * (1. + eZH2_ZuL ) * deltaGL_f(quarks[UP])
                -1.598 * (1. + eZH2_ZuR ) * deltaGR_f(quarks[UP])
                -1.229 * (1. + eZH2_ZdL ) * deltaGL_f(quarks[DOWN])
                +0.227 * (1. + eZH2_ZdR ) * deltaGR_f(quarks[DOWN])
                +3215.38 * (1. + eZH2_HZZ1 ) * deltaG1_hZZ()
                -2922.42 * (1. + eZH2_HZZ2 ) * deltaG2_hZZ()
                +0.059 * (1. + eZH2_HZZ3 ) * deltaG3_hZZ()
                +495.399 * (1. + eZH2_HZA1 ) * deltaG1_hZA()
                -838.743 * (1. + eZH2_HZA2 ) * deltaG2_hZA()
                +5931.99 * (1. + eZH2_HZuL ) * deltaGL_Zffh(quarks[UP])
                -2684.23 * (1. + eZH2_HZuR ) * deltaGR_Zffh(quarks[UP])
                -1878.46 * (1. + eZH2_HZdL ) * deltaGL_Zffh(quarks[DOWN])
                +346.694 * (1. + eZH2_HZdR ) * deltaGR_Zffh(quarks[DOWN]);
        
        if (FlagQuadraticTerms) {
            //Add contributions that are quadratic in the effective coefficients
            //(Only valid under the assumptions of one dim 6 operator at a time)
            mu += +5.126 * pow(deltaGL_f(quarks[UP]),2.0)
                +5.126 * pow(deltaGR_f(quarks[UP]),2.0)
                +1.456 * pow(deltaGL_f(quarks[DOWN]),2.0)
                +1.454 * pow(deltaGR_f(quarks[DOWN]),2.0)
                +3525123.  * pow(deltaG1_hZZ(),2.0)
                +2844179.  * pow(deltaG2_hZZ(),2.0)
                +0.001  * pow(deltaG3_hZZ(),2.0)
                +662397.  * pow(deltaG1_hZA(),2.0)
                +2006248.  * pow(deltaG2_hZA(),2.0)
                +21799545.  * pow(deltaGL_Zffh(quarks[UP]),2.0)
                +21795795.  * pow(deltaGR_Zffh(quarks[UP]),2.0)
                +4723149.  * pow(deltaGL_Zffh(quarks[DOWN]),2.0)
                +4725123.  * pow(deltaGR_Zffh(quarks[DOWN]),2.0);
        }
        
    } else if (sqrt_s == 7.0) {
        mu += +2.583 * (1. + eZH78_ZuL ) * deltaGL_f(quarks[UP])
                -1.17 * (1. + eZH78_ZuR ) * deltaGR_f(quarks[UP])
                -2.127 * (1. + eZH78_ZdL ) * deltaGL_f(quarks[DOWN])
                +0.392 * (1. + eZH78_ZdR ) * deltaGR_f(quarks[DOWN])
                +3269.53 * (1. + eZH78_HZZ1 ) * deltaG1_hZZ()
                -3201.65 * (1. + eZH78_HZZ2 ) * deltaG2_hZZ()
                +0.059 * (1. + eZH78_HZZ3 ) * deltaG3_hZZ()
                +473.267 * (1. + eZH78_HZA1 ) * deltaG1_hZA()
                -873.421 * (1. + eZH78_HZA2 ) * deltaG2_hZA()
                +4763.44 * (1. + eZH78_HZuL ) * deltaGL_Zffh(quarks[UP])
                -2156.99 * (1. + eZH78_HZuR ) * deltaGR_Zffh(quarks[UP])
                -3853.2 * (1. + eZH78_HZdL ) * deltaGL_Zffh(quarks[DOWN])
                +712.124 * (1. + eZH78_HZdR ) * deltaGR_Zffh(quarks[DOWN]);
        
        if (FlagQuadraticTerms) {
            //Add contributions that are quadratic in the effective coefficients
            //(Only valid under the assumptions of one dim 6 operator at a time)
            mu +=  +3.752 * pow(deltaGL_f(quarks[UP]),2.0)
                +3.753 * pow(deltaGR_f(quarks[UP]),2.0)
                +2.519 * pow(deltaGL_f(quarks[DOWN]),2.0)
                +2.517 * pow(deltaGR_f(quarks[DOWN]),2.0)
                +4051505.  * pow(deltaG1_hZZ(),2.0)
                +4597749.  * pow(deltaG2_hZZ(),2.0)
                +0.001  * pow(deltaG3_hZZ(),2.0)
                +610510.  * pow(deltaG1_hZA(),2.0)
                +2766996.  * pow(deltaG2_hZA(),2.0)
                +27425400.  * pow(deltaGL_Zffh(quarks[UP]),2.0)
                +27416894.  * pow(deltaGR_Zffh(quarks[UP]),2.0)
                +17043782.  * pow(deltaGL_Zffh(quarks[DOWN]),2.0)
                +17039528.  * pow(deltaGR_Zffh(quarks[DOWN]),2.0);
        }
        
    } else if (sqrt_s == 8.0) {
        mu += +2.569 * (1. + eZH78_ZuL ) * deltaGL_f(quarks[UP])
                -1.163 * (1. + eZH78_ZuR ) * deltaGR_f(quarks[UP])
                -2.14 * (1. + eZH78_ZdL ) * deltaGL_f(quarks[DOWN])
                +0.395 * (1. + eZH78_ZdR ) * deltaGR_f(quarks[DOWN])
                +3282.79 * (1. + eZH78_HZZ1 ) * deltaG1_hZZ()
                -3262.46 * (1. + eZH78_HZZ2 ) * deltaG2_hZZ()
                +0.059 * (1. + eZH78_HZZ3 ) * deltaG3_hZZ()
                +475.044 * (1. + eZH78_HZA1 ) * deltaG1_hZA()
                -892.243 * (1. + eZH78_HZA2 ) * deltaG2_hZA()
                +4847.78 * (1. + eZH78_HZuL ) * deltaGL_Zffh(quarks[UP])
                -2193.61 * (1. + eZH78_HZuR ) * deltaGR_Zffh(quarks[UP])
                -3960.46 * (1. + eZH78_HZdL ) * deltaGL_Zffh(quarks[DOWN])
                +731.438 * (1. + eZH78_HZdR ) * deltaGR_Zffh(quarks[DOWN]);
        
        if (FlagQuadraticTerms) {
            //Add contributions that are quadratic in the effective coefficients
            //(Only valid under the assumptions of one dim 6 operator at a time)
            mu += +3.732 * pow(deltaGL_f(quarks[UP]),2.0)
                +3.736 * pow(deltaGR_f(quarks[UP]),2.0)
                +2.535 * pow(deltaGL_f(quarks[DOWN]),2.0)
                +2.536 * pow(deltaGR_f(quarks[DOWN]),2.0)
                +4164701.  * pow(deltaG1_hZZ(),2.0)
                +5067698.  * pow(deltaG2_hZZ(),2.0)
                +0.001  * pow(deltaG3_hZZ(),2.0)
                +627966.  * pow(deltaG1_hZA(),2.0)
                +3087745.  * pow(deltaG2_hZA(),2.0)
                +30566228.  * pow(deltaGL_Zffh(quarks[UP]),2.0)
                +30559313.  * pow(deltaGR_Zffh(quarks[UP]),2.0)
                +19107837.  * pow(deltaGL_Zffh(quarks[DOWN]),2.0)
                +19109134.  * pow(deltaGR_Zffh(quarks[DOWN]),2.0);
        }
        
    } else if (sqrt_s == 14.0) {
        mu += +2.477 * deltaGL_f(quarks[UP])
                -1.103 * deltaGR_f(quarks[UP])
                -2.226 * deltaGL_f(quarks[DOWN])
                +0.405 * deltaGR_f(quarks[DOWN])
                +3321.75 * deltaG1_hZZ()
                -3494.38 * deltaG2_hZZ()
                +0.059 * deltaG3_hZZ()
                +481.727 * deltaG1_hZA()
                -967.239 * deltaG2_hZA()
                +5106.92 * deltaGL_Zffh(quarks[UP])
                -2270.81 * deltaGR_Zffh(quarks[UP])
                -4434.64 * deltaGL_Zffh(quarks[DOWN])
                +807.186 * deltaGR_Zffh(quarks[DOWN]);
        if (FlagQuadraticTerms) {
            //Add contributions that are quadratic in the effective coefficients
            //(Only valid under the assumptions of one dim 6 operator at a time)
            mu +=  +3.579 * pow(deltaGL_f(quarks[UP]),2.0)
                +3.58  * pow(deltaGR_f(quarks[UP]),2.0)
                +2.631  * pow(deltaGL_f(quarks[DOWN]),2.0)
                +2.629  * pow(deltaGR_f(quarks[DOWN]),2.0)
                +4609160.  * pow(deltaG1_hZZ(),2.0)
                +7946470.  * pow(deltaG2_hZZ(),2.0)
                +0.001  * pow(deltaG3_hZZ(),2.0)
                +683466.  * pow(deltaG1_hZA(),2.0)
                +5019397.  * pow(deltaG2_hZA(),2.0)
                +50036976.  * pow(deltaGL_Zffh(quarks[UP]),2.0)
                +50008570.  * pow(deltaGR_Zffh(quarks[UP]),2.0)
                +31660707.  * pow(deltaGL_Zffh(quarks[DOWN]),2.0)
                +31666009.  * pow(deltaGR_Zffh(quarks[DOWN]),2.0);
        }
        
    } else
        throw std::runtime_error("Bad argument in NPEffectiveGIMR::muZH()");

    if (mu < 0) return std::numeric_limits<double>::quiet_NaN();
    
    return mu;
}

double NPEffectiveGIMR::mueeZH(const double sqrt_s) const
{
    double mu = 1.0;

    if (sqrt_s == 0.24) {
        mu += 2690.84 * deltaG1_hZZ()
                -1951.93 * deltaG2_hZZ()
                +0.059 * deltaG3_hZZ()
                +147.761 * deltaG1_hZA()
                -185.735 * deltaG2_hZA()
                -4217.73 * deltaGL_Zffh(leptons[ELECTRON])
                +3619.82 * deltaGR_Zffh(leptons[ELECTRON])
                -4.282 * deltaGL_f(leptons[ELECTRON])
                +3.674 * deltaGR_f(leptons[ELECTRON]);
        
        if (FlagQuadraticTerms) {
        //Add contributions that are quadratic in the effective coefficients
        //(Only valid under the assumptions of one dim 6 operator at a time)
            mu += +7.966 * pow(deltaGL_f(leptons[DOWN]),2.0)
                +7.966 * pow(deltaGR_f(leptons[DOWN]),2.0)
                +1841343.  * pow(deltaG1_hZZ(),2.0)
                +952412.  * pow(deltaG2_hZZ(),2.0)
                +0.001  * pow(deltaG3_hZZ(),2.0)
                +961714.  * pow(deltaG1_hZA(),2.0)
                +1520521.  * pow(deltaG2_hZA(),2.0)
                +7731703.  * pow(deltaGL_Zffh(leptons[DOWN]),2.0)
                +7731703.  * pow(deltaGR_Zffh(leptons[DOWN]),2.0);
            }

    } else if (sqrt_s == 0.25) {
        mu +=  2829.45 * deltaG1_hZZ()
                -2097.01 * deltaG2_hZZ()
                +0.059 * deltaG3_hZZ()
                +156.787 * deltaG1_hZA()
                -204.357 * deltaG2_hZA()
                -4635.23 * deltaGL_Zffh(leptons[ELECTRON])
                +3979.84 * deltaGR_Zffh(leptons[ELECTRON])
                -4.282 * deltaGL_f(leptons[ELECTRON])
                +3.674 * deltaGR_f(leptons[ELECTRON]);
        
        if (FlagQuadraticTerms) {
            //Add contributions that are quadratic in the effective coefficients
            //(Only valid under the assumptions of one dim 6 operator at a time)
            mu +=  +0.0;
        }
        
    } else if (sqrt_s == 0.35) {
        mu +=  3893.41 * deltaG1_hZZ()
                -3873.83 * deltaG2_hZZ()
                +0.059 * deltaG3_hZZ()
                +231.963 * deltaG1_hZA()
                -424.266 * deltaG2_hZA()
                -9763.89 * deltaGL_Zffh(leptons[ELECTRON])
                +8387.72 * deltaGR_Zffh(leptons[ELECTRON])
                -4.282 * deltaGL_f(leptons[ELECTRON])
                +3.674 * deltaGR_f(leptons[ELECTRON]);
        
        if (FlagQuadraticTerms) {
            //Add contributions that are quadratic in the effective coefficients
            //(Only valid under the assumptions of one dim 6 operator at a time)
            mu +=  +0.0;
        }
          
    } else if (sqrt_s == 0.5) {
        mu +=  4747.11 * deltaG1_hZZ()
                -7649.28 * deltaG2_hZZ()
                +0.059 * deltaG3_hZZ()
                +291.854 * deltaG1_hZA()
                -902.663 * deltaG2_hZA()
                -20668.9 * deltaGL_Zffh(leptons[ELECTRON])
                +17754.2 * deltaGR_Zffh(leptons[ELECTRON])
                -4.282 * deltaGL_f(leptons[ELECTRON])
                +3.674 * deltaGR_f(leptons[ELECTRON]);
        
        if (FlagQuadraticTerms) {
            //Add contributions that are quadratic in the effective coefficients
            //(Only valid under the assumptions of one dim 6 operator at a time)
            mu +=  +0.0;
        }
        
    } else if (sqrt_s == 1.0) {
        mu +=  5576.18 * deltaG1_hZZ()
                -29856.9 * deltaG2_hZZ()
                +0.059 * deltaG3_hZZ()
                +351.186 * deltaG1_hZA()
                -3727.98 * deltaG2_hZA()
                -84814.3 * deltaGL_Zffh(leptons[ELECTRON])
                +72844.9 * deltaGR_Zffh(leptons[ELECTRON])
                -4.282 * deltaGL_f(leptons[ELECTRON])
                +3.673 * deltaGR_f(leptons[ELECTRON]);
        
        if (FlagQuadraticTerms) {
            //Add contributions that are quadratic in the effective coefficients
            //(Only valid under the assumptions of one dim 6 operator at a time)
            mu +=  +0.0;
        }
        
    } else
        throw std::runtime_error("Bad argument in NPEffectiveGIMR::mueeZH()");
    
    if (mu < 0) return std::numeric_limits<double>::quiet_NaN();
    
    return mu;
}

double NPEffectiveGIMR::muVH(const double sqrt_s) const
{
    double sigmaWH_SM = computeSigmaWH(sqrt_s);
    double sigmaZH_SM = computeSigmaZH(sqrt_s);
    double sigmaWH = muWH(sqrt_s) * sigmaWH_SM;
    double sigmaZH = muZH(sqrt_s) * sigmaZH_SM;
    double mu = ((sigmaWH + sigmaZH) / (sigmaWH_SM + sigmaZH_SM));
    
    if (mu < 0) return std::numeric_limits<double>::quiet_NaN();
    
    return mu;
}

double NPEffectiveGIMR::muVBFpVH(const double sqrt_s) const
{
    double sigmaWH_SM = computeSigmaWH(sqrt_s);
    double sigmaZH_SM = computeSigmaZH(sqrt_s);
    double sigmaVBF_SM = computeSigmaVBF(sqrt_s);
    double sigmaWH = muWH(sqrt_s) * sigmaWH_SM;
    double sigmaZH = muZH(sqrt_s) * sigmaZH_SM;
    double sigmaVBF = muVBF(sqrt_s) * sigmaVBF_SM;
    double mu = ((sigmaWH + sigmaZH + sigmaVBF) / (sigmaWH_SM + sigmaZH_SM + sigmaVBF_SM));
    
    if (mu < 0) return std::numeric_limits<double>::quiet_NaN();
    
    return mu;
}

double NPEffectiveGIMR::muttH(const double sqrt_s) const
{
    double mu = 1.0;
    if (sqrt_s == 1.96) {
        mu += -2.863 * (1. + ettH2_Htt ) * deltaG_hff(quarks[TOP]).real()
            +1737.35 * (1. + ettH2_Hgg ) * deltaG_hgg();
        
        if (FlagQuadraticTerms) {
            //Add contributions that are quadratic in the effective coefficients
            //(Only valid under the assumptions of one dim 6 operator at a time)
            mu += +2.036 * pow(deltaG_hff(quarks[TOP]).real(),2.0)
                +885586. * pow(deltaG_hgg(),2.0);            
        }
        
    } else if (sqrt_s == 7.0) {
        mu += -2.861 * (1. + ettH78_Htt ) * deltaG_hff(quarks[TOP]).real()
            +2583.3 * (1. + ettH78_Hgg ) * deltaG_hgg();
        
        if (FlagQuadraticTerms) {
            //Add contributions that are quadratic in the effective coefficients
            //(Only valid under the assumptions of one dim 6 operator at a time)
            mu += +2.073 * pow(deltaG_hff(quarks[TOP]).real(),2.0)
                +3909554. * pow(deltaG_hgg(),2.0);
        }
        
    } else if (sqrt_s == 8.0) {
        mu += -2.861 * (1. + ettH78_Htt ) * deltaG_hff(quarks[TOP]).real()
            +2636.88 * (1. + ettH78_Hgg ) * deltaG_hgg();
        
        if (FlagQuadraticTerms) {
            //Add contributions that are quadratic in the effective coefficients
            //(Only valid under the assumptions of one dim 6 operator at a time)
            mu += +1.963 * pow(deltaG_hff(quarks[TOP]).real(),2.0)
                +4367338. * pow(deltaG_hgg(),2.0);
        }
        
    } else if (sqrt_s == 14.0) {
        mu += -2.861 * deltaG_hff(quarks[TOP]).real()
            +2769.79 * deltaG_hgg();
        
        if (FlagQuadraticTerms) {
            //Add contributions that are quadratic in the effective coefficients
            //(Only valid under the assumptions of one dim 6 operator at a time)
            mu +=  +2.012 * pow(deltaG_hff(quarks[TOP]).real(),2.0)
            +5689423. * pow(deltaG_hgg(),2.0);
        }
        
    } else
        throw std::runtime_error("Bad argument in NPEffectiveGIMR::muttH()");

    if (mu < 0) return std::numeric_limits<double>::quiet_NaN();
    
    return mu;
}

double NPEffectiveGIMR::muggHpttH(const double sqrt_s) const
{
    double sigmaggH_SM = computeSigmaggH(sqrt_s);
    double sigmattH_SM = computeSigmattH(sqrt_s);
    double sigmaggH = muggH(sqrt_s) * sigmaggH_SM;
    double sigmattH = muttH(sqrt_s) * sigmattH_SM;

    double mu = ((sigmaggH + sigmattH) / (sigmaggH_SM + sigmattH_SM));
    
    if (mu < 0) return std::numeric_limits<double>::quiet_NaN();
    
    return mu;
} 

double NPEffectiveGIMR::mueettH(const double sqrt_s) const
{
    double mu = 1.0;
    if (sqrt_s == 0.5) {
        mu +=  85.139 * deltaG1_hZZ()
                -51.41 * deltaG2_hZZ()
                +0. * deltaG3_hZZ()
                +276.673 * deltaG1_hZA()
                -159.708 * deltaG2_hZA()
                +1017.44 * deltaG_hAA()
                -2.833 * deltaG_hff(quarks[TOP]).real()
                -287.92 * deltaGL_Zffh(leptons[ELECTRON])
                -132.849 * deltaGR_Zffh(leptons[ELECTRON])
                +84.883 * deltaGL_Zffh(quarks[TOP])
                +72.935 * deltaGR_Zffh(quarks[TOP])
                +863716. * deltaG_hAff(quarks[TOP]).real()
                +154393. * deltaG_hZff(quarks[TOP]).real()
                -1.046 * deltaGL_f(leptons[ELECTRON])
                -0.608 * deltaGR_f(leptons[ELECTRON])
                +0.663 * deltaGL_f(quarks[TOP])
                +0.585 * deltaGR_f(quarks[TOP])
                +8833.35 * deltaG_Aff(quarks[TOP]).real()
                +1650.94 * deltaG_Zff(quarks[TOP]).real();
        
        if (FlagQuadraticTerms) {
            //Add contributions that are quadratic in the effective coefficients
            //(Only valid under the assumptions of one dim 6 operator at a time)
            mu +=  +0.0;
        }
        
    } else if (sqrt_s == 1.0) {
        mu +=  446.758 * deltaG1_hZZ()
                -1500.92 * deltaG2_hZZ()
                +0.003 * deltaG3_hZZ()
                +657.283 * deltaG1_hZA()
                -572.102 * deltaG2_hZA()
                +2443.18 * deltaG_hAA()
                -2.701 * deltaG_hff(quarks[TOP]).real()
                -4591.53 * deltaGL_Zffh(leptons[ELECTRON])
                +2945.96 * deltaGR_Zffh(leptons[ELECTRON])
                +251.003 * deltaGL_Zffh(quarks[TOP])
                +49.581 * deltaGR_Zffh(quarks[TOP])
                +3025550. * deltaG_hAff(quarks[TOP]).real()
                +519896. * deltaG_hZff(quarks[TOP]).real()
                -1.426 * deltaGL_f(leptons[ELECTRON])
                -0.041 * deltaGR_f(leptons[ELECTRON])
                +1.066 * deltaGL_f(quarks[TOP])
                -0.038 * deltaGR_f(quarks[TOP])
                +12745.4 * deltaG_Aff(quarks[TOP]).real()
                +2238.61 * deltaG_Zff(quarks[TOP]).real();
        
        if (FlagQuadraticTerms) {
            //Add contributions that are quadratic in the effective coefficients
            //(Only valid under the assumptions of one dim 6 operator at a time)
            mu +=  +0.0;
        }
        
    } else
        throw std::runtime_error("Bad argument in NPEffectiveGIMR::mueettH()");

    if (mu < 0) return std::numeric_limits<double>::quiet_NaN();
    
    return mu;
}

double NPEffectiveGIMR::BrHggRatio() const
{
    double Br = 1.0;
    
    Br += deltaGammaHggRatio1() - deltaGammaTotalRatio1();
    
    if (FlagQuadraticTerms) {
        //Add contributions that are quadratic in the effective coefficients
        //(Only valid under the assumptions of one dim 6 operator at a time)
        Br += - deltaGammaHggRatio1() * deltaGammaTotalRatio1()
                + deltaGammaHggRatio2() - deltaGammaTotalRatio2()
                + pow(deltaGammaTotalRatio1(),2.0);            
        }
    
    if (Br < 0) return std::numeric_limits<double>::quiet_NaN();
    
    return Br;

}

double NPEffectiveGIMR::BrHWWRatio() const
{
    double Br = 1.0;
    
    Br += deltaGammaHWWRatio1() - deltaGammaTotalRatio1();
    
    if (FlagQuadraticTerms) {
        //Add contributions that are quadratic in the effective coefficients
        //(Only valid under the assumptions of one dim 6 operator at a time)
        Br += - deltaGammaHWWRatio1() * deltaGammaTotalRatio1()
                + deltaGammaHWWRatio2() - deltaGammaTotalRatio2()
                + pow(deltaGammaTotalRatio1(),2.0);            
        }
    
    if (Br < 0) return std::numeric_limits<double>::quiet_NaN();
    
    return Br;

}

double NPEffectiveGIMR::BrHZZRatio() const
{
    double Br = 1.0;
    
    Br += deltaGammaHZZRatio1() - deltaGammaTotalRatio1();
    
    if (FlagQuadraticTerms) {
        //Add contributions that are quadratic in the effective coefficients
        //(Only valid under the assumptions of one dim 6 operator at a time)
        Br += - deltaGammaHZZRatio1() * deltaGammaTotalRatio1()
                + deltaGammaHZZRatio2() - deltaGammaTotalRatio2()
                + pow(deltaGammaTotalRatio1(),2.0);            
        }
    
    if (Br < 0) return std::numeric_limits<double>::quiet_NaN();
    
    return Br;

}

double NPEffectiveGIMR::BrHZgaRatio() const
{
    double Br = 1.0;
    
    Br += deltaGammaHZgaRatio1() - deltaGammaTotalRatio1();
    
    if (FlagQuadraticTerms) {
        //Add contributions that are quadratic in the effective coefficients
        //(Only valid under the assumptions of one dim 6 operator at a time)
        Br += - deltaGammaHZgaRatio1() * deltaGammaTotalRatio1()
                + deltaGammaHZgaRatio2() - deltaGammaTotalRatio2()
                + pow(deltaGammaTotalRatio1(),2.0);            
        }
    
    if (Br < 0) return std::numeric_limits<double>::quiet_NaN();
    
    return Br;

}

double NPEffectiveGIMR::BrHgagaRatio() const
{
    double Br = 1.0;
    
    Br += deltaGammaHgagaRatio1() - deltaGammaTotalRatio1();
    
    if (FlagQuadraticTerms) {
        //Add contributions that are quadratic in the effective coefficients
        //(Only valid under the assumptions of one dim 6 operator at a time)
        Br += - deltaGammaHgagaRatio1() * deltaGammaTotalRatio1()
                + deltaGammaHgagaRatio2() - deltaGammaTotalRatio2()
                + pow(deltaGammaTotalRatio1(),2.0);            
        }
    
    if (Br < 0) return std::numeric_limits<double>::quiet_NaN();
    
    return Br;

}

double NPEffectiveGIMR::BrHmumuRatio() const
{
    double Br = 1.0;
    
    Br += deltaGammaHmumuRatio1() - deltaGammaTotalRatio1();
    
    if (FlagQuadraticTerms) {
        //Add contributions that are quadratic in the effective coefficients
        //(Only valid under the assumptions of one dim 6 operator at a time)
        Br += - deltaGammaHmumuRatio1() * deltaGammaTotalRatio1()
                + deltaGammaHmumuRatio2() - deltaGammaTotalRatio2()
                + pow(deltaGammaTotalRatio1(),2.0);            
        }
    
    if (Br < 0) return std::numeric_limits<double>::quiet_NaN();
    
    return Br;

}

double NPEffectiveGIMR::BrHtautauRatio() const
{
    double Br = 1.0;
    
    Br += deltaGammaHtautauRatio1() - deltaGammaTotalRatio1();
    
    if (FlagQuadraticTerms) {
        //Add contributions that are quadratic in the effective coefficients
        //(Only valid under the assumptions of one dim 6 operator at a time)
        Br += - deltaGammaHtautauRatio1() * deltaGammaTotalRatio1()
                + deltaGammaHtautauRatio2() - deltaGammaTotalRatio2()
                + pow(deltaGammaTotalRatio1(),2.0);            
        }
    
    if (Br < 0) return std::numeric_limits<double>::quiet_NaN();
    
    return Br;

}

double NPEffectiveGIMR::BrHccRatio() const
{
    double Br = 1.0;
    
    Br += deltaGammaHccRatio1() - deltaGammaTotalRatio1();
    
    if (FlagQuadraticTerms) {
        //Add contributions that are quadratic in the effective coefficients
        //(Only valid under the assumptions of one dim 6 operator at a time)
        Br += - deltaGammaHccRatio1() * deltaGammaTotalRatio1()
                + deltaGammaHccRatio2() - deltaGammaTotalRatio2()
                + pow(deltaGammaTotalRatio1(),2.0);            
        }
    
    if (Br < 0) return std::numeric_limits<double>::quiet_NaN();
    
    return Br;

}

double NPEffectiveGIMR::BrHbbRatio() const
{
    double Br = 1.0;
    
    Br += deltaGammaHbbRatio1() - deltaGammaTotalRatio1();
    
    if (FlagQuadraticTerms) {
        //Add contributions that are quadratic in the effective coefficients
        //(Only valid under the assumptions of one dim 6 operator at a time)
        Br += - deltaGammaHbbRatio1() * deltaGammaTotalRatio1()
                + deltaGammaHbbRatio2() - deltaGammaTotalRatio2()
                + pow(deltaGammaTotalRatio1(),2.0);            
        }
    
    if (Br < 0) return std::numeric_limits<double>::quiet_NaN();
    
    return Br;

}

double NPEffectiveGIMR::computeGammaTotalRatio() const
{
    return (trueSM.computeBrHtogg() * GammaHggRatio()
            + trueSM.computeBrHtoWW() * GammaHWWRatio()
            + trueSM.computeBrHtoZZ() * GammaHZZRatio()
            + trueSM.computeBrHtoZga() * GammaHZgaRatio()
            + trueSM.computeBrHtogaga() * GammaHgagaRatio()
            + trueSM.computeBrHtomumu() * GammaHmumuRatio()
            + trueSM.computeBrHtotautau() * GammaHtautauRatio()
            + trueSM.computeBrHtocc() * GammaHccRatio()
            + trueSM.computeBrHtobb() * GammaHbbRatio());
}

double NPEffectiveGIMR::deltaGammaTotalRatio1() const
{
    return (trueSM.computeBrHtogg() * deltaGammaHggRatio1()
            + trueSM.computeBrHtoWW() * deltaGammaHWWRatio1()
            + trueSM.computeBrHtoZZ() * deltaGammaHZZRatio1()
            + trueSM.computeBrHtoZga() * deltaGammaHZgaRatio1()
            + trueSM.computeBrHtogaga() * deltaGammaHgagaRatio1()
            + trueSM.computeBrHtomumu() * deltaGammaHmumuRatio1()
            + trueSM.computeBrHtotautau() * deltaGammaHtautauRatio1()
            + trueSM.computeBrHtocc() * deltaGammaHccRatio1()
            + trueSM.computeBrHtobb() * deltaGammaHbbRatio1());
}

double NPEffectiveGIMR::deltaGammaTotalRatio2() const
{
    return (trueSM.computeBrHtogg() * deltaGammaHggRatio2()
            + trueSM.computeBrHtoWW() * deltaGammaHWWRatio2()
            + trueSM.computeBrHtoZZ() * deltaGammaHZZRatio2()
            + trueSM.computeBrHtoZga() * deltaGammaHZgaRatio2()
            + trueSM.computeBrHtogaga() * deltaGammaHgagaRatio2()
            + trueSM.computeBrHtomumu() * deltaGammaHmumuRatio2()
            + trueSM.computeBrHtotautau() * deltaGammaHtautauRatio2()
            + trueSM.computeBrHtocc() * deltaGammaHccRatio2()
            + trueSM.computeBrHtobb() * deltaGammaHbbRatio2());
}

double NPEffectiveGIMR::GammaHggRatio() const
{
    double width = 1.0;

    width += deltaGammaHggRatio1();
    
    if (FlagQuadraticTerms) {
            //Add contributions that are quadratic in the effective coefficients
            //(Only valid under the assumptions of one dim 6 operator at a time)
        width += deltaGammaHggRatio2();            
        }
    
    return width;
    
}

double NPEffectiveGIMR::deltaGammaHggRatio1() const
{
    return ( +151669. * deltaG_hgg()
            -3.006 * deltaG_hff(quarks[TOP]).real()
            +5.853 * deltaG_hff(quarks[BOTTOM]).real()
            +4.71 * deltaG_hff(quarks[CHARM]).real() );
}

double NPEffectiveGIMR::deltaGammaHggRatio2() const
{
    //Contributions that are quadratic in the effective coefficients
    //(Only valid under the assumptions of one dim 6 operator at a time)
    return ( +5879800851. * pow(deltaG_hgg(),2.0)
            +2.284 * pow(deltaG_hff(quarks[TOP]).real(),2.0)
            +40.881 * pow(deltaG_hff(quarks[BOTTOM]).real(),2.0)
            +2.17 * pow(deltaG_hff(quarks[CHARM]).real(),2.0) );            
    
}

double NPEffectiveGIMR::GammaHWWRatio() const
{
    double width = 1.0;

    width += deltaGammaHWWRatio1();
    
    if (FlagQuadraticTerms) {
            //Add contributions that are quadratic in the effective coefficients
            //(Only valid under the assumptions of one dim 6 operator at a time)
        width += deltaGammaHWWRatio2();
        }
    
    return width;
    
}

double NPEffectiveGIMR::deltaGammaHWWRatio1() const
{

    return ( -183.404 * deltaG1_hWW()
            -274.568 * deltaG2_hWW()
            +0.039 * deltaG3_hWW() );
    
}

double NPEffectiveGIMR::deltaGammaHWWRatio2() const
{
    //Contributions that are quadratic in the effective coefficients
    //(Only valid under the assumptions of one dim 6 operator at a time)
    return ( +1267. * pow(deltaG1_hWW(),2.0)
            +868.393 * pow(deltaG2_hWW(),2.0) );
    
}

double NPEffectiveGIMR::GammaHZZRatio() const
{
    double width = 1.0;

    width += deltaGammaHZZRatio1();
    
    if (FlagQuadraticTerms) {
            //Add contributions that are quadratic in the effective coefficients
            //(Only valid under the assumptions of one dim 6 operator at a time)
        width += deltaGammaHZZRatio2();            
        }
    
    return width;
    
}

double NPEffectiveGIMR::deltaGammaHZZRatio1() const
{

    return ( -246.654 * deltaG1_hZZ()
            -240.846 * deltaG2_hZZ()
            +0.059 * deltaG3_hZZ() );
    
}

double NPEffectiveGIMR::deltaGammaHZZRatio2() const
{
    //Contributions that are quadratic in the effective coefficients
    //(Only valid under the assumptions of one dim 6 operator at a time)
    return ( +6391.57 * pow(deltaG1_hZZ(),2.0)
            +2088.67 * pow(deltaG2_hZZ(),2.0)
            +0.001 * pow(deltaG3_hZZ(),2.0) );            
    
}

double NPEffectiveGIMR::GammaHZgaRatio() const
{
    double width = 1.0;

    width += deltaGammaHZgaRatio1();
    
    if (FlagQuadraticTerms) {
            //Add contributions that are quadratic in the effective coefficients
            //(Only valid under the assumptions of one dim 6 operator at a time)
        width += deltaGammaHZgaRatio2();            
        }
    
    return width;
    
}

double NPEffectiveGIMR::deltaGammaHZgaRatio1() const
{

    return ( -71321.5 * deltaG1_hZA()
            +0.041 * deltaG3_hWW()
            +0.172 * deltaG_hff(quarks[TOP]).real()
            -0.301 * deltaG_hff(quarks[BOTTOM]).real()
            +0.196 * deltaG_hff(leptons[TAU]).real()
            +0.232 * deltaG_hff(quarks[CHARM]).real() );
    
}

double NPEffectiveGIMR::deltaGammaHZgaRatio2() const
{
    //Contributions that are quadratic in the effective coefficients
    //(Only valid under the assumptions of one dim 6 operator at a time)
    return ( +1271853409. * pow(deltaG1_hZA(),2.0)
            +0.003 * pow(deltaG_hff(quarks[TOP]).real(),2.0)
            +3.539 * pow(deltaG_hff(quarks[BOTTOM]).real(),2.0)
            -14.568 * pow(deltaG_hff(leptons[TAU]).real(),2.0)
            -31.197 * pow(deltaG_hff(quarks[CHARM]).real(),2.0) );            
    
}

double NPEffectiveGIMR::GammaHgagaRatio() const
{
    double width = 1.0;

    width += deltaGammaHgagaRatio1();
    
    if (FlagQuadraticTerms) {
            //Add contributions that are quadratic in the effective coefficients
            //(Only valid under the assumptions of one dim 6 operator at a time)
        width += deltaGammaHgagaRatio2();            
        }
    
    return width;
    
}

double NPEffectiveGIMR::deltaGammaHgagaRatio1() const
{
    return ( -257366. * deltaG_hAA()
            +0.049 * deltaG3_hWW()
            +0.761 * deltaG_hff(quarks[TOP]).real()
            -0.441 * deltaG_hff(quarks[BOTTOM]).real()
            -1.087 * deltaG_hff(leptons[TAU]).real()
            -0.646 * deltaG_hff(quarks[CHARM]).real() );
    
}

double NPEffectiveGIMR::deltaGammaHgagaRatio2() const
{
    //Contributions that are quadratic in the effective coefficients
    //(Only valid under the assumptions of one dim 6 operator at a time)
    return ( +16479108529. * pow(deltaG_hAA(),2.0)
            +0.001 * pow(deltaG3_hWW(),2.0)
            +0.146 * pow(deltaG_hff(quarks[TOP]).real(),2.0)
            +1.828 * pow(deltaG_hff(quarks[BOTTOM]).real(),2.0)
            +6.672 * pow(deltaG_hff(leptons[TAU]).real(),2.0)
            +9.962 * pow(deltaG_hff(quarks[CHARM]).real(),2.0) );            
    
}

double NPEffectiveGIMR::GammaHmumuRatio() const
{
    double width = 1.0;

    width += deltaGammaHmumuRatio1();
    
    if (FlagQuadraticTerms) {
            //Add contributions that are quadratic in the effective coefficients
            //(Only valid under the assumptions of one dim 6 operator at a time)
        width += deltaGammaHmumuRatio2();
        }
    
    return width;
    
}

double NPEffectiveGIMR::deltaGammaHmumuRatio1() const
{
    return ( -4653.43 * deltaG_hff(leptons[MU]).real() );
        
}

double NPEffectiveGIMR::deltaGammaHmumuRatio2() const
{
    //Contributions that are quadratic in the effective coefficients
    //(Only valid under the assumptions of one dim 6 operator at a time)
    return 0.0;
    
}

double NPEffectiveGIMR::GammaHtautauRatio() const
{
    double width = 1.0;

    width += deltaGammaHtautauRatio1();
    
    if (FlagQuadraticTerms) {
            //Add contributions that are quadratic in the effective coefficients
            //(Only valid under the assumptions of one dim 6 operator at a time)
        width += deltaGammaHtautauRatio2();            
        }
    
    return width;
    
}

double NPEffectiveGIMR::deltaGammaHtautauRatio1() const
{
    return ( -277.458 * deltaG_hff(leptons[TAU]).real() );
        
}

double NPEffectiveGIMR::deltaGammaHtautauRatio2() const
{
    //Contributions that are quadratic in the effective coefficients
    //(Only valid under the assumptions of one dim 6 operator at a time)
    return ( +19223. * pow(deltaG_hff(leptons[TAU]).real(),2.0) );            
    
}

double NPEffectiveGIMR::GammaHccRatio() const
{
    double width = 1.0;

    width += deltaGammaHccRatio1();
    
    if (FlagQuadraticTerms) {
            //Add contributions that are quadratic in the effective coefficients
            //(Only valid under the assumptions of one dim 6 operator at a time)
        width += deltaGammaHccRatio2();            
        }
    
    return width;
    
}

double NPEffectiveGIMR::deltaGammaHccRatio1() const
{
    return ( -383.036 * deltaG_hff(quarks[CHARM]).real() );   
}

double NPEffectiveGIMR::deltaGammaHccRatio2() const
{
    //Contributions that are quadratic in the effective coefficients
    //(Only valid under the assumptions of one dim 6 operator at a time)
    return ( +36709.1 * pow(deltaG_hff(quarks[CHARM]).real(),2.0) );            
    
}

double NPEffectiveGIMR::GammaHbbRatio() const
{
    double width = 1.0;
    
    width += deltaGammaHbbRatio1();
    
    if (FlagQuadraticTerms) {
            //Add contributions that are quadratic in the effective coefficients
            //(Only valid under the assumptions of one dim 6 operator at a time)
        width += deltaGammaHbbRatio2();        
        }
    
    return width;
}

double NPEffectiveGIMR::deltaGammaHbbRatio1() const
{    
    return ( -0.013 * deltaG_hff(quarks[TOP]).real()
            -117.431 * deltaG_hff(quarks[BOTTOM]).real() );
}

double NPEffectiveGIMR::deltaGammaHbbRatio2() const
{
    //Contributions that are quadratic in the effective coefficients
    //(Only valid under the assumptions of one dim 6 operator at a time)
    return ( +3443.96 * pow(deltaG_hff(quarks[BOTTOM]).real(),2.0) );        
    
}

///////////////////////////////////////////////////////////////////////////////

double NPEffectiveGIMR::CLL_mu() const
{
    return ((CLL_1122 + CLL_1221) / LambdaNP2);
}

double NPEffectiveGIMR::CLL_tau() const
{
    return ((CLL_1133 + CLL_1331) / LambdaNP2);
}

double NPEffectiveGIMR::CLL_up() const
{
    return ((CLQ1_1111 - CLQ3_1111) / LambdaNP2);
}

double NPEffectiveGIMR::CLL_down() const
{
    return ((CLQ1_1111 + CLQ3_1111) / LambdaNP2);
}

double NPEffectiveGIMR::CLL_charm() const
{
    return ((CLQ1_1122 - CLQ3_1122) / LambdaNP2);
}

double NPEffectiveGIMR::CLL_strange() const
{
    return ((CLQ1_1122 + CLQ3_1122) / LambdaNP2);
}

double NPEffectiveGIMR::CLL_bottom() const
{
    return ((CLQ1_1133 + CLQ3_1133) / LambdaNP2);
}

double NPEffectiveGIMR::CLR_mu() const
{
    return ((CLe_1122 + CLe_2211) / LambdaNP2);
}

double NPEffectiveGIMR::CLR_tau() const
{
    return ((CLe_1133 + CLe_3311) / LambdaNP2);
}

double NPEffectiveGIMR::CLR_up() const
{
    return (CLu_1111 / LambdaNP2);
}

double NPEffectiveGIMR::CLR_down() const
{
    return (CLd_1111 / LambdaNP2);
}

double NPEffectiveGIMR::CLR_charm() const
{
    return (CLu_1122 / LambdaNP2);
}

double NPEffectiveGIMR::CLR_strange() const
{
    return (CLd_1122 / LambdaNP2);
}

double NPEffectiveGIMR::CLR_bottom() const
{
    return (CLd_1133 / LambdaNP2);
}

double NPEffectiveGIMR::CRL_mu() const
{
    return (CLe_1122 / LambdaNP2);
}

double NPEffectiveGIMR::CRL_tau() const
{
    return (CLe_1133 / LambdaNP2);
}

double NPEffectiveGIMR::CRL_up() const
{
    return (CQe_1111 / LambdaNP2);
}

double NPEffectiveGIMR::CRL_down() const
{
    return (CQe_1111 / LambdaNP2);
}

double NPEffectiveGIMR::CRL_charm() const
{
    return (CQe_1122 / LambdaNP2);
}

double NPEffectiveGIMR::CRL_strange() const
{
    return (CQe_1122 / LambdaNP2);
}

double NPEffectiveGIMR::CRL_bottom() const
{
    return (CQe_1133 / LambdaNP2);
}

double NPEffectiveGIMR::CRR_mu() const
{
    return (Cee_1122 / LambdaNP2);
}

double NPEffectiveGIMR::CRR_tau() const
{
    return (Cee_1133 / LambdaNP2);
}


double NPEffectiveGIMR::CRR_up() const
{
    return (Ceu_1111 / LambdaNP2);
}

double NPEffectiveGIMR::CRR_down() const
{
    return (Ced_1111 / LambdaNP2);
}

double NPEffectiveGIMR::CRR_charm() const
{
    return (Ceu_1122 / LambdaNP2);
}

double NPEffectiveGIMR::CRR_strange() const
{
    return (Ced_1122 / LambdaNP2);
}

double NPEffectiveGIMR::CRR_bottom() const
{
    return (Ced_1133 / LambdaNP2);
}



double NPEffectiveGIMR::LEP2sigmaMu(const double s) const
{
   double sigma_mu;
   double GIMRParam[10] = {CLL_mu(), CLR_mu(), CRL_mu(), CRR_mu(), deltaGamma_Z(), 
                          deltaGR_f(leptons[MU]), 
                          deltaGL_f(leptons[MU]),
                          deltaGR_f(leptons[ELECTRON]),
                          deltaGL_f(leptons[ELECTRON]),
                          deltaMz2()};
   
   sigma_mu = trueSM.LEP2sigmaMu(s) + 
              myLEP2GIMR.sigma_l_LEP2_GIMR(QCD::lepton(MU), s, GIMRParam);
   
   return sigma_mu;
}

double NPEffectiveGIMR::LEP2sigmaTau(const double s) const
{
   double sigma_tau;
   double GIMRParam[10] = {CLL_tau(), CLR_tau(), CRL_tau(), CRR_tau(), deltaGamma_Z(), 
                          deltaGR_f(leptons[TAU]), 
                          deltaGL_f(leptons[TAU]),
                          deltaGR_f(leptons[ELECTRON]),
                          deltaGL_f(leptons[ELECTRON]),
                          deltaMz2()};
   
   sigma_tau = trueSM.LEP2sigmaTau(s) + 
               myLEP2GIMR.sigma_l_LEP2_GIMR(QCD::lepton(TAU), s, GIMRParam);   
   
   return sigma_tau;
}

double NPEffectiveGIMR::LEP2sigmaHadron(const double s) const
{
   double sigma_had;
   double GIMRParamUP[10] = {CLL_up(), CLR_up(), CRL_up(), CRR_up(), deltaGamma_Z(), 
                          deltaGR_f(quarks[UP]), 
                          deltaGL_f(quarks[UP]),
                          deltaGR_f(leptons[ELECTRON]),
                          deltaGL_f(leptons[ELECTRON]),
                          deltaMz2()};
   
   double GIMRParamDOWN[10] = {CLL_down(), CLR_down(), CRL_down(), CRR_down(), deltaGamma_Z(), 
                          deltaGR_f(quarks[DOWN]), 
                          deltaGL_f(quarks[DOWN]),
                          deltaGR_f(leptons[ELECTRON]),
                          deltaGL_f(leptons[ELECTRON]),
                          deltaMz2()};
      
   double GIMRParamSTRANGE[10] = {CLL_strange(), CLR_strange(), CRL_strange(), CRR_strange(), deltaGamma_Z(), 
                          deltaGR_f(quarks[STRANGE]), 
                          deltaGL_f(quarks[STRANGE]),
                          deltaGR_f(leptons[ELECTRON]),
                          deltaGL_f(leptons[ELECTRON]),
                          deltaMz2()};
         
   double GIMRParamCHARM[10] = {CLL_charm(), CLR_charm(), CRL_charm(), CRR_charm(), deltaGamma_Z(), 
                          deltaGR_f(quarks[CHARM]), 
                          deltaGL_f(quarks[CHARM]),
                          deltaGR_f(leptons[ELECTRON]),
                          deltaGL_f(leptons[ELECTRON]),
                          deltaMz2()};
            
   double GIMRParamBOTTOM[10] = {CLL_bottom(), CLR_bottom(), CRL_bottom(), CRR_bottom(), deltaGamma_Z(), 
                          deltaGR_f(quarks[BOTTOM]), 
                          deltaGL_f(quarks[BOTTOM]),
                          deltaGR_f(leptons[ELECTRON]),
                          deltaGL_f(leptons[ELECTRON]),
                          deltaMz2()};
           
           
   sigma_had = trueSM.LEP2sigmaHadron(s) + 
              myLEP2GIMR.sigma_q_LEP2_GIMR(QCD::UP, s, GIMRParamUP) +
              myLEP2GIMR.sigma_q_LEP2_GIMR(QCD::DOWN, s, GIMRParamDOWN) +
              myLEP2GIMR.sigma_q_LEP2_GIMR(QCD::CHARM, s, GIMRParamCHARM) +
              myLEP2GIMR.sigma_q_LEP2_GIMR(QCD::STRANGE, s, GIMRParamSTRANGE) +
              myLEP2GIMR.sigma_q_LEP2_GIMR(QCD::BOTTOM, s, GIMRParamBOTTOM);
   
   return sigma_had;
}


double NPEffectiveGIMR::LEP2AFBmu(const double s) const
{
   double AFB_mu = 0., sigma = 0.0,dFminusB =0.0;
   double GIMRParam[10] = {CLL_mu(), CLR_mu(), CRL_mu(), CRR_mu(), deltaGamma_Z(), 
                          deltaGR_f(leptons[MU]), 
                          deltaGL_f(leptons[MU]),
                          deltaGR_f(leptons[ELECTRON]),
                          deltaGL_f(leptons[ELECTRON]),
                          deltaMz2()};
   
   sigma = trueSM.LEP2sigmaMu(s);
   dFminusB = myLEP2GIMR.sigmaFminusB_l_LEP2_GIMR(QCD::lepton(MU), s, GIMRParam);
   
   AFB_mu = trueSM.LEP2AFBmu(s) * (1. - myLEP2GIMR.sigma_l_LEP2_GIMR(QCD::lepton(MU), s, GIMRParam) / sigma)
                   + dFminusB / sigma;
   
   return AFB_mu;
}

double NPEffectiveGIMR::LEP2AFBtau(const double s) const
{
   double AFB_tau, sigma = 0.0,dFminusB =0.0;
   double GIMRParam[10] = {CLL_tau(), CLR_tau(), CRL_tau(), CRR_tau(), deltaGamma_Z(), 
                          deltaGR_f(leptons[TAU]), 
                          deltaGL_f(leptons[TAU]),
                          deltaGR_f(leptons[ELECTRON]),
                          deltaGL_f(leptons[ELECTRON]),
                          deltaMz2()};
   sigma = trueSM.LEP2sigmaTau(s);
   dFminusB = myLEP2GIMR.sigmaFminusB_l_LEP2_GIMR(QCD::lepton(TAU), s, GIMRParam);

   AFB_tau = trueSM.LEP2AFBtau(s) * (1. - myLEP2GIMR.sigma_l_LEP2_GIMR(QCD::lepton(TAU), s, GIMRParam) / sigma)
                   + dFminusB / sigma;
   
   return AFB_tau;
}

double NPEffectiveGIMR::LEP2AFBcharm(const double s) const
{
   double AFB_charm, sigma = 0.0,dFminusB =0.0;
   double GIMRParam[10] = {CLL_charm(), CLR_charm(), CRL_charm(), CRR_charm(), deltaGamma_Z(), 
                          deltaGR_f(quarks[CHARM]), 
                          deltaGL_f(quarks[CHARM]),
                          deltaGR_f(leptons[ELECTRON]),
                          deltaGL_f(leptons[ELECTRON]),
                          deltaMz2()};
   sigma = trueSM.LEP2sigmaCharm(s);
   dFminusB = myLEP2GIMR.sigmaFminusB_q_LEP2_GIMR(QCD::quark(CHARM), s, GIMRParam);
   
   AFB_charm = trueSM.LEP2AFBcharm(s) * (1. - myLEP2GIMR.sigma_q_LEP2_GIMR(QCD::quark(CHARM), s, GIMRParam) / sigma)
                   + dFminusB / sigma;
   
   return AFB_charm;
}

double NPEffectiveGIMR::LEP2AFBbottom(const double s) const
{
   double AFB_bottom, sigma = 0.0,dFminusB =0.0;;
   double GIMRParam[10] = {CLL_bottom(), CLR_bottom(), CRL_bottom(), CRR_bottom(), deltaGamma_Z(), 
                          deltaGR_f(quarks[BOTTOM]), 
                          deltaGL_f(quarks[BOTTOM]),
                          deltaGR_f(leptons[ELECTRON]),
                          deltaGL_f(leptons[ELECTRON]),
                          deltaMz2()};
   sigma = trueSM.LEP2sigmaBottom(s);
   dFminusB = myLEP2GIMR.sigmaFminusB_q_LEP2_GIMR(QCD::quark(BOTTOM), s, GIMRParam);
   
   AFB_bottom = trueSM.LEP2AFBbottom(s) * (1. - myLEP2GIMR.sigma_q_LEP2_GIMR(QCD::quark(BOTTOM), s, GIMRParam) / sigma)
                   + dFminusB / sigma;
   
   return AFB_bottom;
}


double NPEffectiveGIMR::LEP2Rcharm(const double s) const
{
   double R_charm = 0., dsigma_c = 0.0, dsigma_had =0.0, sigma_hSM = 0.0;
   double GIMRParamUP[10] = {CLL_up(), CLR_up(), CRL_up(), CRR_up(), deltaGamma_Z(), 
                          deltaGR_f(quarks[UP]), 
                          deltaGL_f(quarks[UP]),
                          deltaGR_f(leptons[ELECTRON]),
                          deltaGL_f(leptons[ELECTRON]),
                          deltaMz2()};
   
   double GIMRParamDOWN[10] = {CLL_down(), CLR_down(), CRL_down(), CRR_down(), deltaGamma_Z(), 
                          deltaGR_f(quarks[DOWN]), 
                          deltaGL_f(quarks[DOWN]),
                          deltaGR_f(leptons[ELECTRON]),
                          deltaGL_f(leptons[ELECTRON]),
                          deltaMz2()};
      
   double GIMRParamSTRANGE[10] = {CLL_strange(), CLR_strange(), CRL_strange(), CRR_strange(), deltaGamma_Z(), 
                          deltaGR_f(quarks[STRANGE]), 
                          deltaGL_f(quarks[STRANGE]),
                          deltaGR_f(leptons[ELECTRON]),
                          deltaGL_f(leptons[ELECTRON]),
                          deltaMz2()};
         
   double GIMRParamCHARM[10] = {CLL_charm(), CLR_charm(), CRL_charm(), CRR_charm(), deltaGamma_Z(), 
                          deltaGR_f(quarks[CHARM]), 
                          deltaGL_f(quarks[CHARM]),
                          deltaGR_f(leptons[ELECTRON]),
                          deltaGL_f(leptons[ELECTRON]),
                          deltaMz2()};
            
   double GIMRParamBOTTOM[10] = {CLL_bottom(), CLR_bottom(), CRL_bottom(), CRR_bottom(), deltaGamma_Z(), 
                          deltaGR_f(quarks[BOTTOM]), 
                          deltaGL_f(quarks[BOTTOM]),
                          deltaGR_f(leptons[ELECTRON]),
                          deltaGL_f(leptons[ELECTRON]),
                          deltaMz2()};
   
   sigma_hSM = trueSM.LEP2sigmaHadron(s);
   
   dsigma_had = myLEP2GIMR.sigma_q_LEP2_GIMR(QCD::UP, s, GIMRParamUP) +
                myLEP2GIMR.sigma_q_LEP2_GIMR(QCD::DOWN, s, GIMRParamDOWN) +
                myLEP2GIMR.sigma_q_LEP2_GIMR(QCD::CHARM, s, GIMRParamCHARM) +
                myLEP2GIMR.sigma_q_LEP2_GIMR(QCD::STRANGE, s, GIMRParamSTRANGE) +
                myLEP2GIMR.sigma_q_LEP2_GIMR(QCD::BOTTOM, s, GIMRParamBOTTOM);
   
   dsigma_c = myLEP2GIMR.sigma_q_LEP2_GIMR(QCD::CHARM, s, GIMRParamCHARM);
           
   R_charm = trueSM.LEP2Rcharm(s) * (1 - dsigma_had / sigma_hSM) + dsigma_c / sigma_hSM;
   
   return R_charm;
}

double NPEffectiveGIMR::LEP2Rbottom(const double s) const
{
   double R_bottom = 0., dsigma_b = 0.0, dsigma_had =0.0, sigma_hSM = 0.0;
   double GIMRParamUP[10] = {CLL_up(), CLR_up(), CRL_up(), CRR_up(), deltaGamma_Z(), 
                          deltaGR_f(quarks[UP]), 
                          deltaGL_f(quarks[UP]),
                          deltaGR_f(leptons[ELECTRON]),
                          deltaGL_f(leptons[ELECTRON]),
                          deltaMz2()};
   
   double GIMRParamDOWN[10] = {CLL_down(), CLR_down(), CRL_down(), CRR_down(), deltaGamma_Z(), 
                          deltaGR_f(quarks[DOWN]), 
                          deltaGL_f(quarks[DOWN]),
                          deltaGR_f(leptons[ELECTRON]),
                          deltaGL_f(leptons[ELECTRON]),
                          deltaMz2()};
      
   double GIMRParamSTRANGE[10] = {CLL_strange(), CLR_strange(), CRL_strange(), CRR_strange(), deltaGamma_Z(), 
                          deltaGR_f(quarks[STRANGE]), 
                          deltaGL_f(quarks[STRANGE]),
                          deltaGR_f(leptons[ELECTRON]),
                          deltaGL_f(leptons[ELECTRON]),
                          deltaMz2()};
         
   double GIMRParamCHARM[10] = {CLL_charm(), CLR_charm(), CRL_charm(), CRR_charm(), deltaGamma_Z(), 
                          deltaGR_f(quarks[CHARM]), 
                          deltaGL_f(quarks[CHARM]),
                          deltaGR_f(leptons[ELECTRON]),
                          deltaGL_f(leptons[ELECTRON]),
                          deltaMz2()};
            
   double GIMRParamBOTTOM[10] = {CLL_bottom(), CLR_bottom(), CRL_bottom(), CRR_bottom(), deltaGamma_Z(), 
                          deltaGR_f(quarks[BOTTOM]), 
                          deltaGL_f(quarks[BOTTOM]),
                          deltaGR_f(leptons[ELECTRON]),
                          deltaGL_f(leptons[ELECTRON]),
                          deltaMz2()};
   
   sigma_hSM = trueSM.LEP2sigmaHadron(s);
   
   dsigma_had = myLEP2GIMR.sigma_q_LEP2_GIMR(QCD::UP, s, GIMRParamUP) +
                myLEP2GIMR.sigma_q_LEP2_GIMR(QCD::DOWN, s, GIMRParamDOWN) +
                myLEP2GIMR.sigma_q_LEP2_GIMR(QCD::CHARM, s, GIMRParamCHARM) +
                myLEP2GIMR.sigma_q_LEP2_GIMR(QCD::STRANGE, s, GIMRParamSTRANGE) +
                myLEP2GIMR.sigma_q_LEP2_GIMR(QCD::BOTTOM, s, GIMRParamBOTTOM);
   
   dsigma_b = myLEP2GIMR.sigma_q_LEP2_GIMR(QCD::BOTTOM, s, GIMRParamCHARM);

   R_bottom = trueSM.LEP2Rbottom(s) * (1 - dsigma_had / sigma_hSM) + dsigma_b / sigma_hSM;
   
   return R_bottom;
}


