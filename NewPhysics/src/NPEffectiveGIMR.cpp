/*
 * Copyright (C) 2014 SusyFit Collaboration
 * All rights reserved.
 *
 * For the licensing terms see doc/COPYING.
 */

#include "NPEffectiveGIMR.h"

const std::string NPEffectiveGIMR::NPEffectiveGIMRVars[NNPEffectiveGIMRVars]
        = {"CW", "CHG", "CHW", "CHB", "CHWB", "CHD", "CHbox", "CH",
    "CHL1_11", "CHL1_12", "CHL1_13", "CHL1_22", "CHL1_23", "CHL1_33",
    "CHL3_11", "CHL3_12", "CHL3_13", "CHL3_22", "CHL3_23", "CHL3_33",
    "CHe_11", "CHe_12", "CHe_13", "CHe_22", "CHe_23", "CHe_33",
    "CHQ1_11", "CHQ1_12", "CHQ1_13", "CHQ1_22", "CHQ1_23", "CHQ1_33",
    "CHQ3_11", "CHQ3_12", "CHQ3_13", "CHQ3_22", "CHQ3_23", "CHQ3_33",
    "CHu_11", "CHu_12", "CHu_13", "CHu_22", "CHu_23", "CHu_33",
    "CHd_11", "CHd_12", "CHd_13", "CHd_22", "CHd_23", "CHd_33",
    "CeH_11r", "CeH_12r", "CeH_13r", "CeH_22r", "CeH_23r", "CeH_33r",
    "CeH_11i", "CeH_12i", "CeH_13i", "CeH_22i", "CeH_23i", "CeH_33i",
    "CuH_11r", "CuH_12r", "CuH_13r", "CuH_22r", "CuH_23r", "CuH_33r",
    "CuH_11i", "CuH_12i", "CuH_13i", "CuH_22i", "CuH_23i", "CuH_33i",
    "CdH_11r", "CdH_12r", "CdH_13r", "CdH_22r", "CdH_23r", "CdH_33r",
    "CdH_11i", "CdH_12i", "CdH_13i", "CdH_22i", "CdH_23i", "CdH_33i",
    "CLL_1221", "Lambda_NP"};

const std::string NPEffectiveGIMR::NPEffectiveGIMRVars_LFU_QFU[NNPEffectiveGIMRVars_LFU_QFU]
        = {"CW", "CHG", "CHW", "CHB", "CHWB", "CHD", "CHbox", "CH",
    "CHL1", "CHL3", "CHe", "CHQ1", "CHQ3", "CHu", "CHd",
    "CeH_r", "CeH_i", "CuH_r", "CuH_i", "CdH_r", "CdH_i", "CLL",
    "Lambda_NP"};

NPEffectiveGIMR::NPEffectiveGIMR(const bool FlagLeptonUniversal_in, const bool FlagQuarkUniversal_in)
: NPbase(), FlagLeptonUniversal(FlagLeptonUniversal_in), FlagQuarkUniversal(FlagQuarkUniversal_in)
{
    if ((!FlagLeptonUniversal && !FlagQuarkUniversal)
            || (!FlagLeptonUniversal && FlagQuarkUniversal)
            || (FlagLeptonUniversal && !FlagQuarkUniversal))
        throw std::runtime_error("Invalid arguments for NPEffectiveGIMR::NPEffectiveGIMR()");

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
    } else {
        ModelParamMap.insert(std::pair<std::string, boost::reference_wrapper<const double> >("CHL1_11", boost::cref(CHL1_11)));
        ModelParamMap.insert(std::pair<std::string, boost::reference_wrapper<const double> >("CHL1_12", boost::cref(CHL1_12)));
        ModelParamMap.insert(std::pair<std::string, boost::reference_wrapper<const double> >("CHL1_13", boost::cref(CHL1_13)));
        ModelParamMap.insert(std::pair<std::string, boost::reference_wrapper<const double> >("CHL1_22", boost::cref(CHL1_22)));
        ModelParamMap.insert(std::pair<std::string, boost::reference_wrapper<const double> >("CHL1_23", boost::cref(CHL1_23)));
        ModelParamMap.insert(std::pair<std::string, boost::reference_wrapper<const double> >("CHL1_33", boost::cref(CHL1_33)));
        ModelParamMap.insert(std::pair<std::string, boost::reference_wrapper<const double> >("CHL3_11", boost::cref(CHL3_11)));
        ModelParamMap.insert(std::pair<std::string, boost::reference_wrapper<const double> >("CHL3_12", boost::cref(CHL3_12)));
        ModelParamMap.insert(std::pair<std::string, boost::reference_wrapper<const double> >("CHL3_13", boost::cref(CHL3_13)));
        ModelParamMap.insert(std::pair<std::string, boost::reference_wrapper<const double> >("CHL3_22", boost::cref(CHL3_22)));
        ModelParamMap.insert(std::pair<std::string, boost::reference_wrapper<const double> >("CHL3_23", boost::cref(CHL3_23)));
        ModelParamMap.insert(std::pair<std::string, boost::reference_wrapper<const double> >("CHL3_33", boost::cref(CHL3_33)));
        ModelParamMap.insert(std::pair<std::string, boost::reference_wrapper<const double> >("CHe_11", boost::cref(CHe_11)));
        ModelParamMap.insert(std::pair<std::string, boost::reference_wrapper<const double> >("CHe_12", boost::cref(CHe_12)));
        ModelParamMap.insert(std::pair<std::string, boost::reference_wrapper<const double> >("CHe_13", boost::cref(CHe_13)));
        ModelParamMap.insert(std::pair<std::string, boost::reference_wrapper<const double> >("CHe_22", boost::cref(CHe_22)));
        ModelParamMap.insert(std::pair<std::string, boost::reference_wrapper<const double> >("CHe_23", boost::cref(CHe_23)));
        ModelParamMap.insert(std::pair<std::string, boost::reference_wrapper<const double> >("CHe_33", boost::cref(CHe_33)));
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
        ModelParamMap.insert(std::pair<std::string, boost::reference_wrapper<const double> >("CLL_1221", boost::cref(CLL_1221)));
    }
    if (FlagQuarkUniversal) {
        ModelParamMap.insert(std::pair<std::string, boost::reference_wrapper<const double> >("CHQ1", boost::cref(CHQ1_11)));
        ModelParamMap.insert(std::pair<std::string, boost::reference_wrapper<const double> >("CHQ3", boost::cref(CHQ3_11)));
        ModelParamMap.insert(std::pair<std::string, boost::reference_wrapper<const double> >("CHu", boost::cref(CHu_11)));
        ModelParamMap.insert(std::pair<std::string, boost::reference_wrapper<const double> >("CHd", boost::cref(CHd_11)));
        ModelParamMap.insert(std::pair<std::string, boost::reference_wrapper<const double> >("CuH_r", boost::cref(CuH_11r)));
        ModelParamMap.insert(std::pair<std::string, boost::reference_wrapper<const double> >("CuH_i", boost::cref(CuH_11i)));
        ModelParamMap.insert(std::pair<std::string, boost::reference_wrapper<const double> >("CdH_r", boost::cref(CdH_11r)));
        ModelParamMap.insert(std::pair<std::string, boost::reference_wrapper<const double> >("CdH_i", boost::cref(CdH_11i)));
    } else {
        ModelParamMap.insert(std::pair<std::string, boost::reference_wrapper<const double> >("CHQ1_11", boost::cref(CHQ1_11)));
        ModelParamMap.insert(std::pair<std::string, boost::reference_wrapper<const double> >("CHQ1_12", boost::cref(CHQ1_12)));
        ModelParamMap.insert(std::pair<std::string, boost::reference_wrapper<const double> >("CHQ1_13", boost::cref(CHQ1_13)));
        ModelParamMap.insert(std::pair<std::string, boost::reference_wrapper<const double> >("CHQ1_22", boost::cref(CHQ1_22)));
        ModelParamMap.insert(std::pair<std::string, boost::reference_wrapper<const double> >("CHQ1_23", boost::cref(CHQ1_23)));
        ModelParamMap.insert(std::pair<std::string, boost::reference_wrapper<const double> >("CHQ1_33", boost::cref(CHQ1_33)));
        ModelParamMap.insert(std::pair<std::string, boost::reference_wrapper<const double> >("CHQ3_11", boost::cref(CHQ3_11)));
        ModelParamMap.insert(std::pair<std::string, boost::reference_wrapper<const double> >("CHQ3_12", boost::cref(CHQ3_12)));
        ModelParamMap.insert(std::pair<std::string, boost::reference_wrapper<const double> >("CHQ3_13", boost::cref(CHQ3_13)));
        ModelParamMap.insert(std::pair<std::string, boost::reference_wrapper<const double> >("CHQ3_22", boost::cref(CHQ3_22)));
        ModelParamMap.insert(std::pair<std::string, boost::reference_wrapper<const double> >("CHQ3_23", boost::cref(CHQ3_23)));
        ModelParamMap.insert(std::pair<std::string, boost::reference_wrapper<const double> >("CHQ3_33", boost::cref(CHQ3_33)));
        ModelParamMap.insert(std::pair<std::string, boost::reference_wrapper<const double> >("CHu_11", boost::cref(CHu_11)));
        ModelParamMap.insert(std::pair<std::string, boost::reference_wrapper<const double> >("CHu_12", boost::cref(CHu_12)));
        ModelParamMap.insert(std::pair<std::string, boost::reference_wrapper<const double> >("CHu_13", boost::cref(CHu_13)));
        ModelParamMap.insert(std::pair<std::string, boost::reference_wrapper<const double> >("CHu_22", boost::cref(CHu_22)));
        ModelParamMap.insert(std::pair<std::string, boost::reference_wrapper<const double> >("CHu_23", boost::cref(CHu_23)));
        ModelParamMap.insert(std::pair<std::string, boost::reference_wrapper<const double> >("CHu_33", boost::cref(CHu_33)));
        ModelParamMap.insert(std::pair<std::string, boost::reference_wrapper<const double> >("CHd_11", boost::cref(CHd_11)));
        ModelParamMap.insert(std::pair<std::string, boost::reference_wrapper<const double> >("CHd_12", boost::cref(CHd_12)));
        ModelParamMap.insert(std::pair<std::string, boost::reference_wrapper<const double> >("CHd_13", boost::cref(CHd_13)));
        ModelParamMap.insert(std::pair<std::string, boost::reference_wrapper<const double> >("CHd_22", boost::cref(CHd_22)));
        ModelParamMap.insert(std::pair<std::string, boost::reference_wrapper<const double> >("CHd_23", boost::cref(CHd_23)));
        ModelParamMap.insert(std::pair<std::string, boost::reference_wrapper<const double> >("CHd_33", boost::cref(CHd_33)));

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
    }
    ModelParamMap.insert(std::pair<std::string, boost::reference_wrapper<const double> >("Lambda_NP", boost::cref(Lambda_NP)));
}

bool NPEffectiveGIMR::PostUpdate()
{
    if (!NPbase::PostUpdate()) return (false);

    LambdaNP2 = Lambda_NP * Lambda_NP;
    v2_over_LambdaNP2 = v() * v() / LambdaNP2;
    cW_tree = Mw_tree() / Mz;
    cW2_tree = cW_tree * cW_tree;
    sW2_tree = 1.0 - cW2_tree;
    sW_tree = sqrt(sW2_tree);

    return (true);
}

void NPEffectiveGIMR::setParameter(const std::string name, const double& value)
{
    if (name.compare("CW") == 0)
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
    else if (name.compare("CHL1_12") == 0)
        CHL1_12 = value;
    else if (name.compare("CHL1_13") == 0)
        CHL1_13 = value;
    else if (name.compare("CHL1_22") == 0)
        CHL1_22 = value;
    else if (name.compare("CHL1_23") == 0)
        CHL1_23 = value;
    else if (name.compare("CHL1_33") == 0)
        CHL1_33 = value;
    else if (name.compare("CHL1") == 0) {
        CHL1_11 = value;
        CHL1_12 = 0.0;
        CHL1_13 = 0.0;
        CHL1_22 = value;
        CHL1_23 = 0.0;
        CHL1_33 = value;
    } else if (name.compare("CHL3_11") == 0)
        CHL3_11 = value;
    else if (name.compare("CHL3_12") == 0)
        CHL3_12 = value;
    else if (name.compare("CHL3_13") == 0)
        CHL3_13 = value;
    else if (name.compare("CHL3_22") == 0)
        CHL3_22 = value;
    else if (name.compare("CHL3_23") == 0)
        CHL3_23 = value;
    else if (name.compare("CHL3_33") == 0)
        CHL3_33 = value;
    else if (name.compare("CHL3") == 0) {
        CHL3_11 = value;
        CHL3_12 = 0.0;
        CHL3_13 = 0.0;
        CHL3_22 = value;
        CHL3_23 = 0.0;
        CHL3_33 = value;
    } else if (name.compare("CHe_11") == 0)
        CHe_11 = value;
    else if (name.compare("CHe_12") == 0)
        CHe_12 = value;
    else if (name.compare("CHe_13") == 0)
        CHe_13 = value;
    else if (name.compare("CHe_22") == 0)
        CHe_22 = value;
    else if (name.compare("CHe_23") == 0)
        CHe_23 = value;
    else if (name.compare("CHe_33") == 0)
        CHe_33 = value;
    else if (name.compare("CHe") == 0) {
        CHe_11 = value;
        CHe_12 = 0.0;
        CHe_13 = 0.0;
        CHe_22 = value;
        CHe_23 = 0.0;
        CHe_33 = value;
    } else if (name.compare("CHQ1_11") == 0)
        CHQ1_11 = value;
    else if (name.compare("CHQ1_12") == 0)
        CHQ1_12 = value;
    else if (name.compare("CHQ1_13") == 0)
        CHQ1_13 = value;
    else if (name.compare("CHQ1_22") == 0)
        CHQ1_22 = value;
    else if (name.compare("CHQ1_23") == 0)
        CHQ1_23 = value;
    else if (name.compare("CHQ1_33") == 0)
        CHQ1_33 = value;
    else if (name.compare("CHQ1") == 0) {
        CHQ1_11 = value;
        CHQ1_12 = 0.0;
        CHQ1_13 = 0.0;
        CHQ1_22 = value;
        CHQ1_23 = 0.0;
        CHQ1_33 = value;
    } else if (name.compare("CHQ3_11") == 0)
        CHQ3_11 = value;
    else if (name.compare("CHQ3_12") == 0)
        CHQ3_12 = value;
    else if (name.compare("CHQ3_13") == 0)
        CHQ3_13 = value;
    else if (name.compare("CHQ3_22") == 0)
        CHQ3_22 = value;
    else if (name.compare("CHQ3_23") == 0)
        CHQ3_23 = value;
    else if (name.compare("CHQ3_33") == 0)
        CHQ3_33 = value;
    else if (name.compare("CHQ3") == 0) {
        CHQ3_11 = value;
        CHQ3_12 = 0.0;
        CHQ3_13 = 0.0;
        CHQ3_22 = value;
        CHQ3_23 = 0.0;
        CHQ3_33 = value;
    } else if (name.compare("CHu_11") == 0)
        CHu_11 = value;
    else if (name.compare("CHu_12") == 0)
        CHu_12 = value;
    else if (name.compare("CHu_13") == 0)
        CHu_13 = value;
    else if (name.compare("CHu_22") == 0)
        CHu_22 = value;
    else if (name.compare("CHu_23") == 0)
        CHu_23 = value;
    else if (name.compare("CHu_33") == 0)
        CHu_33 = value;
    else if (name.compare("CHu") == 0) {
        CHu_11 = value;
        CHu_12 = 0.0;
        CHu_13 = 0.0;
        CHu_22 = value;
        CHu_23 = 0.0;
        CHu_33 = value;
    } else if (name.compare("CHd_11") == 0)
        CHd_11 = value;
    else if (name.compare("CHd_12") == 0)
        CHd_12 = value;
    else if (name.compare("CHd_13") == 0)
        CHd_13 = value;
    else if (name.compare("CHd_22") == 0)
        CHd_22 = value;
    else if (name.compare("CHd_23") == 0)
        CHd_23 = value;
    else if (name.compare("CHd_33") == 0)
        CHd_33 = value;
    else if (name.compare("CHd") == 0) {
        CHd_11 = value;
        CHd_12 = 0.0;
        CHd_13 = 0.0;
        CHd_22 = value;
        CHd_23 = 0.0;
        CHd_33 = value;
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
    } else if (name.compare("CLL_1221") == 0) {
        CLL_1221 = value;
        CLL_2112 = value;
    } else if (name.compare("CLL") == 0) {
        CLL_1221 = value;
        CLL_2112 = value;
    } else if (name.compare("Lambda_NP") == 0)
        Lambda_NP = value;
    else
        NPbase::setParameter(name, value);
}

bool NPEffectiveGIMR::CheckParameters(const std::map<std::string, double>& DPars)
{
    if (FlagLeptonUniversal && FlagQuarkUniversal) {
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


////////////////////////////////////////////////////////////////////////

double NPEffectiveGIMR::DeltaGF() const
{
    return ((CHL3_11 + CHL3_22 - 0.5 * (CLL_1221 + CLL_2112)) * v2_over_LambdaNP2);
}

double NPEffectiveGIMR::obliqueS() const
{
    return (4.0 * sW2_tree * cW_tree * CHWB / alphaMz() * v2_over_LambdaNP2);
}

double NPEffectiveGIMR::obliqueT() const
{
    return (-CHD / 2.0 / alphaMz() * v2_over_LambdaNP2);
}

double NPEffectiveGIMR::obliqueU() const
{
    return 0.0;
}

double NPEffectiveGIMR::Mw() const
{
    return (trueSM.Mw() - Mw_tree() / 4.0 / (cW2_tree - sW2_tree)
            *(4.0 * sW_tree * cW_tree * CHWB * v2_over_LambdaNP2
            + cW2_tree * CHD * v2_over_LambdaNP2
            + 2.0 * sW2_tree * DeltaGF()));
}

double NPEffectiveGIMR::GammaW() const
{
    double G0 = GF * pow(Mw(), 3.0) / 6.0 / sqrt(2.0) / M_PI;
    double GammaW_tree = (3.0 + 2.0 * Nc) * G0;

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
    double CHF1, CHF3;
    if (p.is("NEUTRINO_1") || p.is("ELECTRON")) {
        CHF1 = CHL1_11;
        CHF3 = CHL3_11;
    } else if (p.is("NEUTRINO_2") || p.is("MU")) {
        CHF1 = CHL1_22;
        CHF3 = CHL3_22;
    } else if (p.is("NEUTRINO_3") || p.is("TAU")) {
        CHF1 = CHL1_33;
        CHF3 = CHL3_33;
    } else if (p.is("UP") || p.is("DOWN")) {
        CHF1 = CHQ1_11;
        CHF3 = CHQ3_11;
    } else if (p.is("CHARM") || p.is("STRANGE")) {
        CHF1 = CHQ1_22;
        CHF3 = CHQ3_22;
    } else if (p.is("TOP")) {
        return 0.0;
    } else if (p.is("BOTTOM")) {
        CHF1 = CHQ1_33;
        CHF3 = CHQ3_33;
    } else
        throw std::runtime_error("Error in NPEffectiveGIMR::deltaGL_f()");

    double NPdirect = -I3p / 4.0 * (CHD * v2_over_LambdaNP2 + 2.0 * DeltaGF())
            - Qp * sW2_tree / 4.0 / (cW2_tree - sW2_tree)
            *((4.0 * cW_tree / sW_tree * CHWB + CHD) * v2_over_LambdaNP2 + 2.0 * DeltaGF());
    double NPindirect = -0.5 * (CHF1 - 2.0 * I3p * CHF3) * v2_over_LambdaNP2;
    return (NPdirect + NPindirect);
}

double NPEffectiveGIMR::deltaGR_f(const Particle p) const
{
    double Qp = p.getCharge();
    double CHf;
    if (p.is("NEUTRINO_1") || p.is("NEUTRINO_2") || p.is("NEUTRINO_3"))
        return 0.0;
    else if (p.is("ELECTRON"))
        CHf = CHe_11;
    else if (p.is("MU"))
        CHf = CHe_22;
    else if (p.is("TAU"))
        CHf = CHe_33;
    else if (p.is("UP"))
        CHf = CHu_11;
    else if (p.is("CHARM"))
        CHf = CHu_22;
    else if (p.is("TOP"))
        return 0.0;
    else if (p.is("DOWN"))
        CHf = CHd_11;
    else if (p.is("STRANGE"))
        CHf = CHd_22;
    else if (p.is("BOTTOM"))
        CHf = CHd_33;
    else
        throw std::runtime_error("Error in NPEffectiveGIMR::deltaGR_f()");

    double NPdirect = -Qp * sW2_tree / 4.0 / (cW2_tree - sW2_tree)
            *((4.0 * cW_tree / sW_tree * CHWB + CHD) * v2_over_LambdaNP2 + 2.0 * DeltaGF());
    double NPindirect = -0.5 * CHf*v2_over_LambdaNP2;
    return (NPdirect + NPindirect);
}


////////////////////////////////////////////////////////////////////////

double NPEffectiveGIMR::muggH(const double sqrt_s) const
{
    return 1.0;
}

double NPEffectiveGIMR::muVBF(const double sqrt_s) const
{
    return 1.0;
}

double NPEffectiveGIMR::muWH(const double sqrt_s) const
{
    return 1.0;
}

double NPEffectiveGIMR::muZH(const double sqrt_s) const
{
    return 1.0;
}

double NPEffectiveGIMR::muVH(const double sqrt_s) const
{
    return 1.0;
}

double NPEffectiveGIMR::muttH(const double sqrt_s) const
{
    return 1.0;
}

double NPEffectiveGIMR::BrHggRatio() const
{
    return (1.0
            + 252512.0 * CdH_33r / LambdaNP2
            + 198682.0 * CuH_22r / LambdaNP2
            - 129195.0 * CuH_33r / LambdaNP2
            + 3.57e7 * CHG / LambdaNP2
            + 30312.0 * (4.0 * CHbox / LambdaNP2 - CHD / LambdaNP2));
}

double NPEffectiveGIMR::BrHWWRatio() const
{
    return (1.0 
            - 89706.0 * CHW / LambdaNP2
            + 30312.0 * (4.0 * CHbox / LambdaNP2 - CHD / LambdaNP2));
}

double NPEffectiveGIMR::BrHZZRatio() const
{
    return (1.0
            - 14429.0 * CHB / LambdaNP2
            - 45453.0 * CHW / LambdaNP2
            - 27649.0 * CHWB / LambdaNP2
            + 30312.0 * (4.0 * CHbox / LambdaNP2 + CHD / LambdaNP2));
}

double NPEffectiveGIMR::BrHZgaRatio() const
{
    return (1.0
            - 7575.0 * CdH_33r / LambdaNP2
            - 4139.0 * CuH_22r / LambdaNP2
            + 7313.0 * CuH_33r / LambdaNP2
            - 535.0 * CeH_33r / LambdaNP2
            + 1.52e7 * CHB / LambdaNP2
            - 1.32e7 * CHW / LambdaNP2
            + 9.0e6 * CHWB / LambdaNP2
            + 30312.0 * (4.0 * CHbox / LambdaNP2 - CHD / LambdaNP2));
}

double NPEffectiveGIMR::BrHgagaRatio() const
{
    return (1.0
            - 17676.0 * CdH_33r / LambdaNP2
            - 24835.0 * CuH_22r / LambdaNP2
            + 32908.0 * CuH_33r / LambdaNP2
            - 41583.0 * CeH_33r / LambdaNP2
            + 4.82e7 * CHB / LambdaNP2
            - 1.27e7 * CHW / LambdaNP2
            + 2.5e7 * CHWB / LambdaNP2
            + 30312.0 * (4.0 * CHbox / LambdaNP2 - CHD / LambdaNP2));
}

double NPEffectiveGIMR::BrHtautauRatio() const
{
    return (1.0
            - 1.19e7 * CeH_33r / LambdaNP2
            + 30312.0 * (4.0 * CHbox / LambdaNP2 - CHD / LambdaNP2));
}

double NPEffectiveGIMR::BrHccRatio() const
{
    return (1.0
            - 1.64e7 * CuH_22r / LambdaNP2
            - 914.0 * CuH_33r / LambdaNP2
            + 30312.0 * (4.0 * CHbox / LambdaNP2 - CHD / LambdaNP2));
}

double NPEffectiveGIMR::BrHbbRatio() const
{
    return (1.0
            - 5.03e6 * CdH_33r / LambdaNP2
            - 518.0 * CuH_33r / LambdaNP2
            + 30312.0 * (4.0 * CHbox / LambdaNP2 - CHD / LambdaNP2));
}

