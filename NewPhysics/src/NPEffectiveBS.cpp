/* 
 * Copyright (C) 2013 HEPfit Collaboration
 * All rights reserved.
 *
 * For the licensing terms see doc/COPYING.
 */

#include "NPEffectiveBS.h"
#include <stdexcept>

const std::string NPEffectiveBS::NPEffectiveBSVars[NNPEffectiveBSVars]
        = {"cWB", "cH", "cL1L1", "cL1L2", "cL1L3", "cL2L2", "cL2L3", "cL3L3",
    "cHL1p", "cHL2p", "cHL3p", "cHQ1p", "cHQ2p", "cHQ3p",
    "cHL1", "cHL2", "cHL3", "cHQ1", "cHQ2", "cHQ3", "cHE1", "cHE2", "cHE3",
    "cHU1", "cHU2", "cHU3", "cHD1", "cHD2", "cHD3", "Lambda_NP"};

const std::string NPEffectiveBS::NPEffectiveBSVars_LFU[NNPEffectiveBSVars_LFU]
        = {"cWB_NP", "cH_NP", "cLL_NP", "cHLp_NP", "cHL_NP",
    "cHQ1p_NP", "cHQ2p_NP", "cHQ3p_NP", "cHQ1_NP", "cHQ2_NP", "cHQ3_NP",
    "cHU1_NP", "cHU2_NP", "cHU3_NP", "cHD1_NP", "cHD2_NP", "cHD3_NP", "cHE_NP", "Lambda_NP"};

const std::string NPEffectiveBS::NPEffectiveBSVars_QFU[NNPEffectiveBSVars_QFU]
        = {"cWB_NP", "cH_NP", "cL1L1", "cL1L2", "cL1L3", "cL2L2", "cL2L3", "cL3L3",
    "cHL1p", "cHL2p", "cHL3p", "cHQp_NP", "cHL1", "cHL2", "cHL3",
    "cHQ_NP", "cHE1", "cHE2", "cHE3", "cHU_NP", "cHD_NP", "Lambda_NP"};

const std::string NPEffectiveBS::NPEffectiveBSVars_LFU_QFU[NNPEffectiveBSVars_LFU_QFU]
        = {"cWB_NP", "cH_NP", "cLL_NP", "cHLp_NP", "cHQp_NP",
    "cHL_NP", "cHQ_NP", "cHE_NP", "cHU_NP", "cHD_NP", "Lambda_NP"};

NPEffectiveBS::NPEffectiveBS(const bool FlagLeptonUniversal_in, const bool FlagQuarkUniversal_in)
: NPbase(), FlagLeptonUniversal(FlagLeptonUniversal_in), FlagQuarkUniversal(FlagQuarkUniversal_in)
{
    ModelParamMap.insert(std::pair<std::string, boost::reference_wrapper<const double> >("cWB_NP", boost::cref(cWB)));
    ModelParamMap.insert(std::pair<std::string, boost::reference_wrapper<const double> >("cH_NP", boost::cref(cH)));
    if (FlagLeptonUniversal) {
        ModelParamMap.insert(std::pair<std::string, boost::reference_wrapper<const double> >("cLL_NP", boost::cref(cL1L1)));
        ModelParamMap.insert(std::pair<std::string, boost::reference_wrapper<const double> >("cHLp_NP", boost::cref(cHL1p)));
        ModelParamMap.insert(std::pair<std::string, boost::reference_wrapper<const double> >("cHL_NP", boost::cref(cHL1)));
        ModelParamMap.insert(std::pair<std::string, boost::reference_wrapper<const double> >("cHE_NP", boost::cref(cHE1)));
    } else {
        ModelParamMap.insert(std::pair<std::string, boost::reference_wrapper<const double> >("cL1L1_NP", boost::cref(cL1L1)));
        ModelParamMap.insert(std::pair<std::string, boost::reference_wrapper<const double> >("cL1L2_NP", boost::cref(cL1L2)));
        ModelParamMap.insert(std::pair<std::string, boost::reference_wrapper<const double> >("cL1L3_NP", boost::cref(cL1L3)));
        ModelParamMap.insert(std::pair<std::string, boost::reference_wrapper<const double> >("cL2L2_NP", boost::cref(cL2L2)));
        ModelParamMap.insert(std::pair<std::string, boost::reference_wrapper<const double> >("cL2L3_NP", boost::cref(cL2L3)));
        ModelParamMap.insert(std::pair<std::string, boost::reference_wrapper<const double> >("cL3L3_NP", boost::cref(cL3L3)));
        ModelParamMap.insert(std::pair<std::string, boost::reference_wrapper<const double> >("cHL1p_NP", boost::cref(cHL1p)));
        ModelParamMap.insert(std::pair<std::string, boost::reference_wrapper<const double> >("cHL2p_NP", boost::cref(cHL2p)));
        ModelParamMap.insert(std::pair<std::string, boost::reference_wrapper<const double> >("cHL3p_NP", boost::cref(cHL3p)));
        ModelParamMap.insert(std::pair<std::string, boost::reference_wrapper<const double> >("cHL1_NP", boost::cref(cHL1)));
        ModelParamMap.insert(std::pair<std::string, boost::reference_wrapper<const double> >("cHL2_NP", boost::cref(cHL2)));
        ModelParamMap.insert(std::pair<std::string, boost::reference_wrapper<const double> >("cHL3_NP", boost::cref(cHL3)));
        ModelParamMap.insert(std::pair<std::string, boost::reference_wrapper<const double> >("cHE1_NP", boost::cref(cHE1)));
        ModelParamMap.insert(std::pair<std::string, boost::reference_wrapper<const double> >("cHE2_NP", boost::cref(cHE2)));
        ModelParamMap.insert(std::pair<std::string, boost::reference_wrapper<const double> >("cHE3_NP", boost::cref(cHE3)));
    }
    if (FlagQuarkUniversal) {
        ModelParamMap.insert(std::pair<std::string, boost::reference_wrapper<const double> >("cHQp_NP", boost::cref(cHQ1p)));
        ModelParamMap.insert(std::pair<std::string, boost::reference_wrapper<const double> >("cHQ_NP", boost::cref(cHQ1)));
        ModelParamMap.insert(std::pair<std::string, boost::reference_wrapper<const double> >("cHU_NP", boost::cref(cHU1)));
        ModelParamMap.insert(std::pair<std::string, boost::reference_wrapper<const double> >("cHD_NP", boost::cref(cHD1)));
    } else {
        ModelParamMap.insert(std::pair<std::string, boost::reference_wrapper<const double> >("cHQ1p_NP", boost::cref(cHQ1p)));
        ModelParamMap.insert(std::pair<std::string, boost::reference_wrapper<const double> >("cHQ2p_NP", boost::cref(cHQ2p)));
        ModelParamMap.insert(std::pair<std::string, boost::reference_wrapper<const double> >("cHQ3p_NP", boost::cref(cHQ3p)));
        ModelParamMap.insert(std::pair<std::string, boost::reference_wrapper<const double> >("cHQ1_NP", boost::cref(cHQ1)));
        ModelParamMap.insert(std::pair<std::string, boost::reference_wrapper<const double> >("cHQ2_NP", boost::cref(cHQ2)));
        ModelParamMap.insert(std::pair<std::string, boost::reference_wrapper<const double> >("cHQ3_NP", boost::cref(cHQ3)));
        ModelParamMap.insert(std::pair<std::string, boost::reference_wrapper<const double> >("cHU1_NP", boost::cref(cHU1)));
        ModelParamMap.insert(std::pair<std::string, boost::reference_wrapper<const double> >("cHU2_NP", boost::cref(cHU2)));
        ModelParamMap.insert(std::pair<std::string, boost::reference_wrapper<const double> >("cHU3_NP", boost::cref(cHU3)));
        ModelParamMap.insert(std::pair<std::string, boost::reference_wrapper<const double> >("cHD1_NP", boost::cref(cHD1)));
        ModelParamMap.insert(std::pair<std::string, boost::reference_wrapper<const double> >("cHD2_NP", boost::cref(cHD2)));
        ModelParamMap.insert(std::pair<std::string, boost::reference_wrapper<const double> >("cHD3_NP", boost::cref(cHD3)));
    }
    ModelParamMap.insert(std::pair<std::string, boost::reference_wrapper<const double> >("Lambda_NP", boost::cref(Lambda_NP)));
}


////////////////////////////////////////////////////////////////////////

void NPEffectiveBS::setParameter(const std::string name, const double& value)
{
    if (name.compare("cWB_NP") == 0)
        cWB = value;
    else if (name.compare("cH_NP") == 0)
        cH = value;
    else if (name.compare("cL1L1_NP") == 0)
        cL1L1 = value;
    else if (name.compare("cL1L2_NP") == 0)
        cL1L2 = value;
    else if (name.compare("cL1L3_NP") == 0)
        cL1L3 = value;
    else if (name.compare("cL2L2_NP") == 0)
        cL2L2 = value;
    else if (name.compare("cL2L3_NP") == 0)
        cL2L3 = value;
    else if (name.compare("cL3L3_NP") == 0)
        cL3L3 = value;
    else if (name.compare("cLL_NP") == 0) {
        cL1L1 = value;
        cL1L2 = value;
        cL1L3 = value;
        cL2L2 = value;
        cL2L3 = value;
        cL3L3 = value;
    } else if (name.compare("cHL1p_NP") == 0)
        cHL1p = value;
    else if (name.compare("cHL2p_NP") == 0)
        cHL2p = value;
    else if (name.compare("cHL3p_NP") == 0)
        cHL3p = value;
    else if (name.compare("cHLp_NP") == 0) {
        cHL1p = value;
        cHL2p = value;
        cHL3p = value;
    } else if (name.compare("cHQ1p_NP") == 0)
        cHQ1p = value;
    else if (name.compare("cHQ2p_NP") == 0)
        cHQ2p = value;
    else if (name.compare("cHQ3p_NP") == 0)
        cHQ3p = value;
    else if (name.compare("cHQp_NP") == 0) {
        cHQ1p = value;
        cHQ2p = value;
        cHQ3p = value;
    } else if (name.compare("cHL1_NP") == 0)
        cHL1 = value;
    else if (name.compare("cHL2_NP") == 0)
        cHL2 = value;
    else if (name.compare("cHL3_NP") == 0)
        cHL3 = value;
    else if (name.compare("cHL_NP") == 0) {
        cHL1 = value;
        cHL2 = value;
        cHL3 = value;
    } else if (name.compare("cHQ1_NP") == 0)
        cHQ1 = value;
    else if (name.compare("cHQ2_NP") == 0)
        cHQ2 = value;
    else if (name.compare("cHQ3_NP") == 0)
        cHQ3 = value;
    else if (name.compare("cHQ_NP") == 0) {
        cHQ1 = value;
        cHQ2 = value;
        cHQ3 = value;
    } else if (name.compare("cHE1_NP") == 0)
        cHE1 = value;
    else if (name.compare("cHE2_NP") == 0)
        cHE2 = value;
    else if (name.compare("cHE3_NP") == 0)
        cHE3 = value;
    else if (name.compare("cHE_NP") == 0) {
        cHE1 = value;
        cHE2 = value;
        cHE3 = value;
    } else if (name.compare("cHU1_NP") == 0)
        cHU1 = value;
    else if (name.compare("cHU2_NP") == 0)
        cHU2 = value;
    else if (name.compare("cHU3_NP") == 0)
        cHU3 = value;
    else if (name.compare("cHU_NP") == 0) {
        cHU1 = value;
        cHU2 = value;
        cHU3 = value;
    } else if (name.compare("cHD1_NP") == 0)
        cHD1 = value;
    else if (name.compare("cHD2_NP") == 0)
        cHD2 = value;
    else if (name.compare("cHD3_NP") == 0)
        cHD3 = value;
    else if (name.compare("cHD_NP") == 0) {
        cHD1 = value;
        cHD2 = value;
        cHD3 = value;
    } else if (name.compare("Lambda_NP") == 0)
        Lambda_NP = value;
    else
        NPbase::setParameter(name, value);
}

bool NPEffectiveBS::CheckParameters(const std::map<std::string, double>& DPars)
{
    if (FlagLeptonUniversal && FlagQuarkUniversal) {
        for (int i = 0; i < NNPEffectiveBSVars_LFU_QFU; i++) {
            if (DPars.find(NPEffectiveBSVars_LFU_QFU[i]) == DPars.end()) {
                std::cout << "ERROR: Missing mandatory NPEffectiveBS_LFU_QFU parameter "
                        << NPEffectiveBSVars_LFU_QFU[i] << std::endl;
                return false;
            }
        }
    } else if (FlagLeptonUniversal && !FlagQuarkUniversal) {
        for (int i = 0; i < NNPEffectiveBSVars_LFU; i++) {
            if (DPars.find(NPEffectiveBSVars_LFU[i]) == DPars.end()) {
                std::cout << "ERROR: Missing mandatory NPEffectiveBS_LFU parameter "
                        << NPEffectiveBSVars_LFU[i] << std::endl;
                return false;
            }
        }
    } else if (!FlagLeptonUniversal && FlagQuarkUniversal) {
        for (int i = 0; i < NNPEffectiveBSVars_QFU; i++) {
            if (DPars.find(NPEffectiveBSVars_QFU[i]) == DPars.end()) {
                std::cout << "ERROR: Missing mandatory NPEffectiveBS_QFU parameter "
                        << NPEffectiveBSVars_QFU[i] << std::endl;
                return false;
            }
        }
    } else if (!FlagLeptonUniversal && !FlagQuarkUniversal) {
        for (int i = 0; i < NNPEffectiveBSVars; i++) {
            if (DPars.find(NPEffectiveBSVars[i]) == DPars.end()) {
                std::cout << "ERROR: Missing mandatory NPEffectiveBS parameter "
                        << NPEffectiveBSVars[i] << std::endl;
                return false;
            }
        }
    } else
        throw std::runtime_error("Error in NPEffectiveBS::CheckParameters()");

    return (NPbase::CheckParameters(DPars));
}



////////////////////////////////////////////////////////////////////////

double NPEffectiveBS::v() const
{
    //return ( sqrt( (1.0 - (cL1L2 - cHL1p - cHL2p)/sqrt(2.0)/GF/LambdaNP/LambdaNP)
    //               /sqrt(2.0)/GF ) );

    /* use the tree-level relation */
    return trueSM.v();
}

double NPEffectiveBS::Mw_tree() const
{
    double GF0 = GF * (1.0 - DeltaGF());
    double tmp = 4.0 * M_PI * ale / sqrt(2.0) / GF0 / Mz / Mz;
    return ( Mz / sqrt(2.0) * sqrt(1.0 + sqrt(1.0 - tmp)));
}

double NPEffectiveBS::DeltaGF() const
{
    double ratio = v() * v() / Lambda_NP / Lambda_NP;

    return ( -(cL1L2 - cHL1p - cHL2p) * ratio);
}


////////////////////////////////////////////////////////////////////////

double NPEffectiveBS::GammaW() const
{
    double Gamma_W = trueSM.GammaW();

    double alpha = alphaMz();
    double c2 = trueSM.cW2();
    double s2 = trueSM.sW2();
    double ratio = v() * v() / Lambda_NP / Lambda_NP;

    Gamma_W *= 1.0 - 3.0 * alpha / 4.0 / (c2 - s2)
            *(obliqueS() - 2.0 * c2 * obliqueT() - (c2 - s2) * obliqueU() / 2.0 / s2)
            - (1.0 + c2) / 2.0 / (c2 - s2) * DeltaGF()
            //+ (cHL1p + cHL2p + cHL3p + cHQ1p + cHQ2p) * ratio; // incorrect
            + 2.0 / 9.0 * (cHL1p + cHL2p + cHL3p + 3.0 * cHQ1p + 3.0 * cHQ2p) * ratio;

    return Gamma_W;
}


////////////////////////////////////////////////////////////////////////     

double NPEffectiveBS::obliqueS() const
{
    double ratio = v() * v() / Lambda_NP / Lambda_NP;
    double sW_SM = sqrt(trueSM.sW2()); /* This has to be the SM value. */
    double cW_SM = sqrt(trueSM.cW2()); /* This has to be the SM value. */

    return ( 4.0 * sW_SM * cW_SM * cWB / alphaMz() * ratio);
}

double NPEffectiveBS::obliqueT() const
{
    double ratio = v() * v() / Lambda_NP / Lambda_NP;

    return ( -cH / 2.0 / alphaMz() * ratio);
}

double NPEffectiveBS::obliqueU() const
{
    return 0.0;
}


////////////////////////////////////////////////////////////////////////

double NPEffectiveBS::deltaGV_f(const Particle f) const
{
    return ( deltaGL_f_tmp(f) + deltaGR_f_tmp(f) + NPbase::deltaGV_f(f));
}

double NPEffectiveBS::deltaGA_f(const Particle f) const
{
    return ( deltaGL_f_tmp(f) - deltaGR_f_tmp(f) + NPbase::deltaGA_f(f));
}


////////////////////////////////////////////////////////////////////////

double NPEffectiveBS::deltaGL_f_tmp(const Particle f) const
{
    double ratio = v() * v() / Lambda_NP / Lambda_NP;
    if (f.is("NEUTRINO_1"))
        return ( (cHL1p - cHL1) / 2.0 * ratio);
    if (f.is("NEUTRINO_2"))
        return ( (cHL2p - cHL2) / 2.0 * ratio);
    if (f.is("NEUTRINO_3"))
        return ( (cHL3p - cHL3) / 2.0 * ratio);
    if (f.is("ELECTRON"))
        return ( -(cHL1p + cHL1) / 2.0 * ratio);
    if (f.is("MU"))
        return ( -(cHL2p + cHL2) / 2.0 * ratio);
    if (f.is("TAU"))
        return ( -(cHL3p + cHL3) / 2.0 * ratio);
    if (f.is("UP"))
        return ( (cHQ1p - cHQ1) / 2.0 * ratio);
    if (f.is("CHARM"))
        return ( (cHQ2p - cHQ2) / 2.0 * ratio);
    if (f.is("TOP"))
        return 0.0;
    if (f.is("DOWN"))
        return ( -(cHQ1p + cHQ1) / 2.0 * ratio);
    if (f.is("STRANGE"))
        return ( -(cHQ2p + cHQ2) / 2.0 * ratio);
    if (f.is("BOTTOM"))
        return ( -(cHQ3p + cHQ3) / 2.0 * ratio);
    throw std::runtime_error("Error in NPEffectiveBS::deltaGL_f()");
}

double NPEffectiveBS::deltaGR_f_tmp(const Particle f) const
{
    double ratio = v() * v() / Lambda_NP / Lambda_NP;
    if (f.is("NEUTRINO_1"))
        return ( 0.);
    if (f.is("NEUTRINO_2"))
        return ( 0.);
    if (f.is("NEUTRINO_3"))
        return ( 0.);
    if (f.is("ELECTRON"))
        return ( -cHE1 / 2.0 * ratio);
    if (f.is("MU"))
        return ( -cHE2 / 2.0 * ratio);
    if (f.is("TAU"))
        return ( -cHE3 / 2.0 * ratio);
    if (f.is("UP"))
        return ( -cHU1 / 2.0 * ratio);
    if (f.is("CHARM"))
        return ( -cHU2 / 2.0 * ratio);
    if (f.is("TOP"))
        return 0.0;
    if (f.is("DOWN"))
        return ( -cHD1 / 2.0 * ratio);
    if (f.is("STRANGE"))
        return ( -cHD2 / 2.0 * ratio);
    if (f.is("BOTTOM"))
        return ( -cHD3 / 2.0 * ratio);
    throw std::runtime_error("Error in NPEffectiveBS::deltaGL_f()");
}


