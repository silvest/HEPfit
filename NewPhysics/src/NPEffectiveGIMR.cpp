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
    "CLL_1221", "CLL_2112", "LambdaNP"};

NPEffectiveGIMR::NPEffectiveGIMR()
: NPbase()
{
}

bool NPEffectiveGIMR::PostUpdate()
{
    if (!NPbase::PostUpdate()) return (false);

    v2_over_LambdaNP2 = v() * v() / LambdaNP / LambdaNP;
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
    else if (name.compare("CHL3_11") == 0)
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
    else if (name.compare("CHe_11") == 0)
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
    else if (name.compare("CHQ1_11") == 0)
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
    else if (name.compare("CHQ3_11") == 0)
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
    else if (name.compare("CHu_11") == 0)
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
    else if (name.compare("CHd_11") == 0)
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
    else if (name.compare("CeH_11r") == 0)
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
    else if (name.compare("CeH_11i") == 0)
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
    else if (name.compare("CuH_11r") == 0)
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
    else if (name.compare("CuH_11i") == 0)
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
    else if (name.compare("CdH_11r") == 0)
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
    else if (name.compare("CdH_11i") == 0)
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
    else if (name.compare("CLL_1221") == 0)
        CLL_1221 = value;
    else if (name.compare("CLL_2112") == 0)
        CLL_2112 = value;
    else if (name.compare("LambdaNP") == 0)
        LambdaNP = value;
    else
        NPbase::setParameter(name, value);
}

bool NPEffectiveGIMR::CheckParameters(const std::map<std::string, double>& DPars)
{
    for (int i = 0; i < NNPEffectiveGIMRVars; i++) {
        if (DPars.find(NPEffectiveGIMRVars[i]) == DPars.end()) {
            std::cout << "ERROR: Missing mandatory NPEffectiveGIMR parameter"
                    << NPEffectiveGIMRVars[i] << std::endl;
            return false;
        }
    }
    return (NPbase::CheckParameters(DPars));
}

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

