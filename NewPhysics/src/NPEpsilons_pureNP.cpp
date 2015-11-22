/*
 * Copyright (C) 2013 HEPfit Collaboration
 * All rights reserved.
 *
 * For the licensing terms see doc/COPYING.
 */

#include "NPEpsilons_pureNP.h"
#include <stdexcept>

const std::string NPEpsilons_pureNP::EPSILONpureNPvars[NEPSILONpureNPvars]
        = {"delEps_1", "delEps_2", "delEps_3", "delEps_b"};

NPEpsilons_pureNP::NPEpsilons_pureNP()
: NPbase()
{
    ModelParamMap.insert(std::pair<std::string, boost::reference_wrapper<const double> >("delEps_1", boost::cref(deltaEps_1)));
    ModelParamMap.insert(std::pair<std::string, boost::reference_wrapper<const double> >("delEps_2", boost::cref(deltaEps_2)));
    ModelParamMap.insert(std::pair<std::string, boost::reference_wrapper<const double> >("delEps_3", boost::cref(deltaEps_3)));
    ModelParamMap.insert(std::pair<std::string, boost::reference_wrapper<const double> >("delEps_b", boost::cref(deltaEps_b)));
}

void NPEpsilons_pureNP::setParameter(const std::string name, const double& value)
{
    if (name.compare("delEps_1") == 0)
        deltaEps_1 = value;
    else if (name.compare("delEps_2") == 0)
        deltaEps_2 = value;
    else if (name.compare("delEps_3") == 0)
        deltaEps_3 = value;
    else if (name.compare("delEps_b") == 0)
        deltaEps_b = value;
    else
        NPbase::setParameter(name, value);
}

bool NPEpsilons_pureNP::CheckParameters(const std::map<std::string, double>& DPars)
{
    for (int i = 0; i < NEPSILONpureNPvars; i++) {
        if (DPars.find(EPSILONpureNPvars[i]) == DPars.end()) {
            std::cout << "ERROR: Missing mandatory NPEpsilons_pureNP parameter "
                    << EPSILONpureNPvars[i] << std::endl;
            return false;
        }
    }
    return (NPbase::CheckParameters(DPars));
}

////////////////////////////////////////////////////////////////////////

double NPEpsilons_pureNP::epsilon1() const
{
    return (trueSM.epsilon1() + deltaEps_1);
}

double NPEpsilons_pureNP::epsilon2() const
{
    return (trueSM.epsilon2() + deltaEps_2);
}

double NPEpsilons_pureNP::epsilon3() const
{
    return (trueSM.epsilon3() + deltaEps_3);
}

double NPEpsilons_pureNP::epsilonb() const
{
    return (trueSM.epsilonb() + deltaEps_b);
}


////////////////////////////////////////////////////////////////////////

double NPEpsilons_pureNP::Mw() const
{
    double Mw_SM = trueSM.Mw();
    double sW2_SM = trueSM.sW2();
    double cW2_SM = trueSM.cW2();
    return ( Mw_SM * (1.0 - (-c02() * deltaEps_1
            + (c02() - s02()) * deltaEps_2
            + 2.0 * s02() * deltaEps_3) / (cW2_SM - sW2_SM) / 2.0));
}

double NPEpsilons_pureNP::GammaW() const
{
    throw std::runtime_error("NPEpsilons_pureNP::GammaW() is not available.");
}


////////////////////////////////////////////////////////////////////////

double NPEpsilons_pureNP::deltaGV_f(const Particle f) const
{
    if (f.is("TOP")) return 0.0;

    /* SM values */
    double gV_SM = trueSM.gV_f(f).real();
    double gA_SM = trueSM.gA_f(f).real();

    if (f.is("BOTTOM"))
        return ( gV_SM * deltaEps_1 / 2.0 + gV_SM * deltaEps_b
            + (gV_SM - gA_SM)*(deltaEps_3 - c02() * deltaEps_1)
            / (c02() - s02())
            - (gV_SM - gA_SM) * deltaEps_b);
    else
        return ( gV_SM * deltaEps_1 / 2.0
            + (gV_SM - gA_SM)*(deltaEps_3 - c02() * deltaEps_1)
            / (c02() - s02()));
}

double NPEpsilons_pureNP::deltaGA_f(const Particle f) const
{
    if (f.is("TOP")) return 0.0;

    if (f.is("BOTTOM"))
        return ( f.getIsospin()*(deltaEps_1 / 2.0 + deltaEps_b));

    return ( f.getIsospin() * deltaEps_1 / 2.0);
}


