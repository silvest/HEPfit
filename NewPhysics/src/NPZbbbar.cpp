/*
 * Copyright (C) 2013 SusyFit Collaboration
 * All rights reserved.
 *
 * For the licensing terms see doc/COPYING.
 */

#include <stdexcept>
#include "NPZbbbar.h"


const std::string NPZbbbar::ZbbbarVAVars[NZbbbarVars]
        = {"deltaGVb", "deltaGAb"};

const std::string NPZbbbar::ZbbbarLRVars[NZbbbarVars]
        = {"deltaGLb", "deltaGRb"};

NPZbbbar::NPZbbbar(const bool FlagNPZbbbarLR_in)
: NPbase(), FlagNPZbbbarLR(FlagNPZbbbarLR_in) {
    if (FlagNPZbbbarLR) {
        ModelParamMap.insert(std::pair<std::string, boost::reference_wrapper<const double> >("deltaGLb", boost::cref(myDeltaGVb)));
        ModelParamMap.insert(std::pair<std::string, boost::reference_wrapper<const double> >("deltaGRb", boost::cref(myDeltaGAb)));
    } else {
        ModelParamMap.insert(std::pair<std::string, boost::reference_wrapper<const double> >("deltaGVb", boost::cref(myDeltaGVb)));
        ModelParamMap.insert(std::pair<std::string, boost::reference_wrapper<const double> >("deltaGAb", boost::cref(myDeltaGAb)));
    }
}

bool NPZbbbar::Update(const std::map<std::string, double>& DPars) {
    for (std::map<std::string, double>::const_iterator it = DPars.begin(); it != DPars.end(); it++)
        setParameter(it->first, it->second);
    if (!NPbase::Update(DPars)) return (false);
    return (true);
}

void NPZbbbar::setParameter(const std::string name, const double& value) {
    if (FlagNPZbbbarLR) {
        if (name.compare("deltaGLb") == 0)
            myDeltaGVb = value;
        else if (name.compare("deltaGRb") == 0)
            myDeltaGAb = value;
        else
            NPbase::setParameter(name, value);
    } else {
        if (name.compare("deltaGVb") == 0)
            myDeltaGVb = value;
        else if (name.compare("deltaGAb") == 0)
            myDeltaGAb = value;
        else
            NPbase::setParameter(name, value);
    }
}

bool NPZbbbar::CheckParameters(const std::map<std::string, double>& DPars) {
    for (int i = 0; i < NZbbbarVars; i++) {
        if (FlagNPZbbbarLR) {
            if (DPars.find(ZbbbarLRVars[i]) == DPars.end()) {
                std::cout << "ERROR: Missing mandatory NPZbbbarLR parameter "
                        << ZbbbarLRVars[i] << std::endl;
                return false;
            }
        } else {
            if (DPars.find(ZbbbarVAVars[i]) == DPars.end()) {
                std::cout << "ERROR: Missing mandatory NPZbbbar parameter "
                        << ZbbbarVAVars[i] << std::endl;
                return false;
            }
        }
    }
    return (NPbase::CheckParameters(DPars));
}

////////////////////////////////////////////////////////////////////////

complex NPZbbbar::rhoZ_f(const Particle p) const {
    if (p.is("BOTTOM"))
        return ( gA_f(p) * gA_f(p)
            / p.getIsospin() / p.getIsospin());
    else
        return trueSM.rhoZ_f(p);
}

complex NPZbbbar::kappaZ_f(const Particle p) const {
    if (p.is("BOTTOM"))
        return ( (1.0 - gV_f(p) / gA_f(p))
            / (4.0 * fabs(p.getCharge()) * sW2()));
    else
        return trueSM.kappaZ_f(p);
}

complex NPZbbbar::gV_f(const Particle p) const {
    if (p.is("BOTTOM"))
        return ( trueSM.gV_f(p) + deltaGV_f(p));
    else
        return trueSM.gV_f(p);
}

complex NPZbbbar::gA_f(const Particle p) const {
    if (p.is("BOTTOM"))
        return ( trueSM.gA_f(p) + deltaGA_f(p));
    else
        return trueSM.gA_f(p);
}

////////////////////////////////////////////////////////////////////////

double NPZbbbar::deltaGV_f(const Particle p) const {
    if (p.is("BOTTOM"))
        if (FlagNPZbbbarLR)
            // delta g_L^b + delta g_R^b
            return ( myDeltaGVb + NPbase::deltaGV_f(p)
                + myDeltaGAb + NPbase::deltaGA_f(p));
        else
            return ( myDeltaGVb + NPbase::deltaGV_f(p));
    else
        return NPbase::deltaGV_f(p);
}

double NPZbbbar::deltaGA_f(const Particle p) const {
    if (p.is("BOTTOM"))
        if (FlagNPZbbbarLR)
            // delta g_L^b + delta g_R^b
            return ( myDeltaGVb + NPbase::deltaGV_f(p)
                - myDeltaGAb - NPbase::deltaGA_f(p));
        else
            return ( myDeltaGAb + NPbase::deltaGA_f(p));
    else
        return NPbase::deltaGA_f(p);
}

double NPZbbbar::Mw() const {
    return (StandardModel::Mw());
}

double NPZbbbar::GammaW() const {
    return (StandardModel::GammaW());
}

double NPZbbbar::Gamma_Z() const {
    return (StandardModel::Gamma_Z());
}

double NPZbbbar::sigma0_had() const {
    return (StandardModel::sigma0_had());
}

double NPZbbbar::sin2thetaEff(const Particle p) const {
    return (StandardModel::sin2thetaEff(p));
}

double NPZbbbar::A_f(const Particle p) const {
    return (StandardModel::A_f(p));
}

double NPZbbbar::AFB(const Particle p) const {
    return (StandardModel::AFB(p));
}

double NPZbbbar::R0_f(const Particle p) const {
    return (StandardModel::R0_f(p));
}

