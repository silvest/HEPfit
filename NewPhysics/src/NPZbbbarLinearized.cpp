/*
 * Copyright (C) 2013 SusyFit Collaboration
 * All rights reserved.
 *
 * For the licensing terms see doc/COPYING.
 */

#include <stdexcept>
#include "NPZbbbarLinearized.h"


const std::string NPZbbbarLinearized::ZbbbarVAVars[NZbbbarVars]
        = {"deltaGVb", "deltaGAb"};

const std::string NPZbbbarLinearized::ZbbbarLRVars[NZbbbarVars]
        = {"deltaGLb", "deltaGRb"};

NPZbbbarLinearized::NPZbbbarLinearized(const bool FlagNPZbbbarLR_in)
: NPbase(), FlagNPZbbbarLR(FlagNPZbbbarLR_in) {
    if (FlagNPZbbbarLR) {
        ModelParamMap.insert(std::pair<std::string, boost::reference_wrapper<const double> >("deltaGLb", boost::cref(myDeltaGVb)));
        ModelParamMap.insert(std::pair<std::string, boost::reference_wrapper<const double> >("deltaGRb", boost::cref(myDeltaGAb)));
    } else {
        ModelParamMap.insert(std::pair<std::string, boost::reference_wrapper<const double> >("deltaGVb", boost::cref(myDeltaGVb)));
        ModelParamMap.insert(std::pair<std::string, boost::reference_wrapper<const double> >("deltaGAb", boost::cref(myDeltaGAb)));
    }

}

bool NPZbbbarLinearized::Update(const std::map<std::string, double>& DPars) {
    for (std::map<std::string, double>::const_iterator it = DPars.begin(); it != DPars.end(); it++)
        setParameter(it->first, it->second);
    if (!NPbase::Update(DPars)) return (false);
    return (true);
}

void NPZbbbarLinearized::setParameter(const std::string name, const double& value) {
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

bool NPZbbbarLinearized::CheckParameters(const std::map<std::string, double>& DPars) {
    for (int i = 0; i < NZbbbarVars; i++) {
        if (FlagNPZbbbarLR) {
            if (DPars.find(ZbbbarLRVars[i]) == DPars.end()) {
                std::cout << "ERROR: Missing mandatory NPZbbbarLinearizedLR parameter "
                        << ZbbbarLRVars[i] << std::endl;
                return false;
            }
        } else {
            if (DPars.find(ZbbbarVAVars[i]) == DPars.end()) {
                std::cout << "ERROR: Missing mandatory NPZbbbarLinearized parameter "
                        << ZbbbarVAVars[i] << std::endl;
                return false;
            }
        }
    }
    return (NPbase::CheckParameters(DPars));
}

////////////////////////////////////////////////////////////////////////

complex NPZbbbarLinearized::rhoZ_f(const Particle p) const {
    if (p.is("BOTTOM"))
        return ( gA_f(p) * gA_f(p)
            / p.getIsospin() / p.getIsospin());
    else
        return trueSM.rhoZ_f(p);
}

complex NPZbbbarLinearized::kappaZ_f(const Particle p) const {
    if (p.is("BOTTOM"))
        return ( (1.0 - gV_f(p) / gA_f(p))
            / (4.0 * fabs(p.getCharge()) * sW2()));
    else
        return trueSM.kappaZ_f(p);
}

complex NPZbbbarLinearized::gV_f(const Particle p) const {
    if (p.is("BOTTOM"))
        return ( trueSM.gA_f(p)*(1.0 - 4.0 * fabs(p.getCharge())*
            (trueSM.kappaZ_f(p)) * sW2()));
    else
        return trueSM.gV_f(p);
}

complex NPZbbbarLinearized::gA_f(const Particle p) const {
    if (p.is("BOTTOM"))
        return ( sqrt(trueSM.rhoZ_f(p)) * p.getIsospin());
    else
        return trueSM.gA_f(p);
}

////////////////////////////////////////////////////////////////////////

double NPZbbbarLinearized::deltaGV_f(const Particle p) const {
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

double NPZbbbarLinearized::deltaGA_f(const Particle p) const {
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
