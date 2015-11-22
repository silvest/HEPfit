/*
 * Copyright (C) 2013 HEPfit Collaboration
 * All rights reserved.
 *
 * For the licensing terms see doc/COPYING.
 */

#include "NPZbbbarLinearized.h"
#include <stdexcept>


const std::string NPZbbbarLinearized::ZbbbarVAVars[NZbbbarVars]
        = {"deltaGVb", "deltaGAb"};

const std::string NPZbbbarLinearized::ZbbbarLRVars[NZbbbarVars]
        = {"deltaGLb", "deltaGRb"};

NPZbbbarLinearized::NPZbbbarLinearized(const bool FlagNPZbbbarLR_in)
: NPbase(), FlagNPZbbbarLR(FlagNPZbbbarLR_in)
{
    ModelParamMap.insert(std::pair<std::string, boost::reference_wrapper<const double> >("deltaGVb", boost::cref(myDeltaGVb)));
    ModelParamMap.insert(std::pair<std::string, boost::reference_wrapper<const double> >("deltaGAb", boost::cref(myDeltaGAb)));
    ModelParamMap.insert(std::pair<std::string, boost::reference_wrapper<const double> >("deltaGLb", boost::cref(myDeltaGLb)));
    ModelParamMap.insert(std::pair<std::string, boost::reference_wrapper<const double> >("deltaGRb", boost::cref(myDeltaGRb)));
}

bool NPZbbbarLinearized::PostUpdate()
{
    if (!NPbase::PostUpdate()) return (false);

    if (FlagNPZbbbarLR) {
        myDeltaGVb = myDeltaGLb + myDeltaGRb;
        myDeltaGAb = myDeltaGLb - myDeltaGRb;
    } else {
        myDeltaGLb = (myDeltaGVb + myDeltaGAb) / 2.0;
        myDeltaGRb = (myDeltaGVb - myDeltaGAb) / 2.0;
    }

    return (true);
}

void NPZbbbarLinearized::setParameter(const std::string name, const double& value)
{
    if (FlagNPZbbbarLR) {
        if (name.compare("deltaGLb") == 0)
            myDeltaGLb = value;
        else if (name.compare("deltaGRb") == 0)
            myDeltaGRb = value;
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

bool NPZbbbarLinearized::CheckParameters(const std::map<std::string, double>& DPars)
{
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

double NPZbbbarLinearized::deltaGV_f(const Particle f) const
{
    if (f.is("BOTTOM"))
        return myDeltaGVb;
    else
        return 0.0;
}

double NPZbbbarLinearized::deltaGA_f(const Particle f) const
{
    if (f.is("BOTTOM"))
        return myDeltaGAb;
    else
        return 0.0;
}

