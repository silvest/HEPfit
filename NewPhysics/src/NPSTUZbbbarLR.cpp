/*
 * Copyright (C) 2015 HEPfit Collaboration
 *
 *
 * For the licensing terms see doc/COPYING.
 */

#include "NPSTUZbbbarLR.h"
#include <stdexcept>


const std::string NPSTUZbbbarLR::STUZbbbarLRvars[NSTUZbbbarLRvars]
        = {"deltaGLb", "deltaGRb"};

NPSTUZbbbarLR::NPSTUZbbbarLR()
: NPSTU()
{
    ModelParamMap.insert(std::make_pair("deltaGLb", std::cref(myDeltaGLb)));
    ModelParamMap.insert(std::make_pair("deltaGRb", std::cref(myDeltaGRb)));
}

void NPSTUZbbbarLR::setParameter(const std::string name, const double& value)
{
    if (name.compare("deltaGLb") == 0)
        myDeltaGLb = value;
    else if (name.compare("deltaGRb") == 0)
        myDeltaGRb = value;
    else
        NPSTU::setParameter(name, value);
}

bool NPSTUZbbbarLR::CheckParameters(const std::map<std::string, double>& DPars)
{
    for (int i = 0; i < NSTUZbbbarLRvars; i++) {
        if (DPars.find(STUZbbbarLRvars[i]) == DPars.end()) {
            std::cout << "ERROR: Missing mandatory NPSTUZbbbarLR parameter "
                    << STUZbbbarLRvars[i] << std::endl;
            raiseMissingModelParameterCount();
            addMissingModelParameter(STUZbbbarLRvars[i]);
        }
    }
    return (NPSTU::CheckParameters(DPars));
}

////////////////////////////////////////////////////////////////////////

double NPSTUZbbbarLR::deltaGV_f(const Particle f) const
{
    if (f.is("BOTTOM")) {
        double DeltaGVb = myDeltaGLb + myDeltaGRb;
        return (NPSTU::deltaGV_f(f) + DeltaGVb);
    } else
        return NPSTU::deltaGV_f(f);
}

double NPSTUZbbbarLR::deltaGA_f(const Particle f) const
{
    if (f.is("BOTTOM")) {
        double DeltaGAb = myDeltaGLb - myDeltaGRb;
        return (NPSTU::deltaGA_f(f) + DeltaGAb);
    } else
        return NPSTU::deltaGA_f(f);
}

