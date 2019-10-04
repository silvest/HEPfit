/*
 * Copyright (C) 2013 HEPfit Collaboration
 *
 *
 * For the licensing terms see doc/COPYING.
 */

#include "NPZbbbar.h"
#include <stdexcept>


const std::string NPZbbbar::ZbbbarVAVars[NZbbbarVars]
        = {"deltaGVb", "deltaGAb"};

const std::string NPZbbbar::ZbbbarLRVars[NZbbbarVars]
        = {"deltaGLb", "deltaGRb"};

NPZbbbar::NPZbbbar(const bool FlagNPZbbbarLR_in)
: NPbase(), FlagNPZbbbarLR(FlagNPZbbbarLR_in)
{
    setFlagNoApproximateGammaZ(true);

    ModelParamMap.insert(std::make_pair("deltaGVb", std::cref(myDeltaGVb)));
    ModelParamMap.insert(std::make_pair("deltaGAb", std::cref(myDeltaGAb)));
    ModelParamMap.insert(std::make_pair("deltaGLb", std::cref(myDeltaGLb)));
    ModelParamMap.insert(std::make_pair("deltaGRb", std::cref(myDeltaGRb)));
}

bool NPZbbbar::PostUpdate()
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

void NPZbbbar::setParameter(const std::string name, const double& value)
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

bool NPZbbbar::CheckParameters(const std::map<std::string, double>& DPars)
{
    for (int i = 0; i < NZbbbarVars; i++) {
        if (FlagNPZbbbarLR) {
            if (DPars.find(ZbbbarLRVars[i]) == DPars.end()) {
                std::cout << "ERROR: Missing mandatory NPZbbbarLR parameter "
                        << ZbbbarLRVars[i] << std::endl;
                raiseMissingModelParameterCount();
                addMissingModelParameter(ZbbbarLRVars[i]);
            }
        } else {
            if (DPars.find(ZbbbarVAVars[i]) == DPars.end()) {
                std::cout << "ERROR: Missing mandatory NPZbbbar parameter "
                        << ZbbbarVAVars[i] << std::endl;
                raiseMissingModelParameterCount();
                addMissingModelParameter(ZbbbarVAVars[i]);
            }
        }
    }
    return (NPbase::CheckParameters(DPars));
}

////////////////////////////////////////////////////////////////////////

double NPZbbbar::deltaGV_f(const Particle f) const
{
    if (f.is("BOTTOM"))
        return myDeltaGVb;
    else
        return 0.0;
}

double NPZbbbar::deltaGA_f(const Particle f) const
{
    if (f.is("BOTTOM"))
        return myDeltaGAb;
    else
        return 0.0;
}

////////////////////////////////////////////////////////////////////////

double NPZbbbar::Mw() const
{
    return (StandardModel::Mw());
}

double NPZbbbar::GammaW() const
{
    return (StandardModel::GammaW());
}

double NPZbbbar::Gamma_Z() const
{
    return (StandardModel::Gamma_Z());
}

double NPZbbbar::sigma0_had() const
{
    return (StandardModel::sigma0_had());
}

double NPZbbbar::sin2thetaEff(const Particle f) const
{
    return (StandardModel::sin2thetaEff(f));
}

double NPZbbbar::A_f(const Particle f) const
{
    return (StandardModel::A_f(f));
}

double NPZbbbar::AFB(const Particle f) const
{
    return (StandardModel::AFB(f));
}

double NPZbbbar::R0_f(const Particle f) const
{
    return (StandardModel::R0_f(f));
}

