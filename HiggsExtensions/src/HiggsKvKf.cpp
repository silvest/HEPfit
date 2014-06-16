/* 
 * Copyright (C) 2014 SusyFit Collaboration
 * All rights reserved.
 *
 * For the licensing terms see doc/COPYING.
 */

#include "HiggsKvKf.h"

const std::string HiggsKvKf::HKvKfvars[NHKvKfvars] = {
    "Kv", "Kf"
};

HiggsKvKf::HiggsKvKf() : NPbase()
{
    ModelParamMap.insert(std::pair<std::string, boost::reference_wrapper<const double> >("Kv", boost::cref(Kv)));
    ModelParamMap.insert(std::pair<std::string, boost::reference_wrapper<const double> >("Kf", boost::cref(Kf)));
}

void HiggsKvKf::setParameter(const std::string name, const double& value)
{
    if (name.compare("Kv") == 0)
        Kv = value;
    else if (name.compare("Kf") == 0)
        Kf = value;
    else
        StandardModel::setParameter(name, value);
}

bool HiggsKvKf::CheckParameters(const std::map<std::string, double>& DPars)
{
    for (int i = 0; i < NHKvKfvars; i++) {
        if (DPars.find(HKvKfvars[i]) == DPars.end()) {
            std::cout << "missing mandatory HiggsKvKf parameter " << HKvKfvars[i] << std::endl;
            return false;
        }
    }
    return (StandardModel::CheckParameters(DPars));
}

double HiggsKvKf::computeKW() const
{
    return Kv;
}

double HiggsKvKf::computeKZ() const
{
    return Kv;
}

double HiggsKvKf::computeKZga() const
{
    double gtt = computeGammaZgatt();
    double gWW = computeGammaZgaWW();
    double gtW = computeGammaZgatW();
    return sqrt((gtt * Kf * Kf + gWW * Kv * Kv + gtW * Kf * Kv) / (gtt + gWW + gtW));
}

double HiggsKvKf::computeKgaga() const
{
    double gtt = computeGammagagatt();
    double gWW = computeGammagagaWW();
    double gtW = computeGammagagatW();
    return sqrt((gtt * Kf * Kf + gWW * Kv * Kv + gtW * Kf * Kv) / (gtt + gWW + gtW));
}

double HiggsKvKf::computeKb() const
{
    return Kf;
}

double HiggsKvKf::computeKc() const
{
    return Kf;
}

double HiggsKvKf::computeKglgl() const
{
    return Kf;
}

double HiggsKvKf::computeKt() const
{
    return Kf;
}

double HiggsKvKf::computeKtau() const
{
    return Kf;
}

double HiggsKvKf::computeGTotalRatio() const
{
    return computeKW() * computeKW() * computeBRWW() +
            computeKZ() * computeKZ() * computeBRZZ() +
            computeKZga() * computeKZga() * computeBRZga() +
            computeKgaga() * computeKgaga() * computeBRgaga() +
            computeKglgl() * computeKglgl() * computeBRglgl() +
            computeKb() * computeKb() * computeBRbb() +
            computeKc() * computeKc() * computeBRcc() +
            computeKtau() * computeKtau() * computeBRtautau();
}
