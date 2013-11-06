/*
 * Copyright (C) 2013 SusyFit Collaboration
 * All rights reserved.
 *
 * For the licensing terms see doc/COPYING.
 */

#include <stdexcept>
#include <EWSM.h>
#include "NPEpsilons_pureNP.h"

const std::string NPEpsilons_pureNP::EPSILONpureNPvars[NEPSILONpureNPvars]
= {"deltaEpsilon_1", "deltaEpsilon_2", "deltaEpsilon_3", "deltaEpsilon_b"};

//const std::string NPEpsilons_pureNP::EPSILONpureNPflags[NEPSILONpureNPflags]
//= {};


NPEpsilons_pureNP::NPEpsilons_pureNP()
: NPbase()
{
}


bool NPEpsilons_pureNP::Update(const std::map<std::string,double>& DPars)
{
    for (std::map<std::string, double>::const_iterator it = DPars.begin(); it != DPars.end(); it++)
        setParameter(it->first, it->second);
    if(!NPbase::Update(DPars)) return (false);
    return (true);
}


bool NPEpsilons_pureNP::Init(const std::map<std::string, double>& DPars)
{
    Update(DPars);
    return(CheckParameters(DPars));
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
    return(NPbase::CheckParameters(DPars));
}


void NPEpsilons_pureNP::setParameter(const std::string name, const double& value)
{
    if (name.compare("deltaEpsilon_1") == 0)
        deltaEps_1 = value;
    else if (name.compare("deltaEpsilon_2") == 0)
        deltaEps_2 = value;
    else if (name.compare("deltaEpsilon_3") == 0)
        deltaEps_3 = value;
    else if (name.compare("deltaEpsilon_b") == 0)
        deltaEps_b = value;
    else
        NPbase::setParameter(name, value);
}


bool NPEpsilons_pureNP::InitializeModel()
{
    setModelInitialized(NPbase::InitializeModel());
    return (IsModelInitialized());
}


void NPEpsilons_pureNP::setEWSMflags(EWSM& myEWSM)
{
    NPbase::setEWSMflags(myEWSM);
}


bool NPEpsilons_pureNP::setFlag(const std::string name, const bool& value)
{
    bool res = false;
    res = NPbase::setFlag(name,value);
    return(res);
}


bool NPEpsilons_pureNP::CheckFlags() const
{
    return(NPbase::CheckFlags());
}


////////////////////////////////////////////////////////////////////////

double NPEpsilons_pureNP::epsilon1() const
{
    return (epsilon1_SM() + deltaEps_1);
}


double NPEpsilons_pureNP::epsilon2() const
{
    return (epsilon2_SM() + deltaEps_2);
}


double NPEpsilons_pureNP::epsilon3() const
{
    return (epsilon3_SM() + deltaEps_3);
}


double NPEpsilons_pureNP::epsilonb() const
{
    return (epsilonb_SM() + deltaEps_b);
}


////////////////////////////////////////////////////////////////////////

double NPEpsilons_pureNP::Mw() const
{
    double Mw_SM = StandardModel::Mw();
    double sW2_SM = StandardModel::sW2();
    double cW2_SM = StandardModel::cW2();
    return ( Mw_SM*(1.0 - (- c02()*deltaEps_1
                           + (c02() - s02())*deltaEps_2
                           + 2.0*s02()*deltaEps_3)/(cW2_SM - sW2_SM)/2.0) );
}


double NPEpsilons_pureNP::cW2() const
{
    return ( Mw()*Mw()/Mz/Mz );
}


double NPEpsilons_pureNP::sW2() const
{
    return ( 1.0 - cW2() );
}


double NPEpsilons_pureNP::deltaGVl(StandardModel::lepton l) const
{
    /* SM values */
    double gV_SM = StandardModel::gVl(l).real();
    double gA_SM = StandardModel::gAl(l).real();

    return ( gV_SM*deltaEps_1/2.0
            + (gV_SM - gA_SM)*(deltaEps_3 - c02()*deltaEps_1)/(c02() - s02()) );
}


double NPEpsilons_pureNP::deltaGVq(StandardModel::quark q) const
{
    if (q==TOP) return 0.0;

    /* SM values */
    double gV_SM = StandardModel::gVq(q).real();
    double gA_SM = StandardModel::gAq(q).real();

    if (q==BOTTOM)
        return ( gV_SM*deltaEps_1/2.0 + gV_SM*deltaEps_b
                + (gV_SM - gA_SM)*(deltaEps_3 - c02()*deltaEps_1)/(c02() - s02())
                - (gV_SM - gA_SM)*deltaEps_b );
    else
        return ( gV_SM*deltaEps_1/2.0
                + (gV_SM - gA_SM)*(deltaEps_3 - c02()*deltaEps_1)/(c02() - s02()) );
}


double NPEpsilons_pureNP::deltaGAl(StandardModel::lepton l) const
{
    return ( getLeptons(l).getIsospin()*deltaEps_1/2.0 );
}


double NPEpsilons_pureNP::deltaGAq(StandardModel::quark q) const
{
    if (q==TOP)
        return 0.0;
    else if (q==BOTTOM)
        return ( getQuarks(q).getIsospin()*(deltaEps_1/2.0 + deltaEps_b) );
    else
        return ( getQuarks(q).getIsospin()*deltaEps_1/2.0 );
}


double NPEpsilons_pureNP::GammaW() const
{
    throw std::runtime_error("NPEpsilons_pureNP::GammaW() is not implemented.");
}



