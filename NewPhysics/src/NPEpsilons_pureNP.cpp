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
= {"delEps_1", "delEps_2", "delEps_3", "delEps_b"};


NPEpsilons_pureNP::NPEpsilons_pureNP()
: NPbase()
{
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


bool NPEpsilons_pureNP::Init(const std::map<std::string, double>& DPars)
{
    Update(DPars);
    return(CheckParameters(DPars));
}


bool NPEpsilons_pureNP::Update(const std::map<std::string,double>& DPars)
{
    for (std::map<std::string, double>::const_iterator it = DPars.begin(); it != DPars.end(); it++)
        setParameter(it->first, it->second);
    if(!NPbase::Update(DPars)) return (false);
    return (true);
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
    return(NPbase::CheckParameters(DPars));
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
    return (myEWSM->epsilon1_SM() + deltaEps_1);
}


double NPEpsilons_pureNP::epsilon2() const
{
    return (myEWSM->epsilon2_SM() + deltaEps_2);
}


double NPEpsilons_pureNP::epsilon3() const
{
    return (myEWSM->epsilon3_SM() + deltaEps_3);
}


double NPEpsilons_pureNP::epsilonb() const
{
    return (myEWSM->epsilonb_SM() + deltaEps_b);
}


////////////////////////////////////////////////////////////////////////

double NPEpsilons_pureNP::Mw() const
{
    double Mw_SM = myEWSM->Mw_SM();
    double sW2_SM = myEWSM->sW2_SM();
    double cW2_SM = myEWSM->cW2_SM();
    return ( Mw_SM*(1.0 - (- myEWSM->c02()*deltaEps_1
                           + (myEWSM->c02() - myEWSM->s02())*deltaEps_2
                           + 2.0*myEWSM->s02()*deltaEps_3)/(cW2_SM - sW2_SM)/2.0) );
}


double NPEpsilons_pureNP::cW2() const
{
    return ( Mw()*Mw()/Mz/Mz );
}


double NPEpsilons_pureNP::sW2() const
{
    return ( 1.0 - cW2() );
}


double NPEpsilons_pureNP::GammaW() const
{
    throw std::runtime_error("NPEpsilons_pureNP::GammaW() is not implemented.");
}


////////////////////////////////////////////////////////////////////////

double NPEpsilons_pureNP::deltaGVl(StandardModel::lepton l) const
{
    /* SM values */
    double gV_SM = myEWSM->gVl_SM(l).real();
    double gA_SM = myEWSM->gAl_SM(l).real();

    return ( gV_SM*deltaEps_1/2.0
            + (gV_SM - gA_SM)*(deltaEps_3 - myEWSM->c02()*deltaEps_1)
              /(myEWSM->c02() - myEWSM->s02()) );
}


double NPEpsilons_pureNP::deltaGVq(StandardModel::quark q) const
{
    if (q==TOP) return 0.0;

    /* SM values */
    double gV_SM = myEWSM->gVq_SM(q).real();
    double gA_SM = myEWSM->gAq_SM(q).real();

    if (q==BOTTOM)
        return ( gV_SM*deltaEps_1/2.0 + gV_SM*deltaEps_b
                + (gV_SM - gA_SM)*(deltaEps_3 - myEWSM->c02()*deltaEps_1)
                  /(myEWSM->c02() - myEWSM->s02())
                - (gV_SM - gA_SM)*deltaEps_b );
    else
        return ( gV_SM*deltaEps_1/2.0
                + (gV_SM - gA_SM)*(deltaEps_3 - myEWSM->c02()*deltaEps_1)
                  /(myEWSM->c02() - myEWSM->s02()) );
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


