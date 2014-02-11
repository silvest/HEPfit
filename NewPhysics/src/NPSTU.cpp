/*
 * Copyright (C) 2012-2014 SusyFit Collaboration
 * All rights reserved.
 *
 * For the licensing terms see doc/COPYING.
 */

#include <stdexcept>
#include <EWSM.h>
#include "NPSTU.h"


const std::string NPSTU::STUvars[NSTUvars]
= {"obliqueS", "obliqueT", "obliqueU"};


NPSTU::NPSTU()
: NPbase()
{
}


bool NPSTU::InitializeModel()
{
    setModelInitialized(NPbase::InitializeModel());
    return(IsModelInitialized());
}


bool NPSTU::Init(const std::map<std::string, double>& DPars)
{
    Update(DPars);
    return(CheckParameters(DPars));
}


bool NPSTU::Update(const std::map<std::string,double>& DPars)
{
    for (std::map<std::string, double>::const_iterator it = DPars.begin(); it != DPars.end(); it++)
        setParameter(it->first, it->second);
    if(!NPbase::Update(DPars)) return (false);
    return (true);
}


void NPSTU::setParameter(const std::string name, const double& value)
{
    if (name.compare("obliqueS") == 0)
        myObliqueS = value;
    else if (name.compare("obliqueT") == 0)
        myObliqueT = value;
    else if (name.compare("obliqueU") == 0)
        myObliqueU = value;
    else
        NPbase::setParameter(name, value);
}


bool NPSTU::CheckParameters(const std::map<std::string, double>& DPars)
{
    for (int i = 0; i < NSTUvars; i++) {
        if (DPars.find(STUvars[i]) == DPars.end()) {
            std::cout << "ERROR: Missing mandatory NPSTU parameter "
                      << STUvars[i] << std::endl;
            return false;
        }
    }
    return(NPbase::CheckParameters(DPars));
}


bool NPSTU::setFlag(const std::string name, const bool value)
{
    bool res = false;
    res = NPbase::setFlag(name,value);

    return(res);
}

bool NPSTU::CheckFlags() const
{
    return(NPbase::CheckFlags());
}


////////////////////////////////////////////////////////////////////////

double NPSTU::epsilon1() const
{
    double That = alphaMz()*obliqueT();

    return ( myEWSM->epsilon1_SM() + That );
}


double NPSTU::epsilon2() const
{
    double sW2_SM = myEWSM->sW2_SM(); /* This has to be the SM value. */
    double Uhat = - alphaMz()/(4.0*sW2_SM)*obliqueU();

    return ( myEWSM->epsilon2_SM() + Uhat );
}


double NPSTU::epsilon3() const
{
    double sW2_SM = myEWSM->sW2_SM(); /* This has to be the SM value. */
    double Shat = alphaMz()/(4.0*sW2_SM)*obliqueS();

    return ( myEWSM->epsilon3_SM() + Shat );
}


double NPSTU::epsilonb() const
{
    return myEWSM->epsilonb_SM();
}


