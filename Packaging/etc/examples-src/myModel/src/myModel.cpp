/* 
 * Copyright (C) 2015 SusyFit Collaboration
 * All rights reserved.
 *
 * For the licensing terms see doc/COPYING.
 */

#include "myModel.h"

/* Define mandatory model parameters here. */
const std::string myModel::myModelvars[NmyModelvars] = {"ct", "cg", "cV", "cA"};

myModel::myModel()
:   StandardModel()
{
    /* Define all the parameters here and port them as observables too */
    ModelParamMap.insert(std::pair<std::string, boost::reference_wrapper<const double> >("ct", boost::cref(ct)));
    ModelParamMap.insert(std::pair<std::string, boost::reference_wrapper<const double> >("cg", boost::cref(cg)));
    ModelParamMap.insert(std::pair<std::string, boost::reference_wrapper<const double> >("cV", boost::cref(cV)));
    ModelParamMap.insert(std::pair<std::string, boost::reference_wrapper<const double> >("cA", boost::cref(cA)));
}

myModel::~myModel()
{
    if (IsModelInitialized()) {
        /* Destroy whatever you want, e.g. potentially dangerous pointers. */
    }
}

/* Initialize model here */
bool myModel::InitializeModel()
{
    onshell = false;
    setModelInitialized(StandardModel::InitializeModel());
    return(true);
}
    
bool myModel::Init(const std::map<std::string, double>& DPars)
{
    return(StandardModel::Init(DPars));
}

/* Do whatever is necessary before parameters are updated by the MCMC. */
bool myModel::PreUpdate()
{
    if(!StandardModel::PreUpdate()) return (false);
    return (true);
}

/* Model update method used be the MCMC to update the model parameters. */
bool myModel::Update(const std::map<std::string, double>& DPars)
{    
    if(!PreUpdate()) return (false);
    
    UpdateError = false;

    for (std::map<std::string, double>::const_iterator it = DPars.begin(); it != DPars.end(); it++)
        setParameter(it->first, it->second);

    if (UpdateError) return (false);

    if(!PostUpdate()) return (false);

    return (true);

    return (true);
}

/* Postupdate method to update whatever is needed after the model parameters are updated */
bool myModel::PostUpdate()
{
    if(!StandardModel::PostUpdate()) return (false);

    return (true);
}

/* Model parameters and their derived quantities can be set here. */
void myModel::setParameter(const std::string name, const double& value)
{
    if(name.compare("ct") == 0)
        ct = value;
    else if(name.compare("cg") == 0)
        cg = value;
    else if(name.compare("cV") == 0)
        cV = value;
    else if(name.compare("cA") == 0)
        cA = value;
    else
        StandardModel::setParameter(name,value);
}

bool myModel::CheckParameters(const std::map<std::string, double>& DPars)
{
    for (int i = 0; i < NmyModelvars; i++) {
        if (DPars.find(myModelvars[i]) == DPars.end()) {
            std::cout << "missing mandatory myModel parameter " << myModelvars[i] << std::endl;
            return false;
        }
    }
    return(StandardModel::CheckParameters(DPars));
}


/* Model Flags can be set here. */
bool myModel::setFlag(const std::string name, const bool value)
{
    bool res = false;
    
    if(name.compare("onshell") == 0){
        onshell = value;
        res = true;
    } else {
        res = StandardModel::setFlag(name,value);
    }

    return(res);
}
