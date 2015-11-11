/* 
 * Copyright (C) 2015 HEPfit Collaboration
 * All rights reserved.
 *
 * For the licensing terms see doc/COPYING.
 */

#include "myModel.h"

/* Define mandatory model parameters here. */
const std::string myModel::myModelvars[NmyModelvars] = {"c1", "c2", "c3", "c4"};

myModel::myModel()
:   StandardModel()
{
    /* Define all the parameters here and port them as observables too */
    ModelParamMap.insert(std::pair<std::string, boost::reference_wrapper<const double> >("c1", boost::cref(c1)));
    ModelParamMap.insert(std::pair<std::string, boost::reference_wrapper<const double> >("c2", boost::cref(c2)));
    ModelParamMap.insert(std::pair<std::string, boost::reference_wrapper<const double> >("c3", boost::cref(c3)));
    ModelParamMap.insert(std::pair<std::string, boost::reference_wrapper<const double> >("c4", boost::cref(c4)));
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
    condition = false;
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
    if(name.compare("c1") == 0)
        c1 = value;
    else if(name.compare("c2") == 0)
        c2 = value;
    else if(name.compare("c3") == 0)
        c3 = value;
    else if(name.compare("c4") == 0)
        c4 = value;
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
    
    if(name.compare("condition") == 0){
        condition = value;
        res = true;
    } else {
        res = StandardModel::setFlag(name,value);
    }

    return(res);
}
