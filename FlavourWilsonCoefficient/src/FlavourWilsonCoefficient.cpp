/* 
 * Copyright (C) 2015 HEPfit Collaboration
 * All rights reserved.
 *
 * For the licensing terms see doc/COPYING.
 */

#include "FlavourWilsonCoefficient.h"

const std::string FlavourWilsonCoefficient::FlavourWilsonCoefficientvars[NFlavourWilsonCoefficientvars] = {"DC1","DC2","DC3","DC4","DC5","DC6","DC7","DC8","DC9","DC10","DC7g","DC8g"};

FlavourWilsonCoefficient::FlavourWilsonCoefficient() : StandardModel() {   

    
    ModelParamMap.insert(std::pair<std::string, boost::reference_wrapper<const double> >("DC1", boost::cref(DC1)));
    ModelParamMap.insert(std::pair<std::string, boost::reference_wrapper<const double> >("DC2", boost::cref(DC2)));
    ModelParamMap.insert(std::pair<std::string, boost::reference_wrapper<const double> >("DC3", boost::cref(DC3)));
    ModelParamMap.insert(std::pair<std::string, boost::reference_wrapper<const double> >("DC4", boost::cref(DC4)));
    ModelParamMap.insert(std::pair<std::string, boost::reference_wrapper<const double> >("DC5", boost::cref(DC5)));
    ModelParamMap.insert(std::pair<std::string, boost::reference_wrapper<const double> >("DC6", boost::cref(DC6)));
    ModelParamMap.insert(std::pair<std::string, boost::reference_wrapper<const double> >("DC7", boost::cref(DC7)));
    ModelParamMap.insert(std::pair<std::string, boost::reference_wrapper<const double> >("DC8", boost::cref(DC8)));
    ModelParamMap.insert(std::pair<std::string, boost::reference_wrapper<const double> >("DC9", boost::cref(DC9)));
    ModelParamMap.insert(std::pair<std::string, boost::reference_wrapper<const double> >("DC10", boost::cref(DC10)));
    ModelParamMap.insert(std::pair<std::string, boost::reference_wrapper<const double> >("DC7g", boost::cref(DC7g)));
    ModelParamMap.insert(std::pair<std::string, boost::reference_wrapper<const double> >("DC8g", boost::cref(DC8g)));
}

FlavourWilsonCoefficient::~FlavourWilsonCoefficient(){
    if (IsModelInitialized()) {
            if (myFlavourWilsonCoefficientMatching != NULL) delete(myFlavourWilsonCoefficientMatching);
        }
}

///////////////////////////////////////////////////////////////////////////
// Initialization

bool FlavourWilsonCoefficient::InitializeModel()
{
    myFlavourWilsonCoefficientMatching = new FlavourWilsonCoefficientMatching(*this);
    setModelInitialized(StandardModel::InitializeModel());
    return(true);
}
    
bool FlavourWilsonCoefficient::Init(const std::map<std::string, double>& DPars) {
    return(StandardModel::Init(DPars));
}

bool FlavourWilsonCoefficient::PreUpdate()
{    
    if(!StandardModel::PreUpdate()) return (false);

    return (true);
}

bool FlavourWilsonCoefficient::Update(const std::map<std::string, double>& DPars) {
    
    if(!PreUpdate()) return (false);
    
    UpdateError = false;

    for (std::map<std::string, double>::const_iterator it = DPars.begin(); it != DPars.end(); it++)
        setParameter(it->first, it->second);

    if (UpdateError) return (false);

    if(!PostUpdate()) return (false);

    return (true);
}

bool FlavourWilsonCoefficient::PostUpdate()
{
    if(!StandardModel::PostUpdate()) return (false);

    /* Necessary for updating StandardModel parameters in StandardModelMatching,
     * and THDM and THDM-derived parameters in THDMMatching */
    myFlavourWilsonCoefficientMatching->updateFlavourWilsonCoefficientParameters();

    return (true);
}

void FlavourWilsonCoefficient::setParameter(const std::string name, const double& value){    
    if(name.compare("DC1") == 0) 
        DC1 = value;
    else if(name.compare("DC2") == 0) 
        DC2 = value;
    else if(name.compare("DC3") == 0)
        DC3 = value;
    else if(name.compare("DC4") == 0)
        DC4 = value;
    else if(name.compare("DC5") == 0)
        DC5 = value;
    else if(name.compare("DC6") == 0)
        DC6 = value;
    else if(name.compare("DC7") == 0)
        DC7 = value;
    else if(name.compare("DC8") == 0)
        DC8 = value;
    else if(name.compare("DC9") == 0)
        DC9 = value;
    else if(name.compare("DC10") == 0)
        DC10 = value;
    else if(name.compare("DC7g") == 0)
        DC7g = value;
    else if(name.compare("DC8g") == 0)
        DC8g = value;
    else
        StandardModel::setParameter(name,value);
}

bool FlavourWilsonCoefficient::CheckParameters(const std::map<std::string, double>& DPars) {
    for (int i = 0; i < NFlavourWilsonCoefficientvars; i++) {
        if (DPars.find(FlavourWilsonCoefficientvars[i]) == DPars.end()) {
            std::cout << "missing mandatory THDM parameter " << FlavourWilsonCoefficientvars[i] << std::endl;
            return false;
        }
    }
    return(StandardModel::CheckParameters(DPars));
}

///////////////////////////////////////////////////////////////////////////
// Flags

bool FlavourWilsonCoefficient::setFlag(const std::string name, const bool value)
{
    bool res = false;
    
    res = StandardModel::setFlag(name,value);

    return(res);
}
