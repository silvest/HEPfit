/* 
 * Copyright (C) 2015 HEPfit Collaboration
 *
 *
 * For the licensing terms see doc/COPYING.
 */

#include "FlavourWilsonCoefficient.h"

const std::string FlavourWilsonCoefficient::FlavourWilsonCoefficientvars[NFlavourWilsonCoefficientvars] = {"reDC9","reDC10","reDC7g","reDC8g",
                                                                                                           "imDC9","imDC10","imDC7g","imDC8g",
                                                                                                           "reDC9p","reDC10p","reDC7gp","reDC8gp",
                                                                                                           "imDC9p","imDC10p","imDC7gp","imDC8gp",
                                                                                                           "WCscale"};

FlavourWilsonCoefficient::FlavourWilsonCoefficient() : StandardModel(), FWCM(*this) {   

    SMM.setObj((StandardModelMatching&) FWCM.getObj());
    ModelParamMap.insert(std::make_pair("reDC9", std::cref(reDC9)));
    ModelParamMap.insert(std::make_pair("reDC10", std::cref(reDC10)));
    ModelParamMap.insert(std::make_pair("reDC7g", std::cref(reDC7g)));
    ModelParamMap.insert(std::make_pair("reDC8g", std::cref(reDC8g)));
    
    ModelParamMap.insert(std::make_pair("imDC9", std::cref(imDC9)));
    ModelParamMap.insert(std::make_pair("imDC10", std::cref(imDC10)));
    ModelParamMap.insert(std::make_pair("imDC7g", std::cref(imDC7g)));
    ModelParamMap.insert(std::make_pair("imDC8g", std::cref(imDC8g)));
    
    ModelParamMap.insert(std::make_pair("reDC9p", std::cref(reDC9p)));
    ModelParamMap.insert(std::make_pair("reDC10p", std::cref(reDC10p)));
    ModelParamMap.insert(std::make_pair("reDC7gp", std::cref(reDC7gp)));
    ModelParamMap.insert(std::make_pair("reDC8gp", std::cref(reDC8gp)));
    
    ModelParamMap.insert(std::make_pair("imDC9p", std::cref(imDC9p)));
    ModelParamMap.insert(std::make_pair("imDC10p", std::cref(imDC10p)));
    ModelParamMap.insert(std::make_pair("imDC7gp", std::cref(imDC7gp)));
    ModelParamMap.insert(std::make_pair("imDC8gp", std::cref(imDC8gp)));
    
    ModelParamMap.insert(std::make_pair("WCscale", std::cref(WCscale)));
}

FlavourWilsonCoefficient::~FlavourWilsonCoefficient(){
    if (IsModelInitialized()) {
        }
}

///////////////////////////////////////////////////////////////////////////
// Initialization

bool FlavourWilsonCoefficient::InitializeModel()
{
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
    
    DC9 = gslpp::complex(reDC9, imDC9, false);
    DC10 = gslpp::complex(reDC10, imDC10, false);
    DC7g = gslpp::complex(reDC7g, imDC7g, false);
    DC8g = gslpp::complex(reDC8g, imDC8g, false);
    
    DC9p = gslpp::complex(reDC9p, imDC9p, false);
    DC10p = gslpp::complex(reDC10p, imDC10p, false);
    DC7gp = gslpp::complex(reDC7gp, imDC7gp, false);
    DC8gp = gslpp::complex(reDC8gp, imDC8gp, false);

    /* Necessary for updating StandardModel parameters in StandardModelMatching,
     * and FlavourWC and FlavourWC-derived parameters in FlavourWCMatching */
    FWCM.getObj().updateFlavourWilsonCoefficientParameters();

    return (true);
}

void FlavourWilsonCoefficient::setParameter(const std::string name, const double& value){    
    if(name.compare("reDC9") == 0)
        reDC9 = value;
    else if(name.compare("reDC10") == 0)
        reDC10 = value;
    else if(name.compare("reDC7g") == 0)
        reDC7g = value;
    else if(name.compare("reDC8g") == 0)
        reDC8g = value;
    else if(name.compare("imDC9") == 0)
        imDC9 = value;
    else if(name.compare("imDC10") == 0)
        imDC10 = value;
    else if(name.compare("imDC7g") == 0)
        imDC7g = value;
    else if(name.compare("imDC8g") == 0)
        imDC8g = value;
    else if(name.compare("reDC9p") == 0)
        reDC9p = value;
    else if(name.compare("reDC10p") == 0)
        reDC10p = value;
    else if(name.compare("reDC7gp") == 0)
        reDC7gp = value;
    else if(name.compare("reDC8gp") == 0)
        reDC8gp = value;
    else if(name.compare("imDC9p") == 0)
        imDC9p = value;
    else if(name.compare("imDC10p") == 0)
        imDC10p = value;
    else if(name.compare("imDC7gp") == 0)
        imDC7gp = value;
    else if(name.compare("imDC8gp") == 0)
        imDC8gp = value;
    else if(name.compare("WCscale") == 0)
        WCscale = value;
    else
        StandardModel::setParameter(name,value);
}

bool FlavourWilsonCoefficient::CheckParameters(const std::map<std::string, double>& DPars) {
    for (int i = 0; i < NFlavourWilsonCoefficientvars; i++) {
        if (DPars.find(FlavourWilsonCoefficientvars[i]) == DPars.end()) {
            std::cout << "ERROR: missing mandatory FlavourWilsonCoefficient parameter " << FlavourWilsonCoefficientvars[i] << std::endl;
            raiseMissingModelParameterCount();
            addMissingModelParameter(FlavourWilsonCoefficientvars[i]);
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
