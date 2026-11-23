/* 
 * Copyright (C) 2015 HEPfit Collaboration
 *
 *
 * For the licensing terms see doc/COPYING.
 */

#include "FlavourWilsonCoefficient.h"

const std::string FlavourWilsonCoefficient::FlavourWilsonCoefficientvars[NFlavourWilsonCoefficientvars] = {"reDC1","reDC2","reDC3","reDC4","reDC5","reDC6","reDC7","reDC8","reDC9","reDC10","reDC7g","reDC8g",
                                                                                                           "imDC1","imDC2","imDC3","imDC4","imDC5","imDC6","imDC7","imDC8","imDC9","imDC10","imDC7g","imDC8g",
                                                                                                           "reDC7p","reDC8p","reDC9p","reDC10p","reDC7gp","reDC8gp",
                                                                                                           "imDC7p","imDC8p","imDC9p","imDC10p","imDC7gp","imDC8gp",
                                                                                                           "WCscale"};

FlavourWilsonCoefficient::FlavourWilsonCoefficient() : StandardModel(), FWCM(*this) {   

    SMM.setObj((StandardModelMatching&) FWCM.getObj());
    ModelParamMap.insert(std::make_pair("reDC1", std::cref(reDC1)));
    ModelParamMap.insert(std::make_pair("reDC2", std::cref(reDC2)));
    ModelParamMap.insert(std::make_pair("reDC3", std::cref(reDC3)));
    ModelParamMap.insert(std::make_pair("reDC4", std::cref(reDC4)));
    ModelParamMap.insert(std::make_pair("reDC5", std::cref(reDC5)));
    ModelParamMap.insert(std::make_pair("reDC6", std::cref(reDC6)));
    ModelParamMap.insert(std::make_pair("reDC7", std::cref(reDC7)));
    ModelParamMap.insert(std::make_pair("reDC8", std::cref(reDC8)));
    ModelParamMap.insert(std::make_pair("reDC9", std::cref(reDC9)));
    ModelParamMap.insert(std::make_pair("reDC10", std::cref(reDC10)));
    ModelParamMap.insert(std::make_pair("reDC7g", std::cref(reDC7g)));
    ModelParamMap.insert(std::make_pair("reDC8g", std::cref(reDC8g)));
    
    ModelParamMap.insert(std::make_pair("imDC1", std::cref(imDC1)));
    ModelParamMap.insert(std::make_pair("imDC2", std::cref(imDC2)));
    ModelParamMap.insert(std::make_pair("imDC3", std::cref(imDC3)));
    ModelParamMap.insert(std::make_pair("imDC4", std::cref(imDC4)));
    ModelParamMap.insert(std::make_pair("imDC5", std::cref(imDC5)));
    ModelParamMap.insert(std::make_pair("imDC6", std::cref(imDC6)));
    ModelParamMap.insert(std::make_pair("imDC7", std::cref(imDC7)));
    ModelParamMap.insert(std::make_pair("imDC8", std::cref(imDC8)));
    ModelParamMap.insert(std::make_pair("imDC9", std::cref(imDC9)));
    ModelParamMap.insert(std::make_pair("imDC10", std::cref(imDC10)));
    ModelParamMap.insert(std::make_pair("imDC7g", std::cref(imDC7g)));
    ModelParamMap.insert(std::make_pair("imDC8g", std::cref(imDC8g)));
    
    ModelParamMap.insert(std::make_pair("reDC7p", std::cref(reDC7p)));
    ModelParamMap.insert(std::make_pair("reDC8p", std::cref(reDC8p)));
    ModelParamMap.insert(std::make_pair("reDC9p", std::cref(reDC9p)));
    ModelParamMap.insert(std::make_pair("reDC10p", std::cref(reDC10p)));
    ModelParamMap.insert(std::make_pair("reDC7gp", std::cref(reDC7gp)));
    ModelParamMap.insert(std::make_pair("reDC8gp", std::cref(reDC8gp)));
    
    ModelParamMap.insert(std::make_pair("imDC7p", std::cref(imDC7p)));
    ModelParamMap.insert(std::make_pair("imDC8p", std::cref(imDC8p)));
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
    
    DC1 = gslpp::complex(reDC1, imDC1, false);
    DC2 = gslpp::complex(reDC2, imDC2, false);
    DC3 = gslpp::complex(reDC3, imDC3, false);
    DC4 = gslpp::complex(reDC4, imDC4, false);
    DC5 = gslpp::complex(reDC5, imDC5, false);
    DC6 = gslpp::complex(reDC6, imDC6, false);
    DC7 = gslpp::complex(reDC7, imDC7, false);
    DC8 = gslpp::complex(reDC8, imDC8, false);
    DC9 = gslpp::complex(reDC9, imDC9, false);
    DC10 = gslpp::complex(reDC10, imDC10, false);
    DC7g = gslpp::complex(reDC7g, imDC7g, false);
    DC8g = gslpp::complex(reDC8g, imDC8g, false);
    
    DC7p = gslpp::complex(reDC7p, imDC7p, false);
    DC8p = gslpp::complex(reDC8p, imDC8p, false);
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
    if(name.compare("reDC1") == 0) 
        reDC1 = value;
    else if(name.compare("reDC2") == 0) 
        reDC2 = value;
    else if(name.compare("reDC3") == 0)
        reDC3 = value;
    else if(name.compare("reDC4") == 0)
        reDC4 = value;
    else if(name.compare("reDC5") == 0)
        reDC5 = value;
    else if(name.compare("reDC6") == 0)
        reDC6 = value;
    else if(name.compare("reDC7") == 0)
        reDC7 = value;
    else if(name.compare("reDC8") == 0)
        reDC8 = value;
    else if(name.compare("reDC9") == 0)
        reDC9 = value;
    else if(name.compare("reDC10") == 0)
        reDC10 = value;
    else if(name.compare("reDC7g") == 0)
        reDC7g = value;
    else if(name.compare("reDC8g") == 0)
        reDC8g = value;
    else if(name.compare("imDC1") == 0) 
        imDC1 = value;
    else if(name.compare("imDC2") == 0) 
        imDC2 = value;
    else if(name.compare("imDC3") == 0)
        imDC3 = value;
    else if(name.compare("imDC4") == 0)
        imDC4 = value;
    else if(name.compare("imDC5") == 0)
        imDC5 = value;
    else if(name.compare("imDC6") == 0)
        imDC6 = value;
    else if(name.compare("imDC7") == 0)
        imDC7 = value;
    else if(name.compare("imDC8") == 0)
        imDC8 = value;
    else if(name.compare("imDC9") == 0)
        imDC9 = value;
    else if(name.compare("imDC10") == 0)
        imDC10 = value;
    else if(name.compare("imDC7g") == 0)
        imDC7g = value;
    else if(name.compare("imDC8g") == 0)
        imDC8g = value;
    else if(name.compare("reDC7p") == 0)
        reDC7p = value;
    else if(name.compare("reDC8p") == 0)
        reDC8p = value;
    else if(name.compare("reDC9p") == 0)
        reDC9p = value;
    else if(name.compare("reDC10p") == 0)
        reDC10p = value;
    else if(name.compare("reDC7gp") == 0)
        reDC7gp = value;
    else if(name.compare("reDC8gp") == 0)
        reDC8gp = value;
    else if(name.compare("imDC7p") == 0)
        imDC7p = value;
    else if(name.compare("imDC8p") == 0)
        imDC8p = value;
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
