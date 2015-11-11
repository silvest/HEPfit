/* 
 * Copyright (C) 2015 HEPfit Collaboration
 * All rights reserved.
 *
 * For the licensing terms see doc/COPYING.
 */

#include "FlavourWilsonCoefficient.h"

const std::string FlavourWilsonCoefficient::FlavourWilsonCoefficientvars[NFlavourWilsonCoefficientvars] = {"absDC1","absDC2","absDC3","absDC4","absDC5","absDC6","absDC7","absDC8","absDC9","absDC10","absDC7g","absDC8g",
                                                                                                           "argDC1","argDC2","argDC3","argDC4","argDC5","argDC6","argDC7","argDC8","argDC9","argDC10","argDC7g","argDC8g",
                                                                                                           "absDC7p","absDC8p","absDC9p","absDC10p","absDC7gp","absDC8gp",
                                                                                                           "argDC7p","argDC8p","argDC9p","argDC10p","argDC7gp","argDC8gp"};

FlavourWilsonCoefficient::FlavourWilsonCoefficient() : StandardModel() {   

    
    ModelParamMap.insert(std::pair<std::string, boost::reference_wrapper<const double> >("absDC1", boost::cref(absDC1)));
    ModelParamMap.insert(std::pair<std::string, boost::reference_wrapper<const double> >("absDC2", boost::cref(absDC2)));
    ModelParamMap.insert(std::pair<std::string, boost::reference_wrapper<const double> >("absDC3", boost::cref(absDC3)));
    ModelParamMap.insert(std::pair<std::string, boost::reference_wrapper<const double> >("absDC4", boost::cref(absDC4)));
    ModelParamMap.insert(std::pair<std::string, boost::reference_wrapper<const double> >("absDC5", boost::cref(absDC5)));
    ModelParamMap.insert(std::pair<std::string, boost::reference_wrapper<const double> >("absDC6", boost::cref(absDC6)));
    ModelParamMap.insert(std::pair<std::string, boost::reference_wrapper<const double> >("absDC7", boost::cref(absDC7)));
    ModelParamMap.insert(std::pair<std::string, boost::reference_wrapper<const double> >("absDC8", boost::cref(absDC8)));
    ModelParamMap.insert(std::pair<std::string, boost::reference_wrapper<const double> >("absDC9", boost::cref(absDC9)));
    ModelParamMap.insert(std::pair<std::string, boost::reference_wrapper<const double> >("absDC10", boost::cref(absDC10)));
    ModelParamMap.insert(std::pair<std::string, boost::reference_wrapper<const double> >("absDC7g", boost::cref(absDC7g)));
    ModelParamMap.insert(std::pair<std::string, boost::reference_wrapper<const double> >("absDC8g", boost::cref(absDC8g)));
    
    ModelParamMap.insert(std::pair<std::string, boost::reference_wrapper<const double> >("argDC1", boost::cref(argDC1)));
    ModelParamMap.insert(std::pair<std::string, boost::reference_wrapper<const double> >("argDC2", boost::cref(argDC2)));
    ModelParamMap.insert(std::pair<std::string, boost::reference_wrapper<const double> >("argDC3", boost::cref(argDC3)));
    ModelParamMap.insert(std::pair<std::string, boost::reference_wrapper<const double> >("argDC4", boost::cref(argDC4)));
    ModelParamMap.insert(std::pair<std::string, boost::reference_wrapper<const double> >("argDC5", boost::cref(argDC5)));
    ModelParamMap.insert(std::pair<std::string, boost::reference_wrapper<const double> >("argDC6", boost::cref(argDC6)));
    ModelParamMap.insert(std::pair<std::string, boost::reference_wrapper<const double> >("argDC7", boost::cref(argDC7)));
    ModelParamMap.insert(std::pair<std::string, boost::reference_wrapper<const double> >("argDC8", boost::cref(argDC8)));
    ModelParamMap.insert(std::pair<std::string, boost::reference_wrapper<const double> >("argDC9", boost::cref(argDC9)));
    ModelParamMap.insert(std::pair<std::string, boost::reference_wrapper<const double> >("argDC10", boost::cref(argDC10)));
    ModelParamMap.insert(std::pair<std::string, boost::reference_wrapper<const double> >("argDC7g", boost::cref(argDC7g)));
    ModelParamMap.insert(std::pair<std::string, boost::reference_wrapper<const double> >("argDC8g", boost::cref(argDC8g)));
    

    ModelParamMap.insert(std::pair<std::string, boost::reference_wrapper<const double> >("absDC7p", boost::cref(absDC7p)));
    ModelParamMap.insert(std::pair<std::string, boost::reference_wrapper<const double> >("absDC8p", boost::cref(absDC8p)));
    ModelParamMap.insert(std::pair<std::string, boost::reference_wrapper<const double> >("absDC9p", boost::cref(absDC9p)));
    ModelParamMap.insert(std::pair<std::string, boost::reference_wrapper<const double> >("absDC10p", boost::cref(absDC10p)));
    ModelParamMap.insert(std::pair<std::string, boost::reference_wrapper<const double> >("absDC7gp", boost::cref(absDC7gp)));
    ModelParamMap.insert(std::pair<std::string, boost::reference_wrapper<const double> >("absDC8gp", boost::cref(absDC8gp)));
    
    ModelParamMap.insert(std::pair<std::string, boost::reference_wrapper<const double> >("argDC7p", boost::cref(argDC7p)));
    ModelParamMap.insert(std::pair<std::string, boost::reference_wrapper<const double> >("argDC8p", boost::cref(argDC8p)));
    ModelParamMap.insert(std::pair<std::string, boost::reference_wrapper<const double> >("argDC9p", boost::cref(argDC9p)));
    ModelParamMap.insert(std::pair<std::string, boost::reference_wrapper<const double> >("argDC10p", boost::cref(argDC10p)));
    ModelParamMap.insert(std::pair<std::string, boost::reference_wrapper<const double> >("argDC7gp", boost::cref(argDC7gp)));
    ModelParamMap.insert(std::pair<std::string, boost::reference_wrapper<const double> >("argDC8gp", boost::cref(argDC8gp)));
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
    
    DC1 = gslpp::complex(absDC1, argDC1, true);
    DC2 = gslpp::complex(absDC2, argDC2, true);
    DC3 = gslpp::complex(absDC3, argDC3, true);
    DC4 = gslpp::complex(absDC4, argDC4, true);
    DC5 = gslpp::complex(absDC5, argDC5, true);
    DC6 = gslpp::complex(absDC6, argDC6, true);
    DC7 = gslpp::complex(absDC7, argDC7, true);
    DC8 = gslpp::complex(absDC8, argDC8, true);
    DC9 = gslpp::complex(absDC9, argDC8, true);
    DC10 = gslpp::complex(absDC10, argDC10, true);
    DC7g = gslpp::complex(absDC7g, argDC7g, true);
    DC8g = gslpp::complex(absDC8g, argDC8g, true);
    
    DC7p = gslpp::complex(absDC7p, argDC7p, true);
    DC8p = gslpp::complex(absDC8p, argDC8p, true);
    DC9p = gslpp::complex(absDC9p, argDC8p, true);
    DC10p = gslpp::complex(absDC10p, argDC10p, true);
    DC7gp = gslpp::complex(absDC7gp, argDC7gp, true);
    DC8gp = gslpp::complex(absDC8gp, argDC8gp, true);

    /* Necessary for updating StandardModel parameters in StandardModelMatching,
     * and FlavourWC and FlavourWC-derived parameters in FlavourWCMatching */
    myFlavourWilsonCoefficientMatching->updateFlavourWilsonCoefficientParameters();

    return (true);
}

void FlavourWilsonCoefficient::setParameter(const std::string name, const double& value){    
    if(name.compare("absDC1") == 0) 
        absDC1 = value;
    else if(name.compare("absDC2") == 0) 
        absDC2 = value;
    else if(name.compare("absDC3") == 0)
        absDC3 = value;
    else if(name.compare("absDC4") == 0)
        absDC4 = value;
    else if(name.compare("absDC5") == 0)
        absDC5 = value;
    else if(name.compare("absDC6") == 0)
        absDC6 = value;
    else if(name.compare("absDC7") == 0)
        absDC7 = value;
    else if(name.compare("absDC8") == 0)
        absDC8 = value;
    else if(name.compare("absDC9") == 0)
        absDC9 = value;
    else if(name.compare("absDC10") == 0)
        absDC10 = value;
    else if(name.compare("absDC7g") == 0)
        absDC7g = value;
    else if(name.compare("absDC8g") == 0)
        absDC8g = value;
    else if(name.compare("argDC1") == 0) 
        argDC1 = value;
    else if(name.compare("argDC2") == 0) 
        argDC2 = value;
    else if(name.compare("argDC3") == 0)
        argDC3 = value;
    else if(name.compare("argDC4") == 0)
        argDC4 = value;
    else if(name.compare("argDC5") == 0)
        argDC5 = value;
    else if(name.compare("argDC6") == 0)
        argDC6 = value;
    else if(name.compare("argDC7") == 0)
        argDC7 = value;
    else if(name.compare("argDC8") == 0)
        argDC8 = value;
    else if(name.compare("argDC9") == 0)
        argDC9 = value;
    else if(name.compare("argDC10") == 0)
        argDC10 = value;
    else if(name.compare("argDC7g") == 0)
        argDC7g = value;
    else if(name.compare("argDC8g") == 0)
        argDC8g = value;
    else if(name.compare("absDC7p") == 0)
        absDC7p = value;
    else if(name.compare("absDC8p") == 0)
        absDC8p = value;
    else if(name.compare("absDC9p") == 0)
        absDC9p = value;
    else if(name.compare("absDC10p") == 0)
        absDC10p = value;
    else if(name.compare("absDC7gp") == 0)
        absDC7gp = value;
    else if(name.compare("absDC8gp") == 0)
        absDC8gp = value;
    else if(name.compare("argDC7p") == 0)
        argDC7p = value;
    else if(name.compare("argDC8p") == 0)
        argDC8p = value;
    else if(name.compare("argDC9p") == 0)
        argDC9p = value;
    else if(name.compare("argDC10p") == 0)
        argDC10p = value;
    else if(name.compare("argDC7gp") == 0)
        argDC7gp = value;
    else if(name.compare("argDC8gp") == 0)
        argDC8gp = value;
    else
        StandardModel::setParameter(name,value);
}

bool FlavourWilsonCoefficient::CheckParameters(const std::map<std::string, double>& DPars) {
    for (int i = 0; i < NFlavourWilsonCoefficientvars; i++) {
        if (DPars.find(FlavourWilsonCoefficientvars[i]) == DPars.end()) {
            std::cout << "missing mandatory FlavourWilsonCoefficient parameter " << FlavourWilsonCoefficientvars[i] << std::endl;
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
