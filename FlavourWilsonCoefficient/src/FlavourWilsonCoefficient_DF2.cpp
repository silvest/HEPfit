/* 
 * Copyright (C) 2019 HEPfit Collaboration
 *
 *
 * For the licensing terms see doc/COPYING.
 */

#include "FlavourWilsonCoefficient_DF2.h"

const std::string FlavourWilsonCoefficient_DF2::FlavourWilsonCoefficient_DF2vars[NFlavourWilsonCoefficient_DF2vars] = {
    "reC1_s","reC2_s","reC3_s","reC4_s","reC5_s",
"reC1_c","reC2_c","reC3_c","reC4_c","reC5_c","reC1_bd","reC2_bd","reC3_bd","reC4_bd","reC5_bd",
"reC1_bs","reC2_bs","reC3_bs","reC4_bs","reC5_bs","imC1_s","imC2_s","imC3_s","imC4_s","imC5_s",
"imC1_c","imC2_c","imC3_c","imC4_c","imC5_c","imC1_bd","imC2_bd","imC3_bd","imC4_bd","imC5_bd",
"imC1_bs","imC2_bs","imC3_bs","imC4_bs","imC5_bs","WCscale_s","WCscale_c","WCscale_bd","WCscale_bs"};

FlavourWilsonCoefficient_DF2::FlavourWilsonCoefficient_DF2() : StandardModel(), FWCM(*this) {   

    SMM.setObj((StandardModelMatching&) FWCM.getObj());
    ModelParamMap.insert(std::make_pair("reC1_s", std::cref(reC1_s)));
    ModelParamMap.insert(std::make_pair("reC2_s", std::cref(reC2_s)));
    ModelParamMap.insert(std::make_pair("reC3_s", std::cref(reC3_s)));
    ModelParamMap.insert(std::make_pair("reC4_s", std::cref(reC4_s)));
    ModelParamMap.insert(std::make_pair("reC5_s", std::cref(reC5_s)));
    ModelParamMap.insert(std::make_pair("imC1_s", std::cref(imC1_s)));
    ModelParamMap.insert(std::make_pair("imC2_s", std::cref(imC2_s)));
    ModelParamMap.insert(std::make_pair("imC3_s", std::cref(imC3_s)));
    ModelParamMap.insert(std::make_pair("imC4_s", std::cref(imC4_s)));
    ModelParamMap.insert(std::make_pair("imC5_s", std::cref(imC5_s)));    
    ModelParamMap.insert(std::make_pair("WCscale_s", std::cref(WCscale_s)));
    ModelParamMap.insert(std::make_pair("reC1_c", std::cref(reC1_c)));
    ModelParamMap.insert(std::make_pair("reC2_c", std::cref(reC2_c)));
    ModelParamMap.insert(std::make_pair("reC3_c", std::cref(reC3_c)));
    ModelParamMap.insert(std::make_pair("reC4_c", std::cref(reC4_c)));
    ModelParamMap.insert(std::make_pair("reC5_c", std::cref(reC5_c)));
    ModelParamMap.insert(std::make_pair("imC1_c", std::cref(imC1_c)));
    ModelParamMap.insert(std::make_pair("imC2_c", std::cref(imC2_c)));
    ModelParamMap.insert(std::make_pair("imC3_c", std::cref(imC3_c)));
    ModelParamMap.insert(std::make_pair("imC4_c", std::cref(imC4_c)));
    ModelParamMap.insert(std::make_pair("imC5_c", std::cref(imC5_c)));    
    ModelParamMap.insert(std::make_pair("WCscale_c", std::cref(WCscale_c)));
    ModelParamMap.insert(std::make_pair("reC1_bd", std::cref(reC1_bd)));
    ModelParamMap.insert(std::make_pair("reC2_bd", std::cref(reC2_bd)));
    ModelParamMap.insert(std::make_pair("reC3_bd", std::cref(reC3_bd)));
    ModelParamMap.insert(std::make_pair("reC4_bd", std::cref(reC4_bd)));
    ModelParamMap.insert(std::make_pair("reC5_bd", std::cref(reC5_bd)));
    ModelParamMap.insert(std::make_pair("imC1_bd", std::cref(imC1_bd)));
    ModelParamMap.insert(std::make_pair("imC2_bd", std::cref(imC2_bd)));
    ModelParamMap.insert(std::make_pair("imC3_bd", std::cref(imC3_bd)));
    ModelParamMap.insert(std::make_pair("imC4_bd", std::cref(imC4_bd)));
    ModelParamMap.insert(std::make_pair("imC5_bd", std::cref(imC5_bd)));    
    ModelParamMap.insert(std::make_pair("WCscale_bd", std::cref(WCscale_bd)));
    ModelParamMap.insert(std::make_pair("reC1_bs", std::cref(reC1_bs)));
    ModelParamMap.insert(std::make_pair("reC2_bs", std::cref(reC2_bs)));
    ModelParamMap.insert(std::make_pair("reC3_bs", std::cref(reC3_bs)));
    ModelParamMap.insert(std::make_pair("reC4_bs", std::cref(reC4_bs)));
    ModelParamMap.insert(std::make_pair("reC5_bs", std::cref(reC5_bs)));
    ModelParamMap.insert(std::make_pair("imC1_bs", std::cref(imC1_bs)));
    ModelParamMap.insert(std::make_pair("imC2_bs", std::cref(imC2_bs)));
    ModelParamMap.insert(std::make_pair("imC3_bs", std::cref(imC3_bs)));
    ModelParamMap.insert(std::make_pair("imC4_bs", std::cref(imC4_bs)));
    ModelParamMap.insert(std::make_pair("imC5_bs", std::cref(imC5_bs)));    
    ModelParamMap.insert(std::make_pair("WCscale_bs", std::cref(WCscale_bs)));
}


///////////////////////////////////////////////////////////////////////////
// Initialization

bool FlavourWilsonCoefficient_DF2::InitializeModel()
{
    setModelInitialized(StandardModel::InitializeModel());
    return(true);
}
    
bool FlavourWilsonCoefficient_DF2::PostUpdate()
{
    if(!StandardModel::PostUpdate()) return (false);
    
    C1_s = gslpp::complex(reC1_s, imC1_s);
    C2_s = gslpp::complex(reC2_s, imC2_s);
    C3_s = gslpp::complex(reC3_s, imC3_s);
    C4_s = gslpp::complex(reC4_s, imC4_s);
    C5_s = gslpp::complex(reC5_s, imC5_s);
    C1_c = gslpp::complex(reC1_c, imC1_c);
    C2_c = gslpp::complex(reC2_c, imC2_c);
    C3_c = gslpp::complex(reC3_c, imC3_c);
    C4_c = gslpp::complex(reC4_c, imC4_c);
    C5_c = gslpp::complex(reC5_c, imC5_c);
    C1_bd = gslpp::complex(reC1_bd, imC1_bd);
    C2_bd = gslpp::complex(reC2_bd, imC2_bd);
    C3_bd = gslpp::complex(reC3_bd, imC3_bd);
    C4_bd = gslpp::complex(reC4_bd, imC4_bd);
    C5_bd = gslpp::complex(reC5_bd, imC5_bd);
    C1_bs = gslpp::complex(reC1_bs, imC1_bs);
    C2_bs = gslpp::complex(reC2_bs, imC2_bs);
    C3_bs = gslpp::complex(reC3_bs, imC3_bs);
    C4_bs = gslpp::complex(reC4_bs, imC4_bs);
    C5_bs = gslpp::complex(reC5_bs, imC5_bs);

    return (true);
}

void FlavourWilsonCoefficient_DF2::setParameter(const std::string name, const double& value){    
    if(name.compare("reC1_s") == 0) 
        reC1_s = value;
    else if(name.compare("reC2_s") == 0) 
        reC2_s = value;
    else if(name.compare("reC3_s") == 0)
        reC3_s = value;
    else if(name.compare("reC4_s") == 0)
        reC4_s = value;
    else if(name.compare("reC5_s") == 0)
        reC5_s = value;
    else if(name.compare("imC1_s") == 0) 
        imC1_s = value;
    else if(name.compare("imC2_s") == 0) 
        imC2_s = value;
    else if(name.compare("imC3_s") == 0)
        imC3_s = value;
    else if(name.compare("imC4_s") == 0)
        imC4_s = value;
    else if(name.compare("imC5_s") == 0)
        imC5_s = value;
    else if(name.compare("WCscale_s") == 0)
        WCscale_s = value;
    else if(name.compare("reC1_c") == 0) 
        reC1_c = value;
    else if(name.compare("reC2_c") == 0) 
        reC2_c = value;
    else if(name.compare("reC3_c") == 0)
        reC3_c = value;
    else if(name.compare("reC4_c") == 0)
        reC4_c = value;
    else if(name.compare("reC5_c") == 0)
        reC5_c = value;
    else if(name.compare("imC1_c") == 0) 
        imC1_c = value;
    else if(name.compare("imC2_c") == 0) 
        imC2_c = value;
    else if(name.compare("imC3_c") == 0)
        imC3_c = value;
    else if(name.compare("imC4_c") == 0)
        imC4_c = value;
    else if(name.compare("imC5_c") == 0)
        imC5_c = value;
    else if(name.compare("WCscale_c") == 0)
        WCscale_c = value;
    else if(name.compare("reC1_bd") == 0) 
        reC1_bd = value;
    else if(name.compare("reC2_bd") == 0) 
        reC2_bd = value;
    else if(name.compare("reC3_bd") == 0)
        reC3_bd = value;
    else if(name.compare("reC4_bd") == 0)
        reC4_bd = value;
    else if(name.compare("reC5_bd") == 0)
        reC5_bd = value;
    else if(name.compare("imC1_bd") == 0) 
        imC1_bd = value;
    else if(name.compare("imC2_bd") == 0) 
        imC2_bd = value;
    else if(name.compare("imC3_bd") == 0)
        imC3_bd = value;
    else if(name.compare("imC4_bd") == 0)
        imC4_bd = value;
    else if(name.compare("imC5_bd") == 0)
        imC5_bd = value;
    else if(name.compare("WCscale_bd") == 0)
        WCscale_bd = value;
    else if(name.compare("reC1_bs") == 0) 
        reC1_bs = value;
    else if(name.compare("reC2_bs") == 0) 
        reC2_bs = value;
    else if(name.compare("reC3_bs") == 0)
        reC3_bs = value;
    else if(name.compare("reC4_bs") == 0)
        reC4_bs = value;
    else if(name.compare("reC5_bs") == 0)
        reC5_bs = value;
    else if(name.compare("imC1_bs") == 0) 
        imC1_bs = value;
    else if(name.compare("imC2_bs") == 0) 
        imC2_bs = value;
    else if(name.compare("imC3_bs") == 0)
        imC3_bs = value;
    else if(name.compare("imC4_bs") == 0)
        imC4_bs = value;
    else if(name.compare("imC5_bs") == 0)
        imC5_bs = value;
    else if(name.compare("WCscale_bs") == 0)
        WCscale_bs = value;
    else
        StandardModel::setParameter(name,value);
}

bool FlavourWilsonCoefficient_DF2::CheckParameters(const std::map<std::string, double>& DPars) {
    for (int i = 0; i < NFlavourWilsonCoefficient_DF2vars; i++) {
        if (DPars.find(FlavourWilsonCoefficient_DF2vars[i]) == DPars.end()) {
            std::cout << "missing mandatory FlavourWilsonCoefficient_DF2 parameter " << FlavourWilsonCoefficient_DF2vars[i] << std::endl;
            return false;
        }
    }
    return(StandardModel::CheckParameters(DPars));
}
