/* 
 * Copyright (C) 2016 HEPfit Collaboration
 *
 * For the licensing terms see doc/COPYING.
 */

#include "NPDF2.h"

const std::string NPDF2::NPDF2vars[NNPDF2vars] = {"CepsK","CDMK","CBd","PhiBd","CBs","PhiBs"};

NPDF2::NPDF2() : StandardModel() {   

    ModelParamMap.insert(std::pair<std::string, boost::reference_wrapper<const double> >("CepsK", boost::cref(CepsK)));
    ModelParamMap.insert(std::pair<std::string, boost::reference_wrapper<const double> >("CDMK", boost::cref(CDMK)));
    ModelParamMap.insert(std::pair<std::string, boost::reference_wrapper<const double> >("CBd", boost::cref(CBd)));
    ModelParamMap.insert(std::pair<std::string, boost::reference_wrapper<const double> >("PhiBd", boost::cref(phiBd)));
    ModelParamMap.insert(std::pair<std::string, boost::reference_wrapper<const double> >("CBs", boost::cref(CBs)));
    ModelParamMap.insert(std::pair<std::string, boost::reference_wrapper<const double> >("PhiBs", boost::cref(PhiBs)));
}

void NPDF2::setParameter(const std::string name, const double& value){    
    if(name.compare("CepsK") == 0) 
        CepsK = value;
    else if(name.compare("CDMK") == 0) 
        CDMK = value;
    else if(name.compare("CBd") == 0)
        CBd = value;
    else if(name.compare("PhiBd") == 0)
        phiBd = value;
    else if(name.compare("CBs") == 0)
        CBs = value;
    else if(name.compare("PhiBs") == 0)
        PhiBs = value;
    else
        StandardModel::setParameter(name,value);
}

bool NPDF2::CheckParameters(const std::map<std::string, double>& DPars) {
    for (int i = 0; i < NNPDF2vars; i++) {
        if (DPars.find(NPDF2vars[i]) == DPars.end()) {
            std::cout << "missing mandatory NPDF2 parameter " << NPDF2vars[i] << std::endl;
            return false;
        }
    }
    return(StandardModel::CheckParameters(DPars));
}
