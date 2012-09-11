/* 
 * File:   NewPhysicsSTU.cpp
 * Author: mishima
 */

#include "NewPhysicsSTU.h"


const std::string NewPhysicsSTU::STUvars[NSTUvars] 
                  = {"obliqueS", "obliqueT", "obliqueU", 
                     "obliqueV", "obliqueW", "obliqueX", "obliqueY"};


NewPhysicsSTU::NewPhysicsSTU() : StandardModel() {
}


bool NewPhysicsSTU::Update(const std::map<std::string,double>& DPars) {
    for (std::map<std::string, double>::const_iterator it = DPars.begin(); it != DPars.end(); it++)
        SetParameter(it->first, it->second);
    if(!StandardModel::Update(DPars)) return (false);
    
    return (true);
}


bool NewPhysicsSTU::Init(const std::map<std::string, double>& DPars) {
    Update(DPars);
    return(CheckParameters(DPars)); 
}


bool NewPhysicsSTU::CheckParameters(const std::map<std::string, double>& DPars) {
    for (int i = 0; i < NSTUvars; i++) {
        if (DPars.find(STUvars[i]) == DPars.end()) {
            std::cout << "missing mandatory NewPhysicsSTU parameter " 
                      << STUvars[i] << std::endl;
            return false;
        }
    }
    return(StandardModel::CheckParameters(DPars));
}

    
void NewPhysicsSTU::SetParameter(const std::string name, const double& value) {
    if (name.compare("obliqueS") == 0)
        myObliqueS = value;
    else if (name.compare("obliqueT") == 0)
        myObliqueT = value;
    else if (name.compare("obliqueU") == 0)
        myObliqueU = value;    
    else if (name.compare("obliqueV") == 0)
        myObliqueU = value;    
    else if (name.compare("obliqueW") == 0)
        myObliqueU = value;    
    else if (name.compare("obliqueX") == 0)
        myObliqueU = value;    
    else if (name.compare("obliqueY") == 0)
        myObliqueU = value;    
    else
        StandardModel::SetParameter(name, value);       
}



