/* 
 * File:   NewPhysicsSTU.cpp
 * Author: mishima
 */

#include "NewPhysicsSTU.h"


const std::string NewPhysicsSTU::STUvars[NSTUvars] 
                  = {"obliqueS", "obliqueT", "obliqueU"};


NewPhysicsSTU::NewPhysicsSTU() : StandardModel() {
}


void NewPhysicsSTU::Update(const std::map<std::string,double>& DPars) {
    for (std::map<std::string, double>::const_iterator it = DPars.begin(); it != DPars.end(); it++)
        SetParameter(it->first, it->second);
    StandardModel::Update(DPars);
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
        obliqueS = value;
    else if (name.compare("obliqueT") == 0)
        obliqueT = value;
    else if (name.compare("obliqueU") == 0)
        obliqueU = value;    
    else
        StandardModel::SetParameter(name, value);       
}



