/* 
 * Copyright (C) 2012 HEPfit Collaboration
 * All rights reserved.
 *
 * For the licensing terms see doc/COPYING.
 */

#include "StandardModelMatching.h"
#include "THDM.h"

const std::string THDM::THDMvars[NTHDMvars] = {"logtb","bma","mHh2","mHh2mmA2","mHh2mmHp2","m12_2","bsgamma_theoryerror"};

THDM::THDM() : StandardModel() {
    
    ModelParamMap.insert(std::pair<std::string, boost::reference_wrapper<const double> >("logtb", boost::cref(logtb)));
    ModelParamMap.insert(std::pair<std::string, boost::reference_wrapper<const double> >("bma", boost::cref(bma)));
    ModelParamMap.insert(std::pair<std::string, boost::reference_wrapper<const double> >("mHh2", boost::cref(mHh2)));
    ModelParamMap.insert(std::pair<std::string, boost::reference_wrapper<const double> >("mHh2mmA2", boost::cref(mHh2mmA2)));
    ModelParamMap.insert(std::pair<std::string, boost::reference_wrapper<const double> >("mHh2mmHp2", boost::cref(mHh2mmHp2)));
    ModelParamMap.insert(std::pair<std::string, boost::reference_wrapper<const double> >("m12_2", boost::cref(m12_2)));
    ModelParamMap.insert(std::pair<std::string, boost::reference_wrapper<const double> >("bsgamma_theoryerror", boost::cref(bsgamma_theoryerror)));
}

THDM::~THDM(){
    if (IsModelInitialized()) {
            if (myTHDMMatching != NULL) delete(myTHDMMatching);
        }
}

///////////////////////////////////////////////////////////////////////////
// Initialization

bool THDM::InitializeModel()
{
    myTHDMMatching = new THDMMatching(*this);
    setModelInitialized(StandardModel::InitializeModel());
    setModelTHDM();
    return(true);
}

bool THDM::Init(const std::map<std::string, double>& DPars) {
    return(StandardModel::Init(DPars));
}

bool THDM::PreUpdate()
{
    if(!StandardModel::PreUpdate()) return (false);

    return (true);
}

bool THDM::Update(const std::map<std::string, double>& DPars) {

    if(!PreUpdate()) return (false);

    UpdateError = false;

    for (std::map<std::string, double>::const_iterator it = DPars.begin(); it != DPars.end(); it++)
        setParameter(it->first, it->second);

    if (UpdateError) return (false);

    if(!PostUpdate()) return (false);

    return (true);
}

bool THDM::PostUpdate()
{
    if(!StandardModel::PostUpdate()) return (false);

    /* Necessary for updating StandardModel parameters in StandardModelMatching,
     * and THDM and THDM-derived parameters in THDMMatching */
    myTHDMMatching->StandardModelMatching::updateSMParameters();
//    myTHDMMatching->updateTHDMParameters();

    return (true);
}

void THDM::setParameter(const std::string name, const double& value){

    if(name.compare("logtb") == 0) {
        logtb = value;
        tanb = pow(10.,logtb);
        if(tanb > 0.) {
            sinb = tanb / sqrt(1. + tanb*tanb);
            cosb = 1. / sqrt(1. + tanb*tanb);
        }
        else {
            throw std::runtime_error("error in THDM::SetParameter, tanb < 0!"); 
          }
        }
    else if(name.compare("bma") == 0) {
        bma = value;
        sin_ba = sin(bma);
    }
    else if(name.compare("mHh2") == 0)
        mHh2 = value;
    else if(name.compare("mHh2mmA2") == 0)
        mHh2mmA2 = value;
    else if(name.compare("mHh2mmHp2") == 0)
        mHh2mmHp2 = value;
    else if(name.compare("m12_2") == 0)
        m12_2 = value;
    else if(name.compare("bsgamma_theoryerror") == 0)
        bsgamma_theoryerror = value;
    else
        StandardModel::setParameter(name,value);
}

bool THDM::CheckParameters(const std::map<std::string, double>& DPars) {
    for (int i = 0; i < NTHDMvars; i++) {
        if (DPars.find(THDMvars[i]) == DPars.end()) {
            std::cout << "missing mandatory THDM parameter " << THDMvars[i] << std::endl;
            return false;
        }
    }
    return(StandardModel::CheckParameters(DPars));
}

///////////////////////////////////////////////////////////////////////////
// Flags

bool THDM::setFlagStr(const std::string name, const std::string value)
{
    bool res = false;
    if(name.compare("modelTypeflag") == 0)
    {
        if (checkmodelType(value))
        {
            flag_model = value;
            res = true;
        }
        else
        {    
            throw std::runtime_error("THDM::setFlagStr(): Invalid flag "
                + name + "=" + value);
        }
    }
    else
    {
        res = StandardModel::setFlagStr(name,value);
    }

    return(res);
}
