/* 
 * Copyright (C) 2012 HEPfit Collaboration
 * All rights reserved.
 *
 * For the licensing terms see doc/COPYING.
 */

#include <StandardModelMatching.h>
#include "THDM.h"

const std::string THDM::THDMvars[NTHDMvars] = {"logtb","bma","mHh","mA","mHp","m12_2","lambda6","lambda7","modelType"};

THDM::THDM() : StandardModel() {   
//    mycache = new THDMcache();
    
    ModelParamMap.insert(std::pair<std::string, boost::reference_wrapper<const double> >("logtb", boost::cref(logtb)));
    ModelParamMap.insert(std::pair<std::string, boost::reference_wrapper<const double> >("bma", boost::cref(bma)));
    ModelParamMap.insert(std::pair<std::string, boost::reference_wrapper<const double> >("mHh", boost::cref(mHh)));
    ModelParamMap.insert(std::pair<std::string, boost::reference_wrapper<const double> >("mA", boost::cref(mA)));
    ModelParamMap.insert(std::pair<std::string, boost::reference_wrapper<const double> >("mHp", boost::cref(mHp)));
    ModelParamMap.insert(std::pair<std::string, boost::reference_wrapper<const double> >("m12_2", boost::cref(m12_2)));
    ModelParamMap.insert(std::pair<std::string, boost::reference_wrapper<const double> >("lambda6", boost::cref(lambda6)));
    ModelParamMap.insert(std::pair<std::string, boost::reference_wrapper<const double> >("lambda7", boost::cref(lambda7)));
    ModelParamMap.insert(std::pair<std::string, boost::reference_wrapper<const double> >("modelType", boost::cref(modelType)));
}

THDM::~THDM(){
    if (IsModelInitialized()) {
            if (myTHDMMatching != NULL) delete(myTHDMMatching);
//            if (mycache != NULL) delete(mycache);
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
    requireCKM = false;
    requireYe = false;
    requireYn = false;

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

    return (true);
}

bool THDM::PostUpdate()
{
    if(!StandardModel::PostUpdate()) return (false);
    
    //In THDM U couple with v2, D with v1 and L with v1
    if (requireYu || requireCKM) {
        Yu = matrix<complex>::Id(3);
        for (int i = 0; i < 3; i++)
            Yu.assign(i, i, this->quarks[UP + 2 * i].getMass() / v2() * sqrt(2.));
        Yu = VCKM.transpose()*Yu;
    }
    if (requireYd) {
        for (int i = 0; i < 3; i++)
            Yd.assign(i, i, this->QCD::quarks[DOWN + 2 * i].getMass() / v1() * sqrt(2.));
    }
    if (requireYe) {
        for (int i = 0; i < 3; i++)
            Ye.assign(i, i, this->leptons[ELECTRON + 2 * i].getMass() / v1() * sqrt(2.));
    }
    if (requireYn) {
        Yn = matrix<complex>::Id(3);
        for (int i = 0; i < 3; i++)
            Yn.assign(i, i, this->leptons[NEUTRINO_1 + 2 * i].getMass() / v1() * sqrt(2.));
        Yn = Yn * UPMNS.hconjugate();
    }

    /* Necessary for updating StandardModel parameters in StandardModelMatching,
     * and THDM and THDM-derived parameters in THDMMatching */
    myTHDMMatching->StandardModelMatching::updateSMParameters();
    myTHDMMatching->updateTHDMParameters();

    return (true);
}

double THDM::computeCosa() const
{
    return cos(atan(pow(10.,logtb))-bma);
}

double THDM::computeSina() const
{
    return sin(atan(pow(10.,logtb))-bma);
}

void THDM::setParameter(const std::string name, const double& value){

    if(name.compare("logtb") == 0) {
        logtb = value;
//        if(tanb > 0.){
        tanb = pow(10.,logtb);
        sinb = tanb / sqrt(1. + tanb*tanb);
        cosb = 1. / sqrt(1. + tanb*tanb);}
//        else {
//            throw std::runtime_error("error in THDM::SetParameter, tanb < 0!"); 
//          }
//        } 
    else if(name.compare("bma") == 0) {
        bma = value;
        sin_ba = sin(bma);
    }
    else if(name.compare("mHh") == 0)
        mHh = value;
    else if(name.compare("mA") == 0)
        mA = value;
    else if(name.compare("mHp") == 0)
        mHp = value;
    else if(name.compare("m12_2") == 0)
        m12_2 = value;
    else if(name.compare("lambda6") == 0)
        lambda6 = value;
    else if(name.compare("lambda7") == 0)
        lambda7 = value;
    else if(name.compare("modelType") == 0)
        modelType = value;
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

bool THDM::setFlag(const std::string name, const bool value)
{
    bool res = false;

    res = StandardModel::setFlag(name,value);

    return(res);
}

double THDM::v1() const {
    return v() * cosb;
}

double THDM::v2() const {
    return v() * sinb;
}
