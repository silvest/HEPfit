/* 
 * Copyright (C) 2017 HEPfit Collaboration
 *
 *
 * For the licensing terms see doc/COPYING.
 */

#include "THDMW.h"
//#include "THDMWcache.h"

const std::string THDMW::THDMWvars[NTHDMWvars] = {"THDMW_m11_2","THDMW_m22_2","THDMW_m12_2",
                                               "THDMW_lambda1","THDMW_lambda2","THDMW_lambda3","THDMW_lambda4","THDMW_lambda5",
                                               "THDMW_mS2","THDMW_mu1","THDMW_mu2","THDMW_mu3","THDMW_mu4","THDMW_mu5","THDMW_mu6",
                                               "THDMW_nu1","THDMW_nu2","THDMW_nu3","THDMW_nu4","THDMW_nu5",
                                               "THDMW_omega1","THDMW_omega2","THDMW_omega3","THDMW_omega4","THDMW_omega5",
                                               "THDMW_kappa1","THDMW_kappa2","THDMW_kappa3","Q_THDMW"};

THDMW::THDMW() : StandardModel()/*, THDMWM(*this)*/ {

//    SMM.setObj((StandardModelMatching&) THDMWM.getObj());
    ModelParamMap.insert(std::pair<std::string, boost::reference_wrapper<const double> >("THDMW_m11_2", boost::cref(THDMW_m11_2)));
    ModelParamMap.insert(std::pair<std::string, boost::reference_wrapper<const double> >("THDMW_m22_2", boost::cref(THDMW_m22_2)));
    ModelParamMap.insert(std::pair<std::string, boost::reference_wrapper<const double> >("THDMW_m12_2", boost::cref(THDMW_m12_2)));
    ModelParamMap.insert(std::pair<std::string, boost::reference_wrapper<const double> >("THDMW_lambda1", boost::cref(THDMW_lambda1)));
    ModelParamMap.insert(std::pair<std::string, boost::reference_wrapper<const double> >("THDMW_lambda2", boost::cref(THDMW_lambda2)));
    ModelParamMap.insert(std::pair<std::string, boost::reference_wrapper<const double> >("THDMW_lambda3", boost::cref(THDMW_lambda3)));
    ModelParamMap.insert(std::pair<std::string, boost::reference_wrapper<const double> >("THDMW_lambda4", boost::cref(THDMW_lambda4)));
    ModelParamMap.insert(std::pair<std::string, boost::reference_wrapper<const double> >("THDMW_lambda5", boost::cref(THDMW_lambda5)));
    ModelParamMap.insert(std::pair<std::string, boost::reference_wrapper<const double> >("THDMW_mS2", boost::cref(THDMW_mS2)));
    ModelParamMap.insert(std::pair<std::string, boost::reference_wrapper<const double> >("THDMW_mu1", boost::cref(THDMW_mu1)));
    ModelParamMap.insert(std::pair<std::string, boost::reference_wrapper<const double> >("THDMW_mu2", boost::cref(THDMW_mu2)));
    ModelParamMap.insert(std::pair<std::string, boost::reference_wrapper<const double> >("THDMW_mu3", boost::cref(THDMW_mu3)));
    ModelParamMap.insert(std::pair<std::string, boost::reference_wrapper<const double> >("THDMW_mu4", boost::cref(THDMW_mu4)));
    ModelParamMap.insert(std::pair<std::string, boost::reference_wrapper<const double> >("THDMW_mu5", boost::cref(THDMW_mu5)));
    ModelParamMap.insert(std::pair<std::string, boost::reference_wrapper<const double> >("THDMW_mu6", boost::cref(THDMW_mu6)));
    ModelParamMap.insert(std::pair<std::string, boost::reference_wrapper<const double> >("THDMW_nu1", boost::cref(THDMW_nu1)));
    ModelParamMap.insert(std::pair<std::string, boost::reference_wrapper<const double> >("THDMW_nu2", boost::cref(THDMW_nu2)));
    ModelParamMap.insert(std::pair<std::string, boost::reference_wrapper<const double> >("THDMW_nu3", boost::cref(THDMW_nu3)));
    ModelParamMap.insert(std::pair<std::string, boost::reference_wrapper<const double> >("THDMW_nu4", boost::cref(THDMW_nu4)));
    ModelParamMap.insert(std::pair<std::string, boost::reference_wrapper<const double> >("THDMW_nu5", boost::cref(THDMW_nu5)));
    ModelParamMap.insert(std::pair<std::string, boost::reference_wrapper<const double> >("THDMW_omega1", boost::cref(THDMW_omega1)));
    ModelParamMap.insert(std::pair<std::string, boost::reference_wrapper<const double> >("THDMW_omega2", boost::cref(THDMW_omega2)));
    ModelParamMap.insert(std::pair<std::string, boost::reference_wrapper<const double> >("THDMW_omega3", boost::cref(THDMW_omega3)));
    ModelParamMap.insert(std::pair<std::string, boost::reference_wrapper<const double> >("THDMW_omega4", boost::cref(THDMW_omega4)));
    ModelParamMap.insert(std::pair<std::string, boost::reference_wrapper<const double> >("THDMW_omega5", boost::cref(THDMW_omega5)));
    ModelParamMap.insert(std::pair<std::string, boost::reference_wrapper<const double> >("THDMW_kappa1", boost::cref(THDMW_kappa1)));
    ModelParamMap.insert(std::pair<std::string, boost::reference_wrapper<const double> >("THDMW_kappa2", boost::cref(THDMW_kappa2)));
    ModelParamMap.insert(std::pair<std::string, boost::reference_wrapper<const double> >("THDMW_kappa3", boost::cref(THDMW_kappa3)));
    ModelParamMap.insert(std::pair<std::string, boost::reference_wrapper<const double> >("Q_THDMW", boost::cref(Q_THDMW)));
}

THDMW::~THDMW(){
    if (IsModelInitialized()) {
//            if (myTHDMWcache != NULL) delete(myTHDMWcache);
        }
}

///////////////////////////////////////////////////////////////////////////
// Initialization

bool THDMW::InitializeModel()
{
//    myTHDMWcache = new THDMWcache(*this);
    setModelInitialized(StandardModel::InitializeModel());
//    setModelTHDMW();
    return(true);
}

bool THDMW::Init(const std::map<std::string, double>& DPars) {
    return(StandardModel::Init(DPars));
}

bool THDMW::PreUpdate()
{
    if(!StandardModel::PreUpdate()) return (false);

    return (true);
}

bool THDMW::Update(const std::map<std::string, double>& DPars) {

    if(!PreUpdate()) return (false);

    UpdateError = false;

    for (std::map<std::string, double>::const_iterator it = DPars.begin(); it != DPars.end(); it++)
        setParameter(it->first, it->second);

    if (UpdateError) return (false);

    if(!PostUpdate()) return (false);

    return (true);
}

bool THDMW::PostUpdate()
{
    if(!StandardModel::PostUpdate()) return (false);

//    mHl2=myTHDMWcache->updateCache();

    return (true);
}

void THDMW::setParameter(const std::string name, const double& value){

    if(name.compare("THDMW_m11_2") == 0)
        THDMW_m11_2 = value;
    else if(name.compare("THDMW_m22_2") == 0)
        THDMW_m22_2 = value;
    else if(name.compare("THDMW_m12_2") == 0)
        THDMW_m12_2 = value;
    else if(name.compare("THDMW_lambda1") == 0)
        THDMW_lambda1 = value;
    else if(name.compare("THDMW_lambda2") == 0)
        THDMW_lambda2 = value;
    else if(name.compare("THDMW_lambda3") == 0)
        THDMW_lambda3 = value;
    else if(name.compare("THDMW_lambda4") == 0)
        THDMW_lambda4 = value;
    else if(name.compare("THDMW_lambda5") == 0)
        THDMW_lambda5 = value;
    else if(name.compare("THDMW_mS2") == 0)
        THDMW_mS2 = value;
    else if(name.compare("THDMW_mu1") == 0)
        THDMW_mu1 = value;
    else if(name.compare("THDMW_mu2") == 0)
        THDMW_mu2 = value;
    else if(name.compare("THDMW_mu3") == 0)
        THDMW_mu3 = value;
    else if(name.compare("THDMW_mu4") == 0)
        THDMW_mu4 = value;
    else if(name.compare("THDMW_mu5") == 0)
        THDMW_mu5 = value;
    else if(name.compare("THDMW_mu6") == 0)
        THDMW_mu6 = value;
    else if(name.compare("THDMW_nu1") == 0)
        THDMW_nu1 = value;
    else if(name.compare("THDMW_nu2") == 0)
        THDMW_nu2 = value;
    else if(name.compare("THDMW_nu3") == 0)
        THDMW_nu3 = value;
    else if(name.compare("THDMW_nu4") == 0)
        THDMW_nu4 = value;
    else if(name.compare("THDMW_nu5") == 0)
        THDMW_nu5 = value;
    else if(name.compare("THDMW_omega1") == 0)
        THDMW_omega1 = value;
    else if(name.compare("THDMW_omega2") == 0)
        THDMW_omega2 = value;
    else if(name.compare("THDMW_omega3") == 0)
        THDMW_omega3 = value;
    else if(name.compare("THDMW_omega4") == 0)
        THDMW_omega4 = value;
    else if(name.compare("THDMW_omega5") == 0)
        THDMW_omega5 = value;
    else if(name.compare("THDMW_kappa1") == 0)
        THDMW_kappa1 = value;
    else if(name.compare("THDMW_kappa2") == 0)
        THDMW_kappa2 = value;
    else if(name.compare("THDMW_kappa3") == 0)
        THDMW_kappa3 = value;
    else if(name.compare("Q_THDMW") == 0)
        Q_THDMW = value;
    else
        StandardModel::setParameter(name,value);
}

bool THDMW::CheckParameters(const std::map<std::string, double>& DPars) {
    for (int i = 0; i < NTHDMWvars; i++) {
        if (DPars.find(THDMWvars[i]) == DPars.end()) {
            std::cout << "missing mandatory THDMW parameter " << THDMWvars[i] << std::endl;
            return false;
        }
    }
    return(StandardModel::CheckParameters(DPars));
}

///////////////////////////////////////////////////////////////////////////
// Flags
