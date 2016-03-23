/* 
 * Copyright (C) 2016 HEPfit Collaboration
 * All rights reserved.
 *
 * For the licensing terms see doc/COPYING.
 */

#include "GeneralTHDM.h"
//#include "GeneralTHDMcache.h"

const std::string GeneralTHDM::GeneralTHDMvars[NGeneralTHDMvars] = {"logtb","mHp2","alpha1","alpha2","alpha3","m2_2","Rem12_2","Relambda6","Imlambda6","Relambda7","Imlambda7","abssigma_u","phi_u","abssigma_d","phi_d","abssigma_l","phi_l","Q_GTHDM"};

GeneralTHDM::GeneralTHDM() : StandardModel()/*, THDMM(*this)*/ {

//    SMM.setObj((StandardModelMatching&) THDMM.getObj());
    ModelParamMap.insert(std::pair<std::string, boost::reference_wrapper<const double> >("logtb", boost::cref(logtb)));
    ModelParamMap.insert(std::pair<std::string, boost::reference_wrapper<const double> >("mHp2", boost::cref(mHp2)));
    ModelParamMap.insert(std::pair<std::string, boost::reference_wrapper<const double> >("alpha1", boost::cref(alpha1)));
    ModelParamMap.insert(std::pair<std::string, boost::reference_wrapper<const double> >("alpha2", boost::cref(alpha2)));
    ModelParamMap.insert(std::pair<std::string, boost::reference_wrapper<const double> >("alpha3", boost::cref(alpha3)));
    ModelParamMap.insert(std::pair<std::string, boost::reference_wrapper<const double> >("m2_2", boost::cref(m2_2)));
    ModelParamMap.insert(std::pair<std::string, boost::reference_wrapper<const double> >("Rem12_2", boost::cref(Rem12_2)));
    ModelParamMap.insert(std::pair<std::string, boost::reference_wrapper<const double> >("Relambda6", boost::cref(Relambda6)));
    ModelParamMap.insert(std::pair<std::string, boost::reference_wrapper<const double> >("Imlambda6", boost::cref(Imlambda6)));
    ModelParamMap.insert(std::pair<std::string, boost::reference_wrapper<const double> >("Relambda7", boost::cref(Relambda7)));
    ModelParamMap.insert(std::pair<std::string, boost::reference_wrapper<const double> >("Imlambda7", boost::cref(Imlambda7)));
    ModelParamMap.insert(std::pair<std::string, boost::reference_wrapper<const double> >("abssigma_u", boost::cref(abssigma_u)));
    ModelParamMap.insert(std::pair<std::string, boost::reference_wrapper<const double> >("phi_u", boost::cref(phi_u)));
    ModelParamMap.insert(std::pair<std::string, boost::reference_wrapper<const double> >("abssigma_d", boost::cref(abssigma_d)));
    ModelParamMap.insert(std::pair<std::string, boost::reference_wrapper<const double> >("phi_d", boost::cref(phi_d)));
    ModelParamMap.insert(std::pair<std::string, boost::reference_wrapper<const double> >("abssigma_l", boost::cref(abssigma_l)));
    ModelParamMap.insert(std::pair<std::string, boost::reference_wrapper<const double> >("phi_l", boost::cref(phi_l)));
    ModelParamMap.insert(std::pair<std::string, boost::reference_wrapper<const double> >("Q_GTHDM", boost::cref(Q_GTHDM)));
}

GeneralTHDM::~GeneralTHDM(){
//    if (IsModelInitialized()) {
//            if (myTHDMcache != NULL) delete(myTHDMcache);
//        }
}

///////////////////////////////////////////////////////////////////////////
// Initialization

bool GeneralTHDM::InitializeModel()
{
//    myTHDMcache = new THDMcache(*this);
//    setModelInitialized(StandardModel::InitializeModel());
//    setModelTHDM();
    return(true);
}

bool GeneralTHDM::Init(const std::map<std::string, double>& DPars) {
    return(StandardModel::Init(DPars));
}

bool GeneralTHDM::PreUpdate()
{
    if(!StandardModel::PreUpdate()) return (false);

    return (true);
}

bool GeneralTHDM::Update(const std::map<std::string, double>& DPars) {

    if(!PreUpdate()) return (false);

    UpdateError = false;

    for (std::map<std::string, double>::const_iterator it = DPars.begin(); it != DPars.end(); it++)
        setParameter(it->first, it->second);

    if (UpdateError) return (false);

    if(!PostUpdate()) return (false);

    return (true);
}

bool GeneralTHDM::PostUpdate()
{
    if(!StandardModel::PostUpdate()) return (false);

//    myTHDMcache->updateCache();

    return (true);
}

void GeneralTHDM::setParameter(const std::string name, const double& value){

    if(name.compare("logtb") == 0) {
        logtb = value;
        tanb = pow(10.,logtb);
        if(tanb > 0.) {
            sinb = tanb / sqrt(1. + tanb*tanb);
            cosb = 1. / sqrt(1. + tanb*tanb);
        }
        else {
            throw std::runtime_error("error in GeneralTHDM::SetParameter, tanb < 0!");
          }
        }
    else if(name.compare("mHp2") == 0)
        mHp2 = value;
    else if(name.compare("alpha1") == 0)
        alpha1 = value;
    else if(name.compare("alpha2") == 0)
        alpha2 = value;
    else if(name.compare("alpha3") == 0)
        alpha3 = value;
    else if(name.compare("m2_2") == 0)
        m2_2 = value;
    else if(name.compare("Rem12_2") == 0)
        Rem12_2 = value;
    else if(name.compare("Relambda6") == 0)
        Relambda6 = value;
    else if(name.compare("Imlambda6") == 0)
        Imlambda6 = value;
    else if(name.compare("Relambda7") == 0)
        Relambda7 = value;
    else if(name.compare("Imlambda7") == 0)
        Imlambda7 = value;
    else if(name.compare("abssigma_u") == 0)
        abssigma_u = value;
    else if(name.compare("phi_u") == 0)
        phi_u = value;
    else if(name.compare("abssigma_d") == 0)
        abssigma_d = value;
    else if(name.compare("phi_d") == 0)
        phi_d = value;
    else if(name.compare("abssigma_l") == 0)
        abssigma_l = value;
    else if(name.compare("phi_l") == 0)
        phi_l = value;
    else if(name.compare("Q_GTHDM") == 0)
        Q_GTHDM = value;
    else
        StandardModel::setParameter(name,value);
}

bool GeneralTHDM::CheckParameters(const std::map<std::string, double>& DPars) {
    for (int i = 0; i < NGeneralTHDMvars; i++) {
        if (DPars.find(GeneralTHDMvars[i]) == DPars.end()) {
            std::cout << "missing mandatory GeneralTHDM parameter " << GeneralTHDMvars[i] << std::endl;
            return false;
        }
    }
    return(StandardModel::CheckParameters(DPars));
}

///////////////////////////////////////////////////////////////////////////
// Flags
//
//bool THDM::setFlagStr(const std::string name, const std::string value)
//{
//    bool res = false;
//    if(name.compare("modelTypeflag") == 0)
//    {
//        if (checkmodelType(value))
//        {
//            flag_model = value;
//            res = true;
//        }
//        else
//        {
//            throw std::runtime_error("THDM::setFlagStr(): Invalid flag "
//                + name + "=" + value);
//        }
//    }
//    else
//    {
//        res = StandardModel::setFlagStr(name,value);
//    }
//
//    return(res);
//}
//
//bool THDM::setFlag(const std::string name, const bool value)
//{
//    bool res = false;
//    if(name.compare("wavefunctionrenormalization") == 0) {
//        flag_wfr = value;
//        res = true;
//    }
//    else
//        res = StandardModel::setFlag(name,value);
//
//    return(res);
//}
