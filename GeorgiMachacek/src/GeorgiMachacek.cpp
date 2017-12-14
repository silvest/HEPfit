/* 
 * Copyright (C) 2015 SusyFit Collaboration
 * 
 *
 * For the licensing terms see doc/COPYING.
 */

//#include <StandardModelMatching.h>
#include "GeorgiMachacek.h"
#include "GMcache.h"

std::string GeorgiMachacek::GMvars[NGMvars] = {"logtb", "alpha", "mHh", "mA", "mH5", "Mu1", "Mu2", "Q_GM"};

GeorgiMachacek::GeorgiMachacek() : StandardModel(), GMM(*this) {
    SMM.setObj((StandardModelMatching&) GMM.getObj());
    ModelParamMap.insert(std::pair<std::string, boost::reference_wrapper<const double> >("logtb", boost::cref(logtb)));
    ModelParamMap.insert(std::pair<std::string, boost::reference_wrapper<const double> >("alpha", boost::cref(alpha)));
    ModelParamMap.insert(std::pair<std::string, boost::reference_wrapper<const double> >("mHh", boost::cref(mHh)));
    ModelParamMap.insert(std::pair<std::string, boost::reference_wrapper<const double> >("mA", boost::cref(mA)));
    ModelParamMap.insert(std::pair<std::string, boost::reference_wrapper<const double> >("mH5", boost::cref(mH5)));
    ModelParamMap.insert(std::pair<std::string, boost::reference_wrapper<const double> >("Mu1", boost::cref(Mu1)));
    ModelParamMap.insert(std::pair<std::string, boost::reference_wrapper<const double> >("Mu2", boost::cref(Mu2)));
    ModelParamMap.insert(std::pair<std::string, boost::reference_wrapper<const double> >("Q_GM", boost::cref(Q_GM)));
    flag_use_sq_masses=false;
}

GeorgiMachacek::~GeorgiMachacek(){
    if (IsModelInitialized()) {
            if (myGMcache != NULL) delete(myGMcache);
        }
}

///////////////////////////////////////////////////////////////////////////
// Initialization

bool GeorgiMachacek::InitializeModel()
{
    myGMcache = new GMcache(*this);
    setModelInitialized(StandardModel::InitializeModel());
    setModelGeorgiMachacek();
    return(true);
}
    
bool GeorgiMachacek::Init(const std::map<std::string, double>& DPars) {
    return(StandardModel::Init(DPars));
}

bool GeorgiMachacek::PreUpdate()
{    
    if(!StandardModel::PreUpdate()) return (false);

    return (true);
}

bool GeorgiMachacek::Update(const std::map<std::string, double>& DPars) {
    
    if(!PreUpdate()) return (false);
    
    UpdateError = false;

    for (std::map<std::string, double>::const_iterator it = DPars.begin(); it != DPars.end(); it++)
        setParameter(it->first, it->second);

    if (UpdateError) return (false);

    if(!PostUpdate()) return (false);

    return (true);

}

bool GeorgiMachacek::PostUpdate()
{
    if(!StandardModel::PostUpdate()) return (false);

    mHl2=myGMcache->updateCache();

    return (true);
}

void GeorgiMachacek::setParameter(const std::string name, const double& value){    
    if(name.compare("logtb") == 0) {
        logtb = value;
        tanb = pow(10.,logtb);
        if(tanb > 0.) {
            sinb = tanb / sqrt(1. + tanb*tanb);
            cosb = 1. / sqrt(1. + tanb*tanb);
        }
        else {
            throw std::runtime_error("error in GeorgiMachacek::SetParameter, tanb < 0!"); 
          }
        }
    else if(name.compare("alpha") == 0) {
        alpha = value;
        sin_ba = sinb*cos(alpha)-cosb*sin(alpha);
    }
    else if(name.compare("mHh") == 0 && !flag_use_sq_masses)
        mHh = value;
    else if(name.compare("mA") == 0 && !flag_use_sq_masses)
        mA = value;
    else if(name.compare("mH5") == 0 && !flag_use_sq_masses)
        mH5 = value;
    else if(name.compare("mHhsq") == 0 && flag_use_sq_masses)
        mHhsq = value;
    else if(name.compare("mAsq") == 0 && flag_use_sq_masses)
        mAsq = value;
    else if(name.compare("mH5sq") == 0 && flag_use_sq_masses)
        mH5sq = value;
    else if(name.compare("Mu1") == 0)
        Mu1 = value;
    else if(name.compare("Mu2") == 0)
        Mu2 = value;
    else if(name.compare("Q_GM") == 0)
       Q_GM = value;
    else
        StandardModel::setParameter(name,value);
}

bool GeorgiMachacek::CheckParameters(const std::map<std::string, double>& DPars) {
    for (int i = 0; i < NGMvars; i++) {
        if (DPars.find(GMvars[i]) == DPars.end()) {
            std::cout << "missing mandatory GeorgiMachacek parameter " << GMvars[i] << std::endl;
            return false;
        }
    }
    return(StandardModel::CheckParameters(DPars));
}

bool GeorgiMachacek::setFlag(const std::string name, const bool value)
{
    bool res = false;
    if(name.compare("use_sq_masses") == 0) {
        flag_use_sq_masses = value;
        res = true;
        if (flag_use_sq_masses) {
            GMvars[std::distance(GMvars,std::find(GMvars,GMvars+NGMvars,"mHh"))] = "mHhsq";
            GMvars[std::distance(GMvars,std::find(GMvars,GMvars+NGMvars,"mA"))] = "mAsq";
            GMvars[std::distance(GMvars,std::find(GMvars,GMvars+NGMvars,"mH5"))] = "mH5sq";
        }
    }
    else
        res = StandardModel::setFlag(name,value);

    return(res);
}

//
//double GeorgiMachacek::sign(const double x) const {
//    return ( (x > 0) ? 1 : ((x < 0) ? -1 : 0) );
//}
