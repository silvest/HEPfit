/* 
 * Copyright (C) 2015 SusyFit Collaboration
 * All rights reserved.
 *
 * For the licensing terms see doc/COPYING.
 */

//#include <StandardModelMatching.h>
#include "GeorgiMachacek.h"
#include "GMcache.h"

const std::string GeorgiMachacek::GMvars[NGMvars] = {"logttheta", "alpha", "mHh", "mH3", "mH5", "M1", "M2", "Q_GM"};

GeorgiMachacek::GeorgiMachacek() : StandardModel() {   
//    myGMcache = new GMcache();

    ModelParamMap.insert(std::pair<std::string, boost::reference_wrapper<const double> >("logttheta", boost::cref(logttheta)));
    ModelParamMap.insert(std::pair<std::string, boost::reference_wrapper<const double> >("alpha", boost::cref(alpha)));
    ModelParamMap.insert(std::pair<std::string, boost::reference_wrapper<const double> >("mHh", boost::cref(mHh)));
    ModelParamMap.insert(std::pair<std::string, boost::reference_wrapper<const double> >("mH3", boost::cref(mH3)));
    ModelParamMap.insert(std::pair<std::string, boost::reference_wrapper<const double> >("mH5", boost::cref(mH5)));
    ModelParamMap.insert(std::pair<std::string, boost::reference_wrapper<const double> >("M1", boost::cref(M1)));
    ModelParamMap.insert(std::pair<std::string, boost::reference_wrapper<const double> >("M2", boost::cref(M2)));
    ModelParamMap.insert(std::pair<std::string, boost::reference_wrapper<const double> >("Q_GM", boost::cref(Q_GM)));
}

GeorgiMachacek::~GeorgiMachacek(){
    if (IsModelInitialized()) {
//            if (myGMMatching != NULL) delete(myGMMatching);
            if (myGMcache != NULL) delete(myGMcache);
        }
}

///////////////////////////////////////////////////////////////////////////
// Initialization

bool GeorgiMachacek::InitializeModel()
{
    myGMcache = new GMcache(*this);
//    myGMMatching = new GMMatching(*this);
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

    myGMcache->updateCache();

    return (true);
}

double GeorgiMachacek::computeCosa() const
{
    return cos(alpha);
}

double GeorgiMachacek::computeSina() const
{
    return sin(alpha);
}

double GeorgiMachacek::computelambda1() const
{
    return (mHh*mHh+mHl*mHl
            +(mHh*mHh-mHl*mHl)*sign(mHh-mHl)*(cos(alpha)*cos(alpha)-sin(alpha)*sin(alpha)))/
           (16.*v()*v()*costheta*costheta);
}

double GeorgiMachacek::computelambda2() const
{
    return (mH3*mH3 -M1*v()/(2.*sqrt(2.)*sintheta)
            +((mHh*mHh-mHl*mHl)*sign(mHh-mHl)*cos(alpha)*sin(alpha))/(sqrt(6.)*costheta*sintheta))/
           (v()*v());
}

double GeorgiMachacek::computelambda3() const
{
    return (mH5*mH5 -3.*sqrt(2.)*M2*v()*sintheta
            +(-3.*mH3*mH3 + sqrt(2.)*M1*v()/sintheta)*costheta*costheta)/
           (v()*v()*sintheta*sintheta);
}

double GeorgiMachacek::computelambda4() const
{
    return (-2.*mH5*mH5 +9.*sqrt(2.)*M2*sintheta*v()
            +(6.*mH3*mH3 -(3.*sqrt(2.)*M1*v())/sintheta)*costheta*costheta
            +mHh*mHh+mHl*mHl
            -(mHh*mHh-mHl*mHl)*sign(mHh-mHl)*(cos(alpha)*cos(alpha)-sin(alpha)*sin(alpha)))/
           (6.*v()*v()*sintheta*sintheta);
}

double GeorgiMachacek::computelambda5() const
{
    return (2.*mH3*mH3-sqrt(2.)*M1*v()/sintheta)/(v()*v());
}

double GeorgiMachacek::computemu2sq() const
{
    return (3.*sqrt(2.)*M1*v()*sintheta
            -4.*(mHh*mHh+mHl*mHl)
            -2.*(mHh*mHh-mHl*mHl)*sign(mHh-mHl)*(2.*cos(alpha)*cos(alpha) -2.*sin(alpha)*sin(alpha)
                                                 +sqrt(6.)*cos(alpha)*sin(alpha)*sintheta/costheta))/16.;
}

double GeorgiMachacek::computemu3sq() const
{
    return (M1*v()*costheta*costheta +3.*M2*v()*sintheta*sintheta
            -((mHh*mHh+mHl*mHl)*sintheta)/sqrt(2.)
            -((mHh*mHh-mHl*mHl)*sign(mHh-mHl)*(3.*sqrt(2.)*(sin(alpha)*sin(alpha)-cos(alpha)*cos(alpha))*sintheta
                                               +8.*sqrt(3.)*cos(alpha)*sin(alpha)*costheta))/6.)/
           (2.*sqrt(2.)*sintheta);
}

void GeorgiMachacek::setParameter(const std::string name, const double& value){    
    if(name.compare("tantheta") == 0) {
        tantheta = value;
//        if(tantheta > 0.){
        logttheta = log(tantheta);
        sintheta = tantheta / sqrt(1. + tantheta*tantheta);
        costheta = 1. / sqrt(1. + tantheta*tantheta);}
//        else {
//            throw std::runtime_error("error in THDM::SetParameter, tantheta < 0!"); 
//          }
//        } 
    else if(name.compare("alpha") == 0) {
        alpha = value;
        sin_tma = sintheta*cos(alpha)-costheta*sin(alpha);
    }
    else if(name.compare("mHh") == 0)
        mHh = value;
    else if(name.compare("mH3") == 0)
        mH3 = value;
    else if(name.compare("mH5") == 0)
        mH5 = value;
    else if(name.compare("M1") == 0)
        M1 = value;
    else if(name.compare("M2") == 0)
        M2 = value;
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

double GeorgiMachacek::sign(const double x) const {
    return ( (x > 0) ? 1 : ((x < 0) ? -1 : 0) );
}
