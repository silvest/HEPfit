/* 
 * Copyright (C) 2015 SusyFit Collaboration
 * 
 *
 * For the licensing terms see doc/COPYING.
 */

#include <StandardModelMatching.h>
#include "GeorgiMachacek.h"
#include "GMcache.h"

std::string GeorgiMachacek::GMvars[NGMvars] = {"vDelta", "alpha", "mHh", "mA", "mH5", "Mu1", "Mu2", "Q_GM"};

GeorgiMachacek::GeorgiMachacek() : NPbase(), GMM(*this) {
    SMM.setObj((StandardModelMatching&) GMM.getObj());
    ModelParamMap.insert(std::pair<std::string, boost::reference_wrapper<const double> >("vDelta", boost::cref(vDelta)));
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
    if(name.compare("vDelta") == 0) {
        if(vDelta >= 0.) {
            vDelta = value;
        }
        else {
            throw std::runtime_error("error in GeorgiMachacek::SetParameter, vDelta < 0!"); 
          }
        }
    else if(name.compare("alpha") == 0) {
        alpha = value;
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
        NPbase::setParameter(name,value);
}

bool GeorgiMachacek::CheckParameters(const std::map<std::string, double>& DPars) {
    for (int i = 0; i < NGMvars; i++) {
        if (DPars.find(GMvars[i]) == DPars.end()) {
            std::cout << "missing mandatory GeorgiMachacek parameter " << GMvars[i] << std::endl;
            return false;
        }
    }
    return(NPbase::CheckParameters(DPars));
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

double GeorgiMachacek::muggH(const double sqrt_s) const
{
    return getMyGMCache()->rh_gg;
}

double GeorgiMachacek::muVBF(const double sqrt_s) const
{
    return getMyGMCache()->rh_VV;
}

double GeorgiMachacek::mueeWBF(const double sqrt_s) const
{
    return getMyGMCache()->rh_VV;
}

double GeorgiMachacek::muWH(const double sqrt_s) const
{
    return getMyGMCache()->rh_VV;
}

double GeorgiMachacek::muZH(const double sqrt_s) const
{
    return getMyGMCache()->rh_VV;
}

double GeorgiMachacek::mueeZH(const double sqrt_s) const
{
    return getMyGMCache()->rh_VV;
}

double GeorgiMachacek::muVH(const double sqrt_s) const
{
    return getMyGMCache()->rh_VV;
}

double GeorgiMachacek::muVBFpVH(const double sqrt_s) const
{
    return getMyGMCache()->rh_VV;
}

double GeorgiMachacek::muttH(const double sqrt_s) const
{
    return getMyGMCache()->rh_ff;
}

double GeorgiMachacek::computeGammaTotalRatio() const
{
    return getMyGMCache()->sumModBRs;
}

double GeorgiMachacek::GammaTotal() const
{
    return getMyGMCache()->Gamma_h;
}

double GeorgiMachacek::BrHggRatio() const
{
    return getMyGMCache()->rh_gg / computeGammaTotalRatio();
}

double GeorgiMachacek::BrHWWRatio() const
{
    return getMyGMCache()->rh_VV / computeGammaTotalRatio();
}

double GeorgiMachacek::BrHZZRatio() const
{
    return getMyGMCache()->rh_VV / computeGammaTotalRatio();
}

double GeorgiMachacek::BrHZgaRatio() const
{
    return getMyGMCache()->rh_Zga / computeGammaTotalRatio();
}

double GeorgiMachacek::BrHgagaRatio() const
{
    return getMyGMCache()->rh_gaga / computeGammaTotalRatio();
}

double GeorgiMachacek::BrHmumuRatio() const
{
    return getMyGMCache()->rh_ff / computeGammaTotalRatio();
}

double GeorgiMachacek::BrHtautauRatio() const
{
    return getMyGMCache()->rh_ff / computeGammaTotalRatio();
}

double GeorgiMachacek::BrHccRatio() const
{
    return getMyGMCache()->rh_ff / computeGammaTotalRatio();
}

double GeorgiMachacek::BrHbbRatio() const
{
    return getMyGMCache()->rh_ff / computeGammaTotalRatio();
}

double GeorgiMachacek::muggHgaga(const double sqrt_s) const
{
    return getMyGMCache()->rh_gg * getMyGMCache()->rh_gaga / computeGammaTotalRatio();
}

double GeorgiMachacek::muVBFHgaga(const double sqrt_s) const
{
    return getMyGMCache()->rh_VV * getMyGMCache()->rh_gaga / computeGammaTotalRatio();
}

double GeorgiMachacek::muVHgaga(const double sqrt_s) const
{
    return getMyGMCache()->rh_VV * getMyGMCache()->rh_gaga / computeGammaTotalRatio();
}

double GeorgiMachacek::muttHgaga(const double sqrt_s) const
{
    return getMyGMCache()->rh_ff * getMyGMCache()->rh_gaga / computeGammaTotalRatio();
}

double GeorgiMachacek::muggHZZ(const double sqrt_s) const
{
    return getMyGMCache()->rh_gg * getMyGMCache()->rh_VV / computeGammaTotalRatio();
}

double GeorgiMachacek::muVBFHZZ(const double sqrt_s) const
{
    return getMyGMCache()->rh_VV * getMyGMCache()->rh_VV / computeGammaTotalRatio();
}

double GeorgiMachacek::muVHZZ(const double sqrt_s) const
{
    return getMyGMCache()->rh_VV * getMyGMCache()->rh_VV / computeGammaTotalRatio();
}

double GeorgiMachacek::muttHZZ(const double sqrt_s) const
{
    return getMyGMCache()->rh_ff * getMyGMCache()->rh_VV / computeGammaTotalRatio();
}

double GeorgiMachacek::muggHWW(const double sqrt_s) const
{
    return getMyGMCache()->rh_gg * getMyGMCache()->rh_VV / computeGammaTotalRatio();
}

double GeorgiMachacek::muVBFHWW(const double sqrt_s) const
{
    return getMyGMCache()->rh_VV * getMyGMCache()->rh_VV / computeGammaTotalRatio();
}

double GeorgiMachacek::muVHWW(const double sqrt_s) const
{
    return getMyGMCache()->rh_VV * getMyGMCache()->rh_VV / computeGammaTotalRatio();
}

double GeorgiMachacek::muttHWW(const double sqrt_s) const
{
    return getMyGMCache()->rh_ff * getMyGMCache()->rh_VV / computeGammaTotalRatio();
}

double GeorgiMachacek::muggHtautau(const double sqrt_s) const
{
    return getMyGMCache()->rh_gg * getMyGMCache()->rh_ff / computeGammaTotalRatio();
}

double GeorgiMachacek::muVBFHtautau(const double sqrt_s) const
{
    return getMyGMCache()->rh_VV * getMyGMCache()->rh_ff / computeGammaTotalRatio();
}

double GeorgiMachacek::muVHtautau(const double sqrt_s) const
{
    return getMyGMCache()->rh_VV * getMyGMCache()->rh_ff / computeGammaTotalRatio();
}

double GeorgiMachacek::muttHtautau(const double sqrt_s) const
{
    return getMyGMCache()->rh_ff * getMyGMCache()->rh_ff / computeGammaTotalRatio();
}

double GeorgiMachacek::muggHbb(const double sqrt_s) const
{
    return getMyGMCache()->rh_gg * getMyGMCache()->rh_ff / computeGammaTotalRatio();
}

double GeorgiMachacek::muVBFHbb(const double sqrt_s) const
{
    return getMyGMCache()->rh_VV * getMyGMCache()->rh_ff / computeGammaTotalRatio();
}

double GeorgiMachacek::muVHbb(const double sqrt_s) const
{
    return getMyGMCache()->rh_VV * getMyGMCache()->rh_ff / computeGammaTotalRatio();
}

double GeorgiMachacek::muttHbb(const double sqrt_s) const
{
    return getMyGMCache()->rh_ff * getMyGMCache()->rh_ff / computeGammaTotalRatio();
}

double GeorgiMachacek::muppHmumu(const double sqrt_s) const
{
    if(sqrt_s==8)
    {
        return (0.872 * getMyGMCache()->rh_gg + 0.122 * getMyGMCache()->rh_VV + 0.006 * getMyGMCache()->rh_ff) * getMyGMCache()->rh_ff / computeGammaTotalRatio();
    }
    else if(sqrt_s==13)
    {
        return (0.871 * getMyGMCache()->rh_gg + 0.119 * getMyGMCache()->rh_VV + 0.010 * getMyGMCache()->rh_ff) * getMyGMCache()->rh_ff / computeGammaTotalRatio();
    }
    else
    {
        throw std::runtime_error("The observable muppHmumu is only defined for 8 or 13 TeV.");
    }
}

double GeorgiMachacek::muppHZga(const double sqrt_s) const
{
    if(sqrt_s==8)
    {
        return (0.872 * getMyGMCache()->rh_gg + 0.122 * getMyGMCache()->rh_VV + 0.006 * getMyGMCache()->rh_ff) * getMyGMCache()->rh_Zga / computeGammaTotalRatio();
    }
    else if(sqrt_s==13)
    {
        return (0.871 * getMyGMCache()->rh_gg + 0.119 * getMyGMCache()->rh_VV + 0.010 * getMyGMCache()->rh_ff) * getMyGMCache()->rh_Zga / computeGammaTotalRatio();
    }
    else
    {
        throw std::runtime_error("The observable muppHZga is only defined for 8 or 13 TeV.");
    }
}


double GeorgiMachacek::Mw() const{
    double MZ = StandardModel::Mz;
    return ( MZ / sqrt(2.0) * sqrt(1.0 + sqrt(1.0 - 4.0 * M_PI * StandardModel::ale / (sqrt(2.0) * StandardModel::GF * MZ* MZ))));
}
