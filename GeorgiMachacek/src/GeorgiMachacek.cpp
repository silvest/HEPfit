/* 
 * Copyright (C) 2015 HEPfit Collaboration
 * 
 *
 * For the licensing terms see doc/COPYING.
 */

#include <StandardModelMatching.h>
#include <algorithm>
#include "GeorgiMachacek.h"
#include "GMcache.h"

std::string GeorgiMachacek::GMvars[NGMvars] = {"vDelta", "alpha", "mHh", "mA", "mH5", "Mu1", "Mu2", "Q_GM"};

GeorgiMachacek::GeorgiMachacek() : NPbase(), GMM(*this) {
    SMM.setObj((StandardModelMatching&) GMM.getObj());
    ModelParamMap.insert(std::make_pair("vDelta", std::cref(vDelta)));
    ModelParamMap.insert(std::make_pair("alpha", std::cref(alpha)));
    ModelParamMap.insert(std::make_pair("mHh", std::cref(mHh)));
    ModelParamMap.insert(std::make_pair("mA", std::cref(mA)));
    ModelParamMap.insert(std::make_pair("mH5", std::cref(mH5)));
    ModelParamMap.insert(std::make_pair("Mu1", std::cref(Mu1)));
    ModelParamMap.insert(std::make_pair("Mu2", std::cref(Mu2)));
    ModelParamMap.insert(std::make_pair("Q_GM", std::cref(Q_GM)));
    flag_use_sq_masses=false;
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
            std::cout << "ERROR: missing mandatory GeorgiMachacek parameter " << GMvars[i] << std::endl;
            raiseMissingModelParameterCount();
            addMissingModelParameter(GMvars[i]);
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
            
            ModelParamMap.insert(std::make_pair("mHhsq", std::cref(mHhsq)));
            ModelParamMap.insert(std::make_pair("mAsq", std::cref(mAsq)));
            ModelParamMap.insert(std::make_pair("mH5sq", std::cref(mH5sq)));
        }
    }
    else
        res = StandardModel::setFlag(name,value);

    return(res);
}

const double GeorgiMachacek::muggH(const double sqrt_s) const
{
    return getMyGMCache()->rh_gg;
}

const double GeorgiMachacek::muVBF(const double sqrt_s) const
{
    return getMyGMCache()->rh_VV;
}

const double GeorgiMachacek::mueeWBF(const double sqrt_s, const double Pol_em, const double Pol_ep) const
{
    return getMyGMCache()->rh_VV;
}

const double GeorgiMachacek::muWH(const double sqrt_s) const
{
    return getMyGMCache()->rh_VV;
}

const double GeorgiMachacek::muZH(const double sqrt_s) const
{
    return getMyGMCache()->rh_VV;
}

const double GeorgiMachacek::mueeZH(const double sqrt_s, const double Pol_em, const double Pol_ep) const
{
    return getMyGMCache()->rh_VV;
}

const double GeorgiMachacek::muVH(const double sqrt_s) const
{
    return getMyGMCache()->rh_VV;
}

const double GeorgiMachacek::muVBFpVH(const double sqrt_s) const
{
    return getMyGMCache()->rh_VV;
}

const double GeorgiMachacek::muttH(const double sqrt_s) const
{
    return getMyGMCache()->rh_ff;
}

const double GeorgiMachacek::computeGammaTotalRatio() const
{
    return getMyGMCache()->sumModBRs;
}

const double GeorgiMachacek::GammaTotal() const
{
    return getMyGMCache()->Gamma_h;
}

const double GeorgiMachacek::BrHggRatio() const
{
    return getMyGMCache()->rh_gg / computeGammaTotalRatio();
}

const double GeorgiMachacek::BrHWWRatio() const
{
    return getMyGMCache()->rh_VV / computeGammaTotalRatio();
}

const double GeorgiMachacek::BrHZZRatio() const
{
    return getMyGMCache()->rh_VV / computeGammaTotalRatio();
}

const double GeorgiMachacek::BrHZgaRatio() const
{
    return getMyGMCache()->rh_Zga / computeGammaTotalRatio();
}

const double GeorgiMachacek::BrHgagaRatio() const
{
    return getMyGMCache()->rh_gaga / computeGammaTotalRatio();
}

const double GeorgiMachacek::BrHmumuRatio() const
{
    return getMyGMCache()->rh_ff / computeGammaTotalRatio();
}

const double GeorgiMachacek::BrHtautauRatio() const
{
    return getMyGMCache()->rh_ff / computeGammaTotalRatio();
}

const double GeorgiMachacek::BrHccRatio() const
{
    return getMyGMCache()->rh_ff / computeGammaTotalRatio();
}

const double GeorgiMachacek::BrHbbRatio() const
{
    return getMyGMCache()->rh_ff / computeGammaTotalRatio();
}

const double GeorgiMachacek::muggHgaga(const double sqrt_s) const
{
    return getMyGMCache()->rh_gg * getMyGMCache()->rh_gaga / computeGammaTotalRatio();
}

const double GeorgiMachacek::muVBFHgaga(const double sqrt_s) const
{
    return getMyGMCache()->rh_VV * getMyGMCache()->rh_gaga / computeGammaTotalRatio();
}

const double GeorgiMachacek::muVHgaga(const double sqrt_s) const
{
    return getMyGMCache()->rh_VV * getMyGMCache()->rh_gaga / computeGammaTotalRatio();
}

const double GeorgiMachacek::muttHgaga(const double sqrt_s) const
{
    return getMyGMCache()->rh_ff * getMyGMCache()->rh_gaga / computeGammaTotalRatio();
}

const double GeorgiMachacek::muggHZZ(const double sqrt_s) const
{
    return getMyGMCache()->rh_gg * getMyGMCache()->rh_VV / computeGammaTotalRatio();
}

const double GeorgiMachacek::muVBFHZZ(const double sqrt_s) const
{
    return getMyGMCache()->rh_VV * getMyGMCache()->rh_VV / computeGammaTotalRatio();
}

const double GeorgiMachacek::muVHZZ(const double sqrt_s) const
{
    return getMyGMCache()->rh_VV * getMyGMCache()->rh_VV / computeGammaTotalRatio();
}

const double GeorgiMachacek::muttHZZ(const double sqrt_s) const
{
    return getMyGMCache()->rh_ff * getMyGMCache()->rh_VV / computeGammaTotalRatio();
}

const double GeorgiMachacek::muggHWW(const double sqrt_s) const
{
    return getMyGMCache()->rh_gg * getMyGMCache()->rh_VV / computeGammaTotalRatio();
}

const double GeorgiMachacek::muVBFHWW(const double sqrt_s) const
{
    return getMyGMCache()->rh_VV * getMyGMCache()->rh_VV / computeGammaTotalRatio();
}

const double GeorgiMachacek::muVHWW(const double sqrt_s) const
{
    return getMyGMCache()->rh_VV * getMyGMCache()->rh_VV / computeGammaTotalRatio();
}

const double GeorgiMachacek::muttHWW(const double sqrt_s) const
{
    return getMyGMCache()->rh_ff * getMyGMCache()->rh_VV / computeGammaTotalRatio();
}

const double GeorgiMachacek::muggHtautau(const double sqrt_s) const
{
    return getMyGMCache()->rh_gg * getMyGMCache()->rh_ff / computeGammaTotalRatio();
}

const double GeorgiMachacek::muVBFHtautau(const double sqrt_s) const
{
    return getMyGMCache()->rh_VV * getMyGMCache()->rh_ff / computeGammaTotalRatio();
}

const double GeorgiMachacek::muVHtautau(const double sqrt_s) const
{
    return getMyGMCache()->rh_VV * getMyGMCache()->rh_ff / computeGammaTotalRatio();
}

const double GeorgiMachacek::muttHtautau(const double sqrt_s) const
{
    return getMyGMCache()->rh_ff * getMyGMCache()->rh_ff / computeGammaTotalRatio();
}

const double GeorgiMachacek::muggHbb(const double sqrt_s) const
{
    return getMyGMCache()->rh_gg * getMyGMCache()->rh_ff / computeGammaTotalRatio();
}

const double GeorgiMachacek::muVBFHbb(const double sqrt_s) const
{
    return getMyGMCache()->rh_VV * getMyGMCache()->rh_ff / computeGammaTotalRatio();
}

const double GeorgiMachacek::muVHbb(const double sqrt_s) const
{
    return getMyGMCache()->rh_VV * getMyGMCache()->rh_ff / computeGammaTotalRatio();
}

const double GeorgiMachacek::muttHbb(const double sqrt_s) const
{
    return getMyGMCache()->rh_ff * getMyGMCache()->rh_ff / computeGammaTotalRatio();
}

const double GeorgiMachacek::muppHmumu(const double sqrt_s) const
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

const double GeorgiMachacek::muppHZga(const double sqrt_s) const
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


const double GeorgiMachacek::Mw() const{
    double MZ = StandardModel::Mz;
    return ( MZ / sqrt(2.0) * sqrt(1.0 + sqrt(1.0 - 4.0 * M_PI * StandardModel::ale / (sqrt(2.0) * StandardModel::GF * MZ* MZ))));
}
