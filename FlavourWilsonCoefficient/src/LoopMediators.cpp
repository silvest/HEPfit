/* 
 * Copyright (C) 2018 HEPfit Collaboration
 *
 *
 * For the licensing terms see doc/COPYING.
 */

#include "LoopMediators.h"

const std::string LoopMediators::LoopMediatorsvars[NLoopMediatorsvars] = {"GammaL", "GammaR", "Gammamu", "mphi", "yQ", "yL", "WCscale"};

LoopMediators::LoopMediators() : StandardModel(), LoopMediatorsM(*this) {   

    SMM.setObj((StandardModelMatching&) LoopMediatorsM.getObj());
    ModelParamMap.insert(std::pair<std::string, boost::reference_wrapper<const double> >("GammaL", boost::cref(GammaL)));
    ModelParamMap.insert(std::pair<std::string, boost::reference_wrapper<const double> >("GammaR", boost::cref(GammaR)));
    ModelParamMap.insert(std::pair<std::string, boost::reference_wrapper<const double> >("Gammamu", boost::cref(Gammamu)));
    ModelParamMap.insert(std::pair<std::string, boost::reference_wrapper<const double> >("mphi", boost::cref(mphi)));
    ModelParamMap.insert(std::pair<std::string, boost::reference_wrapper<const double> >("yQ", boost::cref(yQ)));
    ModelParamMap.insert(std::pair<std::string, boost::reference_wrapper<const double> >("yL", boost::cref(yL)));
    
    ModelParamMap.insert(std::pair<std::string, boost::reference_wrapper<const double> >("WCscale", boost::cref(WCscale)));
}

LoopMediators::~LoopMediators(){
    if (IsModelInitialized()) {
        }
}

///////////////////////////////////////////////////////////////////////////
// Initialization

bool LoopMediators::InitializeModel()
{
    setModelInitialized(StandardModel::InitializeModel());
    return(true);
}
    
bool LoopMediators::Init(const std::map<std::string, double>& DPars) {
    return(StandardModel::Init(DPars));
}

bool LoopMediators::PreUpdate()
{    
    if(!StandardModel::PreUpdate()) return (false);

    return (true);
}

bool LoopMediators::Update(const std::map<std::string, double>& DPars) {
    
    if(!PreUpdate()) return (false);
    
    UpdateError = false;

    for (std::map<std::string, double>::const_iterator it = DPars.begin(); it != DPars.end(); it++)
        setParameter(it->first, it->second);

    if (UpdateError) return (false);

    if(!PostUpdate()) return (false);

    return (true);
}

bool LoopMediators::PostUpdate()
{
    if(!StandardModel::PostUpdate()) return (false);

    double GammaL2 = GammaL*GammaL;
    double GammaR2 = GammaR*GammaR;
    double GammaLR = GammaL*GammaR;
    double Gammamu2 = Gammamu * Gammamu;
    double mphi2 = mphi * mphi;
    
    double M_PI2 = M_PI*M_PI;
    gslpp::complex Norm = sqrt(2.) / (4. * GF * getCKM().computelamt_s()) / (32. * M_PI * ale);

    C1 = 1. / (128. * M_PI2 * mphi2) * GammaL2 * F9(yQ,yQ);
    C2 = 0.;
    C3 = 0.;
    C4 = 0.;
    C5 = - 1. / (32. * M_PI2 * mphi2) * GammaLR * F9(yQ,yQ);
    
    C1p = 1. / (128. * M_PI2 * mphi2) * GammaR2 * F9(yQ,yQ);
    C2p = 0.;
    C3p = 0.;
    
    C7 = 0.;
    C8 = 0.;
    C9 = - Norm * GammaL*Gammamu2 / mphi2 * F9(yQ,yL);
    C10 =  Norm * GammaL*Gammamu2 / mphi2 * F9(yQ,yL);
    CS = 0.;
    CP = 0.;
    
    C7p = 0.;
    C8p = 0.;
    C9p =  - Norm * GammaR*Gammamu2 / mphi2 * F9(yQ,yL);
    C10p = - Norm * GammaR*Gammamu2 / mphi2 * F9(yQ,yL);
    CSp = 0.;
    CPp = 0.;

    /* Necessary for updating StandardModel parameters in StandardModelMatching,
     * and LoopMediators and LoopMediators-derived parameters in LoopMediatorsMatching */
    LoopMediatorsM.getObj().updateLoopMediatorsParameters();

    return (true);
}

void LoopMediators::setParameter(const std::string name, const double& value){    
    if(name.compare("GammaL") == 0) 
        GammaL = value;  
    else if(name.compare("GammaR") == 0) 
        GammaR = value;  
    else if(name.compare("Gammamu") == 0) 
        Gammamu = value;  
    else if(name.compare("mphi") == 0) 
        mphi = value;
    else if(name.compare("yQ") == 0) 
        yQ = value;
    else if(name.compare("yL") == 0) 
        yL = value;
    else if(name.compare("WCscale") == 0)
        WCscale = value;
    else
        StandardModel::setParameter(name,value);
}

bool LoopMediators::CheckParameters(const std::map<std::string, double>& DPars) {
    for (int i = 0; i < NLoopMediatorsvars; i++) {
        if (DPars.find(LoopMediatorsvars[i]) == DPars.end()) {
            std::cout << "missing mandatory LoopMediators parameter " << LoopMediatorsvars[i] << std::endl;
            return false;
        }
    }
    return(StandardModel::CheckParameters(DPars));
}

///////////////////////////////////////////////////////////////////////////
// Flags

bool LoopMediators::setFlag(const std::string name, const bool value)
{
    bool res = false;
    
    res = StandardModel::setFlag(name,value);

    return(res);
}



double LoopMediators::F9(double x, double y)
{
    double ym1 = y - 1.;
    double xm1 = x - 1.;

    if (x == 1. && y == 1.)
        return 1./3.;
    else if (x == 1.)
        return (-3.*y*y + 2.*y*y*log(y) + 4.*y - 1.)/2./ym1/ym1/ym1;
    else if (y == 1.)
        return (-3.*x*x + 2.*x*x*log(x) + 4.*x - 1.)/2./xm1/xm1/xm1;
    else if (x == y)
        return (y*y - 2.*y*log(y) - 1.)/ym1/ym1/ym1;
    else
        return 1./xm1/ym1 + x*x*log(x)/xm1/xm1/(x-y) + y*y*log(y)/ym1/ym1/(y-x);
}
