/* 
 * Copyright (C) 2013 SusyFit Collaboration
 * All rights reserved.
 *
 * For the licensing terms see doc/COPYING.
 */

#include <stdexcept>
#include "NPEffective.h"


const std::string NPEffective::NPEffectiveVars[NNPEffectiveVars] 
= {"cWB_NP", "cH_NP", "cLL_NP", "cHLp_NP", "cHQp_NP", 
   "cHL_NP", "cHQ_NP", "cHE_NP", "cHU_NP", "cHD_NP", "Lambda_NP"};


NPEffective::NPEffective() 
: StandardModel() 
{
}


bool NPEffective::Update(const std::map<std::string,double>& DPars) 
{
    for (std::map<std::string, double>::const_iterator it = DPars.begin(); it != DPars.end(); it++)
        SetParameter(it->first, it->second);
    if(!StandardModel::Update(DPars)) return (false);

    return (true);
}


bool NPEffective::Init(const std::map<std::string, double>& DPars) 
{
    Update(DPars);
    return(CheckParameters(DPars)); 
}


bool NPEffective::CheckParameters(const std::map<std::string, double>& DPars) 
{
    for (int i = 0; i < NNPEffectiveVars; i++) {
        if (DPars.find(NPEffectiveVars[i]) == DPars.end()) {
            std::cout << "ERROR: Missing mandatory NPEffective parameter" 
                      << NPEffectiveVars[i] << std::endl;
            return false;
        }
    }
    return(StandardModel::CheckParameters(DPars));
}

    
void NPEffective::SetParameter(const std::string name, const double& value) 
{
    if (name.compare("cWB_NP") == 0)
        cWB = value;
    else if (name.compare("cH_NP") == 0)
        cH = value;
    else if (name.compare("cLL_NP") == 0)
        cLL = value;
    else if (name.compare("cHLp_NP") == 0)
        cHLp = value;
    else if (name.compare("cHQp_NP") == 0)
        cHQp = value;
    else if (name.compare("cHL_NP") == 0)
        cHL = value;
    else if (name.compare("cHQ_NP") == 0)
        cHQ = value;
    else if (name.compare("cHE_NP") == 0)
        cHE = value;
    else if (name.compare("cHU_NP") == 0)
        cHU = value;
    else if (name.compare("cHD_NP") == 0)
        cHD = value;
    else if (name.compare("Lambda_NP") == 0)
        LambdaNP = value;
    else
        StandardModel::SetParameter(name, value);       
}


bool NPEffective::InitializeModel() 
{
    SetModelInitialized(StandardModel::InitializeModel());
    return (IsModelInitialized());
}


void NPEffective::SetEWSMflags(EWSM& myEWSM) 
{
    StandardModel::SetEWSMflags(myEWSM);
}


bool NPEffective::SetFlag(const std::string name, const bool& value) 
{
    bool res = false;
    
    if (name.compare("EWABC") == 0) 
        throw std::runtime_error("ERROR: Flag EWABC is not applicable to NPEffective"); 
    else if (name.compare("EWABC2") == 0)
        throw std::runtime_error("ERROR: Flag EWABC2 is not applicable to NPEffective"); 
    else if (name.compare("epsilon1SM") == 0) 
        throw std::runtime_error("ERROR: Flag epsilon1SM is not applicable to NPEffective"); 
    else if (name.compare("epsilon2SM") == 0) 
        throw std::runtime_error("ERROR: Flag epsilon2SM is not applicable to NPEffective"); 
    else if (name.compare("epsilon3SM") == 0) 
        throw std::runtime_error("ERROR: Flag epsilon3SM is not applicable to NPEffective"); 
    else if (name.compare("epsilonbSM") == 0) 
        throw std::runtime_error("ERROR: Flag epsilonbSM is not applicable to NPEffective"); 
    else if (name.compare("NPZbbbarLinearize") == 0) 
        throw std::runtime_error("ERROR: Flag NPZbbbarLinearize is not applicable to NPEffective"); 
    else
        res = StandardModel::SetFlag(name,value);

    return(res);
}


////////////////////////////////////////////////////////////////////////     

double NPEffective::obliqueS() const 
{
    double sW_SM = sqrt(StandardModel::sW2()); /* This has to be the SM value. */
    double cW_SM = sqrt(StandardModel::cW2()); /* This has to be the SM value. */
    double ratio = v()*v()/LambdaNP/LambdaNP;

    return ( 4.0*sW_SM*cW_SM*cWB/alphaMz()*ratio );
}


double NPEffective::obliqueT() const
{
    double ratio = v()*v()/LambdaNP/LambdaNP;

    return ( - cH/2.0/alphaMz()*ratio + (cLL - 2.0*cHLp)/alphaMz()*ratio );    
}


double NPEffective::obliqueU() const
{
    double sW2_SM = StandardModel::sW2(); /* This has to be the SM value. */
    double ratio = v()*v()/LambdaNP/LambdaNP;

    return ( - 4.0*sW2_SM*(cLL - 2.0*cHLp)/alphaMz()*ratio );
}


////////////////////////////////////////////////////////////////////////     

double NPEffective::deltaGLl(StandardModel::lepton l) const
{
    double ratio = v()*v()/LambdaNP/LambdaNP;
    switch (l) {
        case StandardModel::NEUTRINO_1:
        case StandardModel::NEUTRINO_2:
        case StandardModel::NEUTRINO_3:
            return ( (cHLp - cHL)/2.0*ratio );
        case StandardModel::ELECTRON:
        case StandardModel::MU:
        case StandardModel::TAU:
            return ( - (cHLp + cHL)/2.0*ratio );
        default:
            throw std::runtime_error("Error in NPEffective::deltaGLl()");        
    }   
}


double NPEffective::deltaGLq(StandardModel::quark q) const
{
    double ratio = v()*v()/LambdaNP/LambdaNP;
    switch (q) {
        case StandardModel::UP:
        case StandardModel::CHARM:
            return ( (cHQp - cHQ)/2.0*ratio );
        case StandardModel::TOP:
            return 0.0;
        case StandardModel::DOWN:
        case StandardModel::STRANGE:
        case StandardModel::BOTTOM:
            return ( - (cHQp + cHQ)/2.0*ratio );
        default:
            throw std::runtime_error("Error in NPEffective::deltaGLq()");        
    }
}


double NPEffective::deltaGRl(StandardModel::lepton l) const
{
    double ratio = v()*v()/LambdaNP/LambdaNP;
    switch (l) {
        case StandardModel::NEUTRINO_1:
        case StandardModel::NEUTRINO_2:
        case StandardModel::NEUTRINO_3:
            return 0.0;
        case StandardModel::ELECTRON:
        case StandardModel::MU:
        case StandardModel::TAU:
            return ( - cHE/2.0*ratio );
        default:
            throw std::runtime_error("Error in NPEffective::deltaGRl()");        
    } 
}


double NPEffective::deltaGRq(StandardModel::quark q) const
{
    double ratio = v()*v()/LambdaNP/LambdaNP;
    switch (q) {
        case StandardModel::UP:
        case StandardModel::CHARM:
            return ( - cHU/2.0*ratio );
        case StandardModel::TOP:
            return 0.0;
        case StandardModel::DOWN:
        case StandardModel::STRANGE:
        case StandardModel::BOTTOM:
            return ( - cHD/2.0*ratio );
        default:
            throw std::runtime_error("Error in NPEffective::deltaGRq()");        
    }   
}


double NPEffective::deltaGVl(StandardModel::lepton l) const
{
    return ( deltaGLl(l) + deltaGRl(l) );
}


double NPEffective::deltaGVq(StandardModel::quark q) const
{
    return ( deltaGLq(q) + deltaGRq(q) );
}


double NPEffective::deltaGAl(StandardModel::lepton l) const
{
    return ( deltaGLl(l) - deltaGRl(l) );
}


double NPEffective::deltaGAq(StandardModel::quark q) const
{
    return ( deltaGLq(q) - deltaGRq(q) );
}


