/* 
 * Copyright (C) 2013 SusyFit Collaboration
 * All rights reserved.
 *
 * For the licensing terms see doc/COPYING.
 */

#include <stdexcept>
#include "NPHiggsST.h"


const std::string NPHiggsST::NPHIGGSSTvars[NNPHIGGSSTvars] 
= {"a", "b", "c_u", "c_d", "c_e", "d_3", "d_4"};


NPHiggsST::NPHiggsST()
: NPZbbbar()
{
}


bool NPHiggsST::Update(const std::map<std::string,double>& DPars) 
{
    for (std::map<std::string, double>::const_iterator it = DPars.begin(); it != DPars.end(); it++)
        SetParameter(it->first, it->second);
    if(!NPZbbbar::Update(DPars)) return (false);

    return (true);
}


bool NPHiggsST::Init(const std::map<std::string, double>& DPars) 
{
    Update(DPars);
    return(CheckParameters(DPars)); 
}


bool NPHiggsST::CheckParameters(const std::map<std::string, double>& DPars) 
{
    for (int i = 0; i < NNPHIGGSSTvars; i++) {
        if (DPars.find(NPHIGGSSTvars[i]) == DPars.end()) {
            std::cout << "ERROR: Missing mandatory NPHiggsST parameter" 
                      << NPHIGGSSTvars[i] << std::endl;
            return false;
        }
    }
    return(NPZbbbar::CheckParameters(DPars));
}

    
void NPHiggsST::SetParameter(const std::string name, const double& value) 
{
    if (name.compare("a") == 0)
        a = value;
    else if (name.compare("b") == 0)
        b = value;
    else if (name.compare("c_u") == 0)
        c_u = value;
    else if (name.compare("c_d") == 0)
        c_d = value;
    else if (name.compare("c_e") == 0)
        c_e = value;
    else if (name.compare("d_3") == 0)
        d_3 = value;
    else if (name.compare("d_4") == 0)
        d_4 = value;
    else
        NPZbbbar::SetParameter(name, value);       
}


bool NPHiggsST::InitializeModel() 
{
    SetModelInitialized(NPZbbbar::InitializeModel());
    return (IsModelInitialized());
}


void NPHiggsST::SetEWSMflags(EWSM& myEWSM) 
{
    StandardModel::SetEWSMflags(myEWSM);
}


bool NPHiggsST::SetFlag(const std::string name, const bool& value) 
{
    bool res = false;
    if (name.compare("epsilon1SM") == 0) 
        throw std::runtime_error("ERROR: Flag epsilon1SM is not applicable to NPHiggsST"); 
    else if (name.compare("epsilon2SM") == 0) 
        throw std::runtime_error("ERROR: Flag epsilon2SM is not applicable to NPHiggsST"); 
    else if (name.compare("epsilon3SM") == 0) 
        throw std::runtime_error("ERROR: Flag epsilon3SM is not applicable to NPHiggsST"); 
    else if (name.compare("epsilonbSM") == 0) 
        throw std::runtime_error("ERROR: Flag epsilonbSM is not applicable to NPHiggsST"); 
    else
        res = NPZbbbar::SetFlag(name,value);

    return(res);
}


////////////////////////////////////////////////////////////////////////

double NPHiggsST::obliqueS() const
{
    double Lambda;
    if (fabs(1.0-a*a) < pow(10.0, -32.0) ) 
        Lambda = pow(10.0, 19.0);
    else
        Lambda = 4.0*M_PI*v()/sqrt(fabs(1.0 - a*a));

    return ( 1.0/12.0/M_PI*(1.0 - a*a)*log(Lambda*Lambda/mHl/mHl) );
}
        

double NPHiggsST::obliqueT() const
{
    double Lambda;
    double cW2_SM = StandardModel::cW2(); /* This has to be the SM value. */
    if (fabs(1.0-a*a) < pow(10.0, -32.0) ) 
        Lambda = pow(10.0, 19.0);
    else
        Lambda = 4.0*M_PI*v()/sqrt(fabs(1.0 - a*a));
    
    return ( - 3.0/16.0/M_PI/cW2_SM*(1.0 - a*a)*log(Lambda*Lambda/mHl/mHl) );
}
    

double NPHiggsST::obliqueU() const
{
    return 0.0;
}


////////////////////////////////////////////////////////////////////////

double NPHiggsST::epsilon1() const
{
    throw std::runtime_error("ERROR: NPHiggsST::epsilon1() is not implemented");
}


double NPHiggsST::epsilon2() const
{
    throw std::runtime_error("ERROR: NPHiggsST::epsilon2() is not implemented");
}


double NPHiggsST::epsilon3() const
{
    throw std::runtime_error("ERROR: NPHiggsST::epsilon3() is not implemented");
}


double NPHiggsST::epsilonb() const
{
    throw std::runtime_error("ERROR: NPHiggsST::epsilonb() is not implemented");
}

    
