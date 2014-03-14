/*
 * Copyright (C) 2013 SusyFit Collaboration
 * All rights reserved.
 *
 * For the licensing terms see doc/COPYING.
 */

#include <stdexcept>
#include <EWSM.h>
#include "NPHiggs.h"


const std::string NPHiggs::NPHIGGSvars[NNPHIGGSvars]
= {"a", "b", "c_u", "c_d", "c_e", "d_3", "d_4", "LambdaNP"};


NPHiggs::NPHiggs()
: NPbase(), LambdaNP_in(0.0)
{
}


bool NPHiggs::InitializeModel()
{
    setModelInitialized(NPbase::InitializeModel());
    return (IsModelInitialized());
}


bool NPHiggs::Init(const std::map<std::string, double>& DPars)
{
    return(NPbase::Init(DPars));
}


bool NPHiggs::Update(const std::map<std::string,double>& DPars)
{
    for (std::map<std::string, double>::const_iterator it = DPars.begin(); it != DPars.end(); it++)
        setParameter(it->first, it->second);
    if(!NPbase::Update(DPars)) return (false);
    return (true);
}


void NPHiggs::setParameter(const std::string name, const double& value)
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
    else if (name.compare("LambdaNP") == 0)
        LambdaNP_in = value;
    else
        NPbase::setParameter(name, value);
}


bool NPHiggs::CheckParameters(const std::map<std::string, double>& DPars)
{
    for (int i = 0; i < NNPHIGGSvars; i++) {
        if (DPars.find(NPHIGGSvars[i]) == DPars.end()) {
            std::cout << "ERROR: Missing mandatory NPHiggs parameter "
                      << NPHIGGSvars[i] << std::endl;
            return false;
        }
    }
    return(NPbase::CheckParameters(DPars));
}


bool NPHiggs::setFlag(const std::string name, const bool value)
{
    bool res = false;
    res = NPbase::setFlag(name,value);
    return(res);
}


bool NPHiggs::CheckFlags() const
{
    return(NPbase::CheckFlags());
}


////////////////////////////////////////////////////////////////////////

double NPHiggs::obliqueS() const
{
    double Lambda;
    if (LambdaNP_in != 0.0)
        Lambda = LambdaNP_in;
    else if (fabs(1.0-a*a) < pow(10.0, -32.0))
        Lambda = pow(10.0, 19.0);
    else
        Lambda = 4.0*M_PI*v()/sqrt(fabs(1.0 - a*a));

    return ( 1.0/12.0/M_PI*(1.0 - a*a)*log(Lambda*Lambda/mHl/mHl) );
}


double NPHiggs::obliqueT() const
{
    double Lambda;
    double cW2_SM = myEWSM->cW2_SM(); /* This has to be the SM value. */
    if (LambdaNP_in != 0.0)
        Lambda = LambdaNP_in;
    else if (fabs(1.0-a*a) < pow(10.0, -32.0))
        Lambda = pow(10.0, 19.0);
    else
        Lambda = 4.0*M_PI*v()/sqrt(fabs(1.0 - a*a));

    return ( - 3.0/16.0/M_PI/cW2_SM*(1.0 - a*a)*log(Lambda*Lambda/mHl/mHl) );
}


double NPHiggs::obliqueU() const
{
    return 0.0;
}


