/* 
 * Copyright (C) 2013 SusyFit Collaboration
 * All rights reserved.
 *
 * For the licensing terms see doc/COPYING.
 */

#include <stdexcept>
#include "NPHiggsST.h"


const std::string NPHiggsST::NPHIGGSSTvars[NNPHIGGSSTvars] 
= {"a", "b", "c_u", "c_d", "c_e", "d_3", "d_4", "LambdaNP"};


NPHiggsST::NPHiggsST()
: NPZbbbar(), LambdaNP_in(0.0)
{
}


bool NPHiggsST::Update(const std::map<std::string,double>& DPars) 
{
    for (std::map<std::string, double>::const_iterator it = DPars.begin(); it != DPars.end(); it++)
        parseParameters(it->first, it->second);
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

    
void NPHiggsST::parseParameters(const std::string name, const double& value) 
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
        NPZbbbar::parseParameters(name, value);       
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
    if (LambdaNP_in != 0.0)
        Lambda = LambdaNP_in;
    else if (fabs(1.0-a*a) < pow(10.0, -32.0))
        Lambda = pow(10.0, 19.0);
    else
        Lambda = 4.0*M_PI*v()/sqrt(fabs(1.0 - a*a));

    return ( 1.0/12.0/M_PI*(1.0 - a*a)*log(Lambda*Lambda/mHl/mHl) );
}
        

double NPHiggsST::obliqueT() const
{
    double Lambda;
    double cW2_SM = StandardModel::cW2(); /* This has to be the SM value. */
    if (LambdaNP_in != 0.0)
        Lambda = LambdaNP_in;
    else if (fabs(1.0-a*a) < pow(10.0, -32.0))
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


////////////////////////////////////////////////////////////////////////

double NPHiggsST::Mw() const
{
    double myMw = StandardModel::Mw();

    if (!IsFlagNotLinearizedNP() ) {
        double alpha = StandardModel::alphaMz();
        double c2 = StandardModel::cW2();
        double s2 = StandardModel::sW2();

        myMw *= 1.0 - alpha/4.0/(c2-s2)
                *( obliqueS() - 2.0*c2*obliqueT() - (c2-s2)*obliqueU()/2.0/s2 );
    } else
        if (obliqueS()!=0.0 || obliqueT()!=0.0 || obliqueU()!=0.0)
            throw std::runtime_error("NPHiggsST::Mw(): The oblique corrections STU cannot be used with flag NotLinearizedNP=1");

    return myMw;
}


double NPHiggsST::cW2() const
{
    return ( Mw()*Mw()/Mz/Mz );
}


double NPHiggsST::sW2() const
{
    return ( 1.0 - cW2() );
}


double NPHiggsST::GammaW() const
{
    double Gamma_W = StandardModel::GammaW();

    if (!IsFlagNotLinearizedNP() ) {
        double alpha = StandardModel::alphaMz();
        double c2 = StandardModel::cW2();
        double s2 = StandardModel::sW2();

        Gamma_W *= 1.0 - 3.0*alpha/4.0/(c2-s2)
                   *( obliqueS() - 2.0*c2*obliqueT()
                      - (c2-s2)*obliqueU()/2.0/s2 );
        } else
            if (obliqueS()!=0.0 || obliqueT()!=0.0 || obliqueU()!=0.0)
                throw std::runtime_error("NPHiggsST::GammaW(): The oblique corrections STU cannot be used with flag NotLinearizedNP=1");

    return Gamma_W;
}




    
