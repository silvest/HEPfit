/* 
 * Copyright (C) 2013 SusyFit Collaboration
 * All rights reserved.
 *
 * For the licensing terms see doc/COPYING.
 */

#include <stdexcept>
#include "NPHiggsST.h"


bool NPHiggsST::Update(const std::map<std::string,double>& DPars) {
    for (std::map<std::string, double>::const_iterator it = DPars.begin(); it != DPars.end(); it++)
        SetParameter(it->first, it->second);
    if(!NPHiggs::Update(DPars)) return (false);

    return (true);
}


bool NPHiggsST::Init(const std::map<std::string, double>& DPars) {
    Update(DPars);
    return(CheckParameters(DPars)); 
}


bool NPHiggsST::CheckParameters(const std::map<std::string, double>& DPars) {
    return(NPHiggs::CheckParameters(DPars));
}

    
void NPHiggsST::SetParameter(const std::string name, const double& value) {
    NPHiggs::SetParameter(name, value);       
}


bool NPHiggsST::InitializeModel() 
{
    SetModelInitialized(NPHiggs::InitializeModel());
    return (IsModelInitialized());
}


void NPHiggsST::SetEWSMflags(EWSM& myEWSM) 
{
    StandardModel::SetEWSMflags(myEWSM);
}


bool NPHiggsST::SetFlag(const std::string name, const bool& value) {
    bool res = false;
    res = NPHiggs::SetFlag(name,value);
    return(res);
}


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
    if (fabs(1.0-a*a) < pow(10.0, -32.0) ) 
        Lambda = pow(10.0, 19.0);
    else
        Lambda = 4.0*M_PI*v()/sqrt(fabs(1.0 - a*a));
    
    return ( - 3.0/16.0/M_PI/c02()*(1.0 - a*a)*log(Lambda*Lambda/mHl/mHl) );
}
    

double NPHiggsST::obliqueU() const
{
    return 0.0;
}
    

double NPHiggsST::epsilon1() const
{
    return ( epsilon1_SM() + alphaMz()*obliqueT() ); 
}


double NPHiggsST::epsilon2() const
{
    return epsilon2_SM();
}

    
double NPHiggsST::epsilon3() const
{
    return ( epsilon3_SM() + alphaMz()/4.0/s02()*obliqueS() ); 
}

    
double NPHiggsST::epsilonb() const
{
    return epsilonb_SM();
}


double NPHiggsST::Mw() const
{
    return StandardModel::Mw();
}


double NPHiggsST::cW2() const
{
    return StandardModel::cW2();
}


double NPHiggsST::sW2() const
{
    return StandardModel::sW2();
}


complex NPHiggsST::rhoZ_l(const StandardModel::lepton l) const
{
    return NPZbbbar::rhoZ_l(l);
}

    
complex NPHiggsST::rhoZ_q(const StandardModel::quark q) const
{
    return NPZbbbar::rhoZ_q(q);
}


complex NPHiggsST::kappaZ_l(const StandardModel::lepton l) const
{
    return NPZbbbar::kappaZ_l(l);
}


complex NPHiggsST::kappaZ_q(const StandardModel::quark q) const
{
    return NPZbbbar::kappaZ_q(q);
}

    
complex NPHiggsST::gVl(const StandardModel::lepton l) const
{
    return NPZbbbar::gVl(l);
}


complex NPHiggsST::gVq(const StandardModel::quark q) const
{
    return NPZbbbar::gVq(q);
}


complex NPHiggsST::gAl(const StandardModel::lepton l) const
{
    return NPZbbbar::gAl(l);
}


complex NPHiggsST::gAq(const StandardModel::quark q) const
{
    return NPZbbbar::gAq(q);
}

    
double NPHiggsST::GammaW() const
{
    return StandardModel::GammaW();
}




