/*
 * Copyright (C) 2013 SusyFit Collaboration
 * All rights reserved.
 *
 * For the licensing terms see doc/COPYING.
 */

#include <EWSM.h>
#include "NPbase.h"


NPbase::NPbase()
: StandardModel()
{
}


bool NPbase::InitializeModel()
{
    setModelInitialized(StandardModel::InitializeModel());
    return (IsModelInitialized());
}


bool NPbase::Init(const std::map<std::string, double>& DPars)
{
    Update(DPars);
    return(CheckParameters(DPars));
}


bool NPbase::Update(const std::map<std::string,double>& DPars)
{
    for (std::map<std::string, double>::const_iterator it = DPars.begin(); it != DPars.end(); it++)
        setParameter(it->first, it->second);
    if(!StandardModel::Update(DPars)) return (false);
    return (true);
}


void NPbase::setParameter(const std::string name, const double& value)
{
    StandardModel::setParameter(name, value);
}


bool NPbase::CheckParameters(const std::map<std::string, double>& DPars)
{
    return(StandardModel::CheckParameters(DPars));
}


bool NPbase::setFlag(const std::string name, const bool value)
{
    return(StandardModel::setFlag(name,value));
}


bool NPbase::CheckFlags() const
{
    return(StandardModel::CheckFlags());
}


////////////////////////////////////////////////////////////////////////

double NPbase::DeltaGF() const
{
    return 0.0;
}


////////////////////////////////////////////////////////////////////////

double NPbase::obliqueS() const
{
    return 0.0;
}


double NPbase::obliqueT() const
{
    return 0.0;
}


double NPbase::obliqueU() const
{
    return 0.0;
}


////////////////////////////////////////////////////////////////////////

double NPbase::Mw() const
{
    double myMw = myEWSM->Mw_SM();

    double alpha = StandardModel::alphaMz();
    double c2 = myEWSM->cW2_SM();
    double s2 = myEWSM->sW2_SM();

    myMw *= 1.0 - alpha/4.0/(c2-s2)
                  *( obliqueS() - 2.0*c2*obliqueT() - (c2-s2)*obliqueU()/2.0/s2 )
            - s2/2.0/(c2-s2)*DeltaGF();

    //std::cout << "Mw: c_S=" << - alpha/4.0/(c2-s2) << std::endl;
    //std::cout << "Mw: c_T=" << - alpha/4.0/(c2-s2)*(- 2.0*c2) << std::endl;
    //std::cout << "Mw: c_U=" << - alpha/4.0/(c2-s2)*(- (c2-s2)/2.0/s2) << std::endl;

    return myMw;
}


double NPbase::cW2() const
{
    return ( Mw()*Mw()/Mz/Mz );
}


double NPbase::sW2() const
{
    return ( 1.0 - cW2() );
}


double NPbase::GammaW() const
{
    double Gamma_W = myEWSM->GammaW_SM();

    double alpha = StandardModel::alphaMz();
    double c2 = myEWSM->cW2_SM();
    double s2 = myEWSM->sW2_SM();

    Gamma_W *= 1.0 - 3.0*alpha/4.0/(c2-s2)
                     *( obliqueS() - 2.0*c2*obliqueT() - (c2-s2)*obliqueU()/2.0/s2 )
               - (1.0 + c2)/2.0/(c2-s2)*DeltaGF();

    //std::cout << "Gw: c_S=" << - 3.0*alpha/4.0/(c2-s2) << std::endl;
    //std::cout << "Gw: c_T=" << - 3.0*alpha/4.0/(c2-s2)*(- 2.0*c2) << std::endl;
    //std::cout << "Gw: c_U=" << - 3.0*alpha/4.0/(c2-s2)*(- (c2-s2)/2.0/s2) << std::endl;

    return Gamma_W;
}


double NPbase::deltaGVl(StandardModel::lepton l) const
{
    /* SM values */
    double alpha = StandardModel::alphaMz();
    double sW2SM = myEWSM->sW2_SM();
    double cW2SM = myEWSM->cW2_SM();
    double gVSM = myEWSM->gVl_SM(l).real();
    double gASM = myEWSM->gAl_SM(l).real();

    return ( gVSM*(alpha*obliqueT() - DeltaGF())/2.0
             + (gVSM - gASM)/4.0/sW2SM/(cW2SM - sW2SM)
               *( alpha*(obliqueS() - 4.0*cW2SM*sW2SM*obliqueT())
                  + 4.0*cW2SM*sW2SM*DeltaGF() ) );
}


double NPbase::deltaGVq(QCD::quark q) const
{
    if (q==TOP) return 0.0;

    /* SM values */
    double alpha = StandardModel::alphaMz();
    double sW2SM = myEWSM->sW2_SM();
    double cW2SM = myEWSM->cW2_SM();
    double gVSM = myEWSM->gVq_SM(q).real();
    double gASM = myEWSM->gAq_SM(q).real();

    return ( gVSM*(alpha*obliqueT() - DeltaGF())/2.0
             + (gVSM - gASM)/4.0/sW2SM/(cW2SM - sW2SM)
               *( alpha*(obliqueS() - 4.0*cW2SM*sW2SM*obliqueT())
                  + 4.0*cW2SM*sW2SM*DeltaGF() ) );
}


double NPbase::deltaGAl(StandardModel::lepton l) const
{
    /* SM values */
    double alpha = StandardModel::alphaMz();
    double gASM = myEWSM->gAl_SM(l).real();

    return ( gASM*(alpha*obliqueT() - DeltaGF())/2.0 );

}


double NPbase::deltaGAq(QCD::quark q) const
{
    if (q==TOP) return 0.0;

    /* SM values */
    double alpha = StandardModel::alphaMz();
    double gASM = myEWSM->gAq_SM(q).real();

    return ( gASM*(alpha*obliqueT() - DeltaGF())/2.0 );
}


