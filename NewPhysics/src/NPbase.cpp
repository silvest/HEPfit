/*
 * Copyright (C) 2013 SusyFit Collaboration
 * All rights reserved.
 *
 * For the licensing terms see doc/COPYING.
 */

#include <EWSM.h>
#include "NPbase.h"


const std::string NPbase::NPbaseflags[NNPbaseflags]
= {"FixSMcontribution"};


NPbase::NPbase()
: StandardModel()
{
    FlagFixSMcontribution = false;
}


bool NPbase::InitializeModel()
{
    setModelInitialized(StandardModel::InitializeModel());
    return (IsModelInitialized());
}


void NPbase::setEWSMflags(EWSM& myEWSM)
{
    StandardModel::setEWSMflags(myEWSM);
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


bool NPbase::setFlag(const std::string name, const bool& value)
{
    bool res = false;
    if (name.compare("FixSMcontribution") == 0) {
        FlagFixSMcontribution = value;
        res = true;
    } else
        res = StandardModel::setFlag(name,value);

    return(res);
}


bool NPbase::CheckFlags() const
{
    return(StandardModel::CheckFlags());
}


////////////////////////////////////////////////////////////////////////

double NPbase::epsilon1() const
{
    return ( myEWSM->epsilon1_SM() + obliqueThat() );
}


double NPbase::epsilon2() const
{
    return ( myEWSM->epsilon2_SM() + obliqueUhat() );
}


double NPbase::epsilon3() const
{
    return ( myEWSM->epsilon3_SM() + obliqueShat() );
}


double NPbase::epsilonb() const
{
    return myEWSM->epsilonb_SM();
}


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


double NPbase::obliqueShat() const
{
    double sW2_SM = myEWSM->sW2_SM(); /* This has to be the SM value. */
    return ( alphaMz()/(4.0*sW2_SM)*obliqueS() );
}


double NPbase::obliqueThat() const
{
    return ( alphaMz()*obliqueT() );
}


double NPbase::obliqueUhat() const
{
    double sW2_SM = myEWSM->sW2_SM(); /* This has to be the SM value. */
    return ( - alphaMz()/(4.0*sW2_SM)*obliqueU() );
}


double NPbase::obliqueV() const
{
    return 0.0;
}


double NPbase::obliqueW() const
{
    return 0.0;
}


double NPbase::obliqueX() const
{
    return 0.0;
}


double NPbase::obliqueY() const
{
    return 0.0;
}


double NPbase::DeltaGF() const
{
    return 0.0;
}


double NPbase::deltaGVl(StandardModel::lepton l) const
{
    /* SM values */
    double alpha = StandardModel::alphaMz();
    double sW2SM = myEWSM->sW2_SM();
    double cW2SM = myEWSM->cW2_SM();
    double gVSM = myEWSM->gVl_SM(l).real();
    double gASM = myEWSM->gAl_SM(l).real();

    return ( gVSM*alpha*obliqueT()/2.0
            + (gVSM - gASM)*alpha/4.0/sW2SM/(cW2SM - sW2SM)
              *(obliqueS() - 4.0*cW2SM*sW2SM*obliqueT()) );
}

double NPbase::deltaGVq(StandardModel::quark q) const
{
    if (q==TOP) return 0.0;

    /* SM values */
    double alpha = StandardModel::alphaMz();
    double sW2SM = myEWSM->sW2_SM();
    double cW2SM = myEWSM->cW2_SM();
    double gVSM = myEWSM->gVq_SM(q).real();
    double gASM = myEWSM->gAq_SM(q).real();

    return ( gVSM*alpha*obliqueT()/2.0
             + (gVSM - gASM)*alpha/4.0/sW2SM/(cW2SM - sW2SM)
               *(obliqueS() - 4.0*cW2SM*sW2SM*obliqueT()) );
}

double NPbase::deltaGAl(StandardModel::lepton l) const
{
    /* SM values */
    double alpha = StandardModel::alphaMz();
    double gASM = myEWSM->gAl_SM(l).real();

    return ( gASM*alpha*obliqueT()/2.0 );

}

double NPbase::deltaGAq(StandardModel::quark q) const
{
    if (q==TOP) return 0.0;

    /* SM values */
    double alpha = StandardModel::alphaMz();
    double gASM = myEWSM->gAq_SM(q).real();

    return ( gASM*alpha*obliqueT()/2.0 );
}

