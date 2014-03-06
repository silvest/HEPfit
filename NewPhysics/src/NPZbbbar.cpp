/*
 * Copyright (C) 2013 SusyFit Collaboration
 * All rights reserved.
 *
 * For the licensing terms see doc/COPYING.
 */

#include <stdexcept>
#include <EWSM.h>
#include "NPZbbbar.h"
#include "EWNPZbbbar.h"


const std::string NPZbbbar::ZbbbarVAVars[NZbbbarVars]
= {"deltaGVb", "deltaGAb"};

const std::string NPZbbbar::ZbbbarLRVars[NZbbbarVars]
= {"deltaGLb", "deltaGRb"};


NPZbbbar::NPZbbbar(const bool FlagNPZbbbarLR_in)
: NPbase(), FlagNPZbbbarLR(FlagNPZbbbarLR_in)
{
    FlagNotLinearizedNP = false;
}


bool NPZbbbar::InitializeModel()
{
    /* do not use setModelInitialized(NPbase::InitializeModel()) in order to
     use EWNPZbbbar */
    std::cout << "Model: " << ModelName() << std::endl;
    myEWSM = new EWNPZbbbar(*this);
    setModelInitialized(true);
    return(IsModelInitialized());
}


bool NPZbbbar::Init(const std::map<std::string, double>& DPars)
{
    Update(DPars);
    return(CheckParameters(DPars));
}


bool NPZbbbar::Update(const std::map<std::string,double>& DPars)
{
    for (std::map<std::string, double>::const_iterator it = DPars.begin(); it != DPars.end(); it++)
        setParameter(it->first, it->second);
    if(!NPbase::Update(DPars)) return (false);
    return (true);
}


void NPZbbbar::setParameter(const std::string name, const double& value)
{
    if (FlagNPZbbbarLR) {
        if (name.compare("deltaGLb") == 0)
            myDeltaGVb = value;
        else if (name.compare("deltaGRb") == 0)
            myDeltaGAb = value;
        else
            NPbase::setParameter(name, value);
    } else {
        if (name.compare("deltaGVb") == 0)
            myDeltaGVb = value;
        else if (name.compare("deltaGAb") == 0)
            myDeltaGAb = value;
        else
            NPbase::setParameter(name, value);
    }
}


bool NPZbbbar::CheckParameters(const std::map<std::string, double>& DPars)
{
    for (int i = 0; i < NZbbbarVars; i++) {
        if (FlagNPZbbbarLR) {
            if (DPars.find(ZbbbarLRVars[i]) == DPars.end()) {
                std::cout << "ERROR: Missing mandatory NPZbbbar parameter "
                        << ZbbbarLRVars[i] << std::endl;
                return false;
            }
        } else {
            if (DPars.find(ZbbbarVAVars[i]) == DPars.end()) {
                std::cout << "ERROR: Missing mandatory NPZbbbar parameter "
                        << ZbbbarVAVars[i] << std::endl;
                return false;
            }
        }
    }
    return(NPbase::CheckParameters(DPars));
}


bool NPZbbbar::setFlag(const std::string name, const bool value)
{
    bool res = false;
    if (name.compare("NotLinearizedNP") == 0) {
        FlagNotLinearizedNP = value;
        res = true;
    } else
        res = NPbase::setFlag(name,value);

    return(res);
}


bool NPZbbbar::CheckFlags() const
{
    return(NPbase::CheckFlags());
}


////////////////////////////////////////////////////////////////////////


double NPZbbbar::deltaGVl(StandardModel::lepton l) const
{
    return NPbase::deltaGVl(l);
}


double NPZbbbar::deltaGVq(QCD::quark q) const
{
    switch (q) {
        case QCD::UP:
        case QCD::CHARM:
        case QCD::TOP:
        case QCD::DOWN:
        case QCD::STRANGE:
            return NPbase::deltaGVq(q);
        case QCD::BOTTOM:
            if (FlagNPZbbbarLR)
                // delta g_L^b + delta g_R^b
                return ( myDeltaGVb + NPbase::deltaGVq(q)
                         + myDeltaGAb + NPbase::deltaGAq(q));
            else
                return ( myDeltaGVb + NPbase::deltaGVq(q) );
        default:
            throw std::runtime_error("Error in NPZbbbar::deltaGVq()");
    }
}


double NPZbbbar::deltaGAl(StandardModel::lepton l) const
{
    return NPbase::deltaGAl(l);
}


 double NPZbbbar::deltaGAq(QCD::quark q) const
 {
     switch (q) {
         case QCD::UP:
         case QCD::CHARM:
         case QCD::TOP:
         case QCD::DOWN:
         case QCD::STRANGE:
             return NPbase::deltaGAq(q);
         case QCD::BOTTOM:
             if (FlagNPZbbbarLR)
                // delta g_L^b - delta g_R^b
                return ( myDeltaGVb + NPbase::deltaGVq(q)
                         - myDeltaGAb - NPbase::deltaGAq(q));
             else
                 return ( myDeltaGAb + NPbase::deltaGAq(q) );
         default:
             throw std::runtime_error("Error in NPZbbbar::deltaGAq()");
     }
 }


