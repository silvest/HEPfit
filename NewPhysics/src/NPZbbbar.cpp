/* 
 * Copyright (C) 2013-2014 SusyFit Collaboration
 * All rights reserved.
 *
 * For the licensing terms see doc/COPYING.
 */

#include <stdexcept>
#include <EWSM.h>
#include "NPZbbbar.h"
#include "EWNPZbbbar.h"


const std::string NPZbbbar::ZbbbarVars[NZbbbarVars] 
= {"deltaGVb", "deltaGAb"};


NPZbbbar::NPZbbbar() 
: NPbase()
{
    FlagNPZbbbarLR = false;
    FlagNotLinearizedNP = false;
}


bool NPZbbbar::InitializeModel()
{
    /* do not use setModelInitialized(NPbase::InitializeModel()) in order to
     use EWNPZbbbar */
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
    if (name.compare("deltaGVb") == 0)
        myDeltaGVb = value;
    else if (name.compare("deltaGAb") == 0)
        myDeltaGAb = value;
    else
        NPbase::setParameter(name, value);
}


bool NPZbbbar::CheckParameters(const std::map<std::string, double>& DPars)
{
    for (int i = 0; i < NZbbbarVars; i++) {
        if (DPars.find(ZbbbarVars[i]) == DPars.end()) {
            std::cout << "ERROR: Missing mandatory NPZbbbar parameter "
                      << ZbbbarVars[i] << std::endl;
            return false;
        }
    }
    return(NPbase::CheckParameters(DPars));
}


bool NPZbbbar::setFlag(const std::string name, const bool& value) 
{
    bool res = false;
    if (name.compare("NPZbbbarLR") == 0) {
        FlagNPZbbbarLR = value;
        res = true;
    } else if (name.compare("NotLinearizedNP") == 0) {
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


double NPZbbbar::deltaGVq(StandardModel::quark q) const
{
    switch (q) {
        case StandardModel::UP:
        case StandardModel::CHARM:
        case StandardModel::TOP:
        case StandardModel::DOWN:
        case StandardModel::STRANGE:
            return NPbase::deltaGVq(q);
        case StandardModel::BOTTOM:
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


 double NPZbbbar::deltaGAq(StandardModel::quark q) const
 {
     switch (q) {
         case StandardModel::UP:
         case StandardModel::CHARM:
         case StandardModel::TOP:
         case StandardModel::DOWN:
         case StandardModel::STRANGE:
             return NPbase::deltaGAq(q);
         case StandardModel::BOTTOM:
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

 
