/*
 * Copyright (C) 2012 HEPfit Collaboration
 * All rights reserved.
 *
 * For the licensing terms see doc/COPYING.
 */

#include "NPSTU.h"
#include <stdexcept>



const std::string NPSTU::STUvars[NSTUvars]
        = {"obliqueS", "obliqueT", "obliqueU"};

NPSTU::NPSTU()
: NPbase()
{
    ModelParamMap.insert(std::pair<std::string, boost::reference_wrapper<const double> >("obliqueS", boost::cref(myObliqueS)));
    ModelParamMap.insert(std::pair<std::string, boost::reference_wrapper<const double> >("obliqueT", boost::cref(myObliqueT)));
    ModelParamMap.insert(std::pair<std::string, boost::reference_wrapper<const double> >("obliqueU", boost::cref(myObliqueU)));
}

void NPSTU::setParameter(const std::string name, const double& value)
{
    if (name.compare("obliqueS") == 0)
        myObliqueS = value;
    else if (name.compare("obliqueT") == 0)
        myObliqueT = value;
    else if (name.compare("obliqueU") == 0)
        myObliqueU = value;
    else
        NPbase::setParameter(name, value);
}

bool NPSTU::CheckParameters(const std::map<std::string, double>& DPars)
{
    for (int i = 0; i < NSTUvars; i++) {
        if (DPars.find(STUvars[i]) == DPars.end()) {
            std::cout << "ERROR: Missing mandatory NPSTU parameter "
                    << STUvars[i] << std::endl;
            return false;
        }
    }
    return (NPbase::CheckParameters(DPars));
}


////////////////////////////////////////////////////////////////////////

double NPSTU::epsilon1() const
{
    double That = alphaMz() * obliqueT();

    return ( trueSM.epsilon1() + That);
}

double NPSTU::epsilon2() const
{
    double sW2_SM = trueSM.sW2(); /* This has to be the SM value. */
    double Uhat = -alphaMz() / (4.0 * sW2_SM) * obliqueU();

    return ( trueSM.epsilon2() + Uhat);
}

double NPSTU::epsilon3() const
{
    double sW2_SM = trueSM.sW2(); /* This has to be the SM value. */
    double Shat = alphaMz() / (4.0 * sW2_SM) * obliqueS();

    return ( trueSM.epsilon3() + Shat);
}

double NPSTU::epsilonb() const
{
    return trueSM.epsilonb();
}


