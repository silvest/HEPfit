/*
 * Copyright (C) 2012 HEPfit Collaboration
 *
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
    ModelParamMap.insert(std::make_pair("obliqueS", std::cref(myObliqueS)));
    ModelParamMap.insert(std::make_pair("obliqueT", std::cref(myObliqueT)));
    ModelParamMap.insert(std::make_pair("obliqueU", std::cref(myObliqueU)));
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
            raiseMissingModelParameterCount();
            addMissingModelParameter(STUvars[i]);
        }
    }
    return (NPbase::CheckParameters(DPars));
}


////////////////////////////////////////////////////////////////////////

const double NPSTU::epsilon1() const
{
    double That = alphaMz() * obliqueT();

    return ( trueSM.epsilon1() + That);
}

const double NPSTU::epsilon2() const
{
    double sW2_SM = trueSM.sW2(); /* This has to be the SM value. */
    double Uhat = -alphaMz() / (4.0 * sW2_SM) * obliqueU();

    return ( trueSM.epsilon2() + Uhat);
}

const double NPSTU::epsilon3() const
{
    double sW2_SM = trueSM.sW2(); /* This has to be the SM value. */
    double Shat = alphaMz() / (4.0 * sW2_SM) * obliqueS();

    return ( trueSM.epsilon3() + Shat);
}

const double NPSTU::epsilonb() const
{
    return trueSM.epsilonb();
}


