/*
 * Copyright (C) 2013 HEPfit Collaboration
 *
 *
 * For the licensing terms see doc/COPYING.
 */

#include "NPHiggs.h"
#include <stdexcept>


const std::string NPHiggs::NPHIGGSvars[NNPHIGGSvars]
        = {"a", "b", "c_u", "c_d", "c_e", "d_3", "d_4", "LambdaNP"};

NPHiggs::NPHiggs()
: NPbase(), LambdaNP_in(0.0)
{
    ModelParamMap.insert(std::make_pair("a", std::cref(a)));
    ModelParamMap.insert(std::make_pair("b", std::cref(b)));
    ModelParamMap.insert(std::make_pair("c_u", std::cref(c_u)));
    ModelParamMap.insert(std::make_pair("c_d", std::cref(c_d)));
    ModelParamMap.insert(std::make_pair("c_e", std::cref(c_e)));
    ModelParamMap.insert(std::make_pair("d_3", std::cref(d_3)));
    ModelParamMap.insert(std::make_pair("d_4", std::cref(d_4)));
    ModelParamMap.insert(std::make_pair("LambdaNP", std::cref(LambdaNP_in)));

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
            raiseMissingModelParameterCount();
            addMissingModelParameter(NPHIGGSvars[i]);
        }
    }
    return (NPbase::CheckParameters(DPars));
}

////////////////////////////////////////////////////////////////////////

const double NPHiggs::obliqueS() const
{
    double Lambda;
    if (LambdaNP_in != 0.0)
        Lambda = LambdaNP_in;
    else if (fabs(1.0 - a * a) < pow(10.0, -32.0))
        Lambda = pow(10.0, 19.0);
    else
        Lambda = 4.0 * M_PI * v() / sqrt(fabs(1.0 - a * a));

    return ( 1.0 / 12.0 / M_PI * (1.0 - a * a) * log(Lambda * Lambda / mHl / mHl));
}

const double NPHiggs::obliqueT() const
{
    double Lambda;
    double cW2_SM = trueSM.cW2(); /* This has to be the SM value. */
    if (LambdaNP_in != 0.0)
        Lambda = LambdaNP_in;
    else if (fabs(1.0 - a * a) < pow(10.0, -32.0))
        Lambda = pow(10.0, 19.0);
    else
        Lambda = 4.0 * M_PI * v() / sqrt(fabs(1.0 - a * a));

    return ( -3.0 / 16.0 / M_PI / cW2_SM * (1.0 - a * a) * log(Lambda * Lambda / mHl / mHl));
}

const double NPHiggs::obliqueU() const
{
    return 0.0;
}


