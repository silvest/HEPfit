/* 
 * Copyright (C) 2012 HEPfit Collaboration
 * All rights reserved.
 *
 * For the licensing terms see doc/COPYING.
 */

#include "NPEpsilons.h"
#include <stdexcept>

const std::string NPEpsilons::EPSILONvars[NEPSILONvars]
        = {"epsilon_1", "epsilon_2", "epsilon_3", "epsilon_b"};

NPEpsilons::NPEpsilons()
: NPbase()
{
    setFlagNoApproximateGammaZ(true);

    FlagEpsilon1SM = false;
    FlagEpsilon2SM = false;
    FlagEpsilon3SM = false;
    FlagEpsilonbSM = false;
    ModelParamMap.insert(std::pair<std::string, boost::reference_wrapper<const double> >("epsilon_1", boost::cref(myEpsilon_1)));
    ModelParamMap.insert(std::pair<std::string, boost::reference_wrapper<const double> >("epsilon_2", boost::cref(myEpsilon_2)));
    ModelParamMap.insert(std::pair<std::string, boost::reference_wrapper<const double> >("epsilon_3", boost::cref(myEpsilon_3)));
    ModelParamMap.insert(std::pair<std::string, boost::reference_wrapper<const double> >("epsilon_b", boost::cref(myEpsilon_b)));
}

bool NPEpsilons::PostUpdate()
{
    if (!NPbase::PostUpdate()) return (false);
    if (FlagEpsilon1SM)
        myEpsilon_1 = trueSM.epsilon1();
    if (FlagEpsilon2SM)
        myEpsilon_2 = trueSM.epsilon2();
    if (FlagEpsilon3SM)
        myEpsilon_3 = trueSM.epsilon3();
    if (FlagEpsilonbSM)
        myEpsilon_b = trueSM.epsilonb();
    return (true);
}

void NPEpsilons::setParameter(const std::string name, const double& value)
{
    if (name.compare("epsilon_1") == 0)
        myEpsilon_1 = value;
    else if (name.compare("epsilon_2") == 0)
        myEpsilon_2 = value;
    else if (name.compare("epsilon_3") == 0)
        myEpsilon_3 = value;
    else if (name.compare("epsilon_b") == 0)
        myEpsilon_b = value;
    else
        NPbase::setParameter(name, value);
}

bool NPEpsilons::CheckParameters(const std::map<std::string, double>& DPars)
{
    for (int i = 0; i < NEPSILONvars; i++) {
        if (DPars.find(EPSILONvars[i]) == DPars.end()) {
            std::cout << "ERROR: Missing mandatory NPEpsilons parameter "
                    << EPSILONvars[i] << std::endl;
            return false;
        }
    }
    return (NPbase::CheckParameters(DPars));
}

bool NPEpsilons::setFlag(const std::string name, const bool value)
{
    bool res = false;
    if (name.compare("epsilon1SM") == 0) {
        FlagEpsilon1SM = value;
        res = true;
    } else if (name.compare("epsilon2SM") == 0) {
        FlagEpsilon2SM = value;
        res = true;
    } else if (name.compare("epsilon3SM") == 0) {
        FlagEpsilon3SM = value;
        res = true;
    } else if (name.compare("epsilonbSM") == 0) {
        FlagEpsilonbSM = value;
        res = true;
    } else
        res = NPbase::setFlag(name, value);

    return (res);
}

////////////////////////////////////////////////////////////////////////

double NPEpsilons::epsilon1() const
{
    return myEpsilon_1;
}

double NPEpsilons::epsilon2() const
{
    return myEpsilon_2;
}

double NPEpsilons::epsilon3() const
{
    return myEpsilon_3;
}

double NPEpsilons::epsilonb() const
{
    return myEpsilon_b;
}


////////////////////////////////////////////////////////////////////////

double NPEpsilons::Mw_eps(const double eps1, const double eps2, const double eps3) const
{
    double Delta_r = 1.0 - (1.0 - DeltaAlpha())*(1.0 - Delta_rW(eps1, eps2, eps3));

    double tmp = 4.0 * M_PI * ale / sqrt(2.0) / GF / Mz / Mz;
    if (tmp / (1.0 - Delta_r) > 1.0)
        throw std::runtime_error("Error in EWNPEpsilons::Mw()");

    return ( Mz / sqrt(2.0) * sqrt(1.0 + sqrt(1.0 - tmp / (1.0 - Delta_r))));
}

double NPEpsilons::Mw() const
{
    return (Mw_eps(myEpsilon_1, myEpsilon_2, myEpsilon_3));
}

double NPEpsilons::GammaW() const
{
    throw std::runtime_error("NPEpsilons::GammaW() is not available.");
}

double NPEpsilons::Gamma_Z() const
{
    return (StandardModel::Gamma_Z());
}

double NPEpsilons::sigma0_had() const
{
    return (StandardModel::sigma0_had());
}

double NPEpsilons::sin2thetaEff(const Particle p) const
{
    return (StandardModel::sin2thetaEff(p));
}

double NPEpsilons::A_f(const Particle p) const
{
    return (StandardModel::A_f(p));
}

double NPEpsilons::AFB(const Particle p) const
{
    return (StandardModel::AFB(p));
}

double NPEpsilons::R0_f(const Particle p) const
{
    return (StandardModel::R0_f(p));
}


////////////////////////////////////////////////////////////////////////

gslpp::complex NPEpsilons::rhoZ_f_eps(const Particle p, const double eps1, const double epsb) const
{
    double gAe = gA_f_eps(leptons[ELECTRON], eps1).real();
    double rhoZ_e = 4. * gAe*gAe;
    double bfact = (1. + epsb)*(1. + epsb);

    if (p.is("TOP"))
        throw std::runtime_error("Error in NPEpsilons::rhoZ_f_eps()");
    if (IsFlagWithoutNonUniversalVC())
        return ( rhoZ_e * bfact);
    else
        //return ( (rhoZ_e + trueSM.deltaRhoZ_f(p).real()) * bfact); // for test
        return ( (rhoZ_e + deltaRhoZ_f(p).real()) * bfact);
}

gslpp::complex NPEpsilons::rhoZ_f(const Particle p) const
{
    if (p.is("TOP"))
        return gslpp::complex(0.0, 0.0, false);
    else if (p.is("BOTTOM"))
        return rhoZ_f_eps(p, myEpsilon_1, myEpsilon_b);
    else
        return rhoZ_f_eps(p, myEpsilon_1);
}

gslpp::complex NPEpsilons::gV_f_eps(const Particle p, const double eps1,
        const double eps3, double epsb) const
{
    double I3f = p.getIsospin();
    double Qf = p.getCharge();
    if (p.is("TOP"))
        throw std::runtime_error("Error in NPEpsilons::gV_f_eps()");
    else if (p.is("BOTTOM"))
        return ( sqrt(rhoZ_f_eps(p, eps1, epsb).real()) * I3f
            * (1.0 - 4.0 * fabs(Qf) * kappaZ_f_eps(p, eps1, eps3, epsb) * sW2()));
    else if (p.is("ELECTRON"))
        return ( (1.0 - 4.0 * (1.0 + Delta_kappaPrime(eps1, eps3)) * s02()) * gA_f_eps(leptons[ELECTRON], eps1));
    else
        return ( sqrt(rhoZ_f_eps(p, eps1).real()) * I3f
            * (1.0 - 4.0 * fabs(Qf) * kappaZ_f_eps(p, eps1, eps3) * sW2()));
}

gslpp::complex NPEpsilons::gA_f_eps(const Particle p, const double eps1, const double epsb) const
{
    double I3f = p.getIsospin();

    if (p.is("TOP"))
        throw std::runtime_error("Error in NPEpsilons::gA_f_eps()");
    else if (p.is("BOTTOM"))
        return ( sqrt(rhoZ_f_eps(p, eps1, epsb).real()) * I3f);
    else if (p.is("ELECTRON"))
        return gslpp::complex(-(1.0 + eps1 / 2.0) / 2.0, 0.0, false);
    else
        return ( gslpp::complex(sqrt(rhoZ_f_eps(p, eps1).real()) * I3f, 0.0, false));
}

gslpp::complex NPEpsilons::kappaZ_f_eps(const Particle p,
        const double eps1, const double eps3, const double epsb) const
{
    double kappaZ_e = (1.0 - (gV_f_eps(leptons[ELECTRON], eps1, eps3)
            / gA_f_eps(leptons[ELECTRON], eps1)).real()) / (4.0 * sW2());
    double bfact = (1. + epsb);

    if (p.is("TOP"))
        throw std::runtime_error("Error in NPEpsilons::kappaZ_f_eps()");
    if (IsFlagWithoutNonUniversalVC())
        return ( kappaZ_e / bfact);
    else
        //return ( (kappaZ_e + trueSM.deltaKappaZ_f(p).real()) / bfact); // for test
        return ( (kappaZ_e + deltaKappaZ_f(p).real()) / bfact);
}

gslpp::complex NPEpsilons::kappaZ_f(const Particle p) const
{
    if (p.is("TOP"))
        return gslpp::complex(0.0, 0.0, false);
    else if (p.is("BOTTOM"))
        return (kappaZ_f_eps(p, myEpsilon_1, myEpsilon_3, myEpsilon_b));
    else
        return (kappaZ_f_eps(p, myEpsilon_1, myEpsilon_3));
}

gslpp::complex NPEpsilons::gV_f(const Particle p) const
{
    if (p.is("TOP"))
        return gslpp::complex(0.0, 0.0, false);
    else if (p.is("BOTTOM"))
        return (gV_f_eps(p, myEpsilon_1, myEpsilon_3, myEpsilon_b));
    else
        return (gV_f_eps(p, myEpsilon_1, myEpsilon_3));
}

gslpp::complex NPEpsilons::gA_f(const Particle p) const
{
    if (p.is("TOP"))
        return gslpp::complex(0.0, 0.0, false);
    else if (p.is("BOTTOM"))
        return (gA_f_eps(p, myEpsilon_1, myEpsilon_b));
    else
        return (gA_f_eps(p, myEpsilon_1));
}

////////////////////////////////////////////////////////////////////////

double NPEpsilons::Delta_rW(const double eps1, const double eps2, const double eps3) const
{
    return ( (c02() - s02()) / s02()
            *(eps2 - c02() * eps1 + 2.0 * s02() * Delta_kappaPrime(eps1, eps3)));
}

double NPEpsilons::Delta_kappaPrime(const double eps1, const double eps3) const
{
    return ( (eps3 - c02() * eps1) / (c02() - s02()));
}
