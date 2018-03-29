/* 
 * Copyright (C) 2014 HEPfit Collaboration
 *
 *
 * For the licensing terms see doc/COPYING.
 */

#include "HiggsKi.h"

const std::string HiggsKi::HKvKfgenvars[NHKvKfgenvars] = {
    "Kv", "Kg", "Kga", "Ku", "Kd", "Kl", "BrHinv"
};

HiggsKi::HiggsKi()
: NPbase()
{
    ModelParamMap.insert(std::pair<std::string, boost::reference_wrapper<const double> >("Kv", boost::cref(Kv)));
    ModelParamMap.insert(std::pair<std::string, boost::reference_wrapper<const double> >("Kg", boost::cref(Kg)));
    ModelParamMap.insert(std::pair<std::string, boost::reference_wrapper<const double> >("Kga", boost::cref(Kga)));
    ModelParamMap.insert(std::pair<std::string, boost::reference_wrapper<const double> >("Ku", boost::cref(Ku)));
    ModelParamMap.insert(std::pair<std::string, boost::reference_wrapper<const double> >("Kd", boost::cref(Kd)));
    ModelParamMap.insert(std::pair<std::string, boost::reference_wrapper<const double> >("Kl", boost::cref(Kl)));
    ModelParamMap.insert(std::pair<std::string, boost::reference_wrapper<const double> >("BrHinv", boost::cref(BrHinv)));
}

void HiggsKi::setParameter(const std::string name, const double& value)
{
    if (name.compare("Kv") == 0)
        Kv = value;
    else if (name.compare("Kg") == 0)
        Kg = value;
    else if (name.compare("Kga") == 0)
        Kga = value;
    else if (name.compare("Ku") == 0)
        Ku = value;
    else if (name.compare("Kd") == 0)
        Kd = value;
    else if (name.compare("Kl") == 0)
        Kl = value;
    else if (name.compare("BrHinv") == 0)
        BrHinv = value;
    else
        NPbase::setParameter(name, value);
}

bool HiggsKi::CheckParameters(const std::map<std::string, double>& DPars)
{
    for (int i = 0; i < NHKvKfgenvars; i++) {
        if (DPars.find(HKvKfgenvars[i]) == DPars.end()) {
            std::cout << "missing mandatory HiggsKi parameter " << HKvKfgenvars[i] << std::endl;
            return false;
        }
    }
    return (NPbase::CheckParameters(DPars));
}

////////////////////////////////////////////////////////////////////////

double HiggsKi::obliqueS() const
{
    double Lambda;
    if (fabs(1.0 - Kv * Kv) < pow(10.0, -32.0))
        Lambda = pow(10.0, 19.0);
    else
        Lambda = 4.0 * M_PI * v() / sqrt(fabs(1.0 - Kv * Kv));

    return ( 1.0 / 12.0 / M_PI * (1.0 - Kv * Kv) * log(Lambda * Lambda / mHl / mHl));
}

double HiggsKi::obliqueT() const
{
    double Lambda;
    double cW2_SM = trueSM.cW2();
    if (fabs(1.0 - Kv * Kv) < pow(10.0, -32.0))
        Lambda = pow(10.0, 19.0);
    else
        Lambda = 4.0 * M_PI * v() / sqrt(fabs(1.0 - Kv * Kv));

    return ( -3.0 / 16.0 / M_PI / cW2_SM * (1.0 - Kv * Kv) * log(Lambda * Lambda / mHl / mHl));
}

double HiggsKi::obliqueU() const
{
    return 0.0;
}

////////////////////////////////////////////////////////////////////////

double HiggsKi::muggH(const double sqrt_s) const
{
    return (computeKg() * computeKg());
}

double HiggsKi::muVBF(const double sqrt_s) const
{
    double sigmaWF_SM = trueSM.computeSigmaWF(sqrt_s);
    double sigmaZF_SM = trueSM.computeSigmaZF(sqrt_s);
    double sigmaZWF_SM = trueSM.computeSigmaZWF(sqrt_s);
    return (computeKW() * computeKW() * sigmaWF_SM
            + computeKZ() * computeKZ() * sigmaZF_SM
            + computeKW() * computeKZ() * sigmaZWF_SM)
            / (sigmaWF_SM + sigmaZF_SM + sigmaZWF_SM);
}

double HiggsKi::mueeWBF(const double sqrt_s) const
{
    return (computeKW() * computeKW());
}

double HiggsKi::muWH(const double sqrt_s) const
{
    return (computeKW() * computeKW());
}

double HiggsKi::muZH(const double sqrt_s) const
{
    return (computeKZ() * computeKZ());
}

double HiggsKi::mueeZH(const double sqrt_s) const
{
    return (computeKZ() * computeKZ());
}

double HiggsKi::muVH(const double sqrt_s) const
{
    double sigmaWH_SM = trueSM.computeSigmaWH(sqrt_s);
    double sigmaZH_SM = trueSM.computeSigmaZH(sqrt_s);
    return ((computeKW() * computeKW() * sigmaWH_SM
            + computeKZ() * computeKZ() * sigmaZH_SM)
            / (sigmaWH_SM + sigmaZH_SM));
}

double HiggsKi::muVBFpVH(const double sqrt_s) const
{
    double sigmaWH_SM = trueSM.computeSigmaWH(sqrt_s);
    double sigmaZH_SM = trueSM.computeSigmaZH(sqrt_s);
    double sigmaWF_SM = trueSM.computeSigmaWF(sqrt_s);
    double sigmaZF_SM = trueSM.computeSigmaZF(sqrt_s);
    double sigmaZWF_SM = trueSM.computeSigmaZWF(sqrt_s);    
    double sigmaVBF_SM = sigmaWF_SM + sigmaZF_SM + sigmaZWF_SM;
    
    double sigmaWH = muWH(sqrt_s) * sigmaWH_SM;
    double sigmaZH = muZH(sqrt_s) * sigmaZH_SM;
    double sigmaVBF = muVBF(sqrt_s) * sigmaVBF_SM;

    return ((sigmaWH + sigmaZH + sigmaVBF) / (sigmaWH_SM + sigmaZH_SM + sigmaVBF_SM));
}

double HiggsKi::muttH(const double sqrt_s) const
{
    return (computeKt() * computeKt());
}

double HiggsKi::muggHpttH(const double sqrt_s) const
{
    double sigmaggH_SM = trueSM.computeSigmaggH(sqrt_s);
    double sigmattH_SM = trueSM.computeSigmattH(sqrt_s);

    double sigmaggH = muggH(sqrt_s) * sigmaggH_SM;
    double sigmattH = muttH(sqrt_s) * sigmattH_SM;

    return ((sigmaggH + sigmattH) / (sigmaggH_SM + sigmattH_SM));
} 

double HiggsKi::mueettH(const double sqrt_s) const
{
    return (computeKt() * computeKt());
}

double HiggsKi::BrHggRatio() const
{
    return (computeKg() * computeKg() / computeGammaTotalRatio());
}

double HiggsKi::BrHWWRatio() const
{
    return (computeKW() * computeKW() / computeGammaTotalRatio());
}

double HiggsKi::BrHZZRatio() const
{
    return (computeKZ() * computeKZ() / computeGammaTotalRatio());
}

double HiggsKi::BrHZgaRatio() const
{
    return (computeKZga() * computeKZga() / computeGammaTotalRatio());
}

double HiggsKi::BrHgagaRatio() const
{
    return (computeKgaga() * computeKgaga() / computeGammaTotalRatio());
}

double HiggsKi::BrHmumuRatio() const
{
    return (computeKmu() * computeKmu() / computeGammaTotalRatio());
}

double HiggsKi::BrHtautauRatio() const
{
    return (computeKtau() * computeKtau() / computeGammaTotalRatio());
}

double HiggsKi::BrHccRatio() const
{
    return (computeKc() * computeKc() / computeGammaTotalRatio());
}

double HiggsKi::BrHbbRatio() const
{
    return (computeKb() * computeKb() / computeGammaTotalRatio());
}

double HiggsKi::muttHWW(const double sqrt_s) const
{
    return computeKt() * computeKt() * computeKW() * computeKW() / computeGammaTotalRatio();
}

double HiggsKi::muttHtautau(const double sqrt_s) const
{
    return computeKt() * computeKt() * computeKtau() * computeKtau() / computeGammaTotalRatio();
}

double HiggsKi::computeGammaTotalRatio() const
{
    return ((computeKg() * computeKg() * trueSM.computeBrHtogg()
            + computeKW() * computeKW() * trueSM.computeBrHtoWW()
            + computeKZ() * computeKZ() * trueSM.computeBrHtoZZ()
            + computeKZga() * computeKZga() * trueSM.computeBrHtoZga()
            + computeKgaga() * computeKgaga() * trueSM.computeBrHtogaga()
            + computeKmu() * computeKmu() * trueSM.computeBrHtomumu()
            + computeKtau() * computeKtau() * trueSM.computeBrHtotautau()
            + computeKc() * computeKc() * trueSM.computeBrHtocc()
            + computeKb() * computeKb() * trueSM.computeBrHtobb()
            + BrHinv)
            / (trueSM.computeBrHtogg()
            + trueSM.computeBrHtoWW()
            + trueSM.computeBrHtoZZ()
            + trueSM.computeBrHtoZga()
            + trueSM.computeBrHtogaga()
            + trueSM.computeBrHtomumu()
            + trueSM.computeBrHtotautau()
            + trueSM.computeBrHtocc()
            + trueSM.computeBrHtobb()));
}

////////////////////////////////////////////////////////////////////////

double HiggsKi::computeKg() const
{
    return Kg;
}

double HiggsKi::computeKW() const
{
    return Kv;
}

double HiggsKi::computeKZ() const
{
    return Kv;
}

double HiggsKi::computeKZga() const
{
    return 1.;
}

double HiggsKi::computeKgaga() const
{
    return Kga;
}

double HiggsKi::computeKmu() const
{
    return Kl;
}

double HiggsKi::computeKtau() const
{
    return Kl;
}

double HiggsKi::computeKc() const
{
    return Ku;
}

double HiggsKi::computeKt() const
{
    return Ku;
}

double HiggsKi::computeKb() const
{
    return Kd;
}


