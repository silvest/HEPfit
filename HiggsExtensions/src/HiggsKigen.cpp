/* 
 * Copyright (C) 2014 HEPfit Collaboration
 *
 *
 * For the licensing terms see doc/COPYING.
 */

#include "HiggsKigen.h"

const std::string HiggsKigen::HKvKfgenvars[NHKvKfgenvars] = {
    "Kw", "Kz", "Kg", "Kga", "Kzga", "Ku", "Kc", "Kt", "Kd", "Ks", "Kb", "Ke", "Kmu", "Ktau", "BrHinv"
};

HiggsKigen::HiggsKigen()
: NPbase()
{
    ModelParamMap.insert(std::pair<std::string, boost::reference_wrapper<const double> >("Kw", boost::cref(Kw)));
    ModelParamMap.insert(std::pair<std::string, boost::reference_wrapper<const double> >("Kz", boost::cref(Kz)));
    ModelParamMap.insert(std::pair<std::string, boost::reference_wrapper<const double> >("Kg", boost::cref(Kg)));
    ModelParamMap.insert(std::pair<std::string, boost::reference_wrapper<const double> >("Kga", boost::cref(Kga)));
    ModelParamMap.insert(std::pair<std::string, boost::reference_wrapper<const double> >("Kzga", boost::cref(Kzga)));
    ModelParamMap.insert(std::pair<std::string, boost::reference_wrapper<const double> >("Ku", boost::cref(Ku)));
    ModelParamMap.insert(std::pair<std::string, boost::reference_wrapper<const double> >("Kc", boost::cref(Kc)));
    ModelParamMap.insert(std::pair<std::string, boost::reference_wrapper<const double> >("Kt", boost::cref(Kt)));
    ModelParamMap.insert(std::pair<std::string, boost::reference_wrapper<const double> >("Kd", boost::cref(Kd)));
    ModelParamMap.insert(std::pair<std::string, boost::reference_wrapper<const double> >("Ks", boost::cref(Ks)));
    ModelParamMap.insert(std::pair<std::string, boost::reference_wrapper<const double> >("Kb", boost::cref(Kb)));
    ModelParamMap.insert(std::pair<std::string, boost::reference_wrapper<const double> >("Ke", boost::cref(Ke)));
    ModelParamMap.insert(std::pair<std::string, boost::reference_wrapper<const double> >("Kmu", boost::cref(Kmu)));
    ModelParamMap.insert(std::pair<std::string, boost::reference_wrapper<const double> >("Ktau", boost::cref(Ktau)));
    ModelParamMap.insert(std::pair<std::string, boost::reference_wrapper<const double> >("BrHinv", boost::cref(BrHinv)));
}

void HiggsKigen::setParameter(const std::string name, const double& value)
{
    if (name.compare("Kw") == 0)
        Kw = value;
    else if (name.compare("Kz") == 0)
        Kz = value;
    else if (name.compare("Kg") == 0)
        Kg = value;
    else if (name.compare("Kga") == 0)
        Kga = value;
    else if (name.compare("Kzga") == 0)
        Kzga = value;
    else if (name.compare("Ku") == 0)
        Ku = value;
    else if (name.compare("Kc") == 0)
        Kc = value;
    else if (name.compare("Kt") == 0)
        Kt = value;
    else if (name.compare("Kd") == 0)
        Kd = value;
    else if (name.compare("Ks") == 0)
        Ks = value;
    else if (name.compare("Kb") == 0)
        Kb = value;
    else if (name.compare("Ke") == 0)
        Ke = value;
    else if (name.compare("Kmu") == 0)
        Kmu = value;
    else if (name.compare("Ktau") == 0)
        Ktau = value;
    else if (name.compare("BrHinv") == 0)
        BrHinv = value;
    else
        NPbase::setParameter(name, value);
}

bool HiggsKigen::CheckParameters(const std::map<std::string, double>& DPars)
{
    for (int i = 0; i < NHKvKfgenvars; i++) {
        if (DPars.find(HKvKfgenvars[i]) == DPars.end()) {
            std::cout << "missing mandatory HiggsKigen parameter " << HKvKfgenvars[i] << std::endl;
            return false;
        }
    }
    return (NPbase::CheckParameters(DPars));
}

////////////////////////////////////////////////////////////////////////

double HiggsKigen::obliqueS() const
{
    return 0.0;
}

double HiggsKigen::obliqueT() const
{
    return 0.0;
}

double HiggsKigen::obliqueU() const
{
    return 0.0;
}

////////////////////////////////////////////////////////////////////////

double HiggsKigen::muggH(const double sqrt_s) const
{
    return (computeKg() * computeKg());
}

double HiggsKigen::muVBF(const double sqrt_s) const
{
    double sigmaWF_SM = trueSM.computeSigmaWF(sqrt_s);
    double sigmaZF_SM = trueSM.computeSigmaZF(sqrt_s);
    double sigmaZWF_SM = trueSM.computeSigmaZWF(sqrt_s);
    return (computeKW() * computeKW() * sigmaWF_SM
            + computeKZ() * computeKZ() * sigmaZF_SM
            + computeKW() * computeKZ() * sigmaZWF_SM)
            / (sigmaWF_SM + sigmaZF_SM + sigmaZWF_SM);
}

double HiggsKigen::mueeWBF(const double sqrt_s) const
{
    return (computeKW() * computeKW());
}

double HiggsKigen::muepWBF(const double sqrt_s) const
{
    return (computeKW() * computeKW());
}

double HiggsKigen::muepZBF(const double sqrt_s) const
{
    return (computeKZ() * computeKZ());
}

double HiggsKigen::muWH(const double sqrt_s) const
{
    return (computeKW() * computeKW());
}

double HiggsKigen::muZH(const double sqrt_s) const
{
    return (computeKZ() * computeKZ());
}

double HiggsKigen::mueeZH(const double sqrt_s) const
{
    return (computeKZ() * computeKZ());
}

double HiggsKigen::muVH(const double sqrt_s) const
{
    double sigmaWH_SM = trueSM.computeSigmaWH(sqrt_s);
    double sigmaZH_SM = trueSM.computeSigmaZH(sqrt_s);
    return ((computeKW() * computeKW() * sigmaWH_SM
            + computeKZ() * computeKZ() * sigmaZH_SM)
            / (sigmaWH_SM + sigmaZH_SM));
}

double HiggsKigen::muVBFpVH(const double sqrt_s) const
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

double HiggsKigen::muttH(const double sqrt_s) const
{
    return (computeKt() * computeKt());
}

double HiggsKigen::muggHpttH(const double sqrt_s) const
{
    double sigmaggH_SM = trueSM.computeSigmaggH(sqrt_s);
    double sigmattH_SM = trueSM.computeSigmattH(sqrt_s);

    double sigmaggH = muggH(sqrt_s) * sigmaggH_SM;
    double sigmattH = muttH(sqrt_s) * sigmattH_SM;

    return ((sigmaggH + sigmattH) / (sigmaggH_SM + sigmattH_SM));
} 

double HiggsKigen::mueettH(const double sqrt_s) const
{
    return (computeKt() * computeKt());
}

double HiggsKigen::BrHggRatio() const
{
    return (computeKg() * computeKg() / computeGammaTotalRatio());
}

double HiggsKigen::BrHWWRatio() const
{
    return (computeKW() * computeKW() / computeGammaTotalRatio());
}

double HiggsKigen::BrHZZRatio() const
{
    return (computeKZ() * computeKZ() / computeGammaTotalRatio());
}

double HiggsKigen::BrHZgaRatio() const
{
    return (computeKZga() * computeKZga() / computeGammaTotalRatio());
}

double HiggsKigen::BrHgagaRatio() const
{
    return (computeKgaga() * computeKgaga() / computeGammaTotalRatio());
}

double HiggsKigen::BrHmumuRatio() const
{
    return (computeKmu() * computeKmu() / computeGammaTotalRatio());
}

double HiggsKigen::BrHtautauRatio() const
{
    return (computeKtau() * computeKtau() / computeGammaTotalRatio());
}

double HiggsKigen::BrHccRatio() const
{
    return (computeKc() * computeKc() / computeGammaTotalRatio());
}

double HiggsKigen::BrHbbRatio() const
{
    return (computeKb() * computeKb() / computeGammaTotalRatio());
}

double HiggsKigen::muttHZbbboost(const double sqrt_s) const
{    
    return computeKt() * computeKt() * computeKb() * computeKb() / computeGammaTotalRatio();
    
}

double HiggsKigen::muttHWW(const double sqrt_s) const
{
    return computeKt() * computeKt() * computeKW() * computeKW() / computeGammaTotalRatio();
}

double HiggsKigen::muttHtautau(const double sqrt_s) const
{
    return computeKt() * computeKt() * computeKtau() * computeKtau() / computeGammaTotalRatio();
}

double HiggsKigen::computeGammaTotalRatio() const
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

double HiggsKigen::computeKg() const
{
    return Kg;
}

double HiggsKigen::computeKW() const
{
    return Kw;
}

double HiggsKigen::computeKZ() const
{
    return Kz;
}

double HiggsKigen::computeKZga() const
{
    return Kzga;
}

double HiggsKigen::computeKgaga() const
{
    return Kga;
}

double HiggsKigen::computeKmu() const
{
    return Kmu;
}

double HiggsKigen::computeKtau() const
{
    return Ktau;
}

double HiggsKigen::computeKc() const
{
    return Kc;
}

double HiggsKigen::computeKt() const
{
    return Kt;
}

double HiggsKigen::computeKb() const
{
    return Kb;
}


