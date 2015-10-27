/* 
 * Copyright (C) 2014 HEPfit Collaboration
 * All rights reserved.
 *
 * For the licensing terms see doc/COPYING.
 */

#include "HiggsKvgenKf.h"

const std::string HiggsKvgenKf::HKvgenKfvars[NHKvgenKfvars] = {
    "KW", "KZ", "Kf", "BrHinv"

};

HiggsKvgenKf::HiggsKvgenKf()
: NPbase()
{
    ModelParamMap.insert(std::pair<std::string, boost::reference_wrapper<const double> >("KW", boost::cref(KW)));
    ModelParamMap.insert(std::pair<std::string, boost::reference_wrapper<const double> >("KZ", boost::cref(KZ)));
    ModelParamMap.insert(std::pair<std::string, boost::reference_wrapper<const double> >("Kf", boost::cref(Kf)));
    ModelParamMap.insert(std::pair<std::string, boost::reference_wrapper<const double> >("BrHinv", boost::cref(BrHinv)));
}

void HiggsKvgenKf::setParameter(const std::string name, const double& value)
{
    if (name.compare("KW") == 0)
        KW = value;
    else if (name.compare("KZ") == 0)
        KZ = value;
    else if (name.compare("Kf") == 0)
        Kf = value;
    else if (name.compare("BrHinv") == 0)
        BrHinv = value;
    else
        NPbase::setParameter(name, value);
}

bool HiggsKvgenKf::CheckParameters(const std::map<std::string, double>& DPars)
{
    for (int i = 0; i < NHKvgenKfvars; i++) {
        if (DPars.find(HKvgenKfvars[i]) == DPars.end()) {
            std::cout << "missing mandatory HiggsKvgenKf parameter " << HKvgenKfvars[i] << std::endl;
            return false;
        }
    }
    return (NPbase::CheckParameters(DPars));
}

////////////////////////////////////////////////////////////////////////

double HiggsKvgenKf::muggH(const double sqrt_s) const
{
    return (computeKg() * computeKg());
}

double HiggsKvgenKf::muVBF(const double sqrt_s) const
{
    double sigmaWF_SM = trueSM.computeSigmaWF(sqrt_s);
    double sigmaZF_SM = trueSM.computeSigmaZF(sqrt_s);
    double sigmaZWF_SM = trueSM.computeSigmaZWF(sqrt_s);
    return (computeKW() * computeKW() * sigmaWF_SM
            + computeKZ() * computeKZ() * sigmaZF_SM
            + computeKW() * computeKZ() * sigmaZWF_SM)
            / (sigmaWF_SM + sigmaZF_SM + sigmaZWF_SM);
}

double HiggsKvgenKf::muWH(const double sqrt_s) const
{
    return (computeKW() * computeKW());
}

double HiggsKvgenKf::muZH(const double sqrt_s) const
{
    return (computeKZ() * computeKZ());
}

double HiggsKvgenKf::muVH(const double sqrt_s) const
{
    double sigmaWH_SM = trueSM.computeSigmaWH(sqrt_s);
    double sigmaZH_SM = trueSM.computeSigmaZH(sqrt_s);
    return ((computeKW() * computeKW() * sigmaWH_SM
            + computeKZ() * computeKZ() * sigmaZH_SM)
            / (sigmaWH_SM + sigmaZH_SM));
}

double HiggsKvgenKf::muVBFpVH(const double sqrt_s) const
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

double HiggsKvgenKf::muttH(const double sqrt_s) const
{
    return (computeKt() * computeKt());
}

double HiggsKvgenKf::muggHpttH(const double sqrt_s) const
{
    double sigmaggH_SM = trueSM.computeSigmaggH(sqrt_s);
    double sigmattH_SM = trueSM.computeSigmattH(sqrt_s);

    double sigmaggH = muggH(sqrt_s) * sigmaggH_SM;
    double sigmattH = muttH(sqrt_s) * sigmattH_SM;

    return ((sigmaggH + sigmattH) / (sigmaggH_SM + sigmattH_SM));
} 

double HiggsKvgenKf::BrHggRatio() const
{
    return (computeKg() * computeKg() / computeGammaTotalRatio());
}

double HiggsKvgenKf::BrHWWRatio() const
{
    return (computeKW() * computeKW() / computeGammaTotalRatio());
}

double HiggsKvgenKf::BrHZZRatio() const
{
    return (computeKZ() * computeKZ() / computeGammaTotalRatio());
}

double HiggsKvgenKf::BrHZgaRatio() const
{
    return (computeKZga() * computeKZga() / computeGammaTotalRatio());
}

double HiggsKvgenKf::BrHgagaRatio() const
{
    return (computeKgaga() * computeKgaga() / computeGammaTotalRatio());
}

double HiggsKvgenKf::BrHtautauRatio() const
{
    return (computeKtau() * computeKtau() / computeGammaTotalRatio());
}

double HiggsKvgenKf::BrHccRatio() const
{
    return (computeKc() * computeKc() / computeGammaTotalRatio());
}

double HiggsKvgenKf::BrHbbRatio() const
{
    return (computeKb() * computeKb() / computeGammaTotalRatio());
}

double HiggsKvgenKf::computeGammaTotalRatio() const
{
    return ((computeKg() * computeKg() * trueSM.computeBrHtogg()
            + computeKW() * computeKW() * trueSM.computeBrHtoWW()
            + computeKZ() * computeKZ() * trueSM.computeBrHtoZZ()
            + computeKZga() * computeKZga() * trueSM.computeBrHtoZga()
            + computeKgaga() * computeKgaga() * trueSM.computeBrHtogaga()
            + computeKtau() * computeKtau() * trueSM.computeBrHtotautau()
            + computeKc() * computeKc() * trueSM.computeBrHtocc()
            + computeKb() * computeKb() * trueSM.computeBrHtobb())
            / (1.0 - BrHinv));
}

////////////////////////////////////////////////////////////////////////

double HiggsKvgenKf::computeKg() const
{
    return Kf;
}

double HiggsKvgenKf::computeKW() const
{
    return KW;
}

double HiggsKvgenKf::computeKZ() const
{
    return KZ;
}

double HiggsKvgenKf::computeKZga() const
{
    double gtt_SM = trueSM.computeGammaHZga_tt();
    double gWW_SM = trueSM.computeGammaHZga_WW();
    double gtW_SM = trueSM.computeGammaHZga_tW();
    return (sqrt((computeKt() * computeKt() * gtt_SM
            + computeKW() * computeKW() * gWW_SM
            + computeKt() * computeKW() * gtW_SM) / (gtt_SM + gWW_SM + gtW_SM)));
}

double HiggsKvgenKf::computeKgaga() const
{
    double gtt_SM = trueSM.computeGammaHgaga_tt();
    double gWW_SM = trueSM.computeGammaHgaga_WW();
    double gtW_SM = trueSM.computeGammaHgaga_tW();
    return (sqrt((computeKt() * computeKt() * gtt_SM
            + computeKW() * computeKW() * gWW_SM
            + computeKt() * computeKW() * gtW_SM) / (gtt_SM + gWW_SM + gtW_SM)));
}

double HiggsKvgenKf::computeKtau() const
{
    return Kf;
}

double HiggsKvgenKf::computeKc() const
{
    return Kf;
}

double HiggsKvgenKf::computeKt() const
{
    return Kf;
}

double HiggsKvgenKf::computeKb() const
{
    return Kf;
}


