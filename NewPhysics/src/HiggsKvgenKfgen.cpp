/* 
 * Copyright (C) 2014 HEPfit Collaboration
 *
 *
 * For the licensing terms see doc/COPYING.
 */

#include "HiggsKvgenKfgen.h"

const std::string HiggsKvgenKfgen::HKvgenKfgenvars[NHKvgenKfgenvars] = {
    "KW", "KZ", "Ku", "Kd", "Kl", "BrHinv"
};

HiggsKvgenKfgen::HiggsKvgenKfgen()
: NPbase()
{
    ModelParamMap.insert(std::make_pair("KW", std::cref(KW)));
    ModelParamMap.insert(std::make_pair("KZ", std::cref(KZ)));
    ModelParamMap.insert(std::make_pair("Ku", std::cref(Ku)));
    ModelParamMap.insert(std::make_pair("Kd", std::cref(Kd)));
    ModelParamMap.insert(std::make_pair("Kl", std::cref(Kl)));
    ModelParamMap.insert(std::make_pair("BrHinv", std::cref(BrHinv)));
}

void HiggsKvgenKfgen::setParameter(const std::string name, const double& value)
{
    if (name.compare("KW") == 0)
        KW = value;
    else if (name.compare("KZ") == 0)
        KZ = value;
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

bool HiggsKvgenKfgen::CheckParameters(const std::map<std::string, double>& DPars)
{
    for (int i = 0; i < NHKvgenKfgenvars; i++) {
        if (DPars.find(HKvgenKfgenvars[i]) == DPars.end()) {
            std::cout << "ERROR: missing mandatory HiggsKvgenKfgen parameter " << HKvgenKfgenvars[i] << std::endl;
            raiseMissingModelParameterCount();
            addMissingModelParameter(HKvgenKfgenvars[i]);
        }
    }
    return (NPbase::CheckParameters(DPars));
}

////////////////////////////////////////////////////////////////////////

double HiggsKvgenKfgen::muggH(const double sqrt_s) const
{
    return (computeKg(sqrt_s) * computeKg(sqrt_s));
}

double HiggsKvgenKfgen::muVBF(const double sqrt_s) const
{
    double sigmaWF_SM = trueSM.computeSigmaWF(sqrt_s);
    double sigmaZF_SM = trueSM.computeSigmaZF(sqrt_s);
    double sigmaZWF_SM = trueSM.computeSigmaZWF(sqrt_s);
    return (computeKW() * computeKW() * sigmaWF_SM
            + computeKZ() * computeKZ() * sigmaZF_SM
            + computeKW() * computeKZ() * sigmaZWF_SM)
            / (sigmaWF_SM + sigmaZF_SM + sigmaZWF_SM);
}

double HiggsKvgenKfgen::mueeWBF(const double sqrt_s) const
{
    return (computeKW() * computeKW());
}

double HiggsKvgenKfgen::muWH(const double sqrt_s) const
{
    return (computeKW() * computeKW());
}

double HiggsKvgenKfgen::muZH(const double sqrt_s) const
{
    return (computeKZ() * computeKZ());
}

double HiggsKvgenKfgen::mueeZH(const double sqrt_s) const
{
    return (computeKZ() * computeKZ());
}

double HiggsKvgenKfgen::muVH(const double sqrt_s) const
{
    double sigmaWH_SM = trueSM.computeSigmaWH(sqrt_s);
    double sigmaZH_SM = trueSM.computeSigmaZH(sqrt_s);
    return ((computeKW() * computeKW() * sigmaWH_SM
            + computeKZ() * computeKZ() * sigmaZH_SM)
            / (sigmaWH_SM + sigmaZH_SM));
}

double HiggsKvgenKfgen::muVBFpVH(const double sqrt_s) const
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

double HiggsKvgenKfgen::muttH(const double sqrt_s) const
{
    return (computeKt() * computeKt());
}

double HiggsKvgenKfgen::muggHpttH(const double sqrt_s) const
{
    double sigmaggH_SM = trueSM.computeSigmaggH(sqrt_s);
    double sigmattH_SM = trueSM.computeSigmattH(sqrt_s);

    double sigmaggH = muggH(sqrt_s) * sigmaggH_SM;
    double sigmattH = muttH(sqrt_s) * sigmattH_SM;

    return ((sigmaggH + sigmattH) / (sigmaggH_SM + sigmattH_SM));
} 

double HiggsKvgenKfgen::mueettH(const double sqrt_s) const
{
    return (computeKt() * computeKt());
}

double HiggsKvgenKfgen::BrHggRatio() const
{
    double gtt_SM = trueSM.computeGammaHgg_tt();
    double gbb_SM = trueSM.computeGammaHgg_bb();
    double gtb_SM = trueSM.computeGammaHgg_tb();
    double Ggg_over_GggSM = (computeKt() * computeKt() * gtt_SM
            + computeKb() * computeKb() * gbb_SM
            + computeKt() * computeKb() * gtb_SM) / (gtt_SM + gbb_SM + gtb_SM);
    return (Ggg_over_GggSM / computeGammaTotalRatio());
}

double HiggsKvgenKfgen::BrHWWRatio() const
{
    return (computeKW() * computeKW() / computeGammaTotalRatio());
}

double HiggsKvgenKfgen::BrHZZRatio() const
{
    return (computeKZ() * computeKZ() / computeGammaTotalRatio());
}

double HiggsKvgenKfgen::BrHZgaRatio() const
{
    return (computeKZga() * computeKZga() / computeGammaTotalRatio());
}

double HiggsKvgenKfgen::BrHgagaRatio() const
{
    return (computeKgaga() * computeKgaga() / computeGammaTotalRatio());
}

double HiggsKvgenKfgen::BrHmumuRatio() const
{
    return (computeKmu() * computeKmu() / computeGammaTotalRatio());
}

double HiggsKvgenKfgen::BrHtautauRatio() const
{
    return (computeKtau() * computeKtau() / computeGammaTotalRatio());
}

double HiggsKvgenKfgen::BrHccRatio() const
{
    return (computeKc() * computeKc() / computeGammaTotalRatio());
}

double HiggsKvgenKfgen::BrHbbRatio() const
{
    return (computeKb() * computeKb() / computeGammaTotalRatio());
}

double HiggsKvgenKfgen::computeGammaTotalRatio() const
{
    double gtt_SM = trueSM.computeGammaHgg_tt();
    double gbb_SM = trueSM.computeGammaHgg_bb();
    double gtb_SM = trueSM.computeGammaHgg_tb();
    double Ggg_over_GggSM = (computeKt() * computeKt() * gtt_SM
            + computeKb() * computeKb() * gbb_SM
            + computeKt() * computeKb() * gtb_SM) / (gtt_SM + gbb_SM + gtb_SM);
    return ((Ggg_over_GggSM * trueSM.computeBrHtogg()
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

double HiggsKvgenKfgen::computeKg(const double sqrt_s) const
{
    if (sqrt_s == 1.96)
        return computeKt();

    double sigmatt_SM = trueSM.computeSigmaggH_tt(sqrt_s);
    double sigmabb_SM = trueSM.computeSigmaggH_bb(sqrt_s);
    double sigmatb_SM = trueSM.computeSigmaggH_tb(sqrt_s);
    return (sqrt((computeKt() * computeKt() * sigmatt_SM
            + computeKb() * computeKb() * sigmabb_SM
            + computeKt() * computeKb() * sigmatb_SM)
            / (sigmatt_SM + sigmabb_SM + sigmatb_SM)));
}

double HiggsKvgenKfgen::computeKW() const
{
    return KW;
}

double HiggsKvgenKfgen::computeKZ() const
{
    return KZ;
}

double HiggsKvgenKfgen::computeKZga() const
{
    double gtt_SM = trueSM.computeGammaHZga_tt();
    double gWW_SM = trueSM.computeGammaHZga_WW();
    double gtW_SM = trueSM.computeGammaHZga_tW();
    return (sqrt((computeKt() * computeKt() * gtt_SM
            + computeKW() * computeKW() * gWW_SM
            + computeKt() * computeKW() * gtW_SM) / (gtt_SM + gWW_SM + gtW_SM)));
}

double HiggsKvgenKfgen::computeKgaga() const
{
    double gtt_SM = trueSM.computeGammaHgaga_tt();
    double gWW_SM = trueSM.computeGammaHgaga_WW();
    double gtW_SM = trueSM.computeGammaHgaga_tW();
    return (sqrt((computeKt() * computeKt() * gtt_SM
            + computeKW() * computeKW() * gWW_SM
            + computeKt() * computeKW() * gtW_SM) / (gtt_SM + gWW_SM + gtW_SM)));
}

double HiggsKvgenKfgen::computeKmu() const
{
    return Kl;
}

double HiggsKvgenKfgen::computeKtau() const
{
    return Kl;
}

double HiggsKvgenKfgen::computeKc() const
{
    return Ku;
}

double HiggsKvgenKfgen::computeKt() const
{
    return Ku;
}

double HiggsKvgenKfgen::computeKb() const
{
    return Kd;
}


