/* 
 * Copyright (C) 2014 HEPfit Collaboration
 * All rights reserved.
 *
 * For the licensing terms see doc/COPYING.
 */

#include "HiggsKvKfgen.h"

const std::string HiggsKvKfgen::HKvKfgenvars[NHKvKfgenvars] = {
    "Kv", "Ku", "Kd", "Kl", "BrHinv"
};

HiggsKvKfgen::HiggsKvKfgen()
: NPbase()
{
    ModelParamMap.insert(std::pair<std::string, boost::reference_wrapper<const double> >("Kv", boost::cref(Kv)));
    ModelParamMap.insert(std::pair<std::string, boost::reference_wrapper<const double> >("Ku", boost::cref(Ku)));
    ModelParamMap.insert(std::pair<std::string, boost::reference_wrapper<const double> >("Kd", boost::cref(Kd)));
    ModelParamMap.insert(std::pair<std::string, boost::reference_wrapper<const double> >("Kl", boost::cref(Kl)));
    ModelParamMap.insert(std::pair<std::string, boost::reference_wrapper<const double> >("BrHinv", boost::cref(BrHinv)));
}

void HiggsKvKfgen::setParameter(const std::string name, const double& value)
{
    if (name.compare("Kv") == 0)
        Kv = value;
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

bool HiggsKvKfgen::CheckParameters(const std::map<std::string, double>& DPars)
{
    for (int i = 0; i < NHKvKfgenvars; i++) {
        if (DPars.find(HKvKfgenvars[i]) == DPars.end()) {
            std::cout << "missing mandatory HiggsKvKfgen parameter " << HKvKfgenvars[i] << std::endl;
            return false;
        }
    }
    return (NPbase::CheckParameters(DPars));
}

////////////////////////////////////////////////////////////////////////

double HiggsKvKfgen::obliqueS() const
{
    double Lambda;
    if (fabs(1.0 - Kv * Kv) < pow(10.0, -32.0))
        Lambda = pow(10.0, 19.0);
    else
        Lambda = 4.0 * M_PI * v() / sqrt(fabs(1.0 - Kv * Kv));

    return ( 1.0 / 12.0 / M_PI * (1.0 - Kv * Kv) * log(Lambda * Lambda / mHl / mHl));
}

double HiggsKvKfgen::obliqueT() const
{
    double Lambda;
    double cW2_SM = trueSM.cW2();
    if (fabs(1.0 - Kv * Kv) < pow(10.0, -32.0))
        Lambda = pow(10.0, 19.0);
    else
        Lambda = 4.0 * M_PI * v() / sqrt(fabs(1.0 - Kv * Kv));

    return ( -3.0 / 16.0 / M_PI / cW2_SM * (1.0 - Kv * Kv) * log(Lambda * Lambda / mHl / mHl));
}

double HiggsKvKfgen::obliqueU() const
{
    return 0.0;
}

////////////////////////////////////////////////////////////////////////

double HiggsKvKfgen::muggH(const double sqrt_s) const
{
    return (computeKg(sqrt_s) * computeKg(sqrt_s));
}

double HiggsKvKfgen::muVBF(const double sqrt_s) const
{
    double sigmaWF_SM = trueSM.computeSigmaWF(sqrt_s);
    double sigmaZF_SM = trueSM.computeSigmaZF(sqrt_s);
    double sigmaZWF_SM = trueSM.computeSigmaZWF(sqrt_s);
    return (computeKW() * computeKW() * sigmaWF_SM
            + computeKZ() * computeKZ() * sigmaZF_SM
            + computeKW() * computeKZ() * sigmaZWF_SM)
            / (sigmaWF_SM + sigmaZF_SM + sigmaZWF_SM);
}

double HiggsKvKfgen::muWH(const double sqrt_s) const
{
    return (computeKW() * computeKW());
}

double HiggsKvKfgen::muZH(const double sqrt_s) const
{
    return (computeKZ() * computeKZ());
}

double HiggsKvKfgen::muVH(const double sqrt_s) const
{
    double sigmaWH_SM = trueSM.computeSigmaWH(sqrt_s);
    double sigmaZH_SM = trueSM.computeSigmaZH(sqrt_s);
    return ((computeKW() * computeKW() * sigmaWH_SM
            + computeKZ() * computeKZ() * sigmaZH_SM)
            / (sigmaWH_SM + sigmaZH_SM));
}

double HiggsKvKfgen::muVBFpVH(const double sqrt_s) const
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

double HiggsKvKfgen::muttH(const double sqrt_s) const
{
    return (computeKt() * computeKt());
}

double HiggsKvKfgen::muggHpttH(const double sqrt_s) const
{
    double sigmaggH_SM = trueSM.computeSigmaggH(sqrt_s);
    double sigmattH_SM = trueSM.computeSigmattH(sqrt_s);

    double sigmaggH = muggH(sqrt_s) * sigmaggH_SM;
    double sigmattH = muttH(sqrt_s) * sigmattH_SM;

    return ((sigmaggH + sigmattH) / (sigmaggH_SM + sigmattH_SM));
} 

double HiggsKvKfgen::BrHggRatio() const
{
    double gtt_SM = trueSM.computeGammaHgg_tt();
    double gbb_SM = trueSM.computeGammaHgg_bb();
    double gtb_SM = trueSM.computeGammaHgg_tb();
    double Ggg_over_GggSM = (computeKt() * computeKt() * gtt_SM
            + computeKb() * computeKb() * gbb_SM
            + computeKt() * computeKb() * gtb_SM) / (gtt_SM + gbb_SM + gtb_SM);
    return (Ggg_over_GggSM / computeGammaTotalRatio());
}

double HiggsKvKfgen::BrHWWRatio() const
{
    return (computeKW() * computeKW() / computeGammaTotalRatio());
}

double HiggsKvKfgen::BrHZZRatio() const
{
    return (computeKZ() * computeKZ() / computeGammaTotalRatio());
}

double HiggsKvKfgen::BrHZgaRatio() const
{
    return (computeKZga() * computeKZga() / computeGammaTotalRatio());
}

double HiggsKvKfgen::BrHgagaRatio() const
{
    return (computeKgaga() * computeKgaga() / computeGammaTotalRatio());
}

double HiggsKvKfgen::BrHtautauRatio() const
{
    return (computeKtau() * computeKtau() / computeGammaTotalRatio());
}

double HiggsKvKfgen::BrHccRatio() const
{
    return (computeKc() * computeKc() / computeGammaTotalRatio());
}

double HiggsKvKfgen::BrHbbRatio() const
{
    return (computeKb() * computeKb() / computeGammaTotalRatio());
}

double HiggsKvKfgen::computeGammaTotalRatio() const
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
            + computeKtau() * computeKtau() * trueSM.computeBrHtotautau()
            + computeKc() * computeKc() * trueSM.computeBrHtocc()
            + computeKb() * computeKb() * trueSM.computeBrHtobb())
            / (1.0 - BrHinv));
}

////////////////////////////////////////////////////////////////////////

double HiggsKvKfgen::computeKg(const double sqrt_s) const
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

double HiggsKvKfgen::computeKW() const
{
    return Kv;
}

double HiggsKvKfgen::computeKZ() const
{
    return Kv;
}

double HiggsKvKfgen::computeKZga() const
{
    double gtt_SM = trueSM.computeGammaHZga_tt();
    double gWW_SM = trueSM.computeGammaHZga_WW();
    double gtW_SM = trueSM.computeGammaHZga_tW();
    return (sqrt((computeKt() * computeKt() * gtt_SM
            + computeKW() * computeKW() * gWW_SM
            + computeKt() * computeKW() * gtW_SM) / (gtt_SM + gWW_SM + gtW_SM)));
}

double HiggsKvKfgen::computeKgaga() const
{
    double gtt_SM = trueSM.computeGammaHgaga_tt();
    double gWW_SM = trueSM.computeGammaHgaga_WW();
    double gtW_SM = trueSM.computeGammaHgaga_tW();
    return (sqrt((computeKt() * computeKt() * gtt_SM
            + computeKW() * computeKW() * gWW_SM
            + computeKt() * computeKW() * gtW_SM) / (gtt_SM + gWW_SM + gtW_SM)));
}

double HiggsKvKfgen::computeKtau() const
{
    return Kl;
}

double HiggsKvKfgen::computeKc() const
{
    return Ku;
}

double HiggsKvKfgen::computeKt() const
{
    return Ku;
}

double HiggsKvKfgen::computeKb() const
{
    return Kd;
}


