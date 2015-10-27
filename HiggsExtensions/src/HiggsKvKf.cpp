/* 
 * Copyright (C) 2014 HEPfit Collaboration
 * All rights reserved.
 *
 * For the licensing terms see doc/COPYING.
 */

#include "HiggsKvKf.h"

const std::string HiggsKvKf::HKvKfvars[NHKvKfvars] = {
    "Kv", "Kf", "BrHinv"
};

HiggsKvKf::HiggsKvKf()
: NPbase()
{
    ModelParamMap.insert(std::pair<std::string, boost::reference_wrapper<const double> >("Kv", boost::cref(Kv)));
    ModelParamMap.insert(std::pair<std::string, boost::reference_wrapper<const double> >("Kf", boost::cref(Kf)));
    ModelParamMap.insert(std::pair<std::string, boost::reference_wrapper<const double> >("BrHinv", boost::cref(BrHinv)));
}

void HiggsKvKf::setParameter(const std::string name, const double& value)
{
    if (name.compare("Kv") == 0)
        Kv = value;
    else if (name.compare("Kf") == 0)
        Kf = value;
    else if (name.compare("BrHinv") == 0)
        BrHinv = value;
    else
        NPbase::setParameter(name, value);
}

bool HiggsKvKf::CheckParameters(const std::map<std::string, double>& DPars)
{
    for (int i = 0; i < NHKvKfvars; i++) {
        if (DPars.find(HKvKfvars[i]) == DPars.end()) {
            std::cout << "missing mandatory HiggsKvKf parameter " << HKvKfvars[i] << std::endl;
            return false;
        }
    }
    return (NPbase::CheckParameters(DPars));
}

////////////////////////////////////////////////////////////////////////

double HiggsKvKf::obliqueS() const
{
    double Lambda;
    if (fabs(1.0 - Kv * Kv) < pow(10.0, -32.0))
        Lambda = pow(10.0, 19.0);
    else
        Lambda = 4.0 * M_PI * v() / sqrt(fabs(1.0 - Kv * Kv));

    return ( 1.0 / 12.0 / M_PI * (1.0 - Kv * Kv) * log(Lambda * Lambda / mHl / mHl));
}

double HiggsKvKf::obliqueT() const
{
    double Lambda;
    double cW2_SM = trueSM.cW2();
    if (fabs(1.0 - Kv * Kv) < pow(10.0, -32.0))
        Lambda = pow(10.0, 19.0);
    else
        Lambda = 4.0 * M_PI * v() / sqrt(fabs(1.0 - Kv * Kv));

    return ( -3.0 / 16.0 / M_PI / cW2_SM * (1.0 - Kv * Kv) * log(Lambda * Lambda / mHl / mHl));
}

double HiggsKvKf::obliqueU() const
{
    return 0.0;
}

////////////////////////////////////////////////////////////////////////

double HiggsKvKf::muggH(const double sqrt_s) const
{
    return (computeKg() * computeKg());
}

double HiggsKvKf::muVBF(const double sqrt_s) const
{
    double sigmaWF_SM = trueSM.computeSigmaWF(sqrt_s);
    double sigmaZF_SM = trueSM.computeSigmaZF(sqrt_s);
    double sigmaZWF_SM = trueSM.computeSigmaZWF(sqrt_s);
    return (computeKW() * computeKW() * sigmaWF_SM
            + computeKZ() * computeKZ() * sigmaZF_SM
            + computeKW() * computeKZ() * sigmaZWF_SM)
            / (sigmaWF_SM + sigmaZF_SM + sigmaZWF_SM);
}

double HiggsKvKf::muWH(const double sqrt_s) const
{
    return (computeKW() * computeKW());
}

double HiggsKvKf::muZH(const double sqrt_s) const
{
    return (computeKZ() * computeKZ());
}

double HiggsKvKf::muVH(const double sqrt_s) const
{
    double sigmaWH_SM = trueSM.computeSigmaWH(sqrt_s);
    double sigmaZH_SM = trueSM.computeSigmaZH(sqrt_s);
    return ((computeKW() * computeKW() * sigmaWH_SM
            + computeKZ() * computeKZ() * sigmaZH_SM)
            / (sigmaWH_SM + sigmaZH_SM));
}

double HiggsKvKf::muVBFpVH(const double sqrt_s) const
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

double HiggsKvKf::muttH(const double sqrt_s) const
{
    return (computeKt() * computeKt());
}

double HiggsKvKf::muggHpttH(const double sqrt_s) const
{
    double sigmaggH_SM = trueSM.computeSigmaggH(sqrt_s);
    double sigmattH_SM = trueSM.computeSigmattH(sqrt_s);

    double sigmaggH = muggH(sqrt_s) * sigmaggH_SM;
    double sigmattH = muttH(sqrt_s) * sigmattH_SM;

    return ((sigmaggH + sigmattH) / (sigmaggH_SM + sigmattH_SM));
} 

double HiggsKvKf::BrHggRatio() const
{
    return (computeKg() * computeKg() / computeGammaTotalRatio());
}

double HiggsKvKf::BrHWWRatio() const
{
    return (computeKW() * computeKW() / computeGammaTotalRatio());
}

double HiggsKvKf::BrHZZRatio() const
{
    return (computeKZ() * computeKZ() / computeGammaTotalRatio());
}

double HiggsKvKf::BrHZgaRatio() const
{
    return (computeKZga() * computeKZga() / computeGammaTotalRatio());
}

double HiggsKvKf::BrHgagaRatio() const
{
    return (computeKgaga() * computeKgaga() / computeGammaTotalRatio());
}

double HiggsKvKf::BrHtautauRatio() const
{
    return (computeKtau() * computeKtau() / computeGammaTotalRatio());
}

double HiggsKvKf::BrHccRatio() const
{
    return (computeKc() * computeKc() / computeGammaTotalRatio());
}

double HiggsKvKf::BrHbbRatio() const
{
    return (computeKb() * computeKb() / computeGammaTotalRatio());
}

double HiggsKvKf::computeGammaTotalRatio() const
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

double HiggsKvKf::computeKg() const
{
    return Kf;
}

double HiggsKvKf::computeKW() const
{
    return Kv;
}

double HiggsKvKf::computeKZ() const
{
    return Kv;
}

double HiggsKvKf::computeKZga() const
{
    double gtt_SM = trueSM.computeGammaHZga_tt();
    double gWW_SM = trueSM.computeGammaHZga_WW();
    double gtW_SM = trueSM.computeGammaHZga_tW();
    return (sqrt((computeKt() * computeKt() * gtt_SM
            + computeKW() * computeKW() * gWW_SM
            + computeKt() * computeKW() * gtW_SM) / (gtt_SM + gWW_SM + gtW_SM)));
}

double HiggsKvKf::computeKgaga() const
{
    double gtt_SM = trueSM.computeGammaHgaga_tt();
    double gWW_SM = trueSM.computeGammaHgaga_WW();
    double gtW_SM = trueSM.computeGammaHgaga_tW();
    return (sqrt((computeKt() * computeKt() * gtt_SM
            + computeKW() * computeKW() * gWW_SM
            + computeKt() * computeKW() * gtW_SM) / (gtt_SM + gWW_SM + gtW_SM)));
}

double HiggsKvKf::computeKtau() const
{
    return Kf;
}

double HiggsKvKf::computeKc() const
{
    return Kf;
}

double HiggsKvKf::computeKt() const
{
    return Kf;
}

double HiggsKvKf::computeKb() const
{
    return Kf;
}


