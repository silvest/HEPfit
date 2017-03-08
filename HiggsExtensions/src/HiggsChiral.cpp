/* 
 * Copyright (C) 2017 HEPfit Collaboration
 *
 *
 * For the licensing terms see doc/COPYING.
 */

#include "HiggsChiral.h"

const std::string HiggsChiral::HChiralvars[NHChiralvars] = {
    "cv", "ct", "cb", "ctau", "cmu", "cg", "cga", "cZga"
};

HiggsChiral::HiggsChiral()
: NPbase()
{
    ModelParamMap.insert(std::pair<std::string, boost::reference_wrapper<const double> >("cv", boost::cref(cv)));
    ModelParamMap.insert(std::pair<std::string, boost::reference_wrapper<const double> >("ct", boost::cref(ct)));
    ModelParamMap.insert(std::pair<std::string, boost::reference_wrapper<const double> >("cb", boost::cref(cb)));
    ModelParamMap.insert(std::pair<std::string, boost::reference_wrapper<const double> >("ctau", boost::cref(ctau)));
    ModelParamMap.insert(std::pair<std::string, boost::reference_wrapper<const double> >("cmu", boost::cref(cmu)));
    ModelParamMap.insert(std::pair<std::string, boost::reference_wrapper<const double> >("cg", boost::cref(cg)));
    ModelParamMap.insert(std::pair<std::string, boost::reference_wrapper<const double> >("cga", boost::cref(cga)));
    ModelParamMap.insert(std::pair<std::string, boost::reference_wrapper<const double> >("cZga", boost::cref(cZga)));
}

void HiggsChiral::setParameter(const std::string name, const double& value)
{
    if (name.compare("cv") == 0)
        cv = value;
    else if (name.compare("ct") == 0)
        ct = value;
    else if (name.compare("cb") == 0)
        cb = value;
    else if (name.compare("ctau") == 0)
        ctau = value;
    else if (name.compare("cmu") == 0)
        cmu = value;
    else if (name.compare("cg") == 0)
        cg = value;
    else if (name.compare("cga") == 0)
        cga = value;
    else if (name.compare("cZga") == 0)
        cZga = value;
    else
        NPbase::setParameter(name, value);
}

bool HiggsChiral::CheckParameters(const std::map<std::string, double>& DPars)
{
    for (int i = 0; i < NHChiralvars; i++) {
        if (DPars.find(HChiralvars[i]) == DPars.end()) {
            std::cout << "missing mandatory HiggsChiral parameter " << HChiralvars[i] << std::endl;
            return false;
        }
    }
    return (NPbase::CheckParameters(DPars));
}

////////////////////////////////////////////////////////////////////////

double HiggsChiral::muggH(const double sqrt_s) const
{
    return (computeKg(sqrt_s) * computeKg(sqrt_s));
}

double HiggsChiral::muVBF(const double sqrt_s) const
{
    double sigmaWF_SM = trueSM.computeSigmaWF(sqrt_s);
    double sigmaZF_SM = trueSM.computeSigmaZF(sqrt_s);
    double sigmaZWF_SM = trueSM.computeSigmaZWF(sqrt_s);
    return (computeKW() * computeKW() * sigmaWF_SM
            + computeKZ() * computeKZ() * sigmaZF_SM
            + computeKW() * computeKZ() * sigmaZWF_SM)
            / (sigmaWF_SM + sigmaZF_SM + sigmaZWF_SM);
}

double HiggsChiral::mueeWBF(const double sqrt_s) const
{
    return (computeKW() * computeKW());
}

double HiggsChiral::muWH(const double sqrt_s) const
{
    return (computeKW() * computeKW());
}

double HiggsChiral::muZH(const double sqrt_s) const
{
    return (computeKZ() * computeKZ());
}

double HiggsChiral::mueeZH(const double sqrt_s) const
{
    return (computeKZ() * computeKZ());
}

double HiggsChiral::muVH(const double sqrt_s) const
{
    double sigmaWH_SM = trueSM.computeSigmaWH(sqrt_s);
    double sigmaZH_SM = trueSM.computeSigmaZH(sqrt_s);
    return ((computeKW() * computeKW() * sigmaWH_SM
            + computeKZ() * computeKZ() * sigmaZH_SM)
            / (sigmaWH_SM + sigmaZH_SM));
}

double HiggsChiral::muVBFpVH(const double sqrt_s) const
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

double HiggsChiral::muttH(const double sqrt_s) const
{
    return (computeKt() * computeKt());
}

double HiggsChiral::muggHpttH(const double sqrt_s) const
{
    double sigmaggH_SM = trueSM.computeSigmaggH(sqrt_s);
    double sigmattH_SM = trueSM.computeSigmattH(sqrt_s);

    double sigmaggH = muggH(sqrt_s) * sigmaggH_SM;
    double sigmattH = muttH(sqrt_s) * sigmattH_SM;

    return ((sigmaggH + sigmattH) / (sigmaggH_SM + sigmattH_SM));
} 

double HiggsChiral::mueettH(const double sqrt_s) const
{
    return (computeKt() * computeKt());
}

double HiggsChiral::BrHggRatio() const
{
    double gtt_SM = trueSM.computeGammaHgg_tt();
    double gbb_SM = trueSM.computeGammaHgg_bb();
    double gtb_SM = trueSM.computeGammaHgg_tb();
    double Ggg_over_GggSM = (computeKt() * computeKt() * gtt_SM
            + computeKb() * computeKb() * gbb_SM
            + computeKt() * computeKb() * gtb_SM) / (gtt_SM + gbb_SM + gtb_SM);
    return (Ggg_over_GggSM / computeGammaTotalRatio());
}

double HiggsChiral::BrHWWRatio() const
{
    return (computeKW() * computeKW() / computeGammaTotalRatio());
}

double HiggsChiral::BrHZZRatio() const
{
    return (computeKZ() * computeKZ() / computeGammaTotalRatio());
}

double HiggsChiral::BrHZgaRatio() const
{
    return (computeKZga() * computeKZga() / computeGammaTotalRatio());
}

double HiggsChiral::BrHgagaRatio() const
{
    return (computeKgaga() * computeKgaga() / computeGammaTotalRatio());
}

double HiggsChiral::BrHmumuRatio() const
{
    return (computeKmu() * computeKmu() / computeGammaTotalRatio());
}

double HiggsChiral::BrHtautauRatio() const
{
    return (computeKtau() * computeKtau() / computeGammaTotalRatio());
}

double HiggsChiral::BrHccRatio() const
{
    return (computeKc() * computeKc() / computeGammaTotalRatio());
}

double HiggsChiral::BrHbbRatio() const
{
    return (computeKb() * computeKb() / computeGammaTotalRatio());
}

double HiggsChiral::computeGammaTotalRatio() const
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
            + computeKb() * computeKb() * trueSM.computeBrHtobb())
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

double HiggsChiral::computeKg(const double sqrt_s) const
{
//    if (sqrt_s == 1.96)
//        return computeKt();
//
//    double sigmatt_SM = trueSM.computeSigmaggH_tt(sqrt_s);
//    double sigmabb_SM = trueSM.computeSigmaggH_bb(sqrt_s);
//    double sigmatb_SM = trueSM.computeSigmaggH_tb(sqrt_s);
//    return (sqrt((computeKt() * computeKt() * sigmatt_SM
//            + computeKb() * computeKb() * sigmabb_SM
//            + computeKt() * computeKb() * sigmatb_SM)
//            / (sigmatt_SM + sigmabb_SM + sigmatb_SM)));
    return cg;
}

double HiggsChiral::computeKW() const
{
    return cv;
}

double HiggsChiral::computeKZ() const
{
    return cv;
}

double HiggsChiral::computeKZga() const
{
    return cZga;
}

double HiggsChiral::computeKgaga() const
{
    return cga;
}

double HiggsChiral::computeKmu() const
{
    return cmu;
}

double HiggsChiral::computeKtau() const
{
    return ctau;
}

//double HiggsKvKfgen::computeKc() const
//{
//    return ct;
//}

double HiggsChiral::computeKt() const
{
    return ct;
}

double HiggsChiral::computeKb() const
{
    return cb;
}
