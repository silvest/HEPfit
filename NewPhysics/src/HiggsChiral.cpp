/* 
 * Copyright (C) 2017 HEPfit Collaboration
 *
 *
 * For the licensing terms see doc/COPYING.
 */

#include "HiggsChiral.h"

const std::string HiggsChiral::HChiralvars[NHChiralvars] = {
    "cv", "ct", "cb", "cc", "ctau", "cmu", "cg", "cga", "cZga", "obsZgaLimitATLAS13", "obsZgaLimitCMS13", "obsZgaLimitATLAS", "obsZgaLimitCMS", "expZgaLimitATLAS13", "expZgaLimitCMS13", "expZgaLimitATLAS", "expZgaLimitCMS"
};

HiggsChiral::HiggsChiral()
: NPbase()
{    
    FlagUniversalcf = false;
    FlagUniversalcvcf = false;
    
    ModelParamMap.insert(std::make_pair("cv", std::cref(cv)));
    ModelParamMap.insert(std::make_pair("ct", std::cref(ct)));
    ModelParamMap.insert(std::make_pair("cb", std::cref(cb)));
    ModelParamMap.insert(std::make_pair("cc", std::cref(cc)));
    ModelParamMap.insert(std::make_pair("ctau", std::cref(ctau)));
    ModelParamMap.insert(std::make_pair("cmu", std::cref(cmu)));
    ModelParamMap.insert(std::make_pair("cg", std::cref(cg)));
    ModelParamMap.insert(std::make_pair("cga", std::cref(cga)));
    ModelParamMap.insert(std::make_pair("cZga", std::cref(cZga)));
    ModelParamMap.insert(std::make_pair("obsZgaLimitATLAS13", std::cref(obsZgaLimitATLAS13)));
    ModelParamMap.insert(std::make_pair("obsZgaLimitCMS13", std::cref(obsZgaLimitCMS13)));
    ModelParamMap.insert(std::make_pair("obsZgaLimitATLAS", std::cref(obsZgaLimitATLAS)));
    ModelParamMap.insert(std::make_pair("obsZgaLimitCMS", std::cref(obsZgaLimitCMS)));
    ModelParamMap.insert(std::make_pair("expZgaLimitATLAS13", std::cref(expZgaLimitATLAS13)));
    ModelParamMap.insert(std::make_pair("expZgaLimitCMS13", std::cref(expZgaLimitCMS13)));
    ModelParamMap.insert(std::make_pair("expZgaLimitATLAS", std::cref(expZgaLimitATLAS)));
    ModelParamMap.insert(std::make_pair("expZgaLimitCMS", std::cref(expZgaLimitCMS)));
    
    loopComputed = false;
    cg_loop = 0.;
    cga_loop = 0.;
    cZga_loop = 0.;
}

bool HiggsChiral::PostUpdate()
{
    if (!NPbase::PostUpdate()) return (false);

//  Check first the flags that control the values of the parameters    

    if (FlagUniversalcf) {
//  Assign to all cf the value of ct
        cb = ct;
        cc = ct;
        ctau = ct;
        cmu = ct;
    }
    
    if (FlagUniversalcvcf) {
//  Assign to all cv and cf the value of ct
        cv = ct;
        cb = ct;
        cc = ct;
        ctau = ct;
        cmu = ct;
    }

//  Then the add the values of the loops to the cache
    
    loopComputed = false;
    
    cg_loop = computecg();
    cga_loop = computecgaga();
    cZga_loop = computecZga();
    
    loopComputed = true;

    return (true);
}

void HiggsChiral::setParameter(const std::string name, const double& value)
{
    if (name.compare("cv") == 0)
        cv = value;
    else if (name.compare("ct") == 0)
        ct = value;
    else if (name.compare("cb") == 0)
        cb = value;
    else if (name.compare("cc") == 0)
        cc = value;
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
    else if (name.compare("obsZgaLimitATLAS13") == 0)
        obsZgaLimitATLAS13 = value;
    else if (name.compare("obsZgaLimitCMS13") == 0)
        obsZgaLimitCMS13 = value;
    else if (name.compare("obsZgaLimitATLAS") == 0)
        obsZgaLimitATLAS = value;
    else if (name.compare("obsZgaLimitCMS") == 0)
        obsZgaLimitCMS = value;
    else if (name.compare("expZgaLimitATLAS13") == 0)
        expZgaLimitATLAS13 = value;
    else if (name.compare("expZgaLimitCMS13") == 0)
        expZgaLimitCMS13 = value;
    else if (name.compare("expZgaLimitATLAS") == 0)
        expZgaLimitATLAS = value;
    else if (name.compare("expZgaLimitCMS") == 0)
        expZgaLimitCMS = value;
    else
        NPbase::setParameter(name, value);
}

bool HiggsChiral::CheckParameters(const std::map<std::string, double>& DPars)
{
    for (int i = 0; i < NHChiralvars; i++) {
        if (DPars.find(HChiralvars[i]) == DPars.end()) {
            std::cout << "ERROR: missing mandatory HiggsChiral parameter " << HChiralvars[i] << std::endl;
            raiseMissingModelParameterCount();
            addMissingModelParameter(HChiralvars[i]);
        }
    }
    return (NPbase::CheckParameters(DPars));
}

bool HiggsChiral::setFlag(const std::string name, const bool value)
{
    bool res = false;
    if (name.compare("Universalcf") == 0) {
        FlagUniversalcf = value;
        res = true;
    } else if (name.compare("Universalcvcf") == 0) {
        FlagUniversalcvcf = value;
        res = true;
    } else
        res = NPbase::setFlag(name, value);

    return (res);
}

////////////////////////////////////////////////////////////////////////

double HiggsChiral::obliqueS() const
{
    double Lambda;
    if (fabs(1.0 - cv * cv) < pow(10.0, -32.0))
        Lambda = pow(10.0, 19.0);
    else
        Lambda = 4.0 * M_PI * v() / sqrt(fabs(1.0 - cv * cv));

    return ( 1.0 / 12.0 / M_PI * (1.0 - cv * cv) * log(Lambda * Lambda / mHl / mHl));
}

double HiggsChiral::obliqueT() const
{
    double Lambda;
    double cW2_SM = trueSM.cW2();
    if (fabs(1.0 - cv * cv) < pow(10.0, -32.0))
        Lambda = pow(10.0, 19.0);
    else
        Lambda = 4.0 * M_PI * v() / sqrt(fabs(1.0 - cv * cv));

    return ( -3.0 / 16.0 / M_PI / cW2_SM * (1.0 - cv * cv) * log(Lambda * Lambda / mHl / mHl));
}

double HiggsChiral::obliqueU() const
{
    return 0.0;
}

////////////////////////////////////////////////////////////////////////

double HiggsChiral::muggH(const double sqrt_s) const
{
//    return Gammagg() / (trueSM.computeBrHtogg()*trueSM.computeGammaHTotal());
      
    return computecg() * computecg();
}

double HiggsChiral::muVBF(const double sqrt_s) const
{
    return computecV() * computecV();
}

double HiggsChiral::muVBFgamma(const double sqrt_s) const
{
    return computecV() * computecV();
}

double HiggsChiral::mueeWBF(const double sqrt_s) const
{
    return computecV() * computecV();
}

double HiggsChiral::mueeWBFPol(const double sqrt_s, const double Pol_em, const double Pol_ep) const
{
    return computecV() * computecV();
}

double HiggsChiral::mueeHvv(const double sqrt_s) const
{
    return computecV() * computecV();
}

double HiggsChiral::mueeHvvPol(const double sqrt_s, const double Pol_em, const double Pol_ep) const
{
    return computecV() * computecV();
}

double HiggsChiral::mueeZBF(const double sqrt_s) const
{
    return computecV() * computecV();
}

double HiggsChiral::mueeZBFPol(const double sqrt_s, const double Pol_em, const double Pol_ep) const
{
    return computecV() * computecV();
}

double HiggsChiral::muWH(const double sqrt_s) const
{
    return computecV() * computecV();
}

double HiggsChiral::muZH(const double sqrt_s) const
{
    return computecV() * computecV();
}

double HiggsChiral::mueeZH(const double sqrt_s) const
{
    return computecV() * computecV();
}

double HiggsChiral::mueeZllH(const double sqrt_s) const
{
    return computecV() * computecV();
}

double HiggsChiral::mueeZqqH(const double sqrt_s) const
{
    return computecV() * computecV();
}

double HiggsChiral::mueeZHPol(const double sqrt_s, const double Pol_em, const double Pol_ep) const
{
    return computecV() * computecV();
}

double HiggsChiral::mueeZllHPol(const double sqrt_s, const double Pol_em, const double Pol_ep) const
{
    return computecV() * computecV();
}

double HiggsChiral::mueeZqqHPol(const double sqrt_s, const double Pol_em, const double Pol_ep) const
{
    return computecV() * computecV();
}

double HiggsChiral::muVH(const double sqrt_s) const
{
    return computecV() * computecV();
}

double HiggsChiral::muVBFpVH(const double sqrt_s) const
{
    return computecV() * computecV();
}

double HiggsChiral::muttH(const double sqrt_s) const
{
    return computect() * computect();
}

double HiggsChiral::mutHq(const double sqrt_s) const
{
    double mu = 1.0;
    double Kt = computect();
    double Kw = computecV();
    
//  From ATLAS parametrization of the K-framework (at 14 TeV?)     
    mu = 2.63 * Kt * Kt + 3.58 * Kw * Kw - 5.21 * Kt * Kw;
    
    return mu;
}

double HiggsChiral::muggHpttH(const double sqrt_s) const
{
    double sigmaggH_SM = trueSM.computeSigmaggH(sqrt_s);
    double sigmattH_SM = trueSM.computeSigmattH(sqrt_s);

    double sigmaggH = muggH(sqrt_s) * sigmaggH_SM;
    double sigmattH = muttH(sqrt_s) * sigmattH_SM;

    return (sigmaggH + sigmattH) / (sigmaggH_SM + sigmattH_SM);
} 

double HiggsChiral::mueettH(const double sqrt_s) const
{
    return computect() * computect();
}

double HiggsChiral::mueettHPol(const double sqrt_s, const double Pol_em, const double Pol_ep) const
{
    return computect() * computect();
}

double HiggsChiral::UpperLimitZgammaA13(const double sqrt_s) const
{
    double sigmaggH_SM = trueSM.computeSigmaggH(sqrt_s);
    double sigmattH_SM = trueSM.computeSigmattH(sqrt_s);
    double sigmaWH_SM = trueSM.computeSigmaWH(sqrt_s);
    double sigmaZH_SM = trueSM.computeSigmaZH(sqrt_s);
    double sigmaWF_SM = trueSM.computeSigmaWF(sqrt_s);
    double sigmaZF_SM = trueSM.computeSigmaZF(sqrt_s);
    double sigmaZWF_SM = trueSM.computeSigmaZWF(sqrt_s);    
    double sigmaVBF_SM = sigmaWF_SM + sigmaZF_SM + sigmaZWF_SM;

    double sigmaggH = muggH(sqrt_s) * sigmaggH_SM;
    double sigmattH = muttH(sqrt_s) * sigmattH_SM;
    double sigmaWH = muWH(sqrt_s) * sigmaWH_SM;
    double sigmaZH = muZH(sqrt_s) * sigmaZH_SM;
    double sigmaVBF = muVBF(sqrt_s) * sigmaVBF_SM;

    double muppH=((sigmaggH + sigmattH + sigmaWH + sigmaZH + sigmaVBF) / (sigmaggH_SM + sigmattH_SM + sigmaWH_SM + sigmaZH_SM + sigmaVBF_SM));
//    double BrHZgaRatio=;
    double nftos=1.95996398454;
        
    return (1.0+(muppH*BrHZgaRatio()-getobsZgaLimitATLAS13())/getexpZgaLimitATLAS13() ) * nftos;
}

double HiggsChiral::UpperLimitZgammaC13(const double sqrt_s) const
{
    double sigmaggH_SM = trueSM.computeSigmaggH(sqrt_s);
    double sigmattH_SM = trueSM.computeSigmattH(sqrt_s);
    double sigmaWH_SM = trueSM.computeSigmaWH(sqrt_s);
    double sigmaZH_SM = trueSM.computeSigmaZH(sqrt_s);
    double sigmaWF_SM = trueSM.computeSigmaWF(sqrt_s);
    double sigmaZF_SM = trueSM.computeSigmaZF(sqrt_s);
    double sigmaZWF_SM = trueSM.computeSigmaZWF(sqrt_s);    
    double sigmaVBF_SM = sigmaWF_SM + sigmaZF_SM + sigmaZWF_SM;

    double sigmaggH = muggH(sqrt_s) * sigmaggH_SM;
    double sigmattH = muttH(sqrt_s) * sigmattH_SM;
    double sigmaWH = muWH(sqrt_s) * sigmaWH_SM;
    double sigmaZH = muZH(sqrt_s) * sigmaZH_SM;
    double sigmaVBF = muVBF(sqrt_s) * sigmaVBF_SM;

    double muppH=((sigmaggH + sigmattH + sigmaWH + sigmaZH + sigmaVBF) / (sigmaggH_SM + sigmattH_SM + sigmaWH_SM + sigmaZH_SM + sigmaVBF_SM));
//    double BrHZgaRatio=;
    double nftos=1.95996398454;
        
    return (1.0+(muppH*BrHZgaRatio()-getobsZgaLimitCMS13())/getexpZgaLimitCMS13() ) * nftos;
}

double HiggsChiral::UpperLimitZgammaA(const double sqrt_s) const
{
    double sigmaggH_SM = trueSM.computeSigmaggH(sqrt_s);
    double sigmattH_SM = trueSM.computeSigmattH(sqrt_s);
    double sigmaWH_SM = trueSM.computeSigmaWH(sqrt_s);
    double sigmaZH_SM = trueSM.computeSigmaZH(sqrt_s);
    double sigmaWF_SM = trueSM.computeSigmaWF(sqrt_s);
    double sigmaZF_SM = trueSM.computeSigmaZF(sqrt_s);
    double sigmaZWF_SM = trueSM.computeSigmaZWF(sqrt_s);    
    double sigmaVBF_SM = sigmaWF_SM + sigmaZF_SM + sigmaZWF_SM;

    double sigmaggH = muggH(sqrt_s) * sigmaggH_SM;
    double sigmattH = muttH(sqrt_s) * sigmattH_SM;
    double sigmaWH = muWH(sqrt_s) * sigmaWH_SM;
    double sigmaZH = muZH(sqrt_s) * sigmaZH_SM;
    double sigmaVBF = muVBF(sqrt_s) * sigmaVBF_SM;

    double muppH=((sigmaggH + sigmattH + sigmaWH + sigmaZH + sigmaVBF) / (sigmaggH_SM + sigmattH_SM + sigmaWH_SM + sigmaZH_SM + sigmaVBF_SM));
//    double BrHZgaRatio=;
    double nftos=1.95996398454;
        
    return (1.0+(muppH*BrHZgaRatio()-getobsZgaLimitATLAS())/getexpZgaLimitATLAS() ) * nftos;
}

double HiggsChiral::UpperLimitZgammaC(const double sqrt_s) const
{
    double sigmaggH_SM = trueSM.computeSigmaggH(sqrt_s);
    double sigmattH_SM = trueSM.computeSigmattH(sqrt_s);
    double sigmaWH_SM = trueSM.computeSigmaWH(sqrt_s);
    double sigmaZH_SM = trueSM.computeSigmaZH(sqrt_s);
    double sigmaWF_SM = trueSM.computeSigmaWF(sqrt_s);
    double sigmaZF_SM = trueSM.computeSigmaZF(sqrt_s);
    double sigmaZWF_SM = trueSM.computeSigmaZWF(sqrt_s);    
    double sigmaVBF_SM = sigmaWF_SM + sigmaZF_SM + sigmaZWF_SM;

    double sigmaggH = muggH(sqrt_s) * sigmaggH_SM;
    double sigmattH = muttH(sqrt_s) * sigmattH_SM;
    double sigmaWH = muWH(sqrt_s) * sigmaWH_SM;
    double sigmaZH = muZH(sqrt_s) * sigmaZH_SM;
    double sigmaVBF = muVBF(sqrt_s) * sigmaVBF_SM;

    double muppH=((sigmaggH + sigmattH + sigmaWH + sigmaZH + sigmaVBF) / (sigmaggH_SM + sigmattH_SM + sigmaWH_SM + sigmaZH_SM + sigmaVBF_SM));
//    double BrHZgaRatio=;
    double nftos=1.95996398454;
        
    return (1.0+(muppH*BrHZgaRatio()-getobsZgaLimitCMS())/getexpZgaLimitCMS() ) * nftos;
}

double HiggsChiral::cgplusct() const
{
    return cg + ct;
}

double HiggsChiral::cgaplusct() const
{
    return cga + ct;
}

double HiggsChiral::cgminuscga() const
{
    return cg - cga;
}

double HiggsChiral::cVpluscb() const
{
    return cv + cb;
}

double HiggsChiral::cVplusctau() const
{
    return cv + ctau;
}

double HiggsChiral::cbminuscc() const
{
    return cb - cc;
}

double HiggsChiral::cbminusctau() const
{
    return cb - ctau;
}

double HiggsChiral::ccminusctau() const
{
    return cc - ctau;
}

double HiggsChiral::Gammagg() const
{
    double Mt=trueSM.getQuarks(QCD::TOP).getMass();
    double Mb=trueSM.getQuarks(QCD::BOTTOM).getMass();
    double TAUt=4.0*Mt*Mt/(mHl*mHl);
    double TAUb=4.0*Mb*Mb/(mHl*mHl);

    double factor = GF*AlsMz*AlsMz*mHl*mHl*mHl/(sqrt(2.0)*4.0*M_PI*M_PI*M_PI);
      
    double ksoftNLO = 1.0 + (AlsMz/M_PI) * ( 73.0/4.0 -7.0*5.0/6.0);

    return factor * ksoftNLO * ( ( ct * 0.5*TAUt*(1.0+(1.0-TAUt)*f_func(TAUt))
                     +cb * 0.5*TAUb*(1.0+(1.0-TAUb)*f_func(TAUb)) ) * (1.0+11.0*AlsMz/(4.0*M_PI))
                     +cg * 0.5 ).abs2();
}

double HiggsChiral::GammaWW() const
{
    return computecV() * computecV() * trueSM.computeBrHtoWW() * trueSM.computeGammaHTotal();
}

double HiggsChiral::GammaZZ() const
{
    return computecV() * computecV() * trueSM.computeBrHtoZZ() * trueSM.computeGammaHTotal();
}

double HiggsChiral::GammaZga() const
{
    double Mt=trueSM.getQuarks(QCD::TOP).getMass();
    double Mb=trueSM.getQuarks(QCD::BOTTOM).getMass();
    double Mtau=trueSM.getLeptons(StandardModel::TAU).getMass();
    double MW=trueSM.Mw_tree();
    double cW2=trueSM.c02();
    double sW2=1.0-cW2;
    double TAUt=4.0*Mt*Mt/(mHl*mHl);
    double TAUb=4.0*Mb*Mb/(mHl*mHl);
    double TAUtau=4.0*Mtau*Mtau/(mHl*mHl);
    double TAUw=4.0*MW*MW/(mHl*mHl);
    double LAMt=4.0*Mt*Mt/(Mz*Mz);
    double LAMb=4.0*Mb*Mb/(Mz*Mz);
    double LAMtau=4.0*Mtau*Mtau/(Mz*Mz);
    double LAMw=4.0*MW*MW/(Mz*Mz);

    double factor = GF*ale*ale/(sqrt(2.0)*64.0*M_PI*M_PI*M_PI)*pow(mHl*(1.0-Mz*Mz/(mHl*mHl)),3);

    return factor * ((-ct * 4.0*(0.5-4.0/3.0*sW2)*(Int1(TAUt,LAMt)-Int2(TAUt,LAMt)) * (1.0-AlsMz/M_PI)
                      +cb * 2.0*(-0.5+2.0/3.0*sW2)*(Int1(TAUb,LAMb)-Int2(TAUb,LAMb))
                      +ctau * 2.0*(-0.5+2.0*sW2)*(Int1(TAUtau,LAMtau)-Int2(TAUtau,LAMtau)) )/sqrt(sW2*cW2)
                     -cv * sqrt(cW2/sW2)*(4.0*(3.0-sW2/cW2)*Int2(TAUw,LAMw)
                            +((1.0+2.0/TAUw)*sW2/cW2-(5.0+2.0/TAUw))*Int1(TAUw,LAMw))
                     +cZga).abs2(); /*CHANGE THE cZga COEFFICIENT!*/
}

double HiggsChiral::Gammagaga() const
{
    double Mt=trueSM.getQuarks(QCD::TOP).getMass();
    double Mb=trueSM.getQuarks(QCD::BOTTOM).getMass();
    double Mtau=trueSM.getLeptons(StandardModel::TAU).getMass();
    double MW=trueSM.Mw_tree();
    double TAUt=4.0*Mt*Mt/(mHl*mHl);
    double TAUb=4.0*Mb*Mb/(mHl*mHl);
    double TAUtau=4.0*Mtau*Mtau/(mHl*mHl);
    double TAUw=4.0*MW*MW/(mHl*mHl);

    double factor = GF*ale*ale*mHl*mHl*mHl/(sqrt(2.0)*128.0*M_PI*M_PI*M_PI);

    return factor * ( ct * (8./3.)*TAUt*(1.+(1.-TAUt)*f_func(TAUt)) * (1.0-AlsMz/M_PI)
                     +cb * (2./3.)*TAUb*(1.+(1.-TAUb)*f_func(TAUb))
                     +ctau * 2.0*TAUtau*(1.+(1.-TAUtau)*f_func(TAUtau))
                     -cv * (2.0+3.0*TAUw+3.0*TAUw*(2.0-TAUw)*f_func(TAUw))
                     +cga * 2.0 ).abs2();
}

double HiggsChiral::Gammamumu() const
{
    return computecmu() * computecmu() * trueSM.computeBrHtomumu() * trueSM.computeGammaHTotal();
}

double HiggsChiral::Gammatautau() const
{
    return computectau() * computectau() * trueSM.computeBrHtotautau() * trueSM.computeGammaHTotal();
}

double HiggsChiral::Gammacc() const
{
    return computecc() * computecc() * trueSM.computeBrHtocc() * trueSM.computeGammaHTotal();
}

double HiggsChiral::Gammabb() const
{
    return computecb() * computecb() * trueSM.computeBrHtobb() * trueSM.computeGammaHTotal();
}

double HiggsChiral::GammaTotal() const
{
    return Gammagg() + GammaWW() + GammaZZ() + GammaZga() + Gammagaga() 
            + Gammamumu() + Gammatautau() + Gammacc() + Gammabb();
}

double HiggsChiral::BrHggRatio() const
{
//    return Gammagg() / GammaTotal() / trueSM.computeBrHtogg();
    
    return (computecg() * computecg() / computeGammaTotalRatio());
}

double HiggsChiral::BrHWWRatio() const
{
//    return GammaWW() / GammaTotal() / trueSM.computeBrHtoWW();

    return (computecV() * computecV() / computeGammaTotalRatio());
}

double HiggsChiral::BrHWW2l2vRatio() const
{
//    return GammaWW() / GammaTotal() / trueSM.computeBrHtoWW();

    return BrHWWRatio();
}

double HiggsChiral::BrHZZRatio() const
{
//    return GammaZZ() / GammaTotal() / trueSM.computeBrHtoZZ();
    
    return (computecV() * computecV() / computeGammaTotalRatio());
}

double HiggsChiral::BrHZZ4lRatio() const
{
//    return GammaZZ() / GammaTotal() / trueSM.computeBrHtoZZ();
    
    return BrHZZRatio();
}

double HiggsChiral::BrHZZ4eRatio() const
{
//    return GammaZZ() / GammaTotal() / trueSM.computeBrHtoZZ();
    
    return BrHZZRatio();
}

double HiggsChiral::BrHZZ2e2muRatio() const
{
//    return GammaZZ() / GammaTotal() / trueSM.computeBrHtoZZ();
    
    return BrHZZRatio();
}

double HiggsChiral::BrHZZ4muRatio() const
{
//    return GammaZZ() / GammaTotal() / trueSM.computeBrHtoZZ();
    
    return BrHZZRatio();
}

double HiggsChiral::BrHZZ4vRatio() const
{
//    return GammaZZ() / GammaTotal() / trueSM.computeBrHtoZZ();
    
    return BrHZZRatio();
}

double HiggsChiral::BrHVVRatio() const
{
    double BrZZSM = trueSM.computeBrHtoZZ(), BrWWSM = trueSM.computeBrHtoWW();

    return (BrZZSM * BrHZZRatio() + BrWWSM * BrHWWRatio() )/(BrZZSM + BrWWSM);
}

double HiggsChiral::BrHZgaRatio() const
{
//    return GammaZga() / GammaTotal() / trueSM.computeBrHtoZga();
    
    return (computecZga() * computecZga() / computeGammaTotalRatio());
}

double HiggsChiral::BrHZgallRatio() const
{
//    return GammaZga() / GammaTotal() / trueSM.computeBrHtoZga();
    
    return BrHZgaRatio();
}

double HiggsChiral::BrHZgaeeRatio() const
{
//    return GammaZga() / GammaTotal() / trueSM.computeBrHtoZga();
    
    return BrHZgaRatio();
}

double HiggsChiral::BrHZgamumuRatio() const
{
//    return GammaZga() / GammaTotal() / trueSM.computeBrHtoZga();
    
    return BrHZgaRatio();
}

double HiggsChiral::BrHgagaRatio() const
{
//    return Gammagaga() / GammaTotal() / trueSM.computeBrHtogaga();
    
    return (computecgaga() * computecgaga() / computeGammaTotalRatio());
}

double HiggsChiral::BrHmumuRatio() const
{
//    return Gammamumu() / GammaTotal() / trueSM.computeBrHtomumu();
    
    return (computecmu() * computecmu() / computeGammaTotalRatio());
}

double HiggsChiral::BrHtautauRatio() const
{
//    return Gammatautau() / GammaTotal() / trueSM.computeBrHtotautau();
    
    return (computectau() * computectau() / computeGammaTotalRatio());
}

double HiggsChiral::BrHccRatio() const
{
//    return Gammacc() / GammaTotal() / trueSM.computeBrHtocc();
    
    return (computecc() * computecc() / computeGammaTotalRatio());
}

double HiggsChiral::BrHbbRatio() const
{
//    return Gammabb() / GammaTotal() / trueSM.computeBrHtobb();
    
    return (computecb() * computecb() / computeGammaTotalRatio());
}

    // General 4 fermion final states 

double HiggsChiral::BrH2L2LRatio() const
{    
    return BrHZZRatio();
}
    
double HiggsChiral::BrH2e2muRatio() const
{    
    return BrHZZRatio();
}    

double HiggsChiral::BrH2v2vRatio() const
{    
    return BrHZZRatio();
}

double HiggsChiral::BrH2L2vRatio() const
{    
    return BrHZZRatio();
}
    
double HiggsChiral::BrH2L2v2Ratio() const
{    
    return BrHZZRatio();
}
    
double HiggsChiral::BrH2e2vRatio() const
{    
    return BrHZZRatio();
}
    
double HiggsChiral::BrH2mu2vRatio() const
{    
    return BrHZZRatio();
}
    
double HiggsChiral::BrH2u2uRatio() const
{    
    return BrHZZRatio();
}
    
double HiggsChiral::BrH2d2dRatio() const
{    
    return BrHZZRatio();
}
    
double HiggsChiral::BrH2u2dRatio() const
{    
    return BrHZZRatio();
}
    
double HiggsChiral::BrH2L2uRatio() const
{    
    return BrHZZRatio();
}
    
double HiggsChiral::BrH2L2dRatio() const
{    
    return BrHZZRatio();
}
    
double HiggsChiral::BrH2v2uRatio() const
{    
    return BrHZZRatio();
}
    
double HiggsChiral::BrH2v2dRatio() const
{    
    return BrHZZRatio();
}
    
double HiggsChiral::BrH4LRatio() const
{    
    return BrHZZRatio();
}
    
double HiggsChiral::BrH4L2Ratio() const
{    
    return BrHZZRatio();
}
    
double HiggsChiral::BrH4eRatio() const
{    
    return BrHZZRatio();
}
    
double HiggsChiral::BrH4muRatio() const
{    
    return BrHZZRatio();
}
    
double HiggsChiral::BrH4vRatio() const
{    
    return BrHZZRatio();
}
    
double HiggsChiral::BrH4uRatio() const
{    
    return BrHZZRatio();
}
    
double HiggsChiral::BrH4dRatio() const
{    
    return BrHZZRatio();
}
    
double HiggsChiral::BrHLvvLRatio() const
{    
    return BrHWWRatio();
}
    
double HiggsChiral::BrHevmuvRatio() const
{    
    return BrHWWRatio();
}
    
double HiggsChiral::BrHudduRatio() const
{    
    return BrHWWRatio();
}
    
double HiggsChiral::BrHLvudRatio() const
{    
    return BrHWWRatio();
}
    
double HiggsChiral::BrH2udRatio() const
{    
    // Approximate by pure CC
    return BrHWWRatio();
}
    
double HiggsChiral::BrH2LvRatio() const
{    
    // Approximate by pure CC
    return BrHWWRatio();
}   
    
double HiggsChiral::BrH2Lv2Ratio() const
{    
    // Approximate by pure CC
    return BrHWWRatio();
}
    
double HiggsChiral::BrH2evRatio() const
{    
    // Approximate by pure CC
    return BrHWWRatio();
}    
    
double HiggsChiral::BrH2muvRatio() const
{    
    // Approximate by pure CC
    return BrHWWRatio();
}    
    
double HiggsChiral::BrH4fRatio() const
{
    double Br = 0.0;
    
    // SM decay widths (from MG simulations)
    double wH2L2LSM=0.65682e-06, wH2v2vSM=0.28126e-05, wH2L2vSM=0.27224e-05;
    double wH2u2uSM=0.22500e-05, wH2d2dSM=0.11906e-04, wH2u2dSM=0.12361e-04;
    double wH2L2uSM=0.45029e-05, wH2L2dSM=0.85830e-05, wH2v2uSM=0.93233e-05;
    double wH2v2dSM=0.17794e-04, wH4LSM=0.33973e-06, wH4vSM=0.16884e-05;
    double wH4uSM=0.23669e-05, wH4dSM=0.60254e-05;
    double wHLvvLSM=0.58098e-04, wHudduSM=0.13384e-03, wHLvudSM=0.34149e-03 ;
    double wH2udSM=0.13711e-03, wH2LvSM=0.27557e-04;

    // Sum: SM
    double wH4fSM=wH2L2LSM+wH2v2vSM+wH2L2vSM+wH2u2uSM+wH2d2dSM+wH2u2dSM+
    wH2L2uSM+wH2L2dSM+wH2v2uSM+wH2v2dSM+wH4LSM+wH4vSM+wH4uSM+wH4dSM+wHLvvLSM+wHudduSM+
    wHLvudSM+wH2udSM+wH2LvSM;

    // Sum: SM+NP
    double wH4f= (computecV() * computecV()) * (wH2L2LSM+wH2v2vSM+wH2L2vSM+wH2u2uSM+wH2d2dSM+wH2u2dSM
            + wH2L2uSM+wH2L2dSM+wH2v2uSM+wH2v2dSM+wH4LSM+wH4vSM+wH4uSM+wH4dSM)
            + (computecV() * computecV()) * (wHLvvLSM+wHudduSM
            + wHLvudSM+wH2udSM+wH2LvSM);
        
    Br = (wH4f/wH4fSM) / computeGammaTotalRatio();
    
    return Br;
}
    
    // DECAYS INVOLVING ONLY ELECTRONS, MUONS OR NEUTRINOS IN THE FINAL STATES 

double HiggsChiral::BrH4lRatio() const
{    
    return BrHZZRatio();
}
    
double HiggsChiral::BrH2l2vRatio() const
{  
    // Approximate by pure CC
    return BrHWWRatio();
}

double HiggsChiral::BrHtoinvRatio() const
{
//    No extra decays so it scales as the H->ZZ
    
    return BrHZZRatio();
}

double HiggsChiral::muggHgaga(const double sqrt_s) const
{
    return computecg() * computecg() * computecgaga() * computecgaga() / computeGammaTotalRatio();
}

double HiggsChiral::muggHgagaInt(const double sqrt_s) const
{
    double muNWA, GammaRatio, fki;
    double kt,kb,kc,ks,ku,kd,kta,km,ke;
    double kw;
    double kg,kgamma;
    
    muNWA = muggH(sqrt_s)*BrHgagaRatio();
    
    kt = computect();
    kb = computecb();
    kc = computecc();
    ks = 1.0;
    ku = 1.0;
    kd = 1.0;
    kta = computectau();
    km = computecmu();
    ke = 1.0;
    
    kw = computecV();
    
    kg = (3.0/2.0) * computecg();
    kgamma = (3.0/2.0) * computecgaga();
    
    GammaRatio = computeGammaTotalRatio();

    fki = 0.000802422 *kb*kb + 0.000312884 *kb*kc + 0.0000182107 *kc*kc +
            5.94769e-9 *kb*kd + 9.62554e-10 *kc*kd + 6.2785e-15 *kd*kd + 
            5.53251e-10 *kb*ke + 3.51863e-11 *kc*ke + 1.09243e-15 *kd*ke - 
            0.00905016 *kb*kg - 0.00190706 *kc*kg - 5.97591e-9 *kd*kg - 
            6.72288e-10 *ke*kg - 0.0271505 *kb*kgamma - 0.00143029 *kc*kgamma - 
            1.79277e-8 *kd*kgamma - 0.174392 *kg*kgamma + 8.97565e-6 *kb*km + 
            6.21013e-7 *kc*km + 2.37616e-11 *kd*km - 0.0000460022 *kg*km + 
            4.93348e-6 *kb*ks + 8.51176e-7 *kc*ks + 1.29558e-11 *kd*ks + 
            1.16267e-12 *ke*ks - 0.0000123381 *kg*ks - 0.0000370143 *kgamma*ks + 
            2.22544e-8 *km*ks + 6.08665e-9 *ks*ks - 0.0467672 *kb*kt - 
            0.00394193 *kc*kt - 3.08808e-8 *kd*kt - 6.94817e-10 *ke*kt - 
            0.240315 *kg*kt - 0.180236 *kgamma*kt - 0.0000475437 *km*kt - 
            0.0000637578 *ks*kt - 0.248368 *kt*kt + 0.00100168 *kb*kta + 
            0.0000759092 *kc*kta + 3.44671e-9 *kd*kta - 0.00975386 *kg*kta + 
            2.93009e-6 *ks*kta - 0.0100807 *kt*kta + 5.30126e-8 *kb*ku + 
            5.54256e-9 *kc*ku + 1.15815e-13 *kd*ku + 4.05206e-15 *ke*ku - 
            1.03323e-7 *kg*ku - 7.74926e-8 *kgamma*ku + 
            8.62762e-11 *km*ku + 1.17664e-10 *ks*ku - 2.13572e-7 *kt*ku + 
            1.2332e-8 *kta*ku + 3.40922e-13 *ku*ku + 0.169912 *kb*kw + 
            0.00895098 *kc*kw + 1.12194e-7 *kd*kw + 1.09137 *kg*kw + 
            0.000231641 *ks*kw + 1.12795 *kt*kw + 4.8496e-7 *ku*kw;
    
    return (muNWA - 0.022 *GammaRatio * fki)/0.978;

}

double HiggsChiral::muVBFHgaga(const double sqrt_s) const
{
    return computecV() * computecV() * computecgaga() * computecgaga() / computeGammaTotalRatio();
}

double HiggsChiral::muZHgaga(const double sqrt_s) const
{
    return computecV() * computecV() * computecgaga() * computecgaga() / computeGammaTotalRatio();
}

double HiggsChiral::muWHgaga(const double sqrt_s) const
{
    return computecV() * computecV() * computecgaga() * computecgaga() / computeGammaTotalRatio();
}

double HiggsChiral::muVHgaga(const double sqrt_s) const
{
    return computecV() * computecV() * computecgaga() * computecgaga() / computeGammaTotalRatio();
}

double HiggsChiral::muttHgaga(const double sqrt_s) const
{
    return computect() * computect() * computecgaga() * computecgaga() / computeGammaTotalRatio();
}

double HiggsChiral::muggHZga(const double sqrt_s) const
{
    return computecg() * computecg() * computecZga() * computecZga() / computeGammaTotalRatio();
}

double HiggsChiral::muVBFHZga(const double sqrt_s) const
{
    return computecV() * computecV() * computecZga() * computecZga() / computeGammaTotalRatio();
}

double HiggsChiral::muZHZga(const double sqrt_s) const
{
    return computecV() * computecV() * computecZga() * computecZga() / computeGammaTotalRatio();
}

double HiggsChiral::muWHZga(const double sqrt_s) const
{
    return computecV() * computecV() * computecZga() * computecZga() / computeGammaTotalRatio();
}

double HiggsChiral::muVHZga(const double sqrt_s) const
{
    return computecV() * computecV() * computecZga() * computecZga() / computeGammaTotalRatio();
}

double HiggsChiral::muttHZga(const double sqrt_s) const
{
    return computect() * computect() * computecZga() * computecZga() / computeGammaTotalRatio();
}

double HiggsChiral::muggHZZ(const double sqrt_s) const
{
    return computecg() * computecg() * computecV() * computecV() / computeGammaTotalRatio();
}

double HiggsChiral::muVBFHZZ(const double sqrt_s) const
{
    return computecV() * computecV() * computecV() * computecV() / computeGammaTotalRatio();
}

double HiggsChiral::muZHZZ(const double sqrt_s) const
{
    return computecV() * computecV() * computecV() * computecV() / computeGammaTotalRatio();
}

double HiggsChiral::muWHZZ(const double sqrt_s) const
{
    return computecV() * computecV() * computecV() * computecV() / computeGammaTotalRatio();
}

double HiggsChiral::muVHZZ(const double sqrt_s) const
{
    return computecV() * computecV() * computecV() * computecV() / computeGammaTotalRatio();
}

double HiggsChiral::muttHZZ(const double sqrt_s) const
{
    return computect() * computect() * computecV() * computecV() / computeGammaTotalRatio();
}

double HiggsChiral::muggHZZ4l(const double sqrt_s) const
{
    return computecg() * computecg() * computecV() * computecV() / computeGammaTotalRatio();
}

double HiggsChiral::muVBFHZZ4l(const double sqrt_s) const
{
    return computecV() * computecV() * computecV() * computecV() / computeGammaTotalRatio();
}

double HiggsChiral::muZHZZ4l(const double sqrt_s) const
{
    return computecV() * computecV() * computecV() * computecV() / computeGammaTotalRatio();
}

double HiggsChiral::muWHZZ4l(const double sqrt_s) const
{
    return computecV() * computecV() * computecV() * computecV() / computeGammaTotalRatio();
}

double HiggsChiral::muVHZZ4l(const double sqrt_s) const
{
    return computecV() * computecV() * computecV() * computecV() / computeGammaTotalRatio();
}

double HiggsChiral::muttHZZ4l(const double sqrt_s) const
{
    return computect() * computect() * computecV() * computecV() / computeGammaTotalRatio();
}

double HiggsChiral::muggHWW(const double sqrt_s) const
{
    return computecg() * computecg() * computecV() * computecV() / computeGammaTotalRatio();
}

double HiggsChiral::muVBFHWW(const double sqrt_s) const
{
    return computecV() * computecV() * computecV() * computecV() / computeGammaTotalRatio();
}

double HiggsChiral::muZHWW(const double sqrt_s) const
{
    return computecV() * computecV() * computecV() * computecV() / computeGammaTotalRatio();
}

double HiggsChiral::muWHWW(const double sqrt_s) const
{
    return computecV() * computecV() * computecV() * computecV() / computeGammaTotalRatio();
}

double HiggsChiral::muVHWW(const double sqrt_s) const
{
    return computecV() * computecV() * computecV() * computecV() / computeGammaTotalRatio();
}

double HiggsChiral::muttHWW(const double sqrt_s) const
{
    return computect() * computect() * computecV() * computecV() / computeGammaTotalRatio();
}

double HiggsChiral::muggHWW2l2v(const double sqrt_s) const
{
    return computecg() * computecg() * computecV() * computecV() / computeGammaTotalRatio();
}

double HiggsChiral::muVBFHWW2l2v(const double sqrt_s) const
{
    return computecV() * computecV() * computecV() * computecV() / computeGammaTotalRatio();
}

double HiggsChiral::muZHWW2l2v(const double sqrt_s) const
{
    return computecV() * computecV() * computecV() * computecV() / computeGammaTotalRatio();
}

double HiggsChiral::muWHWW2l2v(const double sqrt_s) const
{
    return computecV() * computecV() * computecV() * computecV() / computeGammaTotalRatio();
}

double HiggsChiral::muVHWW2l2v(const double sqrt_s) const
{
    return computecV() * computecV() * computecV() * computecV() / computeGammaTotalRatio();
}

double HiggsChiral::muttHWW2l2v(const double sqrt_s) const
{
    return computect() * computect() * computecV() * computecV() / computeGammaTotalRatio();
}

double HiggsChiral::muggHmumu(const double sqrt_s) const
{
    return computecg() * computecg() * computecmu() * computecmu() / computeGammaTotalRatio();
}

double HiggsChiral::muVBFHmumu(const double sqrt_s) const
{
    return computecV() * computecV() * computecmu() * computecmu() / computeGammaTotalRatio();
}

double HiggsChiral::muZHmumu(const double sqrt_s) const
{
    return computecV() * computecV() * computecmu() * computecmu() / computeGammaTotalRatio();
}

double HiggsChiral::muWHmumu(const double sqrt_s) const
{
    return computecV() * computecV() * computecmu() * computecmu() / computeGammaTotalRatio();
}

double HiggsChiral::muVHmumu(const double sqrt_s) const
{
    return computecV() * computecV() * computecmu() * computecmu() / computeGammaTotalRatio();
}

double HiggsChiral::muttHmumu(const double sqrt_s) const
{
    return computect() * computect() * computecmu() * computecmu() / computeGammaTotalRatio();
}

double HiggsChiral::muggHtautau(const double sqrt_s) const
{
    return computecg() * computecg() * computectau() * computectau() / computeGammaTotalRatio();
}

double HiggsChiral::muVBFHtautau(const double sqrt_s) const
{
    return computecV() * computecV() * computectau() * computectau() / computeGammaTotalRatio();
}

double HiggsChiral::muZHtautau(const double sqrt_s) const
{
    return computecV() * computecV() * computectau() * computectau() / computeGammaTotalRatio();
}

double HiggsChiral::muWHtautau(const double sqrt_s) const
{
    return computecV() * computecV() * computectau() * computectau() / computeGammaTotalRatio();
}

double HiggsChiral::muVHtautau(const double sqrt_s) const
{
    return computecV() * computecV() * computectau() * computectau() / computeGammaTotalRatio();
}

double HiggsChiral::muttHtautau(const double sqrt_s) const
{
    return computect() * computect() * computectau() * computectau() / computeGammaTotalRatio();
}

double HiggsChiral::muggHbb(const double sqrt_s) const
{
    return computecg() * computecg() * computecb() * computecb() / computeGammaTotalRatio();
}

double HiggsChiral::muVBFHbb(const double sqrt_s) const
{
    return computecV() * computecV() * computecb() * computecb() / computeGammaTotalRatio();
}

double HiggsChiral::muZHbb(const double sqrt_s) const
{
    return computecV() * computecV() * computecb() * computecb() / computeGammaTotalRatio();
}

double HiggsChiral::muWHbb(const double sqrt_s) const
{
    return computecV() * computecV() * computecb() * computecb() / computeGammaTotalRatio();
}

double HiggsChiral::muVHbb(const double sqrt_s) const
{
    return computecV() * computecV() * computecb() * computecb() / computeGammaTotalRatio();
}

double HiggsChiral::muttHbb(const double sqrt_s) const
{
    return computect() * computect() * computecb() * computecb() / computeGammaTotalRatio();
}

double HiggsChiral::muppHmumu(const double sqrt_s) const
{
    if(sqrt_s==8)
    {
        return (0.872 * computecg() * computecg() + 0.122 * computecV() * computecV() + 0.006 * computect() * computect()) * computecmu() * computecmu() / computeGammaTotalRatio();
    }
    if(sqrt_s==13)
    {
        return (0.871 * computecg() * computecg() + 0.119 * computecV() * computecV() + 0.010 * computect() * computect()) * computecmu() * computecmu() / computeGammaTotalRatio();
    }
    else
    {
        throw std::runtime_error("The observable muppHmumu is only defined for 8 or 13 TeV.");
    }
}


////////////////////////////////////////////////////////////////////////////////////////////
//-----------------------------------------------------------------------------------------
//-- Special Hadron collider signal strengths with separate full TH unc U(prod x decay) ---
//-----------------------------------------------------------------------------------------
//////////////////////////////////////////////////////////////////////////////////////////// 

double HiggsChiral::muTHUggHgaga(const double sqrt_s) const
{
    return muggH(sqrt_s)*BrHgagaRatio();
}
   
double HiggsChiral::muTHUVBFHgaga(const double sqrt_s) const
{
    return muVBF(sqrt_s)*BrHgagaRatio();
}

double HiggsChiral::muTHUZHgaga(const double sqrt_s) const
{
    return muZH(sqrt_s)*BrHgagaRatio();
}

double HiggsChiral::muTHUWHgaga(const double sqrt_s) const
{
    return muWH(sqrt_s)*BrHgagaRatio();
}

double HiggsChiral::muTHUVHgaga(const double sqrt_s) const
{
    return muVH(sqrt_s)*BrHgagaRatio();
}

double HiggsChiral::muTHUttHgaga(const double sqrt_s) const
{
    return muttH(sqrt_s)*BrHgagaRatio();
}

double HiggsChiral::muTHUggHZga(const double sqrt_s) const
{
    return muggH(sqrt_s)*BrHZgaRatio();
}

double HiggsChiral::muTHUVBFHZga(const double sqrt_s) const
{
    return muVBF(sqrt_s)*BrHZgaRatio();
}

double HiggsChiral::muTHUZHZga(const double sqrt_s) const
{
    return muZH(sqrt_s)*BrHZgaRatio();
}

double HiggsChiral::muTHUWHZga(const double sqrt_s) const
{
    return muWH(sqrt_s)*BrHZgaRatio();
}

double HiggsChiral::muTHUVHZga(const double sqrt_s) const
{
    return muVH(sqrt_s)*BrHZgaRatio();
}

double HiggsChiral::muTHUttHZga(const double sqrt_s) const
{
    return muttH(sqrt_s)*BrHZgaRatio();
}

double HiggsChiral::muTHUggHZZ(const double sqrt_s) const
{
    return muggH(sqrt_s)*BrHZZRatio();
}

double HiggsChiral::muTHUVBFHZZ(const double sqrt_s) const
{
    return muVBF(sqrt_s)*BrHZZRatio();
}

double HiggsChiral::muTHUZHZZ(const double sqrt_s) const
{
    return muZH(sqrt_s)*BrHZZRatio();
}

double HiggsChiral::muTHUWHZZ(const double sqrt_s) const
{
    return muWH(sqrt_s)*BrHZZRatio();
}

double HiggsChiral::muTHUVHZZ(const double sqrt_s) const
{
    return muVH(sqrt_s)*BrHZZRatio();
}

double HiggsChiral::muTHUttHZZ(const double sqrt_s) const
{
    return muttH(sqrt_s)*BrHZZRatio();
}

double HiggsChiral::muTHUggHZZ4l(const double sqrt_s) const
{
    return muggH(sqrt_s)*BrHZZ4lRatio();
}

double HiggsChiral::muTHUVBFHZZ4l(const double sqrt_s) const
{
    return muVBF(sqrt_s)*BrHZZ4lRatio();
}

double HiggsChiral::muTHUZHZZ4l(const double sqrt_s) const
{
    return muZH(sqrt_s)*BrHZZ4lRatio();
}

double HiggsChiral::muTHUWHZZ4l(const double sqrt_s) const
{
    return muWH(sqrt_s)*BrHZZ4lRatio();
}

double HiggsChiral::muTHUVHZZ4l(const double sqrt_s) const
{
    return muVH(sqrt_s)*BrHZZ4lRatio();
}

double HiggsChiral::muTHUttHZZ4l(const double sqrt_s) const
{
    return muttH(sqrt_s)*BrHZZ4lRatio();
}

double HiggsChiral::muTHUggHWW(const double sqrt_s) const
{
    return muggH(sqrt_s)*BrHWWRatio();
}

double HiggsChiral::muTHUVBFHWW(const double sqrt_s) const
{
    return muVBF(sqrt_s)*BrHWWRatio();
}

double HiggsChiral::muTHUZHWW(const double sqrt_s) const
{
    return muZH(sqrt_s)*BrHWWRatio();
}

double HiggsChiral::muTHUWHWW(const double sqrt_s) const
{
    return muWH(sqrt_s)*BrHWWRatio();
}

double HiggsChiral::muTHUVHWW(const double sqrt_s) const
{
    return muVH(sqrt_s)*BrHWWRatio();
}

double HiggsChiral::muTHUttHWW(const double sqrt_s) const
{
    return muttH(sqrt_s)*BrHWWRatio();
}

double HiggsChiral::muTHUggHWW2l2v(const double sqrt_s) const
{
    return muggH(sqrt_s)*BrHWW2l2vRatio();
}

double HiggsChiral::muTHUVBFHWW2l2v(const double sqrt_s) const
{
    return muVBF(sqrt_s)*BrHWW2l2vRatio();
}

double HiggsChiral::muTHUZHWW2l2v(const double sqrt_s) const
{
    return muZH(sqrt_s)*BrHWW2l2vRatio();
}

double HiggsChiral::muTHUWHWW2l2v(const double sqrt_s) const
{
    return muWH(sqrt_s)*BrHWW2l2vRatio();
}

double HiggsChiral::muTHUVHWW2l2v(const double sqrt_s) const
{
    return muVH(sqrt_s)*BrHWW2l2vRatio();
}

double HiggsChiral::muTHUttHWW2l2v(const double sqrt_s) const
{
    return muttH(sqrt_s)*BrHWW2l2vRatio();
}

double HiggsChiral::muTHUggHmumu(const double sqrt_s) const
{
    return muggH(sqrt_s)*BrHmumuRatio();
}

double HiggsChiral::muTHUVBFHmumu(const double sqrt_s) const
{
    return muVBF(sqrt_s)*BrHmumuRatio();
}

double HiggsChiral::muTHUZHmumu(const double sqrt_s) const
{
    return muZH(sqrt_s)*BrHmumuRatio();
}

double HiggsChiral::muTHUWHmumu(const double sqrt_s) const
{
    return muWH(sqrt_s)*BrHmumuRatio();
}

double HiggsChiral::muTHUVHmumu(const double sqrt_s) const
{
    return muVH(sqrt_s)*BrHmumuRatio();
}

double HiggsChiral::muTHUttHmumu(const double sqrt_s) const
{
    return muttH(sqrt_s)*BrHmumuRatio();
}

double HiggsChiral::muTHUggHtautau(const double sqrt_s) const
{
    return muggH(sqrt_s)*BrHtautauRatio();
}

double HiggsChiral::muTHUVBFHtautau(const double sqrt_s) const
{
    return muVBF(sqrt_s)*BrHtautauRatio();
}

double HiggsChiral::muTHUZHtautau(const double sqrt_s) const
{
    return muZH(sqrt_s)*BrHtautauRatio();
}

double HiggsChiral::muTHUWHtautau(const double sqrt_s) const
{
    return muWH(sqrt_s)*BrHtautauRatio();
}

double HiggsChiral::muTHUVHtautau(const double sqrt_s) const
{
    return muVH(sqrt_s)*BrHtautauRatio();
}

double HiggsChiral::muTHUttHtautau(const double sqrt_s) const
{
    return muttH(sqrt_s)*BrHtautauRatio();
}

double HiggsChiral::muTHUggHbb(const double sqrt_s) const
{
    return muggH(sqrt_s)*BrHbbRatio();
}

double HiggsChiral::muTHUVBFHbb(const double sqrt_s) const
{
    return muVBF(sqrt_s)*BrHbbRatio();
}

double HiggsChiral::muTHUZHbb(const double sqrt_s) const
{
    return muZH(sqrt_s)*BrHbbRatio();
}

double HiggsChiral::muTHUWHbb(const double sqrt_s) const
{
    return muWH(sqrt_s)*BrHbbRatio();
}

double HiggsChiral::muTHUVHbb(const double sqrt_s) const
{
    return muVH(sqrt_s)*BrHbbRatio();
}

double HiggsChiral::muTHUttHbb(const double sqrt_s) const
{
    return muttH(sqrt_s)*BrHbbRatio();
}

double HiggsChiral::muTHUVBFBRinv(const double sqrt_s) const
{
    return muVBF(sqrt_s)*Br_H_inv();
}

double HiggsChiral::muTHUVBFHinv(const double sqrt_s) const
{
    return muVBF(sqrt_s)*BrHtoinvRatio();
}

double HiggsChiral::muTHUVHBRinv(const double sqrt_s) const
{
    return muVH(sqrt_s)*Br_H_inv();
}

double HiggsChiral::muTHUVHinv(const double sqrt_s) const
{
    return muVH(sqrt_s)*BrHtoinvRatio();
}

double HiggsChiral::muTHUggHZZ4mu(const double sqrt_s) const
{
    return muggH(sqrt_s)*BrHZZ4muRatio();
}

double HiggsChiral::muTHUggHZgamumu(const double sqrt_s) const
{
    return muggH(sqrt_s)*BrHZgamumuRatio();
}

////////////////////////////////////////////////////////////////////////////////////////////
//-----------------------------------------------------------------------------------------
//-- Decays ---
//-----------------------------------------------------------------------------------------
//////////////////////////////////////////////////////////////////////////////////////////// 


double HiggsChiral::computeGammaTotalRatio() const
{
    
    return ((computecg() * computecg() * trueSM.computeBrHtogg()
            + computecV() * computecV() * trueSM.computeBrHtoWW()
            + computecV() * computecV() * trueSM.computeBrHtoZZ()
            + computecZga() * computecZga() * trueSM.computeBrHtoZga()
            + computecgaga() * computecgaga() * trueSM.computeBrHtogaga()
            + computecmu() * computecmu() * trueSM.computeBrHtomumu()
            + computectau() * computectau() * trueSM.computeBrHtotautau()
            + computecc() * computecc() * trueSM.computeBrHtocc()
            + computecb() * computecb() * trueSM.computeBrHtobb())
            / (trueSM.computeBrHtogg()
            + trueSM.computeBrHtoWW()
            + trueSM.computeBrHtoZZ()
            + trueSM.computeBrHtoZga()
            + trueSM.computeBrHtogaga()
            + trueSM.computeBrHtomumu()
            + trueSM.computeBrHtotautau()
            + trueSM.computeBrHtocc()
            + trueSM.computeBrHtobb()));
    
//    return (GammaTotal() / trueSM.computeGammaHTotal());
}

////////////////////////////////////////////////////////////////////////

double HiggsChiral::computecg() const
{
    if (loopComputed) return cg_loop;
    
    double Mt=trueSM.getQuarks(QCD::TOP).getMass();
    double Mb=trueSM.getQuarks(QCD::BOTTOM).getMass();
    double Mc=trueSM.getQuarks(QCD::CHARM).getMass();
    double TAUt=4.0*Mt*Mt/(mHl*mHl);
    double TAUb=4.0*Mb*Mb/(mHl*mHl);
    double TAUc=4.0*Mc*Mc/(mHl*mHl);
    double cgEff;

    cgEff = ( ( ct * 0.5*TAUt*(1.0+(1.0-TAUt)*f_func(TAUt))
                     +cb * 0.5*TAUb*(1.0+(1.0-TAUb)*f_func(TAUb)) 
                     +cc * 0.5*TAUc*(1.0+(1.0-TAUc)*f_func(TAUc)) ) * (1.0+11.0*AlsMz/(4.0*M_PI))
                     +cg * 0.5 ).abs2();

    cgEff = cgEff / ( ( 0.5*TAUt*(1.0+(1.0-TAUt)*f_func(TAUt))
                     + 0.5*TAUb*(1.0+(1.0-TAUb)*f_func(TAUb))
                     + 0.5*TAUc*(1.0+(1.0-TAUc)*f_func(TAUc)) ) * (1.0+11.0*AlsMz/(4.0*M_PI))).abs2();

    return (sqrt(cgEff));
}

double HiggsChiral::computecV() const
{
    return cv;
}

double HiggsChiral::computecZga() const
{
    if (loopComputed) return cZga_loop;
    
    double Mt=trueSM.getQuarks(QCD::TOP).getMass();
    double Mb=trueSM.getQuarks(QCD::BOTTOM).getMass();
    double Mc=trueSM.getQuarks(QCD::CHARM).getMass();
    double Mtau=trueSM.getLeptons(StandardModel::TAU).getMass();
    double Mmu=trueSM.getLeptons(StandardModel::MU).getMass();
    double MW=trueSM.Mw_tree();
    double cW2=trueSM.c02();
    double sW2=1.0-cW2;
    double TAUt=4.0*Mt*Mt/(mHl*mHl);
    double TAUb=4.0*Mb*Mb/(mHl*mHl);
    double TAUc=4.0*Mc*Mc/(mHl*mHl);
    double TAUtau=4.0*Mtau*Mtau/(mHl*mHl);
    double TAUmu=4.0*Mmu*Mmu/(mHl*mHl);
    double TAUw=4.0*MW*MW/(mHl*mHl);
    double LAMt=4.0*Mt*Mt/(Mz*Mz);
    double LAMb=4.0*Mb*Mb/(Mz*Mz);
    double LAMc=4.0*Mc*Mc/(Mz*Mz);
    double LAMtau=4.0*Mtau*Mtau/(Mz*Mz);
    double LAMmu=4.0*Mmu*Mmu/(Mz*Mz);
    double LAMw=4.0*MW*MW/(Mz*Mz);
    double cZgaEff;

    cZgaEff = ((-ct * 4.0*(0.5-4.0/3.0*sW2)*(Int1(TAUt,LAMt)-Int2(TAUt,LAMt)) * (1.0-AlsMz/M_PI)
                      +cb * 2.0*(-0.5+2.0/3.0*sW2)*(Int1(TAUb,LAMb)-Int2(TAUb,LAMb))
                      -cc * 4.0*(0.5-4.0/3.0*sW2)*(Int1(TAUc,LAMc)-Int2(TAUc,LAMc))            
                      +ctau * 2.0*(-0.5+2.0*sW2)*(Int1(TAUtau,LAMtau)-Int2(TAUtau,LAMtau))
                      +cmu * 2.0*(-0.5+2.0*sW2)*(Int1(TAUmu,LAMmu)-Int2(TAUmu,LAMmu)) )/sqrt(sW2*cW2)
                     -cv * sqrt(cW2/sW2)*(4.0*(3.0-sW2/cW2)*Int2(TAUw,LAMw)
                            +((1.0+2.0/TAUw)*sW2/cW2-(5.0+2.0/TAUw))*Int1(TAUw,LAMw))
                     +cZga).abs2();
    
    cZgaEff = cZgaEff / ((-4.0*(0.5-4.0/3.0*sW2)*(Int1(TAUt,LAMt)-Int2(TAUt,LAMt)) * (1.0-AlsMz/M_PI)
                      + 2.0*(-0.5+2.0/3.0*sW2)*(Int1(TAUb,LAMb)-Int2(TAUb,LAMb))
                      - 4.0*(0.5-4.0/3.0*sW2)*(Int1(TAUc,LAMc)-Int2(TAUc,LAMc))
                      + 2.0*(-0.5+2.0*sW2)*(Int1(TAUtau,LAMtau)-Int2(TAUtau,LAMtau))
                      + 2.0*(-0.5+2.0*sW2)*(Int1(TAUmu,LAMmu)-Int2(TAUmu,LAMmu)) )/sqrt(sW2*cW2)
                      - sqrt(cW2/sW2)*(4.0*(3.0-sW2/cW2)*Int2(TAUw,LAMw)
                            +((1.0+2.0/TAUw)*sW2/cW2-(5.0+2.0/TAUw))*Int1(TAUw,LAMw))).abs2();
    
    return (sqrt(cZgaEff));
}

double HiggsChiral::computecgaga() const
{
    if (loopComputed) return cga_loop;
    
    double Mt=trueSM.getQuarks(QCD::TOP).getMass();
    double Mb=trueSM.getQuarks(QCD::BOTTOM).getMass();
    double Mc=trueSM.getQuarks(QCD::CHARM).getMass();
    double Mtau=trueSM.getLeptons(StandardModel::TAU).getMass();
    double Mmu=trueSM.getLeptons(StandardModel::MU).getMass();
    double MW=trueSM.Mw_tree();
    double TAUt=4.0*Mt*Mt/(mHl*mHl);
    double TAUb=4.0*Mb*Mb/(mHl*mHl);
    double TAUc=4.0*Mc*Mc/(mHl*mHl);
    double TAUtau=4.0*Mtau*Mtau/(mHl*mHl);
    double TAUmu=4.0*Mmu*Mmu/(mHl*mHl);
    double TAUw=4.0*MW*MW/(mHl*mHl);
    double cgagaEff;

    cgagaEff = ( ct * (8./3.)*TAUt*(1.+(1.-TAUt)*f_func(TAUt)) * (1.0-AlsMz/M_PI)
                     +cb * (2./3.)*TAUb*(1.+(1.-TAUb)*f_func(TAUb))
                     +cc * (8./3.)*TAUc*(1.+(1.-TAUc)*f_func(TAUc))
                     +ctau * 2.0*TAUtau*(1.+(1.-TAUtau)*f_func(TAUtau))
                     +cmu * 2.0*TAUmu*(1.+(1.-TAUmu)*f_func(TAUmu))
                     -cv * (2.0+3.0*TAUw+3.0*TAUw*(2.0-TAUw)*f_func(TAUw))
                     +cga * 2.0 ).abs2();
    
    cgagaEff = cgagaEff / ( (8./3.)*TAUt*(1.+(1.-TAUt)*f_func(TAUt)) * (1.0-AlsMz/M_PI)
                     + (2./3.)*TAUb*(1.+(1.-TAUb)*f_func(TAUb))
                     + (8./3.)*TAUc*(1.+(1.-TAUc)*f_func(TAUc))
                     + 2.0*TAUtau*(1.+(1.-TAUtau)*f_func(TAUtau))
                     + 2.0*TAUmu*(1.+(1.-TAUmu)*f_func(TAUmu))
                     - (2.0+3.0*TAUw+3.0*TAUw*(2.0-TAUw)*f_func(TAUw)) ).abs2();
    
    return (sqrt(cgagaEff));
}

double HiggsChiral::computecmu() const
{
    return cmu;
}

double HiggsChiral::computectau() const
{
    return ctau;
}

double HiggsChiral::computecc() const
{
    return cc;
}

double HiggsChiral::computect() const
{
    return ct;
}

double HiggsChiral::computecb() const
{
    return cb;
}


////////////////////////////////////////////////////////////////////////

double HiggsChiral::deltaG_hgg() const
{
    // Not needed/implemented yet
    return 0.0;
}

double HiggsChiral::deltaG_hggRatio() const{
    return (computecg() - 1.0);
}

double HiggsChiral::deltaG1_hWW() const{
    return 0.0;
}

double HiggsChiral::deltaG2_hWW() const{
    return 0.0;
}

double HiggsChiral::deltaG3_hWW() const{
    
    double gSM = 2.0 * (trueSM.Mw_tree())* (trueSM.Mw_tree()) / (trueSM.v());
    
    return ( gSM * (computecV() - 1.0 ) );
}

double HiggsChiral::deltaG1_hZZ() const{
    return 0.0;
}

double HiggsChiral::deltaG2_hZZ() const{
    return 0.0;
}

double HiggsChiral::deltaG3_hZZ() const{
    
    double gSM = (trueSM.getMz()) * (trueSM.getMz()) / (trueSM.v());
    
    return ( gSM * (computecV() - 1.0 ) );
}

double HiggsChiral::deltaG1_hZA() const{
    // Not needed/implemented yet
    return 0.0;
}

double HiggsChiral::deltaG1_hZARatio() const{
    return (computecZga() - 1.0);
}

double HiggsChiral::deltaG2_hZA() const{
    return 0.0;
}

double HiggsChiral::deltaG_hAA() const{
    // Not needed/implemented yet
    return 0.0;
}

double HiggsChiral::deltaG_hAARatio() const{
    return (computecgaga() - 1.0);
}

gslpp::complex HiggsChiral::deltaG_hff(const Particle p) const{
    
    double gSM = -(p.getMass()) / (trueSM.v());
    
    if ( p.is("ELECTRON") ) {
        return ( 0.0 );
    } else if ( p.is("MU") ) {
         return ( gSM * (computecmu() - 1.0 ) );       
    } else if ( p.is("TAU") ) {
        return ( gSM * (computectau() - 1.0 ) );        
    } else if ( p.is("UP") ) {
        return ( 0.0 );        
    } else if ( p.is("CHARM") ) {
        return ( gSM * (computecc() - 1.0 ) );        
    } else if ( p.is("TOP") ) {
        return ( gSM * (computect() - 1.0 ) );        
    } else if ( p.is("DOWN") ) {
        return ( 0.0 );        
    } else if ( p.is("STRANGE") ) {
        return ( 0.0 );        
    } else if ( p.is("BOTTOM") ) {
        return ( gSM * (computecb() - 1.0 ) );        
    } else {
        return 0.0;
    }
}

///////////////////////////////////////////////////////////////////////////////

double HiggsChiral::kappamueff() const
{
      return computecmu();
}

double HiggsChiral::kappataueff() const
{
      return computectau();
}

double HiggsChiral::kappaceff() const
{
      return computecc();
}

double HiggsChiral::kappabeff() const
{
      return computecb();
}

double HiggsChiral::kappaGeff() const
{
      return computecg();
}

double HiggsChiral::kappaZeff() const
{
      return computecV();
}

double HiggsChiral::kappaWeff() const
{
      return computecV();
}

double HiggsChiral::kappaAeff() const
{
      return computecgaga();
}

double HiggsChiral::kappaZAeff() const
{
      return computecZga();
}

//////////////////////////////////////////////////////////////////////

gslpp::complex HiggsChiral::f_func(const double x) const{
    if(x<1) {
    gslpp::complex z = -gslpp::complex::i()*M_PI;
    return -pow(log((1+sqrt(1-x))/(1-sqrt(1-x)))+z,2)/4.0;
    }
    else {
        return pow(asin(sqrt(1.0/x)),2);
    }
}

gslpp::complex HiggsChiral::g_func(const double x) const{
    if(x<1) {
    gslpp::complex z = -gslpp::complex::i()*M_PI;
    gslpp::complex gs1 = sqrt(1.0-x)*(log((1.0+sqrt(1.0-x))/(1.0-sqrt(1.0-x)))+z)/2.0;
    return gs1;
    }
    else {
        gslpp::complex gg1 = sqrt(x-1.0)*asin(sqrt(1.0/x));
        return gg1;
    }
}

gslpp::complex HiggsChiral::Int1(const double tau, const double lambda) const{
    return tau*lambda/(tau-lambda)/2.0+tau*tau*lambda*lambda/((tau-lambda)
           *(tau-lambda))/2.0*(f_func(tau)-f_func(lambda))+tau*tau*lambda/((tau-lambda)
           *(tau-lambda))*(g_func(tau)-g_func(lambda));
}

gslpp::complex HiggsChiral::Int2(const double tau, const double lambda) const{
    return -tau*lambda/(tau-lambda)/2.0*(f_func(tau)-f_func(lambda));
}
