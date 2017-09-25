/* 
 * Copyright (C) 2017 HEPfit Collaboration
 *
 *
 * For the licensing terms see doc/COPYING.
 */

#include "HiggsChiral.h"

const std::string HiggsChiral::HChiralvars[NHChiralvars] = {
    "cv", "ct", "cb", "cc", "ctau", "cmu", "cg", "cga", "cZga", "obsZgaLimitATLAS13", "obsZgaLimitATLAS", "obsZgaLimitCMS", "expZgaLimitATLAS13", "expZgaLimitATLAS", "expZgaLimitCMS"
};

HiggsChiral::HiggsChiral()
: NPbase()
{
    ModelParamMap.insert(std::pair<std::string, boost::reference_wrapper<const double> >("cv", boost::cref(cv)));
    ModelParamMap.insert(std::pair<std::string, boost::reference_wrapper<const double> >("ct", boost::cref(ct)));
    ModelParamMap.insert(std::pair<std::string, boost::reference_wrapper<const double> >("cb", boost::cref(cb)));
    ModelParamMap.insert(std::pair<std::string, boost::reference_wrapper<const double> >("cc", boost::cref(cc)));
    ModelParamMap.insert(std::pair<std::string, boost::reference_wrapper<const double> >("ctau", boost::cref(ctau)));
    ModelParamMap.insert(std::pair<std::string, boost::reference_wrapper<const double> >("cmu", boost::cref(cmu)));
    ModelParamMap.insert(std::pair<std::string, boost::reference_wrapper<const double> >("cg", boost::cref(cg)));
    ModelParamMap.insert(std::pair<std::string, boost::reference_wrapper<const double> >("cga", boost::cref(cga)));
    ModelParamMap.insert(std::pair<std::string, boost::reference_wrapper<const double> >("cZga", boost::cref(cZga)));
    ModelParamMap.insert(std::pair<std::string, boost::reference_wrapper<const double> >("obsZgaLimitATLAS13", boost::cref(obsZgaLimitATLAS13)));
    ModelParamMap.insert(std::pair<std::string, boost::reference_wrapper<const double> >("obsZgaLimitATLAS", boost::cref(obsZgaLimitATLAS)));
    ModelParamMap.insert(std::pair<std::string, boost::reference_wrapper<const double> >("obsZgaLimitCMS", boost::cref(obsZgaLimitCMS)));
    ModelParamMap.insert(std::pair<std::string, boost::reference_wrapper<const double> >("expZgaLimitATLAS13", boost::cref(expZgaLimitATLAS13)));
    ModelParamMap.insert(std::pair<std::string, boost::reference_wrapper<const double> >("expZgaLimitATLAS", boost::cref(expZgaLimitATLAS)));
    ModelParamMap.insert(std::pair<std::string, boost::reference_wrapper<const double> >("expZgaLimitCMS", boost::cref(expZgaLimitCMS)));
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
    else if (name.compare("obsZgaLimitATLAS") == 0)
        obsZgaLimitATLAS = value;
    else if (name.compare("obsZgaLimitCMS") == 0)
        obsZgaLimitCMS = value;
    else if (name.compare("expZgaLimitATLAS13") == 0)
        expZgaLimitATLAS13 = value;
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
            std::cout << "missing mandatory HiggsChiral parameter " << HChiralvars[i] << std::endl;
            return false;
        }
    }
    return (NPbase::CheckParameters(DPars));
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

double HiggsChiral::mueeWBF(const double sqrt_s) const
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

double HiggsChiral::BrHZZRatio() const
{
//    return GammaZZ() / GammaTotal() / trueSM.computeBrHtoZZ();
    
    return (computecV() * computecV() / computeGammaTotalRatio());
}

double HiggsChiral::BrHZgaRatio() const
{
//    return GammaZga() / GammaTotal() / trueSM.computeBrHtoZga();
    
    return (computecZga() * computecZga() / computeGammaTotalRatio());
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
//
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
                     +cZga).abs2(); /*CHANGE THE cZga COEFFICIENT!*/
    
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
