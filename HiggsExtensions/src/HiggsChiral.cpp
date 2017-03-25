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
    return Gammagg() / (trueSM.computeBrHtogg()*trueSM.computeGammaHTotal());
}

double HiggsChiral::muVBF(const double sqrt_s) const
{
    return computecV() * computecV();
}

double HiggsChiral::mueeWBF(const double sqrt_s) const
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

double HiggsChiral::Gammagg() const
{
    double Mt=trueSM.getQuarks(QCD::TOP).getMass();
    double Mb=trueSM.getQuarks(QCD::BOTTOM).getMass();
    double TAUt=4.0*Mt*Mt/(mHl*mHl);
    double TAUb=4.0*Mb*Mb/(mHl*mHl);

    double factor = GF*AlsMz*AlsMz*mHl*mHl*mHl/(sqrt(2.0)*4.0*M_PI*M_PI*M_PI);

    return factor * ( ct * 0.5*TAUt*(1.0+(1.0-TAUt)*f_func(TAUt)) * (1.0+11.0*AlsMz/(4.0*M_PI))
                     +cb * 0.5*TAUb*(1.0+(1.0-TAUb)*f_func(TAUb))
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
                     +cZga * 0.0 ).abs2(); /*CHANGE THE cZga COEFFICIENT!*/
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
    return 0.0;
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

double HiggsChiral::BrHWWRatio() const
{
    return GammaWW() / GammaTotal() / trueSM.computeBrHtoWW();
}

double HiggsChiral::BrHZZRatio() const
{
    return GammaZZ() / GammaTotal() / trueSM.computeBrHtoZZ();
}

double HiggsChiral::BrHZgaRatio() const
{
    return GammaZga() / GammaTotal() / trueSM.computeBrHtoZga();
}

double HiggsChiral::BrHgagaRatio() const
{
    return Gammagaga() / GammaTotal() / trueSM.computeBrHtogaga();
}

double HiggsChiral::BrHmumuRatio() const
{
    return Gammamumu() / GammaTotal() / trueSM.computeBrHtomumu();
}

double HiggsChiral::BrHtautauRatio() const
{
    return Gammatautau() / GammaTotal() / trueSM.computeBrHtotautau();
}

double HiggsChiral::BrHccRatio() const
{
    return Gammacc() / GammaTotal() / trueSM.computeBrHtocc();
}

double HiggsChiral::BrHbbRatio() const
{
    return Gammabb() / GammaTotal() / trueSM.computeBrHtobb();
}
//
//double HiggsChiral::computeGammaTotalRatio() const
//{
//    double gtt_SM = trueSM.computeGammaHgg_tt();
//    double gbb_SM = trueSM.computeGammaHgg_bb();
//    double gtb_SM = trueSM.computeGammaHgg_tb();
//    double Ggg_over_GggSM = (computeKt() * computeKt() * gtt_SM
//            + computeKb() * computeKb() * gbb_SM
//            + computeKt() * computeKb() * gtb_SM) / (gtt_SM + gbb_SM + gtb_SM);
//    return ((Ggg_over_GggSM * trueSM.computeBrHtogg()
//            + computeKW() * computeKW() * trueSM.computeBrHtoWW()
//            + computeKZ() * computeKZ() * trueSM.computeBrHtoZZ()
//            + computeKZga() * computeKZga() * trueSM.computeBrHtoZga()
//            + computeKgaga() * computeKgaga() * trueSM.computeBrHtogaga()
//            + computeKmu() * computeKmu() * trueSM.computeBrHtomumu()
//            + computeKtau() * computeKtau() * trueSM.computeBrHtotautau()
//            + computeKc() * computeKc() * trueSM.computeBrHtocc()
//            + computeKb() * computeKb() * trueSM.computeBrHtobb())
//            / (trueSM.computeBrHtogg()
//            + trueSM.computeBrHtoWW()
//            + trueSM.computeBrHtoZZ()
//            + trueSM.computeBrHtoZga()
//            + trueSM.computeBrHtogaga()
//            + trueSM.computeBrHtomumu()
//            + trueSM.computeBrHtotautau()
//            + trueSM.computeBrHtocc()
//            + trueSM.computeBrHtobb()));
//}

////////////////////////////////////////////////////////////////////////

double HiggsChiral::computecg() const
{
    return cg;
}

double HiggsChiral::computecV() const
{
    return cv;
}

double HiggsChiral::computecZga() const
{
    return cZga;
}

double HiggsChiral::computecga() const
{
    return cga;
}

double HiggsChiral::computecmu() const
{
    return cmu;
}

double HiggsChiral::computectau() const
{
    return ctau;
}

//double HiggsKvKfgen::computeKc() const
//{
//    return cc;
//}

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
