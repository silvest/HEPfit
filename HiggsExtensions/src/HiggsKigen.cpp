/* 
 * Copyright (C) 2014 HEPfit Collaboration
 *
 *
 * For the licensing terms see doc/COPYING.
 */

#include "HiggsKigen.h"

const std::string HiggsKigen::HKvKfgenvars[NHKvKfgenvars] = {
    "Kw", "Kz", "Kg", "Kga", "Kzga", "Ku", "Kc", "Kt", "Kd", "Ks", "Kb", "Ke", "Kmu", "Ktau", "BrHinv", "BrHexo"
};

HiggsKigen::HiggsKigen()
: NPbase()
{
    FlagKiLoop = true;
    FlagCustodial = false;
    
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
    ModelParamMap.insert(std::pair<std::string, boost::reference_wrapper<const double> >("BrHexo", boost::cref(BrHexo)));
}

bool HiggsKigen::PostUpdate()
{
    if (!NPbase::PostUpdate()) return (false);
    
    if (!FlagKiLoop) {
        Kg = computeKgLoop();
        Kga = computeKgagaLoop();
        Kzga = computeKZgaLoop();
    }

    if (FlagCustodial) {
        Kz = Kw;
    }

    return (true);
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
//  Always positive
        BrHinv = fabs(value);
    else if (name.compare("BrHexo") == 0)
//  Always positive
        BrHexo = fabs(value);
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

bool HiggsKigen::setFlag(const std::string name, const bool value)
{
    bool res = false;
    if (name.compare("KiLoop") == 0) {
        FlagKiLoop = value;
        res = true;
    } else if (name.compare("Custodial") == 0) {
        FlagCustodial = value;
        res = true;
    } else
        res = NPbase::setFlag(name, value);

    return (res);
}

////////////////////////////////////////////////////////////////////////

double HiggsKigen::obliqueS() const
{
    if (FlagCustodial) {
        double Lambda, Kv;
    
        Kv=Kw;
    
        if (fabs(1.0 - Kv * Kv) < pow(10.0, -32.0))
            Lambda = pow(10.0, 19.0);
        else
            Lambda = 4.0 * M_PI * v() / sqrt(fabs(1.0 - Kv * Kv));

        return ( 1.0 / 12.0 / M_PI * (1.0 - Kv * Kv) * log(Lambda * Lambda / mHl / mHl));
    
    } else {
        return 0.0;
    }
}

double HiggsKigen::obliqueT() const
{
    if (FlagCustodial) {
        double Lambda, Kv;
        double cW2_SM = trueSM.cW2();
        
        Kv=Kw;
        
        if (fabs(1.0 - Kv * Kv) < pow(10.0, -32.0))
            Lambda = pow(10.0, 19.0);
        else
            Lambda = 4.0 * M_PI * v() / sqrt(fabs(1.0 - Kv * Kv));

        return ( -3.0 / 16.0 / M_PI / cW2_SM * (1.0 - Kv * Kv) * log(Lambda * Lambda / mHl / mHl));
    
    } else {
        return 0.0;
    }
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

double HiggsKigen::muVBFgamma(const double sqrt_s) const
{
    // Assume photon comes from the internal W lines -> no ZBF
    
    return (computeKW() * computeKW());
}

double HiggsKigen::mueeWBF(const double sqrt_s) const
{
    return (computeKW() * computeKW());
}

double HiggsKigen::mueeWBFPol(const double sqrt_s, const double Pol_em, const double Pol_ep) const
{
    return (computeKW() * computeKW());
}

double HiggsKigen::mueeZBF(const double sqrt_s) const
{
    return (computeKZ() * computeKZ());
}

double HiggsKigen::mueeZBFPol(const double sqrt_s, const double Pol_em, const double Pol_ep) const
{
    return (computeKZ() * computeKZ());
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

double HiggsKigen::mueeZHPol(const double sqrt_s, const double Pol_em, const double Pol_ep) const
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

double HiggsKigen::mummH(const double sqrt_s) const
{
    return (computeKmu() * computeKmu());
}

double HiggsKigen::mueettHPol(const double sqrt_s, const double Pol_em, const double Pol_ep) const
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

double HiggsKigen::BrHWW2l2vRatio() const
{
    return BrHWWRatio();
}

double HiggsKigen::BrHZZRatio() const
{
    return (computeKZ() * computeKZ() / computeGammaTotalRatio());
}

double HiggsKigen::BrHZZ4lRatio() const
{
    return BrHZZRatio();
}

double HiggsKigen::BrHZgaRatio() const
{
    return (computeKZga() * computeKZga() / computeGammaTotalRatio());
}

double HiggsKigen::BrHZgallRatio() const
{
    return BrHZgaRatio();
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


double HiggsKigen::muggHgaga(const double sqrt_s) const
{
    return muggH(sqrt_s)*BrHgagaRatio();
}
    
double HiggsKigen::muggHgagaInt(const double sqrt_s) const
{
    double muNWA, GammaRatio, fki;
    double kt,kb,kc,ks,ku,kd,kta,km,ke;
    double kw;
    double kg,kgamma;
    
    muNWA = muggH(sqrt_s)*BrHgagaRatio();
    
    if (!FlagKiLoop) {
    
        kt = computeKt();
        kb = computeKb();
        kc = computeKc();
        ks = computeKs();
        ku = computeKu();
        kd = computeKd();
        kta = computeKtau();
        km = computeKmu();
        ke = computeKe();
    
        kw = computeKW();
    
        kg = 0.0;
        kgamma = 0.0;
    
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
    } else {
        return muNWA;
    }
}
    
double HiggsKigen::muVBFHgaga(const double sqrt_s) const
{
    return muVBF(sqrt_s)*BrHgagaRatio();
}

double HiggsKigen::muZHgaga(const double sqrt_s) const
{
    return muZH(sqrt_s)*BrHgagaRatio();
}

double HiggsKigen::muWHgaga(const double sqrt_s) const
{
    return muWH(sqrt_s)*BrHgagaRatio();
}

double HiggsKigen::muVHgaga(const double sqrt_s) const
{
    return muVH(sqrt_s)*BrHgagaRatio();
}

double HiggsKigen::muttHgaga(const double sqrt_s) const
{
    return muttH(sqrt_s)*BrHgagaRatio();
}

double HiggsKigen::muggHZga(const double sqrt_s) const
{
    return muggH(sqrt_s)*BrHZgaRatio();
}

double HiggsKigen::muVBFHZga(const double sqrt_s) const
{
    return muVBF(sqrt_s)*BrHZgaRatio();
}

double HiggsKigen::muZHZga(const double sqrt_s) const
{
    return muZH(sqrt_s)*BrHZgaRatio();
}

double HiggsKigen::muWHZga(const double sqrt_s) const
{
    return muWH(sqrt_s)*BrHZgaRatio();
}

double HiggsKigen::muVHZga(const double sqrt_s) const
{
    return muVH(sqrt_s)*BrHZgaRatio();
}

double HiggsKigen::muttHZga(const double sqrt_s) const
{
    return muttH(sqrt_s)*BrHZgaRatio();
}

double HiggsKigen::muggHZZ(const double sqrt_s) const
{
    return muggH(sqrt_s)*BrHZZRatio();
}

double HiggsKigen::muVBFHZZ(const double sqrt_s) const
{
    return muVBF(sqrt_s)*BrHZZRatio();
}

double HiggsKigen::muZHZZ(const double sqrt_s) const
{
    return muZH(sqrt_s)*BrHZZRatio();
}

double HiggsKigen::muWHZZ(const double sqrt_s) const
{
    return muWH(sqrt_s)*BrHZZRatio();
}

double HiggsKigen::muVHZZ(const double sqrt_s) const
{
    return muVH(sqrt_s)*BrHZZRatio();
}

double HiggsKigen::muttHZZ(const double sqrt_s) const
{
    return muttH(sqrt_s)*BrHZZRatio();
}

double HiggsKigen::muggHWW(const double sqrt_s) const
{
    return muggH(sqrt_s)*BrHWWRatio();
}

double HiggsKigen::muVBFHWW(const double sqrt_s) const
{
    return muVBF(sqrt_s)*BrHWWRatio();
}

double HiggsKigen::muZHWW(const double sqrt_s) const
{
    return muZH(sqrt_s)*BrHWWRatio();
}

double HiggsKigen::muWHWW(const double sqrt_s) const
{
    return muWH(sqrt_s)*BrHWWRatio();
}

double HiggsKigen::muVHWW(const double sqrt_s) const
{
    return muVH(sqrt_s)*BrHWWRatio();
}

double HiggsKigen::muttHWW(const double sqrt_s) const
{
    return muttH(sqrt_s)*BrHWWRatio();
}

double HiggsKigen::muggHmumu(const double sqrt_s) const
{
    return muggH(sqrt_s)*BrHmumuRatio();
}

double HiggsKigen::muVBFHmumu(const double sqrt_s) const
{
    return muVBF(sqrt_s)*BrHmumuRatio();
}

double HiggsKigen::muZHmumu(const double sqrt_s) const
{
    return muZH(sqrt_s)*BrHmumuRatio();
}

double HiggsKigen::muWHmumu(const double sqrt_s) const
{
    return muWH(sqrt_s)*BrHmumuRatio();
}

double HiggsKigen::muVHmumu(const double sqrt_s) const
{
    return muVH(sqrt_s)*BrHmumuRatio();
}

double HiggsKigen::muttHmumu(const double sqrt_s) const
{
    return muttH(sqrt_s)*BrHmumuRatio();
}

double HiggsKigen::muggHtautau(const double sqrt_s) const
{
    return muggH(sqrt_s)*BrHtautauRatio();
}

double HiggsKigen::muVBFHtautau(const double sqrt_s) const
{
    return muVBF(sqrt_s)*BrHtautauRatio();
}

double HiggsKigen::muZHtautau(const double sqrt_s) const
{
    return muZH(sqrt_s)*BrHtautauRatio();
}

double HiggsKigen::muWHtautau(const double sqrt_s) const
{
    return muWH(sqrt_s)*BrHtautauRatio();
}

double HiggsKigen::muVHtautau(const double sqrt_s) const
{
    return muVH(sqrt_s)*BrHtautauRatio();
}

double HiggsKigen::muttHtautau(const double sqrt_s) const
{
    return muttH(sqrt_s)*BrHtautauRatio();
}

double HiggsKigen::muggHbb(const double sqrt_s) const
{
    return muggH(sqrt_s)*BrHbbRatio();
}

double HiggsKigen::muVBFHbb(const double sqrt_s) const
{
    return muVBF(sqrt_s)*BrHbbRatio();
}

double HiggsKigen::muZHbb(const double sqrt_s) const
{
    return muZH(sqrt_s)*BrHbbRatio();
}

double HiggsKigen::muWHbb(const double sqrt_s) const
{
    return muWH(sqrt_s)*BrHbbRatio();
}

double HiggsKigen::muVHbb(const double sqrt_s) const
{
    return muVH(sqrt_s)*BrHbbRatio();
}

double HiggsKigen::muttHbb(const double sqrt_s) const
{
    return muttH(sqrt_s)*BrHbbRatio();
}


double HiggsKigen::computeGammaTotalRatio() const
{
    double GammaRatioSMmodes;

    GammaRatioSMmodes = ((computeKg() * computeKg() * trueSM.computeBrHtogg()
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
        
    return (GammaRatioSMmodes/(1.0 - BrHinv - BrHexo));

}

double HiggsKigen::Br_H_exo() const
{
    return BrHexo;
};


double HiggsKigen::Br_H_inv() const
{
    return BrHinv;
}


double HiggsKigen::BrHvisRatio() const
{    
        return (1.0 - BrHinv);   
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

double HiggsKigen::computeKe() const
{
    return Ke;
}

double HiggsKigen::computeKmu() const
{
    return Kmu;
}

double HiggsKigen::computeKtau() const
{
    return Ktau;
}

double HiggsKigen::computeKu() const
{
    return Ku;
}

double HiggsKigen::computeKc() const
{
    return Kc;
}

double HiggsKigen::computeKt() const
{
    return Kt;
}

double HiggsKigen::computeKd() const
{
    return Kd;
}

double HiggsKigen::computeKs() const
{
    return Ks;
}

double HiggsKigen::computeKb() const
{
    return Kb;
}

//  USE THESE TO DEFINE KAPPAg,ETC

double HiggsKigen::computeKgLoop() const
{
    double Mt=trueSM.getQuarks(QCD::TOP).getMass();
    double Mb=trueSM.getQuarks(QCD::BOTTOM).getMass();
    double Mc=trueSM.getQuarks(QCD::CHARM).getMass();
    double TAUt=4.0*Mt*Mt/(mHl*mHl);
    double TAUb=4.0*Mb*Mb/(mHl*mHl);
    double TAUc=4.0*Mc*Mc/(mHl*mHl);
    double KgEff;

    KgEff = ( ( Kt * 0.5*TAUt*(1.0+(1.0-TAUt)*f_func(TAUt))
                     +Kb * 0.5*TAUb*(1.0+(1.0-TAUb)*f_func(TAUb)) 
                     +Kc * 0.5*TAUc*(1.0+(1.0-TAUc)*f_func(TAUc)) ) * (1.0+11.0*AlsMz/(4.0*M_PI))).abs2();

    KgEff = KgEff / ( ( 0.5*TAUt*(1.0+(1.0-TAUt)*f_func(TAUt))
                     + 0.5*TAUb*(1.0+(1.0-TAUb)*f_func(TAUb))
                     + 0.5*TAUc*(1.0+(1.0-TAUc)*f_func(TAUc)) ) * (1.0+11.0*AlsMz/(4.0*M_PI))).abs2();

    return (sqrt(KgEff));
}

double HiggsKigen::computeKZgaLoop() const
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
    double KZgaEff;

    KZgaEff = ((-Kt * 4.0*(0.5-4.0/3.0*sW2)*(Int1(TAUt,LAMt)-Int2(TAUt,LAMt)) * (1.0-AlsMz/M_PI)
                      +Kb * 2.0*(-0.5+2.0/3.0*sW2)*(Int1(TAUb,LAMb)-Int2(TAUb,LAMb))
                      -Kc * 4.0*(0.5-4.0/3.0*sW2)*(Int1(TAUc,LAMc)-Int2(TAUc,LAMc))            
                      +Ktau * 2.0*(-0.5+2.0*sW2)*(Int1(TAUtau,LAMtau)-Int2(TAUtau,LAMtau))
                      +Kmu * 2.0*(-0.5+2.0*sW2)*(Int1(TAUmu,LAMmu)-Int2(TAUmu,LAMmu)) )/sqrt(sW2*cW2)
                     -Kw * sqrt(cW2/sW2)*(4.0*(3.0-sW2/cW2)*Int2(TAUw,LAMw)
                            +((1.0+2.0/TAUw)*sW2/cW2-(5.0+2.0/TAUw))*Int1(TAUw,LAMw))).abs2();
    
    KZgaEff = KZgaEff / ((-4.0*(0.5-4.0/3.0*sW2)*(Int1(TAUt,LAMt)-Int2(TAUt,LAMt)) * (1.0-AlsMz/M_PI)
                      + 2.0*(-0.5+2.0/3.0*sW2)*(Int1(TAUb,LAMb)-Int2(TAUb,LAMb))
                      - 4.0*(0.5-4.0/3.0*sW2)*(Int1(TAUc,LAMc)-Int2(TAUc,LAMc))
                      + 2.0*(-0.5+2.0*sW2)*(Int1(TAUtau,LAMtau)-Int2(TAUtau,LAMtau))
                      + 2.0*(-0.5+2.0*sW2)*(Int1(TAUmu,LAMmu)-Int2(TAUmu,LAMmu)) )/sqrt(sW2*cW2)
                      - sqrt(cW2/sW2)*(4.0*(3.0-sW2/cW2)*Int2(TAUw,LAMw)
                            +((1.0+2.0/TAUw)*sW2/cW2-(5.0+2.0/TAUw))*Int1(TAUw,LAMw))).abs2();
    
    return (sqrt(KZgaEff));
}

double HiggsKigen::computeKgagaLoop() const
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
    double KgagaEff;

    KgagaEff = ( Kt * (8./3.)*TAUt*(1.+(1.-TAUt)*f_func(TAUt)) * (1.0-AlsMz/M_PI)
                     +Kb * (2./3.)*TAUb*(1.+(1.-TAUb)*f_func(TAUb))
                     +Kc * (8./3.)*TAUc*(1.+(1.-TAUc)*f_func(TAUc))
                     +Ktau * 2.0*TAUtau*(1.+(1.-TAUtau)*f_func(TAUtau))
                     +Kmu * 2.0*TAUmu*(1.+(1.-TAUmu)*f_func(TAUmu))
                     -Kw * (2.0+3.0*TAUw+3.0*TAUw*(2.0-TAUw)*f_func(TAUw)) ).abs2();
    
    KgagaEff = KgagaEff / ( (8./3.)*TAUt*(1.+(1.-TAUt)*f_func(TAUt)) * (1.0-AlsMz/M_PI)
                     + (2./3.)*TAUb*(1.+(1.-TAUb)*f_func(TAUb))
                     + (8./3.)*TAUc*(1.+(1.-TAUc)*f_func(TAUc))
                     + 2.0*TAUtau*(1.+(1.-TAUtau)*f_func(TAUtau))
                     + 2.0*TAUmu*(1.+(1.-TAUmu)*f_func(TAUmu))
                     - (2.0+3.0*TAUw+3.0*TAUw*(2.0-TAUw)*f_func(TAUw)) ).abs2();
    
    return (sqrt(KgagaEff));
}

//////////////////////////////////////////////////////////////////////

gslpp::complex HiggsKigen::f_func(const double x) const{
    if(x<1) {
    gslpp::complex z = -gslpp::complex::i()*M_PI;
    return -pow(log((1+sqrt(1-x))/(1-sqrt(1-x)))+z,2)/4.0;
    }
    else {
        return pow(asin(sqrt(1.0/x)),2);
    }
}

gslpp::complex HiggsKigen::g_func(const double x) const{
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

gslpp::complex HiggsKigen::Int1(const double tau, const double lambda) const{
    return tau*lambda/(tau-lambda)/2.0+tau*tau*lambda*lambda/((tau-lambda)
           *(tau-lambda))/2.0*(f_func(tau)-f_func(lambda))+tau*tau*lambda/((tau-lambda)
           *(tau-lambda))*(g_func(tau)-g_func(lambda));
}

gslpp::complex HiggsKigen::Int2(const double tau, const double lambda) const{
    return -tau*lambda/(tau-lambda)/2.0*(f_func(tau)-f_func(lambda));
}


