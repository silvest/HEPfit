/*
 * Copyright (C) 2016 HEPfit Collaboration
 *
 * For the licensing terms see doc/COPYING.
 */

#include "SigmaBR.h"

const std::string SigmaBR::SigmaBRVars[NSigmaBRVars]
        = {"ggH", "VBF", "WH", "ZH", "ttH", 
        "eeZH240","eeZH250","eeZH350","eeZH365","eeZH380","eeZH500","eeZH1000","eeZH1500","eeZH3000",
        "eeHvv240","eeHvv250","eeHvv350","eeHvv365","eeHvv380","eeHvv500","eeHvv1000","eeHvv1500","eeHvv3000",
        "eeHee240","eeHee250","eeHee350","eeHee365","eeHee380","eeHee500","eeHee1000","eeHee1500","eeHee3000",
        "eettH500","eettH1000","eettH1500","eettH3000",
        "BrHggRatio", "BrHWWRatio", "BrHZZRatio","BrHZgaRatio", "BrHgagaRatio", "BrHtautauRatio", "BrHmumuRatio", "BrHccRatio", "BrHbbRatio"};

SigmaBR::SigmaBR()
: NPbase()
{
    
    ModelParamMap.insert(std::make_pair("ggH", std::cref(ggh)));
    ModelParamMap.insert(std::make_pair("VBF", std::cref(vbf)));
    ModelParamMap.insert(std::make_pair("WH", std::cref(wh)));
    ModelParamMap.insert(std::make_pair("ZH", std::cref(zh)));
    ModelParamMap.insert(std::make_pair("ttH", std::cref(tth)));
    ModelParamMap.insert(std::make_pair("eeZH240", std::cref(eeZH240)));
    ModelParamMap.insert(std::make_pair("eeZH250", std::cref(eeZH250)));
    ModelParamMap.insert(std::make_pair("eeZH350", std::cref(eeZH350)));
    ModelParamMap.insert(std::make_pair("eeZH365", std::cref(eeZH365)));
    ModelParamMap.insert(std::make_pair("eeZH380", std::cref(eeZH380)));    
    ModelParamMap.insert(std::make_pair("eeZH500", std::cref(eeZH500)));
    ModelParamMap.insert(std::make_pair("eeZH1000", std::cref(eeZH1000)));
    ModelParamMap.insert(std::make_pair("eeZH1500", std::cref(eeZH1500)));
    ModelParamMap.insert(std::make_pair("eeZH3000", std::cref(eeZH3000)));
    ModelParamMap.insert(std::make_pair("eeHvv240", std::cref(eeHvv240)));
    ModelParamMap.insert(std::make_pair("eeHvv250", std::cref(eeHvv250)));
    ModelParamMap.insert(std::make_pair("eeHvv350", std::cref(eeHvv350)));
    ModelParamMap.insert(std::make_pair("eeHvv365", std::cref(eeHvv365)));
    ModelParamMap.insert(std::make_pair("eeHvv380", std::cref(eeHvv380)));    
    ModelParamMap.insert(std::make_pair("eeHvv500", std::cref(eeHvv500)));
    ModelParamMap.insert(std::make_pair("eeHvv1000", std::cref(eeHvv1000)));
    ModelParamMap.insert(std::make_pair("eeHvv1500", std::cref(eeHvv1500)));
    ModelParamMap.insert(std::make_pair("eeHvv3000", std::cref(eeHvv3000)));   
    ModelParamMap.insert(std::make_pair("eeHee240", std::cref(eeHee240)));
    ModelParamMap.insert(std::make_pair("eeHee250", std::cref(eeHee250)));
    ModelParamMap.insert(std::make_pair("eeHee350", std::cref(eeHee350)));
    ModelParamMap.insert(std::make_pair("eeHee365", std::cref(eeHee365)));
    ModelParamMap.insert(std::make_pair("eeHee380", std::cref(eeHee380)));    
    ModelParamMap.insert(std::make_pair("eeHee500", std::cref(eeHee500)));
    ModelParamMap.insert(std::make_pair("eeHee1000", std::cref(eeHee1000)));
    ModelParamMap.insert(std::make_pair("eeHee1500", std::cref(eeHee1500)));
    ModelParamMap.insert(std::make_pair("eeHee3000", std::cref(eeHee3000))); 
    ModelParamMap.insert(std::make_pair("eettH500", std::cref(eettH500)));
    ModelParamMap.insert(std::make_pair("eettH1000", std::cref(eettH1000)));
    ModelParamMap.insert(std::make_pair("eettH1500", std::cref(eettH1500)));
    ModelParamMap.insert(std::make_pair("eettH3000", std::cref(eeHee3000)));
    ModelParamMap.insert(std::make_pair("BrHggRatio", std::cref(brhggratio)));
    ModelParamMap.insert(std::make_pair("BrHWWRatio", std::cref(brhwwratio)));
    ModelParamMap.insert(std::make_pair("BrHZZRatio", std::cref(brhzzratio)));
    ModelParamMap.insert(std::make_pair("BrHZgaRatio", std::cref(brhzgaratio)));
    ModelParamMap.insert(std::make_pair("BrHgagaRatio", std::cref(brhgagaratio)));
    ModelParamMap.insert(std::make_pair("BrHtautauRatio", std::cref(brhtautauratio)));
    ModelParamMap.insert(std::make_pair("BrHmumuRatio", std::cref(brhmumuratio)));
    ModelParamMap.insert(std::make_pair("BrHccRatio", std::cref(brhccratio)));
    ModelParamMap.insert(std::make_pair("BrHbbRatio", std::cref(brhbbratio)));
}

void SigmaBR::setParameter(const std::string name, const double& value)
{
    if (name.compare("ggH") == 0)
        ggh = value;
    else if (name.compare("VBF") == 0)
        vbf = value;
    else if (name.compare("WH") == 0)
        wh = value;
    else if (name.compare("ZH") == 0)
        zh = value;
    else if (name.compare("ttH") == 0)
        tth = value;
    else if (name.compare("eeZH240") == 0)
        eeZH240 = value; 
    else if (name.compare("eeZH250") == 0)
        eeZH250 = value; 
    else if (name.compare("eeZH350") == 0)
        eeZH350 = value;     
    else if (name.compare("eeZH365") == 0)
        eeZH365 = value; 
    else if (name.compare("eeZH380") == 0)
        eeZH380 = value; 
    else if (name.compare("eeZH500") == 0)
        eeZH500 = value; 
    else if (name.compare("eeZH1000") == 0)
        eeZH1000 = value; 
    else if (name.compare("eeZH1500") == 0)
        eeZH1500 = value; 
    else if (name.compare("eeZH3000") == 0)
        eeZH3000 = value; 
    else if (name.compare("eeHvv240") == 0)
        eeHvv240 = value; 
    else if (name.compare("eeHvv250") == 0)
        eeHvv250 = value; 
    else if (name.compare("eeHvv350") == 0)
        eeHvv350 = value;     
    else if (name.compare("eeHvv365") == 0)
        eeHvv365 = value; 
    else if (name.compare("eeHvv380") == 0)
        eeHvv380 = value; 
    else if (name.compare("eeHvv500") == 0)
        eeHvv500 = value; 
    else if (name.compare("eeHvv1000") == 0)
        eeHvv1000 = value; 
    else if (name.compare("eeHvv1500") == 0)
        eeHvv1500 = value; 
    else if (name.compare("eeHvv3000") == 0)
        eeHvv3000 = value;
    else if (name.compare("eeHee240") == 0)
        eeHee240 = value; 
    else if (name.compare("eeHee250") == 0)
        eeHee250 = value; 
    else if (name.compare("eeHee350") == 0)
        eeHee350 = value;     
    else if (name.compare("eeHee365") == 0)
        eeHee365 = value; 
    else if (name.compare("eeHee380") == 0)
        eeHee380 = value; 
    else if (name.compare("eeHee500") == 0)
        eeHee500 = value; 
    else if (name.compare("eeHee1000") == 0)
        eeHee1000 = value; 
    else if (name.compare("eeHee1500") == 0)
        eeHee1500 = value; 
    else if (name.compare("eeHee3000") == 0)
        eeHee3000 = value;
    else if (name.compare("eettH500") == 0)
        eettH500 = value; 
    else if (name.compare("eettH1000") == 0)
        eettH1000 = value; 
    else if (name.compare("eettH1500") == 0)
        eettH1500 = value; 
    else if (name.compare("eettH3000") == 0)
        eettH3000 = value; 
    else if (name.compare("BrHggRatio") == 0)
        brhggratio = value;
    else if (name.compare("BrHWWRatio") == 0)
        brhwwratio = value;
    else if (name.compare("BrHZZRatio") == 0)
        brhzzratio = value;
    else if (name.compare("BrHZgaRatio") == 0)
        brhzgaratio = value;
    else if (name.compare("BrHgagaRatio") == 0)
        brhgagaratio = value;
    else if (name.compare("BrHtautauRatio") == 0)
        brhtautauratio = value;
    else if (name.compare("BrHmumuRatio") == 0)
        brhmumuratio = value;
    else if (name.compare("BrHccRatio") == 0)
        brhccratio = value;
    else if (name.compare("BrHbbRatio") == 0)
        brhbbratio = value;
    else
        NPbase::setParameter(name, value);
}

bool SigmaBR::CheckParameters(const std::map<std::string, double>& DPars)
{
        for (int i = 0; i < NSigmaBRVars; i++) 
            if (DPars.find(SigmaBRVars[i]) == DPars.end()) {
                std::cout << "ERROR: Missing mandatory SigmaBR parameter"
                        << SigmaBRVars[i] << std::endl;
                raiseMissingModelParameterCount();
                addMissingModelParameter(SigmaBRVars[i]);
            }
    return (NPbase::CheckParameters(DPars));
}

double SigmaBR::muggH(const double sqrt_s) const
{    
    return ggh;
}

double SigmaBR::muVBF(const double sqrt_s) const
{
    return vbf;
}

double SigmaBR::muWH(const double sqrt_s) const
{
    return wh;
}

double SigmaBR::muZH(const double sqrt_s) const
{
    return zh;
}

double SigmaBR::muVH(const double sqrt_s) const
{
    double sigmaWH_SM = computeSigmaWH(sqrt_s);
    double sigmaZH_SM = computeSigmaZH(sqrt_s);
    double sigmaWH = muWH(sqrt_s) * sigmaWH_SM;
    double sigmaZH = muZH(sqrt_s) * sigmaZH_SM;
    double mu = ((sigmaWH + sigmaZH) / (sigmaWH_SM + sigmaZH_SM));
        
    return mu;
}

double SigmaBR::muVBFpVH(const double sqrt_s) const
{
    double sigmaWH_SM = computeSigmaWH(sqrt_s);
    double sigmaZH_SM = computeSigmaZH(sqrt_s);
    double sigmaVBF_SM = computeSigmaVBF(sqrt_s);
    double sigmaWH = muWH(sqrt_s) * sigmaWH_SM;
    double sigmaZH = muZH(sqrt_s) * sigmaZH_SM;
    double sigmaVBF = muVBF(sqrt_s) * sigmaVBF_SM;
    double mu = ((sigmaWH + sigmaZH + sigmaVBF) / (sigmaWH_SM + sigmaZH_SM + sigmaVBF_SM));
        
    return mu;
}

double SigmaBR::muttH(const double sqrt_s) const
{
    return tth;
}

double SigmaBR::muggHpttH(const double sqrt_s) const
{
    double sigmaggH_SM = computeSigmaggH(sqrt_s);
    double sigmattH_SM = computeSigmattH(sqrt_s);
    double sigmaggH = muggH(sqrt_s) * sigmaggH_SM;
    double sigmattH = muttH(sqrt_s) * sigmattH_SM;

    double mu = ((sigmaggH + sigmattH) / (sigmaggH_SM + sigmattH_SM));
        
    return mu;
} 

double SigmaBR::mueeZH(const double sqrt_s) const
{
    
    if (sqrt_s == 0.240) {
        
        return eeZH240;
        
    } else if (sqrt_s == 0.250) {
        
        return eeZH250;
        
    } else if (sqrt_s == 0.350) {

        return eeZH350;        
        
    } else if (sqrt_s == 0.365) {
        
        return eeZH365;
        
    } else if (sqrt_s == 0.380) {
        
        return eeZH380;
        
    } else if (sqrt_s == 0.500) {
        
        return eeZH500;
        
    } else if (sqrt_s == 1.0) {
        
        return eeZH1000;
        
    } else if (sqrt_s == 1.5) {
        
        return eeZH1500;
        
    } else if (sqrt_s == 3.0) {
        
        return eeZH3000;
        
    } else
        throw std::runtime_error("Bad argument in SigmaBR::mueeZH()");
    
    return 1.0;

}

double SigmaBR::mueeZHPol(const double sqrt_s, const double Pol_em, const double Pol_ep) const
{
    return mueeZH(sqrt_s);
}


double SigmaBR::mueeHvv(const double sqrt_s) const
{
    if (sqrt_s == 0.240) {
        
        return eeHvv240;
        
    } else if (sqrt_s == 0.250) {
        
        return eeHvv250;
        
    } else if (sqrt_s == 0.350) {

        return eeHvv350;        
        
    } else if (sqrt_s == 0.365) {
        
        return eeHvv365;
        
    } else if (sqrt_s == 0.380) {
        
        return eeHvv380;
        
    } else if (sqrt_s == 0.500) {
        
        return eeHvv500;
        
    } else if (sqrt_s == 1.0) {
        
        return eeHvv1000;
        
    } else if (sqrt_s == 1.5) {
        
        return eeHvv1500;
        
    } else if (sqrt_s == 3.0) {
        
        return eeHvv3000;
        
    } else
        throw std::runtime_error("Bad argument in SigmaBR::mueeHvv()");
    
    return 1.0;
}

double SigmaBR::mueeHvvPol(const double sqrt_s, const double Pol_em, const double Pol_ep) const
{
    return mueeHvv(sqrt_s);
}

double SigmaBR::mueeZBF(const double sqrt_s) const
{
    if (sqrt_s == 0.240) {
        
        return eeHee240;
        
    } else if (sqrt_s == 0.250) {
        
        return eeHee250;
        
    } else if (sqrt_s == 0.350) {

        return eeHee350;        
        
    } else if (sqrt_s == 0.365) {
        
        return eeHee365;
        
    } else if (sqrt_s == 0.380) {
        
        return eeHee380;
        
    } else if (sqrt_s == 0.500) {
        
        return eeHee500;
        
    } else if (sqrt_s == 1.0) {
        
        return eeHee1000;
        
    } else if (sqrt_s == 1.5) {
        
        return eeHee1500;
        
    } else if (sqrt_s == 3.0) {
        
        return eeHee3000;
        
    } else
        throw std::runtime_error("Bad argument in SigmaBR::mueeHee()");
    
    return 1.0;
}

double SigmaBR::mueeZBFPol(const double sqrt_s, const double Pol_em, const double Pol_ep) const
{
    return mueeZBF(sqrt_s);
}

double SigmaBR::mueettH(const double sqrt_s) const
{
    if (sqrt_s == 0.500) {
        
        return eettH500;
        
    } else if (sqrt_s == 1.0) {
        
        return eettH1000;
        
    } else if (sqrt_s == 1.5) {
        
        return eettH1500;
        
    } else if (sqrt_s == 3.0) {
        
        return eettH3000;
        
    } else
        throw std::runtime_error("Bad argument in SigmaBR::mueettH()");
    
    return 1.0;
}

double SigmaBR::mueettHPol(const double sqrt_s, const double Pol_em, const double Pol_ep) const
{
    return mueettH(sqrt_s);
}

double SigmaBR::BrHggRatio() const
{
    return brhggratio;
}

double SigmaBR::BrHWWRatio() const
{
    return brhwwratio;
}

double SigmaBR::BrHWW2l2vRatio() const
{
    return BrHWWRatio();
}

double SigmaBR::BrHZZRatio() const
{
    return brhzzratio;
}

double SigmaBR::BrHZZ4lRatio() const
{
    return BrHZZRatio();
}

double SigmaBR::BrHZZ4eRatio() const
{
    return BrHZZRatio();
}

double SigmaBR::BrHZZ2e2muRatio() const
{
    return BrHZZRatio();
}

double SigmaBR::BrHZZ4muRatio() const
{
    return BrHZZRatio();
}

double SigmaBR::BrHZgaRatio() const
{
    return brhzgaratio;
}

double SigmaBR::BrHZgallRatio() const
{
    return BrHZgaRatio();
}

double SigmaBR::BrHZgaeeRatio() const
{
    return BrHZgaRatio();
}

double SigmaBR::BrHZgamumuRatio() const
{
    return BrHZgaRatio();
}

double SigmaBR::BrHgagaRatio() const
{
    return brhgagaratio;
}

double SigmaBR::BrHtautauRatio() const
{
    return brhtautauratio;
}

double SigmaBR::BrHmumuRatio() const
{
    return brhmumuratio;
}

double SigmaBR::BrHccRatio() const
{
    return brhccratio;
}

double SigmaBR::BrHbbRatio() const
{
    return brhbbratio;
}

double SigmaBR::muggHgaga(const double sqrt_s) const
{
    return muggH(sqrt_s)*BrHgagaRatio();
}
    
double SigmaBR::muVBFHgaga(const double sqrt_s) const
{
    return muVBF(sqrt_s)*BrHgagaRatio();
}

double SigmaBR::muZHgaga(const double sqrt_s) const
{
    return muZH(sqrt_s)*BrHgagaRatio();
}

double SigmaBR::muWHgaga(const double sqrt_s) const
{
    return muWH(sqrt_s)*BrHgagaRatio();
}

double SigmaBR::muVHgaga(const double sqrt_s) const
{
    return muVH(sqrt_s)*BrHgagaRatio();
}

double SigmaBR::muttHgaga(const double sqrt_s) const
{
    return muttH(sqrt_s)*BrHgagaRatio();
}

double SigmaBR::muggHZga(const double sqrt_s) const
{
    return muggH(sqrt_s)*BrHZgaRatio();
}

double SigmaBR::muVBFHZga(const double sqrt_s) const
{
    return muVBF(sqrt_s)*BrHZgaRatio();
}

double SigmaBR::muZHZga(const double sqrt_s) const
{
    return muZH(sqrt_s)*BrHZgaRatio();
}

double SigmaBR::muWHZga(const double sqrt_s) const
{
    return muWH(sqrt_s)*BrHZgaRatio();
}

double SigmaBR::muVHZga(const double sqrt_s) const
{
    return muVH(sqrt_s)*BrHZgaRatio();
}

double SigmaBR::muttHZga(const double sqrt_s) const
{
    return muttH(sqrt_s)*BrHZgaRatio();
}

double SigmaBR::muggHZZ(const double sqrt_s) const
{
    return muggH(sqrt_s)*BrHZZRatio();
}

double SigmaBR::muVBFHZZ(const double sqrt_s) const
{
    return muVBF(sqrt_s)*BrHZZRatio();
}

double SigmaBR::muZHZZ(const double sqrt_s) const
{
    return muZH(sqrt_s)*BrHZZRatio();
}

double SigmaBR::muWHZZ(const double sqrt_s) const
{
    return muWH(sqrt_s)*BrHZZRatio();
}

double SigmaBR::muVHZZ(const double sqrt_s) const
{
    return muVH(sqrt_s)*BrHZZRatio();
}

double SigmaBR::muttHZZ(const double sqrt_s) const
{
    return muttH(sqrt_s)*BrHZZRatio();
}

double SigmaBR::muggHZZ4l(const double sqrt_s) const
{
    return muggH(sqrt_s)*BrHZZ4lRatio();
}

double SigmaBR::muVBFHZZ4l(const double sqrt_s) const
{
    return muVBF(sqrt_s)*BrHZZ4lRatio();
}

double SigmaBR::muZHZZ4l(const double sqrt_s) const
{
    return muZH(sqrt_s)*BrHZZ4lRatio();
}

double SigmaBR::muWHZZ4l(const double sqrt_s) const
{
    return muWH(sqrt_s)*BrHZZ4lRatio();
}

double SigmaBR::muVHZZ4l(const double sqrt_s) const
{
    return muVH(sqrt_s)*BrHZZ4lRatio();
}

double SigmaBR::muttHZZ4l(const double sqrt_s) const
{
    return muttH(sqrt_s)*BrHZZ4lRatio();
}

double SigmaBR::muggHWW(const double sqrt_s) const
{
    return muggH(sqrt_s)*BrHWWRatio();
}

double SigmaBR::muVBFHWW(const double sqrt_s) const
{
    return muVBF(sqrt_s)*BrHWWRatio();
}

double SigmaBR::muZHWW(const double sqrt_s) const
{
    return muZH(sqrt_s)*BrHWWRatio();
}

double SigmaBR::muWHWW(const double sqrt_s) const
{
    return muWH(sqrt_s)*BrHWWRatio();
}

double SigmaBR::muVHWW(const double sqrt_s) const
{
    return muVH(sqrt_s)*BrHWWRatio();
}

double SigmaBR::muttHWW(const double sqrt_s) const
{
    return muttH(sqrt_s)*BrHWWRatio();
}

double SigmaBR::muggHWW2l2v(const double sqrt_s) const
{
    return muggH(sqrt_s)*BrHWW2l2vRatio();
}

double SigmaBR::muVBFHWW2l2v(const double sqrt_s) const
{
    return muVBF(sqrt_s)*BrHWW2l2vRatio();
}

double SigmaBR::muZHWW2l2v(const double sqrt_s) const
{
    return muZH(sqrt_s)*BrHWW2l2vRatio();
}

double SigmaBR::muWHWW2l2v(const double sqrt_s) const
{
    return muWH(sqrt_s)*BrHWW2l2vRatio();
}

double SigmaBR::muVHWW2l2v(const double sqrt_s) const
{
    return muVH(sqrt_s)*BrHWW2l2vRatio();
}

double SigmaBR::muttHWW2l2v(const double sqrt_s) const
{
    return muttH(sqrt_s)*BrHWW2l2vRatio();
}

double SigmaBR::muggHmumu(const double sqrt_s) const
{
    return muggH(sqrt_s)*BrHmumuRatio();
}

double SigmaBR::muVBFHmumu(const double sqrt_s) const
{
    return muVBF(sqrt_s)*BrHmumuRatio();
}

double SigmaBR::muZHmumu(const double sqrt_s) const
{
    return muZH(sqrt_s)*BrHmumuRatio();
}

double SigmaBR::muWHmumu(const double sqrt_s) const
{
    return muWH(sqrt_s)*BrHmumuRatio();
}

double SigmaBR::muVHmumu(const double sqrt_s) const
{
    return muVH(sqrt_s)*BrHmumuRatio();
}

double SigmaBR::muttHmumu(const double sqrt_s) const
{
    return muttH(sqrt_s)*BrHmumuRatio();
}

double SigmaBR::muggHtautau(const double sqrt_s) const
{
    return muggH(sqrt_s)*BrHtautauRatio();
}

double SigmaBR::muVBFHtautau(const double sqrt_s) const
{
    return muVBF(sqrt_s)*BrHtautauRatio();
}

double SigmaBR::muZHtautau(const double sqrt_s) const
{
    return muZH(sqrt_s)*BrHtautauRatio();
}

double SigmaBR::muWHtautau(const double sqrt_s) const
{
    return muWH(sqrt_s)*BrHtautauRatio();
}

double SigmaBR::muVHtautau(const double sqrt_s) const
{
    return muVH(sqrt_s)*BrHtautauRatio();
}

double SigmaBR::muttHtautau(const double sqrt_s) const
{
    return muttH(sqrt_s)*BrHtautauRatio();
}

double SigmaBR::muggHbb(const double sqrt_s) const
{
    return muggH(sqrt_s)*BrHbbRatio();
}

double SigmaBR::muVBFHbb(const double sqrt_s) const
{
    return muVBF(sqrt_s)*BrHbbRatio();
}

double SigmaBR::muZHbb(const double sqrt_s) const
{
    return muZH(sqrt_s)*BrHbbRatio();
}

double SigmaBR::muWHbb(const double sqrt_s) const
{
    return muWH(sqrt_s)*BrHbbRatio();
}

double SigmaBR::muVHbb(const double sqrt_s) const
{
    return muVH(sqrt_s)*BrHbbRatio();
}

double SigmaBR::muttHbb(const double sqrt_s) const
{
    return muttH(sqrt_s)*BrHbbRatio();
}

////////////////////////////////////////////////////////////////////////////////////////////
//-----------------------------------------------------------------------------------------
//-- Special Hadron collider signal strengths with separate full TH unc U(prod x decay) ---
//-----------------------------------------------------------------------------------------
//////////////////////////////////////////////////////////////////////////////////////////// 

double SigmaBR::muTHUggHgaga(const double sqrt_s) const
{
    return muggH(sqrt_s)*BrHgagaRatio();
}
   
double SigmaBR::muTHUVBFHgaga(const double sqrt_s) const
{
    return muVBF(sqrt_s)*BrHgagaRatio();
}

double SigmaBR::muTHUZHgaga(const double sqrt_s) const
{
    return muZH(sqrt_s)*BrHgagaRatio();
}

double SigmaBR::muTHUWHgaga(const double sqrt_s) const
{
    return muWH(sqrt_s)*BrHgagaRatio();
}

double SigmaBR::muTHUVHgaga(const double sqrt_s) const
{
    return muVH(sqrt_s)*BrHgagaRatio();
}

double SigmaBR::muTHUttHgaga(const double sqrt_s) const
{
    return muttH(sqrt_s)*BrHgagaRatio();
}

double SigmaBR::muTHUggHZga(const double sqrt_s) const
{
    return muggH(sqrt_s)*BrHZgaRatio();
}

double SigmaBR::muTHUVBFHZga(const double sqrt_s) const
{
    return muVBF(sqrt_s)*BrHZgaRatio();
}

double SigmaBR::muTHUZHZga(const double sqrt_s) const
{
    return muZH(sqrt_s)*BrHZgaRatio();
}

double SigmaBR::muTHUWHZga(const double sqrt_s) const
{
    return muWH(sqrt_s)*BrHZgaRatio();
}

double SigmaBR::muTHUVHZga(const double sqrt_s) const
{
    return muVH(sqrt_s)*BrHZgaRatio();
}

double SigmaBR::muTHUttHZga(const double sqrt_s) const
{
    return muttH(sqrt_s)*BrHZgaRatio();
}

double SigmaBR::muTHUggHZZ(const double sqrt_s) const
{
    return muggH(sqrt_s)*BrHZZRatio();
}

double SigmaBR::muTHUVBFHZZ(const double sqrt_s) const
{
    return muVBF(sqrt_s)*BrHZZRatio();
}

double SigmaBR::muTHUZHZZ(const double sqrt_s) const
{
    return muZH(sqrt_s)*BrHZZRatio();
}

double SigmaBR::muTHUWHZZ(const double sqrt_s) const
{
    return muWH(sqrt_s)*BrHZZRatio();
}

double SigmaBR::muTHUVHZZ(const double sqrt_s) const
{
    return muVH(sqrt_s)*BrHZZRatio();
}

double SigmaBR::muTHUttHZZ(const double sqrt_s) const
{
    return muttH(sqrt_s)*BrHZZRatio();
}

double SigmaBR::muTHUggHZZ4l(const double sqrt_s) const
{
    return muggH(sqrt_s)*BrHZZ4lRatio();
}

double SigmaBR::muTHUVBFHZZ4l(const double sqrt_s) const
{
    return muVBF(sqrt_s)*BrHZZ4lRatio();
}

double SigmaBR::muTHUZHZZ4l(const double sqrt_s) const
{
    return muZH(sqrt_s)*BrHZZ4lRatio();
}

double SigmaBR::muTHUWHZZ4l(const double sqrt_s) const
{
    return muWH(sqrt_s)*BrHZZ4lRatio();
}

double SigmaBR::muTHUVHZZ4l(const double sqrt_s) const
{
    return muVH(sqrt_s)*BrHZZ4lRatio();
}

double SigmaBR::muTHUttHZZ4l(const double sqrt_s) const
{
    return muttH(sqrt_s)*BrHZZ4lRatio();
}

double SigmaBR::muTHUggHWW(const double sqrt_s) const
{
    return muggH(sqrt_s)*BrHWWRatio();
}

double SigmaBR::muTHUVBFHWW(const double sqrt_s) const
{
    return muVBF(sqrt_s)*BrHWWRatio();
}

double SigmaBR::muTHUZHWW(const double sqrt_s) const
{
    return muZH(sqrt_s)*BrHWWRatio();
}

double SigmaBR::muTHUWHWW(const double sqrt_s) const
{
    return muWH(sqrt_s)*BrHWWRatio();
}

double SigmaBR::muTHUVHWW(const double sqrt_s) const
{
    return muVH(sqrt_s)*BrHWWRatio();
}

double SigmaBR::muTHUttHWW(const double sqrt_s) const
{
    return muttH(sqrt_s)*BrHWWRatio();
}

double SigmaBR::muTHUggHWW2l2v(const double sqrt_s) const
{
    return muggH(sqrt_s)*BrHWW2l2vRatio();
}

double SigmaBR::muTHUVBFHWW2l2v(const double sqrt_s) const
{
    return muVBF(sqrt_s)*BrHWW2l2vRatio();
}

double SigmaBR::muTHUZHWW2l2v(const double sqrt_s) const
{
    return muZH(sqrt_s)*BrHWW2l2vRatio();
}

double SigmaBR::muTHUWHWW2l2v(const double sqrt_s) const
{
    return muWH(sqrt_s)*BrHWW2l2vRatio();
}

double SigmaBR::muTHUVHWW2l2v(const double sqrt_s) const
{
    return muVH(sqrt_s)*BrHWW2l2vRatio();
}

double SigmaBR::muTHUttHWW2l2v(const double sqrt_s) const
{
    return muttH(sqrt_s)*BrHWW2l2vRatio();
}

double SigmaBR::muTHUggHmumu(const double sqrt_s) const
{
    return muggH(sqrt_s)*BrHmumuRatio();
}

double SigmaBR::muTHUVBFHmumu(const double sqrt_s) const
{
    return muVBF(sqrt_s)*BrHmumuRatio();
}

double SigmaBR::muTHUZHmumu(const double sqrt_s) const
{
    return muZH(sqrt_s)*BrHmumuRatio();
}

double SigmaBR::muTHUWHmumu(const double sqrt_s) const
{
    return muWH(sqrt_s)*BrHmumuRatio();
}

double SigmaBR::muTHUVHmumu(const double sqrt_s) const
{
    return muVH(sqrt_s)*BrHmumuRatio();
}

double SigmaBR::muTHUttHmumu(const double sqrt_s) const
{
    return muttH(sqrt_s)*BrHmumuRatio();
}

double SigmaBR::muTHUggHtautau(const double sqrt_s) const
{
    return muggH(sqrt_s)*BrHtautauRatio();
}

double SigmaBR::muTHUVBFHtautau(const double sqrt_s) const
{
    return muVBF(sqrt_s)*BrHtautauRatio();
}

double SigmaBR::muTHUZHtautau(const double sqrt_s) const
{
    return muZH(sqrt_s)*BrHtautauRatio();
}

double SigmaBR::muTHUWHtautau(const double sqrt_s) const
{
    return muWH(sqrt_s)*BrHtautauRatio();
}

double SigmaBR::muTHUVHtautau(const double sqrt_s) const
{
    return muVH(sqrt_s)*BrHtautauRatio();
}

double SigmaBR::muTHUttHtautau(const double sqrt_s) const
{
    return muttH(sqrt_s)*BrHtautauRatio();
}

double SigmaBR::muTHUggHbb(const double sqrt_s) const
{
    return muggH(sqrt_s)*BrHbbRatio();
}

double SigmaBR::muTHUVBFHbb(const double sqrt_s) const
{
    return muVBF(sqrt_s)*BrHbbRatio();
}

double SigmaBR::muTHUZHbb(const double sqrt_s) const
{
    return muZH(sqrt_s)*BrHbbRatio();
}

double SigmaBR::muTHUWHbb(const double sqrt_s) const
{
    return muWH(sqrt_s)*BrHbbRatio();
}

double SigmaBR::muTHUVHbb(const double sqrt_s) const
{
    return muVH(sqrt_s)*BrHbbRatio();
}

double SigmaBR::muTHUttHbb(const double sqrt_s) const
{
    return muttH(sqrt_s)*BrHbbRatio();
}

double SigmaBR::muTHUVBFBRinv(const double sqrt_s) const
{
    return muVBF(sqrt_s)*Br_H_inv();
}

double SigmaBR::muTHUVBFHinv(const double sqrt_s) const
{
    return muVBF(sqrt_s)*BrHtoinvRatio();
}

double SigmaBR::muTHUVHBRinv(const double sqrt_s) const
{
    return muVH(sqrt_s)*Br_H_inv();
}

double SigmaBR::muTHUVHinv(const double sqrt_s) const
{
    return muVH(sqrt_s)*BrHtoinvRatio();
}

double SigmaBR::muTHUggHZZ4mu(const double sqrt_s) const
{
    return muggH(sqrt_s)*BrHZZ4muRatio();
}

double SigmaBR::muTHUggHZgamumu(const double sqrt_s) const
{
    return muggH(sqrt_s)*BrHZgamumuRatio();
}
