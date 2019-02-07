/*
 * Copyright (C) 2016 HEPfit Collaboration
 *
 * For the licensing terms see doc/COPYING.
 */

#include "SigmaBR.h"

const std::string SigmaBR::SigmaBRVars[NSigmaBRVars]
        = {"ggH", "VBF", "WH", "ZH", "ttH", "BrHggRatio", "BrHWWRatio", "BrHZZRatio",
    "BrHZgaRatio", "BrHgagaRatio", "BrHtautauRatio", "BrHccRatio", "BrHbbRatio"};

SigmaBR::SigmaBR()
: NPbase()
{
    
    ModelParamMap.insert(std::make_pair("ggH", std::cref(ggh)));
    ModelParamMap.insert(std::make_pair("VBF", std::cref(vbf)));
    ModelParamMap.insert(std::make_pair("WH", std::cref(wh)));
    ModelParamMap.insert(std::make_pair("ZH", std::cref(zh)));
    ModelParamMap.insert(std::make_pair("ttH", std::cref(tth)));
    ModelParamMap.insert(std::make_pair("BrHggRatio", std::cref(brhggratio)));
    ModelParamMap.insert(std::make_pair("BrHWWRatio", std::cref(brhwwratio)));
    ModelParamMap.insert(std::make_pair("BrHZZRatio", std::cref(brhzzratio)));
    ModelParamMap.insert(std::make_pair("BrHZgaRatio", std::cref(brhzgaratio)));
    ModelParamMap.insert(std::make_pair("BrHgagaRatio", std::cref(brhgagaratio)));
    ModelParamMap.insert(std::make_pair("BrHtautauRatio", std::cref(brhtautauratio)));
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
                return false;
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

double SigmaBR::BrHggRatio() const
{
    
    return brhggratio;

}

double SigmaBR::BrHWWRatio() const
{
    return brhwwratio;
}

double SigmaBR::BrHZZRatio() const
{
    return brhzzratio;
}

double SigmaBR::BrHZgaRatio() const
{
    return brhzgaratio;
}

double SigmaBR::BrHgagaRatio() const
{
    return brhgagaratio;
}

double SigmaBR::BrHtautauRatio() const
{
    return brhtautauratio;
}

double SigmaBR::BrHccRatio() const
{
    return brhccratio;
}

double SigmaBR::BrHbbRatio() const
{
    return brhbbratio;
}


