/* 
 * Copyright (C) 2017 HEPfit Collaboration
 *
 *
 * For the licensing terms see doc/COPYING.
 */

#include "THDMWSignalStrengths.h"

double THDMWSignalStrengths::muggH(const double sqrt_s) const
{
//    double lambda1=myTHDMW.getTHDMW_lambda1();
//    double lambda2=myTHDMW.getTHDMW_lambda2();
//    double lambda3=myTHDMW.getTHDMW_lambda3();
//
//    if (lambda1>0 && lambda2>0)
//    {
//        return lambda3+sqrt(lambda1*lambda2);
//    }
//    else
//    {
//        return -1.;
//    }
    return 1.;
}



//
//
//double HiggsChiral::muggH(const double sqrt_s) const
//{
////    return Gammagg() / (trueSM.computeBrHtogg()*trueSM.computeGammaHTotal());
//      
//    return computecg() * computecg();
//}
//
//double HiggsChiral::muVBF(const double sqrt_s) const
//{
//    return computecV() * computecV();
//}
//
//double HiggsChiral::mueeWBF(const double sqrt_s) const
//{
//    return computecV() * computecV();
//}
//
//double HiggsChiral::muWH(const double sqrt_s) const
//{
//    return computecV() * computecV();
//}
//
//double HiggsChiral::muZH(const double sqrt_s) const
//{
//    return computecV() * computecV();
//}
//
//double HiggsChiral::mueeZH(const double sqrt_s) const
//{
//    return computecV() * computecV();
//}
//
//double HiggsChiral::muVH(const double sqrt_s) const
//{
//    return computecV() * computecV();
//}
//
//double HiggsChiral::muVBFpVH(const double sqrt_s) const
//{
//    return computecV() * computecV();
//}
//
//double HiggsChiral::muttH(const double sqrt_s) const
//{
//    return computect() * computect();
//}
//
//double HiggsChiral::muggHpttH(const double sqrt_s) const
//{
//    double sigmaggH_SM = trueSM.computeSigmaggH(sqrt_s);
//    double sigmattH_SM = trueSM.computeSigmattH(sqrt_s);
//
//    double sigmaggH = muggH(sqrt_s) * sigmaggH_SM;
//    double sigmattH = muttH(sqrt_s) * sigmattH_SM;
//
//    return (sigmaggH + sigmattH) / (sigmaggH_SM + sigmattH_SM);
//} 
