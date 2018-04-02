/* 
 * Copyright (C) 2016 HEPfit Collaboration
 * All rights reserved.
 *
 * For the licensing terms see doc/COPYING.
 */

#include "GeneralTHDMEffectivePot.h"
#include <Math/Minimizer.h>
#include <Math/Factory.h>
#include <Math/Functor.h>

EffectivePot_GTHDM::EffectivePot_GTHDM(const StandardModel& SM_i)
//: myGTHDM(static_cast<const GeneralTHDM&> (SM_i)),
: mat_neutral(4,4,0.)
{}

const double* EffectivePot_GTHDM::potentialminimizer(double S1_start, double S2_start, double S3_start)
{
    ROOT::Math::Minimizer *minimizer = ROOT::Math::Factory::CreateMinimizer("Minuit","Migrad");
    
    minimizer->SetMaxFunctionCalls(10000);
    minimizer->SetMaxIterations(1000);
    minimizer->SetTolerance(0.01);
    
    ROOT::Math::Functor fwrap(this,&EffectivePot_GTHDM::potentialfunction,3);
    
    minimizer->SetFunction(fwrap);
    minimizer->SetVariable(0, "S1", S1_start, 0.1);
    minimizer->SetVariable(1, "S2", S2_start, 0.1);
    minimizer->SetVariable(2, "S3", S3_start, 0.1);
    
    minimizer->Minimize();
    
    const double *minimum = minimizer->X();
    
//    std::cout << "min_S1 = " << minimum[0] << std::endl;
//    std::cout << "min_S2 = " << minimum[1] << std::endl;
//    std::cout << "min_S3 = " << minimum[2] << std::endl << std::endl;
    
    return(minimum);
}

double EffectivePot_GTHDM::potentialfunction(const double *Svec)
{
    const double S1 = Svec[0];
    const double S2 = Svec[1];
    const double S3 = Svec[2];
    
    double Func = 3./2.*S1*S1*S1*S1 + 3.*S1*S1*S1 - 5.*S1 + S2*S2;
    
    double TraceMn2 = (Mneutral_2(S1,S2,S3))(0,0) + (Mneutral_2(S1,S2,S3))(1,1) +
                      (Mneutral_2(S1,S2,S3))(2,2) + (Mneutral_2(S1,S2,S3))(3,3);
    
    return Func + TraceMn2;
}

gslpp::matrix<double> EffectivePot_GTHDM::Mneutral_2(const double S1, const double S2, const double S3)
{
   
    mat_neutral.assign(0,0, S1);
    mat_neutral.assign(1,1, 2.*S2);
    mat_neutral.assign(2,2, 3.*S1);
    mat_neutral.assign(3,3, 4.*S3*S3);
    
    return mat_neutral;
}

EffectivePotMin1_GTHDM::EffectivePotMin1_GTHDM(const StandardModel& SM_i)
: ThObservable(SM_i), myEffectivePot_GTHDM(SM_i)
{}

double EffectivePotMin1_GTHDM::computeThValue()
{
    return (myEffectivePot_GTHDM.potentialminimizer(-0.39,0.,0.))[0];
}

EffectivePotMin2_GTHDM::EffectivePotMin2_GTHDM(const StandardModel& SM_i)
: ThObservable(SM_i), myEffectivePot_GTHDM(SM_i)
{}

double EffectivePotMin2_GTHDM::computeThValue()
{
    return (myEffectivePot_GTHDM.potentialminimizer(-0.38,0.,0.))[0];
}
