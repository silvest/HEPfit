/* 
 * Copyright (C) 2017 HEPfit Collaboration
 *
 *
 * For the licensing terms see doc/COPYING.
 */

#include "GeneralTHDMMatching.h"
#include "GeneralTHDM.h"
#include <gsl/gsl_integration.h>
#include <gsl/gsl_sf_dilog.h>
#include "GeneralTHDMcache.h"
#include <stdexcept>

GeneralTHDMMatching::GeneralTHDMMatching(const GeneralTHDM & GeneralTHDM_i) :

    StandardModelMatching(GeneralTHDM_i),
    myGTHDM(GeneralTHDM_i),
    myCKM(3, 3, 0.),
    mcdbs2(5, NDR, NLO),
    mcbtaunu(4, NDR, LO),
    mccleptonnu(4, NDR, LO),
    mcsleptonnu(4, NDR, LO),
    mculeptonnu(4, NDR, LO),
    mcBMll(13, NDR, NLO),
    mcbsg(8, NDR, NNLO),
    mcgminus2mu(2, NDR, NLO),
    mcbsmm(8, NDR, NNLO, NLO_QED22)
    
{
}
void GeneralTHDMMatching::updateGTHDMParameters()
{
    GF=myGTHDM.getGF();
    mMU=myGTHDM.getLeptons(StandardModel::MU).getMass();
}



//Maybe it's worth to define this at the cache, then we'd need to define one for each scalar
double GeneralTHDMMatching::F1oneloopgm2(const double ratio_sq){
    
    
    //gslpp::complex rsqc=ratio_sq;//We need the functions to be evaluated in the complex space so we define the ratio to be complex
    //gslpp::complex exact_result;
    
    double approx_result;
    
    //The exact expression has numerical problems for very small values of the ratio, better not to use it
    //if(ratio_sq==0.25){
    //    exact_result=-10 + 4*log(16);
    //}
    //else{
    //    exact_result=(-2*sqrt(1 - 4*rsqc)*rsqc + 3*sqrt(1 - 4*rsqc)*pow(rsqc,2) + 
    // 4*rsqc*(-1 + 2*rsqc)*arctanh(1/sqrt(1 - 4*rsqc)) + 
    // 4*rsqc*(-1 + 2*rsqc)*arctanh((-1 + 2*rsqc)/sqrt(1 - 4*rsqc)) + 
    // 3*rsqc*log(1 - sqrt(1 - 4*rsqc)) - 3*rsqc*log(1 + sqrt(1 - 4*rsqc)) - 
    // log((1 + sqrt(1 - 4*rsqc))/(1 - sqrt(1 - 4*rsqc))) - 
    // 3*rsqc*log(1 - sqrt(1 - 4*rsqc) - 2*rsqc) + 
    // 3*rsqc*log(1 + sqrt(1 - 4*rsqc) - 2*rsqc) - sqrt(1 - 4*rsqc)*log(rsqc) + 
    // 3*sqrt(1 - 4*rsqc)*rsqc*log(rsqc))/(2.*sqrt(1 - 4*rsqc)*pow(rsqc,3));
    //}
    
    approx_result=(-7.0/6.0-log(ratio_sq));
    
    return approx_result;
}



double GeneralTHDMMatching::F2oneloopgm2(const double ratio_sq){
    
    
    //gslpp::complex rsqc=ratio_sq;//We need the functions to be evaluated in the complex space so we define the ratio to be complex
    //gslpp::complex exact_result;
    
    double approx_result;
    
    //The exact expression has numerical problems for very small values of the ratio, better not to use it
    //if(ratio_sq==0.25){
    //    exact_result=-34 + 24*log(4);
    //}
    //else{
    //    
    //    exact_result=((-2 + 6*rsqc)*arccot(sqrt(-1 + 4*rsqc)) + 
    //    (-2 + 6*rsqc)*arctan((-1 + 2*rsqc)/sqrt(-1 + 4*rsqc)) + 
    //    sqrt(-1 + 4*rsqc)*(-(rsqc*(2 + rsqc)) + (-1 + rsqc)*log(rsqc)))/(2.*pow(rsqc,3)*sqrt(-1 + 4*rsqc));
    //}
    
    approx_result=(11.0/6.0+log(ratio_sq));
    
    return approx_result;
}



double GeneralTHDMMatching::F3oneloopgm2(const double ratio_sq){
    
    
    //gslpp::complex exact_result;
    double approx_result;
    
    //gslpp::complex rsqc=ratio_sq;//We need the functions to be evaluated in the complex space so we define the ratio to be complex

    //For our range of masses the expansion is really good so we keep only it
    //exact_result = ((-2 + rsqc)*rsqc + 2*(-1 + rsqc)*log(1 - rsqc))/(2.*pow(rsqc,3));
    
    approx_result = -1/6. - ratio_sq/12.;
    
    //std::cout<<"\033[1;31m exact_result = \033[0m "<< exact_result <<std::endl;
    //std::cout<<"\033[1;31m approx_result = \033[0m "<< approx_result <<"\n"<<std::endl;
    
    
    return approx_result;
}



double GeneralTHDMMatching::gminus2muLO() {
    
    //We'll include the contributions from the paper [1502.04199] of V. Ilisie 
    
    updateGTHDMParameters();
    double gminus2muLO;
    
     //add something to note that this is only valid in the aligned case and in the CP-conserving limit
    
    
    
    double pi=M_PI;
    double GF=myGTHDM.getGF();
    double mMU=myGTHDM.getLeptons(StandardModel::MU).getMass();
    
    
 
    
    /*Mass of the physical scalars*/
    double mH1_2 = myGTHDM.getMyGTHDMCache()->mH1sq;
    double mH2_2 = myGTHDM.getMyGTHDMCache()->mH2sq;
    double mH3_2 = myGTHDM.getMyGTHDMCache()->mH3sq;
    double mHp2 = myGTHDM.getMyGTHDMCache()->mHp2;
    
    //gslpp::complex yl1 = myGTHDM.getyl1();
    gslpp::complex yl1 = myGTHDM.getMyGTHDMCache()->yl1;
    gslpp::complex yl2 = myGTHDM.getMyGTHDMCache()->yl2;
    gslpp::complex yl3 = myGTHDM.getMyGTHDMCache()->yl3;
    gslpp::complex sl = myGTHDM.getMyGTHDMCache()->sl;
    
    //double alpha1 = myGTHDM.getalpha1();
    ///////////////////////
    //OLD part of the code, this doesn't make sense anymore
    //double eta = 0.0;
    //if(myGTHDM.getSMHiggs()){
    //    eta = alpha1;
    //}
    //else{
    //    eta = alpha1 - pi/2.0;
    //}
        
    //This was obviously WRONG, here you're expanding for small values of alpha but in the 
    //non "heavy" case this doesn't make any sense! Let's forget about the expansion and consider
    //the most general case. We will take the Yukawas from the cache
    //double yl_h = 1.0 + eta*sl.real();
    //double yl_H =  -sl.real() + eta;
    //double yl_A = - sl.real();  
    ////////////////////////
    

    
    double rmu_hSM, rmu_h, rmu_H, rmu_A, rmu_Hp;
    double part_hSM, part_h, part_H, part_A, part_Hp;
    
    if( mH1_2<100.0 || mH2_2<100.0 || mH3_2<100.0 || mHp2<100.0)
    {
        throw std::runtime_error("The implemented approximation for g-2_\\mu only works for Higgs masses above 10 GeV.");
    }
    
    if (!myGTHDM.getATHDMflag())
    {
        throw std::runtime_error("(g-2) is only available in the A2HDM.");
    }
    
     if (!myGTHDM.getCPconservationflag())
    {
         //Now the CP-violating case is being included, check again if we can remove this
        throw std::runtime_error("(g-2) is only available in the CP-conserving limit.");
    }
     
    
    rmu_hSM=mMU*mMU/mH1_2;
    rmu_h=mMU*mMU/mH1_2;
    rmu_H=mMU*mMU/mH2_2;
    rmu_A=mMU*mMU/mH3_2;
    rmu_Hp=mMU*mMU/mHp2;
    
    part_hSM=rmu_hSM*(F1oneloopgm2(rmu_hSM));
    part_h=rmu_h*(yl1.real()*yl1.real()*F1oneloopgm2(rmu_h)+yl1.imag()*yl1.imag()*F2oneloopgm2(rmu_h));
    part_H=rmu_H*(yl2.real()*yl2.real()*F1oneloopgm2(rmu_H)+yl2.imag()*yl2.imag()*F2oneloopgm2(rmu_H));
    part_A=rmu_A*(yl3.real()*yl3.real()*F1oneloopgm2(rmu_A)+yl3.imag()*yl3.imag()*F2oneloopgm2(rmu_A));;
    part_Hp=rmu_Hp*sl.abs2()*F3oneloopgm2(rmu_Hp);
    
    
   // We're interested only in the NP contribution so we need to remove the SM part
   gminus2muLO=GF*mMU*mMU/(4.0*pi*pi*sqrt(2.0)) * (-part_hSM+part_h+part_H+part_A+part_Hp);
   
       return(gminus2muLO);
}

struct __f_params{
  double a;
};

double __fS2_integrand(double x, void* p){
  __f_params &params= *reinterpret_cast<__f_params *>(p);
  double integ = (2.0*x*(1.0-x)-1.0) * log(x*(1.0-x)/params.a) / (x*(1.0-x)-params.a);
  return integ;
}

double __fPS2_integrand(double x, void* p){
  __f_params &params= *reinterpret_cast<__f_params *>(p);
  double integ = log(x*(1.0-x)/params.a) / (x*(1.0-x)-params.a);
  return integ;
}


double gscalar2(double r)
{
    __f_params params;
    params.a=r;

    double result;
    gsl_integration_glfixed_table * w
        = gsl_integration_glfixed_table_alloc(100);
    gsl_function F;

    F.function = &__fS2_integrand;
    F.params = reinterpret_cast<void *>(&params);

    result = gsl_integration_glfixed (&F, 0, 1, w);

    gsl_integration_glfixed_table_free (w);
    
  return result;
}

double gpseudoscalar2(double r)
{
    __f_params params;
    params.a=r;

    double result;
    gsl_integration_glfixed_table * w
        = gsl_integration_glfixed_table_alloc(100);
    gsl_function F;

    F.function = &__fPS2_integrand;
    F.params = reinterpret_cast<void *>(&params);

    result = gsl_integration_glfixed (&F, 0, 1, w);

    gsl_integration_glfixed_table_free (w);

  return result;
}



gslpp::complex GeneralTHDMMatching::negsquareroot(double x){
    gslpp::complex result;
    if (x > 0)
    {
        result = sqrt(x);
    }
    else{
        result = gslpp::complex::i()*sqrt(-x);
    }
    return result;
}

    gslpp::complex GeneralTHDMMatching::negpow(double basis, double exp){
    gslpp::complex result;
    if (basis > 0)
    {
        result = pow(basis, exp);
    }
    else{
        result = pow(-basis, exp)*(cos(M_PI*exp) + gslpp::complex::i()*sin(M_PI*exp));
    }
    return result;
}
    
        gslpp::complex GeneralTHDMMatching::neglog(gslpp::complex argument){
            
    gslpp::complex result;
    
      result = log(argument.abs()) + gslpp::complex::i()*argument.arg();
   
    
    return result;
}

    /*
    gslpp::complex GeneralTHDMMatching::TF(double m1, double m2, double m3){
            
        double pi=M_PI;
        
        double ml, mm, mh;
        ml=0.0; mm=0.0; mh=0.0;
        if(m1<=m3 && m3<=m2)
    {
        //1<3<2 swap 2 and 3
        ml  = m1;
        mm = m3;
        mh  = m2;
    }
    else if(m3<=m2 && m2<=m1)
    {
        //3<2<1 swap 1 and 3
        ml  = m3;
        mm = m2;
        mh  = m1;
    }
    else if(m2<=m1 && m1<=m3)
    {
        //2<1<3 swap 1 and 2
        ml  = m2;
        mm = m1;
        mh  = m3;
    }
    else if(m2<=m3 && m3<=m1)
    {
        //2<3<1: 3->2, 1->3, 2->1
        ml  = m2;
        mm = m3;
        mh  = m1;
    }
    else if(m3<=m1 && m1<=m2)
    {
        //3<1<2 3->1, 1->2, 2->3
        ml  = m3;
        mm = m1;
        mh  = m2;
    }
    
    else if(m1<=m2 && m2<=m3)
    {
        //1<2<3 ok
        ml  = m1;
        mm = m2;
        mh  = m3;
    }
        
        
        gslpp::complex lambda = negsquareroot(ml*ml*ml*ml + mm*mm*mm*mm + mh*mh*mh*mh 
         - 2.0*ml*ml*mm*mm - 2.0*mm*mm*mh*mh - 2.0*mh*mh*ml*ml);
        gslpp::complex ap = 1.0/(2.0*mh*mh)*( mh*mh + ml*ml - mm*mm - lambda);
        gslpp::complex am = 1.0/(2.0*mh*mh)*( mh*mh - ml*ml + mm*mm - lambda);
            
         gslpp::complex result = (lambda/2.0)*(2.0*neglog(ap)*neglog(am) 
         - neglog(ml*ml/(mh*mh))*neglog(mm*mm/(mh*mh)) - 2.0*PolyLog.Li2(ap) 
         - 2.0*PolyLog.Li2(am) + pi*pi/3.0);
         
        return result;
        }
    */


        
double GeneralTHDMMatching::F1twoloopgm2(const double ratio_sq){
    
    gslpp::complex rsqc=ratio_sq;//We need the functions to be evaluated in the complex space so we define the ratio to be complex
    double final_result;
    
    if(ratio_sq==0.25){
        final_result=-2*ratio_sq;
    }
    else{
        //Let's expand when the function is small to avoid numerical instabilities with the complete function
        if(ratio_sq<0.001){
        final_result=-(1/6.)*(ratio_sq*(12 + M_PI*M_PI - 12*ratio_sq - 9*ratio_sq*ratio_sq + 
        2*M_PI*M_PI*ratio_sq*ratio_sq + 6*(1 + 2*ratio_sq + 3*ratio_sq*ratio_sq)*log(ratio_sq) + 
        (3 + 6*ratio_sq*ratio_sq)*log(ratio_sq)*log(ratio_sq)));
        }
        else{
           final_result = ((rsqc*((pow(M_PI,2) - 12*sqrt(1 - 4*rsqc) - 6*sqrt(1 - 4*rsqc)*log(rsqc) + 
            3*pow(log((-2*rsqc)/(-1 + sqrt(1 - 4*rsqc) + 2*rsqc)),2) - 
            2*rsqc*(pow(M_PI,2) + 6*pow(log(1 - sqrt(1 - 4*rsqc)),2) - 
            6*pow(log(1 + sqrt(1 - 4*rsqc)),2) + 
            2*arctanh(sqrt(1 - 4*rsqc))*log(4096*pow(rsqc,6)) + 
            3*pow(log((-2*rsqc)/(-1 + sqrt(1 - 4*rsqc) + 2*rsqc)),2)))/6. + 
            (2 - 4*rsqc)*PolyLog.Li2(-0.5*(1 + sqrt(1 - 4*rsqc) - 2*rsqc)/rsqc)))/sqrt(1 - 4*rsqc)).real();
        }
        
    }
    
    //std::cout<<"\033[1;31m final_result = \033[0m "<< final_result <<std::endl;
    return (final_result);
}
      



double GeneralTHDMMatching::F1tildetwoloopgm2(const double ratio_sq){
    
    gslpp::complex rsqc=ratio_sq;//We need the functions to be evaluated in the complex space so we define the ratio to be complex
    double final_result;
    
    if(ratio_sq==0.25){
        final_result=log(2);
    }
    else{
        //Let's expand when the function is small to avoid numerical instabilities with the complete function
        if(ratio_sq<0.001){
        final_result=(ratio_sq*(pow(M_PI,2) + 3*pow(log(ratio_sq),2)))/6. + (pow(ratio_sq,2)*(-6 + pow(M_PI,2) 
                + 6*log(ratio_sq) + 3*pow(log(ratio_sq),2)))/3. + (pow(ratio_sq,3)*(-11 + 2*pow(M_PI,2) 
                + 14*log(ratio_sq) + 6*pow(log(ratio_sq),2)))/2.;
        }
        else{
           final_result = ((rsqc*(-PolyLog.Li2(-0.5*(1 + sqrt(1 - 4*rsqc) - 2*rsqc)/rsqc) + 
                   PolyLog.Li2((-1 + sqrt(1 - 4*rsqc) + 2*rsqc)/(2.*rsqc))))/sqrt(1 - 4*rsqc)).real();
        }
        
    }
    
    //std::cout<<"\033[1;31m final_result = \033[0m "<< final_result <<std::endl;
    return (final_result);
}





double GeneralTHDMMatching::F2twoloopgm2(const double ratio_sq){
    
    gslpp::complex rsqc=ratio_sq;//We need the functions to be evaluated in the complex space so we define the ratio to be complex
    double final_result;
    
    if(ratio_sq==0.25){
        final_result=(1-log(4))/4;
    }
    else{
        //Let's expand when the function is small to avoid numerical instabilities with the complete function
        if(ratio_sq<0.001){
        final_result=(ratio_sq*(2 + log(ratio_sq)))/2. + (pow(ratio_sq,2)*(-pow(M_PI,2) - 3*pow(log(ratio_sq),2)))/6. 
                + (pow(ratio_sq,3)*(6 - pow(M_PI,2) - 6*log(ratio_sq) - 3*pow(log(ratio_sq),2)))/3.;
        }
        else{
           final_result = ((rsqc*(6 + 3*log(rsqc) + (rsqc*(pow(M_PI,2) + 6*arctanh(sqrt(1 - 4*rsqc))*log((-2*rsqc)/(-1 + sqrt(1 - 4*rsqc) 
                   + 2*rsqc)) + 12*PolyLog.Li2(-0.5*(1 + sqrt(1 - 4*rsqc) - 2*rsqc)/rsqc)))/sqrt(1 - 4*rsqc)))/6.).real();
        }
        
    }
    
    //std::cout<<"\033[1;31m final_result = \033[0m "<< final_result <<std::endl;
    return (final_result);
}





double GeneralTHDMMatching::F3twoloopgm2(const double ratio_sq){
    
    gslpp::complex rsqc=ratio_sq;//We need the functions to be evaluated in the complex space so we define the ratio to be complex
    double final_result;
    
    if(ratio_sq==0.25){
        final_result=19./16;
    }
    else{
        //Let's expand when the function is small to avoid numerical instabilities with the complete function
        if(ratio_sq<0.001){
        final_result= (ratio_sq*(2 + log(ratio_sq)))/2. + (pow(ratio_sq,3)*(-51 + pow(M_PI,2) + 51*log(ratio_sq) + 3*pow(log(ratio_sq),2)))/3. 
                + (pow(ratio_sq,2)*(180 + 17*pow(M_PI,2) + 90*log(ratio_sq) + 51*pow(log(ratio_sq),2)))/12.;
        }
        else{
           final_result = ((rsqc*(12*sqrt(1 - 4*rsqc)*(1 + 15*rsqc) + pow(M_PI,2)*rsqc*(-17 + 30*rsqc) + 6*sqrt(1 - 4*rsqc)*log(rsqc) 
                   + 3*rsqc*((17 - 30*rsqc)*log(1 + sqrt(1 - 4*rsqc))*log((1 + sqrt(1 - 4*rsqc))/pow(rsqc,2)) + 30*sqrt(1 - 4*rsqc)*log(rsqc) 
                   + (-17 + 30*rsqc)*(8*arctanh(sqrt(1 - 4*rsqc))*log(2) + log(1 - sqrt(1 - 4*rsqc))*(3*log(1 - sqrt(1 - 4*rsqc)) 
                   - 2*log((1 + sqrt(1 - 4*rsqc))*rsqc)))) + 12*rsqc*(-17 + 30*rsqc)*PolyLog.Li2(-0.5*(1 + sqrt(1 - 4*rsqc) - 2*rsqc)/rsqc)))/(12.*sqrt(1 - 4*rsqc))).real();
        }
        
    }
    
    //std::cout<<"\033[1;31m final_result = \033[0m "<< final_result <<std::endl;
    return (final_result);
}





double GeneralTHDMMatching::F4twoloopgm2(const double ratio_sq){
    
    gslpp::complex rsqc=ratio_sq;//We need the functions to be evaluated in the complex space so we define the ratio to be complex
    double final_result;
    
    if(ratio_sq==1.){
        final_result=(-51 + 2*pow(M_PI,2))/36.;
    }
    else{
    //Let's expand when the function is small to avoid numerical instabilities with the complete function
        if(ratio_sq<0.001){
            final_result= (pow(ratio_sq,3)*(13 - 3*log(ratio_sq)))/27. + 
            (-87 - 16*pow(M_PI,2) - 78*log(ratio_sq) - 24*pow(log(ratio_sq),2))/36. + 
            (pow(ratio_sq,2)*(-14 - 2*pow(M_PI,2) + 12*log(ratio_sq) - 3*pow(log(ratio_sq),2)))/
            6. + (ratio_sq*(21 + 8*pow(M_PI,2) - 3*log(ratio_sq) + 12*pow(log(ratio_sq),2)))/9.;
        }
        else{
            final_result = (-29./12 + rsqc + ((-13 + 6*rsqc + 
            2*(4 - 8*rsqc + 3*pow(rsqc,2))*log((-1 + rsqc)/rsqc))*log(rsqc))/6. - 
            ((4 - 8*rsqc + 3*pow(rsqc,2))*PolyLog.Li2(1/rsqc))/3.).real();
        }
    }    
    
    //std::cout<<"\033[1;31m final_result = \033[0m "<< final_result <<std::endl;
    return (final_result);
}

       




double GeneralTHDMMatching::F5twoloopgm2(const double ratio_sq){
    
    gslpp::complex rsqc=ratio_sq;//We need the functions to be evaluated in the complex space so we define the ratio to be complex
    double final_result;
    
    if(ratio_sq==1.){
        final_result=(-21 + 2*pow(M_PI,2))/36.;
    }
    else{
    //Let's expand when the function is small to avoid numerical instabilities with the complete function
        if(ratio_sq<0.001){
            final_result= (5 + 2*log(ratio_sq))/12. + (pow(ratio_sq,3)*(-5 + 4*log(ratio_sq)))/6. + 
            (ratio_sq*(-9 - 2*pow(M_PI,2) - 9*log(ratio_sq) - 3*pow(log(ratio_sq),2)))/9. + 
            (pow(ratio_sq,2)*(4 + 2*pow(M_PI,2) - 4*log(ratio_sq) + 3*pow(log(ratio_sq),2)))/6.;
        }
        else{
            final_result = ((5 - 12*rsqc + 2*(1 - 6*rsqc + 2*(2 - 3*rsqc)*rsqc*log((-1 + rsqc)/rsqc))*
            log(rsqc) + 4*rsqc*(-2 + 3*rsqc)*PolyLog.Li2(1/rsqc))/12.).real();
        }
    }    
    
    //std::cout<<"\033[1;31m final_result = \033[0m "<< final_result <<std::endl;
    return (final_result);
}




double GeneralTHDMMatching::gminus2muNLO() {
    
    updateGTHDMParameters();
    
    if (!myGTHDM.getSMHiggs())
    {
        //throw std::runtime_error("The NLO computation of g-2 is only implemented for SM Higgs since the coupling of charged scalars to neutral scalars is included only for that case.");
        throw std::runtime_error("The NLO computation of g-2 for the heavy higgs scenario must be reviewed although in principle is implemented");
    
    }
    
    
    
    //SM constants
    double Ale = myGTHDM.getAle();
    double vev = myGTHDM.v();
    double sW2 = myGTHDM.sW2();
    double Vtb2=(myGTHDM.getVCKM()(2,2)).abs2();
    double Nc=3.;
    double Qu=2/3;
    double Qd=-1/3;
    double Ql=-1;
    
    //SM masses
    double mMU=myGTHDM.getLeptons(StandardModel::MU).getMass();
    double mTAU=myGTHDM.getLeptons(StandardModel::TAU).getMass();
    double mt=myGTHDM.getQuarks(QCD::TOP).getMass();
    double mb=myGTHDM.getQuarks(QCD::BOTTOM).getMass();
    double Mw=myGTHDM.Mw();
    
    //std::cout<<"\033[1;31m Vtb = \033[0m "<< Vtb <<std::endl;

    //GTHDM Yukawa couplings
    gslpp::complex yu1 = myGTHDM.getMyGTHDMCache()->yu1;
    double yu1R=yu1.real();
    double yu1I=yu1.imag();
    gslpp::complex yu2 = myGTHDM.getMyGTHDMCache()->yu2;
    double yu2R=yu2.real();
    double yu2I=yu2.imag();
    gslpp::complex yu3 = myGTHDM.getMyGTHDMCache()->yu3;
    double yu3R=yu3.real();
    double yu3I=yu3.imag();
    
    gslpp::complex yd1 = myGTHDM.getMyGTHDMCache()->yd1;
    double yd1R=yd1.real();
    double yd1I=yd1.imag();
    gslpp::complex yd2 = myGTHDM.getMyGTHDMCache()->yd2;
    double yd2R=yd2.real();
    double yd2I=yd2.imag();
    gslpp::complex yd3 = myGTHDM.getMyGTHDMCache()->yd3;
    double yd3R=yd3.real();
    double yd3I=yd3.imag();
    
    gslpp::complex yl1 = myGTHDM.getMyGTHDMCache()->yl1;
    double yl1R=yl1.real();
    double yl1I=yl1.imag();
    gslpp::complex yl2 = myGTHDM.getMyGTHDMCache()->yl2;
    double yl2R=yl2.real();
    double yl2I=yl2.imag();
    gslpp::complex yl3 = myGTHDM.getMyGTHDMCache()->yl3;
    double yl3R=yl3.real();
    double yl3I=yl3.imag();
    
    gslpp::complex su = myGTHDM.getMyGTHDMCache()->su;
    gslpp::complex sd = myGTHDM.getMyGTHDMCache()->sd;
    gslpp::complex sl = myGTHDM.getMyGTHDMCache()->sl;
    
    //GTHDM masses
    double mH1_2 = myGTHDM.getMyGTHDMCache()->mH1sq;
    double mH2_2 = myGTHDM.getMyGTHDMCache()->mH2sq;
    double mH3_2 = myGTHDM.getMyGTHDMCache()->mH3sq;
    double mHp2 = myGTHDM.getMyGTHDMCache()->mHp2;
    
    //GTHDM Rotation Couplings
    double R11 = myGTHDM.getMyGTHDMCache()->R11;
    double R12 = myGTHDM.getMyGTHDMCache()->R12;
    double R13 = myGTHDM.getMyGTHDMCache()->R13;
    double R21 = myGTHDM.getMyGTHDMCache()->R21;
    double R22 = myGTHDM.getMyGTHDMCache()->R22;
    double R23 = myGTHDM.getMyGTHDMCache()->R23;
    double R31 = myGTHDM.getMyGTHDMCache()->R31;
    double R32 = myGTHDM.getMyGTHDMCache()->R32;
    double R33 = myGTHDM.getMyGTHDMCache()->R33;
    
    double Rh1;
    double RH1;
    double RA1;
    
    
    if (!myGTHDM.getSMHiggs())
    {
        RH1=R11;
        Rh1=R21;
        RA1=R31;
    }
    else{
        Rh1=R11;
        RH1=R21;
        RA1=R31;
    }
    
    //GTHDM lambda couplings
    double lambda3 = myGTHDM.getlambda3();
    double Relambda7 = myGTHDM.getRelambda7();
    double Imlambda7 = myGTHDM.getImlambda7();
    
    //GTHDM coupling of charged scalar with neutral scalars
    double lambdaHph;
    double lambdaHpH;
    double lambdaHpA;
    
    lambdaHpA= R31*lambda3+R32*Relambda7-R33*Imlambda7;
    
    if (!myGTHDM.getSMHiggs())
    {
        lambdaHpH=R11*lambda3+R12*Relambda7-R13*Imlambda7;
        lambdaHph=R21*lambda3+R22*Relambda7-R23*Imlambda7;
    }
    else{
        lambdaHph=R11*lambda3+R12*Relambda7-R13*Imlambda7;
        lambdaHpH=R21*lambda3+R22*Relambda7-R23*Imlambda7;
    }
    
    
    
    //top-quark ratios square
    double rsqt_h=mt*mt/mH1_2;
    double rsqt_H=mt*mt/mH2_2;
    double rsqt_A=mt*mt/mH3_2;
    double rsqt_Hp=mt*mt/mHp2;
    double rsqt_W=mt*mt/(Mw*Mw);
    //bottom-quark ratios square
    double rsqb_h=mb*mb/mH1_2;
    double rsqb_H=mb*mb/mH2_2;
    double rsqb_A=mb*mb/mH3_2;
    double rsqb_Hp=mb*mb/mHp2;
    //tau-lepton ratios square
    double rsqtau_h=mTAU*mTAU/mH1_2;
    double rsqtau_H=mTAU*mTAU/mH2_2;
    double rsqtau_A=mTAU*mTAU/mH3_2;
    double rsqtau_Hp=mTAU*mTAU/mHp2;
    //charged-higgs ratios square
    double rsqHp_h=mHp2/mH1_2;
    double rsqHp_H=mHp2/mH2_2;
    double rsqHp_A=mHp2/mH3_2;
    //W-boson ratios square
    double rsqW_h=Mw*Mw/mH1_2;
    double rsqW_H=Mw*Mw/mH2_2;
    double rsqW_A=Mw*Mw/mH3_2;
    double rsqW_Hp=Mw*Mw/mHp2;
    
    
    
    //Computation a1
    double a1_const=(Ale*mMU*mMU/(4*M_PI*M_PI*M_PI*vev*vev));
    double a1_SM;
    double a1_h;
    double a1_H;
    double a1_A;    
    double a1_total;
    
    
    a1_SM=a1_const*(Nc*Qu*Qu*F1twoloopgm2(rsqt_h)+Nc*Qd*Qd*F1twoloopgm2(rsqb_h)+F1twoloopgm2(rsqtau_h));
    
    a1_h=a1_const*(Nc*Qu*Qu*(yu1R*yu1R*F1twoloopgm2(rsqt_h)+yu1I*yu1I*F1tildetwoloopgm2(rsqt_h))
            +Nc*Qd*Qd*(yd1R*yd1R*F1twoloopgm2(rsqb_h)+yd1I*yd1I*F1tildetwoloopgm2(rsqb_h))
            +(yl1R*yl1R*F1twoloopgm2(rsqtau_h)+yl1I*yl1I*F1tildetwoloopgm2(rsqtau_h)));

    a1_H=a1_const*(Nc*Qu*Qu*(yu2R*yu2R*F1twoloopgm2(rsqt_H)+yu2I*yu2I*F1tildetwoloopgm2(rsqt_H))
            +Nc*Qd*Qd*(yd2R*yd2R*F1twoloopgm2(rsqb_H)+yd2I*yd2I*F1tildetwoloopgm2(rsqb_H))
            +(yl2R*yl2R*F1twoloopgm2(rsqtau_H)+yl2I*yl2I*F1tildetwoloopgm2(rsqtau_H)));
    
    a1_A=a1_const*(Nc*Qu*Qu*(yu3R*yu3R*F1twoloopgm2(rsqt_A)+yu3I*yu3I*F1tildetwoloopgm2(rsqt_A))
            +Nc*Qd*Qd*(yd3R*yd3R*F1twoloopgm2(rsqb_A)+yd3I*yd3I*F1tildetwoloopgm2(rsqb_A))
            +(yl3R*yl3R*F1twoloopgm2(rsqtau_A)+yl3I*yl3I*F1tildetwoloopgm2(rsqtau_A)));
    
    a1_total=-a1_SM+a1_h+a1_H+a1_A;
    
    
    
    //Computation a2
    double a2_const=(Ale*mMU*mMU/(8*M_PI*M_PI*M_PI));
    double a2_SM=0.;//The SM contribution for this part is zero
    double a2_h;
    double a2_H;
    double a2_A;
    double a2_total;
    
    
    a2_h=a2_const*(yl1R/mH1_2)*lambdaHph*F2twoloopgm2(rsqHp_h);
    a2_H=a2_const*(yl1R/mH2_2)*lambdaHpH*F2twoloopgm2(rsqHp_H);
    a2_A=a2_const*(yl1R/mH3_2)*lambdaHpA*F2twoloopgm2(rsqHp_A);
    
    a2_total=-a2_SM+a2_h+a2_H+a2_A;
    
    
    
    //Computation a3 
    double a3_const=(Ale*mMU*mMU/(8*M_PI*M_PI*M_PI*vev*vev));
    double a3_SM;
    double a3_h;
    double a3_H;
    double a3_A;
    double a3_total;
            
    a3_SM=a3_const*F3twoloopgm2(rsqW_h);
    a3_h=a3_const*yl1R*Rh1*F3twoloopgm2(rsqW_h);
    a3_H=a3_const*yl2R*RH1*F3twoloopgm2(rsqW_H);
    a3_A=a3_const*yl3R*RA1*F3twoloopgm2(rsqW_A);
            
    a3_total=-a3_SM+a3_h+a3_H+a3_A;      
       
    //Computation a4
    double a4_const=(Ale*mMU*mMU*Nc*Vtb2/(32*M_PI*M_PI*M_PI*vev*vev*sW2*(mHp2-Mw*Mw)));
    double a4_Hp_t;
    double a4_Hp_b;
    double a4_total;
            
    a4_Hp_t=a4_const*((su*(sl.conjugate())).real())*mt*mt*(F4twoloopgm2(rsqt_Hp)-F4twoloopgm2(rsqt_W));
    a4_Hp_t=a4_const*((sd*(sl.conjugate())).real())*mb*mb*(F5twoloopgm2(rsqt_Hp)-F5twoloopgm2(rsqt_W));
            
    a4_total=a4_Hp_t+a4_Hp_b;  
    

   
    double gminus2muNLO = a1_total+a2_total+a3_total+a4_total;
  
    return(gminus2muNLO);
    }




std::vector<WilsonCoefficient>& GeneralTHDMMatching::CMgminus2mu() {
    
   /* gslpp::complex su = myGTHDM.getNu_11();
    gslpp::complex sd = myGTHDM.getNd_11();
    gslpp::complex sl = myGTHDM.getNl_11();*/
    
  /*  double Imlambda5=myGTHDM.getImlambda5();
    double Imlambda6=myGTHDM.getMyGTHDMCache()->Imlambda6;
    double Imlambda7=myGTHDM.getMyGTHDMCache()->Imlambda7;
    double sinalpha2=myGTHDM.getsinalpha2();
    double sinalpha3=myGTHDM.getsinalpha3();*/
    
    vmcgminus2mu = StandardModelMatching::CMgminus2mu();

    double gminus2muLOvalue=gminus2muLO();
    double gminus2muNLOvalue=gminus2muNLO();
    //gminus2muNLOvalue=0.; //For the moment we set the NLO contribution to zero

   /* std::cout << " gminus2muLOvalue =  " << gminus2muLOvalue << std::endl;
    std::cout << " gminus2muNLOvalue =  " << gminus2muNLOvalue << std::endl;*/
    

    
    switch (mcgminus2mu.getOrder()) {
        case LO:
            mcgminus2mu.setCoeff(0, gminus2muLOvalue, LO);  //g-2_muR
            mcgminus2mu.setCoeff(1, 0., LO);  //g-2_muL
            break;
        case NLO:
            mcgminus2mu.setCoeff(0, gminus2muLOvalue+gminus2muNLOvalue, NLO);  //g-2_muR
            mcgminus2mu.setCoeff(1, 0., NLO);  //g-2_muL
            break;
        case NNLO:
        default:
            std::stringstream out;
            out << mcgminus2mu.getOrder();
            throw std::runtime_error("GeneralTHDMMatching::CMgminus2mu(): order " + out.str() + " not implemented.\nOnly leading order (LO) or next-to-leading order (NLO) are allowed.");
    }

    vmcgminus2mu.push_back(mcgminus2mu);
    return(vmcgminus2mu);
   
}

/*******************************************************************************
 * Wilson coefficients meson mixing           *  
 * ****************************************************************************/


std::vector<WilsonCoefficient>& GeneralTHDMMatching::CMdbs2() {

    double Mut = myGTHDM.getMut();
    double xt = x_t(Mut);
    double GF=myGTHDM.getGF();
    double MW=myGTHDM.Mw();
    gslpp::complex co = GF / 4. / M_PI * MW * myGTHDM.getCKM().computelamt_s();
    //double tanb = myGTHDM.gettanb();
    double mHp2=myGTHDM.getmHp2();
    double xHW=mHp2/(MW*MW);
    double xtH=xt/xHW;
    double mb=myGTHDM.getQuarks(QCD::BOTTOM).getMass();
    double sd = (myGTHDM.getNd_11()).real();
    double su = (myGTHDM.getNu_11()).real();
    double SWH=xtH*((2.0*xHW-8.0)*log(xtH)/((1.0-xHW)*(1.0-xtH)*(1.0-xtH))+6.0*xHW*log(xt)/((1.0-xHW)*(1.0-xt)*(1.0-xt))-(8.0-2.0*xt)/((1.0-xt)*(1.0-xtH)))*su*su;//su*su = sigu.abs2()
    double SHH=xtH*((1.0+xtH)/((1.0-xtH)*(1.0-xtH))+2.0*xtH*log(xtH)/((1.0-xtH)*(1.0-xtH)*(1.0-xtH)))*su*su*su*su;//su*su*su*su = sigu.abs2()*sigu.abs2()
    double C1bsSRR = 4.0*mb*mb*xt*xtH*sd*su/(mHp2*(xtH-1.0)*(xtH-1.0)*(xtH-1.0)) //sd*su = sigd*sigu.conjugate()
                     * (sd*su*(2.0*(xtH-1.0)-(xtH+1.0)*log(xtH)) //sd*su = sigd*sigu.conjugate()
                        + (2.0*xt*xt*(xtH-1.0)*(xtH-1.0)*(xtH-1.0)*log(xt)/(xt-1.0)
                           +2.0*xt*(xtH-1.0)*((xt-xtH)*(xtH-1.0)+(xtH-xt*xtH)*log(xtH)))/((xt-1.0)*(xt-xtH)*xtH));

    vmcds = StandardModelMatching::CMdbs2();
    mcdbs2.setMu(Mut);

    switch (mcdbs2.getOrder()) {
        case NNLO:
        case NLO:
        case LO:
            mcdbs2.setCoeff(0, co * co * xt * (SWH+SHH), LO);
            break;
        default:
            std::stringstream out;
            out << mcdbs2.getOrder();
            throw std::runtime_error("GeneralTHDMMatching::CMdbs2(): order " + out.str() + "not implemented");
    }

    vmcds.push_back(mcdbs2);
    //The following are the primed coefficients.
    mcdbs2.setMu(Mut);

    switch (mcdbs2.getOrder()) {
        case NNLO:
        case NLO:
        case LO:
            mcdbs2.setCoeff(1, co * co * C1bsSRR, LO);
            break;
        default:
            std::stringstream out;
            out << mcdbs2.getOrder();
            throw std::runtime_error("GeneralTHDMMatching::CMdbs2(): order " + out.str() + "not implemented");
    }

    vmcds.push_back(mcdbs2);

    return(vmcds);
}
//
//std::vector<WilsonCoefficient>& GeneralTHDMMatching::CMdbsp2() {
//
//    double Mut = myGTHDM.getMut();
//    double xt = x_t(Mut);
//    double GF=myGTHDM.getGF();
//    double MW=myGTHDM.Mw();
//    gslpp::complex co = GF / 4. / M_PI * MW * myGTHDM.getCKM().computelamt_s();
//    double mHp2=myGTHDM.getmHp2();
//    double xHW=mHp2/(MW*MW);
//    double xtH=xt/xHW;
//    double mb=myGTHDM.getQuarks(QCD::BOTTOM).getMass();
//    double sd = (myGTHDM.getNd_11()).real();
//    double su = (myGTHDM.getNu_11()).real();
////Taken from 1006.0470:
//    double C1bsSRR = 4.0*mb*mb*xt*xtH*sd*su/(mHp2*(xtH-1.0)*(xtH-1.0)*(xtH-1.0)) //sd*su = sigd*sigu.conjugate()
//                     * (sd*su*(2.0*(xtH-1.0)-(xtH+1.0)*log(xtH)) //sd*su = sigd*sigu.conjugate()
//                        + (2.0*xt*xt*(xtH-1.0)*(xtH-1.0)*(xtH-1.0)*log(xt)/(xt-1.0)
//                           +2.0*xt*(xtH-1.0)*((xt-xtH)*(xtH-1.0)+(xtH-xt*xtH)*log(xtH)))/((xt-1.0)*(xt-xtH)*xtH));
//    
//    vmcdsp = StandardModelMatching::CMdbs2();
//    mcdbsp2.setMu(Mut);
//
//    switch (mcdbsp2.getOrder()) {
//        case NNLO:
//        case NLO:
//        case LO:
//            mcdbsp2.setCoeff(1, co * co * C1bsSRR, LO);
//            break;
//        default:
//            std::stringstream out;
//            out << mcdbsp2.getOrder();
//            throw std::runtime_error("THDMMatching::CMdbs2(): order " + out.str() + "not implemented");
//    }
//
//    vmcdsp.push_back(mcdbsp2);
//    return(vmcdsp);
//}


/*******************************************************************************
 * Wilson coefficients calculus, Bs -> mu+mu-.  From 1404.5865              *  
 * ****************************************************************************/


double GeneralTHDMMatching::C10Bll(double xt, double xHp, gslpp::complex su) {
    
    double C10 = su.abs2()*xt*xt/8*(1/(xHp-xt) + xHp/((xHp-xt)*(xHp-xt))*(log(xt)- log(xHp)));
    return C10;
    }

 gslpp::complex  GeneralTHDMMatching::CSboxBll(double xt, double xHp, gslpp::complex su, gslpp::complex sd, gslpp::complex sl) {
   
      gslpp::complex  CSboxU = xt/(8*(xHp-xt))*(sl*su.conjugate()*(xt/(xt-1)*log(xt)
             - xHp/(xHp-1)*log(xHp)) 
     + su*sl.conjugate()*(1-(xHp-xt*xt)/((xHp-xt)*(xt-1))*log(xt) 
     - (xHp*(xt-1))/((xHp-xt)*(xHp-1))*log(xHp))  + 2*sd*sl.conjugate()*(log(xt)-log(xHp))); 
             
      return(CSboxU);
    }


 gslpp::complex  GeneralTHDMMatching::CPboxBll(double xt, double xHp, gslpp::complex su, gslpp::complex sd, gslpp::complex sl)  {

    
        gslpp::complex  CPboxU = -xt/(8*(xHp-xt))*(-sl*su.conjugate()*(xt/(xt-1)*log(xt) - xHp/(xHp-1)*log(xHp)) 
     + su*sl.conjugate()*(1-(xHp-xt*xt)/((xHp-xt)*(xt-1))*log(xt) 
     - (xHp*(xt-1))/((xHp-xt)*(xHp-1))*log(xHp))  + 2*sd*sl.conjugate()*(log(xt)-log(xHp))); 
      
        return CPboxU;
    }

 
  gslpp::complex  GeneralTHDMMatching:: CPZUBll(double xt, double xHp, double sW2, gslpp::complex su, gslpp::complex sd) {

    
         // CPZF. Z-penguins diagrams. Eq. (52)
       
     gslpp::complex CPZF = (xt/(4*(xHp - xt)*(xHp - xt)))*(sd*su.conjugate()*(-((xt + xHp)/2) 
             +  ((xt*xHp)/(xHp - xt))*(log(xHp) - log(xt))) + 
           su.abs2()*(1/(6*(xHp - xt)))*((xHp*xHp - 8*xHp*xt -  17*xt*xt)/ 6 
             + ((xt*xt*(3*xHp + xt))/(xHp - xt))*(log(xHp) - log(xt)))) 
             + sW2*(xt/(6*(xHp - xt)*(xHp - xt)))*(sd*su.conjugate()*((5*xt - 3*xHp)/2 
             + ((xHp*(2*xHp - 3*xt))/(xHp - xt))*(log(xHp) - log(xt))) - 
            su.abs2()*(1/(6*(xHp - xt)))*(((4*xHp*xHp*xHp - 12*xHp*xHp*xHp*xt + 9*xHp*xt*xt + 3*xt*xt*xt)/(xHp - xt))* (log(xHp) - log(xt))
             - (17*xHp*xHp - 64*xHp*xt + 71.0*xt*xt)/6));
    
     //CPGBF. Goldstone penguin diagrams. Eq. (53)
     
    gslpp::complex  CPGBF = su.abs2()*(1 - sW2)*(xt*xt/(4*(xHp - xt)*(xHp - xt)))*
      (xHp*(log(xHp) - log(xt)) + xt - xHp);
    
      //CPZU. Z-penguin diagrams in the unitary gauge. Eq. (54)
    
    gslpp::complex  CPZU =  CPZF + CPGBF;
    
    return CPZU;
      
  }
 
double GeneralTHDMMatching::f1(double xHp, double xt){


      double f1 = (1.0/(2.0*(xHp - xt)))*(-xHp + xt + xHp*log(xHp) - 
     xt*log(xt));
      return f1;
      
  }


double GeneralTHDMMatching::f2(double xHp, double xt){

    double f2 =  (1.0/(2.0*(xHp - xt)))*(xt - ((xHp*xt)/(xHp - xt))*(log(xHp) - log(xt))); 
         
    return f2 ;
      
  }

double GeneralTHDMMatching::f3(double xHp, double xt){

       double f3= (1.0/(2.0*(xHp - xt)))*
    (xHp - (xHp*xHp*log(xHp))/(xHp - xt) + (xt*(2*xHp - xt)*log(xt))/
      (xHp - xt)); 
 
      return f3;
      
  }

double GeneralTHDMMatching::f4(double xHp, double xt){


    double f4 = (1.0/(4*(xHp - xt)*(xHp - xt)))*((xt*(3.0*xHp - xt))/2.0 - 
     ((xHp*xHp*xt)/(xHp - xt))*(log(xHp) - log(xt)));   
      return f4;
      
  }

double GeneralTHDMMatching::f5(double xHp, double xt){


    double f5 = (1.0/(4*(xHp - xt)*(xHp - xt)))*((xt*(xHp - 3.0*xt))/2.0 - 
    ((xHp*xt*(xHp - 2.0*xt))/(xHp - xt))*(log(xHp) - log(xt)));  
      return f5;
      
  }

double GeneralTHDMMatching::f6(double xHp, double xt){


    double f6 = (1.0/(2.0*(xHp - xt)))*
    ((xt*(xt*xt - 3.0*xHp*xt + 9*xHp - 5*xt - 2.0))/(4*(xt - 1.0)*(xt - 1.0)) + 
     ((xHp*(xHp*xt - 3.0*xHp + 2.0*xt))/(2.0*(xHp - 1.0)*(xHp - xt)))*
      log(xHp) + ((xHp*xHp*(-2.0*xt*xt*xt + 6*xt*xt - 9*xt + 2.0) + 
        3.0*xHp*xt*xt*(xt*xt - 2.0*xt + 3.0) - xt*xt*(2.0*xt*xt*xt - 3.0*xt*xt + 
          3.0*xt + 1.0))/(2.0*(xt - 1.0)*(xt - 1.0)*(xt - 1.0)*(xHp - xt)))*log(xt)); 
 
      return f6;
      
  }

double GeneralTHDMMatching::f7(double xHp, double xt){


     double f7= (1.0/(2.0*(xHp - xt)))*
    (((xt*xt + xt - 8)*(xHp - xt))/(4*(xt - 1.0)*(xt - 1.0)) - 
     ((xHp*(xHp + 2.0))/(2.0*(xHp - 1.0)))*log(xHp) + 
     ((xHp*(xt*xt*xt - 3.0*xt*xt + 3.0*xt + 2.0) + 3.0*(xt - 2.0)*xt*xt)/
       (2.0*(xt - 1.0)*(xt - 1.0)*(xt - 1.0)))*log(xt));  
      return f7;
      
  }

double GeneralTHDMMatching::f8(double xHp, double xt){

     double f8 = (1.0/(4*(xHp - xt)))*((xt*log(xt))/(xt - 1.0) - 
     (xHp*log(xHp))/(xHp - 1.0)); 
 
      return f8;
      
  }

double GeneralTHDMMatching::f9(double xHp, double xt){


     double f9 = (1.0/(8*(xHp - xt)))*(xHp/(xHp - 1.0) + 
     (xt*xt*log(xt))/((xt - 1.0)*(xHp - xt)) - 
     ((xHp*(xHp*xt + xHp - 2.0*xt))/((xHp - 1.0)*(xHp - 1.0)*(xHp - xt)))*
      log(xHp));  
      return f9;
      
  }

double GeneralTHDMMatching::f10(double xHp, double xt){
double f10 = (1.0/(8*(xHp - xt)))*
    ((xHp - xt)/((xHp - 1.0)*(xt - 1.0)) + ((xt*(xt - 2.0))/(xt - 1.0)*(xt - 1.0))*
      log(xt) - ((xHp*(xHp - 2.0))/(xHp - 1.0)*(xHp - 1.0))*log(xHp)); 
     
      return f10;
      
  }
  
  
gslpp::complex  GeneralTHDMMatching::g0(double xHp, double xt, gslpp::complex su, gslpp::complex sd){

    gslpp::complex  g0 = (1.0/(4.0*(xHp-xt)))*((sd*su.conjugate()*(xt/(xHp-xt)*(log(xHp)-log(xt))-1.0))
            + su.abs2()*(xt*xt/(2.0*(xHp-xt)*(xHp-xt))*(log(xHp)-log(xt)) + (xHp-3.0*xt)/(4.0*(xHp-xt))));
    
    return g0;
}

gslpp::complex  GeneralTHDMMatching::g1a(double xHp, double xt, gslpp::complex su, gslpp::complex sd){
    
     gslpp::complex  g1a = -(3.0/4.0) + sd*su.conjugate()*(xt/(xHp - xt))*(1.0 - (xHp/(xHp - xt))*(log(xHp) - log(xt))) 
    + su.abs2()*(xt/(2.0*(xHp - xt)*(xHp - xt)))*((xHp + xt)/2.0 - ((xHp*xt)/(xHp - xt))*(log(xHp) - log(xt)));
    
     return g1a;
    
}

gslpp::complex  GeneralTHDMMatching::g2a(double xHp, double xt, gslpp::complex su, gslpp::complex sd){
    
     gslpp::complex  g2a = sd*sd*su.conjugate()*f1(xHp, xt) + sd*su.conjugate()*su.conjugate()*f2(xHp, xt)+sd*su.abs2()*f3(xHp, xt) 
    + su*su.abs2()*f4(xHp, xt) - su.conjugate()*su.abs2()*f5(xHp, xt) + su*f6(xHp, xt) - su.conjugate()*f7(xHp, xt) +sd*f1(xHp, xt); 
   
     return g2a;
    
}

gslpp::complex  GeneralTHDMMatching::g3a(double xHp, double xt, gslpp::complex su, gslpp::complex sd){
    
    gslpp::complex g3a =sd*sd*su.conjugate()*f1(xHp, xt) - sd*su.conjugate()*su.conjugate()*f2(xHp, xt) 
    + sd*su.abs2()*f3(xHp, xt) + su*su.abs2()*f4(xHp, xt) + su.conjugate()*su.abs2()*f5(xHp, xt) 
    + su* f6(xHp, xt) + su.conjugate()*f7(xHp, xt) + sd*f1(xHp, xt); 
    
    return g3a;
}


    gslpp::complex  GeneralTHDMMatching::lambdaHHphi(double lambda3, double Relambda7,double Imlambda7, double Ri1, double Ri2, double Ri3 ){
         gslpp::complex  lambdaHHphi = lambda3*Ri1 + Relambda7*Ri2 - Imlambda7*Ri3;
        return lambdaHHphi;
    }
            
            


gslpp::complex  GeneralTHDMMatching::CphiU(double xHp, double xt, double vev, double xphi, double mu, double Ri1, double Ri2, double Ri3, double mHi_2, double lambda3, double Relambda7,double Imlambda7, gslpp::complex su, gslpp::complex sd){
     gslpp::complex i = gslpp::complex::i();
     gslpp::complex CphiU = xt*((1/(2*xphi))*(su-sd)*(1 + su.conjugate()*sd)*
      (Ri2 + i*Ri3)*mu + (vev*vev/mHi_2)*lambdaHHphi(lambda3, Relambda7, Imlambda7, Ri1, Ri2, Ri3)*g0(xHp, xt, su, sd) + Ri1*((1/(2.0*xphi))*g1a(xHp, xt, su, sd)) + Ri2*(1/(2.0*xphi)*g2a(xHp, xt, su, sd)) + i*Ri3*(1/(2.0*xphi))*g3a(xHp, xt, su, sd));
    
     return CphiU;
}


std::vector<WilsonCoefficient>& GeneralTHDMMatching::CMBMll(QCD::lepton lepton)
{    
    updateGTHDMParameters();

            
     //From 1404.5865
        
     //complex i 
    gslpp::complex i = gslpp::complex::i();
       

    double Muw = myGTHDM.getMuw();
  //  double Mut = myGTHDM.getMut();
    double mHp2 = myGTHDM.getmHp2();
    double MW = myGTHDM.Mw();
    double Mt_muw = myGTHDM.Mrun(Muw, myGTHDM.getQuarks(QCD::TOP).getMass_scale(), 
                        myGTHDM.getQuarks(QCD::TOP).getMass(), FULLNNLO);
    /*double mt_mt = myGTHDM.Mrun(Mut, myGTHDM.getQuarks(QCD::TOP).getMass_scale(), 
                        myGTHDM.getQuarks(QCD::TOP).getMass(), FULLNNLO);*/
   double mb=myGTHDM.getQuarks(QCD::BOTTOM).getMass();
   
   double ml=myGTHDM.getLeptons(lepton).getMass();

  //  mlep = SM.getLeptons(lep).getMass();

    
    double xt = (Mt_muw*Mt_muw)/(MW*MW);
    double xHp = (mHp2)/(MW*MW);
    double vev = myGTHDM.v();
    double sW2 = myGTHDM.sW2();
    double mHl = myGTHDM.getMHl();
    double mH1_2=mHl*mHl;
   
    
     //mu contains the missalignemtn dependece. It should be mu -> CR(mu0) - log(mu/mu0). Eq (22)

   

    // mu relates the high alignment scale (\Lambda), not know (which can be until the Plank scale), to the electroweak scale
    //double mu = log(Lambda/MW); 
    // Lambda is not a paraemeter of the model so now mu is set to 0 -> the model is aligned at the ew scale
    double mu = 0;
  
    double Imlambda7=myGTHDM.getMyGTHDMCache()->Imlambda7H;
    double Relambda7=myGTHDM.getMyGTHDMCache()->Relambda7H;
    double lambda3=myGTHDM.getMyGTHDMCache()->lambda3H;
    
    double  R11 = myGTHDM.getMyGTHDMCache()->R11_GTHDM;
    double  R12 = myGTHDM.getMyGTHDMCache()->R12_GTHDM;
    double  R13 = myGTHDM.getMyGTHDMCache()->R13_GTHDM;
    double  R21 = myGTHDM.getMyGTHDMCache()->R21_GTHDM;
    double  R22 = myGTHDM.getMyGTHDMCache()->R22_GTHDM;
    double  R23 = myGTHDM.getMyGTHDMCache()->R23_GTHDM;
    double  R31 = myGTHDM.getMyGTHDMCache()->R31_GTHDM;
    double  R32 = myGTHDM.getMyGTHDMCache()->R32_GTHDM;
    double  R33 = myGTHDM.getMyGTHDMCache()->R33_GTHDM;

   gslpp::complex sl = myGTHDM.getNl_11();
   gslpp::complex su = myGTHDM.getNu_11();
   gslpp::complex sd = myGTHDM.getNd_11();
    


    //Mass of the physical scalars

   double mH2_2 = myGTHDM.getmH2sq();
   double mH3_2 = myGTHDM.getmH3sq();
    
   
    double xphi1 = mH1_2/(MW*MW);
    double xphi2 = mH2_2/(MW*MW);
    double xphi3 = mH3_2/(MW*MW);
       
    //Yukawa couplings. Eq. (19)
    
    gslpp::complex yl1 = 0.0;
    gslpp::complex yl2 = 0.0;
    gslpp::complex yl3 = R31  + (R32 + i*R33)*sl;
    
   if(myGTHDM.getSMHiggs()){
        yl1 = R11  + (R12 + i*R13)*sl;  
        yl2 = R21  + (R22 + i*R23)*sl;
   }
    else{
        yl2 = R11  + (R12 + i*R13)*sl;  
        yl1 = R21  + (R22 + i*R23)*sl;
   }
     
      
    gslpp::complex CSboxU = CSboxBll(xt,  xHp, su, sd, sl);
    gslpp::complex CPboxU = CPboxBll(xt,  xHp, su, sd, sl);
    gslpp::complex CPZU = CPZUBll(xt,  xHp,  sW2, su, sd);
    
    gslpp::complex CSphi1U = 0.0;
    gslpp::complex CSphi2U = 0.0;
    gslpp::complex CSphi3U = yl3.real()*CphiU(xHp,  xt,  vev,  xphi3,  mu,  R31,  R32,  R33,  mH3_2,  lambda3,  Relambda7, Imlambda7, su, sd);

    gslpp::complex CPphi1U = 0.0;
    gslpp::complex CPphi2U = 0.0;
    gslpp::complex CPphi3U = i*yl3.imag()*CphiU(xHp,  xt,  vev,  xphi3,  mu,  R31,  R32,  R33,  mH3_2,  lambda3,  Relambda7, Imlambda7, su, sd);
    
    
    if(myGTHDM.getSMHiggs()){
        CSphi1U = yl1.real()*CphiU(xHp,  xt,  vev,  xphi1,  mu,  R11,  R12,  R13,  mH1_2,  lambda3,  Relambda7, Imlambda7, su, sd);
        CSphi2U = yl2.real()*CphiU(xHp,  xt,  vev,  xphi2,  mu,  R21,  R22,  R23,  mH2_2,  lambda3,  Relambda7, Imlambda7, su, sd);
        CPphi1U = i*yl1.imag()*CphiU(xHp,  xt,  vev,  xphi1,  mu,  R11,  R12,  R13,  mH1_2,  lambda3,  Relambda7, Imlambda7, su, sd);
        CPphi2U = i*yl2.imag()*CphiU(xHp,  xt,  vev,  xphi2,  mu,  R21,  R22,  R23,  mH2_2,  lambda3,  Relambda7, Imlambda7, su, sd);
    }
    else{
        CSphi1U = yl1.real()*CphiU(xHp,  xt,  vev,  xphi1,  mu,  R21,  R22,  R23,  mH1_2,  lambda3,  Relambda7, Imlambda7, su, sd);
        CSphi2U = yl2.real()*CphiU(xHp,  xt,  vev,  xphi2,  mu,  R11,  R12,  R13,  mH2_2,  lambda3,  Relambda7, Imlambda7, su, sd);
        CPphi1U = i*yl1.imag()*CphiU(xHp,  xt,  vev,  xphi1,  mu,  R21,  R22,  R23,  mH1_2,  lambda3,  Relambda7, Imlambda7, su, sd);
        CPphi2U = i*yl2.imag()*CphiU(xHp,  xt,  vev,  xphi2,  mu,  R11,  R12,  R13,  mH2_2,  lambda3,  Relambda7, Imlambda7, su, sd);
    }
    
    
    //Total 2HDM Wilson coefficients CS and CP PART. Eq. (31)-(33) without SM part
   
    gslpp::complex CSphiU =  CSboxU + CSphi1U + CSphi2U + CSphi3U;
    gslpp::complex CPphiU =  CPboxU + CPZU + CPphi1U + CPphi2U + CPphi3U;

    
    vmcBMll = StandardModelMatching::CMBMll(lepton);
      switch (mcbsg.getScheme()) {
        case NDR:

            break;
        default:
        std::stringstream out;
       out << mcBMll.getScheme();
     throw std::runtime_error("GeneralTHDMMatching::CMBMll(): scheme " + out.str() + "not implemented");
    }
    mcBMll.setMu(Muw);
    
     switch (mcBMll.getOrder()) {
        case NNLO:
        case NLO:           
        case LO:           
            mcBMll.setCoeff(9 , C10Bll(xt, xHp, su)/(sW2), LO);
            mcBMll.setCoeff(10 , (CSphiU*mb*ml)/(MW*MW*sW2), LO);
            mcBMll.setCoeff(11 , (CPphiU*mb*ml)/(MW*MW*sW2), LO);
            break;
        default:
            std::stringstream out;
            out << mcBMll.getOrder();
            throw std::runtime_error("GeneralTHDMMatching::CMBMll(): order " + out.str() + "not implemeted"); 
            }
    vmcBMll.push_back(mcBMll);
    return (vmcBMll);
}


/*******************************************************************************
 * Wilson coefficients calculus, Misiak base for B -> \tau \nu                  *  
 * ****************************************************************************/

std::vector<WilsonCoefficient>& GeneralTHDMMatching::CMbtaunu(QCD::meson meson_i) {

    
    if (!myGTHDM.getATHDMflag())
    {
        throw std::runtime_error("CMbtaunu is only available in the ATHDM at the moment.");
        return (vmcbsg);
    }
    else{
        
    
    double Muw = myGTHDM.getMuw();
    double GF = myGTHDM.getGF();
    myCKM = myGTHDM.getVCKM();
    double mB = myGTHDM.getMesons(meson_i).getMass();
    gslpp::complex zetad=myGTHDM.getNd_11();
    gslpp::complex zetau=myGTHDM.getNu_11();
    gslpp::complex zetal=myGTHDM.getNl_11();
    double mHp2=myGTHDM.getmHp2();
    
    double mb = myGTHDM.getQuarks(QCD::BOTTOM).getMass();
    double mu = myGTHDM.getQuarks(QCD::UP).getMass();
    double mc = myGTHDM.getQuarks(QCD::CHARM).getMass();
    

    vmcbtaunu = StandardModelMatching::CMbtaunu(meson_i);
    mcbtaunu.setMu(Muw);
 
    switch (mcbtaunu.getOrder()) {
        case NNLO:
        case NLO:
        case LO:
            if (meson_i == QCD::B_P) mcbtaunu.setCoeff(0, -4.*GF * myCKM(0,2) / sqrt(2.)*(mB*mB*zetal.conjugate()*(mb*zetad+mu*zetau)/mHp2/(mb+mu)) , LO);
            else if (meson_i == QCD::B_C) mcbtaunu.setCoeff(0, -4.*GF * myCKM(1,2) / sqrt(2.)*(mB*mB*zetal.conjugate()*(mb*zetad+mc*zetau)/mHp2/(mb+mc)), LO);
            else throw std::runtime_error("GeneralTHDMMatching::CMbtaunu(): meson not implemeted"); 
            break;
        default:
            std::stringstream out;
            out << mcbtaunu.getOrder();
            throw std::runtime_error("GeneralTHDMMatching::CMbtaunu(): order " + out.str() + "not implemented");
    }
    
    vmcbtaunu.push_back(mcbtaunu);
    
    return(vmcbtaunu);
    }
}



/*******************************************************************************
 * Wilson coefficients calculus, LEFT base for D -> \lepton \nu                  *  
 * ****************************************************************************/

std::vector<WilsonCoefficient>& GeneralTHDMMatching::CMcleptonnu(QCD::meson meson_i, QCD::lepton lepton_i) {

    if (!myGTHDM.getATHDMflag())
    {
        throw std::runtime_error("CMcleptonnu is only available in the ATHDM at the moment.");
        return (vmcbsg);
    }
    else{
        
    
    double Muw = myGTHDM.getMuw();
    double GF = myGTHDM.getGF();
    myCKM = myGTHDM.getVCKM();
    double mD = myGTHDM.getMesons(meson_i).getMass();
    gslpp::complex zetad=myGTHDM.getNd_11();
    gslpp::complex zetal=myGTHDM.getNl_11();
    gslpp::complex zetau=myGTHDM.getNu_11();
    double mHp2=myGTHDM.getmHp2();
    
    double mc = myGTHDM.getQuarks(QCD::CHARM).getMass();
    double md = myGTHDM.getQuarks(QCD::DOWN).getMass();
    double ms = myGTHDM.getQuarks(QCD::STRANGE).getMass();
    
    double mtau = myGTHDM.getLeptons(StandardModel::TAU).getMass();
    double mmuon = myGTHDM.getLeptons(StandardModel::MU).getMass();
    
    

    vmccleptonnu = StandardModelMatching::CMcleptonnu(meson_i, lepton_i);
    mccleptonnu.setMu(Muw);
 
    switch (mccleptonnu.getOrder()) {
        case NNLO:
        case NLO:
        case LO:
            
            switch (meson_i){
                case QCD::D_P:
                    mccleptonnu.setCoeff(0, -4.*GF * myCKM(1,0) / sqrt(2.)*(mD*mD*zetal.conjugate()*(md*zetad+mc*zetau)/mHp2/(mc+md)), LO);
                break;
                case QCD::D_S:
                    mccleptonnu.setCoeff(0, -4.*GF * myCKM(1,1) / sqrt(2.)*(mD*mD*zetal.conjugate()*(ms*zetad+mc*zetau)/mHp2/(mc+ms)), LO);
                break;
                // We need to include also D_s, for that we need to include that meson in the SM
                default:
                    throw std::runtime_error("GeneralTHDMMatching::CMcleptonnu(): doesn't include that meson");
                
            }
            
            
            
            break;
        default:
            std::stringstream out;
            out << mccleptonnu.getOrder();
            throw std::runtime_error("GeneralTHDMMatching::CMcleptonnu(): order " + out.str() + "not implemented");
    }
    
    vmccleptonnu.push_back(mccleptonnu);
    
    return(vmccleptonnu);
    
    }
}




 /*******************************************************************************
 * Wilson coefficients calculus, LEFT basis [1709.04486] for (K -> lepton nu) from [1706.00410]                *
 * ****************************************************************************/ 
  std::vector<WilsonCoefficient>& GeneralTHDMMatching::CMsleptonnu(QCD::meson meson_i, QCD::lepton lepton_i) 
{
    
    if (!myGTHDM.getATHDMflag())
    {
        throw std::runtime_error("CMsleptonnu is only available in the ATHDM at the moment.");
        return (vmcbsg);
    }
    else{
        
    
    double Muw = myGTHDM.getMuw();
    double GF = myGTHDM.getGF();
    myCKM = myGTHDM.getVCKM();
    double mK = myGTHDM.getMesons(meson_i).getMass();
    gslpp::complex zetad=myGTHDM.getNd_11();
    gslpp::complex zetal=myGTHDM.getNl_11();
    gslpp::complex zetau=myGTHDM.getNu_11();
    double mHp2=myGTHDM.getmHp2();
    
    double mu = myGTHDM.getQuarks(QCD::UP).getMass();
    double ms = myGTHDM.getQuarks(QCD::STRANGE).getMass();  
      
      
    vmcsleptonnu = StandardModelMatching::CMsleptonnu(meson_i, lepton_i);
    mcsleptonnu.setMu(Muw);
 
    switch (mcsleptonnu.getOrder()) {
        case NNLO:
        case NLO:
        case LO:
            switch (meson_i){
                case QCD::K_P:
                    mcsleptonnu.setCoeff(0, -4.*GF * myCKM(0,1) / sqrt(2.)*(mK*mK*zetal.conjugate()*(ms*zetad+mu*zetau)/mHp2/(ms+mu)), LO);
                break;
                default:
                throw std::runtime_error("GeneralTHDMMatching::CMsleptonnu(): doesn't include that meson");
                
            }
           
        break;
        default:
            std::stringstream out;
            out << mcsleptonnu.getOrder();
            throw std::runtime_error("GeneralTHDMMatching::CMsleptonnu(): order " + out.str() + "not implemented");
    }
    
    
    
    vmcsleptonnu.push_back(mcsleptonnu);
    return(vmcsleptonnu);
    
    }
}     

  
  

  
 /*******************************************************************************
 * Wilson coefficients calculus, LEFT basis [1709.04486] for (\pi -> lepton nu) from [1706.00410]                *
 * ****************************************************************************/ 
  std::vector<WilsonCoefficient>& GeneralTHDMMatching::CMuleptonnu(QCD::meson meson_i, QCD::lepton lepton_i) 
{
      
      
    if (!myGTHDM.getATHDMflag())
    {
        throw std::runtime_error("CMuleptonnu is only available in the ATHDM at the moment.");
        return (vmcbsg);
    }
    else{
        
    
    double Muw = myGTHDM.getMuw();
    double GF = myGTHDM.getGF();
    myCKM = myGTHDM.getVCKM();
    double mP = myGTHDM.getMesons(meson_i).getMass();
    gslpp::complex zetad=myGTHDM.getNd_11();
    gslpp::complex zetal=myGTHDM.getNl_11();
    gslpp::complex zetau=myGTHDM.getNu_11();
    double mHp2=myGTHDM.getmHp2();
    
    double mu = myGTHDM.getQuarks(QCD::UP).getMass();
    double md = myGTHDM.getQuarks(QCD::DOWN).getMass();

    
    
    vmculeptonnu = StandardModelMatching::CMuleptonnu(meson_i, lepton_i);
    
    mculeptonnu.setMu(Muw);
 
    switch (mculeptonnu.getOrder()) {
        case NNLO:
        case NLO:
        case LO:
            switch (meson_i){
                case QCD::P_P:
                    mculeptonnu.setCoeff(0, -4.*GF * myCKM(0,0) / sqrt(2.)*(mP*mP*zetal.conjugate()*(md*zetad+mu*zetau)/mHp2/(md+mu)) , LO);
                break;
                default:
                throw std::runtime_error("GeneralTHDMMatching::CMuleptonnu(): doesn't include that meson");
                
            }
           
        break;
        default:
            std::stringstream out;
            out << mculeptonnu.getOrder();
            throw std::runtime_error("GeneralTHDMMatching::CMuleptonnu(): order " + out.str() + "not implemented");
    }
    
    vmculeptonnu.push_back(mculeptonnu);
    return(vmculeptonnu);
    
    
    }
    
    
}  
  
  
  
  





std::vector<WilsonCoefficient>& GeneralTHDMMatching::CMbsg() 
{
    vmcbsg = StandardModelMatching::CMbsg();

    if (!myGTHDM.getATHDMflag())
    {
        throw std::runtime_error("bsgamma is only available in the ATHDM at the moment.");
        return (vmcbsg);
    }

    double Muw = myGTHDM.getMuw();
    
    double Mut = myGTHDM.getMut();
    double mt = myGTHDM.Mrun(Mut, myGTHDM.getQuarks(QCD::TOP).getMass_scale(), 
                        myGTHDM.getQuarks(QCD::TOP).getMass(), FULLNNLO);
    double mHp=myGTHDM.getmHp();
    
    gslpp::complex sigmau = myGTHDM.getNu_11();
    gslpp::complex sigmad = myGTHDM.getNd_11();
    
    gslpp::complex co = 1.; // (- 4. * GF / sqrt(2)) * SM.computelamt_s(); THIS SHOULD ALREADY BE IMPLEMENTED IN THE OBSERVABLE 
    mcbsg.setMu(Muw);
    
    switch (mcbsg.getOrder()) {
        case NNLO:
            for (int j=0; j<8; j++){
            mcbsg.setCoeff(j, co * myGTHDM.Alstilde5(Muw) * myGTHDM.Alstilde5(Muw) * setWCbsg(j, sigmau, sigmad, mt, mHp, Muw, NNLO), NNLO);
            }
        case NLO:
            for (int j=0; j<8; j++){
            mcbsg.setCoeff(j, co * myGTHDM.Alstilde5(Muw) * setWCbsg(j, sigmau, sigmad, mt, mHp, Muw, NLO), NLO);
            }
        case LO:
            for (int j=0; j<8; j++){
            mcbsg.setCoeff(j, co * setWCbsg(j, sigmau, sigmad, mt, mHp, Muw, LO), LO);
            }
            break;
        default:
            std::stringstream out;
            out << mcbsg.getOrder();
            throw std::runtime_error("THDMMatching::CMbsg(): order " + out.str() + "not implemented"); 
    }
    
    vmcbsg.push_back(mcbsg);
    return(vmcbsg);
}

/*******************************************************************************
 * Wilson coefficients calculus, Misiak base for b -> s gamma                  *  
 * ****************************************************************************/

gslpp::complex GeneralTHDMMatching::setWCbsg(int i, gslpp::complex sigmau, gslpp::complex sigmad, double mt, double mhp, double mu, orders order)
{
    if ( su.abs() == sigmau.abs() && su.arg() == sigmau.arg() &&
         sd.abs() == sigmad.abs() && sd.arg() == sigmad.arg() &&
         mtbsg == mt  && mhpbsg == mhp && mubsg == mu){
        switch (order){
        case NNLO:
            return ( CWbsgArrayNNLO[i] );
        case NLO:
            return ( CWbsgArrayNLO[i] );
            break;
        case LO:
            return ( CWbsgArrayLO[i] );
            break;
        default:
            std::stringstream out;
            out << order;
            throw std::runtime_error("order" + out.str() + "not implemeted");
        }
    }

    su = sigmau; sd = sigmad; mtbsg = mt; mhpbsg = mhp; mubsg = mu;

    double x = mt*mt/mhp/mhp;

    double x2 = x*x;
    double x3 = x2*x;
    double x4 = x3*x;
    double x5 = x4*x;
    double xm = 1. - x;
    double xm2 = xm*xm;
    double xm3 = xm2*xm;
    double xm4 = xm3*xm;
    double xm6 = xm4*xm2;
    double xm8 = xm4*xm4;
    double xo = 1. - 1./x;
    double xo2 = xo*xo;
    double xo4 = xo*xo2*xo;
    double xo6 = xo4*xo2;
    double xo8 = xo4*xo4;
    double Lx = log(x);
    double Lx2 = Lx*Lx;
    double Lx3 = Lx2*Lx;
    double Li2 = gsl_sf_dilog(1.-1./x);
    
    double abssu = su.abs();
    double abssu2 = abssu*abssu;
    gslpp::complex susd = su.conjugate()*sd;

    double lstmu = 2. * log(mu/mt);

    double n70ct = - ( (7. - 5.*x - 8.*x2)/(36.*xm3) + (2.*x - 3.*x2)*Lx/(6.*xm4) ) * x/2.;
    double n70fr = ( (3.*x - 5.*x2)/(6.*xm2) + (2.*x - 3.*x2)*Lx/(3.*xm3) ) / 2.;

    double n80ct = - ( (2. + 5.*x - x2)/(12.*xm3) + (x*Lx)/(2.*xm4) ) * x/2.;
    double n80fr = ( (3.*x - x2)/(2.*xm2) + (x*Lx)/xm3 ) / 2.;

    double n41ct = (-16.*x + 29.*x2 - 7.*x3)/(36.*xm3) + (-2.*x + 3.*x2)*Lx/(6.*xm4);

    double n71ct = (797.*x - 5436.*x2 + 7569.*x3 - 1202.*x4)/(486.*xm4) + 
                    (36.*x2 - 74.*x3 + 16.*x4)*Li2/(9.*xm4) + 
                    ((7.*x - 463.*x2 + 807.*x3 - 63.*x4)*Lx)/(81.*xm3*xm2);
    double cd71ct = (-31.*x - 18.*x2 + 135.*x3 - 14.*x4)/(27.*xm4) + (-28.*x2 + 46.*x3 + 6.*x4)*Lx/(9.*xm3*xm2);
    double n71fr = (28.*x - 52.*x2 + 8.*x3)/(3.*xm3) + (-48.*x + 112.*x2 - 32.*x3)*Li2/(9.*xm3) + 
                    (66.*x - 128.*x2 + 14.*x3)*Lx/(9.*xm4);
    double cd71fr = (42.*x - 94.*x2 + 16.*x3)/(9.*xm3) + (32.*x - 56.*x2 - 12.*x3)*Lx/(9.*xm4);
    
    double n81ct = (1130.*x - 18153.*x2 + 7650.*x3 - 4451.*x4)/(1296.*xm4) + 
                    (30.*x2 - 17.*x3 + 13.*x4)*Li2/(6.*xm4) + 
                    (-2.*x - 2155.*x2 + 321.*x3 - 468.*x4)*Lx/(216.*xm3*xm2);
    double cd81ct = (-38.*x - 261.*x2 + 18.*x3 - 7.*x4)/(36.*xm4) + (-31.*x2 - 17.*x3)*Lx/(6.*xm3*xm2);
    double n81fr = (143.*x - 44.*x2 + 29.*x3)/(8.*xm3) + (-36.*x + 25.*x2 - 17.*x3)*Li2/(6.*xm3) + 
                    (165.*x - 7.*x2 + 34.*x3)*Lx/(12.*xm4);
    double cd81fr = (81.*x - 16.*x2 + 7.*x3)/(6.*xm3) + (19.*x + 17.*x2)*Lx/(3.*xm4);

    double n32ct = (10.*x4 + 30.*x2 - 20.*x)/(27.*xm4) * Li2 + 
                    (30.*x3 - 66.*x2 - 56.*x)/(81.*xm4) * Lx + (6.*x3 - 187.*x2 + 213.*x)/(-81.*xm3);
    double cd32ct = (-30.*x2 + 20.*x)/(27.*xm4)*Lx + (-35.*x3 + 145.*x2 - 80.*x)/(-81.*xm3);

    double n42ct = (515.*x4 - 906.*x3 + 99.*x2 + 182.*x)/(54.*xm4) * Li2 + 
                    (1030.*x4 - 2763.*x3 - 15.*x2 + 980.*x)/(-108.*xm3*xm2)*Lx + 
                    (-29467.*x4 + 68142.*x3 - 6717.*x2 - 18134.*x)/(1944.*xm4);
    double cd42ct = (-375.*x3 - 95.*x2 + 182.*x)/(-54.*xm3*xm2)*Lx + 
                    (133.*x4 - 108.*x3 + 4023.*x2 - 2320.*x)/(324.*xm4);

    double cd72ct = -(x * (67930.*x4 - 470095.*x3 + 1358478.*x2 - 700243.*x + 54970.))/(-2187.*xm3*xm2) + 
                    (x * (10422.*x4 - 84390.*x3 + 322801.*x2 - 146588.*x + 1435.))/(729.*xm6)*Lx +
                    (2.*x2 * (260.*x3 - 1515.*x2 + 3757.*x - 1446.))/(-27. * xm3*xm2) * Li2;
    double ce72ct = (x * (-518.*x4 + 3665.*x3 - 17397.*x2 + 3767.*x + 1843.))/(-162.*xm3*xm2) +
                    (x2 * (-63.*x3 + 532.*x2 + 2089.*x - 1118.))/(27.*xm6)*Lx;
    double cd72fr = -(x * (3790.*x3 - 22511.*x2 + 53614.*x - 21069.))/(81.*xm4) -
                    (2.*x * (-1266.*x3 + 7642.*x2 - 21467.*x + 8179.))/(-81.*xm3*xm2)*Lx +
                    (8.*x * (139.*x3 - 612.*x2 + 1103.*x - 342.))/(27.*xm4) * Li2;
    double ce72fr = -(x * (284.*x3 - 1435.*x2 + 4304.*x - 1425.))/(27.*xm4) -
                    (2.*x * (63.*x3 - 397.*x2 - 970.*x + 440.))/(-27.*xm3*xm2)*Lx;

    double cd82ct = -(x * (522347.*x4 - 2423255.*x3 + 2706021.*x2 - 5930609.*x + 148856))/(-11664.*xm3*xm2) + 
                    (x * (51948.*x4 - 233781.*x3 + 48634.*x2 - 698693.*x + 2452.))/(1944.*xm6)*Lx + 
                    (x2 * (481.*x3 - 1950.*x2 + 1523.*x - 2550.))/(-18.*xm3*xm2) * Li2;
    double ce82ct = (x * (-259.*x4 + 1117.*x3 + 2925.*x2 + 28411.*x + 2366.))/(-216.*xm3*xm2) -
                    (x2 * (139.*x2 + 2938.*x + 2683.))/(36.*xm6)*Lx;
    double cd82fr = -(x * (1463.*x3 - 5794.*x2 + 5543.*x - 15036.))/(27.*xm4) - 
                    (x * (-1887.*x3 + 7115.*x2 + 2519.*x + 19901.))/(-54.*xm3*xm2)*Lx -
                    (x * (-629.*x3 + 2178.*x2 - 1729.*x + 2196.))/(18.*xm4) * Li2;
    double ce82fr = -(x * (259.*x3 - 947.*x2 - 251.*x - 5973.))/(36.*xm4)-
                    (x * (139.*x2 + 2134.*x + 1183.))/(-18.*xm3*xm2)*Lx;
    
    double n72ct = 0.;
    if (mhp < 50.) 
        n72ct = -274.2/x5 - 72.13*Lx/x5 + 24.41/x4 - 168.3*Lx/x4 + 79.15/x3 - 
                103.8*Lx/x3 + 47.09/x2 - 38.12*Lx/x2 + 15.35/x - 8.753*Lx/x + 3.970;
    else if (mhp < mt)
        n72ct = 1.283 + 0.7158 * xo + 0.4119 * xo2 + 0.2629 * xo*xo2 + 0.1825 * xo4 + 
                0.1347 * xo*xo4 + 0.1040 * xo6 + 0.0831 * xo*xo6 + 0.06804 * xo8 + 
                0.05688 * xo*xo8 + 0.04833 * xo2*xo8 + 0.04163 * xo*xo2*xo8 + 0.03625 * xo4*xo8 + 
                0.03188 * xo*xo4*xo8 + 0.02827 * xo6*xo8 + 0.02525 * xo*xo6*xo8 + 0.02269 * xo8*xo8;
    else if (mhp < 520.)
        n72ct = 1.283 - 0.7158 * xm - 0.3039 * xm2 - 0.1549 * xm3 - 0.08625 * xm4 - 
                0.05020 * xm3*xm2 - 0.02970 * xm6 - 0.01740 * xm3*xm4 - 0.009752 * xm8 - 
                0.004877 * xm3*xm6 - 0.001721 * xm2*xm8 + 0.0003378 * xm3*xm8 + 0.001679 * xm4*xm8 + 
                0.002542 * xm*xm4*xm8 + 0.003083 * xm6*xm8 + 0.003404 * xm*xm6*xm8 + 0.003574 * xm8*xm8;
    else n72ct = -823.0*x5 + 42.30*x5*Lx3 - 412.4*x5*Lx2 - 3362*x5*Lx -
                1492*x4 - 23.26*x4*Lx3 - 541.4*x4*Lx2 - 2540*x4*Lx -
                1158*x3 - 34.50*x3*Lx3 - 348.2*x3*Lx2 - 1292*x3*Lx -
                480.9*x2 - 20.73*x2*Lx3 - 112.4*x2*Lx2 - 396.1*x2*Lx -
                8.278*x + 0.9225*x*Lx2 + 4.317*x*Lx;
    
    double n72fr = 0.;
    if (mhp < 50.) 
        n72fr = -( 194.3/x5 + 101.1*Lx/x5 - 24.97/x4 + 168.4*Lx/x4 - 78.90/x3 + 
                106.2*Lx/x3 - 49.32/x2 + 38.43*Lx/x2 - 12.91/x + 9.757*Lx/x + 8.088 );
    else if (mhp < mt)
        n72fr = -( 12.82 - 1.663 * xo - 0.8852 * xo2 - 0.4827 * xo*xo2 - 0.2976 * xo4 - 
                0.2021 * xo*xo4 - 0.1470 * xo6 - 0.1125 * xo*xo6 - 0.08931 * xo8 - 
                0.07291 * xo*xo8 - 0.06083 * xo2*xo8 - 0.05164 * xo*xo2*xo8 - 0.04446 * xo4*xo8 - 
                0.03873 * xo*xo4*xo8 - 0.03407 * xo6*xo8 - 0.03023 * xo*xo6*xo8 - 0.02702 * xo8*xo8 );
    else if (mhp < 400.)
        n72fr = -( 12.82 + 1.663 * xm + 0.7780 * xm2 + 0.3755 * xm3 + 0.1581 * xm4 + 
                0.03021 * xm3*xm2 - 0.04868 * xm6 - 0.09864 * xm3*xm4 - 0.1306 * xm8 - 
                0.1510 * xm3*xm6 - 0.1637 * xm2*xm8 - 0.1712 * xm3*xm8 - 0.1751 * xm4*xm8 - 
                0.1766 * xm*xm4*xm8 - 0.1763 * xm6*xm8 - 0.1748 * xm*xm6*xm8 - 0.1724 * xm8*xm8 );
    else n72fr = -( 2828 * x5 - 66.63 * x5*Lx3 + 469.4 * x5*Lx2 + 1986 * x5*Lx + 
                1480 * x4 + 36.08 * x4*Lx3 + 323.2 * x4*Lx2 + 169.9 * x4*Lx + 
                166.7 * x3 + 19.73 * x3*Lx3 - 46.61 * x3*Lx2 - 826.2 * x3*Lx - 
                524.1 * x2 - 8.889 * x2*Lx3 - 195.7 * x2*Lx2 - 870.3 * x2*Lx - 
                572.2 * x - 20.94 * x*Lx3 - 123.5 * x*Lx2 - 453.5 * x*Lx );
    
    double n82ct = 0.;
    if (mhp < 50.) 
        n82ct = 826.2/x5 - 300.7*Lx/x5 + 96.35/x4 + 91.89*Lx/x4 - 66.39/x3 +
                78.58*Lx/x3 - 39.76/x2 + 20.02*Lx/x2 - 5.214/x + 2.278;
    else if (mhp < mt)
        n82ct = 1.188 + 0.4078 * xo + 0.2002 * xo2 + 0.1190 * xo*xo2 + 0.07861 * xo4 + 
                0.05531 * xo*xo4 + 0.04061 * xo6 + 0.03075 * xo*xo6 + 0.02386 * xo8 + 
                0.01888 * xo*xo8 + 0.01520 * xo2*xo8 + 0.01241 * xo*xo2*xo8 + 0.01026 * xo4*xo8 + 
                0.008575 * xo*xo4*xo8 + 0.007238 * xo6*xo8 + 0.006164 * xo*xo6*xo8 + 0.005290 * xo8*xo8;
    else if (mhp < 600.)
        n82ct = 1.188 - 0.4078 * xm - 0.2076 * xm2 - 0.1265 * xm3 - 0.08570 * xm4 -
                0.06204 * xm3*xm2 - 0.04689 * xm6 - 0.03652 * xm3*xm4 - 0.02907 * xm8 -
                0.02354 * xm3*xm6 - 0.01933 * xm2*xm8 - 0.01605 * xm3*xm8 - 0.01345 * xm4*xm8 -
                0.01137 * xm*xm4*xm8 - 0.009678 * xm6*xm8 - 0.008293 * xm*xm6*xm8 - 0.007148 * xm8*xm8;
    else n82ct = -19606 * x5 - 226.7 * x5*Lx3 - 5251 * x5*Lx2 - 26090 * x5*Lx - 
                9016 * x4 - 143.4 * x4*Lx3 - 2244 * x4*Lx2 - 10102 * x4*Lx - 
                3357 * x3 - 66.32 * x3*Lx3 - 779.6 * x3*Lx2 - 3077 * x3*Lx - 
                805.5 * x2 - 22.98 * x2*Lx3 - 169.1 * x2*Lx2 - 602.7 * x2*Lx + 
                0.7437  * x + 0.6908 * x*Lx2 + 3.238 * x*Lx;
    
    double n82fr = 0.;
    if (mhp < 50.) 
        n82fr = -( -1003/x5 + 476.9*Lx/x5 - 205.7/x4 - 71.62*Lx/x4 + 62.26/x3 -
                110.7*Lx/x3 + 63.74/x2 - 35.42*Lx/x2 + 10.89/x - 3.174 );
    else if (mhp < mt)
        n82fr = -( -0.6110 - 1.095 * xo - 0.4463 * xo2 - 0.2568 * xo*xo2 - 0.1698 * xo4 - 
                0.1197 * xo*xo4 - 0.08761 * xo6 - 0.06595 * xo*xo6 - 0.05079 * xo8 - 
                0.03987 * xo*xo8 - 0.03182 * xo2*xo8 - 0.02577 * xo*xo2*xo8 - 0.02114 * xo4*xo8 - 
                0.01754 * xo*xo4*xo8 - 0.01471 * xo6*xo8 - 0.01244 * xo*xo6*xo8 - 0.01062 * xo8*xo8 );
    else if (mhp < 520.)
        n82fr = -( -0.6110 + 1.095 * xm + 0.6492 * xm2 + 0.4596 * xm3 + 0.3569 * xm4 +
                0.2910 * xm3*xm2 + 0.2438 * xm6 + 0.2075 * xm3*xm4 + 0.1785 * xm8 +
                0.1546 * xm3*xm6 + 0.1347 * xm2*xm8 + 0.1177 * xm3*xm8 + 0.1032 * xm4*xm8 +
                0.09073 * xm*xm4*xm8 + 0.07987 * xm6*xm8 + 0.07040 * xm*xm6*xm8 + 0.06210 * xm8*xm8 );
    else n82fr = -( -15961 * x5 + 1003 * x5*Lx3 - 2627 * x5*Lx2 - 29962 * x5*Lx - 
                11683 * x4 + 54.66 * x4*Lx3 - 2777 * x4*Lx2 - 17770 * x4*Lx - 
                6481 * x3 - 40.68 * x3*Lx3 - 1439 * x3*Lx2 - 7906 * x3*Lx - 
                2943 * x2 - 31.83 * x2*Lx3 - 612.6 * x2*Lx2 - 2770 * x2*Lx -
                929.8 * x - 19.80 * x*Lx3 - 174.7 * x*Lx2 - 658.4 * x*Lx );
        
    
    switch (order){
        case NNLO:
            CWbsgArrayNNLO[2] = ( n32ct + cd32ct * lstmu ) * abssu2;
            CWbsgArrayNNLO[3] = ( n42ct + cd42ct * lstmu ) * abssu2;
            CWbsgArrayNNLO[4] = 2./15. * n41ct * abssu2 - CWbsgArrayNNLO[2] / 10.;
            CWbsgArrayNNLO[5] = 0.25 * n41ct * abssu2 - CWbsgArrayNNLO[2] * 3./16.;
            CWbsgArrayNNLO[6] = -( n72ct + cd72ct * lstmu + ce72ct * lstmu * lstmu) * abssu2
                                +( n72fr + cd72fr * lstmu + ce72fr * lstmu * lstmu) * susd
                     - 1./3.*CWbsgArrayNNLO[2] - 4./9.*CWbsgArrayNNLO[3] - 20./3.*CWbsgArrayNNLO[4] - 80./9.*CWbsgArrayNNLO[5];
            CWbsgArrayNNLO[7] = -( n82ct + cd82ct * lstmu + ce82ct * lstmu * lstmu) * abssu2
                                +( n82fr + cd82fr * lstmu + ce82fr * lstmu * lstmu) * susd
                    + CWbsgArrayNNLO[2] - 1./6.*CWbsgArrayNNLO[3] - 20.*CWbsgArrayNNLO[4] - 10./3.*CWbsgArrayNNLO[5];
        case NLO:
            CWbsgArrayNLO[3] = n41ct * abssu2;
            CWbsgArrayNLO[6] = ( n71ct + cd71ct * lstmu ) * abssu2 - (n71fr + cd71fr * lstmu) * susd - 4./9.*CWbsgArrayNLO[3];
            CWbsgArrayNLO[7] = ( n81ct + cd81ct * lstmu ) * abssu2 - (n81fr + cd81fr * lstmu) * susd - 1./6.*CWbsgArrayNLO[3];
        case LO:
            CWbsgArrayLO[6] = n70ct * abssu2 - n70fr * susd;
            CWbsgArrayLO[7] = n80ct * abssu2 - n80fr * susd;
            break;
        default:
            std::stringstream out;
            out << order;
            throw std::runtime_error("order" + out.str() + "not implemeted"); 
            }
    
    /*std::cout << "CWbsgArrayLO[6] = " << CWbsgArrayLO[6] << std::endl;
    std::cout << "CWbsgArrayLO[7] = " << CWbsgArrayLO[7] << std::endl;
    std::cout << "CWbsgArrayNLO[3] = " << CWbsgArrayNLO[3] << std::endl;
    std::cout << "CWbsgArrayNLO[6] = " << CWbsgArrayNLO[6] << std::endl;
    std::cout << "CWbsgArrayNLO[7] = " << CWbsgArrayNLO[7] << std::endl;
    std::cout << "CWbsgArrayNNLO[2] = " << CWbsgArrayNNLO[2] << std::endl;
    std::cout << "CWbsgArrayNNLO[3] = " << CWbsgArrayNNLO[3] << std::endl;
    std::cout << "CWbsgArrayNNLO[4] = " << CWbsgArrayNNLO[4] << std::endl;
    std::cout << "CWbsgArrayNNLO[5] = " << CWbsgArrayNNLO[5] << std::endl;
    std::cout << "CWbsgArrayNNLO[6] = " << CWbsgArrayNNLO[6] << std::endl;
    std::cout << "CWbsgArrayNNLO[7] = " << CWbsgArrayNNLO[7] << std::endl;*/
    
    
    switch (order){
        case NNLO:
            return ( CWbsgArrayNNLO[i] );
        case NLO:
            return ( CWbsgArrayNLO[i] );
            break;
        case LO:
            return ( CWbsgArrayLO[i] );
            break;
        default:
            std::stringstream out;
            out << order;
            throw std::runtime_error("order" + out.str() + "not implemeted"); 
        }
}

