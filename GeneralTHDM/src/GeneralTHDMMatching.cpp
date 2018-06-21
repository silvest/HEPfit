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
    mcbtaunu(3, NDR, LO),
    mcbsg(8, NDR, NNLO),
    mcgminus2mu(2, NDR, NLO),
    mcbsmm(8, NDR, NNLO, NLO_QED22),
    mcBMll(13, NDR, NLO)

{
}
void GeneralTHDMMatching::updateGTHDMParameters()
{
    GF=myGTHDM.getGF();
    mMU=myGTHDM.getLeptons(StandardModel::MU).getMass();
}

double GeneralTHDMMatching::gminus2muLO() {
  
    updateGTHDMParameters();
    double gminus2muLO;
    
     //add something to note that this is only valid in the aligned case and in the CP-conserving limit
    
    
    
    double pi=M_PI;
    double GF=myGTHDM.getGF();
    double mMU=myGTHDM.getLeptons(StandardModel::MU).getMass();
    
    double mHp2=myGTHDM.getmHp2();
    double mHl=myGTHDM.getMHl();
    //double mH1_2=mHl*mHl;
    

    gslpp::complex sl = myGTHDM.getNl_11();
 
    
    /*Mass of the physical scalars*/
   
   
    
    double mH1_2 = myGTHDM.getMyGTHDMCache()->mH1sq;
    double mH2_2 = myGTHDM.getMyGTHDMCache()->mH2sq;
    double mH3_2 = myGTHDM.getMyGTHDMCache()->mH3sq;

    /*eta is the deviation from the SM limit*/
    double alpha1 = myGTHDM.getalpha1();
    double beta = atan(myGTHDM.gettanb());
    double eta = M_PI/2.0-(beta-alpha1);
    
    
    
 //in the CP-conserving limit su, sd, sl are real. have to fit it
    
 
    double yl_h = 1.0 + eta*sl.real();
    double yl_H = -sl.real() + eta;
    double yl_A = - sl.real();
    
    double rmu_hSM, rmu_h, rmu_H, rmu_A, rmu_Hp;
    double part_hSM, part_h, part_H, part_A, part_Hp;
    
    if( mH1_2<1.0 || mH2_2<1.0 || mH3_2<1.0 || mHp2<1.0)
    {
        throw std::runtime_error("The implemented approximation for g-2_\\mu only works for Higgs masses above 1 GeV.");
    }
    
    if (!myGTHDM.getATHDMflag())
    {
        throw std::runtime_error("(g-2) is only aviable in the A2HDM.");
    }
    
     if (!myGTHDM.getCPconservationflag())
    {
        throw std::runtime_error("(g-2) is only aviable in the CP-conserving limit.");
    }
     if( mH1_2 == mH2_2 || mH1_2==mH3_2 || mH1_2== mHp2 || mH2_2==mH3_2 || mH2_2== mHp2 || mH3_2== mHp2)
    {
        throw std::runtime_error(" Masses are equal. Not valid point.");
    }
    
    
    
    rmu_hSM=mMU*mMU/mH1_2;
    rmu_h=mMU*mMU/mH1_2;
    rmu_H=mMU*mMU/mH2_2;
    rmu_A=mMU*mMU/mH3_2;
    rmu_Hp=mMU*mMU/mHp2;
    
    part_hSM=rmu_hSM*(-7.0/6.0-log(rmu_hSM));
    part_h=yl_h*yl_h*rmu_h*(-7.0/6.0-log(rmu_h));
    part_H=yl_H*yl_H*rmu_H*(-7.0/6.0-log(rmu_H));
    part_A=yl_A*yl_A*rmu_A*(11.0/6.0+log(rmu_A));
    part_Hp=-yl_A*yl_A*rmu_Hp*1.0/6.0;
    
   
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

    
    gslpp::complex GeneralTHDMMatching::TF(double m1, double m2, double m3){
            
        double pi=M_PI;
        
        double ml, mm, mh;
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



double GeneralTHDMMatching::gminus2muNLOF() {
    
    updateGTHDMParameters();
    double  aFNphoton, aFNZ;
    gslpp::complex aFC, gminus2muNLOF; 
    
    
     //add something to note that this is only valid in the aligned case and in the CP-conserving limit
    
    double pi=M_PI;
    double GF=myGTHDM.getGF();
    double aem=myGTHDM.getAle();
    double mMU=myGTHDM.getLeptons(StandardModel::MU).getMass();
    double mTAU=myGTHDM.getLeptons(StandardModel::TAU).getMass();
    double mt=myGTHDM.getQuarks(QCD::TOP).getMass();
    double mb=myGTHDM.getQuarks(QCD::BOTTOM).getMass();
    double MZ = myGTHDM.getMz();
    double MW = myGTHDM.Mw();
    double MW2 = MW*MW;
    double cW2 = myGTHDM.cW2();
    
    double sW2 = myGTHDM.sW2();

    double SW4 = sW2*sW2;
    
    double mHp2=myGTHDM.getmHp2();
    
    gslpp::complex su = myGTHDM.getNu_11();
    gslpp::complex sd = myGTHDM.getNd_11();
    gslpp::complex sl = myGTHDM.getNl_11();

 
    
    /*Mass of the physical scalars*/
   
    
    double mH1_2 = myGTHDM.getMyGTHDMCache()->mH1sq;
    double mH2_2 = myGTHDM.getMyGTHDMCache()->mH2sq;
    double mH3_2 = myGTHDM.getMyGTHDMCache()->mH3sq;

    
   
    /*eta is the deviation from the SM limit*/
    double alpha1 = myGTHDM.getalpha1();
    double beta = atan(myGTHDM.gettanb());
    double eta = M_PI/2.0-(beta-alpha1);
    
    
    double mt2 = mt*mt;
    double mb2 = mb*mb;
    double mTAU2 = mTAU*mTAU;
    
 //in the CP-conserving limit su, sd, sl are real. have to fit it
    
    double yt_h = 1.0 + eta*su.real();
    double yb_h = 1.0 + eta*sd.real();
    double ytau_h = 1.0 + eta*sl.real();
    
    double yt_H = -su.real() + eta;
    double yb_H = -sd.real() + eta;
    double ytau_H = -sl.real() + eta;
    
    double yt_A = su.real();
    double yb_A = - sd.real();
    double ytau_A = - sl.real();
    
    double rtau_hSM, rtau_h, rtau_H, rtau_A, rt_hSM, rt_h, rt_H, rt_A, rb_hSM, rb_h, rb_H, rb_A;
    double stau_hSM, stau_h, stau_H, stau_A, st_hSM, st_h, st_H, st_A, sb_hSM, sb_h, sb_H, sb_A;

    double part_hSM_photon, part_h_photon, part_H_photon, part_A_photon;
    double part_hSM_Z, part_h_Z, part_H_Z, part_A_Z;

    
    if( mH1_2<1.0 || mH2_2<1.0 || mH3_2<1.0 || mHp2<1.0)
    {
        throw std::runtime_error("The implemented approximation for g-2_\\mu only works for Higgs masses above 1 GeV.");
    }
    

    
    rtau_hSM=mTAU*mTAU/mH1_2;
    rtau_h=mTAU*mTAU/mH1_2;
    rtau_H=mTAU*mTAU/mH2_2;
    rtau_A=mTAU*mTAU/mH3_2;
    rt_hSM=mt*mt/mH1_2;
    rt_h=mt*mt/mH1_2;
    rt_H=mt*mt/mH2_2;
    rt_A=mt*mt/mH3_2;
    rb_hSM=mb*mb/mH1_2;
    rb_h=mb*mb/mH1_2;
    rb_H=mb*mb/mH2_2;
    rb_A=mb*mb/mH3_2;
    
    stau_hSM=mTAU*mTAU/(mH1_2 - MZ*MZ);
    stau_h=mTAU*mTAU/(mH1_2-MZ*MZ);
    stau_H=mTAU*mTAU/(mH2_2-MZ*MZ);
    stau_A=mTAU*mTAU/(mH3_2-MZ*MZ);
    st_hSM=mt*mt/(15647.5081-MZ*MZ);
    st_h=mt*mt/(mH1_2-MZ*MZ);
    st_H=mt*mt/(mH2_2-MZ*MZ);
    st_A=mt*mt/(mH3_2-MZ*MZ);
    sb_hSM=mb*mb/(mH1_2 -MZ*MZ);
    sb_h=mb*mb/(mH1_2-MZ*MZ);
    sb_H=mb*mb/(mH2_2-MZ*MZ);
    sb_A=mb*mb/(mH3_2-MZ*MZ);
    
    part_hSM_photon=rtau_hSM*gscalar2(rtau_hSM)+(4.0/3.0)*rt_hSM*gscalar2(rt_hSM)+(1.0/3.0)*rb_hSM*gscalar2(rb_hSM);
    part_h_photon=ytau_h*(ytau_h*rtau_h*gscalar2(rtau_h)+(4.0/3.0)*yt_h*rt_h*gscalar2(rt_h)+(1.0/3.0)*yb_h*rb_h*gscalar2(rb_h));
    part_H_photon=ytau_H*(ytau_H*rtau_H*gscalar2(rtau_H)+(4.0/3.0)*yt_H*rt_H*gscalar2(rt_H)+(1.0/3.0)*yb_H*rb_H*gscalar2(rb_H));
    part_A_photon=ytau_A*(ytau_A*rtau_A*gpseudoscalar2(rtau_A)+(4.0/3.0)*yt_A*rt_A*gpseudoscalar2(rt_A)+(1.0/3.0)*yb_A*rb_A*gpseudoscalar2(rb_A));

    part_hSM_Z = stau_hSM*(3.0/4.0 - cW2)*(3.0/4.0 - cW2)*(gscalar2(rtau_hSM)-gscalar2(mTAU*mTAU/(MZ*MZ)))/(cW2*(1-cW2)) + (-2.0)*st_hSM*(3.0/4.0 - cW2)*(-5.0/12.0 + 2.0/3.0*cW2)*(gscalar2(rt_hSM)-gscalar2(mt*mt/(MZ*MZ)))/(cW2*(1-cW2))+ (1.0)*sb_hSM*(3.0/4.0 - cW2)*(1.0/12.0 + 1.0/3.0*cW2)*(gscalar2(rb_hSM)-gscalar2(mb*mb/(MZ*MZ)))/(cW2*(1-cW2));
    part_h_Z = ytau_h*(ytau_h*stau_h*(3.0/4.0 - cW2)*(3.0/4.0 - cW2)*(gscalar2(rtau_h)-gscalar2(mTAU*mTAU/(MZ*MZ)))/(cW2*(1-cW2)) + (-2.0)*yt_h*st_h*(3.0/4.0 - cW2)*(-5.0/12.0 + 2.0/3.0*cW2)*(gscalar2(rt_h)-gscalar2(mt*mt/(MZ*MZ)))/(cW2*(1-cW2))+ (1.0)*yb_h*sb_h*(3.0/4.0 - cW2)*(1.0/12.0 + 1.0/3.0*cW2)*(gscalar2(rb_h)-gscalar2(mb*mb/(MZ*MZ)))/(cW2*(1-cW2)));
    part_H_Z = ytau_H*(ytau_H*stau_h*(3.0/4.0 - cW2)*(3.0/4.0 - cW2)*(gscalar2(rtau_H)-gscalar2(mTAU*mTAU/(MZ*MZ)))/(cW2*(1-cW2)) + (-2.0)*yt_H*st_H*(3.0/4.0 - cW2)*(-5.0/12.0 + 2.0/3.0*cW2)*(gscalar2(rt_H)-gscalar2(mt*mt/(MZ*MZ)))/(cW2*(1-cW2))+ (1.0)*yb_H*sb_H*(3.0/4.0 - cW2)*(1.0/12.0 + 1.0/3.0*cW2)*(gscalar2(rb_H)-gscalar2(mb*mb/(MZ*MZ)))/(cW2*(1-cW2)));
    part_A_Z = ytau_A*(ytau_A*stau_A*(3.0/4.0 - cW2)*(3.0/4.0 - cW2)*(gpseudoscalar2(rtau_A)-gpseudoscalar2(mTAU*mTAU/(MZ*MZ)))/(cW2*(1-cW2)) + (-2.0)*yt_A*st_A*(3.0/4.0 - cW2)*(-5.0/12.0 + 2.0/3.0*cW2)*(gpseudoscalar2(rt_A)-gpseudoscalar2(mt*mt/(MZ*MZ)))/(cW2*(1-cW2))+ (1.0)*yb_h*sb_A*(3.0/4.0 - cW2)*(1.0/12.0 + 1.0/3.0*cW2)*(gpseudoscalar2(rb_hSM)-gpseudoscalar2(mb*mb/(MZ*MZ)))/(cW2*(1-cW2)));

    
    aFNphoton = GF*mMU*mMU/(4.0*pi*pi*pi*sqrt(2.0)) * aem * (-part_hSM_photon+part_h_photon+part_H_photon+part_A_photon);
    aFNZ = GF*mMU*mMU/(4.0*pi*pi*pi*sqrt(2.0)) * aem * (-part_hSM_Z+part_h_Z+part_H_Z+part_A_Z);
    
   
  aFC = sl*((aem*aem*mTAU*mTAU*mMU*mMU*
      sl*(mTAU*mTAU/mHp2 - 
        mTAU*mTAU/MW2 + (-1.0/2.0 + mTAU*mTAU/mHp2)*
         log (mTAU*mTAU/mHp2) - (-1.0/2.0 + mTAU*mTAU/MW2)*
         log (mTAU*mTAU/MW2) + (mTAU*
           mTAU*(-1.0 + mTAU*mTAU/mHp2)*(-pi*pi/6.0 + 
             PolyLog.Li2 (1.0 - mHp2/mTAU2)))/
         mHp2 - (mTAU*
           mTAU*(-1.0 + mTAU*mTAU/MW2)*(-pi*pi/6.0 + 
             PolyLog.Li2 (1.0 - MW2/mTAU2)))/MW2))/(32.0*
      MW2*(mHp2 - MW2)*pi*pi*SW4)- (3.0*aem*aem*mMU*mMU*mt2*
      su*(mb2/mHp2 - mt2/mHp2 - mb2/MW2 + 
        mt2/MW2 + (13.0/12.0 + mb2/mHp2)*
         log (mb2/mHp2) + (13.0/12.0 - mt2/mHp2)*log (mt2/mHp2) + 
              (-log (mb2/mHp2)*log (mb2/mHp2) +  log (mt2/mHp2)*log (mt2/mHp2))/
         3.0 - (13.0/12.0 + mb2/MW2)*
         log (mb2/MW2) - (13.0/12.0 - mt2/MW2)*
         log (mt2/MW2) + (log (mb2/MW2)*log (mb2/MW2) - 
           log (mt2/MW2)*log (mt2/MW2))/
         3.0 - (4.0*(-1.0 - mb2/mHp2 + mt2/mHp2)*
           TF (sqrt (mb2/mHp2), sqrt (mt2/mHp2), 
            1.0))/(3.0*(1.0 + (-(mb2/mHp2) + mt2/mHp2)*(-(mb2/mHp2) + 
                mt2/mHp2) - 
             2.0*(mb2/mHp2 + 
                mt2/mHp2))) + ((-((mb2*(5.0/3.0 + mb2/mHp2))/
                mHp2) + (mt2*(-8.0/3.0 + mt2/mHp2))/mHp2)*
           TF (sqrt (mb2/mHp2), sqrt (mt2/mHp2), 
            1.0))/(1.0 + (-(mb2/mHp2) + mt2/mHp2)*(-(mb2/mHp2) + 
              mt2/mHp2) - 
           2.0*(mb2/mHp2 + mt2/mHp2))+ ((5.0*mb2)/(3.0*mHp2) - (8.0*
              mt2)/(3.0*mHp2) + (-(mb2/mHp2) + 
              mt2/mHp2)*(-(mb2/mHp2) + 
              mt2/mHp2))*(-(log (mb2/mt2)*log (mt2/mHp2))/2.0 + 
           PolyLog.Li2 (1.0 - mb2/mt2) - ((-(mb2/mHp2) + mt2/mHp2)*
              TF (sqrt (mb2/mHp2), sqrt (mt2/mHp2), 
               1.0))/(1.0 + (-(mb2/mHp2) + mt2/mHp2)*(-(mb2/mHp2) + 
                 mt2/mHp2) - 
              2.0*(mb2/mHp2 + mt2/mHp2))) + (4.0*(-1.0 - mb2/MW2 + 
             mt2/MW2)*
           TF (sqrt (mb2/MW2), sqrt (mt2/MW2), 
            1.0))/(3.0*(1.0 + (-(mb2/MW2) + mt2/MW2)*(-(mb2/MW2) + 
                mt2/MW2) - 
             2.0*(mb2/MW2 + 
                mt2/MW2))) - ((-((mb2*(5.0/3.0 + mb2/MW2))/
                MW2) + (mt2*(-8.0/3.0 + mt2/MW2))/MW2)*
           TF (sqrt (mb2/MW2), sqrt (mt2/MW2), 
            1.0))/(1.0 + (-(mb2/MW2) + mt2/MW2)*(-(mb2/MW2) + 
              mt2/MW2) - 
           2.0*(mb2/MW2 + mt2/MW2))- ((-(mb2/MW2) + 
              mt2/MW2)*(-(mb2/MW2) + mt2/MW2) + (5.0*mb2)/(3.0*
              MW2) - (8.0*mt2)/(3.0*
              MW2))*(-(log (mb2/mt2)*log (mt2/MW2))/2.0 + 
           PolyLog.Li2 (1.0 - mb2/mt2) - ((-(mb2/MW2) + mt2/MW2)*
              TF (sqrt (mb2/MW2), sqrt (mt2/MW2), 
               1.0))/(1.0 + (-(mb2/MW2) + mt2/MW2)*(-(mb2/MW2) + 
                 mt2/MW2) - 2.0*(mb2/MW2 + mt2/MW2)))))/(32.0*
      MW2*(mHp2 - MW2)*pi*pi*SW4) + (3.0*aem*aem*mb2*mMU*mMU*
      sd*(mb2/mHp2 - mt2/mHp2 - mb2/MW2 + 
        mt2/MW2 + (1.0/12.0 + mb2/mHp2)*
         log (mb2/mHp2) + (1.0/12.0 - mt2/mHp2)*log (mt2/mHp2) - (1.0/12.0 + mb2/MW2)*
         log (mb2/MW2) - (1.0/12.0 - mt2/MW2)*
         log (mt2/
           MW2) + ((-((mb2*(-1.0/3.0 + mb2/mHp2))/
                mHp2) + (mt2*(-2.0/3.0 + mt2/mHp2))/mHp2)*
           TF (sqrt (mb2/mHp2), sqrt (mt2/mHp2), 
            1.0))/(1.0 + (-(mb2/mHp2) + mt2/mHp2)*(-(mb2/mHp2) + 
              mt2/mHp2) - 
           2.0*(mb2/mHp2 + mt2/mHp2)) + (-mb2/(3.0*mHp2) - (2.0*
              mt2)/(3.0*mHp2) + (-(mb2/mHp2) + 
              mt2/mHp2)*(-(mb2/mHp2) + 
              
              mt2/mHp2))*(-(log (mb2/mt2)*log (mt2/mHp2))/2.0 + 
           PolyLog.Li2 (1.0 - mb2/mt2) - ((-(mb2/mHp2) + mt2/mHp2)*
              TF (sqrt (mb2/mHp2), sqrt (mt2/mHp2), 
               1.0))/(1.0 + (-(mb2/mHp2) + mt2/mHp2)*(-(mb2/mHp2) + 
                 mt2/mHp2) - 
              2.0*(mb2/mHp2 + 
                 mt2/mHp2))) - ((-((mb2*(-1.0/3.0 + mb2/MW2))/
                MW2) + (mt2*(-2.0/3.0 + mt2/MW2))/MW2)*
           TF (sqrt (mb2/MW2), sqrt (mt2/MW2), 
            1.0))/(1.0 + (-(mb2/MW2) + mt2/MW2)*(-(mb2/MW2) + 
              mt2/MW2) - 
           2.0*(mb2/MW2 + mt2/MW2)) - ((-(mb2/MW2) + 
              mt2/MW2)*(-(mb2/MW2) + mt2/MW2) - 
           mb2/(3.0*MW2) - (2.0*mt2)/(3.0*
              MW2))*(-(log (mb2/mt2)*log (mt2/MW2))/2.0 + 
           PolyLog.Li2 (1.0 - mb2/mt2) - ((-(mb2/MW2) + mt2/MW2)*
              TF (sqrt (mb2/MW2), sqrt (mt2/MW2), 
               1.0))/(1.0 + (-(mb2/MW2) + mt2/MW2)*(-(mb2/MW2) + 
                 mt2/MW2) - 2.0*(mb2/MW2 + mt2/MW2)))))/(32.0*
      MW2*(mHp2 - MW2)*pi*pi*SW4)) ;
   
   gminus2muNLOF = aFNphoton + aFNZ +aFC;
        
    return(gminus2muNLOF.real());
}




double GeneralTHDMMatching::gminus2muNLOB() {
    
    updateGTHDMParameters();
    gslpp::complex aEWadd, aNonYuk, aYuk, gminus2muNLOB;
    
    double pi=M_PI;
    double aem=myGTHDM.getAle();
    double mMU=myGTHDM.getLeptons(StandardModel::MU).getMass();
    double MZ = myGTHDM.getMz();
    double CW2 = myGTHDM.cW2();
    
    double MH=myGTHDM.getMHl();
    double mHp2=myGTHDM.getmHp2();
    double mHp=sqrt(mHp2);

   
    double vev = myGTHDM.v();
    
    double M2 = myGTHDM.getMyGTHDMCache()->M2;
    double tanb = myGTHDM.gettanb();

 
    
    //double mH1_2 = myGTHDM.getMyGTHDMCache()->mH1sq*myGTHDM.getMyGTHDMCache()->mH1sq;
    double mH2_2 = myGTHDM.getMyGTHDMCache()->mH2sq;
    double mH3_2 = myGTHDM.getMyGTHDMCache()->mH3sq;

    //Lambda5 defined as in 1607.06292, Eq. (14)
    
    
    double Lambda5 = (2.0*M2)/(vev*vev);
   
    double sl = (myGTHDM.getNl_11()).real();
    
     //This has to be solved to consider the CP-conserving limtit
   
    /*eta is the deviation from the SM limit*/
    double alpha1 = myGTHDM.getalpha1();
    double beta = atan(myGTHDM.gettanb());
    double eta = M_PI/2.0-(beta-alpha1);
    
    double CW = sqrt(CW2);
    double CW4 = CW2*CW2;
    double CW6 = CW4*CW2;
    double CW8 = CW4*CW4;
    double CW10 = CW8*CW2;
    double CW12 = CW10*CW2;
    double CW14 = CW10*CW4;
    
    double MZ2 = MZ*MZ;
    double MZ4 = MZ2*MZ2;
    double MZ6 = MZ4*MZ2;
    double MZ8 = MZ4*MZ4;
    double MZ10 = MZ8*MZ2;
    double MZ12 = MZ8*MZ4;
    double MZ14 = MZ10*MZ4;
    double MZ16 = MZ8*MZ8;
    
    double MH2 = MH*MH;
    double MH4 = MH2*MH2;
    double MH6 = MH4*MH2;
    double MH8 = MH4*MH4;
    double MH10 = MH8*MH2;
    double MH12 = MH10*MH2;
    double MH14 = MH12*MH2;
    double MH16 = MH10*MH4*MH2;
    
    double mH3 = sqrt(mH3_2);
    double mH3_6 = mH3_2*mH3_2*mH3_2;
    double mH3_8 = mH3_6*mH3_2;
    double mH3_10 = mH3_8*mH3_2;
    
    double mHp4 = mHp2*mHp2;
    double mHp3 = mHp2*sqrt(mHp2);
    double mHp6 = mHp4*mHp2;
    double mHp8  = mHp6*mHp2;
    double mHp10  = mHp8*mHp2;
    
    
    double mH2 = sqrt(mH2_2);
    double mH2_4 = mH2_2*mH2_2;
    double mH2_6 = mH2_4*mH2_2;
    double mH2_8 = mH2_4*mH2_4;
    double mH2_10 = mH2_8*mH2_2;
    double mH2_12 = mH2_10*mH2_2;


    
    

    
 aEWadd = (eta*aem*aem*mMU*mMU*(-24.0*(1.0 - 4.0*CW2)*(1.0 - 4.0*CW2)*(1.0- 3.0*CW2 + 2.0*CW4)*MH16*
         pi*pi - 24.0*(1.0- 4.0*CW2)*(1.0- 4.0*CW2)*(1.0- 3.0*CW2 + 2.0*CW4)*MH14*MZ2*
         (-5.0 - 8.0*CW2 + negsquareroot((MH2*(MH2 - 4.0*MZ2))/MZ4))*pi*pi + 
        256.0*CW4*MZ16*(5.0 - 12.0*CW2 + 8.0*CW4 + 
          9.0*CW8*negsquareroot((MH2*(MH2 - 4.0*CW2*MZ2))/(CW4*MZ4)) - 
          CW6*(-2.0 + 9.0*negsquareroot((MH2*(MH2 - 4.0*CW2*MZ2))/(CW4*MZ4))))*
         (pi - 4.0*CW2*pi)*(pi - 4.0*CW2*pi) - 24.0*(1.0 - 4.0*CW2)*(1.0 - 4.0*CW2)*(1.0 - 3.0*CW2 + 2*CW4)*MH12*
         MZ4*(6 + (4.0 + 16.0*CW4 - 3.0*negsquareroot((MH2*(MH2 - 4.0*MZ2))/MZ4) - 
            8.0*CW2*(-5.0 + negsquareroot((MH2*(MH2 - 4.0*MZ2))/MZ4)))*pi*pi) - 
        MH10*MZ6*(-24.0 + 68.0*negsquareroot((MH2*(MH2 - 4.0*MZ2))/MZ4)*pi*pi + 
          12288.0*CW12*(-5.0 + negsquareroot((MH2*(MH2 - 4.0*MZ2))/MZ4))*pi*pi + 
          16.0*CW8*(4512.0 + (-1968.0 + 15.0*negsquareroot(1.0 - 4.0*CW2) - 1168.0*negsquareroot(
                (MH2*(MH2 - 4.0*MZ2))/MZ4))*pi*pi) - 
          6144.0*CW10*(6.0 + (-16.0 + negsquareroot((MH2*(MH2 - 4.0*MZ2))/MZ4))*pi*pi) - 
          4.0*CW6*(10668.0 + (2784.0 + 75.0*negsquareroot(1.0 - 4.0*CW2) - 3968.0*negsquareroot(
                (MH2*(MH2 - 4.0*MZ2))/MZ4) - 128.0*negsquareroot(
                (MH2*(MH2 - 4.0*CW2*MZ2))/(CW4*MZ4)))*pi*pi) + 
          2.0*CW4*(4860.0 + (3264.0 + 57.0*negsquareroot(1.0 - 4.0*CW2)- 1600.0*negsquareroot(
                (MH2*(MH2 - 4.0*MZ2))/MZ4) - 128.0*negsquareroot(
                (MH2*(MH2 - 4.0*CW2*MZ2))/(CW4*MZ4)))*pi*pi) - 
          CW2*(651.0 + 2.0*(384.0 + 9.0*negsquareroot(1.0 - 4.0*CW2) + 80.0*negsquareroot(
                (MH2*(MH2 - 4.0*MZ2))/MZ4) - 16.0*negsquareroot(
                (MH2*(MH2 - 4.0*CW2*MZ2))/(CW4*MZ4)))*pi*pi)) - 
        2*MH4*MH2*MZ10*(-128.0*CW12*(1824 + (57.0*negsquareroot(1.0 - 4.0*CW2) - 
              128.0*negsquareroot((MH2*(MH2 - 4.0*MZ2))/MZ4))*pi*pi) + 
          10.0*(24.0 + (1.0 - 2*negsquareroot((MH2*(MH2 - 4.0*MZ2))/MZ4))*pi*pi) + 
          8.0*CW2*(-12.0 + (-43.0+ 16.0*negsquareroot((MH2*(MH2 - 4.0*MZ2))/MZ4))*pi*pi) + 
          2*CW6*(36516.0 + (-4664.0 + 2460.0*negsquareroot(1.0 - 4.0*CW2) - 3904.0*negsquareroot(
                (MH2*(MH2 - 4.0*MZ2))/MZ4) - 123.0*negsquareroot(
                (MH2*(MH2 - 4.0*CW2*MZ2))/(CW4*MZ4)))*pi*pi) - 
          CW4*(13836 + (-2852 + 504*negsquareroot(1.0 - 4.0*CW2) - 640.0*negsquareroot(
                (MH2*(MH2 - 4.0*MZ2))/MZ4) - 45.0*negsquareroot(
                (MH2*(MH2 - 4.0*CW2*MZ2))/(CW4*MZ4)))*pi*pi) - 
          32.0*CW8*(6792.0 + (-406.0 + 519.0*negsquareroot(1.0 - 4.0*CW2) - 800.0*negsquareroot(
                (MH2*(MH2 - 4.0*MZ2))/MZ4) + 6.0*negsquareroot((MH2*(MH2 - 
                   4.0*CW2*MZ2))/(CW4*MZ4)))*pi*pi) + 
          32.0*CW10*(13236 + (-232.0 + 627.0*negsquareroot(1.0 - 4.0*CW2) - 1088.0*negsquareroot(
                (MH2*(MH2 - 4.0*MZ2))/MZ4) + 57.0*negsquareroot(
                (MH2*(MH2 - 4.0*CW2*MZ2))/(CW4*MZ4)))*pi*pi)) + 
        MH8*MZ8*(-20.0*(-30.0 + (3.0 + negsquareroot((MH2*(MH2 - 4.0*MZ2))/MZ4))*pi*pi) + 
          12288.0*CW12*(-6.0 + (-4.0 + 3.0*negsquareroot((MH2*(MH2 - 4.0*MZ2))/MZ4))*
             pi*pi) - 128.0*CW10*(-624 + (-768.0 + 3.0*negsquareroot(1.0 - 4.0*CW2) + 
              448.0*negsquareroot((MH2*(MH2 - 4.0*MZ2))/MZ4))*pi*pi) + 
          2*CW4*(13056.0 + (-1884.0 + 561.0*negsquareroot(1.0 - 4.0*CW2) - 2736.0*negsquareroot(
                (MH2*(MH2 - 4.0*MZ2))/MZ4) - 29.0*negsquareroot(
                (MH2*(MH2 - 4.0*CW2*MZ2))/(CW4*MZ4)))*pi*pi) - 
          CW2*(6267.0 + 2.0*(-312.0 + 63.0*negsquareroot(1.0 - 4.0*CW2) - 376.0*negsquareroot(
                (MH2*(MH2 - 4.0*MZ2))/MZ4) - 16.0*negsquareroot(
                (MH2*(MH2 - 4.0*CW2*MZ2))/(CW4*MZ4)))*pi*pi) + 
          16.0*CW8*(3144.0 + (-4152.0 + 201.0*negsquareroot(1.0 - 4.0*CW2) + 944.0*negsquareroot(
                (MH2*(MH2 - 4.0*MZ2))/MZ4) + 198.0*negsquareroot(
                (MH2*(MH2 - 4.0*CW2*MZ2))/(CW4*MZ4)))*pi*pi) - 
          4.0*CW6*(14748.0 + (-5040.0 + 867.0*negsquareroot(1.0 - 4.0*CW2) - 2592.0*negsquareroot(
                (MH2*(MH2 - 4.0*MZ2))/MZ4) + 268.0*negsquareroot(
                (MH2*(MH2 - 4.0*CW2*MZ2))/(CW4*MZ4)))*pi*pi)) - 
        4.0*(negsquareroot(CW2) - 4.0*CW2*negsquareroot(CW2))*(negsquareroot(CW2) - 4.0*CW2*negsquareroot(CW2))*MH2*MZ14*(160.0*pi*pi + 
          64.0*CW4*(-72.0+ (1.0 + 6.0*negsquareroot((MH2*(MH2 - 4.0*MZ2))/MZ4))*pi*pi) - 
          16.0*CW2*(-120.0 + (19.0 + 10.0*negsquareroot((MH2*(MH2 - 4.0*MZ2))/MZ4))*pi*pi) + 
          128.0*CW8*(-108.0 + (1.0 + 6.0*negsquareroot((MH2*(MH2 - 4.0*CW2*MZ2))/
                 (CW4*MZ4)))*pi*pi) - CW6*(-17664.0 + 
            (-96.0 + 256.0*negsquareroot((MH2*(MH2 - 4.0*MZ2))/MZ4) + 
              907.0*negsquareroot((MH2*(MH2 - 4.0*CW2*MZ2))/(CW4*MZ4)))*pi*pi)) - 
        4.0*MH4*MZ12*(-20.0*pi*pi + 9216.0*CW14*(24.0 + negsquareroot(1.0 - 4.0*CW2)*pi*pi) + 
          8.0*CW2*(-120 + (21.0 + 10.0*negsquareroot((MH2*(MH2 - 4.0*MZ2))/MZ4))*pi*pi) - 
          16.0*CW4*(-474.0 + (5.0 + 47.0*negsquareroot((MH2*(MH2 - 4.0*MZ2))/MZ4))*pi*pi) + 
          3.0*CW6*(-6108.0 + (-968.0 + 168.0*negsquareroot(1.0 - 4.0*CW2) + 704.0*negsquareroot(
                (MH2*(MH2 - 4.0*MZ2))/MZ4) + 77.0*negsquareroot(
                (MH2*(MH2 - 4.0*CW2*MZ2))/(CW4*MZ4)))*pi*pi) - 
          16.0*CW12*(18432.0 + (-416.0 + 1404.0*negsquareroot(1.0 - 4.0*CW2) - 
              128.0*negsquareroot((MH2*(MH2 - 4.0*MZ2))/MZ4) + 113.0*negsquareroot(
                (MH2*(MH2 - 4.0*CW2*MZ2))/(CW4*MZ4)))*pi*pi) + 
          8.0*CW10*(5784.0 + (-1616.0 + 2190.0*negsquareroot(1.0 - 4.0*CW2) - 
              256.0*negsquareroot((MH2*(MH2 - 4.0*MZ2))/MZ4) + 575.0*negsquareroot(
                (MH2*(MH2 - 4.0*CW2*MZ2))/(CW4*MZ4)))*pi*pi) - 
          CW8*(-16800.0 + (-10080.0 + 5064.0*negsquareroot(1.0 - 4.0*CW2) + 1152.0*negsquareroot(
                (MH2*(MH2 - 4.0*MZ2))/MZ4) + 1961.0*negsquareroot(
                (MH2*(MH2 - 4.0*CW2*MZ2))/(CW4*MZ4)))*pi*pi)))*sl)/
      (4608.0*CW4*(1.0 - 4.0*CW2)*(1.0 - 4.0*CW2)*(-1.0 + CW2)*(-1.0 + CW2)*MH4*MZ8*(MH2 - MZ2)*
       (MH2 - 4.0*CW2*MZ2)*(MH2 - 4.0*CW2*MZ2)*pi*pi) + 
     (eta*aem*aem*mMU*mMU*((71.0 - 278.0*CW2 + 120.0*CW4)*MH4 + 
        6.0*CW2*(9.0 - 88.0*CW2 + 112.0*CW4)*MH2*MZ2 - 
        128.0*CW4*(10.0 - 49.0*CW2 + 36.0*CW4)*MZ4)*sl*log(CW2))/
      (768.0*CW2*(-1.0 + CW2)*(-1.0 + CW2)*(-1.0 + 4.0*CW2)*MH2*MZ2*(-MH2 + 4.0*CW2*MZ2)*
       pi*pi) - (eta*aem*aem*(-3.0 + 4.0*CW2)*mMU*mMU*((1.0 - 5.0*CW2 + 10.0*CW4)*MH2 + 
        (-7.0 + 61.0*CW2 - 162.0*CW4 + 96.0*CW6)*MZ2)*sl*log(CW2)*log(CW2))/
      (256.0*CW2*negpow((1.0 - 4.0*CW2), 3.0/2.0)*(-1.0 + CW2)*(-1.0 + CW2)*MZ2*(-MH2 + MZ2)*pi*pi) + 
     (eta*aem*aem*(-3.0 + 4.0*CW2)*mMU*mMU*((1.0 - 5.0*CW2 + 10.0*CW4)*MH2 + 
        (-7.0 + 61.0*CW2 - 162.0*CW4 + 96.0*CW6)*MZ2)*sl*log((1.0 - negsquareroot(1.0 - 4.0*CW2))/2.0)*log((1.0 - negsquareroot(1.0 - 4.0*CW2))/2.0))/(128.0*CW2*negpow((1.0 - 4.0*CW2),3.0/2.0)*(-1.0 + CW2)*(-1.0 + CW2)*MZ2*(-MH2 + MZ2)*pi*pi)
    -  (eta*aem*aem*mMU*mMU*negsquareroot((MH4 - 4.0*MH2*MZ2)/MZ4)*
       (6.0*(1.0 - 3.0*CW2 + 2.0*CW4)*MH4*MH2 - 12.0*(1.0 - 3.0*CW2 + 2.0*CW4)*MH4*MZ2 + 
        (5.0 - 12.0*CW2 + 8.0*CW4)*MH2*MZ4 + 2.0*(5.0 - 12.0*CW2 + 8.0*CW4)*MZ6)*
       sl*log((2.0 - negsquareroot(MH4/MZ4 - (4.0*MH2)/MZ2) - MH2/MZ2)/2.0)*
      log((-negsquareroot(MH4/MZ4 - (4.0*MH2)/MZ2) + MH2/MZ2)/2.0))/
      (192.0*CW4*(-1.0 + CW2)*(-1.0 + CW2)*MH2*MZ6*pi*pi) + 
     (eta*aem*aem*mMU*mMU*negsquareroot((MH4 - 4.0*CW2*MH2*MZ2)/(CW4*MZ4))*
       (-16.0*MH10 + (16.0 + 99.0*CW2)*MH8*MZ2 - 3.0*CW2*(15.0 + 38.0*CW2)*MH4*MH2*
         MZ4 + 2.0*CW4*(-231.0 + 113.0*CW2)*MH4*MZ6 - 2.0*CW6*(-907.0 + 768.0*CW2)*
         MH2*MZ8 + 1152.0*CW8*(-1.0 + CW2)*MZ10)*sl*
       log((2.0 - negsquareroot(MH4/(CW4*MZ4) - (4.0*MH2)/(CW2*MZ2)) - 
          MH2/(CW2*MZ2))/2.0)*
       log((-negsquareroot(MH4/(CW4*MZ4) - (4.0*MH2)/(CW2*MZ2)) + 
          MH2/(CW2*MZ2))/2.0))/(384.0*CW2*(-1.0 + CW2)*(-1.0 + CW2)*MH4*MZ2*
       (MH2 - MZ2)*(MH2 - 4.0*CW2*MZ2)*(MH2 - 4.0*CW2*MZ2)*pi*pi) + 
     (eta*aem*aem*mMU*mMU*(-48.0*(1.0 - 3.0*CW2 + 2*CW4)*MH8 + 
        (10.0 + 73.0*CW2 - 512.0*CW4 + 384.0*CW6)*MH4*MH2*MZ2 + 
        (-60.0 + 151.0*CW2 + 258.0*CW4 - 256.0*CW6)*MH4*MZ4 + 
        2*(40 + 24.0*CW2 - 149.0*CW4 + 256.0*CW6)*MH2*MZ6 + 
        64*CW2*(-5.0 + 12.0*CW2 - 28.0*CW4 + 18.0*CW6)*MZ8)*sl*
       log(MH2/MZ2))/(768.0*CW4*(-1.0 + CW2)*(-1.0 + CW2)*MH2*MZ4*(MH2 - MZ2)*
       (MH2 - 4.0*CW2*MZ2)*pi*pi) - (eta*aem*aem*(-1.0 + 2.0*CW2)*MH4*mMU*mMU*
       (MH2 - 4.0*MZ2)*sl*log(MH2/MZ2)*log(MH2/MZ2))/(64*CW4*(-1.0 + CW2)*MZ8*
       pi*pi) - (eta*aem*aem*(-3.0 + 4.0*CW2)*mMU*mMU*((1.0 - 5.0*CW2 + 10.0*CW4)*MH2 + 
        (-7.0 + 61.0*CW2 - 162.0*CW4 + 96.0*CW6)*MZ2)*sl*
       PolyLog.Li2( (1.0 - negsquareroot(1.0 - 4.0*CW2))/2))/(64.0*CW2*negpow((1.0 - 4.0*CW2), 3.0/2.0)*
       (-1.0 + CW2)*
       (-1.0 + CW2)*MZ2*(-MH2 + MZ2)*pi*pi) - 
     (eta*aem*aem*mMU*mMU*(12.0*(1.0 - 3.0*CW2 + 2.0*CW4)*MH10 - 
        48.0*(1.0 - 3.0*CW2 + 2*CW4)*MH8*MZ2 + (17.0 - 48.0*CW2 + 32.0*CW4)*MH4*MH2*
         MZ4 - 3.0*(5.0 - 12.0*CW2 + 8.0*CW4)*MH2*MZ8 - 4.0*(5.0 - 12.0*CW2 + 8.0*CW4)*
         MZ10)*sl*PolyLog.Li2(1.0 - MH2/MZ2))/(192.0*CW4*(-1.0 + CW2)*(-1.0 + CW2)*MH4*
       MZ8*pi*pi) + (eta*aem*aem*mMU*mMU*negsquareroot((MH4 - 4.0*MH2*MZ2)/MZ4)*
       (6.0*(1.0 - 3.0*CW2 + 2.0*CW4)*MH4*MH2 - 12.0*(1.0 - 3.0*CW2 + 2.0*CW4)*MH4*MZ2 + 
        (5.0 - 12.0*CW2 + 8.0*CW4)*MH2*MZ4 + 2.0*(5.0 - 12.0*CW2 + 8.0*CW4)*MZ6)*
       sl*PolyLog.Li2((2.0 - negsquareroot(MH4/MZ4 - (4.0*MH2)/MZ2) - MH2/MZ2)/2.0))/
      (192.0*CW4*(-1.0 + CW2)*(-1.0 + CW2)*MH2*MZ6*pi*pi) + 
     (eta*aem*aem*mMU*mMU*negsquareroot((MH4 - 4.0*MH2*MZ2)/MZ4)*
       (6.0*(1.0 - 3.0*CW2 + 2*CW4)*MH4*MH2 - 12.0*(1.0 - 3.0*CW2 + 2.0*CW4)*MH4*MZ2 + 
        (5.0 - 12.0*CW2 + 8.0*CW4)*MH2*MZ4 + 2.0*(5.0 - 12.0*CW2 + 8.0*CW4)*MZ6)*
       sl*PolyLog.Li2( (-negsquareroot(MH4/MZ4 - (4.0*MH2)/MZ2) + MH2/MZ2)/2.0))/
      (192.0*CW4*(-1.0 + CW2)*(-1.0 + CW2)*MH2*MZ6*pi*pi) + 
     (eta*aem*aem*mMU*mMU*(-16.0*MH4*MH2 + 3.0*CW2*MH4*MZ2 + 12.0*CW4*MH2*MZ4 + 
        16.0*CW6*MZ6)*sl*PolyLog.Li2( 1.0 - MH2/(CW2*MZ2)))/
      (384.0*CW4*(-1.0 + CW2)*(-1.0 + CW2)*MH4*MZ4*pi*pi) + 
     (eta*aem*aem*mMU*mMU*negsquareroot((MH4 - 4.0*CW2*MH2*MZ2)/(CW4*MZ4))*
       (16.0*MH10 - (16.0 + 99.0*CW2)*MH8*MZ2 + 3.0*CW2*(15.0 + 38.0*CW2)*MH4*MH2*
         MZ4 + 2*CW4*(231.0 - 113.0*CW2)*MH4*MZ6 + 2.0*CW6*(-907.0 + 768.0*CW2)*
         MH2*MZ8 - 1152.0*CW8*(-1.0 + CW2)*MZ10)*sl*
      PolyLog.Li2( (2.0 - negsquareroot(MH4/(CW4*MZ4) - (4.0*MH2)/(CW2*MZ2)) - 
          MH2/(CW2*MZ2))/2.0))/(384.0*CW2*(-1.0 + CW2)*(-1.0 + CW2)*MH4*MZ2*
       (MH2 - MZ2)*(MH2 - 4.0*CW2*MZ2)*(MH2 - 4.0*CW2*MZ2)*pi*pi) + 
     (eta*aem*aem*mMU*mMU*negsquareroot((MH4 - 4.0*CW2*MH2*MZ2)/(CW4*MZ4))*
       (16.0*MH10 - (16.0 + 99.0*CW2)*MH8*MZ2 + 3.0*CW2*(15.0 + 38.0*CW2)*MH4*MH2*
         MZ4 + 2.0*CW4*(231.0 - 113.0*CW2)*MH4*MZ6 + 2.0*CW6*(-907.0 + 768.0*CW2)*
         MH2*MZ8 - 1152.0*CW8*(-1.0+ CW2)*MZ10)*sl*
       PolyLog.Li2( (-negsquareroot(MH4/(CW4*MZ4) - (4.0*MH2)/(CW2*MZ2)) + 
          MH2/(CW2*MZ2))/2.0))/(384.0*CW2*(-1.0+ CW2)*(-1.0+ CW2)*MH4*MZ2*
       (MH2 - MZ2)*(MH2 - 4.0*CW2*MZ2)*(MH2 - 4.0*CW2*MZ2)*pi*pi); 
 

 
aNonYuk = (aem*aem*mMU*
     mMU*(-15.0*(-25.0 + 32.0*CW2 + CW4 - 16.0*CW6 + 8.0*CW8)*mH3_2*
        mH3_2*mHp2 - 
       15.0*(-25.0 + 32.0*CW2 + CW4 - 16.0*CW6 + 8.0*CW8)*mH2_4*
        mHp2 - 15.0*mH2_2*
        mHp2*((50.0 - 64.0*CW2 + 8.0*CW4)*mHp2 + 
          CW2*(13.0 - 43.0*CW2 + 50.0*CW4 - 20.0*CW6)*MZ2) + 
       2.0*(15.0*(25.0 - 32.0*CW2 + 4.0*CW4)*mHp6 + 
          30.0*CW2*(-17.0 + 41.0*CW2 - 80.0*CW4 + 184.0*CW6 - 
             192.0*CW8 + 64.0*CW10)*mHp4*MZ2 - 
          80.0*CW4*(-1.0 + CW2)*(-1.0 + CW2)*(-1.0 - 8.0*CW2 + 
             8.0*CW4)*mHp2*MZ4 + 
          24.0*CW6*(-1.0 + CW2)*(-1.0 + CW2)*(-1.0 + CW2)*MZ6) + 
       15.0*mH3_2*
        mHp2*(-50.0*mHp2 + CW2*(64.0*mHp2 - 13.0*MZ2) + 
          4.0*CW8*(4.0*mH2_2 + 5.0*MZ2) - 
          2.0*CW6*(16.0*mH2_2 + 25.0*MZ2) + 
          CW4*(10.0*mH2_2 - 8.0*mHp2 + 43.0*MZ2))))/(17280.0*
     CW6*(-1.0 + CW2)*(-1.0 + CW2)*(-1.0 + CW2)*mHp2*MZ6*pi*pi) + 
  log(CW2)*((aem*aem*mMU*
        mMU*(2.0*mH3_2*mH3_2 + 2.0*mH2_4 + 4.0*mHp4 - 
          2.0*CW2*mHp2*MZ2 + mH3_2*(-4.0*mHp2 + CW2*MZ2) + 
          mH2_2*(-4.0*mHp2 + CW2*MZ2)))/(128.0*
        CW6*(-1.0 + CW2)*(-1.0 + CW2)*MZ6*pi*pi) - (aem*
        aem*(mH3_2 - mHp2)*mMU*
        mMU*(mH3_2*mH3_2 - 2.0*mH3_2*mHp2 + mHp4 - CW2*mHp2*MZ2)*
        log(mH3_2/mHp2))/(128.0*CW8*(-1.0 + CW2)*(-1.0 + CW2)*MZ8*
        pi*pi) - (aem*aem*(mH2_2 - mHp2)*mMU*
        mMU*(mH2_4 - 2.0*mH2_2*mHp2 + mHp4 - CW2*mHp2*MZ2)*
        log(mH2_2/mHp2))/(128.0*CW8*(-1.0 + CW2)*(-1.0 + CW2)*MZ8*
        pi*pi)) + (aem*aem*mMU*
     mMU*((-7.0 + 14.0*CW2 - 4.0*CW4 + 5.0*CW6 - 16.0*CW8 + 8.0*CW10)*
        mH3_10 + 
       mH3_8*((7.0 - 14.0*CW2 + 4.0*CW4 - 20.0*CW6 + 64.0*CW8 - 
             32.0*CW10)*
           mH2_2 + (28.0 - 56.0*CW2 + 16.0*CW4 - 5.0*CW6 + 
             16.0*CW8 - 8.0*CW10)*mHp2 - 
          3.0*CW2*(-19.0 + 26.0*CW2 + CW4 - 16.0*CW6 + 8.0*CW8)*MZ2) +
        mH3_6*(-14.0*mHp2*(2.0*mH2_2 + 3.0*mHp2) + 
          CW2*(84.0*mHp4 - 93.0*mHp2*MZ2 + 
             mH2_2*(56.0*mHp2 - 57*MZ2)) + 
          8.0*CW10*(6.0*mH2_4 + 3.0*MZ2*(mHp2 + MZ2) + 
             mH2_2*(4.0*mHp2 + 3.0*MZ2)) - 
          4.0*CW8*(24.0*mH2_4 + 4.0*mH2_2*(4.0*mHp2 + 3.0*MZ2) + 
             3.0*MZ2*(4.0*mHp2 + 5.0*MZ2)) + 
    CW6*(30.0*mH2_4 + mH2_2*(20.0*mHp2 + 3.0*MZ2) + 
             3.0*MZ2*(mHp2 + 13.0*MZ2)) + 
          CW4*(mH2_2*(-16.0*mHp2 + 78.0*MZ2) - 
             3.0*(8.0*mHp4 - 38*mHp2*MZ2 + MZ4))) + 
       mH3_2*mH3_2*(42*mH2_2*mHp4 + 28.0*mHp6 - 
          4.0*CW10*(8.0*mH2_6 + 6.0*mH2_2*mHp2*MZ2 + 
             6.0*mHp2*MZ4 - MZ6 + 6.0*mH2_4*(2.0*mHp2 - MZ2)) - 
          CW6*(20.0*mH2_6 + 12*mHp4*MZ2 + 69.0*mHp2*MZ4 + 
             2.0*MZ6 + 3.0*mH2_2*MZ2*(mHp2 - 10.0*MZ2) + 
             15.0*mH2_4*(2.0*mHp2 - MZ2)) + 
          CW2*(-56.0*mHp6 + 15.0*mHp4*MZ2 + 
             mH2_2*(-84.0*mHp4 + 93.0*mHp2*MZ2)) + 
          2.0*CW8*(32.0*mH2_6 + 36*mHp2*MZ4 - MZ6 + 
             24.0*mH2_4*(2.0*mHp2 - MZ2) + 
             6.0*mH2_2*(4.0*mHp2*MZ2 - MZ4)) + 
          CW4*(3.0*mH2_2*(8.0*mHp4 - 38*mHp2*MZ2 - 3.0*MZ4) + 
             2.0*mHp2*(8.0*mHp4 + 3.0*mHp2*MZ2 + 6.0*MZ4))) + 
       mH3_2*(-7.0*(4.0*mH2_2*mHp6 + mHp8) - 
          2.0*CW8*(8.0*mH2_8 + 6.0*mHp4*MZ4 - mH2_2*MZ6 - 
             mHp2*MZ6 + 8*mH2_6*(4.0*mHp2 - 3.0*MZ2) + 
             24.0*mH2_4*MZ2*(-mHp2 + MZ2)) - 
          CW4*mHp4*(4.0*mHp4 + 42*mHp2*MZ2 + 21.0*MZ4 + 
             2.0*mH2_2*(8.0*mHp2 + 3.0*MZ2)) + 
          CW2*(7.0*mHp6*(2.0*mHp2 + 3.0*MZ2) + 
             mH2_2*(56.0*mHp6 - 15.0*mHp4*MZ2)) + 
          4.0*CW10*(2.0*mH2_8 - mH2_2*MZ6 - mHp2*MZ6 + 
             mH2_6*(8.0*mHp2 - 6.0*MZ2) + 
             mH2_4*(-6.0*mHp2*MZ2 + 6.0*MZ4)) + 
          CW6*(5.0*mH2_8 + 5.0*mH2_6*(4.0*mHp2 - 3.0*MZ2) + 
             15.0*mH2_4*MZ2*(-mHp2 + MZ2) + 
             2.0*mHp2*MZ2*(6.0*mHp4 + 21.0*mHp2*MZ2 + MZ4) + 
             2.0*mH2_2*(6.0*mHp4*MZ2 + MZ6))) + 
       mH2_2*mHp2*(7.0*mHp6 - 7.0*CW2*(2.0*mHp6 + 3.0*mHp4*MZ2) + 
          CW4*(4.0*mHp6 + 42*mHp4*MZ2 + 21.0*mHp2*MZ4) - 
          4.0*CW10*(2.0*mH2_6 - 6.0*mH2_4*MZ2 + 6.0*mH2_2*MZ4 - 
             MZ6) + 2.0*
           CW8*(8.0*mH2_6 - 24.0*mH2_4*MZ2 + 24.0*mH2_2*MZ4 + 
             6.0*mHp2*MZ4 - MZ6) - 
          CW6*(5.0*mH2_6 - 15.0*mH2_4*MZ2 + 15.0*mH2_2*MZ4 + 
             2.0*(6.0*mHp4*MZ2 + 21.0*mHp2*MZ4 + MZ6))))*
     log(mH3_2/MZ2))/(2304.0*
     CW8*(-1.0 + CW2)*(-1.0 + CW2)*(-1.0 + CW2)*(mH3_2 - 
       mH2_2)*(mH3_2 - mHp2)*MZ8*pi*pi) + (aem*aem*mMU*
     mMU*(-7.0*(mH3_2 - mH2_2)*(mH2_2 - mHp2)*(mH2_2 - 
          mHp2)*(mH2_2 - mHp2)*(mH2_2 - mHp2) + 
       CW2*(mH3_2 - mH2_2)*(mH2_2 - mHp2)*(mH2_2 - 
          mHp2)*(14.0*mH2_4 + 7.0*mHp2*(2.0*mHp2 + 3.0*MZ2) + 
          mH2_2*(-28.0*mHp2 + 57*MZ2)) - 
       4.0*CW10*(mH2_2 - mHp2)*(2.0*mH3_8 - 
          2.0*mH3_6*(4.0*mH2_2 + 3.0*MZ2) + 
          6.0*mH3_2*mH3_2*(2.0*mH2_4 + mH2_2*MZ2 + MZ4) - 
          mH3_2*(8.0*mH2_6 - 6.0*mH2_4*MZ2 + MZ6) + 
          mH2_2*(2.0*mH2_6 - 6.0*mH2_4*MZ2 + 6.0*mH2_2*MZ4 + 
             MZ6)) + 
       CW4*(mH3_2*(-4.0*mH2_8 + 2.0*mH2_6*(8.0*mHp2 - 39*MZ2) + 
             2.0*mH2_2*(8.0*mHp6 + 3.0*mHp4*MZ2) + 
             mH2_4*(-24.0*mHp4 + 114*mHp2*MZ2 + 9*MZ4) - 
             mHp4*(4.0*mHp4 + 42*mHp2*MZ2 + 21.0*MZ4)) + 
          mH2_2*(4.0*mH2_8 + 4.0*mHp8 + 42*mHp6*MZ2 + 
             21.0*mHp4*MZ4 + mH2_6*(-16.0*mHp2 + 78.0*MZ2) + 
             3.0*mH2_4*(8.0*mHp4 - 38*mHp2*MZ2 + MZ4) - 
             2.0*mH2_2*(8.0*mHp6 + 3.0*mHp4*MZ2 + 6.0*mHp2*MZ4))) + 
       2.0*CW8*(8.0*mH3_8*(mH2_2 - mHp2) - 
          8*mH3_6*(mH2_2 - mHp2)*(4.0*mH2_2 + 3.0*MZ2) + 
          24.0*mH3_2*
           mH3_2*(mH2_2 - mHp2)*(2.0*mH2_4 + mH2_2*MZ2 + MZ4) + 
          mH3_2*(-32.0*mH2_8 - 6.0*mHp4*MZ4 - mH2_2*MZ6 + 
             mHp2*MZ6 + 6.0*mH2_4*MZ2*(-4.0*mHp2 + MZ2) + 
             8*mH2_6*(4.0*mHp2 + 3.0*MZ2)) + 
          mH2_2*(8.0*mH2_8 + 6.0*mHp4*MZ4 - mHp2*MZ6 - 
             8*mH2_6*(mHp2 + 3.0*MZ2) + 
             6.0*mH2_4*(4.0*mHp2*MZ2 + 5.0*MZ4) + 
             mH2_2*(-36*mHp2*MZ4 + MZ6))) + 
       CW6*(-5.0*mH3_8*(mH2_2 - mHp2) + 
          5.0*mH3_6*(mH2_2 - mHp2)*(4.0*mH2_2 + 3.0*MZ2) - 
          15.0*mH3_2*
           mH3_2*(mH2_2 - mHp2)*(2.0*mH2_4 + mH2_2*MZ2 + MZ4) + 
          mH3_2*(20.0*mH2_8 + 3.0*mH2_4*MZ2*(mHp2 - 10.0*MZ2) - 
             mH2_6*(20.0*mHp2 + 3.0*MZ2) + 
             2.0*mHp2*MZ2*(6.0*mHp4 + 21.0*mHp2*MZ2 + MZ4) - 
             2.0*mH2_2*(6.0*mHp4*MZ2 + MZ6)) + 
          mH2_2*(-5.0*mH2_8 + mH2_6*(5.0*mHp2 + 3.0*MZ2) - 
             3.0*mH2_4*MZ2*(mHp2 + 13.0*MZ2) - 
             2.0*mHp2*MZ2*(6.0*mHp4 + 21.0*mHp2*MZ2 + MZ4) + 
             mH2_2*(12*mHp4*MZ2 + 69.0*mHp2*MZ4 + 2.0*MZ6))))*
     log(mH2_2/MZ2))/(2304.0*
     CW8*(-1.0 + CW2)*(-1.0 + CW2)*(-1.0 + CW2)*(mH3_2 - 
       mH2_2)*(mH2_2 - mHp2)*MZ8*pi*pi) + (aem*aem*mMU*
     mMU*(-((7.0 - 14.0*CW2 + 4.0*CW4)*
          mH3_8*(mH2_2 - mHp2)) + (7.0 - 14.0*CW2 + 4.0*CW4)*
        mH3_6*(mH2_2 - mHp2)*(4.0*mHp2 + 3.0*CW2*MZ2) - 
       3.0*mH3_2*
        mH3_2*(mH2_2 - mHp2)*(2.0*(7.0 - 14.0*CW2 + 4.0*CW4)*mHp4 + 
          CW2*(-5.0 - 2.0*CW2 + 4.0*CW4)*mHp2*MZ2 + 
          CW4*(7.0 - 14.0*CW2 + 4.0*CW4)*MZ4) + 
       mHp2*((7.0 - 14.0*CW2 + 4.0*CW4)*
           mH2_8 - (7.0 - 14.0*CW2 + 4.0*CW4)*
           mH2_6*(4.0*mHp2 + 3.0*CW2*MZ2) + 
          3.0*mH2_4*(2.0*(7.0 - 14.0*CW2 + 4.0*CW4)*mHp4 + 
             CW2*(-5.0 - 2.0*CW2 + 4.0*CW4)*mHp2*MZ2 + 
             CW4*(7.0 - 14.0*CW2 + 4.0*CW4)*MZ4) + 
          mH2_2*(-5.0*(7.0 - 14.0*CW2 + 4.0*CW4)*mHp6 + 
             6.0*CW2*(25.0 - 32.0*CW2 + 4.0*CW4)*mHp4*MZ2 - 
             3.0*CW4*(29 - 50.0*CW2 + 12*CW4)*mHp2*MZ4 + 
             4.0*CW6*(1.0 + CW2 - 2.0*CW4)*MZ6) + 
          2.0*mHp2*((7.0 - 14.0*CW2 + 4.0*CW4)*mHp6 - 
             3.0*CW2*(19.0 - 26.0*CW2 + 4.0*CW4)*mHp4*MZ2 + 
             3.0*CW4*(13.0 - 20.0*CW2 + 4.0*CW4)*mHp2*MZ4 + 
             2.0*CW6*(-1.0 - CW2 + 2.0*CW4)*MZ6)) + 
       mH3_2*((-7.0 + 14.0*CW2 - 4.0*CW4)*
           mH2_8 + (7.0 - 14.0*CW2 + 4.0*CW4)*
           mH2_6*(4.0*mHp2 + 3.0*CW2*MZ2) - 
          3.0*mH2_4*(2.0*(7.0 - 14.0*CW2 + 4.0*CW4)*mHp4 + 
             CW2*(-5.0 - 2.0*CW2 + 4.0*CW4)*mHp2*MZ2 + 
             CW4*(7.0 - 14.0*CW2 + 4.0*CW4)*MZ4) + 
          mHp2*(-5.0*(7.0 - 14.0*CW2 + 4.0*CW4)*mHp6 + 
             6.0*CW2*(25.0 - 32.0*CW2 + 4.0*CW4)*mHp4*MZ2 - 
             3.0*CW4*(29 - 50.0*CW2 + 12*CW4)*mHp2*MZ4 + 
             4.0*CW6*(1.0 + CW2 - 2.0*CW4)*MZ6) + 
          2.0*mH2_2*(4.0*(7.0 - 14.0*CW2 + 4.0*CW4)*mHp6 - 
             3.0*CW2*(31 - 38*CW2 + 4.0*CW4)*mHp4*MZ2 + 
             6.0*CW4*(8 - 15.0*CW2 + 4.0*CW4)*mHp2*MZ4 + 
             2.0*CW6*(-1.0 - CW2 + 2.0*CW4)*MZ6)))*
     log(mHp2/MZ2))/(2304.0*
     CW8*(-1.0 + CW2)*(-1.0 + CW2)*(-1.0 + CW2)*(mH3_2 - 
       mHp2)*(-mH2_2 + mHp2)*MZ8*pi*pi) + (aem*aem*(mH3_2 - mHp2)*
     mMU*mMU*(mH3_2*mH3_2 - 2.0*mH3_2*mHp2 + mHp4 - CW2*mHp2*MZ2)*
     log(mH3_2/mHp2)*log(mHp2/MZ2))/(128.0*
     CW8*(-1.0 + CW2)*(-1.0 + CW2)*MZ8*pi*pi) + (aem*
     aem*(mH2_2 - mHp2)*mMU*
     mMU*(mH2_4 - 2.0*mH2_2*mHp2 + mHp4 - CW2*mHp2*MZ2)*
     log(mH2_2/mHp2)*log(mHp2/MZ2))/(128.0*
     CW8*(-1.0 + CW2)*(-1.0 + CW2)*MZ8*pi*pi) + (aem*
     aem*(1.0 - 2.0*CW2)*(1.0 - 2.0*CW2)*(5.0 - 16.0*CW2 + 8.0*CW4)*mMU*
     mMU*(-4.0*mHp2 + MZ2)*(MZ4 - MZ2*negsquareroot(-4.0*mHp2*MZ2 + MZ4) + 
       mHp2*(-4.0*MZ2 + 2.0*negsquareroot (-4.0*mHp2*MZ2 + MZ4)))*
     log((2.0*mHp2 - MZ2 + negsquareroot (-4.0*mHp2*MZ2 + MZ4))/(2.0*
         mHp2)))/(1152.0*CW2*(-1.0 + CW2)*(-1.0 + CW2)*(-1.0 + CW2)*
     MZ6*(-2.0*mHp2 + MZ2 - negsquareroot (-4.0*mHp2*MZ2 + MZ4))*pi*pi) - (aem*
     aem*(5.0 - 16.0*CW2 + 8.0*CW4)*mMU*
     mMU*(-4.0*mH3_2*
        mH2_2 + (mH3_2 + mH2_2 - MZ2)*(mH3_2 + mH2_2 - 
          MZ2))*(mH3_2*mH3_2 + 
       mH3_2*(-2.0*mH2_2 - 2.0*MZ2 + 
          negsquareroot(mH3_2*mH3_2 + (mH2_2 - MZ2)*(mH2_2 - MZ2) - 
             2.0*mH3_2*(mH2_2 + MZ2))) + (mH2_2 - MZ2)*(mH2_2 - 
          MZ2 + negsquareroot(mH3_2*mH3_2 + (mH2_2 - MZ2)*(mH2_2 - MZ2) - 
             2.0*mH3_2*(mH2_2 + MZ2))))*
     log((mH3_2 + mH2_2 - MZ2 + 
         negsquareroot (-4.0*mH3_2*       
             mH2_2 + (mH3_2 + mH2_2 - MZ2)*(mH3_2 + mH2_2 - 
               MZ2)))/(2.0*mH3*mH2)))/(1152.0*
     CW2*(-1.0 + CW2)*(-1.0 + CW2)*(-1.0 + CW2)*
     MZ8*(mH3_2 + mH2_2 - MZ2 + 
       negsquareroot (-4.0*mH3_2*
           mH2_2 + (mH3_2 + mH2_2 - MZ2)*(mH3_2 + mH2_2 - MZ2)))*
     pi*pi) + (aem*aem*(7.0 - 14.0*CW2 + 4.0*CW4)*mMU*
     mMU*(-4.0*mH3_2*
        mHp2 + (mH3_2 + mHp2 - CW2*MZ2)*(mH3_2 + mHp2 - 
          CW2*MZ2))*(mH3_2*mH3_2 + 
       mH3_2*(-2.0*mHp2 - 2.0*CW2*MZ2 + 
          negsquareroot (mH3_2*mH3_2 + (mHp2 - CW2*MZ2)*(mHp2 - CW2*MZ2) - 
             2.0*mH3_2*(mHp2 + CW2*MZ2))) + (mHp2 - CW2*MZ2)*(mHp2 - 
          CW2*MZ2 + 
          negsquareroot (mH3_2*mH3_2 + (mHp2 - CW2*MZ2)*(mHp2 - CW2*MZ2) - 
             2.0*mH3_2*(mHp2 + CW2*MZ2))))*
     log((mH3_2 + mHp2 - CW2*MZ2 + 
         negsquareroot (-4.0*mH3_2*
             mHp2 + (mH3_2 + mHp2 - CW2*MZ2)*(mH3_2 + mHp2 - 
               CW2*MZ2)))/(2.0*mH3*mHp)))/(1152.0*
     CW8*(-1.0 + CW2)*(-1.0 + CW2)*(-1.0 + CW2)*
     MZ8*(mH3_2 + mHp2 - CW2*MZ2 + 
       negsquareroot (-4.0*mH3_2*
           mHp2 + (mH3_2 + mHp2 - CW2*MZ2)*(mH3_2 + mHp2 - 
             CW2*MZ2)))*pi*pi) + (aem*aem*(7.0 - 14.0*CW2 + 4.0*CW4)*
     mMU*mMU*(-4.0*mH2_2*
        mHp2 + (mH2_2 + mHp2 - CW2*MZ2)*(mH2_2 + mHp2 - 
          CW2*MZ2))*(mH2_4 + 
       mH2_2*(-2.0*mHp2 - 2.0*CW2*MZ2 + 
          negsquareroot (mH2_4 + (mHp2 - CW2*MZ2)*(mHp2 - CW2*MZ2) - 
             2.0*mH2_2*(mHp2 + CW2*MZ2))) + (mHp2 - CW2*MZ2)*(mHp2 - 
          CW2*MZ2 + 
          negsquareroot (mH2_4 + (mHp2 - CW2*MZ2)*(mHp2 - CW2*MZ2) - 
             2.0*mH2_2*(mHp2 + CW2*MZ2))))*
     log((mH2_2 + mHp2 - CW2*MZ2 + 
         negsquareroot (-4.0*mH2_2*
             mHp2 + (mH2_2 + mHp2 - CW2*MZ2)*(mH2_2 + mHp2 - 
               CW2*MZ2)))/(2.0*mH2*mHp)))/(1152.0*
     CW8*(-1.0 + CW2)*(-1.0 + CW2)*(-1.0 + CW2)*
     MZ8*(mH2_2 + mHp2 - CW2*MZ2 + 
       negsquareroot(-4.0*mH2_2*
           mHp2 + (mH2_2 + mHp2 - CW2*MZ2)*(mH2_2 + mHp2 - 
             CW2*MZ2)))*pi*pi) - (aem*aem*(mH3_2 - mHp2)*mMU*
     mMU*(mH3_2*mH3_2 - 2.0*mH3_2*mHp2 + mHp4 - CW2*mHp2*MZ2)*PolyLog.Li2( 1.0 - mH3_2/mHp2))/(64.0*
     CW8*(-1.0 + CW2)*(-1.0 + CW2)*MZ8*pi*pi) - (aem*
     aem*(mH2_2 - mHp2)*mMU*
     mMU*(mH2_4 - 2.0*mH2_2*mHp2 + mHp4 - CW2*mHp2*MZ2)*
     PolyLog.Li2( 1.0 - mH2_2/mHp2))/(64.0*
     CW8*(-1.0 + CW2)*(-1.0 + CW2)*MZ8*pi*pi) + (aem*
     aem*(mH3_2 - mHp2)*mMU*
     mMU*(-mH3_6 + 
       mH3_2*mH3_2*(3.0*mHp2 + CW2*MZ2) + (mHp3 - 
          CW2*mHp*MZ2)*(mHp3 - CW2*mHp*MZ2) + 
       mH3_2*(-3.0*mHp4 + CW2*mHp2*MZ2))*TF (mH3, mHp, CW*MZ))/(64.0*
     CW8*(-1.0 + CW2)*(-1.0 + CW2)*
     MZ8*(mH3_2*mH3_2 + (mHp2 - CW2*MZ2)*(mHp2 - CW2*MZ2) - 
       2.0*mH3_2*(mHp2 + CW2*MZ2))*pi*pi) + (aem*aem*(mH2_2 - mHp2)*
     mMU*mMU*(-mH2_6 + 
       mH2_4*(3.0*mHp2 + CW2*MZ2) + (mHp3 - CW2*mHp*MZ2)*(mHp3 - 
          CW2*mHp*MZ2) + mH2_2*(-3.0*mHp4 + CW2*mHp2*MZ2))*
     TF(mH2, mHp, CW*MZ))/(64.0*CW8*(-1.0 + CW2)*(-1.0 + CW2)*
     MZ8*(mH2_4 + (mHp2 - CW2*MZ2)*(mHp2 - CW2*MZ2) - 
       2.0*mH2_2*(mHp2 + CW2*MZ2))*pi*pi);


/*std::cout << "exp  = " << (aem*aem*mMU*
     mMU*((-7.0 + 14.0*CW2 - 4.0*CW4 + 5.0*CW6 - 16.0*CW8 + 8.0*CW10)*
        mH3_10 + 
       mH3_8*((7.0 - 14.0*CW2 + 4.0*CW4 - 20.0*CW6 + 64.0*CW8 - 
             32.0*CW10)*
           mH2_2 + (28.0 - 56.0*CW2 + 16.0*CW4 - 5.0*CW6 + 
             16.0*CW8 - 8.0*CW10)*mHp2 - 
          3.0*CW2*(-19.0 + 26.0*CW2 + CW4 - 16.0*CW6 + 8.0*CW8)*MZ2) +
        mH3_6*(-14.0*mHp2*(2.0*mH2_2 + 3.0*mHp2) + 
          CW2*(84.0*mHp4 - 93.0*mHp2*MZ2 + 
             mH2_2*(56.0*mHp2 - 57*MZ2)) + 
          8.0*CW10*(6.0*mH2_4 + 3.0*MZ2*(mHp2 + MZ2) + 
             mH2_2*(4.0*mHp2 + 3.0*MZ2)) - 
          4.0*CW8*(24.0*mH2_4 + 4.0*mH2_2*(4.0*mHp2 + 3.0*MZ2) + 
             3.0*MZ2*(4.0*mHp2 + 5.0*MZ2)) + 
    CW6*(30.0*mH2_4 + mH2_2*(20.0*mHp2 + 3.0*MZ2) + 
             3.0*MZ2*(mHp2 + 13.0*MZ2)) + 
          CW4*(mH2_2*(-16.0*mHp2 + 78.0*MZ2) - 
             3.0*(8.0*mHp4 - 38*mHp2*MZ2 + MZ4))) + 
       mH3_2*mH3_2*(42*mH2_2*mHp4 + 28.0*mHp6 - 
          4.0*CW10*(8.0*mH2_6 + 6.0*mH2_2*mHp2*MZ2 + 
             6.0*mHp2*MZ4 - MZ6 + 6.0*mH2_4*(2.0*mHp2 - MZ2)) - 
          CW6*(20.0*mH2_6 + 12*mHp4*MZ2 + 69.0*mHp2*MZ4 + 
             2.0*MZ6 + 3.0*mH2_2*MZ2*(mHp2 - 10.0*MZ2) + 
             15.0*mH2_4*(2.0*mHp2 - MZ2)) + 
          CW2*(-56.0*mHp6 + 15.0*mHp4*MZ2 + 
             mH2_2*(-84.0*mHp4 + 93.0*mHp2*MZ2)) + 
          2.0*CW8*(32.0*mH2_6 + 36*mHp2*MZ4 - MZ6 + 
             24.0*mH2_4*(2.0*mHp2 - MZ2) + 
             6.0*mH2_2*(4.0*mHp2*MZ2 - MZ4)) + 
          CW4*(3.0*mH2_2*(8.0*mHp4 - 38*mHp2*MZ2 - 3.0*MZ4) + 
             2.0*mHp2*(8.0*mHp4 + 3.0*mHp2*MZ2 + 6.0*MZ4))) + 
       mH3_2*(-7.0*(4.0*mH2_2*mHp6 + mHp8) - 
          2.0*CW8*(8.0*mH2_8 + 6.0*mHp4*MZ4 - mH2_2*MZ6 - 
             mHp2*MZ6 + 8*mH2_6*(4.0*mHp2 - 3.0*MZ2) + 
             24.0*mH2_4*MZ2*(-mHp2 + MZ2)) - 
          CW4*mHp4*(4.0*mHp4 + 42*mHp2*MZ2 + 21.0*MZ4 + 
             2.0*mH2_2*(8.0*mHp2 + 3.0*MZ2)) + 
          CW2*(7.0*mHp6*(2.0*mHp2 + 3.0*MZ2) + 
             mH2_2*(56.0*mHp6 - 15.0*mHp4*MZ2)) + 
          4.0*CW10*(2.0*mH2_8 - mH2_2*MZ6 - mHp2*MZ6 + 
             mH2_6*(8.0*mHp2 - 6.0*MZ2) + 
             mH2_4*(-6.0*mHp2*MZ2 + 6.0*MZ4)) + 
          CW6*(5.0*mH2_8 + 5.0*mH2_6*(4.0*mHp2 - 3.0*MZ2) + 
             15.0*mH2_4*MZ2*(-mHp2 + MZ2) + 
             2.0*mHp2*MZ2*(6.0*mHp4 + 21.0*mHp2*MZ2 + MZ4) + 
             2.0*mH2_2*(6.0*mHp4*MZ2 + MZ6))) + 
       mH2_2*mHp2*(7.0*mHp6 - 7.0*CW2*(2.0*mHp6 + 3.0*mHp4*MZ2) + 
          CW4*(4.0*mHp6 + 42*mHp4*MZ2 + 21.0*mHp2*MZ4) - 
          4.0*CW10*(2.0*mH2_6 - 6.0*mH2_4*MZ2 + 6.0*mH2_2*MZ4 - 
             MZ6) + 2.0*
           CW8*(8.0*mH2_6 - 24.0*mH2_4*MZ2 + 24.0*mH2_2*MZ4 + 
             6.0*mHp2*MZ4 - MZ6) - 
          CW6*(5.0*mH2_6 - 15.0*mH2_4*MZ2 + 15.0*mH2_2*MZ4 + 
             2.0*(6.0*mHp4*MZ2 + 21.0*mHp2*MZ4 + MZ6))))*
     log(mH3_2/MZ2))/(2304.0*
     CW8*(-1.0 + CW2)*(-1.0 + CW2)*(-1.0 + CW2)*(mH3_2 - 
       mH2_2)*(mH3_2 - mHp2)*MZ8*pi*pi) << std::endl;

std::cout << "num = " << (aem*aem*mMU*
     mMU*((-7.0 + 14.0*CW2 - 4.0*CW4 + 5.0*CW6 - 16.0*CW8 + 8.0*CW10)*
        mH3_10 + 
       mH3_8*((7.0 - 14.0*CW2 + 4.0*CW4 - 20.0*CW6 + 64.0*CW8 - 
             32.0*CW10)*
           mH2_2 + (28.0 - 56.0*CW2 + 16.0*CW4 - 5.0*CW6 + 
             16.0*CW8 - 8.0*CW10)*mHp2 - 
          3.0*CW2*(-19.0 + 26.0*CW2 + CW4 - 16.0*CW6 + 8.0*CW8)*MZ2) +
        mH3_6*(-14.0*mHp2*(2.0*mH2_2 + 3.0*mHp2) + 
          CW2*(84.0*mHp4 - 93.0*mHp2*MZ2 + 
             mH2_2*(56.0*mHp2 - 57*MZ2)) + 
          8.0*CW10*(6.0*mH2_4 + 3.0*MZ2*(mHp2 + MZ2) + 
             mH2_2*(4.0*mHp2 + 3.0*MZ2)) - 
          4.0*CW8*(24.0*mH2_4 + 4.0*mH2_2*(4.0*mHp2 + 3.0*MZ2) + 
             3.0*MZ2*(4.0*mHp2 + 5.0*MZ2)) + 
    CW6*(30.0*mH2_4 + mH2_2*(20.0*mHp2 + 3.0*MZ2) + 
             3.0*MZ2*(mHp2 + 13.0*MZ2)) + 
          CW4*(mH2_2*(-16.0*mHp2 + 78.0*MZ2) - 
             3.0*(8.0*mHp4 - 38*mHp2*MZ2 + MZ4))) + 
       mH3_2*mH3_2*(42*mH2_2*mHp4 + 28.0*mHp6 - 
          4.0*CW10*(8.0*mH2_6 + 6.0*mH2_2*mHp2*MZ2 + 
             6.0*mHp2*MZ4 - MZ6 + 6.0*mH2_4*(2.0*mHp2 - MZ2)) - 
          CW6*(20.0*mH2_6 + 12*mHp4*MZ2 + 69.0*mHp2*MZ4 + 
             2.0*MZ6 + 3.0*mH2_2*MZ2*(mHp2 - 10.0*MZ2) + 
             15.0*mH2_4*(2.0*mHp2 - MZ2)) + 
          CW2*(-56.0*mHp6 + 15.0*mHp4*MZ2 + 
             mH2_2*(-84.0*mHp4 + 93.0*mHp2*MZ2)) + 
          2.0*CW8*(32.0*mH2_6 + 36*mHp2*MZ4 - MZ6 + 
             24.0*mH2_4*(2.0*mHp2 - MZ2) + 
             6.0*mH2_2*(4.0*mHp2*MZ2 - MZ4)) + 
          CW4*(3.0*mH2_2*(8.0*mHp4 - 38*mHp2*MZ2 - 3.0*MZ4) + 
             2.0*mHp2*(8.0*mHp4 + 3.0*mHp2*MZ2 + 6.0*MZ4))) + 
       mH3_2*(-7.0*(4.0*mH2_2*mHp6 + mHp8) - 
          2.0*CW8*(8.0*mH2_8 + 6.0*mHp4*MZ4 - mH2_2*MZ6 - 
             mHp2*MZ6 + 8*mH2_6*(4.0*mHp2 - 3.0*MZ2) + 
             24.0*mH2_4*MZ2*(-mHp2 + MZ2)) - 
          CW4*mHp4*(4.0*mHp4 + 42*mHp2*MZ2 + 21.0*MZ4 + 
             2.0*mH2_2*(8.0*mHp2 + 3.0*MZ2)) + 
          CW2*(7.0*mHp6*(2.0*mHp2 + 3.0*MZ2) + 
             mH2_2*(56.0*mHp6 - 15.0*mHp4*MZ2)) + 
          4.0*CW10*(2.0*mH2_8 - mH2_2*MZ6 - mHp2*MZ6 + 
             mH2_6*(8.0*mHp2 - 6.0*MZ2) + 
             mH2_4*(-6.0*mHp2*MZ2 + 6.0*MZ4)) + 
          CW6*(5.0*mH2_8 + 5.0*mH2_6*(4.0*mHp2 - 3.0*MZ2) + 
             15.0*mH2_4*MZ2*(-mHp2 + MZ2) + 
             2.0*mHp2*MZ2*(6.0*mHp4 + 21.0*mHp2*MZ2 + MZ4) + 
             2.0*mH2_2*(6.0*mHp4*MZ2 + MZ6))) + 
       mH2_2*mHp2*(7.0*mHp6 - 7.0*CW2*(2.0*mHp6 + 3.0*mHp4*MZ2) + 
          CW4*(4.0*mHp6 + 42*mHp4*MZ2 + 21.0*mHp2*MZ4) - 
          4.0*CW10*(2.0*mH2_6 - 6.0*mH2_4*MZ2 + 6.0*mH2_2*MZ4 - 
             MZ6) + 2.0*
           CW8*(8.0*mH2_6 - 24.0*mH2_4*MZ2 + 24.0*mH2_2*MZ4 + 
             6.0*mHp2*MZ4 - MZ6) - 
          CW6*(5.0*mH2_6 - 15.0*mH2_4*MZ2 + 15.0*mH2_2*MZ4 + 
             2.0*(6.0*mHp4*MZ2 + 21.0*mHp2*MZ4 + MZ6))))*
     log(mH3_2/MZ2)) << std::endl; 


std::cout << "den = " << (2304.0*
     CW8*(-1.0 + CW2)*(-1.0 + CW2)*(-1.0 + CW2)*(mH3_2 - 
       mH2_2)*(mH3_2 - mHp2)*MZ8*pi*pi)<< std::endl; */






gslpp::complex aYuk2 =  -(aem*aem*(MH2 + 2.0*mHp2)*mMU*mMU)/(8.0*CW2*(-1.0 + CW2)*MH2*MZ2*pi*
    pi) + (aem*aem*(MH2 + 2.0*mHp2)*mMU*
    mMU*((-3.0 + 2.0*CW2)*MH2 - 8.0*CW2*(-1.0 + CW2)*MZ2)*
    log (MH2/MZ2))/(128.0*CW4*(-1.0 + CW2)*(-1.0 + CW2)*MH2*
    MZ2*(MH2 - MZ2)*pi*pi) - (aem*aem*(MH2 + 2.0*mHp2)*mMU*mMU*
    log (mHp2/MZ2))/(16.0*CW2*(-1.0 + CW2)*MH2*MZ2*pi*pi) + (aem*aem*
    mHp2*(MH2 + 2.0*mHp2)*mMU*
    mMU*((-3.0 + 2.0*CW2)*MH2 - 8.0*CW2*(-1.0 + CW2)*MZ2)*
    TF (MH, mHp, mHp))/(64.0*CW4*(-1.0 + CW2)*(-1.0 + CW2)*
    MH4*(MH2 - 4.0*mHp2)*MZ2*(MH2 - MZ2)*pi*pi) - (aem*
    aem*(3.0 - 10.0*CW2 + 8.0*CW4)*mHp2*(MH2 + 2.0*mHp2)*mMU*mMU*
    TF (mHp, mHp, MZ))/(64.0*CW4*(-1.0 + CW2)*(-1.0 + CW2)*
    MZ4*(-MH2 + MZ2)*(-4.0*mHp2 + MZ2)*pi*pi);


gslpp::complex aYuk3 = sl*((aem*aem*mMU*
       mMU*(-mH2_4 + mH2_2*mHp2 + 
         8.0*CW2*(-1.0 + CW2)*mHp2*MZ2)*(-1.0 + tanb*tanb))/(64.0*
       CW4*(-1.0 + CW2)*(-1.0 + CW2)*mHp2*MZ4*pi*pi*tanb) - (aem*aem*
       mH3_2*mH2_2*mHp2*mMU*mMU*(-1.0 + tanb*tanb)*
       log(mH3_2/mHp2)*
       log(mH3_2/mHp2))/(512.0*
       CW2*(-1.0 + CW2)*(-1.0 + CW2)*
       MZ2*(-mHp2 + CW2*MZ2)*(-mHp2 + CW2*MZ2)*(-mHp2 + CW2*MZ2)*pi*
       pi*tanb) + (aem*aem*mH2_4*mMU*
       mMU*(-7.0*mHp6 + 8.0*CW2*mHp4*MZ2 - 4.0*CW4*mHp2*MZ4 + 
         2.0*mH2_2*(mHp2 - CW2*MZ2)*(mHp2 - CW2*MZ2))*(-1.0 + 
         tanb*tanb)*log(mH2_2/mHp2)*log(mH2_2/mHp2))/(512.0*
       CW2*(-1.0 + CW2)*(-1.0 + CW2)*mHp4*
       MZ2*(-mHp2 + CW2*MZ2)*(-mHp2 + CW2*MZ2)*(-mHp2 + CW2*MZ2)*pi*
       pi*tanb) + (aem*aem*mH3_2*mH2_2*mHp2*mMU*
       mMU*(-1.0 + tanb*tanb)*log (mH3_2/MZ2)*log (mH3_2/MZ2))/(512.0*
       CW2*(-1.0 + CW2)*(-1.0 + CW2)*
       MZ2*(-mHp2 + CW2*MZ2)*(-mHp2 + CW2*MZ2)*(-mHp2 + CW2*MZ2)*pi*
       pi*tanb) - (aem*aem*mH2_4*mMU*
       mMU*(-7.0*mHp6 + 8.0*CW2*mHp4*MZ2 - 4.0*CW4*mHp2*MZ4 + 
         2.0*mH2_2*(mHp2 - CW2*MZ2)*(mHp2 - CW2*MZ2))*(-1.0 + 
         tanb*tanb)*log (mH2_2/MZ2)*log (mH2_2/MZ2))/(512.0*
       CW2*(-1.0 + CW2)*(-1.0 + CW2)*mHp4*
       MZ2*(-mHp2 + CW2*MZ2)*(-mHp2 + CW2*MZ2)*(-mHp2 + CW2*MZ2)*pi*
       pi*tanb) + (aem*aem*mMU*
       mMU*(2.0*CW2*mH2_4*MZ2 + 
         8.0*CW2*(-1.0 + CW2)*mHp2*MZ2*(-mHp2 + CW2*MZ2) + 
         mH2_2*(-2.0*mHp4 + CW2*mHp2*MZ2))*(-1.0 + tanb*tanb)*
       log (mHp2/MZ2))/(128.0*CW4*(-1.0 + CW2)*(-1.0 + CW2)*mHp2*
       MZ4*(-mHp2 + CW2*MZ2)*pi*pi*tanb) + (aem*aem*mH2_2*mMU*
       mMU*(mHp6 + mH2_4*(mHp2 + CW2*MZ2) - 
         2.0*mH2_2*(mHp4 + CW2*mHp2*MZ2))*(-1.0 + tanb*tanb)*
       log (mH2_2/mHp2)*log (mHp2/MZ2))/(128.0*
       CW6*(-1.0 + CW2)*(-1.0 + CW2)*mHp4*MZ6*pi*pi*tanb) - (aem*aem*
       mH3_2*mH2_2*mHp2*mMU*mMU*(-1.0 + tanb*tanb)*log (mH3_2/MZ2)*
       log (mHp2/MZ2))/(256.0*CW2*(-1.0 + CW2)*(-1.0 + CW2)*
       MZ2*(-mHp2 + CW2*MZ2)*(-mHp2 + CW2*MZ2)*(-mHp2 + CW2*MZ2)*pi*
       pi*tanb) + (aem*aem*mH2_2*mMU*
       mMU*(mH3_2*mHp6 + 
         2.0*mH2_4*(mHp2 - CW2*MZ2)*(mHp2 - CW2*MZ2) - 
         mH2_2*(mHp6 - 8.0*CW2*mHp4*MZ2 + 4.0*CW4*mHp2*MZ4))*(-1.0 + 
         tanb*tanb)*log (mHp2/MZ2)*log (mHp2/MZ2))/(512.0*
       CW2*(-1.0 + CW2)*(-1.0 + CW2)*mHp4*
       MZ2*(-mHp2 + CW2*MZ2)*(-mHp2 + CW2*MZ2)*(-mHp2 + CW2*MZ2)*pi*
       pi*tanb)+log(mH2_2/
        MZ2)*((aem*aem*mMU*
          mMU*(2.0*mH2_6 - 
            2.0*mH2_4*MZ2 + (-3.0 + 2.0*CW2)*mH2_2*mHp2*MZ2 - 
            8.0*CW2*(-1.0 + CW2)*mHp2*MZ4)*(-1.0 + tanb*tanb))/(128.0*
          CW4*(-1.0 + CW2)*(-1.0 + CW2)*mHp2*MZ4*(-mH2_2 + MZ2)*pi*
          pi*tanb) - (3.0*aem*aem*mH2_4*mHp2*mMU*
          mMU*(-1.0 + tanb*tanb)*log(mHp2/MZ2))/(256.0*
          CW2*(-1.0 + CW2)*(-1.0 + CW2)*
          MZ2*(-mHp2 + CW2*MZ2)*(-mHp2 + CW2*MZ2)*(-mHp2 + CW2*MZ2)*
          pi*pi*tanb))+ log (CW2)*(-(aem*aem*mH2_2*mMU*
           mMU*(2.0*mH2_2 - 2.0*mHp2 + CW2*MZ2)*(-1.0 + 
             tanb*tanb))/(128.0*CW4*(-1.0 + CW2)*(-1.0 + CW2)*
          MZ4*(-mHp2 + CW2*MZ2)*pi*pi*tanb) - (aem*aem*mH2_2*mMU*
          mMU*(mHp6 + mH2_4*(mHp2 + CW2*MZ2) - 
            2.0*mH2_2*(mHp4 + CW2*mHp2*MZ2))*(-1.0 + tanb*tanb)*
          log (mH2_2/mHp2))/(128.0*CW6*(-1.0 + CW2)*(-1.0 + CW2)*
          mHp4*MZ6*pi*pi*tanb) + (aem*aem*mH2_4*(mH2_2 - 2.0*mHp2)*
          mMU*mMU*(-1.0 + tanb*tanb)*log (mH2_2/MZ2))/(128.0*
          CW2*(-1.0 + CW2)*(-1.0 + CW2)*mHp4*MZ2*(-mHp2 + CW2*MZ2)*pi*
          pi*tanb) - (aem*aem*mH2_4*(mH2_2 - 2.0*mHp2)*mMU*
          mMU*(-1.0 + tanb*tanb)*log (mHp2/MZ2))/(128.0*
          CW2*(-1.0 + CW2)*(-1.0 + CW2)*mHp4*MZ2*(-mHp2 + CW2*MZ2)*pi*
          pi*tanb)) - (aem*aem*mH2_2*mMU*
       mMU*(mHp6 + mH2_4*(mHp2 + CW2*MZ2) - 
         2.0*mH2_2*(mHp4 + CW2*mHp2*MZ2))*(-1.0 + tanb*tanb)*
       PolyLog.Li2 (1.0 - mH2_2/mHp2))/(64.0*
       CW6*(-1.0 + CW2)*(-1.0 + CW2)*mHp4*MZ6*pi*pi*tanb)+ (aem*aem*
       mMU*mMU*(-3.0*mH2_2*mHp8 + 8.0*CW6*mHp6*MZ4 - 
         2.0*CW4*mHp6*MZ2*(mH2_2 + 4.0*(mHp2 + MZ2)) + 
         CW2*(-mH2_10 + 2.0*mH2_4*mHp4*MZ2 + 8.0*mHp8*MZ2 + 
            mH2_8*(4.0*mHp2 + MZ2) - 
            2.0*mH2_6*(mHp4 + 2.0*mHp2*MZ2) + 
            mH2_2*(2.0*mHp8 + 3.0*mHp6*MZ2)))*(-1.0 + tanb*tanb)*
       TF (mH2, mHp, mHp))/(64.0*CW4*(-1.0 + CW2)*(-1.0 + CW2)*mH2_2*
       mHp4*(mH2_2 - 4.0*mHp2)*MZ2*(mH2_2 - MZ2)*(-mHp2 + CW2*MZ2)*
       pi*pi*tanb) - (aem*aem*mH2_2*mMU*
       mMU*(-mH2_6 + 
         mH2_4*(3.0*mHp2 + CW2*MZ2) + (mHp3 - CW2*mHp*MZ2)*(mHp3 - 
            CW2*mHp*MZ2) + mH2_2*(-3.0*mHp4 + CW2*mHp2*MZ2))*(-1.0 + 
         tanb*tanb)*TF(mH2, mHp, CW*MZ))/(64.0*
       CW6*(-1.0 + CW2)*(-1.0 + CW2)*
       MZ6*(-mHp2 + 
         CW2*MZ2)*(mH2_4 + (mHp2 - CW2*MZ2)*(mHp2 - CW2*MZ2) - 
         2.0*mH2_2*(mHp2 + CW2*MZ2))*pi*pi*tanb)+  (aem*
       aem*(3.0 - 10.0*CW2 + 8.0*CW4)*mH2_2*mHp2*mMU*
       mMU*(-1.0 + tanb*tanb)*TF (mHp, mHp, MZ))/(64.0*
       CW4*(-1.0 + CW2)*(-1.0 + CW2)*
       MZ4*(-mH2_2 + MZ2)*(-4.0*mHp2 + MZ2)*pi*pi*tanb));


gslpp::complex aYuk4 = eta*((aem*aem*mMU*
       mMU*((-3.0 + 2.0*CW2)*MH2 - 8.0*CW2*(-1.0 + CW2)*MZ2)*(-1.0 + 
         tanb*tanb)*log(MH2/MZ2))/(128.0*
       CW4*(-1.0 + CW2)*(-1.0 + CW2)*MZ2*(-MH2 + MZ2)*pi*pi*
       tanb) + (aem*aem*mMU*
       mMU*((3.0 - 2.0*CW2)*mH2_2 + 8.0*CW2*(-1.0 + CW2)*MZ2)*(-1.0 +
          tanb*tanb)*log(mH2_2/MZ2))/(128.0*
       CW4*(-1.0 + CW2)*(-1.0 + CW2)*MZ2*(-mH2_2 + MZ2)*pi*pi*
       tanb) + (aem*aem*mHp2*mMU*
       mMU*((3.0 - 2.0*CW2)*MH2 + 8.0*CW2*(-1.0 + CW2)*MZ2)*(-1.0 + 
         tanb*tanb)*TF (MH, mHp, mHp))/(64.0*
       CW4*(-1.0 + CW2)*(-1.0 + CW2)*MH2*(MH2 - 4.0*mHp2)*
       MZ2*(MH2 - MZ2)*pi*pi*tanb) + (aem*aem*mHp2*mMU*
       mMU*((-3.0 + 2.0*CW2)*mH2_2 - 
         8.0*CW2*(-1.0 + CW2)*MZ2)*(-1.0 + tanb*tanb)*
       TF (mH2, mHp, mHp))/(64.0*CW4*(-1.0 + CW2)*(-1.0 + CW2)*
       mH2_2*(mH2_2 - 4.0*mHp2)*MZ2*(mH2_2 - MZ2)*pi*pi*
       tanb) + (aem*aem*(3.0 - 10.0*CW2 + 8.0*CW4)*(MH2 - mH2_2)*
       mHp2*mMU*mMU*(-1.0 + tanb*tanb)*TF (mHp, mHp, MZ))/(64.0*
       CW4*(-1.0 + CW2)*(-1.0 + CW2)*
       MZ2*(-MH2 + MZ2)*(-mH2_2 + MZ2)*(-4.0*mHp2 + MZ2)*pi*pi*
       tanb)+ sl*((aem*aem*mMU*
          mMU*(144.0*MH6*mH2_4*MZ4 - 144.0*CW2*MH4*mH2_4*MZ6 - 
            1152.0*CW2*(-1.0 + CW2)*mH2_4*mHp4*MZ6 + 
            MH2*(3.0*(40.0 - 41.0*CW2 + 64.0*CW4)*mH2_4*mHp2*MZ6 + 
               144.0*mH2_6*
                MZ4*((1.0 - 3.0*CW2 + 2.0*CW4)*mHp2 + CW2*MZ2) + 
               24.0*(1.0 - 3.0*CW2 + 2.0*CW4)*mH2_10*mHp2*pi*pi + 
               16.0*(5.0 - 12.0*CW2 + 8.0*CW4 + 2.0*CW6)*mHp2*MZ10*pi*
                pi - 48.0*mH2_8*
                MZ2*(3.0*MZ2 + 
                  2.0*(1.0 - 3.0*CW2 + 2.0*CW4)*mHp2*pi*pi) + 
               12.0*mH2_2*mHp2*
                MZ6*(288*CW6*MZ2 + 5.0*MZ2*(-8.0 + pi*pi) - 
                  12.0*CW2*(8.0*mHp2 + MZ2*(-8.0 + pi*pi)) +  2.0*CW4*(48.0*mHp2 + 
                    MZ2*(-184.0 + 5.0*pi*pi))))))/(4608.0*
          CW4*(-1.0 + CW2)*(-1.0 + CW2)*MH2*mH2_4*mHp2*MZ8*pi*
          pi) + (aem*aem*(MH2 - mH2_2)*mMU*
          mMU*(MH4 + mH2_4 - 3.0*CW2*mHp2*MZ2 + 3.0*CW4*MZ4 + 
            MH2*(mH2_2 - mHp2 - 3.0*CW2*MZ2) - 
            mH2_2*(mHp2 + 3.0*CW2*MZ2))*log (CW2)*log (CW2))/(128.0*
          CW2*(-1.0 + CW2)*(-1.0 + CW2)*mHp4*MZ2*(-mHp2 + CW2*MZ2)*pi*
          pi) +  (aem*aem*mH3_2*(MH2 - mH2_2)*mHp2*mMU*mMU*
          log (mH3_2/mHp2)*
          log (mH3_2/mHp2))/(512.0*
          CW2*(-1.0 + CW2)*(-1.0 + CW2)*
          MZ2*(-mHp2 + CW2*MZ2)*(-mHp2 + CW2*MZ2)*(-mHp2 + CW2*MZ2)*
          pi*pi) + (aem*aem*mMU*
          mMU*(-4.0*MH6*(mHp2 - CW2*MZ2)*(mHp2 - CW2*MZ2) + 
            2.0*CW2*
             MZ2*(-3.0*mHp2 + CW2*MZ2)*(-mHp2 + CW2*MZ2)*(-mHp2 + 
               CW2*MZ2)*(-mHp2 + CW2*MZ2) + 
            MH4*(7.0*mHp6 + 2.0*CW2*mHp4*MZ2 - 10.0*CW4*mHp2*MZ4 + 
               6.0*CW6*MZ6) + 
            MH2*(mH2_2*mHp6 + 
               2.0*(mHp2 - CW2*MZ2)*(mHp2 - CW2*MZ2)*(4.0*mHp4 + 
                  3.0*CW2*mHp2*MZ2 - 3.0*CW4*MZ4)))*
          log (MH2/mHp2)*
          log (MH2/mHp2))/(512.0*
          CW2*(-1.0 + CW2)*(-1.0 + CW2)*mHp4*
          MZ2*(-mHp2 + CW2*MZ2)*(-mHp2 + CW2*MZ2)*(-mHp2 + CW2*MZ2)*
          pi*pi) +   (aem*aem*mMU*
          mMU*(4.0*mH2_6*(mHp2 - CW2*MZ2)*(mHp2 - CW2*MZ2) - 
            2.0*CW2*
             MZ2*(-3.0*mHp2 + CW2*MZ2)*(-mHp2 + CW2*MZ2)*(-mHp2 + 
               CW2*MZ2)*(-mHp2 + CW2*MZ2) - 
            mH2_4*(9.0*mHp6 + 2.0*CW2*mHp4*MZ2 - 10.0*CW4*mHp2*MZ4 + 
               6.0*CW6*MZ6) + 
            mH2_2*(MH2*mHp6 - 
               2.0*(mHp2 - CW2*MZ2)*(mHp2 - CW2*MZ2)*(4.0*mHp4 + 
                  3.0*CW2*mHp2*MZ2 - 3.0*CW4*MZ4)))*
          log (mH2_2/mHp2)*
          log (mH2_2/mHp2))/(512.0*
          CW2*(-1.0 + CW2)*(-1.0 + CW2)*mHp4*
          MZ2*(-mHp2 + CW2*MZ2)*(-mHp2 + CW2*MZ2)*(-mHp2 + CW2*MZ2)*
          pi*pi) +  (aem*aem*mH3_2*(-MH2 + mH2_2)*mHp2*mMU*mMU*
          log (mH3_2/MZ2)*
          log (mH3_2/MZ2))/(512.0*
          CW2*(-1.0 + CW2)*(-1.0 + CW2)*
          MZ2*(-mHp2 + CW2*MZ2)*(-mHp2 + CW2*MZ2)*(-mHp2 + CW2*MZ2)*
          pi*pi) - (aem*aem*mMU*
          mMU*(-4.0*MH6*(mHp2 - CW2*MZ2)*(mHp2 - CW2*MZ2) + 
            2.0*CW2*
             MZ2*(-3.0*mHp2 + CW2*MZ2)*(-mHp2 + CW2*MZ2)*(-mHp2 + 
               CW2*MZ2)*(-mHp2 + CW2*MZ2) + 
            MH4*(7.0*mHp6 + 2.0*CW2*mHp4*MZ2 - 10.0*CW4*mHp2*MZ4 + 
               6.0*CW6*MZ6) + 
            MH2*(mH2_2*mHp6 + 
               2.0*(mHp2 - CW2*MZ2)*(mHp2 - CW2*MZ2)*(4.0*mHp4 + 
                  3.0*CW2*mHp2*MZ2 - 3.0*CW4*MZ4)))*
          log (MH2/MZ2)*
          log (MH2/MZ2))/(512.0*
          CW2*(-1.0 + CW2)*(-1.0 + CW2)*mHp4*
          MZ2*(-mHp2 + CW2*MZ2)*(-mHp2 + CW2*MZ2)*(-mHp2 + CW2*MZ2)*
          pi*pi) +  (aem*aem*mMU*
          mMU*(2.0*CW4*
             MZ8*(-3.0*mHp2 + CW2*MZ2)*(-mHp2 + CW2*MZ2)*(-mHp2 + 
               CW2*MZ2)*(-mHp2 + CW2*MZ2) - 
            4.0*mH2_6*(mHp2 - CW2*MZ2)*(mHp2 - 
               CW2*MZ2)*((2.0 - 6.0*CW2 + 4.0*CW4)*mHp6 - 
               2.0*CW2*(1.0 - 3.0*CW2 + 2.0*CW4)*mHp4*MZ2 + CW2*MZ6) +
             mH2_4*
             MZ2*(32*(1.0 - 3.0*CW2 + 2.0*CW4)*mHp10 - 
               96.0*CW2*(1.0 - 3.0*CW2 + 2.0*CW4)*mHp8*MZ2 + 
               3.0*CW2*(3.0 + 32*CW2 - 96.0*CW4 + 64.0*CW6)*mHp6*
                MZ4 + 2.0*CW4*(1.0 - 16.0*CW2 + 48.0*CW4 - 32*CW6)*
                mHp4*MZ6 - 10.0*CW6*mHp2*MZ8 + 6.0*CW8*MZ10) - 
            CW2*mH2_2*
             MZ6*(MH2*mHp6 - 
               2.0*(mHp2 - CW2*MZ2)*(mHp2 - CW2*MZ2)*(4.0*mHp4 + 
                  3.0*CW2*mHp2*MZ2 - 3.0*CW4*MZ4)))*
          log(mH2_2/MZ2)*
          log(mH2_2/MZ2))/(512.0*
          CW4*(-1.0 + CW2)*(-1.0 + CW2)*mHp4*
          MZ8*(-mHp2 + CW2*MZ2)*(-mHp2 + CW2*MZ2)*(-mHp2 + CW2*MZ2)*
          pi*pi) +   (aem*aem*(MH2 - mH2_2)*mMU*
          mMU*(-2.0*CW2*MH4*mH2_2*MZ2 + 
            8.0*CW2*(-1.0 + CW2)*mHp4*MZ2*(-mHp2 + CW2*MZ2) + 
            MH2*mH2_2*(mHp4 - 2.0*CW2*mH2_2*MZ2 - 
               2.0*CW2*mHp2*MZ2 + 2.0*CW4*MZ4))*log (mHp2/MZ2))/(64.0*
          CW4*(-1.0 + CW2)*(-1.0 + CW2)*MH2*mH2_2*mHp2*
          MZ4*(-mHp2 + CW2*MZ2)*pi*pi)  - (aem*aem*mMU*
          mMU*(2.0*mHp8 + MH6*(mHp2 + CW2*MZ2) - 
            MH2*(3.0*mHp6 + 4.0*CW2*mHp4*MZ2))*log (MH2/mHp2)*
          log (mHp2/MZ2))/(128.0*CW6*(-1.0 + CW2)*(-1.0 + CW2)*mHp4*
          MZ6*pi*pi)+ (aem*aem*mMU*
          mMU*(2.0*mHp8 + mH2_6*(mHp2 + CW2*MZ2) - 
            mH2_2*(3.0*mHp6 + 4.0*CW2*mHp4*MZ2))*log (mH2_2/mHp2)*
          log (mHp2/MZ2))/(128.0*CW6*(-1.0 + CW2)*(-1.0 + CW2)*mHp4*
          MZ6*pi*pi)+ (aem*aem*mH3_2*(MH2 - mH2_2)*mHp2*mMU*mMU*
          log (mH3_2/MZ2)*log (mHp2/MZ2))/(256.0*
          CW2*(-1.0 + CW2)*(-1.0 + CW2)*
          MZ2*(-mHp2 + CW2*MZ2)*(-mHp2 + CW2*MZ2)*(-mHp2 + CW2*MZ2)*
          pi*pi)  - (aem*aem*(MH2 - mH2_2)*mMU*
          mMU*(mH3_2*mHp6 + 9.0*mH2_2*mHp6 - 8.0*mHp8 + 
            2.0*CW2*mH2_2*mHp4*MZ2 + 22.0*CW2*mHp6*MZ2 - 
            10.0*CW4*mH2_2*mHp2*MZ4 - 26*CW4*mHp4*MZ4 + 
            6.0*CW6*mH2_2*MZ6 + 18.0*CW6*mHp2*MZ6 - 6.0*CW8*MZ8 + 
            MH2*(7.0*mHp6 + 2.0*CW2*mHp4*MZ2 - 10.0*CW4*mHp2*MZ4 + 
               6.0*CW6*MZ6))*log(mHp2/MZ2)*log(mHp2/MZ2))/(512.0*
          CW2*(-1.0 + CW2)*(-1.0 + CW2)*mHp4*
          MZ2*(-mHp2 + CW2*MZ2)*(-mHp2 + CW2*MZ2)*(-mHp2 + CW2*MZ2)*
          pi*pi) +  log (mH2_2/
           MZ2)*((aem*aem*mMU*
             mMU*(-6.0*mH2_8 + 
               mH2_4*MZ2*((5.0 + 26*CW2 - 16.0*CW4)*mHp2 - 
                  3.0*CW2*MZ2) + 
               3.0*mH2_6*((3.0 - 12.0*CW2 + 8.0*CW4)*
                   mHp2 + (2.0 + CW2)*MZ2) + 
               mH2_2*mHp2*
                MZ2*((9.0 - 6.0*CW2)*mHp2 + (15 - 59*CW2 + 62*CW4)*
                   MZ2) + 
               4.0*mHp2*
                MZ4*(-5.0*MZ2 + 18.0*CW6*MZ2 + 
                  CW4*(6.0*mHp2 - 28*MZ2) - 
                  6.0*CW2*(mHp2 - 2.0*MZ2))))/(192*
             CW4*(-1.0 + CW2)*(-1.0 + CW2)*mH2_2*mHp2*
             MZ4*(mH2_2 - MZ2)*pi*pi) + (aem*aem*
             mH2_2*(MH2 - 7.0*mH2_2)*mHp2*mMU*mMU*
             log (mHp2/MZ2))/(256.0*CW2*(-1.0 + CW2)*(-1.0 + CW2)*
             MZ2*(-mHp2 + CW2*MZ2)*(-mHp2 + CW2*MZ2)*(-mHp2 + 
               CW2*MZ2)*pi*pi)) + 
       log (MH2/
          MZ2)*(-(aem*aem*mMU*
              mMU*(-8.0*MH10- 128.0*CW4*(-1.0 + CW2)*mHp4*MZ6 + 
                2.0*CW2*MH2*mHp2*
                 MZ4*(8.0*(-5.0 + 4.0*CW2)*mHp2 + 
                   CW2*(33 - 32*CW2)*MZ2) + 
                MH6*MZ2*((10.0 + 11.0*CW2)*mHp2 - 
                   4.0*CW2*(9.0 + 4.0*CW2)*MZ2) - 
                4.0*MH8*(mHp2 - (2.0 + 9.0*CW2)*MZ2) + 
                MH4*MZ2*((12.0 - 8.0*CW2)*mHp4 + 
                   5.0*CW2*(-11.0 + 6.0*CW2)*mHp2*MZ2 + 
                   16.0*CW4*MZ4)))/(256.0*CW4*(-1.0 + CW2)*(-1.0 + CW2)*
             MH2*mHp2*MZ4*(MH2 - MZ2)*(MH2 - 4.0*CW2*MZ2)*pi*
             pi) + (aem*aem*MH2*(5.0*MH2 + mH2_2)*mHp2*mMU*mMU*
             log (mHp2/MZ2))/(256.0*CW2*(-1.0 + CW2)*(-1.0 + CW2)*
             MZ2*(-mHp2 + CW2*MZ2)*(-mHp2 + CW2*MZ2)*(-mHp2 + 
               CW2*MZ2)*pi*pi)) + 
       log (CW2)*((aem*aem*mMU*
             mMU*(24.0*(1.0 - 4.0*CW2)*MH6*mH2_2*mHp2 + 
               12.0*CW2*(-1.0 + 4.0*CW2)*MH4*mH2_2*
                MZ2*(8.0*mHp2 + CW2*MZ2) + 2.0*CW2*
                MZ2*(-48.0*(-1.0 + 4.0*CW2)*mH2_6*mHp2 + 
                  24.0*CW4*(-1.0 + 4.0*CW2)*mH2_4*MZ4 + 
                  CW2*(133 - 520*CW2 + 240*CW4)*mH2_2*mHp2*
                   MZ2*(-mHp2 + CW2*MZ2) + 
                  64.0*CW4*(10.0 - 49*CW2 + 36*CW4)*mHp2*
                   MZ4*(-mHp2 + CW2*MZ2)) + 
               MH2*(24.0*(-1.0 + 4.0*CW2)*mH2_6*mHp2 + 
                  12.0*CW4*(1.0 - 4.0*CW2)*mH2_4*MZ4 - 
                  32*CW4*(10.0 - 49*CW2 + 36*CW4)*mHp2*
                   MZ4*(-mHp2 + CW2*MZ2) + 
                  CW2*mH2_2*
                   MZ2*((71.0 - 278.0*CW2 + 120*CW4)*mHp4 + 
                    CW2*(-71.0 + 278.0*CW2 - 120*CW4)*mHp2*MZ2 + 
                    48.0*CW4*(1.0 - 4.0*CW2)*MZ4))))/(768*
             CW4*(-1.0 + CW2)*(-1.0 + CW2)*(-1.0 + 4.0*CW2)*mH2_2*
             mHp2*MZ4*(-mHp2 + CW2*MZ2)*(-MH2 + 4.0*CW2*MZ2)*pi*
             pi) + (aem*aem*mMU*
             mMU*(2.0*mHp8 + MH6*(mHp2 + CW2*MZ2) - 
               MH2*(3.0*mHp6 + 4.0*CW2*mHp4*MZ2))*
             log (MH2/mHp2))/(128.0*CW6*(-1.0 + CW2)*(-1.0 + CW2)*
             mHp4*MZ6*pi*pi) - (aem*aem*mMU*
             mMU*(2.0*mHp8 + mH2_6*(mHp2 + CW2*MZ2) - 
               mH2_2*(3.0*mHp6 + 4.0*CW2*mHp4*MZ2))*
             log (mH2_2/mHp2))/(128.0*CW6*(-1.0 + CW2)*(-1.0 + CW2)*
             mHp4*MZ6*pi*pi) + (aem*aem*mMU*
             mMU*(-2.0*MH6 + 3.0*CW2*mHp4*MZ2 - 4.0*CW4*mHp2*MZ4 + 
               CW6*MZ6 + MH4*(mHp2 + 3.0*CW2*MZ2) + 
               MH2*(4.0*mHp4 + 3.0*CW2*mHp2*MZ2 - 3.0*CW4*MZ4))*
             log (MH2/MZ2))/(128.0*CW2*(-1.0 + CW2)*(-1.0 + CW2)*mHp4*
             MZ2*(-mHp2 + CW2*MZ2)*pi*pi) - (aem*aem*mMU*
             mMU*(-2.0*mH2_6 + 3.0*CW2*mHp4*MZ2 - 4.0*CW4*mHp2*MZ4 + 
               CW6*MZ6 + mH2_4*(mHp2 + 3.0*CW2*MZ2) + 
               mH2_2*(4.0*mHp4 + 3.0*CW2*mHp2*MZ2 - 3.0*CW4*MZ4))*
             log (mH2_2/MZ2))/(128.0*CW2*(-1.0 + CW2)*(-1.0 + CW2)*
             mHp4*MZ2*(-mHp2 + CW2*MZ2)*pi*pi) + (aem*aem*mMU*
             mMU*(MH4*(mHp2 + 3.0*CW2*MZ2) + 
               MH2*(-4.0*mHp4 + 3.0*CW2*mHp2*MZ2 - 3.0*CW4*MZ4) - 
               mH2_2*(-4.0*mHp4 + 3.0*CW2*mHp2*MZ2 - 3.0*CW4*MZ4 + 
                  mH2_2*(mHp2 + 3.0*CW2*MZ2)))*
             log (mHp2/MZ2))/(128.0*CW2*(-1.0 + CW2)*(-1.0 + CW2)*
             mHp4*MZ2*(-mHp2 + CW2*MZ2)*pi*pi)) + (aem*aem*mMU*
          mMU*(2.0*mHp8 + MH6*(mHp2 + CW2*MZ2) - 
            MH2*(3.0*mHp6 + 4.0*CW2*mHp4*MZ2))*
          PolyLog.Li2(1.0 - MH2/mHp2))/(64.0*
          CW6*(-1.0 + CW2)*(-1.0 + CW2)*mHp4*MZ6*pi*pi)  - (aem*aem*
          mMU*mMU*(2.0*mHp8 + mH2_6*(mHp2 + CW2*MZ2) - 
            mH2_2*(3.0*mHp6 + 4.0*CW2*mHp4*MZ2))*
          PolyLog.Li2(1.0 - mH2_2/mHp2))/(64.0*
          CW6*(-1.0 + CW2)*(-1.0 + CW2)*mHp4*MZ6*pi*pi)+ (aem*aem*
          mMU*mMU*(2.0*MH6*(mHp2 + CW2*MZ2) - 
            CW4*MZ4*(mHp4 - 6.0*CW2*mHp2*MZ2 + 2.0*CW4*MZ4) - 
            2.0*MH4*(mHp4 + 4.0*CW2*mHp2*MZ2 + 3.0*CW4*MZ4) + 
            MH2*(4.0*CW2*mHp4*MZ2 + 6.0*CW6*MZ6))*
          PolyLog.Li2(1.0 - MH2/(CW2*MZ2)))/(128.0*
          CW6*(-1.0 + CW2)*(-1.0 + CW2)*mHp4*MZ6*pi*pi)+ (aem*aem*
          mMU*mMU*(12.0*(1.0 - 3.0*CW2 + 2.0*CW4)*mH2_10 - 
            48.0*(1.0 - 3.0*CW2 + 2.0*CW4)*mH2_8*
             MZ2 + (17.0 - 48.0*CW2 + 32.0*CW4)*mH2_6*MZ4 - 
            3.0*(5.0 - 12.0*CW2 + 8.0*CW4)*mH2_2*MZ8 - 
            4.0*(5.0 - 12.0*CW2 + 8.0*CW4)*MZ10)*
          PolyLog.Li2(  1.0 - mH2_2/MZ2))/(192.0*
          CW4*(-1.0 + CW2)*(-1.0 + CW2)*mH2_4*MZ8*pi*pi) + (aem*aem*
          mMU*mMU*(-6.0*CW6*mH2_2*mHp4*MZ6 - 8.0*CW8*mHp4*MZ8 - 
            3.0*mH2_10*(mHp2 + CW2*MZ2) + 
            3.0*mH2_8*(mHp4 + 4.0*CW2*mHp2*MZ2 + 3.0*CW4*MZ4) + 
            mH2_6*(2.0*CW2*mHp4*MZ2 - 9.0*CW6*MZ6) + 
            3.0*mH2_4*(-3.0*CW6*mHp2*MZ6 + CW8*MZ8))*
          PolyLog.Li2( 1.0 - mH2_2/(CW2*MZ2)))/(192*
          CW6*(-1.0 + CW2)*(-1.0 + CW2)*mH2_4*mHp4*MZ6*pi*pi)   - (aem*
          aem*(-3.0 + 4.0*CW2)*mMU*
          mMU*((1.0 - 5.0*CW2 + 10.0*CW4)*
             mH2_2 + (-7.0 + 61*CW2 - 162*CW4 + 96.0*CW6)*MZ2)*
          TF (CW, CW, 1.0))/(128.0*
          CW2*(1.0 - 4.0*CW2)*(1.0 - 4.0*CW2)*(-1.0 + CW2)*(-1.0 + 
            CW2)*MZ2*(-mH2_2 + MZ2)*pi*pi) + (aem*
          aem*(MH2 + 2.0*mHp2)*mMU*
          mMU*(3.0*MH2*mHp8 - 8.0*CW6*mHp6*MZ4 + 
            2.0*CW4*mHp6*MZ2*(MH2 + 4.0*(mHp2 + MZ2)) + 
            CW2*(MH10 - 2.0*MH4*mHp4*MZ2 - 8.0*mHp8*MZ2 - 
               MH8*(4.0*mHp2 + MZ2) + 
               2.0*MH6*(mHp4 + 2.0*mHp2*MZ2) - 
               MH2*(2.0*mHp8 + 3.0*mHp6*MZ2)))*
          TF (MH, mHp, mHp))/(64.0*CW4*(-1.0 + CW2)*(-1.0 + CW2)*MH4*
          mHp4*(MH2 - 4.0*mHp2)*MZ2*(MH2 - MZ2)*(-mHp2 + CW2*MZ2)*pi*
          pi) + (aem*aem*mMU*
          mMU*(MH8*(-mHp4 + CW4*MZ4) + 
            MH6*(mHp6 + CW2*mHp4*MZ2 - 2.0*CW4*mHp2*MZ4 - 
               4.0*CW6*MZ6) + (mHp2 - CW2*MZ2)*(mHp2 - 
               CW2*MZ2)*(mHp2 - CW2*MZ2)*(2.0*mHp6 + 
               2.0*CW2*mHp4*MZ2 + 2.0*CW4*mHp2*MZ4 - CW6*MZ6) + 
            MH4*(3.0*mHp8 + 3.0*CW2*mHp6*MZ2 + CW4*mHp4*MZ4 - 
               CW6*mHp2*MZ6 + 6.0*CW8*MZ8) + 
            MH2*(-5.0*mHp10 + CW4*mHp6*MZ4 + 8.0*CW8*mHp2*MZ8 - 
               4.0*CW10*MZ10))*TF (MH, mHp, CW*MZ))/(64.0*
          CW6*(-1.0 + CW2)*(-1.0 + CW2)*mHp4*
          MZ6*(-mHp2 + 
            CW2*MZ2)*(MH4 + (mHp2 - CW2*MZ2)*(mHp2 - CW2*MZ2) - 
            2.0*MH2*(mHp2 + CW2*MZ2))*pi*pi) + (aem*aem*mMU*
          mMU*(-2.0*MH10 + 2.0*CW8*MZ8*(-mHp2 + CW2*MZ2) + 
            2.0*MH8*(mHp2 + 9.0*CW2*MZ2) - 
            4.0*MH6*(4.0*CW2*mHp2*MZ2 + 11.0*CW4*MZ4) + 
            MH4*(37*CW4*mHp2*MZ4 + 3.0*CW6*MZ6) + 
            2.0*MH2*(-9.0*CW6*mHp2*MZ6 + 25*CW8*MZ8))*
          TF (MH, CW*MZ, CW*MZ))/(128.0*
          CW6*(-1.0 + CW2)*(-1.0 + CW2)*MH2*
          MZ6*(MH2 - 4.0*CW2*MZ2)*(MH2 - 4.0*CW2*MZ2)*(-mHp2 + 
            CW2*MZ2)*pi*pi) + (aem*aem*(mH2_2 + 2.0*mHp2)*mMU*
          mMU*(-3.0*mH2_2*mHp8 + 8.0*CW6*mHp6*MZ4 - 
            2.0*CW4*mHp6*MZ2*(mH2_2 + 4.0*(mHp2 + MZ2)) + 
            
            CW2*(-mH2_10 + 2.0*mH2_4*mHp4*MZ2 + 8.0*mHp8*MZ2 + 
               mH2_8*(4.0*mHp2 + MZ2) - 
               2.0*mH2_6*(mHp4 + 2.0*mHp2*MZ2) + 
               mH2_2*(2.0*mHp8 + 3.0*mHp6*MZ2)))*
          TF (mH2, mHp, mHp))/(64.0*CW4*(-1.0 + CW2)*(-1.0 + CW2)*
          mH2_4*mHp4*(mH2_2 - 4.0*mHp2)*
          MZ2*(mH2_2 - MZ2)*(-mHp2 + CW2*MZ2)*pi*pi) - (aem*aem*mMU*
          mMU*(mH2_8*(-mHp4 + CW4*MZ4) + 
            mH2_6*(mHp6 + CW2*mHp4*MZ2 - 2.0*CW4*mHp2*MZ4 - 
               4.0*CW6*MZ6) + (mHp2 - CW2*MZ2)*(mHp2 - 
               CW2*MZ2)*(mHp2 - CW2*MZ2)*(2.0*mHp6 + 
               2.0*CW2*mHp4*MZ2 + 2.0*CW4*mHp2*MZ4 - CW6*MZ6) + 
            mH2_4*(3.0*mHp8 + 3.0*CW2*mHp6*MZ2 + CW4*mHp4*MZ4 - 
               CW6*mHp2*MZ6 + 6.0*CW8*MZ8) + 
            mH2_2*(-5.0*mHp10 + CW4*mHp6*MZ4 + 8.0*CW8*mHp2*MZ8 - 
               4.0*CW10*MZ10))*TF (mH2, mHp, CW*MZ))/(64.0*
          CW6*(-1.0 + CW2)*(-1.0 + CW2)*mHp4*
          MZ6*(-mHp2 + 
            CW2*MZ2)*(mH2_4 + (mHp2 - CW2*MZ2)*(mHp2 - CW2*MZ2) - 
            2.0*mH2_2*(mHp2 + CW2*MZ2))*pi*pi) + (aem*aem*mMU*
          mMU*(6.0*(1.0 - 3.0*CW2 + 2.0*CW4)*mH2_6 - 
            12.0*(1.0 - 3.0*CW2 + 2.0*CW4)*mH2_4*
             MZ2 + (5.0 - 12.0*CW2 + 8.0*CW4)*mH2_2*MZ4 + 
            2.0*(5.0 - 12.0*CW2 + 8.0*CW4)*MZ6)*
          TF (mH2, MZ, MZ))/(192*CW4*(-1.0 + CW2)*(-1.0 + CW2)*mH2_2*
          MZ8*pi*pi) + (aem*aem*mMU*
          mMU*(3.0*mH2_12 - 
            2.0*CW6*(-95 + 78*CW2)*mH2_2*MZ8*(-mHp2 + CW2*MZ2) + 
            144*CW8*(-1.0 + CW2)*MZ10*(-mHp2 + CW2*MZ2) + 
            2.0*CW2*mH2_6*
             MZ4*((-2.0 + 5.0*CW2)*mHp2 + CW2*(-7.0 + CW2)*MZ2) - 
            3.0*mH2_10*(mHp2 + (1.0 + 5.0*CW2)*MZ2) + 
            CW4*mH2_4*
             MZ6*((17.0 + 10.0*CW2)*mHp2 - CW2*(29 + 10.0*CW2)*MZ2) + 
            mH2_8*MZ2*((3.0 + 4.0*CW2)*mHp2 + 
               CW2*(15 + 14*CW2)*MZ2))*TF (mH2, CW*MZ, CW*MZ))/(192*
          CW6*(-1.0 + CW2)*(-1.0 + CW2)*mH2_4*
          MZ6*(mH2_2 - MZ2)*(mH2_2 - 4.0*CW2*MZ2)*(-mHp2 + CW2*MZ2)*
          pi*pi) + (aem*aem*(3.0 - 10.0*CW2 + 8.0*CW4)*(MH2 - mH2_2)*
          mHp2*mMU*mMU*(2.0*mHp2 + MZ2)*TF(mHp, mHp, MZ))/(64.0*
          CW4*(-1.0 + CW2)*(-1.0 + CW2)*
          MZ4*(MH2 - MZ2)*(-mH2_2 + MZ2)*(-4.0*mHp2 + MZ2)*pi*pi)));
       

        
      
gslpp::complex aYuk5 = Lambda5*(-(aem*mMU*mMU)/(8.0*MH2*pi*pi*pi) + (aem*mMU*
       mMU*((-3.0 + 2.0*CW2)*MH2 - 8.0*CW2*(-1.0 + CW2)*MZ2)*
       log (MH2/MZ2))/(128.0*CW2*(-1.0 + CW2)*MH2*(MH2 - MZ2)*pi*pi*
       pi) - (aem*mMU*mMU*log(mHp2/MZ2))/(16.0*MH2*pi*pi*pi) + (aem*
       mHp2*mMU*mMU*((-3.0 + 2.0*CW2)*MH2 - 8.0*CW2*(-1.0 + CW2)*MZ2)*
       TF (MH, mHp, mHp))/(64.0*CW2*(-1.0 + CW2)*
       MH4*(MH2 - 4.0*mHp2)*(MH2 - MZ2)*pi*pi*
       pi) - (aem*(3.0 - 10.0*CW2 + 8.0*CW4)*mHp2*mMU*mMU*
       TF (mHp, mHp, MZ))/(64.0*CW2*(-1.0 + CW2)*
       MZ2*(-MH2 + MZ2)*(-4.0*mHp2 + MZ2)*pi*pi*pi) +  sl*((aem*mMU*
          mMU*(-mH2_4 + mH2_2*mHp2 + 
            8.0*CW2*(-1.0 + CW2)*mHp2*MZ2)*(-1.0 + tanb*tanb))/(128.0*
          CW2*(-1.0 + CW2)*mH2_2*mHp2*MZ2*pi*pi*pi*tanb) - (aem*
          mH3_2*mHp2*mMU*mMU*(-1.0 + tanb*tanb)*
          log(mH3_2/mHp2)*log(mH3_2/mHp2))/(1024*(-1.0 + 
            CW2)*(-mHp2 + CW2*MZ2)*(-mHp2 + CW2*MZ2)*(-mHp2 + 
            CW2*MZ2)*pi*pi*pi*tanb) + (aem*mH2_2*mMU*
          mMU*(-7.0*mHp6 + 8.0*CW2*mHp4*MZ2 - 4.0*CW4*mHp2*MZ4 + 
            2.0*mH2_2*(mHp2 - CW2*MZ2)*(mHp2 - CW2*MZ2))*(-1.0 + 
            tanb*tanb)*
          log(mH2_2/mHp2)*
          log(mH2_2/mHp2))/(1024*(-1.0 + CW2)*
          mHp4*(-mHp2 + CW2*MZ2)*(-mHp2 + CW2*MZ2)*(-mHp2 + CW2*MZ2)*
          pi*pi*pi*tanb) + (aem*mH3_2*mHp2*mMU*
          mMU*(-1.0 + tanb*tanb)*
          log (mH3_2/MZ2)* log (mH3_2/MZ2))/(1024*(-1.0 + CW2)*(-mHp2 +
             CW2*MZ2)*(-mHp2 + CW2*MZ2)*(-mHp2 + CW2*MZ2)*pi*pi*pi*
          tanb) - (aem*mH2_2*mMU*
          mMU*(-7.0*mHp6 + 8.0*CW2*mHp4*MZ2 - 4.0*CW4*mHp2*MZ4 + 
            2.0*mH2_2*(mHp2 - CW2*MZ2)*(mHp2 - CW2*MZ2))*(-1.0 + 
            tanb*tanb)*
          log (mH2_2/MZ2)*
          log (mH2_2/MZ2))/(1024*(-1.0 + CW2)*
          mHp4*(-mHp2 + CW2*MZ2)*(-mHp2 + CW2*MZ2)*(-mHp2 + CW2*MZ2)*
          pi*pi*pi*tanb) + (aem*mMU*
          mMU*(2.0*CW2*mH2_4*MZ2 + 
            8.0*CW2*(-1.0 + CW2)*mHp2*MZ2*(-mHp2 + CW2*MZ2) + 
            mH2_2*(-2.0*mHp4 + CW2*mHp2*MZ2))*(-1.0 + tanb*tanb)*
          log (mHp2/MZ2))/(256.0*CW2*(-1.0 + CW2)*mH2_2*mHp2*
          MZ2*(-mHp2 + CW2*MZ2)*pi*pi*pi*tanb) + (aem*mMU*
          mMU*(mHp6 + mH2_4*(mHp2 + CW2*MZ2) - 
            2.0*mH2_2*(mHp4 + CW2*mHp2*MZ2))*(-1.0 + tanb*tanb)*
          log (mH2_2/mHp2)*log (mHp2/MZ2))/(256.0*CW4*(-1.0 + CW2)*
          mHp4*MZ4*pi*pi*pi*tanb) - (aem*mH3_2*mHp2*mMU*
          mMU*(-1.0 + tanb*tanb)*log (mH3_2/MZ2)*
          log (mHp2/MZ2))/(512.0*(-1.0 + CW2)*(-mHp2 + 
            CW2*MZ2)*(-mHp2 + CW2*MZ2)*(-mHp2 + CW2*MZ2)*pi*pi*pi*
          tanb) + (aem*mMU*
          mMU*(mH3_2*mHp6 + 
            2.0*mH2_4*(mHp2 - CW2*MZ2)*(mHp2 - CW2*MZ2) - 
            mH2_2*(mHp6 - 8.0*CW2*mHp4*MZ2 + 
               4.0*CW4*mHp2*MZ4))*(-1.0 + tanb*tanb)*
          log (mHp2/MZ2)*
          log (mHp2/MZ2))/(1024*(-1.0 + CW2)*
          mHp4*(-mHp2 + CW2*MZ2)*(-mHp2 + CW2*MZ2)*(-mHp2 + CW2*MZ2)*
          pi*pi*pi*tanb) + 
       log (mH2_2/
           MZ2)*((aem*mMU*
             mMU*(-2.0*mH2_6 + 
               2.0*mH2_4*MZ2 + (3.0 - 2.0*CW2)*mH2_2*mHp2*MZ2 + 
               8.0*CW2*(-1.0 + CW2)*mHp2*MZ4)*(-1.0 + 
               tanb*tanb))/(256.0*CW2*(-1.0 + CW2)*mH2_2*mHp2*
             MZ2*(mH2_2 - MZ2)*pi*pi*pi*tanb) - (3.0*aem*mH2_2*mHp2*
             mMU*mMU*(-1.0 + tanb*tanb)*
             log (mHp2/MZ2))/(512.0*(-1.0 + CW2)*(-mHp2 + 
               CW2*MZ2)*(-mHp2 + CW2*MZ2)*(-mHp2 + CW2*MZ2)*pi*pi*pi*
             tanb)) + 
       log (CW2)*( -(aem*mMU*
              mMU*(2.0*mH2_2 - 2.0*mHp2 + CW2*MZ2)*(-1.0 + 
                tanb*tanb))/(256.0*CW2*(-1.0 + CW2)*
             MZ2*(-mHp2 + CW2*MZ2)*pi*pi*pi*tanb)
         - (aem*mMU*
             mMU*(mHp6 + mH2_4*(mHp2 + CW2*MZ2) - 
               2.0*mH2_2*(mHp4 + CW2*mHp2*MZ2))*(-1.0 + tanb*tanb)*
             log (mH2_2/mHp2))/(256.0*CW4*(-1.0 + CW2)*mHp4*MZ4*pi*pi*
             pi*tanb)   + (aem*mH2_2*(mH2_2 - 2.0*mHp2)*mMU*
             mMU*(-1.0 + tanb*tanb)*
             log (mH2_2/MZ2))/(256.0*(-1.0 + CW2)*
             mHp4*(-mHp2 + CW2*MZ2)*pi*pi*pi*tanb) - (aem*
             mH2_2*(mH2_2 - 2.0*mHp2)*mMU*mMU*(-1.0 + tanb*tanb)*
             log (mHp2/MZ2))/(256.0*(-1.0 + CW2)*mHp4*(-mHp2 + CW2*MZ2)*
             pi*pi*pi*tanb)  ) -(aem*mMU*
    mMU*(mHp6 + mH2_4*(mHp2 + CW2*MZ2) - 
      2.0*mH2_2*(mHp4 + CW2*mHp2*MZ2))*(-1.0 + tanb*tanb)*
    PolyLog.Li2 (1.0 - mH2_2/mHp2))/(128.0*CW4*(-1.0 + CW2)*mHp4*MZ4*
   pi*pi*pi*tanb)  + (aem*mMU*
          mMU*(-3.0*mH2_2*mHp8 + 8.0*CW6*mHp6*MZ4 - 
            2.0*CW4*mHp6*MZ2*(mH2_2 + 4.0*(mHp2 + MZ2)) + 
            CW2*(-mH2_10 + 2.0*mH2_4*mHp4*MZ2 + 8.0*mHp8*MZ2 + 
               mH2_8*(4.0*mHp2 + MZ2) - 
               2.0*mH2_6*(mHp4 + 2.0*mHp2*MZ2) + 
               mH2_2*(2.0*mHp8 + 3.0*mHp6*MZ2)))*(-1.0 + tanb*tanb)*
          TF (mH2, mHp, mHp))/(128.0*CW2*(-1.0 + CW2)*mH2_4*
          mHp4*(mH2_2 - 4.0*mHp2)*(mH2_2 - MZ2)*(-mHp2 + CW2*MZ2)*
          pi*pi*pi*tanb)  - (aem*mMU*
          mMU*(-mH2_6 + 
            mH2_4*(3.0*mHp2 + CW2*MZ2) + (mHp3 - CW2*mHp*MZ2)*(mHp3 -
                CW2*mHp*MZ2) + 
            mH2_2*(-3.0*mHp4 + CW2*mHp2*MZ2))*(-1.0 + tanb*tanb)*
          TF (mH2, mHp, CW*MZ))/(128.0*CW4*(-1.0 + CW2)*
          MZ4*(-mHp2 + 
            CW2*MZ2)*(mH2_4 + (mHp2 - CW2*MZ2)*(mHp2 - CW2*MZ2) - 
            2.0*mH2_2*(mHp2 + CW2*MZ2))*pi*pi*pi*
          tanb)+(aem*(3.0 - 10.0*CW2 + 8.0*CW4)*mHp2*mMU*
          mMU*(-1.0 + tanb*tanb)*TF (mHp, mHp, MZ))/(128.0*
          CW2*(-1.0 + CW2)*MZ2*(-mH2_2 + MZ2)*(-4.0*mHp2 + MZ2)*pi*
          pi*pi*tanb)) + eta*(-(aem*(MH2 - mH2_2)*mMU*mMU*(-1.0 + tanb*tanb))/(16.0*MH2*
          mH2_2*pi*pi*pi*tanb) + (aem*mMU*
          mMU*((3.0 - 2.0*CW2)*MH2 + 8.0*CW2*(-1.0 + CW2)*MZ2)*(-1.0 +
             tanb*tanb)*log (MH2/MZ2))/(256.0*CW2*(-1.0 + CW2)*
          MH2*(MH2 - MZ2)*pi*pi*pi*tanb) + (aem*mMU*
          mMU*((-3.0 + 2.0*CW2)*mH2_2 - 
            8.0*CW2*(-1.0 + CW2)*MZ2)*(-1.0 + tanb*tanb)*
          log (mH2_2/MZ2))/(256.0*CW2*(-1.0 + CW2)*
          mH2_2*(mH2_2 - MZ2)*pi*pi*pi*tanb) - (aem*(MH2 - mH2_2)*
          mMU*mMU*(-1.0 + tanb*tanb)*log (mHp2/MZ2))/(32*MH2*mH2_2*
          pi*pi*pi*tanb) + (aem*mHp2*mMU*
          mMU*((3.0 - 2.0*CW2)*MH2 + 8.0*CW2*(-1.0 + CW2)*MZ2)*(-1.0 +
             tanb*tanb)*TF (MH, mHp, mHp))/(128.0*CW2*(-1.0 + CW2)*
          MH4*(MH2 - 4.0*mHp2)*(MH2 - MZ2)*pi*pi*pi*tanb)+ (aem*mHp2*
          mMU*mMU*((-3.0 + 2.0*CW2)*mH2_2 - 
            8.0*CW2*(-1.0 + CW2)*MZ2)*(-1.0 + tanb*tanb)*
          TF (mH2, mHp, mHp))/(128.0*CW2*(-1.0 + CW2)*
          mH2_4*(mH2_2 - 4.0*mHp2)*(mH2_2 - MZ2)*pi*pi*pi*
          tanb) + (aem*(3.0 - 10.0*CW2 + 8.0*CW4)*(MH2 - mH2_2)*mHp2*
          mMU*mMU*(-1.0 + tanb*tanb)*TF (mHp, mHp, MZ))/(128.0*
          CW2*(-1.0 + CW2)*
          MZ2*(-MH2 + MZ2)*(-mH2_2 + MZ2)*(-4.0*mHp2 + MZ2)*pi*pi*pi*
          tanb) + sl*((aem*(MH2 - mH2_2)*mMU*
             mMU*(MH2*mH2_2 + 8.0*CW2*(-1.0 + CW2)*mHp2*MZ2))/(64.0*
             CW2*(-1.0 + CW2)*MH2*mH2_2*mHp2*MZ2*pi*pi*pi) - (aem*
             MH2*mMU*
             mMU*(-3.0*mHp6 + 4.0*CW2*mHp4*MZ2 - 2.0*CW4*mHp2*MZ4 + 
               MH2*(mHp2 - CW2*MZ2)*(mHp2 - CW2*MZ2))*
             log (MH2/mHp2)*
             log (MH2/mHp2))/(256.0*(-1.0 + CW2)*
             mHp4*(-mHp2 + CW2*MZ2)*(-mHp2 + CW2*MZ2)*(-mHp2 + 
               CW2*MZ2)*pi*pi*pi) + (aem*mH2_2*mMU*
             mMU*(-3.0*mHp6 + 4.0*CW2*mHp4*MZ2 - 2.0*CW4*mHp2*MZ4 + 
               mH2_2*(mHp2 - CW2*MZ2)*(mHp2 - CW2*MZ2))*
             log (mH2_2/mHp2)*
             log (mH2_2/mHp2))/(256.0*(-1.0 + CW2)*
             mHp4*(-mHp2 + CW2*MZ2)*(-mHp2 + CW2*MZ2)*(-mHp2 + 
               CW2*MZ2)*pi*pi*pi) + (aem*MH2*mMU*
             mMU*(-3.0*mHp6 + 4.0*CW2*mHp4*MZ2 - 2.0*CW4*mHp2*MZ4 + 
               MH2*(mHp2 - CW2*MZ2)*(mHp2 - CW2*MZ2))*
             log (MH2/MZ2)*
             log (MH2/MZ2))/(256.0*(-1.0 + CW2)*
             mHp4*(-mHp2 + CW2*MZ2)*(-mHp2 + CW2*MZ2)*(-mHp2 + 
               CW2*MZ2)*pi*pi*pi) - (aem*mH2_2*mMU*
             mMU*(-3.0*mHp6 + 4.0*CW2*mHp4*MZ2 - 2.0*CW4*mHp2*MZ4 + 
               mH2_2*(mHp2 - CW2*MZ2)*(mHp2 - CW2*MZ2))*
             log (mH2_2/MZ2)*
             log (mH2_2/MZ2))/(256.0*(-1.0 + CW2)*
             
             mHp4*(-mHp2 + CW2*MZ2)*(-mHp2 + CW2*MZ2)*(-mHp2 + 
               CW2*MZ2)*pi*pi*pi)  - (aem*(MH2 - mH2_2)*mMU*
             mMU*(MH2*mH2_2 + 
               4.0*(-1.0 + CW2)*mHp2*(mHp2 - CW2*MZ2))*
             log (mHp2/MZ2))/(64.0*(-1.0 + CW2)*MH2*mH2_2*
             mHp2*(-mHp2 + CW2*MZ2)*pi*pi*pi) - (aem*mMU*
             mMU*(mHp6 + MH4*(mHp2 + CW2*MZ2) - 
               2.0*MH2*(mHp4 + CW2*mHp2*MZ2))*log (MH2/mHp2)*
             log (mHp2/MZ2))/(128.0*CW4*(-1.0 + CW2)*mHp4*MZ4*pi*pi*
             pi) + (aem*mMU*
             mMU*(mHp6 + mH2_4*(mHp2 + CW2*MZ2) - 
               2.0*mH2_2*(mHp4 + CW2*mHp2*MZ2))*log(mH2_2/mHp2)*
             log (mHp2/MZ2))/(128.0*CW4*(-1.0 + CW2)*mHp4*MZ4*pi*pi*
             pi) + (aem*mMU*
             mMU*(-(MH4*(mHp2 - CW2*MZ2)*(mHp2 - CW2*MZ2)) + 
               MH2*(mHp6 - 4.0*CW2*mHp4*MZ2 + 2.0*CW4*mHp2*MZ4) + 
               mH2_2*(mH2_2*(mHp2 - CW2*MZ2)*(mHp2 - CW2*MZ2) - 
                  mHp2*(mHp4 - 4.0*CW2*mHp2*MZ2 + 2.0*CW4*MZ4)))*
             log(mHp2/MZ2)*
             log(mHp2/MZ2))/(256.0*(-1.0 + CW2)*
             mHp4*(-mHp2 + CW2*MZ2)*(-mHp2 + CW2*MZ2)*(-mHp2 + 
               CW2*MZ2)*pi*pi*pi)   + 
          log (MH2/
             MZ2)*((aem*mMU*
                mMU*(2.0*MH6 - 
                  2.0*MH4*MZ2 + (-3.0 + 2.0*CW2)*MH2*mHp2*MZ2 - 
                  8.0*CW2*(-1.0 + CW2)*mHp2*MZ4))/(128.0*
                CW2*(-1.0 + CW2)*MH2*mHp2*MZ2*(MH2 - MZ2)*pi*pi*
                pi) + (aem*MH2*mHp2*mMU*mMU*
                log (mHp2/MZ2))/(128.0*(-1.0 + CW2)*(-mHp2 + 
                  CW2*MZ2)*(-mHp2 + CW2*MZ2)*(-mHp2 + CW2*MZ2)*pi*pi*
                pi))  + log (mH2_2/
              MZ2)*((aem*mMU*
                mMU*(-2.0*mH2_6 + 
                  2.0*mH2_4*MZ2 + (3.0 - 2.0*CW2)*mH2_2*mHp2*MZ2 + 
                  8.0*CW2*(-1.0 + CW2)*mHp2*MZ4))/(128.0*
                CW2*(-1.0 + CW2)*mH2_2*mHp2*MZ2*(mH2_2 - MZ2)*pi*pi*
                pi) - (aem*mH2_2*mHp2*mMU*mMU*
                log (mHp2/MZ2))/(128.0*(-1.0 + CW2)*(-mHp2 + 
                  CW2*MZ2)*(-mHp2 + CW2*MZ2)*(-mHp2 + CW2*MZ2)*pi*pi*
                pi))+ log (CW2)*((aem*(MH2 - mH2_2)*mMU*mMU)/(64.0*
                CW2*(-1.0 + CW2)*MZ2*(-mHp2 + CW2*MZ2)*pi*pi*
                pi) + (aem*mMU*
                mMU*(mHp6 + MH4*(mHp2 + CW2*MZ2) - 
                  2.0*MH2*(mHp4 + CW2*mHp2*MZ2))*
                log (MH2/mHp2))/(128.0*CW4*(-1.0 + CW2)*mHp4*MZ4*pi*
                pi*pi) - (aem*mMU*
                mMU*(mHp6 + mH2_4*(mHp2 + CW2*MZ2) - 
                  2.0*mH2_2*(mHp4 + CW2*mHp2*MZ2))*
                log (mH2_2/mHp2))/(128.0*CW4*(-1.0 + CW2)*mHp4*MZ4*
                pi*pi*pi) - (aem*MH2*(MH2 - 2.0*mHp2)*mMU*mMU*
                log (MH2/MZ2))/(128.0*(-1.0 + CW2)*
                mHp4*(-mHp2 + CW2*MZ2)*pi*pi*pi) + (aem*
                mH2_2*(mH2_2 - 2.0*mHp2)*mMU*mMU*
                log (mH2_2/MZ2))/(128.0*(-1.0 + CW2)*
                mHp4*(-mHp2 + CW2*MZ2)*pi*pi*
                pi) + (aem*(MH2 - mH2_2)*(MH2 + mH2_2 - 2.0*mHp2)*
                mMU*mMU*log (mHp2/MZ2))/(128.0*(-1.0 + CW2)*
                mHp4*(-mHp2 + CW2*MZ2)*pi*pi*pi)) + (aem*mMU*
             mMU*(mHp6 + MH4*(mHp2 + CW2*MZ2) - 
               2.0*MH2*(mHp4 + CW2*mHp2*MZ2))*
             PolyLog.Li2 (1.0 - MH2/mHp2))/(64.0*CW4*(-1.0 + CW2)*
             mHp4*MZ4*pi*pi*pi) - (aem*mMU*
             mMU*(mHp6 + mH2_4*(mHp2 + CW2*MZ2) - 
               2.0*mH2_2*(mHp4 + CW2*mHp2*MZ2))*
             PolyLog.Li2(1.0 - mH2_2/mHp2))/(64.0*CW4*(-1.0 + CW2)*
             mHp4*MZ4*pi*pi*pi) + (aem*mMU*
             mMU*(3.0*MH2*mHp8 - 8.0*CW6*mHp6*MZ4 + 
               2.0*CW4*mHp6*MZ2*(MH2 + 4.0*(mHp2 + MZ2)) + 
               CW2*(MH10 - 2.0*MH4*mHp4*MZ2 - 8.0*mHp8*MZ2 - 
                  MH8*(4.0*mHp2 + MZ2) + 
                  2.0*MH6*(mHp4 + 2.0*mHp2*MZ2) - 
                  MH2*(2.0*mHp8 + 3.0*mHp6*MZ2)))*
             TF (MH, mHp, mHp))/(64.0*CW2*(-1.0 + CW2)*MH4*
             mHp4*(MH2 - 4.0*mHp2)*(MH2 - MZ2)*(-mHp2 + CW2*MZ2)*pi*
             pi*pi) + (aem*mMU*
             mMU*(-MH6 + 
               MH4*(3.0*mHp2 + CW2*MZ2) + (mHp3 - CW2*mHp*MZ2)*(mHp3 -
                   CW2*mHp*MZ2) + MH2*(-3.0*mHp4 + CW2*mHp2*MZ2))*
             TF (MH, mHp, CW*MZ))/(64.0*CW4*(-1.0 + CW2)*
             MZ4*(-mHp2 + 
               CW2*MZ2)*(MH4 + (mHp2 - CW2*MZ2)*(mHp2 - CW2*MZ2) - 
               2.0*MH2*(mHp2 + CW2*MZ2))*pi*pi*pi)  + (aem*mMU*
             mMU*(-3.0*mH2_2*mHp8 + 8.0*CW6*mHp6*MZ4 - 
               2.0*CW4*mHp6*MZ2*(mH2_2 + 4.0*(mHp2 + MZ2)) + 
               CW2*(-mH2_10 + 2.0*mH2_4*mHp4*MZ2 + 8.0*mHp8*MZ2 + 
                  mH2_8*(4.0*mHp2 + MZ2) - 
                  2.0*mH2_6*(mHp4 + 2.0*mHp2*MZ2) + 
                  mH2_2*(2.0*mHp8 + 3.0*mHp6*MZ2)))*
             TF (mH2, mHp, mHp))/(64.0*CW2*(-1.0 + CW2)*mH2_4*
             mHp4*(mH2_2 - 4.0*mHp2)*(mH2_2 - MZ2)*(-mHp2 + 
               CW2*MZ2)*pi*pi*pi) - (aem*mMU*
             mMU*(-mH2_6 + 
               mH2_4*(3.0*mHp2 + CW2*MZ2) + (mHp3 - 
                  CW2*mHp*MZ2)*(mHp3 - CW2*mHp*MZ2) + 
               mH2_2*(-3.0*mHp4 + CW2*mHp2*MZ2))*
             TF (mH2, mHp, CW*MZ))/(64.0*CW4*(-1.0 + CW2)*
             MZ4*(-mHp2 + 
               CW2*MZ2)*(mH2_4 + (mHp2 - CW2*MZ2)*(mHp2 - CW2*MZ2) - 
               2.0*mH2_2*(mHp2 + CW2*MZ2))*pi*pi*
             pi) + (aem*(3.0 - 10.0*CW2 + 8.0*CW4)*(MH2 - mH2_2)*
             mHp2*mMU*mMU*TF (mHp, mHp, MZ))/(64.0*CW2*(-1.0 + CW2)*
             MZ2*(MH2 - MZ2)*(-mH2_2 + MZ2)*(-4.0*mHp2 + MZ2)*pi*pi*
             pi))));
        
 
    aYuk = aYuk2+aYuk3+aYuk4+aYuk5;

    gminus2muNLOB = aEWadd + aNonYuk  + aYuk;
      
    /*std::cout << " aEWadd = " << aEWadd << std::endl;
    std::cout << " aNonYuk = " << aNonYuk << std::endl;
    std::cout << " aYuk = " << aYuk << std::endl;*/

    
     return(gminus2muNLOB.real());
}


double GeneralTHDMMatching::gminus2muNLO() {
    
    updateGTHDMParameters();
     
  //  double gminus2muLOvalue=gminus2muLO();
    double gminus2muNLOFvalue=gminus2muNLOF();
    double gminus2muNLOBvalue = gminus2muNLOB();
    
   /* std::cout << "gminus2muNLOFvalue =" << gminus2muNLOFvalue << std::endl;
    std::cout << "gminus2muNLOBvalue =" << gminus2muNLOBvalue << std::endl;*/
    
   double gminus2muNLO = gminus2muNLOFvalue + gminus2muNLOBvalue;
   
  
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

std::vector<WilsonCoefficient>& GeneralTHDMMatching::CMdbs2() {

    double Mut = myGTHDM.getMut();
    double xt = x_t(Mut);
    double GF=myGTHDM.getGF();
    double MW=myGTHDM.Mw();
    gslpp::complex co = GF / 4. / M_PI * MW * myGTHDM.getCKM().computelamt_s();
    double tanb = myGTHDM.gettanb();
    double mHp2=myGTHDM.getmHp2();
    double xHW=mHp2/(MW*MW);
    double xtH=xt/xHW;
    double SWH=xtH*((2.0*xHW-8.0)*log(xtH)/((1.0-xHW)*(1.0-xtH)*(1.0-xtH))+6.0*xHW*log(xt)/((1.0-xHW)*(1.0-xt)*(1.0-xt))-(8.0-2.0*xt)/((1.0-xt)*(1.0-xtH)))/(tanb*tanb);
    double SHH=xtH*((1.0+xtH)/((1.0-xtH)*(1.0-xtH))+2.0*xtH*log(xtH)/((1.0-xtH)*(1.0-xtH)*(1.0-xtH)))/(tanb*tanb*tanb*tanb);

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
            throw std::runtime_error("THDMMatching::CMdbs2(): order " + out.str() + "not implemented");
    }

    vmcds.push_back(mcdbs2);
    return(vmcds);
}

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

    
         // CPZF. Z-penguins diagrms. Eq. (52)
       
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
    double Mut = myGTHDM.getMut();
    double mHp2 = myGTHDM.getmHp2();
    double MW = myGTHDM.Mw();
    double Mt_muw = myGTHDM.Mrun(Muw, myGTHDM.getQuarks(QCD::TOP).getMass_scale(), 
                        myGTHDM.getQuarks(QCD::TOP).getMass(), FULLNNLO);
    double mt_mt = myGTHDM.Mrun(Mut, myGTHDM.getQuarks(QCD::TOP).getMass_scale(), 
                        myGTHDM.getQuarks(QCD::TOP).getMass(), FULLNNLO);
   double mb=myGTHDM.getQuarks(QCD::BOTTOM).getMass();
   
   double ml=myGTHDM.getLeptons(lepton).getMass();

  //  mlep = SM.getLeptons(lep).getMass();

    
    double xt = (Mt_muw*Mt_muw)/(MW*MW);
    double xHp = (mHp2)/(MW*MW);
    double vev = myGTHDM.v();
    double sW2 = myGTHDM.sW2();
    double mHl=myGTHDM.getMHl();
    double mH1_2=mHl*mHl;
   
    
     //mu contains the missalignemtn dependece. It should be mu -> CR(mu0) - log(mu/mu0). Eq (22)

    double mu = log(MW/myGTHDM.getQ_GTHDM());

  
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
    
    gslpp::complex yl1 = R11  + (R12 + i*R13)*sl;
    gslpp::complex yl2 = R21  + (R22 + i*R23)*sl;
    gslpp::complex yl3 = R31  + (R32 + i*R33)*sl;
    
        
   
    gslpp::complex CSboxU = CSboxBll(xt,  xHp, su, sd, sl);
    gslpp::complex CPboxU = CPboxBll(xt,  xHp, su, sd, sl);
    gslpp::complex CPZU = CPZUBll(xt,  xHp,  sW2, su, sd);
    
    gslpp::complex CSphi1U = yl1.real()*CphiU(xHp,  xt,  vev,  xphi1,  mu,  R11,  R12,  R13,  mH1_2,  lambda3,  Relambda7, Imlambda7, su, sd);
    gslpp::complex CSphi2U = yl2.real()*CphiU(xHp,  xt,  vev,  xphi2,  mu,  R21,  R22,  R23,  mH2_2,  lambda3,  Relambda7, Imlambda7, su, sd);
    gslpp::complex CSphi3U = yl3.real()*CphiU(xHp,  xt,  vev,  xphi3,  mu,  R31,  R32,  R33,  mH3_2,  lambda3,  Relambda7, Imlambda7, su, sd);

    gslpp::complex CPphi1U = i*yl1.imag()*CphiU(xHp,  xt,  vev,  xphi1,  mu,  R11,  R12,  R13,  mH1_2,  lambda3,  Relambda7, Imlambda7, su, sd);
    gslpp::complex CPphi2U = i*yl2.imag()*CphiU(xHp,  xt,  vev,  xphi2,  mu,  R21,  R22,  R23,  mH2_2,  lambda3,  Relambda7, Imlambda7, su, sd);
    gslpp::complex CPphi3U = i*yl3.imag()*CphiU(xHp,  xt,  vev,  xphi3,  mu,  R31,  R32,  R33,  mH3_2,  lambda3,  Relambda7, Imlambda7, su, sd);
    
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


std::vector<WilsonCoefficient>& GeneralTHDMMatching::CMbtaunu() {

    double Muw = myGTHDM.getMuw();
    double GF = myGTHDM.getGF();
    myCKM = myGTHDM.getVCKM();
    double mB = myGTHDM.getMesons(QCD::B_P).getMass();
    double tanb = myGTHDM.gettanb();
    double mHp2=myGTHDM.getmHp2();

    vmcbtaunu = StandardModelMatching::CMbtaunu();
    mcbtaunu.setMu(Muw);
 
    switch (mcbtaunu.getOrder()) {
        case NNLO:
        case NLO:
        case LO:
            mcbtaunu.setCoeff(0, -4.*GF * myCKM(0,2) / negsquareroot(2.) * mB*mB*tanb*tanb/mHp2, LO);
            break;
        default:
            std::stringstream out;
            out << mcbtaunu.getOrder();
            throw std::runtime_error("THDMMatching::CMbtaunu(): order " + out.str() + "not implemented");
    }

    vmcbtaunu.push_back(mcbtaunu);
    return(vmcbtaunu);

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

