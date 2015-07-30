/* 
 * Copyright (C) 2015 SusyFit Collaboration
 * All rights reserved.
 *
 * For the licensing terms see doc/COPYING.
 */

#include "STU.h"
#include "StandardModel.h"

STU::STU(const StandardModel& SM_i, int obsFlag) 
: ThObservable(SM_i), myTHDM(static_cast<const THDM*> (&SM_i))
{
    mycache = new THDMcache();
    if (obsFlag > 0 and obsFlag < 4) obs = obsFlag;
    else throw std::runtime_error("obsFlag in STU() called from ThFactory::ThFactory() can only be 1 (S) or 2 (T) or 3 (U)");
};

double STU::computeThValue()
{
    double mHl=myTHDM->getMHl();
    double mA=myTHDM->getMA();
    double mHh=myTHDM->getMHh();
    double mHp=myTHDM->getMHp();
    double vev=myTHDM->v();
    double sina=myTHDM->computeSina();
    double cosa=myTHDM->computeCosa();
    double tanb=myTHDM->getTanb();
    double sinb=myTHDM->getSinb();
    double cosb=myTHDM->getCosb();
    double sin_ba=myTHDM->getsin_ba();
    double m12_2=myTHDM->getM12_2(); 
    
    double Mz=myTHDM->getMz();


    //////////////////
    
    complex B00prime_Mz_Mz2_mH_mA;
    complex B00prime_Mz_Mz2_mHp_mHp;
    complex B00prime_Mz_Mz2_mh_mA;
    complex B00prime_Mz_Mz2_Mz_mH;
    complex B00prime_Mz_Mz2_Mz_mh;
    
    complex B0prime_Mz_Mz2_Mz_mH;
    complex B0prime_Mz_Mz2_Mz_mh;
    
    double mh = mHl;
    double Mz2 = Mz*Mz;
    double sin2_ba = sin_ba*sin_ba;
    double cos2_ba = 1. - sin2_ba;

    double M_w = myTHDM->Mw();
    double Mw2 = M_w*M_w;
    double s_W2 = myTHDM->sW2(); 
    
    B00prime_Mz_Mz2_mH_mA = - mycache->B00_Mz_Mz2_mH_mA(Mz,mHh,mA) + mycache->B00_Mz_0_mH_mA(Mz,mHh,mA);
    B00prime_Mz_Mz2_mHp_mHp = - mycache->B00_Mz_Mz2_mHp_mHp(Mz,mHp) + mycache->B00_Mz_0_mHp_mHp(Mz,mHp);
    B00prime_Mz_Mz2_mh_mA = - mycache->B00_Mz_Mz2_mh_mA(Mz,mh,mA) + mycache->B00_Mz_0_mh_mA(Mz,mh,mA);
    B00prime_Mz_Mz2_Mz_mH = - mycache->B00_Mz_Mz2_Mz_mH(Mz,mHh) + mycache->B00_Mz_0_Mz_mH(Mz,mHh);
    B00prime_Mz_Mz2_Mz_mh = - mycache->B00_Mz_0_Mz_mh(Mz,mh) + mycache->B00_Mz_0_Mz_mh(Mz,mh);
    B0prime_Mz_Mz2_Mz_mH = mycache->B0_Mz_Mz2_Mz_mH(Mz,mHh) - mycache->B0_Mz_0_Mz_mH(Mz,mHh);
    B0prime_Mz_Mz2_Mz_mh = mycache->B0_Mz_Mz2_Mz_mh(Mz,mh) - mycache->B0_Mz_0_Mz_mh(Mz,mh);
    
    double DeltaS = 1./Mz2/M_PI*(sin2_ba * B00prime_Mz_Mz2_mH_mA.real() - B00prime_Mz_Mz2_mHp_mHp.real()
           + cos2_ba * (B00prime_Mz_Mz2_mh_mA.real() + B00prime_Mz_Mz2_Mz_mH.real()
           - B00prime_Mz_Mz2_Mz_mh.real() - Mz2 * B0prime_Mz_Mz2_Mz_mH.real()
           + Mz2 * B0prime_Mz_Mz2_Mz_mh.real()));

    complex B0_Mz_0_Mz_mH;
    complex B0_Mz_0_Mz_mh;
    complex B0_Mz_0_Mw_mH;
    complex B0_Mz_0_Mw_mh;    
    
    B0_Mz_0_Mw_mH = mycache->B0_Mz_0_Mw_mH(Mz,M_w,mHh);
    B0_Mz_0_Mz_mH = mycache->B0_Mz_0_Mz_mH(Mz,mHh);
    B0_Mz_0_Mz_mh = mycache->B0_Mz_0_Mw_mh(Mz,M_w,mh);
    B0_Mz_0_Mw_mh = mycache->B0_Mz_0_Mw_mh(Mz,M_w,mh); 
    
    double DeltaT = 1. / 16. / M_PI / Mw2 / s_W2 * (F(mHp,mA)
           + sin2_ba * (F(mHp,mHh) - F(mA,mHh)) + cos2_ba * (F(mHp,mh) 
           - F(mA,mh) + F(M_w,mHh) - F(M_w,mh) - F(Mz,mHh) 
           + F(Mz,mh) + 4. * Mz2 * (B0_Mz_0_Mz_mH.real() - B0_Mz_0_Mz_mh.real()) 
           - 4. * Mw2 * (B0_Mz_0_Mw_mH.real() - B0_Mz_0_Mw_mh.real()))); 
     
    complex B00prime_Mz_Mw2_mA_mHp;
    complex B00prime_Mz_Mw2_mHp_mHp;
    complex B00prime_Mz_Mw2_mH_mHp;
    complex B00prime_Mz_Mw2_mh_mHp;
    complex B00prime_Mz_Mw2_Mw_mH;
    complex B00prime_Mz_Mw2_Mw_mh;
    
    complex B0prime_Mz_Mw2_Mw_mH;
    complex B0prime_Mz_Mw2_Mw_mh;
      
    B00prime_Mz_Mw2_mA_mHp = - mycache->B00_Mz_Mw2_mA_mHp(Mz,M_w,mA,mHp) + mycache->B00_Mz_0_mA_mHp(Mz,mA,mHp);
    B00prime_Mz_Mw2_mHp_mHp = - mycache->B00_Mz_Mw2_mHp_mHp(Mz,M_w,mHp) + mycache->B00_Mz_0_mHp_mHp(Mz,mHp);
    B00prime_Mz_Mw2_mh_mHp = - mycache->B00_Mz_Mw2_mh_mHp(Mz,M_w,mh,mHp) + mycache->B00_Mz_0_mh_mHp(Mz,mh,mHp);
    B00prime_Mz_Mw2_Mw_mH = - mycache->B00_Mz_Mw2_Mw_mH(Mz,M_w,mHh) + mycache->B00_Mz_0_Mw_mH(Mz,M_w,mHh);
    B00prime_Mz_Mw2_Mw_mh = - mycache->B00_Mz_Mw2_Mw_mh(Mz,M_w,mh) + mycache->B00_Mz_0_Mw_mh(Mz,M_w,mh);
    B0prime_Mz_Mw2_Mw_mH = mycache->B0_Mz_Mw2_Mw_mH(Mz,M_w,mHh) - mycache->B0_Mz_0_Mw_mH(Mz,M_w,mHh);
    B0prime_Mz_Mw2_Mw_mh = mycache->B0_Mz_Mw2_Mw_mh(Mz,M_w,mh) - mycache->B0_Mz_0_Mw_mh(Mz,M_w,mh);
    B00prime_Mz_Mw2_mH_mHp = - mycache->B00_Mz_Mw2_mH_mHp(Mz,M_w,mHh,mHp) + mycache->B00_Mz_0_mH_mHp(Mz,mHh,mHp);
    
    double DeltaU = - DeltaS + 1. / M_PI / Mz2 * (B00prime_Mz_Mw2_mA_mHp.real()
           - 2. * B00prime_Mz_Mw2_mHp_mHp.real() + sin2_ba * B00prime_Mz_Mw2_mH_mHp.real()
           + cos2_ba * (B00prime_Mz_Mw2_mh_mHp.real() + B00prime_Mz_Mw2_Mw_mH.real()
           - B00prime_Mz_Mw2_Mw_mh.real() - Mw2 * B0prime_Mz_Mw2_Mw_mH.real()
           + Mw2 * B0prime_Mz_Mw2_Mw_mh.real()));    

    ///////////////////
    
    if (obs == 1) return(DeltaS);
    if (obs == 2) return(DeltaT);
    if (obs == 3) return(DeltaU);

    throw std::runtime_error("STU::computeThValue(): Observable type not defined. Can be only any of (1,2,3)");
    return (EXIT_FAILURE);

}


//
//double STU::obliqueS() const {
//  
//    complex B00prime_Mz_Mz2_mH_mA;
//    complex B00prime_Mz_Mz2_mHp_mHp;
//    complex B00prime_Mz_Mz2_mh_mA;
//    complex B00prime_Mz_Mz2_Mz_mH;
//    complex B00prime_Mz_Mz2_Mz_mh;
//    
//    complex B0prime_Mz_Mz2_Mz_mH;
//    complex B0prime_Mz_Mz2_Mz_mh;
//    
//    double mh = mHl;
//    double Mz2 = Mz*Mz;
//    double sin2_ba = sin_ba*sin_ba;
//    double cos2_ba = 1. - sin2_ba;
//    
//    B00prime_Mz_Mz2_mH_mA = - mycache->B00_Mz_Mz2_mH_mA(Mz,mHh,mA) + mycache->B00_Mz_0_mH_mA(Mz,mHh,mA);
//    B00prime_Mz_Mz2_mHp_mHp = - mycache->B00_Mz_Mz2_mHp_mHp(Mz,mHp) + mycache->B00_Mz_0_mHp_mHp(Mz,mHp);
//    B00prime_Mz_Mz2_mh_mA = - mycache->B00_Mz_Mz2_mh_mA(Mz,mh,mA) + mycache->B00_Mz_0_mh_mA(Mz,mh,mA);
//    B00prime_Mz_Mz2_Mz_mH = - mycache->B00_Mz_Mz2_Mz_mH(Mz,mHh) + mycache->B00_Mz_0_Mz_mH(Mz,mHh);
//    B00prime_Mz_Mz2_Mz_mh = - mycache->B00_Mz_0_Mz_mh(Mz,mh) + mycache->B00_Mz_0_Mz_mh(Mz,mh);
//    B0prime_Mz_Mz2_Mz_mH = mycache->B0_Mz_Mz2_Mz_mH(Mz,mHh) - mycache->B0_Mz_0_Mz_mH(Mz,mHh);
//    B0prime_Mz_Mz2_Mz_mh = mycache->B0_Mz_Mz2_Mz_mh(Mz,mh) - mycache->B0_Mz_0_Mz_mh(Mz,mh);
//    
//    double DeltaS = 1./Mz2/M_PI*(sin2_ba * B00prime_Mz_Mz2_mH_mA.real() - B00prime_Mz_Mz2_mHp_mHp.real()
//           + cos2_ba * (B00prime_Mz_Mz2_mh_mA.real() + B00prime_Mz_Mz2_Mz_mH.real()
//           - B00prime_Mz_Mz2_Mz_mh.real() - Mz2 * B0prime_Mz_Mz2_Mz_mH.real()
//           + Mz2 * B0prime_Mz_Mz2_Mz_mh.real()));
//    
//    return DeltaS;
//   
//}
//
//double STU::obliqueT() const {
//    
//    complex B0_Mz_0_Mz_mH;
//    complex B0_Mz_0_Mz_mh;
//    complex B0_Mz_0_Mw_mH;
//    complex B0_Mz_0_Mw_mh;    
//    
//    
//    double M_w = Mw();
//    double mh = mHl;
//    double Mz2 = Mz*Mz;
//    double Mw2 = M_w*M_w;
//    double sin2_ba = sin_ba*sin_ba;
//    double cos2_ba = 1. - sin2_ba;
//    double s_W2 = sW2(); 
//    
//    B0_Mz_0_Mw_mH = mycache->B0_Mz_0_Mw_mH(Mz,M_w,mHh);
//    B0_Mz_0_Mz_mH = mycache->B0_Mz_0_Mz_mH(Mz,mHh);
//    B0_Mz_0_Mz_mh = mycache->B0_Mz_0_Mw_mh(Mz,M_w,mh);
//    B0_Mz_0_Mw_mh = mycache->B0_Mz_0_Mw_mh(Mz,M_w,mh); 
//    
//    double DeltaT = 1. / 16. / M_PI / Mw2 / s_W2 * (F(mHp,mA)
//           + sin2_ba * (F(mHp,mHh) - F(mA,mHh)) + cos2_ba * (F(mHp,mh) 
//           - F(mA,mh) + F(M_w,mHh) - F(M_w,mh) - F(Mz,mHh) 
//           + F(Mz,mh) + 4. * Mz2 * (B0_Mz_0_Mz_mH.real() - B0_Mz_0_Mz_mh.real()) 
//           - 4. * Mw2 * (B0_Mz_0_Mw_mH.real() - B0_Mz_0_Mw_mh.real()))); 
//     
//    return DeltaT;
//}
//
//double STU::obliqueU() const {
//    
//    complex B00prime_Mz_Mw2_mA_mHp;
//    complex B00prime_Mz_Mw2_mHp_mHp;
//    complex B00prime_Mz_Mw2_mH_mHp;
//    complex B00prime_Mz_Mw2_mh_mHp;
//    complex B00prime_Mz_Mw2_Mw_mH;
//    complex B00prime_Mz_Mw2_Mw_mh;
//    
//    complex B0prime_Mz_Mw2_Mw_mH;
//    complex B0prime_Mz_Mw2_Mw_mh;
//    
//    double M_w = Mw();
//    double mh = mHl;
//    double Mz2 = Mz*Mz;
//    double Mw2 = M_w*M_w;//
//    double sin2_ba = sin_ba*sin_ba;
//    double cos2_ba = 1. - sin2_ba;
//      
//    B00prime_Mz_Mw2_mA_mHp = - mycache->B00_Mz_Mw2_mA_mHp(Mz,M_w,mA,mHp) + mycache->B00_Mz_0_mA_mHp(Mz,mA,mHp);
//    B00prime_Mz_Mw2_mHp_mHp = - mycache->B00_Mz_Mw2_mHp_mHp(Mz,M_w,mHp) + mycache->B00_Mz_0_mHp_mHp(Mz,mHp);
//    B00prime_Mz_Mw2_mh_mHp = - mycache->B00_Mz_Mw2_mh_mHp(Mz,M_w,mh,mHp) + mycache->B00_Mz_0_mh_mHp(Mz,mh,mHp);
//    B00prime_Mz_Mw2_Mw_mH = - mycache->B00_Mz_Mw2_Mw_mH(Mz,M_w,mHh) + mycache->B00_Mz_0_Mw_mH(Mz,M_w,mHh);
//    B00prime_Mz_Mw2_Mw_mh = - mycache->B00_Mz_Mw2_Mw_mh(Mz,M_w,mh) + mycache->B00_Mz_0_Mw_mh(Mz,M_w,mh);
//    B0prime_Mz_Mw2_Mw_mH = mycache->B0_Mz_Mw2_Mw_mH(Mz,M_w,mHh) - mycache->B0_Mz_0_Mw_mH(Mz,M_w,mHh);
//    B0prime_Mz_Mw2_Mw_mh = mycache->B0_Mz_Mw2_Mw_mh(Mz,M_w,mh) - mycache->B0_Mz_0_Mw_mh(Mz,M_w,mh);
//    B00prime_Mz_Mw2_mH_mHp = - mycache->B00_Mz_Mw2_mH_mHp(Mz,M_w,mHh,mHp) + mycache->B00_Mz_0_mH_mHp(Mz,mHh,mHp);
//    
//    double DeltaU = - obliqueS() + 1. / M_PI / Mz2 * (B00prime_Mz_Mw2_mA_mHp.real()
//           - 2. * B00prime_Mz_Mw2_mHp_mHp.real() + sin2_ba * B00prime_Mz_Mw2_mH_mHp.real()
//           + cos2_ba * (B00prime_Mz_Mw2_mh_mHp.real() + B00prime_Mz_Mw2_Mw_mH.real()
//           - B00prime_Mz_Mw2_Mw_mh.real() - Mw2 * B0prime_Mz_Mw2_Mw_mH.real()
//           + Mw2 * B0prime_Mz_Mw2_Mw_mh.real()));    
//    
//    return DeltaU;
// 
//}

/////////////////////////////////////////////////////////////////////////////

double STU::F(const double m0, const double m1) const {
    double m12 = m1 * m1;
    double m02 = m0 * m0;
    double F;
    
    if ( m0<=0.0 || m1<0.0 )
        throw std::runtime_error("Invalid argument for THDM::F()\n"); 
    
    if(m0 == 0. && m1 != 0.) {
        F=0.5 * m12;
    } else if(m0 != 0. && m1 == 0.){
        F=0.5 * m02;
    } else if((m0 == 0. && m1 == 0.) || (fabs(m0-m1) < LEPS)){
        F=0.;
    } else if (m0 != 0 && m1 != 0){
        F=0.5 * (m02 + m12) - (m02 * m12) / (m02 - m12) * log(m02 / m12);
    } else
        throw std::runtime_error("Error in THDM::F()");
    return (F);
}
