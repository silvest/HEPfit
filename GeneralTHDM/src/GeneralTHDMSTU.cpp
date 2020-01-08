/*
 * Copyright (C) 2015 HEPfit Collaboration
 *
 *
 * For the licensing terms see doc/COPYING.
 */

#include "GeneralTHDMSTU.h"
#include "GeneralTHDM.h"
#include "GeneralTHDMcache.h"

GeneralTHDMSTU::GeneralTHDMSTU(const StandardModel& SM_i)
: ThObservable(SM_i), myGTHDM(static_cast<const GeneralTHDM*> (&SM_i))

{
    mycache = new GeneralTHDMcache(SM_i);
};

double GeneralTHDMSTU::computeThValue()
{
    return 0.0;
}
/////////////////////////////////////////////////////////////////////////////

double GeneralTHDMSTU::F(const double m02, const double m12) const {
    double F;

    if(m02 == 0. && m12 != 0.) {
        F=0.5 * m12;
    } else if(m02 != 0. && m12 == 0.){
        F=0.5 * m02;
    } else if((m02 == 0. && m12 == 0.) || (fabs(m02-m12) < LEPS)){
        F=0.;
    } else if (m02 != 0 && m12 != 0){
        F=0.5 * (m02 + m12) - (m02 * m12) / (m02 - m12) * log(m02 / m12);
    } else
        throw std::runtime_error("Error in GeneralTHDM::F()");
    return (F);
}

/////////////////////////////////////////////////////////////////////////////
GTHDMDeltaS::GTHDMDeltaS(const StandardModel& SM_i)
: GeneralTHDMSTU(SM_i)
{}

double GTHDMDeltaS::computeThValue()
{
    double mH1_2 = myGTHDM->getMyGTHDMCache()->mH1sq;
    double mH2_2 = myGTHDM->getMyGTHDMCache()->mH2sq;
    double mH3_2 = myGTHDM->getMyGTHDMCache()->mH3sq;
    double mHp2=myGTHDM->getmHp2();
    
    double mHref_2 = myGTHDM->getMyGTHDMCache()->mH1sq;  
    double  R11 = 0.0;
    double  R21 = 0.0;
    double  R31 = myGTHDM->getMyGTHDMCache()->R31_GTHDM;

    
    if(myGTHDM->getSMHiggs()){
           R11 = myGTHDM->getMyGTHDMCache()->R11_GTHDM;
           R21 = myGTHDM->getMyGTHDMCache()->R21_GTHDM;
    }
    else{
           R21 = myGTHDM->getMyGTHDMCache()->R11_GTHDM;
           R11 = myGTHDM->getMyGTHDMCache()->R21_GTHDM;
    }
    
        double  R11_2 = R11*R11;
        double  R21_2 = R21*R21;
        double  R31_2 = R31*R31;



    
    double MZ=myGTHDM->getMz();
    double MZ2 = MZ*MZ;
    
    gslpp::complex B00prime_MZ2_MZ2_mH_2_mH_3;
    gslpp::complex B00prime_MZ2_MZ2_mHp2_mHp2;
    gslpp::complex B00prime_MZ2_MZ2_mH_1_mH_3;
    gslpp::complex B00prime_MZ2_MZ2_mH_1_mH_2;

    
    gslpp::complex B00prime_MZ2_MZ2_MZ2_mH2_2;
    gslpp::complex B00prime_MZ2_MZ2_MZ2_mH1_2;
    gslpp::complex B0prime_MZ2_MZ2_MZ2_mH2_2;
    gslpp::complex B0prime_MZ2_MZ2_MZ2_mH1_2;
    
    gslpp::complex B00prime_MZ2_MZ2_MZ2_mH3_2;
    gslpp::complex B0prime_MZ2_MZ2_MZ2_mH3_2;
    
     gslpp::complex B00prime_MZ2_MZ2_MZ2_mHref_2;
    gslpp::complex B0prime_MZ2_MZ2_MZ2_mHref_2;



    B00prime_MZ2_MZ2_mH_2_mH_3 = mycache->B00_MZ2_MZ2_mHh2_mA2(MZ2,mH2_2,mH3_2) - mycache->B00_MZ2_0_mHh2_mA2(MZ2,mH2_2,mH3_2);
    B00prime_MZ2_MZ2_mHp2_mHp2 = mycache->B00_MZ2_MZ2_mHp2_mHp2(MZ2,mHp2) - mycache->B00_MZ2_0_mHp2_mHp2(MZ2,mHp2);
    B00prime_MZ2_MZ2_mH_1_mH_3 = mycache->B00_MZ2_MZ2_mHl2_mA2(MZ2,mH1_2,mH3_2) - mycache->B00_MZ2_0_mHl2_mA2(MZ2,mH1_2,mH3_2);
    B00prime_MZ2_MZ2_mH_1_mH_2 = mycache->B00_MZ2_MZ2_mHl2_mA2(MZ2,mH1_2,mH2_2) - mycache->B00_MZ2_0_mHl2_mA2(MZ2,mH1_2,mH2_2);

    
    B00prime_MZ2_MZ2_MZ2_mH2_2 = mycache->B00_MZ2_MZ2_MZ2_mHh2(MZ2,mH2_2) - mycache->B00_MZ2_0_MZ2_mHh2(MZ2,mH2_2);
    B00prime_MZ2_MZ2_MZ2_mH1_2 = mycache->B00_MZ2_MZ2_MZ2_mHl2(MZ2,mH1_2) - mycache->B00_MZ2_0_MZ2_mHl2(MZ2,mH1_2);
    B0prime_MZ2_MZ2_MZ2_mH2_2 = mycache->B0_MZ2_MZ2_MZ2_mHh2(MZ2,mH2_2) - mycache->B0_MZ2_0_MZ2_mHh2(MZ2,mH2_2);
    B0prime_MZ2_MZ2_MZ2_mH1_2 = mycache->B0_MZ2_MZ2_MZ2_mHl2(MZ2,mH1_2) - mycache->B0_MZ2_0_MZ2_mHl2(MZ2,mH1_2);
    B00prime_MZ2_MZ2_MZ2_mH3_2 = mycache->B00_MZ2_MZ2_MZ2_mHh2(MZ2,mH3_2) - mycache->B00_MZ2_0_MZ2_mHh2(MZ2,mH3_2);
    B0prime_MZ2_MZ2_MZ2_mH3_2 = mycache->B0_MZ2_MZ2_MZ2_mHh2(MZ2,mH3_2) - mycache->B0_MZ2_0_MZ2_mHh2(MZ2,mH3_2);

    B00prime_MZ2_MZ2_MZ2_mHref_2 = mycache->B00_MZ2_MZ2_MZ2_mHh2(MZ2,mHref_2) - mycache->B00_MZ2_0_MZ2_mHh2(MZ2,mHref_2);
    B0prime_MZ2_MZ2_MZ2_mHref_2 = mycache->B0_MZ2_MZ2_MZ2_mHh2(MZ2,mHref_2) - mycache->B0_MZ2_0_MZ2_mHh2(MZ2,mHref_2);


    
    return 1/MZ2/M_PI*( R11_2*(B00prime_MZ2_MZ2_MZ2_mH1_2.real() - MZ2*B0prime_MZ2_MZ2_MZ2_mH1_2.real())+
            R21_2*(B00prime_MZ2_MZ2_MZ2_mH2_2.real() - MZ2*B0prime_MZ2_MZ2_MZ2_mH2_2.real())+
            R31_2*(B00prime_MZ2_MZ2_MZ2_mH3_2.real() - MZ2*B0prime_MZ2_MZ2_MZ2_mH3_2.real())+
             R11_2*B00prime_MZ2_MZ2_mH_2_mH_3.real() + R21_2*B00prime_MZ2_MZ2_mH_1_mH_3.real()
             +R31_2*B00prime_MZ2_MZ2_mH_1_mH_2.real() - B00prime_MZ2_MZ2_mHp2_mHp2.real()
              - B00prime_MZ2_MZ2_MZ2_mHref_2.real() + MZ2*B0prime_MZ2_MZ2_MZ2_mHref_2.real());
            
                 
}

GTHDMDeltaT::GTHDMDeltaT(const StandardModel& SM_i)
: GeneralTHDMSTU(SM_i)
{}

double GTHDMDeltaT::computeThValue()
{
    
    gslpp::complex I = gslpp::complex::i();
    double mH1_2 = myGTHDM->getMyGTHDMCache()->mH1sq;
    double mH2_2 = myGTHDM->getMyGTHDMCache()->mH2sq;
    double mH3_2 = myGTHDM->getMyGTHDMCache()->mH3sq;
    double mHp2=myGTHDM->getmHp2();
    
    double mHref_2 = myGTHDM->getMyGTHDMCache()->mH1sq;
        
        
        
    double  R11 = 0.0;
    double  R12 = 0.0;
    double  R13 = 0.0;
    double  R21 = 0.0;
    double  R22 = 0.0;
    double  R23 = 0.0;
    double  R31 = myGTHDM->getMyGTHDMCache()->R31_GTHDM;
    double  R32 = myGTHDM->getMyGTHDMCache()->R32_GTHDM;
    double  R33 = myGTHDM->getMyGTHDMCache()->R33_GTHDM;

    
    if(myGTHDM->getSMHiggs()){
           R11 = myGTHDM->getMyGTHDMCache()->R11_GTHDM;
           R12 = myGTHDM->getMyGTHDMCache()->R12_GTHDM;
           R13 = myGTHDM->getMyGTHDMCache()->R13_GTHDM;
           R21 = myGTHDM->getMyGTHDMCache()->R21_GTHDM;
           R22 = myGTHDM->getMyGTHDMCache()->R22_GTHDM;
           R23 = myGTHDM->getMyGTHDMCache()->R23_GTHDM;
    }
    else{
           R21 = myGTHDM->getMyGTHDMCache()->R11_GTHDM;
           R22 = myGTHDM->getMyGTHDMCache()->R12_GTHDM;
           R23 = myGTHDM->getMyGTHDMCache()->R13_GTHDM;
           R11 = myGTHDM->getMyGTHDMCache()->R21_GTHDM;
           R12 = myGTHDM->getMyGTHDMCache()->R22_GTHDM;
           R13 = myGTHDM->getMyGTHDMCache()->R23_GTHDM;
    }
    
        double  R11_2 = R11*R11;
        double  R21_2 = R21*R21;
        double  R31_2 = R31*R31;

    double MZ=myGTHDM->getMz();
    double MZ2 = MZ*MZ;
    
    double MW=mycache->MWGTHDM(myGTHDM->Mw_tree());
    double MW2 = MW*MW;
    double s_W2 = 1.0-mycache->cW2GTHDM(myGTHDM->c02());

    gslpp::complex B0_MZ2_0_MZ2_mH1_2;
    gslpp::complex B0_MZ2_0_MZ2_mH2_2;
    gslpp::complex B0_MZ2_0_MZ2_mH3_2;
    gslpp::complex B0_MZ2_0_MZ2_mHref_2;


    gslpp::complex B0_MZ2_0_MW2_mH1_2;
    gslpp::complex B0_MZ2_0_MW2_mH2_2;
    gslpp::complex B0_MZ2_0_MW2_mH3_2;
    gslpp::complex B0_MZ2_0_MW2_mHref_2;



    B0_MZ2_0_MZ2_mH1_2 = mycache->B0_MZ2_0_MZ2_mHh2(MZ2,mH1_2);
    B0_MZ2_0_MZ2_mH2_2 = mycache->B0_MZ2_0_MZ2_mHl2(MZ2,mH2_2);
    B0_MZ2_0_MZ2_mH3_2 = mycache->B0_MZ2_0_MZ2_mHl2(MZ2,mH3_2);
    B0_MZ2_0_MZ2_mHref_2 = mycache->B0_MZ2_0_MZ2_mHl2(MZ2,mHref_2);


    
    B0_MZ2_0_MW2_mH1_2 = mycache->B0_MZ2_0_MW2_mHh2(MZ2,MW2,mH1_2);
    B0_MZ2_0_MW2_mH2_2 = mycache->B0_MZ2_0_MW2_mHl2(MZ2,MW2,mH2_2);
    B0_MZ2_0_MW2_mH3_2 = mycache->B0_MZ2_0_MW2_mHl2(MZ2,MW2,mH3_2);
    B0_MZ2_0_MW2_mHref_2 = mycache->B0_MZ2_0_MW2_mHl2(MZ2,MW2,mHref_2);

    
    return 1. / 16. / M_PI / MW2 / s_W2*((R12 + I*R13).abs2()*F(mHp2, mH1_2)+(R22 + I*R23).abs2()*F(mHp2, mH2_2)+
             (R32 + I*R33).abs2()*F(mHp2, mH3_2)- R11_2*F(mH2_2, mH3_2) - R21_2*F(mH1_2, mH3_2) - R31_2*F(mH1_2, mH2_2)
            + R11_2*(F(MW2, mH1_2) - F(MZ2, mH1_2) - 4.0*MW2*B0_MZ2_0_MW2_mH1_2.real() + 4.0*MZ2*B0_MZ2_0_MZ2_mH1_2.real())
            + R21_2*(F(MW2, mH2_2) - F(MZ2, mH2_2) - 4.0*MW2*B0_MZ2_0_MW2_mH2_2.real() + 4.0*MZ2*B0_MZ2_0_MZ2_mH2_2.real())
            + R31_2*(F(MW2, mH3_2) - F(MZ2, mH3_2) - 4.0*MW2*B0_MZ2_0_MW2_mH3_2.real() + 4.0*MZ2*B0_MZ2_0_MZ2_mH3_2.real())
            -(F(MW2, mHref_2) - F(MZ2, mHref_2) - 4.0*MW2*B0_MZ2_0_MW2_mHref_2.real() + 4.0*MZ2*B0_MZ2_0_MZ2_mHref_2.real()));

     
}

GTHDMDeltaU::GTHDMDeltaU(const StandardModel& SM_i)
: GeneralTHDMSTU(SM_i)
{
    myDeltaS = new GTHDMDeltaS(SM_i);
}

double GTHDMDeltaU::computeThValue()
{
 
    gslpp::complex I = gslpp::complex::i();
    double mH1_2 = myGTHDM->getMyGTHDMCache()->mH1sq;
    double mH2_2 = myGTHDM->getMyGTHDMCache()->mH2sq;
    double mH3_2 = myGTHDM->getMyGTHDMCache()->mH3sq;
    double mHp2=myGTHDM->getmHp2();
    
    double mHref_2 = myGTHDM->getMyGTHDMCache()->mH1sq;
        
        
    
          
    double  R11 = 0.0;
    double  R12 = 0.0;
    double  R13 = 0.0;
    double  R21 = 0.0;
    double  R22 = 0.0;
    double  R23 = 0.0;
    double  R31 = myGTHDM->getMyGTHDMCache()->R31_GTHDM;
    double  R32 = myGTHDM->getMyGTHDMCache()->R32_GTHDM;
    double  R33 = myGTHDM->getMyGTHDMCache()->R33_GTHDM;

    
    if(myGTHDM->getSMHiggs()){
           R11 = myGTHDM->getMyGTHDMCache()->R11_GTHDM;
           R12 = myGTHDM->getMyGTHDMCache()->R12_GTHDM;
           R13 = myGTHDM->getMyGTHDMCache()->R13_GTHDM;
           R21 = myGTHDM->getMyGTHDMCache()->R21_GTHDM;
           R22 = myGTHDM->getMyGTHDMCache()->R22_GTHDM;
           R23 = myGTHDM->getMyGTHDMCache()->R23_GTHDM;
    }
    else{
           R21 = myGTHDM->getMyGTHDMCache()->R11_GTHDM;
           R22 = myGTHDM->getMyGTHDMCache()->R12_GTHDM;
           R23 = myGTHDM->getMyGTHDMCache()->R13_GTHDM;
           R11 = myGTHDM->getMyGTHDMCache()->R21_GTHDM;
           R12 = myGTHDM->getMyGTHDMCache()->R22_GTHDM;
           R13 = myGTHDM->getMyGTHDMCache()->R23_GTHDM;
    }
    
    double  R11_2 = R11*R11;
    double  R21_2 = R21*R21;
    double  R31_2 = R31*R31;

    double MZ=myGTHDM->getMz();
    double MZ2 = MZ*MZ;
    
    double MW=mycache->MWGTHDM(myGTHDM->Mw_tree());
    double MW2 = MW*MW;


    gslpp::complex B00prime_MZ2_MW2_mH1_2_mHp2;
    gslpp::complex B00prime_MZ2_MW2_mH2_2_mHp2;
    gslpp::complex B00prime_MZ2_MW2_mH3_2_mHp2;
    gslpp::complex B00prime_MZ2_MW2_mHref_2_mHp2;


    
    gslpp::complex B00prime_MZ2_MW2_mHp2_mHp2;
    
    
    gslpp::complex B00prime_MZ2_MW2_MW2_mH1_2;
    gslpp::complex B00prime_MZ2_MW2_MW2_mH2_2;
    gslpp::complex B00prime_MZ2_MW2_MW2_mH3_2;
    gslpp::complex B00prime_MZ2_MW2_MW2_mHref_2;

    
    gslpp::complex B0prime_MZ2_MW2_MW2_mH1_2;
    gslpp::complex B0prime_MZ2_MW2_MW2_mH2_2;
    gslpp::complex B0prime_MZ2_MW2_MW2_mH3_2;
    gslpp::complex B0prime_MZ2_MW2_MW2_mHref_2;



    B00prime_MZ2_MW2_mH1_2_mHp2 = mycache->B00_MZ2_MW2_mHl2_mHp2(MZ2,MW2,mH1_2,mHp2) - mycache->B00_MZ2_0_mHl2_mHp2(MZ2,mH1_2,mHp2);
    B00prime_MZ2_MW2_mH2_2_mHp2 = mycache->B00_MZ2_MW2_mHh2_mHp2(MZ2,MW2,mH2_2,mHp2) - mycache->B00_MZ2_0_mHh2_mHp2(MZ2,mH2_2,mHp2);
    B00prime_MZ2_MW2_mH3_2_mHp2 = mycache->B00_MZ2_MW2_mHh2_mHp2(MZ2,MW2,mH3_2,mHp2) - mycache->B00_MZ2_0_mHh2_mHp2(MZ2,mH3_2,mHp2);
    B00prime_MZ2_MW2_mHref_2_mHp2 = mycache->B00_MZ2_MW2_mHl2_mHp2(MZ2,MW2,mHref_2,mHp2) - mycache->B00_MZ2_0_mHl2_mHp2(MZ2,mHref_2,mHp2);

    B00prime_MZ2_MW2_mHp2_mHp2 = mycache->B00_MZ2_MW2_mHp2_mHp2(MZ2,MW2,mHp2) - mycache->B00_MZ2_0_mHp2_mHp2(MZ2,mHp2);

    B00prime_MZ2_MW2_MW2_mH1_2 = mycache->B00_MZ2_MW2_MW2_mHl2(MZ2,MW2,mH1_2) - mycache->B00_MZ2_0_MW2_mHl2(MZ2,MW2,mH1_2);
    B00prime_MZ2_MW2_MW2_mH2_2 = mycache->B00_MZ2_MW2_MW2_mHh2(MZ2,MW2,mH2_2) - mycache->B00_MZ2_0_MW2_mHh2(MZ2,MW2,mH2_2);
    B00prime_MZ2_MW2_MW2_mH3_2 = mycache->B00_MZ2_MW2_MW2_mHh2(MZ2,MW2,mH3_2) - mycache->B00_MZ2_0_MW2_mHh2(MZ2,MW2,mH3_2);
    B00prime_MZ2_MW2_MW2_mHref_2 = mycache->B00_MZ2_MW2_MW2_mHl2(MZ2,MW2,mHref_2) - mycache->B00_MZ2_0_MW2_mHl2(MZ2,MW2,mHref_2);

    B0prime_MZ2_MW2_MW2_mH1_2 = mycache->B0_MZ2_MW2_MW2_mHl2(MZ2,MW2,mH1_2) - mycache->B0_MZ2_0_MW2_mHl2(MZ2,MW2,mH1_2);
    B0prime_MZ2_MW2_MW2_mH2_2 = mycache->B0_MZ2_MW2_MW2_mHh2(MZ2,MW2,mH2_2) - mycache->B0_MZ2_0_MW2_mHh2(MZ2,MW2,mH2_2);
    B0prime_MZ2_MW2_MW2_mH3_2 = mycache->B0_MZ2_MW2_MW2_mHh2(MZ2,MW2,mH3_2) - mycache->B0_MZ2_0_MW2_mHh2(MZ2,MW2,mH3_2);
    B0prime_MZ2_MW2_MW2_mHref_2 = mycache->B0_MZ2_MW2_MW2_mHl2(MZ2,MW2,mHref_2) - mycache->B0_MZ2_0_MW2_mHl2(MZ2,MW2,mHref_2);

    
    
    return - myDeltaS->computeThValue() + 1. / M_PI / MZ2 * (- R11_2*MW2*B0prime_MZ2_MW2_MW2_mH1_2.real() - R21_2*MW2*B0prime_MZ2_MW2_MW2_mH2_2.real() - R31_2*MW2*B0prime_MZ2_MW2_MW2_mH3_2.real() 
    + MW2*B0prime_MZ2_MW2_MW2_mHref_2.real() - B00prime_MZ2_MW2_MW2_mHref_2.real()
           
   +R11_2* B00prime_MZ2_MW2_MW2_mH1_2.real() +R21_2* B00prime_MZ2_MW2_MW2_mH2_2.real() +R31_2* B00prime_MZ2_MW2_MW2_mH3_2.real()
          
+ (R12+I*R13).abs2()*B00prime_MZ2_MW2_mH1_2_mHp2.real() + (R22+I*R23).abs2()*B00prime_MZ2_MW2_mH2_2_mHp2.real()
 + (R32+I*R33).abs2()*B00prime_MZ2_MW2_mH3_2_mHp2.real() - 2. * B00prime_MZ2_MW2_mHp2_mHp2.real());

   
}