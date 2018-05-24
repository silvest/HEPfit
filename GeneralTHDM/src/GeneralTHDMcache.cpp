/* 
 * Copyright (C) 2016 HEPfit Collaboration
 * All rights reserved.
 *
 * For the licensing terms see doc/COPYING.
 */

#include "GeneralTHDMcache.h"

GeneralTHDMcache::GeneralTHDMcache(const StandardModel& SM_i)
:   Mu_GTHDM(3,3,0.), Md_GTHDM(3,3,0.), Ml_GTHDM(3,3,0.),
    Nu_GTHDM(3,3,0.), Nd_GTHDM(3,3,0.), Nl_GTHDM(3,3,0.),
    Yu1_GTHDM(3,3,0.), Yu2_GTHDM(3,3,0.), Yd1_GTHDM(3,3,0.), Yd2_GTHDM(3,3,0.),
    Yl1_GTHDM(3,3,0.), Yl2_GTHDM(3,3,0.),
    myGTHDM(static_cast<const GeneralTHDM*> (&SM_i))
{}

GeneralTHDMcache::~GeneralTHDMcache()
{}

void GeneralTHDMcache::updateCache()
{
    mHl=myGTHDM->getMHl();
    mH1sq=mHl*mHl;
    mH2sq=myGTHDM->getmH2sq();
    mH3sq=myGTHDM->getmH3sq();
    vev=myGTHDM->v();
    tanb=myGTHDM->gettanb();
    cosb=myGTHDM->getcosb();
    sinb=myGTHDM->getsinb();
    cosa1=myGTHDM->getcosalpha1();
    sina1=myGTHDM->getsinalpha1();
    cosa2=myGTHDM->getcosalpha2();
    sina2=myGTHDM->getsinalpha2();
    cosa3=myGTHDM->getcosalpha3();
    sina3=myGTHDM->getsinalpha3();
    Relambda5=myGTHDM->getRelambda5();
    Imlambda5=myGTHDM->getImlambda5();
    mHpsq=myGTHDM->getmHp2();
    Relambda6=myGTHDM->getRelambda6();
    Relambda7=myGTHDM->getRelambda7();

    /*The Mij_2 are defined such that Msqdiag = -2*RT*M_2*R with the rotation Matrix R
     * and Msqdiag containing the physical mass squares on the diagonal. */

    M11_2 = -0.5*(mH1sq*cosa1*cosa1*cosa2*cosa2
                  +mH2sq*sina1*sina1*cosa2*cosa2 + mH3sq*sina2*sina2);
    M12_2 = 0.5*cosa2*((mH1sq-mH2sq)*cosa1*cosa3*sina1
                       +(-mH3sq+mH1sq*cosa1*cosa1+mH2sq*sina1*sina1)*sina2*sina3);
    M13_2 = 0.5*cosa2*(cosa3*(-mH3sq+mH1sq*cosa1*cosa1+mH2sq*sina1*sina1)*sina2
                       +(mH2sq-mH1sq)*cosa1*sina1*sina3);
    M22_2 = -0.5*(mH3sq*cosa2*cosa2*sina3*sina3
                  +mH1sq*(cosa3*sina1+cosa1*sina2*sina3)*(cosa3*sina1+cosa1*sina2*sina3)
                  +mH2sq*(cosa1*cosa3-sina1*sina2*sina3)*(cosa1*cosa3-sina1*sina2*sina3));
    M23_2 = 0.5*((mH2sq-mH1sq)*cosa1*(cosa3*cosa3-sina3*sina3)*sina1*sina2
                 +cosa1*cosa1*cosa3*(mH2sq-mH1sq*sina2*sina2)*sina3
                 -cosa3*sina3*(mH3sq*cosa2*cosa2+sina1*sina1*(-mH1sq+mH2sq*sina2*sina2)));
    M33_2 = -0.5*(mH3sq*cosa2*cosa2*cosa3*cosa3
                  +mH2sq*(cosa3*sina1*sina2+cosa1*sina3)*(cosa3*sina1*sina2+cosa1*sina3)
                  +mH1sq*(cosa1*cosa3*sina2-sina1*sina3)*(cosa1*cosa3*sina2-sina1*sina3));

    //Remaining general potential parameters
    m11sq     = M11_2 - M33_2 - M12_2*tanb + 0.5*Relambda5*vev*vev
                + (M33_2-0.5*Relambda5*vev*vev)*(cosb*cosb-sinb*sinb)
                + 0.5*vev*vev*((Relambda6-Relambda7)*sinb*cosb+Relambda7*tanb);

    m22sq     = M11_2 - M33_2 + M12_2/tanb + 0.5*Relambda5*vev*vev
                - (M33_2-0.5*Relambda5*vev*vev)*(cosb*cosb-sinb*sinb)
                + 0.25*vev*vev*(Relambda6+Relambda7+(Relambda6-Relambda7)*(cosb*cosb-sinb*sinb))/tanb;

    Rem12sq   = 0.25*vev*vev*(Relambda6+Relambda7+(Relambda6-Relambda7)*(cosb*cosb-sinb*sinb))
                - (2.0*M33_2-Relambda5*vev*vev)*sinb*cosb;

    Imm12sq   = M13_2;

    lambda1   = (-2.0*(M11_2-M22_2+M33_2) + Relambda5*vev*vev
                 - (2.0*M22_2-2.0*M33_2+Relambda5*vev*vev)/(cosb*cosb)
                 + (4.0*M12_2-2.0*Relambda6*vev*vev)*tanb)/(vev*vev);

    lambda2   = (-2.0*(M11_2-M22_2+M33_2) + Relambda5*vev*vev
                 - (2.0*M22_2-2.0*M33_2+Relambda5*vev*vev)/(sinb*sinb)
                 - (4.0*M12_2+2.0*Relambda7*vev*vev)/tanb)/(vev*vev);

    lambda3   = -(2.0*(M11_2-M22_2-M33_2-mHpsq) + Relambda5*vev*vev
                  + (2.0*M12_2+Relambda6*vev*vev)/tanb
                  - (2.0*M12_2-Relambda7*vev*vev)*tanb)/(vev*vev);

    lambda4   = Relambda5 - (2.0*mHpsq+4.0*M33_2)/(vev*vev);

    Imlambda6 = (2.0*M13_2-(2.0*M23_2+0.5*Imlambda5*vev*vev)*tanb)/(vev*vev);

    Imlambda7 = 2.0*M13_2/(vev*vev) + (-0.5*Imlambda5+(2.0*M23_2)/(vev*vev))/tanb;

    //Higgs potential parameters

    m11sqH     = M11_2;
    m22sqH     = M11_2-2.0*M33_2+Relambda5*vev*vev
                 +(M12_2+0.5*(Relambda6*vev*vev))/tanb
                 -(M12_2-0.5*(Relambda7*vev*vev))*tanb;
    Rem12sqH   = -M12_2;
    Imm12sqH   = M13_2;
    lambda1H   = -2.0*M11_2/(vev*vev);
    lambda2H   = -2.0*((2.0*M12_2+Relambda7*vev*vev)/tanb
                       +M11_2-4.0*M22_2+4.0*M33_2-2.0*Relambda5*vev*vev
                       +(M22_2-M33_2+0.5*Relambda5*vev*vev)/(sinb*sinb*cosb*cosb)
                       -2.0*M12_2*tanb+Relambda6*vev*vev*tanb)/(vev*vev);
    lambda3H   = -((2.0*(M11_2-2.0*M33_2-mHpsq+Relambda5*vev*vev)
                    +(2.0*M12_2+Relambda6*vev*vev)/tanb
                    -(2.0*M12_2-Relambda7*vev*vev)*tanb)/(vev*vev));
    lambda4H   = -2.0*(M22_2+M33_2+mHpsq)/(vev*vev);
    Relambda5H = -2.0*(M22_2-M33_2)/(vev*vev);
    Imlambda5H = 4.0*M23_2/(vev*vev);
    Relambda6H = -2.0*M12_2/(vev*vev);
    Imlambda6H = 2.0*M13_2/(vev*vev);
    Relambda7H = (-2.0*M12_2+(Relambda6-Relambda7)*vev*vev
                  +(2.0*M22_2-2.0*M33_2+Relambda5*vev*vev)*(tanb-1.0/tanb))/(vev*vev);
    Imlambda7H = 2.0*(M13_2-M23_2*(tanb-1.0/tanb))/(vev*vev)-0.5*Imlambda5/(sinb*cosb);







//    R11_GTHDM = cosalpha1*cosalpha2;
//    R12_GTHDM = sinalpha1*cosalpha2;
//    R13_GTHDM = -sinalpha2;
//    R21_GTHDM = cosalpha1*sinalpha2*sinalpha3 - sinalpha1*cosalpha3;
//    R22_GTHDM = sinalpha1*sinalpha2*sinalpha3 + cosalpha1*cosalpha3;
//    R23_GTHDM = cosalpha2*sinalpha3;
//    R31_GTHDM = cosalpha1*sinalpha2*cosalpha3 + sinalpha1*sinalpha3;
//    R32_GTHDM = sinalpha1*sinalpha2*cosalpha3 - cosalpha1*sinalpha3;
//    R33_GTHDM = cosalpha2*cosalpha3;

//    M13_2 = -vev*vev*(sinb*cosb*Imlambda5 + cosb*cosb*Imlambda6 + sinb*sinb*Imlambda7);
//    M23_2 = -vev*vev*((cosb*cosb - sinb*sinb)*Imlambda5 + 2.*sinb*cosb*(Imlambda7 - Imlambda6))/2.;

//        std::cout<<"mH1_2 before ordering = "<<mH1_2<<std::endl;
//        std::cout<<"mH2_2 before ordering = "<<mH2_2<<std::endl;
//        std::cout<<"mH3_2 before ordering = "<<mH3_2<<std::endl;

    if(mH1sq<mH3sq && mH3sq<mH2sq)
    {
        //1<3<2 swap 2 and 3
        mHlight_2  = mH1sq;
        mHmedium_2 = mH3sq;
        mHheavy_2  = mH2sq;
    }
    else if(mH3sq<mH2sq && mH2sq<mH1sq)
    {
        //3<2<1 swap 1 and 3
        mHlight_2  = mH3sq;
        mHmedium_2 = mH2sq;
        mHheavy_2  = mH1sq;
    }
    else if(mH2sq<mH1sq && mH1sq<mH3sq)
    {
        //2<1<3 swap 1 and 2
        mHlight_2  = mH2sq;
        mHmedium_2 = mH1sq;
        mHheavy_2  = mH3sq;
    }
    else if(mH2sq<mH3sq && mH3sq<mH1sq)
    {
        //2<3<1: 3->2, 1->3, 2->1
        mHlight_2  = mH2sq;
        mHmedium_2 = mH3sq;
        mHheavy_2  = mH1sq;
    }
    else if(mH3sq<mH1sq && mH1sq<mH2sq)
    {
        //3<1<2 3->1, 1->2, 2->3
        mHlight_2  = mH3sq;
        mHmedium_2 = mH1sq;
        mHheavy_2  = mH2sq;
    }
    
    else if(mH1sq<mH2sq && mH2sq<mH3sq)
    {
        //1<2<3 ok
        mHlight_2  = mH1sq;
        mHmedium_2 = mH2sq;
        mHheavy_2  = mH3sq;
    }

//    M11_2 = (mH1_2*cosalpha1*cosalpha1*cosalpha2*cosalpha2 + mH2_2*sinalpha1*sinalpha1*cosalpha2*cosalpha2 + mH3_2*sinalpha2*sinalpha2);
//
//    M12_2 = (mH1_2*cosalpha1*cosalpha2*(cosalpha1*sinalpha2*sinalpha3 - cosalpha3*sinalpha1)
//             + mH2_2*cosalpha2*sinalpha1*(cosalpha1*cosalpha3 + sinalpha1*sinalpha2*sinalpha3)
//             - mH3_2*cosalpha2*sinalpha2*sinalpha3);
//
////    M13_2 = ?
//
//    M22_2 = (mH1_2*(cosalpha1*sinalpha2*sinalpha3 - cosalpha3*sinalpha1)*(cosalpha1*sinalpha2*sinalpha3 - cosalpha3*sinalpha1)
//             + mH2_2*(cosalpha1*cosalpha3 + sinalpha1*sinalpha2*sinalpha3)*(cosalpha1*cosalpha3 + sinalpha1*sinalpha2*sinalpha3)
//             + mH3_2*cosalpha2*cosalpha2*sinalpha3*sinalpha3);
//
////    M23_2 = ?
//
//    M33_2 = (mH1_2*(cosalpha1*cosalpha3*sinalpha2 + sinalpha1*sinalpha3)*(cosalpha1*cosalpha3*sinalpha2 + sinalpha1*sinalpha3)
//             + mH2_2*(cosalpha3*sinalpha1*sinalpha2 - cosalpha1*sinalpha3)*(cosalpha3*sinalpha1*sinalpha2 - cosalpha1*sinalpha3)
//             + mH3_2*cosalpha2*cosalpha2*cosalpha3*cosalpha3);
//
//    m11_2_GTHDM = M2_GTHDM*(1. - cosb*cosb + 3.*sinb*sinb)/4. + (M12_2*tanb - M11_2)/2.;
//    m22_2_GTHDM = M2_GTHDM*(1. + 3.*cosb*cosb - sinb*sinb)/4. - (M12_2/tanb + M11_2)/2.;
//    Imm12_2_GTHDM = 0.5*(cosb*sinb*Imlambda5 + cosb*cosb*Imlambda6 + sinb*sinb*Imlambda7)*vev*vev;
//    lambda1_GTHDM = (M11_2 + tanb*tanb*(M22_2-M2_GTHDM) - 2.0*tanb*M12_2)/(vev*vev) + tanb*(tanb*tanb*Relambda7 - 3.0*Relambda6)/2.0;
//    lambda2_GTHDM = (M11_2 + (M22_2-M2_GTHDM)/(tanb*tanb) + 2.0*M12_2/tanb)/(vev*vev) + (0.5*Relambda6/(tanb*tanb) - 1.5*Relambda7)/tanb;
//    lambda3_GTHDM = (M11_2 - M22_2 - M2_GTHDM + (1.0/tanb - tanb)*M12_2 + 2.0*mHp2)/(vev*vev) - (Relambda6/tanb + tanb*Relambda7)/2.0;
//    lambda4_GTHDM = (M2_GTHDM + M33_2 - 2.0*mHp2)/(vev*vev) - 0.5*(Relambda6/tanb + tanb*Relambda7);
//    Relambda5_GTHDM = (M2_GTHDM - M33_2)/(vev*vev) - 0.5*(Relambda6/tanb + tanb*Relambda7);
//    
    Mu_GTHDM.assign(0,0, myGTHDM->getQuarks(QCD::UP).getMass());
    Mu_GTHDM.assign(1,1, myGTHDM->getQuarks(QCD::CHARM).getMass());
    Mu_GTHDM.assign(2,2, myGTHDM->getQuarks(QCD::TOP).getMass());
    
    Md_GTHDM.assign(0,0, myGTHDM->getQuarks(QCD::DOWN).getMass());
    Md_GTHDM.assign(1,1, myGTHDM->getQuarks(QCD::STRANGE).getMass());
    Md_GTHDM.assign(2,2, myGTHDM->getQuarks(QCD::BOTTOM).getMass());
    
    Ml_GTHDM.assign(0,0, myGTHDM->getLeptons(StandardModel::ELECTRON).getMass());
    Ml_GTHDM.assign(1,1, myGTHDM->getLeptons(StandardModel::MU).getMass());
    Ml_GTHDM.assign(2,2, myGTHDM->getLeptons(StandardModel::TAU).getMass());
    
    if(myGTHDM->getATHDMflag() == true)
    {
        sigmau_ATHDM = myGTHDM->getNu_33()/myGTHDM->getQuarks(StandardModel::TOP).getMass();
        sigmad_ATHDM = myGTHDM->getNd_33()/myGTHDM->getQuarks(StandardModel::DOWN).getMass();
        sigmal_ATHDM = myGTHDM->getNl_33()/myGTHDM->getLeptons(StandardModel::TAU).getMass();
        
        Nu_GTHDM.assign(0,0, sigmau_ATHDM*Mu_GTHDM(0,0));
        Nu_GTHDM.assign(1,1, sigmau_ATHDM*Mu_GTHDM(1,1));
        Nu_GTHDM.assign(2,2, sigmau_ATHDM*Mu_GTHDM(2,2));
        
        Nd_GTHDM.assign(0,0, sigmad_ATHDM*Mu_GTHDM(0,0));
        Nd_GTHDM.assign(1,1, sigmad_ATHDM*Md_GTHDM(1,1));
        Nd_GTHDM.assign(2,2, sigmad_ATHDM*Mu_GTHDM(2,2));
        
        Nl_GTHDM.assign(0,0, sigmal_ATHDM*Ml_GTHDM(0,0));
        Nl_GTHDM.assign(1,1, sigmal_ATHDM*Ml_GTHDM(1,1));
        Nl_GTHDM.assign(2,2, sigmal_ATHDM*Ml_GTHDM(2,2));
    }
    else
    {
        Nu_GTHDM.assign(0,0, myGTHDM->getNu_11());
        Nu_GTHDM.assign(0,1, myGTHDM->getNu_12());
        Nu_GTHDM.assign(0,2, myGTHDM->getNu_13());
        Nu_GTHDM.assign(1,0, myGTHDM->getNu_21());
        Nu_GTHDM.assign(1,1, myGTHDM->getNu_22());
        Nu_GTHDM.assign(1,2, myGTHDM->getNu_23());
        Nu_GTHDM.assign(2,0, myGTHDM->getNu_31());
        Nu_GTHDM.assign(2,1, myGTHDM->getNu_32());
        Nu_GTHDM.assign(2,2, myGTHDM->getNu_33());
    
        Nd_GTHDM.assign(0,0, myGTHDM->getNd_11());
        Nd_GTHDM.assign(0,1, myGTHDM->getNd_12());
        Nd_GTHDM.assign(0,2, myGTHDM->getNd_13());
        Nd_GTHDM.assign(1,0, myGTHDM->getNd_21());
        Nd_GTHDM.assign(1,1, myGTHDM->getNd_22());
        Nd_GTHDM.assign(1,2, myGTHDM->getNd_23());
        Nd_GTHDM.assign(2,0, myGTHDM->getNd_31());
        Nd_GTHDM.assign(2,1, myGTHDM->getNd_32());
        Nd_GTHDM.assign(2,2, myGTHDM->getNd_33());
    
        Nl_GTHDM.assign(0,0, myGTHDM->getNl_11());
        Nl_GTHDM.assign(0,1, myGTHDM->getNl_12());
        Nl_GTHDM.assign(0,2, myGTHDM->getNl_13());
        Nl_GTHDM.assign(1,0, myGTHDM->getNl_21());
        Nl_GTHDM.assign(1,1, myGTHDM->getNl_22());
        Nl_GTHDM.assign(1,2, myGTHDM->getNl_23());
        Nl_GTHDM.assign(2,0, myGTHDM->getNl_31());
        Nl_GTHDM.assign(2,1, myGTHDM->getNl_32());
        Nl_GTHDM.assign(2,2, myGTHDM->getNl_33());
    }
       
    //Definition of Yukawa matrices
    Yu1_GTHDM = (cosb*Mu_GTHDM - sinb*Nu_GTHDM)*sqrt(2.)/vev;
    Yu2_GTHDM = (cosb*Nu_GTHDM + sinb*Mu_GTHDM)*sqrt(2.)/vev;
    Yd1_GTHDM = (cosb*Md_GTHDM - sinb*Nd_GTHDM)*sqrt(2.)/vev;
    Yd2_GTHDM = (cosb*Nd_GTHDM + sinb*Md_GTHDM)*sqrt(2.)/vev;
    Yl1_GTHDM = (cosb*Ml_GTHDM - sinb*Nl_GTHDM)*sqrt(2.)/vev;
    Yl2_GTHDM = (cosb*Nl_GTHDM + sinb*Ml_GTHDM)*sqrt(2.)/vev;

//    std::cout<<"mH1_2 = "<<mH1_2<<std::endl;
//    std::cout<<"mH2_2 = "<<mH2_2<<std::endl;
//    std::cout<<"mH3_2 = "<<mH3_2<<std::endl;
//    std::cout<<"M11_2 = "<<M11_2<<std::endl;
//    std::cout<<"M12_2 = "<<M12_2<<std::endl;
//    std::cout<<"M22_2 = "<<M22_2<<std::endl;
//    std::cout<<"M33_2 = "<<M33_2<<std::endl;
//    std::cout<<"m11_2_GTHDM = "<<m11_2_GTHDM<<std::endl;
//    std::cout<<"m22_2_GTHDM = "<<m22_2_GTHDM<<std::endl;
//    std::cout<<"Imm12_2_GTHDM = "<<Imm12_2_GTHDM<<std::endl;
//    std::cout<<"lambda1_GTHDM = "<<lambda1_GTHDM<<std::endl;
//    std::cout<<"lambda2_GTHDM = "<<lambda2_GTHDM<<std::endl;
//    std::cout<<"lambda3_GTHDM = "<<lambda3_GTHDM<<std::endl;
//    std::cout<<"lambda4_GTHDM = "<<lambda4_GTHDM<<std::endl;
//    std::cout<<"Relambda5_GTHDM = "<<Relambda5_GTHDM<<std::endl;

//    Q_THDM=myTHDM->getQ_THDM();
////    m12_2=myTHDM->getm12_2();
//    MW=MWTHDM(myTHDM->Mw_tree());
//    cW2=cW2THDM(myTHDM->c02());
//    Ale=myTHDM->getAle();
//    Als=myTHDM->getAlsMz();
//    Mt=myTHDM->getQuarks(QCD::TOP).getMass();
//    Mb=myTHDM->getQuarks(QCD::BOTTOM).getMass();   
//    Mtau=myTHDM->getLeptons(StandardModel::TAU).getMass();
//    Mc=myTHDM->getQuarks(QCD::CHARM).getMass();
//    Ms=myTHDM->getQuarks(QCD::STRANGE).getMass();
//    Mmu=myTHDM->getLeptons(StandardModel::MU).getMass();
//    Mu=myTHDM->getQuarks(QCD::UP).getMass();
//    Md=myTHDM->getQuarks(QCD::DOWN).getMass();
//    Me=myTHDM->getLeptons(StandardModel::ELECTRON).getMass();
//    MZ=myTHDM->getMz();
//    modelflag=myTHDM->getModelTypeflag();
//    WFRflag=myTHDM->getWFRflag();

}
