/* 
 * Copyright (C) 2012 SUSYfit Collaboration
 * All rights reserved.
 *
 * For the licensing terms see doc/COPYING.
 */


#include <stdlib.h>
#include <iostream>
#include <stdexcept>
#include <fstream>
#include <EW.h>
#include <sigmamuLEP2.h>
#include <EWSMOneLoopLEP2.h>
#include <PVfunctions.h>
#include <EWSMcache.h>
#include <ZFitter.h>
#include <ZFEWObservables.h>
//#include <TGraph.h>
//#include <TCanvas.h>
//#include <TAxis.h>
//#include <TH1D.h>
//#include <TVirtualFFT.h>
//#include <TF1.h>
//#include <TH2.h>
//#include <TMath.h>
#define EPSILON 0.00001


using namespace std;

/*
 * Simple C++ Test Suite
 */


void setStandardModelparameters(StandardModel& SM_i) {
    map<std::string, double> Parameters;
    
    // 17 parameters defined in StandardModel
    Parameters["GF"] = 1.16637E-5;
    Parameters["mneutrino_1"] = 0.0;
    Parameters["mneutrino_2"] = 0.0;
    Parameters["mneutrino_3"] = 0.0;
    Parameters["melectron"] = 0.54857990943e-3;
    Parameters["mmu"] = 0.1134289256;
    Parameters["mtau"] = 1.77682;
    Parameters["lambda"] = 0.2253;
    Parameters["A"] = 0.808;
    Parameters["rhob"] = 0.132;
    Parameters["etab"] = 0.341;
    Parameters["ale"] = 1.0/137.035999679;
    Parameters["dAle5Mz"] = 0.02758;
    Parameters["mHl"] = 126;
    Parameters["muw"] = 0.0;
    Parameters["mub"] = 0.0;
    Parameters["muc"] = 0.0;
    
    // 26 parameters defined in QCD    
    Parameters["AlsMz"] = 0.1184;
    Parameters["Mz"] = 91.1876;
    Parameters["mup"] = 0.003;
    Parameters["mdown"] = 0.007;
    Parameters["mcharm"] = 1.5;
    Parameters["mstrange"] = 0.1;
    Parameters["mtop"] = 173.0;
    Parameters["mbottom"] = 4.28;
    Parameters["mut"] = 163.4;
    Parameters["mub"] = 4.21;
    Parameters["muc"] = 1.3;
    Parameters["MBd"] = 0.0;
    Parameters["MBs"] = 0.0;
    Parameters["MBp"] = 0.0;
    Parameters["MK0"] = 0.0;
    Parameters["MKp"] = 0.0;
    Parameters["FBs"] = 0.0;
    Parameters["FBsoFBd"] = 0.0;
    Parameters["BBsoBBd"] = 0.0;
    Parameters["BBs1"] = 0.0;
    Parameters["BBs2"] = 0.0;  
    Parameters["BBs3"] = 0.0;
    Parameters["BBs4"] = 0.0;
    Parameters["BBs5"] = 0.0;
    Parameters["BBsscale"] = 0.0;
    Parameters["BBsscheme"] = 0.0;
    Parameters["phiEpsK"] = 0.0;
    Parameters["DeltaMK"] = 0.0;
    Parameters["KbarEpsK"] = 0.0;
    Parameters["Dmk"] = 0.0;
    Parameters["SM_M12D"] = 0.0;
    Parameters["MD"] = 0.0;
    Parameters["FD"] = 0.0;
    Parameters["BD1"] = 0.0;
    Parameters["BD2"] = 0.0;
    Parameters["BD3"] = 0.0;
    Parameters["BD4"] = 0.0;
    Parameters["BD5"] = 0.0;
    Parameters["BDscale"] = 0.0;
    Parameters["BDscheme"] = 0.0;
    Parameters["BK1"] = 0.0;
    Parameters["BK2"] = 0.0;
    Parameters["BK3"] = 0.0;
    Parameters["BK4"] = 0.0;
    Parameters["BK5"] = 0.0;
    Parameters["BKscale"] = 0.0;
    Parameters["BKscheme"] = 0.0;
    Parameters["BK4"] = 0.0;
    Parameters["BK4"] = 0.0;
    Parameters["BK4"] = 0.0;
    Parameters["BK4"] = 0.0;
    
//    //8 parameters defined in THDM
//    Parameters["mHp"] = 580.;
//    Parameters["mA"] = 900.;
//    Parameters["tanb"] = 35.;
//    Parameters["sin_ba"] = 1.;
//    Parameters["m12_2"] = 300.;
//    Parameters["mH"] = 600.;
//    Parameters["lambda6"] = 300.;
//    Parameters["lambda7"] = 300.;
    
   
    
    
    /** Test for alpha_lep **/
    //Parameters["melectron"] = 0.00051099907;
    //Parameters["mmu"] = 0.105658389;
    //Parameters["mtau"] = 1.777;    
    //Parameters["AlsMz"] = 1.0/137.0359895;
    //Parameters["Mz"] = 91.187;    
    
    
    /** To make comparisons with ZFITTER codes **/
    //Parameters["GF"] = 1.16637E-5;  // for GFER=2
    //Parameters["melectron"] = 0.51099907e-3;
    //Parameters["mmu"] = 0.105658389;
    //Parameters["mtau"] = 1.77705;
    //Parameters["mup"] = 0.062;
    //Parameters["mdown"] = 0.083;
    //Parameters["mcharm"] = 1.50;
    //Parameters["mstrange"] = 0.215;
    //Parameters["mbottom"] = 4.70;
    //Parameters["ale"] = 1.0/137.0359895;


    /* TEST for Table 6.1 in hep-ph/0507146*/
    /* flags: AMT4=6, ALEM=2 */
    /* mcMz = 0.55381685, mbMz = 2.8194352 */
    //Parameters["Mz"] = 91.1875;
    //Parameters["mtop"] = 173.0;    
    //Parameters["mHl"] = 500.0;
    //Parameters["AlsMz"] = 0.118;
    //Parameters["dAle5Mz"] = 0.02758;    
    
    
    SM_i.Init(Parameters);
}



int main(int argc, char** argv) {
    
    
    try{
        
        const double sqrt_s[12] = {130.0, 136.0, 161.0, 172.0, 183.0, 189.0, 
                                   192.0, 196.0, 200.0, 202.0, 205.0, 207.0};
    StandardModel* myNewModel;
    myNewModel = new StandardModel();
    
    sigmamuLEP2* mySigmaLEP2[12];
    //mySigmaLEP2 = new sigmamuLEP2();
    
    EW* myEW = new EW(*myNewModel);
    
    EWSMcache* mycache;
    mycache = new EWSMcache(*myNewModel); 
    
    //EWSMcache* cache = new EWSMcache();
    
    EWSMOneLoopLEP2* myLEP2;
    PVfunctions* PV;
    
  //  EWSM* myEWSM;
//    myEWSM = new EWSM(*myNewModel);
    
    myLEP2 = new EWSMOneLoopLEP2(*mycache, *myNewModel);
    
    ZFitter* ZF;
    ZF = new ZFitter(*myNewModel);
    
    
    ZFsigmaMuLEP2*  myZFmu;
    myZFmu = new ZFsigmaMuLEP2(*ZF,sqrt_s[0]);
//    std::ofstream outfile;
//    outfile.open("TestSigmamuLEP2.txt");
    
    setStandardModelparameters(*myNewModel);
    const double s = 130.*130.;
    const double mu = 100.;
    double mf2 = 0.137*0.137;
    const double m1 = 0.137;
    const double m2 = 0.137;
    const double m3 = 90.;
    const double Mw = 80.;
    const double Mz = 91.;
    const double W = 0.;
    const double X = 0.;
    const double Y = 0.;
    double cos_theta = 1.0;
    double t = -s/2.*(1.-cos_theta);
    double u = -s/2.*(1.+cos_theta);
            
    complex i = complex::i();
    complex M = Mz*Mz-i*Mz*myEW->Gamma_Z();
    
    
    
    cout << "myNewModel->getMz()=" << "\t" << myNewModel->getMz() << endl;
    cout << "myNewModel->getAle()=" << "\t" << myNewModel->getAle() << endl;
    
    cout << "myNewModel->Mw()=" << "\t" << myNewModel->Mw() << endl;
    cout << "myNewModel->getMHl()= " << "\t" << myNewModel->getMHl() << endl;
    cout << "Mmu =\t" << myNewModel->getLeptons(myNewModel->MU).getMass()<< endl;
    
    cout << "t = " << t << endl;
    cout << "u = " << u << endl;
    cout << "M1rhok_M1rhopk_l(s,0.5,0.5,0.5,l,cos_theta) = " << myLEP2->M1rhok_M1rhopk_l(s,0.5,0.5,0.5,myNewModel->MU,cos_theta)<< endl;
    cout << "M1rhok_M1rhopk_l(s,-0.5,-0.5,0.5,l,cos_theta) = " << myLEP2->M1rhok_M1rhopk_l(s,-0.5,-0.5,0.5,myNewModel->MU,cos_theta)<< endl;
    cout << "M1rhok_M1rhopk_l(s,-0.5,0.5,0.5,l,cos_theta) = " << myLEP2->M1rhok_M1rhopk_l(s,-0.5,0.5,0.5,myNewModel->MU,cos_theta)<< endl;
    cout << "M1rhok_M1rhopk_l(s,0.5,-0.5,0.5,l,cos_theta) = " << myLEP2->M1rhok_M1rhopk_l(s,0.5,-0.5,0.5,myNewModel->MU,cos_theta)<< endl;
    
//    complex a1 = 1/s*myLEP2->Al(mu,myNewModel->MU,0.5,0.5,s,Mw,W,X,Y)*
//                 myLEP2->ATOTl(mu,myNewModel->MU,0.5,0.5,s,Mw,cos_theta,W,X,Y,myEW->Gamma_Z()).conjugate()*
//                 myLEP2->M1rhok_M1rhopk_l(s,0.5,0.5,0.5,myNewModel->MU,cos_theta);
//    
//    cout << "myLEP2->Al(mu,myNewModel->MU,0.5,0.5,s,Mw,W,X,Y) = " 
//         << myLEP2->Al(mu,myNewModel->MU,0.5,0.5,s,Mw,W,X,Y) << "\n\n\n"<< endl;
//    
//    cout << "myLEP2->ATOTl(mu,myNewModel->MU,0.5,0.5,s,Mw,cos_theta,W,X,Y,myEW->Gamma_Z()).conjugate() = " 
//         << myLEP2->ATOTl(mu,myNewModel->MU,0.5,0.5,s,Mw,cos_theta,W,X,Y,myEW->Gamma_Z()).conjugate() << "\n\n\n"<< endl;
//    
//    cout << "myLEP2->M1rhok_M1rhopk_l(s,0.5,0.5,0.5,myNewModel->MU,cos_theta) = " 
//         << myLEP2->M1rhok_M1rhopk_l(s,0.5,0.5,0.5,myNewModel->MU,cos_theta) << "\n\n\n"<< endl;
//    
//    cout << "a1 = " << a1 << "\n\n\n"<< endl;
    
    double rho = 0.5;
    double k = 0.5;
    double Mw_i = Mw;
    
    
//    complex ATOTl = 1/s*myLEP2->Bl(mu,myNewModel->MU,rho,k,s,Mw_i,W,X,Y)
//            +myLEP2->Cl(mu,myNewModel->MU,rho,k,s,Mw_i,W,X,Y)
//            +myLEP2->Dl_rho(mu,myNewModel->MU,rho,k,s,Mw_i,W,X,Y)
//            +myLEP2->E1(mu,k,s,Mw_i,W,X,Y)*myLEP2->F1_l(mu,rho,s,Mw_i,myNewModel->MU)
//            +myLEP2->E2(mu,k,s,Mw_i,W,X,Y)*myLEP2->F2_l(mu,rho,s,Mw_i,myNewModel->MU)
//            +1/s*myLEP2->E3(s)*myLEP2->Al(mu,myNewModel->MU,rho,k,s,Mw_i,W,X,Y)
//            +myLEP2->E4l(s,myNewModel->MU,mu,Mw_i,W,X,Y)*myLEP2->Lambdal(s,myNewModel->MU,mu)
//            +myLEP2->E5(s,k,mu,Mw_i,W,X,Y)*myLEP2->F5l(s,rho,mu,Mw_i,myNewModel->MU)
//            +myLEP2->E6l(s,myNewModel->MU,mu,Mw_i,W,X,Y)*myLEP2->F6rhol(s,rho,k,myNewModel->MU,cos_theta)
//            +myLEP2->g_rhoe(k,Mw_i)*myLEP2->E5(s,k,mu,Mw_i,W,X,Y)*myLEP2->F7rhol(s,rho,k,Mw_i,myNewModel->MU,cos_theta,myEW->Gamma_Z());
//    
//   
//    cout << "myLEP2->Cl(mu,myNewModel->MU,rho,k,s,Mw_i,W,X,Y) = " 
//            << myLEP2->Cl(mu,myNewModel->MU,rho,k,s,Mw_i,W,X,Y) << "\n\n\n"<< endl;
//    cout << "myLEP2->Dl_rho(mu,myNewModel->MU,rho,k,s,Mw_i,W,X,Y) = " 
//            << myLEP2->Dl_rho(mu,myNewModel->MU,rho,k,s,Mw_i,W,X,Y) << "\n\n\n"<< endl;
//    
//    cout << "myLEP2->F6rhol(s,rho,k,myNewModel->MU,cos_theta) = " 
//            << myLEP2->F6rhol(s,rho,k,myNewModel->MU,cos_theta) << "\n\n\n"<< endl;
//    
//    
//    complex testlog = log(-s+i*EPSILON);
//    complex testlog2 = log(s)+i*M_PI;
//    
//    cout << "\n\n*******************    TEST LOG     ***********************" <<endl;
//    cout <<"testlog =\t"<< testlog <<endl;
//    cout <<"testlog2=\t"<< testlog2 <<endl;    
//    cout << "*********************************************************\n\n"  << endl;
//    
//    
//    complex xxx = myLEP2->g_rhoe(rho,Mw_i)
//                *myLEP2->Lambda2(Mz,s)+2.*myLEP2->FLZl(s, Mw_i)*myLEP2->Chi_Z(mu,s,Mw_i,W,X,Y)/s*myLEP2->g_rhofl(myNewModel->MU,rho,Mw_i);
//    
////    
////    cout << "\n\n*******************    C0     ***********************" <<endl;
////    cout <<"C0  =\t"<< PV->C0(s,0.,Mz,Mz)<<endl;
////    cout <<"***************************************************************"<<endl;
////    
////    
////    
////    
////    //Dl_rho
////    cout << "\n\n*******************    Dl_rho     ***********************" <<endl;
////    cout <<"myLEP2->g_rhoe(rho,Mw_i) =\t"<< myLEP2->g_rhoe(rho,Mw_i)<<endl;
////    cout <<"myLEP2->Lambda2(Mz,s)=\t"<<myLEP2->Lambda2(Mz,s) <<endl;
////    cout <<"myLEP2->FLZl(s, Mw_i)=\t"<<myLEP2->FLZl(s, Mw_i) <<endl;
////    
////    cout <<"\nLambda3(Mw,s)=\t"<<myLEP2->Lambda3(Mw,s)<<endl;
////    
////    cout <<"\nmyLEP2->Chi_Z(mu,s,Mw_i,W,X,Y)=\t"<< myLEP2->Chi_Z(mu,s,Mw_i,W,X,Y)<<endl;
////    cout <<"myLEP2->g_rhofl(l,rho,Mw_i)=\t"<< myLEP2->g_rhofl(myNewModel->MU,rho,Mw_i)<<endl;       
////    cout << "*********************************************************\n\n"  << endl;     
////    
////    //Cl
////    cout << "\n\n*******************       Cl        ***********************" <<endl;
////    cout <<"myLEP2->FLgammal(s, Mw_i) =\t"<< myLEP2->FLgammal(s, Mw_i)<<endl;
////    cout <<"myLEP2->Chi_Z(mu,s,Mw_i,W,X,Y)=\t"<<myLEP2->Chi_Z(mu,s,Mw_i,W,X,Y) <<endl;        
////    cout << "*************************************************************\n\n"  << endl;      
////            
////            
////    cout << "\n\n*******************       F6rhol        ***********************" <<endl;
////    cout <<"myLEP2->V_gammagammal(s,l,cos_theta) =\t"<< myLEP2->V_gammagammal(s,myNewModel->MU,cos_theta)<<endl;
////    cout <<"myLEP2->A_gammagammal(s,l,cos_theta)=\t"<<myLEP2->A_gammagammal(s,myNewModel->MU,cos_theta) <<endl;  
////     cout <<"\nmyLEP2-> Gfunc(s,t)=\t"<<myLEP2->Gfunc(s,t) <<endl;
////      cout <<"myLEP2-> Gfunc(s,u)=\t"<<myLEP2->Gfunc(s,u) <<endl;
////      cout <<"log(-t)=\t"<<log(-t) <<endl;
////       cout <<"log(-u)=\t"<<log(-u) <<endl;
////    cout << "*************************************************************\n\n"<< endl;
////    
////    
////    cout <<"log(t/(s-M*M)=\t"<<log(u/(s-M*M)) <<endl;
            
            
            
            
//    cout << "desigmaSM = " << myNewModel->getEWSM()->dsigmaLEP2_l(myNewModel->MU,s,Mw,0.5,0.,0.,0.,myEW->Gamma_Z()) << "\n\n\n"<< endl;
    
    cout << "\n\ndesigmaSM = " << myNewModel->DsigmaLEP2_l(myNewModel->MU,s,cos_theta,0.,0.,0.,myEW->Gamma_Z())*myZFmu->GeVminus2_to_nb*1000. << "\n\n\n"<< endl;
//    
    cout << "desigmaEW = " << myEW->dsigma_lLEP2(myNewModel->MU,s,0.,0.,0.,cos_theta)*myZFmu->GeVminus2_to_nb*1000. << "\n\n\n"<< endl;
    
    cout << "desigmaZF = " << myZFmu->getThValue() << "\n\n\n"<< endl;
    
//    cout << "PV.C0(130.,0,80,0) = " << PV->C0(130.*130.,0.,80.,0.) << "\n\n\n" << endl;
//    cout << "PV->C0(130.*130.,100.,80,100.) = " << PV->C0(130.*130.,100.,80,100.) << "\n\n\n" << endl;
//    cout << "PV->C0(130.*130.,1.,80,1.) = " << PV->C0(130.*130.,1.,80,1.) << "\n\n\n" << endl;
//    cout << "PV->C0(130.*130.,80.,0.,80.) = " << PV->C0(130.*130.,80.,0.,80.) << "\n\n\n" << endl;
//    cout << "PV->C0(130.*130.,0.,125.,0.) = " << PV->C0(130.*130.,0.,125.,0.) << "\n\n\n" << endl;
//    cout << "PV->C0(130.*130.,90.,0.,125.) = " << PV->C0(130.*130.,90.,0.,125.) << "\n\n\n" << endl;
//    cout << "PV->C0(130.*130.,125.,0.,90.) = " << PV->C0(130.*130.,125.,0.,90.) << "\n\n\n" << endl;
//    cout << "PV->C0(130.*130.,0.,1.,10.) = " << endl;
//    cout << PV->C0(130.*130.,0.,10.,1.) << "\n\n\n" << endl;
//    cout << "PV->C0(130.*130.,1.,2.,10.) = " << endl;
//    cout << PV->C0(130.*130.,1.,2.,10.) << "\n\n\n" << endl;
//    cout << "The last PV->C0(130.*130.,0.137,0.,0.137) = " << PV->C0(130.*130.,0.137,0.,0.137) << "\n\n\n" << endl;
//    cout << "t = " <<myLEP2->t(0.,130.*130.,-1.) << "\n\n\n" << endl;
//   
//    
//    //cout << " C1minus_l " << myLEP2->C1minus_l(10.,myNewModel->MU,0.137,0.137,90.,130.*130.) << "\n\n\n" << endl;
//    
//    
//    
//    const double s = 130.*130.;
//    const double mu = 100.;
//    double mf2 = 0.137*0.137;
//    const double m1 = 0.137;
//    const double m2 = 0.137;
//    const double m3 = 90.;
//    const double Mw = 80.;
//    const double W = 0.;
//    const double X = 0.;
//    const double Y = 0.;
//    
//    complex C0 = -PV->C0(s,m1,m3,m2);
//    
//    cout << "C0 = " << C0 << "\n\n\n" << endl;
//    cout << "PV->B0(mu,mf2,m3,m1) = " << PV->B0(mu,mf2,m3,m1) << "\n\n\n" << endl;
//    cout << "PV->B0(mu,mf2,m3,m2) = " << PV->B0(mu,mf2,m3,m2) << "\n\n\n" << endl;
//    cout << "(m2*m2-m1*m1) = " << (m2*m2-m1*m1) << "\n\n\n" << endl;
//    //cout << "C0 = " << C0 << "\n\n\n" << endl;
//    
//    complex x = 1./2./s*(PV->B0(mu,mf2,m3,m1)-PV->B0(mu,mf2,m3,m2) + (m2*m2-m1*m1)*C0);
//    complex y = myLEP2->C1minus_l(mu, myNewModel->MU,m1,m2,m3,s);
//    complex z = myLEP2->C1plus_l(mu,myNewModel->MU,m1,m2,m3,s);
//    
//    cout << "Chi_gamma = " << myLEP2->Chi_gamma(mu,s,Mw,W,X,Y) << "\n\n\n" << endl;
//    
//    
//
    
//    cout << "PV->D0(0.,0.,0,0,130.*130.,-t,Mz,0.,Mz,0.137) = " << PV->D0(0.,0.,0,0,130.*130.,-t,Mz,0.,Mz,0.137) << "\n\n\n" << endl;
//    cout << "PV->D0(0.,0.,0,0,130.*130.,-t,Mz,0.,Mz,0.) = " << PV->D0(0.,0.,0,0,130.*130.,-t,Mz,0.,Mz,0.) << "\n\n\n" << endl;
//    
//    
//    
//    cout << "t = " << myLEP2->t(0.137,130.*130.,1.) 
//           //<< "\t" << "x = " << x << "\n\n\n" 
//             << "\t" << "t1 = " << t << "\n\n\n"
//             //<< "\t" << "z = " << z << "\n\n\n"
//            << endl;
//    
    
    
//    cout << "sqrt{s}           sigma(mu)"
//             << endl;
//    for (int i=0; i<1; i++) {
//            mySigmaLEP2[i] = new sigmamuLEP2(*myEW, sqrt_s[i]);
//            
//            cout << "\t" << sqrt_s[i]
//                 << "\t sigmamu130 =  " << mySigmaLEP2[i]->getThValue() 
//                 << endl;
//            
//        }
    
  

    return (EXIT_SUCCESS);
    
    } catch (const runtime_error& e) {
        cerr << e.what() << endl;
        return EXIT_FAILURE;
    }
}

