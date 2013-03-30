/* 
 * Copyright (C) 2012 SusyFit Collaboration
 * All rights reserved.
 *
 * For the licensing terms see doc/COPYING.
 */

#include "Spectrum.h"
#include "SUSY.h"
#include <gslpp.h>
#include <iostream>
#include <fstream>

using namespace gslpp;

Spectrum::Spectrum() {
}

void Spectrum::CalcSpectrum(SUSY & SUSY){
    matrix<complex> CKM(SUSY.getVCKM());
    matrix<double> Mu(3,3,0.);
    matrix<double> Md(3,3,0.);
    
//    Mu(0,0) = SUSY.Mrun(SUSY.GetQ(),SUSY.getQuarks(QCD::UP).getMass_scale(),SUSY.getQuarks(QCD::UP).getMass());
//    Mu(1,1) = SUSY.Mrun(SUSY.GetQ(),SUSY.getQuarks(QCD::CHARM).getMass());
//    Mu(2,2) = SUSY.Mrun(SUSY.GetQ(),SUSY.getQuarks(QCD::TOP).getMass());
//       
//    Md(0, 0) = SUSY.Mrun(SUSY.GetQ(), SUSY.getQuarks(QCD::DOWN).getMass_scale(), SUSY.getQuarks(QCD::DOWN).getMass());
//    Md(1, 1) = SUSY.Mrun(SUSY.GetQ(), SUSY.getQuarks(QCD::STRANGE).getMass_scale(), SUSY.getQuarks(QCD::STRANGE).getMass());
//    Md(2, 2) = SUSY.Mrun(SUSY.GetQ(), SUSY.getQuarks(QCD::BOTTOM).getMass());

      Mu(0,0) = 0;//SUSY.getQuarks(QCD::UP).getMass();
      Mu(1,1) = 0;//SUSY.getQuarks(QCD::CHARM).getMass();
      Mu(2,2) = SUSY.getQuarks(QCD::TOP).getMass();
      
      Md(0,0) = 0;//SUSY.getQuarks(QCD::DOWN).getMass();
      Md(1,1) = 0;//SUSY.getQuarks(QCD::STRANGE).getMass();
      Md(2,2) = SUSY.getQuarks(QCD::BOTTOM).getMass();
    
    //Up-type squark masses
    matrix<complex> uLL(CKM * SUSY.MsQ2 * CKM.hconjugate() + Mu * Mu + matrix<complex>::Id(3)*(.5 - 2. / 3. * SUSY.s02()) * SUSY.getMz() * SUSY.getMz() * cos(2. * atan(SUSY.getTanb())));
    matrix<complex> uRR(SUSY.MsU2 + Mu * Mu + matrix<complex>::Id(3)*(2. / 3. * SUSY.s02()) * SUSY.getMz() * SUSY.getMz() * cos(2. * atan(SUSY.getTanb())));
    matrix<complex> uLR(SUSY.v2() / sqrt(2.) * SUSY.GetTU().hconjugate() - SUSY.getMuH() * Mu / SUSY.getTanb());

    // Down-type squark masses

    matrix<complex> dLL(SUSY.MsQ2 + Md * Md + matrix<complex>::Id(3)*(-.5 + 1. / 3. * SUSY.s02()) * SUSY.getMz() * SUSY.getMz() * cos(2. * atan(SUSY.getTanb())));
    matrix<complex> dRR(SUSY.MsD2 + Md * Md + matrix<complex>::Id(3)*(-1. / 3. * SUSY.s02()) * SUSY.getMz() * SUSY.getMz() * cos(2. * atan(SUSY.getTanb())));
    matrix<complex> dLR(SUSY.v1() / sqrt(2.) * SUSY.GetTD().hconjugate() - SUSY.getMuH() * Md * SUSY.getTanb());
    
    matrix<complex> MUP2(6,6,0.);
    matrix<complex> MDOWN2(6,6,0.);
    
    
    for(int i = 0; i < 3; i++)
        for(int j = 0; j < 3; j++){
            MUP2.assign(i,j,uLL(i,j));
            MUP2.assign(i+3,j,uLR(i,j));
            MUP2.assign(i,j+3,uLR(j,i).conjugate());
            MUP2.assign(i+3,j+3,uRR(i,j));
            
            MDOWN2.assign(i,j,dLL(i,j));
            MDOWN2.assign(i+3,j,dLR(i,j));
            MDOWN2.assign(i,j+3,dLR(j,i).conjugate());
            MDOWN2.assign(i+3,j+3,dRR(i,j));
        }
   
    // test lines
    
    Sflav_Spectrum("\"Up-squark LR mixing matrix, real part:\"");
    Sflav_Spectrum("\"Up-squark LR mixing matrix, imaginary part:\"");
    std::cout << "Matrice uLR = " << uLR << std::endl;
    
    Sflav_Spectrum("\"Down-squark LR mixing matrix, real part:\"");
    Sflav_Spectrum("\"Down-squark LR mixing matrix, imaginary part:\"");
    std::cout << "Matrice dLR = " << dLR << std::endl;
    
    

    
    // end lines
    
    MUP2.eigensystem(SUSY.Ru,SUSY.Msu2);
    MDOWN2.eigensystem(SUSY.Rd,SUSY.Msd2);
    
    
    // Neutralino masses
    
    
    matrix<complex> MN(4,4,0.);
    MN.assign(0,0,SUSY.GetM1());
    MN.assign(1,1,SUSY.GetM2());
    MN.assign(0,2,-SUSY.getCosb()*sqrt(SUSY.s02())*SUSY.getMz());
    MN.assign(2,0,MN(0,2));
    MN.assign(0,3,SUSY.getSinb()*sqrt(SUSY.s02())*SUSY.getMz());
    MN.assign(3,0,MN(0,3));
    MN.assign(1,2,SUSY.getCosb()*sqrt(SUSY.c02())*SUSY.getMz());
    MN.assign(2,1,MN(1,2));
    MN.assign(1,3,-SUSY.getSinb()*sqrt(SUSY.c02()) * SUSY.getMz());
    MN.assign(3,1,MN(1,3));
    MN.assign(2,3,-SUSY.getMuH());
    MN.assign(3,2,MN(2,3));
    
    //std::cout << "MN = " << MN << std::endl;
    //MN.eigensystem(SUSY.N,SUSY.Mneu);
    matrix<complex> Ntemp(4,4,0.);
    Ntemp = SUSY.N.transpose();
    MN.singularvalue(Ntemp, SUSY.N , SUSY.Mneu);
    
    
    
    matrix<complex> MC(2,2,0.);
    MC.assign(0,0,SUSY.GetM2());
    MC.assign(0,1,sqrt(2) * SUSY.getSinb() * SUSY.getMz());
    MC.assign(1,0,sqrt(2) * SUSY.getCosb()*SUSY.getMz());
    MC.assign(1,1,SUSY.getMuH());
    
    vector<double> M2Chi(2,0);
    (MC.hconjugate() * MC).eigensystem(SUSY.V,M2Chi);
    SUSY.Mch(0) = sqrt(M2Chi(0));
    SUSY.Mch(1) = sqrt(M2Chi(1));
    
    (MC * MC.hconjugate()).eigensystem(SUSY.U,M2Chi);
    
    
    SUSY.mh[2] = sqrt(SUSY.mHptree*SUSY.mHptree - SUSY.Mw()*SUSY.Mw());
    double temp = SUSY.mh[2]*SUSY.mh[2] + SUSY.getMz() * SUSY.getMz();
    double temp1 = SUSY.mh[2] * SUSY.getMz();
    SUSY.mh[0] = sqrt(0.5 * ( temp - sqrt(temp * temp - 4. * temp1 * temp1 * cos(2. * atan(SUSY.getTanb())) * cos(2. * atan(SUSY.getTanb())) ) ) );
    SUSY.mh[1] = sqrt(0.5 * ( temp + sqrt(temp * temp - 4. * temp1 * temp1 * cos(2. * atan(SUSY.getTanb())) * cos(2. * atan(SUSY.getTanb())) ) ) );
    SUSY.mHp = SUSY.mHptree;
   
    
}

void Spectrum::Sflav_Spectrum(const std::string block){
    
        
    const std::string com1 = "grep -r -A3 ";
    const std::string file = "mssm_data.txt";
    std::string com;
    
    
    com.append(com1);
    com.append(block);
    com.append(" ");
    com.append(file);
    com.append(" \n ");
    
    const char * c = com.c_str();
    
    system(c);
    
}
void Spectrum::Sflav_init(SUSY& SUSY){
    
   
    system("rm GeneralSUSY_NEW.conf");
    std::ofstream myfile("GeneralSUSY_NEW.conf");
     if (myfile.is_open())
    {
         
    int i,j;
    
    
    for( i = 0; i < 3; i++){
        for( j = 0; j < 3; j++){
            
           
    const std::string temp = "TU_";
    const std::string temp1 = "TD_";
    std::string TUr;
    std::string TUi;
    std::string TDr;
    std::string TDi;
    TUr.append(temp);
    TDr.append(temp1);
    std::ostringstream ss;
    ss << i + 1;
    ss << j + 1;
    TUr.append(ss.str());
    TDr.append(ss.str());
    TUi.append(TUr);
    TDi.append(TDr);
    TUi.append("i");
    TUr.append("r");
    TDr.append("r");
    TDi.append("i");
    
    // test - lines
    
    //std::cout << "TUr(" << i+1 << "," << j + 1 << ")= " << TUr << std::endl;
    //std::cout << "TUi(" << i+1 << "," << j + 1 << ")= " << TUi << std::endl;
    
    // end - test
    
    myfile << TUr;
    myfile << "  ";
    myfile << SUSY.TU(i,j).real();
    myfile << " 0.  0.  "; 
    myfile << "\n" ;
    myfile << TUi;
    myfile << "  ";
    myfile << SUSY.TU(i,j).imag();
    myfile << " 0.  0.  ";
    myfile << "\n" ;
    
    myfile << TDr;
    myfile << "  ";
    myfile << SUSY.TD(i,j).real();
    myfile << " 0.  0.  ";
    myfile << "\n" ;
    myfile << TDi;
    myfile << "  ";
    myfile << SUSY.TD(i,j).imag();
    myfile << " 0.  0.  ";
    myfile << "\n" ;
    
    
    const std::string temp2 = "msQ2_";
    std::string msQ2r;
    std::string msQ2i;
    
    msQ2r.append(temp2);
    msQ2r.append(ss.str());
    msQ2i.append(msQ2r);
    msQ2r.append("r");
    msQ2i.append("i");
    
    myfile << msQ2r;
    myfile << "  ";
    myfile << SUSY.MsQ2(i,j).real();
    myfile << " 0.  0.  ";
    myfile << "\n" ;
    myfile << msQ2i;
    myfile << "  ";
    myfile << SUSY.MsQ2(i,j).imag();
    myfile << " 0.  0.  ";
    myfile << "\n" ;
    
    
    const std::string temp3 = "msU2_";
    std::string msU2r;
    std::string msU2i;
    
    msU2r.append(temp3);
    msU2r.append(ss.str());
    msU2i.append(msU2r);
    msU2r.append("r");
    msU2i.append("i");
    
    myfile << msU2r;
    myfile << "  ";
    myfile << SUSY.MsU2(i,j).real();
    myfile << " 0.  0.  ";
    myfile << "\n" ;
    myfile << msU2i;
    myfile << "  ";
    myfile << SUSY.MsU2(i,j).imag();
    myfile << " 0.  0.  ";
    myfile << "\n" ;
    
    const std::string temp4 = "msD2_";
    std::string msD2r;
    std::string msD2i;
    
    msD2r.append(temp4);
    msD2r.append(ss.str());
    msD2i.append(msD2r);
    msD2r.append("r");
    msD2i.append("i");
    
    myfile << msD2r;
    myfile << "  ";
    myfile << SUSY.MsD2(i,j).real();
    myfile << " 0.  0.  ";
    myfile << "\n" ;
    myfile << msD2i;
    myfile << "  ";
    myfile << SUSY.MsD2(i,j).imag();
    myfile << " 0.  0.  ";
    myfile << "\n" ;  
    
        }
    }
    
 
    myfile << "m1r  ";
    myfile << SUSY.m1.real();
    myfile << " 0.  0.  ";
    myfile << "\n" ;
    myfile << "m1i  ";
    myfile << SUSY.m1.imag();
    myfile << " 0.  0.  ";
    myfile << "\n" ;
    myfile << "m2r  ";
    myfile << SUSY.m2.real();
    myfile << " 0.  0.  ";
    myfile << "\n" ;
    myfile << "m2i  ";
    myfile << SUSY.m2.imag();
    myfile << " 0.  0.  ";
    myfile << "\n" ;
    myfile << "m3  ";    
    myfile << SUSY.m3;
    myfile << " 0.  0.  ";
    myfile << "\n" ;
    myfile << "muHr  ";
    myfile << SUSY.muH.real();
    myfile << " 0.  0.  ";
    myfile << "\n" ;
    myfile << "muHi  ";
    myfile << SUSY.muH.imag();
    myfile << " 0.  0.  ";
    myfile << "\n" ;
    myfile << "tanb  ";
    myfile << SUSY.getTanb();
    myfile << " 0.  0.  ";
    myfile << "\n" ;
    myfile << "etab  ";
    myfile << SUSY.GetEtab();
    myfile << " 0.  0.  ";
    myfile << "\n" ;
    myfile << "rhob  ";
    myfile << SUSY.GetRhob();
    myfile << " 0.  0.  ";
    myfile << "\n" ;
    myfile << "A  ";
    myfile << SUSY.GetA();
    myfile << " 0.  0.  ";
    myfile << "\n" ;
    myfile << "lambda  ";
    myfile << SUSY.GetLambda();
    myfile << " 0.  0.  ";
    myfile << "\n" ;
    myfile << "M_z  ";
    myfile << SUSY.getMz();
    myfile << " 0.  0.  ";
    myfile << "\n" ;
    myfile << "GF  ";
    myfile << SUSY.getGF();
    myfile << " 0.  0.  ";
    myfile << "\n" ;
    
    
    
    myfile.close();
    }
    else std::cout << "Unable to open GeneralSUSY_NEW.conf";
 
    
      system("./../susy_flavor/sflav \n");


    
}

