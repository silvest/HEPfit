/* 
 * Copyright (C) 2018 HEPfit Collaboration
 *
 *
 * For the licensing terms see doc/COPYING.
 */

#include "NPSMEFTd6GeneralMatching.h"
#include "NPSMEFTd6General.h"
#include <stdexcept>

NPSMEFTd6GeneralMatching::NPSMEFTd6GeneralMatching(const NPSMEFTd6General & NPSMEFTd6General_i) :

    StandardModelMatching(NPSMEFTd6General_i),
    mySMEFT(NPSMEFTd6General_i), VuL(3,0.), VuR(3,0.), VdL(3,0.), VdR(3,0.), VeL(3,0.), VeR(3,0.)
{}

void NPSMEFTd6GeneralMatching::updateLEFTGeneralParameters()
{
        
    // Dimension 6 operators with no flavour index are assigned directly here
    
    LambdaNP2=mySMEFT.getLambda_NP()*mySMEFT.getLambda_NP();
    v = mySMEFT.v();
    v2 = v*v;
    double vosq2 = v / sqrt(2.);
    
    CG = mySMEFT.getSMEFTCoeffEW("CG")*LambdaNP2;  
    CW = mySMEFT.getSMEFTCoeffEW("CW")*LambdaNP2; 
    CHG = mySMEFT.getSMEFTCoeffEW("CHG")*LambdaNP2;  
    CHW = mySMEFT.getSMEFTCoeffEW("CHW")*LambdaNP2;  
    CHB = mySMEFT.getSMEFTCoeffEW("CHB")*LambdaNP2;  
    CHWB = mySMEFT.getSMEFTCoeffEW("CHWB")*LambdaNP2;  
    CHD = mySMEFT.getSMEFTCoeffEW("CHD")*LambdaNP2;  
    CHbox = mySMEFT.getSMEFTCoeffEW("CHbox")*LambdaNP2;  
    CH = mySMEFT.getSMEFTCoeffEW("CH")*LambdaNP2;  
    CGtilde = mySMEFT.getSMEFTCoeffEW("CGtilde")*LambdaNP2;  
    CWtilde = mySMEFT.getSMEFTCoeffEW("CWtilde")*LambdaNP2;  
    CHGtilde = mySMEFT.getSMEFTCoeffEW("CHGtilde")*LambdaNP2;  
    CHWtilde = mySMEFT.getSMEFTCoeffEW("CHWtilde")*LambdaNP2;  
    CHBtilde = mySMEFT.getSMEFTCoeffEW("CHBtilde")*LambdaNP2;  
    CHWtildeB = mySMEFT.getSMEFTCoeffEW("CHWtildeB")*LambdaNP2;  
    
    //Now we do not use the SILH basis anymore, we'll set these operators to zero 
    C2B = 0.; 
    C2W = 0.; 
    C2BS = 0.; 
    C2WS = 0.; 
    CDHB = 0.; 
    CDHW = 0.; 
    CDB = 0.; 
    CDW = 0.; 
    CT = 0.;
    
    //For operators with fermionic indices we need to switch to the mass eigenstate basis 
    
    // Let us first define the full mass matrices, including the effect of dimension six operators
    
    gslpp::matrix<complex> MU(3,0.), MD(3,0.), ME(3,0.);
    
    for (int i = 0; i < 3; i++)
        for (int j = 0; j < 3; j++) 
        {
        MU.assignre(i,j,vosq2 * (mySMEFT.getSMEFTCoeffEW("YuR",i,j) * (1 + mySMEFT.getDelta_v()) + mySMEFT.getSMEFTCoeffEW("CuHR",i,j)*v2));
        MU.assignim(i,j,vosq2 * (mySMEFT.getSMEFTCoeffEW("YuI",i,j) * (1 + mySMEFT.getDelta_v()) + mySMEFT.getSMEFTCoeffEW("CuHI",i,j)*v2));
        MD.assignre(i,j,vosq2 * (mySMEFT.getSMEFTCoeffEW("YdR",i,j) * (1 + mySMEFT.getDelta_v()) + mySMEFT.getSMEFTCoeffEW("CdHR",i,j)*v2));
        MD.assignim(i,j,vosq2 * (mySMEFT.getSMEFTCoeffEW("YdI",i,j) * (1 + mySMEFT.getDelta_v()) + mySMEFT.getSMEFTCoeffEW("CdHI",i,j)*v2));
        ME.assignre(i,j,vosq2 * (mySMEFT.getSMEFTCoeffEW("YeR",i,j) * (1 + mySMEFT.getDelta_v()) + mySMEFT.getSMEFTCoeffEW("CeHR",i,j)*v2));
        ME.assignim(i,j,vosq2 * (mySMEFT.getSMEFTCoeffEW("YeI",i,j) * (1 + mySMEFT.getDelta_v()) + mySMEFT.getSMEFTCoeffEW("CeHI",i,j)*v2));
       }
    
    gslpp::vector<double> m2(3);
    
    MU.singularvalue(VuR, VuL, m2);
    mySMEFT.getQuarks(QCD::UP).setMass(mySMEFT.Mrun(mySMEFT.getQuarks(QCD::UP).getMass_scale(),mySMEFT.getMuw(),sqrt(m2(0))));
    mySMEFT.getQuarks(QCD::CHARM).setMass(mySMEFT.Mrun(mySMEFT.getQuarks(QCD::CHARM).getMass_scale(),mySMEFT.getMuw(),sqrt(m2(1))));
    mySMEFT.getQuarks(QCD::TOP).setMass(mySMEFT.Mrun(mySMEFT.getQuarks(QCD::TOP).getMass_scale(),mySMEFT.getMuw(),sqrt(m2(2))));

    MD.singularvalue(VdR, VdL, m2);
    
    mySMEFT.getQuarks(QCD::DOWN).setMass(mySMEFT.Mrun(mySMEFT.getQuarks(QCD::DOWN).getMass_scale(),mySMEFT.getMuw(),sqrt(m2(0))));
    mySMEFT.getQuarks(QCD::STRANGE).setMass(mySMEFT.Mrun(mySMEFT.getQuarks(QCD::STRANGE).getMass_scale(),mySMEFT.getMuw(),sqrt(m2(1))));
    mySMEFT.getQuarks(QCD::BOTTOM).setMass(mySMEFT.Mrun(mySMEFT.getQuarks(QCD::BOTTOM).getMass_scale(),mySMEFT.getMuw(),sqrt(m2(2))));
    
    ME.singularvalue(VeR, VeL, m2);
    
    mySMEFT.getLeptons(QCD::ELECTRON).setMass(sqrt(m2(0)));
    mySMEFT.getLeptons(QCD::MU).setMass(sqrt(m2(1)));
    mySMEFT.getLeptons(QCD::TAU).setMass(sqrt(m2(2)));
    
    //Computing the CKM 
    gslpp::matrix<complex> CKMUnphys = (VuL.hconjugate()) * VdL;

    //std::cout << "CKM unphys = " << CKMUnphys << std::endl;
    
    mySMEFT.getCKM().computeCKM(CKMUnphys(0,1).abs(),CKMUnphys(1,2).abs(),CKMUnphys(0,2).abs(),(-CKMUnphys(0,0)*CKMUnphys(0,2).conjugate()/(CKMUnphys(1,0)*CKMUnphys(1,2).conjugate())).arg());
    
    //std::cout << "computed CKM = " << mySMEFT.getCKM().getCKM() << std::endl;
    
    double a11 = remainder(CKMUnphys(0, 0).arg() - mySMEFT.getCKM().getV_ud().arg(), 2.*M_PI);
    double a12 = remainder(CKMUnphys(0, 1).arg() - mySMEFT.getCKM().getV_us().arg(), 2.*M_PI);
    double a13 = remainder(CKMUnphys(0, 2).arg() - mySMEFT.getCKM().getV_ub().arg(), 2.*M_PI);


    //    double a23 = (gslpp::complex(CKMUnphys(1, 0) / CKM(1, 0))).arg() - a11 + a13;
    //    double a33 = (gslpp::complex(CKMUnphys(2, 0) / CKM(2, 0))).arg() - a11 + a13;
    double a23 = remainder(CKMUnphys(1, 0).arg() - mySMEFT.getCKM().getV_cd().arg(), 2.*M_PI) - a11 + a13;
    double a33 = remainder(CKMUnphys(2, 0).arg() - mySMEFT.getCKM().getV_td().arg(), 2.*M_PI) - a11 + a13;

    gslpp::matrix<gslpp::complex> phi1(3, 3, 0.);
    phi1.assign(0, 0, 1.);
    phi1.assign(1, 1, gslpp::complex(1., a23 - a13, true));
    phi1.assign(2, 2, gslpp::complex(1., a33 - a13, true));

    gslpp::matrix<gslpp::complex> phi2dag(3, 3, 0.);
    phi2dag.assign(0, 0, gslpp::complex(1., -a11, true));
    phi2dag.assign(1, 1, gslpp::complex(1., -a12, true));
    phi2dag.assign(2, 2, gslpp::complex(1., -a13, true));

    gslpp::matrix<gslpp::complex> phie(3, 3, 0.);
    phie.assign(0, 0, gslpp::complex(1., -(VeR(0, 0)).arg(), true));
    phie.assign(1, 1, gslpp::complex(1., -(VeR(1, 1)).arg(), true));
    phie.assign(2, 2, gslpp::complex(1., -(VeR(2, 2)).arg(), true));

    VeR = VeR*phie;
    VeL = VeL*phie;
    VuL = VuL*phi1;
    VuR = VuR*phi1;
    VdL = VdL*phi2dag;
    VdR = VdR*phi2dag;

    //std::cout << "CKM from rotated UfA = " << (VuL.hconjugate()) * VdL << std::endl;

    //std::cout << "has the diagonalization worked? " << VuR.hconjugate()*MU*VuL << std::endl;
    //std::cout << "has the diagonalization worked? " << VdR.hconjugate()*MD*VdL << std::endl;
    //std::cout << "has the diagonalization worked? " << VeR.hconjugate()*ME*VeL << std::endl;
    
    //rotate the WC
    
    
    
    StandardModelMatching::updateSMParameters();
}

NPSMEFTd6GeneralMatching::~NPSMEFTd6GeneralMatching()
{}
