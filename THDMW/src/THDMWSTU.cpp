/*
 * Copyright (C) 2018 HEPfit Collaboration
 *
 *
 * For the licensing terms see doc/COPYING.
 */


#include "THDMWSTU.h"
#include "THDMW.h"
#include "THDMWcache.h"

//Theoretical formulae taken from https://arxiv.org/pdf/0907.2696.pdf

THDMWSTU::THDMWSTU(const StandardModel& SM_i)
: ThObservable(SM_i), myTHDMW(static_cast<const THDMW&> (SM_i))

{
    mycache = new THDMWcache(SM_i);
};

double THDMWSTU::computeThValue()
{
    return 0.0;
}

double THDMWSTU::F(const double m02, const double m12) const {
    double F;

    if(m02 == 0. && m12 != 0.) {
        F= m12;
    } else if(m02 != 0. && m12 == 0.){
        F=m02;
    } else if((m02 == 0. && m12 == 0.) || (fabs(m02-m12) < LEPS)){
        F=0.;
    } else if (m02 != 0 && m12 != 0){
        F=(m02 + m12) - 2*((m02 * m12) / (m02 - m12)) * log(m02/m12);
    } else
        throw std::runtime_error("Error in THDM::F()");
    return (F);
}







THDMWDeltaS::THDMWDeltaS(const StandardModel& SM_i)
: THDMWSTU(SM_i)
{}

double THDMWDeltaS::computeThValue()
{
    double mSp2=myTHDMW.getMyTHDMWCache()->mSpsq;
    double mSR2=myTHDMW.getMyTHDMWCache()->mSRsq;
    double mSI2=myTHDMW.getMyTHDMWCache()->mSIsq;
    double s_W2 = myTHDMW.sW2();
    //std::cout<< "SinW_2=" << s_W2 << std::endl;//Why this value?
    double c_W2 = 1-myTHDMW.sW2();
    double MZ2=pow(myTHDMW.getMz(),2);
    //double MW2=pow(myTHDMW->Mw(),2);
    //We define the self energies removing the coupling because it will cancel
    
    
    
    //gslpp::complex delPiWW_MZ;
    gslpp::complex delPiZZ_MZ;
    //gslpp::complex delPiWW_0;
    gslpp::complex delPiZZ_0;
    gslpp::complex delPiAAprime_0;
    gslpp::complex delPiZAprime_0;
    //gslpp::complex delPiAA_MZ;
    //gslpp::complex delPiZA_MZ;
    
    //g1 on arxiv:0907.2696 is e/s_W
    //delPiWW_MZ=(1./(2.*pow(M_PI,2)*s_W2))*(mycache->B00_MZ2_MZ2_mSi2_mSp2(MZ2,mSI2,mSp2)+mycache->B00_MZ2_MZ2_mSr2_mSp2(MZ2,mSR2,mSp2)-0.5*mycache->A0_MZ2_mSp2(MZ2,mSp2)-0.25*mycache->A0_MZ2_mSp2(MZ2,mSR2)-0.25*mycache->A0_MZ2_mSp2(MZ2,mSI2));
    
        
    //delPiWW_0=(1./(8.*pow(M_PI,2)*s_W2))*(0.5*F(mSp2,mSR2)+0.5*F(mSp2,mSI2));
    
    
    //delPiAA_MZ=(2./pow(M_PI,2))*(mycache->B00_MZ2_MZ2_mSp2_mSp2(MZ2,mSp2)-0.5*mycache->A0_MZ2_mSp2(MZ2,mSp2));
    //delPiZA_MZ=((1.-2.*s_W2)/(pow(M_PI,2)*sqrt(c_W2)*sqrt(s_W2)))*(mycache->B00_MZ2_MZ2_mSp2_mSp2(MZ2,mSp2)-0.5*mycache->A0_MZ2_mSp2(MZ2,mSp2));

    
    //double nu2=myTHDMW->getTHDMW_nu2();
    //double nu3=myTHDMW->getTHDMW_nu3();
    //double v=246.2;
    
    //std::cout<< "prove10=" << delPiAAprime_0 << std::endl;
    //std::cout<< "prove11=" << delPiAA_MZ/MZ2 << std::endl;
    //std::cout<< "prove20=" << delPiZAprime_0 << std::endl;
    //std::cout<< "prove21=" << delPiZA_MZ/MZ2 << std::endl;
    //std::cout<< "DeltaS approx.=" << pow(v,2)*nu2/(6.*M_PI*mSp2) << std::endl;
    //std::cout<< "DeltaS old=" << (16.*M_PI*s_W2*c_W2)*(((delPiZZ_MZ.real()-delPiZZ_0.real())/MZ2)-(c_W2-s_W2)*delPiZAprime_0.real()/(sqrt(s_W2)*sqrt(c_W2))-delPiAAprime_0.real()) << std::endl;

    //std::cout<< "DeltaS fraction old=" << (16.*M_PI*s_W2*c_W2)*(((delPiZZ_MZ.real()-delPiZZ_0.real())/MZ2)-(c_W2-s_W2)*delPiZAprime_0.real()/(sqrt(s_W2)*sqrt(c_W2))-delPiAAprime_0.real())/(pow(v,2)*nu2/(6.*M_PI*mSp2)) << std::endl;
    //std::cout<< "DeltaS fraction new=" << ((16.*M_PI*s_W2*c_W2)*(((delPiZZ_MZ.real()-delPiZZ_0.real())/MZ2)-(c_W2-s_W2)*delPiZA_MZ.real()/(sqrt(s_W2)*sqrt(c_W2)*MZ2)-delPiAA_MZ.real()/MZ2))/(pow(v,2)*nu2/(6.*M_PI*mSp2)) << std::endl;

    
    //std::cout<< "B00_MZ2_MZ2_mSr2_mSi2(MZ2,mSR2,mSI2)=" << mycache->B00_MZ2_MZ2_mSr2_mSi2(MZ2,mSR2,mSI2) << std::endl;
    //std::cout<< "A0_MZ2_mSp2=" << mycache->A0_MZ2_mSp2(MZ2,mSp2) << std::endl;
    //std::cout<< "mSp2=" << mSp2 << std::endl;
    //std::cout<< "mSr2=" << mSR2 << std::endl;
    //std::cout<< "MZ2=" << MZ2 << std::endl;
    //std::cout<< "THDMWpositiveMassSquares=" << myTHDMW->getMyTHDMWCache()->setOtherParameters() << std::endl;
    //std::cout<< "THDMWpositiveMassSquares=" << THDMWpositiveMassSquares << std::endl;
    if(fabs(myTHDMW.getMyTHDMWCache()->setOtherParameters()-1.)<1.e-10)
    {
        delPiZZ_MZ=(1./(2.*pow(M_PI,2)*c_W2*s_W2))*(pow(1.-2.*s_W2,2)*(mycache->B00_MZ2_MZ2_mSp2_mSp2(MZ2,mSp2)-0.5*mycache->A0_MZ2_mSp2(MZ2,mSp2))+mycache->B00_MZ2_MZ2_mSr2_mSi2(MZ2,mSR2,mSI2)-0.25*mycache->A0_MZ2_mSp2(MZ2,mSR2)-0.25*mycache->A0_MZ2_mSp2(MZ2,mSI2));    
        delPiZZ_0=(1./(8.*pow(M_PI,2)*s_W2*c_W2))*(0.5*F(mSI2,mSR2));
        delPiAAprime_0=-(1/(6.*pow(M_PI,2)))*mycache->B0_MZ2_0_mSp2_mSp2(MZ2,mSp2);
        delPiZAprime_0=-((1-2*s_W2)/(12.*pow(M_PI,2)*sqrt(c_W2)*sqrt(s_W2)))*mycache->B0_MZ2_0_mSp2_mSp2(MZ2,mSp2);

        return (16.*M_PI*s_W2*c_W2)*(((delPiZZ_MZ.real()-delPiZZ_0.real())/MZ2)-(c_W2-s_W2)*delPiZAprime_0.real()/(sqrt(s_W2)*sqrt(c_W2))-delPiAAprime_0.real());

    }
    else
    {
        return std::numeric_limits<double>::quiet_NaN();
    }
    }

THDMWDeltaT::THDMWDeltaT(const StandardModel& SM_i)
: THDMWSTU(SM_i)
{}

double THDMWDeltaT::computeThValue()
{
    double mSp2=myTHDMW.getMyTHDMWCache()->mSpsq;
    double mSR2=myTHDMW.getMyTHDMWCache()->mSRsq;
    double mSI2=myTHDMW.getMyTHDMWCache()->mSIsq;
    double s_W2 = myTHDMW.sW2();
    //std::cout<< "SinW_2=" << s_W2 << std::endl;//Why this value?
    double c_W2 = 1-myTHDMW.sW2();
    double MZ2=pow(myTHDMW.getMz(),2);
    double MW2=pow(myTHDMW.Mw(),2);
    //We define the self energies removing the coupling because it will cancel
    

    gslpp::complex delPiWW_0;
    gslpp::complex delPiZZ_0;

    
    //g1 on arxiv:0907.2696 is e/s_W
    
    
    //double nu2=myTHDMW->getTHDMW_nu2();
    //double nu3=myTHDMW->getTHDMW_nu3();
    //double v=246.2;
    //std::cout<< "DeltaT approx.=" << pow(v,4)*(pow(nu2,2)-pow(2.*nu3,2))/(96.*pow(M_PI,1)*mSp2*MW2*s_W2) << std::endl;//typo in arxiv:0907.269, look at arxiv:0606172
    //std::cout<< "DeltaT fraction=" << pow(v,4)*(pow(nu2,2)-pow(2.*nu3,2))/(96.*pow(M_PI,1)*mSp2*MW2*s_W2)/((4.*M_PI)*(delPiWW_0.real()/MW2-delPiZZ_0.real()/MZ2)) << std::endl;
    if(fabs(myTHDMW.getMyTHDMWCache()->setOtherParameters()-1.)<1.e-10)
    {
        delPiWW_0=(1./(8.*pow(M_PI,2)*s_W2))*(0.5*F(mSp2,mSR2)+0.5*F(mSp2,mSI2));
        delPiZZ_0=(1./(8.*pow(M_PI,2)*s_W2*c_W2))*(0.5*F(mSI2,mSR2));
        return (4.*M_PI)*(delPiWW_0.real()/MW2-delPiZZ_0.real()/MZ2);
    }
    else
    {
        return std::numeric_limits<double>::quiet_NaN();
    }
    
    
    
    
    
}


THDMWDeltaU::THDMWDeltaU(const StandardModel& SM_i)
: THDMWSTU(SM_i)
{}

double THDMWDeltaU::computeThValue()
{
    double mSp2=myTHDMW.getMyTHDMWCache()->mSpsq;
    double mSR2=myTHDMW.getMyTHDMWCache()->mSRsq;
    double mSI2=myTHDMW.getMyTHDMWCache()->mSIsq;
    double s_W2 = myTHDMW.sW2();
    //std::cout<< "SinW_2=" << s_W2 << std::endl;//Why this value?
    double c_W2 = 1-myTHDMW.sW2();
    double MZ2=pow(myTHDMW.getMz(),2);
    double MW2=pow(myTHDMW.Mw(),2);
    //We define the self energies removing the coupling because it will cancel
    
    gslpp::complex delPiWW_MW;
    gslpp::complex delPiZZ_MZ;
    gslpp::complex delPiWW_0;
    gslpp::complex delPiZZ_0;
    gslpp::complex delPiAAprime_0;
    gslpp::complex delPiZAprime_0;
    //gslpp::complex delPiAA_MZ;
    //gslpp::complex delPiZA_MZ;
    //gslpp::complex delPiAA_MW;
    //gslpp::complex delPiZA_MW;
    
    
    
    //gslpp::complex S;
    
    //g1 on arxiv:0907.2696 is e/s_W
    
    
    //delPiAA_MZ=(2./pow(M_PI,2))*(mycache->B00_MZ2_MZ2_mSp2_mSp2(MZ2,mSp2)-0.5*mycache->A0_MZ2_mSp2(MZ2,mSp2));
    //delPiZA_MZ=((1.-2.*s_W2)/(pow(M_PI,2)*sqrt(c_W2)*sqrt(s_W2)))*(mycache->B00_MZ2_MZ2_mSp2_mSp2(MZ2,mSp2)-0.5*mycache->A0_MZ2_mSp2(MZ2,mSp2));
    //delPiAA_MW=(2./pow(M_PI,2))*(mycache->B00_MZ2_MZ2_mSp2_mSp2(MW2,mSp2)-0.5*mycache->A0_MZ2_mSp2(MW2,mSp2));
    //delPiZA_MW=((1.-2.*s_W2)/(pow(M_PI,2)*sqrt(c_W2)*sqrt(s_W2)))*(mycache->B00_MZ2_MZ2_mSp2_mSp2(MW2,mSp2)-0.5*mycache->A0_MZ2_mSp2(MW2,mSp2));


    //S=(16.*M_PI*s_W2*c_W2)*(((delPiZZ_MZ.real()-delPiZZ_0.real())/MZ2)-(c_W2-s_W2)*delPiZA_MZ.real()/(sqrt(s_W2)*sqrt(c_W2)*MZ2)-delPiAA_MZ.real()/MZ2);

    //std::cout<< "U new=" << (16.*M_PI*s_W2)*((delPiWW_MW.real()-delPiWW_0.real())/MW2-(sqrt(c_W2)/sqrt(s_W2))*(delPiZA_MZ/MZ2)-delPiAA_MZ/MZ2)-c_W2*S << std::endl;
    //std::cout<< "U new MW=" << (16.*M_PI*s_W2)*((delPiWW_MW.real()-delPiWW_0.real())/MW2-(sqrt(c_W2)/sqrt(s_W2))*(delPiZA_MW/MW2)-delPiAA_MW/MW2)-c_W2*S << std::endl;
    //std::cout<< "U new mix=" << (16.*M_PI*s_W2)*((delPiWW_MW.real()-delPiWW_0.real())/MW2-(sqrt(c_W2)/sqrt(s_W2))*(delPiZA_MW/MZ2)-delPiAA_MW/MZ2)-c_W2*S << std::endl;
    //std::cout<< "Delta MW=" << (delPiWW_MW.real()-delPiWW_0.real())/MW2 << std::endl;
    //std::cout<< "Term with delMZ=" << c_W2*(delPiZZ_MZ.real()-delPiZZ_0.real())/MZ2 << std::endl;
    //std::cout<< "delPiAA=" << s_W2*delPiAA_MZ.real()/MZ2 << std::endl;
    //std::cout<< "delPiAZ=" << 2*sqrt(s_W2)*sqrt(c_W2)*delPiZA_MZ.real()/MZ2 << std::endl;
    
    if(fabs(myTHDMW.getMyTHDMWCache()->setOtherParameters()-1.)<1.e-10)
    {
        delPiWW_MW=(1./(2.*pow(M_PI,2)*s_W2))*(mycache->B00_MZ2_MZ2_mSi2_mSp2(MW2,mSI2,mSp2)+mycache->B00_MZ2_MZ2_mSr2_mSp2(MW2,mSR2,mSp2)-0.5*mycache->A0_MZ2_mSp2(MW2,mSp2)-0.25*mycache->A0_MZ2_mSp2(MW2,mSR2)-0.25*mycache->A0_MZ2_mSp2(MW2,mSI2));
        delPiZZ_MZ=(1./(2.*pow(M_PI,2)*c_W2*s_W2))*(pow(1.-2.*s_W2,2)*(mycache->B00_MZ2_MZ2_mSp2_mSp2(MZ2,mSp2)-0.5*mycache->A0_MZ2_mSp2(MZ2,mSp2))+mycache->B00_MZ2_MZ2_mSr2_mSi2(MZ2,mSR2,mSI2)-0.25*mycache->A0_MZ2_mSp2(MZ2,mSR2)-0.25*mycache->A0_MZ2_mSp2(MZ2,mSI2));    
    
        delPiWW_0=(1./(8.*pow(M_PI,2)*s_W2))*(0.5*F(mSp2,mSR2)+0.5*F(mSp2,mSI2));
        delPiZZ_0=(1./(8.*pow(M_PI,2)*s_W2*c_W2))*(0.5*F(mSI2,mSR2));
        delPiAAprime_0=-(1/(6.*pow(M_PI,2)))*mycache->B0_MZ2_0_mSp2_mSp2(MZ2,mSp2);
        delPiZAprime_0=-((1-2*s_W2)/(12.*pow(M_PI,2)*sqrt(c_W2)*sqrt(s_W2)))*mycache->B0_MZ2_0_mSp2_mSp2(MZ2,mSp2);
        return (16.*M_PI*s_W2)*((delPiWW_MW.real()-delPiWW_0.real())/MW2-c_W2*((delPiZZ_MZ.real()-delPiZZ_0.real())/MZ2)-(s_W2*delPiAAprime_0.real())-(2.*sqrt(s_W2)*sqrt(c_W2)*delPiZAprime_0.real()));

    }
    else
    {
    
        return std::numeric_limits<double>::quiet_NaN();
    }
    
    
    
    }