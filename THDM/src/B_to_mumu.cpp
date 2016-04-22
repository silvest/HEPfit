/* 
 * Copyright (C) 2015 SusyFit Collaboration
 * All rights reserved.
 *
 * For the licensing terms see doc/COPYING.
 */

#include "B_to_mumu.h"
#include "StandardModel.h"
#include "gslpp_complex.h"

B_to_mumu::B_to_mumu(const StandardModel& SM_i): 

        ThObservable(SM_i), 
        myTHDM(static_cast<const THDM*> (&SM_i))
{

}

B_to_mumu::~B_to_mumu()
{}

double B_to_mumu::computeThValue()
{
   return 0.;
}

double B_to_mumu::computeC10()
{
    double mt=myTHDM->getQuarks(QCD::TOP).getMass();
    double MW=myTHDM->Mw_tree();
    double xt=mt*mt/MW/MW;
    double tanb=myTHDM->gettanb();
    double mHp2=myTHDM->getmHp2();
    double xiuA=-1.0/tanb;
    double xHp=mHp2/(MW*MW);

    double C10ZSM = -(xt*(6.0-7.0*xt+xt*xt+(2.0+3.0*xt)*log(xt))) / (8.0*(xt-1.0)*(xt-1.0));
    double C10boxSM = (xt*(1.0-xt+log(xt))) / (4.0*(xt-1.0)*(xt-1.0));
    //C10ZSM+C10boxSM=Y0
    double C10ZTHDM = xiuA*xiuA*xt*xt/8.0*(1.0/(xHp-xt)-xHp*log(xHp/xt)/((xHp-xt)*(xHp-xt)));

    return C10ZSM+C10boxSM+C10ZTHDM;
}

double B_to_mumu::computeCP()
{
    double mt=myTHDM->getQuarks(QCD::TOP).getMass();
    double MW=myTHDM->Mw_tree();
    double sw2=myTHDM->s02();
    double mA2=myTHDM->getmA2();
    double mHp2=myTHDM->getmHp2();
    double tanb=myTHDM->gettanb();
    std::string modelflag=myTHDM->getModelTypeflag();
    double xt=mt*mt/(MW*MW);
    double xA=mA2/(MW*MW);
    double xHp=mHp2/(MW*MW);
    double xiuA=-1.0/tanb;

    double xidA=0.0;
    double xilA=0.0;

    if( modelflag == "type1" ) {
        xidA=1.0/tanb;
        xilA=1.0/tanb;
    }
    else if( modelflag == "type2" ) {
        xidA=-tanb;
        xilA=-tanb;
    }
    else if( modelflag == "typeX" ) {
        xidA=1.0/tanb;
        xilA=-tanb;
    }
    else if( modelflag == "typeY" ) {
        xidA=-tanb;
        xilA=1.0/tanb;
    }
    else {
        throw std::runtime_error("modelflag can be only any of \"type1\", \"type2\", \"typeX\" or \"typeY\"");
    }

    double f1=(xt-xHp+xHp*log(xHp)-xt*log(xt))/(2.0*(xHp-xt));
    double f2=(xt-(xHp*xt*(log(xHp)-log(xt)))/(xHp-xt))/(2.0*(xHp-xt));
    double f3=(xHp-xHp*xHp*log(xHp)/(xHp-xt)+xt*(2.0*xHp-xt)*log(xt)/(xHp-xt))/(2.0*(xHp-xt));
    double f4=(0.5*xt*(3.0*xHp-xt)-(xHp*xHp*xt*(log(xHp)-log(xt)))/(xHp-xt))/(4.0*(xHp-xt)*(xHp-xt));
    double f5=(0.5*xt*(xHp-3.0*xt)-(xHp*xt*(xHp-2.0*xt)*(log(xHp)-log(xt)))/(xHp-xt))/(4.0*(xHp-xt)*(xHp-xt));
    double f6=(xt*(xt*xt-3.0*xHp*xt+9.0*xHp-5.0*xt-2.0))/(8.0*(xt-1.0)*(xt-1.0)*(xHp-xt))
              +(xHp*(xHp*xt-3.0*xHp+2.0*xt)*log(xHp))/(4.0*(xHp-1.0)*(xHp-xt)*(xHp-xt))
              +((xHp*xHp*(-2.0*xt*xt*xt+6.0*xt*xt-9.0*xt+2.0)+3.0*xHp*xt*xt*(xt*xt-2.0*xt+3.0)-xt*xt*(2.0*xt*xt*xt-3.0*xt*xt+3.0*xt+1.0))*log(xt))
               /(4.0*(xt-1.0)*(xt-1.0)*(xt-1.0)*(xHp-xt)*(xHp-xt));
    double f7=(xt*xt+xt-8.0)/(8.0*(xt-1.0)*(xt-1.0))-(xHp*(xHp+2.0)*log(xHp))/(4.0*(xHp-1.0)*(xHp-xt))
              +((xHp*(xt*xt*xt-3.0*xt*xt+3.0*xt+2.0)+3.0*xt*xt*(xt-2.0))*log(xt))/(4.0*(xt-1.0)*(xt-1.0)*(xt-1.0)*(xHp-xt));

    double g3a=-xidA*xidA*xiuA*f1+xidA*xiuA*xiuA*(f3-f2)-xiuA*xiuA*xiuA*(f5+f4)-xiuA*(f7+f6)+xidA*f1;

    double CPboxSM=xt*(71.0*xt*xt-172.0*xt-19.0)/(144.0*(xt-1.0)*(xt-1.0)*(xt-1.0))
                   +(xt*xt*xt*xt-12.0*xt*xt*xt+34.0*xt*xt-xt-2.0)*log(xt)/(24.0*pow(xt-1.0,4));
    double CPZSM=(xt*(18.0*xt*xt*xt-137.0*xt*xt+262.0*xt-95.0)/(6.0*(xt-1.0)*(xt-1.0)*(xt-1.0))
                  +(8.0*xt*xt*xt*xt-11.0*xt*xt*xt-15.0*xt*xt+12.0*xt-2.0)*log(xt)/pow(xt-1.0,4))/12.0
                 -sw2*(xt*(18.0*xt*xt*xt-139.0*xt*xt+274.0*xt-129.0)/(2.0*(xt-1.0)*(xt-1.0)*(xt-1.0))
                       +(24.0*xt*xt*xt*xt-33.0*xt*xt*xt-45.0*xt*xt+50.0*xt-8.0)*log(xt)/pow(xt-1.0,4))/36.0;
    double CPZTHDM=xt*(-xidA*xiuA*(-0.5*(xt+xHp)+xt*xHp*log(xHp/xt)/(xHp-xt))
                       +xiuA*xiuA*((xHp*xHp-8.0*xHp*xt-17.0*xt*xt)/(36.0*(xHp-xt))-xt*(xHp-xt)
                                   +(xt*xt*(3.0*xHp+xt)/(6.0*(xHp-xt)*(xHp-xt))+xt*xHp)*log(xHp/xt)))/(4.0*(xHp-xt)*(xHp-xt))
                   +sw2*xt*(-xidA*xiuA*((2.5*xt-1.5*xHp)+xHp*(2.0*xHp-3.0*xt)*log(xHp/xt)/(xHp-xt))
                            -xiuA*xiuA*(((4.0*xHp*xHp*xHp-12.0*xHp*xHp*xt+9.0*xHp*xt*xt+3.0*xt*xt*xt)/(6.0*(xHp-xt)*(xHp-xt))+1.5*xt*xHp)*log(xHp/xt)
                                        -(17.0*xHp*xHp-64.0*xHp*xt+71.0*xt*xt)/(36.0*(xHp-xt))-1.5*xt*(xHp-xt)))/(6.0*(xHp-xt)*(xHp-xt));
    double CPboxTHDM=xiuA*xilA*xt*(1.0+(2.0*xt*xt-xHp*xt-xHp*xHp)*log(xt)/((xt-1.0)*(xHp-xt))
                                   +xHp*(1.0-2.0*xt+xHp)*log(xHp)/((xHp-1.0)*(xHp-xt)))/(8.0*(xHp-xt))
                     +xidA*xilA*xt*log(xHp/xt)/(4.0*(xHp-xt));
    double CPphiTHDM=-xt*xilA*g3a/(2.0*xA);

    return CPboxSM+CPZSM+CPZTHDM+CPboxTHDM+CPphiTHDM;
}

double B_to_mumu::computeCS()
{
    double mt=myTHDM->getQuarks(QCD::TOP).getMass();
    double MW=myTHDM->Mw_tree();
    double mHh2=myTHDM->getmHh2();
    double mHl=myTHDM->getMHl();
    double mHp2=myTHDM->getmHp2();
    double m12_2=myTHDM->getm12_2();
    double tanb=myTHDM->gettanb();
    double sinb=tanb/sqrt(1.0+tanb*tanb);
    double cosb=1.0/sqrt(1.0+tanb*tanb);
    double bma=myTHDM->getbma();
    double sinbma=sin(bma);
    double cosbma=cos(bma);
    double sina=sinb*cos(bma)-cosb*sin(bma);
    double cosa=cosb*cos(bma)+sinb*sin(bma);
    std::string modelflag=myTHDM->getModelTypeflag();
    double xt=mt*mt/(MW*MW);
    double xh=mHl*mHl/(MW*MW);
    double xHh=mHh2/(MW*MW);
    double xHp=mHp2/(MW*MW);
    double xiuA=-1.0/tanb;

    //These correspond to -vev times the triple Higgs couplings g_hHpHm and g_HhHpHm from THDMquantities.
    double lambdahHpHm=((mHl*mHl-2.0*mHp2)*(cosbma*(cosb*cosb-sinb*sinb)-2.0*sinbma*sinb*cosb)
                        +(-4.0*m12_2/(sinb*cosb)+3.0*mHl*mHl+2.0*mHp2)*(cosbma*(cosb*cosb-sinb*sinb)+2.0*sinbma*sinb*cosb))
                       /(4.0*cosb*sinb);
    double lambdaHhHpHm=((mHh2-2.0*mHp2)*(-sinbma*(cosb*cosb-sinb*sinb)-2.0*cosbma*sinb*cosb)
                         +(-4.0*m12_2/(sinb*cosb)+3.0*mHh2+2.0*mHp2)*(-sinbma*(cosb*cosb-sinb*sinb)+2.0*cosbma*sinb*cosb))
                       /(4.0*cosb*sinb);

    double xidA=0.0;
    double xilh=0.0;
    double xilHh=0.0;
    double xilA=0.0;

    if( modelflag == "type1" ) {
        xidA=1.0/tanb;
        xilh=cosa/sinb;
        xilHh=sina/sinb;
        xilA=1.0/tanb;
    }
    else if( modelflag == "type2" ) {
        xidA=-tanb;
        xilh=-sina/cosb;
        xilHh=cosa/cosb;
        xilA=-tanb;
    }
    else if( modelflag == "typeX" ) {
        xidA=1.0/tanb;
        xilh=-sina/cosb;
        xilHh=cosa/cosb;
        xilA=-tanb;
    }
    else if( modelflag == "typeY" ) {
        xidA=-tanb;
        xilh=cosa/sinb;
        xilHh=sina/sinb;
        xilA=1.0/tanb;
    }
    else {
        throw std::runtime_error("modelflag can be only any of \"type1\", \"type2\", \"typeX\" or \"typeY\"");
    }

    double f1=(xt-xHp+xHp*log(xHp)-xt*log(xt))/(2.0*(xHp-xt));
    double f2=(xt-(xHp*xt*(log(xHp)-log(xt)))/(xHp-xt))/(2.0*(xHp-xt));
    double f3=(xHp-xHp*xHp*log(xHp)/(xHp-xt)+xt*(2.0*xHp-xt)*log(xt)/(xHp-xt))/(2.0*(xHp-xt));
    double f4=(0.5*xt*(3.0*xHp-xt)-(xHp*xHp*xt*(log(xHp)-log(xt)))/(xHp-xt))/(4.0*(xHp-xt)*(xHp-xt));
    double f5=(0.5*xt*(xHp-3.0*xt)-(xHp*xt*(xHp-2.0*xt)*(log(xHp)-log(xt)))/(xHp-xt))/(4.0*(xHp-xt)*(xHp-xt));
    double f6=(xt*(xt*xt-3.0*xHp*xt+9.0*xHp-5.0*xt-2.0))/(8.0*(xt-1.0)*(xt-1.0)*(xHp-xt))
              +(xHp*(xHp*xt-3.0*xHp+2.0*xt)*log(xHp))/(4.0*(xHp-1.0)*(xHp-xt)*(xHp-xt))
              +((xHp*xHp*(-2.0*xt*xt*xt+6.0*xt*xt-9.0*xt+2.0)+3.0*xHp*xt*xt*(xt*xt-2.0*xt+3.0)-xt*xt*(2.0*xt*xt*xt-3.0*xt*xt+3.0*xt+1.0))*log(xt))
               /(4.0*(xt-1.0)*(xt-1.0)*(xt-1.0)*(xHp-xt)*(xHp-xt));
    double f7=(xt*xt+xt-8.0)/(8.0*(xt-1.0)*(xt-1.0))-(xHp*(xHp+2.0)*log(xHp))/(4.0*(xHp-1.0)*(xHp-xt))
              +((xHp*(xt*xt*xt-3.0*xt*xt+3.0*xt+2.0)+3.0*xt*xt*(xt-2.0))*log(xt))/(4.0*(xt-1.0)*(xt-1.0)*(xt-1.0)*(xHp-xt));

    double g0=-(-xiuA*xidA*(f1+f2+f3+1.0)+xiuA*xiuA*(f4-f5-0.25))/(4.0*xHp);
    double g1a=-0.75-xiuA*xidA*(f1+f2+f3)+xiuA*xiuA*(f4-f5);
    double g2a=-xidA*xidA*xiuA*f1+xidA*xiuA*xiuA*(f3+f2)+xiuA*xiuA*xiuA*(f5-f4)+xiuA*(f7-f6)+xidA*f1;

    double CSboxSM=-xt*(xt+1.0)/(48.0*(xt-1.0)*(xt-1.0))-(xt-2.0)*(3.0*xt*xt-3.0*xt+1.0)*log(xt)/(24.0*(xt-1.0)*(xt-1.0)*(xt-1.0));
    double CSboxTHDM=-xiuA*xilA*xt*(1.0-xHp*log(xHp/xt)/(xHp-xt))/(8.0*(xHp-xt)) -xidA*xilA*xt*log(xHp/xt)/(4.0*(xHp-xt));
    double CSphiTHDM=xt*xilh*(sinbma*g1a+cosbma*g2a+2.0/(MW*MW)*lambdahHpHm*g0)/(2.0*xh)
                     +xt*xilHh*(cosbma*g1a-sinbma*g2a+2.0/(MW*MW)*lambdaHhHpHm*g0)/(2.0*xHh);

    return CSboxSM+CSboxTHDM+CSphiTHDM;
}

BR_BsmumuTHDM::BR_BsmumuTHDM(const StandardModel& SM_i)
: B_to_mumu(SM_i)
{}

double BR_BsmumuTHDM::computeThValue()
{

    //Expressions taken from 1511.01829v2 and the SM part

    double Mb=myTHDM->getQuarks(QCD::BOTTOM).getMass();   
    double Ms=myTHDM->getQuarks(QCD::STRANGE).getMass();
    double GF=myTHDM->getGF();
    double mMu=myTHDM->getLeptons(StandardModel::MU).getMass();
    double mBs=myTHDM->getMesons(QCD::B_S).getMass();
    double MW=myTHDM->Mw_tree();
    double GammaBs=myTHDM->getMesons(QCD::B_S).computeWidth();
    double FBs=myTHDM->getMesons(QCD::B_S).getDecayconst();
    double VtbstarVtsAbs2=1.0*1.0*0.04107*0.04107; /*replace this by the CKM matrix elements*/

    double CS=computeCS();      //this is new (in the SM part the box part is neglected)
    double C10=computeC10();    //this contains the SM contribution and a THDM part
    double CP=computeCP();      //this is also new and contains the pseudoscalar contributions
    double beta = sqrt(1.0 - 4. * mMu*mMu / (mBs*mBs));
    double ampSq = VtbstarVtsAbs2 * ( pow(C10+(CP*Mb*mBs)/(2.0*(Mb+Ms)*MW*MW),2) + pow(CS*Mb*mBs*mBs*beta/(2.0*(Mb+Ms)*MW*MW),2) );

    //from here consistent with the SM notation
    double coupling = GF*GF * MW*MW /(M_PI*M_PI);
    double PRF = coupling*coupling * FBs*FBs * mMu*mMu * mBs * beta / (8.0*M_PI*GammaBs);

    return PRF * ampSq;
}

BR_BdmumuTHDM::BR_BdmumuTHDM(const StandardModel& SM_i)
: B_to_mumu(SM_i)
{}

double BR_BdmumuTHDM::computeThValue()
{

    double Mb=myTHDM->getQuarks(QCD::BOTTOM).getMass();   
    double Md=myTHDM->getQuarks(QCD::DOWN).getMass();
    double GF=myTHDM->getGF();
    double mMu=myTHDM->getLeptons(StandardModel::MU).getMass();
    double mBd=myTHDM->getMesons(QCD::B_D).getMass();
    double MW=myTHDM->Mw_tree();
    double GammaBd=myTHDM->getMesons(QCD::B_D).computeWidth();
    double FBd=myTHDM->getMesons(QCD::B_D).getDecayconst();
    double VtbstarVtdAbs2=1.0*1.0*0.008676*0.008676; /*replace this by the CKM matrix elements*/

    double CS=computeCS();
    double C10=computeC10();
    double CP=computeCP();
    double beta = sqrt(1.0 - 4. * mMu*mMu / (mBd*mBd));
    double ampSq = VtbstarVtdAbs2 * ( pow(C10+(CP*Mb*mBd)/(2.0*(Mb+Md)*MW*MW),2) + pow(CS*Mb*mBd*mBd*beta/(2.0*(Mb+Md)*MW*MW),2) );

    //from here consistent with the SM notation
    double coupling = GF*GF * MW*MW /(M_PI*M_PI);
    double PRF = coupling*coupling * FBd*FBd * mMu*mMu * mBd * beta / (8.0*M_PI*GammaBd);

    return PRF * ampSq;
}
