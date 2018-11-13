/*
 * Copyright (C) 2018 HEPfit Collaboration
 *
 *
 * For the licensing terms see doc/COPYING.
 */


#include "THDMWMatching.h"
#include "THDMW.h"
#include <gsl/gsl_integration.h>
#include <gsl/gsl_sf_dilog.h>
#include "THDMWcache.h"
#include <stdexcept>

THDMWMatching::THDMWMatching(const THDMW & THDMW_i) :

    StandardModelMatching(THDMW_i),
    myTHDMW(THDMW_i),
    myCKM(3, 3, 0.),
    mcBMll(13, NDR, NLO),
    mcbsg(8, NDR, NNLO),
    mcdbs2(5, NDR, NLO),
    mcdbsp2(5, NDR, NLO)
{
}





std::vector<WilsonCoefficient>& THDMWMatching::CMBMll(QCD::lepton lepton)
{           
    //From 1504.00839
    double Mut = myTHDMW.getMut();
    double Muw = myTHDMW.getMuw();
    double MW = myTHDMW.Mw();
    double xt = x_t(Mut);
    double mhsq = myTHDMW.getMyTHDMWCache()->mhsq;
    double xh = mhsq/pow(MW,2);
    double sW2 = myTHDMW.sW2();
    double mSp2 =myTHDMW.getMyTHDMWCache()->mSpsq;
    double xS = mSp2/pow(MW,2);
    double nu1 = myTHDMW.getTHDMW_nu1();
    gslpp::complex etaU = myTHDMW.getTHDMW_etaU();
    gslpp::complex etaD = myTHDMW.getTHDMW_etaD();
    
    double mb=myTHDMW.getQuarks(QCD::BOTTOM).getMass();
   
    double ml=myTHDMW.getLeptons(lepton).getMass();
    
    
    
    vmcBMll = StandardModelMatching::CMBMll(lepton);
      switch (mcbsg.getScheme()) {
        case NDR:

            break;
        default:
        std::stringstream out;
       out << mcBMll.getScheme();
     throw std::runtime_error("THDMWMatching::CMBMll(): scheme " + out.str() + "not implemented");
    }
    mcBMll.setMu(Muw);
    
     switch (mcBMll.getOrder()) {
        case NNLO:
        case NLO:           
        case LO:           
            mcBMll.setCoeff(9 , C10NP(xt, xS, etaU)/(sW2), LO);
            mcBMll.setCoeff(10 , CSNP(nu1, xh, xt, xS, etaU, etaD)*mb*ml/(MW*MW*sW2), LO);
            mcBMll.setCoeff(11 , CPNP(xt, xS, etaU, etaD)*(mb*ml)/(MW*MW*sW2), LO);
            break;
        default:
            std::stringstream out;
            out << mcBMll.getOrder();
            throw std::runtime_error("THDMWMatching::CMBMll(): order " + out.str() + "not implemented"); 
            }
    vmcBMll.push_back(mcBMll);
    return (vmcBMll);
}

























gslpp::complex THDMWMatching::C10NP(double xt, double xS , gslpp::complex etaU) {
    
    gslpp::complex C10 = etaU.abs2()*(pow(xt,2)/6)*(1/(xS-xt)+(xS/pow(xS-xt,2))*(log(xt)-log(xS)));
    return C10;
    }

gslpp::complex THDMWMatching::CSNP(double nu1, double xh,double xt, double xS, gslpp::complex etaU, gslpp::complex etaD) {
    double mhsq= myTHDMW.getMyTHDMWCache()->mhsq;
    double vev= myTHDMW.v();
    
    gslpp::complex CS =(pow(vev,2)*xt*nu1/(mhsq*6*(xS-xt)))*(etaD*etaU.conjugate()*(xt*(log(xS)-log(xt))/(xS-xt)-1)+
            etaU.abs2()*(pow(xt,2)*(log(xS)-log(xt))/(2*pow(xS-xt,2))+(xS-3*xt)/(4*(xS-xt))))+
            (2*xt/(3*xh))*(etaD*etaU.conjugate()*xt*(1-xS*(log(xS)-log(xt))/(xS-xt))/(xS-xt)+
            etaU.abs2()*(xt/(2*pow(xS-xt,2)))*((xS+xt)/2-xS*xt*(log(xS)-log(xt))/(xS-xt)));
    return CS;
    }

gslpp::complex THDMWMatching::CPNP(double xt, double xS, gslpp::complex etaU, gslpp::complex etaD) {
    double sW2 = myTHDMW.sW2();
    gslpp::complex CP =(xt/(3*pow(xS-xt,2)))*(etaD*etaU.conjugate()*(-(xt+xS)/2+xt*xS*(log(xS)-log(xt))/(xS-xt))+
    (etaU.abs2()/(6*(xS-xt)))*((pow(xS,2)-8*xS*xt-17*pow(xt,2))/6+pow(xt,2)*(3*xS+xt)*(log(xS)-log(xt))/(xS-xt)))+
    (2*sW2*xt/(9*pow(xS-xt,2)))*(etaD*etaU.conjugate()*((5*xt-3*xS)/2+xS*(2*xS-3*xt)*(log(xS)-log(xt))/(xS-xt))
    +(etaU.abs2()/(6*(xS-xt)))*((17*pow(xS,2)-64*xS*xt+71*pow(xt,2))/6)
    -(4*pow(xS,3)-12*pow(xS,2)*xt+9*xS*pow(xt,2)+3*pow(xt,3))*(log(xS)-log(xt))/(xS-xt))
    +etaU.abs2()*(1-sW2)*pow(xt,2)*(xS*(log(xS)-log(xt))+xt-xS)/(3*pow(xS-xt,2));
    return CP;
    }












//For including the contribution to DelB_s we'll use the result of 1504.00839. 
//The prediction of SM in this paper is the same than HEPfit except for lambdat_s which is conjugated (There are no additional factors)
//Note that in 1504.00839 what is computed is the decay of a B\bar_s^0 to a B_s^0 but usually the matrix elements are define for the conjugated elements.
//We should take the complex conjugate in order to be consistent with the HEPfit calculations.


std::vector<WilsonCoefficient>& THDMWMatching::CMdbs2() {

    double Mut = myTHDMW.getMut();
    double Muc = myTHDMW.getMuc();
    double Mub = myTHDMW.getMub();
    double xt = x_t(Mut);//ratio mass of the top square at the scale of its mass over mass of the W square. Why not Mw scale???
    double xc= x_c(Muc);
    double mb =myTHDMW.Mrun(Mub, myTHDMW.getQuarks(QCD::BOTTOM).getMass_scale(),
                              myTHDMW.getQuarks(QCD::BOTTOM).getMass(), FULLNNLO);
    double mSp2=myTHDMW.getMyTHDMWCache()->mSpsq;
    double MW=myTHDMW.Mw();
    double xb=pow(mb,2)/(MW*MW);
    
    double xS=mSp2/(MW*MW);
    gslpp::complex etaU = myTHDMW.getTHDMW_etaU();
    gslpp::complex etaD = myTHDMW.getTHDMW_etaD();
    double GF=myTHDMW.getGF();
    double co = GF / 4. / M_PI * MW;
    


    vmcds = StandardModelMatching::CMdbs2();
    mcdbs2.setMu(Mut);

    switch (mcdbs2.getOrder()) {
        case NNLO:
        case NLO:
        case LO:
            mcdbs2.setCoeff(0, co * co * CNPVLL(xc, xb, xt, xS, etaU, etaD).conjugate(), LO);
            break;
        default:
            std::stringstream out;
            out << mcdbs2.getOrder();
            throw std::runtime_error("THDMMatching::CMdbs2(): order " + out.str() + "not implemented");
    }

    vmcds.push_back(mcdbs2);
    return(vmcds);
}

std::vector<WilsonCoefficient>& THDMWMatching::CMdbsp2() {

    double Mut = myTHDMW.getMut();
    double Muc = myTHDMW.getMuc();
    double Mub = myTHDMW.getMub();
    double xt = x_t(Mut);//ratio mass of the top square at the scale of its mass over mass of the W square. Why not Mw scale???
    double xc= x_c(Muc);
    double mb =myTHDMW.Mrun(Mub, myTHDMW.getQuarks(QCD::BOTTOM).getMass_scale(),
                              myTHDMW.getQuarks(QCD::BOTTOM).getMass(), FULLNNLO);
    double mSp2=myTHDMW.getMyTHDMWCache()->mSpsq;
    double MW=myTHDMW.Mw();
    double xb=pow(mb,2)/(MW*MW);
    
    double xS=mSp2/(MW*MW);
    gslpp::complex etaU = myTHDMW.getTHDMW_etaU();
    gslpp::complex etaD = myTHDMW.getTHDMW_etaD();
    double GF=myTHDMW.getGF();
    double co = GF / 4. / M_PI * MW;
    
    vmcdsp = StandardModelMatching::CMdbs2();
    mcdbsp2.setMu(Mut);

    switch (mcdbsp2.getOrder()) {
        case NNLO:
        case NLO:
        case LO:
            mcdbsp2.setCoeff(1, co * co * CNPSRR1(xc, xb, xt, xS, etaU, etaD).conjugate()+4*co * co * CNPSRR2(xc, xb, xt, xS, etaU, etaD).conjugate(), LO);
            mcdbsp2.setCoeff(2, 8*co * co * CNPSRR2(xc, xb, xt, xS, etaU, etaD).conjugate(), LO);
            //OSRR2=4*O2+8*O3 -> From Buras basis to SUSY basis
            break;
        default:
            std::stringstream out;
            out << mcdbsp2.getOrder();
            throw std::runtime_error("THDMMatching::CMdbs2(): order " + out.str() + "not implemented");
    }

    vmcdsp.push_back(mcdbsp2);
    return(vmcdsp);
}







double THDMWMatching::CNPVLLctEtadEtasu(double xc, double xb, double xt, double xS){
      return (-2*xb*xc*xt*(-1 - pow(xS,2) + xt + xS*xt))/(3.*(-1 + xS)*(-xS + xc*(1 + xS))*(xS - xt)*(-1 + xt)) + 
   (2*xb*xc*xt*log(xc))/(3.*(-(xS*xt) + xc*(2*xt + xS*(1 + 2*xt)))) + 
   (2*xb*xc*xt*(-2*pow(xS,2) + xt + xS*xt)*log(xS))/(3.*pow(-1 + xS,2)*(-2*xc + xS)*pow(xS - xt,2)) + 
   (2*xb*xc*xt*(xS + (-2 + xt)*xt)*log(xt))/(3.*pow(xS - xt,2)*pow(-1 + xt,2)*(-xc + xt));
  }

double THDMWMatching::CNPVLLctEtau4(double xc, double xb, double xt, double xS){
      return (11*xc*pow(xt,2)*(-18*xS*pow(xS - xt,2) + xb*(5*pow(xS,2) - 22*xS*xt + 5*pow(xt,2))))/(324.*pow(xS - xt,3)*(-(xS*xt) + xc*(2*xS + 3*xt))) - 
   (11*xc*xS*xt*(xb*xS*(xS - 3*xt) - 3*pow(xS - xt,2)*xt)*log(xS))/(54.*(4*xc - xS)*pow(xS - xt,4)) + 
   (11*xc*pow(xt,2)*(xb*xS*(xS - 3*xt) - 3*pow(xS - xt,2)*xt)*log(xt))/(54.*(3*xc - xt)*pow(xS - xt,4));
  }

double THDMWMatching::CNPVLLctEtau2(double xc, double xb, double xt, double xS){
      return -((f10(xt,xS) + f11(xt,xS) + f12(xt,xS))*f9(xc,xb,xt,xS)) + 
   f1(xc,xt,xS)*(f2(xb,xt,xS) + f3(xb,xt,xS) + f4(xb,xt,xS))*log(xS) + 
   f5(xc,xt,xS)*(f6(xb,xt,xS) + f7(xb,xt,xS) + f8(xb,xt))*log(xt);
  }

gslpp::complex THDMWMatching::CNPVLLct(double xc, double xb, double xt, double xS,  gslpp::complex etaU, gslpp::complex etaD){
      return etaD*etaU.conjugate()*CNPVLLctEtadEtasu(xc ,xb ,xt , xS)+pow(etaU.abs2(),2)*CNPVLLctEtau4(xc ,xb ,xt , xS) +
              etaU.abs2()*CNPVLLctEtau2(xc ,xb ,xt , xS);
  }

double THDMWMatching::CNPVLLttEtadEtasu(double xc, double xb, double xt, double xS){
      return (2*xb*xt*(pow(xS,2)*(1 + xt) - xS*pow(1 + xt,2) + xt*(3 - 3*xt + 2*pow(xt,2))))/(3.*(-1 + xS)*pow(xS - xt,2)*pow(-1 + xt,2)) + 
   (4*xb*pow(xt,2)*(-pow(xS,2) + xt)*log(xS))/(3.*pow(-1 + xS,2)*pow(xS - xt,3)) + 
   (4*xb*pow(xt,2)*(-pow(xS,2) + xt + 2*xS*xt - 3*pow(xt,2) + pow(xt,3))*log(xt))/(3.*pow(xS - xt,3)*pow(-1 + xt,3));
  }

double THDMWMatching::CNPVLLttEtau4(double xc, double xb, double xt, double xS){
      return (-11*pow(xt,2)*(-9*pow(xS - xt,2)*(xS + xt) + xb*(5*pow(xS,2) - 22*xS*xt + 5*pow(xt,2))))/(162.*pow(xS - xt,4)) + (11*pow(xt,2)*(-6*xS*pow(xS - xt,2)*xt + 
        xb*(pow(xS,3) - 3*pow(xS,2)*xt - 3*xS*pow(xt,2) + pow(xt,3)))*log(xS))/(54.*pow(xS - xt,5)) - 
   (11*pow(xt,2)*(-6*xS*pow(xS - xt,2)*xt + 
        xb*(pow(xS,3) - 3*pow(xS,2)*xt - 3*xS*pow(xt,2) + pow(xt,3)))*log(xt))/(54.*pow(xS - xt,5));
  }

double THDMWMatching::CNPVLLttEtau2(double xc, double xb, double xt, double xS){
      return  f24(xt,xS)*(f30(xb,xt,xS) + f31(xb,xt,xS)) + 
   f13(xt,xS)*(f14(xb,xt,xS) + f15(xb,xt,xS) + f16(xb,xt,xS))*log(xS) - 
   f17(xt,xS)*(f22(xb,xt,xS) + f23(xb,xt,xS))*log(xt);
  }

gslpp::complex THDMWMatching::CNPVLLtt(double xc, double xb, double xt, double xS,  gslpp::complex etaU, gslpp::complex etaD){
      return etaD*etaU.conjugate()*CNPVLLttEtadEtasu(xc ,xb ,xt , xS)+pow(etaU.abs2(),2)*CNPVLLttEtau4(xc ,xb ,xt , xS) +
              etaU.abs2()*CNPVLLttEtau2(xc ,xb ,xt , xS);
  }

double THDMWMatching::CNPVLLccEtadEtasu(double xc, double xb, double xS){
      return (-2*xb*xc)/(3.*(-xS + 2*xc*(1 + xS)));
  }

double THDMWMatching::CNPVLLccEtau2(double xc, double xb, double xS){
      return (-4*xb*xc)/(9.*(-xS + 3*xc*(1 + xS)));
  }

gslpp::complex THDMWMatching::CNPVLLcc(double xc, double xb, double xS,  gslpp::complex etaU, gslpp::complex etaD){
      return etaD*etaU.conjugate()*CNPVLLccEtadEtasu(xc, xb, xS) +
              etaU.abs2()*CNPVLLccEtau2(xc, xb ,xS);
  }

gslpp::complex THDMWMatching::CNPVLL(double xc, double xb, double xt, double xS,  gslpp::complex etaU, gslpp::complex etaD){
      return pow(myTHDMW.getCKM().computelamt_s().conjugate(),2)*CNPVLLtt(xc, xb, xt, xS, etaU, etaD)+
            2*myTHDMW.getCKM().computelamt_s().conjugate()*myTHDMW.getCKM().computelamc_s().conjugate()*CNPVLLct(xc, xb, xt, xS, etaU, etaD)+
            pow(myTHDMW.getCKM().computelamc_s().conjugate(),2)*CNPVLLcc(xc, xb, xS, etaU, etaD)  ;
  }











double THDMWMatching::CNPSRR1ctEtadEtasuEtau2(double xc, double xb, double xt, double xS){
      return (-19*xb*xc*xt)/(18.*(xc - xS)*(xS - xt)) + (19*xb*xc*xS*xt*log(xS))/(18.*(2*xc - xS)*pow(xS - xt,2)) - 
   (19*xb*xc*pow(xt,2)*log(xt))/(18.*(xc - xt)*pow(xS - xt,2));
  }

double THDMWMatching::CNPSRR1ctEtad2Etasu2(double xc, double xb, double xt, double xS){
      return (19*xb*xc*xt)/(18.*(xc - xS)*(xS - xt)) - (19*xb*xc*xS*xt*log(xS))/(18.*(2*xc - xS)*pow(xS - xt,2)) + 
   (19*xb*xc*pow(xt,2)*log(xt))/(18.*(xc - xt)*pow(xS - xt,2));
  }

double THDMWMatching::CNPSRR1ctEtadEtasu(double xc, double xb, double xt, double xS){
      return (-5*xb*xc*xS*(1 + xS - 2*xt)*xt)/(6.*(-1 + xS)*(-xS + xc*(1 + xS))*(xS - xt)*(-1 + xt)) + 
   (5*xb*xc*xS*xt*(-3*xt + xS*(2 + xt))*log(xS))/(6.*pow(-1 + xS,2)*(-2*xc + xS)*pow(xS - xt,2)) + 
   (5*xb*xc*pow(xt,2)*(xS*(-2 + xt) + (3 - 2*xt)*xt)*log(xt))/(6.*pow(xS - xt,2)*pow(-1 + xt,2)*(-xc + xt));
  }

double THDMWMatching::CNPSRR1ctEtau4(double xc, double xb, double xt, double xS){
      return (19*xb*xc*pow(xt,2)*(5*pow(xS,2) - 22*xS*xt + 5*pow(xt,2)))/(648.*pow(xS - xt,3)*(-(xS*xt) + xc*(2*xS + 3*xt))) - 
   (19*xb*xc*pow(xS,2)*(xS - 3*xt)*xt*log(xS))/(108.*(4*xc - xS)*pow(xS - xt,4)) + 
   (19*xb*xc*xS*(xS - 3*xt)*pow(xt,2)*log(xt))/(108.*(3*xc - xt)*pow(xS - xt,4));
  }

double THDMWMatching::CNPSRR1ctEtau2(double xc, double xb, double xt, double xS){
      return f35(xc,xb,xt,xS)*(f36(xt,xS) + f37(xt,xS)) + f32(xc,xb,xt,xS)*log(xS) - 
   f33(xc,xb,xt,xS)*f34(xt,xS)*log(xt);
  }

gslpp::complex THDMWMatching::CNPSRR1ct(double xc, double xb, double xt, double xS,  gslpp::complex etaU, gslpp::complex etaD){
      return etaD*etaU.conjugate()*etaU.abs2()*CNPSRR1ctEtadEtasuEtau2(xc ,xb ,xt , xS)+pow(etaD*etaU.conjugate(),2)*CNPSRR1ctEtad2Etasu2(xc ,xb ,xt , xS)+
              etaD*etaU.conjugate()*CNPSRR1ctEtadEtasu(xc ,xb ,xt , xS)+pow(etaU.abs2(),2)*CNPSRR1ctEtau4(xc ,xb ,xt , xS) +
              etaU.abs2()*CNPSRR1ctEtau2(xc ,xb ,xt , xS);
  }

double THDMWMatching::CNPSRR1ttEtadEtasuEtau2(double xc, double xb, double xt, double xS){
      return (19*xb*pow(xt,2))/(9.*pow(xS - xt,2)) - (19*xb*pow(xt,2)*(xS + xt)*log(xS))/(18.*pow(xS - xt,3)) + 
   (19*xb*pow(xt,2)*(xS + xt)*log(xt))/(18.*pow(xS - xt,3));
  }

double THDMWMatching::CNPSRR1ttEtad2Etasu2(double xc, double xb, double xt, double xS){
      return (-19*xb*pow(xt,2))/(9.*pow(xS - xt,2)) + (19*xb*pow(xt,2)*(xS + xt)*log(xS))/(18.*pow(xS - xt,3)) - 
   (19*xb*pow(xt,2)*(xS + xt)*log(xt))/(18.*pow(xS - xt,3));
  }

double THDMWMatching::CNPSRR1ttEtadEtasu(double xc, double xb, double xt, double xS){
      return (5*xb*pow(xt,2)*(pow(xS,2)*(-3 + xt) + (-3 + xt)*xt + xS*(1 + 6*xt - 3*pow(xt,2))))/(6.*(-1 + xS)*pow(xS - xt,2)*pow(-1 + xt,2)) + 
   (5*xb*xS*pow(xt,2)*(-2*xt + xS*(1 + xt))*log(xS))/(3.*pow(-1 + xS,2)*pow(xS - xt,3)) + 
   (5*xb*pow(xt,2)*(pow(xS,2) - 2*xS*xt - (-2 + xt)*pow(xt,3))*log(xt))/(3.*pow(xS - xt,3)*pow(-1 + xt,3));
  }

double THDMWMatching::CNPSRR1ttEtau4(double xc, double xb, double xt, double xS){
      return (-19*xb*pow(xt,2)*(5*pow(xS,2) - 22*xS*xt + 5*pow(xt,2)))/(324.*pow(xS - xt,4)) + 
   (19*xb*pow(xt,2)*(pow(xS,3) - 3*pow(xS,2)*xt - 3*xS*pow(xt,2) + pow(xt,3))*log(xS))/(108.*pow(xS - xt,5)) 
    - (19*xb*pow(xt,2)*(pow(xS,3) - 3*pow(xS,2)*xt - 3*xS*pow(xt,2) + pow(xt,3))*log(xt))/(108.*pow(xS - xt,5));
  }

double THDMWMatching::CNPSRR1ttEtau2(double xc, double xb, double xt, double xS){
      return -(f40(xb,xt,xS)*(f41(xt,xS) + f42(xt,xS) + f43(xt,xS))) + f38(xb,xt,xS)*log(xS) + 
   f39(xb,xt,xS)*log(xt);
  }

gslpp::complex THDMWMatching::CNPSRR1tt(double xc, double xb, double xt, double xS,  gslpp::complex etaU, gslpp::complex etaD){
      return etaD*etaU.conjugate()*etaU.abs2()*CNPSRR1ttEtadEtasuEtau2(xc ,xb ,xt , xS)+pow(etaD*etaU.conjugate(),2)*CNPSRR1ttEtad2Etasu2(xc ,xb ,xt , xS)+
              etaD*etaU.conjugate()*CNPSRR1ttEtadEtasu(xc ,xb ,xt , xS)+pow(etaU.abs2(),2)*CNPSRR1ttEtau4(xc ,xb ,xt , xS) +
              etaU.abs2()*CNPSRR1ttEtau2(xc ,xb ,xt , xS);
  }

gslpp::complex THDMWMatching::CNPSRR1(double xc, double xb, double xt, double xS,  gslpp::complex etaU, gslpp::complex etaD){
      return pow(myTHDMW.getCKM().computelamt_s().conjugate(),2)*CNPSRR1tt(xc, xb, xt, xS, etaU, etaD)+
            2*myTHDMW.getCKM().computelamt_s().conjugate()*myTHDMW.getCKM().computelamc_s().conjugate()*CNPSRR1ct(xc, xb, xt, xS, etaU, etaD);
  }











double THDMWMatching::CNPSRR2ctEtadEtasuEtau2(double xc, double xb, double xt, double xS){
      return -21*CNPSRR1ctEtadEtasuEtau2(xc, xb, xt, xS)/76;
  }

double THDMWMatching::CNPSRR2ctEtad2Etasu2(double xc, double xb, double xt, double xS){
      return -21*CNPSRR1ctEtad2Etasu2(xc, xb, xt, xS)/76;
  }

double THDMWMatching::CNPSRR2ctEtadEtasu(double xc, double xb, double xt, double xS){
      return -3*CNPSRR1ctEtadEtasu(xc, xb, xt, xS)/20;
  }

double THDMWMatching::CNPSRR2ctEtau4(double xc, double xb, double xt, double xS){
      return -21*CNPSRR1ctEtau4(xc, xb, xt, xS)/76;
  }

double THDMWMatching::CNPSRR2ctEtau2(double xc, double xb, double xt, double xS){
      return -3*CNPSRR1ctEtau2(xc, xb, xt, xS)/20;
  }

gslpp::complex THDMWMatching::CNPSRR2ct(double xc, double xb, double xt, double xS,  gslpp::complex etaU, gslpp::complex etaD){
      return etaD*etaU.conjugate()*etaU.abs2()*CNPSRR2ctEtadEtasuEtau2(xc ,xb ,xt , xS)+pow(etaD*etaU.conjugate(),2)*CNPSRR2ctEtad2Etasu2(xc ,xb ,xt , xS)+
              etaD*etaU.conjugate()*CNPSRR2ctEtadEtasu(xc ,xb ,xt , xS)+pow(etaU.abs2(),2)*CNPSRR2ctEtau4(xc ,xb ,xt , xS) +
              etaU.abs2()*CNPSRR2ctEtau2(xc ,xb ,xt , xS);
  }




double THDMWMatching::CNPSRR2ttEtadEtasuEtau2(double xc, double xb, double xt, double xS){
      return -21*CNPSRR1ttEtadEtasuEtau2(xc, xb, xt, xS)/76;
  }

double THDMWMatching::CNPSRR2ttEtad2Etasu2(double xc, double xb, double xt, double xS){
      return -21*CNPSRR1ttEtad2Etasu2(xc, xb, xt, xS)/76;
  }

double THDMWMatching::CNPSRR2ttEtadEtasu(double xc, double xb, double xt, double xS){
      return -3*CNPSRR1ttEtadEtasu(xc, xb, xt, xS)/20;
  }

double THDMWMatching::CNPSRR2ttEtau4(double xc, double xb, double xt, double xS){
      return -21*CNPSRR1ttEtau4(xc, xb, xt, xS)/76;
  }

double THDMWMatching::CNPSRR2ttEtau2(double xc, double xb, double xt, double xS){
      return -3*CNPSRR1ttEtau2(xc, xb, xt, xS)/20;
  }

gslpp::complex THDMWMatching::CNPSRR2tt(double xc, double xb, double xt, double xS,  gslpp::complex etaU, gslpp::complex etaD){
      return etaD*etaU.conjugate()*etaU.abs2()*CNPSRR2ttEtadEtasuEtau2(xc ,xb ,xt , xS)+pow(etaD*etaU.conjugate(),2)*CNPSRR2ttEtad2Etasu2(xc ,xb ,xt , xS)+
              etaD*etaU.conjugate()*CNPSRR2ttEtadEtasu(xc ,xb ,xt , xS)+pow(etaU.abs2(),2)*CNPSRR2ttEtau4(xc ,xb ,xt , xS) +
              etaU.abs2()*CNPSRR2ttEtau2(xc ,xb ,xt , xS);
  }



gslpp::complex THDMWMatching::CNPSRR2(double xc, double xb, double xt, double xS,  gslpp::complex etaU, gslpp::complex etaD){
      return pow(myTHDMW.getCKM().computelamt_s().conjugate(),2)*CNPSRR2tt(xc, xb, xt, xS, etaU, etaD)+
            2*myTHDMW.getCKM().computelamt_s().conjugate()*myTHDMW.getCKM().computelamc_s().conjugate()*CNPSRR2ct(xc, xb, xt, xS, etaU, etaD);
  }










//Functions Needed from 1504.00839

double THDMWMatching::f1(double xc, double xt, double xS){
      double f1 = -(xc*xS*xt)/(9.*pow(-1 + xS,3)*(-3*xc + xS)*pow(xS - xt,3));
      return f1;
  }
        
double THDMWMatching::f2(double xb, double xt, double xS){
      double f2 = -6*pow(xS,5) + 12*pow(xS,4)*(3 + xt) - 3*pow(xS,3)*(xb*(-6 + xt) 
      + 2*(9 + 12*xt + pow(xt,2)));
      return f2;
  }

double THDMWMatching::f3(double xb, double xt, double xS){
      double f3 = pow(xS,2)*(xb*(2 - 8*xt + pow(xt,2)) + 12*(2 + 9*xt + 3*pow(xt,2)));
      return f3;
  }

double THDMWMatching::f4(double xb, double xt, double xS){
      double f4 = xS*xt*(xb*(-41 + 6*xt) - 6*(8 + 9*xt)) + xt*(24*xt + xb*(12 + 13*xt));
      return f4;
  }

double THDMWMatching::f5(double xc, double xt, double xS){
      double f5 = (xc*pow(xt,2))/(9.*pow(xS - xt,3)*pow(-1 + xt,3)*(-3*xc + xt));
      return f5;
  }

double THDMWMatching::f6(double xb, double xt, double xS){
      double f6 = pow(xS,2)*(-6*(-4 + xt)*pow(-1 + xt,2) + xb*(2 + 7*xt + pow(xt,2)));
      return f6;
  }

double THDMWMatching::f7(double xb, double xt, double xS){
      double f7 = xS*xt*(12*(-4 + xt)*pow(-1 + xt,2) - xb*(5 + 12*xt + 3*pow(xt,2)));
      return f7;
  }

double THDMWMatching::f8(double xb, double xt){
      double f8 = xt*(-6*(-4 + xt)*pow(-1 + xt,2)*xt + xb*(12 - 23*xt + 21*pow(xt,2)));
      return f8;
  }

double THDMWMatching::f9(double xc,double xb, double xt, double xS){
      double f9 = (2*xb*xc*xt)/(9.*pow(-1 + xS,2)*pow(xS - xt,2)*pow(-1 + xt,2)*(2*xc*xS + 2*xc*xt - xS*xt + 2*xc*xS*xt));
      return f9;
  }

double THDMWMatching::f10(double xt, double xS){
      double f10 = pow(xS,4)*(6 - 15*xt + 4*pow(xt,2)) + pow(xS,3)*(-12 + 19*xt + 17*pow(xt,2) - 9*pow(xt,3));
      return f10;
  }

double THDMWMatching::f11(double xt, double xS){
      double f11 = pow(xS,2)*(6 + 21*xt - 70*pow(xt,2) + 24*pow(xt,3) + 4*pow(xt,4)) ;
      return f11;
  }

double THDMWMatching::f12(double xt, double xS){
      double f12 = 3*pow(-1 + xt,2)*pow(xt,2) + xS*xt*(-15 + 16*xt + 21*pow(xt,2) - 17*pow(xt,3));
      return f12;
  }

double THDMWMatching::f13(double xt, double xS){
      double f13 = (2*xS*pow(xt,2))/(9.*pow(-1 + xS,3)*pow(xS - xt,4));
      return f13;
  }

double THDMWMatching::f14(double xb, double xt, double xS){
      double f14 = 3*pow(xS,5) - 6*pow(xS,4)*(3 + xt) + 3*pow(xS,3)*(9 + xb*(-3 + xt) + 12*xt + pow(xt,2));
      return f14;
  }

double THDMWMatching::f15(double xb, double xt, double xS){
      double f15 = -(pow(xS,2)*(xb*(1 + 10*xt + pow(xt,2)) + 6*(2 + 9*xt + 3*pow(xt,2))));
      return f15;
  }

double THDMWMatching::f16(double xb, double xt, double xS){
      double f16 = -12*xt*(xt + xb*(1 + xt)) + 3*xS*xt*(8 + 9*xt + xb*(13 + xt));
      return f16;
  }

double THDMWMatching::f17( double xt, double xS){
      double f17 = (2*pow(xt,2))/(9.*pow(xS - xt,4)*pow(-1 + xt,4));
      return f17;
  }

double THDMWMatching::f18(double xb, double xt, double xS){
      double f18 = pow(xS,3)*(xb*(1 + 9*xt) + 3*pow(-1 + xt,2)*(4 - 2*xt + pow(xt,2)));
      return f18;
  }

double THDMWMatching::f19(double xb, double xt, double xS){
      double f19 = -3*pow(xS,2)*xt*(xb*(1 + 9*xt) + pow(-1 + xt,2)*(8 - xt + 2*pow(xt,2)));
      return f19;
  }

double THDMWMatching::f20(double xb, double xt, double xS){
      double f20 = 3*xS*xt*(xt*pow(-2 + xt + pow(xt,2),2) + xb*(4 - 12*xt + 24*pow(xt,2) - 7*pow(xt,3) + pow(xt,4)));
      return f20;
  }

double THDMWMatching::f21(double xb, double xt){
      double f21 = -(pow(xt,4)*(9*pow(-1 + xt,2) + xb*(12 - 3*xt + pow(xt,2))));
      return f21;
  }

double THDMWMatching::f22(double xb, double xt, double xS){
      double f22 = f18(xb,xt,xS)+f19(xb,xt,xS);
      return f22;
  }

double THDMWMatching::f23(double xb, double xt, double xS){
      double f23 = f20(xb,xt,xS)+f21(xb,xt);
      return f23;
  }

double THDMWMatching::f24(double xt, double xS){
      double f24 = xt/(27.*pow(-1 + xS,2)*pow(xS - xt,3)*pow(-1 + xt,3));
      return f24;
  }

double THDMWMatching::f25(double xb, double xt, double xS){
      double f25 = xb*pow(xS,4)*(-12 + 65*xt + 2*pow(xt,2) + 5*pow(xt,3)) - 
   2*xb*pow(xS,3)*(-12 + 47*xt + 85*pow(xt,2) - 11*pow(xt,3) + 11*pow(xt,4));
      return f25;
  }

double THDMWMatching::f26(double xb, double xt, double xS){
      double f26 = xb*pow(xS,2)*(-12 - 91*xt + 574*pow(xt,2) - 246*pow(xt,3) + 130*pow(xt,4) + 5*pow(xt,5));
      return f26;
  }

double THDMWMatching::f27(double xb, double xt, double xS){
      double f27 = -2*xb*xS*xt*(-30 + 95*xt + 49*pow(xt,2) - 11*pow(xt,3) + 17*pow(xt,4));
      return f27;
  }

double THDMWMatching::f28(double xb, double xt){
      double f28 = xb*pow(xt,2)*(24 - 43*xt + 110*pow(xt,2) - 31*pow(xt,3));
      return f28;
  }

double THDMWMatching::f29(double xt, double xS){
      double f29 = -18*pow(-1 + xS,2)*pow(xS - xt,2)*(-4 + xt)*pow(-1 + xt,2)*xt;
      return f29;
  }

double THDMWMatching::f30(double xb, double xt, double xS){
      double f30 = f25(xb,xt,xS)+f26(xb,xt,xS)+f26(xb,xt,xS);
      return f30;
  }

double THDMWMatching::f31(double xb, double xt, double xS){
      double f31 = f28(xb,xt)+f29(xt,xS);
      return f31;
  }

double THDMWMatching::f32(double xc, double xb, double xt, double xS){
      double f32 = (5*xb*xc*xS*xt*(4*pow(xt,2) - xS*xt*(5 + 3*xt) + pow(xS,2)*(2 + xt + pow(xt,2))))/(9.*pow(-1 + xS,3)*(-3*xc + xS)*pow(xS - xt,3));
      return f32;
  }

double THDMWMatching::f33(double xc, double xb, double xt, double xS){
      double f33 =  (5*xb*xc*pow(xt,2))/(9.*pow(xS - xt,3)*pow(-1 + xt,3)*(-3*xc + xt));
      return f33;
  }

double THDMWMatching::f34(double xt, double xS){
      double f34 =xS*xt*(-5 + 6*xt - 3*pow(xt,2)) + pow(xS,2)*(2 - 2*xt + pow(xt,2)) + pow(xt,2)*(4 - 6*xt + 3*pow(xt,2));
      return f34;
  }

double THDMWMatching::f35(double xc, double xb, double xt, double xS){
      double f35 =  (5*xb*xc*xS*pow(xt,2))/(18.*pow(-1 + xS,2)*pow(xS - xt,2)*pow(-1 + xt,2)*(-(xS*xt) + 2*xc*(xt + xS*(1 + xt))));
      return f35;
  }

double THDMWMatching::f36(double xt, double xS){
      double f36 = pow(xS,3)*(-3 + xt) + pow(xS,2)*(4 + 5*xt - 3*pow(xt,2)) ;
      return f36;
  }

double THDMWMatching::f37(double xt, double xS){
      double f37 = xt*(-5 + 15*xt - 8*pow(xt,2)) + xS*(3 - 13*xt + 4*pow(xt,3)) ;
      return f37;
  }

double THDMWMatching::f38(double xb, double xt, double xS){
      double f38 = (10*xb*xS*pow(xt,2)*(3*pow(xt,2) - 3*xS*xt*(1 + xt) + pow(xS,2)*(1 + xt + pow(xt,2))))/(9.*pow(-1 + xS,3)*pow(xS - xt,4));
      return f38;
  }

double THDMWMatching::f39(double xb, double xt, double xS){
      double f39 = (10*xb*pow(xt,2)*(pow(xS,3) - 3*pow(xS,2)*xt + 3*xS*pow(xt,2) - pow(xt,4)*(3 - 3*xt + pow(xt,2))))/(9.*pow(xS - xt,4)*pow(-1 + xt,4));
      return f39;
  }

double THDMWMatching::f40(double xb, double xt, double xS){
      double f40 = (5*xb*pow(xt,2))/(27.*pow(-1 + xS,2)*pow(xS - xt,3)*pow(-1 + xt,3));
      return f40;
  }

double THDMWMatching::f41(double xt, double xS){
      double f41 = pow(xS,4)*(11 - 7*xt + 2*pow(xt,2)) - pow(xS,3)*(19 + 17*xt - 19*pow(xt,2) + 7*pow(xt,3));
      return f41;
  }

double THDMWMatching::f42(double xt, double xS){
      double f42 = pow(xS,2)*(2 + 55*xt - 15*pow(xt,2) - 17*pow(xt,3) + 11*pow(xt,4));
      return f42;
  }

double THDMWMatching::f43(double xt, double xS){
      double f43 = pow(xt,2)*(11 - 7*xt + 2*pow(xt,2)) - xS*xt*(7 + 53*xt - 55*pow(xt,2) + 19*pow(xt,3));
      return f43;
  }