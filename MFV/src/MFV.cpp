/* 
 * File:   MFV.cpp
 * Author: silvest
 * 
 * Created on September 24, 2010, 10:53 AM
 */

#include "MFV.h"

MFV::MFV(const gslpp::matrix<gslpp::complex>& VCKM_i, double mu_i, double md_i,
        double mc_i, double ms_i, double mt_i, double mb_i,
        const gslpp::matrix<gslpp::complex>& UPMNS_i, double me_i,
        double mmu_i, double mtau_i, double mnu1_i, double mnu2_i,
        double mnu3_i, double tanb_i, gslpp::complex muH_i) :
        SUSY(VCKM_i, mu_i, md_i, mc_i, ms_i, mt_i, mb_i, UPMNS_i, me_i, mmu_i,
        mtau_i, mnu1_i, mnu2_i, mnu3_i, tanb_i, muH_i)
{
}

MFV::MFV(const SUSY& SUSY_i) : SUSY(SUSY_i)
{

}

MFV::MFV(const MFV& orig) : SUSY(orig.getVCKM(), orig.getMu(),
        orig.getMd(), orig.getMc(), orig.getMs(), orig.getMt(), orig.getMb(),
        orig.getUPMNS(), orig.getMe(), orig.getMmu(), orig.getMtau(),
        orig.getMnu1(), orig.getMnu2(), orig.getMnu3(), orig.getTanb(),
        orig.getMuH())
{
}

MFV::~MFV() {
}

void MFV::setParameters(double a1, double a2, double a3, gslpp::complex a4,
        gslpp::complex a5, double b1, double b2, gslpp::complex b3, double b5,
        double b6, gslpp::complex b7, gslpp::complex b8, double c4,
        gslpp::complex c7, gslpp::complex c10, gslpp::complex c11) {

    gslpp::matrix<gslpp::complex> mQ2a(3,3,0.), mU2a(3,3,0.), mD2a(3,3,0.),
            Aua(3,3,0.), Ada(3,3,0.);
    gslpp::matrix<gslpp::complex> mL2a(3,3,0.), mn2a(3,3,0.), mcha(3,3,0.),
            Ana(3,3,0.), Acha(3,3,0.);
    gslpp::matrix<gslpp::complex> mU2(6,6,0.), mD2(6,6,0.), Yu2(3,3,0.),
            Yd2(3,3,0.);
    gslpp::matrix<double> mqd(3,3,0.), mqd2(3,3,0), mqu(3,3,0.), mqu2(3,3,0),
            I=gslpp::matrix<double>::Id(3);
    double cos2b;
   
    Yu2 = Yu.hconjugate()*Yu;
    Yd2 = Yd.hconjugate()*Yd;

    // Eq. (15) of Colangelo et al., 0807.0801 reduced according to Eq. (21)
    // Notice that Yukawa couplings are LR for Bertolini but RL for Colangelo

    mQ2a = a1*I+b1*Yu2;//+b2*Yd2+(b3*Yd2*Yu2+b3.conjugate()*Yu2*Yd2);
    mU2a = a2*I+b5*Yu*Yu.hconjugate();
    mD2a = a3*I+Yd*(b6+c4*Yu2)*Yd.hconjugate();
    Aua = Yu*(a4*I+c7*Yu2+b7*Yd2);
    Ada = Yd*(a5*I+b8*Yu2+c10*Yd2+c11*Yd2*Yu2);
    mqd(0,0) = md;
    mqd(1,1) = ms;
    mqd(2,2) = mb;
    mqd2=mqd*mqd;
    mqu(0,0) = mu;
    mqu(1,1) = mc;
    mqu(2,2) = mt;
    mqu2=mqu*mqu;
    cos2b = (1.-tanb*tanb)/(1.+tanb*tanb);

    //SUPER-CKM Basis

    mD2.assign(0, 0, mQ2a+mqd2+I*(mZ*mZ*cos2b*(-0.5+1./3.*sin2tw())));
    mD2.assign(0, 3, v1()/sqrt(2.)*Ada.hconjugate()-muH*mqd*tanb);
    mD2.assign(3, 0, v1()/sqrt(2.)*Ada-muH.conjugate()*mqd*tanb);
    mD2.assign(3, 3, mD2a+mqd2+I*(mZ*mZ*cos2b*(-1./3.*sin2tw())));

    mD2.eigensystem(Rd, Msd2);
    
    mU2.assign(0, 0, VCKM*mQ2a*VCKM.hconjugate()+mqu2+I*
    (mZ*mZ*cos2b*(0.5-2./3.*sin2tw())));
    mU2.assign(0, 3, v2()/sqrt(2.)*VCKM*Aua.hconjugate()-muH*mqu/tanb);
    mU2.assign(3, 0, v2()/sqrt(2.)*Aua*VCKM.hconjugate()-
    muH.conjugate()*mqu/tanb);
    mU2.assign(3, 3, mU2a+mqu2+I*(mZ*mZ*cos2b*(2./3.*sin2tw())));

    std::cout<< mU2 << std::endl;

    mU2.eigensystem(Ru, Msu2);
}
