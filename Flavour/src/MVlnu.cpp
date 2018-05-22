/* 
 * Copyright (C) 2014 HEPfit Collaboration
 *
 *
 * For the licensing terms see doc/COPYING.
 */

#include "StandardModel.h"
#include "MVlnu.h"
#include "std_make_vector.h"
#include "gslpp_function_adapter.h"
#include <gsl/gsl_sf_zeta.h>
#include <boost/bind.hpp>
#include <limits>
#include <TFitResult.h>
#include <gsl/gsl_sf_gegenbauer.h>
#include <gsl/gsl_sf_expint.h>

MVlnu::MVlnu(const StandardModel& SM_i, QCD::meson meson_i, QCD::meson vector_i, QCD::lepton lep_i)
: mySM(SM_i)
{
    lep = lep_i;
    meson = meson_i;
    vectorM = vector_i;
    
    w_J = gsl_integration_cquad_workspace_alloc (100); 
}

MVlnu::~MVlnu() {
}

std::vector<std::string> MVlnu::initializeMVlnuParameters()
{
    
    if (vectorM == StandardModel::D_star_P) mvlnuParameters = make_vector<std::string>()
        << "hA1w1" << "rho2" << "R1w1" << "R2w1";
    else {
        std::stringstream out;
        out << vectorM;
        throw std::runtime_error("MVlnu: vector " + out.str() + " not implemented");
    }
    
    mySM.initializeMeson(meson);
    mySM.initializeMeson(vectorM);
    return mvlnuParameters;
}

void MVlnu::updateParameters() 
{
    if (!mySM.getFlavour().getUpdateFlag(meson, vectorM, lep)) return;
  
    GF = mySM.getGF();
    Mlep = mySM.getLeptons(lep).getMass();
    Mnu = 0.; // neutrinos assumed to be massless
    MM = mySM.getMesons(meson).getMass();
    MV = mySM.getMesons(vectorM).getMass();
    width = mySM.getMesons(meson).computeWidth();
    w0 = (MM*MM+MV*MV)/(2.*MM*MV);
    RV = 2.*sqrt(MM*MV)/(MM+MV);
    mu_b = mySM.getMub();
    Mb = mySM.getQuarks(QCD::BOTTOM).getMass(); // add the PS b mass
    mb_pole = mySM.Mbar2Mp(Mb); /* Conversion to pole mass*/
    Mc = mySM.getQuarks(QCD::CHARM).getMass(); // add the PS b mass
    mc_pole = mySM.Mbar2Mp(Mc); /* Conversion to pole mass*/
    Vcb = mySM.getCKM().getV_cb();
    ale_mub = mySM.Ale(mu_b,FULLNLO);
    /* Amplitude propto 4*GF*Vij/sqrt(2) & kinematics requires 1/(2^9 pi^3 MB^3) */
    amplsq_factor = GF*GF*Vcb.abs2()/(64.*M_PI*M_PI*M_PI*MM*MM*MM);
    
    mySM.getFlavour().setUpdateFlag(meson, vectorM, lep, false);
    
    /* SM Wilson coefficients */
    CV_SM = 1./2.*(1.+ale_mub/M_PI*log(mySM.getMz()/mu_b));
    CA_SM = -CV_SM;

    /* SM + NP Wilson coefficients */
    CV = CV_SM;
    CVp = 0.;
    CA = CA_SM;
    CAp = 0.;
    CS = 0.;
    CSp = 0.;
    CP  = 0.;
    CPp  = 0.;
    C7  = 0.;
    C7p  = 0.;
    CT  = 0.;
    CTp  = 0.;
    
    switch (vectorM) {
        case StandardModel::D_star_P:
            hA1w1 = mySM.getOptionalParameter("hA1w1");
            rho2 = mySM.getOptionalParameter("rho2");
            R1w1 = mySM.getOptionalParameter("R1w1");
            R2w1 = mySM.getOptionalParameter("R2w1");
            break;
        default:
            std::stringstream out;
            out << vectorM;
            throw std::runtime_error("MVlnu: vector " + out.str() + " not implemented");
    }

    return;
    
}

/*******************************************************************************
 * Kinematic functions                                                          *
 * ****************************************************************************/

double MVlnu::lambda_half(double a, double b, double c) 
{   
    return sqrt(a*a+b*b+c*c-2.*(a*b+a*c+b*c));
}

/*******************************************************************************
 * Form factors                                                                *
 * ****************************************************************************/

double MVlnu::hA1(double q2) 
{
    double w = w0-q2/(2.*MM*MV);
    double z = (sqrt(w+1.)-sqrt(2.))/(sqrt(w+1.)+sqrt(2.));
    return hA1w1*(1.-8.*rho2*z+(53.*rho2-15.)*z*z-(231.*rho2-91.)*z*z*z);
}

double MVlnu::R1(double q2) 
{
    double w = w0-q2/(2.*MM*MV);
    return R1w1-0.12*(w-1.)+0.05*(w-1.)*(w-1.);
}

double MVlnu::R2(double q2) 
{
    double w = w0-q2/(2.*MM*MV);
    return R2w1+0.11*(w-1.)-0.06*(w-1.)*(w-1.);
}

double MVlnu::R0(double q2) 
{
    /* this can be derived from A0, A1 and A2 q2=0 relations */
    double R2q2at0 = R2(0.);
    double R0q2at0 = (MM+MV-(MM-MV)*R2q2at0)/(2.*MV);
    double R0w1 = R0q2at0+0.11*(w0-1.)-0.01*(w0-1.)*(w0-1.);
    double w = w0-q2/(2.*MM*MV);
    return R0w1-0.11*(w-1.)+0.01*(w-1.)*(w-1.);
}

double MVlnu::V(double q2) 
{
    return R1(q2)/RV*hA1(q2);
}

double MVlnu::A0(double q2) 
{
    return R0(q2)/RV*hA1(q2);
}

double MVlnu::A1(double q2) 
{
    double w = w0-q2/(2.*MM*MV);
    return (w+1.)*RV/2.*hA1(q2);
}

double MVlnu::A2(double q2) 
{
    return R2(q2)/RV*hA1(q2);
}

double MVlnu::T1(double q2) 
{
    double delta_T1 = 0.;
    return (Mb+Mc)/(MM+MV)*V(q2)*(1.+delta_T1);
}

double MVlnu::T2(double q2) 
{
    double delta_T2 = 0.;
    return (Mb-Mc)/(MM-MV)*A1(q2)*(1.+delta_T2);
}

double MVlnu::A12(double q2) 
{
    return (A1(q2)*(MM+MV)*(MM+MV)*(MM*MM-MV*MV-q2)- 
            A2(q2)*(MM*MM*MM*MM+(MV*MV-q2)*(MV*MV-q2)- 
            2.*MM*MM*(MV*MV+q2)))/(16.*MM*MV*MV*(MM+MV));
}

double MVlnu::T23(double q2) 
{
    double delta_T23 = 0.;
    return ((Mb-Mc)*((MM-MV)*(MM-MV)-q2)*((MM+MV)*(MM+MV)-q2)*A0(q2)+ 
            8*MM*MV*(MV*MV-MM*MM)*A12(q2))/(4.*MM*(MV-MM)*MV*q2)*(1.+delta_T23);
}

/*******************************************************************************
 * Helicity amplitudes  (normalization such that all H \propto (mass scale)^-1)                                                       *
 * ****************************************************************************/

gslpp::complex MVlnu::HV0(double q2) 
{
    return 4.*gslpp::complex::i()*MM*MV/(sqrt(q2)*(MM+MV))*((CV-CVp)*(MM+MV)*A12(q2)+Mb*(C7-C7p)*T23(q2));   
}

gslpp::complex MVlnu::HVp(double q2) 
{
    return gslpp::complex::i()*((((CV + CVp)*lambda_half(MM*MM,MV*MV,q2)*V(q2)-(MM+MV)*(MM+MV)*(CV-CVp)*A1(q2)))/(2.*(MM+MV))
            + (Mb/q2)*((C7+C7p)*lambda_half(MM*MM,MV*MV,q2)*T1(q2)-(C7-C7p)*(MM*MM-MV*MV)*T2(q2)));
}

gslpp::complex MVlnu::HVm(double q2) 
{
    return gslpp::complex::i()*(((-(CV + CVp)*lambda_half(MM*MM,MV*MV,q2)*V(q2)-(MM+MV)*(MM+MV)*(CV-CVp)*A1(q2)))/(2.*(MM+MV))
            + (Mb/q2)*(-(C7+C7p)*lambda_half(MM*MM,MV*MV,q2)*T1(q2)-(C7-C7p)*(MM*MM-MV*MV)*T2(q2)));   
}

gslpp::complex MVlnu::HAp(double q2) 
{
    return gslpp::complex::i()*((CA+CAp)*lambda_half(MM*MM,MV*MV,q2)*V(q2)-(MM+MV)*(MM+MV)*(CA-CAp)*A1(q2))/(2.*(MM+MV));   
}

gslpp::complex MVlnu::HAm(double q2) 
{
    return gslpp::complex::i()*(-(CA+CAp)*lambda_half(MM*MM,MV*MV,q2)*V(q2)-(MM+MV)*(MM+MV)*(CA-CAp)*A1(q2))/(2.*(MM+MV));   
}

gslpp::complex MVlnu::HA0(double q2) 
{
    return 4.*gslpp::complex::i()*MV*MM/(sqrt(q2))*(CA-CAp)*A12(q2);   
}

gslpp::complex MVlnu::HP(double q2) 
{
    return gslpp::complex::i()*lambda_half(MM*MM,MV*MV,q2)/2.*((CP-CPp)/(Mb+Mc)+(Mlep+Mnu)/q2*(CA-CAp))*A0(q2);   
}

gslpp::complex MVlnu::HS(double q2) 
{
    return gslpp::complex::i()*lambda_half(MM*MM,MV*MV,q2)/2.*((CS-CSp)/(Mb+Mc)+(Mlep+Mnu)/q2*(CV-CVp))*A0(q2);
}

gslpp::complex MVlnu::HT0(double q2) 
{
    return 2.*sqrt(2.)*(MM*MV)/(MM+MV)*(CT+CTp)*T23(q2);   
}

gslpp::complex MVlnu::HT0t(double q2) 
{
    return 2.*(MM*MV)/(MM+MV)*(CT-CTp)*T23(q2);  
}

gslpp::complex MVlnu::HTp(double q2) 
{
    return ((CT-CTp)*lambda_half(MM*MM,MV*MV,q2)*T1(q2)-(CT+CTp)*(MM*MM-MV*MV)*T2(q2))/(sqrt(2.*q2));   
}

gslpp::complex MVlnu::HTpt(double q2) 
{
    return ((CT+CTp)*lambda_half(MM*MM,MV*MV,q2)*T1(q2)-(CT-CTp)*(MM*MM-MV*MV)*T2(q2))/(sqrt(2.*q2));  
}

gslpp::complex MVlnu::HTm(double q2)
{
    return (-(CT-CTp)*lambda_half(MM*MM,MV*MV,q2)*T1(q2)-(CT+CTp)*(MM*MM-MV*MV)*T2(q2))/(sqrt(2.*q2));   
}

gslpp::complex MVlnu::HTmt(double q2) 
{
    return (-(CT+CTp)*lambda_half(MM*MM,MV*MV,q2)*T1(q2)-(CT-CTp)*(MM*MM-MV*MV)*T2(q2))/(sqrt(2.*q2));  
}

/*******************************************************************************
 * Generalized angular coefficients  (see 1506.03970)                          *
 * ****************************************************************************/

gslpp::complex MVlnu::G000(double q2) 
{
    double lambda_MM = lambda_half(MM*MM,MV*MV,q2);
    double lambda_lep = lambda_half(Mlep*Mlep,Mnu*Mnu,q2);
    double lambda_lep2 = lambda_lep*lambda_lep;
    double Elep = sqrt(Mlep*Mlep+lambda_lep2/(4.*q2));
    double Enu = sqrt(Mnu*Mnu+lambda_lep2/(4.*q2));
    double Gprefactor = lambda_MM*lambda_lep/q2;
    
    return Gprefactor*(4./9.*(3*Elep*Enu+lambda_lep2/(4.*q2))*(HVp(q2).abs2()+HVm(q2).abs2()+HV0(q2).abs2()+HAp(q2).abs2()+HAm(q2).abs2()+HA0(q2).abs2())+ 
            4.*Mlep*Mnu/3.*(HVp(q2).abs2()+HVm(q2).abs2()+HV0(q2).abs2()-HAp(q2).abs2()-HAm(q2).abs2()-HA0(q2).abs2())+
            4./3.*(Elep*Enu-Mlep*Mnu+lambda_lep2/(4.*q2))*HS(q2).abs2()+4/3*(Elep*Enu+Mlep*Mnu+lambda_lep2/(4*q2))*HP(q2).abs2()+
            16./9.*(3.*(Elep*Enu+Mlep*Mnu)-lambda_lep2/(4.*q2))*(HTpt(q2).abs2()+HTmt(q2).abs2()+HT0t(q2).abs2())+ 
            8./9.*(3*(Elep*Enu-Mlep*Mnu)-lambda_lep2/(4.*q2))*(HTpt(q2).abs2()+HTmt(q2).abs2()+HT0t(q2).abs2())+ 
            16./3.*(Mlep*Enu+Mnu*Elep)*(HVp(q2)*HTpt(q2).conjugate()+HVm(q2)*HTmt(q2).conjugate()+HV0(q2)*HT0t(q2).conjugate()).imag()+
            8*sqrt(2.)/3.*(Mlep*Enu-Mnu*Elep)*(HAp(q2)*HTp(q2).conjugate()+HAm(q2)*HTm(q2).conjugate()+HA0(q2)*HT0(q2).conjugate()).imag());
}

gslpp::complex MVlnu::G010(double q2) 
{
    double lambda_MM = lambda_half(MM*MM,MV*MV,q2);
    double lambda_lep = lambda_half(Mlep*Mlep,Mnu*Mnu,q2);
    double Gprefactor = lambda_MM*lambda_lep/q2;
    
    return Gprefactor*4./3.*lambda_lep*((HVp(q2)*HAp(q2).conjugate()-HVm(q2)*HAm(q2).conjugate()).real()+
            (2.*sqrt(2.))/q2*(Mlep*Mlep-Mnu*Mnu)*(HTp(q2)*HTpt(q2).conjugate()-HTm(q2)*HTmt(q2).conjugate()).real()+
            2.*(Mlep+Mnu)/sqrt(q2)*(HAp(q2)*HTpt(q2).conjugate()-HAm(q2)*HTmt(q2).conjugate()).imag()+
            sqrt(2.)*(Mlep-Mnu)/sqrt(q2)*(HVp(q2)*HTp(q2).conjugate()-HVm(q2)*HTm(q2).conjugate()).imag()-
            (Mlep-Mnu)/sqrt(q2)*(HA0(q2)*HP(q2).conjugate()).real()-(Mlep+Mnu)/sqrt(q2)*(HV0(q2)*HS(q2).conjugate()).real()+
            (sqrt(2.)*HT0(q2)*HP(q2).conjugate()+2.*HT0t(q2)*HS(q2).conjugate()).imag());
    
}

gslpp::complex MVlnu::G020(double q2) 
{
    double lambda_MM = lambda_half(MM*MM,MV*MV,q2);
    double lambda_lep = lambda_half(Mlep*Mlep,Mnu*Mnu,q2);
    double lambda_lep2 = lambda_lep*lambda_lep;
    double Gprefactor = lambda_MM*lambda_lep/q2;
    
    return -Gprefactor*2./9.*lambda_lep2/q2*((-HVp(q2).abs2()-HVm(q2).abs2()+2.*HV0(q2).abs2()-HAp(q2).abs2()-HAm(q2).abs2()+2.*HA0(q2).abs2())- 
            2.*(2.*HT0(q2).abs2()-HTp(q2).abs2()-HTm(q2).abs2())-4.*(2.*HT0t(q2).abs2()-HTpt(q2).abs2()-HTmt(q2).abs2()));    
}

gslpp::complex MVlnu::G200(double q2) 
{
    double lambda_MM = lambda_half(MM*MM,MV*MV,q2);
    double lambda_lep = lambda_half(Mlep*Mlep,Mnu*Mnu,q2);
    double lambda_lep2 = lambda_lep*lambda_lep;
    double Elep = sqrt(Mlep*Mlep+lambda_lep2/(4.*q2));
    double Enu = sqrt(Mnu*Mnu+lambda_lep2/(4.*q2));
    double Gprefactor = lambda_MM*lambda_lep/q2;

    return Gprefactor*(-4./9.*(3.*Elep*Enu+lambda_lep2/(4.*q2))*(HVp(q2).abs2()+HVm(q2).abs2()-2.*HV0(q2).abs2()+HAp(q2).abs2()+HAm(q2).abs2()-2.*HA0(q2).abs2())-
            4./3.*Mlep*Mnu*(HVp(q2).abs2()+HVm(q2).abs2()-2.*HV0(q2).abs2()-HAp(q2).abs2()-HAm(q2).abs2()+2.*HA0(q2).abs2())+
            8./3.*(Elep*Enu-Mlep*Mnu+lambda_lep2/(4.*q2))*HS(q2).abs2()+8./3.*(Elep*Enu+Mlep*Mnu+lambda_lep2/(4.*q2))*HP(q2).abs2()-
            16./9.*(3.*(Elep*Enu+Mlep*Mnu)-lambda_lep2/(4.*q2))*(HTpt(q2).abs2()+HTmt(q2).abs2()-2.*HT0t(q2).abs2())-8./9.*(3.*(Elep*Enu-Mlep*Mnu)-lambda_lep2/(4.*q2))*
            (HTp(q2).abs2()+HTm(q2).abs2()-2*HT0(q2).abs2())-16./3.*(Mlep*Enu+Mnu*Elep)*(HVp(q2)*HTpt(q2).conjugate()+HVm(q2)*HTmt(q2).conjugate()-2.*HV0(q2)*HT0t(q2).conjugate()).imag()-
            8.*sqrt(2.)/3.*(Mlep*Enu-Mnu*Elep)*(HAp(q2)*HTp(q2).conjugate()+HAm(q2)*HTm(q2).conjugate()-2.*HA0(q2)*HT0(q2).conjugate()).imag());
}

gslpp::complex MVlnu::G210(double q2) 
{
    double lambda_MM = lambda_half(MM*MM,MV*MV,q2);
    double lambda_lep = lambda_half(Mlep*Mlep,Mnu*Mnu,q2);
    double Gprefactor = lambda_MM*lambda_lep/q2;
    
    return -Gprefactor*4.*lambda_lep/3.*((HVp(q2)*HAp(q2).conjugate()-HVm(q2)*HAm(q2).conjugate()).real()+ 
            2*sqrt(2.)*(Mlep*Mlep-Mnu*Mnu)/q2*(HTp(q2)*HTpt(q2).conjugate()-HTm(q2)*HTmt(q2).conjugate()).real()+ 
            2*(Mlep+Mnu)/sqrt(q2)*(HAp(q2)*HTpt(q2).conjugate()-HAm(q2)*HTmt(q2).conjugate()).imag()+ 
            sqrt(2.)*(Mlep-Mnu)/sqrt(q2)*(HVp(q2)*HTp(q2).conjugate()-HVm(q2)*HTm(q2).conjugate()).imag()+ 
            2*(Mlep-Mnu)/sqrt(q2)*(HA0(q2)*HP(q2).conjugate()).real()+2.*(Mlep+Mnu)/sqrt(q2)*(HV0(q2)*HS(q2).conjugate()).real()- 
            2*(sqrt(2.)*HT0(q2)*HP(q2).conjugate()+2.*HT0t(q2)*HS(q2).conjugate()).imag());
}

gslpp::complex MVlnu::G220(double q2) 
{
    double lambda_MM = lambda_half(MM*MM,MV*MV,q2);
    double lambda_lep = lambda_half(Mlep*Mlep,Mnu*Mnu,q2);
    double lambda_lep2 = lambda_lep*lambda_lep;
    double Gprefactor = lambda_MM*lambda_lep/q2;
    
    return -Gprefactor*2./9.*lambda_lep2/q2*(HVp(q2).abs2()+HVm(q2).abs2()+4.*HV0(q2).abs2()+HAp(q2).abs2()+ 
            HAm(q2).abs2()+4.*HA0(q2).abs2()-2.*(HTp(q2).abs2()+HTm(q2).abs2()+4.*HT0(q2).abs2())-
            4.*(HTpt(q2).abs2()+HTmt(q2).abs2()+4.*HT0t(q2).abs2()));
}

gslpp::complex MVlnu::G211(double q2) 
{
    double lambda_MM = lambda_half(MM*MM,MV*MV,q2);
    double lambda_lep = lambda_half(Mlep*Mlep,Mnu*Mnu,q2);
    double Gprefactor = lambda_MM*lambda_lep/q2;
    
    return Gprefactor*4./sqrt(3.)*lambda_lep*(HVp(q2)*HA0(q2).conjugate()+HAp(q2)*HV0(q2).conjugate()- 
            HV0(q2)*HAm(q2).conjugate()-HA0(q2)*HVm(q2).conjugate()+(Mlep+Mnu)/sqrt(q2)*(HVp(q2)*HS(q2).conjugate()+HS(q2)*HVm(q2).conjugate())- 
            gslpp::complex::i()*sqrt(2.)*(HP(q2)*HTm(q2).conjugate()-HTp(q2)*HP(q2).conjugate()+sqrt(2.)*(HS(q2)*HTmt(q2).conjugate()-HTpt(q2)*HS(q2).conjugate()))+ 
            (Mlep-Mnu)/sqrt(q2)*(HAp(q2)*HP(q2).conjugate()+HP(q2)*HAm(q2).conjugate())- 
            gslpp::complex::i()*2.*(Mlep+Mnu)/sqrt(q2)*(HAp(q2)*HT0t(q2).conjugate()+HT0t(q2)*HAm(q2).conjugate()-HTpt(q2)*HA0(q2).conjugate()-HA0(q2)*HTmt(q2).conjugate())-
            gslpp::complex::i()*sqrt(2.)*(Mlep-Mnu)/sqrt(q2)*(HVp(q2)*HT0(q2).conjugate()+HT0(q2)*HVm(q2).conjugate()-HTp(q2)*HV0(q2).conjugate()-HV0(q2)*HTm(q2).conjugate())+
            2.*sqrt(2.)*(Mlep*Mlep-Mnu*Mnu)/q2*(HTp(q2)*HT0t(q2).conjugate()+HTpt(q2)*HT0(q2).conjugate()-HT0(q2)*HTmt(q2).conjugate()-HT0t(q2)*HTm(q2).conjugate()));
}

gslpp::complex MVlnu::G221(double q2) 
{
    double lambda_MM = lambda_half(MM*MM,MV*MV,q2);
    double lambda_lep = lambda_half(Mlep*Mlep,Mnu*Mnu,q2);
    double lambda_lep2 = lambda_lep*lambda_lep;
    double Gprefactor = lambda_MM*lambda_lep/q2;
    
    return Gprefactor*4./3.*lambda_lep2/q2*(HVp(q2)*HV0(q2).conjugate()+HV0(q2)*HVm(q2).conjugate()+ 
            HAp(q2)*HA0(q2).conjugate()+HA0(q2)*HAm(q2).conjugate()-2.*(HTp(q2)*HT0(q2).conjugate()+ 
            HT0(q2)*HTm(q2).conjugate()+2.*(HTpt(q2)*HT0t(q2).conjugate()+HT0t(q2)*HTmt(q2).conjugate())));   
}

gslpp::complex MVlnu::G222(double q2) 
{
    double lambda_MM = lambda_half(MM*MM,MV*MV,q2);
    double lambda_lep = lambda_half(Mlep*Mlep,Mnu*Mnu,q2);
    double lambda_lep2 = lambda_lep*lambda_lep;
    double Gprefactor = lambda_MM*lambda_lep/q2;
    
    return -Gprefactor*8./3.*lambda_lep2/q2*(HVp(q2)*HVm(q2).conjugate()+HAp(q2)*HAm(q2).conjugate()- 
            2.*(HTp(q2)*HTm(q2).conjugate()+2.*HTpt(q2)*HTmt(q2).conjugate())); 
}

/***************************************************************************
 * 12 independent J angular coefficients  (see again 1506.03970)           *
 * ************************************************************************/

double MVlnu::J1s(double q2) 
{
    return amplsq_factor*(8.*G000(q2)+2.*G020(q2)-4.*G200(q2)-G220(q2)).real()/3.;
}

double MVlnu::J1c(double q2) 
{
    return amplsq_factor*(8.*G000(q2)+2.*G020(q2)+8.*G200(q2)+2.*G220(q2)).real()/3.;
}

double MVlnu::J2s(double q2) 
{
    return amplsq_factor*(2.*G020(q2)-G220(q2)).real();
}
        
double MVlnu::J2c(double q2) 
{
    return amplsq_factor*(2.*(G020(q2)+G220(q2))).real();
}

double MVlnu::J3(double q2) 
{
    return amplsq_factor*(G222(q2).real());
}

double MVlnu::J4(double q2) 
{
    return -amplsq_factor*(G221(q2).real());
}

double MVlnu::J5(double q2) 
{
    return amplsq_factor*(2.*G211(q2).real()/sqrt(3.));
}

double MVlnu::J6s(double q2) 
{
    return -amplsq_factor*(4.*(2.*G010(q2)-G210(q2)).real()/3.);
}

double MVlnu::J6c(double q2) 
{
    return -amplsq_factor*(8.*(G010(q2)+G210(q2)).real()/3.);
}

double MVlnu::J7(double q2) 
{
    return -amplsq_factor*(2.*sqrt(3.)*(G211(q2).imag())/3.);
}

double MVlnu::J8(double q2) 
{
    return amplsq_factor*(G221(q2).imag());
}

double MVlnu::J9(double q2) 
{
    return -amplsq_factor*(G222(q2).imag());
}

/***************************************************************************
 * Integration of angular coefficients Js                                  *
 * ************************************************************************/

double MVlnu::integrateJ(int i, double q2_min, double q2_max) 
{
    updateParameters();
    
    old_handler = gsl_set_error_handler_off();
    
    switch (i) {
        case 0:
                FJ = convertToGslFunction(boost::bind(&MVlnu::J1s, &(*this), _1));
                if (gsl_integration_cquad(&FJ, q2_min, q2_max, 1.e-2, 1.e-1, w_J, &J_res, &J_err, NULL) != 0) std::numeric_limits<double>::quiet_NaN();
                gsl_set_error_handler(old_handler);
                return J_res;
                break;
        case 1:
                FJ = convertToGslFunction(boost::bind(&MVlnu::J1c, &(*this), _1));
                if (gsl_integration_cquad(&FJ, q2_min, q2_max, 1.e-2, 1.e-1, w_J, &J_res, &J_err, NULL) != 0) std::numeric_limits<double>::quiet_NaN();
                gsl_set_error_handler(old_handler);
                return J_res;
                break;
        case 2:
                FJ = convertToGslFunction(boost::bind(&MVlnu::J2s, &(*this), _1));
                if (gsl_integration_cquad(&FJ, q2_min, q2_max, 1.e-2, 1.e-1, w_J, &J_res, &J_err, NULL) != 0) std::numeric_limits<double>::quiet_NaN();
                gsl_set_error_handler(old_handler);
                return J_res;
                break;
        case 3:
                FJ = convertToGslFunction(boost::bind(&MVlnu::J2c, &(*this), _1));
                if (gsl_integration_cquad(&FJ, q2_min, q2_max, 1.e-2, 1.e-1, w_J, &J_res, &J_err, NULL) != 0) std::numeric_limits<double>::quiet_NaN();
                return J_res;
                break;
        case 4:
                FJ = convertToGslFunction(boost::bind(&MVlnu::J3, &(*this), _1));
                if (gsl_integration_cquad(&FJ, q2_min, q2_max, 1.e-2, 1.e-1, w_J, &J_res, &J_err, NULL) != 0) std::numeric_limits<double>::quiet_NaN();
                gsl_set_error_handler(old_handler);
                return J_res;
                break;
        case 5:
                FJ = convertToGslFunction(boost::bind(&MVlnu::J4, &(*this), _1));
                if (gsl_integration_cquad(&FJ, q2_min, q2_max, 1.e-2, 1.e-1, w_J, &J_res, &J_err, NULL) != 0) std::numeric_limits<double>::quiet_NaN();
                return J_res;
                break;
        case 6:
                FJ = convertToGslFunction(boost::bind(&MVlnu::J5, &(*this), _1));
                if (gsl_integration_cquad(&FJ, q2_min, q2_max, 1.e-2, 1.e-1, w_J, &J_res, &J_err, NULL) != 0) std::numeric_limits<double>::quiet_NaN();
                return J_res;
                break;
        case 7:
                FJ = convertToGslFunction(boost::bind(&MVlnu::J6s, &(*this), _1));
                if (gsl_integration_cquad(&FJ, q2_min, q2_max, 1.e-2, 1.e-1, w_J, &J_res, &J_err, NULL) != 0) std::numeric_limits<double>::quiet_NaN();
                gsl_set_error_handler(old_handler);
                return J_res;
                break;
        case 9:
                FJ = convertToGslFunction(boost::bind(&MVlnu::J6c, &(*this), _1));
                if (gsl_integration_cquad(&FJ, q2_min, q2_max, 1.e-2, 1.e-1, w_J, &J_res, &J_err, NULL) != 0) std::numeric_limits<double>::quiet_NaN();
                return J_res;
                break;
        case 10:
                FJ = convertToGslFunction(boost::bind(&MVlnu::J7, &(*this), _1));
                if (gsl_integration_cquad(&FJ, q2_min, q2_max, 1.e-2, 1.e-1, w_J, &J_res, &J_err, NULL) != 0) std::numeric_limits<double>::quiet_NaN();
                gsl_set_error_handler(old_handler);
                return J_res;
                break;
        case 11:
                FJ = convertToGslFunction(boost::bind(&MVlnu::J8, &(*this), _1));
                if (gsl_integration_cquad(&FJ, q2_min, q2_max, 1.e-2, 1.e-1, w_J, &J_res, &J_err, NULL) != 0) std::numeric_limits<double>::quiet_NaN();
                gsl_set_error_handler(old_handler);
                return J_res;
                break;
        case 12:
                FJ = convertToGslFunction(boost::bind(&MVlnu::J9, &(*this), _1));
                if (gsl_integration_cquad(&FJ, q2_min, q2_max, 1.e-2, 1.e-1, w_J, &J_res, &J_err, NULL) != 0) std::numeric_limits<double>::quiet_NaN();
                gsl_set_error_handler(old_handler);
                return J_res;
                break;
        default:
            gsl_set_error_handler(old_handler);
            std::stringstream out;
            out << i;
            throw std::runtime_error("MVlnu::integrateJ: index " + out.str() + " not implemented");
    }
}

double MVlnu::dGammadw(double q2)
{
    updateParameters();
    
    return 3./4.*(2.*J1s(q2)+J1c(q2))-1./4.*(2.*J2s(q2)+J2c(q2));
}

double MVlnu::getDeltaGammaDeltaw(double w_min, double w_max) 
{
    updateParameters();
    
    double q2_min = (2.*MM*MV)*(w0-w_max); // min is Mlep*Mlep;
    double q2_max = (2.*MM*MV)*(w0-w_min); // max is (MM-MV)*(MM-MV);
    
    double intJ1s = integrateJ(1,q2_min,q2_max);
    double intJ1c = integrateJ(2,q2_min,q2_max);
    double intJ2s = integrateJ(3,q2_min,q2_max);
    double intJ2c = integrateJ(4,q2_min,q2_max);

    return 3./4.*(2.*intJ1s+intJ1c)-1./4.*(2.*intJ2s+intJ2c);
}

double MVlnu::dGammadcldq2(double q2, double cl)
{
    updateParameters();
    
    return 3./8.*((J1s(q2)+2.*J1c(q2))+cl*(J6s(q2)+2.*J6c(q2))+(2.*cl*cl-1.)*(J2s(q2)+2.*J2c(q2)));
}

double MVlnu::dGammadcl(double cl)
{
    updateParameters();
    
    double q2_min = Mlep*Mlep;
    double q2_max = (MM-MV)*(MM-MV);
    
    double intJ1s = integrateJ(1,q2_min,q2_max);
    double intJ1c = integrateJ(2,q2_min,q2_max);
    double intJ2s = integrateJ(3,q2_min,q2_max);
    double intJ2c = integrateJ(4,q2_min,q2_max);
    double intJ6s = integrateJ(6,q2_min,q2_max);
    double intJ6c = integrateJ(7,q2_min,q2_max);

    return 3./8.*((intJ1s+2.*intJ1c)+cl*(intJ6s+2.*intJ6c)+(2.*cl*cl-1.)*(intJ2s+2.*intJ2c));
}

double MVlnu::getDeltaGammaDeltacl(double cl_min, double cl_max)
{
    updateParameters();
    
    double q2_min = Mlep*Mlep;
    double q2_max = (MM-MV)*(MM-MV);
    
    double intJ1s = integrateJ(1,q2_min,q2_max);
    double intJ1c = integrateJ(2,q2_min,q2_max);
    double intJ2s = integrateJ(3,q2_min,q2_max);
    double intJ2c = integrateJ(4,q2_min,q2_max);
    double intJ6s = integrateJ(6,q2_min,q2_max);
    double intJ6c = integrateJ(7,q2_min,q2_max);
    
    return 3./8.*((cl_max-cl_min)*(intJ1s+2.*intJ1c)+
            (cl_max*cl_max-cl_min*cl_min)/2.*(intJ6s+2.*intJ6c)+
            (2./3.*(cl_max*cl_max*cl_max-cl_min*cl_min*cl_min)-(cl_max-cl_min))*(intJ2s+2.*intJ2c));
}

double MVlnu::dGammadcVdq2(double q2, double cl)
{
    updateParameters();
    
    return 3./8.*((J1s(q2)+2.*J1c(q2))+cl*(J6s(q2)+2.*J6c(q2))+(2.*cl*cl-1.)*(J2s(q2)+2.*J2c(q2)));
}

double MVlnu::dGammadcV(double cV)
{
    updateParameters();
    
    double q2_min = Mlep*Mlep;
    double q2_max = (MM-MV)*(MM-MV);
    
    double intJ1s = integrateJ(1,q2_min,q2_max);
    double intJ1c = integrateJ(2,q2_min,q2_max);
    double intJ2s = integrateJ(3,q2_min,q2_max);
    double intJ2c = integrateJ(4,q2_min,q2_max);

    return 3./8.*(cV*cV*(3.*intJ1c-intJ2c)+(1.-cV*cV)*(3.*intJ1s-intJ2s));
}

double MVlnu::getDeltaGammaDeltacV(double cV_min, double cV_max)
{
    updateParameters();
    
    double q2_min = Mlep*Mlep;
    double q2_max = (MM-MV)*(MM-MV);
    
    double intJ1s = integrateJ(1,q2_min,q2_max);
    double intJ1c = integrateJ(2,q2_min,q2_max);
    double intJ2s = integrateJ(3,q2_min,q2_max);
    double intJ2c = integrateJ(4,q2_min,q2_max);
    
    return 3./8.*((cV_max*cV_max*cV_max-cV_min*cV_min*cV_min)*(3.*intJ1c-intJ2c)+
            ((cV_max-cV_min)-(cV_max*cV_max*cV_max-cV_min*cV_min*cV_min)/3.)*(3.*intJ1s-intJ2s));
}

double MVlnu::dGammadchidq2(double q2, double chi)
{
    updateParameters();
    
    return (3.*J1c(q2)+6.*J1s(q2)-J2c(q2)-2.*J2s(q2))/8./M_PI+
            cos(2.*chi)/2./M_PI*J3(q2)+sin(2.*chi)/2./M_PI*J9(q2);
}

double MVlnu::dGammadchi(double chi)
{
    updateParameters();
    
    double q2_min = Mlep*Mlep;
    double q2_max = (MM-MV)*(MM-MV);
    
    double intJ1s = integrateJ(1,q2_min,q2_max);
    double intJ1c = integrateJ(2,q2_min,q2_max);
    double intJ2s = integrateJ(3,q2_min,q2_max);
    double intJ2c = integrateJ(4,q2_min,q2_max);
    double intJ3 = integrateJ(5,q2_min,q2_max);
    double intJ9 = integrateJ(12,q2_min,q2_max);

    return ((3.*intJ1c+6.*intJ1s-intJ2c-2.*intJ2s)/4.+
            cos(2.*chi)*intJ3+sin(2.*chi)*intJ9)/2./M_PI;
}

double MVlnu::getDeltaGammaDeltachi(double chi_min, double chi_max)
{
    updateParameters();
    
    double q2_min = Mlep*Mlep;
    double q2_max = (MM-MV)*(MM-MV);
    
    double intJ1s = integrateJ(1,q2_min,q2_max);
    double intJ1c = integrateJ(2,q2_min,q2_max);
    double intJ2s = integrateJ(3,q2_min,q2_max);
    double intJ2c = integrateJ(4,q2_min,q2_max);
    double intJ3 = integrateJ(5,q2_min,q2_max);
    double intJ9 = integrateJ(12,q2_min,q2_max);

    return ((chi_max-chi_min)*(3.*intJ1c+6.*intJ1s-intJ2c-2.*intJ2s)/4.+
            (sin(2.*chi_max)-sin(2.*chi_min))/2.*intJ3-
            (cos(2.*chi_max)-cos(2.*chi_min))/2.*intJ9)/2./M_PI;
}