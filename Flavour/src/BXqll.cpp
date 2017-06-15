/* 
 * Copyright (C) 2016 HEPfit Collaboration
 *
 *
 * For the licensing terms see doc/COPYING.
 */

#include <gsl/gsl_sf.h>
#include <gslpp.h>
#include <boost/bind.hpp>
//#include <limits>
#include "BXqll.h"
#include "StandardModel.h"
#include "F_1.h"
#include "F_2.h"
#include "gslpp_function_adapter.h"


BXqll::BXqll(const StandardModel& SM_i, QCD::quark quark_i, QCD::lepton lep_i)
: mySM(SM_i), myF_1(*(new F_1())), myF_2(*(new F_2()))
{    
    lep = lep_i;
    quark = quark_i;
    w_Rquark = gsl_integration_cquad_workspace_alloc(100);
}

BXqll::~BXqll() 
{
}

void BXqll::updateParameters() 
{
    CF = mySM.getCF();
    GF = mySM.getGF();
    ale = mySM.getAle();
    Mlep = mySM.getLeptons(lep).getMass();
    mu_b = mySM.getMub();
    mu_c = mySM.getMuc();
    Mb = mySM.getQuarks(QCD::BOTTOM).getMass(); // add the PS b mass
    Mc = mySM.getQuarks(QCD::CHARM).getMass();
    Ms = mySM.getQuarks(QCD::STRANGE).getMass();
    MW = mySM.Mw();
    abslambdat_over_Vcb = mySM.computelamt_s().abs()/mySM.getCKM().V_cb().abs();
    Vts_over_Vcb = mySM.getCKM().V_ts().abs()/mySM.getCKM().V_cb().abs();    
    muh = mu_b/Mb;
    alsmu = mySM.Als(mu_b);
    alsmuc = mySM.Als(mu_c);
    Mb_pole = mySM.Mbar2Mp(Mb, FULLNNLO);
    //Mc_pole = mySM.Mbar2Mp(Mc, FULLNNLO); //*** Mbar2Mp does not receive Mc ***/
    Mc_pole = Mc*(1.+4.*alsmuc/3./M_PI+alsmuc*alsmuc/M_PI/M_PI*(-1.0414*(1.-4.*Ms/3.*Mc)+13.4434));
    z = Mc_pole*Mc_pole/Mb_pole/Mb_pole; //****** Must be pole masses ****/
    lambda_2 = (5.32483*5.32483-5.27955*5.27955)/4.;

    //ISIDORI VALUES
    z = 0.29*0.29;
    mu_b = 5.0;
    Mb = 4.9;
    Mtau = 1.77;
    muh = mu_b/Mb;
    ale = 0.0078125;
    abslambdat_over_Vcb = 0.97;
    Vts_over_Vcb = 0.97;
    alsmu = 0.215;

    allcoeff = mySM.getFlavour().ComputeCoeffBMll(mu_b, lep); //check the mass scale, scheme fixed to NDR
    allcoeffprime = mySM.getFlavour().ComputeCoeffprimeBMll(mu_b, lep); //check the mass scale, scheme fixed to NDR

    for(int ord=LO; ord <= NLO; ord++)
    {
        C_1[ord] = (*(allcoeff[ord]))(0);
        C_1L_bar[ord] = (*(allcoeff[ord]))(0) / 2.;
        C_2[ord] = (*(allcoeff[ord]))(1);
        C_2L_bar[ord] = (*(allcoeff[ord]))(1) - (*(allcoeff[LO]))(0) / 6.;
        C_3[ord] = (*(allcoeff[ord]))(2);
        C_4[ord] = (*(allcoeff[ord]))(3);
        C_5[ord] = (*(allcoeff[ord]))(4);
        C_6[ord] = (*(allcoeff[ord]))(5);
        C_7[ord] = (*(allcoeff[ord]))(6);
        C_8L[ord] = (*(allcoeff[ord]))(7);
        C_9[ord] = (*(allcoeff[ord]))(8);
        C_10[ord] = (*(allcoeff[ord]))(9);
        C_S[ord] = (*(allcoeff[ord]))(10);
        C_P[ord] = (*(allcoeff[ord]))(11);

        C_7p[ord] = (*(allcoeffprime[ord]))(6);
        C_9p[ord] = (*(allcoeffprime[ord]))(8);
        C_10p[ord] = (*(allcoeffprime[ord]))(9);
        C_Sp[ord] = (*(allcoeffprime[ord]))(10);
        C_Pp[ord] = (*(allcoeffprime[ord]))(11);
        T_9[ord] = 4./3.*C_1[ord]+C_2[ord]+6.*C_3[ord]+60.*C_5[ord];
        U_9[ord] = -7./2.*C_3[ord]-2./3.*C_4[ord]-38.*C_5[ord]-32./3.*C_6[ord];
        W_9[ord] = -1./2.*C_3[ord]-2./3.*C_4[ord]-8.*C_5[ord]-32./3.*C_6[ord];
   
    }
    
//     ISIDORI VALUES
//    C_1[LO] = -0.487;
//    C_2[LO] = 1.024;
//    C_7[LO] = -0.321;
//    C_7[NLO] = 0.019;
//    C_8L[LO] = -0.148;
//    C_9[LO] = 4.129;
//    C_9[NLO] = 0.013;
//    C_10[LO] = -4.372;
//    C_10[NLO] = 0.135;
//    T_9[LO] = 0.374;
//    T_9[NLO] = 0.251;
//    U_9[LO] = 0.032;
//    U_9[NLO] = 0.016;
//    W_9[LO] = 0.032;
//    W_9[NLO] = 0.012;

    C_1[FULLNLO] = C_1[LO] + C_1[NLO];
    C_1L_bar[FULLNLO] = C_1L_bar[LO] + C_1L_bar[NLO];
    C_2[FULLNLO] = C_2[LO] + C_2[NLO];
    C_2L_bar[FULLNLO] = C_2L_bar[LO] + C_2L_bar[NLO];
    C_3[FULLNLO] = C_3[LO] + C_3[NLO];
    C_4[FULLNLO] = C_4[LO] + C_4[NLO];
    C_5[FULLNLO] = C_5[LO] + C_5[NLO];
    C_6[FULLNLO] = C_6[LO] + C_6[NLO];
    C_7[FULLNLO] = C_7[LO] + C_7[NLO];
    C_8L[FULLNLO] = C_8L[LO] + C_8L[NLO];
    C_9[FULLNLO] = C_9[LO] + C_9[NLO];
    C_10[FULLNLO] = C_10[LO] + C_10[NLO];
    C_S[FULLNLO] = C_S[LO] + C_S[NLO];
    C_P[FULLNLO] = C_P[LO] + C_P[NLO];
    C_7p[FULLNLO] = C_7p[LO] + C_7p[NLO];
    C_9p[FULLNLO] = C_9p[LO] + C_9p[NLO];
    C_10p[FULLNLO] = C_10p[LO] + C_10p[NLO];
    C_Sp[FULLNLO] = C_Sp[LO] + C_Sp[NLO];
    C_Pp[FULLNLO] = C_Pp[LO] + C_Pp[NLO];
    T_9[FULLNLO] = T_9[LO]+T_9[NLO];
    U_9[FULLNLO] = U_9[LO]+U_9[NLO];
    W_9[FULLNLO] = W_9[LO]+W_9[NLO];    
}

double BXqll::integrate_Rquark(double sh_min, double sh_max, q2regions q2region)
{
    updateParameters();

    old_handler = gsl_set_error_handler_off();
    
    switch(q2region)
    {
        case LOWQ2:
            FR = convertToGslFunction(boost::bind(&BXqll::getR_LOWQ2, &(*this), _1));
    
            if (gsl_integration_cquad(&FR, sh_min, sh_max, 1.e-5, 1.e-4, w_Rquark, &avaRquark, &errRquark, NULL) != 0)
                return std::numeric_limits<double>::quiet_NaN();
            return avaRquark;
            break;
        case HIGHQ2:
            FR = convertToGslFunction(boost::bind(&BXqll::getR_HIGHQ2, &(*this), _1));
    
            if (gsl_integration_cquad(&FR, sh_min, sh_max, 1.e-5, 1.e-4, w_Rquark, &avaRquark, &errRquark, NULL) != 0)
                return std::numeric_limits<double>::quiet_NaN();
            return avaRquark;
            break;
        default:
            throw std::runtime_error("BXqll::integrate_R: region of q^2 not implemented");         
    }
    
    gsl_set_error_handler(old_handler);
}


double BXqll::getR_LOWQ2(double sh)
{
    updateParameters();
//   return (R_quark(sh,LOWQ2)/* + deltaMb2_Rquark(sh,LOWQ2)*/);
//    gslpp::matrix<gslpp::complex> test = matH_L(sh,LO);
//    return test(6,6).real();
//    return (H_T(sh)+H_L(sh));
//    return S77_L(sh,NLO);
   return(C_9[LO].real());
}

double BXqll::getR_HIGHQ2(double sh)
{
    updateParameters();
    return (R_quark(sh,HIGHQ2) + deltaMb2_Rquark(sh,HIGHQ2));
}


double BXqll::deltaMb2_Rquark(double sh, q2regions q2region)
{
    double pre;
    double delta_R;
    
    pre = ale/2./M_PI*Vts_over_Vcb;
    
    delta_R = -4.*(6.+3.*sh-5.*sh*sh*sh)/sh*(C7eff(sh, LO).abs2() + 2.*(C7eff(sh, LO).conjugate()*C7eff(sh, NLO)).real());
    delta_R += (1.-15.*sh*sh+10.*sh*sh*sh)*(C9eff(sh, LO).abs2() + 2.*(C9eff(sh, LO).conjugate()*C9eff(sh, NLO)).real());
    delta_R += (1.-15.*sh*sh+10.*sh*sh*sh)*(C10eff(sh, LO).abs2() + 2.*(C10eff(sh, LO).conjugate()*C10eff(sh, NLO)).real());
    delta_R += -4.*(5.+6.*sh-7.*sh*sh)*(C9eff(sh, LO).conjugate()*C7eff(sh, LO)+
            (C9eff(sh, NLO)+alsmu/M_PI*(omega79(sh)-omega9(sh))*C9eff(sh,LO)).conjugate()*C7eff(sh, LO)+
             C9eff(sh, LO).conjugate()*(C7eff(sh, NLO)+alsmu/M_PI*(omega79(sh)-omega7(sh))*C7eff(sh,LO))).real();
    
    if(q2region <= HIGHQ2)
        return (3.*lambda_2/2./Mb/Mb*(pre*pre/f_sl(z)*delta_R + g_lambda(z)/g(z)*R_quark(sh,q2region)));
    else
        throw std::runtime_error("BXqll::deltaMb2_Rquark: q2 region not implemented");
}

double BXqll::R_quark(double sh, q2regions q2region)
{
    double pre;
    double gamma;
   
    pre = ale/2./M_PI*abslambdat_over_Vcb;
 
    gamma = (1.+2.*sh)*(C9eff(sh, LO).abs2() + 2.*(C9eff(sh, LO).conjugate()*C9eff(sh, NLO)).real());
    gamma += (1.+2.*sh)*(C10eff(sh, LO).abs2() + 2.*(C10eff(sh, LO).conjugate()*C10eff(sh, NLO)).real());
    gamma += 4.*(1.+2./sh)*(C7eff(sh, LO).abs2() + 2.*(C7eff(sh, LO).conjugate()*C7eff(sh, NLO)).real());
    // omega_7,9 replaced with omega_79 for the interaction term
    gamma += 12.*(C9eff(sh, LO).conjugate()*C7eff(sh, LO)+
            (C9eff(sh, NLO)+(alsmu/M_PI*(omega79(sh)-omega9(sh))*C9eff(sh,LO))).conjugate()*C7eff(sh, LO)+
             C9eff(sh, LO).conjugate()*(C7eff(sh, NLO)+(alsmu/M_PI*(omega79(sh)-omega7(sh))*C7eff(sh,LO)))).real();

//    Without finite bremsstrahlung corrections
//    return (pre*pre/f_sl*(1.-sh)*(1.-sh)*gamma);
    
    if(q2region <= HIGHQ2)
        return(pre*pre/f_sl(z)*((1.-sh)*(1.-sh)*gamma /*+ alsmu/4./M_PI*R_bremsstrahlung(sh,q2region)*/));
    else
        throw std::runtime_error("BXqll::R_quark: q2 region not implemented");
}

double BXqll::f_sl(double z)
{
    return (g(z)*(1.-2.*alsmu/3./M_PI/g(z)* // according to asatrian g is computed with pole masses
            (-(1.-z*z)*(25./4.-239./3.*z+25./4.*z*z)+z*log(z)*(20.+90.*z-4./3.*z*z+17./3.*z*z*z)+
            z*z*pow(log(z),2.)*(36.+z*z)+(1.-z*z)*(17./3.-64./3.*z+17./3.*z*z)*log(1.-z)-
            4.*(1.+30.*z*z+z*z*z*z)*log(z)*log(1.-z)-(1.+16.*z*z+z*z*z*z)*(6.*dilog((gslpp::complex) z).real()-M_PI*M_PI)-
            32*pow(z,1.5)*(1.+z)*(M_PI*M_PI-4.*dilog((gslpp::complex) sqrt(z)).real()+4.*dilog((gslpp::complex) (-sqrt(z))).real()-2.*
            log(z)*log((1-sqrt(z))/(1+sqrt(z)))))));
}

double BXqll::g(double z)
{
    return (1.-8.*z+8.*z*z*z-z*z*z*z-12.*z*z*log(z));
}

double BXqll::g_lambda(double z)
{
    return (3.-8.*z+24.*z*z-24.*z*z*z+5.*z*z*z*z+12.*z*z*log(z));
}

double BXqll::omega7(double sh)
{
    double umsh = 1.-sh;

    return (-8./3.*log(muh)-4./3.*dilog((gslpp::complex) sh).real()-2./9.*M_PI*M_PI-2./3.*log(sh)*
            log(umsh)-1./3.*(8.+sh)/(2.+sh)*log(umsh)-2./3.*sh*(2.-2.*sh-sh*sh)/
            umsh/umsh/(2.+sh)*log(sh)-1./18.*(16.-11.*sh-17.*sh*sh)/(2.+sh)/umsh);
}

double BXqll::omega9(double sh)
{
    double umsh = 1.-sh;

    return (-4./3.*dilog((gslpp::complex) sh).real()-2./9.*M_PI*M_PI-2./3.*log(umsh)*log(sh)-1./3.*
            (5.+4.*sh)/(1.+2.*sh)*log(umsh)-2./3.*sh*(1.+sh)*(1.-2.*sh)/umsh/umsh/
            (1.+2.*sh)*log(sh)+1./6.*(5.+9.*sh-6.*sh*sh)/umsh/(1.+2.*sh));
}

double BXqll::omega79(double sh)
{
    double umsh = 1.-sh;

    return (-4./3.*log(muh)-4./3.*dilog((gslpp::complex) sh).real()-2./9.*M_PI*M_PI-2./3.*log(umsh)*log(sh)-
            1./9.*(2.+7.*sh)/sh*log(umsh)-2./9.*sh*(3.-2.*sh)/umsh/umsh*log(sh)+
            1./18.*(5.-9.*sh)/umsh);
}

gslpp::complex BXqll::h_z(double zed, double sh)
{
    gslpp::complex i = gslpp::complex::i();
    gslpp::complex h_z;
    
    if (zed == 0.)
    {
        h_z = 8./27.-4./9.*(log(sh)-i*M_PI);
    }
    else
    {
        h_z = -2./9.*(2.+4./sh*zed)*sqrt(std::abs(4.*zed-sh)/sh);
        if(sh > 4.*zed)
            h_z *= log((sqrt(sh)+sqrt(sh-4.*zed))/(sqrt(sh)-sqrt(sh-4.*zed)))-i*M_PI;
        else
            h_z *= 2.*atan(sqrt(sh/(4.*zed-sh)));
        
        h_z += -4./9.*log(zed)+8./27.+16./9.*zed/sh;
   }   
     
    return(h_z);
}

gslpp::complex BXqll::ATUW(double sh, orders order)
{
        if(order == FULLNLO || order <= NLO)
        return(C_9[order]+T_9[order]*h_z(z,sh)+U_9[order]*h_z(1.,sh)+W_9[order]*h_z(0.,sh));
    else
        throw std::runtime_error("BXqll::ATUW: order not implemented");
}

gslpp::complex BXqll::C7eff(double sh, orders order)
{
    gslpp::complex C7eff;
    
    switch(order)
    {
        case LO:
            C7eff = C_7[LO];
            break;
        case NLO:
            C7eff = C_7[NLO]+alsmu/M_PI*omega7(sh)*C_7[LO]-alsmu/4./M_PI*(C_1[LO]*F17(sh)+
            C_2[LO]*F27(sh)+C_8L[LO]*F87(sh));
            break;
        default:
            throw std::runtime_error("BXqll::C7eff: order not implemented"); 
    }
    return (C7eff);
}

gslpp::complex BXqll::C9eff(double sh, orders order)
{
    gslpp::complex C9eff;
    
    switch(order)
    {
        case LO:
            C9eff = ATUW(sh,LO);
            break;
        case NLO:
            C9eff = ATUW(sh,NLO)+alsmu/M_PI*omega9(sh)*ATUW(sh,LO)-alsmu/4./M_PI*(C_1[LO]*F19(sh)+
            C_2[LO]*F29(sh)+C_8L[LO]*F89(sh));
            break;
        default:
            throw std::runtime_error("BXqll::C9eff: order not implemented");   
    }
    return (C9eff);
}

gslpp::complex BXqll::C10eff(double sh, orders order)
{
    gslpp::complex C10eff;
    
    switch(order)
    {
        case LO:
            C10eff = C_10[LO];
            break;
        case NLO:
            C10eff = C_10[NLO]+alsmu/M_PI*omega9(sh)*C_10[LO];
            break;
        default:
            throw std::runtime_error("BXqll::C10eff: order not implemented");         
    }
    return (C10eff);
}

gslpp::complex BXqll::F19(double sh) 
{
    gslpp::complex i = gslpp::complex::i();
    
    return (F_19re(muh, z, sh) + i*F_19im(muh, z, sh));
}

gslpp::complex BXqll::F29(double sh) 
{
    gslpp::complex i = gslpp::complex::i();
    
    return (F_29re(muh, z, sh) + i*F_29im(muh, z, sh));
}

gslpp::complex BXqll::F17(double sh) 
{
    gslpp::complex i = gslpp::complex::i();
    
    return (F_17re(muh, z, sh) + i*F_17im(muh, z, sh));
}

gslpp::complex BXqll::F27(double sh) 
{
    gslpp::complex i = gslpp::complex::i();
    
    return (F_27re(muh, z, sh) + i*F_27im(muh, z, sh));
}

gslpp::complex BXqll::F87(double sh) 
{
    gslpp::complex i = gslpp::complex::i();
    double ash = asin(sqrt(sh)/2.);
    double umsh = 1.-sh;

    return (4.*M_PI*M_PI/27.*(2.+sh)/umsh/umsh/umsh/umsh-4./9.*(11.-16.*sh+8.*sh*sh)/umsh/umsh-
            8./9.*sqrt(sh*(4.-sh))/umsh/umsh/umsh*(9.-5.*sh+2.*sh*sh)*ash-16./3.*(2+sh)/
            umsh/umsh/umsh/umsh*ash*ash-8./9.*sh/umsh*log(sh)-32./9.*log(muh)-i*8./9.*M_PI);
}

double BXqll::F89(double sh) 
{
    double ash = asin(sqrt(sh)/2.);
    double umsh = 1.-sh;

    return (-8.*M_PI*M_PI/27.*(4.-sh)/umsh/umsh/umsh/umsh+8./9.*(5.-2.*sh)/umsh/umsh+
            16./9.*sqrt((4.-sh)/sh)/umsh/umsh/umsh*(4.+3.*sh-sh*sh)*ash+32./3.*(4.-sh)/
            umsh/umsh/umsh/umsh*ash*ash+16./9./umsh*log(sh));
}

double BXqll::R_bremsstrahlung(double sh, q2regions q2region)
{
    gslpp::complex c78, c89, c17, c27, c18, c28, c19, c29;
    double c88, c11, c12, c22;
    double Brem_a;
    double Brem_b;
    double ctau1 = 1./27.;
    double ctau2 = -2./9.;
   
    c78 = CF*C7eff(sh, LO)*C_8L[LO].conjugate();
    c89 = CF*C_8L[LO]*C9eff(sh, LO).conjugate();
    c88 = CF*C_8L[LO].abs2();
    
    c11 = ctau1*C_1[LO].abs2();
    c12 = ctau2*2.*(C_1[LO]*C_2[LO].conjugate()).real();
    c22 = CF*C_2[LO].abs2();
    
    c17 = ctau2*C_1[LO]*C7eff(sh, LO).conjugate();
    c27 = CF*C_2[LO]*C7eff(sh, LO).conjugate();
    c18 = ctau2*C_1[LO]*C_8L[LO].conjugate();
    c28 = CF*C_2[LO]*C_8L[LO].conjugate();
    c19 = ctau2*C_1[LO]*C9eff(sh, LO).conjugate();
    c29 = CF*C_2[LO]*C9eff(sh, LO).conjugate();
    
    Brem_a = 2.*(c78*tau78(sh)+c89*tau89(sh)).real() + c88*tau88(sh);
    
    if(q2region <= HIGHQ2)
        Brem_b = (c11 + c12 + c22)*tau22fit(sh,q2region)+2.*(c17 + c27).real()*tau27fit_Re(sh,q2region)-
            2.*(c17 + c27).imag()*tau27fit_Im(sh,q2region)+2.*(c18 + c28).real()*tau28fit_Re(sh,q2region)-
            2.*(c18 + c28).imag()*tau28fit_Im(sh,q2region)+2.*(c19 + c29).real()*tau29fit_Re(sh,q2region)-
            2.*(c19 + c29).imag()*tau29fit_Im(sh,q2region);
    else
        throw std::runtime_error("BXqll::R_bremsstrahlung: q2 region not implemented");
    
    return (Brem_a + Brem_b);
}

double BXqll::tau78(double sh)
{
    gslpp::complex i = gslpp::complex::i();
    double dmsh = 2.-sh;
    double qmsh = 4.-sh;
    
    return (8./9./sh*(25.-2.*M_PI*M_PI-27.*sh+3.*sh*sh-sh*sh*sh+12.*(sh+sh*sh)*log(sh)+
            6.*(M_PI/2.-atan((2.-4.*sh+sh*sh)/dmsh/sqrt(sh)/sqrt(qmsh)))*
            (M_PI/2.-atan((2.-4.*sh+sh*sh)/dmsh/sqrt(sh)/sqrt(qmsh)))-
            24.*(dilog((sh-i*sqrt(sh)*sqrt(qmsh))/2.)).real()-12.*((1.-sh)*sqrt(sh)*sqrt(qmsh)-
            2.*atan(sqrt(sh)*sqrt(qmsh)/dmsh))*(atan(sqrt(qmsh/sh))-atan(sqrt(sh)*sqrt(qmsh)/dmsh))));
}

double BXqll::tau89(double sh)
{
    gslpp::complex i = gslpp::complex::i();
    double dmsh = 2.-sh;
    double qmsh = 4.-sh;
    
    return (2./3.*(sh*qmsh-3.-4.*log(sh)*(1.-sh-sh*sh)-8.*(dilog(sh/2.+i*sqrt(sh)*sqrt(qmsh)/2.)-
            dilog((sh*qmsh-2.)/2.+i*dmsh*sqrt(sh)*sqrt(qmsh)/2.)).real()+
            4.*(sh*sh*sqrt(qmsh/sh)+2.*atan(sqrt(sh)*sqrt(qmsh)/dmsh))*(atan(sqrt(qmsh/sh))-
            atan(sqrt(sh)*sqrt(qmsh)/dmsh))));
}

double BXqll::tau88(double sh)
{   
    gslpp::complex i = gslpp::complex::i();
    double umsh = 1.-sh;
    double qmsh = 4.-sh;
    
    return (4./27./sh*(-8.*M_PI*M_PI+umsh*(77.-sh-4.*sh*sh)-24.*dilog((gslpp::complex) umsh).real()+
            3.*(10.-4.*sh-9.*sh*sh+8.*log(sqrt(sh)/umsh))*log(sh)+48.*(dilog((3.-sh)/2.+
            i*umsh*sqrt(qmsh)/2./sqrt(sh))).real()-6.*((20.*sh+10.*sh*sh-3.*sh*sh*sh)/sqrt(sh)/sqrt(qmsh)
            -8.*M_PI+8.*atan(sqrt(qmsh/sh)))*(atan(sqrt(qmsh/sh))-atan(sqrt(sh)*sqrt(qmsh)/(2.-sh)))));
}

double BXqll::tau22fit(double sh, q2regions q2region)
{
    double fit;
    
    switch(q2region)
    {
        case LOWQ2:
            fit = -186.96738127 + 1313.45792139*sh - 8975.40399683*sh*sh + 47018.56440838*sh*sh*sh -
                    159603.3217871*sh*sh*sh*sh + 309228.13379963*sh*sh*sh*sh*sh - 258317.14719949*sh*sh*sh*sh*sh*sh -
                    51.2467544*log(sh);
            break;
        case HIGHQ2:
            fit = -322.73989723 + 4.75813656/sh/sh - 80.36414222/sh + 687.70415138*sh - 491.08241967*sh*sh +
                    303.28125994*sh*sh*sh - 132.82124268*sh*sh*sh*sh + 35.63127394*sh*sh*sh*sh*sh -
                    4.36712003*sh*sh*sh*sh*sh*sh - 306.899641*log(sh);
            break;
        default:
            throw std::runtime_error("BXqll::tau22fit: region of q^2 not implemented");
    }
    
    return (fit);
}

double BXqll::tau27fit_Re(double sh, q2regions q2region)
{
    double fit;
    
    switch(q2region)
    {
        case LOWQ2:
            fit = -45.40905903+334.92509492*sh-2404.69181358*sh*sh+12847.93973401*sh*sh*sh-
                    44421.35127703*sh*sh*sh*sh+87786.54536182*sh*sh*sh*sh*sh-75574.96266083*sh*sh*sh*sh*sh*sh-
                    13.79251644*log(sh);
            break;
        case HIGHQ2:
            fit = 87.43391175-196.67646862*sh+219.51106756*sh*sh-184.44868587*sh*sh*sh+
                    103.59892754*sh*sh*sh*sh-34.56056777*sh*sh*sh*sh*sh+5.14181565*sh*sh*sh*sh*sh*sh+
                    38.55667004*log(sh);
            break;
        default:
            throw std::runtime_error("BXqll::tau27fit_Re: region of q^2 not implemented");
    }
    
    return (fit);
}

double BXqll::tau27fit_Im(double sh, q2regions q2region)
{
    double fit;
    
    switch(q2region)
    {
        case LOWQ2:
            fit = -189.61083508+1349.85607262*sh-9198.62227938*sh*sh+48104.40980548*sh*sh*sh-
                    162998.75872037*sh*sh*sh*sh+315224.375522*sh*sh*sh*sh*sh-262649.64320483*sh*sh*sh*sh*sh*sh-
                    52.52183304*log(sh);
            break;
        case HIGHQ2:
            fit = 523.76263422+49.97156504/sh-1120.42920341*sh+1024.46949612*sh*sh-767.28958612*sh*sh*sh+
                    393.62561539*sh*sh*sh*sh-120.74162898*sh*sh*sh*sh*sh+16.63110789*sh*sh*sh*sh*sh*sh+
                    352.74960196*log(sh);
            break;
        default:
            throw std::runtime_error("BXqll::tau27fit_Im: region of q^2 not implemented");
    }
    
    return (fit);
}

double BXqll::tau28fit_Re(double sh, q2regions q2region)
{
    double fit;
    
    switch(q2region)
    {
        case LOWQ2:
            fit = 8.67757227-85.91172547*sh+666.57779497*sh*sh-3619.65234448*sh*sh*sh+
                    12475.74303361*sh*sh*sh*sh-24365.45545631*sh*sh*sh*sh*sh+20446.33269814*sh*sh*sh*sh*sh*sh+
                    1.54278226*log(sh);
            break;
        case HIGHQ2:
            fit = -4.11234905-0.52681762/sh+8.21844628*sh-6.04601107*sh*sh+3.67099354*sh*sh*sh-
                    1.57120958*sh*sh*sh*sh+0.41975346*sh*sh*sh*sh*sh-0.05280596*sh*sh*sh*sh*sh*sh-
                    3.16331567*log(sh);
            break;
        default:
            throw std::runtime_error("BXqll::tau28fit_Re: region of q^2 not implemented");
    }
    
    return (fit);
}

double BXqll::tau28fit_Im(double sh, q2regions q2region)
{
    double fit;
    
    switch(q2region)
    {
        case LOWQ2:
            fit = 57.88258299-430.77957254*sh+3002.9999511*sh*sh-15808.63980887*sh*sh*sh+
                    53787.08410769*sh*sh*sh*sh-104360.60205475*sh*sh*sh*sh*sh+87294.84251167*sh*sh*sh*sh*sh*sh+
                    14.61062129*log(sh);
            break;
        case HIGHQ2:
            fit = -24.92802842+0.3842418/sh/sh-6.38294139/sh+53.15600599*sh-37.59024636*sh*sh+
                    23.04316804*sh*sh*sh-10.03556518*sh*sh*sh*sh+2.68088049*sh*sh*sh*sh*sh-
                    0.32751495*sh*sh*sh*sh*sh*sh-24.01652729*log(sh);
            break;
        default:
            throw std::runtime_error("BXqll::tau28fit_Im: region of q^2 not implemented");
    }
    
    return (fit);
}

double BXqll::tau29fit_Re(double sh, q2regions q2region)
{
    double fit;
    
    switch(q2region)
    {
        case LOWQ2:
            fit = 0.53834924+0.47775224*sh-16.20063604*sh*sh+101.36668267*sh*sh*sh-
                    466.94537092*sh*sh*sh*sh+1224.77742613*sh*sh*sh*sh*sh-1469.41817323*sh*sh*sh*sh*sh*sh-
                    0.01686348*log(sh);
            break; 
        case HIGHQ2:
            fit = 4.46985355-6.16130742*sh+0.84917331*sh*sh+1.7696124*sh*sh*sh-1.14453916*sh*sh*sh*sh+
                    0.24261178*sh*sh*sh*sh*sh-0.02540446*sh*sh*sh*sh*sh*sh+2.67164817*log(sh);
            break;
        default:
            throw std::runtime_error("BXqll::tau29fit_Re: region of q^2 not implemented");
    }
    
    return (fit);
}

double BXqll::tau29fit_Im(double sh, q2regions q2region)
{
    double fit;
    
    switch(q2region)
    {
        case LOWQ2:
            fit = 0.7688748-0.21680402*sh-1.16934757*sh*sh+8.31833871*sh*sh*sh-4.81289468*sh*sh*sh*sh-
                    51.53765482*sh*sh*sh*sh*sh+158.06040784*sh*sh*sh*sh*sh*sh-0.00485643*log(sh);
            break;
        case HIGHQ2:
            fit = -38.80905455+95.60697233*sh-124.04368889*sh*sh+118.64599185*sh*sh*sh-
                    73.76081228*sh*sh*sh*sh+26.55080999*sh*sh*sh*sh*sh-4.19021877*sh*sh*sh*sh*sh*sh-
                    16.02711369*log(sh);
            break;
        default:
            throw std::runtime_error("BXqll::tau29fit_Im: region of q^2 not implemented");
    }
    
    return (fit);
}

double BXqll::F_17re(double muh, double z, double sh, int maxpow)
{
    return myF_1.F_17re(muh, z, sh, maxpow);
};

double BXqll::F_17im(double muh, double z, double sh, int maxpow)
{
    return myF_1.F_17im(muh, z, sh, maxpow);
};

double BXqll::F_19re(double muh, double z, double sh, int maxpow)
{
    return myF_1.F_19re(muh, z, sh, maxpow);
};

double BXqll::F_19im(double muh, double z, double sh, int maxpow)
{
    return myF_1.F_19im(muh, z, sh, maxpow);
};

double BXqll::F_27re(double muh, double z, double sh, int maxpow)
{
    return myF_2.F_27re(muh, z, sh, maxpow);
};

double BXqll::F_27im(double muh, double z, double sh, int maxpow)
{
    return myF_2.F_27im(muh, z, sh, maxpow);
};

double BXqll::F_29re(double muh, double z, double sh, int maxpow)
{
    return myF_2.F_29re(muh, z, sh, maxpow);
};

double BXqll::F_29im(double muh, double z, double sh, int maxpow)
{
    return myF_2.F_29im(muh, z, sh, maxpow);
};

double BXqll::DeltaF_19re(double muh, double z, double sh, int maxpow)
{
    return myF_1.DeltaF_19re(muh, z, sh, maxpow);
};

double BXqll::DeltaF_19im(double muh, double z, double sh, int maxpow)
{
    return myF_1.DeltaF_19im(muh, z, sh, maxpow);
};

double BXqll::DeltaF_29re(double muh, double z, double sh, int maxpow)
{
    return myF_2.DeltaF_29re(muh, z, sh, maxpow);
};

double BXqll::DeltaF_29im(double muh, double z, double sh, int maxpow)
{
    return myF_2.DeltaF_29im(muh, z, sh, maxpow);
};

/*
 * Implementation of the notation of @cite Huber:2015sra
 */

double BXqll::H_T(double sh)
{
    gslpp::matrix<gslpp::complex> HT_LO = matH_T(sh,LO);
    gslpp::matrix<gslpp::complex> HT_NLO = matH_T(sh,NLO);
//    for(int i=0; i<10; i++)
//        for(int j=0; j<10; j++)
//           std::cout << i << "," << j << ": " << (HT_LO(i,j)!=0. ? HT_NLO(i,j)/HT_LO(i,j): 0.) << std::endl; 
    double Phill_T = 0.;
    double pre = 2.*abslambdat_over_Vcb;
    
    for(int j=0; j<10; j++)
    {
        for(int i=0; i<=j; i++)
        {
            Phill_T += ((*(allcoeff[LO]))(i)*(*(allcoeff[LO]))(j).conjugate()*(HT_LO(i,j)/*+HT_NLO(i,j)*/) +
                       ((*(allcoeff[LO]))(i)*(*(allcoeff[NLO]))(j).conjugate()+
                       (*(allcoeff[NLO]))(i)*(*(allcoeff[LO]))(j).conjugate())*HT_LO(i,j)).real();
        }
    }

    return pre*pre*Phill_T/f_sl(z);
}

double BXqll::H_L(double sh)
{
    gslpp::matrix<gslpp::complex> HL_LO = matH_L(sh,LO);
    gslpp::matrix<gslpp::complex> HL_NLO = matH_L(sh,NLO);
    double Phill_L = 0.;
    double pre = 2.*abslambdat_over_Vcb;
    
    for(int j=0; j<10; j++)
    {
        for(int i=0; i<=j; i++)
        {
            Phill_L += ((*(allcoeff[LO]))(i)*((*(allcoeff[LO]))(j)).conjugate()*(HL_LO(i,j)/*+HL_NLO(i,j)*/) +
                       ((*(allcoeff[LO]))(i)*(*(allcoeff[NLO]))(j).conjugate()+
                       (*(allcoeff[NLO]))(i)*(*(allcoeff[LO]))(j).conjugate())*HL_LO(i,j)).real();
        }
    }

    return pre*pre*Phill_L/f_sl(z);
}

gslpp::matrix<gslpp::complex> BXqll::matH_T(double sh, orders order)
{
    gslpp::matrix<gslpp::complex> Hij_T(10,10,0.);
    
    gslpp::vector<gslpp::complex> M7_LO = vecM7(sh,LO);
    gslpp::vector<gslpp::complex> M9_LO = vecM9(sh,LO);
    gslpp::vector<double>         M10_LO = vecM10(sh);
    gslpp::vector<gslpp::complex> M7_NLO = vecM7(sh,NLO);
    gslpp::vector<gslpp::complex> M9_NLO = vecM9(sh,NLO);

    switch(order)
    {
        case LO:
            for(int j=0; j<10; j++)
            {
                for(int i=0; i<=j; i++)
                {
                    if(i==j)
                        Hij_T.assign(i,j, S77_T(sh,LO)*M7_LO(i)*M7_LO(i).conjugate() +
                                          S99_T(sh,LO)*M9_LO(i)*M9_LO(i).conjugate() +
                                          S1010_T(sh,LO)*M10_LO(i)*M10_LO(i) /*+
                                          S79_T(sh,LO)*(M7_LO(i)*M9_LO(i).conjugate()).real()*/);
            
                    else
                        Hij_T.assign(i,j, 2.*S77_T(sh,LO)*M7_LO(i)*M7_LO(j).conjugate() +
                                          2.*S99_T(sh,LO)*M9_LO(i)*M9_LO(j).conjugate() +
                                          2.*S1010_T(sh,LO)*M10_LO(i)*M10_LO(j) /*+
                                          S79_T(sh,LO)*(M7_LO(i)*M9_LO(j).conjugate()+M9_LO(i)*M7_LO(j).conjugate())*/);
                        
                }
            }
            break;
        case NLO:
            for(int j=0; j<10; j++)
            {
                for(int i=0; i<=j; i++)
                {
                    if(i==j)
                        Hij_T.assign(i,j, S77_T(sh,NLO)*M7_LO(i)*M7_LO(i).conjugate()+
                        2*S77_T(sh,LO)*(M7_LO(i)*M7_NLO(i).conjugate()).real() +
                                          S99_T(sh,NLO)*M9_LO(i)*M9_LO(i).conjugate()+
                        2*S99_T(sh,LO)*(M9_LO(i)*M9_NLO(i).conjugate()).real() +
                                          S1010_T(sh,NLO)*M10_LO(i)*M10_LO(i) /*+
                                          S79_T(sh,NLO)*(M7_LO(i)*M9_LO(i).conjugate()).real()+
                        S79_T(sh,LO)*(M7_LO(i)*M9_NLO(i).conjugate()+M7_NLO(i)*M9_LO(i).conjugate()).real()*/);
            
                    else
                        Hij_T.assign(i,j, 2.*(S77_T(sh,NLO)*M7_LO(i)*M7_LO(j).conjugate()+
                        S77_T(sh,LO)*(M7_LO(i)*M7_NLO(j).conjugate() + M7_NLO(i)*M7_LO(j).conjugate())) +
                                          2.*(S99_T(sh,NLO)*M9_LO(i)*M9_LO(j).conjugate()+
                        S99_T(sh,LO)*(M9_LO(i)*M9_NLO(j).conjugate() + M9_NLO(i)*M9_LO(j).conjugate())) +
                                          2.*S1010_T(sh,NLO)*M10_LO(i)*M10_LO(j) /*+
                                          S79_T(sh,NLO)*(M7_LO(i)*M9_LO(j).conjugate()+M9_LO(i)*M7_LO(j).conjugate())+
                        S79_T(sh,LO)*(M7_LO(i)*M9_NLO(j).conjugate() + M7_NLO(i)*M9_LO(j).conjugate()+
                        M9_LO(i)*M7_NLO(j).conjugate() + M9_NLO(i)*M7_LO(j).conjugate())*/);
                }
            }
            break;
        default:
            throw std::runtime_error("BXqll::matH_T: order not implemented");
    }
    return Hij_T;
}

gslpp::matrix<gslpp::complex> BXqll::matH_L(double sh, orders order)
{
    gslpp::matrix<gslpp::complex> Hij_L(10,10,0.);
    
    gslpp::vector<gslpp::complex> M7_LO = vecM7(sh,LO);
    gslpp::vector<gslpp::complex> M7_NLO = vecM7(sh,NLO);
    gslpp::vector<gslpp::complex> M9_LO = vecM9(sh,LO);
    gslpp::vector<gslpp::complex> M9_NLO = vecM9(sh,NLO);
    gslpp::vector<double>         M10_LO = vecM10(sh);
    
    switch(order)
    {
        case LO:
            for(int j=0; j<10; j++)
            {
                for(int i=0; i<=j; i++)
                {
                    if(i==j)
                        Hij_L.assign(i,j, S77_L(sh,LO)*M7_LO(i)*M7_LO(i).conjugate() +
                                          S99_L(sh,LO)*M9_LO(i)*M9_LO(i).conjugate() +
                                          S1010_L(sh,LO)*M10_LO(i)*M10_LO(i) +
                                          S79_L(sh,LO)*(M7_LO(i)*M9_LO(i).conjugate()).real());
            
                    else
                        Hij_L.assign(i,j, 2.*S77_L(sh,LO)*M7_LO(i)*M7_LO(j).conjugate() +
                                          2.*S99_L(sh,LO)*M9_LO(i)*M9_LO(j).conjugate() +
                                          2.*S1010_L(sh,LO)*M10_LO(i)*M10_LO(j) +
                                          S79_L(sh,LO)*(M7_LO(i)*M9_LO(j).conjugate()+M9_LO(i)*M7_LO(j).conjugate()));
                        
                }
            }
            break;
        case NLO:
            for(int j=0; j<10; j++)
            {
                for(int i=0; i<=j; i++)
                {
                    if(i==j)
                        Hij_L.assign(i,j, S77_L(sh,NLO)*M7_LO(i)*M7_LO(i).conjugate()+
                        2*S77_L(sh,LO)*(M7_LO(i)*M7_NLO(i).conjugate()).real() +
                                          S99_L(sh,NLO)*M9_LO(i)*M9_LO(i).conjugate()+
                        2*S99_L(sh,LO)*(M9_LO(i)*M9_NLO(i).conjugate()).real() +
                                          S1010_L(sh,NLO)*M10_LO(i)*M10_LO(i) +
                                          S79_L(sh,NLO)*(M7_LO(i)*M9_LO(i).conjugate()).real()+
                        S79_L(sh,LO)*(M7_LO(i)*M9_NLO(i).conjugate()+M7_NLO(i)*M9_LO(i).conjugate()).real());
            
                    else
                        Hij_L.assign(i,j, 2.*(S77_L(sh,NLO)*M7_LO(i)*M7_LO(j).conjugate()+
                        S77_L(sh,LO)*(M7_LO(i)*M7_NLO(j).conjugate() + M7_NLO(i)*M7_LO(j).conjugate())) +
                                          2.*(S99_L(sh,NLO)*M9_LO(i)*M9_LO(j).conjugate()+
                        S99_L(sh,LO)*(M9_LO(i)*M9_NLO(j).conjugate() + M9_NLO(i)*M9_LO(j).conjugate())) +
                                          2.*S1010_L(sh,NLO)*M10_LO(i)*M10_LO(j) +
                                          S79_L(sh,NLO)*(M7_LO(i)*M9_LO(j).conjugate()+M9_LO(i)*M7_LO(j).conjugate())+
                        S79_L(sh,LO)*(M7_LO(i)*M9_NLO(j).conjugate() + M7_NLO(i)*M9_LO(j).conjugate()+
                        M9_LO(i)*M7_NLO(j).conjugate() + M9_NLO(i)*M7_LO(j).conjugate()));
                }
            }
            break;
        default:
            throw std::runtime_error("BXqll::matH_L: order not implemented");
    }
    
    return Hij_L;
}

gslpp::vector<gslpp::complex> BXqll::vecM7(double sh, orders order)
{
    gslpp::vector<gslpp::complex> M7(10,0.);

    switch(order)
    {
        case LO:
            M7.assign(6, ale/4./M_PI);
            break;
        case NLO:
            M7.assign(0, -ale*alsmu/16./M_PI/M_PI*F17(sh));
            M7.assign(1, -ale*alsmu/16./M_PI/M_PI*F27(sh));
            M7.assign(7, -ale*alsmu/16./M_PI/M_PI*F87(sh));
            break;
        default:
            throw std::runtime_error("BXqll::vecM7: order not implemented");
    }
    
    return M7;
}

gslpp::vector<gslpp::complex> BXqll::vecM9(double sh, orders order)
{
    gslpp::vector<gslpp::complex> M9(10,0.);
    
    switch(order)
    {
        case LO:
            M9.assign(0, ale/4./M_PI*f_Huber(sh, -32./27., 4./3., 0.,      0.,       -16./27.));
            M9.assign(1, ale/4./M_PI*f_Huber(sh, -8./9.,   1.,    0.,      0.,       -4./9.));
            M9.assign(2, ale/4./M_PI*f_Huber(sh, -16./9.,  6.,    -7./2.,  2./9.,    2./27.));
            M9.assign(3, ale/4./M_PI*f_Huber(sh, 32./27.,  0.,    -2./3.,  8./27.,   8./81.));
            M9.assign(4, ale/4./M_PI*f_Huber(sh, -112./9., 60.,   -38.,    32./9.,   -136./27.));
            M9.assign(5, ale/4./M_PI*f_Huber(sh, 512./27., 0.,    -32./3., 128./27., 320./81.));
//            M9.assign(8, 1. + ale/4./M_PI*f9pen_Huber(sh));
            break;
        case NLO:
            M9.assign(0, -ale*alsmu/16./M_PI/M_PI*F19(sh));
            M9.assign(1, -ale*alsmu/16./M_PI/M_PI*F19(sh));
            M9.assign(7, -ale*alsmu/16./M_PI/M_PI*F89(sh));
            break;
        default:
            throw std::runtime_error("BXqll::vecM9: order not implemented");
    }
    
    return M9;
}

gslpp::vector<double> BXqll::vecM10(double sh)
{
    gslpp::vector<double> M10(10,0.);

//    M10(9) = 1.;
    
    return M10;
}

double BXqll::S77_T(double sh, orders order)
{
    double sigma = 8.*(1.-sh)*(1.-sh)/sh;
    
    switch(order)
    {
        case LO:
            return sigma;
            break;
        case NLO:
            return sigma*8.*alsmu*omega77_T(sh)/4./M_PI;
            break;
        default:
            throw std::runtime_error("BXqll::S77_T: order not implemented");
    }
}

double BXqll::S79_T(double sh, orders order)
{
    double sigma = 8.*(1.-sh)*(1.-sh);
    
    switch(order)
    {
        case LO:
            return sigma;
            break;
        case NLO:
            return sigma*8.*alsmu*omega79_T(sh)/4./M_PI;
            break;
        default:
            throw std::runtime_error("BXqll::S79_T: order not implemented");
    }
}

double BXqll::S99_T(double sh, orders order)
{
    double sigma = 2.*sh*(1.-sh)*(1.-sh);
    
    switch(order)
    {
        case LO:
            return sigma;
            break;
        case NLO:
            return sigma*8.*alsmu*omega99_T(sh)/4./M_PI;
            break;
        default:
            throw std::runtime_error("BXqll::S99_T: order not implemented");
    }
}

double BXqll::S1010_T(double sh, orders order)
{
    return S99_T(sh,order);
}

double BXqll::S77_L(double sh, orders order)
{
    double sigma = 4.*(1.-sh)*(1.-sh);
    
    switch(order)
    {
        case LO:
            return sigma;
            break;
        case NLO:
            return sigma*8.*alsmu*omega77_L(sh)/4./M_PI;
            break;
        default:
            throw std::runtime_error("BXqll::S77_L: order not implemented");
    }
}

double BXqll::S79_L(double sh, orders order)
{
    double sigma = 4.*(1.-sh)*(1.-sh);
    
    switch(order)
    {
        case LO:
            return sigma;
            break;
        case NLO:
            return sigma*8.*alsmu*omega79_L(sh)/4./M_PI;
            break;
        default:
            throw std::runtime_error("BXqll::S79_L: order not implemented");
    }
}

double BXqll::S99_L(double sh, orders order)
{
    double sigma = (1.-sh)*(1.-sh);
    
    switch(order)
    {
        case LO:
            return sigma;
            break;
        case NLO:
            return sigma*8.*alsmu*omega99_L(sh)/4./M_PI;
            break;
        default:
            throw std::runtime_error("BXqll::S99_L: order not implemented");
    }
}

double BXqll::S1010_L(double sh, orders order)
{
    return S99_L(sh,order);
}

double BXqll::S710_A(double sh, orders order)
{
    double sigma = -8.*(1.-sh)*(1.-sh);
    
    switch(order)
    {
        case LO:
            return sigma;
            break;
        case NLO:
            return sigma*8.*alsmu*omega710_A(sh)/4./M_PI;
            break;
        default:
            throw std::runtime_error("BXqll::S710_A: order not implemented");
    }
}

double BXqll::S910_A(double sh, orders order)
{
    double sigma = -4.*(1.-sh)*(1.-sh);
    
    switch(order)
    {
        case LO:
            return sigma;
            break;
        case NLO:
            return sigma*8.*alsmu*omega910_A(sh)/4./M_PI;
            break;
        default:
            throw std::runtime_error("BXqll::S910_A: order not implemented");
    }
}


double BXqll::omega77_T(double sh)
{
    double umsh = 1.-sh;
    double umsqrt = 1.-sqrt(sh);
    
    return (-8./3.*log(muh) - (sqrt(sh)+1.)*(sqrt(sh)+1.)*(pow(sh,1.5)-10.*sh+13.*sqrt(sh)-8.)*
            dilog((gslpp::complex) umsh).real()/6./umsh/umsh + 2.*sqrt(sh)*(sh*sh-6.*sh-3.)*
            dilog((gslpp::complex) umsqrt).real()/3./umsh/umsh - M_PI*M_PI*(3.*pow(sh,1.5)+22.*sh+23.*sqrt(sh)+16.)*
            umsqrt*umsqrt/36./umsh/umsh + (5.*sh*sh*sh-54.*sh*sh+57.*sh-8.)/18./umsh/umsh - log(umsh) + sh*(5.*sh+
            1.)*log(sh)/3./umsh/umsh + 2./3.*log(umsh)*log(sh));
}

double BXqll::omega79_T(double sh)
{
    double umsh = 1.-sh;
    double umsqrt = 1.-sqrt(sh);
    
    return (-4./3.*log(muh) - 2.*sqrt(sh)*(sh+3.)*dilog((gslpp::complex) umsqrt).real()/3./umsh/umsh - M_PI*M_PI*
            (16.*sh+29.*sqrt(sh)+19.)*umsqrt*umsqrt/36./umsh/umsh + (sh*sh-6.*sh+5.)/6./umsh/umsh + (sqrt(sh)+1.)*
            (sqrt(sh)+1.)*(8.*sh-15.*sqrt(sh)+9.)*dilog((gslpp::complex) umsh).real()/6./umsh/umsh - (5.*sh+1.)*
            log(umsh)/6./sh + sh*(3.*sh+1.)*log(sh)/6./umsh/umsh + 2./3.*log(umsh)*log(sh));
}

double BXqll::omega99_T(double sh)
{
    double umsh = 1.-sh;
    double umsqrt = 1.-sqrt(sh);
    
    return ((sqrt(sh)+1.)*(sqrt(sh)+1.)*(8.*pow(sh,1.5)-15.*sh+4.*sqrt(sh)-5.)*dilog((gslpp::complex) umsh).real()/
            6./umsh/umsh/sqrt(sh) - 2.*(sh*sh-12.*sh-5.)*dilog((gslpp::complex) umsqrt).real()/3./umsh/umsh/sqrt(sh) -
            M_PI*M_PI*(16.*pow(sh,1.5)+29.*sh+4.*sqrt(sh)+15.)*umsqrt*umsqrt/36./umsh/umsh/sqrt(sh) + (2.*sh*sh-7.*sh-
            5.)*log(sh)/3./umsh/umsh + (sh*sh+18.*sh-19.)/6./umsh/umsh - (2.*sh+1)*log(umsh)/3./sh + 2./3.*log(umsh)*
            log(sh));
}

double BXqll::omega77_L(double sh)
{
    double umsh = 1.-sh;
    double umsqrt = 1.-sqrt(sh);
    
    return (-8./3.*log(muh) + (sqrt(sh)+1.)*(sqrt(sh)+1.)*(4.*pow(sh,1.5)-7.*sh+2.*sqrt(sh)-3.)*
            dilog((gslpp::complex) umsh).real()/3./umsh/umsh/sqrt(sh) - (9.*sh*sh-38.*sh+29.)/6./umsh/umsh -
            4.*(sh*sh-6.*sh-3.)*dilog((gslpp::complex) umsqrt).real()/3./umsh/umsh/sqrt(sh) - M_PI*M_PI*
            (8.*pow(sh,1.5)+13.*sh+2.*sqrt(sh)+9.)*umsqrt*umsqrt/18./umsh/umsh/sqrt(sh) - (sh*sh*sh-3.*sh+2.)*
            log(umsh)/3./umsh/umsh/sh + 2.*(sh*sh-3.*sh-3.)*log(sh)/3./umsh/umsh + 2./3.*log(umsh)*log(sh));
}

double BXqll::omega79_L(double sh)
{
    double umsh = 1.-sh;
    double umsqrt = 1.-sqrt(sh);
    
    return (-4./3.*log(muh) + 4.*sqrt(sh)*(sh+3.)*dilog((gslpp::complex) umsqrt).real()/3./umsh/umsh + (sqrt(sh)+1.)*
            (sqrt(sh)+1.)*(4.*sh-9.*sqrt(sh)+3.)*dilog((gslpp::complex) umsh).real()/3./umsh/umsh + (7.*sh*sh-2.*sh-5.)/
            6./umsh/umsh - M_PI*M_PI*(8.*sh+19.*sqrt(sh)+5.)*umsqrt*umsqrt/18./umsh/umsh - (2.*sh+1.)*log(umsh)/3./sh +
            (sh-7.)*sh*log(sh)/3./umsh/umsh + 2./3.*log(umsh)*log(sh));
}

double BXqll::omega99_L(double sh)
{
    double umsh = 1.-sh;
    double umsqrt = 1.-sqrt(sh);
    
    return (-(sqrt(sh)+1.)*(sqrt(sh)+1.)*(pow(sh,1.5)-8.*sh+3.*sqrt(sh)-4.)*dilog((gslpp::complex) umsh).real()/
            3./umsh/umsh + 4.*sqrt(sh)*(sh*sh-12.*sh-5.)*dilog((gslpp::complex) umsqrt).real()/3./umsh/umsh -
            M_PI*M_PI*(3.*pow(sh,1.5)+20.*sh+sqrt(sh)+8.)*umsqrt*umsqrt/18./umsh/umsh + (4.*sh*sh*sh-51.*sh*sh+42*sh+5.)/
            6./umsh/umsh - log(umsh) + 8.*sh*(2.*sh+1.)*log(sh)/3./umsh/umsh + + 2./3.*log(umsh)*log(sh));
}

double BXqll::omega710_A(double sh)
{
    double umsh = 1.-sh;
    double umsqrt = 1.-sqrt(sh);
    
    return (-4./3.*log(muh) + 2.*(4.*sh*sh-13.*sh-1.)*dilog((gslpp::complex) umsqrt).real()/3./umsh/umsh - (2.*sh*sh-
            9.*sh-3.)*dilog((gslpp::complex) umsh).real()/3./umsh/umsh - (3.*sh*sh-16.*sh+13.)*log(umsqrt)/3./umsh/umsh +
            (4.*sh*sh-13.*sh-1.)*log(umsqrt)*log(sh)/3./umsh/umsh - (2.*sh*sh-9.*sh-3.)*log(umsh)*log(sh)/3./umsh/umsh +
            (sh*sh*sh-23.*sh*sh+23.*sh-1.)*log(umsh)/6./umsh/umsh/sh + (sh-20.*sqrt(sh)+5.)*umsqrt*umsqrt/6./umsh/umsh -
            M_PI*M_PI/3.);
}

double BXqll::omega910_A(double sh)
{
    double umsh = 1.-sh;
    double umsqrt = 1.-sqrt(sh);
    
    return (-2.*(sh*sh-3.*sh-1.)*dilog((gslpp::complex) umsh).real()/3./umsh/umsh - 4.*(5.-2.*sh)*sh*
            dilog((gslpp::complex) umsqrt).real()/3./umsh/umsh - (4.*sqrt(sh)-3.)*umsqrt*umsqrt/3./umsh/umsh -
            2.*(2.*sh*sh-7.*sh+5.)*log(umsqrt)/3./umsh/umsh - 2.*(sh*sh-3.*sh-1.)*log(umsh)*log(sh)/3./umsh/umsh +
            (2.*sh*sh*sh-11.*sh*sh+10.*sh-1.)*log(umsh)/3./umsh/umsh/sh + 2.*sh*(2.*sh-5.)*log(umsqrt)*log(sh)/
            3./umsh/umsh - M_PI*M_PI/3.);
}

gslpp::complex BXqll::f_Huber(double sh, double gamma_9, double rho_c, double rho_b, double rho_0, double rho_num)
{
    gslpp::complex i = gslpp::complex::i();
    
    return (-gamma_9*log(muh) + rho_c*(g_Huber(4.*z/sh)+8./9.*log(Mb/Mc)) + rho_b*g_Huber(4.*1./sh) +
            rho_0*(log(sh)-i*M_PI) + rho_num);
}

gslpp::complex BXqll::f9pen_Huber(double sh)
{
    gslpp::complex i = gslpp::complex::i();
    
    return (-8.*log(muh) - 3.*g_Huber(4.*Mtau*Mtau/Mb/Mb/sh) - 8./3.*log(Mb/Mtau) + 8./3.*(log(sh)-i*M_PI) - 40./9.);
}

gslpp::complex BXqll::g_Huber(double y)
{
    gslpp::complex i = gslpp::complex::i();
    gslpp::complex g_y;
    
    g_y = -2./9.*(2.+y)*sqrt(abs(1.-y));
        if(y < 1.)
            g_y *= log(abs((1.+sqrt(1.-y))/(1.-sqrt(1.-y))))-i*M_PI;
        else
            g_y *= 2.*atan(1./sqrt(y-1.));
        
        g_y += 20./27. + 4./9.*y;
    
    return (g_y);
}