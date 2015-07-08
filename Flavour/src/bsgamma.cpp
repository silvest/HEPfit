/* 
 * Copyright (C) 2015 SusyFit Collaboration
 * All rights reserved.
 *
 * For the licensing terms see doc/COPYING.
 */

#include "bsgamma.h"
#include <gslpp.h>
#include <gsl/gsl_sf_dilog.h>
#include <gsl/gsl_sf_zeta.h>
#include <boost/bind.hpp>

Bsgamma::Bsgamma(const StandardModel& SM_i, int obsFlag)
: ThObservable(SM_i)
{
    if (obsFlag > 0 and obsFlag < 2) obs = obsFlag;
    else throw std::runtime_error("obsFlag in bsgamma can only be 1 (BR)");
    
    w_Phi221_1 = gsl_integration_workspace_alloc (1000);
    w_Phi221_2 = gsl_integration_workspace_alloc (1000);
    w_Phi221_3 = gsl_integration_workspace_alloc (1000);
    w_Phi221_4 = gsl_integration_workspace_alloc (1000);
    w_Phi221_5 = gsl_integration_workspace_alloc (1000);
    
    w_Phi271_1 = gsl_integration_workspace_alloc (1000);
    w_Phi271_2 = gsl_integration_workspace_alloc (1000);
    w_Phi271_3 = gsl_integration_workspace_alloc (1000);
}

double Bsgamma::delta(double E0)
{
    return 1. - 2.*E0/Mb1s;
}

double Bsgamma::zeta()
{
    return Mc*Mc/Mb1s/Mb1s;
}

gslpp::complex Bsgamma::a(double z)
{
    double zeta3 = gsl_sf_zeta_int(3);
    
    double z2=z*z;
    double z3=z2*z;
    double z4=z3*z;
    double z5=z4*z;
    double z6=z5*z;
    
    double L=log(z);
    double L2=L*L;
    double L3=L2*L;
    
    double pi2=M_PI*M_PI;
    
    if (z == 1.) return 4.0859 + 4./9. * M_PI * gslpp::complex::i();
    else return 16./9. * ( ( 5./2. - pi2/3. - 3.*zeta3 + ( 5./2. - 3./4.*pi2 )*L + L2/4. + L3/12. )*z 
            + ( 7./4. + 2./3.*pi2 - pi2*L/2. - L2/4. + L3/12. )*z2 + ( -7./6. - pi2/4. + 2*L - 3./4.*L2 )*z3 
            + ( 457./216. - 5./18*pi2 - L/72. - 5./6.*L2 )*z4 + ( 35101./8640. - 35./72.*pi2 - 185./144.*L - 35./24.*L2 )*z5
            + ( 67801./8000. - 21./20.*pi2 - 3303./800.*L - 63./20.*L2 )*z6 + 
            gslpp::complex::i()*M_PI*( ( 2. - pi2/6. + L/2. + L2/2. )*z + ( 1./2. - pi2/6. - L + L2/2. )*z2
            + z3 + 5./9.*z4 + 49./72.*z5 + 231./200.*z6) );
}

gslpp::complex Bsgamma::b(double z)
{
    
    double z2=z*z;
    double z3=z2*z;
    double z4=z3*z;
    double z5=z4*z;
    double z6=z5*z;
    
    double L=log(z);
    double L2=L*L;
    
    double pi2=M_PI*M_PI;
    
    if (z == 1.) return 0.0316 + 4./81. * M_PI * gslpp::complex::i();
    else return -8./9. * ( ( -3. + pi2/6. - L )*z - 2./3.*pi2*pow(z,3./2.) + ( 1./2. + pi2 -2.*L - L2/2. )*z2 
            + ( -25./12. - pi2/9. - 19./18.*L + 2.*L2 )*z3 + ( -1376./225. + 137./30.*L + 2.*L2 + 2./3.*pi2 )*z4 
            + ( -131317./11760. + 887./84.*L + 5.*L2 + 5./3.*pi2 )*z5 
            + ( -2807617./97200. + 16597./540.*L + 14.*L2 + 14./3.*pi2 )*z6 + 
            gslpp::complex::i()*M_PI*( -z + ( 1 - 2.*L )*z2 + ( -10./9. + 4./3.*L )*z3 + z4 + 2./3.*z5 + 7./9.*z6) );
}

gslpp::complex Bsgamma::Gamma_t(double t)
{
    if (t<4) return -2. * atan( sqrt(t/(4.-t)) ) * atan( sqrt(t/(4.-t)) );
    else return -M_PI*M_PI/2. + 2.*log( ( sqrt(t) + sqrt(t-4.) ) / 2. )*log( ( sqrt(t) + sqrt(t-4.) ) / 2. )
            - 2.*gslpp::complex::i()*M_PI*log( ( sqrt(t) + sqrt(t-4.) ) / 2. );
}

gslpp::complex Bsgamma::r(int i, double z)
{
    
    double Xb = -0.16844083981858157;
    
    switch(i){
        case 1:
            return 833./729. - (a(z) + b(z))/3. + 40./243.*gslpp::complex::i()*M_PI;
        case 2:
            return - 1666./243. + 2.*(a(z) + b(z)) - 80./81.*gslpp::complex::i()*M_PI;
        case 3:
            return 2392./243. + 8.*M_PI/3./sqrt(3.) + 32./9.*Xb - a(1.) + 2.*b(1.) + 56./81.*gslpp::complex::i()*M_PI;
        case 4:
            return -761./729. - 4.*M_PI/9./sqrt(3.) - 16./27.*Xb + a(1.)/6. + 5.*b(1.)/3. + 2.*b(z) - 148./243.*gslpp::complex::i()*M_PI;
        case 5:
            return 56680./243. + 32.*M_PI/3./sqrt(3.) + 128./9.*Xb - 16.*a(1.) + 32.*b(1.) + 896./81.*gslpp::complex::i()*M_PI;
        case 6:
            return 5710./729. - 16.*M_PI/9./sqrt(3.) - 64./27.*Xb - 10./3.*a(1.) + 44./3.*b(1.) + 12.*a(z) + 20.*b(z) 
                    - 2296./243.*gslpp::complex::i()*M_PI;
        default:
            std::stringstream out;
            out << i;
            throw std::runtime_error("Bsgamma::r(): index " + out.str() + " not implemented");
    }
}

double Bsgamma::Phi22_1(double E0)
{
    
    double t1 = (1. - delta(E0))/zeta();
    double t2 = 1./zeta();
    
    FPhi221_1 = convertToGslFunction( boost::bind( &Bsgamma::getPhi221, &(*this), _1 ) );
    gsl_integration_qags (&FPhi221_1, 0, t1, 1.e-5, 1.e-3, 1000, w_Phi221_1, &avaPhi221_1, &errPhi221_1);

    FPhi221_2 = convertToGslFunction( boost::bind( &Bsgamma::getPhi221_t, &(*this), _1 ) );
    gsl_integration_qags (&FPhi221_2, 0, t1, 1.e-5, 1.e-3, 1000, w_Phi221_2, &avaPhi221_2, &errPhi221_2);
    
    FPhi221_3 = convertToGslFunction( boost::bind( &Bsgamma::getPhi221, &(*this), _1 ) );
    gsl_integration_qags (&FPhi221_3, t1, t2, 1.e-5, 1.e-3, 1000, w_Phi221_3, &avaPhi221_3, &errPhi221_3);
    
    FPhi221_4 = convertToGslFunction( boost::bind( &Bsgamma::getPhi221_t, &(*this), _1 ) );
    gsl_integration_qags (&FPhi221_4, t1, t2, 1.e-5, 1.e-3, 1000, w_Phi221_4, &avaPhi221_4, &errPhi221_4);
    
    FPhi221_5 = convertToGslFunction( boost::bind( &Bsgamma::getPhi221_t2, &(*this), _1 ) );
    gsl_integration_qags (&FPhi221_5, t1, t2, 1.e-5, 1.e-3, 1000, w_Phi221_5, &avaPhi221_5, &errPhi221_5);
    
    
    return 16./27. * zeta() * (delta(E0)*(avaPhi221_1 - zeta()*avaPhi221_2) 
            + avaPhi221_3 - 2.*zeta()*avaPhi221_4 + zeta()*zeta()*avaPhi221_5);
}

double Bsgamma::Phi11_1(double E0)
{
    return Phi22_1(E0)/36.;
}

double Bsgamma::Phi12_1(double E0)
{
    return -Phi22_1(E0)/3.;
}

double Bsgamma::Phi27_1(double E0)
{
    
    double t1 = (1. - delta(E0))/zeta();
    double t2 = 1./zeta();
    
    FPhi271_1 = convertToGslFunction( boost::bind( &Bsgamma::getPhi271, &(*this), _1 ) );
    gsl_integration_qags (&FPhi271_1, 0, t1, 1.e-5, 1.e-3, 1000, w_Phi271_1, &avaPhi271_1, &errPhi271_1);

    FPhi271_2 = convertToGslFunction( boost::bind( &Bsgamma::getPhi271, &(*this), _1 ) );
    gsl_integration_qags (&FPhi271_2, t1, t2, 1.e-5, 1.e-3, 1000, w_Phi271_2, &avaPhi271_2, &errPhi271_2);

    FPhi271_3 = convertToGslFunction( boost::bind( &Bsgamma::getPhi271_t, &(*this), _1 ) );
    gsl_integration_qags (&FPhi271_3, t1, t2, 1.e-5, 1.e-3, 1000, w_Phi271_3, &avaPhi271_3, &errPhi271_3);
    
    
    return -8./9. * zeta() * zeta() * (delta(E0)*avaPhi271_1 + avaPhi271_2 - zeta()*avaPhi271_3);
}

double Bsgamma::Phi17_1(double E0)
{
    return -Phi27_1(E0)/6.;
}

double Bsgamma::Phi18_1(double E0)
{
    return Phi27_1(E0)/18.;
}

double Bsgamma::Phi28_1(double E0)
{
    return -Phi27_1(E0)/3.;
}

double Bsgamma::Phi47_1(double E0)
{
    double d=delta(E0);
    double d2=d*d;
    double d3=d2*d;
    
    return 1./54.*M_PI*( 3*sqrt(3) - M_PI ) + 1./81.*d3 - 25./108.*d2 + 5./54.*d 
            + 2./9.*( d2 +2.*d + 3. )*pow(atan(sqrt( (1. - d) / (3. + d) )),2)
            - 1./3.*( d2 +4.*d + 3. )*sqrt( (1. - d) / (3. + d) )*atan(sqrt( (1. - d) / (3. + d) ));
}

double Bsgamma::Phi77_1(double E0)
{
    double d=delta(E0);
    double d2=d*d;
    double d3=d2*d;
    
    return -2./3.*pow(log(d),2.) - 7./3.*log(d) - 31./9. + 10./3.*d + d2/3. - 2./9.*d3 + d*(d - 4.)*log(d)/3.;
}

double Bsgamma::Phi78_1(double E0)
{
    double d=delta(E0);
    double d2=d*d;
    double d3=d2*d;
    
    double Li2 = gsl_sf_dilog(1. - d);
    
    double pi2=M_PI*M_PI;
    
    return 8./9.*( Li2 -  pi2/6. - d*log(d) + 9./4.*d - d2/4. + d3/12.);
}

double Bsgamma::Phi88_1(double E0)
{
    double d=delta(E0);
    double d2=d*d;
    double d3=d2*d;
    
    double Li2 = gsl_sf_dilog(1. - d);
    
    double pi2=M_PI*M_PI;
    
    return 1./27.*( -2.*log(Mb1s/Ms)*( d2 + 2.*d + 4.*log(1. - d) ) + 4.*Li2 - 2./3.*pi2 - d*(2. + d)*log(d) 
            + 8.*log(1. - d) -  2./3.*d3 + 3.*d2 +7*d);
}

double Bsgamma::Kij_1(int i, int j, double E0, double mu)
{
    if (i > j) throw std::runtime_error("Bsgamma::Kij_1(): index i must not be greater than index j");
    
    double gamma_i7[8] = {-208./243., 416./81., -176./81., -152./243., -6272./81., 4624./243., 32./3., -32./9.};
    double Phi_i7[6] = {Phi17_1(E0), Phi27_1(E0), 0., Phi47_1(E0), 0., 0.};
    
    if (i<7 && j==7) return r(i,zeta()).real() - gamma_i7[i-1]*log(mu/Mb1s) + 2.*Phi_i7[i-1];
    
    else if (i==1 && j==1) return 4.*Phi11_1(E0);
    else if (i==1 && j==2) return 2.*Phi12_1(E0);
    else if (i==1 && j==8) return 2.*Phi18_1(E0);
    
    else if (i==2 && j==2) return 4.*Phi22_1(E0);
    else if (i==2 && j==8) return 2.*Phi28_1(E0);
    
    else if (i==4 && j==8) return -2./3.*Phi47_1(E0);
    
    else if (i==7 && j==7) return -182./9. + 8./9.*M_PI*M_PI - gamma_i7[6]*2.*log(mu/Mb1s) + 4.*Phi77_1(E0);
    else if (i==7 && j==8) return 44./9. - 8./27.*M_PI*M_PI - gamma_i7[7]*log(mu/Mb1s) + 2.*Phi78_1(E0);
    
    else if (i==8 && j==8) return 4.*Phi88_1(E0);
    
    else throw std::runtime_error("Bsgamma::Kij_1(): indexes (i,j) not implemented");
}

void Bsgamma::computeCoeff(double mu)
{
    allcoeff = SM.getMyFlavour()->ComputeCoeffsgamma(mu);
    
    C1_0 = (*(allcoeff[LO]))(0);
    C2_0 = (*(allcoeff[LO]))(1);
    C3_0 = (*(allcoeff[LO]))(2);
    C4_0 = (*(allcoeff[LO]))(3);
    C5_0 = (*(allcoeff[LO]))(4);
    C6_0 = (*(allcoeff[LO]))(5);
    C7_0 = (*(allcoeff[LO]))(6);
    C8_0 = (*(allcoeff[LO]))(7);
    
    C1_1 = 4.*M_PI/SM.Als(mu)*(*(allcoeff[NLO]))(0);
    C2_1 = 4.*M_PI/SM.Als(mu)*(*(allcoeff[NLO]))(1);
    C3_1 = 4.*M_PI/SM.Als(mu)*(*(allcoeff[NLO]))(2);
    C4_1 = 4.*M_PI/SM.Als(mu)*(*(allcoeff[NLO]))(3);
    C5_1 = 4.*M_PI/SM.Als(mu)*(*(allcoeff[NLO]))(4);
    C6_1 = 4.*M_PI/SM.Als(mu)*(*(allcoeff[NLO]))(5);
    C7_1 = 4.*M_PI/SM.Als(mu)*(*(allcoeff[NLO]))(6);
    C8_1 = 4.*M_PI/SM.Als(mu)*(*(allcoeff[NLO]))(7);
    
    /*std::cout << "C1_0: " << C1_0 << std::endl;
    std::cout << "C2_0: " << C2_0 << std::endl;
    std::cout << "C3_0: " << C3_0 << std::endl;
    std::cout << "C4_0: " << C4_0 << std::endl;
    std::cout << "C5_0: " << C5_0 << std::endl;
    std::cout << "C6_0: " << C6_0 << std::endl;
    std::cout << "C7_0: " << C7_0 << std::endl;
    std::cout << "C8_0: " << C8_0 << std::endl;
    
    std::cout << "C1_1: " << C1_1 << std::endl;
    std::cout << "C2_1: " << C2_1 << std::endl;
    std::cout << "C3_1: " << C3_1 << std::endl;
    std::cout << "C4_1: " << C4_1 << std::endl;
    std::cout << "C5_1: " << C5_1 << std::endl;
    std::cout << "C6_1: " << C6_1 << std::endl;
    std::cout << "C7_1: " << C7_1 << std::endl;
    std::cout << "C8_1: " << C8_1 << std::endl;*/
    
    /*C1_0 = -0.8411;
    C2_0 = 1.0647;
    C3_0 = -0.0133;
    C4_0 = -0.1276;
    C5_0 = 0.0012;
    C6_0 = 0.0028;
    C7_0 = -0.3736;
    C8_0 = -0.1729;
    
    C1_1 = 15.278;
    C2_1 = -2.124;
    C3_1 = 0.096;
    C4_1 = -0.463;
    C5_1 = -0.021;
    C6_1 = -0.013;
    C7_1 = 2.027;
    C8_1 = -0.617;*/

}

double Bsgamma::P21(double E0, double mu)
{
    return (C1_0*C1_0).real() * Kij_1(1,1,E0,mu) + 2.*(C1_0*C2_0).real() * Kij_1(1,2,E0,mu) 
            + 2.*(C1_0*C7_0).real() * Kij_1(1,7,E0,mu) + 2.*(C1_0*C8_0).real() * Kij_1(1,8,E0,mu)
            + (C2_0*C2_0).real() * Kij_1(2,2,E0,mu) + 2.*(C2_0*C7_0).real() * Kij_1(2,7,E0,mu) 
            + 2.*(C2_0*C8_0).real() * Kij_1(2,8,E0,mu) + 2.*(C3_0*C7_0).real() * Kij_1(3,7,E0,mu)
            + 2.*(C4_0*C7_0).real() * Kij_1(4,7,E0,mu) + 2.*(C4_0*C8_0).real() * Kij_1(4,8,E0,mu)
            + 2.*(C5_0*C7_0).real() * Kij_1(5,7,E0,mu) + 2.*(C6_0*C7_0).real() * Kij_1(6,7,E0,mu) 
            + (C7_0*C7_0).real() * Kij_1(7,7,E0,mu) + 2.*(C7_0*C8_0).real() * Kij_1(7,8,E0,mu) 
            + (C8_0*C8_0).real() * Kij_1(8,8,E0,mu);
}

double Bsgamma::P32(double E0, double mu)
{
    return 2.*( (C1_0*C1_1).real() * Kij_1(1,1,E0,mu) + (C1_0*C2_1 + C1_1*C2_0).real() * Kij_1(1,2,E0,mu) 
            + (C1_0*C7_1 + C1_1*C7_0).real() * Kij_1(1,7,E0,mu) + (C1_0*C8_1 + C1_1*C8_0).real() * Kij_1(1,8,E0,mu)
            + (C2_0*C2_1).real() * Kij_1(2,2,E0,mu) + (C2_0*C7_1 + C2_1*C7_0).real() * Kij_1(2,7,E0,mu) 
            + (C2_0*C8_1 + C2_1*C8_0).real() * Kij_1(2,8,E0,mu) + (C3_0*C7_1 + C3_1*C7_0).real() * Kij_1(3,7,E0,mu)
            + (C4_0*C7_1 + C4_1*C7_0).real() * Kij_1(4,7,E0,mu) + (C4_0*C8_1 + C4_1*C8_0).real() * Kij_1(4,8,E0,mu) 
            + (C5_0*C7_1 + C5_1*C7_0).real() * Kij_1(5,7,E0,mu) + (C6_0*C7_1 + C6_1*C7_0).real() * Kij_1(6,7,E0,mu) 
            + (C7_0*C7_1).real() * Kij_1(7,7,E0,mu) + (C7_0*C8_1 + C7_1*C8_0).real() * Kij_1(7,8,E0,mu) 
            + (C8_0*C8_1).real() * Kij_1(8,8,E0,mu) );
}

double Bsgamma::P(double E0, double mu, orders order)
{
    switch(order) {
        case NLO:
//            std::cout << "p0: " << C7_0.abs2() << std::endl;
//            std::cout << "p11: " << 2.*(C7_0*C7_1).real() << std::endl;
//            std::cout << "p21: " << P21(E0,mu) << std::endl;
//            std::cout << "p32: " << P32(E0,mu) << std::endl;
            return C7_0.abs2() + SM.Als(mu,FULLNLO)/4./M_PI * (2.*(C7_0*C7_1).real() + P21(E0,mu));
            break;
        case LO:
            return C7_0.abs2();
            break;
        default:
            std::stringstream out;
            out << order;
            throw std::runtime_error("Bsgamma::P(): order " + out.str() + " not implemented");
    }
}

void Bsgamma::computeBR(orders order)
{
    ale=SM.getAle();
    E0=SM.getbsgamma_E0();
    Mb1s=4.68;
    mu_b = Mb1s/2.;
    Mc=SM.getQuarks(QCD::CHARM).getMass();
    Ms=SM.getQuarks(QCD::STRANGE).getMass();
    BRsl=SM.getBr_B_Xcenu();
    C=SM.getbsgamma_C();
    lambda_t=SM.computelamt_s();
    V_cb=SM.getCKM().getVcb();
    V_cb=0.04174304221;
    
    computeCoeff(mu_b);
    
    BR = BRsl * (lambda_t/V_cb).abs2() * 6. * ale / (M_PI * C) * P(E0, mu_b, order);
}

double Bsgamma::computeThValue()
{
    computeBR(NLO);
    
    if (obs == 1) return BR;
    
    throw std::runtime_error("Bsgamma::computeThValue(): Observable type not defined. Can be only 1");
}