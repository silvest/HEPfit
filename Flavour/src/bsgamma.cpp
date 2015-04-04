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

Bsgamma::Bsgamma(const StandardModel& SM_i, int obsFlag): ThObservable(SM_i), mySM(SM_i){
    if (obsFlag > 0 and obsFlag < 3) obs = obsFlag;
    else throw std::runtime_error("obsFlag in bsgamma can only be 1 (BR) or 2 (A_{CP})");
    
    w_Phi221_1 = gsl_integration_workspace_alloc (1000);
    w_Phi221_2 = gsl_integration_workspace_alloc (1000);
    w_Phi221_3 = gsl_integration_workspace_alloc (1000);
    w_Phi221_4 = gsl_integration_workspace_alloc (1000);
    w_Phi221_5 = gsl_integration_workspace_alloc (1000);
    
    w_Phi271_1 = gsl_integration_workspace_alloc (1000);
    w_Phi271_2 = gsl_integration_workspace_alloc (1000);
    w_Phi271_3 = gsl_integration_workspace_alloc (1000);
}

/*
double Bsgamma::H1(double E0, double Mu){
    double z=2*E0/Mb;
    double z2=z*z;
    double z3=z2*z;
    double L=log(Mu/Mb);
    
    return -4./9. * (-1. + z) * (16. - 12.*L - 4.*M_PI*M_PI - 30.*z - 3.*z2 + 2.*z3 + 3.*(-3. + 2.*z + z2)*log(1. - z));
}

double Bsgamma::f(double rho){
    double rho2 = rho*rho;
    double Li2 = gsl_sf_dilog(rho);
    double Li2sqrt = gsl_sf_dilog(sqrt(rho));
    Polylogarithms myPolylog;
    
    return - M_PI*M_PI/9.*(162.*sqrt(rho) + 70.*pow(rho,3./2.) - 36.*rho2) + 32./9.*rho*log(rho) + 2./9.*(25. + 27.*rho2)*log(rho)*log(rho)
            + 4./9.*log(rho)*log(rho)*log(rho) - 4./9.*(31. + 27.*rho2)*log(1. - rho)*log(rho)
            - 4./9.*sqrt(rho)*(81. + 35.*rho)*atanh(sqrt(rho))*log(rho) - 2./9.*(62. + 81.*sqrt(rho) + 35.*pow(rho,3./2.) + 54.*rho2)*Li2
            + 8./3.*log(rho)*Li2 + 8./9.*sqrt(rho)*(81. + 35.*rho)*Li2sqrt - 16./3.*myPolylog.Li3(rho) + 5578./81. + 172./9.*rho;
}

double Bsgamma::H2NH(double E0, double Mu){
    double z=2*E0/Mb;
    double z2=z*z;
    double z3=z2*z;
    double L=log(Mu/Mb);
    double zeta3 = gsl_sf_zeta_int(3);
    
    return 4./243. * (-1. + z) * (-3563. - 216.*L*L + 348.*M_PI*M_PI - 36.*L*(-3. + 4.*M_PI*M_PI + 30.*z + 3.*z2 - 2.*z3) 
            + 108.*L*(-3. + 2.*z + z2)*log(1. - z) + 216.*zeta3);
}

double Bsgamma::H2NL(double E0, double Mu){
    double z=2*E0/Mb;
    double z2=z*z;
    double z3=z2*z;
    double L=log(Mu/Mb);
    double lg = log(1. - z);
    double lg2 = lg*log(1. - z);
    double Li2 = gsl_sf_dilog(z);
    Polylogarithms myPolylog;
    double zeta3 = gsl_sf_zeta_int(3);
    
    return -32./1296. * (-1. + z) * (-251 - 72.*L + 144.*L*L + 128.*M_PI*M_PI + 96.*L*M_PI*M_PI + 720*z + 720.*L*z + 246.*z2 
            + 72.*L*z2 - 84.*z3 - 48.*L*z3 + 426.*lg + 216.*L*lg - 24.*M_PI*M_PI*lg - 264.*z*lg - 144.*L*z*lg - 186.*z2*lg
            - 72.*L*z2*lg + 24.*z3*lg - 162.*lg2 + 108.*z*lg2 + 54.*z2*lg2 + 144.*lg2*log(z) + 
            36.*(6. + 2.*z + z2 + 4.*lg)*Li2 + 288.*myPolylog.Li3(1. - z) + 144.*zeta3);
}

double Bsgamma::H2NV(double E0, double Mu){
    double z=2*E0/Mb;
    double z2=z*z;
    double z3=z2*z;
    double rho = Mc*Mc/Mb/Mb;
    double L=log(Mu/Mb);
    
    return 2./81. * (-1. + z) * (-27.*f(rho) + 72.*L - 144.*L*L - 124.*M_PI*M_PI - 96.*L*M_PI*M_PI - 720.*L*z - 72.*L*z2 + 48.*L*z3
            + 72.*L*(-3. + 2.*z + z2)*log(1. - z) + 2.*log(rho)*(-493. + 12.*M_PI*M_PI + 180.*z + 18.*z2 - 12.*z3 
            - 18.*(-3. + 2.*z + z2)*log(1. - z)));
}

double Bsgamma::H2(double E0, double Mu){
    return  NH*H2NH(E0,Mu) + NL*H2NL(E0,Mu) + NV*H2NV(E0,Mu);// + CF*H2a(E0,Mu) + CA*H2na(E0,Mu);
}

double Bsgamma::G77(orders order, double E0, double Mu){
    switch(order){
        case FULLNNLO:
            return 1. + SM.Als(Mu, FULLNNLO)/(4. * M_PI) * H1(E0, Mu) + pow(SM.Als(Mu, FULLNNLO)/(4. * M_PI), 2.) * H2(E0, Mu);
        case FULLNLO:
            return 1. + SM.Als(Mu, FULLNLO)/(4. * M_PI) * H1(E0, Mu);
        case LO:
            return 1.;
        default:
            std::stringstream out;
            out << order;
            throw std::runtime_error("Bsgamma::G77(): order " + out.str() + " not implemented");
    }
}

double Bsgamma::Phi1(double rho){
    double y = (1. - sqrt(1. - 4*rho)) / (1. + sqrt(1. - 4*rho));
    
    if (rho < 0.25)
        return log(y)*log(y) - M_PI*M_PI;
    else return -acos(1. - 1./2./rho) * acos(1. - 1./2./rho);
}

double Bsgamma::Phi2(double rho){
    double y = (1. - sqrt(1. - 4*rho)) / (1. + sqrt(1. - 4*rho));
    
    if (rho < 0.25)
        return sqrt(1. - 4*rho)*log(y);
    else return -sqrt(4*rho - 1.) * acos(1. - 1./2./rho);
}

double Bsgamma::Phi3(double rho){
    double y = (1. - sqrt(1. - 4*rho)) / (1. + sqrt(1. - 4*rho));
    double Li2 = gsl_sf_dilog(-y);
    ClausenFunctions myCl;
    
    if (rho < 0.25)
        return sqrt(1. - 4*rho) * (Li2 + log(y)*log(y)/4. + M_PI*M_PI/12.);
    else return -sqrt(4*rho - 1.) * myCl.Cl2( 2.*asin( 1./2./sqrt(rho) ) );
}

double Bsgamma::Phi4(double rho){
    double y = (1. - sqrt(1. - 4*rho)) / (1. + sqrt(1. - 4*rho));
    Polylogarithms myPolylog;
    ClausenFunctions myCl;
    
    if (rho < 0.25)
        return myPolylog.Li3(-y) + log(y)*log(y)*log(y)/12. + M_PI*M_PI*log(y)/12.;
    else return myCl.Cl3( 2.*asin( 1./2./sqrt(rho) ) );
}

double Bsgamma::dY1(double E0){
    double z0=2*E0/Mb;
    double Li2 = gsl_sf_dilog(z0);
    
    return 2./9.*z0*(z0*z0 + 24.) - 8./3.*(z0 - 1.)*log(1. - z0) - 8./3.*Li2;
}

double Bsgamma::Y2CF(double E0, double Mu){
    double z=2*E0/Mb;
    double L=log(Mu/Mb);
    double Li2 = gsl_sf_dilog(z);
    
    double z2=z*z;
    double z3=z2*z;
    double z4=z3*z;
    double z5=z4*z;
    double z6=z5*z;
    double z7=z6*z;
    double z8=z7*z;
    double z9=z8*z;
    double z10=z9*z;
    double z11=z10*z;
    double z12=z11*z;
    
    return - 35.9375 - 98.4349*L - 106.6667*L*L + 22.1701*z + 45.1015*L*z + 106.6667*L*L*z + 82.3522*z2 + 48*L*z2 
            - 70.3394*z3 + 8.8889*L*z3 - 2.2770*z4 - 3.5556*L*z4 + 10.6755*z5 - 18.2907*z6 + 27.5167*z7 
            - 30.6703*z8 + 23.6558*z9 - 11.8480*z10 + 3.4262*z11 - 0.4333*z12 + 
            ( -25.1395 + 78.9646*z - 82.5106*z2 + 28.6855*z3 + L*(-16 + 26.6667*z - 5.3333*z2 - 5.3333*z3) ) * log(1. - z) 
            + (-6.3944 + 12.5166*z - 5.8499*z2 - 0.2723*z3) * log(1. - z)*log(1. - z) + 
            (-2.8857 + 8.6572*z - 8.6572*z2 + 2.8858*z3) * log(1. - z)*log(1. - z)*log(1. - z) + (-32 + 32*z) * Li2;
}

double Bsgamma::Y2CA(double E0, double Mu){
    double z=2*E0/Mb;
    double L=log(Mu/Mb);
    double Li2 = gsl_sf_dilog(z);
    
    double z2=z*z;
    double z3=z2*z;
    double z4=z3*z;
    double z5=z4*z;
    double z6=z5*z;
    double z7=z6*z;
    double z8=z7*z;
    double z9=z8*z;
    double z10=z9*z;
    double z11=z10*z;
    double z12=z11*z;
    
    return 22.8959 + 50.0741*L + 90.6667*L*L - 19.9497*z + 10.3704*L*z - 90.6667*L*L*z - 16.1984*z2 - 60.4444*L*z2 
            + 6.2832*z3 + 2.5185*L*z3 + 7.087*z4 - 2.5185*L*z4 - 2.1478*z5 + 4.9339*z6 - 6.4431*z7 + 6.0963*z8 
            - 3.4988*z9 + 0.9589*z10 + 0.0424*z11 - 0.0598*z12 + (2.9462 + 30.2222*L*(1. - z)*(1. - z) - 15.3409*z 
            + 21.8433*z2 - 9.4486*z3) * log(1. - z) + 6.6159 * (z - 1.)*(z - 1.)*(z - 1.) * log(1. - z)*log(1. - z) + 
            (0.69995 - 2.0998*z + 2.0999*z2 - 0.69995*z3) * log(1. - z)*log(1. - z)*log(1. - z) + L * 30.2222 * (z - 1.) * Li2;
}

double Bsgamma::Y2NH(double E0, double Mu){
    double L=log(Mu/Mb);
    double zeta3 = gsl_sf_zeta_int(3);
    ClausenFunctions myCl;
    
    return 8./81. * (244. - 27.*sqrt(3.)*M_PI - 61.*M_PI*M_PI) - 64./27. * (18 - M_PI*M_PI) * L - 64./9. * L*L - 64./27. * zeta3 
            + 32.*sqrt(3)*myCl.Cl2(M_PI/3.) + 8./3. * dY1(E0) * L;
}

double Bsgamma::Y2NL(double E0, double Mu){
    double z=2*E0/Mb;
    double L=log(Mu/Mb);
    double zeta3 = gsl_sf_zeta_int(3);
    double Li2 = gsl_sf_dilog(z);
    Polylogarithms myPolylog;
    
    return -16./81. * (328. - 13.*M_PI*M_PI) - 64./27. * (18 - M_PI*M_PI) * L - 64./9. * L * L + 64./3. * zeta3
            + 4./27. * z * (7*z*z - 17*z + 238) + 8./3. * dY1(E0) * L - 8./27. * (z*z*z - 6.*z*z + 80*z - 75 + 6*M_PI*M_PI) * log(1.-z)
            + 16./3. * (z - 1.) * log(1. - z) * log(1. - z) + 16./3. * log(z) * log(1. - z) * log(1. - z) + 32./27. * (3*z - 8) * Li2
            + 32./3. * log(1. - z) * Li2 - 32./9. * myPolylog.Li3(z) + 32./3. * myPolylog.Li3(1. - z) - 32./3. * zeta3;
}

double Bsgamma::Y2NV(double E0, double Mu){
    double L=log(Mu/Mb);
    double rho = Mc*Mc/Mb/Mb;
    double Li2 = gsl_sf_dilog(1. - rho);
    double Li2sqrt = gsl_sf_dilog(1. - sqrt(rho));
    Polylogarithms myPolylog;
    
    return -16./81. * ( 157. - 279.*rho - M_PI*M_PI*( 5. + 9.*rho*rho - 42.*pow(rho,3./2.) ) ) - 64./27. * (18. - M_PI*M_PI) * L
            -64./9. * L*L + 16./27. * (22. - M_PI*M_PI + 10*rho) * log(rho) + 16./27. * (8. + 9.*rho*rho) * log(rho)*log(rho)
            -16./27. * log(rho)*log(rho)*log(rho) - 8./9. * (1. - 6.*rho*rho)*Phi1(rho)  - 8./27. * (19. - 46.*rho)*Phi2(rho)
            -32./27. * (13. + 14.*rho)*Phi3(rho) - 64./9. * Phi4(rho) - 32./9. * log(rho) * Li2 
            + 32./27. * ( 5. + 9.*rho*rho + 14.*pow(rho, 3./2.) ) * Li2 - 1792./27. * pow(rho, 3./2.) * Li2sqrt
            + 64./9. * myPolylog.Li3(1. - rho) + 64./9. * myPolylog.Li3(1. - 1./rho) + 4./3. * dY1(E0) * (2.*L - log(rho));
}

double Bsgamma::Y1(double E0, double Mu){
    double L=log(Mu/Mb);
    
    return 4./9. * (29. - 2. * M_PI * M_PI) + 16./3.*L - dY1(E0);
}

double Bsgamma::Y2(double E0, double Mu){
    return CF*Y2CF(E0,Mu) + CA*Y2CA(E0,Mu) + TR*NH*Y2NH(E0,Mu) + TR*NL*Y2NL(E0,Mu) + TR*NV*Y2NV(E0,Mu);
}

double Bsgamma::G78(orders order, double E0, double Mu){
    switch(order){
        case FULLNNLO:
            return SM.Als(Mu, FULLNNLO)/(4. * M_PI) * CF * Y1(E0, Mu) + pow(SM.Als(Mu, FULLNNLO)/(4. * M_PI), 2.) * CF * Y2(E0, Mu);
        case FULLNLO:
            return SM.Als(Mu, FULLNLO)/(4. * M_PI) * CF * Y1(E0, Mu);
        case LO:
            return 0.;
        default:
            std::stringstream out;
            out << order;
            throw std::runtime_error("Bsgamma::G78(): order " + out.str() + " not implemented");
    }
}*/

double Bsgamma::delta(double E0){
    return 1. - 2.*E0/Mb1s;
}

double Bsgamma::zeta(){
    return Mc*Mc/Mb1s/Mb1s;
}

gslpp::complex Bsgamma::a(double z){
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

gslpp::complex Bsgamma::b(double z){
    
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

gslpp::complex Bsgamma::Gamma_t(double t){
    if (t<4) return -2. * atan( sqrt(t/(4.-t)) ) * atan( sqrt(t/(4.-t)) );
    else return -M_PI*M_PI/2. + 2.*log( ( sqrt(t) + sqrt(t-4.) ) / 2. )*log( ( sqrt(t) + sqrt(t-4.) ) / 2. )
            - 2.*gslpp::complex::i()*M_PI*log( ( sqrt(t) + sqrt(t-4.) ) / 2. );
}

gslpp::complex Bsgamma::r(int i, double z){
    
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

double Bsgamma::Phi22_1(double E0){
    
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

double Bsgamma::Phi11_1(double E0){
    return Phi22_1(E0)/36.;
}

double Bsgamma::Phi12_1(double E0){
    return -Phi22_1(E0)/3.;
}

double Bsgamma::Phi27_1(double E0){
    
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

double Bsgamma::Phi17_1(double E0){
    return -Phi27_1(E0)/6.;
}

double Bsgamma::Phi18_1(double E0){
    return Phi27_1(E0)/18.;
}

double Bsgamma::Phi28_1(double E0){
    return -Phi27_1(E0)/3.;
}

double Bsgamma::Phi47_1(double E0){
    double d=delta(E0);
    double d2=d*d;
    double d3=d2*d;
    
    return 1./54.*M_PI*( 3*sqrt(3) - M_PI ) + 1./81.*d3 - 25./108.*d2 + 5./54.*d 
            + 2./9.*( d2 +2.*d + 3. )*pow(atan(sqrt( (1. - d) / (3. + d) )),2)
            - 1./3.*( d2 +4.*d + 3. )*sqrt( (1. - d) / (3. + d) )*atan(sqrt( (1. - d) / (3. + d) ));
}

double Bsgamma::Phi77_1(double E0){
    double d=delta(E0);
    double d2=d*d;
    double d3=d2*d;
    
    return -2./3.*pow(log(d),2.) - 7./3.*log(d) - 31./9. + 10./3.*d + d2/3. - 2./9.*d3 + d*(d - 4.)*log(d)/3.;
}

double Bsgamma::Phi78_1(double E0){
    double d=delta(E0);
    double d2=d*d;
    double d3=d2*d;
    
    double Li2 = gsl_sf_dilog(1. - d);
    
    double pi2=M_PI*M_PI;
    
    return 8./9.*( Li2 -  pi2/6. - d*log(d) + 9./4.*d - d2/4. + d3/12.);
}

double Bsgamma::Phi88_1(double E0){
    double d=delta(E0);
    double d2=d*d;
    double d3=d2*d;
    
    double Li2 = gsl_sf_dilog(1. - d);
    
    double pi2=M_PI*M_PI;
    
    return 1./27.*( -2.*log(Mb1s/Ms)*( d2 + 2.*d + 4.*log(1. - d) ) + 4.*Li2 - 2./3.*pi2 - d*(2. + d)*log(d) 
            + 8.*log(1. - d) -  2./3.*d3 + 3.*d2 +7*d);
}

double Bsgamma::Kij_1(int i, int j, double E0, double mu){
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

void Bsgamma::computeCoeff(double mu){
    allcoeff = mySM.getMyFlavour()->ComputeCoeffsgamma(mu);
    
    C1_0 = -0.8411;//(*(allcoeff[LO]))(0);
    C2_0 = 1.0647;//(*(allcoeff[LO]))(1);
    C3_0 = -0.0133;//(*(allcoeff[LO]))(2);
    C4_0 = -0.1276;//(*(allcoeff[LO]))(3);
    C5_0 = 0.0012;//(*(allcoeff[LO]))(4);
    C6_0 = 0.0028;//(*(allcoeff[LO]))(5);
    C7_0 = -0.3736;//(*(allcoeff[LO]))(6);
    C8_0 = -0.1729;//(*(allcoeff[LO]))(7);
    
    C1_1 = 15.278;//(*(allcoeff[NLO]))(0);
    C2_1 = -2.124;//(*(allcoeff[NLO]))(1);
    C3_1 = 0.096;//(*(allcoeff[NLO]))(2);
    C4_1 = -0.463;//(*(allcoeff[NLO]))(3);
    C5_1 = -0.021;//(*(allcoeff[NLO]))(4);
    C6_1 = -0.013;//(*(allcoeff[NLO]))(5);
    C7_1 = 2.027;//(*(allcoeff[NLO]))(6);
    C8_1 = -0.617;//(*(allcoeff[NLO]))(7);
    
    //C7_2 = (*(allcoeff[NNLO]))(6);
}

double Bsgamma::P21(double E0, double mu){
    return (C1_0*C1_0).real() * Kij_1(1,1,E0,mu) + 2.*(C1_0*C2_0).real() * Kij_1(1,2,E0,mu) 
            + 2.*(C1_0*C7_0).real() * Kij_1(1,7,E0,mu) + 2.*(C1_0*C8_0).real() * Kij_1(1,8,E0,mu)
            + (C2_0*C2_0).real() * Kij_1(2,2,E0,mu) + 2.*(C2_0*C7_0).real() * Kij_1(2,7,E0,mu) 
            + 2.*(C2_0*C8_0).real() * Kij_1(2,8,E0,mu) + 2.*(C3_0*C7_0).real() * Kij_1(3,7,E0,mu)
            + 2.*(C4_0*C7_0).real() * Kij_1(4,7,E0,mu) + 2.*(C4_0*C8_0).real() * Kij_1(4,8,E0,mu)
            + 2.*(C5_0*C7_0).real() * Kij_1(5,7,E0,mu) + 2.*(C6_0*C7_0).real() * Kij_1(6,7,E0,mu) 
            + (C7_0*C7_0).real() * Kij_1(7,7,E0,mu) + 2.*(C7_0*C8_0).real() * Kij_1(7,8,E0,mu) 
            + (C8_0*C8_0).real() * Kij_1(8,8,E0,mu);
}

double Bsgamma::P32(double E0, double mu){
    return 2.*( (C1_0*C1_1).real() * Kij_1(1,1,E0,mu) + (C1_0*C2_1 + C1_1*C2_0).real() * Kij_1(1,2,E0,mu) 
            + (C1_0*C7_1 + C1_1*C7_0).real() * Kij_1(1,7,E0,mu) + (C1_0*C8_1 + C1_1*C8_0).real() * Kij_1(1,8,E0,mu)
            + (C2_0*C2_1).real() * Kij_1(2,2,E0,mu) + (C2_0*C7_1 + C2_1*C7_0).real() * Kij_1(2,7,E0,mu) 
            + (C2_0*C8_1 + C2_1*C8_0).real() * Kij_1(2,8,E0,mu) + (C3_0*C7_1 + C3_1*C7_0).real() * Kij_1(3,7,E0,mu)
            + (C4_0*C7_1 + C4_1*C7_0).real() * Kij_1(4,7,E0,mu) + (C4_0*C8_1 + C4_1*C8_0).real() * Kij_1(4,8,E0,mu) 
            + (C5_0*C7_1 + C5_1*C7_0).real() * Kij_1(5,7,E0,mu) + (C6_0*C7_1 + C6_1*C7_0).real() * Kij_1(6,7,E0,mu) 
            + (C7_0*C7_1).real() * Kij_1(7,7,E0,mu) + (C7_0*C8_1 + C7_1*C8_0).real() * Kij_1(7,8,E0,mu) 
            + (C8_0*C8_1).real() * Kij_1(8,8,E0,mu) );
}

double Bsgamma::P(double E0, double mu, orders order){
    switch(order) {
        case NNLO:
            /*return C7_0.abs2() + SM.Als(mu,NNLO)/4./M_PI * (2.*(C7_0*C7_1).real() + P21(E0,mu)) + 
                    pow(SM.Als(mu,NNLO)/4./M_PI, 2.) * (C7_1.abs2() + 2.*(C7_0*C7_2).real() + P22(E0,mu) + P32(E0,mu));
            break;*/
        case NLO:
            return C7_0.abs2() + SM.Als(mu,NLO)/4./M_PI * (2.*(C7_0*C7_1).real() + P21(E0,mu));
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

void Bsgamma::computeBR(orders order){
    ale=mySM.getAle();
    E0=mySM.getbsgamma_E0();
    mu_b = mySM.getMub();
    Mb1s=4.68;
    Mc=mySM.getQuarks(QCD::CHARM).getMass();
    Ms=mySM.getQuarks(QCD::STRANGE).getMass();
    BRsl=mySM.getBr_B_Xcenu();
    C=0.568;
    lambda_t=mySM.computelamt_s();
    V_cb=0.04174304221;
    
    computeCoeff(mu_b);
    
    BR = BRsl * (lambda_t/V_cb).abs2() * 6. * ale / (M_PI * C) * P(E0, mu_b, order);
    BR_conj = BRsl * (lambda_t/V_cb).conjugate().abs2() * 6. * ale / (M_PI * C) * P(E0, mu_b, order);
}

double Bsgamma::computeThValue(){
    computeBR(NLO);
    
    if (obs == 1) return BR;
    if (obs == 2) return (BR - BR_conj) / (BR + BR_conj);
    
    throw std::runtime_error("Bsgamma::computeThValue(): Observable type not defined. Can be only any of (1,2)");
}