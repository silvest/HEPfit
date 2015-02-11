/* 
 * Copyright (C) 2015 SusyFit Collaboration
 * All rights reserved.
 *
 * For the licensing terms see doc/COPYING.
 */

#include "bsgamma.h"
#include <gsl/gsl_sf_dilog.h>
#include <gsl/gsl_sf_zeta.h>
#include <gsl/gsl_sf_clausen.h>

Bsgamma::Bsgamma(const StandardModel& SM_i, int obsFlag): ThObservable(SM_i), mySM(SM_i){
    if (obsFlag > 0 and obsFlag < 3) obs = obsFlag;
    else throw std::runtime_error("obsFlag in bsgamma can only be 1 (BR) or 2 (A_{CP})");
}

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
}

void Bsgamma::computeCoeff(orders order){
    allcoeff = mySM.getMyFlavour()->ComputeCoeffsgamma(mu_b);
    E0=1.6; //fixed to experimental cutoff
    
    switch(order) {
        case FULLNNLO: 
            C_7 = (*(allcoeff[LO]))(6) + (*(allcoeff[NLO]))(6);// + (*(allcoeff[NNLO]))(6);
            C_8 = (*(allcoeff[LO]))(7) + (*(allcoeff[NLO]))(7);// + (*(allcoeff[NNLO]))(7);
            coeff = C_7.abs2()*G77(FULLNNLO,E0,mu_b) + (C_7*C_8).real()*G78(FULLNNLO,E0,mu_b);
            break;
        case FULLNLO: 
            C_7 = (*(allcoeff[LO]))(6) + (*(allcoeff[NLO]))(6);
            C_8 = (*(allcoeff[LO]))(7) + (*(allcoeff[NLO]))(7);
            coeff = C_7.abs2()*G77(FULLNLO,E0,mu_b) + (C_7*C_8).real()*G78(FULLNLO,E0,mu_b);
            break;
        case LO: 
            C_7 = (*(allcoeff[LO]))(6);
            coeff = C_7.abs2()*G77(LO,E0,mu_b);
            break;
        default:
            std::stringstream out;
            out << order;
            throw std::runtime_error("Bsgamma::computeCoeff(): order " + out.str() + " not implemented");
    }
}

void Bsgamma::computeGamma(orders order){
    GF = mySM.getGF();
    ale=mySM.getAle();
    Mb=mySM.getQuarks(QCD::BOTTOM).getMass();  // both pole mass and running MSbar mass should be implemented
    Mc=mySM.getQuarks(QCD::CHARM).getMass();
    mu_b = mySM.getMub();
    width = mySM.getMesons(StandardModel::B_D).computeWidth();
    lambda_t=mySM.computelamt_s();
    
    CF=4./3.;
    CA=3.;
    TR=1./2.;
    
    NH=1.;
    NV=1.;
    NL=3.;
    
    computeCoeff(order);
    
    Gamma = GF*GF*pow(Mb,5.)*ale*lambda_t.abs2()/(32. * pow(M_PI,4.)) * coeff;
    Gamma_conj = GF*GF*pow(Mb,5.)*ale*(lambda_t.conjugate()).abs2()/(32. * pow(M_PI,4.)) * coeff;
}

double Bsgamma::computeThValue(){
    computeGamma(FULLNNLO);
    
    if (obs == 1) return Gamma/width;
    if (obs == 2) return (Gamma - Gamma_conj) / (Gamma + Gamma_conj);
    
    throw std::runtime_error("Bsgamma::computeThValue(): Observable type not defined. Can be only any of (1,2)");
}