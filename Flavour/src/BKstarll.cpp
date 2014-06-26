/* 
<<<<<<< HEAD
 * Copyright (C) 2014 SusyFit Collaboration
 * All rights reserved.
 *
 * For the licensing terms see doc/COPYING.
=======
 * File:   BKstarll.cpp
 * Author: marco
 * 
 * Created on 16 giugno 2014, 15.30
>>>>>>> 97e15c2d836f31cf6994885995100219a3d7c67a
 */

#include "BKstarll.h"
#include <gslpp.h>
#include <gslpp_complex.h>





<<<<<<< HEAD
BKstarll::BKstarll(const StandardModel& SM_i, int lep_i) : ThObservable(SM_i), mySM(dynamic_cast<const StandardModel&> (SM_i)) {
    GF = mySM.getGF();    
    ale=mySM.getAle();
    switch(lep_i){
        case '0': Mm=mySM.getLeptons(StandardModel::ELECTRON).getMass();
        case '1': Mm=mySM.getLeptons(StandardModel::MU).getMass();
        default:
            std::stringstream out;
            out << lep_i;
            throw std::runtime_error("Lep_i: index " + out.str() + "not implemented");
    }
    MB=mySM.getMesons(QCD::B_D).getMass();
    Mb=mySM.getQuarks(QCD::BOTTOM).getMass();
    Ms=mySM.getQuarks(QCD::STRANGE).getMass();
    MW=mySM.Mw();
    lambda_t=mySM.computelamt();
    b=1.;                           //please check
=======
BKstarll::BKstarll(const StandardModel& SM_i) : ThObservable(SM_i), mySM(dynamic_cast<const StandardModel&> (SM_i)) {

>>>>>>> 97e15c2d836f31cf6994885995100219a3d7c67a
}


BKstarll::~BKstarll() {
}



/*******************************************************************************
 * Helicity amplitudes                                                         *
 * ****************************************************************************/



<<<<<<< HEAD
complex BKstarll::N(){
    return -(4*GF*MB*ale*lambda_t)/(sqrt(2)*4*M_PI);
=======
double BKstarll::N(){
    return -(4*GF*MB*e*e*lambda_t)/(sqrt(2)*16*M_PI*M_PI);
>>>>>>> 97e15c2d836f31cf6994885995100219a3d7c67a
}


/**
 * polarization   i
 *      0         0
 *      +         1
 *      -         2
 */


gslpp::complex BKstarll::H_V(int i) {
<<<<<<< HEAD
    return -gslpp::complex::i()*N()*( C9*V_L(i) + C9p*V_R(i) + MB*MB/q2*( 2*Mb/MB*( C7*T_L(i) 
=======
    return -gslpp::complex::i()*N()*( C9*V_L(i) + C9p*V_R(i) + MB*MB/q2*( 2*Mb_MSB/MB*( C7*T_L(i) 
>>>>>>> 97e15c2d836f31cf6994885995100219a3d7c67a
            + C7p*T_R(i) ) - 16*M_PI*M_PI*h_lambda ) );
}



gslpp::complex BKstarll::H_A(int i) {
    return -gslpp::complex::i()*N()*( C10*V_L(i) + C10p*V_R(i) );
}



gslpp::complex BKstarll::H_S() {
<<<<<<< HEAD
    return gslpp::complex::i()*N()*Mb/MW*( CS*S_L + CSp*S_R );
=======
    return gslpp::complex::i()*N()*Mb_MSB/MW*( CS*S_L + CSp*S_R );
>>>>>>> 97e15c2d836f31cf6994885995100219a3d7c67a
}



gslpp::complex BKstarll::H_P() {
<<<<<<< HEAD
    return gslpp::complex::i()*N()*( Mb/MW*( CS*S_L + CSp*S_R ) + 2*Mm*Mb/q2*( C10*( S_L - Ms/Mb*S_R ) + C10p*( S_R - Ms/Mb*S_L ) ) );
=======
    return gslpp::complex::i()*N()*( Mb_MSB/MW*( CS*S_L + CSp*S_R ) + 2*Mm*Mb_MSB/q2*( C10*( S_L - Ms/Mb*S_R ) + C10p*( S_R - Ms/Mb*S_L ) ) );
>>>>>>> 97e15c2d836f31cf6994885995100219a3d7c67a
}



/*******************************************************************************
<<<<<<< HEAD
 * Angular coefficients                                                         *
=======
 * Anguar coefficients                                                         *
>>>>>>> 97e15c2d836f31cf6994885995100219a3d7c67a
 * ****************************************************************************/



double BKstarll::beta(double q2) {
    return sqrt(1-4*Mm*Mm/q2);
}



double BKstarll::lambda(double k2) {
    return 4*MB*MB*k2;
}



double BKstarll::F(double q2, double b) {
    return sqrt(lambda(k2))*beta(q2)*q2*b/(96*M_PI*M_PI*M_PI*MB*MB*MB);   
}



double BKstarll::I(int i) {

    double Mm2 = Mm*Mm;
    double beta2 = beta(q2)*beta(q2);
    
    switch (i){
        case '0':
            return F(q2,b)*( ( H_V(0).abs2() + H_A(0).abs2() )/2  +  H_P().abs2()  +  2*Mm2/q2*( H_V(0).abs2() - H_A(0).abs2() )  +  
                beta2*H_S().abs2() ); 
        case '1':
            return F(q2,b)*( (beta2 + 2)/8*( H_V(1).abs2() + H_V(2).abs2() + H_A(1).abs2() + H_A(2).abs2() )  +  
                Mm2/q2*( H_V(1).abs2() + H_V(2).abs2() - H_A(1).abs2() - H_A(2).abs2() ) );
        case '2':
            return -F(q2,b)*beta2/2*( ( H_V(0).abs2() + H_A(0).abs2() )/2 );
        case '3':
            return F(q2,b)*beta2/8*( ( H_V(1).abs2() + H_V(2).abs2() )/2  +  ( H_A(1).abs2() + H_A(2).abs2() )/2 );
        case '4':
            return -F(q2,b)/2*( ( H_V(1)*H_V(2).conjugate() ).real()  +  ( H_A(1)*H_A(2).conjugate() ).real() );
        case '5':
            return F(q2,b)*beta2/4*( ( (H_V(2) + H_V(1))*H_V(0).conjugate() ).real()  +  ( (H_A(2) + H_A(1))*H_A(0).conjugate() ).real() );
        case '6':
            return F(q2,b)*( beta(q2)/2*( ( (H_V(2) - H_V(1))*H_A(0).conjugate() ).real()  +  ( (H_A(2) - H_A(1))*H_V(0).conjugate() ).real() )  -
                beta(q2)*Mm/sqrt(q2)*( H_S().conjugate()*(H_V(1) + H_V(2)) ).real() );
        case '7':
            return F(q2,b)*beta(q2)*( H_V(2)*H_A(2).conjugate() - H_V(1)*H_A(1).conjugate() ).real();
        case '8':
            return 2*F(q2,b)*beta(q2)*Mm/sqrt(q2)*( H_S().conjugate()*H_V(0) ).real();
        case '9':
            return F(q2,b)*( beta(q2)/2*( ( (H_V(2) + H_V(1))*H_A(0).conjugate() ).imag()  +  ( (H_A(2) + H_A(1))*H_V(0).conjugate() ).imag() )  -
                beta(q2)*Mm/sqrt(q2)*( H_S().conjugate()*(H_V(2) - H_V(1)) ).imag() );
        case '10':
            return F(q2,b)*beta2/4*( ( (H_V(2) - H_V(1))*H_V(0).conjugate() ).imag()  +  ( (H_A(2) - H_A(1))*H_A(0).conjugate() ).imag() );
        case '11':
            return F(q2,b)*beta2/2*( ( H_V(1)*H_V(2).conjugate() ).imag()  +  ( H_A(1)*H_A(2).conjugate() ).imag() );
        default:
            std::stringstream out;
            out << i;
            throw std::runtime_error("I: index " + out.str() + "not implemented");
    }
}



double BKstarll::I_bar(int i) {
    
    if ( i==0 | i==1 | i==2 | i==3 | i==4 | i==5 | i==9)
    {
        return I(i);
    }
    else if ( i==6 | i==7 | i==8 | i==10 | i==11)
    {
        return -I(i);
    }
    else{
        std::stringstream out;
        out << i;
        throw std::runtime_error("I: index " + out.str() + "not implemented");
    }   
}




double BKstarll::Sigma(int i) {
    return (I(i) + I_bar(i))/2;
}



/*******************************************************************************
 * Observables                                                                 *
 * ****************************************************************************/



<<<<<<< HEAD
P_1::P_1(const StandardModel& SM_i, int lep_i) : BKstarll(SM_i, lep_i) {  
=======
P_1::P_1(const StandardModel& SM_i) : BKstarll(SM_i) {  
>>>>>>> 97e15c2d836f31cf6994885995100219a3d7c67a
}

double P_1::computeThValue() {
    return Sigma(4)/(2.*Sigma(3));
    
}


<<<<<<< HEAD
P_2::P_2(const StandardModel& SM_i, int lep_i) : BKstarll(SM_i, lep_i) {  
=======
P_2::P_2(const StandardModel& SM_i) : BKstarll(SM_i) {  
>>>>>>> 97e15c2d836f31cf6994885995100219a3d7c67a
}

double P_2::computeThValue() {
    return Sigma(7)/(8.*Sigma(3));
}


<<<<<<< HEAD
P_3::P_3(const StandardModel& SM_i, int lep_i) : BKstarll(SM_i, lep_i) {  
=======
P_3::P_3(const StandardModel& SM_i) : BKstarll(SM_i) {  
>>>>>>> 97e15c2d836f31cf6994885995100219a3d7c67a
}

double P_3::computeThValue() {
    return -Sigma(11)/(4.*Sigma(3));

}


<<<<<<< HEAD
P_4Prime::P_4Prime(const StandardModel& SM_i, int lep_i) : BKstarll(SM_i, lep_i) {  
=======
P_4Prime::P_4Prime(const StandardModel& SM_i) : BKstarll(SM_i) {  
>>>>>>> 97e15c2d836f31cf6994885995100219a3d7c67a
}

double P_4Prime::computeThValue() {
    return Sigma(5)/sqrt(-Sigma(2)*Sigma(3));
   
}


<<<<<<< HEAD
P_5Prime::P_5Prime(const StandardModel& SM_i, int lep_i) : BKstarll(SM_i, lep_i) {  
=======
P_5Prime::P_5Prime(const StandardModel& SM_i) : BKstarll(SM_i) {  
>>>>>>> 97e15c2d836f31cf6994885995100219a3d7c67a
}

double P_5Prime::computeThValue() {
    return Sigma(6)/(2.*sqrt(-Sigma(2)*Sigma(3)));

}


<<<<<<< HEAD
P_6Prime::P_6Prime(const StandardModel& SM_i, int lep_i) : BKstarll(SM_i, lep_i) {  
=======
P_6Prime::P_6Prime(const StandardModel& SM_i) : BKstarll(SM_i) {  
>>>>>>> 97e15c2d836f31cf6994885995100219a3d7c67a
}

double P_6Prime::computeThValue() {
    return -Sigma(9)/(2.*sqrt(-Sigma(2)*Sigma(3)));
 
}


<<<<<<< HEAD
GammaPrime::GammaPrime(const StandardModel& SM_i, int lep_i) : BKstarll(SM_i, lep_i) {  
=======
GammaPrime::GammaPrime(const StandardModel& SM_i) : BKstarll(SM_i) {  
>>>>>>> 97e15c2d836f31cf6994885995100219a3d7c67a
}

double GammaPrime::computeThValue() {
    return ((3.*Sigma(0) - Sigma(2)) + 2.*(3*Sigma(1) - Sigma(3)))/4.;

}


<<<<<<< HEAD
F_L::F_L(const StandardModel& SM_i, int lep_i) : BKstarll(SM_i, lep_i), mySM(SM_i) {
    lep = lep_i;
}

double F_L::computeThValue() {
    GammaPrime myGammaPrime(mySM, lep);
=======
F_L::F_L(const StandardModel& SM_i) : BKstarll(SM_i), mySM(SM_i) {  
}

double F_L::computeThValue() {
    GammaPrime myGammaPrime(mySM);
>>>>>>> 97e15c2d836f31cf6994885995100219a3d7c67a
    return (3.*Sigma(0) - Sigma(2))/(4.*myGammaPrime.computeThValue());

}


<<<<<<< HEAD
M_1Prime::M_1Prime(const StandardModel& SM_i, int lep_i) : BKstarll(SM_i, lep_i) {  
=======
M_1Prime::M_1Prime(const StandardModel& SM_i) : BKstarll(SM_i) {  
>>>>>>> 97e15c2d836f31cf6994885995100219a3d7c67a
}

double M_1Prime::computeThValue() {
    return ( H_V(1).abs2() + H_V(2).abs2() - H_A(1).abs2() - H_A(2).abs2() )/( 2*( H_V(1).abs2() + H_V(2).abs2() + H_A(1).abs2() + H_A(2).abs2() ) );
  
}


<<<<<<< HEAD
M_2Prime::M_2Prime(const StandardModel& SM_i, int lep_i) : BKstarll(SM_i, lep_i) {  
=======
M_2Prime::M_2Prime(const StandardModel& SM_i) : BKstarll(SM_i) {  
>>>>>>> 97e15c2d836f31cf6994885995100219a3d7c67a
}

double M_2Prime::computeThValue() {
    return ( q2/(2*Mm*Mm)*( H_P().abs2() + beta(q2)*beta(q2)*H_S().abs2() ) + H_V(0).abs2() - H_A(0).abs2() )/( H_V(0).abs2() + H_A(0).abs2() );  
}