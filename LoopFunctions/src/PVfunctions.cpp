/* 
 * File:   PVfunctions.cpp
 * Author: mishima
 * 
 * Created on June 10, 2011, 2:47 PM
 */

#include <iostream>
#include <cmath>
#include "PVfunctions.h"

PVfunctions::PVfunctions() {
}

PVfunctions::PVfunctions(const PVfunctions& orig) {
}

PVfunctions::~PVfunctions() {
}

double PVfunctions::A0(const double mu, const double m) {
    if ( mu<=0.0 || m<0.0 ) {
        throw "Invalid argument for PVfunctions::A0()";
    }
    return ( -m*m*(-log(m*m/mu/mu)+1.0) );    
}

complex PVfunctions::B0(const double mu, const double p2, 
                        const double m0, const double m1) {   
    if ( mu<=0.0 || p2<0.0 || m0<0.0 || m1<0.0 ) {
        throw "Invalid argument for PVfunctions::B0()";
    }
    double mu2=mu*mu, m02=m0*m0, m12=m1*m1;
    complex B0(0.0, 0.0, false);
        
    if ( p2!=0.0 && m0!=0.0 && m1!=0.0 ) { 
        double Lambda = sqrt( abs((p2-m02-m12)*(p2-m02-m12) - 4.0*m02*m12) );
        double R;
        if ( p2 > (m0-m1)*(m0-m1) && p2 < (m0+m1)*(m0+m1) ) {
            R = - Lambda/p2*atan(Lambda/(p2-m02-m12));
        } else {
            R = Lambda/p2*log( abs((p2-m02-m12+Lambda)/2.0/m0/m1) );
        }               
        B0 = - log(m0*m1/mu2) - (m02-m12)/2.0/p2*log(m12/m02) - R + 2.0;
        if ( p2>(m0+m1)*(m0+m1) ) { 
            B0 += M_PI*Lambda/p2*complex::i();// imaginary part
        }
    } else if ( p2==0.0 && m0!=0.0 && m1!=0.0 ) {
        if ( m0!=m1 ) {
            B0 = - m02/(m02-m12)*log(m02/mu2) + m12/(m02-m12)*log(m12/mu2) + 1.0;
        } else {
            B0 = - log(m02/mu2);
        }
    } else if ( p2!=0.0 && ((m0==0.0 && m1!=0.0) || (m0!=0.0 && m1==0.0)) ) {
        double M2;
        if (m0!=0.0) M2 = m02;
        if (m1!=0.0) M2 = m12;
        B0 = - log(M2/mu2) + 2.0;
        if ( p2<M2 ) {
            B0 += - (1.0 - M2/p2)*log(1.0 - p2/M2);
        } else if ( p2>M2 ) {
            B0 += - (1.0 - M2/p2)*log(p2/M2 - 1.0);
            B0 += (1.0 - M2/p2)*M_PI*complex::i();// imaginary part        
        } else {
            // the case of p2=M2
            // B0 += 0.0;
            throw "Write codes for B0(M2,0,M2) in PVfunctions::B0().";
        }
    } else if ( p2==0.0 && m0==0.0 && m1!=0.0 ) {
        B0 = - log(m12/mu2) + 1.0;
    } else if ( p2==0.0 && m0!=0.0 && m1==0.0 ) {    
        B0 = - log(m02/mu2) + 1.0;        
    } else if ( p2!=0.0 && m0==0.0 && m1==0.0 ) {
        if ( p2<0.0 ) {
            B0 = - log(-p2/mu2) + 2.0;
        } else {
            B0 = - log(p2/mu2) + 2.0;
            B0 += M_PI*complex::i();// imaginary part    
        }
    } else if ( p2==0.0 && m0==0.0 && m1==0.0 ) {           
        throw "PVfunctions::B0() is IR divergent.";            
    } else {
        throw "There is a missing case in the codes of PVfunctions::B0().";
    }
    return B0;
}

complex PVfunctions::C0(const double p2, 
                        const double m0, const double m1, const double m2) {
    if ( p2<0.0 || m0<0.0 || m1<0.0 || m2<0.0 ) {
        throw "Invalid argument for PVfunctions::C0()";
    }
    double m02=m0*m0, m12=m1*m1, m22=m2*m2;
    complex C0(0.0, 0.0, false);
    
    throw "Write codes for PVfunctions::C0()";

    
    
    
    
    return C0;
}
        