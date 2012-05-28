/* 
 * File:   PVfunctions.cpp
 * Author: mishima
 */

#include <iostream>
#include <cmath>
#include <gsl/gsl_complex.h>
#include <gsl/gsl_sf.h>
#include "PVfunctions.h"
#define LEPS 1.e-5 // tolerance in the limit of masses


PVfunctions::PVfunctions() {
}

double PVfunctions::A0(const double mu, const double m) const {
    if ( mu<=0.0 || m<0.0 ) {
        throw "Invalid argument for PVfunctions::A0()";
    }
    
    if (m==0.0) {
        return ( 0.0 );
    } else {
        return ( -m*m*(-log(m*m/mu/mu)+1.0) ); 
    }
}

complex PVfunctions::B0(const double mu, const double p2, 
                        const double m0, const double m1) const {   
    if ( mu<=0.0 || p2<0.0 || m0<0.0 || m1<0.0 ) {
        throw "Invalid argument for PVfunctions::B0()";
    }
    double mu2=mu*mu, m02=m0*m0, m12=m1*m1;
    complex B0(0.0, 0.0, false);
        
    if ( p2!=0.0 && m0!=0.0 && m1!=0.0 ) { 
        double Lambda = sqrt( fabs((p2-m02-m12)*(p2-m02-m12) - 4.0*m02*m12) );
        double R;
        if ( p2 > (m0-m1)*(m0-m1) && p2 < (m0+m1)*(m0+m1) ) {
            R = - Lambda/p2*atan(Lambda/(p2-m02-m12));
        } else {
            R = Lambda/p2*log( fabs((p2-m02-m12+Lambda)/2.0/m0/m1) );
        }               
        B0 = - log(m0*m1/mu2) + (m02-m12)/2.0/p2*log(m12/m02) - R + 2.0;
        if ( p2>(m0+m1)*(m0+m1) ) { 
            B0 += M_PI*Lambda/p2*complex::i();// imaginary part
        }
    } else if ( p2==0.0 && m0!=0.0 && m1!=0.0 ) {
        if ( fabs(m0 - m1) > LEPS ) {             ///////////////////////////////////////////////
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
            B0 += 0.0;
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
        throw "PVfunctions::B0() is IR divergent. (vanishes in DR)";            
    } else {
        throw "Missing case in the codes of PVfunctions::B0().";
    }
    return B0;
}

complex PVfunctions::B1(const double mu, const double p2, 
                        const double m0, const double m1) const {   
    if ( mu<=0.0 || p2<0.0 || m0<0.0 || m1<0.0 ) {
        throw "Invalid argument for PVfunctions::B1()";
    }
    double mu2=mu*mu, m02=m0*m0, m12=m1*m1;
    double DeltaM2 = m02 - m12; 
    complex B1(0.0, 0.0, false);
    
    if (p2==0.0) {
        if (m02!=0.0 && m12!=0.0) {
            if (fabs(m02 - m12) > LEPS) {////////////////////////////////////////////////////////////////////////
                double F0 = - log(m12/m02);
                double F1 = - 1.0 + m02/DeltaM2*F0;
                double F2 = - 1.0/2.0 + m02/DeltaM2*F1;
                B1.real() = 1.0/2.0*( log(m12/mu2) + F2 );
            } else {
                B1.real() = 1.0/2.0*log(m12/mu2);                
            }
        } else if (m02==0.0 && m12!=0.0) {
            B1.real() = 1.0/2.0*log(m12/mu2) - 1.0/4.0;
        } else if (m02!=0.0 && m12==0.0) { 
            B1.real() = 1.0/2.0*log(m02/mu2) - 1.0/4.0;
        } else {                
            B1 = 0.0;
        }
    } else {
        B1 = -1.0/2.0/p2*(A0(mu,m0) - A0(mu,m1) 
                          + (DeltaM2 + p2)*B0(mu,p2,m0,m1));
    }
    return B1;
}

complex PVfunctions::B21(const double mu, const double p2, 
                         const double m0, const double m1) const {   
    if ( mu<=0.0 || p2<0.0 || m0<0.0 || m1<0.0 ) {
        throw "Invalid argument for PVfunctions::B21()";
    }
    double mu2=mu*mu, m02=m0*m0, m12=m1*m1;
    double DeltaM2 = m02 - m12; 
    complex B21(0.0, 0.0, false);

    if (p2==0.0) {
        if (m02!=0.0 && m12!=0.0) {
            if (fabs(m02 - m12) > LEPS ) {                                  //////////////////////////////
                double F0 = - log(m12/m02);
                double F1 = - 1.0 + m02/DeltaM2*F0;
                double F2 = - 1.0/2.0 + m02/DeltaM2*F1;
                double F3 = - 1.0/3.0 + m02/DeltaM2*F2;
                B21.real() = - 1.0/3.0*( log(m12/mu2) + F3 );        
            } else {
                B21.real() = - 1.0/3.0*log(m12/mu2);
            }
        } else if (m02==0.0 && m12!=0.0) {
            B21.real() = - 1.0/3.0*log(m12/mu2) + 1.0/9.0;
        } else if (m02!=0.0 && m12==0.0) { 
            B21.real() = - 1.0/3.0*log(m02/mu2) + 1.0/9.0;
        } else {                
            throw "PVfunctions::B21() is undefined."; 
        }        
    } else {
        double Lambdabar2 = (p2-m02-m12)*(p2-m02-m12) - 4.0*m02*m12;
        B21 = - (3.0*(m02 + m12) - p2)/18.0/p2 
              + (DeltaM2 + p2)/3.0/p2/p2*A0(mu,m0) 
              - (DeltaM2 + 2.0*p2)/3.0/p2/p2*A0(mu,m1) 
              + (Lambdabar2 + 3.0*p2*m02)/3.0/p2/p2*B0(mu,p2,m0,m1);
    }
    return B21;
}

complex PVfunctions::B22(const double mu, const double p2, 
                         const double m0, const double m1) const {   
    if ( mu<=0.0 || p2<0.0 || m0<0.0 || m1<0.0 ) {
        throw "Invalid argument for PVfunctions::B22()";
    }
    double mu2=mu*mu, m02=m0*m0, m12=m1*m1;
    double DeltaM2 = m02 - m12;
    complex B22(0.0, 0.0, false);

    if(p2 == 0.){
        if(m02 != 0. && m12 != 0.){
            if(fabs(m02 - m12) < LEPS){
                B22 =  m02 / 2. * (- log(m02 / mu2) + 1.);
            } else {
                B22 =  1. / 4. * (m02 + m12) * (- log(m0 * m1 / mu2) + 1.5) 
                      - (m02 * m02 + m12 * m12) / 8. / (m02 - m12) * log(m02 / m12);  
            }             
        } else
            throw "PVfunctions::B22() is undefined."; 
    } else {
        if(m0 != 0. && m1 != 0.){
            if(fabs(m02 - m12) < LEPS){
              B22 =  (6. * m02 - p2) / 18. - A0(mu,m0) /6. 
                    - (p2 - 4. * m02) / 12. * B0(mu,p2,m0,m1);  
            } else {
              B22 =  (3. * (m02 + m12) - p2) / 18. 
                    - (DeltaM2 + p2) / 12. / p2 * A0(mu,m0) 
                    + (DeltaM2 - p2) / 12. / p2 * A0(mu,m1) 
                    - (- m02 - m12 + p2 / 2. 
                       + (DeltaM2 * DeltaM2) / 2. / p2) * B0(mu,p2,m0,m1) / 6.;                
            }
        } else
            throw "PVfunctions::B22() is undefined."; 
    }
    return B22;
}

complex PVfunctions::Bf(const double mu, const double p2, 
                        const double m0, const double m1) const {   
    if ( mu<=0.0 || p2<0.0 || m0<0.0 || m1<0.0 ) {
        throw "Invalid argument for PVfunctions::Bf()";
    }
    complex Bf(0.0, 0.0, false);
    Bf = 2.0*(B21(mu,p2,m0,m1) + B1(mu,p2,m0,m1));
    return Bf;
}

complex PVfunctions::B0p(const double mu, const double p2, 
                         const double m0, const double m1) const {   
    if ( mu<=0.0 || p2<0.0 || m0<0.0 || m1<0.0 ) {
        throw "Invalid argument for PVfunctions::B0p()";
    }
    double mu2=mu*mu, m02=m0*m0, m12=m1*m1;
    complex B0p(0.0, 0.0, false);
        
    if (p2==0.0) {
        if ( m0!=0.0 && m1!=0.0 ) {                 
            double DeltaM2 = m02 - m12;
            B0p = (m02 + m12)/2.0/pow(DeltaM2,2.0)
                   + m02*m12/pow(DeltaM2,3.0)*log(m12/m02);
        } else {        
            throw "PVfunctions::B0p() is undefined.";
        }
    } else {
        if ( m0!=0.0 && m1!=0.0 ) {         
            double Lambda = sqrt( fabs((p2-m02-m12)*(p2-m02-m12) - 4.0*m02*m12) );
            double Rprime;
            if ( p2 > (m0-m1)*(m0-m1) && p2 < (m0+m1)*(m0+m1) ) {
                Rprime = ((p2 - m02 - m12)/Lambda + Lambda/p2)
                          *atan(Lambda/(p2-m02-m12));
            } else {
                Rprime = ((p2 - m02 - m12)/Lambda - Lambda/p2)
                          *log( fabs((p2-m02-m12+Lambda)/2.0/m0/m1) );
            }                  
            B0p = - (m02 - m12)/2.0/p2/p2*log(m12/m02) - (Rprime + 1.0)/p2;
            if ( p2>(m0+m1)*(m0+m1) ) { 
                B0p += M_PI/p2*((p2 - m02 - m12)/Lambda - Lambda/p2)
                      *complex::i();// imaginary part
            }   
        } else if ( (m0==0.0 && m1!=0.0) || (m0!=0.0 && m1==0.0) ) {
            double M2;
            if (m0!=0.0) M2 = m02;
            if (m1!=0.0) M2 = m12;            
            if ( p2<M2 ) {
                B0p = - M2/p2/p2*log(1.0 - p2/M2) - 1.0/p2;
            } else if ( p2>M2 ) {
                B0p = - M2/p2/p2*log(p2/M2 - 1.0) - 1.0/p2;
                B0p += M2/p2/p2*M_PI*complex::i();// imaginary part        
            } else { /* p2=M2 */
                B0p = 1.0/2.0/M2*(log(M2/mu2) - 2.0);
            }
        } else if ( m0==0.0 && m1==0.0 ) {
            B0p = - 1.0/p2;
        } else {
            throw "PVfunctions::B0p() is undefined.";
        }
    }
    return B0p;
}

complex PVfunctions::B1p(const double mu, const double p2, 
                         const double m0, const double m1) const {   
    if ( mu<=0.0 || p2<0.0 || m0<0.0 || m1<0.0 ) {
        throw "Invalid argument for PVfunctions::B1p()";
    }
    complex B1p(0.0, 0.0, false);
    
    if (p2==0.0) {
        throw "PVfunctions::B1p() is undefined.";                     
    } else {
        double DeltaM2 = m0*m0 - m1*m1;
        B1p = - ( 2.0*B1(mu,p2,m0,m1) + B0(mu,p2,m0,m1) 
                 + (DeltaM2 + p2)*B0p(mu,p2,m0,m1) )/2.0/p2;
    }
    return B1p;
}

complex PVfunctions::B21p(const double mu, const double p2, 
                          const double m0, const double m1) const {   
    if ( mu<=0.0 || p2<0.0 || m0<0.0 || m1<0.0 ) {
        throw "Invalid argument for PVfunctions::B21p()";
    }
    double m02=m0*m0, m12=m1*m1;
    double p4 = p2*p2, p6=p2*p2*p2;
    double DeltaM2 = m02 - m12; 
    complex B21p(0.0, 0.0, false);

    if (p2==0.0) {
        throw "PVfunctions::B21p() is undefined.";           
    } else {
        double Lambdabar2 = (p2-m02-m12)*(p2-m02-m12) - 4.0*m02*m12;
        B21p = (m02 + m12)/6.0/p4 - (2.0*DeltaM2 + p2)/3.0/p6*A0(mu,m0)
                + 2.0*(DeltaM2 + p2)/3.0/p6*A0(mu,m1)
//                + (2.0*m12 - m02)/3.0/p4*B0(mu,p2,m0,m1) 
                - (2.0*DeltaM2*DeltaM2 + p2*m02 - 2.0*p2*m12)/3.0/p6*B0(mu,p2,m0,m1) 
                + (Lambdabar2 + 3.0*p2*m02)/3.0/p4*B0p(mu,p2,m0,m1);
    }
    return B21p;
}

//complex PVfunctions::B22p(const double mu, const double p2, 
//                          const double m0, const double m1) {   
//    if ( mu<=0.0 || p2<0.0 || m0<0.0 || m1<0.0 ) {
//        throw "Invalid argument for PVfunctions::B22p()";
//    }
//    double mu2=mu*mu, m02=m0*m0, m12=m1*m1;
//    complex B22p(0.0, 0.0, false);
//
//    return B22p;
//}

complex PVfunctions::Bfp(const double mu, const double p2, 
                         const double m0, const double m1) const {   
    if ( mu<=0.0 || p2<0.0 || m0<0.0 || m1<0.0 ) {
        throw "Invalid argument for PVfunctions::Bfp()";
    }
    complex Bfp(0.0, 0.0, false);
    Bfp = 2.0*(B21p(mu,p2,m0,m1) + B1p(mu,p2,m0,m1));
    return Bfp;
}

complex PVfunctions::C0(const double p2, 
                        const double m0, const double m1, const double m2) const {
    if ( p2<0.0 || m0<0.0 || m1<0.0 || m2<0.0 ) {
        throw "Invalid argument for PVfunctions::C0()";
    }

    complex C0(0.0, 0.0, false);   
    if (p2==0.0) {
        throw "PVfunctions::C0() is undefined.";
    } else {
        if (fabs(m0 - m2) < LEPS && fabs(m0 - m1) > LEPS) {///////////////////////////////////
            double m02 = m0*m0;
            double m12 = m1*m1;
            double epsilon = 1.0e-12;
            gsl_complex tmp = gsl_complex_rect(1.0 - 4.0*m02/p2, epsilon);
            tmp = gsl_complex_sqrt(tmp);
            complex tmp_complex(GSL_REAL(tmp), GSL_IMAG(tmp), false);
            complex x0 = 1.0 - (m02 - m12)/p2;
            complex x1 = (1.0 + tmp_complex)/2.0;
            complex x2 = (1.0 - tmp_complex)/2.0;            
            complex x3 = m02/(m02 - m12);

            if ( x0==x1 || x0==x2 || x0==x3) {///////////////////////////////////////////////????????
                throw "PVfunctions::C0() is undefined.";
            }

            complex arg[6];
            arg[0] = (x0 - 1.0)/(x0 - x1);
            arg[1] = x0/(x0 - x1);
            arg[2] = (x0 - 1.0)/(x0 - x2);
            arg[3] = x0/(x0 - x2);
            arg[4] = (x0 - 1.0)/(x0 - x3);
            arg[5] = x0/(x0 - x3);
            
            complex Li2[6];
            for (int i=0; i<6; i++) {
                gsl_sf_result re, im;
                gsl_sf_complex_dilog_xy_e(arg[i].real(), arg[i].imag(), &re, &im);
                Li2[i].real() = re.val;
                Li2[i].imag() = im.val;
                
                /* Check the sizes of errors */
                //std::cout << "re.val=" << re.val << "  re.err=" << re.err << std::endl;
                //std::cout << "im.val=" << im.val << "  im.err=" << im.err << std::endl;                
            }
            C0 = - 1.0/p2*( Li2[0] - Li2[1] + Li2[2] - Li2[3] - Li2[4] + Li2[5]);        
        } else if (m0!=0.0 && m2!=0.0 && fabs(m0 - m2) > LEPS && m1==0.0) {              ////////////////
            double m02 = m0*m0;
            double m22 = m2*m2;
            double epsilon = 1.0e-12;
            double tmp_real = pow((m02+m22-p2),2.0) - 4.0*m02*m22;
            gsl_complex tmp = gsl_complex_rect(tmp_real, epsilon);
            tmp = gsl_complex_sqrt(tmp);
            complex tmp_complex(GSL_REAL(tmp), GSL_IMAG(tmp), false);
            complex x1 = (p2 - m02 + m22 + tmp_complex)/2.0/p2;
            complex x2 = (p2 - m02 + m22 - tmp_complex)/2.0/p2;            

            if ( x1==0.0 || x1==1.0 || x2==0.0 || x2==1.0 ) { //////////////////////////////////???
                throw "PVfunctions::C0() is undefined.";
            }            
            
            complex arg1 = (x1 - 1.0)/x1;
            complex arg2 = x2/(x2 - 1.0);
            gsl_complex arg1_tmp = gsl_complex_rect(arg1.real(), arg1.imag());
            gsl_complex arg2_tmp = gsl_complex_rect(arg2.real(), arg2.imag());            
            
            C0.real() = - 1.0/p2*( GSL_REAL(gsl_complex_log(arg1_tmp))
                                    *GSL_REAL(gsl_complex_log(arg2_tmp))
                                   - GSL_IMAG(gsl_complex_log(arg1_tmp))
                                     *GSL_IMAG(gsl_complex_log(arg2_tmp)) );
            C0.imag() = - 1.0/p2*( GSL_REAL(gsl_complex_log(arg1_tmp))
                                    *GSL_IMAG(gsl_complex_log(arg2_tmp))
                                   + GSL_IMAG(gsl_complex_log(arg1_tmp))
                                      *GSL_REAL(gsl_complex_log(arg2_tmp)) );            
        } else {
            throw "PVfunctions::C0() is undefined.";            
        }
    }
    return C0;
}

double PVfunctions::F(const double m0, const double m1) const {
    double m12 = m1 * m1;
    double m02 = m0 * m0;
    
    if ( m0<=0.0 || m1<0.0 ) {
        throw "Invalid argument for PVfunctions::F()";
    }
    
    if(m0 == 0. && m1 != 0.) {
        return (0.5 * m12);
    } else if(m0 != 0. && m1 == 0.){
        return (0.5 * m02);
    } else if((m0 == 0. && m1 == 0.) || (fabs(m0-m1) < LEPS)){
        return (0.);
    } else if (m0 != 0 && m1 != 0){
        return (0.5 * (m02 + m12) - (m02 * m12) / (m02 - m12) * log(m02 / m12));
    }
}
