/* 
 * File:   PVfunctions.cpp
 * Author: mishima
 */

#include <iostream>
#include <cmath>
#include <gsl/gsl_complex.h>
#include <gsl/gsl_sf.h>
#include "PVfunctions.h"
#include <stdexcept>

#define LEPS 1.e-7 // tolerance in the limit of masses


PVfunctions::PVfunctions() {
}

double PVfunctions::A0(const double mu, const double m) const {
    if ( mu<=0.0 || m<0.0 ) {
        throw std::runtime_error("Invalid argument for PVfunctions::A0()"); 
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
        throw std::runtime_error("Invalid argument for PVfunctions::B0()");    
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
        if ( fabs(m0 - m1) > LEPS ) {
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
        throw std::runtime_error("PVfunctions::B0() is IR divergent. (vanishes in DR)");             
    } else {
        throw std::runtime_error("Missing case in the codes of PVfunctions::B0()."); 
    }
    return B0;
}

complex PVfunctions::B1(const double mu, const double p2, 
                        const double m0, const double m1) const {   
    if ( mu<=0.0 || p2<0.0 || m0<0.0 || m1<0.0 ) {
        throw std::runtime_error("Invalid argument for PVfunctions::B1()"); 
    }
    double mu2=mu*mu, m02=m0*m0, m12=m1*m1;
    double DeltaM2 = m02 - m12; 
    complex B1(0.0, 0.0, false);
    
    if (p2==0.0) {
        if (m02!=0.0 && m12!=0.0) {
            if (fabs(m02 - m12) > LEPS) {
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
        throw std::runtime_error("Invalid argument for PVfunctions::B21()"); 
    }
    double mu2=mu*mu, m02=m0*m0, m12=m1*m1;
    double DeltaM2 = m02 - m12; 
    complex B21(0.0, 0.0, false);

    if (p2==0.0) {
        if (m02!=0.0 && m12!=0.0) {
            if (fabs(m02 - m12) > LEPS ) {
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
            throw std::runtime_error("PVfunctions::B21() is undefined.");  
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
        throw std::runtime_error("Invalid argument for PVfunctions::B22()"); 
    }
    double mu2=mu*mu, m02=m0*m0, m12=m1*m1;
    double DeltaM2 = m02 - m12;
    complex B22(0.0, 0.0, false);

    if(p2 == 0.){
        if(m02 != 0. && m12 != 0.){
            if(fabs(m02 - m12) < LEPS){
                B22 = - m02 / 2. * (- log(m02 / mu2) + 1.);
            } else {
                B22 = - 1. / 4. * (m02 + m12) * (- log(m0 * m1 / mu2) + 1.5) 
                      + (m02 * m02 + m12 * m12) / 8. / (m02 - m12) * log(m02 / m12);  
            }             
        } else
            throw std::runtime_error("PVfunctions::B22() is undefined.");  
    } else {
        if(m0 != 0. && m1 != 0.){
            if(fabs(m02 - m12) < LEPS){
              B22 = - (6. * m02 - p2) / 18. + A0(mu,m0) /6. 
                    + (p2 - 4. * m02) / 12. * B0(mu,p2,m0,m1);  
            } else {
              B22 = - (3. * (m02 + m12) - p2) / 18. 
                    + (DeltaM2 + p2) / 12. / p2 * A0(mu,m0) 
                    - (DeltaM2 - p2) / 12. / p2 * A0(mu,m1) 
                    + (- m02 - m12 + p2 / 2. 
                       + (DeltaM2 * DeltaM2) / 2. / p2) * B0(mu,p2,m0,m1) / 6.
                      ;                
            }
        } else
            throw std::runtime_error("PVfunctions::B22() is undefined.");  
    }
    return B22;
}

complex PVfunctions::Bf(const double mu, const double p2, 
                        const double m0, const double m1) const {   
    if ( mu<=0.0 || p2<0.0 || m0<0.0 || m1<0.0 ) {
        throw std::runtime_error("Invalid argument for PVfunctions::Bf()"); 
    }
    complex Bf(0.0, 0.0, false);
    Bf = 2.0*(B21(mu,p2,m0,m1) + B1(mu,p2,m0,m1));
    return Bf;
}

complex PVfunctions::B0p(const double mu, const double p2, 
                         const double m0, const double m1) const {   
    if ( mu<=0.0 || p2<0.0 || m0<0.0 || m1<0.0 )
        throw std::runtime_error("Invalid argument for PVfunctions::B0p()"); 
    
    double mu2=mu*mu, m02=m0*m0, m12=m1*m1;
    complex B0p(0.0, 0.0, false);
        
    if (p2==0.0) {
        if ( m0!=0.0 && m1!=0.0 ) {
            if ( fabs(m02 - m12) < LEPS )
                B0p = 1.0/6.0/m02;
            else {
                double DeltaM2 = m02 - m12;
                B0p = (m02 + m12)/2.0/pow(DeltaM2,2.0)
                        + m02*m12/pow(DeltaM2,3.0)*log(m12/m02);
            } 
        } else if ( m0!=0.0 && m1==0.0 )
            B0p = 1.0/2.0/m02;
        else if ( m0==0.0 && m1!=0.0 )
            B0p = 1.0/2.0/m12;
        else
            throw std::runtime_error("PVfunctions::B0p() is undefined."); 
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
            throw std::runtime_error("PVfunctions::B0p() is undefined."); 
        }
    }
    return B0p;
}

complex PVfunctions::B1p(const double mu, const double p2, 
                         const double m0, const double m1) const {   
    if ( mu<=0.0 || p2<0.0 || m0<0.0 || m1<0.0 ) {
        throw std::runtime_error("Invalid argument for PVfunctions::B1p()"); 
    }
    complex B1p(0.0, 0.0, false);
    
    if (p2==0.0) {
        throw std::runtime_error("PVfunctions::B1p() is undefined.");                      
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
        throw std::runtime_error("Invalid argument for PVfunctions::B21p()"); 
    }
    double m02=m0*m0, m12=m1*m1;
    double p4 = p2*p2, p6=p2*p2*p2;
    double DeltaM2 = m02 - m12; 
    complex B21p(0.0, 0.0, false);

    if (p2==0.0) {
        throw std::runtime_error("PVfunctions::B21p() is undefined.");            
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
//        throw std::runtime_error("Invalid argument for PVfunctions::B22p()"); 
//    }
//    double mu2=mu*mu, m02=m0*m0, m12=m1*m1;
//    complex B22p(0.0, 0.0, false);
//
//    return B22p;
//}

complex PVfunctions::Bfp(const double mu, const double p2, 
                         const double m0, const double m1) const {   
    if ( mu<=0.0 || p2<0.0 || m0<0.0 || m1<0.0 ) {
        throw std::runtime_error("Invalid argument for PVfunctions::Bfp()"); 
    }
    complex Bfp(0.0, 0.0, false);
    Bfp = 2.0*(B21p(mu,p2,m0,m1) + B1p(mu,p2,m0,m1));
    return Bfp;
}

complex PVfunctions::C0(const double p2, 
                        const double m0, const double m1, const double m2) const {
    if ( p2<0.0 || m0<0.0 || m1<0.0 || m2<0.0 ) {
        throw std::runtime_error("Invalid argument for PVfunctions::C0()"); 
    }

    complex C0(0.0, 0.0, false);   
    if (p2==0.0) {
        throw std::runtime_error("\nPVfunctions::C0() is undefined-1.\n"); 
    } else {
        if (fabs(m0 - m2) < LEPS && fabs(m0 - m1) > LEPS) {
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

            if ( x0==x1 || x0==x2 || x0==x3)
                throw std::runtime_error("\nPVfunctions::C0() is undefined-2.\n"); 

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
        } else if (m0!=0.0 && m2!=0.0 && fabs(m0 - m2) > LEPS && m1==0.0) { 
            double m02 = m0*m0;
            double m22 = m2*m2;
            double epsilon = 1.0e-12;
            double tmp_real = pow((m02+m22-p2),2.0) - 4.0*m02*m22;
            gsl_complex tmp = gsl_complex_rect(tmp_real, epsilon);
            tmp = gsl_complex_sqrt(tmp);
            complex tmp_complex(GSL_REAL(tmp), GSL_IMAG(tmp), false);
            complex x1 = (p2 - m02 + m22 + tmp_complex)/2.0/p2;
            complex x2 = (p2 - m02 + m22 - tmp_complex)/2.0/p2;            

            if ( x1==0.0 || x1==1.0 || x2==0.0 || x2==1.0 )
                throw std::runtime_error("\nPVfunctions::C0() is undefined-3.\n"); 
            
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
        } else if (m0 == 0. && m1 != 0. && m2 != 0.) { 
            complex arg[2];
            arg[0] = 1.-m2*m2/m1/m1;
            arg[1] = 1.-(-p2+m2*m2)/m1/m1;
            complex Li2[2];
            for (int i=0; i<2; i++) {
                gsl_sf_result re, im;
                gsl_sf_complex_dilog_xy_e(arg[i].real(), arg[i].imag(), &re, &im);
                Li2[i].real() = re.val;
                Li2[i].imag() = im.val;
            }
            C0 = 1./(-p2)*(Li2[0]-Li2[1]);
        } else if(m0 != 0 && m1 != 0 && m2 != 0 && fabs(m0-m1) > LEPS && fabs(m1-m2) > LEPS){
            double x0 = 1.-(m0*m0-m1*m1)/p2;
            double x1 = -(-p2+m0*m0-m2*m2-sqrt(fabs((m0*m0+m2*m2-p2)*(m0*m0+m2*m2-p2)-4.*m0*m0*m2*m2)))/p2/2.;
            double x2 = -(-p2+m0*m0-m2*m2+sqrt(fabs((m0*m0+m2*m2-p2)*(m0*m0+m2*m2-p2)-4.*m0*m0*m2*m2)))/p2/2.;
            double x3 = m2*m2/(m2*m2-m1*m1);   
            
            complex arg[6];
            arg[0] = (x0-1.)/(x0-x1);
            arg[1] = x0/(x0-x1);
            arg[2] = (x0-1.)/(x0-x2);
            arg[3] = x0/(x0-x2);
            arg[4] = (x0-1.)/(x0-x3);
            arg[5] = x0/(x0-x3);
            
            complex Li2[6];
            for (int i=0; i<2; i++) {
                gsl_sf_result re, im;
                gsl_sf_complex_dilog_xy_e(arg[i].real(), arg[i].imag(), &re, &im);
                Li2[i].real() = re.val;
                Li2[i].imag() = im.val;
            }
            
            C0 = -1./p2*(Li2[0]-Li2[1]+Li2[2]-Li2[3]-Li2[4]+Li2[5]);
            
        } else {
            throw std::runtime_error("\nPVfunctions::C0() is undefined-4.\n");             
        }
        
        
    }
    return C0;
}

complex PVfunctions::C11(const double mu,const double p12, const double p22,const double s, 
                        const double m1, const double m2, const double m3) const{
    
    complex C02 = B0(mu,p12,m1,m2);
    complex C01 = B0(mu,s,m1,m3);
    double f2c = s+p12+m2*m2-m3*m3;
    
    return (-2./s*(C02-C01+f2c*C0(s,m1,m2,m3)));
    
}

complex PVfunctions::C12(const double mu,const double p12, const double p22,const double s, 
                        const double m1, const double m2, const double m3) const{
    
    complex C00 = B0(mu,p22,m2,m3);
    complex C01 = B0(mu,s,m1,m3);
    double f1c = -p12+m1*m1-m2*m2;
    
    return (-2./s*(C01-C00+f1c*C0(s,m1,m2,m3)));
    
}

double PVfunctions::F(const double m0, const double m1) const {
    double m12 = m1 * m1;
    double m02 = m0 * m0;
    double F;
    
    if ( m0<=0.0 || m1<0.0 )
        throw std::runtime_error("Invalid argument for PVfunctions::F()\n"); 
    
    if(m0 == 0. && m1 != 0.) {
        F=0.5 * m12;
    } else if(m0 != 0. && m1 == 0.){
        F=0.5 * m02;
    } else if((m0 == 0. && m1 == 0.) || (fabs(m0-m1) < LEPS)){
        F=0.;
    } else if (m0 != 0 && m1 != 0){
        F=0.5 * (m02 + m12) - (m02 * m12) / (m02 - m12) * log(m02 / m12);
    }
    return (F);
}


complex PVfunctions::D0(const double p12, const double p22, const double p32,
                        const double p42, const double s, const double t, 
                        const double m1, const double m2, const double m3, 
                        const double m4) const{

    complex D0;

    if(p12==0. && p22==0. && p32==0. && p42==0.) {
        if(m2 == 0. && fabs(m1-m3) < LEPS && m1 != 0.){

            double d4 =1.+(4.*m1*m1*t*(t+m4*m4-m1*m1))/(s*(t+m4*m4)*(t+m4*m4));
            double x1 = 0.5*(1.-sqrt(1.+4.*m1*m1/s));
            double x2 = 0.5*(1.+sqrt(1.+4.*m1*m1/s));
            double x3 = m4*m4/(m4*m4-m1*m1);
            double x4 =(t+m4*m4)/(t+m4*m4-m1*m1);

            if(m4 != 0.){
                
                double xbar1 = x4/2*(1.-sqrt(d4));
                double xbar2 = x4/2*(1.+sqrt(d4));
                
                
                complex arg[16];
                arg[0] = xbar1/(xbar1-x1);
                arg[1] = (xbar1-1.)/(xbar1-x1);
                arg[2] = xbar2/(xbar2-x1);
                arg[3] = (xbar2-1.)/(xbar2-x1);
                arg[4] = xbar1/(xbar1-x2);
                arg[5] = (xbar1-1.)/(xbar1-x2);
                arg[6] = xbar2/(xbar2-x2);
                arg[7] = (xbar2-1.)/(xbar2-x2);
                arg[8] = xbar1/(xbar1-x3);
                arg[9] = (xbar1-1.)/(xbar1-x3);
                arg[10] = xbar2/(xbar2-x3);
                arg[11] = (xbar2-1.)/(xbar2-x3);
                arg[12] = xbar1/(xbar1-x4);
                arg[13] = (xbar1-1.)/(xbar1-x4);
                arg[14] = xbar2/(xbar2-x4);
                arg[15] = (xbar2-1.)/(xbar2-x4);
                complex Li2[16];
                for (int i=0; i<16; i++) {
                    gsl_sf_result re, im;
                    gsl_sf_complex_dilog_xy_e(arg[i].real(), arg[i].imag(), &re, &im);
                    Li2[i].real() = re.val;
                    Li2[i].imag() = im.val;
                }
                
                double A = 1./(s*(t+m4*m4)*sqrt(d4));
                
                D0=A*(Li2[0]-Li2[1]-Li2[2]+Li2[3]+Li2[4]-Li2[5]-Li2[6]
                        +Li2[7]-Li2[8]+Li2[9]+Li2[10]-Li2[11]+Li2[12]-Li2[13]
                        -Li2[14]+Li2[15]);
                
            } else if (m4 == 0.){
              
                double xbar1 = 0.5*x4*(1.-sqrt(d4));
                double xbar2 = 0.5*x4*(1.+sqrt(d4));
                complex arg[4];
                arg[0] = xbar1/(xbar1-x1);
                arg[1] = xbar1/(xbar1-x2);
                arg[2] = xbar2/(xbar2-x1);
                arg[3] = xbar2/(xbar2-x2);
                complex Li2[4];
                for (int i=0; i<4; i++) {
                    gsl_sf_result re, im;
                    gsl_sf_complex_dilog_xy_e(arg[i].real(), arg[i].imag(), &re, &im);
                    Li2[i].real() = re.val;
                    Li2[i].imag() = im.val;
                }
                double A = 2./(s*t*sqrt(d4));
                
                D0=A*(Li2[0]-Li2[1]+Li2[2]-Li2[3]);                
            }
        
        } else {
            throw std::runtime_error("Invalid argument m2 != 0. for PVfunctions::D0()!"); 
        }
        
    }else {
        throw std::runtime_error("Invalid argument for PVfunctions::D0(). The initial and final momenta must be 0."); 
    }
    
    return (D0);
    
}
    
    

complex PVfunctions::D11(const double p12, const double p22,
                         const double p32,const double p42, const double s,
                         const double m1, const double m2, const double m3, 
                         const double m4,const double cos_theta) const{
    
    complex D00 = C0(+s/2.*(1.+cos_theta),m2,m3,m4);
    complex D01 = C0(s,m3,m4,m1);
    complex D02 = C0(+s/2.*(1.+cos_theta),m4,m1,m2);
    complex D03 = C0(s,m1,m2,m3);
    double f1d = m1*m1-m2*m2-p12;
    double f2d = m2*m2-m3*m3+p12+s;
    double f3d = m3*m3-m4*m4-p42-s;
    double t = -s/2.*(1.-cos_theta);
    
    complex a =D01-D00+f1d*D0(p12,p22,p32,p42,s,-t,m1,m2,m3,m4);
    complex b = D02-D01+f2d*D0(p12,p22,p32,p42,s,-t,m1,m2,m3,m4);       
    complex c = D03-D02+f3d*D0(p12,p22,p32,p42,s,-t,m1,m2,m3,m4);
    
    double X11=0.25*(1.-cos_theta)*(1.-cos_theta);
    double X12=-.25*(1.-cos_theta*cos_theta);
    double X13=0.5*(1.-cos_theta);
    double g = -1./s/(1.-cos_theta*cos_theta);
            
    return (g*(X11*a+X12*b+X13*c));
}


complex PVfunctions::D12(const double p12, const double p22,
                         const double p32,const double p42, const double s,
                         const double m1, const double m2, const double m3, 
                         const double m4,const double cos_theta) const{
    
    complex D00 = C0(+s/2.*(1.+cos_theta),m2,m3,m4);
    complex D01 = C0(s,m3,m4,m1);
    complex D02 = C0(+s/2.*(1.+cos_theta),m4,m1,m2);
    complex D03 = C0(s,m1,m2,m3);
    double f1d = m1*m1-m2*m2-p12;
    double f2d = m2*m2-m3*m3+p12+s;
    double f3d = m3*m3-m4*m4-p42-s;
    double t = -s/2.*(1.-cos_theta);
    
    complex a =D01-D00+f1d*D0(p12,p22,p32,p42,s,-t,m1,m2,m3,m4);
    complex b = D02-D01+f2d*D0(p12,p22,p32,p42,s,-t,m1,m2,m3,m4);       
    complex c = D03-D02+f3d*D0(p12,p22,p32,p42,s,-t,m1,m2,m3,m4);
    
    double X22=0.25*(1.+cos_theta)*(1.+cos_theta);
    double X12=-0.25*(1.-cos_theta*cos_theta);
    double X23=-0.5*(1.+cos_theta);
    double g = -1./s/(1.-cos_theta*cos_theta);
            
    return (g*(X12*a+X22*b+X23*c));
}


complex PVfunctions::D13(const double p12, const double p22,
                         const double p32,const double p42, const double s,
                         const double m1, const double m2, const double m3, 
                         const double m4,const double cos_theta) const{
    
    complex D00 = C0(+s/2.*(1.+cos_theta),m2,m3,m4);
    complex D01 = C0(s,m3,m4,m1);
    complex D02 = C0(+s/2.*(1.+cos_theta),m4,m1,m2);
    complex D03 = C0(s,m1,m2,m3);
    double f1d = m1*m1-m2*m2-p12;
    double f2d = m2*m2-m3*m3+p12+s;
    double f3d = m3*m3-m4*m4-p42-s;
    double t = -s/2.*(1.-cos_theta);
    
    complex a =D01-D00+f1d*D0(p12,p22,p32,p42,s,-t,m1,m2,m3,m4);
    complex b = D02-D01+f2d*D0(p12,p22,p32,p42,s,-t,m1,m2,m3,m4);       
    complex c = D03-D02+f3d*D0(p12,p22,p32,p42,s,-t,m1,m2,m3,m4);
    
    double X13=0.5*(1.-cos_theta);
    double X23=-0.5*(1.+cos_theta);
    double g = -1./s/(1.-cos_theta*cos_theta);
            
    return (g*(X13*a+X23*b+c));
}

complex PVfunctions::D27(const double p12, const double p22,
                         const double p32,const double p42, const double s,
                         const double m1, const double m2, const double m3, 
                         const double m4,const double cos_theta) const{
    
    complex D00 = C0(+s/2.*(1.+cos_theta),m2,m3,m4);
    complex D01 = C0(s,m3,m4,m1);
    complex D02 = C0(+s/2.*(1.+cos_theta),m4,m1,m2);
    complex D03 = C0(s,m1,m2,m3);
    double f1d = m1*m1-m2*m2-p12;
    double f2d = m2*m2-m3*m3+p12+s;
    double f3d = m3*m3-m4*m4-p42-s;
    double t = -s/2.*(1.-cos_theta);
    
//    complex a =D01-D00+f1d*D0(p12,p22,p32,p42,s,-t,m1,m2,m3,m4);
//    complex b = D02-D01+f2d*D0(p12,p22,p32,p42,s,-t,m1,m2,m3,m4);       
//    complex c = D03-D02+f3d*D0(p12,p22,p32,p42,s,-t,m1,m2,m3,m4);
//    
//    double X13=0.5*(1.-cos_theta);
//    double X23=-0.5*(1.+cos_theta);
//    double g = -1./s/(1.-cos_theta*cos_theta);
            
    return (-m1*m1*D0(p12,p22,p32,p42,s,-t,m1,m2,m3,m4)
            +0.5*(D00-f1d*D11(p12,p22,p32,p42,s,m1,m2,m3,m4,cos_theta)
            -f2d*D12(p12,p22,p32,p42,s,m1,m2,m3,m4,cos_theta)
            -f3d*D13(p12,p22,p32,p42,s,m1,m2,m3,m4,cos_theta)));
}


complex PVfunctions::D21(const double mu,const double p12, const double p22,
                         const double p32,const double p42, const double s,
                         const double m1, const double m2, const double m3, 
                         const double m4,const double cos_theta) const{
    
    complex D00 = C0(+s/2.*(1.+cos_theta),m2,m3,m4);
    complex D01 = C0(s,m3,m4,m1);
    complex D02 = C0(+s/2.*(1.+cos_theta),m4,m1,m2);
    complex D03 = C0(s,m1,m2,m3);
    double f1d = m1*m1-m2*m2-p12;
    double f2d = m2*m2-m3*m3+p12+s;
    double f3d = m3*m3-m4*m4-p42-s;
    double t = -s/2.*(1.-cos_theta);
    complex D11_0= C11(mu,p22,p32,-t,m2,m3,m4);
    complex D11_1=-C0(s,m3,m4,m1)-C12(mu,p32,p42,s,m3,m4,m1);
    complex D11_2=-C0(-t,m4,m1,m2)-C11(mu,p42,p12,-t,m4,m1,m2)
                  +C12(mu,p42,p12,-t,m4,m1,m2);
    complex D11_3=C11(mu,p12,p22,s,m1,m2,m3);
    complex D12_0=C12(mu,p22,p32,-t,m2,m3,m4);
    complex D12_1=C11(mu,p32,p42,s,m3,m4,m1)-C12(mu,p32,p42,s,m3,m4,m1);
    complex D12_2=-C0(-t,m4,m1,m2)-C11(mu,p42,p12,s,m4,m1,m2);
    complex D12_3=C12(mu,p12,p22,s,m1,m2,m3);
    
    double X11=0.25*(1.-cos_theta)*(1.-cos_theta);
    double X12=-.25*(1.-cos_theta*cos_theta);
    double X13=0.5*(1.-cos_theta);
    double g = -1./s/(1.-cos_theta*cos_theta);
    complex a = 0.5*(D11_1+D00+f1d*D11(p12,p22,p32,p42,s,m1,m2,m3,m4,cos_theta)
                -2.*D27(p12,p22,p32,p42,s,m1,m2,m3,m4,cos_theta));
    complex b = 0.5*(D11_2-D11_1+f2d*D11(p12,p22,p32,p42,s,m1,m2,m3,m4,cos_theta));
    complex c = 0.5*(D11_3-D11_2+f3d*D11(p12,p22,p32,p42,s,m1,m2,m3,m4,cos_theta));
            
    return (g*(X11*a+X12*b+X13*c));
}

complex PVfunctions::D24(const double mu,const double p12, const double p22,
                         const double p32,const double p42, const double s,
                         const double m1, const double m2, const double m3, 
                         const double m4,const double cos_theta) const{
    
    complex D00 = C0(+s/2.*(1.+cos_theta),m2,m3,m4);
    complex D01 = C0(s,m3,m4,m1);
    complex D02 = C0(+s/2.*(1.+cos_theta),m4,m1,m2);
    complex D03 = C0(s,m1,m2,m3);
    double f1d = m1*m1-m2*m2-p12;
    double f2d = m2*m2-m3*m3+p12+s;
    double f3d = m3*m3-m4*m4-p42-s;
    double t = -s/2.*(1.-cos_theta);
    complex D11_0= C11(mu,p22,p32,-t,m2,m3,m4);
    complex D11_1=-C0(s,m3,m4,m1)-C12(mu,p32,p42,s,m3,m4,m1);
    complex D11_2=-C0(-t,m4,m1,m2)-C11(mu,p42,p12,-t,m4,m1,m2)
                  +C12(mu,p42,p12,-t,m4,m1,m2);
    complex D11_3=C11(mu,p12,p22,s,m1,m2,m3);
    complex D12_0=C12(mu,p22,p32,-t,m2,m3,m4);
    complex D12_1=C11(mu,p32,p42,s,m3,m4,m1)-C12(mu,p32,p42,s,m3,m4,m1);
    complex D12_2=-C0(-t,m4,m1,m2)-C11(mu,p42,p12,s,m4,m1,m2);
    complex D12_3=C12(mu,p12,p22,s,m1,m2,m3);
    
    double X22=0.25*(1.+cos_theta)*(1.+cos_theta);
    double X12=-0.25*(1.-cos_theta*cos_theta);
    double X23=-0.5*(1.+cos_theta);
    double g = -1./s/(1.-cos_theta*cos_theta);
            
    
    complex a = 0.5*(D11_1+D00+f1d*D11(p12,p22,p32,p42,s,m1,m2,m3,m4,cos_theta)
                -2.*D27(p12,p22,p32,p42,s,m1,m2,m3,m4,cos_theta));
    complex b = 0.5*(D11_2-D11_1+f2d*D11(p12,p22,p32,p42,s,m1,m2,m3,m4,cos_theta));
    complex c = 0.5*(D11_3-D11_2+f3d*D11(p12,p22,p32,p42,s,m1,m2,m3,m4,cos_theta));
            
    return (g*(X12*a+X22*b+X23*c));
}

complex PVfunctions::D25(const double mu,const double p12, const double p22,
                         const double p32,const double p42, const double s,
                         const double m1, const double m2, const double m3, 
                         const double m4,const double cos_theta) const{
    
    complex D00 = C0(+s/2.*(1.+cos_theta),m2,m3,m4);
    complex D01 = C0(s,m3,m4,m1);
    complex D02 = C0(+s/2.*(1.+cos_theta),m4,m1,m2);
    complex D03 = C0(s,m1,m2,m3);
    double f1d = m1*m1-m2*m2-p12;
    double f2d = m2*m2-m3*m3+p12+s;
    double f3d = m3*m3-m4*m4-p42-s;
    double t = -s/2.*(1.-cos_theta);
    complex D11_0= C11(mu,p22,p32,-t,m2,m3,m4);
    complex D11_1=-C0(s,m3,m4,m1)-C12(mu,p32,p42,s,m3,m4,m1);
    complex D11_2=-C0(-t,m4,m1,m2)-C11(mu,p42,p12,-t,m4,m1,m2)
                  +C12(mu,p42,p12,-t,m4,m1,m2);
    complex D11_3=C11(mu,p12,p22,s,m1,m2,m3);
    complex D12_0=C12(mu,p22,p32,-t,m2,m3,m4);
    complex D12_1=C11(mu,p32,p42,s,m3,m4,m1)-C12(mu,p32,p42,s,m3,m4,m1);
    complex D12_2=-C0(-t,m4,m1,m2)-C11(mu,p42,p12,s,m4,m1,m2);
    complex D12_3=C12(mu,p12,p22,s,m1,m2,m3);
    
    
    complex a = 0.5*(D11_1+D00+f1d*D11(p12,p22,p32,p42,s,m1,m2,m3,m4,cos_theta)
                -2.*D27(p12,p22,p32,p42,s,m1,m2,m3,m4,cos_theta));
    complex b = 0.5*(D11_2-D11_1+f2d*D11(p12,p22,p32,p42,s,m1,m2,m3,m4,cos_theta));
    complex c = 0.5*(D11_3-D11_2+f3d*D11(p12,p22,p32,p42,s,m1,m2,m3,m4,cos_theta));
    
    double X13=0.5*(1.-cos_theta);
    double X23=-0.5*(1.+cos_theta);
    double g = -1./s/(1.-cos_theta*cos_theta);
            
    return (g*(X13*a+X23*b+c));
}

complex PVfunctions::D22(const double mu,const double p12, const double p22,
                         const double p32,const double p42, const double s,
                         const double m1, const double m2, const double m3, 
                         const double m4,const double cos_theta) const{
    
    complex D00 = C0(+s/2.*(1.+cos_theta),m2,m3,m4);
    complex D01 = C0(s,m3,m4,m1);
    complex D02 = C0(+s/2.*(1.+cos_theta),m4,m1,m2);
    complex D03 = C0(s,m1,m2,m3);
    double f1d = m1*m1-m2*m2-p12;
    double f2d = m2*m2-m3*m3+p12+s;
    double f3d = m3*m3-m4*m4-p42-s;
    double t = -s/2.*(1.-cos_theta);
    complex D11_0= C11(mu,p22,p32,-t,m2,m3,m4);
    complex D11_1=-C0(s,m3,m4,m1)-C12(mu,p32,p42,s,m3,m4,m1);
    complex D11_2=-C0(-t,m4,m1,m2)-C11(mu,p42,p12,-t,m4,m1,m2)
                  +C12(mu,p42,p12,-t,m4,m1,m2);
    complex D11_3=C11(mu,p12,p22,s,m1,m2,m3);
    complex D12_0=C12(mu,p22,p32,-t,m2,m3,m4);
    complex D12_1=C11(mu,p32,p42,s,m3,m4,m1)-C12(mu,p32,p42,s,m3,m4,m1);
    complex D12_2=-C0(-t,m4,m1,m2)-C11(mu,p42,p12,s,m4,m1,m2);
    complex D12_3=C12(mu,p12,p22,s,m1,m2,m3);
    
    
    complex d = 0.5*(D11_1-D11_0+f1d*D12(p12,p22,p32,p42,s,m1,m2,m3,m4,cos_theta));
    complex e = 0.5*(D12_2-D11_1+f2d*D12(p12,p22,p32,p42,s,m1,m2,m3,m4,cos_theta)-
            2.*D27(p12,p22,p32,p42,s,m1,m2,m3,m4,cos_theta));
    complex f = 0.5*(D12_3-D12_2+f3d*D12(p12,p22,p32,p42,s,m1,m2,m3,m4,cos_theta));
    
    double X22=0.25*(1.+cos_theta)*(1.+cos_theta);
    double X12=-0.25*(1.-cos_theta*cos_theta);
    double X23=-0.5*(1.+cos_theta);
    double g = -1./s/(1.-cos_theta*cos_theta);
            
    return (g*(X12*d+X22*e+X23*f));
}

complex PVfunctions::D26(const double mu,const double p12, const double p22,
                         const double p32,const double p42, const double s,
                         const double m1, const double m2, const double m3, 
                         const double m4,const double cos_theta) const{
    
    complex D00 = C0(+s/2.*(1.+cos_theta),m2,m3,m4);
    complex D01 = C0(s,m3,m4,m1);
    complex D02 = C0(+s/2.*(1.+cos_theta),m4,m1,m2);
    complex D03 = C0(s,m1,m2,m3);
    double f1d = m1*m1-m2*m2-p12;
    double f2d = m2*m2-m3*m3+p12+s;
    double f3d = m3*m3-m4*m4-p42-s;
    double t = -s/2.*(1.-cos_theta);
    complex D11_0= C11(mu,p22,p32,-t,m2,m3,m4);
    complex D11_1=-C0(s,m3,m4,m1)-C12(mu,p32,p42,s,m3,m4,m1);
    complex D11_2=-C0(-t,m4,m1,m2)-C11(mu,p42,p12,-t,m4,m1,m2)
                  +C12(mu,p42,p12,-t,m4,m1,m2);
    complex D11_3=C11(mu,p12,p22,s,m1,m2,m3);
    complex D12_0=C12(mu,p22,p32,-t,m2,m3,m4);
    complex D12_1=C11(mu,p32,p42,s,m3,m4,m1)-C12(mu,p32,p42,s,m3,m4,m1);
    complex D12_2=-C0(-t,m4,m1,m2)-C11(mu,p42,p12,s,m4,m1,m2);
    complex D12_3=C12(mu,p12,p22,s,m1,m2,m3);
    
    
    complex d = 0.5*(D11_1-D11_0+f1d*D12(p12,p22,p32,p42,s,m1,m2,m3,m4,cos_theta));
    complex e = 0.5*(D12_2-D11_1+f2d*D12(p12,p22,p32,p42,s,m1,m2,m3,m4,cos_theta)-
            2.*D27(p12,p22,p32,p42,s,m1,m2,m3,m4,cos_theta));
    complex f = 0.5*(D12_3-D12_2+f3d*D12(p12,p22,p32,p42,s,m1,m2,m3,m4,cos_theta));
    
    double X13=0.5*(1.-cos_theta);
    double X23=-0.5*(1.+cos_theta);
    double g = -1./s/(1.-cos_theta*cos_theta);
            
    return (g*(X13*d+X23*e+f));
}

complex PVfunctions::D23(const double mu,const double p12, const double p22,
                         const double p32,const double p42, const double s,
                         const double m1, const double m2, const double m3, 
                         const double m4,const double cos_theta) const{
    
    complex D00 = C0(+s/2.*(1.+cos_theta),m2,m3,m4);
    complex D01 = C0(s,m3,m4,m1);
    complex D02 = C0(+s/2.*(1.+cos_theta),m4,m1,m2);
    complex D03 = C0(s,m1,m2,m3);
    double f1d = m1*m1-m2*m2-p12;
    double f2d = m2*m2-m3*m3+p12+s;
    double f3d = m3*m3-m4*m4-p42-s;
    double t = -s/2.*(1.-cos_theta);
    complex D11_0= C11(mu,p22,p32,-t,m2,m3,m4);
    complex D11_1=-C0(s,m3,m4,m1)-C12(mu,p32,p42,s,m3,m4,m1);
    complex D11_2=-C0(-t,m4,m1,m2)-C11(mu,p42,p12,-t,m4,m1,m2)
                  +C12(mu,p42,p12,-t,m4,m1,m2);
    complex D11_3=C11(mu,p12,p22,s,m1,m2,m3);
    complex D12_0=C12(mu,p22,p32,-t,m2,m3,m4);
    complex D12_1=C11(mu,p32,p42,s,m3,m4,m1)-C12(mu,p32,p42,s,m3,m4,m1);
    complex D12_2=-C0(-t,m4,m1,m2)-C11(mu,p42,p12,s,m4,m1,m2);
    complex D12_3=C12(mu,p12,p22,s,m1,m2,m3);
    
    
    complex h = 0.5*(D12_1-D12_0+f1d*D13(p12,p22,p32,p42,s,m1,m2,m3,m4,cos_theta));
    complex i = 0.5*(D12_2-D12_1+f2d*D13(p12,p22,p32,p42,s,m1,m2,m3,m4,cos_theta));
    complex l = 0.5*(-D12_2+f3d*D13(p12,p22,p32,p42,s,m1,m2,m3,m4,cos_theta)-
                2.*D27(p12,p22,p32,p42,s,m1,m2,m3,m4,cos_theta));
    
    double X13=0.5*(1.-cos_theta);
    double X23=-0.5*(1.+cos_theta);
    double g = -1./s/(1.-cos_theta*cos_theta);
        
    return (g*(X13*h+X23*i+l));
} 


