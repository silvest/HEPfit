/* 
 * Copyright (C) 2012 HEPfit Collaboration
 * All rights reserved.
 *
 * For the licensing terms see doc/COPYING.
 */

#include <iostream>
#include <stdexcept>
#include <cmath>
#include <algorithm>
#include <gsl/gsl_complex.h>
#include <gsl/gsl_sf.h>
#include "PVfunctions.h"


PVfunctions::PVfunctions(const bool bExtraMinusSign)
{
    if (bExtraMinusSign) ExtraMinusSign = -1.0;
    else ExtraMinusSign = 1.0;
}

double PVfunctions::A0(const double mu2, const double m2) const
{
#ifdef USE_LOOPTOOLS
    return ( ExtraMinusSign * myLT.PV_A0(mu2, m2) );
#else    
    if ( mu2<=0.0 || m2<0.0 )
        throw std::runtime_error("PVfunctions::A0(): Invalid argument!");

    const double Tolerance = 1.0e-10;
    const bool m2zero = (m2 <= mu2*Tolerance);
    
    if ( m2zero )
        return ( 0.0 );
    else
        return ( ExtraMinusSign * m2*(-log(m2/mu2)+1.0) );
#endif
}

gslpp::complex PVfunctions::B0(const double mu2, const double p2,
                        const double m02, const double m12) const
{   
#ifdef USE_LOOPTOOLS
    return myLT.PV_B0(mu2, p2, m02, m12);
#else    
    if ( mu2<=0.0 || p2<0.0 || m02<0.0 || m12<0.0 )
        throw std::runtime_error("PVfunctions::B0(): Invalid argument!");

    const double Tolerance = 1.0e-8;
    const double maxM2 = std::max(p2, std::max(m02, m12));
    const bool p2zero = (p2 <= maxM2*Tolerance);
    const bool m02zero = (m02 <= maxM2*Tolerance);
    const bool m12zero = (m12 <= maxM2*Tolerance);
    const bool m02_eq_m12 = (fabs(m02 - m12) <= (m02 + m12)*Tolerance);

    gslpp::complex B0(0.0, 0.0, false);
    if ( p2zero ) {
        if ( !m02zero && !m12zero ) {
            if ( m02_eq_m12 )
                B0 = - log(m02/mu2);
            else
                B0 = - m02/(m02-m12)*log(m02/mu2) + m12/(m02-m12)*log(m12/mu2) + 1.0;
        } else if ( m02zero && !m12zero )
            B0 = - log(m12/mu2) + 1.0;
        else if ( !m02zero && m12zero )
            B0 = - log(m02/mu2) + 1.0;
        else if ( m02zero && m12zero )
            throw std::runtime_error("PVfunctions::B0(): IR divergent! (vanishes in DR)");
        else
            throw std::runtime_error("PVfunctions::B0(): Undefined!");
    } else {
        if ( !m02zero && !m12zero ) {
            double m0 = sqrt(m02), m1 = sqrt(m12);
            double Lambda = sqrt( fabs((p2-m02-m12)*(p2-m02-m12) - 4.0*m02*m12) );
            double R;
            if ( p2 > (m0-m1)*(m0-m1) && p2 < (m0+m1)*(m0+m1) ) {
                if ( p2-m02-m12 > 0.0 ) {
                    if ( p2-m02-m12 > Lambda*Tolerance )
                        R = - Lambda/p2*(atan(Lambda/(p2-m02-m12)) - M_PI);
                    else
                        R = - Lambda/p2*(M_PI/2.0 - M_PI);
                } else {
                    if ( - (p2-m02-m12) > Lambda*Tolerance )
                        R = - Lambda/p2*atan(Lambda/(p2-m02-m12));
                    else
                        R = - Lambda/p2*( -M_PI/2.0 );
                }
            } else
                R = Lambda/p2*log( fabs((p2-m02-m12+Lambda)/2.0/m0/m1) );
            B0 = - log(m0*m1/mu2) + (m02-m12)/2.0/p2*log(m12/m02) - R + 2.0;
            if ( p2 > (m0+m1)*(m0+m1) )
                B0 += M_PI*Lambda/p2*gslpp::complex::i();// imaginary part
        } else if ( (m02zero && !m12zero) || (!m02zero && m12zero) ) {
            double M2;
            if (!m02zero) M2 = m02;
            else M2 = m12;
            B0 = - log(M2/mu2) + 2.0;
            if ( p2 < M2 )
                B0 += - (1.0 - M2/p2)*log(1.0 - p2/M2);
            else if ( p2 > M2 ) {
                B0 += - (1.0 - M2/p2)*log(p2/M2 - 1.0);
                B0 += (1.0 - M2/p2)*M_PI*gslpp::complex::i();// imaginary part
            } else
                B0 += 0.0;
        } else if ( m02zero && m12zero ) {
            if ( p2 < 0.0 )
                B0 = - log(-p2/mu2) + 2.0;
            else {
                B0 = - log(p2/mu2) + 2.0;
                B0 += M_PI*gslpp::complex::i();// imaginary part
            }
        } else
            throw std::runtime_error("PVfunctions::B0(): Undefined!");
    }
    return B0;
#endif    
}

gslpp::complex PVfunctions::B1(const double mu2, const double p2,
                        const double m02, const double m12) const
{
#ifdef USE_LOOPTOOLS
    return myLT.PV_B1(mu2, p2, m02, m12);
#else
    if ( mu2<=0.0 || p2<0.0 || m02<0.0 || m12<0.0 )
        throw std::runtime_error("PVfunctions::B1(): Invalid argument!");

    const double Tolerance = 1.0e-8;
    const double maxM2 = std::max(p2, std::max(m02, m12));
    const bool p2zero = (p2 <= maxM2*Tolerance);
    const bool m02zero = (m02 <= maxM2*Tolerance);
    const bool m12zero = (m12 <= maxM2*Tolerance);
    const bool m02_eq_m12 = (fabs(m02 - m12) <= (m02 + m12)*Tolerance);

    gslpp::complex B1(0.0, 0.0, false);    
    double DeltaM2 = m02 - m12;
    if ( p2zero ) {
        if ( !m02zero && !m12zero ) {
            if ( m02_eq_m12 ) 
                B1.real() = 1.0/2.0*log(m12/mu2);
            else {
                double F0 = - log(m12/m02);
                double F1 = - 1.0 + m02/DeltaM2*F0;
                double F2 = - 1.0/2.0 + m02/DeltaM2*F1;
                B1.real() = 1.0/2.0*( log(m12/mu2) + F2 );
            }
        } else if ( m02zero && !m12zero )
            B1.real() = 1.0/2.0*log(m12/mu2) - 1.0/4.0;
        else if ( !m02zero && m12zero )
            B1.real() = 1.0/2.0*log(m02/mu2) - 1.0/4.0;
        else
            B1 = 0.0;
    } else
        B1 = -1.0/2.0/p2*(- ExtraMinusSign*A0(mu2,m02)
                          + ExtraMinusSign*A0(mu2,m12)
                          + (DeltaM2 + p2)*B0(mu2,p2,m02,m12));
    return B1;
#endif
}

gslpp::complex PVfunctions::B11(const double mu2, const double p2,
                         const double m02, const double m12) const
{   
#ifdef USE_LOOPTOOLS
    return myLT.PV_B11(mu2, p2, m02, m12);
#else
    if ( mu2<=0.0 || p2<0.0 || m02<0.0 || m12<0.0 )
        throw std::runtime_error("PVfunctions::B11(): Invalid argument!");

    const double Tolerance = 1.0e-8;
    const double maxM2 = std::max(p2, std::max(m02, m12));
    const bool p2zero = (p2 <= maxM2*Tolerance);
    const bool m02zero = (m02 <= maxM2*Tolerance);
    const bool m12zero = (m12 <= maxM2*Tolerance);
    const bool m02_eq_m12 = (fabs(m02 - m12) <= (m02 + m12)*Tolerance);

    gslpp::complex B11(0.0, 0.0, false);
    double DeltaM2 = m02 - m12;
    if ( p2zero ) {
        if ( !m02zero && !m12zero ) {
            if ( m02_eq_m12 )
                B11.real() = - 1.0/3.0*log(m12/mu2);
            else {
                double F0 = - log(m12/m02);
                double F1 = - 1.0 + m02/DeltaM2*F0;
                double F2 = - 1.0/2.0 + m02/DeltaM2*F1;
                double F3 = - 1.0/3.0 + m02/DeltaM2*F2;
                B11.real() = - 1.0/3.0*( log(m12/mu2) + F3 );
            }
        } else if ( m02zero && !m12zero )
            B11.real() = - 1.0/3.0*log(m12/mu2) + 1.0/9.0;
        else if ( !m02zero && m12zero )
            B11.real() = - 1.0/3.0*log(m02/mu2) + 1.0/9.0;
        else
            throw std::runtime_error("PVfunctions::B11(): Undefined!");
    } else {
        double Lambdabar2 = (p2-m02-m12)*(p2-m02-m12) - 4.0*m02*m12;
        B11 = - (3.0*(m02 + m12) - p2)/18.0/p2
              + (DeltaM2 + p2)/3.0/p2/p2*(-ExtraMinusSign)*A0(mu2,m02)
              - (DeltaM2 + 2.0*p2)/3.0/p2/p2*(-ExtraMinusSign)*A0(mu2,m12)
              + (Lambdabar2 + 3.0*p2*m02)/3.0/p2/p2*B0(mu2,p2,m02,m12);
    }
    return B11;
#endif
}

gslpp::complex PVfunctions::B00(const double mu2, const double p2,
                         const double m02, const double m12) const
{   
#ifdef USE_LOOPTOOLS
    return myLT.PV_B00(mu2, p2, m02, m12);
#else
    if ( mu2<=0.0 || p2<0.0 || m02<0.0 || m12<0.0 )
        throw std::runtime_error("PVfunctions::B00(): Invalid argument!");

    const double Tolerance = 1.0e-8;
    const double maxM2 = std::max(p2, std::max(m02, m12));
    const bool p2zero = (p2 <= maxM2*Tolerance);
    const bool m02zero = (m02 <= maxM2*Tolerance);
    const bool m12zero = (m12 <= maxM2*Tolerance);
    const bool m02_eq_m12 = (fabs(m02 - m12) <= (m02 + m12)*Tolerance);

    gslpp::complex B00(0.0, 0.0, false);
    if ( p2zero ) {
        if ( !m02zero && !m12zero ) {
            if( m02_eq_m12 )
                B00 = m02/2.0*(- log(m02/mu2) + 1.0);
            else
                B00 = 1.0/4.0*(m02 + m12)*(- log(sqrt(m02)*sqrt(m12)/mu2) + 3.0/2.0)
                      - (m02*m02 + m12*m12)/8.0/(m02 - m12)*log(m02/m12);
        } else if ( (!m02zero && m12zero) || (m02zero && !m12zero) ) {
            double M2;
            if ( !m02zero ) M2 = m02;
            else M2 = m12;
            B00 = M2/4.0*(- log(M2/mu2) + 3.0/2.0);
        } else
            B00 = 0.0;
    } else {
        if ( !m02zero && !m12zero ) {
            if ( m02_eq_m12 )
                B00 = (6.0*m02 - p2)/18.0 + ExtraMinusSign*A0(mu2,m02)/6.0
                       - (p2 - 4.0*m02)/12.0*B0(mu2,p2,m02,m12);
            else {
                double DeltaM2 = m02 - m12;
                double Lambdabar2 = (p2-m02-m12)*(p2-m02-m12) - 4.0*m02*m12;
                B00 = (3.0*(m02 + m12) - p2)/18.0
                       - (DeltaM2 + p2)/12.0/p2*(-ExtraMinusSign)*A0(mu2,m02)
                       + (DeltaM2 - p2)/12.0/p2*(-ExtraMinusSign)*A0(mu2,m12)
                       - Lambdabar2*B0(mu2,p2,m02,m12)/12.0/p2;
            }
        } else if ( (!m02zero && m12zero) || (m02zero && !m12zero) ) {
            double M2;
            if ( !m02zero ) M2 = m02;
            else M2 = m12;
            B00 = (3.0*M2 - p2)/18.0 - (M2 + p2)/12.0/p2*(-ExtraMinusSign)*A0(mu2,M2)
                  - (M2 - p2)*(M2 - p2)/12.0/p2*B0(mu2,p2,M2,0.0);
        } else
            B00 = - p2/18.0 - p2/12.0*B0(mu2,p2,0.0,0.0);
    }
    return B00;
#endif
}

gslpp::complex PVfunctions::Bf(const double mu2, const double p2,
                        const double m02, const double m12) const
{   
    if ( mu2<=0.0 || p2<0.0 || m02<0.0 || m12<0.0 )
        throw std::runtime_error("PVfunctions::Bf(): Invalid argument!");

    gslpp::complex Bf(0.0, 0.0, false);
    Bf = 2.0*(B11(mu2,p2,m02,m12) + B1(mu2,p2,m02,m12));
    return Bf;
}

gslpp::complex PVfunctions::B0p(const double muIR2, const double p2,
                         const double m02, const double m12) const
{   
#ifdef USE_LOOPTOOLS
    return myLT.PV_B0p(muIR2, p2, m02, m12);
#else
    if ( muIR2<=0.0 || p2<0.0 || m02<0.0 || m12<0.0 )
        throw std::runtime_error("PVfunctions::B0p(): Invalid argument!");

    const double Tolerance = 1.0e-8;
    const double maxM2 = std::max(p2, std::max(m02, m12));
    const bool p2zero = (p2 <= maxM2*Tolerance);
    const bool m02zero = (m02 <= maxM2*Tolerance);
    const bool m12zero = (m12 <= maxM2*Tolerance);
    const bool m02_eq_m12 = (fabs(m02 - m12) <= (m02 + m12)*Tolerance);

    gslpp::complex B0p(0.0, 0.0, false);        
    if ( p2zero ) {
        if ( !m02zero && !m12zero ) {
            if ( m02_eq_m12 )
                B0p = 1.0/6.0/m02;
            else {
                double DeltaM2 = m02 - m12;
                B0p = (m02 + m12)/2.0/pow(DeltaM2,2.0)
                        + m02*m12/pow(DeltaM2,3.0)*log(m12/m02);
            } 
        } else if ( !m02zero && m12zero )
            B0p = 1.0/2.0/m02;
        else if ( m02zero && !m12zero )
            B0p = 1.0/2.0/m12;
        else
            throw std::runtime_error("PVfunctions::B0p(): Undefined!");
    } else {
        if ( !m02zero && !m12zero ) {
            double m0 = sqrt(m02), m1 = sqrt(m12);
            double Lambda = sqrt( fabs((p2-m02-m12)*(p2-m02-m12) - 4.0*m02*m12) );
            double Rprime;
            if ( p2 > (m0-m1)*(m0-m1) && p2 < (m0+m1)*(m0+m1) ) {
                 if ( p2-m02-m12 > 0.0 ) {
                    if ( p2-m02-m12 > Lambda*Tolerance )
                        Rprime = ((p2 - m02 - m12)/Lambda + Lambda/p2)
                                  *(atan(Lambda/(p2-m02-m12)) - M_PI);
                    else
                        Rprime = ((p2 - m02 - m12)/Lambda + Lambda/p2)
                                  *(M_PI/2.0 - M_PI);
                } else {
                    if ( - (p2-m02-m12) > Lambda*Tolerance )
                        Rprime = ((p2 - m02 - m12)/Lambda + Lambda/p2)
                                  *atan(Lambda/(p2-m02-m12));
                    else
                        Rprime = ((p2 - m02 - m12)/Lambda + Lambda/p2)
                                  *( -M_PI/2.0 );
                }
            } else
                Rprime = ((p2 - m02 - m12)/Lambda - Lambda/p2)
                          *log( fabs((p2-m02-m12+Lambda)/2.0/m0/m1) );
            B0p = - (m02 - m12)/2.0/p2/p2*log(m12/m02) - (Rprime + 1.0)/p2;
            if ( p2 > (m0+m1)*(m0+m1) )
                B0p += M_PI/p2*((p2 - m02 - m12)/Lambda - Lambda/p2)
                      *gslpp::complex::i();// imaginary part
        } else if ( (m02zero && !m12zero) || (!m02zero && m12zero) ) {
            double M2;
            if ( !m02zero ) M2 = m02;
            else M2 = m12;
            if ( p2 < M2 )
                B0p = - M2/p2/p2*log(1.0 - p2/M2) - 1.0/p2;
            else if ( p2 > M2 ) {
                B0p = - M2/p2/p2*log(p2/M2 - 1.0) - 1.0/p2;
                B0p += M2/p2/p2*M_PI*gslpp::complex::i();// imaginary part        
            } else /* p2=M2 */
                B0p = 1.0/2.0/M2*(log(M2/muIR2) - 2.0);
        } else if ( m02zero && m12zero )
            B0p = - 1.0/p2;
        else
            throw std::runtime_error("PVfunctions::B0p(): Undefined!");
    }
    return B0p;
#endif
}

gslpp::complex PVfunctions::B1p(const double mu2, const double p2,
                         const double m02, const double m12) const
{   
#ifdef USE_LOOPTOOLS
    return myLT.PV_B1p(mu2, p2, m02, m12);
#else
    if ( mu2<=0.0 || p2<0.0 || m02<0.0 || m12<0.0 )
        throw std::runtime_error("PVfunctions::B1p(): Invalid argument!");

    const double Tolerance = 1.0e-8;
    const double maxM2 = std::max(p2, std::max(m02, m12));
    const bool p2zero = (p2 <= maxM2*Tolerance);

    gslpp::complex B1p(0.0, 0.0, false);    
    if ( p2zero )
        throw std::runtime_error("PVfunctions::B1p(): Undefined!");
    else {
        double DeltaM2 = m02 - m12;
        B1p = - ( 2.0*B1(mu2,p2,m02,m12) + B0(mu2,p2,m02,m12)
                 + (DeltaM2 + p2)*B0p(mu2,p2,m02,m12) )/2.0/p2;
    }
    return B1p;
#endif
}

gslpp::complex PVfunctions::B11p(const double mu2, const double p2,
                          const double m02, const double m12) const
{   
#ifdef USE_LOOPTOOLS
    return myLT.PV_B11p(mu2, p2, m02, m12);
#else
    if ( mu2<=0.0 || p2<0.0 || m02<0.0 || m12<0.0 )
        throw std::runtime_error("PVfunctions::B11p(): Invalid argument!");

    const double Tolerance = 1.0e-8;
    const double maxM2 = std::max(p2, std::max(m02, m12));
    const bool p2zero = (p2 <= maxM2*Tolerance);

    double p4 = p2*p2, p6=p2*p2*p2;
    double DeltaM2 = m02 - m12; 
    gslpp::complex B11p(0.0, 0.0, false);
    if ( p2zero )
        throw std::runtime_error("PVfunctions::B11p(): Undefined!");
    else {
        double Lambdabar2 = (p2-m02-m12)*(p2-m02-m12) - 4.0*m02*m12;
        B11p = (m02 + m12)/6.0/p4 - (2.0*DeltaM2 + p2)/3.0/p6*(-ExtraMinusSign)*A0(mu2,m02)
                + 2.0*(DeltaM2 + p2)/3.0/p6*(-ExtraMinusSign)*A0(mu2,m12)
                - (2.0*DeltaM2*DeltaM2 + p2*m02 - 2.0*p2*m12)/3.0/p6*B0(mu2,p2,m02,m12)
                + (Lambdabar2 + 3.0*p2*m02)/3.0/p4*B0p(mu2,p2,m02,m12);
    }
    return B11p;
#endif
}

gslpp::complex PVfunctions::B00p(const double mu2, const double p2,
                          const double m02, const double m12) const
{   
#ifdef USE_LOOPTOOLS
    return myLT.PV_B00p(mu2, p2, m02, m12);
#else
    if ( mu2<=0.0 || p2<0.0 || m02<0.0 || m12<0.0 )
        throw std::runtime_error("PVfunctions::B00p(): Invalid argument!");

    const double Tolerance = 1.0e-8;
    const double maxM2 = std::max(p2, std::max(m02, m12));
    const bool p2zero = (p2 <= maxM2*Tolerance);
    const bool m02zero = (m02 <= maxM2*Tolerance);
    const bool m12zero = (m12 <= maxM2*Tolerance);
    const bool m02_eq_m12 = (fabs(m02 - m12) <= (m02 + m12)*Tolerance);

    gslpp::complex B00p(0.0, 0.0, false);
    double DeltaM2 = m02 - m12;
    if ( p2zero ) {
        if ( m02_eq_m12 )
            B00p = - 1.0/18.0 - 1.0/12.0*B0(mu2,0.0,m02,m12)
                   + (m02 + m12)/6.0*B0p(mu2,0.0,m02,m12);
        else if ( m02zero || m12zero )
            B00p = - 1.0/18.0 - 1.0/12.0*B0(mu2,0.0,m02,m12)
                   + (m02 + m12)/6.0*B0p(mu2,0.0,m02,m12) - 1.0/72.0;
        else
            B00p = - 1.0/18.0 - 1.0/12.0*B0(mu2,0.0,m02,m12)
                   + (m02 + m12)/6.0*B0p(mu2,0.0,m02,m12)
                   - 1.0/24.0
                   *( (m02*m02 + 10.0*m02*m12 + m12*m12)/3.0/DeltaM2/DeltaM2
                       + 2.0*m02*m12*(m02 + m12)/DeltaM2/DeltaM2/DeltaM2*log(m12/m02) );
    } else {
        double Lambdabar2 = (p2-m02-m12)*(p2-m02-m12) - 4.0*m02*m12;
        if ( m02_eq_m12 )
            B00p = - 1.0/18.0 - B0(mu2,p2,m02,m12)/12.0
                   - Lambdabar2/12.0/p2*B0p(mu2,p2,m02,m12);
        else
            B00p = - 1.0/18.0
                   + DeltaM2/12.0/p2/p2*(- ExtraMinusSign*A0(mu2,m02)
                                         + ExtraMinusSign*A0(mu2,m12)
                                         + DeltaM2*B0(mu2,p2,m02,m12))
                   - B0(mu2,p2,m02,m12)/12.0
                   - Lambdabar2/12.0/p2*B0p(mu2,p2,m02,m12);
    }
    return B00p;
#endif
}

gslpp::complex PVfunctions::Bfp(const double mu2, const double p2,
                         const double m02, const double m12) const
{   
    if ( mu2<=0.0 || p2<0.0 || m02<0.0 || m12<0.0 )
        throw std::runtime_error("PVfunctions::Bfp(): Invalid argument!");

    gslpp::complex Bfp(0.0, 0.0, false);
    Bfp = 2.0*(B11p(mu2,p2,m02,m12) + B1p(mu2,p2,m02,m12));
    return Bfp;
}

gslpp::complex PVfunctions::C0(const double p2, 
                        const double m02, const double m12, const double m22) const
{
#ifdef USE_LOOPTOOLS
    return ( ExtraMinusSign * myLT.PV_C0(p2, m02, m12, m22) );
#else 
    if ( p2<0.0 || m02<0.0 || m12<0.0 || m22<0.0 )
        throw std::runtime_error("PVfunctions::C0(): Invalid argument!");

    const double Tolerance = 1.0e-8;
    const double maxM2 = std::max(p2, std::max(m02, std::max(m12, m22)));
    const bool p2zero = (p2 <= maxM2*Tolerance);
    const bool m02zero = (m02 <= maxM2*Tolerance);
    const bool m12zero = (m12 <= maxM2*Tolerance);
    const bool m22zero = (m22 <= maxM2*Tolerance);
    bool diff01 = (fabs(m02 - m12) > (m02 + m12)*Tolerance);
    bool diff12 = (fabs(m12 - m22) > (m12 + m22)*Tolerance);
    bool diff20 = (fabs(m22 - m02) > (m22 + m02)*Tolerance);

    gslpp::complex C0(0.0, 0.0, false);   
    if ( p2zero ) {
        if ( !m02zero && !m12zero && !m22zero ) {
            if ( diff01 && diff12 && diff20 )
                return ( - ( m12/(m02 - m12)*log(m12/m02)
                            - m22/(m02 - m22)*log(m22/m02) )/(m12 - m22) );
            else if ( !diff01 && diff12 && diff20 )
                return ( - (- m02 + m22 - m22*log(m22/m02))/(m02 - m22)/(m02 - m22) );
            else if ( diff01 && !diff12 && diff20 )
                return ( - (  m02 - m12 + m02*log(m12/m02))/(m02 - m12)/(m02 - m12) );
            else if ( diff01 && diff12 && !diff20 )
                return ( - (- m02 + m12 - m12*log(m12/m02))/(m02 - m12)/(m02 - m12) );
            else
                return ( 1.0/2.0/m02 );
        }
    } else {
        if ( !diff20 && diff01 ) {
            double epsilon = 1.0e-12;
            gsl_complex tmp = gsl_complex_rect(1.0 - 4.0*m02/p2, epsilon);
            tmp = gsl_complex_sqrt(tmp);
            gslpp::complex tmp_complex(GSL_REAL(tmp), GSL_IMAG(tmp), false);
            gslpp::complex x0 = 1.0 - (m02 - m12)/p2;
            gslpp::complex x1 = (1.0 + tmp_complex)/2.0;
            gslpp::complex x2 = (1.0 - tmp_complex)/2.0;
            gslpp::complex x3 = m02/(m02 - m12);

            if ( x0==x1 || x0==x2 || x0==x3)
                throw std::runtime_error("PVfunctions::C0(): Undefined-2!");

            gslpp::complex arg[6];
            arg[0] = (x0 - 1.0)/(x0 - x1);
            arg[1] = x0/(x0 - x1);
            arg[2] = (x0 - 1.0)/(x0 - x2);
            arg[3] = x0/(x0 - x2);
            arg[4] = (x0 - 1.0)/(x0 - x3);
            arg[5] = x0/(x0 - x3);
            
            gslpp::complex Li2[6];
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
        } else if ( !m02zero && !m22zero && diff20 && m12zero ) {
            double epsilon = 1.0e-12;
            double tmp_real = pow((m02+m22-p2),2.0) - 4.0*m02*m22;
            gsl_complex tmp = gsl_complex_rect(tmp_real, epsilon);
            tmp = gsl_complex_sqrt(tmp);
            gslpp::complex tmp_complex(GSL_REAL(tmp), GSL_IMAG(tmp), false);
            gslpp::complex x1 = (p2 - m02 + m22 + tmp_complex)/2.0/p2;
            gslpp::complex x2 = (p2 - m02 + m22 - tmp_complex)/2.0/p2;            

            if ( x1==0.0 || x1==1.0 || x2==0.0 || x2==1.0 )
                throw std::runtime_error("PVfunctions::C0(): Undefined-3!");
            
            gslpp::complex arg1 = (x1 - 1.0)/x1;
            gslpp::complex arg2 = x2/(x2 - 1.0);
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
        } else if ( m02zero && !m12zero && !m22zero ) {
            gslpp::complex arg[2];
            arg[0] = 1.0 - m22/m12;
            arg[1] = 1.0 - (-p2+m22)/m12;
            gslpp::complex Li2[2];
            for (int i=0; i<2; i++) {
                gsl_sf_result re, im;
                gsl_sf_complex_dilog_xy_e(arg[i].real(), arg[i].imag(), &re, &im);
                Li2[i].real() = re.val;
                Li2[i].imag() = im.val;
            }
            C0 = 1./(-p2)*(Li2[0]-Li2[1]);
        } else if ( !m02zero && !m12zero && !m22zero && diff01 && diff12 ) {
            double x0 = 1.0 - (m02-m12)/p2;
            double x1 = -(-p2+m02-m22-sqrt(fabs((m02+m22-p2)*(m02+m22-p2) - 4.*m02*m22)))/p2/2.0;
            double x2 = -(-p2+m02-m22+sqrt(fabs((m02+m22-p2)*(m02+m22-p2) - 4.*m02*m22)))/p2/2.0;
            double x3 = m22/(m22-m12);
            
            gslpp::complex arg[6];
            arg[0] = (x0-1.0)/(x0-x1);
            arg[1] = x0/(x0-x1);
            arg[2] = (x0-1.0)/(x0-x2);
            arg[3] = x0/(x0-x2);
            arg[4] = (x0-1.0)/(x0-x3);
            arg[5] = x0/(x0-x3);
            
            gslpp::complex Li2[6];
            for (int i=0; i<2; i++) {
                gsl_sf_result re, im;
                gsl_sf_complex_dilog_xy_e(arg[i].real(), arg[i].imag(), &re, &im);
                Li2[i].real() = re.val;
                Li2[i].imag() = im.val;
            }
            
            C0 = -1.0/p2*(Li2[0] - Li2[1] + Li2[2] - Li2[3] - Li2[4] + Li2[5]);
        } else
            throw std::runtime_error("PVfunctions::C0(): Undefined-4!");
    }
    return ( - ExtraMinusSign * C0 );
#endif    
}

double PVfunctions::C11(const double m12, const double m22, const double m32) const
{
    if ( m12<=0.0 || m22<=0.0 || m32<=0.0 )
        throw std::runtime_error("PVfunctions::C11(): Argument is not positive!");

    const double Tolerance = 2.5e-3;
    double C11;

    if ( 2.0*fabs(m12-m22) > (m12+m22)*Tolerance && 2.0*fabs(m12-m32) > (m12+m32)*Tolerance && 2.0*fabs(m22-m32) > (m22+m32)*Tolerance ) {
        C11=(-m12*(m12-m22)*(m12-m32)*(m22-m32)
             +m12*m12*m22*(2.0*m12-m22)*log(m12/m22)
             +m12*m12*m32*(-2.0*m12+m32)*log(m12/m32)
             +m22*m32*(-2.0*m12+m22)*(-2.0*m12+m32)*log(m22/m32))
            /(2.0*(m12-m22)*(m12-m22)*(m12-m32)*(m12-m32)*(m22-m32));
    }
    else if ( 2.0*fabs(m12-m22) > (m12+m22)*Tolerance && 2.0*fabs(m12-m32) > (m12+m32)*Tolerance && 2.0*fabs(m22-m32) <= (m22+m32)*Tolerance ) {
        C11=-((-12.0*m12*m12+8.0*m12*(m22+m32)-(m22+m32)*(m22+m32)+8.0*m12*m12*log((2.0*m12)/(m22+m32)))
              /pow(-2.0*m12+m22+m32,3));
    }
    else if ( 2.0*fabs(m12-m22) > (m12+m22)*Tolerance && 2.0*fabs(m12-m32) <= (m12+m32)*Tolerance && 2.0*fabs(m22-m32) > (m22+m32)*Tolerance ) {
        C11=(3.0*m12*m12-8.0*m12*m22+4.0*m22*m22+6.0*m12*m32-8.0*m22*m32+3.0*m32*m32+8.0*m22*(m12-m22+m32)*log((2.0*m22)/(m12+m32)))
            /(2.0*pow(m12-2.0*m22+m32,3));
    }
    else if ( 2.0*fabs(m12-m22) <= (m12+m22)*Tolerance && 2.0*fabs(m12-m32) <= (m12+m32)*Tolerance && 2.0*fabs(m22-m32) <= (m22+m32)*Tolerance ) {
        C11=1/(m12+m22+m32);
    }
    else {
        C11=(3.0*m12*m12+6.0*m12*m22+3.0*m22*m22-8.0*m12*m32-8.0*m22*m32+4.0*m32*m32-8.0*m32*(m12+m22-m32)*log((m12+m22)/(2.0*m32)))
            /(2.0*pow(m12+m22-2.0*m32,3));
    }

    return C11;
}

double PVfunctions::C12(const double m12, const double m22, const double m32) const
{
    if ( m12<=0.0 || m22<=0.0 || m32<=0.0 )
        throw std::runtime_error("PVfunctions::C12(): Argument is not positive!");

    const double Tolerance = 2.5e-3;
    double C12;

    if ( 2.0*fabs(m12-m22) > (m12+m22)*Tolerance && 2.0*fabs(m12-m32) > (m12+m32)*Tolerance && 2.0*fabs(m22-m32) > (m22+m32)*Tolerance ) {
        C12=((m12-m22)*(m12-m32)*(m22-m32)*m32
             +m12*m12*(m22*m22*log(m12/m22) +m32*(-2.0*m22+m32)*log(m12/m32))
             +m22*m22*(2.0*m12-m32)*m32*log(m22/m32))
            /(2.0*(m12-m22)*(m12-m32)*(m12-m32)*(m22-m32)*(m22-m32));
    }
    else if ( 2.0*fabs(m12-m22) > (m12+m22)*Tolerance && 2.0*fabs(m12-m32) > (m12+m32)*Tolerance && 2.0*fabs(m22-m32) <= (m22+m32)*Tolerance ) {
        C12=-(-12.0*m12*m12+8.0*m12*(m22+m32)-(m22+m32)*(m22+m32)+8.0*m12*m12*log((2.0*m12)/(m22+m32)))
             /(2.0*pow(-2.0*m12+m22+m32,3));
    }
    else if ( 2.0*fabs(m12-m22) > (m12+m22)*Tolerance && 2.0*fabs(m12-m32) <= (m12+m32)*Tolerance && 2.0*fabs(m22-m32) > (m22+m32)*Tolerance ) {
        C12=(m12*m12-8.0*m12*m22+12.0*m22*m22+2.0*m12*m32-8.0*m22*m32+m32*m32-8.0*m22*m22*log((2.0*m22)/(m12+m32)))
            /(2.0*pow(m12-2.0*m22+m32,3));
    }
    else if ( 2.0*fabs(m12-m22) <= (m12+m22)*Tolerance && 2.0*fabs(m12-m32) <= (m12+m32)*Tolerance && 2.0*fabs(m22-m32) <= (m22+m32)*Tolerance ) {
        C12=1.0/(2.0*(m12+m22+m32));
    }
    else {
        C12=(m12*m12+2.0*m12*m22+m22*m22-4.0*m32*m32-4.0*(m12+m22)*m32*log((m12+m22)/(2.0*m32)))
            /pow(m12+m22-2.0*m32,3);
    }

    return C12;
}

gslpp::complex PVfunctions::D0(const double s, const double t, const double m02,
                        const double m12, const double m22, const double m32) const
{
    if ( m02<0.0 || m12<0.0 || m22<0.0 || m32<0.0 )
        throw std::runtime_error("PVfunctions::D0(): Invalid argument!");

    const double Tolerance = 1.0e-8;
    const double maxM2 = std::max(s, std::max(t, std::max(m02, std::max(m12, std::max(m22, m32)))));
    const bool szero = (s <= maxM2*Tolerance);
    const bool tzero = (t <= maxM2*Tolerance);
    const bool m02zero = (m02 <= maxM2*Tolerance);
    const bool m12zero = (m12 <= maxM2*Tolerance);
    const bool m22zero = (m22 <= maxM2*Tolerance);
    const bool m32zero = (m32 <= maxM2*Tolerance);
    bool diff01 = (fabs(m02 - m12) > (m02 + m12)*Tolerance);
    bool diff02 = (fabs(m02 - m22) > (m02 + m22)*Tolerance);
    bool diff03 = (fabs(m02 - m32) > (m02 + m32)*Tolerance);
    bool diff23 = (fabs(m22 - m32) > (m22 + m32)*Tolerance);

    if ( szero && tzero ) {
        if ( diff01 )
            return ( ( ExtraMinusSign * C0(0.0, m02, m22, m32)
                       - ExtraMinusSign * C0(0.0, m12, m22, m32) ) / (m02 - m12) );
        else {
            if ( !m02zero && !m12zero && !m22zero && !m32zero ) {
                if ( diff02 && diff03 && diff23 )
                    return ( - 1.0/(m02 - m22)/(m02 - m32)
                             + m22/(m02 - m22)/(m02 - m22)/(m32 - m22)*log(m22/m02)
                             + m32/(m02 - m32)/(m02 - m32)/(m22 - m32)*log(m32/m02) );
                else if ( !diff02 && diff03 && diff23 )
                    return ( (m02*m02 - m32*m32 + 2.0*m02*m32*log(m32/m02))
                             /2.0/m02/(m02 - m32)/(m02 - m32)/(m02 - m32) );
                else if ( diff02 && !diff03 && diff23 )
                    return ( (m02*m02 - m22*m22 + 2.0*m02*m22*log(m22/m02))
                             /2.0/m02/(m02 - m22)/(m02 - m22)/(m02 - m22) );
                else if ( diff02 && diff03 && !diff23 )
                    return ( ( - 2.0*m02 + 2.0*m22 - (m02 + m22)*log(m22/m02))
                             /(m02 - m22)/(m02 - m22)/(m02 - m22) );
                else
                    return ( 1.0/6.0/m02/m02 );
            } else
                throw std::runtime_error("PVfunctions::D0(): Undefined!");
        }
    } 

#ifdef USE_LOOPTOOLS
    return myLT.PV_D0(s, t, m02, m12, m22, m32);
#else
    gslpp::complex D0(0.0, 0.0, false);

    if ( s>0.0 && t<0.0 && !m02zero && m12zero && !diff02 && !m32zero
            && m02!=m32 && t-m32+m02!=0.0 && t-m32!=0.0 ) {
        //D0(s,t; m02, 0.0, m02, m32)

        /*
        double x1, x2;
        if ( s >= 4.0*m02 ) {
            x1 = (1.0 - sqrt(1.0 - 4.0*m02/s))/2.0;
            x2 = (1.0 + sqrt(1.0 - 4.0*m02/s))/2.0;        
        } else {
            throw std::runtime_error("PVfunctions::D0(): Undefined!");
        }
        double x3 = m32/(m32 - m02);
        double x4 = (t - m32)/(t - m32 + m02);
        double d4 = 1.0 - 4.0*m02*t*(t - m32 + m02)/(s*(t - m32)*(t - m32));
        double x1tilde, x2tilde;
        if ( d4 >= 0.0 ) {
            x1tilde = x4/2.0*(1.0 - sqrt(d4));
            x2tilde = x4/2.0*(1.0 + sqrt(d4));
        } else {
            throw std::runtime_error("PVfunctions::D0(): Undefined!");
        }
         */

        /* Write codes! */

        throw std::runtime_error("PVfunctions::D0(): Undefined!");
    } else if ( s>0.0 && t<0.0 && !m02zero && m12zero && !diff02 && m32zero ) {
        //D0(s,t; m02, 0.0, m02, 0.0)
        
        throw std::runtime_error("PVfunctions::D0(): Undefined!");
    } else 
        throw std::runtime_error("PVfunctions::D0(): Undefined!");
        
    return D0; 
#endif
}

gslpp::complex PVfunctions::D00(const double s, const double t, const double m02,
                         const double m12, const double m22, const double m32) const
{
    if ( m02<0.0 || m12<0.0 || m22<0.0 || m32<0.0 )
        throw std::runtime_error("PVfunctions::D00(): Invalid argument!");

    const double Tolerance = 1.0e-8;
    const double maxM2 = std::max(s, std::max(t, std::max(m02, std::max(m12, std::max(m22, m32)))));
    const bool szero = (s <= maxM2*Tolerance);
    const bool tzero = (t <= maxM2*Tolerance);
    const bool m02zero = (m02 <= maxM2*Tolerance);
    const bool m12zero = (m12 <= maxM2*Tolerance);
    const bool m22zero = (m22 <= maxM2*Tolerance);
    const bool m32zero = (m32 <= maxM2*Tolerance);
    bool diff01 = (fabs(m02 - m12) > (m02 + m12)*Tolerance);
    bool diff02 = (fabs(m02 - m22) > (m02 + m22)*Tolerance);
    bool diff03 = (fabs(m02 - m32) > (m02 + m32)*Tolerance);
    bool diff23 = (fabs(m22 - m32) > (m22 + m32)*Tolerance);

    if ( szero && tzero ) {
        if ( diff01 )
            return ( 0.25/(m02 - m12)*(m02*ExtraMinusSign*C0(0.0, m02, m22, m32)
                     - m12*ExtraMinusSign*C0(0.0, m12, m22, m32)) );
        else {
            if ( !m02zero && !m12zero && !m22zero && !m32zero ) {
                if ( diff02 && diff03 && diff23 )
                    return ( m02/4.0/(m02 - m22)/(m32 - m02)
                             + m22*m22/4.0/(m02 - m22)/(m02 - m22)/(m32 - m22)*log(m22/m02)
                             + m32*m32/4.0/(m02 - m32)/(m02 - m32)/(m22 - m32)*log(m32/m02) );
                else if ( !diff02 && diff03 && diff23 )
                    return ( ( - m02*m02 + 4.0*m02*m32 - 3.0*m32*m32
                               + 2.0*m32*m32*log(m32/m02) )
                             /8.0/(m02 - m32)/(m02 - m32)/(m02 - m32) );
                else if ( diff02 && !diff03 && diff23 )
                    return ( ( - m02*m02 + 4.0*m02*m22 - 3.0*m22*m22
                               + 2.0*m22*m22*log(m22/m02) )
                             /8.0/(m02 - m22)/(m02 - m22)/(m02 - m22) );
                else if ( diff02 && diff03 && !diff23 )
                    return ( ( - m02*m02 + m22*m22 - 2.0*m02*m22*log(m22/m02) )
                             /4.0/(m02 - m22)/(m02 - m22)/(m02 - m22) );
                else
                    return ( - 1.0/12.0/m02 );
            } else
                throw std::runtime_error("PVfunctions::D00(): Undefined!");
        }
    }

#ifdef USE_LOOPTOOLS
    return myLT.PV_D00(s, t, m02, m12, m22, m32);
#else
    gslpp::complex D00(0.0, 0.0, false);

    throw std::runtime_error("PVfunctions::D00(): Undefined!");

    return D00;
#endif
}

