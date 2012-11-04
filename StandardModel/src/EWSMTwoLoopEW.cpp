 /* 
 * Copyright (C) 2012 SUSYfit Collaboration
 * All rights reserved.
 *
 * For the licensing terms see doc/COPYING.
 */

#include <TMath.h>
#include "EWSMTwoLoopEW.h"

/* include O(alpha^2 M_t^2/M_Z^2) in addition to O(alpha^2 M_t^4/M_Z^4) */
//#define EW_SUBLEADING_ALPHA2


EWSMTwoLoopEW::EWSMTwoLoopEW(const EWSMcache& cache_i) : cache(cache_i), myOneLoopEW(cache_i) {
}


////////////////////////////////////////////////////////////////////////

double EWSMTwoLoopEW::DeltaAlpha_l(const double s) const {
    double xl[3] = { s/cache.ml(StandardModel::ELECTRON)/cache.ml(StandardModel::ELECTRON), 
                     s/cache.ml(StandardModel::MU)/cache.ml(StandardModel::MU), 
                     s/cache.ml(StandardModel::TAU)/cache.ml(StandardModel::TAU) };
    double log_l[3];
    if (s==cache.Mz()*cache.Mz()) {
        log_l[0] = 2.0*cache.logMZtoME();
        log_l[1] = 2.0*cache.logMZtoMMU();
        log_l[2] = 2.0*cache.logMZtoMTAU();
    } else {
        log_l[0] = log(xl[0]);
        log_l[1] = log(xl[1]);
        log_l[2] = log(xl[2]);
    }    
    
    double twoLoop[3];
    for (int i = 0; i < 3; i++) {
        twoLoop[i] = - 5.0/24.0 + cache.GetZeta3() + log_l[i]/4.0 
                     + 3.0/xl[i]*log_l[i];
    }
            
    return ( pow(cache.ale()/M_PI, 2.0)
             *(twoLoop[0] + twoLoop[1] + twoLoop[2]) );
}    


double EWSMTwoLoopEW::DeltaAlpha_t(const double s) const {   
    return (0.0);
}


double EWSMTwoLoopEW::DeltaRho(const double Mw_i) const {
    double Mz = cache.Mz();
    double Mw = cache.Mw(Mw_i);
    double sW2 = cache.sW2(Mw);
    double cW2 = cache.cW2(Mw);
   
    double DeltaRho = 0.0;
    
    #ifndef EW_SUBLEADING_ALPHA2
    /* O(\alpha^2 Mt^4/Mz^4) */
    DeltaRho += 3.0*rho_2();
    DeltaRho *= pow(cache.Xt_alpha(Mw), 2.0);
    #else
    /* O(\alpha^2 Mt^4/Mz^4 + \alpha^2 Mt^2/Mz^2) */
    double zt = Mz*Mz/cache.Mt()/cache.Mt();
    DeltaRho += 3.0*pow(cache.Xt_alpha(Mw), 2.0)
                *( DeltaRho2(Mw) + 4.0*zt*cW2*DeltaRho2Add(Mw) );
    #endif
    
    /* add O(alpha^2) contribution from the Z-gamma mixing */
    DeltaRho += - pow(cache.ale()/4.0/M_PI, 2.0)*cW2/sW2
                  *pow(myOneLoopEW.PiZgamma_fer(Mz,Mz*Mz,Mw).real(), 2.0);

    return DeltaRho;
}


double EWSMTwoLoopEW::DeltaR_rem(const double Mw_i) const {
    double DeltaRrem = 0.0;
    
    #ifdef EW_SUBLEADING_ALPHA2
    /* O(\alpha^2 Mt^4/Mz^4 + \alpha^2 Mt^2/Mz^2) */
    double Mz = cache.Mz();
    double Mw = cache.Mw(Mw_i);
    double sW2 = cache.sW2(Mw);
    DeltaRrem += 3.0*pow(cache.ale()*cache.Mt()/4.0/M_PI/sW2/Mw, 2.0)
                 *( DeltaRw2(Mw) + sW2*deltaEoverE2(Mw) + f2Add(Mw)/4.0 );
    #endif    
    
    return DeltaRrem;     
}


complex EWSMTwoLoopEW::deltaRho_rem_l(const StandardModel::lepton l, 
                                      const double Mw_i) const {
    complex dRho = complex(0.0,0.0,false);
    
    #ifdef EW_SUBLEADING_ALPHA2
    /* O(\alpha^2 Mt^4/Mz^4 + \alpha^2 Mt^2/Mz^2) */
    double Mz = cache.Mz();
    double Mw = cache.Mw(Mw_i);
    double cW2 = cache.cW2(Mw);
    double zt = Mz*Mz/cache.Mt()/cache.Mt();
    dRho += 3.0*pow(cache.Xt_alpha(Mw), 2.0)
            *( 16.0*zt*cW2*DeltaEta2(Mw) + 4.0*zt*cW2*DeltaEta2Add_l(l, Mw) );
    #endif 
    
    return dRho;
}


complex EWSMTwoLoopEW::deltaRho_rem_q(const StandardModel::quark q, 
                                      const double Mw_i) const {
    if(q==StandardModel::TOP) return ( complex(0.0,0.0,false) );

    complex dRho = complex(0.0,0.0,false);
    
    #ifdef EW_SUBLEADING_ALPHA2
    /* O(\alpha^2 Mt^4/Mz^4 + \alpha^2 Mt^2/Mz^2) */
    double Mz = cache.Mz();
    double Mw = cache.Mw(Mw_i);
    double cW2 = cache.cW2(Mw);
    double zt = Mz*Mz/cache.Mt()/cache.Mt();
    dRho += 3.0*pow(cache.Xt_alpha(Mw), 2.0)
            *( 16.0*zt*cW2*DeltaEta2(Mw) + 4.0*zt*cW2*DeltaEta2Add_q(q, Mw) );
    #endif 
    
    return dRho;
}


complex EWSMTwoLoopEW::deltaKappa_rem_l(const StandardModel::lepton l, 
                                        const double Mw_i) const {
    complex dKappa = complex(0.0,0.0,false);
    
    #ifdef EW_SUBLEADING_ALPHA2
    /* O(\alpha^2 Mt^4/Mz^4 + \alpha^2 Mt^2/Mz^2) */
    double Mz = cache.Mz();
    double Mw = cache.Mw(Mw_i);
    double cW2 = cache.cW2(Mw);
    double zt = Mz*Mz/cache.Mt()/cache.Mt();
    dKappa += 3.0*pow(cache.Xt_alpha(Mw), 2.0)
              *( 16.0*zt*cW2*DeltaKappa2(Mw) + 4.0*zt*cW2*DeltaKappa2Add_l(l, Mw) );
    #endif 
    
    return dKappa;
}


complex EWSMTwoLoopEW::deltaKappa_rem_q(const StandardModel::quark q, 
                                        const double Mw_i) const {
    if(q==StandardModel::TOP) return ( complex(0.0,0.0,false) );

    complex dKappa = complex(0.0,0.0,false);
    
    #ifdef EW_SUBLEADING_ALPHA2
    /* O(\alpha^2 Mt^4/Mz^4 + \alpha^2 Mt^2/Mz^2) */
    double Mz = cache.Mz();
    double Mw = cache.Mw(Mw_i);
    double cW2 = cache.cW2(Mw);
    double zt = Mz*Mz/cache.Mt()/cache.Mt();
    dKappa += 3.0*pow(cache.Xt_alpha(Mw), 2.0)
              *( 16.0*zt*cW2*DeltaKappa2(Mw) + 4.0*zt*cW2*DeltaKappa2Add_q(q, Mw) );
    #endif 
    
    return dKappa;
}


////////////////////////////////////////////////////////////////////////   `

double EWSMTwoLoopEW::rho_2() const {
    double a = cache.mh()*cache.mh()/cache.Mt()/cache.Mt();
    if (a<=0.0) throw std::runtime_error("a is out of range in EWSMTwoLoopEW::rho_2"); 
    double g_a = g(a);
    double f_a_0 = f0(a);// f(a,0)
    double f_a_1 = f1(a);// f(a,1)
    double log_a = - 2.0*cache.logMTOPtoMH();
    return ( 25.0 - 4.0*a + 0.5*(a*a - 12.0*a - 12.0)*log_a
             + (a - 2.0)/2.0/a*M_PI*M_PI + 0.5*(a - 4.0)*sqrt(a)*g_a
             - 3.0/a*(a - 1.0)*(a - 1.0)*(a - 2.0)*f_a_0
             + 3.0*(a*a - 6.0*a + 10.0)*f_a_1 );
}


double EWSMTwoLoopEW::tau_2() const {
    double a = cache.mh()*cache.mh()/cache.Mt()/cache.Mt();
    if (a<=0.0) throw std::runtime_error("a is out of range in EWSMTwoLoopEW::tau_2"); 
    double g_a = g(a);
    double f_a_0 = f0(a);// f(a,0)
    double f_a_1 = f1(a);// f(a,1)
    double log_a = - 2.0*cache.logMTOPtoMH();
    return ( 9.0 - 13.0/4.0*a - 2.0*a*a - a/4.0*(19.0 + 6.0*a)*log_a
             - a*a/4.0*(7.0 - 6.0*a)*log_a*log_a
             - (1.0/4.0 + 7.0/2.0*a*a - 3.0*a*a*a)*M_PI*M_PI/6.0
             + (a/2.0 - 2.0)*sqrt(a)*g_a 
             + (a - 1.0)*(a - 1.0)*(4.0*a - 7.0/4.0)*f_a_0
             - (a*a*a - 33.0/4.0*a*a + 18.0*a - 7.0)*f_a_1 );
}


double EWSMTwoLoopEW::g(const double a) const {
    if (a >= 0.0 && a <= 4.0 ) {
        double phi = 2.0*asin(sqrt(a/4.0));
        return ( sqrt(4.0 - a)*(M_PI - phi) );
    } else if (a > 4.0) {
        double y = 4.0/a;
        double xi = (sqrt(1.0-y) - 1.0)/(sqrt(1.0-y) + 1.0);
        return ( sqrt(a - 4.0)*log(-xi) );
    } else
        throw std::runtime_error("Out of range in EWSMTwoLoopEW::g()"); 
}


double EWSMTwoLoopEW::f0(const double a) const {
    if (a >= 0.0 )
        return ( cache.getPolyLog().Li2(1.0-a).real() ); // 1-a<1
    else
        throw std::runtime_error("Out of range in EWSMTwoLoopEW::f0()"); 
}


double EWSMTwoLoopEW::f1(const double a) const {
    if (a >= 0.0 && a <= 4.0 ) {
        double y = 4.0/a;
        double phi = 2.0*asin(sqrt(a/4.0));
        return ( -2.0/sqrt(y-1.0)*cache.getClausen().Cl2(phi) );
    } else if (a > 4.0) {
        double y = 4.0/a; // 0<y<1
        double xi = (sqrt(1.0-y) - 1.0)/(sqrt(1.0-y) + 1.0); // -1<xi<0
        return ( -1.0/sqrt(1.0-y)*(cache.getPolyLog().Li2(xi).real()
                                   - cache.getPolyLog().Li2(1.0/xi).real()) );
    } else
        throw std::runtime_error("Out of range in EWSMTwoLoopEW::f1()"); 
}


////////////////////////////////////////////////////////////////////////        
// O(alpha^2 Mt^4/M_Z^4 + alpha^2 Mt^2/M_Z^2) contributions

double EWSMTwoLoopEW::DeltaRho2(const double Mw_i) const {
    
    /* Write codes!! */

    return 0.0;
}


double EWSMTwoLoopEW::DeltaRho2Add(const double Mw_i) const {
    
    /* Write codes!! */

    return 0.0;
}


double EWSMTwoLoopEW::DeltaRw2(const double Mw_i) const {
    
    /* Write codes!! */

    return 0.0;
}


double EWSMTwoLoopEW::deltaEoverE2(const double Mw_i) const {
    
    /* Write codes!! */

    return 0.0;
}


double EWSMTwoLoopEW::f2Add(const double Mw_i) const {
    
    /* Write codes!! */

    return 0.0;
}


double EWSMTwoLoopEW::DeltaEta2(const double Mw_i) const {
    
    /* Write codes!! */

    return 0.0;
}


complex EWSMTwoLoopEW::DeltaEta2Add_l(const StandardModel::lepton l, 
                                      const double Mw_i) const {
    
    /* Write codes!! */

    return complex(0.0, 0.0, false);
}


complex EWSMTwoLoopEW::DeltaEta2Add_q(const StandardModel::quark q, 
                                      const double Mw_i) const {
    
    /* Write codes!! */

    return complex(0.0, 0.0, false);
}


double EWSMTwoLoopEW::DeltaKappa2(const double Mw_i) const {
    
    /* Write codes!! */

    return 0.0;
}


complex EWSMTwoLoopEW::DeltaKappa2Add_l(const StandardModel::lepton l, 
                                        const double Mw_i) const {
    
    /* Write codes!! */

    return complex(0.0, 0.0, false);
}


complex EWSMTwoLoopEW::DeltaKappa2Add_q(const StandardModel::quark q, 
                                        const double Mw_i) const {
    
    /* Write codes!! */

    return complex(0.0, 0.0, false);
}


complex EWSMTwoLoopEW::Vadd_l(const StandardModel::lepton l, 
                              const double Mw_i) const {
    
    /* Write codes!! */

    return complex(0.0, 0.0, false);
}


complex EWSMTwoLoopEW::Vadd_q(const StandardModel::quark q, 
                              const double Mw_i) const {
    
    /* Write codes!! */

    return complex(0.0, 0.0, false);
}


complex EWSMTwoLoopEW::DeltaEtaf1_l(const StandardModel::lepton l, 
                                    const double Mw_i) const {
    
    /* Write codes!! */

    return complex(0.0, 0.0, false);
}


complex EWSMTwoLoopEW::DeltaEtaf1_q(const StandardModel::quark q, 
                                    const double Mw_i) const {
    
    /* Write codes!! */

    return complex(0.0, 0.0, false);
}


complex EWSMTwoLoopEW::Vfi_l(const StandardModel::lepton l, 
                             const double Mw_i) const {
    
    /* Write codes!! */

    return complex(0.0, 0.0, false);
}


complex EWSMTwoLoopEW::Vfi_q(const StandardModel::quark q, 
                             const double Mw_i) const {
    
    /* Write codes!! */

    return complex(0.0, 0.0, false);
}


double EWSMTwoLoopEW::Lambda(const double x) const {
    
    /* Write codes!! */

    return 0.0;
}


double EWSMTwoLoopEW::phi(const double x) const {
    
    /* Write codes!! */

    return 0.0;
}


complex EWSMTwoLoopEW::FV(const double x) const {
    
    /* Write codes!! */

    return complex(0.0, 0.0, false);
}


complex EWSMTwoLoopEW::GV(const double x) const {
    
    /* Write codes!! */

    return complex(0.0, 0.0, false);
}

  


