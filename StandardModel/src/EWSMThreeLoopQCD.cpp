/* 
 * File:   EWSMThreeLoopQCD.cpp
 * Author: mishima
 */

#include "EWSMThreeLoopQCD.h"


EWSMThreeLoopQCD::EWSMThreeLoopQCD(const EWSMcache& cache_i) : cache(cache_i) {
}


////////////////////////////////////////////////////////////////////////

double EWSMThreeLoopQCD::DeltaAlpha_l(const double s) const {
    return (0.0);
}    


double EWSMThreeLoopQCD::DeltaAlpha_t(const double s) const {   
    double xt = s/cache.Mt()/cache.Mt();
    double log_t, als;
    if (s==cache.Mz()*cache.Mz()) {
        log_t = 2.0*cache.logMZtoMTOP();
        als = cache.alsMz();
    } else {
        log_t = log(s/cache.mq(StandardModel::TOP)/cache.mq(StandardModel::TOP));
        als = cache.getSM().Als(sqrt(s),FULLNNLO);
    }
    double tmp = ( (28.220 + 9.702*log_t) 
                   + xt*(6.924 + 1.594*log_t) )
                 *pow(als/M_PI, 2.0);
    tmp *= -4.0/45.0*cache.ale()/M_PI*xt;
    return tmp;
}


double EWSMThreeLoopQCD::DeltaRho(const double Mw_i) const {
    double Mw = cache.Mw(Mw_i);
    return ( 3.0*cache.Xt_alpha(Mw)*pow(cache.alsMt()/M_PI,2.0)*deltaQCD_3(Mw));     
}


double EWSMThreeLoopQCD::DeltaR_rem(const double Mw_i) const {
    double Mw = cache.Mw(Mw_i);
    double sW2 = cache.sW2(Mw);
    double cW2 = cache.cW2(Mw);
    
    /* Logarithm */
    double log_cW2 = cache.log_cW2(Mw);     
    
    // O(alpha_s) correction to Delta r^{ud} of O(alpha alpha_s). 
    double DeltaR;
    DeltaR = - log_cW2;
    DeltaR *= (cW2 - sW2)/4.0/sW2/sW2;
    DeltaR *= cache.ale()*cache.alsMz()/M_PI/M_PI;
    DeltaR *= 1.4097*cache.alsMz()/M_PI;
    return DeltaR;     
}


complex EWSMThreeLoopQCD::deltaRho_rem_l(const StandardModel::lepton l, 
                                         const double Mw_i) const {
    return ( complex(0.0,0.0,false) );
}


complex EWSMThreeLoopQCD::deltaRho_rem_q(const StandardModel::quark q, 
                                         const double Mw_i) const {
    if(q==StandardModel::TOP) return ( complex(0.0,0.0,false) );
    return ( complex(0.0,0.0,false) );
}


complex EWSMThreeLoopQCD::deltaKappa_rem_l(const StandardModel::lepton l, 
                                           const double Mw_i) const {
    double Mw = cache.Mw(Mw_i);
    return ( - 3.0*cache.Xt_alpha(Mw)*cache.cW2(Mw)/cache.sW2(Mw)
               *pow(cache.alsMt()/M_PI,2.0)
               *(deltaQCD_3(Mw)+deltaQCD_kappa3(Mw).real()) ); 
}


complex EWSMThreeLoopQCD::deltaKappa_rem_q(const StandardModel::quark q, 
                                           const double Mw_i) const {
    if(q==StandardModel::TOP) return ( complex(0.0,0.0,false) );
    double Mw = cache.Mw(Mw_i);
    return ( - 3.0*cache.Xt_alpha(Mw)*cache.cW2(Mw)/cache.sW2(Mw)
               *pow(cache.alsMt()/M_PI,2.0)
               *(deltaQCD_3(Mw)+deltaQCD_kappa3(Mw).real()) ); 
}


////////////////////////////////////////////////////////////////////////

double EWSMThreeLoopQCD::deltaQCD_3(const double Mw_i) const {
    double dQCD3;
    double lZ = 2.0*cache.logMZtoMTOP();
    double Mw = cache.Mw(Mw_i);
    double sW2 = cache.sW2(Mw);
    double log2 = cache.GetLog2();
    double zeta2 = cache.GetZeta2();
    double zeta3 = cache.GetZeta3();    
    double zeta4 = cache.GetZeta4();
    double S2 = cache.GetS2(), D3 = cache.GetD3(), B4 = cache.GetB4();
    double MZtoMT = cache.Mz()/cache.Mt();
    double nf = 6.0;
    dQCD3 = 157.0/648.0 - 3313.0/162.0*zeta2 - 308.0/27.0*zeta3 
            + 143.0/18.0*zeta4 - 4.0/3.0*zeta2*log2 
            + 441.0/8.0*S2 - B4/9.0 - D3/18.0
            - ( 1.0/18.0 - 13.0/9.0*zeta2 + 4.0/9.0*zeta3 )*nf 
            + pow(MZtoMT, 2.0)
              *( - 17.224 + 0.08829*lZ + 0.4722*lZ*lZ
                 + ( 22.6367 + 1.2527*lZ - 0.8519*lZ*lZ )*sW2 )
            + pow(MZtoMT, 4.0)
              *( - 7.7781 - 0.07226*lZ + 0.004938*lZ*lZ
                 + ( 21.497 + 0.05794*lZ - 0.006584*lZ*lZ )*sW2
                 - 21.0799*sW2*sW2 );
    return dQCD3;
}


complex EWSMThreeLoopQCD::deltaQCD_kappa3(const double Mw_i) const {
    complex dQCDk3;
    double lZ = 2.0*cache.logMZtoMTOP();
    double Mw = cache.Mw(Mw_i);
    double sW2 = cache.sW2(Mw);
    double MZtoMT = cache.Mz()/cache.Mt();
    dQCDk3.real() = - deltaQCD_3(Mw)
                    + pow(MZtoMT, 2.0)
                      *( ( 22.6367 + 1.2527*lZ - 0.8519*lZ*lZ )*sW2
                         + ( - 11.3184 - 0.6263*lZ + 0.4259*lZ*lZ )*sW2 )
                    + pow(MZtoMT, 4.0)
                      *( ( 21.497 + 0.05794*lZ - 0.006584*lZ*lZ )*sW2
                         + ( - 16.0186 - 0.02897*lZ + 0.003292*lZ*lZ )*sW2
                         - 21.0799*sW2*sW2 + 10.54*sW2*sW2 );
    dQCDk3.imag() = pow(MZtoMT, 2.0)
                      *( ( - 1.968 + 2.676*lZ )*sW2 
                         + ( 2.6235 - 3.5682*lZ )*sW2*sW2 )
                    + pow(MZtoMT, 4.0)
                      *( ( - 0.09102 + 0.02069*lZ )*sW2 
                         + ( 0.1214 - 0.02758*lZ )*sW2*sW2 );
    return dQCDk3;
}







