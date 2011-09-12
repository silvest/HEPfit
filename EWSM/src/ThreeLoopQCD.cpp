/* 
 * File:   ThreeLoopQCD.cpp
 * Author: mishima
 */

#include "ThreeLoopQCD.h"


ThreeLoopQCD::ThreeLoopQCD(const EWSMcommon& EWSMC_i) : EWSMC(EWSMC_i) {
}

//ThreeLoopQCD::ThreeLoopQCD(const ThreeLoopQCD& orig) {
//}

ThreeLoopQCD::~ThreeLoopQCD() {
}


////////////////////////////////////////////////////////////////////////

double ThreeLoopQCD::DeltaAlpha_l() const {
    return (0.0);
}    

double ThreeLoopQCD::DeltaAlpha_t() const {   
    double xt = pow(EWSMC.GetSM().getMz()
                /EWSMC.GetSM().getQuarks(EWSMC.GetSM().TOP).getMass(), 2.0);
    double log_t = 2.0*EWSMC.GetLogMZtoMTOP();

    double tmp = ( (28.220 + 9.702*log_t) 
                   + xt*(6.924 + 1.594*log_t) )
                 *pow(EWSMC.GetSM().getAlsMz()/M_PI, 2.0);
    tmp *= -4.0/45.0*EWSMC.GetSM().getAle()/M_PI*xt;
    return tmp;
}

double ThreeLoopQCD::DeltaRho() const {
    return ( 3.0*EWSMC.GetXt_alpha()*pow(EWSMC.GetAlsMt()/M_PI,2.0)*deltaQCD_3());     
}

double ThreeLoopQCD::DeltaR_rem() const {
    return (0.0);     
}

complex ThreeLoopQCD::deltaRho_rem_l(const StandardModel::lepton l) const {
    complex zero(0.0,0.0,false);
    return zero;      
}

complex ThreeLoopQCD::deltaRho_rem_q(const StandardModel::quark q) const {
    complex zero(0.0,0.0,false);
    return zero;        
}

complex ThreeLoopQCD::deltaKappa_rem_l(const StandardModel::lepton l) const {
    return ( - 3.0*EWSMC.GetXt_alpha()*EWSMC.GetCW2()/EWSMC.GetSW2()
               *pow(EWSMC.GetAlsMt()/M_PI,2.0)
               *(deltaQCD_3()+deltaQCD_kappa3()) ); 
}

complex ThreeLoopQCD::deltaKappa_rem_q(const StandardModel::quark q) const {
    return ( - 3.0*EWSMC.GetXt_alpha()*EWSMC.GetCW2()/EWSMC.GetSW2()
               *pow(EWSMC.GetAlsMt()/M_PI,2.0)
               *(deltaQCD_3()+deltaQCD_kappa3()) ); 
}


////////////////////////////////////////////////////////////////////////

double ThreeLoopQCD::deltaQCD_3() const {
    double dQCD3;
    double lZ = 2.0*EWSMC.GetLogMZtoMTOP();
    double sW2 = EWSMC.GetSW2();
    double log2 = EWSMC.GetLog2();
    double zeta2 = EWSMC.GetZeta2();
    double zeta3 = EWSMC.GetZeta3();    
    double zeta4 = EWSMC.GetZeta4();
    double S2 = EWSMC.GetS2(), D3 = EWSMC.GetD3(), B4 = EWSMC.GetB4();
    double MZtoMT = EWSMC.GetSM().getMz()
                      /EWSMC.GetSM().getQuarks(EWSMC.GetSM().TOP).getMass();
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

complex ThreeLoopQCD::deltaQCD_kappa3() const {
    complex dQCDk3;
    double lZ = 2.0*EWSMC.GetLogMZtoMTOP();
    double sW2 = EWSMC.GetSW2();
    double MZtoMT = EWSMC.GetSM().getMz()
                    /EWSMC.GetSM().getQuarks(EWSMC.GetSM().TOP).getMass();
    dQCDk3.real() = - deltaQCD_3()
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





