/* 
 * File:   ThreeLoopEW2QCD.cpp
 * Author: mishima
 */

#include <TError.h>

#include "ThreeLoopEW2QCD.h"


ThreeLoopEW2QCD::ThreeLoopEW2QCD(const EWSMcommon& EWSMC_i) : EWSMC(EWSMC_i) {
}


////////////////////////////////////////////////////////////////////////

double ThreeLoopEW2QCD::DeltaAlpha_l() const {
    return (0.0);
}    

double ThreeLoopEW2QCD::DeltaAlpha_t() const {   
    return (0.0);
}

double ThreeLoopEW2QCD::DeltaRho() const {
    double mh = EWSMC.GetSM().getMHl();
    double Mt = EWSMC.GetSM().getQuarks(EWSMC.GetSM().TOP).getMass();
    double DeltaRho;
    if (mh==0.0) {
        DeltaRho = 4.0*( 185.0/3.0 + 729.0/4.0*EWSMC.GetS2() 
                         - 48.0*EWSMC.GetZeta2()*EWSMC.GetLog2()
                         - 151.0/6.0*EWSMC.GetZeta2() + 29.0*EWSMC.GetZeta3()
                         - 24.0*EWSMC.GetZeta4() + 12.0*EWSMC.GetB4() );
    } else if (mh > 0.0 && mh <= 2.5*Mt) {
        double delta = mh/Mt -1.0;
        DeltaRho = 157.295 + 112.00*delta - 24.73*delta*delta
                   + 7.39*pow(delta, 3.0) - 3.52*pow(delta, 4.0) 
                   + 2.06*pow(delta, 5.0);        
    } else if (mh > 2.5*Mt) {
        double Y = 4.0*pow(Mt/mh,2.0);
        double logY = 2.0*(EWSMC.GetLog2() + EWSMC.GetLogMTOPtoMH());
        double logY2 = logY*logY;
        double logY3 = logY2*logY;
        DeltaRho = 79.73 - 47.77*logY + 42.07*logY2 + 9.00*logY3
                   +  Y*( 225.16 - 179.74*logY + 70.22*logY2 - 19.22*logY3 ) 
                   +  Y*Y*( -76.07 + 25.33*logY - 9.17*logY2 - 5.57*logY3 ) 
                   +  Y*Y*Y*( -10.10 - 24.69*logY - 0.30*logY2 - 5.46*logY3 ) 
                   +  Y*Y*Y*Y*( -4.52 - 32.85*logY + 0.72*logY2 - 5.25*logY3 ) 
                   +  Y*Y*Y*Y*Y*( -2.55 - 36.61*logY + 1.06*logY2 - 5.14*logY3 );
    } else {
        throw "Higgs mass is out of range in ThreeLoopEW2QCD::DeltaRho()";
    }
    DeltaRho *= pow(EWSMC.GetXt_alpha(), 2.0) * EWSMC.GetAlsMt()/M_PI;
    return DeltaRho;     
}

double ThreeLoopEW2QCD::DeltaR_rem() const {
    return (0.0);     
}

complex ThreeLoopEW2QCD::deltaRho_rem_l(const StandardModel::lepton l) const {
    complex zero(0.0,0.0,false);
    return zero;      
}

complex ThreeLoopEW2QCD::deltaRho_rem_q(const StandardModel::quark q) const {
    if(q==StandardModel::TOP) return ( complex(0.0,0.0,false) );
    complex zero(0.0,0.0,false);
    return zero;      
}

complex ThreeLoopEW2QCD::deltaKappa_rem_l(const StandardModel::lepton l) const {
    complex zero(0.0,0.0,false);
    return zero;      
}

complex ThreeLoopEW2QCD::deltaKappa_rem_q(const StandardModel::quark q) const {
    if(q==StandardModel::TOP) return ( complex(0.0,0.0,false) );
    complex zero(0.0,0.0,false);
    return zero;       
}









