/* 
 * File:   AFBbottom.cpp
 * Author: mishima
 */

#include "AFBbottom.h"


AFBbottom::AFBbottom(const EW& EW_i) : ThObservable(EW_i) {
    std::string Model = EW_i.getSM().ModelName();
    if (Model=="StandardModel" || Model=="SUSY")
        AFB_b = 3.0/4.0*EW_i.A_l(SM.ELECTRON)*EW_i.A_q(SM.BOTTOM);
    else if (Model=="THDM") {
        double S = EW_i.getSM().obliqueS();
        double T = EW_i.getSM().obliqueT();        
        double U = EW_i.getSM().obliqueU();
        double AFB_b_SM = 0.0;
        // Write codes!
        AFB_b = AFB_b_SM;
    } else 
        throw "Error";
}

double AFBbottom::getThValue() {   
    return AFB_b;
}
        

