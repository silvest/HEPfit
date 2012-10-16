/* 
 * File:   EWSMTwoFermionsLEP2.cpp
 * Author: mishima
 */

#include <cmath>
#include <stdexcept>
#include "EWSMTwoFermionsLEP2.h"


EWSMTwoFermionsLEP2::EWSMTwoFermionsLEP2(const StandardModel& SM_i) : SM(SM_i), 
        myOneLoopEW(SM_i) {

}


//////////////////////////////////////////////////////////////////////// 

double EWSMTwoFermionsLEP2::sigma_l(const StandardModel::lepton l, 
                                    const double s, const double Mw, const double GammaZ, 
                                    const bool bWEAK, const bool bWEAKBOX, const bool bQED) const {
    return 0.0;
}


double EWSMTwoFermionsLEP2::sigma_q(const StandardModel::quark q, 
                                    const double s, const double Mw, const double GammaZ, 
                                    const bool bWEAK, const bool bWEAKBOX, const bool bQED) const {
    return 0.0;
}


double EWSMTwoFermionsLEP2::AFB_l(const StandardModel::lepton l, 
                                  const double s, const double Mw, const double GammaZ, 
                                  const bool bWEAK, const bool bWEAKBOX, const bool bQED) const {
    return 0.0;
}


double EWSMTwoFermionsLEP2::AFB_q(const StandardModel::quark q, 
                                  const double s, const double Mw, const double GammaZ, 
                                  const bool bWEAK, const bool bWEAKBOX, const bool bQED) const {
    return 0.0;
}


////////////////////////////////////////////////////////////////////////


