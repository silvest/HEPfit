/* 
 * File:   EWSMTwoFermionsLEP2.cpp
 * Author: mishima
 */

#include <cmath>
#include <stdexcept>
#include "EWSMTwoFermionsLEP2.h"


EWSMTwoFermionsLEP2::EWSMTwoFermionsLEP2(const StandardModel& SM_i, 
                                         const bool bDebug_i) : 
                                         SM(SM_i), myOneLoopEW(SM_i) {
    bDebug = bDebug_i;
}


//////////////////////////////////////////////////////////////////////// 

double EWSMTwoFermionsLEP2::sigma_l(const StandardModel::lepton l, 
                                    const double s, const double Mw, const double GammaZ, 
                                    const std::map<std::string, bool>& flags) const {
    return 0.0;
}


double EWSMTwoFermionsLEP2::sigma_q(const StandardModel::quark q, 
                                    const double s, const double Mw, const double GammaZ, 
                                    const std::map<std::string, bool>& flags) const {
    return 0.0;
}


double EWSMTwoFermionsLEP2::AFB_l(const StandardModel::lepton l, 
                                  const double s, const double Mw, const double GammaZ, 
                                  const std::map<std::string, bool>& flags) const {
    return 0.0;
}


double EWSMTwoFermionsLEP2::AFB_q(const StandardModel::quark q, 
                                  const double s, const double Mw, const double GammaZ, 
                                  const std::map<std::string, bool>& flags) const {
    return 0.0;
}


////////////////////////////////////////////////////////////////////////

complex EWSMTwoFermionsLEP2::chi_Z(const double s, const double Mw, 
                                   const double GammaZ) const {
    double Mz = SM.getMz();
    //double cW2 = Mw*Mw/Mz/Mz, sW2 = 1.0 - cW2;
    complex denom = complex(s - Mz*Mz, GammaZ/Mz*s, false);
    double prefactor = SM.getGF()*Mz*Mz/(sqrt(2.0)*8.0*M_PI*SM.getAle());
    
    return ( prefactor*s/denom );
}



////////////////////////////////////////////////////////////////////////




