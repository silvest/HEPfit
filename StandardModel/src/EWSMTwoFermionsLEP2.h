/* 
 * File:   EWSMTwoFermionsLEP2.h
 * Author: mishima
 */

#ifndef EWSMTWOFERMIONSLEP2_H
#define	EWSMTWOFERMIONSLEP2_H

#include <string>
#include <gslpp.h>
#include <Polylogarithms.h>
#include <PVfunctions.h>
#include "EWSMOneLoopEW.h"
using namespace gslpp;


/**
 * @class EWSMTwoFermionsLEP2
 * @brief Form factors G_1, G_2 and G_3 and the contributions from the box diagrams 
 * in the processes e^+e^- -> f fbar at LEP-II
 * 
 * Formulae used in the current class are calculated in the unitary gauge. 
 */
class EWSMTwoFermionsLEP2 {
public:

    /**
     * @brief EWSMTwoFermionsLEP2 constructor
     * @param[in] SM_i reference to a StandardModel object
     */
    EWSMTwoFermionsLEP2(const StandardModel& SM_i);

    
    ////////////////////////////////////////////////////////////////////////  

    complex Vpol_inv(const double s) const;
    
    complex chi_Z(const double s, const double Mw, const double GammaZ) const;

    complex rho_ef(const double s, const double Mw, const double I3f, 
                   const double Qf, const double mfp) const;
    complex kappa_e(const double s, const double Mw, const double I3f, 
                    const double Qf, const double mfp) const;
    complex kappa_f(const double s, const double Mw, const double I3f, 
                    const double Qf, const double mfp) const;
    complex kappa_ef(const double s, const double Mw, const double I3f, 
                     const double Qf, const double mfp) const;
    
    complex G_e(const double s, const double Mw, const double I3f, 
                const double Qf, const double mfp) const;
    complex G_f(const double s, const double Mw, const double I3f, 
                const double Qf, const double mfp) const;
    complex G_ef(const double s, const double Mw, const double I3f, 
                 const double Qf, const double mfp) const;
    
    double G_1(const double s, const double Mw, const double GammaZ, 
               const double I3f, const double Qf, const double mfp,
               const bool bWeak) const;
    double G_2(const double s, const double Mw, const double GammaZ, 
               const double I3f, const double Qf, const double mfp,
               const bool bWeak) const;
    double G_3(const double s, const double Mw, const double GammaZ, 
               const double I3f, const double Qf, const double mfp,
               const bool bWeak) const; 

    ////////////////////////////////////////////////////////////////////////  
    // Weak corrections

    complex I2e(const double s, const double Mw) const;
    complex I2f(const double s, const double Mw) const;
    
    complex DeltaRhobar(const double mu, const double Mw) const;
    complex DeltaRhobarZ(const double mu, const double Mw) const;
    complex D_Z_hat(const double mu, const double s, const double Mw) const;
    complex Pibar_Zgamma_hat(const double mu, const double s, const double Mw) const;
    complex Pibar_gg_bos_hat(const double mu, const double s, const double Mw) const;

    complex F_za_0(const double s, const double Mw) const;

    complex F_Wa_0(const double s, const double Mw) const;    
    complex F_Wa_t(const double s, const double Mw) const;    
    complex F_Wn_0_hat(const double mu, const double s, const double Mw) const;    
    complex F_Wn_t_hat(const double mu, const double s, const double Mw) const;    
    complex F_W_0_hat(const double mu, const double s, const double Mw) const;    
    complex F_W_t_hat(const double mu, const double s, const double Mw) const;    
    
    
    ////////////////////////////////////////////////////////////////////////  
private:
    bool bDebug; // true for debugging    
    
    const StandardModel& SM;
    const EWSMcache myCache;
    const EWSMOneLoopEW myOneLoopEW;
    const PVfunctions PV; 
    
};

#endif	/* EWSMTWOFERMIONSLEP2_H */

