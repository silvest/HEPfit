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
 * for the processes e^+e^- -> f fbar at LEP-II
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
                    const double Qf) const;
    complex kappa_f(const double s, const double Mw, const double I3f, 
                    const double Qf, const double mfp) const;
    complex kappa_ef(const double s, const double Mw, const double I3f, 
                     const double Qf, const double mfp) const;

    complex I2e(const double s, const double Mw, const double I3f, 
                const double Qf) const;
    complex I2f(const double s, const double Mw, const double I3f, 
                const double Qf, const double mfp) const;
    
    complex G_e(const double s, const double Mw, const double I3f, 
                const double Qf) const;
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
private:
    bool bDebug; // true for debugging    
    
    const StandardModel& SM;
    const EWSMcache myCache;
    const EWSMOneLoopEW myOneLoopEW;
    const PVfunctions PV; 
    
};

#endif	/* EWSMTWOFERMIONSLEP2_H */

