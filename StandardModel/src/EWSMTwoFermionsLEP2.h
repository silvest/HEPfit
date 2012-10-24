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
 * @brief Form factors G_1, G_2 and G_3 in the processes e^+ e^- -> f fbar at LEP-II
 */
class EWSMTwoFermionsLEP2 {
public:

    /**
     * @brief EWSMTwoFermionsLEP2 constructor
     * @param[in] SM_i reference to a StandardModel object
     * @param[in] bKeepNonUnitary_i true if keeping non-unitary terms
     */
    EWSMTwoFermionsLEP2(const StandardModel& SM_i, const bool bKeepNonUnitary_i=false);

    ////////////////////////////////////////////////////////////////////////  
    
    double G_1(const double s, const double t, const double Mw, 
               const double GammaZ, const double I3f, const double Qf, 
               const double mf, const double mfp, const bool bWeak, 
               const bool bWWbox, const bool bZZbox) const;
    double G_2(const double s, const double t, const double Mw, 
               const double GammaZ, const double I3f, const double Qf, 
               const double mf, const double mfp, const bool bWeak, 
               const bool bWWbox, const bool bZZbox) const;
    double G_3(const double s, const double t, const double Mw, 
               const double GammaZ, const double I3f, const double Qf, 
               const double mf, const double mfp, const bool bWeak, 
               const bool bWWbox, const bool bZZbox) const;

    double G_1_noBox(const double s, const double Mw, const double GammaZ, 
                     const double I3f, const double Qf, const double mf, 
                     const double mfp, const bool bWeak) const;
    double G_2_noBox(const double s, const double Mw, const double GammaZ, 
                     const double I3f, const double Qf, const double mf, 
                     const double mfp, const bool bWeak) const;
    double G_3_noBox(const double s, const double Mw, const double GammaZ, 
                     const double I3f, const double Qf, const double mf, 
                     const double mfp, const bool bWeak) const; 
    
    double G_1_box(const double s, const double t, const double Mw, 
                   const double GammaZ, const double I3f, const double Qf, 
                   const double mf, const double mfp, const bool bWWbox=true, 
                   const bool bZZbox=true) const;
    double G_2_box(const double s, const double t, const double Mw, 
                   const double GammaZ, const double I3f, const double Qf, 
                   const double mf, const double mfp, const bool bWWbox=true, 
                   const bool bZZbox=true) const;
    double G_3_box(const double s, const double t, const double Mw, 
                   const double GammaZ, const double I3f, const double Qf, 
                   const double mf, const double mfp, const bool bWWbox=true, 
                   const bool bZZbox=true) const;    
    
    ////////////////////////////////////////////////////////////////////////  

    complex V_pol(const double s) const;
    
    complex chi_Z(const double s, const double Mw, const double GammaZ) const;

    complex G_e(const double s, const double t, const double Mw, 
                const double I3f, const double Qf, const double mf, 
                const double mfp, const bool bWeak, const bool bWWbox, 
                const bool bZZbox) const;
    complex G_f(const double s, const double t, const double Mw, 
                const double I3f, const double Qf, const double mf, 
                const double mfp, const bool bWeak, const bool bWWbox, 
                const bool bZZbox) const;
    complex G_ef(const double s, const double t, const double Mw, 
                 const double I3f, const double Qf, const double mf, 
                 const double mfp, const bool bWeak, const bool bWWbox, 
                 const bool bZZbox) const;
       
    complex rho_ef(const double s, const double t, const double Mw, 
                   const double I3f, const double Qf, const double mf, 
                   const double mfp, const bool bWeak, const bool bWWbox, 
                   const bool bZZbox) const;
    complex kappa_e(const double s, const double t, const double Mw, 
                    const double I3f, const double Qf, const double mf, 
                    const double mfp, const bool bWeak, const bool bWWbox, 
                    const bool bZZbox) const;
    complex kappa_f(const double s, const double t, const double Mw, 
                    const double I3f, const double Qf, const double mf, 
                    const double mfp, const bool bWeak, const bool bWWbox, 
                    const bool bZZbox) const;
    complex kappa_ef(const double s, const double t, const double Mw, 
                     const double I3f, const double Qf, const double mf, 
                     const double mfp, const bool bWeak, const bool bWWbox, 
                     const bool bZZbox) const;
    
    complex Delta_rho_ef_TOP(const double s, const double t, const double u, 
                             const double Mw, const bool bWWbox) const;
    complex Delta_kappa_e_TOP(const double s, const double t, const double u, 
                              const double Mw, const bool bWWbox) const;
    complex Delta_kappa_f_TOP(const double s, const double t, const double u, 
                              const double Mw, const bool bWWbox) const;
    complex Delta_kappa_ef_TOP(const double s, const double t, const double u, 
                               const double Mw, const bool bWWbox) const;

    complex Delta_rho_ef_WW_hat(const double s, const double t, const double u, 
                                const double Mw, const double I3f) const;
    complex Delta_kappa_e_WW_hat(const double s, const double t, const double u, 
                                 const double Mw, const double I3f) const;
    complex Delta_kappa_f_WW_hat(const double s, const double t, const double u, 
                                 const double Mw, const double I3f) const;
    complex Delta_kappa_ef_WW_hat(const double s, const double t, const double u, 
                                  const double Mw, const double I3f) const;

    complex Delta_rho_ef_WW_TOP_hat(const double s, const double t, const double u, 
                                    const double Mw) const;
    complex Delta_kappa_e_WW_TOP_hat(const double s, const double t, const double u, 
                                     const double Mw) const;
    complex Delta_kappa_f_WW_TOP_hat(const double s, const double t, const double u, 
                                     const double Mw) const;
    complex Delta_kappa_ef_WW_TOP_hat(const double s, const double t, const double u, 
                                      const double Mw) const;
    
    complex Delta_rho_ef_ZZ(const double mu, const double s, const double t, 
                            const double u, const double Mw, const double I3f, 
                            const double Qf) const;
    complex Delta_kappa_e_ZZ(const double mu, const double s, const double t, 
                             const double u, const double Mw, const double I3f, 
                             const double Qf) const;
    complex Delta_kappa_f_ZZ(const double mu, const double s, const double t, 
                             const double u, const double Mw, const double I3f, 
                             const double Qf) const;
    complex Delta_kappa_ef_ZZ(const double mu, const double s, const double t, 
                              const double u, const double Mw, const double I3f, 
                              const double Qf) const;
    
    ////////////////////////////////////////////////////////////////////////  

    // Weak corrections
    //complex I2e(const double s, const double Mw, const bool bWeak) const;
    //complex I2f(const double s, const double Mw, const bool bWeak) const;
    complex DeltaRhobar(const double mu, const double Mw) const;
    complex DeltaRhobarZ(const double mu, const double Mw) const;
    complex D_Z(const double mu, const double s, const double Mw) const;
    complex Pibar_Zgamma(const double mu, const double s, const double Mw) const;
    complex Pibar_gg_bos(const double mu, const double s, const double Mw) const;
    complex F_za_0(const double s, const double Mw) const;
    complex F_Wa_0(const double s, const double Mw) const;    
    complex F_Wa_t(const double s, const double Mw) const;    
    complex F_Wn_0(const double s, const double Mw) const;    
    complex F_Wn_t(const double s, const double Mw) const;    
    complex F_W_0(const double s, const double Mw) const;    
    complex F_W_t(const double s, const double Mw) const;    

    // WW box
    complex B_WW_d_0(const double mu, const double s, const double t, 
                     const double u, const double Mw) const;
    complex B_WW_d(const double mu, const double s, const double t, 
                   const double u, const double Mw) const;
    complex Delta_B_WW_d(const double mu, const double s, const double t, 
                         const double u, const double Mw) const;
    complex B_WW_c_0(const double mu, const double s, const double t, 
                     const double u, const double Mw) const;
    
    // ZZ box
    complex B_ZZ_0(const double mu, const double s, const double t, 
                   const double u) const;
    
    // weak corrections without non-unitary terms
    complex D_Z_hat(const double s, const double Mw) const;
    complex Pibar_Zgamma_hat(const double s, const double Mw) const;
    complex Pibar_gg_bos_hat(const double s, const double Mw) const;
    complex F_Wn_0_hat(const double s, const double Mw) const;    
    complex F_Wn_t_hat(const double s, const double Mw) const;    
    complex F_W_0_hat(const double s, const double Mw) const;    
    complex F_W_t_hat(const double s, const double Mw) const;    
    complex B_WW_d_0_hat(const double s, const double t, const double u, 
                         const double Mw) const;
    complex B_WW_d_0_hat_TEST(const double s, const double t, const double u, 
                              const double Mw) const;
    complex Delta_B_WW_d_hat(const double s, const double t, const double u, 
                             const double Mw) const;
    complex B_WW_c_0_hat(const double s, const double t, const double u, 
                         const double Mw) const;

    
    ////////////////////////////////////////////////////////////////////////  
private:
    bool bDebug; // true for debugging    
    bool bKeepNonUnitary; // true if keeping non-unitary terms
    
    const StandardModel& SM;
    const EWSMcache myCache;
    const EWSMOneLoopEW myOneLoopEW;
    const PVfunctions PV; 
    
};

#endif	/* EWSMTWOFERMIONSLEP2_H */

