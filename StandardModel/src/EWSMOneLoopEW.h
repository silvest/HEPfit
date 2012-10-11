/* 
 * File:   EWSMOneLoopEW.h
 * Author: mishima
 */

#ifndef EWSMONELOOPEW_H
#define	EWSMONELOOPEW_H

#include "EWSMcache.h"
using namespace gslpp;


class EWSMOneLoopEW {
    
public:

    /**
     * @brief EWSMOneLoopEW constructor
     * @param[in] cache_i reference to an EWSMcache object
     */
    EWSMOneLoopEW(const EWSMcache& cache_i);


    ////////////////////////////////////////////////////////////////////////

    /**
     * @brief leptonic contribution to alpha
     * @param[in] s invariant mass squared 
     * @return Delta alpha_{lept}^{alpha}
     */
    double DeltaAlpha_l(const double s) const;

    /**
     * @brief top-quark contribution to alpha
     * @param[in] s invariant mass squared 
     * @return Delta alpha_{top}^{alpha}
     */
    double DeltaAlpha_t(const double s) const;

    /**
     * @brief leading contribution to Delta r
     * @param[in] Mw_i the W-boson mass
     * @return Delta rho^{alpha}
     */
    double DeltaRho(const double Mw_i) const;

    /**
     * @brief remainder contribution to Delta r
     * @param[in] Mw_i the W-boson mass
     * @return Delta r_{rem}^{alpha}
     */
    double DeltaR_rem(const double Mw_i) const;

    /**
     * @brief remainder contribution for rho_Z^f and kappa_Z^f
     * @param[in] Mw_i the W-boson mass
     * @return Delta rbar_{rem}^{alpha}
     */
    double DeltaRbar_rem(const double Mw_i) const;

    /**
     * @brief remainder contribution to rho_Z^f for a given u_f
     * @param[in] uf a combination of the unified form factors
     * @param[in] Mw_i the W-boson mass
     * @return delta rho_{rem}^{f, alpha}(u_f)
     */
    complex deltaRho_rem_tmp(const complex u_f, const double Mw_i) const;
    
    /**
     * @brief remainder contribution to rho_Z^l
     * @param[in] l name of lepton
     * @param[in] Mw_i the W-boson mass
     * @return delta rho_{rem}^{l, alpha}
     */
    complex deltaRho_rem_l(const StandardModel::lepton l, const double Mw_i) const;

    /**
     * @brief remainder contribution to rho_Z^q
     * @param[in] q name of quark
     * @param[in] Mw_i the W-boson mass
     * @return delta rho_{rem}^{q, alpha}
     */
    complex deltaRho_rem_q(const StandardModel::quark q, const double Mw_i) const;

    /**
     * @brief remainder contribution to kappa_Z^f for given delta_f and u_f
     * @param[in] deltaf a combination of the effective couplings 
     * @param[in] uf a combination of the unified form factors
     * @param[in] Mw_i the W-boson mass
     * @return delta kappa_{rem}^{f, alpha}(delta_f, u_f)
     */
    complex deltaKappa_rem_tmp(const double deltaf, const complex uf,
                               const double Mw_i) const;    
    
    /**
     * @brief remainder contribution to kappa_Z^l
     * @param[in] l name of lepton
     * @param[in] Mw_i the W-boson mass
     * @return delta kappa_{rem}^{l, alpha}
     */
    complex deltaKappa_rem_l(const StandardModel::lepton l, const double Mw_i) const;
                                                  
    /**
     * @brief remainder contribution to kappa_Z^q
     * @param[in] q name of quark
     * @param[in] Mw_i the W-boson mass
     * @return delta kappa_{rem}^{q, alpha}
     */
    complex deltaKappa_rem_q(const StandardModel::quark q, const double Mw_i) const;
    
    /**
     * @param[in] Qi the electric charge of f_i
     * @param[in] Qj the electric charge of f_j
     * @param[in] Mw_i the W-boson mass
     * @return O(alpha) contribution to the width of W -> f_i bar{f}_j for given Q_i and Q_j
     * @attention masses for virtual fermions are neglected. 
     */
    double rho_GammaW_tmp(const double Qi, const double Qj, 
                          const double Mw_i) const;    
    
    /**
     * @param[in] li name of a neutrino
     * @param[in] lj name of a charged lepton
     * @param[in] Mw_i the W-boson mass
     * @return O(alpha) contribution to rho_ij^W for Gamma_W 
     */
    double rho_GammaW_l(const StandardModel::lepton li, 
                        const StandardModel::lepton lj, 
                        const double Mw_i) const;

    /**
     * @param[in] qi name of a up-type quark
     * @param[in] qj name of a down-type quark
     * @param[in] Mw_i the W-boson mass
     * @return O(alpha) contribution to rho_ij^W for Gamma_W 
     */
    double rho_GammaW_q(const StandardModel::quark qi, 
                        const StandardModel::quark qj, 
                        const double Mw_i) const;

    
    ////////////////////////////////////////////////////////////////////////    
       
    /**
     * @param[in] mu renormalization scale
     * @param[in] s momentum-squared
     * @param[in] Mw_i the W-boson mass
     * @return bosonic contribution to the self-energy function of W boson
     */
    complex SigmaWW_bos(const double mu, const double s, const double Mw_i) const;
 
    /**
     * @param[in] mu renormalization scale
     * @param[in] s momentum-squared
     * @param[in] Mw_i the W-boson mass
     * @return fermionic contribution to the self-energy function of W boson
     */
    complex SigmaWW_fer(const double mu, const double s, const double Mw_i) const;
    
    /**
     * @param[in] mu renormalization scale
     * @param[in] s momentum-squared
     * @param[in] Mw_i the W-boson mass
     * @return bosonic contribution to the self-energy function of Z boson
     */
    complex SigmaZZ_bos(const double mu, const double s, const double Mw_i) const;    
    
    /**
     * @param[in] mu renormalization scale
     * @param[in] s momentum-squared
     * @param[in] Mw_i the W-boson mass
     * @return fermionic contribution to the self-energy function of Z boson
     */
    complex SigmaZZ_fer(const double mu, const double s, const double Mw_i) const;
    
    /**
     * @param[in] mu renormalization scale
     * @param[in] s momentum-squared
     * @param[in] Mw_i the W-boson mass
     * @return bosonic contribution to the self-energy function of photon 
     */
    complex PiGammaGamma_bos(const double mu, const double s, const double Mw_i) const;

     /**
     * @param[in] mu renormalization scale
     * @param[in] s momentum-squared
     * @param[in] l name of lepton
     * @return contribution to the self-energy function of photon from lepton l
     */
    complex PiGammaGamma_fer_l(const double mu, const double s, const StandardModel::lepton l) const;

     /**
     * @param[in] mu renormalization scale
     * @param[in] s momentum-squared
     * @param[in] q name of quark
     * @return contribution to the self-energy function of photon from quark q
     */
    complex PiGammaGamma_fer_q(const double mu, const double s, const StandardModel::quark q) const;    
    
    /**
     * @param[in] mu renormalization scale
     * @param[in] s momentum-squared
     * @return fermionic contribution to the self-energy function of photon 
     */
    complex PiGammaGamma_fer(const double mu, const double s) const;

    /**
     * @param[in] mu renormalization scale
     * @param[in] s momentum-squared
     * @param[in] Mw_i the W-boson mass
     * @return bosonic contribution to the self-energy function of the Z-gamma mixing
     */
    complex PiZgamma_bos(const double mu, const double s, const double Mw_i) const;

    /**
     * @param[in] mu renormalization scale
     * @param[in] s momentum-squared
     * @param[in] Mw_i the W-boson mass
     * @return fermionic contribution to the self-energy function of the Z-gamma mixing
     */    
    complex PiZgamma_fer(const double mu, const double s, const double Mw_i) const;

    
    ////////////////////////////////////////////////////////////////////////   
    
    /**
     * @param[in] mu renormalization scale
     * @param[in] Mw_i the W-boson mass
     * @return bosonic contribution to the wave-function renormalization of W boson with s=Mw^2
     */   
    complex SigmaPrime_WW_bos_Mw2(const double mu, const double Mw_i) const;
    
     /**
     * @param[in] mu renormalization scale
     * @param[in] Mw_i the W-boson mass
     * @return fermionic contribution to the wave-function renormalization of W boson with s=Mw^2
     */   
    complex SigmaPrime_WW_fer_Mw2(const double mu, const double Mw_i) const;   
    
     /**
     * @param[in] mu renormalization scale
     * @param[in] Mw_i the W-boson mass
     * @return bosonic contribution to the wave-function renormalization of Z boson with s=Mz^2
     */   
    complex SigmaPrime_ZZ_bos_Mz2(const double mu, const double Mw_i) const;
    
    /**
     * @param[in] mu renormalization scale
     * @param[in] Mw_i the W-boson mass
     * @return fermionic contribution to the wave-function renormalization of Z boson with s=Mz^2
     */   
    complex SigmaPrime_ZZ_fer_Mz2(const double mu, const double Mw_i) const;     
    
    ////////////////////////////////////////////////////////////////////////       
    
    /**
     * @param[in] mu renormalization scale
     * @param[in] Mw_i the W-boson mass
     * @return Delta Rhobar^{F} for the renormalization scale mu
     */
    double DeltaRhobar(const double mu, const double Mw_i) const;
    
    /**
     * @param[in] mu renormalization scale
     * @param[in] Mw_i the W-boson mass
     * @return Delta Rhobar_W^{F} for the renormalization scale mu
     */
    double DeltaRhobarW(const double mu, const double Mw_i) const;
    
    
    ////////////////////////////////////////////////////////////////////////   
    
    /**
     * @brief TEST function
     * @param[in] Mw_i the W-boson mass
     * @return Delta Rhobar^{bos,F} for mu=Mw
     * @attention The renormalization scale is fixed to be mu=Mw.
     */
    double TEST_DeltaRhobar_bos(const double Mw_i) const;

    /**
     * @brief TEST function
     * @param[in] Mw_i the W-boson mass
     * @return Delta Rhobar_W^{bos,F}
     * @attention The renormalization scale is fixed to be mu=Mw.
     */
    double TEST_DeltaRhobarW_bos(const double Mw_i) const;    


    ////////////////////////////////////////////////////////////////////////    

    /**
     * @param[in] s momentum-squared
     * @param[in] Mw_i the W-boson mass
     * @return abelian-type vertex corrections with virtual Z in the chiral limit
     */
    complex FZa_0(const double s, const double Mw_i) const;

    /**
     * @param[in] s momentum-squared
     * @param[in] Mw_i the W-boson mass
     * @return abelian-type vertex corrections with virtual W in the chiral limit
     */    
    complex FWa_0(const double s, const double Mw_i) const;

    /**
     * @param[in] s momentum-squared
     * @return non-abelian-type vertex corrections with virtual W in the chiral limit
     */     
    complex FbarWa_0(const double s) const;

    /**
     * @param[in] s momentum-squared
     * @param[in] Mw_i the W-boson mass
     * @return non-abelian-type vertex corrections with virtual W in the chiral limit
     */        
    complex FWn_0(const double s, const double Mw_i) const;
    
    /**
     * @param[in] s momentum-squared
     * @param[in] Mw_i the W-boson mass
     * @return abelian-type vertex corrections with virtual W and top-quark for Zbb
     */  
    complex FWa_t(const double s, const double Mw_i) const;

    /**
     * @param[in] s momentum-squared
     * @param[in] Mw_i the W-boson mass
     * @return abelian-type vertex corrections with virtual W and top-quark for Zbb
     */      
    complex FbarWa_t(const double s, const double Mw_i) const;

    /**
     * @param[in] s momentum-squared
     * @param[in] Mw_i the W-boson mass
     * @return non-abelian-type vertex corrections with virtual W and top-quark for Zbb
     */     
    complex FWn_t(const double s, const double Mw_i) const;

    /**
     * @param[in] s momentum-squared
     * @param[in] Mw_i the W-boson mass
     * @return Unified form factor F_Z
     */      
    complex FZ(const double s, const double Mw_i) const;
    
    /**
     * @param[in] s momentum-squared
     * @param[in] l name of lepton 
     * @param[in] Mw_i the W-boson mass
     * @return Unified form factor F_W for the l-lbar channel
     */  
    complex FW_l(const double s, const StandardModel::lepton l, const double Mw_i) const;

    /**
     * @param[in] s momentum-squared
     * @param[in] q name of quark
     * @param[in] Mw_i the W-boson mass
     * @return Unified form factor F_W for the q-qbar channel
     */  
    complex FW_q(const double s, const StandardModel::quark q, const double Mw_i) const;    
    
    ////////////////////////////////////////////////////////////////////////        
    
    /**
     * @brief TEST function
     * @param[in] s momentum-squared
     * @param[in] mf the mass of the fermion in the loop
     * @param[in] Mw_i the W-boson mass
     * @return FWn_0 + FWn_t   
     */
    complex TEST_FWn(const double s, const double mf, const double Mw_i) const;
    
    
    ////////////////////////////////////////////////////////////////////////    

private:
    const EWSMcache& cache;    
    
    
};

#endif	/* EWSMONELOOPEW_H */

